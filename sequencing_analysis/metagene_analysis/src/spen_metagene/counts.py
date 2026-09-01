from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def prepare_counts(expr_dir: str, late_dir: str, output_dir: str) -> None:
    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    samples = []
    for condition in ["ctrl", "Aux1h", "Aux2h", "Dox", "DoxAux24h"]:
        for rep in ["A", "B"]:
            name = f"total__SHA_{condition}_rep{rep}"
            base = Path(expr_dir if condition in {"ctrl", "Aux1h", "Aux2h"} else late_dir)
            f = base / f"{name}.tsv"
            if not f.exists(): raise FileNotFoundError(f)
            samples.append((name, condition, rep, f))
    annotation = None
    for experiment, conds in [("early", {"ctrl", "Aux1h", "Aux2h"}), ("late", {"Dox", "DoxAux24h"})]:
        subset = [x for x in samples if x[1] in conds]
        meta = pd.DataFrame([(n, c, r) for n,c,r,_ in subset], columns=["sample", "condition", "replicate"])
        meta.to_csv(out/f"{experiment}_metadata.tsv", sep="\t", index=False)
        for readout, col in [("total", "all_frag"), ("nascent", "nascent_frag")]:
            mats=[]
            for n,_,_,f in subset:
                d=pd.read_csv(f,sep="\t")
                if annotation is None: annotation=d[["gene_id","gene_name","chrom","strand"]].drop_duplicates("gene_id")
                mats.append(d[["gene_id",col]].rename(columns={col:n}))
            m=mats[0]
            for z in mats[1:]: m=m.merge(z,on="gene_id",how="outer")
            m.fillna(0).to_csv(out/f"{experiment}_{readout}_counts.tsv",sep="\t",index=False)
    annotation.to_csv(out/"gene_annotation.tsv",sep="\t",index=False)


def merge_de_annotation(de_tsv: str, annotation_tsv: str, output_tsv: str) -> None:
    d=pd.read_csv(de_tsv,sep="\t"); a=pd.read_csv(annotation_tsv,sep="\t")
    d.merge(a,on="gene_id",how="left").to_csv(output_tsv,sep="\t",index=False)


def write_manifest(config: dict, output: str) -> None:
    paths = {"gtf": config["gtf"], "refgene_bed": config["refgene_bed"], "tc_mask": config["tc_mask"], "prior_report": config["prior_report"]}
    paths.update({f"late::{k}":v for k,v in config["late_bams"].items()})
    paths.update({f"clap::{k}":v["bam"] for k,v in config["clap_samples"].items()})
    for condition in ["ctrl","Aux1h","Aux2h"]:
        for rep in ["A","B"]:
            label=f"total__SHA_{condition}_rep{rep}"
            paths[f"early_count::{label}"]=str(Path(config["prior_expr_dir"])/f"{label}.tsv")
    import hashlib, os
    rows=[]
    for label,p in paths.items():
        st=os.stat(p); h=hashlib.sha256(); is_bam=str(p).endswith(".bam")
        with open(p,"rb") as fh:
            if is_bam: h.update(fh.read(1024*1024))
            else:
                for block in iter(lambda:fh.read(8*1024*1024),b""): h.update(block)
        rows.append({"label":label,"path":p,"size":st.st_size,"mtime":st.st_mtime,
                     "digest_scope":"first_1MiB" if is_bam else "full_file","sha256":h.hexdigest()})
    pd.DataFrame(rows).to_csv(output,sep="\t",index=False)
