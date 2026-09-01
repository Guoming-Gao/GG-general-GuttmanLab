from __future__ import annotations
import json, os
from pathlib import Path
import pandas as pd
import pysam


def _bam_status(path):
    status="ok"; detail=[]
    try:
        with pysam.AlignmentFile(path) as b:
            refs=set(b.references); detail.append(f"contigs={len(refs)}")
            if "chr1" not in refs: status="failed"; detail.append("chr1_missing")
            paired=0; n=0
            for r in b.fetch(until_eof=True):
                n+=1; paired+=int(r.is_paired)
                if n>=1000: break
            frac=paired/n if n else 0; detail.append(f"paired_fraction={frac:.3f}")
            if n==0 or frac<.9: status="failed"
        if not (Path(path+".bai").exists() or Path(path).with_suffix(".bai").exists()):
            status="failed"; detail.append("index_missing")
    except Exception as e: status="failed"; detail.append(repr(e))
    return status,";".join(detail)


def validate_config(config, output):
    rows=[]
    expected_cols={"gene_id","gene_name","chrom","strand","all_frag","nascent_frag"}
    for condition in ["ctrl","Aux1h","Aux2h"]:
        for rep in ["A","B"]:
            label=f"total__SHA_{condition}_rep{rep}"; path=Path(config["prior_expr_dir"])/f"{label}.tsv"
            try:
                d=pd.read_csv(path,sep="\t")
                numeric=d[["all_frag","nascent_frag"]]
                ok=expected_cols.issubset(d.columns) and d.gene_id.is_unique and numeric.ge(0).all().all() and (numeric%1==0).all().all()
                detail=f"genes={len(d)};all_frag={int(d.all_frag.sum())};nascent_frag={int(d.nascent_frag.sum())}"
            except Exception as e: ok=False; detail=repr(e)
            rows.append({"type":"early_count","label":label,"path":str(path),"status":"ok" if ok else "failed","detail":detail})
    for label,entry in config["clap_samples"].items():
        path=entry["bam"]; status,detail=_bam_status(path)
        rows.append({"type":"clap_bam","label":label,"path":path,"status":status,"detail":detail})
    for label,path in config["late_bams"].items():
        status,detail=_bam_status(path)
        rows.append({"type":"late_bam","label":label,"path":path,"status":status,"detail":detail})
    for label in ["gtf","refgene_bed","tc_mask","prior_report"]:
        p=config[label]; rows.append({"type":"reference","label":label,"path":p,"status":"ok" if Path(p).exists() else "failed","detail":""})
    pairs={(v["context"],v["replicate"]):set() for v in config["clap_samples"].values()}
    for v in config["clap_samples"].values(): pairs[(v["context"],v["replicate"])].add(v["assay"])
    for key,assays in sorted(pairs.items()):
        rows.append({"type":"sample_pairing","label":"::".join(key),"path":"","status":"ok" if assays=={"CLAP","Input"} else "failed","detail":",".join(sorted(assays))})
    late_labels=set(config["late_bams"])
    expected={f"total__SHA_{c}_rep{r}" for c in ["Dox","DoxAux24h"] for r in ["A","B"]}
    rows.append({"type":"sample_labels","label":"late_design","path":"","status":"ok" if late_labels==expected else "failed","detail":f"observed={sorted(late_labels)}"})
    try:
        with Path(config["gtf"]).open() as fh:
            header="".join(next(fh,"") for _ in range(5))
        version_ok="vM25" in header or "version M25" in header
    except Exception: version_ok=False
    rows.append({"type":"annotation_version","label":"GENCODE_mm10_vM25","path":config["gtf"],"status":"ok" if version_ok else "failed","detail":"header audit"})
    d=pd.DataFrame(rows); Path(output).parent.mkdir(parents=True,exist_ok=True); d.to_csv(output,sep="\t",index=False)
    if d.status.ne("ok").any(): raise RuntimeError("Input validation failed; see "+output)
