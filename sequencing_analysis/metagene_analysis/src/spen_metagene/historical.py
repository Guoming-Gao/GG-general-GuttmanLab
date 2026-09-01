from __future__ import annotations

from pathlib import Path
import pandas as pd


def summarize(enrichment_dir: str, output: str) -> None:
    root=Path(enrichment_dir); rows=[]
    files=sorted(root.glob("*endospen*_CA.tsv"))+sorted(root.glob("*endospen*_EBS.tsv"))
    for f in files:
        method="CA" if f.name.endswith("_CA.tsv") else "EBS"
        if method=="CA":
            d=pd.read_csv(f,sep="\t")
            q=d[(d["sample count"]>=5)&(d["enrichment (intra)"]>2)&(d["enrichment (inter)"]>2)&(d["p-val (intra)"]<.01)&(d["p-val (inter)"]<.01)]
            genes=q.loc[q["feature name"]!="intergenic","feature name"].astype(str).str.split("|").str[-1].str.replace(r"_(exon|intron)$","",regex=True)
        else:
            d=pd.read_csv(f,sep="\t",header=None); q=d
            genes=q[0].astype(str).str.split("|").str[-1].str.replace(r"_(exon|intron)$","",regex=True)
        rows.append({
            "file":str(f),"method":method,"window_bp":100,
            "permutations":10 if method=="EBS" else pd.NA,
            "filter":("CLAP>=5; intra>2; inter>2; intra_p<0.01; inter_p<0.01"
                      if method=="CA" else "caller output (unfiltered)"),
            "all_rows":len(d),"retained_rows":len(q),"unique_genes":genes.nunique(),
        })
    pd.DataFrame(rows).to_csv(output,sep="\t",index=False)
