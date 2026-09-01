from __future__ import annotations
import json
from pathlib import Path
import pandas as pd


def build(class_files, output_dir):
    out=Path(output_dir); out.mkdir(parents=True,exist_ok=True); counts=[]; changes=[]; concord=[]
    loaded={k:pd.read_csv(v,sep="\t") for k,v in class_files.items()}
    for label,d in loaded.items():
        contrast,readout=label.split("__",1)
        for tag in ["fdr05","fdr10"]:
            vc=d[f"class_{tag}"].value_counts(dropna=False)
            counts.extend({"contrast":contrast,"readout":readout,"threshold":tag,"class":c,"n":int(n)} for c,n in vc.items())
        ch=d.class_fdr05.ne(d.class_fdr10)
        changes.extend({"contrast":contrast,"readout":readout,"gene_id":r.gene_id,"gene_name":getattr(r,"gene_name",None),"class_fdr05":r.class_fdr05,"class_fdr10":r.class_fdr10} for r in d[ch].itertuples())
    for contrast in sorted({k.split("__",1)[0] for k in loaded}):
        a=loaded[f"{contrast}__nascent"][["gene_id","class_fdr05","class_fdr10"]].add_suffix("_nascent").rename(columns={"gene_id_nascent":"gene_id"})
        b=loaded[f"{contrast}__total"][["gene_id","class_fdr05","class_fdr10"]].add_suffix("_total").rename(columns={"gene_id_total":"gene_id"})
        x=a.merge(b,on="gene_id")
        for tag in ["fdr05","fdr10"]:
            tab=pd.crosstab(x[f"class_{tag}_nascent"],x[f"class_{tag}_total"])
            for nrow,row in tab.iterrows():
                for ncol,n in row.items(): concord.append({"contrast":contrast,"threshold":tag,"nascent_class":nrow,"total_class":ncol,"n":int(n)})
    pd.DataFrame(counts).to_csv(out/"class_counts.tsv",sep="\t",index=False)
    pd.DataFrame(changes).to_csv(out/"fdr05_to_fdr10_changes.tsv",sep="\t",index=False)
    pd.DataFrame(concord).to_csv(out/"nascent_total_concordance.tsv",sep="\t",index=False)
