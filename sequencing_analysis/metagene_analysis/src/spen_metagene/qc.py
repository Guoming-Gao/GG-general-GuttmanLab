from __future__ import annotations

from pathlib import Path
import itertools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def build(root):
    root=Path(root); tables=root/"tables"; figures=root/"figures"/"qc"; tables.mkdir(parents=True,exist_ok=True); figures.mkdir(parents=True,exist_ok=True)
    clap=[]
    for context in ["nodox","dox"]:
        for geometry in ["fixed","scaled","features"]:
            d=pd.read_csv(root/"clap_matrices"/"normalized"/f"{context}.{geometry}.tsv",sep="\t")
            x=d[["enrichment_repA","enrichment_repB"]].replace([np.inf,-np.inf],np.nan).dropna()
            clap.append({"context":context,"geometry":geometry,"n_bins":len(x),"pearson":x.corr(method="pearson").iloc[0,1],"spearman":x.corr(method="spearman").iloc[0,1]})
    clap=pd.DataFrame(clap); clap.to_csv(tables/"clap_replicate_qc.tsv",sep="\t",index=False)
    expression=[]
    sources={
      "early_nascent":"Aux1h_vs_ctrl__nascent.tsv", "early_total":"Aux1h_vs_ctrl__total.tsv",
      "late_nascent":"DoxAux24h_vs_Dox__nascent.tsv", "late_total":"DoxAux24h_vs_Dox__total.tsv",
    }
    for label,name in sources.items():
        d=pd.read_csv(root/"de_results"/"normalized"/name,sep="\t").set_index("gene_id"); x=np.log1p(d)
        for a,b in itertools.combinations(x.columns,2):
            expression.append({"dataset":label,"sample_a":a,"sample_b":b,"pearson_log1p":x[a].corr(x[b]),"spearman":x[a].corr(x[b],method="spearman")})
    expression=pd.DataFrame(expression); expression.to_csv(tables/"expression_sample_correlations.tsv",sep="\t",index=False)
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    labels=(clap.context+"/"+clap.geometry); axes[0].bar(np.arange(len(clap)),clap.pearson,color="#0072B2"); axes[0].set_xticks(np.arange(len(clap)),labels,rotation=35,ha="right"); axes[0].set_ylim(-1,1); axes[0].set(ylabel="Pearson r",title="CLAP/Input enrichment replicate concordance")
    means=expression.groupby("dataset").pearson_log1p.mean(); axes[1].bar(np.arange(len(means)),means,color="#D55E00"); axes[1].set_xticks(np.arange(len(means)),means.index,rotation=35,ha="right"); axes[1].set_ylim(0,1); axes[1].set(ylabel="mean pairwise Pearson r",title="Expression sample concordance (log1p normalized counts)")
    fig.tight_layout(); fig.savefig(figures/"replicate_qc.png",dpi=250,bbox_inches="tight"); plt.close(fig)
