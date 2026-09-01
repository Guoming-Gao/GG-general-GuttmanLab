from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .classes import CLASS_ORDER,COLORS


def aggregate(signal_tsv, classes_tsv, class_col, coordinate, subset, bootstrap, seed, output_tsv, output_png, title):
    sig=pd.read_csv(signal_tsv,sep="\t"); cls=pd.read_csv(classes_tsv,sep="\t")[["gene_id",class_col]]
    x=sig.merge(cls,on="gene_id",how="inner")
    if subset=="autosomes": x=x[~x.chrom.isin(["chrX","chrY","chrM"])]
    elif subset=="chrX": x=x[x.chrom.eq("chrX")]
    elif subset=="all": x=x[~x.chrom.eq("chrY")]
    rng=np.random.default_rng(seed); rows=[]
    fig,ax=plt.subplots(figsize=(8,5.5))
    for c in CLASS_ORDER:
        z=x[x[class_col].eq(c)]; genes=z.gene_id.unique(); pivot=z.pivot_table(index="gene_id",columns=coordinate,values="enrichment",aggfunc="mean")
        if len(genes)==0: continue
        center=pivot.median(axis=0,skipna=True)
        boots=np.empty((bootstrap,pivot.shape[1]))
        arr=pivot.to_numpy()
        for i in range(bootstrap): boots[i]=np.nanmedian(arr[rng.integers(0,len(arr),len(arr))],axis=0)
        lo,hi=np.nanquantile(boots,[.025,.975],axis=0)
        coord=np.asarray(pivot.columns)
        ax.plot(coord,center,color=COLORS[c],label=f"{c.replace('_',' ')} (n={len(genes):,})")
        ax.fill_between(coord,lo,hi,color=COLORS[c],alpha=.16)
        rows.extend({"class":c,"coordinate":q,"median":m,"ci_low":l,"ci_high":h,"n_genes":len(genes)} for q,m,l,h in zip(coord,center,lo,hi))
    ax.set(xlabel=coordinate,ylabel="median log2(CLAP/Input)",title=title); ax.axhline(0,color="black",lw=.7); ax.legend(fontsize=7,frameon=False,bbox_to_anchor=(1.01,1),loc="upper left")
    fig.tight_layout(); Path(output_png).parent.mkdir(parents=True,exist_ok=True); fig.savefig(output_png,dpi=300,bbox_inches="tight"); plt.close(fig)
    pd.DataFrame(rows).to_csv(output_tsv,sep="\t",index=False)

