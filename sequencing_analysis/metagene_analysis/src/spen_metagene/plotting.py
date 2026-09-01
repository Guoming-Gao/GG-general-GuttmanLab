from __future__ import annotations

import math
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from .classes import CLASS_ORDER, COLORS


def paired_volcano(input_tsv: str, output_png: str, class_col: str, title: str, delta: float, fdr: float) -> None:
    """Plot MLE and apeglm LFC against BH padj with invariant MLE-class colors."""
    df=pd.read_csv(input_tsv,sep="\t")
    y=-np.log10(pd.to_numeric(df["padj"],errors="coerce").clip(lower=np.nextafter(0,1)))
    fig,axes=plt.subplots(1,2,figsize=(13.5,5.8),sharey=True)
    specs=[("log2FoldChange_mle","DESeq2 maximum-likelihood LFC"),("log2FoldChange_apeglm","apeglm-shrunken LFC")]
    for ai,(ax,(col,label)) in enumerate(zip(axes,specs)):
        x=pd.to_numeric(df[col],errors="coerce")
        for cls in ["not_tested"]+CLASS_ORDER:
            m=df[class_col].eq(cls)
            ax.scatter(x[m],y[m],s=7 if cls!="not_tested" else 4,alpha=.6,c=COLORS[cls],linewidths=0,rasterized=True,
                       label=f"{cls.replace('_',' ')} (n={m.sum():,})" if ai==0 else None)
        ax.axvline(-delta,color="black",ls="--",lw=.8); ax.axvline(delta,color="black",ls="--",lw=.8)
        ax.axhline(-math.log10(fdr),color="black",ls=":",lw=.8); ax.set_xlabel(label); ax.grid(alpha=.12)
        candidates=df.assign(_x=x,_y=y).dropna(subset=["_x","_y"])
        fixed=candidates[candidates.get("gene_name",pd.Series(index=candidates.index,dtype=str)).isin(["Spen","Xist"])]
        top=candidates.assign(_abs=candidates._x.abs()).sort_values(["padj","_abs"],ascending=[True,False]).head(6)
        labels=pd.concat([fixed,top]).drop_duplicates("gene_id")
        for j,(_,r) in enumerate(labels.iterrows()):
            ax.annotate(str(r.get("gene_name",r["gene_id"])),(r._x,r._y),fontsize=7,xytext=(3,3+7*(j%3)),textcoords="offset points")
    axes[0].set_ylabel("−log10(BH-adjusted p-value)")
    axes[0].legend(frameon=False,fontsize=7,loc="upper left",bbox_to_anchor=(0,-.15),ncol=3)
    fig.suptitle(title+" | BH FDR 0.05"); fig.tight_layout()
    out=Path(output_png); out.parent.mkdir(parents=True,exist_ok=True); fig.savefig(out,dpi=300,bbox_inches="tight"); plt.close(fig)
    plotted=df.copy(); plotted["plotted_y"] = y; plotted["class_color_source"]="MLE LFC + BH padj"
    plotted.to_csv(out.with_suffix(".tsv"),sep="\t",index=False)


def ma_plot(input_tsv: str, output_png: str, class_col: str, title: str, lfc_column: str="log2FoldChange_mle") -> None:
    df=pd.read_csv(input_tsv,sep="\t"); fig,ax=plt.subplots(figsize=(7.5,5.5))
    for cls in ["not_tested"]+CLASS_ORDER:
        m=df[class_col].eq(cls); ax.scatter(df.loc[m,"baseMean"],df.loc[m,lfc_column],s=6,alpha=.55,c=COLORS[cls],linewidths=0,rasterized=True)
    ax.axhline(0,color="black",lw=.8); ax.set_xscale("symlog",linthresh=1)
    ax.set(xlabel="DESeq2 baseMean",ylabel="DESeq2 maximum-likelihood log2 fold change",title=title)
    Path(output_png).parent.mkdir(parents=True,exist_ok=True); fig.tight_layout(); fig.savefig(output_png,dpi=250); plt.close(fig)
