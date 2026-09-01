from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

from .classes import CLASS_ORDER, COLORS

FEATURES = [
    ("five_prime_utr", "5′ UTR"),
    ("first_exon", "First Exon"),
    ("cds_exon", "CDS Exons"),
    ("later_exons", "Later Exons"),
    ("first_intron", "First Intron"),
    ("later_introns", "Later Introns"),
    ("all_introns", "All Introns"),
    ("three_prime_utr", "3′ UTR"),
]


def _density(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Scott-bandwidth KDE whose numerical area is one on the plotted grid."""
    y=gaussian_kde(values,bw_method="scott")(grid)
    area=np.trapezoid(y,grid)
    return y/area if area>0 else y


def build_distribution_data(feature_tsv: str, class_tsv: str, class_col: str="class_fdr05", min_kde_n: int=30) -> pd.DataFrame:
    f=pd.read_csv(feature_tsv,sep="\t",compression="infer")
    c=pd.read_csv(class_tsv,sep="\t")[["gene_id","gene_name","baseMean",class_col]]
    # The reused feature export may contain earlier derived class/name fields.
    # The DE class table is authoritative for this plot.
    f=f.drop(columns=["class",class_col,"gene_name","baseMean"],errors="ignore")
    z=f.merge(c,on="gene_id",how="inner").rename(columns={class_col:"class"})
    z=z[z["class"].isin(CLASS_ORDER)&z["feature"].isin(dict(FEATURES))].copy()
    z["panel"]=z["feature"].map(dict(FEATURES)); z["value"]=pd.to_numeric(z["enrichment"],errors="coerce")
    z["value_type"]="stabilized_log2_clap_input"
    eligible=z[["gene_id","gene_name","class","baseMean"]].drop_duplicates("gene_id")
    base=eligible.copy(); base["feature"]="deseq2_basemean"; base["panel"]="Expression"; base["value"]=np.log10(pd.to_numeric(base["baseMean"],errors="coerce").clip(lower=0)+1); base["value_type"]="log10_basemean_plus_1"
    cols=["gene_id","gene_name","class","feature","panel","value","value_type","baseMean"]
    out=pd.concat([z[cols],base[cols]],ignore_index=True).dropna(subset=["value"])
    n=out.groupby(["panel","class"],observed=True)["gene_id"].transform("nunique")
    out["n_genes"]=n.astype(int); out["encoding"]=np.where(out.n_genes>=min_kde_n,"kde","exact_impulse")
    return out


def plot_distribution_data(data: pd.DataFrame, output_prefix: str, title: str, min_kde_n: int=30,
                           enrich_xlim=None, no_change_label="No Change") -> None:
    panels=[x[1] for x in FEATURES]+["Expression"]
    enrich=data[data.value_type.eq("stabilized_log2_clap_input")].value
    # Pooled limits are shared by all eight features. For the review dataset,
    # these quantiles are approximately -4 and +5 and omit only ~0.22% of
    # extreme observations, retaining substantially more tail than 1%/99%.
    enrich_xlim=tuple(np.quantile(enrich,[.0015,.9994])) if enrich_xlim is None else tuple(enrich_xlim)
    fig=plt.figure(figsize=(12,19))
    gs=fig.add_gridspec(10,1,height_ratios=[1]*8+[.55,1.12],hspace=.16)
    axes=[fig.add_subplot(gs[i,0]) for i in range(8)]+[fig.add_subplot(gs[9,0])]
    for i,(ax,panel) in enumerate(zip(axes,panels)):
        d=data[data.panel.eq(panel)]; is_base=panel=="Expression"
        lo,hi=((d.value.min(),d.value.max()) if is_base else enrich_xlim); span=max(hi-lo,1e-6); grid=np.linspace(lo-.03*span,hi+.03*span,600)
        ymax=0.0
        densities=[]
        for cls in CLASS_ORDER:
            v=d.loc[d["class"].eq(cls),"value"].to_numpy(float); n=len(v)
            if n>=min_kde_n and np.unique(v).size>1:
                y=_density(v,grid); densities.append((cls,y,n)); ymax=max(ymax,float(y.max()))
        for cls,y,n in densities:
            ax.fill_between(grid,0,y,color=COLORS[cls],alpha=.20,lw=0)
            ax.plot(grid,y,color=COLORS[cls],lw=1.8)
        ymax=max(ymax,.1)
        for j,cls in enumerate(CLASS_ORDER):
            v=d.loc[(d["class"].eq(cls))&(d["encoding"].eq("exact_impulse")),"value"].to_numpy(float)
            if len(v):
                height=ymax*(.12+.035*j)
                ax.vlines(v,0,height,color=COLORS[cls],lw=2.2,alpha=.95,zorder=5)
                ax.scatter(v,np.full(len(v),height),s=19,color=COLORS[cls],edgecolor="white",linewidth=.4,zorder=6)
        title_y=1.14 if is_base else .88
        ax.text(.012,title_y,panel,transform=ax.transAxes,ha="left",va="top",fontsize=15,fontweight="semibold",clip_on=False)
        ax.set_ylabel("Density",fontsize=15); ax.tick_params(labelsize=15,length=4,width=.8)
        ax.spines[["top","right"]].set_visible(False); ax.grid(axis="x",color="#dddddd",lw=.55,alpha=.65); ax.set_ylim(bottom=0)
        if not is_base:
            ax.set_xlim(enrich_xlim)
            if i<7: ax.tick_params(labelbottom=False)
            else: ax.set_xlabel("SPEN CLAP enrichment over Input, log₂(CLAP/Input)",fontsize=15)
        else:
            ax.set_xlim(grid.min(),grid.max()); ax.set_xlabel("log₁₀(DESeq2 baseMean + 1)",fontsize=15)
    counts=data[data.panel.eq("All Introns")].groupby("class").gene_id.nunique()
    def display(c): return no_change_label if c=="no_change" else c.replace("_"," ").title()
    handles=[Line2D([0],[0],color=COLORS[c],lw=4,label=f"{display(c)} (n={int(counts.get(c,0)):,})") for c in CLASS_ORDER if counts.get(c,0)>0]
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.5,.948),ncol=3,frameon=False,fontsize=15,handlelength=1.8,columnspacing=1.5)
    fig.suptitle(title,fontsize=15,fontweight="bold",y=.975)
    fig.subplots_adjust(left=.12,right=.985,bottom=.05,top=.89)
    prefix=Path(output_prefix); prefix.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(prefix.with_suffix(".png"),dpi=300,bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".pdf"),bbox_inches="tight")
    plt.close(fig)
    data.to_csv(prefix.with_suffix(".tsv.gz"),sep="\t",index=False,compression="gzip")


def plot_distribution_prototype(feature_tsv: str, class_tsv: str, output_prefix: str, title: str="Aux 2 h vs control — total RNA", min_kde_n: int=30) -> None:
    data=build_distribution_data(feature_tsv,class_tsv,"class_fdr05",min_kde_n)
    plot_distribution_data(data,output_prefix,title,min_kde_n)
