from __future__ import annotations

from pathlib import Path
import json
import warnings
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .classes import CLASS_ORDER,COLORS
from .statistics import gene_metrics, stratified_tests

FEATURE_ORDER=["five_prime_utr","first_exon","cds_exon","later_exons","first_intron","later_introns","all_introns","three_prime_utr"]
SHAPE_FEATURES=["five_prime_utr","cds_exon","all_introns","three_prime_utr"]


def _feature_widths(models_json):
    models=json.loads(Path(models_json).read_text()); rows=[]
    for gid,t in models.items():
        ex=sorted(t["exons"]); ordered=ex if t["strand"]=="+" else list(reversed(ex)); cds=sorted(t["cds"])
        intr=[[ex[i][1],ex[i+1][0]] for i in range(len(ex)-1)]; oi=intr if t["strand"]=="+" else list(reversed(intr))
        widths={"first_exon":ordered[0][1]-ordered[0][0] if ordered else 0,
                "later_exons":sum(e-s for s,e in ordered[1:]),"cds_exon":sum(e-s for s,e in cds),
                "all_introns":sum(e-s for s,e in intr),"first_intron":oi[0][1]-oi[0][0] if oi else 0,
                "later_introns":sum(e-s for s,e in oi[1:])}
        if cds:
            lo=min(s for s,_ in cds); hi=max(e for _,e in cds)
            low=sum(max(0,min(e,lo)-s) for s,e in ex if s<lo); high=sum(max(0,e-max(s,hi)) for s,e in ex if e>hi)
            widths["five_prime_utr"],widths["three_prime_utr"]=(low,high) if t["strand"]=="+" else (high,low)
        else: widths["five_prime_utr"]=widths["three_prime_utr"]=0
        rows.extend({"gene_id":gid,"feature":f,"width_bp":w} for f,w in widths.items() if w>0)
    return pd.DataFrame(rows)


def _scaled_width_kb(df):
    b=pd.to_numeric(df["bin"]); return np.where((b>=20)&(b<120),df["gene_length"]/100000.0,0.1)


def _violin_box(ax,data,y,order,hue_order,ylabel,title):
    counts=data.groupby(["feature","class"],observed=True).size(); keep=set(counts[counts>=30].index)
    vd=data[data.apply(lambda r:(r["feature"],r["class"]) in keep,axis=1)]
    if len(vd): sns.violinplot(data=vd,x="feature",y=y,hue="class",order=order,hue_order=hue_order,palette=COLORS,cut=0,inner=None,density_norm="width",common_norm=False,linewidth=.5,alpha=.45,ax=ax,legend=False)
    sns.boxplot(data=data,x="feature",y=y,hue="class",order=order,hue_order=hue_order,palette=COLORS,width=.62,dodge=True,gap=.12,showfliers=False,fill=False,linewidth=.65,ax=ax,legend=False)
    small=data[data.apply(lambda r:counts.get((r["feature"],r["class"]),0)<30,axis=1)]
    if len(small): sns.stripplot(data=small,x="feature",y=y,hue="class",order=order,hue_order=hue_order,palette=COLORS,dodge=True,jitter=.12,size=2.2,alpha=.75,ax=ax,legend=False)
    ax.set(xlabel="",ylabel=ylabel,title=title); ax.tick_params(axis="x",rotation=32); ax.grid(axis="y",alpha=.12)


def _feature_figure(x,classes,class_col,feature_widths,stem,title):
    z=x.merge(classes[["gene_id",class_col]],on="gene_id",how="inner").rename(columns={class_col:"class"})
    if feature_widths.empty: feature_widths=z[["gene_id","feature"]].drop_duplicates().assign(width_bp=1000.0)
    z=z.merge(feature_widths,on=["gene_id","feature"],how="inner"); z=z[z["class"].isin(CLASS_ORDER)].copy()
    z["clap_fpkm"]=z["clap_cpm"]/(z["width_bp"]/1000); z["input_fpkm"]=z["input_cpm"]/(z["width_bp"]/1000)
    piv=z[z.feature.isin(SHAPE_FEATURES)].pivot_table(index="gene_id",columns="feature",values="enrichment",aggfunc="mean")
    pos=piv.clip(lower=0); mass=pos.div(pos.sum(axis=1).replace(0,np.nan),axis=0).stack().rename("positive_mass").reset_index()
    z=z.merge(mass,on=["gene_id","feature"],how="left")
    z["signed_component_density"]=np.log2(1+z["clap_fpkm"])
    inp=z.copy(); inp["signed_component_density"]=-np.log2(1+inp["input_fpkm"]); inp["component"]="Input"
    comp=pd.concat([z.assign(component="CLAP"),inp],ignore_index=True)
    present=[f for f in FEATURE_ORDER if f in set(z.feature)]; shape_present=[f for f in SHAPE_FEATURES if f in set(z.feature)]
    fig,axes=plt.subplots(3,1,figsize=(17,15)); hue=[c for c in CLASS_ORDER if c in set(z["class"])]
    _violin_box(axes[0],z,"enrichment",present,hue,"mean paired log2(CLAP/Input)","SPEN occupancy enrichment")
    _violin_box(axes[1],comp,"signed_component_density",present,hue,"signed log2(1 + coverage density)","CLAP and Input coverage density")
    axes[1].axhline(0,color="black",lw=.6)
    shaped=z[z.feature.isin(shape_present)&z.positive_mass.notna()]
    _violin_box(axes[2],shaped,"positive_mass",shape_present,hue,"fraction of within-gene positive enrichment","Within-gene distribution of positive SPEN enrichment")
    handles=[plt.Line2D([0],[0],marker="s",color="none",markerfacecolor=COLORS[c],label=c.replace("_"," "),markersize=8) for c in hue]
    fig.legend(handles=handles,loc="lower center",ncol=max(1,len(hue)),frameon=False); fig.suptitle(title); fig.tight_layout(rect=(0,.04,1,.98))
    png=Path(str(stem)+".png"); fig.savefig(png,dpi=250,bbox_inches="tight"); plt.close(fig)
    z.to_csv(str(stem)+".tsv.gz",sep="\t",index=False,compression="gzip")
    return png,Path(str(stem)+".tsv.gz")


def _subset(x, subset):
    if subset=="autosomes": return x[~x.chrom.isin(["chrX","chrY","chrM"])]
    if subset=="chrX": return x[x.chrom.eq("chrX")]
    return x[~x.chrom.eq("chrY")]


def _curve(pivot, iterations, rng):
    a=pivot.to_numpy(dtype=float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",category=RuntimeWarning)
        center=np.nanmean(a,axis=0)
    if len(a)<2: return center,center,center
    # A deterministic cap prevents very large no-change classes from making
    # 1,000 gene-level resamples intractable while preserving equal gene weight.
    if len(a)>2500: a=a[rng.choice(len(a),2500,replace=False)]
    boot=np.empty((iterations,a.shape[1]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",category=RuntimeWarning)
        for i in range(iterations): boot[i]=np.nanmean(a[rng.integers(0,len(a),len(a))],axis=0)
    lo,hi=np.nanquantile(boot,[.025,.975],axis=0); return center,lo,hi


def _plot_coordinates(geometry, columns):
    """Return publication coordinates while retaining raw bin IDs in tables."""
    if geometry == "fixed":
        # 100-bp bins span -2 kb through +20 kb; use bin midpoints.
        return np.asarray(columns, dtype=float) * 100 - 1950
    if geometry == "features":
        return np.arange(len(columns))
    return np.asarray(columns, dtype=float)


def _decorate_position_axis(axes, geometry):
    if geometry == "fixed":
        for ax in axes:
            ax.axvline(0, color="black", lw=.6, ls=":")
            ax.set_xlabel("position relative to TSS (bp)")
    elif geometry == "scaled":
        for ax in axes:
            ax.axvline(20, color="black", lw=.6, ls=":")
            ax.axvline(120, color="black", lw=.6, ls=":")
            ax.set_xticks([0, 20, 70, 120, 140], ["-2 kb", "TSS", "50%", "TES", "+2 kb"])
            ax.set_xlabel("scaled transcript position")


def make_suite(signals, class_files, output_dir, bootstrap=100, seed=20260812, agreement_files=None, covariates_file=None, models_json=None):
    """Build all contrast/readout/FDR metagenes and machine-readable curve tables."""
    out=Path(output_dir); out.mkdir(parents=True,exist_ok=True); manifest=[]
    annotation_cov=pd.read_csv(covariates_file,sep="\t") if covariates_file else pd.DataFrame()
    feature_widths=_feature_widths(models_json) if models_json else pd.DataFrame()
    for label,class_file in class_files.items():
        contrast,readout=label.split("__",1); late=contrast=="DoxAux24h_vs_Dox"
        context="dox" if late else "nodox"; subsets=["autosomes","chrX","all"] if late else ["all"]
        classes=pd.read_csv(class_file,sep="\t")
        context_signals={g:pd.read_csv(signals[f"{context}_{g}"],sep="\t") for g in ["fixed","scaled","features"]}
        # FDR10 classifications remain in machine-readable gene tables, but
        # displayed metagenes are deliberately restricted to primary FDR05.
        for fdrtag in ["fdr05"]:
            class_col=f"class_{fdrtag}"
            for geometry in ["fixed","scaled","features"]:
                sig=pd.read_csv(signals[f"{context}_{geometry}"],sep="\t"); coord="feature" if geometry=="features" else "bin"
                if geometry=="fixed": sig=sig[sig.gene_length>=20000]
                if geometry=="scaled":
                    width=_scaled_width_kb(sig); sig["clap_fpkm"]=sig["clap_cpm"]/width; sig["input_fpkm"]=sig["input_cpm"]/width
                joined=sig.merge(classes[["gene_id",class_col]],on="gene_id",how="inner")
                for subset in subsets:
                    stem=out/f"{contrast}__{readout}__{fdrtag}__{subset}__{geometry}"
                    if geometry=="features":
                        sx=_subset(sig,subset); png,tsv=_feature_figure(sx,classes,class_col,feature_widths,stem,f"{contrast} | {readout} | {fdrtag} | {subset} | categorical features")
                        manifest.append({"contrast":contrast,"readout":readout,"fdr":fdrtag,"subset":subset,"geometry":geometry,"figure":str(png),"table":str(tsv)}); continue
                    x=_subset(joined,subset); fig,axes=plt.subplots(1,3,figsize=(16,4.8),sharex=True)
                    table=[]; rng=np.random.default_rng(seed)
                    for cls in CLASS_ORDER:
                        z=x[x[class_col].eq(cls)]
                        if z.empty: continue
                        pivots={
                            "enrichment":z.pivot_table(index="gene_id",columns=coord,values="enrichment",aggfunc="mean",fill_value=0),
                            "components":z.pivot_table(index="gene_id",columns=coord,values="clap_cpm",aggfunc="mean",fill_value=0),
                        }
                        p=pivots["enrichment"]
                        if geometry=="features": p=p.reindex(columns=["five_prime_utr","first_exon","cds_exon","later_exons","first_intron","later_introns","all_introns","three_prime_utr"],fill_value=0)
                        center,lo,hi=_curve(p,bootstrap,rng); q=_plot_coordinates(geometry,p.columns)
                        axes[0].plot(q,center,color=COLORS[cls],label=f"{cls.replace('_',' ')} (n={len(p):,})"); axes[0].fill_between(q,lo,hi,color=COLORS[cls],alpha=.12)
                        component_clap="clap_fpkm" if geometry=="scaled" else "clap_cpm"; component_input="input_fpkm" if geometry=="scaled" else "input_cpm"
                        cp=z.pivot_table(index="gene_id",columns=coord,values=component_clap,aggfunc="mean",fill_value=0).reindex(columns=p.columns,fill_value=0); ip=z.pivot_table(index="gene_id",columns=coord,values=component_input,aggfunc="mean",fill_value=0).reindex(columns=p.columns,fill_value=0)
                        axes[1].plot(q,np.nanmean(cp,axis=0),color=COLORS[cls]); axes[1].plot(q,np.nanmean(ip,axis=0),color=COLORS[cls],ls="--",alpha=.8)
                        positive=np.clip(p.to_numpy(dtype=float),0,None); denom=positive.sum(axis=1,keepdims=True); shape=np.divide(positive,denom,out=np.full_like(positive,np.nan),where=denom>0)
                        shape_df=pd.DataFrame(shape,index=p.index,columns=p.columns); sc,sl,sh=_curve(shape_df,bootstrap,rng); axes[2].plot(q,sc,color=COLORS[cls]); axes[2].fill_between(q,sl,sh,color=COLORS[cls],alpha=.12)
                        for qi,m,l,h in zip(p.columns,center,lo,hi): table.append({"class":cls,coord:qi,"mean_enrichment":m,"ci_low":l,"ci_high":h,"n_genes":len(p),"n_positive_shape":int(np.isfinite(shape).any(axis=1).sum())})
                    if geometry=="features":
                        for ax in axes: ax.set_xticks(range(len(p.columns)),list(p.columns),rotation=35,ha="right")
                    _decorate_position_axis(axes,geometry)
                    axes[0].set_ylabel("gene-mean log2(CLAP/Input)"); axes[0].axhline(0,color="black",lw=.6); axes[0].legend(frameon=False,fontsize=6)
                    axes[1].set_ylabel("gene-mean FPKM" if geometry=="scaled" else "gene-mean normalized CPM"); axes[1].set_title("CLAP solid; Input dashed")
                    axes[2].set_ylabel("gene-mean within-gene positive mass"); axes[2].set_title("shape-normalized")
                    fig.suptitle(f"{contrast} | {readout} | {fdrtag} | {subset} | {geometry}")
                    fig.tight_layout(); png=Path(str(stem)+".png"); tsv=Path(str(stem)+".tsv")
                    fig.savefig(png,dpi=250,bbox_inches="tight"); plt.close(fig); pd.DataFrame(table).to_csv(tsv,sep="\t",index=False)
                    manifest.append({"contrast":contrast,"readout":readout,"fdr":fdrtag,"subset":subset,"geometry":geometry,"figure":str(png),"table":str(tsv)})
            for subset in subsets:
                sx={g:_subset(d,subset) for g,d in context_signals.items()}
                metrics=gene_metrics(sx["fixed"],sx["scaled"],sx["features"])
                cov=sx["fixed"][["gene_id","gene_length","intron_count"]].drop_duplicates("gene_id")
                c=classes.merge(cov,on="gene_id",how="left")
                if len(annotation_cov): c=c.merge(annotation_cov.drop(columns=["chrom"],errors="ignore"),on="gene_id",how="left")
                stats=stratified_tests(metrics,c,class_col,bootstrap,seed)
                stats.to_csv(out/f"{contrast}__{readout}__{fdrtag}__{subset}__stratified_tests.tsv",sep="\t",index=False)
    agreement_colors={"concordant_decrease":"#0072B2","stable_no_change":"#999999","concordant_elevation":"#D55E00","transient_or_weak":"#E69F00","discordant_direction":"#CC79A7"}
    for label,path in (agreement_files or {}).items():
        readout,fdrtag=label.split("__",1); classes=pd.read_csv(path,sep="\t")[["gene_id","agreement"]]
        if fdrtag != "fdr05":
            continue
        for geometry in ["scaled"]:
            sig=pd.read_csv(signals[f"nodox_{geometry}"],sep="\t"); coord="feature" if geometry=="features" else "bin"
            if geometry=="fixed": sig=sig[sig.gene_length>=20000]
            x=sig.merge(classes,on="gene_id")
            fig,ax=plt.subplots(figsize=(8,5.2)); table=[]; rng=np.random.default_rng(seed)
            for cls,color in agreement_colors.items():
                z=x[x.agreement.eq(cls)]; p=z.pivot_table(index="gene_id",columns=coord,values="enrichment",aggfunc="mean",fill_value=0)
                if p.empty: continue
                center,lo,hi=_curve(p,bootstrap,rng); q=_plot_coordinates(geometry,p.columns)
                ax.plot(q,center,color=color,label=f"{cls.replace('_',' ')} (n={len(p):,})"); ax.fill_between(q,lo,hi,color=color,alpha=.14)
                table.extend({"agreement":cls,coord:q0,"mean_enrichment":m,"ci_low":l,"ci_high":h,"n_genes":len(p)} for q0,m,l,h in zip(p.columns,center,lo,hi))
            if geometry=="features" and table:
                cols=list(x[coord].drop_duplicates()); ax.set_xticks(range(len(cols)),cols,rotation=35,ha="right")
            _decorate_position_axis([ax],geometry)
            ax.axhline(0,color="black",lw=.6); ax.set(ylabel="gene-mean log2(CLAP/Input)",title=f"Early 1 h/2 h agreement | {readout} | {fdrtag} | {geometry}"); ax.legend(frameon=False,fontsize=7,bbox_to_anchor=(1.01,1),loc="upper left")
            fig.tight_layout(); stem=f"early_agreement__{readout}__{fdrtag}__all__{geometry}"; png=out/f"{stem}.png"; tsv=out/f"{stem}.tsv"; fig.savefig(png,dpi=250,bbox_inches="tight"); plt.close(fig); pd.DataFrame(table).to_csv(tsv,sep="\t",index=False)
            manifest.append({"contrast":"early_agreement","readout":readout,"fdr":fdrtag,"subset":"all","geometry":geometry,"figure":str(png),"table":str(tsv)})
    pd.DataFrame(manifest).to_csv(out/"metagene_manifest.tsv",sep="\t",index=False)
