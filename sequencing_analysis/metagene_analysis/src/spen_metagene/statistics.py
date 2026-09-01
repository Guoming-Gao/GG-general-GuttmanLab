from __future__ import annotations

import numpy as np
import pandas as pd

from .classes import CLASS_ORDER


def _bh(p):
    p=np.asarray(p,float); order=np.argsort(p); ranked=p[order]; q=ranked*len(p)/np.arange(1,len(p)+1); q=np.minimum.accumulate(q[::-1])[::-1].clip(max=1); out=np.empty_like(q); out[order]=q; return out


def gene_metrics(fixed, scaled, features):
    f=fixed.pivot_table(index="gene_id",columns="bin",values="enrichment",aggfunc="mean",fill_value=0)
    s=scaled.pivot_table(index="gene_id",columns="bin",values="enrichment",aggfunc="mean",fill_value=0)
    z=features.pivot_table(index="gene_id",columns="feature",values="enrichment",aggfunc="mean",fill_value=0)
    out=pd.DataFrame(index=sorted(set(f.index)|set(s.index)|set(z.index)))
    out["first_5kb"]=f.reindex(out.index).loc[:,[c for c in f.columns if 20<=int(c)<70]].mean(axis=1)
    out["first20_minus_remaining80"]=s.reindex(out.index).loc[:,[c for c in s.columns if 20<=int(c)<40]].mean(axis=1)-s.reindex(out.index).loc[:,[c for c in s.columns if 40<=int(c)<120]].mean(axis=1)
    out["tss20_minus_tes20"]=s.reindex(out.index).loc[:,[c for c in s.columns if 20<=int(c)<40]].mean(axis=1)-s.reindex(out.index).loc[:,[c for c in s.columns if 100<=int(c)<120]].mean(axis=1)
    out["first_intron_minus_later"]=z.reindex(out.index).get("first_intron",0)-z.reindex(out.index).get("later_introns",0)
    exon=(z.reindex(out.index).get("first_exon",0)+z.reindex(out.index).get("later_exons",0))/2
    out["intronic_minus_exonic"]=z.reindex(out.index).get("all_introns",0)-exon
    return out.rename_axis("gene_id").reset_index()


def stratified_tests(metrics, classes, class_col, iterations=1000, seed=20260812):
    keep=[c for c in ["gene_id",class_col,"baseMean","gene_length","intron_count","intron_fraction"] if c in classes]
    cov=classes[keep].merge(metrics,on="gene_id")
    # Geometry-derived covariates are joined by the caller when available.
    for col in ["baseMean","gene_length","intron_count","intron_fraction"]:
        if col not in cov: cov[col]=0
        try: cov[col+"_q"]=pd.qcut(cov[col].rank(method="first"),4,labels=False,duplicates="drop")
        except ValueError: cov[col+"_q"]=0
    cov["stratum"]=cov[["baseMean_q","gene_length_q","intron_count_q","intron_fraction_q"]].astype(str).agg("|".join,axis=1)
    rng=np.random.default_rng(seed); rows=[]; metrics_cols=[c for c in metrics if c!="gene_id"]
    for metric in metrics_cols:
        for cls in [c for c in CLASS_ORDER if c!="no_change"]:
            use=cov[cov[class_col].isin([cls,"no_change"])].dropna(subset=[metric]).copy().reset_index(drop=True)
            if not (use[class_col].eq(cls).any() and use[class_col].eq("no_change").any()): continue
            class_values=use.loc[use[class_col].eq(cls),metric].to_numpy()
            ref_values=use.loc[use[class_col].eq("no_change"),metric].to_numpy()
            obs=np.nanmedian(class_values)-np.nanmedian(ref_values); null=[]; boot=[]
            labels=use[class_col].to_numpy().copy()
            for _ in range(iterations):
                perm=labels.copy()
                for _,idx in use.groupby("stratum").indices.items(): perm[idx]=rng.permutation(perm[idx])
                null.append(np.nanmedian(use.loc[perm==cls,metric])-np.nanmedian(use.loc[perm=="no_change",metric]))
                boot.append(np.nanmedian(rng.choice(class_values,len(class_values),replace=True))-np.nanmedian(rng.choice(ref_values,len(ref_values),replace=True)))
            p=(1+np.sum(np.abs(null)>=abs(obs)))/(iterations+1)
            ci=np.nanquantile(boot,[.025,.975])
            rows.append({"metric":metric,"class":cls,"reference":"no_change","median_difference":obs,"bootstrap_ci_low":ci[0],"bootstrap_ci_high":ci[1],"permutation_pvalue":p,"n_class":int(use[class_col].eq(cls).sum()),"n_reference":int(use[class_col].eq("no_change").sum()),"iterations":iterations})
    out=pd.DataFrame(rows)
    if len(out): out["padj_bh"]= _bh(out.permutation_pvalue)
    return out
