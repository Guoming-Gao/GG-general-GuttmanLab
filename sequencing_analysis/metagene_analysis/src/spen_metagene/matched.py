from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

from .distributions import build_distribution_data, plot_distribution_data, FEATURES


def _smd(a, b):
    den=np.sqrt((np.var(a,ddof=1)+np.var(b,ddof=1))/2)
    return (np.mean(a)-np.mean(b))/den if den>0 else np.nan


def select_no_change(base: pd.DataFrame, bin_width=.1, seed=20260812):
    """Select a pooled-direction expression-matched no-change subset without replacement."""
    b=base.drop_duplicates("gene_id").copy()
    b["expression_bin"]=np.floor(b["value"]/bin_width).astype(int)
    target=b[b["class"].ne("no_change")].copy(); controls=b[b["class"].eq("no_change")].copy()
    rng=np.random.default_rng(seed); chosen=[]; rows=[]
    bins=sorted(set(target.expression_bin)|set(controls.expression_bin))
    for k in bins:
        t=target[target.expression_bin.eq(k)]; c=controls[controls.expression_bin.eq(k)].sort_values("gene_id")
        need=len(t); take=min(need,len(c))
        if take:
            ix=rng.choice(len(c),size=take,replace=False); s=c.iloc[np.sort(ix)].copy(); chosen.append(s)
        rows.append({"expression_bin":k,"bin_lower":k*bin_width,"bin_upper":(k+1)*bin_width,
                     "n_directional":need,"n_no_change_available":len(c),"n_no_change_selected":take,
                     "shortage":max(0,need-len(c))})
    selected=pd.concat(chosen,ignore_index=True) if chosen else controls.iloc[0:0].copy()
    selected["selection_seed"]=seed; selected["bin_width"]=bin_width
    return selected,pd.DataFrame(rows),target,controls


def balance_table(target, original, matched):
    rows=[]
    for label,x in [("directional_target",target),("no_change_unmatched",original),("no_change_matched",matched)]:
        v=x.value.to_numpy(float)
        rows.append({"group":label,"n_genes":len(x),"median_log10_basemean_plus_1":np.median(v),
                     "mean_log10_basemean_plus_1":np.mean(v),"sd_log10_basemean_plus_1":np.std(v,ddof=1)})
    out=pd.DataFrame(rows); tv=target.value.to_numpy(float)
    out["smd_vs_directional"]=np.nan; out["ks_vs_directional"]=np.nan
    for label,x in [("no_change_unmatched",original),("no_change_matched",matched)]:
        i=out.index[out.group.eq(label)][0]; xv=x.value.to_numpy(float)
        out.loc[i,"smd_vs_directional"]=_smd(xv,tv); out.loc[i,"ks_vs_directional"]=ks_2samp(xv,tv).statistic
    return out


def feature_effects(full, matched):
    rows=[]
    for feature,label in FEATURES:
        a=full[full.feature.eq(feature)]; m=matched[matched.feature.eq(feature)]
        for direction in ["decrease","elevation"]:
            x=a[a["class"].str.contains(direction)].value; nc=a[a["class"].eq("no_change")].value
            xm=m[m["class"].str.contains(direction)].value; ncm=m[m["class"].eq("no_change")].value
            u=float(x.median()-nc.median()) if len(x) and len(nc) else np.nan
            mm=float(xm.median()-ncm.median()) if len(xm) and len(ncm) else np.nan
            if not np.isfinite(u) or not np.isfinite(mm): status="unavailable"
            elif np.sign(u)!=np.sign(mm) and abs(u)>1e-9 and abs(mm)>1e-9: status="reversed"
            elif abs(u)<1e-9: status="new_or_unchanged"
            else:
                ratio=abs(mm/u)
                status="persisted" if ratio>=.75 else ("weakened" if ratio>=.25 else "attenuated")
            rows.append({"feature":feature,"feature_label":label,"direction":direction,
                         "median_difference_unmatched":u,"median_difference_matched":mm,"trend_status":status,
                         "n_directional":len(xm),"n_matched_no_change":len(ncm)})
    return pd.DataFrame(rows)


def run(feature_tsv, class_tsv, output_prefix, title, bin_width=.1, seed=20260812, min_kde_n=30):
    full=build_distribution_data(feature_tsv,class_tsv,"class_fdr05",min_kde_n)
    base=full[full.value_type.eq("log10_basemean_plus_1")]
    selected,bins,target,controls=select_no_change(base,bin_width,seed)
    keep=set(selected.gene_id)|set(target.gene_id)
    matched=full[full.gene_id.isin(keep)].copy()
    n=matched.groupby(["panel","class"],observed=True)["gene_id"].transform("nunique")
    matched["n_genes"]=n.astype(int); matched["encoding"]=np.where(n>=min_kde_n,"kde","exact_impulse")
    enrich=full[full.value_type.eq("stabilized_log2_clap_input")].value
    xlim=np.quantile(enrich,[.0015,.9994])
    matched["enrichment_xlim_low"]=xlim[0]; matched["enrichment_xlim_high"]=xlim[1]
    plot_distribution_data(matched,output_prefix,title,min_kde_n,xlim,"Expression-Matched No Change")
    prefix=Path(output_prefix); selected.to_csv(str(prefix)+"__selected_no_change.tsv",sep="\t",index=False)
    bins.to_csv(str(prefix)+"__bin_diagnostics.tsv",sep="\t",index=False)
    balance_table(target,controls,selected).to_csv(str(prefix)+"__balance.tsv",sep="\t",index=False)
    feature_effects(full,matched).to_csv(str(prefix)+"__feature_effects.tsv",sep="\t",index=False)
