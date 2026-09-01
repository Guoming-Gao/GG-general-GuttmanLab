from __future__ import annotations

import base64
import html
from io import BytesIO
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import gaussian_kde, ks_2samp, spearmanr, wasserstein_distance

from .distributions import FEATURES, _density

TERTILES = ["Low", "Middle", "High"]
COLORS = {"Low": "#56B4E9", "Middle": "#999999", "High": "#D55E00"}


def baseline_expression(normalized_tsv: str, columns: list[str]) -> pd.DataFrame:
    d = pd.read_csv(normalized_tsv, sep="\t")
    missing = [c for c in columns if c not in d]
    if missing:
        raise ValueError(f"Missing untreated total-RNA columns: {missing}")
    out = d[["gene_id", *columns]].copy()
    out["untreated_mean_normalized_count"] = out[columns].mean(axis=1)
    out = out[out.untreated_mean_normalized_count > 0].copy()
    out["log10_untreated_expression_plus_1"] = np.log10(out.untreated_mean_normalized_count + 1)
    return out


def assign_tertiles(expression: pd.DataFrame) -> tuple[pd.DataFrame, tuple[float, float]]:
    x = expression.untreated_mean_normalized_count
    q1, q2 = (float(x.quantile(1/3)), float(x.quantile(2/3)))
    out = expression.copy()
    out["expression_tertile"] = np.select([x <= q1, x <= q2], ["Low", "Middle"], default="High")
    return out, (q1, q2)


def feature_data(feature_tsv: str, expression: pd.DataFrame, subset: str) -> pd.DataFrame:
    f = pd.read_csv(feature_tsv, sep="\t", compression="infer")
    keep = ["gene_id", "gene_name", "chrom", "feature", "enrichment"]
    f = f[[c for c in keep if c in f]].drop_duplicates(["gene_id", "feature"])
    d = f.merge(expression, on="gene_id", how="inner")
    if subset == "autosomes": d = d[d.chrom.ne("chrX")]
    elif subset == "chrX": d = d[d.chrom.eq("chrX")]
    d = d[d.feature.isin(dict(FEATURES))].copy()
    d["feature_label"] = d.feature.map(dict(FEATURES))
    d["enrichment"] = pd.to_numeric(d.enrichment, errors="coerce")
    return d.dropna(subset=["enrichment"])


def eligible_gene_ids(feature_tsvs: list[str | Path]) -> set[str]:
    ids=set()
    for path in feature_tsvs:
        d=pd.read_csv(path,sep="\t",compression="infer",usecols=["gene_id"])
        ids.update(d.gene_id.dropna().astype(str))
    return ids


def _limits(d: pd.DataFrame) -> tuple[float, float]:
    return tuple(np.quantile(d.enrichment, [.0015, .9994]))


def plot_tertiles(d: pd.DataFrame, prefix: Path, title: str, limits: tuple[float, float], min_kde_n: int = 30):
    fig = plt.figure(figsize=(12, 19))
    gs = fig.add_gridspec(10, 1, height_ratios=[1]*8+[.55, 1.12], hspace=.16)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(8)] + [fig.add_subplot(gs[9, 0])]
    for i, ((feature, label), ax) in enumerate(zip(FEATURES, axes[:8])):
        z = d[d.feature.eq(feature)]; lo, hi = limits; grid = np.linspace(lo, hi, 600); ymax = .1
        for tertile in TERTILES:
            v = z.loc[z.expression_tertile.eq(tertile), "enrichment"].to_numpy(float)
            if len(v) >= min_kde_n and np.unique(v).size > 1:
                y = _density(v, grid); ymax = max(ymax, y.max())
                ax.fill_between(grid, 0, y, color=COLORS[tertile], alpha=.20, lw=0)
                ax.plot(grid, y, color=COLORS[tertile], lw=1.8)
            elif len(v):
                ax.vlines(v, 0, ymax*.18, color=COLORS[tertile], lw=1.7)
        ax.text(.012, .88, label, transform=ax.transAxes, fontsize=15, fontweight="semibold", va="top")
        ax.set_xlim(limits); ax.set_ylim(bottom=0); ax.set_ylabel("Density", fontsize=15)
        ax.tick_params(labelsize=15); ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", color="#dddddd", lw=.55, alpha=.65)
        if i < 7: ax.tick_params(labelbottom=False)
        else: ax.set_xlabel("SPEN CLAP enrichment over Input, log₂(CLAP/Input)", fontsize=15)
    ax = axes[-1]
    xall = d[["gene_id", "expression_tertile", "log10_untreated_expression_plus_1"]].drop_duplicates("gene_id")
    lo, hi = xall.log10_untreated_expression_plus_1.min(), xall.log10_untreated_expression_plus_1.max()
    grid = np.linspace(lo, hi, 600)
    for tertile in TERTILES:
        v = xall.loc[xall.expression_tertile.eq(tertile), "log10_untreated_expression_plus_1"].to_numpy()
        if len(v) >= min_kde_n and np.unique(v).size > 1:
            y = _density(v, grid); ax.fill_between(grid, 0, y, color=COLORS[tertile], alpha=.2); ax.plot(grid, y, color=COLORS[tertile], lw=1.8)
        elif len(v): ax.vlines(v, 0, .2, color=COLORS[tertile], lw=1.7)
    ax.text(.012, 1.14, "Untreated Total-RNA Expression", transform=ax.transAxes, fontsize=15, fontweight="semibold", va="top")
    ax.set(xlabel="log₁₀(mean untreated normalized count + 1)", ylabel="Density")
    ax.xaxis.label.set_size(15); ax.yaxis.label.set_size(15); ax.tick_params(labelsize=15)
    ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="x", color="#dddddd", lw=.55, alpha=.65)
    counts = xall.groupby("expression_tertile").gene_id.nunique()
    handles = [Line2D([0], [0], color=COLORS[t], lw=4, label=f"{t} Expression (n={counts.get(t,0):,})") for t in TERTILES]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .948), ncol=3, frameon=False, fontsize=15)
    fig.suptitle(title, fontsize=15, fontweight="bold", y=.975)
    fig.subplots_adjust(left=.12, right=.985, bottom=.05, top=.89)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".pdf"), bbox_inches="tight"); plt.close(fig)


def _trend(x, y, grid):
    edges = np.quantile(x, np.linspace(0, 1, 31)); edges = np.unique(edges)
    b = np.clip(np.digitize(x, edges[1:-1]), 0, len(edges)-2)
    bx, by = [], []
    for k in np.unique(b):
        m = b == k
        if m.sum() >= 5: bx.append(np.median(x[m])); by.append(np.median(y[m]))
    if len(bx) < 3: return np.full_like(grid, np.nan)
    return gaussian_filter1d(np.interp(grid, bx, by), 1.2)


def plot_joint(d: pd.DataFrame, prefix: Path, title: str, seed: int, bootstrap: int):
    rng = np.random.default_rng(seed); fig = plt.figure(figsize=(13, 18)); outer = fig.add_gridspec(4, 2, hspace=.28, wspace=.24)
    stat_rows = []
    for n, (feature, label) in enumerate(FEATURES):
        z = d[d.feature.eq(feature)].dropna(subset=["enrichment", "log10_untreated_expression_plus_1"])
        x = z.log10_untreated_expression_plus_1.to_numpy(); y = z.enrichment.to_numpy(); rho = spearmanr(x, y).statistic
        grid = np.linspace(x.min(), x.max(), 160); fitted = _trend(x, y, grid); boots=[]; rhos=[]
        for _ in range(bootstrap):
            ix = rng.integers(0, len(x), len(x)); boots.append(_trend(x[ix], y[ix], grid)); rhos.append(spearmanr(x[ix], y[ix]).statistic)
        boots=np.asarray(boots); lo, hi=np.nanpercentile(boots,[2.5,97.5],axis=0); rlo,rhi=np.nanpercentile(rhos,[2.5,97.5])
        sub=outer[n//2,n%2].subgridspec(2,2,height_ratios=[1,4],width_ratios=[4,1],hspace=.03,wspace=.03)
        top=fig.add_subplot(sub[0,0]); ax=fig.add_subplot(sub[1,0],sharex=top); right=fig.add_subplot(sub[1,1],sharey=ax)
        ax.scatter(x,y,s=4,alpha=.16,color="#4c78a8",linewidths=0,rasterized=True)
        ax.fill_between(grid,lo,hi,color="#D55E00",alpha=.18,lw=0); ax.plot(grid,fitted,color="#D55E00",lw=2)
        gx=np.linspace(x.min(),x.max(),300); gy=np.linspace(y.min(),y.max(),300)
        if len(x)>=30 and np.unique(x).size>1: top.fill_between(gx,0,_density(x,gx),color="#777777",alpha=.35)
        if len(y)>=30 and np.unique(y).size>1: right.fill_betweenx(gy,0,_density(y,gy),color="#777777",alpha=.35)
        top.axis("off"); right.axis("off"); ax.spines[["top","right"]].set_visible(False); ax.tick_params(labelsize=11)
        ax.set_title(f"{label}\nn={len(z):,}; Spearman ρ={rho:.2f} [{rlo:.2f}, {rhi:.2f}]",fontsize=12,loc="left")
        if n//2==3: ax.set_xlabel("log₁₀(mean untreated normalized count + 1)",fontsize=12)
        if n%2==0: ax.set_ylabel("log₂(CLAP/Input)",fontsize=12)
        stat_rows.append({"feature":feature,"feature_label":label,"n_genes":len(z),"spearman_rho":rho,"spearman_ci_low":rlo,"spearman_ci_high":rhi})
    fig.suptitle(title,fontsize=16,fontweight="bold",y=.995); prefix.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(prefix.with_suffix(".png"),dpi=300,bbox_inches="tight"); fig.savefig(prefix.with_suffix(".pdf"),bbox_inches="tight"); plt.close(fig)
    return pd.DataFrame(stat_rows)


def summaries(d: pd.DataFrame) -> tuple[pd.DataFrame,pd.DataFrame]:
    rows=[]; dist=[]
    for (feature,label), z in d.groupby(["feature","feature_label"]):
        vals={t:z.loc[z.expression_tertile.eq(t),"enrichment"].to_numpy() for t in TERTILES}
        for t,v in vals.items():
            rows.append({"feature":feature,"feature_label":label,"tertile":t,"n":len(v),"median":np.median(v) if len(v) else np.nan,"q1":np.quantile(v,.25) if len(v) else np.nan,"q3":np.quantile(v,.75) if len(v) else np.nan,"fraction_positive":np.mean(v>0) if len(v) else np.nan,"fraction_gt_1":np.mean(v>1) if len(v) else np.nan})
        for a,b in [("Low","Middle"),("Low","High"),("Middle","High")]:
            ok=len(vals[a]) and len(vals[b])
            dist.append({"feature":feature,"feature_label":label,"comparison":f"{a} vs {b}","wasserstein_distance":wasserstein_distance(vals[a],vals[b]) if ok else np.nan,"ks_distance":ks_2samp(vals[a],vals[b]).statistic if ok else np.nan})
    return pd.DataFrame(rows),pd.DataFrame(dist)


def _embedded(path: Path) -> str:
    data=base64.b64encode(path.read_bytes()).decode(); return f"data:image/png;base64,{data}"


def run(root: str, output: str, seed: int=20260812, bootstrap: int=100, min_kde_n: int=30):
    root=Path(root); outdir=root/"figures/expression_tertiles"; tabdir=root/"tables/expression_tertiles"; tabdir.mkdir(parents=True,exist_ok=True)
    nodox_feat=root/"metagene/Aux2h_vs_ctrl__total__fdr05__all__features.tsv.gz"
    early=baseline_expression(root/"de_results/normalized/Aux2h_vs_ctrl__total.tsv",["total__SHA_ctrl_repA","total__SHA_ctrl_repB"])
    # Percentiles are defined only after restricting to the No-Dox
    # CLAP-eligible universe.
    early=early[early.gene_id.isin(eligible_gene_ids([nodox_feat]))].copy()
    early,(e1,e2)=assign_tertiles(early)
    contexts=[("nodox_all",early,nodox_feat,"all","Endogenous SPEN CLAP — untreated-expression tertiles")]
    sections=[]; allstats=[]
    for i,(name,expr,feat,subset,title) in enumerate(contexts):
        d=feature_data(str(feat),expr,subset); lim=_limits(d); prefix=outdir/f"{name}__tertile_distributions"; joint=outdir/f"{name}__expression_enrichment_joint"
        plot_tertiles(d,prefix,title,lim,min_kde_n); corr=plot_joint(d,joint,title+" — gene-level relationship",seed+i,bootstrap)
        summ,dist=summaries(d); d.to_csv(tabdir/f"{name}__gene_feature_data.tsv.gz",sep="\t",index=False,compression="gzip")
        summ.to_csv(tabdir/f"{name}__tertile_summary.tsv",sep="\t",index=False); dist.to_csv(tabdir/f"{name}__distribution_distances.tsv",sep="\t",index=False); corr.to_csv(tabdir/f"{name}__correlations.tsv",sep="\t",index=False)
        corr.insert(0,"context",name); allstats.append(corr)
        counts=d[["gene_id","expression_tertile"]].drop_duplicates().expression_tertile.value_counts()
        sections.append((title,prefix.with_suffix(".png"),joint.with_suffix(".png"),counts,lim))
    early.to_csv(tabdir/"nodox_baseline_expression_and_tertiles.tsv.gz",sep="\t",index=False,compression="gzip")
    pd.DataFrame([{"context":"Endogenous No-Dox","lower_cutoff":e1,"upper_cutoff":e2}]).to_csv(tabdir/"tertile_cutoffs.tsv",sep="\t",index=False)
    cross=pd.concat(allstats,ignore_index=True); cross.to_csv(tabdir/"feature_correlations.tsv",sep="\t",index=False)
    def rho(context, feature):
        q=cross[(cross.context.eq(context))&(cross.feature.eq(feature))]
        return float(q.spearman_rho.iloc[0]) if len(q) else np.nan
    finding=(f"In the No-Dox data, untreated expression is inversely associated with CLAP/Input enrichment: "
             f"the association is modest for 5′ UTR (ρ={rho('nodox_all','five_prime_utr'):.2f}) and first exon (ρ={rho('nodox_all','first_exon'):.2f}), "
             f"but stronger for later exons (ρ={rho('nodox_all','later_exons'):.2f}), later introns (ρ={rho('nodox_all','later_introns'):.2f}), and all introns (ρ={rho('nodox_all','all_introns'):.2f}). "
             "Thus the right-shifted enrichment component is concentrated among lower-expression genes and is most expression-associated after the first exon. This supports expression dependence, but the inverse direction and the CLAP/Input pseudocount make a low-coverage ratio effect an important alternative explanation.")
    cards=[]
    for title,p1,p2,counts,lim in sections:
        cards.append(f"<h2>{html.escape(title)}</h2><p><b>Gene counts:</b> "+", ".join(f"{t}={int(counts.get(t,0)):,}" for t in TERTILES)+f". Enrichment display range: {lim[0]:.2f} to {lim[1]:.2f}.</p><figure><img src='{_embedded(p1)}'><figcaption>Rows 1–8 show independently area-normalized gene-level distributions of paired, pseudocount-adjusted SPEN CLAP enrichment over Input. Density height does not encode group abundance. Row 9 verifies tertile separation using only untreated total-RNA normalized counts. Tertiles are expression groups, not differential-expression response classes.</figcaption></figure><figure><img src='{_embedded(p2)}'><figcaption>Each point is one eligible gene. Orange curves are nonlinear binned-median smooths; bands are exploratory 95% intervals from 100 seeded gene bootstraps. Marginal KDEs summarize each axis. Spearman correlations quantify monotonic association and do not establish causality.</figcaption></figure>")
    table=cross[["feature_label","n_genes","spearman_rho","spearman_ci_low","spearman_ci_high"]].round(3).to_html(index=False,classes="data")
    body=f"""<!doctype html><html><head><meta charset='utf-8'><title>Endogenous SPEN CLAP and untreated expression</title><style>body{{font:15px system-ui;max-width:1280px;margin:auto;padding:30px;color:#222}}h1,h2{{color:#17365d}}img{{max-width:100%;border:1px solid #ddd}}figure{{margin:24px 0}}figcaption{{line-height:1.5}}.note{{background:#fff4cc;padding:14px;border-left:5px solid #e69f00}}.result{{background:#eaf4ea;padding:14px;border-left:5px solid #24733c}}.data{{border-collapse:collapse}}.data th,.data td{{border:1px solid #ddd;padding:5px 8px}}</style></head><body><h1>Endogenous SPEN CLAP enrichment by untreated total-RNA expression</h1><h2>Executive summary and hypothesis</h2><p>This focused analysis tests whether the downstream-feature right shift in endogenous SPEN CLAP enrichment is associated with untreated RNA abundance. It contains no Dox data, Aux-response classes, volcano plots, nascent readouts, or individual treatment contrasts.</p><div class='result'><b>Observed result.</b> {finding}</div><div class='note'>The earlier Expression row used DESeq2 baseMean, which averages normalized counts across every sample in a fitted experiment. Here, expression is exclusively the arithmetic mean of untreated control total-RNA replicates A and B. Every Aux-treated sample is excluded.</div><h2>Definitions and validation</h2><p>Genes with zero untreated expression are excluded before percentile calculation. Tertiles are calculated after restricting to the endogenous No-Dox CLAP-eligible universe. Low, Middle, and High use value cutoffs at the 33.3rd and 66.7th percentiles; ties remain together. Cutoffs are {e1:.3f} and {e2:.3f} normalized counts. CLAP enrichment is log₂[(TMM-normalized CLAP CPM + 0.1)/(TMM-normalized Input CPM + 0.1)], averaged across matched No-Dox replicates.</p>{''.join(cards)}<h2>Feature correlation summary</h2>{table}<h2>Interpretation and limitations</h2><p>The result supports expression dependence, specifically an inverse relationship: lower-expression genes carry more of the positive/right-shifted CLAP/Input component after the first exon. It does not by itself prove stronger SPEN binding at low expression. Expression correlates with gene length, transcript architecture, mappability, RNA processing, and Input behavior; moreover, a fixed pseudocount affects ratios most strongly at low coverage. Separate CLAP and Input component modeling is therefore needed to distinguish biology from denominator or low-count effects. KDE modes and bootstrap bands are descriptive, not causal evidence.</p><h2>Reproducibility</h2><p>Exact gene assignments, untreated replicate values, feature measurements, summaries, distances, correlations, PNGs, and vector PDFs are stored under <code>{html.escape(str(outdir))}</code> and <code>{html.escape(str(tabdir))}</code>. Seed: {seed}; bootstrap iterations: {bootstrap}.</p></body></html>"""
    Path(output).parent.mkdir(parents=True,exist_ok=True); Path(output).write_text(body)
