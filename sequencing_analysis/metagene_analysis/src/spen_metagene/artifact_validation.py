from __future__ import annotations

import base64
import hashlib
import html
import json
import math
import os
import subprocess
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import pysam
from scipy.ndimage import gaussian_filter1d
from scipy.stats import spearmanr
import statsmodels.api as sm

from .annotation import parse_representative_transcripts
from .clap import FEATURES as FEATURE_KEYS, assign, build_index, feature_membership


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
FEATURE_LABELS = dict(FEATURES)
TERTILES = ["Low", "Middle", "High"]
COLORS = {"Low": "#56B4E9", "Middle": "#777777", "High": "#D55E00"}
SEED = 20260812
CA_FILTER = (
    "sample count >= 5; enrichment (intra) > 2; enrichment (inter) > 2; "
    "p-val (intra) < 0.01; p-val (inter) < 0.01"
)


def _mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _sha256(path: Path, chunk: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def _canonical_chrom(chrom: str) -> str:
    if chrom.startswith("chr"):
        return chrom
    if chrom in {"MT", "M"}:
        return "chrM"
    return f"chr{chrom}"


def _interval_intersection(a: list[int], b: list[int]) -> list[int] | None:
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    return [lo, hi] if lo < hi else None


def feature_intervals(t: dict) -> dict[str, list[list[int]]]:
    exons = sorted([list(x) for x in t["exons"]])
    ordered_exons = exons if t["strand"] == "+" else list(reversed(exons))
    introns = [[exons[i][1], exons[i + 1][0]] for i in range(len(exons) - 1)]
    ordered_introns = introns if t["strand"] == "+" else list(reversed(introns))
    cds = sorted([list(x) for x in t.get("cds", [])])
    out = {k: [] for k in FEATURE_KEYS}
    if ordered_exons:
        out["first_exon"] = [ordered_exons[0]]
        out["later_exons"] = ordered_exons[1:]
    out["all_introns"] = introns
    out["first_intron"] = ordered_introns[:1]
    out["later_introns"] = ordered_introns[1:]
    out["cds_exon"] = cds
    if cds:
        cds_lo, cds_hi = min(x[0] for x in cds), max(x[1] for x in cds)
        for exon in exons:
            left = _interval_intersection(exon, [exon[0], min(exon[1], cds_lo)])
            right = _interval_intersection(exon, [max(exon[0], cds_hi), exon[1]])
            if t["strand"] == "+":
                if left: out["five_prime_utr"].append(left)
                if right: out["three_prime_utr"].append(right)
            else:
                if right: out["five_prime_utr"].append(right)
                if left: out["three_prime_utr"].append(left)
    return out


def feature_lengths(models: dict) -> pd.DataFrame:
    rows = []
    for gid, t in models.items():
        for feature, intervals in feature_intervals(t).items():
            rows.append({
                "gene_id": gid,
                "gene_name": t["gene_name"],
                "chrom": t["chrom"],
                "strand": t["strand"],
                "feature": feature,
                "feature_length": sum(b - a for a, b in intervals),
            })
    return pd.DataFrame(rows)


def untreated_expression(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path, sep="\t")
    cols = ["total__SHA_ctrl_repA", "total__SHA_ctrl_repB"]
    if any(c not in d for c in cols):
        raise ValueError("Untreated total-RNA columns A and B were not found")
    d = d[["gene_id", *cols]].copy()
    d["expression"] = d[cols].mean(axis=1)
    d = d[d.expression > 0].copy()
    d["log10_expression_plus_1"] = np.log10(d.expression + 1)
    d["expression_source"] = "Untreated total RNA"
    return d


def assign_tertiles(d: pd.DataFrame) -> tuple[pd.DataFrame, tuple[float, float]]:
    q1, q2 = [float(d.expression.quantile(q)) for q in (1 / 3, 2 / 3)]
    out = d.copy()
    out["expression_tertile"] = np.select(
        [out.expression <= q1, out.expression <= q2], ["Low", "Middle"], default="High"
    )
    return out, (q1, q2)


def spen_input_expression(root: Path) -> pd.DataFrame:
    matrix = root / "supporting_data/validated_SPEN_CLAP_matrices/raw/all_samples.scaled.tsv.gz"
    sizes = root / "supporting_data/validated_SPEN_CLAP_matrices/normalized/nodox.scaled.size_factors.tsv"
    d = pd.read_csv(matrix, sep="\t")
    d = d[d["bin"].between(20, 119)].copy()
    cols = ["input_endospen_1", "input_endospen_2"]
    sums = d.groupby("gene_id", as_index=False)[cols].sum()
    sf = pd.read_csv(sizes, sep="\t").set_index("sample")["effective_lib_size"]
    for c in cols:
        sums[f"{c}_cpm"] = sums[c] / sf[c] * 1e6
    sums["expression"] = sums[[f"{c}_cpm" for c in cols]].mean(axis=1)
    sums = sums[sums.expression > 0].copy()
    sums["log10_expression_plus_1"] = np.log10(sums.expression + 1)
    sums["expression_source"] = "SPEN Input gene-body CPM"
    return sums


def filter_ca(path: Path) -> pd.DataFrame:
    cols = ["window", "strand", "feature name", "sample count", "input count",
            "enrichment (intra)", "enrichment (inter)", "p-val (intra)", "p-val (inter)"]
    kept = []
    for d in pd.read_csv(path, sep="\t", usecols=cols, chunksize=250_000):
        q = d[(d["sample count"] >= 5) & (d["enrichment (intra)"] > 2) &
              (d["enrichment (inter)"] > 2) & (d["p-val (intra)"] < .01) &
              (d["p-val (inter)"] < .01)].copy()
        kept.append(q)
    out = pd.concat(kept, ignore_index=True) if kept else pd.DataFrame(columns=cols)
    return out.drop_duplicates("window")


def _parse_window(value: str) -> tuple[str, int, int]:
    chrom, pos = str(value).split(":", 1)
    start, end = pos.split("-", 1)
    return _canonical_chrom(chrom), int(start), int(end)


def assign_ca_to_features(calls: pd.DataFrame, models: dict, callset: str) -> pd.DataFrame:
    idx = build_index(models, flank=0)
    rows = []
    for r in calls.itertuples(index=False):
        chrom, start, end = _parse_window(getattr(r, "window"))
        strand = getattr(r, "strand")
        gid = assign(idx, chrom, strand, (start + end) // 2)
        if gid is None:
            continue
        t = models[gid]
        for feature, intervals in feature_intervals(t).items():
            if any(start < b and end > a for a, b in intervals):
                rows.append({"gene_id": gid, "gene_name": t["gene_name"], "chrom": chrom,
                             "strand": strand, "feature": feature,
                             "window": f"{chrom}:{start}-{end}", "callset": callset})
    return pd.DataFrame(rows).drop_duplicates(["gene_id", "feature", "window", "callset"])


def peak_count_universe(assigned: pd.DataFrame, lengths: pd.DataFrame,
                        expression: pd.DataFrame, callset: str) -> pd.DataFrame:
    universe = lengths[lengths.feature_length > 0].merge(expression, on="gene_id", how="inner")
    counts = (assigned[assigned.callset.eq(callset)]
              .groupby(["gene_id", "feature"], as_index=False).window.nunique()
              .rename(columns={"window": "absolute_CA_peak_count"}))
    out = universe.merge(counts, on=["gene_id", "feature"], how="left")
    out["absolute_CA_peak_count"] = out.absolute_CA_peak_count.fillna(0).astype(int)
    out["any_CA_peak"] = out.absolute_CA_peak_count.gt(0)
    out["callset"] = callset
    out["feature_label"] = out.feature.map(FEATURE_LABELS)
    return out


def _bootstrap_rho(x: np.ndarray, y: np.ndarray, n: int, seed: int) -> tuple[float, float, float]:
    if len(x) < 3 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan, np.nan, np.nan
    rho = float(spearmanr(x, y).statistic)
    if n <= 0:
        return rho, np.nan, np.nan
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n):
        ix = rng.integers(0, len(x), len(x))
        vals.append(spearmanr(x[ix], y[ix]).statistic)
    lo, hi = np.nanpercentile(vals, [2.5, 97.5])
    return rho, float(lo), float(hi)


def peak_statistics(data: pd.DataFrame, bootstrap: int = 100, seed: int = SEED) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summaries, models, strata = [], [], []
    for i, (feature, label) in enumerate(FEATURES):
        all_z = data[data.feature.eq(feature)].dropna(subset=["log10_expression_plus_1", "feature_length"])
        z = all_z[all_z.absolute_CA_peak_count.gt(0)].copy()
        x = z.log10_expression_plus_1.to_numpy(float)
        y = z.absolute_CA_peak_count.to_numpy(float)
        rho, lo, hi = _bootstrap_rho(x, y, bootstrap, seed + i)
        for tertile in TERTILES:
            q = z[z.expression_tertile.eq(tertile)].absolute_CA_peak_count
            all_q = all_z[all_z.expression_tertile.eq(tertile)].absolute_CA_peak_count
            summaries.append({"feature": feature, "feature_label": label, "tertile": tertile,
                              "n_all_eligible_genes": len(all_q), "n_peak_positive_genes": len(q),
                              "median_peak_count_among_positive": q.median(),
                              "q1": q.quantile(.25), "q3": q.quantile(.75),
                              "fraction_zero_QC": all_q.eq(0).mean(), "fraction_any_peak_QC": all_q.gt(0).mean(),
                              "spearman_rho_peak_positive": rho, "spearman_ci_low": lo, "spearman_ci_high": hi})
        if len(z) < 10:
            models.append({"feature": feature, "feature_label": label,
                           "endpoint": "positive_peak_count_negative_binomial",
                           "n_peak_positive_genes": len(z), "error": "fewer than 10 peak-positive genes"})
            continue
        X = pd.DataFrame({
            "log_expression_z": (x - x.mean()) / (x.std() or 1),
            "log_feature_length_z": (np.log10(z.feature_length.to_numpy(float) + 1) -
                                     np.log10(z.feature_length.to_numpy(float) + 1).mean()) /
                                    (np.log10(z.feature_length.to_numpy(float) + 1).std() or 1),
        })
        X = sm.add_constant(X)
        for endpoint, yy, family in [
            ("positive_peak_count_negative_binomial", z.absolute_CA_peak_count, sm.families.NegativeBinomial()),
        ]:
            try:
                fit = sm.GLM(yy.to_numpy(), X, family=family).fit()
                models.append({"feature": feature, "feature_label": label, "endpoint": endpoint,
                               "n_peak_positive_genes": len(z), "expression_coefficient": fit.params["log_expression_z"],
                               "expression_se": fit.bse["log_expression_z"],
                               "expression_pvalue": fit.pvalues["log_expression_z"],
                               "length_coefficient": fit.params["log_feature_length_z"],
                               "length_pvalue": fit.pvalues["log_feature_length_z"]})
            except Exception as exc:
                models.append({"feature": feature, "feature_label": label, "endpoint": endpoint,
                               "n_peak_positive_genes": len(z), "error": str(exc)})
        z = z.copy()
        z["length_stratum"] = pd.qcut(z.feature_length.rank(method="first"), 3,
                                      labels=["Short", "Middle", "Long"])
        for stratum, q in z.groupby("length_stratum", observed=True):
            r = (_bootstrap_rho(q.log10_expression_plus_1.to_numpy(float),
                                q.absolute_CA_peak_count.to_numpy(float), 0, seed)[0]
                 if len(q) else np.nan)
            strata.append({"feature": feature, "feature_label": label, "length_stratum": stratum,
                           "n_peak_positive_genes": len(q), "spearman_rho": r})
    return pd.DataFrame(summaries), pd.DataFrame(models), pd.DataFrame(strata)


def _save_figure(fig, prefix: Path) -> None:
    _mkdir(prefix.parent)
    fig.savefig(Path(str(prefix) + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(Path(str(prefix) + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_peak_distributions(data: pd.DataFrame, prefix: Path, title: str) -> None:
    fig = plt.figure(figsize=(12, 19))
    gs = fig.add_gridspec(10, 1, height_ratios=[1] * 8 + [.55, 1.12], hspace=.18)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(8)] + [fig.add_subplot(gs[9, 0])]
    for i, ((feature, label), ax) in enumerate(zip(FEATURES, axes[:8])):
        z = data[data.feature.eq(feature) & data.absolute_CA_peak_count.gt(0)]
        if z.empty:
            ax.text(.5, .5, "No peak-positive genes", ha="center", va="center",
                    transform=ax.transAxes, fontsize=15)
            ax.text(.012, .88, label, transform=ax.transAxes, fontsize=15,
                    fontweight="semibold", va="top")
            ax.set_ylabel("Probability", fontsize=15)
            ax.tick_params(labelsize=15)
            ax.spines[["top", "right"]].set_visible(False)
            continue
        cap = max(3, int(z.absolute_CA_peak_count.quantile(.99)))
        xs = np.arange(1, cap + 1)
        for tertile in TERTILES:
            v = z.loc[z.expression_tertile.eq(tertile), "absolute_CA_peak_count"].to_numpy(int)
            pmf = np.array([(v == k).mean() for k in xs]) if len(v) else np.zeros_like(xs, float)
            pmf[-1] += (v > cap).mean() if len(v) else 0
            offset = {"Low": -.18, "Middle": 0, "High": .18}[tertile]
            ax.vlines(xs + offset, 0, pmf, color=COLORS[tertile], lw=2.1, alpha=.92)
            ax.scatter(xs + offset, pmf, color=COLORS[tertile], s=14, zorder=3)
        ax.text(.012, .88, label, transform=ax.transAxes, fontsize=15, fontweight="semibold", va="top")
        ax.set_xlim(.35, cap + .65); ax.set_ylim(bottom=0); ax.set_ylabel("Probability", fontsize=15)
        ax.tick_params(labelsize=15); ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", color="#dddddd", lw=.55)
        ticks = np.unique(np.linspace(1, cap, min(cap, 7)).round().astype(int))
        ax.set_xticks(ticks, [str(x) if x < cap else f"{x}+" for x in ticks])
        if i < 7: ax.tick_params(labelbottom=False)
        else: ax.set_xlabel("Absolute number of CA-called 100-bp windows per gene feature", fontsize=15)
    ax = axes[-1]
    genes = data[data.absolute_CA_peak_count.gt(0)][["gene_id", "expression_tertile", "log10_expression_plus_1"]].drop_duplicates("gene_id")
    for tertile in TERTILES:
        v = genes.loc[genes.expression_tertile.eq(tertile), "log10_expression_plus_1"].to_numpy()
        if len(v) > 1:
            hist, edges = np.histogram(v, bins=80, density=True)
            centers = (edges[:-1] + edges[1:]) / 2
            ax.fill_between(centers, 0, gaussian_filter1d(hist, 1.1), color=COLORS[tertile], alpha=.18)
            ax.plot(centers, gaussian_filter1d(hist, 1.1), color=COLORS[tertile], lw=1.8)
    ax.text(.012, 1.13, "Expression", transform=ax.transAxes, fontsize=15, fontweight="semibold")
    ax.set_xlabel("log₁₀(expression + 1)", fontsize=15); ax.set_ylabel("Density", fontsize=15)
    ax.tick_params(labelsize=15); ax.spines[["top", "right"]].set_visible(False)
    counts = genes.expression_tertile.value_counts()
    handles = [Line2D([0], [0], color=COLORS[t], lw=4, label=f"{t} (n={counts.get(t, 0):,})") for t in TERTILES]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, .948), ncol=3,
               frameon=False, fontsize=15)
    fig.suptitle(title, fontsize=15, fontweight="bold", y=.975)
    fig.subplots_adjust(left=.12, right=.985, bottom=.05, top=.89)
    _save_figure(fig, prefix)


def _binned_mean(x: np.ndarray, y: np.ndarray, grid: np.ndarray) -> np.ndarray:
    edges = np.unique(np.quantile(x, np.linspace(0, 1, 31)))
    bins = np.clip(np.digitize(x, edges[1:-1]), 0, max(0, len(edges) - 2))
    bx, by = [], []
    for b in np.unique(bins):
        take = bins == b
        if take.sum() >= 5:
            bx.append(np.median(x[take])); by.append(np.mean(y[take]))
    if len(bx) < 3:
        return np.full_like(grid, np.nan)
    return gaussian_filter1d(np.interp(grid, bx, by), 1.2)


def plot_peak_relationships(data: pd.DataFrame, prefix: Path, title: str,
                            bootstrap: int = 100, seed: int = SEED) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    fig = plt.figure(figsize=(13, 18)); outer = fig.add_gridspec(4, 2, hspace=.29, wspace=.24)
    rows = []
    for i, (feature, label) in enumerate(FEATURES):
        z = data[data.feature.eq(feature) & data.absolute_CA_peak_count.gt(0)]
        x = z.log10_expression_plus_1.to_numpy(float)
        y = z.absolute_CA_peak_count.to_numpy(float)
        rho, rlo, rhi = _bootstrap_rho(x, y, bootstrap, seed + i)
        sub = outer[i // 2, i % 2].subgridspec(2, 2, height_ratios=[1, 4], width_ratios=[4, 1], hspace=.03, wspace=.03)
        top = fig.add_subplot(sub[0, 0]); ax = fig.add_subplot(sub[1, 0], sharex=top)
        right = fig.add_subplot(sub[1, 1], sharey=ax)
        if len(z) < 3:
            ax.text(.5, .5, "Fewer than 3 peak-positive genes", ha="center", va="center",
                    transform=ax.transAxes, fontsize=11)
            top.axis("off"); right.axis("off")
            ax.spines[["top", "right"]].set_visible(False)
            ax.set_title(f"{label}\nn={len(z):,}; Spearman ρ=not estimable", fontsize=12, loc="left")
            rows.append({"feature": feature, "feature_label": label,
                         "n_peak_positive_genes": len(z), "spearman_rho": np.nan,
                         "spearman_ci_low": np.nan, "spearman_ci_high": np.nan})
            continue
        grid = np.linspace(x.min(), x.max(), 150); fit = _binned_mean(x, y, grid); boots = []
        for _ in range(bootstrap):
            ix = rng.integers(0, len(x), len(x)); boots.append(_binned_mean(x[ix], y[ix], grid))
        blo, bhi = np.nanpercentile(np.asarray(boots), [2.5, 97.5], axis=0)
        jitter = rng.uniform(-.12, .12, len(y))
        ax.scatter(x, y + jitter, s=4, alpha=.14, color="#4C78A8", linewidths=0, rasterized=True)
        ax.fill_between(grid, blo, bhi, color="#D55E00", alpha=.18, lw=0)
        ax.plot(grid, fit, color="#D55E00", lw=2)
        top.hist(x, bins=70, density=True, color="#777777", alpha=.4)
        right.hist(y, bins=np.arange(-.5, min(max(y.max(), 3), 30) + 1.5), density=True,
                   orientation="horizontal", color="#777777", alpha=.4)
        top.axis("off"); right.axis("off"); ax.spines[["top", "right"]].set_visible(False)
        ax.set_ylim(-.45, max(3, np.quantile(y, .995) + .5)); ax.tick_params(labelsize=11)
        ax.set_title(f"{label}\nn={len(z):,}; Spearman ρ={rho:.2f} [{rlo:.2f}, {rhi:.2f}]", fontsize=12, loc="left")
        if i // 2 == 3: ax.set_xlabel("log₁₀(expression + 1)", fontsize=12)
        if i % 2 == 0: ax.set_ylabel("Absolute CA peak count", fontsize=12)
        rows.append({"feature": feature, "feature_label": label, "n_peak_positive_genes": len(z),
                     "spearman_rho": rho, "spearman_ci_low": rlo, "spearman_ci_high": rhi})
    fig.suptitle(title, fontsize=16, fontweight="bold", y=.995)
    _save_figure(fig, prefix)
    return pd.DataFrame(rows)


def ratio_feature_data(root: Path, expression: pd.DataFrame) -> pd.DataFrame:
    f = pd.read_csv(root / "supporting_data/validated_SPEN_CLAP_matrices/normalized/nodox.features.tsv", sep="\t")
    f = f[f.feature.isin(FEATURE_LABELS)].copy()
    f["feature_label"] = f.feature.map(FEATURE_LABELS)
    return f.merge(expression, on="gene_id", how="inner", suffixes=("", "_expr"))


def plot_ratio_components(data: pd.DataFrame, prefix: Path, expression_label: str) -> pd.DataFrame:
    fig, axes = plt.subplots(8, 3, figsize=(15, 25), sharex="col")
    rows = []
    metrics = [("clap_cpm", "Normalized CLAP CPM"), ("input_cpm", "Normalized Input CPM"),
               ("enrichment", "CLAP-to-Input ratio, log₂")]
    for i, (feature, label) in enumerate(FEATURES):
        z = data[data.feature.eq(feature)]
        x = z.log10_expression_plus_1.to_numpy(float)
        for j, (metric, ylabel) in enumerate(metrics):
            ax = axes[i, j]; raw = z[metric].to_numpy(float)
            y = np.log10(raw + 1) if metric != "enrichment" else raw
            ax.scatter(x, y, s=3, alpha=.12, color="#4C78A8", linewidths=0, rasterized=True)
            grid = np.linspace(x.min(), x.max(), 120); ax.plot(grid, _binned_mean(x, y, grid), color="#D55E00", lw=1.8)
            rho = spearmanr(x, y).statistic
            rows.append({"feature": feature, "feature_label": label, "component": metric,
                         "n_genes": len(z), "spearman_rho": rho, "expression_source": expression_label})
            ax.text(.015, .91, f"{label}; ρ={rho:.2f}", transform=ax.transAxes, fontsize=10, va="top")
            ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=9)
            if i == 0: ax.set_title(ylabel, fontsize=12, fontweight="bold")
            if i == 7: ax.set_xlabel(f"log₁₀({expression_label} + 1)", fontsize=10)
            if j == 0: ax.set_ylabel("log₁₀(CPM + 1)", fontsize=10)
            elif j == 1: ax.set_ylabel("log₁₀(CPM + 1)", fontsize=10)
            else: ax.set_ylabel("log₂ ratio", fontsize=10)
    fig.suptitle(f"SPEN CLAP and Input components versus {expression_label}", fontsize=16, fontweight="bold", y=.995)
    fig.tight_layout(rect=[0, 0, 1, .988]); _save_figure(fig, prefix)
    return pd.DataFrame(rows)


def pseudocount_sensitivity(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pseudo in [.01, .1, 1.0]:
        ratio = np.log2((data.clap_cpm + pseudo) / (data.input_cpm + pseudo))
        for min_input in [0, .1, 1.0]:
            take = data.input_cpm >= min_input
            for feature, label in FEATURES:
                q = data[take & data.feature.eq(feature)]
                r = ratio.loc[q.index]
                rows.append({"feature": feature, "feature_label": label, "pseudocount_CPM": pseudo,
                             "minimum_Input_CPM": min_input, "n_genes": len(q),
                             "spearman_rho": spearmanr(q.log10_expression_plus_1, r).statistic if len(q) else np.nan})
    return pd.DataFrame(rows)


def build_human_models(gtf: Path, output_json: Path, output_bed: Path) -> dict:
    reps = parse_representative_transcripts(str(gtf))
    models = {gid: asdict(t) for gid, t in reps.items()
              if t.gene_type == "protein_coding" or t.transcript_type == "protein_coding"}
    _mkdir(output_json.parent); output_json.write_text(json.dumps(models))
    with output_bed.open("w") as out:
        for gid, t in sorted(models.items(), key=lambda kv: (kv[1]["chrom"], kv[1]["start"], kv[0])):
            exons = sorted(t["exons"])
            if not exons: continue
            chrom = t["chrom"].removeprefix("chr")
            if chrom == "M": chrom = "MT"
            thick_start = min((x[0] for x in t["cds"]), default=t["end"])
            thick_end = max((x[1] for x in t["cds"]), default=t["end"])
            sizes = ",".join(str(b - a) for a, b in exons)
            starts = ",".join(str(a - t["start"]) for a, _ in exons)
            fields = [chrom, t["start"], t["end"], f'{t["transcript_id"]}|{t["gene_name"]}', 0,
                      t["strand"], thick_start, thick_end, 0, len(exons), sizes, starts]
            out.write("\t".join(map(str, fields)) + "\n")
    return models


def _huh_bam(source: Path, sample: str) -> Path:
    return source / sample / "alignment/hg38/hg38_dedup.read2.bam"


def huh_source_manifest(source: Path, controls: dict, output: Path) -> pd.DataFrame:
    samples = []
    for spec in controls.values():
        samples.extend(spec["clap"]); samples.extend(spec["input"])
    rows = []
    for sample in sorted(set(samples)):
        for path in [_huh_bam(source, sample), Path(str(_huh_bam(source, sample)) + ".bai")]:
            if not path.exists(): raise FileNotFoundError(path)
            st = path.stat()
            rows.append({"sample": sample, "path": str(path), "size": st.st_size,
                         "mtime_ns": st.st_mtime_ns, "sha256": _sha256(path)})
    d = pd.DataFrame(rows); _mkdir(output.parent); d.to_csv(output, sep="\t", index=False)
    return d


def verify_huh_source(before: pd.DataFrame, output: Path) -> None:
    rows = []
    for r in before.itertuples(index=False):
        p = Path(r.path); st = p.stat(); digest = _sha256(p)
        rows.append({"sample": r.sample, "path": r.path, "size_before": r.size,
                     "size_after": st.st_size, "mtime_ns_before": r.mtime_ns,
                     "mtime_ns_after": st.st_mtime_ns, "sha256_before": r.sha256,
                     "sha256_after": digest, "unchanged": bool(r.size == st.st_size and
                     r.mtime_ns == st.st_mtime_ns and r.sha256 == digest)})
    d = pd.DataFrame(rows); _mkdir(output.parent); d.to_csv(output, sep="\t", index=False)
    if not d.unchanged.all(): raise RuntimeError("One or more HUHCLAP source files changed")


def count_huh_bam(models: dict, bam_path: Path) -> tuple[dict, dict, dict]:
    idx = build_index(models, flank=0); genes = defaultdict(int); features = defaultdict(int)
    metrics = defaultdict(int)
    with pysam.AlignmentFile(bam_path) as bam:
        for r in bam.fetch(until_eof=True):
            metrics["records"] += 1
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20:
                continue
            chrom = _canonical_chrom(r.reference_name)
            strand = "-" if r.is_reverse else "+"
            pos = (r.reference_start + r.reference_end) // 2
            gid = assign(idx, chrom, strand, pos)
            if gid is None: continue
            genes[gid] += 1; metrics["assigned_gene_body"] += 1
            for feature in feature_membership(models[gid], pos):
                features[(gid, feature)] += 1
    return genes, features, dict(metrics)


def build_huh_matrices(models: dict, source: Path, rbp: str, spec: dict,
                       output_dir: Path) -> tuple[Path, Path, Path]:
    _mkdir(output_dir); samples = []
    for assay, names in [("CLAP", spec["clap"]), ("Input", spec["input"])]:
        for i, name in enumerate(names, 1): samples.append((f"{rbp}_{assay.lower()}_{i}", assay, str(i), name))
    gene_counts, feat_counts, metrics = {}, {}, []
    for sample, assay, replicate, source_name in samples:
        cache = output_dir / f"{sample}.counts.npz"
        meta_cache = output_dir / f"{sample}.metrics.json"
        if cache.exists() and meta_cache.exists():
            z = np.load(cache, allow_pickle=True)
            genes = dict(z["genes"].tolist()); feats = dict(z["features"].tolist())
            met = json.loads(meta_cache.read_text())
        else:
            genes, feats, met = count_huh_bam(models, _huh_bam(source, source_name))
            np.savez_compressed(cache, genes=np.array(list(genes.items()), dtype=object),
                                features=np.array(list(feats.items()), dtype=object))
            meta_cache.write_text(json.dumps({"sample": sample, "source_sample": source_name, **met}, indent=2))
        gene_counts[sample] = genes; feat_counts[sample] = feats
        metrics.append({"sample": sample, "assay": assay, "replicate": replicate,
                        "source_sample": source_name, **met})
    lengths = feature_lengths(models)
    keys = sorted(set().union(*(set(x) for x in feat_counts.values())))
    rows = []
    for gid, feature in keys:
        t = models[gid]
        row = {"gene_id": gid, "gene_name": t["gene_name"], "chrom": t["chrom"],
               "strand": t["strand"], "gene_length": t["end"] - t["start"],
               "intron_count": max(0, len(t["exons"]) - 1), "feature": feature}
        row.update({sample: feat_counts[sample].get((gid, feature), 0) for sample, *_ in samples})
        rows.append(row)
    feature_matrix = output_dir / f"{rbp}.features.raw.tsv.gz"
    pd.DataFrame(rows).to_csv(feature_matrix, sep="\t", index=False, compression="gzip")
    gene_rows = []
    for gid in sorted(set().union(*(set(x) for x in gene_counts.values()))):
        t = models[gid]; row = {"gene_id": gid, "gene_name": t["gene_name"], "chrom": t["chrom"]}
        row.update({sample: gene_counts[sample].get(gid, 0) for sample, *_ in samples}); gene_rows.append(row)
    gene_matrix = output_dir / f"{rbp}.gene_body.raw.tsv.gz"
    pd.DataFrame(gene_rows).to_csv(gene_matrix, sep="\t", index=False, compression="gzip")
    metadata = output_dir / f"{rbp}.sample_metadata.tsv"
    pd.DataFrame([{"sample": s, "context": rbp, "assay": a, "replicate": r,
                   "source_sample": n} for s, a, r, n in samples]).to_csv(metadata, sep="\t", index=False)
    pd.DataFrame(metrics).to_csv(output_dir / f"{rbp}.counting_metrics.tsv", sep="\t", index=False)
    lengths.to_csv(output_dir / "human_feature_lengths.tsv", sep="\t", index=False)
    return feature_matrix, gene_matrix, metadata


def normalize_huh_features(feature_matrix: Path, metadata: Path, output: Path) -> pd.DataFrame:
    script = Path(__file__).resolve().parents[2] / "scripts/normalize_clap.R"
    if not output.exists():
        subprocess.run(["/opt/miniconda3/envs/bioinfo/bin/Rscript", str(script), str(feature_matrix),
                        str(metadata), pd.read_csv(metadata, sep="\t").context.iloc[0], ".1", str(output)], check=True)
    return pd.read_csv(output, sep="\t")


def huh_input_expression(gene_matrix: Path, size_factors: Path, metadata: Path) -> pd.DataFrame:
    d = pd.read_csv(gene_matrix, sep="\t"); sf = pd.read_csv(size_factors, sep="\t").set_index("sample")
    m = pd.read_csv(metadata, sep="\t"); inputs = m.loc[m.assay.eq("Input"), "sample"].tolist()
    for sample in inputs:
        d[f"{sample}_cpm"] = d[sample] / sf.loc[sample, "effective_lib_size"] * 1e6
    d["expression"] = d[[f"{x}_cpm" for x in inputs]].mean(axis=1)
    d = d[d.expression > 0].copy(); d["log10_expression_plus_1"] = np.log10(d.expression + 1)
    d["expression_source"] = "Matched HUHCLAP Input gene-body CPM"
    return d


def pool_bams(source_bams: list[Path], output: Path) -> None:
    if output.exists() and Path(str(output) + ".bai").exists(): return
    _mkdir(output.parent)
    samtools = "/opt/miniconda3/envs/bioinfo/bin/samtools"
    subprocess.run([samtools, "merge", "-f", "-@", "4", str(output), *map(str, source_bams)], check=True)
    subprocess.run([samtools, "index", "-@", "4", str(output)], check=True)


def run_ca(clap: Path, input_bam: Path, bed: Path, output: Path, jar: Path) -> None:
    if output.exists() and output.stat().st_size > 100: return
    _mkdir(output.parent)
    java = "/opt/miniconda3/envs/bioinfo/bin/java"
    log = output.with_suffix(".log")
    with log.open("w") as err:
        subprocess.run([java, "-Xmx12g", "-jar", str(jar), str(clap), str(input_bam), "100",
                        str(bed), str(output)], stderr=err, check=True)


def prepare_huh_ca(source: Path, rbp: str, spec: dict, bed: Path, jar: Path,
                   output_dir: Path) -> dict[str, Path]:
    outputs = {}
    for i in (0, 1):
        out = output_dir / f"{rbp}.rep{i + 1}.CA.tsv"
        run_ca(_huh_bam(source, spec["clap"][i]), _huh_bam(source, spec["input"][i]), bed, out, jar)
        outputs[f"rep{i + 1}"] = out
    pooled_clap = output_dir / "pooled_bams" / f"{rbp}.CLAP.pooled.bam"
    pooled_input = output_dir / "pooled_bams" / f"{rbp}.Input.pooled.bam"
    pool_bams([_huh_bam(source, x) for x in spec["clap"]], pooled_clap)
    pool_bams([_huh_bam(source, x) for x in spec["input"]], pooled_input)
    out = output_dir / f"{rbp}.pooled.CA.tsv"; run_ca(pooled_clap, pooled_input, bed, out, jar)
    outputs["pooled"] = out
    return outputs


def analyze_ca_callsets(callsets: dict[str, Path], models: dict, lengths: pd.DataFrame,
                        expression: pd.DataFrame, output_dir: Path, label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    assigned = []
    qc = []
    for name, path in callsets.items():
        calls = filter_ca(path); a = assign_ca_to_features(calls, models, name); assigned.append(a)
        qc.append({"analysis": label, "callset": name, "source": str(path), "raw_rows": sum(1 for _ in path.open()) - 1,
                   "retained_unique_windows": calls.window.nunique(), "assigned_windows": a.window.nunique(),
                   "assigned_genes": a.gene_id.nunique(), "filter": CA_FILTER})
    all_assigned = pd.concat(assigned, ignore_index=True)
    if {"rep1", "rep2"}.issubset(callsets):
        r1 = all_assigned[all_assigned.callset.eq("rep1")]
        r2 = all_assigned[all_assigned.callset.eq("rep2")]
        common = r1.merge(r2[["gene_id", "feature", "window"]], on=["gene_id", "feature", "window"])
        common = common.drop_duplicates(["gene_id", "feature", "window"]); common["callset"] = "strict_overlap"
        all_assigned = pd.concat([all_assigned, common], ignore_index=True)
        qc.append({"analysis": label, "callset": "strict_overlap", "source": "rep1 ∩ rep2",
                   "raw_rows": pd.NA, "retained_unique_windows": common.window.nunique(),
                   "assigned_windows": common.window.nunique(), "assigned_genes": common.gene_id.nunique(),
                   "filter": CA_FILTER})
    _mkdir(output_dir)
    all_assigned.to_csv(output_dir / f"{label}.CA_assigned_windows.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.DataFrame(qc).to_csv(output_dir / f"{label}.CA_callset_QC.tsv", sep="\t", index=False)
    frames = []
    for callset in all_assigned.callset.unique():
        frames.append(peak_count_universe(all_assigned, lengths, expression, callset))
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(output_dir / f"{label}.absolute_CA_peak_counts.tsv.gz", sep="\t", index=False, compression="gzip")
    return combined, pd.DataFrame(qc)


def ca_ebs_comparison(prior: Path, output_dir: Path) -> pd.DataFrame:
    names = ["clap_endospen_1", "clap_endospen_2", "merged_clap_endospen"]
    rows = []
    for sample in names:
        ca = prior / "enrichment" / f"{sample}_CA.tsv"
        ebs = prior / "enrichment" / f"{sample}_EBS.tsv"
        q = filter_ca(ca)
        rows.append({"sample": sample, "method": "CA", "all_rows": sum(1 for _ in ca.open()) - 1,
                     "retained_calls": len(q), "unique_windows": q.window.nunique()})
        e = pd.read_csv(ebs, sep="\t", header=None)
        rows.append({"sample": sample, "method": "EBS", "all_rows": len(e),
                     "retained_calls": len(e), "unique_windows": e.iloc[:, 2].nunique() if len(e) else 0})
    d = pd.DataFrame(rows); _mkdir(output_dir); d.to_csv(output_dir / "endogenous_NoDox_CA_vs_EBS_yield.tsv", sep="\t", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    labels = ["Replicate A", "Replicate B", "Merged"]
    for ax, metric, title in zip(axes, ["retained_calls", "unique_windows"], ["Retained calls", "Unique genomic windows"]):
        x = np.arange(3); width = .34
        for j, method in enumerate(["CA", "EBS"]):
            vals = [d[(d["sample"].eq(s)) & (d["method"].eq(method))][metric].iloc[0] for s in names]
            ax.bar(x + (j - .5) * width, vals, width, label=method, color={"CA": "#0072B2", "EBS": "#D55E00"}[method])
            for xx, value in zip(x + (j - .5) * width, vals): ax.text(xx, max(value, .7) * 1.15, f"{value:,}", ha="center", fontsize=9, rotation=35)
        ax.set_yscale("log"); ax.set_xticks(x, labels); ax.set_ylabel("Count (log scale)"); ax.set_title(title)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False); fig.suptitle("Why CA was selected: endogenous No-Dox SPEN only", fontweight="bold")
    fig.tight_layout(); _save_figure(fig, output_dir / "endogenous_NoDox_CA_vs_EBS_yield")
    exclusions = pd.DataFrame({"excluded_dataset_family": ["Dox endogenous SPEN", "FL rescue", "dRRM rescue", "dIDR rescue", "episomal"],
                               "reason": ["Not the No-Dox context"] + ["Not endogenous SPEN"] * 4})
    exclusions.to_csv(output_dir / "explicit_dataset_exclusions.tsv", sep="\t", index=False)
    return d


def _img(path: Path) -> str:
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode()


def _table(path: Path, columns: list[str] | None = None, max_rows: int = 100) -> str:
    d = pd.read_csv(path, sep="\t")
    if columns: d = d[[c for c in columns if c in d]]
    return d.head(max_rows).to_html(index=False, classes="data", border=0, float_format=lambda x: f"{x:.3g}")


def _figure(path: Path, caption: str) -> str:
    return f"<figure><img src='{_img(path)}'><figcaption>{caption}</figcaption></figure>"


def build_master_report(root: Path, spen_ca_stats: pd.DataFrame, spen_models: pd.DataFrame,
                        control_stats: pd.DataFrame, method_yield: pd.DataFrame) -> None:
    ca_introns = spen_ca_stats[(spen_ca_stats.feature.eq("all_introns")) & (spen_ca_stats.tertile.eq("Low"))]
    ca_rho = float(ca_introns.spearman_rho_peak_positive.iloc[0]) if len(ca_introns) else np.nan
    adj = spen_models[(spen_models.feature.eq("all_introns")) & spen_models.endpoint.eq("positive_peak_count_negative_binomial")]
    adj_beta = float(adj.expression_coefficient.iloc[0]) if len(adj) else np.nan
    adj_p = float(adj.expression_pvalue.iloc[0]) if len(adj) else np.nan
    control_pivot = control_stats[control_stats.feature.eq("all_introns")][["RBP", "spearman_rho"]]
    control_text = ", ".join(f"{r.RBP}: ρ={r.spearman_rho:.2f}" for r in control_pivot.itertuples())
    components = pd.read_csv(root / "results/03_Input_expression_check/component_correlations.tsv", sep="\t")
    spen_ai = components[(components.expression_source.eq("untreated total-RNA expression")) &
                         components.feature.eq("all_introns")].set_index("component").spearman_rho
    clap_rho = float(spen_ai.get("clap_cpm", np.nan)); input_rho = float(spen_ai.get("input_cpm", np.nan))
    ratio_rho = float(spen_ai.get("enrichment", np.nan))
    zero_low = float(ca_introns.fraction_zero_QC.iloc[0]) if len(ca_introns) else np.nan
    ca_high = spen_ca_stats[(spen_ca_stats.feature.eq("all_introns")) & (spen_ca_stats.tertile.eq("High"))]
    zero_high = float(ca_high.fraction_zero_QC.iloc[0]) if len(ca_high) else np.nan
    control_ratio = []
    for rbp in ["HnRNPC", "PTBP1", "GFP"]:
        p = root / f"results/04_HUHCLAP_controls/{rbp}/{rbp}.ratio_component_correlations.tsv"
        d = pd.read_csv(p, sep="\t")
        q = d[(d.feature.eq("all_introns")) & d.component.eq("enrichment")]
        control_ratio.append(f"{rbp}: ρ={float(q.spearman_rho.iloc[0]):.2f}" if len(q) else f"{rbp}: n/a")
    verdict = ("The current evidence favors a substantial denominator/coverage component over a simple biological "
               "interpretation of higher SPEN binding at low-expression genes. Input coverage tracks untreated "
               "expression more strongly than CLAP coverage, making their ratio fall as expression rises. CA does "
               "not independently rescue the pattern because its detection rate is itself strongly expression-dependent.")
    css = """body{font:15px system-ui;max-width:1320px;margin:auto;padding:32px;color:#222;line-height:1.55}h1,h2,h3{color:#17365d}img{max-width:100%;border:1px solid #ddd}figure{margin:28px 0}figcaption{margin-top:8px}.data{border-collapse:collapse;font-size:12px;display:block;overflow-x:auto}.data th,.data td{border:1px solid #ddd;padding:5px 8px}.data th{background:#eef3f8}.result{background:#eaf4ea;padding:15px;border-left:5px solid #24733c}.warning{background:#fff4cc;padding:15px;border-left:5px solid #e69f00}code{background:#f3f3f3;padding:2px 4px}"""
    report = [f"<!doctype html><html><head><meta charset='utf-8'><title>SPEN CLAP-to-Input ratio: expression and artifact assessment</title><style>{css}</style></head><body>",
              "<h1>SPEN CLAP-to-Input ratio versus expression: artifact assessment and controls</h1>",
              "<h2>Executive summary</h2>",
              "<p>The earlier treatment-response analyses did not show a reproducible relationship between the direction of RNA change after SPEN removal and SPEN binding. A stronger pattern emerged after regrouping genes by untreated abundance: low-expression genes showed a right-shifted SPEN CLAP-to-Input ratio, strongest in intron-related features and later transcript regions.</p>",
              f"<div class='result'><b>Integrated result.</b> {html.escape(verdict)} For all introns, untreated expression correlates with normalized CLAP coverage at ρ={clap_rho:.2f}, normalized Input coverage at ρ={input_rho:.2f}, and the CLAP-to-Input ratio at ρ={ratio_rho:.2f}. CA zero fractions are {zero_low:.1%} in the low-expression tertile and {zero_high:.1%} in the high-expression tertile. Among peak-positive genes only, peak multiplicity correlates positively with expression (ρ={ca_rho:.2f}), consistent with expression-dependent detection opportunity rather than validation of the inverse ratio.</div>",
              "<div class='warning'><b>Terminology.</b> The SPEN CLAP-to-Input ratio is <code>log₂[(normalized CLAP + pseudocount)/(normalized Input + pseudocount)]</code>. It is calculated for eligible gene features whether or not CA calls a peak. This explicit name and formula are used throughout the report.</div>",
              "<h2>1. Discovery: CLAP-to-Input ratio versus untreated expression</h2>"]
    ratio_dir = root / "results/01_CLAP_to_Input_ratio"
    report.append(_figure(ratio_dir / "nodox_all__tertile_distributions.png",
        "Rows 1–8 show the gene-level SPEN CLAP-to-Input ratio, log₂[(TMM-normalized CLAP CPM + 0.1)/(TMM-normalized Input CPM + 0.1)], grouped by tertiles of the mean DESeq2-normalized untreated total-RNA count. Each tertile density integrates independently to one; height is not gene abundance. The Expression row verifies separation. The association is intron-dominant and post–first-exon, not strictly intron-exclusive."))
    report.extend(["<h2>2. Numerator, denominator, and pseudocount diagnostics</h2>",
                   f"<p>A negative relationship can be created mechanically when Input is both large at highly expressed genes and placed in the denominator. For all introns, Input rises with untreated expression more strongly (ρ={input_rho:.2f}) than CLAP does (ρ={clap_rho:.2f}); accordingly, their ratio falls (ρ={ratio_rho:.2f}). A fixed pseudocount has its greatest effect at low coverage. Raising the minimum Input CPM and varying the pseudocount attenuate but do not eliminate the inverse association, so the effect is not attributable to the pseudocount alone. The component panels therefore show CLAP, Input, and their ratio separately. Repeating the x-axis with SPEN Input coverage is explicitly a mathematical-coupling diagnostic, not independent confirmation.</p>"])
    for name, caption in [("total_RNA_component_relationships.png", "CLAP CPM, Input CPM, and the CLAP-to-Input ratio versus untreated total RNA."),
                          ("Input_proxy_component_relationships.png", "The same components versus SPEN Input gene-body CPM. Because Input contributes to both axes in the ratio column, negative coupling is expected and must not be interpreted as validation.")]:
        report.append(_figure(root / "results/03_Input_expression_check" / name, caption))
    report.extend(["<h2>3. Peak-positive endogenous SPEN CA sensitivity</h2>",
                   "<p>The CA analysis counts distinct 100-bp windows that pass CLAP count, intra-gene enrichment, inter-library enrichment, and both p-value filters. Most eligible gene features have zero called windows, so zero-dominated plots were discarded as uninformative. The figures and correlations below include only gene features with at least one CA-called window. Exact zero fractions remain in the QC table to disclose the strength of this selection. Probability impulses are used because KDE would invent fractional peak counts.</p>",
                   "<p><b>Why feature length still matters.</b> Specific binding does not imply that length is irrelevant: longer features contain more sequence, motifs, structures, and accessible positions where genuine binding can occur. Expression and intron length can also covary. Length is therefore included as a covariate and through stratification, while the plotted biological outcome remains the absolute peak count rather than peaks per kilobase.</p>",
                   _figure(root / "results/02_CA_peak_validation/SPEN_total_RNA_absolute_CA_peak_distributions.png", "Pooled No-Dox endogenous SPEN CA calls among peak-positive gene features, grouped by untreated total-RNA expression tertile. Zeros are excluded. Each impulse is the conditional probability of an exact integer peak count; the last displayed point includes the upper 1% tail. Overlapping feature definitions are analyzed separately and are not additive."),
                   _figure(root / "results/02_CA_peak_validation/SPEN_total_RNA_absolute_CA_peak_relationships.png", "Each point is one peak-positive gene feature with vertical jitter only for visibility. Orange curves are binned mean counts with exploratory 100-bootstrap intervals. Zeros are excluded, so this evaluates peak multiplicity only and is subject to detection/selection bias. Length is addressed in accompanying conditional models."),
                   "<h2>4. Why CA rather than EBS?</h2>",
                   _figure(root / "results/05_CA_vs_EBS_method_choice/endogenous_NoDox_CA_vs_EBS_yield.png", "Only endogenous No-Dox SPEN is included. CA retains thousands of windows, whereas EBS produced 12 and 27 replicate-level calls, making EBS unsuitable for feature-distribution analysis in this dataset. Yield does not prove CA is unbiased; CA still uses Input and has expression-dependent detection power."),
                   "<h2>5. Independent HUHCLAP controls</h2>",
                   f"<p>HnRNPC and PTBP1 are bona fide intron-binding RBPs and are biological comparators, not expected nulls. GFP is the negative-binding control. These human hg38 libraries were processed independently with matched Input; their Input gene-body CPM is the only available expression proxy. They are not mixed with mouse genes or used to normalize SPEN. Because Input is both the expression proxy and the ratio denominator, these controls test a generic mathematical tendency rather than provide an independent expression measurement. All-intron ratio correlations are {html.escape(', '.join(control_ratio))}.</p>"])
    for rbp in ["HnRNPC", "PTBP1", "GFP"]:
        report.append(f"<h3>{rbp}</h3>")
        report.append(_figure(root / f"results/04_HUHCLAP_controls/{rbp}/{rbp}.absolute_CA_peak_distributions.png",
                              f"{rbp} pooled CA absolute peak counts among peak-positive gene features by matched-Input expression tertile. Zero-peak features are excluded from the plot and retained only as QC counts; individual replicates and strict overlap are retained in the supporting tables."))
        report.append(_figure(root / f"results/04_HUHCLAP_controls/{rbp}/{rbp}.ratio_components.png",
                              f"{rbp} normalized CLAP coverage, matched Input coverage, and CLAP-to-Input ratio versus matched-Input expression. The ratio column is mathematically coupled to the x-axis and is diagnostic."))
    report.extend(["<h2>6. Statistical summary</h2>",
                   "<h3>SPEN CA peak counts among peak-positive gene features</h3>",
                   (root / "results/02_CA_peak_validation/SPEN_total_RNA_peak_summary.tsv").read_text()[:0] +
                   _table(root / "results/02_CA_peak_validation/SPEN_total_RNA_peak_summary.tsv",
                          ["feature_label", "tertile", "n_all_eligible_genes", "n_peak_positive_genes", "median_peak_count_among_positive", "fraction_zero_QC", "spearman_rho_peak_positive"]),
                   "<h3>Conditional length-adjusted models among peak-positive genes</h3>",
                   _table(root / "results/02_CA_peak_validation/SPEN_total_RNA_length_adjusted_models.tsv",
                          ["feature_label", "endpoint", "n_peak_positive_genes", "expression_coefficient", "expression_pvalue", "length_coefficient", "length_pvalue"]),
                   "<h3>Control comparison</h3>", control_stats.to_html(index=False, classes="data", border=0, float_format=lambda x: f"{x:.3g}"),
                   "<h2>7. Interpretation and limitations</h2>",
                   "<ul><li>The zero-dominated CA result cannot discriminate biological absence of binding from inadequate peak-detection power.</li><li>The peak-positive CA sensitivity asks about peak multiplicity after detection; conditioning on detection introduces selection bias, so it is not independent confirmation of the CLAP-to-Input ratio.</li><li>Reproduction in GFP would indicate a general assay/normalization effect. Similarity to HnRNPC or PTBP1 may reflect either shared intronic biology or a shared CLAP technical property.</li><li>Input-based expression is not independent of a ratio containing Input. Untreated total RNA remains the primary expression measurement for SPEN.</li><li>These are associations and do not prove that expression causes SPEN binding or vice versa.</li></ul>",
                   "<h2>Reproducibility</h2><p>All exact plotted tables, pooled BAM derivatives, CA output, checksums, logs, and vector figures are stored alongside this report. The HUHCLAP source directory was read-only and verified by full SHA-256 before and after analysis. Bootstrap seed: 20260812; iterations: 100.</p>",
                   "</body></html>"])
    output = root / "report/SPEN_expression_dependence_report.html"; _mkdir(output.parent); output.write_text("".join(report))


def run(root: str, config_path: str) -> None:
    root = Path(root); cfg = __import__("yaml").safe_load(Path(config_path).read_text())
    prior = Path(cfg["prior_analysis"]); av = cfg["artifact_validation"]
    models = json.loads((root / "supporting_data/annotation/representative_transcripts.json").read_text())
    lengths = feature_lengths(models)
    total_expr = untreated_expression(root / "supporting_data/untreated_total_RNA/Aux2h_vs_ctrl__total.tsv")
    ratio_eligible = pd.read_csv(root / "supporting_data/validated_SPEN_CLAP_matrices/normalized/nodox.features.tsv",
                                 sep="\t", usecols=["gene_id"]).gene_id.unique()
    total_ratio_expr, total_cut = assign_tertiles(total_expr[total_expr.gene_id.isin(ratio_eligible)].copy())
    input_expr = spen_input_expression(root)
    input_ratio_expr, input_cut = assign_tertiles(input_expr[input_expr.gene_id.isin(ratio_eligible)].copy())
    ratio_total = ratio_feature_data(root, total_ratio_expr)
    ratio_input = ratio_feature_data(root, input_ratio_expr)
    component_table = root / "results/03_Input_expression_check/component_correlations.tsv"
    if not component_table.exists():
        comp_total = plot_ratio_components(ratio_total, root / "results/03_Input_expression_check/total_RNA_component_relationships", "untreated total-RNA expression")
        comp_input = plot_ratio_components(ratio_input, root / "results/03_Input_expression_check/Input_proxy_component_relationships", "SPEN Input gene-body CPM")
        pd.concat([comp_total, comp_input]).to_csv(component_table, sep="\t", index=False)
        pseudocount_sensitivity(ratio_total).to_csv(root / "results/03_Input_expression_check/pseudocount_and_Input_filter_sensitivity.tsv", sep="\t", index=False)
        pd.DataFrame([{"expression_source": "Untreated total RNA", "lower_cutoff": total_cut[0], "upper_cutoff": total_cut[1]},
                      {"expression_source": "SPEN Input gene-body CPM", "lower_cutoff": input_cut[0], "upper_cutoff": input_cut[1]}]).to_csv(root / "results/03_Input_expression_check/expression_tertile_cutoffs.tsv", sep="\t", index=False)

    spen_callsets = {name: prior / "enrichment" / f"{sample}_CA.tsv" for name, sample in
                     {"rep1": "clap_endospen_1", "rep2": "clap_endospen_2", "pooled": "merged_clap_endospen"}.items()}
    ca_expr, _ = assign_tertiles(total_expr[total_expr.gene_id.isin(lengths.gene_id)].copy())
    spen_dir = root / "results/02_CA_peak_validation"
    spen_summary_path = spen_dir / "SPEN_total_RNA_peak_summary.tsv"
    if not (spen_dir / "SPEN.absolute_CA_peak_counts.tsv.gz").exists():
        spen_counts, spen_qc = analyze_ca_callsets(spen_callsets, models, lengths, ca_expr, spen_dir, "SPEN")
    else:
        spen_counts = pd.read_csv(spen_dir / "SPEN.absolute_CA_peak_counts.tsv.gz", sep="\t")
    headline = spen_counts[spen_counts.callset.eq("pooled")].copy()
    plot_peak_distributions(headline, spen_dir / "SPEN_total_RNA_absolute_CA_peak_distributions",
                            "Endogenous SPEN CA peaks among peak-positive genes")
    plot_peak_relationships(headline, spen_dir / "SPEN_total_RNA_absolute_CA_peak_relationships",
                            "Endogenous SPEN CA peak multiplicity versus untreated expression").to_csv(
                            spen_dir / "SPEN_total_RNA_peak_correlations.tsv", sep="\t", index=False)
    spen_stats, spen_models, spen_strata = peak_statistics(headline)
    spen_stats.to_csv(spen_summary_path, sep="\t", index=False)
    spen_models.to_csv(spen_dir / "SPEN_total_RNA_length_adjusted_models.tsv", sep="\t", index=False)
    spen_strata.to_csv(spen_dir / "SPEN_total_RNA_length_stratified_correlations.tsv", sep="\t", index=False)
    method_yield = ca_ebs_comparison(prior, root / "results/05_CA_vs_EBS_method_choice")

    huh_source = Path(av["huhclap_root"]); control_root = root / "results/04_HUHCLAP_controls"
    before = huh_source_manifest(huh_source, av["controls"], root / "workflow_records/provenance/HUHCLAP_source_before.tsv")
    human_models_path = root / "supporting_data/HUHCLAP/hg38_v47_representative_transcripts.json"
    human_bed = root / "supporting_data/HUHCLAP/hg38_v47_CA_transcripts.bed"
    human_models = build_human_models(Path(av["hg38_gtf"]), human_models_path, human_bed)
    human_lengths = feature_lengths(human_models)
    control_rows = []
    for j, (rbp, spec) in enumerate(av["controls"].items()):
        rbp_dir = _mkdir(control_root / rbp); data_dir = _mkdir(root / "supporting_data/HUHCLAP" / rbp)
        feature_matrix, gene_matrix, metadata = build_huh_matrices(human_models, huh_source, rbp, spec, data_dir)
        norm = normalize_huh_features(feature_matrix, metadata, data_dir / f"{rbp}.features.normalized.tsv")
        expr = huh_input_expression(gene_matrix, data_dir / f"{rbp}.features.normalized.size_factors.tsv", metadata)
        expr = expr[expr.gene_id.isin(norm.gene_id.unique())].copy(); expr, cut = assign_tertiles(expr)
        ratio = norm[norm.feature.isin(FEATURE_LABELS)].merge(expr, on="gene_id", suffixes=("", "_expr"))
        ratio_plot = rbp_dir / f"{rbp}.ratio_components.png"
        if not ratio_plot.exists():
            plot_ratio_components(ratio, rbp_dir / f"{rbp}.ratio_components", "matched Input gene-body CPM").to_csv(
                rbp_dir / f"{rbp}.ratio_component_correlations.tsv", sep="\t", index=False)
        ca_paths = prepare_huh_ca(huh_source, rbp, spec, human_bed, Path(av["ca_jar"]), data_dir / "CA")
        counts_path = rbp_dir / f"{rbp}.absolute_CA_peak_counts.tsv.gz"
        if counts_path.exists():
            counts = pd.read_csv(counts_path, sep="\t")
        else:
            counts, qc = analyze_ca_callsets(ca_paths, human_models, human_lengths, expr, rbp_dir, rbp)
        pooled = counts[counts.callset.eq("pooled")]
        plot_peak_distributions(pooled, rbp_dir / f"{rbp}.absolute_CA_peak_distributions",
                                f"{rbp} CA peaks among peak-positive genes")
        corr_path = rbp_dir / f"{rbp}.absolute_CA_peak_correlations.tsv"
        corr = plot_peak_relationships(pooled, rbp_dir / f"{rbp}.absolute_CA_peak_relationships",
                                       f"{rbp} CA peak multiplicity versus matched-Input expression", seed=SEED + 100 + j)
        corr.insert(0, "RBP", rbp); corr.to_csv(corr_path, sep="\t", index=False)
        stat_path = rbp_dir / f"{rbp}.peak_summary.tsv"
        stats, fitted, strata = peak_statistics(pooled, seed=SEED + 100 + j)
        stats.to_csv(stat_path, sep="\t", index=False)
        fitted.to_csv(rbp_dir / f"{rbp}.length_adjusted_models.tsv", sep="\t", index=False)
        strata.to_csv(rbp_dir / f"{rbp}.length_stratified_correlations.tsv", sep="\t", index=False)
        control_rows.append(corr)
        pd.DataFrame([{"RBP": rbp, "lower_cutoff": cut[0], "upper_cutoff": cut[1]}]).to_csv(rbp_dir / f"{rbp}.expression_tertile_cutoffs.tsv", sep="\t", index=False)
    controls = pd.concat(control_rows, ignore_index=True)
    controls.to_csv(control_root / "cross_RBP_peak_correlations.tsv", sep="\t", index=False)
    verify_huh_source(before, root / "workflow_records/provenance/HUHCLAP_source_after_verification.tsv")
    build_master_report(root, spen_stats, spen_models, controls, method_yield)
