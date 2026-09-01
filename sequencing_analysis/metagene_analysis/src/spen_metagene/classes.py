from __future__ import annotations

import numpy as np
import pandas as pd

CLASS_ORDER = [
    "significant_decrease",
    "nonsignificant_decrease",
    "no_change",
    "nonsignificant_elevation",
    "significant_elevation",
]

COLORS = {
    "significant_decrease": "#0072B2",
    "nonsignificant_decrease": "#56B4E9",
    "no_change": "#999999",
    "nonsignificant_elevation": "#E69F00",
    "significant_elevation": "#D55E00",
    "not_tested": "#D9D9D9",
}


def classify_frame(df: pd.DataFrame, fdr: float, delta: float, lfc_column: str = "log2FoldChange_mle") -> pd.Series:
    """Assign MLE-effect volcano bins; NA padj is explicitly not_tested."""
    lfc = pd.to_numeric(df[lfc_column], errors="coerce")
    padj = pd.to_numeric(df["padj"], errors="coerce")
    out = pd.Series("not_tested", index=df.index, dtype="object")
    tested = lfc.notna() & padj.notna()
    out.loc[tested & (lfc.abs() <= delta)] = "no_change"
    out.loc[tested & (lfc > delta) & (padj < fdr)] = "significant_elevation"
    out.loc[tested & (lfc > delta) & (padj >= fdr)] = "nonsignificant_elevation"
    out.loc[tested & (lfc < -delta) & (padj < fdr)] = "significant_decrease"
    out.loc[tested & (lfc < -delta) & (padj >= fdr)] = "nonsignificant_decrease"
    return out


def classify_file(input_tsv: str, output_tsv: str, delta: float, fdrs=(0.05, 0.10)) -> None:
    df = pd.read_csv(input_tsv, sep="\t")
    for fdr in fdrs:
        label = f"class_fdr{int(round(fdr * 100)):02d}"
        df[label] = classify_frame(df, fdr, delta, "log2FoldChange_mle")
        df[label + "_apeglm_sensitivity"] = classify_frame(df, fdr, delta, "log2FoldChange_apeglm")
    df["minus_log10_padj"] = -np.log10(pd.to_numeric(df["padj"], errors="coerce").clip(lower=np.nextafter(0, 1)))
    df.to_csv(output_tsv, sep="\t", index=False)


def early_agreement(one: pd.DataFrame, two: pd.DataFrame, fdr: float, delta: float) -> pd.DataFrame:
    cols = ["gene_id", "gene_name", "chrom", "log2FoldChange_mle", "padj", f"class_fdr{int(fdr*100):02d}"]
    a = one[[c for c in cols if c in one]].copy().add_suffix("_1h").rename(columns={"gene_id_1h": "gene_id"})
    b = two[[c for c in cols if c in two]].copy().add_suffix("_2h").rename(columns={"gene_id_2h": "gene_id"})
    x = a.merge(b, on="gene_id", how="outer")
    l1, l2 = x["log2FoldChange_mle_1h"], x["log2FoldChange_mle_2h"]
    tested = l1.notna() & l2.notna() & x["padj_1h"].notna() & x["padj_2h"].notna()
    x["agreement"] = "not_tested"
    x.loc[tested, "agreement"] = "transient_or_weak"
    x.loc[tested & (l1.abs() <= delta) & (l2.abs() <= delta), "agreement"] = "stable_no_change"
    x.loc[tested & (l1 > delta) & (l2 > delta), "agreement"] = "concordant_elevation"
    x.loc[tested & (l1 < -delta) & (l2 < -delta), "agreement"] = "concordant_decrease"
    x.loc[tested & (((l1 > delta) & (l2 < -delta)) | ((l1 < -delta) & (l2 > delta))), "agreement"] = "discordant_direction"
    s1 = x["padj_1h"] < fdr
    s2 = x["padj_2h"] < fdr
    x["significance_support"] = "neither"
    x.loc[s1 ^ s2, "significance_support"] = "one_timepoint"
    x.loc[s1 & s2, "significance_support"] = "both_timepoints"
    return x
