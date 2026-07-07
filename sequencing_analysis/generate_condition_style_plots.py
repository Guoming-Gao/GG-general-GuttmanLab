#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


BASE_DIR = Path("/Volumes/guttman/users/gmgao/Data_seq/20260629-AmpSeq-KOveri on polyclonal-cross-junction SNPs-aware on WT vs Dox")
FIG_DIR = BASE_DIR / "figures"
SUMMARY_CSV = BASE_DIR / "quantification" / "allele_quantification_summary.csv"

GDNA_ORDER = ["RNF12_Koveri_P2", "Tsix_Koveri_F1_R1", "Tsix_Koveri_F2_R1", "Tsix_Koveri_F3_R1"]
RNA_ORDER = ["esc_Xist_F1_R1", "esc_Xist_F2_R2", "esc_Kdm5c", "esc_Kdm6a", "si_Tsix", "si_RNF12_Rlim", "si_Rbmx"]
CONDITION_ORDER = {
    "gDNA_WT": GDNA_ORDER,
    "gDNA_polyclonal_dRNF12": GDNA_ORDER,
    "gDNA_polyclonal_dTsix": GDNA_ORDER,
    "RNA_WT_diff": RNA_ORDER,
    "RNA_Dox_24h": RNA_ORDER,
}
PRIMER_LABELS = {
    "RNF12_Koveri_P2": "RNF12_Koveri_P2\nF: RNF12_Koveri-P2-F / R: RNF12_Koveri-P2-R",
    "Tsix_Koveri_F1_R1": "Tsix_Koveri_F1_R1\nF: Tsix_Koveri-0616F1 / R: Tsix_Koveri-0616R1",
    "Tsix_Koveri_F2_R1": "Tsix_Koveri_F2_R1\nF: Tsix_Koveri-0616F2 / R: Tsix_Koveri-0616R1",
    "Tsix_Koveri_F3_R1": "Tsix_Koveri_F3_R1\nF: Tsix_Koveri-0616F3 / R: Tsix_Koveri-0616R1",
    "esc_Xist_F1_R1": "esc_Xist_F1_R1\nF: esc-Xist-0616F1 / R: esc-Xist-0616R1",
    "esc_Xist_F2_R2": "esc_Xist_F2_R2\nF: esc-Xist-0616F2 / R: esc-Xist-0616R2",
    "esc_Kdm5c": "esc_Kdm5c\nF: esc-Kdm5c-0616F / R: esc-Kdm5c-0616R",
    "esc_Kdm6a": "esc_Kdm6a\nF: esc-Kdm6a-0616F / R: esc-Kdm6a-0616R",
    "si_Tsix": "si_Tsix\nF: si-Tsix-0616F / R: si-Tsix-0616R",
    "si_RNF12_Rlim": "si_RNF12_Rlim\nF: si-RNF12-0616F / R: si-RNF12-0616R",
    "si_Rbmx": "si_Rbmx\nF: si-Rbmx-0616F / R: si-Rbmx-0616R",
}


def plot_condition(df: pd.DataFrame, condition: str, fontsize: int = 17) -> Path:
    amplicon_order = CONDITION_ORDER[condition]
    plot_df = df[df["sample"] == condition].set_index("amplicon_id").reindex(amplicon_order).reset_index()
    plot_df["Label"] = plot_df["amplicon_id"].map(PRIMER_LABELS)
    plot_df["sort_order"] = plot_df["amplicon_id"].apply(lambda x: amplicon_order.index(x))
    plot_df = plot_df.sort_values("sort_order", ascending=False)

    fig, ax = plt.subplots(figsize=(10.5, 0.85 * len(plot_df) + 2.2))

    labels = plot_df["Label"].tolist()
    cast_pct = plot_df["cast_percent"].fillna(0).tolist()
    b6_pct = plot_df["b6_percent"].fillna(0).tolist()

    ax.barh(labels, cast_pct, label="Cast Reads", color="#9a3324", alpha=0.8)
    ax.barh(labels, b6_pct, left=cast_pct, label="B6 Reads", color="#00274c", alpha=0.8)

    for i, (_, row) in enumerate(plot_df.iterrows()):
        cast = int(row["cast_count"]) if not pd.isna(row["cast_count"]) else 0
        b6 = int(row["b6_count"]) if not pd.isna(row["b6_count"]) else 0
        aligned = int(row["aligned_reads"]) if not pd.isna(row["aligned_reads"]) else 0
        informative = int(row["informative_observations"]) if not pd.isna(row["informative_observations"]) else 0

        if row["status"] == "OK":
            ax.text(101, i, f"{row.cast_percent:.1f}%", va="center", fontweight="bold", fontsize=fontsize)
            ax.text(1, i, f"{cast:,}", va="center", ha="left", color="white", fontsize=fontsize - 7)
            ax.text(99, i, f"{b6:,}", va="center", ha="right", color="white", fontsize=fontsize - 7)
            if row.cast_percent < 4 or row.b6_percent < 4:
                ax.text(112, i, f"Cast={cast:,}; B6={b6:,}", va="center", fontsize=fontsize - 9)
        else:
            ax.text(101, i, "NA", va="center", fontweight="bold", fontsize=fontsize)
            ax.text(
                112,
                i,
                f"{row.status}; aligned={aligned:,}; informative={informative:,}; Cast={cast:,}; B6={b6:,}",
                va="center",
                fontsize=fontsize - 9,
            )

    ax.text(101, len(plot_df) - 0.5, "Cast%", va="bottom", fontweight="bold", fontsize=fontsize)
    ax.set_title(condition, fontsize=fontsize + 3, fontweight="bold", pad=12)
    ax.set_xlabel("Percent of Reads", fontsize=fontsize)
    ax.tick_params(axis="both", which="major", labelsize=fontsize)

    for spine in ax.spines.values():
        spine.set_linewidth(2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0, 150)

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=False, fontsize=fontsize)
    out = FIG_DIR / f"{condition}_allele_barplot_by_amplicon.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def write_report(df: pd.DataFrame, outputs: list[Path]) -> None:
    lines = [
        "# 20260629 AmpSeq Figure Report",
        "",
        "These figures use the same plotting style as the earlier pipeline notebook: one condition per horizontal Cast/B6 stacked bar plot. The y-axis rows are the amplicon primer pairs used for that sample group. Raw Cast and B6 read counts are printed inside the bars using the smaller count font; the Cast percentage is printed to the right.",
        "",
        "Uninformative gDNA rows are not missing amplicons. They are shown as `NA` on the right side, with aligned read counts and zero informative SNP observations.",
        "",
        "## Figure Files",
        "",
    ]
    for path in outputs:
        lines.append(f"- `{path.name}`")
    lines.extend(["", "## Per-Condition Notes", ""])
    for condition in CONDITION_ORDER:
        lines.append(f"### {condition}")
        sub = df[df["sample"] == condition].set_index("amplicon_id").reindex(CONDITION_ORDER[condition])
        for amp, row in sub.iterrows():
            cast = int(row["cast_count"]) if not pd.isna(row["cast_count"]) else 0
            b6 = int(row["b6_count"]) if not pd.isna(row["b6_count"]) else 0
            aligned = int(row["aligned_reads"]) if not pd.isna(row["aligned_reads"]) else 0
            informative = int(row["informative_observations"]) if not pd.isna(row["informative_observations"]) else 0
            if row["status"] == "OK":
                lines.append(f"- `{amp}`: Cast {row.cast_percent:.1f}% with Cast={cast:,}, B6={b6:,}, informative={informative:,}, aligned={aligned:,}.")
            else:
                lines.append(f"- `{amp}`: `{row.status}`; aligned={aligned:,}, informative={informative:,}, Cast={cast:,}, B6={b6:,}.")
        lines.append("")
    (FIG_DIR / "20260629_ampseq_figure_report.md").write_text("\n".join(lines))


def notebook_source() -> str:
    return Path(__file__).read_text()


def write_notebook() -> None:
    nb = {
        "cells": [
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "# 20260629 AmpSeq Condition Plots\n",
                    "\n",
                    "One condition per plot, matching the earlier pipeline notebook style. Rows are amplicon primer pairs.\n",
                ],
            },
            {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": notebook_source().splitlines(True)},
        ],
        "metadata": {
            "kernelspec": {"display_name": "bioinfo", "language": "python", "name": "bioinfo"},
            "language_info": {"name": "python"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    (FIG_DIR / "plot_20260629_ampseq_figures.ipynb").write_text(json.dumps(nb, indent=2))


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(SUMMARY_CSV)
    outputs = [plot_condition(df, condition) for condition in CONDITION_ORDER]
    write_report(df, outputs)
    write_notebook()


if __name__ == "__main__":
    main()
