#!/usr/bin/env python3
"""Create a consolidated B6/Cast allele stacked-bar plot."""

import argparse
import csv
import os

import matplotlib.pyplot as plt
import pandas as pd


def sanitize_sample_name(name):
    import re

    name = name.strip()
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    return name.strip("._-")


def sample_order_from_map(sample_map):
    if not sample_map:
        return None
    with open(sample_map, newline="") as handle:
        reader = csv.DictReader(handle)
        if "SampleName" not in (reader.fieldnames or []):
            raise ValueError("sample map must include SampleName")
        order = [sanitize_sample_name(row["SampleName"]) for row in reader if row.get("SampleName", "").strip()]
    return order or None


def load_summary(summary_csv):
    df = pd.read_csv(summary_csv)
    required = {"Sample", "B6", "Cast", "Cast_Ratio"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"summary is missing required columns: {', '.join(sorted(missing))}")
    df = df.copy()
    df["Cast_percent"] = df["Cast_Ratio"] * 100
    return df


def plot_allele_summary(df, output_path, sample_order=None, fontsize=25, show_raw_counts=True):
    if sample_order:
        df = df[df["Sample"].isin(sample_order)].copy()
        df["sort_order"] = df["Sample"].apply(lambda x: sample_order.index(x))
        df = df.sort_values("sort_order", ascending=False)
    else:
        df = df.sort_values("Sample", ascending=False)

    if df.empty:
        raise ValueError("no samples available to plot")

    fig, ax = plt.subplots(figsize=(8, 0.8 * len(df) + 2))
    labels = df["Sample"].tolist()
    cast_pct = df["Cast_percent"].tolist()
    b6_pct = [100 - x for x in cast_pct]

    ax.barh(labels, cast_pct, label="Cast Reads", color="#9a3324", alpha=0.8)
    ax.barh(labels, b6_pct, left=cast_pct, label="B6 Reads", color="#00274c", alpha=0.8)

    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(
            101,
            i,
            f"{row.Cast_percent:.1f}%",
            va="center",
            fontweight="bold",
            fontsize=fontsize,
        )
        if show_raw_counts:
            ax.text(
                1,
                i,
                f"{int(row.Cast)}",
                va="center",
                ha="left",
                color="white",
                fontsize=fontsize - 7,
            )
            ax.text(
                99,
                i,
                f"{int(row.B6)}",
                va="center",
                ha="right",
                color="white",
                fontsize=fontsize - 7,
            )

    ax.text(
        101,
        len(df) - 0.5,
        "Cast%",
        va="bottom",
        fontweight="bold",
        fontsize=fontsize,
    )

    ax.set_xlabel("Percent of Reads", fontsize=fontsize)
    ax.tick_params(axis="both", which="major", labelsize=fontsize)

    for spine in ax.spines.values():
        spine.set_linewidth(2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_xlim(0, 100)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.1),
        ncol=2,
        frameon=False,
        fontsize=fontsize,
    )

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, help="Nanopore results directory, e.g. results/mode_1")
    parser.add_argument("--summary-csv", help="Allele summary CSV. Default: results-dir/quantification/allele_quantification_summary.csv")
    parser.add_argument("--sample-map", default=None, help="Optional sample map whose SampleName order controls row order")
    parser.add_argument("--output", help="Output PNG. Default: results-dir/allele_barplot_consolidated.png")
    parser.add_argument("--fontsize", type=int, default=25)
    parser.add_argument("--hide-raw-counts", action="store_true", help="Hide raw Cast/B6 read counts inside the bars")
    args = parser.parse_args()

    summary_csv = args.summary_csv or os.path.join(args.results_dir, "quantification", "allele_quantification_summary.csv")
    output = args.output or os.path.join(args.results_dir, "allele_barplot_consolidated.png")
    order = sample_order_from_map(args.sample_map) if args.sample_map else None
    df = load_summary(summary_csv)
    plot_allele_summary(
        df,
        output,
        sample_order=order,
        fontsize=args.fontsize,
        show_raw_counts=not args.hide_raw_counts,
    )
    print(f"Saved plot: {output}")


if __name__ == "__main__":
    main()
