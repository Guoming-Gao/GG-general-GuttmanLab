"""Step 7: generate per-dataset outputs and the combined biological comparison report."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from pbsa_shared import atomic_csv, condition_metadata, load_config, matplotlib_setup, output_dirs, stage_timer
from step01_inspect_inputs import load_manifest
from step05_quickpbsa_threshold_pilot import pooled_pilot_root, selected_thresholds_path
from step06_count_photobleaching_steps import count_paths


CONDITION_ORDER = ["SHA noDox", "SHA Dox", "dSPEN FL", "dSPEN dRRM"]
ROI_CLASSES = ["point_punctum", "extended_condensate"]


def combined_report_root(cfg: dict[str, Any]) -> Path:
    return Path(cfg["output_root"]) / "07_summary_reports"


def condition_from_filename(filename: str) -> str:
    return str(condition_metadata(filename)["condition"])


def _true(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return values.astype(str).str.lower().isin({"true", "1", "yes"})


def load_dataset_counts(dataset: dict[str, Any], manifest: pd.DataFrame) -> pd.DataFrame:
    tables = []
    for row in manifest[manifest.status == "accepted"].to_dict("records"):
        path = count_paths(dataset["output_dir"], row["fov"])["canonical"]
        if not path.exists():
            raise FileNotFoundError(f"Run photobleaching step counting first: {path}")
        table = pd.read_csv(path)
        metadata = condition_metadata(row["filename"])
        table["dataset"] = dataset["name"]
        table["condition"] = metadata["condition"]
        table["condition_key"] = metadata["condition_key"]
        table["analysis_set"] = metadata["analysis_set"]
        table["is_primary_comparison"] = metadata["is_primary_comparison"]
        tables.append(table)
    return pd.concat(tables, ignore_index=True, sort=False) if tables else pd.DataFrame()


def summary_table(counts: pd.DataFrame) -> pd.DataFrame:
    records = []
    group_columns = ["analysis_set", "condition", "roi_class"]
    for keys, group in counts.groupby(group_columns, dropna=False):
        accepted = group[_true(group.accepted_count)]
        valid = accepted[accepted.validated_range_status == "within_validated_range"]
        records.append({
            "analysis_set": keys[0], "condition": keys[1], "roi_class": keys[2],
            "total_rois": len(group), "accepted_rois": len(accepted),
            "accepted_fraction": len(accepted) / len(group) if len(group) else 0,
            "within_validated_range": len(valid),
            "median_count": float(valid.photobleaching_step_count.median()) if len(valid) else np.nan,
            "mean_count": float(valid.photobleaching_step_count.mean()) if len(valid) else np.nan,
            "outside_validated_range": int((accepted.validated_range_status == "outside_validated_range").sum()),
        })
    return pd.DataFrame.from_records(records)


def per_fov_summary(counts: pd.DataFrame) -> pd.DataFrame:
    records = []
    for keys, group in counts.groupby(["dataset", "filename", "fov", "analysis_set", "condition", "roi_class"], dropna=False):
        accepted = group[_true(group.accepted_count)]
        records.append({
            "dataset": keys[0], "filename": keys[1], "fov": keys[2], "analysis_set": keys[3],
            "condition": keys[4], "roi_class": keys[5], "total_rois": len(group),
            "accepted_rois": len(accepted),
            "median_count": float(accepted.photobleaching_step_count.median()) if len(accepted) else np.nan,
            "outside_validated_range": int((accepted.validated_range_status == "outside_validated_range").sum()),
        })
    return pd.DataFrame.from_records(records)


def make_condition_histogram(counts: pd.DataFrame, roi_class: str):
    matplotlib_setup(str(counts.attrs.get("output_dir", "/tmp")))
    import matplotlib.pyplot as plt
    selected = counts[
        _true(counts.is_primary_comparison)
        & _true(counts.accepted_count)
        & (counts.roi_class == roi_class)
    ].copy()
    numeric = pd.to_numeric(selected.photobleaching_step_count, errors="coerce")
    selected = selected[numeric.notna()].copy()
    selected["photobleaching_step_count"] = numeric[numeric.notna()]
    maximum = int(max(41, selected.photobleaching_step_count.max() if len(selected) else 41))
    bins = np.arange(-0.5, maximum + 1.5, 1)
    fig, ax = plt.subplots(figsize=(10, 6.5))
    plotted = 0
    for condition in CONDITION_ORDER:
        values = selected.loc[selected.condition == condition, "photobleaching_step_count"].to_numpy(float)
        if not len(values):
            continue
        weights = np.full(len(values), 100.0 / len(values))
        ax.hist(values, bins=bins, weights=weights, histtype="step", linewidth=2,
                label=f"{condition} (n={len(values):,})")
        plotted += 1
    ax.axvline(40.5, color="black", ls="--", lw=1.2, label="Published validated range (≤40)")
    ax.set_xlabel("Photobleaching step count")
    ax.set_ylabel("Accepted ROIs per integer bin (%)")
    ax.set_title("Point puncta" if roi_class == "point_punctum" else "Extended condensates")
    ax.grid(axis="y", alpha=.2)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(.5, .015), ncol=4, frameon=False)
    fig.subplots_adjust(bottom=.22)
    return fig, ax


def save_comparison_figures(counts: pd.DataFrame, cfg: dict[str, Any]) -> dict[str, Path]:
    import matplotlib.pyplot as plt
    root = combined_report_root(cfg); root.mkdir(parents=True, exist_ok=True)
    counts.attrs["output_dir"] = cfg["output_root"]
    paths = {
        "point_punctum": root / "point_puncta_condition_step_count_histogram.png",
        "extended_condensate": root / "extended_condensates_condition_step_count_histogram.png",
    }
    for roi_class, path in paths.items():
        fig, _ = make_condition_histogram(counts, roi_class)
        fig.savefig(path, dpi=300, bbox_inches="tight"); plt.close(fig)
    return paths


def _stage_tables(cfg: dict[str, Any], stage_key: str, filename: str) -> pd.DataFrame:
    tables = []
    for dataset in cfg["datasets"]:
        path = output_dirs(dataset["output_dir"])[stage_key] / filename
        if path.exists():
            table = pd.read_csv(path)
            if "dataset" not in table.columns: table.insert(0, "dataset", dataset["name"])
            tables.append(table)
    return pd.concat(tables, ignore_index=True, sort=False) if tables else pd.DataFrame()


def save_combined_report(counts: pd.DataFrame, summary: pd.DataFrame, fov_summary: pd.DataFrame,
                         figures: dict[str, Path], cfg: dict[str, Any]) -> Path:
    matplotlib_setup(cfg["output_root"])
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    root = combined_report_root(cfg); report = root / "combined_stepwise_photobleach_report.pdf"
    thresholds = pd.read_csv(selected_thresholds_path(cfg))
    inventory = _stage_tables(cfg, "input_inspection", "input_manifest.csv")
    stage_specs = [
        ("drift_correction", "drift_correction_QC_summary.csv", "Drift correction QC",
         ["filename", "status", "frames", "median_residual_shift_px", "p95_residual_shift_px", "invalid_block_estimates", "crop_margin_px"]),
        ("roi_detection", "roi_detection_and_growth_QC_summary.csv", "ROI detection and growth QC",
         ["filename", "seeds", "rejected_seeds", "rois", "point_puncta", "extended_condensates", "growth_limit_warnings"]),
        ("traces", "background_corrected_trace_QC_summary.csv", "Background-corrected trace QC",
         ["filename", "accepted_rois", "rejected_rois", "point_puncta", "extended_condensates", "median_late_baseline", "median_late_slope"]),
        ("step_counts", "photobleaching_step_counting_QC_summary.csv", "quickPBSA counting QC",
         ["filename", "rois", "accepted", "accepted_fraction", "point_puncta_accepted", "extended_condensates_accepted", "outside_validated_range"]),
    ]
    stage_pages = []
    for stage_key, filename, title, columns in stage_specs:
        table = _stage_tables(cfg, stage_key, filename)
        if not table.empty:
            stage_pages.append((title, table[[column for column in columns if column in table]]))
    primary = summary[summary.analysis_set == "primary_standardized"]
    point = primary[primary.roi_class == "point_punctum"].set_index("condition")
    extended = primary[primary.roi_class == "extended_condensate"].set_index("condition")
    key_results = [
        "Key results", "",
        f"SHA noDox point puncta: n={int(point.loc['SHA noDox','accepted_rois'])}, median={point.loc['SHA noDox','median_count']:.1f} steps.",
        f"SHA Dox point puncta: n={int(point.loc['SHA Dox','accepted_rois'])}, median={point.loc['SHA Dox','median_count']:.1f} steps.",
        f"SHA noDox extended condensates: n={int(extended.loc['SHA noDox','accepted_rois'])}, median={extended.loc['SHA noDox','median_count']:.1f} steps.",
        f"SHA Dox extended condensates: n={int(extended.loc['SHA Dox','accepted_rois'])}, median={extended.loc['SHA Dox','median_count']:.1f} steps.",
        "dSPEN FL: no ROIs passed quickPBSA acceptance filters.",
        "dSPEN dRRM: one extended ROI passed (14 steps); no point-punctum comparison is available.",
        "", "Interpretation limits",
        "Counts above 40 are retained but lie outside quickPBSA's published validated range.",
        "The primary comparison excludes five nonstandardized SHA-noDox acquisition-test FOVs.",
        "Histograms pool accepted ROIs and normalize within condition; per-FOV summaries expose replicate structure.",
        "Large grown masks frequently retained late signal and were usually rejected by quickPBSA; extended-condensate results are underpowered.",
    ]
    with PdfPages(report) as pdf:
        pages = [
            ("Combined stepwise photobleaching analysis\n\nPrimary condition summary", summary[summary.analysis_set == "primary_standardized"]),
            ("Acquisition inventory", inventory[[column for column in ["filename", "frames", "acquisition_profile", "status"] if column in inventory]]),
            ("Automatically selected quickPBSA thresholds", thresholds),
            ("Per-FOV summary", fov_summary[[column for column in ["filename", "condition", "roi_class", "total_rois", "accepted_rois", "median_count", "outside_validated_range"] if column in fov_summary]]),
            ("Acquisition-test appendix", summary[summary.analysis_set == "acquisition_test"]),
        ]
        total_pages = len(pages) + len(ROI_CLASSES) + len(stage_pages) + 1
        page_number = 0
        def footer(fig):
            nonlocal page_number
            page_number += 1
            fig.text(.5, .02, f"Page {page_number} of {total_pages}", ha="center", fontsize=8, color="0.4")
        for title, table in pages:
            fig, ax = plt.subplots(figsize=(11, 8.5)); ax.axis("off")
            ax.text(.02, .96, title + "\n\n" + table.to_string(index=False), va="top", family="monospace", fontsize=7)
            footer(fig)
            pdf.savefig(fig); plt.close(fig)
        for roi_class in ROI_CLASSES:
            fig, ax = plt.subplots(figsize=(11, 8.5)); ax.imshow(plt.imread(figures[roi_class])); ax.axis("off")
            footer(fig)
            pdf.savefig(fig); plt.close(fig)
        for title, table in stage_pages:
            fig, ax = plt.subplots(figsize=(11, 8.5)); ax.axis("off")
            ax.text(.02, .96, title + "\n\n" + table.to_string(index=False), va="top", family="monospace", fontsize=6.5)
            footer(fig)
            pdf.savefig(fig); plt.close(fig)
        fig, ax = plt.subplots(figsize=(11, 8.5)); ax.axis("off")
        wrapped_results = "\n".join(textwrap.fill(line, width=95) if line else "" for line in key_results)
        ax.text(.05, .94, wrapped_results, va="top", fontsize=11, linespacing=1.5)
        footer(fig)
        pdf.savefig(fig); plt.close(fig)
    return report


def reporting_stage(cfg: dict[str, Any], dataset: dict[str, Any]) -> dict[str, str]:
    manifest = load_manifest(dataset); root = output_dirs(dataset["output_dir"])["reports"]; root.mkdir(parents=True, exist_ok=True)
    with stage_timer(dataset["output_dir"], "07_summary_reports"):
        counts = load_dataset_counts(dataset, manifest); summary = summary_table(counts)
        atomic_csv(counts, root / "photobleaching_step_counts_all_ROIs.csv")
        atomic_csv(summary, root / "photobleaching_step_count_summary.csv")
    return {"counts": str(root / "photobleaching_step_counts_all_ROIs.csv"), "summary": str(root / "photobleaching_step_count_summary.csv")}


def combined_reporting_stage(cfg: dict[str, Any]) -> dict[str, str]:
    root = combined_report_root(cfg); root.mkdir(parents=True, exist_ok=True)
    with stage_timer(cfg["output_root"], "07_combined_summary_reports"):
        tables = [load_dataset_counts(dataset, load_manifest(dataset)) for dataset in cfg["datasets"]]
        counts = pd.concat(tables, ignore_index=True, sort=False)
        summary = summary_table(counts); fov_summary = per_fov_summary(counts)
        atomic_csv(counts, root / "combined_photobleaching_step_counts_all_ROIs.csv")
        atomic_csv(summary, root / "combined_photobleaching_step_count_summary.csv")
        atomic_csv(fov_summary, root / "combined_photobleaching_step_count_per_FOV.csv")
        figures = save_comparison_figures(counts, cfg)
        report = save_combined_report(counts, summary, fov_summary, figures, cfg)
    return {"counts": str(root / "combined_photobleaching_step_counts_all_ROIs.csv"),
            "summary": str(root / "combined_photobleaching_step_count_summary.csv"),
            "per_fov": str(root / "combined_photobleaching_step_count_per_FOV.csv"),
            "point_figure": str(figures["point_punctum"]), "condensate_figure": str(figures["extended_condensate"]),
            "report": str(report)}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", required=True)
    args = parser.parse_args(argv); cfg = load_config(args.config)
    for dataset in cfg["datasets"]:
        reporting_stage(cfg, dataset)
    for key, value in combined_reporting_stage(cfg).items(): print(f"combined {key}: {value}")


if __name__ == "__main__":
    main()
