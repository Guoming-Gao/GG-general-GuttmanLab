"""Step 7: generate condition-level AIO/saSPT summaries and notebook-style reports."""

from __future__ import annotations

import argparse
import ast
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from spt_shared import atomic_csv, output_dirs, stage_timer
from step01_inspect_inputs import inspect_inputs


METRICS = [
    ("N_steps", "Track length, locations", (5, 50), 30),
    ("displacement_nm", "End-to-end displacement, nm", (0, 2000), 40),
    ("mean_stepsize_nm", "Mean Step Size, nm", (0, 600), 40),
    ("max_d_anytwo_nm", "Maximum Excursion, nm", (0, 2500), 40),
    ("linear_fit_sigma", "Localization Error, nm", (0, 300), 30),
    ("linear_fit_log10D", r"log$_{10}$D$_{app}$, $\mu$m$^2$/s", (-2.5, 1), 25),
    ("alpha", r"$\alpha$ Component", (0, 2), 25),
]


def _style() -> None:
    import seaborn as sns
    sns.set(color_codes=True, style="white")


def _filtered(data: pd.DataFrame, metric: str, cfg: dict[str, Any]) -> pd.DataFrame:
    analysis = cfg["analysis"]; mobile = float(analysis["immobile_stepsize_nm"]); r2 = float(analysis["fit_r2_threshold"])
    if metric == "linear_fit_sigma": return data[(data["mean_stepsize_nm"] > mobile) & (data["linear_fit_R2"] > r2)]
    if metric == "linear_fit_log10D": return data[(data["mean_stepsize_nm"] > mobile) & (data["linear_fit_R2"] > r2) & (data["alpha"] > 0.5)]
    if metric == "alpha": return data[(data["mean_stepsize_nm"] > mobile) & (data["loglog_fit_R2"] > r2) & (data["alpha"] > 0)]
    return data


def _format_axis(ax: Any, xlabel: str) -> None:
    ax.set_xlabel(xlabel, fontsize=15); ax.set_ylabel("Probability", fontsize=15)
    ax.spines[:].set_linewidth(1); ax.tick_params(axis="both", which="major", labelsize=15, direction="in", bottom=True, left=True, length=5, width=1)
    legend = ax.get_legend()
    if legend is not None: legend.set_title(None); legend.set_frame_on(False)


def _save_figure(fig: Any, name: str, report_dir: Path, pdf: Any) -> Path:
    path = report_dir/"plots"/f"{name}.png"; path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight"); pdf.savefig(fig, bbox_inches="tight"); return path


def load_condition_tables(cfg: dict[str, Any]) -> pd.DataFrame:
    paths = sorted((output_dirs(cfg["output_dir"])["analysis"]/"by_condition").glob("SPT_results_AIO_concat-*.csv"))
    if not paths: raise FileNotFoundError("No condition AIO tables; run step06 first")
    tables = []
    for path in paths:
        table = pd.read_csv(path)
        table["plot_label"] = table["condition"].astype(str) + np.where(table["channel"].astype(str) == "spt", "", " · " + table["channel"].astype(str))
        tables.append(table)
    return pd.concat(tables, ignore_index=True, sort=False)


def fraction_tables(data: pd.DataFrame, cfg: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    threshold = float(cfg["analysis"]["immobile_stepsize_nm"]); alpha_threshold = float(cfg["analysis"]["alpha_threshold"])
    rows = []
    for (condition, channel, filename), group in data.groupby(["condition", "channel", "filename"]):
        total = len(group); mobile = group["mean_stepsize_nm"] >= threshold
        constrained = mobile & (group["alpha"] <= alpha_threshold); normal = mobile & (group["alpha"] > alpha_threshold)
        rows.append({
            "condition": condition, "channel": channel, "filename": filename, "N_total": total,
            "N_immobile": int((~mobile).sum()), "N_constrained": int(constrained.sum()), "N_normal": int(normal.sum()),
            "immobile_fraction": float((~mobile).mean()), "constrained_fraction": float(constrained.mean()),
            "normal_fraction": float(normal.mean()),
        })
    replicate = pd.DataFrame.from_records(rows)
    summary = replicate.groupby(["condition", "channel"], as_index=False).agg(
        replicates=("filename", "nunique"), total_trajectories=("N_total", "sum"),
        immobile_fraction_mean=("immobile_fraction", "mean"), immobile_fraction_sem=("immobile_fraction", "sem"),
        constrained_fraction_mean=("constrained_fraction", "mean"), constrained_fraction_sem=("constrained_fraction", "sem"),
        normal_fraction_mean=("normal_fraction", "mean"), normal_fraction_sem=("normal_fraction", "sem"),
    )
    return replicate, summary


def metric_summary(data: pd.DataFrame) -> pd.DataFrame:
    metrics = [item[0] for item in METRICS]
    rows = []
    for (condition, channel), group in data.groupby(["condition", "channel"]):
        for metric in metrics:
            values = pd.to_numeric(group[metric], errors="coerce").dropna()
            rows.append({"condition": condition, "channel": channel, "metric": metric, "N": len(values),
                         "mean": values.mean(), "median": values.median(), "std": values.std(), "sem": values.sem()})
    return pd.DataFrame.from_records(rows)


def cell_count_summary(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (condition, channel, filename), group in data.groupby(["condition", "channel", "filename"]):
        ids = pd.to_numeric(group.get("cell_id", pd.Series(dtype=float)), errors="coerce")
        rows.append({
            "condition": condition, "channel": channel, "filename": filename,
            "cells_with_trajectories": int(ids[ids > 0].nunique()),
            "assigned_trajectories": int((ids > 0).sum()),
            "background_trajectories": int((ids.fillna(0) == 0).sum()),
            "total_trajectories": len(group),
        })
    return pd.DataFrame.from_records(rows)


def generate_report(cfg: dict[str, Any]) -> dict[str, str]:
    cache = Path(cfg["output_dir"])/".cache"/"matplotlib"; cache.mkdir(parents=True, exist_ok=True); os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn as sns
    _style(); report_dir = output_dirs(cfg["output_dir"])["reports"]; report_dir.mkdir(parents=True, exist_ok=True)
    data = load_condition_tables(cfg); replicate, fractions = fraction_tables(data, cfg); metrics = metric_summary(data)
    cells = cell_count_summary(data)
    atomic_csv(replicate, report_dir/"replicate_fractions.csv"); atomic_csv(fractions, report_dir/"condition_summary.csv")
    atomic_csv(metrics, report_dir/"diffusion_metric_summary.csv"); atomic_csv(cells, report_dir/"replicate_cell_counts.csv")
    labels = list(dict.fromkeys(data["plot_label"])); palette = dict(zip(labels, sns.color_palette("colorblind", len(labels))))
    pdf_path = report_dir/"diffusion_comparison_report.pdf"; plot_paths = []
    with PdfPages(pdf_path) as pdf:
        for index, (metric, xlabel, limits, bins) in enumerate(METRICS, 1):
            current = _filtered(data, metric, cfg)
            fig, ax = plt.subplots(figsize=(4, 3))
            sns.histplot(data=current, x=metric, hue="plot_label", hue_order=labels, palette=palette, bins=bins,
                         binrange=limits, stat="probability", common_norm=False, lw=2, element="step", fill=False, ax=ax)
            ax.set_xlim(*limits); _format_axis(ax, xlabel)
            if metric == "mean_stepsize_nm": ax.axvline(float(cfg["analysis"]["immobile_stepsize_nm"]), ls="--", color="gray")
            if metric == "alpha": ax.axvline(float(cfg["analysis"]["alpha_threshold"]), ls="--", color="gray")
            plot_paths.append(str(_save_figure(fig, f"{index:02d}_{metric}_comparison", report_dir, pdf))); plt.close(fig)
        angle_rows = []
        for row in data[data["mean_stepsize_nm"] > float(cfg["analysis"]["immobile_stepsize_nm"])].itertuples():
            try: angles = np.abs(np.asarray(ast.literal_eval(row.list_of_angles), dtype=float))
            except (ValueError, SyntaxError): angles = np.array([])
            angle_rows.extend({"angle": value, "plot_label": row.plot_label} for value in angles)
        angle_data = pd.DataFrame.from_records(angle_rows)
        fig, ax = plt.subplots(figsize=(4, 3))
        if not angle_data.empty:
            sns.histplot(data=angle_data, x="angle", hue="plot_label", hue_order=labels, palette=palette, bins=30,
                         binrange=(0, 180), stat="probability", common_norm=False, lw=2, element="step", fill=False, ax=ax)
        ax.axhline(1/30, color="gray", ls="--", lw=2); ax.set_xlim(0, 180); ax.set_xticks([0, 90, 180]); _format_axis(ax, "Angle, °")
        plot_paths.append(str(_save_figure(fig, "08_turning_angle_comparison", report_dir, pdf))); plt.close(fig)
        saspt_paths = sorted((output_dirs(cfg["output_dir"])["analysis"]/"saspt").glob("saSPT-*.csv")); saspt_rows = []
        for path in saspt_paths:
            table = pd.read_csv(path)
            if table.empty or "diff_coef" not in table: continue
            label = str(table["condition"].iloc[0]) + ("" if str(table["channel"].iloc[0]) == "spt" else f" · {table['channel'].iloc[0]}")
            for diffusion, group in table.groupby("diff_coef"):
                saspt_rows.append({"plot_label": label, "log10D": np.log10(float(diffusion)), "Probability": group["mean_posterior_occupation"].sum()})
        saspt = pd.DataFrame.from_records(saspt_rows); fig, ax = plt.subplots(figsize=(4, 3))
        if not saspt.empty: sns.lineplot(data=saspt, x="log10D", y="Probability", hue="plot_label", palette=palette, lw=3, ax=ax)
        ax.set_xlim(-2, 1); ax.set_xticks([-2, -1, 0, 1]); _format_axis(ax, r"log$_{10}$D$_{app}$, $\mu$m$^2$/s"); ax.set_ylabel("SA Occupation", fontsize=15)
        plot_paths.append(str(_save_figure(fig, "09_saSPT_occupation_comparison", report_dir, pdf))); plt.close(fig)
    return {"pdf": str(pdf_path), "condition_summary": str(report_dir/"condition_summary.csv"), "plots": str(report_dir/"plots")}


def report_stage(cfg: dict[str, Any]) -> dict[str, str]:
    with stage_timer(cfg, "07_report", {}): return generate_report(cfg)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", required=True)
    args = parser.parse_args(argv); cfg, _ = inspect_inputs(args.config)
    for key, value in report_stage(cfg).items(): print(f"{key}: {value}")


if __name__ == "__main__":
    main()
