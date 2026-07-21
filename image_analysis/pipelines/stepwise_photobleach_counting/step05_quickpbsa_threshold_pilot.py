"""Step 5: pool traces across datasets, sweep quickPBSA thresholds, and select half-step thresholds."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter

from pbsa_shared import atomic_csv, load_config, matplotlib_setup, output_dirs, stage_timer
from quickpbsa_compat import result_summary, run_quickpbsa
from step01_inspect_inputs import load_manifest
from step04_extract_background_corrected_traces import trace_paths


def pooled_pilot_root(cfg: dict[str, Any]) -> Path:
    return Path(cfg["output_root"]) / "05_quickpbsa_threshold_pilot"


def selected_thresholds_path(cfg: dict[str, Any]) -> Path:
    return pooled_pilot_root(cfg) / "selected_quickpbsa_thresholds.csv"


def frame_columns(table: pd.DataFrame) -> list[str]:
    return sorted([column for column in table.columns if str(column).isdigit()], key=lambda value: int(value))


def load_sampled_traces(path: Path, maximum: int) -> pd.DataFrame:
    table = pd.read_csv(path)
    points = table[table.roi_class == "point_punctum"]
    selected = points if not points.empty else table
    source = "point_punctum" if not points.empty else "all_accepted_rois_fallback"
    if len(selected) > maximum:
        selected = selected.iloc[np.linspace(0, len(selected) - 1, maximum).round().astype(int)]
    selected = selected.copy().reset_index(drop=True)
    selected["threshold_estimation_source"] = source
    return selected


def estimate_step_scale(values: np.ndarray) -> tuple[float, float, np.ndarray]:
    if not values.size:
        return float("nan"), float("nan"), np.array([])
    smoothed = median_filter(values.astype(np.float32), size=(1, 3), mode="nearest")
    start = smoothed.shape[1] // 2
    differences = -np.diff(smoothed[:, start:], axis=1).ravel()
    center = float(np.median(differences))
    noise = float(1.4826 * np.median(np.abs(differences - center)))
    jumps = differences[differences > center + 4 * max(noise, 1e-6)]
    if jumps.size:
        cutoff = np.percentile(jumps, 50)
        estimate = float(np.median(jumps[jumps <= cutoff]))
    else:
        estimate = float(max(8 * noise, 1.0))
    return estimate, noise, jumps


def select_half_step_threshold(step_scale: float) -> float:
    if not np.isfinite(step_scale) or step_scale <= 0:
        raise RuntimeError(f"Cannot select quickPBSA threshold from invalid step scale: {step_scale}")
    return float(0.5 * step_scale)


def write_quickpbsa_input(table: pd.DataFrame, path: Path) -> None:
    metadata = [column for column in ["dataset", "filename", "fov", "roi_id", "roi_class"] if column in table]
    atomic_csv(table[metadata + frame_columns(table)], path)


def pilot_profile(profile: str, sample: pd.DataFrame, cfg: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any], Path]:
    root = pooled_pilot_root(cfg) / profile
    root.mkdir(parents=True, exist_ok=True)
    # Profiles may contain 1,000- to 3,000-frame movies. Use their shared full
    # prefix so quickPBSA receives a rectangular, NaN-free pilot matrix.
    frames = [column for column in frame_columns(sample) if sample[column].notna().all()]
    if not frames:
        raise RuntimeError(f"No shared finite frames for pooled acquisition profile {profile}")
    metadata = [column for column in sample.columns if not str(column).isdigit()]
    sample = sample[metadata + frames].copy()
    values = sample[frames].to_numpy(float)
    step_scale, noise, jumps = estimate_step_scale(values)
    selected_threshold = select_half_step_threshold(step_scale)
    multipliers = [float(value) for value in cfg["quickpbsa"]["candidate_threshold_multipliers"]]
    infile = root / f"{profile}__pooled_pilot_integrated_difference.csv"
    write_quickpbsa_input(sample, infile)
    records = []
    for multiplier in multipliers:
        threshold = float(step_scale * multiplier)
        label = f"threshold_{threshold:.4g}".replace(".", "p")
        result = run_quickpbsa(infile, root / label, threshold, cfg["quickpbsa"], maxiter=min(100, int(cfg["quickpbsa"]["maxiter"])))
        summary = result_summary(result)
        accepted = summary[summary.flag == 1] if "flag" in summary else summary.iloc[0:0]
        records.append({
            "acquisition_profile": profile, "threshold_multiplier": multiplier,
            "candidate_threshold": threshold, "estimated_single_step": step_scale,
            "late_difference_noise": noise, "traces": len(summary), "accepted_traces": len(accepted),
            "accepted_fraction": len(accepted) / len(summary) if len(summary) else 0,
            "median_accepted_count": float(accepted.photobleaching_step_count.median()) if len(accepted) else np.nan,
            "selected": bool(np.isclose(multiplier, 0.5)),
            "result_file": str(root / label / f"{infile.stem}_result.csv"),
        })
    sweep = pd.DataFrame.from_records(records)
    atomic_csv(sweep, root / f"{profile}__threshold_sweep.csv")
    source_counts = sample.threshold_estimation_source.value_counts().to_dict()
    explicit = cfg["quickpbsa"].get("thresholds", {}).get(profile)
    selection = {
        "acquisition_profile": profile,
        "estimated_single_step": step_scale,
        "selected_threshold": float(explicit) if explicit is not None else selected_threshold,
        "selection_source": "config_override" if explicit is not None else "automatic_half_step",
        "pilot_traces": len(sample), "datasets": int(sample.dataset.nunique()),
        "point_punctum_traces": int(source_counts.get("point_punctum", 0)),
        "fallback_traces": int(source_counts.get("all_accepted_rois_fallback", 0)),
        "late_difference_noise": noise,
    }

    matplotlib_setup(cfg["output_root"])
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for trace in values[:int(cfg["qc"]["max_trace_panels_per_class"])]:
        axes[0, 0].plot(trace, alpha=.35)
    axes[0, 0].set_title(f"Pooled pilot traces — {profile}")
    axes[0, 0].set_xlabel("Frame"); axes[0, 0].set_ylabel("Integrated difference")
    axes[0, 1].hist(jumps, bins=50, histtype="step")
    axes[0, 1].axvline(step_scale, color="red", label="estimated single step")
    axes[0, 1].legend(); axes[0, 1].set_title("Late downward jumps")
    axes[1, 0].plot(sweep.candidate_threshold, sweep.accepted_fraction, "o-")
    axes[1, 0].axvline(selection["selected_threshold"], color="red", ls="--")
    axes[1, 0].set_xlabel("quickPBSA threshold"); axes[1, 0].set_ylabel("Accepted fraction")
    axes[1, 1].plot(sweep.candidate_threshold, sweep.median_accepted_count, "o-")
    axes[1, 1].axvline(selection["selected_threshold"], color="red", ls="--")
    axes[1, 1].set_xlabel("quickPBSA threshold"); axes[1, 1].set_ylabel("Median accepted count")
    fig.tight_layout()
    plot = root / f"{profile}__quickpbsa_threshold_pilot.png"
    fig.savefig(plot, dpi=200); plt.close(fig)
    return sweep, selection, plot


def collect_pilot_traces(cfg: dict[str, Any], max_files: int | None = None) -> dict[str, list[pd.DataFrame]]:
    grouped: dict[str, list[pd.DataFrame]] = {}
    maximum = max(1, int(cfg["qc"].get("threshold_pilot_traces_per_fov", 24)))
    for dataset in cfg["datasets"]:
        accepted = load_manifest(dataset)
        accepted = accepted[accepted.status == "accepted"]
        if max_files is not None:
            accepted = accepted.iloc[:max_files]
        for row in accepted.to_dict("records"):
            path = trace_paths(dataset["output_dir"], row["fov"])["integrated_difference"]
            if not path.exists():
                raise FileNotFoundError(f"Run trace extraction first: {path}")
            table = load_sampled_traces(path, maximum)
            if table.empty:
                continue
            table.insert(0, "dataset", dataset["name"])
            table.insert(1, "filename", row["filename"])
            table.insert(2, "fov", row["fov"])
            grouped.setdefault(row["acquisition_profile"], []).append(table)
    return grouped


def save_stage_report(sweeps: pd.DataFrame, selections: pd.DataFrame, plots: list[Path], cfg: dict[str, Any]) -> Path:
    matplotlib_setup(cfg["output_root"])
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    report = pooled_pilot_root(cfg) / "quickpbsa_threshold_pilot_QC_report.pdf"
    with PdfPages(report) as pdf:
        fig, ax = plt.subplots(figsize=(11, 8.5)); ax.axis("off")
        text = "Pooled quickPBSA threshold pilot\n\nSelected thresholds:\n" + selections.to_string(index=False)
        text += "\n\nSensitivity sweep:\n" + sweeps.to_string(index=False)
        ax.text(.02, .98, text, va="top", family="monospace", fontsize=6)
        pdf.savefig(fig); plt.close(fig)
        for path in plots:
            fig, ax = plt.subplots(figsize=(11, 8.5)); ax.imshow(plt.imread(path)); ax.axis("off")
            pdf.savefig(fig); plt.close(fig)
    return report


def pooled_threshold_pilot_stage(cfg: dict[str, Any], *, max_files: int | None = None) -> pd.DataFrame:
    root = pooled_pilot_root(cfg); root.mkdir(parents=True, exist_ok=True)
    grouped = collect_pilot_traces(cfg, max_files=max_files)
    sweeps, selections, plots = [], [], []
    with stage_timer(cfg["output_root"], "05_quickpbsa_threshold_pilot", {"profiles": len(grouped)}):
        for profile, tables in sorted(grouped.items()):
            print(f"[pooled quickPBSA threshold pilot] {profile}", flush=True)
            sweep, selection, plot = pilot_profile(profile, pd.concat(tables, ignore_index=True, sort=False), cfg)
            sweeps.append(sweep); selections.append(selection); plots.append(plot)
        if not selections:
            raise RuntimeError("No usable traces were available for quickPBSA threshold selection")
        sweep_table = pd.concat(sweeps, ignore_index=True)
        selection_table = pd.DataFrame.from_records(selections)
        atomic_csv(sweep_table, root / "quickpbsa_threshold_pilot_QC_summary.csv")
        atomic_csv(selection_table, selected_thresholds_path(cfg))
        save_stage_report(sweep_table, selection_table, plots, cfg)
    return selection_table


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True); parser.add_argument("--max-files", type=int)
    args = parser.parse_args(argv)
    print(pooled_threshold_pilot_stage(load_config(args.config), max_files=args.max_files).to_string(index=False))


if __name__ == "__main__":
    main()
