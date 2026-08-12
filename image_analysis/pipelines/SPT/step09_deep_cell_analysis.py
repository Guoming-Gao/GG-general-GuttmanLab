"""Step 9: selected-cell diffusion analysis for staged 30/100 ms review.

This stage reads an existing ONI SPT result tree.  It never re-runs detection,
tracking, or the original AIO calculation and never modifies their artifacts.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import yaml
from scipy import stats
from scipy.optimize import linear_sum_assignment
from tifffile import TiffFile, imread

from spt_shared import atomic_csv, atomic_json, atomic_text
from step06_calculate_diffusion import run_saspt
from step07_generate_reports import METRICS, _filtered


PHASE_ORDER = ["Before", "FVP 2 h", "Recovery 2 h"]
PHASE_COLORS = {"Before": "#4d4d4d", "FVP 2 h": "#d55e00", "Recovery 2 h": "#0072b2"}
STATE_METRICS = ["immobile_fraction", "constrained_fraction", "normal_fraction"]
CELL_METRICS = [item[0] for item in METRICS] + STATE_METRICS
PAIR_METRICS = ["mean_stepsize_nm", "linear_fit_sigma", "linear_fit_log10D", "alpha"] + STATE_METRICS
DEEP_PLOT_METRICS = [item for item in METRICS if item[0] in {
    "mean_stepsize_nm", "linear_fit_sigma", "linear_fit_log10D", "alpha",
}]
DEFAULT_CONFIG: dict[str, Any] = {
    "schema_version": 3,
    "result_root": "/Volumes/guttman/users/gmgao/Imaging_ProcessedData/SPEN/SPT/20260617_ONI-gmgao-SPEN_SPT-SHA_FVP_bef_2h_rec",
    "trajectory_cutoff": 50,
    "cell_bbox_padding_px": 5,
    "scale_bar_um": 1.0,
    "overlay_dpi": 1200,
    "trajectory_linewidth_pt": 0.1,
    "pairing_minimum_iou": 0.5,
    "bootstrap_iterations": 5000,
    "random_seed": 20260722,
    "output_folders": {
        "30ms": "30ms_comparison",
        "100ms": "100ms_comparison",
        "timing_comparison": "30ms_vs_100ms_comparison",
        "selected_cells": "selected_cells",
    },
    "selected_overlay": {
        "colormap": "plasma",
        "alpha_min": 0.3,
        "alpha_max": 1.0,
        "missing_alpha_color": "#9e9e9e",
        "distinguish_gap_closing": False,
        "legend_frame": False,
    },
    "step_size_distribution": {
        "link_policy": "adjacent_only",
        "minimum_nm": 0.0,
        "maximum_nm": 600.0,
        "histogram_bins": 60,
        "cdf_grid_points": 601,
        "timing_comparison_cohort": "matched_fvp_recovery_independent_before",
    },
    "plots": {
        "raster_format": "png",
        "vector_format": "svg",
        "retain_report_pdf": True,
        "legend_include_cells": True,
        "legend_include_tracks": True,
        "legend_frame": False,
    },
    "phases": {
        "Before": "SHA-before",
        "FVP 2 h": "SHA-FVP2h",
        "Recovery 2 h": "SHA-FVP2h_recover2h",
    },
    "analysis": {
        "immobile_stepsize_nm": 30.0,
        "alpha_threshold": 0.7,
        "fit_r2_threshold": 0.7,
        "focal_depth_um": 0.7,
        "pixel_size_um": 0.117,
    },
}


def _deep_merge(base: dict[str, Any], update: dict[str, Any]) -> dict[str, Any]:
    result = json.loads(json.dumps(base))
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def load_deep_config(path: str | Path) -> dict[str, Any]:
    source = Path(path).expanduser().resolve()
    with source.open() as handle:
        user_config = yaml.safe_load(handle) or {}
    # Schema 2 is the immediately preceding local format.  Promote it before
    # merging so existing overrides survive the new adjacent-step defaults.
    if int(user_config.get("schema_version", 2)) == 2:
        user_config["schema_version"] = 3
    cfg = _deep_merge(DEFAULT_CONFIG, user_config)
    root = Path(cfg["result_root"]).expanduser()
    if not root.is_absolute():
        root = (source.parent / root).resolve()
    cfg["result_root"] = str(root)
    cfg["_config_path"] = str(source)
    if int(cfg["trajectory_cutoff"]) < 1:
        raise ValueError("trajectory_cutoff must be positive")
    if int(cfg["cell_bbox_padding_px"]) < 0:
        raise ValueError("cell_bbox_padding_px cannot be negative")
    if float(cfg["pairing_minimum_iou"]) <= 0 or float(cfg["pairing_minimum_iou"]) > 1:
        raise ValueError("pairing_minimum_iou must be in (0, 1]")
    if int(cfg["schema_version"]) != 3:
        raise ValueError("This implementation requires deep-analysis config schema_version 3")
    alpha_min = float(cfg["selected_overlay"]["alpha_min"])
    alpha_max = float(cfg["selected_overlay"]["alpha_max"])
    if not alpha_min < alpha_max:
        raise ValueError("selected_overlay alpha_min must be below alpha_max")
    step_cfg = cfg["step_size_distribution"]
    if step_cfg["link_policy"] != "adjacent_only":
        raise ValueError("Only adjacent_only step-size distributions are supported")
    if float(step_cfg["minimum_nm"]) >= float(step_cfg["maximum_nm"]):
        raise ValueError("step-size minimum_nm must be below maximum_nm")
    if int(step_cfg["histogram_bins"]) < 1 or int(step_cfg["cdf_grid_points"]) < 2:
        raise ValueError("step-size histogram_bins and cdf_grid_points are invalid")
    for folder in ("00_run_metadata", "02_segmentation", "04_trajectories", "05_diffusion_analysis"):
        if not (root / folder).is_dir():
            raise FileNotFoundError(f"Missing existing pipeline folder: {root / folder}")
    return cfg


def _public_config(cfg: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in cfg.items() if not key.startswith("_")}


def comparison_dir(cfg: dict[str, Any], timing_ms: int) -> Path:
    if timing_ms not in {30, 100}:
        raise ValueError("timing_ms must be 30 or 100")
    return Path(cfg["result_root"]) / cfg["output_folders"][f"{timing_ms}ms"]


def timing_comparison_dir(cfg: dict[str, Any]) -> Path:
    return Path(cfg["result_root"]) / cfg["output_folders"]["timing_comparison"]


def selected_cells_dir(cfg: dict[str, Any], timing_ms: int) -> Path:
    return Path(cfg["result_root"]) / cfg["output_folders"]["selected_cells"] / f"{timing_ms}ms"


def _safe(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_") or "value"


def _phase(condition: str, phases: dict[str, str]) -> str:
    # Check Recovery before the broader FVP substring.
    ordered = sorted(phases.items(), key=lambda item: len(item[1]), reverse=True)
    hits = [(label, token) for label, token in ordered if token in condition]
    if not hits:
        raise ValueError(f"Condition {condition!r} does not map to a configured phase")
    longest = len(hits[0][1])
    winners = [label for label, token in hits if len(token) == longest]
    if len(winners) != 1:
        raise ValueError(f"Condition {condition!r} maps ambiguously to {winners}")
    return winners[0]


def _timing_ms(frame_interval_s: float) -> int:
    value = int(round(float(frame_interval_s) * 1000))
    if value not in {30, 100}:
        raise ValueError(f"Unsupported frame interval {frame_interval_s} s")
    return value


def _manifest(root: Path, cfg: dict[str, Any]) -> pd.DataFrame:
    table = pd.read_csv(root / "00_run_metadata" / "input_manifest.csv")
    table = table[table.status == "accepted"].copy()
    table["phase"] = table.condition.map(lambda value: _phase(str(value), cfg["phases"]))
    table["timing_ms"] = table.frame_interval_s.map(_timing_ms)
    return table


def _trajectory_path(root: Path, fov: str) -> Path:
    return root / "04_trajectories" / f"{fov}__spt__trajectories.csv"


def _edge_path(root: Path, fov: str) -> Path:
    return root / "04_trajectories" / f"{fov}__spt__trajectory_edges.csv"


def _aio_path(root: Path, fov: str) -> Path:
    return root / "05_diffusion_analysis" / "per_fov" / f"SPT_results_AIO-{fov}__spt.csv"


def _mask_path(root: Path, fov: str) -> Path:
    return root / "02_segmentation" / f"{fov}__cellposeSAM_masks.tif"


def _dog_path(root: Path, fov: str) -> Path:
    paths = list((root / "01_preprocessed" / "bandpass_float32").glob(f"{fov}__spt__DoG_*_float32.tif"))
    if len(paths) != 1:
        raise FileNotFoundError(f"Expected one SPT DoG stack for {fov}, found {len(paths)}")
    return paths[0]


def build_cell_inventory(cfg: dict[str, Any]) -> pd.DataFrame:
    root = Path(cfg["result_root"])
    manifest = _manifest(root, cfg)
    cells = pd.concat(
        [pd.read_csv(root / "02_segmentation" / f"{row.fov}__cells.csv") for row in manifest.itertuples()],
        ignore_index=True,
    )
    cells = cells.merge(
        manifest[["fov", "condition", "phase", "timing_ms", "frame_interval_s", "pixel_size_um"]],
        on="fov", how="inner", validate="many_to_one",
    )
    retained_parts, aio_parts = [], []
    for row in manifest.itertuples():
        retained = pd.read_csv(_trajectory_path(root, row.fov), usecols=["fov", "cell_id", "trackID"])
        retained_parts.append(retained.drop_duplicates(["fov", "cell_id", "trackID"]))
        aio_parts.append(pd.read_csv(_aio_path(root, row.fov), usecols=["fov", "cell_id", "trackID"]))
    retained = pd.concat(retained_parts, ignore_index=True)
    aio = pd.concat(aio_parts, ignore_index=True)
    retained_counts = (
        retained[retained.cell_id > 0].groupby(["fov", "cell_id"]).trackID.nunique()
        .rename("retained_trajectories")
    )
    aio_counts = (
        aio[aio.cell_id > 0].groupby(["fov", "cell_id"]).trackID.nunique()
        .rename("aio_valid_trajectories")
    )
    cells = cells.join(retained_counts, on=["fov", "cell_id"]).join(aio_counts, on=["fov", "cell_id"])
    for column in ("retained_trajectories", "aio_valid_trajectories"):
        cells[column] = cells[column].fillna(0).astype(int)
    cells["aio_excluded_trajectories"] = cells.retained_trajectories - cells.aio_valid_trajectories
    cells["selected"] = cells.retained_trajectories >= int(cfg["trajectory_cutoff"])
    columns = [
        "phase", "timing_ms", "condition", "fov", "cell_id", "selected",
        "retained_trajectories", "aio_valid_trajectories", "aio_excluded_trajectories",
        "segmentation_source", "centroid_y", "centroid_x", "bbox_min_y", "bbox_min_x",
        "bbox_max_y", "bbox_max_x", "area_px", "frame_interval_s", "pixel_size_um",
    ]
    return cells[columns].sort_values(["timing_ms", "phase", "fov", "cell_id"]).reset_index(drop=True)


def selection_counts(inventory: pd.DataFrame) -> pd.DataFrame:
    return inventory.groupby(["timing_ms", "phase"], sort=False).agg(
        segmented_cells=("cell_id", "size"),
        cells_with_trajectories=("retained_trajectories", lambda values: int((values > 0).sum())),
        selected_cells=("selected", "sum"),
        total_retained_trajectories=("retained_trajectories", "sum"),
        median_retained_trajectories=("retained_trajectories", "median"),
    ).reset_index()


def _config_fingerprint(cfg: dict[str, Any]) -> str:
    payload = yaml.safe_dump(_public_config(cfg), sort_keys=True).encode()
    return hashlib.sha256(payload).hexdigest()


def write_analysis_metadata(cfg: dict[str, Any]) -> None:
    root = Path(cfg["result_root"])
    config_path = root / "deep_analysis_config.yaml"
    resolved = yaml.safe_dump(_public_config(cfg), sort_keys=False)
    # The loaded configuration already contains every user override.  Refreshing
    # the resolved local copy atomically adds new schema defaults without losing
    # those overrides and makes the exact run configuration self-contained.
    atomic_text(resolved, config_path)
    provenance = {
        "created_unix": time.time(),
        "created_local": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        "config_source": cfg["_config_path"],
        "config_fingerprint": _config_fingerprint(cfg),
        "python": sys.version,
        "platform": platform.platform(),
        "command": " ".join(sys.argv),
        "historical_output_dir": (yaml.safe_load((root / "00_run_metadata" / "resolved_config.yaml").read_text()) or {}).get("output_dir"),
        "active_result_root": str(root),
    }
    atomic_json(provenance, root / "deep_analysis_provenance.json")


def inventory_stage(cfg: dict[str, Any]) -> pd.DataFrame:
    root = Path(cfg["result_root"])
    write_analysis_metadata(cfg)
    inventory = build_cell_inventory(cfg)
    atomic_csv(inventory, root / "cell_selection_summary.csv")
    atomic_csv(selection_counts(inventory), root / "cell_selection_counts.csv")
    atomic_json({"inventory": {"status": "complete", "finished_unix": time.time()}}, root / "deep_analysis_stage_status.json")
    return inventory


def _dog_mip(path: Path) -> np.ndarray:
    result = None
    with TiffFile(path) as tif:
        for page in tif.pages:
            frame = np.maximum(page.asarray().astype(np.float32), 0)
            result = frame.copy() if result is None else np.maximum(result, frame)
    if result is None:
        raise ValueError(f"Empty TIFF: {path}")
    return result


def _expanded_bbox(row: Any, shape: tuple[int, int], padding: int) -> tuple[int, int, int, int]:
    height, width = shape
    return (
        max(0, int(row.bbox_min_x) - padding), max(0, int(row.bbox_min_y) - padding),
        min(width, int(row.bbox_max_x) + padding), min(height, int(row.bbox_max_y) + padding),
    )


def selected_overlay_name(row: Any) -> str:
    return (
        f"{int(row.retained_trajectories):05d}_trajectories__{int(row.timing_ms)}ms__"
        f"{_safe(row.phase)}__{_safe(row.fov)}__cell-{int(row.cell_id)}.png"
    )


def trajectory_segment_style(alpha: float, frame_difference: int, cfg: dict[str, Any]) -> dict[str, Any]:
    """Return one consistent segment style; frame gaps are intentionally not distinguished."""
    import matplotlib
    from matplotlib import colors

    overlay = cfg["selected_overlay"]
    if np.isfinite(alpha):
        norm = colors.Normalize(
            vmin=float(overlay["alpha_min"]), vmax=float(overlay["alpha_max"]), clip=True,
        )
        color = matplotlib.colormaps[overlay["colormap"]](norm(float(alpha)))
    else:
        color = overlay["missing_alpha_color"]
    return {
        "color": color,
        "linewidth": float(cfg["trajectory_linewidth_pt"]),
        "linestyle": "-",
        "solid_capstyle": "round",
    }


def render_selected_cell_overlay(
    cfg: dict[str, Any], row: Any, mip: np.ndarray, mask: np.ndarray,
    tracks: pd.DataFrame, edges: pd.DataFrame, track_alpha: dict[int, float], destination: Path,
) -> dict[str, Any]:
    cache = Path(cfg["result_root"]) / ".cache" / "matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from matplotlib.lines import Line2D

    padding = int(cfg["cell_bbox_padding_px"])
    x0, y0, x1, y1 = _expanded_bbox(row, mask.shape, padding)
    crop = mip[y0:y1, x0:x1]
    cell_mask = mask[y0:y1, x0:x1] == int(row.cell_id)
    selected = tracks[tracks.cell_id == int(row.cell_id)].copy()
    selected_ids = set(selected.trackID.unique())
    selected_edges = edges[edges.trackID.isin(selected_ids)].copy()
    finite = crop[np.isfinite(crop)]
    positive = finite[finite > 0]
    upper = float(np.percentile(positive, 99.8)) if positive.size else 1.0
    height, width = crop.shape
    fig, ax = plt.subplots(figsize=(max(2.2, width / 80), max(2.2, height / 80)))
    ax.imshow(crop, cmap="gray", vmin=0, vmax=upper, interpolation="nearest")
    # A direct binary contour is one geometric line, unlike contouring a thick
    # precomputed boundary band, which creates an inner and outer cyan line.
    if np.any(cell_mask):
        ax.contour(cell_mask.astype(float), levels=[0.5], colors="cyan", linewidths=0.35)
    missing_alpha_tracks = 0
    for track_id, group in selected.groupby("trackID", sort=True):
        alpha = float(track_alpha.get(int(track_id), np.nan))
        if not np.isfinite(alpha):
            missing_alpha_tracks += 1
        ordered = group.sort_values("frame")
        values = list(ordered.itertuples())
        for source, target in zip(values[:-1], values[1:]):
            ax.plot(
                [source.x - x0, target.x - x0], [source.y - y0, target.y - y0],
                **trajectory_segment_style(alpha, int(target.frame - source.frame), cfg),
            )
    bar_um = float(cfg["scale_bar_um"])
    bar_px = bar_um / float(row.pixel_size_um)
    bar_y = max(2, height - 7)
    bar_x1 = max(3, width - 5)
    bar_x0 = max(2, bar_x1 - bar_px)
    ax.plot([bar_x0, bar_x1], [bar_y, bar_y], color="white", lw=1.6, solid_capstyle="butt")
    ax.text((bar_x0 + bar_x1) / 2, bar_y - 2.5, f"{bar_um:g} µm", color="white", fontsize=4,
            ha="center", va="bottom", path_effects=[])
    missing_color = cfg["selected_overlay"]["missing_alpha_color"]
    handles = [
        Line2D([0], [0], color="cyan", lw=0.8, label="Cell boundary"),
        Line2D([0], [0], color=missing_color, lw=0.8, label="α unavailable"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=3.5, frameon=False,
              labelcolor="white", handlelength=1.4)
    overlay = cfg["selected_overlay"]
    norm = colors.Normalize(vmin=float(overlay["alpha_min"]), vmax=float(overlay["alpha_max"]), clip=True)
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=overlay["colormap"])
    colorbar = fig.colorbar(scalar, ax=ax, fraction=0.035, pad=0.015)
    colorbar.set_label("α", fontsize=4, rotation=0, labelpad=3)
    colorbar.set_ticks([float(overlay["alpha_min"]), float(overlay["alpha_max"])])
    colorbar.set_ticklabels([f"≤{float(overlay['alpha_min']):g}", f"≥{float(overlay['alpha_max']):g}"])
    colorbar.ax.tick_params(labelsize=3.5, direction="in", length=1.5)
    colorbar.outline.set_visible(False)
    ax.set_title(
        f"{row.phase} · {int(row.timing_ms)} ms · {row.fov} · cell {int(row.cell_id)}\n"
        f"N retained={int(row.retained_trajectories):,}; n={int(row.aio_valid_trajectories):,} trajectories",
        fontsize=5,
    )
    ax.set_xlim(-0.5, width - 0.5)
    ax.set_ylim(height - 0.5, -0.5)
    ax.axis("off")
    fig.tight_layout(pad=0.08)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=int(cfg["overlay_dpi"]), bbox_inches="tight", pad_inches=0.01)
    plt.close(fig)
    return {
        "output": str(destination), "fov": row.fov, "cell_id": int(row.cell_id),
        "timing_ms": int(row.timing_ms), "phase": row.phase,
        "retained_trajectories": int(row.retained_trajectories),
        "aio_valid_trajectories": int(row.aio_valid_trajectories),
        "crop_x0": x0, "crop_y0": y0, "crop_x1": x1, "crop_y1": y1,
        "padding_px": padding, "scale_bar_um": bar_um,
        "plotted_trajectories": int(selected.trackID.nunique()),
        "alpha_available_trajectories": int(selected.trackID.nunique() - missing_alpha_tracks),
        "alpha_unavailable_trajectories": int(missing_alpha_tracks),
        "gap_closing_edges": int((selected_edges.edge_type == "gap_closing").sum()) if not selected_edges.empty else 0,
        "single_boundary_contour": True,
        "gap_edges_styled_identically": True,
    }


def render_selected_cells(
    cfg: dict[str, Any], inventory: pd.DataFrame, timing_ms: int,
    *, resume: bool, force: bool, progress: bool = True,
) -> pd.DataFrame:
    root = Path(cfg["result_root"])
    selected = inventory[(inventory.selected) & (inventory.timing_ms == timing_ms)]
    records: list[dict[str, Any]] = []
    for fov_index, (fov, group) in enumerate(selected.groupby("fov", sort=True), 1):
        mip = _dog_mip(_dog_path(root, fov))
        mask = np.asarray(imread(_mask_path(root, fov)))
        tracks = pd.read_csv(_trajectory_path(root, fov), low_memory=False)
        edges = pd.read_csv(_edge_path(root, fov), low_memory=False)
        aio = pd.read_csv(_aio_path(root, fov), usecols=["trackID", "alpha"], low_memory=False)
        aio["alpha"] = pd.to_numeric(aio.alpha, errors="coerce")
        track_alpha = aio.drop_duplicates("trackID").set_index("trackID").alpha.to_dict()
        if progress:
            print(f"[selected cells] {timing_ms} ms FOV {fov_index}/{selected.fov.nunique()}: {fov}", flush=True)
        for row in group.itertuples(index=False):
            destination = selected_cells_dir(cfg, timing_ms) / selected_overlay_name(row)
            if destination.exists() and resume and not force:
                record = {
                    "output": str(destination), "fov": row.fov, "cell_id": int(row.cell_id),
                    "timing_ms": timing_ms, "phase": row.phase,
                    "retained_trajectories": int(row.retained_trajectories),
                    "aio_valid_trajectories": int(row.aio_valid_trajectories),
                    "plotted_trajectories": int(row.retained_trajectories), "status": "existing",
                }
            elif destination.exists() and not force:
                raise FileExistsError(f"{destination} exists; use --resume or --force")
            else:
                record = render_selected_cell_overlay(cfg, row, mip, mask, tracks, edges, track_alpha, destination)
                record["status"] = "rendered"
            records.append(record)
    result = pd.DataFrame.from_records(records)
    atomic_csv(result, comparison_dir(cfg, timing_ms) / f"{timing_ms}ms_selected_cell_overlays.csv")
    return result


def load_selected_aio(cfg: dict[str, Any], inventory: pd.DataFrame, timing_ms: int) -> pd.DataFrame:
    root = Path(cfg["result_root"])
    selected = inventory[(inventory.selected) & (inventory.timing_ms == timing_ms)]
    tables = []
    for fov, cells in selected.groupby("fov", sort=True):
        table = pd.read_csv(_aio_path(root, fov), low_memory=False)
        table = table[table.cell_id.isin(cells.cell_id)].copy()
        metadata = cells.set_index("cell_id")
        table["phase"] = table.cell_id.map(metadata.phase)
        table["timing_ms"] = timing_ms
        table["pixel_size_um"] = float(cells.pixel_size_um.iloc[0])
        table["frame_interval_s"] = float(cells.frame_interval_s.iloc[0])
        table["retained_trajectories_in_cell"] = table.cell_id.map(metadata.retained_trajectories)
        table["aio_valid_trajectories_in_cell"] = table.cell_id.map(metadata.aio_valid_trajectories)
        tables.append(table)
    if not tables:
        raise ValueError(f"No selected {timing_ms} ms AIO trajectories")
    return pd.concat(tables, ignore_index=True, sort=False)


def _analysis_cfg(cfg: dict[str, Any]) -> dict[str, Any]:
    return {"analysis": cfg["analysis"]}


def cell_metrics(aio: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    rows = []
    step_threshold = float(cfg["analysis"]["immobile_stepsize_nm"])
    alpha_threshold = float(cfg["analysis"]["alpha_threshold"])
    for (timing, phase, fov, cell_id), group in aio.groupby(["timing_ms", "phase", "fov", "cell_id"], sort=True):
        mobile = group.mean_stepsize_nm >= step_threshold
        row: dict[str, Any] = {
            "timing_ms": int(timing), "phase": phase, "fov": fov, "cell_id": int(cell_id),
            "aio_valid_trajectories": len(group),
            "retained_trajectories": int(group.retained_trajectories_in_cell.iloc[0]),
            "immobile_fraction": float((~mobile).mean()),
            "constrained_fraction": float((mobile & (group.alpha <= alpha_threshold)).mean()),
            "normal_fraction": float((mobile & (group.alpha > alpha_threshold)).mean()),
        }
        for metric, *_ in METRICS:
            current = _filtered(group, metric, _analysis_cfg(cfg))
            values = pd.to_numeric(current[metric], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
            row[metric] = float(values.median()) if not values.empty else np.nan
            row[f"{metric}_n"] = len(values)
        rows.append(row)
    return pd.DataFrame.from_records(rows)


def holm_adjust(pvalues: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(pvalues), dtype=float)
    result = np.full(values.shape, np.nan)
    valid = np.flatnonzero(np.isfinite(values))
    if not len(valid):
        return result
    order = valid[np.argsort(values[valid])]
    adjusted_sorted = np.maximum.accumulate(
        np.asarray([(len(order) - rank) * values[index] for rank, index in enumerate(order)])
    )
    adjusted_sorted = np.minimum(adjusted_sorted, 1.0)
    result[order] = adjusted_sorted
    return result


def cliffs_delta(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if not len(left) or not len(right):
        return np.nan
    total = 0
    for value in left:
        total += int(np.sum(value > right)) - int(np.sum(value < right))
    return float(total / (len(left) * len(right)))


def _cluster_arrays(group: pd.DataFrame, metric: str) -> list[np.ndarray]:
    return [
        current[metric].dropna().to_numpy(float)
        for _, current in group.groupby("fov", sort=False)
        if current[metric].notna().any()
    ]


def _cluster_bootstrap_values(clusters: list[np.ndarray], rng: np.random.Generator) -> np.ndarray:
    if not clusters:
        return np.asarray([], dtype=float)
    sampled_clusters = rng.integers(0, len(clusters), size=len(clusters))
    values = []
    for index in sampled_clusters:
        current = clusters[index]
        values.extend(rng.choice(current, current.size, replace=True))
    return np.asarray(values, dtype=float)


def grouped_metric_summary(metrics: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    rng = np.random.default_rng(int(cfg["random_seed"]))
    iterations = int(cfg["bootstrap_iterations"])
    rows = []
    for (timing, phase), group in metrics.groupby(["timing_ms", "phase"], sort=False):
        for metric in CELL_METRICS:
            values = group[metric].dropna().to_numpy(float)
            clusters = _cluster_arrays(group, metric)
            estimates = []
            for _ in range(iterations):
                sampled = _cluster_bootstrap_values(clusters, rng)
                if sampled.size:
                    estimates.append(float(np.median(sampled)))
            low, high = np.quantile(estimates, [0.025, 0.975]) if estimates else (np.nan, np.nan)
            rows.append({
                "timing_ms": int(timing), "phase": phase, "metric": metric,
                "cells": len(values), "fovs": group.fov.nunique(),
                "median": float(np.median(values)) if len(values) else np.nan,
                "q25": float(np.quantile(values, 0.25)) if len(values) else np.nan,
                "q75": float(np.quantile(values, 0.75)) if len(values) else np.nan,
                "cluster_bootstrap_ci_low": float(low), "cluster_bootstrap_ci_high": float(high),
            })
    return pd.DataFrame.from_records(rows)


def phase_statistics(metrics: pd.DataFrame, cfg: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(int(cfg["random_seed"]) + 1)
    iterations = int(cfg["bootstrap_iterations"])
    omnibus_rows, pair_rows = [], []
    for timing, timing_data in metrics.groupby("timing_ms"):
        for metric in CELL_METRICS:
            # Use FOV medians as the inferential unit; cell-level values remain
            # visible in plots and drive the clustered bootstrap intervals.
            fov_values = (
                timing_data.groupby(["phase", "fov"])[metric].median().dropna().rename("value").reset_index()
            )
            groups = [fov_values[fov_values.phase == phase].value.to_numpy(float) for phase in PHASE_ORDER]
            if all(len(values) for values in groups):
                statistic, pvalue = stats.kruskal(*groups)
            else:
                statistic = pvalue = np.nan
            omnibus_rows.append({
                "timing_ms": int(timing), "metric": metric, "test": "Kruskal-Wallis on FOV medians",
                "statistic": statistic, "p_raw": pvalue,
                **{f"n_fov_{_safe(phase)}": len(values) for phase, values in zip(PHASE_ORDER, groups)},
            })
            metric_pairs = []
            for left_phase, right_phase in combinations(PHASE_ORDER, 2):
                left_fov = fov_values[fov_values.phase == left_phase].value.to_numpy(float)
                right_fov = fov_values[fov_values.phase == right_phase].value.to_numpy(float)
                if len(left_fov) and len(right_fov):
                    statistic, pvalue = stats.mannwhitneyu(left_fov, right_fov, alternative="two-sided")
                else:
                    statistic = pvalue = np.nan
                left_cells = timing_data[timing_data.phase == left_phase]
                right_cells = timing_data[timing_data.phase == right_phase]
                left_clusters = _cluster_arrays(left_cells, metric)
                right_clusters = _cluster_arrays(right_cells, metric)
                differences = []
                for _ in range(iterations):
                    left = _cluster_bootstrap_values(left_clusters, rng)
                    right = _cluster_bootstrap_values(right_clusters, rng)
                    if left.size and right.size:
                        differences.append(float(np.median(right) - np.median(left)))
                low, high = np.quantile(differences, [0.025, 0.975]) if differences else (np.nan, np.nan)
                metric_pairs.append({
                    "timing_ms": int(timing), "metric": metric, "left_phase": left_phase,
                    "right_phase": right_phase, "test": "Mann-Whitney U on FOV medians",
                    "statistic": statistic, "p_raw": pvalue,
                    "cliffs_delta_left_vs_right": cliffs_delta(left_fov, right_fov),
                    "median_difference_right_minus_left": (
                        float(np.nanmedian(right_cells[metric]) - np.nanmedian(left_cells[metric]))
                    ),
                    "cluster_bootstrap_difference_ci_low": float(low),
                    "cluster_bootstrap_difference_ci_high": float(high),
                    "n_fov_left": len(left_fov), "n_fov_right": len(right_fov),
                })
            adjusted = holm_adjust([row["p_raw"] for row in metric_pairs])
            for row, value in zip(metric_pairs, adjusted):
                row["p_holm_within_metric"] = value
            pair_rows.extend(metric_pairs)
    return pd.DataFrame.from_records(omnibus_rows), pd.DataFrame.from_records(pair_rows)


def _style() -> None:
    import seaborn as sns
    sns.set_theme(style="white", context="notebook")


def _format_axis(ax: Any, ylabel: str | None = None) -> None:
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=11)
    ax.tick_params(direction="in", bottom=True, left=True, labelsize=9)
    for spine in ax.spines.values():
        spine.set_linewidth(1)
    legend = ax.get_legend()
    if legend is not None:
        legend.set_title(None)
        legend.set_frame_on(False)


def _save(fig: Any, output_dir: Path, stem: str, pdf: Any | None = None) -> list[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    png = output_dir / f"{stem}.png"
    vector = output_dir / f"{stem}.svg"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(vector, bbox_inches="tight")
    if pdf is not None:
        pdf.savefig(fig, bbox_inches="tight")
    return [str(png), str(vector)]


def phase_sample_sizes(metrics: pd.DataFrame, aio: pd.DataFrame) -> dict[str, dict[str, int]]:
    sizes: dict[str, dict[str, int]] = {}
    for phase in PHASE_ORDER:
        phase_metrics = metrics[metrics.phase == phase]
        phase_aio = aio[aio.phase == phase]
        sizes[phase] = {
            "cells": int(phase_metrics[["fov", "cell_id"]].drop_duplicates().shape[0]),
            "trajectories": int(phase_aio[["fov", "trackID"]].drop_duplicates().shape[0]),
        }
    return sizes


def phase_sample_label(phase: str, sizes: dict[str, dict[str, int]]) -> str:
    current = sizes[phase]
    return f"{phase}\nN={current['cells']:,} cells\nn={current['trajectories']:,} trajectories"


def _trajectory_histogram(aio: pd.DataFrame, metric: str, cfg: dict[str, Any], bins: int,
                          limits: tuple[float, float]) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    result = {}
    for phase in PHASE_ORDER:
        group = _filtered(aio[aio.phase == phase], metric, _analysis_cfg(cfg)).copy()
        group[metric] = pd.to_numeric(group[metric], errors="coerce")
        group = group[np.isfinite(group[metric]) & group[metric].between(*limits)]
        if group.empty:
            result[phase] = (np.linspace(*limits, bins + 1), np.zeros(bins))
            continue
        hist, edges = np.histogram(group[metric], bins=bins, range=limits)
        hist = hist / hist.sum() if hist.sum() else hist
        result[phase] = (edges, hist)
    return result


def trajectory_metric_cdf(
    aio: pd.DataFrame, metric: str, cfg: dict[str, Any], limits: tuple[float, float],
    *, grid_points: int = 601,
) -> pd.DataFrame:
    """Return the pooled empirical CDF, giving every valid trajectory one vote."""
    grid = np.linspace(float(limits[0]), float(limits[1]), int(grid_points))
    rows: list[dict[str, Any]] = []
    for phase in PHASE_ORDER:
        group = _filtered(aio[aio.phase == phase], metric, _analysis_cfg(cfg)).copy()
        group[metric] = pd.to_numeric(group[metric], errors="coerce")
        group = group[np.isfinite(group[metric]) & group[metric].between(*limits)]
        if group.empty:
            cumulative = np.zeros_like(grid)
            cells = trajectories = 0
        else:
            ordered = group.sort_values(metric)
            values = ordered[metric].to_numpy(float)
            cumulative = np.searchsorted(values, grid, side="right") / len(values)
            cells = group[["fov", "cell_id"]].drop_duplicates().shape[0]
            trajectories = len(group)
        rows.extend({
            "phase": phase, "metric": metric, "weighting": "pooled_trajectories",
            "x_nm": float(x_value), "cumulative_probability": float(y_value),
            "cells": int(cells), "trajectories": int(trajectories),
        } for x_value, y_value in zip(grid, cumulative))
    return pd.DataFrame.from_records(rows)


def _numeric_list(value: object, *, dtype: type = float) -> np.ndarray:
    return np.fromstring(str(value).strip()[1:-1], sep=",", dtype=dtype)


def extract_adjacent_steps(aio: pd.DataFrame) -> pd.DataFrame:
    """Expand one-frame displacements from diffusion-valid trajectories."""
    records: list[pd.DataFrame] = []
    required = {"timing_ms", "phase", "fov", "cell_id", "trackID", "list_of_t", "list_of_x", "list_of_y", "pixel_size_um"}
    missing = required - set(aio.columns)
    if missing:
        raise ValueError(f"AIO table is missing adjacent-step columns: {sorted(missing)}")
    for row in aio.itertuples(index=False):
        frame = _numeric_list(row.list_of_t, dtype=int)
        x = _numeric_list(row.list_of_x)
        y = _numeric_list(row.list_of_y)
        if not (len(frame) == len(x) == len(y)):
            raise ValueError(f"Coordinate-list length mismatch for {row.fov}/track {row.trackID}")
        adjacent = np.diff(frame) == 1
        if not adjacent.any():
            continue
        distances = np.hypot(np.diff(x), np.diff(y)) * float(row.pixel_size_um) * 1000
        source = frame[:-1][adjacent]
        records.append(pd.DataFrame({
            "timing_ms": int(row.timing_ms), "phase": row.phase, "fov": row.fov,
            "cell_id": int(row.cell_id), "trackID": row.trackID,
            "source_frame": source, "destination_frame": source + 1,
            "frame_difference": 1, "step_size_nm": distances[adjacent],
        }))
    columns = ["timing_ms", "phase", "fov", "cell_id", "trackID", "source_frame",
               "destination_frame", "frame_difference", "step_size_nm"]
    return pd.concat(records, ignore_index=True)[columns] if records else pd.DataFrame(columns=columns)


def step_distribution_curves(aio: pd.DataFrame, cfg: dict[str, Any], *, cohort: str) -> pd.DataFrame:
    steps = extract_adjacent_steps(aio)
    settings = cfg["step_size_distribution"]
    minimum = float(settings["minimum_nm"])
    maximum = float(settings["maximum_nm"])
    edges = np.linspace(minimum, maximum, int(settings["histogram_bins"]) + 1)
    grid = np.linspace(minimum, maximum, int(settings["cdf_grid_points"]))
    rows: list[dict[str, Any]] = []
    for (timing_ms, phase), trajectories in aio.groupby(["timing_ms", "phase"], sort=False):
        current = steps[(steps.timing_ms == timing_ms) & (steps.phase == phase)]
        values = current.step_size_nm.to_numpy(float)
        values = values[np.isfinite(values)]
        counts, _ = np.histogram(values, bins=edges)
        probability = counts / counts.sum() if counts.sum() else counts.astype(float)
        cells = trajectories[["fov", "cell_id"]].drop_duplicates().shape[0]
        n_trajectories = trajectories[["fov", "trackID"]].drop_duplicates().shape[0]
        common = {
            "cohort": cohort, "timing_ms": int(timing_ms), "phase": phase,
            "cells": int(cells), "trajectories": int(n_trajectories), "steps": int(len(values)),
        }
        for left, right, value in zip(edges[:-1], edges[1:], probability):
            rows.append({**common, "curve": "histogram", "x_nm": (left + right) / 2,
                         "y": float(value), "bin_left_nm": left, "bin_right_nm": right})
        ordered = np.sort(values)
        cumulative = np.searchsorted(ordered, grid, side="right") / len(ordered) if len(ordered) else np.zeros_like(grid)
        for x_value, value in zip(grid, cumulative):
            rows.append({**common, "curve": "cdf", "x_nm": x_value, "y": float(value),
                         "bin_left_nm": np.nan, "bin_right_nm": np.nan})
    return pd.DataFrame.from_records(rows)


def _run_selected_saspt(cfg: dict[str, Any], aio: pd.DataFrame, timing_ms: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    output_dir = comparison_dir(cfg, timing_ms)
    rng = np.random.default_rng(int(cfg["random_seed"]) + timing_ms)
    all_rows, balanced_rows = [], []
    for phase in PHASE_ORDER:
        phase_data = aio[aio.phase == phase].copy()
        all_result = run_saspt(
            phase_data, frame_interval_s=timing_ms / 1000,
            pixel_size_um=float(cfg["analysis"]["pixel_size_um"]),
            focal_depth_um=float(cfg["analysis"]["focal_depth_um"]),
        )
        all_result.insert(0, "mode", "all_selected")
        all_result.insert(0, "input_trajectories", len(phase_data))
        all_result.insert(0, "input_cells", phase_data[["fov", "cell_id"]].drop_duplicates().shape[0])
        all_result.insert(0, "phase", phase)
        all_result.insert(0, "timing_ms", timing_ms)
        all_rows.append(all_result)

        balanced_parts = []
        for _, cell in phase_data.groupby(["fov", "cell_id"], sort=True):
            if len(cell) < int(cfg["trajectory_cutoff"]):
                continue
            indices = rng.choice(cell.index.to_numpy(), int(cfg["trajectory_cutoff"]), replace=False)
            balanced_parts.append(cell.loc[indices])
        balanced_input = pd.concat(balanced_parts, ignore_index=True) if balanced_parts else phase_data.iloc[:0]
        balanced_result = run_saspt(
            balanced_input, frame_interval_s=timing_ms / 1000,
            pixel_size_um=float(cfg["analysis"]["pixel_size_um"]),
            focal_depth_um=float(cfg["analysis"]["focal_depth_um"]),
        )
        balanced_result.insert(0, "input_cells", len(balanced_parts))
        balanced_result.insert(1, "input_trajectories", len(balanced_input))
        balanced_result.insert(0, "mode", "balanced_50_aio_per_cell")
        balanced_result.insert(0, "phase", phase)
        balanced_result.insert(0, "timing_ms", timing_ms)
        balanced_rows.append(balanced_result)
    all_table = pd.concat(all_rows, ignore_index=True, sort=False)
    balanced_table = pd.concat(balanced_rows, ignore_index=True, sort=False)
    atomic_csv(all_table, output_dir / f"{timing_ms}ms_selected_cells_saSPT_all.csv")
    atomic_csv(balanced_table, output_dir / f"{timing_ms}ms_selected_cells_saSPT_balanced.csv")
    return all_table, balanced_table


def _plot_timing_report(
    cfg: dict[str, Any], timing_ms: int, inventory: pd.DataFrame, aio: pd.DataFrame,
    metrics: pd.DataFrame, saspt_all: pd.DataFrame, saspt_balanced: pd.DataFrame,
    mean_step_cdf: pd.DataFrame,
) -> list[str]:
    root = Path(cfg["result_root"])
    output_dir = comparison_dir(cfg, timing_ms)
    output_dir.mkdir(parents=True, exist_ok=True)
    cache = root / ".cache" / "matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn as sns
    _style()
    outputs: list[str] = []
    report = output_dir / f"{timing_ms}ms_selected_cell_diffusion_report.pdf"
    current_inventory = inventory[inventory.timing_ms == timing_ms]
    sample_sizes = phase_sample_sizes(metrics, aio)
    sample_labels = [phase_sample_label(phase, sample_sizes) for phase in PHASE_ORDER]
    with PdfPages(report) as pdf:
        counts = selection_counts(current_inventory)
        fig, ax = plt.subplots(figsize=(5.2, 3.6))
        positions = np.arange(len(PHASE_ORDER))
        total = counts.set_index("phase").reindex(PHASE_ORDER).segmented_cells
        selected = counts.set_index("phase").reindex(PHASE_ORDER).selected_cells
        ax.bar(positions, total, color="#d9d9d9", edgecolor="black", label="All segmented")
        ax.bar(positions, selected, color=[PHASE_COLORS[p] for p in PHASE_ORDER], edgecolor="black", label="≥50 tracks")
        for x, value in zip(positions, selected):
            ax.text(x, value + max(total) * 0.02, str(int(value)), ha="center", fontsize=9)
        ax.set_xticks(positions, sample_labels)
        ax.set_title(f"{timing_ms} ms selected-cell yield")
        _format_axis(ax, "Cells")
        ax.legend(frameon=False, loc="upper left")
        outputs += _save(fig, output_dir, f"{timing_ms}ms_01_selected_cell_yield", pdf)
        plt.close(fig)

        display_metrics = [item[0] for item in DEEP_PLOT_METRICS]
        labels = {item[0]: item[1] for item in DEEP_PLOT_METRICS}
        fig, axes = plt.subplots(2, 3, figsize=(10, 6.5))
        for ax, metric in zip(axes.flat, display_metrics):
            sns.boxplot(data=metrics, x="phase", y=metric, hue="phase", order=PHASE_ORDER,
                        palette=PHASE_COLORS, legend=False, width=0.55, showfliers=False, linewidth=1, ax=ax)
            sns.stripplot(data=metrics, x="phase", y=metric, order=PHASE_ORDER, color="black",
                          size=1.8, alpha=0.55, jitter=0.22, ax=ax)
            ax.set_xlabel("")
            ax.set_xticks(np.arange(len(PHASE_ORDER)), PHASE_ORDER, rotation=0, fontsize=7)
            _format_axis(ax, labels[metric])
        from matplotlib.lines import Line2D
        sample_handles = [
            Line2D([0], [0], color=PHASE_COLORS[phase], lw=3,
                   label=phase_sample_label(phase, sample_sizes))
            for phase in PHASE_ORDER
        ]
        axes.flat[-1].axis("off")
        axes.flat[-2].axis("off")
        axes.flat[-2].legend(handles=sample_handles, loc="center", frameon=False, fontsize=9)
        fig.suptitle(f"{timing_ms} ms · cell-level diffusion summaries", fontsize=14)
        fig.tight_layout()
        outputs += _save(fig, output_dir, f"{timing_ms}ms_02_cell_level_diffusion_metrics", pdf)
        plt.close(fig)

        state_long = metrics.melt(
            id_vars=["phase", "fov", "cell_id"], value_vars=STATE_METRICS,
            var_name="state", value_name="fraction",
        )
        state_long["state"] = state_long.state.str.replace("_fraction", "", regex=False).str.title()
        fig, ax = plt.subplots(figsize=(6.8, 4.1))
        sns.boxplot(data=state_long, x="phase", y="fraction", hue="state", order=PHASE_ORDER,
                    showfliers=False, linewidth=1, ax=ax)
        sns.stripplot(data=state_long, x="phase", y="fraction", hue="state", order=PHASE_ORDER,
                      dodge=True, palette={"Immobile": "#333333", "Constrained": "#333333", "Normal": "#333333"},
                      size=1.3, alpha=0.35, ax=ax)
        handles, labels_state = ax.get_legend_handles_labels()
        ax.legend(handles[:3], labels_state[:3], frameon=False, title=None)
        ax.set_ylim(0, 1)
        ax.set_xlabel("")
        ax.set_xticks(np.arange(len(PHASE_ORDER)), sample_labels, fontsize=7)
        ax.set_title(f"{timing_ms} ms · per-cell diffusion-state fractions")
        _format_axis(ax, "Fraction of trajectories")
        outputs += _save(fig, output_dir, f"{timing_ms}ms_03_cell_level_state_fractions", pdf)
        plt.close(fig)

        for index, (metric, xlabel, limits, bins) in enumerate(DEEP_PLOT_METRICS, 4):
            fig, ax = plt.subplots(figsize=(4.5, 3.4))
            histograms = _trajectory_histogram(aio, metric, cfg, bins, limits)
            for phase in PHASE_ORDER:
                edges, hist = histograms[phase]
                centers = (edges[:-1] + edges[1:]) / 2
                ax.step(centers, hist, where="mid", color=PHASE_COLORS[phase], lw=2,
                        label=phase_sample_label(phase, sample_sizes))
            ax.set_xlim(*limits)
            ax.set_xlabel(xlabel, fontsize=11)
            ax.set_title(f"{timing_ms} ms · pooled trajectories")
            _format_axis(ax, "Probability")
            ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))
            fig.subplots_adjust(right=0.67)
            outputs += _save(fig, output_dir, f"{timing_ms}ms_{index:02d}_{metric}_trajectory_histogram", pdf)
            plt.close(fig)

        fig, ax = plt.subplots(figsize=(4.5, 3.4))
        for phase in PHASE_ORDER:
            current = mean_step_cdf[mean_step_cdf.phase == phase]
            ax.step(current.x_nm, current.cumulative_probability, where="post",
                    color=PHASE_COLORS[phase], lw=2, label=phase_sample_label(phase, sample_sizes))
        mean_step_metric = next(item for item in DEEP_PLOT_METRICS if item[0] == "mean_stepsize_nm")
        ax.set_xlim(*mean_step_metric[2])
        ax.set_ylim(0, 1)
        ax.set_xlabel("Trajectory mean step size, nm", fontsize=11)
        ax.set_title(f"{timing_ms} ms · pooled trajectory means")
        _format_axis(ax, "Cumulative probability")
        ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))
        fig.subplots_adjust(right=0.67)
        outputs += _save(fig, output_dir, f"{timing_ms}ms_08_mean_stepsize_nm_trajectory_cdf", pdf)
        plt.close(fig)

        angle_rows = []
        mobile = aio[aio.mean_stepsize_nm > float(cfg["analysis"]["immobile_stepsize_nm"])]
        for (phase, fov, cell_id), cell in mobile.groupby(["phase", "fov", "cell_id"]):
            values = []
            for value in cell.list_of_angles:
                try:
                    values.extend(np.abs(np.asarray(ast.literal_eval(str(value)), dtype=float)))
                except (ValueError, SyntaxError):
                    continue
            if values:
                angle_rows.extend({"phase": phase, "angle": item} for item in values)
        fig, ax = plt.subplots(figsize=(4.5, 3.4))
        angles = pd.DataFrame.from_records(angle_rows)
        for phase in PHASE_ORDER:
            group = angles[angles.phase == phase]
            hist, edges = np.histogram(group.angle, bins=30, range=(0, 180))
            hist = hist / hist.sum() if hist.sum() else hist
            ax.step((edges[:-1] + edges[1:]) / 2, hist, where="mid", color=PHASE_COLORS[phase], lw=2,
                    label=phase_sample_label(phase, sample_sizes))
        ax.set_xlim(0, 180)
        ax.set_xticks([0, 90, 180])
        ax.set_xlabel("Turning angle, °")
        ax.set_title(f"{timing_ms} ms · pooled turning angles")
        _format_axis(ax, "Probability")
        ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))
        fig.subplots_adjust(right=0.67)
        outputs += _save(fig, output_dir, f"{timing_ms}ms_09_turning_angle_pooled_histogram", pdf)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(9, 3.5), sharey=True)
        for ax, table, title in zip(
            axes, [saspt_all, saspt_balanced], ["All selected trajectories", "Balanced: 50 trajectories/cell"],
        ):
            for phase in PHASE_ORDER:
                group = table[table.phase == phase]
                if group.empty:
                    continue
                curve = group.groupby("diff_coef", as_index=False).mean_posterior_occupation.sum()
                input_cells = int(group.input_cells.iloc[0])
                input_tracks = int(group.input_trajectories.iloc[0])
                ax.plot(np.log10(curve.diff_coef), curve.mean_posterior_occupation,
                        color=PHASE_COLORS[phase], lw=2,
                        label=f"{phase}\nN={input_cells:,} cells\nn={input_tracks:,} trajectories")
            ax.set_xlim(-2, 1)
            ax.set_xlabel(r"log$_{10}$D, µm$^2$/s")
            ax.set_title(title)
            _format_axis(ax, "SA occupation" if ax is axes[0] else None)
        for ax in axes:
            ax.legend(frameon=False, fontsize=6.5, loc="upper left", bbox_to_anchor=(1.01, 1))
        fig.suptitle(f"{timing_ms} ms · selected-cell saSPT")
        fig.tight_layout(rect=(0, 0, 0.84, 0.95))
        outputs += _save(fig, output_dir, f"{timing_ms}ms_10_selected_cell_saSPT", pdf)
        plt.close(fig)
    outputs.append(str(report))
    return outputs


def _update_stage_status(root: Path, stage: str, value: dict[str, Any]) -> None:
    path = root / "deep_analysis_stage_status.json"
    state = json.loads(path.read_text()) if path.exists() else {}
    state[stage] = value
    atomic_json(state, path)


def timing_plot_stems(timing_ms: int) -> list[str]:
    stems = [
        f"{timing_ms}ms_01_selected_cell_yield",
        f"{timing_ms}ms_02_cell_level_diffusion_metrics",
        f"{timing_ms}ms_03_cell_level_state_fractions",
    ]
    stems.extend(
        f"{timing_ms}ms_{index:02d}_{metric}_trajectory_histogram"
        for index, (metric, *_rest) in enumerate(DEEP_PLOT_METRICS, 4)
    )
    stems.extend([
        f"{timing_ms}ms_08_mean_stepsize_nm_trajectory_cdf",
        f"{timing_ms}ms_09_turning_angle_pooled_histogram",
        f"{timing_ms}ms_10_selected_cell_saSPT",
    ])
    return stems


def remove_obsolete_timing_plots(cfg: dict[str, Any], timing_ms: int) -> list[str]:
    """Remove superseded plots only after every replacement exists."""
    output_dir = comparison_dir(cfg, timing_ms)
    expected = set(timing_plot_stems(timing_ms))
    for stem in expected:
        for suffix in (".png", ".svg"):
            if not (output_dir / f"{stem}{suffix}").exists():
                raise AssertionError(f"Replacement plot is missing: {stem}{suffix}")
    report = output_dir / f"{timing_ms}ms_selected_cell_diffusion_report.pdf"
    if not report.exists():
        raise AssertionError(f"Replacement report is missing: {report}")
    removed: list[str] = []
    for path in sorted(output_dir.glob(f"{timing_ms}ms_*")):
        if path.suffix in {".png", ".svg"} and path.stem not in expected:
            path.unlink(); removed.append(str(path))
        elif path.suffix == ".pdf" and path != report:
            path.unlink(); removed.append(str(path))
    return removed


def validate_timing_outputs(cfg: dict[str, Any], timing_ms: int, expected_cells: int) -> dict[str, int]:
    output_dir = comparison_dir(cfg, timing_ms)
    overlay_dir = selected_cells_dir(cfg, timing_ms)
    pngs = sorted(output_dir.glob(f"{timing_ms}ms_*.png"))
    svgs = sorted(output_dir.glob(f"{timing_ms}ms_*.svg"))
    pdfs = sorted(output_dir.glob("*.pdf"))
    report = output_dir / f"{timing_ms}ms_selected_cell_diffusion_report.pdf"
    expected_csvs = {
        f"{timing_ms}ms_selected_cell_overlays.csv",
        f"{timing_ms}ms_selected_trajectories_AIO.csv",
        f"{timing_ms}ms_selected_cell_metrics.csv",
        f"{timing_ms}ms_selected_cell_metric_summary.csv",
        f"{timing_ms}ms_phase_omnibus_statistics.csv",
        f"{timing_ms}ms_phase_pairwise_statistics.csv",
        f"{timing_ms}ms_selected_cells_saSPT_all.csv",
        f"{timing_ms}ms_selected_cells_saSPT_balanced.csv",
        f"{timing_ms}ms_adjacent_step_size_curves.csv",
        f"{timing_ms}ms_mean_stepsize_nm_trajectory_cdf.csv",
    }
    observed_csvs = {path.name for path in output_dir.glob("*.csv")}
    if {path.stem for path in pngs} != set(timing_plot_stems(timing_ms)):
        raise AssertionError(f"Unexpected PNG plot set in {output_dir}")
    if {path.stem for path in svgs} != set(timing_plot_stems(timing_ms)):
        raise AssertionError(f"Unexpected SVG plot set in {output_dir}")
    if pdfs != [report]:
        raise AssertionError(f"Only the multipage report PDF may remain in {output_dir}: {pdfs}")
    if observed_csvs != expected_csvs:
        raise AssertionError(f"Unexpected CSV set in {output_dir}: {sorted(observed_csvs ^ expected_csvs)}")
    for svg in svgs:
        if ET.parse(svg).getroot().tag.rsplit("}", 1)[-1] != "svg":
            raise AssertionError(f"Invalid SVG: {svg}")
    pdf_info = subprocess.run(["pdfinfo", str(report)], check=True, capture_output=True, text=True).stdout
    page_match = re.search(r"^Pages:\s+(\d+)\s*$", pdf_info, flags=re.MULTILINE)
    report_pages = int(page_match.group(1)) if page_match else -1
    if report_pages != 10:
        raise AssertionError(f"Expected a 10-page report: {report}")
    overlays = sorted(overlay_dir.glob("*.png"))
    if len(overlays) != expected_cells:
        raise AssertionError(f"Expected {expected_cells} overlays, observed {len(overlays)}")
    table = pd.read_csv(output_dir / f"{timing_ms}ms_selected_cell_overlays.csv")
    if len(table) != expected_cells or not (table.plotted_trajectories == table.retained_trajectories).all():
        raise AssertionError("Overlay table count validation failed")
    if not (table.padding_px == int(cfg["cell_bbox_padding_px"])).all():
        raise AssertionError("Overlay padding validation failed")
    for row in table.itertuples(index=False):
        path = Path(row.output)
        if path.parent != overlay_dir or not path.exists():
            raise AssertionError(f"Overlay path is missing or routed incorrectly: {path}")
        if int(path.name.split("_", 1)[0]) != int(row.retained_trajectories):
            raise AssertionError(f"Overlay filename count mismatch: {path.name}")
    curves = pd.read_csv(output_dir / f"{timing_ms}ms_adjacent_step_size_curves.csv")
    expected_steps = {30: 366678, 100: 179828}[timing_ms]
    observed_steps = int(curves.groupby(["timing_ms", "phase"]).steps.first().sum())
    if observed_steps != expected_steps:
        raise AssertionError(f"Expected {expected_steps} adjacent steps, observed {observed_steps}")
    for phase in PHASE_ORDER:
        cdf = curves[(curves.phase == phase) & (curves.curve == "cdf")].sort_values("x_nm")
        if cdf.empty or np.any(np.diff(cdf.y) < -1e-12) or not np.isclose(cdf.y.iloc[-1], 1):
            raise AssertionError(f"Invalid adjacent-step CDF for {phase}")
    mean_cdf = pd.read_csv(output_dir / f"{timing_ms}ms_mean_stepsize_nm_trajectory_cdf.csv")
    for phase in PHASE_ORDER:
        current = mean_cdf[mean_cdf.phase == phase].sort_values("x_nm")
        if (current.empty or np.any(np.diff(current.cumulative_probability) < -1e-12)
                or not np.isclose(current.cumulative_probability.iloc[-1], 1)):
            raise AssertionError(f"Invalid trajectory mean-step-size CDF for {phase}")
    return {"png_plots": len(pngs), "svg_plots": len(svgs), "report_pages": 10,
            "overlays": len(overlays), "adjacent_steps": observed_steps}


def remove_legacy_timing_outputs(cfg: dict[str, Any], timing_ms: int) -> list[str]:
    """Remove old flat artifacts only after their reorganized replacements validate."""
    root = Path(cfg["result_root"])
    removed: list[str] = []
    for path in sorted(root.glob(f"{timing_ms}ms_*")):
        if path.is_file():
            path.unlink()
            removed.append(str(path))
    flat_overlay_dir = root / cfg["output_folders"]["selected_cells"]
    for path in sorted(flat_overlay_dir.glob(f"*__{timing_ms}ms__*.png")):
        path.unlink()
        removed.append(str(path))
    return removed


def analyze_timing_stage(cfg: dict[str, Any], timing_ms: int, *, resume: bool, force: bool) -> dict[str, Any]:
    if timing_ms not in {30, 100}:
        raise ValueError("timing_ms must be 30 or 100")
    root = Path(cfg["result_root"])
    output_dir = comparison_dir(cfg, timing_ms)
    output_dir.mkdir(parents=True, exist_ok=True)
    inventory_path = root / "cell_selection_summary.csv"
    inventory = pd.read_csv(inventory_path) if inventory_path.exists() else inventory_stage(cfg)
    expected = int(((inventory.timing_ms == timing_ms) & inventory.selected).sum())
    report_path = output_dir / f"{timing_ms}ms_selected_cell_diffusion_report.pdf"
    overlay_table_path = output_dir / f"{timing_ms}ms_selected_cell_overlays.csv"
    if resume and report_path.exists() and overlay_table_path.exists() and not force:
        overlay_table = pd.read_csv(overlay_table_path)
        if len(overlay_table) == expected and all(Path(path).exists() for path in overlay_table.output):
            validation = validate_timing_outputs(cfg, timing_ms, expected)
            removed_legacy = remove_legacy_timing_outputs(cfg, timing_ms)
            _update_stage_status(root, f"analyze_{timing_ms}ms", {
                "status": "complete", "revalidated_unix": time.time(), "selected_cells": expected,
                "report": str(report_path), "validation": validation,
            })
            print(f"{timing_ms} ms stage already complete; skipping", flush=True)
            return {"status": "skipped", "selected_cells": expected, "report": str(report_path),
                    "validation": validation, "removed_legacy_artifacts": removed_legacy}
    _update_stage_status(root, f"analyze_{timing_ms}ms", {"status": "running", "started_unix": time.time()})
    overlays = render_selected_cells(cfg, inventory, timing_ms, resume=resume, force=force)
    aio = load_selected_aio(cfg, inventory, timing_ms)
    metrics = cell_metrics(aio, cfg)
    summary = grouped_metric_summary(metrics, cfg)
    omnibus, pairwise = phase_statistics(metrics, cfg)
    atomic_csv(aio, output_dir / f"{timing_ms}ms_selected_trajectories_AIO.csv")
    atomic_csv(metrics, output_dir / f"{timing_ms}ms_selected_cell_metrics.csv")
    atomic_csv(summary, output_dir / f"{timing_ms}ms_selected_cell_metric_summary.csv")
    atomic_csv(omnibus, output_dir / f"{timing_ms}ms_phase_omnibus_statistics.csv")
    atomic_csv(pairwise, output_dir / f"{timing_ms}ms_phase_pairwise_statistics.csv")
    step_curves = step_distribution_curves(aio, cfg, cohort="all_selected_cells")
    atomic_csv(step_curves, output_dir / f"{timing_ms}ms_adjacent_step_size_curves.csv")
    mean_step_metric = next(item for item in DEEP_PLOT_METRICS if item[0] == "mean_stepsize_nm")
    mean_step_cdf = trajectory_metric_cdf(
        aio, "mean_stepsize_nm", cfg, mean_step_metric[2],
        grid_points=int(cfg["step_size_distribution"]["cdf_grid_points"]),
    )
    atomic_csv(mean_step_cdf, output_dir / f"{timing_ms}ms_mean_stepsize_nm_trajectory_cdf.csv")
    saspt_all, saspt_balanced = _run_selected_saspt(cfg, aio, timing_ms)
    plots = _plot_timing_report(cfg, timing_ms, inventory, aio, metrics, saspt_all, saspt_balanced, mean_step_cdf)
    if len(overlays) != expected or not (overlays.plotted_trajectories == overlays.retained_trajectories).all():
        raise AssertionError(f"{timing_ms} ms selected-cell overlay validation failed")
    removed_obsolete_plots = remove_obsolete_timing_plots(cfg, timing_ms)
    validation = validate_timing_outputs(cfg, timing_ms, expected)
    removed_legacy = remove_legacy_timing_outputs(cfg, timing_ms)
    result = {
        "status": "complete", "finished_unix": time.time(), "selected_cells": expected,
        "overlays": len(overlays), "aio_trajectories": len(aio), "plots": plots,
        "validation": validation, "removed_obsolete_plots": removed_obsolete_plots,
        "removed_legacy_artifacts": removed_legacy,
    }
    _update_stage_status(root, f"analyze_{timing_ms}ms", result)
    return result


def refresh_timing_report(cfg: dict[str, Any], timing_ms: int) -> dict[str, Any]:
    """Rebuild timing plots/reports from existing tables without rerendering overlays."""
    root = Path(cfg["result_root"])
    output_dir = comparison_dir(cfg, timing_ms)
    inventory = pd.read_csv(root / "cell_selection_summary.csv")
    expected = int(((inventory.timing_ms == timing_ms) & inventory.selected).sum())
    aio = pd.read_csv(output_dir / f"{timing_ms}ms_selected_trajectories_AIO.csv")
    metrics = pd.read_csv(output_dir / f"{timing_ms}ms_selected_cell_metrics.csv")
    saspt_all = pd.read_csv(output_dir / f"{timing_ms}ms_selected_cells_saSPT_all.csv")
    saspt_balanced = pd.read_csv(output_dir / f"{timing_ms}ms_selected_cells_saSPT_balanced.csv")
    mean_step_metric = next(item for item in DEEP_PLOT_METRICS if item[0] == "mean_stepsize_nm")
    mean_step_cdf = trajectory_metric_cdf(
        aio, "mean_stepsize_nm", cfg, mean_step_metric[2],
        grid_points=int(cfg["step_size_distribution"]["cdf_grid_points"]),
    )
    atomic_csv(mean_step_cdf, output_dir / f"{timing_ms}ms_mean_stepsize_nm_trajectory_cdf.csv")
    plots = _plot_timing_report(
        cfg, timing_ms, inventory, aio, metrics, saspt_all, saspt_balanced, mean_step_cdf,
    )
    removed_obsolete_plots = remove_obsolete_timing_plots(cfg, timing_ms)
    validation = validate_timing_outputs(cfg, timing_ms, expected)
    result = {
        "status": "complete", "finished_unix": time.time(), "selected_cells": expected,
        "aio_trajectories": len(aio), "plots": plots, "validation": validation,
        "removed_obsolete_plots": removed_obsolete_plots, "overlays_rerendered": False,
    }
    _update_stage_status(root, f"analyze_{timing_ms}ms", result)
    return result


def _timing_counterpart_fov(fov: str, target_ms: int) -> str:
    return re.sub(r"-(30|100)ms-", f"-{target_ms}ms-", fov)


def match_timing_cells(cfg: dict[str, Any], inventory: pd.DataFrame) -> pd.DataFrame:
    root = Path(cfg["result_root"])
    minimum_iou = float(cfg["pairing_minimum_iou"])
    records: list[dict[str, Any]] = []
    for phase in ("FVP 2 h", "Recovery 2 h"):
        left_cells = inventory[(inventory.phase == phase) & (inventory.timing_ms == 30)]
        for fov30, cells30 in left_cells.groupby("fov", sort=True):
            fov100 = _timing_counterpart_fov(fov30, 100)
            cells100 = inventory[(inventory.fov == fov100) & (inventory.phase == phase)]
            if cells100.empty:
                for row in cells30.itertuples(index=False):
                    records.append({
                        "phase": phase, "fov_30ms": fov30, "fov_100ms": fov100,
                        "cell_id_30ms": int(row.cell_id), "cell_id_100ms": np.nan,
                        "mask_iou": np.nan, "retained_30ms": int(row.retained_trajectories),
                        "retained_100ms": np.nan, "both_pass_cutoff": False,
                        "accepted_pair": False, "rejection_reason": "missing_100ms_fov",
                    })
                continue
            mask30 = np.asarray(imread(_mask_path(root, fov30)))
            mask100 = np.asarray(imread(_mask_path(root, fov100)))
            ids30 = cells30.cell_id.astype(int).to_numpy()
            ids100 = cells100.cell_id.astype(int).to_numpy()
            scores = np.zeros((len(ids30), len(ids100)), dtype=float)
            for i, cell30 in enumerate(ids30):
                left = mask30 == cell30
                for j, cell100 in enumerate(ids100):
                    right = mask100 == cell100
                    intersection = np.sum(left & right)
                    if intersection:
                        scores[i, j] = intersection / np.sum(left | right)
            rows, columns = linear_sum_assignment(-scores)
            assigned30, assigned100 = set(), set()
            for i, j in zip(rows, columns):
                assigned30.add(int(ids30[i])); assigned100.add(int(ids100[j]))
                row30 = cells30[cells30.cell_id == ids30[i]].iloc[0]
                row100 = cells100[cells100.cell_id == ids100[j]].iloc[0]
                passes = bool(row30.selected and row100.selected)
                accepted = passes and scores[i, j] >= minimum_iou
                reason = "" if accepted else "below_iou" if passes else "cell_below_trajectory_cutoff"
                records.append({
                    "phase": phase, "fov_30ms": fov30, "fov_100ms": fov100,
                    "cell_id_30ms": int(ids30[i]), "cell_id_100ms": int(ids100[j]),
                    "mask_iou": float(scores[i, j]),
                    "retained_30ms": int(row30.retained_trajectories),
                    "retained_100ms": int(row100.retained_trajectories),
                    "both_pass_cutoff": passes, "accepted_pair": accepted, "rejection_reason": reason,
                })
            for cell30 in set(ids30) - assigned30:
                row30 = cells30[cells30.cell_id == cell30].iloc[0]
                records.append({
                    "phase": phase, "fov_30ms": fov30, "fov_100ms": fov100,
                    "cell_id_30ms": int(cell30), "cell_id_100ms": np.nan, "mask_iou": np.nan,
                    "retained_30ms": int(row30.retained_trajectories), "retained_100ms": np.nan,
                    "both_pass_cutoff": False, "accepted_pair": False, "rejection_reason": "unmatched_30ms_cell",
                })
            for cell100 in set(ids100) - assigned100:
                row100 = cells100[cells100.cell_id == cell100].iloc[0]
                records.append({
                    "phase": phase, "fov_30ms": fov30, "fov_100ms": fov100,
                    "cell_id_30ms": np.nan, "cell_id_100ms": int(cell100), "mask_iou": np.nan,
                    "retained_30ms": np.nan, "retained_100ms": int(row100.retained_trajectories),
                    "both_pass_cutoff": False, "accepted_pair": False, "rejection_reason": "unmatched_100ms_cell",
                })
    return pd.DataFrame.from_records(records)


def _paired_bootstrap_difference(left: np.ndarray, right: np.ndarray, iterations: int,
                                 rng: np.random.Generator) -> tuple[float, float]:
    valid = np.isfinite(left) & np.isfinite(right)
    differences = right[valid] - left[valid]
    if not len(differences):
        return np.nan, np.nan
    estimates = [np.median(rng.choice(differences, len(differences), replace=True)) for _ in range(iterations)]
    return tuple(map(float, np.quantile(estimates, [0.025, 0.975])))


def timing_comparison_aio(aio30: pd.DataFrame, aio100: pd.DataFrame,
                          matches: pd.DataFrame) -> pd.DataFrame:
    """Use independent Before cells and accepted matched FVP/Recovery cells."""
    parts = [aio30[aio30.phase == "Before"].copy(), aio100[aio100.phase == "Before"].copy()]
    accepted = matches[matches.accepted_pair]
    for phase in ("FVP 2 h", "Recovery 2 h"):
        current = accepted[accepted.phase == phase]
        keys30 = pd.MultiIndex.from_frame(current[["fov_30ms", "cell_id_30ms"]].rename(
            columns={"fov_30ms": "fov", "cell_id_30ms": "cell_id"}))
        keys100 = pd.MultiIndex.from_frame(current[["fov_100ms", "cell_id_100ms"]].rename(
            columns={"fov_100ms": "fov", "cell_id_100ms": "cell_id"}))
        index30 = pd.MultiIndex.from_frame(aio30[["fov", "cell_id"]])
        index100 = pd.MultiIndex.from_frame(aio100[["fov", "cell_id"]])
        parts.append(aio30[(aio30.phase == phase) & index30.isin(keys30)].copy())
        parts.append(aio100[(aio100.phase == phase) & index100.isin(keys100)].copy())
    return pd.concat(parts, ignore_index=True, sort=False)


def timing_comparison_plot_stems() -> list[str]:
    return ["30ms_vs_100ms_01_pairing_diagnostics"] + [
        f"30ms_vs_100ms_{index:02d}_{metric}" for index, metric in enumerate(PAIR_METRICS, 2)
    ] + [
        "30ms_vs_100ms_09_adjacent_step_size_histogram",
        "30ms_vs_100ms_10_adjacent_step_size_cdf",
    ]


def validate_timing_comparison_outputs(cfg: dict[str, Any]) -> dict[str, int]:
    output_dir = timing_comparison_dir(cfg)
    expected_stems = set(timing_comparison_plot_stems())
    pngs = sorted(output_dir.glob("30ms_vs_100ms_*.png"))
    svgs = sorted(output_dir.glob("30ms_vs_100ms_*.svg"))
    report = output_dir / "30ms_vs_100ms_selected_cell_comparison_report.pdf"
    if {path.stem for path in pngs} != expected_stems or {path.stem for path in svgs} != expected_stems:
        raise AssertionError("Timing-comparison PNG/SVG plot set is incomplete")
    for svg in svgs:
        if ET.parse(svg).getroot().tag.rsplit("}", 1)[-1] != "svg":
            raise AssertionError(f"Invalid SVG: {svg}")
    pdfs = sorted(output_dir.glob("*.pdf"))
    if pdfs != [report]:
        raise AssertionError("Only the timing-comparison report PDF may remain")
    pdf_info = subprocess.run(["pdfinfo", str(report)], check=True, capture_output=True, text=True).stdout
    page_match = re.search(r"^Pages:\s+(\d+)\s*$", pdf_info, flags=re.MULTILINE)
    if not page_match or int(page_match.group(1)) != 10:
        raise AssertionError("Timing-comparison report must contain 10 pages")
    expected_csvs = {
        "30ms_vs_100ms_cell_matching.csv", "30ms_vs_100ms_accepted_cell_pairs.csv",
        "30ms_vs_100ms_statistics.csv", "30ms_vs_100ms_adjacent_step_size_curves.csv",
    }
    if {path.name for path in output_dir.glob("*.csv")} != expected_csvs:
        raise AssertionError("Timing-comparison CSV set is incomplete")
    matches = pd.read_csv(output_dir / "30ms_vs_100ms_cell_matching.csv")
    accepted = matches[matches.accepted_pair]
    counts = accepted.groupby("phase").size().to_dict()
    if counts != {"FVP 2 h": 32, "Recovery 2 h": 38}:
        raise AssertionError(f"Unexpected accepted-pair counts: {counts}")
    curves = pd.read_csv(output_dir / "30ms_vs_100ms_adjacent_step_size_curves.csv")
    for (timing, phase), group in curves[curves.curve == "cdf"].groupby(["timing_ms", "phase"]):
        ordered = group.sort_values("x_nm").y.to_numpy(float)
        if np.any(np.diff(ordered) < -1e-12) or not np.isclose(ordered[-1], 1):
            raise AssertionError(f"Invalid timing CDF for {timing}/{phase}")
    return {"png_plots": len(pngs), "svg_plots": len(svgs), "report_pages": 10,
            "accepted_fvp_pairs": 32, "accepted_recovery_pairs": 38}


def compare_timings_stage(cfg: dict[str, Any], *, resume: bool, force: bool) -> dict[str, Any]:
    root = Path(cfg["result_root"])
    output_dir = timing_comparison_dir(cfg)
    required = [comparison_dir(cfg, timing) / f"{timing}ms_selected_cell_metrics.csv" for timing in (30, 100)]
    if not all(path.exists() for path in required):
        raise RuntimeError("Complete and review both timing stages before compare-timings")
    report_path = output_dir / "30ms_vs_100ms_selected_cell_comparison_report.pdf"
    if resume and report_path.exists() and not force:
        validation = validate_timing_comparison_outputs(cfg)
        return {"status": "skipped", "report": str(report_path), "validation": validation}
    _update_stage_status(root, "compare_timings", {"status": "running", "started_unix": time.time()})
    output_dir.mkdir(parents=True, exist_ok=True)
    inventory = pd.read_csv(root / "cell_selection_summary.csv")
    metrics30 = pd.read_csv(required[0])
    metrics100 = pd.read_csv(required[1])
    aio30 = pd.read_csv(comparison_dir(cfg, 30) / "30ms_selected_trajectories_AIO.csv", low_memory=False)
    aio100 = pd.read_csv(comparison_dir(cfg, 100) / "100ms_selected_trajectories_AIO.csv", low_memory=False)
    matches = match_timing_cells(cfg, inventory)
    atomic_csv(matches, output_dir / "30ms_vs_100ms_cell_matching.csv")
    accepted = matches[matches.accepted_pair].copy()
    comparison_aio = timing_comparison_aio(aio30, aio100, matches)
    step_curves = step_distribution_curves(
        comparison_aio, cfg, cohort="matched_fvp_recovery_independent_before",
    )
    atomic_csv(step_curves, output_dir / "30ms_vs_100ms_adjacent_step_size_curves.csv")
    left = metrics30.rename(columns={column: f"{column}_30ms" for column in PAIR_METRICS})
    right = metrics100.rename(columns={column: f"{column}_100ms" for column in PAIR_METRICS})
    paired = accepted.merge(
        left, left_on=["fov_30ms", "cell_id_30ms"], right_on=["fov", "cell_id"], how="left"
    ).drop(columns=["fov", "cell_id"])
    paired = paired.merge(
        right, left_on=["fov_100ms", "cell_id_100ms"], right_on=["fov", "cell_id"], how="left",
        suffixes=("", "_right"),
    ).drop(columns=["fov", "cell_id"])
    atomic_csv(paired, output_dir / "30ms_vs_100ms_accepted_cell_pairs.csv")

    rng = np.random.default_rng(int(cfg["random_seed"]) + 200)
    stat_rows = []
    for phase in PHASE_ORDER:
        for metric in PAIR_METRICS:
            if phase == "Before":
                a = metrics30[metrics30.phase == phase][metric].dropna().to_numpy(float)
                b = metrics100[metrics100.phase == phase][metric].dropna().to_numpy(float)
                statistic, pvalue = stats.mannwhitneyu(a, b, alternative="two-sided") if len(a) and len(b) else (np.nan, np.nan)
                stat_rows.append({
                    "phase": phase, "metric": metric, "paired": False,
                    "test": "Mann-Whitney U; Before fields are independent", "n_30ms": len(a), "n_100ms": len(b),
                    "statistic": statistic, "p_raw": pvalue,
                    "median_difference_100ms_minus_30ms": float(np.median(b) - np.median(a)) if len(a) and len(b) else np.nan,
                    "difference_ci_low": np.nan, "difference_ci_high": np.nan,
                })
            else:
                current = paired[paired.phase == phase]
                a = current[f"{metric}_30ms"].to_numpy(float)
                b = current[f"{metric}_100ms"].to_numpy(float)
                valid = np.isfinite(a) & np.isfinite(b)
                statistic, pvalue = stats.wilcoxon(a[valid], b[valid]) if valid.sum() else (np.nan, np.nan)
                low, high = _paired_bootstrap_difference(a, b, int(cfg["bootstrap_iterations"]), rng)
                stat_rows.append({
                    "phase": phase, "metric": metric, "paired": True,
                    "test": "Wilcoxon signed-rank on IoU-matched cells", "n_30ms": int(valid.sum()),
                    "n_100ms": int(valid.sum()), "statistic": statistic, "p_raw": pvalue,
                    "median_difference_100ms_minus_30ms": float(np.nanmedian(b[valid] - a[valid])) if valid.sum() else np.nan,
                    "difference_ci_low": low, "difference_ci_high": high,
                })
    statistics = pd.DataFrame.from_records(stat_rows)
    statistics["p_holm_within_phase"] = statistics.groupby("phase").p_raw.transform(holm_adjust)
    atomic_csv(statistics, output_dir / "30ms_vs_100ms_statistics.csv")

    cache = root / ".cache" / "matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn as sns
    _style()
    outputs = []
    with PdfPages(report_path) as pdf:
        fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.5))
        sns.histplot(data=matches[matches.both_pass_cutoff], x="mask_iou", hue="phase", hue_order=PHASE_ORDER[1:],
                     palette=PHASE_COLORS, bins=20, element="step", fill=False, lw=2, ax=axes[0])
        axes[0].axvline(float(cfg["pairing_minimum_iou"]), color="black", ls="--")
        _format_axis(axes[0], "Cell matches")
        axes[0].set_title("Mask-overlap matching")
        pair_counts = accepted.groupby("phase").size().reindex(PHASE_ORDER[1:]).fillna(0)
        axes[1].bar(pair_counts.index, pair_counts.values, color=[PHASE_COLORS[p] for p in pair_counts.index])
        axes[1].set_title("Accepted 30/100 ms pairs")
        _format_axis(axes[1], "Cell pairs")
        fig.tight_layout()
        outputs += _save(fig, output_dir, "30ms_vs_100ms_01_pairing_diagnostics", pdf)
        plt.close(fig)

        for index, metric in enumerate(PAIR_METRICS, 2):
            fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.3), sharey=True)
            for ax, phase in zip(axes, PHASE_ORDER):
                if phase == "Before":
                    values = pd.concat([
                        metrics30[metrics30.phase == phase][[metric]].assign(timing="30 ms"),
                        metrics100[metrics100.phase == phase][[metric]].assign(timing="100 ms"),
                    ])
                    sns.boxplot(data=values, x="timing", y=metric, color="white", showfliers=False, ax=ax)
                    sns.stripplot(data=values, x="timing", y=metric, color=PHASE_COLORS[phase], size=2, alpha=0.55, ax=ax)
                    n30 = len(metrics30[metrics30.phase == phase])
                    n100 = len(metrics100[metrics100.phase == phase])
                    tracks30 = int(metrics30[metrics30.phase == phase].aio_valid_trajectories.sum())
                    tracks100 = int(metrics100[metrics100.phase == phase].aio_valid_trajectories.sum())
                    ax.set_title(
                        f"Before · independent\n30 ms: N={n30} cells, n={tracks30:,} trajectories\n"
                        f"100 ms: N={n100} cells, n={tracks100:,} trajectories", fontsize=8,
                    )
                else:
                    current = paired[paired.phase == phase]
                    for row in current.itertuples():
                        y30, y100 = getattr(row, f"{metric}_30ms"), getattr(row, f"{metric}_100ms")
                        if np.isfinite(y30) and np.isfinite(y100):
                            ax.plot([0, 1], [y30, y100], color=PHASE_COLORS[phase], alpha=0.35, lw=0.6)
                            ax.scatter([0, 1], [y30, y100], color=PHASE_COLORS[phase], s=7)
                    ax.set_xticks([0, 1], ["30 ms", "100 ms"])
                    ax.set_title(f"{phase} · matched\nN={len(current)} cells")
                ax.set_xlabel("")
                _format_axis(ax, metric if ax is axes[0] else None)
            fig.suptitle(metric)
            fig.tight_layout()
            outputs += _save(fig, output_dir, f"30ms_vs_100ms_{index:02d}_{metric}", pdf)
            plt.close(fig)

        for index, curve_name, ylabel, title in (
            (9, "histogram", "Probability", "Adjacent-frame step-size distribution"),
            (10, "cdf", "Cumulative probability", "Adjacent-frame step-size CDF"),
        ):
            fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.3), sharey=True)
            for ax, phase in zip(axes, PHASE_ORDER):
                for timing, color in ((30, "#6a3d9a"), (100, "#1b9e77")):
                    current = step_curves[
                        (step_curves.phase == phase) & (step_curves.timing_ms == timing) &
                        (step_curves.curve == curve_name)
                    ]
                    if current.empty:
                        continue
                    info = current.iloc[0]
                    ax.step(current.x_nm, current.y,
                            where="mid" if curve_name == "histogram" else "post",
                            color=color, lw=2,
                            label=f"{timing} ms\nN={int(info.cells):,} cells\n"
                                  f"n={int(info.trajectories):,} trajectories")
                ax.set_xlim(float(cfg["step_size_distribution"]["minimum_nm"]),
                            float(cfg["step_size_distribution"]["maximum_nm"]))
                if curve_name == "cdf":
                    ax.set_ylim(0, 1)
                ax.set_title(f"{phase}\n{'independent' if phase == 'Before' else 'matched cells'}")
                ax.set_xlabel("Step size, nm")
                _format_axis(ax, ylabel if ax is axes[0] else None)
                ax.legend(frameon=False, fontsize=6.5, loc="lower right" if curve_name == "cdf" else "upper right")
            fig.suptitle(f"30 ms versus 100 ms · {title}")
            fig.tight_layout()
            outputs += _save(fig, output_dir, f"30ms_vs_100ms_{index:02d}_adjacent_step_size_{curve_name}", pdf)
            plt.close(fig)
    outputs.append(str(report_path))
    validation = validate_timing_comparison_outputs(cfg)
    result = {
        "status": "complete", "finished_unix": time.time(),
        "accepted_FVP_pairs": int(((accepted.phase == "FVP 2 h")).sum()),
        "accepted_recovery_pairs": int(((accepted.phase == "Recovery 2 h")).sum()),
        "report": str(report_path), "plots": outputs, "validation": validation,
    }
    _update_stage_status(root, "compare_timings", result)
    return result


def validate_inventory(inventory: pd.DataFrame, cfg: dict[str, Any]) -> None:
    if len(inventory) != 388:
        raise AssertionError(f"Expected 388 segmented cells, observed {len(inventory)}")
    if int(inventory.selected.sum()) != 269:
        raise AssertionError(f"Expected 269 selected cells, observed {int(inventory.selected.sum())}")
    expected = {(30, "Before"): 57, (30, "FVP 2 h"): 58, (30, "Recovery 2 h"): 50,
                (100, "Before"): 29, (100, "FVP 2 h"): 34, (100, "Recovery 2 h"): 41}
    observed = inventory[inventory.selected].groupby(["timing_ms", "phase"]).size().to_dict()
    if observed != expected:
        raise AssertionError(f"Unexpected selected-cell groups: {observed}")
    if int(cfg["trajectory_cutoff"]) != 50:
        raise AssertionError("Dataset acceptance counts are defined for a 50-trajectory cutoff")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=(
        "inventory", "analyze-30ms", "analyze-100ms", "refresh-30ms-report",
        "refresh-100ms-report", "compare-timings", "run-all",
    ))
    parser.add_argument("--config", required=True)
    policy = parser.add_mutually_exclusive_group()
    policy.add_argument("--resume", action="store_true")
    policy.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    cfg = load_deep_config(args.config)
    write_analysis_metadata(cfg)
    inventory_path = Path(cfg["result_root"]) / "cell_selection_summary.csv"
    inventory = inventory_stage(cfg) if args.action == "inventory" or not inventory_path.exists() else pd.read_csv(inventory_path)
    validate_inventory(inventory, cfg)
    if args.action == "inventory":
        print(selection_counts(inventory).to_string(index=False))
        return
    if args.action == "refresh-30ms-report":
        print(json.dumps(refresh_timing_report(cfg, 30), indent=2, default=str))
        return
    if args.action == "refresh-100ms-report":
        print(json.dumps(refresh_timing_report(cfg, 100), indent=2, default=str))
        return
    if args.action in {"analyze-30ms", "run-all"}:
        print(json.dumps(analyze_timing_stage(cfg, 30, resume=args.resume, force=args.force), indent=2, default=str))
    if args.action in {"analyze-100ms", "run-all"}:
        print(json.dumps(analyze_timing_stage(cfg, 100, resume=args.resume, force=args.force), indent=2, default=str))
    if args.action in {"compare-timings", "run-all"}:
        print(json.dumps(compare_timings_stage(cfg, resume=args.resume, force=args.force), indent=2, default=str))


if __name__ == "__main__":
    main()
