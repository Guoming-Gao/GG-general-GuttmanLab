"""Step 5: link Spotiflow detections with LapTrack, assign cells, and render trajectory QC."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from skimage.segmentation import find_boundaries
from tifffile import TiffFile, imread

from spt_shared import atomic_csv, output_action, output_dirs, stage_timer
from spt_video import render_overlay_video
from step01_inspect_inputs import accepted_rows, inspect_inputs
from step02_preprocess_videos import artifact_paths, channel_specs
from step03_segment_cells import segmentation_paths
from step04_detect_spots import detection_paths


LEGACY_COLUMNS = [
    "spotID", "trackID", "QUALITY", "x", "y", "t", "R", "meanIntensity",
    "medianIntensity", "minIntensity", "maxIntensity", "totalIntensity",
    "stdIntensity", "contrast", "SNR",
]


def tracking_paths(output_dir: str | Path, fov: str, channel: str) -> dict[str, Path]:
    dirs = output_dirs(output_dir); root = dirs["trajectories"]
    return {
        "canonical": root / f"{fov}__{channel}__trajectories.csv",
        "legacy": root / f"{fov}__{channel}__tracks_TrackMate_compatible.csv",
        "summary": root / f"{fov}__{channel}__track_summary.csv",
        "edges": root / f"{fov}__{channel}__trajectory_edges.csv",
        "overlay_png": dirs["qc"] / "trajectory_overlays" / f"{fov}__{channel}__trajectory_overlay_1200dpi.png",
        "overlay_pdf": dirs["qc"] / "trajectory_overlays" / f"{fov}__{channel}__trajectory_overlay_vector.pdf",
        "video": dirs["qc"] / "trajectory_videos" / f"{fov}__{channel}__trajectory_QC.mp4",
        "dashboard": dirs["qc"] / "secondary_metrics" / f"{fov}__{channel}__tracking_metrics.png",
    }


def assign_track_labels(track: pd.DataFrame) -> dict[str, Any]:
    labels = track["detection_cell_id"].fillna(0).astype(int).to_numpy()
    nonzero = labels[labels > 0]
    if nonzero.size:
        values, counts = np.unique(nonzero, return_counts=True); cell_id = int(values[np.argmax(counts)])
    else: cell_id = 0
    confidence = float(np.mean(labels == cell_id)) if labels.size else 0.0
    encountered = sorted(set(map(int, labels)))
    ambiguous = len([value for value in encountered if value > 0]) > 1 or confidence < 1
    return {"cell_id": cell_id, "assignment_confidence": confidence,
            "encountered_cell_ids": ";".join(map(str, encountered)), "assignment_ambiguous": ambiguous}


def _edge_table(track: pd.DataFrame, trajectory_uid: str) -> pd.DataFrame:
    ordered = track.sort_values("frame")
    rows = []
    values = list(ordered.itertuples())
    for source, target in zip(values[:-1], values[1:]):
        frame_difference = int(target.frame - source.frame)
        rows.append({
            "dataset": source.dataset, "condition": source.condition, "fov": source.fov,
            "channel": source.channel, "trackID": int(source.trackID), "trajectory_uid": trajectory_uid,
            "source_spotID": int(source.spotID), "target_spotID": int(target.spotID),
            "source_frame": int(source.frame), "target_frame": int(target.frame),
            "frame_difference": frame_difference,
            "distance_px": float(np.hypot(target.x-source.x, target.y-source.y)),
            "edge_type": "adjacent" if frame_difference == 1 else "gap_closing",
        })
    return pd.DataFrame.from_records(rows)


def link_detections(detections: pd.DataFrame, cfg: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if detections.empty:
        empty = detections.copy(); empty["trackID"] = pd.Series(dtype=int)
        return empty, pd.DataFrame(), pd.DataFrame()
    detections = detections.copy()
    if "condition" not in detections:
        detections["condition"] = ""
    from laptrack import LapTrack
    settings = cfg["laptrack"]; distance = float(settings["max_link_distance_px"])
    tracker = LapTrack(
        metric=settings.get("metric", "sqeuclidean"), cutoff=distance**2,
        gap_closing_cutoff=distance**2 if settings.get("gap_closing", True) else False,
        gap_closing_max_frame_count=int(settings.get("gap_closing_max_frame_count", 2)),
        splitting_cutoff=False, merging_cutoff=False,
    )
    tracked, _, _ = tracker.predict_dataframe(
        detections.reset_index(drop=True), coordinate_cols=["x", "y"], frame_col="frame", only_coordinate_cols=False,
    )
    tracked = tracked.rename(columns={"track_id": "trackID"}); tracked["trackID"] = tracked["trackID"].astype(int)
    lengths = tracked.groupby("trackID").size(); keep = lengths[lengths >= int(settings["min_track_length"])].index
    tracked = tracked[tracked["trackID"].isin(keep)].copy()
    enriched, summaries, edge_tables = [], [], []
    for track_id, group in tracked.groupby("trackID", sort=True):
        assignment = assign_track_labels(group); group = group.copy()
        for key, value in assignment.items(): group[key] = value
        uid = f"{cfg['dataset']}/{group['fov'].iloc[0]}/{group['channel'].iloc[0]}/{int(track_id)}"
        group["trajectory_uid"] = uid; edges = _edge_table(group, uid)
        maximum = float(edges["distance_px"].max()) if not edges.empty else 0.0
        summaries.append({
            "dataset": cfg["dataset"], "condition": group["condition"].iloc[0], "fov": group["fov"].iloc[0],
            "channel": group["channel"].iloc[0], "trackID": int(track_id), "trajectory_uid": uid,
            "n_detections": len(group), "first_frame": int(group["frame"].min()), "last_frame": int(group["frame"].max()),
            "max_step_px": maximum,
            "gap_closing_edges": int((edges["edge_type"] == "gap_closing").sum()) if not edges.empty else 0,
            **assignment,
        })
        enriched.append(group); edge_tables.append(edges)
    canonical = pd.concat(enriched, ignore_index=True) if enriched else tracked
    canonical["t"] = canonical["frame"].astype(int)
    edges = pd.concat(edge_tables, ignore_index=True) if edge_tables else pd.DataFrame()
    return canonical.sort_values(["trackID", "frame"]), pd.DataFrame.from_records(summaries), edges


def legacy_table(canonical: pd.DataFrame) -> pd.DataFrame:
    if canonical.empty: return pd.DataFrame(columns=LEGACY_COLUMNS)
    result = canonical.copy(); result["QUALITY"] = result["spotiflow_probability"]
    return result.reindex(columns=LEGACY_COLUMNS).astype(float)


def _dog_mip(path: Path) -> np.ndarray:
    result = None
    with TiffFile(path) as tif:
        for page in tif.pages:
            positive = np.maximum(page.asarray().astype(np.float32), 0)
            result = positive.copy() if result is None else np.maximum(result, positive)
    return result


def save_trajectory_overlay(dog_path: Path, mask: np.ndarray, tracks: pd.DataFrame, edges: pd.DataFrame,
                            png_path: Path, pdf_path: Path, cfg: dict[str, Any], title: str) -> None:
    png_path.parent.mkdir(parents=True, exist_ok=True)
    cache = Path(cfg["output_dir"]) / ".cache" / "matplotlib"; cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mip = _dog_mip(dog_path); finite = mip[np.isfinite(mip)]
    upper = float(np.percentile(finite[finite > 0], 99.8)) if np.any(finite > 0) else 1.0
    height, width = mip.shape; figsize = (max(3, width/100), max(3, height/100))
    fig, ax = plt.subplots(figsize=figsize)
    ax.imshow(mip, cmap="gray", vmin=0, vmax=upper, interpolation="nearest")
    ax.contour(find_boundaries(mask), levels=[0.5], colors="cyan", linewidths=0.1)
    cmap = plt.get_cmap("turbo")
    for index, (track_id, group) in enumerate(tracks.groupby("trackID")):
        ordered = group.sort_values("frame"); color = cmap((index * 0.61803398875) % 1)
        values = list(ordered.itertuples())
        for source, target in zip(values[:-1], values[1:]):
            is_gap = int(target.frame-source.frame) > 1
            ax.plot([source.x, target.x], [source.y, target.y], color="magenta" if is_gap else color,
                    linewidth=float(cfg["qc"]["track_linewidth_pt"]), linestyle="--" if is_gap else "-")
    ax.set_title(title, fontsize=7); ax.axis("off"); fig.tight_layout(pad=0.05)
    fig.savefig(png_path, dpi=int(cfg["qc"]["track_overlay_dpi"]), bbox_inches="tight", pad_inches=0.01)
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight", pad_inches=0.01); plt.close(fig)


def save_tracking_dashboard(detections: pd.DataFrame, tracks: pd.DataFrame, summary: pd.DataFrame,
                            edges: pd.DataFrame, path: Path, title: str) -> None:
    """Save compact secondary metrics; overlays remain the primary track QC."""
    path.parent.mkdir(parents=True, exist_ok=True)
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(8, 6))
    counts = detections.groupby("frame").size()
    axes[0, 0].plot(counts.index, counts.values, color="black", lw=0.8)
    axes[0, 0].set(xlabel="Frame", ylabel="Detections")
    if not summary.empty:
        axes[0, 1].hist(summary["n_detections"], bins=30, histtype="step", lw=1.5)
    axes[0, 1].set(xlabel="Locations per retained track", ylabel="Tracks")
    if not edges.empty:
        for edge_type, group in edges.groupby("edge_type"):
            axes[1, 0].hist(group["distance_px"], bins=np.linspace(0, 5, 31), histtype="step", lw=1.5, label=edge_type)
        axes[1, 0].legend(frameon=False)
    axes[1, 0].set(xlabel="Linked distance, px", ylabel="Edges", xlim=(0, 5))
    if not summary.empty:
        assigned = int((summary["cell_id"] > 0).sum()); unassigned = len(summary) - assigned
        axes[1, 1].bar(["Assigned", "Background"], [assigned, unassigned], color=["#4c72b0", "#b0b0b0"])
    axes[1, 1].set(ylabel="Tracks")
    for ax in axes.flat:
        ax.tick_params(direction="in"); ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(title); fig.tight_layout(); fig.savefig(path, dpi=300, bbox_inches="tight"); plt.close(fig)


def track_file(row: dict[str, Any], cfg: dict[str, Any], *, resume: bool = False, force: bool = False) -> list[dict[str, str]]:
    mask = np.asarray(imread(segmentation_paths(cfg["output_dir"], row["fov"])["mask"])); results = []
    for spec in channel_specs(cfg):
        if spec["role"] != "spt": continue
        paths = tracking_paths(cfg["output_dir"], row["fov"], spec["name"])
        if output_action(list(paths.values()), resume=resume, force=force) == "skip":
            results.append({key: str(value) for key, value in paths.items()}); continue
        detections = pd.read_csv(detection_paths(cfg["output_dir"], row["fov"], spec["name"])["table"])
        canonical, summary, edges = link_detections(detections, cfg)
        atomic_csv(canonical, paths["canonical"]); atomic_csv(legacy_table(canonical), paths["legacy"])
        atomic_csv(summary, paths["summary"]); atomic_csv(edges, paths["edges"])
        dog_path = artifact_paths(cfg["output_dir"], row["fov"], spec, cfg)["filtered"]
        save_trajectory_overlay(dog_path, mask, canonical, edges, paths["overlay_png"], paths["overlay_pdf"], cfg,
                                f"{row['fov']} — {spec['name']} — {len(summary)} tracks")
        save_tracking_dashboard(detections, canonical, summary, edges, paths["dashboard"],
                                f"{row['fov']} — {spec['name']} secondary tracking metrics")
        qc, spot = cfg["qc"], cfg["spotiflow"]
        threshold = spot.get("probability_threshold"); threshold_label = "0.5 (model)" if threshold is None else f"{float(threshold):g}"
        render_overlay_video(
            dog_path, paths["video"], detections, tracks=canonical, mask=mask,
            frame_interval_s=float(row["frame_interval_s"]), pixel_size_um=float(row["pixel_size_um"]),
            model_label=f"Spotiflow {spot.get('pretrained_model', 'general')}", threshold_label=threshold_label,
            link_label=f"LapTrack: <5 px, gap ≤2 frames, ≥5 locations", maximum_frames=int(qc["max_video_frames"]),
            scale=int(qc["render_scale"]), fps=int(qc["video_fps"]), crf=int(qc["video_crf"]),
            tail_frames=int(qc["track_tail_frames"]), scale_bar_um=float(qc["scale_bar_um"]),
        )
        results.append({key: str(value) for key, value in paths.items()})
    return results


def track_stage(cfg: dict[str, Any], manifest: pd.DataFrame, *, max_files: int | None = None,
                resume: bool = False, force: bool = False) -> None:
    rows = accepted_rows(manifest, max_files)
    with stage_timer(cfg, "05_track", {"files": len(rows), **cfg["laptrack"]}):
        for index, row in enumerate(rows, 1):
            print(f"[track {index}/{len(rows)}] {row['filename']}", flush=True)
            track_file(row, cfg, resume=resume, force=force)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", required=True)
    parser.add_argument("--max-files", type=int); policy = parser.add_mutually_exclusive_group()
    policy.add_argument("--resume", action="store_true"); policy.add_argument("--force", action="store_true")
    args = parser.parse_args(argv); cfg, manifest = inspect_inputs(args.config)
    track_stage(cfg, manifest, max_files=args.max_files, resume=args.resume, force=args.force)


if __name__ == "__main__":
    main()
