from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from skimage.segmentation import find_boundaries
from tifffile import imread

from .detection import detection_path
from .preprocess import artifact_paths, channel_specs
from .segmentation import segmentation_paths
from .utils import atomic_csv, output_action


LEGACY_COLUMNS = [
    "spotID",
    "trackID",
    "QUALITY",
    "x",
    "y",
    "t",
    "R",
    "meanIntensity",
    "medianIntensity",
    "minIntensity",
    "maxIntensity",
    "totalIntensity",
    "stdIntensity",
    "contrast",
    "SNR",
]


def assign_track_labels(track: pd.DataFrame) -> dict[str, Any]:
    labels = track["detection_cell_id"].fillna(0).astype(int).to_numpy()
    nonzero = labels[labels > 0]
    if nonzero.size:
        values, counts = np.unique(nonzero, return_counts=True)
        cell_id = int(values[np.argmax(counts)])
    else:
        cell_id = 0
    confidence = float(np.mean(labels == cell_id)) if labels.size else 0.0
    encountered = sorted(set(int(value) for value in labels))
    nonzero_labels = [value for value in encountered if value > 0]
    ambiguous = len(nonzero_labels) > 1 or confidence < 1.0
    return {
        "cell_id": cell_id,
        "assignment_confidence": confidence,
        "encountered_cell_ids": ";".join(map(str, encountered)),
        "assignment_ambiguous": bool(ambiguous),
    }


def _max_step(track: pd.DataFrame) -> float:
    ordered = track.sort_values("frame")
    xy = ordered[["x", "y"]].to_numpy(dtype=float)
    if len(xy) < 2:
        return 0.0
    return float(np.sqrt(np.sum(np.diff(xy, axis=0) ** 2, axis=1)).max())


def link_detections(detections: pd.DataFrame, cfg: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    if detections.empty:
        empty = detections.copy()
        empty["trackID"] = pd.Series(dtype=int)
        return empty, pd.DataFrame()
    from laptrack import LapTrack

    link = cfg["laptrack"]
    max_distance = float(link["max_link_distance_px"])
    tracker = LapTrack(
        metric=link.get("metric", "sqeuclidean"),
        cutoff=max_distance**2,
        gap_closing_cutoff=False if not link.get("gap_closing", False) else max_distance**2,
        splitting_cutoff=False,
        merging_cutoff=False,
    )
    tracked, _, _ = tracker.predict_dataframe(
        detections.reset_index(drop=True),
        coordinate_cols=["x", "y"],
        frame_col="frame",
        only_coordinate_cols=False,
    )
    tracked = tracked.rename(columns={"track_id": "trackID"})
    tracked["trackID"] = tracked["trackID"].astype(int)
    minimum = int(link["min_track_length"])
    lengths = tracked.groupby("trackID").size()
    keep_ids = lengths[lengths >= minimum].index
    tracked = tracked[tracked["trackID"].isin(keep_ids)].copy()
    rows = []
    enriched = []
    for track_id, group in tracked.groupby("trackID", sort=True):
        assignment = assign_track_labels(group)
        group = group.copy()
        for key, value in assignment.items():
            group[key] = value
        uid = f"{cfg['dataset']}/{group['fov'].iloc[0]}/{group['channel'].iloc[0]}/{int(track_id)}"
        group["trajectory_uid"] = uid
        enriched.append(group)
        rows.append(
            {
                "dataset": cfg["dataset"],
                "fov": group["fov"].iloc[0],
                "channel": group["channel"].iloc[0],
                "trackID": int(track_id),
                "trajectory_uid": uid,
                "n_detections": len(group),
                "first_frame": int(group["frame"].min()),
                "last_frame": int(group["frame"].max()),
                "max_step_px": _max_step(group),
                **assignment,
            }
        )
    canonical = pd.concat(enriched, ignore_index=True) if enriched else tracked
    canonical["t"] = canonical["frame"].astype(int)
    return canonical.sort_values(["trackID", "frame"]), pd.DataFrame.from_records(rows)


def legacy_table(canonical: pd.DataFrame) -> pd.DataFrame:
    if canonical.empty:
        return pd.DataFrame(columns=LEGACY_COLUMNS)
    mapping = canonical.copy()
    mapping["QUALITY"] = mapping["spotiflow_probability"]
    return mapping.reindex(columns=LEGACY_COLUMNS).astype(float)


def tracking_paths(output_dir: str | Path, fov: str, channel: str) -> dict[str, Path]:
    root = Path(output_dir) / "tracks"
    return {
        "canonical": root / f"{fov}__{channel}__trajectories.csv",
        "legacy": root / f"{fov}__{channel}__tracks_legacy.csv",
        "summary": root / f"{fov}__{channel}__track_summary.csv",
        "qc": Path(output_dir) / "qc" / f"{fov}__{channel}__tracks_qc.png",
    }


def _all_steps(tracks: pd.DataFrame) -> np.ndarray:
    steps: list[float] = []
    for _, group in tracks.groupby("trackID"):
        xy = group.sort_values("frame")[["x", "y"]].to_numpy(float)
        if len(xy) > 1:
            steps.extend(np.sqrt(np.sum(np.diff(xy, axis=0) ** 2, axis=1)))
    return np.asarray(steps)


def save_tracking_qc(
    mip: np.ndarray,
    mask: np.ndarray,
    tracks: pd.DataFrame,
    summary: pd.DataFrame,
    path: Path,
    title: str,
    max_tracks: int,
) -> None:
    cache = path.parents[1] / ".cache" / "matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes[0, 0].imshow(mip, cmap="gray", vmin=np.percentile(mip, 1), vmax=np.percentile(mip, 99.8))
    axes[0, 0].contour(find_boundaries(mask), levels=[0.5], colors="cyan", linewidths=0.35)
    for _, group in list(tracks.groupby("trackID"))[:max_tracks]:
        ordered = group.sort_values("frame")
        axes[0, 0].plot(ordered["x"], ordered["y"], linewidth=0.65)
    axes[0, 0].set_title(title)
    axes[0, 0].axis("off")
    if not summary.empty:
        axes[0, 1].hist(summary["n_detections"], bins=30)
        assigned = int((summary["cell_id"] > 0).sum())
        total = len(summary)
        axes[1, 1].bar(
            ["assigned", "cell_id=0"],
            [assigned / total, (total - assigned) / total],
        )
    axes[0, 1].set_title("Track lengths")
    axes[0, 1].set_xlabel("Detections")
    steps = _all_steps(tracks)
    if steps.size:
        axes[1, 0].hist(steps, bins=30)
    axes[1, 0].set_title("Linked step distances")
    axes[1, 0].set_xlabel("Pixels")
    axes[1, 1].set_title("Cell assignment fraction")
    axes[1, 1].set_ylim(0, 1)
    axes[1, 1].set_ylabel("Fraction of tracks")
    counts = tracks.groupby("frame").size() if not tracks.empty else pd.Series(dtype=int)
    if not counts.empty:
        axes[0, 2].plot(counts.index, counts.values, linewidth=1)
    axes[0, 2].set_title("Linked detections per frame")
    axes[0, 2].set_xlabel("Frame")
    axes[0, 2].set_ylabel("Detections")
    if not summary.empty:
        axes[1, 2].hist(summary["assignment_confidence"], bins=np.linspace(0, 1, 21))
    axes[1, 2].set_title("Assignment confidence")
    axes[1, 2].set_xlabel("Modal-label fraction")
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def track_file(
    row: dict[str, Any],
    cfg: dict[str, Any],
    *,
    resume: bool = False,
    force: bool = False,
) -> list[dict[str, str]]:
    mask = np.asarray(imread(segmentation_paths(cfg["output_dir"], row["fov"])["mask"]))
    results = []
    for spec in channel_specs(cfg):
        if spec["role"] != "spt":
            continue
        paths = tracking_paths(cfg["output_dir"], row["fov"], spec["name"])
        if output_action(list(paths.values()), resume=resume, force=force) == "skip":
            results.append({key: str(value) for key, value in paths.items()})
            continue
        for path in paths.values():
            path.parent.mkdir(parents=True, exist_ok=True)
        detections = pd.read_csv(detection_path(cfg["output_dir"], row["fov"], spec["name"]))
        canonical, summary = link_detections(detections, cfg)
        atomic_csv(canonical, paths["canonical"])
        atomic_csv(legacy_table(canonical), paths["legacy"])
        atomic_csv(summary, paths["summary"])
        mip = np.asarray(imread(artifact_paths(cfg["output_dir"], row["fov"], spec, cfg)["mip"]))
        save_tracking_qc(
            mip,
            mask,
            canonical,
            summary,
            paths["qc"],
            f"{row['fov']} — {spec['name']} — {len(summary)} tracks",
            int(cfg["qc"]["max_tracks_draw"]),
        )
        results.append({key: str(value) for key, value in paths.items()})
    return results
