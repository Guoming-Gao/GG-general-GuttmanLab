from __future__ import annotations

import os
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
from skimage.segmentation import find_boundaries
from tifffile import TiffFile, imread

from .preprocess import artifact_paths, channel_specs
from .segmentation import normalize01, segmentation_paths
from .utils import atomic_csv, output_action


DETECTION_COLUMNS = [
    "dataset",
    "fov",
    "channel",
    "spotID",
    "frame",
    "t",
    "x",
    "y",
    "spotiflow_probability",
    "spotiflow_intensity",
    "raw_center_intensity",
    "raw_mean_intensity",
    "R",
    "meanIntensity",
    "medianIntensity",
    "minIntensity",
    "maxIntensity",
    "totalIntensity",
    "stdIntensity",
    "contrast",
    "SNR",
    "detection_cell_id",
]


def resolve_torch_device(name: str) -> str:
    import torch

    name = str(name).lower()
    if name in {"auto", "gpu"}:
        if torch.backends.mps.is_available():
            return "mps"
        if torch.cuda.is_available():
            return "cuda"
        return "cpu"
    if name == "mps" and not torch.backends.mps.is_available():
        raise RuntimeError("MPS requested but unavailable")
    if name == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    if name not in {"cpu", "mps", "cuda"}:
        raise ValueError("Spotiflow device must be auto, cpu, mps, or cuda")
    return name


def load_spotiflow(cfg: dict[str, Any]) -> SimpleNamespace:
    from spotiflow.model import Spotiflow

    spot = cfg["spotiflow"]
    cache = Path(spot["model_cache_dir"])
    cache.mkdir(parents=True, exist_ok=True)
    device = resolve_torch_device(spot.get("device", "auto"))
    model = Spotiflow.from_pretrained(
        spot.get("pretrained_model", "general"),
        cache_dir=cache,
        map_location=device,
    )
    return SimpleNamespace(model=model, device=device)


def disk_statistics(image: np.ndarray, y: float, x: float, radius: float) -> dict[str, float]:
    y0 = max(0, int(np.floor(y - radius)))
    y1 = min(image.shape[0], int(np.ceil(y + radius)) + 1)
    x0 = max(0, int(np.floor(x - radius)))
    x1 = min(image.shape[1], int(np.ceil(x + radius)) + 1)
    yy, xx = np.ogrid[y0:y1, x0:x1]
    keep = (yy - y) ** 2 + (xx - x) ** 2 <= radius**2
    values = image[y0:y1, x0:x1][keep].astype(float)
    if values.size == 0:
        values = np.array([np.nan])
    mean = float(np.nanmean(values))
    std = float(np.nanstd(values))
    minimum = float(np.nanmin(values))
    maximum = float(np.nanmax(values))
    return {
        "meanIntensity": mean,
        "medianIntensity": float(np.nanmedian(values)),
        "minIntensity": minimum,
        "maxIntensity": maximum,
        "totalIntensity": float(np.nansum(values)),
        "stdIntensity": std,
        "contrast": (maximum - minimum) / (maximum + minimum) if maximum + minimum else 0.0,
        "SNR": mean / std if std > 0 else 0.0,
    }


def _detail_value(details: Any, name: str, index: int) -> float:
    value = getattr(details, name, None)
    if value is None or len(value) <= index:
        return float("nan")
    # Spotiflow 0.6 returns probabilities as (N,) but intensities as (N, C).
    # ONI SPT inference is always performed on one grayscale channel, so accept
    # either representation and take that channel's sole value.
    item = np.asarray(value[index], dtype=float)
    if item.size == 0:
        return float("nan")
    return float(item.reshape(-1)[0])


def detect_stack(
    raw_path: Path,
    filtered_path: Path,
    mask: np.ndarray,
    cfg: dict[str, Any],
    context: SimpleNamespace,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    spot = cfg["spotiflow"]
    radius = float(cfg["legacy"]["spot_radius_px"])
    spot_id = 0
    with TiffFile(raw_path) as raw_tif, TiffFile(filtered_path) as filtered_tif:
        if len(raw_tif.pages) != len(filtered_tif.pages):
            raise ValueError("Raw and filtered stacks have different frame counts")
        for frame_index, (raw_page, filtered_page) in enumerate(
            zip(raw_tif.pages, filtered_tif.pages)
        ):
            raw = raw_page.asarray()
            filtered = filtered_page.asarray()
            points, details = context.model.predict(
                filtered,
                prob_thresh=spot.get("probability_threshold"),
                min_distance=int(spot.get("min_distance", 1)),
                exclude_border=int(spot.get("exclude_border", 1)),
                normalizer=spot.get("normalizer", "auto"),
                subpix=True,
                device=context.device,
                verbose=False,
            )
            for index, (y, x) in enumerate(np.asarray(points, dtype=float)):
                stats = disk_statistics(filtered, float(y), float(x), radius)
                raw_stats = disk_statistics(raw, float(y), float(x), radius)
                yi = int(np.clip(np.rint(y), 0, mask.shape[0] - 1))
                xi = int(np.clip(np.rint(x), 0, mask.shape[1] - 1))
                row = {
                    "spotID": spot_id,
                    "frame": frame_index,
                    "t": frame_index,
                    "x": float(x),
                    "y": float(y),
                    "spotiflow_probability": _detail_value(details, "prob", index),
                    "spotiflow_intensity": _detail_value(details, "intens", index),
                    "raw_center_intensity": float(raw[yi, xi]),
                    "raw_mean_intensity": raw_stats["meanIntensity"],
                    "R": radius,
                    "detection_cell_id": int(mask[yi, xi]),
                    **stats,
                }
                rows.append(row)
                spot_id += 1
    return pd.DataFrame.from_records(rows)


def detection_path(output_dir: str | Path, fov: str, channel: str) -> Path:
    return Path(output_dir) / "detections" / f"{fov}__{channel}__detections.csv"


def save_detection_overlay(
    mip: np.ndarray,
    mask: np.ndarray,
    detections: pd.DataFrame,
    path: Path,
    title: str,
    max_points: int,
) -> None:
    cache = path.parents[1] / ".cache" / "matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(mip, cmap="gray", vmin=np.percentile(mip, 1), vmax=np.percentile(mip, 99.8))
    boundary = find_boundaries(mask, mode="outer")
    ax.contour(boundary, levels=[0.5], colors="cyan", linewidths=0.4)
    points = detections.head(max_points)
    if not points.empty:
        ax.scatter(points["x"], points["y"], s=8, facecolors="none", edgecolors="yellow", linewidths=0.4)
    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def detect_file(
    row: dict[str, Any],
    cfg: dict[str, Any],
    context: SimpleNamespace,
    *,
    resume: bool = False,
    force: bool = False,
) -> list[str]:
    mask = np.asarray(imread(segmentation_paths(cfg["output_dir"], row["fov"])["mask"]))
    outputs = []
    for spec in channel_specs(cfg):
        if spec["role"] != "spt":
            continue
        paths = artifact_paths(cfg["output_dir"], row["fov"], spec, cfg)
        output = detection_path(cfg["output_dir"], row["fov"], spec["name"])
        overlay = Path(cfg["output_dir"]) / "qc" / f"{row['fov']}__{spec['name']}__detections.png"
        if output_action([output, overlay], resume=resume, force=force) == "skip":
            outputs.append(str(output))
            continue
        output.parent.mkdir(parents=True, exist_ok=True)
        overlay.parent.mkdir(parents=True, exist_ok=True)
        detections = detect_stack(paths["raw"], paths["filtered"], mask, cfg, context)
        for column, value in [("channel", spec["name"]), ("fov", row["fov"]), ("dataset", cfg["dataset"])]:
            detections.insert(0, column, value)
        detections = detections.reindex(columns=DETECTION_COLUMNS)
        atomic_csv(detections, output)
        save_detection_overlay(
            np.asarray(imread(paths["mip"])),
            mask,
            detections,
            overlay,
            f"{row['fov']} — {spec['name']} — {len(detections)} detections",
            int(cfg["qc"]["max_detections_draw"]),
        )
        outputs.append(str(output))
    return outputs
