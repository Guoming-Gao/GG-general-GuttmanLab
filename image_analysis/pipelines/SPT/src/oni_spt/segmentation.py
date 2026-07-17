from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from skimage import measure
from skimage.segmentation import find_boundaries
from tifffile import imread, imwrite

from .preprocess import artifact_paths, channel_specs
from .utils import atomic_csv, output_action


def normalize01(image: np.ndarray, low: float = 1, high: float = 99.8) -> np.ndarray:
    arr = image.astype(np.float32, copy=False)
    lo, hi = np.percentile(arr[np.isfinite(arr)], (low, high))
    if hi <= lo:
        return np.zeros_like(arr)
    return np.clip((arr - lo) / (hi - lo), 0, 1)


def choose_segmentation_image(row: dict[str, Any], cfg: dict[str, Any]) -> tuple[np.ndarray, str]:
    specs = channel_specs(cfg)
    requested = cfg["segmentation"].get("source", "auto")
    by_name = {spec["name"]: spec for spec in specs}
    if requested not in {"auto", "combined_spt_mip"}:
        if requested not in by_name:
            raise ValueError(f"Unknown segmentation source channel: {requested}")
        selected = [by_name[requested]]
    else:
        markers = [spec for spec in specs if spec["role"] == "marker"]
        spt = [spec for spec in specs if spec["role"] == "spt"]
        selected = markers or spt
        if requested == "auto" and markers:
            selected = markers[:1]
        elif requested == "auto" and len(spt) == 1:
            selected = spt
    images = []
    for spec in selected:
        paths = artifact_paths(cfg["output_dir"], row["fov"], spec, cfg)
        images.append(np.asarray(imread(paths["mip"])))
    if len(images) == 1:
        return images[0], selected[0]["name"]
    combined = np.maximum.reduce([normalize01(image) for image in images])
    return combined.astype(np.float32), "combined_spt_mip"


def resolve_cellpose_gpu(device: str) -> bool:
    if device == "cpu":
        return False
    from cellpose import core

    return bool(core.use_gpu())


def segment_image(image: np.ndarray, cfg: dict[str, Any], model: Any | None = None) -> np.ndarray:
    seg = cfg["segmentation"]
    if model is None:
        from cellpose import models

        model = models.CellposeModel(
            gpu=resolve_cellpose_gpu(str(seg.get("device", "auto"))),
            pretrained_model=seg.get("model", "cpsam"),
        )
    masks, _, _ = model.eval(
        image,
        diameter=float(seg["diameter"]),
        flow_threshold=float(seg["flow_threshold"]),
        cellprob_threshold=float(seg["cellprob_threshold"]),
        min_size=int(seg.get("min_size", 15)),
    )
    return np.asarray(masks, dtype=np.uint16)


def region_table(masks: np.ndarray) -> pd.DataFrame:
    rows = []
    for prop in measure.regionprops(masks):
        rows.append(
            {
                "cell_id": int(prop.label),
                "centroid_y": float(prop.centroid[0]),
                "centroid_x": float(prop.centroid[1]),
                "bbox_min_y": int(prop.bbox[0]),
                "bbox_min_x": int(prop.bbox[1]),
                "bbox_max_y": int(prop.bbox[2]),
                "bbox_max_x": int(prop.bbox[3]),
                "area_px": int(prop.area),
            }
        )
    return pd.DataFrame.from_records(rows)


def save_segmentation_overlay(image: np.ndarray, masks: np.ndarray, path: Path, title: str) -> None:
    cache = path.parents[1] / ".cache" / "matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    display = normalize01(image)
    boundary = find_boundaries(masks, mode="outer")
    rgb = np.repeat(display[..., None], 3, axis=-1)
    rgb[boundary] = (1, 0.1, 0.1)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(rgb)
    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def segmentation_paths(output_dir: str | Path, fov: str) -> dict[str, Path]:
    root = Path(output_dir) / "segmentation"
    return {
        "source": root / f"{fov}__segmentation_source.tif",
        "mask": root / f"{fov}__cellposeSAM_masks.tif",
        "regions": root / f"{fov}__cells.csv",
        "overlay": root / f"{fov}__segmentation_overlay.png",
    }


def segment_file(
    row: dict[str, Any],
    cfg: dict[str, Any],
    *,
    resume: bool = False,
    force: bool = False,
    model: Any | None = None,
) -> dict[str, str]:
    paths = segmentation_paths(cfg["output_dir"], row["fov"])
    if output_action(list(paths.values()), resume=resume, force=force) == "skip":
        return {key: str(value) for key, value in paths.items()}
    paths["mask"].parent.mkdir(parents=True, exist_ok=True)
    image, source = choose_segmentation_image(row, cfg)
    masks = segment_image(image, cfg, model=model)
    if masks.shape != image.shape:
        raise ValueError(f"Cellpose mask shape {masks.shape} != source shape {image.shape}")
    imwrite(paths["source"], image, metadata=None)
    imwrite(paths["mask"], masks, metadata=None)
    table = region_table(masks)
    table.insert(0, "segmentation_source", source)
    table.insert(0, "fov", row["fov"])
    table.insert(0, "dataset", cfg["dataset"])
    atomic_csv(table, paths["regions"])
    save_segmentation_overlay(
        image, masks, paths["overlay"], f"{row['fov']} — {source} — {int(masks.max())} masks"
    )
    return {key: str(value) for key, value in paths.items()}
