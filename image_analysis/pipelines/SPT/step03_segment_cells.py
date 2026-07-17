"""Step 3: segment cells/nuclei from a marker MIP or SPT MIP with CellposeSAM."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from skimage import measure
from skimage.segmentation import find_boundaries
from tifffile import imread, imwrite

from spt_shared import atomic_csv, output_action, output_dirs, stage_timer
from step01_inspect_inputs import accepted_rows, inspect_inputs
from step02_preprocess_videos import artifact_paths, channel_specs


def normalize01(image: np.ndarray, low: float = 1, high: float = 99.8) -> np.ndarray:
    array = image.astype(np.float32, copy=False)
    finite = array[np.isfinite(array)]
    if not finite.size:
        return np.zeros_like(array)
    lower, upper = np.percentile(finite, (low, high))
    return np.zeros_like(array) if upper <= lower else np.clip((array - lower) / (upper - lower), 0, 1)


def choose_segmentation_image(row: dict[str, Any], cfg: dict[str, Any]) -> tuple[np.ndarray, str]:
    specs = channel_specs(cfg)
    requested = cfg["segmentation"].get("source", "auto")
    by_name = {spec["name"]: spec for spec in specs}
    if requested not in {"auto", "combined_spt_mip"}:
        if requested not in by_name:
            raise ValueError(f"Unknown segmentation source {requested}")
        selected = [by_name[requested]]
    else:
        markers = [spec for spec in specs if spec["role"] == "marker"]
        spt = [spec for spec in specs if spec["role"] == "spt"]
        selected = markers[:1] if requested == "auto" and markers else spt
    images = [np.asarray(imread(artifact_paths(cfg["output_dir"], row["fov"], spec, cfg)["mip"])) for spec in selected]
    if len(images) == 1:
        return images[0], selected[0]["name"]
    return np.maximum.reduce([normalize01(image) for image in images]).astype(np.float32), "combined_spt_mip"


def segmentation_paths(output_dir: str | Path, fov: str) -> dict[str, Path]:
    root = output_dirs(output_dir)["segmentation"]
    return {
        "source": root / f"{fov}__segmentation_source.tif",
        "mask": root / f"{fov}__cellposeSAM_masks.tif",
        "regions": root / f"{fov}__cells.csv",
        "overlay": root / f"{fov}__segmentation_overlay.png",
    }


def segment_image(image: np.ndarray, cfg: dict[str, Any], model: Any | None = None) -> np.ndarray:
    settings = cfg["segmentation"]
    if model is None:
        from cellpose import core, models
        gpu = False if settings.get("device") == "cpu" else bool(core.use_gpu())
        model = models.CellposeModel(gpu=gpu, pretrained_model=settings.get("model", "cpsam"))
    masks, _, _ = model.eval(
        image, diameter=float(settings["diameter"]), flow_threshold=float(settings["flow_threshold"]),
        cellprob_threshold=float(settings["cellprob_threshold"]), min_size=int(settings.get("min_size", 15)),
    )
    return np.asarray(masks, dtype=np.uint16)


def region_table(mask: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame.from_records([{
        "cell_id": int(prop.label), "centroid_y": float(prop.centroid[0]), "centroid_x": float(prop.centroid[1]),
        "bbox_min_y": int(prop.bbox[0]), "bbox_min_x": int(prop.bbox[1]),
        "bbox_max_y": int(prop.bbox[2]), "bbox_max_x": int(prop.bbox[3]), "area_px": int(prop.area),
    } for prop in measure.regionprops(mask)])


def save_overlay(image: np.ndarray, mask: np.ndarray, path: Path, title: str) -> None:
    cache = output_dirs(path.parents[1])["root"] / ".cache" / "matplotlib"
    cache.mkdir(parents=True, exist_ok=True); os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rgb = np.repeat(normalize01(image)[..., None], 3, axis=-1)
    rgb[find_boundaries(mask, mode="outer")] = (1, 0.1, 0.1)
    fig, ax = plt.subplots(figsize=(6, 8)); ax.imshow(rgb); ax.set_title(title); ax.axis("off")
    fig.tight_layout(); fig.savefig(path, dpi=300, bbox_inches="tight"); plt.close(fig)


def segment_file(row: dict[str, Any], cfg: dict[str, Any], *, resume: bool = False,
                 force: bool = False, model: Any | None = None) -> dict[str, str]:
    paths = segmentation_paths(cfg["output_dir"], row["fov"])
    if output_action(list(paths.values()), resume=resume, force=force) == "skip":
        return {key: str(value) for key, value in paths.items()}
    paths["mask"].parent.mkdir(parents=True, exist_ok=True)
    image, source = choose_segmentation_image(row, cfg)
    mask = segment_image(image, cfg, model)
    if mask.shape != image.shape:
        raise ValueError("Cellpose mask shape differs from source")
    imwrite(paths["source"], image, metadata=None); imwrite(paths["mask"], mask, metadata=None)
    table = region_table(mask)
    for column, value in reversed([("dataset", cfg["dataset"]), ("fov", row["fov"]), ("segmentation_source", source)]):
        table.insert(0, column, value)
    atomic_csv(table, paths["regions"])
    save_overlay(image, mask, paths["overlay"], f"{row['fov']} — {source} — {int(mask.max())} masks")
    return {key: str(value) for key, value in paths.items()}


def segment_stage(cfg: dict[str, Any], manifest: pd.DataFrame, *, max_files: int | None = None,
                  resume: bool = False, force: bool = False, model: Any | None = None) -> None:
    rows = accepted_rows(manifest, max_files)
    with stage_timer(cfg, "03_segment", {"files": len(rows)}):
        for index, row in enumerate(rows, 1):
            print(f"[segment {index}/{len(rows)}] {row['filename']}", flush=True)
            segment_file(row, cfg, resume=resume, force=force, model=model)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", required=True)
    parser.add_argument("--max-files", type=int); policy = parser.add_mutually_exclusive_group()
    policy.add_argument("--resume", action="store_true"); policy.add_argument("--force", action="store_true")
    args = parser.parse_args(argv); cfg, manifest = inspect_inputs(args.config)
    segment_stage(cfg, manifest, max_files=args.max_files, resume=args.resume, force=args.force)


if __name__ == "__main__":
    main()

