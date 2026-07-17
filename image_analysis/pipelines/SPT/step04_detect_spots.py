"""Step 4: detect frame-resolved SPT spots with Spotiflow and render detection QC MP4s."""

from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
from tifffile import TiffFile, imread

from spt_shared import atomic_csv, output_action, output_dirs, stage_timer
from spt_video import render_overlay_video
from step01_inspect_inputs import accepted_rows, inspect_inputs
from step02_preprocess_videos import artifact_paths, channel_specs
from step03_segment_cells import segmentation_paths


DETECTION_COLUMNS = [
    "dataset", "condition", "fov", "channel", "spotID", "frame", "t", "x", "y",
    "spotiflow_probability", "spotiflow_intensity", "raw_center_intensity", "raw_mean_intensity",
    "R", "meanIntensity", "medianIntensity", "minIntensity", "maxIntensity", "totalIntensity",
    "stdIntensity", "contrast", "SNR", "signed_dog_mean", "signed_dog_median",
    "signed_dog_min", "signed_dog_max", "signed_dog_std", "detection_cell_id",
]


def resolve_torch_device(name: str) -> str:
    import torch
    name = str(name).lower()
    if name in {"auto", "gpu"}:
        if torch.backends.mps.is_available(): return "mps"
        if torch.cuda.is_available(): return "cuda"
        return "cpu"
    if name == "mps" and not torch.backends.mps.is_available(): raise RuntimeError("MPS unavailable")
    if name == "cuda" and not torch.cuda.is_available(): raise RuntimeError("CUDA unavailable")
    if name not in {"cpu", "mps", "cuda"}: raise ValueError("device must be auto/cpu/mps/cuda")
    return name


def load_spotiflow(cfg: dict[str, Any]) -> SimpleNamespace:
    from spotiflow.model import Spotiflow
    settings = cfg["spotiflow"]; cache = Path(settings["model_cache_dir"]); cache.mkdir(parents=True, exist_ok=True)
    device = resolve_torch_device(settings.get("device", "auto"))
    model = Spotiflow.from_pretrained(settings.get("pretrained_model", "general"), cache_dir=cache, map_location=device)
    stored = getattr(model, "_prob_thresh", [0.5])
    resolved = float(stored[0] if isinstance(stored, (list, tuple, np.ndarray)) else stored)
    if settings.get("probability_threshold") is not None:
        resolved = float(settings["probability_threshold"])
    return SimpleNamespace(model=model, device=device, resolved_threshold=resolved)


def disk_values(image: np.ndarray, y: float, x: float, radius: float) -> np.ndarray:
    y0, y1 = max(0, int(np.floor(y-radius))), min(image.shape[0], int(np.ceil(y+radius))+1)
    x0, x1 = max(0, int(np.floor(x-radius))), min(image.shape[1], int(np.ceil(x+radius))+1)
    yy, xx = np.ogrid[y0:y1, x0:x1]
    return image[y0:y1, x0:x1][(yy-y)**2 + (xx-x)**2 <= radius**2].astype(float)


def legacy_statistics(image: np.ndarray, y: float, x: float, radius: float) -> dict[str, float]:
    values = disk_values(image, y, x, radius)
    if not values.size: values = np.array([np.nan])
    mean, std = float(np.nanmean(values)), float(np.nanstd(values))
    minimum, maximum = float(np.nanmin(values)), float(np.nanmax(values))
    return {
        "meanIntensity": mean, "medianIntensity": float(np.nanmedian(values)), "minIntensity": minimum,
        "maxIntensity": maximum, "totalIntensity": float(np.nansum(values)), "stdIntensity": std,
        "contrast": (maximum-minimum)/(maximum+minimum) if maximum+minimum else 0.0,
        "SNR": mean/std if std > 0 else 0.0,
    }


def _detail_value(details: Any, name: str, index: int) -> float:
    values = getattr(details, name, None)
    if values is None or len(values) <= index: return float("nan")
    item = np.asarray(values[index], dtype=float)
    return float(item.reshape(-1)[0]) if item.size else float("nan")


def detection_paths(output_dir: str | Path, fov: str, channel: str) -> dict[str, Path]:
    dirs = output_dirs(output_dir)
    return {
        "table": dirs["detections"] / f"{fov}__{channel}__detections.csv",
        "video": dirs["qc"] / "detection_videos" / f"{fov}__{channel}__spot_detection_QC.mp4",
    }


def detect_stack(raw_path: Path, filtered_path: Path, mask: np.ndarray, cfg: dict[str, Any],
                 context: SimpleNamespace) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []; spot_id = 0
    radius = float(cfg["legacy"]["spot_radius_px"]); settings = cfg["spotiflow"]
    with TiffFile(raw_path) as raw_tif, TiffFile(filtered_path) as dog_tif:
        if len(raw_tif.pages) != len(dog_tif.pages): raise ValueError("Raw/DoG frame mismatch")
        for frame_index, (raw_page, dog_page) in enumerate(zip(raw_tif.pages, dog_tif.pages)):
            raw = raw_page.asarray(); dog = dog_page.asarray().astype(np.float32)
            points, details = context.model.predict(
                dog, prob_thresh=settings.get("probability_threshold"), min_distance=int(settings.get("min_distance", 1)),
                exclude_border=int(settings.get("exclude_border", 1)), normalizer=settings.get("normalizer", "auto"),
                subpix=True, device=context.device, verbose=False,
            )
            positive_dog = np.maximum(dog, 0)
            for index, (y, x) in enumerate(np.asarray(points, dtype=float)):
                legacy = legacy_statistics(positive_dog, y, x, radius)
                raw_values = disk_values(raw, y, x, radius); signed = disk_values(dog, y, x, radius)
                yi, xi = int(np.clip(np.rint(y), 0, mask.shape[0]-1)), int(np.clip(np.rint(x), 0, mask.shape[1]-1))
                rows.append({
                    "spotID": spot_id, "frame": frame_index, "t": frame_index, "x": float(x), "y": float(y),
                    "spotiflow_probability": _detail_value(details, "prob", index),
                    "spotiflow_intensity": _detail_value(details, "intens", index),
                    "raw_center_intensity": float(raw[yi, xi]),
                    "raw_mean_intensity": float(np.mean(raw_values)) if raw_values.size else np.nan,
                    "R": radius, **legacy,
                    "signed_dog_mean": float(np.mean(signed)), "signed_dog_median": float(np.median(signed)),
                    "signed_dog_min": float(np.min(signed)), "signed_dog_max": float(np.max(signed)),
                    "signed_dog_std": float(np.std(signed)), "detection_cell_id": int(mask[yi, xi]),
                }); spot_id += 1
    return pd.DataFrame.from_records(rows)


def detect_file(row: dict[str, Any], cfg: dict[str, Any], context: SimpleNamespace, *,
                resume: bool = False, force: bool = False) -> list[str]:
    mask = np.asarray(imread(segmentation_paths(cfg["output_dir"], row["fov"])["mask"]))
    outputs = []
    for spec in channel_specs(cfg):
        if spec["role"] != "spt": continue
        artifacts = artifact_paths(cfg["output_dir"], row["fov"], spec, cfg)
        paths = detection_paths(cfg["output_dir"], row["fov"], spec["name"])
        if output_action(list(paths.values()), resume=resume, force=force) == "skip":
            outputs.append(str(paths["table"])); continue
        detections = detect_stack(artifacts["raw"], artifacts["filtered"], mask, cfg, context)
        for column, value in reversed([
            ("dataset", cfg["dataset"]), ("condition", row["condition"]), ("fov", row["fov"]), ("channel", spec["name"]),
        ]): detections.insert(0, column, value)
        detections = detections.reindex(columns=DETECTION_COLUMNS); atomic_csv(detections, paths["table"])
        qc = cfg["qc"]
        render_overlay_video(
            artifacts["filtered"], paths["video"], detections, tracks=None, mask=mask,
            frame_interval_s=float(row["frame_interval_s"]), pixel_size_um=float(row["pixel_size_um"]),
            model_label=settings_label(cfg), threshold_label=f"{context.resolved_threshold:g}",
            maximum_frames=int(qc["max_video_frames"]), scale=int(qc["render_scale"]), fps=int(qc["video_fps"]),
            crf=int(qc["video_crf"]), scale_bar_um=float(qc["scale_bar_um"]),
        )
        outputs.append(str(paths["table"]))
    return outputs


def settings_label(cfg: dict[str, Any]) -> str:
    return f"Spotiflow {cfg['spotiflow'].get('pretrained_model', 'general')}"


def detect_stage(cfg: dict[str, Any], manifest: pd.DataFrame, *, max_files: int | None = None,
                 resume: bool = False, force: bool = False, context: Any | None = None) -> None:
    rows = accepted_rows(manifest, max_files); context = context or load_spotiflow(cfg)
    with stage_timer(cfg, "04_detect", {"files": len(rows), "threshold": context.resolved_threshold}):
        for index, row in enumerate(rows, 1):
            print(f"[detect {index}/{len(rows)}] {row['filename']}", flush=True)
            detect_file(row, cfg, context, resume=resume, force=force)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", required=True)
    parser.add_argument("--max-files", type=int); policy = parser.add_mutually_exclusive_group()
    policy.add_argument("--resume", action="store_true"); policy.add_argument("--force", action="store_true")
    args = parser.parse_args(argv); cfg, manifest = inspect_inputs(args.config)
    detect_stage(cfg, manifest, max_files=args.max_files, resume=args.resume, force=args.force)


if __name__ == "__main__":
    main()
