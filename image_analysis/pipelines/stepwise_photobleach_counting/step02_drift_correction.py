"""Step 2: estimate and apply subpixel cross-correlation drift correction."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import gaussian_filter, shift as ndi_shift
from skimage.registration import phase_cross_correlation

from pbsa_shared import atomic_csv, load_config, matplotlib_setup, output_action, output_dirs, require_free_space, stage_timer
from step01_inspect_inputs import load_manifest


def drift_correction_paths(output_dir: str | Path, fov: str) -> dict[str, Path]:
    root = output_dirs(output_dir)["drift_correction"]
    return {
        "stack": root / "corrected_tiffs" / f"{fov}__drift_corrected.tif",
        "table": root / "shift_tables" / f"{fov}__drift_correction.csv",
        "plot": root / "plots" / f"{fov}__drift_correction.png",
        "overlay": root / "overlays" / f"{fov}__drift_correction_overlay.png",
    }


def registration_image(image: np.ndarray, sigma: float) -> np.ndarray:
    array = image.astype(np.float32, copy=False)
    filtered = array - gaussian_filter(array, sigma)
    filtered -= np.mean(filtered)
    scale = float(np.std(filtered))
    return filtered / scale if scale > 0 else filtered


def read_block(tif: tifffile.TiffFile, start: int, stop: int) -> np.ndarray:
    images = [tif.pages[index].asarray().astype(np.float32) for index in range(start, stop)]
    return np.mean(images, axis=0, dtype=np.float32)


def estimate_block_shifts(path: str | Path, settings: dict[str, Any], max_frames: int | None = None) -> tuple[pd.DataFrame, list[np.ndarray]]:
    block_size = int(settings["block_frames"]); sigma = float(settings["highpass_sigma_px"])
    with tifffile.TiffFile(path) as tif:
        total = min(len(tif.pages), max_frames) if max_frames else len(tif.pages)
        starts = list(range(0, total, block_size)); blocks = [read_block(tif, start, min(total, start + block_size)) for start in starts]
    if not blocks: raise ValueError(f"No frames in {path}")
    reference = registration_image(blocks[0], sigma); records = []
    for index, (start, block) in enumerate(zip(starts, blocks)):
        moving = registration_image(block, sigma)
        shift, error, phase = phase_cross_correlation(
            reference, moving, upsample_factor=int(settings["upsample_factor"]), normalization=None,
        )
        valid = bool(np.all(np.isfinite(shift)) and np.max(np.abs(shift)) <= float(settings["max_abs_shift_px"]))
        records.append({"block": index, "frame_start": start, "frame_center": start + (min(total, start + block_size) - start - 1) / 2,
                        "estimated_y_shift_px": float(shift[0]), "estimated_x_shift_px": float(shift[1]),
                        "registration_error": float(error), "phase_difference": float(phase), "estimate_valid": valid})
    table = pd.DataFrame.from_records(records)
    valid = table.estimate_valid.to_numpy(bool)
    if valid.sum() == 0: raise RuntimeError(f"No valid drift estimates for {path}")
    x = table.frame_center.to_numpy(float); xv = x[valid]
    for column in ["estimated_y_shift_px", "estimated_x_shift_px"]:
        values = table[column].to_numpy(float)
        table[column.replace("estimated_", "interpolated_")] = np.interp(x, xv, values[valid])
    return table, blocks


def frame_shifts(table: pd.DataFrame, frame_count: int) -> np.ndarray:
    centers = table.frame_center.to_numpy(float); frames = np.arange(frame_count, dtype=float)
    return np.column_stack([
        np.interp(frames, centers, table.interpolated_y_shift_px, left=table.interpolated_y_shift_px.iloc[0], right=table.interpolated_y_shift_px.iloc[-1]),
        np.interp(frames, centers, table.interpolated_x_shift_px, left=table.interpolated_x_shift_px.iloc[0], right=table.interpolated_x_shift_px.iloc[-1]),
    ])


def write_corrected_stack(source: str | Path, destination: Path, shifts: np.ndarray, settings: dict[str, Any]) -> tuple[int, float]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.partial")
    margin = int(math.ceil(float(np.max(np.abs(shifts))))) + 1
    before_means, after_means = [], []
    with tifffile.TiffFile(source) as tif, tifffile.TiffWriter(temporary, bigtiff=True) as writer:
        first_description = tif.pages[0].description
        for index in range(len(shifts)):
            raw = tif.pages[index].asarray()
            corrected = ndi_shift(raw.astype(np.float32), shifts[index], order=int(settings["interpolation_order"]),
                                  mode="constant", cval=0.0, prefilter=False)
            if margin:
                corrected = corrected[margin:-margin, margin:-margin]
                original_valid = raw[margin:-margin, margin:-margin]
            else: original_valid = raw
            before_means.append(float(np.mean(original_valid))); after_means.append(float(np.mean(corrected)))
            corrected = np.clip(np.rint(corrected), 0, np.iinfo(raw.dtype).max).astype(raw.dtype)
            description = None
            if index == 0:
                try: meta = json.loads(first_description or "{}")
                except json.JSONDecodeError: meta = {"original_description": first_description}
                meta.update({"PBSA_DriftCorrected": True, "PBSA_DriftCropMargin_px": margin,
                             "PBSA_DriftMethod": "skimage.phase_cross_correlation"})
                description = json.dumps(meta)
            writer.write(corrected, compression=settings["compression"],
                         compressionargs={"level": int(settings["compression_level"])}, metadata=None, description=description)
    os.replace(temporary, destination)
    preservation_error = float(np.median(np.abs(np.asarray(after_means) - np.asarray(before_means))))
    return margin, preservation_error


def residual_shifts(blocks: list[np.ndarray], table: pd.DataFrame, settings: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    sigma = float(settings["highpass_sigma_px"]); reference = registration_image(blocks[0], sigma)
    ry, rx = [], []
    for block, row in zip(blocks, table.itertuples()):
        corrected = ndi_shift(block, [row.interpolated_y_shift_px, row.interpolated_x_shift_px], order=1, mode="nearest", prefilter=False)
        shift, _, _ = phase_cross_correlation(reference, registration_image(corrected, sigma),
                                               upsample_factor=int(settings["upsample_factor"]), normalization=None)
        ry.append(float(shift[0])); rx.append(float(shift[1]))
    return np.asarray(ry), np.asarray(rx)


def save_file_qc(source: Path, paths: dict[str, Path], table: pd.DataFrame, blocks: list[np.ndarray], settings: dict[str, Any], dataset: dict[str, Any]) -> None:
    matplotlib_setup(dataset["output_dir"]); import matplotlib.pyplot as plt
    paths["plot"].parent.mkdir(parents=True, exist_ok=True); paths["overlay"].parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axes[0].plot(table.frame_center, table.interpolated_x_shift_px, label="X"); axes[0].plot(table.frame_center, table.interpolated_y_shift_px, label="Y")
    axes[0].set_ylabel("Applied shift (px)"); axes[0].legend(); axes[0].set_title(source.name)
    axes[1].plot(table.frame_center, table.residual_x_shift_px, label="X residual"); axes[1].plot(table.frame_center, table.residual_y_shift_px, label="Y residual")
    axes[1].axhline(0, color="black", lw=.5); axes[1].set_ylabel("Residual (px)"); axes[1].set_xlabel("Frame"); axes[1].legend()
    fig.tight_layout(); fig.savefig(paths["plot"], dpi=200); plt.close(fig)
    indices = sorted(set([0, len(blocks)//2, len(blocks)-1])); colors = [(1,0,0),(0,1,0),(0,0,1)]
    before = np.zeros((*blocks[0].shape, 3), np.float32); after = np.zeros_like(before)
    def norm(image: np.ndarray) -> np.ndarray:
        lo, hi = np.percentile(image, (1, 99.8)); return np.clip((image-lo)/(hi-lo), 0, 1) if hi > lo else np.zeros_like(image)
    for idx, color in zip(indices, colors):
        raw = norm(blocks[idx]); shifted = norm(ndi_shift(blocks[idx], [table.interpolated_y_shift_px.iloc[idx], table.interpolated_x_shift_px.iloc[idx]], order=1, mode="nearest"))
        before += raw[..., None] * np.asarray(color); after += shifted[..., None] * np.asarray(color)
    fig, axes = plt.subplots(1, 2, figsize=(10, 7)); axes[0].imshow(np.clip(before,0,1)); axes[0].set_title("Before correction"); axes[1].imshow(np.clip(after,0,1)); axes[1].set_title("After correction")
    for ax in axes: ax.axis("off")
    fig.tight_layout(); fig.savefig(paths["overlay"], dpi=200); plt.close(fig)


def correct_file(row: dict[str, Any], cfg: dict[str, Any], dataset: dict[str, Any], *, resume: bool, force: bool, max_frames: int | None) -> dict[str, Any]:
    paths = drift_correction_paths(dataset["output_dir"], row["fov"])
    if output_action(list(paths.values()), resume=resume, force=force) == "skip":
        table = pd.read_csv(paths["table"])
        residual = np.hypot(table.residual_y_shift_px, table.residual_x_shift_px)
        median = float(np.median(residual)); p95 = float(np.percentile(residual, 95))
        status = "effective" if median <= float(cfg["qc"]["drift_median_residual_target_px"]) and p95 <= float(cfg["qc"]["drift_p95_residual_target_px"]) else "warning"
        return {"filename": row["filename"], "fov": row["fov"], "status": status,
                "frames": int(row["frames"]), "median_residual_shift_px": median,
                "p95_residual_shift_px": p95, "invalid_block_estimates": int((~table.estimate_valid.astype(bool)).sum()),
                "crop_margin_px": int(table.crop_margin_px.iloc[0]), "corrected_tiff": str(paths["stack"]),
                "qc_plot": str(paths["plot"]), "qc_overlay": str(paths["overlay"])}
    settings = cfg["drift_correction"]; table, blocks = estimate_block_shifts(row["path"], settings, max_frames=max_frames)
    total_frames = min(int(row["frames"]), max_frames) if max_frames else int(row["frames"])
    shifts = frame_shifts(table, total_frames); ry, rx = residual_shifts(blocks, table, settings)
    table["residual_y_shift_px"] = ry; table["residual_x_shift_px"] = rx
    margin, intensity_error = write_corrected_stack(row["path"], paths["stack"], shifts, settings)
    table.insert(0, "filename", row["filename"]); table.insert(1, "fov", row["fov"]); table["crop_margin_px"] = margin
    table["registration_block_frames"] = int(settings["block_frames"])
    table["registration_upsample_factor"] = int(settings["upsample_factor"])
    table["registration_highpass_sigma_px"] = float(settings["highpass_sigma_px"])
    atomic_csv(table, paths["table"]); save_file_qc(Path(row["path"]), paths, table, blocks, settings, dataset)
    residual = np.hypot(ry, rx); median = float(np.median(residual)); p95 = float(np.percentile(residual,95))
    status = "effective" if median <= float(cfg["qc"]["drift_median_residual_target_px"]) and p95 <= float(cfg["qc"]["drift_p95_residual_target_px"]) else "warning"
    return {"filename": row["filename"], "fov": row["fov"], "status": status, "frames": total_frames,
            "maximum_applied_shift_px": float(np.max(np.hypot(shifts[:,0], shifts[:,1]))), "median_residual_shift_px": median,
            "p95_residual_shift_px": p95, "invalid_block_estimates": int((~table.estimate_valid).sum()),
            "crop_margin_px": margin, "median_absolute_intensity_change_counts": intensity_error,
            "corrected_tiff": str(paths["stack"]), "qc_plot": str(paths["plot"]), "qc_overlay": str(paths["overlay"])}


def save_stage_report(summary: pd.DataFrame, dataset: dict[str, Any]) -> None:
    matplotlib_setup(dataset["output_dir"]); import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    root = output_dirs(dataset["output_dir"])["drift_correction"]
    report = root / "drift_correction_QC_report.pdf"
    with PdfPages(report) as pdf:
        fig, ax = plt.subplots(figsize=(11,8.5)); ax.axis("off"); ax.text(.02,.98, f"Drift correction QC — {dataset['name']}\n\n" + summary.to_string(index=False), va="top", family="monospace", fontsize=7); pdf.savefig(fig); plt.close(fig)
        for path in summary.qc_plot:
            image = plt.imread(path); fig, ax = plt.subplots(figsize=(11,8.5)); ax.imshow(image); ax.axis("off"); pdf.savefig(fig); plt.close(fig)
        for path in summary.qc_overlay:
            image = plt.imread(path); fig, ax = plt.subplots(figsize=(11,8.5)); ax.imshow(image); ax.axis("off"); pdf.savefig(fig); plt.close(fig)


def drift_correction_stage(cfg: dict[str, Any], dataset: dict[str, Any], manifest: pd.DataFrame, *, resume: bool = False, force: bool = False, max_files: int | None = None, max_frames: int | None = None) -> pd.DataFrame:
    rows = manifest[manifest.status == "accepted"].to_dict("records")[:max_files]
    root = output_dirs(dataset["output_dir"])["drift_correction"]; root.mkdir(parents=True, exist_ok=True); records = []
    minimum_free = float(cfg["qc"]["minimum_free_space_gib"])
    require_free_space(cfg["output_root"], minimum_free)
    with stage_timer(dataset["output_dir"], "02_drift_correction", {"files": len(rows), "max_frames": max_frames}):
        for index, row in enumerate(rows,1):
            print(f"[drift correction {index}/{len(rows)}] {row['filename']}", flush=True)
            records.append(correct_file(row, cfg, dataset, resume=resume, force=force, max_frames=max_frames))
        summary = pd.DataFrame.from_records(records); atomic_csv(summary, root / "drift_correction_QC_summary.csv"); save_stage_report(summary, dataset)
        warnings = summary[summary.status == "warning"]
        if len(warnings):
            raise RuntimeError("Drift correction QC failed for: " + ", ".join(warnings.filename.astype(str)))
        require_free_space(cfg["output_root"], minimum_free)
    return summary


def main(argv: list[str] | None = None) -> None:
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--config",required=True); parser.add_argument("--max-files",type=int); parser.add_argument("--max-frames",type=int)
    policy=parser.add_mutually_exclusive_group(); policy.add_argument("--resume",action="store_true"); policy.add_argument("--force",action="store_true")
    args=parser.parse_args(argv); cfg=load_config(args.config)
    for dataset in cfg["datasets"]: drift_correction_stage(cfg,dataset,load_manifest(dataset),resume=args.resume,force=args.force,max_files=args.max_files,max_frames=args.max_frames)


if __name__ == "__main__": main()
