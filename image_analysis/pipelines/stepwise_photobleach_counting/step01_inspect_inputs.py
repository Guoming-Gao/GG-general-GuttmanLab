"""Step 1: inspect TIFF inputs and create acquisition-profile and storage QC."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Any

import pandas as pd
import tifffile

from pbsa_shared import acquisition_profile, atomic_csv, load_config, matplotlib_setup, output_dirs, stage_timer, write_run_metadata


def inspect_tiff(path: Path, cfg: dict[str, Any], dataset_name: str) -> dict[str, Any]:
    try:
        with tifffile.TiffFile(path) as tif:
            series = tif.series[0]
            shape = tuple(series.shape); dtype = str(series.dtype); axes = series.axes
            description = tif.pages[0].description or "{}"
            try: metadata = json.loads(description)
            except json.JSONDecodeError: metadata = {}
        pixel_size = metadata.get("PixelSize_um")
        expected_yx = tuple(cfg["expected_shape_yx"])
        reasons = []
        if len(shape) != 3 or tuple(shape[-2:]) != expected_yx: reasons.append("unexpected_shape")
        if dtype != cfg["expected_dtype"]: reasons.append("unexpected_dtype")
        if pixel_size is None: reasons.append("missing_pixel_size")
        elif abs(float(pixel_size) - float(cfg["expected_pixel_size_um"])) > float(cfg["pixel_size_tolerance_um"]):
            reasons.append("unexpected_pixel_size")
        return {
            "dataset": dataset_name, "filename": path.name, "fov": path.stem, "path": str(path),
            "status": "accepted" if not reasons else "rejected", "reason": ";".join(reasons),
            "frames": int(shape[0]) if len(shape) == 3 else 0,
            "height": int(shape[-2]), "width": int(shape[-1]), "axes": axes, "dtype": dtype,
            "size_GiB": path.stat().st_size / 2**30, "pixel_size_um": pixel_size,
            "exposure_ms": metadata.get("Exposure_ms"), "frames_per_second": metadata.get("FramesPerSecond"),
            "illumination_angle_deg": metadata.get("IlluminationAngle_deg"),
            "laser_wavelengths_nm": json.dumps(metadata.get("LaserWavelength_nm", [])),
            "laser_active": json.dumps(metadata.get("LaserActive", [])),
            "laser_power_percent": json.dumps(metadata.get("LaserPowerPercent", [])),
            "camera_binning": metadata.get("cameraBinning"), "acquisition_profile": acquisition_profile(metadata),
            "source_hash": metadata.get("Hash"), "start_time_ms": metadata.get("startTime_ms"),
        }
    except Exception as exc:
        return {"dataset": dataset_name, "filename": path.name, "fov": path.stem, "path": str(path),
                "status": "rejected", "reason": f"{type(exc).__name__}: {exc}", "frames": 0, "size_GiB": path.stat().st_size / 2**30}


def build_manifest(cfg: dict[str, Any], dataset: dict[str, Any]) -> pd.DataFrame:
    root = Path(dataset["input_dir"])
    files = sorted(root.glob(cfg["file_glob"]))
    if dataset.get("expected_files") is not None and len(files) != int(dataset["expected_files"]):
        raise RuntimeError(f"{dataset['name']}: expected {dataset['expected_files']} TIFFs, observed {len(files)}")
    return pd.DataFrame.from_records([inspect_tiff(path, cfg, dataset["name"]) for path in files])


def save_qc_report(manifest: pd.DataFrame, dataset: dict[str, Any]) -> Path:
    matplotlib_setup(dataset["output_dir"])
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    out = output_dirs(dataset["output_dir"])["input_inspection"]
    report = out / "input_inspection_QC_report.pdf"
    accepted = manifest[manifest.status == "accepted"]
    profiles = accepted.groupby("acquisition_profile", as_index=False).agg(
        files=("filename", "count"), frames=("frames", "sum"), input_GiB=("size_GiB", "sum"),
        exposure_ms=("exposure_ms", "first"), pixel_size_um=("pixel_size_um", "first"))
    atomic_csv(profiles, out / "acquisition_profiles.csv")
    free_gib = shutil.disk_usage(Path(dataset["output_dir"]).parent).free / 2**30
    with PdfPages(report) as pdf:
        fig, ax = plt.subplots(figsize=(11, 8.5)); ax.axis("off")
        text = [f"Input inspection QC — {dataset['name']}", "",
                f"TIFFs: {len(manifest)} ({len(accepted)} accepted, {(manifest.status != 'accepted').sum()} rejected)",
                f"Frames: {int(accepted.frames.sum()) if not accepted.empty else 0:,}",
                f"Input size: {accepted.size_GiB.sum():.2f} GiB", f"Free output space: {free_gib:.1f} GiB", "",
                "Acquisition profiles:", profiles.to_string(index=False)]
        if (manifest.status != "accepted").any(): text += ["", "Rejected:", manifest.loc[manifest.status != "accepted", ["filename", "reason"]].to_string(index=False)]
        ax.text(0.02, 0.98, "\n".join(text), va="top", family="monospace", fontsize=9); pdf.savefig(fig); plt.close(fig)
    return report


def inspect_dataset(cfg: dict[str, Any], dataset: dict[str, Any]) -> pd.DataFrame:
    write_run_metadata(cfg, dataset)
    out = output_dirs(dataset["output_dir"])["input_inspection"]; out.mkdir(parents=True, exist_ok=True)
    with stage_timer(dataset["output_dir"], "01_input_inspection"):
        manifest = build_manifest(cfg, dataset); atomic_csv(manifest, out / "input_manifest.csv"); save_qc_report(manifest, dataset)
    return manifest


def load_manifest(dataset: dict[str, Any]) -> pd.DataFrame:
    path = output_dirs(dataset["output_dir"])["input_inspection"] / "input_manifest.csv"
    if not path.exists(): raise FileNotFoundError(f"Run input inspection first: {path}")
    return pd.read_csv(path)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", required=True)
    args = parser.parse_args(argv); cfg = load_config(args.config)
    for dataset in cfg["datasets"]:
        table = inspect_dataset(cfg, dataset); print(f"{dataset['name']}: {(table.status == 'accepted').sum()}/{len(table)} accepted")


if __name__ == "__main__": main()
