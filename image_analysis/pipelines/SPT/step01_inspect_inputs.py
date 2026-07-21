"""Step 1: inspect ONI TIFFs, validate layout/metadata, and write manifests."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd
from tifffile import TiffFile

from spt_shared import atomic_csv, condition_from_fov, load_config, output_dirs, write_run_metadata


MANIFEST_COLUMNS = [
    "dataset", "path", "filename", "fov", "condition", "actual_frames",
    "expected_frames", "height", "width", "dtype", "pixel_size_um",
    "frame_interval_s", "exposure_ms", "layout", "status", "reason", "warnings",
]


def read_oni_metadata(tif: TiffFile) -> dict[str, Any]:
    description = tif.pages[0].description
    if not description:
        return {}
    try:
        value = json.loads(description)
        return value if isinstance(value, dict) else {}
    except json.JSONDecodeError:
        return {}


def inspect_tiff(path: str | Path, cfg: dict[str, Any]) -> dict[str, Any]:
    path = Path(path)
    warnings: list[str] = []
    with TiffFile(path) as tif:
        if not tif.pages:
            raise ValueError("TIFF contains no pages")
        page = tif.pages[0]
        height, width = page.shape[-2:]
        actual_frames = len(tif.pages)
        dtype = str(page.dtype)
        metadata = read_oni_metadata(tif)
    expected_frames = int(metadata.get("Frames", actual_frames))
    pixel_size = float(metadata.get("PixelSize_um", float("nan")))
    exposure_ms = float(metadata.get("Exposure_ms", float("nan")))
    fps = float(metadata.get("FramesPerSecond", 0) or 0)
    frame_interval = 1 / fps if fps > 0 else exposure_ms / 1000
    expected_width = int(cfg["layout"]["expected_width"])
    status = "accepted" if width == expected_width else "skipped"
    reason = "" if status == "accepted" else cfg["layout"].get("unexpected_width_reason", "unexpected_width")
    if actual_frames != expected_frames:
        warnings.append(f"page_count_mismatch:{actual_frames}!={expected_frames}")
    if abs(pixel_size - float(cfg["expected_pixel_size_um"])) > float(cfg["pixel_size_tolerance_um"]):
        warnings.append(f"pixel_size_unexpected:{pixel_size:g}")
    if not any(
        abs(frame_interval - float(value)) <= float(cfg["frame_interval_tolerance_s"])
        for value in cfg["allowed_frame_intervals_s"]
    ):
        warnings.append(f"frame_interval_unexpected:{frame_interval:g}")
    overrides = cfg["analysis"].get("condition_overrides", {})
    override = overrides.get(path.name, overrides.get(path.stem))
    if override is not None:
        condition, condition_warning = str(override), ""
    else:
        condition, condition_warning = condition_from_fov(path.stem, cfg["analysis"].get("condition_regex"))
    if condition_warning:
        warnings.append(condition_warning)
    return {
        "dataset": cfg["dataset"], "path": str(path.resolve()), "filename": path.name,
        "fov": path.stem, "condition": condition, "actual_frames": actual_frames,
        "expected_frames": expected_frames, "height": int(height), "width": int(width),
        "dtype": dtype, "pixel_size_um": pixel_size, "frame_interval_s": frame_interval,
        "exposure_ms": exposure_ms, "layout": cfg["layout"]["mode"], "status": status,
        "reason": reason, "warnings": ";".join(warnings),
    }


def build_manifest(cfg: dict[str, Any]) -> pd.DataFrame:
    records = []
    for path in sorted(Path(cfg["input_dir"]).glob(cfg.get("file_glob", "*.tif"))):
        try:
            records.append(inspect_tiff(path, cfg))
        except Exception as exc:
            records.append({
                "dataset": cfg["dataset"], "path": str(path.resolve()), "filename": path.name,
                "fov": path.stem, "status": "error", "reason": type(exc).__name__,
                "warnings": str(exc),
            })
    return pd.DataFrame.from_records(records, columns=MANIFEST_COLUMNS)


def inspect_inputs(
    config_path: str | Path, *, resume: bool = False, force: bool = False,
    batch_config_path: str | Path | None = None, command: str | None = None,
) -> tuple[dict[str, Any], pd.DataFrame]:
    cfg = load_config(config_path)
    write_run_metadata(
        cfg, resume=resume, force=force, batch_config_path=batch_config_path, command=command,
    )
    manifest = build_manifest(cfg)
    metadata = output_dirs(cfg["output_dir"])["metadata"]
    atomic_csv(manifest, metadata / "input_manifest.csv")
    atomic_csv(manifest[manifest["status"] != "accepted"], metadata / "unused_files.csv")
    return cfg, manifest


def accepted_rows(manifest: pd.DataFrame, max_files: int | None = None) -> list[dict[str, Any]]:
    accepted = manifest[manifest["status"] == "accepted"]
    if max_files is not None:
        accepted = accepted.head(max_files)
    return accepted.to_dict(orient="records")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)
    cfg, manifest = inspect_inputs(args.config)
    accepted = int((manifest["status"] == "accepted").sum())
    print(f"Accepted: {accepted}; skipped/errors: {len(manifest) - accepted}")
    print(output_dirs(cfg["output_dir"])["metadata"] / "input_manifest.csv")


if __name__ == "__main__":
    main()
