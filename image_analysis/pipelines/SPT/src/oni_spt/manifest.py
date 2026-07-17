from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
from tifffile import TiffFile


MANIFEST_COLUMNS = [
    "dataset",
    "path",
    "filename",
    "fov",
    "actual_frames",
    "expected_frames",
    "height",
    "width",
    "dtype",
    "pixel_size_um",
    "frame_interval_s",
    "exposure_ms",
    "layout",
    "status",
    "reason",
    "warnings",
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
            raise ValueError(f"TIFF contains no pages: {path}")
        page = tif.pages[0]
        height, width = page.shape[-2:]
        actual_frames = len(tif.pages)
        dtype = str(page.dtype)
        meta = read_oni_metadata(tif)
    expected_frames = int(meta.get("Frames", actual_frames))
    pixel_size = float(meta.get("PixelSize_um", float("nan")))
    exposure_ms = float(meta.get("Exposure_ms", float("nan")))
    fps = float(meta.get("FramesPerSecond", 0) or 0)
    frame_interval = 1.0 / fps if fps > 0 else exposure_ms / 1000.0
    expected_width = int(cfg["layout"]["expected_width"])
    if width != expected_width:
        status = "skipped"
        reason = cfg["layout"].get("unexpected_width_reason", "unexpected_width")
    else:
        status = "accepted"
        reason = ""
    if actual_frames != expected_frames:
        warnings.append(f"page_count_mismatch:{actual_frames}!={expected_frames}")
    px_expected = float(cfg["expected_pixel_size_um"])
    px_tol = float(cfg["pixel_size_tolerance_um"])
    if not (abs(pixel_size - px_expected) <= px_tol):
        warnings.append(f"pixel_size_unexpected:{pixel_size:g}")
    allowed_intervals = [float(v) for v in cfg["allowed_frame_intervals_s"]]
    interval_tol = float(cfg["frame_interval_tolerance_s"])
    if not any(abs(frame_interval - value) <= interval_tol for value in allowed_intervals):
        warnings.append(f"frame_interval_unexpected:{frame_interval:g}")
    return {
        "dataset": cfg["dataset"],
        "path": str(path.resolve()),
        "filename": path.name,
        "fov": path.stem,
        "actual_frames": actual_frames,
        "expected_frames": expected_frames,
        "height": int(height),
        "width": int(width),
        "dtype": dtype,
        "pixel_size_um": pixel_size,
        "frame_interval_s": frame_interval,
        "exposure_ms": exposure_ms,
        "layout": cfg["layout"]["mode"],
        "status": status,
        "reason": reason,
        "warnings": ";".join(warnings),
    }


def build_manifest(cfg: dict[str, Any]) -> pd.DataFrame:
    input_dir = Path(cfg["input_dir"])
    paths = sorted(input_dir.glob(cfg.get("file_glob", "*.tif")))
    records: list[dict[str, Any]] = []
    for path in paths:
        try:
            records.append(inspect_tiff(path, cfg))
        except Exception as exc:
            records.append(
                {
                    "dataset": cfg["dataset"],
                    "path": str(path.resolve()),
                    "filename": path.name,
                    "fov": path.stem,
                    "status": "error",
                    "reason": type(exc).__name__,
                    "warnings": str(exc),
                }
            )
    return pd.DataFrame.from_records(records, columns=MANIFEST_COLUMNS)


def write_manifest(manifest: pd.DataFrame, output_dir: str | Path) -> None:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(out / "input_manifest.csv", index=False)
    unused = manifest[manifest["status"] != "accepted"].copy()
    unused.to_csv(out / "unused_files.csv", index=False)


def accepted_rows(manifest: pd.DataFrame, max_files: int | None = None) -> list[dict[str, Any]]:
    rows = manifest[manifest["status"] == "accepted"]
    if max_files is not None:
        rows = rows.head(max_files)
    return rows.to_dict(orient="records")
