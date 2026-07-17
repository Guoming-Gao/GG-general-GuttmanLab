from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

import pandas as pd

from .config import load_config, write_resolved_config
from .detection import detect_file, load_spotiflow
from .manifest import accepted_rows, build_manifest, write_manifest
from .preprocess import preprocess_file
from .segmentation import segment_file
from .tracking import track_file
from .utils import atomic_json, provenance, stage_timer


def initialize_run(config_path: str | Path) -> tuple[dict[str, Any], pd.DataFrame]:
    cfg = load_config(config_path)
    output = Path(cfg["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    write_resolved_config(cfg, output / "run_config.resolved.yaml")
    manifest = build_manifest(cfg)
    write_manifest(manifest, output)
    atomic_json(provenance(cfg), output / "provenance.json")
    return cfg, manifest


def _rows(manifest: pd.DataFrame, max_files: int | None) -> list[dict[str, Any]]:
    rows = accepted_rows(manifest, max_files=max_files)
    if not rows:
        raise ValueError("No accepted TIFF inputs were found")
    return rows


def preprocess_stage(
    cfg: dict[str, Any],
    manifest: pd.DataFrame,
    *,
    max_files: int | None = None,
    max_frames: int | None = None,
    resume: bool = False,
    force: bool = False,
) -> None:
    status = Path(cfg["output_dir"]) / "stage_status.json"
    rows = _rows(manifest, max_files)
    with stage_timer(status, "preprocess", {"files": len(rows), "max_frames": max_frames}):
        for index, row in enumerate(rows, 1):
            print(f"[preprocess {index}/{len(rows)}] {row['filename']}", flush=True)
            preprocess_file(row, cfg, max_frames=max_frames, resume=resume, force=force)


def segment_stage(
    cfg: dict[str, Any],
    manifest: pd.DataFrame,
    *,
    max_files: int | None = None,
    resume: bool = False,
    force: bool = False,
    model: Any | None = None,
) -> None:
    status = Path(cfg["output_dir"]) / "stage_status.json"
    rows = _rows(manifest, max_files)
    with stage_timer(status, "segment", {"files": len(rows)}):
        for index, row in enumerate(rows, 1):
            print(f"[segment {index}/{len(rows)}] {row['filename']}", flush=True)
            segment_file(row, cfg, resume=resume, force=force, model=model)


def detect_stage(
    cfg: dict[str, Any],
    manifest: pd.DataFrame,
    *,
    max_files: int | None = None,
    resume: bool = False,
    force: bool = False,
    context: Any | None = None,
) -> None:
    status = Path(cfg["output_dir"]) / "stage_status.json"
    rows = _rows(manifest, max_files)
    context = context or load_spotiflow(cfg)
    with stage_timer(status, "detect", {"files": len(rows)}):
        for index, row in enumerate(rows, 1):
            print(f"[detect {index}/{len(rows)}] {row['filename']}", flush=True)
            detect_file(row, cfg, context, resume=resume, force=force)


def track_stage(
    cfg: dict[str, Any],
    manifest: pd.DataFrame,
    *,
    max_files: int | None = None,
    resume: bool = False,
    force: bool = False,
) -> None:
    status = Path(cfg["output_dir"]) / "stage_status.json"
    rows = _rows(manifest, max_files)
    with stage_timer(status, "track", {"files": len(rows)}):
        for index, row in enumerate(rows, 1):
            print(f"[track {index}/{len(rows)}] {row['filename']}", flush=True)
            track_file(row, cfg, resume=resume, force=force)


def run_all(
    config_path: str | Path,
    *,
    max_files: int | None = None,
    max_frames: int | None = None,
    resume: bool = False,
    force: bool = False,
) -> tuple[dict[str, Any], pd.DataFrame]:
    cfg, manifest = initialize_run(config_path)
    preprocess_stage(
        cfg, manifest, max_files=max_files, max_frames=max_frames, resume=resume, force=force
    )
    segment_stage(cfg, manifest, max_files=max_files, resume=resume, force=force)
    detect_stage(cfg, manifest, max_files=max_files, resume=resume, force=force)
    track_stage(cfg, manifest, max_files=max_files, resume=resume, force=force)
    return cfg, manifest
