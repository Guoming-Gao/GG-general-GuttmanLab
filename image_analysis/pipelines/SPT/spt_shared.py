"""Shared configuration, paths, atomic I/O, and run bookkeeping for ONI SPT."""

from __future__ import annotations

import copy
import hashlib
import json
import os
import platform
import re
import sys
import time
from contextlib import contextmanager
from importlib import metadata
from pathlib import Path
from typing import Any, Iterator

import pandas as pd
import yaml


DEFAULT_CONFIG: dict[str, Any] = {
    "file_glob": "*.tif",
    "expected_pixel_size_um": 0.117,
    "pixel_size_tolerance_um": 0.002,
    "allowed_frame_intervals_s": [0.03, 0.1],
    "frame_interval_tolerance_s": 0.003,
    "bandpass": {"sigma_low": 1.0, "sigma_high": 3.0},
    "segmentation": {
        "model": "cpsam",
        "source": "auto",
        "diameter": 80.0,
        "flow_threshold": 0.4,
        "cellprob_threshold": 0.0,
        "min_size": 15,
        "device": "auto",
    },
    "spotiflow": {
        "pretrained_model": "general",
        "probability_threshold": None,
        "min_distance": 1,
        "exclude_border": 1,
        "normalizer": "auto",
        "device": "auto",
    },
    "laptrack": {
        "metric": "sqeuclidean",
        "max_link_distance_px": 5.0,
        "gap_closing": True,
        "gap_closing_max_frame_count": 2,
        "min_track_length": 5,
    },
    "legacy": {"spot_radius_px": 2.5},
    "qc": {
        "max_video_frames": 100,
        "render_scale": 2,
        "video_fps": 10,
        "video_crf": 10,
        "track_tail_frames": 20,
        "track_overlay_dpi": 1200,
        "track_linewidth_pt": 0.1,
        "scale_bar_um": 1.0,
        "crop_padding_px": 12,
    },
    "analysis": {
        "condition_regex": None,
        "condition_overrides": {},
        "immobile_stepsize_nm": 30.0,
        "alpha_threshold": 0.7,
        "fit_r2_threshold": 0.7,
        "focal_depth_um": 0.7,
        "run_saspt": True,
    },
}


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def _resolve_path(value: str | Path, config_dir: Path) -> str:
    path = Path(value).expanduser()
    return str(path if path.is_absolute() else (config_dir / path).resolve())


def validate_config(cfg: dict[str, Any]) -> None:
    required = ("dataset", "input_dir", "output_dir", "layout", "channels")
    missing = [key for key in required if key not in cfg]
    if missing:
        raise ValueError(f"Missing required configuration keys: {', '.join(missing)}")
    mode = cfg["layout"].get("mode")
    if mode not in {"single", "dual"}:
        raise ValueError("layout.mode must be single or dual")
    expected = {"full"} if mode == "single" else {"left", "right"}
    if set(cfg["channels"]) != expected:
        raise ValueError(f"{mode} layout requires channel keys {sorted(expected)}")
    roles = {value.get("role") for value in cfg["channels"].values()}
    if not roles <= {"spt", "marker", "ignore"} or "spt" not in roles:
        raise ValueError("Channels must use spt/marker/ignore roles and include SPT")
    if float(cfg["bandpass"]["sigma_low"]) >= float(cfg["bandpass"]["sigma_high"]):
        raise ValueError("bandpass sigma_low must be smaller than sigma_high")
    if int(cfg["laptrack"]["gap_closing_max_frame_count"]) < 1:
        raise ValueError("gap_closing_max_frame_count must be at least one")


def load_config(path: str | Path) -> dict[str, Any]:
    config_path = Path(path).expanduser().resolve()
    with config_path.open() as handle:
        raw = yaml.safe_load(handle) or {}
    cfg = _deep_merge(DEFAULT_CONFIG, raw)
    cfg["input_dir"] = _resolve_path(cfg["input_dir"], config_path.parent)
    cfg["output_dir"] = _resolve_path(cfg["output_dir"], config_path.parent)
    cache = cfg["spotiflow"].get("model_cache_dir")
    cfg["spotiflow"]["model_cache_dir"] = _resolve_path(
        cache, config_path.parent
    ) if cache else str(Path(cfg["output_dir"]) / ".cache" / "spotiflow_models")
    cfg["_config_path"] = str(config_path)
    validate_config(cfg)
    return cfg


def public_config(cfg: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in cfg.items() if not key.startswith("_")}


def output_dirs(output_dir: str | Path) -> dict[str, Path]:
    root = Path(output_dir)
    return {
        "root": root,
        "metadata": root / "00_run_metadata",
        "raw": root / "01_preprocessed" / "raw_channels",
        "bandpass": root / "01_preprocessed" / "bandpass_float32",
        "mips": root / "01_preprocessed" / "mips",
        "segmentation": root / "02_segmentation",
        "detections": root / "03_spot_detection",
        "trajectories": root / "04_trajectories",
        "analysis": root / "05_diffusion_analysis",
        "reports": root / "06_reports",
        "qc": root / "07_quality_control",
    }


def ensure_output_dirs(output_dir: str | Path) -> dict[str, Path]:
    dirs = output_dirs(output_dir)
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def atomic_csv(table: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    table.to_csv(temporary, index=False)
    os.replace(temporary, path)


def atomic_json(value: Any, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, default=str)
    os.replace(temporary, path)


def output_action(paths: list[Path], *, resume: bool, force: bool) -> str:
    existing = [path for path in paths if path.exists()]
    if not existing:
        return "run"
    if force:
        return "run"
    if resume and len(existing) == len(paths):
        return "skip"
    raise FileExistsError(
        f"Output already exists: {existing[0]}. Use --resume for complete outputs "
        "or --force to replace this stage's outputs."
    )


def condition_from_fov(fov: str, condition_regex: str | None = None) -> tuple[str, str]:
    if condition_regex:
        match = re.search(condition_regex, fov)
        if not match:
            return fov, "condition_regex_no_match"
        value = match.groupdict().get("condition") or match.group(1)
        return value.rstrip("-_ "), ""
    match = re.search(r"(?i)(.*?)[-_ ]*FOV(?:[-_ ]*\d+)?", fov)
    if not match:
        return fov, "filename_has_no_FOV_token"
    return match.group(1).rstrip("-_ "), ""


def dependency_versions() -> dict[str, str | None]:
    names = [
        "numpy", "pandas", "scipy", "scikit-image", "tifffile", "spotiflow",
        "cellpose", "laptrack", "saspt", "imageio", "imageio-ffmpeg",
    ]
    result: dict[str, str | None] = {}
    for name in names:
        try:
            result[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            result[name] = None
    return result


def config_fingerprint(cfg: dict[str, Any]) -> str:
    payload = json.dumps(public_config(cfg), sort_keys=True, default=str).encode()
    return hashlib.sha256(payload).hexdigest()[:16]


def write_run_metadata(cfg: dict[str, Any]) -> None:
    dirs = ensure_output_dirs(cfg["output_dir"])
    with (dirs["metadata"] / "resolved_config.yaml").open("w") as handle:
        yaml.safe_dump(public_config(cfg), handle, sort_keys=False)
    atomic_json(
        {
            "created_unix": time.time(),
            "python": sys.version,
            "platform": platform.platform(),
            "config_fingerprint": config_fingerprint(cfg),
            "dependencies": dependency_versions(),
        },
        dirs["metadata"] / "provenance.json",
    )


@contextmanager
def stage_timer(cfg: dict[str, Any], stage: str, details: dict[str, Any]) -> Iterator[None]:
    path = output_dirs(cfg["output_dir"])["metadata"] / "stage_status.json"
    try:
        status = json.loads(path.read_text()) if path.exists() else {}
    except json.JSONDecodeError:
        status = {}
    entry = {"status": "running", "started_unix": time.time(), **details}
    status[stage] = entry
    atomic_json(status, path)
    try:
        yield
    except Exception as exc:
        entry.update(status="failed", finished_unix=time.time(), error=f"{type(exc).__name__}: {exc}")
        atomic_json(status, path)
        raise
    else:
        entry.update(status="complete", finished_unix=time.time())
        atomic_json(status, path)
