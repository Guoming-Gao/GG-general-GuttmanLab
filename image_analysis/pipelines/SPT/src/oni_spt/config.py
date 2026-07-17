from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml


DEFAULTS: dict[str, Any] = {
    "file_glob": "*.tif",
    "expected_pixel_size_um": 0.117,
    "pixel_size_tolerance_um": 0.002,
    "allowed_frame_intervals_s": [0.03, 0.1],
    "frame_interval_tolerance_s": 0.003,
    "save": {"raw_channels": True, "filtered_spt": True},
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
        "gap_closing": False,
        "min_track_length": 5,
    },
    "legacy": {"spot_radius_px": 2.5},
    "qc": {"max_tracks_draw": 500, "max_detections_draw": 5000},
}


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def _resolve_path(value: str | Path, config_dir: Path) -> str:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (config_dir / path).resolve()
    return str(path)


def validate_config(cfg: dict[str, Any]) -> None:
    required = ("dataset", "input_dir", "output_dir", "layout", "channels")
    missing = [key for key in required if key not in cfg]
    if missing:
        raise ValueError(f"Missing required configuration keys: {', '.join(missing)}")
    layout = cfg["layout"]
    if layout.get("mode") not in {"single", "dual"}:
        raise ValueError("layout.mode must be 'single' or 'dual'")
    if not isinstance(layout.get("expected_width"), int) or layout["expected_width"] <= 0:
        raise ValueError("layout.expected_width must be a positive integer")
    channels = cfg["channels"]
    expected_keys = {"full"} if layout["mode"] == "single" else {"left", "right"}
    if set(channels) != expected_keys:
        raise ValueError(
            f"{layout['mode']} layout requires channel keys {sorted(expected_keys)}"
        )
    roles = {spec.get("role") for spec in channels.values()}
    if not roles <= {"spt", "marker", "ignore"} or "spt" not in roles:
        raise ValueError("Channel roles must be spt/marker/ignore and include at least one spt")
    low = float(cfg["bandpass"]["sigma_low"])
    high = float(cfg["bandpass"]["sigma_high"])
    if not 0 < low < high:
        raise ValueError("bandpass requires 0 < sigma_low < sigma_high")
    if int(cfg["laptrack"]["min_track_length"]) < 2:
        raise ValueError("laptrack.min_track_length must be at least 2")


def load_config(path: str | Path) -> dict[str, Any]:
    config_path = Path(path).expanduser().resolve()
    with config_path.open() as handle:
        raw = yaml.safe_load(handle) or {}
    cfg = _deep_merge(DEFAULTS, raw)
    cfg["input_dir"] = _resolve_path(cfg["input_dir"], config_path.parent)
    cfg["output_dir"] = _resolve_path(cfg["output_dir"], config_path.parent)
    cache = cfg["spotiflow"].get("model_cache_dir")
    if cache:
        cfg["spotiflow"]["model_cache_dir"] = _resolve_path(cache, config_path.parent)
    else:
        cfg["spotiflow"]["model_cache_dir"] = str(
            Path(cfg["output_dir"]) / ".cache" / "spotiflow_models"
        )
    cfg["_config_path"] = str(config_path)
    validate_config(cfg)
    return cfg


def public_config(cfg: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in cfg.items() if not key.startswith("_")}


def write_resolved_config(cfg: dict[str, Any], output_path: str | Path) -> None:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        yaml.safe_dump(public_config(cfg), handle, sort_keys=False)
