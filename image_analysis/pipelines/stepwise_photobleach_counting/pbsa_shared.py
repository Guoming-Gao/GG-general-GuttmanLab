"""Shared configuration, paths, atomic I/O, and bookkeeping for the PBSA pipeline."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shutil
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
    "expected_shape_yx": [684, 428],
    "expected_dtype": "uint16",
    "expected_pixel_size_um": 0.117,
    "pixel_size_tolerance_um": 0.002,
    "drift_correction": {
        "block_frames": 5,
        "highpass_sigma_px": 8.0,
        "upsample_factor": 10,
        "max_abs_shift_px": 10.0,
        "interpolation_order": 1,
        "compression": "zlib",
        "compression_level": 1,
    },
    "spotiflow": {
        "pretrained_model": "general",
        "probability_threshold": 0.4,
        "min_distance": 1,
        "normalizer": "auto",
        "device": "auto",
        "detection_frames": 5,
        "model_cache_dir": "../BulkFluo_RDF/cache/spotiflow_models",
    },
    "roi_growth": {
        "smooth_sigma_px": 1.0,
        "foreground_sigma": 3.0,
        "maximum_radius_px": 50,
        "outer_background_width_px": 10,
        "point_max_equivalent_diameter_px": 5.0,
        "minimum_area_px": 2,
    },
    "trace_extraction": {
        "point_radius_px": 2.0,
        "point_background_inner_radius_px": 3.5,
        "point_background_outer_radius_px": 5.0,
        "point_minimum_neighbor_distance_px": 3.5,
        "mask_background_gap_px": 2.0,
        "mask_background_width_px": 2.0,
    },
    "quickpbsa": {
        "thresholds": {},
        "maxiter": 200,
        "num_cores": 4,
        "candidate_threshold_multipliers": [0.3, 0.4, 0.5, 0.6, 0.7],
        "maximum_validated_count": 40,
        "preliminary": {"crop": True, "bg_frames": 500},
        "filter": {"subtracted": True, "percentile_step": 90, "length_laststep": 20},
        "refinement": {},
    },
    "qc": {
        "max_montage_fovs": 12,
        "max_trace_panels_per_class": 24,
        "drift_median_residual_target_px": 0.2,
        "drift_p95_residual_target_px": 0.5,
        "minimum_free_space_gib": 50.0,
        "threshold_pilot_traces_per_fov": 24,
    },
}


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    result = json.loads(json.dumps(base))
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def load_config(path: str | Path) -> dict[str, Any]:
    config_path = Path(path).expanduser().resolve()
    with config_path.open() as handle:
        raw = yaml.safe_load(handle) or {}
    cfg = _deep_merge(DEFAULT_CONFIG, raw)
    if not cfg.get("output_root") or not cfg.get("datasets"):
        raise ValueError("Configuration requires output_root and datasets")
    cfg["output_root"] = str(Path(cfg["output_root"]).expanduser().resolve())
    cache = Path(cfg["spotiflow"]["model_cache_dir"]).expanduser()
    if not cache.is_absolute():
        cache = (config_path.parent / cache).resolve()
    cfg["spotiflow"]["model_cache_dir"] = str(cache)
    seen = set()
    datasets = []
    for item in cfg["datasets"]:
        dataset = dict(item)
        dataset["input_dir"] = str(Path(dataset["input_dir"]).expanduser().resolve())
        dataset.setdefault("name", Path(dataset["input_dir"]).name)
        if dataset["name"] in seen:
            raise ValueError(f"Duplicate dataset name: {dataset['name']}")
        seen.add(dataset["name"])
        dataset["output_dir"] = str(Path(cfg["output_root"]) / dataset["name"])
        datasets.append(dataset)
    cfg["datasets"] = datasets
    cfg["_config_path"] = str(config_path)
    validate_config(cfg)
    return cfg


def validate_config(cfg: dict[str, Any]) -> None:
    if float(cfg["spotiflow"]["probability_threshold"]) != 0.4:
        raise ValueError("spotiflow.probability_threshold must be 0.4 for this production analysis")
    radii = cfg["trace_extraction"]
    if not (float(radii["point_radius_px"]) < float(radii["point_background_inner_radius_px"])
            < float(radii["point_background_outer_radius_px"])):
        raise ValueError("Point ROI and background radii must increase monotonically")
    if int(cfg["drift_correction"]["block_frames"]) < 1:
        raise ValueError("drift_correction.block_frames must be positive")


def public_config(cfg: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in cfg.items() if not key.startswith("_")}


def config_fingerprint(cfg: dict[str, Any]) -> str:
    payload = json.dumps(public_config(cfg), sort_keys=True, default=str).encode()
    return hashlib.sha256(payload).hexdigest()[:16]


def output_dirs(output_dir: str | Path) -> dict[str, Path]:
    root = Path(output_dir)
    return {
        "root": root,
        "metadata": root / "00_run_metadata",
        "input_inspection": root / "01_input_inspection",
        "drift_correction": root / "02_drift_correction",
        "roi_detection": root / "03_roi_detection_and_growth",
        "traces": root / "04_background_corrected_traces",
        "threshold_pilot": root / "05_quickpbsa_threshold_pilot",
        "step_counts": root / "06_photobleaching_step_counts",
        "reports": root / "07_summary_reports",
    }


def ensure_output_dirs(output_dir: str | Path) -> dict[str, Path]:
    dirs = output_dirs(output_dir)
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def require_free_space(path: str | Path, minimum_gib: float) -> float:
    target = Path(path)
    target.mkdir(parents=True, exist_ok=True)
    free_gib = shutil.disk_usage(target).free / 2**30
    if free_gib < float(minimum_gib):
        raise RuntimeError(f"Only {free_gib:.1f} GiB free at {target}; at least {minimum_gib:.1f} GiB is required")
    return free_gib


def atomic_csv(table: pd.DataFrame, path: str | Path) -> None:
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    table.to_csv(temporary, index=False); os.replace(temporary, path)


def atomic_json(value: Any, path: str | Path) -> None:
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str))
    os.replace(temporary, path)


def atomic_text(value: str, path: str | Path) -> None:
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp"); temporary.write_text(value)
    os.replace(temporary, path)


def output_action(paths: list[Path], *, resume: bool, force: bool) -> str:
    existing = [path for path in paths if path.exists()]
    if not existing or force:
        return "run"
    if resume and len(existing) == len(paths):
        return "skip"
    raise FileExistsError(
        f"Output already exists: {existing[0]}. Use --resume for complete outputs or --force for this stage."
    )


def acquisition_profile(metadata_dict: dict[str, Any]) -> str:
    wavelengths = metadata_dict.get("LaserWavelength_nm", [])
    active = metadata_dict.get("LaserActive", [])
    powers = metadata_dict.get("LaserPowerPercent", [])
    selected = [(wavelengths[i], powers[i]) for i in range(min(len(wavelengths), len(active), len(powers))) if active[i]]
    wavelength, power = selected[0] if selected else ("unknown", "unknown")
    exposure = metadata_dict.get("Exposure_ms", "unknown")
    binning = metadata_dict.get("cameraBinning", "unknown")
    def compact(value: Any) -> str:
        try: return f"{float(value):g}".replace(".", "p")
        except (TypeError, ValueError): return re.sub(r"\W+", "", str(value))
    return f"{compact(wavelength)}nm_{compact(power)}pct_{compact(exposure)}ms_{compact(binning)}x"


CONDITION_LABELS = {
    "SHA_noDox": "SHA noDox",
    "SHA_Dox": "SHA Dox",
    "dSPEN_FL": "dSPEN FL",
    "dSPEN_dRRM": "dSPEN dRRM",
}


def condition_metadata(filename: str) -> dict[str, str | bool]:
    """Return biological condition and primary/appendix assignment from a TIFF name."""
    stem = Path(filename).stem
    standardized = stem.startswith("standarized-")
    if standardized:
        stem = stem[len("standarized-"):]
    condition_key = re.sub(r"[-_]*FOV(?:[-_]\d+)?$", "", stem, flags=re.IGNORECASE).rstrip("-_")
    return {
        "condition_key": condition_key,
        "condition": CONDITION_LABELS.get(condition_key, condition_key.replace("_", " ")),
        "analysis_set": "primary_standardized" if standardized else "acquisition_test",
        "is_primary_comparison": standardized,
    }


def dependency_versions() -> dict[str, str | None]:
    result = {}
    for name in ["numpy", "pandas", "scipy", "scikit-image", "tifffile", "imagecodecs", "spotiflow", "quickpbsa", "matplotlib", "PyYAML"]:
        try: result[name] = metadata.version(name)
        except metadata.PackageNotFoundError: result[name] = None
    return result


def write_run_metadata(cfg: dict[str, Any], dataset: dict[str, Any], command: str | None = None) -> None:
    dirs = ensure_output_dirs(dataset["output_dir"])
    source = Path(cfg["_config_path"])
    atomic_text(source.read_text(), dirs["metadata"] / "pipeline_config.yaml")
    resolved = public_config(cfg)
    atomic_text(yaml.safe_dump(resolved, sort_keys=False), dirs["metadata"] / "resolved_config.yaml")
    atomic_json({
        "config_fingerprint": config_fingerprint(cfg),
        "dataset": dataset["name"],
        "command": command or " ".join(sys.argv),
        "python": sys.version,
        "platform": platform.platform(),
        "dependencies": dependency_versions(),
        "created_unix": time.time(),
    }, dirs["metadata"] / "provenance.json")


@contextmanager
def stage_timer(output_dir: str | Path, stage_name: str, details: dict[str, Any] | None = None) -> Iterator[None]:
    path = output_dirs(output_dir)["metadata"] / "stage_status" / f"{stage_name}.json"
    started = time.time(); atomic_json({"stage": stage_name, "status": "running", "started_unix": started, **(details or {})}, path)
    try:
        yield
    except Exception as exc:
        atomic_json({"stage": stage_name, "status": "failed", "started_unix": started,
                     "finished_unix": time.time(), "error": f"{type(exc).__name__}: {exc}", **(details or {})}, path)
        raise
    atomic_json({"stage": stage_name, "status": "complete", "started_unix": started,
                 "finished_unix": time.time(), "elapsed_seconds": time.time() - started, **(details or {})}, path)


def matplotlib_setup(output_dir: str | Path) -> None:
    cache = Path(output_dir) / ".cache" / "matplotlib"; cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib
    matplotlib.use("Agg")
