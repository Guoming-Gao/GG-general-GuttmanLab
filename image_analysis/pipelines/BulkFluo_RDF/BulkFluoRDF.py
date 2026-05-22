#!/usr/bin/env python
"""Bulk fluorescence per-hub RDF pipeline for paired SACD images.

The main analysis calculates H3K27ac radial intensity independently around each
SPEN hub using overlapping physical annular bins. Older Sofi-style summed-source
helpers are retained as reference utilities, but are not the default analysis
because a nucleus can contain many SPEN hubs and H3K27ac is not assumed to be
fully predicted by summed SPEN hub locations.

- object channel: SPEN `*-SACD-left.tif`, used to detect hub centers
- intensity channel: H3K27ac `*-SACD-right.tif`, used for nucleus masks and hub RDF
- nucleus: labeled mask defining valid pixels and nucleus-level reference intensity

Run with:
    conda run -n smlm python BulkFluoRDF.py run --config config.yaml
"""

from __future__ import annotations

import argparse
import copy
import os
import re
import shutil
from dataclasses import dataclass
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Iterable

import numpy as np
import pandas as pd
import tifffile
import yaml
from scipy import ndimage
from scipy.optimize import curve_fit
from scipy.spatial.distance import cdist
from skimage import exposure, feature, measure, morphology, segmentation

plt = None
pe = None
ListedColormap = None


DEFAULT_CONFIG = {
    "input_dir": "data",
    "output_dir": "results",
    "object_pattern": "*-SACD-left.tif",
    "intensity_pattern": "*-SACD-right.tif",
    "pixel_size_nm": 58.5,
    "radius_nm": 3000,
    "radius_px": 51,
    "max_fovs": 1,
    "cellpose": {
        "use_existing_masks": True,
        "diameter": 50,
        "flow_threshold": 0.4,
        "cellprob_threshold": 0.0,
        "device": "auto",
        "downsample": 1.0,
        "preprocess": {
            "method": "sacd_percentile",
            "lower_percentile": 1,
            "upper_percentile": 99.8,
            "background_percentile": 1,
            "gamma": 1.0,
            "gaussian_sigma": 0,
        },
    },
    "spotiflow": {
        "use_existing_spots": True,
        "pretrained_model": "general",
        "model_cache_dir": None,
        "probability_threshold": 0.4,
        "min_distance": 1,
        "exclude_border": 1,
        "normalizer": "auto",
        "device": "auto",
    },
    "nucleus_filter": {
        "enabled": True,
        "min_nucleus_area_px": 15000,
        "min_edge_distance_px": 20,
    },
    "rdf": {
        "mode": "per_hub_annular",
        "radius_nm": 3000,
        "bin_width_nm": 100,
        "bin_step_nm": 50,
        "normalization": "local_intensity_mean",
        "tail_normalization": {
            "enabled": False,
            "last_n_bins": 5,
        },
        "aggregation": {
            "plot_column": "h3k27ac_rdf_local_norm",
        },
    },
    "hub_filter": {
        "enabled": True,
        "metric": "spotiflow_intensity",
        "threshold_source": "nucleus_spen_median_plus_std",
        "std_multiplier": 2.0,
        "rule": "greater_equal",
        "hard_threshold": {
            "enabled": False,
            "source": "manual",
            "value": None,
            "control_dirs": [],
            "metric": "spotiflow_intensity",
            "statistic": "quantile",
            "quantile": 0.95,
        },
    },
    "hub_filter_comparison": {
        "enabled": True,
        "std_multipliers": [1, 2, 3],
        "output_root": "results/hub_filter_comparison",
    },
    "qc": {
        "make_overlays": True,
        "max_points_draw": 5000,
        "contrast_lower_percentile": 1,
        "contrast_upper_percentile": 99.8,
        "nucleus_crop_padding_px": 25,
        "hub_circle_size": 36,
        "scale_bar_um": 1.0,
    },
    "progress": {
        "enabled": True,
    },
    "detector_comparison": {
        "input_dir": "data",
        "object_pattern": "*-SACD-left.tif",
        "intensity_pattern": "*-SACD-right.tif",
        "max_fovs": 1,
        "domain": "nuclear",
        "output_dir": "results/detector_comparison",
        "ranking_metric": "area_integrated",
        "local_window_px": 15,
        "component_threshold_mad": 3,
        "top_n": 50,
        "spotiflow": {
            "probability_threshold": [0.5, 0.6, 0.7, 0.8, 0.9],
            "min_distance": [3],
        },
        "dog": {
            "min_sigma": [1, 2],
            "max_sigma": [4, 6, 8],
            "threshold_rel": [0.05, 0.1, 0.2, 0.3],
            "overlap": 0.5,
        },
    },
    "spotiflow_screen": {
        "input_dir": "data",
        "object_pattern": "*-SACD-left.tif",
        "intensity_pattern": "*-SACD-right.tif",
        "max_fovs": 1,
        "output_dir": "results/detector_comparison/spotiflow_center_intensity_screen",
        "thresholds": [0.2, 0.3, 0.4, 0.5],
        "min_distance": 1,
        "intensity_metric": "spotiflow_intensity",
        "selection_rule": "none",
        "selection_std_multiplier": 0.0,
        "min_nucleus_area_px": 15000,
        "min_edge_distance_px": 20,
        "top_n_per_nucleus": 50,
        "local_max_radii": [1, 2, 3],
        "circle_size": 85,
    },
}


@dataclass(frozen=True)
class ImagePair:
    fov: str
    object_path: Path
    intensity_path: Path


def load_config(path: str | Path = "config.yaml") -> dict:
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    path = Path(path)
    if path.exists():
        with path.open() as fh:
            user_cfg = yaml.safe_load(fh) or {}
        cfg = _deep_update(cfg, user_cfg)
    return cfg


def _deep_update(base: dict, update: dict) -> dict:
    out = dict(base)
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_update(out[key], value)
        else:
            out[key] = value
    return out


def fov_key(path: str | Path) -> str:
    name = Path(path).name
    name = re.sub(r"-(SACD|SACDpy)-(left|right)\.tiff?$", "", name, flags=re.IGNORECASE)
    return name


def pair_sacd_files(
    input_dir: str | Path,
    object_pattern: str = "*-SACD-left.tif",
    intensity_pattern: str = "*-SACD-right.tif",
    max_fovs: int | None = None,
) -> list[ImagePair]:
    input_dir = Path(input_dir)
    object_files = {fov_key(p): p for p in sorted(input_dir.glob(object_pattern))}
    intensity_files = {fov_key(p): p for p in sorted(input_dir.glob(intensity_pattern))}
    keys = sorted(set(object_files) & set(intensity_files))
    if max_fovs is not None:
        keys = keys[: int(max_fovs)]
    return [ImagePair(key, object_files[key], intensity_files[key]) for key in keys]


def raise_if_no_pairs(
    pairs: list[ImagePair],
    input_dir: str | Path,
    object_pattern: str,
    intensity_pattern: str,
) -> None:
    """Raise a clear pairing error before downstream empty DataFrames obscure it."""
    if pairs:
        return
    input_path = Path(input_dir)
    object_matches = sorted(input_path.glob(object_pattern))
    intensity_matches = sorted(input_path.glob(intensity_pattern))
    object_keys = {fov_key(p) for p in object_matches}
    intensity_keys = {fov_key(p) for p in intensity_matches}
    shared_keys = sorted(object_keys & intensity_keys)
    key_hint = ""
    if object_matches and intensity_matches and not shared_keys:
        object_examples = ", ".join(sorted(object_keys)[:3])
        intensity_examples = ", ".join(sorted(intensity_keys)[:3])
        key_hint = (
            " Matched files were found, but their inferred FOV keys did not overlap. "
            f"Example object keys: {object_examples or 'none'}. "
            f"Example intensity keys: {intensity_examples or 'none'}."
        )
    raise FileNotFoundError(
        "No paired SACD files found. "
        f"input_dir={input_path!s}; "
        f"object_pattern={object_pattern!r} matched {len(object_matches)} file(s); "
        f"intensity_pattern={intensity_pattern!r} matched {len(intensity_matches)} file(s). "
        f"{key_hint} "
        "Check the notebook parameter cell or use paired patterns like "
        "'*-SACD-left.tif'/'*-SACD-right.tif' or '*-SACDpy-left.tif'/'*-SACDpy-right.tif'."
    )


def resolve_cache_dir(out_dir: str | Path, explicit_cache_dir: str | Path | None = None, name: str = "cache") -> Path:
    """Return an explicit cache path or an output-local hidden cache path."""
    if explicit_cache_dir:
        return Path(explicit_cache_dir)
    return Path(out_dir) / ".cache" / name


def ensure_matplotlib() -> None:
    """Import Matplotlib after MPLCONFIGDIR has been set for the active run."""
    global plt, pe, ListedColormap
    if plt is not None:
        return
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as _plt
    from matplotlib import patheffects as _pe
    from matplotlib.colors import ListedColormap as _ListedColormap

    plt = _plt
    pe = _pe
    ListedColormap = _ListedColormap


def resolve_spotiflow_cache_dir(cfg: dict, out_dir: str | Path) -> Path:
    spotiflow_cfg = cfg.get("spotiflow", {})
    explicit_cache_dir = spotiflow_cfg.get("model_cache_dir") or spotiflow_cfg.get("cache_dir")
    if explicit_cache_dir:
        return Path(explicit_cache_dir)
    return Path(__file__).resolve().parent / "cache" / "spotiflow_models"


def configure_runtime_caches(cfg: dict, out_dir: str | Path) -> None:
    """Create output-local caches for libraries that need writable cache dirs."""
    out_path = Path(out_dir)
    mpl_cache = resolve_cache_dir(out_path, cfg.get("matplotlib_cache_dir"), "matplotlib")
    mpl_cache.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(mpl_cache)
    resolve_spotiflow_cache_dir(cfg, out_path).mkdir(parents=True, exist_ok=True)
    ensure_matplotlib()


def make_progress(enabled: bool = True):
    """Return a Rich progress bar when available; callers can fall back to prints."""
    if not enabled:
        return None
    try:
        from rich.progress import (
            BarColumn,
            MofNCompleteColumn,
            Progress,
            TextColumn,
            TimeElapsedColumn,
            TimeRemainingColumn,
        )
    except ImportError:
        return None
    return Progress(
        TextColumn("[bold blue]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        transient=False,
    )


def read_image(path: str | Path) -> np.ndarray:
    image = tifffile.imread(path)
    image = np.asarray(image)
    if image.ndim != 2:
        raise ValueError(f"Expected 2D SACD image for {path}, got shape {image.shape}")
    return image.astype(np.float32, copy=False)


def normalize01(image: np.ndarray, lower_percentile: float = 1, upper_percentile: float = 99.8) -> np.ndarray:
    """Robust display normalization for QC plots."""
    img = image.astype(np.float32, copy=False)
    finite = np.isfinite(img)
    if not finite.any():
        return np.zeros_like(img, dtype=np.float32)
    p_low, p_high = np.percentile(img[finite], (lower_percentile, upper_percentile))
    if p_high <= p_low:
        p_low, p_high = float(np.min(img[finite])), float(np.max(img[finite]))
    if p_high <= p_low:
        return np.zeros_like(img, dtype=np.float32)
    return np.clip((img - p_low) / (p_high - p_low), 0, 1).astype(np.float32)


def resolve_radius(cfg: dict) -> tuple[int, float, float]:
    """Return RDF radius in px/nm plus configured pixel size."""
    pixel_size_nm = float(cfg.get("pixel_size_nm", 58.5))
    if pixel_size_nm <= 0:
        raise ValueError("pixel_size_nm must be > 0")

    radius_nm = cfg.get("radius_nm")
    if radius_nm is not None:
        radius_nm = float(radius_nm)
        if radius_nm <= 0:
            raise ValueError("radius_nm must be > 0 when provided")
        radius_px = int(round(radius_nm / pixel_size_nm))
    else:
        radius_px = int(cfg.get("radius_px", 26))
        if radius_px <= 0:
            raise ValueError("radius_px must be > 0")
        radius_nm = radius_px * pixel_size_nm

    if radius_px <= 0:
        raise ValueError("Resolved radius_px must be > 0")
    return radius_px, radius_nm, pixel_size_nm


def normalize_img(image: np.ndarray) -> np.ndarray:
    """Match the CellposeSAM GUI normalization: float32 divided by max."""
    img = image.astype(np.float32)
    max_val = float(np.max(img))
    if max_val > 0:
        img /= max_val
    return img


def preprocess_cellpose_image(image: np.ndarray, preprocess_cfg: dict | None = None) -> np.ndarray:
    """Prepare SACD or conventional fluorescence image for CellposeSAM.

    `method=max` matches the existing CellposeSAM GUI. `method=sacd_percentile`
    is intended for SACD images with small floating-point values and occasional
    bright outliers.
    """
    cfg = preprocess_cfg or {}
    method = cfg.get("method", "max")
    img = image.astype(np.float32)

    sigma = float(cfg.get("gaussian_sigma", 0) or 0)
    if sigma > 0:
        img = ndimage.gaussian_filter(img, sigma=sigma).astype(np.float32)

    if method == "max":
        out = normalize_img(img)
    elif method == "sacd_percentile":
        finite = img[np.isfinite(img)]
        if finite.size == 0:
            out = np.zeros_like(img, dtype=np.float32)
        else:
            bg = np.percentile(finite, float(cfg.get("background_percentile", 1)))
            lo = np.percentile(finite, float(cfg.get("lower_percentile", 1)))
            hi = np.percentile(finite, float(cfg.get("upper_percentile", 99.8)))
            lo = min(lo, bg)
            img = img - bg
            hi = hi - bg
            lo = lo - bg
            if hi <= lo:
                out = normalize_img(img)
            else:
                out = np.clip((img - lo) / (hi - lo), 0, 1).astype(np.float32)
    else:
        raise ValueError("cellpose.preprocess.method must be 'max' or 'sacd_percentile'")

    gamma = float(cfg.get("gamma", 1.0) or 1.0)
    if gamma != 1.0:
        out = np.power(np.clip(out, 0, 1), gamma).astype(np.float32)
    return out


def run_cellpose_nuclei(
    intensity_image: np.ndarray,
    diameter: float = 50,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    device: str = "auto",
    preprocess: dict | None = None,
    downsample: float = 1.0,
) -> np.ndarray:
    try:
        from cellpose import core, models
    except ImportError as exc:
        raise RuntimeError("Cellpose is not installed in the active environment.") from exc

    gpu = False if device == "cpu" else core.use_gpu()
    print(
        f"  CellposeSAM nuclei: diameter={diameter}, gpu={gpu}, "
        f"preprocess={preprocess.get('method', 'max') if preprocess else 'max'}, "
        f"downsample={downsample}",
        flush=True,
    )
    model = models.CellposeModel(gpu=gpu)
    seg_image = preprocess_cellpose_image(intensity_image, preprocess)
    original_shape = seg_image.shape
    ds = float(downsample or 1.0)
    if ds <= 0:
        raise ValueError("cellpose.downsample must be > 0")
    if ds != 1.0:
        seg_image = ndimage.zoom(seg_image, zoom=ds, order=1)
        diameter = max(1.0, float(diameter) * ds)
    img = np.zeros((2,) + seg_image.shape, dtype=np.float32)
    img[0, :, :] = seg_image
    masks, _, _ = model.eval(
        img,
        batch_size=32,
        diameter=diameter,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
        normalize={"tile_norm_blocksize": 0},
    )
    masks = np.asarray(masks, dtype=np.uint16)
    if ds != 1.0 and masks.shape != original_shape:
        zoom = (original_shape[0] / masks.shape[0], original_shape[1] / masks.shape[1])
        masks = ndimage.zoom(masks, zoom=zoom, order=0).astype(np.uint16)
    return masks


def load_or_run_nuclei(pair: ImagePair, cfg: dict, out_dir: Path) -> np.ndarray:
    mask_rel = Path("nucleus_masks") / f"{pair.fov}_nuclei.tif"
    mask_path = out_dir / mask_rel
    if cfg["cellpose"].get("use_existing_masks", True):
        candidate_roots = [out_dir]
        if cfg.get("reuse_output_dir"):
            candidate_roots.append(Path(cfg["reuse_output_dir"]))
        seen = set()
        for root in candidate_roots:
            root = Path(root)
            if root in seen:
                continue
            seen.add(root)
            candidate = root / mask_rel
            if candidate.exists():
                return np.asarray(tifffile.imread(candidate))

    intensity = read_image(pair.intensity_path)
    masks = run_cellpose_nuclei(
        intensity,
        diameter=cfg["cellpose"]["diameter"],
        flow_threshold=cfg["cellpose"]["flow_threshold"],
        cellprob_threshold=cfg["cellpose"]["cellprob_threshold"],
        device=cfg["cellpose"].get("device", "auto"),
        preprocess=cfg["cellpose"].get("preprocess", {}),
        downsample=cfg["cellpose"].get("downsample", 1.0),
    )
    mask_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(mask_path, masks)
    return masks


def run_spotiflow_spots(
    object_image: np.ndarray,
    pretrained_model: str = "general",
    cache_dir: str | Path | None = None,
    probability_threshold: float | None = None,
    min_distance: int = 1,
    exclude_border: int = 1,
    normalizer: str = "auto",
    device: str = "auto",
) -> pd.DataFrame:
    cfg = {
        "spotiflow": {
            "pretrained_model": pretrained_model,
            "model_cache_dir": cache_dir,
            "device": device,
        }
    }
    context = load_spotiflow_context(cfg, Path("."))
    return predict_spotiflow_spots(
        object_image,
        context,
        probability_threshold=probability_threshold,
        min_distance=min_distance,
        exclude_border=exclude_border,
        normalizer=normalizer,
    )


def load_spotiflow_context(cfg: dict, out_dir: str | Path) -> SimpleNamespace:
    try:
        from spotiflow.model import Spotiflow
    except ImportError as exc:
        raise RuntimeError("Spotiflow is not installed in the active environment.") from exc

    spot_cfg = cfg.get("spotiflow", {})
    pretrained_model = spot_cfg.get("pretrained_model", "general")
    torch_device = resolve_torch_device(spot_cfg.get("device", "auto"))
    cache_path = resolve_spotiflow_cache_dir(cfg, out_dir)
    cache_path.mkdir(parents=True, exist_ok=True)
    print(
        f"  Spotiflow model: model={pretrained_model}, device={torch_device}, cache={cache_path}",
        flush=True,
    )
    model = Spotiflow.from_pretrained(
        pretrained_model,
        cache_dir=cache_path,
        map_location=str(torch_device),
    )
    return SimpleNamespace(
        model=model,
        device=torch_device,
        cache_dir=cache_path,
        pretrained_model=pretrained_model,
    )


def predict_spotiflow_spots(
    object_image: np.ndarray,
    context: SimpleNamespace,
    probability_threshold: float | None = None,
    min_distance: int = 1,
    exclude_border: int = 1,
    normalizer: str = "auto",
) -> pd.DataFrame:
    points, details = context.model.predict(
        object_image,
        prob_thresh=probability_threshold,
        min_distance=min_distance,
        exclude_border=exclude_border,
        normalizer=normalizer,
        device=str(context.device),
    )
    spots = pd.DataFrame(np.round(points, 4), columns=["y", "x"])
    if hasattr(details, "prob") and details.prob is not None:
        spots["probability"] = np.round(details.prob, 4)
    if hasattr(details, "intens") and details.intens is not None:
        spots["object_intensity"] = np.round(details.intens, 4)
    return spots


def resolve_torch_device(device: str):
    import torch

    device = str(device).lower()
    if device == "gpu":
        device = "auto"
    if device == "auto":
        if torch.backends.mps.is_available():
            return torch.device("mps")
        if torch.cuda.is_available():
            return torch.device("cuda")
        return torch.device("cpu")
    if device == "mps":
        if not torch.backends.mps.is_available():
            raise RuntimeError(
                "MPS/GPU was requested, but PyTorch reports MPS is unavailable. "
                "Check `python -c \"import torch; print(torch.backends.mps.is_available())\"` "
                "inside the smlm environment."
            )
        return torch.device("mps")
    if device == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA was requested, but PyTorch reports CUDA is unavailable.")
        return torch.device("cuda")
    if device == "cpu":
        return torch.device("cpu")
    raise ValueError("device must be one of: auto, mps, cuda, cpu")


def load_or_run_spots(
    pair: ImagePair,
    cfg: dict,
    out_dir: Path,
    spotiflow_context: SimpleNamespace | None = None,
    spotiflow_context_factory=None,
) -> pd.DataFrame:
    spot_dir = out_dir / "spotiflow_spots"
    spot_rel = Path("spotiflow_spots") / f"{pair.fov}_spots.csv"
    spot_path = out_dir / spot_rel
    if cfg["spotiflow"].get("use_existing_spots", True):
        candidate_roots = [out_dir]
        if cfg.get("reuse_output_dir"):
            candidate_roots.append(Path(cfg["reuse_output_dir"]))
        seen = set()
        for root in candidate_roots:
            root = Path(root)
            if root in seen:
                continue
            seen.add(root)
            candidate = root / spot_rel
            if candidate.exists():
                spots = pd.read_csv(candidate)
                spots.attrs["source_path"] = str(candidate)
                spots.attrs["loaded_from_cache"] = True
                return spots

    object_image = read_image(pair.object_path)
    if spotiflow_context is None:
        spotiflow_context = spotiflow_context_factory() if spotiflow_context_factory else load_spotiflow_context(cfg, out_dir)
    spots = predict_spotiflow_spots(
        object_image,
        spotiflow_context,
        probability_threshold=cfg["spotiflow"]["probability_threshold"],
        min_distance=cfg["spotiflow"]["min_distance"],
        exclude_border=cfg["spotiflow"]["exclude_border"],
        normalizer=cfg["spotiflow"]["normalizer"],
    )
    spots.insert(0, "fov", pair.fov)
    spot_dir.mkdir(parents=True, exist_ok=True)
    spots.to_csv(spot_path, index=False)
    spots.attrs["source_path"] = str(spot_path)
    spots.attrs["loaded_from_cache"] = False
    return spots


def assign_spots_to_nuclei(spots: pd.DataFrame, nuclei: np.ndarray) -> pd.DataFrame:
    if spots.empty:
        out = spots.copy()
        out["nucleus_id"] = pd.Series(dtype=int)
        return out
    yy = np.clip(np.rint(spots["y"].to_numpy()).astype(int), 0, nuclei.shape[0] - 1)
    xx = np.clip(np.rint(spots["x"].to_numpy()).astype(int), 0, nuclei.shape[1] - 1)
    out = spots.copy()
    out["nucleus_id"] = nuclei[yy, xx].astype(int)
    return out[out["nucleus_id"] > 0].reset_index(drop=True)


def _mad(values: np.ndarray) -> float:
    values = values.astype(float, copy=False)
    med = float(np.median(values))
    return float(np.median(np.abs(values - med)))


def measure_punctum_features(
    image: np.ndarray,
    y: float,
    x: float,
    pixel_size_nm: float,
    local_window_px: int = 15,
    component_threshold_mad: float = 3,
) -> dict:
    half = max(2, int(local_window_px))
    yc = int(np.clip(round(float(y)), 0, image.shape[0] - 1))
    xc = int(np.clip(round(float(x)), 0, image.shape[1] - 1))
    y0, y1 = max(0, yc - half), min(image.shape[0], yc + half + 1)
    x0, x1 = max(0, xc - half), min(image.shape[1], xc + half + 1)
    crop = image[y0:y1, x0:x1].astype(float, copy=False)

    local_background = float(np.median(crop))
    mad = _mad(crop)
    robust_sigma = 1.4826 * mad
    threshold = local_background + float(component_threshold_mad) * robust_sigma
    binary = crop > threshold

    cy, cx = yc - y0, xc - x0
    if not binary[cy, cx]:
        peak_rel = np.unravel_index(int(np.argmax(crop)), crop.shape)
        cy, cx = int(peak_rel[0]), int(peak_rel[1])
        binary[cy, cx] = True

    labels = measure.label(binary, connectivity=2)
    component_id = int(labels[cy, cx])
    component = labels == component_id if component_id > 0 else np.zeros_like(binary, dtype=bool)
    if not component.any():
        component = np.zeros_like(binary, dtype=bool)
        component[cy, cx] = True

    component_values = crop[component]
    area_px = int(component.sum())
    peak_intensity = float(np.max(component_values))
    integrated_intensity = float(np.sum(np.maximum(component_values - local_background, 0)))
    area_um2 = area_px * (float(pixel_size_nm) / 1000.0) ** 2
    equiv_diameter_px = float(2.0 * np.sqrt(area_px / np.pi))
    equiv_diameter_nm = equiv_diameter_px * float(pixel_size_nm)
    return {
        "peak_intensity": peak_intensity,
        "local_background": local_background,
        "background_corrected_integrated_intensity": integrated_intensity,
        "punctum_area_px": area_px,
        "punctum_area_um2": area_um2,
        "equivalent_diameter_px": equiv_diameter_px,
        "equivalent_diameter_nm": equiv_diameter_nm,
        "component_threshold": threshold,
        "local_window_px": half,
    }


def add_shared_spot_features(
    spots: pd.DataFrame,
    image: np.ndarray,
    nuclei: np.ndarray,
    pixel_size_nm: float,
    local_window_px: int,
    component_threshold_mad: float,
    domain: str = "nuclear",
) -> pd.DataFrame:
    if spots.empty:
        out = spots.copy()
        for col in [
            "nucleus_id",
            "peak_intensity",
            "local_background",
            "background_corrected_integrated_intensity",
            "punctum_area_px",
            "punctum_area_um2",
            "equivalent_diameter_px",
            "equivalent_diameter_nm",
            "component_threshold",
            "local_window_px",
        ]:
            out[col] = pd.Series(dtype=float)
        return out

    assigned = assign_spots_to_nuclei(spots, nuclei)
    if domain != "nuclear":
        yy = np.clip(np.rint(spots["y"].to_numpy()).astype(int), 0, nuclei.shape[0] - 1)
        xx = np.clip(np.rint(spots["x"].to_numpy()).astype(int), 0, nuclei.shape[1] - 1)
        assigned = spots.copy()
        assigned["nucleus_id"] = nuclei[yy, xx].astype(int)

    rows = []
    for row in assigned.to_dict("records"):
        feats = measure_punctum_features(
            image,
            row["y"],
            row["x"],
            pixel_size_nm,
            local_window_px=local_window_px,
            component_threshold_mad=component_threshold_mad,
        )
        row.update(feats)
        rows.append(row)
    return pd.DataFrame(rows)


def rdf_fit_original(
    intensity_image: np.ndarray,
    nucleus_mask: np.ndarray,
    centers_yx: np.ndarray,
    radius_px: int,
) -> tuple[np.ndarray, np.ndarray]:
    radii = np.arange(radius_px + 1, dtype=int)
    if centers_yx.size == 0:
        return np.full_like(radii, np.nan, dtype=float), np.zeros_like(radii, dtype=float)

    mask = nucleus_mask.copy()
    mask_inds = np.stack(np.where(mask), axis=-1)
    if mask_inds.size == 0:
        return np.full_like(radii, np.nan, dtype=float), np.zeros_like(radii, dtype=float)

    p_dists = cdist(mask_inds, centers_yx)
    p_dists[p_dists > radius_px] = -1
    p_dists = p_dists.round().astype(int)

    no_info = np.all(p_dists == -1, axis=-1)
    mask[mask] = ~no_info
    p_dists = p_dists[~no_info]

    if int(mask.sum()) < radius_px * 10:
        return np.full_like(radii, np.nan, dtype=float), np.zeros_like(radii, dtype=float)

    values, raw_counts = np.unique(p_dists, return_counts=True)
    counts = np.zeros(radius_px + 1, dtype=float)
    valid_values = values >= 0
    counts[values[valid_values]] = raw_counts[valid_values]

    img = intensity_image[mask].astype(float)
    im_mean = float(img.mean())
    im_std = float(img.std())
    if im_std == 0:
        return np.full(radius_px + 1, im_mean, dtype=float), counts
    img_scaled = (img - im_mean) / im_std

    def rdf_translate(x, *rdf_est):
        dists = p_dists[x.astype(int), :]
        params = np.append(rdf_est, 0.0)
        return params[dists].sum(axis=1)

    try:
        rdf, _ = curve_fit(
            rdf_translate,
            np.arange(p_dists.shape[0]),
            img_scaled,
            p0=np.zeros(radius_px + 1, dtype=float),
            maxfev=10000,
        )
        intensity = rdf * im_std + im_mean
    except Exception:
        design = np.zeros((p_dists.shape[0], radius_px + 1), dtype=float)
        rows = np.repeat(np.arange(p_dists.shape[0]), p_dists.shape[1])
        cols = p_dists.ravel()
        keep = cols >= 0
        np.add.at(design, (rows[keep], cols[keep]), 1.0)
        beta, *_ = np.linalg.lstsq(design, img_scaled, rcond=None)
        intensity = beta * im_std + im_mean
    return intensity, counts


def compute_rdf_for_fov(
    fov: str,
    intensity_image: np.ndarray,
    nuclei: np.ndarray,
    nuclear_spots: pd.DataFrame,
    radius_px: int,
    pixel_size_nm: float,
    nucleus_ids: Iterable[int] | None = None,
) -> pd.DataFrame:
    frames = []
    radii = np.arange(radius_px + 1, dtype=int)
    radii_nm = radii.astype(float) * float(pixel_size_nm)
    ids = sorted(int(v) for v in np.unique(nuclei) if v > 0) if nucleus_ids is None else sorted(int(v) for v in nucleus_ids)
    for nucleus_id in ids:
        nucleus_mask = nuclei == nucleus_id
        spots_nuc = nuclear_spots[nuclear_spots["nucleus_id"] == nucleus_id]
        centers = spots_nuc[["y", "x"]].to_numpy(float)
        intensity, counts = rdf_fit_original(intensity_image, nucleus_mask, centers, radius_px)
        frames.append(
            pd.DataFrame(
                {
                    "fov": fov,
                    "nucleus_id": nucleus_id,
                    "radius_px": radii,
                    "radius_nm": radii_nm,
                    "radius_um": radii_nm / 1000.0,
                    "h3k27ac_rdf_intensity": intensity,
                    "counts": counts,
                    "n_spen_hubs": len(spots_nuc),
                    "n_nucleus_pixels": int(nucleus_mask.sum()),
                    "status": "ok" if len(spots_nuc) > 0 else "no_spen_hubs",
                }
            )
        )
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def make_physical_rdf_bins(radius_nm: float, bin_width_nm: float, bin_step_nm: float) -> pd.DataFrame:
    if radius_nm <= 0 or bin_width_nm <= 0 or bin_step_nm <= 0:
        raise ValueError("radius_nm, bin_width_nm, and bin_step_nm must all be > 0")
    starts = np.arange(0, radius_nm - bin_width_nm + 1e-9, bin_step_nm, dtype=float)
    ends = starts + float(bin_width_nm)
    return pd.DataFrame(
        {
            "bin_index": np.arange(len(starts), dtype=int),
            "radius_start_nm": starts,
            "radius_end_nm": ends,
            "radius_center_nm": (starts + ends) / 2.0,
        }
    )


def _safe_series_stats(values: np.ndarray, prefix: str) -> dict:
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return {
            f"{prefix}_mean": np.nan,
            f"{prefix}_median": np.nan,
            f"{prefix}_std": np.nan,
            f"{prefix}_q90": np.nan,
            f"{prefix}_max": np.nan,
        }
    vals = values[finite]
    return {
        f"{prefix}_mean": float(np.mean(vals)),
        f"{prefix}_median": float(np.median(vals)),
        f"{prefix}_std": float(np.std(vals)),
        f"{prefix}_q90": float(np.quantile(vals, 0.9)),
        f"{prefix}_max": float(np.max(vals)),
    }


def validate_spotiflow_intensity_scale(
    hub_properties: pd.DataFrame,
    fov: str,
    spot_source_path: str | Path | None = None,
    min_median_ratio: float = 0.01,
    max_median_ratio: float = 100.0,
) -> None:
    """Catch stale or normalized Spotiflow intensity values before hub filtering silently fails."""
    required = {"spotiflow_intensity", "spen_center_intensity"}
    if hub_properties.empty or not required.issubset(hub_properties.columns):
        return
    metric = hub_properties["spotiflow_intensity"].astype(float).to_numpy()
    center = hub_properties["spen_center_intensity"].astype(float).to_numpy()
    valid = np.isfinite(metric) & np.isfinite(center) & (metric > 0) & (center > 0)
    if valid.sum() < 20:
        return
    ratio = float(np.median(metric[valid]) / np.median(center[valid]))
    if min_median_ratio <= ratio <= max_median_ratio:
        return
    source = f" Source spot CSV: {spot_source_path}." if spot_source_path else ""
    raise ValueError(
        f"Spotiflow intensity scale looks inconsistent for {fov}: "
        f"median(spotiflow_intensity)/median(raw SPEN center intensity) = {ratio:.4g}. "
        "This usually means a stale spot CSV from a differently scaled run was reused. "
        "Delete/regenerate that FOV's spot CSV or set spotiflow.use_existing_spots=false."
        f"{source}"
    )


def _resolve_relative_path(path: str | Path, config_path: str | Path | None = None) -> Path:
    out = Path(path).expanduser()
    if out.is_absolute():
        return out
    if config_path is not None:
        return Path(config_path).expanduser().resolve().parent / out
    return out


def _read_control_metric_values(control_dirs: Iterable[str | Path], metric: str, config_path: str | Path | None = None) -> np.ndarray:
    values = []
    missing = []
    for control_dir in control_dirs:
        control_path = _resolve_relative_path(control_dir, config_path)
        csv_path = control_path if control_path.is_file() else control_path / "hub_properties.csv"
        if not csv_path.exists():
            missing.append(str(csv_path))
            continue
        df = pd.read_csv(csv_path)
        if metric not in df.columns:
            raise ValueError(f"Hard hub threshold metric `{metric}` is missing from {csv_path}")
        metric_values = pd.to_numeric(df[metric], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        if not metric_values.empty:
            values.append(metric_values.to_numpy(float))
    if missing:
        raise FileNotFoundError("Missing hard-threshold control hub_properties.csv file(s): " + "; ".join(missing))
    if not values:
        return np.array([], dtype=float)
    return np.concatenate(values)


def resolve_hub_filter_hard_threshold(hub_filter_cfg: dict, config_path: str | Path | None = None) -> dict:
    hard_cfg = hub_filter_cfg.get("hard_threshold", {}) or {}
    enabled = bool(hard_cfg.get("enabled", False))
    metric = hard_cfg.get("metric", hub_filter_cfg.get("metric", "spotiflow_intensity"))
    source = hard_cfg.get("source", "manual")
    statistic = hard_cfg.get("statistic", "quantile")
    quantile = hard_cfg.get("quantile", np.nan)
    resolved = {
        "enabled": enabled,
        "value": np.nan,
        "source": source,
        "metric": metric,
        "statistic": statistic,
        "quantile": float(quantile) if quantile is not None and np.isfinite(float(quantile)) else np.nan,
        "n_control_values": 0,
        "control_dirs": list(hard_cfg.get("control_dirs", []) or hard_cfg.get("negative_control_dirs", []) or []),
    }
    if not enabled:
        return resolved

    manual_value = hard_cfg.get("value", hard_cfg.get("threshold", None))
    if manual_value is not None:
        value = float(manual_value)
        if not np.isfinite(value):
            raise ValueError("hub_filter.hard_threshold.value must be finite when provided")
        resolved.update({"value": value, "source": source or "manual", "statistic": "manual", "n_control_values": 0})
        return resolved

    if source != "negative_control_result_dirs":
        raise ValueError(
            "hub_filter.hard_threshold.source must be `negative_control_result_dirs` "
            "or provide hub_filter.hard_threshold.value"
        )

    control_dirs = resolved["control_dirs"]
    if not control_dirs:
        raise ValueError("hub_filter.hard_threshold.control_dirs must list negative-control result directories")
    values = _read_control_metric_values(control_dirs, metric, config_path)
    if values.size == 0:
        raise ValueError("No finite values found for hard hub threshold controls")

    if statistic == "quantile":
        q = float(hard_cfg.get("quantile", 0.95))
        if q < 0 or q > 1:
            raise ValueError("hub_filter.hard_threshold.quantile must be between 0 and 1")
        value = float(np.quantile(values, q))
        resolved["quantile"] = q
    elif statistic == "max":
        value = float(np.max(values))
    else:
        raise ValueError("hub_filter.hard_threshold.statistic must be `quantile` or `max`")

    resolved.update({"value": value, "n_control_values": int(values.size)})
    return resolved


def combine_hub_filter_threshold(local_threshold: float, hard_threshold_info: dict) -> float:
    thresholds = []
    if np.isfinite(local_threshold):
        thresholds.append(float(local_threshold))
    if hard_threshold_info.get("enabled", False):
        hard_value = hard_threshold_info.get("value", np.nan)
        if np.isfinite(hard_value):
            thresholds.append(float(hard_value))
    return float(max(thresholds)) if thresholds else np.nan


def resolve_hub_filter_threshold(row: dict, hub_filter_cfg: dict) -> float:
    source = hub_filter_cfg.get("threshold_source", "nucleus_spen_median_plus_std")
    if source == "nucleus_spen_median_plus_std":
        median = row.get("nucleus_spen_median", np.nan)
        std = row.get("nucleus_spen_std", np.nan)
        multiplier = float(hub_filter_cfg.get("std_multiplier", 1.0))
        return float(median + multiplier * std) if np.isfinite(median) and np.isfinite(std) else np.nan
    return row.get(source, np.nan)


def compute_hub_properties_and_annular_rdf(
    fov: str,
    object_image: np.ndarray,
    intensity_image: np.ndarray,
    nuclei: np.ndarray,
    nuclear_spots: pd.DataFrame,
    keep_nucleus_ids: Iterable[int],
    bins: pd.DataFrame,
    pixel_size_nm: float,
    hub_filter_cfg: dict,
    rdf_cfg: dict | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Measure per-hub SPEN properties and H3K27ac annular RDF.

    This is the main RDF mode for SPEN/H3K27ac. It does not use the Sofi
    summed-source model. Each Spotiflow point is treated independently, and each
    annular bin is clipped to the hub's parent nucleus.
    """
    keep_ids = set(int(v) for v in keep_nucleus_ids)
    spots = nuclear_spots[nuclear_spots["nucleus_id"].isin(keep_ids)].copy().reset_index(drop=True)
    if spots.empty:
        return pd.DataFrame(), pd.DataFrame()

    if "object_intensity" in spots.columns and "spotiflow_intensity" not in spots.columns:
        spots["spotiflow_intensity"] = spots["object_intensity"].astype(float)
    elif "spotiflow_intensity" not in spots.columns:
        spots["spotiflow_intensity"] = np.nan

    filter_enabled = bool(hub_filter_cfg.get("enabled", True))
    filter_metric = hub_filter_cfg.get("metric", "spotiflow_intensity")
    threshold_source = hub_filter_cfg.get("threshold_source", "nucleus_spen_median_plus_std")
    hard_threshold_info = hub_filter_cfg.get("_resolved_hard_threshold") or resolve_hub_filter_hard_threshold(hub_filter_cfg)
    hard_threshold_enabled = bool(hard_threshold_info.get("enabled", False))
    hard_threshold_value = float(hard_threshold_info.get("value", np.nan))
    rdf_cfg = rdf_cfg or {}
    tail_cfg = rdf_cfg.get("tail_normalization", {})
    tail_enabled = bool(tail_cfg.get("enabled", False))
    tail_last_n = max(1, int(tail_cfg.get("last_n_bins", 5)))
    local_radius_nm = float(bins["radius_end_nm"].iloc[0]) if not bins.empty else 200.0
    max_radius_nm = float(bins["radius_end_nm"].max()) if not bins.empty else float(rdf_cfg.get("radius_nm", 1000))

    hub_rows = []
    rdf_rows = []
    hub_counter = 0
    for nucleus_id in sorted(spots["nucleus_id"].dropna().astype(int).unique()):
        nucleus_mask = nuclei == nucleus_id
        yy, xx = np.where(nucleus_mask)
        if yy.size == 0:
            continue
        spen_vals = object_image[nucleus_mask].astype(float)
        h3k_vals = intensity_image[nucleus_mask].astype(float)
        nucleus_stats = {
            **_safe_series_stats(spen_vals, "nucleus_spen"),
            **_safe_series_stats(h3k_vals, "nucleus_h3k27ac"),
            "nucleus_area_px": int(yy.size),
            "nucleus_area_um2": float(yy.size * (pixel_size_nm / 1000.0) ** 2),
        }
        nucleus_h3k_mean = nucleus_stats["nucleus_h3k27ac_mean"]
        spots_nuc = spots[spots["nucleus_id"] == nucleus_id].reset_index(drop=True)

        for _, spot in spots_nuc.iterrows():
            hub_counter += 1
            hub_y = float(spot["y"])
            hub_x = float(spot["x"])
            dist_nm = np.sqrt((yy.astype(float) - hub_y) ** 2 + (xx.astype(float) - hub_x) ** 2) * float(pixel_size_nm)
            local_mask = dist_nm <= local_radius_nm
            local_spen = spen_vals[local_mask]
            hub_stats = _safe_series_stats(local_spen, f"spen_local_r{int(round(local_radius_nm))}nm")
            local_reference_mask = dist_nm < max_radius_nm
            local_h3k_vals = h3k_vals[local_reference_mask]
            finite_local_h3k = np.isfinite(local_h3k_vals)
            local_h3k_mean = float(np.mean(local_h3k_vals[finite_local_h3k])) if finite_local_h3k.any() else np.nan
            local_spen_vals = spen_vals[local_reference_mask]
            finite_local_spen = np.isfinite(local_spen_vals)
            local_spen_mean = float(np.mean(local_spen_vals[finite_local_spen])) if finite_local_spen.any() else np.nan
            local_reference_pixel_count = int(finite_local_h3k.sum())
            center_y = int(np.clip(round(hub_y), 0, object_image.shape[0] - 1))
            center_x = int(np.clip(round(hub_x), 0, object_image.shape[1] - 1))
            row = {
                "fov": fov,
                "hub_id": hub_counter,
                "hub_uid": f"{fov}_hub_{hub_counter:05d}",
                "nucleus_id": int(nucleus_id),
                "hub_y": hub_y,
                "hub_x": hub_x,
                "spotiflow_probability": float(spot["probability"]) if "probability" in spot and pd.notna(spot["probability"]) else np.nan,
                "spotiflow_intensity": float(spot["spotiflow_intensity"]) if pd.notna(spot["spotiflow_intensity"]) else np.nan,
                "spen_center_intensity": float(object_image[center_y, center_x]),
                **nucleus_stats,
                **hub_stats,
                "spen_local_area_px": int(local_mask.sum()),
                "local_h3k27ac_mean": local_h3k_mean,
                "local_spen_mean": local_spen_mean,
                "local_reference_pixel_count": local_reference_pixel_count,
            }
            local_threshold_value = resolve_hub_filter_threshold(row, hub_filter_cfg)
            threshold_value = combine_hub_filter_threshold(local_threshold_value, hard_threshold_info)
            metric_value = row.get(filter_metric, np.nan)
            if filter_enabled:
                row["hub_filter_metric"] = filter_metric
                row["hub_filter_threshold_source"] = threshold_source
                row["hub_filter_std_multiplier"] = float(hub_filter_cfg.get("std_multiplier", np.nan))
                row["hub_filter_local_threshold"] = local_threshold_value
                row["hub_filter_hard_threshold_enabled"] = hard_threshold_enabled
                row["hub_filter_hard_threshold"] = hard_threshold_value
                row["hub_filter_hard_threshold_source"] = hard_threshold_info.get("source", "")
                row["hub_filter_hard_threshold_metric"] = hard_threshold_info.get("metric", "")
                row["hub_filter_hard_threshold_statistic"] = hard_threshold_info.get("statistic", "")
                row["hub_filter_hard_threshold_quantile"] = hard_threshold_info.get("quantile", np.nan)
                row["hub_filter_hard_threshold_n_control_values"] = int(hard_threshold_info.get("n_control_values", 0))
                row["hub_filter_threshold"] = threshold_value
                row["hub_filter_final_threshold"] = threshold_value
                row["keep_hub"] = bool(np.isfinite(metric_value) and np.isfinite(threshold_value) and metric_value >= threshold_value)
            else:
                row["hub_filter_metric"] = filter_metric
                row["hub_filter_threshold_source"] = "filter_disabled"
                row["hub_filter_std_multiplier"] = np.nan
                row["hub_filter_local_threshold"] = np.nan
                row["hub_filter_hard_threshold_enabled"] = False
                row["hub_filter_hard_threshold"] = np.nan
                row["hub_filter_hard_threshold_source"] = "filter_disabled"
                row["hub_filter_hard_threshold_metric"] = ""
                row["hub_filter_hard_threshold_statistic"] = ""
                row["hub_filter_hard_threshold_quantile"] = np.nan
                row["hub_filter_hard_threshold_n_control_values"] = 0
                row["hub_filter_threshold"] = np.nan
                row["hub_filter_final_threshold"] = np.nan
                row["keep_hub"] = True
            hub_rows.append(row)

            if not row["keep_hub"]:
                continue
            hub_bin_rows = []
            for bin_row in bins.to_dict("records"):
                in_bin = (dist_nm >= float(bin_row["radius_start_nm"])) & (dist_nm < float(bin_row["radius_end_nm"]))
                pixel_count = int(in_bin.sum())
                h3k_bin_vals = h3k_vals[in_bin]
                spen_bin_vals = spen_vals[in_bin]
                if pixel_count > 0:
                    h3k_sum = float(np.sum(h3k_bin_vals))
                    expected_h3k_sum = float(local_h3k_mean * pixel_count) if np.isfinite(local_h3k_mean) else np.nan
                    h3k_rdf = h3k_sum / expected_h3k_sum if expected_h3k_sum and np.isfinite(expected_h3k_sum) else np.nan
                    h3k_mean = float(np.mean(h3k_bin_vals))
                    h3k_median = float(np.median(h3k_bin_vals))
                    spen_sum = float(np.sum(spen_bin_vals))
                    expected_spen_sum = float(local_spen_mean * pixel_count) if np.isfinite(local_spen_mean) else np.nan
                    spen_rdf = spen_sum / expected_spen_sum if expected_spen_sum and np.isfinite(expected_spen_sum) else np.nan
                    spen_mean = float(np.mean(spen_bin_vals))
                    spen_median = float(np.median(spen_bin_vals))
                else:
                    h3k_sum = np.nan
                    expected_h3k_sum = np.nan
                    h3k_rdf = np.nan
                    h3k_mean = np.nan
                    h3k_median = np.nan
                    spen_sum = np.nan
                    expected_spen_sum = np.nan
                    spen_rdf = np.nan
                    spen_mean = np.nan
                    spen_median = np.nan
                hub_bin_rows.append(
                    {
                        "fov": fov,
                        "hub_id": hub_counter,
                        "hub_uid": row["hub_uid"],
                        "nucleus_id": int(nucleus_id),
                        "radius_start_nm": float(bin_row["radius_start_nm"]),
                        "radius_end_nm": float(bin_row["radius_end_nm"]),
                        "radius_center_nm": float(bin_row["radius_center_nm"]),
                        "bin_index": int(bin_row["bin_index"]),
                        "h3k27ac_sum_intensity": h3k_sum,
                        "h3k27ac_mean_intensity": h3k_mean,
                        "h3k27ac_median_intensity": h3k_median,
                        "expected_h3k27ac_sum_intensity": expected_h3k_sum,
                        "local_h3k27ac_mean": local_h3k_mean,
                        "spen_sum_intensity": spen_sum,
                        "spen_mean_intensity": spen_mean,
                        "spen_median_intensity": spen_median,
                        "expected_spen_sum_intensity": expected_spen_sum,
                        "local_spen_mean": local_spen_mean,
                        "local_reference_pixel_count": local_reference_pixel_count,
                        "h3k27ac_rdf_local_norm": h3k_rdf,
                        "spen_rdf_local_norm": spen_rdf,
                        "pixel_count": pixel_count,
                        "spotiflow_intensity": row["spotiflow_intensity"],
                        "nucleus_spen_median": row["nucleus_spen_median"],
                        "nucleus_h3k27ac_mean": nucleus_h3k_mean,
                    }
                )
            if tail_enabled:
                h3k_values = np.array([r["h3k27ac_rdf_local_norm"] for r in hub_bin_rows], dtype=float)
                spen_values = np.array([r["spen_rdf_local_norm"] for r in hub_bin_rows], dtype=float)
                finite_h3k = np.isfinite(h3k_values)
                finite_spen = np.isfinite(spen_values)
                h3k_tail_values = h3k_values[finite_h3k][-tail_last_n:] if finite_h3k.any() else np.array([], dtype=float)
                spen_tail_values = spen_values[finite_spen][-tail_last_n:] if finite_spen.any() else np.array([], dtype=float)
                h3k_tail_baseline = float(np.mean(h3k_tail_values)) if h3k_tail_values.size > 0 else np.nan
                spen_tail_baseline = float(np.mean(spen_tail_values)) if spen_tail_values.size > 0 else np.nan
                if not np.isfinite(h3k_tail_baseline) or h3k_tail_baseline == 0:
                    h3k_tail_baseline = np.nan
                if not np.isfinite(spen_tail_baseline) or spen_tail_baseline == 0:
                    spen_tail_baseline = np.nan
                for rdf_row in hub_bin_rows:
                    rdf_row["h3k27ac_tail_baseline_last5"] = h3k_tail_baseline
                    rdf_row["spen_tail_baseline_last5"] = spen_tail_baseline
                    h3k_raw = rdf_row["h3k27ac_rdf_local_norm"]
                    spen_raw = rdf_row["spen_rdf_local_norm"]
                    rdf_row["h3k27ac_rdf_tail_norm"] = h3k_raw / h3k_tail_baseline if np.isfinite(h3k_raw) and np.isfinite(h3k_tail_baseline) else np.nan
                    rdf_row["spen_rdf_tail_norm"] = spen_raw / spen_tail_baseline if np.isfinite(spen_raw) and np.isfinite(spen_tail_baseline) else np.nan
            rdf_rows.extend(hub_bin_rows)
    return pd.DataFrame(hub_rows), pd.DataFrame(rdf_rows)


def aggregate_hub_rdf(
    hub_rdf: pd.DataFrame,
    value_column: str = "h3k27ac_rdf_local_norm",
) -> pd.DataFrame:
    if hub_rdf.empty:
        return pd.DataFrame()
    if value_column not in hub_rdf.columns:
        value_column = "h3k27ac_rdf_local_norm"
    if value_column.startswith("spen_"):
        prefix = "spen_rdf"
    elif value_column.startswith("h3k27ac_"):
        prefix = "h3k27ac_rdf"
    else:
        prefix = "rdf"
    rows = []
    group_cols = ["bin_index", "radius_start_nm", "radius_end_nm", "radius_center_nm"]
    for name, group in hub_rdf.groupby(group_cols, sort=True):
        values = group[value_column].replace([np.inf, -np.inf], np.nan).dropna().astype(float)
        mean = float(values.mean()) if not values.empty else np.nan
        std = float(values.std(ddof=1)) if len(values) > 1 else np.nan
        rows.append(
            {
                "bin_index": int(name[0]),
                "radius_start_nm": float(name[1]),
                "radius_end_nm": float(name[2]),
                "radius_center_nm": float(name[3]),
                "aggregation_value_column": value_column,
                f"{value_column}_mean": mean,
                f"{value_column}_std": std,
                f"{prefix}_mean": mean,
                f"{prefix}_std": std,
                "n_hubs": int(values.size),
                "total_pixel_count": int(group["pixel_count"].sum()),
                f"{prefix}_bin_intensity_mean": float(group["h3k27ac_mean_intensity"].mean()) if value_column.startswith("h3k27ac_") else float(group["spen_mean_intensity"].mean()),
            }
        )
    return pd.DataFrame(rows)


def aggregate_dual_channel_rdf(
    hub_rdf: pd.DataFrame,
    h3k27ac_column: str = "h3k27ac_rdf_local_norm",
    spen_column: str = "spen_rdf_local_norm",
) -> pd.DataFrame:
    h3 = aggregate_hub_rdf(hub_rdf, h3k27ac_column)
    spen = aggregate_hub_rdf(hub_rdf, spen_column)
    if h3.empty:
        return spen
    if spen.empty:
        return h3
    keys = ["bin_index", "radius_start_nm", "radius_end_nm", "radius_center_nm"]
    spen_keep = keys + [c for c in spen.columns if c.startswith("spen_rdf")]
    return h3.merge(spen[spen_keep], on=keys, how="left")


def clamp_rdf_ylim(ax, lower: float = 0.5, upper: float = 2.0) -> None:
    ymin, ymax = ax.get_ylim()
    if ymin < lower or ymax > upper:
        ax.set_ylim(lower, upper)


def plot_rdf_mean_with_std(
    ax,
    x: np.ndarray,
    mean: np.ndarray,
    std: np.ndarray,
    color,
    label: str,
    zorder: int = 4,
) -> None:
    ax.errorbar(
        x,
        mean,
        yerr=std,
        fmt=".",
        color=color,
        ecolor=color,
        elinewidth=1,
        capsize=2,
        alpha=0.5,
        label=label,
        zorder=zorder,
    )
    ax.plot(x, mean, color=color, lw=2, zorder=zorder)


def make_hub_aggregate_rdf_plot(
    aggregate_summary: pd.DataFrame,
    out_dir: Path,
    hub_rdf: pd.DataFrame | None = None,
    value_column: str = "h3k27ac_rdf_local_norm",
) -> Path | None:
    ensure_matplotlib()
    if aggregate_summary.empty:
        return None
    out_path = out_dir / "rdf_aggregate.png"
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True, sharex=True)
    plot_specs = [
        (axes[0], "h3k27ac_rdf_local_norm", "h3k27ac_rdf", "local-intensity normalized H3K27ac RDF", "H3K27ac RDF"),
        (axes[1], "spen_rdf_local_norm", "spen_rdf", "local-intensity normalized SPEN RDF", "SPEN positive-control RDF"),
    ]
    cmap = plt.get_cmap("viridis")
    for ax, plot_column, prefix, ylabel, title in plot_specs:
        if hub_rdf is not None and not hub_rdf.empty and plot_column in hub_rdf.columns:
            draw = hub_rdf[np.isfinite(hub_rdf[plot_column])].copy()
            intens = draw.groupby("hub_uid")["spotiflow_intensity"].first()
            if not intens.empty:
                vmin = float(intens.min())
                vmax = float(intens.max())
                denom = vmax - vmin if vmax > vmin else 1.0
            else:
                vmin = 0.0
                denom = 1.0
            for hub_uid, curve in draw.groupby("hub_uid"):
                hub_intensity = float(curve["spotiflow_intensity"].iloc[0])
                color = cmap((hub_intensity - vmin) / denom)
                ax.plot(
                    curve["radius_start_nm"],
                    curve[plot_column],
                    color=color,
                    alpha=0.08,
                    lw=0.55,
                    zorder=1,
                )
        x = aggregate_summary["radius_start_nm"].to_numpy(float)
        y = aggregate_summary[f"{prefix}_mean"].to_numpy(float)
        std = aggregate_summary[f"{prefix}_std"].to_numpy(float)
        plot_rdf_mean_with_std(ax, x, y, std, "black", "Mean +/- STD", zorder=4)
        ax.axhline(1.0, color="0.45", lw=1, ls="--")
        ax.set_xlim(float(aggregate_summary["radius_start_nm"].min()), float(aggregate_summary["radius_start_nm"].max()))
        ax.set_xlabel("Distance from SPEN hub center (nm)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        clamp_rdf_ylim(ax)
        n_hubs = int(aggregate_summary["n_hubs"].max()) if "n_hubs" in aggregate_summary else 0
        ax.plot([], [], linestyle="None", label=f"N = {n_hubs} SPEN hubs")
        ax.legend(frameon=False, fontsize=8, loc="best")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def add_micrometer_scale_bar(ax, crop_shape: tuple[int, int], pixel_size_nm: float, scale_bar_um: float = 1.0) -> None:
    ensure_matplotlib()
    bar_px = float(scale_bar_um) * 1000.0 / float(pixel_size_nm)
    margin = max(8, int(min(crop_shape) * 0.05))
    x1 = crop_shape[1] - margin
    x0 = x1 - bar_px
    y = crop_shape[0] - margin
    ax.plot([x0, x1], [y, y], color="white", lw=3, solid_capstyle="butt")
    ax.text(
        (x0 + x1) / 2,
        y - max(4, int(margin * 0.35)),
        f"{scale_bar_um:g} um",
        color="white",
        ha="center",
        va="bottom",
        fontsize=8,
        fontweight="bold",
        path_effects=[pe.withStroke(linewidth=2, foreground="black")],
    )


def make_single_hub_qc_overlay(
    pair: ImagePair,
    object_display: np.ndarray,
    intensity_display: np.ndarray,
    boundary_crop: np.ndarray,
    hubs_nuc: pd.DataFrame,
    kept: pd.DataFrame,
    rejected: pd.DataFrame,
    color,
    crop_bounds: tuple[int, int, int, int],
    out_dir: Path,
    circle_size: float,
    pixel_size_nm: float,
    scale_bar_um: float,
) -> Path:
    ensure_matplotlib()
    y0, y1, x0, x1 = crop_bounds
    nucleus_id = int(hubs_nuc["nucleus_id"].iloc[0]) if not hubs_nuc.empty else 0
    out_path = out_dir / f"{pair.fov}_nucleus_{nucleus_id:03d}_qc_overlay.png"
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.2), constrained_layout=True)

    for ax, display, title, cmap in (
        (axes[0], object_display, "SPEN", "magma"),
        (axes[1], intensity_display, "H3K27ac", "gray"),
    ):
        ax.imshow(display[y0:y1, x0:x1], cmap=cmap)
        ax.imshow(np.ma.masked_where(~boundary_crop, boundary_crop), cmap=ListedColormap([color]), alpha=0.9)
        if not rejected.empty:
            ax.scatter(
                rejected["hub_x"] - x0,
                rejected["hub_y"] - y0,
                s=circle_size,
                facecolors="none",
                edgecolors="0.65",
                linewidths=0.7,
                alpha=0.65,
            )
        if not kept.empty:
            ax.scatter(
                kept["hub_x"] - x0,
                kept["hub_y"] - y0,
                s=circle_size,
                facecolors="none",
                edgecolors=[color],
                linewidths=1.2,
            )
        add_micrometer_scale_bar(ax, boundary_crop.shape, pixel_size_nm, scale_bar_um)
        ax.set_title(title, fontsize=10)
        ax.set_axis_off()

    fig.suptitle(f"{pair.fov} nucleus {nucleus_id:03d}: retained={len(kept)}, rejected={len(rejected)}", fontsize=10)
    fig.savefig(out_path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    return out_path


def make_hub_qc_overlay(
    pair: ImagePair,
    object_image: np.ndarray,
    intensity_image: np.ndarray,
    nuclei: np.ndarray,
    hub_properties: pd.DataFrame,
    hub_rdf: pd.DataFrame,
    out_dir: Path,
    max_points_draw: int = 5000,
    contrast_lower_percentile: float = 1,
    contrast_upper_percentile: float = 99.8,
    nucleus_crop_padding_px: int = 25,
    hub_circle_size: float = 36,
    pixel_size_nm: float = 58.5,
    scale_bar_um: float = 1.0,
    rdf_plot_column: str = "h3k27ac_rdf_local_norm",
) -> Path:
    ensure_matplotlib()
    qc_dir = out_dir / "qc_overlays"
    qc_dir.mkdir(parents=True, exist_ok=True)
    single_overlay_dir = qc_dir / "overlays"
    single_overlay_dir.mkdir(parents=True, exist_ok=True)
    out_path = qc_dir / f"{pair.fov}_qc.png"

    object_display = normalize01(object_image, contrast_lower_percentile, contrast_upper_percentile)
    intensity_display = normalize01(intensity_image, contrast_lower_percentile, contrast_upper_percentile)
    nucleus_ids = sorted(int(v) for v in hub_properties["nucleus_id"].dropna().unique()) if not hub_properties.empty else []
    if not nucleus_ids:
        nucleus_ids = sorted(int(v) for v in np.unique(nuclei) if v > 0)

    n_tiles = max(1, len(nucleus_ids))
    montage_cols = max(1, int(math.ceil(math.sqrt(n_tiles))))
    montage_rows = int(math.ceil(n_tiles / montage_cols))
    fig_width = max(14, 4.2 * montage_cols + 6)
    fig_height = max(6, 3.0 * montage_rows)
    fig = plt.figure(figsize=(fig_width, fig_height), constrained_layout=True)
    outer = fig.add_gridspec(1, 2, width_ratios=[montage_cols, 1.35])
    tile_gs = outer[0, 0].subgridspec(montage_rows, montage_cols)
    rdf_ax = fig.add_subplot(outer[0, 1])
    colors = plt.get_cmap("tab20")(np.linspace(0, 1, max(1, len(nucleus_ids))))
    color_by_nucleus = {nid: colors[i % len(colors)] for i, nid in enumerate(nucleus_ids)}
    padding = int(nucleus_crop_padding_px)
    circle_size = float(hub_circle_size)

    for row_idx, nucleus_id in enumerate(nucleus_ids):
        tile_row = row_idx // montage_cols
        tile_col = row_idx % montage_cols
        pair_gs = tile_gs[tile_row, tile_col].subgridspec(1, 2, wspace=0.02)
        spen_ax = fig.add_subplot(pair_gs[0, 0])
        h3k_ax = fig.add_subplot(pair_gs[0, 1])
        mask = nuclei == nucleus_id
        yy, xx = np.where(mask)
        if yy.size == 0:
            for ax in (spen_ax, h3k_ax):
                ax.set_axis_off()
            continue
        y0 = max(0, int(yy.min()) - padding)
        y1 = min(nuclei.shape[0], int(yy.max()) + padding + 1)
        x0 = max(0, int(xx.min()) - padding)
        x1 = min(nuclei.shape[1], int(xx.max()) + padding + 1)
        boundary_crop = segmentation.find_boundaries(mask, mode="outer")[y0:y1, x0:x1]
        hubs_nuc = hub_properties[hub_properties["nucleus_id"] == nucleus_id].head(max_points_draw)
        kept = hubs_nuc[hubs_nuc["keep_hub"]]
        rejected = hubs_nuc[~hubs_nuc["keep_hub"]]
        color = color_by_nucleus[nucleus_id]
        make_single_hub_qc_overlay(
            pair,
            object_display,
            intensity_display,
            boundary_crop,
            hubs_nuc,
            kept,
            rejected,
            color,
            (y0, y1, x0, x1),
            single_overlay_dir,
            circle_size,
            pixel_size_nm,
            scale_bar_um,
        )

        for ax, display, title in (
            (spen_ax, object_display, "SPEN"),
            (h3k_ax, intensity_display, "H3K27ac"),
        ):
            ax.imshow(display[y0:y1, x0:x1], cmap="magma" if title == "SPEN" else "gray")
            ax.imshow(np.ma.masked_where(~boundary_crop, boundary_crop), cmap=ListedColormap([color]), alpha=0.9)
            if not rejected.empty:
                ax.scatter(
                    rejected["hub_x"] - x0,
                    rejected["hub_y"] - y0,
                    s=circle_size,
                    facecolors="none",
                    edgecolors="0.65",
                    linewidths=0.7,
                    alpha=0.65,
                )
            if not kept.empty:
                ax.scatter(
                    kept["hub_x"] - x0,
                    kept["hub_y"] - y0,
                    s=circle_size,
                    facecolors="none",
                    edgecolors=[color],
                    linewidths=0.9,
                )
            ax.set_title(f"N{nucleus_id} {title} ({len(kept)}/{len(hubs_nuc)} hubs)", fontsize=9)
            ax.set_axis_off()

    for empty_idx in range(len(nucleus_ids), montage_rows * montage_cols):
        tile_row = empty_idx // montage_cols
        tile_col = empty_idx % montage_cols
        ax = fig.add_subplot(tile_gs[tile_row, tile_col])
        ax.set_axis_off()

    if not hub_rdf.empty:
        plot_column = rdf_plot_column if rdf_plot_column in hub_rdf.columns else "h3k27ac_rdf_local_norm"
        draw = hub_rdf[np.isfinite(hub_rdf[plot_column])].copy()
        for hub_uid, curve in draw.groupby("hub_uid"):
            nid = int(curve["nucleus_id"].iloc[0])
            rdf_ax.plot(
                curve["radius_start_nm"],
                curve[plot_column],
                color=color_by_nucleus.get(nid, "0.5"),
                alpha=0.18,
                lw=0.8,
            )
        summary = aggregate_hub_rdf(draw, value_column=plot_column)
        if not summary.empty:
            std = summary["h3k27ac_rdf_std"].to_numpy(float)
            mean = summary["h3k27ac_rdf_mean"].to_numpy(float)
            x = summary["radius_start_nm"].to_numpy(float)
            plot_rdf_mean_with_std(rdf_ax, x, mean, std, "black", "FOV mean +/- STD", zorder=4)
    rdf_ax.axhline(1.0, color="0.45", lw=1, ls="--")
    if not hub_rdf.empty:
        rdf_ax.set_xlim(float(hub_rdf["radius_start_nm"].min()), float(hub_rdf["radius_start_nm"].max()))
    rdf_ax.set_xlabel("Distance from SPEN hub center (nm)")
    rdf_ax.set_ylabel("local-intensity normalized H3K27ac RDF")
    rdf_ax.set_title("Per-hub RDF in retained nuclei")
    clamp_rdf_ylim(rdf_ax)
    rdf_ax.plot([], [], linestyle="None", label=f"N = {int(hub_rdf['hub_uid'].nunique()) if not hub_rdf.empty else 0} SPEN hubs")
    rdf_ax.legend(loc="best", fontsize=8, frameon=False)
    fig.suptitle(pair.fov)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def write_pairs_csv(pairs: Iterable[ImagePair], out_dir: Path) -> pd.DataFrame:
    df = pd.DataFrame(
        [{"fov": p.fov, "object_path": str(p.object_path), "intensity_path": str(p.intensity_path)} for p in pairs]
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "paired_inputs.csv", index=False)
    return df


def run_pipeline(config_path: str | Path = "config.yaml", cfg: dict | None = None) -> SimpleNamespace:
    cfg = load_config(config_path) if cfg is None else cfg
    out_dir = Path(cfg["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    configure_runtime_caches(cfg, out_dir)
    radius_px, radius_nm, pixel_size_nm = resolve_radius(cfg)
    rdf_cfg = cfg.get("rdf", {})
    rdf_radius_nm = float(rdf_cfg.get("radius_nm", radius_nm))
    aggregation_cfg = rdf_cfg.get("aggregation", {})
    rdf_plot_column = aggregation_cfg.get("plot_column", "h3k27ac_rdf_local_norm")
    hub_filter_cfg = copy.deepcopy(cfg.get("hub_filter", {}))
    hard_threshold_info = resolve_hub_filter_hard_threshold(hub_filter_cfg, config_path)
    hub_filter_cfg["_resolved_hard_threshold"] = hard_threshold_info
    bins = make_physical_rdf_bins(
        rdf_radius_nm,
        float(rdf_cfg.get("bin_width_nm", 100)),
        float(rdf_cfg.get("bin_step_nm", 50)),
    )

    pairs = pair_sacd_files(
        cfg["input_dir"],
        cfg["object_pattern"],
        cfg["intensity_pattern"],
        cfg.get("max_fovs"),
    )
    raise_if_no_pairs(pairs, cfg["input_dir"], cfg["object_pattern"], cfg["intensity_pattern"])
    pair_df = write_pairs_csv(pairs, out_dir)
    make_qc = bool(cfg.get("qc", {}).get("make_overlays", True))
    per_fov_steps = 6 if make_qc else 5
    progress = make_progress(bool(cfg.get("progress", {}).get("enabled", True)))
    progress_task = None
    spotiflow_model_cache = resolve_spotiflow_cache_dir(cfg, out_dir)
    spotiflow_context = None

    def log(message: str) -> None:
        if progress is not None:
            progress.console.print(message)
        else:
            print(message, flush=True)

    def set_progress(description: str) -> None:
        if progress is not None and progress_task is not None:
            progress.update(progress_task, description=description)

    def advance_progress(description: str | None = None) -> None:
        if progress is not None and progress_task is not None:
            if description is not None:
                progress.update(progress_task, description=description)
            progress.advance(progress_task)

    def get_spotiflow_context() -> SimpleNamespace:
        nonlocal spotiflow_context
        if spotiflow_context is None:
            spotiflow_context = load_spotiflow_context(cfg, out_dir)
        return spotiflow_context

    if progress is not None:
        progress.start()
        progress_task = progress.add_task("Starting pipeline", total=len(pairs) * per_fov_steps + 1)

    try:
        log(f"Paired {len(pairs)} FOV(s)")
        log(
            f"Per-hub RDF radius: {rdf_radius_nm:g} nm; "
            f"{len(bins)} overlapping bins "
            f"({float(rdf_cfg.get('bin_width_nm', 100)):g} nm width, "
            f"{float(rdf_cfg.get('bin_step_nm', 50)):g} nm step)"
        )
        log(
            f"Normalization: {rdf_cfg.get('normalization', 'local_intensity_mean')}; "
            f"aggregate value={rdf_plot_column}"
        )
        if hard_threshold_info.get("enabled", False):
            log(
                "Hub hard threshold: "
                f"{hard_threshold_info.get('metric')} >= {float(hard_threshold_info.get('value')):g} "
                f"({hard_threshold_info.get('statistic')}, n={int(hard_threshold_info.get('n_control_values', 0))})"
            )
        log(f"Spotiflow model cache: {spotiflow_model_cache}")

        all_hub_properties = []
        all_hub_rdf = []
        all_nucleus_qc = []
        qc_paths = []
        for index, pair in enumerate(pairs, start=1):
            log(f"Processing {pair.fov} ({index}/{len(pairs)})")
            set_progress(f"{pair.fov}: loading images")
            object_image = read_image(pair.object_path)
            intensity_image = read_image(pair.intensity_path)
            log("  Loaded object/intensity images")
            advance_progress()

            set_progress(f"{pair.fov}: nuclei")
            nuclei = load_or_run_nuclei(pair, cfg, out_dir)
            log(f"  Nuclei: {max(0, len(np.unique(nuclei)) - 1)}")
            advance_progress()

            set_progress(f"{pair.fov}: nucleus filter")
            filter_cfg = cfg.get("nucleus_filter", {})
            if filter_cfg.get("enabled", True):
                nucleus_qc = compute_nucleus_qc(
                    nuclei,
                    pixel_size_nm,
                    int(filter_cfg.get("min_nucleus_area_px", 15000)),
                    int(filter_cfg.get("min_edge_distance_px", 20)),
                )
                nucleus_qc.insert(0, "fov", pair.fov)
                keep_ids = nucleus_qc.loc[nucleus_qc["keep_nucleus"], "nucleus_id"].astype(int).tolist()
            else:
                keep_ids = sorted(int(v) for v in np.unique(nuclei) if v > 0)
                nucleus_qc = pd.DataFrame(
                    {
                        "fov": pair.fov,
                        "nucleus_id": keep_ids,
                        "keep_nucleus": True,
                        "exclude_reason": "filter_disabled",
                    }
                )
            all_nucleus_qc.append(nucleus_qc)
            log(f"  Retained nuclei: {len(keep_ids)}")
            advance_progress()

            set_progress(f"{pair.fov}: Spotiflow hubs")
            spots = load_or_run_spots(pair, cfg, out_dir, spotiflow_context_factory=get_spotiflow_context)
            spot_source_path = spots.attrs.get("source_path")
            log(f"  Raw SPEN hubs: {len(spots)}")
            nuclear_spots = assign_spots_to_nuclei(spots, nuclei)
            nuclear_spots = nuclear_spots[nuclear_spots["nucleus_id"].isin(keep_ids)].reset_index(drop=True)
            log(f"  Nuclear SPEN hubs: {len(nuclear_spots)}")
            (out_dir / "spotiflow_spots").mkdir(parents=True, exist_ok=True)
            nuclear_spots.to_csv(out_dir / "spotiflow_spots" / f"{pair.fov}_nuclear_spots.csv", index=False)
            advance_progress()

            set_progress(f"{pair.fov}: RDF")
            hub_properties, hub_rdf = compute_hub_properties_and_annular_rdf(
                pair.fov,
                object_image,
                intensity_image,
                nuclei,
                nuclear_spots,
                keep_ids,
                bins,
                pixel_size_nm,
                hub_filter_cfg,
                rdf_cfg,
            )
            validate_spotiflow_intensity_scale(hub_properties, pair.fov, spot_source_path)
            threshold_source = hub_filter_cfg.get("threshold_source", "nucleus_spen_median_plus_std")
            std_multiplier = hub_filter_cfg.get("std_multiplier")
            filter_label = threshold_source if std_multiplier is None else f"{threshold_source}; std_multiplier={float(std_multiplier):g}"
            if hard_threshold_info.get("enabled", False):
                filter_label = f"{filter_label}; hard>={float(hard_threshold_info.get('value')):g}"
            log(
                f"  Retained SPEN hubs after filter ({filter_label}): "
                f"{int(hub_properties['keep_hub'].sum()) if not hub_properties.empty else 0}"
            )
            log(f"  Hub RDF rows: {len(hub_rdf)}")
            all_hub_properties.append(hub_properties)
            all_hub_rdf.append(hub_rdf)
            advance_progress()

            if make_qc:
                set_progress(f"{pair.fov}: QC overlay")
                qc_paths.append(
                    make_hub_qc_overlay(
                        pair,
                        object_image,
                        intensity_image,
                        nuclei,
                        hub_properties,
                        hub_rdf,
                        out_dir,
                        max_points_draw=int(cfg.get("qc", {}).get("max_points_draw", 5000)),
                        contrast_lower_percentile=float(cfg.get("qc", {}).get("contrast_lower_percentile", 1)),
                        contrast_upper_percentile=float(cfg.get("qc", {}).get("contrast_upper_percentile", 99.8)),
                        nucleus_crop_padding_px=int(cfg.get("qc", {}).get("nucleus_crop_padding_px", 25)),
                        hub_circle_size=float(cfg.get("qc", {}).get("hub_circle_size", 36)),
                        pixel_size_nm=pixel_size_nm,
                        scale_bar_um=float(cfg.get("qc", {}).get("scale_bar_um", 1.0)),
                        rdf_plot_column=rdf_plot_column,
                    )
                )
                advance_progress()

        set_progress("Writing outputs")
        hub_properties_df = pd.concat(all_hub_properties, ignore_index=True) if all_hub_properties else pd.DataFrame()
        if "keep_hub" not in hub_properties_df.columns:
            hub_properties_df["keep_hub"] = pd.Series(dtype=bool)
        hub_rdf_df = pd.concat(all_hub_rdf, ignore_index=True) if all_hub_rdf else pd.DataFrame()
        rdf_aggregate_df = aggregate_dual_channel_rdf(hub_rdf_df, h3k27ac_column=rdf_plot_column, spen_column="spen_rdf_local_norm")
        nucleus_qc_df = pd.concat(all_nucleus_qc, ignore_index=True) if all_nucleus_qc else pd.DataFrame()
        nucleus_qc_df.to_csv(out_dir / "nucleus_qc.csv", index=False)
        hub_properties_df.to_csv(out_dir / "hub_properties.csv", index=False)
        hub_rdf_df.to_csv(out_dir / "hub_rdf_results.csv", index=False)
        rdf_aggregate_df.to_csv(out_dir / "aggregated_rdf_summary.csv", index=False)
        for legacy_name in ["rdf_results.csv", "rdf_summary.csv", "rdf_aggregate_summary.csv"]:
            legacy_path = out_dir / legacy_name
            if legacy_path.exists():
                legacy_path.unlink()
        aggregate_plot_path = make_hub_aggregate_rdf_plot(rdf_aggregate_df, out_dir, hub_rdf_df, rdf_plot_column)
        advance_progress("Finished")
    finally:
        if progress is not None:
            progress.stop()

    return SimpleNamespace(
        config=cfg,
        pairs=pair_df,
        nucleus_qc=nucleus_qc_df,
        hub_properties=hub_properties_df,
        hub_rdf=hub_rdf_df,
        rdf_aggregate=rdf_aggregate_df,
        aggregate_plot_path=aggregate_plot_path,
        qc_paths=qc_paths,
    )


def format_std_multiplier(multiplier: float) -> str:
    value = float(multiplier)
    if value.is_integer():
        return str(int(value))
    return f"{value:g}".replace(".", "p")


def hub_filter_config_for_std_multiplier(base_filter: dict, multiplier: float) -> dict:
    cfg = dict(base_filter or {})
    cfg.update(
        {
            "enabled": True,
            "metric": cfg.get("metric", "spotiflow_intensity"),
            "threshold_source": "nucleus_spen_median_plus_std",
            "std_multiplier": float(multiplier),
            "rule": "greater_equal",
        }
    )
    return cfg


def make_filter_comparison_rdf_overlay(records: list[dict], out_dir: Path) -> Path | None:
    ensure_matplotlib()
    if not records:
        return None
    out_path = out_dir / "filter_comparison_rdf_overlay.png"
    fig, ax = plt.subplots(figsize=(6, 4.2), constrained_layout=True)
    plotted = False
    color_idx = 0
    for record in records:
        aggregate = record.get("rdf_aggregate")
        if aggregate is None or aggregate.empty:
            continue
        x = aggregate["radius_start_nm"].to_numpy(float)
        y = aggregate["h3k27ac_rdf_mean"].to_numpy(float)
        std = aggregate["h3k27ac_rdf_std"].to_numpy(float)
        label = f"{record['filter_label']} ({record['retained_hubs']} hubs)"
        color = f"C{color_idx}"
        plot_rdf_mean_with_std(ax, x, y, std, color, label, zorder=4)
        color_idx += 1
        plotted = True
    if not plotted:
        plt.close(fig)
        return None
    ax.axhline(1.0, color="0.45", lw=1, ls="--")
    ax.set_xlim(float(aggregate["radius_start_nm"].min()), float(aggregate["radius_start_nm"].max()))
    ax.set_xlabel("Distance from SPEN hub center (nm)")
    ax.set_ylabel("local-intensity normalized H3K27ac RDF")
    ax.set_title("Hub-filter comparison RDF")
    clamp_rdf_ylim(ax)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def run_hub_filter_comparison(config_path: str | Path = "config.yaml") -> SimpleNamespace:
    base_cfg = load_config(config_path)
    compare_cfg = base_cfg.get("hub_filter_comparison", {})
    output_root = Path(compare_cfg.get("output_root", Path(base_cfg.get("output_dir", "results")) / "hub_filter_comparison"))
    output_root.mkdir(parents=True, exist_ok=True)
    multipliers = [float(v) for v in _as_list(compare_cfg.get("std_multipliers", [1, 2, 3]))]
    base_output_dir = Path(base_cfg.get("output_dir", "results"))

    summary_rows = []
    records = []
    for multiplier in multipliers:
        mult_label = format_std_multiplier(multiplier)
        filter_label = f"median_plus_{mult_label}std"
        work_dir = output_root / f"_working_{filter_label}"
        if work_dir.exists():
            shutil.rmtree(work_dir)

        run_cfg = copy.deepcopy(base_cfg)
        run_cfg["reuse_output_dir"] = str(base_output_dir)
        run_cfg["output_dir"] = str(work_dir)
        run_cfg.setdefault("cellpose", {})["use_existing_masks"] = True
        run_cfg.setdefault("spotiflow", {})["use_existing_spots"] = True
        run_cfg["hub_filter"] = hub_filter_config_for_std_multiplier(base_cfg.get("hub_filter", {}), multiplier)

        print(f"Running hub filter comparison: {filter_label}", flush=True)
        result = run_pipeline(config_path, cfg=run_cfg)
        retained_hubs = int(result.hub_properties["keep_hub"].sum()) if not result.hub_properties.empty else 0
        final_dir = output_root / f"{filter_label}_{retained_hubs}_hubs"
        if final_dir.exists():
            shutil.rmtree(final_dir)
        work_dir.rename(final_dir)

        row = {
            "filter_label": filter_label,
            "std_multiplier": multiplier,
            "output_dir": str(final_dir),
            "n_fovs": int(len(result.pairs)),
            "nuclear_detections": int(len(result.hub_properties)),
            "retained_hubs": retained_hubs,
            "hub_rdf_rows": int(len(result.hub_rdf)),
        }
        if not result.rdf_aggregate.empty:
            row["first_bin_h3k27ac_rdf_mean"] = float(result.rdf_aggregate["h3k27ac_rdf_mean"].iloc[0])
            row["last_bin_h3k27ac_rdf_mean"] = float(result.rdf_aggregate["h3k27ac_rdf_mean"].iloc[-1])
            row["first_bin_spen_rdf_mean"] = float(result.rdf_aggregate["spen_rdf_mean"].iloc[0])
            row["last_bin_spen_rdf_mean"] = float(result.rdf_aggregate["spen_rdf_mean"].iloc[-1])
        summary_rows.append(row)
        records.append({**row, "rdf_aggregate": result.rdf_aggregate.copy()})

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(output_root / "filter_comparison_summary.csv", index=False)
    overlay_path = make_filter_comparison_rdf_overlay(records, output_root)
    return SimpleNamespace(summary=summary_df, output_root=output_root, overlay_path=overlay_path, records=records)


def _as_list(value) -> list:
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    return [value]


def load_or_run_comparison_nuclei(pair: ImagePair, cfg: dict, compare_dir: Path) -> np.ndarray:
    candidates = [
        Path("results") / "nucleus_masks" / f"{pair.fov}_nuclei.tif",
        Path(cfg.get("output_dir", "results")) / "nucleus_masks" / f"{pair.fov}_nuclei.tif",
        compare_dir / "nucleus_masks" / f"{pair.fov}_nuclei.tif",
    ]
    for mask_path in candidates:
        if mask_path.exists():
            return np.asarray(tifffile.imread(mask_path))

    intensity = read_image(pair.intensity_path)
    masks = run_cellpose_nuclei(
        intensity,
        diameter=cfg["cellpose"]["diameter"],
        flow_threshold=cfg["cellpose"]["flow_threshold"],
        cellprob_threshold=cfg["cellpose"]["cellprob_threshold"],
        device=cfg["cellpose"].get("device", "auto"),
        preprocess=cfg["cellpose"].get("preprocess", {}),
        downsample=cfg["cellpose"].get("downsample", 1.0),
    )
    out_path = compare_dir / "nucleus_masks" / f"{pair.fov}_nuclei.tif"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(out_path, masks)
    return masks


def run_spotiflow_detector_sweep(object_image: np.ndarray, cfg: dict, sweep_cfg: dict, out_dir: str | Path) -> pd.DataFrame:
    try:
        from spotiflow.model import Spotiflow
    except ImportError as exc:
        raise RuntimeError("Spotiflow is not installed in the active environment.") from exc

    torch_device = resolve_torch_device(cfg["spotiflow"].get("device", "auto"))
    cache_path = resolve_spotiflow_cache_dir(cfg, out_dir)
    cache_path.mkdir(parents=True, exist_ok=True)
    model = Spotiflow.from_pretrained(
        cfg["spotiflow"].get("pretrained_model", "general"),
        cache_dir=cache_path,
        map_location=str(torch_device),
    )

    frames = []
    probability_thresholds = sorted(float(v) for v in _as_list(sweep_cfg.get("probability_threshold", [0.5])))
    for min_distance in _as_list(sweep_cfg.get("min_distance", [3])):
        base_threshold = min(probability_thresholds)
        points, details = model.predict(
            object_image,
            prob_thresh=base_threshold,
            min_distance=int(min_distance),
            exclude_border=int(cfg["spotiflow"].get("exclude_border", 1)),
            normalizer=cfg["spotiflow"].get("normalizer", "auto"),
            device=str(torch_device),
        )
        base_spots = pd.DataFrame(np.round(points, 4), columns=["y", "x"])
        if hasattr(details, "prob") and details.prob is not None:
            base_spots["detector_score"] = np.round(details.prob, 4)
            base_spots["probability"] = base_spots["detector_score"]
        else:
            base_spots["detector_score"] = np.nan
            base_spots["probability"] = np.nan
        if hasattr(details, "intens") and details.intens is not None:
            base_spots["spotiflow_intensity"] = np.round(details.intens, 4)

        for prob_thresh in probability_thresholds:
            if "probability" in base_spots:
                spots = base_spots[base_spots["probability"] >= prob_thresh].copy()
            else:
                spots = base_spots.copy()
            spots["method"] = "spotiflow"
            spots["parameter_set"] = f"prob={prob_thresh:.2f};min_distance={int(min_distance)}"
            spots["spotiflow_probability_threshold"] = prob_thresh
            spots["spotiflow_min_distance"] = int(min_distance)
            frames.append(spots)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _dog_response_at(response: np.ndarray, y: float, x: float) -> float:
    yy = int(np.clip(round(float(y)), 0, response.shape[0] - 1))
    xx = int(np.clip(round(float(x)), 0, response.shape[1] - 1))
    return float(response[yy, xx])


def run_dog_detector_sweep(object_image: np.ndarray, cfg: dict, sweep_cfg: dict) -> pd.DataFrame:
    qc_cfg = cfg.get("qc", {})
    normalized = normalize01(
        object_image,
        float(qc_cfg.get("contrast_lower_percentile", 1)),
        float(qc_cfg.get("contrast_upper_percentile", 99.8)),
    )
    frames = []
    overlap = float(sweep_cfg.get("overlap", 0.5))
    for min_sigma in _as_list(sweep_cfg.get("min_sigma", [1])):
        for max_sigma in _as_list(sweep_cfg.get("max_sigma", [6])):
            if float(max_sigma) < float(min_sigma):
                continue
            for threshold_rel in _as_list(sweep_cfg.get("threshold_rel", [0.1])):
                blobs = feature.blob_dog(
                    normalized,
                    min_sigma=float(min_sigma),
                    max_sigma=float(max_sigma),
                    threshold=0.0,
                    threshold_rel=float(threshold_rel),
                    overlap=overlap,
                    exclude_border=int(cfg["spotiflow"].get("exclude_border", 1)),
                )
                spots = pd.DataFrame(blobs, columns=["y", "x", "dog_sigma"]) if len(blobs) else pd.DataFrame(columns=["y", "x", "dog_sigma"])
                if not spots.empty:
                    response = ndimage.gaussian_filter(normalized, sigma=float(min_sigma)) - ndimage.gaussian_filter(
                        normalized,
                        sigma=float(max_sigma),
                    )
                    spots["detector_score"] = [
                        _dog_response_at(response, row.y, row.x) for row in spots.itertuples(index=False)
                    ]
                else:
                    spots["detector_score"] = pd.Series(dtype=float)
                spots["method"] = "dog"
                spots["parameter_set"] = (
                    f"min_sigma={float(min_sigma):g};max_sigma={float(max_sigma):g};"
                    f"threshold_rel={float(threshold_rel):g}"
                )
                spots["dog_min_sigma"] = float(min_sigma)
                spots["dog_max_sigma"] = float(max_sigma)
                spots["dog_threshold_rel"] = float(threshold_rel)
                spots["dog_overlap"] = overlap
                frames.append(spots)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def summarize_detector_comparison(detections: pd.DataFrame) -> pd.DataFrame:
    if detections.empty:
        return pd.DataFrame()
    grouped = detections.groupby(["method", "parameter_set"], dropna=False)
    return grouped.agg(
        n_detections=("y", "count"),
        median_peak_intensity=("peak_intensity", "median"),
        q90_peak_intensity=("peak_intensity", lambda x: float(np.quantile(x, 0.9))),
        median_integrated_intensity=("background_corrected_integrated_intensity", "median"),
        q90_integrated_intensity=("background_corrected_integrated_intensity", lambda x: float(np.quantile(x, 0.9))),
        median_area_px=("punctum_area_px", "median"),
        q90_area_px=("punctum_area_px", lambda x: float(np.quantile(x, 0.9))),
        median_equivalent_diameter_nm=("equivalent_diameter_nm", "median"),
        q90_equivalent_diameter_nm=("equivalent_diameter_nm", lambda x: float(np.quantile(x, 0.9))),
    ).reset_index()


def make_detector_comparison_qc(detections: pd.DataFrame, pairs: list[ImagePair], out_dir: Path, cfg: dict) -> list[Path]:
    ensure_matplotlib()
    paths = []
    if detections.empty:
        return paths
    qc_cfg = cfg.get("qc", {})
    lower = float(qc_cfg.get("contrast_lower_percentile", 1))
    upper = float(qc_cfg.get("contrast_upper_percentile", 99.8))

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    summary = detections.groupby(["method", "parameter_set"]).size().reset_index(name="count")
    for method, method_df in summary.groupby("method"):
        axes[0, 0].plot(np.arange(len(method_df)), method_df["count"], marker="o", label=method)
    axes[0, 0].set_title("Detection count per parameter set")
    axes[0, 0].set_xlabel("Parameter set index")
    axes[0, 0].set_ylabel("Nuclear detections")
    axes[0, 0].legend()

    detections.boxplot(column="background_corrected_integrated_intensity", by="method", ax=axes[0, 1])
    axes[0, 1].set_title("Integrated intensity")
    axes[0, 1].set_xlabel("")
    fig.suptitle("")

    detections.boxplot(column="equivalent_diameter_nm", by="method", ax=axes[1, 0])
    axes[1, 0].set_title("Equivalent diameter")
    axes[1, 0].set_xlabel("")

    for method, method_df in detections.groupby("method"):
        axes[1, 1].scatter(
            method_df["punctum_area_px"],
            method_df["peak_intensity"],
            s=8,
            alpha=0.35,
            label=method,
        )
    axes[1, 1].set_xlabel("Punctum area (px)")
    axes[1, 1].set_ylabel("Peak SPEN intensity")
    axes[1, 1].set_title("Peak intensity vs area")
    axes[1, 1].legend()
    out_path = out_dir / "detector_metric_summary.png"
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    paths.append(out_path)

    for pair in pairs:
        fov_df = detections[detections["fov"] == pair.fov]
        if fov_df.empty:
            continue
        image = read_image(pair.object_path)
        display = normalize01(image, lower, upper)
        fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
        axes[0].imshow(display, cmap="magma")
        axes[0].set_title("SPEN SACD")
        for ax, method in zip(axes[1:], ["spotiflow", "dog"]):
            ax.imshow(display, cmap="magma")
            top = select_top_unique_hubs(fov_df[fov_df["method"] == method], top_n=50)
            if not top.empty:
                sizes = np.clip(top["equivalent_diameter_px"].to_numpy(float) * 8, 12, 90)
                ax.scatter(top["x"], top["y"], s=sizes, facecolors="none", edgecolors="cyan", linewidths=0.8)
            ax.set_title(f"Top {method} hubs")
        for ax in axes:
            ax.set_axis_off()
        out_path = out_dir / f"{pair.fov}_top_hub_overlay.png"
        fig.savefig(out_path, dpi=180)
        plt.close(fig)
        paths.append(out_path)
    return paths


def select_top_unique_hubs(detections: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    if detections.empty:
        return detections.copy()
    cols = ["fov", "method", "y", "x"]
    df = detections.copy()
    df["center_y_px"] = np.rint(df["y"].astype(float)).astype(int)
    df["center_x_px"] = np.rint(df["x"].astype(float)).astype(int)
    sort_cols = ["area_integrated_score", "background_corrected_integrated_intensity", "punctum_area_px"]
    df = df.sort_values(sort_cols, ascending=[False, False, False])
    df = df.drop_duplicates(["fov", "method", "center_y_px", "center_x_px"], keep="first")
    keep_cols = [c for c in df.columns if c not in set(cols) or c in detections.columns]
    return df[keep_cols].groupby("method", as_index=False).head(top_n)


def run_detector_comparison(config_path: str | Path = "config.yaml") -> SimpleNamespace:
    cfg = load_config(config_path)
    comp_cfg = cfg.get("detector_comparison", {})
    out_dir = Path(comp_cfg.get("output_dir", "results/detector_comparison"))
    out_dir.mkdir(parents=True, exist_ok=True)
    configure_runtime_caches(cfg, out_dir)
    pixel_size_nm = float(cfg.get("pixel_size_nm", 58.5))
    max_fovs = comp_cfg.get("max_fovs", 1)
    input_dir = comp_cfg.get("input_dir", cfg["input_dir"])
    object_pattern = comp_cfg.get("object_pattern", cfg["object_pattern"])
    intensity_pattern = comp_cfg.get("intensity_pattern", cfg["intensity_pattern"])
    pairs = pair_sacd_files(
        input_dir,
        object_pattern,
        intensity_pattern,
        max_fovs,
    )
    raise_if_no_pairs(pairs, input_dir, object_pattern, intensity_pattern)
    print(f"Detector comparison FOVs: {len(pairs)}", flush=True)

    all_detections = []
    for pair in pairs:
        print(f"Comparing detectors for {pair.fov}", flush=True)
        object_image = read_image(pair.object_path)
        nuclei = load_or_run_comparison_nuclei(pair, cfg, out_dir)
        spotiflow = run_spotiflow_detector_sweep(object_image, cfg, comp_cfg.get("spotiflow", {}), out_dir)
        dog = run_dog_detector_sweep(object_image, cfg, comp_cfg.get("dog", {}))
        raw = pd.concat([spotiflow, dog], ignore_index=True)
        raw.insert(0, "fov", pair.fov)
        featured = add_shared_spot_features(
            raw,
            object_image,
            nuclei,
            pixel_size_nm,
            int(comp_cfg.get("local_window_px", 15)),
            float(comp_cfg.get("component_threshold_mad", 3)),
            domain=comp_cfg.get("domain", "nuclear"),
        )
        all_detections.append(featured)
        print(f"  Nuclear detections: {len(featured)}", flush=True)

    detections = pd.concat(all_detections, ignore_index=True) if all_detections else pd.DataFrame()
    if not detections.empty:
        detections["area_integrated_score"] = (
            detections["punctum_area_px"].astype(float)
            * detections["background_corrected_integrated_intensity"].astype(float)
        )
        detections["area_integrated_rank"] = detections.groupby("method")["area_integrated_score"].rank(
            ascending=False,
            method="first",
        )
    summary = summarize_detector_comparison(detections)
    top_n = int(comp_cfg.get("top_n", 50))
    top_hubs = select_top_unique_hubs(detections, top_n=top_n) if not detections.empty else pd.DataFrame()

    detections.to_csv(out_dir / "detections.csv", index=False)
    summary.to_csv(out_dir / "summary.csv", index=False)
    top_hubs.to_csv(out_dir / "top_hubs.csv", index=False)
    qc_paths = make_detector_comparison_qc(detections, pairs, out_dir, cfg)
    return SimpleNamespace(config=cfg, pairs=pairs, detections=detections, summary=summary, top_hubs=top_hubs, qc_paths=qc_paths)


def add_raw_coordinate_intensity_features(
    spots: pd.DataFrame,
    image: np.ndarray,
    nuclei: np.ndarray,
    local_max_radii: Iterable[int] = (1, 2, 3),
    domain: str = "nuclear",
) -> pd.DataFrame:
    assigned = assign_spots_to_nuclei(spots, nuclei) if domain == "nuclear" else spots.copy()
    if domain != "nuclear" and not assigned.empty:
        yy = np.clip(np.rint(assigned["y"].to_numpy()).astype(int), 0, nuclei.shape[0] - 1)
        xx = np.clip(np.rint(assigned["x"].to_numpy()).astype(int), 0, nuclei.shape[1] - 1)
        assigned["nucleus_id"] = nuclei[yy, xx].astype(int)

    out = assigned.copy()
    if out.empty:
        out["center_y_px"] = pd.Series(dtype=int)
        out["center_x_px"] = pd.Series(dtype=int)
        out["center_pixel_intensity"] = pd.Series(dtype=float)
        return out

    center_y = np.clip(np.rint(out["y"].to_numpy()).astype(int), 0, image.shape[0] - 1)
    center_x = np.clip(np.rint(out["x"].to_numpy()).astype(int), 0, image.shape[1] - 1)
    out["center_y_px"] = center_y
    out["center_x_px"] = center_x
    out["center_pixel_intensity"] = image[center_y, center_x].astype(float)

    for radius in local_max_radii:
        r = int(radius)
        values = []
        for y, x in zip(center_y, center_x):
            crop = image[max(0, y - r) : min(image.shape[0], y + r + 1), max(0, x - r) : min(image.shape[1], x + r + 1)]
            values.append(float(np.max(crop)))
        out[f"local_max_r{r}"] = values
    return out.drop_duplicates(["center_y_px", "center_x_px"], keep="first").reset_index(drop=True)


def compute_nucleus_qc(nuclei: np.ndarray, pixel_size_nm: float, min_area_px: int, min_edge_distance_px: int) -> pd.DataFrame:
    rows = []
    height, width = nuclei.shape
    for nucleus_id in sorted(int(v) for v in np.unique(nuclei) if v > 0):
        yy, xx = np.where(nuclei == nucleus_id)
        area_px = int(yy.size)
        y0, y1, x0, x1 = int(yy.min()), int(yy.max()), int(xx.min()), int(xx.max())
        edge_distance = int(min(y0, height - 1 - y1, x0, width - 1 - x1))
        reasons = []
        if area_px < int(min_area_px):
            reasons.append("small_area")
        if edge_distance < int(min_edge_distance_px):
            reasons.append("near_edge")
        rows.append(
            {
                "nucleus_id": nucleus_id,
                "area_px": area_px,
                "area_um2": area_px * (float(pixel_size_nm) / 1000.0) ** 2,
                "bbox_y0": y0,
                "bbox_y1": y1,
                "bbox_x0": x0,
                "bbox_x1": x1,
                "min_edge_distance_px": edge_distance,
                "keep_nucleus": len(reasons) == 0,
                "exclude_reason": "keep" if not reasons else ";".join(reasons),
            }
        )
    return pd.DataFrame(rows)


def add_spotiflow_selection(detections: pd.DataFrame, nucleus_qc: pd.DataFrame, std_multiplier: float = 1.0) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = detections.copy()
    if "object_intensity" not in out.columns:
        raise ValueError("Spotiflow screen requires Spotiflow `object_intensity`.")
    out["spotiflow_intensity"] = out["object_intensity"].astype(float)
    kept_ids = set(nucleus_qc.loc[nucleus_qc["keep_nucleus"], "nucleus_id"].astype(int))
    out["keep_nucleus"] = out["nucleus_id"].astype(int).isin(kept_ids)

    stats = (
        out[out["keep_nucleus"]]
        .groupby(["fov", "spotiflow_probability_threshold", "nucleus_id"])
        .agg(
            n_candidates=("spotiflow_intensity", "count"),
            median_spotiflow_intensity=("spotiflow_intensity", "median"),
            std_spotiflow_intensity=("spotiflow_intensity", "std"),
        )
        .reset_index()
    )
    stats["std_spotiflow_intensity"] = stats["std_spotiflow_intensity"].fillna(0.0)
    stats["selection_std_multiplier"] = float(std_multiplier)
    stats["spotiflow_intensity_threshold"] = (
        stats["median_spotiflow_intensity"] + float(std_multiplier) * stats["std_spotiflow_intensity"]
    )
    out = out.merge(
        stats[
            [
                "fov",
                "spotiflow_probability_threshold",
                "nucleus_id",
                "median_spotiflow_intensity",
                "std_spotiflow_intensity",
                "selection_std_multiplier",
                "spotiflow_intensity_threshold",
            ]
        ],
        on=["fov", "spotiflow_probability_threshold", "nucleus_id"],
        how="left",
    )
    out["selected_hub"] = out["keep_nucleus"] & (out["spotiflow_intensity"] >= out["spotiflow_intensity_threshold"])
    selected = out[out["selected_hub"]].copy()
    if not selected.empty:
        selected = selected.sort_values(["fov", "spotiflow_probability_threshold", "nucleus_id", "spotiflow_intensity"], ascending=[True, True, True, False])
        selected["rank_in_nucleus"] = selected.groupby(["fov", "spotiflow_probability_threshold", "nucleus_id"]).cumcount() + 1
    return out, selected


def keep_all_spotiflow_candidates(detections: pd.DataFrame, nucleus_qc: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = detections.copy()
    if "object_intensity" not in out.columns:
        raise ValueError("Spotiflow screen requires Spotiflow `object_intensity`.")
    out["spotiflow_intensity"] = out["object_intensity"].astype(float)
    kept_ids = set(nucleus_qc.loc[nucleus_qc["keep_nucleus"], "nucleus_id"].astype(int))
    out["keep_nucleus"] = out["nucleus_id"].astype(int).isin(kept_ids)
    out["median_spotiflow_intensity"] = np.nan
    out["std_spotiflow_intensity"] = np.nan
    out["selection_std_multiplier"] = 0.0
    out["spotiflow_intensity_threshold"] = np.nan
    out["selected_hub"] = out["keep_nucleus"]
    selected = out[out["selected_hub"]].copy()
    if not selected.empty:
        selected = selected.sort_values(
            ["fov", "spotiflow_probability_threshold", "nucleus_id", "spotiflow_intensity"],
            ascending=[True, True, True, False],
        )
        selected["rank_in_nucleus"] = selected.groupby(["fov", "spotiflow_probability_threshold", "nucleus_id"]).cumcount() + 1
    return out, selected


def make_spotiflow_screen_plots(
    detections: pd.DataFrame,
    selected_hubs: pd.DataFrame,
    nucleus_qc: pd.DataFrame,
    pairs: list[ImagePair],
    out_dir: Path,
    cfg: dict,
) -> list[Path]:
    ensure_matplotlib()
    paths = []
    if detections.empty:
        return paths

    screen_cfg = cfg.get("spotiflow_screen", {})
    intensity_metric = screen_cfg.get("intensity_metric", "spotiflow_intensity")
    circle_size = float(screen_cfg.get("circle_size", 85))
    qc_cfg = cfg.get("qc", {})
    lower = float(qc_cfg.get("contrast_lower_percentile", 1))
    upper = float(qc_cfg.get("contrast_upper_percentile", 99.8))
    kept_ids = sorted(int(v) for v in nucleus_qc.loc[nucleus_qc["keep_nucleus"], "nucleus_id"])

    for threshold, df in detections.groupby("spotiflow_probability_threshold"):
        threshold_label = str(threshold).replace(".", "p")
        th_dir = out_dir / f"threshold_{threshold_label}"
        th_dir.mkdir(parents=True, exist_ok=True)

        candidate_df = df[df["keep_nucleus"]].copy()
        selected_df = selected_hubs[selected_hubs["spotiflow_probability_threshold"] == threshold].copy()
        values = candidate_df[intensity_metric].replace([np.inf, -np.inf], np.nan).dropna().to_numpy(float)
        fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
        if values.size:
            hi = np.quantile(values, 0.995)
            lo = values.min()
            bins = np.linspace(lo, hi, 60) if hi > lo else 30
            ax.hist(values[values <= hi], bins=bins, color="#4c78a8", alpha=0.85)
            ax.axvline(np.median(values), color="black", linestyle="--", lw=1.5, label=f"median={np.median(values):.3g}")
            ax.axvline(np.quantile(values, 0.9), color="#f28e2b", linestyle="--", lw=1.5, label=f"q90={np.quantile(values, 0.9):.3g}")
            ax.axvline(values.max(), color="#d62728", linestyle="-", lw=1.5, label=f"max={values.max():.3g}")
        ax.set_xlabel("Spotiflow intensity")
        ax.set_ylabel("Nuclear Spotiflow detections")
        ax.set_title(f"Spotiflow min_distance={int(screen_cfg.get('min_distance', 1))}, threshold={threshold:g}")
        ax.legend(frameon=False)
        hist_path = th_dir / f"spotiflow_threshold_{threshold_label}_spotiflow_intensity_histogram.png"
        fig.savefig(hist_path, dpi=220)
        plt.close(fig)
        paths.append(hist_path)

        fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=True)
        ax.scatter(candidate_df["center_pixel_intensity"], candidate_df[intensity_metric], s=8, alpha=0.35)
        ax.set_xlabel("Raw center-pixel intensity")
        ax.set_ylabel("Spotiflow intensity")
        ax.set_title(f"Intensity comparison, threshold={threshold:g}")
        scatter_path = th_dir / f"spotiflow_threshold_{threshold_label}_spotiflow_vs_center_intensity.png"
        fig.savefig(scatter_path, dpi=220)
        plt.close(fig)
        paths.append(scatter_path)

        hist_dir = th_dir / "per_nucleus_histograms"
        hist_dir.mkdir(exist_ok=True)
        for (fov, nucleus_id), nuc_df in candidate_df.groupby(["fov", "nucleus_id"]):
            vals = nuc_df[intensity_metric].dropna().to_numpy(float)
            if vals.size == 0:
                continue
            median = float(nuc_df["median_spotiflow_intensity"].iloc[0])
            sd = float(nuc_df["std_spotiflow_intensity"].iloc[0])
            fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
            ax.hist(vals, bins=30, color="#4c78a8", alpha=0.85)
            if np.isfinite(median):
                threshold_value = float(nuc_df["spotiflow_intensity_threshold"].iloc[0])
                ax.axvline(median, color="black", linestyle="--", lw=1.2, label=f"median={median:.3g}")
                ax.axvline(threshold_value, color="#d62728", linestyle="-", lw=1.5, label=f"threshold={threshold_value:.3g}")
            ax.set_xlabel("Spotiflow intensity")
            ax.set_ylabel("Candidate spots")
            sd_text = f"SD={sd:.3g}, " if np.isfinite(sd) else ""
            ax.set_title(f"{fov} nucleus {int(nucleus_id)}\n{sd_text}shown={int((nuc_df['selected_hub']).sum())}")
            ax.legend(frameon=False, fontsize=8)
            p = hist_dir / f"{fov}_nucleus_{int(nucleus_id):03d}_spotiflow_intensity_histogram.png"
            fig.savefig(p, dpi=180)
            plt.close(fig)
            paths.append(p)

        for pair in pairs:
            fov_candidates = candidate_df[candidate_df["fov"] == pair.fov]
            fov_selected = selected_df[selected_df["fov"] == pair.fov]
            image = read_image(pair.object_path)
            mask_candidates = [
                Path("results") / "nucleus_masks" / f"{pair.fov}_nuclei.tif",
                out_dir / "nucleus_masks" / f"{pair.fov}_nuclei.tif",
            ]
            for mask_path in mask_candidates:
                if mask_path.exists():
                    nuclei = np.asarray(tifffile.imread(mask_path))
                    break
            else:
                nuclei = load_or_run_comparison_nuclei(pair, cfg, out_dir)

            display = normalize01(image, lower, upper)
            nuc_ids = [n for n in kept_ids if n in set(int(v) for v in np.unique(nuclei) if v > 0)]
            cols = 4
            rows = max(1, int(np.ceil(len(nuc_ids) / cols)))
            fig, axes = plt.subplots(rows, cols, figsize=(cols * 4.8, rows * 4.8), constrained_layout=True)
            axes = np.atleast_1d(axes).ravel()
            for ax, nucleus_id in zip(axes, nuc_ids):
                ys, xs = np.where(nuclei == nucleus_id)
                pad = 35
                y0, y1 = max(0, ys.min() - pad), min(image.shape[0], ys.max() + pad + 1)
                x0, x1 = max(0, xs.min() - pad), min(image.shape[1], xs.max() + pad + 1)
                ax.imshow(display[y0:y1, x0:x1], cmap="magma")
                boundary = segmentation.find_boundaries(nuclei[y0:y1, x0:x1] == nucleus_id, mode="outer")
                ax.contour(boundary, levels=[0.5], colors="cyan", linewidths=0.7)
                nuc_candidates = fov_candidates[fov_candidates["nucleus_id"] == nucleus_id]
                nuc_selected = fov_selected[fov_selected["nucleus_id"] == nucleus_id]
                if not nuc_selected.empty:
                    ax.scatter(
                        nuc_selected["x"] - x0,
                        nuc_selected["y"] - y0,
                        s=circle_size,
                        facecolors="none",
                        edgecolors="lime",
                        linewidths=1.0,
                    )
                    for row in nuc_selected.itertuples(index=False):
                        ax.text(row.x - x0 + 3, row.y - y0 + 3, str(int(row.rank_in_nucleus)), color="white", fontsize=6.5)
                ax.set_title(f"Nucleus {nucleus_id}: spots {len(nuc_selected)}")
                ax.set_axis_off()
            for ax in axes[len(nuc_ids) :]:
                ax.set_axis_off()
            fig.suptitle(
                f"Spotiflow threshold={threshold:g}, min_distance={int(screen_cfg.get('min_distance', 1))}: "
                f"lime=all retained-nucleus Spotiflow detections"
            )
            montage_path = th_dir / f"{pair.fov}_selected_hubs_by_nucleus.png"
            fig.savefig(montage_path, dpi=220)
            plt.close(fig)
            paths.append(montage_path)

            per_nucleus_dir = th_dir / "per_nucleus"
            per_nucleus_dir.mkdir(exist_ok=True)
            for nucleus_id in nuc_ids:
                ys, xs = np.where(nuclei == nucleus_id)
                pad = 45
                y0, y1 = max(0, ys.min() - pad), min(image.shape[0], ys.max() + pad + 1)
                x0, x1 = max(0, xs.min() - pad), min(image.shape[1], xs.max() + pad + 1)
                fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)
                ax.imshow(display[y0:y1, x0:x1], cmap="magma")
                boundary = segmentation.find_boundaries(nuclei[y0:y1, x0:x1] == nucleus_id, mode="outer")
                ax.contour(boundary, levels=[0.5], colors="cyan", linewidths=0.8)
                nuc_candidates = fov_candidates[fov_candidates["nucleus_id"] == nucleus_id]
                nuc_selected = fov_selected[fov_selected["nucleus_id"] == nucleus_id]
                if not nuc_selected.empty:
                    ax.scatter(
                        nuc_selected["x"] - x0,
                        nuc_selected["y"] - y0,
                        s=circle_size,
                        facecolors="none",
                        edgecolors="lime",
                        linewidths=1.2,
                    )
                    for row in nuc_selected.itertuples(index=False):
                        ax.text(row.x - x0 + 3, row.y - y0 + 3, str(int(row.rank_in_nucleus)), color="white", fontsize=8)
                ax.set_title(f"Nucleus {nucleus_id}: spots {len(nuc_selected)}")
                ax.set_axis_off()
                panel_path = per_nucleus_dir / f"nucleus_{nucleus_id:03d}_selected_hubs.png"
                fig.savefig(panel_path, dpi=220)
                plt.close(fig)
                paths.append(panel_path)

        counts = (
            selected_df.groupby(["fov", "nucleus_id"]).agg(n_selected_hubs=("y", "count")).reset_index()
            if not selected_df.empty
            else pd.DataFrame(columns=["fov", "nucleus_id", "n_selected_hubs"])
        )
        counts.to_csv(th_dir / f"spotiflow_threshold_{threshold_label}_selected_hub_counts_by_nucleus.csv", index=False)
    return paths


def run_spotiflow_screen(config_path: str | Path = "config.yaml") -> SimpleNamespace:
    cfg = load_config(config_path)
    screen_cfg = cfg.get("spotiflow_screen", {})
    out_dir = Path(screen_cfg.get("output_dir", "results/detector_comparison/spotiflow_center_intensity_screen"))
    out_dir.mkdir(parents=True, exist_ok=True)
    configure_runtime_caches(cfg, out_dir)
    thresholds = sorted(float(v) for v in _as_list(screen_cfg.get("thresholds", [0.1, 0.2])))
    min_distance = int(screen_cfg.get("min_distance", 1))
    input_dir = screen_cfg.get("input_dir", cfg["input_dir"])
    object_pattern = screen_cfg.get("object_pattern", cfg["object_pattern"])
    intensity_pattern = screen_cfg.get("intensity_pattern", cfg["intensity_pattern"])
    pairs = pair_sacd_files(
        input_dir,
        object_pattern,
        intensity_pattern,
        screen_cfg.get("max_fovs", 1),
    )
    raise_if_no_pairs(pairs, input_dir, object_pattern, intensity_pattern)
    print(f"Spotiflow screen FOVs: {len(pairs)}", flush=True)

    all_detections = []
    all_nucleus_qc = []
    for pair in pairs:
        print(f"Screening {pair.fov}", flush=True)
        image = read_image(pair.object_path)
        nuclei = load_or_run_comparison_nuclei(pair, cfg, out_dir)
        nucleus_qc = compute_nucleus_qc(
            nuclei,
            float(cfg.get("pixel_size_nm", 58.5)),
            int(screen_cfg.get("min_nucleus_area_px", 27000)),
            int(screen_cfg.get("min_edge_distance_px", 20)),
        )
        nucleus_qc.insert(0, "fov", pair.fov)
        all_nucleus_qc.append(nucleus_qc)
        base_spots = run_spotiflow_spots(
            image,
            pretrained_model=cfg["spotiflow"]["pretrained_model"],
            cache_dir=resolve_spotiflow_cache_dir(cfg, out_dir),
            probability_threshold=min(thresholds),
            min_distance=min_distance,
            exclude_border=cfg["spotiflow"].get("exclude_border", 1),
            normalizer=cfg["spotiflow"].get("normalizer", "auto"),
            device=cfg["spotiflow"].get("device", "cpu"),
        )
        base_spots.insert(0, "fov", pair.fov)
        base_spots["method"] = "spotiflow"
        base_spots["spotiflow_min_distance"] = min_distance
        if "probability" in base_spots.columns:
            base_spots["detector_score"] = base_spots["probability"]
        for threshold in thresholds:
            spots = base_spots[base_spots.get("probability", pd.Series(1.0, index=base_spots.index)) >= threshold].copy()
            spots["spotiflow_probability_threshold"] = threshold
            featured = add_raw_coordinate_intensity_features(
                spots,
                image,
                nuclei,
                local_max_radii=screen_cfg.get("local_max_radii", [1, 2, 3]),
                domain="nuclear",
            )
            all_detections.append(featured)
            print(f"  threshold={threshold:g}: nuclear detections={len(featured)}", flush=True)

    detections = pd.concat(all_detections, ignore_index=True) if all_detections else pd.DataFrame()
    nucleus_qc_df = pd.concat(all_nucleus_qc, ignore_index=True) if all_nucleus_qc else pd.DataFrame()
    if screen_cfg.get("selection_rule", "none") == "median_plus_std":
        detections, selected_hubs = add_spotiflow_selection(
            detections,
            nucleus_qc_df,
            float(screen_cfg.get("selection_std_multiplier", 1.0)),
        )
    else:
        detections, selected_hubs = keep_all_spotiflow_candidates(detections, nucleus_qc_df)
    nucleus_qc_df.to_csv(out_dir / "nucleus_qc.csv", index=False)
    detections.to_csv(out_dir / "spotiflow_screen_detections.csv", index=False)
    selected_hubs.to_csv(out_dir / "spotiflow_selected_hubs.csv", index=False)
    summary = (
        detections.groupby(["fov", "spotiflow_probability_threshold", "nucleus_id"])
        .agg(
            n_detections=("y", "count"),
            keep_nucleus=("keep_nucleus", "first"),
            n_selected_hubs=("selected_hub", "sum"),
            median_spotiflow_intensity=("spotiflow_intensity", "median"),
            std_spotiflow_intensity=("spotiflow_intensity", "std"),
            spotiflow_intensity_threshold=("spotiflow_intensity_threshold", "first"),
            q90_spotiflow_intensity=("spotiflow_intensity", lambda x: float(np.quantile(x, 0.9))),
            max_spotiflow_intensity=("spotiflow_intensity", "max"),
        )
        .reset_index()
        if not detections.empty
        else pd.DataFrame()
    )
    summary.to_csv(out_dir / "spotiflow_screen_summary_by_nucleus.csv", index=False)
    qc_paths = make_spotiflow_screen_plots(detections, selected_hubs, nucleus_qc_df, pairs, out_dir, cfg)
    return SimpleNamespace(
        config=cfg,
        pairs=pairs,
        nucleus_qc=nucleus_qc_df,
        detections=detections,
        selected_hubs=selected_hubs,
        summary=summary,
        qc_paths=qc_paths,
    )


def summed_source_reference_check(radius_px: int = 18, seed: int = 7) -> dict:
    rng = np.random.default_rng(seed)
    shape = (96, 96)
    centers = np.array([[24, 28], [38, 65], [68, 33], [72, 72]], dtype=float)
    radii = np.arange(radius_px + 1)
    true_rdf = 80 + 170 * np.exp(-(radii / 5.0) ** 2) + 30 * np.exp(-((radii - 12) / 3.5) ** 2)
    yy, xx = np.indices(shape)
    image = np.zeros(shape, dtype=float)
    for cy, cx in centers:
        d = np.rint(np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)).astype(int)
        keep = d <= radius_px
        image[keep] += true_rdf[d[keep]]
    image += rng.normal(0, 3, shape)
    nucleus = image > 0
    fit, counts = rdf_fit_original(image, nucleus, centers, radius_px)
    rmse = float(np.sqrt(np.nanmean((fit - true_rdf) ** 2)))
    corr = float(np.corrcoef(fit, true_rdf)[0, 1])
    return {"rmse": rmse, "pearson_r": corr, "min_counts": float(np.min(counts))}


def synthetic_rdf_check(seed: int = 7) -> dict:
    rng = np.random.default_rng(seed)
    shape = (96, 96)
    yy, xx = np.indices(shape)
    nucleus = (((yy - 48) ** 2 + (xx - 48) ** 2) <= 42**2).astype(np.uint16)
    centers = np.array([[38.0, 39.0], [59.0, 57.0]], dtype=float)

    object_image = rng.normal(0.02, 0.003, shape).astype(float)
    intensity_image = np.full(shape, 1.0, dtype=float)
    nearest_nm = np.full(shape, np.inf, dtype=float)
    for cy, cx in centers:
        d_px = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)
        nearest_nm = np.minimum(nearest_nm, d_px * 50.0)
        object_image += 0.6 * np.exp(-(d_px / 1.8) ** 2)
    intensity_image += 0.7 * np.exp(-(nearest_nm / 120.0) ** 2)
    intensity_image += rng.normal(0, 0.01, shape)
    intensity_image[nucleus == 0] = 0

    spots = pd.DataFrame(
        {
            "y": centers[:, 0],
            "x": centers[:, 1],
            "probability": [0.9, 0.9],
            "object_intensity": [0.7, 1.0],
            "nucleus_id": [1, 1],
        }
    )
    bins = make_physical_rdf_bins(radius_nm=600, bin_width_nm=200, bin_step_nm=50)
    hub_props, hub_rdf = compute_hub_properties_and_annular_rdf(
        "synthetic",
        object_image,
        intensity_image,
        nucleus,
        spots,
        [1],
        bins,
        pixel_size_nm=50.0,
        hub_filter_cfg={
            "enabled": True,
            "metric": "spotiflow_intensity",
            "threshold_source": "nucleus_spen_median_plus_std",
            "std_multiplier": 2.0,
            "rule": "greater_equal",
            "hard_threshold": {
                "enabled": True,
                "source": "manual",
                "value": 0.8,
                "metric": "spotiflow_intensity",
            },
        },
        rdf_cfg={
            "normalization": "local_intensity_mean",
            "tail_normalization": {"enabled": False, "last_n_bins": 5},
            "aggregation": {"plot_column": "h3k27ac_rdf_local_norm"},
        },
    )
    agg = aggregate_hub_rdf(hub_rdf)
    first = float(agg.loc[agg["bin_index"].idxmin(), "h3k27ac_rdf_mean"])
    last = float(agg.loc[agg["bin_index"].idxmax(), "h3k27ac_rdf_mean"])
    return {
        "mode": "per_hub_annular",
        "retained_hubs": int(hub_props["keep_hub"].sum()),
        "hard_threshold": float(hub_props["hub_filter_hard_threshold"].dropna().iloc[0]),
        "rdf_rows": int(len(hub_rdf)),
        "bins": int(len(agg)),
        "first_bin_rdf": first,
        "last_bin_rdf": last,
        "first_gt_last": bool(first > last),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Bulk fluorescence SPEN-H3K27ac per-hub RDF pipeline")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Run the full pipeline")
    run.add_argument("--config", default="config.yaml")
    sub.add_parser("check", help="Run the synthetic RDF check")
    sub.add_parser("pair", help="Pair inputs and write paired_inputs.csv").add_argument("--config", default="config.yaml")
    compare = sub.add_parser("compare-detectors", help="Compare Spotiflow and DoG SPEN hub detection")
    compare.add_argument("--config", default="config.yaml")
    hub_compare = sub.add_parser("compare-hub-filters", help="Compare median plus STD hub filters")
    hub_compare.add_argument("--config", default="config.yaml")
    screen = sub.add_parser("screen-spotiflow", help="Screen Spotiflow thresholds using raw SPEN center-pixel intensity")
    screen.add_argument("--config", default="config.yaml")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "check":
        print(synthetic_rdf_check())
    elif args.command == "pair":
        cfg = load_config(args.config)
        out_dir = Path(cfg["output_dir"])
        pairs = pair_sacd_files(cfg["input_dir"], cfg["object_pattern"], cfg["intensity_pattern"], cfg.get("max_fovs"))
        raise_if_no_pairs(pairs, cfg["input_dir"], cfg["object_pattern"], cfg["intensity_pattern"])
        df = write_pairs_csv(pairs, out_dir)
        print(df.to_string(index=False))
    elif args.command == "run":
        result = run_pipeline(args.config)
        print(f"FOVs: {len(result.pairs)}")
        print(f"Nuclear SPEN detections: {len(result.hub_properties)}")
        print(f"Retained SPEN hubs: {int(result.hub_properties['keep_hub'].sum()) if not result.hub_properties.empty else 0}")
        print(f"Hub RDF rows: {len(result.hub_rdf)}")
    elif args.command == "compare-detectors":
        result = run_detector_comparison(args.config)
        print(f"FOVs: {len(result.pairs)}")
        print(f"Detections: {len(result.detections)}")
        print(f"Parameter summaries: {len(result.summary)}")
    elif args.command == "compare-hub-filters":
        result = run_hub_filter_comparison(args.config)
        print(f"Output root: {result.output_root}")
        print(result.summary.to_string(index=False))
        print(f"Comparison RDF overlay: {result.overlay_path}")
    elif args.command == "screen-spotiflow":
        result = run_spotiflow_screen(args.config)
        print(f"FOVs: {len(result.pairs)}")
        print(f"Detections: {len(result.detections)}")
        print(f"Selected hubs: {len(result.selected_hubs)}")
        print(f"QC files: {len(result.qc_paths)}")


if __name__ == "__main__":
    main()
