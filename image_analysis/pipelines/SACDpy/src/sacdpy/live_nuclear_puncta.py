"""Frame-2 nucleus segmentation and PunctaTools analysis for SACD live MIPs."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import shutil
import sys
import time
import uuid
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import tifffile
from scipy import ndimage, stats

from .phase_diagram import preprocess_cellpose_image
from .puncta import PunctaTools2DParams, punctatools_segment_2d


DEFAULT_DATASET_ROOT = Path(
    "/Volumes/guttman/users/gmgao/Imaging_ProcessedData/SPEN/SACDlive/"
    "20260617_ONI-gmgao-SPEN_SACDlive-SHA_bef_FVP_rec"
)
CONDITIONS = ("before", "FVP_2h", "recover_2h")
CONDITION_COUNTS = {"before": 9, "FVP_2h": 6, "recover_2h": 6}
DEFAULT_PUNCTA_PARAMS = PunctaTools2DParams()


@dataclass(frozen=True)
class LiveNuclearPunctaConfig:
    dataset_root: Path = DEFAULT_DATASET_ROOT
    output_folder: str = "selected_nucleus"
    frame_index: int = 1
    pixel_size_um: float = 0.0585
    cellpose_model: str = "cpsam"
    cellpose_diameter_px: float = 140.0
    cellpose_flow_threshold: float = 0.4
    cellpose_cellprob_threshold: float = 0.0
    cellpose_lower_percentile: float = 1.0
    cellpose_upper_percentile: float = 99.8
    cellpose_device: str = "cpu"
    crop_padding_px: int = 5
    nucleus_area_quantile: float = 10.0
    pilot_seed: int = 20260722
    pilot_per_condition_per_stratum: int = 2
    display_lower_percentile: float = 5.0
    display_upper_percentile: float = 95.0
    png_dpi: int = 1200
    png_max_edge_inches: float = 3.0
    scale_bar_um: float = 2.0
    bootstrap_iterations: int = 10_000
    fixed_threshold_enabled: bool = False
    fixed_threshold_multipliers: tuple[float, ...] = (1.25, 1.5, 2.0, 3.0)
    punctatools: PunctaTools2DParams = field(default_factory=PunctaTools2DParams)

    @property
    def output_root(self) -> Path:
        return self.dataset_root / self.output_folder

    @property
    def staging_root(self) -> Path:
        return self.dataset_root / f".{self.output_folder}.staging"


@dataclass(frozen=True)
class MIPPlan:
    path: Path
    condition: str
    fov: str


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _progress_eta(label: str, index: int, total: int, started: float) -> None:
    elapsed = max(time.perf_counter() - started, 0.0)
    remaining = elapsed / index * (total - index) if index else float("nan")
    print(
        f"{label} {index}/{total} elapsed={elapsed:.1f}s eta={remaining:.1f}s",
        flush=True,
    )


def parse_condition(name: str) -> str:
    if "recover2h" in name:
        return "recover_2h"
    if "SHA-FVP2h-" in name:
        return "FVP_2h"
    if "SHA-before-" in name:
        return "before"
    raise ValueError(f"Cannot assign condition from filename: {name}")


def fov_from_mip(path: str | Path) -> str:
    name = Path(path).name
    suffix = "-SACDpy-647-posXY0-MIP-TYX.tif"
    if not name.endswith(suffix):
        raise ValueError(f"Not a SACD MIP filename: {name}")
    return name[: -len(suffix)]


def discover_mips(dataset_root: str | Path) -> list[MIPPlan]:
    root = Path(dataset_root)
    if not root.is_dir():
        raise FileNotFoundError(root)
    paths = sorted(root.glob("*-SACDpy-647-posXY0-MIP-TYX.tif"))
    plans = [MIPPlan(path=p, condition=parse_condition(p.name), fov=fov_from_mip(p)) for p in paths]
    counts = {condition: sum(plan.condition == condition for plan in plans) for condition in CONDITIONS}
    if counts != CONDITION_COUNTS:
        raise ValueError(f"Expected condition counts {CONDITION_COUNTS}, found {counts}")
    if len({plan.fov for plan in plans}) != len(plans):
        raise ValueError("Duplicate MIP-derived FOV identifiers")
    return plans


def bbox_with_padding(
    mask: np.ndarray,
    padding: int,
    shape: tuple[int, int] | None = None,
) -> tuple[int, int, int, int]:
    binary = np.asarray(mask, dtype=bool)
    if binary.ndim != 2 or not np.any(binary):
        raise ValueError("Bounding boxes require a nonempty 2D mask")
    height, width = shape or binary.shape
    ys, xs = np.where(binary)
    return (
        max(0, int(ys.min()) - padding),
        min(height, int(ys.max()) + 1 + padding),
        max(0, int(xs.min()) - padding),
        min(width, int(xs.max()) + 1 + padding),
    )


def touches_fov_boundary(mask: np.ndarray) -> bool:
    binary = np.asarray(mask, dtype=bool)
    return bool(
        np.any(binary[0])
        or np.any(binary[-1])
        or np.any(binary[:, 0])
        or np.any(binary[:, -1])
    )


def pooled_area_cutoff(areas: Sequence[int], quantile: float = 10.0) -> float:
    values = np.asarray(areas, dtype=np.float64)
    if not len(values):
        raise ValueError("Cannot calculate an area cutoff without nuclei")
    return float(np.percentile(values, quantile, method="linear"))


def default_puncta_params() -> PunctaTools2DParams:
    """Return a fresh immutable copy of the PunctaTools-compatible defaults."""

    return PunctaTools2DParams(**asdict(DEFAULT_PUNCTA_PARAMS))


def fixed_threshold_params(multiplier: float) -> PunctaTools2DParams:
    if multiplier <= 0:
        raise ValueError("Fixed-threshold multiplier must be positive")
    return PunctaTools2DParams(
        minsize_um=DEFAULT_PUNCTA_PARAMS.minsize_um,
        maxsize_um=DEFAULT_PUNCTA_PARAMS.maxsize_um,
        num_sigma=DEFAULT_PUNCTA_PARAMS.num_sigma,
        overlap=DEFAULT_PUNCTA_PARAMS.overlap,
        threshold_detection=0.0,
        threshold_background=float(multiplier),
        threshold_segmentation=float(multiplier),
        segmentation_mode=2,
        maxrad_um=DEFAULT_PUNCTA_PARAMS.maxrad_um,
    )


def puncta_pilot_safety_metrics(
    puncta_counts: Sequence[float],
    puncta_area_fractions: Sequence[float],
) -> dict[str, Any]:
    counts = np.asarray(puncta_counts, dtype=float)
    area_fractions = np.asarray(puncta_area_fractions, dtype=float)
    if not len(counts) or len(counts) != len(area_fractions):
        raise ValueError("Pilot counts and area fractions must be nonempty and aligned")
    result: dict[str, Any] = {
        "median_puncta_area_fraction": float(np.median(area_fractions)),
        "maximum_puncta_area_fraction": float(np.max(area_fractions)),
        "median_puncta_count": float(np.median(counts)),
    }
    result["passes_safety_gates"] = bool(
        np.all(np.isfinite(area_fractions))
        and np.all(np.isfinite(counts))
        and result["median_puncta_area_fraction"] <= 0.35
        and result["maximum_puncta_area_fraction"] <= 0.80
        and result["median_puncta_count"] < 200
    )
    return result


def pooled_nucleus_background(rows: Sequence[Mapping[str, Any]]) -> float:
    values = np.asarray(
        [float(row["nucleus_median_sacd"]) for row in rows],
        dtype=float,
    )
    if not len(values) or not np.all(np.isfinite(values)):
        raise ValueError("Full-nucleus background values must be finite and nonempty")
    return float(np.median(values))


def run_cellpose_frame(image: np.ndarray, model: Any, config: LiveNuclearPunctaConfig) -> np.ndarray:
    normalized = preprocess_cellpose_image(
        image,
        config.cellpose_lower_percentile,
        config.cellpose_upper_percentile,
    )
    model_input = np.zeros((2,) + normalized.shape, dtype=np.float32)
    model_input[0] = normalized
    result = model.eval(
        model_input,
        batch_size=8,
        diameter=config.cellpose_diameter_px,
        flow_threshold=config.cellpose_flow_threshold,
        cellprob_threshold=config.cellpose_cellprob_threshold,
        normalize={"tile_norm_blocksize": 0},
    )
    labels = np.asarray(result[0] if isinstance(result, tuple) else result)
    if labels.shape != image.shape:
        raise ValueError(f"Cellpose returned {labels.shape}, expected {image.shape}")
    maximum = int(np.max(labels)) if labels.size else 0
    dtype = np.uint16 if maximum <= np.iinfo(np.uint16).max else np.uint32
    return labels.astype(dtype, copy=False)


def _write_imagej(
    path: Path,
    image: np.ndarray,
    *,
    axes: str,
    pixel_size_um: float,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    array = np.asarray(image)
    # ImageJ's legacy TIFF format cannot encode uint32. Instance-label images
    # may legitimately use that dtype, so keep their IDs exact in a standard
    # TIFF with tifffile's axes metadata. Float crops and uint16 labels retain
    # ImageJ compatibility.
    imagej_compatible = array.dtype != np.dtype(np.uint32)
    tifffile.imwrite(
        temporary,
        array,
        imagej=imagej_compatible,
        metadata={"axes": axes, "unit": "um"},
        resolution=(1.0 / pixel_size_um, 1.0 / pixel_size_um),
        photometric="minisblack",
    )
    os.replace(temporary, path)


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    os.replace(temporary, path)


def _csv_value(value: Any) -> Any:
    if isinstance(value, (list, tuple, dict)):
        return json.dumps(value, sort_keys=True)
    return value


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row}) if rows else []
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        if fields:
            writer.writeheader()
            writer.writerows(
                {field: _csv_value(row.get(field, "")) for field in fields}
                for row in rows
            )
    os.replace(temporary, path)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _safe_stem(value: str) -> str:
    return "".join(character if character.isalnum() or character in "-_" else "_" for character in value)


def nucleus_key(condition: str, fov: str, nucleus_id: int) -> str:
    return f"{condition}__{_safe_stem(fov)}__nucleus-{nucleus_id:04d}"


def crop_filename(condition: str, fov: str, nucleus_id: int) -> str:
    return f"{nucleus_key(condition, fov, nucleus_id)}__SACD-mask-CYX.tif"


def _condition_counts(plans: Sequence[MIPPlan]) -> dict[str, int]:
    return {condition: sum(plan.condition == condition for plan in plans) for condition in CONDITIONS}


def discovery_summary(config: LiveNuclearPunctaConfig) -> dict[str, Any]:
    plans = discover_mips(config.dataset_root)
    shapes: dict[str, list[int]] = {}
    for plan in plans:
        with tifffile.TiffFile(plan.path) as tif:
            series = tif.series[0]
            if series.axes != "TYX":
                raise ValueError(f"{plan.path.name} has axes {series.axes}, expected TYX")
            if series.shape[0] <= config.frame_index:
                raise ValueError(f"{plan.path.name} lacks frame {config.frame_index}")
            if str(series.dtype) != "float32":
                raise ValueError(f"{plan.path.name} has dtype {series.dtype}, expected float32")
            shapes[plan.path.name] = list(series.shape)
    return {
        "dataset_root": str(config.dataset_root),
        "output_root": str(config.output_root),
        "frame_index": config.frame_index,
        "mip_count": len(plans),
        "condition_counts": _condition_counts(plans),
        "shapes": shapes,
    }


def _load_cellpose(config: LiveNuclearPunctaConfig) -> Any:
    from cellpose import core as cellpose_core
    from cellpose import models

    gpu = config.cellpose_device == "auto" and bool(cellpose_core.use_gpu())
    print(f"Loading CellposeSAM model={config.cellpose_model} gpu={gpu}", flush=True)
    return models.CellposeModel(gpu=gpu, pretrained_model=config.cellpose_model)


def _read_frame(path: Path, frame_index: int) -> np.ndarray:
    image = tifffile.imread(path, key=frame_index)
    image = np.asarray(image, dtype=np.float32)
    if image.ndim != 2 or not np.all(np.isfinite(image)):
        raise ValueError(f"Invalid frame from {path}: {image.shape}")
    return image


def _save_segmentation_record(
    plan: MIPPlan,
    labels: np.ndarray,
    config: LiveNuclearPunctaConfig,
    analysis_root: Path,
) -> tuple[list[dict[str, Any]], Path]:
    label_path = _segmentation_label_path(plan, config, analysis_root)
    _write_imagej(label_path, labels, axes="YX", pixel_size_um=config.pixel_size_um)
    rows: list[dict[str, Any]] = []
    for label_id in (int(value) for value in np.unique(labels) if value > 0):
        mask = labels == label_id
        ys, xs = np.where(mask)
        raw_bbox = bbox_with_padding(mask, 0)
        crop_bbox = bbox_with_padding(mask, config.crop_padding_px)
        rows.append(
            {
                "nucleus_key": nucleus_key(plan.condition, plan.fov, label_id),
                "condition": plan.condition,
                "fov": plan.fov,
                "source_mip": str(plan.path),
                "frame_index_zero_based": config.frame_index,
                "frame_number_one_based": config.frame_index + 1,
                "nucleus_id": label_id,
                "area_px": int(mask.sum()),
                "centroid_y_px": float(np.mean(ys)),
                "centroid_x_px": float(np.mean(xs)),
                "bbox_y0": raw_bbox[0],
                "bbox_y1": raw_bbox[1],
                "bbox_x0": raw_bbox[2],
                "bbox_x1": raw_bbox[3],
                "crop_y0": crop_bbox[0],
                "crop_y1": crop_bbox[1],
                "crop_x0": crop_bbox[2],
                "crop_x1": crop_bbox[3],
                "touches_fov_boundary": touches_fov_boundary(mask),
                "segmentation_label_path": str(label_path),
            }
        )
    return rows, label_path


def _segmentation_label_path(
    plan: MIPPlan,
    config: LiveNuclearPunctaConfig,
    analysis_root: Path,
) -> Path:
    return (
        analysis_root
        / "segmentation_labels"
        / plan.condition
        / f"{_safe_stem(plan.fov)}__frame-{config.frame_index + 1:02d}__CellposeSAM-labels-YX.tif"
    )


def _retained_and_excluded(
    rows: Sequence[Mapping[str, Any]],
    area_cutoff: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    retained: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    for source in rows:
        row = dict(source)
        small = int(row["area_px"]) < area_cutoff
        border = bool(row["touches_fov_boundary"])
        reasons: list[str] = []
        if small:
            reasons.append("bottom_10_percent_area")
        if border:
            reasons.append("touches_fov_boundary")
        row["area_cutoff_px"] = area_cutoff
        row["small_nucleus"] = small
        row["retained"] = not reasons
        row["exclusion_reasons"] = reasons
        (excluded if reasons else retained).append(row)
    return retained, excluded


def _export_crop(
    row: dict[str, Any],
    config: LiveNuclearPunctaConfig,
    staging_root: Path,
) -> tuple[dict[str, Any], np.ndarray]:
    image = _read_frame(Path(row["source_mip"]), config.frame_index)
    labels = tifffile.imread(row["segmentation_label_path"])
    label_id = int(row["nucleus_id"])
    mask = labels == label_id
    y0, y1 = int(row["crop_y0"]), int(row["crop_y1"])
    x0, x1 = int(row["crop_x0"]), int(row["crop_x1"])
    image_crop = image[y0:y1, x0:x1]
    mask_crop = mask[y0:y1, x0:x1]
    if not np.any(mask_crop):
        raise ValueError(f"Empty nucleus mask for {row['nucleus_key']}")
    stack = np.stack((image_crop, mask_crop.astype(np.float32)), axis=0).astype(np.float32)
    crop_path = staging_root / row["condition"] / crop_filename(
        row["condition"],
        row["fov"],
        label_id,
    )
    _write_imagej(crop_path, stack, axes="CYX", pixel_size_um=config.pixel_size_um)
    result = dict(row)
    result.update(
        {
            "crop_tif": str(crop_path),
            "crop_shape_y": int(image_crop.shape[0]),
            "crop_shape_x": int(image_crop.shape[1]),
            "nucleus_area_px": int(mask_crop.sum()),
            "nucleus_mean_sacd": float(np.mean(image_crop[mask_crop])),
            "nucleus_median_sacd": float(np.median(image_crop[mask_crop])),
        }
    )
    return result, np.asarray(image_crop[mask_crop], dtype=np.float32)


def select_blinded_pilot(
    rows: Sequence[Mapping[str, Any]],
    *,
    seed: int,
    per_condition_per_stratum: int,
) -> list[dict[str, Any]]:
    rng = np.random.default_rng(seed)
    selected: list[dict[str, Any]] = []
    for condition in CONDITIONS:
        candidates = sorted(
            (dict(row) for row in rows if row["condition"] == condition),
            key=lambda row: (float(row["nucleus_mean_sacd"]), row["nucleus_key"]),
        )
        if len(candidates) < 3 * per_condition_per_stratum:
            raise ValueError(f"Not enough retained nuclei for pilot condition {condition}")
        strata = np.array_split(np.asarray(candidates, dtype=object), 3)
        for stratum_index, stratum in enumerate(strata):
            choices = rng.choice(len(stratum), size=per_condition_per_stratum, replace=False)
            for choice in sorted(int(value) for value in choices):
                row = dict(stratum[choice])
                row["pilot_stratum"] = ("low", "medium", "high")[stratum_index]
                selected.append(row)
    rng.shuffle(selected)
    for index, row in enumerate(selected, start=1):
        row["blinded_pilot_id"] = f"pilot-{index:03d}"
    return selected


def quantify_puncta(
    image: np.ndarray,
    nucleus_mask: np.ndarray,
    labels: np.ndarray,
    *,
    pixel_size_um: float,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    values_image = np.asarray(image, dtype=np.float64)
    binary_nucleus = np.asarray(nucleus_mask, dtype=bool)
    instances = np.asarray(labels, dtype=np.uint32)
    if values_image.shape != binary_nucleus.shape or values_image.shape != instances.shape:
        raise ValueError("Puncta image/nucleus/labels are not aligned")
    if not np.any(binary_nucleus):
        raise ValueError("Nucleus mask is empty")
    if np.any(instances[~binary_nucleus] > 0):
        raise ValueError("Puncta labels extend outside the Cellpose nucleus")
    puncta_mask = instances > 0
    background = binary_nucleus & ~puncta_mask
    nonpuncta_mean = (
        float(np.mean(values_image[background]))
        if np.any(background)
        else float("nan")
    )
    nucleus_boundary = binary_nucleus & ndimage.binary_dilation(~binary_nucleus)
    puncta_rows: list[dict[str, Any]] = []
    for punctum_id in (int(value) for value in np.unique(instances) if value > 0):
        punctum = instances == punctum_id
        values = values_image[punctum]
        corrected = (
            np.maximum(values - nonpuncta_mean, 0.0)
            if np.isfinite(nonpuncta_mean)
            else np.full_like(values, np.nan)
        )
        ys, xs = np.where(punctum)
        puncta_rows.append(
            {
                "punctum_id": punctum_id,
                "area_px": int(values.size),
                "area_um2": float(values.size * pixel_size_um**2),
                "raw_mean_sacd": float(np.mean(values)),
                "corrected_mean_sacd": float(np.mean(corrected)),
                "max_sacd": float(np.max(values)),
                "raw_integrated_sacd": float(np.sum(values)),
                "corrected_integrated_sacd": float(np.sum(corrected)),
                "centroid_y_crop_px": float(np.mean(ys)),
                "centroid_x_crop_px": float(np.mean(xs)),
                "touches_nucleus_boundary": bool(np.any(punctum & nucleus_boundary)),
            }
        )
    count = len(puncta_rows)
    corrected_means = [float(row["corrected_mean_sacd"]) for row in puncta_rows]
    areas = [float(row["area_um2"]) for row in puncta_rows]
    summary = {
        "puncta_count": count,
        "mean_corrected_puncta_brightness": (
            float(np.mean(corrected_means)) if corrected_means else float("nan")
        ),
        "median_punctum_area_um2": float(np.median(areas)) if areas else float("nan"),
        "puncta_area_fraction": float(np.sum(puncta_mask) / np.sum(binary_nucleus)),
        "total_corrected_puncta_intensity": float(
            np.nansum([row["corrected_integrated_sacd"] for row in puncta_rows])
        ),
        "nonpuncta_nucleus_mean_sacd": (
            nonpuncta_mean if count else float(np.mean(values_image[binary_nucleus]))
        ),
    }
    return summary, puncta_rows


def _run_default_puncta(
    image: np.ndarray,
    mask: np.ndarray,
    config: LiveNuclearPunctaConfig,
) -> np.ndarray:
    return punctatools_segment_2d(
        image,
        mask,
        pixel_size_um=config.pixel_size_um,
        params=default_puncta_params(),
    )


def _make_overlay(
    path: Path,
    image: np.ndarray,
    nucleus_mask: np.ndarray,
    puncta_labels: np.ndarray,
    *,
    title: str,
    display_limits: tuple[float, float],
    config: LiveNuclearPunctaConfig,
    dpi: int | None = None,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.lines import Line2D

    height, width = image.shape
    maximum = max(height, width)
    figure_size = (
        config.png_max_edge_inches * width / maximum,
        config.png_max_edge_inches * height / maximum,
    )
    render_dpi = dpi or config.png_dpi
    fig, ax = plt.subplots(figsize=figure_size, dpi=render_dpi)
    lo, hi = display_limits
    ax.imshow(
        np.clip(image, lo, hi),
        cmap="magma",
        norm=LogNorm(vmin=lo, vmax=hi),
        interpolation="nearest",
    )
    ax.contour(nucleus_mask.astype(np.uint8), [0.5], colors=["cyan"], linewidths=0.8)
    if np.any(puncta_labels):
        ax.contour(
            (puncta_labels > 0).astype(np.uint8),
            [0.5],
            colors=["white"],
            linewidths=0.65,
        )
    bar_pixels = config.scale_bar_um / config.pixel_size_um
    bar_y = height - max(5, int(round(height * 0.06)))
    bar_x1 = width - max(5, int(round(width * 0.05)))
    bar_x0 = bar_x1 - bar_pixels
    ax.plot([bar_x0, bar_x1], [bar_y, bar_y], color="white", linewidth=2.0)
    ax.text(
        (bar_x0 + bar_x1) / 2,
        bar_y - max(3, height * 0.025),
        f"{config.scale_bar_um:g} µm",
        color="white",
        ha="center",
        va="bottom",
        fontsize=5,
    )
    handles = [
        Line2D([0], [0], color="cyan", lw=1, label="Nucleus"),
        Line2D([0], [0], color="white", lw=1, label="Puncta"),
    ]
    ax.legend(
        handles=handles,
        loc="lower left",
        frameon=False,
        fontsize=4,
        labelcolor="white",
        handlelength=1.4,
        borderaxespad=0.2,
    )
    ax.set_title(title, fontsize=5, pad=2, color="white")
    ax.set_axis_off()
    fig.subplots_adjust(left=0, right=1, bottom=0, top=0.94)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.png")
    fig.savefig(temporary, dpi=render_dpi, facecolor="black", bbox_inches=None)
    plt.close(fig)
    os.replace(temporary, path)


def _compact_overlay_title(condition: str, fov: str, nucleus_id: int) -> str:
    display_condition = {
        "before": "Before",
        "FVP_2h": "FVP 2 h",
        "recover_2h": "Recovery 2 h",
    }.get(condition, condition)
    suffix = fov.rsplit("-FOV", 1)[-1]
    short_fov = f"FOV{suffix}" if suffix else "FOV"
    return f"{display_condition} | {short_fov} | nucleus {nucleus_id}"


def _make_pilot_montage(
    path: Path,
    entries: Sequence[Mapping[str, Any]],
    display_limits: tuple[float, float],
    *,
    title: str = "Blinded puncta pilot: cyan nucleus, white puncta",
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    columns = 6
    rows = int(np.ceil(len(entries) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(18, 3 * rows), constrained_layout=True)
    lo, hi = display_limits
    for ax, entry in zip(np.ravel(axes), entries, strict=False):
        image = entry["image"]
        mask = entry["mask"]
        labels = entry["labels"]
        ax.imshow(np.clip(image, lo, hi), cmap="magma", norm=LogNorm(lo, hi))
        ax.contour(mask, [0.5], colors=["cyan"], linewidths=0.8)
        if np.any(labels):
            ax.contour(labels > 0, [0.5], colors=["white"], linewidths=0.7)
        nucleus_area = int(np.sum(mask))
        nucleus_area_fraction = (
            float(np.sum(labels > 0) / nucleus_area) if nucleus_area else float("nan")
        )
        ax.set_title(
            f"{entry['blinded_pilot_id']} | n={int(np.max(labels))} | "
            f"nucleus AF={nucleus_area_fraction:.3f}",
            fontsize=8,
        )
        ax.axis("off")
    for ax in np.ravel(axes)[len(entries) :]:
        ax.axis("off")
    fig.suptitle(title, fontsize=13)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.png")
    fig.savefig(temporary, dpi=250, facecolor="white")
    plt.close(fig)
    os.replace(temporary, path)


def prepare_analysis(
    config: LiveNuclearPunctaConfig,
    *,
    resume_staging: bool = False,
    reuse_segmentation_from: Path | None = None,
    allow_existing_output: bool = False,
) -> dict[str, Any]:
    reuse_root = reuse_segmentation_from.resolve() if reuse_segmentation_from else None
    if config.output_root.exists() and not allow_existing_output:
        raise FileExistsError(f"Final output already exists: {config.output_root}")
    if allow_existing_output and reuse_root != config.output_root.resolve():
        raise ValueError(
            "Existing output may only be retained during preparation when it is the "
            "explicit reusable-segmentation source"
        )
    if config.staging_root.exists() and not resume_staging:
        raise FileExistsError(
            f"Staging output already exists: {config.staging_root}; inspect before resuming or removing"
        )
    plans = discover_mips(config.dataset_root)
    staging = config.staging_root
    analysis = staging / "_analysis"
    staging.mkdir(parents=True, exist_ok=resume_staging)
    for condition in CONDITIONS:
        (staging / condition).mkdir(exist_ok=resume_staging)

    if resume_staging:
        saved_config_path = analysis / "config.json"
        if not saved_config_path.exists():
            raise FileNotFoundError(f"Cannot resume without {saved_config_path}")
        saved_config = json.loads(saved_config_path.read_text())
        if saved_config != config_to_json(config):
            raise ValueError("Staged configuration does not exactly match the requested configuration")
    else:
        _write_json(analysis / "config.json", config_to_json(config))
        _write_json(analysis / "discovery.json", discovery_summary(config))

    all_rows: list[dict[str, Any]] = []
    segmentation_provenance: list[dict[str, Any]] = []
    if resume_staging:
        for index, plan in enumerate(plans, start=1):
            label_path = _segmentation_label_path(plan, config, analysis)
            if not label_path.exists():
                raise FileNotFoundError(f"Cannot resume; missing staged label map: {label_path}")
            print(f"RESUME_SEGMENTATION {index}/{len(plans)} {plan.path.name}", flush=True)
            labels = tifffile.imread(label_path)
            expected_shape = _read_frame(plan.path, config.frame_index).shape
            if labels.shape != expected_shape or not np.issubdtype(labels.dtype, np.integer):
                raise ValueError(f"Invalid staged label map {label_path}: {labels.shape} {labels.dtype}")
            rows, _ = _save_segmentation_record(plan, labels, config, analysis)
            all_rows.extend(rows)
    elif reuse_root is not None:
        source_analysis = reuse_root / "_analysis"
        started = time.perf_counter()
        for index, plan in enumerate(plans, start=1):
            source_label = _segmentation_label_path(plan, config, source_analysis)
            if not source_label.is_file():
                raise FileNotFoundError(
                    f"Reusable Cellpose label map is missing: {source_label}"
                )
            labels = tifffile.imread(source_label)
            expected_shape = _read_frame(plan.path, config.frame_index).shape
            if labels.shape != expected_shape or not np.issubdtype(labels.dtype, np.integer):
                raise ValueError(
                    f"Invalid reusable label map {source_label}: "
                    f"{labels.shape} {labels.dtype}; expected {expected_shape} integer"
                )
            rows, destination = _save_segmentation_record(
                plan,
                labels,
                config,
                analysis,
            )
            source_hash = _sha256(source_label)
            destination_hash = _sha256(destination)
            if source_hash != destination_hash:
                raise ValueError(f"Reusable label hash mismatch for {plan.fov}")
            segmentation_provenance.append(
                {
                    "condition": plan.condition,
                    "fov": plan.fov,
                    "source_label": str(source_label),
                    "staged_label": str(destination),
                    "sha256": source_hash,
                    "shape": list(labels.shape),
                    "dtype": str(labels.dtype),
                    "label_count": int(np.count_nonzero(np.unique(labels) > 0)),
                }
            )
            all_rows.extend(rows)
            _progress_eta("REUSE_SEGMENTATION", index, len(plans), started)
        _write_json(
            analysis / "segmentation_reuse.json",
            {
                "reused": True,
                "cellpose_inference_rerun": False,
                "source_output_root": str(reuse_root),
                "source_run_state_sha256": (
                    _sha256(source_analysis / "run_state.json")
                    if (source_analysis / "run_state.json").is_file()
                    else None
                ),
                "label_maps": segmentation_provenance,
                "validated_at": utc_now(),
            },
        )
    else:
        model = _load_cellpose(config)
        for index, plan in enumerate(plans, start=1):
            print(f"SEGMENT {index}/{len(plans)} {plan.path.name}", flush=True)
            image = _read_frame(plan.path, config.frame_index)
            labels = run_cellpose_frame(image, model, config)
            rows, _ = _save_segmentation_record(plan, labels, config, analysis)
            all_rows.extend(rows)
    cutoff = pooled_area_cutoff([int(row["area_px"]) for row in all_rows], config.nucleus_area_quantile)
    retained, excluded = _retained_and_excluded(all_rows, cutoff)

    mask_pixels: list[np.ndarray] = []
    exported: list[dict[str, Any]] = []
    crop_started = time.perf_counter()
    for index, row in enumerate(retained, start=1):
        if index % 25 == 0 or index == len(retained):
            _progress_eta("EXPORT_CROPS", index, len(retained), crop_started)
        updated, pixels = _export_crop(row, config, staging)
        exported.append(updated)
        mask_pixels.append(pixels)
    pooled = np.concatenate(mask_pixels).astype(np.float64, copy=False)
    positive = pooled[np.isfinite(pooled) & (pooled > 0)]
    if not len(positive):
        raise ValueError("No positive retained nucleus pixels for logarithmic display")
    display_limits = tuple(
        float(value)
        for value in np.percentile(
            positive,
            [config.display_lower_percentile, config.display_upper_percentile],
        )
    )
    if display_limits[1] <= display_limits[0]:
        raise ValueError(f"Invalid global display limits: {display_limits}")

    pilot = select_blinded_pilot(
        exported,
        seed=config.pilot_seed,
        per_condition_per_stratum=config.pilot_per_condition_per_stratum,
    )
    pilot_entries: list[dict[str, Any]] = []
    pilot_rows: list[dict[str, Any]] = []
    pilot_mask_dir = analysis / "qc" / "default_pilot_masks"
    pilot_started = time.perf_counter()
    for row in pilot:
        crop = tifffile.imread(row["crop_tif"]).astype(np.float32, copy=False)
        image, mask = crop[0], crop[1] > 0.5
        labels = _run_default_puncta(image, mask, config)
        label_path = pilot_mask_dir / f"{row['blinded_pilot_id']}__labels-YX.tif"
        _write_imagej(label_path, labels, axes="YX", pixel_size_um=config.pixel_size_um)
        summary, _ = quantify_puncta(
            image,
            mask,
            labels,
            pixel_size_um=config.pixel_size_um,
        )
        pilot_rows.append(
            {
                "blinded_pilot_id": row["blinded_pilot_id"],
                "nucleus_key": row["nucleus_key"],
                "condition_hidden_from_montage": row["condition"],
                "pilot_stratum": row["pilot_stratum"],
                "puncta_count": summary["puncta_count"],
                "puncta_area_fraction": summary["puncta_area_fraction"],
                "mask_path": str(label_path),
            }
        )
        pilot_entries.append(
            {
                "blinded_pilot_id": row["blinded_pilot_id"],
                "image": image,
                "mask": mask,
                "labels": labels,
            }
        )
        _progress_eta(
            "DEFAULT_PILOT",
            len(pilot_rows),
            len(pilot),
            pilot_started,
        )
    montage = analysis / "qc" / "default_punctatools_pilot_montage.png"
    _make_pilot_montage(
        montage,
        pilot_entries,
        display_limits,
        title="Blinded PunctaTools-default pilot: cyan nucleus, white puncta",
    )
    area_fractions = np.asarray([row["puncta_area_fraction"] for row in pilot_rows], dtype=float)
    counts = np.asarray([row["puncta_count"] for row in pilot_rows], dtype=float)
    default_pilot_metrics = puncta_pilot_safety_metrics(counts, area_fractions)
    objective_pass = bool(default_pilot_metrics["passes_safety_gates"])
    _write_csv(analysis / "nucleus_manifest_pre_puncta.csv", exported)
    _write_csv(analysis / "excluded_nuclei.csv", excluded)
    _write_csv(analysis / "qc" / "default_pilot_worklist.csv", pilot_rows)
    state = {
        "status": "pilot_ready",
        "updated_at": utc_now(),
        "source_mip_count": len(plans),
        "segmented_nucleus_count": len(all_rows),
        "retained_nucleus_count": len(exported),
        "excluded_nucleus_count": len(excluded),
        "area_cutoff_px": cutoff,
        "display_limits_sacd": list(display_limits),
        "default_punctatools_params": asdict(default_puncta_params()),
        "fixed_threshold_enabled": config.fixed_threshold_enabled,
        "analysis_roi": "full_CellposeSAM_nucleus_mask",
        "nuclear_erosion_applied": False,
        "segmentation_reused": reuse_root is not None,
        "cellpose_inference_rerun": reuse_root is None and not resume_staging,
        "segmentation_label_map_count": len(plans),
        "segmentation_provenance": str(analysis / "segmentation_reuse.json")
        if reuse_root is not None
        else None,
        "pilot_nucleus_count": len(pilot_rows),
        "pilot_objective_pass": objective_pass,
        "default_pilot_safety_metrics": default_pilot_metrics,
        "pilot_montage": str(montage),
    }
    _write_json(analysis / "run_state.json", state)
    return state


def run_fixed_threshold_pilot(config: LiveNuclearPunctaConfig) -> dict[str, Any]:
    """Generate optional fixed-threshold pilot masks without changing primary results."""

    staging = config.staging_root
    analysis = staging / "_analysis"
    if not staging.is_dir():
        raise FileNotFoundError(staging)
    state_path = analysis / "run_state.json"
    state = json.loads(state_path.read_text())
    if state.get("status") != "pilot_ready":
        raise ValueError("Fixed-threshold pilot requires a pilot-ready staging run")
    nuclei = {row["nucleus_key"]: row for row in _read_csv(analysis / "nucleus_manifest_pre_puncta.csv")}
    pilot_rows = _read_csv(analysis / "qc" / "default_pilot_worklist.csv")
    global_background = pooled_nucleus_background(list(nuclei.values()))
    display_limits = tuple(float(value) for value in state["display_limits_sacd"])
    output_rows: list[dict[str, Any]] = []
    montage_paths: dict[str, str] = {}
    safety_by_multiplier: dict[str, dict[str, Any]] = {}
    pilot_started = time.perf_counter()
    completed = 0
    total = len(config.fixed_threshold_multipliers) * len(pilot_rows)
    for multiplier in config.fixed_threshold_multipliers:
        entries: list[dict[str, Any]] = []
        token = f"{multiplier:g}x".replace(".", "p")
        for pilot_row in pilot_rows:
            row = nuclei[pilot_row["nucleus_key"]]
            crop = tifffile.imread(row["crop_tif"]).astype(np.float32, copy=False)
            image, mask = crop[0], crop[1] > 0.5
            labels = punctatools_segment_2d(
                image,
                mask,
                pixel_size_um=config.pixel_size_um,
                params=fixed_threshold_params(multiplier),
                global_background=global_background,
            )
            label_path = (
                analysis
                / "qc"
                / "fixed_threshold_pilot"
                / token
                / f"{pilot_row['blinded_pilot_id']}__labels-YX.tif"
            )
            _write_imagej(label_path, labels, axes="YX", pixel_size_um=config.pixel_size_um)
            summary, _ = quantify_puncta(
                image,
                mask,
                labels,
                pixel_size_um=config.pixel_size_um,
            )
            output_rows.append(
                {
                    "blinded_pilot_id": pilot_row["blinded_pilot_id"],
                    "nucleus_key": row["nucleus_key"],
                    "condition_hidden_from_montage": row["condition"],
                    "multiplier": multiplier,
                    "global_background_sacd": global_background,
                    "puncta_count": summary["puncta_count"],
                    "puncta_area_fraction": summary["puncta_area_fraction"],
                    "mask_path": str(label_path),
                }
            )
            entries.append(
                {
                    "blinded_pilot_id": pilot_row["blinded_pilot_id"],
                    "image": image,
                    "mask": mask,
                    "labels": labels,
                }
            )
            completed += 1
            _progress_eta("FIXED_THRESHOLD_PILOT", completed, total, pilot_started)
        montage = analysis / "qc" / "fixed_threshold_pilot" / f"{token}__montage.png"
        _make_pilot_montage(
            montage,
            entries,
            display_limits,
            title=(
                "Blinded fixed-threshold sensitivity pilot "
                f"({multiplier:g}× pooled background)"
            ),
        )
        montage_paths[f"{multiplier:g}"] = str(montage)
        multiplier_rows = [
            row for row in output_rows if float(row["multiplier"]) == float(multiplier)
        ]
        area_fractions = np.asarray(
            [float(row["puncta_area_fraction"]) for row in multiplier_rows],
            dtype=float,
        )
        counts = np.asarray(
            [float(row["puncta_count"]) for row in multiplier_rows],
            dtype=float,
        )
        gate = puncta_pilot_safety_metrics(counts, area_fractions)
        safety_by_multiplier[f"{multiplier:g}"] = gate
    passing = [
        float(multiplier)
        for multiplier in config.fixed_threshold_multipliers
        if safety_by_multiplier[f"{multiplier:g}"]["passes_safety_gates"]
    ]
    selected_fallback = min(passing) if passing else None
    _write_csv(analysis / "qc" / "fixed_threshold_pilot" / "pilot_metrics.csv", output_rows)
    result = {
        "generated_at": utc_now(),
        "primary_results_changed": False,
        "global_background_sacd": global_background,
        "background_definition": (
            "median across retained cells of each full Cellpose nucleus median SACD"
        ),
        "multipliers": list(config.fixed_threshold_multipliers),
        "safety_gates": {
            "median_puncta_area_fraction_max": 0.35,
            "maximum_puncta_area_fraction_max": 0.80,
            "median_puncta_count_strict_max": 200,
        },
        "pilot_metrics_by_multiplier": safety_by_multiplier,
        "selected_fallback_multiplier": selected_fallback,
        "selection_rule": "lowest tested multiplier passing all safety gates",
        "analysis_roi": "full_CellposeSAM_nucleus_mask",
        "montages": montage_paths,
    }
    _write_json(analysis / "qc" / "fixed_threshold_pilot" / "summary.json", result)
    state["fixed_threshold_pilot"] = result
    state["updated_at"] = utc_now()
    _write_json(state_path, state)
    return result


def _fdr_bh(values: Sequence[float]) -> list[float]:
    pvalues = np.asarray(values, dtype=float)
    order = np.argsort(pvalues)
    adjusted = np.empty_like(pvalues)
    ranked = pvalues[order] * len(pvalues) / np.arange(1, len(pvalues) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted[order] = np.clip(ranked, 0, 1)
    return adjusted.tolist()


def _holm(values: Sequence[float]) -> list[float]:
    pvalues = np.asarray(values, dtype=float)
    order = np.argsort(pvalues)
    adjusted_sorted = np.maximum.accumulate(
        pvalues[order] * (len(pvalues) - np.arange(len(pvalues)))
    )
    adjusted = np.empty_like(pvalues)
    adjusted[order] = np.clip(adjusted_sorted, 0, 1)
    return adjusted.tolist()


def _rank_biserial(first: np.ndarray, second: np.ndarray, statistic: float) -> float:
    return float(2.0 * statistic / (len(first) * len(second)) - 1.0)


def _bootstrap_difference(
    first: np.ndarray,
    second: np.ndarray,
    *,
    seed: int,
    iterations: int,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    differences = np.empty(iterations, dtype=np.float64)
    for index in range(iterations):
        a = rng.choice(first, len(first), replace=True)
        b = rng.choice(second, len(second), replace=True)
        differences[index] = np.median(a) - np.median(b)
    low, high = np.percentile(differences, [2.5, 97.5])
    return float(low), float(high)


METRICS = (
    ("fov_mean_puncta_count", "Puncta count per nucleus"),
    ("fov_median_corrected_brightness", "Corrected puncta brightness (a.u.)"),
    ("fov_median_punctum_area_um2", "Punctum area (µm²)"),
)

CELL_METRICS = (
    ("puncta_count", "Puncta count per nucleus"),
    ("mean_corrected_puncta_brightness", "Corrected puncta brightness (a.u.)"),
    ("median_punctum_area_um2", "Median punctum area per nucleus (µm²)"),
)
CONDITION_COLORS = {
    "before": "#4C78A8",
    "FVP_2h": "#E45756",
    "recover_2h": "#54A24B",
}
CONDITION_DISPLAY_NAMES = {
    "before": "Before",
    "FVP_2h": "FVP 2 h",
    "recover_2h": "Recovery 2 h",
}
PRESENTATION_CONDITION_NAMES = {
    "before": "untreated",
    "FVP_2h": "transcription\ninhibited",
    "recover_2h": "recovered",
}
PLOT_FONT_SIZE = 11
PRESENTATION_FONT_SIZE = 22
PRESENTATION_FIGSIZE_INCHES = (8.5, 6.0)
PRESENTATION_BOXPLOT_METRICS = (
    "puncta_count",
    "mean_corrected_puncta_brightness",
)
PRESENTATION_Y_LABELS = {
    "puncta_count": "Puncta count\nper nucleus",
    "mean_corrected_puncta_brightness": "Puncta intensity, A.U.",
}
HISTOGRAM_BINS = 11


def _numeric_or_nan(value: Any) -> float:
    return float(value) if value not in ("", None) else float("nan")


def histogram_edges(
    values: Sequence[float],
    bins: int = HISTOGRAM_BINS,
    *,
    logarithmic: bool = False,
) -> np.ndarray:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if not len(finite):
        raise ValueError("Cannot construct histogram edges without finite values")
    minimum = float(np.min(finite))
    maximum = float(np.max(finite))
    if maximum == minimum:
        padding = max(abs(minimum) * 0.01, 0.5)
        minimum -= padding
        maximum += padding
    if logarithmic:
        if minimum <= 0:
            raise ValueError("Logarithmic histogram edges require strictly positive values")
        return np.geomspace(minimum, maximum, int(bins) + 1)
    return np.linspace(minimum, maximum, int(bins) + 1)


def aggregate_fovs(nucleus_rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], list[Mapping[str, Any]]] = {}
    for row in nucleus_rows:
        grouped.setdefault((str(row["condition"]), str(row["fov"])), []).append(row)
    output: list[dict[str, Any]] = []
    for (condition, fov), rows in sorted(grouped.items()):
        count_values = np.asarray([float(row["puncta_count"]) for row in rows])
        brightness = np.asarray(
            [float(row["mean_corrected_puncta_brightness"]) for row in rows],
            dtype=float,
        )
        area = np.asarray([float(row["median_punctum_area_um2"]) for row in rows], dtype=float)
        output.append(
            {
                "condition": condition,
                "fov": fov,
                "retained_nuclei": len(rows),
                "fov_mean_puncta_count": float(np.mean(count_values)),
                "fov_median_corrected_brightness": (
                    float(np.nanmedian(brightness)) if np.any(np.isfinite(brightness)) else float("nan")
                ),
                "fov_median_punctum_area_um2": (
                    float(np.nanmedian(area)) if np.any(np.isfinite(area)) else float("nan")
                ),
            }
        )
    return output


def calculate_statistics(
    fov_rows: Sequence[Mapping[str, Any]],
    *,
    seed: int,
    bootstrap_iterations: int,
) -> list[dict[str, Any]]:
    omnibus: list[dict[str, Any]] = []
    for metric, label in METRICS:
        arrays = []
        for condition in CONDITIONS:
            values = np.asarray(
                [
                    float(row[metric])
                    for row in fov_rows
                    if row["condition"] == condition and np.isfinite(float(row[metric]))
                ],
                dtype=float,
            )
            arrays.append(values)
        result = stats.kruskal(*arrays)
        omnibus.append(
            {
                "test_family": "omnibus",
                "metric": metric,
                "metric_label": label,
                "comparison": "before_vs_FVP_2h_vs_recover_2h",
                "statistic": float(result.statistic),
                "p_value": float(result.pvalue),
            }
        )
    omnibus_q = _fdr_bh([float(row["p_value"]) for row in omnibus])
    for row, value in zip(omnibus, omnibus_q, strict=True):
        row["adjusted_p_value"] = value

    pairs = (("before", "FVP_2h"), ("before", "recover_2h"), ("FVP_2h", "recover_2h"))
    significant_metrics = {
        str(row["metric"])
        for row in omnibus
        if float(row["adjusted_p_value"]) < 0.05
    }
    pairwise: list[dict[str, Any]] = []
    for metric, label in METRICS:
        if metric not in significant_metrics:
            continue
        metric_rows: list[dict[str, Any]] = []
        for pair_index, (first_condition, second_condition) in enumerate(pairs):
            first = np.asarray(
                [
                    float(row[metric])
                    for row in fov_rows
                    if row["condition"] == first_condition and np.isfinite(float(row[metric]))
                ]
            )
            second = np.asarray(
                [
                    float(row[metric])
                    for row in fov_rows
                    if row["condition"] == second_condition and np.isfinite(float(row[metric]))
                ]
            )
            test = stats.mannwhitneyu(first, second, alternative="two-sided", method="auto")
            ci_low, ci_high = _bootstrap_difference(
                first,
                second,
                seed=seed + pair_index,
                iterations=bootstrap_iterations,
            )
            metric_rows.append(
                {
                    "test_family": "pairwise",
                    "metric": metric,
                    "metric_label": label,
                    "comparison": f"{first_condition}_vs_{second_condition}",
                    "statistic": float(test.statistic),
                    "p_value": float(test.pvalue),
                    "rank_biserial_first_minus_second": _rank_biserial(
                        first,
                        second,
                        float(test.statistic),
                    ),
                    "median_difference_first_minus_second": float(
                        np.median(first) - np.median(second)
                    ),
                    "bootstrap_95_ci_low": ci_low,
                    "bootstrap_95_ci_high": ci_high,
                    "n_fov_first": len(first),
                    "n_fov_second": len(second),
                }
            )
        holm = _holm([float(row["p_value"]) for row in metric_rows])
        for row, value in zip(metric_rows, holm, strict=True):
            row["adjusted_p_value"] = value
        pairwise.extend(metric_rows)
    return omnibus + pairwise


def calculate_cell_statistics(
    nucleus_rows: Sequence[Mapping[str, Any]],
    *,
    seed: int,
    bootstrap_iterations: int,
) -> list[dict[str, Any]]:
    """Cell-level ANOVA gate followed by Holm-corrected Mann–Whitney tests."""

    import pandas as pd
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm

    omnibus: list[dict[str, Any]] = []
    metric_values: dict[str, dict[str, np.ndarray]] = {}
    for metric, label in CELL_METRICS:
        by_condition: dict[str, np.ndarray] = {}
        records: list[dict[str, Any]] = []
        for condition in CONDITIONS:
            values = np.asarray(
                [
                    _numeric_or_nan(row[metric])
                    for row in nucleus_rows
                    if row["condition"] == condition
                    and np.isfinite(_numeric_or_nan(row[metric]))
                ],
                dtype=float,
            )
            if len(values) < 2:
                raise ValueError(f"ANOVA requires at least two cells for {metric}/{condition}")
            by_condition[condition] = values
            records.extend({"condition": condition, "value": value} for value in values)
        metric_values[metric] = by_condition
        frame = pd.DataFrame.from_records(records)
        model = smf.ols("value ~ C(condition)", data=frame).fit()
        table = anova_lm(model, typ=2)
        effect = table.loc["C(condition)"]
        omnibus.append(
            {
                "independent_unit": "cell",
                "test_family": "omnibus",
                "test_name": "one_way_ANOVA_type_II_sums_of_squares",
                "metric": metric,
                "metric_label": label,
                "comparison": "before_vs_FVP_2h_vs_recover_2h",
                "statistic": float(effect["F"]),
                "degrees_of_freedom_between": float(effect["df"]),
                "degrees_of_freedom_residual": float(table.loc["Residual", "df"]),
                "p_value": float(effect["PR(>F)"]),
                "n_cell_before": len(by_condition["before"]),
                "n_cell_FVP_2h": len(by_condition["FVP_2h"]),
                "n_cell_recover_2h": len(by_condition["recover_2h"]),
            }
        )
    omnibus_q = _fdr_bh([float(row["p_value"]) for row in omnibus])
    for row, value in zip(omnibus, omnibus_q, strict=True):
        row["adjusted_p_value_BH_across_properties"] = value
        row["pairwise_triggered_by_raw_ANOVA_p_lt_0_05"] = bool(
            float(row["p_value"]) < 0.05
        )

    pairs = (("before", "FVP_2h"), ("before", "recover_2h"), ("FVP_2h", "recover_2h"))
    pairwise: list[dict[str, Any]] = []
    for omnibus_row in omnibus:
        if not bool(omnibus_row["pairwise_triggered_by_raw_ANOVA_p_lt_0_05"]):
            continue
        metric = str(omnibus_row["metric"])
        label = str(omnibus_row["metric_label"])
        metric_rows: list[dict[str, Any]] = []
        for pair_index, (first_condition, second_condition) in enumerate(pairs):
            first = metric_values[metric][first_condition]
            second = metric_values[metric][second_condition]
            test = stats.mannwhitneyu(
                first,
                second,
                alternative="two-sided",
                method="auto",
            )
            ci_low, ci_high = _bootstrap_difference(
                first,
                second,
                seed=seed + pair_index,
                iterations=bootstrap_iterations,
            )
            metric_rows.append(
                {
                    "independent_unit": "cell",
                    "test_family": "pairwise",
                    "test_name": "Mann-Whitney_two-sided",
                    "metric": metric,
                    "metric_label": label,
                    "comparison": f"{first_condition}_vs_{second_condition}",
                    "statistic": float(test.statistic),
                    "p_value": float(test.pvalue),
                    "rank_biserial_first_minus_second": _rank_biserial(
                        first,
                        second,
                        float(test.statistic),
                    ),
                    "median_difference_first_minus_second": float(
                        np.median(first) - np.median(second)
                    ),
                    "bootstrap_95_ci_low": ci_low,
                    "bootstrap_95_ci_high": ci_high,
                    "n_cell_first": len(first),
                    "n_cell_second": len(second),
                    "multiple_comparison_correction": "Holm",
                }
            )
        holm = _holm([float(row["p_value"]) for row in metric_rows])
        for row, value in zip(metric_rows, holm, strict=True):
            row["adjusted_p_value_Holm"] = value
        pairwise.extend(metric_rows)
    return omnibus + pairwise


def _significance_text(value: float) -> str:
    if value < 0.001:
        return "***"
    if value < 0.01:
        return "**"
    if value < 0.05:
        return "*"
    return "ns"


def make_distribution_histograms(
    nucleus_rows: Sequence[Mapping[str, Any]],
    path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    style = {
        "font.size": PLOT_FONT_SIZE,
        "axes.titlesize": PLOT_FONT_SIZE,
        "axes.labelsize": PLOT_FONT_SIZE,
        "xtick.labelsize": PLOT_FONT_SIZE,
        "ytick.labelsize": PLOT_FONT_SIZE,
        "legend.fontsize": PLOT_FONT_SIZE,
        "figure.titlesize": PLOT_FONT_SIZE,
    }
    with matplotlib.rc_context(style):
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
        for ax, (metric, label) in zip(axes, CELL_METRICS, strict=True):
            pooled = np.asarray(
                [
                    _numeric_or_nan(row[metric])
                    for row in nucleus_rows
                    if np.isfinite(_numeric_or_nan(row[metric]))
                ]
            )
            edges = histogram_edges(pooled, bins=HISTOGRAM_BINS)
            for condition in CONDITIONS:
                values = np.asarray(
                    [
                        _numeric_or_nan(row[metric])
                        for row in nucleus_rows
                        if row["condition"] == condition
                        and np.isfinite(_numeric_or_nan(row[metric]))
                    ]
                )
                ax.hist(
                    values,
                    bins=edges,
                    density=True,
                    histtype="step",
                    linewidth=1.8,
                    color=CONDITION_COLORS[condition],
                    label=f"{CONDITION_DISPLAY_NAMES[condition]} (n={len(values)})",
                )
            ax.set_xlim(float(edges[0]), float(edges[-1]))
            ax.margins(x=0)
            ax.set_xlabel(label)
            ax.set_ylabel("Density")
            ax.grid(alpha=0.2)
            ax.legend(frameon=False, fontsize=PLOT_FONT_SIZE)
            ax.tick_params(labelsize=PLOT_FONT_SIZE)
        fig.suptitle(
            "Frame-2 nucleus-level puncta distributions (descriptive)",
            fontsize=PLOT_FONT_SIZE,
        )
        temporary = path.with_name(path.name + ".tmp.png")
        fig.savefig(temporary, dpi=600, facecolor="white")
        plt.close(fig)
    os.replace(temporary, path)


def make_cell_boxplots(
    nucleus_rows: Sequence[Mapping[str, Any]],
    statistics_rows: Sequence[Mapping[str, Any]],
    path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import pandas as pd
    import seaborn as sns
    from statannotations.Annotator import Annotator

    style = {
        "font.size": PLOT_FONT_SIZE,
        "axes.titlesize": PLOT_FONT_SIZE,
        "axes.labelsize": PLOT_FONT_SIZE,
        "xtick.labelsize": PLOT_FONT_SIZE,
        "ytick.labelsize": PLOT_FONT_SIZE,
        "legend.fontsize": PLOT_FONT_SIZE,
        "figure.titlesize": PLOT_FONT_SIZE,
    }
    pair_order = [
        ("before", "FVP_2h"),
        ("before", "recover_2h"),
        ("FVP_2h", "recover_2h"),
    ]
    with matplotlib.rc_context(style):
        fig, axes = plt.subplots(1, 3, figsize=(15, 5.7), constrained_layout=True)
        for ax, (metric, label) in zip(axes, CELL_METRICS, strict=True):
            records = [
                {
                    "condition": str(row["condition"]),
                    "value": _numeric_or_nan(row[metric]),
                }
                for row in nucleus_rows
                if np.isfinite(_numeric_or_nan(row[metric]))
            ]
            frame = pd.DataFrame.from_records(records)
            sns.boxplot(
                data=frame,
                x="condition",
                y="value",
                order=list(CONDITIONS),
                palette=CONDITION_COLORS,
                hue="condition",
                hue_order=list(CONDITIONS),
                legend=False,
                showfliers=False,
                width=0.55,
                saturation=0.65,
                ax=ax,
            )
            sns.stripplot(
                data=frame,
                x="condition",
                y="value",
                order=list(CONDITIONS),
                hue="condition",
                hue_order=list(CONDITIONS),
                palette=CONDITION_COLORS,
                legend=False,
                jitter=0.20,
                alpha=0.72,
                edgecolor="black",
                linewidth=0.3,
                size=4,
                ax=ax,
            )
            relevant = [
                row
                for row in statistics_rows
                if row["test_family"] == "pairwise" and row["metric"] == metric
            ]
            if relevant:
                adjusted_by_comparison = {
                    str(row["comparison"]): float(row["adjusted_p_value_Holm"])
                    for row in relevant
                }
                adjusted = [
                    adjusted_by_comparison[f"{first}_vs_{second}"]
                    for first, second in pair_order
                ]
                annotator = Annotator(
                    ax,
                    pair_order,
                    data=frame,
                    x="condition",
                    y="value",
                    order=list(CONDITIONS),
                )
                annotator.configure(
                    text_format="star",
                    loc="inside",
                    show_test_name=False,
                    verbose=0,
                    fontsize=PLOT_FONT_SIZE,
                    line_width=1.0,
                )
                annotator.set_pvalues_and_annotate(adjusted)
            omnibus = next(
                row
                for row in statistics_rows
                if row["test_family"] == "omnibus" and row["metric"] == metric
            )
            ax.set_title(
                "One-way ANOVA "
                f"p={float(omnibus['p_value']):.3g}; "
                f"BH q={float(omnibus['adjusted_p_value_BH_across_properties']):.3g}",
                fontsize=PLOT_FONT_SIZE,
            )
            ax.set_xlabel("")
            ax.set_ylabel(PRESENTATION_Y_LABELS.get(metric, label) if presentation else label)
            ax.set_xticks(range(len(CONDITIONS)))
            ax.set_xticklabels(
                [CONDITION_DISPLAY_NAMES[condition] for condition in CONDITIONS],
                fontsize=PLOT_FONT_SIZE,
            )
            counts = frame.groupby("condition", observed=False).size().to_dict()
            handles = [
                Patch(
                    facecolor=CONDITION_COLORS[condition],
                    edgecolor="black",
                    alpha=0.65,
                    label=f"{CONDITION_DISPLAY_NAMES[condition]} (n={counts.get(condition, 0)})",
                )
                for condition in CONDITIONS
            ]
            ax.legend(
                handles=handles,
                frameon=False,
                fontsize=PLOT_FONT_SIZE,
                loc="center right",
            )
            ax.grid(axis="y", alpha=0.2)
            ax.tick_params(labelsize=PLOT_FONT_SIZE)
        fig.suptitle(
            "Cell-level puncta comparisons; Holm-corrected Mann–Whitney annotations",
            fontsize=PLOT_FONT_SIZE,
        )
        temporary = path.with_name(path.name + ".tmp.png")
        fig.savefig(temporary, dpi=600, facecolor="white")
        plt.close(fig)
    os.replace(temporary, path)


def make_individual_distribution_histograms(
    nucleus_rows: Sequence[Mapping[str, Any]],
    paths: Mapping[str, Path],
) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    style = {
        "font.size": PLOT_FONT_SIZE,
        "axes.titlesize": PLOT_FONT_SIZE,
        "axes.labelsize": PLOT_FONT_SIZE,
        "xtick.labelsize": PLOT_FONT_SIZE,
        "ytick.labelsize": PLOT_FONT_SIZE,
        "legend.fontsize": PLOT_FONT_SIZE,
        "figure.titlesize": PLOT_FONT_SIZE,
    }
    generated: list[Path] = []
    for metric, label in CELL_METRICS:
        path = Path(paths[metric])
        logarithmic = metric == "mean_corrected_puncta_brightness"
        pooled = np.asarray(
            [
                _numeric_or_nan(row[metric])
                for row in nucleus_rows
                if np.isfinite(_numeric_or_nan(row[metric]))
            ],
            dtype=float,
        )
        edges = histogram_edges(
            pooled,
            bins=HISTOGRAM_BINS,
            logarithmic=logarithmic,
        )
        with matplotlib.rc_context(style):
            fig, ax = plt.subplots(figsize=(8.0, 6.0))
            for condition in CONDITIONS:
                values = np.asarray(
                    [
                        _numeric_or_nan(row[metric])
                        for row in nucleus_rows
                        if row["condition"] == condition
                        and np.isfinite(_numeric_or_nan(row[metric]))
                    ],
                    dtype=float,
                )
                ax.hist(
                    values,
                    bins=edges,
                    density=True,
                    histtype="step",
                    linewidth=1.8,
                    color=CONDITION_COLORS[condition],
                    label=f"{CONDITION_DISPLAY_NAMES[condition]} (n={len(values)})",
                )
            if logarithmic:
                ax.set_xscale("log")
            ax.set_xlim(float(edges[0]), float(edges[-1]))
            ax.margins(x=0)
            ax.set_xlabel(label)
            ax.set_ylabel("Density")
            ax.set_title(f"Distribution of {label}", fontsize=PLOT_FONT_SIZE)
            ax.grid(alpha=0.2)
            ax.tick_params(labelsize=PLOT_FONT_SIZE)
            handles, labels = ax.get_legend_handles_labels()
            fig.legend(
                handles,
                labels,
                loc="lower center",
                bbox_to_anchor=(0.5, 0.015),
                ncol=3,
                frameon=False,
                fontsize=PLOT_FONT_SIZE,
            )
            fig.subplots_adjust(left=0.14, right=0.98, top=0.91, bottom=0.22)
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_name(path.name + ".tmp.png")
            fig.savefig(temporary, dpi=600, facecolor="white")
            plt.close(fig)
        os.replace(temporary, path)
        generated.append(path)
    return generated


def make_individual_cell_boxplots(
    nucleus_rows: Sequence[Mapping[str, Any]],
    statistics_rows: Sequence[Mapping[str, Any]],
    paths: Mapping[str, Path],
    *,
    metrics: Sequence[str] | None = None,
) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    import pandas as pd
    import seaborn as sns
    from statannotations.Annotator import Annotator

    style = {
        "font.size": PLOT_FONT_SIZE,
        "axes.titlesize": PLOT_FONT_SIZE,
        "axes.labelsize": PLOT_FONT_SIZE,
        "xtick.labelsize": PLOT_FONT_SIZE,
        "ytick.labelsize": PLOT_FONT_SIZE,
        "legend.fontsize": PLOT_FONT_SIZE,
        "figure.titlesize": PLOT_FONT_SIZE,
    }
    pair_order = [
        ("before", "FVP_2h"),
        ("before", "recover_2h"),
        ("FVP_2h", "recover_2h"),
    ]
    selected_metrics = set(metrics) if metrics is not None else None
    generated: list[Path] = []
    for metric, label in CELL_METRICS:
        if selected_metrics is not None and metric not in selected_metrics:
            continue
        path = Path(paths[metric])
        presentation = metric in PRESENTATION_BOXPLOT_METRICS
        font_size = PRESENTATION_FONT_SIZE if presentation else PLOT_FONT_SIZE
        condition_names = (
            PRESENTATION_CONDITION_NAMES if presentation else CONDITION_DISPLAY_NAMES
        )
        metric_style = {
            **style,
            "font.size": font_size,
            "axes.titlesize": font_size,
            "axes.labelsize": font_size,
            "xtick.labelsize": font_size,
            "ytick.labelsize": font_size,
            "legend.fontsize": font_size,
            "figure.titlesize": font_size,
        }
        records = [
            {
                "condition": str(row["condition"]),
                "value": _numeric_or_nan(row[metric]),
            }
            for row in nucleus_rows
            if np.isfinite(_numeric_or_nan(row[metric]))
        ]
        frame = pd.DataFrame.from_records(records)
        with matplotlib.rc_context(metric_style):
            fig, ax = plt.subplots(
                figsize=PRESENTATION_FIGSIZE_INCHES if presentation else (8.0, 6.2)
            )
            sns.boxplot(
                data=frame,
                x="condition",
                y="value",
                order=list(CONDITIONS),
                palette=CONDITION_COLORS,
                hue="condition",
                hue_order=list(CONDITIONS),
                legend=False,
                showfliers=False,
                width=0.55,
                saturation=0.65,
                ax=ax,
            )
            sns.stripplot(
                data=frame,
                x="condition",
                y="value",
                order=list(CONDITIONS),
                hue="condition",
                hue_order=list(CONDITIONS),
                palette=CONDITION_COLORS,
                legend=False,
                jitter=0.20,
                alpha=0.72,
                edgecolor="black",
                linewidth=0.3,
                size=4,
                ax=ax,
            )
            if metric == "mean_corrected_puncta_brightness":
                ax.set_yscale("log")
            relevant = [
                row
                for row in statistics_rows
                if row["test_family"] == "pairwise" and row["metric"] == metric
            ]
            if relevant:
                adjusted_by_comparison = {
                    str(row["comparison"]): float(row["adjusted_p_value_Holm"])
                    for row in relevant
                }
                adjusted = [
                    adjusted_by_comparison[f"{first}_vs_{second}"]
                    for first, second in pair_order
                ]
                annotator = Annotator(
                    ax,
                    pair_order,
                    data=frame,
                    x="condition",
                    y="value",
                    order=list(CONDITIONS),
                )
                annotator.configure(
                    text_format="star",
                    loc="inside",
                    show_test_name=False,
                    verbose=0,
                    fontsize=font_size,
                    line_width=1.0,
                )
                annotator.set_pvalues_and_annotate(adjusted)
            omnibus = next(
                row
                for row in statistics_rows
                if row["test_family"] == "omnibus" and row["metric"] == metric
            )
            if not presentation:
                ax.set_title(
                    "One-way ANOVA "
                    f"p={float(omnibus['p_value']):.3g}; "
                    f"BH q={float(omnibus['adjusted_p_value_BH_across_properties']):.3g}",
                    fontsize=font_size,
                )
            ax.set_xlabel("")
            ax.set_ylabel(
                PRESENTATION_Y_LABELS.get(metric, label) if presentation else label,
                fontsize=font_size,
            )
            ax.set_xticks(range(len(CONDITIONS)))
            ax.set_xticklabels(
                [condition_names[condition] for condition in CONDITIONS],
                fontsize=font_size,
            )
            ax.grid(axis="y", alpha=0.2)
            ax.tick_params(axis="both", labelsize=font_size)
            if presentation:
                fig.subplots_adjust(left=0.22, right=0.98, top=0.96, bottom=0.24)
            else:
                counts = frame.groupby("condition", observed=False).size().to_dict()
                handles = [
                    Patch(
                        facecolor=CONDITION_COLORS[condition],
                        edgecolor="black",
                        alpha=0.65,
                        label=(
                            f"{condition_names[condition]} "
                            f"(n={counts.get(condition, 0)})"
                        ),
                    )
                    for condition in CONDITIONS
                ]
                fig.legend(
                    handles=handles,
                    loc="lower center",
                    bbox_to_anchor=(0.5, 0.015),
                    ncol=3,
                    frameon=False,
                    fontsize=font_size,
                )
                fig.subplots_adjust(left=0.15, right=0.98, top=0.91, bottom=0.22)
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_name(path.name + ".tmp.png")
            fig.savefig(temporary, dpi=600, facecolor="white")
            plt.close(fig)
        os.replace(temporary, path)
        generated.append(path)
    return generated


def make_fov_boxplots(
    fov_rows: Sequence[Mapping[str, Any]],
    statistics_rows: Sequence[Mapping[str, Any]],
    path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"before": "#4C78A8", "FVP_2h": "#E45756", "recover_2h": "#54A24B"}
    rng = np.random.default_rng(20260722)
    display_names = {
        "before": "Before",
        "FVP_2h": "FVP 2 h",
        "recover_2h": "Recovery 2 h",
    }
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4), constrained_layout=True)
    for ax, (metric, label) in zip(axes, METRICS, strict=True):
        arrays = [
            np.asarray(
                [
                    float(row[metric])
                    for row in fov_rows
                    if row["condition"] == condition and np.isfinite(float(row[metric]))
                ]
            )
            for condition in CONDITIONS
        ]
        box = ax.boxplot(
            arrays,
            tick_labels=[display_names[condition] for condition in CONDITIONS],
            patch_artist=True,
            showfliers=False,
        )
        for patch, condition in zip(box["boxes"], CONDITIONS, strict=True):
            patch.set_facecolor(colors[condition])
            patch.set_alpha(0.45)
        for index, (values, condition) in enumerate(zip(arrays, CONDITIONS, strict=True), start=1):
            jitter = rng.uniform(-0.07, 0.07, size=len(values))
            ax.scatter(
                index + jitter,
                values,
                color=colors[condition],
                edgecolor="black",
                linewidth=0.35,
                s=28,
                zorder=3,
                label=f"{condition} n={len(values)}",
            )
        relevant = [
            row
            for row in statistics_rows
            if row["test_family"] == "pairwise" and row["metric"] == metric
        ]
        omnibus = next(
            row
            for row in statistics_rows
            if row["test_family"] == "omnibus" and row["metric"] == metric
        )
        pairwise_summary = "\n".join(
            f"{str(row['comparison']).replace('before', 'Before').replace('FVP_2h', 'FVP 2 h').replace('recover_2h', 'Recovery 2 h').replace('_vs_', ' vs ')} "
            f"{_significance_text(float(row['adjusted_p_value']))} "
            f"(q={float(row['adjusted_p_value']):.3g})"
            for row in relevant
        )
        title = (
            f"Kruskal–Wallis, BH q={float(omnibus['adjusted_p_value']):.3g}"
            + (f"\n{pairwise_summary}" if pairwise_summary else "\nPairwise tests not triggered")
        )
        ax.set_title(title, fontsize=7.5, pad=8)
        ax.set_ylabel(label)
        ax.grid(axis="y", alpha=0.2)
        ax.tick_params(axis="x", rotation=15)
    fig.suptitle(
        "FOV-level puncta comparisons (points are independent FOVs)",
        fontsize=13,
    )
    temporary = path.with_name(path.name + ".tmp.png")
    fig.savefig(temporary, dpi=600, facecolor="white")
    plt.close(fig)
    os.replace(temporary, path)


def _environment_manifest() -> dict[str, Any]:
    packages = (
        "cellpose",
        "numpy",
        "scipy",
        "scikit-image",
        "tifffile",
        "matplotlib",
        "pandas",
        "seaborn",
        "statannotations",
        "statsmodels",
    )
    versions: dict[str, str] = {}
    for package in packages:
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return {
        "captured_at": utc_now(),
        "python": sys.version,
        "platform": platform.platform(),
        "packages": versions,
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _snapshot_sources(
    analysis_root: Path,
    *,
    config_path: Path | None,
    script_path: Path | None,
) -> None:
    snapshot = analysis_root / "source_snapshot"
    snapshot.mkdir(parents=True, exist_ok=True)
    candidates = [Path(__file__), Path(__file__).with_name("puncta.py")]
    if config_path:
        candidates.append(config_path)
    if script_path:
        candidates.append(script_path)
    rows = []
    for source in candidates:
        if not source.is_file():
            continue
        destination = snapshot / source.name
        shutil.copy2(source, destination)
        rows.append(
            {
                "source": str(source),
                "snapshot": str(destination),
                "sha256": _sha256(destination),
            }
        )
    _write_json(snapshot / "snapshot_manifest.json", rows)


def refresh_published_qc(
    config: LiveNuclearPunctaConfig,
    *,
    config_path: Path | None = None,
    script_path: Path | None = None,
) -> dict[str, Any]:
    """Regenerate render-only artifacts without changing measurements or masks."""

    root = config.output_root
    analysis = root / "_analysis"
    if not root.is_dir():
        raise FileNotFoundError(root)
    selected = json.loads((analysis / "selected_detector.json").read_text())
    result_token = str(selected["result_token"])
    state_path = analysis / "run_state.json"
    state = json.loads(state_path.read_text())
    display_limits = tuple(float(value) for value in state["display_limits_sacd"])
    nucleus_path = root / f"nucleus_manifest__{result_token}.csv"
    fov_path = root / f"fov_summary__{result_token}.csv"
    statistics_path = root / f"statistical_tests__{result_token}.csv"
    nuclei = _read_csv(nucleus_path)
    for index, row in enumerate(nuclei, start=1):
        if index % 25 == 0 or index == len(nuclei):
            print(f"REFRESH_OVERLAY {index}/{len(nuclei)}", flush=True)
        crop = tifffile.imread(row["crop_tif"]).astype(np.float32, copy=False)
        labels = tifffile.imread(row["puncta_mask_path"])
        image, mask = crop[0], crop[1] > 0.5
        _make_overlay(
            Path(row["overlay_png"]),
            image,
            mask,
            labels,
            title=_compact_overlay_title(
                str(row["condition"]),
                str(row["fov"]),
                int(row["nucleus_id"]),
            ),
            display_limits=display_limits,
            config=config,
        )
    replot = replot_published_cell_statistics(
        config,
        config_path=config_path,
        script_path=script_path,
    )
    _snapshot_sources(analysis, config_path=config_path, script_path=script_path)
    result = {
        "refreshed_at": utc_now(),
        "render_only": True,
        "nucleus_overlays_refreshed": len(nuclei),
        "summary_plots_refreshed": 6,
        "cell_level_plot_refresh": replot,
        "measurements_changed": False,
        "masks_changed": False,
    }
    _write_json(analysis / "qc_refresh.json", result)
    state["qc_refresh"] = result
    _write_json(state_path, state)
    return result


def replot_published_cell_statistics(
    config: LiveNuclearPunctaConfig,
    *,
    config_path: Path | None = None,
    script_path: Path | None = None,
    presentation_only: bool = False,
) -> dict[str, Any]:
    """Publish cell-level plots, optionally limiting writes to presentation plots."""

    root = config.output_root
    analysis = root / "_analysis"
    if not root.is_dir():
        raise FileNotFoundError(root)
    selected = json.loads((analysis / "selected_detector.json").read_text())
    result_token = str(selected["result_token"])
    state_path = analysis / "run_state.json"
    state = json.loads(state_path.read_text())
    nucleus_path = root / f"nucleus_manifest__{result_token}.csv"
    nuclei = _read_csv(nucleus_path)
    statistics_path = root / f"statistical_tests__{result_token}__cell_level.csv"
    if presentation_only:
        if not statistics_path.is_file():
            raise FileNotFoundError(statistics_path)
        statistics_rows = _read_csv(statistics_path)
    else:
        statistics_rows = calculate_cell_statistics(
            nuclei,
            seed=config.pilot_seed,
            bootstrap_iterations=config.bootstrap_iterations,
        )
    filename_tokens = {
        "puncta_count": "puncta_count",
        "mean_corrected_puncta_brightness": "corrected_puncta_brightness",
        "median_punctum_area_um2": "median_punctum_area",
    }
    distribution_paths = {
        metric: root / f"{token}__distribution__{result_token}.png"
        for metric, token in filename_tokens.items()
    }
    boxplot_paths = {
        metric: root / f"{token}__cell_boxplot_statistics__{result_token}.png"
        for metric, token in filename_tokens.items()
    }
    if presentation_only:
        generated_distributions: list[Path] = []
    else:
        _write_csv(statistics_path, _sanitize_rows(statistics_rows))
        generated_distributions = make_individual_distribution_histograms(
            nuclei,
            distribution_paths,
        )
    generated_boxplots = make_individual_cell_boxplots(
        nuclei,
        statistics_rows,
        boxplot_paths,
        metrics=PRESENTATION_BOXPLOT_METRICS if presentation_only else None,
    )

    legacy_candidates = [
        root / f"puncta_property_FOV_boxplots_statistics__{result_token}.png",
        root / f"puncta_property_distributions__{result_token}.png",
        root / f"puncta_property_cell_boxplots_statistics__{result_token}.png",
    ]
    legacy_destinations: list[str] = []
    if not presentation_only:
        for legacy_plot in legacy_candidates:
            legacy_destination = analysis / "legacy_outputs" / legacy_plot.name
            if legacy_plot.exists():
                legacy_destination.parent.mkdir(parents=True, exist_ok=True)
                os.replace(legacy_plot, legacy_destination)
            if legacy_destination.exists():
                legacy_destinations.append(str(legacy_destination))

    _write_json(analysis / "environment.json", _environment_manifest())
    _snapshot_sources(analysis, config_path=config_path, script_path=script_path)
    result = {
        "replotted_at": utc_now(),
        "render_scope": (
            "presentation_boxplots_only" if presentation_only else "all_summary_plots"
        ),
        "independent_unit": "cell",
        "retained_cell_count": len(nuclei),
        "histogram_bins": HISTOGRAM_BINS,
        "plot_font_size_points": {
            "histograms_and_area_boxplot": PLOT_FONT_SIZE,
            "presentation_count_and_brightness_boxplots": PRESENTATION_FONT_SIZE,
        },
        "histogram_paths": [str(path) for path in generated_distributions],
        "boxplot_paths": [str(path) for path in generated_boxplots],
        "intensity_axis_scale": "logarithmic",
        "count_and_area_axis_scale": "linear",
        "legend_layout": {
            "presentation_boxplots": "none",
            "other_plots": "bottom, three columns",
        },
        "presentation_boxplot_style": {
            "metrics": list(PRESENTATION_BOXPLOT_METRICS),
            "figsize_inches": list(PRESENTATION_FIGSIZE_INCHES),
            "condition_labels": PRESENTATION_CONDITION_NAMES,
            "y_axis_labels": PRESENTATION_Y_LABELS,
            "legend_removed_as_redundant": True,
            "statistical_number_title_removed": True,
            "significance_brackets_and_stars_retained": True,
        },
        "statistics_path": str(statistics_path),
        "omnibus_test": "one-way ANOVA with type-II sums of squares",
        "pairwise_gate": "raw omnibus ANOVA p < 0.05 for each property",
        "pairwise_test": "two-sided Mann-Whitney",
        "pairwise_correction": "Holm within each triggered property",
        "annotations_package": "statannotations",
        "measurements_changed": False,
        "masks_changed": False,
        "statistics_changed": not presentation_only,
        "legacy_combined_plots": legacy_destinations,
    }
    _write_json(analysis / "cell_level_replot.json", result)
    if presentation_only:
        all_summary_paths = list(distribution_paths.values()) + list(boxplot_paths.values())
        state["active_summary_plots"] = [
            str(path) for path in all_summary_paths if path.is_file()
        ]
    else:
        state["active_summary_plots"] = [
            str(path) for path in generated_distributions + generated_boxplots
        ]
    state["active_statistics"] = str(statistics_path)
    state["cell_level_replot"] = result
    _write_json(state_path, state)
    return result


def _rewrite_staging_paths(root: Path, old_root: Path, new_root: Path) -> None:
    """Rewrite staged absolute paths in text manifests immediately before publication."""

    old = str(old_root)
    new = str(new_root)
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in {".csv", ".json"}:
            continue
        original = path.read_text()
        updated = original.replace(old, new)
        if updated == original:
            continue
        temporary = path.with_name(path.name + ".rewrite.tmp")
        temporary.write_text(updated)
        os.replace(temporary, path)


def _finite_or_none(value: Any) -> Any:
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _sanitize_rows(rows: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    return [
        {key: _finite_or_none(value) for key, value in row.items()}
        for row in rows
    ]


def finalize_analysis(
    config: LiveNuclearPunctaConfig,
    *,
    approve_pilot: bool,
    detector_mode: str = "default",
    fixed_threshold_multiplier: float | None = None,
    config_path: Path | None = None,
    script_path: Path | None = None,
    replace_existing: bool = False,
) -> dict[str, Any]:
    staging = config.staging_root
    analysis = staging / "_analysis"
    if config.output_root.exists() and not replace_existing:
        raise FileExistsError(config.output_root)
    if not staging.is_dir():
        raise FileNotFoundError(staging)
    state = json.loads((analysis / "run_state.json").read_text())
    if state.get("status") != "pilot_ready":
        raise ValueError(f"Staging is not pilot-ready: {state.get('status')}")
    if not approve_pilot:
        raise ValueError("Full processing requires explicit --approve-pilot")
    if detector_mode not in {"default", "fixed-threshold"}:
        raise ValueError(f"Unknown detector mode: {detector_mode}")
    global_background: float | None = None
    if detector_mode == "default":
        if not bool(state.get("pilot_objective_pass")):
            raise ValueError("Default PunctaTools pilot failed objective safety gates")
        detector_token = "default"
        result_token = "default_full_nucleus"
        detector_name = "PunctaTools_2D_v0.2.0_compatible_defaults_full_nucleus"
        detector_parameters = {
            **asdict(default_puncta_params()),
            "analysis_roi": "full_CellposeSAM_nucleus_mask",
            "nuclear_erosion_applied": False,
        }
    else:
        fixed_pilot = state.get("fixed_threshold_pilot")
        if not isinstance(fixed_pilot, dict):
            raise ValueError("Fixed-threshold full analysis requires saved pilot diagnostics")
        if fixed_threshold_multiplier is None:
            raise ValueError("Fixed-threshold full analysis requires one frozen multiplier")
        candidates = np.asarray(config.fixed_threshold_multipliers, dtype=float)
        if not np.any(np.isclose(candidates, fixed_threshold_multiplier)):
            raise ValueError(
                f"Multiplier {fixed_threshold_multiplier} was not in the pilot set "
                f"{config.fixed_threshold_multipliers}"
            )
        selected_fallback = fixed_pilot.get("selected_fallback_multiplier")
        if selected_fallback is None:
            raise ValueError("No fixed-threshold multiplier passed the pilot safety gates")
        pilot_metrics = fixed_pilot.get("pilot_metrics_by_multiplier", {})
        requested_key = f"{fixed_threshold_multiplier:g}"
        requested_metrics = pilot_metrics.get(requested_key)
        if not isinstance(requested_metrics, dict) or not bool(
            requested_metrics.get("passes_safety_gates")
        ):
            raise ValueError(
                f"Requested multiplier {fixed_threshold_multiplier:g} did not pass "
                "the predeclared pilot safety gates"
            )
        global_background = float(fixed_pilot["global_background_sacd"])
        multiplier_token = f"{fixed_threshold_multiplier:g}x".replace(".", "p")
        detector_token = f"fixed_threshold/{multiplier_token}"
        result_token = f"fixed_threshold_{multiplier_token}_full_nucleus_sensitivity"
        detector_name = (
            "PunctaTools_2D_v0.2.0_compatible_fixed_global_threshold_"
            "full_nucleus_sensitivity"
        )
        detector_parameters = {
            **asdict(fixed_threshold_params(fixed_threshold_multiplier)),
            "fixed_threshold_multiplier": float(fixed_threshold_multiplier),
            "pooled_global_background_sacd": global_background,
            "absolute_intensity_threshold_sacd": (
                global_background * float(fixed_threshold_multiplier)
            ),
            "analysis_roi": "full_CellposeSAM_nucleus_mask",
            "nuclear_erosion_applied": False,
            "pilot_safety_metrics": requested_metrics,
            "automatic_lowest_passing_multiplier": float(selected_fallback),
            "selection_basis": "explicit_user_selection_among_passing_multipliers",
        }
    rows = _read_csv(analysis / "nucleus_manifest_pre_puncta.csv")
    display_limits = tuple(float(value) for value in state["display_limits_sacd"])
    puncta_dir = analysis / "puncta_masks" / detector_token
    nucleus_output: list[dict[str, Any]] = []
    puncta_output: list[dict[str, Any]] = []
    puncta_started = time.perf_counter()
    for index, source in enumerate(rows, start=1):
        crop_path = Path(source["crop_tif"])
        crop = tifffile.imread(crop_path).astype(np.float32, copy=False)
        image, mask = crop[0], crop[1] > 0.5
        if detector_mode == "default":
            labels = _run_default_puncta(image, mask, config)
        else:
            labels = punctatools_segment_2d(
                image,
                mask,
                pixel_size_um=config.pixel_size_um,
                params=fixed_threshold_params(float(fixed_threshold_multiplier)),
                global_background=global_background,
            )
        label_path = puncta_dir / source["condition"] / (
            f"{source['nucleus_key']}__PunctaTools-{result_token}-labels-YX.tif"
        )
        _write_imagej(label_path, labels, axes="YX", pixel_size_um=config.pixel_size_um)
        summary, puncta_rows = quantify_puncta(
            image,
            mask,
            labels,
            pixel_size_um=config.pixel_size_um,
        )
        png_path = crop_path.with_suffix(".png")
        _make_overlay(
            png_path,
            image,
            mask,
            labels,
            title=_compact_overlay_title(
                str(source["condition"]),
                str(source["fov"]),
                int(source["nucleus_id"]),
            ),
            display_limits=display_limits,
            config=config,
        )
        nucleus_row: dict[str, Any] = {
            **source,
            **summary,
            "puncta_detector": detector_name,
            "puncta_result_mode": result_token,
            "puncta_mask_path": str(label_path),
            "overlay_png": str(png_path),
        }
        nucleus_output.append(nucleus_row)
        for punctum in puncta_rows:
            puncta_output.append(
                {
                    "nucleus_key": source["nucleus_key"],
                    "condition": source["condition"],
                    "fov": source["fov"],
                    "puncta_result_mode": result_token,
                    **punctum,
                }
            )
        _progress_eta("PUNCTA", index, len(rows), puncta_started)
    fov_rows = aggregate_fovs(nucleus_output)
    cell_stats_rows = calculate_cell_statistics(
        nucleus_output,
        seed=config.pilot_seed,
        bootstrap_iterations=config.bootstrap_iterations,
    )
    nucleus_manifest_path = staging / f"nucleus_manifest__{result_token}.csv"
    puncta_manifest_path = staging / f"puncta_manifest__{result_token}.csv"
    fov_summary_path = staging / f"fov_summary__{result_token}.csv"
    statistical_tests_path = (
        staging / f"statistical_tests__{result_token}__cell_level.csv"
    )
    _write_csv(nucleus_manifest_path, _sanitize_rows(nucleus_output))
    _write_csv(staging / "excluded_nuclei.csv", _read_csv(analysis / "excluded_nuclei.csv"))
    _write_csv(puncta_manifest_path, _sanitize_rows(puncta_output))
    _write_csv(fov_summary_path, _sanitize_rows(fov_rows))
    _write_csv(statistical_tests_path, _sanitize_rows(cell_stats_rows))
    filename_tokens = {
        "puncta_count": "puncta_count",
        "mean_corrected_puncta_brightness": "corrected_puncta_brightness",
        "median_punctum_area_um2": "median_punctum_area",
    }
    distribution_paths = {
        metric: staging / f"{token}__distribution__{result_token}.png"
        for metric, token in filename_tokens.items()
    }
    boxplot_paths = {
        metric: staging / f"{token}__cell_boxplot_statistics__{result_token}.png"
        for metric, token in filename_tokens.items()
    }
    plot_started = time.perf_counter()
    generated_distributions = make_individual_distribution_histograms(
        nucleus_output,
        distribution_paths,
    )
    _progress_eta("SUMMARY_PLOTS", 3, 6, plot_started)
    generated_boxplots = make_individual_cell_boxplots(
        nucleus_output,
        cell_stats_rows,
        boxplot_paths,
    )
    _progress_eta("SUMMARY_PLOTS", 6, 6, plot_started)
    plot_paths = generated_distributions + generated_boxplots
    _write_json(
        analysis / "selected_detector.json",
        {
            "selected_at": utc_now(),
            "mode": detector_mode,
            "result_token": result_token,
            "name": detector_name,
            "parameters": detector_parameters,
            "default_pilot_objective_pass": bool(state.get("pilot_objective_pass")),
            "analysis_roi": "full_CellposeSAM_nucleus_mask",
            "nuclear_erosion_applied": False,
            "selection_note": (
                "Default mode retained after blinded pilot."
                if detector_mode == "default"
                else (
                    "User-selected passing fixed-threshold sensitivity mode after the "
                    "default pilot labeled nearly all full-nucleus foreground."
                )
            ),
        },
    )
    _write_json(analysis / "environment.json", _environment_manifest())
    _snapshot_sources(analysis, config_path=config_path, script_path=script_path)
    validation = validate_staging(
        config,
        nucleus_manifest_path=nucleus_manifest_path,
        fov_summary_path=fov_summary_path,
        expected_plot_paths=plot_paths,
    )
    completed = {
        **state,
        "status": "complete",
        "completed_at": utc_now(),
        "pilot_approved": True,
        "puncta_detector": detector_name,
        "puncta_result_mode": result_token,
        "detector_parameters": detector_parameters,
        "analysis_roi": "full_CellposeSAM_nucleus_mask",
        "nuclear_erosion_applied": False,
        "nucleus_count": len(nucleus_output),
        "punctum_count": len(puncta_output),
        "fov_count": len(fov_rows),
        "active_statistics": str(statistical_tests_path),
        "active_summary_plots": [str(path) for path in plot_paths],
        "cell_level_plotting": {
            "independent_unit": "cell",
            "histogram_bins": HISTOGRAM_BINS,
            "plot_font_size_points": PLOT_FONT_SIZE,
            "intensity_axis_scale": "logarithmic",
            "count_and_area_axis_scale": "linear",
            "legend_layout": "bottom, three columns",
            "omnibus_test": "one-way ANOVA with type-II sums of squares",
            "pairwise_gate": "raw omnibus ANOVA p < 0.05 for each property",
            "pairwise_test": "two-sided Mann-Whitney",
            "pairwise_correction": "Holm within each triggered property",
            "annotations_package": "statannotations",
        },
        "validation": validation,
    }
    _write_json(analysis / "run_state.json", completed)
    _rewrite_staging_paths(staging, staging, config.output_root)
    replacement_backup: Path | None = None
    if config.output_root.exists():
        replacement_backup = config.dataset_root / (
            f".{config.output_folder}.replacing-old-{uuid.uuid4().hex}"
        )
        os.replace(config.output_root, replacement_backup)
    try:
        os.replace(staging, config.output_root)
        published_validation = validate_result_root(
            config,
            root=config.output_root,
            nucleus_manifest_path=(
                config.output_root / nucleus_manifest_path.name
            ),
            fov_summary_path=config.output_root / fov_summary_path.name,
            expected_plot_paths=[
                config.output_root / path.name for path in plot_paths
            ],
        )
    except Exception:
        if config.output_root.exists():
            os.replace(config.output_root, staging)
        if replacement_backup is not None and replacement_backup.exists():
            os.replace(replacement_backup, config.output_root)
        raise
    if replacement_backup is not None:
        shutil.rmtree(replacement_backup)
    published_state_path = config.output_root / "_analysis" / "run_state.json"
    published_state = json.loads(published_state_path.read_text())
    published_state["published_validation"] = published_validation
    published_state["legacy_erosion_result_deleted"] = replacement_backup is not None
    published_state["published_at"] = utc_now()
    _write_json(published_state_path, published_state)
    return published_state


def validate_staging(
    config: LiveNuclearPunctaConfig,
    *,
    nucleus_manifest_path: Path | None = None,
    fov_summary_path: Path | None = None,
    expected_plot_paths: Sequence[Path] = (),
) -> dict[str, Any]:
    staging = config.staging_root
    return validate_result_root(
        config,
        root=staging,
        nucleus_manifest_path=nucleus_manifest_path or staging / "nucleus_manifest.csv",
        fov_summary_path=fov_summary_path or staging / "fov_summary.csv",
        expected_plot_paths=expected_plot_paths,
    )


def validate_result_root(
    config: LiveNuclearPunctaConfig,
    *,
    root: Path,
    nucleus_manifest_path: Path,
    fov_summary_path: Path,
    expected_plot_paths: Sequence[Path] = (),
) -> dict[str, Any]:
    analysis = root / "_analysis"
    nuclei = _read_csv(nucleus_manifest_path)
    fovs = _read_csv(fov_summary_path)
    failures: list[str] = []
    if len(fovs) != 21:
        failures.append(f"Expected 21 FOV summaries, found {len(fovs)}")
    for row in nuclei:
        crop_path = Path(row["crop_tif"])
        png_path = Path(row["overlay_png"])
        puncta_path = Path(row["puncta_mask_path"])
        if not crop_path.is_file() or not png_path.is_file() or not puncta_path.is_file():
            failures.append(f"Missing nucleus artifact for {row['nucleus_key']}")
            continue
        crop = tifffile.imread(crop_path)
        puncta = tifffile.imread(puncta_path)
        if crop.ndim != 3 or crop.shape[0] != 2 or crop.dtype != np.float32:
            failures.append(f"Invalid CYX crop for {row['nucleus_key']}: {crop.shape} {crop.dtype}")
        nucleus_mask = crop[1] > 0.5
        if puncta.shape != crop.shape[1:] or np.any(puncta[~nucleus_mask] > 0):
            failures.append(f"Invalid puncta-mask alignment for {row['nucleus_key']}")
        forbidden = [
            key
            for key in row
            if "core" in key.lower() or "erod" in key.lower()
        ]
        if forbidden:
            failures.append(
                f"Erosion-specific manifest fields for {row['nucleus_key']}: {forbidden}"
            )
    missing_plots = [str(path) for path in expected_plot_paths if not path.is_file()]
    failures.extend(f"Missing summary plot: {path}" for path in missing_plots)
    reuse_path = analysis / "segmentation_reuse.json"
    if not reuse_path.is_file():
        failures.append(f"Missing segmentation reuse provenance: {reuse_path}")
        reused_labels: list[dict[str, Any]] = []
    else:
        reuse = json.loads(reuse_path.read_text())
        reused_labels = list(reuse.get("label_maps", []))
        if reuse.get("cellpose_inference_rerun") is not False:
            failures.append("Segmentation provenance does not confirm Cellpose reuse")
        if len(reused_labels) != 21:
            failures.append(
                f"Expected 21 reused Cellpose label maps, found {len(reused_labels)}"
            )
        for record in reused_labels:
            label_path = Path(str(record["staged_label"]))
            if not label_path.is_file() or _sha256(label_path) != record["sha256"]:
                failures.append(f"Reusable label hash validation failed: {label_path}")
    if failures:
        raise ValueError("Staging validation failed:\n" + "\n".join(failures[:20]))
    return {
        "nucleus_artifacts_validated": len(nuclei),
        "fov_summaries_validated": len(fovs),
        "segmentation_label_maps_validated": len(reused_labels),
        "summary_plots_validated": len(expected_plot_paths),
        "analysis_roi": "full_CellposeSAM_nucleus_mask",
        "nuclear_erosion_applied": False,
        "condition_fov_counts": {
            condition: sum(row["condition"] == condition for row in fovs)
            for condition in CONDITIONS
        },
        "failures": 0,
        "validated_at": utc_now(),
    }


def config_to_json(config: LiveNuclearPunctaConfig) -> dict[str, Any]:
    value = asdict(config)
    value["dataset_root"] = str(config.dataset_root)
    value["fixed_threshold_multipliers"] = list(config.fixed_threshold_multipliers)
    return value


def config_from_json(path: str | Path) -> LiveNuclearPunctaConfig:
    value = json.loads(Path(path).read_text())
    value["dataset_root"] = Path(value["dataset_root"])
    value.pop("core_radius_fraction", None)
    value["fixed_threshold_multipliers"] = tuple(value.get("fixed_threshold_multipliers", (1.25, 1.5, 2.0, 3.0)))
    value["punctatools"] = PunctaTools2DParams(**value.get("punctatools", {}))
    return LiveNuclearPunctaConfig(**value)
