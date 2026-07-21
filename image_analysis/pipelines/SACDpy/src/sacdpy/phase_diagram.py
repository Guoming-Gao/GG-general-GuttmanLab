from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import uuid
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable

import numpy as np
import tifffile
from scipy.ndimage import distance_transform_edt

from .discrete_timelapse import read_oni_metadata
from .params import SACDParams
from .reconstruction import reconstruct


DEFAULT_DATASET_ROOT = Path(
    "/Volumes/guttman/primarydata/imaging_rawdata/"
    "20260701_ONI-gmgao-Dox_phase_diagram_FL_dRRM"
)
DEFAULT_OUTPUT_ROOT = Path(
    "/Volumes/guttman/users/gmgao/Imaging_ProcessedData/SPEN/"
    "SACDlive-phase diagram"
)

_ONI_SINGLE_TIME_RE = re.compile(
    r"^(?P<prefix>.+)_posXY(?P<pos_xy>\d+)_channels_t(?P<time_index>\d+)_posZ(?P<z_index>\d+)\.tiff?$",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class PhaseDiagramConfig:
    dataset_root: Path = DEFAULT_DATASET_ROOT
    output_root: Path = DEFAULT_OUTPUT_ROOT
    position_folder: str = "pos_0"
    glob_pattern: str = "*.tif"
    hoechst_wavelength_nm: float = 461.0
    spen_wavelength_nm: float = 647.0
    fallback_pixel_nm: float = 117.0
    fallback_na: float = 1.45
    mag: int = 2
    iter1: int = 7
    iter2: int = 8
    ac_order: int = 2
    subfactor: float = 0.8
    ifbackground: bool = False
    backgroundfactor: float = 2.0
    ifregistration: bool = False
    ifsparsedecon: bool = False
    fidelity: float = 100.0
    tcontinuity: float = 0.1
    sparsity: float = 1.0
    sparse_iterations: int = 100
    cellpose_diameter_px: float = 140.0
    cellpose_flow_threshold: float = 0.4
    cellpose_cellprob_threshold: float = 0.0
    cellpose_lower_percentile: float = 1.0
    cellpose_upper_percentile: float = 99.8
    cellpose_device: str = "auto"
    crop_padding_px: int = 10
    review_scale_bar_um: float = 2.0
    review_lower_percentile: float = 1.0
    review_upper_percentile: float = 99.8
    core_erosion_radius_fraction: float = 1.0 / 3.0
    overwrite: bool = False


@dataclass(frozen=True)
class PhaseFOVPlan:
    fov_folder: Path
    fov_name: str
    condition: str
    prefix: str
    pos_xy: int
    files: dict[int, Path]

    @property
    def z_indices(self) -> tuple[int, ...]:
        return tuple(sorted(self.files))


def parse_condition(fov_name: str) -> str:
    """Return the biological condition portion of an observed FOV folder name."""

    if "-Hoechst" in fov_name:
        return fov_name.split("-Hoechst", 1)[0]
    return fov_name


def parse_review_condition(condition: str) -> tuple[str, str]:
    """Split an acquisition condition into pooled genotype and dose."""

    match = re.fullmatch(r"(dSPEN_(?:FL|dRRM))-(1x|p1x)", condition)
    if match is None:
        raise ValueError(f"Unsupported phase-diagram condition: {condition}")
    return match.group(1), match.group(2)


def parse_fov_token(fov_name: str) -> str:
    """Return the compact terminal FOV token used in review filenames."""

    match = re.search(r"(FOV(?:-[A-Za-z0-9_]+)*)$", fov_name)
    if match is None:
        raise ValueError(f"FOV name does not end in a safe FOV token: {fov_name}")
    return match.group(1)


def nucleus_review_stem(
    plan: PhaseFOVPlan,
    nucleus_id: int,
    rounded_mean_intensity: int,
) -> str:
    pooled_condition, dose = parse_review_condition(plan.condition)
    return (
        f"{rounded_mean_intensity}-{pooled_condition}-{dose}-"
        f"{parse_fov_token(plan.fov_name)}-Hoechst_SPEN_mask-nucleus-{nucleus_id:04d}"
    )


def discover_phase_fovs(
    dataset_root: str | Path,
    *,
    position_folder: str = "pos_0",
    glob_pattern: str = "*.tif",
) -> list[PhaseFOVPlan]:
    root = Path(dataset_root)
    if not root.is_dir():
        raise FileNotFoundError(f"Dataset root does not exist: {root}")

    plans: list[PhaseFOVPlan] = []
    for position_path in sorted(path for path in root.rglob(position_folder) if path.is_dir()):
        fov_folder = position_path.parent
        parsed: list[tuple[str, int, int, int, Path]] = []
        for path in sorted(position_path.glob(glob_pattern)):
            if not path.is_file():
                continue
            match = _ONI_SINGLE_TIME_RE.match(path.name)
            if match is None:
                continue
            parsed.append(
                (
                    match.group("prefix"),
                    int(match.group("pos_xy")),
                    int(match.group("time_index")),
                    int(match.group("z_index")),
                    path,
                )
            )
        if not parsed:
            raise FileNotFoundError(f"No ONI z-stack TIFFs found in {position_path}")

        prefixes = {item[0] for item in parsed}
        positions = {item[1] for item in parsed}
        times = {item[2] for item in parsed}
        if len(prefixes) != 1 or len(positions) != 1:
            raise ValueError(f"Expected one prefix and position in {position_path}; got {prefixes}, {positions}")
        if times != {0}:
            raise ValueError(f"Expected a single t0 acquisition in {position_path}; got time indices {sorted(times)}")

        files: dict[int, Path] = {}
        for _, _, _, z_index, path in parsed:
            if z_index in files:
                raise ValueError(f"Duplicate z{z_index} in {position_path}: {files[z_index]} and {path}")
            files[z_index] = path
        expected = tuple(range(min(files), max(files) + 1))
        if tuple(sorted(files)) != expected:
            raise ValueError(f"Non-contiguous z grid in {position_path}: {sorted(files)}")

        fov_name = fov_folder.name
        plans.append(
            PhaseFOVPlan(
                fov_folder=fov_folder,
                fov_name=fov_name,
                condition=parse_condition(fov_name),
                prefix=next(iter(prefixes)),
                pos_xy=next(iter(positions)),
                files=files,
            )
        )
    return plans


def split_dual_view(stack: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Split an ONI dual-view stack into left Hoechst and right SPEN movies."""

    arr = np.asarray(stack)
    if arr.ndim != 3:
        raise ValueError(f"Expected a frames,Y,X movie, got shape {arr.shape}")
    if arr.shape[-1] % 2:
        raise ValueError(f"Dual-view width must be even, got {arr.shape[-1]}")
    midpoint = arr.shape[-1] // 2
    return arr[..., :midpoint], arr[..., midpoint:]


def get_bbox_with_padding(
    mask: np.ndarray,
    padding: int,
    image_shape: tuple[int, int] | None = None,
) -> tuple[int, int, int, int]:
    if padding < 0:
        raise ValueError("padding must be nonnegative")
    coords = np.where(mask)
    if coords[0].size == 0:
        raise ValueError("Cannot calculate a bounding box for an empty mask")
    height, width = image_shape or mask.shape
    y0 = max(0, int(coords[0].min()) - padding)
    y1 = min(height, int(coords[0].max()) + padding + 1)
    x0 = max(0, int(coords[1].min()) - padding)
    x1 = min(width, int(coords[1].max()) + padding + 1)
    return y0, y1, x0, x1


def make_intensity_core(
    mask: np.ndarray,
    erosion_radius_fraction: float = 1.0 / 3.0,
) -> tuple[np.ndarray, float, float, float]:
    """Return an equivalent-radius eroded core and its geometry in pixels."""

    binary = np.asarray(mask, dtype=bool)
    area = int(binary.sum())
    if area == 0:
        raise ValueError("Cannot erode an empty nucleus mask")
    if not 0 < erosion_radius_fraction < 1:
        raise ValueError("erosion_radius_fraction must be between zero and one")
    equivalent_radius_px = float(np.sqrt(area / np.pi))
    erosion_distance_px = equivalent_radius_px * erosion_radius_fraction
    distances = distance_transform_edt(binary)
    core = binary & (distances > erosion_distance_px)
    return core, equivalent_radius_px, erosion_distance_px, float(np.max(distances))


def calculate_dispersion_scores(values: np.ndarray) -> tuple[float, float, float]:
    """Return Gini, coefficient of variation, and positive-median percentile dispersion."""

    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if not finite.size or np.any(finite < 0):
        raise ValueError("Dispersion scores require finite, nonnegative SPEN intensities")
    mean = float(np.mean(finite))
    if mean <= 0:
        raise ValueError("Dispersion scores require a positive mean")
    ordered = np.sort(finite)
    indices = np.arange(1, ordered.size + 1, dtype=np.float64)
    gini = float(
        2.0 * np.dot(indices, ordered) / (ordered.size * np.sum(ordered))
        - (ordered.size + 1.0) / ordered.size
    )
    coefficient_of_variation = float(np.std(finite) / mean)
    positive = finite[finite > 0]
    positive_median = float(np.median(positive))
    p10, p90 = (float(value) for value in np.percentile(finite, [10.0, 90.0]))
    robust_percentile_dispersion = float((p90 - p10) / positive_median)
    return gini, coefficient_of_variation, robust_percentile_dispersion


def preprocess_cellpose_image(image: np.ndarray, lower: float = 1.0, upper: float = 99.8) -> np.ndarray:
    arr = np.asarray(image, dtype=np.float32)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return np.zeros_like(arr, dtype=np.float32)
    lo, hi = np.percentile(finite, [lower, upper])
    if hi <= lo:
        maximum = float(np.max(finite))
        return arr / maximum if maximum > 0 else np.zeros_like(arr, dtype=np.float32)
    return np.clip((arr - lo) / (hi - lo), 0, 1).astype(np.float32)


def run_cellpose_nuclei(image: np.ndarray, model: Any, config: PhaseDiagramConfig) -> np.ndarray:
    normalized = preprocess_cellpose_image(
        image,
        config.cellpose_lower_percentile,
        config.cellpose_upper_percentile,
    )
    model_input = np.zeros((2,) + normalized.shape, dtype=np.float32)
    model_input[0] = normalized
    masks, _, _ = model.eval(
        model_input,
        batch_size=8,
        diameter=config.cellpose_diameter_px,
        flow_threshold=config.cellpose_flow_threshold,
        cellprob_threshold=config.cellpose_cellprob_threshold,
        normalize={"tile_norm_blocksize": 0},
    )
    maximum = int(np.max(masks)) if np.size(masks) else 0
    dtype = np.uint16 if maximum <= np.iinfo(np.uint16).max else np.uint32
    return np.asarray(masks, dtype=dtype)


def _sacd_params(config: PhaseDiagramConfig, *, pixel_nm: float, na: float, wavelength_nm: float) -> SACDParams:
    return SACDParams(
        pixel_nm=pixel_nm,
        wavelength_nm=wavelength_nm,
        na=na,
        mag=config.mag,
        iter1=config.iter1,
        iter2=config.iter2,
        ac_order=config.ac_order,
        subfactor=config.subfactor,
        frames_per_sacd=None,
        ifbackground=config.ifbackground,
        backgroundfactor=config.backgroundfactor,
        ifregistration=config.ifregistration,
        ifsparsedecon=config.ifsparsedecon,
        fidelity=config.fidelity,
        tcontinuity=config.tcontinuity,
        sparsity=config.sparsity,
        sparse_iterations=config.sparse_iterations,
    )


def _stage_z_summary(files: Iterable[Path]) -> tuple[list[float | None], float | None, list[str]]:
    positions: list[float | None] = []
    for path in files:
        stage = read_oni_metadata(path).get("StagePos_um")
        positions.append(float(stage[2]) if isinstance(stage, list) and len(stage) >= 3 else None)
    warnings: list[str] = []
    valid = [value for value in positions if value is not None]
    if len(valid) < 2:
        return positions, None, ["z_spacing_unavailable"]
    signed = np.diff(valid)
    nonzero = signed[np.abs(signed) > 0.01]
    if len(nonzero) and np.any(np.sign(nonzero) != np.sign(nonzero[0])):
        warnings.append("non_monotonic_stage_z")
    if np.any(np.abs(signed) <= 0.01):
        warnings.append("repeated_stage_z")
    spacing = float(np.median(np.abs(signed)))
    return positions, spacing, warnings


def _write_imagej(
    path: Path,
    image: np.ndarray,
    *,
    axes: str,
    pixel_size_um: float,
    z_spacing_um: float | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    arr = np.asarray(image)
    metadata: dict[str, Any] = {"axes": axes, "unit": "um"}
    if "Z" in axes and z_spacing_um is not None:
        metadata["spacing"] = float(z_spacing_um)
    tifffile.imwrite(
        temporary,
        arr,
        imagej=True,
        metadata=metadata,
        resolution=(1.0 / pixel_size_um, 1.0 / pixel_size_um),
    )
    os.replace(temporary, path)


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    fieldnames = list(rows[0]) if rows else []
    with temporary.open("w", newline="") as handle:
        if fieldnames:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
    os.replace(temporary, path)


def _fov_output_paths(config: PhaseDiagramConfig, plan: PhaseFOVPlan) -> dict[str, Path]:
    recon = config.output_root / "reconstructions" / plan.condition / plan.fov_name
    pooled_condition, _ = parse_review_condition(plan.condition)
    nuclei = config.output_root / "nuclei" / pooled_condition
    qc = config.output_root / "qc" / plan.condition
    stem = plan.fov_name
    return {
        "recon_dir": recon,
        "nuclei_dir": nuclei,
        "qc_dir": qc,
        "hoechst_stack": recon / f"{stem}__Hoechst-SACD-ZYX.tif",
        "hoechst_mip": recon / f"{stem}__Hoechst-SACD-MIP-YX.tif",
        "spen_stack": recon / f"{stem}__SPEN-SACD-ZYX.tif",
        "spen_mip": recon / f"{stem}__SPEN-SACD-MIP-YX.tif",
        "raw_hoechst_mip": recon / f"{stem}__Hoechst-raw-mean-MIP-YX.tif",
        "raw_spen_mip": recon / f"{stem}__SPEN-raw-mean-MIP-YX.tif",
        "nuclei_labels": recon / f"{stem}__CellposeSAM-nuclei-labels-YX.tif",
        "qc_overlay": qc / f"{stem}__CellposeSAM-overlay.png",
        "complete": recon / "_COMPLETE.json",
    }


def _export_nucleus_crops(
    plan: PhaseFOVPlan,
    config: PhaseDiagramConfig,
    nuclei_dir: Path,
    labels: np.ndarray,
    hoechst_mip: np.ndarray,
    spen_mip: np.ndarray,
    *,
    pixel_size_um: float | None = None,
    manual_review: dict[tuple[str, int], tuple[str, str]] | None = None,
    reuse_existing: bool = False,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    records: list[dict[str, Any]] = []
    excluded_records: list[dict[str, Any]] = []
    output_pixel_um = pixel_size_um or config.fallback_pixel_nm / 1000.0 / config.mag
    pooled_condition, dose = parse_review_condition(plan.condition)
    fov_token = parse_fov_token(plan.fov_name)
    nuclei_dir.mkdir(parents=True, exist_ok=True)

    for nucleus_id in (int(value) for value in np.unique(labels) if value > 0):
        nucleus_mask = labels == nucleus_id
        y0, y1, x0, x1 = get_bbox_with_padding(
            nucleus_mask,
            config.crop_padding_px,
            labels.shape,
        )
        mask_crop = nucleus_mask[y0:y1, x0:x1]
        hoechst_crop = np.asarray(hoechst_mip[y0:y1, x0:x1], dtype=np.float32)
        spen_crop = np.asarray(spen_mip[y0:y1, x0:x1], dtype=np.float32)
        full_mask_values = np.asarray(spen_mip[nucleus_mask], dtype=np.float64)
        if not full_mask_values.size or not np.all(np.isfinite(full_mask_values)):
            raise ValueError(f"Nucleus {nucleus_id} in {plan.fov_name} has no finite SPEN SACD pixels")
        core_mask, equivalent_radius_px, erosion_distance_px, max_inscribed_radius_px = make_intensity_core(
            nucleus_mask,
            config.core_erosion_radius_fraction,
        )
        manual_use, manual_note = (manual_review or {}).get((plan.fov_name, nucleus_id), ("", ""))
        key = f"{plan.fov_name}__nucleus-{nucleus_id:04d}"
        exclusion_reason = ""
        if not np.any(core_mask):
            exclusion_reason = "empty_core_after_equivalent_radius_erosion"
            core_values = np.asarray([], dtype=np.float64)
        else:
            core_values = np.asarray(spen_mip[core_mask], dtype=np.float64)
            if not np.all(np.isfinite(core_values)):
                exclusion_reason = "nonfinite_core_intensity"
            elif float(np.mean(core_values)) <= 0:
                exclusion_reason = "nonpositive_core_mean"
        if exclusion_reason:
            excluded_records.append(
                {
                    "nucleus_key": key,
                    "condition": plan.condition,
                    "pooled_condition": pooled_condition,
                    "dose": dose,
                    "fov": plan.fov_name,
                    "fov_token": fov_token,
                    "nucleus_id": nucleus_id,
                    "exclusion_reason": exclusion_reason,
                    "area_sacd_px": int(nucleus_mask.sum()),
                    "equivalent_radius_sacd_px": equivalent_radius_px,
                    "erosion_distance_sacd_px": erosion_distance_px,
                    "max_inscribed_radius_sacd_px": max_inscribed_radius_px,
                    "core_area_sacd_px": int(core_mask.sum()),
                    "manual_use": manual_use,
                    "manual_note": manual_note,
                }
            )
            continue

        mean_intensity = float(np.mean(core_values))
        max_intensity = float(np.max(core_values))
        gini, coefficient_of_variation, robust_percentile_dispersion = calculate_dispersion_scores(core_values)
        rounded_mean = int(np.rint(mean_intensity))
        stem = nucleus_review_stem(plan, nucleus_id, rounded_mean)
        review_tif = nuclei_dir / f"{stem}.tif"
        review_png = nuclei_dir / f"{stem}.png"
        review_stack = np.stack(
            (hoechst_crop, spen_crop, mask_crop.astype(np.float32)),
            axis=0,
        ).astype(np.float32, copy=False)
        if not (reuse_existing and review_tif.is_file() and review_png.is_file()):
            _write_imagej(review_tif, review_stack, axes="CYX", pixel_size_um=output_pixel_um)
            _make_nucleus_review_png(
                review_png,
                spen_crop,
                mask_crop,
                core_mask[y0:y1, x0:x1],
                pixel_size_um=output_pixel_um,
                scale_bar_um=config.review_scale_bar_um,
                lower_percentile=config.review_lower_percentile,
                upper_percentile=config.review_upper_percentile,
            )

        ys, xs = np.where(nucleus_mask)
        touches_border = bool(
            np.any(ys == 0)
            or np.any(xs == 0)
            or np.any(ys == labels.shape[0] - 1)
            or np.any(xs == labels.shape[1] - 1)
        )
        records.append(
            {
                "nucleus_key": key,
                "condition": plan.condition,
                "pooled_condition": pooled_condition,
                "dose": dose,
                "fov": plan.fov_name,
                "fov_token": fov_token,
                "nucleus_id": nucleus_id,
                "spen_sacd_core_mean": mean_intensity,
                "spen_sacd_core_mean_rounded": rounded_mean,
                "spen_sacd_core_max": max_intensity,
                "spen_sacd_core_gini": gini,
                "spen_sacd_core_coefficient_of_variation": coefficient_of_variation,
                "spen_sacd_core_robust_percentile_dispersion": robust_percentile_dispersion,
                "spen_sacd_full_mask_mean": float(np.mean(full_mask_values)),
                "area_sacd_px": int(nucleus_mask.sum()),
                "equivalent_radius_sacd_px": equivalent_radius_px,
                "erosion_distance_sacd_px": erosion_distance_px,
                "max_inscribed_radius_sacd_px": max_inscribed_radius_px,
                "core_area_sacd_px": int(core_mask.sum()),
                "core_area_fraction": float(core_mask.sum() / nucleus_mask.sum()),
                "centroid_y_sacd_px": float(np.mean(ys)),
                "centroid_x_sacd_px": float(np.mean(xs)),
                "bbox_y0_sacd_px": y0,
                "bbox_y1_sacd_px": y1,
                "bbox_x0_sacd_px": x0,
                "bbox_x1_sacd_px": x1,
                "touches_border": touches_border,
                "review_tif": str(review_tif),
                "review_png": str(review_png),
                "manual_use": manual_use,
                "manual_note": manual_note,
            }
        )
    return records, excluded_records


def _make_nucleus_review_png(
    path: Path,
    spen_crop: np.ndarray,
    mask_crop: np.ndarray,
    core_crop: np.ndarray,
    *,
    pixel_size_um: float,
    scale_bar_um: float,
    lower_percentile: float,
    upper_percentile: float,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    image = np.asarray(spen_crop, dtype=np.float32)
    mask = np.asarray(mask_crop, dtype=bool)
    core = np.asarray(core_crop, dtype=bool)
    finite = image[mask & np.isfinite(image)]
    if finite.size == 0:
        raise ValueError("Cannot render a nucleus review PNG without finite masked SPEN pixels")
    vmin, vmax = (float(value) for value in np.percentile(finite, [lower_percentile, upper_percentile]))
    if vmax <= vmin:
        vmin = float(np.min(finite))
        vmax = float(np.max(finite))
    if vmax <= vmin:
        vmax = vmin + 1.0

    fig, ax = plt.subplots(figsize=(5.4, 5.1))
    shown = ax.imshow(image, cmap="magma", vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.contour(mask.astype(np.uint8), levels=[0.5], colors=["cyan"], linewidths=0.9)
    ax.contour(core.astype(np.uint8), levels=[0.5], colors=["lime"], linewidths=1.1, linestyles=":")
    bar_px = scale_bar_um / pixel_size_um
    x0 = max(4.0, image.shape[1] * 0.08)
    y0 = image.shape[0] * 0.90
    ax.plot([x0, x0 + bar_px], [y0, y0], color="white", linewidth=3, solid_capstyle="butt")
    ax.text(
        x0 + bar_px / 2,
        y0 - max(2.0, image.shape[0] * 0.025),
        f"{scale_bar_um:g} µm",
        color="white",
        ha="center",
        va="bottom",
        fontsize=8,
    )
    ax.set_axis_off()
    colorbar = fig.colorbar(shown, ax=ax, fraction=0.046, pad=0.025)
    colorbar.set_label("SPEN SACD intensity (a.u.)")
    fig.legend(
        handles=[
            Line2D([0], [0], color="cyan", linewidth=1.2, linestyle="-", label="CellposeSAM boundary"),
            Line2D([0], [0], color="lime", linewidth=1.4, linestyle=":", label="Intensity core (r_eq/3 erosion)"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.01),
        ncol=2,
        frameon=False,
        fontsize=8,
    )
    fig.subplots_adjust(bottom=0.13)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.png")
    fig.savefig(temporary, dpi=160, facecolor="white")
    plt.close(fig)
    os.replace(temporary, path)


def _make_qc_overlay(path: Path, image: np.ndarray, labels: np.ndarray) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from skimage.segmentation import find_boundaries

    normalized = preprocess_cellpose_image(image)
    rgb = np.repeat(normalized[..., None], 3, axis=2)
    boundary = find_boundaries(labels, mode="outer")
    rgb[boundary] = (1.0, 0.0, 0.0)
    fig, ax = plt.subplots(figsize=(8, 12))
    ax.imshow(rgb)
    for nucleus_id in (int(value) for value in np.unique(labels) if value > 0):
        ys, xs = np.where(labels == nucleus_id)
        ax.text(float(np.mean(xs)), float(np.mean(ys)), str(nucleus_id), color="yellow", fontsize=7, ha="center", va="center")
    ax.set_title(f"CellposeSAM nuclei: {int(np.max(labels)) if labels.size else 0}")
    ax.axis("off")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.png")
    fig.savefig(temporary, dpi=150, bbox_inches="tight")
    plt.close(fig)
    os.replace(temporary, path)


def make_dispersion_phase_diagram(records: list[dict[str, Any]], path: Path) -> None:
    """Compare three core-intensity dispersion scores against core mean."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not records:
        raise ValueError("Cannot plot a phase diagram without nucleus records")
    required = {
        "pooled_condition",
        "spen_sacd_core_mean",
        "spen_sacd_core_gini",
        "spen_sacd_core_coefficient_of_variation",
        "spen_sacd_core_robust_percentile_dispersion",
    }
    missing = required - set(records[0])
    if missing:
        raise ValueError(f"Nucleus records lack phase-diagram columns: {sorted(missing)}")

    palette = {"dSPEN_FL": "#2468b4", "dSPEN_dRRM": "#e4572e"}
    score_columns = (
        ("spen_sacd_core_gini", "Gini coefficient", "Gini coefficient"),
        (
            "spen_sacd_core_coefficient_of_variation",
            "Coefficient of variation",
            "Standard deviation / mean",
        ),
        (
            "spen_sacd_core_robust_percentile_dispersion",
            "Robust percentile dispersion",
            "(P90−P10) / median(positive pixels)",
        ),
    )
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 4.9), sharex=True, constrained_layout=True)
    for condition in ("dSPEN_FL", "dSPEN_dRRM"):
        selected = [record for record in records if record["pooled_condition"] == condition]
        x = np.asarray([float(record["spen_sacd_core_mean"]) for record in selected])
        for ax, (column, _, _) in zip(axes, score_columns, strict=True):
            y = np.asarray([float(record[column]) for record in selected])
            ax.scatter(
                x,
                y,
                s=13,
                alpha=0.55,
                linewidths=0,
                color=palette[condition],
                label=f"{condition} (n={len(selected)})",
                rasterized=True,
            )
    for ax, (_, title, ylabel) in zip(axes, score_columns, strict=True):
        ax.set_title(title)
        ax.set_xscale("log")
        ax.set_yscale("linear")
        ax.set_xlabel("Mean SPEN SACD intensity in eroded core (a.u.)")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.2, linewidth=0.6)
    axes[0].legend(frameon=False)
    fig.suptitle("SPEN phase-diagram dispersion comparison: retained nuclei")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.png")
    fig.savefig(temporary, dpi=180, facecolor="white")
    plt.close(fig)
    os.replace(temporary, path)


def _validate_fov_outputs(
    paths: dict[str, Path],
    nucleus_records: list[dict[str, Any]],
    excluded_records: list[dict[str, Any]],
) -> None:
    required = [
        "hoechst_stack",
        "hoechst_mip",
        "spen_stack",
        "spen_mip",
        "raw_hoechst_mip",
        "raw_spen_mip",
        "nuclei_labels",
        "qc_overlay",
    ]
    missing = [name for name in required if not paths[name].exists()]
    if missing:
        raise FileNotFoundError(f"Missing required FOV outputs: {missing}")
    hoechst_stack = tifffile.imread(paths["hoechst_stack"])
    spen_stack = tifffile.imread(paths["spen_stack"])
    if not np.array_equal(tifffile.imread(paths["hoechst_mip"]), np.max(hoechst_stack, axis=0)):
        raise AssertionError("Saved Hoechst MIP does not equal max(saved stack, axis=Z)")
    if not np.array_equal(tifffile.imread(paths["spen_mip"]), np.max(spen_stack, axis=0)):
        raise AssertionError("Saved SPEN MIP does not equal max(saved stack, axis=Z)")
    labels = tifffile.imread(paths["nuclei_labels"])
    if len(nucleus_records) + len(excluded_records) != len([value for value in np.unique(labels) if value > 0]):
        raise AssertionError("Retained plus excluded nucleus counts do not match the label image")
    for record in nucleus_records:
        for column in ("review_tif", "review_png"):
            if not Path(record[column]).exists():
                raise FileNotFoundError(record[column])
        review = tifffile.imread(record["review_tif"])
        if review.ndim != 3 or review.shape[0] != 3 or review.dtype != np.float32:
            raise AssertionError(f"Invalid review TIFF shape/dtype: {record['review_tif']} {review.shape} {review.dtype}")
        if not np.all(np.isin(np.unique(review[2]), (0.0, 1.0))):
            raise AssertionError(f"Review TIFF mask is not binary: {record['review_tif']}")


def _relative_records(records: list[dict[str, Any]], staging: Path, final: Path) -> list[dict[str, Any]]:
    path_columns = ("review_tif", "review_png")
    converted: list[dict[str, Any]] = []
    for record in records:
        item = dict(record)
        for column in path_columns:
            item[column] = str(final / Path(item[column]).relative_to(staging))
        converted.append(item)
    return converted


def _remove_superseded_nucleus_files(
    records: Iterable[dict[str, Any]],
    output_root: Path,
    *,
    keep: set[Path] | None = None,
) -> None:
    nuclei_root = (output_root / "nuclei").resolve()
    keep_resolved = {path.resolve() for path in (keep or set())}
    path_columns = (
        "review_tif",
        "review_png",
        "hoechst_sacd_crop",
        "spen_sacd_crop",
        "raw_hoechst_crop",
        "raw_spen_crop",
        "mask_crop",
    )
    for record in records:
        for column in path_columns:
            value = record.get(column)
            if not value:
                continue
            path = Path(value)
            resolved = path.resolve()
            if not resolved.is_relative_to(nuclei_root):
                raise ValueError(f"Refusing to remove nucleus output outside {nuclei_root}: {path}")
            if resolved not in keep_resolved and path.is_file():
                path.unlink()


def process_phase_fov(plan: PhaseFOVPlan, config: PhaseDiagramConfig, model: Any) -> dict[str, Any]:
    final_paths = _fov_output_paths(config, plan)
    if final_paths["complete"].exists() and not config.overwrite:
        completion = json.loads(final_paths["complete"].read_text())
        completion["status"] = "skipped_existing"
        return completion
    previous_completion = (
        json.loads(final_paths["complete"].read_text())
        if final_paths["complete"].exists()
        else None
    )

    attempt = config.output_root / ".work" / f"{plan.fov_name}-{uuid.uuid4().hex[:8]}"
    stage_recon = attempt / "reconstruction"
    stage_nuclei = attempt / "nuclei"
    stage_qc = attempt / "qc" / final_paths["qc_overlay"].name
    stage_paths = {
        **final_paths,
        "recon_dir": stage_recon,
        "nuclei_dir": stage_nuclei,
        "qc_dir": stage_qc.parent,
        "hoechst_stack": stage_recon / final_paths["hoechst_stack"].name,
        "hoechst_mip": stage_recon / final_paths["hoechst_mip"].name,
        "spen_stack": stage_recon / final_paths["spen_stack"].name,
        "spen_mip": stage_recon / final_paths["spen_mip"].name,
        "raw_hoechst_mip": stage_recon / final_paths["raw_hoechst_mip"].name,
        "raw_spen_mip": stage_recon / final_paths["raw_spen_mip"].name,
        "nuclei_labels": stage_recon / final_paths["nuclei_labels"].name,
        "qc_overlay": stage_qc,
    }
    stage_recon.mkdir(parents=True, exist_ok=True)

    started = perf_counter()
    first_metadata = read_oni_metadata(plan.files[plan.z_indices[0]])
    pixel_nm = float(first_metadata.get("PixelSize_um", config.fallback_pixel_nm / 1000.0)) * 1000.0
    na = float(first_metadata.get("Objective_NA", config.fallback_na))
    output_pixel_um = pixel_nm / 1000.0 / config.mag
    raw_pixel_um = pixel_nm / 1000.0
    stage_positions, z_spacing_um, metadata_warnings = _stage_z_summary(plan.files[z] for z in plan.z_indices)

    hoechst_params = _sacd_params(config, pixel_nm=pixel_nm, na=na, wavelength_nm=config.hoechst_wavelength_nm)
    spen_params = _sacd_params(config, pixel_nm=pixel_nm, na=na, wavelength_nm=config.spen_wavelength_nm)
    hoechst_z: list[np.ndarray] = []
    spen_z: list[np.ndarray] = []
    raw_hoechst_mip: np.ndarray | None = None
    raw_spen_mip: np.ndarray | None = None
    input_shape: tuple[int, ...] | None = None

    for ordinal, z_index in enumerate(plan.z_indices, start=1):
        input_file = plan.files[z_index]
        raw = tifffile.imread(input_file)
        input_shape = tuple(raw.shape)
        left, right = split_dual_view(raw)
        left_mean = np.mean(left, axis=0, dtype=np.float32)
        right_mean = np.mean(right, axis=0, dtype=np.float32)
        raw_hoechst_mip = left_mean if raw_hoechst_mip is None else np.maximum(raw_hoechst_mip, left_mean)
        raw_spen_mip = right_mean if raw_spen_mip is None else np.maximum(raw_spen_mip, right_mean)
        print(f"  {plan.fov_name}: z{z_index} ({ordinal}/{len(plan.z_indices)}) Hoechst", flush=True)
        hoechst_z.append(reconstruct(left, hoechst_params))
        print(f"  {plan.fov_name}: z{z_index} ({ordinal}/{len(plan.z_indices)}) SPEN", flush=True)
        spen_z.append(reconstruct(right, spen_params))

    assert raw_hoechst_mip is not None and raw_spen_mip is not None
    hoechst_stack = np.stack(hoechst_z).astype(np.float32, copy=False)
    spen_stack = np.stack(spen_z).astype(np.float32, copy=False)
    hoechst_mip = np.max(hoechst_stack, axis=0).astype(np.float32, copy=False)
    spen_mip = np.max(spen_stack, axis=0).astype(np.float32, copy=False)

    _write_imagej(stage_paths["hoechst_stack"], hoechst_stack, axes="ZYX", pixel_size_um=output_pixel_um, z_spacing_um=z_spacing_um)
    _write_imagej(stage_paths["hoechst_mip"], hoechst_mip, axes="YX", pixel_size_um=output_pixel_um)
    _write_imagej(stage_paths["spen_stack"], spen_stack, axes="ZYX", pixel_size_um=output_pixel_um, z_spacing_um=z_spacing_um)
    _write_imagej(stage_paths["spen_mip"], spen_mip, axes="YX", pixel_size_um=output_pixel_um)
    _write_imagej(stage_paths["raw_hoechst_mip"], raw_hoechst_mip, axes="YX", pixel_size_um=raw_pixel_um)
    _write_imagej(stage_paths["raw_spen_mip"], raw_spen_mip, axes="YX", pixel_size_um=raw_pixel_um)

    print(f"  {plan.fov_name}: CellposeSAM", flush=True)
    labels = run_cellpose_nuclei(hoechst_mip, model, config)
    _write_imagej(stage_paths["nuclei_labels"], labels, axes="YX", pixel_size_um=output_pixel_um)
    nucleus_records, excluded_records = _export_nucleus_crops(
        plan,
        config,
        stage_nuclei,
        labels,
        hoechst_mip,
        spen_mip,
        pixel_size_um=output_pixel_um,
    )
    _make_qc_overlay(stage_paths["qc_overlay"], hoechst_mip, labels)
    _validate_fov_outputs(stage_paths, nucleus_records, excluded_records)

    final_paths["recon_dir"].mkdir(parents=True, exist_ok=True)
    final_paths["qc_dir"].mkdir(parents=True, exist_ok=True)
    for name in (
        "hoechst_stack",
        "hoechst_mip",
        "spen_stack",
        "spen_mip",
        "raw_hoechst_mip",
        "raw_spen_mip",
        "nuclei_labels",
    ):
        os.replace(stage_paths[name], final_paths[name])
    os.replace(stage_paths["qc_overlay"], final_paths["qc_overlay"])
    final_paths["nuclei_dir"].mkdir(parents=True, exist_ok=True)
    for staged_file in stage_nuclei.iterdir():
        os.replace(staged_file, final_paths["nuclei_dir"] / staged_file.name)
    final_records = _relative_records(nucleus_records, stage_nuclei, final_paths["nuclei_dir"])
    if previous_completion is not None:
        keep = {
            Path(record[column])
            for record in final_records
            for column in ("review_tif", "review_png")
        }
        _remove_superseded_nucleus_files(
            previous_completion.get("nucleus_records", []),
            config.output_root,
            keep=keep,
        )
        legacy_dir = config.output_root / "nuclei" / plan.condition / plan.fov_name
        if legacy_dir.is_dir():
            shutil.rmtree(legacy_dir)

    fov_record = {
        "condition": plan.condition,
        "fov": plan.fov_name,
        "source_folder": str(plan.fov_folder),
        "source_files": [str(plan.files[z]) for z in plan.z_indices],
        "z_indices": list(plan.z_indices),
        "stage_z_um": stage_positions,
        "z_spacing_um": z_spacing_um,
        "metadata_warnings": metadata_warnings,
        "pixel_nm": pixel_nm,
        "output_pixel_nm": pixel_nm / config.mag,
        "na": na,
        "input_shape": list(input_shape or ()),
        "sacd_stack_shape": list(hoechst_stack.shape),
        "n_nuclei": len(final_records),
        "n_segmented_nuclei": len(final_records) + len(excluded_records),
        "n_excluded_nuclei": len(excluded_records),
        "runtime_s": perf_counter() - started,
        "status": "written",
        "outputs": {name: str(final_paths[name]) for name in final_paths if name not in {"recon_dir", "nuclei_dir", "qc_dir", "complete"}},
    }
    completion = {
        "fov_record": fov_record,
        "nucleus_records": final_records,
        "excluded_nucleus_records": excluded_records,
        "status": "written",
    }
    _write_json(final_paths["complete"], completion)
    shutil.rmtree(attempt, ignore_errors=True)
    return completion


def _config_json(config: PhaseDiagramConfig) -> dict[str, Any]:
    value = asdict(config)
    value["dataset_root"] = str(config.dataset_root)
    value["output_root"] = str(config.output_root)
    return value


def _load_cellpose_model(config: PhaseDiagramConfig) -> Any:
    try:
        from cellpose import core, models
    except ImportError as exc:
        raise RuntimeError("Cellpose is required; run with the conda-smlm kernel/environment") from exc
    gpu = False if config.cellpose_device == "cpu" else bool(core.use_gpu())
    print(f"Loading CellposeSAM (gpu={gpu})", flush=True)
    return models.CellposeModel(gpu=gpu)


def consolidate_manifests(config: PhaseDiagramConfig, completions: list[dict[str, Any]], failures: list[dict[str, Any]]) -> None:
    manifest_dir = config.output_root / "manifests"
    fov_rows = [item["fov_record"] for item in completions]
    flat_fov_rows: list[dict[str, Any]] = []
    for row in fov_rows:
        flat_fov_rows.append(
            {
                **{key: value for key, value in row.items() if key not in {"source_files", "stage_z_um", "metadata_warnings", "input_shape", "sacd_stack_shape", "outputs"}},
                "z_indices": json.dumps(row["z_indices"]),
                "stage_z_um": json.dumps(row["stage_z_um"]),
                "metadata_warnings": json.dumps(row["metadata_warnings"]),
                "input_shape": json.dumps(row["input_shape"]),
                "sacd_stack_shape": json.dumps(row["sacd_stack_shape"]),
                "source_files": json.dumps(row["source_files"]),
                "outputs": json.dumps(row["outputs"], sort_keys=True),
            }
        )
    nucleus_rows = [record for item in completions for record in item["nucleus_records"]]
    excluded_rows = [record for item in completions for record in item.get("excluded_nucleus_records", [])]
    _write_csv(manifest_dir / "fov_manifest.csv", flat_fov_rows)
    _write_csv(manifest_dir / "nucleus_manifest.csv", nucleus_rows)
    _write_csv(manifest_dir / "excluded_nuclei.csv", excluded_rows)
    if nucleus_rows and all("spen_sacd_core_gini" in row for row in nucleus_rows):
        make_dispersion_phase_diagram(
            nucleus_rows,
            config.output_root / "qc" / "phase_diagram-SPEN_core_mean-vs-dispersion.png",
        )
    _write_json(manifest_dir / "failures.json", failures)
    summary = {
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "dataset_root": str(config.dataset_root),
        "output_root": str(config.output_root),
        "completed_fovs": len(completions),
        "failed_fovs": len(failures),
        "total_nuclei": len(nucleus_rows),
        "total_segmented_nuclei": len(nucleus_rows) + len(excluded_rows),
        "total_retained_nuclei": len(nucleus_rows),
        "total_excluded_nuclei": len(excluded_rows),
        "total_source_movies": sum(len(item["fov_record"]["source_files"]) for item in completions),
        "conditions": {
            condition: sum(1 for row in fov_rows if row["condition"] == condition)
            for condition in sorted({row["condition"] for row in fov_rows})
        },
    }
    _write_json(manifest_dir / "run_summary.json", summary)


def _load_manual_review(path: Path) -> dict[tuple[str, int], tuple[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        return {
            (row["fov"], int(row["nucleus_id"])): (
                row.get("manual_use", ""),
                row.get("manual_note", ""),
            )
            for row in csv.DictReader(handle)
        }


def repackage_existing_nuclei(
    config: PhaseDiagramConfig,
    *,
    expected_nuclei: int | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Rebuild review files from completed SACD MIPs and Cellpose labels only."""

    completion_paths = sorted((config.output_root / "reconstructions").rglob("_COMPLETE.json"))
    if not completion_paths:
        raise FileNotFoundError(f"No completed FOV checkpoints found under {config.output_root}")
    original_text = {path: path.read_text() for path in completion_paths}
    original_completions = [json.loads(original_text[path]) for path in completion_paths]
    manual_review = _load_manual_review(config.output_root / "manifests" / "nucleus_manifest.csv")
    failures_path = config.output_root / "manifests" / "failures.json"
    failures = json.loads(failures_path.read_text()) if failures_path.exists() else []

    work_root = config.output_root / ".work"
    resumable_attempts = sorted(
        path for path in work_root.glob("nuclei-repackage-*")
        if (path / "nuclei").is_dir()
    )
    attempt = (
        resumable_attempts[-1]
        if resumable_attempts
        else work_root / f"nuclei-repackage-{uuid.uuid4().hex[:8]}"
    )
    if resumable_attempts:
        print(f"Resuming staged nucleus repackaging from {attempt}", flush=True)
    stage_nuclei = attempt / "nuclei"
    stage_plot = attempt / "phase_diagram-SPEN_core_mean-vs-dispersion.png"
    staged_records: list[dict[str, Any]] = []
    staged_excluded: list[dict[str, Any]] = []
    updated_completions: list[dict[str, Any]] = []
    updates_by_path: list[tuple[Path, dict[str, Any]]] = []
    seen_names: set[tuple[str, str]] = set()

    for index, (completion_path, completion) in enumerate(
        zip(completion_paths, original_completions, strict=True),
        start=1,
    ):
        fov_record = completion["fov_record"]
        condition = str(fov_record["condition"])
        fov_name = str(fov_record["fov"])
        plan = PhaseFOVPlan(
            fov_folder=Path(fov_record["source_folder"]),
            fov_name=fov_name,
            condition=condition,
            prefix="",
            pos_xy=0,
            files={},
        )
        outputs = fov_record["outputs"]
        hoechst_mip = tifffile.imread(outputs["hoechst_mip"])
        spen_mip = tifffile.imread(outputs["spen_mip"])
        labels = tifffile.imread(outputs["nuclei_labels"])
        if hoechst_mip.shape != spen_mip.shape or hoechst_mip.shape != labels.shape:
            raise AssertionError(f"MIP/label shape mismatch in {fov_name}")
        pooled_condition, _ = parse_review_condition(condition)
        stage_condition = stage_nuclei / pooled_condition
        pixel_size_um = float(fov_record["output_pixel_nm"]) / 1000.0
        records, excluded = _export_nucleus_crops(
            plan,
            config,
            stage_condition,
            labels,
            hoechst_mip,
            spen_mip,
            pixel_size_um=pixel_size_um,
            manual_review=manual_review,
            reuse_existing=False,
        )
        expected_for_fov = len([value for value in np.unique(labels) if value > 0])
        if len(records) + len(excluded) != expected_for_fov:
            raise AssertionError(f"Review export count mismatch in {fov_name}")
        for record in records:
            for column in ("review_tif", "review_png"):
                key = (column, Path(record[column]).name)
                if key in seen_names:
                    raise AssertionError(f"Duplicate review filename: {key[1]}")
                seen_names.add(key)
        final_condition = config.output_root / "nuclei" / pooled_condition
        final_records = _relative_records(records, stage_condition, final_condition)
        updated = json.loads(json.dumps(completion))
        updated["nucleus_records"] = final_records
        updated["excluded_nucleus_records"] = excluded
        updated["status"] = "repackaged_existing"
        updated["fov_record"]["n_nuclei"] = len(final_records)
        updated["fov_record"]["n_segmented_nuclei"] = len(final_records) + len(excluded)
        updated["fov_record"]["n_excluded_nuclei"] = len(excluded)
        updated_completions.append(updated)
        updates_by_path.append((completion_path, updated))
        staged_records.extend(records)
        staged_excluded.extend(excluded)
        print(
            f"Repackaged FOV {index}/{len(completion_paths)}: {fov_name} "
            f"({len(records)} retained, {len(excluded)} excluded)",
            flush=True,
        )

    if expected_nuclei is not None and len(staged_records) + len(staged_excluded) != expected_nuclei:
        raise AssertionError(
            f"Expected {expected_nuclei} segmented nuclei, staged "
            f"{len(staged_records)} retained plus {len(staged_excluded)} excluded"
        )
    expected_conditions = {"dSPEN_FL", "dSPEN_dRRM"}
    observed_conditions = {path.name for path in stage_nuclei.iterdir() if path.is_dir()}
    if observed_conditions != expected_conditions:
        raise AssertionError(f"Expected pooled folders {sorted(expected_conditions)}, got {sorted(observed_conditions)}")
    tif_count = len(list(stage_nuclei.rglob("*.tif")))
    png_count = len(list(stage_nuclei.rglob("*.png")))
    if tif_count != len(staged_records) or png_count != len(staged_records):
        raise AssertionError(
            f"Staged review count mismatch: records={len(staged_records)}, TIFF={tif_count}, PNG={png_count}"
        )
    for record in staged_records:
        review = tifffile.imread(record["review_tif"])
        if review.ndim != 3 or review.shape[0] != 3 or review.dtype != np.float32:
            raise AssertionError(f"Invalid staged review TIFF: {record['review_tif']}")
        if not np.all(np.isin(np.unique(review[2]), (0.0, 1.0))):
            raise AssertionError(f"Nonbinary staged mask: {record['review_tif']}")

    final_records = [record for completion in updated_completions for record in completion["nucleus_records"]]
    make_dispersion_phase_diagram(final_records, stage_plot)
    if not stage_plot.is_file() or stage_plot.stat().st_size == 0:
        raise AssertionError("The staged rough phase diagram was not written")

    final_nuclei = config.output_root / "nuclei"
    backup_nuclei = config.output_root / f"nuclei.legacy-{datetime.now().strftime('%Y%m%d-%H%M%S')}"
    swapped = False
    try:
        if final_nuclei.exists():
            os.replace(final_nuclei, backup_nuclei)
        os.replace(stage_nuclei, final_nuclei)
        swapped = True
        for completion_path, updated in updates_by_path:
            _write_json(completion_path, updated)
        consolidate_manifests(config, updated_completions, failures)
        if len(list(final_nuclei.rglob("*.tif"))) != len(final_records):
            raise AssertionError("Final nucleus TIFF count changed after the atomic swap")
        if len(list(final_nuclei.rglob("*.png"))) != len(final_records):
            raise AssertionError("Final nucleus PNG count changed after the atomic swap")
    except Exception:
        for completion_path, text in original_text.items():
            temporary = completion_path.with_name(completion_path.name + ".rollback.tmp")
            temporary.write_text(text)
            os.replace(temporary, completion_path)
        if swapped and final_nuclei.exists():
            failed_nuclei = attempt / "failed_nuclei"
            os.replace(final_nuclei, failed_nuclei)
        if backup_nuclei.exists():
            os.replace(backup_nuclei, final_nuclei)
        consolidate_manifests(config, original_completions, failures)
        raise
    else:
        old_plot = config.output_root / "qc" / "rough_phase_diagram-SPEN_SACD-mean-vs-max.png"
        if old_plot.is_file():
            old_plot.unlink()
        if backup_nuclei.exists():
            shutil.rmtree(backup_nuclei, ignore_errors=True)
        shutil.rmtree(attempt, ignore_errors=True)
    return updated_completions, failures


def rebuild_manifests_from_completions(config: PhaseDiagramConfig) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Rebuild dataset-wide manifests from independently checkpointed FOVs."""

    completions: list[dict[str, Any]] = []
    for path in sorted((config.output_root / "reconstructions").rglob("_COMPLETE.json")):
        completions.append(json.loads(path.read_text()))
    failures_path = config.output_root / "manifests" / "failures.json"
    failures = json.loads(failures_path.read_text()) if failures_path.exists() else []
    consolidate_manifests(config, completions, failures)
    return completions, failures


def run_phase_dataset(
    config: PhaseDiagramConfig,
    *,
    fov_names: Iterable[str] | None = None,
    stop_on_error: bool = False,
    model: Any | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    config = PhaseDiagramConfig(**{**asdict(config), "dataset_root": Path(config.dataset_root), "output_root": Path(config.output_root)})
    plans = discover_phase_fovs(
        config.dataset_root,
        position_folder=config.position_folder,
        glob_pattern=config.glob_pattern,
    )
    selected = set(fov_names or ())
    if selected:
        plans = [plan for plan in plans if plan.fov_name in selected]
        missing = selected - {plan.fov_name for plan in plans}
        if missing:
            raise ValueError(f"Unknown requested FOV names: {sorted(missing)}")
    if not plans:
        raise ValueError("No FOVs selected")

    config.output_root.mkdir(parents=True, exist_ok=True)
    _write_json(config.output_root / "manifests" / "run_config.json", _config_json(config))
    pending = [plan for plan in plans if config.overwrite or not _fov_output_paths(config, plan)["complete"].exists()]
    if pending and model is None:
        model = _load_cellpose_model(config)

    completions: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    total_started = perf_counter()
    for index, plan in enumerate(plans, start=1):
        print(f"\nFOV {index}/{len(plans)}: {plan.fov_name}", flush=True)
        try:
            completion = process_phase_fov(plan, config, model)
            completions.append(completion)
            consolidate_manifests(config, completions, failures)
            elapsed = perf_counter() - total_started
            print(
                f"Completed {plan.fov_name}: status={completion['status']}, "
                f"nuclei={len(completion['nucleus_records'])}, elapsed={elapsed / 60:.1f} min",
                flush=True,
            )
        except Exception as exc:
            failure = {"fov": plan.fov_name, "condition": plan.condition, "error_type": type(exc).__name__, "error": str(exc)}
            failures.append(failure)
            consolidate_manifests(config, completions, failures)
            print(f"FAILED {plan.fov_name}: {type(exc).__name__}: {exc}", flush=True)
            if stop_on_error:
                raise
    consolidate_manifests(config, completions, failures)
    return completions, failures


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the Hoechst-SPEN SACDlive phase-diagram preparation pipeline")
    parser.add_argument("--dataset-root", type=Path, default=DEFAULT_DATASET_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--fov", action="append", dest="fov_names", help="Process only this exact FOV name; repeat as needed")
    parser.add_argument("--start-index", type=int, help="Process FOVs from this 1-based sorted index")
    parser.add_argument("--end-index", type=int, help="Process FOVs through this inclusive 1-based sorted index")
    parser.add_argument("--cellpose-device", choices=("auto", "cpu"), default="auto")
    parser.add_argument("--consolidate-only", action="store_true", help="Rebuild manifests from FOV completion records and exit")
    parser.add_argument(
        "--repackage-nuclei-only",
        action="store_true",
        help="Rebuild nucleus review TIFF/PNG files from completed SACD and Cellpose outputs",
    )
    parser.add_argument("--expected-nuclei", type=int, help="Require this many nuclei during repackaging")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--stop-on-error", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    config = PhaseDiagramConfig(
        dataset_root=args.dataset_root,
        output_root=args.output_root,
        overwrite=args.overwrite,
        cellpose_device=args.cellpose_device,
    )
    if args.consolidate_only:
        completions, failures = rebuild_manifests_from_completions(config)
        print(f"Consolidated: {len(completions)} completed, {len(failures)} failed", flush=True)
        return 1 if failures else 0
    if args.repackage_nuclei_only:
        completions, failures = repackage_existing_nuclei(
            config,
            expected_nuclei=args.expected_nuclei,
        )
        total_nuclei = sum(len(item["nucleus_records"]) for item in completions)
        print(
            f"Repackaged: {len(completions)} FOVs, {total_nuclei} nuclei, {len(failures)} failures",
            flush=True,
        )
        return 1 if failures else 0
    fov_names = args.fov_names
    if args.start_index is not None or args.end_index is not None:
        if fov_names:
            raise ValueError("Use either --fov or --start-index/--end-index, not both")
        plans = discover_phase_fovs(args.dataset_root)
        start = 1 if args.start_index is None else args.start_index
        end = len(plans) if args.end_index is None else args.end_index
        if start < 1 or end < start or end > len(plans):
            raise ValueError(f"Invalid FOV index range {start}..{end}; dataset has {len(plans)} FOVs")
        fov_names = [plan.fov_name for plan in plans[start - 1 : end]]
    completions, failures = run_phase_dataset(
        config,
        fov_names=fov_names,
        stop_on_error=args.stop_on_error,
    )
    print(f"Finished: {len(completions)} completed, {len(failures)} failed", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
