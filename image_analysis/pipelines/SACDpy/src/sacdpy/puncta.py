"""2D SPEN-puncta analysis for SACD phase-diagram data.

Spotiflow is used only to detect seed coordinates. Final punctum masks are
constructed on unnormalized SACD values in the configured intensity space with deterministic grouping,
bounded local growth, and a fixed-footprint fallback. Logarithms are confined
to review rendering and the phase diagram's x-axis.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import os
import shutil
import uuid
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import tifffile
from scipy import ndimage
from scipy.spatial import cKDTree

from .progress import report_progress, progress_items, progress_message, stage
from skimage import filters
from skimage.feature import blob_log
from skimage.measure import label as connected_components
from skimage.segmentation import relabel_sequential, watershed


PUNCTATOOLS_REFERENCE = {
    "package": "stjude/punctatools",
    "version": "0.2.0",
    "commit": "dfc2c33e7e0b18dfd5aa859ae589ad2ed078cf1c",
    "algorithm": "2D LoG centers, fixed background filtering, thresholding, watershed",
}

SPOTIFLOW_GENERAL_WEIGHT_SHA256 = (
    "1c3575464d621924b27f4deb66495b807f175a0ccd995d3533403f00daf806f2"
)
DEFAULT_SPOTIFLOW_CACHE = Path(
    "/Users/gmgao/GGscripts/GG-general-GuttmanLab/image_analysis/pipelines/"
    "BulkFluo_RDF/cache/spotiflow_models"
)


@dataclass(frozen=True)
class PunctaTools2DParams:
    """Legacy PunctaTools contract retained for the separate live-cell pipeline."""

    minsize_um: float = 0.2
    maxsize_um: float = 2.0
    num_sigma: int = 5
    overlap: float = 1.0
    threshold_detection: float = 0.001
    threshold_background: float = 0.0
    threshold_segmentation: float = 0.001
    segmentation_mode: int = 0
    maxrad_um: float | None = None


@dataclass(frozen=True)
class SpotiflowHybridParams:
    """Frozen Spotiflow seed and constrained-growth settings."""

    pretrained_model: str = "general"
    probability_threshold: float = 0.5
    min_distance_px: int = 1
    exclude_border: bool = False
    normalizer: str = "auto"
    seed_radius_um: float = 0.234
    grouped_link_distance_um: float = 0.468
    growth_smoothing_sigma_um: float = 0.117
    growth_foreground_sigma: float = 3.0
    maximum_growth_distance_um: float = 1.17
    background_annulus_width_um: float = 0.234
    maximum_group_area_um2: float = 13.69
    minimum_background_pixels: int = 10

    def validate(self) -> None:
        if self.pretrained_model != "general":
            raise ValueError("The phase-diagram pipeline requires Spotiflow general")
        if self.probability_threshold != 0.5:
            raise ValueError("The Spotiflow probability threshold must remain 0.5")
        if self.min_distance_px != 1 or self.exclude_border:
            raise ValueError("Spotiflow default peak settings were changed")
        if self.normalizer != "auto":
            raise ValueError("Spotiflow general requires its default auto normalizer")
        positive = (
            self.seed_radius_um,
            self.grouped_link_distance_um,
            self.growth_smoothing_sigma_um,
            self.growth_foreground_sigma,
            self.maximum_growth_distance_um,
            self.background_annulus_width_um,
            self.maximum_group_area_um2,
        )
        if any(not np.isfinite(value) or value <= 0 for value in positive):
            raise ValueError("Spotiflow hybrid physical parameters must be positive")
        if self.minimum_background_pixels < 1:
            raise ValueError("minimum_background_pixels must be positive")


def make_core_crop(mask: np.ndarray, radius_fraction: float = 1.0 / 3.0) -> np.ndarray:
    """Recreate the area-equivalent-radius erosion on a cropped nucleus mask."""

    mask = np.asarray(mask, dtype=bool)
    if mask.ndim != 2:
        raise ValueError(f"Expected a 2D nucleus mask, got {mask.shape}")
    if not np.any(mask):
        return np.zeros_like(mask)
    equivalent_radius = float(np.sqrt(np.sum(mask) / np.pi))
    return ndimage.distance_transform_edt(mask) > equivalent_radius * radius_fraction


def _laplace_scale_space(
    image: np.ndarray,
    minsize_um: float,
    maxsize_um: float,
    num_sigma: int,
    pixel_size_um: float,
) -> np.ndarray:
    response = np.zeros(image.shape, dtype=np.float32)
    for sigma_um in np.linspace(minsize_um, maxsize_um, int(num_sigma), endpoint=True):
        gaussian = filters.gaussian(image, sigma=float(sigma_um / pixel_size_um))
        response = np.maximum(
            response,
            filters.laplace(gaussian).astype(np.float32, copy=False),
        )
    return response


def punctatools_segment_2d(
    image: np.ndarray,
    core_mask: np.ndarray,
    *,
    pixel_size_um: float,
    params: PunctaTools2DParams,
    global_background: float | None = None,
    global_log_background: float | None = None,
) -> np.ndarray:
    """Legacy PunctaTools implementation used by another SACD workflow."""

    image = np.asarray(image, dtype=np.float32)
    core = np.asarray(core_mask, dtype=bool)
    if image.ndim != 2 or image.shape != core.shape:
        raise ValueError("PunctaTools image and core mask must be aligned")
    if not np.all(np.isfinite(image)) or pixel_size_um <= 0:
        raise ValueError("PunctaTools input must be finite with positive pixel size")
    blobs = blob_log(
        image,
        min_sigma=params.minsize_um / pixel_size_um,
        max_sigma=params.maxsize_um / pixel_size_um,
        num_sigma=int(params.num_sigma),
        overlap=params.overlap,
        threshold=params.threshold_detection,
    )
    markers = np.zeros(image.shape, dtype=np.uint8)
    if len(blobs):
        indices = np.rint(blobs[:, :2]).astype(int)
        indices[:, 0] = np.clip(indices[:, 0], 0, image.shape[0] - 1)
        indices[:, 1] = np.clip(indices[:, 1], 0, image.shape[1] - 1)
        keep = core[indices[:, 0], indices[:, 1]]
        if params.threshold_background > 0:
            if global_background is None:
                raise ValueError("global_background is required")
            keep &= (
                image[indices[:, 0], indices[:, 1]]
                > global_background * params.threshold_background
            )
        indices = indices[keep]
        markers[indices[:, 0], indices[:, 1]] = 1
    marker_labels = ndimage.label(markers)[0]
    if params.segmentation_mode in (0, 1):
        response = _laplace_scale_space(
            image,
            params.minsize_um,
            params.maxsize_um,
            params.num_sigma,
            pixel_size_um,
        )
        threshold = params.threshold_segmentation
        if params.segmentation_mode == 1:
            if global_log_background is None:
                raise ValueError("global_log_background is required")
            threshold *= global_log_background
        foreground = response > threshold
    elif params.segmentation_mode == 2:
        if global_background is None:
            raise ValueError("global_background is required")
        foreground = image > params.threshold_segmentation * global_background
    else:
        raise ValueError("segmentation_mode must be 0, 1, or 2")
    foreground &= core
    distance = ndimage.distance_transform_edt(
        foreground,
        sampling=(pixel_size_um, pixel_size_um),
    )
    labels = watershed(-distance, marker_labels, mask=foreground)
    if params.maxrad_um is not None and np.any(labels):
        ids = np.unique(labels)[1:]
        areas = ndimage.sum(labels > 0, labels, ids) * pixel_size_um**2
        too_large = ids[np.asarray(areas) > np.pi * params.maxrad_um**2]
        labels[np.isin(labels, too_large)] = 0
    return relabel_sequential(labels)[0].astype(np.uint32, copy=False)


class _UnionFind:
    def __init__(self, size: int):
        self.parent = list(range(size))

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> None:
        a, b = self.find(left), self.find(right)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def robust_local_background(values: np.ndarray) -> tuple[float, float]:
    """Return upper-clipped median and MAD noise in the configured SACD intensity space."""

    data = np.asarray(values, dtype=np.float32)
    data = data[np.isfinite(data)]
    if not data.size:
        return float("nan"), float("nan")
    for _ in range(3):
        median = float(np.median(data))
        noise = float(1.4826 * np.median(np.abs(data - median)))
        if noise <= 0:
            break
        clipped = data[data <= median + 3.0 * noise]
        if not clipped.size or clipped.size == data.size:
            break
        data = clipped
    median = float(np.median(data))
    noise = float(1.4826 * np.median(np.abs(data - median)))
    return median, max(noise, float(np.finfo(np.float32).eps))


def _details_value(details: Any, name: str, index: int) -> float:
    value = getattr(details, name, None)
    if value is None or len(value) <= index:
        return float("nan")
    current = np.asarray(value[index]).squeeze()
    return float(current) if current.size == 1 else float("nan")


def spotiflow_seed_table(
    image: np.ndarray,
    core_mask: np.ndarray,
    model: Any,
    *,
    params: SpotiflowHybridParams,
    device: str = "auto",
) -> list[dict[str, Any]]:
    """Detect default-p0.5 Spotiflow seeds and retain centers inside the core."""

    params.validate()
    values = np.asarray(image, dtype=np.float32)
    core = np.asarray(core_mask, dtype=bool)
    if values.ndim != 2 or values.shape != core.shape:
        raise ValueError("Spotiflow image and core must be aligned 2D arrays")
    if not np.all(np.isfinite(values)):
        raise ValueError("Spotiflow input contains nonfinite SACD values")
    points, details = model.predict(
        values,
        prob_thresh=params.probability_threshold,
        min_distance=params.min_distance_px,
        exclude_border=params.exclude_border,
        normalizer=params.normalizer,
        device=device,
        verbose=False,
    )
    retained: list[dict[str, Any]] = []
    for source_index, point in enumerate(np.asarray(points)):
        y, x = float(point[0]), float(point[1])
        yi, xi = int(np.rint(y)), int(np.rint(x))
        if not (0 <= yi < core.shape[0] and 0 <= xi < core.shape[1]):
            continue
        if not core[yi, xi]:
            continue
        retained.append(
            {
                "source_seed_index": source_index,
                "y_sacd_px": y,
                "x_sacd_px": x,
                "rounded_y_sacd_px": yi,
                "rounded_x_sacd_px": xi,
                "spotiflow_probability": _details_value(
                    details,
                    "prob",
                    source_index,
                ),
                "spotiflow_model_intensity": _details_value(
                    details,
                    "intens",
                    source_index,
                ),
            }
        )
    retained.sort(
        key=lambda row: (
            row["y_sacd_px"],
            row["x_sacd_px"],
            row["source_seed_index"],
        )
    )
    for seed_id, row in enumerate(retained, start=1):
        row["seed_id"] = seed_id
    return retained


def group_spotiflow_seeds(
    seeds: list[dict[str, Any]],
    *,
    link_distance_px: float,
) -> list[list[int]]:
    """Return deterministic seed-index groups linked in Euclidean distance."""

    if link_distance_px <= 0:
        raise ValueError("link_distance_px must be positive")
    if not seeds:
        return []
    coordinates = np.asarray(
        [[row["y_sacd_px"], row["x_sacd_px"]] for row in seeds],
        dtype=float,
    )
    union = _UnionFind(len(seeds))
    if len(seeds) > 1:
        for left, right in sorted(cKDTree(coordinates).query_pairs(link_distance_px)):
            union.union(int(left), int(right))
    grouped: dict[int, list[int]] = {}
    for index in range(len(seeds)):
        grouped.setdefault(union.find(index), []).append(index)
    return [grouped[key] for key in sorted(grouped)]


def _seed_distance(
    shape: tuple[int, int],
    seeds: list[dict[str, Any]],
    group: list[int],
    *,
    pixel_size_um: float,
) -> np.ndarray:
    markers = np.zeros(shape, dtype=bool)
    for index in group:
        row = seeds[index]
        markers[row["rounded_y_sacd_px"], row["rounded_x_sacd_px"]] = True
    return ndimage.distance_transform_edt(
        ~markers,
        sampling=(pixel_size_um, pixel_size_um),
    )


def construct_spotiflow_hybrid_mask(
    image: np.ndarray,
    core_mask: np.ndarray,
    seeds: list[dict[str, Any]],
    *,
    pixel_size_um: float,
    params: SpotiflowHybridParams,
) -> tuple[np.ndarray, list[dict[str, Any]], dict[int, dict[str, Any]]]:
    """Group seeds and form bounded, core-confined punctum masks."""

    params.validate()
    values = np.asarray(image, dtype=np.float32)
    core = np.asarray(core_mask, dtype=bool)
    if values.shape != core.shape or values.ndim != 2:
        raise ValueError("Hybrid image and core must be aligned 2D arrays")
    if pixel_size_um <= 0:
        raise ValueError("pixel_size_um must be positive")
    if not seeds:
        return np.zeros(values.shape, dtype=np.uint32), [], {}

    groups = group_spotiflow_seeds(
        seeds,
        link_distance_px=params.grouped_link_distance_um / pixel_size_um,
    )
    smoothed = ndimage.gaussian_filter(
        values,
        sigma=params.growth_smoothing_sigma_um / pixel_size_um,
        mode="nearest",
    )
    union = np.zeros(values.shape, dtype=bool)
    group_rows: list[dict[str, Any]] = []
    group_masks: dict[int, np.ndarray] = {}
    for group_id, group in enumerate(groups, start=1):
        distance = _seed_distance(
            values.shape,
            seeds,
            group,
            pixel_size_um=pixel_size_um,
        )
        base = (distance <= params.seed_radius_um) & core
        search = (distance <= params.maximum_growth_distance_um) & core
        inner = max(
            params.seed_radius_um,
            params.maximum_growth_distance_um - params.background_annulus_width_um,
        )
        annulus = (
            (distance > inner)
            & (distance <= params.maximum_growth_distance_um)
            & core
        )
        background, noise = robust_local_background(smoothed[annulus])
        threshold = (
            background + params.growth_foreground_sigma * noise
            if np.isfinite(background) and np.isfinite(noise)
            else float("nan")
        )
        status = "grown"
        candidate = np.zeros_like(core)
        if int(np.sum(annulus)) < params.minimum_background_pixels:
            status = "fallback_insufficient_background"
        elif not np.isfinite(threshold):
            status = "fallback_invalid_background"
        else:
            foreground = (smoothed > threshold) & search
            components, _ = ndimage.label(foreground)
            component_ids = {
                int(components[
                    seeds[index]["rounded_y_sacd_px"],
                    seeds[index]["rounded_x_sacd_px"],
                ])
                for index in group
            }
            component_ids.discard(0)
            if component_ids:
                candidate = np.isin(components, sorted(component_ids))
            if not np.any(candidate):
                status = "fallback_seed_below_local_threshold"
            elif float(np.sum(candidate) * pixel_size_um**2) > params.maximum_group_area_um2:
                status = "fallback_growth_area_exceeded"
                candidate[:] = False
        final_group = base | candidate
        union |= final_group
        group_masks[group_id] = final_group
        member_seed_ids = [int(seeds[index]["seed_id"]) for index in group]
        group_rows.append(
            {
                "group_id": group_id,
                "member_seed_ids": ";".join(map(str, member_seed_ids)),
                "member_seed_count": len(member_seed_ids),
                "local_background_sacd": background,
                "local_noise_sacd": noise,
                "growth_threshold_sacd": threshold,
                "background_pixel_count": int(np.sum(annulus)),
                "base_area_sacd_px": int(np.sum(base)),
                "candidate_growth_area_sacd_px": int(np.sum(candidate)),
                "final_group_area_sacd_px": int(np.sum(final_group)),
                "growth_status": status,
            }
        )
        for index in group:
            seeds[index]["group_id"] = group_id
            seeds[index]["group_member_count"] = len(group)
            seeds[index]["growth_status"] = status

    labels = connected_components(union, connectivity=2).astype(np.uint32)
    labels = relabel_sequential(labels)[0].astype(np.uint32, copy=False)
    provenance: dict[int, dict[str, Any]] = {}
    group_by_id = {int(row["group_id"]): row for row in group_rows}
    for row in seeds:
        punctum_id = int(
            labels[row["rounded_y_sacd_px"], row["rounded_x_sacd_px"]]
        )
        if punctum_id <= 0:
            raise AssertionError("A retained Spotiflow seed lost its fallback footprint")
        row["punctum_id"] = punctum_id
        record = provenance.setdefault(
            punctum_id,
            {"seed_ids": [], "group_ids": set(), "growth_statuses": set()},
        )
        record["seed_ids"].append(int(row["seed_id"]))
        record["group_ids"].add(int(row["group_id"]))
        record["growth_statuses"].add(str(row["growth_status"]))
    for punctum_id, record in provenance.items():
        group_ids = sorted(record["group_ids"])
        record["seed_ids"] = sorted(record["seed_ids"])
        record["group_ids"] = group_ids
        record["growth_statuses"] = sorted(record["growth_statuses"])
        record["member_seed_count"] = len(record["seed_ids"])
        record["merged_group_count"] = len(group_ids)
        record["growth_status"] = ";".join(record["growth_statuses"])
        record["group_member_seed_ids"] = ";".join(
            group_by_id[group_id]["member_seed_ids"] for group_id in group_ids
        )
    if labels.dtype != np.uint32 or np.any(labels[~core]):
        raise AssertionError("Hybrid labels are not uint32 and core-confined")
    return labels, group_rows, provenance


def quantify_nucleus_puncta(
    image: np.ndarray,
    core_mask: np.ndarray,
    puncta_labels: np.ndarray,
    *,
    pixel_size_um: float,
    provenance: Mapping[int, Mapping[str, Any]] | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Measure final puncta on unnormalized values in the configured SACD intensity space."""

    image = np.asarray(image, dtype=np.float64)
    core = np.asarray(core_mask, dtype=bool)
    labels = np.asarray(puncta_labels, dtype=np.uint32)
    if image.shape != core.shape or image.shape != labels.shape:
        raise ValueError("Image, core, and puncta labels must be aligned")
    if np.any(labels[~core] > 0):
        raise ValueError("Puncta labels extend outside the eroded core")
    if not np.any(core):
        raise ValueError("Cannot quantify an empty core")
    puncta_mask = labels > 0
    nonpuncta = core & ~puncta_mask
    nonpuncta_mean = (
        float(np.mean(image[nonpuncta]))
        if np.any(nonpuncta)
        else float("nan")
    )
    boundary = core & ~ndimage.binary_erosion(core)
    rows: list[dict[str, Any]] = []
    for punctum_id in (int(value) for value in np.unique(labels) if value > 0):
        mask = labels == punctum_id
        values = image[mask]
        corrected = (
            np.maximum(values - nonpuncta_mean, 0.0)
            if np.isfinite(nonpuncta_mean)
            else np.full_like(values, np.nan)
        )
        ys, xs = np.where(mask)
        row: dict[str, Any] = {
            "punctum_id": punctum_id,
            "area_sacd_px": int(values.size),
            "area_um2": float(values.size * pixel_size_um**2),
            "mean_sacd_intensity": float(np.mean(values)),
            "max_sacd_intensity": float(np.max(values)),
            "raw_integrated_sacd_intensity": float(np.sum(values)),
            "background_corrected_integrated_sacd_intensity": float(
                np.sum(corrected)
            ),
            "centroid_y_sacd_px": float(np.mean(ys)),
            "centroid_x_sacd_px": float(np.mean(xs)),
            "touches_core_boundary": bool(np.any(mask & boundary)),
        }
        source = dict((provenance or {}).get(punctum_id, {}))
        if source:
            row.update(
                {
                    "member_seed_count": source["member_seed_count"],
                    "member_seed_ids": ";".join(map(str, source["seed_ids"])),
                    "merged_group_count": source["merged_group_count"],
                    "member_group_ids": ";".join(map(str, source["group_ids"])),
                    "growth_status": source["growth_status"],
                }
            )
        rows.append(row)
    count = len(rows)
    summary = {
        "puncta_count": count,
        "puncta_pixel_mean_sacd": (
            float(np.mean(image[puncta_mask])) if count else float("nan")
        ),
        "nonpuncta_core_mean_sacd": (
            nonpuncta_mean if count else float(np.mean(image[core]))
        ),
        "max_raw_integrated_punctum_sacd": max(
            (float(row["raw_integrated_sacd_intensity"]) for row in rows),
            default=0.0,
        ),
        "max_background_corrected_integrated_punctum_sacd": max(
            (
                float(row["background_corrected_integrated_sacd_intensity"])
                for row in rows
            ),
            default=0.0,
        ),
        "puncta_area_fraction": float(np.sum(puncta_mask) / np.sum(core)),
    }
    return summary, rows


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with temporary.open("w", newline="") as handle:
        if fieldnames:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
    os.replace(temporary, path)


def _write_json(path: Path, value: Any) -> None:
    def json_safe(current: Any) -> Any:
        if isinstance(current, Mapping):
            return {str(key): json_safe(item) for key, item in current.items()}
        if isinstance(current, (list, tuple)):
            return [json_safe(item) for item in current]
        if isinstance(current, np.generic):
            current = current.item()
        if isinstance(current, float) and not np.isfinite(current):
            return None
        return current

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(
            json_safe(value),
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n"
    )
    os.replace(temporary, path)


def _pixel_size_um_from_tiff(path: str | Path) -> float:
    with tifffile.TiffFile(path) as tif:
        value = tif.pages[0].tags["XResolution"].value
    pixels_per_um = (
        float(value[0]) / float(value[1])
        if isinstance(value, tuple)
        else float(value)
    )
    if pixels_per_um <= 0:
        raise ValueError(f"Invalid TIFF resolution in {path}")
    return 1.0 / pixels_per_um


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _aggregate_file_digest(paths: Iterable[Path], root: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    count = 0
    for path in sorted(paths):
        count += 1
        digest.update(str(path.relative_to(root)).encode())
        digest.update(bytes.fromhex(_sha256_file(path)))
    return {"count": count, "aggregate_sha256": digest.hexdigest()}


def phase_dataset_invariants(
    output_root: str | Path,
    *,
    nucleus_paths: Iterable[str | Path] | None = None,
) -> dict[str, dict[str, Any]]:
    """Hash protected TIFF families and FOV-level QC overlays."""

    root = Path(output_root)
    protected_nuclei = (
        [Path(path) for path in nucleus_paths]
        if nucleus_paths is not None
        else list((root / "nuclei").rglob("*.tif"))
    )
    return {
        "nucleus_tiffs": _aggregate_file_digest(
            protected_nuclei,
            root,
        ),
        "reconstruction_tiffs": _aggregate_file_digest(
            (root / "reconstructions").rglob("*.tif"),
            root,
        ),
        "fov_qc_overlays": _aggregate_file_digest(
            (root / "qc").glob("*/*__CellposeSAM-overlay.png"),
            root,
        ),
    }


def _assign_expression_deciles(
    records: list[dict[str, Any]],
) -> dict[str, int]:
    ordered = sorted(
        records,
        key=lambda row: (
            float(row["spen_sacd_core_mean"]),
            str(row["nucleus_key"]),
        ),
    )
    return {
        str(row["nucleus_key"]): min(9, int(index * 10 / len(ordered)))
        for index, row in enumerate(ordered)
    }


def select_comparison_sample(
    records: list[dict[str, Any]],
    *,
    seed: int = 20260701,
) -> list[dict[str, Any]]:
    """Select 35 deterministic cells per condition across expression deciles."""

    rng = np.random.default_rng(seed)
    chosen: list[dict[str, Any]] = []
    targets = (
        list(range(10))
        + [0, 2, 4, 6, 8]
        + [1, 3, 5, 7, 9]
        + list(range(10))
        + [1, 3, 5, 7, 9]
    )
    for condition in ("dSPEN_FL", "dSPEN_dRRM"):
        condition_rows = [
            dict(row) for row in records
            if row["pooled_condition"] == condition
        ]
        deciles = _assign_expression_deciles(condition_rows)
        groups = {value: [] for value in range(10)}
        for row in condition_rows:
            decile = deciles[str(row["nucleus_key"])]
            row["expression_decile"] = decile
            groups[decile].append(row)
        for group in groups.values():
            rng.shuffle(group)
        for target in targets:
            available = [value for value, group in groups.items() if group]
            decile = min(available, key=lambda value: (abs(value - target), value))
            chosen.append(groups[decile].pop())
    for order, row in enumerate(chosen, start=1):
        row["comparison_order"] = order
    return chosen


def _draw_review_panel(
    ax: Any,
    image: np.ndarray,
    nucleus: np.ndarray,
    core: np.ndarray,
    labels: np.ndarray | None,
    *,
    vmin: float,
    vmax: float,
    title: str,
) -> None:
    from matplotlib.colors import LogNorm

    ax.imshow(
        image,
        cmap="magma",
        norm=LogNorm(vmin=vmin, vmax=vmax, clip=True),
        interpolation="nearest",
    )
    ax.contour(nucleus.astype(np.uint8), [0.5], colors=["cyan"], linewidths=0.6)
    ax.contour(
        core.astype(np.uint8),
        [0.5],
        colors=["lime"],
        linewidths=0.7,
        linestyles=":",
    )
    if labels is not None and np.any(labels):
        ax.contour(
            (labels > 0).astype(np.uint8),
            [0.5],
            colors=["white"],
            linewidths=0.6,
        )
    ax.set_title(title, fontsize=6.5)
    ax.axis("off")


def _render_spotiflow_qc_montages(
    stage_qc: Path,
    qc_cells: list[dict[str, Any]],
    *,
    review_vmin: float,
    review_vmax: float,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for page_index in range(7):
        page = qc_cells[page_index * 10 : (page_index + 1) * 10]
        fig, axes = plt.subplots(10, 2, figsize=(7.2, 24), constrained_layout=True)
        for row_index, cell in enumerate(page):
            record = cell["record"]
            base = (
                f"{record['pooled_condition']} {record['dose']} "
                f"{record['fov_token']} n{int(record['nucleus_id']):04d}\n"
                f"core mean={float(record['spen_sacd_core_mean']):.3g}"
            )
            _draw_review_panel(
                axes[row_index, 0],
                cell["image"],
                cell["nucleus"],
                cell["core"],
                None,
                vmin=review_vmin,
                vmax=review_vmax,
                title=base + "\nSPEN SACD",
            )
            _draw_review_panel(
                axes[row_index, 1],
                cell["image"],
                cell["nucleus"],
                cell["core"],
                cell["labels"],
                vmin=review_vmin,
                vmax=review_vmax,
                title=(
                    base
                    + f"\nSpotiflow seeds={cell['seed_count']}, "
                    f"puncta={cell['puncta_count']}, "
                    f"area={cell['area_fraction']:.3f}"
                ),
            )
        fig.suptitle(
            f"Spotiflow p=0.5 constrained-hybrid QC — page {page_index + 1}/7",
            fontsize=12,
        )
        path = stage_qc / "montages" / (
            f"spotiflow-hybrid-page-{page_index + 1:02d}-of-07.png"
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=160, facecolor="white")
        plt.close(fig)


def _build_puncta_phase_figure(
    records: list[dict[str, Any]],
) -> tuple[Any, np.ndarray]:
    """Build shared-scale condition facets for the three retained metrics."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    panels = (
        ("puncta_count", "Puncta count"),
        ("puncta_pixel_mean_sacd", "Mean puncta-pixel intensity (a.u.)"),
        ("nonpuncta_core_mean_sacd", "Mean non-puncta core intensity (a.u.)"),
    )
    conditions = ("dSPEN_FL", "dSPEN_dRRM")
    palette = {"dSPEN_FL": "#2468b4", "dSPEN_dRRM": "#e4572e"}
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(15.5, 9.5),
        sharex=True,
        sharey="col",
    )
    for row_index, condition in enumerate(conditions):
        selected = [
            row for row in records if row["pooled_condition"] == condition
        ]
        x = np.asarray(
            [float(row["spen_sacd_core_mean"]) for row in selected]
        )
        for column_index, (column, _) in enumerate(panels):
            ax = axes[row_index, column_index]
            y = np.asarray(
                [float(row[column]) for row in selected],
                dtype=np.float64,
            )
            finite = np.isfinite(x) & np.isfinite(y) & (x > 0)
            if column_index > 0:
                finite &= y > 0
            ax.scatter(
                x[finite],
                y[finite],
                s=13,
                alpha=0.55,
                linewidths=0,
                color=palette[condition],
                rasterized=True,
            )
            if column == "puncta_pixel_mean_sacd":
                displayed = int(np.sum(finite))
                omitted = len(selected) - displayed
                ax.text(
                    0.97,
                    0.04,
                    f"n={displayed}\n{omitted} zero-puncta omitted",
                    transform=ax.transAxes,
                    ha="right",
                    va="bottom",
                    fontsize=8,
                )
        axes[row_index, 0].annotate(
            condition,
            xy=(-0.17, 0.5),
            xycoords="axes fraction",
            rotation=90,
            ha="center",
            va="center",
            fontsize=12,
        )

    for column_index, (_, ylabel) in enumerate(panels):
        axes[0, column_index].set_title(ylabel)
        for row_index in range(2):
            ax = axes[row_index, column_index]
            ax.set_xscale("log")
            ax.set_yscale("linear" if column_index == 0 else "log")
            ax.set_ylabel(ylabel)
            ax.grid(alpha=0.2, linewidth=0.6)
    count_max = max(float(row["puncta_count"]) for row in records)
    axes[0, 0].set_ylim(0, max(1.0, count_max * 1.05))
    for ax in axes[1]:
        ax.set_xlabel(
            "Mean linear SPEN SACD intensity in eroded core (a.u.; log axis)"
        )

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markersize=6,
            markerfacecolor=palette[condition],
            markeredgewidth=0,
            label=(
                f"{condition} "
                f"(n={sum(row['pooled_condition'] == condition for row in records)})"
            )
        )
        for condition in conditions
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.012),
        ncol=2,
        frameon=False,
    )
    fig.suptitle("SPEN 2D puncta phase diagram — Spotiflow hybrid")
    fig.subplots_adjust(
        left=0.09,
        right=0.985,
        top=0.91,
        bottom=0.13,
        wspace=0.25,
        hspace=0.23,
    )
    return fig, axes


def make_puncta_phase_diagram(
    records: list[dict[str, Any]],
    path: str | Path,
) -> None:
    """Plot three puncta metrics as shared-scale condition facets."""

    import matplotlib.pyplot as plt

    fig, _ = _build_puncta_phase_figure(records)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp.png")
    fig.savefig(temporary, dpi=180, facecolor="white")
    plt.close(fig)
    os.replace(temporary, path)


def load_spotiflow_general(
    *,
    cache_dir: str | Path = DEFAULT_SPOTIFLOW_CACHE,
) -> tuple[Any, dict[str, Any]]:
    """Load and verify the installed pretrained Spotiflow general model."""

    from spotiflow import __version__ as spotiflow_version
    from spotiflow.model import Spotiflow

    cache_dir = Path(cache_dir)
    weight_path = cache_dir / "general" / "best.pt"
    if not weight_path.is_file():
        raise FileNotFoundError(f"Missing cached Spotiflow weights: {weight_path}")
    weight_sha = _sha256_file(weight_path)
    if weight_sha != SPOTIFLOW_GENERAL_WEIGHT_SHA256:
        raise ValueError(
            f"Unexpected Spotiflow general weight hash: {weight_sha}"
        )
    model = Spotiflow.from_pretrained(
        "general",
        cache_dir=cache_dir,
        map_location="cpu",
        verbose=False,
    )
    metadata = {
        "spotiflow_version": str(spotiflow_version),
        "spotiflow_distribution_version": importlib.metadata.version("spotiflow"),
        "pretrained_model": "general",
        "model_weight_path": str(weight_path),
        "model_weight_sha256": weight_sha,
    }
    return model, metadata


def _move_to_backup(
    final: Path,
    backup_root: Path,
    output_root: Path,
) -> Path | None:
    if not final.exists():
        return None
    backup = backup_root / final.relative_to(output_root)
    backup.parent.mkdir(parents=True, exist_ok=True)
    os.replace(final, backup)
    return backup


@report_progress
def apply_spotiflow_hybrid(
    records: list[dict[str, Any]],
    output_root: str | Path,
    *,
    expected_nuclei: int = 1231,
    comparison_seed: int = 20260701,
    device: str = "auto",
    cache_dir: str | Path = DEFAULT_SPOTIFLOW_CACHE,
    params: SpotiflowHybridParams | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Run, validate, and atomically publish the complete hybrid analysis."""

    from .phase_diagram import _make_nucleus_review_png

    params = params or SpotiflowHybridParams()
    params.validate()
    if len(records) != expected_nuclei:
        raise AssertionError(
            f"Expected {expected_nuclei:,} retained nuclei, found {len(records):,}"
        )
    output_root = Path(output_root)
    protected_nucleus_paths = [row["review_tif"] for row in records]
    if len(set(protected_nucleus_paths)) != expected_nuclei:
        raise AssertionError("Nucleus manifest does not contain unique review TIFF paths")
    stage("Validate protected inputs")
    protected_before = phase_dataset_invariants(
        output_root,
        nucleus_paths=protected_nucleus_paths,
    )
    expected_protected = {
        "nucleus_tiffs": 1231,
        "reconstruction_tiffs": 574,
        "fov_qc_overlays": 82,
    }
    for family, expected in expected_protected.items():
        if protected_before[family]["count"] != expected:
            raise AssertionError(
                f"Protected {family} count is "
                f"{protected_before[family]['count']}, expected {expected}"
            )

    stage("Load Spotiflow model")
    model, model_metadata = load_spotiflow_general(cache_dir=cache_dir)
    attempt = output_root / ".work" / f"spotiflow-hybrid-{uuid.uuid4().hex[:8]}"
    stage_masks = attempt / "puncta_masks"
    stage_pngs = attempt / "nucleus_pngs"
    stage_manifests = attempt / "manifests"
    stage_qc = attempt / "qc" / "spotiflow_hybrid"
    sample_keys = {
        row["nucleus_key"]
        for row in select_comparison_sample(records, seed=comparison_seed)
    }
    nucleus_rows: list[dict[str, Any]] = []
    puncta_rows: list[dict[str, Any]] = []
    seed_rows: list[dict[str, Any]] = []
    qc_cells: list[dict[str, Any]] = []
    staged_mask_paths: list[Path] = []
    staged_png_paths: list[Path] = []
    started = datetime.now(timezone.utc)
    try:
        for index, source in enumerate(progress_items(records, "Spotiflow nuclei"), start=1):
            review = tifffile.imread(source["review_tif"])
            if review.ndim != 3 or review.shape[0] != 3:
                raise ValueError(f"Invalid review TIFF: {source['review_tif']}")
            image = np.asarray(review[1], dtype=np.float32)
            nucleus = np.asarray(review[2] > 0.5, dtype=bool)
            core = make_core_crop(nucleus)
            pixel_size_um = _pixel_size_um_from_tiff(source["review_tif"])
            seeds = spotiflow_seed_table(
                image,
                core,
                model,
                params=params,
                device=device,
            )
            labels, groups, provenance = construct_spotiflow_hybrid_mask(
                image,
                core,
                seeds,
                pixel_size_um=pixel_size_um,
                params=params,
            )
            summary, objects = quantify_nucleus_puncta(
                image,
                core,
                labels,
                pixel_size_um=pixel_size_um,
                provenance=provenance,
            )
            mask_name = (
                Path(source["review_tif"]).stem
                + "__SPEN-puncta-labels-YX.tif"
            )
            staged_mask = (
                stage_masks / source["pooled_condition"] / mask_name
            )
            final_mask = (
                output_root
                / "puncta_masks"
                / source["pooled_condition"]
                / mask_name
            )
            staged_mask.parent.mkdir(parents=True, exist_ok=True)
            tifffile.imwrite(
                staged_mask,
                labels,
                metadata={"axes": "YX", "unit": "um"},
                resolution=(1.0 / pixel_size_um, 1.0 / pixel_size_um),
                photometric="minisblack",
            )
            staged_mask_paths.append(staged_mask)
            updated = dict(source)
            updated.update(summary)
            updated.update(
                {
                    "puncta_detector": "spotiflow_general_p0.5_constrained_hybrid",
                    "puncta_detector_parameters_json": json.dumps(
                        asdict(params),
                        sort_keys=True,
                    ),
                    "spotiflow_version": model_metadata["spotiflow_version"],
                    "spotiflow_model_weight_sha256": model_metadata[
                        "model_weight_sha256"
                    ],
                    "spotiflow_seed_count": len(seeds),
                    "spotiflow_group_count": len(groups),
                    "spotiflow_grown_group_count": sum(
                        row["growth_status"] == "grown" for row in groups
                    ),
                    "spotiflow_fallback_group_count": sum(
                        row["growth_status"] != "grown" for row in groups
                    ),
                    "puncta_mask_path": str(final_mask),
                }
            )
            nucleus_rows.append(updated)
            for seed_row in seeds:
                seed_rows.append(
                    {
                        "nucleus_key": source["nucleus_key"],
                        "pooled_condition": source["pooled_condition"],
                        "condition": source["condition"],
                        "dose": source["dose"],
                        "fov": source["fov"],
                        "nucleus_id": source["nucleus_id"],
                        **seed_row,
                        "puncta_mask_path": str(final_mask),
                    }
                )
            for obj in objects:
                puncta_rows.append(
                    {
                        "nucleus_key": source["nucleus_key"],
                        "pooled_condition": source["pooled_condition"],
                        "condition": source["condition"],
                        "dose": source["dose"],
                        "fov": source["fov"],
                        "nucleus_id": source["nucleus_id"],
                        **obj,
                        "puncta_mask_path": str(final_mask),
                        "detector": "spotiflow_general_p0.5_constrained_hybrid",
                    }
                )
            staged_png = (
                stage_pngs
                / source["pooled_condition"]
                / Path(source["review_png"]).name
            )
            _make_nucleus_review_png(
                staged_png,
                image,
                nucleus,
                core,
                pixel_size_um=pixel_size_um,
                scale_bar_um=2.0,
                lower_percentile=3.0,
                upper_percentile=99.5,
                display_limits=(
                    float(source["review_vmin_sacd"]),
                    float(source["review_vmax_sacd"]),
                ),
                puncta_labels=labels,
            )
            staged_png_paths.append(staged_png)
            if source["nucleus_key"] in sample_keys:
                qc_cells.append(
                    {
                        "record": source,
                        "image": image,
                        "nucleus": nucleus,
                        "core": core,
                        "labels": labels,
                        "seed_count": len(seeds),
                        "puncta_count": summary["puncta_count"],
                        "area_fraction": summary["puncta_area_fraction"],
                    }
                )
            if index % 25 == 0 or index == len(records):
                progress_message(
                    f"Spotiflow hybrid {index}/{len(records)}: "
                    f"{len(seed_rows):,} retained seeds, "
                    f"{len(puncta_rows):,} puncta",
                    flush=True,
                    update_stage=False,
                )

        stage("Render Spotiflow QC and prepare manifests")
        if len(qc_cells) != 70:
            raise AssertionError(f"Expected 70 QC cells, found {len(qc_cells)}")
        qc_order = {
            row["nucleus_key"]: row["comparison_order"]
            for row in select_comparison_sample(records, seed=comparison_seed)
        }
        qc_cells.sort(key=lambda cell: qc_order[cell["record"]["nucleus_key"]])
        _render_spotiflow_qc_montages(
            stage_qc,
            qc_cells,
            review_vmin=float(records[0]["review_vmin_sacd"]),
            review_vmax=float(records[0]["review_vmax_sacd"]),
        )
        _write_csv(stage_manifests / "nucleus_manifest.csv", nucleus_rows)
        _write_csv(stage_manifests / "puncta_manifest.csv", puncta_rows)
        _write_csv(
            stage_manifests / "spotiflow_seed_manifest.csv",
            seed_rows,
        )
        metadata = {
            "schema_version": 1,
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "started_at": started.isoformat(),
            "detector": "spotiflow_general_p0.5_constrained_hybrid",
            "parameters": asdict(params),
            "model": model_metadata,
            "retained_nuclei": len(records),
            "retained_seeds": len(seed_rows),
            "final_puncta": len(puncta_rows),
            "zero_puncta_nuclei": sum(
                int(row["puncta_count"]) == 0 for row in nucleus_rows
            ),
            "detection_input": (
                "linear SPEN SACD passed to Spotiflow default auto normalizer"
            ),
            "growth_and_measurement_intensity_space": "native_linear_SACD",
            "log_usage": "review_rendering_and_phase_diagram_x_axis_only",
            "protected_before": protected_before,
        }
        _write_json(
            stage_manifests / "spotiflow_hybrid_metadata.json",
            metadata,
        )
        _write_json(stage_qc / "summary.json", metadata)
        make_puncta_phase_diagram(
            nucleus_rows,
            attempt / "qc" / "phase_diagram-SPEN-puncta-metrics.png",
        )

        if (
            len(staged_mask_paths) != expected_nuclei
            or not all(path.is_file() for path in staged_mask_paths)
        ):
            raise AssertionError("Staged puncta-mask count mismatch")
        if (
            len(staged_png_paths) != expected_nuclei
            or not all(path.is_file() for path in staged_png_paths)
        ):
            raise AssertionError("Staged review-PNG count mismatch")
        if len(list((stage_qc / "montages").glob("*.png"))) != 7:
            raise AssertionError("Spotiflow QC must contain seven montage pages")
        rows_by_key = {row["nucleus_key"]: row for row in nucleus_rows}
        for source in progress_items(records, "Validate Spotiflow outputs"):
            row = rows_by_key[source["nucleus_key"]]
            staged = (
                stage_masks
                / row["pooled_condition"]
                / Path(row["puncta_mask_path"]).name
            )
            labels = tifffile.imread(staged)
            review = tifffile.imread(row["review_tif"])
            image = np.asarray(review[1], dtype=np.float32)
            core = make_core_crop(review[2] > 0.5)
            if (
                labels.dtype != np.uint32
                or labels.shape != core.shape
                or np.any(labels[~core])
            ):
                raise AssertionError(
                    f"Invalid staged puncta mask for {row['nucleus_key']}"
                )
            recomputed, _ = quantify_nucleus_puncta(
                image,
                core,
                labels,
                pixel_size_um=_pixel_size_um_from_tiff(row["review_tif"]),
            )
            for key, expected in recomputed.items():
                observed = float(row[key])
                if np.isnan(expected) and np.isnan(observed):
                    continue
                if not np.isclose(observed, expected, rtol=1e-10, atol=1e-6):
                    raise AssertionError(
                        f"Metric mismatch {row['nucleus_key']} {key}"
                    )

        updated_by_key = {row["nucleus_key"]: row for row in nucleus_rows}
        completion_pairs: list[tuple[Path, Path]] = []
        for completion_path in sorted(
            (output_root / "reconstructions").rglob("_COMPLETE.json")
        ):
            completion = json.loads(completion_path.read_text())
            completion["nucleus_records"] = [
                updated_by_key.get(row["nucleus_key"], row)
                for row in completion.get("nucleus_records", [])
            ]
            completion["puncta_detector_metadata"] = metadata
            staged_completion = (
                attempt
                / "checkpoints"
                / completion_path.parent.name
                / completion_path.name
            )
            _write_json(staged_completion, completion)
            completion_pairs.append((staged_completion, completion_path))

        replacements: list[tuple[Path, Path]] = [
            (
                stage_manifests / "nucleus_manifest.csv",
                output_root / "manifests" / "nucleus_manifest.csv",
            ),
            (
                stage_manifests / "puncta_manifest.csv",
                output_root / "manifests" / "puncta_manifest.csv",
            ),
            (
                stage_manifests / "spotiflow_seed_manifest.csv",
                output_root / "manifests" / "spotiflow_seed_manifest.csv",
            ),
            (
                stage_manifests / "spotiflow_hybrid_metadata.json",
                output_root / "manifests" / "spotiflow_hybrid_metadata.json",
            ),
            (
                attempt / "qc" / "phase_diagram-SPEN-puncta-metrics.png",
                output_root / "qc" / "phase_diagram-SPEN-puncta-metrics.png",
            ),
        ]
        replacements.extend(completion_pairs)
        replacements.extend(
            (
                stage_pngs
                / row["pooled_condition"]
                / Path(row["review_png"]).name,
                Path(row["review_png"]),
            )
            for row in nucleus_rows
        )
        directory_replacements = [
            (stage_masks, output_root / "puncta_masks"),
            (stage_qc, output_root / "qc" / "spotiflow_hybrid"),
        ]
        removals = [
            output_root / "qc" / "puncta_benchmark",
            output_root
            / "qc"
            / "phase_diagram-SPEN_core_mean-vs-dispersion.png",
        ]
        stage("Publish Spotiflow outputs")
        backup_root = attempt / "rollback"
        committed_files: list[tuple[Path, Path | None]] = []
        committed_dirs: list[tuple[Path, Path | None]] = []
        removed: list[tuple[Path, Path]] = []
        try:
            for staged, final in directory_replacements:
                old = _move_to_backup(final, backup_root, output_root)
                final.parent.mkdir(parents=True, exist_ok=True)
                os.replace(staged, final)
                committed_dirs.append((final, old))
            for staged, final in replacements:
                old = _move_to_backup(final, backup_root, output_root)
                final.parent.mkdir(parents=True, exist_ok=True)
                os.replace(staged, final)
                committed_files.append((final, old))
            for final in removals:
                old = _move_to_backup(final, backup_root, output_root)
                if old is not None:
                    removed.append((final, old))
            protected_after = phase_dataset_invariants(
                output_root,
                nucleus_paths=protected_nucleus_paths,
            )
            if protected_after != protected_before:
                raise AssertionError(
                    "Protected TIFF or FOV-QC hashes changed during publication"
                )
        except Exception:
            for final, old in reversed(committed_files):
                if final.exists():
                    final.unlink()
                if old is not None and old.exists():
                    final.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(old, final)
            for final, old in reversed(committed_dirs):
                if final.exists():
                    shutil.rmtree(final)
                if old is not None and old.exists():
                    final.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(old, final)
            for final, old in reversed(removed):
                if old.exists():
                    final.parent.mkdir(parents=True, exist_ok=True)
                    os.replace(old, final)
            raise
        else:
            shutil.rmtree(attempt, ignore_errors=True)
        return nucleus_rows, puncta_rows
    except Exception:
        shutil.rmtree(attempt, ignore_errors=True)
        raise


def detector_environment_manifest() -> dict[str, Any]:
    """Return the Spotiflow hybrid provenance contract."""

    return {
        "working_smlm_environment_mutated": False,
        "custom_training": False,
        "detector": "spotiflow_general_p0.5_constrained_hybrid",
        "parameters": asdict(SpotiflowHybridParams()),
        "model_weight_sha256": SPOTIFLOW_GENERAL_WEIGHT_SHA256,
        "intensity_space": {
            "spotiflow_inference": "linear_input_with_model_default_auto_normalizer",
            "growth_and_measurement": "native_linear_SACD",
            "log_usage": "rendering_and_phase_diagram_x_axis_only",
        },
    }
