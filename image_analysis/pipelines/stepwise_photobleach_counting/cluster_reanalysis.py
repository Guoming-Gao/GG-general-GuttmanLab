"""Cluster-aware analysis units and local-background trace extraction.

The revised workflow deliberately separates morphology from eligibility:
every Spotiflow seed receives a seed-level unit and a deterministic grouped
unit whether or not legacy local growth succeeded.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import distance_transform_edt, gaussian_filter
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from skimage.draw import disk
from skimage.segmentation import find_boundaries, watershed

from pbsa_shared import atomic_csv, condition_metadata, matplotlib_setup, require_free_space


STREAMS = ("seed_level_cluster", "grouped_cluster", "extended_cluster_candidate")
TRACE_REPRESENTATIONS = (
    "signal_mean", "background_mean", "mean_difference",
    "integrated_difference", "baseline_centered_integrated_difference",
)
EXTENDED_COLUMNS = [
    "analysis_unit_id", "analysis_stream", "source_seed_id", "source_growth_status",
    "area_px", "equivalent_diameter_px", "centroid_y_px", "centroid_x_px",
    "background_pixels", "extended_segmentation_valid", "extended_rejection_reason",
    "condensate_candidate", "condensate_validation_status", "label_id",
]


def load_reanalysis_config(path: str | Path) -> dict[str, Any]:
    import yaml
    config_path = Path(path).expanduser().resolve()
    with config_path.open() as handle:
        cfg = yaml.safe_load(handle)
    cfg["_config_path"] = str(config_path)
    cfg["source_results_root"] = str(Path(cfg["source_results_root"]).expanduser().resolve())
    cfg["output_root"] = str(Path(cfg["output_root"]).expanduser().resolve())
    if float(cfg["spotiflow"]["probability_threshold"]) != 0.4:
        raise ValueError("The primary Spotiflow threshold must remain 0.4")
    return cfg


def roots(cfg: dict[str, Any]) -> dict[str, Path]:
    root = Path(cfg["output_root"])
    return {
        "root": root,
        "report": root / "01_report",
        "figures": root / "02_figures",
        "tables": root / "03_tables",
        "experiments": root / "04_experiments",
        "technical": root / "99_technical",
    }


def ensure_layout(cfg: dict[str, Any]) -> dict[str, Path]:
    paths = roots(cfg)
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    for name in ("seed_level_clusters", "grouped_clusters", "extended_clusters", "method_QC"):
        (paths["figures"] / name).mkdir(parents=True, exist_ok=True)
    for name in ("provenance", "sensitivity_analysis", "quickpbsa_native", "exploratory_segmentation"):
        (paths["technical"] / name).mkdir(parents=True, exist_ok=True)
    require_free_space(paths["root"], float(cfg["qc"]["minimum_free_space_gib"]))
    return paths


def experiment_root(cfg: dict[str, Any], dataset: str) -> Path:
    return roots(cfg)["experiments"] / dataset


def source_dataset_root(cfg: dict[str, Any], dataset: str) -> Path:
    return Path(cfg["source_results_root"]) / dataset


def source_manifest(cfg: dict[str, Any], dataset: str) -> pd.DataFrame:
    path = source_dataset_root(cfg, dataset) / "01_input_inspection" / "input_manifest.csv"
    if not path.exists(): path = experiment_root(cfg, dataset) / "input_QC" / "input_manifest.csv"
    table = pd.read_csv(path)
    return table[table.status == "accepted"].copy()


def source_seed_table(cfg: dict[str, Any], dataset: str, fov: str) -> pd.DataFrame:
    path = source_dataset_root(cfg, dataset) / "03_roi_detection_and_growth" / "tables" / f"{fov}__spotiflow_seeds.csv"
    if not path.exists(): path = roots(cfg)["technical"] / "exploratory_segmentation" / "legacy_growth_inputs" / dataset / "tables" / f"{fov}__spotiflow_seeds.csv"
    table = pd.read_csv(path)
    return table.rename(columns={"roi_id": "legacy_roi_id", "growth_status": "legacy_growth_status"})


def source_corrected_stack(cfg: dict[str, Any], dataset: str, fov: str) -> Path:
    path = source_dataset_root(cfg, dataset) / "02_drift_correction" / "corrected_tiffs" / f"{fov}__drift_corrected.tif"
    return path if path.exists() else experiment_root(cfg, dataset) / "drift_correction" / "corrected_tiffs" / f"{fov}__drift_corrected.tif"


def source_detection_image(cfg: dict[str, Any], dataset: str, fov: str) -> Path:
    path = source_dataset_root(cfg, dataset) / "03_roi_detection_and_growth" / "detection_images" / f"{fov}__first_frames_mean.tif"
    return path if path.exists() else roots(cfg)["technical"] / "exploratory_segmentation" / "legacy_growth_inputs" / dataset / "detection_images" / f"{fov}__first_frames_mean.tif"


def source_legacy_mask(cfg: dict[str, Any], dataset: str, fov: str) -> Path:
    path = source_dataset_root(cfg, dataset) / "03_roi_detection_and_growth" / "masks" / f"{fov}__roi_labels.tif"
    return path if path.exists() else roots(cfg)["technical"] / "exploratory_segmentation" / "legacy_growth_inputs" / dataset / "masks" / f"{fov}__roi_labels.tif"


class UnionFind:
    def __init__(self, values: Iterable[int]):
        self.parent = {int(value): int(value) for value in values}

    def find(self, value: int) -> int:
        value = int(value)
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> None:
        a, b = self.find(left), self.find(right)
        if a == b:
            return
        lo, hi = sorted((a, b))
        self.parent[hi] = lo


def circle_overlap_fraction(distance: float, radius: float) -> float:
    if not np.isfinite(distance) or distance >= 2 * radius:
        return 0.0
    if distance <= 0:
        return 1.0
    area = 2 * radius * radius * math.acos(distance / (2 * radius))
    area -= 0.5 * distance * math.sqrt(max(0.0, 4 * radius * radius - distance * distance))
    return float(area / (math.pi * radius * radius))


def build_seed_and_group_units(seeds: pd.DataFrame, fov: str, link_distance: float, radius: float) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    seeds = seeds.sort_values("seed_id").reset_index(drop=True).copy()
    coordinates = seeds[["y_px", "x_px"]].to_numpy(float)
    tree = cKDTree(coordinates) if len(coordinates) else None
    nearest = np.full(len(seeds), np.inf)
    pairs: set[tuple[int, int]] = set()
    if tree is not None and len(seeds) > 1:
        distances, _ = tree.query(coordinates, k=2)
        nearest = distances[:, 1]
        pairs = set(tree.query_pairs(float(link_distance)))
    union = UnionFind(seeds.seed_id.astype(int))
    for left_index, right_index in sorted(pairs):
        union.union(int(seeds.iloc[left_index].seed_id), int(seeds.iloc[right_index].seed_id))
    roots_by_seed = {int(seed_id): union.find(int(seed_id)) for seed_id in seeds.seed_id}
    root_order = {root: index + 1 for index, root in enumerate(sorted(set(roots_by_seed.values())))}
    seeds["nearest_neighbor_distance_px"] = nearest
    seeds["crowded_seed"] = nearest < float(link_distance)
    seeds["maximum_signal_overlap_fraction"] = [circle_overlap_fraction(value, radius) for value in nearest]
    seeds["seed_level_unit_id"] = [f"{fov}__seed_{int(value):05d}" for value in seeds.seed_id]
    seeds["group_index"] = [root_order[roots_by_seed[int(value)]] for value in seeds.seed_id]
    seeds["grouped_unit_id"] = [f"{fov}__group_{int(value):05d}" for value in seeds.group_index]
    groups = []
    for grouped_unit_id, group in seeds.groupby("grouped_unit_id", sort=True):
        ids = sorted(group.seed_id.astype(int).tolist())
        groups.append({
            "analysis_unit_id": grouped_unit_id,
            "analysis_stream": "grouped_cluster",
            "member_seed_ids": ";".join(map(str, ids)),
            "member_seed_count": len(ids),
            "centroid_y_px": float(group.y_px.mean()),
            "centroid_x_px": float(group.x_px.mean()),
            "crowded_group": len(ids) > 1,
        })
    group_table = pd.DataFrame.from_records(groups)
    seed_units = pd.DataFrame({
        "analysis_unit_id": seeds.seed_level_unit_id,
        "analysis_stream": "seed_level_cluster",
        "member_seed_ids": seeds.seed_id.astype(int).astype(str),
        "member_seed_count": 1,
        "centroid_y_px": seeds.y_px,
        "centroid_x_px": seeds.x_px,
        "crowded_group": seeds.crowded_seed,
    })
    mapping = seeds[["seed_id", "seed_level_unit_id", "grouped_unit_id", "group_index"]].copy()
    return seeds, seed_units, group_table, mapping


def extended_watershed(cfg: dict[str, Any], dataset: str, fov: str, seeds: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    image = tifffile.imread(source_detection_image(cfg, dataset, fov)).astype(np.float32)
    legacy = tifffile.imread(source_legacy_mask(cfg, dataset, fov)) > 0
    successful = seeds[seeds.legacy_growth_status.isin(["grown", "grown_with_warning"])].copy()
    markers = np.zeros(image.shape, np.int32)
    marker_rows = []
    occupied: set[tuple[int, int]] = set()
    next_marker = 1
    for row in successful.sort_values(["spotiflow_probability", "seed_id"], ascending=[False, True]).itertuples():
        y, x = int(round(row.y_px)), int(round(row.x_px))
        if not (0 <= y < image.shape[0] and 0 <= x < image.shape[1]) or (y, x) in occupied or not legacy[y, x]:
            continue
        occupied.add((y, x)); markers[y, x] = next_marker
        marker_rows.append({"marker": next_marker, "seed_id": int(row.seed_id), "source_growth_status": row.legacy_growth_status})
        next_marker += 1
    if next_marker == 1:
        return (np.zeros_like(markers, dtype=np.uint16), pd.DataFrame(columns=EXTENDED_COLUMNS),
                pd.DataFrame(columns=["seed_id", "extended_unit_id", "extended_segmentation_status", "extended_rejection_reason"]))
    labels = watershed(-gaussian_filter(image, 1.0), markers=markers, mask=legacy, watershed_line=True)
    foreground = labels > 0
    settings = cfg["analysis_units"]
    gap = float(settings["extended_background_gap_px"]); width = float(settings["extended_background_width_px"])
    minimum = int(settings.get("minimum_extended_area_px", 20)); maximum = int(settings["maximum_extended_area_px"]); minimum_bg = int(settings["minimum_background_pixels"])
    records, mappings = [], []
    for marker in marker_rows:
        label_id = marker["marker"]; mask = labels == label_id; area = int(mask.sum())
        if not area:
            continue
        yy, xx = np.where(mask); touches_image = bool((yy == 0).any() or (xx == 0).any() or (yy == image.shape[0] - 1).any() or (xx == image.shape[1] - 1).any())
        distance = distance_transform_edt(~mask)
        background = (distance > gap) & (distance <= gap + width) & (~foreground)
        reasons = []
        if marker["source_growth_status"] == "grown_with_warning": reasons.append("search_boundary_touched")
        if touches_image: reasons.append("image_boundary_touched")
        if area < minimum: reasons.append("not_extended")
        if area > maximum: reasons.append("cell_or_nucleus_scale")
        if int(background.sum()) < minimum_bg: reasons.append("insufficient_background_pixels")
        valid = not reasons
        unit_id = f"{fov}__extended_{label_id:05d}"
        records.append({
            "analysis_unit_id": unit_id, "analysis_stream": "extended_cluster_candidate",
            "source_seed_id": marker["seed_id"], "source_growth_status": marker["source_growth_status"],
            "area_px": area, "equivalent_diameter_px": float(2 * math.sqrt(area / math.pi)),
            "centroid_y_px": float(yy.mean()), "centroid_x_px": float(xx.mean()),
            "background_pixels": int(background.sum()), "extended_segmentation_valid": valid,
            "extended_rejection_reason": ";".join(reasons), "condensate_candidate": False,
            "condensate_validation_status": "not_manually_validated",
            "label_id": label_id,
        })
        mappings.append({"seed_id": marker["seed_id"], "extended_unit_id": unit_id if valid else "", "extended_segmentation_status": "valid" if valid else "rejected", "extended_rejection_reason": ";".join(reasons)})
    return (labels.astype(np.uint16), pd.DataFrame.from_records(records, columns=EXTENDED_COLUMNS),
            pd.DataFrame.from_records(mappings, columns=["seed_id", "extended_unit_id", "extended_segmentation_status", "extended_rejection_reason"]))


def _save_unit_overlay(cfg: dict[str, Any], dataset: str, fov: str, seeds: pd.DataFrame, labels: np.ndarray, extended: pd.DataFrame) -> Path:
    matplotlib_setup(cfg["output_root"])
    import matplotlib.pyplot as plt
    image = tifffile.imread(source_detection_image(cfg, dataset, fov)).astype(float)
    lo, hi = np.percentile(image, [1, 99.8]); display = np.clip((image - lo) / max(hi - lo, 1e-6), 0, 1)
    rgb = np.repeat(display[..., None], 3, axis=2)
    valid_ids = set(extended.loc[extended.extended_segmentation_valid.astype(bool), "label_id"].astype(int)) if not extended.empty else set()
    invalid_ids = set(extended.label_id.astype(int)) - valid_ids if not extended.empty else set()
    valid_boundary = find_boundaries(np.isin(labels, list(valid_ids)), mode="outer") if valid_ids else np.zeros(labels.shape, bool)
    invalid_boundary = find_boundaries(np.isin(labels, list(invalid_ids)), mode="outer") if invalid_ids else np.zeros(labels.shape, bool)
    rgb[valid_boundary] = (1, 1, 0); rgb[invalid_boundary] = (.45, .45, .45)
    fig, ax = plt.subplots(figsize=(7, 10)); ax.imshow(rgb)
    ax.scatter(seeds.x_px, seeds.y_px, s=9, facecolors="none", edgecolors="cyan", linewidths=.45)
    crowded = seeds[seeds.crowded_seed.astype(bool)]
    if len(crowded): ax.scatter(crowded.x_px, crowded.y_px, s=12, marker="+", c="magenta", linewidths=.5)
    ax.set_title(f"{fov}: all {len(seeds):,} seeds retained; {len(valid_ids)} valid extended candidates")
    ax.axis("off"); fig.tight_layout()
    path = roots(cfg)["figures"] / "method_QC" / f"{dataset}__{fov}__analysis_units.png"
    fig.savefig(path, dpi=220, bbox_inches="tight"); plt.close(fig)
    return path


def build_analysis_units(cfg: dict[str, Any], *, resume: bool = False) -> dict[str, Path]:
    paths = ensure_layout(cfg); seed_ledgers, mappings, sensitivity = [], [], []
    radius = float(cfg["analysis_units"]["signal_radius_px"]); link = float(cfg["analysis_units"]["grouped_link_distance_px"])
    for dataset_cfg in cfg["datasets"]:
        dataset = dataset_cfg["name"]; manifest = source_manifest(cfg, dataset)
        for row in manifest.to_dict("records"):
            fov = row["fov"]; out = experiment_root(cfg, dataset) / "spot_detection" / fov
            complete = out / f"{fov}__seed_to_analysis_unit_mapping.csv"
            if resume and complete.exists():
                seed_ledgers.append(pd.read_csv(out / f"{fov}__complete_seed_ledger.csv")); mappings.append(pd.read_csv(complete)); continue
            out.mkdir(parents=True, exist_ok=True)
            seeds = source_seed_table(cfg, dataset, fov)
            seeds, seed_units, groups, mapping = build_seed_and_group_units(seeds, fov, link, radius)
            labels, extended, extended_map = extended_watershed(cfg, dataset, fov, seeds)
            mapping = mapping.merge(extended_map, on="seed_id", how="left")
            mapping["extended_segmentation_status"] = mapping.extended_segmentation_status.fillna("not_successfully_grown")
            mapping["extended_rejection_reason"] = mapping.extended_rejection_reason.fillna("")
            mapping["extended_unit_id"] = mapping.extended_unit_id.fillna("")
            seeds["dataset"] = dataset; seeds["filename"] = row["filename"]; seeds["fov"] = fov
            seeds["acquisition_profile"] = row["acquisition_profile"]
            for key, value in condition_metadata(row["filename"]).items(): seeds[key] = value
            seeds["growth_is_morphology_only"] = True
            mapping.insert(0, "dataset", dataset); mapping.insert(1, "filename", row["filename"]); mapping.insert(2, "fov", fov)
            for table in (seed_units, groups, extended):
                if not table.empty:
                    table.insert(0, "dataset", dataset); table.insert(1, "filename", row["filename"]); table.insert(2, "fov", fov)
                    table["acquisition_profile"] = row["acquisition_profile"]
                    for key, value in condition_metadata(row["filename"]).items(): table[key] = value
            atomic_csv(seeds, out / f"{fov}__complete_seed_ledger.csv")
            atomic_csv(seed_units, out / f"{fov}__seed_level_units.csv")
            atomic_csv(groups, out / f"{fov}__grouped_cluster_units.csv")
            atomic_csv(extended, out / f"{fov}__extended_cluster_candidates.csv")
            atomic_csv(mapping, complete)
            tifffile.imwrite(out / f"{fov}__extended_watershed_labels.tif", labels, metadata={"axes": "YX"})
            _save_unit_overlay(cfg, dataset, fov, seeds, labels, extended)
            for threshold in [0.4] + list(cfg["spotiflow"]["sensitivity_thresholds"]):
                sensitivity.append({"dataset": dataset, "filename": row["filename"], "fov": fov, "spotiflow_threshold": threshold, "seed_count": int((seeds.spotiflow_probability >= threshold).sum())})
            seed_ledgers.append(seeds); mappings.append(mapping)
    ledger = pd.concat(seed_ledgers, ignore_index=True, sort=False); mapping_table = pd.concat(mappings, ignore_index=True, sort=False)
    atomic_csv(ledger, paths["tables"] / "complete_seed_ledger.csv")
    atomic_csv(mapping_table, paths["tables"] / "seed_to_analysis_unit_mapping.csv")
    if sensitivity:
        atomic_csv(pd.DataFrame.from_records(sensitivity), paths["technical"] / "sensitivity_analysis" / "spotiflow_detection_threshold_sensitivity.csv")
    return {"ledger": paths["tables"] / "complete_seed_ledger.csv", "mapping": paths["tables"] / "seed_to_analysis_unit_mapping.csv"}


def _disk_indices(shape: tuple[int, int], y: float, x: float, radius: float) -> np.ndarray:
    rr, cc = disk((float(y), float(x)), radius, shape=shape)
    return np.ravel_multi_index((rr, cc), shape)


def _annulus_indices(shape: tuple[int, int], y: float, x: float, inner: float, outer: float, excluded: np.ndarray) -> np.ndarray:
    y0 = max(0, int(math.floor(y - outer))); y1 = min(shape[0], int(math.ceil(y + outer)) + 1)
    x0 = max(0, int(math.floor(x - outer))); x1 = min(shape[1], int(math.ceil(x + outer)) + 1)
    yy, xx = np.ogrid[y0:y1, x0:x1]; distance = np.hypot(yy - y, xx - x)
    selector = (distance >= inner) & (distance <= outer) & (~excluded[y0:y1, x0:x1])
    rr, cc = np.where(selector); return np.ravel_multi_index((rr + y0, cc + x0), shape)


@dataclass
class TraceEntry:
    metadata: dict[str, Any]
    signal: np.ndarray
    background: np.ndarray
    status: str
    reason: str


def build_trace_entries(cfg: dict[str, Any], dataset: str, fov: str, stream: str) -> list[TraceEntry]:
    out = experiment_root(cfg, dataset) / "spot_detection" / fov
    seeds = pd.read_csv(out / f"{fov}__complete_seed_ledger.csv")
    stack = source_corrected_stack(cfg, dataset, fov)
    with tifffile.TiffFile(stack) as tif: shape = tif.pages[0].shape
    settings = cfg["analysis_units"]; radius = float(settings["signal_radius_px"])
    inner = float(settings["background_inner_radius_px"]); outer = float(settings["background_outer_radius_px"])
    min_fraction = float(settings["minimum_signal_fraction"]); minimum_bg = int(settings["minimum_background_pixels"])
    seed_signals = {int(row.seed_id): _disk_indices(shape, row.y_px, row.x_px, radius) for row in seeds.itertuples()}
    all_seed_foreground = np.zeros(shape[0] * shape[1], bool)
    for indices in seed_signals.values(): all_seed_foreground[indices] = True
    all_seed_foreground = all_seed_foreground.reshape(shape)
    expected_signal = len(_disk_indices((21, 21), 10, 10, radius))
    entries: list[TraceEntry] = []
    if stream == "seed_level_cluster":
        for row in seeds.itertuples():
            signal = seed_signals[int(row.seed_id)]
            background = _annulus_indices(shape, row.y_px, row.x_px, inner, outer, all_seed_foreground)
            reasons = []
            if len(signal) < math.ceil(min_fraction * expected_signal): reasons.append("insufficient_signal_pixels")
            if len(background) < minimum_bg: reasons.append("insufficient_background_pixels")
            entries.append(TraceEntry({"analysis_unit_id": row.seed_level_unit_id, "analysis_stream": stream,
                                       "member_seed_ids": str(int(row.seed_id)), "member_seed_count": 1,
                                       "crowded_seed": bool(row.crowded_seed), "nearest_neighbor_distance_px": row.nearest_neighbor_distance_px,
                                       "signal_overlap_fraction": row.maximum_signal_overlap_fraction}, signal, background,
                                      "trace_ready" if not reasons else "trace_rejected", ";".join(reasons)))
    elif stream == "grouped_cluster":
        groups = pd.read_csv(out / f"{fov}__grouped_cluster_units.csv")
        group_masks: dict[str, np.ndarray] = {}
        for row in groups.itertuples():
            member_ids = [int(value) for value in str(row.member_seed_ids).split(";")]
            mask = np.zeros(shape[0] * shape[1], bool)
            for seed_id in member_ids: mask[seed_signals[seed_id]] = True
            group_masks[row.analysis_unit_id] = mask.reshape(shape)
        group_foreground = np.logical_or.reduce(list(group_masks.values())) if group_masks else np.zeros(shape, bool)
        for row in groups.itertuples():
            mask = group_masks[row.analysis_unit_id]; signal = np.flatnonzero(mask)
            distance = distance_transform_edt(~mask)
            background = np.flatnonzero((distance > inner) & (distance <= outer) & (~group_foreground))
            reasons = []
            if len(signal) < math.ceil(min_fraction * expected_signal): reasons.append("insufficient_signal_pixels")
            if len(background) < minimum_bg: reasons.append("insufficient_background_pixels")
            entries.append(TraceEntry({"analysis_unit_id": row.analysis_unit_id, "analysis_stream": stream,
                                       "member_seed_ids": row.member_seed_ids, "member_seed_count": int(row.member_seed_count),
                                       "crowded_seed": bool(row.crowded_group)}, signal, background,
                                      "trace_ready" if not reasons else "trace_rejected", ";".join(reasons)))
    elif stream == "extended_cluster_candidate":
        units = pd.read_csv(out / f"{fov}__extended_cluster_candidates.csv")
        labels = tifffile.imread(out / f"{fov}__extended_watershed_labels.tif")
        valid = units[units.extended_segmentation_valid.astype(bool)].copy() if len(units) else units
        extended_foreground = labels > 0
        excluded = extended_foreground | all_seed_foreground
        gap = float(settings["extended_background_gap_px"]); width = float(settings["extended_background_width_px"])
        for row in valid.itertuples():
            mask = labels == int(row.label_id); signal = np.flatnonzero(mask); distance = distance_transform_edt(~mask)
            background = np.flatnonzero((distance > gap) & (distance <= gap + width) & (~excluded))
            reasons = []
            if len(background) < minimum_bg: reasons.append("insufficient_background_pixels")
            entries.append(TraceEntry({"analysis_unit_id": row.analysis_unit_id, "analysis_stream": stream,
                                       "member_seed_ids": str(int(row.source_seed_id)), "member_seed_count": 1,
                                       "crowded_seed": False, "area_px": int(row.area_px)}, signal, background,
                                      "trace_ready" if not reasons else "trace_rejected", ";".join(reasons)))
    else:
        raise ValueError(stream)
    return entries


def _weight_matrix(entries: list[TraceEntry], pixels: int, which: str) -> csr_matrix:
    rows, columns, values = [], [], []
    for row_index, entry in enumerate(entries):
        indices = entry.signal if which == "signal" else entry.background
        if not len(indices): continue
        rows.extend([row_index] * len(indices)); columns.extend(indices.tolist()); values.extend([1.0 / len(indices)] * len(indices))
    return csr_matrix((values, (rows, columns)), shape=(len(entries), pixels), dtype=np.float32)


def extract_entry_traces(stack: Path, entries: list[TraceEntry]) -> dict[str, np.ndarray]:
    ready = [entry for entry in entries if entry.status == "trace_ready"]
    if not ready:
        return {name: np.empty((0, 0), np.float32) for name in TRACE_REPRESENTATIONS}
    with tifffile.TiffFile(stack) as tif:
        frames = len(tif.pages); shape = tif.pages[0].shape; pixels = int(np.prod(shape))
        signal_weights = _weight_matrix(ready, pixels, "signal"); background_weights = _weight_matrix(ready, pixels, "background")
        signal_mean = np.empty((len(ready), frames), np.float32); background_mean = np.empty_like(signal_mean)
        for frame_index, page in enumerate(tif.pages):
            image = page.asarray().reshape(-1).astype(np.float32, copy=False)
            signal_mean[:, frame_index] = signal_weights @ image
            background_mean[:, frame_index] = background_weights @ image
    mean_difference = signal_mean - background_mean
    signal_pixels = np.asarray([len(entry.signal) for entry in ready], np.float32)[:, None]
    integrated = mean_difference * signal_pixels
    tail = min(500, max(20, frames // 10)); centered = integrated - np.median(integrated[:, -tail:], axis=1, keepdims=True)
    return {"signal_mean": signal_mean, "background_mean": background_mean, "mean_difference": mean_difference,
            "integrated_difference": integrated, "baseline_centered_integrated_difference": centered}


def _write_trace_table(metadata: pd.DataFrame, values: np.ndarray, path: Path) -> None:
    if metadata.empty and not len(metadata.columns):
        metadata = pd.DataFrame(columns=["dataset", "filename", "fov", "analysis_unit_id", "analysis_stream",
                                                "member_seed_ids", "member_seed_count", "signal_pixels", "background_pixels",
                                                "acquisition_profile", "condition", "analysis_set", "is_primary_comparison"])
    frames = pd.DataFrame(values, columns=[str(index) for index in range(values.shape[1])])
    atomic_csv(pd.concat([metadata.reset_index(drop=True), frames], axis=1), path)


def extract_all_traces(cfg: dict[str, Any], *, resume: bool = False) -> None:
    ensure_layout(cfg)
    for dataset_cfg in cfg["datasets"]:
        dataset = dataset_cfg["name"]; manifest = source_manifest(cfg, dataset)
        for row in manifest.to_dict("records"):
            fov = row["fov"]
            for stream in STREAMS:
                out = experiment_root(cfg, dataset) / "puncta_traces" / stream / fov
                qc_path = out / f"{fov}__trace_QC.csv"
                if resume and qc_path.exists() and qc_path.stat().st_size > 1: continue
                out.mkdir(parents=True, exist_ok=True)
                entries = build_trace_entries(cfg, dataset, fov, stream)
                traces = extract_entry_traces(source_corrected_stack(cfg, dataset, fov), entries)
                ready = [entry for entry in entries if entry.status == "trace_ready"]
                metadata = pd.DataFrame.from_records([entry.metadata for entry in ready])
                if len(metadata):
                    metadata.insert(0, "dataset", dataset); metadata.insert(1, "filename", row["filename"]); metadata.insert(2, "fov", fov)
                    metadata["signal_pixels"] = [len(entry.signal) for entry in ready]
                    metadata["background_pixels"] = [len(entry.background) for entry in ready]
                    metadata["acquisition_profile"] = row["acquisition_profile"]
                    for key, value in condition_metadata(row["filename"]).items(): metadata[key] = value
                for representation, values in traces.items(): _write_trace_table(metadata, values, out / f"{fov}__{representation}.csv")
                qc_records = []
                ready_index = 0
                for entry in entries:
                    record = {**entry.metadata, "dataset": dataset, "filename": row["filename"], "fov": fov,
                              "acquisition_profile": row["acquisition_profile"], "trace_status": entry.status,
                              "trace_rejection_reason": entry.reason, "signal_pixels": len(entry.signal), "background_pixels": len(entry.background)}
                    if entry.status == "trace_ready":
                        for key, value in condition_metadata(row["filename"]).items(): record[key] = value
                        values = traces["integrated_difference"][ready_index]; mean_values = traces["mean_difference"][ready_index]
                        tail = min(500, max(20, len(values) // 10))
                        record.update({"late_integrated_baseline": float(np.median(values[-tail:])),
                                       "late_mean_baseline": float(np.median(mean_values[-tail:])),
                                       "late_integrated_noise": float(1.4826 * np.median(np.abs(np.diff(values[-tail:]) - np.median(np.diff(values[-tail:]))))) if tail > 2 else np.nan})
                        ready_index += 1
                    qc_records.append(record)
                qc_columns = ["analysis_unit_id", "analysis_stream", "member_seed_ids", "member_seed_count", "crowded_seed",
                              "dataset", "filename", "fov", "acquisition_profile", "trace_status", "trace_rejection_reason",
                              "signal_pixels", "background_pixels", "condition", "analysis_set", "is_primary_comparison",
                              "late_integrated_baseline", "late_mean_baseline", "late_integrated_noise"]
                atomic_csv(pd.DataFrame.from_records(qc_records).reindex(columns=qc_columns), qc_path)


def verify_unit_invariants(cfg: dict[str, Any]) -> dict[str, int]:
    ledger = pd.read_csv(roots(cfg)["tables"] / "complete_seed_ledger.csv")
    mapping = pd.read_csv(roots(cfg)["tables"] / "seed_to_analysis_unit_mapping.csv")
    keys = ["dataset", "fov", "seed_id"]
    if ledger.duplicated(keys).any() or mapping.duplicated(keys).any(): raise AssertionError("Seed identifiers are not unique")
    if len(ledger) != len(mapping): raise AssertionError("Every seed must map to both co-primary streams")
    if mapping.seed_level_unit_id.isna().any() or mapping.grouped_unit_id.isna().any(): raise AssertionError("Missing co-primary analysis-unit mapping")
    failed = ledger.legacy_growth_status.isin(["seed_below_local_threshold", "region_too_small"])
    if ledger.loc[failed, "seed_level_unit_id"].isna().any(): raise AssertionError("Failed-growth seed was discarded")
    grouped_units = mapping[["dataset", "fov", "grouped_unit_id"]].drop_duplicates()
    return {"seeds": len(ledger), "failed_growth_seeds_retained": int(failed.sum()), "grouped_units": len(grouped_units)}
