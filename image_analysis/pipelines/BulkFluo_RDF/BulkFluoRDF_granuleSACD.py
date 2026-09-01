#!/usr/bin/env python
"""Four-channel SACD granule segmentation and size-normalized fluorescence RDF.

The 405 channel is retained for exported crops.  Granules are segmented from an
equal-weight aggregate of independently normalized 488/561/647 images and must
be supported by at least two RNA channels.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import os
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import tifffile
import yaml
from scipy import ndimage
from skimage import feature, filters, measure, morphology, segmentation


CHANNELS = ("405", "488", "561", "647")
RNA_CHANNELS = ("488", "561", "647")
FILE_RE = re.compile(r"^(?P<fov>.+)__(?P<channel>405|488|561|647)-SACD-MIP-YX\.tiff?$", re.I)

DEFAULT_CONFIG = {
    "input_dir": ".",
    "output_dir": "BulkFluoRDF_granuleSACD_results",
    "expected_fovs": 10,
    "pixel_size_nm": 58.5,
    "crop_padding_px": 10,
    "segmentation": {
        "background_sigma_px": 12.0,
        "smooth_sigma_px": 1.5,
        "normalization_upper_percentile": 99.7,
        "minimum_raw_area_px": 2,
        "watershed_min_distance_px": 5,
        "minimum_support_channels": 2,
        "minimum_channel_support_fraction": 0.10,
        "minimum_equivalent_diameter_px": 15.0,
    },
    "rdf": {"maximum_normalized_radius": 1.3, "bin_width": 0.10, "bin_step": 0.05,
            "subpixel_samples_per_axis": 8},
    "qc": {"contrast_upper_percentile": 99.7, "make_overlays": True},
}


@dataclass(frozen=True)
class FourChannelFOV:
    fov: str
    paths: dict[str, Path]


@dataclass
class PipelineResult:
    pairs: list[FourChannelFOV]
    properties: pd.DataFrame
    rdf: pd.DataFrame
    aggregate: pd.DataFrame
    rdf_wide: pd.DataFrame
    correlations: pd.DataFrame
    output_dir: Path


def _deep_update(base: dict, update: dict) -> dict:
    result = copy.deepcopy(base)
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_update(result[key], value)
        else:
            result[key] = value
    return result


def load_config(config: str | Path | dict) -> dict:
    if isinstance(config, dict):
        raw = config
    else:
        with Path(config).open() as handle:
            raw = yaml.safe_load(handle) or {}
    cfg = _deep_update(DEFAULT_CONFIG, raw)
    cfg["input_dir"] = str(Path(cfg["input_dir"]).expanduser().resolve())
    output = Path(cfg["output_dir"]).expanduser()
    if not output.is_absolute():
        output = Path(cfg["input_dir"]) / output
    cfg["output_dir"] = str(output.resolve())
    return cfg


def pair_four_channel_files(input_dir: str | Path) -> list[FourChannelFOV]:
    root = Path(input_dir)
    found: dict[str, dict[str, Path]] = {}
    duplicates: list[str] = []
    for path in sorted(root.glob("*SACD-MIP-YX.tif*")):
        match = FILE_RE.match(path.name)
        if not match:
            continue
        fov, channel = match.group("fov"), match.group("channel")
        if channel in found.setdefault(fov, {}):
            duplicates.append(f"{fov}:{channel}")
        found[fov][channel] = path.resolve()
    if duplicates:
        raise ValueError(f"Duplicate FOV/channel inputs: {duplicates}")
    incomplete = {fov: sorted(set(CHANNELS) - set(paths)) for fov, paths in found.items() if set(paths) != set(CHANNELS)}
    if incomplete:
        raise ValueError(f"Incomplete four-channel FOVs: {incomplete}")
    if not found:
        raise FileNotFoundError(f"No four-channel SACD-MIP-YX inputs found in {root}")
    return [FourChannelFOV(fov, {ch: found[fov][ch] for ch in CHANNELS}) for fov in sorted(found)]


def read_fov(pair: FourChannelFOV) -> dict[str, np.ndarray]:
    images = {ch: np.asarray(tifffile.imread(pair.paths[ch]), dtype=np.float32) for ch in CHANNELS}
    shapes = {image.shape for image in images.values()}
    if len(shapes) != 1 or any(image.ndim != 2 for image in images.values()):
        raise ValueError(f"{pair.fov}: channels are not aligned 2D images: {[x.shape for x in images.values()]}")
    if not all(np.isfinite(image).all() for image in images.values()):
        raise ValueError(f"{pair.fov}: non-finite input intensities")
    return images


def background_correct_and_normalize(image: np.ndarray, background_sigma_px: float, upper_percentile: float) -> np.ndarray:
    corrected = np.maximum(image.astype(np.float32) - ndimage.gaussian_filter(image, background_sigma_px), 0)
    scale = float(np.percentile(corrected, upper_percentile))
    if not np.isfinite(scale) or scale <= 0:
        return np.zeros_like(corrected)
    return np.clip(corrected / scale, 0, 1).astype(np.float32)


def _otsu_or_inf(image: np.ndarray) -> float:
    finite = image[np.isfinite(image)]
    return float(filters.threshold_otsu(finite)) if finite.size and np.ptp(finite) > 0 else float("inf")


def segment_granules(images: dict[str, np.ndarray], config: dict) -> tuple[np.ndarray, np.ndarray, pd.DataFrame, dict]:
    cfg = config["segmentation"]
    normalized = {
        ch: background_correct_and_normalize(
            images[ch], float(cfg["background_sigma_px"]), float(cfg["normalization_upper_percentile"])
        )
        for ch in RNA_CHANNELS
    }
    smoothed = {ch: ndimage.gaussian_filter(normalized[ch], float(cfg["smooth_sigma_px"])) for ch in RNA_CHANNELS}
    aggregate = np.mean(np.stack([smoothed[ch] for ch in RNA_CHANNELS]), axis=0).astype(np.float32)
    aggregate_threshold = _otsu_or_inf(aggregate)
    channel_thresholds = {ch: _otsu_or_inf(smoothed[ch]) for ch in RNA_CHANNELS}
    foreground = aggregate > aggregate_threshold
    foreground = morphology.remove_small_objects(foreground, min_size=int(cfg["minimum_raw_area_px"]))
    foreground = ndimage.binary_fill_holes(foreground)
    distance = ndimage.distance_transform_edt(foreground)
    coordinates = feature.peak_local_max(
        distance, min_distance=int(cfg["watershed_min_distance_px"]), labels=foreground, exclude_border=False
    )
    markers = np.zeros(foreground.shape, dtype=np.int32)
    for index, (y, x) in enumerate(coordinates, 1):
        markers[y, x] = index
    candidates = (
        segmentation.watershed(-distance, markers, mask=foreground).astype(np.int32)
        if len(coordinates)
        else measure.label(foreground).astype(np.int32)
    )

    retained = np.zeros_like(candidates, dtype=np.uint16)
    records: list[dict] = []
    next_id = 1
    h, w = foreground.shape
    maximum_rho = float(config["rdf"]["maximum_normalized_radius"])
    min_support = int(cfg["minimum_support_channels"])
    min_fraction = float(cfg["minimum_channel_support_fraction"])
    min_diameter = float(cfg["minimum_equivalent_diameter_px"])
    for region in measure.regionprops(candidates):
        mask = candidates == region.label
        fractions = {ch: float(np.mean(smoothed[ch][mask] > channel_thresholds[ch])) for ch in RNA_CHANNELS}
        supported = [ch for ch, fraction in fractions.items() if fraction >= min_fraction]
        cy, cx = map(float, region.centroid)
        area = int(region.area)
        equivalent_radius = math.sqrt(area / math.pi)
        equivalent_diameter = 2 * equivalent_radius
        touches_edge = bool(region.bbox[0] == 0 or region.bbox[1] == 0 or region.bbox[2] == h or region.bbox[3] == w)
        analysis_radius = maximum_rho * equivalent_radius
        complete_rdf = bool(cy - analysis_radius >= 0 and cx - analysis_radius >= 0 and cy + analysis_radius < h and cx + analysis_radius < w)
        reasons = []
        if len(supported) < min_support:
            reasons.append("insufficient_channel_support")
        if not equivalent_diameter > min_diameter:
            reasons.append("below_minimum_granule_size")
        if touches_edge:
            reasons.append("touches_image_edge")
        if not complete_rdf:
            reasons.append("incomplete_rdf_radius")
        keep = not reasons
        granule_id = next_id if keep else 0
        if keep:
            retained[mask] = next_id
            next_id += 1
        records.append(
            {
                "candidate_id": int(region.label), "granule_id": int(granule_id), "keep_granule": keep,
                "rejection_reason": ";".join(reasons), "area_px": area,
                "equivalent_radius_px": equivalent_radius, "equivalent_diameter_px": equivalent_diameter,
                "centroid_y_px": cy, "centroid_x_px": cx,
                "bbox_min_y": int(region.bbox[0]), "bbox_min_x": int(region.bbox[1]),
                "bbox_max_y": int(region.bbox[2]), "bbox_max_x": int(region.bbox[3]),
                "touches_image_edge": touches_edge, "complete_rdf_radius": complete_rdf,
                "n_support_channels": len(supported), "support_channels": ",".join(supported),
                **{f"support_fraction_{ch}": fractions[ch] for ch in RNA_CHANNELS},
            }
        )
    diagnostics = {
        "aggregate_threshold": aggregate_threshold,
        **{f"threshold_{ch}": channel_thresholds[ch] for ch in RNA_CHANNELS},
        "normalized": normalized,
    }
    return retained, candidates, pd.DataFrame.from_records(records), {"aggregate": aggregate, **diagnostics}


def rdf_bins(config: dict) -> list[tuple[float, float]]:
    cfg = config["rdf"]
    starts = np.arange(0, float(cfg["maximum_normalized_radius"]) - float(cfg["bin_width"]) + 1e-9,
                       float(cfg["bin_step"]))
    return [(float(start), float(start + float(cfg["bin_width"]))) for start in starts]


def calculate_granule_rdf(fov: str, images: dict[str, np.ndarray], labels: np.ndarray, properties: pd.DataFrame,
                          config: dict) -> pd.DataFrame:
    rows: list[dict] = []
    kept = properties[properties["keep_granule"]]
    for prop in kept.itertuples():
        mask = labels == int(prop.granule_id)
        radius = float(prop.equivalent_radius_px)
        maximum_radius = float(config["rdf"]["maximum_normalized_radius"]) * radius
        margin = maximum_radius + math.sqrt(2) / 2
        y0 = max(0, int(math.floor(float(prop.centroid_y_px) - margin)))
        y1 = min(labels.shape[0], int(math.ceil(float(prop.centroid_y_px) + margin)) + 1)
        x0 = max(0, int(math.floor(float(prop.centroid_x_px) - margin)))
        x1 = min(labels.shape[1], int(math.ceil(float(prop.centroid_x_px) + margin)) + 1)
        local_y, local_x = np.indices((y1 - y0, x1 - x0), dtype=np.float32)
        local_y += y0
        local_x += x0
        samples = int(config["rdf"].get("subpixel_samples_per_axis", 8))
        offsets = (np.arange(samples, dtype=np.float32) + 0.5) / samples - 0.5
        subpixel_rho = np.hypot(
            local_y[None, None, :, :] + offsets[:, None, None, None] - float(prop.centroid_y_px),
            local_x[None, None, :, :] + offsets[None, :, None, None] - float(prop.centroid_x_px),
        ) / radius
        for channel in RNA_CHANNELS:
            inside_mean = float(np.mean(images[channel][mask]))
            local_image = images[channel][y0:y1, x0:x1]
            for start, end in rdf_bins(config):
                weights = np.mean((subpixel_rho >= start) & (subpixel_rho < end), axis=(0, 1))
                effective_area = float(weights.sum())
                radial_mean = float(np.sum(local_image * weights) / effective_area) if effective_area > 0 else float("nan")
                normalized = radial_mean / inside_mean if effective_area > 0 and inside_mean > 0 else float("nan")
                rows.append({
                    "fov": fov, "granule_id": int(prop.granule_id), "channel": channel,
                    "radius_start_r_over_R": start, "radius_end_r_over_R": end,
                    "radius_mid_r_over_R": (start + end) / 2, "effective_pixel_area": effective_area,
                    "inside_granule_mean": inside_mean, "annular_mean": radial_mean,
                    "rdf_normalized": normalized,
                })
    return pd.DataFrame.from_records(rows)


def aggregate_rdf(rdf: pd.DataFrame) -> pd.DataFrame:
    columns = ["channel", "radius_start_r_over_R", "radius_end_r_over_R", "radius_mid_r_over_R",
               "rdf_mean", "rdf_std", "rdf_sem", "n_granules", "n_fovs"]
    if rdf.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    group_cols = ["channel", "radius_start_r_over_R", "radius_end_r_over_R", "radius_mid_r_over_R"]
    for keys, group in rdf.groupby(group_cols, sort=True):
        valid = group[np.isfinite(group["rdf_normalized"])].copy()
        values = valid["rdf_normalized"].to_numpy(float)
        n = len(values)
        std = float(np.std(values, ddof=1)) if n > 1 else 0.0 if n == 1 else float("nan")
        rows.append(dict(zip(group_cols, keys)) | {
            "rdf_mean": float(np.mean(values)) if n else float("nan"), "rdf_std": std,
            "rdf_sem": std / math.sqrt(n) if n else float("nan"), "n_granules": n,
            "n_fovs": int(valid["fov"].nunique()),
        })
    return pd.DataFrame.from_records(rows, columns=columns)


def rdf_to_wide(rdf: pd.DataFrame) -> pd.DataFrame:
    """Return one row per granule/bin with raw and normalized values for all RNA channels."""
    index = ["fov", "granule_id", "radius_start_r_over_R", "radius_end_r_over_R", "radius_mid_r_over_R"]
    if rdf.empty:
        columns = index + ["effective_pixel_area"]
        for channel in RNA_CHANNELS:
            columns += [f"inside_granule_mean_{channel}", f"annular_mean_{channel}", f"rdf_normalized_{channel}"]
        return pd.DataFrame(columns=columns)
    base = rdf[rdf["channel"] == RNA_CHANNELS[0]][index + ["effective_pixel_area"]].copy()
    for channel in RNA_CHANNELS:
        current = rdf[rdf["channel"] == channel][index + ["inside_granule_mean", "annular_mean", "rdf_normalized"]].copy()
        current = current.rename(columns={
            "inside_granule_mean": f"inside_granule_mean_{channel}",
            "annular_mean": f"annular_mean_{channel}",
            "rdf_normalized": f"rdf_normalized_{channel}",
        })
        base = base.merge(current, on=index, how="inner", validate="one_to_one")
    return base.sort_values(index).reset_index(drop=True)


CORRELATION_PAIRS = (("488", "561"), ("488", "647"), ("561", "647"))


def calculate_rdf_correlations(rdf_wide: pd.DataFrame, properties: pd.DataFrame) -> pd.DataFrame:
    columns = ["fov", "granule_id", "equivalent_diameter_px", "equivalent_diameter_nm", "n_bins"] + [
        f"pearson_r_{left}_{right}" for left, right in CORRELATION_PAIRS
    ]
    if rdf_wide.empty:
        return pd.DataFrame(columns=columns)
    diameters = properties[properties["keep_granule"]][
        ["fov", "granule_id", "equivalent_diameter_px", "equivalent_diameter_nm"]
    ]
    rows = []
    for (fov, granule_id), group in rdf_wide.groupby(["fov", "granule_id"], sort=True):
        group = group.sort_values("radius_mid_r_over_R")
        row = {"fov": fov, "granule_id": int(granule_id), "n_bins": len(group)}
        for left, right in CORRELATION_PAIRS:
            x = group[f"rdf_normalized_{left}"].to_numpy(float)
            y = group[f"rdf_normalized_{right}"].to_numpy(float)
            row[f"pearson_r_{left}_{right}"] = float(np.corrcoef(x, y)[0, 1])
        rows.append(row)
    result = pd.DataFrame.from_records(rows).merge(diameters, on=["fov", "granule_id"], validate="one_to_one")
    return result[columns].sort_values(["fov", "granule_id"]).reset_index(drop=True)


def export_granule_crops(pair: FourChannelFOV, images: dict[str, np.ndarray], properties: pd.DataFrame,
                         output_dir: Path, padding: int) -> list[Path]:
    paths = []
    h, w = next(iter(images.values())).shape
    for prop in properties[properties["keep_granule"]].itertuples():
        y0, x0 = max(0, int(prop.bbox_min_y) - padding), max(0, int(prop.bbox_min_x) - padding)
        y1, x1 = min(h, int(prop.bbox_max_y) + padding), min(w, int(prop.bbox_max_x) + padding)
        stack = np.stack([images[ch][y0:y1, x0:x1] for ch in CHANNELS]).astype(np.float32)
        path = output_dir / "granule_crops" / f"{pair.fov}__granule-{int(prop.granule_id):03d}.tif"
        path.parent.mkdir(parents=True, exist_ok=True)
        tifffile.imwrite(path, stack, imagej=True, metadata={"axes": "CYX", "Labels": list(CHANNELS), "unit": "pixel"})
        paths.append(path)
    return paths


def _setup_matplotlib(output_dir: Path):
    cache = output_dir / ".cache" / "matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib.pyplot as plt
    return plt


def save_qc_overlay(pair: FourChannelFOV, aggregate: np.ndarray, retained: np.ndarray, candidates: np.ndarray,
                    properties: pd.DataFrame, output_dir: Path, upper_percentile: float) -> Path:
    plt = _setup_matplotlib(output_dir)
    high = float(np.percentile(aggregate, upper_percentile))
    display = np.clip(aggregate / high, 0, 1) if high > 0 else np.zeros_like(aggregate)
    rgb = np.repeat(display[..., None], 3, axis=-1)
    kept_boundary = segmentation.find_boundaries(retained, mode="outer")
    rejected_ids = properties.loc[~properties["keep_granule"], "candidate_id"].astype(int).tolist()
    rejected_mask = np.isin(candidates, rejected_ids)
    rejected_boundary = segmentation.find_boundaries(rejected_mask, mode="outer")
    rgb[rejected_boundary] = (1, 0.75, 0)
    rgb[kept_boundary] = (1, 0, 1)
    fig, ax = plt.subplots(figsize=(7, 10))
    ax.imshow(rgb)
    for row in properties[properties["keep_granule"]].itertuples():
        ax.text(row.centroid_x_px, row.centroid_y_px, str(int(row.granule_id)), color="cyan", fontsize=5,
                ha="center", va="center")
    ax.set_title(f"{pair.fov}: retained {int(properties.keep_granule.sum())}/{len(properties)}\nmagenta retained; yellow rejected")
    ax.axis("off")
    path = output_dir / "qc_overlays" / f"{pair.fov}__granule_segmentation.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def _primary_rejection_reason(reason: str) -> str:
    reasons = set(str(reason).split(";"))
    for candidate in ("below_minimum_granule_size", "insufficient_channel_support", "touches_image_edge", "incomplete_rdf_radius"):
        if candidate in reasons:
            return candidate
    return "other"


def _significance_label(result) -> str:
    """Format a statannotations result as stars plus p, or only ns after correction."""
    if not result.is_significant:
        return "ns"
    pvalue = float(result.pvalue)
    stars = "****" if pvalue <= 1e-4 else "***" if pvalue <= 1e-3 else "**" if pvalue <= 1e-2 else "*"
    return f"{stars}\np={pvalue:.2e}"


def _save_correlation_violin_box(correlation_long: pd.DataFrame, order: list[str], pairs: list[tuple[str, str]],
                                 output_path: Path, test: str, title: str) -> None:
    import seaborn as sns
    from statannotations.Annotator import Annotator

    plt = _setup_matplotlib(output_path.parent)
    fig, ax = plt.subplots(figsize=(3, 3))
    sns.violinplot(data=correlation_long, x="channel_pair", y="pearson_r", order=order,
                   inner=None, cut=0, color="0.78", linewidth=0.8, ax=ax)
    sns.boxplot(data=correlation_long, x="channel_pair", y="pearson_r", order=order, width=0.22,
                showfliers=False, linewidth=0.8, boxprops={"facecolor": "white", "zorder": 3},
                medianprops={"color": "black", "linewidth": 1}, ax=ax)
    ax.set_ylim(-0.75, 1.50)
    annotator = Annotator(ax, pairs, data=correlation_long, x="channel_pair", y="pearson_r", order=order)
    annotator.configure(test=test, comparisons_correction="Benjamini-Hochberg", text_format="star",
                        loc="inside", fontsize=5.5, line_height=0.015, text_offset=0, verbose=0)
    annotator.apply_test()
    label_by_pair = {
        frozenset(struct["label"] for struct in annotation.structs): _significance_label(annotation.data)
        for annotation in annotator.annotations
    }
    labels = [label_by_pair[frozenset(pair)] for pair in pairs]
    annotator.set_custom_annotations(labels).annotate(line_offset_to_group=0.02)
    ax.axhline(0, color="0.6", linestyle="--", linewidth=0.7)
    ax.set(xlabel="RDF channel pair", ylabel="Pearson r", title=title)
    ax.set_title(title, fontsize=7, pad=3)
    ax.set_xlabel("RDF channel pair", fontsize=7)
    ax.set_ylabel("Pearson r", fontsize=7)
    ax.tick_params(axis="both", labelsize=6)
    fig.tight_layout(pad=0.5)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_summary_plots(properties: pd.DataFrame, rdf: pd.DataFrame, aggregate: pd.DataFrame,
                       correlations: pd.DataFrame, output_dir: Path, config: dict) -> None:
    plt = _setup_matplotlib(output_dir)
    rdf = rdf.copy()
    aggregate = aggregate.copy()
    rdf["channel"] = rdf["channel"].astype(str)
    aggregate["channel"] = aggregate["channel"].astype(str)
    cutoff = float(config["segmentation"]["minimum_equivalent_diameter_px"])
    pixel_size_nm = float(config["pixel_size_nm"])
    cutoff_um = cutoff * pixel_size_nm / 1000.0
    properties = properties.copy()
    properties["equivalent_diameter_um"] = properties["equivalent_diameter_px"] * pixel_size_nm / 1000.0
    fig, ax = plt.subplots(figsize=(7, 4))
    all_min = float(properties["equivalent_diameter_um"].min())
    all_max = float(properties["equivalent_diameter_um"].max())
    all_bins = np.linspace(all_min, all_max, 41)
    for keep, group in properties.groupby("keep_granule"):
        ax.hist(group["equivalent_diameter_um"], bins=all_bins, alpha=0.6,
                label="retained" if keep else "rejected")
    ax.axvline(cutoff_um, color="black", linestyle="--", label=f"size cutoff >{cutoff_um:.3f} µm")
    ax.set(xlabel="Equivalent diameter (µm)", ylabel="Candidates", title="Granule size distribution",
           xlim=(all_min, all_max))
    ax.margins(x=0)
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(output_dir / "granule_size_distribution.png", dpi=200); plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharex=True, sharey=True)
    for ax, channel in zip(axes, RNA_CHANNELS):
        for (_, _), curve in rdf[rdf.channel == channel].groupby(["fov", "granule_id"]):
            ax.plot(curve.radius_mid_r_over_R, curve.rdf_normalized, color="0.75", alpha=0.2, lw=0.6)
        curve = aggregate[aggregate.channel == channel]
        if not curve.empty:
            ax.errorbar(curve.radius_mid_r_over_R, curve.rdf_mean, yerr=curve.rdf_sem, color="black", lw=1.6)
        ax.axvline(1, color="magenta", linestyle="--", lw=1)
        ax.set(title=channel, xlabel="Distance from centroid (r/R)")
    axes[0].set_ylabel("Annular mean / granule mean")
    fig.suptitle("Size-normalized granule RDF")
    fig.tight_layout(); fig.savefig(output_dir / "rdf_aggregate.png", dpi=220); plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharex=True, sharey=True)
    for ax, (left, right) in zip(axes, CORRELATION_PAIRS):
        column = f"pearson_r_{left}_{right}"
        ax.hist(correlations[column], bins=np.linspace(-1, 1, 31), color="0.3", alpha=0.8)
        ax.axvline(float(correlations[column].mean()), color="magenta", linestyle="--", label="mean")
        ax.set(title=f"{left}–{right}", xlabel="Pearson r", xlim=(-1, 1))
        ax.legend(frameon=False)
    axes[0].set_ylabel("Granules")
    fig.suptitle("Within-granule RDF correlation distributions (0–1.3 R)")
    fig.tight_layout(); fig.savefig(output_dir / "granule_rdf_correlation_distributions.png", dpi=220); plt.close(fig)

    order = [f"{left}–{right}" for left, right in CORRELATION_PAIRS]
    correlation_long = correlations.melt(
        id_vars=["fov", "granule_id"],
        value_vars=[f"pearson_r_{left}_{right}" for left, right in CORRELATION_PAIRS],
        var_name="channel_pair", value_name="pearson_r",
    )
    correlation_long["channel_pair"] = correlation_long["channel_pair"].map({
        f"pearson_r_{left}_{right}": f"{left}–{right}" for left, right in CORRELATION_PAIRS
    })
    pairs = [(order[0], order[1]), (order[0], order[2]), (order[1], order[2])]
    _save_correlation_violin_box(
        correlation_long, order, pairs,
        output_dir / "granule_rdf_correlation_violin_box_statannotations.png",
        "Wilcoxon", "Paired Wilcoxon (BH correction)",
    )
    _save_correlation_violin_box(
        correlation_long, order, pairs,
        output_dir / "granule_rdf_correlation_violin_box_statannotations_paired_ttest.png",
        "t-test_paired", "Paired t-test (BH correction)",
    )


def export_correlation_representatives(correlations: pd.DataFrame, rdf_wide: pd.DataFrame,
                                       output_dir: Path) -> pd.DataFrame:
    """Export true pairwise max/min selections; duplicate biological granules are allowed."""
    plt = _setup_matplotlib(output_dir)
    representative_root = output_dir / "correlation_representatives"
    representative_root.mkdir(parents=True, exist_ok=True)
    rows = []
    for left, right in CORRELATION_PAIRS:
        score_column = f"pearson_r_{left}_{right}"
        for extreme in ("max", "min"):
            ascending = extreme == "min"
            selected = correlations.sort_values(
                [score_column, "fov", "granule_id"], ascending=[ascending, True, True], kind="mergesort"
            ).iloc[0]
            slot = f"{extreme}_{left}-{right}"
            slot_dir = representative_root / slot
            slot_dir.mkdir(parents=True, exist_ok=True)
            fov, granule_id = str(selected.fov), int(selected.granule_id)
            source_crop = output_dir / "granule_crops" / f"{fov}__granule-{granule_id:03d}.tif"
            destination_crop = slot_dir / f"{slot}__{fov}__granule-{granule_id:03d}.tif"
            shutil.copy2(source_crop, destination_crop)
            curves = rdf_wide[(rdf_wide["fov"] == fov) & (rdf_wide["granule_id"] == granule_id)].sort_values(
                "radius_mid_r_over_R"
            )
            fig, ax = plt.subplots(figsize=(7, 5))
            colors = {"488": "#00a651", "561": "#f28e2b", "647": "#d62728"}
            for channel in RNA_CHANNELS:
                ax.plot(curves["radius_mid_r_over_R"], curves[f"rdf_normalized_{channel}"],
                        lw=2, color=colors[channel], label=channel)
            ax.axvline(1, color="0.35", linestyle="--", lw=1, label="r/R = 1")
            score_text = "\n".join(
                f"{a}–{b}: r = {float(selected[f'pearson_r_{a}_{b}']):.3f}" for a, b in CORRELATION_PAIRS
            )
            ax.text(0.98, 0.98, score_text, transform=ax.transAxes, ha="right", va="top",
                    bbox={"facecolor": "white", "edgecolor": "0.7", "alpha": 0.9})
            ax.set(xlabel="Distance from centroid (r/R)", ylabel="Annular mean / granule mean",
                   title=f"{slot}: {fov}, granule {granule_id}\nselected r = {float(selected[score_column]):.3f}")
            ax.legend(frameon=False, loc="lower left")
            fig.tight_layout()
            plot_path = slot_dir / f"{slot}__{fov}__granule-{granule_id:03d}__rdf.png"
            fig.savefig(plot_path, dpi=220); plt.close(fig)
            rows.append({
                "selection_slot": slot, "selection_pair": f"{left}-{right}", "selection_extreme": extreme,
                "selection_score": float(selected[score_column]), "fov": fov, "granule_id": granule_id,
                "equivalent_diameter_px": float(selected.equivalent_diameter_px),
                **{f"pearson_r_{a}_{b}": float(selected[f"pearson_r_{a}_{b}"]) for a, b in CORRELATION_PAIRS},
                "representative_tif": str(destination_crop.relative_to(output_dir)),
                "representative_plot": str(plot_path.relative_to(output_dir)),
            })
    index = pd.DataFrame.from_records(rows)
    index.to_csv(representative_root / "representative_index.csv", index=False)
    return index


def synthetic_validation() -> None:
    shape = (176, 176)
    yy, xx = np.indices(shape)
    images = {ch: np.zeros(shape, np.float32) for ch in CHANNELS}
    # Diffraction-limited object and a shared object <=15 px must be rejected.
    small = np.zeros(shape, dtype=bool)
    small[28:32, 29:32] = True  # 12 px; equivalent diameter 3.91 px.
    threshold_sized = np.hypot(yy - 65, xx - 65) <= 7
    # Shared larger granule: must be retained by the strict >15 px rule.
    large = np.hypot(yy - 105, xx - 105) <= 10
    # Large single-channel object: must be rejected for support.
    one_channel = np.hypot(yy - 130, xx - 35) <= 10
    for ch in ("488", "561"):
        images[ch][small | threshold_sized | large] = 100
    images["647"][large | one_channel] = 100
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["segmentation"].update({"background_sigma_px": 8, "smooth_sigma_px": 0.5,
                                 "watershed_min_distance_px": 4, "minimum_channel_support_fraction": 0.05})
    _, _, props, _ = segment_granules(images, cfg)
    near_small = props.iloc[((props.centroid_y_px - 30) ** 2 + (props.centroid_x_px - 30) ** 2).argmin()]
    near_threshold = props.iloc[((props.centroid_y_px - 65) ** 2 + (props.centroid_x_px - 65) ** 2).argmin()]
    near_large = props.iloc[((props.centroid_y_px - 105) ** 2 + (props.centroid_x_px - 105) ** 2).argmin()]
    near_one = props.iloc[((props.centroid_y_px - 130) ** 2 + (props.centroid_x_px - 35) ** 2).argmin()]
    assert not bool(near_small.keep_granule), "diffraction-limited synthetic molecule was retained"
    assert not bool(near_threshold.keep_granule), "synthetic object <=15 px was retained"
    assert bool(near_large.keep_granule), "shared synthetic granule was rejected"
    assert not bool(near_one.keep_granule), "single-channel synthetic object was retained"


def run_pipeline(config: str | Path | dict) -> PipelineResult:
    cfg = load_config(config)
    output_dir = Path(cfg["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "run_config.generated.yaml").open("w") as handle:
        yaml.safe_dump(cfg, handle, sort_keys=False)
    pairs = pair_four_channel_files(cfg["input_dir"])
    manifest = pd.DataFrame([{ "fov": pair.fov, **{f"path_{ch}": str(pair.paths[ch]) for ch in CHANNELS}} for pair in pairs])
    manifest.to_csv(output_dir / "paired_inputs.csv", index=False)
    all_properties, all_rdf, status_rows = [], [], []
    for index, pair in enumerate(pairs, 1):
        print(f"[{index}/{len(pairs)}] {pair.fov}", flush=True)
        try:
            images = read_fov(pair)
            retained, candidates, properties, diagnostics = segment_granules(images, cfg)
            properties.insert(0, "fov", pair.fov)
            properties["equivalent_radius_nm"] = properties.equivalent_radius_px * float(cfg["pixel_size_nm"])
            properties["equivalent_diameter_nm"] = properties.equivalent_diameter_px * float(cfg["pixel_size_nm"])
            (output_dir / "granule_masks").mkdir(parents=True, exist_ok=True)
            tifffile.imwrite(output_dir / "granule_masks" / f"{pair.fov}__granule_labels.tif", retained,
                             metadata={"axes": "YX"})
            (output_dir / "aggregate_images").mkdir(parents=True, exist_ok=True)
            tifffile.imwrite(output_dir / "aggregate_images" / f"{pair.fov}__normalized_RNA_aggregate.tif",
                             diagnostics["aggregate"].astype(np.float32), metadata={"axes": "YX"})
            crops = export_granule_crops(pair, images, properties, output_dir, int(cfg["crop_padding_px"]))
            fov_rdf = calculate_granule_rdf(pair.fov, images, retained, properties, cfg)
            if cfg["qc"].get("make_overlays", True):
                save_qc_overlay(pair, diagnostics["aggregate"], retained, candidates, properties, output_dir,
                                float(cfg["qc"]["contrast_upper_percentile"]))
            all_properties.append(properties); all_rdf.append(fov_rdf)
            status_rows.append({"fov": pair.fov, "status": "success", "input_files": 4,
                                "candidates": len(properties), "retained_granules": int(properties.keep_granule.sum()),
                                "rejected_candidates": int((~properties.keep_granule).sum()), "crops": len(crops),
                                "rdf_rows": len(fov_rdf), "error": ""})
        except Exception as exc:
            status_rows.append({"fov": pair.fov, "status": "failed", "input_files": 4, "candidates": 0,
                                "retained_granules": 0, "rejected_candidates": 0, "crops": 0, "rdf_rows": 0,
                                "error": f"{type(exc).__name__}: {exc}"})
            pd.DataFrame(status_rows).to_csv(output_dir / "fov_run_status.csv", index=False)
            raise
    properties = pd.concat(all_properties, ignore_index=True) if all_properties else pd.DataFrame()
    rdf = pd.concat(all_rdf, ignore_index=True) if all_rdf else pd.DataFrame()
    aggregate = aggregate_rdf(rdf)
    rdf_wide = rdf_to_wide(rdf)
    correlations = calculate_rdf_correlations(rdf_wide, properties)
    properties.to_csv(output_dir / "granule_properties.csv", index=False)
    rdf.to_csv(output_dir / "granule_rdf_results.csv", index=False)
    rdf_wide.to_csv(output_dir / "granule_rdf_profiles_wide.csv", index=False)
    correlations.to_csv(output_dir / "granule_rdf_correlations.csv", index=False)
    aggregate.to_csv(output_dir / "aggregated_rdf_summary.csv", index=False)
    status = pd.DataFrame(status_rows); status.to_csv(output_dir / "fov_run_status.csv", index=False)
    save_summary_plots(properties, rdf, aggregate, correlations, output_dir, cfg)
    representative_index = export_correlation_representatives(correlations, rdf_wide, output_dir)
    reasons = properties.loc[~properties.keep_granule, "rejection_reason"].str.split(";").explode().value_counts().to_dict()
    retained_count = int(properties.keep_granule.sum())
    correlation_columns = [f"pearson_r_{left}_{right}" for left, right in CORRELATION_PAIRS]
    nonfinite_correlations = int((~np.isfinite(correlations[correlation_columns])).sum().sum())
    crop_count = len(list((output_dir / "granule_crops").glob("*.tif")))
    summary = {
        "status": "success" if len(pairs) == int(cfg["expected_fovs"]) and (status.status == "success").all() else "incomplete",
        "fovs_processed": len(pairs), "input_files": len(pairs) * 4, "candidates": len(properties),
        "retained_granules": retained_count, "rejected_candidates": int((~properties.keep_granule).sum()),
        "rejection_reason_counts": {str(k): int(v) for k, v in reasons.items()},
        "crop_count": crop_count, "rdf_rows": len(rdf), "rdf_wide_rows": len(rdf_wide),
        "rdf_nonfinite_values": int((~np.isfinite(rdf.rdf_normalized)).sum()) if not rdf.empty else 0,
        "correlation_rows": len(correlations), "correlation_nonfinite_values": nonfinite_correlations,
        "representative_slots": len(representative_index),
        "rdf_channels": sorted(rdf.channel.unique().tolist()) if not rdf.empty else [],
        "diameter_px": properties.equivalent_diameter_px.describe().to_dict() if not properties.empty else {},
    }
    with (output_dir / "full_dataset_run_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
    pd.DataFrame([summary | {"rejection_reason_counts": json.dumps(summary["rejection_reason_counts"]),
                             "rdf_channels": ",".join(summary["rdf_channels"]),
                             "diameter_px": json.dumps(summary["diameter_px"]) }]).to_csv(
        output_dir / "full_dataset_run_summary.csv", index=False)
    expected_bins = len(rdf_bins(cfg))
    if (summary["status"] != "success" or summary["rdf_nonfinite_values"] or nonfinite_correlations
            or crop_count != retained_count or len(rdf_wide) != retained_count * expected_bins
            or len(correlations) != retained_count or len(representative_index) != 6):
        raise RuntimeError(f"Production acceptance checks failed: {summary}")
    return PipelineResult(pairs, properties, rdf, aggregate, rdf_wide, correlations, output_dir)


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("check", "run"))
    parser.add_argument("--config", required=False)
    args = parser.parse_args(argv)
    if args.command == "check":
        synthetic_validation(); print("Synthetic validation passed")
    else:
        if not args.config:
            parser.error("run requires --config")
        result = run_pipeline(args.config)
        print(f"Completed {len(result.pairs)} FOVs; retained {int(result.properties.keep_granule.sum())} granules")


if __name__ == "__main__":
    main()
