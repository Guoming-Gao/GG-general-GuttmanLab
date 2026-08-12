"""Generate presentation-ready grouped-cluster PBSA figures.

The plotted idealization preserves quickPBSA's state sequence and transition
locations.  Only the display height of each contiguous plateau is refit to the
median of the exact corrected trace that was supplied to quickPBSA.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

from cluster_counting import _native_result_path, frame_columns, trace_file
from cluster_reanalysis import load_reanalysis_config, roots
from pbsa_shared import atomic_csv


STREAM_FILES = {
    "grouped_cluster": "grouped_cluster_counts.csv",
    "seed_level_cluster": "seed_level_counts.csv",
}
CONDITION_SPECS = (
    ("SHA noDox", "endogenous SPEN", "#0072B2"),
    ("SHA Dox", "endogenous SPEN +Dox", "#D55E00"),
    ("dSPEN FL", "over-expressed SPEN", "#009E73"),
)
CONDITION_ORDER = tuple(item[0] for item in CONDITION_SPECS)
DISPLAY_LABELS = {item[0]: item[1] for item in CONDITION_SPECS}
COLORS = {item[0]: item[2] for item in CONDITION_SPECS}
FONT_SIZE = 11.0


def _bool(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return values.astype(str).str.lower().isin({"true", "1", "yes"})


def configure_style() -> None:
    plt.rcParams.update({
        "font.size": FONT_SIZE,
        "axes.titlesize": FONT_SIZE,
        "axes.labelsize": FONT_SIZE,
        "xtick.labelsize": FONT_SIZE,
        "ytick.labelsize": FONT_SIZE,
        "legend.fontsize": FONT_SIZE,
        "figure.titlesize": FONT_SIZE,
        "figure.labelsize": FONT_SIZE,
        "svg.fonttype": "none",
    })


def load_presentation_counts(cfg: dict[str, Any], stream: str = "grouped_cluster") -> pd.DataFrame:
    if stream not in STREAM_FILES:
        raise ValueError(f"Unsupported presentation stream: {stream}")
    table = pd.read_csv(roots(cfg)["tables"] / STREAM_FILES[stream], low_memory=False)
    selected = table[
        table.analysis_set.eq("primary_standardized")
        & _bool(table.accepted_count)
        & table.condition.isin(CONDITION_ORDER)
    ].copy()
    selected["photobleaching_step_count"] = pd.to_numeric(selected.photobleaching_step_count, errors="coerce")
    selected = selected[selected.photobleaching_step_count.notna()].copy()
    selected["display_condition"] = selected.condition.map(DISPLAY_LABELS)
    if set(selected.condition.unique()) != set(CONDITION_ORDER):
        raise RuntimeError("One or more presentation conditions have no accepted grouped clusters")
    if not selected.analysis_set.eq("primary_standardized").all():
        raise AssertionError("Acquisition-test data entered the presentation table")
    return selected


def histogram_percentages(table: pd.DataFrame) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    maximum = int(table.photobleaching_step_count.max())
    bins = np.arange(0.5, maximum + 1.5, 1.0)
    percentages: dict[str, np.ndarray] = {}
    for condition in CONDITION_ORDER:
        values = table.loc[table.condition.eq(condition), "photobleaching_step_count"].to_numpy(float)
        weights = np.full(len(values), 100.0 / len(values))
        percentages[condition] = np.histogram(values, bins=bins, weights=weights)[0]
    return bins, percentages


def empirical_cdf(values: Iterable[float]) -> tuple[np.ndarray, np.ndarray]:
    ordered = np.sort(np.asarray(list(values), dtype=float))
    if not len(ordered):
        raise ValueError("Cannot calculate an ECDF from no values")
    return np.r_[0.5, ordered], np.r_[0.0, np.arange(1, len(ordered) + 1) / len(ordered)]


def make_histogram_figure(table: pd.DataFrame) -> plt.Figure:
    configure_style()
    bins, percentages = histogram_percentages(table)
    figure, axis = plt.subplots(figsize=(7.5, 5.25))
    for condition in CONDITION_ORDER:
        count = int(table.condition.eq(condition).sum())
        axis.stairs(
            percentages[condition], bins, color=COLORS[condition], linewidth=2.0,
            label=f"{DISPLAY_LABELS[condition]} (n={count:,})",
        )
    axis.set_xlim(bins[0], bins[-1])
    axis.set_ylim(bottom=0)
    axis.set_xlabel("Photobleaching step count")
    axis.set_ylabel("Accepted clusters per bin (%)")
    axis.set_title("Photobleaching step-count distribution")
    axis.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=12))
    axis.yaxis.set_major_locator(MaxNLocator(nbins=6))
    axis.grid(axis="y", color="0.85", linewidth=0.8)
    axis.legend(frameon=False)
    figure.tight_layout()
    return figure


def make_cdf_figure(table: pd.DataFrame) -> plt.Figure:
    configure_style()
    maximum = int(table.photobleaching_step_count.max())
    figure, axis = plt.subplots(figsize=(7.5, 5.25))
    for condition in CONDITION_ORDER:
        values = table.loc[table.condition.eq(condition), "photobleaching_step_count"].to_numpy(float)
        x_values, y_values = empirical_cdf(values)
        axis.step(
            x_values, y_values, where="post", color=COLORS[condition], linewidth=2.0,
            label=f"{DISPLAY_LABELS[condition]} (n={len(values):,})",
        )
    axis.set_xlim(0.5, maximum + 0.5)
    axis.set_ylim(0, 1.01)
    axis.set_xlabel("Photobleaching step count")
    axis.set_ylabel("Cumulative fraction")
    axis.set_title("Cumulative step-count distribution")
    axis.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=12))
    axis.set_yticks(np.linspace(0, 1, 6))
    axis.grid(color="0.85", linewidth=0.8)
    axis.legend(frameon=False, loc="lower right")
    figure.tight_layout()
    return figure


def height_refit_idealization(corrected: np.ndarray, states: np.ndarray) -> np.ndarray:
    """Refit each contiguous quickPBSA state plateau to the trace median."""
    corrected = np.asarray(corrected, dtype=float)
    states = np.asarray(states, dtype=float)
    if corrected.ndim != 1 or states.ndim != 1 or len(corrected) != len(states):
        raise ValueError("Corrected trace and quickPBSA states must be aligned one-dimensional arrays")
    if not len(corrected) or not np.isfinite(states).all():
        raise ValueError("Trace idealization requires finite quickPBSA states")
    changes = np.flatnonzero(np.r_[True, states[1:] != states[:-1], True])
    idealized = np.full(len(corrected), np.nan, dtype=float)
    for start, end in zip(changes[:-1], changes[1:]):
        segment = corrected[start:end]
        finite = segment[np.isfinite(segment)]
        if not len(finite):
            raise ValueError("A quickPBSA plateau contains no finite corrected values")
        idealized[start:end] = float(np.median(finite))
    return idealized


def normalized_fit_residual(corrected: np.ndarray, idealized: np.ndarray) -> float:
    corrected = np.asarray(corrected, dtype=float)
    idealized = np.asarray(idealized, dtype=float)
    finite = np.isfinite(corrected) & np.isfinite(idealized)
    if not finite.any():
        return float("inf")
    scale = float(np.percentile(corrected[finite], 95) - np.percentile(corrected[finite], 5))
    if not np.isfinite(scale) or scale <= 0:
        scale = max(float(np.std(corrected[finite])), 1.0)
    return float(np.sqrt(np.mean(np.square(corrected[finite] - idealized[finite]))) / scale)


@dataclass(frozen=True)
class TraceFit:
    corrected: np.ndarray
    states: np.ndarray
    idealized: np.ndarray
    residual: float


def _trace_fit_from_loaded(traces: pd.DataFrame, native: pd.DataFrame, row: pd.Series) -> TraceFit:
    frames = frame_columns(traces)
    trace_match = traces[traces.analysis_unit_id.astype(str).eq(str(row.analysis_unit_id))]
    if len(trace_match) != 1:
        raise RuntimeError(f"Expected one corrected trace for {row.analysis_unit_id}, found {len(trace_match)}")
    corrected = trace_match.iloc[0][frames].to_numpy(float)

    fit_match = native[
        native.analysis_unit_id.astype(str).eq(str(row.analysis_unit_id))
        & native.type.astype(str).eq("fluors_full")
    ]
    if len(fit_match) != 1:
        raise RuntimeError(f"Expected one quickPBSA state trace for {row.analysis_unit_id}, found {len(fit_match)}")
    missing = [frame for frame in frames if frame not in fit_match.columns]
    if missing:
        raise RuntimeError(f"Native state trace is missing {len(missing)} corrected-trace frames")
    states = fit_match.iloc[0][frames].to_numpy(float)
    idealized = height_refit_idealization(corrected, states)
    return TraceFit(corrected, states, idealized, normalized_fit_residual(corrected, idealized))


def load_trace_fit(cfg: dict[str, Any], row: pd.Series, stream: str) -> TraceFit:
    representation = str(row.trace_representation)
    traces = pd.read_csv(trace_file(cfg, row.dataset, stream, row.fov, representation), low_memory=False)
    native = pd.read_csv(_native_result_path(cfg, row.dataset, stream, row.fov), skiprows=1, low_memory=False)
    return _trace_fit_from_loaded(traces, native, row)


def score_candidates(cfg: dict[str, Any], candidates: pd.DataFrame, stream: str) -> pd.DataFrame:
    scored_groups = []
    keys = ["dataset", "fov", "trace_representation"]
    for (dataset, fov, representation), group in candidates.groupby(keys, sort=False):
        traces = pd.read_csv(trace_file(cfg, dataset, stream, fov, representation), low_memory=False)
        native = pd.read_csv(_native_result_path(cfg, dataset, stream, fov), skiprows=1, low_memory=False)
        group = group.copy()
        group["fit_residual"] = [
            _trace_fit_from_loaded(traces, native, pd.Series(row._asdict())).residual
            for row in group.itertuples(index=False)
        ]
        scored_groups.append(group)
    return pd.concat(scored_groups, ignore_index=True)


def select_monomer_candidates(scored: pd.DataFrame, per_condition: int = 4) -> pd.DataFrame:
    selected = []
    for condition in CONDITION_ORDER:
        group = scored[scored.condition.eq(condition)].copy()
        group["amplitude_error"] = pd.to_numeric(group.single_step_amplitude_relative_error, errors="coerce").fillna(np.inf)
        group = group.sort_values(
            ["fit_residual", "amplitude_error", "dataset", "fov", "analysis_unit_id"],
            kind="stable",
        ).head(per_condition)
        if len(group) != per_condition:
            raise RuntimeError(f"Only {len(group)} monomer candidates are available for {condition}")
        group["selection_rank"] = np.arange(1, per_condition + 1)
        selected.append(group)
    return pd.concat(selected, ignore_index=True)


def highest_step_candidate_pool(table: pd.DataFrame) -> pd.DataFrame:
    records = []
    for (_, _), group in table.groupby(["condition", "fov"], sort=False):
        maximum = group.photobleaching_step_count.max()
        records.append(group[group.photobleaching_step_count.eq(maximum)])
    return pd.concat(records, ignore_index=True)


def select_highest_step_candidates(scored_pool: pd.DataFrame, per_condition: int = 4) -> pd.DataFrame:
    selected = []
    for condition in CONDITION_ORDER:
        group = scored_pool[scored_pool.condition.eq(condition)].copy()
        per_fov = (
            group.sort_values(["photobleaching_step_count", "fit_residual", "analysis_unit_id"], ascending=[False, True, True], kind="stable")
            .groupby("fov", sort=False, as_index=False).head(1)
            .sort_values(["photobleaching_step_count", "fit_residual", "fov", "analysis_unit_id"], ascending=[False, True, True, True], kind="stable")
            .head(per_condition)
        )
        if len(per_fov) < per_condition:
            remaining = group[~group.analysis_unit_id.isin(per_fov.analysis_unit_id)].sort_values(
                ["photobleaching_step_count", "fit_residual", "analysis_unit_id"], ascending=[False, True, True], kind="stable"
            ).head(per_condition - len(per_fov))
            per_fov = pd.concat([per_fov, remaining], ignore_index=True)
        if len(per_fov) != per_condition:
            raise RuntimeError(f"Only {len(per_fov)} highest-step candidates are available for {condition}")
        per_fov = per_fov.sort_values(
            ["photobleaching_step_count", "fit_residual", "fov", "analysis_unit_id"],
            ascending=[False, True, True, True], kind="stable",
        )
        per_fov["selection_rank"] = np.arange(1, per_condition + 1)
        selected.append(per_fov)
    return pd.concat(selected, ignore_index=True)


def build_panel_records(cfg: dict[str, Any], selected: pd.DataFrame, stream: str) -> list[dict[str, Any]]:
    records = []
    for row in selected.itertuples(index=False):
        series = pd.Series(row._asdict())
        fit = load_trace_fit(cfg, series, stream)
        record = series.to_dict()
        record.update({"corrected": fit.corrected, "idealized": fit.idealized, "fit_residual": fit.residual})
        records.append(record)
    return records


def make_trace_panel(records: list[dict[str, Any]], title: str) -> plt.Figure:
    configure_style()
    figure, axes = plt.subplots(4, 3, figsize=(10.5, 9.5), squeeze=False)
    lookup = {(record["condition"], int(record["selection_rank"])): record for record in records}
    for column, condition in enumerate(CONDITION_ORDER):
        for row_index in range(4):
            axis = axes[row_index, column]
            record = lookup[(condition, row_index + 1)]
            corrected = np.asarray(record["corrected"], dtype=float)
            idealized = np.asarray(record["idealized"], dtype=float)
            frames = np.arange(len(corrected))
            axis.plot(frames, corrected, color=COLORS[condition], linewidth=0.7, alpha=0.9)
            axis.plot(frames, idealized, color="black", linewidth=1.25)
            steps = int(record["photobleaching_step_count"])
            axis.text(0.02, 0.96, f"{steps} step" if steps == 1 else f"{steps} steps", transform=axis.transAxes, ha="left", va="top")
            axis.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
            axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
    figure.suptitle(title, y=0.99)
    for column, condition in enumerate(CONDITION_ORDER):
        figure.text((column + 0.5) / 3, 0.945, DISPLAY_LABELS[condition], ha="center", va="center")
    figure.supxlabel("Frame", y=0.055)
    figure.supylabel("Corrected intensity (a.u.)", x=0.025)
    legend = [
        Line2D([0], [0], color="0.35", linewidth=1.2, label="Corrected trace (condition color)"),
        Line2D([0], [0], color="black", linewidth=1.4, label="quickPBSA states, height-refit"),
    ]
    figure.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, 0.002), ncol=2, frameon=False)
    figure.subplots_adjust(left=0.09, right=0.99, bottom=0.13, top=0.90, wspace=0.30, hspace=0.42)
    return figure


def save_figure(figure: plt.Figure, root: Path, stem: str) -> list[Path]:
    paths = [root / f"{stem}.png", root / f"{stem}.svg"]
    figure.savefig(paths[0], dpi=300, bbox_inches="tight", facecolor="white")
    figure.savefig(paths[1], bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return paths


def generate_presentation_figures(cfg: dict[str, Any], stream: str = "grouped_cluster") -> dict[str, Path]:
    table = load_presentation_counts(cfg, stream)
    root = roots(cfg)["root"]
    outputs: dict[str, Path] = {}

    for key, figure, stem in (
        ("histogram", make_histogram_figure(table), "presentation_step_count_histogram"),
        ("cdf", make_cdf_figure(table), "presentation_step_count_cdf"),
    ):
        png, svg = save_figure(figure, root, stem)
        outputs[f"{key}_png"] = png; outputs[f"{key}_svg"] = svg

    monomer_pool = table[_bool(table.monomer_candidate)].copy()
    top_pool = highest_step_candidate_pool(table)
    scoring_pool = pd.concat([monomer_pool, top_pool], ignore_index=True).drop_duplicates(
        ["dataset", "fov", "analysis_unit_id"]
    )
    scored = score_candidates(cfg, scoring_pool, stream)
    score_columns = ["dataset", "fov", "analysis_unit_id", "fit_residual"]
    monomer_scored = monomer_pool.merge(scored[score_columns], on=["dataset", "fov", "analysis_unit_id"], how="left", validate="one_to_one")
    top_scored = top_pool.merge(scored[score_columns], on=["dataset", "fov", "analysis_unit_id"], how="left", validate="one_to_one")
    monomers = select_monomer_candidates(monomer_scored)
    highest = select_highest_step_candidates(top_scored)

    for key, selected, title, stem in (
        ("monomer", monomers, "Monomer-candidate traces", "presentation_monomer_candidate_traces"),
        ("highest", highest, "Highest-step representative traces", "presentation_highest_step_traces"),
    ):
        records = build_panel_records(cfg, selected, stream)
        png, svg = save_figure(make_trace_panel(records, title), root, stem)
        outputs[f"{key}_png"] = png; outputs[f"{key}_svg"] = svg

    manifest_fields = [
        "condition", "display_condition", "dataset", "fov", "analysis_unit_id",
        "photobleaching_step_count", "trace_representation", "fit_residual", "selection_rank",
    ]
    manifest = pd.concat([
        monomers.assign(figure="monomer_candidate_traces"),
        highest.assign(figure="highest_step_traces"),
    ], ignore_index=True)
    manifest["display_condition"] = manifest.condition.map(DISPLAY_LABELS)
    manifest["analysis_stream"] = stream
    manifest_path = roots(cfg)["tables"] / "presentation_trace_selection.csv"
    atomic_csv(manifest[["figure", "analysis_stream", *manifest_fields]], manifest_path)
    outputs["selection_manifest"] = manifest_path
    return outputs


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config_reanalysis.yaml")
    parser.add_argument("--stream", choices=tuple(STREAM_FILES), default="grouped_cluster")
    args = parser.parse_args(argv)
    cfg = load_reanalysis_config(args.config)
    outputs = generate_presentation_figures(cfg, args.stream)
    for name, path in outputs.items():
        print(f"{name}: {path}", flush=True)


if __name__ == "__main__":
    main()
