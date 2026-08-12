"""Consolidated figures, ledgers, and PDF audit report for the revised analysis."""

from __future__ import annotations

import textwrap
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from cluster_counting import _native_result_path, frame_columns, trace_file
from cluster_reanalysis import STREAMS, experiment_root, roots, source_manifest
from pbsa_shared import atomic_csv, atomic_text, condition_metadata, matplotlib_setup


CONDITION_ORDER = ["SHA noDox", "SHA Dox", "dSPEN FL", "dSPEN dRRM"]
STREAM_LABELS = {
    "seed_level_cluster": "Seed-level clusters",
    "grouped_cluster": "Grouped clusters",
    "extended_cluster_candidate": "Extended cluster candidates",
}


def _bool(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values): return values.fillna(False)
    return values.astype(str).str.lower().isin({"true", "1", "yes"})


def load_counts(cfg: dict[str, Any]) -> dict[str, pd.DataFrame]:
    names = {"seed_level_cluster": "seed_level_counts.csv", "grouped_cluster": "grouped_cluster_counts.csv", "extended_cluster_candidate": "extended_cluster_counts.csv"}
    return {stream: pd.read_csv(roots(cfg)["tables"] / filename) for stream, filename in names.items()}


def load_trace_qc(cfg: dict[str, Any], stream: str) -> pd.DataFrame:
    tables = []
    for dataset_cfg in cfg["datasets"]:
        dataset = dataset_cfg["name"]
        for row in source_manifest(cfg, dataset).to_dict("records"):
            path = experiment_root(cfg, dataset) / "puncta_traces" / stream / row["fov"] / f"{row['fov']}__trace_QC.csv"
            if path.exists(): tables.append(pd.read_csv(path))
    return pd.concat(tables, ignore_index=True, sort=False) if tables else pd.DataFrame()


def build_inventory_tables(cfg: dict[str, Any], ledger: pd.DataFrame, counts: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    paths = roots(cfg); manifests = []; drift = []
    for dataset_cfg in cfg["datasets"]:
        dataset = dataset_cfg["name"]; manifest = source_manifest(cfg, dataset).copy(); manifest["dataset"] = dataset
        normalized = manifest.filename.map(condition_metadata).apply(pd.Series)
        for column in normalized.columns: manifest[column] = normalized[column].to_numpy()
        manifests.append(manifest)
        drift_path = Path(cfg["source_results_root"]) / dataset / "02_drift_correction" / "drift_correction_QC_summary.csv"
        if not drift_path.exists(): drift_path = experiment_root(cfg, dataset) / "drift_correction" / "drift_correction_QC_summary.csv"
        table = pd.read_csv(drift_path); table["dataset"] = dataset
        table["corrected_tiff_relative_path"] = table.fov.map(lambda fov: f"04_experiments/{dataset}/drift_correction/corrected_tiffs/{fov}__drift_corrected.tif")
        drift.append(table.drop(columns=[column for column in ["corrected_tiff", "qc_plot", "qc_overlay"] if column in table]))
    inventory = pd.concat(manifests, ignore_index=True, sort=False)
    drift_qc = pd.concat(drift, ignore_index=True, sort=False)
    atomic_csv(inventory, paths["tables"] / "input_and_acquisition_inventory.csv")
    atomic_csv(drift_qc, paths["tables"] / "drift_correction_QC_summary.csv")
    trace_records = []
    for stream in STREAMS:
        qc = load_trace_qc(cfg, stream)
        for keys, group in qc.groupby(["dataset", "fov"], dropna=False):
            trace_records.append({"analysis_stream": stream, "dataset": keys[0], "fov": keys[1], "analysis_units": len(group),
                                  "trace_ready": int((group.trace_status == "trace_ready").sum()),
                                  "trace_rejected": int((group.trace_status != "trace_ready").sum()),
                                  "insufficient_signal_pixels": int((group.trace_rejection_reason == "insufficient_signal_pixels").sum()),
                                  "insufficient_background_pixels": int((group.trace_rejection_reason == "insufficient_background_pixels").sum()),
                                  "nonfinite_trace": int((group.trace_rejection_reason == "nonfinite_trace").sum())})
    trace_attrition = pd.DataFrame.from_records(trace_records)
    atomic_csv(trace_attrition, paths["tables"] / "trace_background_QC_and_attrition.csv")
    mapping = pd.read_csv(paths["tables"] / "seed_to_analysis_unit_mapping.csv")
    growth = (ledger.groupby(["dataset", "fov", "legacy_growth_status"], dropna=False).size().reset_index(name="seeds"))
    atomic_csv(growth, paths["tables"] / "growth_status_routing.csv")
    spot_records = []
    for threshold in [cfg["spotiflow"]["probability_threshold"], *cfg["spotiflow"]["sensitivity_thresholds"]]:
        for keys, group in ledger.groupby(["dataset", "fov", "condition", "analysis_set"], dropna=False):
            spot_records.append({"dataset": keys[0], "fov": keys[1], "condition": keys[2], "analysis_set": keys[3],
                                 "spotiflow_threshold": threshold, "retained_seeds": int((pd.to_numeric(group.spotiflow_probability, errors="coerce") >= threshold).sum())})
    spot_sensitivity = pd.DataFrame.from_records(spot_records)
    atomic_csv(spot_sensitivity, paths["technical"] / "sensitivity_analysis" / "spotiflow_detection_threshold_sensitivity.csv")
    per_fov = []
    flags = []
    for stream, table in counts.items():
        for keys, group in table.groupby(["dataset", "fov", "condition", "analysis_set"], dropna=False):
            accepted = group[_bool(group.accepted_count)]
            per_fov.append({"analysis_stream": stream, "dataset": keys[0], "fov": keys[1], "condition": keys[2], "analysis_set": keys[3],
                            "analysis_units": len(group), "quickpbsa_accepted": len(accepted), "monomer_candidates": int(_bool(group.monomer_candidate).sum()),
                            "median_step_count": float(accepted.photobleaching_step_count.median()) if len(accepted) else np.nan,
                            "outside_validated_range": int((accepted.validated_range_status == "outside_validated_range").sum())})
        flag_table = group_flag_distribution(table, stream); flags.append(flag_table)
    per_fov_table = pd.DataFrame.from_records(per_fov); flag_table = pd.concat(flags, ignore_index=True, sort=False)
    atomic_csv(per_fov_table, paths["tables"] / "per_FOV_count_summary.csv")
    atomic_csv(flag_table, paths["tables"] / "quickpbsa_flag_distribution.csv")
    crowding = []
    for stream in ("seed_level_cluster", "grouped_cluster"):
        table = counts[stream]
        for keys, group in table.groupby(["condition", "analysis_set", "crowded_seed"], dropna=False):
            accepted = group[_bool(group.accepted_count)]
            crowding.append({"analysis_stream": stream, "condition": keys[0], "analysis_set": keys[1], "crowded_seed": keys[2],
                             "accepted": len(accepted), "monomer_candidates": int(_bool(accepted.monomer_candidate).sum()),
                             "monomer_fraction": float(_bool(accepted.monomer_candidate).mean()) if len(accepted) else np.nan})
    crowding_table = pd.DataFrame.from_records(crowding); atomic_csv(crowding_table, paths["tables"] / "monomer_candidate_crowding_summary.csv")
    source_root = Path(cfg["source_results_root"])
    old_seed_files = list(source_root.glob("*/03_roi_detection_and_growth/tables/*__spotiflow_seeds.csv"))
    old_roi_files = list(source_root.glob("*/03_roi_detection_and_growth/tables/*__roi_properties.csv"))
    if not old_seed_files:
        legacy = paths["technical"] / "exploratory_segmentation" / "legacy_growth_inputs"
        old_seed_files = list(legacy.glob("*/tables/*__spotiflow_seeds.csv")); old_roi_files = list(legacy.glob("*/tables/*__roi_properties.csv"))
    revised_seed_qc = load_trace_qc(cfg, "seed_level_cluster"); revised_group_qc = load_trace_qc(cfg, "grouped_cluster")
    existing_funnel = paths["tables"] / "original_vs_revised_analysis_funnel.csv"
    if (source_root / "07_summary_reports" / "combined_photobleaching_step_counts_all_ROIs.csv").exists():
        old_seeds = pd.concat([pd.read_csv(path) for path in old_seed_files], ignore_index=True, sort=False)
        old_rois = pd.concat([pd.read_csv(path) for path in old_roi_files], ignore_index=True, sort=False)
        old_counts = pd.read_csv(source_root / "07_summary_reports" / "combined_photobleaching_step_counts_all_ROIs.csv")
        failed = old_seeds.growth_status.isin(["seed_below_local_threshold", "region_too_small"])
        original_records = [
            {"pipeline": "superseded", "transition": "Spotiflow 0.4 seeds", "units": len(old_seeds), "note": "starting denominator"},
            {"pipeline": "superseded", "transition": "failed-growth seeds", "units": int(failed.sum()), "note": "incorrectly treated as ineligible"},
            {"pipeline": "superseded", "transition": "failed-growth seeds nevertheless assigned to grown masks", "units": int((failed & (old_seeds.roi_id > 0)).sum()), "note": "assignment defect"},
            {"pipeline": "superseded", "transition": "grown/merged ROIs", "units": len(old_rois), "note": f"{int((old_rois.area_px > 1000).sum())} exceeded 1000 px"},
            {"pipeline": "superseded", "transition": "traces passed to quickPBSA", "units": len(old_counts), "note": "one ROI lost before counting"},
            {"pipeline": "superseded", "transition": "quickPBSA accepted", "units": int(_bool(old_counts.accepted_count).sum()), "note": "64/409"},
        ]
    else:
        original_records = pd.read_csv(existing_funnel).query("pipeline == 'superseded'").to_dict("records") if existing_funnel.exists() else []
    old_new_funnel = pd.DataFrame(original_records + [
        {"pipeline": "revised", "transition": "Spotiflow 0.4 seeds", "units": len(ledger), "note": "all entered both co-primary streams"},
        {"pipeline": "revised", "transition": "seed-level traces ready", "units": int((revised_seed_qc.trace_status == "trace_ready").sum()), "note": "failed growth retained"},
        {"pipeline": "revised", "transition": "grouped-cluster units", "units": mapping[["dataset", "fov", "grouped_unit_id"]].drop_duplicates().shape[0], "note": "nonduplicated connected components"},
        {"pipeline": "revised", "transition": "grouped-cluster traces ready", "units": int((revised_group_qc.trace_status == "trace_ready").sum()), "note": "neighbor footprints excluded from background"},
    ])
    atomic_csv(old_new_funnel, paths["tables"] / "original_vs_revised_analysis_funnel.csv")
    return {"inventory": inventory, "drift": drift_qc, "trace_attrition": trace_attrition, "growth": growth,
            "spot_sensitivity": spot_sensitivity, "per_fov": per_fov_table, "flags": flag_table, "crowding": crowding_table,
            "old_new_funnel": old_new_funnel}


def group_flag_distribution(table: pd.DataFrame, stream: str) -> pd.DataFrame:
    data = table.copy(); data["pbsa_flag_label"] = data.pbsa_flag_meaning.fillna("calibration_blocked")
    return (data.groupby(["analysis_set", "condition", "pbsa_flag", "pbsa_flag_label"], dropna=False).size()
            .reset_index(name="analysis_units").assign(analysis_stream=stream))


def enrich_seed_ledger(cfg: dict[str, Any], counts: dict[str, pd.DataFrame]) -> pd.DataFrame:
    ledger = pd.read_csv(roots(cfg)["tables"] / "complete_seed_ledger.csv")
    mapping = pd.read_csv(roots(cfg)["tables"] / "seed_to_analysis_unit_mapping.csv")
    # The report command is intentionally resumable.  Remove fields written by a
    # previous report pass before rebuilding them, otherwise pandas suffixes the
    # mapping keys and the stream joins can no longer find ``extended_unit_id``.
    ledger = ledger.drop(columns=[
        column for column in ["extended_unit_id", "extended_segmentation_status", "extended_rejection_reason"]
        if column in ledger
    ])
    ledger = ledger.merge(mapping[["dataset", "fov", "seed_id", "extended_unit_id", "extended_segmentation_status", "extended_rejection_reason"]], on=["dataset", "fov", "seed_id"], how="left")
    for stream, unit_column, prefix in [
        ("seed_level_cluster", "seed_level_unit_id", "seed_level"),
        ("grouped_cluster", "grouped_unit_id", "grouped"),
        ("extended_cluster_candidate", "extended_unit_id", "extended"),
    ]:
        qc = load_trace_qc(cfg, stream)
        if len(qc):
            qc = qc[["dataset", "fov", "analysis_unit_id", "trace_status", "trace_rejection_reason", "signal_pixels", "background_pixels"]].drop_duplicates(["dataset", "fov", "analysis_unit_id"])
            qc = qc.rename(columns={column: f"{prefix}_{column}" for column in ["trace_status", "trace_rejection_reason", "signal_pixels", "background_pixels"]})
            ledger = ledger.drop(columns=[column for column in qc.columns if column.startswith(f"{prefix}_") and column in ledger])
            ledger = ledger.merge(qc, left_on=["dataset", "fov", unit_column], right_on=["dataset", "fov", "analysis_unit_id"], how="left").drop(columns=["analysis_unit_id"])
        count = counts[stream]
        if len(count):
            fields = ["dataset", "fov", "analysis_unit_id", "pbsa_flag", "pbsa_flag_meaning", "photobleaching_step_count", "accepted_count", "monomer_candidate", "validated_range_status", "calibration_status"]
            count = count[[column for column in fields if column in count]].drop_duplicates(["dataset", "fov", "analysis_unit_id"])
            count = count.rename(columns={column: f"{prefix}_{column}" for column in fields[3:] if column in count})
            ledger = ledger.drop(columns=[column for column in count.columns if column.startswith(f"{prefix}_") and column in ledger])
            ledger = ledger.merge(count, left_on=["dataset", "fov", unit_column], right_on=["dataset", "fov", "analysis_unit_id"], how="left").drop(columns=["analysis_unit_id"])
    atomic_csv(ledger, roots(cfg)["tables"] / "complete_seed_ledger.csv")
    return ledger


def count_summary(counts: dict[str, pd.DataFrame]) -> pd.DataFrame:
    records = []
    for stream, table in counts.items():
        if table.empty: continue
        for keys, group in table.groupby(["analysis_set", "condition"], dropna=False):
            accepted = group[_bool(group.accepted_count)]
            monomers = accepted[_bool(accepted.monomer_candidate)] if len(accepted) else accepted
            records.append({"analysis_stream": stream, "analysis_set": keys[0], "condition": keys[1],
                            "analysis_units": len(group), "quickpbsa_accepted": len(accepted),
                            "accepted_fraction": len(accepted) / len(group) if len(group) else 0,
                            "monomer_candidates": len(monomers), "monomer_fraction_of_accepted": len(monomers) / len(accepted) if len(accepted) else np.nan,
                            "median_step_count": float(accepted.photobleaching_step_count.median()) if len(accepted) else np.nan,
                            "outside_validated_range": int((accepted.validated_range_status == "outside_validated_range").sum())})
    return pd.DataFrame.from_records(records)


def make_histogram(cfg: dict[str, Any], table: pd.DataFrame, stream: str) -> Path:
    matplotlib_setup(cfg["output_root"]); import matplotlib.pyplot as plt
    selected = table[_bool(table.is_primary_comparison) & _bool(table.accepted_count)].copy()
    selected["photobleaching_step_count"] = pd.to_numeric(selected.photobleaching_step_count, errors="coerce")
    selected = selected[selected.photobleaching_step_count.notna()]
    maximum = int(max(41, selected.photobleaching_step_count.max() if len(selected) else 41)); bins = np.arange(.5, maximum + 1.5, 1)
    fig, ax = plt.subplots(figsize=(10, 6.5))
    for condition in CONDITION_ORDER:
        values = selected.loc[selected.condition == condition, "photobleaching_step_count"].to_numpy(float)
        if not len(values): continue
        ax.hist(values, bins=bins, weights=np.full(len(values), 100 / len(values)), histtype="step", lw=2, label=f"{condition} (n={len(values):,})")
    ax.axvline(40.5, color="black", ls="--", lw=1.2)
    ax.text(40.5, .98, "validated range ends at 40", transform=ax.get_xaxis_transform(), ha="right", va="top", fontsize=8)
    ax.set_xlim(left=.5); ax.set_xlabel("Photobleaching step count"); ax.set_ylabel("Accepted analysis units per integer bin (%)")
    ax.set_title(STREAM_LABELS[stream]); ax.grid(axis="y", alpha=.2)
    handles, labels = ax.get_legend_handles_labels(); fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(.5, .015), ncol=4, frameon=False)
    fig.subplots_adjust(bottom=.22)
    folder = {"seed_level_cluster": "seed_level_clusters", "grouped_cluster": "grouped_clusters", "extended_cluster_candidate": "extended_clusters"}[stream]
    path = roots(cfg)["figures"] / folder / f"{stream}__condition_step_count_histogram.png"
    fig.savefig(path, dpi=300, bbox_inches="tight"); plt.close(fig); return path


def make_monomer_figure(cfg: dict[str, Any], counts: dict[str, pd.DataFrame]) -> Path:
    matplotlib_setup(cfg["output_root"]); import matplotlib.pyplot as plt
    rows = []
    for stream in ("seed_level_cluster", "grouped_cluster"):
        table = counts[stream]; table = table[_bool(table.is_primary_comparison) & _bool(table.accepted_count)]
        for condition in CONDITION_ORDER:
            group = table[table.condition == condition]
            rows.append({"stream": STREAM_LABELS[stream], "condition": condition, "accepted": len(group),
                         "monomer_percent": 100 * _bool(group.monomer_candidate).mean() if len(group) else 0})
    data = pd.DataFrame.from_records(rows); x = np.arange(len(CONDITION_ORDER)); width = .36
    fig, ax = plt.subplots(figsize=(10, 6))
    for index, stream in enumerate([STREAM_LABELS["seed_level_cluster"], STREAM_LABELS["grouped_cluster"]]):
        group = data[data.stream == stream].set_index("condition").reindex(CONDITION_ORDER)
        labels = [f"{value:.1f}%\n(n={int(n)})" for value, n in zip(group.monomer_percent, group.accepted)]
        bars = ax.bar(x + (index - .5) * width, group.monomer_percent, width, label=stream)
        ax.bar_label(bars, labels=labels, fontsize=8, padding=3)
    ax.set_xticks(x, CONDITION_ORDER); ax.set_ylabel("Monomer candidates among accepted units (%)"); ax.set_title("Single-step population")
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.12), ncol=2, frameon=False); ax.set_ylim(0, max(10, data.monomer_percent.max() * 1.25 if len(data) else 10)); fig.subplots_adjust(bottom=.22)
    path = roots(cfg)["figures"] / "method_QC" / "monomer_candidate_comparison.png"; fig.savefig(path, dpi=300, bbox_inches="tight"); plt.close(fig); return path


def make_funnel_figure(cfg: dict[str, Any], ledger: pd.DataFrame, counts: dict[str, pd.DataFrame]) -> Path:
    matplotlib_setup(cfg["output_root"]); import matplotlib.pyplot as plt
    mapping = pd.read_csv(roots(cfg)["tables"] / "seed_to_analysis_unit_mapping.csv")
    values = {
        "Raw Spotiflow seeds": len(ledger),
        "Seed traces formed": int((ledger.seed_level_trace_status == "trace_ready").sum()),
        "Seed quickPBSA accepted": int(_bool(counts["seed_level_cluster"].accepted_count).sum()),
        "Grouped units": len(mapping[["dataset", "fov", "grouped_unit_id"]].drop_duplicates()),
        "Grouped traces formed": int((load_trace_qc(cfg, "grouped_cluster").trace_status == "trace_ready").sum()),
        "Grouped quickPBSA accepted": int(_bool(counts["grouped_cluster"].accepted_count).sum()),
        "Valid extended units": int(mapping.loc[mapping.extended_segmentation_status == "valid", ["dataset", "fov", "extended_unit_id"]].drop_duplicates().shape[0]),
        "Extended quickPBSA accepted": int(_bool(counts["extended_cluster_candidate"].accepted_count).sum()),
    }
    fig, ax = plt.subplots(figsize=(11, 6)); labels = list(values); y = np.arange(len(labels)); bars = ax.barh(y, list(values.values()), color=["#2563EB"] * 3 + ["#7C3AED"] * 3 + ["#D97706"] * 2)
    ax.set_yticks(y, labels); ax.invert_yaxis(); ax.set_xlabel("Analysis units (log scale)"); ax.set_xscale("symlog", linthresh=1); ax.bar_label(bars, labels=[f"{value:,}" for value in values.values()], padding=4)
    ax.set_title("Revised analysis funnel - denominators shown explicitly"); fig.tight_layout()
    path = roots(cfg)["figures"] / "method_QC" / "revised_analysis_funnel.png"; fig.savefig(path, dpi=300, bbox_inches="tight"); plt.close(fig); return path


def make_growth_routing_figure(cfg: dict[str, Any], ledger: pd.DataFrame) -> Path:
    matplotlib_setup(cfg["output_root"]); import matplotlib.pyplot as plt
    data = ledger.legacy_growth_status.fillna("unrecorded").value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(10, 6)); bars = ax.barh(data.index, data.values, color="#0F766E")
    ax.bar_label(bars, labels=[f"{value:,}" for value in data.values], padding=3); ax.set_xlabel("Spotiflow seeds")
    ax.set_title("Growth is morphology routing, not eligibility filtering"); ax.text(.99, .02, "All statuses retain seed-level and grouped-cluster mappings", transform=ax.transAxes, ha="right", fontsize=9)
    fig.tight_layout(); path = roots(cfg)["figures"] / "method_QC" / "growth_status_routing.png"; fig.savefig(path, dpi=300); plt.close(fig); return path


def make_spotiflow_sensitivity_figure(cfg: dict[str, Any], sensitivity: pd.DataFrame) -> Path:
    matplotlib_setup(cfg["output_root"]); import matplotlib.pyplot as plt
    primary = sensitivity[sensitivity.analysis_set == "primary_standardized"].groupby(["condition", "spotiflow_threshold"]).retained_seeds.sum().reset_index()
    fig, ax = plt.subplots(figsize=(10, 6))
    for condition in CONDITION_ORDER:
        group = primary[primary.condition == condition].sort_values("spotiflow_threshold")
        if len(group): ax.plot(group.spotiflow_threshold, group.retained_seeds, marker="o", lw=2, label=condition)
    ax.set_xlabel("Spotiflow probability threshold"); ax.set_ylabel("Retained seeds across primary FOVs"); ax.set_title("Spotiflow threshold sensitivity")
    ax.set_xticks(sorted(sensitivity.spotiflow_threshold.unique())); ax.grid(alpha=.2); ax.legend(loc="upper center", bbox_to_anchor=(.5, -.13), ncol=4, frameon=False); fig.subplots_adjust(bottom=.22)
    path = roots(cfg)["figures"] / "method_QC" / "spotiflow_threshold_sensitivity.png"; fig.savefig(path, dpi=300); plt.close(fig); return path


def make_mask_comparison(cfg: dict[str, Any]) -> Path:
    matplotlib_setup(cfg["output_root"]); import matplotlib.pyplot as plt
    examples = [
        (cfg["datasets"][0]["name"], "standarized-SHA_noDox-FOV-1"),
        (cfg["datasets"][1]["name"], "standarized-dSPEN_FL-FOV-1"),
    ]
    fig, axes = plt.subplots(len(examples), 2, figsize=(11, 12), squeeze=False)
    for row_index, (dataset, fov) in enumerate(examples):
        old = Path(cfg["source_results_root"]) / dataset / "03_roi_detection_and_growth" / "overlays" / f"{fov}__roi_overlay.png"
        if not old.exists():
            candidates = list((old.parent).glob(f"{fov}*.png")); old = candidates[0] if candidates else old
        if not old.exists():
            candidates = list((roots(cfg)["technical"] / "exploratory_segmentation" / "legacy_growth_inputs" / dataset / "overlays").glob(f"{fov}*.png"))
            old = candidates[0] if candidates else old
        revised = roots(cfg)["figures"] / "method_QC" / f"{dataset}__{fov}__analysis_units.png"
        for column, (image_path, title) in enumerate([(old, "Superseded unrestricted growth/union"), (revised, "Revised seed-preserving watershed")]):
            axes[row_index, column].axis("off"); axes[row_index, column].set_title(f"{fov}\n{title}", fontsize=9)
            if image_path.exists(): axes[row_index, column].imshow(plt.imread(image_path))
    fig.suptitle("Representative mask correction: seeds remain eligible even when extended segmentation fails")
    fig.tight_layout(rect=[0, 0, 1, .97]); path = roots(cfg)["figures"] / "method_QC" / "superseded_vs_revised_masks.png"
    fig.savefig(path, dpi=220); plt.close(fig); return path


def monomer_trace_montage(cfg: dict[str, Any], counts: pd.DataFrame) -> Path | None:
    candidates = counts[_bool(counts.monomer_candidate)].sort_values(["dataset", "fov", "analysis_unit_id"])
    maximum = int(cfg["qc"]["maximum_review_traces"]); half = maximum // 2
    isolated_pool = candidates[~_bool(candidates.crowded_seed)]
    crowded_pool = candidates[_bool(candidates.crowded_seed)]
    isolated = isolated_pool.sample(n=min(half, len(isolated_pool)), random_state=20260722)
    crowded = crowded_pool.sample(n=min(maximum - len(isolated), len(crowded_pool)), random_state=20260723)
    selected = pd.concat([isolated, crowded], ignore_index=True)
    if len(selected) < maximum:
        remaining = candidates[~candidates.analysis_unit_id.isin(selected.analysis_unit_id)]
        selected = pd.concat([selected, remaining.sample(n=min(maximum - len(selected), len(remaining)), random_state=20260724)], ignore_index=True)
    if selected.empty: return None
    review_columns = ["dataset", "fov", "analysis_unit_id", "condition", "analysis_set", "crowded_seed",
                      "fitted_single_step_amplitude", "calibrated_single_step_amplitude", "single_step_amplitude_within_20pct"]
    review = selected[[column for column in review_columns if column in selected]].copy()
    review["visual_review_status"] = "pending"
    review["visual_review_note"] = ""
    review_path = roots(cfg)["tables"] / "monomer_candidate_visual_review.csv"
    if review_path.exists():
        prior = pd.read_csv(review_path).set_index(["dataset", "fov", "analysis_unit_id"])
        keyed = review.set_index(["dataset", "fov", "analysis_unit_id"])
        common = keyed.index.intersection(prior.index)
        for column in ("visual_review_status", "visual_review_note"):
            if column in prior:
                keyed.loc[common, column] = prior.loc[common, column]
        review = keyed.reset_index()
    atomic_csv(review, review_path)
    matplotlib_setup(cfg["output_root"]); import matplotlib.pyplot as plt
    fig, axes = plt.subplots(6, 5, figsize=(14, 12), squeeze=False)
    for ax in axes.ravel(): ax.axis("off")
    for plot_index, row in enumerate(selected.itertuples()):
        representation = getattr(row, "trace_representation", "mean_difference")
        table = pd.read_csv(trace_file(cfg, row.dataset, "seed_level_cluster", row.fov, representation)); frames = frame_columns(table)
        matched = table[table.analysis_unit_id == row.analysis_unit_id]
        if matched.empty: continue
        values = matched.iloc[0][frames].to_numpy(float); ax = axes.ravel()[plot_index]; ax.axis("on"); ax.plot(values, lw=.65, label="corrected")
        native_path = _native_result_path(cfg, row.dataset, "seed_level_cluster", row.fov)
        if native_path.exists():
            native = pd.read_csv(native_path, skiprows=1, low_memory=False)
            fitted = native[(native.analysis_unit_id == row.analysis_unit_id) & (native.type == "fluors_full")]
            if len(fitted):
                occupancy = fitted.iloc[0][frames].to_numpy(float); amplitude = float(fitted.iloc[0]["laststep"]); baseline = float(fitted.iloc[0]["bg"])
                ax.plot(baseline + occupancy * amplitude, lw=1, color="black", alpha=.8, label="quickPBSA fit")
        crowd = "crowded" if bool(row.crowded_seed) else "isolated"
        ax.set_title(f"{row.analysis_unit_id}\n{row.condition}; {crowd}", fontsize=7); ax.set_xlabel("Frame", fontsize=7); ax.tick_params(labelsize=6)
    axes.ravel()[0].legend(fontsize=6, frameon=False)
    fig.suptitle("Random review sample of 30 seed-level monomer-candidate traces and fits (fixed sampling seed)"); fig.tight_layout(rect=[0, 0, 1, .97])
    path = roots(cfg)["figures"] / "seed_level_clusters" / "monomer_candidate_trace_review.png"; fig.savefig(path, dpi=220); plt.close(fig); return path


def representative_trace_components(cfg: dict[str, Any], counts: pd.DataFrame) -> Path | None:
    selected = counts[_bool(counts.monomer_candidate)].sort_values(["dataset", "fov", "analysis_unit_id"]).head(6)
    if selected.empty: return None
    matplotlib_setup(cfg["output_root"]); import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 2, figsize=(12, 10), squeeze=False)
    for ax, row in zip(axes.ravel(), selected.itertuples()):
        traces = {}
        for representation in ("signal_mean", "background_mean", "mean_difference"):
            table = pd.read_csv(trace_file(cfg, row.dataset, "seed_level_cluster", row.fov, representation)); frames = frame_columns(table)
            match = table[table.analysis_unit_id == row.analysis_unit_id]
            if len(match): traces[representation] = match.iloc[0][frames].to_numpy(float)
        if "signal_mean" in traces: ax.plot(traces["signal_mean"], lw=.55, alpha=.8, label="signal mean")
        if "background_mean" in traces: ax.plot(traces["background_mean"], lw=.55, alpha=.8, label="background mean")
        if "mean_difference" in traces: ax.plot(traces["mean_difference"], lw=.8, label="corrected")
        native_path = _native_result_path(cfg, row.dataset, "seed_level_cluster", row.fov)
        if native_path.exists() and "mean_difference" in traces:
            native = pd.read_csv(native_path, skiprows=1, low_memory=False); fitted = native[(native.analysis_unit_id == row.analysis_unit_id) & (native.type == "fluors_full")]
            if len(fitted):
                occupancy = fitted.iloc[0][frames].to_numpy(float); ax.plot(float(fitted.iloc[0]["bg"]) + occupancy * float(fitted.iloc[0]["laststep"]), color="black", lw=1.2, label="quickPBSA fit")
        ax.set_title(f"{row.analysis_unit_id} | {row.condition}", fontsize=8); ax.set_xlabel("Frame"); ax.set_ylabel("Intensity")
    axes.ravel()[0].legend(ncol=2, fontsize=7, frameon=False); fig.suptitle("Representative monomer-candidate trace components"); fig.tight_layout(rect=[0, 0, 1, .96])
    path = roots(cfg)["figures"] / "seed_level_clusters" / "representative_monomer_trace_components.png"; fig.savefig(path, dpi=220); plt.close(fig); return path


def _table_page(pdf, title: str, table: pd.DataFrame, page: int, total: int, fontsize: float = 7) -> None:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(11, 8.5)); ax.axis("off")
    text = table.to_string(index=False, max_rows=45) if len(table) else "No rows"
    ax.text(.03, .95, title + "\n\n" + text, va="top", family="monospace", fontsize=fontsize)
    fig.text(.5, .02, f"Page {page} of {total}", ha="center", fontsize=8, color=".4"); pdf.savefig(fig); plt.close(fig)


def generate_report(cfg: dict[str, Any]) -> dict[str, Path]:
    paths = roots(cfg); counts = load_counts(cfg); ledger = enrich_seed_ledger(cfg, counts); summary = count_summary(counts)
    atomic_csv(summary, paths["tables"] / "analysis_stream_condition_summary.csv")
    audit = build_inventory_tables(cfg, ledger, counts)
    histograms = {stream: make_histogram(cfg, table, stream) for stream, table in counts.items()}
    monomer = make_monomer_figure(cfg, counts); funnel = make_funnel_figure(cfg, ledger, counts); growth_figure = make_growth_routing_figure(cfg, ledger)
    spot_figure = make_spotiflow_sensitivity_figure(cfg, audit["spot_sensitivity"]); mask_comparison = make_mask_comparison(cfg); montage = monomer_trace_montage(cfg, counts["seed_level_cluster"])
    components = representative_trace_components(cfg, counts["seed_level_cluster"])
    thresholds = pd.read_csv(paths["tables"] / "quickpbsa_thresholds.csv")
    production = thresholds[thresholds.selected_representation.astype(bool)][["analysis_stream", "acquisition_profile", "trace_representation", "trace_count", "single_step_events", "estimated_single_step", "selected_threshold", "late_baseline_pass_fraction", "calibration_status"]]
    original_growth = ledger.legacy_growth_status.value_counts().rename_axis("legacy_growth_status").reset_index(name="seeds")
    inventory_summary = (audit["inventory"].groupby(["analysis_set", "condition", "acquisition_profile"], dropna=False).size().reset_index(name="FOVs"))
    drift_summary = (audit["drift"].groupby("status").agg(FOVs=("fov", "size"), maximum_median_residual_px=("median_residual_shift_px", "max"), maximum_p95_residual_px=("p95_residual_shift_px", "max")).reset_index())
    trace_summary = (audit["trace_attrition"].groupby("analysis_stream")[["analysis_units", "trace_ready", "trace_rejected", "insufficient_signal_pixels", "insufficient_background_pixels", "nonfinite_trace"]].sum().reset_index())
    flag_summary = (audit["flags"].groupby(["analysis_stream", "pbsa_flag_label"], dropna=False).analysis_units.sum().reset_index())
    acquisition_appendix = audit["per_fov"][audit["per_fov"].analysis_set == "acquisition_test"]
    monomer_validation = []
    for stream in ("seed_level_cluster", "grouped_cluster"):
        table = counts[stream]; candidates = table[_bool(table.monomer_candidate)]
        for crowded, group in candidates.groupby("crowded_seed", dropna=False):
            fitted = pd.to_numeric(group.get("fitted_single_step_amplitude"), errors="coerce"); calibrated = pd.to_numeric(group.get("calibrated_single_step_amplitude"), errors="coerce")
            fitted_median = float(fitted.median()) if fitted.notna().any() else np.nan; calibrated_median = float(calibrated.median()) if calibrated.notna().any() else np.nan
            median_relative_error = abs(fitted_median - calibrated_median) / calibrated_median if np.isfinite(fitted_median) and np.isfinite(calibrated_median) and calibrated_median else np.nan
            monomer_validation.append({"analysis_stream": stream, "crowded_seed": crowded, "monomer_candidates": len(group),
                                       "with_fitted_amplitude": int(fitted.notna().sum()), "median_fitted_step_amplitude": fitted_median,
                                       "median_calibrated_step_amplitude": calibrated_median, "population_median_relative_error": median_relative_error,
                                       "population_median_within_20pct": bool(median_relative_error <= .2) if np.isfinite(median_relative_error) else False,
                                       "within_20pct_of_calibration": int(_bool(group.get("single_step_amplitude_within_20pct", pd.Series(False, index=group.index))).sum()),
                                       "within_20pct_fraction": float(_bool(group.get("single_step_amplitude_within_20pct", pd.Series(False, index=group.index))).mean()) if len(group) else np.nan})
    monomer_validation = pd.DataFrame.from_records(monomer_validation); atomic_csv(monomer_validation, paths["tables"] / "monomer_candidate_amplitude_validation.csv")
    sweep_path = paths["technical"] / "sensitivity_analysis" / "quickpbsa_threshold_sweep_results.csv"
    sweep = pd.read_csv(sweep_path) if sweep_path.exists() else pd.DataFrame()
    narrative = [
        "Pipeline audit and cluster-aware PBSA reanalysis", "",
        f"All {len(ledger):,} Spotiflow detections are preserved in the canonical seed ledger.",
        "The superseded pipeline compressed 14,712 seeds into 410 merged masks, formed 409 traces, and accepted 64 quickPBSA results. It also assigned 5,520 failed-growth seeds to grown masks, explaining inflated mask seed counts while simultaneously excluding those seeds as independent clusters.",
        f"Failed local growth is treated as morphology only: {int(ledger.legacy_growth_status.isin(['seed_below_local_threshold','region_too_small']).sum()):,} failed-growth seeds still receive seed-level and grouped-cluster analysis mappings.",
        "Seed-level and grouped-cluster outputs are co-primary but are never pooled because overlapping seed measurements are not independent.",
        "Extended masks are marker-controlled watershed candidates; cell-scale, boundary-touching, and background-deficient masks are rejected only from the extended stream.",
        "Background subtraction precedes threshold calibration because quickPBSA thresholds are defined in corrected-trace intensity units.",
        "A monomer candidate is an accepted one-step trace, not proof of biochemical monomeric state or complete labeling.",
        "Production calibration requires at least 30 persistent late-step events, baseline and ROI-area checks, and visual montage review. Extended masks failed independent calibration and remain exploratory/uninterpreted.",
        "Primary biological histograms use percentages of accepted units and keep the seed-level and grouped-cluster counting units separate. Acquisition-test FOVs are reported only in the appendix.",
    ]
    image_pages = [funnel, growth_figure, spot_figure, mask_comparison, histograms["seed_level_cluster"], histograms["grouped_cluster"], histograms["extended_cluster_candidate"], monomer]
    if montage is not None: image_pages.append(montage)
    if components is not None: image_pages.append(components)
    review_path = paths["tables"] / "monomer_candidate_visual_review.csv"
    visual_review = pd.read_csv(review_path) if review_path.exists() else pd.DataFrame()
    if len(visual_review):
        reviewed = visual_review.visual_review_status.ne("pending").sum()
        review_counts = visual_review.visual_review_status.value_counts().to_dict()
        narrative.append(f"Visual review covered {reviewed}/{len(visual_review)} randomly sampled one-step traces; classifications: {review_counts}.")
    profile_labels = {
        "640nm_12pct_100p006ms_1x": "12%@100ms",
        "640nm_14p9929pct_100p006ms_1x": "14.993%@100ms",
        "640nm_20pct_100p006ms_1x": "20%@100ms",
        "640nm_30pct_30p0029ms_1x": "30%@30ms",
    }
    stream_labels = {"seed_level_cluster": "seed", "grouped_cluster": "grouped", "extended_cluster_candidate": "extended"}
    sweep_display = sweep[["analysis_stream", "acquisition_profile", "threshold_multiplier", "threshold", "pilot_traces", "accepted", "one_step_accepted", "median_accepted_count"]].copy() if len(sweep) else sweep
    if len(sweep_display):
        sweep_display["analysis_stream"] = sweep_display.analysis_stream.map(stream_labels)
        sweep_display["acquisition_profile"] = sweep_display.acquisition_profile.map(profile_labels).fillna(sweep_display.acquisition_profile)
        sweep_display = sweep_display.rename(columns={"analysis_stream": "stream", "acquisition_profile": "profile", "threshold_multiplier": "mult", "pilot_traces": "n", "one_step_accepted": "one_step", "median_accepted_count": "median_count"})
    visual_review_display = visual_review.copy()
    if len(visual_review_display):
        visual_review_display["amplitude_ratio"] = (pd.to_numeric(visual_review_display.fitted_single_step_amplitude, errors="coerce") / pd.to_numeric(visual_review_display.calibrated_single_step_amplitude, errors="coerce")).round(2)
        visual_review_display = visual_review_display[["analysis_unit_id", "analysis_set", "crowded_seed", "amplitude_ratio", "single_step_amplitude_within_20pct", "visual_review_status", "visual_review_note"]]
        visual_review_display = visual_review_display.rename(columns={"analysis_unit_id": "unit", "analysis_set": "set", "crowded_seed": "crowded", "single_step_amplitude_within_20pct": "within20", "visual_review_status": "review", "visual_review_note": "note"})
        visual_review_display["set"] = visual_review_display["set"].replace({"primary_standardized": "primary", "acquisition_test": "test"})
    acquisition_display = acquisition_appendix.copy()
    if len(acquisition_display):
        acquisition_display["analysis_stream"] = acquisition_display.analysis_stream.map(stream_labels)
        acquisition_display = acquisition_display[["analysis_stream", "fov", "analysis_units", "quickpbsa_accepted", "monomer_candidates", "median_step_count", "outside_validated_range"]]
        acquisition_display = acquisition_display.rename(columns={"analysis_stream": "stream", "analysis_units": "units", "quickpbsa_accepted": "accepted", "monomer_candidates": "one_step", "median_step_count": "median", "outside_validated_range": ">40"})
    tables = [("Original versus revised analysis funnel", audit["old_new_funnel"], 7),
              ("Input and acquisition inventory", inventory_summary, 7), ("Drift-correction effectiveness", drift_summary, 7),
              ("Legacy growth-status inventory", original_growth, 8), ("Trace-background QC and attrition", trace_summary, 7),
              ("Selected calibration records", production, 6.5), ("quickPBSA threshold sensitivity", sweep_display, 7),
              ("quickPBSA flag distribution", flag_summary, 7), ("Condition summary", summary, 7),
              ("Monomer-candidate amplitude validation", monomer_validation, 7),
              ("Random monomer-candidate visual review", visual_review_display, 5.8),
              ("Acquisition-test appendix", acquisition_display, 7)]
    total = 1 + len(tables) + len(image_pages); report = paths["report"] / "pipeline_audit_and_reanalysis_report.pdf"; report.parent.mkdir(parents=True, exist_ok=True)
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(report) as pdf:
        fig, ax = plt.subplots(figsize=(11, 8.5)); ax.axis("off")
        wrapped = "\n".join(textwrap.fill(line, 105) if line else "" for line in narrative)
        ax.text(.05, .94, wrapped, va="top", fontsize=11, linespacing=1.5); fig.text(.5, .02, f"Page 1 of {total}", ha="center", fontsize=8, color=".4"); pdf.savefig(fig); plt.close(fig)
        page = 2
        for title, table, fontsize in tables: _table_page(pdf, title, table, page, total, fontsize); page += 1
        for image_path in image_pages:
            fig, ax = plt.subplots(figsize=(11, 8.5)); ax.imshow(plt.imread(image_path)); ax.axis("off"); fig.text(.5, .02, f"Page {page} of {total}", ha="center", fontsize=8, color=".4"); pdf.savefig(fig); plt.close(fig); page += 1
    readme = f"""# Cluster-aware stepwise photobleaching results

Start with `01_report/pipeline_audit_and_reanalysis_report.pdf`.

- `02_figures/`: separate seed-level, grouped-cluster, extended-cluster, and method-QC figures.
- `03_tables/complete_seed_ledger.csv`: one row for every Spotiflow 0.4 seed and its complete routing/status.
- `03_tables/seed_level_counts.csv` and `grouped_cluster_counts.csv`: co-primary counting units; do not pool them.
- `03_tables/extended_cluster_counts.csv`: separately qualified extended candidates.
- `03_tables/input_and_acquisition_inventory.csv`, `drift_correction_QC_summary.csv`, and `trace_background_QC_and_attrition.csv`: stage-by-stage QC.
- `03_tables/per_FOV_count_summary.csv` and `quickpbsa_flag_distribution.csv`: explicit per-FOV denominators and every quickPBSA outcome.
- `04_experiments/`: exact experiment folder names with per-FOV detection, traces, and counts.
- `99_technical/`: sensitivity analyses, native quickPBSA files, and provenance.

Failed local growth is a morphology result, not a seed rejection. A `monomer_candidate` is an accepted one-step trace and is not by itself proof of biochemical monomeric state.
"""
    atomic_text(readme, paths["root"] / "00_README_RESULTS.md")
    return {"report": report, "funnel": funnel, "monomer": monomer, "growth_routing": growth_figure, "spotiflow_sensitivity": spot_figure, "mask_comparison": mask_comparison,
            **{f"{stream}_histogram": path for stream, path in histograms.items()}}
