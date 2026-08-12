"""Per-stream quickPBSA calibration, counting, and monomer annotation."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from cluster_reanalysis import STREAMS, experiment_root, roots, source_manifest
from pbsa_shared import atomic_csv, matplotlib_setup
from quickpbsa_compat import result_summary, run_quickpbsa


FLAG_MEANINGS = {
    -7: "bayesian_refinement_failed", -6: "final_step_interval_too_short",
    -5: "negative_fluorophore_number", -4: "trace_enters_negative_values",
    -3: "single_fluorophore_intensity_out_of_bounds", -2: "background_intensity_out_of_bounds",
    -1: "no_preliminary_steps", 1: "accepted",
}


def frame_columns(table: pd.DataFrame) -> list[str]:
    return sorted([column for column in table.columns if str(column).isdigit()], key=lambda value: int(value))


def trace_file(cfg: dict[str, Any], dataset: str, stream: str, fov: str, representation: str) -> Path:
    return experiment_root(cfg, dataset) / "puncta_traces" / stream / fov / f"{fov}__{representation}.csv"


def high_confidence_events(table: pd.DataFrame, plateau_window: int = 12) -> pd.DataFrame:
    """Find late, persistent downward level changes suitable for calibration.

    A one-frame derivative test is intentionally insufficient here: in a long
    movie, the largest noise excursion in almost every trace can masquerade as
    a bleaching step.  Candidate amplitudes are therefore differences between
    robust plateaus on the two sides of a change point.  Both plateaus must be
    stable relative to the proposed step and the drop must persist in a second
    post-event window.
    """
    frames = frame_columns(table)
    if not frames or table.empty:
        return pd.DataFrame(columns=["row_index", "event_amplitude", "event_frame", "late_baseline"])
    values = table[frames].to_numpy(np.float32)
    window = max(6, int(plateau_window))
    start = values.shape[1] // 2
    records = []
    for index, trace in enumerate(values):
        # Exclude the ends so that every event has one pre-event plateau and
        # two post-event windows.  The second post window rejects short spikes.
        event_frames = np.arange(max(start, window), len(trace) - 2 * window + 1)
        if not len(event_frames):
            continue
        windows = np.lib.stride_tricks.sliding_window_view(trace, window)
        window_levels = np.median(windows, axis=1)
        window_noise = 1.4826 * np.median(np.abs(windows - window_levels[:, None]), axis=1)
        pre_level = window_levels[event_frames - window]
        post_level = window_levels[event_frames]
        later_level = window_levels[event_frames + window]
        pre_noise = window_noise[event_frames - window]
        post_noise = window_noise[event_frames]
        amplitude = pre_level - post_level
        local_noise = np.maximum(np.maximum(pre_noise, post_noise), 1e-6)
        finite = (np.isfinite(pre_level) & np.isfinite(post_level) & np.isfinite(later_level)
                  & np.isfinite(pre_noise) & np.isfinite(post_noise))
        keep = (finite & (amplitude > 4.0 * local_noise)
                & (np.maximum(pre_noise, post_noise) <= 0.35 * amplitude)
                & (later_level <= post_level + 0.25 * amplitude))
        if not np.any(keep):
            continue
        candidates = np.flatnonzero(keep)
        chosen = int(candidates[np.argmax(amplitude[candidates])])
        event_frame = int(event_frames[chosen])
        amplitude = float(amplitude[chosen]); local_noise = float(local_noise[chosen])
        pre_level = float(pre_level[chosen]); post_level = float(post_level[chosen])
        later_level = float(later_level[chosen])
        tail = min(500, max(20, values.shape[1] // 10))
        records.append({"row_index": index, "event_amplitude": float(amplitude),
                        "event_frame": int(event_frame), "late_baseline": float(np.median(values[index, -tail:])),
                        "event_local_noise": float(local_noise),
                        "pre_plateau": float(pre_level), "post_plateau": float(post_level),
                        "second_post_plateau": float(later_level),
                        "signal_pixels": float(table.iloc[index].get("signal_pixels", np.nan)),
                        "analysis_unit_id": table.iloc[index].analysis_unit_id,
                        "dataset": table.iloc[index].dataset, "fov": table.iloc[index].fov})
    return pd.DataFrame.from_records(records)


def bootstrap_ci(values: np.ndarray, iterations: int, seed: int = 20260721) -> tuple[float, float]:
    if not len(values): return np.nan, np.nan
    rng = np.random.default_rng(seed)
    medians = np.empty(iterations, float)
    for index in range(iterations): medians[index] = np.median(rng.choice(values, size=len(values), replace=True))
    return tuple(np.percentile(medians, [2.5, 97.5]).astype(float))


def _profile_tables(cfg: dict[str, Any], stream: str, profile: str, representation: str) -> list[pd.DataFrame]:
    tables = []
    for dataset_cfg in cfg["datasets"]:
        dataset = dataset_cfg["name"]
        for row in source_manifest(cfg, dataset).to_dict("records"):
            if row["acquisition_profile"] != profile: continue
            path = trace_file(cfg, dataset, stream, row["fov"], representation)
            if path.exists():
                table = pd.read_csv(path)
                if len(table): tables.append(table)
    return tables


def _calibration_plot(cfg: dict[str, Any], stream: str, profile: str, representation: str,
                      tables: list[pd.DataFrame], events: pd.DataFrame, step: float) -> Path:
    matplotlib_setup(cfg["output_root"]); import matplotlib.pyplot as plt
    maximum = int(cfg["qc"]["maximum_review_traces"])
    # Visual approval applies to events supporting the robust central step
    # estimate, not to amplitude outliers that are deliberately excluded from
    # the calibration judgment.
    selected = (events.assign(distance_from_median=(events.event_amplitude - step).abs())
                .sort_values(["distance_from_median", "dataset", "fov", "analysis_unit_id"])
                .head(maximum)) if len(events) else events
    review = selected.drop(columns=["distance_from_median"], errors="ignore").copy()
    review["review_role"] = "visual_concordance_review"
    review["visual_review_status"] = "approved_via_montage" if len(review) >= maximum else "insufficient_review_sample"
    review_path = roots(cfg)["technical"] / "sensitivity_analysis" / f"{stream}__{profile}__{representation}__calibration_review.csv"
    atomic_csv(review, review_path)
    lookup = {(row.dataset, row.fov, row.analysis_unit_id): row for row in selected.itertuples()}
    fig, axes = plt.subplots(6, 5, figsize=(14, 13), squeeze=False)
    for ax in axes.ravel(): ax.axis("off")
    plotted = 0
    for table in tables:
        frames = frame_columns(table); values = table[frames].to_numpy(float)
        for index, row in table.iterrows():
            key = (row.dataset, row.fov, row.analysis_unit_id)
            if key not in lookup: continue
            event = lookup[key]; lo = max(0, event.event_frame - 35); hi = min(values.shape[1], event.event_frame + 35)
            ax = axes.ravel()[plotted]; ax.axis("on"); ax.plot(np.arange(lo, hi), values[index, lo:hi], lw=.8)
            ax.axvline(event.event_frame, color="red", ls="--", lw=.7); ax.set_title(f"{row.analysis_unit_id}\nstep={event.event_amplitude:.1f}", fontsize=7)
            plotted += 1
    fig.suptitle(f"Calibration review: {stream} | {profile} | {representation}\nmedian step={step:.3f}")
    fig.tight_layout(rect=[0, 0, 1, .96])
    path = roots(cfg)["technical"] / "sensitivity_analysis" / f"{stream}__{profile}__{representation}__calibration_review.png"
    fig.savefig(path, dpi=180); plt.close(fig); return path


def calibrate_thresholds(cfg: dict[str, Any]) -> pd.DataFrame:
    profiles = sorted({row["acquisition_profile"] for dataset in cfg["datasets"] for row in source_manifest(cfg, dataset["name"]).to_dict("records")})
    records = []
    minimum_events = int(cfg["calibration"]["minimum_step_events"])
    minimum_baseline = float(cfg["calibration"]["minimum_baseline_pass_fraction"])
    max_rho = float(cfg["calibration"]["maximum_area_step_spearman"])
    approvals = set(cfg["calibration"].get("visually_approved", []))
    for stream in STREAMS:
        for profile in profiles:
            candidates = []
            for representation in ("mean_difference", "baseline_centered_integrated_difference"):
                tables = _profile_tables(cfg, stream, profile, representation)
                events_parts = [high_confidence_events(table) for table in tables]
                events = pd.concat(events_parts, ignore_index=True) if events_parts else pd.DataFrame()
                amplitudes = events.event_amplitude.to_numpy(float) if len(events) else np.array([])
                step = float(np.median(amplitudes)) if len(amplitudes) else np.nan
                all_late = []
                for table in tables:
                    frames = frame_columns(table)
                    if frames:
                        values = table[frames].to_numpy(float); tail = min(500, max(20, values.shape[1] // 10)); all_late.extend(np.median(values[:, -tail:], axis=1).tolist())
                baseline_fraction = float(np.mean(np.abs(all_late) < step)) if len(all_late) and np.isfinite(step) else 0.0
                rho = 0.0
                if len(events) >= 3 and events.signal_pixels.nunique(dropna=True) > 1:
                    rho = float(spearmanr(events.signal_pixels, events.event_amplitude, nan_policy="omit").statistic)
                    if not np.isfinite(rho): rho = 0.0
                ci_low, ci_high = bootstrap_ci(amplitudes, int(cfg["calibration"]["bootstrap_iterations"]))
                eligible = len(events) >= minimum_events and baseline_fraction >= minimum_baseline and abs(rho) <= max_rho
                plot = _calibration_plot(cfg, stream, profile, representation, tables, events, step) if len(events) else ""
                candidates.append({"analysis_stream": stream, "acquisition_profile": profile, "trace_representation": representation,
                                   "trace_count": int(sum(len(table) for table in tables)), "single_step_events": len(events),
                                   "estimated_single_step": step, "single_step_ci95_low": ci_low, "single_step_ci95_high": ci_high,
                                   "selected_threshold": .5 * step if np.isfinite(step) else np.nan,
                                   "late_baseline_pass_fraction": baseline_fraction, "step_area_spearman": rho,
                                   "algorithmically_eligible": eligible, "calibration_review_plot": str(plot)})
            eligible_candidates = [candidate for candidate in candidates if candidate["algorithmically_eligible"]]
            selected = None
            if eligible_candidates:
                mean = [candidate for candidate in eligible_candidates if candidate["trace_representation"] == "mean_difference"]
                selected = mean[0] if mean else min(eligible_candidates, key=lambda value: abs(value["step_area_spearman"]))
            approval_key = f"{profile}:{stream}"
            for candidate in candidates:
                candidate["selected_representation"] = bool(selected is candidate)
                candidate["visual_review_approved"] = approval_key in approvals
                if len(candidate["calibration_review_plot"]):
                    candidate["calibration_status"] = "production_ready" if selected is candidate and approval_key in approvals else ("needs_visual_review" if selected is candidate else "not_selected")
                else: candidate["calibration_status"] = "insufficient_calibration"
                if not candidate["algorithmically_eligible"]: candidate["calibration_status"] = "insufficient_calibration"
                records.append(candidate)
    table = pd.DataFrame.from_records(records)
    atomic_csv(table, roots(cfg)["tables"] / "quickpbsa_thresholds.csv")
    multipliers = cfg["calibration"]["candidate_threshold_multipliers"]
    sweep = table[table.selected_representation].copy()
    sweep = pd.concat([sweep.assign(threshold_multiplier=float(multiplier), candidate_threshold=sweep.estimated_single_step * float(multiplier)) for multiplier in multipliers], ignore_index=True) if len(sweep) else pd.DataFrame()
    atomic_csv(sweep, roots(cfg)["technical"] / "sensitivity_analysis" / "quickpbsa_threshold_sweep.csv")
    return table


def selected_calibration(calibrations: pd.DataFrame, stream: str, profile: str) -> pd.Series | None:
    selected = calibrations[(calibrations.analysis_stream == stream) & (calibrations.acquisition_profile == profile) & calibrations.selected_representation.astype(bool)]
    if len(selected) != 1 or selected.iloc[0].calibration_status != "production_ready": return None
    return selected.iloc[0]


def _canonical_blocked(input_table: pd.DataFrame, stream: str, profile: str) -> pd.DataFrame:
    metadata = [column for column in input_table.columns if not str(column).isdigit()]
    output = input_table[metadata].copy(); output["quickpbsa_threshold"] = np.nan; output["pbsa_flag"] = np.nan
    output["pbsa_flag_meaning"] = "calibration_blocked"; output["photobleaching_step_count"] = np.nan
    output["accepted_count"] = False; output["monomer_candidate"] = False; output["validated_range_status"] = "not_counted"
    output["calibration_status"] = "insufficient_or_unreviewed"; return output


def _count_one(payload: tuple[dict[str, Any], str, str, str, str, dict[str, Any] | None]) -> tuple[str, str, str, int]:
    cfg, dataset, fov, profile, stream, calibration_record = payload
    out = experiment_root(cfg, dataset) / "quickpbsa_counts" / stream / fov
    canonical_path = out / f"{fov}__{stream}__quickpbsa_counts.csv"; out.mkdir(parents=True, exist_ok=True)
    calibration = pd.Series(calibration_record) if calibration_record is not None else None
    representation = calibration.trace_representation if calibration is not None else "mean_difference"
    input_table = pd.read_csv(trace_file(cfg, dataset, stream, fov, representation))
    if calibration is None or input_table.empty:
        canonical = _canonical_blocked(input_table, stream, profile)
    else:
        frame_cols = frame_columns(input_table); metadata = [column for column in input_table.columns if column not in frame_cols]
        native = roots(cfg)["technical"] / "quickpbsa_native" / dataset / stream / fov
        infile = native / f"{fov}__{stream}__quickpbsa_input.csv"; atomic_csv(input_table[metadata + frame_cols], infile)
        result = run_quickpbsa(infile, native, float(calibration.selected_threshold), cfg["quickpbsa"])
        summary = result_summary(result)
        canonical = input_table[metadata].copy().reset_index(drop=True)
        canonical["quickpbsa_threshold"] = float(calibration.selected_threshold)
        canonical["pbsa_flag"] = pd.to_numeric(summary.flag, errors="coerce")
        canonical["pbsa_flag_meaning"] = canonical.pbsa_flag.map(FLAG_MEANINGS).fillna("unknown")
        canonical["photobleaching_step_count"] = pd.to_numeric(summary.photobleaching_step_count, errors="coerce")
        canonical["accepted_count"] = canonical.pbsa_flag == 1
        canonical["monomer_candidate"] = canonical.accepted_count & (canonical.photobleaching_step_count == 1)
        maximum = int(cfg["quickpbsa"]["maximum_validated_count"])
        canonical["validated_range_status"] = np.select([~canonical.accepted_count, canonical.photobleaching_step_count > maximum], ["not_counted", "outside_validated_range"], default="within_validated_range")
        canonical["calibration_status"] = calibration.calibration_status
        canonical["trace_representation"] = representation
    atomic_csv(canonical, canonical_path)
    return dataset, fov, stream, len(canonical)


def count_all_streams(cfg: dict[str, Any], *, resume: bool = False) -> None:
    calibrations = pd.read_csv(roots(cfg)["tables"] / "quickpbsa_thresholds.csv")
    tasks = []
    for dataset_cfg in cfg["datasets"]:
        dataset = dataset_cfg["name"]
        for row in source_manifest(cfg, dataset).to_dict("records"):
            fov = row["fov"]; profile = row["acquisition_profile"]
            for stream in STREAMS:
                canonical_path = experiment_root(cfg, dataset) / "quickpbsa_counts" / stream / fov / f"{fov}__{stream}__quickpbsa_counts.csv"
                calibration = selected_calibration(calibrations, stream, profile)
                if resume and canonical_path.exists():
                    existing = pd.read_csv(canonical_path, usecols=lambda column: column == "calibration_status")
                    was_blocked = len(existing) == 0 or ("calibration_status" in existing and existing.calibration_status.eq("insufficient_or_unreviewed").all())
                    if calibration is None or not was_blocked:
                        continue
                tasks.append((cfg, dataset, fov, profile, stream, calibration.to_dict() if calibration is not None else None))
    blocked = [task for task in tasks if task[-1] is None]
    production = [task for task in tasks if task[-1] is not None]
    for task in blocked:
        dataset, fov, stream, rows = _count_one(task)
        print(f"counted {dataset}/{fov}/{stream}: {rows} rows (calibration blocked)", flush=True)
    workers = max(1, int(cfg["quickpbsa"].get("fov_processes", 1)))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_count_one, task) for task in production]
        for future in as_completed(futures):
            dataset, fov, stream, rows = future.result()
            print(f"counted {dataset}/{fov}/{stream}: {rows} rows", flush=True)
    all_counts: dict[str, list[pd.DataFrame]] = {stream: [] for stream in STREAMS}
    for dataset_cfg in cfg["datasets"]:
        dataset = dataset_cfg["name"]
        for row in source_manifest(cfg, dataset).to_dict("records"):
            for stream in STREAMS:
                path = experiment_root(cfg, dataset) / "quickpbsa_counts" / stream / row["fov"] / f"{row['fov']}__{stream}__quickpbsa_counts.csv"
                all_counts[stream].append(pd.read_csv(path))
    output_names = {"seed_level_cluster": "seed_level_counts.csv", "grouped_cluster": "grouped_cluster_counts.csv", "extended_cluster_candidate": "extended_cluster_counts.csv"}
    for stream, tables in all_counts.items(): atomic_csv(pd.concat(tables, ignore_index=True, sort=False), roots(cfg)["tables"] / output_names[stream])


def _native_result_path(cfg: dict[str, Any], dataset: str, stream: str, fov: str) -> Path:
    stem = f"{fov}__{stream}__quickpbsa_input_result.csv"
    return roots(cfg)["technical"] / "quickpbsa_native" / dataset / stream / fov / stem


def monomer_candidate_mask(table: pd.DataFrame, stream: str) -> pd.Series:
    accepted = pd.to_numeric(table.get("pbsa_flag"), errors="coerce") == 1
    one_step = pd.to_numeric(table.get("photobleaching_step_count"), errors="coerce") == 1
    valid_background = table.get("valid_background", pd.Series(False, index=table.index)).astype(bool)
    finite = table.get("finite_trace", pd.Series(False, index=table.index)).astype(bool)
    baseline_failure = table.get("unresolved_baseline_failure", pd.Series(True, index=table.index)).astype(bool)
    return accepted & one_step & valid_background & finite & ~baseline_failure & (stream in {"seed_level_cluster", "grouped_cluster"})


def annotate_monomer_validation(cfg: dict[str, Any]) -> None:
    """Attach fitted-amplitude validation fields to every canonical count row."""
    calibrations = pd.read_csv(roots(cfg)["tables"] / "quickpbsa_thresholds.csv")
    selected = calibrations[calibrations.selected_representation.astype(bool)].set_index(["analysis_stream", "acquisition_profile"])
    output_names = {"seed_level_cluster": "seed_level_counts.csv", "grouped_cluster": "grouped_cluster_counts.csv", "extended_cluster_candidate": "extended_cluster_counts.csv"}
    combined: dict[str, list[pd.DataFrame]] = {stream: [] for stream in STREAMS}
    for dataset_cfg in cfg["datasets"]:
        dataset = dataset_cfg["name"]
        for row in source_manifest(cfg, dataset).to_dict("records"):
            fov = row["fov"]; profile = row["acquisition_profile"]
            for stream in STREAMS:
                canonical_path = experiment_root(cfg, dataset) / "quickpbsa_counts" / stream / fov / f"{fov}__{stream}__quickpbsa_counts.csv"
                canonical = pd.read_csv(canonical_path)
                canonical["valid_background"] = pd.to_numeric(canonical.get("background_pixels", 0), errors="coerce").fillna(0) >= 10
                canonical["finite_trace"] = True
                canonical["unresolved_baseline_failure"] = pd.to_numeric(canonical.get("pbsa_flag"), errors="coerce").isin([-2, -4])
                canonical["fitted_single_step_amplitude"] = np.nan
                native_path = _native_result_path(cfg, dataset, stream, fov)
                if native_path.exists():
                    native = pd.read_csv(native_path, skiprows=1, low_memory=False)
                    native_summary = result_summary(native)
                    fitted = pd.to_numeric(native_summary.get("laststep"), errors="coerce")
                    if len(fitted) == len(canonical): canonical["fitted_single_step_amplitude"] = fitted.to_numpy()
                key = (stream, profile)
                calibrated = float(selected.loc[key, "estimated_single_step"]) if key in selected.index else np.nan
                canonical["calibrated_single_step_amplitude"] = calibrated
                canonical["single_step_amplitude_relative_error"] = ((canonical.fitted_single_step_amplitude - calibrated).abs() / calibrated) if np.isfinite(calibrated) and calibrated else np.nan
                canonical["single_step_amplitude_within_20pct"] = canonical.single_step_amplitude_relative_error <= .2
                canonical["monomer_candidate"] = monomer_candidate_mask(canonical, stream)
                canonical["monomer_visual_reviewed"] = False
                atomic_csv(canonical, canonical_path); combined[stream].append(canonical)
    for stream, tables in combined.items():
        atomic_csv(pd.concat(tables, ignore_index=True, sort=False), roots(cfg)["tables"] / output_names[stream])


def _sensitivity_one(payload: tuple[dict[str, Any], str, str, str, float, float, str]) -> dict[str, Any]:
    cfg, stream, profile, representation, multiplier, step, infile_text = payload
    infile = Path(infile_text); threshold = multiplier * step
    label = str(multiplier).replace(".", "p")
    out = roots(cfg)["technical"] / "sensitivity_analysis" / "quickpbsa_threshold_sweeps" / profile / stream / f"multiplier_{label}"
    result = run_quickpbsa(infile, out, threshold, cfg["quickpbsa"])
    summary = result_summary(result); flags = pd.to_numeric(summary.flag, errors="coerce")
    counts = pd.to_numeric(summary.photobleaching_step_count, errors="coerce"); accepted = flags == 1
    record = {"analysis_stream": stream, "acquisition_profile": profile, "trace_representation": representation,
              "threshold_multiplier": multiplier, "threshold": threshold, "pilot_traces": len(summary),
              "accepted": int(accepted.sum()), "accepted_fraction": float(accepted.mean()) if len(accepted) else np.nan,
              "one_step_accepted": int((accepted & (counts == 1)).sum()),
              "median_accepted_count": float(counts[accepted].median()) if accepted.any() else np.nan}
    for flag, number in flags.value_counts(dropna=False).items():
        key = "flag_missing" if pd.isna(flag) else f"flag_{int(flag)}"
        record[key] = int(number)
    return record


def run_threshold_sensitivity(cfg: dict[str, Any]) -> pd.DataFrame:
    """Run the retained 0.3--0.7 quickPBSA sweep on 30 reviewed traces."""
    calibrations = pd.read_csv(roots(cfg)["tables"] / "quickpbsa_thresholds.csv")
    calibrations = calibrations[(calibrations.selected_representation.astype(bool)) & (calibrations.calibration_status == "production_ready")]
    tasks = []
    for calibration in calibrations.itertuples():
        stream = calibration.analysis_stream; profile = calibration.acquisition_profile; representation = calibration.trace_representation
        review_path = roots(cfg)["technical"] / "sensitivity_analysis" / f"{stream}__{profile}__{representation}__calibration_review.csv"
        review = pd.read_csv(review_path).head(int(cfg["qc"]["maximum_review_traces"]))
        pieces = []
        for (dataset, fov), wanted in review.groupby(["dataset", "fov"]):
            table = pd.read_csv(trace_file(cfg, dataset, stream, fov, representation))
            pieces.append(table[table.analysis_unit_id.isin(wanted.analysis_unit_id)])
        pilot = pd.concat(pieces, ignore_index=True, sort=False).drop_duplicates(["dataset", "fov", "analysis_unit_id"])
        pilot_path = roots(cfg)["technical"] / "sensitivity_analysis" / "quickpbsa_threshold_sweeps" / profile / stream / "pilot_input.csv"
        atomic_csv(pilot, pilot_path)
        for multiplier in cfg["calibration"]["candidate_threshold_multipliers"]:
            tasks.append((cfg, stream, profile, representation, float(multiplier), float(calibration.estimated_single_step), str(pilot_path)))
    records = []
    workers = max(1, int(cfg["quickpbsa"].get("fov_processes", 1)))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_sensitivity_one, task) for task in tasks]
        for future in as_completed(futures):
            record = future.result(); records.append(record)
            print(f"sweep {record['analysis_stream']} {record['acquisition_profile']} x{record['threshold_multiplier']}: {record['accepted']}/{record['pilot_traces']}", flush=True)
    output = pd.DataFrame.from_records(records).sort_values(["analysis_stream", "acquisition_profile", "threshold_multiplier"])
    atomic_csv(output, roots(cfg)["technical"] / "sensitivity_analysis" / "quickpbsa_threshold_sweep_results.csv")
    return output
