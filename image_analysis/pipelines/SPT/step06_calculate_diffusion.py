"""Step 6: calculate historical AIO diffusion metrics and pooled saSPT by condition."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from scipy import stats

from spt_shared import atomic_csv, output_dirs, stage_timer
from step01_inspect_inputs import accepted_rows, inspect_inputs
from step02_preprocess_videos import channel_specs
from step05_link_trajectories import tracking_paths


ANGLE_BINS = np.linspace(0, 180, 7).astype(int)
ANGLE_COLUMNS = [f"({ANGLE_BINS[index]},{ANGLE_BINS[index+1]}]" for index in range(6)]
AIO_COLUMNS = [
    "trackID", "list_of_t", "list_of_x", "list_of_y", "N_steps", "displacement_nm",
    "mean_stepsize_nm", "max_d_anytwo_nm", "mean_x_pxl", "mean_y_pxl",
    "mean_spot_intensity_max_in_track", "list_of_MSD_um2", "list_of_tau_s",
    "linear_fit_slope", "linear_fit_R2", "linear_fit_sigma", "linear_fit_D_um2s",
    "linear_fit_log10D", "loglog_fit_R2", "loglog_fit_log10D", "alpha", "list_of_angles",
] + ANGLE_COLUMNS
PROPAGATED_COLUMNS = [
    "dataset", "condition", "fov", "channel", "cell_id", "assignment_confidence",
    "assignment_ambiguous", "trajectory_uid",
]


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_") or "condition"


def _msd(track: pd.DataFrame, lags: np.ndarray) -> np.ndarray:
    ordered = track.sort_values("t"); first, last = int(ordered["t"].min()), int(ordered["t"].max())
    full = ordered.set_index(ordered["t"].astype(int)).reindex(range(first, last+1))
    x, y = full["x"].to_numpy(float), full["y"].to_numpy(float)
    return np.asarray([np.nanmean((x[lag:]-x[:-lag])**2 + (y[lag:]-y[:-lag])**2) for lag in lags])


def _angles(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    vector = np.stack([np.diff(x), np.diff(y)]); headings = np.angle(vector[0] + 1j*vector[1], deg=True)
    result = np.diff(headings); result[result < -180] += 360; result[result > 180] -= 360
    return result


def _constant_metadata(track: pd.DataFrame) -> dict[str, object]:
    result = {}
    for column in PROPAGATED_COLUMNS:
        if column in track and track[column].notna().any(): result[column] = track[column].dropna().iloc[0]
    return result


def calculate_aio_table(tracks: pd.DataFrame, *, frame_interval_s: float, pixel_size_um: float = 0.117,
                        min_track_length: int = 5) -> pd.DataFrame:
    required = {"trackID", "x", "y", "t", "meanIntensity"}; missing = required-set(tracks)
    if missing: raise ValueError(f"Trajectory table missing {sorted(missing)}")
    tracks = tracks.copy()
    for column in required: tracks[column] = pd.to_numeric(tracks[column], errors="raise")
    rows = []
    for track_id, track in tracks.groupby("trackID", sort=True):
        track = track.sort_values("t")
        if len(track) < min_track_length: continue
        x, y = track["x"].to_numpy(float), track["y"].to_numpy(float)
        if np.unique(x, return_counts=True)[1].max() > 2 or np.unique(y, return_counts=True)[1].max() > 2: continue
        lags = np.arange(1, len(track)); msd = _msd(track, lags)*pixel_size_um**2; tau = lags*frame_interval_s
        fit_count = min(max(3, round(len(lags)/2)), len(lags)); fit_tau, fit_msd = tau[:fit_count], msd[:fit_count]
        slope, intercept, linear_r, _, _ = stats.linregress(fit_tau, fit_msd)
        if slope > 0:
            diffusion = slope/(8/3); log10_diffusion = np.log10(diffusion)
            sigma_nm = np.sqrt(intercept/4)*1000 if intercept >= 0 else np.nan
        else: diffusion = log10_diffusion = sigma_nm = np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            alpha, log_intercept, log_r, _, _ = stats.linregress(np.log10(fit_tau), np.log10(fit_msd))
        angles = _angles(x, y); densities, _ = np.histogram(np.abs(angles), ANGLE_BINS, density=True)
        fractions = densities*(ANGLE_BINS[1]-ANGLE_BINS[0]); steps = np.hypot(np.diff(x), np.diff(y))
        pairwise = np.hypot(x[:, None]-x[None, :], y[:, None]-y[None, :])
        rows.append({
            "trackID": track_id, "list_of_t": track["t"].astype(int).tolist(), "list_of_x": x.tolist(),
            "list_of_y": y.tolist(), "N_steps": len(track),
            "displacement_nm": np.hypot(x[-1]-x[0], y[-1]-y[0])*pixel_size_um*1000,
            "mean_stepsize_nm": float(steps.mean()*pixel_size_um*1000),
            "max_d_anytwo_nm": float(pairwise.max()*pixel_size_um*1000),
            "mean_x_pxl": float(x.mean()), "mean_y_pxl": float(y.mean()),
            "mean_spot_intensity_max_in_track": float(track["meanIntensity"].max()),
            "list_of_MSD_um2": msd.tolist(), "list_of_tau_s": tau.tolist(), "linear_fit_slope": slope,
            "linear_fit_R2": linear_r**2, "linear_fit_sigma": sigma_nm, "linear_fit_D_um2s": diffusion,
            "linear_fit_log10D": log10_diffusion, "loglog_fit_R2": log_r**2,
            "loglog_fit_log10D": log_intercept-np.log10(4), "alpha": alpha,
            "list_of_angles": angles.tolist(), **dict(zip(ANGLE_COLUMNS, fractions)), **_constant_metadata(track),
        })
    columns = AIO_COLUMNS + [column for column in PROPAGATED_COLUMNS if any(column in row for row in rows)]
    return pd.DataFrame.from_records(rows, columns=columns)


def _parse_array(value: object) -> np.ndarray:
    return np.fromstring(str(value).strip()[1:-1], sep=",", dtype=float)


def reformat_aio_for_saspt(aio: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for trajectory, (_, row) in enumerate(aio.iterrows()):
        x, y, frame = _parse_array(row["list_of_x"]), _parse_array(row["list_of_y"]), _parse_array(row["list_of_t"])
        rows.append(pd.DataFrame({"x": x, "y": y, "trajectory": trajectory, "frame": frame}, dtype=float))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["x", "y", "trajectory", "frame"], dtype=float)


def run_saspt(aio: pd.DataFrame, *, frame_interval_s: float, pixel_size_um: float,
              focal_depth_um: float) -> pd.DataFrame:
    from saspt import RBME, StateArray
    filtered = aio[(aio["mean_stepsize_nm"] > 30) & (aio["alpha"] > 0.7)]
    detections = reformat_aio_for_saspt(filtered)
    if detections.empty: return pd.DataFrame(columns=["diff_coef", "loc_error", "mean_posterior_occupation"])
    state = StateArray.from_detections(
        detections, likelihood_type=RBME, pixel_size_um=pixel_size_um,
        frame_interval=frame_interval_s, focal_depth=focal_depth_um, progress_bar=True,
    )
    return state.occupations_dataframe


def analyze_aio_files(paths: Iterable[str | Path], *, frame_interval_s: float, pixel_size_um: float = 0.117,
                      output_dir: str | Path | None = None) -> list[Path]:
    outputs = []
    for value in paths:
        path = Path(value); table = calculate_aio_table(pd.read_csv(path), frame_interval_s=frame_interval_s, pixel_size_um=pixel_size_um)
        output = (Path(output_dir) if output_dir else path.parent)/f"SPT_results_AIO-{path.name}"
        atomic_csv(table, output); outputs.append(output)
    return outputs


def concat_aio_files(paths: Iterable[str | Path], output: str | Path) -> Path:
    tables = []
    for value in paths:
        path = Path(value); table = pd.read_csv(path); table.insert(0, "filename", path.name); tables.append(table)
    if not tables: raise ValueError("No AIO inputs")
    atomic_csv(pd.concat(tables, ignore_index=True, sort=False), output); return Path(output)


def diffusion_stage(cfg: dict[str, Any], manifest: pd.DataFrame, *, max_files: int | None = None) -> dict[str, list[str]]:
    dirs = output_dirs(cfg["output_dir"]); per_fov_dir = dirs["analysis"]/"per_fov"
    condition_dir, saspt_dir = dirs["analysis"]/"by_condition", dirs["analysis"]/"saspt"
    rows = accepted_rows(manifest, max_files); outputs: dict[str, list[str]] = {"per_fov": [], "conditions": [], "saspt": []}
    grouped: dict[tuple[str, str], list[tuple[dict[str, Any], Path]]] = {}
    with stage_timer(cfg, "06_diffusion", {"files": len(rows), "saspt": bool(cfg["analysis"]["run_saspt"])}):
        for row in rows:
            for spec in channel_specs(cfg):
                if spec["role"] != "spt": continue
                tracks = pd.read_csv(tracking_paths(cfg["output_dir"], row["fov"], spec["name"])["canonical"])
                aio = calculate_aio_table(tracks, frame_interval_s=float(row["frame_interval_s"]), pixel_size_um=float(row["pixel_size_um"]))
                path = per_fov_dir/f"SPT_results_AIO-{row['fov']}__{spec['name']}.csv"; atomic_csv(aio, path)
                outputs["per_fov"].append(str(path)); grouped.setdefault((row["condition"], spec["name"]), []).append((row, path))
        for (condition, channel), members in grouped.items():
            intervals = np.asarray([float(row["frame_interval_s"]) for row, _ in members]); pixels = np.asarray([float(row["pixel_size_um"]) for row, _ in members])
            if np.ptp(intervals) > float(cfg["frame_interval_tolerance_s"]): raise ValueError(f"Mixed frame intervals in {condition}/{channel}")
            if np.ptp(pixels) > float(cfg["pixel_size_tolerance_um"]): raise ValueError(f"Mixed pixel sizes in {condition}/{channel}")
            tables = []
            for row, path in members:
                table = pd.read_csv(path); table.insert(0, "filename", row["filename"]); tables.append(table)
            combined = pd.concat(tables, ignore_index=True, sort=False); stem = f"{safe_name(condition)}__{safe_name(channel)}"
            condition_path = condition_dir/f"SPT_results_AIO_concat-{stem}.csv"; atomic_csv(combined, condition_path)
            outputs["conditions"].append(str(condition_path))
            if cfg["analysis"]["run_saspt"]:
                # AIO uses each acquisition's validated metadata value.  The
                # historical pooled saSPT contract is explicitly 0.117 µm/px.
                saspt = run_saspt(
                    combined,
                    frame_interval_s=float(intervals[0]),
                    pixel_size_um=float(cfg["expected_pixel_size_um"]),
                    focal_depth_um=float(cfg["analysis"]["focal_depth_um"]),
                )
                saspt.insert(0, "channel", channel); saspt.insert(0, "condition", condition)
                saspt_path = saspt_dir/f"saSPT-pooled-mobile-{stem}.csv"; atomic_csv(saspt, saspt_path); outputs["saspt"].append(str(saspt_path))
    return outputs


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", required=True); parser.add_argument("--max-files", type=int)
    args = parser.parse_args(argv); cfg, manifest = inspect_inputs(args.config); result = diffusion_stage(cfg, manifest, max_files=args.max_files)
    for values in result.values(): print("\n".join(values))


if __name__ == "__main__":
    main()
