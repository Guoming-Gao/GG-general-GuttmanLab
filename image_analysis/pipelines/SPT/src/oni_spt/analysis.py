from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats

from .utils import atomic_csv


ANGLE_BINS = np.linspace(0, 180, 7).astype(int)
ANGLE_COLUMNS = [f"({ANGLE_BINS[i]},{ANGLE_BINS[i + 1]}]" for i in range(6)]
AIO_COLUMNS = [
    "trackID",
    "list_of_t",
    "list_of_x",
    "list_of_y",
    "N_steps",
    "displacement_nm",
    "mean_stepsize_nm",
    "max_d_anytwo_nm",
    "mean_x_pxl",
    "mean_y_pxl",
    "mean_spot_intensity_max_in_track",
    "list_of_MSD_um2",
    "list_of_tau_s",
    "linear_fit_slope",
    "linear_fit_R2",
    "linear_fit_sigma",
    "linear_fit_D_um2s",
    "linear_fit_log10D",
    "loglog_fit_R2",
    "loglog_fit_log10D",
    "alpha",
    "list_of_angles",
] + ANGLE_COLUMNS
PROPAGATED_COLUMNS = [
    "dataset",
    "fov",
    "channel",
    "cell_id",
    "assignment_confidence",
    "assignment_ambiguous",
    "trajectory_uid",
]


def _msd(track: pd.DataFrame, lags: np.ndarray) -> np.ndarray:
    ordered = track.sort_values("t")
    first = int(ordered["t"].min())
    last = int(ordered["t"].max())
    full = ordered.set_index(ordered["t"].astype(int)).reindex(range(first, last + 1))
    x = full["x"].to_numpy(float)
    y = full["y"].to_numpy(float)
    return np.asarray(
        [np.nanmean((x[lag:] - x[:-lag]) ** 2 + (y[lag:] - y[:-lag]) ** 2) for lag in lags]
    )


def _angles(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    vector = np.stack([np.diff(x), np.diff(y)])
    headings = np.angle(vector[0] + 1j * vector[1], deg=True)
    angles = np.diff(headings)
    angles[angles < -180] += 360
    angles[angles > 180] -= 360
    return angles


def _constant_metadata(track: pd.DataFrame) -> dict[str, object]:
    values: dict[str, object] = {}
    for column in PROPAGATED_COLUMNS:
        if column in track.columns and track[column].notna().any():
            unique = track[column].dropna().unique()
            values[column] = unique[0] if len(unique) else np.nan
    return values


def calculate_aio_table(
    tracks: pd.DataFrame,
    *,
    frame_interval_s: float,
    pixel_size_um: float = 0.117,
    min_track_length: int = 5,
) -> pd.DataFrame:
    required = {"trackID", "x", "y", "t", "meanIntensity"}
    missing = required - set(tracks.columns)
    if missing:
        raise ValueError(f"Trajectory table missing columns: {sorted(missing)}")
    numeric = ["trackID", "x", "y", "t", "meanIntensity"]
    tracks = tracks.copy()
    for column in numeric:
        tracks[column] = pd.to_numeric(tracks[column], errors="raise")
    rows: list[dict[str, object]] = []
    for track_id, track in tracks.groupby("trackID", sort=True):
        track = track.sort_values("t")
        if len(track) < min_track_length:
            continue
        x = track["x"].to_numpy(float)
        y = track["y"].to_numpy(float)
        # Preserve the original dead-pixel rejection behavior.
        if np.unique(x, return_counts=True)[1].max() > 2 or np.unique(y, return_counts=True)[1].max() > 2:
            continue
        lags = np.arange(1, len(track))
        msd_um2 = _msd(track, lags) * pixel_size_um**2
        tau_s = lags * frame_interval_s
        n_fit = max(3, round(len(lags) / 2))
        n_fit = min(n_fit, len(lags))
        fit_tau = tau_s[:n_fit]
        fit_msd = msd_um2[:n_fit]
        slope, intercept, r_value, _, _ = stats.linregress(fit_tau, fit_msd)
        if slope > 0:
            diffusion = slope / (8 / 3)
            log10_diffusion = np.log10(diffusion)
            sigma_nm = np.sqrt(intercept / 4) * 1000 if intercept >= 0 else np.nan
        else:
            diffusion = log10_diffusion = sigma_nm = np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            log_slope, log_intercept, log_r, _, _ = stats.linregress(
                np.log10(fit_tau), np.log10(fit_msd)
            )
        angles = _angles(x, y)
        densities, _ = np.histogram(np.abs(angles), ANGLE_BINS, density=True)
        fractions = densities * (ANGLE_BINS[1] - ANGLE_BINS[0])
        steps = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)
        pairwise = np.sqrt((x[:, None] - x[None, :]) ** 2 + (y[:, None] - y[None, :]) ** 2)
        row = {
            "trackID": track_id,
            "list_of_t": track["t"].astype(int).tolist(),
            "list_of_x": x.tolist(),
            "list_of_y": y.tolist(),
            "N_steps": len(track),
            "displacement_nm": np.hypot(x[-1] - x[0], y[-1] - y[0]) * pixel_size_um * 1000,
            "mean_stepsize_nm": float(steps.mean() * pixel_size_um * 1000),
            "max_d_anytwo_nm": float(pairwise.max() * pixel_size_um * 1000),
            "mean_x_pxl": float(x.mean()),
            "mean_y_pxl": float(y.mean()),
            "mean_spot_intensity_max_in_track": float(track["meanIntensity"].max()),
            "list_of_MSD_um2": msd_um2.tolist(),
            "list_of_tau_s": tau_s.tolist(),
            "linear_fit_slope": slope,
            "linear_fit_R2": r_value**2,
            "linear_fit_sigma": sigma_nm,
            "linear_fit_D_um2s": diffusion,
            "linear_fit_log10D": log10_diffusion,
            "loglog_fit_R2": log_r**2,
            "loglog_fit_log10D": log_intercept - np.log10(4),
            "alpha": log_slope,
            "list_of_angles": angles.tolist(),
            **dict(zip(ANGLE_COLUMNS, fractions)),
            **_constant_metadata(track),
        }
        rows.append(row)
    columns = AIO_COLUMNS + [column for column in PROPAGATED_COLUMNS if any(column in row for row in rows)]
    return pd.DataFrame.from_records(rows, columns=columns)


def analyze_aio_files(
    paths: Iterable[str | Path],
    *,
    frame_interval_s: float,
    pixel_size_um: float = 0.117,
    output_dir: str | Path | None = None,
) -> list[Path]:
    outputs = []
    for value in paths:
        path = Path(value)
        table = calculate_aio_table(
            pd.read_csv(path), frame_interval_s=frame_interval_s, pixel_size_um=pixel_size_um
        )
        destination = Path(output_dir) if output_dir else path.parent
        output = destination / f"SPT_results_AIO-{path.name}"
        atomic_csv(table, output)
        outputs.append(output)
    return outputs


def concat_aio_files(paths: Iterable[str | Path], output: str | Path) -> Path:
    tables = []
    for value in paths:
        path = Path(value)
        table = pd.read_csv(path)
        table.insert(0, "filename", path.name)
        tables.append(table)
    if not tables:
        raise ValueError("No AIO files provided")
    destination = Path(output)
    atomic_csv(pd.concat(tables, ignore_index=True, sort=False), destination)
    return destination


def _parse_array(value: object) -> np.ndarray:
    text = str(value).strip()[1:-1]
    return np.fromstring(text, sep=",", dtype=float)


def reformat_aio_for_saspt(aio: pd.DataFrame) -> pd.DataFrame:
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    frames: list[np.ndarray] = []
    trajectories: list[np.ndarray] = []
    for trajectory, (_, row) in enumerate(aio.iterrows()):
        x = _parse_array(row["list_of_x"])
        y = _parse_array(row["list_of_y"])
        frame = _parse_array(row["list_of_t"])
        xs.append(x)
        ys.append(y)
        frames.append(frame)
        trajectories.append(np.full(len(frame), trajectory, dtype=float))
    if not xs:
        return pd.DataFrame(columns=["x", "y", "trajectory", "frame"], dtype=float)
    return pd.DataFrame(
        {
            "x": np.concatenate(xs),
            "y": np.concatenate(ys),
            "trajectory": np.concatenate(trajectories),
            "frame": np.concatenate(frames),
        },
        dtype=float,
    )


def analyze_saspt_file(
    path: str | Path,
    *,
    frame_interval_s: float,
    pixel_size_um: float = 0.117,
    focal_depth_um: float = 0.7,
    output: str | Path | None = None,
) -> Path:
    from saspt import RBME, StateArray

    path = Path(path)
    aio = pd.read_csv(path)
    filtered = aio[(aio["mean_stepsize_nm"] > 30) & (aio["alpha"] > 0.7)]
    detections = reformat_aio_for_saspt(filtered)
    if detections.empty:
        raise ValueError("No mobile, non-confined trajectories remain for saSPT")
    state = StateArray.from_detections(
        detections,
        likelihood_type=RBME,
        pixel_size_um=pixel_size_um,
        frame_interval=frame_interval_s,
        focal_depth=focal_depth_um,
        progress_bar=True,
    )
    destination = Path(output) if output else path.with_name(f"saSPT-pooled-mobile-{path.name}")
    atomic_csv(state.occupations_dataframe, destination)
    return destination
