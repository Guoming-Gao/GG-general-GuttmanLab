"""Compatibility helpers for quickPBSA 2021.0.1 on the current smlm environment."""

from __future__ import annotations

import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def _patch_numpy2_refinement() -> None:
    """Patch two quickPBSA 2021 routines that relied on ragged NumPy arrays.

    NumPy 2 no longer turns a list of unequal-length trace segments into an
    implicitly broadcastable object array.  The calculations below are the
    original quickPBSA equations, evaluated segment-by-segment instead.
    """
    from scipy.special import factorial
    from quickpbsa.steps_refinement import improve_steps as improve_module
    from quickpbsa.steps_refinement import refinement_lowlevel as refinement

    if getattr(refinement, "_pbsa_numpy2_patched", False):
        return

    def posterior_single(data, jpos, means, variances, df, i_in=0, lamb=0.1, gamma0=0.5):
        del i_in  # retained for API compatibility with the published package
        mb = means[0]
        vb = variances[0]
        mf = means[1] - means[0]
        vf = variances[1] - variances[0]
        if vf < 0:
            vf = variances[1]
        occupancy = np.cumsum([0] + list(df))
        if not np.all(occupancy >= 0):
            return float("inf")
        segments = np.split(data, jpos)
        nphi = np.asarray([len(segment) for segment in segments])
        varphi = occupancy * vf + vb
        residual = np.asarray([
            np.sum((segment - level * mf - mb) ** 2)
            for segment, level in zip(segments, occupancy)
        ])
        k_steps = len(jpos)
        multiplicity = sum(map(abs, df))
        _, counts = np.unique(np.abs(df), return_counts=True)
        sic = np.sum(nphi * np.log(varphi) + residual / varphi)
        sic += 2 * (
            -k_steps * np.log(lamb)
            - np.log(factorial(multiplicity - k_steps) * factorial(k_steps) / factorial(multiplicity - 1))
            + np.sum(np.log(factorial(counts)))
        )
        sic += (
            2 * gamma0 * (multiplicity - k_steps + 1) / k_steps
            + np.log(multiplicity - k_steps + 2)
            + np.log(multiplicity - k_steps + 1)
            - np.log(
                multiplicity - k_steps + 2
                - (multiplicity - k_steps + 1) * np.exp(-gamma0 / k_steps)
            )
        )
        return sic

    def add_step(data, jpos, means, variances, steps, lamb=0.1, gamma0=0.5):
        mb = means[0]
        vb = variances[0]
        mf = means[1] - means[0]
        vf = variances[1] - variances[0]
        if vf < 0:
            vf = variances[1]
        size = len(data)
        sicmin = posterior_single(data, jpos, means, variances, steps)
        _, counts = np.unique(np.abs(steps), return_counts=True)
        counts[0] += 1
        step_out = steps
        jpos_out = jpos
        occupancy = np.cumsum(np.hstack((0, steps)))
        step_locations = np.hstack((0, jpos, size))
        segment_sse = np.asarray([
            np.sum((segment - level * mf - mb) ** 2)
            for segment, level in zip(np.split(data, jpos), occupancy)
        ])
        rows = size - jpos[1]
        diffar = np.tile(segment_sse, [rows, 1])
        diffar = np.hstack((diffar, np.zeros([rows, 1])))
        varphi = np.tile(np.hstack((occupancy * vf + vb, 0)), [rows, 1])
        nphi = np.tile(np.hstack((np.diff(step_locations), 0)), [rows, 1])

        sics = []
        step_positions = []
        for step_sign in (-1, 1):
            for index in range(2, len(step_locations) - 1):
                segment_length = int(step_locations[index + 1] - step_locations[index])
                if index < len(step_locations) - 2:
                    tail = data[step_locations[index + 1]:]
                    tail_positions = jpos[index + 1:] - step_locations[index + 1]
                    tail_levels = occupancy[index + 1:] + step_sign
                    tail_sse = np.asarray([
                        np.sum((segment - level * mf - mb) ** 2)
                        for segment, level in zip(np.split(tail, tail_positions), tail_levels)
                    ])
                    row_slice = slice(
                        step_locations[index] - jpos[1],
                        step_locations[index + 1] - jpos[1],
                    )
                    diffar[row_slice, index + 2:] = np.tile(tail_sse, [segment_length, 1])
                    varphi[row_slice, index + 2:] = np.tile(
                        tail_levels * vf + vb, [segment_length, 1]
                    )
                    nphi[row_slice, index + 2:] = np.tile(
                        np.diff(step_locations[index + 1:]), [segment_length, 1]
                    )

                current = np.tile(
                    data[step_locations[index]:step_locations[index + 1]],
                    [segment_length, 1],
                )
                current[np.triu(np.ones([segment_length, segment_length]), 1).astype(bool)] = (
                    occupancy[index] * mf + mb
                )
                row_slice = slice(
                    step_locations[index] - jpos[1],
                    step_locations[index + 1] - jpos[1],
                )
                diffar[row_slice, index] = np.sum(
                    (current - occupancy[index] * mf - mb) ** 2, axis=1
                )

                current = np.tile(
                    data[step_locations[index]:step_locations[index + 1]],
                    [segment_length, 1],
                )
                current[np.tril(np.ones([segment_length, segment_length]), 0).astype(bool)] = (
                    (occupancy[index] + step_sign) * mf + mb
                )
                diffar[row_slice, index + 1] = np.sum(
                    (current - (occupancy[index] + step_sign) * mf - mb) ** 2, axis=1
                )
                varphi[row_slice, index] = occupancy[index] * vf + vb
                varphi[row_slice, index + 1] = (occupancy[index] + step_sign) * vf + vb
                nphi[row_slice, index] = np.arange(segment_length) + 1
                nphi[row_slice, index + 1] = segment_length - np.arange(segment_length) - 1

            multiplicity = np.sum(np.abs(steps)) + 1
            k_steps = len(jpos) + 1
            with np.errstate(divide="ignore", invalid="ignore"):
                sic = np.sum(nphi * np.log(varphi) + diffar / varphi, axis=1)
                sic += 2 * (
                    -k_steps * np.log(lamb)
                    - np.log(factorial(multiplicity - k_steps) * factorial(k_steps) / factorial(multiplicity - 1))
                    + np.sum(np.log(factorial(counts)))
                )
                sic += (
                    2 * gamma0 * (multiplicity - k_steps + 1) / k_steps
                    + np.log(multiplicity - k_steps + 2)
                    + np.log(multiplicity - k_steps + 1)
                    - np.log(
                        multiplicity - k_steps + 2
                        - (multiplicity - k_steps + 1) * np.exp(-gamma0 / k_steps)
                    )
                )
            ignored = np.ravel(np.tile(jpos, [5, 1]).T + np.arange(-2, 3)) - jpos[1]
            ignored = ignored[(ignored > 0) * (ignored < len(sic) - 2)]
            included = np.delete(np.arange(len(sic) - 2), ignored)
            candidate = sic[included]
            sics.append(np.min(candidate))
            step_positions.append(included[np.argmin(candidate)] + jpos[1] + 1)

        if np.nanmin(sics) < sicmin:
            sicmin = np.nanmin(sics)
            new_step = step_positions[np.nanargmin(sics)]
            step_sign = np.nanargmin(sics) * 2 - 1
            insert_at = np.max(np.where(jpos < new_step)[0]) + 1
            step_out = np.hstack((steps[:insert_at], step_sign, steps[insert_at:]))
            jpos_out = sorted(np.hstack((jpos, new_step)))
            variances = [np.var(segment) for segment in np.split(data, jpos_out)]
            means = [np.mean(segment) for segment in np.split(data, jpos_out)]
        return sicmin, jpos_out, step_out, means, variances

    refinement.posterior_single = posterior_single
    refinement.add_step = add_step
    # improve_steps imported add_step directly, so update that module reference too.
    improve_module.add_step = add_step
    refinement._pbsa_numpy2_patched = True


def prepare_quickpbsa() -> Any:
    # quickPBSA 2021.0.1 calls the NumPy alias removed in NumPy 2.
    if not hasattr(np, "float"):
        np.float = float  # type: ignore[attr-defined]
    if not hasattr(np, "in1d"):
        np.in1d = np.isin  # type: ignore[attr-defined]
    try:
        import quickpbsa
    except ImportError as exc:
        raise RuntimeError(
            "quickPBSA is not installed in smlm. Run: conda run -n smlm pip install quickpbsa==2021.0.1"
        ) from exc
    _patch_numpy2_refinement()
    return quickpbsa


def run_quickpbsa(infile: str | Path, outfolder: str | Path, threshold: float, settings: dict[str, Any], *, maxiter: int | None = None) -> pd.DataFrame:
    """Run quickPBSA's published low-level algorithms without its socket-based Manager.

    quickPBSA's high-level file functions require ``multiprocessing.Manager``, which is
    unavailable in some managed/macOS environments. This adapter executes the same KV
    and Bayesian functions in a thread pool and serializes their results with the
    package's own helpers. NumPy-2 compatibility changes preserve the original
    posterior equations while evaluating ragged trace segments explicitly.
    """
    pbsa = prepare_quickpbsa(); outfolder = Path(outfolder); outfolder.mkdir(parents=True, exist_ok=True)
    from quickpbsa.helpers import export_csv
    from quickpbsa.steps_preliminary.helpers import crop_traces, kv_to_json, kvresult_from_json, read_tracedf
    from quickpbsa.steps_preliminary.kv_lowlevel import kv_single, kv_single_fast
    from quickpbsa.steps_refinement.helpers import read_kvresult, step2_to_json, step2result_from_json
    from quickpbsa.steps_refinement.improve_steps import improve_steps_single
    from quickpbsa.steps_refinement.refinement_file import generate_flags

    infile = Path(infile); stem = infile.stem
    kvjson = outfolder / f"{stem}_kvout.json"; kvfile = outfolder / f"{stem}_kvout.csv"
    result_json = outfolder / f"{stem}_result.json"; result_file = outfolder / f"{stem}_result.csv"
    for path in (kvjson, kvfile, result_json, result_file):
        if path.exists(): path.unlink()

    preliminary = {"norm": 1, "max_memory": 4.0, "crop": True, "bg_frames": 500}
    preliminary.update(settings.get("preliminary", {})); max_iterations = int(maxiter or settings["maxiter"])
    traces, basedf, comment, parameters, n_traces, _ = read_tracedf(str(infile))
    parameters.update({"KV_threshold": float(threshold), "maxiter": max_iterations,
                       "norm": preliminary["norm"], "crop": preliminary["crop"], "bg_frames": preliminary["bg_frames"]})
    traces = np.fliplr(traces) / float(preliminary["norm"])
    crop_index = crop_traces(traces, float(threshold) / 2, int(preliminary["bg_frames"])) if preliminary["crop"] else np.zeros(n_traces, dtype=int)

    def preliminary_one(index: int) -> tuple[int, tuple[Any, ...], float]:
        started = time.time(); trace = traces[index, crop_index[index]:]
        algorithm = kv_single_fast if len(trace) ** 2 * 4 / 1e9 < float(preliminary["max_memory"]) else kv_single
        return index, algorithm(trace, float(threshold), max_iterations), time.time() - started

    with ThreadPoolExecutor(max_workers=max(1, int(settings.get("num_cores", 1)) - 1)) as pool:
        preliminary_results = list(pool.map(preliminary_one, range(n_traces)))
    for index, result, elapsed in preliminary_results:
        steppos, means, variances, posbysic, kv_iter = result
        steppos = steppos + crop_index[index]; posbysic = posbysic + crop_index[index]
        if len(steppos):
            means_sub = means - means[0]; fluors_intensity = np.round(means_sub / means_sub[1]).astype(int)
            fluors_kv = np.cumsum(np.hstack((0, np.sign(np.diff(means)))))
        else:
            fluors_intensity = np.array([0]); fluors_kv = np.array([0])
        kv_to_json(str(kvjson), parameters, index, steppos, means, variances, posbysic,
                   fluors_intensity, fluors_kv, 0, elapsed, kv_iter)
    kvout, kvarray = kvresult_from_json(str(kvjson), traces, basedf)
    export_csv(kvout, kvarray, str(kvfile), parameters, comment)

    filter_settings = {"subtracted": True, "percentile_step": 90, "length_laststep": 20}
    filter_settings.update(settings.get("filter", {}))
    refinement = {"maxmult": 5, "lambda": 0.1, "gamma0": 0.5, "multstep_fraction": 0.5,
                  "nonegatives": False, "mult_threshold": 1.0, "maxadded": 10,
                  "splitcomb": 30000, "combcutoff": 2000000}
    refinement.update(settings.get("refinement", {}))
    kvresult, original_traces, _, result_comment, result_parameters, n_traces, _ = read_kvresult(str(kvfile))
    original_traces = np.fliplr(original_traces)
    kvdata = json.loads(kvjson.read_text())
    flags, avg_laststep = generate_flags(kvdata, float(threshold), bool(filter_settings["subtracted"]),
                                         filter_settings["percentile_step"], int(filter_settings["length_laststep"]))
    result_parameters.update({"subtracted": filter_settings["subtracted"], "percentile_step": filter_settings["percentile_step"],
                              "length_laststep": filter_settings["length_laststep"], **refinement})

    eligible = [index for index in range(n_traces) if flags[index] == 0]
    def refine_one(index: int) -> tuple[int, tuple[Any, ...], float]:
        location = int(np.where(np.asarray(kvdata["trace_index"]) == index)[0][0]); started = time.time()
        result = improve_steps_single(
            original_traces[index], np.asarray(kvdata["steppos"][location]), np.asarray(kvdata["means"][location]),
            np.asarray(kvdata["variances"][location]), np.asarray(kvdata["posbysic"][location]),
            int(refinement["maxmult"]), float(refinement["lambda"]), float(refinement["gamma0"]),
            float(refinement["multstep_fraction"]), bool(refinement["nonegatives"]), float(refinement["mult_threshold"]),
            int(refinement["maxadded"]), int(refinement["splitcomb"]), int(refinement["combcutoff"]),
        )
        return index, result, time.time() - started
    with ThreadPoolExecutor(max_workers=max(1, int(settings.get("num_cores", 1)) - 1)) as pool:
        refined = list(pool.map(refine_one, eligible))
    if not refined:
        result_json.write_text(json.dumps({"trace_index": [], "steppos": [], "steps": [], "step2_time": [], "sic_final": [], "flag": [], "parameters": list(result_parameters.values())}))
    for index, result, elapsed in refined:
        success, sic_final, steppos, steps = result
        step2_to_json(str(result_json), result_parameters, index, steppos, steps, elapsed, sic_final, 1 if success else -7)
    result_out, result_array = step2result_from_json(str(result_json), str(kvfile), flags, avg_laststep)
    return export_csv(result_out, result_array, str(result_file), result_parameters, result_comment)


def result_summary(result: pd.DataFrame) -> pd.DataFrame:
    full = result[result["type"] == "fluors_full"].copy().reset_index(drop=True)
    frame_columns = [column for column in full.columns if str(column).isdigit()]
    if frame_columns:
        counts = full[frame_columns].apply(pd.to_numeric, errors="coerce").max(axis=1)
    else:
        counts = pd.Series(np.nan, index=full.index)
    metadata_columns = [column for column in full.columns if column not in frame_columns]
    output = full[metadata_columns].copy()
    # quickPBSA leaves the ``fluors_full`` row at zero when a trace is rejected
    # by its quality filters.  Zero is a plausible biological result, so expose
    # rejected values as missing and preserve the package flag for diagnosis.
    accepted = pd.to_numeric(output.get("flag"), errors="coerce") == 1
    output["photobleaching_step_count"] = counts.where(accepted)
    return output
