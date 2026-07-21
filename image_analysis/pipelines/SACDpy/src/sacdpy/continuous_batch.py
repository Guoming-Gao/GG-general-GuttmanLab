from __future__ import annotations

import json
import platform
import shutil
import sys
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from time import perf_counter
from typing import Callable, Iterable

import numpy as np
import tifffile

from .discrete_timelapse import (
    TimelapseGroup,
    apply_incomplete_grid_policy,
    discover_folder_plans,
    discover_fov_folders,
    infer_na,
    infer_pixel_nm,
    infer_time_interval_s,
    infer_z_spacing_um,
    output_paths_for_group,
    write_timelapse_tiff,
)
from .params import SACDParams
from .reconstruction import reconstruct
from .tiffio import read_tiff_stack


@dataclass(frozen=True)
class FOVPlan:
    dataset_name: str
    raw_root: Path
    result_dir: Path
    fov_folder: Path
    relative_fov: str
    group: TimelapseGroup
    stack_output: Path
    mip_output: Path
    pixel_nm: float
    na: float
    z_spacing_um: float | None
    time_interval_s: float | None

    @property
    def movie_count(self) -> int:
        return len(self.group.files)


@dataclass(frozen=True)
class BatchPlan:
    fovs: tuple[FOVPlan, ...]
    exclusions: tuple[dict[str, str], ...]
    legacy_datasets: tuple[dict, ...]

    @property
    def movie_count(self) -> int:
        return sum(fov.movie_count for fov in self.fovs)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_config(path: str | Path) -> dict:
    with Path(path).open() as handle:
        config = json.load(handle)
    if config.get("version") != 1:
        raise ValueError("Batch config version must be 1")
    if not config.get("datasets"):
        raise ValueError("Batch config must contain at least one dataset")
    return config


def build_batch_plan(config: dict) -> BatchPlan:
    defaults = config["processing"]
    output_root = Path(config["output_root"])
    fovs: list[FOVPlan] = []
    exclusions: list[dict[str, str]] = []
    legacy: list[dict] = []
    planned_paths: set[Path] = set()

    for dataset in config["datasets"]:
        raw_root = Path(dataset["raw_root"])
        result_dir = output_root / dataset.get("result_name", raw_root.name)
        if dataset.get("mode", "process") == "legacy_provenance":
            for relative, reason in dataset.get("exclude_fovs", {}).items():
                exclusions.append(_exclusion(raw_root, relative, reason))
            legacy.append({**dataset, "result_dir": str(result_dir)})
            continue

        discovery = discover_fov_folders(
            raw_root,
            position_folder=defaults.get("position_folder", "pos_0"),
            include_folder_patterns=dataset.get("include_folder_patterns", ()),
            exclude_folder_patterns=(),
        )
        selected = set(dataset.get("selected_fov_folders", ()))
        excluded = dataset.get("exclude_fovs", {})
        aliases = dataset.get("output_prefix_aliases", {})

        for folder in discovery.discovered:
            relative = folder.relative_to(raw_root).as_posix()
            if selected and relative not in selected:
                exclusions.append(_exclusion(raw_root, relative, "not selected by explicit FOV list"))
                continue
            if relative in excluded:
                exclusions.append(_exclusion(raw_root, relative, excluded[relative]))
                continue

            plans = discover_folder_plans(
                [folder / defaults.get("position_folder", "pos_0")],
                glob_pattern=defaults.get("glob_pattern", "*.tif"),
                channel_name=str(defaults.get("channel_name", "647")),
                recursive=False,
                output_dir=result_dir,
            )
            if len(plans) != 1 or len(plans[0].groups) != 1:
                raise ValueError(f"Expected one ONI group in {folder}, found {len(plans[0].groups)}")
            policy_result = apply_incomplete_grid_policy(
                plans[0].groups[0], defaults.get("incomplete_grid_policy", "error")
            )
            if policy_result.group is None:
                exclusions.append(_exclusion(raw_root, relative, policy_result.status))
                continue
            group = policy_result.group
            if relative in aliases:
                group = replace(group, prefix=aliases[relative])

            first_file = next(iter(group.files.values()))
            pixel_nm = defaults.get("pixel_nm")
            if pixel_nm is None:
                pixel_nm = infer_pixel_nm(first_file, defaults.get("fallback_pixel_nm"))
            na = defaults.get("na")
            if na is None:
                na = infer_na(first_file, defaults.get("fallback_na"))
            if pixel_nm is None or na is None:
                raise ValueError(f"Could not resolve pixel size or NA for {folder}")

            stack_output, mip_output = output_paths_for_group(group, result_dir)
            for path in (stack_output, mip_output):
                if path in planned_paths:
                    raise ValueError(f"Output path collision: {path}")
                planned_paths.add(path)
            fovs.append(
                FOVPlan(
                    dataset_name=raw_root.name,
                    raw_root=raw_root,
                    result_dir=result_dir,
                    fov_folder=folder,
                    relative_fov=relative,
                    group=group,
                    stack_output=stack_output,
                    mip_output=mip_output,
                    pixel_nm=float(pixel_nm),
                    na=float(na),
                    z_spacing_um=infer_z_spacing_um(group),
                    time_interval_s=infer_time_interval_s(
                        group,
                        source=defaults.get("time_interval_source", "metadata"),
                        nominal_time_interval_s=defaults.get("nominal_time_interval_s"),
                    ),
                )
            )
    return BatchPlan(tuple(fovs), tuple(exclusions), tuple(legacy))


def preflight_summary(plan: BatchPlan) -> dict:
    by_dataset: dict[str, dict[str, int]] = {}
    for fov in plan.fovs:
        item = by_dataset.setdefault(fov.dataset_name, {"fovs": 0, "movies": 0, "outputs": 0})
        item["fovs"] += 1
        item["movies"] += fov.movie_count
        item["outputs"] += 2
    return {
        "fovs": len(plan.fovs),
        "movies": plan.movie_count,
        "outputs": len(plan.fovs) * 2,
        "exclusions": len(plan.exclusions),
        "legacy_datasets": len(plan.legacy_datasets),
        "by_dataset": by_dataset,
    }


def run_batch(
    config: dict,
    config_path: str | Path,
    *,
    max_new_fovs: int | None = None,
    progress_callback: Callable[[dict], None] | None = None,
) -> list[dict]:
    plan = build_batch_plan(config)
    repo_root = Path(config.get("repo_root", Path.cwd()))
    config_path = Path(config_path)
    _prepare_provenance(config, config_path, plan, repo_root)

    defaults = config["processing"]
    initial_seconds_per_movie = float(config.get("initial_seconds_per_movie", 10.09))
    completed_movies = 0
    elapsed_processing = 0.0
    newly_written = 0
    results: list[dict] = []

    for index, fov in enumerate(plan.fovs, start=1):
        manifest = _load_manifest(fov.result_dir)
        previous = next(
            (item for item in manifest["fovs"] if item.get("relative_fov") == fov.relative_fov), None
        )
        if previous and previous.get("status") in {"written", "skipped_existing"}:
            completed_movies += fov.movie_count
            result = {**previous, "status": "resumed_manifest"}
            results.append(result)
            _emit(progress_callback, result)
            continue
        if max_new_fovs is not None and newly_written >= max_new_fovs:
            break

        try:
            result = process_fov(fov, defaults)
            newly_written += result["status"] == "written"
        except Exception as exc:
            result = _result_record(fov, status="failed", runtime_s=0.0, error=repr(exc))
        results.append(result)
        _upsert_manifest(fov.result_dir, result)

        completed_movies += fov.movie_count
        elapsed_processing += float(result.get("runtime_s", 0.0))
        measured_movies = sum(
            item["movie_count"] for item in results if item.get("status") == "written"
        )
        seconds_per_movie = (
            elapsed_processing / measured_movies if measured_movies else initial_seconds_per_movie
        )
        remaining_movies = plan.movie_count - completed_movies
        status = {
            "updated_at": utc_now(),
            "current_dataset": fov.dataset_name,
            "current_fov": fov.relative_fov,
            "completed_fovs": index,
            "total_fovs": len(plan.fovs),
            "completed_movies": completed_movies,
            "total_movies": plan.movie_count,
            "remaining_movies": remaining_movies,
            "seconds_per_movie": seconds_per_movie,
            "eta_seconds": remaining_movies * seconds_per_movie,
            "last_result": dict(result),
        }
        _write_json(fov.result_dir / "_processing" / "run_status.json", status)
        _append_log(fov.result_dir, f"FOV_DONE {json.dumps(status, default=str)}")
        result["batch_status"] = status
        _emit(progress_callback, result)
    return results


def process_fov(fov: FOVPlan, defaults: dict) -> dict:
    final_paths = (fov.stack_output, fov.mip_output)
    existing = [path for path in final_paths if path.exists()]
    if existing:
        if len(existing) != 2:
            raise FileExistsError(f"Partial final output pair exists: {existing}")
        shapes = validate_output_pair(*final_paths)
        return _result_record(fov, "skipped_existing", 0.0, shapes=shapes)

    fov.result_dir.mkdir(parents=True, exist_ok=True)
    partial_stack = _partial_path(fov.stack_output)
    partial_mip = _partial_path(fov.mip_output)
    for path in (partial_stack, partial_mip):
        if path.exists():
            path.unlink()

    params = SACDParams(
        pixel_nm=fov.pixel_nm,
        wavelength_nm=float(defaults["wavelength_nm"]),
        na=fov.na,
        mag=int(defaults["mag"]),
        iter1=int(defaults["iter1"]),
        iter2=int(defaults["iter2"]),
        ac_order=int(defaults["ac_order"]),
        subfactor=float(defaults["subfactor"]),
        frames_per_sacd=None,
        ifbackground=bool(defaults.get("ifbackground", False)),
        backgroundfactor=float(defaults.get("backgroundfactor", 2.0)),
        ifregistration=bool(defaults.get("ifregistration", False)),
        ifsparsedecon=bool(defaults.get("ifsparsedecon", False)),
        fidelity=float(defaults.get("fidelity", 100.0)),
        tcontinuity=float(defaults.get("tcontinuity", 0.1)),
        sparsity=float(defaults.get("sparsity", 1.0)),
        sparse_iterations=int(defaults.get("sparse_iterations", 100)),
    )
    started = perf_counter()
    sacd_times: list[np.ndarray] = []
    input_shape = None
    for time_index in fov.group.time_indices:
        sacd_z: list[np.ndarray] = []
        for z_index in fov.group.z_indices:
            raw = read_tiff_stack(fov.group.files[(time_index, z_index)])
            input_shape = tuple(raw.shape)
            sacd = reconstruct(raw, params)
            if sacd.ndim != 2:
                raise ValueError(f"Expected a 2D reconstruction, got {sacd.shape}")
            sacd_z.append(sacd)
        sacd_times.append(np.stack(sacd_z, axis=0))

    stack = np.stack(sacd_times, axis=0).astype(np.float32, copy=False)
    mip = np.max(stack, axis=1).astype(np.float32, copy=False)
    pixel_um = fov.pixel_nm / 1000.0 / params.mag
    write_timelapse_tiff(
        partial_stack,
        stack,
        axes="TZYX",
        pixel_size_um=pixel_um,
        z_spacing_um=fov.z_spacing_um,
        time_interval_s=fov.time_interval_s,
    )
    write_timelapse_tiff(
        partial_mip,
        mip,
        axes="TYX",
        pixel_size_um=pixel_um,
        time_interval_s=fov.time_interval_s,
    )
    shapes = validate_output_pair(partial_stack, partial_mip)
    partial_stack.replace(fov.stack_output)
    partial_mip.replace(fov.mip_output)
    return _result_record(
        fov,
        "written",
        perf_counter() - started,
        shapes=shapes,
        input_shape=input_shape,
    )


def validate_output_pair(stack_path: str | Path, mip_path: str | Path) -> dict:
    stack_path, mip_path = Path(stack_path), Path(mip_path)
    with tifffile.TiffFile(stack_path) as tif:
        stack_shape, stack_axes, stack_dtype = tif.series[0].shape, tif.series[0].axes, tif.series[0].dtype
    with tifffile.TiffFile(mip_path) as tif:
        mip_shape, mip_axes, mip_dtype = tif.series[0].shape, tif.series[0].axes, tif.series[0].dtype
    if stack_axes != "TZYX" or mip_axes != "TYX":
        raise ValueError(f"Unexpected axes: {stack_axes}, {mip_axes}")
    if stack_dtype != np.dtype("float32") or mip_dtype != np.dtype("float32"):
        raise ValueError(f"Unexpected dtypes: {stack_dtype}, {mip_dtype}")
    stack = tifffile.imread(stack_path)
    mip = tifffile.imread(mip_path)
    if not np.array_equal(mip, np.max(stack, axis=1)):
        raise ValueError("Saved MIP does not equal max(saved TZYX, axis=1)")
    return {"stack_shape": list(stack_shape), "mip_shape": list(mip_shape)}


def prepare_legacy_provenance(config: dict, config_path: str | Path) -> None:
    plan = build_batch_plan(config)
    _prepare_provenance(config, Path(config_path), plan, Path(config.get("repo_root", Path.cwd())))


def _prepare_provenance(config: dict, config_path: Path, plan: BatchPlan, repo_root: Path) -> None:
    result_dirs = {fov.result_dir for fov in plan.fovs}
    result_dirs.update(Path(item["result_dir"]) for item in plan.legacy_datasets)
    for result_dir in result_dirs:
        processing = result_dir / "_processing"
        processing.mkdir(parents=True, exist_ok=True)
        _write_json(processing / "config.json", config)
        manifest_path = processing / "manifest.json"
        if not manifest_path.exists():
            mode = next(
                (item.get("mode") for item in plan.legacy_datasets if Path(item["result_dir"]) == result_dir),
                "process",
            )
            existing = sorted(str(path) for path in result_dir.glob("*SACDpy-647-posXY0-*.tif"))
            _write_json(
                manifest_path,
                {
                    "created_at": utc_now(),
                    "mode": mode,
                    "existing_outputs": existing if mode == "legacy_provenance" else [],
                    "exclusions": [x for x in plan.exclusions if x["dataset_name"] == result_dir.name],
                    "fovs": [],
                },
            )
        environment = processing / "environment.txt"
        if not environment.exists():
            environment.write_text(_environment_text())
        snapshot = processing / "source_snapshot"
        if not snapshot.exists():
            snapshot.mkdir()
            for name in (
                "SACDpy_pipeline-batch_timelapse_zstack_SACD.ipynb",
                "pyproject.toml",
            ):
                source = repo_root / name
                if source.exists():
                    shutil.copy2(source, snapshot / name)
            for directory in ("src", "scripts"):
                source = repo_root / directory
                if source.exists():
                    shutil.copytree(source, snapshot / directory)


def _result_record(
    fov: FOVPlan,
    status: str,
    runtime_s: float,
    *,
    shapes: dict | None = None,
    input_shape: Iterable[int] | None = None,
    error: str | None = None,
) -> dict:
    return {
        "updated_at": utc_now(),
        "dataset_name": fov.dataset_name,
        "relative_fov": fov.relative_fov,
        "source_folder": str(fov.fov_folder),
        "movie_count": fov.movie_count,
        "time_indices": list(fov.group.time_indices),
        "z_indices": list(fov.group.z_indices),
        "input_shape": list(input_shape) if input_shape else None,
        **(shapes or {}),
        "pixel_nm": fov.pixel_nm,
        "na": fov.na,
        "z_spacing_um": fov.z_spacing_um,
        "time_interval_s": fov.time_interval_s,
        "stack_output": str(fov.stack_output),
        "mip_output": str(fov.mip_output),
        "runtime_s": runtime_s,
        "status": status,
        "error": error,
    }


def _load_manifest(result_dir: Path) -> dict:
    path = result_dir / "_processing" / "manifest.json"
    with path.open() as handle:
        return json.load(handle)


def _upsert_manifest(result_dir: Path, result: dict) -> None:
    manifest = _load_manifest(result_dir)
    manifest["fovs"] = [
        item for item in manifest["fovs"] if item.get("relative_fov") != result["relative_fov"]
    ]
    manifest["fovs"].append({k: v for k, v in result.items() if k != "batch_status"})
    manifest["updated_at"] = utc_now()
    _write_json(result_dir / "_processing" / "manifest.json", manifest)


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    temporary.replace(path)


def _append_log(result_dir: Path, line: str) -> None:
    with (result_dir / "_processing" / "run.log").open("a") as handle:
        handle.write(f"{utc_now()} {line}\n")


def _partial_path(path: Path) -> Path:
    return path.with_name(f"{path.stem}.partial{path.suffix}")


def _exclusion(raw_root: Path, relative: str, reason: str) -> dict[str, str]:
    return {"dataset_name": raw_root.name, "relative_fov": relative, "reason": reason}


def _environment_text() -> str:
    packages = ["numpy", "scipy", "scikit-image", "tifffile", "rich", "PyWavelets"]
    lines = [f"created_at={utc_now()}", f"python={sys.version}", f"platform={platform.platform()}"]
    for package in packages:
        try:
            package_version = version(package)
        except PackageNotFoundError:
            package_version = "not installed"
        lines.append(f"{package}={package_version}")
    return "\n".join(lines) + "\n"


def _emit(callback: Callable[[dict], None] | None, result: dict) -> None:
    if callback:
        callback(result)
    else:
        print("FOV_DONE " + json.dumps(result, default=str), flush=True)
