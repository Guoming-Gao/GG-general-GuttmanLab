from __future__ import annotations

import json
import platform
import re
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from statistics import median
from time import perf_counter
from typing import Callable, Iterable

import numpy as np
import tifffile

from .params import SACDParams
from .reconstruction import reconstruct


_ONI_ZSTACK_RE = re.compile(
    r"^(?P<prefix>.+)_posXY(?P<pos_xy>\d+)_channels_t(?P<time_index>\d+)_posZ(?P<z_index>\d+)\.tiff?$",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class ChannelSpec:
    name: str
    wavelength_nm: float
    camera_half: str
    metadata_wavelengths_nm: tuple[float, ...]


@dataclass(frozen=True)
class MovieInfo:
    path: Path
    prefix: str
    pos_xy: int
    time_index: int
    z_index: int
    shape: tuple[int, int, int]
    dtype: str
    metadata: dict
    frame_ranges: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class FOVPlan:
    raw_root: Path
    output_root: Path
    fov_folder: Path
    relative_fov: str
    prefix: str
    pos_xy: int
    movies: tuple[MovieInfo, ...]
    channels: tuple[ChannelSpec, ...]
    stack_output: Path
    mip_output: Path
    pixel_nm: float
    na: float
    z_spacing_um: float | None

    @property
    def z_indices(self) -> tuple[int, ...]:
        return tuple(movie.z_index for movie in self.movies)

    @property
    def reconstruction_count(self) -> int:
        return len(self.movies) * len(self.channels)


@dataclass(frozen=True)
class BatchPlan:
    raw_root: Path
    output_root: Path
    fovs: tuple[FOVPlan, ...]
    exclusions: tuple[dict[str, str], ...]

    @property
    def reconstruction_count(self) -> int:
        return sum(fov.reconstruction_count for fov in self.fovs)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_config(path: str | Path) -> dict:
    with Path(path).open() as handle:
        config = json.load(handle)
    if config.get("version") != 1:
        raise ValueError("Multicolor z-stack config version must be 1")
    if not config.get("raw_root") or not config.get("output_root"):
        raise ValueError("Config must define raw_root and output_root")
    _channel_specs(config)
    return config


def _channel_specs(config: dict) -> tuple[ChannelSpec, ...]:
    entries = config.get("processing", {}).get("channels", ())
    if not entries:
        raise ValueError("processing.channels must contain at least one channel")
    channels: list[ChannelSpec] = []
    names: set[str] = set()
    for entry in entries:
        name = str(entry["name"])
        if name in names:
            raise ValueError(f"Duplicate channel name: {name}")
        names.add(name)
        camera_half = str(entry["camera_half"]).lower()
        if camera_half not in {"left", "right"}:
            raise ValueError(f"camera_half for channel {name} must be left or right")
        wavelength = float(entry["wavelength_nm"])
        accepted = tuple(float(value) for value in entry.get("metadata_wavelengths_nm", [wavelength]))
        channels.append(ChannelSpec(name, wavelength, camera_half, accepted))
    return tuple(channels)


def parse_movie_filename(path: str | Path) -> tuple[str, int, int, int]:
    path = Path(path)
    match = _ONI_ZSTACK_RE.match(path.name)
    if match is None:
        raise ValueError(f"Filename does not match ONI multicolor z-stack grammar: {path.name}")
    return (
        match.group("prefix"),
        int(match.group("pos_xy")),
        int(match.group("time_index")),
        int(match.group("z_index")),
    )


def read_movie_info(path: str | Path, channels: Iterable[ChannelSpec]) -> MovieInfo:
    path = Path(path)
    prefix, pos_xy, time_index, z_index = parse_movie_filename(path)
    with tifffile.TiffFile(path) as tif:
        page_count = len(tif.pages)
        if not page_count:
            raise ValueError(f"TIFF has no pages: {path}")
        height, width = tif.pages[0].shape
        dtype = str(tif.pages[0].dtype)
        description = tif.pages[0].description
    try:
        metadata = json.loads(description) if description else {}
    except json.JSONDecodeError as exc:
        raise ValueError(f"First TIFF page does not contain valid ONI JSON metadata: {path}") from exc
    if not isinstance(metadata, dict):
        raise ValueError(f"ONI metadata must be a JSON object: {path}")

    channel_tuple = tuple(channels)
    wavelengths = metadata.get("LaserWavelength_nm")
    steps = metadata.get("laserProgram", {}).get("steps")
    if not isinstance(wavelengths, list) or not isinstance(steps, list):
        raise ValueError(f"Missing LaserWavelength_nm or laserProgram.steps metadata: {path}")
    if len(wavelengths) != len(channel_tuple) or len(steps) != len(channel_tuple):
        raise ValueError(
            f"Expected {len(channel_tuple)} metadata channel(s), got "
            f"{len(wavelengths)} wavelength(s) and {len(steps)} laser step(s): {path}"
        )

    frame_ranges: list[tuple[int, int]] = []
    start = 0
    for index, (spec, observed, step) in enumerate(zip(channel_tuple, wavelengths, steps)):
        if not any(np.isclose(float(observed), allowed) for allowed in spec.metadata_wavelengths_nm):
            raise ValueError(
                f"Channel {index} ({spec.name}) expected metadata wavelength "
                f"{spec.metadata_wavelengths_nm}, got {observed}: {path}"
            )
        repeats = step.get("nRepeats") if isinstance(step, dict) else None
        if not isinstance(repeats, int) or repeats <= 0:
            raise ValueError(f"Invalid nRepeats for channel {spec.name}: {path}")
        frame_ranges.append((start, start + repeats))
        start += repeats
    if start != page_count:
        raise ValueError(
            f"Laser-program repeats total {start}, but TIFF contains {page_count} page(s): {path}"
        )
    if width % 2:
        raise ValueError(f"Camera-split TIFF width must be even, got {width}: {path}")

    return MovieInfo(
        path=path,
        prefix=prefix,
        pos_xy=pos_xy,
        time_index=time_index,
        z_index=z_index,
        shape=(page_count, height, width),
        dtype=dtype,
        metadata=metadata,
        frame_ranges=tuple(frame_ranges),
    )


def build_batch_plan(config: dict) -> BatchPlan:
    processing = config["processing"]
    channels = _channel_specs(config)
    raw_root = Path(config["raw_root"])
    output_root = Path(config["output_root"])
    if not raw_root.is_dir():
        raise FileNotFoundError(f"Raw dataset root does not exist: {raw_root}")

    position_folder = str(processing.get("position_folder", "pos_0"))
    glob_pattern = str(processing.get("glob_pattern", "*.tif"))
    selected = set(config.get("selected_fov_folders", ()))
    excluded_config = config.get("exclude_fovs", {})
    aliases = config.get("output_prefix_aliases", {})
    position_dirs = sorted(path for path in raw_root.rglob(position_folder) if path.is_dir())
    fovs: list[FOVPlan] = []
    exclusions: list[dict[str, str]] = []
    planned_paths: set[Path] = set()

    for position_dir in position_dirs:
        fov_folder = position_dir.parent
        relative_fov = fov_folder.relative_to(raw_root).as_posix()
        if selected and relative_fov not in selected:
            exclusions.append(_exclusion(relative_fov, "not selected by explicit FOV list"))
            continue
        if relative_fov in excluded_config:
            exclusions.append(_exclusion(relative_fov, str(excluded_config[relative_fov])))
            continue

        candidates = sorted(
            path for path in position_dir.glob(glob_pattern)
            if path.is_file() and "sacdpy" not in path.name.lower()
        )
        parsed: list[MovieInfo] = []
        for candidate in candidates:
            try:
                parsed.append(read_movie_info(candidate, channels))
            except ValueError as exc:
                if _ONI_ZSTACK_RE.match(candidate.name):
                    raise
                continue
        if not parsed:
            exclusions.append(_exclusion(relative_fov, "no matching ONI multicolor z-stack TIFFs"))
            continue

        keys = [(item.prefix, item.pos_xy, item.time_index) for item in parsed]
        if len(set(keys)) != 1:
            raise ValueError(f"Expected one prefix/position/time group in {position_dir}, found {sorted(set(keys))}")
        prefix, pos_xy, time_index = keys[0]
        if time_index != 0:
            raise ValueError(f"Multicolor z-stack pipeline expects t0 only, got t{time_index}: {position_dir}")
        z_indices = [item.z_index for item in parsed]
        if len(set(z_indices)) != len(z_indices):
            raise ValueError(f"Duplicate z-plane index in {position_dir}: {z_indices}")
        expected_z = list(range(min(z_indices), max(z_indices) + 1))
        if sorted(z_indices) != expected_z:
            missing = sorted(set(expected_z) - set(z_indices))
            raise ValueError(f"Missing z-plane index/indices {missing} in {position_dir}")
        movies = tuple(sorted(parsed, key=lambda item: item.z_index))
        _validate_consistent_movies(movies)

        output_prefix = str(aliases.get(relative_fov, prefix))
        stack_output, mip_output = output_paths(output_root, output_prefix, pos_xy)
        for path in (stack_output, mip_output):
            if path in planned_paths:
                raise ValueError(f"Output path collision: {path}")
            planned_paths.add(path)

        pixel_nm = processing.get("pixel_nm")
        if pixel_nm is None:
            value = movies[0].metadata.get("PixelSize_um")
            pixel_nm = float(value) * 1000.0 if value is not None else processing.get("fallback_pixel_nm")
        na = processing.get("na")
        if na is None:
            na = movies[0].metadata.get("Objective_NA", processing.get("fallback_na"))
        if pixel_nm is None or na is None:
            raise ValueError(f"Could not resolve pixel size or NA for {position_dir}")

        fovs.append(
            FOVPlan(
                raw_root=raw_root,
                output_root=output_root,
                fov_folder=fov_folder,
                relative_fov=relative_fov,
                prefix=output_prefix,
                pos_xy=pos_xy,
                movies=movies,
                channels=channels,
                stack_output=stack_output,
                mip_output=mip_output,
                pixel_nm=float(pixel_nm),
                na=float(na),
                z_spacing_um=_infer_z_spacing_um(movies),
            )
        )

    return BatchPlan(raw_root, output_root, tuple(fovs), tuple(exclusions))


def _validate_consistent_movies(movies: tuple[MovieInfo, ...]) -> None:
    first = movies[0]
    for movie in movies[1:]:
        if movie.shape != first.shape or movie.dtype != first.dtype:
            raise ValueError(
                f"Inconsistent TIFF shape/dtype: {first.path.name}={first.shape}/{first.dtype}, "
                f"{movie.path.name}={movie.shape}/{movie.dtype}"
            )
        if movie.frame_ranges != first.frame_ranges:
            raise ValueError(f"Inconsistent laser-program frame ranges: {movie.path}")
        for key in ("LaserWavelength_nm", "PixelSize_um", "Objective_NA", "ROI"):
            if movie.metadata.get(key) != first.metadata.get(key):
                raise ValueError(f"Inconsistent {key} metadata: {movie.path}")


def _infer_z_spacing_um(movies: Iterable[MovieInfo]) -> float | None:
    positions: list[float] = []
    for movie in movies:
        stage_pos = movie.metadata.get("StagePos_um")
        if isinstance(stage_pos, list) and len(stage_pos) >= 3:
            positions.append(float(stage_pos[2]))
    spacings = [abs(right - left) for left, right in zip(positions, positions[1:])]
    return float(median(spacings)) if spacings else None


def output_paths(output_root: str | Path, prefix: str, pos_xy: int) -> tuple[Path, Path]:
    root = Path(output_root)
    stem = f"{prefix}-SACDpy-4color-posXY{pos_xy}"
    return root / f"{stem}-CZYX.ome.tif", root / f"{stem}-MIP-CYX.ome.tif"


def preflight_summary(plan: BatchPlan) -> dict:
    return {
        "raw_root": str(plan.raw_root),
        "output_root": str(plan.output_root),
        "fovs": len(plan.fovs),
        "movies": sum(len(fov.movies) for fov in plan.fovs),
        "channels": [channel.name for channel in plan.fovs[0].channels] if plan.fovs else [],
        "reconstructions": plan.reconstruction_count,
        "outputs": len(plan.fovs) * 2,
        "exclusions": len(plan.exclusions),
    }


def select_channel_frames(raw: np.ndarray, frame_range: tuple[int, int], camera_half: str) -> np.ndarray:
    arr = np.asarray(raw)
    if arr.ndim != 3:
        raise ValueError(f"Expected a TYX movie, got shape {arr.shape}")
    if arr.shape[2] % 2:
        raise ValueError(f"Camera-split movie width must be even, got {arr.shape[2]}")
    start, end = frame_range
    if start < 0 or end <= start or end > arr.shape[0]:
        raise ValueError(f"Invalid zero-based frame range [{start}, {end}) for {arr.shape[0]} frames")
    halfwidth = arr.shape[2] // 2
    if camera_half == "left":
        selected = arr[start:end, :, :halfwidth]
    elif camera_half == "right":
        selected = arr[start:end, :, halfwidth:]
    else:
        raise ValueError(f"camera_half must be left or right, got {camera_half!r}")
    return np.ascontiguousarray(selected)


def _params_for_channel(fov: FOVPlan, channel: ChannelSpec, processing: dict) -> SACDParams:
    return SACDParams(
        pixel_nm=fov.pixel_nm,
        wavelength_nm=channel.wavelength_nm,
        na=fov.na,
        mag=int(processing.get("mag", 2)),
        iter1=int(processing.get("iter1", 7)),
        iter2=int(processing.get("iter2", 8)),
        ac_order=int(processing.get("ac_order", 2)),
        subfactor=float(processing.get("subfactor", 0.8)),
        frames_per_sacd=None,
        ifbackground=bool(processing.get("ifbackground", False)),
        backgroundfactor=float(processing.get("backgroundfactor", 2.0)),
        ifregistration=bool(processing.get("ifregistration", False)),
        ifsparsedecon=bool(processing.get("ifsparsedecon", False)),
        fidelity=float(processing.get("fidelity", 100.0)),
        tcontinuity=float(processing.get("tcontinuity", 0.1)),
        sparsity=float(processing.get("sparsity", 1.0)),
        sparse_iterations=int(processing.get("sparse_iterations", 100)),
    )


def process_fov(
    fov: FOVPlan,
    processing: dict,
    *,
    progress_callback: Callable[[dict], None] | None = None,
) -> dict:
    existing = [path for path in (fov.stack_output, fov.mip_output) if path.exists()]
    if existing:
        if len(existing) != 2:
            raise FileExistsError(f"Partial final output pair exists: {existing}")
        shapes = validate_output_pair(fov.stack_output, fov.mip_output, fov.channels)
        return _result_record(fov, "skipped_existing", 0.0, shapes=shapes)

    fov.output_root.mkdir(parents=True, exist_ok=True)
    partial_stack = _partial_path(fov.stack_output)
    partial_mip = _partial_path(fov.mip_output)
    for partial in (partial_stack, partial_mip):
        if partial.exists():
            partial.unlink()

    started = perf_counter()
    channel_planes: list[list[np.ndarray]] = [[] for _ in fov.channels]
    input_shape: tuple[int, ...] | None = None
    completed_reconstructions = 0
    for movie in fov.movies:
        raw = tifffile.imread(movie.path)
        input_shape = tuple(raw.shape)
        for channel_index, channel in enumerate(fov.channels):
            reconstruction_started = perf_counter()
            selected = select_channel_frames(raw, movie.frame_ranges[channel_index], channel.camera_half)
            sacd = reconstruct(selected, _params_for_channel(fov, channel, processing))
            if sacd.ndim != 2:
                raise ValueError(f"Expected one 2D SACD image, got {sacd.shape} for {movie.path}")
            channel_planes[channel_index].append(np.asarray(sacd, dtype=np.float32))
            completed_reconstructions += 1
            if progress_callback is not None:
                progress_callback(
                    {
                        "event": "reconstruction_done",
                        "relative_fov": fov.relative_fov,
                        "z_index": movie.z_index,
                        "channel": channel.name,
                        "completed_in_fov": completed_reconstructions,
                        "total_in_fov": fov.reconstruction_count,
                        "runtime_s": perf_counter() - reconstruction_started,
                    }
                )

    stack = np.stack([np.stack(planes, axis=0) for planes in channel_planes], axis=0)
    mip = np.max(stack, axis=1).astype(np.float32, copy=False)
    pixel_um = fov.pixel_nm / 1000.0 / int(processing.get("mag", 2))
    write_ome_tiff(
        partial_stack,
        stack,
        axes="CZYX",
        channels=fov.channels,
        pixel_size_um=pixel_um,
        z_spacing_um=fov.z_spacing_um,
    )
    write_ome_tiff(
        partial_mip,
        mip,
        axes="CYX",
        channels=fov.channels,
        pixel_size_um=pixel_um,
    )
    shapes = validate_output_pair(partial_stack, partial_mip, fov.channels)
    partial_stack.replace(fov.stack_output)
    partial_mip.replace(fov.mip_output)
    return _result_record(
        fov,
        "written",
        perf_counter() - started,
        shapes=shapes,
        input_shape=input_shape,
    )


def write_ome_tiff(
    path: str | Path,
    image: np.ndarray,
    *,
    axes: str,
    channels: Iterable[ChannelSpec],
    pixel_size_um: float,
    z_spacing_um: float | None = None,
) -> None:
    arr = np.asarray(image, dtype=np.float32)
    channel_tuple = tuple(channels)
    if arr.ndim != len(axes):
        raise ValueError(f"Image ndim {arr.ndim} does not match axes {axes!r}")
    if "C" not in axes or arr.shape[axes.index("C")] != len(channel_tuple):
        raise ValueError("OME image channel axis does not match channel metadata")
    metadata: dict[str, object] = {
        "axes": axes,
        "Channel": {
            "Name": [channel.name for channel in channel_tuple],
            "ExcitationWavelength": [channel.wavelength_nm for channel in channel_tuple],
            "ExcitationWavelengthUnit": ["nm"] * len(channel_tuple),
        },
        "PhysicalSizeX": float(pixel_size_um),
        "PhysicalSizeXUnit": "µm",
        "PhysicalSizeY": float(pixel_size_um),
        "PhysicalSizeYUnit": "µm",
    }
    if "Z" in axes and z_spacing_um is not None:
        metadata["PhysicalSizeZ"] = float(z_spacing_um)
        metadata["PhysicalSizeZUnit"] = "µm"
    tifffile.imwrite(
        path,
        arr,
        ome=True,
        bigtiff=arr.nbytes >= 4_000_000_000,
        photometric="minisblack",
        metadata=metadata,
    )


def validate_output_pair(
    stack_path: str | Path,
    mip_path: str | Path,
    channels: Iterable[ChannelSpec] | None = None,
) -> dict:
    stack_path, mip_path = Path(stack_path), Path(mip_path)
    with tifffile.TiffFile(stack_path) as tif:
        stack_shape = tuple(tif.series[0].shape)
        stack_axes = tif.series[0].axes
        stack_dtype = tif.series[0].dtype
        stack_xml = tif.ome_metadata or ""
    with tifffile.TiffFile(mip_path) as tif:
        mip_shape = tuple(tif.series[0].shape)
        mip_axes = tif.series[0].axes
        mip_dtype = tif.series[0].dtype
        mip_xml = tif.ome_metadata or ""
    if stack_axes != "CZYX" or mip_axes != "CYX":
        raise ValueError(f"Unexpected output axes: {stack_axes}, {mip_axes}")
    if stack_dtype != np.dtype("float32") or mip_dtype != np.dtype("float32"):
        raise ValueError(f"Unexpected output dtypes: {stack_dtype}, {mip_dtype}")
    channel_tuple = tuple(channels or ())
    for channel in channel_tuple:
        marker = f'Name="{channel.name}"'
        if marker not in stack_xml or marker not in mip_xml:
            raise ValueError(f"Missing OME channel name {channel.name!r}")
    stack = tifffile.imread(stack_path)
    mip = tifffile.imread(mip_path)
    if not np.array_equal(mip, np.max(stack, axis=1)):
        raise ValueError("Saved CYX MIP does not equal max(saved CZYX, axis=Z)")
    return {"stack_shape": list(stack_shape), "mip_shape": list(mip_shape)}


def run_batch(
    config: dict,
    config_path: str | Path,
    *,
    max_new_fovs: int | None = None,
    progress_callback: Callable[[dict], None] | None = None,
) -> list[dict]:
    plan = build_batch_plan(config)
    _prepare_provenance(plan, config, Path(config_path))
    processing = config["processing"]
    manifest = _load_manifest(plan.output_root)
    completed_reconstructions = 0
    elapsed_processing = 0.0
    newly_written = 0
    results: list[dict] = []

    for index, fov in enumerate(plan.fovs, start=1):
        previous = next(
            (item for item in manifest["fovs"] if item.get("relative_fov") == fov.relative_fov),
            None,
        )
        if previous and previous.get("status") in {"written", "skipped_existing"}:
            try:
                validate_output_pair(fov.stack_output, fov.mip_output, fov.channels)
            except (FileNotFoundError, ValueError):
                previous = None
        if previous is not None:
            result = {**previous, "status": "resumed_manifest"}
            results.append(result)
            completed_reconstructions += fov.reconstruction_count
            _emit(progress_callback, result)
            continue
        if max_new_fovs is not None and newly_written >= max_new_fovs:
            break

        try:
            step_callback = progress_callback
            if step_callback is None:
                step_callback = lambda event: print(
                    "SACD_STEP " + json.dumps(event, default=str), flush=True
                )
            result = process_fov(fov, processing, progress_callback=step_callback)
            newly_written += result["status"] == "written"
        except Exception as exc:
            result = _result_record(fov, "failed", 0.0, error=repr(exc))
        results.append(result)
        manifest = _upsert_manifest(plan.output_root, manifest, result, plan.exclusions)

        completed_reconstructions += fov.reconstruction_count
        elapsed_processing += float(result.get("runtime_s", 0.0))
        measured = sum(
            item["reconstruction_count"] for item in results if item.get("status") == "written"
        )
        seconds_per_reconstruction = (
            elapsed_processing / measured
            if measured
            else float(config.get("initial_seconds_per_reconstruction", 10.0))
        )
        remaining = plan.reconstruction_count - completed_reconstructions
        status = {
            "updated_at": utc_now(),
            "current_fov": fov.relative_fov,
            "completed_fovs": index,
            "total_fovs": len(plan.fovs),
            "completed_reconstructions": completed_reconstructions,
            "total_reconstructions": plan.reconstruction_count,
            "remaining_reconstructions": remaining,
            "seconds_per_reconstruction": seconds_per_reconstruction,
            "eta_seconds": remaining * seconds_per_reconstruction,
            "last_result": dict(result),
        }
        _write_json(plan.output_root / "_processing" / "run_status.json", status)
        _append_log(plan.output_root, f"FOV_DONE {json.dumps(status, default=str)}")
        result["batch_status"] = status
        _emit(progress_callback, result)
    return results


def _prepare_provenance(plan: BatchPlan, config: dict, config_path: Path) -> None:
    processing_dir = plan.output_root / "_processing"
    processing_dir.mkdir(parents=True, exist_ok=True)
    _write_json(processing_dir / "config.json", config)
    manifest_path = processing_dir / "manifest.json"
    if not manifest_path.exists():
        _write_json(
            manifest_path,
            {
                "created_at": utc_now(),
                "config_source": str(config_path),
                "exclusions": list(plan.exclusions),
                "fovs": [],
            },
        )
    environment_path = processing_dir / "environment.txt"
    if not environment_path.exists():
        environment_path.write_text(_environment_text())


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
        "relative_fov": fov.relative_fov,
        "source_folder": str(fov.fov_folder),
        "z_indices": list(fov.z_indices),
        "channels": [asdict(channel) for channel in fov.channels],
        "frame_ranges_zero_based_end_exclusive": [list(value) for value in fov.movies[0].frame_ranges],
        "reconstruction_count": fov.reconstruction_count,
        "input_shape": list(input_shape) if input_shape else None,
        **(shapes or {}),
        "pixel_nm": fov.pixel_nm,
        "na": fov.na,
        "z_spacing_um": fov.z_spacing_um,
        "stack_output": str(fov.stack_output),
        "mip_output": str(fov.mip_output),
        "runtime_s": runtime_s,
        "status": status,
        "error": error,
    }


def _partial_path(path: str | Path) -> Path:
    path = Path(path)
    return path.with_name(f"{path.stem}.partial{path.suffix}")


def _exclusion(relative_fov: str, reason: str) -> dict[str, str]:
    return {"relative_fov": relative_fov, "reason": reason}


def _load_manifest(output_root: Path) -> dict:
    with (output_root / "_processing" / "manifest.json").open() as handle:
        return json.load(handle)


def _upsert_manifest(output_root: Path, manifest: dict, result: dict, exclusions: tuple[dict, ...]) -> dict:
    updated = dict(manifest)
    updated["fovs"] = [
        item for item in manifest.get("fovs", ()) if item.get("relative_fov") != result["relative_fov"]
    ]
    updated["fovs"].append({key: value for key, value in result.items() if key != "batch_status"})
    updated["exclusions"] = list(exclusions)
    updated["updated_at"] = utc_now()
    _write_json(output_root / "_processing" / "manifest.json", updated)
    return updated


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    temporary.replace(path)


def _append_log(output_root: Path, line: str) -> None:
    with (output_root / "_processing" / "run.log").open("a") as handle:
        handle.write(f"{utc_now()} {line}\n")


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
