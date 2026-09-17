from __future__ import annotations

import json
import os
import platform
import re
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from statistics import median
from time import perf_counter
from typing import Callable, Iterable

import numpy as np
import tifffile

from .params import SACDParams
from .execution import execute_tasks, initialize_worker, pooled_batch, worker_count
from .progress import FOVEvents
from .reconstruction import apply_intensity_transform, reconstruct


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
class ChannelOutputPaths:
    channel: str
    stack: Path
    mip: Path
    raw_max_mip: Path

    @property
    def all(self) -> tuple[Path, Path, Path]:
        return self.stack, self.mip, self.raw_max_mip


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
    outputs: tuple[ChannelOutputPaths, ...]
    pixel_nm: float
    na: float
    z_spacing_um: float | None

    @property
    def z_indices(self) -> tuple[int, ...]:
        return tuple(movie.z_index for movie in self.movies)

    @property
    def reconstruction_count(self) -> int:
        return len(self.movies) * len(self.channels)

    @property
    def all_outputs(self) -> tuple[Path, ...]:
        return tuple(path for output in self.outputs for path in output.all)

    def output_for(self, channel: str) -> ChannelOutputPaths:
        return next(output for output in self.outputs if output.channel == channel)


@dataclass(frozen=True)
class BatchPlan:
    raw_root: Path
    output_root: Path
    fovs: tuple[FOVPlan, ...]
    exclusions: tuple[dict[str, str], ...]

    @property
    def reconstruction_count(self) -> int:
        return sum(fov.reconstruction_count for fov in self.fovs)


@dataclass(frozen=True)
class ReconstructionTask:
    movie_path: Path
    frame_range: tuple[int, int]
    camera_half: str
    params: SACDParams
    z_position: int
    z_index: int
    channel_index: int
    channel_name: str


@dataclass(frozen=True)
class ReconstructionResult:
    z_position: int
    z_index: int
    channel_index: int
    channel_name: str
    input_shape: tuple[int, ...]
    sacd: np.ndarray
    frame_max: np.ndarray
    runtime_s: float


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


def read_movie_info(
    path: str | Path,
    channels: Iterable[ChannelSpec],
    *,
    frame_mode: str = "auto",
) -> MovieInfo:
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
    if frame_mode not in {"auto", "sequential", "simultaneous"}:
        raise ValueError("processing.frame_mode must be auto, sequential, or simultaneous")
    resolved_mode = frame_mode
    if resolved_mode == "auto":
        resolved_mode = "sequential" if steps else "simultaneous"

    frame_ranges: list[tuple[int, int]] = []
    if resolved_mode == "sequential":
        if len(steps) != len(channel_tuple):
            raise ValueError(
                f"Sequential mode expected {len(channel_tuple)} laser step(s), got "
                f"{len(steps)}: {path}"
            )
        start = 0
        for index, (spec, step) in enumerate(zip(channel_tuple, steps, strict=True)):
            observed = wavelengths[index]
            if not any(
                np.isclose(float(observed), allowed)
                for allowed in spec.metadata_wavelengths_nm
            ):
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
    else:
        if steps:
            raise ValueError(f"Simultaneous mode requires an empty laser program: {path}")
        active = metadata.get("LaserActive")
        if not isinstance(active, list) or len(active) != len(wavelengths):
            raise ValueError(f"Simultaneous mode requires LaserActive metadata: {path}")
        powers = metadata.get("LaserPowerPercent")
        if not isinstance(powers, list) or len(powers) != len(wavelengths):
            raise ValueError(f"Simultaneous mode requires LaserPowerPercent metadata: {path}")
        for spec in channel_tuple:
            matching_indices = [
                index
                for index, observed in enumerate(wavelengths)
                if any(
                    np.isclose(float(observed), allowed)
                    for allowed in spec.metadata_wavelengths_nm
                )
            ]
            if not matching_indices:
                raise ValueError(
                    f"Channel {spec.name} wavelength {spec.metadata_wavelengths_nm} is absent "
                    f"from LaserWavelength_nm={wavelengths}: {path}"
                )
            if not any(
                bool(active[index]) or float(powers[index]) > 0
                for index in matching_indices
            ):
                raise ValueError(
                    f"Configured channel {spec.name} has neither an active flag nor positive "
                    f"laser power: {path}"
                )
            frame_ranges.append((0, page_count))
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
    frame_mode = str(processing.get("frame_mode", "auto"))
    prefix_source = str(processing.get("output_prefix_source", "fov_folder"))
    if prefix_source not in {"fov_folder", "movie_prefix"}:
        raise ValueError("processing.output_prefix_source must be fov_folder or movie_prefix")
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
                parsed.append(read_movie_info(candidate, channels, frame_mode=frame_mode))
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

        default_prefix = fov_folder.name if prefix_source == "fov_folder" else prefix
        output_prefix = str(aliases.get(relative_fov, default_prefix))
        outputs = output_paths(output_root, output_prefix, channels)
        for path in (path for output in outputs for path in output.all):
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
                outputs=outputs,
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


def output_paths(
    output_root: str | Path,
    prefix: str,
    channels: Iterable[ChannelSpec],
) -> tuple[ChannelOutputPaths, ...]:
    root = Path(output_root)
    return tuple(
        ChannelOutputPaths(
            channel=channel.name,
            stack=root / f"{prefix}__{channel.name}-SACD-ZYX.tif",
            mip=root / f"{prefix}__{channel.name}-SACD-MIP-YX.tif",
            raw_max_mip=root / f"{prefix}__{channel.name}-raw-max-MIP-YX.tif",
        )
        for channel in channels
    )


def legacy_ome_paths(fov: FOVPlan) -> tuple[Path, Path]:
    stem = f"{fov.movies[0].prefix}-SACDpy-4color-posXY{fov.pos_xy}"
    return (
        fov.output_root / f"{stem}-CZYX.ome.tif",
        fov.output_root / f"{stem}-MIP-CYX.ome.tif",
    )


def preflight_summary(plan: BatchPlan) -> dict:
    return {
        "raw_root": str(plan.raw_root),
        "output_root": str(plan.output_root),
        "fovs": len(plan.fovs),
        "movies": sum(len(fov.movies) for fov in plan.fovs),
        "channels": [channel.name for channel in plan.fovs[0].channels] if plan.fovs else [],
        "reconstructions": plan.reconstruction_count,
        "outputs": len(plan.fovs) * len(plan.fovs[0].channels) * 3 if plan.fovs else 0,
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
        intensity_transform=str(processing.get("intensity_transform", "order_root")),
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


def _initialize_reconstruction_worker() -> None:
    """Keep each SACD worker single-threaded so processes do not oversubscribe the CPU."""
    initialize_worker()


def _reconstruct_task(task: ReconstructionTask) -> ReconstructionResult:
    started = perf_counter()
    raw = tifffile.imread(task.movie_path)
    selected = select_channel_frames(raw, task.frame_range, task.camera_half)
    sacd = reconstruct(selected, task.params)
    if sacd.ndim != 2:
        raise ValueError(f"Expected one 2D SACD image, got {sacd.shape} for {task.movie_path}")
    return ReconstructionResult(
        z_position=task.z_position,
        z_index=task.z_index,
        channel_index=task.channel_index,
        channel_name=task.channel_name,
        input_shape=tuple(raw.shape),
        sacd=np.asarray(sacd, dtype=np.float32),
        frame_max=np.max(selected, axis=0),
        runtime_s=perf_counter() - started,
    )


def _task_sequence(fov: FOVPlan, processing: dict) -> tuple[ReconstructionTask, ...]:
    return tuple(
        ReconstructionTask(
            movie_path=movie.path,
            frame_range=movie.frame_ranges[channel_index],
            camera_half=channel.camera_half,
            params=_params_for_channel(fov, channel, processing),
            z_position=z_position,
            z_index=movie.z_index,
            channel_index=channel_index,
            channel_name=channel.name,
        )
        for z_position, movie in enumerate(fov.movies)
        for channel_index, channel in enumerate(fov.channels)
    )


def _execute_reconstruction_tasks(
    tasks: tuple[ReconstructionTask, ...],
    *,
    executor: ProcessPoolExecutor | None,
    max_workers: int,
) -> Iterable[ReconstructionResult]:
    yield from execute_tasks(
        _reconstruct_task, tasks, executor=executor, max_workers=max_workers,
    )


def process_fov(
    fov: FOVPlan,
    processing: dict,
    *,
    progress_callback: Callable[[dict], None] | None = None,
    executor: ProcessPoolExecutor | None = None,
) -> dict:
    intensity_transform = str(processing.get("intensity_transform", "order_root"))
    existing = [path for path in fov.all_outputs if path.exists()]
    if existing:
        if len(existing) != len(fov.all_outputs):
            raise FileExistsError(f"Partial final multicolor output set exists: {existing}")
        shapes = validate_fov_outputs(
            fov,
            intensity_transform=intensity_transform,
            mag=int(processing.get("mag", 2)),
        )
        return _result_record(
            fov,
            "skipped_existing",
            0.0,
            shapes=shapes,
            intensity_transform=intensity_transform,
        )
    legacy_existing = [path for path in legacy_ome_paths(fov) if path.exists()]
    if legacy_existing:
        raise ValueError(
            f"Legacy raw-cumulant OME outputs exist for {fov.relative_fov}. "
            "Run migrate_legacy_ome_dataset before normal batch processing."
        )

    fov.output_root.mkdir(parents=True, exist_ok=True)
    partial_outputs = tuple(
        ChannelOutputPaths(
            output.channel,
            _partial_path(output.stack),
            _partial_path(output.mip),
            _partial_path(output.raw_max_mip),
        )
        for output in fov.outputs
    )
    for partial in (path for output in partial_outputs for path in output.all):
        if partial.exists():
            partial.unlink()

    started = perf_counter()
    max_workers = worker_count(processing)
    channel_planes: list[list[np.ndarray | None]] = [
        [None for _ in fov.movies] for _ in fov.channels
    ]
    raw_max_mips: list[np.ndarray | None] = [None for _ in fov.channels]
    input_shape: tuple[int, ...] | None = None
    completed_reconstructions = 0
    tasks = _task_sequence(fov, processing)
    for result in _execute_reconstruction_tasks(
        tasks,
        executor=executor,
        max_workers=max_workers,
    ):
        input_shape = result.input_shape
        channel_planes[result.channel_index][result.z_position] = result.sacd
        raw_max_mips[result.channel_index] = (
            result.frame_max
            if raw_max_mips[result.channel_index] is None
            else np.maximum(raw_max_mips[result.channel_index], result.frame_max)
        )
        completed_reconstructions += 1
        if progress_callback is not None:
            progress_callback(
                {
                    "event": "reconstruction_done",
                    "relative_fov": fov.relative_fov,
                    "z_index": result.z_index,
                    "channel": result.channel_name,
                    "completed_in_fov": completed_reconstructions,
                    "total_in_fov": fov.reconstruction_count,
                    "runtime_s": result.runtime_s,
                }
            )

    if any(plane is None for planes in channel_planes for plane in planes):
        raise RuntimeError(f"Incomplete parallel reconstruction result set for {fov.relative_fov}")
    stacks = [
        np.stack(planes, axis=0).astype(np.float32, copy=False)
        for planes in channel_planes
    ]
    raw_references = [_require_uint16(image) for image in raw_max_mips]
    if progress_callback is not None:
        progress_callback({"event": "validation_started", "relative_fov": fov.relative_fov})
    output_pixel_um = fov.pixel_nm / 1000.0 / int(processing.get("mag", 2))
    raw_pixel_um = fov.pixel_nm / 1000.0
    for channel, paths, stack, raw_reference in zip(
        fov.channels, partial_outputs, stacks, raw_references, strict=True
    ):
        write_imagej_tiff(
            paths.stack,
            stack,
            axes="ZYX",
            pixel_size_um=output_pixel_um,
            z_spacing_um=fov.z_spacing_um,
            intensity_transform=intensity_transform,
        )
        write_imagej_tiff(
            paths.mip,
            np.max(stack, axis=0).astype(np.float32, copy=False),
            axes="YX",
            pixel_size_um=output_pixel_um,
            intensity_transform=intensity_transform,
        )
        write_imagej_tiff(
            paths.raw_max_mip,
            raw_reference,
            axes="YX",
            pixel_size_um=raw_pixel_um,
            reference_projection="frame_max_then_z_max",
        )
    shapes = validate_output_set(
        partial_outputs,
        fov.channels,
        intensity_transform=intensity_transform,
        expected_stacks=stacks,
        expected_raw_max_mips=raw_references,
    )
    for partial, final in zip(partial_outputs, fov.outputs, strict=True):
        for partial_path, final_path in zip(partial.all, final.all, strict=True):
            partial_path.replace(final_path)
    return _result_record(
        fov,
        "written",
        perf_counter() - started,
        shapes=shapes,
        input_shape=input_shape,
        intensity_transform=intensity_transform,
    )


def write_imagej_tiff(
    path: str | Path,
    image: np.ndarray,
    *,
    axes: str,
    pixel_size_um: float,
    z_spacing_um: float | None = None,
    intensity_transform: str | None = None,
    reference_projection: str | None = None,
) -> None:
    arr = np.asarray(image)
    if arr.ndim != len(axes):
        raise ValueError(f"Image ndim {arr.ndim} does not match axes {axes!r}")
    metadata: dict[str, object] = {
        "axes": axes,
        "unit": "um",
    }
    if intensity_transform is not None:
        metadata["intensity_transform"] = intensity_transform
    if reference_projection is not None:
        metadata["reference_projection"] = reference_projection
    if "Z" in axes and z_spacing_um is not None:
        metadata["spacing"] = float(z_spacing_um)
    tifffile.imwrite(
        path,
        arr,
        imagej=True,
        bigtiff=arr.nbytes >= 4_000_000_000,
        photometric="minisblack",
        resolution=(1.0 / pixel_size_um, 1.0 / pixel_size_um),
        metadata=metadata,
    )


def validate_output_set(
    outputs: Iterable[ChannelOutputPaths],
    channels: Iterable[ChannelSpec],
    *,
    intensity_transform: str = "order_root",
    expected_stacks: Iterable[np.ndarray] | None = None,
    expected_raw_max_mips: Iterable[np.ndarray] | None = None,
    output_pixel_um: float | None = None,
    raw_pixel_um: float | None = None,
    z_spacing_um: float | None = None,
) -> dict:
    output_tuple = tuple(outputs)
    channel_tuple = tuple(channels)
    expected_stack_tuple = tuple(expected_stacks) if expected_stacks is not None else None
    expected_raw_tuple = (
        tuple(expected_raw_max_mips) if expected_raw_max_mips is not None else None
    )
    if len(output_tuple) != len(channel_tuple):
        raise ValueError("Output set does not match configured channels")
    shapes: dict[str, dict[str, list[int]]] = {}
    for index, (paths, channel) in enumerate(zip(output_tuple, channel_tuple, strict=True)):
        if paths.channel != channel.name:
            raise ValueError(f"Output channel {paths.channel!r} does not match {channel.name!r}")
        with tifffile.TiffFile(paths.stack) as tif:
            stack_shape, stack_axes, stack_dtype = tif.series[0].shape, tif.series[0].axes, tif.series[0].dtype
            stack_metadata = tif.imagej_metadata or {}
            stack_pixel_um = _pixel_size_um(tif)
        with tifffile.TiffFile(paths.mip) as tif:
            mip_shape, mip_axes, mip_dtype = tif.series[0].shape, tif.series[0].axes, tif.series[0].dtype
            mip_metadata = tif.imagej_metadata or {}
            mip_pixel_um = _pixel_size_um(tif)
        with tifffile.TiffFile(paths.raw_max_mip) as tif:
            raw_shape, raw_axes, raw_dtype = tif.series[0].shape, tif.series[0].axes, tif.series[0].dtype
            raw_metadata = tif.imagej_metadata or {}
            raw_observed_pixel_um = _pixel_size_um(tif)
        if stack_axes != "ZYX" or mip_axes != "YX" or raw_axes != "YX":
            raise ValueError(f"Unexpected axes for {channel.name}: {stack_axes}, {mip_axes}, {raw_axes}")
        if stack_dtype != np.dtype("float32") or mip_dtype != np.dtype("float32"):
            raise ValueError(f"Unexpected SACD dtype for {channel.name}: {stack_dtype}, {mip_dtype}")
        if raw_dtype != np.dtype("uint16"):
            raise ValueError(f"Unexpected raw reference dtype for {channel.name}: {raw_dtype}")
        if (
            stack_metadata.get("intensity_transform") != intensity_transform
            or mip_metadata.get("intensity_transform") != intensity_transform
        ):
            raise ValueError(f"Incompatible intensity transform for channel {channel.name}")
        if raw_metadata.get("reference_projection") != "frame_max_then_z_max":
            raise ValueError(f"Incorrect raw reference provenance for channel {channel.name}")
        if output_pixel_um is not None and (
            not np.isclose(stack_pixel_um, output_pixel_um)
            or not np.isclose(mip_pixel_um, output_pixel_um)
        ):
            raise ValueError(f"Incorrect SACD pixel calibration for channel {channel.name}")
        if raw_pixel_um is not None and not np.isclose(raw_observed_pixel_um, raw_pixel_um):
            raise ValueError(f"Incorrect raw pixel calibration for channel {channel.name}")
        if z_spacing_um is not None and not np.isclose(
            float(stack_metadata.get("spacing", np.nan)), z_spacing_um
        ):
            raise ValueError(f"Incorrect z spacing for channel {channel.name}")
        stack = tifffile.imread(paths.stack)
        mip = tifffile.imread(paths.mip)
        raw_mip = tifffile.imread(paths.raw_max_mip)
        if not np.all(np.isfinite(stack)) or np.any(stack < 0):
            raise ValueError(f"Invalid SACD values for channel {channel.name}")
        if not np.array_equal(mip, np.max(stack, axis=0)):
            raise ValueError(f"Saved SACD MIP is not max(Z) for channel {channel.name}")
        if expected_stack_tuple is not None and not np.array_equal(stack, expected_stack_tuple[index]):
            raise ValueError(f"Saved SACD stack differs from expected values for channel {channel.name}")
        if expected_raw_tuple is not None and not np.array_equal(raw_mip, expected_raw_tuple[index]):
            raise ValueError(f"Saved raw reference differs from expected values for channel {channel.name}")
        shapes[channel.name] = {
            "stack_shape": list(stack_shape),
            "mip_shape": list(mip_shape),
            "raw_max_mip_shape": list(raw_shape),
        }
    return {"output_shapes": shapes}


def validate_fov_outputs(
    fov: FOVPlan,
    *,
    intensity_transform: str = "order_root",
    mag: int = 2,
    expected_raw_max_mips: Iterable[np.ndarray] | None = None,
) -> dict:
    return validate_output_set(
        fov.outputs,
        fov.channels,
        intensity_transform=intensity_transform,
        expected_raw_max_mips=expected_raw_max_mips,
        output_pixel_um=fov.pixel_nm / 1000.0 / mag,
        raw_pixel_um=fov.pixel_nm / 1000.0,
        z_spacing_um=fov.z_spacing_um,
    )


def validate_dataset_outputs(
    config: dict,
    *,
    verify_raw_references: bool = False,
    progress_callback: Callable[[dict], None] | None = None,
) -> list[dict]:
    """Independently validate every new-format file in a configured dataset."""
    plan = build_batch_plan(config)
    processing = config["processing"]
    intensity_transform = str(processing.get("intensity_transform", "order_root"))
    mag = int(processing.get("mag", 2))
    validations: list[dict] = []
    for index, fov in enumerate(plan.fovs, start=1):
        raw_references: list[np.ndarray] | None = None
        if verify_raw_references:
            maxima: list[np.ndarray | None] = [None for _ in fov.channels]
            for movie in fov.movies:
                raw = tifffile.imread(movie.path)
                for channel_index, channel in enumerate(fov.channels):
                    selected = select_channel_frames(
                        raw, movie.frame_ranges[channel_index], channel.camera_half
                    )
                    frame_max = np.max(selected, axis=0)
                    maxima[channel_index] = (
                        frame_max
                        if maxima[channel_index] is None
                        else np.maximum(maxima[channel_index], frame_max)
                    )
            raw_references = [_require_uint16(image) for image in maxima]
        validation = {
            "relative_fov": fov.relative_fov,
            **validate_fov_outputs(
                fov,
                intensity_transform=intensity_transform,
                mag=mag,
                expected_raw_max_mips=raw_references,
            ),
        }
        validations.append(validation)
        _emit(
            progress_callback,
            {
                "event": "validation_complete",
                "relative_fov": fov.relative_fov,
                "completed_fovs": index,
                "total_fovs": len(plan.fovs),
                "validated_files": index * len(fov.channels) * 3,
            },
        )
    return validations


def _pixel_size_um(tif: tifffile.TiffFile) -> float:
    metadata = tif.imagej_metadata or {}
    if metadata.get("unit") not in {"um", "micron"}:
        raise ValueError(f"TIFF calibration unit is not micrometers: {metadata.get('unit')!r}")

    def resolution(tag_name: str) -> float:
        value = tif.pages[0].tags[tag_name].value
        return float(value[0]) / float(value[1]) if isinstance(value, tuple) else float(value)

    x_pixels_per_um = resolution("XResolution")
    y_pixels_per_um = resolution("YResolution")
    if not np.isclose(x_pixels_per_um, y_pixels_per_um):
        raise ValueError("TIFF X/Y pixel calibrations differ")
    return 1.0 / x_pixels_per_um


def validate_legacy_ome_pair(
    stack_path: str | Path,
    mip_path: str | Path,
    channels: Iterable[ChannelSpec],
) -> tuple[np.ndarray, np.ndarray]:
    """Validate one historical raw-cumulant OME pair before explicit migration."""
    channel_tuple = tuple(channels)
    with tifffile.TiffFile(stack_path) as tif:
        stack_axes = tif.series[0].axes
        stack = tif.asarray()
        ome_metadata = tif.ome_metadata or ""
    with tifffile.TiffFile(mip_path) as tif:
        mip_axes = tif.series[0].axes
        mip = tif.asarray()
    if stack_axes != "CZYX" or mip_axes != "CYX":
        raise ValueError(f"Unexpected legacy OME axes: {stack_axes}, {mip_axes}")
    if stack.dtype != np.float32 or mip.dtype != np.float32:
        raise ValueError(f"Unexpected legacy OME dtypes: {stack.dtype}, {mip.dtype}")
    if stack.shape[0] != len(channel_tuple) or mip.shape[0] != len(channel_tuple):
        raise ValueError("Legacy OME channel count does not match the current configuration")
    if not np.all(np.isfinite(stack)) or np.any(stack < 0):
        raise ValueError("Legacy OME stack contains nonfinite or negative values")
    if not np.array_equal(mip, np.max(stack, axis=1)):
        raise ValueError("Legacy OME MIP is not exactly max(Z) of its stack")
    missing_names = [channel.name for channel in channel_tuple if channel.name not in ome_metadata]
    if missing_names:
        raise ValueError(f"Legacy OME metadata is missing channel name(s): {missing_names}")
    return stack, mip


def migrate_legacy_ome_dataset(
    config: dict,
    config_path: str | Path,
    *,
    progress_callback: Callable[[dict], None] | None = None,
) -> list[dict]:
    """Explicitly migrate one configured raw-cumulant OME dataset without rerunning SACD."""
    plan = build_batch_plan(config)
    processing = config["processing"]
    intensity_transform = str(processing.get("intensity_transform", "order_root"))
    if intensity_transform != "order_root":
        raise ValueError("Legacy OME migration requires intensity_transform='order_root'")
    if any(path.exists() for fov in plan.fovs for path in fov.all_outputs):
        raise FileExistsError("New-format outputs already exist; migration will not overwrite them")
    legacy_pairs = [(fov, legacy_ome_paths(fov)) for fov in plan.fovs]
    missing = [str(path) for _, pair in legacy_pairs for path in pair if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing legacy OME input(s): {missing}")

    staging_root = plan.output_root / ".order_root_migration"
    if staging_root.exists():
        raise FileExistsError(
            f"Migration staging directory already exists: {staging_root}. "
            "Inspect or remove it before retrying."
        )
    staging_root.mkdir(parents=True)
    staged_sets: list[tuple[FOVPlan, tuple[ChannelOutputPaths, ...], list[np.ndarray], list[np.ndarray]]] = []
    results: list[dict] = []
    try:
        for fov_index, (fov, (legacy_stack_path, legacy_mip_path)) in enumerate(
            legacy_pairs, start=1
        ):
            raw_stack, _ = validate_legacy_ome_pair(
                legacy_stack_path, legacy_mip_path, fov.channels
            )
            params = SACDParams(
                ac_order=int(processing.get("ac_order", 2)),
                intensity_transform="order_root",
            )
            params.validate_core()
            transformed = apply_intensity_transform(raw_stack, params)
            stacks = [transformed[index] for index in range(len(fov.channels))]
            raw_max_mips: list[np.ndarray | None] = [None for _ in fov.channels]
            for movie in fov.movies:
                raw = tifffile.imread(movie.path)
                for channel_index, channel in enumerate(fov.channels):
                    selected = select_channel_frames(
                        raw, movie.frame_ranges[channel_index], channel.camera_half
                    )
                    frame_max = np.max(selected, axis=0)
                    raw_max_mips[channel_index] = (
                        frame_max
                        if raw_max_mips[channel_index] is None
                        else np.maximum(raw_max_mips[channel_index], frame_max)
                    )
            raw_references = [_require_uint16(image) for image in raw_max_mips]
            staged_outputs = output_paths(staging_root, fov.prefix, fov.channels)
            output_pixel_um = fov.pixel_nm / 1000.0 / int(processing.get("mag", 2))
            raw_pixel_um = fov.pixel_nm / 1000.0
            for paths, stack, raw_reference in zip(
                staged_outputs, stacks, raw_references, strict=True
            ):
                write_imagej_tiff(
                    paths.stack,
                    stack,
                    axes="ZYX",
                    pixel_size_um=output_pixel_um,
                    z_spacing_um=fov.z_spacing_um,
                    intensity_transform="order_root",
                )
                write_imagej_tiff(
                    paths.mip,
                    np.max(stack, axis=0).astype(np.float32, copy=False),
                    axes="YX",
                    pixel_size_um=output_pixel_um,
                    intensity_transform="order_root",
                )
                write_imagej_tiff(
                    paths.raw_max_mip,
                    raw_reference,
                    axes="YX",
                    pixel_size_um=raw_pixel_um,
                    reference_projection="frame_max_then_z_max",
                )
            shapes = validate_output_set(
                staged_outputs,
                fov.channels,
                intensity_transform="order_root",
                expected_stacks=stacks,
                expected_raw_max_mips=raw_references,
                output_pixel_um=output_pixel_um,
                raw_pixel_um=raw_pixel_um,
                z_spacing_um=fov.z_spacing_um,
            )
            staged_sets.append((fov, staged_outputs, stacks, raw_references))
            result = _result_record(
                fov,
                "migrated_order_root",
                0.0,
                shapes=shapes,
                input_shape=fov.movies[0].shape,
                intensity_transform="order_root",
            )
            results.append(result)
            _emit(
                progress_callback,
                {
                    "event": "migration_staged",
                    "relative_fov": fov.relative_fov,
                    "completed_fovs": fov_index,
                    "total_fovs": len(plan.fovs),
                },
            )

        # Publication starts only after every staged file has passed exact validation.
        for fov, staged_outputs, _, _ in staged_sets:
            for staged, final in zip(staged_outputs, fov.outputs, strict=True):
                for staged_path, final_path in zip(staged.all, final.all, strict=True):
                    os.replace(staged_path, final_path)
        for fov, _, stacks, raw_references in staged_sets:
            validate_output_set(
                fov.outputs,
                fov.channels,
                intensity_transform="order_root",
                expected_stacks=stacks,
                expected_raw_max_mips=raw_references,
                output_pixel_um=fov.pixel_nm / 1000.0 / int(processing.get("mag", 2)),
                raw_pixel_um=fov.pixel_nm / 1000.0,
                z_spacing_um=fov.z_spacing_um,
            )

        _prepare_provenance(plan, config, Path(config_path))
        manifest = _load_manifest(plan.output_root)
        manifest["fovs"] = results
        manifest["exclusions"] = list(plan.exclusions)
        manifest["intensity_transform"] = "order_root"
        manifest["output_format"] = "separate_imagej_tiff_v1"
        manifest["migration"] = {
            "completed_at": utc_now(),
            "source_intensity_transform": "raw_cumulant",
            "transform": f"power_1_over_{int(processing.get('ac_order', 2))}",
            "legacy_ome_files_removed": len(legacy_pairs) * 2,
        }
        manifest["updated_at"] = utc_now()
        _write_json(plan.output_root / "_processing" / "manifest.json", manifest)
        for _, pair in legacy_pairs:
            for path in pair:
                path.unlink()
        shutil.rmtree(staging_root)
        return results
    except Exception:
        # Keep staged files for diagnosis; legacy inputs are retained until the success path above.
        raise


def migrate_output_names(
    config: dict,
    config_path: str | Path,
    *,
    progress_callback: Callable[[dict], None] | None = None,
) -> dict:
    """Atomically shorten existing output names and update their manifest records."""
    plan = build_batch_plan(config)
    aliases = config.get("output_prefix_aliases", {})
    processing_dir = plan.output_root / "_processing"
    manifest_path = processing_dir / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Missing output manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text())
    manifest_by_fov = {
        item.get("relative_fov"): item for item in manifest.get("fovs", ())
    }

    pairs: list[tuple[Path, Path, str]] = []
    for fov in plan.fovs:
        old_prefix = str(aliases.get(fov.relative_fov, fov.movies[0].prefix))
        old_outputs = output_paths(plan.output_root, old_prefix, fov.channels)
        for old_group, new_group in zip(old_outputs, fov.outputs, strict=True):
            for old_path, new_path in zip(old_group.all, new_group.all, strict=True):
                if old_path != new_path:
                    pairs.append((old_path, new_path, fov.relative_fov))

    missing = [str(old) for old, _, _ in pairs if not old.is_file()]
    collisions = [str(new) for old, new, _ in pairs if new.exists() and new != old]
    new_paths = [new for _, new, _ in pairs]
    if len(new_paths) != len(set(new_paths)):
        raise ValueError("Short-name migration contains duplicate destination paths")
    if missing or collisions:
        raise ValueError(
            f"Output rename preflight failed; missing={missing}, collisions={collisions}"
        )

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup_path = processing_dir / f"manifest.before_short_names.{timestamp}.json"
    mapping_path = processing_dir / f"short_name_mapping.{timestamp}.json"
    checksums = {str(old): _sha256_file(old) for old, _, _ in pairs}
    mapping = {
        "created_at": utc_now(),
        "config_source": str(config_path),
        "files": [
            {
                "relative_fov": relative_fov,
                "old": str(old),
                "new": str(new),
                "sha256": checksums[str(old)],
            }
            for old, new, relative_fov in pairs
        ],
    }
    shutil.copy2(manifest_path, backup_path)
    _write_json(mapping_path, mapping)

    moved: list[tuple[Path, Path]] = []
    try:
        for index, (old, new, relative_fov) in enumerate(pairs, start=1):
            os.replace(old, new)
            moved.append((old, new))
            observed = _sha256_file(new)
            if observed != checksums[str(old)]:
                raise ValueError(f"Checksum changed while renaming {old} to {new}")
            _emit(
                progress_callback,
                {
                    "event": "output_renamed",
                    "relative_fov": relative_fov,
                    "completed_files": index,
                    "total_files": len(pairs),
                },
            )

        for fov in plan.fovs:
            record = manifest_by_fov.get(fov.relative_fov)
            if record is None:
                raise ValueError(f"Manifest has no FOV record for {fov.relative_fov}")
            record["outputs"] = _outputs_record(fov)
            record["updated_at"] = utc_now()
            validate_fov_outputs(
                fov,
                intensity_transform=str(
                    config.get("processing", {}).get("intensity_transform", "order_root")
                ),
                mag=int(config.get("processing", {}).get("mag", 2)),
            )
        manifest["naming"] = {
            "prefix_source": "fov_folder_name",
            "migrated_at": utc_now(),
            "mapping": str(mapping_path),
            "previous_manifest": str(backup_path),
        }
        manifest["updated_at"] = utc_now()
        _write_json(manifest_path, manifest)
        _write_json(processing_dir / "config.json", config)
    except BaseException:
        for old, new in reversed(moved):
            if new.exists() and not old.exists():
                os.replace(new, old)
        shutil.copy2(backup_path, manifest_path)
        raise
    return {
        "renamed_files": len(pairs),
        "manifest": str(manifest_path),
        "previous_manifest": str(backup_path),
        "mapping": str(mapping_path),
    }


@pooled_batch
def run_batch(
    config: dict,
    config_path: str | Path,
    *,
    max_new_fovs: int | None = None,
    progress_callback: Callable[[dict], None] | None = None,
    executor: ProcessPoolExecutor | None = None,
) -> list[dict]:
    plan = build_batch_plan(config)
    _prepare_provenance(plan, config, Path(config_path))
    processing = config["processing"]
    intensity_transform = str(processing.get("intensity_transform", "order_root"))
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
        if previous and previous.get("status") in {
            "written", "skipped_existing", "migrated_order_root"
        }:
            if previous.get("intensity_transform") != intensity_transform:
                raise ValueError(
                    f"Existing FOV {fov.relative_fov} uses intensity_transform="
                    f"{previous.get('intensity_transform')!r}, expected {intensity_transform!r}. "
                    "Use a new output location; untagged and raw-cumulant outputs cannot be resumed."
                )
            try:
                validate_fov_outputs(
                    fov,
                    intensity_transform=intensity_transform,
                    mag=int(processing.get("mag", 2)),
                )
            except (FileNotFoundError, ValueError):
                previous = None
        if previous and previous.get("status") in {
            "written", "skipped_existing", "migrated_order_root"
        }:
            result = {**previous, "status": "resumed_manifest"}
            results.append(result)
            completed_reconstructions += fov.reconstruction_count
            _emit(progress_callback, result)
            continue
        if max_new_fovs is not None and newly_written >= max_new_fovs:
            break

        step_callback = FOVEvents(progress_callback)
        try:
            result = process_fov(
                fov,
                processing,
                progress_callback=step_callback,
                executor=executor,
            )
            newly_written += result["status"] == "written"
        except Exception as exc:
            result = _result_record(
                fov,
                "failed",
                0.0,
                error=repr(exc),
                intensity_transform=intensity_transform,
            )
        if result["status"] == "failed":
            result["failure_stage"] = step_callback.failure_stage
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
                "intensity_transform": str(
                    config.get("processing", {}).get("intensity_transform", "order_root")
                ),
                "output_format": "separate_imagej_tiff_v1",
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
    intensity_transform: str | None = None,
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
        "intensity_transform": intensity_transform,
        "output_format": "separate_imagej_tiff_v1",
        "outputs": _outputs_record(fov),
        "runtime_s": runtime_s,
        "status": status,
        "error": error,
    }


def _outputs_record(fov: FOVPlan) -> dict[str, dict[str, str]]:
    return {
        output.channel: {
            "stack": str(output.stack),
            "mip": str(output.mip),
            "raw_max_mip": str(output.raw_max_mip),
        }
        for output in fov.outputs
    }


def _sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _partial_path(path: str | Path) -> Path:
    path = Path(path)
    return path.with_name(f"{path.stem}.partial{path.suffix}")


def _require_uint16(image: np.ndarray | None) -> np.ndarray:
    if image is None:
        raise ValueError("Raw reference was not constructed")
    arr = np.asarray(image)
    if arr.dtype != np.uint16:
        raise ValueError(f"Raw reference source must be uint16, got {arr.dtype}")
    return arr


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
