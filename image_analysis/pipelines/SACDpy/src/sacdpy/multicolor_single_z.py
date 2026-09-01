from __future__ import annotations

import json
import platform
import re
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from time import perf_counter
from typing import Callable, Iterable

import numpy as np
import tifffile

from .multicolor_zstack import select_channel_frames, write_imagej_tiff
from .params import SACDParams
from .reconstruction import reconstruct


_ONI_SINGLE_Z_RE = re.compile(
    r"^(?P<prefix>.+)_posXY(?P<pos_xy>\d+)_channels_t(?P<time_index>\d+)_posZ(?P<z_index>\d+)\.tiff?$",
    re.IGNORECASE,
)
_OUTPUT_FORMAT = "separate_imagej_tiff_single_z_v1"
_PIPELINE = "multicolor_single_z"


@dataclass(frozen=True)
class ChannelSpec:
    name: str
    wavelength_nm: float
    metadata_laser_nm: float
    camera_half: str


@dataclass(frozen=True)
class MovieInfo:
    path: Path
    file_prefix: str
    pos_xy: int
    time_index: int
    z_index: int
    shape: tuple[int, int, int]
    dtype: str
    metadata: dict
    frame_ranges: tuple[tuple[int, int], ...]
    active_lasers_nm: tuple[float, ...]


@dataclass(frozen=True)
class ChannelOutputPaths:
    channel: str
    sacd: Path
    mip: Path

    @property
    def all(self) -> tuple[Path, Path]:
        return self.sacd, self.mip


@dataclass(frozen=True)
class FOVPlan:
    raw_root: Path
    output_root: Path
    fov_folder: Path
    relative_fov: str
    prefix: str
    movie: MovieInfo
    channels: tuple[ChannelSpec, ...]
    outputs: tuple[ChannelOutputPaths, ...]
    pixel_nm: float
    na: float

    @property
    def reconstruction_count(self) -> int:
        return len(self.channels)

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


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _channel_specs(config: dict) -> tuple[ChannelSpec, ...]:
    entries = config.get("processing", {}).get("channels", ())
    if not entries:
        raise ValueError("processing.channels must contain at least one channel")
    channels: list[ChannelSpec] = []
    names: set[str] = set()
    for entry in entries:
        name = str(entry["name"])
        if not name:
            raise ValueError("Channel names must not be empty")
        if name in names:
            raise ValueError(f"Duplicate channel name: {name}")
        names.add(name)
        camera_half = str(entry["camera_half"]).lower()
        if camera_half not in {"left", "right"}:
            raise ValueError(f"camera_half for channel {name} must be left or right")
        channels.append(
            ChannelSpec(
                name=name,
                wavelength_nm=float(entry["wavelength_nm"]),
                metadata_laser_nm=float(entry["metadata_laser_nm"]),
                camera_half=camera_half,
            )
        )
    return tuple(channels)


def validate_config(config: dict) -> dict:
    if config.get("version") != 1:
        raise ValueError("Multicolor single-z config version must be 1")
    if not config.get("raw_root") or not config.get("output_root"):
        raise ValueError("Config must define raw_root and output_root")
    _channel_specs(config)
    return config


def parse_movie_filename(path: str | Path) -> tuple[str, int, int, int]:
    path = Path(path)
    match = _ONI_SINGLE_Z_RE.match(path.name)
    if match is None:
        raise ValueError(f"Filename does not match ONI multicolor single-z grammar: {path.name}")
    return (
        match.group("prefix"),
        int(match.group("pos_xy")),
        int(match.group("time_index")),
        int(match.group("z_index")),
    )


def _active_laser_nm(step: object, wavelengths: tuple[float, ...], path: Path) -> float:
    if not isinstance(step, dict):
        raise ValueError(f"Laser-program step is not an object: {path}")
    states = step.get("states")
    if not isinstance(states, list) or not states:
        raise ValueError(f"Laser-program step has no states: {path}")
    active_indices: set[int] = set()
    for state in states:
        values = state.get("values") if isinstance(state, dict) else None
        if not isinstance(values, list) or len(values) != len(wavelengths):
            raise ValueError(f"Laser-program state does not match LaserWavelength_nm: {path}")
        active_indices.update(index for index, value in enumerate(values) if float(value) > 0)
    if len(active_indices) != 1:
        raise ValueError(
            f"Expected exactly one active laser per acquisition step, found {sorted(active_indices)}: {path}"
        )
    return wavelengths[next(iter(active_indices))]


def read_movie_info(path: str | Path, channels: Iterable[ChannelSpec]) -> MovieInfo:
    path = Path(path)
    file_prefix, pos_xy, time_index, z_index = parse_movie_filename(path)
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
    wavelength_values = metadata.get("LaserWavelength_nm")
    steps = metadata.get("laserProgram", {}).get("steps")
    if not isinstance(wavelength_values, list) or not isinstance(steps, list):
        raise ValueError(f"Missing LaserWavelength_nm or laserProgram.steps metadata: {path}")
    wavelengths = tuple(float(value) for value in wavelength_values)
    if len(steps) != len(channel_tuple):
        raise ValueError(
            f"Expected {len(channel_tuple)} acquisition step(s), got {len(steps)}: {path}"
        )

    frame_ranges: list[tuple[int, int]] = []
    active_lasers: list[float] = []
    start = 0
    for index, (channel, step) in enumerate(zip(channel_tuple, steps, strict=True)):
        repeats = step.get("nRepeats") if isinstance(step, dict) else None
        if not isinstance(repeats, int) or repeats <= 0:
            raise ValueError(f"Invalid nRepeats for channel {channel.name}: {path}")
        observed_laser = _active_laser_nm(step, wavelengths, path)
        if not np.isclose(observed_laser, channel.metadata_laser_nm):
            raise ValueError(
                f"Acquisition step {index} ({channel.name}) expected active metadata laser "
                f"{channel.metadata_laser_nm:g} nm, got {observed_laser:g} nm: {path}"
            )
        frame_ranges.append((start, start + repeats))
        active_lasers.append(observed_laser)
        start += repeats
    if start != page_count:
        raise ValueError(
            f"Laser-program repeats total {start}, but TIFF contains {page_count} page(s): {path}"
        )
    if width % 2:
        raise ValueError(f"Camera-split TIFF width must be even, got {width}: {path}")

    return MovieInfo(
        path=path,
        file_prefix=file_prefix,
        pos_xy=pos_xy,
        time_index=time_index,
        z_index=z_index,
        shape=(page_count, height, width),
        dtype=dtype,
        metadata=metadata,
        frame_ranges=tuple(frame_ranges),
        active_lasers_nm=tuple(active_lasers),
    )


def output_paths(
    output_root: str | Path,
    prefix: str,
    channels: Iterable[ChannelSpec],
) -> tuple[ChannelOutputPaths, ...]:
    root = Path(output_root)
    return tuple(
        ChannelOutputPaths(
            channel=channel.name,
            sacd=root / f"{prefix}__{channel.name}-SACD.tif",
            mip=root / f"{prefix}__{channel.name}-MIP.tif",
        )
        for channel in channels
    )


def build_batch_plan(config: dict) -> BatchPlan:
    validate_config(config)
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
    fovs: list[FOVPlan] = []
    exclusions: list[dict[str, str]] = []
    planned_paths: set[Path] = set()

    for position_dir in sorted(path for path in raw_root.rglob(position_folder) if path.is_dir()):
        fov_folder = position_dir.parent
        relative_fov = fov_folder.relative_to(raw_root).as_posix()
        if selected and relative_fov not in selected:
            exclusions.append(_exclusion(relative_fov, "not selected by explicit FOV list"))
            continue
        if relative_fov in excluded_config:
            exclusions.append(_exclusion(relative_fov, str(excluded_config[relative_fov])))
            continue

        matching = [
            path
            for path in sorted(position_dir.glob(glob_pattern))
            if path.is_file() and _ONI_SINGLE_Z_RE.match(path.name) and "sacdpy" not in path.name.lower()
        ]
        if not matching:
            exclusions.append(_exclusion(relative_fov, "no matching ONI multicolor single-z TIFF"))
            continue
        if len(matching) != 1:
            raise ValueError(
                f"Expected exactly one ONI TIFF/z-plane in {position_dir}, found {len(matching)}"
            )
        movie = read_movie_info(matching[0], channels)
        if movie.time_index != 0:
            raise ValueError(f"Multicolor single-z pipeline expects t0 only, got t{movie.time_index}: {position_dir}")

        prefix = fov_folder.name
        outputs = output_paths(output_root, prefix, channels)
        for path in (path for output in outputs for path in output.all):
            if path in planned_paths:
                raise ValueError(f"Output path collision: {path}")
            planned_paths.add(path)

        pixel_nm = processing.get("pixel_nm")
        if pixel_nm is None:
            value = movie.metadata.get("PixelSize_um")
            pixel_nm = float(value) * 1000.0 if value is not None else processing.get("fallback_pixel_nm")
        na = processing.get("na")
        if na is None:
            na = movie.metadata.get("Objective_NA", processing.get("fallback_na"))
        if pixel_nm is None or na is None:
            raise ValueError(f"Could not resolve pixel size or NA for {position_dir}")

        fovs.append(
            FOVPlan(
                raw_root=raw_root,
                output_root=output_root,
                fov_folder=fov_folder,
                relative_fov=relative_fov,
                prefix=prefix,
                movie=movie,
                channels=channels,
                outputs=outputs,
                pixel_nm=float(pixel_nm),
                na=float(na),
            )
        )
    return BatchPlan(raw_root, output_root, tuple(fovs), tuple(exclusions))


def preflight_summary(plan: BatchPlan) -> dict:
    return {
        "raw_root": str(plan.raw_root),
        "output_root": str(plan.output_root),
        "fovs": len(plan.fovs),
        "movies": len(plan.fovs),
        "channels": [channel.name for channel in plan.fovs[0].channels] if plan.fovs else [],
        "reconstructions": plan.reconstruction_count,
        "outputs": sum(len(fov.all_outputs) for fov in plan.fovs),
        "exclusions": len(plan.exclusions),
    }


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


def process_fov(
    fov: FOVPlan,
    processing: dict,
    *,
    progress_callback: Callable[[dict], None] | None = None,
) -> dict:
    intensity_transform = str(processing.get("intensity_transform", "order_root"))
    existing = [path for path in fov.all_outputs if path.exists()]
    if existing:
        if len(existing) != len(fov.all_outputs):
            raise FileExistsError(f"Partial final multicolor single-z output set exists: {existing}")
        shapes = validate_fov_outputs(
            fov,
            intensity_transform=intensity_transform,
            mag=int(processing.get("mag", 2)),
        )
        return _result_record(fov, "skipped_existing", 0.0, shapes=shapes, intensity_transform=intensity_transform)

    fov.output_root.mkdir(parents=True, exist_ok=True)
    partial_outputs = tuple(
        ChannelOutputPaths(output.channel, _partial_path(output.sacd), _partial_path(output.mip))
        for output in fov.outputs
    )
    for partial in (path for output in partial_outputs for path in output.all):
        if partial.exists():
            partial.unlink()

    started = perf_counter()
    raw = tifffile.imread(fov.movie.path)
    sacd_images: list[np.ndarray] = []
    time_mips: list[np.ndarray] = []
    for index, channel in enumerate(fov.channels):
        reconstruction_started = perf_counter()
        selected = select_channel_frames(raw, fov.movie.frame_ranges[index], channel.camera_half)
        time_mips.append(_require_uint16(np.max(selected, axis=0)))
        sacd = reconstruct(selected, _params_for_channel(fov, channel, processing))
        if sacd.ndim != 2:
            raise ValueError(f"Expected one 2D SACD image, got {sacd.shape} for {fov.movie.path}")
        sacd_images.append(np.asarray(sacd, dtype=np.float32))
        _emit(
            progress_callback,
            {
                "event": "reconstruction_done",
                "relative_fov": fov.relative_fov,
                "z_index": fov.movie.z_index,
                "channel": channel.name,
                "completed_in_fov": index + 1,
                "total_in_fov": fov.reconstruction_count,
                "runtime_s": perf_counter() - reconstruction_started,
            },
        )

    output_pixel_um = fov.pixel_nm / 1000.0 / int(processing.get("mag", 2))
    raw_pixel_um = fov.pixel_nm / 1000.0
    for paths, sacd, time_mip in zip(partial_outputs, sacd_images, time_mips, strict=True):
        write_imagej_tiff(
            paths.sacd,
            sacd,
            axes="YX",
            pixel_size_um=output_pixel_um,
            intensity_transform=intensity_transform,
        )
        write_imagej_tiff(
            paths.mip,
            time_mip,
            axes="YX",
            pixel_size_um=raw_pixel_um,
            reference_projection="frame_max",
        )
    shapes = validate_output_set(
        partial_outputs,
        fov.channels,
        intensity_transform=intensity_transform,
        expected_sacd=sacd_images,
        expected_mips=time_mips,
        output_pixel_um=output_pixel_um,
        raw_pixel_um=raw_pixel_um,
    )
    for partial, final in zip(partial_outputs, fov.outputs, strict=True):
        for partial_path, final_path in zip(partial.all, final.all, strict=True):
            partial_path.replace(final_path)
    return _result_record(
        fov,
        "written",
        perf_counter() - started,
        shapes=shapes,
        input_shape=raw.shape,
        intensity_transform=intensity_transform,
    )


def validate_output_set(
    outputs: Iterable[ChannelOutputPaths],
    channels: Iterable[ChannelSpec],
    *,
    intensity_transform: str = "order_root",
    expected_sacd: Iterable[np.ndarray] | None = None,
    expected_mips: Iterable[np.ndarray] | None = None,
    output_pixel_um: float | None = None,
    raw_pixel_um: float | None = None,
) -> dict:
    output_tuple = tuple(outputs)
    channel_tuple = tuple(channels)
    expected_sacd_tuple = tuple(expected_sacd) if expected_sacd is not None else None
    expected_mip_tuple = tuple(expected_mips) if expected_mips is not None else None
    if len(output_tuple) != len(channel_tuple):
        raise ValueError("Output set does not match configured channels")
    shapes: dict[str, dict[str, list[int]]] = {}
    for index, (paths, channel) in enumerate(zip(output_tuple, channel_tuple, strict=True)):
        if paths.channel != channel.name:
            raise ValueError(f"Output channel {paths.channel!r} does not match {channel.name!r}")
        with tifffile.TiffFile(paths.sacd) as tif:
            sacd_shape, sacd_axes, sacd_dtype = tif.series[0].shape, tif.series[0].axes, tif.series[0].dtype
            sacd_metadata = tif.imagej_metadata or {}
            sacd_pixel_um = _pixel_size_um(tif)
        with tifffile.TiffFile(paths.mip) as tif:
            mip_shape, mip_axes, mip_dtype = tif.series[0].shape, tif.series[0].axes, tif.series[0].dtype
            mip_metadata = tif.imagej_metadata or {}
            mip_pixel_um = _pixel_size_um(tif)
        if sacd_axes != "YX" or mip_axes != "YX":
            raise ValueError(f"Unexpected axes for {channel.name}: {sacd_axes}, {mip_axes}")
        if sacd_dtype != np.dtype("float32") or mip_dtype != np.dtype("uint16"):
            raise ValueError(f"Unexpected output dtypes for {channel.name}: {sacd_dtype}, {mip_dtype}")
        if sacd_metadata.get("intensity_transform") != intensity_transform:
            raise ValueError(f"Incompatible intensity transform for channel {channel.name}")
        if mip_metadata.get("reference_projection") != "frame_max":
            raise ValueError(f"Incorrect time-MIP provenance for channel {channel.name}")
        if output_pixel_um is not None and not np.isclose(sacd_pixel_um, output_pixel_um):
            raise ValueError(f"Incorrect SACD pixel calibration for channel {channel.name}")
        if raw_pixel_um is not None and not np.isclose(mip_pixel_um, raw_pixel_um):
            raise ValueError(f"Incorrect MIP pixel calibration for channel {channel.name}")
        sacd = tifffile.imread(paths.sacd)
        mip = tifffile.imread(paths.mip)
        if not np.all(np.isfinite(sacd)) or np.any(sacd < 0):
            raise ValueError(f"Invalid SACD values for channel {channel.name}")
        if expected_sacd_tuple is not None and not np.array_equal(sacd, expected_sacd_tuple[index]):
            raise ValueError(f"Saved SACD differs from expected values for channel {channel.name}")
        if expected_mip_tuple is not None and not np.array_equal(mip, expected_mip_tuple[index]):
            raise ValueError(f"Saved MIP differs from expected values for channel {channel.name}")
        shapes[channel.name] = {"sacd_shape": list(sacd_shape), "mip_shape": list(mip_shape)}
    return {"output_shapes": shapes}


def validate_fov_outputs(
    fov: FOVPlan,
    *,
    intensity_transform: str = "order_root",
    mag: int = 2,
) -> dict:
    return validate_output_set(
        fov.outputs,
        fov.channels,
        intensity_transform=intensity_transform,
        output_pixel_um=fov.pixel_nm / 1000.0 / mag,
        raw_pixel_um=fov.pixel_nm / 1000.0,
    )


def run_batch(
    config: dict,
    config_source: str = "notebook settings cell",
    *,
    max_new_fovs: int | None = None,
    progress_callback: Callable[[dict], None] | None = None,
) -> list[dict]:
    plan = build_batch_plan(config)
    _prepare_provenance(plan, config, config_source)
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
        if previous and previous.get("status") in {"written", "skipped_existing"}:
            expected_provenance = {
                "channels": [asdict(channel) for channel in fov.channels],
                "frame_ranges_zero_based_end_exclusive": [
                    list(value) for value in fov.movie.frame_ranges
                ],
                "active_lasers_nm": list(fov.movie.active_lasers_nm),
                "pixel_nm": fov.pixel_nm,
                "na": fov.na,
                "intensity_transform": intensity_transform,
                "pipeline": _PIPELINE,
                "output_format": _OUTPUT_FORMAT,
            }
            if any(previous.get(key) != value for key, value in expected_provenance.items()):
                raise ValueError(
                    f"Existing FOV {fov.relative_fov} uses incompatible provenance and cannot be resumed"
                )
            try:
                validate_fov_outputs(fov, intensity_transform=intensity_transform, mag=int(processing.get("mag", 2)))
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
                step_callback = lambda event: print("SACD_STEP " + json.dumps(event, default=str), flush=True)
            result = process_fov(fov, processing, progress_callback=step_callback)
            newly_written += result["status"] == "written"
        except Exception as exc:
            result = _result_record(
                fov,
                "failed",
                0.0,
                error=repr(exc),
                intensity_transform=intensity_transform,
            )
        results.append(result)
        manifest = _upsert_manifest(plan.output_root, manifest, result, plan.exclusions)

        completed_reconstructions += fov.reconstruction_count
        elapsed_processing += float(result.get("runtime_s", 0.0))
        measured = sum(item["reconstruction_count"] for item in results if item.get("status") == "written")
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


def _prepare_provenance(plan: BatchPlan, config: dict, config_source: str) -> None:
    processing_dir = plan.output_root / "_processing"
    processing_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = processing_dir / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("pipeline") != _PIPELINE or manifest.get("output_format") != _OUTPUT_FORMAT:
            raise ValueError(
                f"Output folder contains an incompatible manifest: {manifest_path}. Use a new output location."
            )
        expected_transform = str(config.get("processing", {}).get("intensity_transform", "order_root"))
        if manifest.get("intensity_transform") != expected_transform:
            raise ValueError(
                f"Output folder contains incompatible intensity provenance: {manifest_path}. "
                "Use a new output location."
            )
    else:
        _write_json(
            manifest_path,
            {
                "created_at": utc_now(),
                "config_source": config_source,
                "pipeline": _PIPELINE,
                "intensity_transform": str(config.get("processing", {}).get("intensity_transform", "order_root")),
                "output_format": _OUTPUT_FORMAT,
                "exclusions": list(plan.exclusions),
                "fovs": [],
            },
        )
    _write_json(processing_dir / "config.json", config)
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
        "source_file": str(fov.movie.path),
        "file_prefix": fov.movie.file_prefix,
        "output_prefix": fov.prefix,
        "z_index": fov.movie.z_index,
        "channels": [asdict(channel) for channel in fov.channels],
        "frame_ranges_zero_based_end_exclusive": [list(value) for value in fov.movie.frame_ranges],
        "active_lasers_nm": list(fov.movie.active_lasers_nm),
        "reconstruction_count": fov.reconstruction_count,
        "input_shape": list(input_shape) if input_shape is not None else None,
        **(shapes or {}),
        "pixel_nm": fov.pixel_nm,
        "na": fov.na,
        "intensity_transform": intensity_transform,
        "pipeline": _PIPELINE,
        "output_format": _OUTPUT_FORMAT,
        "outputs": {
            output.channel: {"sacd": str(output.sacd), "mip": str(output.mip)}
            for output in fov.outputs
        },
        "runtime_s": runtime_s,
        "status": status,
        "error": error,
    }


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


def _partial_path(path: str | Path) -> Path:
    path = Path(path)
    return path.with_name(f"{path.stem}.partial{path.suffix}")


def _require_uint16(image: np.ndarray) -> np.ndarray:
    arr = np.asarray(image)
    if arr.dtype != np.uint16:
        raise ValueError(f"Raw time-MIP source must be uint16, got {arr.dtype}")
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


def _emit(callback: Callable[[dict], None] | None, event: dict) -> None:
    if callback is not None:
        callback(event)
