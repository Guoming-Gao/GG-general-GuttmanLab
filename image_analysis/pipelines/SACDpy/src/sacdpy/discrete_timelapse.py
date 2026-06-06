from __future__ import annotations

import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Iterable

import numpy as np
import tifffile


_ONI_TIMELAPSE_RE = re.compile(
    r"^(?P<prefix>.+)_posXY(?P<pos_xy>\d+)_channels_t(?P<time_index>\d+)_posZ(?P<z_index>\d+)\.tiff?$",
    re.IGNORECASE,
)
_SECONDS_RE = re.compile(r"(?<!m)(?P<seconds>\d+(?:\.\d+)?)s(?:[-_/]|$)", re.IGNORECASE)


@dataclass(frozen=True)
class TimelapseFile:
    path: Path
    prefix: str
    channel: str
    pos_xy: int
    time_index: int
    z_index: int


@dataclass(frozen=True)
class TimelapseGroup:
    prefix: str
    channel: str
    pos_xy: int
    files: dict[tuple[int, int], Path]

    @property
    def time_indices(self) -> tuple[int, ...]:
        return tuple(sorted({time_index for time_index, _ in self.files}))

    @property
    def z_indices(self) -> tuple[int, ...]:
        return tuple(sorted({z_index for _, z_index in self.files}))

    @property
    def missing_entries(self) -> tuple[tuple[int, int], ...]:
        return tuple(
            (time_index, z_index)
            for time_index in self.time_indices
            for z_index in self.z_indices
            if (time_index, z_index) not in self.files
        )

    @property
    def is_complete(self) -> bool:
        return not self.missing_entries


@dataclass(frozen=True)
class FolderTimelapsePlan:
    input_folder: Path
    output_dir: Path
    files: tuple[TimelapseFile, ...]
    groups: tuple[TimelapseGroup, ...]


def parse_timelapse_filename(path: str | Path, channel_name: str = "647") -> TimelapseFile:
    """Parse ONI timelapse-z filenames ending in `_posXY#_channels_t#_posZ#.tif`."""

    path = Path(path)
    match = _ONI_TIMELAPSE_RE.match(path.name)
    if match is None:
        raise ValueError(f"Filename does not match ONI timelapse-z grammar: {path.name}")

    return TimelapseFile(
        path=path,
        prefix=match.group("prefix"),
        channel=str(channel_name),
        pos_xy=int(match.group("pos_xy")),
        time_index=int(match.group("time_index")),
        z_index=int(match.group("z_index")),
    )


def find_timelapse_files(
    input_path: str | Path,
    glob_pattern: str = "*.tif",
    *,
    channel_name: str = "647",
    recursive: bool = True,
    exclude_name_contains: tuple[str, ...] = ("SACDpy",),
) -> list[TimelapseFile]:
    """Find and parse matching timelapse-z TIFF files from a file or folder."""

    path = Path(input_path)
    if path.is_file():
        candidates = [path]
    elif path.exists():
        globber = path.rglob if recursive else path.glob
        candidates = sorted(p for p in globber(glob_pattern) if p.is_file())
    else:
        raise FileNotFoundError(f"Input path does not exist: {path}")

    excluded = tuple(token.lower() for token in exclude_name_contains)
    parsed: list[TimelapseFile] = []
    for candidate in candidates:
        if any(token in candidate.name.lower() for token in excluded):
            continue
        try:
            parsed.append(parse_timelapse_filename(candidate, channel_name=channel_name))
        except ValueError:
            continue
    return parsed


def group_timelapse_files(files: Iterable[TimelapseFile]) -> list[TimelapseGroup]:
    grouped: dict[tuple[str, str, int], dict[tuple[int, int], Path]] = defaultdict(dict)
    for item in files:
        key = (item.prefix, item.channel, item.pos_xy)
        grouped[key][(item.time_index, item.z_index)] = item.path

    return [
        TimelapseGroup(prefix=prefix, channel=channel, pos_xy=pos_xy, files=files_by_tz)
        for (prefix, channel, pos_xy), files_by_tz in sorted(grouped.items())
    ]


def discover_folder_plans(
    input_folders: Iterable[str | Path],
    *,
    glob_pattern: str = "*.tif",
    channel_name: str = "647",
    recursive: bool = True,
    output_dir: str | Path | None = None,
) -> list[FolderTimelapsePlan]:
    """Discover timelapse-z groups for one or more FOV folders."""

    plans: list[FolderTimelapsePlan] = []
    for input_folder in input_folders:
        folder = Path(input_folder)
        files = find_timelapse_files(
            folder,
            glob_pattern,
            channel_name=channel_name,
            recursive=recursive,
        )
        groups = group_timelapse_files(files)
        folder_output_dir = Path(output_dir) if output_dir is not None else folder
        plans.append(
            FolderTimelapsePlan(
                input_folder=folder,
                output_dir=folder_output_dir,
                files=tuple(files),
                groups=tuple(groups),
            )
        )
    return plans


def read_oni_metadata(path: str | Path) -> dict:
    """Read the ONI JSON metadata stored in the first TIFF page description."""

    with tifffile.TiffFile(path) as tif:
        description = tif.pages[0].description
    if not description:
        return {}
    try:
        metadata = json.loads(description)
    except json.JSONDecodeError:
        return {}
    return metadata if isinstance(metadata, dict) else {}


def metadata_value(path: str | Path, key: str):
    return read_oni_metadata(path).get(key)


def infer_pixel_nm(path: str | Path, fallback_nm: float | None = None) -> float | None:
    pixel_size_um = metadata_value(path, "PixelSize_um")
    if pixel_size_um is None:
        return fallback_nm
    return float(pixel_size_um) * 1000.0


def infer_na(path: str | Path, fallback_na: float | None = None) -> float | None:
    objective_na = metadata_value(path, "Objective_NA")
    if objective_na is None:
        return fallback_na
    return float(objective_na)


def infer_z_spacing_um(group: TimelapseGroup) -> float | None:
    spacings: list[float] = []
    for time_index in group.time_indices:
        positions: list[float] = []
        for z_index in group.z_indices:
            path = group.files.get((time_index, z_index))
            if path is None:
                continue
            stage_pos = metadata_value(path, "StagePos_um")
            if isinstance(stage_pos, list) and len(stage_pos) >= 3:
                positions.append(float(stage_pos[2]))
        spacings.extend(abs(b - a) for a, b in zip(positions, positions[1:]))
    return float(median(spacings)) if spacings else None


def infer_time_interval_s(
    group: TimelapseGroup,
    *,
    source: str = "metadata",
    nominal_time_interval_s: float | None = None,
) -> float | None:
    """Infer the effective time spacing between reconstructed time points."""

    source = source.lower()
    if source == "nominal":
        return nominal_time_interval_s

    if source in {"metadata", "auto"}:
        interval = _metadata_time_interval_s(group)
        if interval is not None:
            return interval
        if nominal_time_interval_s is not None:
            return nominal_time_interval_s
        return parse_time_interval_from_path(next(iter(group.files.values())))

    if source == "folder":
        interval = parse_time_interval_from_path(next(iter(group.files.values())))
        return interval if interval is not None else nominal_time_interval_s

    raise ValueError("time interval source must be 'metadata', 'folder', 'nominal', or 'auto'.")


def parse_time_interval_from_path(path: str | Path) -> float | None:
    path = Path(path)
    for part in [path.stem, *(parent.name for parent in path.parents)]:
        match = _SECONDS_RE.search(part)
        if match is not None:
            return float(match.group("seconds"))
    return None


def output_paths_for_group(
    group: TimelapseGroup,
    output_dir: str | Path,
    *,
    output_label: str = "SACDpy",
) -> tuple[Path, Path]:
    output_dir = Path(output_dir)
    stem = f"{group.prefix}-{output_label}-{group.channel}-posXY{group.pos_xy}"
    return output_dir / f"{stem}-TZYX.tif", output_dir / f"{stem}-MIP-TYX.tif"


def write_timelapse_tiff(
    path: str | Path,
    image: np.ndarray,
    *,
    axes: str,
    pixel_size_um: float | None = None,
    z_spacing_um: float | None = None,
    time_interval_s: float | None = None,
) -> None:
    """Write assembled SACD timelapse outputs as float32 ImageJ TIFFs."""

    arr = np.asarray(image, dtype=np.float32)
    if arr.ndim != len(axes):
        raise ValueError(f"Image ndim {arr.ndim} does not match axes {axes!r}.")

    metadata: dict[str, object] = {"axes": axes}
    if pixel_size_um is not None:
        metadata["unit"] = "um"
    if "Z" in axes and z_spacing_um is not None:
        metadata["spacing"] = float(z_spacing_um)
    if "T" in axes and time_interval_s is not None:
        metadata["finterval"] = float(time_interval_s)

    kwargs: dict[str, object] = {"imagej": True, "metadata": metadata}
    if pixel_size_um is not None and pixel_size_um > 0:
        pixels_per_um = 1.0 / float(pixel_size_um)
        kwargs["resolution"] = (pixels_per_um, pixels_per_um)

    tifffile.imwrite(path, arr, **kwargs)


def _metadata_time_interval_s(group: TimelapseGroup) -> float | None:
    intervals: list[float] = []
    for z_index in group.z_indices:
        timestamps: list[tuple[int, int]] = []
        for time_index in group.time_indices:
            path = group.files.get((time_index, z_index))
            if path is None:
                continue
            timestamp_us = metadata_value(path, "timestamp_us")
            if timestamp_us is not None:
                timestamps.append((time_index, int(timestamp_us)))
        timestamps.sort()
        intervals.extend((b - a) / 1_000_000.0 for (_, a), (_, b) in zip(timestamps, timestamps[1:]))
    return float(median(intervals)) if intervals else None
