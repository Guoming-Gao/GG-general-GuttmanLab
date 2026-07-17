from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import numpy as np
from skimage.filters import difference_of_gaussians
from skimage.util import img_as_uint
from tifffile import TiffFile, TiffWriter, imwrite

from .utils import atomic_json, output_action


def channel_specs(cfg: dict[str, Any]) -> list[dict[str, str]]:
    specs = []
    for position, value in cfg["channels"].items():
        if value["role"] == "ignore":
            continue
        specs.append(
            {
                "position": position,
                "name": value.get("name", position),
                "role": value["role"],
            }
        )
    return specs


def split_frame(frame: np.ndarray, mode: str) -> dict[str, np.ndarray]:
    if frame.ndim != 2:
        raise ValueError(f"Expected 2D ONI frame, received shape {frame.shape}")
    if mode == "single":
        return {"full": frame}
    if frame.shape[1] % 2:
        raise ValueError(f"Dual-channel frame width must be even, got {frame.shape[1]}")
    midpoint = frame.shape[1] // 2
    return {"left": frame[:, :midpoint], "right": frame[:, midpoint:]}


def bandpass_frame(frame: np.ndarray, sigma_low: float, sigma_high: float) -> np.ndarray:
    filtered = difference_of_gaussians(frame, sigma_low, sigma_high)
    return img_as_uint(filtered)


def artifact_paths(output_dir: str | Path, fov: str, spec: dict[str, str], cfg: dict[str, Any]) -> dict[str, Path]:
    out = Path(output_dir)
    name = spec["name"]
    low = cfg["bandpass"]["sigma_low"]
    high = cfg["bandpass"]["sigma_high"]
    paths = {
        "raw": out / "channels" / f"{fov}__{name}__raw.tif",
        "mip": out / "mips" / f"{fov}__{name}__raw_mip.tif",
        "metadata": out / "channels" / f"{fov}__{name}__metadata.json",
    }
    if spec["role"] == "spt":
        paths["filtered"] = (
            out
            / "bandpass"
            / f"{fov}__{name}__dog_sigma1-{low:g}_sigma2-{high:g}.tif"
        )
    return paths


def _tmp_tiff(path: Path) -> Path:
    return path.with_name(f".{path.stem}.tmp{path.suffix}")


def preprocess_file(
    row: dict[str, Any],
    cfg: dict[str, Any],
    *,
    max_frames: int | None = None,
    resume: bool = False,
    force: bool = False,
) -> dict[str, dict[str, str]]:
    specs = channel_specs(cfg)
    all_paths = {
        spec["name"]: artifact_paths(cfg["output_dir"], row["fov"], spec, cfg)
        for spec in specs
    }
    required = [path for paths in all_paths.values() for key, path in paths.items() if key != "metadata"]
    if output_action(required, resume=resume, force=force) == "skip":
        return {
            name: {key: str(path) for key, path in paths.items()}
            for name, paths in all_paths.items()
        }
    for path in required:
        path.parent.mkdir(parents=True, exist_ok=True)
    writers: dict[tuple[str, str], TiffWriter] = {}
    tmp_paths: dict[tuple[str, str], Path] = {}
    mips: dict[str, np.ndarray] = {}
    source = Path(row["path"])
    try:
        with TiffFile(source) as tif:
            n_frames = len(tif.pages) if max_frames is None else min(len(tif.pages), max_frames)
            if n_frames <= 0:
                raise ValueError("max_frames selected zero frames")
            description = tif.pages[0].description
            try:
                oni_metadata = json.loads(description) if description else {}
            except json.JSONDecodeError:
                oni_metadata = {"ImageDescription": description}
            for spec in specs:
                paths = all_paths[spec["name"]]
                if cfg["save"]["raw_channels"]:
                    tmp = _tmp_tiff(paths["raw"])
                    writers[(spec["name"], "raw")] = TiffWriter(tmp, bigtiff=True)
                    tmp_paths[(spec["name"], "raw")] = tmp
                if spec["role"] == "spt" and cfg["save"]["filtered_spt"]:
                    tmp = _tmp_tiff(paths["filtered"])
                    writers[(spec["name"], "filtered")] = TiffWriter(tmp, bigtiff=True)
                    tmp_paths[(spec["name"], "filtered")] = tmp
            low = float(cfg["bandpass"]["sigma_low"])
            high = float(cfg["bandpass"]["sigma_high"])
            for frame_index in range(n_frames):
                frame = tif.pages[frame_index].asarray()
                split = split_frame(frame, cfg["layout"]["mode"])
                for spec in specs:
                    image = np.asarray(split[spec["position"]])
                    name = spec["name"]
                    if name not in mips:
                        mips[name] = image.copy()
                    else:
                        np.maximum(mips[name], image, out=mips[name])
                    raw_writer = writers.get((name, "raw"))
                    if raw_writer is not None:
                        raw_writer.write(
                            image,
                            contiguous=True,
                            description=description if frame_index == 0 else None,
                            metadata=None,
                        )
                    filtered_writer = writers.get((name, "filtered"))
                    if filtered_writer is not None:
                        filtered_writer.write(
                            bandpass_frame(image, low, high),
                            contiguous=True,
                            description=description if frame_index == 0 else None,
                            metadata=None,
                        )
            for writer in writers.values():
                writer.close()
            writers.clear()
            for (name, kind), tmp in tmp_paths.items():
                os.replace(tmp, all_paths[name][kind])
            for spec in specs:
                paths = all_paths[spec["name"]]
                imwrite(paths["mip"], mips[spec["name"]], metadata=None)
                atomic_json(
                    {
                        "source_path": str(source),
                        "dataset": cfg["dataset"],
                        "fov": row["fov"],
                        "channel": spec["name"],
                        "role": spec["role"],
                        "frames_written": n_frames,
                        "source_actual_frames": len(tif.pages),
                        "partial": n_frames != len(tif.pages),
                        "oni_metadata": oni_metadata,
                    },
                    paths["metadata"],
                )
    finally:
        for writer in writers.values():
            writer.close()
        for tmp in tmp_paths.values():
            if tmp.exists():
                tmp.unlink()
    return {
        name: {key: str(path) for key, path in paths.items()}
        for name, paths in all_paths.items()
    }
