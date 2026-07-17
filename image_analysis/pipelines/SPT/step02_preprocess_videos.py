"""Step 2: split/copy ONI channels, save signed float32 DoG, and make raw MIPs."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from skimage.filters import difference_of_gaussians
from tifffile import TiffFile, TiffWriter, imwrite

from spt_shared import atomic_csv, atomic_json, output_action, output_dirs, stage_timer
from step01_inspect_inputs import accepted_rows, inspect_inputs


def channel_specs(cfg: dict[str, Any]) -> list[dict[str, str]]:
    return [
        {"position": position, "name": value.get("name", position), "role": value["role"]}
        for position, value in cfg["channels"].items() if value["role"] != "ignore"
    ]


def split_frame(frame: np.ndarray, mode: str) -> dict[str, np.ndarray]:
    if frame.ndim != 2:
        raise ValueError(f"Expected 2D ONI frame, got {frame.shape}")
    if mode == "single":
        return {"full": frame}
    if frame.shape[1] % 2:
        raise ValueError("Dual-channel width must be even")
    midpoint = frame.shape[1] // 2
    return {"left": frame[:, :midpoint], "right": frame[:, midpoint:]}


def bandpass_frame(frame: np.ndarray, sigma_low: float, sigma_high: float) -> np.ndarray:
    return difference_of_gaussians(frame.astype(np.float32), sigma_low, sigma_high).astype(np.float32)


def artifact_paths(output_dir: str | Path, fov: str, spec: dict[str, str], cfg: dict[str, Any]) -> dict[str, Path]:
    dirs = output_dirs(output_dir)
    name = spec["name"]
    paths = {
        "raw": dirs["raw"] / f"{fov}__{name}__raw.tif",
        "mip": dirs["mips"] / f"{fov}__{name}__raw_mip.tif",
        "metadata": dirs["raw"] / f"{fov}__{name}__metadata.json",
    }
    if spec["role"] == "spt":
        low, high = cfg["bandpass"]["sigma_low"], cfg["bandpass"]["sigma_high"]
        paths["filtered"] = dirs["bandpass"] / f"{fov}__{name}__DoG_sigma-{low:g}-{high:g}_float32.tif"
        paths["stats"] = dirs["bandpass"] / f"{fov}__{name}__DoG_statistics.csv"
    return paths


def _tmp(path: Path) -> Path:
    return path.with_name(f".{path.stem}.tmp{path.suffix}")


def preprocess_file(row: dict[str, Any], cfg: dict[str, Any], *, max_frames: int | None = None,
                    resume: bool = False, force: bool = False) -> dict[str, dict[str, str]]:
    specs = channel_specs(cfg)
    paths_by_name = {spec["name"]: artifact_paths(cfg["output_dir"], row["fov"], spec, cfg) for spec in specs}
    required = [path for paths in paths_by_name.values() for key, path in paths.items() if key != "metadata"]
    if output_action(required, resume=resume, force=force) == "skip":
        return {name: {key: str(path) for key, path in paths.items()} for name, paths in paths_by_name.items()}
    for path in required:
        path.parent.mkdir(parents=True, exist_ok=True)
    writers: dict[tuple[str, str], TiffWriter] = {}
    temporary: dict[tuple[str, str], Path] = {}
    mips: dict[str, np.ndarray] = {}
    stats: dict[str, dict[str, float]] = {}
    percentile_samples: dict[str, list[np.ndarray]] = {}
    try:
        with TiffFile(row["path"]) as tif:
            frame_count = len(tif.pages) if max_frames is None else min(max_frames, len(tif.pages))
            if frame_count < 1:
                raise ValueError("No frames selected")
            description = tif.pages[0].description
            try:
                oni_metadata = json.loads(description) if description else {}
            except json.JSONDecodeError:
                oni_metadata = {"ImageDescription": description}
            for spec in specs:
                paths = paths_by_name[spec["name"]]
                raw_tmp = _tmp(paths["raw"])
                writers[(spec["name"], "raw")] = TiffWriter(raw_tmp, bigtiff=True)
                temporary[(spec["name"], "raw")] = raw_tmp
                if spec["role"] == "spt":
                    dog_tmp = _tmp(paths["filtered"])
                    writers[(spec["name"], "filtered")] = TiffWriter(dog_tmp, bigtiff=True)
                    temporary[(spec["name"], "filtered")] = dog_tmp
                    stats[spec["name"]] = {"count": 0, "negative": 0, "zero": 0, "finite": 0, "minimum": np.inf, "maximum": -np.inf}
                    percentile_samples[spec["name"]] = []
            low, high = float(cfg["bandpass"]["sigma_low"]), float(cfg["bandpass"]["sigma_high"])
            for frame_index in range(frame_count):
                pieces = split_frame(tif.pages[frame_index].asarray(), cfg["layout"]["mode"])
                for spec in specs:
                    image = np.asarray(pieces[spec["position"]])
                    name = spec["name"]
                    mips[name] = image.copy() if name not in mips else np.maximum(mips[name], image)
                    writers[(name, "raw")].write(image, contiguous=True, description=description if frame_index == 0 else None, metadata=None)
                    if spec["role"] == "spt":
                        dog = bandpass_frame(image, low, high)
                        writers[(name, "filtered")].write(dog, contiguous=True, description=description if frame_index == 0 else None, metadata=None)
                        current = stats[name]
                        current["count"] += dog.size
                        current["negative"] += int(np.count_nonzero(dog < 0))
                        current["zero"] += int(np.count_nonzero(dog == 0))
                        current["finite"] += int(np.count_nonzero(np.isfinite(dog)))
                        current["minimum"] = min(current["minimum"], float(np.nanmin(dog)))
                        current["maximum"] = max(current["maximum"], float(np.nanmax(dog)))
                        # A deterministic bounded sample gives useful run-level
                        # percentiles without ever materializing the stack.
                        stride = max(1, dog.size // 20_000)
                        percentile_samples[name].append(dog.ravel()[::stride][:20_000])
            for writer in writers.values():
                writer.close()
            writers.clear()
            for (name, kind), temp_path in temporary.items():
                os.replace(temp_path, paths_by_name[name][kind])
            for spec in specs:
                name, paths = spec["name"], paths_by_name[spec["name"]]
                imwrite(paths["mip"], mips[name], metadata=None)
                atomic_json({
                    "source_path": row["path"], "dataset": cfg["dataset"], "fov": row["fov"],
                    "condition": row.get("condition"), "channel": name, "role": spec["role"],
                    "frames_written": frame_count, "source_actual_frames": len(tif.pages),
                    "partial": frame_count != len(tif.pages), "oni_metadata": oni_metadata,
                }, paths["metadata"])
                if spec["role"] == "spt":
                    current = stats[name]
                    count = current.pop("count")
                    sampled = np.concatenate(percentile_samples[name])
                    sampled = sampled[np.isfinite(sampled)]
                    percentiles = np.percentile(sampled, [0.1, 1, 50, 99, 99.8, 99.9])
                    atomic_csv(pd.DataFrame([{
                        "dataset": cfg["dataset"], "fov": row["fov"], "channel": name,
                        "dtype": "float32", "minimum": current["minimum"], "maximum": current["maximum"],
                        "negative_fraction": current["negative"] / count,
                        "zero_fraction": current["zero"] / count,
                        "finite_fraction": current["finite"] / count,
                        "percentile_0_1": percentiles[0], "percentile_1": percentiles[1],
                        "percentile_50": percentiles[2], "percentile_99": percentiles[3],
                        "percentile_99_8": percentiles[4], "percentile_99_9": percentiles[5],
                        "sigma_low": low, "sigma_high": high,
                    }]), paths["stats"])
    finally:
        for writer in writers.values():
            writer.close()
        for path in temporary.values():
            if path.exists():
                path.unlink()
    return {name: {key: str(path) for key, path in paths.items()} for name, paths in paths_by_name.items()}


def preprocess_stage(cfg: dict[str, Any], manifest: pd.DataFrame, *, max_files: int | None = None,
                     max_frames: int | None = None, resume: bool = False, force: bool = False) -> None:
    rows = accepted_rows(manifest, max_files)
    if not rows:
        raise ValueError("No accepted TIFFs")
    with stage_timer(cfg, "02_preprocess", {"files": len(rows), "max_frames": max_frames}):
        for index, row in enumerate(rows, 1):
            print(f"[preprocess {index}/{len(rows)}] {row['filename']}", flush=True)
            preprocess_file(row, cfg, max_frames=max_frames, resume=resume, force=force)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True); parser.add_argument("--max-files", type=int)
    parser.add_argument("--max-frames", type=int)
    policy = parser.add_mutually_exclusive_group(); policy.add_argument("--resume", action="store_true"); policy.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    cfg, manifest = inspect_inputs(args.config)
    preprocess_stage(cfg, manifest, max_files=args.max_files, max_frames=args.max_frames, resume=args.resume, force=args.force)


if __name__ == "__main__":
    main()
