"""Read-only acceptance audit for a completed production ONI SPT batch.

The optional JSON reports are deliberately stored beside each dataset's resolved
configuration so a future reader can see exactly what was checked after the run.
"""

from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path
import sys
from typing import Any

import imageio_ffmpeg
import numpy as np
import pandas as pd
from tifffile import TiffFile
import yaml

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from run_ONI_SPT_batch import load_batch
from spt_shared import atomic_json, load_config, output_dirs


REQUIRED_STAGES = {
    "02_preprocess", "03_segment", "04_detect", "05_track",
    "06_diffusion", "07_report",
}
DETECTION_COLUMNS = {
    "dataset", "condition", "fov", "channel", "frame", "x", "y",
    "spotiflow_probability", "raw_mean_intensity", "signed_dog_mean",
    "detection_cell_id",
}
TRAJECTORY_COLUMNS = {
    "trackID", "frame", "x", "y", "cell_id", "assignment_confidence",
    "assignment_ambiguous", "trajectory_uid",
}
LEGACY_COLUMNS = {
    "spotID", "trackID", "x", "y", "t", "R", "meanIntensity",
    "medianIntensity", "minIntensity", "maxIntensity", "totalIntensity",
    "stdIntensity", "contrast", "SNR",
}


def _assert_columns(path: Path, required: set[str]) -> pd.DataFrame:
    table = pd.read_csv(path, low_memory=False)
    missing = required - set(table.columns)
    assert not missing, f"{path}: missing columns {sorted(missing)}"
    return table


def _validate_tiff(path: Path, frames: int, height: int, width: int, dtype: str) -> dict[str, Any]:
    with TiffFile(path) as tif:
        assert len(tif.pages) == frames, f"{path}: {len(tif.pages)} pages != {frames}"
        assert tif.pages[0].shape == (height, width), f"{path}: wrong page shape"
        assert tif.pages[0].dtype == np.dtype(dtype), f"{path}: wrong dtype"
        indices = sorted({0, frames // 2, frames - 1})
        samples = np.concatenate([tif.pages[index].asarray().ravel()[::97] for index in indices])
    assert np.isfinite(samples).all(), f"{path}: non-finite TIFF samples"
    return {"minimum": float(samples.min()), "maximum": float(samples.max())}


def _validate_dataset(cfg: dict[str, Any], expected: dict[str, Any], check_videos: bool) -> dict[str, Any]:
    dirs = output_dirs(cfg["output_dir"])
    metadata = dirs["metadata"]
    for name in (
        "pipeline_config.yaml", "resolved_config.yaml", "batch_config.yaml",
        "provenance.json", "run_command.txt", "input_manifest.csv",
        "unused_files.csv", "stage_status.json",
    ):
        assert (metadata / name).exists(), f"{cfg['dataset']}: missing {name}"
    assert len(list((metadata / "config_sources").glob("*.yaml"))) >= 2

    with (metadata / "resolved_config.yaml").open() as handle:
        resolved = yaml.safe_load(handle)
    assert resolved["spotiflow"]["pretrained_model"] == "general"
    assert float(resolved["spotiflow"]["probability_threshold"]) == 0.4
    assert resolved["laptrack"]["metric"] == "sqeuclidean"
    assert float(resolved["laptrack"]["max_link_distance_px"]) == 5
    assert bool(resolved["laptrack"]["gap_closing"])
    assert int(resolved["laptrack"]["gap_closing_max_frame_count"]) == 2
    assert int(resolved["laptrack"]["min_track_length"]) == 5

    with (metadata / "stage_status.json").open() as handle:
        stages = json.load(handle)
    assert REQUIRED_STAGES <= set(stages)
    assert all(stages[name]["status"] == "complete" for name in REQUIRED_STAGES)

    manifest = pd.read_csv(metadata / "input_manifest.csv")
    accepted = manifest[manifest.status == "accepted"].copy()
    unused = pd.read_csv(metadata / "unused_files.csv")
    assert len(accepted) == int(expected["expected_accepted_files"])
    assert int(accepted.actual_frames.sum()) == int(expected["expected_accepted_frames"])
    assert len(unused) == int(expected["expected_skipped_files"])

    output_channels = [value["name"] for value in cfg["channels"].values() if value["role"] != "ignore"]
    spt_channels = [value["name"] for value in cfg["channels"].values() if value["role"] == "spt"]
    marker_channels = [value["name"] for value in cfg["channels"].values() if value["role"] == "marker"]
    expected_spt_outputs = len(accepted) * len(spt_channels)
    expected_raw_outputs = len(accepted) * len(output_channels)
    assert len(list(dirs["raw"].glob("*__raw.tif"))) == expected_raw_outputs
    assert len(list(dirs["bandpass"].glob("*__DoG_*_float32.tif"))) == expected_spt_outputs
    assert len(list(dirs["detections"].glob("*__detections.csv"))) == expected_spt_outputs
    assert len(list(dirs["trajectories"].glob("*__trajectories.csv"))) == expected_spt_outputs
    for marker in marker_channels:
        assert not list(dirs["detections"].glob(f"*__{marker}__*")), f"detections found on marker {marker}"

    dog_minimum = float("inf")
    dog_maximum = float("-inf")
    for row in accepted.itertuples(index=False):
        for channel in output_channels:
            raw_path = dirs["raw"] / f"{row.fov}__{channel}__raw.tif"
            raw_check = _validate_tiff(raw_path, int(row.actual_frames), int(row.height), 428, "uint16")
            assert raw_check["maximum"] > raw_check["minimum"]
        for channel in spt_channels:
            dog_path = dirs["bandpass"] / f"{row.fov}__{channel}__DoG_sigma-1-3_float32.tif"
            dog_check = _validate_tiff(dog_path, int(row.actual_frames), int(row.height), 428, "float32")
            assert dog_check["minimum"] < 0 < dog_check["maximum"], f"{dog_path}: signed range lost"
            dog_minimum = min(dog_minimum, dog_check["minimum"])
            dog_maximum = max(dog_maximum, dog_check["maximum"])
            stats = pd.read_csv(dirs["bandpass"] / f"{row.fov}__{channel}__DoG_statistics.csv")
            assert set(["minimum", "maximum", "negative_fraction", "finite_fraction", "dtype"]) <= set(stats.columns)
            assert (stats["dtype"] == "float32").all()
            assert (stats["minimum"] < 0).all() and (stats["maximum"] > 0).all()
            assert np.allclose(stats["finite_fraction"], 1.0)

    detections = retained_tracks = trajectory_points = gap_edges = 0
    assigned_points = 0
    for path in sorted(dirs["detections"].glob("*__detections.csv")):
        table = _assert_columns(path, DETECTION_COLUMNS)
        assert set(table.channel.astype(str)) <= set(spt_channels)
        assert table.frame.min() >= 0
        assert table.spotiflow_probability.min() >= 0.4 - 1e-7
        detections += len(table)
    for path in sorted(dirs["trajectories"].glob("*__trajectories.csv")):
        table = _assert_columns(path, TRAJECTORY_COLUMNS)
        if not table.empty:
            assert np.allclose(table.cell_id, table.cell_id.astype(int))
            assert table.assignment_confidence.between(0, 1).all()
            lengths = table.groupby("trackID").size()
            assert int(lengths.min()) >= 5
            retained_tracks += len(lengths)
            trajectory_points += len(table)
            assigned_points += int((table.cell_id > 0).sum())
        edge_path = path.with_name(path.name.replace("__trajectories.csv", "__trajectory_edges.csv"))
        edges = _assert_columns(edge_path, {"distance_px", "frame_difference", "edge_type"})
        if not edges.empty:
            assert edges.distance_px.max() <= 5 + 1e-9
            assert edges.frame_difference.min() >= 1 and edges.frame_difference.max() <= 2
            assert set(edges.edge_type) <= {"adjacent", "gap_closing"}
            gap_edges += int((edges.edge_type == "gap_closing").sum())
        legacy_path = path.with_name(path.name.replace("__trajectories.csv", "__tracks_TrackMate_compatible.csv"))
        legacy = _assert_columns(legacy_path, LEGACY_COLUMNS)
        assert all(np.issubdtype(legacy[name].dtype, np.number) for name in LEGACY_COLUMNS)

    per_fov = list((dirs["analysis"] / "per_fov").glob("*.csv"))
    pooled = list((dirs["analysis"] / "by_condition").glob("*.csv"))
    saspt = list((dirs["analysis"] / "saspt").glob("*.csv"))
    assert len(per_fov) == expected_spt_outputs
    assert pooled and saspt
    for path in per_fov:
        assert {"cell_id", "fov", "channel", "trajectory_uid"} <= set(pd.read_csv(path, nrows=0).columns)
    assert (dirs["reports"] / "condition_summary.csv").exists()
    assert (dirs["reports"] / "diffusion_comparison_report.pdf").exists()
    assert (dirs["reports"] / "replicate_cell_counts.csv").exists()
    assert len(list((dirs["reports"] / "plots").glob("*.png"))) == 9
    assert len(list((dirs["qc"] / "trajectory_overlays").glob("*1200dpi.png"))) == expected_spt_outputs
    assert len(list((dirs["qc"] / "trajectory_overlays").glob("*vector.pdf"))) == expected_spt_outputs
    assert len(list((dirs["qc"] / "secondary_metrics").glob("*.png"))) == expected_spt_outputs

    video_count = 0
    if check_videos:
        videos = sorted((dirs["qc"] / "detection_videos").glob("*.mp4"))
        videos += sorted((dirs["qc"] / "trajectory_videos").glob("*.mp4"))
        assert len(videos) == expected_spt_outputs * 2
        for path in videos:
            frames, _ = imageio_ffmpeg.count_frames_and_secs(path)
            assert 1 <= frames <= 100, f"{path}: {frames} frames"
            video_count += 1

    return {
        "dataset": cfg["dataset"],
        "validated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "accepted_files": len(accepted),
        "accepted_frames": int(accepted.actual_frames.sum()),
        "skipped_files": len(unused),
        "detections": detections,
        "retained_tracks": retained_tracks,
        "trajectory_points": trajectory_points,
        "assigned_point_fraction": assigned_points / trajectory_points if trajectory_points else 0.0,
        "gap_closing_edges": gap_edges,
        "dog_sampled_minimum": dog_minimum,
        "dog_sampled_maximum": dog_maximum,
        "per_fov_aio_tables": len(per_fov),
        "pooled_condition_tables": len(pooled),
        "saspt_tables": len(saspt),
        "qc_videos_checked": video_count,
        "status": "passed",
    }


def validate_batch(batch_path: str | Path, *, check_videos: bool, write_reports: bool) -> dict[str, Any]:
    batch = load_batch(batch_path)
    status_path = Path(batch["base_output_dir"]) / "_batch_status" / f"{batch['name']}.json"
    with status_path.open() as handle:
        status = json.load(handle)
    assert status["status"] == "complete" and not status.get("failures")

    results = []
    for expected, config_path in zip(batch["datasets"], batch["_configs"]):
        cfg = load_config(config_path)
        result = _validate_dataset(cfg, expected, check_videos)
        results.append(result)
        if write_reports:
            atomic_json(result, output_dirs(cfg["output_dir"])["metadata"] / "production_validation.json")
        print(
            f"PASS {cfg['dataset']}: {result['accepted_files']} files, "
            f"{result['detections']:,} detections, {result['retained_tracks']:,} tracks",
            flush=True,
        )
    summary = {
        "batch": batch["name"],
        "validated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "status": "passed",
        "datasets": results,
        "totals": {
            "accepted_files": sum(value["accepted_files"] for value in results),
            "accepted_frames": sum(value["accepted_frames"] for value in results),
            "skipped_files": sum(value["skipped_files"] for value in results),
            "detections": sum(value["detections"] for value in results),
            "retained_tracks": sum(value["retained_tracks"] for value in results),
            "trajectory_points": sum(value["trajectory_points"] for value in results),
            "gap_closing_edges": sum(value["gap_closing_edges"] for value in results),
        },
    }
    if write_reports:
        atomic_json(summary, status_path.with_name(f"{batch['name']}__production_validation.json"))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-config", required=True)
    parser.add_argument("--skip-video-check", action="store_true")
    parser.add_argument("--write-reports", action="store_true")
    args = parser.parse_args()
    result = validate_batch(
        args.batch_config,
        check_videos=not args.skip_video_check,
        write_reports=args.write_reports,
    )
    print(json.dumps(result["totals"], indent=2))


if __name__ == "__main__":
    main()
