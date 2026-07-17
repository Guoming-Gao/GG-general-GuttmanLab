"""Read-only acceptance audit for a completed revised ONI SPT pilot."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import imageio_ffmpeg
import numpy as np
import pandas as pd
from tifffile import TiffFile

from spt_shared import load_config, output_dirs


def validate(config_path: str | Path, expected_skipped: int) -> dict[str, object]:
    cfg = load_config(config_path)
    dirs = output_dirs(cfg["output_dir"])
    manifest = pd.read_csv(dirs["metadata"] / "input_manifest.csv")
    unused = pd.read_csv(dirs["metadata"] / "unused_files.csv")
    assert len(unused) == expected_skipped
    accepted = manifest[manifest.status == "accepted"]
    assert len(accepted) >= 1

    dogs = sorted(dirs["bandpass"].glob("*__DoG_*_float32.tif"))
    assert dogs
    dog_checks = []
    for path in dogs:
        with TiffFile(path) as tif:
            assert len(tif.pages) == 50
            assert tif.pages[0].shape == (684, 428)
            assert tif.pages[0].dtype == np.dtype("float32")
            values = np.concatenate([page.asarray().ravel()[::100] for page in tif.pages])
        assert values.min() < 0 < values.max()
        dog_checks.append({"file": path.name, "minimum": float(values.min()), "maximum": float(values.max())})

    detections = sorted(dirs["detections"].glob("*__detections.csv"))
    assert detections and not list(dirs["detections"].glob("*hoechst*"))
    detection_count = 0
    for path in detections:
        table = pd.read_csv(path)
        detection_count += len(table)
        assert set(["spotiflow_probability", "signed_dog_mean", "raw_mean_intensity", "detection_cell_id"]) <= set(table)

    tracks = sorted(dirs["trajectories"].glob("*__trajectories.csv"))
    assert tracks
    track_count = gap_edges = 0
    for path in tracks:
        table = pd.read_csv(path)
        assert set(["cell_id", "assignment_confidence", "trajectory_uid"]) <= set(table)
        assert np.allclose(table["cell_id"], table["cell_id"].astype(int))
        lengths = table.groupby("trackID").size()
        assert lengths.empty or int(lengths.min()) >= 5
        track_count += len(lengths)
        edge_path = path.with_name(path.name.replace("__trajectories.csv", "__trajectory_edges.csv"))
        edges = pd.read_csv(edge_path)
        if not edges.empty:
            assert edges["distance_px"].max() < 5 + 1e-9
            assert edges["frame_difference"].max() <= 2
            assert set(edges["edge_type"]) <= {"adjacent", "gap_closing"}
            gap_edges += int((edges.edge_type == "gap_closing").sum())

    videos = sorted((dirs["qc"] / "detection_videos").glob("*.mp4"))
    videos += sorted((dirs["qc"] / "trajectory_videos").glob("*.mp4"))
    assert len(videos) >= 2
    video_checks = []
    for path in videos:
        frames, seconds = imageio_ffmpeg.count_frames_and_secs(path)
        reader = imageio_ffmpeg.read_frames(path); metadata = next(reader); reader.close()
        assert frames == 50
        assert metadata["size"] == (856, 1508)
        video_checks.append({"file": path.name, "frames": frames, "size": metadata["size"], "seconds": seconds})
    for path in sorted((dirs["qc"] / "crop_videos").glob("*.mp4")):
        frames, _ = imageio_ffmpeg.count_frames_and_secs(path)
        assert frames <= 100

    assert list((dirs["analysis"] / "per_fov").glob("*.csv"))
    assert list((dirs["analysis"] / "by_condition").glob("*.csv"))
    assert list((dirs["analysis"] / "saspt").glob("*.csv"))
    assert (dirs["reports"] / "diffusion_comparison_report.pdf").exists()
    assert len(list((dirs["reports"] / "plots").glob("*.png"))) == 9
    assert (dirs["reports"] / "replicate_cell_counts.csv").exists()
    assert list((dirs["qc"] / "trajectory_overlays").glob("*1200dpi.png"))
    assert list((dirs["qc"] / "trajectory_overlays").glob("*vector.pdf"))
    assert list((dirs["qc"] / "secondary_metrics").glob("*.png"))
    return {
        "dataset": cfg["dataset"], "accepted": len(accepted), "skipped": len(unused),
        "detections": detection_count, "retained_tracks": track_count, "gap_edges": gap_edges,
        "dog": dog_checks, "videos": video_checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--expected-skipped", type=int, default=0)
    args = parser.parse_args()
    print(json.dumps(validate(args.config, args.expected_skipped), indent=2))


if __name__ == "__main__":
    main()
