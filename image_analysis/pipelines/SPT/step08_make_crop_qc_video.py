"""Step 8: regenerate detailed crop QC videos from existing detections and trajectories."""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from tifffile import TiffFile, TiffWriter, imread

from spt_shared import load_config, output_dirs
from spt_video import render_overlay_video
from step01_inspect_inputs import build_manifest
from step02_preprocess_videos import artifact_paths, bandpass_frame, channel_specs
from step03_segment_cells import segmentation_paths
from step04_detect_spots import detection_paths
from step05_link_trajectories import tracking_paths


def _spec(cfg: dict, channel: str) -> dict[str, str]:
    matches = [spec for spec in channel_specs(cfg) if spec["name"] == channel and spec["role"] == "spt"]
    if not matches: raise ValueError(f"Unknown SPT channel {channel}")
    return matches[0]


def _expand_bbox(bbox: tuple[int, int, int, int], shape: tuple[int, int], padding: int) -> tuple[int, int, int, int]:
    x0, y0, x1, y1 = bbox; height, width = shape
    return max(0, x0-padding), max(0, y0-padding), min(width, x1+padding), min(height, y1+padding)


def select_crop(mask: np.ndarray, tracks: pd.DataFrame, *, cell_id: int | None, track_id: int | None,
                bbox: tuple[int, int, int, int] | None, padding: int) -> tuple[tuple[int, int, int, int], pd.DataFrame]:
    if cell_id is not None:
        yy, xx = np.where(mask == cell_id)
        if not len(xx): raise ValueError(f"cell_id {cell_id} is absent from the mask")
        crop = (int(xx.min()), int(yy.min()), int(xx.max())+1, int(yy.max())+1)
        selected = tracks[tracks["cell_id"] == cell_id]
    elif track_id is not None:
        selected = tracks[tracks["trackID"] == track_id]
        if selected.empty: raise ValueError(f"trackID {track_id} not found")
        crop = (int(np.floor(selected.x.min())), int(np.floor(selected.y.min())),
                int(np.ceil(selected.x.max()))+1, int(np.ceil(selected.y.max()))+1)
    elif bbox is not None:
        crop = bbox; x0, y0, x1, y1 = bbox
        selected_ids = tracks[(tracks.x >= x0) & (tracks.x < x1) & (tracks.y >= y0) & (tracks.y < y1)].trackID.unique()
        selected = tracks[tracks.trackID.isin(selected_ids)]
    else: raise ValueError("Select a cell, track, or bounding box")
    return _expand_bbox(crop, mask.shape, padding), selected


def _manual_dog(source: Path, destination: Path, cfg: dict) -> tuple[int, int]:
    low, high = float(cfg["bandpass"]["sigma_low"]), float(cfg["bandpass"]["sigma_high"])
    with TiffFile(source) as tif, TiffWriter(destination, bigtiff=True) as writer:
        shape = tif.pages[0].shape
        for page in tif.pages: writer.write(bandpass_frame(page.asarray(), low, high), contiguous=True, metadata=None)
    return shape


def make_crop_qc_video(config_path: str | Path, fov: str, channel: str, *, cell_id: int | None = None,
                       track_id: int | None = None, bbox: tuple[int, int, int, int] | None = None,
                       cropped_tif: str | Path | None = None, origin_x: int | None = None,
                       origin_y: int | None = None, start_frame: int = 0, end_frame: int | None = None,
                       maximum_frames: int = 100, output: str | Path | None = None) -> Path:
    cfg = load_config(config_path); manifest = build_manifest(cfg); rows = manifest[manifest.fov == fov]
    if rows.empty: raise ValueError(f"FOV {fov} not found in input manifest")
    row = rows.iloc[0]; spec = _spec(cfg, channel)
    mask = np.asarray(imread(segmentation_paths(cfg["output_dir"], fov)["mask"]))
    detections = pd.read_csv(detection_paths(cfg["output_dir"], fov, channel)["table"])
    tracks = pd.read_csv(tracking_paths(cfg["output_dir"], fov, channel)["canonical"])
    crop, selected_tracks = select_crop(mask, tracks, cell_id=cell_id, track_id=track_id, bbox=bbox,
                                        padding=int(cfg["qc"]["crop_padding_px"]))
    label = f"cell-{cell_id}" if cell_id is not None else f"track-{track_id}" if track_id is not None else "bbox"
    destination = Path(output) if output else output_dirs(cfg["output_dir"])["qc"]/"crop_videos"/f"{fov}__{channel}__{label}_QC.mp4"
    qc, spot = cfg["qc"], cfg["spotiflow"]
    if cropped_tif is None:
        dog_path = artifact_paths(cfg["output_dir"], fov, spec, cfg)["filtered"]; video_bbox = crop
        video_mask, video_detections, video_tracks = mask, detections, selected_tracks
        render_overlay_video(
            dog_path, destination, video_detections, tracks=video_tracks, mask=video_mask,
            frame_interval_s=float(row.frame_interval_s), pixel_size_um=float(row.pixel_size_um),
            model_label=f"Spotiflow {spot.get('pretrained_model', 'general')}",
            threshold_label="0.5 (model)" if spot.get("probability_threshold") is None else str(spot["probability_threshold"]),
            link_label="Existing LapTrack results: <5 px, gap ≤2 frames, ≥5 locations",
            maximum_frames=min(maximum_frames, int(qc["max_video_frames"])), scale=int(qc["render_scale"]),
            fps=int(qc["video_fps"]), crf=int(qc["video_crf"]), tail_frames=int(qc["track_tail_frames"]),
            scale_bar_um=float(qc["scale_bar_um"]), crop_bbox=video_bbox,
            start_frame=start_frame, end_frame=end_frame,
        )
    else:
        if origin_x is None or origin_y is None: raise ValueError("Manual cropped TIFF requires origin_x and origin_y")
        with tempfile.TemporaryDirectory(prefix="oni_spt_crop_qc_") as temporary:
            dog_path = Path(temporary)/"crop_DoG_float32.tif"; shape = _manual_dog(Path(cropped_tif), dog_path, cfg)
            transformed_detections = detections.copy(); transformed_tracks = selected_tracks.copy()
            for table in (transformed_detections, transformed_tracks):
                table["x"] -= origin_x; table["y"] -= origin_y
            cropped_mask = mask[origin_y:origin_y+shape[0], origin_x:origin_x+shape[1]]
            render_overlay_video(
                dog_path, destination, transformed_detections, tracks=transformed_tracks, mask=cropped_mask,
                frame_interval_s=float(row.frame_interval_s), pixel_size_um=float(row.pixel_size_um),
                model_label=f"Spotiflow {spot.get('pretrained_model', 'general')}",
                threshold_label="0.5 (model)" if spot.get("probability_threshold") is None else str(spot["probability_threshold"]),
                link_label="Existing full-FOV results", maximum_frames=min(maximum_frames, int(qc["max_video_frames"])),
                scale=int(qc["render_scale"]), fps=int(qc["video_fps"]), crf=int(qc["video_crf"]),
                tail_frames=int(qc["track_tail_frames"]), scale_bar_um=float(qc["scale_bar_um"]),
                start_frame=start_frame, end_frame=end_frame,
            )
    return destination


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", required=True)
    parser.add_argument("--fov", required=True); parser.add_argument("--channel", default="spt")
    selector = parser.add_mutually_exclusive_group(required=True); selector.add_argument("--cell-id", type=int)
    selector.add_argument("--track-id", type=int); selector.add_argument("--bbox", nargs=4, type=int, metavar=("X0", "Y0", "X1", "Y1"))
    parser.add_argument("--cropped-tif"); parser.add_argument("--origin-x", type=int); parser.add_argument("--origin-y", type=int)
    parser.add_argument("--start-frame", type=int, default=0); parser.add_argument("--end-frame", type=int)
    parser.add_argument("--max-frames", type=int, default=100); parser.add_argument("--output")
    args = parser.parse_args(argv)
    path = make_crop_qc_video(args.config, args.fov, args.channel, cell_id=args.cell_id, track_id=args.track_id,
                              bbox=tuple(args.bbox) if args.bbox else None, cropped_tif=args.cropped_tif,
                              origin_x=args.origin_x, origin_y=args.origin_y, start_frame=args.start_frame,
                              end_frame=args.end_frame, maximum_frames=args.max_frames, output=args.output)
    print(path)


if __name__ == "__main__":
    main()

