"""Loss-aware MP4 rendering helpers shared by detection, tracking, and crop QC."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from skimage.segmentation import find_boundaries
from tifffile import TiffFile


def selected_frames(total: int, maximum: int, start: int = 0, end: int | None = None) -> np.ndarray:
    stop = total if end is None else min(total, end + 1)
    values = np.arange(max(0, start), stop, dtype=int)
    if len(values) <= maximum:
        return values
    return np.unique(np.rint(np.linspace(values[0], values[-1], maximum)).astype(int))


def _font(size: int) -> ImageFont.ImageFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for candidate in candidates:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size)
    return ImageFont.load_default()


def _contrast(path: Path, indices: Iterable[int]) -> tuple[float, float]:
    samples = []
    with TiffFile(path) as tif:
        for index in indices:
            frame = tif.pages[int(index)].asarray().astype(np.float32)
            positive = frame[np.isfinite(frame) & (frame > 0)]
            if positive.size:
                stride = max(1, positive.size // 20000)
                samples.append(positive[::stride])
    if not samples:
        return 0.0, 1.0
    values = np.concatenate(samples)
    upper = float(np.percentile(values, 99.8))
    return 0.0, upper if upper > 0 else 1.0


def _base_canvas(frame: np.ndarray, *, scale: int, lower: float, upper: float,
                 header_px: int, footer_px: int, canvas_width: int
                 ) -> tuple[Image.Image, ImageDraw.ImageDraw, int, int]:
    gray = np.clip((frame.astype(np.float32) - lower) / (upper - lower), 0, 1)
    image = Image.fromarray(np.uint8(gray * 255), mode="L").resize(
        (frame.shape[1] * scale, frame.shape[0] * scale), Image.Resampling.BILINEAR
    ).convert("RGB")
    image_left = (canvas_width - image.width) // 2
    canvas = Image.new("RGB", (canvas_width, image.height + header_px + footer_px), "black")
    canvas.paste(image, (image_left, header_px))
    return canvas, ImageDraw.Draw(canvas), image_left, header_px


def _draw_mask(draw: ImageDraw.ImageDraw, boundary_yx: np.ndarray, scale: int, y_offset: int,
               crop_bbox: tuple[int, int, int, int], image_left: int) -> None:
    x0, y0, x1, y1 = crop_bbox
    for y, x in boundary_yx[::2]:
        if not (x0 <= x < x1 and y0 <= y < y1):
            continue
        xx, yy = (int(x) - x0) * scale + image_left, (int(y) - y0) * scale + y_offset
        draw.point((xx, yy), fill=(0, 255, 255))


def _draw_time_and_scale(draw: ImageDraw.ImageDraw, width: int, image_bottom: int,
                         frame: int, last_frame: int, frame_interval_s: float,
                         pixel_size_um: float, scale: int, scale_bar_um: float,
                         font: ImageFont.ImageFont, image_left: int, image_width: int) -> None:
    bar_px = max(1, round(scale_bar_um / pixel_size_um * scale))
    x2, y = image_left + image_width - 16, image_bottom - 18
    x1 = x2 - bar_px
    draw.line((x1, y, x2, y), fill="white", width=max(2, scale))
    draw.text((x1, y - 22), f"{scale_bar_um:g} µm", fill="white", font=font)
    progress_left, progress_right = 12, width - 12
    progress_y = image_bottom + 13
    draw.line((progress_left, progress_y, progress_right, progress_y), fill=(90, 90, 90), width=4)
    fraction = 0 if last_frame <= 0 else frame / last_frame
    draw.line((progress_left, progress_y, progress_left + int((progress_right-progress_left)*fraction), progress_y), fill="white", width=4)
    draw.text((12, image_bottom + 18), f"t = {frame * frame_interval_s:.3f} s    Δt = {frame_interval_s:.3f} s", fill="white", font=font)


def _writer(path: Path, size: tuple[int, int], fps: int, crf: int):
    try:
        import imageio_ffmpeg
    except ImportError as exc:
        raise RuntimeError("MP4 QC requires imageio-ffmpeg; install the environment specification") from exc
    path.parent.mkdir(parents=True, exist_ok=True)
    generator = imageio_ffmpeg.write_frames(
        str(path), size, fps=fps, codec="libx264", pix_fmt_in="rgb24", pix_fmt_out="yuv420p",
        macro_block_size=1,
        output_params=["-crf", str(crf), "-preset", "slow", "-movflags", "+faststart"],
    )
    generator.send(None)
    return generator


def render_overlay_video(
    bandpass_path: str | Path,
    output_path: str | Path,
    detections: pd.DataFrame,
    *,
    tracks: pd.DataFrame | None,
    mask: np.ndarray | None,
    frame_interval_s: float,
    pixel_size_um: float,
    model_label: str,
    threshold_label: str,
    link_label: str = "",
    maximum_frames: int = 100,
    scale: int = 2,
    fps: int = 10,
    crf: int = 10,
    tail_frames: int = 20,
    scale_bar_um: float = 1.0,
    crop_bbox: tuple[int, int, int, int] | None = None,
    start_frame: int = 0,
    end_frame: int | None = None,
) -> list[int]:
    bandpass_path, output_path = Path(bandpass_path), Path(output_path)
    with TiffFile(bandpass_path) as tif:
        total = len(tif.pages); shape = tif.pages[0].shape
    indices = selected_frames(total, maximum_frames, start_frame, end_frame)
    lower, upper = _contrast(bandpass_path, indices)
    x0, y0, x1, y1 = crop_bbox or (0, 0, shape[1], shape[0])
    image_width, height = (x1 - x0) * scale, (y1 - y0) * scale
    width = max(856, image_width)
    header, footer = 42 * scale, 28 * scale
    font = _font(11 * scale); small = _font(9 * scale)
    boundary = np.argwhere(find_boundaries(mask, mode="outer")) if mask is not None else np.empty((0, 2), int)
    writer = _writer(output_path, (width, height + header + footer), fps, crf)
    colors = [(255, 70, 70), (80, 220, 80), (80, 150, 255), (255, 180, 40), (210, 80, 255)]
    try:
        with TiffFile(bandpass_path) as tif:
            for frame_index in indices:
                frame = tif.pages[int(frame_index)].asarray()[y0:y1, x0:x1]
                canvas, draw, image_left, y_offset = _base_canvas(
                    frame, scale=scale, lower=lower, upper=upper,
                    header_px=header, footer_px=footer, canvas_width=width,
                )
                current = detections[
                    (detections["frame"] == frame_index)
                    & (detections["x"] >= x0) & (detections["x"] < x1)
                    & (detections["y"] >= y0) & (detections["y"] < y1)
                ]
                draw.text((10, 6), f"Frame {frame_index}/{total-1}  |  detections: {len(current)}  |  {model_label} threshold={threshold_label}", fill="white", font=font)
                if link_label:
                    draw.text((10, 8 + 17 * scale), link_label, fill=(200, 200, 200), font=small)
                _draw_mask(draw, boundary, scale, y_offset, (x0, y0, x1, y1), image_left)
                radius = 2.5 * scale
                for row in current.itertuples():
                    xx = (float(row.x) - x0) * scale + image_left
                    yy = (float(row.y) - y0) * scale + y_offset
                    draw.ellipse((xx-radius, yy-radius, xx+radius, yy+radius), outline=(255, 255, 0), width=max(1, scale))
                if tracks is not None and not tracks.empty:
                    active_ids = sorted(
                        tracks.loc[tracks["frame"] == frame_index, "trackID"].astype(int).unique().tolist()
                    )
                    active = tracks[
                        tracks["trackID"].isin(active_ids)
                        & (tracks["frame"] <= frame_index)
                        & (tracks["frame"] >= frame_index-tail_frames)
                    ]
                    if active_ids:
                        shown = active_ids[:5]
                        ids = ",".join(map(str, shown))
                        if len(active_ids) > len(shown):
                            ids += f" +{len(active_ids)-len(shown)}"
                        draw.text((width - 10, 8 + 17 * scale), f"active tracks: {ids}",
                                  anchor="ra", fill=(200, 200, 200), font=small)
                    for track_id, group in active.groupby("trackID"):
                        ordered = group.sort_values("frame")
                        points = [
                            ((float(row.x)-x0)*scale + image_left, (float(row.y)-y0)*scale+y_offset)
                            for row in ordered.itertuples()
                        ]
                        if len(points) > 1:
                            draw.line(points, fill=colors[int(track_id) % len(colors)], width=max(1, scale))
                _draw_time_and_scale(draw, width, header + height, int(frame_index), total-1,
                                     frame_interval_s, pixel_size_um, scale, scale_bar_um, small,
                                     image_left, image_width)
                writer.send(np.asarray(canvas, dtype=np.uint8).tobytes())
    finally:
        writer.close()
    return indices.tolist()
