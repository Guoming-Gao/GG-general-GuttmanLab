import numpy as np
import pandas as pd
from tifffile import TiffWriter

from spt_video import render_overlay_video, selected_frames
from step08_make_crop_qc_video import select_crop


def test_qc_frame_selection_is_even_and_limited_to_100():
    frames = selected_frames(1000, 100)
    assert len(frames) == 100
    assert frames[0] == 0 and frames[-1] == 999
    assert np.all(np.diff(frames) > 0)


def test_crop_selection_by_cell_and_track_transforms_bbox():
    mask = np.zeros((20, 30), np.uint16)
    mask[5:10, 8:15] = 3
    tracks = pd.DataFrame({"trackID": [7, 7], "cell_id": [3, 3], "x": [9.0, 14.0], "y": [6.0, 9.0]})
    bbox, selected = select_crop(mask, tracks, cell_id=3, track_id=None, bbox=None, padding=2)
    assert bbox == (6, 3, 17, 12)
    assert set(selected.trackID) == {7}


def test_mp4_renderer_limits_frames_and_writes_h264(tmp_path):
    stack_path, video_path = tmp_path / "dog.tif", tmp_path / "qc.mp4"
    with TiffWriter(stack_path) as writer:
        for frame in range(12):
            image = np.zeros((16, 18), np.float32)
            image[4:7, 5:8] = frame + 1
            image -= 0.25
            writer.write(image, contiguous=True, metadata=None)
    detections = pd.DataFrame({"frame": [0, 6, 11], "x": [6., 6., 6.], "y": [5., 5., 5.]})
    indices = render_overlay_video(
        stack_path, video_path, detections, tracks=None, mask=np.zeros((16, 18), np.uint16),
        frame_interval_s=0.03, pixel_size_um=0.117, model_label="Spotiflow general",
        threshold_label="0.5", maximum_frames=5, scale=2, fps=10, crf=20,
    )
    assert indices == [0, 3, 6, 8, 11]
    assert video_path.exists() and video_path.stat().st_size > 1000
    import imageio_ffmpeg
    frames, seconds = imageio_ffmpeg.count_frames_and_secs(video_path)
    assert frames == 5
    assert seconds > 0
    reader = imageio_ffmpeg.read_frames(video_path)
    metadata = next(reader)
    reader.close()
    assert metadata["size"] == (856, 172)
