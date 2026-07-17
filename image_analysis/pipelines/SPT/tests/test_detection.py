from types import SimpleNamespace

import numpy as np
from tifffile import TiffFile, TiffWriter

from step04_detect_spots import _detail_value, detect_stack


def test_spotiflow_detail_value_accepts_scalar_and_single_channel_arrays():
    details = SimpleNamespace(
        prob=np.array([0.7, 0.8]),
        intens=np.array([[10.0], [12.0]]),
    )
    assert _detail_value(details, "prob", 1) == 0.8
    assert _detail_value(details, "intens", 1) == 12.0
    assert np.isnan(_detail_value(details, "missing", 0))


def test_spotiflow_receives_exact_saved_float32_frames(tmp_path):
    raw_path, dog_path = tmp_path / "raw.tif", tmp_path / "dog.tif"
    raw = np.arange(2 * 8 * 9, dtype=np.uint16).reshape(2, 8, 9)
    dog = np.linspace(-7.25, 12.75, raw.size, dtype=np.float32).reshape(raw.shape)
    for path, stack in ((raw_path, raw), (dog_path, dog)):
        with TiffWriter(path) as writer:
            for frame in stack:
                writer.write(frame, contiguous=True, metadata=None)

    class RecordingModel:
        def __init__(self):
            self.frames = []

        def predict(self, frame, **kwargs):
            self.frames.append(frame.copy())
            return np.empty((0, 2)), SimpleNamespace(prob=[], intens=[])

    model = RecordingModel()
    cfg = {
        "legacy": {"spot_radius_px": 2.5},
        "spotiflow": {"probability_threshold": None, "min_distance": 1,
                       "exclude_border": 1, "normalizer": "auto"},
    }
    detect_stack(raw_path, dog_path, np.zeros(raw.shape[1:], np.uint16), cfg,
                 SimpleNamespace(model=model, device="cpu"))
    assert len(model.frames) == 2
    with TiffFile(dog_path) as tif:
        for recorded, page in zip(model.frames, tif.pages):
            np.testing.assert_array_equal(recorded, page.asarray())
            assert recorded.dtype == np.float32
    assert np.any(np.concatenate([frame.ravel() for frame in model.frames]) < 0)
