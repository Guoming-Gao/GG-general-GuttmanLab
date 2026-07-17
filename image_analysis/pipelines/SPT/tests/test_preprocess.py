from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from skimage.filters import difference_of_gaussians
from tifffile import imread

from spt_shared import load_config
from step01_inspect_inputs import build_manifest
from step02_preprocess_videos import artifact_paths, bandpass_frame, channel_specs, preprocess_file, split_frame

from conftest import write_config, write_oni_tiff


def test_split_and_bandpass_exact(tmp_path):
    frame = np.arange(8 * 16, dtype=np.uint16).reshape(8, 16)
    halves = split_frame(frame, "dual")
    assert halves["left"].shape == (8, 8)
    np.testing.assert_array_equal(halves["right"], frame[:, 8:])
    with pytest.raises(ValueError):
        split_frame(frame[:, :-1], "dual")
    expected = difference_of_gaussians(frame.astype(np.float32), 1, 3).astype(np.float32)
    np.testing.assert_array_equal(bandpass_frame(frame, 1, 3), expected)
    assert expected.dtype == np.float32
    assert expected.min() < 0


def test_preprocess_streams_raw_filtered_and_mip_and_resume(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    rng = np.random.default_rng(4)
    data = rng.integers(0, 2000, (4, 10, 12), dtype=np.uint16)
    source = raw / "movie.tif"
    write_oni_tiff(source, data)
    cfg = load_config(write_config(tmp_path / "config.yaml", raw, tmp_path / "out", mode="single", width=12))
    row = build_manifest(cfg).iloc[0].to_dict()
    preprocess_file(row, cfg, max_frames=3)
    spec = channel_specs(cfg)[0]
    paths = artifact_paths(cfg["output_dir"], row["fov"], spec, cfg)
    raw_out = imread(paths["raw"])
    filtered = imread(paths["filtered"])
    assert raw_out.shape == (3, 10, 12)
    np.testing.assert_array_equal(raw_out, data[:3])
    np.testing.assert_array_equal(filtered[0], bandpass_frame(data[0], 1, 3))
    assert filtered.dtype == np.float32
    assert filtered.min() < 0
    assert np.any(filtered != np.rint(filtered))
    statistics = pd.read_csv(paths["stats"])
    assert {"percentile_0_1", "percentile_50", "percentile_99_9"} <= set(statistics)
    assert statistics.loc[0, "negative_fraction"] > 0
    np.testing.assert_array_equal(imread(paths["mip"]), data[:3].max(axis=0))
    mtime = paths["filtered"].stat().st_mtime_ns
    preprocess_file(row, cfg, max_frames=3, resume=True)
    assert paths["filtered"].stat().st_mtime_ns == mtime
    with pytest.raises(FileExistsError):
        preprocess_file(row, cfg, max_frames=3)
    preprocess_file(row, cfg, max_frames=3, force=True)
    np.testing.assert_array_equal(imread(paths["filtered"])[0], bandpass_frame(data[0], 1, 3))
