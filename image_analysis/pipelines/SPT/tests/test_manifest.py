from __future__ import annotations

import numpy as np

from oni_spt.config import load_config
from oni_spt.manifest import build_manifest

from conftest import write_config, write_oni_tiff


def test_manifest_accepts_profile_and_reports_wrong_width_and_truncation(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    write_oni_tiff(raw / "accepted.tif", np.zeros((3, 8, 16), np.uint16), expected_frames=5)
    write_oni_tiff(raw / "wrong.tif", np.zeros((3, 8, 8), np.uint16))
    cfg = load_config(write_config(tmp_path / "config.yaml", raw, tmp_path / "out", mode="dual", width=16))
    manifest = build_manifest(cfg).set_index("filename")
    assert manifest.loc["accepted.tif", "status"] == "accepted"
    assert "page_count_mismatch:3!=5" in manifest.loc["accepted.tif", "warnings"]
    assert manifest.loc["wrong.tif", "status"] == "skipped"
    assert manifest.loc["wrong.tif", "reason"] == "unexpected_single_channel_width"
    assert manifest.loc["accepted.tif", "frame_interval_s"] == 0.03
