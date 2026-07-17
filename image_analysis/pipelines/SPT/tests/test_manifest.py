from __future__ import annotations

import numpy as np

from spt_shared import load_config
from step01_inspect_inputs import build_manifest

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


def test_manifest_condition_override_accepts_filename_or_stem(tmp_path):
    raw = tmp_path / "raw"; raw.mkdir()
    write_oni_tiff(raw / "unstructured_name.tif", np.zeros((2, 8, 8), np.uint16))
    cfg = load_config(write_config(tmp_path / "config.yaml", raw, tmp_path / "out", mode="single", width=8))
    cfg["analysis"]["condition_overrides"] = {"unstructured_name.tif": "manual-condition"}
    manifest = build_manifest(cfg)
    assert manifest.loc[0, "condition"] == "manual-condition"
    assert "filename_has_no_FOV_token" not in manifest.loc[0, "warnings"]
