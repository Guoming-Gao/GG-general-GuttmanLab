from __future__ import annotations

import numpy as np
from tifffile import imwrite

from oni_spt.preprocess import artifact_paths, channel_specs
from oni_spt.segmentation import choose_segmentation_image, normalize01, region_table


def test_marker_is_preferred_and_regions_are_exportable(tmp_path):
    cfg = {
        "output_dir": str(tmp_path),
        "layout": {"mode": "dual"},
        "channels": {
            "left": {"name": "hoechst", "role": "marker"},
            "right": {"name": "spt", "role": "spt"},
        },
        "segmentation": {"source": "auto"},
        "bandpass": {"sigma_low": 1, "sigma_high": 3},
    }
    row = {"fov": "sample"}
    for spec in channel_specs(cfg):
        path = artifact_paths(tmp_path, "sample", spec, cfg)["mip"]
        path.parent.mkdir(parents=True, exist_ok=True)
        imwrite(path, np.full((5, 5), 3 if spec["role"] == "marker" else 7, np.uint16))
    image, source = choose_segmentation_image(row, cfg)
    assert source == "hoechst"
    assert np.all(image == 3)
    masks = np.zeros((5, 5), np.uint16)
    masks[1:3, 1:4] = 1
    table = region_table(masks)
    assert table.loc[0, "cell_id"] == 1
    assert table.loc[0, "area_px"] == 6


def test_dual_spt_builds_shared_percentile_normalized_mip(tmp_path):
    cfg = {
        "output_dir": str(tmp_path),
        "layout": {"mode": "dual"},
        "channels": {
            "left": {"name": "spt_left", "role": "spt"},
            "right": {"name": "spt_right", "role": "spt"},
        },
        "segmentation": {"source": "combined_spt_mip"},
        "bandpass": {"sigma_low": 1, "sigma_high": 3},
    }
    row = {"fov": "sample"}
    images = [np.arange(25, dtype=np.uint16).reshape(5, 5), np.fliplr(np.arange(25, dtype=np.uint16).reshape(5, 5)) * 10]
    for spec, image in zip(channel_specs(cfg), images):
        path = artifact_paths(tmp_path, "sample", spec, cfg)["mip"]
        path.parent.mkdir(parents=True, exist_ok=True)
        imwrite(path, image)
    combined, source = choose_segmentation_image(row, cfg)
    assert source == "combined_spt_mip"
    np.testing.assert_allclose(combined, np.maximum(normalize01(images[0]), normalize01(images[1])))
