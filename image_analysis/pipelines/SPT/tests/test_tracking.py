from __future__ import annotations

import importlib.util

import numpy as np
import pandas as pd
import pytest

from oni_spt.tracking import LEGACY_COLUMNS, assign_track_labels, legacy_table, link_detections


def test_track_assignment_retains_background_and_flags_mixed_masks():
    group = pd.DataFrame({"detection_cell_id": [0, 2, 2, 2, 3]})
    result = assign_track_labels(group)
    assert result["cell_id"] == 2
    assert result["assignment_confidence"] == 0.6
    assert result["encountered_cell_ids"] == "0;2;3"
    assert result["assignment_ambiguous"]
    background = assign_track_labels(pd.DataFrame({"detection_cell_id": [0, 0]}))
    assert background["cell_id"] == 0
    assert background["assignment_confidence"] == 1.0


@pytest.mark.skipif(importlib.util.find_spec("laptrack") is None, reason="laptrack not installed")
def test_laptrack_cutoff_minimum_length_and_legacy_contract():
    rows = []
    for frame in range(6):
        rows.append(
            {
                "dataset": "d",
                "fov": "f",
                "channel": "spt",
                "spotID": frame,
                "frame": frame,
                "t": frame,
                "x": float(frame),
                "y": float(frame) / 2,
                "spotiflow_probability": 0.8,
                "R": 2.5,
                "meanIntensity": 10,
                "medianIntensity": 9,
                "minIntensity": 2,
                "maxIntensity": 15,
                "totalIntensity": 100,
                "stdIntensity": 3,
                "contrast": 0.5,
                "SNR": 3,
                "detection_cell_id": 1,
            }
        )
    cfg = {
        "dataset": "d",
        "laptrack": {
            "metric": "sqeuclidean",
            "max_link_distance_px": 5,
            "gap_closing": False,
            "min_track_length": 5,
        },
    }
    tracks, summary = link_detections(pd.DataFrame(rows), cfg)
    assert len(summary) == 1
    assert summary.loc[0, "max_step_px"] <= 5
    legacy = legacy_table(tracks)
    assert list(legacy.columns) == LEGACY_COLUMNS
    assert all(np.issubdtype(dtype, np.floating) for dtype in legacy.dtypes)
