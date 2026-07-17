from __future__ import annotations

import numpy as np
import pandas as pd

from oni_spt.analysis import calculate_aio_table, concat_aio_files, reformat_aio_for_saspt


def test_aio_accepts_legacy_and_propagates_new_cell_metadata():
    frame = np.arange(6)
    tracks = pd.DataFrame(
        {
            "trackID": np.zeros(6),
            "x": [0, 1, 2, 3, 4, 5],
            "y": [0, 1, 0, 2, 1, 3],
            "t": frame,
            "meanIntensity": [10, 11, 12, 13, 14, 15],
            "dataset": "test",
            "fov": "fov1",
            "channel": "spt",
            "cell_id": 7,
            "assignment_confidence": 1.0,
            "assignment_ambiguous": False,
            "trajectory_uid": "test/fov1/spt/0",
        }
    )
    aio = calculate_aio_table(tracks, frame_interval_s=0.03)
    assert len(aio) == 1
    assert aio.loc[0, "cell_id"] == 7
    assert aio.loc[0, "N_steps"] == 6
    saspt = reformat_aio_for_saspt(aio)
    assert list(saspt.columns) == ["x", "y", "trajectory", "frame"]
    assert len(saspt) == 6


def test_old_trackmate_table_and_aio_concat_remain_supported(tmp_path):
    legacy = pd.DataFrame(
        {
            "spotID": np.arange(6, dtype=float),
            "trackID": np.zeros(6, dtype=float),
            "x": [0, 1, 2, 3, 4, 5],
            "y": [0, 1, 0, 2, 1, 3],
            "t": np.arange(6, dtype=float),
            "R": 2.5,
            "meanIntensity": np.arange(10, 16, dtype=float),
        }
    )
    aio = calculate_aio_table(legacy, frame_interval_s=0.03)
    assert len(aio) == 1
    assert "cell_id" not in aio.columns
    first = tmp_path / "first.csv"
    second = tmp_path / "second.csv"
    aio.to_csv(first, index=False)
    aio.to_csv(second, index=False)
    output = concat_aio_files([first, second], tmp_path / "condition.csv")
    combined = pd.read_csv(output)
    assert len(combined) == 2
    assert set(combined["filename"]) == {"first.csv", "second.csv"}
