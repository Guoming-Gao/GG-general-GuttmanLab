from __future__ import annotations

import numpy as np
import pandas as pd

from pbsa_shared import condition_metadata
from step05_quickpbsa_threshold_pilot import select_half_step_threshold
from step07_generate_reports import CONDITION_ORDER, make_condition_histogram


def test_condition_normalization_and_primary_assignment():
    primary = condition_metadata("standarized-dSPEN_dRRM-FOV-1.tif")
    acquisition_test = condition_metadata("SHA_noDox-FOV-2.tif")
    assert primary == {
        "condition_key": "dSPEN_dRRM", "condition": "dSPEN dRRM",
        "analysis_set": "primary_standardized", "is_primary_comparison": True,
    }
    assert acquisition_test["condition"] == "SHA noDox"
    assert acquisition_test["analysis_set"] == "acquisition_test"
    assert acquisition_test["is_primary_comparison"] is False


def test_automatic_threshold_is_half_the_estimated_step():
    assert select_half_step_threshold(124.0) == 62.0
    for invalid in [0.0, -1.0, np.nan]:
        try:
            select_half_step_threshold(invalid)
        except RuntimeError:
            pass
        else:
            raise AssertionError("Invalid step scales must block automatic counting")


def test_condition_histogram_has_bottom_multicolumn_legend(tmp_path):
    records = []
    for condition_index, condition in enumerate(CONDITION_ORDER):
        for count in range(1, 5):
            records.append({
                "condition": condition, "roi_class": "point_punctum",
                "photobleaching_step_count": count + condition_index,
                "accepted_count": True, "is_primary_comparison": True,
            })
    counts = pd.DataFrame.from_records(records); counts.attrs["output_dir"] = str(tmp_path)
    fig, ax = make_condition_histogram(counts, "point_punctum")
    assert len(fig.legends) == 1
    assert fig.legends[0]._ncols == 4
    assert fig.subplotpars.bottom >= 0.2
    assert len(ax.patches) > 0
