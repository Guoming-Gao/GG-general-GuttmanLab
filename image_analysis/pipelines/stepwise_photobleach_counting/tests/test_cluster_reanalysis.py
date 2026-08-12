from __future__ import annotations

import numpy as np
import pandas as pd

from cluster_counting import bootstrap_ci, high_confidence_events, monomer_candidate_mask
from cluster_reanalysis import _annulus_indices, _disk_indices, build_seed_and_group_units, circle_overlap_fraction
from cluster_reporting import count_summary, make_histogram


def sample_seeds() -> pd.DataFrame:
    return pd.DataFrame({
        "seed_id": [1, 2, 3, 4],
        "y_px": [10.0, 10.0, 30.0, 50.0],
        "x_px": [10.0, 13.0, 30.0, 50.0],
        "spotiflow_probability": [.9, .8, .7, .6],
        "legacy_growth_status": ["seed_below_local_threshold", "grown", "region_too_small", "grown_with_warning"],
    })


def test_every_seed_gets_both_coprimary_units_even_when_growth_failed():
    seeds, seed_units, groups, mapping = build_seed_and_group_units(sample_seeds(), "FOV", 4.0, 2.0)
    assert len(seeds) == len(seed_units) == len(mapping) == 4
    assert mapping.seed_level_unit_id.notna().all() and mapping.grouped_unit_id.notna().all()
    failed = seeds.legacy_growth_status.isin(["seed_below_local_threshold", "region_too_small"])
    assert seeds.loc[failed, "seed_level_unit_id"].notna().all()
    assert len(groups) == 3


def test_overlapping_seed_groups_are_deterministic_and_preserve_members():
    first = build_seed_and_group_units(sample_seeds(), "FOV", 4.0, 2.0)
    shuffled = build_seed_and_group_units(sample_seeds().sample(frac=1, random_state=3), "FOV", 4.0, 2.0)
    first_groups = first[2].sort_values("analysis_unit_id").reset_index(drop=True)
    second_groups = shuffled[2].sort_values("analysis_unit_id").reset_index(drop=True)
    pd.testing.assert_series_equal(first_groups.member_seed_ids, second_groups.member_seed_ids)
    assert first_groups.iloc[0].member_seed_ids == "1;2"


def test_signal_overlap_fraction_has_expected_limits():
    assert circle_overlap_fraction(4.0, 2.0) == 0
    assert circle_overlap_fraction(0.0, 2.0) == 1
    assert 0 < circle_overlap_fraction(2.0, 2.0) < 1


def test_high_confidence_events_find_late_single_steps():
    rng = np.random.default_rng(8); frames = 200; records = []
    for index in range(35):
        values = np.r_[np.full(130, 100.0), np.full(70, 0.0)] + rng.normal(0, 1, frames)
        records.append({"dataset": "D", "fov": "F", "analysis_unit_id": f"u{index}", "signal_pixels": 9, **{str(i): value for i, value in enumerate(values)}})
    events = high_confidence_events(pd.DataFrame.from_records(records))
    assert len(events) >= 30
    assert np.median(events.event_amplitude) > 90


def test_high_confidence_events_reject_isolated_noise_spike():
    rng = np.random.default_rng(12)
    values = rng.normal(0, 1, 240)
    values[180] = -80
    table = pd.DataFrame([{"dataset": "D", "fov": "F", "analysis_unit_id": "spike",
                           "signal_pixels": 9, **{str(i): value for i, value in enumerate(values)}}])
    assert high_confidence_events(table).empty


def test_bootstrap_step_interval_contains_median():
    values = np.arange(80, 121, dtype=float)
    low, high = bootstrap_ci(values, 400)
    assert low < np.median(values) < high


def test_background_annulus_excludes_neighbor_signal_footprint():
    shape = (25, 25)
    excluded = np.zeros(shape, dtype=bool)
    neighbor = _disk_indices(shape, 12, 16, 2)
    excluded.ravel()[neighbor] = True
    background = _annulus_indices(shape, 12, 12, 3.5, 5, excluded)
    assert not set(background).intersection(set(neighbor))


def test_monomer_candidates_require_all_validation_fields_and_cluster_stream():
    table = pd.DataFrame({"pbsa_flag": [1, 1, -2], "photobleaching_step_count": [1, 2, 1],
                          "valid_background": [True, True, True], "finite_trace": [True, True, True],
                          "unresolved_baseline_failure": [False, False, True]})
    assert monomer_candidate_mask(table, "seed_level_cluster").tolist() == [True, False, False]
    assert not monomer_candidate_mask(table, "extended_cluster_candidate").any()


def test_coprimary_streams_remain_separate_in_summary():
    base = pd.DataFrame({"analysis_set": ["primary"], "condition": ["SHA noDox"],
                         "accepted_count": [True], "monomer_candidate": [True],
                         "photobleaching_step_count": [1], "validated_range_status": ["within_validated_range"]})
    summary = count_summary({"seed_level_cluster": base.copy(), "grouped_cluster": base.copy()})
    assert set(summary.analysis_stream) == {"seed_level_cluster", "grouped_cluster"}
    assert len(summary) == 2


def test_revised_histogram_uses_bottom_four_column_legend(tmp_path, monkeypatch):
    from matplotlib.figure import Figure
    captured = {}; original = Figure.legend
    def spy(self, *args, **kwargs):
        captured.update(kwargs); return original(self, *args, **kwargs)
    monkeypatch.setattr(Figure, "legend", spy)
    folder = tmp_path / "02_figures" / "seed_level_clusters"; folder.mkdir(parents=True)
    table = pd.DataFrame({"condition": ["SHA noDox", "SHA Dox", "dSPEN FL", "dSPEN dRRM"],
                          "is_primary_comparison": [True] * 4, "accepted_count": [True] * 4,
                          "photobleaching_step_count": [1, 2, 3, 4]})
    make_histogram({"output_root": str(tmp_path)}, table, "seed_level_cluster")
    assert captured["ncol"] == 4
    assert captured["loc"] == "lower center"
