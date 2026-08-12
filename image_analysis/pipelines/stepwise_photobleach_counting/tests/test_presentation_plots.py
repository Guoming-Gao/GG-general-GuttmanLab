from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from presentation_plots import (
    CONDITION_ORDER,
    DISPLAY_LABELS,
    FONT_SIZE,
    empirical_cdf,
    height_refit_idealization,
    histogram_percentages,
    load_presentation_counts,
    load_trace_fit,
    make_cdf_figure,
    make_histogram_figure,
    make_trace_panel,
    select_highest_step_candidates,
    select_monomer_candidates,
)


def test_presentation_filtering_and_display_labels(tmp_path: Path):
    tables = tmp_path / "03_tables"; tables.mkdir()
    rows = []
    for condition in (*CONDITION_ORDER, "dSPEN dRRM"):
        rows.append({
            "condition": condition, "analysis_set": "primary_standardized", "accepted_count": True,
            "photobleaching_step_count": 2, "analysis_unit_id": condition,
        })
    rows.extend([
        {"condition": "SHA noDox", "analysis_set": "acquisition_test", "accepted_count": True, "photobleaching_step_count": 1, "analysis_unit_id": "test"},
        {"condition": "SHA Dox", "analysis_set": "primary_standardized", "accepted_count": False, "photobleaching_step_count": 1, "analysis_unit_id": "rejected"},
    ])
    pd.DataFrame(rows).to_csv(tables / "grouped_cluster_counts.csv", index=False)
    selected = load_presentation_counts({"output_root": str(tmp_path)})
    assert tuple(selected.condition) == CONDITION_ORDER
    assert set(selected.display_condition) == set(DISPLAY_LABELS.values())
    assert "dSPEN dRRM" not in selected.condition.values
    assert selected.analysis_set.eq("primary_standardized").all()


def test_histogram_percentages_and_ecdf_are_normalized():
    table = pd.DataFrame({
        "condition": np.repeat(CONDITION_ORDER, [3, 4, 5]),
        "photobleaching_step_count": [1, 2, 3, 1, 1, 2, 4, 2, 2, 3, 4, 5],
    })
    _, percentages = histogram_percentages(table)
    for condition in CONDITION_ORDER:
        assert np.isclose(percentages[condition].sum(), 100.0)
        _, y_values = empirical_cdf(table.loc[table.condition.eq(condition), "photobleaching_step_count"])
        assert y_values[0] == 0
        assert y_values[-1] == 1

    for figure in (make_histogram_figure(table), make_cdf_figure(table)):
        visible_text = [text for text in figure.findobj(match=plt.Text) if text.get_text()]
        assert all(np.isclose(text.get_fontsize(), FONT_SIZE) for text in visible_text)
        plt.close(figure)


def test_height_refit_matches_each_contiguous_plateau_median():
    corrected = np.array([11, 13, 12, 7, 9, 8, 2, 4], dtype=float)
    states = np.array([2, 2, 2, 1, 1, 1, 0, 0], dtype=float)
    idealized = height_refit_idealization(corrected, states)
    np.testing.assert_allclose(idealized, [12, 12, 12, 8, 8, 8, 3, 3])


def test_trace_loader_uses_the_counted_representation(tmp_path: Path):
    dataset = "dataset"; fov = "fov"; stream = "grouped_cluster"
    trace_dir = tmp_path / "04_experiments" / dataset / "puncta_traces" / stream / fov
    native_dir = tmp_path / "99_technical" / "quickpbsa_native" / dataset / stream / fov
    trace_dir.mkdir(parents=True); native_dir.mkdir(parents=True)
    representation = "baseline_centered_integrated_difference"
    pd.DataFrame([{
        "analysis_unit_id": "unit-1", "0": 102.0, "1": 98.0, "2": 11.0, "3": 9.0,
    }]).to_csv(trace_dir / f"{fov}__{representation}.csv", index=False)
    native_table = pd.DataFrame([{
        "analysis_unit_id": "unit-1", "type": "fluors_full", "0": 1.0, "1": 1.0, "2": 0.0, "3": 0.0,
    }])
    native_path = native_dir / f"{fov}__{stream}__quickpbsa_input_result.csv"
    native_path.write_text("metadata\n" + native_table.to_csv(index=False))
    row = pd.Series({
        "dataset": dataset, "fov": fov, "analysis_unit_id": "unit-1", "trace_representation": representation,
    })
    fit = load_trace_fit({"output_root": str(tmp_path)}, row, stream)
    np.testing.assert_allclose(fit.corrected, [102, 98, 11, 9])
    np.testing.assert_allclose(fit.idealized, [100, 100, 10, 10])


def _selection_table() -> pd.DataFrame:
    rows = []
    for condition_index, condition in enumerate(CONDITION_ORDER):
        for index in range(6):
            rows.append({
                "condition": condition, "dataset": "dataset", "fov": f"{condition_index}-fov-{index}",
                "analysis_unit_id": f"{condition_index}-{index}", "fit_residual": 0.6 - index * 0.05,
                "single_step_amplitude_relative_error": index / 10,
                "photobleaching_step_count": 10 + index,
            })
    return pd.DataFrame(rows)


def test_trace_selection_is_deterministic_and_balanced():
    table = _selection_table()
    first = select_monomer_candidates(table)
    second = select_monomer_candidates(table.sample(frac=1, random_state=4))
    assert first.analysis_unit_id.tolist() == second.analysis_unit_id.tolist()
    assert first.groupby("condition").size().eq(4).all()

    highest = select_highest_step_candidates(table)
    assert highest.groupby("condition").size().eq(4).all()
    assert highest.groupby("condition").fov.nunique().eq(4).all()
    assert highest.groupby("condition").photobleaching_step_count.max().eq(15).all()


def test_trace_panel_is_four_by_three_with_11_point_text():
    records = []
    for condition in CONDITION_ORDER:
        for rank in range(1, 5):
            records.append({
                "condition": condition, "selection_rank": rank, "photobleaching_step_count": rank,
                "corrected": np.array([10, 9, 3, 2], dtype=float),
                "idealized": np.array([9.5, 9.5, 2.5, 2.5], dtype=float),
            })
    figure = make_trace_panel(records, "Test panel")
    assert len(figure.axes) == 12
    visible_text = [text for text in figure.findobj(match=plt.Text) if text.get_text()]
    assert visible_text
    assert all(np.isclose(text.get_fontsize(), FONT_SIZE) for text in visible_text)
    plt.close(figure)
