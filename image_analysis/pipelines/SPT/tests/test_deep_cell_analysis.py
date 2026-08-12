from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from tifffile import imwrite

from step09_deep_cell_analysis import (
    DEFAULT_CONFIG,
    _save,
    _expanded_bbox,
    _phase,
    cell_metrics,
    cliffs_delta,
    comparison_dir,
    extract_adjacent_steps,
    holm_adjust,
    load_deep_config,
    match_timing_cells,
    phase_sample_label,
    phase_sample_sizes,
    remove_legacy_timing_outputs,
    selected_cells_dir,
    selected_overlay_name,
    step_distribution_curves,
    timing_comparison_aio,
    timing_comparison_dir,
    trajectory_metric_cdf,
    trajectory_segment_style,
)


def test_phase_mapping_checks_recovery_before_fvp():
    phases = {"Before": "SHA-before", "FVP 2 h": "SHA-FVP2h", "Recovery 2 h": "SHA-FVP2h_recover2h"}
    assert _phase("SHA-FVP2h_recover2h-SPEN_JFX650-30ms", phases) == "Recovery 2 h"
    assert _phase("SHA-FVP2h-SPEN_JFX650-30ms", phases) == "FVP 2 h"


def test_holm_and_cliffs_delta():
    adjusted = holm_adjust([0.01, 0.04, 0.03])
    assert np.allclose(adjusted, [0.03, 0.06, 0.06])
    assert cliffs_delta(np.array([3, 4]), np.array([1, 2])) == 1.0
    assert cliffs_delta(np.array([1, 2]), np.array([3, 4])) == -1.0


def test_bbox_filename_and_count_prefix():
    row = SimpleNamespace(
        bbox_min_x=1, bbox_min_y=2, bbox_max_x=8, bbox_max_y=9,
        retained_trajectories=57, timing_ms=30, phase="FVP 2 h", fov="sample-FOV-1", cell_id=4,
    )
    assert _expanded_bbox(row, (10, 10), 5) == (0, 0, 10, 10)
    assert selected_overlay_name(row).startswith("00057_trajectories__30ms__FVP_2_h__")


def test_cell_metrics_uses_existing_quality_filters():
    table = pd.DataFrame({
        "timing_ms": [30] * 4, "phase": ["Before"] * 4, "fov": ["fov"] * 4,
        "cell_id": [1] * 4, "retained_trajectories_in_cell": [60] * 4,
        "N_steps": [5, 6, 7, 8], "displacement_nm": [10, 20, 30, 40],
        "mean_stepsize_nm": [10, 40, 50, 60], "max_d_anytwo_nm": [20, 30, 40, 50],
        "linear_fit_sigma": [1, 2, 3, 4], "linear_fit_log10D": [-2, -1.5, -1, -.5],
        "linear_fit_R2": [1, 1, .1, 1], "loglog_fit_R2": [1, 1, 1, 1],
        "alpha": [.4, .6, .8, 1.0],
    })
    cfg = {"analysis": {"immobile_stepsize_nm": 30, "alpha_threshold": .7, "fit_r2_threshold": .7}}
    result = cell_metrics(table, cfg).iloc[0]
    assert result.retained_trajectories == 60
    assert result.immobile_fraction == .25
    assert result.constrained_fraction == .25
    assert result.normal_fraction == .5
    # D keeps mobile, R2-valid, alpha > 0.5 rows: indices 1 and 3.
    assert result.linear_fit_log10D == -1.0


def test_mask_matching_accepts_only_selected_high_iou_pairs(tmp_path: Path):
    segmentation = tmp_path / "02_segmentation"
    segmentation.mkdir()
    fov30 = "SHA-FVP2h-SPEN_JFX650-30ms-FOV-1"
    fov100 = "SHA-FVP2h-SPEN_JFX650-100ms-FOV-1"
    left = np.zeros((12, 12), np.uint16); left[2:8, 2:8] = 1; left[9:11, 9:11] = 2
    right = np.zeros((12, 12), np.uint16); right[2:8, 2:8] = 5; right[9:11, 9:11] = 6
    imwrite(segmentation / f"{fov30}__cellposeSAM_masks.tif", left)
    imwrite(segmentation / f"{fov100}__cellposeSAM_masks.tif", right)
    inventory = pd.DataFrame([
        {"phase": "FVP 2 h", "timing_ms": 30, "fov": fov30, "cell_id": 1, "selected": True, "retained_trajectories": 60},
        {"phase": "FVP 2 h", "timing_ms": 30, "fov": fov30, "cell_id": 2, "selected": False, "retained_trajectories": 10},
        {"phase": "FVP 2 h", "timing_ms": 100, "fov": fov100, "cell_id": 5, "selected": True, "retained_trajectories": 70},
        {"phase": "FVP 2 h", "timing_ms": 100, "fov": fov100, "cell_id": 6, "selected": True, "retained_trajectories": 80},
    ])
    result = match_timing_cells({"result_root": str(tmp_path), "pairing_minimum_iou": .5}, inventory)
    accepted = result[result.accepted_pair]
    assert len(accepted) == 1
    assert accepted.iloc[0].cell_id_30ms == 1
    assert accepted.iloc[0].cell_id_100ms == 5
    assert accepted.iloc[0].mask_iou == 1.0


def test_output_folder_routing_is_review_gated(tmp_path: Path):
    cfg = {**DEFAULT_CONFIG, "result_root": str(tmp_path)}
    assert comparison_dir(cfg, 30) == tmp_path / "30ms_comparison"
    assert comparison_dir(cfg, 100) == tmp_path / "100ms_comparison"
    assert timing_comparison_dir(cfg) == tmp_path / "30ms_vs_100ms_comparison"
    assert selected_cells_dir(cfg, 30) == tmp_path / "selected_cells" / "30ms"


def test_schema_two_local_config_migrates_without_losing_overrides(tmp_path: Path):
    for folder in ("00_run_metadata", "02_segmentation", "04_trajectories", "05_diffusion_analysis"):
        (tmp_path / folder).mkdir()
    config = tmp_path / "config.yaml"
    config.write_text(f"schema_version: 2\nresult_root: {tmp_path}\ntrajectory_cutoff: 77\n")
    loaded = load_deep_config(config)
    assert loaded["schema_version"] == 3
    assert loaded["trajectory_cutoff"] == 77
    assert loaded["step_size_distribution"]["link_policy"] == "adjacent_only"


def test_alpha_style_clips_uses_plasma_and_ignores_gap_type():
    from matplotlib import colormaps, colors

    cfg = {**DEFAULT_CONFIG, "trajectory_linewidth_pt": 0.1}
    low_adjacent = trajectory_segment_style(0.1, 1, cfg)
    low_gap = trajectory_segment_style(0.1, 2, cfg)
    high = trajectory_segment_style(4.0, 1, cfg)
    missing = trajectory_segment_style(np.nan, 2, cfg)
    assert colors.to_rgba(low_adjacent["color"]) == colors.to_rgba(colormaps["plasma"](0.0))
    assert colors.to_rgba(high["color"]) == colors.to_rgba(colormaps["plasma"](1.0))
    assert low_adjacent == low_gap
    assert low_gap["linestyle"] == "-"
    assert colors.to_rgba(missing["color"]) == colors.to_rgba("#9e9e9e")


def test_sample_labels_include_multiline_cell_and_trajectory_counts():
    metrics = pd.DataFrame([
        {"phase": "Before", "fov": "a", "cell_id": 1},
        {"phase": "Before", "fov": "a", "cell_id": 2},
        {"phase": "FVP 2 h", "fov": "b", "cell_id": 1},
        {"phase": "Recovery 2 h", "fov": "c", "cell_id": 1},
    ])
    aio = pd.DataFrame([
        {"phase": "Before", "fov": "a", "trackID": 1},
        {"phase": "Before", "fov": "a", "trackID": 2},
        {"phase": "FVP 2 h", "fov": "b", "trackID": 1},
        {"phase": "Recovery 2 h", "fov": "c", "trackID": 1},
    ])
    label = phase_sample_label("Before", phase_sample_sizes(metrics, aio))
    assert label == "Before\nN=2 cells\nn=2 trajectories"
    assert "AIO" not in label


def test_adjacent_steps_exclude_gaps_and_convert_pixels_to_nm():
    aio = pd.DataFrame([{
        "timing_ms": 30, "phase": "Before", "fov": "fov", "cell_id": 1, "trackID": 7,
        "list_of_t": "[0, 1, 3, 4]", "list_of_x": "[0, 1, 4, 4]",
        "list_of_y": "[0, 0, 0, 2]", "pixel_size_um": 0.1,
    }])
    steps = extract_adjacent_steps(aio)
    assert steps.frame_difference.tolist() == [1, 1]
    assert steps.source_frame.tolist() == [0, 3]
    assert np.allclose(steps.step_size_nm, [100, 200])


def test_step_curves_pool_every_adjacent_step_and_make_valid_cdf():
    aio = pd.DataFrame([
        {"timing_ms": 30, "phase": "Before", "fov": "a", "cell_id": 1, "trackID": 1,
         "list_of_t": "[0, 1, 2]", "list_of_x": "[0, 1, 2]", "list_of_y": "[0, 0, 0]",
         "pixel_size_um": .1},
        {"timing_ms": 30, "phase": "Before", "fov": "a", "cell_id": 2, "trackID": 2,
         "list_of_t": "[0, 1]", "list_of_x": "[0, 3]", "list_of_y": "[0, 0]",
         "pixel_size_um": .1},
    ])
    curves = step_distribution_curves(aio, DEFAULT_CONFIG, cohort="test")
    histogram = curves[curves.curve == "histogram"]
    cdf = curves[curves.curve == "cdf"].sort_values("x_nm")
    assert histogram.steps.iloc[0] == 3
    assert histogram.trajectories.iloc[0] == 2
    assert histogram.cells.iloc[0] == 2
    assert np.isclose(histogram.y.sum(), 1)
    assert np.all(np.diff(cdf.y) >= 0)
    assert cdf.y.iloc[-1] == 1


def test_trajectory_mean_step_cdf_gives_every_trajectory_equal_weight():
    aio = pd.DataFrame([
        {"phase": "Before", "fov": "a", "cell_id": 1, "trackID": 1, "mean_stepsize_nm": 100},
        {"phase": "Before", "fov": "a", "cell_id": 1, "trackID": 2, "mean_stepsize_nm": 200},
        {"phase": "Before", "fov": "a", "cell_id": 2, "trackID": 3, "mean_stepsize_nm": 300},
    ])
    curves = trajectory_metric_cdf(
        aio, "mean_stepsize_nm", DEFAULT_CONFIG, (0, 600), grid_points=7,
    )
    before = curves[curves.phase == "Before"].set_index("x_nm")
    assert np.isclose(before.loc[100, "cumulative_probability"], 1 / 3)
    assert np.isclose(before.loc[200, "cumulative_probability"], 2 / 3)
    assert np.isclose(before.loc[300, "cumulative_probability"], 1.00)
    assert before.trajectories.iloc[0] == 3
    assert before.cells.iloc[0] == 2
    assert (before.weighting == "pooled_trajectories").all()
    assert np.all(np.diff(before.cumulative_probability) >= 0)


def test_timing_step_cohort_uses_matched_treated_cells_and_independent_before():
    columns = ["timing_ms", "phase", "fov", "cell_id", "trackID"]
    aio30 = pd.DataFrame([
        [30, "Before", "before30", 1, 1],
        [30, "FVP 2 h", "fov30", 1, 2],
        [30, "FVP 2 h", "fov30", 2, 3],
    ], columns=columns)
    aio100 = pd.DataFrame([
        [100, "Before", "before100", 1, 1],
        [100, "FVP 2 h", "fov100", 5, 2],
        [100, "FVP 2 h", "fov100", 6, 3],
    ], columns=columns)
    matches = pd.DataFrame([
        {"phase": "FVP 2 h", "fov_30ms": "fov30", "cell_id_30ms": 1,
         "fov_100ms": "fov100", "cell_id_100ms": 5, "accepted_pair": True},
        {"phase": "FVP 2 h", "fov_30ms": "fov30", "cell_id_30ms": 2,
         "fov_100ms": "fov100", "cell_id_100ms": 6, "accepted_pair": False},
    ])
    result = timing_comparison_aio(aio30, aio100, matches)
    assert set(result[result.phase == "Before"].fov) == {"before30", "before100"}
    assert set(result[result.phase == "FVP 2 h"].cell_id) == {1, 5}


def test_save_writes_valid_svg_and_no_individual_pdf(tmp_path: Path):
    import xml.etree.ElementTree as ET
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    paths = _save(fig, tmp_path, "plot")
    plt.close(fig)
    assert {Path(path).suffix for path in paths} == {".png", ".svg"}
    assert not (tmp_path / "plot.pdf").exists()
    assert ET.parse(tmp_path / "plot.svg").getroot().tag.rsplit("}", 1)[-1] == "svg"


def test_safe_legacy_migration_only_removes_flat_timing_artifacts(tmp_path: Path):
    cfg = {**DEFAULT_CONFIG, "result_root": str(tmp_path)}
    legacy_plot = tmp_path / "30ms_old_plot.pdf"
    legacy_plot.write_text("old")
    flat = tmp_path / "selected_cells"
    flat.mkdir()
    legacy_overlay = flat / "00050_trajectories__30ms__Before__fov__cell-1.png"
    legacy_overlay.write_text("old")
    future_overlay = flat / "00050_trajectories__100ms__Before__fov__cell-1.png"
    future_overlay.write_text("future")
    replacement = flat / "30ms"
    replacement.mkdir()
    replacement_file = replacement / legacy_overlay.name
    replacement_file.write_text("new")
    removed = remove_legacy_timing_outputs(cfg, 30)
    assert set(removed) == {str(legacy_plot), str(legacy_overlay)}
    assert replacement_file.exists()
    assert future_overlay.exists()
