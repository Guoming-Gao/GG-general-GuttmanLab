from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import tifffile

from sacdpy.live_nuclear_puncta import (
    CONDITION_COUNTS,
    CONDITIONS,
    LiveNuclearPunctaConfig,
    PRESENTATION_BOXPLOT_METRICS,
    PRESENTATION_CONDITION_NAMES,
    PRESENTATION_FIGSIZE_INCHES,
    PRESENTATION_FONT_SIZE,
    PRESENTATION_Y_LABELS,
    _fdr_bh,
    _holm,
    _retained_and_excluded,
    _write_imagej,
    aggregate_fovs,
    bbox_with_padding,
    calculate_cell_statistics,
    calculate_statistics,
    config_from_json,
    config_to_json,
    default_puncta_params,
    discover_mips,
    fixed_threshold_params,
    histogram_edges,
    make_individual_cell_boxplots,
    parse_condition,
    pooled_area_cutoff,
    pooled_nucleus_background,
    puncta_pilot_safety_metrics,
    quantify_puncta,
    run_cellpose_frame,
    select_blinded_pilot,
    touches_fov_boundary,
)


class FakeCellpose:
    def __init__(self) -> None:
        self.image = None
        self.kwargs = None

    def eval(self, image, **kwargs):
        self.image = np.asarray(image)
        self.kwargs = kwargs
        labels = np.zeros(image.shape[1:], dtype=np.uint16)
        labels[2:7, 3:9] = 1
        return labels, None, None


class LiveNuclearPunctaTests(unittest.TestCase):
    def test_condition_parsing_order(self) -> None:
        self.assertEqual(parse_condition("SHA-before-SPEN-FOV"), "before")
        self.assertEqual(parse_condition("SHA-FVP2h-SPEN-FOV"), "FVP_2h")
        self.assertEqual(parse_condition("SHA-FVP2h_recover2h-SPEN-FOV"), "recover_2h")
        with self.assertRaises(ValueError):
            parse_condition("unknown")

    def test_discovery_requires_exact_condition_inventory_and_tyx(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for condition, count in CONDITION_COUNTS.items():
                token = {
                    "before": "SHA-before",
                    "FVP_2h": "SHA-FVP2h",
                    "recover_2h": "SHA-FVP2h_recover2h",
                }[condition]
                for index in range(count):
                    name = (
                        f"ONI-{token}-SPEN-FOV-{index}-"
                        "SACDpy-647-posXY0-MIP-TYX.tif"
                    )
                    tifffile.imwrite(
                        root / name,
                        np.zeros((7, 8, 9), dtype=np.float32),
                        imagej=True,
                        metadata={"axes": "TYX"},
                    )
            plans = discover_mips(root)
            self.assertEqual(len(plans), 21)
            self.assertEqual(
                {condition: sum(p.condition == condition for p in plans) for condition in CONDITION_COUNTS},
                CONDITION_COUNTS,
            )

    def test_bbox_padding_is_clipped(self) -> None:
        mask = np.zeros((10, 12), dtype=bool)
        mask[1:4, 2:6] = True
        self.assertEqual(bbox_with_padding(mask, 5), (0, 9, 0, 11))

    def test_border_and_pooled_area_filter_are_independent(self) -> None:
        mask = np.zeros((5, 5), dtype=bool)
        mask[0, 2] = True
        self.assertTrue(touches_fov_boundary(mask))
        cutoff = pooled_area_cutoff([10, 20, 30, 40, 50], 10)
        self.assertAlmostEqual(cutoff, 14.0)
        rows = [
            {"nucleus_key": "a", "area_px": 10, "touches_fov_boundary": False},
            {"nucleus_key": "b", "area_px": 20, "touches_fov_boundary": True},
            {"nucleus_key": "c", "area_px": 30, "touches_fov_boundary": False},
        ]
        retained, excluded = _retained_and_excluded(rows, 15.0)
        self.assertEqual([row["nucleus_key"] for row in retained], ["c"])
        self.assertEqual(
            {row["nucleus_key"]: row["exclusion_reasons"] for row in excluded},
            {"a": ["bottom_10_percent_area"], "b": ["touches_fov_boundary"]},
        )

    def test_cellpose_uses_frame_image_and_locked_segmentation_normalization(self) -> None:
        image = np.arange(120, dtype=np.float32).reshape(10, 12)
        model = FakeCellpose()
        config = LiveNuclearPunctaConfig()
        labels = run_cellpose_frame(image, model, config)
        self.assertEqual(labels.shape, image.shape)
        self.assertEqual(model.image.shape, (2, 10, 12))
        self.assertEqual(model.kwargs["normalize"], {"tile_norm_blocksize": 0})
        self.assertEqual(model.kwargs["diameter"], 140.0)

    def test_default_and_fixed_threshold_modes_do_not_mix(self) -> None:
        default = default_puncta_params()
        fixed = fixed_threshold_params(1.5)
        self.assertEqual(default.threshold_detection, 0.001)
        self.assertEqual(default.threshold_background, 0.0)
        self.assertEqual(default.threshold_segmentation, 0.001)
        self.assertEqual(default.segmentation_mode, 0)
        self.assertEqual(fixed.threshold_detection, 0.0)
        self.assertEqual(fixed.threshold_background, 1.5)
        self.assertEqual(fixed.threshold_segmentation, 1.5)
        self.assertEqual(fixed.segmentation_mode, 2)

    def test_fixed_threshold_pilot_safety_gates(self) -> None:
        passing = puncta_pilot_safety_metrics(
            [10, 20, 30],
            [0.10, 0.20, 0.70],
        )
        self.assertTrue(passing["passes_safety_gates"])
        failing = puncta_pilot_safety_metrics(
            [10, 20, 30],
            [0.10, 0.36, 0.81],
        )
        self.assertFalse(failing["passes_safety_gates"])

    def test_pooled_background_uses_full_nucleus_medians(self) -> None:
        rows = [
            {"nucleus_median_sacd": 10.0},
            {"nucleus_median_sacd": 30.0},
            {"nucleus_median_sacd": 20.0},
        ]
        self.assertEqual(pooled_nucleus_background(rows), 20.0)

    def test_uint32_instance_labels_roundtrip_without_losing_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "labels.tif"
            labels = np.zeros((8, 9), dtype=np.uint32)
            labels[2:6, 3:7] = 70000
            _write_imagej(path, labels, axes="YX", pixel_size_um=0.117)
            observed = tifffile.imread(path)
            np.testing.assert_array_equal(observed, labels)
            self.assertEqual(observed.dtype, np.uint32)

    def test_pilot_is_balanced_across_condition_and_stratum(self) -> None:
        rows = []
        for condition in CONDITION_COUNTS:
            for index in range(12):
                rows.append(
                    {
                        "condition": condition,
                        "nucleus_key": f"{condition}-{index}",
                        "nucleus_mean_sacd": index + 1,
                    }
                )
        pilot = select_blinded_pilot(rows, seed=4, per_condition_per_stratum=2)
        self.assertEqual(len(pilot), 18)
        for condition in CONDITION_COUNTS:
            selected = [row for row in pilot if row["condition"] == condition]
            self.assertEqual(len(selected), 6)
            self.assertEqual(
                {stratum: sum(row["pilot_stratum"] == stratum for row in selected)
                 for stratum in ("low", "medium", "high")},
                {"low": 2, "medium": 2, "high": 2},
            )

    def test_quantification_keeps_zero_and_corrected_metrics(self) -> None:
        image = np.full((9, 9), 10.0)
        nucleus_mask = np.zeros_like(image, dtype=bool)
        nucleus_mask[1:8, 1:8] = True
        labels = np.zeros_like(image, dtype=np.uint32)
        labels[3:5, 3:5] = 1
        image[labels == 1] = [20, 30, 40, 50]
        summary, rows = quantify_puncta(
            image,
            nucleus_mask,
            labels,
            pixel_size_um=0.1,
        )
        self.assertEqual(summary["puncta_count"], 1)
        self.assertAlmostEqual(summary["mean_corrected_puncta_brightness"], 25.0)
        self.assertAlmostEqual(summary["median_punctum_area_um2"], 0.04)
        self.assertAlmostEqual(rows[0]["corrected_integrated_sacd"], 100.0)
        self.assertFalse(rows[0]["touches_nucleus_boundary"])
        self.assertIn("nonpuncta_nucleus_mean_sacd", summary)
        self.assertNotIn("nonpuncta_core_mean_sacd", summary)
        empty, empty_rows = quantify_puncta(
            image,
            nucleus_mask,
            np.zeros_like(labels),
            pixel_size_um=0.1,
        )
        self.assertEqual(empty_rows, [])
        self.assertEqual(empty["puncta_count"], 0)
        self.assertTrue(np.isnan(empty["mean_corrected_puncta_brightness"]))
        self.assertTrue(np.isnan(empty["median_punctum_area_um2"]))

    def test_quantification_rejects_puncta_outside_full_nucleus(self) -> None:
        image = np.ones((7, 7), dtype=np.float32)
        nucleus_mask = np.zeros_like(image, dtype=bool)
        nucleus_mask[1:6, 1:6] = True
        labels = np.zeros_like(image, dtype=np.uint32)
        labels[0, 0] = 1
        with self.assertRaisesRegex(ValueError, "outside the Cellpose nucleus"):
            quantify_puncta(image, nucleus_mask, labels, pixel_size_um=0.1)

    def test_full_nucleus_boundary_contact_is_recorded(self) -> None:
        image = np.full((7, 7), 10.0, dtype=np.float32)
        nucleus_mask = np.zeros_like(image, dtype=bool)
        nucleus_mask[1:6, 1:6] = True
        labels = np.zeros_like(image, dtype=np.uint32)
        labels[1:3, 2:4] = 1
        image[labels == 1] = 20.0
        _, rows = quantify_puncta(image, nucleus_mask, labels, pixel_size_um=0.1)
        self.assertTrue(rows[0]["touches_nucleus_boundary"])

    def test_fov_aggregation_prevents_nucleus_pseudoreplication(self) -> None:
        rows = []
        for condition in CONDITION_COUNTS:
            for fov_index in range(2):
                for nucleus_index in range(3):
                    rows.append(
                        {
                            "condition": condition,
                            "fov": f"{condition}-{fov_index}",
                            "puncta_count": nucleus_index,
                            "mean_corrected_puncta_brightness": 10 + nucleus_index,
                            "median_punctum_area_um2": 0.2 + nucleus_index / 10,
                        }
                    )
        aggregated = aggregate_fovs(rows)
        self.assertEqual(len(aggregated), 6)
        self.assertTrue(all(row["fov_mean_puncta_count"] == 1 for row in aggregated))

    def test_multiple_testing_corrections_and_statistics_schema(self) -> None:
        self.assertEqual(_fdr_bh([0.01, 0.02, 0.5]), [0.03, 0.03, 0.5])
        self.assertEqual(_holm([0.01, 0.02, 0.5]), [0.03, 0.04, 0.5])
        rows = []
        for condition_index, condition in enumerate(CONDITION_COUNTS):
            for index in range(6):
                value = float(condition_index * 10 + index)
                rows.append(
                    {
                        "condition": condition,
                        "fov": f"{condition}-{index}",
                        "fov_mean_puncta_count": value,
                        "fov_median_corrected_brightness": value + 1,
                        "fov_median_punctum_area_um2": value + 2,
                    }
                )
        output = calculate_statistics(rows, seed=1, bootstrap_iterations=100)
        self.assertEqual(len(output), 12)
        self.assertTrue(all("adjusted_p_value" in row for row in output))

    def test_cell_level_anova_gate_and_mann_whitney_schema(self) -> None:
        rows = []
        for condition_index, condition in enumerate(CONDITION_COUNTS):
            for index in range(12):
                value = float(condition_index * 100 + index)
                rows.append(
                    {
                        "condition": condition,
                        "puncta_count": value,
                        "mean_corrected_puncta_brightness": value + 10,
                        "median_punctum_area_um2": value / 100 + 0.1,
                    }
                )
        output = calculate_cell_statistics(rows, seed=2, bootstrap_iterations=50)
        omnibus = [row for row in output if row["test_family"] == "omnibus"]
        pairwise = [row for row in output if row["test_family"] == "pairwise"]
        self.assertEqual(len(omnibus), 3)
        self.assertEqual(len(pairwise), 9)
        self.assertTrue(all(row["independent_unit"] == "cell" for row in output))
        self.assertTrue(
            all(row["test_name"] == "Mann-Whitney_two-sided" for row in pairwise)
        )
        self.assertTrue(all("adjusted_p_value_Holm" in row for row in pairwise))

    def test_histogram_edges_use_exact_range_and_eleven_bins(self) -> None:
        edges = histogram_edges([2.0, 3.0, 8.0])
        self.assertEqual(len(edges), 12)
        self.assertEqual(edges[0], 2.0)
        self.assertEqual(edges[-1], 8.0)
        logarithmic = histogram_edges([10.0, 100.0, 1000.0], logarithmic=True)
        self.assertEqual(len(logarithmic), 12)
        self.assertEqual(logarithmic[0], 10.0)
        self.assertEqual(logarithmic[-1], 1000.0)
        np.testing.assert_allclose(
            logarithmic[1:] / logarithmic[:-1],
            np.repeat(logarithmic[1] / logarithmic[0], 11),
        )

    def test_presentation_boxplots_use_requested_style_and_metric_scope(self) -> None:
        self.assertEqual(
            PRESENTATION_Y_LABELS["mean_corrected_puncta_brightness"],
            "Puncta intensity, A.U.",
        )
        rows = []
        for condition_index, condition in enumerate(CONDITIONS):
            for index in range(6):
                rows.append(
                    {
                        "condition": condition,
                        "puncta_count": float(condition_index * 10 + index + 1),
                        "mean_corrected_puncta_brightness": float(
                            10 ** (condition_index + 1) + index + 1
                        ),
                        "median_punctum_area_um2": 0.1 + index / 100,
                    }
                )
        statistics_rows = []
        for metric in PRESENTATION_BOXPLOT_METRICS:
            statistics_rows.append(
                {
                    "test_family": "omnibus",
                    "metric": metric,
                    "p_value": 0.001,
                    "adjusted_p_value_BH_across_properties": 0.003,
                }
            )
            for first, second in (
                ("before", "FVP_2h"),
                ("before", "recover_2h"),
                ("FVP_2h", "recover_2h"),
            ):
                statistics_rows.append(
                    {
                        "test_family": "pairwise",
                        "metric": metric,
                        "comparison": f"{first}_vs_{second}",
                        "adjusted_p_value_Holm": 0.01,
                    }
                )

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = {
                "puncta_count": root / "count.png",
                "mean_corrected_puncta_brightness": root / "brightness.png",
                "median_punctum_area_um2": root / "area.png",
            }
            captured = []
            import matplotlib.figure

            original_savefig = matplotlib.figure.Figure.savefig

            def inspect_and_save(figure, *args, **kwargs):
                captured.append(figure)
                return original_savefig(figure, *args, **kwargs)

            with patch.object(matplotlib.figure.Figure, "savefig", inspect_and_save):
                generated = make_individual_cell_boxplots(
                    rows,
                    statistics_rows,
                    paths,
                    metrics=PRESENTATION_BOXPLOT_METRICS,
                )

            self.assertEqual(
                generated,
                [paths[metric] for metric in PRESENTATION_BOXPLOT_METRICS],
            )
            self.assertFalse(paths["median_punctum_area_um2"].exists())
            self.assertEqual(len(captured), 2)
            for metric, figure in zip(PRESENTATION_BOXPLOT_METRICS, captured, strict=True):
                np.testing.assert_allclose(
                    figure.get_size_inches(),
                    PRESENTATION_FIGSIZE_INCHES,
                )
                self.assertLess(
                    PRESENTATION_FIGSIZE_INCHES[0] / PRESENTATION_FIGSIZE_INCHES[1],
                    1.5,
                )
                axis = figure.axes[0]
                self.assertEqual(axis.get_title(), "")
                self.assertEqual(axis.get_ylabel(), PRESENTATION_Y_LABELS[metric])
                self.assertEqual(
                    [tick.get_text() for tick in axis.get_xticklabels()],
                    [PRESENTATION_CONDITION_NAMES[condition] for condition in CONDITIONS],
                )
                self.assertTrue(axis.texts)
                self.assertTrue(
                    all(
                        float(text.get_fontsize()) == PRESENTATION_FONT_SIZE
                        for text in axis.texts
                    )
                )
                self.assertEqual(
                    float(axis.yaxis.label.get_fontsize()),
                    PRESENTATION_FONT_SIZE,
                )
                self.assertTrue(
                    all(
                        float(tick.get_fontsize()) == PRESENTATION_FONT_SIZE
                        for tick in axis.get_xticklabels() + axis.get_yticklabels()
                    )
                )
                self.assertEqual(
                    axis.get_yscale(),
                    "log" if metric.endswith("brightness") else "linear",
                )
                self.assertEqual(figure.legends, [])

    def test_cyx_tiff_and_config_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "crop.tif"
            _write_imagej(
                path,
                np.zeros((2, 7, 9), dtype=np.float32),
                axes="CYX",
                pixel_size_um=0.0585,
            )
            with tifffile.TiffFile(path) as tif:
                self.assertEqual(tif.series[0].axes, "CYX")
                self.assertEqual(str(tif.series[0].dtype), "float32")
            config = LiveNuclearPunctaConfig(dataset_root=root)
            config_path = root / "config.json"
            config_path.write_text(json.dumps(config_to_json(config)))
            restored = config_from_json(config_path)
            self.assertEqual(restored, config)
            serialized = config_to_json(config)
            self.assertNotIn("core_radius_fraction", serialized)
            self.assertFalse(hasattr(restored, "core_radius_fraction"))


if __name__ == "__main__":
    unittest.main()
