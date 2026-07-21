from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import tifffile

from sacdpy.phase_diagram import (
    PhaseDiagramConfig,
    PhaseFOVPlan,
    _export_nucleus_crops,
    _stage_z_summary,
    calculate_dispersion_scores,
    discover_phase_fovs,
    get_bbox_with_padding,
    make_dispersion_phase_diagram,
    make_intensity_core,
    nucleus_review_stem,
    parse_condition,
    parse_fov_token,
    parse_review_condition,
    preprocess_cellpose_image,
    run_cellpose_nuclei,
    split_dual_view,
)


class FakeCellposeModel:
    def eval(self, image, **kwargs):
        mask = np.zeros(image.shape[-2:], dtype=np.uint16)
        mask[2:6, 3:8] = 1
        return mask, None, None


class PhaseDiagramTests(unittest.TestCase):
    def test_parse_condition_keeps_biological_prefix(self) -> None:
        self.assertEqual(
            parse_condition("dSPEN_dRRM-p1x-Hoechst-SPEN_JFX650-200ms-FOV-mysterycell"),
            "dSPEN_dRRM-p1x",
        )

    def test_review_condition_and_filename_are_pooled_but_keep_dose(self) -> None:
        plan = PhaseFOVPlan(
            Path("."),
            "dSPEN_FL-1x-Hoechst-SPEN_JFX650-200ms-FOV-1",
            "dSPEN_FL-1x",
            "prefix",
            0,
            {},
        )
        self.assertEqual(parse_review_condition(plan.condition), ("dSPEN_FL", "1x"))
        self.assertEqual(parse_fov_token(plan.fov_name), "FOV-1")
        self.assertEqual(
            parse_fov_token("dSPEN_dRRM-p1x-Hoechst-SPEN_JFX650-200ms-FOV-cytoplasmic"),
            "FOV-cytoplasmic",
        )
        self.assertEqual(
            nucleus_review_stem(plan, 1, 14370),
            "14370-dSPEN_FL-1x-FOV-1-Hoechst_SPEN_mask-nucleus-0001",
        )

    def test_discover_phase_fovs_requires_single_t0_contiguous_z_grid(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pos = root / "dSPEN_FL-1x-Hoechst-SPEN_JFX650-200ms-FOV" / "pos_0"
            pos.mkdir(parents=True)
            for z in range(3):
                (pos / f"sample_posXY0_channels_t0_posZ{z}.tif").touch()

            plans = discover_phase_fovs(root)

        self.assertEqual(len(plans), 1)
        self.assertEqual(plans[0].z_indices, (0, 1, 2))
        self.assertEqual(plans[0].condition, "dSPEN_FL-1x")

    def test_split_dual_view(self) -> None:
        movie = np.arange(2 * 3 * 8).reshape(2, 3, 8)
        left, right = split_dual_view(movie)
        self.assertEqual(left.shape, (2, 3, 4))
        self.assertEqual(right.shape, (2, 3, 4))
        np.testing.assert_array_equal(left, movie[..., :4])
        np.testing.assert_array_equal(right, movie[..., 4:])
        with self.assertRaisesRegex(ValueError, "even"):
            split_dual_view(np.zeros((2, 3, 7)))

    def test_bbox_padding_clips_to_boundaries(self) -> None:
        mask = np.zeros((10, 12), dtype=bool)
        mask[0:3, 8:12] = True
        self.assertEqual(get_bbox_with_padding(mask, 2), (0, 5, 6, 12))

    def test_equivalent_radius_core_erodes_by_one_third_radius(self) -> None:
        yy, xx = np.ogrid[:101, :101]
        mask = (yy - 50) ** 2 + (xx - 50) ** 2 <= 30**2
        core, equivalent_radius, erosion_distance, max_inscribed = make_intensity_core(mask)
        self.assertAlmostEqual(erosion_distance, equivalent_radius / 3.0)
        self.assertAlmostEqual(max_inscribed, 30.0, delta=1.1)
        self.assertAlmostEqual(core.sum() / mask.sum(), 4.0 / 9.0, delta=0.04)

    def test_dispersion_scores_are_scale_independent(self) -> None:
        values = np.array([1, 1, 2, 6], dtype=np.float64)
        first = calculate_dispersion_scores(values)
        second = calculate_dispersion_scores(values * 100)
        np.testing.assert_allclose(first, second)

    def test_percentile_preprocessing_is_bounded(self) -> None:
        image = np.arange(100, dtype=np.float32).reshape(10, 10)
        result = preprocess_cellpose_image(image, 1, 99)
        self.assertEqual(result.dtype, np.float32)
        self.assertGreaterEqual(float(result.min()), 0)
        self.assertLessEqual(float(result.max()), 1)

    def test_cellpose_wrapper_returns_label_image(self) -> None:
        config = PhaseDiagramConfig()
        labels = run_cellpose_nuclei(np.ones((10, 12), dtype=np.float32), FakeCellposeModel(), config)
        self.assertEqual(labels.shape, (10, 12))
        self.assertEqual(labels.dtype, np.uint16)
        self.assertEqual(int(labels.max()), 1)

    def test_stage_z_summary_flags_nonmonotonic_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            paths = []
            for index, z_value in enumerate((0.0, -0.7, 0.0, -2.1)):
                path = Path(tmp) / f"z{index}.tif"
                tifffile.imwrite(path, np.zeros((3, 4), dtype=np.uint16), description=f'{{"StagePos_um": [0, 0, {z_value}]}}')
                paths.append(path)
            positions, spacing, warnings = _stage_z_summary(paths)

        self.assertEqual(positions, [0.0, -0.7, 0.0, -2.1])
        self.assertAlmostEqual(spacing, 0.7)
        self.assertIn("non_monotonic_stage_z", warnings)

    def test_export_crops_writes_combined_review_files_and_masked_metrics(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            plan = PhaseFOVPlan(
                root,
                "dSPEN_dRRM-p1x-Hoechst-SPEN_JFX650-200ms-FOV",
                "dSPEN_dRRM-p1x",
                "prefix",
                0,
                {0: root / "z0.tif"},
            )
            config = PhaseDiagramConfig(output_root=root, crop_padding_px=2, mag=2)
            labels = np.zeros((12, 16), dtype=np.uint16)
            labels[3:8, 5:11] = 1
            hoechst = np.arange(labels.size, dtype=np.float32).reshape(labels.shape)
            spen = np.full(labels.shape, 1_000_000, dtype=np.float32)
            spen[labels == 1] = np.arange(30, dtype=np.float32) + 100

            records, excluded = _export_nucleus_crops(
                plan,
                config,
                root / "crops",
                labels,
                hoechst,
                spen,
                pixel_size_um=0.0585,
            )

            self.assertEqual(len(records), 1)
            self.assertEqual(excluded, [])
            review = tifffile.imread(records[0]["review_tif"])
            self.assertEqual(review.shape, (3, 9, 10))
            self.assertEqual(review.dtype, np.float32)
            np.testing.assert_array_equal(review[0], hoechst[1:10, 3:13])
            np.testing.assert_array_equal(review[1], spen[1:10, 3:13])
            self.assertEqual(set(np.unique(review[2])), {0.0, 1.0})
            self.assertTrue(Path(records[0]["review_png"]).is_file())
            full_mask = labels == 1
            core, _, _, _ = make_intensity_core(full_mask)
            expected = spen[core].astype(np.float64)
            self.assertAlmostEqual(records[0]["spen_sacd_core_mean"], float(expected.mean()))
            self.assertAlmostEqual(records[0]["spen_sacd_core_max"], float(expected.max()))
            self.assertEqual(records[0]["spen_sacd_core_mean_rounded"], int(np.rint(expected.mean())))
            self.assertGreater(records[0]["spen_sacd_core_gini"], 0)
            self.assertEqual(records[0]["pooled_condition"], "dSPEN_dRRM")
            self.assertEqual(records[0]["dose"], "p1x")
            self.assertEqual(records[0]["bbox_y0_sacd_px"], 1)
            self.assertEqual(records[0]["bbox_x0_sacd_px"], 3)

    def test_empty_and_nonpositive_cores_are_excluded(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            plan = PhaseFOVPlan(root, "dSPEN_FL-1x-Hoechst-SPEN-FOV", "dSPEN_FL-1x", "", 0, {})
            labels = np.zeros((40, 50), dtype=np.uint16)
            labels[5, 5:35] = 1
            labels[15:35, 10:30] = 2
            spen = np.ones(labels.shape, dtype=np.float32)
            spen[labels == 2] = 0
            records, excluded = _export_nucleus_crops(
                plan, PhaseDiagramConfig(crop_padding_px=2), root / "crops", labels, spen, spen
            )
            self.assertEqual(records, [])
            self.assertEqual(
                {row["exclusion_reason"] for row in excluded},
                {"empty_core_after_equivalent_radius_erosion", "nonpositive_core_mean"},
            )

    def test_dispersion_phase_diagram_has_three_score_panels(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "phase.png"
            records = [
                {
                    "pooled_condition": "dSPEN_FL",
                    "spen_sacd_core_mean": 100.0,
                    "spen_sacd_core_gini": 0.2,
                    "spen_sacd_core_coefficient_of_variation": 0.5,
                    "spen_sacd_core_robust_percentile_dispersion": 0.8,
                },
                {
                    "pooled_condition": "dSPEN_dRRM",
                    "spen_sacd_core_mean": 200.0,
                    "spen_sacd_core_gini": 0.3,
                    "spen_sacd_core_coefficient_of_variation": 0.7,
                    "spen_sacd_core_robust_percentile_dispersion": 1.1,
                },
            ]
            make_dispersion_phase_diagram(records, output)
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 1_000)


if __name__ == "__main__":
    unittest.main()
