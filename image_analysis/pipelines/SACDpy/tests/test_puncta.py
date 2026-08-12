from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from sacdpy.puncta import (
    SPOTIFLOW_GENERAL_WEIGHT_SHA256,
    SpotiflowHybridParams,
    _build_puncta_phase_figure,
    _write_json,
    construct_spotiflow_hybrid_mask,
    detector_environment_manifest,
    group_spotiflow_seeds,
    make_puncta_phase_diagram,
    quantify_nucleus_puncta,
    robust_local_background,
    spotiflow_seed_table,
)


def seed(seed_id: int, y: float, x: float) -> dict[str, float | int]:
    return {
        "seed_id": seed_id,
        "source_seed_index": seed_id - 1,
        "y_sacd_px": y,
        "x_sacd_px": x,
        "rounded_y_sacd_px": int(round(y)),
        "rounded_x_sacd_px": int(round(x)),
        "spotiflow_probability": 0.8,
        "spotiflow_model_intensity": 1.0,
    }


class FakeSpotiflow:
    def __init__(self) -> None:
        self.kwargs: dict[str, object] = {}

    def predict(self, image: np.ndarray, **kwargs: object):
        self.kwargs = kwargs
        points = np.asarray([[10.0, 10.0], [2.0, 2.0], [50.0, 50.0]])
        details = SimpleNamespace(
            prob=np.asarray([0.9, 0.8, 0.7]),
            intens=np.asarray([3.0, 2.0, 1.0]),
        )
        return points, details


class SpotiflowHybridTests(unittest.TestCase):
    def test_contract_and_physical_pixel_conversion(self) -> None:
        params = SpotiflowHybridParams()
        params.validate()
        pixel_size = 0.0585
        self.assertAlmostEqual(params.seed_radius_um / pixel_size, 4.0)
        self.assertAlmostEqual(params.grouped_link_distance_um / pixel_size, 8.0)
        self.assertAlmostEqual(
            params.growth_smoothing_sigma_um / pixel_size,
            2.0,
        )
        self.assertAlmostEqual(
            params.maximum_growth_distance_um / pixel_size,
            20.0,
        )
        manifest = detector_environment_manifest()
        self.assertFalse(manifest["custom_training"])
        self.assertEqual(
            manifest["model_weight_sha256"],
            SPOTIFLOW_GENERAL_WEIGHT_SHA256,
        )

    def test_spotiflow_defaults_and_core_center_filter(self) -> None:
        image = np.ones((32, 32), dtype=np.float32)
        core = np.zeros_like(image, dtype=bool)
        core[5:20, 5:20] = True
        model = FakeSpotiflow()
        rows = spotiflow_seed_table(
            image,
            core,
            model,
            params=SpotiflowHybridParams(),
        )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["seed_id"], 1)
        self.assertEqual(model.kwargs["prob_thresh"], 0.5)
        self.assertEqual(model.kwargs["min_distance"], 1)
        self.assertFalse(model.kwargs["exclude_border"])
        self.assertEqual(model.kwargs["normalizer"], "auto")
        self.assertIs(model.kwargs["device"], "auto")

    def test_grouping_is_deterministic_and_transitive(self) -> None:
        rows = [seed(1, 10, 10), seed(2, 10, 17), seed(3, 10, 24), seed(4, 40, 40)]
        groups = group_spotiflow_seeds(rows, link_distance_px=8.0)
        self.assertEqual(groups, [[0, 1, 2], [3]])

    def test_growth_and_fallback_are_core_confined(self) -> None:
        yy, xx = np.mgrid[:80, :80]
        core = (yy - 40) ** 2 + (xx - 40) ** 2 <= 30**2
        image = np.full(core.shape, 10.0, dtype=np.float32)
        image += (
            150.0
            * np.exp(-((yy - 40) ** 2 + (xx - 40) ** 2) / (2 * 5.0**2))
        ).astype(np.float32)
        rows = [seed(1, 40, 40)]
        labels, groups, provenance = construct_spotiflow_hybrid_mask(
            image,
            core,
            rows,
            pixel_size_um=0.0585,
            params=SpotiflowHybridParams(),
        )
        self.assertEqual(labels.dtype, np.uint32)
        self.assertFalse(np.any(labels[~core]))
        self.assertEqual(groups[0]["growth_status"], "grown")
        base_area = np.pi * (0.234 / 0.0585) ** 2
        self.assertGreater(np.sum(labels > 0), base_area)
        self.assertEqual(provenance[1]["member_seed_count"], 1)

        flat = np.full(core.shape, 10.0, dtype=np.float32)
        fallback_labels, fallback_groups, _ = construct_spotiflow_hybrid_mask(
            flat,
            core,
            [seed(1, 40, 40)],
            pixel_size_um=0.0585,
            params=SpotiflowHybridParams(),
        )
        self.assertEqual(
            fallback_groups[0]["growth_status"],
            "fallback_seed_below_local_threshold",
        )
        self.assertGreater(np.sum(fallback_labels > 0), 0)
        self.assertFalse(np.any(fallback_labels[~core]))

    def test_overlapping_seed_footprints_become_one_punctum(self) -> None:
        core = np.ones((64, 64), dtype=bool)
        image = np.full((64, 64), 10.0, dtype=np.float32)
        rows = [seed(1, 30, 30), seed(2, 30, 37)]
        labels, groups, provenance = construct_spotiflow_hybrid_mask(
            image,
            core,
            rows,
            pixel_size_um=0.0585,
            params=SpotiflowHybridParams(),
        )
        self.assertEqual(len(groups), 1)
        self.assertEqual(int(labels.max()), 1)
        self.assertEqual(provenance[1]["member_seed_count"], 2)

    def test_robust_background_clips_bright_outlier(self) -> None:
        values = np.asarray([10.0] * 100 + [1_000_000.0], dtype=np.float32)
        background, noise = robust_local_background(values)
        self.assertAlmostEqual(background, 10.0)
        self.assertGreater(noise, 0)

    def test_quantification_and_zero_puncta_semantics(self) -> None:
        image = np.full((8, 8), 10.0)
        image[2:4, 2:4] = np.asarray([[20.0, 30.0], [40.0, 50.0]])
        core = np.ones_like(image, dtype=bool)
        labels = np.zeros_like(image, dtype=np.uint32)
        labels[2:4, 2:4] = 1
        summary, rows = quantify_nucleus_puncta(
            image,
            core,
            labels,
            pixel_size_um=0.1,
            provenance={
                1: {
                    "member_seed_count": 2,
                    "seed_ids": [1, 2],
                    "merged_group_count": 1,
                    "group_ids": [1],
                    "growth_status": "grown",
                }
            },
        )
        self.assertEqual(summary["puncta_count"], 1)
        self.assertAlmostEqual(summary["puncta_pixel_mean_sacd"], 35.0)
        self.assertAlmostEqual(summary["nonpuncta_core_mean_sacd"], 10.0)
        self.assertAlmostEqual(
            summary["max_background_corrected_integrated_punctum_sacd"],
            100.0,
        )
        self.assertEqual(rows[0]["member_seed_count"], 2)

        empty, empty_rows = quantify_nucleus_puncta(
            image,
            core,
            np.zeros_like(labels),
            pixel_size_um=0.1,
        )
        self.assertEqual(empty["puncta_count"], 0)
        self.assertEqual(empty["max_raw_integrated_punctum_sacd"], 0.0)
        self.assertTrue(np.isnan(empty["puncta_pixel_mean_sacd"]))
        self.assertAlmostEqual(empty["nonpuncta_core_mean_sacd"], image.mean())
        self.assertEqual(empty_rows, [])

    def test_phase_diagram_is_created(self) -> None:
        rows = []
        for condition in ("dSPEN_FL", "dSPEN_dRRM"):
            for index in range(3):
                rows.append(
                    {
                        "pooled_condition": condition,
                        "spen_sacd_core_mean": 100.0 * (index + 1),
                        "puncta_count": index,
                        "puncta_pixel_mean_sacd": (
                            float("nan") if index == 0 else 500.0
                        ),
                        "nonpuncta_core_mean_sacd": 80.0,
                        "max_background_corrected_integrated_punctum_sacd": 1000.0,
                    }
                )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "phase.png"
            make_puncta_phase_diagram(rows, path)
            self.assertTrue(path.is_file())
            self.assertGreater(path.stat().st_size, 1000)

    def test_phase_diagram_facets_scales_counts_and_legend(self) -> None:
        import matplotlib.pyplot as plt

        rows = []
        specifications = (
            ("dSPEN_FL", 932, 410),
            ("dSPEN_dRRM", 299, 173),
        )
        for condition, count, zero_count in specifications:
            for index in range(count):
                has_puncta = index >= zero_count
                rows.append(
                    {
                        "pooled_condition": condition,
                        "spen_sacd_core_mean": 2_000.0 * 1.01**index,
                        "puncta_count": (index % 41) + 1 if has_puncta else 0,
                        "puncta_pixel_mean_sacd": (
                            5_000.0 * 1.008**index
                            if has_puncta
                            else float("nan")
                        ),
                        "nonpuncta_core_mean_sacd": 1_500.0 * 1.009**index,
                    }
                )
        fig, axes = _build_puncta_phase_figure(rows)
        try:
            self.assertEqual(axes.shape, (2, 3))
            self.assertEqual(len(fig.axes), 6)
            for ax in axes.flat:
                self.assertEqual(ax.get_xscale(), "log")
            self.assertEqual(
                [[ax.get_yscale() for ax in row] for row in axes],
                [["linear", "log", "log"], ["linear", "log", "log"]],
            )
            x_limits = [ax.get_xlim() for ax in axes.flat]
            self.assertTrue(all(limits == x_limits[0] for limits in x_limits))
            for column in range(3):
                self.assertEqual(
                    axes[0, column].get_ylim(),
                    axes[1, column].get_ylim(),
                )
            point_counts = [
                [
                    len(axes[row, column].collections[0].get_offsets())
                    for column in range(3)
                ]
                for row in range(2)
            ]
            self.assertEqual(
                point_counts,
                [[932, 522, 932], [299, 126, 299]],
            )
            self.assertEqual(
                [axes[0, column].get_title() for column in range(3)],
                [
                    "Puncta count",
                    "Mean puncta-pixel intensity (a.u.)",
                    "Mean non-puncta core intensity (a.u.)",
                ],
            )
            self.assertEqual(len(fig.legends), 1)
            self.assertEqual(fig.legends[0]._ncols, 2)
            self.assertEqual(fig.legends[0]._loc, 8)
        finally:
            plt.close(fig)

    def test_checkpoint_json_maps_nan_to_null(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "checkpoint.json"
            _write_json(path, {"defined": 1.0, "undefined": float("nan")})
            value = json.loads(path.read_text())
            self.assertEqual(value["defined"], 1.0)
            self.assertIsNone(value["undefined"])


if __name__ == "__main__":
    unittest.main()
