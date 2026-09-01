import copy
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

import BulkFluoRDF_granuleSACD as granule


class GranuleSACDTests(unittest.TestCase):
    def test_synthetic_acceptance_rules(self):
        granule.synthetic_validation()

    def test_four_channel_pairing_and_crop_order(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for index, channel in enumerate(granule.CHANNELS):
                tifffile.imwrite(root / f"sample__{channel}-SACD-MIP-YX.tif", np.full((20, 20), index, np.float32))
            pair = granule.pair_four_channel_files(root)[0]
            images = granule.read_fov(pair)
            props = granule.pd.DataFrame([{
                "keep_granule": True, "granule_id": 1, "bbox_min_y": 5, "bbox_min_x": 6,
                "bbox_max_y": 10, "bbox_max_x": 12,
            }])
            paths = granule.export_granule_crops(pair, images, props, root / "out", 3)
            crop = tifffile.imread(paths[0])
            self.assertEqual(crop.shape, (4, 11, 12))
            self.assertEqual([float(crop[i, 0, 0]) for i in range(4)], [0, 1, 2, 3])

    def test_size_normalized_rdf_alignment(self):
        cfg = copy.deepcopy(granule.DEFAULT_CONFIG)
        cfg["rdf"].update({"maximum_normalized_radius": 1.3, "bin_width": 0.2, "bin_step": 0.1})
        shape = (160, 160)
        yy, xx = np.indices(shape)
        images = {ch: np.zeros(shape, np.float32) for ch in granule.CHANNELS}
        labels = np.zeros(shape, np.uint16)
        records = []
        for gid, (cy, cx, radius) in enumerate(((45, 45, 10), (110, 110, 20)), 1):
            rho = np.hypot(yy - cy, xx - cx) / radius
            mask = rho <= 1
            labels[mask] = gid
            profile = np.maximum(0, 2 - rho)
            for channel in granule.RNA_CHANNELS:
                images[channel] += profile.astype(np.float32)
            area = int(mask.sum())
            records.append({"keep_granule": True, "granule_id": gid, "centroid_y_px": cy,
                            "centroid_x_px": cx, "equivalent_radius_px": np.sqrt(area / np.pi)})
        props = granule.pd.DataFrame(records)
        rdf = granule.calculate_granule_rdf("synthetic", images, labels, props, cfg)
        curves = [g.rdf_normalized.to_numpy() for _, g in rdf[rdf.channel == "488"].groupby("granule_id")]
        self.assertLess(float(np.nanmean(np.abs(curves[0] - curves[1]))), 0.08)

    def test_wide_profiles_and_correlations(self):
        radius = np.arange(25, dtype=float)
        profiles = {"488": radius, "561": 2 * radius + 5, "647": -radius + 30}
        rows = []
        for channel, values in profiles.items():
            for index, value in enumerate(values):
                rows.append({
                    "fov": "fov", "granule_id": 1, "channel": channel,
                    "radius_start_r_over_R": index * 0.05, "radius_end_r_over_R": index * 0.05 + 0.1,
                    "radius_mid_r_over_R": index * 0.05 + 0.05, "effective_pixel_area": 10.0,
                    "inside_granule_mean": 2.0, "annular_mean": value * 2 + 10,
                    "rdf_normalized": value,
                })
        rdf = pd.DataFrame(rows)
        wide = granule.rdf_to_wide(rdf)
        props = pd.DataFrame([{"fov": "fov", "granule_id": 1, "keep_granule": True,
                               "equivalent_diameter_px": 20.0, "equivalent_diameter_nm": 1170.0}])
        correlations = granule.calculate_rdf_correlations(wide, props)
        self.assertEqual(len(wide), 25)
        self.assertAlmostEqual(correlations.iloc[0].pearson_r_488_561, 1.0)
        self.assertAlmostEqual(correlations.iloc[0].pearson_r_488_647, -1.0)
        raw_r = np.corrcoef(wide.annular_mean_488, wide.annular_mean_561)[0, 1]
        self.assertAlmostEqual(raw_r, correlations.iloc[0].pearson_r_488_561)

    def test_representative_extremes_and_duplicate_slots(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "granule_crops").mkdir()
            correlation_rows, wide_rows = [], []
            values = [
                ("a", 1, 0.9, 0.2, -0.3),
                ("b", 2, -0.8, 0.8, 0.1),
                ("c", 3, 0.1, -0.7, 0.95),
            ]
            for fov, gid, r1, r2, r3 in values:
                correlation_rows.append({"fov": fov, "granule_id": gid, "equivalent_diameter_px": 20,
                                         "equivalent_diameter_nm": 1170, "n_bins": 25,
                                         "pearson_r_488_561": r1, "pearson_r_488_647": r2,
                                         "pearson_r_561_647": r3})
                tifffile.imwrite(root / "granule_crops" / f"{fov}__granule-{gid:03d}.tif",
                                 np.zeros((4, 12, 12), np.float32), imagej=True,
                                 metadata={"axes": "CYX", "Labels": list(granule.CHANNELS)})
                for index in range(25):
                    wide_rows.append({"fov": fov, "granule_id": gid, "radius_mid_r_over_R": index * 0.05 + 0.05,
                                      **{f"rdf_normalized_{channel}": index / 25 for channel in granule.RNA_CHANNELS}})
            index = granule.export_correlation_representatives(
                pd.DataFrame(correlation_rows), pd.DataFrame(wide_rows), root
            )
            self.assertEqual(len(index), 6)
            self.assertEqual(index.set_index("selection_slot").loc["max_488-561", "fov"], "a")
            self.assertEqual(index.set_index("selection_slot").loc["min_488-561", "fov"], "b")
            for row in index.itertuples():
                self.assertTrue((root / row.representative_tif).exists())
                self.assertTrue((root / row.representative_plot).exists())
            self.assertEqual(len(list((root / "correlation_representatives").glob("*/*.tif"))), 6)
            self.assertEqual(len(list((root / "correlation_representatives").glob("*/*.png"))), 6)


if __name__ == "__main__":
    unittest.main()
