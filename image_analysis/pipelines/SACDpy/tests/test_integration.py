from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import unittest

from sacdpy import SACDParams, reconstruct


DATA = Path(__file__).resolve().parent / "testdata"


class IntegrationTests(unittest.TestCase):
    def test_reconstruct_matches_fixture_shape_and_similarity(self) -> None:
        cases = [
        (
            "TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-left.tif",
            "TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-SACD-left.tif",
            560.0,
        ),
        (
            "TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-right.tif",
            "TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-SACD-right.tif",
            647.0,
        ),
        ]
        for raw_name, ref_name, wavelength in cases:
            with self.subTest(raw_name=raw_name):
                raw = tifffile.imread(DATA / raw_name)
                ref = tifffile.imread(DATA / ref_name).astype(np.float64)
                result = reconstruct(
                    raw,
                    SACDParams(pixel_nm=117.0, wavelength_nm=wavelength, na=1.45, iter1=1, iter2=1),
                ).astype(np.float64)

                self.assertEqual(result.shape, ref.shape)
                self.assertTrue(np.isfinite(result).all())
                self.assertGreaterEqual(float(result.min()), 0.0)

                result_norm = result / max(float(result.max()), 1.0)
                ref_norm = ref / max(float(ref.max()), 1.0)
                corr = np.corrcoef(result_norm.ravel(), ref_norm.ravel())[0, 1]
                self.assertGreater(corr, 0.25)


if __name__ == "__main__":
    unittest.main()
