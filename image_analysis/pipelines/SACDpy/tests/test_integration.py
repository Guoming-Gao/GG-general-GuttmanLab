from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile
import unittest
from skimage.metrics import normalized_root_mse, structural_similarity

from sacdpy import SACDParams, reconstruct


DATA = Path(__file__).resolve().parent / "testdata"
MATLAB_REF = DATA.parents[1] / "validation_report_SACDm" / "matlab_reference"


class IntegrationTests(unittest.TestCase):
    def test_reconstruct_matches_matlab_sacdm_reference(self) -> None:
        cases = [
        (
            "TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-left.tif",
            MATLAB_REF / "sacdm-left.tif",
            560.0,
            0.999,
            0.999,
            0.03,
        ),
        (
            "TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-right.tif",
            MATLAB_REF / "sacdm-right.tif",
            647.0,
            0.9999,
            0.9999,
            0.005,
        ),
        ]
        for raw_name, ref_path, wavelength, min_corr, min_ssim, max_nrmse in cases:
            with self.subTest(raw_name=raw_name):
                raw = tifffile.imread(DATA / raw_name)
                ref = tifffile.imread(ref_path).astype(np.float64)
                result = reconstruct(
                    raw,
                    SACDParams(pixel_nm=117.0, wavelength_nm=wavelength, na=1.45),
                ).astype(np.float64)

                self.assertEqual(result.shape, ref.shape)
                self.assertTrue(np.isfinite(result).all())
                self.assertGreaterEqual(float(result.min()), 0.0)

                result_norm = result / max(float(result.max()), 1.0)
                ref_norm = ref / max(float(ref.max()), 1.0)
                corr = np.corrcoef(result_norm.ravel(), ref_norm.ravel())[0, 1]
                data_range = max(float(ref_norm.max() - ref_norm.min()), 1e-12)
                ssim = structural_similarity(ref_norm, result_norm, data_range=data_range)
                nrmse = normalized_root_mse(ref_norm, result_norm)
                self.assertGreater(corr, min_corr)
                self.assertGreater(ssim, min_ssim)
                self.assertLess(nrmse, max_nrmse)


if __name__ == "__main__":
    unittest.main()
