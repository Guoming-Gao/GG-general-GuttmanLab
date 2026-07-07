from __future__ import annotations

from pathlib import Path
import unittest

import numpy as np
from scipy.io import loadmat
from skimage.metrics import normalized_root_mse

from sacdpy.background import background_estimation
from sacdpy.registration import register_image
from sacdpy.sparse_hessian import sparse_hessian_core


FIXTURES = Path(__file__).resolve().parents[1] / "validation_report_SACDm" / "optional_reference"


def normalized_metrics(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    aa = np.asarray(a, dtype=np.float64)
    bb = np.asarray(b, dtype=np.float64)
    aa = aa / max(float(np.max(np.abs(aa))), 1e-12)
    bb = bb / max(float(np.max(np.abs(bb))), 1e-12)
    corr = float(np.corrcoef(aa.ravel(), bb.ravel())[0, 1])
    nrmse = float(normalized_root_mse(bb, aa))
    mae = float(np.mean(np.abs(aa - bb)))
    return corr, nrmse, mae


class OptionalMatlabParityTests(unittest.TestCase):
    def test_background_estimation_matches_matlab_fixture_with_wavelet_tolerance(self) -> None:
        fixture = loadmat(FIXTURES / "background_fixture.mat")
        result = background_estimation(fixture["bg_input"], decomposition_level=4, iterations=2)
        corr, nrmse, mae = normalized_metrics(result, fixture["bg_output"])
        self.assertGreater(corr, 0.99)
        self.assertLess(nrmse, 0.05)
        self.assertLess(mae, 0.01)

    def test_registration_matches_matlab_fixture(self) -> None:
        fixture = loadmat(FIXTURES / "registration_fixture.mat")
        result = register_image(fixture["reg_moving"], fixture["reg_reference"])
        corr, nrmse, mae = normalized_metrics(result, fixture["reg_output"])
        self.assertGreater(corr, 0.99999)
        self.assertLess(nrmse, 1e-4)
        self.assertLess(mae, 1e-5)

    def test_sparse_hessian_matches_matlab_fixture(self) -> None:
        fixture = loadmat(FIXTURES / "sparse_fixture.mat")
        result = sparse_hessian_core(fixture["sparse_input"], iterations=6)
        corr, nrmse, mae = normalized_metrics(result, fixture["sparse_output"])
        self.assertGreater(corr, 0.99999)
        self.assertLess(nrmse, 1e-5)
        self.assertLess(mae, 1e-6)


if __name__ == "__main__":
    unittest.main()
