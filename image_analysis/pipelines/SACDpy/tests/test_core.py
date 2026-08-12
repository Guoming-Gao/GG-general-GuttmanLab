from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from sacdpy.background import background_estimation
from sacdpy.cumulant import cumulant
from sacdpy.fourier import fourier_interpolate
from sacdpy.params import SACDParams
from sacdpy.psf import airy_kernel, generate_rsf
from sacdpy.registration import register_stack
from sacdpy.reconstruction import apply_intensity_transform, as_yxt, reconstruct
from sacdpy.sparse_hessian import sparse_hessian_core
from sacdpy.tiffio import read_tiff_stack, write_tiff_image


class CoreTests(unittest.TestCase):
    def test_as_yxt_converts_tyx_stack(self) -> None:
        stack = np.zeros((5, 8, 9), dtype=np.uint16)
        self.assertEqual(as_yxt(stack).shape, (8, 9, 5))

    def test_generate_rsf_is_normalized(self) -> None:
        rsf = generate_rsf(2.5, 32)
        self.assertEqual(rsf.ndim, 2)
        self.assertTrue(np.isclose(rsf.sum(), 1.0))
        self.assertTrue(np.all(rsf >= 0))

    def test_airy_kernel_is_normalized(self) -> None:
        psf = airy_kernel(117e-9, 560e-9, 1.45, n=428)
        self.assertEqual(psf.shape, (129, 129))
        self.assertTrue(np.isclose(psf.sum(), 1.0))
        self.assertTrue(np.all(psf >= 0))

    def test_fourier_interpolate_lateral_3d_shape(self) -> None:
        rng = np.random.default_rng(0)
        stack = rng.random((8, 6, 3))
        out = fourier_interpolate(stack, (2, 2, 1), mirror_mode="lateral")
        self.assertEqual(out.shape, (16, 12, 3))
        self.assertTrue(np.isfinite(out).all())

    def test_cumulant_order2_matches_matlab_formula(self) -> None:
        stack = np.arange(2 * 3 * 5, dtype=np.float64).reshape(2, 3, 5)
        expected = np.mean(stack[:, :, :-1] * stack[:, :, 1:], axis=2)
        self.assertTrue(np.allclose(cumulant(stack, 2), expected))

    def test_cumulant_order4_matches_matlab_formula(self) -> None:
        stack = np.arange(2 * 3 * 6, dtype=np.float64).reshape(2, 3, 6) + 1
        a = stack[:, :, :-3]
        b = stack[:, :, 1:-2]
        c = stack[:, :, 2:-1]
        d = stack[:, :, 3:]
        expected = (
            np.mean(a * b * c * d, axis=2)
            - np.mean(a * d, axis=2) * np.mean(b * c, axis=2)
            - np.mean(a * c, axis=2) * np.mean(b * d, axis=2)
            - np.mean(a * b, axis=2) * np.mean(c * d, axis=2)
        )
        self.assertTrue(np.allclose(cumulant(stack, 4), expected))

    def test_reconstruct_without_frame_chunk_returns_single_image(self) -> None:
        stack = np.random.default_rng(1).random((6, 5, 50))
        result = reconstruct(stack, self._fast_params())
        self.assertEqual(result.shape, (12, 10))

    def test_reconstruct_with_frame_chunks_returns_tyx_stack(self) -> None:
        stack = np.random.default_rng(2).random((6, 5, 50))
        result = reconstruct(stack, self._fast_params(frames_per_sacd=25))
        self.assertEqual(result.shape, (2, 12, 10))

    def test_order_root_transform_supports_orders_two_through_six(self) -> None:
        values = np.array([[0.0, 1.0, 64.0, 729.0]], dtype=np.float32)
        for order in range(2, 7):
            params = SACDParams(ac_order=order, intensity_transform="order_root")
            transformed = apply_intensity_transform(values, params)
            np.testing.assert_allclose(transformed, values ** (1.0 / order), rtol=1e-6)
            self.assertEqual(transformed.dtype, np.float32)

    def test_raw_cumulant_transform_is_unchanged(self) -> None:
        values = np.array([[0.0, 4.0, 9.0]], dtype=np.float32)
        params = SACDParams(ac_order=2, intensity_transform="raw_cumulant")
        np.testing.assert_array_equal(apply_intensity_transform(values, params), values)

    def test_order_root_uses_ac_order_independently_of_post_psf_scale(self) -> None:
        values = np.array([64.0], dtype=np.float32)
        params = SACDParams(ac_order=3, scale=2, intensity_transform="order_root")
        np.testing.assert_allclose(apply_intensity_transform(values, params), [4.0])

    def test_transform_rejects_negative_and_nonfinite_values(self) -> None:
        params = SACDParams(ac_order=2)
        with self.assertRaisesRegex(ValueError, "finite"):
            apply_intensity_transform(np.array([np.inf]), params)
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            apply_intensity_transform(np.array([-1.0]), params)

    def test_default_reconstruction_is_root_of_raw_for_single_and_chunked(self) -> None:
        stack = np.random.default_rng(12).random((6, 5, 50))
        for frames_per_sacd in (None, 25):
            default_params = self._fast_params(frames_per_sacd=frames_per_sacd)
            raw_params = self._fast_params(frames_per_sacd=frames_per_sacd)
            raw_params.intensity_transform = "raw_cumulant"
            default = reconstruct(stack, default_params)
            raw = reconstruct(stack, raw_params)
            np.testing.assert_allclose(default, np.sqrt(raw), rtol=2e-6, atol=1e-6)

    def test_reconstruct_drops_partial_frame_chunk(self) -> None:
        stack = np.random.default_rng(3).random((6, 5, 55))
        result = reconstruct(stack, self._fast_params(frames_per_sacd=25))
        self.assertEqual(result.shape, (2, 12, 10))

    def test_reconstruct_rejects_invalid_frame_chunks(self) -> None:
        stack = np.random.default_rng(4).random((6, 5, 10))
        with self.assertRaises(ValueError):
            reconstruct(stack, self._fast_params(frames_per_sacd=0))
        with self.assertRaises(ValueError):
            reconstruct(stack, self._fast_params(frames_per_sacd=1))

    def test_background_estimation_preserves_shape_and_finite_values(self) -> None:
        stack = np.random.default_rng(6).random((32, 24, 3)).astype(np.float32)
        background = background_estimation(stack, decomposition_level=2, iterations=1)
        self.assertEqual(background.shape, stack.shape)
        self.assertTrue(np.isfinite(background).all())

    def test_registration_preserves_shape_and_normalizes_registered_frames(self) -> None:
        stack = np.zeros((16, 16, 3), dtype=np.float64)
        stack[6:10, 6:10, 0] = 1
        stack[7:11, 5:9, 1] = 1
        stack[5:9, 7:11, 2] = 1
        registered = register_stack(stack, upsample_factor=10)
        self.assertEqual(registered.shape, stack.shape)
        self.assertTrue(np.isfinite(registered).all())
        self.assertTrue(np.allclose(registered[:, :, 1:].max(axis=(0, 1)), 1.0))

    def test_sparse_hessian_core_accepts_2d_images(self) -> None:
        image = np.random.default_rng(7).random((12, 10)).astype(np.float32)
        out = sparse_hessian_core(image, iterations=2)
        self.assertEqual(out.shape, image.shape)
        self.assertTrue(np.isfinite(out).all())
        self.assertGreaterEqual(float(out.min()), 0.0)

    def test_reconstruct_advanced_branches_execute_on_small_stack(self) -> None:
        stack = np.random.default_rng(8).random((32, 24, 4))
        params = self._fast_params()
        params.ifbackground = True
        params.ifregistration = True
        params.ifsparsedecon = True
        params.sparse_iterations = 2
        result = reconstruct(stack, params)
        self.assertEqual(result.shape, (64, 48))
        self.assertTrue(np.isfinite(result).all())

    def test_cli_help_includes_frames_per_sacd(self) -> None:
        result = subprocess.run(
            [sys.executable, "-m", "sacdpy", "--help"],
            check=True,
            capture_output=True,
            text=True,
            env={
                **os.environ,
                "PYTHONPATH": str(Path(__file__).parents[1] / "src"),
            },
        )
        self.assertIn("--frames-per-sacd", result.stdout)
        self.assertIn("--intensity-transform", result.stdout)

    def test_write_tiff_image_preserves_3d_tyx_shape(self) -> None:
        stack = np.random.default_rng(5).random((2, 8, 9)).astype(np.float32)
        with tempfile.TemporaryDirectory() as tmp:
            path = f"{tmp}/chunked_sacd.tif"
            write_tiff_image(path, stack)
            loaded = read_tiff_stack(path)
        self.assertEqual(loaded.shape, stack.shape)

    @staticmethod
    def _fast_params(frames_per_sacd: int | None = None) -> SACDParams:
        psf = np.ones((3, 3), dtype=np.float64)
        psf /= psf.sum()
        return SACDParams(
            pixel_nm=117.0,
            wavelength_nm=560.0,
            mag=2,
            iter1=0,
            iter2=0,
            ac_order=2,
            frames_per_sacd=frames_per_sacd,
            psf=psf,
        )


if __name__ == "__main__":
    unittest.main()
