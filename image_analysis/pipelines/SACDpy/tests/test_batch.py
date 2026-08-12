from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from sacdpy.batch import (
    find_input_files,
    params_for_file,
    run_batch_reconstruction,
    sacdpy_output_path,
    select_frame_range,
    wavelength_for_file,
)


class BatchTests(unittest.TestCase):
    def test_sacdpy_output_path_places_label_before_channel(self) -> None:
        input_file = Path(
            "/data/TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-left_frames_1-50.tif"
        )

        output = sacdpy_output_path(input_file)

        self.assertEqual(
            output,
            Path("/data/TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-SACDpy-left.tif"),
        )

    def test_sacdpy_output_path_handles_right_channel(self) -> None:
        input_file = Path("/data/sample-right_frames_1-50.tif")

        output = sacdpy_output_path(input_file)

        self.assertEqual(output, Path("/data/sample-SACDpy-right.tif"))

    def test_find_input_files_uses_requested_pattern(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            wanted = folder / "a-left_frames_1-50.tif"
            ignored = folder / "a-left_frames_1-100.tif"
            wanted.touch()
            ignored.touch()

            self.assertEqual(find_input_files(folder), [wanted])

    def test_find_input_files_can_exclude_generated_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            wanted = folder / "sample-DAPI.tif"
            generated = folder / "sample-DAPI-SACDpy.tif"
            wanted.touch()
            generated.touch()

            self.assertEqual(
                find_input_files(folder, "*.tif", exclude_name_contains=("SACDpy",)),
                [wanted],
            )

    def test_wavelength_for_file_uses_filename_tokens(self) -> None:
        self.assertEqual(
            wavelength_for_file("sample-right_frames_1-50.tif", 561, {"right": 640}),
            640,
        )
        self.assertEqual(wavelength_for_file("sample-other.tif", 561, {"right": 640}), 561)

    def test_params_for_file_passes_optional_sacdm_branches(self) -> None:
        params = params_for_file(
            "sample-right_frames_1-50.tif",
            pixel_nm=117.0,
            na=1.45,
            default_wavelength_nm=561,
            wavelength_by_name={"right": 640},
            ifbackground=True,
            backgroundfactor=3.0,
            ifregistration=True,
            ifsparsedecon=True,
            fidelity=50.0,
            tcontinuity=0.25,
            sparsity=2.0,
            sparse_iterations=5,
        )

        self.assertEqual(params.wavelength_nm, 640)
        self.assertTrue(params.ifbackground)
        self.assertEqual(params.backgroundfactor, 3.0)
        self.assertTrue(params.ifregistration)
        self.assertTrue(params.ifsparsedecon)
        self.assertEqual(params.fidelity, 50.0)
        self.assertEqual(params.tcontinuity, 0.25)
        self.assertEqual(params.sparsity, 2.0)
        self.assertEqual(params.sparse_iterations, 5)

    def test_select_frame_range_uses_one_based_inclusive_bounds_for_tyx(self) -> None:
        frame_values = np.arange(500, dtype=np.uint16)[:, None, None]
        stack = np.broadcast_to(frame_values, (500, 684, 428))

        selected = select_frame_range(stack, frame_start=1, frame_end=25)

        self.assertEqual(selected.shape, (25, 684, 428))
        np.testing.assert_array_equal(selected, stack[0:25])

    def test_select_frame_range_supports_alternate_and_all_frame_ranges(self) -> None:
        stack = np.arange(10 * 4 * 5).reshape(10, 4, 5)

        selected = select_frame_range(stack, frame_start=2, frame_end=4)

        np.testing.assert_array_equal(selected, stack[1:4])
        np.testing.assert_array_equal(select_frame_range(stack), stack)

    def test_select_frame_range_rejects_invalid_bounds(self) -> None:
        stack = np.zeros((10, 4, 5))

        with self.assertRaisesRegex(ValueError, "frame_start must be >= 1"):
            select_frame_range(stack, frame_start=0, frame_end=5)
        with self.assertRaisesRegex(ValueError, "frame_end must be >= frame_start"):
            select_frame_range(stack, frame_start=6, frame_end=5)

    def test_select_frame_range_rejects_frames_beyond_input(self) -> None:
        stack = np.zeros((10, 4, 5))

        with self.assertRaisesRegex(ValueError, "contains only 10 frame"):
            select_frame_range(stack, frame_start=1, frame_end=25)

    def test_select_frame_range_requires_a_tyx_stack(self) -> None:
        with self.assertRaisesRegex(ValueError, "3D TYX"):
            select_frame_range(np.zeros((4, 5)), frame_start=1, frame_end=1)

    def test_run_batch_passes_selected_frames_to_reconstruct_as_yxt(self) -> None:
        input_file = Path("sample-cropped-left.tif")
        stack = np.arange(6 * 8 * 10, dtype=np.uint16).reshape(6, 8, 10)
        expected_yxt = np.moveaxis(stack[1:4], 0, -1)

        with (
            tempfile.TemporaryDirectory() as output_dir,
            patch("sacdpy.batch.read_tiff_stack", return_value=stack),
            patch("sacdpy.batch.reconstruct", return_value=np.zeros((16, 20))) as reconstruct_mock,
            patch("sacdpy.batch.write_tiff_image"),
        ):
            results = run_batch_reconstruction(
                [input_file],
                output_dir=output_dir,
                pixel_nm=117.0,
                na=1.45,
                default_wavelength_nm=561,
                wavelength_by_name={"cropped-left": 549},
                frame_start=2,
                frame_end=4,
            )

        received_yxt = reconstruct_mock.call_args.args[0]
        self.assertEqual(received_yxt.shape, (8, 10, 3))
        np.testing.assert_array_equal(received_yxt, expected_yxt)
        self.assertEqual(results[0].input_shape, (3, 8, 10))
        self.assertEqual(results[0].wavelength_nm, 549)


if __name__ == "__main__":
    unittest.main()
