from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from sacdpy.batch import find_input_files, sacdpy_output_path, wavelength_for_file


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

    def test_wavelength_for_file_uses_filename_tokens(self) -> None:
        self.assertEqual(
            wavelength_for_file("sample-right_frames_1-50.tif", 561, {"right": 640}),
            640,
        )
        self.assertEqual(wavelength_for_file("sample-other.tif", 561, {"right": 640}), 561)


if __name__ == "__main__":
    unittest.main()
