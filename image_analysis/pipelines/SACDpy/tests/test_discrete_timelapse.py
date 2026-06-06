from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import tifffile

from sacdpy.discrete_timelapse import (
    discover_folder_plans,
    find_timelapse_files,
    group_timelapse_files,
    infer_na,
    infer_pixel_nm,
    infer_time_interval_s,
    infer_z_spacing_um,
    parse_time_interval_from_path,
    parse_timelapse_filename,
    read_oni_metadata,
    write_timelapse_tiff,
)


class DiscreteTimelapseTests(unittest.TestCase):
    def test_parse_timelapse_filename_uses_observed_oni_grammar(self) -> None:
        item = parse_timelapse_filename(
            "ONI-gmgao-SPEN_SACDlive_sample_posXY0_channels_t14_posZ4.tif",
            channel_name="647",
        )

        self.assertEqual(item.prefix, "ONI-gmgao-SPEN_SACDlive_sample")
        self.assertEqual(item.channel, "647")
        self.assertEqual(item.pos_xy, 0)
        self.assertEqual(item.time_index, 14)
        self.assertEqual(item.z_index, 4)

    def test_group_timelapse_files_reports_complete_and_missing_grids(self) -> None:
        files = [
            parse_timelapse_filename(f"sample_posXY0_channels_t{time}_posZ{z}.tif")
            for time in range(2)
            for z in range(3)
        ]
        complete = group_timelapse_files(files)[0]
        self.assertTrue(complete.is_complete)
        self.assertEqual(complete.time_indices, (0, 1))
        self.assertEqual(complete.z_indices, (0, 1, 2))

        missing = group_timelapse_files(files[:-1])[0]
        self.assertFalse(missing.is_complete)
        self.assertEqual(missing.missing_entries, ((1, 2),))

    def test_find_timelapse_files_recurses_and_ignores_unmatched_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            pos = folder / "pos_0"
            pos.mkdir()
            wanted = pos / "sample_posXY0_channels_t0_posZ0.tif"
            wanted.touch()
            (pos / "sample-SACDpy-647-TZYX.tif").touch()
            (pos / "notes.txt").touch()

            found = find_timelapse_files(folder)

        self.assertEqual([item.path for item in found], [wanted])

    def test_discover_folder_plans_keeps_one_plan_per_fov_folder(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fov_a = root / "fov_a"
            fov_b = root / "fov_b"
            fov_a.mkdir()
            fov_b.mkdir()
            (fov_a / "sampleA_posXY0_channels_t0_posZ0.tif").touch()
            (fov_b / "sampleB_posXY0_channels_t0_posZ0.tif").touch()

            plans = discover_folder_plans([fov_a, fov_b], output_dir=root / "out")

        self.assertEqual([plan.input_folder.name for plan in plans], ["fov_a", "fov_b"])
        self.assertEqual([len(plan.files) for plan in plans], [1, 1])
        self.assertEqual([len(plan.groups) for plan in plans], [1, 1])
        self.assertEqual({plan.output_dir.name for plan in plans}, {"out"})

    def test_read_oni_metadata_and_infer_scalar_values(self) -> None:
        metadata = {
            "PixelSize_um": 0.117,
            "Objective_NA": 1.45,
            "timestamp_us": 1_000_000,
            "StagePos_um": [0.0, 0.0, -5.0],
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "meta.tif"
            tifffile.imwrite(path, np.zeros((2, 3), dtype=np.uint16), description=json.dumps(metadata))

            self.assertEqual(read_oni_metadata(path), metadata)
            self.assertEqual(infer_pixel_nm(path), 117.0)
            self.assertEqual(infer_na(path), 1.45)

    def test_infer_z_spacing_and_metadata_time_interval(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            paths: list[Path] = []
            for time in range(3):
                for z in range(2):
                    path = folder / f"sample_30s-FOV_posXY0_channels_t{time}_posZ{z}.tif"
                    metadata = {
                        "timestamp_us": 1_000_000 + time * 30_000_000 + z * 1_000_000,
                        "StagePos_um": [0.0, 0.0, -5.0 - z * 0.75],
                    }
                    tifffile.imwrite(path, np.zeros((2, 3), dtype=np.uint16), description=json.dumps(metadata))
                    paths.append(path)
            group = group_timelapse_files(parse_timelapse_filename(path) for path in paths)[0]

            self.assertEqual(infer_z_spacing_um(group), 0.75)
            self.assertEqual(infer_time_interval_s(group), 30.0)

    def test_parse_time_interval_from_path_skips_ms_tokens(self) -> None:
        path = Path("/data/sample_JFX650-200ms15f_30s-FOV/pos_0/file.tif")

        self.assertEqual(parse_time_interval_from_path(path), 30.0)

    def test_write_timelapse_tiff_round_trips_tzyx_and_tyx_shapes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            stack = np.random.default_rng(1).random((2, 3, 4, 5)).astype(np.float32)
            mip = stack.max(axis=1)
            stack_path = Path(tmp) / "stack.tif"
            mip_path = Path(tmp) / "mip.tif"

            write_timelapse_tiff(
                stack_path,
                stack,
                axes="TZYX",
                pixel_size_um=0.0585,
                z_spacing_um=0.75,
                time_interval_s=30.0,
            )
            write_timelapse_tiff(
                mip_path,
                mip,
                axes="TYX",
                pixel_size_um=0.0585,
                time_interval_s=30.0,
            )

            with tifffile.TiffFile(stack_path) as tif:
                self.assertEqual(tif.series[0].shape, stack.shape)
                self.assertEqual(tif.series[0].axes, "TZYX")
                self.assertEqual(tif.imagej_metadata["spacing"], 0.75)
                self.assertEqual(tif.imagej_metadata["finterval"], 30.0)
            with tifffile.TiffFile(mip_path) as tif:
                self.assertEqual(tif.series[0].shape, mip.shape)
                self.assertEqual(tif.series[0].axes, "TYX")


if __name__ == "__main__":
    unittest.main()
