from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import tifffile

from sacdpy.multicolor_zstack import (
    _channel_specs,
    _partial_path,
    build_batch_plan,
    load_config,
    preflight_summary,
    process_fov,
    read_movie_info,
    run_batch,
    select_channel_frames,
    validate_output_pair,
    write_ome_tiff,
)


class MulticolorZStackTests(unittest.TestCase):
    def _config(self, raw_root: Path, output_root: Path) -> dict:
        return {
            "version": 1,
            "raw_root": str(raw_root),
            "output_root": str(output_root),
            "processing": {
                "position_folder": "pos_0",
                "glob_pattern": "*.tif",
                "fallback_pixel_nm": 117.0,
                "fallback_na": 1.45,
                "mag": 2,
                "iter1": 1,
                "iter2": 1,
                "ac_order": 2,
                "subfactor": 0.8,
                "channels": [
                    {
                        "name": "405",
                        "wavelength_nm": 405.0,
                        "metadata_wavelengths_nm": [405.0],
                        "camera_half": "left",
                    },
                    {
                        "name": "488",
                        "wavelength_nm": 488.0,
                        "metadata_wavelengths_nm": [488.0],
                        "camera_half": "left",
                    },
                    {
                        "name": "561",
                        "wavelength_nm": 561.0,
                        "metadata_wavelengths_nm": [561.0],
                        "camera_half": "left",
                    },
                    {
                        "name": "647",
                        "wavelength_nm": 647.0,
                        "metadata_wavelengths_nm": [640.0, 647.0],
                        "camera_half": "right",
                    },
                ],
            },
        }

    def _metadata(self, z_index: int, repeats: tuple[int, ...] = (2, 2, 2, 2)) -> dict:
        return {
            "LaserWavelength_nm": [405, 488, 561, 640],
            "laserProgram": {"steps": [{"nRepeats": value} for value in repeats]},
            "PixelSize_um": 0.117,
            "Objective_NA": 1.45,
            "ROI": [0, 0, 8, 3],
            "StagePos_um": [0.0, 0.0, -0.667 * z_index],
        }

    def _write_movie(
        self,
        path: Path,
        z_index: int,
        *,
        width: int = 8,
        repeats: tuple[int, ...] = (2, 2, 2, 2),
    ) -> np.ndarray:
        raw = np.empty((sum(repeats), 3, width), dtype=np.uint16)
        halfwidth = width // 2
        start = 0
        for channel_index, count in enumerate(repeats):
            raw[start : start + count, :, :halfwidth] = 10 + channel_index
            raw[start : start + count, :, halfwidth:] = 20 + channel_index
            start += count
        path.parent.mkdir(parents=True, exist_ok=True)
        metadata = self._metadata(z_index, repeats)
        metadata["ROI"] = [0, 0, width, 3]
        tifffile.imwrite(path, raw, photometric="minisblack", description=json.dumps(metadata))
        return raw

    def _write_fov(self, raw_root: Path, name: str = "FOV", z_indices=(0, 1)) -> None:
        for z_index in z_indices:
            self._write_movie(
                raw_root / name / "pos_0" / f"sample_{name}_posXY0_channels_t0_posZ{z_index}.tif",
                z_index,
            )

    def test_load_config_validates_version_and_channel_camera(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.json"
            config = self._config(Path(tmp) / "raw", Path(tmp) / "out")
            path.write_text(json.dumps(config))
            self.assertEqual(load_config(path)["version"], 1)
            config["processing"]["channels"][0]["camera_half"] = "middle"
            path.write_text(json.dumps(config))
            with self.assertRaisesRegex(ValueError, "left or right"):
                load_config(path)

    def test_plan_discovers_complete_z_stack_and_metadata_frame_ranges(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            raw_root = Path(tmp) / "raw"
            output_root = raw_root / "SACDpy_results"
            self._write_fov(raw_root)

            plan = build_batch_plan(self._config(raw_root, output_root))

            self.assertEqual(preflight_summary(plan)["fovs"], 1)
            self.assertEqual(preflight_summary(plan)["reconstructions"], 8)
            self.assertEqual(plan.fovs[0].z_indices, (0, 1))
            self.assertEqual(
                plan.fovs[0].movies[0].frame_ranges,
                ((0, 2), (2, 4), (4, 6), (6, 8)),
            )
            self.assertEqual(plan.fovs[0].z_spacing_um, 0.667)
            self.assertEqual(plan.fovs[0].stack_output.parent, output_root)

    def test_plan_rejects_missing_z_plane(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            raw_root = Path(tmp) / "raw"
            self._write_fov(raw_root, z_indices=(0, 2))
            with self.assertRaisesRegex(ValueError, "Missing z-plane"):
                build_batch_plan(self._config(raw_root, raw_root / "out"))

    def test_metadata_accepts_640_for_647_but_rejects_odd_width(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            config = self._config(root, root / "out")
            channels = _channel_specs(config)
            odd_path = root / "odd_posXY0_channels_t0_posZ0.tif"
            self._write_movie(odd_path, 0, width=7)
            with self.assertRaisesRegex(ValueError, "width must be even"):
                read_movie_info(odd_path, channels)
            self.assertEqual(config["processing"]["channels"][-1]["wavelength_nm"], 647.0)

    def test_metadata_rejects_laser_repeat_total_that_differs_from_page_count(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "bad_posXY0_channels_t0_posZ0.tif"
            raw = np.zeros((7, 3, 8), dtype=np.uint16)
            tifffile.imwrite(
                path,
                raw,
                photometric="minisblack",
                description=json.dumps(self._metadata(0, (2, 2, 2, 2))),
            )
            with self.assertRaisesRegex(ValueError, "repeats total 8.*7 page"):
                read_movie_info(path, _channel_specs(self._config(root, root / "out")))

    def test_select_channel_frames_matches_camera_split_midpoint(self) -> None:
        raw = np.zeros((8, 3, 8), dtype=np.uint16)
        raw[0:2, :, :4] = 11
        raw[6:8, :, 4:] = 23
        left = select_channel_frames(raw, (0, 2), "left")
        right = select_channel_frames(raw, (6, 8), "right")
        self.assertEqual(left.shape, (2, 3, 4))
        self.assertTrue(np.all(left == 11))
        self.assertTrue(np.all(right == 23))

    def test_process_routes_frames_camera_halves_and_psf_wavelengths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            raw_root = Path(tmp) / "raw"
            output_root = raw_root / "SACDpy_results"
            self._write_fov(raw_root)
            config = self._config(raw_root, output_root)
            fov = build_batch_plan(config).fovs[0]
            calls: list[tuple[float, float, tuple[int, ...]]] = []
            progress_events: list[dict] = []

            def fake_reconstruct(stack: np.ndarray, params) -> np.ndarray:
                calls.append((float(stack.mean()), params.wavelength_nm, tuple(stack.shape)))
                return np.full(stack.shape[1:], stack.mean(), dtype=np.float32)

            with patch("sacdpy.multicolor_zstack.reconstruct", side_effect=fake_reconstruct):
                result = process_fov(
                    fov,
                    config["processing"],
                    progress_callback=progress_events.append,
                )

            self.assertEqual(result["status"], "written")
            expected_one_z = [
                (10.0, 405.0, (2, 3, 4)),
                (11.0, 488.0, (2, 3, 4)),
                (12.0, 561.0, (2, 3, 4)),
                (23.0, 647.0, (2, 3, 4)),
            ]
            self.assertEqual(calls, expected_one_z * 2)
            self.assertEqual(result["stack_shape"], [4, 2, 3, 4])
            self.assertEqual(result["mip_shape"], [4, 3, 4])
            self.assertEqual(len(progress_events), 8)
            self.assertEqual(progress_events[-1]["completed_in_fov"], 8)
            self.assertEqual(progress_events[-1]["channel"], "647")
            self.assertFalse(_partial_path(fov.stack_output).exists())

            resumed = process_fov(fov, config["processing"])
            self.assertEqual(resumed["status"], "skipped_existing")

    def test_ome_pair_has_channel_names_axes_physical_sizes_and_exact_mip(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            raw_root = Path(tmp) / "raw"
            self._write_fov(raw_root)
            fov = build_batch_plan(self._config(raw_root, raw_root / "out")).fovs[0]
            stack = np.arange(4 * 2 * 3 * 4, dtype=np.float32).reshape(4, 2, 3, 4)
            mip = stack.max(axis=1)
            stack_path = Path(tmp) / "stack.ome.tif"
            mip_path = Path(tmp) / "mip.ome.tif"
            write_ome_tiff(
                stack_path,
                stack,
                axes="CZYX",
                channels=fov.channels,
                pixel_size_um=0.0585,
                z_spacing_um=0.667,
            )
            write_ome_tiff(
                mip_path,
                mip,
                axes="CYX",
                channels=fov.channels,
                pixel_size_um=0.0585,
            )

            shapes = validate_output_pair(stack_path, mip_path, fov.channels)
            self.assertEqual(shapes["stack_shape"], [4, 2, 3, 4])
            with tifffile.TiffFile(stack_path) as tif:
                xml = tif.ome_metadata or ""
                self.assertIn('PhysicalSizeX="0.0585"', xml)
                self.assertIn('PhysicalSizeZ="0.667"', xml)
                self.assertIn('ExcitationWavelength="647.0"', xml)

    def test_partial_final_pair_is_not_silently_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            raw_root = Path(tmp) / "raw"
            self._write_fov(raw_root)
            config = self._config(raw_root, raw_root / "out")
            fov = build_batch_plan(config).fovs[0]
            fov.stack_output.parent.mkdir(parents=True)
            fov.stack_output.touch()
            with self.assertRaisesRegex(FileExistsError, "Partial final output pair"):
                process_fov(fov, config["processing"])

    def test_pilot_manifest_is_validated_and_resumed_by_full_run(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw_root = root / "raw"
            output_root = raw_root / "SACDpy_results"
            self._write_fov(raw_root)
            config = self._config(raw_root, output_root)
            config_path = root / "config.json"
            config_path.write_text(json.dumps(config))

            def fake_reconstruct(stack: np.ndarray, params) -> np.ndarray:
                return np.full(stack.shape[1:], params.wavelength_nm, dtype=np.float32)

            with patch("sacdpy.multicolor_zstack.reconstruct", side_effect=fake_reconstruct):
                pilot = run_batch(
                    config, config_path, max_new_fovs=1, progress_callback=lambda result: None
                )
                resumed = run_batch(config, config_path, progress_callback=lambda result: None)

            self.assertEqual(pilot[0]["status"], "written")
            self.assertEqual(resumed[0]["status"], "resumed_manifest")
            manifest = json.loads((output_root / "_processing" / "manifest.json").read_text())
            self.assertEqual(manifest["fovs"][0]["status"], "written")


if __name__ == "__main__":
    unittest.main()
