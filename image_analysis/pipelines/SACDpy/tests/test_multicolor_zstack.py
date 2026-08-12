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
    legacy_ome_paths,
    load_config,
    migrate_legacy_ome_dataset,
    preflight_summary,
    process_fov,
    read_movie_info,
    run_batch,
    select_channel_frames,
    validate_fov_outputs,
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
                "intensity_transform": "order_root",
                "subfactor": 0.8,
                "channels": [
                    {"name": "405", "wavelength_nm": 405.0, "metadata_wavelengths_nm": [405.0], "camera_half": "left"},
                    {"name": "488", "wavelength_nm": 488.0, "metadata_wavelengths_nm": [488.0], "camera_half": "left"},
                    {"name": "561", "wavelength_nm": 561.0, "metadata_wavelengths_nm": [561.0], "camera_half": "left"},
                    {"name": "647", "wavelength_nm": 647.0, "metadata_wavelengths_nm": [640.0, 647.0], "camera_half": "right"},
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

    def _write_movie(self, path: Path, z_index: int, *, width: int = 8, repeats=(2, 2, 2, 2)) -> np.ndarray:
        raw = np.empty((sum(repeats), 3, width), dtype=np.uint16)
        halfwidth = width // 2
        start = 0
        for channel_index, count in enumerate(repeats):
            # Frame 2 is brighter so the raw-reference frame-max behavior is testable.
            for frame_index in range(count):
                raw[start + frame_index, :, :halfwidth] = 10 + channel_index + frame_index + z_index
                raw[start + frame_index, :, halfwidth:] = 20 + channel_index + frame_index + z_index
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

    def test_plan_discovers_z_stack_and_builds_twelve_named_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            raw_root = Path(tmp) / "raw"
            output_root = raw_root / "SACDpy_results"
            self._write_fov(raw_root)
            plan = build_batch_plan(self._config(raw_root, output_root))
            fov = plan.fovs[0]
            self.assertEqual(preflight_summary(plan)["reconstructions"], 8)
            self.assertEqual(preflight_summary(plan)["outputs"], 12)
            self.assertEqual(fov.z_indices, (0, 1))
            self.assertEqual(fov.movies[0].frame_ranges, ((0, 2), (2, 4), (4, 6), (6, 8)))
            self.assertEqual(fov.z_spacing_um, 0.667)
            self.assertEqual(fov.output_for("405").stack.name, "sample_FOV__405-SACD-ZYX.tif")
            self.assertEqual(fov.output_for("647").raw_max_mip.name, "sample_FOV__647-raw-max-MIP-YX.tif")

    def test_plan_rejects_missing_z_plane(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            raw_root = Path(tmp) / "raw"
            self._write_fov(raw_root, z_indices=(0, 2))
            with self.assertRaisesRegex(ValueError, "Missing z-plane"):
                build_batch_plan(self._config(raw_root, raw_root / "out"))

    def test_metadata_accepts_640_for_647_but_rejects_odd_width(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            channels = _channel_specs(self._config(root, root / "out"))
            odd_path = root / "odd_posXY0_channels_t0_posZ0.tif"
            self._write_movie(odd_path, 0, width=7)
            with self.assertRaisesRegex(ValueError, "width must be even"):
                read_movie_info(odd_path, channels)

    def test_select_channel_frames_matches_camera_split_midpoint(self) -> None:
        raw = np.zeros((8, 3, 8), dtype=np.uint16)
        raw[0:2, :, :4] = 11
        raw[6:8, :, 4:] = 23
        np.testing.assert_array_equal(select_channel_frames(raw, (0, 2), "left"), 11)
        np.testing.assert_array_equal(select_channel_frames(raw, (6, 8), "right"), 23)

    def test_process_routes_channels_and_writes_calibrated_imagej_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            raw_root = Path(tmp) / "raw"
            output_root = raw_root / "SACDpy_results"
            self._write_fov(raw_root)
            config = self._config(raw_root, output_root)
            fov = build_batch_plan(config).fovs[0]
            calls: list[tuple[float, float, str]] = []

            def fake_reconstruct(stack: np.ndarray, params) -> np.ndarray:
                calls.append((float(stack.mean()), params.wavelength_nm, params.intensity_transform))
                return np.full((6, 8), stack.mean(), dtype=np.float32)

            with patch("sacdpy.multicolor_zstack.reconstruct", side_effect=fake_reconstruct):
                result = process_fov(fov, config["processing"])

            self.assertEqual(result["status"], "written")
            self.assertEqual([value[1] for value in calls[:4]], [405.0, 488.0, 561.0, 647.0])
            self.assertEqual([value[2] for value in calls], ["order_root"] * 8)
            shapes = result["output_shapes"]["405"]
            self.assertEqual(shapes["stack_shape"], [2, 6, 8])
            self.assertEqual(shapes["raw_max_mip_shape"], [3, 4])
            self.assertEqual(len(list(output_root.glob("*.tif"))), 12)
            validate_fov_outputs(fov)
            raw_405 = tifffile.imread(fov.output_for("405").raw_max_mip)
            raw_647 = tifffile.imread(fov.output_for("647").raw_max_mip)
            self.assertTrue(np.all(raw_405 == 12))
            self.assertTrue(np.all(raw_647 == 25))
            for output in fov.outputs:
                self.assertFalse(_partial_path(output.stack).exists())
            self.assertEqual(process_fov(fov, config["processing"])["status"], "skipped_existing")

    def test_partial_final_set_is_not_silently_overwritten(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            raw_root = Path(tmp) / "raw"
            self._write_fov(raw_root)
            config = self._config(raw_root, raw_root / "out")
            fov = build_batch_plan(config).fovs[0]
            fov.output_root.mkdir(parents=True)
            fov.output_for("405").stack.touch()
            with self.assertRaisesRegex(FileExistsError, "Partial final multicolor output set"):
                process_fov(fov, config["processing"])

    def test_manifest_resume_rejects_missing_intensity_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw_root = root / "raw"
            output_root = raw_root / "SACDpy_results"
            self._write_fov(raw_root)
            config = self._config(raw_root, output_root)
            config_path = root / "config.json"
            config_path.write_text(json.dumps(config))
            with patch("sacdpy.multicolor_zstack.reconstruct", return_value=np.ones((6, 8), np.float32)):
                first = run_batch(config, config_path, progress_callback=lambda _: None)
            self.assertEqual(first[0]["status"], "written")
            manifest_path = output_root / "_processing" / "manifest.json"
            manifest = json.loads(manifest_path.read_text())
            manifest["fovs"][0].pop("intensity_transform")
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "cannot be resumed"):
                run_batch(config, config_path, progress_callback=lambda _: None)

    def test_explicit_legacy_migration_roots_values_and_removes_only_ome_pair(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw_root = root / "raw"
            output_root = raw_root / "SACDpy_results"
            self._write_fov(raw_root)
            config = self._config(raw_root, output_root)
            config_path = root / "config.json"
            config_path.write_text(json.dumps(config))
            fov = build_batch_plan(config).fovs[0]
            output_root.mkdir(parents=True)
            legacy_stack, legacy_mip = legacy_ome_paths(fov)
            values = np.arange(1, 4 * 2 * 6 * 8 + 1, dtype=np.float32).reshape(4, 2, 6, 8)
            raw_cumulants = values**2
            metadata = {"axes": "CZYX", "Channel": {"Name": ["405", "488", "561", "647"]}}
            tifffile.imwrite(legacy_stack, raw_cumulants, ome=True, metadata=metadata)
            tifffile.imwrite(
                legacy_mip,
                raw_cumulants.max(axis=1),
                ome=True,
                metadata={"axes": "CYX", "Channel": {"Name": ["405", "488", "561", "647"]}},
            )

            results = migrate_legacy_ome_dataset(config, config_path, progress_callback=lambda _: None)

            self.assertEqual(results[0]["status"], "migrated_order_root")
            self.assertFalse(legacy_stack.exists())
            self.assertFalse(legacy_mip.exists())
            self.assertEqual(len(list(output_root.glob("*.tif"))), 12)
            np.testing.assert_allclose(tifffile.imread(fov.output_for("405").stack), values[0])
            manifest = json.loads((output_root / "_processing" / "manifest.json").read_text())
            self.assertEqual(manifest["intensity_transform"], "order_root")
            self.assertEqual(manifest["migration"]["legacy_ome_files_removed"], 2)


if __name__ == "__main__":
    unittest.main()
