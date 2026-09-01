from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import tifffile

from sacdpy.multicolor_single_z import (
    _channel_specs,
    _partial_path,
    build_batch_plan,
    preflight_summary,
    process_fov,
    read_movie_info,
    run_batch,
    validate_fov_outputs,
)


class MulticolorSingleZTests(unittest.TestCase):
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
                # Acquisition-step order, not ascending wavelength order.
                "channels": [
                    {"name": "AF647", "wavelength_nm": 647.0, "metadata_laser_nm": 640.0, "camera_half": "right"},
                    {"name": "AF488", "wavelength_nm": 488.0, "metadata_laser_nm": 488.0, "camera_half": "left"},
                    {"name": "Hoechst", "wavelength_nm": 405.0, "metadata_laser_nm": 405.0, "camera_half": "left"},
                ],
            },
            "selected_fov_folders": [],
            "exclude_fovs": {},
        }

    def _metadata(self, repeats: tuple[int, ...] = (2, 2, 2)) -> dict:
        active_indices = (3, 1, 0)
        steps = []
        for count, active_index in zip(repeats, active_indices, strict=True):
            values = [0.0, 0.0, 0.0, 0.0]
            values[active_index] = 10.0
            steps.append(
                {"nRepeats": count, "states": [{"group": 0, "record": True, "values": values}]}
            )
        return {
            "LaserWavelength_nm": [405, 488, 561, 640],
            "laserProgram": {"nLights": 4, "steps": steps},
            "PixelSize_um": 0.117,
            "Objective_NA": 1.45,
            "ROI": [0, 0, 8, 3],
        }

    def _write_movie(
        self,
        path: Path,
        *,
        width: int = 8,
        repeats: tuple[int, ...] = (2, 2, 2),
        metadata: dict | None = None,
    ) -> np.ndarray:
        raw = np.empty((sum(repeats), 3, width), dtype=np.uint16)
        halfwidth = width // 2
        start = 0
        for channel_index, count in enumerate(repeats):
            for frame_index in range(count):
                raw[start + frame_index, :, :halfwidth] = 10 + 10 * channel_index + frame_index
                raw[start + frame_index, :, halfwidth:] = 100 + 10 * channel_index + frame_index
            start += count
        path.parent.mkdir(parents=True, exist_ok=True)
        md = metadata if metadata is not None else self._metadata(repeats)
        md["ROI"] = [0, 0, width, 3]
        tifffile.imwrite(path, raw, photometric="minisblack", description=json.dumps(md))
        return raw

    def _write_fov(self, raw_root: Path, folder: str = "FOV-10", file_prefix: str = "wrong-FOV-1") -> Path:
        path = raw_root / folder / "pos_0" / f"{file_prefix}_posXY0_channels_t0_posZ0.tif"
        self._write_movie(path)
        return path

    def test_plan_uses_folder_prefix_and_metadata_driven_step_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            raw_root = Path(tmp) / "raw"
            output_root = raw_root / "SACDpy_results"
            self._write_fov(raw_root)
            plan = build_batch_plan(self._config(raw_root, output_root))
            fov = plan.fovs[0]

            self.assertEqual(fov.prefix, "FOV-10")
            self.assertEqual(fov.movie.active_lasers_nm, (640.0, 488.0, 405.0))
            self.assertEqual(fov.movie.frame_ranges, ((0, 2), (2, 4), (4, 6)))
            self.assertEqual(preflight_summary(plan)["reconstructions"], 3)
            self.assertEqual(preflight_summary(plan)["outputs"], 6)
            self.assertEqual(fov.output_for("AF647").sacd.name, "FOV-10__AF647-SACD.tif")
            self.assertEqual(fov.output_for("AF647").mip.name, "FOV-10__AF647-MIP.tif")

    def test_process_routes_channels_and_writes_only_two_calibrated_yx_outputs(self) -> None:
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

            with patch("sacdpy.multicolor_single_z.reconstruct", side_effect=fake_reconstruct):
                result = process_fov(fov, config["processing"])

            self.assertEqual(result["status"], "written")
            self.assertEqual([call[1] for call in calls], [647.0, 488.0, 405.0])
            self.assertEqual([call[2] for call in calls], ["order_root"] * 3)
            self.assertEqual([call[0] for call in calls], [100.5, 20.5, 30.5])
            self.assertEqual(len(list(output_root.glob("*.tif"))), 6)
            shapes = validate_fov_outputs(fov)["output_shapes"]["AF647"]
            self.assertEqual(shapes, {"sacd_shape": [6, 8], "mip_shape": [3, 4]})
            np.testing.assert_array_equal(tifffile.imread(fov.output_for("AF647").mip), 101)
            with tifffile.TiffFile(fov.output_for("AF647").mip) as tif:
                self.assertEqual(tif.imagej_metadata["reference_projection"], "frame_max")
            for output in fov.outputs:
                self.assertFalse(_partial_path(output.sacd).exists())
                self.assertNotIn("ZYX", output.sacd.name)
            self.assertEqual(process_fov(fov, config["processing"])["status"], "skipped_existing")

    def test_rejects_multiple_z_planes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            raw_root = Path(tmp) / "raw"
            self._write_fov(raw_root, file_prefix="sample")
            second = raw_root / "FOV-10" / "pos_0" / "sample_posXY0_channels_t0_posZ1.tif"
            self._write_movie(second)
            with self.assertRaisesRegex(ValueError, "exactly one ONI TIFF/z-plane"):
                build_batch_plan(self._config(raw_root, raw_root / "out"))

    def test_rejects_wrong_active_laser_and_ambiguous_active_lasers(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "sample_posXY0_channels_t0_posZ0.tif"
            config = self._config(root, root / "out")
            channels = _channel_specs(config)
            metadata = self._metadata()
            metadata["laserProgram"]["steps"][0]["states"][0]["values"] = [0, 0, 10, 0]
            self._write_movie(path, metadata=metadata)
            with self.assertRaisesRegex(ValueError, "expected active metadata laser 640"):
                read_movie_info(path, channels)

            metadata = self._metadata()
            metadata["laserProgram"]["steps"][0]["states"][0]["values"] = [0, 10, 0, 10]
            self._write_movie(path, metadata=metadata)
            with self.assertRaisesRegex(ValueError, "exactly one active laser"):
                read_movie_info(path, channels)

    def test_rejects_duplicate_channel_names_and_partial_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            raw_root = Path(tmp) / "raw"
            self._write_fov(raw_root)
            config = self._config(raw_root, raw_root / "out")
            config["processing"]["channels"][1]["name"] = "AF647"
            with self.assertRaisesRegex(ValueError, "Duplicate channel name"):
                build_batch_plan(config)

            config = self._config(raw_root, raw_root / "out")
            fov = build_batch_plan(config).fovs[0]
            fov.output_root.mkdir(parents=True)
            fov.output_for("AF647").sacd.touch()
            with self.assertRaisesRegex(FileExistsError, "Partial final"):
                process_fov(fov, config["processing"])

    def test_manifest_resume_rejects_incompatible_intensity_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw_root = root / "raw"
            output_root = raw_root / "SACDpy_results"
            self._write_fov(raw_root)
            config = self._config(raw_root, output_root)
            with patch("sacdpy.multicolor_single_z.reconstruct", return_value=np.ones((6, 8), np.float32)):
                first = run_batch(config, progress_callback=lambda _: None)
            self.assertEqual(first[0]["status"], "written")
            manifest_path = output_root / "_processing" / "manifest.json"
            manifest = json.loads(manifest_path.read_text())
            manifest["fovs"][0]["intensity_transform"] = "raw_cumulant"
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "cannot be resumed"):
                run_batch(config, progress_callback=lambda _: None)

    def test_incompatible_manifest_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            raw_root = Path(tmp) / "raw"
            output_root = raw_root / "SACDpy_results"
            self._write_fov(raw_root)
            processing_dir = output_root / "_processing"
            processing_dir.mkdir(parents=True)
            (processing_dir / "manifest.json").write_text(json.dumps({"pipeline": "multicolor_zstack"}))
            with self.assertRaisesRegex(ValueError, "incompatible manifest"):
                run_batch(self._config(raw_root, output_root), progress_callback=lambda _: None)


if __name__ == "__main__":
    unittest.main()
