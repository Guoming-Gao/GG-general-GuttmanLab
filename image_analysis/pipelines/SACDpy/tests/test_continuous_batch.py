from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import tifffile

from sacdpy.continuous_batch import (
    FOVPlan,
    _partial_path,
    build_batch_plan,
    load_config,
    preflight_summary,
    process_fov,
    validate_output_pair,
)
from sacdpy.discrete_timelapse import TimelapseGroup


class ContinuousBatchTests(unittest.TestCase):
    def _config(self, root: Path, output_root: Path) -> dict:
        return {
            "version": 1,
            "output_root": str(output_root),
            "processing": {
                "position_folder": "pos_0",
                "glob_pattern": "*.tif",
                "channel_name": "647",
                "incomplete_grid_policy": "error",
                "wavelength_nm": 647.0,
                "fallback_pixel_nm": 117.0,
                "fallback_na": 1.45,
                "mag": 2,
                "iter1": 1,
                "iter2": 1,
                "ac_order": 2,
                "subfactor": 0.8,
            },
            "datasets": [{"raw_root": str(root)}],
        }

    def _write_grid(self, fov: Path, prefix: str = "sample") -> None:
        pos = fov / "pos_0"
        pos.mkdir(parents=True)
        metadata = {"PixelSize_um": 0.117, "Objective_NA": 1.45}
        for time in range(2):
            for z in range(2):
                tifffile.imwrite(
                    pos / f"{prefix}_posXY0_channels_t{time}_posZ{z}.tif",
                    np.ones((3, 4, 5), dtype=np.uint16),
                    description=json.dumps(metadata),
                )

    def test_load_config_requires_supported_version(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.json"
            path.write_text(json.dumps({"version": 2, "datasets": [{}]}))
            with self.assertRaisesRegex(ValueError, "version"):
                load_config(path)

    def test_plan_applies_exclusion_selection_and_alias(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "raw"
            output = Path(tmp) / "out"
            self._write_grid(root / "keep", "raw-prefix")
            self._write_grid(root / "exclude", "excluded-prefix")
            config = self._config(root, output)
            config["datasets"][0].update(
                {
                    "selected_fov_folders": ["keep", "exclude"],
                    "exclude_fovs": {"exclude": "failed acquisition"},
                    "output_prefix_aliases": {"keep": "clean-prefix"},
                }
            )

            plan = build_batch_plan(config)

        self.assertEqual(preflight_summary(plan)["fovs"], 1)
        self.assertEqual(plan.movie_count, 4)
        self.assertEqual(plan.fovs[0].group.prefix, "clean-prefix")
        self.assertEqual(plan.exclusions[0]["reason"], "failed acquisition")
        self.assertEqual(plan.fovs[0].stack_output.parent, output / root.name)

    def test_plan_rejects_output_alias_collisions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "raw"
            output = Path(tmp) / "out"
            self._write_grid(root / "one", "one")
            self._write_grid(root / "two", "two")
            config = self._config(root, output)
            config["datasets"][0]["output_prefix_aliases"] = {"one": "same", "two": "same"}
            with self.assertRaisesRegex(ValueError, "collision"):
                build_batch_plan(config)

    def test_process_fov_uses_partial_files_and_resumes_pair(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw_path = root / "raw.tif"
            tifffile.imwrite(raw_path, np.ones((3, 4, 5), dtype=np.uint16))
            group = TimelapseGroup(
                "sample", "647", 0,
                {(time, z): raw_path for time in range(2) for z in range(2)},
            )
            stack_path = root / "sample-TZYX.tif"
            mip_path = root / "sample-MIP-TYX.tif"
            fov = FOVPlan(
                "dataset", root, root, root, ".", group, stack_path, mip_path,
                117.0, 1.45, 0.75, 30.0,
            )
            defaults = self._config(root, root)["processing"]
            with patch("sacdpy.continuous_batch.reconstruct", side_effect=lambda raw, params: raw[0]):
                written = process_fov(fov, defaults)
            resumed = process_fov(fov, defaults)

            self.assertEqual(written["status"], "written")
            self.assertEqual(resumed["status"], "skipped_existing")
            self.assertFalse(_partial_path(stack_path).exists())
            self.assertEqual(validate_output_pair(stack_path, mip_path)["stack_shape"], [2, 2, 4, 5])

    def test_eta_status_can_embed_a_result_without_a_cycle(self) -> None:
        result = {"status": "written"}
        status = {"last_result": dict(result)}
        result["batch_status"] = status
        self.assertIn('"batch_status"', json.dumps(result))


if __name__ == "__main__":
    unittest.main()
