"""Filename/resume fixtures only: no SACD reconstruction is executed."""
import copy
import json

import numpy as np
import pytest

from sacdpy import continuous_batch as continuous
from sacdpy import multicolor_zstack as multicolor
from sacdpy.discrete_timelapse import write_timelapse_tiff
from sacdpy.output_naming import resolve_recorded_paths
from tests import test_continuous_batch as continuous_fixtures
from tests import test_multicolor_zstack as multicolor_fixtures


def fixture_dataset(tmp_path, kind):
    raw = tmp_path / "20260917_ONI-dataset"
    output = tmp_path / "out"
    if kind == "timelapse":
        helper = continuous_fixtures.ContinuousBatchTests()
        helper._write_grid(raw / "FOV", "ONI-dataset_FOV")
        config = helper._config(raw, output)
        module = continuous
    else:
        helper = multicolor_fixtures.MulticolorZStackTests()
        helper._write_fov(raw)
        config = helper._config(raw, output)
        module = multicolor
    config["processing"]["max_workers"] = 1
    legacy_config = copy.deepcopy(config)
    legacy_config["processing"]["output_prefix_source"] = "movie_prefix"
    old = module.build_batch_plan(legacy_config).fovs[0]
    if kind == "timelapse":
        old.result_dir.mkdir(parents=True)
        stack = np.ones((2, 2, 8, 10), np.float32)
        for path, data, axes in ((old.stack_output, stack, "TZYX"), (old.mip_output, stack.max(axis=1), "TYX")):
            write_timelapse_tiff(path, data, axes=axes, pixel_size_um=0.0585,
                                z_spacing_um=old.z_spacing_um, time_interval_s=old.time_interval_s,
                                intensity_transform="order_root")
        paths = (old.stack_output, old.mip_output)
        root = old.result_dir
    else:
        old.output_root.mkdir(parents=True)
        stack = np.ones((2, 6, 8), np.float32)
        for output_paths in old.outputs:
            multicolor.write_imagej_tiff(output_paths.stack, stack, axes="ZYX", pixel_size_um=0.0585,
                                        z_spacing_um=old.z_spacing_um, intensity_transform="order_root")
            multicolor.write_imagej_tiff(output_paths.mip, stack.max(axis=0), axes="YX", pixel_size_um=0.0585,
                                        intensity_transform="order_root")
            multicolor.write_imagej_tiff(output_paths.raw_max_mip, np.ones((3, 4), np.uint16), axes="YX",
                                        pixel_size_um=0.117, reference_projection="frame_max_then_z_max")
        paths = old.all_outputs
        root = old.output_root
    record = module._result_record(old, "written", 0.0, intensity_transform="order_root")
    manifest = root / "_processing" / "manifest.json"
    manifest.parent.mkdir()
    manifest.write_text(json.dumps({"fovs": [record]}))
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    return module, config, config_path, old, paths, manifest


@pytest.mark.parametrize("kind", ["timelapse", "multicolor"])
def test_legacy_manifest_paths_resume_without_reconstruction_or_renaming(tmp_path, monkeypatch, kind):
    module, config, config_path, old, paths, manifest = fixture_dataset(tmp_path, kind)
    monkeypatch.setattr(module, "process_fov", lambda *a, **kw: pytest.fail("must resume, not reconstruct"))
    before = {p: p.read_bytes() for p in (*paths, manifest)}
    plan = module.resolve_resume_plan(module.build_batch_plan(config), config)
    selected = plan.fovs[0]
    selected_paths = (selected.stack_output, selected.mip_output) if kind == "timelapse" else selected.all_outputs
    assert selected_paths == paths
    results = module.run_batch(config, config_path, progress_callback=lambda _: None)
    assert results[0]["status"] == "resumed_manifest"
    assert before == {p: p.read_bytes() for p in before}


@pytest.mark.parametrize("kind", ["timelapse", "multicolor"])
@pytest.mark.parametrize("problem", ["orphan", "both_names", "identity", "intensity", "missing_member", "missing_all"])
def test_unsafe_legacy_resume_stops_before_provenance_writes(tmp_path, monkeypatch, kind, problem):
    module, config, config_path, old, paths, manifest = fixture_dataset(tmp_path, kind)
    canonical = module.build_batch_plan(config).fovs[0]
    if problem == "orphan":
        manifest.unlink()
    elif problem == "both_names":
        path = canonical.stack_output if kind == "timelapse" else canonical.all_outputs[0]
        path.write_bytes(paths[0].read_bytes())
    elif problem in {"identity", "intensity"}:
        data = json.loads(manifest.read_text())
        data["fovs"][0]["source_folder" if problem == "identity" else "intensity_transform"] = "wrong"
        manifest.write_text(json.dumps(data))
    elif problem == "missing_member":
        paths[0].unlink()
    else:
        for path in paths:
            path.unlink()
    monkeypatch.setattr(module, "_prepare_provenance", lambda *a, **kw: pytest.fail("must stop before writing"))
    with pytest.raises(ValueError):
        module.run_batch(config, config_path, progress_callback=lambda _: None)


@pytest.mark.parametrize("kind", ["timelapse", "multicolor"])
def test_failed_entry_without_files_is_retried_under_short_names(tmp_path, monkeypatch, kind):
    module, config, config_path, old, paths, manifest = fixture_dataset(tmp_path, kind)
    for path in paths:
        path.unlink()
    data = json.loads(manifest.read_text())
    data["fovs"][0]["status"] = "failed"
    manifest.write_text(json.dumps(data))
    calls = []
    def fake_process(fov, settings, **kwargs):
        calls.append(fov)
        return module._result_record(fov, "written", 0.0, intensity_transform="order_root")
    monkeypatch.setattr(module, "process_fov", fake_process)
    module.run_batch(config, config_path, progress_callback=lambda _: None)
    assert len(calls) == 1
    assert (calls[0].group.prefix if kind == "timelapse" else calls[0].prefix) == "FOV"


@pytest.mark.parametrize("movie_prefix", ["long_dataset_SHA-FOV-10", "SHA-FOV-10"])
def test_timelapse_short_default_and_alias_override(tmp_path, movie_prefix):
    helper = continuous_fixtures.ContinuousBatchTests()
    raw = tmp_path / "raw"
    helper._write_grid(raw / "SHA-FOV-10", movie_prefix)
    config = helper._config(raw, tmp_path / "out")
    fov = continuous.build_batch_plan(config).fovs[0]
    assert fov.stack_output.name == "SHA-FOV-10-SACDpy-647-posXY0-TZYX.tif"
    assert fov.mip_output.name == "SHA-FOV-10-SACDpy-647-posXY0-MIP-TYX.tif"
    config["datasets"][0]["output_prefix_aliases"] = {"SHA-FOV-10": "explicit"}
    assert continuous.build_batch_plan(config).fovs[0].stack_output.name.startswith("explicit-")


def test_orphan_dataset_prefixed_alias_is_detected_without_suffix_false_positive(tmp_path):
    canonical = (tmp_path / "FOV.tif",)
    (tmp_path / "another-FOV.tif").touch()
    kwargs = dict(canonical=canonical, legacy=(), record=None, recorded_paths=(),
                  expected={"relative_fov": "FOV"}, dataset_name="20260917_ONI-dataset")
    assert resolve_recorded_paths(**kwargs) == canonical
    (tmp_path / "ONI-dataset_FOV.tif").touch()
    with pytest.raises(ValueError, match="lack successful manifest"):
        resolve_recorded_paths(**kwargs)
