from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from conftest import write_oni_tiff
from run_ONI_SPT_batch import BatchProgress, preflight_batch
from spt_shared import load_config, output_dirs, write_run_metadata


def _configs(tmp_path: Path) -> tuple[Path, Path]:
    raw = tmp_path / "raw"; raw.mkdir()
    output = tmp_path / "processed"; output.mkdir()
    common = tmp_path / "common.yaml"
    common.write_text(yaml.safe_dump({
        "spotiflow": {"pretrained_model": "general", "probability_threshold": 0.4},
        "bandpass": {"sigma_low": 1, "sigma_high": 3},
    }))
    dataset = tmp_path / "dataset.yaml"
    dataset.write_text(yaml.safe_dump({
        "extends": "common.yaml", "dataset": "test-production", "input_dir": str(raw),
        "output_dir": str(output), "layout": {"mode": "single", "expected_width": 8},
        "channels": {"full": {"name": "spt", "role": "spt"}},
        "allowed_frame_intervals_s": [0.03],
    }, sort_keys=False))
    return common, dataset


def test_inherited_config_is_snapshotted_and_resume_is_fingerprint_safe(tmp_path):
    common, dataset = _configs(tmp_path)
    cfg = load_config(dataset)
    assert cfg["spotiflow"]["probability_threshold"] == 0.4
    write_run_metadata(cfg, command="run production")
    metadata = output_dirs(cfg["output_dir"])["metadata"]
    assert (metadata / "pipeline_config.yaml").read_text() == dataset.read_text()
    assert (metadata / "config_sources" / "00_common.yaml").read_text() == common.read_text()
    assert (metadata / "config_sources" / "01_dataset.yaml").read_text() == dataset.read_text()
    assert (metadata / "run_command.txt").read_text() == "run production\n"
    changed = copy.deepcopy(cfg); changed["spotiflow"]["probability_threshold"] = 0.5
    with pytest.raises(RuntimeError, match="changed processing parameters"):
        write_run_metadata(changed, resume=True)


def test_batch_preflight_uses_top_level_files_and_updates_eta_state(tmp_path):
    _, dataset = _configs(tmp_path)
    cfg = load_config(dataset); raw = Path(cfg["input_dir"])
    write_oni_tiff(raw / "condition-FOV-1.tif", np.zeros((3, 6, 8), np.uint16))
    nested = raw / "SACD"; nested.mkdir()
    write_oni_tiff(nested / "ignored-FOV-1.tif", np.zeros((9, 6, 8), np.uint16))
    batch_path = tmp_path / "batch.yaml"
    batch_path.write_text(yaml.safe_dump({
        "name": "test-batch", "base_output_dir": str(tmp_path), "minimum_free_GiB": 0,
        "datasets": [{"config": str(dataset), "expected_accepted_files": 1,
                      "expected_accepted_frames": 3, "expected_skipped_files": 0}],
    }))
    batch, prepared, table = preflight_batch(batch_path)
    assert table.loc[0, "accepted_files"] == 1
    assert table.loc[0, "accepted_frames"] == 3
    status = tmp_path / "status.json"
    progress = BatchProgress(batch, prepared, status)
    progress.set_dataset("test-production", "running")
    row = prepared[0]["manifest"].iloc[0].to_dict()
    progress.callback("preprocess", "finish", 1, 1, row)
    saved = json.loads(status.read_text())
    assert saved["current_stage"] == "preprocess"
    assert saved["estimated_remaining_seconds"] >= 0
    assert (output_dirs(cfg["output_dir"])["metadata"] / "batch_progress.json").exists()


def test_production_notebook_has_one_parameter_cell_and_no_pilot_execution():
    import nbformat
    notebook = nbformat.read(Path(__file__).parents[1] / "ONI_SPT_pipeline.ipynb", as_version=4)
    nbformat.validate(notebook)
    parameter_cells = [cell for cell in notebook.cells if "parameters" in cell.metadata.get("tags", [])]
    assert len(parameter_cells) == 1
    text = "".join("".join(cell.get("source", [])) for cell in notebook.cells)
    assert "pilot_max_frames" not in text
    assert "--max-frames" not in text
