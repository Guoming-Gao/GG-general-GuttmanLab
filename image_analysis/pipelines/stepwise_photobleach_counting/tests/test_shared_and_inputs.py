from __future__ import annotations

import json
from pathlib import Path

import nbformat
import numpy as np
import yaml

from conftest import write_stack
from pbsa_shared import acquisition_profile, load_config, output_dirs
from step01_inspect_inputs import build_manifest


def test_acquisition_profile_is_stable():
    profile = acquisition_profile({"LaserWavelength_nm":[405,640],"LaserActive":[False,True],"LaserPowerPercent":[0,15],"Exposure_ms":100,"cameraBinning":1})
    assert profile == "640nm_15pct_100ms_1x"


def test_single_config_resolves_dataset_outputs(tmp_path: Path):
    raw=tmp_path/"raw"; raw.mkdir(); output=tmp_path/"out"; output.mkdir()
    path=tmp_path/"config.yaml"; path.write_text(yaml.safe_dump({"output_root":str(output),"datasets":[{"input_dir":str(raw),"expected_files":0}],"spotiflow":{"probability_threshold":0.4}}))
    cfg=load_config(path)
    assert cfg["datasets"][0]["name"]=="raw"
    assert Path(cfg["datasets"][0]["output_dir"])==output/"raw"
    assert output_dirs(cfg["datasets"][0]["output_dir"])["drift_correction"].name=="02_drift_correction"


def test_manifest_reads_oni_metadata(tmp_path: Path):
    raw=tmp_path/"dataset"; raw.mkdir(); write_stack(raw/"sample-FOV.tif",np.zeros((6,20,24),np.uint16),laser_power=20)
    cfg={"file_glob":"*.tif","expected_shape_yx":[20,24],"expected_dtype":"uint16","expected_pixel_size_um":.117,"pixel_size_tolerance_um":.002}
    dataset={"name":"dataset","input_dir":str(raw),"expected_files":1}
    manifest=build_manifest(cfg,dataset)
    assert manifest.loc[0,"status"]=="accepted"
    assert manifest.loc[0,"acquisition_profile"]=="640nm_20pct_100ms_1x"


def test_notebook_has_one_parameter_cell():
    notebook=nbformat.read(Path(__file__).parents[1]/"PBSA_HILO_pipeline.ipynb",as_version=4); nbformat.validate(notebook)
    assert sum("parameters" in cell.metadata.get("tags",[]) for cell in notebook.cells)==1

