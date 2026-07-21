# SACDpy

Python implementation of the core SACD reconstruction pipeline from the MATLAB
SACDm package.

The first implementation covers the default 2D SACD workflow:

1. per-frame background offset removal
2. pre Richardson-Lucy deconvolution
3. lateral Fourier interpolation
4. autocumulant calculation
5. post Richardson-Lucy deconvolution

Advanced MATLAB options such as registration, wavelet background subtraction,
GPU execution, and sparse Hessian deconvolution are intentionally not included
in this first core port.

## CLI

From a checkout, either install the package into the `smlm` environment:

```bash
conda run -n smlm python -m pip install -e .
```

or run commands with `PYTHONPATH=src`.

```bash
PYTHONPATH=src conda run -n smlm python -m sacdpy tests/testdata/input.tif output.tif \
  --pixel 117 --wavelength 560 --na 1.45
```

The CLI reads TIFF stacks as `TYX` and writes float32 TIFF outputs. Use
`--frames-per-sacd 25` to reconstruct non-overlapping 25-frame chunks; multiple
SACD frames are saved as a `TYX` TIFF stack.

The same option is available in Python:

```python
from sacdpy import SACDParams, reconstruct

params = SACDParams(pixel_nm=117, wavelength_nm=560, na=1.45, frames_per_sacd=25)
sacd = reconstruct(raw_stack, params)
```

## Repository layout

The package uses the standard Python `src/` layout. Importable SACDpy code
lives in `src/sacdpy`, while fixtures, validation outputs, and test helpers live
under `tests/`.

## Validation and Notebooks

- Validation report: `SACDpy_validation_report.md`
- Batch single-image notebook: `SACDpy_pipeline-batch_sinlgeSACD.ipynb`
- Canonical timelapse z-stack notebook: `SACDpy_pipeline-batch_timelapse_zstack_SACD.ipynb`
- Reproducible multi-dataset config: `configs/sacd_live_batch_2026.json`
- Resumable runner: `scripts/run_sacd_live_batch.py`
- Reproduce test-data validation:

```bash
PYTHONPATH=src conda run -n smlm python tests/validate_testdata.py
```
