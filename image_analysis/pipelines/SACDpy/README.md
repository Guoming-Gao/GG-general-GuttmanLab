# SACDpy

Python implementation of the core SACD reconstruction pipeline from the MATLAB
SACDm package.

The first implementation covers the default 2D SACD workflow:

1. per-frame background offset removal
2. pre Richardson-Lucy deconvolution
3. lateral Fourier interpolation
4. autocumulant calculation
5. post Richardson-Lucy deconvolution
6. autocumulant-order root conversion to the public SACD intensity space

`reconstruct()` returns float32 order-root SACD intensity by default: an order
`n` raw cumulant is raised to `1 / n` after post-deconvolution. This follows the
original SACDm visualization convention without thresholding, percentile
stretching, max normalization, or integer conversion. The result is
intensity-like but remains in arbitrary units. Use
`SACDParams(intensity_transform="raw_cumulant")` (or CLI
`--intensity-transform raw-cumulant`) only for MATLAB core-parity diagnostics
and explicit backward-compatibility work; raw-cumulant magnitude is not a
calibrated camera intensity.

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

The multicolor z-stack, multicolor single-plane, and timelapse z-stack notebooks
default to `processing.max_workers = 2`. Set the visible `max_workers` notebook
setting to `1` for serial execution. FOVs remain sequential; independent
z/channel, channel, or time/z reconstructions within each FOV share a persistent
spawn-based process pool. Each worker uses one numerical-library thread and the
queue is bounded to twice the worker count. Memory demand grows with the number
of concurrent reconstructions; core count alone is not a suitable worker limit.

All five notebooks use Rich progress with counts, elapsed time, and ETA where
the stage has a known total. Single-image batch and Dox analysis remain serial.
The three parallel runners expose `progress_callback`; notebooks connect it to
`sacdpy.progress.PipelineProgress`. Dox reconstruction, repackaging, PNG refresh,
and Spotiflow entry points accept the same optional callback. Model loading and
other stages without measurable totals use indeterminate progress.

Reconstruction completion and output validation are reported separately.
Resumed FOVs advance immediately; failed work is not counted as successful.
Changing worker count does not invalidate compatible existing results. On
interruption, queued work is cancelled and running jobs are allowed to settle
before the pool exits; this can take up to the duration of a reconstruction.

Timelapse and multicolor z-stack outputs default to the descriptive FOV folder
name, without the parent dataset prefix. Explicit `output_prefix_aliases` still
override this name; channel, position, and axis suffixes are preserved. For
example, `long_dataset_SHA-FOV-10-SACDpy-647-posXY0-TZYX.tif` becomes
`SHA-FOV-10-SACDpy-647-posXY0-TZYX.tif` for newly processed FOVs.

Both runners and their notebook preflights use `resolve_resume_plan` to retain
completed, manifest-recorded legacy filenames. This read-only check verifies
FOV identity, intensity provenance, and TIFF outputs. It rejects conflicting
old/new variants, orphan legacy outputs, and incomplete legacy sets before
writing provenance or reconstructing anything. No automatic rename or migration
is performed. `build_batch_plan` alone describes the canonical new filenames;
wrap it with `resolve_resume_plan(plan, config)` when inspecting resumable paths.

- Validation report: `SACDpy_validation_report.md`
- Batch single-image notebook: `SACDpy_pipeline-batch_sinlgeSACD.ipynb`
- Canonical timelapse z-stack notebook: `SACDpy_pipeline-batch_timelapse_zstack_SACD.ipynb`
- Reproducible multi-dataset config: `configs/sacd_live_batch_2026.json`
- Resumable runner: `scripts/run_sacd_live_batch.py`
- Reproduce test-data validation:

```bash
PYTHONPATH=src conda run -n smlm python tests/validate_testdata.py
```
