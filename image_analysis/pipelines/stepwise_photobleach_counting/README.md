# Stepwise Photobleach Counting

Flat, stage-by-stage HILO photobleaching pipeline following the neighboring ONI SPT workflow style. It performs image-level cross-correlation drift correction, Spotiflow seed detection, local seed growth, quickPBSA-compatible local background extraction, and quickPBSA SIC/Bayesian step counting.

## Environment

Use the existing `smlm` Conda environment. Install the one missing analysis dependency once:

```bash
conda run -n smlm pip install quickpbsa==2021.0.1
```

Spotiflow is already installed from the local checkout and its `general` model is shared with `../BulkFluo_RDF/cache/spotiflow_models`.

## Workflow

The primary interface is `PBSA_HILO_pipeline.ipynb`. Edit its single parameter cell, then run one stage at a time and inspect the linked QC report before continuing.

Terminal equivalents:

```bash
conda run -n smlm python run_stepwise_pbsa.py inspect-inputs --config config.yaml
conda run -n smlm python run_stepwise_pbsa.py drift-correction --config config.yaml --resume
conda run -n smlm python run_stepwise_pbsa.py detect-and-grow-rois --config config.yaml --resume
conda run -n smlm python run_stepwise_pbsa.py extract-background-corrected-traces --config config.yaml --resume
conda run -n smlm python run_stepwise_pbsa.py quickpbsa-threshold-pilot --config config.yaml
conda run -n smlm python run_stepwise_pbsa.py count-photobleaching-steps --config config.yaml --resume
conda run -n smlm python run_stepwise_pbsa.py generate-reports --config config.yaml
```

Stage 05 pools pilot traces by acquisition profile across all configured datasets and automatically selects `0.5 × estimated single-step intensity`. It also retains the 0.3–0.7 sensitivity sweep. Explicit values under `quickpbsa.thresholds` in `config.yaml` override the automatic selections. `run --resume` executes the full workflow in order while preserving completed per-FOV products.

## Output organization

Each input dataset keeps its source folder name beneath the configured result root:

```text
00_run_metadata/
01_input_inspection/
02_drift_correction/
03_roi_detection_and_growth/
04_background_corrected_traces/
05_quickpbsa_threshold_pilot/
06_photobleaching_step_counts/
07_summary_reports/
```

Every processing stage creates machine-readable QC records and a stage report. The pooled threshold products are stored at the result root under `05_quickpbsa_threshold_pilot/`; the combined figures, tables, and final PDF are under `07_summary_reports/`. Point puncta and extended condensates are always reported separately. Counts above 40 are retained with `outside_validated_range` rather than silently mixed into the validated result range.

## Important implementation details

- Drift correction uses scikit-image's published [`phase_cross_correlation`](https://scikit-image.org/docs/stable/api/skimage.registration.html#skimage.registration.phase_cross_correlation) implementation (Guizar-Sicairos et al., 2008) with five-frame registration blocks and 0.1-pixel precision, then applies interpolated shifts to every original frame.
- Spotiflow uses the `general` model and a fixed probability threshold of `0.4`.
- Every ROI must contain a Spotiflow seed. One-seed regions no wider than 5 pixels are point puncta; larger or multi-seed regions are extended condensates.
- Point puncta use a 2-pixel signal circle and 3.5–5-pixel background ring. Condensates use their grown mask and a mask-shaped local background band.
- quickPBSA receives integrated local-background-corrected traces, avoiding artificial dilution of a single-emitter step in larger masks.
- Native quickPBSA artifacts are retained alongside simpler canonical count tables.

## Tests

```bash
conda run -n smlm python -m pytest -q
```
