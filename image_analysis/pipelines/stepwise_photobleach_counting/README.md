# Stepwise Photobleach Counting

Compact HILO photobleaching pipeline following the neighboring ONI SPT workflow style. It performs cross-correlation drift correction, Spotiflow detection, cluster-aware local-background extraction, independently calibrated quickPBSA Bayesian counting, and stage-by-stage QC.

## Environment

Use the existing `smlm` Conda environment. Install the one missing analysis dependency once:

```bash
conda run -n smlm pip install quickpbsa==2021.0.1
```

Spotiflow is already installed from the local checkout and its `general` model is shared with `../BulkFluo_RDF/cache/spotiflow_models`.

## Workflow

The validated drift products were produced with the original stage runner. The corrected downstream reanalysis is driven by `run_cluster_reanalysis.py` and `config_reanalysis.yaml`:

Terminal equivalents:

```bash
conda run -n smlm python run_cluster_reanalysis.py build-units --config config_reanalysis.yaml --resume
conda run -n smlm python run_cluster_reanalysis.py extract-traces --config config_reanalysis.yaml --resume
conda run -n smlm python run_cluster_reanalysis.py calibrate-thresholds --config config_reanalysis.yaml
conda run -n smlm python run_cluster_reanalysis.py count --config config_reanalysis.yaml --resume
conda run -n smlm python run_cluster_reanalysis.py threshold-sensitivity --config config_reanalysis.yaml
conda run -n smlm python run_cluster_reanalysis.py validate-monomers --config config_reanalysis.yaml
conda run -n smlm python run_cluster_reanalysis.py report --config config_reanalysis.yaml
```

Calibration is independent by acquisition profile and analysis stream. A production threshold is `0.5 × median persistent late-step amplitude`, with bootstrap confidence intervals, baseline and ROI-area checks, a 30-trace visual-review montage, and a 0.3–0.7 quickPBSA sensitivity sweep. `run --resume` executes every revised stage and preserves completed per-FOV products.

## Output organization

The final result root is deliberately reader-oriented:

```text
00_README_RESULTS.md
01_report/
02_figures/{seed_level_clusters,grouped_clusters,extended_clusters,method_QC}/
03_tables/
04_experiments/<original experiment name>/{input_QC,drift_correction,spot_detection,puncta_traces,quickpbsa_counts}/
99_technical/{provenance,sensitivity_analysis,quickpbsa_native,exploratory_segmentation}/
```

Every processing stage creates machine-readable QC. Seed-level and grouped-cluster counts are co-primary but never pooled. Extended candidates remain a separately qualified exploratory stream. Counts above 40 are retained and marked `outside_validated_range`.

## Presentation figures

Generate the presentation-ready grouped-cluster histogram, ECDF, monomer-candidate traces, and highest-step traces with:

```bash
conda run -n smlm python presentation_plots.py --config config_reanalysis.yaml
```

The four PNG/SVG figure pairs are written directly to the result root. `03_tables/presentation_trace_selection.csv` records every trace chosen for the two 4 × 3 panels. Plot idealizations retain quickPBSA transition locations while refitting each displayed plateau height to the median of the exact corrected trace representation used for counting.

## Important implementation details

- Drift correction uses scikit-image's published [`phase_cross_correlation`](https://scikit-image.org/docs/stable/api/skimage.registration.html#skimage.registration.phase_cross_correlation) implementation (Guizar-Sicairos et al., 2008) with five-frame registration blocks and 0.1-pixel precision, then applies interpolated shifts to every original frame.
- Spotiflow uses the `general` model and a fixed probability threshold of `0.4`.
- Every Spotiflow 0.4 seed receives a fixed circular seed-level ROI even if local growth fails.
- Seeds with overlapping 2-pixel signal footprints form one deterministic grouped-cluster unit; member IDs are preserved.
- Seed signal uses a 2-pixel circle and a clipped 3.5–5-pixel ring that excludes every neighboring seed footprint.
- Marker-controlled watershed creates extended candidates only from successfully grown seeds. Failure affects morphology classification, not seed/group eligibility.
- Mean-difference and baseline-centered integrated representations are evaluated independently. The representation must pass late-baseline, ROI-area, synthetic, and visual checks before counting.
- A `monomer_candidate` is an accepted, valid-background, finite one-step seed or single-seed group trace. It is not proof of biochemical monomeric state.
- Native quickPBSA artifacts are retained alongside simpler canonical count tables.

## Tests

```bash
conda run -n smlm python -m pytest -q
```
