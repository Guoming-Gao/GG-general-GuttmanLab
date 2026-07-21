# ONI Single-Particle Tracking

This folder is the complete ONI SPT pipeline. The primary interface is
`ONI_SPT_pipeline.ipynb`; the ordered `stepXX_*.py` files are the readable,
independently runnable implementation.

The workflow replaces TrackMate detection/tracking with Spotiflow and LapTrack,
assigns trajectories to CellposeSAM masks, calculates the historical AIO
diffusion metrics, runs pooled saSPT by condition, and creates comparison reports.

## Quick start

Use the existing `smlm` environment. Its required dependencies are documented in
`environment.yml`; Spotiflow and LapTrack currently resolve to the local checkouts
under `/Users/gmgao/GGscripts/`.

Open `ONI_SPT_pipeline.ipynb`, edit only its parameter cell, and use **Run All**.
The notebook performs a production preflight and then runs or safely resumes the
configured batch; there is no pilot stage.

Terminal equivalent:

```bash
conda run -n smlm python run_ONI_SPT_batch.py preflight --config configs/production/batch.yaml
conda run -n smlm python run_ONI_SPT_batch.py run --config configs/production/batch.yaml --resume
```

Run one dataset independently, or explicitly replace one reviewed stage:

```bash
conda run -n smlm python run_ONI_SPT.py run --config configs/production/fvp.yaml --resume
conda run -n smlm python step04_detect_spots.py --config CONFIG.yaml --force
```

## Ordered stages

1. `step01_inspect_inputs.py` — validate ONI layout and metadata.
2. `step02_preprocess_videos.py` — split/copy channels and save signed float32 DoG.
3. `step03_segment_cells.py` — CellposeSAM segmentation from marker or SPT MIP.
4. `step04_detect_spots.py` — frame-resolved Spotiflow detections and MP4 QC.
5. `step05_link_trajectories.py` — LapTrack, cell assignment, edge tables, track QC.
6. `step06_calculate_diffusion.py` — historical AIO metrics and pooled saSPT.
7. `step07_generate_reports.py` — condition comparison PDF, PNGs, and summaries.
8. `step08_make_crop_qc_video.py` — regenerate a detailed cell/track crop video from existing results.

The default LapTrack settings are `<5 px` between adjacent frames, `<5 px` gap
closing over at most two frames, no splitting/merging, and at least five detected
locations per retained trajectory.

Conditions default to the filename prefix before the first case-insensitive
`FOV` token. Exceptional names can be mapped explicitly under
`analysis.condition_overrides` using either the TIFF filename or stem as the key.

## Production configuration

- `configs/production/common.yaml` contains shared processing parameters,
  including the explicit Spotiflow threshold.
- Each production dataset has a short YAML that inherits `common.yaml` and
  declares only its paths, timing, layout, and channel roles.
- `configs/production/batch.yaml` defines execution order and exact expected
  accepted/skipped file and frame counts.

Dataset configs are independently runnable. The two mixed SACD/SPT profiles use
top-level `*.tif` matching only, so acquisition subfolders are excluded. dRRM
accepts only 856-pixel inputs and uses left Hoechst/right SPT.

## Output folders

```text
00_run_metadata/                 manifests, source/resolved configs, provenance, ETA/status
01_preprocessed/raw_channels/    standardized raw channel BigTIFFs
01_preprocessed/bandpass_float32/ signed float32 DoG BigTIFFs and statistics
01_preprocessed/mips/            raw channel MIPs
02_segmentation/                 masks, regions, source image, overlay
03_spot_detection/               rich per-frame detections
04_trajectories/                 canonical/legacy tracks, summaries, edge tables
05_diffusion_analysis/           per-FOV AIO, condition AIO, pooled saSPT
06_reports/                      comparison PDF, 300-dpi plots, summary CSVs
07_quality_control/              detection/track/crop MP4s and trajectory overlays
```

The tracking-metrics dashboards in `07_quality_control/secondary_metrics/` are
secondary summaries. Use the frame-resolved MP4 and the 1200-dpi/vector overlay
to judge whether detections were linked correctly.

DoG is saved as signed `float32` in camera-intensity units. Negative values are
not clipped and frames are not integer-quantized. Legacy TrackMate intensity
columns use the positive DoG component without rounding; rich detection tables
also retain signed-DoG and raw-intensity measurements.

Every result folder contains the exact source configuration chain, standalone
resolved processing parameters, batch configuration, command, dependency/model
versions, and configuration fingerprint. Resume is refused if that fingerprint
changes unless `--force` is explicitly supplied.

## Crop QC video

Render existing results for one segmented cell:

```bash
conda run -n smlm python step08_make_crop_qc_video.py \
  --config CONFIG.yaml --fov FOV_NAME --channel spt --cell-id 3
```

Selectors may instead be `--track-id ID` or `--bbox X0 Y0 X1 Y1`. The video uses
at most 100 source frames and includes a 1 µm scale bar, source frame, physical
timestamp, `Δt`, a physical-time progress bar, detections, track IDs, and resolved
Spotiflow/LapTrack settings.

A manually cropped TIFF can be used only with its full-FOV coordinate origin:

```bash
conda run -n smlm python step08_make_crop_qc_video.py \
  --config CONFIG.yaml --fov FOV_NAME --channel spt --cell-id 3 \
  --cropped-tif CELL_CROP.tif --origin-x 120 --origin-y 80
```

## Historical tools

Original GUI scripts and plotting notebooks are preserved under `legacy/`.
The original root CLI filenames remain as non-GUI compatibility wrappers around
the new shared AIO, concatenation, and saSPT implementations.
