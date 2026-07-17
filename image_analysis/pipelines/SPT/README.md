# ONI SPT

`oni-spt` is a stand-alone, resumable pipeline for ONI single-particle tracking.
It replaces TrackMate detection/tracking with Spotiflow and LapTrack while
retaining compatibility with the existing AIO and saSPT analyses.

## Pipeline

1. Inspect ONI metadata and validate the configured acquisition layout.
2. Stream the input into standardized channel BigTIFFs.
3. Save DoG-filtered (`sigma=1,3`) TIFFs for SPT channels.
4. Segment a marker MIP, or an SPT MIP when no marker is present, with
   CellposeSAM.
5. Detect subpixel spots independently in every filtered SPT frame with
   Spotiflow.
6. Link detections with LapTrack and assign each trajectory to its modal cell
   mask.
7. Export rich CSVs, exact numeric TrackMate-compatible CSVs, summaries, and QC.

The source TIFFs are never modified. Existing outputs are not silently deleted:
use `--resume` to reuse a complete stage or `--force` to replace its individual
artifacts.

## Installation

Create the combined environment:

```bash
conda env create -f environment.yml
conda activate oni-spt
```

For development in the existing `smlm` environment:

```bash
conda run -n smlm pip install -e /Users/gmgao/GGscripts/spotiflow
conda run -n smlm pip install -e /Users/gmgao/GGscripts/laptrack
conda run -n smlm pip install -e ".[analysis,test]"
```

## Commands

Inspect inputs without running models:

```bash
oni-spt inspect --config configs/fvp-validation.yaml
oni-spt inspect --config configs/drrm-validation.yaml
```

Run a 50-frame validation:

```bash
oni-spt run --config configs/fvp-validation.yaml --max-files 1 --max-frames 50
oni-spt run --config configs/drrm-validation.yaml --max-files 1 --max-frames 50
```

Each stage can be run independently:

```bash
oni-spt preprocess --config CONFIG.yaml --max-frames 50
oni-spt segment --config CONFIG.yaml
oni-spt detect --config CONFIG.yaml
oni-spt track --config CONFIG.yaml
```

Run the classical analysis on a new trajectory file:

```bash
oni-spt analyze-aio tracks.csv --frame-interval 0.03
oni-spt concat-aio SPT_results_AIO-*.csv --output condition.csv
oni-spt analyze-saspt condition.csv --frame-interval 0.03
```

The historical `cli_reformat-TrackMate-tracks.py` remains available for old
TrackMate exports.

## Acquisition profiles

- `fvp-*.yaml`: the 428-pixel frame is one full-width SPT channel.
- `drrm-*.yaml`: only 856-pixel frames are accepted and split into a 428-pixel
  left Hoechst channel and 428-pixel right SPT channel. The five mistaken
  428-pixel files are left in primary data and listed in `unused_files.csv`.
- `dual-spt-example.yaml`: both halves are independently filtered, detected,
  and tracked; normalized channel MIPs are combined for a shared mask.

## Output layout

```text
output/
  input_manifest.csv
  unused_files.csv
  run_config.resolved.yaml
  provenance.json
  stage_status.json
  channels/       # raw standardized/split stacks and metadata sidecars
  bandpass/       # filtered SPT stacks
  mips/           # raw per-channel MIPs
  segmentation/   # Cellpose source, masks, region table, overlay
  detections/     # canonical per-spot tables
  tracks/         # canonical trajectories, legacy CSVs, track summaries
  qc/             # detection and tracking overlays/histograms
```

Coordinates remain in channel-local pixels. Pixel size and frame interval are
recorded from ONI metadata; the supplied profiles validate 0.117 um/pixel and
30/100 ms acquisitions. Channel registration is assumed and no transform is
applied.
