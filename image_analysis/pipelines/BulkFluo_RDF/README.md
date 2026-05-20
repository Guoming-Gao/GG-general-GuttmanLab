# BulkFluoRDF

Lightweight per-hub RDF pipeline for paired SPEN/H3K27ac SACD images.

## Terms

- `object channel`: SPEN `*-SACD-left.tif`; Spotiflow detects SPEN hub centers.
- `intensity channel`: H3K27ac `*-SACD-right.tif`; Cellpose segments nuclei and RDF measures H3K27ac intensity around hubs.
- `nucleus`: labeled mask defining valid pixels for local per-hub RDF.

The main RDF method is per-hub annular RDF. Each retained SPEN hub gets its own
H3K27ac radial profile using overlapping physical bins clipped to the hub's
parent nucleus. The older Sofi summed-source model was tested but is not the
default because it assumes H3K27ac pixels are explained by summed radial
contributions from all nearby SPEN hubs, which is not appropriate when one
nucleus contains many SPEN hubs and H3K27ac has other biological determinants.

## Environment

Run everything in `smlm`:

```bash
conda run -n smlm python BulkFluoRDF.py check
conda run -n smlm python BulkFluoRDF.py pair --config config.yaml
conda run -n smlm python BulkFluoRDF.py run --config config.yaml
conda run -n smlm python BulkFluoRDF.py compare-hub-filters --config config.yaml
conda run -n smlm python BulkFluoRDF.py compare-detectors --config config.yaml
conda run -n smlm python BulkFluoRDF.py screen-spotiflow --config config.yaml
```

Installed/required key packages: `cellpose`, `spotiflow`, `torch`, `numpy`, `pandas`, `scipy`, `scikit-image`, `matplotlib`, `tifffile`, `pyyaml`.

`config.yaml` currently points at the full `undiff` SACD dataset and uses saved
masks/spots when available. CellposeSAM is called in the same nuclei-only style as
`/Users/gmgao/GGscripts/GG-general-GuttmanLab/image_analysis/segmentation/gui_CellposeSAM_run.py`.
For SACD images, the default Cellpose input preprocessing uses percentile
background subtraction/rescaling rather than max-only normalization, because the
H3K27ac SACD intensity channel has small floating-point values and bright
outliers. Set `cellpose.preprocess.method: max` to reproduce the GUI
normalization exactly.

Tune `cellpose.diameter` first if nuclei are missed or fragmented. The default
local smoke-test config uses `diameter: 120` with `downsample: 0.5`; the effective
diameter passed to Cellpose is scaled with the downsampled image.

Distance defaults are `pixel_size_nm: 58.5`, `rdf.radius_nm: 1000`,
`rdf.bin_width_nm: 100`, and `rdf.bin_step_nm: 50`. Bins are overlapping, e.g.
`0-100 nm`, `50-150 nm`, `100-200 nm`.

Hub filtering compares Spotiflow's `object_intensity` value, reported as
`spotiflow_intensity`, against the median plus one STD of SPEN pixel intensities
inside that hub's own nucleus. This uses the raw SPEN nucleus pixel distribution,
not the detected spot intensity distribution. Use `compare-hub-filters` to
generate median + 1/2/3 STD comparison result folders.

QC overlays use robust percentile contrast controlled by
`qc.contrast_lower_percentile` and `qc.contrast_upper_percentile`. One PNG is
written per SACD pair. Each overlay shows retained nuclei as crop pairs: SPEN
SACD with nucleus outline and hub circles, H3K27ac SACD with the same outline
and hubs, plus a right-side per-hub RDF panel. Gray circles are rejected
Spotiflow detections; colored circles are retained SPEN hubs.
Individual per-nucleus versions of the same SPEN/H3K27ac crop overlays are
also written to `results/qc_overlays/overlays/`, with titles and 1 um scale
bars.

## Files

- `BulkFluoRDF.py`: reusable functions and CLI.
- `BulkFluoRDF_pipeline.ipynb`: notebook wrapper for running and QC.
- `config.yaml`: paths and parameters.
- `results/`: generated outputs.

Main outputs:

- `results/paired_inputs.csv`
- `results/nucleus_masks/`
- `results/nucleus_qc.csv`
- `results/spotiflow_cache/`
- `results/hub_properties.csv`
- `results/hub_rdf_results.csv`
- `results/aggregated_rdf_summary.csv`
- `results/rdf_aggregate.png`
- `results/qc_overlays/`
- `results/qc_overlays/overlays/`

Set `cellpose.use_existing_masks` or `spotiflow.use_existing_spots` to `true` to reuse saved outputs.

Current RDF defaults use Spotiflow `probability_threshold: 0.4`,
`min_distance: 1`, per-hub annular RDF, local-intensity H3K27ac normalization,
plain per-bin mean/STD aggregation across retained hubs, nucleus area/edge
filtering, and median + 2 STD SPEN hub filtering.

## Detector Comparison

`compare-detectors` sweeps Spotiflow and DoG parameters on SPEN SACD images and
measures every retained nuclear detection with the same raw-image features:
peak intensity, local background, background-corrected integrated intensity,
punctum area, and equivalent diameter. This is intended to test whether
Spotiflow enriches for PSF-like detections rather than the largest and strongest
SPEN hub puncta.
The default Spotiflow comparison keeps `min_distance: 3` fixed and sweeps only
probability threshold.

Outputs are written to `results/detector_comparison/`:

- `detections.csv`
- `summary.csv`
- `top_hubs.csv`
- QC plots for distributions and top-hub overlays

## Spotiflow Screen

`screen-spotiflow` runs a focused Spotiflow threshold screen for visual QC.
The main RDF pipeline now uses the selected screen result directly:
Spotiflow threshold `0.4`, nucleus area/edge filtering, and a per-nucleus hub
filter where `spotiflow_intensity >= median(SPEN pixels inside that nucleus) + 2 * std(SPEN pixels inside that nucleus)`.

## Hub Filter Comparison

Compare median + 1/2/3 STD hub filters without rerunning Cellpose or Spotiflow:

```bash
conda run -n smlm python BulkFluoRDF.py compare-hub-filters --config config.yaml
```

Outputs are written under `results/hub_filter_comparison/`. Each condition
subfolder name includes the runtime retained hub count, and the comparison root
also contains `filter_comparison_summary.csv` and
`filter_comparison_rdf_overlay.png`.

## Per-Hub RDF Definition

For each retained hub `h` in nucleus `n`, BulkFluoRDF measures H3K27ac pixels
inside overlapping annular bins clipped to the same nucleus mask. For a bin
`[r0, r1)`, the output is:

```text
local_reference(h) = pixels in nucleus n with distance(pixel, h) < rdf.radius_nm
local_mean(h) = mean H3K27ac intensity in local_reference(h)
observed_sum(h, r0, r1) = sum H3K27ac pixels in nucleus n with r0 <= distance(pixel, h) < r1
expected_sum(h, r0, r1) = local_mean(h) * pixel_count(h, r0, r1)
rdf_local_norm(h, r0, r1) = observed_sum / expected_sum
```

`pixel_count` is the actual number of valid nucleus pixels in that clipped
annular bin. This provides ring-area/boundary normalization while the denominator
uses the hub-local max-radius intensity reference instead of the whole nucleus. Aggregated
plots show the per-bin mean and STD across retained hubs, with individual hub
curves in the background. The same RDF calculation is also run on the SPEN object
channel as a positive control.
