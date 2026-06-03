# BulkFluoRDF Final Main Workflow

## Summary
BulkFluoRDF now uses a per-hub, annular, physical-bin RDF instead of the earlier
Sofi summed-source nucleus fit. The summed-source model is retained only as
historical context because it assumes H3K27ac pixel intensity can be explained by
the summed radial contribution of all SPEN hubs in a nucleus. That assumption is
not appropriate for this SPEN/H3K27ac use case.

## Default Parameters
```yaml
rdf:
  mode: per_hub_annular
  radius_nm: 3000
  bin_width_nm: 100
  bin_step_nm: 50
  normalization: local_intensity_mean
  tail_normalization:
    enabled: false
  aggregation:
    plot_column: h3k27ac_rdf_local_norm

hub_filter:
  enabled: true
  metric: spotiflow_intensity
  threshold_source: nucleus_spen_median_plus_std
  std_multiplier: 2.0
  hard_threshold:
    enabled: true
    source: negative_control_result_dirs
    control_dirs:
      - /Volumes/guttman/users/gmgao/Imaging_ProcessedData/SPEN/20260413_ONI-gmgao-SPEN_H3K27ac_dSTORM_SACD/processed-SACD/Aux1h
      - /Volumes/guttman/users/gmgao/Imaging_ProcessedData/SPEN/20260413_ONI-gmgao-SPEN_H3K27ac_dSTORM_SACD/processed-SACD/Aux24h
    metric: spotiflow_intensity
    statistic: quantile
    quantile: 0.95
```

## RDF Calculation
For each retained SPEN hub, annular bins are generated in physical units:
`0-100 nm`, `50-150 nm`, `100-200 nm`, continuing until `3000 nm`.

For each hub and annular bin:
- pixels are clipped to the same retained nucleus mask
- `pixel_count` is the number of valid nucleus pixels in that annulus
- `local_h3k27ac_mean` is the mean H3K27ac intensity in the max-radius local
  reference region around that hub, clipped to the same nucleus
- `h3k27ac_rdf_local_norm = h3k27ac_sum_intensity / (local_h3k27ac_mean * pixel_count)`

The same calculation is also performed on the SPEN object channel:
- `local_spen_mean`
- `spen_rdf_local_norm`

SPEN-channel RDF is a positive control and should generally show enrichment near
SPEN hub centers.

## Hub Filtering
All nuclear Spotiflow detections are written to `hub_properties.csv`. The default
main hub filter is:

```text
keep_hub = spotiflow_intensity >= nucleus_spen_median + 2 * nucleus_spen_std
```

When negative-control hard filtering is enabled, the final threshold is:

```text
final_threshold = max(nucleus_spen_median + 2 * nucleus_spen_std, pooled_Aux1h_Aux24h_spotiflow_intensity_q95)
keep_hub = spotiflow_intensity >= final_threshold
```

The `compare-hub-filters` entrypoint still produces median + 1/2/3 STD
comparison folders under `results/hub_filter_comparison/`.

## Aggregation
`aggregated_rdf_summary.csv` is aggregated directly across retained hubs, not
first across nuclei. It includes:
- H3K27ac RDF mean, STD, and SEM
- SPEN positive-control RDF mean, STD, and SEM
- SEM N as the number of contributing SPEN hubs per bin

`rdf_aggregate.png` plots:
- H3K27ac RDF and SPEN positive-control RDF
- all individual hub RDF curves in the background using a shared intensity color
  palette
- aggregate mean curves with SEM error bars, where SEM is hub-level STD divided
  by the square root of the number of contributing SPEN hubs
- y-axis clamped to `[0.5, 2]` if the plotted range would otherwise exceed it
- x coordinates set to each bin's left bound, with the x-axis ending at the last
  left bound

## QC
The per-FOV QC overlay keeps the established montage layout:
- each retained nucleus is shown as a paired SPEN/H3K27ac crop
- nucleus outline is shown on both channels
- retained hubs are hollow circles
- the RDF panel shows per-hub H3K27ac RDF curves plus the mean curve

## Test Plan
```bash
conda run -n smlm python -m py_compile BulkFluoRDF.py
conda run -n smlm python BulkFluoRDF.py check
conda run -n smlm python BulkFluoRDF.py run --config config.yaml
conda run -n smlm python BulkFluoRDF.py compare-hub-filters --config config.yaml
```

Validate that:
- max bin end is `3000 nm`
- bin width is `100 nm`, bin step is `50 nm`
- main outputs use median + 2 STD hub filtering plus the configured Aux control
  hard threshold
- `hub_rdf_results.csv` contains H3K27ac and SPEN local-intensity normalized RDF
  columns
- `aggregated_rdf_summary.csv` contains H3K27ac/SPEN mean, STD, SEM, `n_hubs`,
  `n_nuclei`, and `sem_n` aggregate columns, with `sem_n == n_hubs`
- aggregate and QC plots use RDF terminology and the final axis ranges
