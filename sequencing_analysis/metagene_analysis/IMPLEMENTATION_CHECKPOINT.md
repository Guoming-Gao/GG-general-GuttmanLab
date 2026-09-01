# SPEN artifact-validation implementation checkpoint

Last updated: 2026-08-20

## Completed

- Renamed the external results root to `SPEN-2026Aug-metagene_analysis`.
- Verified SHA-256 identity for all 511 pre-existing files across the rename.
- Reorganized prior outputs under `prior_exploration/` and created the approved
  reader-facing results structure.
- Moved validated No-Dox SPEN CLAP matrices, representative-transcript models,
  and untreated total-RNA normalized counts into `supporting_data/`.
- Updated `config/config.yaml` and `README.md` to the new root.
- Installed OpenJDK 25 in the `bioinfo` environment for compatibility with the
  historical `CLAPAnalysis.jar`.
- Began `artifact_validation.py`: expression inputs, transcript feature geometry,
  CA filtering/assignment, absolute peak counts, bootstrapped correlations,
  length-adjusted logistic/NB models, and length-stratified sensitivities exist.
- Completed and cached all SPEN-only artifact diagnostics, CA peak-count figures,
  length-adjusted tables, and Input/pseudocount diagnostics under `results/02_*`
  and `results/03_*`.
- The initial recovery test suite passed 25 tests; the final expanded suite
  passes 26 tests after adding zero-exclusion coverage.
- Completed and cached all HUHCLAP BAM counts, TMM-normalized matrices, and
  replicate/pooled CA calls for HnRNPC, PTBP1, and GFP. These expensive files
  must not be rebuilt.
- User rejected zero-dominated CA plots. All such PNG/PDF files were deleted.
  Code now excludes zero-peak gene features from CA distributions,
  correlations, and conditional count models, while retaining zero fractions
  only as QC in tables. CA is explicitly secondary rather than headline proof.
- Regenerated the sole authoritative self-contained report at
  `report/SPEN_expression_dependence_report.html`. The executive conclusion now
  reports the observed component correlations and treats denominator/coverage
  bias as the leading explanation, with appropriate limitations.
- Full automated test suite passes: 26 tests. HUHCLAP full-file SHA-256
  verification reports all 24 BAM/index source files unchanged. Final Snakemake
  dry run reports that all requested files are present and up to date.

## Next exact step

Implementation is complete. If future revision is requested, start from the
cached matrices and CA files; do not rerun BAM counting or CA scoring. Modify
the plotting/report code, force only `artifact_validation`, run tests, and end
with a dry run.

## Reuse / do not rebuild

- Rename provenance:
  `SPEN-2026Aug-metagene_analysis/workflow_records/provenance/`
- Validated SPEN matrices:
  `SPEN-2026Aug-metagene_analysis/supporting_data/validated_SPEN_CLAP_matrices/`
- Untreated normalized counts:
  `SPEN-2026Aug-metagene_analysis/supporting_data/untreated_total_RNA/`
- Representative annotation:
  `SPEN-2026Aug-metagene_analysis/supporting_data/annotation/`

No HUHCLAP source file has been modified. Cached HUHCLAP artifacts are under
`supporting_data/HUHCLAP/{HnRNPC,PTBP1,GFP}/`; pooled CA outputs are as large as
~1.2 GB and are intentionally retained for reproducibility.
