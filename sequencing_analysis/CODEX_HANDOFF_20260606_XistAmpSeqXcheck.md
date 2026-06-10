# Codex Handoff: 20260606 Xist Amplicon Sequencing Work

Generated from the prior Codex chat on 2026-06-10 after the workspace was renamed to `sequencing_analysis`.

## Workspace State

The sequencing pipeline workspace now lives at:

`/Users/gmgao/GGscripts/GG-general-GuttmanLab/sequencing_analysis`

The two Dox-enabled amplicon packages were renamed to:

- `allelic_specific_amplicon_nanoporeseq`
- `allelic_specific_amplicon_shortreadseq`

Both packages are documented to run in the `bioinfo` conda environment. Use:

```bash
conda activate bioinfo
```

For command-line scripts that need an explicit interpreter, the validated interpreter was:

```bash
/opt/miniconda3/envs/bioinfo/bin/python
```

## Main Dataset Processed

Target dataset:

`/Volumes/guttman/users/gmgao/Data_seq/20260606-XistAmpSeqXcheck`

Raw FASTQ:

`/Volumes/guttman/users/gmgao/Data_seq/20260606-XistAmpSeqXcheck/DRPD44_1_XistAmpSeq-Xcheck-8wControl.fastq`

Sample map:

`/Volumes/guttman/users/gmgao/Data_seq/20260606-XistAmpSeqXcheck/barcode_sample_map.csv`

The sample map was changed by the user to the lab naming convention with `SampleName,i5,i7`. The Nanopore demultiplexer treats `i5` as the P5-side barcode and reverse-complements `i7` into the internal P7 orientation.

## Nanopore Package Changes

Package path:

`/Users/gmgao/GGscripts/GG-general-GuttmanLab/sequencing_analysis/allelic_specific_amplicon_nanoporeseq`

Added/updated workflow pieces:

- `00_demultiplex_inline_barcodes.py`
  - Accepts `SampleName,i5,i7` and `SampleName,P7_Barcode,P5_Barcode` maps.
  - Supports forward and reverse read orientations.
  - Uses Illumina adapter anchors around inline barcodes.
  - Searches small barcode offset windows.
  - Allows barcode edit distance via `--max-mismatches`, validated with `--max-mismatches 1`.
  - Writes demultiplexed sample FASTQs, `demux_summary.csv`, and `demux_audit/unassigned.fastq`.

- `run_nanopore_pipeline.py`
  - Supports multiplexed mode via `--fastq` plus `--sample-map`.
  - Supports existing per-sample FASTQ mode via `--fastq-dir` or `--data-dir`.
  - Uses output layout:
    - `fastq_by_sample/`
    - `demux_summary.csv`
    - `demux_audit/`
    - `results/mode_1/`

- `05_plot_allele_summary.py`
  - Creates `results/mode_1/allele_barplot_consolidated.png`.
  - Uses sample-map order when `--sample-map` is provided.
  - Matches the short-read visualization style: horizontal stacked bars, Cast `#9a3324`, B6 `#00274c`, `fontsize=25`, bold Cast percent labels at x=101, hidden top/right spines, legend below.
  - Raw Cast and B6 read counts are shown inside bars by default.
  - Use `--hide-raw-counts` for a percent-only figure.

- `README.md`, `pipeline.ipynb`, and `visualization.ipynb`
  - Updated for `bioinfo` environment.
  - Document multiplexed and already-demultiplexed workflows.
  - Document the `i5/i7` convention.
  - Visualization notebook exposes `SHOW_RAW_COUNTS = True`.

## Short-Read Package Changes

Package path:

`/Users/gmgao/GGscripts/GG-general-GuttmanLab/sequencing_analysis/allelic_specific_amplicon_shortreadseq`

Updates:

- Added `README.md` stating the `bioinfo` environment requirement.
- Updated `pipeline.ipynb` title and environment wording.
- Updated `visualization.ipynb`:
  - Uses `bioinfo` kernel metadata.
  - Adds `SHOW_RAW_COUNTS = True`.
  - Turns on the previously commented raw Cast/B6 count labels inside horizontal stacked bars.

Existing caveat:

- Some short-read scripts/notebooks still contain pre-existing hard-coded paths for `/Volumes/guttman/users/gmgao/Data_seq/20260207-DoxSeqAllRepsMultiplexed`. This was inspected but not generalized. The hard-coded large short-read alignment was not rerun.

## Validated 20260606 Commands

Regenerate the final Nanopore plot with raw read counts:

```bash
cd /Users/gmgao/GGscripts/GG-general-GuttmanLab/sequencing_analysis
/opt/miniconda3/envs/bioinfo/bin/python allelic_specific_amplicon_nanoporeseq/05_plot_allele_summary.py   --results-dir /Volumes/guttman/users/gmgao/Data_seq/20260606-XistAmpSeqXcheck/results/mode_1   --sample-map /Volumes/guttman/users/gmgao/Data_seq/20260606-XistAmpSeqXcheck/barcode_sample_map.csv
```

The generated plot is:

`/Volumes/guttman/users/gmgao/Data_seq/20260606-XistAmpSeqXcheck/results/mode_1/allele_barplot_consolidated.png`

## Validation Results From Prior Run

Demultiplexing on the 20260606 FASTQ reproduced the expected counts:

- Total reads: `10774`
- Assigned reads: `5765`
- Ambiguous reads: `0`
- Blank sample matches: `0`
- Unrecognized structure: `4700`
- Unmatched barcode: `309`

Pipeline output checks:

- Demultiplexed FASTQs: `9`
- Sorted BAMs: `9`
- BAM indexes: `9`
- `xist_snps.txt`: `102` lines total, meaning `101` SNPs plus header
- Quantification files: `9`
- Allele quantification summary: nonempty

The final plot was visually inspected and contains raw Cast/B6 read counts inside each bar. Examples visible in the 20260606 plot include:

- `dRNF12_1A7_gDNA`: Cast `541`, B6 `0`, Cast `100.0%`
- `dTsix_1A10_gDNA`: Cast `2`, B6 `782`, Cast `0.3%`
- `Tx_gDNA_positive_ctrl`: Cast `71`, B6 `75`, Cast `48.6%`

Smoke checks passed:

```bash
/opt/miniconda3/envs/bioinfo/bin/python -m py_compile   allelic_specific_amplicon_nanoporeseq/*.py   allelic_specific_amplicon_shortreadseq/*.py

/opt/miniconda3/envs/bioinfo/bin/python   allelic_specific_amplicon_nanoporeseq/05_plot_allele_summary.py --help
```

Note: Matplotlib may warn that `/Users/gmgao/.matplotlib` is not writable and create a temporary cache. This did not block plotting.

## 20260606 Data Folder Cleanup

The local one-off script copies were removed from:

`/Volumes/guttman/users/gmgao/Data_seq/20260606-XistAmpSeqXcheck`

Removed local copies included:

- `00_demultiplex_inline_barcodes.py`
- `01_fastq_quality_metrics.py`
- `02_align_to_genome.py`
- `03_extract_snps_from_vcf.py`
- `04_quantify_alleles.py`
- `run_20260606_pipeline.sh`
- `.DS_Store` files

Kept outputs:

- Raw FASTQ
- `barcode_sample_map.csv`
- `demux_summary.csv`
- `demux_audit/unassigned.fastq`
- `fastq_by_sample/`
- `results/mode_1/`

## Reference Paths Used

Genome FASTA:

`/Volumes/guttman/genomes/mm10/fasta/mm10.fa`

SNP VCF:

`/Volumes/guttman/genomes/mm10/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz`

Xist region:

`chrX:103460373-103483233`

Strains:

- B6: `C57BL_6NJ`
- Cast: `CAST_EiJ`

## Remaining Follow-Up Candidates

- Add a tiny synthetic FASTQ test for Nanopore demux covering forward orientation, reverse orientation, one barcode mismatch, ambiguous barcode, and unassigned read.
- Generalize old short-read hard-coded dataset paths if this package will be reused broadly.
- Consider setting `MPLCONFIGDIR` to a writable directory for cleaner Matplotlib startup in batch runs.
