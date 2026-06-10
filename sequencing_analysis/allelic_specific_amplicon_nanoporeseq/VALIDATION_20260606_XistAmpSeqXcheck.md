# Validation Note: 20260606 XistAmpSeqXcheck

This note captures the validated Nanopore demultiplexing and plotting behavior from the prior Codex session.

## Dataset

Data folder:

`/Volumes/guttman/users/gmgao/Data_seq/20260606-XistAmpSeqXcheck`

Raw FASTQ:

`DRPD44_1_XistAmpSeq-Xcheck-8wControl.fastq`

Sample map:

`barcode_sample_map.csv`

The sample map uses `SampleName,i5,i7`. In `00_demultiplex_inline_barcodes.py`, `i5` is treated as the P5-side barcode and `i7` is reverse-complemented into the internal P7 orientation.

## Validated Demux Outcome

Using `--max-mismatches 1`:

- Total reads: `10774`
- Assigned reads: `5765`
- Ambiguous reads: `0`
- Blank sample matches: `0`
- Unrecognized structure: `4700`
- Unmatched barcode: `309`

Expected generated files:

- `fastq_by_sample/*.fastq`: 9 files
- `demux_summary.csv`
- `demux_audit/unassigned.fastq`

## Validated Downstream Outputs

Under `results/mode_1/`:

- `aligned/*.sorted.bam`: 9 files
- `aligned/*.sorted.bam.bai`: 9 files
- `xist_snps.txt`: 101 SNPs plus header
- `quantification/*_quant.csv`: 9 files
- `quantification/allele_quantification_summary.csv`: nonempty
- `allele_barplot_consolidated.png`: regenerated and visually confirmed to show raw Cast/B6 read counts inside bars

## Plot Regeneration Command

```bash
cd /Users/gmgao/GGscripts/GG-general-GuttmanLab/sequencing_analysis
/opt/miniconda3/envs/bioinfo/bin/python allelic_specific_amplicon_nanoporeseq/05_plot_allele_summary.py   --results-dir /Volumes/guttman/users/gmgao/Data_seq/20260606-XistAmpSeqXcheck/results/mode_1   --sample-map /Volumes/guttman/users/gmgao/Data_seq/20260606-XistAmpSeqXcheck/barcode_sample_map.csv
```

Output:

`/Volumes/guttman/users/gmgao/Data_seq/20260606-XistAmpSeqXcheck/results/mode_1/allele_barplot_consolidated.png`

The plot defaults to raw read counts inside bars. Add `--hide-raw-counts` only for percent-only output.

## Environment

Run from the `bioinfo` conda environment:

```bash
conda activate bioinfo
```

Validated explicit interpreter:

```bash
/opt/miniconda3/envs/bioinfo/bin/python
```
