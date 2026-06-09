# Allele-Specific Amplicon Short-Read Sequencing Pipeline

Short-read amplicon sequencing workflow for SNP-based B6/Cast allelic quantification.

## Requirements

- **Conda environment**: `bioinfo`
- Run all commands and notebooks from this environment: `conda activate bioinfo`
- **Key tools**: `fastp`, `bowtie2`, `samtools`, `pysam`, `pandas`, `matplotlib`, `numpy`, `seaborn`

The pipeline and visualization notebooks use the `bioinfo` kernel. If Jupyter opens with a different kernel, switch to `bioinfo` before running cells.

## Visualization

`visualization.ipynb` generates the short-read horizontal stacked allele plot and SNP-match heatmaps. The allele plot now shows raw Cast and B6 read counts inside the bars by default via `SHOW_RAW_COUNTS = True`.
