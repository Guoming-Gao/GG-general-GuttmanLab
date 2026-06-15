# SNPaware RNA Primer Designer

Design SNP-aware PCR primers for mature RNA-derived cDNA from a gene symbol.

The pipeline resolves local mm10 refGene transcript models, builds a strand-aware
spliced cDNA template, and designs primers around exon-exon junctions that are
common to curated multi-exon RefSeq transcripts. It preserves the KOveri primer
quality workflow: Primer3 design, B6/Cast informative SNP coverage, primer SNP
exclusion, BLAST specificity annotation, 3-prime GC clamp, thermodynamic checks,
and strict-to-relaxed output tiers.

## Quick Start

For interactive use, open `SNPaware_RNA_primer_pipeline.ipynb`. The notebook
keeps all user-facing setup parameters in the first code cell, including gene
list, output folders, local resource paths, Primer3 settings, SNP thresholds,
junction thresholds, and BLAST settings. It also runs multi-gene designs with a
`rich.progress` progress bar. There is no separate `config.py` file; notebook
runs are configured directly in the notebook.

For command-line use, run from this directory with the `bioinfo` conda
environment. CLI defaults are embedded in `design_snpaware_rna_primers.py` and
can be overridden with command-line arguments:

```bash
conda run -n bioinfo python design_snpaware_rna_primers.py \
  --gene Xist \
  --output-dir "/Users/gmgao/Dropbox/Caltech_PostDoc_GuttmanLab/constructs_primers_FISHprobes/qPCR LibAmp PCR primers/AmpliconSeq-Xist-crossjunction"
```

Batch command-line runs are also supported:

```bash
conda run -n bioinfo python design_snpaware_rna_primers.py \
  --genes Xist,Tsix \
  --output-dir "/Users/gmgao/Dropbox/Caltech_PostDoc_GuttmanLab/constructs_primers_FISHprobes/qPCR LibAmp PCR primers"
```

Batch runs write one subfolder per gene and a `SNPaware_RNA_batch_summary.csv`
file in the output root.

## Inputs

- Gene symbol, for example `Xist`.
- Local mm10 refGene GTF.
- mm10/GRCm38 genome FASTA and BLAST database.
- Mouse Genomes Project VCF with B6 and Cast samples.

Default local resources:

- `/Volumes/guttman/genomes/mm10/annotation/mm10.refGene.gtf.gz`
- `/Volumes/guttman/genomes/mm10/fasta/mm10.fa`
- `/Volumes/guttman/genomes/mm10/fasta/blastdb/mm10_blastdb`
- `/Volumes/guttman/genomes/mm10/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz`

External tools default to `/opt/miniconda3/envs/bioinfo/bin`.

## Notebook Setup

The first notebook code cell is the main user setup surface. Edit these values
there before running:

- `GENE_NAMES`: one or more genes to design.
- `OUTPUT_ROOT`: root folder for batch outputs.
- `OUTPUT_DIR_BY_GENE`: optional exact output folders for specific genes.
- local GTF, FASTA, BLAST database, and VCF paths.
- junction, SNP, amplicon, Primer3, thermodynamic, and BLAST thresholds.

For multiple genes, the notebook writes each gene to its own folder unless an
exact folder is listed in `OUTPUT_DIR_BY_GENE`. It also writes
`SNPaware_RNA_batch_summary.csv` in `OUTPUT_ROOT`.

## Design Logic

1. Parse all GTF transcript and exon records for the requested gene.
2. Keep curated `NM_`/`NR_` multi-exon transcripts for common-junction discovery.
3. Select the longest curated transcript as the spliced cDNA design template.
4. Compute exon-exon junctions common to all curated multi-exon transcripts.
5. Build the mature transcript sequence in 5-prime-to-3-prime orientation,
   including minus-strand genes such as `Xist`.
6. Query B6/Cast-different SNPs in exonic sequence represented in the template.
7. First design primer pairs requiring at least one primer to cross a common
   junction.
8. If no orderable crossing-primer pair passes, fall back to cDNA amplicons that
   span a common junction.
9. Exclude SNP positions from primer binding sites and post-filter any primer
   pair that overlaps an informative SNP.
10. Annotate primer specificity using short-query BLAST.
11. Apply the same strict-to-relaxed KOveri quality tiers and report only the
    strictest available non-failed tier.

## Junction Rules

For the preferred crossing-primer mode, a primer must overlap a common exon-exon
junction with at least 6 bases on each side of the junction and at least 4 bases
anchoring the primer 3-prime side. The 3-prime side is orientation-aware for
left and right primers. These thresholds can be tuned with
`--min-junction-overlap` and `--min-junction-3p-anchor`.

For fallback mode, the cDNA amplicon must span a common junction. Fallback rows
are labeled as `junction_spanning_fallback`.

## Outputs

- `SNPaware_RNA_transcripts.csv`: transcript models, selected template, and
  excluded transcript reasons.
- `SNPaware_RNA_common_junctions.csv`: common junctions with genomic and cDNA
  coordinates.
- `SNPaware_RNA_primer_candidates.csv`: all deduplicated primer pairs.
- `SNPaware_RNA_primer_top_candidates.csv`: ranked top pairs for ordering.
- `SNPaware_RNA_primers.fasta`: primer sequences.
- `SNPaware_RNA_amplicons.bed`: common junctions, amplicon exon segments, and
  primer genomic segments.
- `SNPaware_RNA_report.html`: consolidated report.
- `SNPaware_RNA_summary.txt`: run summary and best ranked primer pair.
- `SNPaware_RNA_batch_summary.csv`: batch-level summary written only for
  multi-gene notebook or CLI runs.

## Coordinates

cDNA coordinates are 1-indexed on the mature spliced transcript. Genomic
coordinates are 1-indexed inclusive except the BED file, which uses standard
0-based starts and 1-based ends. Primer sequences are reported 5-prime to
3-prime. The right primer is the reverse primer sequence to order.
