# KOveri Primer Designer

Design PCR primers for CRISPR KO genotyping and amplicon sequencing in TX cells.

The pipeline takes user-provided CRISPR gRNA sequences, finds their exact mm10
genomic locations, and designs primer pairs that flank the full edited interval
plus B6/Cast informative SNPs. Primers are rejected if they overlap informative
B6/Cast SNPs, so amplification should not be biased toward one allele.

## Quick Start

Run from this directory with the `bioinfo` conda environment:

```bash
conda activate bioinfo
python design_koveri_primers.py \
  --output-dir "/Users/gmgao/Dropbox/Caltech_PostDoc_GuttmanLab/constructs_primers_FISHprobes/qPCR LibAmp PCR primers/AmpliconSeq-KOveri-RNF12"
```

By default, the input FASTA is expected to be in the same folder as the result
files. For the RNF12 run above, the package reads:

```text
/Users/gmgao/Dropbox/Caltech_PostDoc_GuttmanLab/constructs_primers_FISHprobes/qPCR LibAmp PCR primers/AmpliconSeq-KOveri-RNF12/RNF12_gRNAs.fa
```

The same workflow is available in `KOveri_primer_pipeline.ipynb` for interactive
tuning.

## Inputs

- gRNA FASTA file. RNA `U` bases are converted to DNA `T`.
- mm10/GRCm38 genome FASTA and BLAST database.
- Mouse Genomes Project VCF with B6 and Cast samples.

Default local resources match `../smfish-like-rt-probe-designer`:

- `/Volumes/guttman/genomes/mm10/fasta/mm10.fa`
- `/Volumes/guttman/genomes/mm10/fasta/blastdb/mm10_blastdb`
- `/Volumes/guttman/genomes/mm10/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz`

## Design Logic

1. Locate each gRNA with `blastn -task blastn-short`.
2. Keep only full-length exact mm10 hits.
3. Select the smallest same-chromosome cluster containing one hit for each gRNA.
4. Treat the min/max guide coordinates as the protected edit interval.
5. Query B6/Cast-different SNPs around that interval with `pysam`.
6. Try increasing flanking windows until Primer3 finds products covering all
   guides and the requested SNP count.
7. Exclude SNP positions from primer binding sites and post-filter any primer
   pair that overlaps an informative SNP.
8. Annotate primer specificity using short-query BLAST.
9. Apply post-Primer3 genotyping filters and select top candidates from the
   strictest passing tier only.

If the ideal short product cannot include at least three SNPs, the pipeline
expands the amplicon window before relaxing the SNP count.

## Post-Primer3 Filters

`KOveri_primer_top_candidates.csv`, `KOveri_primers.fasta`, `KOveri_amplicons.bed`,
`KOveri_summary.txt`, and `KOveri_report.html` use only the strictest available
non-failed filter tier. If one strict primer pair passes, the package reports
that one pair instead of padding the top list with relaxed lower-quality pairs.
Candidates that fail hard selectable rules are never orderable.

Relaxation order is fixed:

1. strict
2. relaxed amplicon size
3. relaxed specificity
4. relaxed GC/Tm
5. failed candidates retained only in the full candidates CSV

Hard selectable filters require GC 40-60%, Tm 55-65 C, pair Tm delta <=5 C,
3-prime G/C clamp for both primers, no primer SNP overlap, guide coverage,
minimum SNP count, no homopolymer runs of 4 or more, no quad-G, no alternating
dinucleotide repeat of 6 or more bases, no inter-primer complement run over
3 bases, ntthal dG > -9 kcal/mol for hairpins and dimers, and configured
maximum amplicon size. The 3-prime G/C clamp is not relaxed.

Soft preferences are relaxed only if no stricter tier exists: 400-700 bp
amplicons relax first to the configured maximum, strict BLAST specificity relaxes
second to annotated non-strict specificity, and GC/Tm relax modestly only as a
final rescue. Primer length is only checked against the broad valid 18-30 nt
range and is not used as a tier preference.

## Outputs

- `KOveri_guide_hits.csv`: all exact guide hits and the selected target hits.
- `KOveri_primer_candidates.csv`: all deduplicated primer pairs.
- `KOveri_primer_top_candidates.csv`: ranked top pairs for ordering.
- `KOveri_primers.fasta`: primer sequences.
- `KOveri_amplicons.bed`: protected interval, amplicons, and primer coordinates.
- `KOveri_report.html`: consolidated report with a genomic coordinate figure.
- `KOveri_summary.txt`: run summary and best ranked primer pair.

## Important Coordinates

All reported genomic coordinates are 1-indexed inclusive except the BED file,
which uses standard 0-based starts and 1-based ends. Primer sequences are
reported 5' to 3'. The right primer is the reverse primer sequence to order.
