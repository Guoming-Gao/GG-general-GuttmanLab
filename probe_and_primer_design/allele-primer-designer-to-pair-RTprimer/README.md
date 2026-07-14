# Allele Primer Designer To Pair RT Primer

Design one RNA-sense, gene-specific primer insert for targeted single-cell RNA
library amplification. The insert is intended to pair with an existing
poly-dT/RT-handle primer. Each accepted insert is emitted with one selected
2P handle; the default is `2PUNI`:

- `2PBC`: `CAGACGTGTGCTCTTCCGATCT`
- `2PUNI`: `CCTACACGACGCTCTTCCGATCT`

The workflow is strand-aware and works in mature spliced RNA coordinates. For
each informative B6/Cast SNP, Primer3 searches the 100 nt window immediately
upstream of the SNP in RNA 5' to 3' order. A passing primer must have the target
SNP downstream toward the transcript/poly-dT 3' end and within the configured
window.

## Environment

Use the `bioinfo` conda environment, which must provide Python 3.12, pandas,
Rich, Primer3, samtools, BLAST+, pysam, and Jupyter. The default resource paths
point to the laboratory mm10 RefSeq GTF, genome FASTA, BLAST database, and Mouse
Genomes Project B6/Cast VCF under `/Volumes/guttman/genomes/mm10`.

For the notebook, select the Python kernel from the `bioinfo` environment and
start Jupyter from this package directory.

## Notebook Workflow

Open `allele_primer_design_pipeline.ipynb`. Its first code cell contains all
configuration, including genes, output path, resources, samples, tools, handle,
primer QC limits, SNP search window, BLAST thresholds, and progress settings.
The defaults reproduce the seven-gene X-chromosome production run with `2PUNI`.

Run the cells in order to:

1. Validate resources and external tools.
2. Run the shared package implementation with Rich progress.
3. Load consolidated CSV outputs.
4. Assert handle, sequence, SNP direction/window, overlap, and specificity QC.
5. Display order-ready oligos and a link to the consolidated HTML report.

The notebook overwrites only the package's named output files. It does not
remove unrelated files from the selected output directory.

## CLI Workflow

Run from this directory with the `bioinfo` environment:

```bash
conda run -n bioinfo python design_allele_primers_to_pair_rtprimer.py \
  --genes Xist,Kdm5c \
  --output-dir ./allele_primer_results \
  --handle 2PUNI \
  --skip-blast-specificity
```

Omit `--skip-blast-specificity` for a production run. Use `--handle 2PBC` to
select the alternative handle.

## Parameters

- `genes` / `--genes`: ordered gene symbols. Default CLI genes are
  `Xist,Kdm5c`; the notebook defaults to the seven-gene production set.
- `output directory` / `--output-dir`: root for consolidated files and one
  subfolder per gene.
- `handle` / `--handle`: exactly one of `2PUNI` or `2PBC`; default `2PUNI`.
- `top_snps` / `--top-snps`: maximum accepted SNP targets per gene; default 5.
- `snp_window_toward_polydt` / `--snp-window-toward-polydt`: maximum spliced
  RNA distance from primer 3-prime end toward the SNP; default 100 nt.
- Primer size, Tm, and GC parameters control both Primer3 and post-filtering.
  The notebook additionally exposes GC clamp, self-complementarity,
  homopolymer, quad-G, and Primer3 return-count settings.
- Resource/sample/tool settings select the GTF, FASTA, VCF, BLAST database,
  B6/Cast VCF samples, samtools, Primer3, and BLAST executables.
- BLAST settings control minimum unique alignment length and identity, noise
  alignment length, and the maximum candidates screened per SNP.

## Progress and BLAST Performance

Rich prints `[Gene i/N]` start/completion messages and a progress bar for SNP
targets within the active gene. The SNP bar shows scanned targets, accepted
targets, percentage, and elapsed time, so a long BLAST step is not mistaken for
a stalled run.

With specificity enabled, SNPs are considered from closest to farthest from
the mature RNA 3-prime end. For each SNP, the five best QC-passing Primer3
candidates are checked in one batched BLAST call by default. The workflow stops
once five distinct SNP targets pass or all informative SNPs are exhausted.
Change the bound with `--blast-candidates-per-snp`; larger values improve search
depth but increase runtime. BLAST-disabled runs are useful for smoke testing,
not final ordering decisions.

## Outputs

For batch runs, each gene writes to `<output>/<gene>_allele_primers/`.

The output root also contains `allele_primer_batch_summary.csv`,
`allele_primer_consolidated_top_oligos.csv`, and
`allele_primer_consolidated_report.html` for reviewing all genes together.

- `allele_primer_transcripts.csv`
- `allele_primer_informative_snps.csv`
- `allele_primer_candidate_inserts.csv`
- `allele_primer_top_oligos.csv`
- `allele_primer_oligos.fasta`
- `allele_primer_report.html`
- `allele_primer_summary.txt`

`SNP_Distance_To_RNA_3p` is always the spliced RNA/cDNA distance to the
transcript 3' end, not genomic distance.

A gene is considered successfully analyzed even if fewer than `top_snps`
primers pass. Zero or fewer-than-five oligos means the available informative
SNPs were exhausted under the configured primer and specificity filters; it is
not necessarily a software failure. The batch summary and consolidated report
retain these genes and their counts.

Open `allele_primer_consolidated_report.html` in a browser to compare all genes
and order-ready sequences. Each gene subfolder also contains its own HTML
report and detailed candidate tables.

## Tests

Run the regression suite from this package directory:

```bash
PYTHONDONTWRITEBYTECODE=1 /opt/miniconda3/envs/bioinfo/bin/python \
  -m unittest discover -s tests
```

The tests cover plus/minus strand coordinate mapping, SNP direction toward
poly-dT, spliced RNA 3-prime distances, handle selection and sequence assembly,
and consolidated report behavior.
