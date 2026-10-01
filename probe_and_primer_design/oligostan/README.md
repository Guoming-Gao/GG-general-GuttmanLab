# Integrated Oligostan Python port

This directory contains the Python Oligostan port and a mouse exon/intron
design workflow. The importable package is `oligostan`; the notebook calls
its core directly.

Use the sibling notebook
[`Mouse_smiFISH_exon_intron_pipeline.ipynb`](../Mouse_smiFISH_exon_intron_pipeline.ipynb)
from a Python environment with Biopython, NumPy, pandas, Rich and Matplotlib. Its BLAST
step invokes the `blastn` executable and mm10 database configured locally;
the notebook can run in any kernel containing those Python packages.
The selected probe sets, hit audit, genomic coverage PDF and manifest are saved
under a new run folder in the configured `output_parent` directory.

## BLAST-first minimum-span add-on

`minimum_span.py` adds compact set selection after the Oligostan design and
BLAST stages. The complete quality-tier-eligible candidate pool is generated
across all annotated tiles, with no default 600-probe intron stop. Every
retained candidate is BLAST checked in batches. Before BLAST, the workflow
rejects any probe sequence with five or more consecutive A, T, C, or G bases,
matching the neighboring RT-probe designer's max-four-base rule. This rule
applies to every quality tier; each candidate records its longest run, and
the run manifest records the limit. The add-on receives the **full**
BLAST-passing table directly, not a 40-probe preselection. For each available
gene/region, it compares all consecutive 30-probe windows in genomic order
and takes the smallest inclusive genomic span. Filter tier and Oligostan
score only break equal-span ties. Fewer than 30 passing probes are reported as
a shortfall; unverified probes are never used. The core Oligostan R parity
claim does not extend to this Python-only selection mode.

`coverage_report.py` writes `coverage_summary.csv`, `coverage_report.pdf`, and
one PNG per set with mm10 genomic-coordinate ticks. The seven-gene workflow
contains 12 SPEN-target exon/intron sets and a Malat1 exon-only control;
Malat1 has no intron in the chosen mm10 RefGene annotation.

Copy `mouse_smifish.example.json` to `.smifish-local.json` in the directory
containing this package, then replace every placeholder with a local path.
The ignored config supplies `gtf`, `fasta`, `blast_db`, `blastn`, and
`output_parent`; `spen_fasta` is optional for the validator. The reference
FASTA needs a `.fai` index and the BLAST database needs a `.njs` file.
The notebook, CLI (`python -m oligostan.run_mouse_smifish`), and audit command
(`python -m oligostan.audit_run RUN_DIR`) read this config by default.
Each command also accepts `--config`.
The CLI's `--max-homopolymer-length` defaults to 4; a new run records the
chosen value in `manifest.json` and the audit checks it against every probe.

The validation script checks the supplied human RNU1 fixture and can compare
mouse Spen to R output. The optional `--source-checkout` parameter compares a
second Python checkout when needed:

```bash
PYTHONPATH=probe_and_primer_design python -m oligostan.validate_spen \
  --spen-fasta /path/to/Spen_mm10.fa \
  --report /path/to/spen_parity.json
```

The matched-settings mouse Spen run passed the original R comparison for all
1,821 candidate rows and all 725 probes passing the active GC/PNAS filters.
See `PARITY_PLAN.md` for the comparison settings and remaining gates. The
original R script's automatic dG37 optimization and RepeatMasker behavior have
not been established as equivalent.

## Design choices

- Mouse GRCm38/mm10 RefGene annotation and genome are supplied through the
  local config; every output contains the transcript accession and
  genomic position. The longest spliced RefGene isoform is selected per gene.
- Exon candidates stay within one selected-transcript exon. Intron regions
  exclude all annotated exons for the gene, including other isoforms.
- The Python core reverse-complements the target interval to generate an
  antisense oligo. The workflow independently checks that every candidate is
  the reverse complement of its recorded target sequence.
- The project package's GC 0.4–0.6 and PNAS 1/2/4 settings are the baseline.
  Additional candidates receive documented relaxed quality tiers; all tiers
  retain the max-four-base homopolymer rule and are BLAST checked for the
  minimum-span mode. Each selected oligo records its
  quality tier. Repeat masking is off for the R parity baseline; dustmasker
  is not equivalent to RepeatMasker.
- Local BLAST uses `blastn-short`, word size 7, both strands, E-value 1,
  masking off, and the mm10 database. A selected probe has an exact
  full-length hit at its intended locus and no other alignment covering at
  least 80% of the oligo at at least 90% identity. The hit audit is retained,
  including weaker hits returned by BLAST.
- The output selects 30 probes when possible. Any shortfall
  is visible in `set_summary.csv`; unverified candidates never pad a set.
  Order sheets are marked `REVIEW` because BLAST and sequence filters do not
  establish experimental hybridization performance.
