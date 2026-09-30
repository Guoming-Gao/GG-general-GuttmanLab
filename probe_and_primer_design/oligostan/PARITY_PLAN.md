# Oligostan parity and SPEN validation plan

## Reference and scope

The reference is Tsanov et al.'s `Oligostan.r` in FISH-quant:
https://bitbucket.org/muellerflorian/fish_quant/src/master/Oligostan/Oligostan.r
The local reference copy has SHA-256
`3f4a75da1403a77617194de327d72fa67b637102373bbdf9da3cca26f2ca253e`.
Its settings are customized (GC 0.2–0.6 and no PNAS rules), so the SPEN
comparison uses a temporary copy with matched settings: dG37 −32, GC 0.4–0.6,
PNAS rules 1/2/4, score ≥0.9, 26–32 nt, 2 nt spacing and masking off.
The Python source was copied from the standalone `oligostan-python` repository
(commit `dc0f579`) into this directory. The RT-probe designer has a modified
Oligostan-like algorithm and is **not** an independent R reference.

## Parity gate, in order

1. Pin one mouse genome build and input: mm10/GRCm38, the exact FASTA bases,
   strand and annotation transcript for `Spen`. Record SHA-256 of input and
   reference code. Disable repeat masking in both implementations first, then
   test masking as a separate, explicitly non-equivalent method (R
   RepeatMasker versus Python dustmasker).
2. Run the same SPEN FASTA and parameters through R and Python. Compare the
   *ordered* candidate list, positions, lengths, sequences, scores, dG37,
   GC and each PNAS flag before comparing the final filtered set and FLAPs.
   Float tolerance: 1e-9; all sequence/coordinate/filter values exact.
3. Isolate discrepancies by stage: reverse complement, 16 dinucleotide
   values/salt correction, length and score tie handling, greedy spacing,
   one-based coordinates, filter boundaries and output ordering. Add small
   regression fixtures for each corrected case, plus full SPEN R and Python
   TSV outputs and a machine-readable difference report.
4. Benchmark wall time and peak RSS on SPEN with both tools after output
   parity is achieved. Repeat three times and report median, environment and
   versions. "Equal performance" means equivalent design results; speed and
   memory are reported separately, not assumed equal.
5. Repeat on a short known fixture and opposite-strand and boundary cases.
   Only then label the port R-equivalent. Recheck whenever parameters, genome
   build or masking strategy changes.

## Current status

The matched-settings full mouse Spen run passed the R comparison: 1,821
candidate rows, all sequences/positions/lengths, scores, dG37, GC and PNAS
flags, and 725 GC+PNAS-filtered probes matched. The comparison ignored only
file naming and R-only masking fields. It used the temporary R source with
the explicit settings above and the same `Spen.fa` bases. The Python port's
fifth composition flag was corrected to mirror the original R script's
first-six-window behavior. The supplied human RNU1 fixture also passed.

This establishes parity for the **fixed −32, masking-off settings used by this
notebook**. The original R script's automatic dG37 optimization, RepeatMasker
behavior and runtime/peak-memory benchmark remain separate gates and are not
claimed equivalent. The user requested that the two Python package copies be
kept identical instead of repeating a Python-vs-Python comparison after sync.
