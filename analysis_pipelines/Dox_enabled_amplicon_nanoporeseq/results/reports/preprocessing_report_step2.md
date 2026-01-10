# Preprocessing Report

This report summarizes the results of **Step 2: Preprocessing & QC**.

## Reference Conversion
- **Source**: `ref_seq/xist-ref-seq-with-snps-annot-1712-2295.gb` (Verified in Step 1)
- **Output**: `ref_seq/xist_amplicon_ref.fa`
- **Reference ID**: `Xist_Amplicon`
- **Length**: 584 bp

## FASTQ Quality Control
We analyzed 4 experimental conditions. Overall, the data quality is exceptionally high for Nanopore sequencing.

### Summary Table
| Condition / File | Total Reads | Mean Length | Median Length | Mean Quality (Q) |
| :--- | :--- | :--- | :--- | :--- |
| **WT (diff)** | 1258 | 522 bp | 570 bp | 36.6 |
| **WT (Dox 72h)** | 1240 | 539 bp | 575 bp | 36.5 |
| **dTsix (Dox 72h)** | 946 | 516 bp | 573 bp | 36.0 |
| **dTsix dSPEN (Dox 72h)** | 158 | 521 bp | 571 bp | 35.1 |

### Observations
- **Read lengths**: The median read lengths (570-575 bp) are very close to the expected 584 bp amplicon size, suggesting minimal fragmentation.
- **Quality**: Mean Q-scores around 36 indicate >99.9% base-call accuracy, which is ideal for high-confidence SNP calling.
- **Throughput**: Sample D has significantly fewer reads (158) compared to the others, but still sufficient for a preliminary stoichiometry analysis.

## Conclusion
The data is high-quality and ready for alignment.

**I am ready to proceed to Step 3: Alignment once this report is approved.**
