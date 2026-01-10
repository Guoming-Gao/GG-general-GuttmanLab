# Alignment Report

This report summarizes the results of **Step 3: Alignment**.

## Alignment Strategy
- **Reference**: `ref_seq/xist_amplicon_ref.fa` (584 bp)
- **Tool**: `minimap2` (v2.30) with `-ax map-ont` settings for Nanopore reads.
- **Processing**: Alignments were converted to BAM, sorted, and indexed using `samtools`.

## Mapping Statistics
We achieved high-confidence mapping across all four samples.

| Sample / Condition | Primary Alignments | Mapping Rate (%) | Mean MQ | Mean Coverage |
| :--- | :--- | :--- | :--- | :--- |
| **WT (diff)** | 1081 | 85.9% | 59.8 | 1035x |
| **WT (Dox 72h)** | 1172 | 88.3% | 59.9 | 1104x |
| **dTsix (Dox 72h)** | 837 | 81.1% | 59.7 | 784x |
| **dTsix dSPEN (Dox 72h)** | 144 | 84.2% | 60.0 | 135x |

### Observations
- **Mapping Quality (MQ)**: The mean MQ is consistently ~60, indicating that almost all primary alignments are unique and highly confident.
- **Coverage**: Even the sample with the fewest reads (Sample D) has >100x coverage, which is more than sufficient for allele quantification.
- **Supplementary Alignments**: We observed a small percentage (~5-8%) of supplementary alignments, typical for long-read data representing split mappings or technical artifacts.

## Conclusion
The alignment is successful and provides a robust foundation for identifying SNP-bearing reads and calculating allele ratios.

**I am ready to proceed to Step 4: SNP Quantification & Allele Assignment once this report is approved.**
