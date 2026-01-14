
# Snapshot Report: De Novo Reference Generation & 5-SNP Analysis

This report documents the fix implemented to transition the pipeline from fixed GenBank annotations to dynamic, de novo reference generation based on user-provided primers.

## 1. De Novo Initialization
We extracted the reference sequence and SNP list directly from the lab's shared resources using the primer pair:
- **Forward**: `AC_XistExAmp_5SNPs-F` (TCTGGTGCCTGTGTGGTCTGCT)
- **Reverse**: `AC_XistExAmp_5SNPs-R` (GCAGTCTGGTGTGAGGAACGGC)

**Genomic Locus**: `chrX:103466603-103467186` (mm10, Plus Strand)

## 2. Updated SNP Profile
Scanning the shared VCF (`mgp.v5.merged.snps_all.dbSNP142.vcf.gz`) revealed **5 B6/Cast SNPs** within the amplicon, correcting the previous count of 3.

| Local Pos (BP) | Genomic Pos (chrX) | B6 Allele | Cast Allele |
| :--- | :--- | :--- | :--- |
| 28 | 103466630 | **G** | **A** |
| 67 | 103466669 | **T** | **C** |
| 71 | 103466673 | **C** | **G** |
| 334 | 103466936 | **G** | **A** |
| 553 | 103467155 | **G** | **A** |

## 3. Revised Allele Quantification
Assignment now follows a majority-rule threshold of **3/5 matching SNPs** for high-confidence calls.

| Sample / Condition | Total Reads | B6 Reads | Cast Reads | Cast Ratio |
| :--- | :---: | :---: | :---: | :---: |
| **WT (diff)** | 1081 | 509 | 557 | **0.522** |
| **WT (Dox 72h)** | 1172 | 1121 | 38 | **0.033** |
| **dTsix (Dox 72h)** | 837 | 719 | 101 | **0.123** |
| **dTsix dSPEN (Dox 72h)** | 144 | 81 | 59 | **0.421** |

## 4. Stoichiometry Validation (5-SNP)
Across the major B6 populations, **~57-60%** of all reads carry **all 5 SNPs** simultaneously, confirming the high fidelity of the de novo reference and the robustness of the allelic assignments.

---
**Status**: The pipeline is now fully de novo and relies exclusively on the mm10 genome, shared VCF, and the target primer pair.

**I am ready for the next "suspicious thing" you've identified.**
