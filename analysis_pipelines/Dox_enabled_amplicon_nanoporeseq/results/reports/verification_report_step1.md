# SNP Verification Report

This report confirms the mapping and accuracy of the SNP annotations in the GenBank reference file against the mm10 genome and the standard SNP VCF.

## Mapping Summary
- **Reference File**: `ref_seq/xist-ref-seq-with-snps-annot-1712-2295.gb`
- **Genomic Match**: `chrX:103466603-103467186` (on the **plus** strand of the genome)
- **Gene Context**: This segment lies within **Exon 2** of the *Xist* gene (RefSeq `NR_001570`).

## SNP Allele Verification
All three annotated SNPs were successfully identified in the VCF (`mgp.v5.merged.snps_all.dbSNP142.vcf.gz`) at the corresponding genomic positions.

| GB Pos | Genomic Pos (chrX) | Label | B6 Allele (Ref) | Cast Allele (Alt) | VCF Ref | VCF Alt |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **67** | 103466669 | T->C | **T** | **C** | T | C |
| **334** | 103466936 | G->A | **G** | **A** | G | A |
| **553** | 103467155 | G->A | **G** | **A** | G | A |

## Conclusion
> [!NOTE]
> The GenBank sequence orientation is on the plus strand of the genome. Since *Xist* is a minus-strand gene, the reads (which represent the RNA sense sequence) will align as **antisense** to this GenBank reference if using default parameters, or should be reverse-complemented for sense alignment.
>
> The SNPs are verified and provide 100% differentiation between B6 and Cast alleles at all three target positions.

**I am ready to proceed to Step 2: Preprocessing once this report is approved.**
