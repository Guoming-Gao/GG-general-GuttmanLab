# Research & Package Selection (High-Speed Mode)

Based on the requirements for error-prone Nanopore amplicon sequencing, strand-specific random priming, and allele quantification, the following peer-reviewed tools and packages have been selected:

| Tool/Package | Purpose | Justification |
| :--- | :--- | :--- |
| **Cutadapt** | Primer & Barcode Trimming | **Standard & Robust**: Peer-reviewed (Martin, 2011). It is the most versatile tool for trimming linked adapters (Forward Primer ... 2PBC Barcode) and supports high error rates (defaulting to 10-20% for Nanopore) and detailed statistical reporting. |
| **Minimap2** | Genomic Alignment | **Gold Standard**: Peer-reviewed (Li, 2018). Optimized specifically for long, error-prone reads. It handles genomic alignment (mm10) significantly better than local amplicon alignment for varying read lengths. |
| **Pysam** | SAM/BAM/VCF Parsing | **Industry Standard**: Essential for programmatically interacting with alignment files and performing high-precision single-read allelic Quantification. |
| **Pandas / Seaborn** | Stats & Visualization | **Data Integrity**: Used for generating the requested statistics on head/tail recovery and plotting length distributions and allele counts. |
| **Nanoplot / NanoComp** | General QC | **Visualization**: Excellent peer-reviewed tools for visualizing Nanopore read quality and length distributions before and after trimming. |

## Strategy for Barcode/Primer Logic
*   **Fuzzy Matching**: We will use `cutadapt`'s `--overlap` and `--error-rate` parameters to handle the inherent errors in Nanopore sequencing at the head (Forward Primer) and tail (RC of 2PBC-9mer).
*   **Strandness**: By defining the head and tail, we enforce a "linked-adapter" check. Only reads containing both (or prioritized parts) will be analyzed, ensuring defined strandness relative to the gene direction.
