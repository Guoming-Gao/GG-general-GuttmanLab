# Walkthrough of Current Pipeline Logic (Legacy)

This report details the logic and parameters of the existing `Dox_enabled_amplicon_nanoporeseq` project before the planned updates.

## 1. Reference Initialization (`scripts/initialize_reference.py`)
*   **Purpose**: To create a localized reference sequence for the targeted amplicon.
*   **Logic**:
    *   Loads primers from `ValidatedPrimers.fa`.
    *   Searches the *Xist* locus (`chrX:103450373-103493233`) in `mm10.fa` for exact matches or RCs of the primers.
    *   Defines the amplicon range as the sequence between the two mapped primers.
    *   Fetches B6/Cast SNPs from the VCF within this range.
*   **Parameters**:
    *   Primers: `AC_XistExAmp_5SNPs-F` and `AC_XistExAmp_5SNPs-R`.
    *   Genomic Samples: `C57BL_6NJ` (B6) and `CAST_EiJ` (Cast).
*   **Outputs**: `target_amplicon.fa` and `snps.json` (containing local coordinates).

## 2. QC & Alignment (`scripts/fastq_qc.py` & `scripts/align_reads.py`)
*   **QC Logic**: Calculates basic N50/Length stats from FASTQ files using `pysam`.
*   **Alignment Logic**: Uses `minimap2` with the `-ax sr` (short read) preset to align Nanopore reads to the `target_amplicon.fa`.
    *   **Note**: Using `-ax sr` for Nanopore amplicons is sometimes done for high-accuracy reads, but `-ax map-ont` is more standard for default Nanopore data.
*   **Outputs**: `aligned_reads.sam` -> `aligned_reads.bam`.

## 3. Allele Quantification (`scripts/quantify_alleles.py`)
*   **Purpose**: Assign reads to B6 or Cast alleles.
*   **Logic**:
    *   Iterates through each aligned read in the BAM.
    *   For each SNP position defined in `snps.json`, it queries the BAM for the base at that position (`pysam.AlignedSegment.query_sequence`).
    *   Assigns a "vote" to B6 or Cast based on the base match.
    *   Final assignment is based on majority rule.
*   **Outputs**: `allele_counts.csv` and `per_read_assignments.csv`.

## 4. Stoichiometry & Reliability (`scripts/analyze_stoichiometry.py`)
*   **Purpose**: Assess the internal consistency of SNP calls within single reads.
*   **Logic**: Calculates the percentage of SNP positions within a read that agree with the final allele assignment. If a read matches B6 at 3 SNPs but Cast at 2, it indicates a potential chimeric read or sequencing error.

## 5. Summary & Consolidation (`scripts/generate_reports.py`)
*   **Purpose**: Merges tables and generates a Markdown summary of the results.

---

## Identified Limitations for Current Upgrade
1.  **Reference Dependency**: The current system requires a pre-defined amplicon. The new plan moves to genomic alignment, which is more robust to variable read lengths.
2.  **Trimming**: There is no explicit primer/barcode trimming in the current pipeline; it relies on `minimap2` soft-clipping. The new plan introduces rigorous `cutadapt` trimming.
3.  **Strandness**: The current logic assumes reads are oriented correctly by the reference. The new plan will explicitly validate strandness based on the Forward Primer vs. 2PBC orientation.
