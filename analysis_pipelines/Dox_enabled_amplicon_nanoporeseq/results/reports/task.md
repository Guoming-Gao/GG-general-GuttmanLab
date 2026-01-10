# Task: Xist SNP Analysis Pipeline

## Phase 1: Research & Package Selection
- [x] Explore local resources and config files <!-- id: 0 -->
- [x] Research standard ONT SNP analysis tools <!-- id: 1 -->
- [x] Select and justify tools/packages <!-- id: 2 -->

## Phase 2: Agentic Planning & Orchestration
- [x] Update Implementation Plan based on feedback <!-- id: 3 -->
- [x] Incorporate SNP Verification step <!-- id: 12 -->
- [x] Define Sub-Agent strategy <!-- id: 4 -->
- [x] Define Data Validation strategy <!-- id: 5 -->

## Phase 3: Implementation (Pending Approval)
- [x] Setup virtual environment (`bioinfo`) and install dependencies <!-- id: 6 -->
- [x] Verify SNP annotations with mm10 and VCF <!-- id: 13 -->
    - [x] Map Xist coordinates in mm10 <!-- id: 14 -->
    - [x] Write and run `verify_snps.py` <!-- id: 15 -->
    - [x] Generate SNP Verification Report <!-- id: 16 -->
- [/] Implement data preprocessing (QC/Filtering) <!-- id: 7 -->
    - [x] Convert GenBank to FASTA <!-- id: 17 -->
    - [x] Perform FASTQ QC and stats <!-- id: 18 -->
    - [x] Generate Preprocessing Report <!-- id: 19 -->
- [x] Implement alignment to reference <!-- id: 8 -->
    - [x] Align reads with minimap2 <!-- id: 20 -->
    - [x] Sort and index BAM files <!-- id: 21 -->
    - [x] Generate Alignment Report <!-- id: 22 -->
- [x] Implement SNP quantification and allele assignment <!-- id: 9 -->
    - [x] Extract SNP bases per read <!-- id: 23 -->
    - [x] Assign alleles and calculate ratios <!-- id: 24 -->
    - [x] Generate Quantification Report <!-- id: 25 -->
- [x] Implement stoichiometry and co-occurrence analysis <!-- id: 10 -->
    - [x] Calculate SNP hit distributions (1, 2, 3 SNPs) <!-- id: 26 -->
    - [x] Calculate co-occurrence frequencies <!-- id: 27 -->
- [x] Final validation, data organization, and reporting <!-- id: 11 -->
    - [x] Consolidate results into `results/` folder <!-- id: 28 -->
    - [x] Generate Final Cumulative Report <!-- id: 29 -->
