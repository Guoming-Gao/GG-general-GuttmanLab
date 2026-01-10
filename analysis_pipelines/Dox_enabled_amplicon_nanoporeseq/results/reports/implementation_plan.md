# Xist SNP Analysis Pipeline Implementation Plan

This project implements a pipeline to analyze Nanopore sequencing reads of an Xist amplicon. The goal is to quantify allele ratios (Cast vs B16/129) and investigate the stoichiometry of SNPs per read, distinguishing biological signal from sequencing noise.

## Proposed Changes

### 1. Environment & Tools
We will use the existing `bioinfo` environment. If new packages are required (e.g., `pycoQC`, `pysam`, `biopython`), we will install them using `mamba` to ensure fast and reliable dependency resolution.

Selected tools:
- **minimap2**: For fast and accurate alignment of long Nanopore reads.
- **pysam**: To interact with alignment data (SAM/BAM) and perform pileup/allele-counting.
- **BioPython**: To parse the annotated GenBank reference sequence and extract SNP metadata.
- **Pandas & Matplotlib**: For data aggregation, ratio calculation, and visualization.

### 2. Pipeline Workflow

#### Step 1: SNP Annotation Verification (CRITICAL)
Before any alignment, we must verify the SNP annotations in the provided GenBank file (`ref_seq/xist-ref-seq-with-snps-annot-1712-2295.gb`):
- Cross-check the annotated positions (67, 334, 553) against the mm10 reference genome and the SNP VCF (`/Volumes/guttman/genomes/mm10/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz`).
- Confirm the specific alleles for `CAST_EiJ` and `C57BL_6NJ` (B16/129) at these positions.
- Ensure the GenBank sequence segment correctly matches the genomic coordinates in mm10.
- **Wait for USER approval of the Verification Report before proceeding.**

#### Step 2: Preprocessing
- Convert verified GenBank reference to FASTA for alignment.
- Basic QC of input FASTQ files.
- **Wait for USER approval of the Preprocessing Report before proceeding.**

#### Step 3: Alignment
- Align reads to the amplicon reference using `minimap2`.
- Sort and index the resulting BAM files using `samtools`.
- **Wait for USER approval of the Alignment Report before proceeding.**

#### Step 4: SNP Quantification & Allele Assignment
- Parse SNP annotations from the verified reference.
- Assign reads to Cast or B16 alleles based on SNP matches.
- **Wait for USER approval of the Quantification Report before proceeding.**

#### Step 5: SNP Stoichiometry & Co-occurrence Analysis
Instead of general noise analysis, we will focus on the biological variation within allelic reads:
- For both Cast and B6 reads, calculate the distribution of "SNP hits" (how many reads have 1, 2, or 3 matching SNPs).
- Calculate co-occurrence frequencies between pairs and triplets of SNPs (e.g., SNP1 & SNP2, SNP1 & SNP3, etc.).
- Export detailed co-occurrence tables per sample.
- **Wait for USER approval of the Final Analysis Report.**

### 3. Data Organization & Export
To ensure re-analysis and plotting capabilities later, we will:
- Consolidate all intermediate CSVs (QC, Alignment, Quantification, Stoichiometry) into a `results/` directory.
- Export a "Final Cumulative Report" as a markdown file that synthesizes all verification reports.
- Ensure the `scripts/` directory contains all Python tools used in the pipeline for reproducibility.

### Automated Tests
- **Unit Tests**: Test the SNP extraction logic using a dummy BAM file with known SNP variations.
- **Integration Test**: Run the full pipeline on a subset of reads and verify that the output CSV contains expected columns (ReadID, SNP1_Base, SNP2_Base, SNP3_Base, Allele_Call).

### Manual Verification
- Review the generated ratio plot.
- Inspect a few "noisy" reads in a genome browser (IGV) if possible.
