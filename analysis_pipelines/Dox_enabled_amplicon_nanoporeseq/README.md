# Nanopore Amplicon Sequencing Pipeline

Whole-genome alignment pipeline for Nanopore amplicon sequencing with SNP-based allelic quantification.

## Requirements

- **Conda environment**: `bioinfo`
- **Key tools**: `minimap2`, `samtools`, `pysam`, `pandas`, `matplotlib`, `numpy`
- **Reference genome**: mm10 FASTA
- **VCF file**: MGP v5 merged SNPs (for B6/Cast strain discrimination)

## Installation

```bash
# Activate conda environment
conda activate bioinfo

# Verify required tools
minimap2 --version
samtools --version
python -c "import pysam, pandas, matplotlib; print('All packages installed')"
```

## Workflow: Mode 1 (No UMI)

### Step 1: Quality Control

Calculate read length distributions, N50, and quality scores.

```bash
python scripts/01_fastq_quality_metrics.py \
  --data_dir /path/to/fastq \
  --output_dir results/qc
```

**Outputs**: `qc_summary.csv`, length distribution plots

---

### Step 2: Genome Alignment

Align reads to mm10 genome using minimap2 with ONT preset.

```bash
python scripts/02_align_to_genome.py \
  --data_dir /path/to/fastq \
  --results_dir results/mode_1 \
  --genome /Volumes/guttman/genomes/mm10/fasta/mm10.fa \
  --threads 8
```

**Outputs**:
- `*.sorted.bam` - Genome-aligned reads (sorted and indexed)
- `*.xist.primary.bam` - Xist locus reads only (chrX:103,460,373-103,483,233)

---

### Step 3: Extract SNPs

Extract B6/Cast discriminating SNPs from VCF for the Xist region.

```bash
python scripts/03_extract_snps_from_vcf.py \
  --vcf /Volumes/guttman/genomes/mm10/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz \
  --region chrX:103460373-103483233 \
  --b6 C57BL_6NJ \
  --cast CAST_EiJ \
  --output results/mode_1/xist_snps.txt
```

**Outputs**: `xist_snps.txt` - TSV file with 101 B6/Cast SNP positions

---

### Step 4: Quantify Alleles

Assign each read to B6 or Cast allele based on SNP match scoring.

```bash
python scripts/04_quantify_alleles.py \
  --bam_dir results/mode_1/aligned \
  --snp_file results/mode_1/xist_snps.txt \
  --output_dir results/mode_1/quantification
```

**Outputs**:
- `*_quant.csv` - Per-read allelic assignments with SNP calls
- `allele_quantification_summary.csv` - Sample-level summary (B6/Cast counts and ratios)

---

### Step 5: Visualize Stoichiometry

Generate 2D heatmaps showing B6 vs Cast SNP match distributions.

```bash
python scripts/05_plot_stoichiometry_heatmap.py \
  --quant_dir results/mode_1/quantification \
  --output_dir results/mode_1/quantification/plots
```

**Outputs**: `*.png` - Heatmaps for each sample

---

### Step 6: Compare Multiple Datasets (Optional)

Compare allelic ratios across multiple experimental datasets.

```bash
python scripts/06_compare_multiple_datasets.py \
  --results_dirs "results/exon_rep1,results/exon_rep2,results/intron" \
  --labels "ExonRep1,ExonRep2,Intron" \
  --output_dir results/comparative_analysis
```

**Outputs**:
- `comparison_summary.csv` - Cast ratios across conditions
- `Comparative_Analysis_Report.md` - Markdown report with plots

---

### Step 7: Generate Summary Report (Optional)

Auto-generate a consolidated markdown report.

```bash
python scripts/07_generate_summary_report.py \
  --results_dir results/mode_1 \
  --omit_snps "2,5"  # Optional: comma-separated SNP indices to exclude
```

**Outputs**: `Automated_Summary_Report.md`

---

## UMI-Based Analysis (Experimental)

For UMI-containing libraries, use `trim_umi_adapters.py` to extract 6bp UMIs and trim adapters before alignment.

```bash
python scripts/trim_umi_adapters.py \
  --fastq /path/to/sample.fastq \
  --output_dir results/mode_2/categorized \
  --f_primer AGCAGACCACACAGGCACCAGA \
  --r_primer GCAGTCTGGTGTGAGGAACGGC
```

**Note**: Current UMI libraries show <2% on-target rate. This workflow is included for reference but **not recommended for production use** until library prep issues are resolved.

---

## Output File Structure

```
results/
├── qc/
│   ├── qc_summary.csv
│   └── *_length_dist.png
├── mode_1/
│   ├── aligned/
│   │   ├── *.sorted.bam
│   │   ├── *.sorted.bam.bai
│   │   └── *.xist.primary.bam
│   ├── quantification/
│   │   ├── *_quant.csv
│   │   ├── allele_quantification_summary.csv
│   │   └── plots/*.png
│   └── xist_snps.txt
└── mode_2/ (experimental)
    ├── categorized/
    │   └── *.trimmed.fastq
    └── reports/
```

---

## Key Features

- **Whole-genome alignment**: Avoids local amplicon reference issues
- **Robust SNP calling**: 101 B6/Cast SNPs for high-confidence allelic assignment
- **Modular design**: Each script is standalone and reusable
- **Comprehensive QC**: Length distributions, mapping rates, stoichiometry validation
- **Efficient processing**: Parallel alignment (~5 min/sample with 8 threads)

---

## Validated Results (Mode 1)

- **8 samples processed**: A-H (WT, dTsix, dSPEN, dRex1, TsixOE variants)
- **High on-target rate**: 80-95% of reads map to Xist locus
- **Robust allelic purity**: Clear B6/Cast separation with minimal ambiguous reads
- **Validation**: SNP match heatmaps confirm expected allelic stoichiometry

---

## Citation

[Add publication info when available]

---

## Contact

Guoming Gao (Guttman Lab)
