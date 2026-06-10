# smfish-like-rt-probe-designer

Design focused RT primers for SWIFT-seq using smFISH probe design principles (Oligostan), with SNP coverage analysis for B6xCast allelic expression studies.

## Quick Start

```bash
conda activate bioinfo
python main.py --genes Nanog Pou5f1 Xist --output ./my_output
```

## Pipeline Overview

1. **Phase 1: Candidate Generation** (`1_design_candidate_probes.py`)
   - Fetch sequences from GTF + FASTA
   - Design primers using Oligostan (dG37 = -32.0)
   - Filter by GC, PNAS rules, and Dustmasker
   - Analyze SNPs in RT coverage region
   - Output: `FISH_RT_probes_CANDIDATES.csv/.fasta`

2. **Phase 2: Specificity Validation** (`2_validate_probe_specificity.py`)
   - Automated `blastn` specificity check
   - HSP-aware hit counting logic
   - Merges BLAST results with candidate data
   - Output: `FISH_RT_probes_FINAL_SELECTION.csv`

---

## RT Coverage and SNP Detection Logic

The key to understanding SNP coverage is the directionality of RT extension:

### For **+ strand** genes (e.g., Nanog on chr6):

```
                        gene direction →
                        mRNA 5'→3'
                        ─────────────────────────────────────────▶
Genomic DNA:  5'═══════════════════════════════════════════════════3'
                        ▲           ▲
                        │           │
              RT region │  probe    │
              (200 bp)  │  binds    │
              ◀─────────┤  here     │
                        │           │
RT extends this way ◀───┘           │
(3'→5' on mRNA)                     │
                                    ▼
                        ┌───────────┐
Probe (DNA):            │ 5'────3'  │  (complementary to mRNA)
                        └───────────┘
                              │
                              ▼ RT extends
                        ┌─────────────────────────┐
RT Product (cDNA):      │ 5'────────────────────3'│  (extends ~200bp UPSTREAM)
                        └─────────────────────────┘
                        ◀─────────────────────────▶
                        RT_Region_Start    Probe_End
                        (lower coords)     (higher coords)
```

**+ strand summary:**
- Probe binds at `Probe_Start..Probe_End` (higher genomic coordinates)
- RT extends **UPSTREAM** (toward lower genomic coordinates)
- `RT_Region = [Probe_Start - 200, Probe_Start - 1]`
- SNPs checked in RT region for allelic distinction

---

### For **- strand** genes (e.g., Xist on chrX):

```
                        ◀ gene direction
                          mRNA 5'→3'
                        ◀─────────────────────────────────────────
Genomic DNA:  5'═══════════════════════════════════════════════════3'
                                    ▲           ▲
                                    │           │
                                 probe      RT region
                                 binds      (200 bp)
                                 here       ─────────▶
                                    │
                                    └───▶ RT extends this way
                                          (3'→5' on mRNA)

                        ┌───────────┐
Probe (DNA):            │ 5'────3'  │  (complementary to mRNA)
                        └───────────┘
                              │
                              ▼ RT extends
                        ┌─────────────────────────┐
RT Product (cDNA):      │ 5'────────────────────3'│  (extends ~200bp DOWNSTREAM)
                        └─────────────────────────┘
                        ◀─────────────────────────▶
                        Probe_Start       RT_Region_End
                        (lower coords)    (higher coords)
```

**- strand summary:**
- Probe binds at `Probe_Start..Probe_End` (lower genomic coordinates)
- RT extends **DOWNSTREAM** (toward higher genomic coordinates)
- `RT_Region = [Probe_End + 1, Probe_End + 200]`
- SNPs checked in RT region for allelic distinction

---

### Key Insight

**RT always extends 3'→5' on the mRNA template**, which means:
- `+ strand gene`: RT goes toward **lower** genomic coordinates (upstream)
- `- strand gene`: RT goes toward **higher** genomic coordinates (downstream)

This is critical for correct SNP coverage calculation!

---

## Design Principles and Internal Logic

To ensure accurate probe design and region classification (intron/exon), the pipeline follows these steps:

### 1. Transcript Selection
By default, the pipeline selects the **"longest spliced" transcript** for each gene (sum of all exon lengths). This ensures that probes are designed against the most comprehensive representation of the gene's mature product.
- **Customizable**: Users can override this by providing a `transcript_id_override` map in `FISH_RT_CONFIG`.
- **Ambiguity Handling**: If multiple transcripts share the same maximum length, the first one encountered in the GTF is used.

### 2. Sequence Extraction (Strand-Aware)
All sequences (gene genomic sequence, RT product sequences, and probe sequences) are stored in the **RNA sense (5'→3')** direction.
- For **negative strand genes**, this means the genomic region is automatically reverse-complemented before any design steps occur.
- Probes are designed to be **complementary** to this RNA sense sequence.

### 3. Coordinate System (1-indexed inclusive)
The pipeline uses a unified **1-indexed inclusive** system for all relative coordinates:
- The first base of the extracted gene sequence is **Position 1**.
- Exon and intron regions are defined as inclusive ranges `[start, end]`.
- A probe binding at `[70, 100]` is 31nt long and covers those exact positions.

### 4. Region Classification (Intron vs Exon)
A probe's `RegionType` is determined by its overlap with the selected transcript's exons and introns:
- **exon**: All bases of the probe overlap with one or more exons.
- **intron**: All bases of the probe overlap with one or more introns.
- **exon-intron**: The probe spans a junction between an exon and an intron.
- **Robustness**: As a safety measure, a probe is labeled **exon** if it overlaps with an exon in *any* known transcript for that gene, preventing accidental targeting of splicing variants.

---

## Output Files

| Step | File | Description |
|------|------|-------------|
| 1 | `FISH_RT_probes_CANDIDATES.csv` | All high-quality candidates before BLAST |
| 2 | `FISH_RT_probes_FINAL_SELECTION.csv` | **Main Output**: Final probe set after BLAST filtering |
| 2 | `FISH_RT_probes_FINAL_SELECTION.fasta` | FASTA format for synthesis |
| 2 | `FISH_RT_probes_FINAL_SELECTION_summary.txt` | Design statistics and quality summary |

## Workflow

```bash
# Phase 1: Generate candidates (Design + SNPs)
conda activate bioinfo
python 1_design_candidate_probes.py --genes Nanog Xist --output ./results

# Phase 2: Run Specificity Validation (BLAST + Filtering)
python 2_validate_probe_specificity.py --candidates ./results/FISH_RT_probes_CANDIDATES.csv --output-dir ./results
```

## Configuration

Edit `config.py`:

```python
FISH_RT_CONFIG = {
    # Input files (mm10)
    "local_gtf_path": "/Volumes/guttman/genomes/mm10/annotation/mm10.refGene.gtf.gz",
    "local_genome_fasta_path": "/Volumes/guttman/genomes/mm10/fasta/mm10.fa",
    "snp_file_path": "/Volumes/guttman/genomes/mm10/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz",

    # VCF samples for allelic SNP detection
    "vcf_b6_sample": "C57BL_6NJ",
    "vcf_cast_sample": "CAST_EiJ",

    # Parameters
    "rt_coverage_downstream": 100,      # nt downstream for SNP detection
    "max_probes_per_gene": 200,         # Maximum probes per gene
    "min_snps_for_selection": 3,        # Min SNPs to include probe in selection
    "fixed_dg37_value": -32.0,          # Target Gibbs free energy
    "probe_length_min": 26,
    "probe_length_max": 32,
}
```

## Requirements

- Python 3.7+
- **pysam** (for VCF/tabix queries)
- biopython, pandas, numpy, rich
- samtools (for FASTA indexing)

**One-time setup** (if FASTA not indexed):
```bash
cd /Volumes/guttman/genomes/mm10/fasta/
gunzip -k mm10.fa.gz && samtools faidx mm10.fa
```

## File Formats Supported

| File | Format | Notes |
|------|--------|-------|
| GTF | Gzipped refGene or Ensembl | Auto-infers genes from transcripts if no "gene" features |
| FASTA | Indexed (`.fai` required) | Use `samtools faidx` to create index |
| SNP | VCF with tabix index (`.tbi`) | Filters for B6 vs Cast **different** genotypes only |

## Test Results (21 genes, mm10)

| Metric | Value |
|--------|-------|
| Total probes designed | 41,118 |
| After quality filtering | 12,603 (30.7%) |
| Probes with ≥5 SNPs | 129 |
| SNP analysis time | ~2 min (tabix) |

### Validation

Probes at identical genomic positions produce identical results:
- Same sequences, dG37 scores, and filter outcomes
- Example: Mecp2 @ chr74027166 → `ATCTTTCTCTCTGCCTTGTCTGCCTGCTCCC`, dG37=-42.978

### SNP Detection Comparison (vs old custom .snps format)

| SNP Count | Old Format | New VCF |
|-----------|------------|---------|
| ≥5 SNPs | 78 probes | **129 probes** (+65%) |

The VCF-based approach finds more high-SNP probes by properly filtering for B6 vs Cast genotype differences.

## Key Files

| File | Purpose |
|------|---------|
| `1_design_candidate_probes.py` | Phase 1: Candidate Design & SNPs |
| `2_validate_probe_specificity.py` | Phase 2: BLAST & Hit Filtering |
| `config.py` | All configuration parameters |
| `gene_fetcher.py` | GTF/FASTA parsing (gzip support) |
| `snp_analyzer.py` | VCF-based SNP analysis with pysam/tabix |
| `4_optional_design_forward_primers.py` | Automated Primer3-based primer design |
| `3_add_rtbc_barcode.py` | Post-processing with synthesis barcodes |

## Citation

Thermodynamic design based on:
> Tsanov et al. "smiFISH and FISH-Quant" *Nucleic Acids Research* 44(22):e165, 2016. https://doi.org/10.1093/NAR/GKW784

## License

MIT License
