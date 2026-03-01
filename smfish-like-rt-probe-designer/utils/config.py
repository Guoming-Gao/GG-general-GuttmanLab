# config.py - smfish-like-rt-probe-designer Configuration (Local Files + Stringent Filtering)
import os

# =============================================================================
# MAIN CONFIGURATION (LOCAL FILES + STRINGENT FILTERING)
# =============================================================================

FISH_RT_CONFIG = {
    # INPUT SETTINGS
    "transcript_selection": "longest",
    "output_directory": "/Users/gmgao/Dropbox/Caltech_PostDoc_GuttmanLab/constructs_and_smiFISH/smFISH_like_focusedRT-XCI",
    # LOCAL FILE PATHS (YOUR SERVER FILES)
    "local_gtf_path": "/Volumes/guttman/genomes/mm10/annotation/mm10.refGene.gtf.gz",
    "local_genome_fasta_path": "/Volumes/guttman/genomes/mm10/fasta/mm10.fa",  # Requires one-time: gunzip -k mm10.fa.gz && samtools faidx mm10.fa
    "snp_file_path": "/Volumes/guttman/genomes/mm10/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz",
    # VCF SAMPLE NAMES (for B6 x Cast SNP filtering)
    "vcf_b6_sample": "C57BL_6NJ",  # B6 reference strain
    "vcf_cast_sample": "CAST_EiJ",  # Cast strain for allelic analysis
    # RT COVERAGE SETTINGS - RT extends 3'→5' on mRNA (upstream direction)
    "rt_coverage_downstream": 500,  # nt that RT extends (based on experimental data showing >200nt products)
    "include_probe_in_coverage": False,  # Strictly downstream
    "extract_rt_sequence": True,  # Extract RT product sequence for output
    # RTBC BARCODE SETTINGS
    "add_rtbc_barcode": False,  # Changed: RTBC now added as final post-processing step
    "rtbc_sequence": "/5Phos/TGACTTGAGGAT",
    # OLIGOSTAN PARAMETERS
    "fixed_dg37_value": -32.0,
    "score_min": 0.9,
    "probe_length_max": 32,
    "probe_length_min": 26,
    "probe_spacing_min": 2,
    "gc_content_min": 0.4,
    "gc_content_max": 0.6,
    "salt_concentration": 0.115,
    # STRINGENT FILTERING (FIXED - all features restored)
    "pnas_filter_rules": [1, 2, 4],  # 3 rules is often enough
    "use_dustmasker": True,  # ENABLED (FIXED)
    "max_masked_percent": 0.1,
    # SNP COVERAGE FILTERING
    # "min_snp_coverage_for_final": 5,  # DEPRECATED: Threshold for HIGH_SNP output (now redundant)

    "max_homopolymer_length": 4,  # Consecutive bases allowed (e.g. 4 means AAAAA is blocked)
    # PROBE SELECTION SETTINGS
    "max_probes_per_gene": 500,       # Maximum probes per gene (was hardcoded 50)
    "min_snps_for_selection": 3,      # Minimum SNPs to include probe in selection
    "generate_all_probes_file": False,  # FIXED: Skip ALL.csv generation
    "focus_on_high_quality_only": True,  # FIXED: Only process filtered probes
    # FASTA OUTPUT SETTINGS
    "generate_blast_fasta": True,  # Generate FASTA for manual BLAST
    "fasta_description_format": "detailed",  # 'simple' or 'detailed'

    # =========================================================================
    # SYNTHESIS OUTPUT SETTINGS
    # =========================================================================
    # Default number of probes per gene for final synthesis output
    "default_probes_per_gene": 3,

    # Per-gene probe count overrides (use this for custom counts per gene)
    # Example: {"Xist": 20, "Tsix": 10} means 20 probes for Xist, 10 for Tsix, default for others
    "probes_per_gene_override": {},

    # =========================================================================
    # LOCAL BLAST SETTINGS
    # =========================================================================
    "run_local_blast": True,  # Set to False to skip BLAST (faster, but no specificity check)
    "blast_database_path": "/Volumes/guttman/genomes/mm10/fasta/blastdb/mm10_blastdb",  # Path to local BLAST database
    "blast_min_identity": 90.0,  # Minimum % identity for a hit to count
    "blast_max_evalue": 0.001,  # Maximum E-value for a hit
    "filter_unique_blast_hits": True,  # Only keep probes with exactly 1 BLAST hit
}

# =============================================================================
# THERMODYNAMIC PARAMETERS (unchanged)
# =============================================================================

DG37_VALUES = {
    "AA": -0.2,
    "AC": -1.5,
    "AG": -0.9,
    "AT": -1.0,
    "CA": -1.0,
    "CC": -2.2,
    "CG": -1.2,
    "CT": -1.4,
    "GA": -0.8,
    "GC": -2.4,
    "GG": -1.5,
    "GT": -1.0,
    "TA": -0.3,
    "TC": -1.4,
    "TG": -1.0,
    "TT": -0.4,
}

# =============================================================================
# TEST GENE SETS (REMOVED - Use --genes argument)
# =============================================================================

# =============================================================================
# DEMO GENE COORDINATES (mm10) - For when local files aren't accessible
# =============================================================================

DEMO_GENE_COORDS = {
    "Nanog": {
        "chromosome": "6",
        "start": 122707489,
        "end": 122714633,
        "strand": 1,
        "biotype": "protein_coding",
    },
    "Pou5f1": {
        "chromosome": "17",
        "start": 35506018,
        "end": 35510772,
        "strand": 1,
        "biotype": "protein_coding",
    },
    "Sox2": {
        "chromosome": "3",
        "start": 34650840,
        "end": 34652882,
        "strand": -1,
        "biotype": "protein_coding",
    },
}
