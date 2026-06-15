"""Configuration for KOveri primer design."""

from pathlib import Path


KOVERI_CONFIG = {
    # Local mm10 / GRCm38 resources, matching smfish-like-rt-probe-designer.
    "genome_fasta": "/Volumes/guttman/genomes/mm10/fasta/mm10.fa",
    "blast_database": "/Volumes/guttman/genomes/mm10/fasta/blastdb/mm10_blastdb",
    "snp_vcf": "/Volumes/guttman/genomes/mm10/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz",
    "vcf_b6_sample": "C57BL_6NJ",
    "vcf_cast_sample": "CAST_EiJ",
    # bioinfo conda environment tools.
    "samtools": "/opt/miniconda3/envs/bioinfo/bin/samtools",
    "blastn": "/opt/miniconda3/envs/bioinfo/bin/blastn",
    "primer3_core": "/opt/miniconda3/envs/bioinfo/bin/primer3_core",
    "ntthal": "/opt/miniconda3/envs/bioinfo/bin/ntthal",
    # Design behavior.
    "min_informative_snps": 3,
    "flank_steps": [250, 500, 750, 1000, 1500, 2500, 5000],
    "min_amplicon_size": 150,
    "ideal_amplicon_max": 700,
    "max_amplicon_size": 5000,
    "max_cluster_span": 100000,
    "max_target_spans_per_flank": 30,
    "primer3_num_return": 250,
    "top_n_output": 20,
    # Post-Primer3 filter settings.
    "valid_primer_length_min": 18,
    "valid_primer_length_max": 30,
    "post_filter_gc_min": 40.0,
    "post_filter_gc_max": 60.0,
    "relaxed_post_filter_gc_min": 35.0,
    "relaxed_post_filter_gc_max": 65.0,
    "post_filter_tm_min": 55.0,
    "post_filter_tm_max": 65.0,
    "relaxed_post_filter_tm_min": 52.0,
    "relaxed_post_filter_tm_max": 68.0,
    "strict_tm_delta_max": 5.0,
    "strict_amplicon_min": 400,
    "strict_amplicon_max": 700,
    "max_homopolymer_run": 3,
    "max_quad_g_run": 3,
    "max_alternating_dinuc_bases": 5,
    "max_interprimer_complement_bases": 3,
    "min_thermo_dg_kcal": -9.0,
    # Primer3 settings.
    "primer_min_size": 18,
    "primer_opt_size": 22,
    "primer_max_size": 30,
    "primer_min_tm": 55.0,
    "primer_opt_tm": 60.0,
    "primer_max_tm": 65.0,
    "primer_min_gc": 40.0,
    "primer_max_gc": 60.0,
    "primer_gc_clamp": 1,
    "primer_max_poly_x": 4,
    "primer_max_self_any": 8.0,
    "primer_max_self_end": 3.0,
    "primer_pair_max_compl_any": 8.0,
    "primer_pair_max_compl_end": 3.0,
    # BLAST specificity scoring for individual primers.
    "blast_specificity": True,
    "specificity_candidates": 60,
    "blast_min_unique_len": 18,
    "blast_min_unique_identity": 95.0,
    "blast_noise_len": 12,
}


DEFAULT_OUTPUT_DIR = Path(
    "/Users/gmgao/Dropbox/Caltech_PostDoc_GuttmanLab/constructs_primers_FISHprobes/"
    "qPCR LibAmp PCR primers/AmpliconSeq-KOveri-Tsix"
)
DEFAULT_INPUT_FASTA = DEFAULT_OUTPUT_DIR / "Tsix_gRNAs.fa"
