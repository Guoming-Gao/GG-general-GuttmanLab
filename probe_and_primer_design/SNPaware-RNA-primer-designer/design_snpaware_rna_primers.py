#!/usr/bin/env python3
"""Design SNP-aware mature-RNA PCR primers from a gene symbol."""

from __future__ import annotations

import argparse
from pathlib import Path

from rich.console import Console

from utils.designer import design_snpaware_rna_primers, design_snpaware_rna_primers_for_genes


console = Console()

DEFAULT_GENE = "Xist"
DEFAULT_OUTPUT_DIR = Path(
    "/Users/gmgao/Dropbox/Caltech_PostDoc_GuttmanLab/constructs_primers_FISHprobes/"
    "qPCR LibAmp PCR primers/AmpliconSeq-Xist-crossjunction"
)

DEFAULT_CONFIG = {
    "gtf_path": "/Volumes/guttman/genomes/mm10/annotation/mm10.refGene.gtf.gz",
    "genome_fasta": "/Volumes/guttman/genomes/mm10/fasta/mm10.fa",
    "blast_database": "/Volumes/guttman/genomes/mm10/fasta/blastdb/mm10_blastdb",
    "snp_vcf": "/Volumes/guttman/genomes/mm10/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz",
    "vcf_b6_sample": "C57BL_6NJ",
    "vcf_cast_sample": "CAST_EiJ",
    "conda_env": "bioinfo",
    "samtools": "/opt/miniconda3/envs/bioinfo/bin/samtools",
    "blastn": "/opt/miniconda3/envs/bioinfo/bin/blastn",
    "primer3_core": "/opt/miniconda3/envs/bioinfo/bin/primer3_core",
    "ntthal": "/opt/miniconda3/envs/bioinfo/bin/ntthal",
    "curated_transcript_prefixes": ("NM_", "NR_"),
    "require_common_junctions": True,
    "min_junction_overlap_bases": 6,
    "min_junction_3p_anchor_bases": 4,
    "min_informative_snps": 3,
    "flank_steps": [250, 500, 750, 1000, 1500, 2500, 5000],
    "min_amplicon_size": 150,
    "ideal_amplicon_max": 700,
    "max_amplicon_size": 5000,
    "max_target_spans_per_flank": 30,
    "primer3_num_return": 250,
    "top_n_output": 20,
    "target_pass_field": "Junction_Target_Pass",
    "target_failure_reason": "junction_requirement_failed",
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
    "blast_specificity": True,
    "specificity_candidates": 60,
    "blast_min_unique_len": 18,
    "blast_min_unique_identity": 95.0,
    "blast_noise_len": 12,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Design SNP-aware cross-junction PCR primers for mature RNA-derived cDNA."
    )
    parser.add_argument("--gene", default=DEFAULT_GENE, help="Single gene symbol to design primers for.")
    parser.add_argument(
        "--genes",
        default="",
        help="Comma-separated gene symbols for batch design. Overrides --gene when provided.",
    )
    parser.add_argument(
        "--genes-file",
        default="",
        help="Optional text file with one gene symbol per line for batch design. Overrides --gene and --genes.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory for SNP-aware RNA primer output files.",
    )
    parser.add_argument("--gtf", default=DEFAULT_CONFIG["gtf_path"], help="mm10 refGene GTF.")
    parser.add_argument("--genome", default=DEFAULT_CONFIG["genome_fasta"], help="mm10 FASTA.")
    parser.add_argument("--blast-db", default=DEFAULT_CONFIG["blast_database"], help="mm10 BLAST DB prefix.")
    parser.add_argument("--snp-vcf", default=DEFAULT_CONFIG["snp_vcf"], help="B6/Cast SNP VCF.")
    parser.add_argument("--min-snps", type=int, default=DEFAULT_CONFIG["min_informative_snps"])
    parser.add_argument("--max-amplicon", type=int, default=DEFAULT_CONFIG["max_amplicon_size"])
    parser.add_argument(
        "--min-junction-overlap",
        type=int,
        default=DEFAULT_CONFIG["min_junction_overlap_bases"],
        help="Minimum primer bases required on each side of a crossed exon-exon junction.",
    )
    parser.add_argument(
        "--min-junction-3p-anchor",
        type=int,
        default=DEFAULT_CONFIG["min_junction_3p_anchor_bases"],
        help="Minimum primer bases required on the 3-prime side of a crossed exon-exon junction.",
    )
    parser.add_argument(
        "--flank-steps",
        default=",".join(str(x) for x in DEFAULT_CONFIG["flank_steps"]),
        help="Comma-separated cDNA flanking windows to try around target junction/SNP spans.",
    )
    parser.add_argument(
        "--skip-blast-specificity",
        action="store_true",
        help="Skip primer BLAST specificity annotation.",
    )
    parser.add_argument(
        "--specificity-candidates",
        type=int,
        default=DEFAULT_CONFIG["specificity_candidates"],
        help="Number of ranked deduplicated pairs to BLAST-annotate for specificity.",
    )
    parser.add_argument("--top-n", type=int, default=DEFAULT_CONFIG["top_n_output"])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = dict(DEFAULT_CONFIG)
    config.update(
        {
            "gtf_path": args.gtf,
            "genome_fasta": args.genome,
            "blast_database": args.blast_db,
            "snp_vcf": args.snp_vcf,
            "min_informative_snps": args.min_snps,
            "max_amplicon_size": args.max_amplicon,
            "min_junction_overlap_bases": args.min_junction_overlap,
            "min_junction_3p_anchor_bases": args.min_junction_3p_anchor,
            "flank_steps": [int(x.strip()) for x in args.flank_steps.split(",") if x.strip()],
            "blast_specificity": not args.skip_blast_specificity,
            "specificity_candidates": args.specificity_candidates,
            "top_n_output": args.top_n,
        }
    )

    genes = _resolve_genes(args)
    console.print("[cyan]Running SNP-aware RNA primer design[/cyan]")
    console.print(f"Genes: {', '.join(genes)}")
    console.print(f"Output: {args.output_dir}")
    console.print(f"Default conda environment: {config['conda_env']}")
    if len(genes) > 1:
        results = design_snpaware_rna_primers_for_genes(genes, args.output_dir, config)
        successes = sum(row["Status"] == "success" for row in results)
        console.print(f"[green]Batch complete.[/green] Successful genes: {successes}/{len(results)}")
        console.print(f"[green]Batch summary:[/green] {args.output_dir}/SNPaware_RNA_batch_summary.csv")
        return

    result = design_snpaware_rna_primers(genes[0], args.output_dir, config)
    selected = result["selected_transcript"]
    top_rows = result["top_rows"]
    console.print(
        f"[green]Selected transcript:[/green] {selected.transcript_id} "
        f"({selected.chrom}:{selected.start}-{selected.end} {selected.strand}, "
        f"{selected.spliced_length} nt spliced)"
    )
    console.print(f"[green]Common junctions:[/green] {len(result['common_junctions'])}")
    console.print(f"[green]Candidate primer pairs:[/green] {len(result['candidate_rows'])}")
    console.print(f"[green]Selected filter tier:[/green] {result['filter_info']['selected_tier']}")
    console.print(f"[green]Fallback used:[/green] {result['fallback_used']}")
    if top_rows:
        best = top_rows[0]
        console.print(
            "[green]Best pair:[/green] "
            f"{best['Left_Primer_Seq']} / {best['Right_Primer_Seq']} "
            f"({best['Amplicon_Size']} bp cDNA, {best['SNP_Count_In_Amplicon']} SNPs, "
            f"{best['Design_Mode']})"
        )
    console.print(f"[green]Done.[/green] Files written to {args.output_dir}")


def _resolve_genes(args: argparse.Namespace) -> list[str]:
    if args.genes_file:
        with open(args.genes_file) as handle:
            genes = [line.strip() for line in handle if line.strip() and not line.startswith("#")]
    elif args.genes:
        genes = [gene.strip() for gene in args.genes.split(",") if gene.strip()]
    else:
        genes = [args.gene]
    if not genes:
        raise ValueError("No genes provided.")
    return list(dict.fromkeys(genes))


if __name__ == "__main__":
    main()
