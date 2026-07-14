#!/usr/bin/env python3
"""Design allele-aware gene-specific primers to pair with an RT/poly-dT handle."""

from __future__ import annotations

import argparse
from pathlib import Path

from rich.console import Console

from utils.designer import design_allele_primers, design_allele_primers_for_genes


console = Console()

DEFAULT_GENES = "Xist,Kdm5c"
DEFAULT_OUTPUT_DIR = Path("./allele_primer_results")

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
    "curated_transcript_prefixes": ("NM_", "NR_"),
    "top_snps": 5,
    "snp_window_toward_polydt": 100,
    "primer3_num_return": 50,
    "valid_primer_length_min": 18,
    "valid_primer_length_max": 30,
    "post_filter_gc_min": 40.0,
    "post_filter_gc_max": 60.0,
    "post_filter_tm_min": 55.0,
    "post_filter_tm_max": 65.0,
    "max_homopolymer_run": 3,
    "max_quad_g_run": 3,
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
    "blast_specificity": True,
    "blast_min_unique_len": 18,
    "blast_min_unique_identity": 95.0,
    "blast_noise_len": 12,
    "blast_candidates_per_snp": 5,
    "handle": "2PUNI",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Design allele-aware primers to pair with an RT/poly-dT handle.")
    parser.add_argument("--gene", default="", help="Single gene symbol. Ignored when --genes or --genes-file is set.")
    parser.add_argument("--genes", default="", help="Comma-separated gene symbols. Defaults to Xist,Kdm5c.")
    parser.add_argument("--genes-file", default="", help="Optional text file with one gene per line.")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="Output directory or output root.")
    parser.add_argument("--gtf", default=DEFAULT_CONFIG["gtf_path"], help="mm10 refGene GTF.")
    parser.add_argument("--genome", default=DEFAULT_CONFIG["genome_fasta"], help="mm10 FASTA.")
    parser.add_argument("--blast-db", default=DEFAULT_CONFIG["blast_database"], help="mm10 BLAST DB prefix.")
    parser.add_argument("--snp-vcf", default=DEFAULT_CONFIG["snp_vcf"], help="B6/Cast SNP VCF.")
    parser.add_argument("--samtools", default=DEFAULT_CONFIG["samtools"])
    parser.add_argument("--blastn", default=DEFAULT_CONFIG["blastn"])
    parser.add_argument("--primer3-core", default=DEFAULT_CONFIG["primer3_core"])
    parser.add_argument("--top-snps", type=int, default=DEFAULT_CONFIG["top_snps"])
    parser.add_argument("--snp-window-toward-polydt", type=int, default=DEFAULT_CONFIG["snp_window_toward_polydt"])
    parser.add_argument("--primer-min-size", type=int, default=DEFAULT_CONFIG["primer_min_size"])
    parser.add_argument("--primer-opt-size", type=int, default=DEFAULT_CONFIG["primer_opt_size"])
    parser.add_argument("--primer-max-size", type=int, default=DEFAULT_CONFIG["primer_max_size"])
    parser.add_argument("--primer-min-tm", type=float, default=DEFAULT_CONFIG["primer_min_tm"])
    parser.add_argument("--primer-opt-tm", type=float, default=DEFAULT_CONFIG["primer_opt_tm"])
    parser.add_argument("--primer-max-tm", type=float, default=DEFAULT_CONFIG["primer_max_tm"])
    parser.add_argument("--primer-min-gc", type=float, default=DEFAULT_CONFIG["primer_min_gc"])
    parser.add_argument("--primer-max-gc", type=float, default=DEFAULT_CONFIG["primer_max_gc"])
    parser.add_argument("--skip-blast-specificity", action="store_true", help="Skip short-primer BLAST specificity.")
    parser.add_argument("--handle", choices=("2PUNI", "2PBC"), default="2PUNI", help="Common 5-prime handle to add to every accepted insert (default: 2PUNI).")
    parser.add_argument("--blast-candidates-per-snp", type=int, default=5, help="Maximum top QC-passing candidates to BLAST per SNP target (default: 5).")
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
            "samtools": args.samtools,
            "blastn": args.blastn,
            "primer3_core": args.primer3_core,
            "top_snps": args.top_snps,
            "snp_window_toward_polydt": args.snp_window_toward_polydt,
            "primer_min_size": args.primer_min_size,
            "primer_opt_size": args.primer_opt_size,
            "primer_max_size": args.primer_max_size,
            "valid_primer_length_min": args.primer_min_size,
            "valid_primer_length_max": args.primer_max_size,
            "primer_min_tm": args.primer_min_tm,
            "primer_opt_tm": args.primer_opt_tm,
            "primer_max_tm": args.primer_max_tm,
            "post_filter_tm_min": args.primer_min_tm,
            "post_filter_tm_max": args.primer_max_tm,
            "primer_min_gc": args.primer_min_gc,
            "primer_max_gc": args.primer_max_gc,
            "post_filter_gc_min": args.primer_min_gc,
            "post_filter_gc_max": args.primer_max_gc,
            "blast_specificity": not args.skip_blast_specificity,
            "handle": args.handle,
            "blast_candidates_per_snp": args.blast_candidates_per_snp,
        }
    )

    genes = _resolve_genes(args)
    console.print("[cyan]Running allele primer design to pair RT/poly-dT primer[/cyan]")
    console.print(f"Genes: {', '.join(genes)}")
    console.print(f"Output: {args.output_dir}")
    console.print(f"SNP window toward poly-dT: {config['snp_window_toward_polydt']} nt")
    console.print(f"Handle: {config['handle']}")

    if len(genes) > 1:
        results = design_allele_primers_for_genes(genes, args.output_dir, config)
        successes = sum(row["Status"] == "success" for row in results)
        console.print(f"[green]Batch complete.[/green] Successful genes: {successes}/{len(results)}")
        console.print(f"[green]Batch summary:[/green] {args.output_dir}/allele_primer_batch_summary.csv")
        console.print(f"[green]Consolidated oligos:[/green] {args.output_dir}/allele_primer_consolidated_top_oligos.csv")
        console.print(f"[green]Consolidated report:[/green] {args.output_dir}/allele_primer_consolidated_report.html")
        return

    result = design_allele_primers(genes[0], args.output_dir, config)
    selected = result["selected_transcript"]
    console.print(
        f"[green]Selected transcript:[/green] {selected.transcript_id} "
        f"({selected.chrom}:{selected.start}-{selected.end} {selected.strand}, "
        f"{selected.spliced_length} nt spliced)"
    )
    console.print(f"[green]Informative SNPs:[/green] {len(result['snps'])}")
    console.print(f"[green]Candidate inserts:[/green] {len(result['candidate_rows'])}")
    console.print(f"[green]Top inserts:[/green] {len(result['top_insert_rows'])}")
    console.print(f"[green]Top oligos:[/green] {len(result['top_oligo_rows'])}")
    console.print(f"[green]Done.[/green] Files written to {args.output_dir}")


def _resolve_genes(args: argparse.Namespace) -> list[str]:
    if args.genes_file:
        with open(args.genes_file) as handle:
            genes = [line.strip() for line in handle if line.strip() and not line.startswith("#")]
    elif args.genes:
        genes = [gene.strip() for gene in args.genes.split(",") if gene.strip()]
    elif args.gene:
        genes = [args.gene]
    else:
        genes = [gene.strip() for gene in DEFAULT_GENES.split(",")]
    if not genes:
        raise ValueError("No genes provided.")
    return list(dict.fromkeys(genes))


if __name__ == "__main__":
    main()
