#!/usr/bin/env python3
"""Design KOveri genotyping primers from CRISPR gRNA FASTA input."""

from __future__ import annotations

import argparse
from pathlib import Path

from rich.console import Console

from config import DEFAULT_INPUT_FASTA, DEFAULT_OUTPUT_DIR, KOVERI_CONFIG
from utils.designer import design_koveri_primers


console = Console()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Design unbiased B6/Cast genotyping PCR primers around CRISPR gRNAs."
    )
    parser.add_argument(
        "--input",
        default=None,
        help=(
            "Input gRNA FASTA file. Defaults to RNF12_gRNAs.fa inside --output-dir, "
            "so inputs and results live in the same folder."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory for KOveri output files.",
    )
    parser.add_argument("--genome", default=KOVERI_CONFIG["genome_fasta"], help="mm10 FASTA.")
    parser.add_argument("--blast-db", default=KOVERI_CONFIG["blast_database"], help="mm10 BLAST DB prefix.")
    parser.add_argument("--snp-vcf", default=KOVERI_CONFIG["snp_vcf"], help="B6/Cast SNP VCF.")
    parser.add_argument("--min-snps", type=int, default=KOVERI_CONFIG["min_informative_snps"])
    parser.add_argument("--max-amplicon", type=int, default=KOVERI_CONFIG["max_amplicon_size"])
    parser.add_argument(
        "--flank-steps",
        default=",".join(str(x) for x in KOVERI_CONFIG["flank_steps"]),
        help="Comma-separated flanking windows to try around all gRNAs.",
    )
    parser.add_argument(
        "--skip-blast-specificity",
        action="store_true",
        help="Skip primer BLAST specificity annotation.",
    )
    parser.add_argument(
        "--specificity-candidates",
        type=int,
        default=KOVERI_CONFIG["specificity_candidates"],
        help="Number of ranked deduplicated pairs to BLAST-annotate for specificity.",
    )
    parser.add_argument("--top-n", type=int, default=KOVERI_CONFIG["top_n_output"])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = dict(KOVERI_CONFIG)
    config.update(
        {
            "genome_fasta": args.genome,
            "blast_database": args.blast_db,
            "snp_vcf": args.snp_vcf,
            "min_informative_snps": args.min_snps,
            "max_amplicon_size": args.max_amplicon,
            "flank_steps": [int(x.strip()) for x in args.flank_steps.split(",") if x.strip()],
            "blast_specificity": not args.skip_blast_specificity,
            "specificity_candidates": args.specificity_candidates,
            "top_n_output": args.top_n,
        }
    )

    input_path = Path(args.input) if args.input else Path(args.output_dir) / DEFAULT_INPUT_FASTA.name
    if not input_path.is_absolute():
        input_path = Path(__file__).resolve().parent / input_path

    console.print("[cyan]Running KOveri primer design[/cyan]")
    console.print(f"Input: {input_path}")
    console.print(f"Output: {args.output_dir}")
    result = design_koveri_primers(input_path, args.output_dir, config)
    target = result["target"]
    top_rows = result["top_rows"]
    console.print(
        f"[green]Selected target interval:[/green] {target.chrom}:{target.start}-{target.end}"
    )
    console.print(f"[green]Candidate primer pairs:[/green] {len(result['candidate_rows'])}")
    console.print(f"[green]Selected filter tier:[/green] {result['filter_info']['selected_tier']}")
    if top_rows:
        best = top_rows[0]
        console.print(
            "[green]Best pair:[/green] "
            f"{best['Left_Primer_Seq']} / {best['Right_Primer_Seq']} "
            f"({best['Amplicon_Size']} bp, {best['SNP_Count_In_Amplicon']} SNPs)"
        )
    console.print(f"[green]Done.[/green] Files written to {args.output_dir}")


if __name__ == "__main__":
    main()
