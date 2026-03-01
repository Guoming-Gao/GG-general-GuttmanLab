#!/usr/bin/env python3

# main.py - smfish-like-rt-probe-designer

import os
import sys
import argparse
import pandas as pd
from pathlib import Path
from rich.progress import track
from rich.console import Console

# Import our modules
from utils.config import FISH_RT_CONFIG
from utils.gene_fetcher import GeneSequenceFetcher
from utils.oligostan_core import design_fish_probes
from utils.snp_analyzer import SNPCoverageAnalyzer
from utils.output_generator import OutputGenerator
from utils.filters import is_ok_4_homopolymer

console = Console()


def setup_output_directory(output_path):
    """Create output directory structure"""
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    (output_path / "gene_sequences").mkdir(exist_ok=True)
    return output_path


def apply_basic_quality_filtering(probes, config):
    """
    Apply basic quality filtering (GC, dustmasker) without selection.
    PNAS filtering is adaptive and happens after SNP analysis.
    """
    console.print(
        f"[cyan]Applying basic quality filtering to {len(probes)} probes...[/cyan]"
    )

    filtered_probes = []

    for probe in probes:
        gc_pass = probe["GCFilter"] == 1
        dustmasker_pass = probe.get("MaskedFilter", 1) == 1

        # Homo-polymer filter
        homopolymer_pass = is_ok_4_homopolymer(
            probe["Probe_Seq"],
            config.get("max_homopolymer_length", 4)
        )

        if gc_pass and dustmasker_pass and homopolymer_pass:
            filtered_probes.append(probe)

    console.print(f"[green]Basic Filtering Results:[/green]")
    console.print(f" Total probes: {len(probes)}")
    console.print(f" Passed GC + dustmasker: {len(filtered_probes)}")
    console.print(f" Retention rate: {100*len(filtered_probes)/len(probes):.1f}%")

    return filtered_probes


def select_top_probes_by_snp(probes, config, max_per_gene=None, min_snps=None):
    """
    After SNP analysis, select probes per gene using the logic:
    min(top N probes, all probes with >= M SNPs)

    Sorted by:
    1. SNP_Count (desc) - most important for allelic discrimination
    2. NbOfPNAS (desc) - probe quality
    3. dGScore (desc) - thermodynamic score

    Also applies adaptive PNAS filtering if needed.
    """
    # Use config values or provided overrides
    max_per_gene = max_per_gene or config.get("max_probes_per_gene", 200)
    min_snps = min_snps if min_snps is not None else config.get("min_snps_for_selection", 3)

    console.print(
        f"[cyan]Selecting probes: min(top {max_per_gene}, probes with ≥{min_snps} SNPs)...[/cyan]"
    )

    # Group probes by gene
    probes_by_gene = {}
    for probe in probes:
        gene = probe["GeneName"]
        if gene not in probes_by_gene:
            probes_by_gene[gene] = []
        probes_by_gene[gene].append(probe)

    # PNAS rule progression: start strict, progressively relax
    pnas_rule_sets = [
        [1, 2, 4],       # Strict: 3 rules
        [1, 2, 3, 4],    # Medium: 4 rules
        [1, 2, 3, 4, 5], # Relaxed: 5 rules (all)
        [],              # Last resort: no PNAS filter
    ]

    selected_probes = []
    gene_stats = {}

    for gene, gene_probes in probes_by_gene.items():
        best_probes = []
        rules_used = None

        for pnas_rules in pnas_rule_sets:
            # Filter with current PNAS rule set
            passed = []
            for probe in gene_probes:
                if len(pnas_rules) > 0:
                    pnas_pass = probe["NbOfPNAS"] >= len(pnas_rules)
                else:
                    pnas_pass = True  # No PNAS filter

                if pnas_pass:
                    passed.append(probe)

            if len(passed) >= max_per_gene:
                best_probes = passed
                rules_used = pnas_rules
                break
            elif len(passed) > len(best_probes):
                best_probes = passed
                rules_used = pnas_rules

        # Sort by SNP_Count (primary), then NbOfPNAS, then dGScore
        best_probes.sort(
            key=lambda p: (p.get("SNP_Count", 0), p["NbOfPNAS"], p["dGScore"]),
            reverse=True
        )

        # Apply new selection logic: min(top N, probes with >= M SNPs)
        # First, get probes meeting SNP threshold
        probes_with_snps = [p for p in best_probes if p.get("SNP_Count", 0) >= min_snps]
        top_n = best_probes[:max_per_gene]

        # Take whichever is smaller (but at least get probes with sufficient SNPs)
        if len(probes_with_snps) <= len(top_n):
            selected = probes_with_snps
        else:
            selected = top_n

        selected_probes.extend(selected)

        # Stats for reporting
        avg_snps = sum(p.get("SNP_Count", 0) for p in selected) / len(selected) if selected else 0
        gene_stats[gene] = {
            "total": len(gene_probes),
            "passed_pnas": len(best_probes),
            "selected": len(selected),
            "rules_used": rules_used,
            "avg_snps": avg_snps,
        }

    # Report statistics
    console.print(f"[green]Probe Selection Results:[/green]")
    console.print(f" Max per gene: {max_per_gene}, Min SNPs: {min_snps}")
    console.print(f" Total genes: {len(probes_by_gene)}")
    console.print(f" Total selected: {len(selected_probes)}")
    console.print(f"\n Per-gene breakdown:")
    for gene in sorted(gene_stats.keys()):
        stats = gene_stats[gene]
        rules_str = ",".join(map(str, stats["rules_used"])) if stats["rules_used"] else "none"
        console.print(f"   {gene}: {stats['selected']}/{stats['passed_pnas']} (PNAS [{rules_str}], avgSNP: {stats['avg_snps']:.1f})")

    return selected_probes


def main():
    """Main FISH-RT probe design pipeline"""
    console.print("[bold blue]🧬 smfish-like-rt-probe-designer[/bold blue]")
    console.print("=" * 60)
    console.print("[cyan]Mode: Local GTF + FASTA files only[/cyan]")

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Design high-quality FISH probes using local files"
    )
    parser.add_argument(
        "--output", help="Output directory path (optional, uses config default)"
    )
    parser.add_argument("--genes", nargs="+", help="List of gene names to process")
    parser.add_argument(
        "--max-probes", type=int, default=None,
        help="Max probes per gene (default: config value or 200)"
    )
    parser.add_argument(
        "--min-snps", type=int, default=None,
        help="Min SNPs for probe selection (default: config value or 3)"
    )
    parser.add_argument(
        "--rt-coverage", type=int, default=None,
        help="RT coverage length in nt (default: config value or 200)"
    )
    args = parser.parse_args()

    # Load configuration
    config = FISH_RT_CONFIG.copy()

    # Use config default if --output not provided
    if args.output:
        config["output_directory"] = args.output
        console.print(f"[yellow]Using custom output directory: {args.output}[/yellow]")
    if args.rt_coverage:
        config["rt_coverage_downstream"] = args.rt_coverage
        console.print(f"[yellow]Using custom RT coverage: {args.rt_coverage} nt[/yellow]")

    if not args.output and not args.rt_coverage:
        console.print(f"[green]Using config default settings[/green]")

    # Determine gene list
    if args.genes:
        gene_list = args.genes
        console.print(
            f"[green]Processing {len(gene_list)} specified genes[/green]"
        )
    else:
        # Default fallback if no genes provided
        gene_list = ["Mecp2"]
        console.print(f"[yellow]No genes specified. Using default: {gene_list}[/yellow]")
        console.print(f"[dim]Tip: Use --genes GENE1 GENE2 ... to specify your targets[/dim]")

    # Setup output directory
    output_path = setup_output_directory(config["output_directory"])
    console.print(f"[green]Output directory: {output_path}[/green]")

    # Show configuration
    console.print(f"\n[bold]Configuration:[/bold]")
    console.print(f" PNAS rules: {config['pnas_filter_rules']}")
    console.print(
        f" Dustmasker: {'✅ Enabled' if config['use_dustmasker'] else '❌ Disabled'}"
    )
    console.print(f" RT coverage: {config['rt_coverage_downstream']} nt downstream")
    # console.print(f" Min SNP coverage: {config.get('min_snp_coverage_for_final', 'N/A')}")

    # Initialize components
    console.print("\n[bold]Initializing pipeline components...[/bold]")
    gene_fetcher = GeneSequenceFetcher(config)
    snp_analyzer = SNPCoverageAnalyzer(config) if config["snp_file_path"] else None
    output_generator = OutputGenerator(config)

    # Step 1: Fetch gene sequences
    console.print("\n[bold cyan]Step 1: Fetching gene sequences...[/bold cyan]")
    all_gene_data = []
    failed_genes = []

    for gene_name in track(gene_list, description="Fetching sequences..."):
        try:
            gene_data = gene_fetcher.fetch_gene_sequence(gene_name)
            if gene_data:
                all_gene_data.append(gene_data)
                console.print(
                    f"✅ {gene_name}: {len(gene_data['sequence'])} bp ({gene_data['source']}) strand {gene_data['strand']}"
                )
            else:
                failed_genes.append(gene_name)
                console.print(f"❌ {gene_name}: Failed to fetch")
        except Exception as e:
            failed_genes.append(gene_name)
            console.print(f"❌ {gene_name}: Error - {str(e)}")

    if failed_genes:
        console.print(
            f"[red]Failed to fetch {len(failed_genes)} genes: {failed_genes}[/red]"
        )

    if not all_gene_data:
        console.print("[red]No genes successfully fetched. Exiting.[/red]")
        return

    # Save gene sequences (for validation)
    gene_fetcher.save_gene_sequences(all_gene_data, output_path / "gene_sequences")

    # Step 2: Design FISH probes
    console.print("\n[bold cyan]Step 2: Designing FISH probes...[/bold cyan]")
    all_probes = []

    for gene_data in track(all_gene_data, description="Designing probes..."):
        try:
            probes = design_fish_probes(gene_data, config)
            if probes:
                all_probes.extend(probes)
                console.print(f"✅ {gene_data['gene_name']}: {len(probes)} probes")
            else:
                console.print(f"⚠️ {gene_data['gene_name']}: No probes generated")
        except Exception as e:
            console.print(f"❌ {gene_data['gene_name']}: Error - {str(e)}")

    if not all_probes:
        console.print("[red]No probes generated. Exiting.[/red]")
        return

    console.print(f"[green]Total probes designed: {len(all_probes)}[/green]")

    # Step 3: Basic quality filtering (GC + dustmasker only)
    console.print("\n[bold cyan]Step 3: Basic quality filtering...[/bold cyan]")
    filtered_probes = apply_basic_quality_filtering(all_probes, config)

    if not filtered_probes:
        console.print(
            "[red]No probes passed basic filtering. Consider relaxing filter parameters.[/red]"
        )
        return

    console.print(f"[green]Filtered probes: {len(filtered_probes)}[/green]")

    # Step 4: SNP analysis on ALL filtered probes (before selection!)
    if snp_analyzer:
        console.print(
            f"\n[bold cyan]Step 4: SNP analysis (all filtered probes)...[/bold cyan]"
        )
        try:
            filtered_probes = snp_analyzer.analyze_probes(filtered_probes)

            # Report SNP coverage statistics
            snp_counts = [p.get("SNP_Count", 0) for p in filtered_probes]
            avg_snps = sum(snp_counts) / len(snp_counts) if snp_counts else 0
            max_snps = max(snp_counts) if snp_counts else 0
            high_snp_count = sum(
                1
                for count in snp_counts
                if count >= config["min_snp_coverage_for_final"]
            )

            console.print(f"✅ SNP coverage analysis completed")
            console.print(f" Average SNPs per probe: {avg_snps:.1f}")
            console.print(f" Maximum SNPs per probe: {max_snps}")
            console.print(
                f" Probes with ≥{config['min_snp_coverage_for_final']} SNPs: {high_snp_count}"
            )

        except Exception as e:
            console.print(f"❌ SNP analysis failed: {str(e)}")
    else:
        console.print(
            "\n[yellow]Step 4: SNP analysis skipped (no SNP file configured)[/yellow]"
        )

    # Step 4b: Select top probes per gene (sorted by SNP_Count!)
    console.print("\n[bold cyan]Step 4b: Selecting top probes by SNP coverage...[/bold cyan]")
    high_quality_probes = select_top_probes_by_snp(
        filtered_probes, config,
        max_per_gene=args.max_probes,
        min_snps=args.min_snps
    )
    console.print(f"[green]High-quality probes selected: {len(high_quality_probes)}[/green]")

    # Step 5: Generate Candidate output files
    console.print(
        f"\n[bold cyan]Step 5: Generating candidate output files...[/bold cyan]"
    )

    try:
        # Generate Candidate Probes
        candidate_files = output_generator.generate_focused_outputs(
            high_quality_probes, output_path, file_prefix="FISH_RT_probes_CANDIDATES"
        )

        console.print("[bold green]✅ Phase 1: Candidate Design completed successfully![/bold green]")
        console.print(f"[green]Candidate files generated:[/green]")
        for file_path in candidate_files:
            console.print(f" 📄 {file_path}")

        # Final summary
        if candidate_files:
            csv_files = [f for f in candidate_files if str(f).endswith(".csv")]
            if csv_files:
                df_cand = pd.read_csv(csv_files[0])
                console.print(f"\n[bold]📊 Design Summary:[/bold]")
                console.print(f" Candidate probes: {len(df_cand)}")
                console.print(f" Genes processed: {len(set(df_cand['GeneName']))}")
                console.print(f" Output directory: {output_path}")

                if "SNP_Count" in df_cand.columns:
                    avg_snps = df_cand["SNP_Count"].mean()
                    console.print(f" Average SNPs covered: {avg_snps:.1f}")

        console.print("\n[bold yellow]Next Step:[/bold yellow] Run BLAST analysis on the candidates FASTA file to ensure specificity.")

    except Exception as e:
        console.print(f"[red]❌ Output generation failed: {str(e)}[/red]")


if __name__ == "__main__":
    main()
