#!/usr/bin/env python3
# validate_probe_specificity.py - Phase 2 of the FISH-RT probe pipeline

import os
import sys
import argparse
import subprocess
import pandas as pd
import re
import traceback
from pathlib import Path
from rich.console import Console
from rich.progress import track
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Try to import config, but provide defaults if not found
try:
    sys.path.append(os.getcwd())
    from utils.config import FISH_RT_CONFIG
except ImportError:
    FISH_RT_CONFIG = {}

console = Console()

def run_blast(query_fasta, output_file, config):
    """Run local BLAST specificity check"""
    blast_db = config.get("blast_database_path", "")
    if not blast_db:
        # Check if it was provided in config but we had import error
        raise ValueError("BLAST database path not configured in config.py or config.py not found.")

    evalue = config.get("blast_max_evalue", 0.01)
    word_size = config.get("blast_word_size", 7)
    max_target_seqs = config.get("max_target_seqs", 100)

    console.print(f"[cyan]Running BLASTn against {blast_db}...[/cyan]")

    # Verify blastn exists
    try:
        subprocess.run(["blastn", "-version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        raise RuntimeError("BLAST+ ('blastn') not found in PATH. Please ensure it is installed and the environment is active.")

    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-db", blast_db,
        "-out", str(output_file),
        "-outfmt", "0",
        "-task", "blastn-short",
        "-max_target_seqs", str(max_target_seqs),
        "-evalue", str(evalue),
        "-word_size", str(word_size)
    ]

    console.print(f"[dim]Command: {' '.join(cmd)}[/dim]")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        console.print(f"[red]BLAST execution failed with exit code {e.returncode}[/red]")
        console.print(f"[red]Error output: {e.stderr}[/red]")
        raise

    console.print("[green]✅ BLAST completed successfully.[/green]")

def parse_blast_results(blast_text, min_alignment_length=15):
    """
    Parse BLAST result text and count High Scoring Pairs (HSPs).
    Only counts hits that meet the min_alignment_length threshold.
    """
    # Split by queries
    queries = re.split(r"Query #\d+:|Query=\s*", blast_text)[1:]
    records = []

    if not queries:
        console.print("[yellow]Warning: No query results found in BLAST output.[/yellow]")
        return pd.DataFrame()

    for query in track(queries, description="Parsing BLAST results..."):
        header_match = re.search(r"^(\S+)", query.strip())
        full_id = header_match.group(1) if header_match else "Unknown"
        gene_name = full_id.split('_')[0] if '_' in full_id else "Unknown"

        # Find all HSP segments and check their alignment length
        # Format: "Identities = 21/29 (72%), Gaps = 0/29 (0%)"
        hsp_identities = re.findall(r"Identities\s*=\s*(\d+)/(\d+)\s*\((\d+)%\)", query)

        # Filter for "good" hits only
        valid_hits = []
        for match_count, total_count, perc in hsp_identities:
            if int(total_count) >= min_alignment_length:
                valid_hits.append({
                    "length": int(total_count),
                    "identity": int(perc)
                })

        num_hits = len(valid_hits)

        # Reporting primary and secondary identities for transparency
        primary_id = valid_hits[0]["identity"] if num_hits > 0 else None
        secondary_id = valid_hits[1]["identity"] if num_hits > 1 else None

        records.append({
            "ProbeFullID": full_id,
            "GeneName": gene_name,
            "BLAST_Hits": num_hits,
            "Primary_Identity": primary_id,
            "Secondary_Identity": secondary_id,
            "BLAST_Unique": num_hits == 1
        })

    return pd.DataFrame(records)

def main():
    parser = argparse.ArgumentParser(description="Validate probe specificity using BLAST")
    parser.add_argument("-c", "--candidates", required=True, help="Candidates CSV file from Phase 1")
    parser.add_argument("-f", "--fasta", help="Candidates FASTA file (optional, inferred from CSV)")
    parser.add_argument("-o", "--output-dir", help="Output directory (optional, same as CSV)")
    parser.add_argument("--skip-blast", action="store_true", help="Skip running BLAST, use existing results file")

    args = parser.parse_args()

    console.print("[bold blue]🧬 Probe Specificity Validation (Phase 2)[/bold blue]")
    console.print("=" * 60)

    try:
        csv_path = Path(args.candidates)
        if not csv_path.exists():
            console.print(f"[red]Error: Candidates CSV not found: {csv_path}[/red]")
            return 1

        fasta_path = Path(args.fasta) if args.fasta else csv_path.with_suffix(".fasta")
        output_dir = Path(args.output_dir) if args.output_dir else csv_path.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        blast_results_path = output_dir / "blast_results.txt"

        # 1. Run BLAST
        if not args.skip_blast:
            if not fasta_path.exists():
                console.print(f"[red]Error: Candidates FASTA not found at {fasta_path}[/red]")
                return 1

            # Research-Validated parameters for short probes
            # word_size=11 and evalue=0.1 capture dangerous 21bp matches
            # while avoiding noise from marginal 11-14bp alignments.
            config_copy = FISH_RT_CONFIG.copy()
            config_copy["blast_max_evalue"] = 0.1
            config_copy["blast_word_size"] = 11
            config_copy["max_target_seqs"] = 5

            run_blast(fasta_path, blast_results_path, config_copy)
        else:
            console.print("[yellow]Skipping BLAST execution, using existing results.[/yellow]")

        # 2. Parse Results
        if not blast_results_path.exists():
            console.print(f"[red]Error: BLAST results file not found at {blast_results_path}[/red]")
            return 1

        console.print("[cyan]Reading BLAST results...[/cyan]")
        with open(blast_results_path, "r") as f:
            blast_text = f.read()

        blast_results_df = parse_blast_results(blast_text)
        if blast_results_df.empty:
            console.print("[red]Error: Failed to parse any data from BLAST results.[/red]")
            return 1

        console.print(f"✅ Parsed {len(blast_results_df)} results.")

        # 3. Load candidates and Merge
        console.print(f"[cyan]Merging with {csv_path.name}...[/cyan]")
        candidates_df = pd.read_csv(csv_path)

        # DROPPING existing BLAST columns to avoid suffixing (_x, _y)
        # This is the fix for the KeyError: 'BLAST_Unique'
        # DROPPING existing BLAST columns to avoid suffixing (_x, _y)
        cols_to_drop = [
            "BLAST_Hits", "BLAST_Identity", "BLAST_Unique", "BLAST_Hit_Name",
            "Primary_Identity", "Secondary_Identity", "BLAST_WordSize",
            "BLAST_EValue", "BLAST_MinAlignment"
        ]
        for col in cols_to_drop:
            if col in candidates_df.columns:
                console.print(f"[dim]Dropping existing '{col}' from candidates...[/dim]")
                candidates_df = candidates_df.drop(columns=[col])

        # Recreate ProbeID matching Phase 1 output_generator.py if not present
        if "ProbeID" not in candidates_df.columns:
            candidates_df["ProbeID"] = candidates_df["GeneName"] + "_probe_" + candidates_df.index.astype(str)

        # Merge only on ProbeID. Drop GeneName from blast results to avoid duplicate
        if "GeneName" in blast_results_df.columns:
            blast_results_df = blast_results_df.drop(columns=["GeneName"])

        # Debug: Check column presence
        console.print(f"[dim]Candidates columns: {list(candidates_df.columns)}[/dim]")
        console.print(f"[dim]BLAST results columns: {list(blast_results_df.columns)}[/dim]")

        merged_df = candidates_df.merge(
            blast_results_df[["ProbeFullID", "BLAST_Hits", "Primary_Identity", "Secondary_Identity", "BLAST_Unique"]],
            left_on="ProbeID",
            right_on="ProbeFullID",
            how="left"
        )

        # Add parameters to the results dataframe for documentation
        # This fulfills the user request to include parameters in the CSV
        merged_df["BLAST_WordSize"] = 11
        merged_df["BLAST_EValue"] = 0.1
        merged_df["BLAST_MinAlignment"] = 15

        # Handle cases where no hits were found at all (NaN)
        # Old behavior: only keep probes with exactly 1 hit
        # So we fill NaN with False (non-specific/not found)
        num_no_results = merged_df["BLAST_Unique"].isna().sum()
        if num_no_results > 0:
            console.print(f"[yellow]Warning: {num_no_results} probes had no BLAST results found (will be filtered out).[/yellow]")

        merged_df["BLAST_Unique"] = merged_df["BLAST_Unique"].fillna(False)
        merged_df["BLAST_Hits"] = merged_df["BLAST_Hits"].fillna(0)

        # Filter for unique
        final_df = merged_df[merged_df["BLAST_Unique"] == True].copy()

        # 4. Save BLAST results summary for all candidates (requested)
        results_csv = output_dir / "FISH_RT_probes_CANDIDATES_BLASTresults.csv"
        # Using merged_df ensures we have ProbeID, GeneName, and the BLAST parameters
        merged_df.to_csv(results_csv, index=False)
        console.print(f"✅ [green]Detailed BLAST results saved to: {results_csv}[/green]")

        # 5. Save Final Selection
        final_csv = output_dir / "FISH_RT_probes_FINAL_SELECTION.csv"
        final_df.to_csv(final_csv, index=False)
        console.print(f"✅ [green]Final selection saved to: {final_csv}[/green]")

        # 5. Save Summary
        summary_txt = output_dir / "FISH_RT_probes_FINAL_SELECTION_summary.txt"
        with open(summary_txt, "w") as f:
            f.write("Probe Specificity Validation Summary\n")
            f.write("=" * 40 + "\n")
            f.write(f"Total Candidates: {len(candidates_df)}\n")
            f.write(f"BLAST Unique: {len(final_df)}\n")
            f.write(f"Retention Rate: {100*len(final_df)/len(candidates_df):.1f}%\n")

        console.print(f"\n[bold]Summary:[/bold]")
        console.print(f"  Total Candidates: {len(candidates_df)}")
        console.print(f"  BLAST Unique: {len(final_df)} ({100*len(final_df)/len(candidates_df):.1f}%)")
        console.print(f"  Summary saved to: {summary_txt.name}")

    except Exception as e:
        console.print(f"[bold red]An unexpected error occurred:[/bold red]")
        console.print(f"[red]{str(e)}[/red]")
        console.print(traceback.format_exc())
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())
