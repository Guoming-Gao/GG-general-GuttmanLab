#!/usr/bin/env python3
"""
Step 1: FASTQ Quality Control with fastp
- Run fastp on all paired-end FASTQ files
- Generate HTML and JSON reports
- Create consolidated summary CSV
"""

import os
import subprocess
import json
import pandas as pd
from pathlib import Path
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn
from rich.table import Table

console = Console()

# Paths
DATA_DIR = "/Volumes/guttman/users/gmgao/Data_seq/20260207-DoxSeqAllRepsMultiplexed"
RESULTS_DIR = os.path.join(DATA_DIR, "results")
QC_DIR = os.path.join(RESULTS_DIR, "qc", "fastp")

def find_fastq_pairs():
    """Find all R1/R2 FASTQ pairs"""
    r1_files = sorted(Path(DATA_DIR).rglob("*_R1.fastq.gz"))
    pairs = []

    for r1 in r1_files:
        r2 = r1.parent / r1.name.replace("_R1.fastq.gz", "_R2.fastq.gz")
        if r2.exists():
            sample_name = r1.name.replace("_R1.fastq.gz", "")
            pairs.append({
                'sample': sample_name,
                'r1': str(r1),
                'r2': str(r2)
            })

    return pairs

def run_fastp(sample_info, progress, task_id):
    """Run fastp on a single sample"""
    sample = sample_info['sample']
    r1 = sample_info['r1']
    r2 = sample_info['r2']

    # Output files
    html_out = os.path.join(QC_DIR, f"{sample}_fastp.html")
    json_out = os.path.join(QC_DIR, f"{sample}_fastp.json")

    # Run fastp (no trimming, just QC)
    cmd = [
        "fastp",
        "-i", r1,
        "-I", r2,
        "--disable_adapter_trimming",
        "--disable_quality_filtering",
        "--disable_length_filtering",
        "-h", html_out,
        "-j", json_out,
        "--thread", "4"
    ]

    progress.update(task_id, description=f"[cyan]Processing {sample}...[/cyan]")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        console.print(f"[red]✗ Error processing {sample}[/red]")
        console.print(f"[dim]{result.stderr}[/dim]")
        return None

    progress.update(task_id, advance=1)
    return json_out

def parse_fastp_json(json_file):
    """Extract key metrics from fastp JSON output"""
    with open(json_file) as f:
        data = json.load(f)

    sample_name = Path(json_file).name.replace("_fastp.json", "")

    summary = data.get('summary', {})
    before = summary.get('before_filtering', {})
    after = summary.get('after_filtering', {})

    return {
        'Sample': sample_name,
        'Total_Reads': before.get('total_reads', 0),
        'Total_Bases': before.get('total_bases', 0),
        'Q20_Rate': before.get('q20_rate', 0),
        'Q30_Rate': before.get('q30_rate', 0),
        'GC_Content': before.get('gc_content', 0),
        'Read1_Mean_Length': data.get('read1_before_filtering', {}).get('total_bases', 0) // max(data.get('read1_before_filtering', {}).get('total_reads', 1), 1),
        'Read2_Mean_Length': data.get('read2_before_filtering', {}).get('total_bases', 0) // max(data.get('read2_before_filtering', {}).get('total_reads', 1), 1),
        'Duplication_Rate': data.get('duplication', {}).get('rate', 0)
    }

def main():
    console.print("\n[bold cyan]═══ Step 1: FASTQ Quality Control (fastp) ═══[/bold cyan]\n")

    # Find FASTQ pairs
    console.print("[yellow]Searching for FASTQ files...[/yellow]")
    pairs = find_fastq_pairs()
    console.print(f"  [green]✓[/green] Found {len(pairs)} sample pairs\n")

    if len(pairs) == 0:
        console.print("[red]✗ No FASTQ files found![/red]")
        return

    # Create output directory
    os.makedirs(QC_DIR, exist_ok=True)

    # Run fastp with progress bar
    json_files = []
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
        console=console
    ) as progress:
        task = progress.add_task("[cyan]Running fastp...", total=len(pairs))

        for pair in pairs:
            json_file = run_fastp(pair, progress, task)
            if json_file:
                json_files.append(json_file)

    console.print(f"\n[green]✓ Completed fastp QC for {len(json_files)}/{len(pairs)} samples[/green]\n")

    # Parse all JSON files and create summary
    console.print("[yellow]Creating summary report...[/yellow]")
    summaries = [parse_fastp_json(jf) for jf in json_files]
    df = pd.DataFrame(summaries)

    # Save summary CSV
    summary_csv = os.path.join(QC_DIR, "fastp_summary.csv")
    df.to_csv(summary_csv, index=False)
    console.print(f"  [green]✓[/green] Summary saved to fastp_summary.csv\n")

    # Print summary table
    table = Table(title="FASTQ QC Summary", show_header=True, header_style="bold magenta")
    table.add_column("Sample", style="cyan")
    table.add_column("Total Reads", justify="right")
    table.add_column("Q30 Rate", justify="right")
    table.add_column("GC%", justify="right")
    table.add_column("R1 Length", justify="right")
    table.add_column("R2 Length", justify="right")

    for _, row in df.iterrows():
        table.add_row(
            row['Sample'],
            f"{row['Total_Reads']:,}",
            f"{row['Q30_Rate']:.1%}",
            f"{row['GC_Content']:.1%}",
            f"{row['Read1_Mean_Length']:.0f}bp",
            f"{row['Read2_Mean_Length']:.0f}bp"
        )

    console.print(table)

    console.print(f"\n[bold green]✓ FASTQ QC complete![/bold green]")
    console.print(f"\n[dim]Output directory: {QC_DIR}[/dim]")
    console.print(f"[dim]  - {len(pairs)} HTML reports[/dim]")
    console.print(f"[dim]  - {len(pairs)} JSON reports[/dim]")
    console.print(f"[dim]  - 1 summary CSV[/dim]\n")

if __name__ == "__main__":
    main()
