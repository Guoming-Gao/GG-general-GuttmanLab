#!/usr/bin/env python3
"""
Step 2: Bowtie2 Alignment (Dual-Mode)
- Align all samples to original reference
- Align all samples to N-masked reference
- Generate sorted BAM files with indices
- Use default bowtie2 parameters (balanced filtering)
"""

import os
import subprocess
import pandas as pd
from pathlib import Path
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn, TaskID

console = Console()

# Paths
DATA_DIR = "/Volumes/guttman/users/gmgao/Data_seq/20260207-DoxSeqAllRepsMultiplexed"
RESULTS_DIR = os.path.join(DATA_DIR, "results")
REF_DIR = os.path.join(RESULTS_DIR, "references")

# Bowtie2 indices
ORIGINAL_INDEX = os.path.join(REF_DIR, "xist_amplicon_original")
MASKED_INDEX = os.path.join(REF_DIR, "xist_amplicon_masked")

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

def run_bowtie2_alignment(sample_info, mode, progress, task_id):
    """
    Run bowtie2 alignment for a single sample

    Args:
        sample_info: dict with 'sample', 'r1', 'r2'
        mode: 'original' or 'masked'
        progress: Rich progress object
        task_id: Task ID for progress bar
    """
    sample = sample_info['sample']
    r1 = sample_info['r1']
    r2 = sample_info['r2']

    # Select index based on mode
    index = ORIGINAL_INDEX if mode == 'original' else MASKED_INDEX

    # Output directory
    out_dir = os.path.join(RESULTS_DIR, f"aligned_{mode}")
    os.makedirs(out_dir, exist_ok=True)

    # Output BAM file
    bam_file = os.path.join(out_dir, f"{sample}.bam")

    progress.update(task_id, description=f"[cyan]{mode}: {sample}...[/cyan]")

    # Bowtie2 command with default parameters (local mode only)
    bowtie2_cmd = [
        "bowtie2",
        "-x", index,
        "-1", r1,
        "-2", r2,
        "--local",  # Use default scoring
        "-p", "8",  # 8 threads
        "--no-unal"  # Don't output unaligned reads
    ]

    # Pipe to samtools for sorting
    samtools_view_cmd = ["samtools", "view", "-bS", "-"]
    samtools_sort_cmd = ["samtools", "sort", "-@", "4", "-o", bam_file, "-"]

    try:
        # Run bowtie2
        bowtie2_proc = subprocess.Popen(
            bowtie2_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Pipe to samtools view
        view_proc = subprocess.Popen(
            samtools_view_cmd,
            stdin=bowtie2_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        # Pipe to samtools sort
        sort_proc = subprocess.Popen(
            samtools_sort_cmd,
            stdin=view_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        # Close upstream pipes
        bowtie2_proc.stdout.close()
        view_proc.stdout.close()

        # Wait for completion
        sort_proc.communicate()
        bowtie2_stderr = bowtie2_proc.stderr.read()

        if sort_proc.returncode != 0:
            console.print(f"[red]✗ Error sorting {sample} ({mode})[/red]")
            return None

        # Index BAM file
        index_cmd = ["samtools", "index", bam_file]
        subprocess.run(index_cmd, check=True, capture_output=True)

        # Parse alignment rate from bowtie2 stderr
        alignment_rate = parse_bowtie2_stats(bowtie2_stderr)

        progress.update(task_id, advance=1)

        return {
            'sample': sample,
            'mode': mode,
            'bam_file': bam_file,
            'alignment_rate': alignment_rate
        }

    except Exception as e:
        console.print(f"[red]✗ Error processing {sample} ({mode}): {e}[/red]")
        return None

def parse_bowtie2_stats(stderr_text):
    """Extract alignment rate from bowtie2 stderr"""
    for line in stderr_text.split('\n'):
        if '% overall alignment rate' in line:
            rate_str = line.strip().split('%')[0].strip()
            try:
                return float(rate_str)
            except ValueError:
                return None
    return None

def main():
    console.print("\n[bold cyan]═══ Step 2: Bowtie2 Alignment (Dual-Mode) ═══[/bold cyan]\n")

    # Find FASTQ pairs
    console.print("[yellow]Finding FASTQ files...[/yellow]")
    pairs = find_fastq_pairs()
    console.print(f"  [green]✓[/green] Found {len(pairs)} samples\n")

    if len(pairs) == 0:
        console.print("[red]✗ No FASTQ files found![/red]")
        return

    # Run alignments for both modes
    results = []

    for mode in ['original', 'masked']:
        console.print(f"\n[bold yellow]Aligning to {mode} reference...[/bold yellow]\n")

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TextColumn("({task.completed}/{task.total})"),
            TimeElapsedColumn(),
            console=console
        ) as progress:
            task = progress.add_task(f"[cyan]Aligning ({mode})...", total=len(pairs))

            for pair in pairs:
                result = run_bowtie2_alignment(pair, mode, progress, task)
                if result:
                    results.append(result)

        console.print(f"\n[green]✓ Completed {mode} alignments[/green]\n")

    # Create summary
    console.print("\n[yellow]Creating alignment summary...[/yellow]")
    df = pd.DataFrame(results)

    # Pivot to compare original vs masked
    if len(df) > 0:
        summary = df.pivot_table(
            index='sample',
            columns='mode',
            values='alignment_rate',
            aggfunc='first'
        ).reset_index()

        summary.columns.name = None
        summary.to_csv(os.path.join(RESULTS_DIR, "alignment_summary.csv"), index=False)

        console.print(f"  [green]✓[/green] Summary saved to alignment_summary.csv\n")

        # Print statistics
        console.print("[bold]Alignment Statistics:[/bold]")
        for mode in ['original', 'masked']:
            mode_df = df[df['mode'] == mode]
            if len(mode_df) > 0:
                rates = mode_df['alignment_rate'].dropna()
                console.print(f"\n  [cyan]{mode.capitalize()} Reference:[/cyan]")
                console.print(f"    Mean: {rates.mean():.1f}%")
                console.print(f"    Median: {rates.median():.1f}%")
                console.print(f"    Range: {rates.min():.1f}% - {rates.max():.1f}%")

    console.print(f"\n[bold green]✓ Dual-mode alignment complete![/bold green]")
    console.print(f"\n[dim]Output:[/dim]")
    console.print(f"[dim]  - {len(pairs)} BAM files in aligned_original/[/dim]")
    console.print(f"[dim]  - {len(pairs)} BAM files in aligned_masked/[/dim]")
    console.print(f"[dim]  - {len(pairs)*2} BAM index files (.bai)[/dim]")
    console.print(f"[dim]  - 1 alignment summary CSV[/dim]\n")

if __name__ == "__main__":
    main()
