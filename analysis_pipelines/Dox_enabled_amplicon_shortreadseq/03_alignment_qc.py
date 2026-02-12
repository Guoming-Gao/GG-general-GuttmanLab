#!/usr/bin/env python3
"""
Step 3: Alignment QC with samtools
- Run samtools flagstat on all BAM files
- Run samtools stats on all BAM files
- Generate summary CSVs for original and masked alignments
"""

import os
import subprocess
import pandas as pd
from pathlib import Path
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn

console = Console()

# Paths
RESULTS_DIR = "/Volumes/guttman/users/gmgao/Data_seq/20260207-DoxSeqAllRepsMultiplexed/results"
QC_DIR = os.path.join(RESULTS_DIR, "qc", "samtools")

def find_bam_files(mode):
    """Find all BAM files for a given mode (original or masked)"""
    bam_dir = os.path.join(RESULTS_DIR, f"aligned_{mode}")
    bam_files = sorted(Path(bam_dir).glob("*.bam"))
    return [str(f) for f in bam_files]

def run_flagstat(bam_file, mode, progress, task_id):
    """Run samtools flagstat on a BAM file"""
    sample = Path(bam_file).stem
    out_dir = os.path.join(QC_DIR, mode)
    os.makedirs(out_dir, exist_ok=True)

    flagstat_file = os.path.join(out_dir, f"{sample}_flagstat.txt")

    progress.update(task_id, description=f"[cyan]flagstat: {sample} ({mode})...[/cyan]")

    cmd = ["samtools", "flagstat", bam_file]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        with open(flagstat_file, 'w') as f:
            f.write(result.stdout)

        progress.update(task_id, advance=1)
        return flagstat_file
    except subprocess.CalledProcessError as e:
        console.print(f"[red]✗ Error running flagstat on {sample}[/red]")
        return None

def run_stats(bam_file, mode, progress, task_id):
    """Run samtools stats on a BAM file"""
    sample = Path(bam_file).stem
    out_dir = os.path.join(QC_DIR, mode)
    os.makedirs(out_dir, exist_ok=True)

    stats_file = os.path.join(out_dir, f"{sample}_stats.txt")

    progress.update(task_id, description=f"[cyan]stats: {sample} ({mode})...[/cyan]")

    cmd = ["samtools", "stats", bam_file]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        with open(stats_file, 'w') as f:
            f.write(result.stdout)

        progress.update(task_id, advance=1)
        return stats_file
    except subprocess.CalledProcessError as e:
        console.print(f"[red]✗ Error running stats on {sample}[/red]")
        return None

def parse_flagstat(flagstat_file):
    """Parse samtools flagstat output"""
    with open(flagstat_file) as f:
        lines = f.readlines()

    sample = Path(flagstat_file).stem.replace("_flagstat", "")

    data = {'sample': sample}

    for line in lines:
        if 'in total' in line:
            data['total_reads'] = int(line.split()[0])
        elif 'mapped (' in line and 'primary' not in line:
            data['mapped_reads'] = int(line.split()[0])
            # Extract percentage
            pct = line.split('(')[1].split('%')[0]
            data['mapped_pct'] = float(pct)
        elif 'properly paired' in line:
            data['properly_paired'] = int(line.split()[0])
            pct = line.split('(')[1].split('%')[0]
            data['properly_paired_pct'] = float(pct)
        elif 'with itself and mate mapped' in line:
            data['both_mapped'] = int(line.split()[0])
        elif 'singletons' in line:
            data['singletons'] = int(line.split()[0])
            pct = line.split('(')[1].split('%')[0]
            data['singletons_pct'] = float(pct)

    return data

def parse_stats(stats_file):
    """Parse samtools stats output to get coverage info"""
    with open(stats_file) as f:
        lines = f.readlines()

    sample = Path(stats_file).stem.replace("_stats", "")

    data = {'sample': sample}

    for line in lines:
        if line.startswith('SN') and 'average length:' in line:
            data['avg_read_length'] = float(line.split('\t')[2].strip())
        elif line.startswith('SN') and 'average quality:' in line:
            data['avg_quality'] = float(line.split('\t')[2].strip())
        elif line.startswith('SN') and 'insert size average:' in line:
            data['insert_size_avg'] = float(line.split('\t')[2].strip())
        elif line.startswith('SN') and 'insert size standard deviation:' in line:
            data['insert_size_std'] = float(line.split('\t')[2].strip())
        elif line.startswith('SN') and 'bases mapped:' in line:
            data['bases_mapped'] = int(line.split('\t')[2].strip())
        elif line.startswith('SN') and 'reads mapped:' in line:
            data['reads_mapped_stats'] = int(line.split('\t')[2].strip())

    return data

def main():
    console.print("\n[bold cyan]═══ Step 3: Alignment QC (samtools) ═══[/bold cyan]\n")

    for mode in ['original', 'masked']:
        console.print(f"\n[bold yellow]Processing {mode} alignments...[/bold yellow]\n")

        # Find BAM files
        bam_files = find_bam_files(mode)
        console.print(f"  [green]✓[/green] Found {len(bam_files)} BAM files\n")

        if len(bam_files) == 0:
            console.print(f"[red]✗ No BAM files found for {mode}![/red]")
            continue

        # Run flagstat
        console.print(f"[yellow]Running samtools flagstat...[/yellow]")
        flagstat_files = []
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=console
        ) as progress:
            task = progress.add_task(f"[cyan]flagstat ({mode})...", total=len(bam_files))

            for bam in bam_files:
                flagstat_file = run_flagstat(bam, mode, progress, task)
                if flagstat_file:
                    flagstat_files.append(flagstat_file)

        console.print(f"[green]✓ Completed flagstat for {len(flagstat_files)} files[/green]\n")

        # Run stats
        console.print(f"[yellow]Running samtools stats...[/yellow]")
        stats_files = []
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=console
        ) as progress:
            task = progress.add_task(f"[cyan]stats ({mode})...", total=len(bam_files))

            for bam in bam_files:
                stats_file = run_stats(bam, mode, progress, task)
                if stats_file:
                    stats_files.append(stats_file)

        console.print(f"[green]✓ Completed stats for {len(stats_files)} files[/green]\n")

        # Parse and create summary
        console.print(f"[yellow]Creating summary CSV...[/yellow]")

        flagstat_data = [parse_flagstat(f) for f in flagstat_files]
        stats_data = {d['sample']: d for d in [parse_stats(f) for f in stats_files]}

        # Merge flagstat and stats data
        for fs_data in flagstat_data:
            sample = fs_data['sample']
            if sample in stats_data:
                fs_data.update({k: v for k, v in stats_data[sample].items() if k != 'sample'})

        df = pd.DataFrame(flagstat_data)

        # Calculate mean coverage (bases mapped / amplicon length)
        amplicon_length = 584
        df['mean_coverage'] = df['bases_mapped'] / amplicon_length

        summary_file = os.path.join(QC_DIR, mode, f"alignment_summary_{mode}.csv")
        df.to_csv(summary_file, index=False)

        console.print(f"  [green]✓[/green] Summary saved to alignment_summary_{mode}.csv\n")

        # Print statistics
        console.print(f"[bold]{mode.capitalize()} Reference Statistics:[/bold]")
        console.print(f"  Mean mapped reads: {df['mapped_reads'].mean():,.0f}")
        console.print(f"  Mean coverage: {df['mean_coverage'].mean():.1f}x")
        console.print(f"  Mean insert size: {df['insert_size_avg'].mean():.1f} bp")
        console.print(f"  Mean properly paired: {df['properly_paired_pct'].mean():.1f}%\n")

    console.print(f"[bold green]✓ Alignment QC complete![/bold green]")
    console.print(f"\n[dim]Output:[/dim]")
    console.print(f"[dim]  - 66 flagstat files[/dim]")
    console.print(f"[dim]  - 66 stats files[/dim]")
    console.print(f"[dim]  - 2 summary CSV files[/dim]\n")

if __name__ == "__main__":
    main()
