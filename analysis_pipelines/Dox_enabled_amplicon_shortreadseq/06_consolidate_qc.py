#!/usr/bin/env python3
"""
Step 6: QC Consolidation
- Run MultiQC to aggregate fastp and samtools reports
- Generate comprehensive QC report
"""

import os
import subprocess
from rich.console import Console

console = Console()

# Paths
RESULTS_DIR = "/Volumes/guttman/users/gmgao/Data_seq/20260207-DoxSeqAllRepsMultiplexed/results"
QC_DIR = os.path.join(RESULTS_DIR, "qc")
MULTIQC_DIR = os.path.join(QC_DIR, "multiqc")

def main():
    console.print("\n[bold cyan]═══ Step 6: QC Consolidation (MultiQC) ═══[/bold cyan]\n")

    os.makedirs(MULTIQC_DIR, exist_ok=True)

    # Input directories
    fastp_dir = os.path.join(QC_DIR, "fastp")
    samtools_dir = os.path.join(QC_DIR, "samtools")

    console.print(f"[yellow]Aggregating reports from:[/yellow]")
    console.print(f"  - {fastp_dir}")
    console.print(f"  - {samtools_dir}")

    # MultiQC command
    cmd = [
        "multiqc",
        fastp_dir,
        samtools_dir,
        "--outdir", MULTIQC_DIR,
        "--filename", "multiqc_report.html",
        "--force"  # Overwrite
    ]

    console.print("\n[yellow]Running MultiQC...[/yellow]")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            console.print(f"\n[green]✓ MultiQC complete![/green]")
            console.print(f"[dim]Output: {os.path.join(MULTIQC_DIR, 'multiqc_report.html')}[/dim]\n")

            # Print parsed stats if available in stdout
            for line in result.stdout.splitlines():
                if "report" in line.lower() or "search" in line.lower():
                    console.print(f"  {line.strip()}")

        else:
            console.print(f"\n[red]✗ MultiQC failed:[/red]")
            console.print(result.stderr)

    except FileNotFoundError:
        console.print("\n[red]✗ MultiQC not found. Please install it with 'conda install multiqc'[/red]")

if __name__ == "__main__":
    main()
