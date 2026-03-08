#!/usr/bin/env python3
"""
Step 0: Reference Preparation
- Map amplicon to genome
- Extract SNPs from VCF
- Create N-masked reference
- Build bowtie2 indices
"""

import pysam
import os
import subprocess
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn

console = Console()

# Paths
AMPLICON_FA = "/Volumes/guttman/users/gmgao/Data_seq/20260207-DoxSeqAllRepsMultiplexed/Xist_amplicon_ref_seq.fa"
GENOME_FA = "/Volumes/guttman/genomes/mm10/fasta/mm10.fa"
VCF_FILE = "/Volumes/guttman/genomes/mm10/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
RESULTS_DIR = "/Volumes/guttman/users/gmgao/Data_seq/20260207-DoxSeqAllRepsMultiplexed/results"
REF_DIR = os.path.join(RESULTS_DIR, "references")

B6_SAMPLE = "C57BL_6NJ"
CAST_SAMPLE = "CAST_EiJ"

def main():
    console.print("\n[bold cyan]═══ Step 0: Reference Preparation ═══[/bold cyan]\n")

    # Step 1: Map amplicon to genome
    console.print("[yellow]1. Mapping amplicon to genome...[/yellow]")
    amplicon_coords = map_amplicon_to_genome()

    # Step 2: Extract SNPs
    console.print("\n[yellow]2. Extracting SNPs from VCF...[/yellow]")
    snps = extract_snps(amplicon_coords)

    # Step 3: Create N-masked reference
    console.print("\n[yellow]3. Creating N-masked reference...[/yellow]")
    create_masked_reference(snps)

    # Step 4: Build bowtie2 indices
    console.print("\n[yellow]4. Building bowtie2 indices...[/yellow]")
    build_bowtie2_indices()

    console.print("\n[bold green]✓ Reference preparation complete![/bold green]\n")

    # Print summary
    print_summary(amplicon_coords, snps)

def map_amplicon_to_genome():
    """Map amplicon sequence to mm10 genome"""
    # Read amplicon
    with open(AMPLICON_FA) as f:
        lines = f.readlines()
        amplicon_seq = ''.join(line.strip() for line in lines if not line.startswith('>')).upper()

    console.print(f"  Amplicon length: {len(amplicon_seq)}bp")

    # Search in Xist locus
    fasta = pysam.FastaFile(GENOME_FA)
    xist_start = 103460373 - 1  # 0-based
    xist_end = 103483233
    xist_seq = fasta.fetch("chrX", xist_start, xist_end).upper()

    # Forward search
    pos = xist_seq.find(amplicon_seq)
    if pos != -1:
        genomic_start = xist_start + pos + 1  # 1-based
        genomic_end = genomic_start + len(amplicon_seq) - 1
        strand = "+"
    else:
        # Reverse complement search
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        amplicon_rc = ''.join(complement.get(base, base) for base in reversed(amplicon_seq))
        pos = xist_seq.find(amplicon_rc)
        if pos != -1:
            genomic_start = xist_start + pos + 1
            genomic_end = genomic_start + len(amplicon_seq) - 1
            strand = "-"
        else:
            console.print(f"  [red]✗ Amplicon not found in Xist region[/red]")
            fasta.close()
            raise ValueError("Amplicon not found")

    fasta.close()

    console.print(f"  [green]✓[/green] Amplicon found at chrX:{genomic_start}-{genomic_end} (strand: {strand})")

    return {
        'chrom': 'chrX',
        'start': genomic_start,
        'end': genomic_end,
        'strand': strand,
        'length': len(amplicon_seq)
    }

def extract_snps(amplicon_coords):
    """Extract B6/Cast SNPs within amplicon region"""
    vcf = pysam.VariantFile(VCF_FILE)

    # Fetch SNPs (VCF uses chromosome names without 'chr')
    chrom = 'X'
    snps = []

    for record in vcf.fetch(chrom, amplicon_coords['start'], amplicon_coords['end']):
        b6_gt = record.samples[B6_SAMPLE]['GT']
        cast_gt = record.samples[CAST_SAMPLE]['GT']

        if b6_gt != cast_gt:
            b6_base = record.ref if b6_gt == (0, 0) else record.alleles[b6_gt[0]]
            cast_base = record.ref if cast_gt == (0, 0) else record.alleles[cast_gt[0]]

            if b6_base != cast_base:
                # Calculate position within amplicon (0-based)
                if amplicon_coords['strand'] == "+":
                    amplicon_pos = record.pos - amplicon_coords['start']
                else:
                    # For minus strand, pos 0 in FASTA corresponds to amplicon_coords['end']
                    # Indexing is reversed: amplicon_pos = end - genomic_pos
                    amplicon_pos = amplicon_coords['end'] - record.pos

                snps.append({
                    'genomic_pos': record.pos,
                    'amplicon_pos': amplicon_pos,
                    'b6': b6_base,
                    'cast': cast_base
                })

    vcf.close()

    console.print(f"  [green]✓[/green] Found {len(snps)} B6/Cast SNPs in amplicon")

    # Filter to accessible SNPs (outside 120-404bp gap)
    accessible_snps = [s for s in snps if s['amplicon_pos'] < 120 or s['amplicon_pos'] >= 404]
    inaccessible_snps = [s for s in snps if 120 <= s['amplicon_pos'] < 404]

    console.print(f"  [cyan]→[/cyan] {len(accessible_snps)} SNPs accessible to short reads")
    console.print(f"  [dim]→[/dim] {len(inaccessible_snps)} SNPs in coverage gap (will be ignored)")

    # Save SNP file
    snp_file = os.path.join(REF_DIR, "xist_snps.txt")
    with open(snp_file, 'w') as f:
        f.write("POS\tAMPLICON_POS\tB6_BASE\tCAST_BASE\tACCESSIBLE\n")
        for snp in snps:
            accessible = "YES" if snp in accessible_snps else "NO"
            f.write(f"{snp['genomic_pos']}\t{snp['amplicon_pos']}\t{snp['b6']}\t{snp['cast']}\t{accessible}\n")

    console.print(f"  [green]✓[/green] SNPs saved to xist_snps.txt")

    return accessible_snps

def create_masked_reference(snps):
    """Create N-masked reference with SNPs replaced by N"""
    # Read original amplicon
    with open(AMPLICON_FA) as f:
        header = f.readline().strip()
        amplicon_seq = ''.join(line.strip() for line in f).upper()

    # Create masked version
    seq_list = list(amplicon_seq)
    for snp in snps:
        seq_list[snp['amplicon_pos']] = 'N'

    masked_seq = ''.join(seq_list)
    n_count = masked_seq.count('N')

    # Copy original to ref dir
    orig_ref = os.path.join(REF_DIR, "Xist_amplicon_ref_seq.fa")
    with open(orig_ref, 'w') as f:
        f.write(header + '\n')
        # Write in 80-char lines
        for i in range(0, len(amplicon_seq), 80):
            f.write(amplicon_seq[i:i+80] + '\n')

    # Save masked version
    masked_ref = os.path.join(REF_DIR, "Xist_amplicon_ref_seq-N_masked.fa")
    with open(masked_ref, 'w') as f:
        f.write(header + f" (N-masked at {n_count} SNP positions)\n")
        # Write in 80-char lines
        for i in range(0, len(masked_seq), 80):
            f.write(masked_seq[i:i+80] + '\n')

    console.print(f"  [green]✓[/green] Original reference: Xist_amplicon_ref_seq.fa")
    console.print(f"  [green]✓[/green] N-masked reference: Xist_amplicon_ref_seq-N_masked.fa ({n_count} positions masked)")

def build_bowtie2_indices():
    """Build bowtie2 indices for both references"""
    orig_ref = os.path.join(REF_DIR, "Xist_amplicon_ref_seq.fa")
    masked_ref = os.path.join(REF_DIR, "Xist_amplicon_ref_seq-N_masked.fa")

    with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
        # Build original index
        task = progress.add_task("Building original index...", total=None)
        result = subprocess.run(
            ["bowtie2-build", orig_ref, os.path.join(REF_DIR, "xist_amplicon_original")],
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            console.print(f"[red]Error building original index:[/red] {result.stderr}")
            raise RuntimeError("bowtie2-build failed")
        progress.update(task, completed=True)
        console.print(f"  [green]✓[/green] Original index built: xist_amplicon_original")

        # Build masked index
        task = progress.add_task("Building N-masked index...", total=None)
        result = subprocess.run(
            ["bowtie2-build", masked_ref, os.path.join(REF_DIR, "xist_amplicon_masked")],
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            console.print(f"[red]Error building masked index:[/red] {result.stderr}")
            raise RuntimeError("bowtie2-build failed")
        progress.update(task, completed=True)
        console.print(f"  [green]✓[/green] N-masked index built: xist_amplicon_masked")

def print_summary(coords, snps):
    """Print summary table"""
    console.print("\n[bold]Summary:[/bold]")
    console.print(f"  Amplicon location: {coords['chrom']}:{coords['start']}-{coords['end']}")
    console.print(f"  Amplicon length: {coords['length']}bp")
    console.print(f"  Strand: {coords['strand']}")
    console.print(f"  Accessible SNPs: {len(snps)}")
    console.print(f"  Coverage gap: 120-404bp (284bp)")
    console.print(f"\n  Files created in {REF_DIR}:")
    console.print(f"    - xist_snps.txt")
    console.print(f"    - Xist_amplicon_ref_seq.fa")
    console.print(f"    - Xist_amplicon_ref_seq-N_masked.fa")
    console.print(f"    - xist_amplicon_original.*.bt2 (6 files)")
    console.print(f"    - xist_amplicon_masked.*.bt2 (6 files)")

if __name__ == "__main__":
    main()
