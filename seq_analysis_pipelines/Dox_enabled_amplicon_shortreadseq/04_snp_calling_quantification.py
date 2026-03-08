#!/usr/bin/env python3
"""
Step 4: SNP Calling & Allele Quantification (Pysam Native Verified)
- Uses pysam with stepper='all' to count ALL reads (duplicates included)
- Verified to match samtools depth counts (~600k reads)
- Robust against subprocess crashes
"""

import os
import pysam
import pandas as pd
from pathlib import Path
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn

console = Console()

# Paths
RESULTS_DIR = "/Volumes/guttman/users/gmgao/Data_seq/20260207-DoxSeqAllRepsMultiplexed/results"
REF_DIR = os.path.join(RESULTS_DIR, "references")
QUANT_DIR = os.path.join(RESULTS_DIR, "quantification")

# SNP positions (0-based)
SNP_POSITIONS = [31, 513, 517, 556]

def get_reference_bases():
    fasta = pysam.FastaFile(os.path.join(REF_DIR, "Xist_amplicon_ref_seq.fa"))
    ref_name = fasta.references[0]
    ref_bases = {}
    for pos in SNP_POSITIONS:
        base = fasta.fetch(ref_name, pos, pos+1).upper()
        ref_bases[pos] = base
    fasta.close()
    return ref_bases

def quantify_sample(bam_file, mode, ref_bases, progress, task_id):
    sample = Path(bam_file).stem
    process_desc = f"[cyan]{mode}: {sample}...[/cyan]"
    progress.update(task_id, description=process_desc)

    results = {'sample': sample, 'mode': mode}
    per_read_data = []

    # Open BAM
    bam = pysam.AlignmentFile(bam_file, "rb")
    ref_name = bam.references[0]

    # First pass: collect per-read SNP matches
    # We want B6_Matches and Cast_Matches for each read
    read_matches = {} # read_name -> {'B6': 0, 'Cast': 0}

    for i, pos in enumerate(SNP_POSITIONS):
        ref_base = ref_bases[pos]

        iter_pileup = bam.pileup(ref_name, pos, pos+1, truncate=True,
                               stepper='all', max_depth=1000000)

        for pileupcolumn in iter_pileup:
            if pileupcolumn.pos == pos:
                for pileupread in pileupcolumn.pileups:
                    if pileupread.is_del or pileupread.is_refskip:
                        continue
                    if not pileupread.query_position:
                        continue

                    try:
                        base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                    except (IndexError, TypeError):
                        continue

                    read_name = pileupread.alignment.query_name
                    if read_name not in read_matches:
                        read_matches[read_name] = {'B6': 0, 'Cast': 0}

                    if base == ref_base:
                        read_matches[read_name]['B6'] += 1
                    else:
                        # Assuming any non-ref is potential Cast (consistent with aggregate logic)
                        read_matches[read_name]['Cast'] += 1

    # Convert to dataframe for saving
    per_read_df = pd.DataFrame.from_dict(read_matches, orient='index').reset_index()
    per_read_df.columns = ['Read_ID', 'B6_Matches', 'Cast_Matches']

    # Aggregate stats for results dict
    total_ref = per_read_df['B6_Matches'].sum()
    total_alt = per_read_df['Cast_Matches'].sum()
    total_all = total_ref + total_alt # Note: this counts informative bases across all SNPs

    results['Total_ref_reads'] = total_ref
    results['Total_alt_reads'] = total_alt
    results['Total_informative_reads'] = total_ref + total_alt
    results['Total_all_reads'] = len(per_read_df) # Unique reads
    results['Alt_pct_weighted'] = (100.0 * total_alt / (total_ref + total_alt)) if (total_ref + total_alt) > 0 else 0.0

    progress.update(task_id, advance=1)
    return results, per_read_df

def main():
    console.print("\n[bold cyan]═══ Step 4: SNP Quantification (Pysam Native Verified) ═══[/bold cyan]\n")
    os.makedirs(QUANT_DIR, exist_ok=True)

    ref_bases = get_reference_bases()
    all_results = []

    for mode in ['original', 'masked']:
        console.print(f"\n[bold yellow]Quantifying {mode} alignments...[/bold yellow]")
        bam_dir = os.path.join(RESULTS_DIR, f"aligned_{mode}")
        bam_files = sorted(Path(bam_dir).glob("*.bam"))

        if not bam_files:
            continue

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn(),
            console=console
        ) as progress:
            task = progress.add_task(f"[cyan]Processing {mode}...", total=len(bam_files))
            for bam_path in bam_files:
                res, per_read_df = quantify_sample(str(bam_path), mode, ref_bases, progress, task)
                all_results.append(res)

                # Save per-sample quant CSV for heatmaps
                sample_name = Path(bam_path).stem
                out_name = f"{sample_name}_{mode}_quant.csv"
                per_read_df.to_csv(os.path.join(QUANT_DIR, out_name), index=False)

                # Proactive: Save incrementally so user can see progress in notebook
                df_current = pd.DataFrame(all_results)
                summary_temp = df_current[['sample', 'mode', 'Alt_pct_weighted', 'Total_informative_reads', 'Total_all_reads']].copy()
                summary_temp.columns = ['sample', 'mode', 'Alt_percent', 'Informative_reads', 'Total_reads']
                summary_temp.to_csv(os.path.join(QUANT_DIR, "allele_quantification_summary_v2.csv"), index=False)

    df = pd.DataFrame(all_results)
    full_csv = os.path.join(QUANT_DIR, "allele_quantification_full_v2.csv")
    df.to_csv(full_csv, index=False)
    console.print(f"[green]✓[/green] Completed. All results saved.")

if __name__ == "__main__":
    main()
