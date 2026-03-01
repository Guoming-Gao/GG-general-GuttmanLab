#!/usr/bin/env python3
"""
Design forward primers for RT probes using Primer3.

Forward primers are designed 200-250bp upstream of the RT probe to ensure
the PCR amplicon fully covers the 200bp RT coverage region containing SNPs.

Usage:
    python design_forward_primers.py input.csv output.csv [options]
"""

import argparse
import subprocess
import tempfile
import pandas as pd
import re
import os
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rich.console import Console
from rich.progress import track

# Initialize Rich console
console = Console()

# Try to import config for BLAST database path
try:
    sys.path.append(os.getcwd())
    from utils.config import FISH_RT_CONFIG
except ImportError:
    FISH_RT_CONFIG = {}


def extract_target_sequence(chrom: str, start: int, end: int, strand: int,
                           genome_fasta: str) -> str:
    """
    Extract sequence from genome using samtools.

    For forward primer design:
    - Plus strand genes: return plus strand sequence
    - Minus strand genes: return minus strand sequence (reverse complement)

    This ensures the forward primer extends toward the RT region.
    """
    try:
        # Find samtools
        samtools_exec = "samtools"
        try:
            subprocess.run(["samtools", "--version"], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            search_paths = [
                "/opt/miniconda3/envs/bioinfo/bin/samtools",
                "/usr/local/bin/samtools",
                "/usr/bin/samtools"
            ]
            for p in search_paths:
                if os.path.exists(p):
                    samtools_exec = p
                    break

        region = f"{chrom}:{start}-{end}"
        result = subprocess.run(
            [samtools_exec, "faidx", genome_fasta, region],
            capture_output=True, text=True, check=True
        )
        lines = result.stdout.strip().split('\n')
        seq = ''.join(lines[1:])  # Skip header

        # For minus strand genes, reverse complement the sequence
        # This gives us a template where Primer3's LEFT primer will
        # match the minus strand (bind plus strand, extend leftward)
        if strand == -1:
            complement = str.maketrans('ATCGatcg', 'TAGCtagc')
            seq = seq.translate(complement)[::-1]

        return seq.upper()
    except Exception as e:
        print(f"Warning: Failed to extract sequence for {region}: {e}")
        return ""


def run_primer3(sequence: str, target_start: int, target_length: int,
                primer_min_size: int = 18, primer_opt_size: int = 20,
                primer_max_size: int = 25, primer_min_tm: float = 57.0,
                primer_opt_tm: float = 60.0, primer_max_tm: float = 63.0,
                primer_min_gc: float = 40.0, primer_opt_gc: float = 50.0,
                primer_max_gc: float = 60.0, gc_clamp: int = 1) -> dict:
    """
    Run Primer3 to design a forward primer.

    Returns list of primer info dicts.
    """
    # Primer3 input format
    # We remove SEQUENCE_TARGET to allow picking a primer anywhere in the template
    primer3_input = f"""SEQUENCE_ID=target
SEQUENCE_TEMPLATE={sequence}
PRIMER_TASK=pick_left_only
PRIMER_PICK_LEFT_ONLY=1
PRIMER_MIN_SIZE={primer_min_size}
PRIMER_OPT_SIZE={primer_opt_size}
PRIMER_MAX_SIZE={primer_max_size}
PRIMER_MIN_TM={primer_min_tm}
PRIMER_OPT_TM={primer_opt_tm}
PRIMER_MAX_TM={primer_max_tm}
PRIMER_MIN_GC={primer_min_gc}
PRIMER_OPT_GC_PERCENT={primer_opt_gc}
PRIMER_MAX_GC={primer_max_gc}
PRIMER_GC_CLAMP={gc_clamp}
PRIMER_NUM_RETURN=20
PRIMER_MAX_POLY_X=4
PRIMER_MAX_SELF_ANY=8.0
PRIMER_MAX_SELF_END=3.0
=
"""

    # Find primer3_core
    p3_exec = "primer3_core"
    try:
        subprocess.run(["primer3_core", "--version"], capture_output=True, check=False)
    except (subprocess.CalledProcessError, FileNotFoundError):
        search_paths = [
            "/opt/miniconda3/envs/bioinfo/bin/primer3_core",
            "/usr/local/bin/primer3_core",
            "/usr/bin/primer3_core"
        ]
        for p in search_paths:
            if os.path.exists(p):
                p3_exec = p
                break

    try:
        result = subprocess.run([p3_exec], input=primer3_input, capture_output=True, text=True, check=True)
        output = result.stdout

        primers = []
        for i in range(20):
            seq_key = f'PRIMER_LEFT_{i}_SEQUENCE'
            if seq_key in output:
                # Find the value in the output block
                seq = ""
                tm = 0.0
                gc = 0.0
                pos = 0
                for line in output.strip().split('\n'):
                    if '=' not in line: continue
                    key, val = line.split('=', 1)
                    if key == f"PRIMER_LEFT_{i}_SEQUENCE": seq = val
                    if key == f"PRIMER_LEFT_{i}_TM": tm = float(val)
                    if key == f"PRIMER_LEFT_{i}_GC_PERCENT": gc = float(val)
                    if key == f"PRIMER_LEFT_{i}":
                        # Format is pos,len
                        pos = int(val.split(',')[0])

                primers.append({
                    'sequence': seq,
                    'tm': tm,
                    'gc': gc,
                    'position': pos
                })
        return primers
    except Exception as e:
        console.print(f"[yellow]Warning: Primer3 failed: {e}[/yellow]")
        return []

def is_low_complexity(sequence: str, threshold: float = 0.5) -> bool:
    """
    Check if a sequence has low complexity (repetitive).
    Uses the ratio of unique 3-mers to total possible 3-mers.
    """
    if len(sequence) < 10: return False

    # Check for simple repeats (e.g. TGTGTGTG)
    kmers = set()
    for i in range(len(sequence) - 2):
        kmers.add(sequence[i:i+3])

    # 20bp sequence has 18 possible 3-mers
    ratio = len(kmers) / (len(sequence) - 2)
    return ratio < threshold


def validate_primer_specificity(sequence: str, genome_fasta_or_db: str) -> bool:
    """
    Validation step: Run BLAST and strictly parse results.
    We want:
    1. Exactly 1 PERFECT or near-perfect hit (Length >= 18, Identity >= 95%).
    2. ZERO other hits that meet a noise threshold (Length >= 12).
    """
    if not sequence or is_low_complexity(sequence):
        return False

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        tmp.write(f">primer\n{sequence}\n")
        tmp_path = tmp.name

    try:
        # Find blastn in common locations
        blastn_exec = "blastn"
        search_paths = [
            "/opt/miniconda3/envs/bioinfo/bin/blastn",
            "/usr/local/bin/blastn",
            "/usr/bin/blastn"
        ]
        for p in search_paths:
            if os.path.exists(p):
                blastn_exec = p
                break

        cmd = [
            blastn_exec,
            "-query", tmp_path,
            "-db", genome_fasta_or_db,
            "-outfmt", "0",
            "-evalue", "10.0",  # High sensitivity to find even weak off-targets
            "-word_size", "7",   # Very sensitive for short primers
            "-max_target_seqs", "30",
            "-dust", "no"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        blast_text = result.stdout

        # Find all HSPs and their alignment stats
        # Format: "Identities = 21/29 (72%), Gaps = 0/29 (0%)"
        hsp_identities = re.findall(r"Identities\s*=\s*(\d+)/(\d+)\s*\((\d+)%\)", blast_text)

        perfect_hits = 0
        noise_hits = 0

        for matched, total, perc in hsp_identities:
            matched = int(matched)
            total = int(total)
            perc = int(perc)

            # A "Perfect" hit for a 20bp primer: at least 18bp and high identity
            if total >= 18 and perc >= 95:
                perfect_hits += 1
            # A "Noise" hit: anything significant enough to cause off-target issues
            elif total >= 12:
                noise_hits += 1

        # Specificity requirements:
        # Must have exactly 1 perfect match (the target)
        # Must have ZERO significant off-target matches
        is_specific = (perfect_hits == 1 and noise_hits == 0)

        # Specificity requirements:
        if is_specific:
            # console.print(f"[dim]Primer {sequence} is specific.[/dim]")
            pass
        else:
            # console.print(f"[dim]Primer {sequence} rejected: perfect={perfect_hits}, noise={noise_hits}[/dim]")
            pass

        return is_specific
    except Exception as e:
        print(f"Warning: Specificity check failed: {e}")
        # If specificity check fails due to environment, we'd rather be safe
        # but here we'll return False to avoid designing non-specific primers
        return False
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def design_forward_primers(input_file: str, output_file: str,
                          genome_fasta: str,
                          target_distance: int = 225,
                          primer_min_size: int = 18,
                          primer_opt_size: int = 20,
                          primer_max_size: int = 25,
                          primer_min_tm: float = 57.0,
                          primer_opt_tm: float = 60.0,
                          primer_max_tm: float = 63.0,
                          primer_min_gc: float = 40.0,
                          primer_max_gc: float = 60.0,
                          args: argparse.Namespace = None):
    """
    Design forward primers for RT probes.

    Args:
        input_file: CSV with RT probes (needs Chromosome, theStartPos, theEndPos, Target_Strand)
        output_file: Output CSV with forward primers added
        genome_fasta: Path to genome FASTA file
        target_distance: Distance upstream of RT probe to design primer (default 225bp)
    """
    df = pd.read_csv(input_file)
    console.print(f"📄 [cyan]Loaded {len(df)} probes from {input_file}[/cyan]")

    # Result columns
    forward_primers = []
    forward_tms = []
    forward_gcs = []
    pcr_sizes = []
    rt_covered = []

    # Use track for progress bar
    for idx, row in track(df.iterrows(), total=len(df), description="Designing forward primers..."):
        chrom = row['Chromosome']
        # Support both column naming conventions
        probe_start = row.get('Probe_Start', row.get('theStartPos', 0))
        probe_end = row.get('Probe_End', row.get('theEndPos', 0))
        strand_raw = row['Target_Strand']
        # Parse strand - handle both string ('-'/'+') and integer (-1/1)
        if strand_raw in [-1, '-', '-1']:
            strand = -1
        else:
            strand = 1
        rt_start = row.get('RT_Region_Start', probe_start)
        rt_end = row.get('RT_Region_End', probe_end)

        # Calculate region to search for forward primer
        # Forward primer should be upstream of RT probe.
        # We expand search AWAY from the probe to ensure the PCR product covers the RT region.
        # User specified: search region follows rt_coverage_downstream + 50
        min_dist = int(args.rt_coverage) + 50
        max_dist = int(args.max_distance)

        if strand == 1:  # Plus strand: primer upstream (lower coords)
            search_start = probe_start - max_dist
            search_end = probe_start - min_dist
        else:  # Minus strand
            # Forward primer region is downstream (larger coordinates)
            search_start = probe_end + min_dist
            search_end = probe_end + max_dist

        target_in_seq = 100 # Not strictly used with pick_left_only and open template

        # Ensure valid coordinates
        search_start = max(1, search_start)

        # Extract sequence for Primer3
        seq = extract_target_sequence(chrom, search_start, search_end, strand, genome_fasta)

        if not seq or len(seq) < 50:
            forward_primers.append('')
            forward_tms.append(0)
            forward_gcs.append(0)
            pcr_sizes.append(0)
            rt_covered.append(False)
            continue

        # Run Primer3 to get up to 20 candidates
        primers = run_primer3(
            seq, 0, len(seq), # Search across the entire template
            primer_min_size, primer_opt_size, primer_max_size,
            primer_min_tm, primer_opt_tm, primer_max_tm,
            primer_min_gc, primer_max_gc
        )

        # Validate each candidate with BLAST
        selected_primer = None
        blast_db = FISH_RT_CONFIG.get("blast_database_path", genome_fasta)

        for p in primers:
            if validate_primer_specificity(p['sequence'], blast_db):
                selected_primer = p
                break

        if selected_primer:
            forward_primers.append(selected_primer['sequence'])
            forward_tms.append(selected_primer['tm'])
            forward_gcs.append(selected_primer['gc'])
            # Calculate PCR size correctly: |ProbeSite - PrimerSite|
            # ProbeSite is probe_start (for + strand) or probe_end (for - strand)
            # Primer site is relative to search_start
            p_pos = selected_primer['position'] # 0-based index in template string
            if strand == 1:
                # Plus strand: Template is [search_start, search_end]
                primer_genomic_pos = search_start + p_pos
                pcr_size = probe_start - primer_genomic_pos
            else:
                # Minus strand: Template is RC of [search_start, search_end]
                # Index 0 is search_end
                primer_genomic_pos = search_end - p_pos
                pcr_size = primer_genomic_pos - probe_end

            pcr_sizes.append(pcr_size)
            rt_covered.append(pcr_size >= int(args.rt_coverage))
        else:
            forward_primers.append('')
            forward_tms.append(0)
            forward_gcs.append(0)
            pcr_sizes.append(0)
            rt_covered.append(False)


    # Add columns to dataframe
    df['Forward_Primer_Seq'] = forward_primers
    df['Forward_Primer_Tm'] = forward_tms
    df['Forward_Primer_GC'] = forward_gcs
    df['PCR_Amplicon_Size'] = pcr_sizes
    df['RT_Region_Covered'] = rt_covered

    # Save output
    df.to_csv(output_file, index=False)

    # Summary
    success_count = sum(1 for p in forward_primers if p)
    console.print(f"\n✅ [green]Forward primer design complete![/green]")
    console.print(f"   Successful designs: [bold]{success_count}/{len(df)}[/bold] ({100*success_count/len(df):.1f}%)")
    console.print(f"   Average Tm: {sum(forward_tms)/max(success_count,1):.1f}°C")
    console.print(f"   Average GC: {sum(forward_gcs)/max(success_count,1):.1f}%")
    console.print(f"   Output saved to: [yellow]{output_file}[/yellow]")

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Design forward primers for RT probes using Primer3"
    )
    parser.add_argument("input", help="Input CSV with RT probes")
    parser.add_argument("output", help="Output CSV with forward primers")
    parser.add_argument(
        "--genome",
        default="/Volumes/guttman/genomes/mm10/fasta/mm10.fa",
        help="Path to genome FASTA file"
    )
    parser.add_argument(
        "--distance",
        type=int, default=225,
        help="Target distance from RT probe (default: 225bp)"
    )
    parser.add_argument(
        "--max-distance",
        type=int, default=750,
        help="Maximum upstream search distance (default: 750bp)"
    )
    parser.add_argument(
        "--rt-coverage",
        type=int, default=None,
        help="RT coverage length in nt (default: config value or 200)"
    )
    parser.add_argument(
        "--primer-length-min", type=int, default=18,
        help="Minimum primer length (default: 18)"
    )
    parser.add_argument(
        "--primer-length-opt", type=int, default=20,
        help="Optimal primer length (default: 20)"
    )
    parser.add_argument(
        "--primer-length-max", type=int, default=25,
        help="Maximum primer length (default: 25)"
    )
    parser.add_argument(
        "--tm-min", type=float, default=57.0,
        help="Minimum Tm (default: 57°C)"
    )
    parser.add_argument(
        "--tm-opt", type=float, default=60.0,
        help="Optimal Tm (default: 60°C)"
    )
    parser.add_argument(
        "--tm-max", type=float, default=63.0,
        help="Maximum Tm (default: 63°C)"
    )
    parser.add_argument(
        "--gc-min", type=float, default=40.0,
        help="Minimum GC%% (default: 40)"
    )
    parser.add_argument(
        "--gc-max", type=float, default=60.0,
        help="Maximum GC%% (default: 60)"
    )

    args = parser.parse_args()

    # Load RT coverage from config if not provided
    rt_coverage = args.rt_coverage if args.rt_coverage is not None else FISH_RT_CONFIG.get("rt_coverage_downstream", 200)
    args.rt_coverage = rt_coverage  # Update args for use in functions

    design_forward_primers(
        args.input, args.output,
        args.genome,
        args.distance,
        args.primer_length_min, args.primer_length_opt, args.primer_length_max,
        args.tm_min, args.tm_opt, args.tm_max,
        args.gc_min, args.gc_max,
        args
    )


if __name__ == "__main__":
    main()
