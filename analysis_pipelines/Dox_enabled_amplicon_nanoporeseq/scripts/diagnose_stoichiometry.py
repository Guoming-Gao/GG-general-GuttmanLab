import pysam
import os
import glob
import pandas as pd
import json
import argparse
from collections import Counter

def get_base(read, pos_0):
    for q_pos, r_pos in read.get_aligned_pairs():
        if r_pos == pos_0:
            if q_pos is None: return "-"
            return read.query_sequence[q_pos]
    return "N"

def diagnostic(bam_path, snps):
    sample = os.path.basename(bam_path)
    print(f"\n--- Diagnostic for {sample} ---")

    # Track hits for ALL primary reads
    all_hits_b6 = []
    all_hits_cast = []

    with pysam.AlignmentFile(bam_path, "rb") as sam:
        count = 0
        for read in sam.fetch("Xist_Amplicon"):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue

            bases = [get_base(read, s["local_pos"] - 1) for s in snps]
            b6_vals = [s["b6"].split("/")[0] for s in snps]
            cast_vals = [s["cast"].split("/")[0] for s in snps]

            b6_hits = sum(1 for i, b in enumerate(bases) if b == b6_vals[i])
            cast_hits = sum(1 for i, b in enumerate(bases) if b == cast_vals[i])

            all_hits_b6.append(b6_hits)
            all_hits_cast.append(cast_hits)

            if count < 5: # Print first 5 reads
                print(f"Read: {read.query_name[:10]} | Strand: {'-' if read.is_reverse else '+'} | Bases: {''.join(bases)} | B6 Hits: {b6_hits} | Cast Hits: {cast_hits}")
            count += 1

    # Frequency of hits
    print("\nB6 Hit Distribution (Sample-wide):")
    print(dict(sorted(Counter(all_hits_b6).items())))
    print("\nCast Hit Distribution (Sample-wide):")
    print(dict(sorted(Counter(all_hits_cast).items())))

def main():
    parser = argparse.ArgumentParser(description="Run diagnostics on stoichiometry.")
    parser.add_argument("--results_dir", default=None, help="Results directory (default: ./results)")
    parser.add_argument("--snp_file", default=None, help="Path to SNP JSON (default: results/ref_seq/snps.json)")
    args = parser.parse_args()

    results_dir = args.results_dir if args.results_dir else os.path.join(".", "results")
    align_dir = os.path.join(results_dir, "aligned")
    snp_file = args.snp_file if args.snp_file else os.path.join(results_dir, "ref_seq", "snps.json")

    if not os.path.exists(snp_file):
        print(f"Error: {snp_file} not found.")
        return

    with open(snp_file, "r") as f:
        snps = json.load(f)

    bam_files = glob.glob(os.path.join(align_dir, "*.sorted.bam"))
    if not bam_files:
        print(f"No sorted BAM files found in {align_dir}")
        return

    # Test on WT (diff) and WT (Dox)
    for b in bam_files[:2]:
        diagnostic(b, snps)

if __name__ == "__main__":
    main()
