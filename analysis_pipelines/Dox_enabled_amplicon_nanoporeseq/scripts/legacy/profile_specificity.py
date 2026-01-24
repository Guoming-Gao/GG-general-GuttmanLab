import os
import pysam
import pandas as pd
import argparse
import glob
from collections import Counter

# Constants
XIST_CHROM = "chrX"
XIST_START = 103460373
XIST_END = 103483233

def get_specificity_stats(bam_path, cat_csv, sample_name):
    """Profiles specificity (Xist vs Off-target) per read category."""

    # Load categorization
    df_cat = pd.read_csv(cat_csv)
    read_to_cat = dict(zip(df_cat["ReadID"], df_cat["Category"]))

    # Initialize stats
    # categories: Both, Head-only, Tail-only, Neither
    stats = {cat: {"Total": 0, "Mapped": 0, "On-Target": 0, "Off-Target": 0, "Unmapped": 0} for cat in ["Both", "Head-only", "Tail-only", "Neither"]}

    # Off-target tracking: category -> Counter(locus)
    off_target_loci = {cat: Counter() for cat in ["Both", "Head-only", "Tail-only", "Neither"]}

    # Track which reads we've seen in BAM to identify unmapped
    seen_reads = set()

    samfile = pysam.AlignmentFile(bam_path, "rb")
    for read in samfile.fetch(until_eof=True):
        rid = read.query_name
        if rid not in read_to_cat:
            continue

        cat = read_to_cat[rid]

        # Only process each read once (pick primary or first alignment found)
        if rid in seen_reads:
            continue
        seen_reads.add(rid)

        stats[cat]["Total"] += 1

        if read.is_unmapped:
            stats[cat]["Unmapped"] += 1
            continue

        # We focus on primary alignments for locus counting
        if read.is_secondary or read.is_supplementary:
            # We already marked it as seen/mapped from the primary alignment logic
            # If we see it here first (no primary), we count it.
            pass

        stats[cat]["Mapped"] += 1

        # Check if On-Target (Xist)
        is_xist = False
        if read.reference_name == XIST_CHROM:
            # Overlap check
            if read.reference_start < XIST_END and read.reference_end > XIST_START:
                is_xist = True

        if is_xist:
            stats[cat]["On-Target"] += 1
        else:
            stats[cat]["Off-Target"] += 1
            # Track off-target locus (roughly by chromosome and genomic bin)
            locus = f"{read.reference_name}:{read.reference_start // 1000000}Mb"
            off_target_loci[cat][locus] += 1

    samfile.close()

    # Identify unmapped reads that never appeared in BAM at all (safeguard)
    for rid, cat in read_to_cat.items():
        if rid not in seen_reads:
            stats[cat]["Total"] += 1
            stats[cat]["Unmapped"] += 1

    # Format per-category summary
    summary_rows = []
    for cat, data in stats.items():
        # Get Top 5 off-target
        top_5 = off_target_loci[cat].most_common(5)
        top_5_str = "; ".join([f"{l} ({c})" for l, c in top_5])

        row = {
            "Sample": sample_name,
            "Category": cat,
            "Total": data["Total"],
            "Mapped": data["Mapped"],
            "On-Target_Xist": data["On-Target"],
            "Off-Target": data["Off-Target"],
            "Unmapped": data["Unmapped"],
            "On-Target_Rate": (data["On-Target"] / data["Total"] * 100) if data["Total"] > 0 else 0,
            "Top5_OffTarget": top_5_str
        }
        summary_rows.append(row)

    return summary_rows

def main():
    parser = argparse.ArgumentParser(description="Profile specificity across categories and conditions.")
    parser.add_argument("--bam_dir", required=True)
    parser.add_argument("--cat_dir", required=True)
    parser.add_argument("--results_dir", required=True)
    args = parser.parse_args()

    spec_dir = os.path.join(args.results_dir, "specificity")
    os.makedirs(spec_dir, exist_ok=True)

    bam_files = glob.glob(os.path.join(args.bam_dir, "*.bam"))

    all_rows = []
    for bam in bam_files:
        sample = os.path.basename(bam).replace(".genomic.sorted.bam", "")
        cat_csv = os.path.join(args.cat_dir, f"{sample}.categorization.csv")

        if not os.path.exists(cat_csv):
            print(f"Warning: Categorization CSV not found for {sample}")
            continue

        print(f"Profiling {sample}...")
        rows = get_specificity_stats(bam, cat_csv, sample)
        all_rows.extend(rows)

    df_results = pd.DataFrame(all_rows)
    output_path = os.path.join(spec_dir, "specificity_metrics.csv")
    df_results.to_csv(output_path, index=False)

    print(f"\nSpecificity Profiling Complete. Summary saved to {output_path}")

if __name__ == "__main__":
    main()
