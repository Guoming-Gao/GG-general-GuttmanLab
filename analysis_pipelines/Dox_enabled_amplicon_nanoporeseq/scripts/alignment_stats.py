import pysam
import os
import glob
import pandas as pd
import numpy as np
import argparse

def analyze_bam(file_path):
    stats = {
        "Sample": os.path.basename(file_path).replace(".sorted.bam", ""),
        "Total_Records": 0,
        "Primary_Alignments": 0,
        "Secondary_Alignments": 0,
        "Supplementary_Alignments": 0,
        "Unmapped_Reads": 0,
        "Mean_MQ": 0,
        "Median_MQ": 0,
        "Mean_Coverage": 0
    }

    mqs = []

    with pysam.AlignmentFile(file_path, "rb") as sam:
        # Check reference name
        if not sam.references:
            return stats

        ref_name = sam.references[0] # Xist_Amplicon
        ref_len = sam.lengths[0]

        for read in sam.fetch(until_eof=True):
            stats["Total_Records"] += 1
            if read.is_unmapped:
                stats["Unmapped_Reads"] += 1
            elif read.is_secondary:
                stats["Secondary_Alignments"] += 1
            elif read.is_supplementary:
                stats["Supplementary_Alignments"] += 1
            else:
                stats["Primary_Alignments"] += 1
                mqs.append(read.mapping_quality)

        if mqs:
            stats["Mean_MQ"] = np.mean(mqs)
            stats["Median_MQ"] = np.median(mqs)

        # Coverage
        depths = []
        for pileupcolumn in sam.pileup(ref_name):
            depths.append(pileupcolumn.n)

        if depths:
            stats["Mean_Coverage"] = sum(depths) / ref_len

    return stats

def main():
    parser = argparse.ArgumentParser(description="Analyze alignment statistics from BAM files.")
    parser.add_argument("--output_dir", default=".", help="Project output directory (default: current directory)")
    args = parser.parse_args()

    results_dir = os.path.join(args.output_dir, "results")
    align_dir = os.path.join(results_dir, "aligned")
    csv_dir = os.path.join(results_dir, "csv")

    bam_files = glob.glob(os.path.join(align_dir, "*.sorted.bam"))

    if not bam_files:
        print(f"No sorted BAM files found in {align_dir}")
        return

    all_stats = []
    for f in bam_files:
        print(f"Analyzing {f}...")
        all_stats.append(analyze_bam(f))

    df = pd.DataFrame(all_stats)
    os.makedirs(csv_dir, exist_ok=True)
    df.to_csv(os.path.join(csv_dir, "alignment_stats.csv"), index=False)
    print("\nAlignment Summary:")
    print(df.to_string())

if __name__ == "__main__":
    main()
