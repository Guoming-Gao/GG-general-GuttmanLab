
import pysam
import os
import glob
import pandas as pd
import numpy as np

ALIGN_DIR = "aligned"
BAM_FILES = glob.glob(os.path.join(ALIGN_DIR, "*.sorted.bam"))

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
        # pileup() result depends on reads being present
        for pileupcolumn in sam.pileup(ref_name):
            depths.append(pileupcolumn.n)

        if depths:
            # We want mean coverage across the WHOLE reference
            # padding with zeros if no reads at some positions
            full_depths = np.zeros(ref_len)
            # pileupcolumn.pos is 0-indexed position
            for i, d in enumerate(depths):
                # This is actually complex if positions are skipped
                pass

            # Simple average of what we have, but really should be sum/ref_len
            stats["Mean_Coverage"] = sum(depths) / ref_len

    return stats

def main():
    all_stats = []
    for f in BAM_FILES:
        print(f"Analyzing {f}...")
        all_stats.append(analyze_bam(f))

    df = pd.DataFrame(all_stats)
    df.to_csv("alignment_stats.csv", index=False)
    print("\nAlignment Summary:")
    print(df.to_string())

if __name__ == "__main__":
    main()
