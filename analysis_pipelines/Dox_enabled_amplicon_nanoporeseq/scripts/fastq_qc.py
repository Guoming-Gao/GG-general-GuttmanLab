
import pysam
import os
import glob
import pandas as pd
import numpy as np

DATA_DIR = "data"
FASTQ_FILES = glob.glob(os.path.join(DATA_DIR, "*.fastq"))

def analyze_fastq(file_path):
    stats = {
        "File": os.path.basename(file_path),
        "Total_Reads": 0,
        "Total_Bases": 0,
        "Mean_Read_Length": 0,
        "Median_Read_Length": 0,
        "Mean_Quality": 0
    }

    lengths = []
    qualities = []

    # Nanopore reads can be many, let's sample or process all if small
    # Since these are amplicons, they shouldn't be too huge in terms of total reads
    with pysam.FastxFile(file_path) as fh:
        for entry in fh:
            stats["Total_Reads"] += 1
            l = len(entry.sequence)
            stats["Total_Bases"] += l
            lengths.append(l)
            # entry.get_quality_array() returns list of integers
            qualities.append(np.mean(entry.get_quality_array()))

    if stats["Total_Reads"] > 0:
        stats["Mean_Read_Length"] = np.mean(lengths)
        stats["Median_Read_Length"] = np.median(lengths)
        stats["Mean_Quality"] = np.mean(qualities)

    return stats

def main():
    all_stats = []
    for f in FASTQ_FILES:
        print(f"Analyzing {f}...")
        all_stats.append(analyze_fastq(f))

    df = pd.DataFrame(all_stats)
    df.to_csv("preprocessing_stats.csv", index=False)
    print("\nQC Summary:")
    print(df.to_string())

if __name__ == "__main__":
    main()
