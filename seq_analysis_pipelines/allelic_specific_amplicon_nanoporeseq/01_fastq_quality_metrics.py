
import os
import argparse
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def calculate_n50(lengths):
    lengths = sorted(lengths, reverse=True)
    total_len = sum(lengths)
    acc_len = 0
    for l in lengths:
        acc_len += l
        if acc_len >= total_len / 2:
            return l
    return 0

def run_qc(fastq_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    qc_data = []

    fastq_files = [f for f in os.listdir(fastq_dir) if f.endswith(".fastq") or f.endswith(".fastq.gz")]

    for f in fastq_files:
        path = os.path.join(fastq_dir, f)
        sample = f.split(".")[0]
        lengths = []
        qualities = []

        with pysam.FastxFile(path) as fx:
            for entry in fx:
                lengths.append(len(entry.sequence))
                # Average quality
                if entry.quality:
                    quals = [ord(q) - 33 for q in entry.quality]
                    qualities.append(np.mean(quals))

        if not lengths:
            continue

        n50 = calculate_n50(lengths)
        mean_len = np.mean(lengths)
        median_len = np.median(lengths)
        min_len = np.min(lengths)
        max_len = np.max(lengths)
        mean_qual = np.mean(qualities) if qualities else 0

        qc_data.append({
            "Sample": sample,
            "TotalReads": len(lengths),
            "N50": n50,
            "MeanLength": mean_len,
            "MedianLength": median_len,
            "MinLength": min_len,
            "MaxLength": max_len,
            "MeanQuality": mean_qual
        })

        # Plot length distribution
        plt.figure(figsize=(10, 6))
        plt.hist(lengths, bins=50, color='skyblue', edgecolor='black')
        plt.title(f"Length Distribution: {sample}")
        plt.xlabel("Length (bp)")
        plt.ylabel("Count")
        plt.savefig(os.path.join(output_dir, f"{sample}_length_dist.png"))
        plt.close()

    df = pd.DataFrame(qc_data)
    df.to_csv(os.path.join(output_dir, "qc_summary.csv"), index=False)
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", required=True)
    parser.add_argument("--results_dir", required=True)
    args = parser.parse_args()

    output_dir = os.path.join(args.results_dir, "qc")
    print(f"Saving QC results to: {output_dir}")
    run_qc(args.data_dir, output_dir)
