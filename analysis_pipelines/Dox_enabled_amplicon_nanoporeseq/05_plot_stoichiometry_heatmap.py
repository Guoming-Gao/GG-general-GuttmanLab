
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
import argparse
import numpy as np

def plot_snp_heatmap(quant_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    quant_files = glob.glob(os.path.join(quant_dir, "*_quant.csv"))

    for f in quant_files:
        sample = os.path.basename(f).replace("_quant.csv", "")
        df = pd.read_csv(f)

        plt.figure(figsize=(10, 8))

        # We use hist2d to create the heatmap with a colorbar
        # range [[0, 6], [0, 6]] with 6 bins captures 0, 1, 2, 3, 4, 5 discrete counts
        h = plt.hist2d(df['B6_Matches'], df['Cast_Matches'], bins=6,
                       cmap='viridis', cmin=1, range=[[0, 6], [0, 6]])

        plt.colorbar(label='Number of Reads')

        plt.xlabel("B6 SNP Matches")
        plt.ylabel("Cast SNP Matches")
        plt.title(f"SNP Match Heatmap: {sample}")
        plt.grid(True, linestyle='--', alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{sample}_snp_heatmap.png"), dpi=300)
        plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--quant_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()
    plot_snp_heatmap(args.quant_dir, args.output_dir)
