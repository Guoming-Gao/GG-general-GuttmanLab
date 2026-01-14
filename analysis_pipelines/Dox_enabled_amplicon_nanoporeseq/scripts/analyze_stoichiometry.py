
import pandas as pd
import os
import glob
import json
from collections import Counter
import itertools

QUANT_DIR = "quantification"
QUANT_FILES = glob.glob(os.path.join(QUANT_DIR, "*_quant.csv"))
OUTPUT_DIR = "stoichiometry"
SNP_FILE = "ref_seq/snps.json"

# Load SNPs dynamically
if os.path.exists(SNP_FILE):
    with open(SNP_FILE, "r") as f:
        SNP_DATA = json.load(f)
    # Extracts B6 and Cast bases for each SNP
    B6_BASES = [s["b6"].split("/")[0] for s in SNP_DATA]
    CAST_BASES = [s["cast"].split("/")[0] for s in SNP_DATA]
    NUM_SNPS = len(SNP_DATA)
else:
    print(f"Warning: {SNP_FILE} not found.")
    B6_BASES = []
    CAST_BASES = []
    NUM_SNPS = 0

def analyze_stoichiometry(file_path):
    sample_name = os.path.basename(file_path).replace("_quant.csv", "")
    df = pd.read_csv(file_path)

    results = []

    for allele_group in ["B6", "Cast"]:
        group_df = df[df["Allele"].str.contains(allele_group, na=False)]
        if group_df.empty:
            continue

        group_size = len(group_df)
        target_bases = B6_BASES if allele_group == "B6" else CAST_BASES

        # 1. SNP Hit Distribution
        def count_hits(row):
            hits = 0
            for i in range(NUM_SNPS):
                if row[f"SNP{i+1}"] == target_bases[i]:
                    hits += 1
            return hits

        group_df = group_df.copy()
        group_df["SNP_Hits"] = group_df.apply(count_hits, axis=1)
        hit_dist = Counter(group_df["SNP_Hits"])

        # 2. Co-occurrence frequencies
        co_occ = {}
        # Individual bits
        for i in range(NUM_SNPS):
            count = len(group_df[group_df[f"SNP{i+1}"] == target_bases[i]])
            co_occ[f"SNP{i+1}_only"] = count / group_size

        # Pairs (top 3 pairs for brevity if many SNPs, but let's do all if few)
        if NUM_SNPS >= 2:
            pairs = list(itertools.combinations(range(NUM_SNPS), 2))
            for p in pairs:
                idx1, idx2 = p
                mask = (group_df[f"SNP{idx1+1}"] == target_bases[idx1]) & (group_df[f"SNP{idx2+1}"] == target_bases[idx2])
                count = len(group_df[mask])
                co_occ[f"SNP{idx1+1}_&_SNP{idx2+1}"] = count / group_size

        # All
        mask_all = True
        for i in range(NUM_SNPS):
            mask_all &= (group_df[f"SNP{i+1}"] == target_bases[i])
        count_all = len(group_df[mask_all])
        co_occ["ALL_SNPS"] = count_all / group_size

        res = {
            "Sample": sample_name,
            "Allele": allele_group,
            "Total_Reads": group_size,
        }
        for hit in range(1, NUM_SNPS + 1):
            res[f"{hit}_SNP_Hits_%"] = (hit_dist.get(hit, 0) / group_size) * 100
        res.update(co_occ)
        results.append(res)

    return results

def main():
    if NUM_SNPS == 0:
        print("No SNPs to analyze.")
        return

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    all_results = []

    for f in QUANT_FILES:
        print(f"Analyzing stoichiometry for {f}...")
        all_results.extend(analyze_stoichiometry(f))

    summary_df = pd.DataFrame(all_results)
    summary_df.to_csv(os.path.join(OUTPUT_DIR, "stoichiometry_summary.csv"), index=False)
    print("\nStoichiometry & Co-occurrence Summary:")
    print(summary_df.to_string())

if __name__ == "__main__":
    main()
