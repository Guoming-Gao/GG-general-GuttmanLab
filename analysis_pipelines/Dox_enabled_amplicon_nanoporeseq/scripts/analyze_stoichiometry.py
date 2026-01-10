
import pandas as pd
import os
import glob
from collections import Counter
import itertools

QUANT_DIR = "quantification"
QUANT_FILES = glob.glob(os.path.join(QUANT_DIR, "*_quant.csv"))
OUTPUT_DIR = "stoichiometry"

# SNP Positions (from quantify_alleles.py)
SNPS = [
    {"pos": 67, "b6": "T", "cast": "C", "label": "SNP1_67"},
    {"pos": 334, "b6": "G", "cast": "A", "label": "SNP2_334"},
    {"pos": 553, "b6": "G", "cast": "A", "label": "SNP3_553"}
]

def analyze_stoichiometry(file_path):
    sample_name = os.path.basename(file_path).replace("_quant.csv", "")
    df = pd.read_csv(file_path)

    results = []

    # Analyze B6 and Cast groups separately
    for allele_group in ["B6", "Cast"]:
        # Include low_conf for stoichiometry depth
        group_df = df[df["Allele"].str.contains(allele_group, na=False)]

        if group_df.empty:
            continue

        group_size = len(group_df)

        # 1. SNP Hit Distribution (1, 2, 3 hits)
        def count_hits(row):
            hits = 0
            if allele_group == "B6":
                if row["SNP1"] == "T": hits += 1
                if row["SNP2"] == "G": hits += 1
                if row["SNP3"] == "G": hits += 1
            else:
                if row["SNP1"] == "C": hits += 1
                if row["SNP2"] == "A": hits += 1
                if row["SNP3"] == "A": hits += 1
            return hits

        group_df = group_df.copy()
        group_df["SNP_Hits"] = group_df.apply(count_hits, axis=1)
        hit_dist = Counter(group_df["SNP_Hits"])

        # 2. Co-occurrence frequencies
        # combinations of SNP1, SNP2, SNP3
        co_occ = {}
        target_bases = ["T", "G", "G"] if allele_group == "B6" else ["C", "A", "A"]

        # Individual
        for i in range(3):
            hit_col = f"SNP{i+1}"
            count = len(group_df[group_df[hit_col] == target_bases[i]])
            co_occ[f"SNP{i+1}_only"] = count / group_size

        # Pairs
        pairs = list(itertools.combinations(range(3), 2))
        for p in pairs:
            idx1, idx2 = p
            mask = (group_df[f"SNP{idx1+1}"] == target_bases[idx1]) & (group_df[f"SNP{idx2+1}"] == target_bases[idx2])
            count = len(group_df[mask])
            co_occ[f"SNP{idx1+1}_&_SNP{idx2+1}"] = count / group_size

        # Triplet
        mask_all = (group_df["SNP1"] == target_bases[0]) & (group_df["SNP2"] == target_bases[1]) & (group_df["SNP3"] == target_bases[2])
        count_all = len(group_df[mask_all])
        co_occ["SNP1_&_SNP2_&_SNP3"] = count_all / group_size

        results.append({
            "Sample": sample_name,
            "Allele": allele_group,
            "Total_Reads": group_size,
            "1_SNP_Hit_%": (hit_dist.get(1, 0) / group_size) * 100,
            "2_SNP_Hits_%": (hit_dist.get(2, 0) / group_size) * 100,
            "3_SNP_Hits_%": (hit_dist.get(3, 0) / group_size) * 100,
            **co_occ
        })

    return results

def main():
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
