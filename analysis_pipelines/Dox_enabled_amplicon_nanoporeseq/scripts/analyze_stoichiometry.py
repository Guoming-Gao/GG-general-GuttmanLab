import pandas as pd
import os
import glob
import json
import argparse
from collections import Counter
import itertools

def analyze_stoichiometry(file_path, num_snps, target_bases):
    sample_name = os.path.basename(file_path).replace("_quant.csv", "")
    df = pd.read_csv(file_path)

    results = []

    for allele_group in ["B6", "Cast"]:
        group_df = df[df["Allele"].str.contains(allele_group, na=False)]
        if group_df.empty:
            continue

        group_size = len(group_df)
        group_target_bases = target_bases[allele_group]

        # 1. SNP Hit Distribution
        def count_hits(row):
            hits = 0
            for i in range(num_snps):
                if row[f"SNP{i+1}"] == group_target_bases[i]:
                    hits += 1
            return hits

        group_df = group_df.copy()
        group_df["SNP_Hits"] = group_df.apply(count_hits, axis=1)
        hit_dist = Counter(group_df["SNP_Hits"])

        # 2. Co-occurrence frequencies
        co_occ = {}
        # Individual bits
        for i in range(num_snps):
            count = len(group_df[group_df[f"SNP{i+1}"] == group_target_bases[i]])
            co_occ[f"SNP{i+1}_only"] = count / group_size

        # Pairs (top 3 pairs for brevity if many SNPs, but let's do all if few)
        if num_snps >= 2:
            pairs = list(itertools.combinations(range(num_snps), 2))
            for p in pairs:
                idx1, idx2 = p
                mask = (group_df[f"SNP{idx1+1}"] == group_target_bases[idx1]) & (group_df[f"SNP{idx2+1}"] == group_target_bases[idx2])
                count = len(group_df[mask])
                co_occ[f"SNP{idx1+1}_&_SNP{idx2+1}"] = count / group_size

        # All
        mask_all = True
        for i in range(num_snps):
            mask_all &= (group_df[f"SNP{i+1}"] == group_target_bases[i])
        count_all = len(group_df[mask_all])
        co_occ["ALL_SNPS"] = count_all / group_size

        res = {
            "Sample": sample_name,
            "Allele": allele_group,
            "Total_Reads": group_size,
        }
        for hit in range(1, num_snps + 1):
            res[f"{hit}_SNP_Hits_%"] = (hit_dist.get(hit, 0) / group_size) * 100
        res.update(co_occ)
        results.append(res)

    return results

def main():
    parser = argparse.ArgumentParser(description="Analyze stoichiometry and co-occurrence from quantification results.")
    parser.add_argument("--output_dir", default=".", help="Project output directory (default: current directory)")
    parser.add_argument("--snp_file", default=None, help="Path to SNP JSON (default: results/ref_seq/snps.json)")
    args = parser.parse_args()

    results_dir = os.path.join(args.output_dir, "results")
    quant_dir = os.path.join(results_dir, "quantification")
    stoic_dir = os.path.join(results_dir, "stoichiometry")
    snp_file = args.snp_file if args.snp_file else os.path.join(results_dir, "ref_seq", "snps.json")

    quant_files = glob.glob(os.path.join(quant_dir, "*_quant.csv"))

    if not quant_files:
        print(f"No quantification files found in {quant_dir}")
        return

    # Load SNPs dynamically
    if os.path.exists(snp_file):
        with open(snp_file, "r") as f:
            snp_data = json.load(f)
        b6_bases = [s["b6"].split("/")[0] for s in snp_data]
        cast_bases = [s["cast"].split("/")[0] for s in snp_data]
        num_snps = len(snp_data)
        target_bases = {"B6": b6_bases, "Cast": cast_bases}
    else:
        print(f"Error: {snp_file} not found.")
        return

    if num_snps == 0:
        print("No SNPs to analyze.")
        return

    os.makedirs(stoic_dir, exist_ok=True)
    all_results = []

    for f in quant_files:
        print(f"Analyzing stoichiometry for {f}...")
        all_results.extend(analyze_stoichiometry(f, num_snps, target_bases))

    summary_df = pd.DataFrame(all_results)
    print("\nStoichiometry & Co-occurrence Summary:")
    print(summary_df.to_string())

if __name__ == "__main__":
    main()
