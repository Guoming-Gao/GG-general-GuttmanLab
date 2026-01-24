import pandas as pd
import os
import glob
import json
import argparse
from collections import Counter
import itertools

def analyze_stoichiometry(file_path, snps, target_bases):
    sample_name = os.path.basename(file_path).replace("_quant.csv", "")
    df = pd.read_csv(file_path)
    num_snps = len(snps)

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
            for i, s in enumerate(snps):
                if row[s["label"]] == group_target_bases[i]:
                    hits += 1
            return hits

        group_df = group_df.copy()
        group_df["SNP_Hits"] = group_df.apply(count_hits, axis=1)
        hit_dist = Counter(group_df["SNP_Hits"])

        # 2. Co-occurrence frequencies
        co_occ = {}
        # Individual bits
        for i, s in enumerate(snps):
            count = len(group_df[group_df[s["label"]] == group_target_bases[i]])
            co_occ[f"{s['label']}_only"] = count / group_size

        # Pairs
        if num_snps >= 2:
            pairs = list(itertools.combinations(range(num_snps), 2))
            for p in pairs:
                idx1, idx2 = p
                label1, label2 = snps[idx1]["label"], snps[idx2]["label"]
                mask = (group_df[label1] == group_target_bases[idx1]) & (group_df[label2] == group_target_bases[idx2])
                count = len(group_df[mask])
                co_occ[f"{label1}_&_{label2}"] = count / group_size

        # All
        mask_all = True
        for i, s in enumerate(snps):
            mask_all &= (group_df[s["label"]] == group_target_bases[i])
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
    parser.add_argument("--results_dir", default=None, help="Results directory (default: ./results)")
    parser.add_argument("--snp_file", default=None, help="Path to SNP JSON (default: results/ref_seq/snps.json)")
    parser.add_argument("--omit_snps", type=str, default=None, help="Comma-separated SNP indices to omit (1-based, e.g., '2')")
    args = parser.parse_args()

    results_dir = args.results_dir if args.results_dir else os.path.join(".", "results")
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

        omit = []
        if args.omit_snps:
            omit = [int(i.strip()) for i in args.omit_snps.split(",")]

        snps = []
        for i, s in enumerate(snp_data):
            idx = i + 1
            if idx in omit:
                continue
            snps.append({
                "label": f"SNP{idx}_{s['local_pos']}",
                "b6": s["b6"].split("/")[0],
                "cast": s["cast"].split("/")[0]
            })

        b6_bases = [s["b6"] for s in snps]
        cast_bases = [s["cast"] for s in snps]
        num_snps = len(snps)
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
        all_results.extend(analyze_stoichiometry(f, snps, target_bases))

    summary_df = pd.DataFrame(all_results)
    summary_path = os.path.join(stoic_dir, "stoichiometry_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"\nStoichiometry & Co-occurrence Summary saved to {summary_path}")
    print(summary_df.to_string())

if __name__ == "__main__":
    main()
