import pysam
import os
import glob
import pandas as pd
import json
import argparse
from collections import Counter

def get_base_at_pos(read, ref_pos):
    """Extract the base from a read at a specific reference position."""
    for q_pos, r_pos in read.get_aligned_pairs():
        if r_pos == ref_pos:
            if q_pos is None:
                return "-"
            return read.query_sequence[q_pos]
    return "N"

def assign_allele(snp_bases, snps, min_matches=None):
    """Assign allele based on SNP matches using a configurable threshold."""
    if not snps:
        return "Unknown"

    b6_matches = 0
    cast_matches = 0

    for i, base in enumerate(snp_bases):
        if base == snps[i]["b6"]:
            b6_matches += 1
        elif base == snps[i]["cast"]:
            cast_matches += 1

    # Majority rule threshold or custom
    if min_matches is not None:
        threshold = min_matches
    else:
        threshold = (len(snps) // 2) + 1

    if b6_matches > cast_matches:
        if b6_matches >= threshold:
             return "B6"
        return "B6_low_conf"
    elif cast_matches > b6_matches:
        if cast_matches >= threshold:
            return "Cast"
        return "Cast_low_conf"
    else:
        if b6_matches > 0:
            return "Ambiguous"
        return "Noise/Unknown"

def process_bam(file_path, output_dir, snps, min_matches=None):
    sample_name = os.path.basename(file_path).replace(".sorted.bam", "")
    results = []

    with pysam.AlignmentFile(file_path, "rb") as sam:
        for read in sam.fetch("Xist_Amplicon"):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue

            bases = [get_base_at_pos(read, s["pos"] - 1) for s in snps]
            allele = assign_allele(bases, snps, min_matches=min_matches)

            res = {"ReadID": read.query_name, "Allele": allele}
            for i, b in enumerate(bases):
                label = snps[i].get("label", f"SNP{i+1}")
                res[label] = b
            results.append(res)

    if not results:
        print(f"No valid reads found in {file_path}")
        return {
            "Sample": sample_name,
            "Total_Primary_Reads": 0,
            "B6": 0,
            "Cast": 0,
            "Ambiguous": 0,
            "Noise": 0,
            "Cast_Ratio": 0
        }

    df = pd.DataFrame(results)
    os.makedirs(output_dir, exist_ok=True)
    df.to_csv(os.path.join(output_dir, f"{sample_name}_quant.csv"), index=False)

    summary = Counter(df["Allele"])
    summary_stats = {
        "Sample": sample_name,
        "Total_Primary_Reads": len(df),
        "B6": summary.get("B6", 0) + summary.get("B6_low_conf", 0),
        "Cast": summary.get("Cast", 0) + summary.get("Cast_low_conf", 0),
        "Ambiguous": summary.get("Ambiguous", 0),
        "Noise": summary.get("Noise/Unknown", 0)
    }

    denom = summary_stats["B6"] + summary_stats["Cast"]
    summary_stats["Cast_Ratio"] = summary_stats["Cast"] / denom if denom > 0 else 0
    return summary_stats

def main():
    parser = argparse.ArgumentParser(description="Quantify alleles from aligned BAM files.")
    parser.add_argument("--results_dir", default=None, help="Results directory (default: ./results)")
    parser.add_argument("--snp_file", default=None, help="Path to SNP JSON (default: results/ref_seq/snps.json)")
    parser.add_argument("--omit_snps", type=str, default=None, help="Comma-separated SNP indices to omit (1-based, e.g., '2')")
    parser.add_argument("--min_matches", type=int, default=None, help="Minimum SNP matches for high-confidence assignment")
    args = parser.parse_args()

    results_dir = args.results_dir if args.results_dir else os.path.join(".", "results")
    align_dir = os.path.join(results_dir, "aligned")
    quant_dir = os.path.join(results_dir, "quantification")
    snp_file = args.snp_file if args.snp_file else os.path.join(results_dir, "ref_seq", "snps.json")

    bam_files = glob.glob(os.path.join(align_dir, "*.sorted.bam"))

    if not bam_files:
        print(f"No sorted BAM files found in {align_dir}")
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
                print(f"Omit SNP {idx} at local pos {s['local_pos']}")
                continue

            snps.append({
                "pos": s["local_pos"],
                "b6": s["b6"].split("/")[0], # GT 0/0 -> 0
                "cast": s["cast"].split("/")[0],
                "label": f"SNP{idx}_{s['local_pos']}"
            })
    else:
        print(f"Error: {snp_file} not found.")
        return

    all_summaries = []
    for f in bam_files:
        print(f"Processing {f}...")
        all_summaries.append(process_bam(f, quant_dir, snps, min_matches=args.min_matches))

    summary_df = pd.DataFrame(all_summaries)
    print("\nAllele Quantification Summary:")
    print(summary_df.to_string())

if __name__ == "__main__":
    main()
