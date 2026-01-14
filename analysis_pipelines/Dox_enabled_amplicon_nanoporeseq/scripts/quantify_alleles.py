
import pysam
import os
import glob
import pandas as pd
import json
from collections import Counter

ALIGN_DIR = "aligned"
BAM_FILES = glob.glob(os.path.join(ALIGN_DIR, "*.sorted.bam"))
OUTPUT_DIR = "quantification"
SNP_FILE = "ref_seq/snps.json"

# Load SNPs dynamically
if os.path.exists(SNP_FILE):
    with open(SNP_FILE, "r") as f:
        SNP_DATA = json.load(f)
    SNPS = []
    for i, s in enumerate(SNP_DATA):
        SNPS.append({
            "pos": s["local_pos"],
            "b6": s["b6"].split("/")[0], # GT 0/0 -> 0
            "cast": s["cast"].split("/")[0],
            "label": f"SNP{i+1}_{s['local_pos']}"
        })
else:
    print(f"Warning: {SNP_FILE} not found. Using empty SNP list.")
    SNPS = []

def get_base_at_pos(read, ref_pos):
    """Extract the base from a read at a specific reference position."""
    for q_pos, r_pos in read.get_aligned_pairs():
        if r_pos == ref_pos:
            if q_pos is None:
                return "-"
            return read.query_sequence[q_pos]
    return "N"

def assign_allele(snp_bases):
    """Assign allele based on SNP matches using majority rule."""
    if not SNPS:
        return "Unknown"

    b6_matches = 0
    cast_matches = 0

    for i, base in enumerate(snp_bases):
        if base == SNPS[i]["b6"]:
            b6_matches += 1
        elif base == SNPS[i]["cast"]:
            cast_matches += 1

    # Majority rule threshold
    threshold = (len(SNPS) // 2) + 1

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

def process_bam(file_path):
    sample_name = os.path.basename(file_path).replace(".sorted.bam", "")
    results = []

    with pysam.AlignmentFile(file_path, "rb") as sam:
        for read in sam.fetch("Xist_Amplicon"):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue

            bases = [get_base_at_pos(read, s["pos"] - 1) for s in SNPS]
            allele = assign_allele(bases)

            res = {"ReadID": read.query_name, "Allele": allele}
            for i, b in enumerate(bases):
                res[f"SNP{i+1}"] = b
            results.append(res)

    df = pd.DataFrame(results)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df.to_csv(os.path.join(OUTPUT_DIR, f"{sample_name}_quant.csv"), index=False)

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
    if not SNPS:
        print("No SNPs defined. Exiting.")
        return

    all_summaries = []
    for f in BAM_FILES:
        print(f"Processing {f}...")
        all_summaries.append(process_bam(f))

    summary_df = pd.DataFrame(all_summaries)
    summary_df.to_csv("allele_quantification_summary.csv", index=False)
    print("\nAllele Quantification Summary:")
    print(summary_df.to_string())

if __name__ == "__main__":
    main()
