
import pysam
import os
import glob
import pandas as pd
from collections import Counter

ALIGN_DIR = "aligned"
BAM_FILES = glob.glob(os.path.join(ALIGN_DIR, "*.sorted.bam"))
OUTPUT_DIR = "quantification"

# SNP Definitions (from Step 1 Verification)
# Positions are 1-indexed in GenBank/Reference FASTA
SNPS = [
    {"pos": 67, "b6": "T", "cast": "C", "label": "SNP1_67"},
    {"pos": 334, "b6": "G", "cast": "A", "label": "SNP2_334"},
    {"pos": 553, "b6": "G", "cast": "A", "label": "SNP3_553"}
]

def get_base_at_pos(read, ref_pos):
    """
    Extract the base from a read at a specific reference position.
    ref_pos is 0-indexed relative to the reference start.
    """
    # get_aligned_pairs(with_seq=True) returns (query_pos, ref_pos, ref_base)
    # query_pos is None if it's a deletion
    # ref_pos is None if it's an insertion
    for q_pos, r_pos in read.get_aligned_pairs():
        if r_pos == ref_pos:
            if q_pos is None:
                return "-" # Deletion
            return read.query_sequence[q_pos]
    return "N" # Not covered

def assign_allele(snp_bases):
    """
    Assign allele based on SNP matches.
    """
    b6_matches = 0
    cast_matches = 0

    for i, base in enumerate(snp_bases):
        if base == SNPS[i]["b6"]:
            b6_matches += 1
        elif base == SNPS[i]["cast"]:
            cast_matches += 1

    if b6_matches > cast_matches:
        if b6_matches >= 2: # Majority rule
             return "B6"
        return "B6_low_conf"
    elif cast_matches > b6_matches:
        if cast_matches >= 2:
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
        # We only care about primary alignments for allele assignment
        for read in sam.fetch("Xist_Amplicon"):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue

            # Extract bases at SNP positions (Convert 1-indexed to 0-indexed)
            bases = [get_base_at_pos(read, s["pos"] - 1) for s in SNPS]
            allele = assign_allele(bases)

            results.append({
                "ReadID": read.query_name,
                "SNP1": bases[0],
                "SNP2": bases[1],
                "SNP3": bases[2],
                "Allele": allele
            })

    df = pd.DataFrame(results)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df.to_csv(os.path.join(OUTPUT_DIR, f"{sample_name}_quant.csv"), index=False)

    # Summary
    summary = Counter(df["Allele"])
    total = len(df)
    summary_stats = {
        "Sample": sample_name,
        "Total_Primary_Reads": total,
        "B6": summary.get("B6", 0) + summary.get("B6_low_conf", 0),
        "Cast": summary.get("Cast", 0) + summary.get("Cast_low_conf", 0),
        "Ambiguous": summary.get("Ambiguous", 0),
        "Noise": summary.get("Noise/Unknown", 0)
    }

    # Calculate ratio (Cast / (B6 + Cast))
    denom = summary_stats["B6"] + summary_stats["Cast"]
    summary_stats["Cast_Ratio"] = summary_stats["Cast"] / denom if denom > 0 else 0

    return summary_stats

def main():
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
