
import pysam
import os
import pandas as pd
import argparse
from collections import Counter

def get_base_at_pos(read, target_ref_pos):
    """
    Extract the base from a read at a specific genomic position.
    target_ref_pos is 0-based.
    """
    # get_aligned_pairs returns (query_pos, ref_pos)
    # query_pos is 0-based index into read.query_sequence
    # ref_pos is 0-based index into the reference
    for q_pos, r_pos in read.get_aligned_pairs():
        if r_pos == target_ref_pos:
            if q_pos is None:
                # Deletion at this position
                return "-"
            return read.query_sequence[q_pos]
    return "N" # Not covered by alignment

def assign_allele(snp_bases, snp_list):
    b6_matches = 0
    cast_matches = 0

    for i, base in enumerate(snp_bases):
        if base == snp_list[i]['b6']:
            b6_matches += 1
        elif base == snp_list[i]['cast']:
            cast_matches += 1

    if b6_matches > cast_matches:
        return "B6", b6_matches, cast_matches
    elif cast_matches > b6_matches:
        return "Cast", b6_matches, cast_matches
    else:
        if b6_matches > 0:
            return "Ambiguous", b6_matches, cast_matches
        return "Noise", b6_matches, cast_matches

def quantify_alleles(bam_dir, snp_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Load SNPs
    snps_df = pd.read_csv(snp_file, sep="\t")
    snp_list = []
    for _, row in snps_df.iterrows():
        snp_list.append({
            'pos': row['POS'] - 1, # Convert to 0-based
            'b6': row['B6_BASE'],
            'cast': row['CAST_BASE']
        })

    bam_files = [f for f in os.listdir(bam_dir) if f.endswith(".xist.primary.bam")]

    all_summaries = []

    for f in bam_files:
        sample = f.replace(".xist.primary.bam", "")
        print(f"Quantifying alleles for {sample}...")
        results = []

        with pysam.AlignmentFile(os.path.join(bam_dir, f), "rb") as sam:
            for read in sam.fetch():
                if read.is_secondary or read.is_supplementary or read.is_unmapped:
                    continue

                snp_bases = [get_base_at_pos(read, s['pos']) for s in snp_list]
                allele, b6_cnt, cast_cnt = assign_allele(snp_bases, snp_list)

                res = {
                    "ReadID": read.query_name,
                    "Allele": allele,
                    "B6_Matches": b6_cnt,
                    "Cast_Matches": cast_cnt,
                    "Total_SNP_Coverage": sum(1 for b in snp_bases if b != "N")
                }
                # Add individual SNP bases for transparency
                for i, b in enumerate(snp_bases):
                    res[f"SNP_{snp_list[i]['pos']+1}"] = b

                results.append(res)

        if results:
            df = pd.DataFrame(results)
            df.to_csv(os.path.join(output_dir, f"{sample}_quant.csv"), index=False)

            summary = Counter(df["Allele"])
            total = len(df)
            all_summaries.append({
                "Sample": sample,
                "Total_Reads": total,
                "B6": summary.get("B6", 0),
                "Cast": summary.get("Cast", 0),
                "Ambiguous": summary.get("Ambiguous", 0),
                "Noise": summary.get("Noise", 0),
                "Cast_Ratio": summary.get("Cast", 0) / (summary.get("B6", 0) + summary.get("Cast", 0)) if (summary.get("B6", 0) + summary.get("Cast", 0)) > 0 else 0
            })

    summary_df = pd.DataFrame(all_summaries)
    summary_df.to_csv(os.path.join(output_dir, "allele_quantification_summary.csv"), index=False)
    return summary_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam_dir", required=True)
    parser.add_argument("--snp_file", required=True)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    quantify_alleles(args.bam_dir, args.snp_file, args.output_dir)
