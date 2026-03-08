
import os
import argparse
import pysam
from Bio.Seq import Seq
import pandas as pd
from collections import Counter

def trim_mode2(fastq_in, output_dir, f_primer, r_primer):
    os.makedirs(output_dir, exist_ok=True)
    sample = os.path.basename(fastq_in).split(".")[0]
    out_fastq = os.path.join(output_dir, f"{sample}.trimmed.fastq")
    out_stats = os.path.join(output_dir, f"{sample}.trimming_stats.csv")

    # Landmarks
    ADAPTER_F = "CAGACGTGTGCTCTTCCGATCT"
    ADAPTER_R = "GATCGGAAGAGCACACGTCTGA" # Universal tail

    PRIMER_F = f_primer.upper()
    PRIMER_R = r_primer.upper()
    PRIMER_F_RC = str(Seq(PRIMER_F).reverse_complement()).upper()
    PRIMER_R_RC = str(Seq(PRIMER_R).reverse_complement()).upper()

    counts = Counter()
    stats = []

    print(f"Trimming {sample}...")

    with pysam.FastxFile(fastq_in) as fx, open(out_fastq, "w") as out:
        for entry in fx:
            counts["Total"] += 1
            seq = entry.sequence.upper()

            # Detect landmarks
            has_head_adapter = ADAPTER_F in seq
            has_tail_adapter = ADAPTER_R in seq

            # Since orientation can be flipped, check for primers
            # Orientation 1: ADP_F - UMI - PRIMER_R ... PRIMER_F_RC - ADP_R
            # Orientation 2: ADP_F - UMI - PRIMER_F ... PRIMER_R_RC - ADP_R

            has_p_f = (PRIMER_F in seq) or (PRIMER_F_RC in seq)
            has_p_r = (PRIMER_R in seq) or (PRIMER_R_RC in seq)

            counts["HasHeadAdapter"] += 1 if has_head_adapter else 0
            counts["HasTailAdapter"] += 1 if has_tail_adapter else 0
            counts["HasPrimerF"] += 1 if has_p_f else 0
            counts["HasPrimerR"] += 1 if has_p_r else 0

            # Logic for extraction (focused on the major pattern observed)
            f_idx = seq.find(ADAPTER_F)
            b_idx = seq.find(ADAPTER_R)

            umi = ""
            if has_head_adapter:
                umi_start = f_idx + len(ADAPTER_F)
                umi = seq[umi_start:umi_start+6]
                counts["HasUMI"] += 1

            if has_head_adapter and has_tail_adapter:
                # Trim sequence: from end of UMI to start of Back Adapter
                new_seq = entry.sequence[f_idx + len(ADAPTER_F) + 6:b_idx]
                new_qual = entry.quality[f_idx + len(ADAPTER_F) + 6:b_idx]

                if len(new_seq) > 50:
                    new_name = f"{entry.name}_{umi}"
                    out.write(f"@{new_name}\n{new_seq}\n+\n{new_qual}\n")
                    counts["Success"] += 1
                    status = "Success"
                else:
                    counts["TooShort"] += 1
                    status = "TooShort"
            else:
                status = "MissingAdapters"

            stats.append({
                "ReadID": entry.name,
                "HeadAdapter": has_head_adapter,
                "UMI": umi != "",
                "PrimerF": has_p_f,
                "PrimerR": has_p_r,
                "TailAdapter": has_tail_adapter,
                "Status": status
            })

    # Summary per sample
    summary = {
        "Sample": sample,
        "TotalReads": counts["Total"],
        "HeadAdapter_Count": counts["HasHeadAdapter"],
        "HeadAdapter_Pct": counts["HasHeadAdapter"] / counts["Total"],
        "UMI_Count": counts["HasUMI"],
        "UMI_Pct": counts["HasUMI"] / counts["Total"],
        "PrimerF_Count": counts["HasPrimerF"],
        "PrimerF_Pct": counts["HasPrimerF"] / counts["Total"],
        "PrimerR_Count": counts["HasPrimerR"],
        "PrimerR_Pct": counts["HasPrimerR"] / counts["Total"],
        "TailAdapter_Count": counts["HasTailAdapter"],
        "TailAdapter_Pct": counts["HasTailAdapter"] / counts["Total"],
        "FullStructure_Count": counts["Success"],
        "FullStructure_Pct": counts["Success"] / counts["Total"]
    }

    pd.DataFrame([summary]).to_csv(out_stats.replace(".csv", "_summary.csv"), index=False)
    # pd.DataFrame(stats).to_csv(out_stats, index=False) # Large file, skip for now unless needed

    return summary

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--f_primer", required=True)
    parser.add_argument("--r_primer", required=True)
    args = parser.parse_args()
    trim_mode2(args.fastq, args.output_dir, args.f_primer, args.r_primer)
