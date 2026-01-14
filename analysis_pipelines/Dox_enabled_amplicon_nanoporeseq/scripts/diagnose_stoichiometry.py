
import pysam
import os
import glob
import pandas as pd
import json

ALIGN_DIR = "/Users/gmgao/GGscripts/GG-general-GuttmanLab/analysis_pipelines/Dox_enabled_amplicon_nanoporeseq/aligned"
SNP_FILE = "/Users/gmgao/GGscripts/GG-general-GuttmanLab/analysis_pipelines/Dox_enabled_amplicon_nanoporeseq/ref_seq/snps.json"
BAM_FILES = glob.glob(os.path.join(ALIGN_DIR, "*.sorted.bam"))

with open(SNP_FILE, "r") as f:
    SNPS = json.load(f)

def get_base(read, pos_0):
    for q_pos, r_pos in read.get_aligned_pairs():
        if r_pos == pos_0:
            if q_pos is None: return "-"
            return read.query_sequence[q_pos]
    return "N"

def diagnostic(bam_path):
    sample = os.path.basename(bam_path)
    print(f"\n--- Diagnostic for {sample} ---")

    # Track hits for ALL primary reads
    all_hits_b6 = []
    all_hits_cast = []

    with pysam.AlignmentFile(bam_path, "rb") as sam:
        count = 0
        for read in sam.fetch("Xist_Amplicon"):
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue

            bases = [get_base(read, s["local_pos"] - 1) for s in SNPS]
            b6_vals = [s["b6"].split("/")[0] for s in SNPS]
            cast_vals = [s["cast"].split("/")[0] for s in SNPS]

            b6_hits = sum(1 for i, b in enumerate(bases) if b == b6_vals[i])
            cast_hits = sum(1 for i, b in enumerate(bases) if b == cast_vals[i])

            all_hits_b6.append(b6_hits)
            all_hits_cast.append(cast_hits)

            if count < 5: # Print first 5 reads
                print(f"Read: {read.query_name[:10]} | Strand: {'-' if read.is_reverse else '+'} | Bases: {''.join(bases)} | B6 Hits: {b6_hits} | Cast Hits: {cast_hits}")
            count += 1

    # Frequency of hits
    from collections import Counter
    print("\nB6 Hit Distribution (Sample-wide):")
    print(dict(sorted(Counter(all_hits_b6).items())))
    print("\nCast Hit Distribution (Sample-wide):")
    print(dict(sorted(Counter(all_hits_cast).items())))

if __name__ == "__main__":
    # Test on WT (diff) and WT (Dox)
    for b in BAM_FILES[:2]:
        diagnostic(b)
