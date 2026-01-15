
import pysam
import pandas as pd
import numpy as np
import os
from scipy import stats

# Paths
base_data_dir = "/Volumes/guttman/users/gmgao/Data_seq/20251214_vs_1218-focusedSwift-RTprimerconc_UMIsplintconc/"
# Result paths
results_base = os.path.join(base_data_dir, "results")
processed_probes_path = os.path.join(results_base, "intermediate/processed_probes.csv")
output_dir = os.path.join(results_base, "intermediate")
os.makedirs(output_dir, exist_ok=True)

# Load probes
probes_df = pd.read_csv(processed_probes_path)

datasets = {
    "1214_10pmol": "1214-C-50CRT_wRT_10pmolPrimer_1uMsplint",
    "1218_10pmol": "1218-A-50CRT_wRT_10pmolPrimer_1uMsplint-Gelcut"
}

def get_start_counts(bam_path, chrom, start, end, gene_strand):
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if chrom not in bam.references and "chr"+chrom in bam.references:
            chrom = "chr" + chrom
        start_counts = np.zeros(max(0, end - start + 1))
        try:
            for read in bam.fetch(chrom, max(0, start), end):
                if gene_strand == '+':
                    if read.is_reverse:
                        pos = read.reference_end - 1
                        if start <= pos <= end:
                            start_counts[pos - start] += 1
                elif gene_strand == '-':
                    if not read.is_reverse:
                        pos = read.reference_start
                        if start <= pos <= end:
                            start_counts[pos - start] += 1
        except:
            pass
    return start_counts

probe_stats = []

for idx, row in probes_df.iterrows():
    gene = row['GeneName']
    chrom = str(row['Chromosome'])
    strand = row['Strand']
    p_start, p_end = int(row['theStartPos']), int(row['theEndPos'])
    bin_start, bin_end = int(row['WindowStart']), int(row['WindowEnd'])

    bg_start = p_start - 500
    bg_end = p_end + 500

    stats_entry = {"ProbeID": idx, "GeneName": gene, "Strand": strand}

    for label, ds_folder in datasets.items():
        raw_bam = os.path.join(base_data_dir, ds_folder, "aligned.mouse.sorted.bam")
        if not os.path.exists(raw_bam): continue

        bg_starts = get_start_counts(raw_bam, chrom, bg_start, bg_end, strand)
        rel_start = bin_start - bg_start
        rel_end = bin_end - bg_start
        bin_starts = bg_starts[max(0, rel_start) : min(len(bg_starts), rel_end + 1)]

        peak_val = np.max(bin_starts) if len(bin_starts) > 0 else 0
        mode_val = stats.mode(bg_starts, keepdims=True).mode[0] if len(bg_starts) > 0 else 0

        stats_entry[f"{label}_Peak_Starts"] = peak_val
        stats_entry[f"{label}_BG_Starts_Mode"] = mode_val

    accepted = False
    for label in datasets.keys():
        peak = stats_entry.get(f"{label}_Peak_Starts", 0)
        mode = stats_entry.get(f"{label}_BG_Starts_Mode", 0)
        if peak >= 3 and peak > mode + 2:
            accepted = True
            break

    stats_entry["Status"] = "Accepted" if accepted else "Rejected"
    probe_stats.append(stats_entry)

stats_df = pd.DataFrame(probe_stats)
stats_df.to_csv(os.path.join(output_dir, "probe_specificity_stats.csv"), index=False)
print("Step 6 Complete.")
