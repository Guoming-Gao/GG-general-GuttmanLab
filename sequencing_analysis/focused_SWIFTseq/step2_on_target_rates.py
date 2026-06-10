
import pysam
import pandas as pd
import os

# Paths
base_data_dir = "/Volumes/guttman/users/gmgao/Data_seq/20251214_vs_1218-focusedSwift-RTprimerconc_UMIsplintconc/"
# Result paths
results_base = os.path.join(base_data_dir, "results")
processed_probes_path = os.path.join(results_base, "intermediate/processed_probes.csv")
reports_dir = os.path.join(results_base, "reports")
output_dir = os.path.join(results_base, "intermediate")

# Load probes
probes_df = pd.read_csv(processed_probes_path)

# Dataset mapping
dataset_map = {
    "1214-A-50CRT_wRT_400pmolPrimer_1uMsplint": "400pmol_Primer_1uM_Splint",
    "1214-B-50CRT_wRT_40pmolPrimer_1uMsplint": "40pmol_Primer_1uM_Splint",
    "1214-C-50CRT_wRT_10pmolPrimer_1uMsplint": "10pmol_Primer_1uM_Splint_1214",
    "1218-A-50CRT_wRT_10pmolPrimer_1uMsplint-Gelcut": "10pmol_Primer_1uM_Splint_1218"
}

datasets = list(dataset_map.keys())

results = []

def count_strict_on_target(bam_path, probes_df):
    on_target_qnames = set()
    total_reads = 0
    if not os.path.exists(bam_path):
        return 0, 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        total_reads = bam.mapped + bam.unmapped

        for index, row in probes_df.iterrows():
            chrom = str(row['Chromosome'])
            if chrom not in bam.references and "chr"+chrom in bam.references:
                chrom = "chr" + chrom

            start = int(row['WindowStart'])
            end = int(row['WindowEnd'])
            strand = row['Strand']

            try:
                for read in bam.fetch(chrom, max(0, start), end):
                    # Strict 5' end check
                    # RNA (+): cDNA is (-). 5' end of read is reference_end - 1.
                    # RNA (-): cDNA is (+). 5' end of read is reference_start.

                    if strand == '+':
                        if read.is_reverse and start <= (read.reference_end - 1) <= end:
                            on_target_qnames.add(read.query_name)
                    elif strand == '-':
                        if not read.is_reverse and start <= read.reference_start <= end:
                            on_target_qnames.add(read.query_name)

            except ValueError:
                pass

    return len(on_target_qnames), total_reads

for ds in datasets:
    raw_bam = os.path.join(base_data_dir, ds, "aligned.mouse.sorted.bam")
    dedup_bam = os.path.join(base_data_dir, ds, "aligned.mouse.sorted.bam.dedup.bam")

    print(f"Processing {ds} (Strict On-Target)...")

    on_target_raw, total_raw = count_strict_on_target(raw_bam, probes_df)
    on_target_dedup, total_dedup = count_strict_on_target(dedup_bam, probes_df)

    results.append({
        "Dataset": dataset_map[ds],
        "Original_Name": ds,
        "TotalReads_Raw": total_raw,
        "TotalReads_Dedup": total_dedup,
        "OnTarget_Raw": on_target_raw,
        "OnTarget_Dedup": on_target_dedup,
        "OnTargetRateRaw_Percent": (on_target_raw / total_raw * 100) if total_raw > 0 else 0,
        "OnTargetRateDedup_Percent": (on_target_dedup / total_dedup * 100) if total_dedup > 0 else 0
    })

results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(output_dir, "on_target_rates.csv"), index=False)
print("Step 2 Complete.")
