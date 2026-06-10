
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Paths
base_data_dir = "/Volumes/guttman/users/gmgao/Data_seq/20251214_vs_1218-focusedSwift-RTprimerconc_UMIsplintconc/"
# Result paths
results_base = os.path.join(base_data_dir, "results")
processed_probes_path = os.path.join(results_base, "intermediate/processed_probes.csv")
stats_path = os.path.join(results_base, "intermediate/probe_specificity_stats.csv")
spot_check_dir = os.path.join(results_base, "reports/spot_checks")
output_dir = os.path.join(results_base, "intermediate")
os.makedirs(spot_check_dir, exist_ok=True)

# Load data
probes_df = pd.read_csv(processed_probes_path)
stats_df = pd.read_csv(stats_path)

# Datasets to compare (Standardized names mapping to folders)
ds_c_folder = "1214-C-50CRT_wRT_10pmolPrimer_1uMsplint"
ds_a_folder = "1218-A-50CRT_wRT_10pmolPrimer_1uMsplint-Gelcut"

def get_coverage(bam_path, chrom, start, end):
    if not os.path.exists(bam_path):
        return None
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if chrom not in bam.references and "chr"+chrom in bam.references:
            chrom = "chr" + chrom
        cov = np.zeros(max(0, end - start + 1))
        try:
            for read in bam.fetch(chrom, max(0, start), end):
                r_start = max(start, read.reference_start)
                r_end = min(end, read.reference_end)
                if r_start < r_end:
                    cov[r_start - start : r_end - start] += 1
        except:
            pass
    return cov

report_plots = []

for idx, row in probes_df.iterrows():
    gene = row['GeneName']
    chrom = str(row['Chromosome'])
    strand = row['Strand']
    probe_start, probe_end = int(row['theStartPos']), int(row['theEndPos'])
    bin_start, bin_end = int(row['WindowStart']), int(row['WindowEnd'])

    # Get status
    status = stats_df[stats_df['ProbeID'] == idx]['Status'].values[0]

    # View window (+/- 500 bp around probe for context)
    w_start = probe_start - 500
    w_end = probe_end + 500

    print(f"Plotting probe {idx}: {gene} ({status})...")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 7), sharey=True)
    x = np.arange(w_start, w_end + 1)

    # Dataset 1: 10pmol (1214)
    c_raw = get_coverage(os.path.join(base_data_dir, ds_c_folder, "aligned.mouse.sorted.bam"), chrom, w_start, w_end)
    c_dedup = get_coverage(os.path.join(base_data_dir, ds_c_folder, "aligned.mouse.sorted.bam.dedup.bam"), chrom, w_start, w_end)

    if c_raw is not None:
        xc = x[:len(c_raw)]
        ax1.fill_between(xc, c_raw, color='blue', alpha=0.1, label='Raw')
        ax1.plot(xc, c_raw, color='blue', linewidth=0.8, alpha=0.5)
        ax1.fill_between(xc, c_dedup, color='red', alpha=0.3, label='Dedup')
        ax1.plot(xc, c_dedup, color='red', linewidth=1)
        ax1.axvspan(probe_start, probe_end, color='green', alpha=0.2, label='Probe')
        ax1.axvspan(bin_start, bin_end, color='orange', alpha=0.3, label='Init Bin (5-15nt)')
        ax1.set_title("10pmol_Primer_1uM_Splint_1214", fontsize=11)
        ax1.set_ylabel("Read Depth")
        ax1.legend(loc='upper right', fontsize=8)
        ax1.grid(True, linestyle='--', alpha=0.3)

    # Dataset 2: 10pmol (1218)
    a_raw = get_coverage(os.path.join(base_data_dir, ds_a_folder, "aligned.mouse.sorted.bam"), chrom, w_start, w_end)
    a_dedup = get_coverage(os.path.join(base_data_dir, ds_a_folder, "aligned.mouse.sorted.bam.dedup.bam"), chrom, w_start, w_end)

    if a_raw is not None:
        xa = x[:len(a_raw)]
        ax2.fill_between(xa, a_raw, color='blue', alpha=0.1, label='Raw')
        ax2.plot(xa, a_raw, color='blue', linewidth=0.8, alpha=0.5)
        ax2.fill_between(xa, a_dedup, color='red', alpha=0.3, label='Dedup')
        ax2.plot(xa, a_dedup, color='red', linewidth=1)
        ax2.axvspan(probe_start, probe_end, color='green', alpha=0.2, label='Probe')
        ax2.axvspan(bin_start, bin_end, color='orange', alpha=0.3, label='Init Bin (5-15nt)')
        ax2.set_title("10pmol_Primer_1uM_Splint_1218", fontsize=11)
        ax2.legend(loc='upper right', fontsize=8)
        ax2.grid(True, linestyle='--', alpha=0.3)

    plt.suptitle(f"Probe {idx}: {gene} ({strand}) | Status: {status}\nWindow: {chrom}:{w_start}-{w_end}", fontsize=14, y=1.02)
    plt.tight_layout()

    plot_fn = f"probe_{idx}_{gene}.png"
    plt.savefig(os.path.join(spot_check_dir, plot_fn), bbox_inches='tight')
    plt.close()
    report_plots.append((idx, gene, status, plot_fn))

# Output a simple manifest for the summary script
plots_manifest = pd.DataFrame(report_plots, columns=["ID", "Gene", "Status", "Filename"])
plots_manifest.to_csv(os.path.join(output_dir, "plots_manifest.csv"), index=False)

print(f"Generated {len(report_plots)} comparative plots in {spot_check_dir}")
