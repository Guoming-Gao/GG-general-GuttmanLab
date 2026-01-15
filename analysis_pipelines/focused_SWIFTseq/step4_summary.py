
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Result paths
base_data_dir = "/Volumes/guttman/users/gmgao/Data_seq/20251214_vs_1218-focusedSwift-RTprimerconc_UMIsplintconc/"
results_base = os.path.join(base_data_dir, "results")
intermediate_dir = os.path.join(results_base, "intermediate")
reports_dir = os.path.join(results_base, "reports")
spot_checks_dir = os.path.join(reports_dir, "spot_checks")

# Load data
rates_df = pd.read_csv(os.path.join(intermediate_dir, "on_target_rates.csv"))
complexity_df = pd.read_csv(os.path.join(intermediate_dir, "on_target_complexity.csv"))

# Merge results on common columns to avoid suffixing
merged_df = pd.merge(rates_df, complexity_df, on=["Dataset", "Original_Name", "OnTarget_Raw", "OnTarget_Dedup"])

# Sort datasets
ds_order = [
    "400pmol_Primer_1uM_Splint",
    "40pmol_Primer_1uM_Splint",
    "10pmol_Primer_1uM_Splint_1214",
    "10pmol_Primer_1uM_Splint_1218"
]
merged_df['Dataset'] = pd.Categorical(merged_df['Dataset'], categories=ds_order, ordered=True)
final_df = merged_df.sort_values("Dataset")

# Load Specificity Stats
stats_df = pd.read_csv(stats_path)
acc_count = len(stats_df[stats_df['Status'] == 'Accepted'])
total_count = len(stats_df)

# Plotting Settings
sns.set_theme(style="whitegrid")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# Plot 1: On-Target Rate (Deduplicated)
sns.barplot(data=final_df, x="Dataset", y="OnTargetRateDedup_Percent", ax=ax1, palette="viridis", hue="Dataset", legend=False)
ax1.set_title("Strict On-Target Rate (Dedup BAM)", fontsize=14)
ax1.set_ylabel("On-Target Reads (%)")
ax1.tick_params(axis='x', rotation=30)

# Plot 2: Complexity (Ratio)
sns.barplot(data=final_df, x="Dataset", y="Complexity_NRF", ax=ax2, palette="magma", hue="Dataset", legend=False)
ax2.set_title("Library Complexity (Dedup / Raw On-Target Ratio)", fontsize=14)
ax2.set_ylabel("Complexity Ratio")
ax2.tick_params(axis='x', rotation=30)

plt.tight_layout()
plt.savefig(os.path.join(reports_dir, "comparative_analysis.png"), bbox_inches='tight')
plt.close()

# Generate Final Report
report_path = os.path.join(reports_dir, "final_summary_report.md")
with open(report_path, "w") as f:
    f.write("# Final Summary Report: Refined Focused-SWIFTseq Analysis\n\n")

    f.write("## 1. Logic Revision (Master Changelog)\n")
    f.write("> [!IMPORTANT]\n")
    f.write("> **Major Correction (Jan 15, 2:50 PM)**: Fixed strand-specific RT initiation direction and implemented strict 5'-end on-target definition.\n\n")
    f.write("- **Strict On-Target Definition**: A read is ONLY counted if its 5' end (initiation site) falls within the **1-15 nt initiation bin**.\n")
    f.write("- **Corrected RT Direction**: \n")
    f.write("  - `+` Genes: RT initiates at `Start-1` and moves upstream (LOWER coords). cDNA is `-` strand.\n")
    f.write("  - `-` Genes: RT initiates at `End+1` and moves downstream (HIGHER coords). cDNA is `+` strand.\n")
    f.write("- **Complexity (NRF)**: Ratio of On-Target reads in **Deduplicated vs Raw** BAMs. This measures the fraction of unique initiation events correctly captured.\n\n")

    f.write("## 2. Comparative Metrics\n")
    f.write(final_df[['Dataset', 'TotalReads_Raw', 'OnTarget_Raw', 'OnTarget_Dedup', 'OnTargetRateDedup_Percent', 'Complexity_NRF']].to_markdown(index=False) + "\n\n")

    f.write("## 3. Probe Specificity Audit\n")
    f.write(f"- **Accepted Probes**: {acc_count} / {total_count} ({acc_count/total_count:.1%})\n")
    f.write("- **Criteria**: Peak start counts in 1-15nt initiation bin > (Mode of ±500bp background counts + 2).\n\n")

    f.write("## 4. Visual Analysis\n")
    f.write("### Global Performance Summary\n")
    f.write("![Comparative Analysis Plots](./comparative_analysis.png)\n\n")

    f.write("### Complete Probe Gallery (61 Probes)\n")
    f.write("Comparison of 10pmol replicates: `10pmol_Primer_1uM_Splint_1214` vs `..._1218`.\n\n")

    if os.path.exists(manifest_path):
        plots_manifest = pd.read_csv(manifest_path)
        for _, prow in plots_manifest.iterrows():
            f.write(f"#### Probe {prow['ID']}: {prow['Gene']} (Status: {prow['Status']})\n")
            f.write(f"![Probe {prow['ID']}](./spot_checks/{prow['Filename']})\n\n")

print(f"Step 4 Complete. Final summary report saved to {report_path}")
