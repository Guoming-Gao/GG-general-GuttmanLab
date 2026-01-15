import os
import pandas as pd
import argparse
import glob
import matplotlib.pyplot as plt
import numpy as np

def get_summary_stats(quant_dir):
    quant_files = glob.glob(os.path.join(quant_dir, "*_quant.csv"))
    data = []
    for f in quant_files:
        df = pd.read_csv(f)
        sample = os.path.basename(f).replace("_quant.csv", "")
        # Consistent with quantify_alleles.py: include low_conf in B6/Cast counts
        b6 = len(df[df['Allele'].str.startswith('B6', na=False)])
        cast = len(df[df['Allele'].str.startswith('Cast', na=False)])
        ratio = cast / (b6 + cast) if (b6 + cast) > 0 else 0
        data.append({
            "Sample": sample,
            "Total": len(df),
            "B6": b6,
            "Cast": cast,
            "Cast_Ratio": ratio
        })
    return pd.DataFrame(data)

def map_sample_to_condition(sample_name):
    """
    Standard conditions:
    1. WT (diff)
    2. WT_diffDox72h
    3. dTsix_diffDox72h
    4. dTsixdSPEN_diffDox72h
    """
    s = sample_name.lower()
    if "dtsixdspen" in s:
        return "dTsixdSPEN_Dox72h"
    if "dtsix" in s:
        return "dTsix_Dox72h"
    if "dox72h" in s:
        return "WT_Dox72h"
    if "diff" in s:
        return "WT_Untreated"
    return "Unknown"

def generate_stoich_plot(stoich_df, label, output_dir):
    """Generate a bar plot for SNP hit distribution."""
    # Filter for columns matching '[0-9]_SNP_Hits_%'
    hit_cols = [c for c in stoich_df.columns if "SNP_Hits_%" in c]
    if not hit_cols:
        return None

    # Sort columns numerically
    hit_cols.sort(key=lambda x: int(x.split("_")[0]))
    labels = [c.split("_")[0] for c in hit_cols]

    # Group by Condition and Allele
    # For simplicity, we aggregate across samples per condition
    # But since conditions are often 1 sample per rep, we just take the mean
    agg_df = stoich_df.groupby(['Condition', 'Allele'])[hit_cols].mean().reset_index()

    conditions = agg_df['Condition'].unique()
    num_conditions = len(conditions)

    fig, axes = plt.subplots(num_conditions, 1, figsize=(8, 4 * num_conditions), sharex=True)
    if num_conditions == 1:
        axes = [axes]

    for i, cond in enumerate(conditions):
        ax = axes[i]
        cond_df = agg_df[agg_df['Condition'] == cond]

        x = np.arange(len(labels))
        width = 0.35

        b6_vals = cond_df[cond_df['Allele'] == 'B6'][hit_cols].values
        cast_vals = cond_df[cond_df['Allele'] == 'Cast'][hit_cols].values

        if b6_vals.size > 0:
            ax.bar(x - width/2, b6_vals[0], width, label='B6', color='blue', alpha=0.6)
        if cast_vals.size > 0:
            ax.bar(x + width/2, cast_vals[0], width, label='Cast', color='red', alpha=0.6)

        ax.set_ylabel('Percentage of Reads (%)')
        ax.set_title(f'{label} - {cond}')
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend()
        ax.set_ylim(0, 100)
        ax.grid(axis='y', linestyle='--', alpha=0.7)

    plt.xlabel('Number of SNP Matches')
    plt.tight_layout()

    img_dir = os.path.join(output_dir, "images")
    os.makedirs(img_dir, exist_ok=True)
    plot_path = os.path.join(img_dir, f"{label}_stoichiometry.png")
    plt.savefig(plot_path)
    plt.close()
    return os.path.join("images", f"{label}_stoichiometry.png")

def main():
    parser = argparse.ArgumentParser(description="Compare multiple Xist datasets.")
    parser.add_argument("--results_dirs", required=True, help="Comma-separated results directories")
    parser.add_argument("--labels", required=True, help="Comma-separated labels for the datasets (e.g. ExonRep1,ExonRep2,Intron)")
    parser.add_argument("--output_dir", required=True, help="Directory to save comparison results")
    args = parser.parse_args()

    results_dirs = [d.strip() for d in args.results_dirs.split(",")]
    labels = [l.strip() for l in args.labels.split(",")]

    if len(results_dirs) != len(labels):
        print("Error: Number of results directories must match number of labels.")
        return

    os.makedirs(args.output_dir, exist_ok=True)

    summary_dfs = []
    stoich_data = {} # label -> dataframe

    for d, label in zip(results_dirs, labels):
        quant_dir = os.path.join(d, "quantification")
        stoich_file = os.path.join(d, "stoichiometry", "stoichiometry_summary.csv")

        if not os.path.exists(quant_dir):
            print(f"Warning: {quant_dir} not found. Skipping {label}.")
            continue

        # 1. Allele Stats
        df = get_summary_stats(quant_dir)
        if df.empty:
            print(f"Warning: No stats for {label}.")
            continue

        df['Condition'] = df['Sample'].apply(map_sample_to_condition)
        df_cond = df.groupby('Condition').agg({'Cast_Ratio': 'mean', 'Total': 'sum'}).reset_index()
        df_cond = df_cond.rename(columns={'Cast_Ratio': f'Cast_Ratio_{label}', 'Total': f'Total_{label}'})
        summary_dfs.append(df_cond)

        # 2. Stoichiometry
        if os.path.exists(stoich_file):
            s_df = pd.read_csv(stoich_file)
            s_df['Condition'] = s_df['Sample'].apply(map_sample_to_condition)
            stoich_data[label] = s_df

    if not summary_dfs:
        print("Error: No datasets to compare.")
        return

    # Merge all allele dataframes on Condition
    merged = summary_dfs[0]
    for next_df in summary_dfs[1:]:
        merged = pd.merge(merged, next_df, on='Condition', how='outer')

    order = ["WT_Untreated", "WT_Dox72h", "dTsix_Dox72h", "dTsixdSPEN_Dox72h"]
    merged['Condition'] = pd.Categorical(merged['Condition'], categories=order, ordered=True)
    merged = merged.sort_values('Condition')
    merged.to_csv(os.path.join(args.output_dir, "comparison_summary.csv"), index=False)

    # Generate Report
    report_path = os.path.join(args.output_dir, "Comparative_Analysis_Report.md")
    with open(report_path, "w") as f:
        f.write("# Xist Allele-Specific Expression Comparison\n\n")
        f.write(f"This report compares {len(labels)} datasets ({', '.join(labels)}) across 4 biological conditions.\n\n")

        f.write("## 1. Comparative Summary Table (Cast Ratios)\n\n")
        f.write(merged.to_markdown(index=False) + "\n\n")

        f.write("## 2. Stoichiometry & Allele Confidence\n")
        f.write("High confidence assignments (ALL_SNPS) indicate reads matching all targeted SNPs for an allele.\n\n")

        for label in labels:
            if label in stoich_data:
                f.write(f"### {label} Stoichiometry\n\n")
                # Show top level stoichiometry summary
                s_df = stoich_data[label]
                # Filter columns for display (Total_Reads, hit %, ALL_SNPS)
                cols = ['Condition', 'Allele', 'Total_Reads'] + [c for c in s_df.columns if "SNP_Hits_%" in c] + ['ALL_SNPS']
                display_df = s_df[cols].copy()
                display_df['Condition'] = pd.Categorical(display_df['Condition'], categories=order, ordered=True)
                display_df = display_df.sort_values(['Condition', 'Allele'])

                f.write(display_df.to_markdown(index=False) + "\n\n")

                # Bar Plot
                img_path = generate_stoich_plot(s_df, label, args.output_dir)
                if img_path:
                    f.write(f"![{label} Stoichiometry Plot]({img_path})\n\n")

        f.write("## 3. Condition Mapping Rationale\n")
        f.write("- **WT_Untreated**: Differentiated cells without Dox induction.\n")
        f.write("- **WT_Dox72h**: Differentiated cells with 72h Dox induction.\n")
        f.write("- **dTsix_Dox72h**: Tsix deletion with 72h Dox induction.\n")
        f.write("- **dTsixdSPEN_Dox72h**: Tsix and SPEN deletion with 72h Dox induction.\n\n")
        f.write("---\n")
        f.write("*Report generated for Guoming Gao (Guttman Lab).*")

    print(f"Comparison report generated at: {report_path}")

if __name__ == "__main__":
    main()
