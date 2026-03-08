
import pandas as pd
import os

# Paths
results_base = "/Volumes/guttman/users/gmgao/Data_seq/20251214_vs_1218-focusedSwift-RTprimerconc_UMIsplintconc/results"
processed_probes_path = os.path.join(results_base, "intermediate/processed_probes.csv")
reports_dir = os.path.join(results_base, "reports")
output_dir = os.path.join(results_base, "intermediate")

rates_path = os.path.join(output_dir, "on_target_rates.csv")
output_path = os.path.join(output_dir, "on_target_complexity.csv")

# Load rates calculated in Step 2
rates_df = pd.read_csv(rates_path)

# Complexity = OnTarget_Dedup / OnTarget_Raw
# This represents the fraction of unique biological fragments that were sequenced.
rates_df['Complexity_NRF'] = (rates_df['OnTarget_Dedup'] / rates_df['OnTarget_Raw']).fillna(0)

# Select relevant columns for complexity report
complexity_df = rates_df[['Dataset', 'Original_Name', 'OnTarget_Raw', 'OnTarget_Dedup', 'Complexity_NRF']]

complexity_df.to_csv(output_path, index=False)

print(f"Step 3 (Complexity Ratio) Complete. Saved to {output_path}")
print(complexity_df[['Dataset', 'Complexity_NRF']].to_string(index=False))
