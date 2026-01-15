
import pandas as pd
import gzip
import os

# Input paths
probe_csv_path = "/Volumes/guttman/users/gmgao/Data_seq/20251214_vs_1218-focusedSwift-RTprimerconc_UMIsplintconc/FISH_RT_probes-final-reduced_3pergene_probeswithRT.csv"
gtf_path = "/Volumes/guttman/genomes/mm10/annotation/mm10.refGene.gtf.gz"
# Result paths (Centralized in dataset directory)
results_base = os.path.join(os.path.dirname(probe_csv_path), "results")
output_dir = os.path.join(results_base, "intermediate")
reports_dir = os.path.join(results_base, "reports")

os.makedirs(output_dir, exist_ok=True)
os.makedirs(reports_dir, exist_ok=True)

# Load probes
df = pd.read_csv(probe_csv_path)
genes = df['GeneName'].unique()

# Extract strands from GTF - Robustly
gene_strands = {}
with gzip.open(gtf_path, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.split('\t')
        if parts[2] == 'gene' or parts[2] == 'transcript':
            info = parts[8]
            if 'gene_name "' in info:
                gene_name = info.split('gene_name "')[1].split('"')[0]
                if gene_name in genes:
                    strand = parts[6]
                    # Prioritize finding the transcript/gene strand
                    gene_strands[gene_name] = strand

# Merge strand info
df['Strand'] = df['GeneName'].map(gene_strands)

# Calculate Windows (Corrected RT Direction)
# + strand genes: RT moves high -> low. Initiation Bin: [Start-15, Start-1]
# - strand genes: RT moves low -> high. Initiation Bin: [End+1, End+15]
def calculate_window(row):
    if row['Strand'] == '+':
        return row['theStartPos'] - 15, row['theStartPos'] - 1
    elif row['Strand'] == '-':
        return row['theEndPos'] + 1, row['theEndPos'] + 15
    else:
        return None, None

df['WindowStart'], df['WindowEnd'] = zip(*df.apply(calculate_window, axis=1))

# Save processed probes
processed_probe_path = os.path.join(output_dir, "processed_probes.csv")
df.to_csv(processed_probe_path, index=False)

# Validation Report
report_path = os.path.join(reports_dir, "step1_validation.md")
with open(report_path, "w") as f:
    f.write("# Step 1 Validation Report: Reference Annotation & Probe Processing\n\n")
    f.write("## Strand Assignment\n")
    f.write(df[['GeneName', 'Strand']].drop_duplicates().to_markdown(index=False) + "\n\n")
    f.write("## On-Target Window Calculation (Sample)\n")
    f.write(df[['GeneName', 'Chromosome', 'theStartPos', 'theEndPos', 'Strand', 'WindowStart', 'WindowEnd']].head(10).to_markdown(index=False) + "\n\n")

print(f"Step 1 Complete. Processed probes saved to {processed_probe_path}")
