import pandas as pd
import argparse
import os

def convert_csv_to_bed(input_csv, output_bed, score_col=None):
    """
    Converts a probe design CSV to an extended BED format file (BED6+).
    Adding a browser track line to help IGV display the extra columns if possible,
    but IGV's support for arbitrary columns in BED files is limited.
    """
    if not os.path.exists(input_csv):
        raise FileNotFoundError(f"Input file not found: {input_csv}")

    print(f"Reading {input_csv}...")
    df = pd.read_csv(input_csv)

    # BED standard: 0-based start, half-open
    df['chrom'] = df['Chromosome']
    df['chromStart'] = df['Probe_Start'] - 1
    df['chromEnd'] = df['Probe_End']
    df['name'] = df['ProbeID']

    # Strand should be opposite of Target_Strand
    strand_map = {'+': '-', '-': '+'}
    df['strand'] = df['Target_Strand'].map(strand_map).fillna('.')

    if score_col and score_col in df.columns:
        df['score'] = df[score_col]
    else:
        df['score'] = 0

    # Standard BED6 columns
    bed_cols = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']

    # BED files don't naturally support arbitrary named metadata columns in standard viewers like IGV.
    # To see metadata in standard BED, we can pack it into the 'name' field,
    # or create a custom format. Let's create an extended name field for IGV.

    # Example: Name | GeneName | RegionType | SNP_Count
    # igv will show this single string when you hover/click.

    df['extended_name'] = df['ProbeID'] + "|" + df['GeneName'] + "|" + df['RegionType'] + "|SNPs:" + df['SNP_Count'].astype(str)

    # Let's replace the standard 'name' with our extended name
    df['name'] = df['extended_name']

    # Write just the BED6 columns, because IGV mostly ignores anything after column 12
    # and doesn't display arbitrary columns well.
    bed_df = df[bed_cols]

    # Write the BED file
    with open(output_bed, 'w') as f:
        # Optional: Add a track line for IGV
        f.write('track name="FISH-like RT Probes" description="smFISH Probes" useScore=1\n')

    bed_df.to_csv(output_bed, sep='\t', header=False, index=False, mode='a')
    print(f"Successfully wrote {len(bed_df)} records to {output_bed}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert probe sequences CSV to BED format.")
    parser.add_argument("-i", "--input", required=True, help="Input CSV file path")
    parser.add_argument("-o", "--output", required=True, help="Output BED file path")
    parser.add_argument("-s", "--score_col", default=None, help="Optional CSV column to use as the BED score (e.g., dGScore, SNP_Count)")

    args = parser.parse_args()

    convert_csv_to_bed(args.input, args.output, args.score_col)
