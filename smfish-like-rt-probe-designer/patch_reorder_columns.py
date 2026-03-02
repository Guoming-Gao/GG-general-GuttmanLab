#!/usr/bin/env python3
"""
patch_reorder_columns.py - Reorder columns and discard unexpected ones to match the current pipeline format.

Usage:
    python patch_reorder_columns.py input.csv

This will generate input-reordered.csv with exactly the ordered list of known columns,
and discard any extra columns that are not recognized.
"""

import argparse
import pandas as pd
from pathlib import Path
import sys

# The exact desired ordered list of columns for the pipeline
TARGET_COLUMNS_ORDER = [
    # Top 9 priority columns requested
    "ProbeID",
    "GeneName",
    "Chromosome",
    "RegionType",
    "SNP_Count",
    "Probe_Start",
    "Probe_End",
    "Probe_Seq",
    "Full_Oligo_Seq",

    # Original candidate probe columns
    "Species",
    "Target_Seq",
    "Probe_Length",
    "Target_Strand",
    "SNPs_Positions",
    "SNPs_Types",
    "RT_Region_Start",
    "RT_Region_End",
    "RT_Product_Seq",
    "dG37",
    "dGScore",
    "dGOpt",
    "GC_Content",
    "GCFilter",
    "PNASFilter",
    "aCompFilter",
    "aStackFilter",
    "cCompFilter",
    "cStackFilter",
    "cSpecStackFilter",
    "NbOfPNAS",
    "MaskedFilter",
    "RepeatMaskerPC",

    # BLAST Specificity validation columns
    "BLAST_Hits",
    "Primary_Identity",
    "Secondary_Identity",
    "BLAST_Unique",
    "BLAST_WordSize",
    "BLAST_EValue",
    "BLAST_MinAlignment",

    # RTBC Barcode columns
    "RTBC_Sequence",
    "Full_Oligo_Length",

    # Optional Forward Primer columns
    "Forward_Primer_Seq",
    "Forward_Primer_Tm",
    "Forward_Primer_GC",
    "PCR_Amplicon_Size",
    "RT_Region_Covered"
]

def format_csv(input_csv):
    """
    Load CSV, reorder columns, and discard any extra columns not in the known list.
    """
    in_path = Path(input_csv)
    if not in_path.exists():
        print(f"Error: Could not find input file '{input_csv}'")
        sys.exit(1)

    print(f"📄 Loading '{in_path.name}'...")
    df = pd.read_csv(in_path)
    original_cols = list(df.columns)

    # Find columns that are present in both the CSV and the target order
    kept_columns = [col for col in TARGET_COLUMNS_ORDER if col in df.columns]

    # Identify unrecognized columns to discard
    discarded_columns = [col for col in original_cols if col not in kept_columns]

    if discarded_columns:
        print(f"⚠️  Discarding unrecognized or extra columns: {', '.join(discarded_columns)}")
    else:
        print("✅ No unrecognized columns to discard.")

    # Reorder the dataframe
    df_reordered = df[kept_columns]

    # Construct output filename
    out_path = in_path.parent / f"{in_path.stem}-reordered.csv"

    # Save the file
    df_reordered.to_csv(out_path, index=False)
    print(f"✅ Generated properly formatted file at: {out_path}")

def main():
    parser = argparse.ArgumentParser(description="Reorder and filter columns of old FISH-RT probe pipeline CSVs.")
    parser.add_argument("input_csv", help="Input CSV file to process")
    args = parser.parse_args()

    format_csv(args.input_csv)

if __name__ == "__main__":
    main()
