#!/usr/bin/env python3
"""
Add RTBC barcode to probe sequences.

This is a standalone post-processing step that can be applied to any
probe CSV file. Run this AFTER the main pipeline to add custom RTBC
barcodes without re-running the entire analysis.

Usage:
    python add_rtbc_barcode.py input.csv output.csv --rtbc "/5Phos/TGACTTGAGGAT"
"""

import argparse
import pandas as pd
from pathlib import Path


def add_rtbc_barcode(input_file: str, output_file: str, rtbc_sequence: str):
    """
    Add RTBC barcode to probe sequences.

    Args:
        input_file: Path to input CSV with probe sequences
        output_file: Path to output CSV with RTBC-modified sequences
        rtbc_sequence: RTBC barcode sequence (e.g., "/5Phos/TGACTTGAGGAT")
    """
    # Read input
    df = pd.read_csv(input_file)

    print(f"📄 Loaded {len(df)} probes from {input_file}")
    print(f"🧬 Adding RTBC barcode: {rtbc_sequence}")

    # Find the probe sequence column
    seq_col = None
    for col in ["Probe_Seq", "Seq", "Target_Seq"]:
        if col in df.columns:
            seq_col = col
            break

    if seq_col is None:
        raise ValueError("No sequence column found (Probe_Seq, Seq, or Target_Seq)")

    # Add RTBC to sequences
    df["RTBC_Sequence"] = rtbc_sequence
    df["Full_Oligo_Seq"] = rtbc_sequence + df[seq_col].astype(str)

    # Calculate new lengths
    df["Full_Oligo_Length"] = df["Full_Oligo_Seq"].str.len()

    # Save output
    df.to_csv(output_file, index=False)
    print(f"✅ Saved {len(df)} probes to {output_file}")

    # Summary
    print(f"\n📊 Summary:")
    print(f"   Original probe length: {df[seq_col].str.len().mean():.1f} nt")
    print(f"   RTBC length: {len(rtbc_sequence.replace('/5Phos/', ''))} nt")
    print(f"   Full oligo length: {df['Full_Oligo_Length'].mean():.1f} nt")

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Add RTBC barcode to probe sequences"
    )
    parser.add_argument("input", help="Input CSV file with probe sequences")
    parser.add_argument("output", help="Output CSV file with RTBC-modified sequences")
    parser.add_argument(
        "--rtbc",
        default="/5Phos/TGACTTGAGGAT",
        help="RTBC barcode sequence (default: /5Phos/TGACTTGAGGAT)"
    )

    args = parser.parse_args()

    add_rtbc_barcode(args.input, args.output, args.rtbc)


if __name__ == "__main__":
    main()
