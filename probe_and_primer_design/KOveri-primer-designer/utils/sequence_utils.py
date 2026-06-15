"""Reference sequence helpers."""

from __future__ import annotations

import subprocess


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1].upper()


def extract_reference_sequence(
    chrom: str,
    start: int,
    end: int,
    genome_fasta: str,
    samtools: str,
) -> str:
    """Extract plus-strand genomic sequence as 5' to 3' increasing coordinates."""
    if start < 1 or end < start:
        raise ValueError(f"Invalid genomic interval: {chrom}:{start}-{end}")
    region = f"{chrom}:{start}-{end}"
    result = subprocess.run(
        [samtools, "faidx", genome_fasta, region],
        capture_output=True,
        text=True,
        check=True,
    )
    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError(f"No sequence returned for {region}")
    return "".join(lines[1:]).upper()
