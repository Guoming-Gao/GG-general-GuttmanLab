"""Pure coordinate rules for allele-primer design."""

from __future__ import annotations


def primer_to_snp_distance_toward_polydt(primer_cdna_end: int, snp_cdna: int, window: int) -> int | None:
    """Return positive primer-3'-end-to-SNP distance only toward the RNA 3' end."""
    distance = snp_cdna - primer_cdna_end
    if 1 <= distance <= window:
        return distance
    return None


def snp_distance_to_rna_3p(template_length: int, snp_cdna: int) -> int:
    """Distance in spliced RNA/cDNA bases from SNP to the transcript 3' end."""
    return template_length - snp_cdna

