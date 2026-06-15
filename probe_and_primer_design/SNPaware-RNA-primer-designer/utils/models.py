"""Shared data models for SNP-aware RNA primer design."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Exon:
    chrom: str
    start: int
    end: int
    exon_number: int | None = None

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class TranscriptModel:
    gene_name: str
    transcript_id: str
    chrom: str
    strand: str
    start: int
    end: int
    exons: tuple[Exon, ...]
    is_curated: bool
    exclusion_reason: str = ""

    @property
    def exon_count(self) -> int:
        return len(self.exons)

    @property
    def spliced_length(self) -> int:
        return sum(exon.length for exon in self.exons)


@dataclass(frozen=True)
class CdnaExonSegment:
    exon: Exon
    cdna_start: int
    cdna_end: int


@dataclass(frozen=True)
class CommonJunction:
    chrom: str
    strand: str
    upstream_exon_end_genomic: int
    downstream_exon_start_genomic: int
    cdna_left: int
    cdna_right: int
    upstream_exon_number: int | None
    downstream_exon_number: int | None
    transcript_count: int

    @property
    def key(self) -> tuple[str, str, int, int]:
        return (
            self.chrom,
            self.strand,
            self.upstream_exon_end_genomic,
            self.downstream_exon_start_genomic,
        )


@dataclass(frozen=True)
class TranscriptTemplate:
    gene_name: str
    transcript_id: str
    chrom: str
    strand: str
    sequence: str
    segments: tuple[CdnaExonSegment, ...]

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class SnpRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    b6_gt: str
    cast_gt: str

    @property
    def genotype(self) -> str:
        return f"{self.b6_gt},{self.cast_gt}"


@dataclass(frozen=True)
class SpecificityResult:
    sequence: str
    perfect_or_near_hits: int
    noise_hits: int
    pass_specificity: bool
