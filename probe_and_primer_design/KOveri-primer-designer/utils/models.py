"""Shared data models."""

from dataclasses import dataclass


@dataclass(frozen=True)
class GuideRecord:
    guide_id: str
    sequence: str


@dataclass(frozen=True)
class GuideHit:
    guide_id: str
    sequence: str
    chrom: str
    start: int
    end: int
    strand: str
    pident: float
    length: int
    evalue: str
    bitscore: float


@dataclass(frozen=True)
class TargetRegion:
    chrom: str
    start: int
    end: int
    guide_hits: tuple[GuideHit, ...]

    @property
    def length(self) -> int:
        return self.end - self.start + 1


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
