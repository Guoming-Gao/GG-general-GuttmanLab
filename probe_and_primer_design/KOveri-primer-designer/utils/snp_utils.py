"""B6/Cast SNP querying."""

from __future__ import annotations

from pathlib import Path

import pysam

from .models import SnpRecord


class SnpDatabase:
    """On-demand tabix queries for B6/Cast informative SNPs."""

    def __init__(self, vcf_path: str, b6_sample: str, cast_sample: str):
        if not Path(vcf_path).exists():
            raise FileNotFoundError(f"VCF not found: {vcf_path}")
        if not Path(f"{vcf_path}.tbi").exists():
            raise FileNotFoundError(f"Tabix index not found: {vcf_path}.tbi")

        self.vcf_path = vcf_path
        self.b6_sample = b6_sample
        self.cast_sample = cast_sample
        self.vcf = pysam.TabixFile(vcf_path)
        self.samples = self._samples_from_header()
        self.b6_idx = self.samples.index(b6_sample)
        self.cast_idx = self.samples.index(cast_sample)

    def _samples_from_header(self) -> list[str]:
        for line in self.vcf.header:
            text = line.decode() if isinstance(line, bytes) else line
            if text.startswith("#CHROM"):
                return text.strip().split("\t")[9:]
        raise ValueError("No #CHROM header line found in VCF")

    def query_informative(self, chrom: str, start: int, end: int) -> list[SnpRecord]:
        """Return SNPs where B6 and Cast genotypes differ."""
        vcf_chrom = chrom[3:] if chrom.startswith("chr") else chrom
        rows: list[SnpRecord] = []
        try:
            iterator = self.vcf.fetch(vcf_chrom, start - 1, end)
        except ValueError:
            return rows

        for line in iterator:
            text = line.decode() if isinstance(line, bytes) else line
            fields = text.rstrip().split("\t")
            if len(fields) < 10:
                continue
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            if len(ref) != 1:
                continue
            sample_fields = fields[9:]
            b6_gt = self._parse_gt(sample_fields[self.b6_idx], ref, alt)
            cast_gt = self._parse_gt(sample_fields[self.cast_idx], ref, alt)
            if b6_gt is None or cast_gt is None or b6_gt == cast_gt:
                continue
            rows.append(
                SnpRecord(
                    chrom=f"chr{vcf_chrom}",
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    b6_gt=b6_gt,
                    cast_gt=cast_gt,
                )
            )
        return rows

    @staticmethod
    def _parse_gt(sample_data: str, ref: str, alt: str) -> str | None:
        gt_field = sample_data.split(":")[0]
        if "." in gt_field:
            return None
        alleles = [ref] + alt.split(",")
        try:
            return "/".join(alleles[int(i)] for i in gt_field.replace("|", "/").split("/"))
        except (ValueError, IndexError):
            return None

    def close(self) -> None:
        self.vcf.close()
