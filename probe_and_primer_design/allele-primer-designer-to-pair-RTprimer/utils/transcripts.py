"""Transcript annotation and spliced cDNA helpers."""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path

from .models import CdnaExonSegment, CommonJunction, Exon, TranscriptModel, TranscriptTemplate
from .sequence_utils import extract_reference_sequence, reverse_complement


def load_gene_transcripts(
    gtf_path: str | Path,
    gene_name: str,
    curated_prefixes: tuple[str, ...],
) -> tuple[list[TranscriptModel], list[TranscriptModel]]:
    """Load all transcript models for one gene and return all plus designable curated models."""
    transcript_meta: dict[str, dict] = {}
    exons: dict[str, list[Exon]] = defaultdict(list)
    gtf_path = Path(gtf_path)
    opener = gzip.open if gtf_path.suffix == ".gz" else open

    with opener(gtf_path, "rt") as handle:
        for raw_line in handle:
            if raw_line.startswith("#") or not raw_line.strip():
                continue
            fields = raw_line.rstrip().split("\t")
            if len(fields) < 9:
                continue
            chrom, _source, feature, start, end, _score, strand, _frame, attrs = fields
            attr = _parse_attributes(attrs)
            if attr.get("gene_name") != gene_name and attr.get("gene_id") != gene_name:
                continue
            transcript_id = attr.get("transcript_id")
            if not transcript_id:
                continue
            transcript_meta.setdefault(
                transcript_id,
                {
                    "gene_name": attr.get("gene_name", gene_name),
                    "chrom": chrom,
                    "strand": strand,
                    "start": int(start),
                    "end": int(end),
                },
            )
            transcript_meta[transcript_id]["start"] = min(transcript_meta[transcript_id]["start"], int(start))
            transcript_meta[transcript_id]["end"] = max(transcript_meta[transcript_id]["end"], int(end))
            if feature == "exon":
                exon_number = _parse_int(attr.get("exon_number"))
                exons[transcript_id].append(Exon(chrom, int(start), int(end), exon_number))

    transcripts: list[TranscriptModel] = []
    for transcript_id, meta in transcript_meta.items():
        exon_tuple = tuple(sorted(exons.get(transcript_id, []), key=lambda exon: (exon.start, exon.end)))
        is_curated = transcript_id.startswith(curated_prefixes)
        reason = ""
        if not is_curated:
            reason = "not_curated_refseq"
        elif len(exon_tuple) < 2:
            reason = "single_exon_or_no_exons"
        transcripts.append(
            TranscriptModel(
                gene_name=meta["gene_name"],
                transcript_id=transcript_id,
                chrom=meta["chrom"],
                strand=meta["strand"],
                start=meta["start"],
                end=meta["end"],
                exons=exon_tuple,
                is_curated=is_curated and len(exon_tuple) >= 2,
                exclusion_reason=reason,
            )
        )

    transcripts.sort(key=lambda tx: (not tx.is_curated, tx.transcript_id))
    curated = [tx for tx in transcripts if tx.is_curated]
    return transcripts, curated


def select_template_transcript(curated: list[TranscriptModel]) -> TranscriptModel:
    if not curated:
        raise ValueError("No curated multi-exon NM_/NR_ transcript models found.")
    return max(curated, key=lambda tx: (tx.spliced_length, tx.exon_count, tx.transcript_id))


def build_transcript_template(
    transcript: TranscriptModel,
    genome_fasta: str,
    samtools: str,
) -> TranscriptTemplate:
    """Build a strand-aware spliced cDNA template for one transcript."""
    pieces: list[str] = []
    segments: list[CdnaExonSegment] = []
    cdna_pos = 1
    for exon in transcript_order_exons(transcript):
        seq = extract_reference_sequence(exon.chrom, exon.start, exon.end, genome_fasta, samtools)
        if transcript.strand == "-":
            seq = reverse_complement(seq)
        cdna_start = cdna_pos
        cdna_end = cdna_pos + len(seq) - 1
        pieces.append(seq)
        segments.append(CdnaExonSegment(exon, cdna_start, cdna_end))
        cdna_pos = cdna_end + 1
    return TranscriptTemplate(
        gene_name=transcript.gene_name,
        transcript_id=transcript.transcript_id,
        chrom=transcript.chrom,
        strand=transcript.strand,
        sequence="".join(pieces),
        segments=tuple(segments),
    )


def common_junctions_for_template(
    curated: list[TranscriptModel],
    template_transcript: TranscriptModel,
    template: TranscriptTemplate,
) -> list[CommonJunction]:
    """Return junctions shared by every curated model and present in the selected template."""
    if not curated:
        return []
    per_tx_keys = [set(_junction_keys(tx)) for tx in curated]
    common_keys = set.intersection(*per_tx_keys) if per_tx_keys else set()
    if not common_keys:
        return []

    segment_by_exon = {segment.exon: segment for segment in template.segments}
    junctions: list[CommonJunction] = []
    ordered = transcript_order_exons(template_transcript)
    for prev_exon, next_exon in zip(ordered, ordered[1:]):
        key = _junction_key(template_transcript.strand, prev_exon, next_exon)
        if key not in common_keys:
            continue
        prev_segment = segment_by_exon[prev_exon]
        next_segment = segment_by_exon[next_exon]
        junctions.append(
            CommonJunction(
                chrom=template_transcript.chrom,
                strand=template_transcript.strand,
                upstream_exon_end_genomic=key[2],
                downstream_exon_start_genomic=key[3],
                cdna_left=prev_segment.cdna_end,
                cdna_right=next_segment.cdna_start,
                upstream_exon_number=prev_exon.exon_number,
                downstream_exon_number=next_exon.exon_number,
                transcript_count=len(curated),
            )
        )
    return junctions


def transcript_order_exons(transcript: TranscriptModel) -> list[Exon]:
    return sorted(transcript.exons, key=lambda exon: exon.start, reverse=(transcript.strand == "-"))


def map_cdna_pos_to_genomic(template: TranscriptTemplate, cdna_pos: int) -> int | None:
    for segment in template.segments:
        if segment.cdna_start <= cdna_pos <= segment.cdna_end:
            offset = cdna_pos - segment.cdna_start
            if template.strand == "+":
                return segment.exon.start + offset
            return segment.exon.end - offset
    return None


def map_genomic_pos_to_cdna(template: TranscriptTemplate, genomic_pos: int) -> int | None:
    for segment in template.segments:
        exon = segment.exon
        if exon.start <= genomic_pos <= exon.end:
            if template.strand == "+":
                return segment.cdna_start + genomic_pos - exon.start
            return segment.cdna_start + exon.end - genomic_pos
    return None


def cdna_interval_to_genomic_segments(
    template: TranscriptTemplate,
    cdna_start: int,
    cdna_end: int,
) -> list[tuple[str, int, int, str]]:
    segments = []
    for segment in template.segments:
        overlap_start = max(cdna_start, segment.cdna_start)
        overlap_end = min(cdna_end, segment.cdna_end)
        if overlap_end < overlap_start:
            continue
        start_genomic = map_cdna_pos_to_genomic(template, overlap_start)
        end_genomic = map_cdna_pos_to_genomic(template, overlap_end)
        if start_genomic is None or end_genomic is None:
            continue
        segments.append((segment.exon.chrom, min(start_genomic, end_genomic), max(start_genomic, end_genomic), template.strand))
    return segments


def format_segments(segments: list[tuple[str, int, int, str]]) -> str:
    return ";".join(f"{chrom}:{start}-{end}:{strand}" for chrom, start, end, strand in segments)


def _junction_keys(transcript: TranscriptModel) -> list[tuple[str, str, int, int]]:
    ordered = transcript_order_exons(transcript)
    return [_junction_key(transcript.strand, prev_exon, next_exon) for prev_exon, next_exon in zip(ordered, ordered[1:])]


def _junction_key(strand: str, prev_exon: Exon, next_exon: Exon) -> tuple[str, str, int, int]:
    if strand == "+":
        return (prev_exon.chrom, strand, prev_exon.end, next_exon.start)
    return (prev_exon.chrom, strand, prev_exon.start, next_exon.end)


def _parse_attributes(attributes: str) -> dict[str, str]:
    parsed = {}
    for key, value in re.findall(r'(\S+)\s+"([^"]*)"', attributes):
        parsed[key] = value
    return parsed


def _parse_int(value: str | None) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except ValueError:
        return None
