"""Locate guide RNAs in the reference genome."""

from __future__ import annotations

import itertools
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from .io_utils import write_fasta
from .models import GuideHit, GuideRecord, TargetRegion


def locate_guides_exact(
    guides: list[GuideRecord],
    blast_db: str,
    blastn: str,
) -> tuple[list[GuideHit], dict[str, int]]:
    """Find full-length exact guide matches with blastn-short."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        write_fasta([(g.guide_id, g.sequence) for g in guides], tmp_path)
        cmd = [
            blastn,
            "-task",
            "blastn-short",
            "-query",
            str(tmp_path),
            "-db",
            blast_db,
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand",
            "-word_size",
            "7",
            "-dust",
            "no",
            "-perc_identity",
            "100",
            "-qcov_hsp_perc",
            "100",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    finally:
        tmp_path.unlink(missing_ok=True)

    guide_by_id = {g.guide_id: g for g in guides}
    hits: list[GuideHit] = []
    partial_counts: dict[str, int] = {g.guide_id: 0 for g in guides}

    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        (
            qseqid,
            sseqid,
            pident,
            length,
            mismatch,
            gapopen,
            qstart,
            qend,
            sstart,
            send,
            evalue,
            bitscore,
            sstrand,
        ) = fields
        guide = guide_by_id[qseqid]
        length_int = int(length)
        is_full_exact = (
            float(pident) == 100.0
            and length_int == len(guide.sequence)
            and int(mismatch) == 0
            and int(gapopen) == 0
            and int(qstart) == 1
            and int(qend) == len(guide.sequence)
        )
        if not is_full_exact:
            partial_counts[qseqid] += 1
            continue
        start = min(int(sstart), int(send))
        end = max(int(sstart), int(send))
        hits.append(
            GuideHit(
                guide_id=qseqid,
                sequence=guide.sequence,
                chrom=sseqid,
                start=start,
                end=end,
                strand="+" if sstrand == "plus" else "-",
                pident=float(pident),
                length=length_int,
                evalue=evalue,
                bitscore=float(bitscore),
            )
        )

    missing = [g.guide_id for g in guides if not any(h.guide_id == g.guide_id for h in hits)]
    if missing:
        raise ValueError(f"No full-length exact mm10 hit found for: {', '.join(missing)}")

    return hits, partial_counts


def select_coherent_target(
    guides: list[GuideRecord],
    hits: list[GuideHit],
    max_cluster_span: int,
) -> TargetRegion:
    """Select the smallest same-chromosome guide hit cluster covering all guides."""
    hits_by_guide: dict[str, list[GuideHit]] = defaultdict(list)
    for hit in hits:
        hits_by_guide[hit.guide_id].append(hit)

    best_combo: tuple[GuideHit, ...] | None = None
    best_span: int | None = None
    guide_ids = [g.guide_id for g in guides]

    chroms = sorted({hit.chrom for hit in hits})
    for chrom in chroms:
        per_guide = []
        for guide_id in guide_ids:
            guide_hits = [h for h in hits_by_guide[guide_id] if h.chrom == chrom]
            if not guide_hits:
                break
            per_guide.append(guide_hits)
        else:
            for combo in itertools.product(*per_guide):
                start = min(hit.start for hit in combo)
                end = max(hit.end for hit in combo)
                span = end - start + 1
                if span > max_cluster_span:
                    continue
                if best_span is None or span < best_span:
                    best_span = span
                    best_combo = tuple(combo)

    if best_combo is None:
        raise ValueError(
            "Could not place all guides in one coherent genomic cluster. "
            f"Increase max_cluster_span if this is expected; current value is {max_cluster_span}."
        )

    return TargetRegion(
        chrom=best_combo[0].chrom,
        start=min(hit.start for hit in best_combo),
        end=max(hit.end for hit in best_combo),
        guide_hits=best_combo,
    )
