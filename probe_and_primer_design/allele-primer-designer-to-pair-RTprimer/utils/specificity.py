"""Short-primer BLAST specificity checks."""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from .io_utils import write_fasta
from .models import SpecificityResult


def check_primer_specificity(
    sequence: str,
    blast_db: str,
    blastn: str,
    config: dict,
) -> SpecificityResult:
    """Count strong BLAST hits for one primer sequence."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        write_fasta([("primer", sequence)], tmp_path)
        cmd = [
            blastn,
            "-task",
            "blastn-short",
            "-query",
            str(tmp_path),
            "-db",
            blast_db,
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-word_size",
            "7",
            "-dust",
            "no",
            "-evalue",
            "10",
            "-max_target_seqs",
            "100",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    finally:
        tmp_path.unlink(missing_ok=True)

    perfect_or_near = 0
    noise = 0
    unique_len = min(int(config["blast_min_unique_len"]), len(sequence))
    unique_identity = float(config["blast_min_unique_identity"])
    noise_len = int(config["blast_noise_len"])

    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        pident = float(fields[2])
        aln_len = int(fields[3])
        if aln_len >= unique_len and pident >= unique_identity:
            perfect_or_near += 1
        elif aln_len >= noise_len:
            noise += 1

    return SpecificityResult(
        sequence=sequence,
        perfect_or_near_hits=perfect_or_near,
        noise_hits=noise,
        pass_specificity=(perfect_or_near == 1 and noise == 0),
    )


def check_primer_specificity_batch(
    sequences: list[str],
    blast_db: str,
    blastn: str,
    config: dict,
) -> list[SpecificityResult]:
    """Check multiple primers in one BLAST process, preserving input order."""
    if not sequences:
        return []
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tmp:
        tmp_path = Path(tmp.name)
    try:
        write_fasta([(f"primer_{index}", sequence) for index, sequence in enumerate(sequences)], tmp_path)
        cmd = [
            blastn,
            "-task", "blastn-short",
            "-query", str(tmp_path),
            "-db", blast_db,
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-word_size", "7",
            "-dust", "no",
            "-evalue", "10",
            "-max_target_seqs", "100",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    finally:
        tmp_path.unlink(missing_ok=True)

    counts = [[0, 0] for _ in sequences]
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        index = int(fields[0].removeprefix("primer_"))
        pident = float(fields[2])
        aln_len = int(fields[3])
        unique_len = min(int(config["blast_min_unique_len"]), len(sequences[index]))
        if aln_len >= unique_len and pident >= float(config["blast_min_unique_identity"]):
            counts[index][0] += 1
        elif aln_len >= int(config["blast_noise_len"]):
            counts[index][1] += 1

    return [
        SpecificityResult(sequence, perfect, noise, perfect == 1 and noise == 0)
        for sequence, (perfect, noise) in zip(sequences, counts)
    ]
