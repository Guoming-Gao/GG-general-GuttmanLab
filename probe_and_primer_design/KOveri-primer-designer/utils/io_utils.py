"""Input and output helpers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .models import GuideRecord


def read_guide_fasta(path: str | Path) -> list[GuideRecord]:
    """Read gRNA records from FASTA without changing sequence orientation."""
    records: list[GuideRecord] = []
    current_id: str | None = None
    current_seq: list[str] = []
    path = Path(path)

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append(
                        GuideRecord(current_id, "".join(current_seq).upper().replace("U", "T"))
                    )
                current_id = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_id is not None:
        records.append(GuideRecord(current_id, "".join(current_seq).upper().replace("U", "T")))

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    for record in records:
        invalid = sorted(set(record.sequence) - set("ACGTN"))
        if invalid:
            raise ValueError(f"{record.guide_id} contains unsupported bases: {''.join(invalid)}")

    return records


def write_fasta(records: list[tuple[str, str]], path: str | Path) -> None:
    path = Path(path)
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def write_dataframe(rows: list[dict], path: str | Path) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    return df
