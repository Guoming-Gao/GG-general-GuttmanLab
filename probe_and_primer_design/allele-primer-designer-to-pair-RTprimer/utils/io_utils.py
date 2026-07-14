"""Input and output helpers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


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
