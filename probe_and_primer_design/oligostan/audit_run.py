"""Independently audit an exported mouse smiFISH run before ordering."""

import argparse
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

from .mouse_smifish import FastaIndex, GENES, load_models, load_run_config, quality_tier


def audit(root, config_path=None):
    root = Path(root)
    config = load_run_config(config_path)
    selected = pd.read_csv(root / "all_selected_blast_verified.csv")
    hits = pd.read_csv(root / "blast_hits.tsv", sep="\t")
    summary = pd.read_csv(root / "set_summary.csv")
    reference = FastaIndex(config["fasta"])
    models = load_models(config["gtf"])
    assert len(summary) == 12
    assert set(summary.gene) == set(GENES)
    assert set(summary.region) == {"exon", "intron"}
    assert selected.probe_id.is_unique and selected.probe_seq.is_unique
    assert selected.blast_verified.all()
    assert selected.blast_target_exact.all()
    assert not selected.blast_strong_offtarget.any()
    for row in selected.itertuples(index=False):
        plus = reference.fetch(row.chrom, int(row.start), int(row.end))
        sense = plus if row.target_strand == "+" else str(Seq(plus).reverse_complement())
        assert row.target_seq == sense, row.probe_id
        assert row.probe_seq == str(Seq(sense).reverse_complement()), row.probe_id
        assert any(row.start >= start and row.end <= end for start, end in models[row.gene][row.region]), row.probe_id
        assert row.quality_tier == quality_tier(row._asdict()), row.probe_id
        intended = hits[(hits.probe_id == row.probe_id) & hits.expected_locus & hits.full_exact]
        assert not intended.empty, row.probe_id
        assert hits[(hits.probe_id == row.probe_id) & hits.strong_offtarget].empty, row.probe_id
    for s in summary.itertuples(index=False):
        table = pd.read_csv(root / s.gene / s.region / "selected_blast_verified.csv")
        assert len(table) == s.selected_count
        assert s.selected_count == int(((selected.gene == s.gene) & (selected.region == s.region)).sum())
        assert s.meets_minimum == (s.selected_count >= s.minimum_requested)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--config", type=Path, help="Local JSON config (default: sibling .smifish-local.json)")
    args = parser.parse_args()
    summary = audit(args.run_dir, config_path=args.config)
    print(summary.to_string(index=False))
    print("Sequence, coordinate, BLAST and export audit passed")


if __name__ == "__main__":
    main()
