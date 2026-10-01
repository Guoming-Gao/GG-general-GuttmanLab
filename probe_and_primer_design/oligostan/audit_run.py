"""Independently audit an exported mouse smiFISH run before ordering."""

import argparse
import json
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from .filters import is_ok_4_homopolymer, longest_homopolymer_run
from .mouse_smifish import FastaIndex, load_models, load_run_config, quality_tier
from .minimum_span import expected_sets, select_minimum_span_sets


def audit_homopolymers(frame, max_len, table_name):
    """Reject a table containing probes beyond the recorded run-length limit."""
    if "max_homopolymer_run" not in frame:
        raise AssertionError(f"{table_name} lacks max_homopolymer_run")
    for row in frame.itertuples(index=False):
        actual = longest_homopolymer_run(row.probe_seq)
        if (actual != row.max_homopolymer_run or actual > max_len or
                not is_ok_4_homopolymer(row.probe_seq, max_len)):
            raise AssertionError(f"{table_name}: homopolymer check failed for {row.probe_id}")


def audit(root, config_path=None):
    root = Path(root)
    config = load_run_config(config_path)
    selected = pd.read_csv(root / "all_selected_blast_verified.csv")
    candidates = pd.read_csv(root / "all_design_candidates.csv")
    verified = pd.read_csv(root / "all_candidates_blast_status.csv")
    hits = pd.read_csv(root / "blast_hits.tsv", sep="\t")
    summary = pd.read_csv(root / "set_summary.csv")
    manifest = json.loads((root / "manifest.json").read_text())
    max_homopolymer_length = int(manifest["max_homopolymer_length"])
    assert max_homopolymer_length >= 1
    for name, frame in (("candidates", candidates), ("BLAST status", verified),
                        ("selected", selected)):
        audit_homopolymers(frame, max_homopolymer_length, name)
    assert set(candidates.probe_id) == set(verified.probe_id)
    reference = FastaIndex(config["fasta"])
    models = load_models(config["gtf"])
    assert set(zip(summary.gene, summary.region)) == set(expected_sets(models))
    assert len(summary) == len(expected_sets(models))
    assert selected.probe_id.is_unique and selected.probe_seq.is_unique
    assert selected.blast_verified.all()
    assert selected.blast_target_exact.all()
    assert not selected.blast_strong_offtarget.any()
    exact_target_ids = set(hits.loc[hits.expected_locus & hits.full_exact, "probe_id"])
    strong_offtarget_ids = set(hits.loc[hits.strong_offtarget, "probe_id"])
    expected_selected, expected_summary = select_minimum_span_sets(verified, models, probes_per_set=30)
    assert set(selected.probe_id) == set(expected_selected.probe_id)
    for col in ("selected_count", "blast_tested_count", "blast_verified_count", "genomic_span_bp"):
        assert summary[col].fillna(-1).tolist() == expected_summary[col].fillna(-1).tolist(), col
    for row in selected.itertuples(index=False):
        plus = reference.fetch(row.chrom, int(row.start), int(row.end))
        sense = plus if row.target_strand == "+" else str(Seq(plus).reverse_complement())
        assert row.target_seq == sense, row.probe_id
        assert row.probe_seq == str(Seq(sense).reverse_complement()), row.probe_id
        assert any(row.start >= start and row.end <= end for start, end in models[row.gene][row.region]), row.probe_id
        assert row.quality_tier == quality_tier(row._asdict()), row.probe_id
        assert row.probe_id in exact_target_ids, row.probe_id
        assert row.probe_id not in strong_offtarget_ids, row.probe_id
    for s in summary.itertuples(index=False):
        table = pd.read_csv(root / s.gene / s.region / "selected_blast_verified.csv")
        assert len(table) == s.selected_count
        assert s.selected_count == int(((selected.gene == s.gene) & (selected.region == s.region)).sum())
        assert s.meets_minimum == (s.selected_count >= s.minimum_requested)
        order = pd.read_csv(root / s.gene / s.region / "order_FlapX_REVIEW.csv")
        assert order.Sequence.tolist() == table.sort_values("start").HybFlpX.tolist()
        fasta = list(SeqIO.parse(root / s.gene / s.region / "selected_probes.fa", "fasta"))
        assert len(fasta) == s.selected_count
        if s.selected_count:
            subset = selected[(selected.gene == s.gene) & (selected.region == s.region)]
            assert s.genomic_span_bp == int(subset.end.max() - subset.start.min() + 1)
    coverage = pd.read_csv(root / "coverage_summary.csv")
    assert len(coverage) == len(summary)
    assert coverage["genomic_span_bp"].fillna(-1).tolist() == summary["genomic_span_bp"].fillna(-1).tolist()
    assert (root / "coverage_report.pdf").stat().st_size > 0
    assert len(list((root / "coverage_plots").glob("*.png"))) == len(summary)
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
