"""Command-line equivalent of the sibling mouse smiFISH notebook."""

import argparse
from pathlib import Path

from .mouse_smifish import (
    adaptive_blast_verify, generate_candidates, load_models, load_run_config,
    select_sets, write_outputs,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="Local JSON config (default: sibling .smifish-local.json)")
    parser.add_argument("--output-parent", type=Path, help="Override configured output directory")
    parser.add_argument("--intron-pool", type=int, default=600)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    config = load_run_config(args.config)
    models = load_models(config["gtf"])
    candidates = generate_candidates(models, fasta=config["fasta"], intron_pool=args.intron_pool)
    verified, hits, command = adaptive_blast_verify(
        candidates, blastn=config["blastn"], database=config["blast_db"], threads=args.threads)
    selected, summary = select_sets(verified)
    output = write_outputs(models, candidates, verified, hits, selected, summary, command,
                           output_parent=args.output_parent or config["output_parent"],
                           reference_paths=config)
    print(summary.to_string(index=False))
    print(f"Saved run: {output}")
    if not summary.meets_minimum.all():
        print("Some sets have fewer than 30 BLAST-verified probes; see set_summary.csv")


if __name__ == "__main__":
    main()
