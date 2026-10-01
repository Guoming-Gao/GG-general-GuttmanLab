"""Command-line equivalent of the sibling mouse smiFISH notebook."""

import argparse
from pathlib import Path

from .config import DEFAULT_SETTINGS
from .mouse_smifish import (
    adaptive_blast_verify, generate_candidates, load_models, load_run_config,
    write_outputs,
)
from .minimum_span import select_minimum_span_sets


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="Local JSON config (default: sibling .smifish-local.json)")
    parser.add_argument("--output-parent", type=Path, help="Override configured output directory")
    parser.add_argument("--run-name", help="Name of a new output subfolder")
    parser.add_argument("--intron-pool", type=int, help="Exploratory intron candidate cap; omit for complete search")
    parser.add_argument("--blast-batch-size", type=int, default=3000)
    parser.add_argument("--max-homopolymer-length", type=int,
                        default=DEFAULT_SETTINGS["max_homopolymer_length"],
                        help="Maximum identical-base run in a probe sequence (default: 4)")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    config = load_run_config(args.config)
    models = load_models(config["gtf"])
    candidates = generate_candidates(models, fasta=config["fasta"], intron_pool=args.intron_pool,
                                     max_homopolymer_length=args.max_homopolymer_length)
    verified, hits, command = adaptive_blast_verify(
        candidates, blastn=config["blastn"], database=config["blast_db"],
        threads=args.threads, batch_size=args.blast_batch_size)
    selected, summary = select_minimum_span_sets(verified, models, probes_per_set=30)
    output = write_outputs(models, candidates, verified, hits, selected, summary, command,
                           output_parent=args.output_parent or config["output_parent"],
                           reference_paths=config, intron_pool=args.intron_pool,
                           run_name=args.run_name,
                           max_homopolymer_length=args.max_homopolymer_length)
    print(summary.to_string(index=False))
    print(f"Saved run: {output}")
    if not summary.meets_minimum.all():
        print("Some sets have fewer than 30 BLAST-verified probes; see set_summary.csv")


if __name__ == "__main__":
    main()
