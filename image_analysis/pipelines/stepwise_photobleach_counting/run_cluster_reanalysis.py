"""Run the revised cluster-, condensate-, and monomer-aware PBSA analysis."""

from __future__ import annotations

import argparse

from cluster_counting import annotate_monomer_validation, calibrate_thresholds, count_all_streams, run_threshold_sensitivity
from cluster_reanalysis import build_analysis_units, extract_all_traces, load_reanalysis_config, verify_unit_invariants
from cluster_reporting import generate_report


STAGES = ("build-units", "extract-traces", "calibrate-thresholds", "count", "threshold-sensitivity", "validate-monomers", "report", "run")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("stage", choices=STAGES); parser.add_argument("--config", default="config_reanalysis.yaml"); parser.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv); cfg = load_reanalysis_config(args.config)
    stages = STAGES[:-1] if args.stage == "run" else (args.stage,)
    for stage in stages:
        print(f"\n[{stage}]", flush=True)
        if stage == "build-units":
            build_analysis_units(cfg, resume=args.resume); print(verify_unit_invariants(cfg), flush=True)
        elif stage == "extract-traces": extract_all_traces(cfg, resume=args.resume)
        elif stage == "calibrate-thresholds": print(calibrate_thresholds(cfg).to_string(index=False), flush=True)
        elif stage == "count": count_all_streams(cfg, resume=args.resume)
        elif stage == "threshold-sensitivity": print(run_threshold_sensitivity(cfg).to_string(index=False), flush=True)
        elif stage == "validate-monomers": annotate_monomer_validation(cfg)
        elif stage == "report": print(generate_report(cfg), flush=True)


if __name__ == "__main__": main()
