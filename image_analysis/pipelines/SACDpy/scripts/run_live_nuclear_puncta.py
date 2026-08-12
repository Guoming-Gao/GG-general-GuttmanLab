#!/usr/bin/env python3
"""Run the frame-2 SHA/FVP/recovery nucleus and puncta analysis."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from sacdpy.live_nuclear_puncta import (  # noqa: E402
    config_from_json,
    discovery_summary,
    finalize_analysis,
    prepare_analysis,
    refresh_published_qc,
    replot_published_cell_statistics,
    run_fixed_threshold_pilot,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    parser.add_argument(
        "--stage",
        required=True,
        choices=(
            "discover",
            "prepare",
            "fixed-threshold-pilot",
            "finalize",
            "refresh-qc",
            "replot-cell-statistics",
        ),
    )
    parser.add_argument(
        "--approve-pilot",
        action="store_true",
        help="Required for finalize after the blinded default-parameter pilot is inspected.",
    )
    parser.add_argument(
        "--resume-staging",
        action="store_true",
        help="Resume prepare from a matching complete set of staged Cellpose label maps.",
    )
    parser.add_argument(
        "--reuse-segmentation-from",
        type=Path,
        help="Published analysis root whose validated Cellpose label maps should be reused.",
    )
    parser.add_argument(
        "--replace-existing",
        action="store_true",
        help="Permit transactional replacement of the existing final output after validation.",
    )
    parser.add_argument(
        "--detector-mode",
        choices=("default", "fixed-threshold"),
        default="default",
        help="Full-run detector; fixed-threshold remains a separately labeled sensitivity mode.",
    )
    parser.add_argument(
        "--fixed-threshold-multiplier",
        type=float,
        help="One multiplier frozen from the saved fixed-threshold pilot.",
    )
    parser.add_argument(
        "--presentation-only",
        action="store_true",
        help=(
            "With --stage replot-cell-statistics, overwrite only the presentation "
            "puncta-count and corrected-brightness boxplots."
        ),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    config = config_from_json(args.config)
    if args.stage == "discover":
        result = discovery_summary(config)
    elif args.stage == "prepare":
        result = prepare_analysis(
            config,
            resume_staging=args.resume_staging,
            reuse_segmentation_from=args.reuse_segmentation_from,
            allow_existing_output=args.replace_existing,
        )
    elif args.stage == "fixed-threshold-pilot":
        result = run_fixed_threshold_pilot(config)
    elif args.stage == "refresh-qc":
        result = refresh_published_qc(
            config,
            config_path=args.config.resolve(),
            script_path=Path(__file__).resolve(),
        )
    elif args.stage == "replot-cell-statistics":
        result = replot_published_cell_statistics(
            config,
            config_path=args.config.resolve(),
            script_path=Path(__file__).resolve(),
            presentation_only=args.presentation_only,
        )
    else:
        result = finalize_analysis(
            config,
            approve_pilot=args.approve_pilot,
            detector_mode=args.detector_mode,
            fixed_threshold_multiplier=args.fixed_threshold_multiplier,
            config_path=args.config.resolve(),
            script_path=Path(__file__).resolve(),
            replace_existing=args.replace_existing,
        )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
