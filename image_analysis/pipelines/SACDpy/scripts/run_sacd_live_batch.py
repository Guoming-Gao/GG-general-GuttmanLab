#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from sacdpy.continuous_batch import (
    build_batch_plan,
    load_config,
    preflight_summary,
    prepare_legacy_provenance,
    run_batch,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run resumable SACD time/z batch reconstruction")
    parser.add_argument("config", type=Path)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--prepare-provenance", action="store_true")
    parser.add_argument("--max-new-fovs", type=int)
    args = parser.parse_args()

    config = load_config(args.config)
    plan = build_batch_plan(config)
    print(json.dumps(preflight_summary(plan), indent=2), flush=True)
    if args.preflight:
        return 0
    if args.prepare_provenance:
        prepare_legacy_provenance(config, args.config)
        return 0
    run_batch(config, args.config, max_new_fovs=args.max_new_fovs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
