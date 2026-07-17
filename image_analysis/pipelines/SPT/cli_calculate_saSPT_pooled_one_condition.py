"""Compatibility entry point for pooled saSPT from a condition AIO table."""

from __future__ import annotations

import argparse

import pandas as pd

from spt_shared import atomic_csv
from step06_calculate_diffusion import run_saspt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("input")
    parser.add_argument("--frame-interval", type=float, required=True); parser.add_argument("--pixel-size", type=float, default=0.117)
    parser.add_argument("--focal-depth", type=float, default=0.7); parser.add_argument("--output", required=True)
    args = parser.parse_args(); table = pd.read_csv(args.input)
    result = run_saspt(table, frame_interval_s=args.frame_interval, pixel_size_um=args.pixel_size, focal_depth_um=args.focal_depth)
    atomic_csv(result, args.output); print(args.output)


if __name__ == "__main__":
    main()
