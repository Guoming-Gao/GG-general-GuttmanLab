"""Compatibility entry point for the historical AIO calculator."""

from __future__ import annotations

import argparse

from step06_calculate_diffusion import analyze_aio_files


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--frame-interval", type=float, required=True)
    parser.add_argument("--pixel-size", type=float, default=0.117)
    parser.add_argument("--output-dir")
    args = parser.parse_args()
    for path in analyze_aio_files(args.inputs, frame_interval_s=args.frame_interval,
                                  pixel_size_um=args.pixel_size, output_dir=args.output_dir):
        print(path)


if __name__ == "__main__":
    main()
