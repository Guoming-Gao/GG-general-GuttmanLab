"""Compatibility entry point for concatenating historical AIO files."""

from __future__ import annotations

import argparse

from step06_calculate_diffusion import concat_aio_files


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+"); parser.add_argument("--output", required=True)
    args = parser.parse_args(); print(concat_aio_files(args.inputs, args.output))


if __name__ == "__main__":
    main()
