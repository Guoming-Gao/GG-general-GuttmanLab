from __future__ import annotations

import argparse
from pathlib import Path

from .analysis import analyze_aio_files, analyze_saspt_file, concat_aio_files
from .pipeline import (
    detect_stage,
    initialize_run,
    preprocess_stage,
    run_all,
    segment_stage,
    track_stage,
)


def _stage_options(parser: argparse.ArgumentParser, include_frames: bool = False) -> None:
    parser.add_argument("--config", required=True, help="YAML acquisition configuration")
    parser.add_argument("--max-files", type=int, default=None)
    if include_frames:
        parser.add_argument("--max-frames", type=int, default=None)
    policy = parser.add_mutually_exclusive_group()
    policy.add_argument("--resume", action="store_true", help="Reuse complete existing stage outputs")
    policy.add_argument("--force", action="store_true", help="Replace existing individual outputs")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="oni-spt", description="Stand-alone ONI SPT pipeline")
    sub = parser.add_subparsers(dest="command", required=True)
    inspect = sub.add_parser("inspect", help="Validate inputs and write manifests")
    inspect.add_argument("--config", required=True)
    _stage_options(sub.add_parser("preprocess", help="Split/copy channels and apply DoG"), True)
    _stage_options(sub.add_parser("segment", help="Run CellposeSAM segmentation"))
    _stage_options(sub.add_parser("detect", help="Run Spotiflow per SPT frame"))
    _stage_options(sub.add_parser("track", help="Run LapTrack and cell assignment"))
    _stage_options(sub.add_parser("run", help="Run all stages"), True)

    aio = sub.add_parser("analyze-aio", help="Calculate classical SPT AIO tables")
    aio.add_argument("inputs", nargs="+")
    aio.add_argument("--frame-interval", type=float, required=True)
    aio.add_argument("--pixel-size", type=float, default=0.117)
    aio.add_argument("--output-dir", default=None)

    concat = sub.add_parser("concat-aio", help="Concatenate AIO files for one condition")
    concat.add_argument("inputs", nargs="+")
    concat.add_argument("--output", required=True)

    saspt = sub.add_parser("analyze-saspt", help="Run pooled saSPT from one AIO file")
    saspt.add_argument("input")
    saspt.add_argument("--frame-interval", type=float, required=True)
    saspt.add_argument("--pixel-size", type=float, default=0.117)
    saspt.add_argument("--focal-depth", type=float, default=0.7)
    saspt.add_argument("--output", default=None)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "inspect":
        cfg, manifest = initialize_run(args.config)
        accepted = int((manifest["status"] == "accepted").sum())
        print(f"Accepted: {accepted}; skipped/errors: {len(manifest) - accepted}")
        print(f"Manifest: {Path(cfg['output_dir']) / 'input_manifest.csv'}")
    elif args.command == "run":
        run_all(
            args.config,
            max_files=args.max_files,
            max_frames=args.max_frames,
            resume=args.resume,
            force=args.force,
        )
    elif args.command in {"preprocess", "segment", "detect", "track"}:
        cfg, manifest = initialize_run(args.config)
        function = {
            "preprocess": preprocess_stage,
            "segment": segment_stage,
            "detect": detect_stage,
            "track": track_stage,
        }[args.command]
        kwargs = {
            "max_files": args.max_files,
            "resume": args.resume,
            "force": args.force,
        }
        if args.command == "preprocess":
            kwargs["max_frames"] = args.max_frames
        function(cfg, manifest, **kwargs)
    elif args.command == "analyze-aio":
        outputs = analyze_aio_files(
            args.inputs,
            frame_interval_s=args.frame_interval,
            pixel_size_um=args.pixel_size,
            output_dir=args.output_dir,
        )
        print("\n".join(map(str, outputs)))
    elif args.command == "concat-aio":
        print(concat_aio_files(args.inputs, args.output))
    elif args.command == "analyze-saspt":
        print(
            analyze_saspt_file(
                args.input,
                frame_interval_s=args.frame_interval,
                pixel_size_um=args.pixel_size,
                focal_depth_um=args.focal_depth,
                output=args.output,
            )
        )


if __name__ == "__main__":
    main()
