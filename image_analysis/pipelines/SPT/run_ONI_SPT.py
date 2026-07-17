"""Notebook-friendly terminal runner for the complete ONI SPT workflow."""

from __future__ import annotations

import argparse

from step01_inspect_inputs import inspect_inputs
from step02_preprocess_videos import preprocess_stage
from step03_segment_cells import segment_stage
from step04_detect_spots import detect_stage
from step05_link_trajectories import track_stage
from step06_calculate_diffusion import diffusion_stage
from step07_generate_reports import report_stage


STAGES = ("inspect", "preprocess", "segment", "detect", "track", "diffusion", "report", "run")


def run_pipeline(config: str, *, max_files: int | None = None, max_frames: int | None = None,
                 resume: bool = False, force: bool = False, through_report: bool = True):
    cfg, manifest = inspect_inputs(config)
    preprocess_stage(cfg, manifest, max_files=max_files, max_frames=max_frames, resume=resume, force=force)
    segment_stage(cfg, manifest, max_files=max_files, resume=resume, force=force)
    detect_stage(cfg, manifest, max_files=max_files, resume=resume, force=force)
    track_stage(cfg, manifest, max_files=max_files, resume=resume, force=force)
    diffusion = diffusion_stage(cfg, manifest, max_files=max_files)
    report = report_stage(cfg) if through_report else {}
    return {"config": cfg, "manifest": manifest, "diffusion": diffusion, "report": report}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("stage", choices=STAGES)
    parser.add_argument("--config", required=True); parser.add_argument("--max-files", type=int); parser.add_argument("--max-frames", type=int)
    policy = parser.add_mutually_exclusive_group(); policy.add_argument("--resume", action="store_true"); policy.add_argument("--force", action="store_true")
    args = parser.parse_args(argv); cfg, manifest = inspect_inputs(args.config)
    if args.stage == "inspect": print(manifest); return
    if args.stage == "run":
        run_pipeline(args.config, max_files=args.max_files, max_frames=args.max_frames, resume=args.resume, force=args.force); return
    functions = {"preprocess": preprocess_stage, "segment": segment_stage, "detect": detect_stage,
                 "track": track_stage, "diffusion": diffusion_stage, "report": report_stage}
    function = functions[args.stage]
    if args.stage == "preprocess": function(cfg, manifest, max_files=args.max_files, max_frames=args.max_frames, resume=args.resume, force=args.force)
    elif args.stage in {"segment", "detect", "track"}: function(cfg, manifest, max_files=args.max_files, resume=args.resume, force=args.force)
    elif args.stage == "diffusion": function(cfg, manifest, max_files=args.max_files)
    else: function(cfg)


if __name__ == "__main__":
    main()
