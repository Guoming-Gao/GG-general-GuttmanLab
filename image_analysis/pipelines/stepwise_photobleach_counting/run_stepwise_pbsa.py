"""Run individual stages or the complete HILO stepwise-photobleaching workflow."""

from __future__ import annotations

import argparse
from typing import Any

from pbsa_shared import load_config


STAGES = (
    "inspect-inputs", "drift-correction", "detect-and-grow-rois",
    "extract-background-corrected-traces", "quickpbsa-threshold-pilot",
    "count-photobleaching-steps", "generate-reports", "run",
)


def run_stage(stage: str, cfg: dict[str, Any], *, resume: bool = False, force: bool = False,
              max_files: int | None = None, max_frames: int | None = None) -> None:
    from step01_inspect_inputs import inspect_dataset, load_manifest
    if stage == "quickpbsa-threshold-pilot":
        from step05_quickpbsa_threshold_pilot import pooled_threshold_pilot_stage
        pooled_threshold_pilot_stage(cfg, max_files=max_files)
        return
    if stage == "generate-reports":
        from step07_generate_reports import combined_reporting_stage, reporting_stage
        for dataset in cfg["datasets"]:
            reporting_stage(cfg, dataset)
        combined_reporting_stage(cfg)
        return
    for dataset in cfg["datasets"]:
        print(f"\n[{stage}] {dataset['name']}", flush=True)
        if stage == "inspect-inputs": inspect_dataset(cfg, dataset); continue
        manifest = load_manifest(dataset)
        if stage == "drift-correction":
            from step02_drift_correction import drift_correction_stage
            drift_correction_stage(cfg, dataset, manifest, resume=resume, force=force, max_files=max_files, max_frames=max_frames)
        elif stage == "detect-and-grow-rois":
            from step03_detect_and_grow_rois import roi_detection_stage
            roi_detection_stage(cfg, dataset, manifest, resume=resume, force=force, max_files=max_files)
        elif stage == "extract-background-corrected-traces":
            from step04_extract_background_corrected_traces import trace_extraction_stage
            trace_extraction_stage(cfg, dataset, manifest, resume=resume, force=force, max_files=max_files, max_frames=max_frames)
        elif stage == "count-photobleaching-steps":
            from step06_count_photobleaching_steps import step_counting_stage
            step_counting_stage(cfg, dataset, manifest, resume=resume, force=force, max_files=max_files)
        else:
            raise ValueError(stage)


def run_pipeline(cfg: dict[str, Any], **kwargs: Any) -> None:
    for stage in STAGES[:-2]:
        run_stage(stage, cfg, **kwargs)
    run_stage("generate-reports", cfg, **kwargs)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("stage", choices=STAGES); parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--max-files", type=int); parser.add_argument("--max-frames", type=int)
    policy = parser.add_mutually_exclusive_group(); policy.add_argument("--resume", action="store_true"); policy.add_argument("--force", action="store_true")
    args = parser.parse_args(argv); cfg = load_config(args.config)
    kwargs = {"resume": args.resume, "force": args.force, "max_files": args.max_files, "max_frames": args.max_frames}
    if args.stage == "run": run_pipeline(cfg, **kwargs)
    else: run_stage(args.stage, cfg, **kwargs)


if __name__ == "__main__": main()
