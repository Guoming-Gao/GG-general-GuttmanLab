"""Preflight, run, resume, and monitor a production batch of ONI SPT datasets."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from run_ONI_SPT import run_pipeline
from spt_shared import atomic_json, atomic_text, config_fingerprint, load_config, output_dirs
from step01_inspect_inputs import build_manifest


STAGE_COST = {
    "preprocess": (0.031, 0.0),
    "segment": (0.0, 4.5),
    "detect": (0.22, 2.0),
    "track": (0.025, 7.0),
    "diffusion": (0.0, 30.0),
    "report": (0.0, 10.0),
}


def load_batch(path: str | Path) -> dict[str, Any]:
    batch_path = Path(path).expanduser().resolve()
    with batch_path.open() as handle:
        batch = yaml.safe_load(handle) or {}
    if not batch.get("name") or not batch.get("datasets"):
        raise ValueError("Batch config requires name and datasets")
    configs = []
    for item in batch["datasets"]:
        value = item["config"] if isinstance(item, dict) else item
        config_path = Path(value).expanduser()
        if not config_path.is_absolute():
            config_path = batch_path.parent / config_path
        configs.append(config_path.resolve())
    batch["_path"] = batch_path
    batch["_configs"] = configs
    return batch


def _estimated_output_bytes(cfg: dict[str, Any], manifest: pd.DataFrame) -> int:
    accepted = manifest[manifest.status == "accepted"]
    input_bytes = sum(Path(path).stat().st_size for path in accepted.path)
    multiplier = 3 if cfg["layout"]["mode"] == "single" else 2
    qc_and_tables = len(accepted) * 60 * 2**20
    return int(input_bytes * multiplier + qc_and_tables)


def preflight_batch(path: str | Path, *, require_space: bool = True) -> tuple[dict[str, Any], list[dict[str, Any]], pd.DataFrame]:
    batch = load_batch(path)
    records, prepared = [], []
    for item, config_path in zip(batch["datasets"], batch["_configs"]):
        expected = item if isinstance(item, dict) else {}
        cfg = load_config(config_path)
        manifest = build_manifest(cfg)
        accepted = manifest[manifest.status == "accepted"]
        skipped = manifest[manifest.status != "accepted"]
        frames = int(accepted.actual_frames.sum())
        checks = {
            "accepted_files": len(accepted), "accepted_frames": frames, "skipped_files": len(skipped),
        }
        for key, actual in checks.items():
            expected_value = expected.get(f"expected_{key}")
            if expected_value is not None and int(expected_value) != int(actual):
                raise RuntimeError(f"{cfg['dataset']}: expected {key}={expected_value}, observed {actual}")
        output_bytes = _estimated_output_bytes(cfg, manifest)
        records.append({
            "dataset": cfg["dataset"], **checks,
            "input_GiB": sum(Path(value).stat().st_size for value in accepted.path) / 2**30,
            "estimated_output_GiB": output_bytes / 2**30,
            "output_dir": cfg["output_dir"], "config": str(config_path),
        })
        prepared.append({"cfg": cfg, "manifest": manifest, "config_path": config_path})
    table = pd.DataFrame.from_records(records)
    base_output = Path(batch["base_output_dir"]).expanduser().resolve()
    required_gib = max(float(batch.get("minimum_free_GiB", 200)), float(table.estimated_output_GiB.sum()) * 1.25)
    free_gib = shutil.disk_usage(base_output).free / 2**30
    if require_space and free_gib < required_gib:
        raise RuntimeError(f"Insufficient free space: {free_gib:.1f} GiB available, {required_gib:.1f} GiB required")
    table.attrs.update(free_GiB=free_gib, required_GiB=required_gib)
    return batch, prepared, table


def _work(stage: str, row: dict[str, Any]) -> float:
    per_frame, per_file = STAGE_COST[stage]
    return per_frame * int(row.get("actual_frames", 0) or 0) + per_file


class BatchProgress:
    def __init__(self, batch: dict[str, Any], prepared: list[dict[str, Any]], status_path: Path):
        self.batch = batch
        self.prepared = prepared
        self.status_path = status_path
        self.started = time.time()
        self.completed_work = 0.0
        self.completed_keys: set[tuple[str, str, str]] = set()
        self.current_dataset = ""
        self.total_work = 0.0
        for value in prepared:
            rows = value["manifest"][value["manifest"].status == "accepted"].to_dict("records")
            for row in rows:
                for stage in ("preprocess", "segment", "detect", "track"):
                    self.total_work += _work(stage, row)
            self.total_work += STAGE_COST["diffusion"][1] + STAGE_COST["report"][1]
        self.state: dict[str, Any] = {
            "batch": batch["name"], "status": "running", "started_unix": self.started,
            "nominal_seconds": self.total_work, "datasets": {},
        }

    def set_dataset(self, dataset: str, status: str, error: str | None = None) -> None:
        self.current_dataset = dataset
        entry = self.state["datasets"].setdefault(dataset, {})
        entry.update(status=status, updated_unix=time.time())
        if error: entry["error"] = error
        self._write()

    def callback(self, stage: str, event: str, index: int, total: int, row: dict[str, Any]) -> None:
        filename = str(row.get("filename", ""))
        key = (self.current_dataset, stage, filename)
        if event == "finish" and key not in self.completed_keys:
            self.completed_keys.add(key)
            self.completed_work += _work(stage, row)
        elapsed = time.time() - self.started
        if self.completed_work >= 20:
            scale = min(3.0, max(0.5, elapsed / self.completed_work))
        else:
            scale = 1.0
        remaining = max(0.0, self.total_work - self.completed_work) * scale
        eta = datetime.now() + timedelta(seconds=remaining)
        self.state.update({
            "current_dataset": self.current_dataset, "current_stage": stage,
            "current_file": filename, "file_index": index, "files_in_stage": total,
            "event": event, "elapsed_seconds": elapsed, "estimated_remaining_seconds": remaining,
            "estimated_completion": eta.isoformat(timespec="seconds"),
            "completed_fraction": self.completed_work / self.total_work if self.total_work else 1.0,
        })
        self._write()
        if event == "finish":
            print(
                f"[batch ETA] {self.current_dataset} | {stage} {index}/{total} | "
                f"elapsed {timedelta(seconds=int(elapsed))} | ETA {eta:%Y-%m-%d %H:%M}",
                flush=True,
            )

    def finish(self, failures: dict[str, str]) -> None:
        self.state.update(
            status="complete_with_errors" if failures else "complete",
            finished_unix=time.time(), failures=failures,
        )
        self._write()

    def _write(self) -> None:
        atomic_json(self.state, self.status_path)
        if self.current_dataset:
            match = next((value for value in self.prepared if value["cfg"]["dataset"] == self.current_dataset), None)
            if match:
                atomic_json(self.state, output_dirs(match["cfg"]["output_dir"])["metadata"] / "batch_progress.json")


def run_batch(path: str | Path, *, resume: bool = False, force: bool = False) -> dict[str, str]:
    batch, prepared, table = preflight_batch(path)
    print(table.to_string(index=False), flush=True)
    nominal = sum(_work(stage, row) for value in prepared
                  for row in value["manifest"][value["manifest"].status == "accepted"].to_dict("records")
                  for stage in ("preprocess", "segment", "detect", "track")) + len(prepared) * 40
    print(
        f"Initial runtime estimate: {timedelta(seconds=int(nominal * 0.85))}–"
        f"{timedelta(seconds=int(nominal * 1.45))}; free space {table.attrs['free_GiB']:.1f} GiB",
        flush=True,
    )
    status_dir = Path(batch["base_output_dir"]).expanduser().resolve() / "_batch_status"
    status_dir.mkdir(parents=True, exist_ok=True)
    status_path = status_dir / f"{batch['name']}.json"
    atomic_text(Path(batch["_path"]).read_text(), status_dir / f"{batch['name']}.yaml")
    progress = BatchProgress(batch, prepared, status_path)
    failures: dict[str, str] = {}
    command = " ".join(sys.argv)
    for value in prepared:
        cfg, config_path = value["cfg"], value["config_path"]
        dataset = cfg["dataset"]
        progress.set_dataset(dataset, "running")
        try:
            run_pipeline(
                str(config_path), resume=resume, force=force, progress=progress.callback,
                batch_config_path=str(batch["_path"]), command=command,
            )
        except Exception as exc:
            message = f"{type(exc).__name__}: {exc}"
            failures[dataset] = message
            progress.set_dataset(dataset, "failed", message)
            print(f"[batch error] {dataset}: {message}", flush=True)
            if not bool(batch.get("continue_on_dataset_error", True)):
                break
        else:
            progress.set_dataset(dataset, "complete")
    progress.finish(failures)
    if failures:
        raise RuntimeError(f"Batch completed with {len(failures)} failed dataset(s); see {status_path}")
    return failures


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("preflight", "run")); parser.add_argument("--config", required=True)
    policy = parser.add_mutually_exclusive_group(); policy.add_argument("--resume", action="store_true"); policy.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    if args.action == "preflight":
        _, _, table = preflight_batch(args.config)
        print(table.to_string(index=False))
        print(f"Total accepted files: {int(table.accepted_files.sum())}")
        print(f"Total accepted frames: {int(table.accepted_frames.sum())}")
        print(f"Estimated outputs: {table.estimated_output_GiB.sum():.1f} GiB")
        print(f"Free space: {table.attrs['free_GiB']:.1f} GiB")
        return
    run_batch(args.config, resume=args.resume, force=args.force)


if __name__ == "__main__":
    main()
