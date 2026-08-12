"""Safely consolidate and atomically install the validated production result tree."""

from __future__ import annotations

import argparse
import importlib.metadata
import os
from pathlib import Path
import shutil
import yaml


EXPERIMENTS = (
    "20260624_ONI-gmgao-SPEN_SHA-stepwise_photobleach",
    "20260625_ONI-gmgao-SPEN_wOSS-stepwise_photobleach",
)


def load_paths(config_path: Path) -> tuple[Path, Path]:
    config = yaml.safe_load(config_path.read_text())
    staging = Path(config["output_root"]).expanduser().resolve()
    current = Path(config["source_results_root"]).expanduser().resolve()
    if staging == current or staging.parent != current.parent:
        raise RuntimeError("staging and current roots must be distinct siblings")
    return staging, current


def move_once(source: Path, destination: Path) -> None:
    if destination.exists():
        if source.exists():
            raise FileExistsError(f"both source and destination exist: {source} -> {destination}")
        return
    if not source.exists():
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(source), str(destination))


def prepare(config_path: Path, final_root: Path) -> None:
    staging, current = load_paths(config_path)
    if final_root.resolve() != current:
        raise RuntimeError("final root must exactly match source_results_root")
    for dataset in EXPERIMENTS:
        old_experiment = current / dataset
        new_experiment = staging / "04_experiments" / dataset
        move_once(old_experiment / "01_input_inspection", new_experiment / "input_QC")
        move_once(old_experiment / "02_drift_correction", new_experiment / "drift_correction")
        move_once(
            old_experiment / "03_roi_detection_and_growth",
            staging / "99_technical" / "exploratory_segmentation" / "legacy_growth_inputs" / dataset,
        )

    moved_inputs = 0
    for source in sorted((staging / "04_experiments").glob("*/quickpbsa_counts/*/*/*__quickpbsa_input.csv")):
        dataset = source.parents[3].name
        stream = source.parents[1].name
        fov = source.parent.name
        destination = staging / "99_technical" / "quickpbsa_native" / dataset / stream / fov / source.name
        if destination.exists():
            raise FileExistsError(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), str(destination))
        moved_inputs += 1
    if moved_inputs not in (0, 54):
        raise RuntimeError(f"expected 54 or 0 quickPBSA inputs, found {moved_inputs}")

    provenance = staging / "99_technical" / "provenance"
    provenance.mkdir(parents=True, exist_ok=True)
    workspace = config_path.parent.resolve()
    for pattern in ("*.py", "*.yaml", "*.yml", "README.md"):
        for source in workspace.glob(pattern):
            destination = provenance / source.name
            if source.resolve() == config_path.resolve():
                destination.write_text(source.read_text().replace(str(staging), str(final_root.resolve())))
            else:
                shutil.copy2(source, destination)
    tests_destination = provenance / "tests"
    tests_destination.mkdir(exist_ok=True)
    for source in (workspace / "tests").glob("*.py"):
        shutil.copy2(source, tests_destination / source.name)
    versions = sorted(f"{dist.metadata['Name']}=={dist.version}" for dist in importlib.metadata.distributions() if dist.metadata.get("Name"))
    (provenance / "conda_smlm_package_versions.txt").write_text("\n".join(versions) + "\n")

    for path in sorted(staging.rglob("*"), reverse=True):
        if path.name in {".DS_Store", ".cache", ".matplotlib", ".pytest_cache", "__pycache__"}:
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()
    for path in sorted((item for item in staging.rglob("*") if item.is_dir()), key=lambda item: len(item.parts), reverse=True):
        if not any(path.iterdir()):
            path.rmdir()

    if list(staging.glob("04_experiments/*/quickpbsa_counts/*/*/*__quickpbsa_input.csv")):
        raise RuntimeError("quickPBSA input duplicates remain in reader-facing experiment folders")
    if list(staging.rglob(".DS_Store")) or list(staging.rglob(".cache")):
        raise RuntimeError("cache artifacts remain")
    empty = [path for path in staging.rglob("*") if path.is_dir() and not any(path.iterdir())]
    if empty:
        raise RuntimeError(f"empty directories remain: {empty[:3]}")
    print({"prepared": str(staging), "moved_quickpbsa_inputs": moved_inputs})


def swap(config_path: Path, obsolete_root: Path) -> None:
    staging, current = load_paths(config_path)
    obsolete = obsolete_root.resolve()
    if obsolete.parent != current.parent or "obsolete_" not in obsolete.name:
        raise RuntimeError("obsolete root must be a clearly named sibling")
    if not current.exists() or not staging.exists() or obsolete.exists():
        raise RuntimeError("swap preconditions not met")
    os.rename(current, obsolete)
    try:
        os.rename(staging, current)
    except Exception:
        os.rename(obsolete, current)
        raise
    print({"installed": str(current), "obsolete": str(obsolete)})


def delete_obsolete(config_path: Path, obsolete_root: Path) -> None:
    config = yaml.safe_load(config_path.read_text())
    current = Path(config["source_results_root"]).expanduser().resolve()
    obsolete = obsolete_root.resolve()
    if obsolete.parent != current.parent or obsolete.name != "stepwise_photobleach_counting_obsolete_20260722":
        raise RuntimeError("refusing to delete an unexpected path")
    if not current.exists() or not (current / "01_report" / "pipeline_audit_and_reanalysis_report.pdf").exists():
        raise RuntimeError("validated replacement is not installed")
    if obsolete.exists():
        shutil.rmtree(obsolete)
    print({"permanently_deleted": str(obsolete)})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("prepare", "swap", "delete-obsolete"))
    parser.add_argument("--config", type=Path, default=Path("config_reanalysis.yaml"))
    parser.add_argument("--final-root", type=Path)
    parser.add_argument("--obsolete-root", type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    if args.action == "prepare":
        if args.final_root is None: parser.error("--final-root is required for prepare")
        prepare(config_path, args.final_root)
    elif args.action == "swap":
        if args.obsolete_root is None: parser.error("--obsolete-root is required for swap")
        swap(config_path, args.obsolete_root)
    else:
        if args.obsolete_root is None: parser.error("--obsolete-root is required for delete-obsolete")
        delete_obsolete(config_path, args.obsolete_root)


if __name__ == "__main__":
    main()
