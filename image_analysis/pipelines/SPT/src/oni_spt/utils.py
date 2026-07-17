from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
import time
from contextlib import contextmanager
from importlib import metadata
from pathlib import Path
from typing import Any, Iterator

import pandas as pd


def atomic_csv(df: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def atomic_json(value: Any, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    with tmp.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, default=str)
    os.replace(tmp, path)


def output_action(paths: list[Path], *, resume: bool, force: bool) -> str:
    existing = [path for path in paths if path.exists()]
    if not existing:
        return "run"
    if force:
        return "run"
    if resume and len(existing) == len(paths):
        return "skip"
    names = ", ".join(str(path) for path in existing[:3])
    raise FileExistsError(
        f"Output already exists ({names}). Use --resume to reuse a complete stage "
        "or --force to replace individual outputs."
    )


def config_fingerprint(cfg: dict[str, Any]) -> str:
    public = {k: v for k, v in cfg.items() if not k.startswith("_")}
    payload = json.dumps(public, sort_keys=True, default=str).encode()
    return hashlib.sha256(payload).hexdigest()[:16]


def dependency_versions(names: list[str]) -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for name in names:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def provenance(cfg: dict[str, Any]) -> dict[str, Any]:
    return {
        "created_unix": time.time(),
        "python": sys.version,
        "platform": platform.platform(),
        "config_fingerprint": config_fingerprint(cfg),
        "dependencies": dependency_versions(
            ["oni-spt", "numpy", "pandas", "scipy", "scikit-image", "tifffile", "spotiflow", "cellpose", "laptrack", "saspt"]
        ),
    }


@contextmanager
def stage_timer(status_path: str | Path, stage: str, details: dict[str, Any] | None = None) -> Iterator[None]:
    path = Path(status_path)
    if path.exists():
        try:
            status = json.loads(path.read_text())
        except json.JSONDecodeError:
            status = {}
    else:
        status = {}
    entry = {"status": "running", "started_unix": time.time()}
    if details:
        entry.update(details)
    status[stage] = entry
    atomic_json(status, path)
    try:
        yield
    except Exception as exc:
        entry.update(
            status="failed",
            finished_unix=time.time(),
            error=f"{type(exc).__name__}: {exc}",
        )
        status[stage] = entry
        atomic_json(status, path)
        raise
    else:
        entry.update(status="complete", finished_unix=time.time())
        status[stage] = entry
        atomic_json(status, path)
