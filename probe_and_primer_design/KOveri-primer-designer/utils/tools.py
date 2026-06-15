"""External tool helpers."""

from __future__ import annotations

import shutil
from pathlib import Path


def resolve_tool(config: dict, key: str, executable: str) -> str:
    """Return an executable path from config, PATH, or common bioinfo locations."""
    configured = config.get(key)
    if configured and Path(configured).exists():
        return str(configured)

    found = shutil.which(executable)
    if found:
        return found

    candidates = [
        Path("/opt/miniconda3/envs/bioinfo/bin") / executable,
        Path("/usr/local/bin") / executable,
        Path("/usr/bin") / executable,
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)

    raise FileNotFoundError(
        f"Could not find {executable}. Activate conda env 'bioinfo' or set config['{key}']."
    )


def require_file(path: str, label: str) -> str:
    """Validate that a required resource exists."""
    if not path or not Path(path).exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    return path
