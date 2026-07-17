from __future__ import annotations

import json
from pathlib import Path
import sys

# The revised workflow intentionally keeps its readable stage modules at the
# repository root instead of installing a package.  Make that root explicit
# for both ``pytest`` and ``python -m pytest`` launch styles.
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import numpy as np
import yaml
from tifffile import TiffWriter


def write_oni_tiff(
    path: Path,
    data: np.ndarray,
    *,
    expected_frames: int | None = None,
    pixel_size_um: float = 0.117,
    exposure_ms: float = 30.0,
) -> None:
    metadata = {
        "Frames": int(expected_frames if expected_frames is not None else len(data)),
        "FramesPerSecond": 1000.0 / exposure_ms,
        "Exposure_ms": exposure_ms,
        "PixelSize_um": pixel_size_um,
        "ROI": [0, 0, int(data.shape[2]), int(data.shape[1])],
    }
    with TiffWriter(path) as writer:
        for index, frame in enumerate(data):
            writer.write(
                frame,
                contiguous=True,
                description=json.dumps(metadata) if index == 0 else None,
                metadata=None,
            )


def write_config(
    path: Path,
    input_dir: Path,
    output_dir: Path,
    *,
    mode: str,
    width: int,
) -> Path:
    channels = (
        {"full": {"name": "spt", "role": "spt"}}
        if mode == "single"
        else {
            "left": {"name": "hoechst", "role": "marker"},
            "right": {"name": "spt", "role": "spt"},
        }
    )
    value = {
        "dataset": "test",
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "layout": {
            "mode": mode,
            "expected_width": width,
            "unexpected_width_reason": "unexpected_single_channel_width",
        },
        "channels": channels,
    }
    path.write_text(yaml.safe_dump(value))
    return path
