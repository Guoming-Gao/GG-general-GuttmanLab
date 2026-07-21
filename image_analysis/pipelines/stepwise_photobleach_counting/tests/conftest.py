from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import tifffile


def write_stack(path: Path, data: np.ndarray, *, pixel_size_um: float = 0.117, laser_power: float = 15.0) -> None:
    metadata = {
        "Frames": int(data.shape[0]), "FramesPerSecond": 10.0, "Exposure_ms": 100.0,
        "PixelSize_um": pixel_size_um, "LaserWavelength_nm": [405, 488, 561, 640],
        "LaserActive": [False, False, False, True], "LaserPowerPercent": [0, 0, 0, laser_power],
        "cameraBinning": 1,
    }
    with tifffile.TiffWriter(path) as writer:
        for index, frame in enumerate(data):
            writer.write(frame, metadata=None, description=json.dumps(metadata) if index == 0 else None)

