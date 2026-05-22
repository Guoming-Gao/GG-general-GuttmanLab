from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile


def read_tiff_stack(path: str | Path) -> np.ndarray:
    return tifffile.imread(path)


def write_tiff_image(path: str | Path, image: np.ndarray) -> None:
    """Write a 2D SACD image or a TYX SACD stack as float32 TIFF."""

    tifffile.imwrite(path, np.asarray(image, dtype=np.float32), imagej=True)
