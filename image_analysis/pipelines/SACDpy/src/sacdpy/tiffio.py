from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile


def read_tiff_stack(path: str | Path) -> np.ndarray:
    return tifffile.imread(path)


def write_tiff_image(
    path: str | Path,
    image: np.ndarray,
    *,
    intensity_transform: str = "order_root",
) -> None:
    """Write a 2D SACD image or a TYX SACD stack as float32 TIFF."""

    arr = np.asarray(image, dtype=np.float32)
    axes = "YX" if arr.ndim == 2 else "TYX"
    tifffile.imwrite(
        path,
        arr,
        imagej=True,
        metadata={"axes": axes, "intensity_transform": intensity_transform},
    )


def read_intensity_transform(path: str | Path) -> str | None:
    """Read SACD intensity-space provenance from an ImageJ TIFF."""

    with tifffile.TiffFile(path) as tif:
        metadata = tif.imagej_metadata or {}
    value = metadata.get("intensity_transform")
    return str(value) if value is not None else None
