from __future__ import annotations

import numpy as np

from .cumulant import cumulant
from .deconvolution import richardson_lucy_image, richardson_lucy_stack
from .fourier import fourier_interpolate
from .params import SACDParams
from .psf import make_psfs


def reconstruct(stack: np.ndarray, params: SACDParams | None = None) -> np.ndarray:
    """Run the core SACD reconstruction.

    Input stacks are accepted as `TYX` by default when the first axis is much
    smaller than the two spatial axes. Otherwise, 3D input is treated as `YXT`.
    When `params.frames_per_sacd` is set, full non-overlapping chunks are
    reconstructed and multiple results are returned as a `TYX` stack.
    """

    params = params or SACDParams()
    params.validate_core()

    yxt = as_yxt(stack)
    psfs = make_psfs(
        yxt.shape[:2],
        pixel_nm=params.pixel_nm,
        wavelength_nm=params.wavelength_nm,
        na=params.na,
        mag=params.mag,
        psf=params.psf,
        resolution_nm=params.resolution_nm,
    )

    if params.frames_per_sacd is None:
        return _reconstruct_single_window(yxt, params, psfs)

    window = params.frames_per_sacd
    assert window is not None
    n_windows = yxt.shape[2] // window
    if n_windows < 1:
        raise ValueError(
            f"Input has {yxt.shape[2]} frame(s), fewer than frames_per_sacd={window}; "
            "no full SACD chunks can be reconstructed."
        )

    results = [
        _reconstruct_single_window(yxt[:, :, start : start + window], params, psfs)
        for start in range(0, n_windows * window, window)
    ]
    if len(results) == 1:
        return results[0]
    return np.stack(results, axis=0).astype(np.float32, copy=False)


def _reconstruct_single_window(
    yxt: np.ndarray,
    params: SACDParams,
    psfs: tuple[np.ndarray, np.ndarray],
) -> np.ndarray:
    psf, psf_high = psfs
    work = yxt.astype(np.float64, copy=False)
    work = work - np.min(work, axis=(0, 1), keepdims=True)
    work[work < 0] = 0.0

    decon = richardson_lucy_stack(work, psf, params.iter1)
    interp = fourier_interpolate(decon, (params.mag, params.mag, 1), mirror_mode="lateral")
    interp[interp < 0] = 0.0

    subtraction = params.subfactor * np.mean(interp, axis=2)
    interp = np.abs(interp - subtraction[:, :, None])
    ac = np.abs(cumulant(interp, params.ac_order))

    post_kernel = psf_high ** params.resolved_scale()
    post_kernel = post_kernel / post_kernel.sum()
    result = richardson_lucy_image(ac, post_kernel, params.iter2)
    return result.astype(np.float32, copy=False)


def as_yxt(stack: np.ndarray) -> np.ndarray:
    arr = np.asarray(stack)
    if arr.ndim == 2:
        return arr[:, :, None]
    if arr.ndim != 3:
        raise ValueError("Input stack must be 2D or 3D.")
    if arr.shape[0] < arr.shape[1] and arr.shape[0] < arr.shape[2]:
        return np.moveaxis(arr, 0, 2)
    return arr
