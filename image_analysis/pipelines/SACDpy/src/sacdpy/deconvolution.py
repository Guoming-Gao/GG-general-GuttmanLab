from __future__ import annotations

import numpy as np
from skimage.restoration import richardson_lucy


def richardson_lucy_image(image: np.ndarray, psf: np.ndarray, iterations: int) -> np.ndarray:
    if iterations == 0:
        return np.asarray(image, dtype=np.float64)
    kernel = np.asarray(psf, dtype=np.float64)
    kernel_sum = kernel.sum()
    if kernel_sum <= 0:
        raise ValueError("PSF must have positive sum.")
    kernel = kernel / kernel_sum
    result = richardson_lucy(
        np.asarray(image, dtype=np.float64),
        kernel,
        num_iter=int(iterations),
        clip=False,
        filter_epsilon=1e-12,
    )
    result = np.asarray(result, dtype=np.float64)
    result[~np.isfinite(result)] = 0.0
    result[result < 0] = 0.0
    return result


def richardson_lucy_stack(stack: np.ndarray, psf: np.ndarray, iterations: int) -> np.ndarray:
    arr = np.asarray(stack, dtype=np.float64)
    if arr.ndim == 2:
        return richardson_lucy_image(arr, psf, iterations)
    if arr.ndim != 3:
        raise ValueError("RL input must be 2D or YXT 3D.")
    out = np.empty(arr.shape, dtype=np.float64)
    for frame in range(arr.shape[2]):
        out[:, :, frame] = richardson_lucy_image(arr[:, :, frame], psf, iterations)
    return out
