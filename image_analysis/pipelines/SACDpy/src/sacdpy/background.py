from __future__ import annotations

import numpy as np
import pywt


def background_estimation(
    stack: np.ndarray,
    threshold: float = 1.0,
    decomposition_level: int = 7,
    wavelet_name: str = "db6",
    iterations: int = 3,
) -> np.ndarray:
    """Wavelet background estimation ported from SACDm's `background_estimation.m`."""

    imgs = np.asarray(stack, dtype=np.float32)
    if imgs.ndim == 2:
        imgs = imgs[:, :, None]
    if imgs.ndim != 3:
        raise ValueError("background_estimation expects a 2D image or YXT stack.")
    if decomposition_level < 1:
        raise ValueError("decomposition_level must be positive.")
    if iterations < 1:
        raise ValueError("iterations must be positive.")

    original_shape = imgs.shape
    x, y, _ = imgs.shape
    if x < y:
        pad_x = max(x, y) - x
        pad_y = max(x, y) - y
        imgs = np.pad(imgs, ((0, pad_x), (0, pad_y), (0, 0)), mode="constant")

    background = np.zeros_like(imgs, dtype=np.float32)
    for frame in range(imgs.shape[2]):
        initial = imgs[:, :, frame].astype(np.float64, copy=False)
        residual = initial.copy()
        b_iter = np.zeros_like(initial)
        level = min(decomposition_level, pywt.dwtn_max_level(initial.shape, wavelet_name))
        if level < 1:
            background[:, :, frame] = initial
            continue

        for _ in range(iterations):
            b_iter = _approximation_only_reconstruction(residual, wavelet_name, level)
            if threshold > 0:
                eps = np.sqrt(np.abs(residual)) / 2.0
                mask = initial > (b_iter + eps)
                residual[mask] = b_iter[mask] + eps[mask]
                b_iter = _approximation_only_reconstruction(residual, wavelet_name, level)
        background[:, :, frame] = b_iter.astype(np.float32, copy=False)

    background = background[: original_shape[0], : original_shape[1], :]
    if stack.ndim == 2:
        return background[:, :, 0]
    return background


def _approximation_only_reconstruction(image: np.ndarray, wavelet_name: str, level: int) -> np.ndarray:
    coeffs = pywt.wavedec2(image, wavelet_name, mode="symmetric", level=level)
    approx = coeffs[0]
    zero_details = tuple(
        tuple(np.zeros_like(detail, dtype=np.float64) for detail in detail_level)
        for detail_level in coeffs[1:]
    )
    reconstructed = pywt.waverec2((approx, *zero_details), wavelet_name, mode="symmetric")
    return reconstructed[: image.shape[0], : image.shape[1]]
