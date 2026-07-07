from __future__ import annotations

import numpy as np


def richardson_lucy_image(image: np.ndarray, psf: np.ndarray, iterations: int) -> np.ndarray:
    if iterations == 0:
        return np.asarray(image, dtype=np.float64)
    kernel = np.asarray(psf, dtype=np.float64)
    kernel_sum = kernel.sum()
    if kernel_sum <= 0:
        raise ValueError("PSF must have positive sum.")
    result = _matlab_deconvlucy(np.asarray(image, dtype=np.float64), kernel / kernel_sum, int(iterations))
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


def _matlab_deconvlucy(image: np.ndarray, psf: np.ndarray, iterations: int) -> np.ndarray:
    """Default MATLAB `deconvlucy(I, PSF, NUMIT)` algorithm."""

    if image.ndim != 2 or psf.ndim != 2:
        raise ValueError("MATLAB deconvlucy port expects 2D image and PSF arrays.")
    if any(p > s for p, s in zip(psf.shape, image.shape)):
        raise ValueError("PSF must not exceed image size.")

    h = _psf2otf(psf, image.shape)
    j1 = np.asarray(image, dtype=np.float64)
    j2 = j1.copy()
    j3 = np.zeros_like(j2)
    previous = np.zeros((j2.size, 2), dtype=np.float64)
    weight = np.ones_like(j2)
    wi = np.maximum(weight * j1, 0.0)
    axes = tuple(range(j2.ndim))
    scale = np.real(np.fft.ifftn(np.conj(h) * np.fft.fftn(weight, axes=axes), axes=axes)) + np.sqrt(
        np.finfo(np.float64).eps
    )

    lambda_value = 2.0 if np.any(previous != 0) else 0.0
    for k in range(int(lambda_value) + 1, int(lambda_value) + iterations + 1):
        if k > 2:
            numerator = float(previous[:, 0].T @ previous[:, 1])
            denominator = float(previous[:, 1].T @ previous[:, 1] + np.finfo(np.float64).eps)
            lambda_value = max(min(numerator / denominator, 1.0), 0.0)
        y = np.maximum(j2 + lambda_value * (j2 - j3), 0.0)
        cc = _corelucy(y, h, wi)
        j3 = j2
        j2 = np.maximum(y * np.real(np.fft.ifftn(np.conj(h) * cc, axes=axes)) / scale, 0.0)
        previous = np.column_stack((j2.ravel() - y.ravel(), previous[:, 0]))
    return j2


def _corelucy(y: np.ndarray, h: np.ndarray, weighted_image: np.ndarray) -> np.ndarray:
    axes = tuple(range(y.ndim))
    reblurred = np.real(np.fft.ifftn(h * np.fft.fftn(y, axes=axes), axes=axes))
    reblurred[reblurred == 0] = np.finfo(np.float64).eps
    image_ratio = weighted_image / reblurred + np.finfo(np.float64).eps
    return np.fft.fftn(image_ratio, axes=axes)


def _psf2otf(psf: np.ndarray, out_shape: tuple[int, int]) -> np.ndarray:
    padded = np.zeros(out_shape, dtype=np.float64)
    padded[: psf.shape[0], : psf.shape[1]] = psf
    for axis, size in enumerate(psf.shape):
        padded = np.roll(padded, shift=-int(np.floor(size / 2)), axis=axis)
    return np.fft.fftn(padded)
