from __future__ import annotations

import numpy as np


def register_image(
    moving: np.ndarray,
    reference: np.ndarray,
    *,
    upsample_factor: int = 100,
    normalize: bool = True,
) -> np.ndarray:
    """Register `moving` to `reference` with SACDm-style Fourier subpixel shifting."""

    moving_arr = np.asarray(moving, dtype=np.float64)
    reference_arr = np.asarray(reference, dtype=np.float64)
    if moving_arr.shape != reference_arr.shape or moving_arr.ndim != 2:
        raise ValueError("register_image expects two same-shaped 2D images.")

    _, shift = drift_detect(reference_arr, moving_arr)
    shifted = np.abs(_shift_image(moving_arr, (-shift[0], -shift[1])))
    if normalize:
        max_value = float(np.max(shifted))
        if max_value > 0:
            shifted = shifted / max_value
    return shifted


def register_stack(stack: np.ndarray, *, upsample_factor: int = 100) -> np.ndarray:
    arr = np.asarray(stack, dtype=np.float64)
    if arr.ndim != 3:
        raise ValueError("register_stack expects a YXT stack.")
    out = arr.copy()
    reference = out[:, :, 0]
    for frame in range(1, out.shape[2]):
        out[:, :, frame] = register_image(out[:, :, frame], reference, upsample_factor=upsample_factor)
    return out


def drift_detect(
    reference: np.ndarray,
    moving: np.ndarray,
    *,
    max_iterations: int = 200,
    roi_radius: int = 5,
) -> tuple[np.ndarray, np.ndarray]:
    """Port of SACDm's `DriftDetect.m`."""

    im1 = np.asarray(reference, dtype=np.float64)
    im2 = np.asarray(moving, dtype=np.float64)
    if im1.shape != im2.shape or im1.ndim != 2:
        raise ValueError("drift_detect expects two same-shaped 2D images.")

    h, w = im1.shape
    nr = h * 2 - 1
    nc = w * 2 - 1
    corr = np.fft.fftshift(np.fft.ifft2(np.fft.fft2(im1, (nr, nc)) * np.conj(np.fft.fft2(im2, (nr, nc)))))
    ij, ji = np.unravel_index(int(np.argmax(np.abs(corr))), corr.shape)
    shift_integer = np.array([h - (ji + 1), w - (ij + 1)], dtype=np.float64)

    corr = np.real(corr)
    shifted_corr = _shift_image(corr, shift_integer)
    corr_h, corr_w = corr.shape
    row_start = int(np.floor(corr_w / 2) - roi_radius - 1)
    row_stop = int(np.floor(corr_w / 2) + roi_radius)
    col_start = int(np.floor(corr_h / 2) - roi_radius - 1)
    col_stop = int(np.floor(corr_h / 2) + roi_radius)
    sub_corr = shifted_corr[row_start:row_stop, col_start:col_stop]
    params = _g2d_fit(sub_corr, max_iterations=max_iterations)
    shift = np.array(
        [
            roi_radius - params["ux"] + shift_integer[0] + 1,
            roi_radius - params["uy"] + shift_integer[1] + 1,
        ],
        dtype=np.float64,
    )
    return shift_integer, shift


def _shift_image(image: np.ndarray, shift: tuple[float, float] | np.ndarray) -> np.ndarray:
    im = np.asarray(image, dtype=np.float64)
    h, w = im.shape
    xx, yy = _frequency_grid(im.shape)
    fft_im = np.fft.fftshift(np.fft.fft2(im))
    phase = np.exp(-1j * 2.0 * np.pi * ((xx * float(shift[0])) / h + (yy * float(shift[1])) / w))
    return np.real(np.fft.ifft2(np.fft.ifftshift(fft_im * phase)))


def _frequency_grid(shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    h, w = shape
    x = np.arange(-np.floor(w / 2), -np.floor(w / 2) + w, dtype=np.float64)
    y = np.arange(-np.floor(h / 2), -np.floor(h / 2) + h, dtype=np.float64)
    return np.meshgrid(x, y)


def _g2d_fit(data: np.ndarray, max_iterations: int, pixel_size: float = 1.0) -> dict[str, float]:
    if max_iterations < 1:
        max_iterations = 1
    if pixel_size <= 0:
        pixel_size = 1.0

    coincidence_level = 0.01
    minimum_cond_number = 1e-20
    arr = np.asarray(data, dtype=np.float64)
    height, width = arr.shape
    background_level = max(float(np.min(arr)), 0.0)
    coeff = np.zeros(6, dtype=np.float64)
    old_coeff = np.zeros(6, dtype=np.float64)
    model_data = np.zeros_like(arr, dtype=np.float64)
    current_iteration = 1

    for current_iteration in range(1, max_iterations + 1):
        amat = np.zeros((6, 6), dtype=np.float64)
        bvec = np.zeros(6, dtype=np.float64)
        for y_idx in range(height):
            for x_idx in range(width):
                value = max(float(arr[y_idx, x_idx]) - background_level, 0.0)
                x = x_idx * pixel_size
                y = y_idx * pixel_size
                if current_iteration > 1:
                    old_value = np.exp(np.dot(coeff, [1, x, y, x**2, y**2, x * y]))
                else:
                    old_value = value
                if value > 0:
                    da = np.array(
                        [
                            [1, x, y, x**2, y**2, x * y],
                            [x, x**2, x * y, x**3, x * y**2, x**2 * y],
                            [y, x * y, y**2, x**2 * y, y**3, x * y**2],
                            [x**2, x**3, x**2 * y, x**4, x**2 * y**2, x**3 * y],
                            [y**2, x * y**2, y**3, x**2 * y**2, y**4, x * y**3],
                            [x * y, x**2 * y, x * y**2, x**3 * y, x * y**3, x**2 * y**2],
                        ],
                        dtype=np.float64,
                    )
                    db = np.array(
                        [
                            np.log(value),
                            x * np.log(value),
                            y * np.log(value),
                            x**2 * np.log(value),
                            y**2 * np.log(value),
                            x * y * np.log(value),
                        ],
                        dtype=np.float64,
                    )
                    amat += da * old_value**2
                    bvec += db * old_value**2

        cond_num = 1.0 / np.linalg.cond(amat)
        if not np.isfinite(cond_num) or cond_num < minimum_cond_number:
            return _invalid_gaussian_parameters()

        coeff = np.linalg.solve(amat, bvec)
        tail_distance_sqr = -4.5 / min(coeff[3], coeff[4])
        hyp_x = -0.5 * coeff[1] / coeff[3]
        hyp_y = -0.5 * coeff[2] / coeff[4]
        mean_error = 0.0
        divisor = 0.0
        for y_idx in range(height):
            for x_idx in range(width):
                x = x_idx * pixel_size
                y = y_idx * pixel_size
                model_data[y_idx, x_idx] = np.exp(np.dot(coeff, [1, x, y, x**2, y**2, x * y]))
                if (x - hyp_x) ** 2 + (y - hyp_y) ** 2 >= tail_distance_sqr:
                    mean_error += arr[y_idx, x_idx] - model_data[y_idx, x_idx]
                    divisor += 1
        if divisor > 0:
            mean_error /= divisor
        background_level = max(mean_error, 0.0)
        if float(np.max(np.abs(coeff - old_coeff))) <= coincidence_level:
            break
        old_coeff = coeff.copy()

    if coeff[1] <= 0 or coeff[2] <= 0 or coeff[3] >= 0 or coeff[4] >= 0:
        return _invalid_gaussian_parameters()

    theta = 0.5 * np.arctan(coeff[5] / (coeff[4] - coeff[3]))
    sin2theta = np.sin(theta) ** 2
    cos2theta = np.cos(theta) ** 2
    s2x = 0.5 * np.cos(2 * theta) / (coeff[4] * sin2theta - coeff[3] * cos2theta)
    s2y = 0.5 * np.cos(2 * theta) / (coeff[3] * sin2theta - coeff[4] * cos2theta)
    ux = coeff[1] * (s2y * sin2theta + s2x * cos2theta) + 0.5 * coeff[2] * (s2y - s2x) * np.sin(2 * theta)
    uy = coeff[2] * (s2y * cos2theta + s2x * sin2theta) + 0.5 * coeff[1] * (s2y - s2x) * np.sin(2 * theta)
    sx = np.sqrt(s2x)
    sy = np.sqrt(s2y)
    peak = np.exp(np.dot(coeff, [1, ux, uy, ux**2, uy**2, ux * uy]))
    data_pos = np.maximum(arr - background_level, 0)
    rss = float(np.sum((data_pos - model_data) ** 2))
    tss = float(np.sum((data_pos - np.mean(data_pos)) ** 2))
    return {
        "sx": float(sx),
        "sy": float(sy),
        "ux": float(ux),
        "uy": float(uy),
        "peak": float(peak),
        "rSquare": float(1 - rss / tss) if tss else -1.0,
        "background": float(background_level),
        "iterations": float(current_iteration),
        "theta": float(theta),
    }


def _invalid_gaussian_parameters() -> dict[str, float]:
    return {
        "sx": -1.0,
        "sy": -1.0,
        "ux": -1.0,
        "uy": -1.0,
        "peak": -1.0,
        "rSquare": -1.0,
        "background": -1.0,
        "iterations": -1.0,
        "theta": -1.0,
    }
