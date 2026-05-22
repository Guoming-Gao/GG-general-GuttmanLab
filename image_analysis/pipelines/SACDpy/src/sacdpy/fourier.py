from __future__ import annotations

import numpy as np


def fourier_interpolate(
    image: np.ndarray,
    factors: int | tuple[int, ...],
    mirror_mode: str = "none",
) -> np.ndarray:
    """Fourier-domain interpolation ported from SACDm's `fourierInterpolation.m`."""

    arr = np.asarray(image)
    if isinstance(factors, int):
        fac = np.full(arr.ndim, factors, dtype=int)
    else:
        fac = np.asarray(factors, dtype=int)
    if fac.size == 1:
        fac = np.full(arr.ndim, int(fac[0]), dtype=int)
    if fac.size != arr.ndim:
        raise ValueError("Need one interpolation factor or one per image dimension.")
    if np.any(fac < 1):
        raise ValueError("Interpolation factors must be >= 1.")
    if np.all(fac == 1):
        return arr.copy()

    input_shape = np.array(arr.shape, dtype=int)
    noip = fac == 1

    if mirror_mode == "none":
        return _finterp(arr, np.rint(fac * input_shape).astype(int))
    if mirror_mode not in {"lateral", "axial", "both"}:
        raise ValueError(f"Unknown mirror_mode {mirror_mode!r}.")
    if arr.ndim not in {2, 3}:
        raise ValueError("Only 2D and 3D arrays are supported.")
    if arr.ndim == 2 and mirror_mode != "lateral":
        raise ValueError("2D data only supports mirror_mode='lateral'.")

    pad_size = np.zeros(arr.ndim, dtype=float)
    if mirror_mode in {"lateral", "both"}:
        pad_size[:2] = np.array(arr.shape[:2], dtype=float) / 2.0
    if arr.ndim == 3 and mirror_mode in {"axial", "both"}:
        pad_size[2] = arr.shape[2] / 2.0
    pad_size[noip] = 0

    pad_width = [(int(np.ceil(p)), int(np.floor(p))) for p in pad_size]
    padded = np.pad(arr, pad_width, mode="symmetric")
    padded_shape = np.array(padded.shape, dtype=int)

    if mirror_mode == "lateral" and arr.ndim == 3:
        new_shape = np.array(
            [
                fac[0] * padded_shape[0] - (fac[0] - 1),
                fac[1] * padded_shape[1] - (fac[1] - 1),
                fac[2] * padded_shape[2],
            ],
            dtype=int,
        )
    elif mirror_mode == "axial" and arr.ndim == 3:
        new_shape = np.array(
            [
                fac[0] * padded_shape[0],
                fac[1] * padded_shape[1],
                fac[2] * padded_shape[2] - (fac[2] - 1),
            ],
            dtype=int,
        )
    else:
        new_shape = np.rint(fac * padded_shape - (fac - 1)).astype(int)

    interp = _finterp(padded, new_shape)
    return _valid_part(interp, input_shape, fac, noip, mirror_mode)


def _valid_part(
    arr: np.ndarray,
    input_shape: np.ndarray,
    fac: np.ndarray,
    noip: np.ndarray,
    mirror_mode: str,
) -> np.ndarray:
    sz = input_shape.copy()
    sz = sz - (sz % 2 == 0)
    idx_1based = np.ceil(sz / 2).astype(int) + 1 + (fac - 1) * np.floor(sz / 2).astype(int)
    start = idx_1based - 1
    stop = start + fac * input_shape

    slices: list[slice] = []
    for dim in range(arr.ndim):
        applies = not noip[dim]
        if arr.ndim == 3 and mirror_mode == "lateral" and dim == 2:
            applies = False
        if applies:
            slices.append(slice(int(start[dim]), int(stop[dim])))
        else:
            slices.append(slice(None))
    return arr[tuple(slices)]


def _finterp(image: np.ndarray, new_shape: np.ndarray) -> np.ndarray:
    if image.ndim == 2:
        return _finterp_2d(image, new_shape)
    if image.ndim == 3:
        return _finterp_3d(image, new_shape)
    raise ValueError("Only 2D and 3D arrays are supported.")


def _resize_for_downsampling(old_shape: np.ndarray, new_shape: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    target = np.array(new_shape, dtype=int)
    incr = np.ones_like(target)
    for dim in range(target.size):
        if target[dim] < old_shape[dim]:
            incr[dim] = int(np.floor(old_shape[dim] / target[dim]) + 1)
            target[dim] = incr[dim] * target[dim]
    return target, incr


def _finterp_2d(image: np.ndarray, new_shape: np.ndarray) -> np.ndarray:
    old = np.array(image.shape, dtype=int)
    target, incr = _resize_for_downsampling(old, new_shape)
    if np.any(target == 0):
        return np.empty(tuple(target), dtype=np.float64)

    spec = np.fft.fft2(image) * (target[0] / old[0]) * (target[1] / old[1])
    out = np.zeros(tuple(target), dtype=np.complex128)
    nyq = np.ceil((old + 1) / 2).astype(int)

    out[: nyq[0], : nyq[1]] = spec[: nyq[0], : nyq[1]]
    out[-(old[0] - nyq[0]) :, : nyq[1]] = spec[nyq[0] :, : nyq[1]]
    out[: nyq[0], -(old[1] - nyq[1]) :] = spec[: nyq[0], nyq[1] :]
    out[-(old[0] - nyq[0]) :, -(old[1] - nyq[1]) :] = spec[nyq[0] :, nyq[1] :]

    if old[0] % 2 == 0 and target[0] != old[0]:
        out[nyq[0] - 1, :] /= 2.0
        out[nyq[0] + target[0] - old[0] - 1, :] = out[nyq[0] - 1, :]
    if old[1] % 2 == 0 and target[1] != old[1]:
        out[:, nyq[1] - 1] /= 2.0
        out[:, nyq[1] + target[1] - old[1] - 1] = out[:, nyq[1] - 1]

    img_ip = np.fft.ifft2(out).real
    return img_ip[:: incr[0], :: incr[1]]


def _finterp_3d(image: np.ndarray, new_shape: np.ndarray) -> np.ndarray:
    old = np.array(image.shape, dtype=int)
    target, incr = _resize_for_downsampling(old, new_shape)
    if np.any(target == 0):
        return np.empty(tuple(target), dtype=np.float64)

    spec = np.fft.fftn(image) * (target[0] / old[0]) * (target[1] / old[1]) * (target[2] / old[2])
    out = np.zeros(tuple(target), dtype=np.complex128)
    nyq = np.ceil((old + 1) / 2).astype(int)

    xlo, ylo, zlo = (slice(0, int(n)) for n in nyq)
    xhi = slice(-(old[0] - nyq[0]), None)
    yhi = slice(-(old[1] - nyq[1]), None)
    zhi = slice(-(old[2] - nyq[2]), None)
    sxhi = slice(int(nyq[0]), old[0])
    syhi = slice(int(nyq[1]), old[1])
    szhi = slice(int(nyq[2]), old[2])

    out[xlo, ylo, zlo] = spec[xlo, ylo, zlo]
    out[xhi, ylo, zlo] = spec[sxhi, ylo, zlo]
    out[xlo, yhi, zlo] = spec[xlo, syhi, zlo]
    out[xlo, ylo, zhi] = spec[xlo, ylo, szhi]
    out[xhi, yhi, zlo] = spec[sxhi, syhi, zlo]
    out[xhi, ylo, zhi] = spec[sxhi, ylo, szhi]
    out[xlo, yhi, zhi] = spec[xlo, syhi, szhi]
    out[xhi, yhi, zhi] = spec[sxhi, syhi, szhi]

    for dim in range(3):
        if old[dim] % 2 == 0 and target[dim] != old[dim]:
            src = [slice(None)] * 3
            dst = [slice(None)] * 3
            src[dim] = nyq[dim] - 1
            dst[dim] = nyq[dim] + target[dim] - old[dim] - 1
            out[tuple(src)] /= 2.0
            out[tuple(dst)] = out[tuple(src)]

    img_ip = np.fft.ifftn(out).real
    return img_ip[:: incr[0], :: incr[1], :: incr[2]]
