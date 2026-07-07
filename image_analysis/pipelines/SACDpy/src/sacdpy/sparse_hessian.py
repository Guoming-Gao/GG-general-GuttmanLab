from __future__ import annotations

import numpy as np


def sparse_hessian_core(
    image: np.ndarray,
    fidelity: float = 100.0,
    continuity: float = 0.1,
    sparsity: float = 1.0,
    iterations: int = 100,
    mu: float = 1.0,
) -> np.ndarray:
    """CPU port of SACDm's `SparseHessian_core.m`."""

    if fidelity <= 0 or sparsity < 0 or continuity < 0 or iterations < 1 or mu <= 0:
        raise ValueError("Invalid sparse Hessian parameters.")

    f = np.asarray(image, dtype=np.float32)
    input_was_2d = f.ndim == 2
    if input_was_2d:
        f = f[:, :, None]
    if f.ndim != 3:
        raise ValueError("sparse_hessian_core expects a 2D image or 3D stack.")

    f_flag = f.shape[2]
    contiz = np.float32(np.sqrt(continuity))
    if f_flag < 3:
        contiz = np.float32(0.0)
        repeats = 3 - f.shape[2]
        f = np.concatenate([f, np.repeat(f[:, :, -1:], repeats, axis=2)], axis=2)

    max_value = float(np.max(f))
    if max_value <= 0:
        return np.zeros_like(image, dtype=np.float32)
    f = (f / max_value).astype(np.float32, copy=False)
    sizeg = f.shape

    operationfft = (
        _operation_xx(sizeg)
        + _operation_yy(sizeg)
        + (contiz**2) * _operation_zz(sizeg)
        + 2 * _operation_xy(sizeg)
        + 2 * contiz * _operation_xz(sizeg)
        + 2 * contiz * _operation_yz(sizeg)
    )
    normalize = (fidelity / mu) + sparsity**2 + operationfft

    bxx = np.zeros(sizeg, dtype=np.float32)
    byy = np.zeros(sizeg, dtype=np.float32)
    bzz = np.zeros(sizeg, dtype=np.float32)
    bxy = np.zeros(sizeg, dtype=np.float32)
    bxz = np.zeros(sizeg, dtype=np.float32)
    byz = np.zeros(sizeg, dtype=np.float32)
    bl1 = np.zeros(sizeg, dtype=np.float32)

    g_update = (fidelity / mu) * f
    for iteration in range(iterations):
        g_update_fft = np.fft.fftn(g_update)
        if iteration > 0:
            g = np.real(np.fft.ifftn(g_update_fft / normalize)).astype(np.float32, copy=False)
        else:
            g = np.real(np.fft.ifftn(g_update_fft / (fidelity / mu))).astype(np.float32, copy=False)
        g_update = (fidelity / mu) * f

        lxx, bxx = _iter_xx(g, bxx, 1.0, mu)
        g_update = g_update + lxx
        lyy, byy = _iter_yy(g, byy, 1.0, mu)
        g_update = g_update + lyy
        lzz, bzz = _iter_zz(g, bzz, float(contiz**2), mu)
        g_update = g_update + lzz
        lxy, bxy = _iter_xy(g, bxy, 2.0, mu)
        g_update = g_update + lxy
        lxz, bxz = _iter_xz(g, bxz, float(2 * contiz), mu)
        g_update = g_update + lxz
        lyz, byz = _iter_yz(g, byz, float(2 * contiz), mu)
        g_update = g_update + lyz
        lsparse, bl1 = _iter_sparse(g, bl1, sparsity, mu)
        g_update = g_update + lsparse

    g[g < 0] = 0
    if f_flag < 3:
        return g[:, :, 1].astype(np.float32, copy=False)
    return g.astype(np.float32, copy=False)


def _kernel_fft(kernel: np.ndarray, shape: tuple[int, int, int]) -> np.ndarray:
    full = np.zeros(shape, dtype=np.float32)
    slices = tuple(slice(0, size) for size in kernel.shape)
    full[slices] = kernel
    spec = np.fft.fftn(full)
    return spec * np.conj(spec)


def _operation_xx(shape: tuple[int, int, int]) -> np.ndarray:
    return _kernel_fft(np.array([[[1.0], [-2.0], [1.0]]], dtype=np.float32), shape)


def _operation_yy(shape: tuple[int, int, int]) -> np.ndarray:
    return _kernel_fft(np.array([[[1.0]], [[-2.0]], [[1.0]]], dtype=np.float32), shape)


def _operation_zz(shape: tuple[int, int, int]) -> np.ndarray:
    return _kernel_fft(np.array([[[1.0, -2.0, 1.0]]], dtype=np.float32), shape)


def _operation_xy(shape: tuple[int, int, int]) -> np.ndarray:
    return _kernel_fft(np.array([[[1.0], [-1.0]], [[-1.0], [1.0]]], dtype=np.float32), shape)


def _operation_xz(shape: tuple[int, int, int]) -> np.ndarray:
    return _kernel_fft(np.array([[[1.0, -1.0], [-1.0, 1.0]]], dtype=np.float32), shape)


def _operation_yz(shape: tuple[int, int, int]) -> np.ndarray:
    return _kernel_fft(np.array([[[1.0, -1.0]], [[-1.0, 1.0]]], dtype=np.float32), shape)


def _shrink(value: np.ndarray, mu: float) -> np.ndarray:
    out = np.abs(value) - 1.0 / mu
    out[out < 0] = 0
    return out * np.sign(value)


def _forward_diff(data: np.ndarray, step: float, dim: int) -> np.ndarray:
    shift = np.array(data.shape)
    position = np.zeros(3, dtype=int)
    shiftdata1 = np.zeros(tuple(shift + 1), dtype=data.dtype)
    shiftdata2 = np.zeros(tuple(shift + 1), dtype=data.dtype)
    shift[dim] += 1
    position[dim] += 1
    target = tuple(slice(position[d], shift[d]) for d in range(3))
    shiftdata1[target] = data
    shiftdata2[target] = data
    shift[dim] -= 1
    source = tuple(slice(0, shift[d]) for d in range(3))
    shiftdata2[source] = data
    diff = (shiftdata1 - shiftdata2) / step
    shift[dim] += 1
    return -diff[target]


def _back_diff(data: np.ndarray, step: float, dim: int) -> np.ndarray:
    shift = np.array(data.shape)
    position = np.zeros(3, dtype=int)
    shiftdata1 = np.zeros(tuple(shift + 1), dtype=data.dtype)
    shiftdata2 = np.zeros(tuple(shift + 1), dtype=data.dtype)
    source = tuple(slice(position[d], shift[d]) for d in range(3))
    shiftdata1[source] = data
    shiftdata2[source] = data
    shift[dim] += 1
    position[dim] += 1
    target = tuple(slice(position[d], shift[d]) for d in range(3))
    shiftdata2[target] = data
    diff = (shiftdata1 - shiftdata2) / step
    shift[dim] -= 1
    out = tuple(slice(0, shift[d]) for d in range(3))
    return diff[out]


def _iter_sparse(g: np.ndarray, b: np.ndarray, para: float, mu: float) -> tuple[np.ndarray, np.ndarray]:
    d = _shrink(g + b, mu)
    b = b + (g - d)
    return para * (d - b), b


def _iter_xx(g: np.ndarray, b: np.ndarray, para: float, mu: float) -> tuple[np.ndarray, np.ndarray]:
    gxx = _back_diff(_forward_diff(g, 1.0, 0), 1.0, 0)
    d = _shrink(gxx + b, mu)
    b = b + (gxx - d)
    return para * _back_diff(_forward_diff(d - b, 1.0, 0), 1.0, 0), b


def _iter_yy(g: np.ndarray, b: np.ndarray, para: float, mu: float) -> tuple[np.ndarray, np.ndarray]:
    gyy = _back_diff(_forward_diff(g, 1.0, 1), 1.0, 1)
    d = _shrink(gyy + b, mu)
    b = b + (gyy - d)
    return para * _back_diff(_forward_diff(d - b, 1.0, 1), 1.0, 1), b


def _iter_zz(g: np.ndarray, b: np.ndarray, para: float, mu: float) -> tuple[np.ndarray, np.ndarray]:
    gzz = _back_diff(_forward_diff(g, 1.0, 2), 1.0, 2)
    d = _shrink(gzz + b, mu)
    b = b + (gzz - d)
    return para * _back_diff(_forward_diff(d - b, 1.0, 2), 1.0, 2), b


def _iter_xy(g: np.ndarray, b: np.ndarray, para: float, mu: float) -> tuple[np.ndarray, np.ndarray]:
    gxy = _forward_diff(_forward_diff(g, 1.0, 0), 1.0, 1)
    d = _shrink(gxy + b, mu)
    b = b + (gxy - d)
    return para * _back_diff(_back_diff(d - b, 1.0, 1), 1.0, 0), b


def _iter_xz(g: np.ndarray, b: np.ndarray, para: float, mu: float) -> tuple[np.ndarray, np.ndarray]:
    gxz = _forward_diff(_forward_diff(g, 1.0, 0), 1.0, 2)
    d = _shrink(gxz + b, mu)
    b = b + (gxz - d)
    return para * _back_diff(_back_diff(d - b, 1.0, 2), 1.0, 0), b


def _iter_yz(g: np.ndarray, b: np.ndarray, para: float, mu: float) -> tuple[np.ndarray, np.ndarray]:
    gyz = _forward_diff(_forward_diff(g, 1.0, 1), 1.0, 2)
    d = _shrink(gyz + b, mu)
    b = b + (gyz - d)
    return para * _back_diff(_back_diff(d - b, 1.0, 2), 1.0, 1), b
