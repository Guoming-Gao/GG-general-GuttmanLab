from __future__ import annotations

import math

import numpy as np


def cumulant(stack: np.ndarray, order: int) -> np.ndarray:
    """Calculate SACDm temporal autocumulants for orders 2 through 6."""

    im = np.asarray(stack)
    if im.ndim != 3:
        raise ValueError("cumulant expects a YXT stack.")
    if not 2 <= order <= 6:
        raise ValueError("order must be in the MATLAB-supported range 2..6.")
    if im.shape[2] < order:
        raise ValueError(f"Need at least {order} frames for cumulant order {order}.")

    windows = [im[:, :, i : im.shape[2] - order + i + 1] for i in range(order)]

    if order <= 3:
        return _moment(windows, range(order))

    # Joint cumulant from the partition formula. This matches the explicit
    # MATLAB expressions for orders 4-6 while keeping the implementation small.
    result = np.zeros(im.shape[:2], dtype=np.result_type(im.dtype, np.float64))
    for partition in _set_partitions(tuple(range(order))):
        if any(len(block) == 1 for block in partition):
            continue
        k = len(partition)
        coeff = ((-1) ** (k - 1)) * math.factorial(k - 1)
        term = None
        for block in partition:
            block_moment = _moment(windows, block)
            term = block_moment if term is None else term * block_moment
        result = result + coeff * term
    return result


def _moment(windows: list[np.ndarray], indices: tuple[int, ...] | range) -> np.ndarray:
    prod = None
    for idx in indices:
        prod = windows[idx] if prod is None else prod * windows[idx]
    return np.mean(prod, axis=2)


def _set_partitions(items: tuple[int, ...]) -> list[tuple[tuple[int, ...], ...]]:
    if len(items) == 1:
        return [((items[0],),)]
    first, rest = items[0], items[1:]
    partitions = []
    for partition in _set_partitions(rest):
        partitions.append(((first,),) + partition)
        for i in range(len(partition)):
            merged = list(partition)
            merged[i] = tuple(sorted((first,) + merged[i]))
            partitions.append(tuple(merged))
    normalized = []
    seen = set()
    for partition in partitions:
        key = tuple(sorted(tuple(block) for block in partition))
        if key not in seen:
            seen.add(key)
            normalized.append(key)
    return normalized
