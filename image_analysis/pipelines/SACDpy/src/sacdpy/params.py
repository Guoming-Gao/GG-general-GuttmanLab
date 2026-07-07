from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(slots=True)
class SACDParams:
    """Parameters for the core SACD reconstruction pipeline."""

    pixel_nm: float = 65.0
    wavelength_nm: float = 525.0
    na: float = 1.45
    mag: int = 2
    iter1: int = 7
    iter2: int = 8
    ac_order: int = 2
    frames_per_sacd: int | None = None
    scale: float | None = None
    subfactor: float = 0.8
    psf: np.ndarray | None = None
    resolution_nm: float | None = None
    ifregistration: bool = False
    ifbackground: bool = False
    backgroundfactor: float = 2.0
    ifsparsedecon: bool = False
    fidelity: float = 100.0
    tcontinuity: float = 0.1
    sparsity: float = 1.0
    sparse_iterations: int = 100

    def resolved_scale(self) -> float:
        return float(self.ac_order if self.scale is None else self.scale)

    def validate_core(self) -> None:
        if self.mag < 1:
            raise ValueError("mag must be >= 1.")
        if self.iter1 < 0 or self.iter2 < 0:
            raise ValueError("RL iteration counts must be nonnegative.")
        if not 2 <= self.ac_order <= 6:
            raise ValueError("ac_order must be in the MATLAB-supported range 2..6.")
        if self.frames_per_sacd is not None:
            if self.frames_per_sacd <= 0:
                raise ValueError("frames_per_sacd must be positive when provided.")
            if self.frames_per_sacd < self.ac_order:
                raise ValueError("frames_per_sacd must be >= ac_order.")
        if self.pixel_nm <= 0:
            raise ValueError("pixel_nm must be positive.")
        if self.psf is None:
            if self.resolution_nm is None and (self.wavelength_nm <= 0 or self.na <= 0):
                raise ValueError("wavelength_nm and na must be positive when psf/resolution_nm are not supplied.")
            if self.resolution_nm is not None and self.resolution_nm <= 0:
                raise ValueError("resolution_nm must be positive.")
        if self.backgroundfactor <= 0:
            raise ValueError("backgroundfactor must be positive.")
        if self.fidelity <= 0:
            raise ValueError("fidelity must be positive.")
        if self.tcontinuity < 0:
            raise ValueError("tcontinuity must be nonnegative.")
        if self.sparsity < 0:
            raise ValueError("sparsity must be nonnegative.")
        if self.sparse_iterations <= 0:
            raise ValueError("sparse_iterations must be positive.")
