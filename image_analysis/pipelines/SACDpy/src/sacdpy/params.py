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
    ifsparsedecon: bool = False

    def resolved_scale(self) -> float:
        return float(self.ac_order if self.scale is None else self.scale)

    def validate_core(self) -> None:
        if self.ifregistration:
            raise NotImplementedError("Registration is not part of the core SACDpy port yet.")
        if self.ifbackground:
            raise NotImplementedError("Wavelet background subtraction is not part of the core SACDpy port yet.")
        if self.ifsparsedecon:
            raise NotImplementedError("Sparse Hessian deconvolution is not part of the core SACDpy port yet.")
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
