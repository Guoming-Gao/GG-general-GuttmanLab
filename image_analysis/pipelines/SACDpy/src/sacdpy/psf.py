from __future__ import annotations

import numpy as np
from scipy import special


def generate_rsf(gamma: float, n: int | None = None) -> np.ndarray:
    """Generate the Gaussian RSF used by SACDm's `generate_rsf.m`."""

    if gamma <= 0:
        raise ValueError("gamma must be positive.")
    sigma = gamma / np.sqrt(8.0 * np.log(2.0))
    if n is None:
        n = int(np.ceil(sigma * np.sqrt(-2.0 * np.log(0.0002))) + 1)
    kernel_radius = min(
        int(np.ceil(sigma * np.sqrt(-2.0 * np.log(0.0002))) + 1),
        int(np.floor(n / 2)),
    )
    ii = np.arange(-kernel_radius, kernel_radius + 1, dtype=np.float64)
    rsf_x = 0.5 * (
        special.erf((ii + 0.5) / (np.sqrt(2.0) * sigma))
        - special.erf((ii - 0.5) / (np.sqrt(2.0) * sigma))
    )
    kernel = np.outer(rsf_x, rsf_x)
    return (kernel / kernel.sum()).astype(np.float64)


def airy_kernel(pixel_m: float, wavelength_m: float, na: float, z: float = 0.0, n: int = 17) -> np.ndarray:
    """Generate the SACDm scalar diffraction PSF approximation.

    The MATLAB code evaluates an integral independently for each radius. This
    uses fixed Gauss-Legendre quadrature over the same integrand, vectorized over
    all image radii.
    """

    if n < 33 and n > 17:
        nn = 8
    elif n < 65 and n > 33:
        nn = 16
    elif n < 129 and n > 65:
        nn = 32
    elif n < 257 and n > 129:
        nn = 64
    else:
        nn = 64

    x = np.arange(-nn, nn + 1, dtype=np.float64) * pixel_m
    xx, yy = np.meshgrid(x, x, indexing="xy")
    radius = np.hypot(xx, yy)
    inside = radius <= 1.0

    sin2 = (1.0 - (1.0 - na**2)) / 2.0
    u = 8.0 * np.pi * z * sin2 / wavelength_m

    nodes, weights = np.polynomial.legendre.leggauss(96)
    p = (nodes + 1.0) / 2.0
    w = weights / 2.0
    phase = np.exp((1j * u * p**2) / 2.0)

    arg = (2.0 * np.pi * radius[..., None] * na / wavelength_m) * p
    field = np.sum(2.0 * phase * special.j0(arg) * w, axis=-1)
    intensity = np.zeros_like(radius, dtype=np.float64)
    intensity[inside] = np.abs(field[inside] ** 2)
    total = intensity.sum()
    if total <= 0:
        raise ValueError("Generated PSF has nonpositive sum.")
    return intensity / total


def make_psfs(
    image_shape: tuple[int, int],
    pixel_nm: float,
    wavelength_nm: float,
    na: float,
    mag: int,
    psf: np.ndarray | None = None,
    resolution_nm: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Create pre- and post-deconvolution PSFs following `SACDm.m`."""

    min_dim = min(image_shape)
    if psf is not None:
        base = np.asarray(psf, dtype=np.float64)
        from .fourier import fourier_interpolate

        high = fourier_interpolate(base, (mag, mag), mirror_mode="lateral")
    elif resolution_nm is not None:
        base = generate_rsf(resolution_nm / pixel_nm, min_dim)
        high = generate_rsf(mag * resolution_nm / pixel_nm, min_dim)
    else:
        pixel_m = pixel_nm * 1e-9
        wavelength_m = wavelength_nm * 1e-9
        base = airy_kernel(pixel_m, wavelength_m, na, 0.0, min_dim)
        high = airy_kernel(pixel_m / mag, wavelength_m, na, 0.0, min_dim)
    return base / base.sum(), high / high.sum()
