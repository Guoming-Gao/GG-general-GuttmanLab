from __future__ import annotations

import argparse

from .params import SACDParams
from .reconstruction import reconstruct
from .tiffio import read_tiff_stack, write_tiff_image


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run core SACD reconstruction on a TIFF stack.")
    parser.add_argument("input", help="Input TIFF stack, read as TYX when loaded from tifffile.")
    parser.add_argument("output", help="Output float32 TIFF. Chunked outputs are written as TYX stacks.")
    parser.add_argument("--pixel", type=float, default=65.0, help="Pixel size in nm.")
    parser.add_argument("--wavelength", type=float, default=525.0, help="Emission wavelength in nm.")
    parser.add_argument("--na", type=float, default=1.45, help="Objective numerical aperture.")
    parser.add_argument("--mag", type=int, default=2, help="Lateral Fourier interpolation factor.")
    parser.add_argument("--iter1", type=int, default=7, help="Pre-RL iterations.")
    parser.add_argument("--iter2", type=int, default=8, help="Post-RL iterations.")
    parser.add_argument("--ac-order", type=int, default=2, help="Autocumulant order.")
    parser.add_argument(
        "--frames-per-sacd",
        type=int,
        default=None,
        help="Number of input frames per SACD reconstruction. Defaults to the full stack.",
    )
    parser.add_argument("--scale", type=float, default=None, help="Post-RL PSF exponent; defaults to AC order.")
    parser.add_argument("--subfactor", type=float, default=0.8, help="Mean subtraction factor before cumulant.")
    parser.add_argument("--background", action="store_true", help="Enable SACDm wavelet background subtraction.")
    parser.add_argument("--backgroundfactor", type=float, default=2.0, help="SACDm background weight.")
    parser.add_argument("--registration", action="store_true", help="Enable subpixel registration to the first frame.")
    parser.add_argument("--sparse-decon", action="store_true", help="Enable sparse Hessian post-deconvolution.")
    parser.add_argument("--fidelity", type=float, default=100.0, help="Sparse Hessian fidelity weight.")
    parser.add_argument("--tcontinuity", type=float, default=0.1, help="Sparse Hessian t/z continuity weight.")
    parser.add_argument("--sparsity", type=float, default=1.0, help="Sparse Hessian sparsity weight.")
    parser.add_argument("--sparse-iterations", type=int, default=100, help="Sparse Hessian iterations.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    params = SACDParams(
        pixel_nm=args.pixel,
        wavelength_nm=args.wavelength,
        na=args.na,
        mag=args.mag,
        iter1=args.iter1,
        iter2=args.iter2,
        ac_order=args.ac_order,
        frames_per_sacd=args.frames_per_sacd,
        scale=args.scale,
        subfactor=args.subfactor,
        ifbackground=args.background,
        backgroundfactor=args.backgroundfactor,
        ifregistration=args.registration,
        ifsparsedecon=args.sparse_decon,
        fidelity=args.fidelity,
        tcontinuity=args.tcontinuity,
        sparsity=args.sparsity,
        sparse_iterations=args.sparse_iterations,
    )
    stack = read_tiff_stack(args.input)
    result = reconstruct(stack, params)
    write_tiff_image(args.output, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
