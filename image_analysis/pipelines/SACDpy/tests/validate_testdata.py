from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import tifffile
from skimage.metrics import normalized_root_mse, structural_similarity

from sacdpy import SACDParams, reconstruct, write_tiff_image


def validation_metrics(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    aa = a.astype(np.float64)
    bb = b.astype(np.float64)
    aa = aa / max(float(aa.max()), 1.0)
    bb = bb / max(float(bb.max()), 1.0)
    data_range = max(float(bb.max() - bb.min()), 1e-12)
    return {
        "corr": float(np.corrcoef(aa.ravel(), bb.ravel())[0, 1]),
        "nrmse": float(normalized_root_mse(bb, aa)),
        "mae": float(np.mean(np.abs(bb - aa))),
        "ssim": float(structural_similarity(bb, aa, data_range=data_range)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate SACDpy against bundled ImageJ SACD TIFFs.")
    test_dir = Path(__file__).resolve().parent
    parser.add_argument("--data-dir", type=Path, default=test_dir / "testdata")
    parser.add_argument("--out-dir", type=Path, default=test_dir / "outputs")
    parser.add_argument("--na", type=float, default=1.45)
    parser.add_argument("--iter1", type=int, default=7)
    parser.add_argument("--iter2", type=int, default=8)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    cases = [
        ("left", 560.0),
        ("right", 647.0),
    ]
    for side, wavelength in cases:
        raw = args.data_dir / f"TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-{side}.tif"
        ref = args.data_dir / f"TXSHA-undiff-SPEN_JF549_H3K27ac_DL650-SACD-FOV-1-SACD-{side}.tif"
        result = reconstruct(
            tifffile.imread(raw),
            SACDParams(
                pixel_nm=117.0,
                wavelength_nm=wavelength,
                na=args.na,
                iter1=args.iter1,
                iter2=args.iter2,
            ),
        )
        out = args.out_dir / f"sacdpy-{side}.tif"
        write_tiff_image(out, result)
        metrics = validation_metrics(result, tifffile.imread(ref))
        print(
            f"{side}: output={out} shape={result.shape} "
            f"corr={metrics['corr']:.4f} ssim={metrics['ssim']:.4f} "
            f"nrmse={metrics['nrmse']:.4f} mae={metrics['mae']:.4f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
