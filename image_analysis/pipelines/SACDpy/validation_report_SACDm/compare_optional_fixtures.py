from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from skimage.metrics import normalized_root_mse

REPORT_DIR = Path(__file__).resolve().parent
REPO_ROOT = REPORT_DIR.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from sacdpy.background import background_estimation  # noqa: E402
from sacdpy.registration import register_image  # noqa: E402
from sacdpy.sparse_hessian import sparse_hessian_core  # noqa: E402


def metrics(a: np.ndarray, b: np.ndarray) -> dict[str, float | str]:
    aa = np.asarray(a, dtype=np.float64)
    bb = np.asarray(b, dtype=np.float64)
    aa = aa / max(float(np.max(np.abs(aa))), 1e-12)
    bb = bb / max(float(np.max(np.abs(bb))), 1e-12)
    diff = aa - bb
    return {
        "shape_python": "x".join(map(str, a.shape)),
        "shape_matlab": "x".join(map(str, b.shape)),
        "pearson": float(np.corrcoef(aa.ravel(), bb.ravel())[0, 1]),
        "nrmse": float(normalized_root_mse(bb, aa)),
        "mae": float(np.mean(np.abs(diff))),
        "max_abs_diff": float(np.max(np.abs(diff))),
    }


def main() -> int:
    fixture_dir = REPORT_DIR / "optional_reference"
    rows: list[dict[str, float | str]] = []

    bg = loadmat(fixture_dir / "background_fixture.mat")
    bg_py = background_estimation(bg["bg_input"], decomposition_level=4, iterations=2)
    rows.append({"branch": "background_estimation", **metrics(bg_py, bg["bg_output"])})

    reg = loadmat(fixture_dir / "registration_fixture.mat")
    reg_py = register_image(reg["reg_moving"], reg["reg_reference"])
    rows.append({"branch": "registration", **metrics(reg_py, reg["reg_output"])})

    sparse = loadmat(fixture_dir / "sparse_fixture.mat")
    sparse_py = sparse_hessian_core(sparse["sparse_input"], iterations=6)
    rows.append({"branch": "sparse_hessian", **metrics(sparse_py, sparse["sparse_output"])})

    out = REPORT_DIR / "metrics" / "optional_branch_metrics.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"{row['branch']}: pearson={row['pearson']:.6g} "
            f"nrmse={row['nrmse']:.6g} mae={row['mae']:.6g} max_abs_diff={row['max_abs_diff']:.6g}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
