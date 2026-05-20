# SMLM PCF-FFT: Fast Fourier Transform Pair Correlation Function for SMLM

This package provides a highly optimized, boundary-corrected **Pair Correlation Function (PCF)** calculation pipeline for Single Molecule Localization Microscopy (SMLM) data. 

By replacing slow geometric intersection calculations (e.g., using Shapely) with a **2D Fast Fourier Transform (FFT)** autocorrelation boundary correction method (inspired by the Veatch Lab SMLM package), this implementation achieves speedups of **100x to 10,000x** on large localization datasets.

---

## Features

- **Fast Boundary Correction**: Implements the Ohser-Stoyan isotropic boundary estimator utilizing 2D FFT autocorrelation of binary ROI masks.
- **Fast Pair Binning**: Uses `scipy.spatial.KDTree` with chunking to quickly count and bin distance pairs without memory bottlenecks.
- **Sliding Bins**: Supports overlapping (sliding) distance bins for smooth PCF curve profiling.
- **Automatic Nucleus Segmentation**: Integrates Cellpose to automatically detect cell/nucleus boundaries from a Gaussian-blurred 2D localization density image.
- **Auto and Cross-Correlation**: Performs both auto-correlation (single channel spatial pattern) and cross-correlation (spatial colocalization between two channels).

---

## Directory Structure

```
SMLM_PCF-FFT/
│
├── smlm_pcf_fft.py           # Core module containing FFT boundary correction and pair counting
├── SMLM_PCF_Pipeline.ipynb   # Interactive Jupyter Notebook for parameter setup and execution
│
├── validation/               # Folder containing validation resources
│   ├── generate_synthetic.py    # Script to generate synthetic point distributions
│   ├── run_synthetic_pcf.py     # Script to verify PCFs on synthetic data
│   ├── compare_pipelines.py     # Benchmark comparing FFT pipeline vs. Shapely pipeline
│   ├── test_synthetic_data/     # Synthetic CSV datasets and validation curves
│   └── test_dSTORM_data/        # Folder containing real sample dSTORM datasets and output results
│
└── README.md                 # Package documentation (this file)
```

---

## Prerequisites & Installation

The package runs in the `smlm` conda environment. To activate the environment:
```bash
conda activate smlm
```

Key dependencies:
- `numpy`, `pandas`, `scipy`
- `matplotlib`, `tqdm`
- `cellpose` (with GPU support if available)
- `shapely` (required only for running the old pipeline comparison benchmark)

---

## How to Use

### Running the Main Pipeline (Interactive Notebook)
Open and run [SMLM_PCF_Pipeline.ipynb](file:///Users/gmgao/GGscripts/GG-general-GuttmanLab/image_analysis/pipelines/SMLM_PCF-FFT/SMLM_PCF_Pipeline.ipynb) to configure parameters, run nuclei segmentation, and visualize results interactively. This is the main interface for adjusting parameters and plotting.
This script will:
- Load localization coordinate CSV files for the left (SPEN) and right (H3K27ac) channels.
- Reconstruct a 2D density image from the high-density channel and segment nucleus boundaries using Cellpose.
- Crop individual nuclei, extract localizations within each mask, and compute the Auto-Left, Auto-Right, and Cross PCF profiles.
- Save individual CSV tables for each nucleus and generate a pooled profile plot (`dstorm_pcf_curves.png`) saved under `validation/test_dSTORM_data/`.

### 2. Running Validation & Benchmarks
To regenerate synthetic datasets:
```bash
python validation/generate_synthetic.py
```

To run PCF calculations on synthetic datasets and verify matching physical curves (Random, Clustered, and Exclusion distributions):
```bash
python validation/run_synthetic_pcf.py
```

To run the speed comparison benchmark against the Shapely-based pipeline:
```bash
python validation/compare_pipelines.py
```

---

## Mathematical Summary

### Isotropic Edge Correction via 2D FFT
For an arbitrary binary mask $W$, the autocorrelation function $C(\mathbf{r}) = \int W(\mathbf{x}) W(\mathbf{x} + \mathbf{r}) d\mathbf{x}$ represents the overlap area of the mask with its shifted copy. The Ohser-Stoyan estimator weights each pair of points separated by displacement vector $\mathbf{r}$ by $1 / C(\mathbf{r})$.

We calculate $C(\mathbf{r})$ globally in Fourier space using the Wiener-Khinchin theorem:
$$C = \text{Re}\left(\mathcal{F}^{-1}\left( |\mathcal{F}(W)|^2 \right)\right)$$

To get the isotropic correction factor at distance $r$, we average $C(\mathbf{r})$ over multiple angles:
$$C_{\text{iso}}(r) = \frac{1}{2\pi} \int_{0}^{2\pi} C(r \cos\theta, r \sin\theta) d\theta$$
This is sampled efficiently in 2D using bilinear grid interpolation.
