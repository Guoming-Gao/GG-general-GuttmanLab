import os
import sys
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
from shapely.vectorized import contains

# Ensure parent directory is in path to import smlm_pcf_fft
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from smlm_pcf_fft import calculate_pcf_fft

# Parameters for comparison
nm_per_pxl = 23.4
pixel_size = 100.0  # 100 nm pixel size for mask
width, height = 20000.0, 20000.0  # nm

# Distance bins: 0 to 1000 nm with 100 nm bin width (similar to old pipeline config)
ringwidth_nm = 100.0
dr_slidingrings_nm = 20.0
bin_starts = np.arange(0, 1000.0 - ringwidth_nm, dr_slidingrings_nm)
bin_ends = bin_starts + ringwidth_nm
r_centers = 0.5 * (bin_starts + bin_ends)

# Precompute ring areas (in nm^2 and pixel^2 for old pipeline)
ring_areas_nm2 = np.pi * (bin_ends**2 - bin_starts**2)
ring_areas_pxl2 = ring_areas_nm2 / (nm_per_pxl**2)

# Create a polygon representation of a 20 um x 20 um box in pixels (for Shapely)
poly_coords_pxl = [
    (0, 0),
    (width / nm_per_pxl, 0),
    (width / nm_per_pxl, height / nm_per_pxl),
    (0, height / nm_per_pxl),
    (0, 0)
]
cell_polygon = Polygon(poly_coords_pxl)

# Create a mask representation of a 20 um x 20 um box (for FFT)
ny = int(height / pixel_size)
nx = int(width / pixel_size)
mask = np.ones((ny, nx), dtype=bool)

# Load synthetic data
print("Loading random synthetic data...")
script_dir = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(script_dir, 'test_synthetic_data/random.csv'))

# Use a small subset of points because Shapely is extremely slow
n_compare = 200
points_subset = df[['x [nm]', 'y [nm]']].values[:n_compare]

print(f"Comparing pipelines with {n_compare} points...")

# =============================================================================
# 1. OLD PIPELINE (SHAPELY)
# =============================================================================
print("Running Old Shapely-based pipeline...")
start_time_old = time.time()

# Coordinates in pixel units for Shapely
x_pxl = points_subset[:, 0] / nm_per_pxl
y_pxl = points_subset[:, 1] / nm_per_pxl

hist_results = []
# Loop over reference points (auto-correlation)
for i in range(len(x_pxl)):
    # Create rings for this point in pixel units
    rings = [Point(x_pxl[i], y_pxl[i]).buffer(end / nm_per_pxl).difference(Point(x_pxl[i], y_pxl[i]).buffer(start / nm_per_pxl))
             for start, end in zip(bin_starts, bin_ends)]
    
    # Calculate overlap areas with cell polygon
    intersect_areas = np.array([cell_polygon.intersection(ring).area for ring in rings])
    
    # Avoid division by zero
    valid_area = intersect_areas > 0
    edge_correction_factors = np.zeros_like(intersect_areas)
    edge_correction_factors[valid_area] = 1.0 / (intersect_areas[valid_area] / ring_areas_pxl2[valid_area])
    
    # Calculate distances to all other points
    distances = np.sqrt((x_pxl[i] - x_pxl) ** 2 + (y_pxl[i] - y_pxl) ** 2) * nm_per_pxl
    
    # Bin distances (excluding self-interaction)
    # Exclude distance == 0 (self-interaction)
    valid_dist = distances > 0
    hist_matrix = (bin_starts[:, np.newaxis] <= distances[valid_dist]) & (bin_ends[:, np.newaxis] > distances[valid_dist])
    hist_per_point = np.sum(hist_matrix, axis=1)
    
    # Apply edge correction
    hist_results.append(hist_per_point * edge_correction_factors)

# Total normalization
area_nm2 = cell_polygon.area * (nm_per_pxl**2)
rho_interest_per_nm2 = len(points_subset) / area_nm2
norm_factors = len(points_subset) * ring_areas_nm2 * rho_interest_per_nm2

pcf_old = np.sum(hist_results, axis=0) / norm_factors
elapsed_old = time.time() - start_time_old
print(f"Old pipeline finished in {elapsed_old:.3f} seconds.")

# =============================================================================
# 2. NEW PIPELINE (FFT + KDTREE)
# =============================================================================
print("Running New FFT-based pipeline...")
start_time_new = time.time()

# calculate_pcf_fft expects bin_starts and bin_ends
pcf_new, _, _, _ = calculate_pcf_fft(
    points1=points_subset,
    points2=None,
    mask=mask,
    pixel_size=pixel_size,
    x_min=0.0,
    y_min=0.0,
    bin_starts=bin_starts,
    bin_ends=bin_ends,
    mode='auto',
    show_progress=False
)

elapsed_new = time.time() - start_time_new
print(f"New pipeline finished in {elapsed_new:.5f} seconds.")
speedup = elapsed_old / elapsed_new
print(f"Speedup factor: {speedup:.1f}x")

# =============================================================================
# 3. PLOT COMPARISON
# =============================================================================
plt.figure(figsize=(10, 5))
plt.plot(r_centers, pcf_old, 'o-', label=f"Old Shapely Pipeline ({elapsed_old:.2f}s)", color='red', alpha=0.7)
plt.plot(r_centers, pcf_new, 's--', label=f"New FFT Pipeline ({elapsed_new:.4f}s)", color='blue', alpha=0.7)
plt.axhline(1.0, color='gray', linestyle=':')
plt.xlabel("Distance r (nm)")
plt.ylabel("G(r)")
plt.title(f"Comparison of Shapely vs FFT Pipelines ({n_compare} points)")
plt.legend()
plt.grid(True, alpha=0.3)

comparison_plot = os.path.join(script_dir, 'test_synthetic_data/pipeline_comparison.png')
plt.savefig(comparison_plot, dpi=150)
print(f"Comparison plot saved to {comparison_plot}")
