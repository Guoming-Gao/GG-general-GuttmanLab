import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Ensure parent directory is in path to import smlm_pcf_fft
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from smlm_pcf_fft import calculate_pcf_fft

# Parameters
pixel_size = 100.0  # nm
x_min, y_min = 0.0, 0.0
width, height = 20000.0, 20000.0  # nm

# Create the logical mask (rectangular area of 200 x 200 pixels)
ny = int(height / pixel_size)
nx = int(width / pixel_size)
mask = np.ones((ny, nx), dtype=bool)

# Radial bins: 0 to 1000 nm with 20 nm bin size
r_edges = np.arange(0, 1001, 20)
bin_starts = r_edges[:-1]
bin_ends = r_edges[1:]

script_dir = os.path.dirname(os.path.abspath(__file__))
datasets = {
    'Random': os.path.join(script_dir, 'test_synthetic_data/random.csv'),
    'Enrichment (Clustered)': os.path.join(script_dir, 'test_synthetic_data/enrichment.csv'),
    'Exclusion (Hard-Core)': os.path.join(script_dir, 'test_synthetic_data/exclusion.csv')
}

plt.figure(figsize=(10, 6))

for name, path in datasets.items():
    print(f"Calculating PCF for {name}...")
    df = pd.read_csv(path)
    points = df[['x [nm]', 'y [nm]']].values
    
    # Calculate auto-correlation PCF
    G, r_centers, N_pairs, edge_cor = calculate_pcf_fft(
        points1=points,
        points2=None,
        mask=mask,
        pixel_size=pixel_size,
        x_min=x_min,
        y_min=y_min,
        bin_starts=bin_starts,
        bin_ends=bin_ends,
        mode='auto',
        show_progress=False
    )
    
    # Plot curve
    plt.plot(r_centers, G, label=f"{name}", linewidth=2)
    
    # Print some statistics
    print(f"  • Max G(r): {np.nanmax(G):.3f} at r={r_centers[np.nanargmax(G)]:.1f} nm")
    print(f"  • Mean G(r) for r > 500 nm: {np.nanmean(G[r_centers > 500]):.3f}")

plt.axhline(1.0, color='gray', linestyle='--', label='G(r) = 1 (No Correlation)')
plt.xlabel('Distance r (nm)', fontsize=12)
plt.ylabel('Pair Correlation Function G(r)', fontsize=12)
plt.title('Validation of FFT-based PCF on Synthetic Data', fontsize=14)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.3)

plot_filename = os.path.join(script_dir, 'test_synthetic_data/synthetic_pcf_curves.png')
plt.savefig(plot_filename, dpi=150)
print(f"Validation plot saved to {plot_filename}")
