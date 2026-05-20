import numpy as np
from scipy.fft import fft2, ifft2, fftshift
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import KDTree
from tqdm import tqdm

def compute_fft_boundary_correction(mask, pixel_size, r_centers, ntheta=60):
    """
    Computes the isotropic boundary correction for a given logical mask using 2D FFT.
    
    Parameters:
    -----------
    mask : 2D array of bool
        Logical mask of the ROI (True/1 inside, False/0 outside).
    pixel_size : float
        Physical size of each pixel in nanometers (or same units as coordinates).
    r_centers : 1D array of float
        Radial distances (in nanometers) at which to evaluate the correction.
    ntheta : int
        Number of angular samples to average over for isotropic correction.
        
    Returns:
    --------
    edge_cor : 1D array of float
        Isotropic edge correction factor for each radial distance.
    """
    ny, nx = mask.shape
    sza = 2 * ny + 1
    szb = 2 * nx + 1
    
    # Compute 2D FFT autocorrelation of the binary mask
    mask_fft = fft2(mask, s=(sza, szb))
    C2 = np.real(fftshift(ifft2(np.abs(mask_fft)**2)))
    
    # Normalize by the total mask area in pixels
    mask_area_pixels = np.sum(mask)
    if mask_area_pixels > 0:
        C2 = C2 / mask_area_pixels
    else:
        return np.zeros_like(r_centers)
        
    # Set up coordinates for interpolation (centered at shift = 0, which is index (ny, nx))
    y_coords = np.arange(-ny, ny + 1) * pixel_size
    x_coords = np.arange(-nx, nx + 1) * pixel_size
    
    # Create the regular grid interpolator
    # First dimension is rows (y_coords), second is columns (x_coords)
    interp = RegularGridInterpolator((y_coords, x_coords), C2, method='linear', bounds_error=False, fill_value=0.0)
    
    # Sample at ntheta angles for each r
    thetas = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    
    # Vectorized calculation:
    # Outer product: rows are r, columns are theta
    dy = np.outer(r_centers, np.sin(thetas))
    dx = np.outer(r_centers, np.cos(thetas))
    
    # Flatten grid and query interpolator
    query_points = np.stack([dy.ravel(), dx.ravel()], axis=1)
    vals = interp(query_points)
    
    # Reshape back and average across theta
    edge_cor = np.mean(vals.reshape(len(r_centers), ntheta), axis=1)
    return edge_cor

def filter_points_in_mask(points, mask, x_min, y_min, pixel_size):
    """
    Filters coordinates to keep only those lying inside the True region of the mask.
    """
    col_indices = ((points[:, 0] - x_min) / pixel_size).astype(int)
    row_indices = ((points[:, 1] - y_min) / pixel_size).astype(int)
    
    in_bounds = (col_indices >= 0) & (col_indices < mask.shape[1]) & (row_indices >= 0) & (row_indices < mask.shape[0])
    
    valid = np.zeros(len(points), dtype=bool)
    if len(points[in_bounds]) > 0:
        valid[in_bounds] = mask[row_indices[in_bounds], col_indices[in_bounds]]
        
    return points[valid]

def calculate_pcf_fft(points1, points2, mask, pixel_size, x_min, y_min, bin_starts, bin_ends, mode='auto', chunk_size=20000, show_progress=False):
    """
    Calculates the Pair Correlation Function G(r) with FFT-based edge correction.
    Supports overlapping/sliding bins.
    
    Parameters:
    -----------
    points1 : 2D array of shape (N1, 2)
        Coordinates of set 1 [x, y] in nanometers.
    points2 : 2D array of shape (N2, 2) or None
        Coordinates of set 2 [x, y] in nanometers (only needed if mode='cross').
    mask : 2D array of bool
        Logical mask of the ROI (True/1 inside, False/0 outside).
    pixel_size : float
        Physical pixel size of the mask in nanometers.
    x_min : float
        Minimum x coordinate of the mask bounding box.
    y_min : float
        Minimum y coordinate of the mask bounding box.
    bin_starts : 1D array of float
        Start distances for each bin in nanometers.
    bin_ends : 1D array of float
        End distances for each bin in nanometers.
    mode : str
        'auto' for auto-correlation, 'cross' for cross-correlation.
    chunk_size : int
        Chunk size for KDTree queries.
    show_progress : bool
        If True, displays progress bars.
        
    Returns:
    --------
    G : 1D array of float
        The calculated Pair Correlation Function values for each bin.
    r_centers : 1D array of float
        The centers of the distance bins in nanometers.
    N_pairs : 1D array of int
        Raw pair counts in each bin.
    edge_cor : 1D array of float
        Isotropic edge correction factor for each bin.
    """
    r_centers = 0.5 * (bin_starts + bin_ends)
    r_max = np.max(bin_ends)
    
    # 1. Filter points to keep only those within the mask
    pts1_filtered = filter_points_in_mask(points1, mask, x_min, y_min, pixel_size)
    if mode == 'cross':
        if points2 is None:
            raise ValueError("points2 must be provided for cross-correlation mode")
        pts2_filtered = filter_points_in_mask(points2, mask, x_min, y_min, pixel_size)
    
    N1 = len(pts1_filtered)
    N2 = len(pts2_filtered) if mode == 'cross' else N1
    
    if N1 == 0 or N2 == 0:
        return np.zeros(len(r_centers)), r_centers, np.zeros(len(r_centers), dtype=int), np.zeros(len(r_centers))
        
    # 2. Compute 2D FFT-based boundary correction
    if show_progress:
        print("Computing FFT boundary correction...")
    edge_cor = compute_fft_boundary_correction(mask, pixel_size, r_centers)
    
    # 3. Compute raw binned pair counts
    N_pairs = np.zeros(len(bin_starts), dtype=np.int64)
    chunk_size = 5000
    
    if mode == 'auto':
        tree = KDTree(pts1_filtered)
        for start_idx in range(0, len(pts1_filtered), chunk_size):
            end_idx = min(start_idx + chunk_size, len(pts1_filtered))
            chunk = pts1_filtered[start_idx:end_idx]
            indices = tree.query_ball_point(chunk, r_max)
            
            lengths = np.array([len(x) for x in indices])
            flat_i = np.repeat(np.arange(len(chunk)), lengths)
            non_empty = [np.array(x, dtype=np.int64) for x in indices if len(x) > 0]
            if non_empty:
                flat_neighbors = np.concatenate(non_empty)
                actual_i_arr = start_idx + flat_i
                valid = flat_neighbors > actual_i_arr
                if np.any(valid):
                    v_neighbors = flat_neighbors[valid]
                    v_i = flat_i[valid]
                    diff = pts1_filtered[v_neighbors] - chunk[v_i]
                    dists = np.linalg.norm(diff, axis=1)
                    for b_idx, (bs, be) in enumerate(zip(bin_starts, bin_ends)):
                        N_pairs[b_idx] += np.sum((dists >= bs) & (dists < be))
                        
        # Multiply by 2 to account for both directions (i to j and j to i)
        N_pairs = N_pairs * 2
        
    else:
        tree2 = KDTree(pts2_filtered)
        for start_idx in range(0, len(pts1_filtered), chunk_size):
            end_idx = min(start_idx + chunk_size, len(pts1_filtered))
            chunk = pts1_filtered[start_idx:end_idx]
            indices = tree2.query_ball_point(chunk, r_max)
            
            lengths = np.array([len(x) for x in indices])
            flat_i = np.repeat(np.arange(len(chunk)), lengths)
            non_empty = [np.array(x, dtype=np.int64) for x in indices if len(x) > 0]
            if non_empty:
                flat_neighbors = np.concatenate(non_empty)
                diff = pts2_filtered[flat_neighbors] - chunk[flat_i]
                dists = np.linalg.norm(diff, axis=1)
                for b_idx, (bs, be) in enumerate(zip(bin_starts, bin_ends)):
                    N_pairs[b_idx] += np.sum((dists >= bs) & (dists < be))
                    
    # 4. Calculate normalization factor
    area = np.sum(mask) * (pixel_size ** 2)
    ring_areas_nm2 = np.pi * (bin_ends**2 - bin_starts**2)
    
    if mode == 'auto':
        density = N1 / area
        basic_norm = area * (density ** 2) * ring_areas_nm2
    else:
        density1 = N1 / area
        density2 = N2 / area
        basic_norm = area * density1 * density2 * ring_areas_nm2
        
    # Combine basic normalization and edge correction
    Norm = basic_norm * edge_cor
    
    # Compute G(r)
    G = np.zeros_like(r_centers, dtype=float)
    valid = Norm > 0
    G[valid] = N_pairs[valid] / Norm[valid]
    G[~valid] = np.nan
    
    return G, r_centers, N_pairs, edge_cor

def is_valid_nucleus(mask, label, min_area_pixels=6000, edge_margin_pixels=5):
    """
    Checks if a nucleus is valid based on area and distance from the mask boundary.
    
    Parameters:
    -----------
    mask : 2D array of int
        Segmented label mask.
    label : int
        Label ID of the nucleus to check.
    min_area_pixels : int
        Minimum area threshold in pixels.
    edge_margin_pixels : int
        Minimum boundary margin in pixels (nucleus must not touch this close to the border).
        
    Returns:
    --------
    is_valid : bool
        True if the nucleus is valid, False otherwise.
    """
    y_indices, x_indices = np.where(mask == label)
    if len(y_indices) == 0:
        return False
        
    # Size check
    area = len(y_indices)
    if area < min_area_pixels:
        return False
        
    # Edge margin check
    y_min, y_max = y_indices.min(), y_indices.max()
    x_min, x_max = x_indices.min(), x_indices.max()
    
    ny, nx = mask.shape
    if (y_min < edge_margin_pixels or 
        y_max >= ny - edge_margin_pixels or 
        x_min < edge_margin_pixels or 
        x_max >= nx - edge_margin_pixels):
        return False
        
    return True

def process_single_nucleus_pcf(args):
    """
    Worker function for parallel processing of a single nucleus PCF calculation.
    
    Parameters:
    -----------
    args : tuple
        (fov_name, label, pts_left, pts_right, nuc_mask, pixel_size, 
         nuc_x_min_phys, nuc_y_min_phys, bin_starts, bin_ends, 
         min_points_threshold, output_dir)
         
    Returns:
    --------
    result : dict
        Summary dictionary of the nucleus execution results.
    """
    (fov_name, label, pts_left, pts_right, nuc_mask, pixel_size, 
     nuc_x_min_phys, nuc_y_min_phys, bin_starts, bin_ends, 
     min_points_threshold, output_dir) = args
     
    # Filter points to keep only those inside the cropped mask
    nuc_pts_left = filter_points_in_mask(pts_left, nuc_mask, nuc_x_min_phys, nuc_y_min_phys, pixel_size)
    nuc_pts_right = filter_points_in_mask(pts_right, nuc_mask, nuc_x_min_phys, nuc_y_min_phys, pixel_size)
    
    n_left = len(nuc_pts_left)
    n_right = len(nuc_pts_right)
    
    if n_left < min_points_threshold or n_right < min_points_threshold:
        return {
            'fov_name': fov_name,
            'label': label,
            'status': 'skipped_insufficient_points',
            'n_left': n_left,
            'n_right': n_right
        }
        
    # Auto-Correlation Left
    G_left, r_centers, _, _ = calculate_pcf_fft(
        nuc_pts_left, None, nuc_mask, pixel_size, nuc_x_min_phys, nuc_y_min_phys, bin_starts, bin_ends, mode='auto'
    )
    
    # Auto-Correlation Right
    G_right, _, _, _ = calculate_pcf_fft(
        nuc_pts_right, None, nuc_mask, pixel_size, nuc_x_min_phys, nuc_y_min_phys, bin_starts, bin_ends, mode='auto'
    )
    
    # Cross-Correlation
    G_cross, _, _, _ = calculate_pcf_fft(
        nuc_pts_left, nuc_pts_right, nuc_mask, pixel_size, nuc_x_min_phys, nuc_y_min_phys, bin_starts, bin_ends, mode='cross'
    )
    
    # Save individual results to CSV
    import pandas as pd
    import os
    
    nuc_df = pd.DataFrame({
        'r_center_nm': r_centers,
        'bin_start_nm': bin_starts,
        'bin_end_nm': bin_ends,
        'G_auto_left': G_left,
        'G_auto_right': G_right,
        'G_cross': G_cross,
    })
    
    csv_filename = f"{fov_name}_nucleus_{label}_pcf.csv"
    nuc_csv_path = os.path.join(output_dir, csv_filename)
    nuc_df.to_csv(nuc_csv_path, index=False)
    
    return {
        'fov_name': fov_name,
        'label': label,
        'status': 'success',
        'n_left': n_left,
        'n_right': n_right,
        'G_left': G_left.tolist(),
        'G_right': G_right.tolist(),
        'G_cross': G_cross.tolist(),
        'csv_path': nuc_csv_path
    }

