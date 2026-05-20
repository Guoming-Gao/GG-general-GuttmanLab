import os
import numpy as np
import pandas as pd

# Parameters for synthetic data
width = 20000.0   # nm (20 um)
height = 20000.0  # nm (20 um)
num_points = 5000 # target number of points

script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(script_dir, 'test_synthetic_data')
os.makedirs(output_dir, exist_ok=True)

# 1. Random Distribution (Uniform Poisson Point Process)
print("Generating Random distribution...")
x_random = np.random.uniform(0, width, num_points)
y_random = np.random.uniform(0, height, num_points)

df_random = pd.DataFrame({'x [nm]': x_random, 'y [nm]': y_random})
df_random.to_csv(os.path.join(output_dir, 'random.csv'), index=False)

# 2. Clustered (Enrichment) Distribution (Neyman-Scott Process)
print("Generating Clustered (Enrichment) distribution...")
num_parents = 100
parent_x = np.random.uniform(0, width, num_parents)
parent_y = np.random.uniform(0, height, num_parents)

# Mean number of offspring per parent
mean_offspring = num_points // num_parents
# standard deviation of cluster (in nm)
sigma_cluster = 150.0  # 150 nm standard deviation

offspring_x = []
offspring_y = []
for px, py in zip(parent_x, parent_y):
    n_offspring = np.random.poisson(mean_offspring)
    ox = np.random.normal(px, sigma_cluster, n_offspring)
    oy = np.random.normal(py, sigma_cluster, n_offspring)
    
    # Clip to boundary
    ox = np.clip(ox, 0, width)
    oy = np.clip(oy, 0, height)
    
    offspring_x.extend(ox)
    offspring_y.extend(oy)

df_enrichment = pd.DataFrame({'x [nm]': offspring_x, 'y [nm]': offspring_y})
df_enrichment.to_csv(os.path.join(output_dir, 'enrichment.csv'), index=False)

# 3. Exclusion Distribution (Matérn Hard-Core Process)
print("Generating Exclusion distribution...")
core_radius = 100.0 # nm exclusion zone
# Generate parent points with higher density, then keep only those that don't violate exclusion
factor = 5
large_num = num_points * factor
x_large = np.random.uniform(0, width, large_num)
y_large = np.random.uniform(0, height, large_num)

keep_x = []
keep_y = []

# Using a KDTree for fast exclusion checking
from scipy.spatial import KDTree
points_large = np.stack([x_large, y_large], axis=1)
tree = KDTree(points_large)

# Greedy selection
selected = np.zeros(large_num, dtype=bool)
for i in range(large_num):
    if selected[i]:
        continue
    # Find all neighbors within exclusion distance
    neighbors = tree.query_ball_point(points_large[i], core_radius)
    # Mark all neighbors as discarded (except itself)
    for n in neighbors:
        if n != i:
            selected[n] = True
    keep_x.append(points_large[i, 0])
    keep_y.append(points_large[i, 1])
    
    if len(keep_x) >= num_points:
        break

df_exclusion = pd.DataFrame({'x [nm]': keep_x, 'y [nm]': keep_y})
df_exclusion.to_csv(os.path.join(output_dir, 'exclusion.csv'), index=False)

print("Synthetic datasets generated and saved:")
print(f"  • Random: {len(df_random)} points")
print(f"  • Enrichment: {len(df_enrichment)} points")
print(f"  • Exclusion: {len(df_exclusion)} points")
