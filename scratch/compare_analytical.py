import sys
import os
import numpy as np
from sklearn.cluster import KMeans

# Add selective_inference directory to path
sys.path.append("scripts/selective_inference")
from kmeans_inference import KMeansSelectiveInference
from kmeans_selective_inferenece_analytically_and_simulation import compute_analytical_limits, generate_mixture_data

# Generate data
X, _ = generate_mixture_data(n_samples_comp1=75, n_samples_comp2=75, mean_diff=0.0, sigma=1.0, random_state=42)

# Fit KMeans
kmeans = KMeans(n_clusters=2, random_state=42, n_init=10)
kmeans.fit(X)
labels = kmeans.labels_
centers = kmeans.cluster_centers_

# 1. Repository implementation (kmeans_inference.py)
kmeans_si = KMeansSelectiveInference(n_clusters=2, random_state=42)
kmeans_si.labels_ = labels
kmeans_si.cluster_centers_ = centers
V_min_rep, V_max_rep, phi_obs_rep, nu_norm_rep = kmeans_si._compute_v_min_max(X, c1=0, c2=1)

# 2. Our implementation
V_min_our, V_max_our, phi_obs_our, v_norm_our = compute_analytical_limits(X, labels, centers, c1=0, c2=1, sigma=1.0)

print("Repository Implementation (kmeans_inference.py):")
print(f"  phi_obs:  {phi_obs_rep:.6f}")
print(f"  nu_norm:  {nu_norm_rep:.6f}")
print(f"  V_min:    {V_min_rep:.6f}")
print(f"  V_max:    {V_max_rep:.6f}")

print("\nOur Implementation:")
print(f"  phi_obs:  {phi_obs_our:.6f}")
print(f"  v_norm:   {v_norm_our:.6f}")
print(f"  V_min:    {V_min_our:.6f}")
print(f"  V_max:    {V_max_our:.6f}")

# Check scaling conversion
# In the repository:
#   std_phi = sigma / nu_norm
#   z_val = phi_obs / std_phi = phi_obs * nu_norm / sigma
# Let's see what p-values each computes:
naive_p_rep, selective_p_rep, _ = kmeans_si.test_difference_in_means(X, c1=0, c2=1, sigma=1.0)
print(f"\nRepository p-values: naive={naive_p_rep:.6f}, selective={selective_p_rep:.6f}")

from kmeans_selective_inferenece_analytically_and_simulation import compute_naive_and_analytical_p
naive_p_our, selective_p_our, _, _, _, _ = compute_naive_and_analytical_p(X, labels, centers, c1=0, c2=1, sigma=1.0)
print(f"Our p-values:        naive={naive_p_our:.6f}, selective={selective_p_our:.6f}")
