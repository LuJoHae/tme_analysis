#!/usr/bin/env python3
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr, spearmanr

# Set aesthetic style
sns.set_theme(style="whitegrid")
plt.rcParams.update({"font.size": 11, "figure.autolayout": True})

OUTPUT_DIR = Path("output") / "visualize_hsic_dependency"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def compute_rbf_kernel(X, bandwidth_factor=1.0, gamma=None):
    """Compute RBF kernel matrix using median distance heuristic for gamma scaled by bandwidth_factor."""
    X = np.atleast_2d(X)
    if X.shape[0] == 1:
        X = X.T
    dists_sq = squareform(pdist(X, metric="sqeuclidean"))
    if gamma is None:
        dists = np.sqrt(dists_sq)
        median_dist = np.median(dists[dists > 0])
        if median_dist == 0 or np.isnan(median_dist):
            median_dist = 1.0
        sigma = median_dist * bandwidth_factor
        gamma = 1.0 / (2.0 * sigma**2)
    K = np.exp(-gamma * dists_sq)
    return K, gamma


def compute_hsic_details(X, Y, bandwidth_factor=1.0):
    """
    Compute total HSIC, sample-wise influence scores h_i, centered kernels Kc and Lc,
    and the pairwise interaction matrix C = Kc * Lc under a given bandwidth_factor.
    """
    X = np.ravel(X)[:, np.newaxis]
    Y = np.ravel(Y)[:, np.newaxis]
    N = len(X)

    K, gamma_x = compute_rbf_kernel(X, bandwidth_factor=bandwidth_factor)
    L, gamma_y = compute_rbf_kernel(Y, bandwidth_factor=bandwidth_factor)

    H = np.eye(N) - (1.0 / N) * np.ones((N, N))
    Kc = H @ K @ H
    Lc = H @ L @ H

    # Pairwise interaction matrix
    C = Kc * Lc

    # Total HSIC statistic
    hsic_stat = float(np.trace(Kc @ Lc) / ((N - 1) ** 2))

    # Sample-wise influence score: h_i = sum_j (Kc_ij * Lc_ij) / (N - 1)
    sample_influence = np.sum(C, axis=1) / (N - 1)

    return hsic_stat, sample_influence, Kc, Lc, C, gamma_x, gamma_y


def compute_rkhs_conditional_moments(X, Y, n_grid=200, bandwidth_factor=0.25, gamma_x=None):
    """
    Compute non-parametric RKHS conditional expectation E[Y | X = x]
    and conditional standard deviation SD(Y | X = x) across a smooth grid.
    bandwidth_factor scales the median heuristic bandwidth sigma = bandwidth_factor * sigma_median.
    """
    X_vec = np.ravel(X)
    Y_vec = np.ravel(Y)
    
    x_grid = np.linspace(np.min(X_vec), np.max(X_vec), n_grid)
    
    if gamma_x is None:
        dists_sq = squareform(pdist(X_vec[:, np.newaxis], metric="sqeuclidean"))
        dists = np.sqrt(dists_sq)
        med = np.median(dists[dists > 0])
        sigma = (med if med > 0 else 1.0) * bandwidth_factor
        gamma_x = 1.0 / (2.0 * sigma**2)

    cond_mean = np.zeros(n_grid)
    cond_sd = np.zeros(n_grid)

    for idx, x_val in enumerate(x_grid):
        w = np.exp(-gamma_x * (x_val - X_vec) ** 2)
        sum_w = np.sum(w)
        if sum_w == 0:
            sum_w = 1e-12
        w_norm = w / sum_w

        mu = np.sum(w_norm * Y_vec)
        var = np.sum(w_norm * (Y_vec - mu) ** 2)

        cond_mean[idx] = mu
        cond_sd[idx] = np.sqrt(var)

    return x_grid, cond_mean, cond_sd


def get_synthetic_datasets(N=800):
    """Generate 3 synthetic benchmarks: Rotated Diamond, Parabola, Sinusoid."""
    np.random.seed(42)
    # 1. Rotated Diamond (theta = 45 deg)
    Z1 = np.random.uniform(-1, 1, size=N)
    Z2 = np.random.uniform(-1, 1, size=N)
    theta = np.radians(45.0)
    R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    XY_rot = R @ np.vstack([Z1, Z2])
    X1, Y1 = XY_rot[0, :], XY_rot[1, :]

    # 2. Parabola Y = X^2
    np.random.seed(101)
    X2 = np.random.uniform(-1.5, 1.5, size=N)
    Y2 = X2**2 + np.random.normal(0, 0.25, size=N)

    # 3. Sinusoid Y = sin(2*pi*X)
    np.random.seed(202)
    X3 = np.random.uniform(-1.0, 1.0, size=N)
    Y3 = np.sin(2.0 * np.pi * X3) + np.random.normal(0, 0.2, size=N)

    datasets = [
        ("Rotated Diamond (θ=45°)", X1, Y1),
        ("Parabola (Y = X² + ε)", X2, Y2),
        ("Sinusoid (Y = sin(2πX) + ε)", X3, Y3)
    ]
    return datasets


def plot_sample_influence_scatter():
    """
    Method 1: Scatter plot color-coded by sample-wise HSIC influence scores h_i.
    4 Bandwidths (1.0x, 0.5x, 0.25x, 0.125x) x 3 Datasets (Diamond, Parabola, Sinusoid).
    """
    datasets = get_synthetic_datasets(N=800)
    bandwidth_factors = [1.0, 0.5, 0.25, 0.125]
    bw_labels = ["Standard (1.0 × σ_med)", "Intermediate (0.5 × σ_med)", "Fine Local (0.25 × σ_med)", "Ultra-Fine (1/8 × σ_med)"]

    fig, axes = plt.subplots(4, 3, figsize=(17, 16), sharex="col")

    for col_idx, (name, X, Y) in enumerate(datasets):
        r_val, _ = pearsonr(X, Y)
        for row_idx, (bw, bw_lbl) in enumerate(zip(bandwidth_factors, bw_labels)):
            ax = axes[row_idx, col_idx]
            hsic_val, h_influence, _, _, _, _, _ = compute_hsic_details(X, Y, bandwidth_factor=bw)
            
            sc = ax.scatter(X, Y, c=h_influence, cmap="magma", s=20, alpha=0.85)
            cbar = fig.colorbar(sc, ax=ax)
            cbar.ax.tick_params(labelsize=8)
            cbar.set_label("hᵢ", fontweight="bold", fontsize=9)

            if row_idx == 0:
                ax.set_title(f"{name}\nPearson r = {r_val:.3f} | {bw_lbl}\nHSIC = {hsic_val:.4f}", fontsize=11, fontweight="bold")
            else:
                ax.set_title(f"{bw_lbl} | HSIC = {hsic_val:.4f}", fontsize=10, fontweight="bold")

            if col_idx == 0:
                ax.set_ylabel("Y", fontweight="bold")
            if row_idx == 3:
                ax.set_xlabel("X", fontweight="bold")

    fig.suptitle("Sample-wise HSIC Influence Scores (hᵢ) across 4 Bandwidths and 3 Datasets\nHigh-contrast yellow/white points indicate observations driving statistical non-linear dependency", fontsize=14, fontweight="bold")
    out_path = OUTPUT_DIR / "hsic_influence_scatter_rotated.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_rkhs_conditional_curves():
    """
    Method 2: Non-parametric RKHS conditional expectation E[Y|X] & variance bands.
    4 Bandwidths (1.0x, 0.5x, 0.25x, 0.125x) x 3 Datasets (Diamond, Parabola, Sinusoid).
    """
    datasets = get_synthetic_datasets(N=800)
    bandwidth_factors = [1.0, 0.5, 0.25, 0.125]
    bw_labels = ["Standard (1.0 × σ_med)", "Intermediate (0.5 × σ_med)", "Fine Local (0.25 × σ_med)", "Ultra-Fine (1/8 × σ_med)"]
    colors = ["#e7298a", "#7570b3", "#d95f02", "#1b9e77"]

    fig, axes = plt.subplots(4, 3, figsize=(17, 16), sharex="col")

    for col_idx, (name, X, Y) in enumerate(datasets):
        for row_idx, (bw, bw_lbl, col) in enumerate(zip(bandwidth_factors, bw_labels, colors)):
            ax = axes[row_idx, col_idx]
            x_grid, mu, sd = compute_rkhs_conditional_moments(X, Y, bandwidth_factor=bw)

            ax.scatter(X, Y, color="#2b5c8f", alpha=0.3, s=15, label="Data")
            ax.plot(x_grid, mu, color=col, linewidth=2.5, label=f"E[Y|X] ({bw}x)")
            ax.fill_between(x_grid, mu - sd, mu + sd, color=col, alpha=0.22, label="±1 SD Envelope")

            if row_idx == 0:
                ax.set_title(f"{name}\n{bw_lbl}", fontsize=11, fontweight="bold")
            else:
                ax.set_title(f"{bw_lbl}", fontsize=10, fontweight="bold")

            if col_idx == 0:
                ax.set_ylabel("Y", fontweight="bold")
            if row_idx == 3:
                ax.set_xlabel("X", fontweight="bold")

            ax.legend(loc="upper right", fontsize=8)

    fig.suptitle("RKHS Kernel Conditional Expectation E[Y|X] & Variance Envelopes across 4 Bandwidths\nRows: Standard (1.0 × σ_med) | Intermediate (0.5 × σ_med) | Fine (0.25 × σ_med) | Ultra-Fine (1/8 × σ_med)", fontsize=14, fontweight="bold")
    out_path = OUTPUT_DIR / "hsic_rkhs_curves_synthetic.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_kernel_matrix_heatmap():
    """
    Method 3: Pairwise Centered Kernel Interaction Matrix C_ij = Kc_ij * Lc_ij Heatmap.
    4 Bandwidths (1.0x, 0.5x, 0.25x, 0.125x) x 3 Datasets (Diamond, Parabola, Sinusoid).
    """
    datasets = get_synthetic_datasets(N=120)  # N=120 for clear matrix resolution
    bandwidth_factors = [1.0, 0.5, 0.25, 0.125]
    bw_labels = ["Standard (1.0 × σ_med)", "Intermediate (0.5 × σ_med)", "Fine Local (0.25 × σ_med)", "Ultra-Fine (1/8 × σ_med)"]

    fig, axes = plt.subplots(4, 3, figsize=(16, 16))

    for col_idx, (name, X, Y) in enumerate(datasets):
        # Sort samples by X for structured matrix visualization
        sort_idx = np.argsort(X)
        X_s, Y_s = X[sort_idx], Y[sort_idx]

        for row_idx, (bw, bw_lbl) in enumerate(zip(bandwidth_factors, bw_labels)):
            ax = axes[row_idx, col_idx]
            hsic_val, _, _, _, C, _, _ = compute_hsic_details(X_s, Y_s, bandwidth_factor=bw)

            im = ax.imshow(C, cmap="magma", origin="lower")
            cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.ax.tick_params(labelsize=8)

            if row_idx == 0:
                ax.set_title(f"{name}\n{bw_lbl}\nHSIC = {hsic_val:.4f}", fontsize=11, fontweight="bold")
            else:
                ax.set_title(f"{bw_lbl} | HSIC = {hsic_val:.4f}", fontsize=10, fontweight="bold")

            if col_idx == 0:
                ax.set_ylabel("Samples (sorted by X)", fontweight="bold")
            if row_idx == 3:
                ax.set_xlabel("Samples (sorted by X)", fontweight="bold")

    fig.suptitle("Pairwise RKHS Kernel Alignment Matrix C = K_c ⊙ L_c across 4 Bandwidths and 3 Datasets\nHeatmaps display sample-pair RKHS similarity interactions (sum / (N-1)² = HSIC)", fontsize=14, fontweight="bold")
    out_path = OUTPUT_DIR / "hsic_kernel_matrix_heatmap.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved {out_path}")


def main():
    print("Starting HSIC dependency visualization script...")
    print(f"Output directory: {OUTPUT_DIR}")

    # Method 1: Sample-wise HSIC influence scatter plots (4x3)
    plot_sample_influence_scatter()

    # Method 2: RKHS kernel conditional expectation curves (4x3)
    plot_rkhs_conditional_curves()

    # Method 3: Kernel interaction matrix heatmaps (4x3)
    plot_kernel_matrix_heatmap()

    print("All tasks completed successfully!")


if __name__ == "__main__":
    main()
