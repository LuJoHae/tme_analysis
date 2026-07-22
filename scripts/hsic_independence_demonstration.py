#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr, spearmanr, gamma

# Set aesthetic style
sns.set_theme(style="whitegrid")
plt.rcParams.update({"font.size": 11, "figure.autolayout": True})

OUTPUT_DIR = os.path.join("output", "hsic_independence_demonstration")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def compute_rbf_kernel(X, gamma=None):
    """
    Compute RBF (Gaussian) Kernel matrix for 1D or 2D inputs.
    Uses median heuristic for gamma if not provided.
    """
    X = np.atleast_2d(X)
    if X.shape[0] == 1:
        X = X.T

    dists_sq = squareform(pdist(X, metric="sqeuclidean"))

    if gamma is None:
        dists = np.sqrt(dists_sq)
        median_dist = np.median(dists[dists > 0])
        if median_dist == 0 or np.isnan(median_dist):
            median_dist = 1.0
        sigma = median_dist
        gamma = 1.0 / (2.0 * sigma**2)

    K = np.exp(-gamma * dists_sq)
    return K, gamma


def compute_hsic(X, Y, n_permutations=500, seed=None):
    """
    Compute Hilbert-Schmidt Independence Criterion (HSIC) between X and Y.
    Computes empirical p-value via permutation test and asymptotic p-value via Gamma approximation.
    """
    if seed is not None:
        np.random.seed(seed)

    X = np.ravel(X)[:, np.newaxis]
    Y = np.ravel(Y)[:, np.newaxis]
    N = len(X)

    K, _ = compute_rbf_kernel(X)
    L, _ = compute_rbf_kernel(Y)

    # Centering matrix H = I - 1/N * 11^T
    H = np.eye(N) - (1.0 / N) * np.ones((N, N))

    # Centered kernels
    Kc = H @ K @ H
    Lc = H @ L @ H

    # Empirical HSIC statistic
    hsic_stat = float(np.trace(Kc @ Lc) / ((N - 1) ** 2))

    # Permutation test using fast matrix trace property Tr(Kc L_perm)
    if n_permutations > 0:
        perm_stats = np.zeros(n_permutations)
        for b in range(n_permutations):
            perm_idx = np.random.permutation(N)
            perm_stats[b] = np.sum(Kc * L[perm_idx, :][:, perm_idx]) / ((N - 1) ** 2)

        p_val_perm = float((1 + np.sum(perm_stats >= hsic_stat)) / (n_permutations + 1))
        
        # Asymptotic Gamma approximation under H0
        mu_null = float(np.mean(perm_stats))
        var_null = float(np.var(perm_stats))
        if var_null > 0 and hsic_stat > 0:
            shape_k = (mu_null ** 2) / var_null
            scale_theta = var_null / mu_null
            p_val_gamma = float(gamma.sf(hsic_stat, a=shape_k, scale=scale_theta))
        else:
            p_val_gamma = 0.0
    else:
        p_val_perm = np.nan
        p_val_gamma = np.nan

    return hsic_stat, p_val_perm, p_val_gamma


def generate_rotated_data(n_samples=500, angle_degrees=45.0, seed=42):
    """
    Generate independent uniform variables Z1, Z2 ~ U(-1, 1) and rotate them by angle_degrees.
    """
    np.random.seed(seed)
    Z1 = np.random.uniform(-1, 1, size=n_samples)
    Z2 = np.random.uniform(-1, 1, size=n_samples)
    Z = np.vstack([Z1, Z2])

    theta = np.radians(angle_degrees)
    R = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta),  np.cos(theta)]
    ])

    XY = R @ Z
    X, Y = XY[0, :], XY[1, :]
    return X, Y


def plot_single_angle_comparison(n_samples=3000, seed=42):
    """
    Plot scatter plots comparing independent unrotated data vs rotated data (theta=45 deg).
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # Panel A: Unrotated independent square (theta = 0)
    X0, Y0 = generate_rotated_data(n_samples=n_samples, angle_degrees=0.0, seed=seed)
    r0, p_r0 = pearsonr(X0, Y0)
    rho0, p_rho0 = spearmanr(X0, Y0)
    hsic0, p_perm0, p_gam0 = compute_hsic(X0, Y0, n_permutations=300, seed=seed)

    axes[0].scatter(X0, Y0, color="#2b5c8f", alpha=0.35, edgecolors="none", s=12)
    axes[0].set_title(f"Independent Uniform Variables (θ = 0°, N={n_samples})\nStrictly Independent (P(X,Y) = P(X)P(Y))", fontsize=12, fontweight="bold")
    axes[0].set_xlabel("X")
    axes[0].set_ylabel("Y")
    axes[0].set_xlim(-1.6, 1.6)
    axes[0].set_ylim(-1.6, 1.6)
    
    stats_text0 = (
        f"Pearson r = {r0:.3f} (p={p_r0:.3f})\n"
        f"Spearman ρ = {rho0:.3f} (p={p_rho0:.3f})\n"
        f"HSIC = {hsic0:.4f} (p_perm={p_perm0:.3f}, p_gamma={p_gam0:.3f})"
    )
    axes[0].text(
        0.04, 0.94, stats_text0, transform=axes[0].transAxes,
        verticalalignment="top", fontsize=10.5,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.85, edgecolor="#cccccc")
    )

    # Panel B: Rotated Diamond (theta = 45 degrees)
    X45, Y45 = generate_rotated_data(n_samples=n_samples, angle_degrees=45.0, seed=seed)
    r45, p_r45 = pearsonr(X45, Y45)
    rho45, p_rho45 = spearmanr(X45, Y45)
    hsic45, p_perm45, p_gam45 = compute_hsic(X45, Y45, n_permutations=300, seed=seed)

    axes[1].scatter(X45, Y45, color="#d95f02", alpha=0.35, edgecolors="none", s=12)
    axes[1].set_title(f"Rotated Uniform Variables (θ = 45°, N={n_samples})\nNon-linearly Dependent (Pearson r = 0, HSIC > 0)", fontsize=12, fontweight="bold")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")
    axes[1].set_xlim(-1.6, 1.6)
    axes[1].set_ylim(-1.6, 1.6)

    gam_str = f"{p_gam45:.1e}" if p_gam45 > 0 else "< 1e-15"
    stats_text45 = (
        f"Pearson r = {r45:.3f} (p={p_r45:.3f})  [FAILS]\n"
        f"Spearman ρ = {rho45:.3f} (p={p_rho45:.3f})  [FAILS]\n"
        f"HSIC = {hsic45:.4f}  [DETECTS DEPENDENCE!]\n"
        f"  • Permutation p = {p_perm45:.3f} (floor for B=300)\n"
        f"  • Asymptotic p = {gam_str}"
    )
    axes[1].text(
        0.04, 0.94, stats_text45, transform=axes[1].transAxes,
        verticalalignment="top", fontsize=10.0, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#ffffcc", alpha=0.9, edgecolor="#d95f02")
    )

    fig.suptitle("Linear/Monotonic Correlation vs. Hilbert-Schmidt Independence Criterion (HSIC)", fontsize=14, fontweight="bold")
    out_path = os.path.join(OUTPUT_DIR, "single_angle_dependence_scatter.png")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_angle_sweep(n_samples=2500, n_angles=91, seed=42):
    """
    Vary rotation angle theta from 0 to 180 degrees and record Pearson r, Spearman rho, and HSIC.
    """
    angles = np.linspace(0, 180, n_angles)
    pearson_r_vals = []
    spearman_rho_vals = []
    hsic_vals = []

    for angle in angles:
        X, Y = generate_rotated_data(n_samples=n_samples, angle_degrees=angle, seed=seed)
        r, _ = pearsonr(X, Y)
        rho, _ = spearmanr(X, Y)
        hsic, _, _ = compute_hsic(X, Y, n_permutations=0)  # No permutation needed for smooth trend curve

        pearson_r_vals.append(r)
        spearman_rho_vals.append(rho)
        hsic_vals.append(hsic)

    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot Correlations on Primary Y Axis
    ax1.plot(angles, pearson_r_vals, "-", color="#2b5c8f", linewidth=2.5, label="Pearson Correlation r")
    ax1.plot(angles, spearman_rho_vals, "--", color="#7570b3", linewidth=2.2, label="Spearman Correlation ρ")
    ax1.axhline(0, color="gray", linestyle=":", linewidth=1)
    ax1.axvline(45, color="#d95f02", linestyle="--", linewidth=1.5, alpha=0.7)
    ax1.axvline(135, color="#d95f02", linestyle="--", linewidth=1.5, alpha=0.7)
    
    ax1.set_xlabel("Rotation Angle θ (degrees)")
    ax1.set_ylabel("Linear / Monotonic Correlation (r, ρ)", color="#2b5c8f", fontweight="bold")
    ax1.tick_params(axis="y", labelcolor="#2b5c8f")
    ax1.set_xlim(0, 180)
    ax1.set_xticks(np.arange(0, 181, 30))

    # Plot HSIC on Secondary Y Axis
    ax2 = ax1.twinx()
    ax2.plot(angles, hsic_vals, "-", color="#d95f02", linewidth=2.8, label="HSIC Statistic")
    ax2.set_ylabel("Hilbert-Schmidt Independence Criterion (HSIC)", color="#d95f02", fontweight="bold")
    ax2.tick_params(axis="y", labelcolor="#d95f02")

    # Combined Legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right", frameon=True)

    plt.title(f"Dependence Metrics vs. Rotation Angle θ of Uniform Bivariate Data (N={n_samples})\nNotice HSIC peaks at θ=45° and 135° where Correlations are EXACTLY ZERO", fontsize=13, fontweight="bold")
    
    out_path = os.path.join(OUTPUT_DIR, "hsic_vs_correlations_angle_sweep.png")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_grid_across_angles(n_samples=2000, seed=42):
    """
    Plot a grid of scatter plots across representative rotation angles (0, 15, 30, 45, 60, 90 deg).
    """
    angles = [0.0, 15.0, 30.0, 45.0, 60.0, 90.0]
    fig, axes = plt.subplots(2, 3, figsize=(15, 9.5), sharex=True, sharey=True)
    axes = axes.ravel()

    for idx, angle in enumerate(angles):
        ax = axes[idx]
        X, Y = generate_rotated_data(n_samples=n_samples, angle_degrees=angle, seed=seed)
        r, _ = pearsonr(X, Y)
        rho, _ = spearmanr(X, Y)
        hsic, p_perm, p_gam = compute_hsic(X, Y, n_permutations=150, seed=seed)

        color = "#d95f02" if abs(r) < 0.1 and hsic > 0.005 else "#2b5c8f"
        ax.scatter(X, Y, color=color, alpha=0.35, edgecolors="none", s=12)
        ax.set_title(f"θ = {angle:.0f}° (N={n_samples})", fontsize=12, fontweight="bold")
        ax.set_xlim(-1.6, 1.6)
        ax.set_ylim(-1.6, 1.6)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")

        p_str = f"{p_gam:.1e}" if p_gam > 0 else f"{p_perm:.3f}"
        annot = f"r = {r:+.2f}\nρ = {rho:+.2f}\nHSIC = {hsic:.4f}\np_asympt = {p_str}"
        ax.text(
            0.05, 0.93, annot, transform=ax.transAxes,
            verticalalignment="top", fontsize=9.5,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.85, edgecolor="#cccccc")
        )

    fig.suptitle(f"Scatter Plots and Independence Metrics across Rotation Angles θ (N={n_samples})", fontsize=15, fontweight="bold")
    out_path = os.path.join(OUTPUT_DIR, "multiple_angles_scatter_grid.png")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Saved {out_path}")


def main():
    print("Starting Hilbert-Schmidt Independence Criterion (HSIC) demonstration...")
    print(f"Output directory: {OUTPUT_DIR}")

    # Plot 1: Single angle comparison (0 vs 45 degrees)
    plot_single_angle_comparison(n_samples=3000, seed=42)

    # Plot 2: Angle sweep (0 to 180 degrees)
    plot_angle_sweep(n_samples=2500, n_angles=91, seed=42)

    # Plot 3: Scatter plot grid across key rotation angles
    plot_grid_across_angles(n_samples=2000, seed=42)

    print("All tasks completed successfully!")


if __name__ == "__main__":
    main()
