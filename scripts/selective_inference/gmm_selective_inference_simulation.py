import os
import time
import json
import csv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
from matplotlib.patches import Ellipse
from tqdm import tqdm
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

def draw_ellipse(position, covariance, ax=None, **kwargs):
    """Draw a 2D covariance ellipse."""
    ax = ax or plt.gca()
    
    # Calculate eigenvalues and eigenvectors
    v, w = np.linalg.eigh(covariance)
    u = w[0] / np.linalg.norm(w[0])
    angle = np.arctan2(u[1], u[0])
    angle = 180 * angle / np.pi  # Convert to degrees
    
    # 2 standard deviations ellipse size
    v = 2.0 * np.sqrt(2.0) * np.sqrt(v)
    
    color = kwargs.pop('color', 'blue')
    alpha = kwargs.pop('alpha', 0.2)
    
    ell = Ellipse(position, v[0], v[1], angle=angle, facecolor='none', edgecolor=color, alpha=alpha, linewidth=2, **kwargs)
    ax.add_patch(ell)

def generate_mixture_data(n_samples_comp1=75, n_samples_comp2=75, n_features=2, random_state=42):
    """Generate data from a mixture of two Gaussians (centered at 0, with standard deviations 1.0 and 2.0)."""
    np.random.seed(random_state)
    # Component 1: Centered at 0 with standard deviation 1
    X1 = np.random.normal(loc=0.0, scale=1.0, size=(n_samples_comp1, n_features))
    # Component 2: Centered at 0 with standard deviation 2
    X2 = np.random.normal(loc=0.0, scale=2.0, size=(n_samples_comp2, n_features))
    X = np.vstack([X1, X2])
    return X

def test_difference_in_means_gmm(X, labels, c1, c2, n_components=2, sigma=1.0, max_mcmc_steps=1000, burn_in=200, random_state=42):
    n = X.shape[0]
    d = X.shape[1]
    
    c1_mask = labels == c1
    c2_mask = labels == c2
    
    n1 = np.sum(c1_mask)
    n2 = np.sum(c2_mask)
    
    if n1 == 0 or n2 == 0:
        return 1.0, 1.0, {}
        
    mean1 = np.mean(X[c1_mask], axis=0)
    mean2 = np.mean(X[c2_mask], axis=0)
    
    diff = mean1 - mean2
    diff_norm = np.linalg.norm(diff)
    
    if diff_norm == 0:
        return 1.0, 1.0, {}
        
    dir_vec = diff / diff_norm
    
    V = np.zeros_like(X)
    V[c1_mask] = dir_vec / n1
    V[c2_mask] = -dir_vec / n2
    
    phi_obs = np.sum(V * X)
    v_norm_sq = np.sum(V**2)
    std_phi = sigma * np.sqrt(v_norm_sq)
    
    # Naive p-value
    z_val_naive = phi_obs / std_phi
    naive_p = 2 * (1 - norm.cdf(abs(z_val_naive)))
    
    # MCMC Hit-and-Run on the 1D line Z
    accepted_z = []
    z_current = phi_obs
    tau = std_phi * 0.5 
    accept_count = 0
    
    for step in range(max_mcmc_steps):
        z_prop = np.random.normal(loc=z_current, scale=tau)
        X_prop = X + (z_prop - phi_obs) * V
        
        # Fit GMM with tied covariance
        gmm_prop = GaussianMixture(n_components=n_components, covariance_type='tied', random_state=random_state)
        gmm_prop.fit(X_prop)
        labels_prop = gmm_prop.predict(X_prop)
        
        # Check if assignments match up to label permutation
        if np.array_equal(labels_prop, labels) or np.array_equal(labels_prop, 1 - labels):
            log_ratio = -0.5 * ((z_prop/std_phi)**2 - (z_current/std_phi)**2)
            if np.log(np.random.uniform()) < log_ratio:
                z_current = z_prop
                accept_count += 1
                
        if step >= burn_in:
            accepted_z.append(z_current)
            
    accepted_z = np.array(accepted_z)
    
    if len(accepted_z) == 0:
        selective_p = 1.0
    else:
        selective_p = 2 * min(np.mean(accepted_z >= phi_obs), np.mean(accepted_z <= phi_obs))
        
    plot_data = {
        'accepted_z': accepted_z,
        'phi_obs': phi_obs,
        'std_phi': std_phi,
        'acceptance_rate': accept_count / max_mcmc_steps
    }
    
    return naive_p, selective_p, plot_data

def main():
    # Set style for publication-quality plots
    sns.set_theme(style="whitegrid")
    plt.rcParams.update({
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'figure.titlesize': 18
    })

    out_dir = "output/gmm_selective_inference_sim_plots"
    os.makedirs(out_dir, exist_ok=True)
    
    logging.info(f"Created output directory: {out_dir}")

    # Set parameters
    n_samples_comp1 = 75
    n_samples_comp2 = 75
    n_features = 2
    
    # 1. Generate null data (mixture of two Gaussians)
    logging.info("Generating mixture of two Gaussians (representing a null with unequal variances)...")
    X = generate_mixture_data(n_samples_comp1=n_samples_comp1, n_samples_comp2=n_samples_comp2, n_features=n_features, random_state=42)
    
    plt.figure(figsize=(8, 6))
    plt.scatter(X[:n_samples_comp1, 0], X[:n_samples_comp1, 1], color='coral', alpha=0.7, label='Comp 1 (sigma=1)')
    plt.scatter(X[n_samples_comp1:, 0], X[n_samples_comp1:, 1], color='royalblue', alpha=0.7, label='Comp 2 (sigma=2)')
    plt.title("Generated Data (Mixture of Two Gaussians - Null)")
    plt.xlabel("$X_1$")
    plt.ylabel("$X_2$")
    plt.legend()
    plt.tight_layout()
    plot_path = os.path.join(out_dir, "null_data.png")
    plt.savefig(plot_path, dpi=150)
    plt.close()
    logging.info(f"Saved null data plot to {plot_path}")

    # 2. Fit GMM (tied covariance) and run single test
    logging.info("Fitting GMM (k=2, tied covariance) to null data...")
    gmm = GaussianMixture(n_components=2, covariance_type='tied', random_state=42)
    gmm.fit(X)
    pred_labels = gmm.predict(X)
    
    # Plot data with indication of the GMM components and covariance
    logging.info("Plotting data with GMM cluster assignments and covariance ellipses...")
    plt.figure(figsize=(9, 7))
    sns.scatterplot(x=X[:, 0], y=X[:, 1], hue=pred_labels, palette="tab10", alpha=0.7, s=50, edgecolor='k')
    plt.scatter(gmm.means_[:, 0], gmm.means_[:, 1], color='black', marker='X', s=200, label='GMM Centers', zorder=10)
    
    # Covariances for tied covariance GMM has shape (n_features, n_features) -> (2, 2)
    cov = gmm.covariances_
    for k in range(gmm.n_components):
        mean = gmm.means_[k]
        # Draw 1 std dev ellipse
        draw_ellipse(mean, cov, color='darkblue', alpha=0.4, label='1 Std Dev Contour' if k == 0 else "")
        # Draw 2 std dev ellipse (scaled by 4 since variance scales by factor^2)
        draw_ellipse(mean, cov * 4, color='darkblue', alpha=0.15, label='2 Std Dev Contour' if k == 0 else "")
        
    plt.title("GMM Clustering Fit on Null Data (Tied Covariance)")
    plt.xlabel("$X_1$")
    plt.ylabel("$X_2$")
    plt.legend(loc='upper right')
    plt.tight_layout()
    gmm_fit_path = os.path.join(out_dir, "gmm_fitted_clusters.png")
    plt.savefig(gmm_fit_path, dpi=150)
    plt.close()
    logging.info(f"Saved GMM fitted clusters plot to {gmm_fit_path}")

    # Dynamically estimate sigma
    sigma_est = np.std(X)
    logging.info(f"Estimated sigma for single run: {sigma_est:.4f}")

    logging.info("Running single GMM MCMC Selective Inference test (10k steps)...")
    naive_p, selective_p, plot_data = test_difference_in_means_gmm(
        X, pred_labels, c1=0, c2=1, n_components=2, sigma=sigma_est, max_mcmc_steps=10000, burn_in=2000, random_state=42
    )
    
    logging.info(f"Single test: Naive p = {naive_p:.6f}, Selective p = {selective_p:.6f}")
    logging.info(f"MCMC Acceptance Rate: {plot_data['acceptance_rate']:.2%}")

    # Export Single Run Results CSV
    csv_path_single = os.path.join(out_dir, "single_run_results.csv")
    with open(csv_path_single, mode="w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["naive_p", "selective_p", "acceptance_rate", "phi_obs", "std_phi"])
        writer.writerow([naive_p, selective_p, plot_data["acceptance_rate"], plot_data["phi_obs"], plot_data["std_phi"]])
    logging.info(f"Saved single run results to {csv_path_single}")

    # Export MCMC Samples CSV
    csv_path_samples = os.path.join(out_dir, "mcmc_samples.csv")
    with open(csv_path_samples, mode="w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample_index", "z_value"])
        for idx, val in enumerate(plot_data["accepted_z"]):
            writer.writerow([idx, val])
    logging.info(f"Saved MCMC samples to {csv_path_samples}")

    # 3. Generate Diagnostics Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))

    # Trace Plot
    ax1.plot(plot_data['accepted_z'], color='teal', alpha=0.6)
    ax1.axhline(plot_data['phi_obs'], color='red', linestyle='--', linewidth=2, label=f'Observed $\\phi$ ({plot_data["phi_obs"]:.3f})')
    ax1.set_title("MCMC Trace Plot")
    ax1.set_xlabel("Steps (post burn-in)")
    ax1.set_ylabel("z value")
    ax1.legend()

    # Density Plot
    sns.histplot(plot_data['accepted_z'], kde=True, color='teal', ax=ax2, stat="density")
    ax2.axvline(plot_data['phi_obs'], color='red', linestyle='--', linewidth=2, label=r'Observed $\phi$')
    z_range = np.linspace(np.min(plot_data['accepted_z']) - 1, np.max(plot_data['accepted_z']) + 1, 200)
    ax2.plot(z_range, norm.pdf(z_range, 0, plot_data['std_phi']), color='darkorange', linestyle=':', linewidth=2, label='Unconditioned Null')
    ax2.set_title("Empirical Null Distribution")
    ax2.set_xlabel("z value")
    ax2.legend()

    plt.tight_layout()
    diag_path = os.path.join(out_dir, "mcmc_diagnostics.png")
    plt.savefig(diag_path, dpi=150)
    plt.close()
    logging.info(f"Saved diagnostics plot to {diag_path}")

    # 4. Monte Carlo Simulation
    n_simulations = 20
    naive_p_vals = []
    selective_p_vals = []
    
    logging.info(f"Running Monte Carlo simulation of {n_simulations} runs under the null (5k steps)...")
    for run in tqdm(range(n_simulations), desc="Running Simulations"):
        X_sim = generate_mixture_data(n_samples_comp1=n_samples_comp1, n_samples_comp2=n_samples_comp2, n_features=n_features, random_state=run)
        
        # Fit GMM with tied covariance
        gmm_sim = GaussianMixture(n_components=2, covariance_type='tied', random_state=run)
        gmm_sim.fit(X_sim)
        labels_sim = gmm_sim.predict(X_sim)
        
        if len(np.unique(labels_sim)) < 2:
            continue
            
        # Dynamically estimate sigma
        sigma_sim_est = np.std(X_sim)
        
        np_val, sp_val, _ = test_difference_in_means_gmm(
            X_sim, labels_sim, c1=0, c2=1, n_components=2, sigma=sigma_sim_est,
            max_mcmc_steps=5000, burn_in=1000, random_state=run
        )
        
        naive_p_vals.append(np_val)
        selective_p_vals.append(sp_val)
        
    naive_p_vals = np.array(naive_p_vals)
    selective_p_vals = np.array(selective_p_vals)
    
    logging.info(f"Completed simulation. Naive p-value mean: {np.mean(naive_p_vals):.4f}, Selective p-value mean: {np.mean(selective_p_vals):.4f}")

    # Export Simulation Results CSV
    csv_path_sim = os.path.join(out_dir, "simulation_results.csv")
    with open(csv_path_sim, mode="w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["run", "naive_p", "selective_p"])
        for idx, (np_val, sp_val) in enumerate(zip(naive_p_vals, selective_p_vals)):
            writer.writerow([idx, np_val, sp_val])
    logging.info(f"Saved simulation results to {csv_path_sim}")

    # 5. Plot ECDF Comparison
    plt.figure(figsize=(10, 7))
    plt.plot([0, 1], [0, 1], color='gray', linestyle='--', label='Uniform (Ideal Calibration)')
    
    x_naive = np.sort(naive_p_vals)
    y_naive = np.arange(1, len(x_naive) + 1) / len(x_naive)
    plt.step(x_naive, y_naive, label=f'Naive p-values (mean={np.mean(x_naive):.3f})', color='crimson', linewidth=2.5)
    
    x_sel = np.sort(selective_p_vals)
    y_sel = np.arange(1, len(x_sel) + 1) / len(x_sel)
    plt.step(x_sel, y_sel, label=f'Selective p-values (mean={np.mean(x_sel):.3f})', color='darkgreen', linewidth=2.5)
    
    plt.title("Empirical Cumulative Distribution Function (ECDF) of P-values")
    plt.xlabel("p-value")
    plt.ylabel("Fraction of trials <= p-value")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.legend(loc='lower right')
    plt.tight_layout()
    ecdf_path = os.path.join(out_dir, "ecdf_comparison.png")
    plt.savefig(ecdf_path, dpi=150)
    plt.close()
    logging.info(f"Saved ECDF comparison plot to {ecdf_path}")
    
    print("\nSimulation, CSV exports, and plotting successfully completed!")
    print(f"Results saved in: {out_dir}/")

if __name__ == "__main__":
    main()
