import os
import time
import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from scipy.stats import norm
from tqdm import tqdm
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def generate_mixture_data(n_samples_comp1=75, n_samples_comp2=75, n_features=2, 
                          mean_diff=1.5, sigma=1.0, random_state=42):
    """
    Generate data from a mixture of two Gaussians with same variance (sigma^2 * I)
    but different means along the first dimension.
    """
    np.random.seed(random_state)
    # Component 1 mean: [-mean_diff / 2, 0, ...]
    mean1 = np.zeros(n_features)
    mean1[0] = -mean_diff / 2.0
    
    # Component 2 mean: [mean_diff / 2, 0, ...]
    mean2 = np.zeros(n_features)
    mean2[0] = mean_diff / 2.0
    
    X1 = np.random.normal(loc=mean1, scale=sigma, size=(n_samples_comp1, n_features))
    X2 = np.random.normal(loc=mean2, scale=sigma, size=(n_samples_comp2, n_features))
    X = np.vstack([X1, X2])
    
    # True component labels
    true_labels = np.hstack([np.zeros(n_samples_comp1, dtype=int), np.ones(n_samples_comp2, dtype=int)])
    return X, true_labels

def compute_analytical_limits(X, labels, centers, c1, c2, sigma=1.0):
    """
    Compute the truncation limits [V_min, V_max] for the contrast testing
    difference in means between cluster c1 and c2 under the final partition constraint.
    """
    n_samples, n_features = X.shape
    n_clusters = len(centers)
    
    diff = centers[c1] - centers[c2]
    diff_norm = np.linalg.norm(diff)
    if diff_norm == 0:
        return -np.inf, np.inf, 0.0, 0.0
        
    dir_vec = diff / diff_norm
    
    n1 = np.sum(labels == c1)
    n2 = np.sum(labels == c2)
    
    V = np.zeros_like(X)
    if n1 > 0:
        V[labels == c1] = dir_vec / n1
    if n2 > 0:
        V[labels == c2] = -dir_vec / n2
        
    v_norm_sq = np.sum(V**2)
    phi_obs = np.sum(V * X)  # Observed difference in means along delta
    
    # Projection of X orthogonal to V
    X_perp = X - V * (phi_obs / v_norm_sq)
    
    V_min = -np.inf
    V_max = np.inf
    
    # centers as functions of z: mu_k(z) = mean(X_perp[C_k]) + z * mean(V[C_k])
    centers_perp = np.zeros((n_clusters, n_features))
    centers_v = np.zeros((n_clusters, n_features))
    
    for k in range(n_clusters):
        mask = labels == k
        if np.sum(mask) > 0:
            centers_perp[k] = np.mean(X_perp[mask], axis=0)
            centers_v[k] = np.mean(V[mask], axis=0)
            
    for i in range(n_samples):
        k = labels[i]
        x_i_perp = X_perp[i]
        x_i_v = V[i]
        
        for j in range(n_clusters):
            if k == j:
                continue
                
            # Condition: || x_i(z) - mu_k(z) ||^2 <= || x_i(z) - mu_j(z) ||^2
            # Which expands to c2 * z^2 + c1 * z + c0 <= 0
            a_perp = x_i_perp - centers_perp[k]
            a_v = x_i_v - centers_v[k]
            
            b_perp = x_i_perp - centers_perp[j]
            b_v = x_i_v - centers_v[j]
            
            c2 = np.sum(a_v**2) - np.sum(b_v**2)
            c1 = 2 * (np.sum(a_perp * a_v) - np.sum(b_perp * b_v))
            c0 = np.sum(a_perp**2) - np.sum(b_perp**2)
            
            if np.abs(c2) < 1e-12:
                # Linear case: c1 * z + c0 <= 0
                if c1 > 1e-12:
                    V_max = min(V_max, -c0 / c1)
                elif c1 < -1e-12:
                    V_min = max(V_min, -c0 / c1)
            else:
                discriminant = c1**2 - 4 * c2 * c0
                if discriminant >= 0:
                    root1 = (-c1 - np.sqrt(discriminant)) / (2 * c2)
                    root2 = (-c1 + np.sqrt(discriminant)) / (2 * c2)
                    r_min = min(root1, root2)
                    r_max = max(root1, root2)
                    
                    if c2 > 0:
                        # Parabola opens upwards: <= 0 between roots
                        V_min = max(V_min, r_min)
                        V_max = min(V_max, r_max)
                    else:
                        # Parabola opens downwards: <= 0 outside roots (union)
                        # Keep the interval containing the observed statistic
                        if phi_obs <= r_min:
                            V_max = min(V_max, r_min)
                        elif phi_obs >= r_max:
                            V_min = max(V_min, r_max)
                            
    return V_min, V_max, phi_obs, np.sqrt(v_norm_sq)

def compute_naive_and_analytical_p(X, labels, centers, c1, c2, sigma=1.0):
    """
    Calculate naive and correct analytical selective p-values.
    """
    V_min, V_max, phi_obs, v_norm = compute_analytical_limits(X, labels, centers, c1, c2, sigma)
    
    if v_norm == 0:
        return 1.0, 1.0, V_min, V_max, phi_obs, 0.0
        
    std_phi = sigma * v_norm
    
    # Naive p-value (assuming normally distributed difference without selection)
    z_obs = np.abs(phi_obs) / std_phi
    naive_p = 2 * (1 - norm.cdf(z_obs))
    
    # Correct Analytical Selective p-value (using truncated normal)
    alpha_min = V_min / std_phi
    alpha_max = V_max / std_phi
    z_val = phi_obs / std_phi
    
    denom = norm.cdf(alpha_max) - norm.cdf(alpha_min)
    if denom <= 1e-15:
        # Fallback if numerical range is extremely restricted
        analytical_p = 1.0
    else:
        cdf_z = (norm.cdf(z_val) - norm.cdf(alpha_min)) / denom
        analytical_p = 2 * min(cdf_z, 1.0 - cdf_z)
        
    return naive_p, analytical_p, V_min, V_max, phi_obs, std_phi

def run_mcmc_metropolis_hastings(X, labels, c1, c2, sigma=1.0, max_mcmc_steps=2000, 
                                 burn_in=500, random_state=42):
    """
    Estimate the selective p-value using a Metropolis-Hastings Hit-and-Run sampler.
    """
    n_samples, n_features = X.shape
    n1 = np.sum(labels == c1)
    n2 = np.sum(labels == c2)
    
    if n1 == 0 or n2 == 0:
        return 1.0, [], 0.0
        
    # Fit original centers
    centers = np.zeros((2, n_features))
    centers[0] = np.mean(X[labels == 0], axis=0)
    centers[1] = np.mean(X[labels == 1], axis=0)
    
    diff = centers[c1] - centers[c2]
    diff_norm = np.linalg.norm(diff)
    if diff_norm == 0:
        return 1.0, [], 0.0
        
    dir_vec = diff / diff_norm
    
    V = np.zeros_like(X)
    V[labels == c1] = dir_vec / n1
    V[labels == c2] = -dir_vec / n2
    
    v_norm_sq = np.sum(V**2)
    phi_obs = np.sum(V * X)
    std_phi = sigma * np.sqrt(v_norm_sq)
    
    # Metropolis-Hastings sampler
    if max_mcmc_steps <= burn_in:
        new_burn_in = max(0, max_mcmc_steps // 2)
        logging.warning(
            f"max_mcmc_steps ({max_mcmc_steps}) is less than or equal to burn_in ({burn_in}). "
            f"Adjusting burn_in to {new_burn_in} to collect samples."
        )
        burn_in = new_burn_in

    np.random.seed(random_state)
    accepted_z = []
    z_current = phi_obs
    tau = std_phi * 0.5  # Proposal step size
    accept_count = 0
    
    for step in range(max_mcmc_steps):
        z_prop = np.random.normal(loc=z_current, scale=tau)
        
        # Correctly scale the perturbation to keep projection equal to z_prop
        X_prop = X + (z_prop - phi_obs) * V / v_norm_sq
        
        # Fit K-Means on perturbed data
        kmeans_prop = KMeans(n_clusters=2, random_state=random_state, n_init=5)
        kmeans_prop.fit(X_prop)
        labels_prop = kmeans_prop.labels_
        
        # Check if partition matches (up to cluster label swap)
        if np.array_equal(labels_prop, labels) or np.array_equal(labels_prop, 1 - labels):
            log_ratio = -0.5 * ((z_prop / std_phi)**2 - (z_current / std_phi)**2)
            if np.log(np.random.uniform()) < log_ratio:
                z_current = z_prop
                accept_count += 1
                
        if step >= burn_in:
            accepted_z.append(z_current)
            
    accepted_z = np.array(accepted_z)
    
    if len(accepted_z) == 0:
        mcmc_p = 1.0
    else:
        mcmc_p = 2 * min(np.mean(accepted_z >= phi_obs), np.mean(accepted_z <= phi_obs))
        
    acceptance_rate = accept_count / max_mcmc_steps
    return mcmc_p, accepted_z, acceptance_rate

def plot_generated_data(X, pred_labels, centers, args, out_dir):
    # Plot generated data
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    sns.scatterplot(x=X[:, 0], y=X[:, 1], hue=pred_labels, palette="tab10", alpha=0.8, s=60, edgecolor='k', ax=ax)
    ax.scatter(centers[:, 0], centers[:, 1], color='black', marker='X', s=250, label='K-Means Centers', zorder=10)
    ax.set_title(f"K-Means Clustering (True Mean Difference = {args.mean_diff})")
    ax.set_xlabel("$X_1$")
    ax.set_ylabel("$X_2$")
    ax.legend()
    fig.tight_layout()
    data_plot_path = os.path.join(out_dir, "generated_data.png")
    fig.savefig(data_plot_path, dpi=150)
    logging.info(f"Saved generated data plot to {data_plot_path}")
    return fig, ax

def plot_mcmc_diagnostics(accepted_z, phi_obs, V_min, V_max, std_phi, out_dir):
    if len(accepted_z) == 0:
        logging.warning("No accepted MCMC samples to plot for diagnostics.")
        return None, (None, None)

    # Plot MCMC diagnostics vs Analytical Truncated Normal
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Trace Plot
    ax1.plot(accepted_z, color='teal', alpha=0.7)
    ax1.axhline(phi_obs, color='red', linestyle='--', linewidth=2, label=f'Observed $\\phi$ ({phi_obs:.3f})')
    if V_min > -np.inf and V_min > np.min(accepted_z):
        ax1.axhline(V_min, color='darkorange', linestyle=':', linewidth=2, label=f'Analytical $V_{{min}}$ ({V_min:.3f})')
    if V_max < np.inf and V_max < np.max(accepted_z):
        ax1.axhline(V_max, color='darkorange', linestyle=':', linewidth=2, label=f'Analytical $V_{{max}}$ ({V_max:.3f})')
    ax1.set_title("MCMC Trace Plot")
    ax1.set_xlabel("Steps (post burn-in)")
    ax1.set_ylabel("Difference in means value ($z$)")
    ax1.legend(loc="upper right")

    # Density plot comparing empirical MCMC and analytical truncated normal
    sns.histplot(accepted_z, kde=False, color='teal', ax=ax2, stat="density", alpha=0.5, label="Empirical MCMC")
    ax2.axvline(phi_obs, color='red', linestyle='--', linewidth=2, label='Observed $\\phi$')

    # Plot analytical truncated normal PDF
    if V_min > -np.inf or V_max < np.inf:
        z_min_plot = max(V_min if V_min > -np.inf else np.min(accepted_z) - 0.5, np.min(accepted_z) - 0.5)
        z_max_plot = min(V_max if V_max < np.inf else np.max(accepted_z) + 0.5, np.max(accepted_z) + 0.5)
        z_grid = np.linspace(z_min_plot, z_max_plot, 200)

        # Truncated normal density
        denom = norm.cdf(V_max / std_phi) - norm.cdf(V_min / std_phi)
        if denom > 0:
            density = norm.pdf(z_grid / std_phi) / (std_phi * denom)
            ax2.plot(z_grid, density, color='darkorange', linewidth=2.5, label='Analytical Truncated Normal')

    # Unconditioned normal PDF (for comparison)
    z_full_grid = np.linspace(np.min(accepted_z) - 1, np.max(accepted_z) + 1, 200)
    ax2.plot(z_full_grid, norm.pdf(z_full_grid / std_phi) / std_phi, color='gray', linestyle=':', linewidth=1.5, label='Unconditioned Null')

    ax2.set_title("Distribution of Test Statistic")
    ax2.set_xlabel("Difference in means value ($z$)")
    ax2.set_ylabel("Density")
    ax2.legend()

    plt.tight_layout()
    diag_plot_path = os.path.join(out_dir, "mcmc_diagnostics.png")
    plt.savefig(diag_plot_path, dpi=150)
    plt.close()
    logging.info(f"Saved MCMC diagnostics plot to {diag_plot_path}")
    return fig, (ax1, ax2)

def run_single_run_diagnostics(args, out_dir):
    # -----------------------------------------------------------------
    # Step 1: Single Run Diagnostics (Visualizing null or alternative)
    # -----------------------------------------------------------------
    logging.info("Generating dataset for single run diagnostics...")
    X, true_labels = generate_mixture_data(n_samples_comp1=75, n_samples_comp2=75, 
                                           mean_diff=args.mean_diff, sigma=1.0, random_state=args.seed)
    
    kmeans = KMeans(n_clusters=2, random_state=args.seed, n_init=10)
    kmeans.fit(X)
    pred_labels = kmeans.labels_
    centers = kmeans.cluster_centers_
    
    _, ax1 = plot_generated_data(X, pred_labels, centers, args, out_dir)

    
    # Run p-value calculations
    naive_p, analytical_p, V_min, V_max, phi_obs, std_phi = compute_naive_and_analytical_p(
        X, pred_labels, centers, c1=0, c2=1, sigma=1.0
    )
    
    logging.info(f"Analytical Single Run Results:")
    logging.info(f"  Observed Test Statistic (phi_obs): {phi_obs:.4f}")
    logging.info(f"  Standard Deviation (std_phi):       {std_phi:.4f}")
    logging.info(f"  Truncation Interval:                 [{V_min:.4f}, {V_max:.4f}]")
    logging.info(f"  Naive p-value:                       {naive_p:.6f}")
    logging.info(f"  Analytical Selective p-value:        {analytical_p:.6f}")
    
    logging.info("Running MCMC Metropolis-Hastings sampler for single run...")
    mcmc_p, accepted_z, acceptance_rate = run_mcmc_metropolis_hastings(
        X, pred_labels, c1=0, c2=1, sigma=1.0, max_mcmc_steps=args.mcmc_steps, 
        burn_in=args.burn_in, random_state=args.seed
    )
    
    logging.info(f"MCMC Single Run Results:")
    logging.info(f"  MCMC Selective p-value:              {mcmc_p:.6f}")
    logging.info(f"  MCMC Acceptance Rate:                {acceptance_rate:.2%}")
    
    _, (ax2, ax3) = plot_mcmc_diagnostics(accepted_z, phi_obs, V_min, V_max, std_phi, out_dir)

    
    # Write Single Run Results to CSV
    csv_single_path = os.path.join(out_dir, "single_run_results.csv")
    with open(csv_single_path, mode="w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        writer.writerow(["mean_diff_setting", args.mean_diff])
        writer.writerow(["phi_obs", phi_obs])
        writer.writerow(["std_phi", std_phi])
        writer.writerow(["V_min", V_min])
        writer.writerow(["V_max", V_max])
        writer.writerow(["naive_p", naive_p])
        writer.writerow(["analytical_p", analytical_p])
        writer.writerow(["mcmc_p", mcmc_p])
        writer.writerow(["mcmc_acceptance_rate", acceptance_rate])
    logging.info(f"Saved single run CSV results to {csv_single_path}")

    return (ax1, ax2, ax3)

def run_ecdf_simulation(args, out_dir):
    # -----------------------------------------------------------------
    # Step 2: Multi-Iteration ECDF Simulation under the Null (D = 0)
    # -----------------------------------------------------------------
    logging.info(f"Starting Multi-Iteration Calibration Simulation under the NULL (D=0, N={args.n_simulations} iterations)...")
    naive_null_pvals = []
    analytical_null_pvals = []
    mcmc_null_pvals = []
    
    for sim in tqdm(range(args.n_simulations), desc="Simulating Null"):
        # Generate Null data (mean difference = 0)
        X_sim, _ = generate_mixture_data(n_samples_comp1=75, n_samples_comp2=75, 
                                         mean_diff=args.mean_diff, sigma=1.0, random_state=sim)
        
        # Fit K-Means
        kmeans_sim = KMeans(n_clusters=2, random_state=sim, n_init=5)
        kmeans_sim.fit(X_sim)
        labels_sim = kmeans_sim.labels_
        centers_sim = kmeans_sim.cluster_centers_
        
        if len(np.unique(labels_sim)) < 2:
            continue
            
        # Analytical p-values
        n_p, a_p, _, _, _, _ = compute_naive_and_analytical_p(
            X_sim, labels_sim, centers_sim, c1=0, c2=1, sigma=1.0
        )
        
        # MCMC p-values
        m_p, _, _ = run_mcmc_metropolis_hastings(
            X_sim, labels_sim, c1=0, c2=1, sigma=1.0, max_mcmc_steps=1000, 
            burn_in=200, random_state=sim
        )
        
        naive_null_pvals.append(n_p)
        analytical_null_pvals.append(a_p)
        mcmc_null_pvals.append(m_p)
        
    naive_null_pvals = np.array(naive_null_pvals)
    analytical_null_pvals = np.array(analytical_null_pvals)
    mcmc_null_pvals = np.array(mcmc_null_pvals)
    
    # Save Null Simulation Results to CSV
    csv_sim_path = os.path.join(out_dir, "null_simulation_results.csv")
    with open(csv_sim_path, mode="w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["iteration", "naive_p", "analytical_p", "mcmc_p"])
        for idx, (n_p, a_p, m_p) in enumerate(zip(naive_null_pvals, analytical_null_pvals, mcmc_null_pvals)):
            writer.writerow([idx, n_p, a_p, m_p])
    logging.info(f"Saved null simulation CSV results to {csv_sim_path}")
    
    # Plot ECDF comparison
    plt.figure(figsize=(10, 7))
    plt.plot([0, 1], [0, 1], color='gray', linestyle='--', label='Uniform (Ideal Calibration)')
    
    # Plot Naive
    x_naive = np.sort(naive_null_pvals)
    y_naive = np.arange(1, len(x_naive) + 1) / len(x_naive)
    plt.step(x_naive, y_naive, label=f'Naive p-values (mean={np.mean(x_naive):.3f}, Type I={np.mean(x_naive <= 0.05):.3f})', 
             color='crimson', linewidth=2.5)
    
    # Plot Analytical
    x_ana = np.sort(analytical_null_pvals)
    y_ana = np.arange(1, len(x_ana) + 1) / len(x_ana)
    plt.step(x_ana, y_ana, label=f'Analytical Selective p-values (mean={np.mean(x_ana):.3f}, Type I={np.mean(x_ana <= 0.05):.3f})', 
             color='darkorange', linewidth=2.5)
    
    # Plot MCMC
    x_mcmc = np.sort(mcmc_null_pvals)
    y_mcmc = np.arange(1, len(x_mcmc) + 1) / len(x_mcmc)
    plt.step(x_mcmc, y_mcmc, label=f'MCMC Selective p-values (mean={np.mean(x_mcmc):.3f}, Type I={np.mean(x_mcmc <= 0.05):.3f})', 
             color='teal', linewidth=2.5)
    
    plt.title("Empirical Cumulative Distribution Function (ECDF) under the Null")
    plt.xlabel("p-value")
    plt.ylabel("Fraction of trials <= p-value")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.legend(loc='lower right')
    plt.tight_layout()
    ecdf_plot_path = os.path.join(out_dir, "ecdf_comparison.png")
    plt.savefig(ecdf_plot_path, dpi=150)
    plt.close()
    logging.info(f"Saved ECDF calibration plot to {ecdf_plot_path}")
