from dataclasses import dataclass, field
import os
import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from scipy.stats import norm, uniform
from tqdm import tqdm
import logging
import pandas as pd
import altair as alt

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

@dataclass
class AnalyticalStats:
    naive_p: float
    analytical_p: float
    V_min: float
    V_max: float
    phi_obs: float
    std_phi: float

@dataclass
class MCMCStats:
    mcmc_p: float
    acceptance_rate: float
    accepted_z: np.ndarray = field(repr=False)


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
        
    return AnalyticalStats(
        naive_p=naive_p,
        analytical_p=analytical_p,
        V_min=V_min,
        V_max=V_max,
        phi_obs=phi_obs,
        std_phi=std_phi
    )

def run_mcmc_metropolis_hastings(X, labels, c1, c2, sigma=1.0, max_mcmc_steps=2000, 
                                 burn_in=500, random_state=42, model_type="kmeans", n_init=None):
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
        
        # Fit clustering model on perturbed data
        if model_type == "kmeans":
            n_init_km = n_init if n_init is not None else 5
            kmeans_prop = KMeans(n_clusters=2, random_state=random_state, n_init=n_init_km)
            kmeans_prop.fit(X_prop)
            labels_prop = kmeans_prop.labels_
        elif model_type in ["gmm_spherical", "gmm_diag", "gmm_full", "gmm_tied"]:
            n_init_gmm = n_init if n_init is not None else 5
            gmm_prop = GaussianMixture(n_components=2, covariance_type=model_type.split("_")[1], random_state=random_state, n_init=n_init_gmm)
            gmm_prop.fit(X_prop)
            labels_prop = gmm_prop.predict(X_prop)
        else:
            raise ValueError(f"Unknown model_type: {model_type}")
        
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
        mcmc_p = min(1.0, 2 * min(np.mean(accepted_z >= phi_obs), np.mean(accepted_z <= phi_obs)))
        
    acceptance_rate = accept_count / max_mcmc_steps

    return MCMCStats(mcmc_p=mcmc_p, accepted_z=accepted_z, acceptance_rate=acceptance_rate)

def plot_generated_data(X, pred_labels, centers, args, ax):
    model_type = getattr(args, 'model_type', 'kmeans')
    model_name = "K-Means" if model_type == "kmeans" else "GMM"
    center_label = f"{model_name} Centers"
    # Plot generated data
    sns.scatterplot(x=X[:, 0], y=X[:, 1], hue=pred_labels, palette="tab10", alpha=0.8, s=60, edgecolor='k', ax=ax)
    ax.scatter(centers[:, 0], centers[:, 1], color='black', marker='X', s=250, label=center_label, zorder=10)
    ax.set_title(f"{model_name} Clustering (True Mean Difference = {args.mean_diff})")
    ax.set_xlabel("$X_1$")
    ax.set_ylabel("$X_2$")
    ax.legend()
        
    return ax

def plot_mcmc_diagnostics(accepted_z, ax):
    if len(accepted_z) == 0:
        logging.warning("No accepted MCMC samples to plot for diagnostics.")
        return (None, None)

    # Plot MCMC diagnostics vs Analytical Truncated Normal
    ax1, ax2 = ax[0], ax[1]
        
    # Trace Plot
    ax1.plot(accepted_z, color='teal', alpha=0.7)
    ax1.set_title("MCMC Trace Plot")
    ax1.set_xlabel("Steps (post burn-in)")
    ax1.set_ylabel("Difference in means value ($z$)")
    ax1.legend(loc="upper right")

    # Density plot comparing empirical MCMC and analytical truncated normal
    sns.histplot(accepted_z, kde=False, color='teal', ax=ax2, stat="density", alpha=0.5, label="Empirical MCMC")
    ax2.set_title("Distribution of Test Statistic")
    ax2.set_xlabel("Difference in means value ($z$)")
    ax2.set_ylabel("Density")
    ax2.legend()
        
    return (ax1, ax2)

def plot_generated_data_altair(X, pred_labels, centers, args):
    df_points = pd.DataFrame({
        'X_1': X[:, 0],
        'X_2': X[:, 1],
        'Cluster': pred_labels.astype(str)
    })
    
    df_centers = pd.DataFrame({
        'X_1': centers[:, 0],
        'X_2': centers[:, 1]
    })
    
    scatter = alt.Chart(df_points).mark_circle(size=60, opacity=0.8, stroke='black', strokeWidth=0.5).encode(
        x=alt.X('X_1:Q', title='X1'),
        y=alt.Y('X_2:Q', title='X2'),
        color=alt.Color('Cluster:N', scale=alt.Scale(scheme='tableau10'), title='Cluster')
    )
    
    centers_plot = alt.Chart(df_centers).mark_point(
        shape='cross', 
        size=250, 
        color='black', 
        strokeWidth=3,
        opacity=1.0
    ).encode(
        x='X_1:Q',
        y='X_2:Q'
    )
    
    model_type = getattr(args, 'model_type', 'kmeans')
    model_name = model_type
    chart = (scatter + centers_plot).properties(
        title=f"{model_name} Clustering (True Mean Difference = {args.mean_diff})",
        width=400,
        height=300
    ).interactive()
    
    return chart

def plot_mcmc_diagnostics_altair(accepted_z, phi_obs):
    if len(accepted_z) == 0:
        logging.warning("No accepted MCMC samples to plot for diagnostics.")
        return None
        
    df_trace = pd.DataFrame({
        'Steps': np.arange(len(accepted_z)),
        'z': accepted_z
    })
    
    trace_chart = alt.Chart(df_trace).mark_line(color='teal', opacity=0.7).encode(
        x=alt.X('Steps:Q', title='Steps (post burn-in)'),
        y=alt.Y('z:Q', title='Difference in means value (z)', scale=alt.Scale(zero=False))
    ).properties(
        title='MCMC Trace Plot',
        width=400,
        height=300
    ).interactive(name='trace_select')
    
    hline = alt.Chart(pd.DataFrame({'y': [phi_obs]})).mark_rule(
        color='red',
        strokeDash=[4, 4]
    ).encode(
        y='y:Q'
    )

    trace_chart = trace_chart + hline

    df_density = pd.DataFrame({
        'z': accepted_z
    })
    
    density_chart = alt.Chart(df_density).mark_bar(opacity=0.5, color='teal').encode(
        x=alt.X('z:Q', bin=alt.Bin(maxbins=30), title='Difference in means value (z)'),
        y=alt.Y('count():Q', title='Frequency')
    ).properties(
        title='Distribution of Test Statistic',
        width=400,
        height=300
    ).interactive(name='density_select')
    
    return trace_chart, density_chart

def run_single_run_diagnostics(args):
    logging.info("Generating dataset for single run diagnostics...")
    X, true_labels = generate_mixture_data(n_samples_comp1=args.n_samples, n_samples_comp2=args.n_samples, 
                                           mean_diff=args.mean_diff, sigma=1.0, random_state=args.seed)
    
    model_type = getattr(args, 'model_type', 'kmeans')
    if model_type == "kmeans":
        kmeans = KMeans(n_clusters=2, random_state=args.seed, n_init=args.n_init)
        kmeans.fit(X)
        pred_labels = kmeans.labels_
        centers = kmeans.cluster_centers_
    elif model_type in ["gmm_spherical", "gmm_diag", "gmm_full", "gmm_tied"]:
        gmm = GaussianMixture(n_components=2, covariance_type=model_type.split("_")[1], random_state=args.seed, n_init=args.n_init)
        gmm.fit(X)
        pred_labels = gmm.predict(X)
        centers = gmm.means_
    else:
        raise ValueError(f"Unknown model_type: {model_type}")
    
   
    logging.info("Running MCMC Metropolis-Hastings sampler for single run...")
    mcmc_stats = run_mcmc_metropolis_hastings(
        X, pred_labels, c1=0, c2=1, sigma=1.0, max_mcmc_steps=args.mcmc_steps, 
        burn_in=args.burn_in, random_state=args.seed, model_type=args.model_type,
        n_init=args.n_init
    )
    logging.info(mcmc_stats)
    
    analytical_stats = compute_naive_and_analytical_p(X, pred_labels, centers, c1=0, c2=1, sigma=1.0)
    logging.info(analytical_stats)

    # Use Altair function
    chart_data = plot_generated_data_altair(X, pred_labels, centers, args)
    trace_chart, density_chart = plot_mcmc_diagnostics_altair(mcmc_stats.accepted_z, analytical_stats.phi_obs)

    return (chart_data, trace_chart, density_chart), (mcmc_stats, analytical_stats)

def run_ecdf_simulation(args, out_dir, ax=None):
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
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 7))
        fig_created = True
    # else:
    #     fig = ax.get_figure()
    #     fig_created = False

    ax.plot([0, 1], [0, 1], color='gray', linestyle='--', label='Uniform (Ideal Calibration)')
    
    # # Plot Naive
    # x_naive = np.sort(naive_null_pvals)
    # y_naive = np.arange(1, len(x_naive) + 1) / len(x_naive)
    # ax.step(x_naive, y_naive, label=f'Naive p-values (mean={np.mean(x_naive):.3f}, Type I={np.mean(x_naive <= 0.05):.3f})', 
    #          color='crimson', linewidth=2.5)
    
    # # Plot Analytical
    # x_ana = np.sort(analytical_null_pvals)
    # y_ana = np.arange(1, len(x_ana) + 1) / len(x_ana)
    # ax.step(x_ana, y_ana, label=f'Analytical Selective p-values (mean={np.mean(x_ana):.3f}, Type I={np.mean(x_ana <= 0.05):.3f})', 
    #          color='darkorange', linewidth=2.5)
    
    # Plot MCMC
    x_mcmc = np.sort(mcmc_null_pvals)
    y_mcmc = np.arange(1, len(x_mcmc) + 1) / len(x_mcmc)
    ax.step(x_mcmc, y_mcmc, label=f'MCMC Selective p-values (mean={np.mean(x_mcmc):.3f}, Type I={np.mean(x_mcmc <= 0.05):.3f})', 
             color='teal', linewidth=2.5)
    
    ax.set_title("Empirical Cumulative Distribution Function (ECDF) under the Null")
    ax.set_xlabel("p-value")
    ax.set_ylabel("Fraction of trials <= p-value")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc='lower right')
    
    # if fig_created:
    #     fig.tight_layout()
    #     ecdf_plot_path = os.path.join(out_dir, "ecdf_comparison.png")
    #     fig.savefig(ecdf_plot_path, dpi=150)
    #     plt.close(fig)
    #     logging.info(f"Saved ECDF calibration plot to {ecdf_plot_path}")

    return ax

def qq_plot(data: np.ndarray, dist=uniform, params=()):
    """
    Creates an interactive Q-Q plot using Altair and NumPy.
    """
    # 1. Sort data and compute empirical quantiles
    sorted_data = np.sort(data)
    n = len(sorted_data)
    probs = (np.arange(n) + 0.5) / n
    
    # 2. Compute theoretical quantiles
    theoretical_quantiles = dist.ppf(probs, *params)
    
    # 3. Prepare data for Altair
    source = alt.Data(values=[
        {"sample": float(s), "theoretical": float(t)} 
        for s, t in zip(sorted_data, theoretical_quantiles)
    ])
    
    # 4. Create base plot
    points = alt.Chart(source).mark_point().encode(
        x=alt.X('theoretical:Q', title='Theoretical Quantiles'),
        y=alt.Y('sample:Q', title='Sample Quantiles')
    )
    
    # 5. Add reference line (y = x) without Pandas
    # Casting to float is required for Altair's JSON serialization
    min_val = float(min(sorted_data[0], theoretical_quantiles[0]))
    max_val = float(max(sorted_data[-1], theoretical_quantiles[-1]))
    
    line_data = alt.Data(values=[
        {'x': min_val, 'y': min_val},
        {'x': max_val, 'y': max_val}
    ])
    
    line = alt.Chart(line_data).mark_line(
        color='red', strokeDash=[5, 5]
    ).encode(
        x='x:Q', 
        y='y:Q'
    )
    
    return points + line