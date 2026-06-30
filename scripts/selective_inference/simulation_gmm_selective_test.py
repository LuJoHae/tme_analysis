import sys
import os
import time
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.datasets import make_blobs
from sklearn.metrics import adjusted_rand_score
from sklearn.decomposition import PCA
import logging
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from gmm_mcmc_inference import GMMSelectiveInference

def plot_mcmc_trace(accepted_z, phi_obs, out_path, title):
    if len(accepted_z) == 0:
        return
    plt.figure(figsize=(10, 5))
    plt.plot(accepted_z, alpha=0.5, color='blue', label='MCMC Trace')
    plt.axhline(phi_obs, color='red', linestyle='--', label=f'Observed $\phi$ ({phi_obs:.2f})')
    plt.title(title)
    plt.xlabel('MCMC Step (post burn-in)')
    plt.ylabel('z value')
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

def plot_pca_comparison(X, true_labels, pred_labels, out_path, title):
    # Run PCA to 2 components for plotting
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5))
    
    # Raw Data
    sns.scatterplot(x=X_pca[:,0], y=X_pca[:,1], color='gray', ax=ax1, s=40, alpha=0.5)
    ax1.set_title("Raw Data")
    ax1.set_xlabel("PC1")
    ax1.set_ylabel("PC2")
    
    # True Labels
    sns.scatterplot(x=X_pca[:,0], y=X_pca[:,1], hue=true_labels, palette="tab10", ax=ax2, s=40, alpha=0.8, legend=False)
    ax2.set_title("True Labels (Ground Truth)")
    ax2.set_xlabel("PC1")
    ax2.set_ylabel("PC2")
    
    # Inferred Labels
    sns.scatterplot(x=X_pca[:,0], y=X_pca[:,1], hue=pred_labels, palette="Set1", ax=ax3, s=40, alpha=0.8, legend=False)
    ax3.set_title("Inferred Labels (GMM)")
    ax3.set_xlabel("PC1")
    ax3.set_ylabel("PC2")
    
    plt.suptitle(title, fontsize=14)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

def run_scenario(scenario_name, n_iterations, n_samples, n_features, sigma, centers, out_dir, start_time, time_limit):
    results = {'naive': [], 'selective': [], 'acceptance_rates': []}
    example_plotted = False
    
    logging.info(f"Running Scenario: {scenario_name} (sigma={sigma})")
    
    for i in tqdm(range(n_iterations)):
        if time.time() - start_time > time_limit:
            logging.warning("Time limit exceeded! Stopping simulation early.")
            break
            
        if centers == 1:
            # Null Hypothesis: 1 Blob
            X = np.random.normal(loc=0, scale=sigma, size=(n_samples, n_features))
            true_labels = np.zeros(n_samples, dtype=int)
        else:
            # Alternative Hypothesis: 2 Blobs
            X, true_labels = make_blobs(n_samples=n_samples, n_features=n_features, centers=centers, cluster_std=sigma, random_state=i)
            
        gmm_si = GMMSelectiveInference(n_components=2, random_state=i, max_mcmc_steps=1000)
        gmm_si.fit(X)
        pred_labels = gmm_si.labels_
        
        c1, c2 = 0, 1
        n1 = np.sum(pred_labels == c1)
        n2 = np.sum(pred_labels == c2)
        
        if n1 == 0 or n2 == 0:
            continue
            
        np_val, sp_val, plot_data = gmm_si.test_difference_in_means(X, c1, c2, sigma=sigma, burn_in=200)
        
        results['naive'].append(np_val)
        results['selective'].append(sp_val)
        if 'acceptance_rate' in plot_data:
            results['acceptance_rates'].append(plot_data['acceptance_rate'])
        
        if not example_plotted and 'accepted_z' in plot_data:
            ari = adjusted_rand_score(true_labels, pred_labels)
            prefix = scenario_name.lower().replace(" ", "_")
            plot_mcmc_trace(plot_data['accepted_z'], plot_data['phi_obs'],
                            os.path.join(out_dir, f"{prefix}_mcmc_trace.png"),
                            f"MCMC Trace: {scenario_name}")
            plot_pca_comparison(X, true_labels, pred_labels,
                                os.path.join(out_dir, f"{prefix}_pca.png"),
                                f"PCA: {scenario_name} (ARI={ari:.2f})")
            example_plotted = True
            
    for k in ['naive', 'selective']:
        results[k] = np.array(results[k])
        
    return results

def main():
    start_time = time.time()
    time_limit = 3600 # 1 hour
    
    out_dir = "output/gmm_mcmc_plots"
    os.makedirs(out_dir, exist_ok=True)
    
    n_iterations = 20
    n_samples = 100
    n_features = 10 # High dim so PCA makes sense
    
    report = {
        'parameters': {
            'n_iterations': n_iterations,
            'n_samples': n_samples,
            'n_features': n_features,
            'mcmc_steps': 1000,
            'burn_in': 200,
            'time_limit_seconds': time_limit
        },
        'results': {}
    }
    
    # Define Scenarios
    scenarios = [
        {"name": "Null_Hypothesis_LowVar", "centers": 1, "sigma": 1.0},
        {"name": "Null_Hypothesis_HighVar", "centers": 1, "sigma": 5.0},
        {"name": "Overlap_1_Tiny", "centers": 2, "sigma": 1.0},
        {"name": "Overlap_2_Low", "centers": 2, "sigma": 2.0},
        {"name": "Overlap_3_Moderate", "centers": 2, "sigma": 3.0},
        {"name": "Overlap_4_High", "centers": 2, "sigma": 4.0},
        {"name": "Overlap_5_Significant", "centers": 2, "sigma": 5.0},
        {"name": "Overlap_6_Severe", "centers": 2, "sigma": 6.0},
        {"name": "Overlap_7_Extreme", "centers": 2, "sigma": 8.0},
        {"name": "Overlap_8_NearNull", "centers": 2, "sigma": 10.0}
    ]
    
    for s in scenarios:
        res = run_scenario(s['name'], n_iterations, n_samples, n_features, s['sigma'], s['centers'], out_dir, start_time, time_limit)
        
        naive_metric = np.mean(res['naive'] <= 0.05) if len(res['naive']) > 0 else 0
        selective_metric = np.mean(res['selective'] <= 0.05) if len(res['selective']) > 0 else 0
        mean_acc = np.mean(res['acceptance_rates']) if len(res['acceptance_rates']) > 0 else 0
        
        report['results'][s['name']] = {
            'naive_metric': float(naive_metric),
            'selective_metric': float(selective_metric),
            'mean_acceptance_rate': float(mean_acc)
        }
        
        logging.info(f"--- {s['name']} Results ---")
        metric_name = "Type I Error" if s['centers'] == 1 else "Power"
        logging.info(f"  Naive {metric_name}:     {naive_metric:.3f}")
        logging.info(f"  Selective {metric_name}: {selective_metric:.3f}")
        logging.info(f"  Mean Acceptance Rate: {mean_acc:.3f}")

    report['total_script_seconds'] = time.time() - start_time
    
    with open(os.path.join(out_dir, "simulation_report.json"), "w") as f:
        json.dump(report, f, indent=4)
        
    logging.info("Simulation complete. Report and PCA plots saved.")

if __name__ == "__main__":
    main()
