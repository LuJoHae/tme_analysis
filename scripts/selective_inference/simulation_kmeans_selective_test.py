import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.datasets import make_blobs
from sklearn.metrics import adjusted_rand_score
from scipy.stats import uniform
import logging
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

# Add the scripts directory to path to import local modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from kmeans_inference import KMeansSelectiveInference

def plot_qq(p_values, out_path, title):
    """Plot Q-Q plot for p-values vs Uniform(0,1) distribution."""
    if len(p_values) == 0:
        return
    sorted_p = np.sort(p_values)
    n = len(p_values)
    expected_p = np.arange(1, n + 1) / n
    
    plt.figure(figsize=(6, 6))
    plt.scatter(expected_p, sorted_p, color='blue', alpha=0.7)
    plt.plot([0, 1], [0, 1], color='red', linestyle='--')
    plt.xlabel('Expected Uniform Quantiles')
    plt.ylabel('Observed p-value Quantiles')
    plt.title(title)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

def plot_histograms(naive_p, selective_p, out_path, title):
    """Plot histograms of naive and selective p-values."""
    if len(naive_p) == 0:
        return
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    sns.histplot(naive_p, bins=20, ax=ax1, color='red', stat='density')
    ax1.set_title("Naive p-values")
    ax1.set_xlabel("p-value")
    ax1.set_xlim(0, 1)
    
    sns.histplot(selective_p, bins=20, ax=ax2, color='blue', stat='density')
    ax2.set_title("Selective p-values")
    ax2.set_xlabel("p-value")
    ax2.set_xlim(0, 1)
    
    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

def plot_dataset_example(X, true_labels, pred_labels, out_path, title):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    sns.scatterplot(x=X[:,0], y=X[:,1], hue=true_labels, palette="tab10", ax=ax1, s=30, alpha=0.8)
    ax1.set_title("True Labels")
    
    sns.scatterplot(x=X[:,0], y=X[:,1], hue=pred_labels, palette="Set1", ax=ax2, s=30, alpha=0.8)
    ax2.set_title("K-Means Labels")
    
    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

def run_simulation(n_iterations=100, n_samples=100, n_features=2, sigma=1.0, alternative=False, out_dir="", prefix=""):
    results = {
        'correct': {'naive': [], 'selective': []},
        'incorrect': {'naive': [], 'selective': []},
        'all': {'naive': [], 'selective': []}
    }
    
    example_plotted = False
    incorrect_example_plotted = False
    
    for i in tqdm(range(n_iterations)):
        if not alternative:
            # Null Hypothesis: 1 Blob (No true clusters)
            X = np.random.normal(loc=0, scale=sigma, size=(n_samples, n_features))
            true_labels = np.zeros(n_samples, dtype=int)
        else:
            # Alternative Hypothesis: 2 Blobs (True clusters)
            X, true_labels = make_blobs(n_samples=n_samples, n_features=n_features, centers=2, cluster_std=sigma, random_state=i)
            
        kmeans_si = KMeansSelectiveInference(n_clusters=2, random_state=i)
        kmeans_si.fit(X)
        pred_labels = kmeans_si.labels_
        
        # We always test the difference between cluster 0 and cluster 1
        c1, c2 = 0, 1
        n1 = np.sum(pred_labels == c1)
        n2 = np.sum(pred_labels == c2)
        
        if n1 == 0 or n2 == 0:
            continue
            
        np_val, sp_val, _ = kmeans_si.test_difference_in_means(X, c1, c2, sigma=sigma)
        
        results['all']['naive'].append(np_val)
        results['all']['selective'].append(sp_val)
        
        ari = adjusted_rand_score(true_labels, pred_labels)
        is_correct = ari == 1.0 if alternative else False
        
        if alternative:
            if is_correct:
                results['correct']['naive'].append(np_val)
                results['correct']['selective'].append(sp_val)
                if not example_plotted:
                    plot_dataset_example(X, true_labels, pred_labels, 
                                         os.path.join(out_dir, f"{prefix}_example_correct.png"),
                                         "Example Dataset (Correct Clustering)")
                    example_plotted = True
            else:
                results['incorrect']['naive'].append(np_val)
                results['incorrect']['selective'].append(sp_val)
                if not incorrect_example_plotted:
                    plot_dataset_example(X, true_labels, pred_labels, 
                                         os.path.join(out_dir, f"{prefix}_example_incorrect.png"),
                                         f"Example Dataset (Incorrect Clustering, ARI={ari:.2f})")
                    incorrect_example_plotted = True
        else:
            if not example_plotted:
                plot_dataset_example(X, true_labels, pred_labels, 
                                     os.path.join(out_dir, f"{prefix}_example_null.png"),
                                     "Example Dataset (Null Hypothesis)")
                example_plotted = True
                
    for k in results:
        results[k]['naive'] = np.array(results[k]['naive'])
        results[k]['selective'] = np.array(results[k]['selective'])
        
    return results

def main():
    out_dir = "output/selective_inference_plots"
    os.makedirs(out_dir, exist_ok=True)
    
    n_iterations = 200
    
    # -----------------------------------------------------
    # Experiment 1: Null Hypothesis (1 Blob)
    # -----------------------------------------------------
    logging.info(f"Running Null Hypothesis Simulation ({n_iterations} iterations)...")
    res_null = run_simulation(n_iterations=n_iterations, n_samples=50, n_features=2, alternative=False, out_dir=out_dir, prefix="null")
    naive_null = res_null['all']['naive']
    selective_null = res_null['all']['selective']
    
    type1_naive = np.mean(naive_null <= 0.05) if len(naive_null) > 0 else 0
    type1_selective = np.mean(selective_null <= 0.05) if len(selective_null) > 0 else 0
    
    logging.info(f"Null Hypothesis Results:")
    logging.info(f"  Naive Type I Error Rate:     {type1_naive:.3f}")
    logging.info(f"  Selective Type I Error Rate: {type1_selective:.3f}")
    
    plot_histograms(naive_null, selective_null, 
                    os.path.join(out_dir, "simulation_null_histograms.png"),
                    "P-values under Null Hypothesis")
    plot_qq(selective_null, 
            os.path.join(out_dir, "simulation_null_qq.png"),
            "Q-Q Plot: Selective P-values vs Uniform(0,1)")

    # -----------------------------------------------------
    # Experiment 2: Alternative Hypothesis (2 Blobs)
    # -----------------------------------------------------
    logging.info(f"\nRunning Alternative Hypothesis Simulation ({n_iterations} iterations)...")
    # Increase cluster_std to create overlapping blobs so k-means sometimes gets it wrong
    res_alt = run_simulation(n_iterations=n_iterations, n_samples=50, n_features=2, sigma=2.5, alternative=True, out_dir=out_dir, prefix="alt")
    
    # All
    naive_alt = res_alt['all']['naive']
    selective_alt = res_alt['all']['selective']
    logging.info(f"Alternative Hypothesis Results (All):")
    logging.info(f"  Naive Power:     {np.mean(naive_alt <= 0.05):.3f}")
    logging.info(f"  Selective Power: {np.mean(selective_alt <= 0.05):.3f}")
    plot_histograms(naive_alt, selective_alt, os.path.join(out_dir, "simulation_alt_all_histograms.png"), "P-values under Alt Hypothesis (All)")

    # Correct
    naive_corr = res_alt['correct']['naive']
    selective_corr = res_alt['correct']['selective']
    logging.info(f"\nAlternative Hypothesis (Correct Clustering, N={len(naive_corr)}):")
    if len(naive_corr) > 0:
        logging.info(f"  Naive Power:     {np.mean(naive_corr <= 0.05):.3f}")
        logging.info(f"  Selective Power: {np.mean(selective_corr <= 0.05):.3f}")
        plot_histograms(naive_corr, selective_corr, os.path.join(out_dir, "simulation_alt_correct_histograms.png"), "P-values under Alt Hypothesis (Correct Clustering)")

    # Incorrect
    naive_inc = res_alt['incorrect']['naive']
    selective_inc = res_alt['incorrect']['selective']
    logging.info(f"\nAlternative Hypothesis (Incorrect Clustering, N={len(naive_inc)}):")
    if len(naive_inc) > 0:
        logging.info(f"  Naive Power:     {np.mean(naive_inc <= 0.05):.3f}")
        logging.info(f"  Selective Power: {np.mean(selective_inc <= 0.05):.3f}")
        plot_histograms(naive_inc, selective_inc, os.path.join(out_dir, "simulation_alt_incorrect_histograms.png"), "P-values under Alt Hypothesis (Incorrect Clustering)")

if __name__ == "__main__":
    main()
