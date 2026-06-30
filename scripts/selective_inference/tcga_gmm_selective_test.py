import sys
import os
import time
import json
import numpy as np
import pandas as pd
import datalair
import scanpy as sc
import logging
import importlib

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from gmm_mcmc_inference import GMMSelectiveInference
from ici_datasets.other_datasets import load_tcga

def main():
    start_time = time.time()
    time_limit = 3600 # 1 hour
    
    out_dir = "output/gmm_mcmc_plots"
    os.makedirs(out_dir, exist_ok=True)
    
    report = {
        'parameters': {
            'n_samples_subsample': 500, # Heavy subsampling
            'mcmc_steps': 1000,
            'burn_in': 200,
            'time_limit_seconds': time_limit
        },
        'runtimes': {}
    }
    
    t0 = time.time()
    logging.info("Initializing Lair and loading TCGA data...")
    lair = datalair.Lair("/storage/halu/lair")
    adata_tcga = load_tcga(lair)
    adata_tcga.var_names_make_unique()
    adata_tcga.X = adata_tcga.X.astype(np.float32)
    report['runtimes']['data_load_seconds'] = time.time() - t0
    
    t0 = time.time()
    logging.info(f"Subsampling to {report['parameters']['n_samples_subsample']} samples for MCMC computational feasibility...")
    np.random.seed(42)
    indices = np.random.choice(adata_tcga.shape[0], report['parameters']['n_samples_subsample'], replace=False)
    adata_tcga = adata_tcga[indices].copy()
    
    sc.pp.filter_genes(adata_tcga, min_cells=1)
    sc.pp.normalize_total(adata_tcga, target_sum=1e4)
    sc.pp.log1p(adata_tcga)
    sc.pp.highly_variable_genes(adata_tcga, n_top_genes=1000)
    adata_tcga = adata_tcga[:, adata_tcga.var.highly_variable].copy()
    sc.pp.scale(adata_tcga, max_value=10)
    sc.tl.pca(adata_tcga, n_comps=20)
    
    X = adata_tcga.obsm['X_pca']
    report['runtimes']['preprocessing_seconds'] = time.time() - t0
    
    t0 = time.time()
    logging.info("Running GMM Clustering (k=5)...")
    gmm_si = GMMSelectiveInference(n_components=5, random_state=42, max_mcmc_steps=report['parameters']['mcmc_steps'])
    gmm_si.fit(X)
    pred_labels = gmm_si.labels_
    
    unique_clusters, counts = np.unique(pred_labels, return_counts=True)
    sorted_clusters = unique_clusters[np.argsort(-counts)]
    c1, c2 = sorted_clusters[0], sorted_clusters[1]
    
    logging.info(f"Testing difference in means between largest cluster {c1} and {c2}...")
    sigma_est = np.std(X)
    
    naive_p, selective_p, plot_data = gmm_si.test_difference_in_means(X, c1, c2, sigma=sigma_est, burn_in=report['parameters']['burn_in'])
    report['runtimes']['gmm_mcmc_inference_seconds'] = time.time() - t0
    
    report['results'] = {
        'naive_p': float(naive_p),
        'selective_p': float(selective_p),
        'cluster_1_size': int(counts[np.where(unique_clusters == c1)][0]),
        'cluster_2_size': int(counts[np.where(unique_clusters == c2)][0]),
        'acceptance_rate': float(plot_data.get('acceptance_rate', 0.0))
    }
    
    logging.info(f"Naive p-value: {naive_p:.3e}")
    logging.info(f"Selective p-value: {selective_p:.3e}")
    logging.info(f"MCMC Acceptance Rate: {plot_data.get('acceptance_rate', 0.0):.3f}")
    
    report['runtimes']['total_script_seconds'] = time.time() - start_time
    
    with open(os.path.join(out_dir, "tcga_report.json"), "w") as f:
        json.dump(report, f, indent=4)
        
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # 1. PCA by GMM Clusters
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=X[:, 0], y=X[:, 1], hue=pred_labels, palette="Set1", s=40, alpha=0.9)
    plt.title(f"PCA of TCGA Data (Subsampled N={report['parameters']['n_samples_subsample']}) colored by GMM Clusters")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend(title="GMM Cluster", bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "tcga_gmm_pca_clusters.png"), dpi=150)
    plt.close()
    
    # 2. MCMC Trace Plot
    if 'accepted_z' in plot_data and len(plot_data['accepted_z']) > 0:
        plt.figure(figsize=(10, 5))
        plt.plot(plot_data['accepted_z'], alpha=0.5, color='purple', label='MCMC Trace')
        plt.axhline(plot_data['phi_obs'], color='red', linestyle='--', label=f'Observed $\phi$ ({plot_data["phi_obs"]:.2f})')
        plt.title(f"MCMC Trace: Testing Cluster {c1} vs {c2} (TCGA Data)")
        plt.xlabel("MCMC Step (post burn-in)")
        plt.ylabel("z value")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, "tcga_gmm_mcmc_trace.png"), dpi=150)
        plt.close()
        
    logging.info("TCGA GMM MCMC complete. Report and plots saved.")

if __name__ == "__main__":
    main()
