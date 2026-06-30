import sys
import os
import numpy as np
import pandas as pd
import datalair
import scanpy as sc
from collections import Counter
from sklearn.metrics import adjusted_rand_score
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

# Add the scripts directory to path to import local modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from kmeans_inference import KMeansSelectiveInference

from ici_datasets.other_datasets import load_tcga
import importlib

def main():
    tcga_bg = importlib.import_module("TCGA-background")
    
    logging.info("Initializing Lair...")
    lair = datalair.Lair("/storage/halu/lair")
    
    logging.info("Loading TCGA data...")
    adata_tcga = load_tcga(lair)
    adata_tcga.var_names_make_unique()
    adata_tcga.X = adata_tcga.X.astype(np.float32)
    
    logging.info("Subsampling TCGA data for computational feasibility (e.g., 2000 random samples)...")
    # Subsample data to make custom k-means fast
    np.random.seed(42)
    indices = np.random.choice(adata_tcga.shape[0], 2000, replace=False)
    adata_tcga = adata_tcga[indices].copy()
    
    sc.pp.filter_genes(adata_tcga, min_cells=1)
    sc.pp.normalize_total(adata_tcga, target_sum=1e4)
    sc.pp.log1p(adata_tcga)
    
    logging.info("Selecting highly variable genes...")
    sc.pp.highly_variable_genes(adata_tcga, n_top_genes=1000)
    adata_tcga = adata_tcga[:, adata_tcga.var.highly_variable].copy()
    
    logging.info("Standardizing data...")
    sc.pp.scale(adata_tcga, max_value=10)
    
    # We will use PCA to reduce to top 20 components to speed up clustering and constraint solving
    logging.info("Running PCA...")
    sc.tl.pca(adata_tcga, n_comps=20)
    
    X = adata_tcga.obsm['X_pca']
    
    # TCGA usually has project or cancer_type metadata
    project_col = None
    for col in ["project", "project_id", "cancer_type", "type", "dataset"]:
        if col in adata_tcga.obs.columns:
            project_col = col
            break
            
    if project_col:
        logging.info(f"Using column '{project_col}' for ground truth cancer types.")
        true_labels = adata_tcga.obs[project_col].values
    else:
        logging.warning("No cancer type column found. Here are the columns: " + str(adata_tcga.obs.columns))
        true_labels = np.zeros(X.shape[0])
    
    logging.info("Running Selective K-Means Inference (k=10)...")
    kmeans_si = KMeansSelectiveInference(n_clusters=10, random_state=42)
    kmeans_si.fit(X)
    
    pred_labels = kmeans_si.labels_
    
    if project_col:
        ari = adjusted_rand_score(true_labels, pred_labels)
        logging.info(f"Adjusted Rand Index against true cancer types: {ari:.3f}")
        
    # Analyze the clusters
    unique_clusters, counts = np.unique(pred_labels, return_counts=True)
    logging.info(f"Cluster sizes: {dict(zip(unique_clusters, counts))}")
    
    # Find two largest clusters to test
    sorted_clusters = unique_clusters[np.argsort(-counts)]
    c1, c2 = sorted_clusters[0], sorted_clusters[1]
    
    logging.info(f"Testing difference in means between cluster {c1} (size {np.sum(pred_labels==c1)}) and cluster {c2} (size {np.sum(pred_labels==c2)})...")
    
    # Estimate noise variance (sigma) using standard sample variance approximation
    # For PCA space, data is already standardized and pca transformed, variance of noise is approx 1
    # We'll use a conservative estimate
    sigma_est = np.std(X)
    
    naive_p, selective_p, plot_data = kmeans_si.test_difference_in_means(X, c1, c2, sigma=sigma_est)
    
    logging.info(f"Naive p-value:     {naive_p:.3e}")
    logging.info(f"Selective p-value: {selective_p:.3e}")
    
    # Let's see what cancer types are in these two clusters
    if project_col:
        logging.info(f"--- Cancer types in Cluster {c1} ---")
        counts_c1 = Counter(true_labels[pred_labels == c1])
        for k, v in counts_c1.most_common(5):
            logging.info(f"  {k}: {v}")
            
        logging.info(f"--- Cancer types in Cluster {c2} ---")
        counts_c2 = Counter(true_labels[pred_labels == c2])
        for k, v in counts_c2.most_common(5):
            logging.info(f"  {k}: {v}")

    # ==========================================
    # PLOTTING
    # ==========================================
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy.stats import norm
    
    out_dir = "output/selective_inference_plots"
    os.makedirs(out_dir, exist_ok=True)
    
    # 1. PCA by True Cancer Type
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=X[:, 0], y=X[:, 1], hue=true_labels, palette="tab20", s=15, alpha=0.8)
    plt.title("PCA of TCGA Data colored by True Cancer Type")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize='small', ncol=2)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "pca_true_cancer_types.png"), dpi=150)
    plt.close()
    
    # 2. PCA by K-Means Clusters
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=X[:, 0], y=X[:, 1], hue=pred_labels, palette="Set1", s=15, alpha=0.8)
    plt.title("PCA of TCGA Data colored by K-Means Clusters")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "pca_kmeans_clusters.png"), dpi=150)
    plt.close()
    
    # 3. Cluster Composition
    # Create a dataframe for counting
    df_counts = pd.DataFrame({"Cluster": pred_labels, "CancerType": true_labels})
    ct = pd.crosstab(df_counts["Cluster"], df_counts["CancerType"], normalize="index") * 100
    
    plt.figure(figsize=(12, 8))
    sns.heatmap(ct, cmap="Blues", annot=False)
    plt.title("Cluster Composition: Percentage of Cancer Types per Cluster")
    plt.ylabel("K-Means Cluster")
    plt.xlabel("True Cancer Type")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "cluster_composition.png"), dpi=150)
    plt.close()
    
    # 4. Selective Inference Truncation Plot
    phi_obs = plot_data['phi_obs']
    std_phi = plot_data['std_phi']
    V_min = plot_data['V_min']
    V_max = plot_data['V_max']
    
    plt.figure(figsize=(10, 6))
    x_vals = np.linspace(-4 * std_phi, 4 * std_phi, 1000)
    pdf_vals = norm.pdf(x_vals, loc=0, scale=std_phi)
    
    plt.plot(x_vals, pdf_vals, label=r"Null Distribution N(0, $\sigma^2$)", color="black")
    
    # Highlight the truncated region
    trunc_mask = (x_vals >= V_min) & (x_vals <= V_max)
    plt.fill_between(x_vals[trunc_mask], 0, pdf_vals[trunc_mask], alpha=0.3, color="blue", label=r"Conditioning Region $[\mathcal{V}^-, \mathcal{V}^+]$")
    
    plt.axvline(phi_obs, color="red", linestyle="--", label=r"Observed $\phi$: " + f"{phi_obs:.2f}")
    plt.axvline(V_min, color="blue", linestyle=":", label=r"$\mathcal{V}^-$: " + f"{V_min:.2f}")
    if V_max < np.inf:
        plt.axvline(V_max, color="blue", linestyle=":", label=r"$\mathcal{V}^+$: " + f"{V_max:.2f}")
    
    plt.title("Selective Inference: Truncated Normal Distribution")
    plt.xlabel(r"Contrast Value ($\phi$)")
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "selective_inference_truncation.png"), dpi=150)
    plt.close()
    
    logging.info(f"Plots saved to {out_dir}")

if __name__ == "__main__":
    main()
