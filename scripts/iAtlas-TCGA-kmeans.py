#!/usr/bin/env python3
"""
TCGA-iAtlas Cell Fraction K-Means (K=100) Clustering and Classification Analysis.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import confusion_matrix, accuracy_score, classification_report
import umap

# Ensure local packages are in path
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))
import datalair
from single_cell_datasets._single_cell_reference import SingleCellDeconvolution


def main():
    lair_path = "/storage/halu/lair"
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output" / "single-cell-exploration"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Load SingleCellDeconvolution filepaths
    ds_deconv = SingleCellDeconvolution()
    deconv_paths = lair.get_dataset_filepaths(ds_deconv)
    
    resolution = "leiden_res_0.5"
    
    # TCGA cohorts and true cancer type mapping
    TCGA_MAP = {
        'TCGA-SKCM': 'SKCM',
        'TCGA-BLCA': 'BLCA',
        'TCGA-PAAD': 'PAAD',
        'TCGA-BRCA': 'BRCA',
        'TCGA-KIRC': 'KIRC'
    }
    
    # iAtlas trials and actual cancer type mapping
    IATLAS_MAP = {
        'Hugo-iAtlas': 'SKCM',
        'Riaz-iAtlas': 'SKCM',
        'Liu-iAtlas': 'SKCM',
        'Gide-iAtlas': 'SKCM',
        'Rosenberg-iAtlas': 'BLCA',
        'Padron-iAtlas': 'PAAD',
        'Anders-iAtlas': 'BRCA',
        'McDermott-iAtlas': 'KIRC',
        'Choueiri-iAtlas': 'KIRC'
    }
    
    # Color palette for the 5 cancer types
    CANCER_COLORS = {
        'SKCM': '#E74C3C', # Red
        'BLCA': '#2ECC71', # Green
        'PAAD': '#3498DB', # Blue
        'BRCA': '#E67E22', # Orange
        'KIRC': '#9B59B6'  # Purple
    }
    
    # 2. Collect TCGA data
    tcga_dfs = []
    for cohort, cancer_type in TCGA_MAP.items():
        key = f"deconv_{cohort}_{resolution}.csv"
        path = deconv_paths.get(key)
        if not path:
            print(f"Warning: Missing file for {cohort}. Skipping.")
            continue
        df = pd.read_csv(path, index_col=0)
        df.index = df.index.astype(str)
        df['cancer_type'] = cancer_type
        df['cohort'] = cohort
        tcga_dfs.append(df)
        
    if not tcga_dfs:
        print("Error: No TCGA data found.")
        return
        
    df_tcga_all = pd.concat(tcga_dfs, axis=0)
    features = [c for c in df_tcga_all.columns if c not in ['cancer_type', 'cohort']]
    X_tcga = df_tcga_all[features].values
    
    print(f"\nLoaded {len(df_tcga_all)} TCGA samples with {len(features)} cell type fractions.")
    
    # 3. Fit K-Means on TCGA (K=100)
    K = 100
    print(f"Running K-Means (K={K}) clustering on TCGA baseline samples...")
    kmeans = KMeans(n_clusters=K, random_state=42, n_init=10)
    kmeans_labels = kmeans.fit_predict(X_tcga)
    df_tcga_all['cluster'] = kmeans_labels
    
    # 4. Analyze cluster centroids
    centroids = kmeans.cluster_centers_
    cluster_metadata = []
    
    for c in range(K):
        mask = df_tcga_all['cluster'] == c
        df_cluster = df_tcga_all[mask]
        size = len(df_cluster)
        
        if size > 0:
            counts = df_cluster['cancer_type'].value_counts()
            dominant_type = counts.index[0]
            purity = counts.iloc[0] / size
        else:
            dominant_type = 'Unknown'
            purity = 0.0
            
        cluster_metadata.append({
            'cluster': c,
            'size': size,
            'dominant_cancer_type': dominant_type,
            'purity': purity
        })
        
    df_clusters = pd.DataFrame(cluster_metadata)
    
    print("\nCentroid properties (top 5 largest clusters):")
    print(df_clusters.sort_values(by='size', ascending=False).head(5))
    
    # Check if there are empty or very small clusters
    min_size = df_clusters['size'].min()
    print(f"Min cluster size: {min_size} samples. Max cluster size: {df_clusters['size'].max()} samples.")
    
    # 5. Plot 1: K-Means cluster sizes and dominant cancer types
    print("Plotting cluster size distribution...")
    plt.figure(figsize=(15, 6), dpi=300)
    df_clusters_sorted = df_clusters.sort_values(by='size', ascending=False)
    
    bar_colors = [CANCER_COLORS.get(t, '#95A5A6') for t in df_clusters_sorted['dominant_cancer_type']]
    
    plt.bar(
        range(K), df_clusters_sorted['size'],
        color=bar_colors, edgecolor='black', linewidth=0.5
    )
    plt.xticks(range(K), df_clusters_sorted['cluster'], rotation=90, fontsize=6)
    plt.xlabel("K-Means Cluster ID (Sorted by Size)", fontsize=10, fontweight='bold')
    plt.ylabel("Number of Samples", fontsize=10, fontweight='bold')
    plt.title(f"TCGA Cluster Size Distribution & Dominant Cancer Type (K={K})", fontsize=12, fontweight='bold', pad=15)
    
    # Add custom legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=color, edgecolor='black', label=f"{name} (Dominant)")
        for name, color in CANCER_COLORS.items()
    ]
    plt.legend(handles=legend_elements, loc='upper right', fontsize=9)
    plt.grid(True, linestyle='--', alpha=0.3, axis='y')
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir / "kmeans_cluster_sizes.svg", format='svg', bbox_inches='tight')
    plt.close()
    
    # 6. Plot 2: Heatmap of centroids
    print("Plotting centroids heatmap...")
    df_centroids = pd.DataFrame(centroids, columns=features)
    df_centroids['dominant_cancer_type'] = df_clusters['dominant_cancer_type']
    
    # Sort centroids by dominant cancer type to group them in the heatmap
    df_centroids_sorted = df_centroids.sort_values(by='dominant_cancer_type')
    row_labels = [f"C{r} ({t})" for r, t in zip(df_centroids_sorted.index, df_centroids_sorted['dominant_cancer_type'])]
    
    plt.figure(figsize=(12, 18), dpi=300)
    sns.heatmap(
        df_centroids_sorted[features],
        yticklabels=row_labels,
        cmap='viridis',
        cbar_kws={'label': 'Mean Cell Type Fraction'}
    )
    plt.title(f"Heatmap of K-Means Cluster Centroids (K={K}, Sorted by Dominant Cancer Type)", fontsize=13, fontweight='bold', pad=15)
    plt.xlabel("Cell Type (Resolution 0.5)", fontsize=10, fontweight='bold')
    plt.ylabel("Cluster Centroid (Dominant Type)", fontsize=10, fontweight='bold')
    plt.tick_params(axis='y', labelsize=6)
    plt.tight_layout()
    plt.savefig(output_dir / "kmeans_centroids_heatmap.svg", format='svg', bbox_inches='tight')
    plt.close()
    
    # 7. Plot 3: TCGA UMAP colored by True Cancer Type and K-Means Cluster
    print("Computing UMAP on TCGA cell fractions...")
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
    embedding_tcga = reducer.fit_transform(X_tcga)
    
    fig, axes = plt.subplots(1, 2, figsize=(18, 8), dpi=300)
    
    # True Cancer Type
    for name, color in CANCER_COLORS.items():
        mask = df_tcga_all['cancer_type'] == name
        axes[0].scatter(
            embedding_tcga[mask, 0], embedding_tcga[mask, 1],
            c=color, label=name, s=15, alpha=0.6, edgecolors='none'
        )
    axes[0].set_title("TCGA Samples by True Cancer Type", fontsize=12, fontweight='bold')
    axes[0].legend(loc='upper right', markerscale=1.5)
    axes[0].grid(True, linestyle='--', alpha=0.3)
    
    # K-Means Cluster ID
    scatter = axes[1].scatter(
        embedding_tcga[:, 0], embedding_tcga[:, 1],
        c=kmeans_labels, cmap='tab20', s=15, alpha=0.6, edgecolors='none'
    )
    axes[1].set_title(f"TCGA Samples by K-Means Cluster ID (K={K})", fontsize=12, fontweight='bold')
    axes[1].grid(True, linestyle='--', alpha=0.3)
    
    plt.suptitle("UMAP Projection of TCGA Deconvolution Cell Fractions", fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(output_dir / "tcga_umap.svg", format='svg', bbox_inches='tight')
    plt.close()
    
    # 8. Load iAtlas trials and classify samples
    print("\nLoading and classifying iAtlas trials...")
    iatlas_dfs = []
    for cohort, cancer_type in IATLAS_MAP.items():
        key = f"deconv_{cohort}_{resolution}.csv"
        path = deconv_paths.get(key)
        if not path:
            print(f"Warning: Missing file for {cohort}. Skipping.")
            continue
        df = pd.read_csv(path, index_col=0)
        df.index = df.index.astype(str)
        df['actual_cancer_type'] = cancer_type
        df['cohort'] = cohort
        iatlas_dfs.append(df)
        
    if not iatlas_dfs:
        print("Error: No iAtlas data found.")
        return
        
    df_iatlas_all = pd.concat(iatlas_dfs, axis=0)
    X_iatlas = df_iatlas_all[features].values
    
    print(f"Loaded {len(df_iatlas_all)} iAtlas trial samples.")
    
    # Map each iAtlas sample to the nearest K-means centroid
    # kmeans.predict assigns it to the nearest centroid
    iatlas_clusters = kmeans.predict(X_iatlas)
    df_iatlas_all['predicted_cluster'] = iatlas_clusters
    
    # Look up dominant cancer type of predicted cluster
    cluster_to_dominant = df_clusters.set_index('cluster')['dominant_cancer_type'].to_dict()
    df_iatlas_all['predicted_cancer_type'] = [cluster_to_dominant[c] for c in iatlas_clusters]
    
    # 9. Evaluate performance
    y_true = df_iatlas_all['actual_cancer_type'].values
    y_pred = df_iatlas_all['predicted_cancer_type'].values
    
    acc = accuracy_score(y_true, y_pred)
    print(f"\niAtlas Classification Accuracy: {acc:.4f}")
    print("\nClassification Report:")
    print(classification_report(y_true, y_pred, zero_division=0))
    
    # 10. Plot 4: iAtlas UMAP colored by Actual vs Predicted Cancer Type
    print("Computing UMAP on iAtlas cell fractions...")
    reducer_iatlas = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
    embedding_iatlas = reducer_iatlas.fit_transform(X_iatlas)
    
    fig, axes = plt.subplots(1, 2, figsize=(18, 8), dpi=300)
    
    # Actual Cancer Type
    for name, color in CANCER_COLORS.items():
        mask = df_iatlas_all['actual_cancer_type'] == name
        axes[0].scatter(
            embedding_iatlas[mask, 0], embedding_iatlas[mask, 1],
            c=color, label=name, s=25, alpha=0.6, edgecolors='black', linewidths=0.2
        )
    axes[0].set_title("iAtlas Trial Samples by Actual Cancer Type", fontsize=12, fontweight='bold')
    axes[0].legend(loc='upper right', markerscale=1.5)
    axes[0].grid(True, linestyle='--', alpha=0.3)
    
    # Predicted Cancer Type
    for name, color in CANCER_COLORS.items():
        mask = df_iatlas_all['predicted_cancer_type'] == name
        axes[1].scatter(
            embedding_iatlas[mask, 0], embedding_iatlas[mask, 1],
            c=color, label=name, s=25, alpha=0.6, edgecolors='black', linewidths=0.2
        )
    axes[1].set_title(f"iAtlas Trial Samples by Predicted Cancer Type (Accuracy: {acc:.2%})", fontsize=12, fontweight='bold')
    axes[1].legend(loc='upper right', markerscale=1.5)
    axes[1].grid(True, linestyle='--', alpha=0.3)
    
    plt.suptitle("UMAP Projection of iAtlas Deconvolution Cell Fractions: Actual vs Predicted", fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(output_dir / "iatlas_classification_umap.svg", format='svg', bbox_inches='tight')
    plt.close()
    
    # 11. Plot 5: Confusion Matrix
    print("Plotting confusion matrix...")
    labels = sorted(list(CANCER_COLORS.keys()))
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    
    # Normalize rows (actuals) to sum to 100%
    cm_norm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
    
    plt.figure(figsize=(8, 7), dpi=300)
    
    # Format cell text with counts and percentages
    cell_texts = []
    for r in range(len(labels)):
        row_text = []
        for c in range(len(labels)):
            val = cm[r, c]
            pct = cm_norm[r, c] * 100
            row_text.append(f"{val}\n({pct:.1f}%)")
        cell_texts.append(row_text)
        
    sns.heatmap(
        cm_norm, annot=np.array(cell_texts), fmt="", cmap="Blues",
        xticklabels=labels, yticklabels=labels, cbar_kws={'label': 'Classification Rate'}
    )
    plt.title(f"iAtlas Cancer Type Classification Confusion Matrix\n(K-Means Centroid Classifier, Accuracy: {acc:.2%})", fontsize=12, fontweight='bold', pad=15)
    plt.xlabel("Predicted Cancer Type", fontsize=10, fontweight='bold')
    plt.ylabel("Actual Cancer Type", fontsize=10, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_dir / "kmeans_confusion_matrix.svg", format='svg', bbox_inches='tight')
    plt.close()
    
    # Save classification summary
    df_iatlas_all[['cohort', 'actual_cancer_type', 'predicted_cancer_type', 'predicted_cluster']].to_csv(
        output_dir / "iatlas_classification_predictions.csv"
    )
    
    print("\nClassification predictions saved to iatlas_classification_predictions.csv.")
    print("All K-Means clustering and classification diagnostic plots successfully generated!")


if __name__ == "__main__":
    main()
