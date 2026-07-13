#!/usr/bin/env python3
"""
Exploration of single-cell RNA-seq reference data in SingleCellSubclustered.
Triggers sub-clustering, cluster means, and bulk deconvolution derivation,
and plots fraction distribution dashboards for all cohorts.
Includes splits by response and TMB.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad

# Ensure local packages are in python path if running directly
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))

import datalair
from single_cell_datasets._single_cell_reference import (
    SingleCellSubclustered,
    SingleCellClusterMeans,
    SingleCellDeconvolution
)


def load_clinical_and_tmb(cohort_name: str, lair):
    """
    Load clinical and TMB data for a given cohort directly to bypass library bugs.
    """
    import ici_datasets
    from gene_utils import calculate_maf_tmb
    
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    try:
        data_dir = ici_datasets.cbioportal_datasets.get_dataset_dir(lair, dataset_class, cohort_name)
        
        filepath_clinical_sample = data_dir / "data_clinical_sample.txt"
        df_clinical = pd.read_csv(filepath_clinical_sample, sep="\t", skiprows=4)
        if "SAMPLE_ID" in df_clinical.columns:
            df_clinical = df_clinical.set_index("SAMPLE_ID")
        elif "sample_id" in df_clinical.columns:
            df_clinical = df_clinical.set_index("sample_id")
        else:
            df_clinical.index = df_clinical.index.astype(str)
            
        # Clean index: convert to string, remove duplicates
        df_clinical.index = df_clinical.index.astype(str)
        df_clinical = df_clinical[~df_clinical.index.duplicated(keep='first')]
        
        # Resolve response column
        RESPONSE_COL_CANDIDATES = [
            "RESPONSE", "Response", "response",
            "BEST_RESPONSE", "Best_Response", "best_response",
            "RECIST", "recist",
        ]
        response_col = None
        for col in df_clinical.columns:
            if col.strip().lower() in [c.lower() for c in RESPONSE_COL_CANDIDATES]:
                response_col = col
                break
                
        if response_col is None:
            return None, None
            
        RESPONSE_VALUE_MAP = {
            "Complete Response": "R", "Partial Response": "R",
            "Stable Disease": "NR", "Progressive Disease": "NR",
            "CR": "R", "PR": "R", "SD": "NR", "PD": "NR",
            "R": "R", "NR": "NR"
        }
        df_clinical['response'] = df_clinical[response_col].map(RESPONSE_VALUE_MAP)
        df_clinical = df_clinical[df_clinical['response'].isin(['R', 'NR'])]
        
        mut_file = data_dir / "data_mutations.txt"
        if not mut_file.exists():
            return df_clinical, pd.Series(0.0, index=df_clinical.index)
            
        df_mut = pd.read_csv(mut_file, sep="\t", low_memory=False)
        df_mut['Tumor_Sample_Barcode'] = df_mut['Tumor_Sample_Barcode'].astype(str)
        
        tmb_matrix = calculate_maf_tmb(df_mut, capture_size_mb=30.0, vaf_threshold=0.05)
        tmb_matrix = tmb_matrix.set_index('Tumor_Sample_Barcode')
        tmb_matrix = tmb_matrix[~tmb_matrix.index.duplicated(keep='first')]
        tmb_scores = tmb_matrix['TMB_Score'].reindex(df_clinical.index, fill_value=0.0)
        
        return df_clinical, tmb_scores
    except Exception as e:
        print(f"Error loading clinical/TMB for {cohort_name}: {e}")
        return None, None


def plot_grouped_boxes(ax, df_merged, cell_types, group_col, group_labels, group_colors):
    """
    Helper to plot custom side-by-side grouped boxplots for cell fractions.
    """
    y_positions = np.arange(len(cell_types))
    
    for g_idx, group in enumerate(group_labels):
        df_group = df_merged[df_merged[group_col] == group]
        data_list = [df_group[ct].values if len(df_group) > 0 else np.array([0.0]) for ct in cell_types]
        positions = y_positions + (g_idx - 0.5) * 0.4
        
        bp = ax.boxplot(
            data_list,
            positions=positions,
            vert=False,
            widths=0.3,
            patch_artist=True,
            manage_ticks=False,
            flierprops=dict(markersize=1, marker='o', alpha=0.3, markeredgecolor='none')
        )
        
        for patch in bp['boxes']:
            patch.set_facecolor(group_colors[g_idx])
            patch.set_alpha(0.7)
            patch.set_edgecolor('black')
            patch.set_linewidth(0.8)
            
        for median in bp['medians']:
            median.set_color('#2C3E50')
            median.set_linewidth(1.2)
            
    ax.set_yticks(y_positions)
    ax.set_yticklabels(cell_types)
    
    from matplotlib.patches import Patch
    legend_patches = [Patch(facecolor=group_colors[i], edgecolor='black', alpha=0.7, label=str(group_labels[i])) for i in range(len(group_labels))]
    ax.legend(handles=legend_patches, loc='upper right', fontsize=8)


def plot_metadata_distributions(obs_df, output_dir):
    """
    Plot cell counts and QC distributions across cancer types.
    """
    print("Generating QC and metadata distribution plots...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), dpi=300)
    
    cancer_counts = obs_df['cancer_code'].value_counts()
    sns.barplot(x=cancer_counts.values, y=cancer_counts.index, ax=axes[0], hue=cancer_counts.index, palette='viridis', legend=False)
    axes[0].set_title("Cell Count per Cancer Type", fontsize=11, fontweight='bold')
    axes[0].set_xlabel("Number of Cells")
    axes[0].set_ylabel("Cancer Code")
    axes[0].grid(True, linestyle='--', alpha=0.5, axis='x')
    
    obs_df['log_total_counts'] = np.log10(obs_df['total_counts'] + 1)
    sns.boxplot(data=obs_df, y='cancer_code', x='log_total_counts', ax=axes[1], hue='cancer_code', palette='plasma', legend=False)
    axes[1].set_title("Total UMI Counts per Cell", fontsize=11, fontweight='bold')
    axes[1].set_xlabel("Log10(Total Counts)")
    axes[1].set_ylabel("")
    axes[1].grid(True, linestyle='--', alpha=0.5, axis='x')
    
    sns.boxplot(data=obs_df, y='cancer_code', x='pct_counts_mt', ax=axes[2], hue='cancer_code', palette='inferno', legend=False)
    axes[2].set_title("Percent Mitochondrial Counts", fontsize=11, fontweight='bold')
    axes[2].set_xlabel("Mitochondrial %")
    axes[2].set_ylabel("")
    axes[2].grid(True, linestyle='--', alpha=0.5, axis='x')
    
    plt.suptitle("SingleCellReference QC and Cohort Distributions", fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(output_dir / "qc_metadata_distributions.svg", format='svg', bbox_inches='tight')
    plt.close()


def plot_umap_granularities(adata, output_dir):
    """
    Plot UMAP coordinates colored by Leiden clusters at 4 resolutions.
    """
    print("Generating UMAP cluster granularity plots...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 12), dpi=300)
    resolutions = [0.5, 1.0, 1.5, 2.0]
    
    for idx, res in enumerate(resolutions):
        ax = axes[idx // 2, idx % 2]
        res_col = f"leiden_res_{res}"
        sc.pl.umap(adata, color=res_col, ax=ax, show=False, title=f"Leiden (Resolution = {res})", legend_loc='none')
        
    plt.suptitle("Single-Cell Clustering at Different Granularities", fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(output_dir / "umap_granularities.svg", format='svg', bbox_inches='tight')
    plt.close()


def plot_umap_tumor_inference(adata, output_dir):
    """
    Plot UMAP coordinates showing inferred Tumor/Normal classifications.
    """
    print("Generating UMAP tumor inference plots...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 12), dpi=300)
    resolutions = [0.5, 1.0, 1.5, 2.0]
    
    for idx, res in enumerate(resolutions):
        ax = axes[idx // 2, idx % 2]
        res_col = f"tumor_status_res_{res}"
        sc.pl.umap(adata, color=res_col, ax=ax, show=False, title=f"Inferred Tumor (Res = {res})", palette={'Tumor': '#E74C3C', 'Normal': '#3498DB'})
        
    plt.suptitle("Inferred Tumor Clusters across Resolutions", fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(output_dir / "umap_tumor_inference.svg", format='svg', bbox_inches='tight')
    plt.close()


def plot_umap_annotations(adata, output_dir):
    """
    Plot UMAP coordinates colored by automatic subtype annotations.
    """
    print("Generating UMAP cell type annotation plots...")
    fig, axes = plt.subplots(2, 2, figsize=(16, 14), dpi=300)
    resolutions = [0.5, 1.0, 1.5, 2.0]
    
    for idx, res in enumerate(resolutions):
        ax = axes[idx // 2, idx % 2]
        res_col = f"cell_type_subtype_res_{res}"
        sc.pl.umap(adata, color=res_col, ax=ax, show=False, title=f"Lineage Subtypes (Res = {res})", legend_loc='none')
        
    plt.suptitle("Hierarchical Cell Subtype Annotations across Resolutions", fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(output_dir / "umap_cell_type_annotations.svg", format='svg', bbox_inches='tight')
    plt.close()


def plot_umap_subclusters(adata, output_dir):
    """
    Plot UMAP coordinates showing K-means sub-clusters for resolutions 0.5 and 1.0.
    """
    print("Generating UMAP sub-clustering plots...")
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), dpi=300)
    
    sc.pl.umap(adata, color="kmeans_subcluster_res_0.5", ax=axes[0], show=False, title="K-Means Sub-clusters (Leiden Res = 0.5)", legend_loc='none')
    sc.pl.umap(adata, color="kmeans_subcluster_res_1.0", ax=axes[1], show=False, title="K-Means Sub-clusters (Leiden Res = 1.0)", legend_loc='none')
    
    plt.suptitle("K-Means Sub-Clustering across Resolutions", fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(output_dir / "umap_subclusters.svg", format='svg', bbox_inches='tight')
    plt.close()


def plot_subclustering_metrics(df_metrics, output_dir):
    """
    Plot sub-clustering quality metrics and Witten selective p-values across resolutions.
    """
    print("Generating sub-clustering metrics distribution plots...")
    df_metrics = df_metrics.dropna(subset=['Silhouette_Score', 'Calinski_Harabasz_Index', 'Davies_Bouldin_Index', 'Witten_p_value']).copy()
    df_metrics['neg_log_p'] = -np.log10(df_metrics['Witten_p_value'] + 1e-30)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=300)
    
    sns.boxplot(data=df_metrics, x='Resolution', y='Silhouette_Score', ax=axes[0, 0], hue='Resolution', palette='Set2', legend=False)
    sns.stripplot(data=df_metrics, x='Resolution', y='Silhouette_Score', ax=axes[0, 0], color='black', alpha=0.3, jitter=0.2)
    axes[0, 0].set_title("K-Means Silhouette Scores", fontsize=11, fontweight='bold')
    axes[0, 0].set_xlabel("Leiden Resolution")
    axes[0, 0].set_ylabel("Silhouette Coefficient")
    axes[0, 0].grid(True, linestyle='--', alpha=0.5, axis='y')
    
    sns.boxplot(data=df_metrics, x='Resolution', y='Davies_Bouldin_Index', ax=axes[0, 1], hue='Resolution', palette='Set2', legend=False)
    sns.stripplot(data=df_metrics, x='Resolution', y='Davies_Bouldin_Index', ax=axes[0, 1], color='black', alpha=0.3, jitter=0.2)
    axes[0, 1].set_title("Davies-Bouldin Index (Lower is Better)", fontsize=11, fontweight='bold')
    axes[0, 1].set_xlabel("Leiden Resolution")
    axes[0, 1].set_ylabel("Davies-Bouldin Index")
    axes[0, 1].grid(True, linestyle='--', alpha=0.5, axis='y')
    
    sns.boxplot(data=df_metrics, x='Resolution', y='Calinski_Harabasz_Index', ax=axes[1, 0], hue='Resolution', palette='Set2', legend=False)
    sns.stripplot(data=df_metrics, x='Resolution', y='Calinski_Harabasz_Index', ax=axes[1, 0], color='black', alpha=0.3, jitter=0.2)
    axes[1, 0].set_title("Calinski-Harabasz Index (Higher is Better)", fontsize=11, fontweight='bold')
    axes[1, 0].set_xlabel("Leiden Resolution")
    axes[1, 0].set_ylabel("Calinski-Harabasz Index")
    axes[1, 0].grid(True, linestyle='--', alpha=0.5, axis='y')
    
    sns.violinplot(data=df_metrics, x='Resolution', y='neg_log_p', ax=axes[1, 1], hue='Resolution', palette='Set2', inner='quartile', legend=False)
    axes[1, 1].axhline(y=-np.log10(0.05), color='#E74C3C', linestyle='--', label='Alpha = 0.05')
    axes[1, 1].set_title("Witten Selective Inference Significance", fontsize=11, fontweight='bold')
    axes[1, 1].set_xlabel("Leiden Resolution")
    axes[1, 1].set_ylabel("-Log10(Selective p-value)")
    axes[1, 1].legend(loc='upper right', fontsize=8.5)
    axes[1, 1].grid(True, linestyle='--', alpha=0.5, axis='y')
    
    plt.suptitle("K-Means Sub-Clustering Evaluation Metrics across Resolutions", fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(output_dir / "subclustering_metrics_distributions.svg", format='svg', bbox_inches='tight')
    plt.close()


def plot_cluster_means_summary(means_stats, output_dir):
    """
    Plot the number of clusters in each generated cluster means file.
    """
    print("Generating cluster means summary plot...")
    plt.figure(figsize=(10, 5), dpi=300)
    
    df_stats = pd.DataFrame(means_stats)
    df_stats['label'] = df_stats['file_name'].str.replace('means_', '').str.replace('.h5ad', '')
    
    colors = ['#1ABC9C', '#2ECC71', '#3498DB', '#9B59B6', '#E67E22', '#E74C3C', '#34495E', '#16A085']
    bars = plt.barh(df_stats['label'], df_stats['num_clusters'], color=colors, edgecolor='black', height=0.6)
    
    plt.xlabel("Number of Clusters / Profiles", fontsize=10, fontweight='bold')
    plt.title("Summary: Number of Exported Profiles in Cluster Means Datasets", fontsize=12, fontweight='bold', pad=15)
    
    for bar in bars:
        width = bar.get_width()
        plt.text(width + 0.5, bar.get_y() + bar.get_height()/2.0, f"{int(width)}", 
                 va='center', ha='left', fontsize=9, fontweight='bold', color='#2C3E50')
                 
    plt.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5, axis='x')
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir / "cluster_means_summary.svg", format='svg', bbox_inches='tight')
    plt.close()


def plot_deconvolution_distributions(deconv_paths, clustering_name, output_dir):
    """
    Plot full cell type fraction distributions for each cohort in a grid layout,
    comparing each iAtlas cohort side-by-side with its matching TCGA cohort.
    """
    print(f"Generating cell fraction distribution plots for {clustering_name} (with TCGA comparison)...")
    
    COHORT_TO_TCGA = {
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
    
    cohorts = list(COHORT_TO_TCGA.keys())
    n_cohorts = len(cohorts)
    cols = 3
    rows = int(np.ceil(n_cohorts / cols))
    
    # Using 12 * rows for extra vertical breathing room
    fig, axes = plt.subplots(rows, cols, figsize=(24, 12 * rows), dpi=300)
    axes_flat = axes.flatten() if n_cohorts > 1 else [axes]
    
    for idx, cohort_name in enumerate(cohorts):
        tcga_project = COHORT_TO_TCGA[cohort_name]
        ax = axes_flat[idx]
        
        # Load iAtlas and TCGA files
        iatlas_key = f"deconv_{cohort_name}_{clustering_name}.csv"
        tcga_key = f"deconv_TCGA-{tcga_project}_{clustering_name}.csv"
        
        iatlas_path = deconv_paths.get(iatlas_key)
        tcga_path = deconv_paths.get(tcga_key)
        
        if not iatlas_path:
            ax.text(0.5, 0.5, f"Missing iAtlas data for {cohort_name}", ha='center', va='center')
            continue
            
        df_iatlas = pd.read_csv(iatlas_path, index_col=0)
        df_iatlas.index = df_iatlas.index.astype(str)
        df_iatlas['cohort_group'] = cohort_name
        
        if tcga_path:
            df_tcga = pd.read_csv(tcga_path, index_col=0)
            df_tcga.index = df_tcga.index.astype(str)
            df_tcga['cohort_group'] = f"TCGA-{tcga_project}"
            df_merged = pd.concat([df_iatlas, df_tcga], axis=0)
            group_labels = [cohort_name, f"TCGA-{tcga_project}"]
            group_colors = ['#2ECC71', '#9B59B6']  # Green vs Purple
        else:
            df_merged = df_iatlas
            group_labels = [cohort_name]
            group_colors = ['#2ECC71']
            
        cell_types = df_iatlas.columns.drop('cohort_group', errors='ignore')
        
        plot_grouped_boxes(
            ax=ax,
            df_merged=df_merged,
            cell_types=cell_types,
            group_col='cohort_group',
            group_labels=group_labels,
            group_colors=group_colors
        )
        
        ax.set_title(f"{cohort_name} vs TCGA-{tcga_project}", fontsize=11, fontweight='bold', color='#2C3E50')
        ax.set_xlabel("Fractions")
        ax.set_ylabel("Clusters")
        ax.tick_params(axis='y', labelsize=8)
        ax.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
        
    for ax in axes_flat[idx + 1:]:
        ax.axis('off')
        
    plt.suptitle(f"Comparison of Inferred Cell Fractions: iAtlas Trial vs TCGA Baseline ({clustering_name})", fontsize=15, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(output_dir / f"deconv_distributions_{clustering_name}.svg", format='svg', bbox_inches='tight')
    plt.close()


def plot_deconvolution_split_response(deconv_paths, clustering_name, output_dir, lair):
    """
    Plot cell fraction distributions split by Response (R vs NR).
    """
    print(f"Generating cell fraction distributions split by response for {clustering_name}...")
    
    filtered_paths = {k: v for k, v in deconv_paths.items() if clustering_name in k}
    if not filtered_paths:
        return
        
    n_cohorts = len(filtered_paths)
    cols = 3
    rows = int(np.ceil(n_cohorts / cols))
    
    fig, axes = plt.subplots(rows, cols, figsize=(24, 12 * rows), dpi=300)
    axes_flat = axes.flatten() if n_cohorts > 1 else [axes]
    
    for idx, (file_key, file_path) in enumerate(sorted(list(filtered_paths.items()))):
        cohort_name = file_key.split('_')[1]
        df_fracs = pd.read_csv(file_path, index_col=0)
        df_fracs.index = df_fracs.index.astype(str)
        
        # Load clinical response
        df_clin, _ = load_clinical_and_tmb(cohort_name, lair)
        ax = axes_flat[idx]
        
        if df_clin is None:
            ax.text(0.5, 0.5, f"No clinical data for {cohort_name}", ha='center', va='center')
            continue
            
        common_samples = df_fracs.index.intersection(df_clin.index)
        if len(common_samples) == 0:
            ax.text(0.5, 0.5, f"No overlapping samples for {cohort_name}", ha='center', va='center')
            continue
            
        df_fracs_aligned = df_fracs.loc[common_samples]
        df_clin_aligned = df_clin.loc[common_samples]
        df_merged = pd.concat([df_fracs_aligned, df_clin_aligned['response']], axis=1)
        
        cell_types = df_fracs.columns
        plot_grouped_boxes(
            ax=ax,
            df_merged=df_merged,
            cell_types=cell_types,
            group_col='response',
            group_labels=['R', 'NR'],
            group_colors=['#2ECC71', '#E74C3C'] # Green vs Red
        )
        
        ax.set_title(f"{cohort_name} Split by Response", fontsize=11, fontweight='bold', color='#2C3E50')
        ax.set_xlabel("Fractions")
        ax.set_ylabel("Clusters")
        ax.tick_params(axis='y', labelsize=8)
        ax.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
        
    for ax in axes_flat[idx + 1:]:
        ax.axis('off')
        
    plt.suptitle(f"Cell Fractions Split by Response (R vs NR) - {clustering_name}", fontsize=15, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(output_dir / f"deconv_split_response_{clustering_name}.svg", format='svg', bbox_inches='tight')
    plt.close()


def plot_deconvolution_split_tmb(deconv_paths, clustering_name, output_dir, lair):
    """
    Plot cell fraction distributions split by TMB status (High vs Low TMB, cohort median split).
    """
    print(f"Generating cell fraction distributions split by TMB for {clustering_name}...")
    
    filtered_paths = {k: v for k, v in deconv_paths.items() if clustering_name in k}
    if not filtered_paths:
        return
        
    n_cohorts = len(filtered_paths)
    cols = 3
    rows = int(np.ceil(n_cohorts / cols))
    
    fig, axes = plt.subplots(rows, cols, figsize=(24, 12 * rows), dpi=300)
    axes_flat = axes.flatten() if n_cohorts > 1 else [axes]
    
    for idx, (file_key, file_path) in enumerate(sorted(list(filtered_paths.items()))):
        cohort_name = file_key.split('_')[1]
        df_fracs = pd.read_csv(file_path, index_col=0)
        df_fracs.index = df_fracs.index.astype(str)
        
        # Load clinical and TMB
        df_clin, tmb_scores = load_clinical_and_tmb(cohort_name, lair)
        ax = axes_flat[idx]
        
        if df_clin is None or tmb_scores is None or tmb_scores.sum() == 0:
            ax.text(0.5, 0.5, f"No TMB/mutation data for {cohort_name}", ha='center', va='center')
            continue
            
        common_samples = df_fracs.index.intersection(tmb_scores.index)
        if len(common_samples) == 0:
            ax.text(0.5, 0.5, f"No overlapping samples for {cohort_name}", ha='center', va='center')
            continue
            
        df_fracs_aligned = df_fracs.loc[common_samples]
        tmb_aligned = tmb_scores.loc[common_samples]
        
        # Define High vs Low TMB based on median split
        median_val = tmb_aligned.median()
        tmb_status = pd.Series(
            np.where(tmb_aligned >= median_val, 'High TMB', 'Low TMB'),
            index=common_samples
        )
        
        df_merged = pd.concat([df_fracs_aligned, tmb_status.rename('tmb_group')], axis=1)
        cell_types = df_fracs.columns
        
        plot_grouped_boxes(
            ax=ax,
            df_merged=df_merged,
            cell_types=cell_types,
            group_col='tmb_group',
            group_labels=['High TMB', 'Low TMB'],
            group_colors=['#3498DB', '#E67E22'] # Blue vs Orange
        )
        
        ax.set_title(f"{cohort_name} Split by TMB (Median={median_val:.2f})", fontsize=11, fontweight='bold', color='#2C3E50')
        ax.set_xlabel("Fractions")
        ax.set_ylabel("Clusters")
        ax.tick_params(axis='y', labelsize=8)
        ax.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
        
    for ax in axes_flat[idx + 1:]:
        ax.axis('off')
        
    plt.suptitle(f"Cell Fractions Split by TMB (High vs Low TMB) - {clustering_name}", fontsize=15, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(output_dir / f"deconv_split_tmb_{clustering_name}.svg", format='svg', bbox_inches='tight')
    plt.close()


def plot_marker_genes(adata, output_dir):
    """
    Plot lineage marker expression across Leiden clusters (Res = 0.5) using dotplot.
    """
    print("Plotting lineage marker gene expression dotplot...")
    lineage_markers = {
        'T/NK cells': ['CD3D', 'CD3E', 'CD8A', 'NKG7'],
        'B/Plasma cells': ['MS4A1', 'MZB1'],
        'Myeloid cells': ['CD68', 'LYZ'],
        'Epithelial/Tumor': ['EPCAM', 'KRT18'],
        'Endothelial cells': ['PECAM1'],
        'Fibroblasts': ['COL1A1']
    }
    
    sc.pl.dotplot(adata, lineage_markers, groupby='leiden_res_0.5', show=False, use_raw=True)
    plt.savefig(output_dir / "marker_gene_dotplot.svg", format='svg', bbox_inches='tight')
    plt.close()


def main():
    lair_path = "/storage/halu/lair"
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output" / "single-cell-exploration"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Derive Sub-clustered AnnData
    ds = SingleCellSubclustered()
    print("Loading/deriving SingleCellSubclustered dataset in Lair...")
    lair.safe_derive(ds, overwrite=False)
    filepaths = lair.get_dataset_filepaths(ds)
    
    print("Loading derived sub-clustered AnnData...")
    adata = ad.read_h5ad(filepaths["adata.h5ad"])
    print(f"AnnData successfully loaded. Dimensions: {adata.shape}")
    
    # 2. Derive Cluster Means AnnData
    ds_means = SingleCellClusterMeans()
    print("\nLoading/deriving SingleCellClusterMeans dataset in Lair...")
    lair.safe_derive(ds_means, overwrite=False)
    means_paths = lair.get_dataset_filepaths(ds_means)
    
    means_stats = []
    for file_key, file_path in sorted(list(means_paths.items())):
        adata_mean = ad.read_h5ad(file_path)
        means_stats.append({
            'file_name': file_key,
            'num_clusters': adata_mean.n_obs,
            'num_genes': adata_mean.n_vars
        })
        
    plot_cluster_means_summary(means_stats, output_dir)
    
    # 3. Derive Single-Cell Deconvolution fractions
    ds_deconv = SingleCellDeconvolution()
    print("\nLoading/deriving SingleCellDeconvolution dataset in Lair...")
    # Overwrite to derive TCGA and iAtlas deconvolution
    lair.safe_derive(ds_deconv, overwrite=True)
    deconv_paths = lair.get_dataset_filepaths(ds_deconv)
    print(f"Deconvolution resolved files: {len(deconv_paths)} CSV files.")
    
    # Plot fraction distributions for the main resolutions and subclusters
    plot_deconvolution_distributions(deconv_paths, "leiden_res_0.5", output_dir)
    plot_deconvolution_distributions(deconv_paths, "kmeans_subcluster_res_0.5", output_dir)
    
    # Plot split distributions
    plot_deconvolution_split_response(deconv_paths, "leiden_res_0.5", output_dir, lair)
    plot_deconvolution_split_response(deconv_paths, "kmeans_subcluster_res_0.5", output_dir, lair)
    
    plot_deconvolution_split_tmb(deconv_paths, "leiden_res_0.5", output_dir, lair)
    plot_deconvolution_split_tmb(deconv_paths, "kmeans_subcluster_res_0.5", output_dir, lair)
    
    # Load sub-clustering metrics CSV
    metrics_csv = filepaths["adata.h5ad"].parent / "subclustering_metrics.csv"
    if metrics_csv.exists():
        df_metrics = pd.read_csv(metrics_csv)
        plot_subclustering_metrics(df_metrics, output_dir)
        df_metrics.to_csv(output_dir / "subclustering_metrics.csv", index=False)
        print(f"Subclustering metrics table saved to output directory.")
        
    # Generate distribution plots
    plot_metadata_distributions(adata.obs, output_dir)
    
    # Generate UMAP granularities
    plot_umap_granularities(adata, output_dir)
    
    # Generate UMAP tumor inference
    plot_umap_tumor_inference(adata, output_dir)
    
    # Generate UMAP cell type annotations
    plot_umap_annotations(adata, output_dir)
    
    # Generate UMAP K-means sub-clusters
    plot_umap_subclusters(adata, output_dir)
    
    # Generate marker gene dotplot
    plot_marker_genes(adata, output_dir)
    
    # Save coordinate mapping CSV
    umap_df = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names, columns=['UMAP1', 'UMAP2'])
    metadata_cols = ['cancer_code', 'staging', 'dataset', 'cnv_score',
                     'leiden_res_0.5', 'cell_type_subtype_res_0.5', 'tumor_status_res_0.5', 'kmeans_subcluster_res_0.5',
                     'leiden_res_1.0', 'cell_type_subtype_res_1.0', 'tumor_status_res_1.0', 'kmeans_subcluster_res_1.0',
                     'leiden_res_1.5', 'cell_type_subtype_res_1.5', 'tumor_status_res_1.5', 'kmeans_subcluster_res_1.5',
                     'leiden_res_2.0', 'cell_type_subtype_res_2.0', 'tumor_status_res_2.0', 'kmeans_subcluster_res_2.0']
    
    metadata_df = pd.concat([adata.obs[metadata_cols], umap_df], axis=1)
    metadata_csv_path = output_dir / "subsampled_cell_metadata_embeddings.csv"
    metadata_df.to_csv(metadata_csv_path)
    print(f"Saved complete cell metadata, clusterings, and annotations to {metadata_csv_path.name}")
    
    print("\nSingle-cell clustered exploration visualization complete!")


if __name__ == "__main__":
    main()
