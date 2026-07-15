#!/usr/bin/env python3
"""
Investigation of Cohort Clustering in Bagaev-Deconvoluted Space.
Analyzes expression normalization differences, biological lineages, and predictive feature rankings.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
import umap

# Ensure local packages are in path
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))
import datalair
from single_cell_datasets import SingleCellBagaevDeconvolution, SingleCellDeconvolution, SingleCellClusterMeans
from ici_datasets.bagaev_datasets import Signature as BagaevSignature
from gene_utils import read_gene_sets

# Setup output dir
output_dir = Path("output/bagaev-constrained-deconv/investigation")
output_dir.mkdir(parents=True, exist_ok=True)


def load_clinical(cohort_name: str, lair):
    """Helper to load clinical response from CBioPortal."""
    import ici_datasets
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
            
        df_clinical.index = df_clinical.index.astype(str)
        df_clinical = df_clinical[~df_clinical.index.duplicated(keep='first')]
        
        RESPONSE_COL_CANDIDATES = ["RESPONSE", "Response", "response", "BEST_RESPONSE", "RECIST", "recist"]
        response_col = None
        for col in df_clinical.columns:
            if col.strip().lower() in [c.lower() for c in RESPONSE_COL_CANDIDATES]:
                response_col = col
                break
        if response_col is None:
            return None
            
        RESPONSE_VALUE_MAP = {
            "Complete Response": "R", "Partial Response": "R",
            "Stable Disease": "NR", "Progressive Disease": "NR",
            "CR": "R", "PR": "R", "SD": "NR", "PD": "NR",
            "R": "R", "NR": "NR"
        }
        df_clinical['response'] = df_clinical[response_col].map(RESPONSE_VALUE_MAP)
        df_clinical = df_clinical[df_clinical['response'].isin(['R', 'NR'])]
        return df_clinical[['response']]
    except Exception as e:
        print(f"Error loading clinical data for {cohort_name}: {e}")
        return None


def main():
    print("Initializing Lair...")
    lair = datalair.Lair("/storage/halu/lair")
    
    ds_bagaev = SingleCellBagaevDeconvolution()
    lair.safe_derive(ds_bagaev)
    bagaev_deconv_paths = lair.get_dataset_filepaths(ds_bagaev)
    
    cohorts = [
        "Hugo-iAtlas",
        "Riaz-iAtlas",
        "Liu-iAtlas",
        "Gide-iAtlas",
        "Rosenberg-iAtlas",
        "Padron-iAtlas",
        "Anders-iAtlas",
        "McDermott-iAtlas",
        "Choueiri-iAtlas"
    ]
    
    COHORT_TYPES = {
        'Hugo-iAtlas': 'Native TPM',
        'Riaz-iAtlas': 'Native TPM',
        'Liu-iAtlas': 'Native TPM',
        'Gide-iAtlas': 'Custom Transformed',
        'Rosenberg-iAtlas': 'Native TPM',
        'Padron-iAtlas': 'Custom Transformed',
        'Anders-iAtlas': 'Native TPM',
        'McDermott-iAtlas': 'Native TPM',
        'Choueiri-iAtlas': 'Custom Transformed'
    }
    
    ref_name = "kmeans_subcluster_res_0.5"
    lineages = ['T_NK', 'B_Plasma', 'Myeloid', 'Endothelial', 'Fibroblasts', 'Tumor_Epithelial']
    
    # Load cluster means annotations to map subclusters to major lineages
    import anndata as ad
    ds_means = SingleCellClusterMeans()
    means_paths = lair.get_dataset_filepaths(ds_means)
    adata_means = ad.read_h5ad(means_paths[f"means_{ref_name}.h5ad"])
    cluster_to_major = adata_means.obs['cell_type_major_res_0.5'].to_dict()
    
    # 1. TECHNICAL CHECK: Expression summaries
    print("\nPhase 1: Analyzing Expression Summaries...")
    expr_stats = []
    
    for cohort in cohorts:
        import ici_datasets
        dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
        ds_cb = dataset_class(name=cohort)
        filepaths_cb = lair.get_dataset_filepaths(ds_cb)
        unpacked_key = next(f for f in filepaths_cb.keys() if not f.endswith('.tar.gz'))
        p_dir = filepaths_cb[unpacked_key] / filepaths_cb[unpacked_key].name
        
        target_file = None
        for f in ["data_mrna_seq_expression.txt", "data_mrna_seq_tpm.txt", "data_mrna_seq_rpkm.txt"]:
            if (p_dir / f).is_file():
                target_file = p_dir / f
                break
        if not target_file:
            continue
            
        df_expr = pd.read_csv(target_file, sep='\t', index_col=0)
        df_expr.index = df_expr.index.str.upper()
        df_expr = df_expr.groupby(level=0).mean()
        
        # We calculate total sum and mean expression across all genes
        sample_sums = df_expr.sum(axis=0)
        sample_means = df_expr.mean(axis=0)
        sample_vars = df_expr.var(axis=0)
        
        for sid in df_expr.columns:
            expr_stats.append({
                'Cohort': cohort,
                'Type': COHORT_TYPES[cohort],
                'Sample_Sum': sample_sums[sid],
                'Sample_Mean': sample_means[sid],
                'Sample_Var': sample_vars[sid]
            })
            
    df_expr_stats = pd.DataFrame(expr_stats)
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    sns.boxplot(data=df_expr_stats, x='Cohort', y='Sample_Sum', hue='Type', ax=axes[0], palette='Set2')
    axes[0].set_title("Distribution of Sample Expression Sums")
    axes[0].set_yscale('log')
    axes[0].tick_params(axis='x', rotation=45)
    
    sns.boxplot(data=df_expr_stats, x='Cohort', y='Sample_Mean', hue='Type', ax=axes[1], palette='Set2')
    axes[1].set_title("Distribution of Sample Expression Means")
    axes[1].set_yscale('log')
    axes[1].tick_params(axis='x', rotation=45)
    
    sns.boxplot(data=df_expr_stats, x='Cohort', y='Sample_Var', hue='Type', ax=axes[2], palette='Set2')
    axes[2].set_title("Distribution of Sample Expression Variances")
    axes[2].set_yscale('log')
    axes[2].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_dir / "cohort_expression_stats.svg", bbox_inches='tight')
    plt.close()
    print("Saved cohort expression statistics comparison.")

    # 2. BIOLOGICAL CHECK: Mean Lineage Composition
    print("\nPhase 2: Analyzing Cell Lineage Compositions...")
    lineage_means = []
    all_cohort_fracs = []
    native_fracs_list = []
    
    for cohort in cohorts:
        bagaev_deconv_key = f"deconv_bagaev_{cohort}_{ref_name}.csv"
        deconv_path = bagaev_deconv_paths.get(bagaev_deconv_key)
        if not deconv_path:
            continue
            
        df_fracs = pd.read_csv(deconv_path, index_col=0)
        df_fracs.index = df_fracs.index.astype(str)
        
        # Calculate major lineages sums using cluster_to_major mapping
        df_lin = pd.DataFrame(index=df_fracs.index)
        for lin in lineages:
            matching_cols = [c for c, maj in cluster_to_major.items() if maj == lin and c in df_fracs.columns]
            df_lin[lin] = df_fracs[matching_cols].sum(axis=1)
            
        df_lin['Cohort'] = cohort
        df_lin['Type'] = COHORT_TYPES[cohort]
        all_cohort_fracs.append(df_lin)
        
        if COHORT_TYPES[cohort] == 'Native TPM':
            native_fracs_list.append(df_lin.drop(columns=['Cohort', 'Type']))
            
        # Calculate mean composition
        mean_comp = df_lin[lineages].mean()
        for lin in lineages:
            lineage_means.append({
                'Cohort': cohort,
                'Type': COHORT_TYPES[cohort],
                'Lineage': lin,
                'Mean_Fraction': mean_comp[lin]
            })
            
    df_lin_means = pd.DataFrame(lineage_means)
    
    # Pivot to plot stacked bar chart
    df_pivot = df_lin_means.pivot(index='Cohort', columns='Lineage', values='Mean_Fraction')[lineages]
    # Reorder index to group custom transformed together
    df_pivot = df_pivot.reindex(sorted(cohorts, key=lambda x: COHORT_TYPES[x]))
    
    df_pivot.plot(kind='bar', stacked=True, figsize=(12, 8), colormap='Set3', edgecolor='black')
    plt.title("Average Cell Lineage Compositions Across Cohorts")
    plt.ylabel("Fraction")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_dir / "cohort_lineage_abundances.svg", bbox_inches='tight')
    plt.close()
    print("Saved cohort lineage compositions stacked bar chart.")

    # 3. UMAP OF NATIVE TPM COHORTS ONLY
    print("\nPhase 3: Generating UMAP of Native TPM Cohorts Only...")
    if native_fracs_list:
        df_native_all = pd.concat(native_fracs_list, axis=0)
        native_cohort_labels = []
        for cohort in cohorts:
            if COHORT_TYPES[cohort] == 'Native TPM':
                bagaev_deconv_key = f"deconv_bagaev_{cohort}_{ref_name}.csv"
                deconv_path = bagaev_deconv_paths.get(bagaev_deconv_key)
                if deconv_path:
                    n_samples = len(pd.read_csv(deconv_path, index_col=0))
                    native_cohort_labels.extend([cohort] * n_samples)
                    
        X_native_scaled = StandardScaler().fit_transform(df_native_all.values)
        reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
        emb_native = reducer.fit_transform(X_native_scaled)
        
        plt.figure(figsize=(10, 8))
        sns.scatterplot(x=emb_native[:, 0], y=emb_native[:, 1], hue=native_cohort_labels, palette='tab10', alpha=0.8)
        plt.title("UMAP of Cell Fractions (Native TPM Cohorts Only)")
        plt.savefig(output_dir / "umap_native_tpm_only.svg", bbox_inches='tight')
        plt.close()
        print("Saved UMAP of Native TPM cohorts only.")

    # 4. PREDICTIVE DRIVER CHECK: Feature Importances
    print("\nPhase 4: Evaluating Prediction Drivers...")
    feature_importances = {}
    
    for cohort in cohorts:
        df_clin = load_clinical(cohort, lair)
        if df_clin is None:
            continue
            
        bagaev_deconv_key = f"deconv_bagaev_{cohort}_{ref_name}.csv"
        deconv_path = bagaev_deconv_paths.get(bagaev_deconv_key)
        if not deconv_path:
            continue
            
        df_fracs = pd.read_csv(deconv_path, index_col=0)
        df_fracs.index = df_fracs.index.astype(str)
        
        # Construct lineage-level features using cluster_to_major mapping
        df_lin = pd.DataFrame(index=df_fracs.index)
        for lin in lineages:
            matching_cols = [c for c, maj in cluster_to_major.items() if maj == lin and c in df_fracs.columns]
            df_lin[lin] = df_fracs[matching_cols].sum(axis=1)
            
        common_samples = df_lin.index.intersection(df_clin.index)
        if len(common_samples) < 10:
            continue
            
        X = df_lin.loc[common_samples, lineages].values
        y = df_clin.loc[common_samples, 'response'].map({'R': 1, 'NR': 0}).values
        
        if len(np.unique(y)) < 2:
            continue
            
        # Fit Random Forest over 5 seeds and average importances
        importances_seeds = []
        for seed in [42, 101, 2023, 7, 888]:
            skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
            min_class = min(np.bincount(y))
            if min_class < 5:
                from sklearn.model_selection import KFold
                skf = KFold(n_splits=5, shuffle=True, random_state=seed)
                
            for train_idx, _ in skf.split(X, y):
                model = RandomForestClassifier(n_estimators=50, max_depth=4, random_state=seed, n_jobs=-1)
                model.fit(X[train_idx], y[train_idx])
                importances_seeds.append(model.feature_importances_)
                
        mean_importances = np.mean(importances_seeds, axis=0)
        feature_importances[cohort] = mean_importances
        
    df_importances = pd.DataFrame(feature_importances, index=lineages).T
    # Reorder rows to group custom transformed cohorts together
    df_importances = df_importances.reindex(sorted(df_importances.index, key=lambda x: COHORT_TYPES[x]))
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(df_importances, annot=True, cmap='YlGnBu', fmt=".3f", cbar_kws={'label': 'Mean Gini Importance'})
    plt.title("Importance of Major Cell Lineages in Response Prediction")
    plt.ylabel("Cohort")
    plt.savefig(output_dir / "predictive_features_comparison.svg", bbox_inches='tight')
    plt.close()
    print("Saved predictive features comparison heatmap.")
    print("Done!")


if __name__ == "__main__":
    main()
