#!/usr/bin/env python3
"""
Advanced Immunotherapy Response Prediction using Random Forest.
Compares different feature engineering strategies across multiple random seeds
and plots the average performance with standard deviation error bars.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import (
    roc_auc_score, average_precision_score, accuracy_score,
    precision_score, recall_score, f1_score, matthews_corrcoef
)

# Ensure local packages are in path
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))
import datalair
from single_cell_datasets._single_cell_reference import SingleCellDeconvolution


def load_clinical_and_tmb(cohort_name: str, lair):
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
        return df_clinical
    except Exception as e:
        print(f"Error loading clinical for {cohort_name}: {e}")
        return None


def run_cross_val_rf_advanced(X, y, feature_eng_method, seed):
    """Runs 5-fold CV Random Forest with training-set feature engineering to prevent target leakage."""
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    oof_preds = np.zeros(len(y))
    
    for train_idx, val_idx in cv.split(X, y):
        X_train, X_val = X[train_idx], X[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]
        
        # Apply Feature Engineering
        if feature_eng_method == 'raw':
            pass
        elif feature_eng_method == 'filter':
            # Drop features with mean abundance < 0.005 (0.5%)
            mean_vals = X_train.mean(axis=0)
            valid_cols = mean_vals >= 0.005
            if not np.any(valid_cols):
                valid_cols = np.ones(X_train.shape[1], dtype=bool)
            X_train = X_train[:, valid_cols]
            X_val = X_val[:, valid_cols]
            
        elif feature_eng_method == 'select':
            # Univariate ANOVA feature selection (k=10)
            k = min(10, X_train.shape[1])
            selector = SelectKBest(score_func=f_classif, k=k)
            X_train = selector.fit_transform(X_train, y_train)
            X_val = selector.transform(X_val)
            
        elif feature_eng_method == 'pca':
            # PCA dimensionality reduction (5 components)
            n_comp = min(5, X_train.shape[1])
            pca = PCA(n_components=n_comp, random_state=seed)
            X_train = pca.fit_transform(X_train)
            X_val = pca.transform(X_val)
            
        elif feature_eng_method == 'filter_select':
            # Filter low-abundance features first, then select top 10
            mean_vals = X_train.mean(axis=0)
            valid_cols = mean_vals >= 0.005
            if not np.any(valid_cols):
                valid_cols = np.ones(X_train.shape[1], dtype=bool)
            X_train = X_train[:, valid_cols]
            X_val = X_val[:, valid_cols]
            
            k = min(10, X_train.shape[1])
            selector = SelectKBest(score_func=f_classif, k=k)
            X_train = selector.fit_transform(X_train, y_train)
            X_val = selector.transform(X_val)
            
        elif feature_eng_method == 'filter_pca':
            # Filter low-abundance features first, then PCA to 5 components
            mean_vals = X_train.mean(axis=0)
            valid_cols = mean_vals >= 0.005
            if not np.any(valid_cols):
                valid_cols = np.ones(X_train.shape[1], dtype=bool)
            X_train = X_train[:, valid_cols]
            X_val = X_val[:, valid_cols]
            
            n_comp = min(5, X_train.shape[1])
            pca = PCA(n_components=n_comp, random_state=seed)
            X_train = pca.fit_transform(X_train)
            X_val = pca.transform(X_val)
            
        # Fit Random Forest
        clf = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=seed)
        clf.fit(X_train, y_train)
        oof_preds[val_idx] = clf.predict_proba(X_val)[:, 1]
        
    return oof_preds


def main():
    lair_path = "/storage/halu/lair"
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output" / "single-cell-exploration"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get deconvolution dataset filepaths
    ds_deconv = SingleCellDeconvolution()
    deconv_paths = lair.get_dataset_filepaths(ds_deconv)
    
    cohorts = [
        "Hugo-iAtlas",
        "Riaz-iAtlas",
        "Liu-iAtlas",
        "Gide-iAtlas",
        "Rosenberg-iAtlas",
        "Padron-iAtlas",
        "Anders-iAtlas",
        "McDermott-iAtlas",
        "Choueiri-iAtlas",
        "Combined-Melanoma"
    ]
    
    resolutions = ["leiden_res_0.5", "kmeans_subcluster_res_0.5"]
    methods = ["raw", "filter", "select", "pca", "filter_select", "filter_pca"]
    seeds = [42, 101, 2023, 7, 888]
    
    # Pre-load clinical data
    clinical_data = {}
    for cohort in cohorts:
        if cohort != 'Combined-Melanoma':
            df_clin = load_clinical_and_tmb(cohort, lair)
            if df_clin is not None:
                clinical_data[cohort] = df_clin
                
    results_rows = []
    
    # Loop over resolutions, seeds, cohorts, and methods
    for resolution in resolutions:
        print(f"\n==========================================")
        print(f"Running experiments for {resolution}")
        print(f"==========================================")
        
        for seed in seeds:
            print(f"  Running random seed {seed}...")
            
            for cohort in cohorts:
                # Load aligned X and y
                if cohort == 'Combined-Melanoma':
                    mel_features_list = []
                    mel_targets_list = []
                    for c in ['Hugo-iAtlas', 'Riaz-iAtlas', 'Liu-iAtlas', 'Gide-iAtlas']:
                        key = f"deconv_{c}_{resolution}.csv"
                        path = deconv_paths.get(key)
                        df_clin = clinical_data.get(c)
                        if path and df_clin is not None:
                            df_fracs = pd.read_csv(path, index_col=0)
                            df_fracs.index = df_fracs.index.astype(str)
                            common_samples = df_fracs.index.intersection(df_clin.index)
                            if len(common_samples) >= 10:
                                mel_features_list.append(df_fracs.loc[common_samples])
                                mel_targets_list.append(df_clin.loc[common_samples, 'response'].map({'R': 1, 'NR': 0}))
                    if not mel_features_list:
                        continue
                    X_df = pd.concat(mel_features_list, axis=0)
                    y_series = pd.concat(mel_targets_list, axis=0)
                else:
                    key = f"deconv_{cohort}_{resolution}.csv"
                    path = deconv_paths.get(key)
                    df_clin = clinical_data.get(cohort)
                    if not path or df_clin is None:
                        continue
                    df_fracs = pd.read_csv(path, index_col=0)
                    df_fracs.index = df_fracs.index.astype(str)
                    common_samples = df_fracs.index.intersection(df_clin.index)
                    if len(common_samples) < 10:
                        continue
                    X_df = df_fracs.loc[common_samples]
                    y_series = df_clin.loc[common_samples, 'response'].map({'R': 1, 'NR': 0})
                    
                X = X_df.values
                y = y_series.values
                p = y.mean()
                
                # Append Random Baseline (same for all seeds)
                results_rows.append({
                    'Resolution': resolution,
                    'Cohort': cohort,
                    'Method': 'random_baseline',
                    'ROC_AUC': 0.5,
                    'PR_AUC': p,
                    'Accuracy': 0.5,
                    'MCC': 0.0,
                    'Seed': seed,
                    'Num_Samples': len(y)
                })
                
                # Loop over methods
                for method in methods:
                    oof_preds = run_cross_val_rf_advanced(X, y, method, seed)
                    
                    # Compute probability metrics
                    roc_auc = roc_auc_score(y, oof_preds)
                    pr_auc = average_precision_score(y, oof_preds)
                    
                    # Compute classification metrics (threshold 0.5)
                    y_pred = (oof_preds >= 0.5).astype(int)
                    acc = accuracy_score(y, y_pred)
                    mcc = matthews_corrcoef(y, y_pred)
                    
                    results_rows.append({
                        'Resolution': resolution,
                        'Cohort': cohort,
                        'Method': method,
                        'ROC_AUC': roc_auc,
                        'PR_AUC': pr_auc,
                        'Accuracy': acc,
                        'MCC': mcc,
                        'Seed': seed,
                        'Num_Samples': len(y)
                    })
                    
    # Save statistics table
    df_results = pd.DataFrame(results_rows)
    df_results.to_csv(output_dir / "deconv_prediction_advanced_metrics.csv", index=False)
    print(f"\nSaved advanced prediction metrics to deconv_prediction_advanced_metrics.csv")
    
    # --- Plotting ---
    # We will generate separate plots for resolution: leiden_res_0.5 and kmeans_subcluster_res_0.5
    for resolution in resolutions:
        df_res = df_results[df_results['Resolution'] == resolution]
        
        # Re-order cohorts so Combined-Melanoma is last
        cohort_order = [c for c in cohorts if c != 'Combined-Melanoma'] + ['Combined-Melanoma']
        
        fig, axes = plt.subplots(2, 2, figsize=(18, 14), dpi=300)
        axes_flat = axes.flatten()
        
        metrics = [
            ('ROC_AUC', 'ROC AUC', 0.3, 0.95),
            ('PR_AUC', 'PR AUC', 0.0, 1.0),
            ('Accuracy', 'Accuracy', 0.3, 0.95),
            ('MCC', 'Matthews Correlation Coefficient (MCC)', -0.3, 0.8)
        ]
        
        palette = {
            'raw': '#3498DB',              # Blue
            'filter': '#E67E22',           # Orange
            'select': '#2ECC71',           # Green
            'pca': '#9B59B6',              # Purple
            'filter_select': '#E74C3C',    # Red
            'filter_pca': '#F1C40F',       # Yellow
            'random_baseline': '#7F8C8D'   # Gray
        }
        
        method_labels = {
            'raw': 'Raw (Baseline)',
            'filter': 'Filtered (Mean > 0.5%)',
            'select': 'SelectKBest (ANOVA k=10)',
            'pca': 'PCA (n_components=5)',
            'filter_select': 'Filtered + ANOVA',
            'filter_pca': 'Filtered + PCA',
            'random_baseline': 'Random Baseline'
        }
        
        for idx, (col_name, title, ymin, ymax) in enumerate(metrics):
            ax = axes_flat[idx]
            
            # Map method names to pretty labels
            df_plot = df_res.copy()
            df_plot['Method_Label'] = df_plot['Method'].map(method_labels)
            
            # Grouped barplot with errorbar='sd' (standard deviation)
            sns.barplot(
                data=df_plot, x='Cohort', y=col_name, hue='Method', ax=ax,
                order=cohort_order,
                palette=palette,
                edgecolor='black', linewidth=0.5,
                errorbar='sd', capsize=0.08, err_kws={'linewidth': 0.8}
            )
            
            if col_name in ['ROC_AUC', 'Accuracy']:
                ax.axhline(0.5, color='#E74C3C', linestyle=':', alpha=0.8)
            elif col_name == 'MCC':
                ax.axhline(0.0, color='#E74C3C', linestyle=':', alpha=0.8)
                
            ax.set_title(title + " (Mean ± SD)", fontsize=12, fontweight='bold')
            ax.set_xlabel("Trial Cohort", fontsize=9)
            ax.set_ylabel(title, fontsize=9)
            ax.set_ylim(ymin, ymax)
            ax.grid(True, linestyle='--', alpha=0.3, axis='y')
            ax.tick_params(axis='x', rotation=15)
            ax.legend().remove()
            
        # Draw a unified legend on the right side
        handles, labels = axes_flat[0].get_legend_handles_labels()
        pretty_labels = [method_labels.get(l, l) for l in labels]
        fig.legend(handles, pretty_labels, loc='upper right', bbox_to_anchor=(0.98, 0.98), title="Feature Engineering Strategy", fontsize=10)
        
        plt.suptitle(f"Feature Engineering & Dimensionality Reduction Comparison ({resolution})\n(Averaged across 5 Random Seeds)", 
                     fontsize=15, fontweight='bold', y=0.99)
        plt.tight_layout()
        plt.savefig(output_dir / f"deconv_prediction_advanced_grid_{resolution}.svg", format='svg', bbox_inches='tight')
        plt.close()
        print(f"Saved advanced metrics plot for {resolution}.")


if __name__ == "__main__":
    main()
