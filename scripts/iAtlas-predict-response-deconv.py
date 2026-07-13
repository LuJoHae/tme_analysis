#!/usr/bin/env python3
"""
Predict response to immunotherapy from deconvolution cell fractions using Random Forest.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score

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


def run_cross_val_rf(X, y):
    """Runs 5-fold CV Random Forest and returns predicted probabilities and importances."""
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    oof_preds = np.zeros(len(y))
    importances = np.zeros(X.shape[1])
    
    for train_idx, val_idx in cv.split(X, y):
        X_train, X_val = X[train_idx], X[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]
        
        clf = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42, n_jobs=-1)
        clf.fit(X_train, y_train)
        
        oof_preds[val_idx] = clf.predict_proba(X_val)[:, 1]
        importances += clf.feature_importances_ / 5.0
        
    return oof_preds, importances


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
        "Choueiri-iAtlas"
    ]
    
    resolutions = ["leiden_res_0.5", "kmeans_subcluster_res_0.5"]
    
    results_rows = []
    importance_dfs = []
    
    # Pre-load clinical data
    clinical_data = {}
    for cohort in cohorts:
        df_clin = load_clinical_and_tmb(cohort, lair)
        if df_clin is not None:
            clinical_data[cohort] = df_clin
            
    for resolution in resolutions:
        print(f"\n--- Running Random Forest predictions for {resolution} ---")
        
        melanoma_features_list = []
        melanoma_targets_list = []
        
        for cohort in cohorts:
            key = f"deconv_{cohort}_{resolution}.csv"
            path = deconv_paths.get(key)
            df_clin = clinical_data.get(cohort)
            
            if not path or df_clin is None:
                continue
                
            # Load deconvolution fractions
            df_fracs = pd.read_csv(path, index_col=0)
            df_fracs.index = df_fracs.index.astype(str)
            
            # Align
            common_samples = df_fracs.index.intersection(df_clin.index)
            if len(common_samples) < 10:
                print(f"Skipping {cohort}: Too few samples ({len(common_samples)})")
                continue
                
            X_df = df_fracs.loc[common_samples]
            y_series = df_clin.loc[common_samples, 'response'].map({'R': 1, 'NR': 0})
            
            X = X_df.values
            y = y_series.values
            
            oof_preds, importances = run_cross_val_rf(X, y)
            
            roc_auc = roc_auc_score(y, oof_preds)
            pr_auc = average_precision_score(y, oof_preds)
            
            print(f"  {cohort}: ROC AUC = {roc_auc:.4f}, PR AUC = {pr_auc:.4f} (n={len(y)})")
            
            results_rows.append({
                'Resolution': resolution,
                'Cohort': cohort,
                'ROC_AUC': roc_auc,
                'PR_AUC': pr_auc,
                'Num_Samples': len(y)
            })
            
            # Save feature importances
            imp_df = pd.DataFrame({
                'Feature': X_df.columns,
                'Importance': importances,
                'Cohort': cohort,
                'Resolution': resolution
            })
            importance_dfs.append(imp_df)
            
            # Collect for Combined Melanoma
            if cohort in ['Hugo-iAtlas', 'Riaz-iAtlas', 'Liu-iAtlas', 'Gide-iAtlas']:
                melanoma_features_list.append(X_df)
                melanoma_targets_list.append(y_series)
                
        # Run Combined Melanoma
        if melanoma_features_list:
            X_mel = pd.concat(melanoma_features_list, axis=0)
            y_mel = pd.concat(melanoma_targets_list, axis=0)
            
            oof_preds, importances = run_cross_val_rf(X_mel.values, y_mel.values)
            roc_auc = roc_auc_score(y_mel.values, oof_preds)
            pr_auc = average_precision_score(y_mel.values, oof_preds)
            
            print(f"  Combined-Melanoma: ROC AUC = {roc_auc:.4f}, PR AUC = {pr_auc:.4f} (n={len(y_mel)})")
            
            results_rows.append({
                'Resolution': resolution,
                'Cohort': 'Combined-Melanoma',
                'ROC_AUC': roc_auc,
                'PR_AUC': pr_auc,
                'Num_Samples': len(y_mel)
            })
            
            imp_df = pd.DataFrame({
                'Feature': X_mel.columns,
                'Importance': importances,
                'Cohort': 'Combined-Melanoma',
                'Resolution': resolution
            })
            importance_dfs.append(imp_df)
            
    # Save performance metrics table
    df_perf = pd.DataFrame(results_rows)
    df_perf.to_csv(output_dir / "deconv_prediction_performance.csv", index=False)
    print("\nSaved prediction performance table.")
    
    # Save all importances
    df_imp_all = pd.concat(importance_dfs, axis=0)
    df_imp_all.to_csv(output_dir / "deconv_feature_importances.csv", index=False)
    
    # --- Plotting ---
    # Plot 1: Performance comparison bar plot
    plt.figure(figsize=(12, 6), dpi=300)
    # Filter out Combined-Melanoma for trial cohort comparison
    df_trials = df_perf[df_perf['Cohort'] != 'Combined-Melanoma']
    
    sns.barplot(
        data=df_trials, x='Cohort', y='ROC_AUC', hue='Resolution',
        palette={'leiden_res_0.5': '#3498DB', 'kmeans_subcluster_res_0.5': '#E67E22'},
        edgecolor='black', linewidth=0.8
    )
    plt.axhline(0.5, color='#7F8C8D', linestyle='--', label='Chance (AUC=0.5)')
    plt.title("Immunotherapy Response Prediction (ROC AUC) using Deconvolution Fractions", fontsize=12, fontweight='bold', pad=15)
    plt.xlabel("Trial Cohort", fontsize=10, fontweight='bold')
    plt.ylabel("Cross-Validated ROC AUC", fontsize=10, fontweight='bold')
    plt.ylim(0.3, 0.8)
    plt.legend(loc='upper right', title="Reference Clustering")
    plt.grid(True, linestyle='--', alpha=0.3, axis='y')
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir / "deconv_prediction_performance.svg", format='svg', bbox_inches='tight')
    plt.close()
    
    # Plot 2: Feature Importances for Combined-Melanoma & Rosenberg-iAtlas
    fig, axes = plt.subplots(2, 2, figsize=(16, 12), dpi=300)
    
    # Leiden 0.5 Combined-Melanoma
    df_sub1 = df_imp_all[(df_imp_all['Cohort'] == 'Combined-Melanoma') & (df_imp_all['Resolution'] == 'leiden_res_0.5')]
    df_sub1 = df_sub1.sort_values(by='Importance', ascending=False).head(15)
    sns.barplot(data=df_sub1, x='Importance', y='Feature', ax=axes[0, 0], palette='viridis', hue='Feature', legend=False, edgecolor='black', linewidth=0.5)
    axes[0, 0].set_title("Combined Melanoma Top Features (Leiden 0.5)", fontsize=11, fontweight='bold')
    axes[0, 0].set_xlabel("Mean Feature Importance")
    axes[0, 0].set_ylabel("Cluster")
    
    # K-means 0.5 Combined-Melanoma
    df_sub2 = df_imp_all[(df_imp_all['Cohort'] == 'Combined-Melanoma') & (df_imp_all['Resolution'] == 'kmeans_subcluster_res_0.5')]
    df_sub2 = df_sub2.sort_values(by='Importance', ascending=False).head(15)
    sns.barplot(data=df_sub2, x='Importance', y='Feature', ax=axes[0, 1], palette='plasma', hue='Feature', legend=False, edgecolor='black', linewidth=0.5)
    axes[0, 1].set_title("Combined Melanoma Top Features (K-Means 0.5)", fontsize=11, fontweight='bold')
    axes[0, 1].set_xlabel("Mean Feature Importance")
    axes[0, 1].set_ylabel("")
    
    # Leiden 0.5 Rosenberg-iAtlas
    df_sub3 = df_imp_all[(df_imp_all['Cohort'] == 'Rosenberg-iAtlas') & (df_imp_all['Resolution'] == 'leiden_res_0.5')]
    df_sub3 = df_sub3.sort_values(by='Importance', ascending=False).head(15)
    sns.barplot(data=df_sub3, x='Importance', y='Feature', ax=axes[1, 0], palette='inferno', hue='Feature', legend=False, edgecolor='black', linewidth=0.5)
    axes[1, 0].set_title("Rosenberg Bladder Top Features (Leiden 0.5)", fontsize=11, fontweight='bold')
    axes[1, 0].set_xlabel("Mean Feature Importance")
    axes[1, 0].set_ylabel("Cluster")
    
    # K-means 0.5 Rosenberg-iAtlas
    df_sub4 = df_imp_all[(df_imp_all['Cohort'] == 'Rosenberg-iAtlas') & (df_imp_all['Resolution'] == 'kmeans_subcluster_res_0.5')]
    df_sub4 = df_sub4.sort_values(by='Importance', ascending=False).head(15)
    sns.barplot(data=df_sub4, x='Importance', y='Feature', ax=axes[1, 1], palette='magma', hue='Feature', legend=False, edgecolor='black', linewidth=0.5)
    axes[1, 1].set_title("Rosenberg Bladder Top Features (K-Means 0.5)", fontsize=11, fontweight='bold')
    axes[1, 1].set_xlabel("Mean Feature Importance")
    axes[1, 1].set_ylabel("")
    
    plt.suptitle("Top 15 Predictive Cell Type Features for Immunotherapy Response", fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(output_dir / "deconv_feature_importance_top.svg", format='svg', bbox_inches='tight')
    plt.close()
    
    print("All deconvolution prediction plots successfully generated!")


if __name__ == "__main__":
    main()
