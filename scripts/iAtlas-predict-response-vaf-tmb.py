#!/usr/bin/env python3
"""
Predict response to immunotherapy using TMB scores computed at different VAF thresholds as features.
Evaluates 6 classifiers (Random Forest, Adaline, MLP-4 ReLU, MLP-4 Sigmoid, Linear Model, Logistic Regression)
across 13 feature configurations inside a Stratified 5-Fold CV loop.
Generates model fit diagnostics (over/underfitting) and ROC curves.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.exceptions import ConvergenceWarning

warnings.filterwarnings('ignore', category=ConvergenceWarning)

# Ensure local packages are in python path if running directly
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))

import datalair
import ici_datasets
from gene_utils import calculate_maf_tmb


def oversample_to_balance(X_tr, y_tr):
    """
    Oversample the minority class in the training split to balance classes.
    """
    y_tr = pd.Series(y_tr)
    c0 = X_tr[y_tr == 0]
    c1 = X_tr[y_tr == 1]
    if len(c0) == len(c1):
        return X_tr, y_tr
    if len(c0) > len(c1):
        maj_df, min_df = c0, c1
        maj_val, min_val = 0, 1
    else:
        maj_df, min_df = c1, c0
        maj_val, min_val = 1, 0
        
    rng = np.random.default_rng(42)
    min_indices = rng.choice(min_df.index, size=len(maj_df), replace=True)
    X_resampled = pd.concat([maj_df, X_tr.loc[min_indices]], axis=0)
    y_resampled = pd.concat([pd.Series(maj_val, index=maj_df.index), pd.Series(min_val, index=min_indices)], axis=0)
    return X_resampled, y_resampled


def get_cohort_vaf_tmb_data(lair, ds_names, vaf_thresholds):
    """
    Load mutations and clinical responses for the given cohorts and compute TMB for different VAF thresholds.
    """
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    all_X = []
    all_y = []
    
    for ds_name in ds_names:
        try:
            data_dir = ici_datasets.cbioportal_datasets.get_dataset_dir(lair, dataset_class, ds_name)
            df_clinical = ici_datasets.cbioportal_datasets.load_data_clinical(data_dir)
            df_clinical = df_clinical[df_clinical['response'].isin(['R', 'NR'])]
            if len(df_clinical) == 0:
                continue
                
            mut_file = data_dir / "data_mutations.txt"
            if not mut_file.exists():
                continue
                
            df_mut = pd.read_csv(mut_file, sep="\t", low_memory=False)
            
            # Compute TMB for each VAF threshold
            tmb_features = pd.DataFrame(index=df_clinical.index)
            for thresh in vaf_thresholds:
                tmb_matrix = calculate_maf_tmb(df_mut, capture_size_mb=30.0, vaf_threshold=thresh)
                tmb_matrix = tmb_matrix.set_index('Tumor_Sample_Barcode')
                tmb_features[f'TMB_VAF_{thresh:.2f}'] = tmb_matrix['TMB_Score'].reindex(df_clinical.index, fill_value=0.0)
                
            tmb_features.index = f"{ds_name}_" + tmb_features.index.astype(str)
            y = df_clinical['response'].copy()
            y.index = f"{ds_name}_" + y.index.astype(str)
            
            all_X.append(tmb_features)
            all_y.append(y)
        except Exception as e:
            print(f"  Error loading dataset {ds_name}: {e}")
            
    if not all_X:
        return None, None
        
    combined_X = pd.concat(all_X, axis=0, sort=False).fillna(0.0)
    combined_y = pd.concat(all_y, axis=0)
    
    return combined_X, combined_y


def plot_overfitting_underfitting_grid(cohort_name, train_means, val_means, plot_path):
    """
    Generate a 2x3 grid of subplots showing Training vs. Validation AUC
    for each of the 6 classifiers across the 13 configurations.
    """
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), dpi=300)
    axes = axes.flatten()
    
    clf_keys = ['rf', 'adaline', 'mlp4_relu', 'mlp4_sigmoid', 'linear_model', 'logistic_regression']
    clf_titles = {
        'rf': 'Random Forest Classifier',
        'adaline': 'Adaline (Single-Layer NN)',
        'mlp4_relu': 'MLP-4 (ReLU Activation)',
        'mlp4_sigmoid': 'MLP-4 (Sigmoid Activation)',
        'linear_model': 'Linear Model (OLS)',
        'logistic_regression': 'Logistic Regression'
    }
    
    x_labels = [
        'VAF 0.01', 'VAF 0.02', 'VAF 0.05', 'VAF 0.10', 'VAF 0.15',
        'VAF 0.20', 'VAF 0.25', 'VAF 0.30', 'VAF 0.40', 'VAF 0.50',
        'Low VAF TMBs', 'High VAF TMBs', 'All VAF TMBs'
    ]
    
    x_indices = np.arange(len(x_labels))
    
    for idx, clf in enumerate(clf_keys):
        ax = axes[idx]
        
        tr_aucs = [train_means[(clf, x)] for x in x_labels]
        val_aucs = [val_means[(clf, x)] for x in x_labels]
        
        ax.plot(x_indices, tr_aucs, 'o-', label='Training AUC', color='#E74C3C', linewidth=1.8, markersize=4)
        ax.plot(x_indices, val_aucs, 'o-', label='Validation AUC', color='#3498DB', linewidth=1.8, markersize=4)
        
        # Shading overfitting gap
        ax.fill_between(x_indices, tr_aucs, val_aucs, color='#E74C3C', alpha=0.1, label='Overfitting Gap')
        
        ax.set_title(clf_titles[clf], fontsize=11, fontweight='bold', color='#2C3E50')
        ax.set_xlabel("VAF TMB Configurations", fontsize=9, fontweight='semibold')
        ax.set_ylabel("ROC AUC Score", fontsize=9, fontweight='semibold')
        ax.set_xticks(x_indices)
        ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=8)
        ax.set_ylim(0.2, 1.05)
        ax.legend(loc='lower left', fontsize=8)
        ax.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
        
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
            
    plt.suptitle(f"Model Fit Diagnostics (TMB at Multiple VAFs): {cohort_name}", 
                 fontsize=14, fontweight='bold', color='#2C3E50', y=0.99)
                 
    plt.tight_layout()
    plt.savefig(plot_path, format='svg', bbox_inches='tight')
    plt.close(fig)


def analyze_cohort_vaf_tmb(X, y, name, output_dir):
    """
    Run Stratified 5-Fold CV on TMB features calculated at multiple VAF thresholds.
    """
    cohort_dir = output_dir / name
    cohort_dir.mkdir(parents=True, exist_ok=True)
    
    n_samples = len(y)
    r_count = sum(y == 'R')
    nr_count = sum(y == 'NR')
    
    if n_samples < 10 or r_count < 3 or nr_count < 3:
        print(f"Skipping {name}: Too few samples/classes (n={n_samples}, R={r_count}, NR={nr_count}).")
        return None
        
    print(f"Analyzing {name} (samples = {n_samples})...")
    
    # Classifiers (raw feature values, no scaling)
    classifiers = {
        'rf': RandomForestClassifier(n_estimators=100, random_state=42),
        'adaline': MLPClassifier(hidden_layer_sizes=(), activation='identity', solver='adam', max_iter=1000, random_state=42),
        'mlp4_relu': MLPClassifier(hidden_layer_sizes=(4,), activation='relu', solver='adam', max_iter=1000, random_state=42),
        'mlp4_sigmoid': MLPClassifier(hidden_layer_sizes=(4,), activation='logistic', solver='adam', max_iter=1000, random_state=42),
        'linear_model': LinearRegression(),
        'logistic_regression': LogisticRegression(random_state=42)
    }
    
    y_encoded = (y == 'R').astype(int)
    cv = StratifiedKFold(n_splits=min(5, n_samples), shuffle=True, random_state=42)
    
    x_labels = [
        'VAF 0.01', 'VAF 0.02', 'VAF 0.05', 'VAF 0.10', 'VAF 0.15',
        'VAF 0.20', 'VAF 0.25', 'VAF 0.30', 'VAF 0.40', 'VAF 0.50',
        'Low VAF TMBs', 'High VAF TMBs', 'All VAF TMBs'
    ]
    
    # Define features mapping for each config
    config_features = {
        'VAF 0.01': ['TMB_VAF_0.01'],
        'VAF 0.02': ['TMB_VAF_0.02'],
        'VAF 0.05': ['TMB_VAF_0.05'],
        'VAF 0.10': ['TMB_VAF_0.10'],
        'VAF 0.15': ['TMB_VAF_0.15'],
        'VAF 0.20': ['TMB_VAF_0.20'],
        'VAF 0.25': ['TMB_VAF_0.25'],
        'VAF 0.30': ['TMB_VAF_0.30'],
        'VAF 0.40': ['TMB_VAF_0.40'],
        'VAF 0.50': ['TMB_VAF_0.50'],
        'Low VAF TMBs': ['TMB_VAF_0.01', 'TMB_VAF_0.02', 'TMB_VAF_0.05'],
        'High VAF TMBs': ['TMB_VAF_0.10', 'TMB_VAF_0.15', 'TMB_VAF_0.20', 'TMB_VAF_0.25', 'TMB_VAF_0.30', 'TMB_VAF_0.40', 'TMB_VAF_0.50'],
        'All VAF TMBs': [f'TMB_VAF_{t:.2f}' for t in [0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]]
    }
    
    # Pre-allocate dictionaries to store fold AUCs
    train_aucs = { (clf, cfg): [] for clf in classifiers for cfg in x_labels }
    val_aucs = { (clf, cfg): [] for clf in classifiers for cfg in x_labels }
    
    try:
        for train_idx, val_idx in cv.split(X, y_encoded):
            X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
            y_train = y_encoded.iloc[train_idx]
            y_val = y_encoded.iloc[val_idx]
            
            for cfg in x_labels:
                features = config_features[cfg]
                
                X_tr_sub = X_train[features]
                X_val_sub = X_val[features]
                
                # Oversample training fold to balance classes
                X_tr_res, y_tr_res = oversample_to_balance(X_tr_sub, y_train)
                
                for clf_key, clf_model in classifiers.items():
                    clf = clone(clf_model)
                    clf.fit(X_tr_res, y_tr_res)
                    
                    # Predict Training
                    if hasattr(clf, 'predict_proba'):
                        tr_prob = clf.predict_proba(X_tr_res)[:, 1]
                    else:
                        tr_prob = clf.predict(X_tr_res)
                    tr_auc = roc_auc_score(y_tr_res, tr_prob)
                    train_aucs[(clf_key, cfg)].append(tr_auc)
                    
                    # Predict Validation
                    if hasattr(clf, 'predict_proba'):
                        val_prob = clf.predict_proba(X_val_sub)[:, 1]
                    else:
                        val_prob = clf.predict(X_val_sub)
                    val_auc = roc_auc_score(y_val, val_prob)
                    val_aucs[(clf_key, cfg)].append(val_auc)
                    
        # Average AUCs across folds
        train_means = { (clf, cfg): np.mean(train_aucs[(clf, cfg)]) for clf in classifiers for cfg in x_labels }
        val_means = { (clf, cfg): np.mean(val_aucs[(clf, cfg)]) for clf in classifiers for cfg in x_labels }
        
        # Plot over/underfitting grid
        plot_path = cohort_dir / "overfitting_underfitting_vaf_tmb.svg"
        print(f"  -> Plotting overfitting diagnostic curves to {plot_path.name}...")
        plot_overfitting_underfitting_grid(name, train_means, val_means, plot_path)
        
        # Return validation performance dictionary for the table
        perf_metrics_list = []
        for clf_key in classifiers:
            perf_row = {
                'Cohort': name,
                'Classifier': clf_key
            }
            for cfg in x_labels:
                # Format header column name, e.g. "VAF 0.05" -> "AUC_VAF_0.05"
                col_name = cfg.lower().replace(' ', '_')
                perf_row[f"AUC_{col_name}"] = val_means[(clf_key, cfg)]
            perf_metrics_list.append(perf_row)
            
        return perf_metrics_list
        
    except Exception as e:
        print(f"  Warning: Cross-validation failed for {name}: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    lair_path = "/storage/halu/lair"
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    # Resolve directories
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output" / "iAtlas-vaf-tmb"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    vaf_thresholds = [0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]
    
    # Individual cohorts with mutation data
    individual_cohorts = [
        'Rosenberg-iAtlas',
        'Riaz-iAtlas',
        'Liu-iAtlas',
        'Hugo-iAtlas',
        'Padron-iAtlas',
        'Choueiri-iAtlas',
        'Anders-iAtlas',
        'McDermott-iAtlas'
    ]
    
    # Cancer-type combinations
    cancer_type_combinations = {
        'Combined-Melanoma': ['Riaz-iAtlas', 'Liu-iAtlas', 'Hugo-iAtlas'],
        'Combined-Kidney': ['Choueiri-iAtlas', 'McDermott-iAtlas']
    }
    
    performances_list = []
    
    # 1. Run Individual Cohort Analyses
    print("\n--- Running Individual Cohort Analyses ---")
    for cohort in individual_cohorts:
        X, y = get_cohort_vaf_tmb_data(lair, [cohort], vaf_thresholds)
        if X is not None and y is not None:
            perf_metrics = analyze_cohort_vaf_tmb(X, y, cohort, output_dir)
            if perf_metrics is not None:
                performances_list.extend(perf_metrics)
                
    # 2. Run Combined Cohort Analyses
    print("\n--- Running Combined Cohort Analyses ---")
    for comb_name, cohorts in cancer_type_combinations.items():
        X, y = get_cohort_vaf_tmb_data(lair, cohorts, vaf_thresholds)
        if X is not None and y is not None:
            perf_metrics = analyze_cohort_vaf_tmb(X, y, comb_name, output_dir)
            if perf_metrics is not None:
                performances_list.extend(perf_metrics)
                
    # 3. Save Summary Performances CSV
    if performances_list:
        df_perf = pd.DataFrame(performances_list)
        perf_csv_path = output_dir / "vaf_tmb_model_performances.csv"
        df_perf.to_csv(perf_csv_path, index=False)
        print(f"\nSaved VAF TMB model performances comparison to {perf_csv_path}")
        
        # Display subset of performances
        print("\nModel Performances Comparison (VAF TMB ROC AUC Subset):")
        subset_cols = ['Cohort', 'Classifier', 'AUC_vaf_0.05', 'AUC_low_vaf_tmbs', 'AUC_high_vaf_tmbs', 'AUC_all_vaf_tmbs']
        print(df_perf[subset_cols].to_string(index=False))
    else:
        print("\nNo cohorts could be analyzed.")


if __name__ == "__main__":
    main()
