#!/usr/bin/env python3
"""
Infer the best gene mutations to predict immunotherapeutic response across iAtlas cohorts.
Performs univariate Fisher's Exact Tests and multivariate classifiers (Random Forest, Adaline, MLP-4 ReLU, MLP-4 Sigmoid)
using multiple feature selection strategies (all, top 5, top 10, top 20, top 100, top 1000)
both with and without TMB as a feature, evaluated strictly inside the CV loop.
Saves diagnostic plots for over/underfitting and separate ROC curve grids for each classifier.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from scipy.stats import fisher_exact
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.exceptions import ConvergenceWarning

warnings.filterwarnings('ignore', category=ConvergenceWarning)

# Ensure local packages are in python path if running directly
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))

import datalair
import ici_datasets


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


def plot_roc_grid(name, rocs_data, n_samples, r_count, nr_count, plot_path, clf_title, clf_key):
    """
    Generate a 4x3 grid of subplots comparing the ROC curves of the 11 strategies.
    Each subplot shows: Mutations Only, Mutations + TMB, and the TMB-only baseline.
    """
    fig, axes = plt.subplots(4, 3, figsize=(15, 18), dpi=300)
    axes = axes.flatten()
    
    strategies = [
        'all',
        'mut_5', 'mut_10', 'mut_20', 'mut_100', 'mut_1000',
        'fish_5', 'fish_10', 'fish_20', 'fish_100', 'fish_1000'
    ]
    
    labels_map = {
        'all': 'All Qualifying Genes',
        'mut_5': 'Top 5 Mutated',
        'mut_10': 'Top 10 Mutated',
        'mut_20': 'Top 20 Mutated',
        'mut_100': 'Top 100 Mutated',
        'mut_1000': 'Top 1000 Mutated',
        'fish_5': 'Top 5 Fisher',
        'fish_10': 'Top 10 Fisher',
        'fish_20': 'Top 20 Fisher',
        'fish_100': 'Top 100 Fisher',
        'fish_1000': 'Top 1000 Fisher'
    }
    
    fpr_tmb, tpr_tmb, auc_tmb = rocs_data[(clf_key, 'tmb_only')]
    
    for idx, strat in enumerate(strategies):
        ax = axes[idx]
        
        # 1. TMB Only (baseline)
        if fpr_tmb is not None:
            ax.plot(fpr_tmb, tpr_tmb, label=f"TMB Only ({auc_tmb:.3f})", color='#7F8C8D', linestyle='--', linewidth=1.2)
            
        # 2. Mutations Only (without TMB)
        fpr_no, tpr_no, auc_no = rocs_data[(clf_key, strat, False)]
        if fpr_no is not None:
            ax.plot(fpr_no, tpr_no, label=f"Mutations Only ({auc_no:.3f})", color='#E67E22', linewidth=1.8)
            
        # 3. Mutations + TMB (with TMB)
        fpr_yes, tpr_yes, auc_yes = rocs_data[(clf_key, strat, True)]
        if fpr_yes is not None:
            ax.plot(fpr_yes, tpr_yes, label=f"Mutations + TMB ({auc_yes:.3f})", color='#2980B9', linewidth=1.8)
            
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)
        ax.set_title(labels_map[strat], fontsize=11, fontweight='bold', color='#2C3E50')
        ax.set_xlabel("False Positive Rate", fontsize=8)
        ax.set_ylabel("True Positive Rate", fontsize=8)
        ax.legend(loc="lower right", fontsize=7)
        ax.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
        ax.tick_params(labelsize=8)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
            
    # Hide the 12th subplot (empty)
    axes[11].axis('off')
    
    plt.suptitle(f"ROC Curves Comparison ({clf_title}): {name}\n(n = {n_samples}, R={r_count}, NR={nr_count})", 
                 fontsize=15, fontweight='bold', color='#2C3E50', y=0.995)
                 
    plt.tight_layout()
    plt.savefig(plot_path, format='svg', bbox_inches='tight')
    plt.close(fig)


def plot_overfitting_underfitting(cohort_name, train_means, val_means, plot_path):
    """
    Generate a 4-panel figure (1 row, 4 columns) showing Training AUC vs. Validation AUC
    for each classifier (Random Forest, Adaline, MLP-4 ReLU, MLP-4 Sigmoid).
    """
    fig, axes = plt.subplots(1, 4, figsize=(24, 5.5), dpi=300)
    
    clf_keys = ['rf', 'adaline', 'mlp4_relu', 'mlp4_sigmoid']
    clf_titles = {
        'rf': 'Random Forest Classifier',
        'adaline': 'Adaline (Single-Layer NN)',
        'mlp4_relu': 'MLP-4 (ReLU Activation)',
        'mlp4_sigmoid': 'MLP-4 (Sigmoid Activation)'
    }
    
    x_labels = [
        'TMB Only',
        'All', 'All+TMB',
        'Mut 5', 'Mut 5+TMB',
        'Mut 10', 'Mut 10+TMB',
        'Mut 20', 'Mut 20+TMB',
        'Mut 100', 'Mut 100+TMB',
        'Mut 1000', 'Mut 1000+TMB',
        'Fish 5', 'Fish 5+TMB',
        'Fish 10', 'Fish 10+TMB',
        'Fish 20', 'Fish 20+TMB',
        'Fish 100', 'Fish 100+TMB',
        'Fish 1000', 'Fish 1000+TMB'
    ]
    
    x_indices = np.arange(len(x_labels))
    
    for idx, clf in enumerate(clf_keys):
        ax = axes[idx]
        
        tr_aucs = [train_means[(clf, x)] for x in x_labels]
        val_aucs = [val_means[(clf, x)] for x in x_labels]
        
        ax.plot(x_indices, tr_aucs, 'o-', label='Training AUC', color='#E74C3C', linewidth=1.8, markersize=4)
        ax.plot(x_indices, val_aucs, 'o-', label='Validation AUC', color='#3498DB', linewidth=1.8, markersize=4)
        
        # Shading the overfitting gap
        ax.fill_between(x_indices, tr_aucs, val_aucs, color='#E74C3C', alpha=0.1, label='Overfitting Gap')
        
        ax.set_title(clf_titles[clf], fontsize=11, fontweight='bold', color='#2C3E50')
        ax.set_xlabel("Model Configurations", fontsize=9, fontweight='semibold')
        ax.set_ylabel("ROC AUC Score", fontsize=9, fontweight='semibold')
        ax.set_xticks(x_indices)
        ax.set_xticklabels(x_labels, rotation=90, fontsize=8)
        ax.set_ylim(0.2, 1.05)
        ax.legend(loc='lower left', fontsize=8)
        ax.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
        
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
            
    plt.suptitle(f"Model Fit Diagnostics (Over vs. Underfitting): {cohort_name}", 
                 fontsize=14, fontweight='bold', color='#2C3E50', y=0.985)
                 
    plt.tight_layout()
    plt.savefig(plot_path, format='svg', bbox_inches='tight')
    plt.close(fig)


def plot_strategy_importances(cohort_name, strat, name_label, importances_no, importances_yes, plot_path):
    """
    Generate side-by-side bar plots of RF feature importances with and without TMB.
    Highlights TMB_Score in red in the right panel.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5), dpi=300)
    
    # Left subplot: Without TMB
    ax_no = axes[0]
    df_no = pd.DataFrame(list(importances_no.items()), columns=['Feature', 'Importance'])
    df_no = df_no.sort_values(by='Importance', ascending=False).head(15)
    
    if len(df_no) > 0:
        df_no = df_no.iloc[::-1]
        colors = sns.color_palette("mako", len(df_no))
        bars = ax_no.barh(df_no['Feature'], df_no['Importance'], color=colors, height=0.6)
        ax_no.set_title("Model WITHOUT TMB", fontsize=11, fontweight='bold', color='#2C3E50')
        ax_no.set_xlabel("Average Feature Importance", fontsize=9, fontweight='semibold')
        ax_no.xaxis.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
        
        max_imp = df_no['Importance'].max()
        for bar in bars:
            width = bar.get_width()
            ax_no.text(width + (max_imp * 0.01), bar.get_y() + bar.get_height()/2, 
                       f"{width:.3f}", va='center', ha='left', fontsize=8, color='#34495E', fontweight='bold')
    else:
        ax_no.text(0.5, 0.5, "No features available", ha='center', va='center', fontsize=11, color='#7F8C8D')
    ax_no.tick_params(labelsize=9)
    for label in ax_no.get_yticklabels():
        label.set_fontweight('semibold')
        label.set_color('#34495E')
    for spine in ['top', 'right', 'left']:
        ax_no.spines[spine].set_visible(False)
        
    # Right subplot: With TMB
    ax_yes = axes[1]
    df_yes = pd.DataFrame(list(importances_yes.items()), columns=['Feature', 'Importance'])
    df_yes = df_yes.sort_values(by='Importance', ascending=False).head(15)
    
    if len(df_yes) > 0:
        df_yes = df_yes.iloc[::-1]
        colors = []
        for feat in df_yes['Feature']:
            if feat == 'TMB_Score':
                colors.append('#E74C3C')  # Red highlight for TMB
            else:
                colors.append('#2980B9')  # Blue for mutations
                
        bars = ax_yes.barh(df_yes['Feature'], df_yes['Importance'], color=colors, height=0.6)
        ax_yes.set_title("Model WITH TMB (TMB Highlighted)", fontsize=11, fontweight='bold', color='#2C3E50')
        ax_yes.set_xlabel("Average Feature Importance", fontsize=9, fontweight='semibold')
        ax_yes.xaxis.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
        
        max_imp = df_yes['Importance'].max()
        for bar in bars:
            width = bar.get_width()
            ax_yes.text(width + (max_imp * 0.01), bar.get_y() + bar.get_height()/2, 
                       f"{width:.3f}", va='center', ha='left', fontsize=8, color='#34495E', fontweight='bold')
    else:
        ax_yes.text(0.5, 0.5, "No features available", ha='center', va='center', fontsize=11, color='#7F8C8D')
    ax_yes.tick_params(labelsize=9)
    for label in ax_yes.get_yticklabels():
        label.set_fontweight('semibold')
        label.set_color('#34495E')
    for spine in ['top', 'right', 'left']:
        ax_yes.spines[spine].set_visible(False)
        
    plt.suptitle(f"{cohort_name} - Random Forest Feature Importances ({name_label})", 
                 fontsize=13, fontweight='bold', color='#2C3E50', y=0.985)
                 
    plt.tight_layout()
    plt.savefig(plot_path, format='svg', bbox_inches='tight')
    plt.close(fig)


def plot_strategy_fisher_pvalues(cohort_name, strat, name_label, fisher_pvals, plot_path):
    """
    Generate a bar chart of the average p-values from Fisher's exact tests in the CV loop.
    """
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    
    df_p = pd.DataFrame(list(fisher_pvals.items()), columns=['Gene', 'p_value'])
    df_p['Neg_Log_P'] = -np.log10(df_p['p_value'])
    df_p = df_p.sort_values(by='Neg_Log_P', ascending=False).head(15)
    
    if len(df_p) > 0:
        df_p = df_p.iloc[::-1]
        colors = sns.color_palette("flare", len(df_p))
        bars = ax.barh(df_p['Gene'], df_p['Neg_Log_P'], color=colors, height=0.6)
        
        ax.set_title(f"Fisher's Exact Test p-values in CV loop\n(Top 15 Genes for {name_label})", 
                     fontsize=11, fontweight='bold', color='#2C3E50', pad=15)
        ax.set_xlabel("-log10(Average p-value)", fontsize=9, fontweight='semibold')
        ax.xaxis.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
        
        max_nlp = df_p['Neg_Log_P'].max()
        for bar, (_, row) in zip(bars, df_p.iterrows()):
            width = bar.get_width()
            pval = row['p_value']
            ax.text(width + (max_nlp * 0.01), bar.get_y() + bar.get_height()/2, 
                    f"p={pval:.3g}", va='center', ha='left', fontsize=8, color='#34495E', fontweight='bold')
    else:
        ax.text(0.5, 0.5, "No genes available", ha='center', va='center', fontsize=11, color='#7F8C8D')
        
    ax.tick_params(labelsize=9)
    for label in ax.get_yticklabels():
        label.set_fontweight('semibold')
        label.set_color('#34495E')
    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)
        
    plt.tight_layout()
    plt.savefig(plot_path, format='svg', bbox_inches='tight')
    plt.close(fig)


def get_cohort_data(lair, ds_names):
    """
    Load mutations, clinical responses, and TMB for the given cohorts and merge them.
    Aligns mutated genes and pads wild-type entries with 0.
    """
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    
    all_X = []
    all_y = []
    all_tmb = []
    
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
            df_mut['VAF'] = np.where(df_mut['t_depth'] > 0, df_mut['t_alt_count'] / df_mut['t_depth'], 0)
            
            qualifying_classifications = {
                'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation',
                'Frame_Shift_Ins', 'Frame_Shift_Del', 'In_Frame_Ins',
                'In_Frame_Del', 'Translation_Start_Site', 'Splice_Site'
            }
            
            qualifying_subset = df_mut[
                (df_mut['Variant_Classification'].isin(qualifying_classifications)) &
                (df_mut['VAF'] >= 0.05)
            ].copy()
            
            qualifying_subset = qualifying_subset[qualifying_subset['Tumor_Sample_Barcode'].isin(df_clinical.index)]
            
            sample_gene_pairs = qualifying_subset[['Tumor_Sample_Barcode', 'Hugo_Symbol']].drop_duplicates()
            X_mut = pd.crosstab(sample_gene_pairs['Tumor_Sample_Barcode'], sample_gene_pairs['Hugo_Symbol'])
            X_mut = (X_mut > 0).astype(int)
            
            X_mut = X_mut.reindex(df_clinical.index, fill_value=0)
            X_mut.index = f"{ds_name}_" + X_mut.index.astype(str)
            
            y = df_clinical['response'].copy()
            y.index = f"{ds_name}_" + y.index.astype(str)
            
            tmb_score = df_clinical['TMB_Score'].copy()
            tmb_score.index = f"{ds_name}_" + tmb_score.index.astype(str)
            
            all_X.append(X_mut)
            all_y.append(y)
            all_tmb.append(tmb_score)
            
        except Exception as e:
            print(f"  Error loading dataset {ds_name}: {e}")
            
    if not all_X:
        return None, None, None
        
    combined_X = pd.concat(all_X, axis=0, sort=False).fillna(0).astype(int)
    combined_y = pd.concat(all_y, axis=0)
    combined_tmb = pd.concat(all_tmb, axis=0)
    
    return combined_X, combined_y, combined_tmb


def get_config_key(label):
    """
    Map x-axis labels to internal dictionary keys.
    """
    if label == 'TMB Only':
        return 'tmb_only'
    if '+' in label:
        base = label.split('+')[0]
        base_key = base.lower().replace(' ', '_')
        return (base_key, True)
    else:
        base_key = label.lower().replace(' ', '_')
        return (base_key, False)


def analyze_cohort_mutations(X, y, tmb, name, output_dir):
    """
    Perform statistical associations (Fisher's exact test) and predictive classification
    evaluating Random Forest, Adaline, MLP-4 ReLU, and MLP-4 Sigmoid in Stratified CV.
    Collects Train and Val AUCs across all configurations and generates overfitting plots.
    """
    # Create cohort specific directories
    cohort_dir = output_dir / name
    imp_dir = cohort_dir / "feature_importances"
    fish_dir = cohort_dir / "fisher_pvalues"
    
    imp_dir.mkdir(parents=True, exist_ok=True)
    fish_dir.mkdir(parents=True, exist_ok=True)

    gene_counts = X.sum(axis=0)
    selected_genes = gene_counts[gene_counts >= 3].index.tolist()
    
    if len(selected_genes) == 0:
        print(f"Skipping {name}: No genes mutated in >= 3 samples.")
        return None, None
        
    X_filtered = X[selected_genes]
    
    if len(np.unique(y)) < 2:
        print(f"Skipping {name}: Only one class present in response labels.")
        return None, None
        
    n_samples = len(y)
    r_count = sum(y == 'R')
    nr_count = sum(y == 'NR')
    if n_samples < 10 or r_count < 3 or nr_count < 3:
        print(f"Skipping {name}: Too few samples/classes (n={n_samples}, R={r_count}, NR={nr_count}).")
        return None, None
        
    print(f"Analyzing {name} (samples = {n_samples}, genes = {len(selected_genes)})...")
    
    # Univariate Fisher's Exact Test on full cohort
    fisher_results = []
    for gene in selected_genes:
        mut_status = X_filtered[gene]
        mut_R = sum((mut_status == 1) & (y == 'R'))
        wt_R = sum((mut_status == 0) & (y == 'R'))
        mut_NR = sum((mut_status == 1) & (y == 'NR'))
        mut_wt_NR = sum((mut_status == 0) & (y == 'NR'))
        
        odds_ratio, p_value = fisher_exact([[mut_R, mut_NR], [wt_R, mut_wt_NR]], alternative='two-sided')
        
        fisher_results.append({
            'Gene': gene,
            'Mutated_R': mut_R,
            'WT_R': wt_R,
            'Mutated_NR': mut_NR,
            'WT_NR': mut_wt_NR,
            'Odds_Ratio': odds_ratio,
            'p_value': p_value
        })
        
    df_fisher = pd.DataFrame(fisher_results)
    
    # Define Classifiers (No feature scaling will be applied)
    classifiers = {
        'rf': RandomForestClassifier(n_estimators=100, random_state=42),
        'adaline': MLPClassifier(hidden_layer_sizes=(), activation='identity', solver='adam', max_iter=1000, random_state=42),
        'mlp4_relu': MLPClassifier(hidden_layer_sizes=(4,), activation='relu', solver='adam', max_iter=1000, random_state=42),
        'mlp4_sigmoid': MLPClassifier(hidden_layer_sizes=(4,), activation='logistic', solver='adam', max_iter=1000, random_state=42)
    }
    
    y_encoded = (y == 'R').astype(int)
    cv = StratifiedKFold(n_splits=min(5, n_samples), shuffle=True, random_state=42)
    
    strategies = [
        'all',
        'mut_5', 'mut_10', 'mut_20', 'mut_100', 'mut_1000',
        'fish_5', 'fish_10', 'fish_20', 'fish_100', 'fish_1000'
    ]
    
    # Pre-allocate for validation probabilities (out-of-fold)
    y_probas_val = {}
    for clf_key in classifiers:
        for strat in strategies:
            for has_tmb in [False, True]:
                y_probas_val[(clf_key, strat, has_tmb)] = np.zeros(n_samples)
        y_probas_val[(clf_key, 'tmb_only')] = np.zeros(n_samples)
        
    # Pre-allocate lists of fold AUCs
    train_aucs = {}
    val_aucs = {}
    for clf_key in classifiers:
        for strat in strategies:
            for has_tmb in [False, True]:
                train_aucs[(clf_key, strat, has_tmb)] = []
                val_aucs[(clf_key, strat, has_tmb)] = []
        train_aucs[(clf_key, 'tmb_only')] = []
        val_aucs[(clf_key, 'tmb_only')] = []
        
    # Pre-allocate for feature importances and Fisher p-values (RF specific or model independent)
    cv_importances = { (strat, has_tmb): {} for strat in strategies for has_tmb in [False, True] }
    cv_fisher_pvals = { strat: {} for strat in strategies }
    
    try:
        for train_idx, val_idx in cv.split(X_filtered, y_encoded):
            X_train, X_val = X_filtered.iloc[train_idx], X_filtered.iloc[val_idx]
            y_train = y_encoded.iloc[train_idx]
            y_val = y_encoded.iloc[val_idx]
            tmb_train, tmb_val = tmb.iloc[train_idx], tmb.iloc[val_idx]
            
            # 1. Baseline TMB Only Model
            X_tr_tmb = tmb_train.to_frame()
            X_val_tmb = tmb_val.to_frame()
            X_tr_tmb_res, y_tr_res = oversample_to_balance(X_tr_tmb, y_train)
            
            for clf_key, clf_model in classifiers.items():
                clf = clone(clf_model)
                clf.fit(X_tr_tmb_res, y_tr_res)
                
                tr_prob = clf.predict_proba(X_tr_tmb_res)[:, 1]
                train_auc = roc_auc_score(y_tr_res, tr_prob)
                train_aucs[(clf_key, 'tmb_only')].append(train_auc)
                
                val_prob = clf.predict_proba(X_val_tmb)[:, 1]
                y_probas_val[(clf_key, 'tmb_only')][val_idx] = val_prob
                val_auc = roc_auc_score(y_val, val_prob)
                val_aucs[(clf_key, 'tmb_only')].append(val_auc)
                
            # Find genes mutated in >= 3 samples in the training split
            train_gene_counts = X_train.sum(axis=0)
            train_selected_genes = train_gene_counts[train_gene_counts >= 3].index.tolist()
            if len(train_selected_genes) == 0:
                train_selected_genes = train_gene_counts.sort_values(ascending=False).head(5).index.tolist()
                
            X_train_filtered = X_train[train_selected_genes]
            sorted_mutated_genes = train_gene_counts.sort_values(ascending=False)
            
            # Run Fisher test on training split
            train_fisher_pvals = {}
            for gene in train_selected_genes:
                mut_status = X_train_filtered[gene]
                mut_R = sum((mut_status == 1) & (y_train == 1))
                wt_R = sum((mut_status == 0) & (y_train == 1))
                mut_NR = sum((mut_status == 1) & (y_train == 0))
                wt_NR = sum((mut_status == 0) & (y_train == 0))
                _, p_val = fisher_exact([[mut_R, mut_NR], [wt_R, wt_NR]], alternative='two-sided')
                train_fisher_pvals[gene] = p_val
                
            sorted_fisher_genes = sorted(train_fisher_pvals.keys(), key=lambda g: train_fisher_pvals[g])
            
            # Define features mapping for each strategy in this fold
            strat_features = {
                'all': train_selected_genes,
                'mut_5': sorted_mutated_genes.head(5).index.tolist(),
                'mut_10': sorted_mutated_genes.head(10).index.tolist(),
                'mut_20': sorted_mutated_genes.head(20).index.tolist(),
                'mut_100': sorted_mutated_genes.head(100).index.tolist(),
                'mut_1000': sorted_mutated_genes.head(1000).index.tolist(),
                'fish_5': sorted_fisher_genes[:5],
                'fish_10': sorted_fisher_genes[:10],
                'fish_20': sorted_fisher_genes[:20],
                'fish_100': sorted_fisher_genes[:100],
                'fish_1000': sorted_fisher_genes[:1000]
            }
            
            # Accumulate Fisher p-values for this fold
            for strat in strategies:
                features = strat_features[strat]
                for gene in features:
                    pval = train_fisher_pvals.get(gene, 1.0)
                    cv_fisher_pvals[strat].setdefault(gene, []).append(pval)
                    
            # Train and predict for all strategies (with/without TMB)
            for strat in strategies:
                features = strat_features[strat]
                
                # A. WITHOUT TMB
                X_tr_no = X_train[features]
                X_val_no = X_val[features]
                X_tr_no_res, y_tr_res = oversample_to_balance(X_tr_no, y_train)
                
                for clf_key, clf_model in classifiers.items():
                    clf = clone(clf_model)
                    clf.fit(X_tr_no_res, y_tr_res)
                    
                    tr_prob = clf.predict_proba(X_tr_no_res)[:, 1]
                    train_auc = roc_auc_score(y_tr_res, tr_prob)
                    train_aucs[(clf_key, strat, False)].append(train_auc)
                    
                    val_prob = clf.predict_proba(X_val_no)[:, 1]
                    y_probas_val[(clf_key, strat, False)][val_idx] = val_prob
                    val_auc = roc_auc_score(y_val, val_prob)
                    val_aucs[(clf_key, strat, False)].append(val_auc)
                    
                    # Accumulate RF feature importances
                    if clf_key == 'rf':
                        importances_no = clf.feature_importances_
                        for gene, imp in zip(features, importances_no):
                            cv_importances[(strat, False)].setdefault(gene, []).append(imp)
                            
                # B. WITH TMB
                X_tr_yes = X_train[features].copy()
                X_tr_yes['TMB_Score'] = tmb_train
                X_val_yes = X_val[features].copy()
                X_val_yes['TMB_Score'] = tmb_val
                X_tr_yes_res, y_tr_res = oversample_to_balance(X_tr_yes, y_train)
                
                for clf_key, clf_model in classifiers.items():
                    clf = clone(clf_model)
                    clf.fit(X_tr_yes_res, y_tr_res)
                    
                    tr_prob = clf.predict_proba(X_tr_yes_res)[:, 1]
                    train_auc = roc_auc_score(y_tr_res, tr_prob)
                    train_aucs[(clf_key, strat, True)].append(train_auc)
                    
                    val_prob = clf.predict_proba(X_val_yes)[:, 1]
                    y_probas_val[(clf_key, strat, True)][val_idx] = val_prob
                    val_auc = roc_auc_score(y_val, val_prob)
                    val_aucs[(clf_key, strat, True)].append(val_auc)
                    
                    # Accumulate RF feature importances
                    if clf_key == 'rf':
                        importances_yes = clf.feature_importances_
                        features_yes = features + ['TMB_Score']
                        for feat, imp in zip(features_yes, importances_yes):
                            cv_importances[(strat, True)].setdefault(feat, []).append(imp)
                            
        # Calculate summary ROCs and mean AUCs
        rocs_data = {}
        performance_metrics_list = []
        
        train_means = {}
        val_means = {}
        
        x_labels = [
            'TMB Only',
            'All', 'All+TMB',
            'Mut 5', 'Mut 5+TMB',
            'Mut 10', 'Mut 10+TMB',
            'Mut 20', 'Mut 20+TMB',
            'Mut 100', 'Mut 100+TMB',
            'Mut 1000', 'Mut 1000+TMB',
            'Fish 5', 'Fish 5+TMB',
            'Fish 10', 'Fish 10+TMB',
            'Fish 20', 'Fish 20+TMB',
            'Fish 100', 'Fish 100+TMB',
            'Fish 1000', 'Fish 1000+TMB'
        ]
        
        for clf_key in classifiers:
            # TMB Baseline
            mean_tr_tmb = np.mean(train_aucs[(clf_key, 'tmb_only')])
            mean_val_tmb = np.mean(val_aucs[(clf_key, 'tmb_only')])
            train_means[(clf_key, 'TMB Only')] = mean_tr_tmb
            val_means[(clf_key, 'TMB Only')] = mean_val_tmb
            
            fpr_tmb, tpr_tmb, _ = roc_curve(y_encoded, y_probas_val[(clf_key, 'tmb_only')])
            rocs_data[(clf_key, 'tmb_only')] = (fpr_tmb, tpr_tmb, mean_val_tmb)
            
            perf_row = {
                'Cohort': name,
                'Classifier': clf_key,
                'AUC_TMB_only': mean_val_tmb
            }
            
            # Strategies
            for strat in strategies:
                for has_tmb in [False, True]:
                    suffix = "+TMB" if has_tmb else ""
                    lbl = f"{strat.capitalize().replace('_', ' ')}{suffix}"
                    
                    mean_tr = np.mean(train_aucs[(clf_key, strat, has_tmb)])
                    mean_val = np.mean(val_aucs[(clf_key, strat, has_tmb)])
                    
                    train_means[(clf_key, lbl)] = mean_tr
                    val_means[(clf_key, lbl)] = mean_val
                    
                    # Store OOF ROC curves
                    y_prob = y_probas_val[(clf_key, strat, has_tmb)]
                    fpr, tpr, _ = roc_curve(y_encoded, y_prob)
                    rocs_data[(clf_key, strat, has_tmb)] = (fpr, tpr, mean_val)
                    
                    col_name = f"AUC_{strat}_with_TMB" if has_tmb else f"AUC_{strat}"
                    perf_row[col_name] = mean_val
                    
            performance_metrics_list.append(perf_row)
            
        # Average the RF feature importances and Fisher p-values
        avg_importances_no = { strat: { gene: np.mean(imps) for gene, imps in cv_importances[(strat, False)].items() } for strat in strategies }
        avg_importances_yes = { strat: { feat: np.mean(imps) for feat, imps in cv_importances[(strat, True)].items() } for strat in strategies }
        avg_fisher_pvals = { strat: { gene: np.mean(pvals) for gene, pvals in cv_fisher_pvals[strat].items() } for strat in strategies }
        
    except Exception as e:
        print(f"  Warning: Cross-validation failed for {name}: {e}")
        import traceback
        traceback.print_exc()
        return None, None
        
    # Fit final global RF model on all data for feature importances (trained on all qualifying genes without TMB)
    clf_global = RandomForestClassifier(n_estimators=100, random_state=42)
    clf_global.fit(X_filtered, y_encoded)
    importances = clf_global.feature_importances_
    df_importances = pd.DataFrame({
        'Gene': selected_genes,
        'RF_Importance': importances
    })
    
    df_cohort = pd.merge(df_fisher, df_importances, on='Gene')
    df_cohort['Cohort'] = name
    cols = ['Cohort', 'Gene', 'Mutated_R', 'WT_R', 'Mutated_NR', 'WT_NR', 'Odds_Ratio', 'p_value', 'RF_Importance']
    df_cohort = df_cohort[cols]
    
    # 1. Plot comparison ROC grids for all 4 classifiers
    clf_titles = {
        'rf': 'Random Forest Classifier',
        'adaline': 'Adaline (Single-Layer NN)',
        'mlp4_relu': 'MLP-4 (ReLU Activation)',
        'mlp4_sigmoid': 'MLP-4 (Sigmoid Activation)'
    }
    
    for clf_key, clf_title in clf_titles.items():
        roc_grid_path = cohort_dir / f"{clf_key}_roc_curves_comparison.svg"
        plot_roc_grid(name, rocs_data, n_samples, r_count, nr_count, roc_grid_path, clf_title, clf_key)
        
    # 2. Plot Overfitting/Underfitting comparison
    overfit_plot_path = cohort_dir / "overfitting_underfitting_comparison.svg"
    print(f"  -> Plotting overfitting diagnostic curves to {overfit_plot_path.name}...")
    plot_overfitting_underfitting(name, train_means, val_means, overfit_plot_path)
    
    # 3. Plot individual RF feature importances and Fisher p-values per strategy
    labels_map = {
        'all': 'All Qualifying Genes',
        'mut_5': 'Top 5 Mutated',
        'mut_10': 'Top 10 Mutated',
        'mut_20': 'Top 20 Mutated',
        'mut_100': 'Top 100 Mutated',
        'mut_1000': 'Top 1000 Mutated',
        'fish_5': 'Top 5 Fisher',
        'fish_10': 'Top 10 Fisher',
        'fish_20': 'Top 20 Fisher',
        'fish_100': 'Top 100 Fisher',
        'fish_1000': 'Top 1000 Fisher'
    }
    
    for strat in strategies:
        imp_plot_path = imp_dir / f"{strat}.svg"
        plot_strategy_importances(name, strat, labels_map[strat], avg_importances_no[strat], avg_importances_yes[strat], imp_plot_path)
        
        fish_plot_path = fish_dir / f"{strat}.svg"
        plot_strategy_fisher_pvalues(name, strat, labels_map[strat], avg_fisher_pvals[strat], fish_plot_path)
        
    return df_cohort, performance_metrics_list


def main():
    lair_path = "/storage/halu/lair"
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    # Resolve directories
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output"
    resp_mut_output_dir = output_dir / "iAtlas-response-mutations"
    resp_mut_output_dir.mkdir(parents=True, exist_ok=True)
    
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
    
    all_results = []
    performances_list = []
    
    # 1. Run Individual Cohort Analyses
    print("\n--- Running Individual Cohort Analyses ---")
    for cohort in individual_cohorts:
        X, y, tmb = get_cohort_data(lair, [cohort])
        if X is not None and y is not None and tmb is not None:
            df_res, perf_metrics = analyze_cohort_mutations(X, y, tmb, cohort, resp_mut_output_dir)
            if df_res is not None:
                all_results.append(df_res)
                performances_list.extend(perf_metrics)
                
    # 2. Run Combined Cohort Analyses (Melanoma and Kidney)
    print("\n--- Running Combined Cohort Analyses ---")
    for comb_name, cohorts in cancer_type_combinations.items():
        X, y, tmb = get_cohort_data(lair, cohorts)
        if X is not None and y is not None and tmb is not None:
            df_res, perf_metrics = analyze_cohort_mutations(X, y, tmb, comb_name, resp_mut_output_dir)
            if df_res is not None:
                all_results.append(df_res)
                performances_list.extend(perf_metrics)
                
    # 3. Combine and Save Results
    if all_results:
        df_summary = pd.concat(all_results, axis=0).reset_index(drop=True)
        summary_path = resp_mut_output_dir / "predictor_mutations_summary.csv"
        df_summary.to_csv(summary_path, index=False)
        print(f"\nSaved combined prediction results summary to {summary_path}")
        
        # Save model performances comparison to CSV
        df_perf = pd.DataFrame(performances_list)
        perf_path = resp_mut_output_dir / "model_performances.csv"
        df_perf.to_csv(perf_path, index=False)
        print(f"Saved model performances comparison to {perf_path}")
        
        # Print a concise subset of performances for readability
        print("\nModel Performances Comparison (ROC AUC Subset):")
        subset_cols = ['Cohort', 'Classifier', 'AUC_TMB_only', 'AUC_all', 'AUC_all_with_TMB', 'AUC_mut_10', 'AUC_mut_10_with_TMB', 'AUC_fish_10', 'AUC_fish_10_with_TMB']
        print(df_perf[subset_cols].to_string(index=False))
        
        # Display top predictor gene for each analyzed dataset
        print("\nTop Predictor Gene (by Random Forest Importance) per Cohort:")
        for name in df_summary['Cohort'].unique():
            sub_df = df_summary[df_summary['Cohort'] == name]
            top_row = sub_df.sort_values(by='RF_Importance', ascending=False).iloc[0]
            print(f"  {name}: {top_row['Gene']} (Importance={top_row['RF_Importance']:.4f}, Fisher p-value={top_row['p_value']:.4f}, OR={top_row['Odds_Ratio']:.2f})")
    else:
        print("\nNo cohorts could be analyzed.")


if __name__ == "__main__":
    main()
