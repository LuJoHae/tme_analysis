#!/usr/bin/env python3
"""
Bagaev-Constrained Deconvolution and Comparative Response Prediction.
Compares predicting immunotherapy response from Bagaev gene expression directly
vs. inferred cell fractions from Bagaev-constrained deconvolution.
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
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score, matthews_corrcoef
from sklearn.manifold import TSNE
import umap

# Ensure local packages are in path
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))
import datalair
from single_cell_datasets import SingleCellBagaevDeconvolution, SingleCellDeconvolution
from ici_datasets.bagaev_datasets import Signature as BagaevSignature
from gene_utils import read_gene_sets

# Setup output dir
output_dir = Path("output/bagaev-constrained-deconv")
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


def run_cross_val_rf(X, y, seed):
    """Run Stratified 5-Fold CV Random Forest and return OOF predictions."""
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    oof_preds = np.zeros(len(y))
    
    # Handle tiny classes
    n_classes = len(np.unique(y))
    if n_classes < 2:
        return oof_preds
    
    min_class_size = min(np.bincount(y))
    if min_class_size < 5:
        # Fallback to standard KFold if stratified is impossible
        from sklearn.model_selection import KFold
        skf = KFold(n_splits=5, shuffle=True, random_state=seed)
        
    for train_idx, val_idx in skf.split(X, y):
        X_train, X_val = X[train_idx], X[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]
        
        # Standardize features
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_val = scaler.transform(X_val)
        
        model = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=seed)
        model.fit(X_train, y_train)
        
        # Return probability of positive class
        oof_preds[val_idx] = model.predict_proba(X_val)[:, 1]
        
    return oof_preds


def main():
    print("Initializing Lair...")
    lair = datalair.Lair("/storage/halu/lair")
    
    # 1. safe_derive SingleCellBagaevDeconvolution
    print("\nDeriving Bagaev-Constrained Deconvolution Dataset...")
    ds_bagaev = SingleCellBagaevDeconvolution()
    lair.safe_derive(ds_bagaev)
    bagaev_deconv_paths = lair.get_dataset_filepaths(ds_bagaev)
    
    # Load also standard deconvolution to compare
    ds_standard = SingleCellDeconvolution()
    standard_deconv_paths = lair.get_dataset_filepaths(ds_standard)
    
    # Load Bagaev signature genes once
    ds_sig = BagaevSignature()
    lair.safe_derive(ds_sig)
    filepaths_sig = lair.get_dataset_filepaths(ds_sig)
    bagaev_signature = read_gene_sets(filepaths_sig["gene_signatures.gmt"])
    bagaev_genes = sorted(list(set.union(*[x.genes for x in bagaev_signature.values()])))
    bagaev_genes = [g.upper() for g in bagaev_genes]
    
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
    
    # We will use the kmeans_subcluster_res_0.5 reference
    ref_name = "kmeans_subcluster_res_0.5"
    
    # We will store results here
    results_rows = []
    
    # Collect expression and deconv data for UMAPs
    expression_matrices = {}
    deconv_matrices = {}
    clinical_labels = {}
    
    # Loop over cohorts
    for cohort in cohorts:
        print(f"\nProcessing cohort: {cohort}...")
        df_clin = load_clinical(cohort, lair)
        if df_clin is None:
            print(f"  No clinical data found for {cohort}. Skipping.")
            continue
            
        clinical_labels[cohort] = df_clin
        
        # A. Load aligned deconv cell fractions
        bagaev_deconv_key = f"deconv_bagaev_{cohort}_{ref_name}.csv"
        deconv_path = bagaev_deconv_paths.get(bagaev_deconv_key)
        if not deconv_path:
            print(f"  Bagaev deconv file not found for {cohort}. Skipping.")
            continue
            
        df_fracs = pd.read_csv(deconv_path, index_col=0)
        df_fracs.index = df_fracs.index.astype(str)
        
        # B. Load and reconstruct bulk expression for Bagaev genes
        import ici_datasets
        dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
        ds_cb = dataset_class(name=cohort)
        filepaths_cb = lair.get_dataset_filepaths(ds_cb)
        unpacked_key = next(f for f in filepaths_cb.keys() if not f.endswith('.tar.gz'))
        p_dir = filepaths_cb[unpacked_key] / filepaths_cb[unpacked_key].name
        
        # Determine target file
        target_file = None
        for f in ["data_mrna_seq_expression.txt", "data_mrna_seq_tpm.txt", "data_mrna_seq_rpkm.txt"]:
            if (p_dir / f).is_file():
                target_file = p_dir / f
                break
        if not target_file:
            print(f"  Expression file not found for {cohort}. Skipping.")
            continue
            
        df_expr = pd.read_csv(target_file, sep='\t', index_col=0)
        import numpy as np
        if df_expr.max().max() < 50:
            print(f"    Detected log-transformed expression for {cohort} (max={df_expr.max().max():.2f}). Linearizing...")
            df_expr = np.power(2, df_expr) - 1
        df_expr.index = df_expr.index.str.upper()
        df_expr = df_expr.groupby(level=0).mean()
        
        # Align with common samples
        common_samples = df_fracs.index.intersection(df_clin.index).intersection(df_expr.columns)
        if len(common_samples) < 10:
            print(f"  Too few common samples ({len(common_samples)}) for {cohort}. Skipping.")
            continue
            
        print(f"  Found {len(common_samples)} common samples.")
        
        # Subset and align
        y_series = df_clin.loc[common_samples, 'response'].map({'R': 1, 'NR': 0})
        df_fracs_aligned = df_fracs.loc[common_samples]
        
        common_bagaev = [g for g in bagaev_genes if g in df_expr.index]
        df_expr_aligned = df_expr.loc[common_bagaev, common_samples].T
        
        expression_matrices[cohort] = df_expr_aligned
        deconv_matrices[cohort] = df_fracs_aligned
        
        # Run Symmetric predictions across 5 seeds
        seeds = [42, 101, 2023, 7, 888]
        
        # Prepare arrays
        X_expr = df_expr_aligned.values
        X_fracs = df_fracs_aligned.values
        y = y_series.values
        
        for seed in seeds:
            # Mode A: Gene Expression Directly
            oof_expr = run_cross_val_rf(X_expr, y, seed)
            auc_expr = roc_auc_score(y, oof_expr)
            pr_expr = average_precision_score(y, oof_expr)
            pred_class_expr = (oof_expr >= 0.5).astype(int)
            acc_expr = accuracy_score(y, pred_class_expr)
            mcc_expr = matthews_corrcoef(y, pred_class_expr)
            
            results_rows.append({
                'Cohort': cohort,
                'Mode': 'Direct Expression',
                'Seed': seed,
                'ROC_AUC': auc_expr,
                'PR_AUC': pr_expr,
                'Accuracy': acc_expr,
                'MCC': mcc_expr
            })
            
            # Mode B: Deconvoluted Fractions
            oof_fracs = run_cross_val_rf(X_fracs, y, seed)
            auc_fracs = roc_auc_score(y, oof_fracs)
            pr_fracs = average_precision_score(y, oof_fracs)
            pred_class_fracs = (oof_fracs >= 0.5).astype(int)
            acc_fracs = accuracy_score(y, pred_class_fracs)
            mcc_fracs = matthews_corrcoef(y, pred_class_fracs)
            
            results_rows.append({
                'Cohort': cohort,
                'Mode': 'Deconvoluted Fractions',
                'Seed': seed,
                'ROC_AUC': auc_fracs,
                'PR_AUC': pr_fracs,
                'Accuracy': acc_fracs,
                'MCC': mcc_fracs
            })
    # Run Combined Cohorts
    combined_specs = [
        ('Combined-Melanoma', ['Hugo-iAtlas', 'Riaz-iAtlas', 'Liu-iAtlas', 'Gide-iAtlas']),
        ('Combined-RCC', ['McDermott-iAtlas', 'Choueiri-iAtlas']),
        ('Combined-All', ['Hugo-iAtlas', 'Riaz-iAtlas', 'Liu-iAtlas', 'Gide-iAtlas', 'Rosenberg-iAtlas', 'Padron-iAtlas', 'Anders-iAtlas', 'McDermott-iAtlas', 'Choueiri-iAtlas'])
    ]
    
    for comb_name, sub_cohorts in combined_specs:
        print(f"\nProcessing combined cohort: {comb_name}...")
        expr_list = []
        fracs_list = []
        y_list = []
        
        for c in sub_cohorts:
            if c in expression_matrices:
                expr_list.append(expression_matrices[c])
                fracs_list.append(deconv_matrices[c])
                df_clin = clinical_labels[c]
                common = expression_matrices[c].index
                y_list.append(df_clin.loc[common, 'response'].map({'R': 1, 'NR': 0}))
                
        if not expr_list:
            continue
            
        df_expr_comb = pd.concat(expr_list, axis=0)
        df_fracs_comb = pd.concat(fracs_list, axis=0)
        y_series_comb = pd.concat(y_list, axis=0)
        
        X_expr = df_expr_comb.values
        X_fracs = df_fracs_comb.values
        y = y_series_comb.values
        
        seeds = [42, 101, 2023, 7, 888]
        for seed in seeds:
            # Mode A: Gene Expression Directly
            oof_expr = run_cross_val_rf(X_expr, y, seed)
            auc_expr = roc_auc_score(y, oof_expr)
            pr_expr = average_precision_score(y, oof_expr)
            pred_class_expr = (oof_expr >= 0.5).astype(int)
            acc_expr = accuracy_score(y, pred_class_expr)
            mcc_expr = matthews_corrcoef(y, pred_class_expr)
            
            results_rows.append({
                'Cohort': comb_name,
                'Mode': 'Direct Expression',
                'Seed': seed,
                'ROC_AUC': auc_expr,
                'PR_AUC': pr_expr,
                'Accuracy': acc_expr,
                'MCC': mcc_expr
            })
            
            # Mode B: Deconvoluted Fractions
            oof_fracs = run_cross_val_rf(X_fracs, y, seed)
            auc_fracs = roc_auc_score(y, oof_fracs)
            pr_fracs = average_precision_score(y, oof_fracs)
            pred_class_fracs = (oof_fracs >= 0.5).astype(int)
            acc_fracs = accuracy_score(y, pred_class_fracs)
            mcc_fracs = matthews_corrcoef(y, pred_class_fracs)
            
            results_rows.append({
                'Cohort': comb_name,
                'Mode': 'Deconvoluted Fractions',
                'Seed': seed,
                'ROC_AUC': auc_fracs,
                'PR_AUC': pr_fracs,
                'Accuracy': acc_fracs,
                'MCC': mcc_fracs
            })

    # Save comparative metrics to CSV
    df_results = pd.DataFrame(results_rows)
    df_results.to_csv(output_dir / "bagaev_prediction_results.csv", index=False)
    print(f"\nSaved prediction results to {output_dir / 'bagaev_prediction_results.csv'}")
    
    # 2. Plot performance comparison with error bars
    df_summary = df_results.groupby(['Cohort', 'Mode']).agg({
        'ROC_AUC': ['mean', 'std'],
        'PR_AUC': ['mean', 'std'],
        'Accuracy': ['mean', 'std'],
        'MCC': ['mean', 'std']
    }).reset_index()
    df_summary.columns = ['Cohort', 'Mode', 'ROC_AUC_mean', 'ROC_AUC_std', 'PR_AUC_mean', 'PR_AUC_std', 'Accuracy_mean', 'Accuracy_std', 'MCC_mean', 'MCC_std']
    
    combined_names = ['Combined-Melanoma', 'Combined-RCC', 'Combined-All']
    cohort_order = [c for c in cohorts if c not in combined_names] + combined_names

    fig, axes = plt.subplots(2, 2, figsize=(20, 12))
    metrics = ['ROC_AUC', 'PR_AUC', 'Accuracy', 'MCC']
    
    for idx, metric in enumerate(metrics):
        ax = axes[idx // 2, idx % 2]
        sns.barplot(
            data=df_summary,
            x='Cohort',
            y=f'{metric}_mean',
            hue='Mode',
            ax=ax,
            palette='Set2',
            order=cohort_order
        )
        
        # Add error bars manually
        # Get x coordinates of bars
        hue_labels = [t.get_text() for t in ax.get_legend().get_texts()]
        n_cohorts = len(ax.get_xticks())
        for patch_idx, bar in enumerate(ax.patches):
            x_pos = bar.get_x() + bar.get_width() / 2.0
            
            # Find closest x-tick to identify cohort
            tick_positions = np.array(ax.get_xticks())
            closest_tick_idx = np.argmin(np.abs(tick_positions - x_pos))
            cohort_name = ax.get_xticklabels()[closest_tick_idx].get_text()
            
            # Identify mode based on patch index
            mode_idx = patch_idx // n_cohorts
            if mode_idx < len(hue_labels):
                mode_name = hue_labels[mode_idx]
                row = df_summary[(df_summary['Cohort'] == cohort_name) & (df_summary['Mode'] == mode_name)]
                if not row.empty:
                    mean_val = row[f'{metric}_mean'].values[0]
                    std_val = row[f'{metric}_std'].values[0]
                    ax.errorbar(x_pos, mean_val, yerr=std_val, fmt='none', c='black', capsize=4)
                    
        ax.set_title(f"Comparative {metric.replace('_', ' ')}")
        ax.set_ylabel(metric.replace('_', ' '))
        ax.tick_params(axis='x', rotation=30)
        ax.set_ylim(bottom=0.0 if metric != 'MCC' else -0.5, top=1.0)
        
    plt.tight_layout()
    plt.savefig(output_dir / "bagaev_prediction_comparison_grid.svg", bbox_inches='tight')
    plt.close()
    print("Saved comparative metrics grid.")

    # 3. Create UMAP of expressions
    print("Generating UMAPs...")
    all_expr_list = []
    all_fracs_list = []
    expr_cohort_labels = []
    frac_cohort_labels = []
    
    for cohort in expression_matrices:
        all_expr_list.append(expression_matrices[cohort])
        expr_cohort_labels.extend([cohort] * len(expression_matrices[cohort]))
        
    for cohort in deconv_matrices:
        all_fracs_list.append(deconv_matrices[cohort])
        frac_cohort_labels.extend([cohort] * len(deconv_matrices[cohort]))
        
    if all_expr_list:
        X_expr_all = pd.concat(all_expr_list, axis=0)
        # Run UMAP on log-transformed expressions
        X_expr_scaled = StandardScaler().fit_transform(np.log1p(X_expr_all.values))
        reducer_expr = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
        emb_expr = reducer_expr.fit_transform(X_expr_scaled)
        
        plt.figure(figsize=(10, 8))
        sns.scatterplot(x=emb_expr[:, 0], y=emb_expr[:, 1], hue=expr_cohort_labels, palette='tab10', alpha=0.8)
        plt.title("UMAP of Bagaev-Constrained Gene Expressions (Bulk)")
        plt.savefig(output_dir / "umap_bagaev_gene_expressions.svg", bbox_inches='tight')
        plt.close()
        
    if all_fracs_list:
        X_fracs_all = pd.concat(all_fracs_list, axis=0)
        X_fracs_scaled = StandardScaler().fit_transform(X_fracs_all.values)
        reducer_fracs = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
        emb_fracs = reducer_fracs.fit_transform(X_fracs_scaled)
        
        plt.figure(figsize=(10, 8))
        sns.scatterplot(x=emb_fracs[:, 0], y=emb_fracs[:, 1], hue=frac_cohort_labels, palette='tab10', alpha=0.8)
        plt.title("UMAP of Bagaev-Deconvoluted Cell Fractions")
        plt.savefig(output_dir / "umap_bagaev_cell_fractions.svg", bbox_inches='tight')
        plt.close()

    # 4. Compare standard deconvolution with Bagaev deconvolution fractions
    print("Comparing standard vs. Bagaev cell fractions...")
    correlations_list = []
    
    for cohort in deconv_matrices:
        bagaev_df = deconv_matrices[cohort]
        
        # Load standard deconv fractions
        std_deconv_key = f"deconv_{cohort}_{ref_name}.csv"
        std_path = standard_deconv_paths.get(std_deconv_key)
        if std_path:
            std_df = pd.read_csv(std_path, index_col=0)
            std_df.index = std_df.index.astype(str)
            common_s = std_df.index.intersection(bagaev_df.index)
            if len(common_s) > 5:
                # Compute Pearson correlation for each cell type
                std_aligned = std_df.loc[common_s]
                bagaev_aligned = bagaev_df.loc[common_s]
                
                for col in std_aligned.columns:
                    if col in bagaev_aligned.columns:
                        r = std_aligned[col].corr(bagaev_aligned[col])
                        correlations_list.append({
                            'Cohort': cohort,
                            'Cell Type': col,
                            'Correlation': r
                        })
                        
    if correlations_list:
        df_corr = pd.DataFrame(correlations_list)
        plt.figure(figsize=(14, 8))
        # Plot correlation distribution by cell type
        sns.boxplot(data=df_corr, x='Cell Type', y='Correlation', palette='vlag')
        plt.xticks(rotation=90)
        plt.title("Pearson Correlation of Cell Type Fractions: Standard vs. Bagaev-Constrained Deconvolution")
        plt.ylabel("Pearson Correlation")
        plt.ylim(-1.0, 1.0)
        plt.savefig(output_dir / "deconv_correlation_standard_vs_bagaev.svg", bbox_inches='tight')
        plt.close()
        
    # Plot detailed cell type fractions and major lineages
    plot_bagaev_deconvolution_distributions(bagaev_deconv_paths, ref_name, output_dir)
    plot_bagaev_major_lineages(bagaev_deconv_paths, ref_name, output_dir)
        
    print("Done!")


def plot_grouped_boxes(ax, df_merged, cell_types, group_col, group_labels, group_colors):
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


def plot_bagaev_deconvolution_distributions(bagaev_deconv_paths, ref_name, output_dir):
    print("Generating Bagaev cell fraction distributions...")
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
    
    first_key = f"deconv_bagaev_Hugo-iAtlas_{ref_name}.csv"
    first_path = bagaev_deconv_paths.get(first_key)
    if not first_path:
        return
    df_temp = pd.read_csv(first_path, index_col=0)
    cell_types = [c for c in df_temp.columns if c != 'cohort_group']
    
    fig, axes = plt.subplots(rows, cols, figsize=(24, 12 * rows), dpi=150)
    axes_flat = axes.flatten() if n_cohorts > 1 else [axes]
    
    for idx, cohort_name in enumerate(cohorts):
        tcga_project = COHORT_TO_TCGA[cohort_name]
        ax = axes_flat[idx]
        
        iatlas_key = f"deconv_bagaev_{cohort_name}_{ref_name}.csv"
        tcga_key = f"deconv_bagaev_TCGA-{tcga_project}_{ref_name}.csv"
        
        iatlas_path = bagaev_deconv_paths.get(iatlas_key)
        tcga_path = bagaev_deconv_paths.get(tcga_key)
        
        if not iatlas_path:
            ax.text(0.5, 0.5, f"Missing data for {cohort_name}", ha='center', va='center')
            continue
            
        df_iatlas = pd.read_csv(iatlas_path, index_col=0)
        df_iatlas['cohort_group'] = cohort_name
        
        if tcga_path:
            df_tcga = pd.read_csv(tcga_path, index_col=0)
            df_tcga['cohort_group'] = f"TCGA-{tcga_project}"
            df_merged = pd.concat([df_iatlas, df_tcga], axis=0)
            group_labels = [cohort_name, f"TCGA-{tcga_project}"]
        else:
            df_merged = df_iatlas
            group_labels = [cohort_name]
            
        plot_grouped_boxes(ax, df_merged, cell_types, 'cohort_group', group_labels, ['#2ECC71', '#3498DB'])
        ax.set_title(f"{cohort_name} vs TCGA-{tcga_project} (Bagaev-Deconv)")
        ax.set_xlabel("Fraction")
        
    plt.tight_layout()
    plt.savefig(output_dir / f"bagaev_deconv_distributions_{ref_name}.svg", bbox_inches='tight')
    plt.close()


def plot_bagaev_major_lineages(bagaev_deconv_paths, ref_name, output_dir):
    print("Generating Bagaev major lineage distributions...")
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
    
    lineages = ['T_NK', 'B_Plasma', 'Myeloid', 'Endothelial', 'Fibroblasts', 'Tumor_Epithelial']
    
    fig, axes = plt.subplots(rows, cols, figsize=(24, 6 * rows), dpi=150)
    axes_flat = axes.flatten() if n_cohorts > 1 else [axes]
    
    for idx, cohort_name in enumerate(cohorts):
        tcga_project = COHORT_TO_TCGA[cohort_name]
        ax = axes_flat[idx]
        
        iatlas_key = f"deconv_bagaev_{cohort_name}_{ref_name}.csv"
        tcga_key = f"deconv_bagaev_TCGA-{tcga_project}_{ref_name}.csv"
        
        iatlas_path = bagaev_deconv_paths.get(iatlas_key)
        tcga_path = bagaev_deconv_paths.get(tcga_key)
        
        if not iatlas_path:
            ax.text(0.5, 0.5, f"Missing data for {cohort_name}", ha='center', va='center')
            continue
            
        df_iatlas = pd.read_csv(iatlas_path, index_col=0)
        df_iatlas_lin = pd.DataFrame(index=df_iatlas.index)
        for lin in lineages:
            matching_cols = [c for c in df_iatlas.columns if c.startswith(lin)]
            df_iatlas_lin[lin] = df_iatlas[matching_cols].sum(axis=1)
        df_iatlas_lin['cohort_group'] = cohort_name
        
        if tcga_path:
            df_tcga = pd.read_csv(tcga_path, index_col=0)
            df_tcga_lin = pd.DataFrame(index=df_tcga.index)
            for lin in lineages:
                matching_cols = [c for c in df_tcga.columns if c.startswith(lin)]
                df_tcga_lin[lin] = df_tcga[matching_cols].sum(axis=1)
            df_tcga_lin['cohort_group'] = f"TCGA-{tcga_project}"
            df_merged = pd.concat([df_iatlas_lin, df_tcga_lin], axis=0)
            group_labels = [cohort_name, f"TCGA-{tcga_project}"]
        else:
            df_merged = df_iatlas_lin
            group_labels = [cohort_name]
            
        plot_grouped_boxes(ax, df_merged, lineages, 'cohort_group', group_labels, ['#2ECC71', '#3498DB'])
        ax.set_title(f"{cohort_name} vs TCGA-{tcga_project} (Bagaev-Deconv Major Lineages)")
        ax.set_xlabel("Fraction")
        
    plt.tight_layout()
    plt.savefig(output_dir / "bagaev_deconv_major_lineages_distributions.svg", bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    main()
