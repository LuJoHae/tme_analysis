#!/usr/bin/env python3
"""
Evaluate the statistical association of TMB and response as a function of the VAF threshold.
Calculates TMB for VAF thresholds from 0.01 to 0.50 in steps of 0.01.
Runs Mann-Whitney U tests, Cohen's d effect sizes, ROC AUCs, and PR AUCs for each threshold.
Saves results and generates cohort-specific 3-panel dashboards.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc

# Ensure local packages are in python path if running directly
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))

import datalair
import ici_datasets
from gene_utils import calculate_maf_tmb


def get_cohort_mutations_and_clinical(lair, ds_name):
    """
    Load mutations and response labels for a cohort.
    """
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    try:
        data_dir = ici_datasets.cbioportal_datasets.get_dataset_dir(lair, dataset_class, ds_name)
        df_clinical = ici_datasets.cbioportal_datasets.load_data_clinical(data_dir)
        df_clinical = df_clinical[df_clinical['response'].isin(['R', 'NR'])]
        
        mut_file = data_dir / "data_mutations.txt"
        if not mut_file.exists():
            return None, None
            
        df_mut = pd.read_csv(mut_file, sep="\t", low_memory=False)
        return df_mut, df_clinical
    except Exception as e:
        print(f"Error loading data for {ds_name}: {e}")
        return None, None


def run_vaf_tmb_statistics(df_mut, df_clinical, ds_name, vaf_thresholds):
    """
    Compute TMB at different VAF thresholds, run Mann-Whitney U, Cohen's d, ROC AUC, and PR AUC.
    """
    clinical_samples = df_clinical.index
    y_encoded = (df_clinical['response'] == 'R').astype(int)
    
    stats_list = []
    
    for thresh in vaf_thresholds:
        # Calculate TMB Score at this threshold
        tmb_matrix = calculate_maf_tmb(df_mut, capture_size_mb=30.0, vaf_threshold=thresh)
        tmb_matrix = tmb_matrix.set_index('Tumor_Sample_Barcode')
        
        # Align with clinical samples
        tmb_scores = tmb_matrix['TMB_Score'].reindex(clinical_samples, fill_value=0.0)
        
        # Split by response status
        tmb_r = tmb_scores[df_clinical['response'] == 'R']
        tmb_nr = tmb_scores[df_clinical['response'] == 'NR']
        
        # Mann-Whitney U test
        if len(tmb_r) > 1 and len(tmb_nr) > 1:
            stat, pval = mannwhitneyu(tmb_r, tmb_nr, alternative='two-sided')
        else:
            pval = np.nan
            
        # Cohen's d effect size
        n_r = len(tmb_r)
        n_nr = len(tmb_nr)
        if n_r > 1 and n_nr > 1:
            std_r = tmb_r.std()
            std_nr = tmb_nr.std()
            pooled_std = np.sqrt(((n_r - 1) * (std_r ** 2) + (n_nr - 1) * (std_nr ** 2)) / (n_r + n_nr - 2))
            if pooled_std > 0:
                cohens_d = (tmb_r.mean() - tmb_nr.mean()) / pooled_std
            else:
                cohens_d = 0.0
        else:
            cohens_d = np.nan
            
        # Single-feature ROC AUC and PR AUC
        if len(np.unique(y_encoded)) > 1:
            roc_auc = roc_auc_score(y_encoded, tmb_scores)
            precision, recall, _ = precision_recall_curve(y_encoded, tmb_scores)
            pr_auc = auc(recall, precision)
        else:
            roc_auc = np.nan
            pr_auc = np.nan
            
        stats_list.append({
            'Cohort': ds_name,
            'VAF_Threshold': thresh,
            'R_Mean_TMB': tmb_r.mean(),
            'R_Median_TMB': tmb_r.median(),
            'NR_Mean_TMB': tmb_nr.mean(),
            'NR_Median_TMB': tmb_nr.median(),
            'Mann_Whitney_U_p_value': pval,
            'Neg_Log_P': -np.log10(pval) if not pd.isna(pval) and pval > 0 else np.nan,
            'Cohens_d': cohens_d,
            'TMB_ROC_AUC': roc_auc,
            'TMB_PR_AUC': pr_auc
        })
        
    return pd.DataFrame(stats_list)


def plot_vaf_tmb_association(df_stats, cohort_name, plot_path, r_count, total_samples):
    """
    Plot a 3-panel line chart showing:
    1. -log10(p-value) from Mann-Whitney U test.
    2. Cohen's d effect size.
    3. ROC AUC & PR AUC with their respective baseline references.
    """
    fig, axes = plt.subplots(1, 3, figsize=(21, 5.5), dpi=300)
    
    vaf = df_stats['VAF_Threshold']
    neg_log_p = df_stats['Neg_Log_P']
    cohens_d = df_stats['Cohens_d']
    roc_auc = df_stats['TMB_ROC_AUC']
    pr_auc = df_stats['TMB_PR_AUC']
    
    # 1. Panel 1: Statistical Significance
    ax_p = axes[0]
    ax_p.plot(vaf, neg_log_p, color='#E74C3C', linewidth=2.0, marker='o', markersize=3, label='-log10(p-val)')
    ax_p.axhline(y=-np.log10(0.05), color='#7F8C8D', linestyle='--', linewidth=1.2, label='p = 0.05')
    
    # Highlight max significance
    max_idx_p = neg_log_p.idxmax()
    if not pd.isna(max_idx_p):
        best_vaf_p = vaf.iloc[max_idx_p]
        best_p_val = df_stats['Mann_Whitney_U_p_value'].iloc[max_idx_p]
        ax_p.axvline(x=best_vaf_p, color='#2C3E50', linestyle=':', alpha=0.7)
        ax_p.plot(best_vaf_p, neg_log_p.iloc[max_idx_p], 'o', color='#2C3E50', markersize=6)
        ax_p.text(best_vaf_p + 0.01, neg_log_p.iloc[max_idx_p], f"Best VAF={best_vaf_p:.2f}\np={best_p_val:.3g}",
                  va='center', ha='left', fontsize=8.5, fontweight='bold', color='#2C3E50')
                  
    ax_p.set_title("Statistical Significance vs. VAF", fontsize=11, fontweight='bold', color='#2C3E50')
    ax_p.set_xlabel("VAF Threshold for TMB calculation", fontsize=9.5)
    ax_p.set_ylabel("-log10(Mann-Whitney U p-value)", fontsize=9.5)
    ax_p.set_xlim(0.0, 0.51)
    ax_p.legend(loc='upper right', fontsize=8.5)
    ax_p.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
    for spine in ['top', 'right']:
        ax_p.spines[spine].set_visible(False)
        
    # 2. Panel 2: Cohen's d Effect Size
    ax_d = axes[1]
    ax_d.plot(vaf, cohens_d, color='#F39C12', linewidth=2.0, marker='o', markersize=3, label="Cohen's d")
    ax_d.axhline(y=0.0, color='#7F8C8D', linestyle='--', linewidth=1.2)
    
    # Highlight max magnitude of Cohen's d (can be positive or negative)
    abs_d = cohens_d.abs()
    max_idx_d = abs_d.idxmax()
    if not pd.isna(max_idx_d):
        best_vaf_d = vaf.iloc[max_idx_d]
        best_d = cohens_d.iloc[max_idx_d]
        ax_d.axvline(x=best_vaf_d, color='#2C3E50', linestyle=':', alpha=0.7)
        ax_d.plot(best_vaf_d, best_d, 'o', color='#2C3E50', markersize=6)
        ax_d.text(best_vaf_d + 0.01, best_d, f"Best VAF={best_vaf_d:.2f}\nd={best_d:.2f}",
                  va='center', ha='left', fontsize=8.5, fontweight='bold', color='#2C3E50')
                  
    ax_d.set_title("Cohen's d Effect Size vs. VAF", fontsize=11, fontweight='bold', color='#2C3E50')
    ax_d.set_xlabel("VAF Threshold for TMB calculation", fontsize=9.5)
    ax_d.set_ylabel("Cohen's d", fontsize=9.5)
    ax_d.set_xlim(0.0, 0.51)
    ax_d.legend(loc='upper right', fontsize=8.5)
    ax_d.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
    for spine in ['top', 'right']:
        ax_d.spines[spine].set_visible(False)
        
    # 3. Panel 3: Performance (ROC AUC & PR AUC)
    ax_a = axes[2]
    ax_a.plot(vaf, roc_auc, color='#3498DB', linewidth=2.0, marker='o', markersize=3, label='TMB ROC AUC')
    ax_a.plot(vaf, pr_auc, color='#2ECC71', linewidth=2.0, marker='s', markersize=3, label='TMB PR AUC')
    
    # Reference baselines
    ax_a.axhline(y=0.5, color='#3498DB', linestyle='--', linewidth=1.0, alpha=0.7, label='ROC Baseline (0.5)')
    baseline_pr = r_count / total_samples
    ax_a.axhline(y=baseline_pr, color='#2ECC71', linestyle='--', linewidth=1.0, alpha=0.7, label=f'PR Baseline ({baseline_pr:.2f})')
    
    # Highlight max ROC AUC
    max_idx_a = roc_auc.idxmax()
    if not pd.isna(max_idx_a):
        best_vaf_a = vaf.iloc[max_idx_a]
        best_auc = roc_auc.iloc[max_idx_a]
        ax_a.axvline(x=best_vaf_a, color='#2C3E50', linestyle=':', alpha=0.7)
        ax_a.plot(best_vaf_a, best_auc, 'o', color='#2C3E50', markersize=6)
        ax_a.text(best_vaf_a + 0.01, best_auc, f"Best VAF={best_vaf_a:.2f}\nROC={best_auc:.3f}",
                  va='center', ha='left', fontsize=8.5, fontweight='bold', color='#2C3E50')
                  
    ax_a.set_title("ROC AUC & PR AUC vs. VAF", fontsize=11, fontweight='bold', color='#2C3E50')
    ax_a.set_xlabel("VAF Threshold for TMB calculation", fontsize=9.5)
    ax_a.set_ylabel("AUC Score", fontsize=9.5)
    ax_a.set_xlim(0.0, 0.51)
    ax_a.set_ylim(0.1, 1.05)
    ax_a.legend(loc='lower right', fontsize=8)
    ax_a.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
    for spine in ['top', 'right']:
        ax_a.spines[spine].set_visible(False)
        
    plt.suptitle(f"TMB Association & Effect Size by VAF Threshold: {cohort_name}", 
                 fontsize=13, fontweight='bold', color='#2C3E50', y=0.985)
                 
    plt.tight_layout()
    plt.savefig(plot_path, format='svg', bbox_inches='tight')
    plt.close(fig)


def main():
    lair_path = "/storage/halu/lair"
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    # Resolve output directory
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output" / "iAtlas-vaf-tmb-stats"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Grid of VAF thresholds
    vaf_thresholds = np.arange(0.01, 0.51, 0.01)
    
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
    
    all_stats_dfs = []
    
    # 1. Run Individual Cohorts
    print("\n--- Running Individual Cohort Statistical Analyses ---")
    for cohort in individual_cohorts:
        df_mut, df_clin = get_cohort_mutations_and_clinical(lair, cohort)
        if df_mut is not None and df_clin is not None:
            df_c_stats = run_vaf_tmb_statistics(df_mut, df_clin, cohort, vaf_thresholds)
            all_stats_dfs.append(df_c_stats)
            
            # Plot
            plot_path = output_dir / f"{cohort}_vaf_tmb_association.svg"
            r_count = sum(df_clin['response'] == 'R')
            total = len(df_clin)
            plot_vaf_tmb_association(df_c_stats, cohort, plot_path, r_count, total)
            print(f"  Saved plot for {cohort} to {plot_path.name}")
            
    # 2. Run Combined Cohorts
    print("\n--- Running Combined Cohort Statistical Analyses ---")
    for comb_name, cohorts in cancer_type_combinations.items():
        # Merge mutations and clinical datasets
        merged_mut_list = []
        merged_clin_list = []
        for c in cohorts:
            df_mut, df_clin = get_cohort_mutations_and_clinical(lair, c)
            if df_mut is not None and df_clin is not None:
                df_clin_c = df_clin.copy()
                df_clin_c.index = f"{c}_" + df_clin_c.index.astype(str)
                df_mut_c = df_mut.copy()
                df_mut_c['Tumor_Sample_Barcode'] = f"{c}_" + df_mut_c['Tumor_Sample_Barcode'].astype(str)
                
                merged_mut_list.append(df_mut_c)
                merged_clin_list.append(df_clin_c)
                
        if merged_mut_list:
            df_mut_comb = pd.concat(merged_mut_list, axis=0, ignore_index=True)
            df_clin_comb = pd.concat(merged_clin_list, axis=0)
            
            df_comb_stats = run_vaf_tmb_statistics(df_mut_comb, df_clin_comb, comb_name, vaf_thresholds)
            all_stats_dfs.append(df_comb_stats)
            
            # Plot
            plot_path = output_dir / f"{comb_name}_vaf_tmb_association.svg"
            r_count = sum(df_clin_comb['response'] == 'R')
            total = len(df_clin_comb)
            plot_vaf_tmb_association(df_comb_stats, comb_name, plot_path, r_count, total)
            print(f"  Saved plot for {comb_name} to {plot_path.name}")
            
    if all_stats_dfs:
        df_all_stats = pd.concat(all_stats_dfs, axis=0).reset_index(drop=True)
        csv_path = output_dir / "vaf_tmb_association_stats.csv"
        df_all_stats.to_csv(csv_path, index=False)
        print(f"\nSaved statistical summaries to {csv_path}")
        
        # Print top threshold per cohort
        print("\nOptimal TMB VAF Threshold per Cohort (by p-value):")
        for cohort in df_all_stats['Cohort'].unique():
            sub = df_all_stats[df_all_stats['Cohort'] == cohort]
            best_row = sub.sort_values(by='Mann_Whitney_U_p_value').iloc[0]
            print(f"  {cohort}: Best VAF = {best_row['VAF_Threshold']:.2f} (p = {best_row['Mann_Whitney_U_p_value']:.4g}, Cohen's d = {best_row['Cohens_d']:.3f}, ROC AUC = {best_row['TMB_ROC_AUC']:.3f}, PR AUC = {best_row['TMB_PR_AUC']:.3f})")
    else:
        print("No cohorts could be analyzed.")


if __name__ == "__main__":
    main()
