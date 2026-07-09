#!/usr/bin/env python3
"""
Grid search optimization for VAF and TMB thresholds in iAtlas datasets.
Binarizes cohorts into high/low TMB sets to maximize MCC, F1, or Accuracy.
Applies Log-Normal priors defined on the positive real axis (centered at VAF = 0.05 and TMB = 10.0).
Computes Fisher's Exact Test, Odds Ratios, and 95% Confidence Intervals.
Generates forest plots of Odds Ratios across cohorts and 4-panel dashboards.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from sklearn.metrics import matthews_corrcoef, f1_score, accuracy_score, roc_auc_score, precision_recall_curve, auc

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
        print(f"Error loading cohort {ds_name}: {e}")
        return None, None


def optimize_thresholds(df_mut, df_clinical, metric_name, lambda_val=0.20):
    """
    Grid search optimization for VAF and TMB thresholds with Log-Normal priors.
    Runs Fisher's exact test for binarized subgroups and calculates 95% Confidence Intervals.
    """
    clinical_samples = df_clinical.index
    y_true = (df_clinical['response'] == 'R').astype(int)
    
    # Grid definitions (strictly positive real axis)
    vaf_grid = np.linspace(0.01, 0.50, 50)
    tmb_grid = np.linspace(0.5, 50.0, 100)  # Avoid TMB = 0 to prevent log(0)
    
    best_score = -999.0
    best_v = 0.05
    best_t = 10.0
    best_raw_metric = 0.0
    
    # Pre-calculate TMB matrices for speed
    tmb_cache = {}
    for v in vaf_grid:
        tmb_matrix = calculate_maf_tmb(df_mut, capture_size_mb=30.0, vaf_threshold=v)
        tmb_matrix = tmb_matrix.set_index('Tumor_Sample_Barcode')
        tmb_scores = tmb_matrix['TMB_Score'].reindex(clinical_samples, fill_value=0.0)
        tmb_cache[v] = tmb_scores
        
    for v in vaf_grid:
        tmb_scores = tmb_cache[v]
        for t in tmb_grid:
            # Predict Responder (1) if TMB >= t
            y_pred = (tmb_scores >= t).astype(int)
            
            # Compute raw metric
            if metric_name == 'MCC':
                metric_val = matthews_corrcoef(y_true, y_pred)
            elif metric_name == 'F1':
                metric_val = f1_score(y_true, y_pred, zero_division=0)
            elif metric_name == 'Accuracy':
                metric_val = accuracy_score(y_true, y_pred)
            else:
                metric_val = accuracy_score(y_true, y_pred)
                
            # Log-Normal prior penalty (Gaussian in log-space)
            v_val = max(v, 1e-4)
            t_val = max(t, 1e-4)
            
            penalty = lambda_val * (
                ((np.log(v_val) - np.log(0.05)) / 0.5) ** 2 +
                ((np.log(t_val) - np.log(10.0)) / 0.5) ** 2
            )
            reg_score = metric_val - penalty
            
            if reg_score > best_score:
                best_score = reg_score
                best_raw_metric = metric_val
                best_v = v
                best_t = t
                
    # Calculate performance for the optimized threshold TMB
    opt_tmb = tmb_cache[best_v]
    if len(np.unique(y_true)) > 1:
        opt_roc = roc_auc_score(y_true, opt_tmb)
        precision, recall, _ = precision_recall_curve(y_true, opt_tmb)
        opt_pr = auc(recall, precision)
    else:
        opt_roc = np.nan
        opt_pr = np.nan
        
    # Calculate Fisher's exact test and Odds Ratio for binarized groups
    y_pred_opt = (opt_tmb >= best_t).astype(int)
    
    tp = sum((y_true == 1) & (y_pred_opt == 1))  # High TMB Responder
    fp = sum((y_true == 0) & (y_pred_opt == 1))  # High TMB Non-Responder
    fn = sum((y_true == 1) & (y_pred_opt == 0))  # Low TMB Responder
    tn = sum((y_true == 0) & (y_pred_opt == 0))  # Low TMB Non-Responder
    
    odds_ratio, fisher_pval = fisher_exact([[tp, fp], [fn, tn]])
    
    # Compute 95% Confidence Interval for Odds Ratio
    if odds_ratio > 0 and (tp > 0 or fp > 0) and (fn > 0 or tn > 0):
        # Wald confidence interval in log space (with 0.5 correction to prevent division by zero)
        se_log_or = np.sqrt(1.0/max(tp, 0.5) + 1.0/max(fp, 0.5) + 1.0/max(fn, 0.5) + 1.0/max(tn, 0.5))
        ci_lower = np.exp(np.log(odds_ratio) - 1.96 * se_log_or)
        ci_upper = np.exp(np.log(odds_ratio) + 1.96 * se_log_or)
    else:
        ci_lower = 0.0
        ci_upper = np.nan
        
    return {
        'Optimized_VAF': best_v,
        'Optimized_TMB': best_t,
        'Regularized_Score': best_score,
        'Raw_Metric': best_raw_metric,
        'ROC_AUC': opt_roc,
        'PR_AUC': opt_pr,
        'Fisher_p_value': fisher_pval,
        'Odds_Ratio': odds_ratio,
        'Odds_Ratio_Lower_CI': ci_lower,
        'Odds_Ratio_Upper_CI': ci_upper,
        'TMB_Scores': opt_tmb
    }


def plot_threshold_dashboard(df_mut, df_clinical, opt_res, cohort_name, metric_name, plot_path):
    """
    Generate the 4-panel threshold dashboard with Log-Normal priors and independent X-axes.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=300)
    
    # 1. Panel 1: Priors (VAF prior on bottom, TMB prior on top)
    ax_prior_v = axes[0, 0]
    
    # Plot VAF Log-Normal Prior (bottom X-axis)
    v_vals = np.linspace(0.005, 0.25, 200)
    v_dens = (1.0 / v_vals) * np.exp(-((np.log(v_vals) - np.log(0.05)) ** 2) / (2 * 0.5 ** 2))
    v_dens /= v_dens.max()
    line_v = ax_prior_v.plot(v_vals, v_dens, color='#E74C3C', linewidth=2.0, label='VAF Prior (Mode=0.05)')
    ax_prior_v.set_xlabel("VAF Threshold (Bottom Axis)", color='#E74C3C', fontsize=9.5)
    ax_prior_v.tick_params(axis='x', labelcolor='#E74C3C')
    ax_prior_v.set_ylabel("Normalized Prior Density", fontsize=9.5)
    ax_prior_v.set_xlim(0.0, 0.25)
    
    # Plot TMB Log-Normal Prior (top X-axis)
    ax_prior_t = ax_prior_v.twiny()
    t_vals = np.linspace(1.0, 50.0, 200)
    t_dens = (1.0 / t_vals) * np.exp(-((np.log(t_vals) - np.log(10.0)) ** 2) / (2 * 0.5 ** 2))
    t_dens /= t_dens.max()
    line_t = ax_prior_t.plot(t_vals, t_dens, color='#3498DB', linewidth=2.0, linestyle='--', label='TMB Prior (Mode=10.0)')
    ax_prior_t.set_xlabel("TMB Threshold (Top Axis)", color='#3498DB', fontsize=9.5)
    ax_prior_t.tick_params(axis='x', labelcolor='#3498DB')
    ax_prior_t.set_xlim(0.0, 50.0)
    
    # Combine legends
    lines = line_v + line_t
    labels = [l.get_label() for l in lines]
    ax_prior_v.legend(lines, labels, loc='upper right', fontsize=8.5)
    ax_prior_v.set_title("Log-Normal Priors on Positive Real Axis", fontsize=11, fontweight='bold', color='#2C3E50', pad=25)
    ax_prior_v.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
    
    # 2. Panel 2: VAF Distribution
    ax_vaf = axes[0, 1]
    df_mut_filtered = df_mut[df_mut['t_depth'] > 0].copy()
    df_mut_filtered['VAF'] = df_mut_filtered['t_alt_count'] / df_mut_filtered['t_depth']
    df_mut_filtered = df_mut_filtered[df_mut_filtered['Tumor_Sample_Barcode'].isin(df_clinical.index)]
    
    sns.histplot(df_mut_filtered['VAF'], bins=40, kde=True, ax=ax_vaf, color='#34495E', stat='density')
    ax_vaf.axvline(x=opt_res['Optimized_VAF'], color='#E74C3C', linestyle='-', linewidth=2.0, 
                   label=f"Optimized Threshold ({opt_res['Optimized_VAF']:.2f})")
    ax_vaf.set_xlim(0, 0.6)
    ax_vaf.set_xlabel("Somatic Mutation VAF", fontsize=9.5)
    ax_vaf.set_ylabel("Density", fontsize=9.5)
    ax_vaf.set_title("Somatic Mutation VAF Distribution", fontsize=11, fontweight='bold', color='#2C3E50')
    ax_vaf.legend(loc='upper right', fontsize=8.5)
    ax_vaf.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
    for spine in ['top', 'right']:
        ax_vaf.spines[spine].set_visible(False)
        
    # 3. Panel 3: TMB Distribution (Fixed Y-axis at [0, 40])
    ax_tmb = axes[1, 0]
    sns.histplot(opt_res['TMB_Scores'], bins=20, ax=ax_tmb, color='#9B59B6', stat='count')
    ax_tmb.axvline(x=opt_res['Optimized_TMB'], color='#2ECC71', linestyle='-', linewidth=2.0, 
                   label=f"Optimized Split ({opt_res['Optimized_TMB']:.1f})")
    ax_tmb.set_ylim(0, 40)  # FIXED Y-AXIS
    ax_tmb.set_xlabel("TMB Score (Mutations/Mb)", fontsize=9.5)
    ax_tmb.set_ylabel("Patient Count", fontsize=9.5)
    ax_tmb.set_title("TMB Distribution (Calculated at Optimized VAF)", fontsize=11, fontweight='bold', color='#2C3E50')
    ax_tmb.legend(loc='upper right', fontsize=8.5)
    ax_tmb.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
    for spine in ['top', 'right']:
        ax_tmb.spines[spine].set_visible(False)
        
    # 4. Panel 4: Performance Comparison Bar Plot
    ax_perf = axes[1, 1]
    perf_metrics = {
        f"Raw {metric_name}": opt_res['Raw_Metric'],
        "ROC AUC": opt_res['ROC_AUC'],
        "PR AUC": opt_res['PR_AUC']
    }
    
    colors = ['#E67E22', '#3498DB', '#2ECC71']
    bars = ax_perf.bar(perf_metrics.keys(), perf_metrics.values(), color=colors, edgecolor='black', width=0.5)
    ax_perf.set_ylim(0.0, 1.05)
    ax_perf.set_ylabel("Score / Performance", fontsize=9.5)
    ax_perf.set_title(f"Optimized TMB Split Performance ({metric_name})", fontsize=11, fontweight='bold', color='#2C3E50')
    
    # Add values on top of bars
    for bar in bars:
        height = bar.get_height()
        if not np.isnan(height):
            ax_perf.text(bar.get_x() + bar.get_width()/2.0, height + 0.02, f"{height:.3f}", 
                        ha='center', va='bottom', fontsize=9.5, fontweight='bold', color='#2C3E50')
                        
    ax_perf.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5, axis='y')
    for spine in ['top', 'right']:
        ax_perf.spines[spine].set_visible(False)
        
    plt.suptitle(f"Strongly Regularized VAF & TMB Log-Normal Optimization ({metric_name}): {cohort_name}\n"
                 f"Fisher p = {opt_res['Fisher_p_value']:.4g}, Odds Ratio = {opt_res['Odds_Ratio']:.3f} "
                 f"(95% CI: [{opt_res['Odds_Ratio_Lower_CI']:.2f}, {opt_res['Odds_Ratio_Upper_CI']:.2f}])", 
                 fontsize=12, fontweight='bold', color='#2C3E50', y=0.985)
                 
    plt.tight_layout()
    plt.savefig(plot_path, format='svg', bbox_inches='tight')
    plt.close(fig)


def plot_odds_ratios_forest(df_res, metric_name, plot_path):
    """
    Generate a Forest Plot of Odds Ratios across the cohorts for a given metric.
    """
    sub_df = df_res[df_res['Optimization_Metric'] == metric_name].copy()
    
    # Filter out entries where CI upper is NaN or infinite
    sub_df = sub_df.dropna(subset=['Odds_Ratio', 'Odds_Ratio_Lower_CI', 'Odds_Ratio_Upper_CI'])
    sub_df = sub_df[np.isfinite(sub_df['Odds_Ratio_Upper_CI'])]
    
    if len(sub_df) == 0:
        return
        
    plt.figure(figsize=(10, 6), dpi=300)
    
    # Reindex for clean y-axis order
    sub_df = sub_df.reset_index(drop=True)
    y_pos = np.arange(len(sub_df))
    
    # Determine significance colors: Green if significant positive (OR > 1 and p < 0.05), else Gray
    colors = []
    for _, row in sub_df.iterrows():
        if row['Fisher_p_value'] < 0.05 and row['Odds_Ratio'] > 1.0:
            colors.append('#2ECC71')  # Significant benefit
        else:
            colors.append('#7F8C8D')  # Not significant / negative
            
    # Plot error bars
    for idx, row in sub_df.iterrows():
        err_lower = row['Odds_Ratio'] - row['Odds_Ratio_Lower_CI']
        err_upper = row['Odds_Ratio_Upper_CI'] - row['Odds_Ratio']
        plt.errorbar(row['Odds_Ratio'], idx, xerr=[[err_lower], [err_upper]], 
                     fmt='o', color=colors[idx], ecolor=colors[idx], capsize=5, 
                     markersize=8, markeredgecolor='black', elinewidth=1.8)
                     
    # Add cohort labels on Y-axis
    cohort_labels = [f"{row['Cohort'].replace('-iAtlas', '')}\n(p={row['Fisher_p_value']:.3g})" for _, row in sub_df.iterrows()]
    plt.yticks(y_pos, cohort_labels, fontsize=9, fontweight='semibold', color='#34495E')
    
    # Logarithmic scale for Odds Ratio
    plt.xscale('log')
    plt.axvline(x=1.0, color='#E74C3C', linestyle='--', linewidth=1.5, label='No Association (OR=1)')
    
    plt.title(f"Clinical Association of Binarized TMB: Forest Plot of Odds Ratios ({metric_name})", 
              fontsize=12, fontweight='bold', color='#2C3E50', pad=15)
    plt.xlabel("Odds Ratio (Log Scale, 95% Confidence Interval)", fontsize=10)
    plt.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5, axis='x')
    plt.legend(loc='lower right', fontsize=8.5)
    
    # Annotate details
    for idx, row in sub_df.iterrows():
        plt.text(row['Odds_Ratio_Upper_CI'] * 1.05 if not pd.isna(row['Odds_Ratio_Upper_CI']) else row['Odds_Ratio'] * 1.5, 
                 idx, f"OR={row['Odds_Ratio']:.2f}\n[{row['Odds_Ratio_Lower_CI']:.2f}, {row['Odds_Ratio_Upper_CI']:.2f}]", 
                 va='center', ha='left', fontsize=8, color='#34495E', fontweight='bold')
                 
    sns.despine()
    plt.tight_layout()
    plt.savefig(plot_path, format='svg', bbox_inches='tight')
    plt.close()


def main():
    lair_path = "/storage/halu/lair"
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    # Resolve output directory
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output" / "iAtlas-vaf-tmb-optimization"
    output_dir.mkdir(parents=True, exist_ok=True)
    
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
    
    metrics = ['MCC', 'F1', 'Accuracy']
    lambda_val = 0.20  # Strong clinical prior
    
    optimization_results = []
    
    print("\n--- Running Strong Log-Normal VAF/TMB Threshold Optimizations (with Fisher's Test) ---")
    for cohort in individual_cohorts:
        df_mut, df_clin = get_cohort_mutations_and_clinical(lair, cohort)
        if df_mut is not None and df_clin is not None:
            for metric in metrics:
                print(f"  Optimizing {cohort} for {metric}...")
                opt_res = optimize_thresholds(df_mut, df_clin, metric, lambda_val)
                
                # Save plot
                plot_path = output_dir / f"{cohort}_opt_{metric}.svg"
                plot_threshold_dashboard(df_mut, df_clin, opt_res, cohort, metric, plot_path)
                
                optimization_results.append({
                    'Cohort': cohort,
                    'Optimization_Metric': metric,
                    'Optimized_VAF': opt_res['Optimized_VAF'],
                    'Optimized_TMB': opt_res['Optimized_TMB'],
                    'Regularized_Score': opt_res['Regularized_Score'],
                    'Raw_Metric_Score': opt_res['Raw_Metric'],
                    'Continuous_ROC_AUC': opt_res['ROC_AUC'],
                    'Continuous_PR_AUC': opt_res['PR_AUC'],
                    'Fisher_p_value': opt_res['Fisher_p_value'],
                    'Odds_Ratio': opt_res['Odds_Ratio'],
                    'Odds_Ratio_Lower_CI': opt_res['Odds_Ratio_Lower_CI'],
                    'Odds_Ratio_Upper_CI': opt_res['Odds_Ratio_Upper_CI']
                })
                
    if optimization_results:
        df_res = pd.DataFrame(optimization_results)
        csv_path = output_dir / "threshold_optimization_results.csv"
        df_res.to_csv(csv_path, index=False)
        print(f"\nSaved optimization statistics to {csv_path}")
        
        # Generate Forest Plots of Odds Ratios for each metric
        for metric in metrics:
            forest_plot_path = output_dir / f"odds_ratios_forest_plot_{metric}.svg"
            plot_odds_ratios_forest(df_res, metric, forest_plot_path)
            print(f"  Saved Forest Plot for {metric} to {forest_plot_path.name}")
            
        # Print summary of optimization results
        print("\nOptimized Threshold Summary (Strong Log-Normal Prior + Fisher's Stats):")
        for metric in metrics:
            print(f"\nMetric: {metric}")
            sub = df_res[df_res['Optimization_Metric'] == metric]
            for _, row in sub.iterrows():
                print(f"  {row['Cohort']}: VAF = {row['Optimized_VAF']:.2f}, TMB = {row['Optimized_TMB']:.1f} "
                      f"(Metric: {row['Raw_Metric_Score']:.3f}, ROC AUC: {row['Continuous_ROC_AUC']:.3f}, "
                      f"Fisher p = {row['Fisher_p_value']:.4g}, Odds Ratio = {row['Odds_Ratio']:.3f} "
                      f"CI: [{row['Odds_Ratio_Lower_CI']:.2f}, {row['Odds_Ratio_Upper_CI']:.2f}])")
    else:
        print("No cohorts could be optimized.")


if __name__ == "__main__":
    main()
