#!/usr/bin/env python3
"""
Evaluate the impact of TGM6 somatic mutations in the combined melanoma cohort
on response to immunotherapy.
Calculates Fisher's Exact Test, Mann-Whitney U test, and runs stratified CV
comparing predictive power of TMB vs. TGM6 mutation vs. TMB + TGM6.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact, mannwhitneyu
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, roc_auc_score

# Ensure local packages are in python path if running directly
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))

import datalair
import ici_datasets


def get_melanoma_combined_data(lair):
    """
    Load and combine clinical response, mutations, and TMB for melanoma cohorts.
    """
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    mel_cohorts = ['Riaz-iAtlas', 'Liu-iAtlas', 'Hugo-iAtlas']
    
    all_X = []
    all_y = []
    all_tmb = []
    
    for ds_name in mel_cohorts:
        try:
            data_dir = ici_datasets.cbioportal_datasets.get_dataset_dir(lair, dataset_class, ds_name)
            df_clinical = ici_datasets.cbioportal_datasets.load_data_clinical(data_dir)
            df_clinical = df_clinical[df_clinical['response'].isin(['R', 'NR'])]
            
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
            
            # Filter TGM6 mutations
            tgm6_muts = qualifying_subset[qualifying_subset['Hugo_Symbol'] == 'TGM6']
            mutated_samples = tgm6_muts['Tumor_Sample_Barcode'].unique()
            
            # Create a Series of TGM6 mutation status for all clinical samples
            tgm6_status = pd.Series(0, index=df_clinical.index)
            tgm6_status.loc[tgm6_status.index.isin(mutated_samples)] = 1
            
            # Align indices
            tgm6_status.index = f"{ds_name}_" + tgm6_status.index.astype(str)
            y = df_clinical['response'].copy()
            y.index = f"{ds_name}_" + y.index.astype(str)
            tmb_score = df_clinical['TMB_Score'].copy()
            tmb_score.index = f"{ds_name}_" + tmb_score.index.astype(str)
            
            all_X.append(tgm6_status)
            all_y.append(y)
            all_tmb.append(tmb_score)
            
        except Exception as e:
            print(f"Error loading {ds_name}: {e}")
            
    if not all_X:
        return None
        
    combined_df = pd.DataFrame({
        'TGM6_Mutation': pd.concat(all_X, axis=0),
        'Response': pd.concat(all_y, axis=0),
        'TMB_Score': pd.concat(all_tmb, axis=0)
    })
    
    return combined_df


def main():
    lair_path = "/storage/halu/lair"
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    # Resolve directories
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output" / "TGM6-analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("\nLoading and combining melanoma cohort data...")
    df = get_melanoma_combined_data(lair)
    if df is None:
        print("Error: Could not retrieve combined melanoma cohort data.")
        return
        
    n_samples = len(df)
    print(f"Combined Melanoma Cohort Size: n = {n_samples}")
    
    # 1. Fisher's Exact Test for TGM6 Mutation vs. Response
    # Contingency Table:
    #                 Responder  Non-Responder
    # TGM6 Mutated       a            b
    # TGM6 Wildtype      c            d
    
    mut_r = sum((df['TGM6_Mutation'] == 1) & (df['Response'] == 'R'))
    mut_nr = sum((df['TGM6_Mutation'] == 1) & (df['Response'] == 'NR'))
    wt_r = sum((df['TGM6_Mutation'] == 0) & (df['Response'] == 'R'))
    wt_nr = sum((df['TGM6_Mutation'] == 0) & (df['Response'] == 'NR'))
    
    odds_ratio, fisher_p = fisher_exact([[mut_r, mut_nr], [wt_r, wt_nr]], alternative='two-sided')
    
    print("\n--- Somatic Mutation Response Association (Fisher's Exact Test) ---")
    print(f"Contingency Table:")
    print(f"                  Responder (R)    Non-Responder (NR)")
    print(f"  TGM6 Mutated:        {mut_r}                 {mut_nr}")
    print(f"  TGM6 Wildtype:       {wt_r}                {wt_nr}")
    print(f"Odds Ratio: {odds_ratio:.4f}")
    print(f"Fisher's p-value: {fisher_p:.4g}")
    
    # 2. TMB Comparison: TGM6 Mutated vs. Wild-type
    tmb_mut = df[df['TGM6_Mutation'] == 1]['TMB_Score']
    tmb_wt = df[df['TGM6_Mutation'] == 0]['TMB_Score']
    
    mwu_stat, mwu_p = mannwhitneyu(tmb_mut, tmb_wt, alternative='two-sided')
    print("\n--- Tumor Mutational Burden Comparison (Mann-Whitney U Test) ---")
    print(f"TGM6 Mutated TMB: Mean = {tmb_mut.mean():.2f}, Median = {tmb_mut.median():.2f}")
    print(f"TGM6 Wildtype TMB: Mean = {tmb_wt.mean():.2f}, Median = {tmb_wt.median():.2f}")
    print(f"Mann-Whitney U p-value: {mwu_p:.4g}")
    
    # 3. Multivariate Predictive Analysis using Stratified 5-Fold Cross-Validation
    y_encoded = (df['Response'] == 'R').astype(int)
    
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    # Pre-allocate for predictions
    y_proba_tmb = np.zeros(n_samples)
    y_proba_tgm6 = np.zeros(n_samples)
    y_proba_combined = np.zeros(n_samples)
    
    for train_idx, val_idx in cv.split(df, y_encoded):
        # Splits
        train_df, val_df = df.iloc[train_idx], df.iloc[val_idx]
        y_train, y_val = y_encoded.iloc[train_idx], y_encoded.iloc[val_idx]
        
        # A. TMB Only
        clf_tmb = LogisticRegression(class_weight='balanced', random_state=42)
        clf_tmb.fit(train_df[['TMB_Score']], y_train)
        y_proba_tmb[val_idx] = clf_tmb.predict_proba(val_df[['TMB_Score']])[:, 1]
        
        # B. TGM6 Mutation Only
        clf_tgm6 = LogisticRegression(class_weight='balanced', random_state=42)
        clf_tgm6.fit(train_df[['TGM6_Mutation']], y_train)
        y_proba_tgm6[val_idx] = clf_tgm6.predict_proba(val_df[['TGM6_Mutation']])[:, 1]
        
        # C. Combined: TMB + TGM6
        clf_comb = LogisticRegression(class_weight='balanced', random_state=42)
        clf_comb.fit(train_df[['TMB_Score', 'TGM6_Mutation']], y_train)
        y_proba_combined[val_idx] = clf_comb.predict_proba(val_df[['TMB_Score', 'TGM6_Mutation']])[:, 1]
        
    auc_tmb = roc_auc_score(y_encoded, y_proba_tmb)
    auc_tgm6 = roc_auc_score(y_encoded, y_proba_tgm6)
    auc_combined = roc_auc_score(y_encoded, y_proba_combined)
    
    print("\n--- Predictive Performances (ROC AUC from 5-Fold Stratified CV) ---")
    print(f"  TMB-only Model ROC AUC:       {auc_tmb:.4f}")
    print(f"  TGM6-only Model ROC AUC:      {auc_tgm6:.4f}")
    print(f"  Combined (TMB+TGM6) ROC AUC:  {auc_combined:.4f}")
    
    # 4. Generate 3-Panel Visualizations
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), dpi=300)
    
    # --- Subplot A: Response Rates Stacked Bar ---
    ax_a = axes[0]
    # Build df for plotting proportions
    resp_rates = pd.DataFrame([
        {'TGM6_Status': 'Wildtype\n(n=205)', 'Response': 'R', 'Count': wt_r, 'Percentage': wt_r / (wt_r + wt_nr) * 100},
        {'TGM6_Status': 'Wildtype\n(n=205)', 'Response': 'NR', 'Count': wt_nr, 'Percentage': wt_nr / (wt_r + wt_nr) * 100},
        {'TGM6_Status': 'Mutated\n(n=14)', 'Response': 'R', 'Count': mut_r, 'Percentage': mut_r / (mut_r + mut_nr) * 100},
        {'TGM6_Status': 'Mutated\n(n=14)', 'Response': 'NR', 'Count': mut_nr, 'Percentage': mut_nr / (mut_r + mut_nr) * 100}
    ])
    
    # Pivot to plot stacked bars
    pivot_df = resp_rates.pivot(index='TGM6_Status', columns='Response', values='Percentage')
    pivot_df = pivot_df[['R', 'NR']] # order
    
    pivot_df.plot(kind='bar', stacked=True, color=['#2ECC71', '#E74C3C'], ax=ax_a, width=0.5, edgecolor='none')
    ax_a.set_title("Response Rate Comparison", fontsize=11, fontweight='bold', color='#2C3E50', pad=10)
    ax_a.set_ylabel("Percentage (%)", fontsize=9, fontweight='semibold')
    ax_a.set_xlabel("TGM6 Mutation Status", fontsize=9, fontweight='semibold')
    ax_a.set_xticklabels(ax_a.get_xticklabels(), rotation=0)
    ax_a.legend(['Responder (R)', 'Non-Responder (NR)'], loc='lower left', fontsize=8.5)
    
    # Annotate details on bars
    # Wildtype
    # ax_a.text(0, wt_r/2, f"{wt_r} ({(wt_r/(wt_r+wt_nr)*100):.1f}%)", ha='center', va='center', color='white', fontweight='bold', fontsize=9)
    # ax_a.text(0, wt_r + wt_nr/2, f"{wt_nr} ({(wt_nr/(wt_r+wt_nr)*100):.1f}%)", ha='center', va='center', color='white', fontweight='bold', fontsize=9)
    # Mutated
    # ax_a.text(1, mut_r/2, f"{mut_r} ({(mut_r/(mut_r+mut_nr)*100):.1f}%)", ha='center', va='center', color='white', fontweight='bold', fontsize=9)
    # ax_a.text(1, mut_r + mut_nr/2, f"{mut_nr} ({(mut_nr/(mut_r+mut_nr)*100):.1f}%)", ha='center', va='center', color='white', fontweight='bold', fontsize=9)
    
    # Annotate stats
    ax_a.text(0.5, 105, f"Fisher's p = {fisher_p:.3g}\nOdds Ratio = {odds_ratio:.2f}", ha='center', va='bottom', fontsize=9, fontweight='bold', color='#2980B9')
    ax_a.set_ylim(0, 120)
    for spine in ['top', 'right']:
        ax_a.spines[spine].set_visible(False)
        
    # --- Subplot B: TMB Comparison Boxplot ---
    ax_b = axes[1]
    df_plot_tmb = df.copy()
    df_plot_tmb['TGM6_Status'] = np.where(df_plot_tmb['TGM6_Mutation'] == 1, 'Mutated\n(n=14)', 'Wildtype\n(n=205)')
    
    sns.boxplot(x='TGM6_Status', y='TMB_Score', data=df_plot_tmb, palette=['#BDC3C7', '#F1C40F'], ax=ax_b, width=0.4, showfliers=False)
    sns.stripplot(x='TGM6_Status', y='TMB_Score', data=df_plot_tmb, color='#2C3E50', alpha=0.4, jitter=0.15, ax=ax_b, size=4)
    
    ax_b.set_title("TMB Distribution by TGM6 Status", fontsize=11, fontweight='bold', color='#2C3E50', pad=10)
    ax_b.set_ylabel("Tumor Mutational Burden (TMB Score)", fontsize=9, fontweight='semibold')
    ax_b.set_xlabel("TGM6 Mutation Status", fontsize=9, fontweight='semibold')
    
    # Annotate stats
    y_max = df['TMB_Score'].max()
    ax_b.text(0.5, y_max * 0.95, f"Mann-Whitney U p = {mwu_p:.3g}", ha='center', va='bottom', fontsize=9, fontweight='bold', color='#D35400')
    for spine in ['top', 'right']:
        ax_b.spines[spine].set_visible(False)
        
    # --- Subplot C: ROC Curves ---
    ax_c = axes[2]
    
    fpr_tmb, tpr_tmb, _ = roc_curve(y_encoded, y_proba_tmb)
    fpr_tgm6, tpr_tgm6, _ = roc_curve(y_encoded, y_proba_tgm6)
    fpr_comb, tpr_comb, _ = roc_curve(y_encoded, y_proba_combined)
    
    ax_c.plot(fpr_tmb, tpr_tmb, label=f"TMB Only (AUC={auc_tmb:.3f})", color='#E67E22', linewidth=2)
    ax_c.plot(fpr_tgm6, tpr_tgm6, label=f"TGM6 Only (AUC={auc_tgm6:.3f})", color='#8E44AD', linewidth=2)
    ax_c.plot(fpr_comb, tpr_comb, label=f"TMB + TGM6 (AUC={auc_combined:.3f})", color='#2980B9', linewidth=2)
    ax_c.plot([0, 1], [0, 1], 'k--', alpha=0.5)
    
    ax_c.set_title("ROC Curves: Immunotherapy Response", fontsize=11, fontweight='bold', color='#2C3E50', pad=10)
    ax_c.set_xlabel("False Positive Rate", fontsize=9, fontweight='semibold')
    ax_c.set_ylabel("True Positive Rate", fontsize=9, fontweight='semibold')
    ax_c.legend(loc="lower right", fontsize=8.5)
    ax_c.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
    for spine in ['top', 'right']:
        ax_c.spines[spine].set_visible(False)
        
    plt.suptitle("Impact of TGM6 Mutation in Combined Melanoma Cohort (Riaz + Liu + Hugo)", 
                 fontsize=14, fontweight='bold', color='#2C3E50', y=0.985)
                 
    plt.tight_layout()
    plot_path = output_dir / "tgm6_melanoma_impact.svg"
    plt.savefig(plot_path, format='svg', bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved premium 3-panel figure to {plot_path}")
    
    # 5. Save stats to CSV
    stats_df = pd.DataFrame([{
        'TGM6_Mut_R': mut_r,
        'TGM6_Mut_NR': mut_nr,
        'TGM6_WT_R': wt_r,
        'TGM6_WT_NR': wt_nr,
        'Fisher_Odds_Ratio': odds_ratio,
        'Fisher_p_value': fisher_p,
        'Mann_Whitney_U_p_value': mwu_p,
        'AUC_TMB_only': auc_tmb,
        'AUC_TGM6_only': auc_tgm6,
        'AUC_TMB_TGM6_Combined': auc_combined
    }])
    
    stats_path = output_dir / "tgm6_melanoma_stats.csv"
    stats_df.to_csv(stats_path, index=False)
    print(f"Saved stats summary to {stats_path}")


if __name__ == "__main__":
    main()
