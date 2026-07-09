#!/usr/bin/env python3
"""
Evaluate the impact of TGM6 expression on response to immunotherapy
across Riaz-iAtlas, Liu-iAtlas, Hugo-iAtlas, and the Combined Melanoma cohort.
Splits the analysis for samples with TGM6 wildtype and TGM6 mutated.
Generates 2x3 diagnostic plots and saves a compiled stats summary CSV.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact, mannwhitneyu

# Ensure local packages are in python path if running directly
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))

import datalair
import ici_datasets


def load_cohort_data(lair, ds_name):
    """
    Load clinical response, mutations, and mRNA expression data for a single dataset.
    """
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    try:
        data_dir = ici_datasets.cbioportal_datasets.get_dataset_dir(lair, dataset_class, ds_name)
        df_clin, df_mrna = ici_datasets.cbioportal_datasets.load_and_process_data(data_dir)
        df_clin = df_clin[df_clin['response'].isin(['R', 'NR'])]
        
        common_samples = df_clin.index.intersection(df_mrna.columns)
        if len(common_samples) == 0:
            return None
            
        expr = df_mrna.loc['TGM6', common_samples]
        resp = df_clin.loc[common_samples, 'response']
        
        # Load mutations to find mutated samples
        mutated_samples = []
        mut_file = data_dir / "data_mutations.txt"
        if mut_file.exists():
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
            
            tgm6_muts = qualifying_subset[qualifying_subset['Hugo_Symbol'] == 'TGM6']
            mutated_samples = tgm6_muts['Tumor_Sample_Barcode'].unique()
            
        df_res = pd.DataFrame({
            'TGM6_Expression': expr,
            'Response': resp,
            'TGM6_Mutation': np.where(expr.index.isin(mutated_samples), 1, 0)
        })
        df_res['Cohort'] = ds_name
        return df_res
    except Exception as e:
        print(f"Error loading {ds_name}: {e}")
        return None


def plot_subplot_boxplot(ax, df_sub, title):
    """
    Plot continuous expression boxplot on a specific subplot.
    """
    if len(df_sub) == 0:
        ax.text(0.5, 0.5, "No samples", ha='center', va='center', fontsize=10, color='#7F8C8D')
        ax.set_title(title, fontsize=10, fontweight='bold', color='#2C3E50')
        ax.axis('off')
        return
        
    r_sub = df_sub[df_sub['Response'] == 'R']['TGM6_Expression']
    nr_sub = df_sub[df_sub['Response'] == 'NR']['TGM6_Expression']
    
    # Handle single class plotting by adjusting palette dynamically
    colors = []
    if len(r_sub) > 0: colors.append('#2ECC71')
    if len(nr_sub) > 0: colors.append('#E74C3C')
    
    sns.boxplot(x='Response', y='TGM6_Expression', data=df_sub, palette=colors, ax=ax, width=0.4, showfliers=False)
    sns.stripplot(x='Response', y='TGM6_Expression', data=df_sub, color='#2C3E50', alpha=0.5, jitter=0.15, ax=ax, size=4)
    
    ax.set_title(title, fontsize=10, fontweight='bold', color='#2C3E50')
    ax.set_ylabel("TGM6 Expression (TPM)", fontsize=8, fontweight='semibold')
    ax.set_xlabel("Immunotherapy Response", fontsize=8, fontweight='semibold')
    
    # Custom x tick labels with counts
    ax.set_xticklabels([f"R (n={len(r_sub)})", f"NR (n={len(nr_sub)})"])
    
    # Calculate MWU test
    if len(r_sub) > 1 and len(nr_sub) > 1:
        _, pval = mannwhitneyu(r_sub, nr_sub, alternative='two-sided')
        y_max = max(df_sub['TGM6_Expression'].max(), 0.1)
        ax.text(0.5, y_max * 0.92, f"MWU p={pval:.3g}", ha='center', va='bottom', fontsize=8.5, fontweight='bold', color='#D35400')
        ax.set_ylim(-0.05, y_max * 1.1)
    else:
        y_max = max(df_sub['TGM6_Expression'].max(), 0.1)
        ax.text(0.5, y_max * 0.5, "MWU: N/A\n(Too few classes/samples)", ha='center', va='center', fontsize=8.5, fontweight='bold', color='#7F8C8D')
        
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)


def plot_subplot_stackedbar(ax, df_sub, title):
    """
    Plot stacked bar chart of response rate on a specific subplot.
    """
    if len(df_sub) == 0:
        ax.text(0.5, 0.5, "No samples", ha='center', va='center', fontsize=10, color='#7F8C8D')
        ax.set_title(title, fontsize=10, fontweight='bold', color='#2C3E50')
        ax.axis('off')
        return
        
    df_binary = df_sub.copy()
    df_binary['Expressed'] = (df_binary['TGM6_Expression'] > 0).astype(int)
    
    exp_r = sum((df_binary['Expressed'] == 1) & (df_binary['Response'] == 'R'))
    exp_nr = sum((df_binary['Expressed'] == 1) & (df_binary['Response'] == 'NR'))
    unexp_r = sum((df_binary['Expressed'] == 0) & (df_binary['Response'] == 'R'))
    unexp_nr = sum((df_binary['Expressed'] == 0) & (df_binary['Response'] == 'NR'))
    
    n_exp = exp_r + exp_nr
    n_unexp = unexp_r + unexp_nr
    
    pct_exp_r = exp_r / n_exp * 100 if n_exp > 0 else 0
    pct_exp_nr = exp_nr / n_exp * 100 if n_exp > 0 else 0
    pct_unexp_r = unexp_r / n_unexp * 100 if n_unexp > 0 else 0
    pct_unexp_nr = unexp_nr / n_unexp * 100 if n_unexp > 0 else 0
    
    resp_rates = pd.DataFrame([
        {'TGM6_Status': f'Unexpressed\n(n={n_unexp})', 'Response': 'R', 'Percentage': pct_unexp_r},
        {'TGM6_Status': f'Unexpressed\n(n={n_unexp})', 'Response': 'NR', 'Percentage': pct_unexp_nr},
        {'TGM6_Status': f'Expressed\n(n={n_exp})', 'Response': 'R', 'Percentage': pct_exp_r},
        {'TGM6_Status': f'Expressed\n(n={n_exp})', 'Response': 'NR', 'Percentage': pct_exp_nr}
    ])
    
    pivot_df = resp_rates.pivot(index='TGM6_Status', columns='Response', values='Percentage')
    for c in ['R', 'NR']:
        if c not in pivot_df.columns:
            pivot_df[c] = 0.0
    pivot_df = pivot_df[['R', 'NR']]
    
    pivot_df.plot(kind='bar', stacked=True, color=['#2ECC71', '#E74C3C'], ax=ax, width=0.5, edgecolor='none')
    ax.set_title(title, fontsize=10, fontweight='bold', color='#2C3E50')
    ax.set_ylabel("Percentage (%)", fontsize=8, fontweight='semibold')
    ax.set_xlabel("TGM6 Expression (TPM > 0)", fontsize=8, fontweight='semibold')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
    
    # Fisher test
    if n_exp > 0 and n_unexp > 0:
        try:
            odds_ratio, fisher_p = fisher_exact([[exp_r, exp_nr], [unexp_r, unexp_nr]], alternative='two-sided')
            ax.text(0.5, 105, f"Fisher p={fisher_p:.3g}\nOR={odds_ratio:.2f}", ha='center', va='bottom', fontsize=8.5, fontweight='bold', color='#2980B9')
        except Exception:
            ax.text(0.5, 105, "Fisher: N/A", ha='center', va='bottom', fontsize=8.5, fontweight='bold', color='#7F8C8D')
    else:
        ax.text(0.5, 105, "Fisher: N/A\n(No variation)", ha='center', va='bottom', fontsize=8.5, fontweight='bold', color='#7F8C8D')
        
    ax.set_ylim(0, 120)
    ax.legend(['R', 'NR'], loc='lower left', fontsize=7.5)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)


def plot_expression_impact_split(df, name, output_path):
    """
    Generate a 2x3 figure layout showing:
    - Row 1: Boxplots of TGM6 expression vs response (All, Wildtype, Mutated)
    - Row 2: Stacked bars of response rates (All, Wildtype, Mutated)
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10), dpi=300)
    
    df_wt = df[df['TGM6_Mutation'] == 0]
    df_mut = df[df['TGM6_Mutation'] == 1]
    
    # Row 1: Boxplots
    plot_subplot_boxplot(axes[0, 0], df, "All Samples (Continuous)")
    plot_subplot_boxplot(axes[0, 1], df_wt, "TGM6 Wildtype (Continuous)")
    plot_subplot_boxplot(axes[0, 2], df_mut, "TGM6 Mutated (Continuous)")
    
    # Row 2: Stacked Bars
    plot_subplot_stackedbar(axes[1, 0], df, "All Samples (Binarized)")
    plot_subplot_stackedbar(axes[1, 1], df_wt, "TGM6 Wildtype (Binarized)")
    plot_subplot_stackedbar(axes[1, 2], df_mut, "TGM6 Mutated (Binarized)")
    
    plt.suptitle(f"Impact of TGM6 mRNA Expression Split by TGM6 Mutation: {name}", 
                 fontsize=14, fontweight='bold', color='#2C3E50', y=0.98)
                 
    plt.tight_layout()
    plt.savefig(output_path, format='svg', bbox_inches='tight')
    plt.close(fig)


def calculate_metrics_for_split(df_sub, cohort_name, split_name):
    """
    Perform statistical tests and return dictionary of results for a specific split.
    """
    n_samples = len(df_sub)
    if n_samples == 0:
        return {
            'Cohort': cohort_name,
            'Split': split_name,
            'Total_Samples': 0,
            'Responders': 0,
            'Non_Responders': 0,
            'TGM6_Expressed_R': 0,
            'TGM6_Expressed_NR': 0,
            'TGM6_Unexpressed_R': 0,
            'TGM6_Unexpressed_NR': 0,
            'Fisher_Odds_Ratio': np.nan,
            'Fisher_p_value': np.nan,
            'Mann_Whitney_U_p_value': np.nan
        }
        
    r_expr = df_sub[df_sub['Response'] == 'R']['TGM6_Expression']
    nr_expr = df_sub[df_sub['Response'] == 'NR']['TGM6_Expression']
    
    if len(r_expr) > 1 and len(nr_expr) > 1:
        _, mwu_p = mannwhitneyu(r_expr, nr_expr, alternative='two-sided')
    else:
        mwu_p = np.nan
        
    df_binary = df_sub.copy()
    df_binary['Expressed'] = (df_binary['TGM6_Expression'] > 0).astype(int)
    
    exp_r = sum((df_binary['Expressed'] == 1) & (df_binary['Response'] == 'R'))
    exp_nr = sum((df_binary['Expressed'] == 1) & (df_binary['Response'] == 'NR'))
    unexp_r = sum((df_binary['Expressed'] == 0) & (df_binary['Response'] == 'R'))
    unexp_nr = sum((df_binary['Expressed'] == 0) & (df_binary['Response'] == 'NR'))
    
    n_exp = exp_r + exp_nr
    n_unexp = unexp_r + unexp_nr
    
    if n_exp > 0 and n_unexp > 0:
        try:
            odds_ratio, fisher_p = fisher_exact([[exp_r, exp_nr], [unexp_r, unexp_nr]], alternative='two-sided')
        except Exception:
            odds_ratio, fisher_p = np.nan, np.nan
    else:
        odds_ratio, fisher_p = np.nan, np.nan
        
    return {
        'Cohort': cohort_name,
        'Split': split_name,
        'Total_Samples': n_samples,
        'Responders': len(r_expr),
        'Non_Responders': len(nr_expr),
        'TGM6_Expressed_R': exp_r,
        'TGM6_Expressed_NR': exp_nr,
        'TGM6_Unexpressed_R': unexp_r,
        'TGM6_Unexpressed_NR': unexp_nr,
        'Fisher_Odds_Ratio': odds_ratio,
        'Fisher_p_value': fisher_p,
        'Mann_Whitney_U_p_value': mwu_p
    }


def analyze_cohort_split_expression(df, name, output_dir):
    """
    Run statistics and save split plots for a cohort or combined cohort.
    """
    df_wt = df[df['TGM6_Mutation'] == 0]
    df_mut = df[df['TGM6_Mutation'] == 1]
    
    print(f"\n--- Running Split Analysis for {name} ---")
    print(f"  All:      n = {len(df)}")
    print(f"  Wildtype: n = {len(df_wt)}")
    print(f"  Mutated:  n = {len(df_mut)}")
    
    m_all = calculate_metrics_for_split(df, name, 'All')
    m_wt = calculate_metrics_for_split(df_wt, name, 'Wildtype')
    m_mut = calculate_metrics_for_split(df_mut, name, 'Mutated')
    
    plot_path = output_dir / f"{name}_tgm6_expression_impact.svg"
    plot_expression_impact_split(df, name, plot_path)
    print(f"  Saved split plot to {plot_path.name}")
    
    return [m_all, m_wt, m_mut]


def main():
    lair_path = "/storage/halu/lair"
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    # Resolve directories
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output" / "TGM6-analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    mel_cohorts = ['Riaz-iAtlas', 'Liu-iAtlas', 'Hugo-iAtlas']
    
    all_dfs = []
    stats_list = []
    
    # 1. Analyze each dataset individually
    for ds_name in mel_cohorts:
        df_cohort = load_cohort_data(lair, ds_name)
        if df_cohort is not None:
            all_dfs.append(df_cohort)
            metrics_splits = analyze_cohort_split_expression(df_cohort, ds_name, output_dir)
            stats_list.extend(metrics_splits)
            
    # 2. Analyze combined melanoma cohort
    if all_dfs:
        combined_df = pd.concat(all_dfs, axis=0)
        metrics_splits = analyze_cohort_split_expression(combined_df, "Combined-Melanoma", output_dir)
        stats_list.extend(metrics_splits)
        
        # Save summary CSV
        df_stats = pd.DataFrame(stats_list)
        stats_path = output_dir / "tgm6_expression_stats.csv"
        df_stats.to_csv(stats_path, index=False)
        print(f"\nSaved combined stats summary to {stats_path}")


if __name__ == "__main__":
    main()
