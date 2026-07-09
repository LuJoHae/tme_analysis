#!/usr/bin/env python3
"""
Explore Variant Allele Frequencies (VAF) across iAtlas cohorts.
Calculates summary stats of VAF per cohort and per variant classification,
compares VAF in Responders vs. Non-Responders, and generates diagnostic plots.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# Ensure local packages are in python path if running directly
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))

import datalair
import ici_datasets


def load_cohort_mutations_vaf(lair, ds_name):
    """
    Load mutations and response labels for a single cohort, calculate VAF.
    """
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    try:
        data_dir = ici_datasets.cbioportal_datasets.get_dataset_dir(lair, dataset_class, ds_name)
        df_clinical = ici_datasets.cbioportal_datasets.load_data_clinical(data_dir)
        df_clinical = df_clinical[df_clinical['response'].isin(['R', 'NR'])]
        
        mut_file = data_dir / "data_mutations.txt"
        if not mut_file.exists():
            return None
            
        df_mut = pd.read_csv(mut_file, sep="\t", low_memory=False)
        
        # Filter for valid depth and counts
        df_mut = df_mut[df_mut['t_depth'] > 0].copy()
        df_mut['VAF'] = df_mut['t_alt_count'] / df_mut['t_depth']
        
        # Keep only mutations in clinical samples with response
        df_mut = df_mut[df_mut['Tumor_Sample_Barcode'].isin(df_clinical.index)].copy()
        
        # Merge response status
        df_mut['Response'] = df_mut['Tumor_Sample_Barcode'].map(df_clinical['response'])
        df_mut['Cohort'] = ds_name
        
        cols = ['Cohort', 'Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification', 't_alt_count', 't_depth', 'VAF', 'Response']
        return df_mut[cols]
    except Exception as e:
        print(f"Error loading mutations for {ds_name}: {e}")
        return None


def main():
    lair_path = "/storage/halu/lair"
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    # Resolve directories
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output" / "iAtlas-VAF"
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
    
    all_muts_list = []
    print("\n--- Loading mutations and calculating VAF across cohorts ---")
    for cohort in individual_cohorts:
        df_m = load_cohort_mutations_vaf(lair, cohort)
        if df_m is not None:
            print(f"  {cohort}: Loaded {len(df_m)} mutations across {df_m['Tumor_Sample_Barcode'].nunique()} samples.")
            all_muts_list.append(df_m)
            
    if not all_muts_list:
        print("No mutations could be loaded.")
        return
        
    df_all = pd.concat(all_muts_list, axis=0).reset_index(drop=True)
    
    # 1. Cohort-level summary stats
    print("\n--- Generating VAF summary statistics per cohort ---")
    cohort_stats = []
    for cohort in df_all['Cohort'].unique():
        sub = df_all[df_all['Cohort'] == cohort]
        cohort_stats.append({
            'Cohort': cohort,
            'Total_Mutations': len(sub),
            'Mean_VAF': sub['VAF'].mean(),
            'Median_VAF': sub['VAF'].median(),
            'SD_VAF': sub['VAF'].std(),
            'Q25_VAF': sub['VAF'].quantile(0.25),
            'Q75_VAF': sub['VAF'].quantile(0.75)
        })
    df_cohort_stats = pd.DataFrame(cohort_stats)
    print(df_cohort_stats.to_string(index=False))
    
    # 2. VAF by Variant Classification
    print("\n--- Generating VAF summary statistics by Variant Classification ---")
    class_stats = []
    for vc in df_all['Variant_Classification'].unique():
        sub = df_all[df_all['Variant_Classification'] == vc]
        class_stats.append({
            'Variant_Classification': vc,
            'Total_Mutations': len(sub),
            'Mean_VAF': sub['VAF'].mean(),
            'Median_VAF': sub['VAF'].median(),
            'SD_VAF': sub['VAF'].std()
        })
    df_class_stats = pd.DataFrame(class_stats)
    
    # 3. VAF Responders vs Non-Responders comparison
    print("\n--- Comparing VAF between Responders and Non-Responders ---")
    response_stats = []
    for cohort in df_all['Cohort'].unique():
        sub = df_all[df_all['Cohort'] == cohort]
        vaf_r = sub[sub['Response'] == 'R']['VAF']
        vaf_nr = sub[sub['Response'] == 'NR']['VAF']
        
        if len(vaf_r) > 5 and len(vaf_nr) > 5:
            stat, pval = mannwhitneyu(vaf_r, vaf_nr, alternative='two-sided')
        else:
            pval = np.nan
            
        response_stats.append({
            'Cohort': cohort,
            'R_Mutations': len(vaf_r),
            'R_Mean_VAF': vaf_r.mean(),
            'R_Median_VAF': vaf_r.median(),
            'NR_Mutations': len(vaf_nr),
            'NR_Mean_VAF': vaf_nr.mean(),
            'NR_Median_VAF': vaf_nr.median(),
            'Mann_Whitney_U_p_value': pval
        })
    df_response_stats = pd.DataFrame(response_stats)
    print(df_response_stats.to_string(index=False))
    
    # Save combined statistics to CSV files
    df_cohort_stats.to_csv(output_dir / "vaf_cohort_summary.csv", index=False)
    df_class_stats.to_csv(output_dir / "vaf_classification_summary.csv", index=False)
    df_response_stats.to_csv(output_dir / "vaf_response_comparison.csv", index=False)
    print(f"\nSaved statistical summaries to {output_dir}/")
    
    # --- PLOT 1: VAF Distribution across Cohorts ---
    plt.figure(figsize=(12, 6), dpi=300)
    sns.violinplot(x='Cohort', y='VAF', data=df_all, inner='quartile', palette='Set3', cut=0)
    plt.title("Somatic Mutation Variant Allele Frequency (VAF) across iAtlas Cohorts", fontsize=13, fontweight='bold', color='#2C3E50', pad=15)
    plt.xlabel("Cohort", fontsize=10, fontweight='semibold')
    plt.ylabel("Variant Allele Frequency (VAF)", fontsize=10, fontweight='semibold')
    plt.xticks(rotation=30, ha='right')
    plt.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5, axis='y')
    sns.despine(left=True, bottom=True)
    plt.tight_layout()
    plot_cohorts_path = output_dir / "vaf_distribution_cohorts.svg"
    plt.savefig(plot_cohorts_path, format='svg', bbox_inches='tight')
    plt.close()
    
    # --- PLOT 2: VAF by Variant Classification (Top 10 classes) ---
    plt.figure(figsize=(12, 6), dpi=300)
    top_vcs = df_all['Variant_Classification'].value_counts().head(10).index
    df_top_vc = df_all[df_all['Variant_Classification'].isin(top_vcs)]
    sns.boxplot(x='Variant_Classification', y='VAF', data=df_top_vc, palette='Set2', showfliers=False)
    plt.title("VAF Distribution by Variant Classification", fontsize=13, fontweight='bold', color='#2C3E50', pad=15)
    plt.xlabel("Variant Classification", fontsize=10, fontweight='semibold')
    plt.ylabel("Variant Allele Frequency (VAF)", fontsize=10, fontweight='semibold')
    plt.xticks(rotation=30, ha='right')
    plt.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5, axis='y')
    sns.despine(left=True, bottom=True)
    plt.tight_layout()
    plot_vc_path = output_dir / "vaf_by_classification.svg"
    plt.savefig(plot_vc_path, format='svg', bbox_inches='tight')
    plt.close()
    
    # --- PLOT 3: VAF Responders vs Non-Responders ---
    plt.figure(figsize=(14, 7), dpi=300)
    # Exclude cohorts with extremely small mutation counts if any
    sns.boxplot(x='Cohort', y='VAF', hue='Response', data=df_all, palette=['#2ECC71', '#E74C3C'], showfliers=False, width=0.6)
    plt.title("Mutation VAF in Responders (R) vs. Non-Responders (NR) across Cohorts", fontsize=14, fontweight='bold', color='#2C3E50', pad=15)
    plt.xlabel("Cohort", fontsize=10, fontweight='semibold')
    plt.ylabel("Variant Allele Frequency (VAF)", fontsize=10, fontweight='semibold')
    plt.xticks(rotation=30, ha='right')
    plt.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5, axis='y')
    
    # Annotate p-values on plot
    for i, cohort in enumerate(individual_cohorts):
        p_row = df_response_stats[df_response_stats['Cohort'] == cohort]
        if not p_row.empty:
            pval = p_row.iloc[0]['Mann_Whitney_U_p_value']
            if not pd.isna(pval):
                text_lbl = f"p={pval:.2g}" if pval >= 0.01 else f"p={pval:.1e}"
                # Plot above the cohort boxplots
                plt.text(i, 0.95, text_lbl, ha='center', va='bottom', fontsize=8.5, fontweight='bold', color='#2980B9')
                
    plt.ylim(-0.02, 1.05)
    plt.legend(title="Immunotherapy Response", loc="upper right")
    sns.despine(left=True, bottom=True)
    plt.tight_layout()
    plot_resp_path = output_dir / "vaf_by_response.svg"
    plt.savefig(plot_resp_path, format='svg', bbox_inches='tight')
    plt.close()
    
    print(f"Successfully generated all VAF exploration plots in {output_dir}/")


if __name__ == "__main__":
    main()
