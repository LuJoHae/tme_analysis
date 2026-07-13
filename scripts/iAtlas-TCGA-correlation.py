#!/usr/bin/env python3
"""
Analyze correlation between inferred cell type fractions in iAtlas clinical trials
and their corresponding TCGA baseline cohorts.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

# Ensure local packages are in path
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))
import datalair
from single_cell_datasets._single_cell_reference import SingleCellDeconvolution


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
    
    resolutions = ["leiden_res_0.5", "kmeans_subcluster_res_0.5"]
    
    correlation_rows = []
    
    for resolution in resolutions:
        print(f"\n--- Analyzing correlations for {resolution} ---")
        
        fig, axes = plt.subplots(3, 3, figsize=(18, 16), dpi=300)
        axes_flat = axes.flatten()
        
        for idx, (trial, tcga_proj) in enumerate(COHORT_TO_TCGA.items()):
            ax = axes_flat[idx]
            
            trial_key = f"deconv_{trial}_{resolution}.csv"
            tcga_key = f"deconv_TCGA-{tcga_proj}_{resolution}.csv"
            
            trial_path = deconv_paths.get(trial_key)
            tcga_path = deconv_paths.get(tcga_key)
            
            if not trial_path or not tcga_path:
                ax.text(0.5, 0.5, f"Missing data for {trial} / TCGA-{tcga_proj}", 
                        ha='center', va='center', fontsize=10, color='gray')
                continue
                
            # Load dataframes
            df_trial = pd.read_csv(trial_path, index_col=0)
            df_tcga = pd.read_csv(tcga_path, index_col=0)
            
            # Compute mean cell type fractions across samples
            mean_trial = df_trial.mean(axis=0)
            mean_tcga = df_tcga.mean(axis=0)
            
            # Align on common cell types
            common_cts = mean_trial.index.intersection(mean_tcga.index)
            y_trial = mean_trial.loc[common_cts].values
            x_tcga = mean_tcga.loc[common_cts].values
            
            # Calculate Pearson and Spearman correlations
            pears_r, pears_p = pearsonr(x_tcga, y_trial)
            spear_r, spear_p = spearmanr(x_tcga, y_trial)
            
            # Save results
            correlation_rows.append({
                'Resolution': resolution,
                'Trial_Cohort': trial,
                'TCGA_Cohort': f"TCGA-{tcga_proj}",
                'Pearson_r': pears_r,
                'Pearson_p': pears_p,
                'Spearman_rho': spear_r,
                'Spearman_p': spear_p,
                'Num_Cell_Types': len(common_cts)
            })
            
            # Scatter plot
            sns.regplot(
                x=x_tcga, y=y_trial, ax=ax,
                scatter_kws={'s': 50, 'alpha': 0.7, 'color': '#8E44AD', 'edgecolors': 'black', 'linewidths': 0.6},
                line_kws={'color': '#E74C3C', 'linewidth': 1.5, 'linestyle': '--'}
            )
            
            # Add identity line
            max_val = max(x_tcga.max(), y_trial.max()) * 1.15
            ax.plot([0, max_val], [0, max_val], color='#BDC3C7', linestyle=':', linewidth=1)
            
            # Label points with cell type labels for key outliers
            for ct, x_val, y_val in zip(common_cts, x_tcga, y_trial):
                # Annotate top cell types
                if x_val > 0.05 or y_val > 0.05:
                    ax.annotate(str(ct), (x_val, y_val), textcoords="offset points", 
                                xytext=(0, 5), ha='center', fontsize=7, weight='bold', color='#2C3E50')
            
            ax.set_xlim(-0.02, max_val)
            ax.set_ylim(-0.02, max_val)
            
            # Add statistics text box
            stats_text = (
                f"Pearson $r = {pears_r:.3f}$ ($p = {pears_p:.1e}$)\n"
                f"Spearman $\\rho = {spear_r:.3f}$ ($p = {spear_p:.1e}$)"
            )
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                    verticalalignment='top', fontsize=9, bbox=dict(boxstyle='round', facecolor='white', alpha=0.85, edgecolor='#BDC3C7'))
            
            ax.set_title(f"{trial} vs TCGA-{tcga_proj}", fontsize=11, fontweight='bold', color='#2C3E50')
            ax.set_xlabel(f"TCGA-{tcga_proj} Mean Fractions", fontsize=9)
            ax.set_ylabel(f"{trial} Mean Fractions", fontsize=9)
            ax.grid(True, linestyle='--', alpha=0.5)
            
        plt.suptitle(f"Cell Type Fraction Correlation: iAtlas Trials vs TCGA Baseline ({resolution})", 
                     fontsize=15, fontweight='bold', y=0.995)
        plt.tight_layout()
        plt.savefig(output_dir / f"deconv_tcga_correlation_{resolution}.svg", format='svg', bbox_inches='tight')
        plt.close()
        print(f"Saved correlation plot for {resolution}.")
        
    # Save statistics table
    df_corrs = pd.DataFrame(correlation_rows)
    corr_csv_path = output_dir / "deconv_tcga_correlation_summary.csv"
    df_corrs.to_csv(corr_csv_path, index=False)
    print(f"Saved correlation summary table to {corr_csv_path.name}")


if __name__ == "__main__":
    main()
