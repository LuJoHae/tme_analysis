#!/usr/bin/env python3
"""
Develop and test methods to estimate the best VAF threshold and model TMB reliability
using genomic mutations and bulk RNA sequencing data.
Compares 3 unsupervised VAF threshold estimation methods and evaluates Expressed RNA-TMB.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, pearsonr, spearmanr, gaussian_kde
from scipy.signal import find_peaks
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc
from sklearn.mixture import GaussianMixture
import warnings
from sklearn.exceptions import ConvergenceWarning

warnings.filterwarnings('ignore', category=ConvergenceWarning)

# Ensure local packages are in python path if running directly
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))

import datalair
import ici_datasets
from gene_utils import calculate_maf_tmb


def otsu_threshold_1d(vafs):
    """
    Otsu's thresholding method adapted for 1D VAF values.
    Finds the threshold that minimizes intra-class variance.
    """
    vafs = np.sort(vafs)
    n = len(vafs)
    if n < 5:
        return np.nan
    best_thresh = np.nan
    min_variance = np.inf
    
    # Evaluate split points
    for i in range(2, n - 2, max(1, n // 200)):
        c0 = vafs[:i]
        c1 = vafs[i:]
        w0 = i / n
        w1 = (n - i) / n
        
        var0 = np.var(c0)
        var1 = np.var(c1)
        
        intra_class_var = w0 * var0 + w1 * var1
        if intra_class_var < min_variance:
            min_variance = intra_class_var
            best_thresh = (vafs[i-1] + vafs[i]) / 2.0
            
    return best_thresh


def gmm_threshold(vafs):
    """
    Gaussian Mixture Model thresholding.
    Fits a 2-component GMM and finds the decision boundary (equal posterior probabilities).
    """
    if len(vafs) < 10:
        return np.nan
    try:
        vafs_2d = np.array(vafs).reshape(-1, 1)
        gmm = GaussianMixture(n_components=2, random_state=42)
        gmm.fit(vafs_2d)
        
        means = gmm.means_.flatten()
        idx_high = np.argmax(means)
        
        # Grid search to find decision boundary in [0.01, 0.60]
        test_vafs = np.linspace(0.01, 0.60, 1000).reshape(-1, 1)
        probs = gmm.predict_proba(test_vafs)
        
        diff = np.abs(probs[:, idx_high] - 0.5)
        best_idx = np.argmin(diff)
        return test_vafs[best_idx][0]
    except Exception:
        return np.nan


def find_kde_valley_and_peak(vafs):
    """
    Kernel Density Estimation valley detection.
    Fits a KDE in [0.02, 0.60], finds the clonal peak, and detects the first valley to its left.
    """
    if len(vafs) < 10:
        return np.nan, np.nan
        
    try:
        # Fit KDE
        kde = gaussian_kde(vafs)
        x_grid = np.linspace(0.02, 0.60, 500)
        y_grid = kde(x_grid)
        
        # Find peaks
        peaks, _ = find_peaks(y_grid, distance=20)
        if len(peaks) == 0:
            peak_x = x_grid[np.argmax(y_grid)]
            return peak_x / 2.0, peak_x
            
        # Major peak (highest density in clonal range)
        major_peak_idx = peaks[np.argmax(y_grid[peaks])]
        peak_x = x_grid[major_peak_idx]
        
        # Find first valley (local minimum) to the left of the peak
        y_left = y_grid[:major_peak_idx]
        x_left = x_grid[:major_peak_idx]
        
        valleys, _ = find_peaks(-y_left, distance=10)
        if len(valleys) > 0:
            valley_x = x_left[valleys[-1]]  # Rightmost valley to the left of peak
        else:
            valley_x = peak_x / 2.0
            
        return valley_x, peak_x
    except Exception:
        return np.nan, np.nan


def get_cohort_data(lair, ds_name):
    """
    Load mutations, clinical response, and mRNA expression for a single cohort.
    Deduplicates mRNA index and columns to prevent duplicate key series conversion.
    """
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    try:
        data_dir = ici_datasets.cbioportal_datasets.get_dataset_dir(lair, dataset_class, ds_name)
        df_clinical = ici_datasets.cbioportal_datasets.load_data_clinical(data_dir)
        df_clinical = df_clinical[df_clinical['response'].isin(['R', 'NR'])]
        df_clinical.index = df_clinical.index.astype(str)
        
        mut_file = data_dir / "data_mutations.txt"
        if not mut_file.exists():
            return None, None, None
            
        df_mut = pd.read_csv(mut_file, sep="\t", low_memory=False)
        df_mut['Tumor_Sample_Barcode'] = df_mut['Tumor_Sample_Barcode'].astype(str)
        
        # Load mRNA expression
        _, df_mrna = ici_datasets.cbioportal_datasets.load_and_process_data(data_dir)
        df_mrna.columns = df_mrna.columns.astype(str)
        df_mrna = df_mrna.loc[~df_mrna.index.duplicated(keep='first')]
        df_mrna = df_mrna.loc[:, ~df_mrna.columns.duplicated(keep='first')]
        
        return df_mut, df_clinical, df_mrna
    except Exception as e:
        print(f"Error loading cohort {ds_name}: {e}")
        return None, None, None


def main():
    lair_path = "/storage/halu/lair"
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    # Resolve output directory
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output" / "iAtlas-vaf-tmb-reliability"
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
    
    # Stromal microenvironment genes for biological purity validation
    stromal_genes = ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'VIM', 'FN1']
    
    cohort_stats = []
    vaf_density_curves = []
    improvement_stats = []
    sample_tmb_distributions = []
    tmb_correlation_stats = []
    paired_tmb_samples = []
    
    print("\n--- Processing Cohorts & Estimating VAF Thresholds ---")
    for cohort in individual_cohorts:
        df_mut, df_clin, df_mrna = get_cohort_data(lair, cohort)
        if df_mut is None or df_clin is None or df_mrna is None:
            continue
            
        n_samples = len(df_clin)
        r_count = sum(df_clin['response'] == 'R')
        nr_count = sum(df_clin['response'] == 'NR')
        
        if n_samples < 10 or r_count < 3 or nr_count < 3:
            continue
            
        # Calculate mutation depth and VAF values
        df_mut = df_mut[df_mut['t_depth'] > 0].copy()
        df_mut['VAF'] = df_mut['t_alt_count'] / df_mut['t_depth']
        
        # Align mutations with clinical samples
        df_mut_filtered = df_mut[df_mut['Tumor_Sample_Barcode'].isin(df_clin.index)].copy()
        vafs_all = df_mut_filtered['VAF'].values
        
        # Filter out very low VAF noise for threshold estimators
        vafs_for_est = vafs_all[(vafs_all >= 0.02) & (vafs_all <= 0.80)]
        
        # 1. Run VAF threshold estimators
        valley_val, peak_val = find_kde_valley_and_peak(vafs_for_est)
        otsu_val = otsu_threshold_1d(vafs_for_est)
        gmm_val = gmm_threshold(vafs_for_est)
        
        # Save KDE curves for plotting
        if len(vafs_for_est) > 10:
            kde = gaussian_kde(vafs_for_est)
            x_grid = np.linspace(0.0, 0.60, 100)
            y_grid = kde(x_grid)
            for x, y in zip(x_grid, y_grid):
                vaf_density_curves.append({
                    'Cohort': cohort,
                    'VAF': x,
                    'Density': y
                })
                
        # 2. Integrate RNA-seq data for reliability
        avail_stromal = [g for g in stromal_genes if g in df_mrna.index]
        if avail_stromal:
            microenv_scores = df_mrna.loc[avail_stromal].mean(axis=0)
            median_microenv = microenv_scores.reindex(df_clin.index, fill_value=0.0).median()
        else:
            median_microenv = np.nan
            
        # Mutation expression rate: somatic mutations (VAF >= 0.05) in expressed genes (TPM >= 1.0)
        df_somatic = df_mut_filtered[df_mut_filtered['VAF'] >= 0.05].copy()
        expressed_count = 0
        total_somatic = len(df_somatic)
        
        for _, row in df_somatic.iterrows():
            gene = row['Hugo_Symbol']
            sample = row['Tumor_Sample_Barcode']
            if gene in df_mrna.index and sample in df_mrna.columns:
                if df_mrna.loc[gene, sample] >= 1.0:
                    expressed_count += 1
                    
        mer = expressed_count / total_somatic if total_somatic > 0 else np.nan
        median_depth = df_mut_filtered['t_depth'].median()
        
        # TMB Reliability Score (TRS)
        trs = mer * np.log10(median_depth) if not pd.isna(mer) and median_depth > 0 else np.nan
        
        # 3. Compute True Optimal VAF Threshold & Max TMB ROC AUC from clinical labels
        grid_vafs = np.arange(0.01, 0.51, 0.01)
        best_pval = 1.0
        best_vaf_true = 0.05
        best_auc = 0.5
        
        y_encoded = (df_clin['response'] == 'R').astype(int)
        
        for thresh in grid_vafs:
            tmb_matrix = calculate_maf_tmb(df_mut, capture_size_mb=30.0, vaf_threshold=thresh)
            tmb_matrix = tmb_matrix.set_index('Tumor_Sample_Barcode')
            tmb_scores = tmb_matrix['TMB_Score'].reindex(df_clin.index, fill_value=0.0)
            
            tmb_r = tmb_scores[df_clin['response'] == 'R']
            tmb_nr = tmb_scores[df_clin['response'] == 'NR']
            
            if len(tmb_r) > 1 and len(tmb_nr) > 1:
                _, pval = mannwhitneyu(tmb_r, tmb_nr, alternative='two-sided')
                auc_score = roc_auc_score(y_encoded, tmb_scores)
            else:
                pval = 1.0
                auc_score = 0.5
                
            if pval < best_pval:
                best_pval = pval
                best_vaf_true = thresh
                best_auc = auc_score
                
        # Save sample TMB values at the True Optimal VAF threshold
        opt_tmb_matrix = calculate_maf_tmb(df_mut, capture_size_mb=30.0, vaf_threshold=best_vaf_true)
        opt_tmb_matrix = opt_tmb_matrix.set_index('Tumor_Sample_Barcode')
        opt_tmb_scores = opt_tmb_matrix['TMB_Score'].reindex(df_clin.index, fill_value=0.0)
        
        for sample, tmb_val in opt_tmb_scores.items():
            sample_tmb_distributions.append({
                'Cohort': cohort,
                'Sample': sample,
                'Response': df_clin.loc[sample, 'response'],
                'TMB_at_Optimal_VAF': tmb_val,
                'Optimal_VAF_Threshold': best_vaf_true
            })
                
        # Estimated Purity from clonal peak
        est_purity = 2 * peak_val if not pd.isna(peak_val) else np.nan
        
        cohort_stats.append({
            'Cohort': cohort,
            'Est_Optimal_Threshold_KDE': valley_val,
            'Est_Optimal_Threshold_Otsu': otsu_val,
            'Est_Optimal_Threshold_GMM': gmm_val,
            'True_Optimal_Threshold': best_vaf_true,
            'Est_Clonal_Peak': peak_val,
            'Est_Purity': est_purity,
            'Stromal_Score': median_microenv,
            'Mutation_Expression_Rate': mer,
            'Median_Sequencing_Depth': median_depth,
            'TMB_Reliability_Score': trs,
            'Max_TMB_ROC_AUC': best_auc,
            'Max_TMB_p_value': best_pval
        })
        
        # 4. Compare DNA-TMB vs Expressed RNA-TMB (DNA VAF >= 0.05 and RNA TPM >= 1.0)
        # DNA-TMB
        dna_tmb_matrix = calculate_maf_tmb(df_mut, capture_size_mb=30.0, vaf_threshold=0.05)
        dna_tmb_matrix = dna_tmb_matrix.set_index('Tumor_Sample_Barcode')
        dna_tmb = dna_tmb_matrix['TMB_Score'].reindex(df_clin.index, fill_value=0.0)
        
        # RNA-TMB (Only mutations in expressed genes)
        df_expressed_somatic = df_somatic.copy()
        df_expressed_somatic['Is_Expressed'] = 0
        for idx, row in df_expressed_somatic.iterrows():
            gene = row['Hugo_Symbol']
            sample = row['Tumor_Sample_Barcode']
            if gene in df_mrna.index and sample in df_mrna.columns:
                if df_mrna.loc[gene, sample] >= 1.0:
                    df_expressed_somatic.at[idx, 'Is_Expressed'] = 1
                    
        df_expressed_only = df_expressed_somatic[df_expressed_somatic['Is_Expressed'] == 1]
        
        # Re-run TMB calculation using only expressed somatic mutations
        if len(df_expressed_only) > 0:
            rna_tmb_matrix = calculate_maf_tmb(df_expressed_only, capture_size_mb=30.0, vaf_threshold=0.0)
            rna_tmb_matrix = rna_tmb_matrix.set_index('Tumor_Sample_Barcode')
            rna_tmb = rna_tmb_matrix['TMB_Score'].reindex(df_clin.index, fill_value=0.0)
        else:
            rna_tmb = pd.Series(0.0, index=df_clin.index)
            
        # Statistical comparisons
        dna_tmb_r = dna_tmb[df_clin['response'] == 'R']
        dna_tmb_nr = dna_tmb[df_clin['response'] == 'NR']
        if len(dna_tmb_r) > 1 and len(dna_tmb_nr) > 1:
            _, dna_pval = mannwhitneyu(dna_tmb_r, dna_tmb_nr, alternative='two-sided')
            dna_auc = roc_auc_score(y_encoded, dna_tmb)
            dna_precision, dna_recall, _ = precision_recall_curve(y_encoded, dna_tmb)
            dna_pr_auc = auc(dna_recall, dna_precision)
        else:
            dna_pval, dna_auc, dna_pr_auc = np.nan, np.nan, np.nan
            
        rna_tmb_r = rna_tmb[df_clin['response'] == 'R']
        rna_tmb_nr = rna_tmb[df_clin['response'] == 'NR']
        if len(rna_tmb_r) > 1 and len(rna_tmb_nr) > 1:
            _, rna_pval = mannwhitneyu(rna_tmb_r, rna_tmb_nr, alternative='two-sided')
            rna_auc = roc_auc_score(y_encoded, rna_tmb)
            rna_precision, rna_recall, _ = precision_recall_curve(y_encoded, rna_tmb)
            rna_pr_auc = auc(rna_recall, rna_precision)
        else:
            rna_pval, rna_auc, rna_pr_auc = np.nan, np.nan, np.nan
            
        improvement_stats.append({
            'Cohort': cohort,
            'DNA_TMB_ROC_AUC': dna_auc,
            'RNA_TMB_ROC_AUC': rna_auc,
            'DNA_TMB_PR_AUC': dna_pr_auc,
            'RNA_TMB_PR_AUC': rna_pr_auc,
            'DNA_TMB_p_value': dna_pval,
            'RNA_TMB_p_value': rna_pval
        })
        
        # DNA-TMB vs Expressed RNA-TMB correlation analysis
        if len(dna_tmb) > 2:
            p_corr, p_pval = pearsonr(dna_tmb, rna_tmb)
            s_corr, s_pval = spearmanr(dna_tmb, rna_tmb)
        else:
            p_corr, p_pval, s_corr, s_pval = np.nan, np.nan, np.nan, np.nan
            
        tmb_correlation_stats.append({
            'Cohort': cohort,
            'Pearson_r': p_corr,
            'Pearson_p': p_pval,
            'Spearman_rho': s_corr,
            'Spearman_p': s_pval
        })
        
        for sample in df_clin.index:
            paired_tmb_samples.append({
                'Cohort': cohort,
                'Sample': sample,
                'Response': df_clin.loc[sample, 'response'],
                'DNA_TMB': dna_tmb.loc[sample],
                'RNA_TMB': rna_tmb.loc[sample]
            })
        
        print(f"  {cohort} processed. KDE Valley: {valley_val:.2f}, Otsu: {otsu_val:.2f}, GMM: {gmm_val:.2f}, True Optimal: {best_vaf_true:.2f}")

    df_cohort_stats = pd.DataFrame(cohort_stats)
    df_curves = pd.DataFrame(vaf_density_curves)
    df_improvement = pd.DataFrame(improvement_stats)
    df_sample_tmb = pd.DataFrame(sample_tmb_distributions)
    df_tmb_corr = pd.DataFrame(tmb_correlation_stats)
    df_paired_samples = pd.DataFrame(paired_tmb_samples)
    
    # Save files
    stats_csv_path = output_dir / "cohort_reliability_stats.csv"
    curves_csv_path = output_dir / "vaf_density_curves.csv"
    improvement_csv_path = output_dir / "tmb_expression_improvement_comparison.csv"
    sample_tmb_csv_path = output_dir / "cohort_sample_tmb_distributions.csv"
    tmb_corr_csv_path = output_dir / "dna_vs_rna_tmb_correlation_stats.csv"
    paired_samples_csv_path = output_dir / "dna_vs_rna_tmb_samples.csv"
    
    df_cohort_stats.to_csv(stats_csv_path, index=False)
    df_curves.to_csv(curves_csv_path, index=False)
    df_improvement.to_csv(improvement_csv_path, index=False)
    df_sample_tmb.to_csv(sample_tmb_csv_path, index=False)
    df_tmb_corr.to_csv(tmb_corr_csv_path, index=False)
    df_paired_samples.to_csv(paired_samples_csv_path, index=False)
    
    print(f"\nSaved cohort stats to {stats_csv_path.name}")
    print(f"Saved VAF density curves to {curves_csv_path.name}")
    print(f"Saved improvement comparison stats to {improvement_csv_path.name}")
    print(f"Saved sample TMB distributions to {sample_tmb_csv_path.name}")
    print(f"Saved DNA vs RNA TMB correlation stats to {tmb_corr_csv_path.name}")
    print(f"Saved paired DNA/RNA TMB samples to {paired_samples_csv_path.name}")

    # --- PLOT 1: VAF Densities with Peaks and Valleys ---
    plt.figure(figsize=(14, 9), dpi=300)
    for idx, cohort in enumerate(df_cohort_stats['Cohort'].unique()):
        plt.subplot(3, 3, idx + 1)
        sub_c = df_curves[df_curves['Cohort'] == cohort]
        stats_c = df_cohort_stats[df_cohort_stats['Cohort'] == cohort].iloc[0]
        
        plt.plot(sub_c['VAF'], sub_c['Density'], color='#34495E', linewidth=1.8)
        
        # Draw estimated thresholds
        plt.axvline(x=stats_c['Est_Optimal_Threshold_KDE'], color='#E74C3C', linestyle='-', label=f"KDE Valley ({stats_c['Est_Optimal_Threshold_KDE']:.2f})")
        plt.axvline(x=stats_c['Est_Optimal_Threshold_Otsu'], color='#F1C40F', linestyle='--', label=f"Otsu ({stats_c['Est_Optimal_Threshold_Otsu']:.2f})")
        plt.axvline(x=stats_c['Est_Optimal_Threshold_GMM'], color='#2ECC71', linestyle=':', label=f"GMM ({stats_c['Est_Optimal_Threshold_GMM']:.2f})")
        
        plt.title(cohort, fontsize=9.5, fontweight='bold', color='#2C3E50')
        plt.xlabel("VAF", fontsize=8)
        plt.ylabel("Density", fontsize=8)
        plt.xlim(0.0, 0.6)
        plt.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
        plt.legend(fontsize=7, loc='upper right')
        
    plt.suptitle("DNA-VAF Density Valley Detection & Threshold Comparison", fontsize=13, fontweight='bold', color='#2C3E50', y=0.995)
    plt.tight_layout()
    plt.savefig(output_dir / "vaf_density_valley_detection.svg", format='svg', bbox_inches='tight')
    plt.close()

    # --- PLOT 2: Compare Estimated vs. True VAF Thresholds ---
    plt.figure(figsize=(8, 6), dpi=300)
    plt.scatter(df_cohort_stats['Est_Optimal_Threshold_KDE'], df_cohort_stats['True_Optimal_Threshold'], 
                color='#E74C3C', s=70, label='KDE Valley', alpha=0.8, edgecolors='black')
    plt.scatter(df_cohort_stats['Est_Optimal_Threshold_Otsu'], df_cohort_stats['True_Optimal_Threshold'], 
                color='#F1C40F', s=70, label='Otsu Method', alpha=0.8, edgecolors='black')
    plt.scatter(df_cohort_stats['Est_Optimal_Threshold_GMM'], df_cohort_stats['True_Optimal_Threshold'], 
                color='#2ECC71', s=70, label='GMM Decision Boundary', alpha=0.8, edgecolors='black')
    plt.plot([0, 0.5], [0, 0.5], color='#7F8C8D', linestyle='--', linewidth=1.2, label='Identity Line')
    
    clean_stats = df_cohort_stats.dropna(subset=['Est_Optimal_Threshold_KDE', 'True_Optimal_Threshold'])
    if len(clean_stats) > 2:
        r_val, p_val = pearsonr(clean_stats['Est_Optimal_Threshold_KDE'], clean_stats['True_Optimal_Threshold'])
        plt.text(0.05, 0.45, f"KDE Valley vs True:\nPearson r = {r_val:.2f}\np = {p_val:.3f}", 
                 fontsize=9, fontweight='bold', color='#E74C3C')
                 
    plt.title("Comparison of Unsupervised Threshold Methods vs. True Optimal", fontsize=11, fontweight='bold', color='#2C3E50', pad=15)
    plt.xlabel("Estimated Optimal VAF Threshold", fontsize=9.5)
    plt.ylabel("True Optimal VAF Threshold (Clinical Response)", fontsize=9.5)
    plt.xlim(-0.02, 0.52)
    plt.ylim(-0.02, 0.52)
    plt.legend(loc='lower right', fontsize=8.5)
    plt.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir / "estimated_vs_true_thresholds.svg", format='svg', bbox_inches='tight')
    plt.close()

    # --- PLOT 3: TMB Reliability Score vs. Absolute ROC AUC Deviation (|AUC - 0.5|) ---
    plt.figure(figsize=(8, 6), dpi=300)
    
    # Calculate absolute deviation from chance
    df_cohort_stats['Abs_AUC_Dev'] = (df_cohort_stats['Max_TMB_ROC_AUC'] - 0.5).abs()
    
    # Only keep cohorts with active RNA-seq (TRS > 0)
    df_rel_plot = df_cohort_stats[df_cohort_stats['TMB_Reliability_Score'] > 0].dropna(subset=['TMB_Reliability_Score', 'Abs_AUC_Dev'])
    
    plt.scatter(df_rel_plot['TMB_Reliability_Score'], df_rel_plot['Abs_AUC_Dev'], 
                color='#3498DB', s=80, alpha=0.8, edgecolors='black')
                
    for _, row in df_rel_plot.iterrows():
        plt.text(row['TMB_Reliability_Score'] + 0.01, row['Abs_AUC_Dev'], row['Cohort'].replace('-iAtlas', ''),
                 va='center', ha='left', fontsize=8, fontweight='semibold', color='#34495E')
                 
    if len(df_rel_plot) > 2:
        r_val, p_val = pearsonr(df_rel_plot['TMB_Reliability_Score'], df_rel_plot['Abs_AUC_Dev'])
        sns.regplot(x='TMB_Reliability_Score', y='Abs_AUC_Dev', data=df_rel_plot, 
                    scatter=False, color='#2980B9', line_kws={'linestyle':'--', 'linewidth':1.2})
        plt.text(df_rel_plot['TMB_Reliability_Score'].min() + 0.02, df_rel_plot['Abs_AUC_Dev'].max() - 0.05, 
                 f"Pearson r = {r_val:.3f}\np = {p_val:.3f}", 
                 fontsize=9.5, fontweight='bold', color='#2980B9')
                 
    plt.title("TRS vs. Absolute Predictive Strength of TMB", fontsize=11, fontweight='bold', color='#2C3E50', pad=15)
    plt.xlabel("TMB Reliability Score (TRS)", fontsize=9.5)
    plt.ylabel("Absolute TMB ROC AUC Deviation: |AUC - 0.5|", fontsize=9.5)
    plt.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir / "reliability_score_vs_performance.svg", format='svg', bbox_inches='tight')
    plt.close()

    # --- PLOT 4: Purity Validation ---
    plt.figure(figsize=(8, 6), dpi=300)
    plt.scatter(df_cohort_stats['Est_Purity'], df_cohort_stats['Stromal_Score'], 
                color='#9B59B6', s=80, alpha=0.8, edgecolors='black')
                
    for _, row in df_cohort_stats.iterrows():
        if not pd.isna(row['Est_Purity']) and not pd.isna(row['Stromal_Score']):
            plt.text(row['Est_Purity'] + 0.01, row['Stromal_Score'], row['Cohort'].replace('-iAtlas', ''),
                     va='center', ha='left', fontsize=8, fontweight='semibold', color='#34495E')
                     
    clean_stats_pur = df_cohort_stats.dropna(subset=['Est_Purity', 'Stromal_Score'])
    if len(clean_stats_pur) > 2:
        r_val, p_val = pearsonr(clean_stats_pur['Est_Purity'], clean_stats_pur['Stromal_Score'])
        sns.regplot(x='Est_Purity', y='Stromal_Score', data=clean_stats_pur, 
                    scatter=False, color='#8E44AD', line_kws={'linestyle':'--', 'linewidth':1.2})
        plt.text(df_cohort_stats['Est_Purity'].min() + 0.01, df_cohort_stats['Stromal_Score'].max() - 0.2, 
                 f"Pearson r = {r_val:.2f}\np = {p_val:.3f}", 
                 fontsize=9.5, fontweight='bold', color='#8E44AD')
                 
    plt.title("Biological Purity Validation: Est Purity vs. Stromal Score", fontsize=11, fontweight='bold', color='#2C3E50', pad=15)
    plt.xlabel("DNA-VAF Estimated Purity (2 * Clonal Peak)", fontsize=9.5)
    plt.ylabel("RNA Stromal Microenvironment Score", fontsize=9.5)
    plt.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir / "purity_vs_stromal_score.svg", format='svg', bbox_inches='tight')
    plt.close()

    # --- PLOT 5: DNA-TMB vs. RNA-TMB ROC AUC ---
    plt.figure(figsize=(10, 6), dpi=300)
    df_plot = df_improvement.dropna(subset=['DNA_TMB_ROC_AUC', 'RNA_TMB_ROC_AUC']).copy()
    
    df_melt = df_plot.melt(id_vars='Cohort', value_vars=['DNA_TMB_ROC_AUC', 'RNA_TMB_ROC_AUC'], 
                          var_name='TMB_Type', value_name='ROC_AUC')
    df_melt['TMB_Type'] = df_melt['TMB_Type'].replace({
        'DNA_TMB_ROC_AUC': 'Standard DNA-TMB (VAF >= 0.05)',
        'RNA_TMB_ROC_AUC': 'Expressed RNA-TMB (DNA >= 0.05 + RNA >= 1.0 TPM)'
    })
    
    sns.barplot(x='ROC_AUC', y='Cohort', hue='TMB_Type', data=df_melt, palette=['#7F8C8D', '#2ECC71'])
    plt.axvline(x=0.5, color='#7F8C8D', linestyle='--', linewidth=1.0)
    plt.title("Performance Comparison: Standard DNA-TMB vs. Expressed RNA-TMB", fontsize=12, fontweight='bold', color='#2C3E50')
    plt.xlabel("Out-of-Fold ROC AUC Score", fontsize=10)
    plt.ylabel("Cohort", fontsize=10)
    plt.xlim(0.3, 0.85)
    plt.legend(title="TMB Type", loc="lower right", fontsize=8.5)
    plt.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5, axis='x')
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_dir / "dna_vs_rna_tmb_performance.svg", format='svg', bbox_inches='tight')
    plt.close()

    # --- PLOT 6: Cohort TMB Distributions by Response (For Each Cohort Individually) ---
    plt.figure(figsize=(15, 10), dpi=300)
    for idx, cohort in enumerate(df_cohort_stats['Cohort'].unique()):
        plt.subplot(3, 3, idx + 1)
        sub_tmb = df_sample_tmb[df_sample_tmb['Cohort'] == cohort]
        stats_c = df_cohort_stats[df_cohort_stats['Cohort'] == cohort].iloc[0]
        
        # Draw box plot overlayed with strip plot
        sns.boxplot(x='Response', y='TMB_at_Optimal_VAF', data=sub_tmb, 
                    palette={'R': '#2ECC71', 'NR': '#E74C3C'}, width=0.5, fliersize=0)
        sns.stripplot(x='Response', y='TMB_at_Optimal_VAF', data=sub_tmb, 
                      color='#2C3E50', size=4, alpha=0.6, jitter=0.15)
                      
        plt.title(f"{cohort}\n(VAF={stats_c['True_Optimal_Threshold']:.2f}, TRS={stats_c['TMB_Reliability_Score']:.2f})", 
                  fontsize=9.5, fontweight='bold', color='#2C3E50')
        plt.xlabel("Response", fontsize=8.5)
        plt.ylabel("TMB Score", fontsize=8.5)
        plt.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5, axis='y')
        
    plt.suptitle("Clinical Separation: TMB at Optimal VAF Threshold by Immunotherapy Response", 
                 fontsize=13, fontweight='bold', color='#2C3E50', y=0.995)
    plt.tight_layout()
    plt.savefig(output_dir / "cohort_tmb_response_distributions.svg", format='svg', bbox_inches='tight')
    plt.close()
    
    # Helper to format p-values cleanly
    def format_p_val(pval):
        if pd.isna(pval):
            return "N/A"
        return f"{pval:.3f}" if pval >= 0.001 else f"{pval:.1e}"

    # --- PLOT 7: Correlation between DNA-TMB and Expressed RNA-TMB ---
    plt.figure(figsize=(15, 10), dpi=300)
    for idx, cohort in enumerate(df_cohort_stats['Cohort'].unique()):
        plt.subplot(3, 3, idx + 1)
        sub_samples = df_paired_samples[df_paired_samples['Cohort'] == cohort]
        corr_c = df_tmb_corr[df_tmb_corr['Cohort'] == cohort].iloc[0]
        
        # Scatter plot colored by response
        sns.scatterplot(
            x='DNA_TMB', y='RNA_TMB', hue='Response', data=sub_samples,
            palette={'R': '#2ECC71', 'NR': '#E74C3C'}, alpha=0.8, edgecolor='black', s=40
        )
        
        # Regression line
        if len(sub_samples) > 2:
            sns.regplot(
                x='DNA_TMB', y='RNA_TMB', data=sub_samples,
                scatter=False, color='#34495E', line_kws={'linestyle': '--', 'linewidth': 1.2}
            )
            
        pearson_p_str = format_p_val(corr_c['Pearson_p'])
        spearman_p_str = format_p_val(corr_c['Spearman_p'])
        
        plt.title(f"{cohort}\nPearson r={corr_c['Pearson_r']:.2f} (p={pearson_p_str})\nSpearman rho={corr_c['Spearman_rho']:.2f} (p={spearman_p_str})", 
                  fontsize=9, fontweight='bold', color='#2C3E50')
        plt.xlabel("Standard DNA-TMB (VAF >= 0.05)", fontsize=8)
        plt.ylabel("Expressed RNA-TMB", fontsize=8)
        plt.grid(True, color='#E5E5E5', linestyle='-', linewidth=0.5)
        plt.legend(fontsize=7, loc='upper left')
        
    plt.suptitle("Correlation: Standard DNA-TMB vs. Expressed RNA-TMB", fontsize=13, fontweight='bold', color='#2C3E50', y=0.995)
    plt.tight_layout()
    plt.savefig(output_dir / "dna_vs_rna_tmb_correlation.svg", format='svg', bbox_inches='tight')
    plt.close()
    
    print(f"\nSuccessfully generated cohort individual distributions, correlation plots, and saved CSV files!")


if __name__ == "__main__":
    main()
