#!/usr/bin/env python3
"""
Deconvolution Outlier Detection & Multi-Test Response Sensitivity Analysis
for LM22_no_IG Cell Type Fractions across ALL 10 iAtlas Datasets & Combined Cohorts.

1. Includes ALL 10 iAtlas datasets:
   - Hugo-iAtlas, Riaz-iAtlas, Liu-iAtlas, Gide-iAtlas (Melanoma)
   - Rosenberg-iAtlas (Bladder Cancer BLCA)
   - Padron-iAtlas (Pancreatic Cancer PAAD)
   - Anders-iAtlas (Breast Cancer BRCA)
   - McDermott-iAtlas, Choueiri-iAtlas (Renal Cell Carcinoma RCC)
   - Cloughesy-iAtlas (Glioblastoma GBM)

2. Implements 3 complementary outlier detection methods using cohort-relative 1.5 * IQR fences:
   - Method 1: Reconstruction Quality (rho_recon < Q1 - 1.5 * IQR)
   - Method 2: CLR Robust Mahalanobis Distance (D_M^2 > Q3 + 1.5 * IQR)
   - Method 3: Single-State Dominance Purity Index (Purity > max(0.60, Q3 + 1.5 * IQR))

3. Generates Outlier Diagnostic Visualizations per cohort & multi-cohort breakdown plots:
   - Multi-panel diagnostic dashboard (rho_recon, D_M^2, Purity, 2D CLR PCA scatter)
   - Overview breakdown stacked bar plot across all cohorts

4. Generates Actual Logistic Regression Curve Plots for every cell type and cohort:
   - X-axis: Cell type fraction
   - Y-axis: Binary response variable (0 = Non-Responder NR, 1 = Responder R)
   - Fitted Sigmoid Curve P(Responder) with 95% Confidence Interval band
   - Jittered scatter points for observations

5. Evaluates 8 combinatorial filtering strategies and re-runs 4 statistical hypothesis tests:
   - Likelihood Ratio Test (LRT)
   - Wald Test (for Odds Ratio & 95% CI)
   - Pearson Linear Correlation (r, p-value)
   - Spearman Rank Correlation (rho, p-value)
"""

import sys
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2, norm, pearsonr, spearmanr
from sklearn.covariance import MinCovDet
from sklearn.decomposition import PCA
import statsmodels.api as sm

# Ensure local packages are in path
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))
import datalair
import ici_datasets

# Set publication plotting style
sns.set_theme(style="whitegrid")
plt.rcParams.update({
    'font.size': 10,
    'figure.autolayout': True,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9
})

RESPONSE_COLORS = {
    'R': '#2ca02c',   # Green for Responders
    'NR': '#d62728'   # Red for Non-Responders
}

OUTLIER_COLORS = {
    'Retained': '#2ca02c',            # Green
    'Recon_Only': '#d62728',          # Red
    'CLR_Only': '#1f77b4',            # Blue
    'Purity_Only': '#ff7f0e',         # Orange
    'Multiple_Flags': '#9467bd'       # Purple
}


def is_ig_gene(symbol: str) -> bool:
    """Check if a gene symbol represents an Immunoglobulin heavy/light chain or J-chain transcript."""
    s = str(symbol).strip().upper()
    return s.startswith(('IGH', 'IGK', 'IGL', 'IGLL')) or s in ['JCHAIN', 'MZB1']


def load_clinical_response(cohort_name: str, lair):
    """Helper to load clinical response labels from CBioPortal datasets."""
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    try:
        data_dir = ici_datasets.cbioportal_datasets.get_dataset_dir(lair, dataset_class, cohort_name)
        filepath_clinical_sample = data_dir / "data_clinical_sample.txt"
        if not filepath_clinical_sample.exists():
            return None
            
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
        print(f"  Error loading clinical response for {cohort_name}: {e}")
        return None


def compute_outlier_metrics(df_fractions: pd.DataFrame, cohort_name: str, lair, lm22_df: pd.DataFrame):
    """
    Computes 3 outlier metrics for each sample in a cohort:
    1. rho_recon: Log1p Pearson correlation between bulk TPM and reconstructed TPM (\sum_k p_ik S_k)
    2. mahalanobis_clr: Robust Mahalanobis distance D_M^2 in Centered Log-Ratio (CLR) space
    3. purity_index: Single-state dominance purity index \sum_k p_ik^2
    """
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    data_dir = ici_datasets.cbioportal_datasets.get_dataset_dir(lair, dataset_class, cohort_name)
    
    df_bulk_tpm = None
    for fname in ['data_mrna_seq_expression.txt', 'data_mrna_seq_tpm.txt', 'data_mrna_seq_rpkm.txt']:
        p_file = data_dir / fname
        if p_file.exists():
            df_b = pd.read_csv(p_file, sep='\t', index_col=0)
            df_b.index = df_b.index.astype(str).str.upper()
            if df_b.max().max() < 50:
                df_b = np.power(2, df_b) - 1
            df_bulk_tpm = df_b.div(df_b.sum(axis=0), axis=1) * 1e6
            break
            
    lm22_norm = lm22_df / lm22_df.sum(axis=0)
    common_genes = df_bulk_tpm.index.intersection(lm22_norm.index) if df_bulk_tpm is not None else []
    common_genes = [g for g in common_genes if not is_ig_gene(g)]
    
    sample_ids = list(df_fractions.index)
    N = len(sample_ids)
    K = df_fractions.shape[1]
    
    # 1. Method 1: Bulk Reconstruction Log1p Correlation (rho_recon)
    rho_recon_dict = {}
    if df_bulk_tpm is not None and len(common_genes) > 10:
        bulk_sub = df_bulk_tpm.loc[common_genes]
        ref_sub = lm22_norm.loc[common_genes]
        ref_mat = ref_sub.values # (G, K)
        
        for sid in sample_ids:
            if sid in bulk_sub.columns:
                y_t = np.log1p(bulk_sub[sid].values)
                y_p = np.log1p(ref_mat @ df_fractions.loc[sid].values * 1e4)
                r, _ = pearsonr(y_t, y_p)
                rho_recon_dict[sid] = float(r) if not np.isnan(r) else 1.0
            else:
                rho_recon_dict[sid] = 1.0
    else:
        for sid in sample_ids:
            rho_recon_dict[sid] = 1.0
            
    # 2. Method 2: Centered Log-Ratio (CLR) & Robust Mahalanobis Distance
    eps = 1e-4
    P_mat = df_fractions.values + eps
    geom_mean = np.exp(np.mean(np.log(P_mat), axis=1, keepdims=True))
    CLR_mat = np.log(P_mat / geom_mean)
    
    mahalanobis_dict = {}
    try:
        if N > K + 2:
            mcd = MinCovDet(random_state=42).fit(CLR_mat)
            d2 = mcd.dist_
        else:
            cov = np.cov(CLR_mat, rowvar=False)
            inv_cov = np.linalg.pinv(cov)
            diff = CLR_mat - np.mean(CLR_mat, axis=0)
            d2 = np.sum((diff @ inv_cov) * diff, axis=1)
            
        for idx, sid in enumerate(sample_ids):
            mahalanobis_dict[sid] = float(d2[idx])
    except Exception:
        for sid in sample_ids:
            mahalanobis_dict[sid] = 0.0
            
    # 3. Method 3: Single-State Purity Index (\sum_k p_ik^2)
    purity_dict = {}
    for sid in sample_ids:
        p_vec = df_fractions.loc[sid].values
        purity_dict[sid] = float(np.sum(p_vec ** 2))
        
    df_outliers = pd.DataFrame({
        'Sample_ID': sample_ids,
        'rho_recon': [rho_recon_dict[s] for s in sample_ids],
        'mahalanobis_clr': [mahalanobis_dict[s] for s in sample_ids],
        'purity_index': [purity_dict[s] for s in sample_ids]
    }).set_index('Sample_ID')
    
    # 4. Cohort-Relative 1.5 * IQR Fence Outlier Flags
    q25_r, q75_r = np.percentile(df_outliers['rho_recon'], [25, 75])
    iqr_r = q75_r - q25_r
    low_fence_r = q25_r - 1.5 * iqr_r
    df_outliers['flag_recon'] = df_outliers['rho_recon'] < low_fence_r
    df_outliers['fence_recon'] = low_fence_r
    
    q25_m, q75_m = np.percentile(df_outliers['mahalanobis_clr'], [25, 75])
    iqr_m = q75_m - q25_m
    high_fence_m = q75_m + 1.5 * iqr_m
    df_outliers['flag_clr'] = df_outliers['mahalanobis_clr'] > high_fence_m
    df_outliers['fence_clr'] = high_fence_m
    
    q25_p, q75_p = np.percentile(df_outliers['purity_index'], [25, 75])
    iqr_p = q75_p - q25_p
    high_fence_p = max(0.60, q75_p + 1.5 * iqr_p)
    df_outliers['flag_purity'] = df_outliers['purity_index'] > high_fence_p
    df_outliers['fence_purity'] = high_fence_p
    
    # Category reason string
    reasons = []
    for sid in sample_ids:
        r = df_outliers.loc[sid, 'flag_recon']
        c = df_outliers.loc[sid, 'flag_clr']
        p = df_outliers.loc[sid, 'flag_purity']
        flags = [r, c, p]
        if sum(flags) == 0:
            reasons.append('Retained')
        elif sum(flags) > 1:
            reasons.append('Multiple_Flags')
        elif r:
            reasons.append('Recon_Only')
        elif c:
            reasons.append('CLR_Only')
        elif p:
            reasons.append('Purity_Only')
    df_outliers['outlier_category'] = reasons
    
    return df_outliers, CLR_mat


def generate_cohort_outlier_diagnostic_plot(df_fractions: pd.DataFrame, df_outliers: pd.DataFrame, CLR_mat: np.ndarray, cohort_name: str, output_path: Path):
    """Generates a 4-panel diagnostic dashboard for outlier detection in a cohort."""
    fig = plt.figure(figsize=(20, 12), dpi=300)
    
    low_fence_r = df_outliers['fence_recon'].iloc[0]
    high_fence_m = df_outliers['fence_clr'].iloc[0]
    high_fence_p = df_outliers['fence_purity'].iloc[0]
    
    # Panel A: Reconstruction Correlation Distribution
    ax1 = fig.add_subplot(2, 2, 1)
    sns.histplot(data=df_outliers, x='rho_recon', hue='flag_recon', palette={False: '#2ca02c', True: '#d62728'}, bins=25, kde=False, ax=ax1)
    ax1.axvline(low_fence_r, color='black', linestyle='--', linewidth=2, label=f'Fence: {low_fence_r:.3f}')
    ax1.set_title(f"A. Reconstruction Correlation (rho_recon)\nLower Fence (Q1 - 1.5*IQR) = {low_fence_r:.3f}", fontsize=11, fontweight='bold')
    ax1.set_xlabel("Log1p Reconstruction Correlation")
    ax1.set_ylabel("Sample Count")
    ax1.legend()
    
    # Panel B: CLR Robust Mahalanobis Distance Scatter
    ax2 = fig.add_subplot(2, 2, 2)
    sample_ranks = np.arange(len(df_outliers))
    colors_b = ['#d62728' if f else '#1f77b4' for f in df_outliers['flag_clr']]
    ax2.scatter(sample_ranks, df_outliers['mahalanobis_clr'], c=colors_b, alpha=0.8, edgecolors='none', s=40)
    ax2.axhline(high_fence_m, color='black', linestyle='--', linewidth=2, label=f'Fence: {high_fence_m:.2f}')
    ax2.set_title(f"B. CLR Simplex Mahalanobis Distance (D_M^2)\nUpper Fence (Q3 + 1.5*IQR) = {high_fence_m:.2f}", fontsize=11, fontweight='bold')
    ax2.set_xlabel("Sample Index")
    ax2.set_ylabel("Mahalanobis Distance D_M^2")
    ax2.legend()
    
    # Panel C: Single-State Dominance Purity Index
    ax3 = fig.add_subplot(2, 2, 3)
    sns.histplot(data=df_outliers, x='purity_index', hue='flag_purity', palette={False: '#2ca02c', True: '#ff7f0e'}, bins=25, kde=False, ax=ax3)
    ax3.axvline(high_fence_p, color='black', linestyle='--', linewidth=2, label=f'Fence: {high_fence_p:.3f}')
    ax3.set_title(f"C. Single-State Dominance Purity Index (Sum p_ik^2)\nUpper Fence = {high_fence_p:.3f}", fontsize=11, fontweight='bold')
    ax3.set_xlabel("Purity Index (\sum_k p_ik^2)")
    ax3.set_ylabel("Sample Count")
    ax3.legend()
    
    # Panel D: 2D Compositional PCA (CLR Simplex Space)
    ax4 = fig.add_subplot(2, 2, 4)
    pca = PCA(n_components=2, random_state=42)
    pcs = pca.fit_transform(CLR_mat)
    df_outliers['PC1'] = pcs[:, 0]
    df_outliers['PC2'] = pcs[:, 1]
    
    sns.scatterplot(
        data=df_outliers,
        x='PC1',
        y='PC2',
        hue='outlier_category',
        palette=OUTLIER_COLORS,
        style='outlier_category',
        s=60,
        alpha=0.9,
        ax=ax4
    )
    var1, var2 = pca.explained_variance_ratio_[:2] * 100
    ax4.set_title(f"D. 2D CLR Compositional PCA (Included vs Excluded Outlier Reasons)\nExpl. Var: PC1={var1:.1f}%, PC2={var2:.1f}%", fontsize=11, fontweight='bold')
    ax4.set_xlabel(f"PC1 ({var1:.1f}%)")
    ax4.set_ylabel(f"PC2 ({var2:.1f}%)")
    ax4.legend(title="Outlier Status", bbox_to_anchor=(1.02, 1), loc='upper left')
    
    n_ret = sum(df_outliers['outlier_category'] == 'Retained')
    n_ex = len(df_outliers) - n_ret
    fig.suptitle(
        f"Deconvolution Outlier Diagnostics & Inclusion Classification: {cohort_name}\n"
        f"Total Samples N={len(df_outliers)} | Included/Retained = {n_ret} ({n_ret/len(df_outliers)*100:.1f}%), Excluded = {n_ex} ({n_ex/len(df_outliers)*100:.1f}%)",
        fontsize=15, fontweight='bold', y=1.02
    )
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    fig.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
    plt.close(fig)


def generate_logistic_regression_curve_grid(df_cohort: pd.DataFrame, df_sub_res: pd.DataFrame, cell_types: list, c_name: str, combo_name: str, output_path: Path):
    """
    Generates a multi-panel grid of actual logistic regression curves:
    - X-axis: Cell type fraction
    - Y-axis: Binary immunotherapy response (0 = Non-Responder NR, 1 = Responder R)
    - Data points: Jittered scatter points for NR (red) and R (green)
    - Fitted Sigmoid Curve: P(Responder) = 1 / (1 + exp(-(beta_0 + beta_1 * x)))
    - 95% Confidence Interval Band around the fitted probability curve
    """
    n_ct = len(cell_types)
    n_cols = 4
    n_rows = int(np.ceil(n_ct / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(24, 4.8 * n_rows), dpi=300)
    axes_flat = axes.flatten()
    
    Y = (df_cohort['response'] == 'R').astype(int).values
    
    for idx, ct in enumerate(cell_types):
        ax = axes_flat[idx]
        X = df_cohort[ct].values
        
        mask = ~np.isnan(X) & ~np.isnan(Y)
        X_val, Y_val = X[mask], Y[mask]
        
        # Scatter points with jitter for visibility
        np.random.seed(42)
        jitter_y = Y_val + np.random.uniform(-0.03, 0.03, size=len(Y_val))
        colors_pts = [RESPONSE_COLORS['R'] if y == 1 else RESPONSE_COLORS['NR'] for y in Y_val]
        
        ax.scatter(X_val, jitter_y, c=colors_pts, alpha=0.6, edgecolors='none', s=25, zorder=2)
        
        # Fit Logit and plot Sigmoid Curve + 95% CI
        x_max = max(0.02, float(np.percentile(X_val, 99) * 1.15)) if len(X_val) > 0 else 1.0
        x_grid = np.linspace(0, x_max, 200)
        
        if len(np.unique(Y_val)) == 2 and len(X_val) >= 10:
            try:
                X_full = sm.add_constant(X_val)
                full_mod = sm.Logit(Y_val, X_full)
                full_res = full_mod.fit(disp=False, maxiter=100)
                
                X_grid_design = np.column_stack([np.ones(len(x_grid)), x_grid])
                eta_grid = X_grid_design @ full_res.params
                
                cov_beta = full_res.cov_params()
                se_eta = np.sqrt(np.sum((X_grid_design @ cov_beta) * X_grid_design, axis=1))
                
                p_grid = 1.0 / (1.0 + np.exp(-eta_grid))
                p_low = 1.0 / (1.0 + np.exp(-(eta_grid - 1.96 * se_eta)))
                p_high = 1.0 / (1.0 + np.exp(-(eta_grid + 1.96 * se_eta)))
                
                ax.plot(x_grid, p_grid, color='#1f77b4', linewidth=2.5, label='P(Responder)', zorder=3)
                ax.fill_between(x_grid, p_low, p_high, color='#1f77b4', alpha=0.2, zorder=1)
            except Exception:
                pass
                
        ax.axhline(0.5, color='gray', linestyle=':', linewidth=1.0, alpha=0.7)
        ax.axhline(0.0, color='gray', linestyle='-', linewidth=0.5)
        ax.axhline(1.0, color='gray', linestyle='-', linewidth=0.5)
        
        row_info = df_sub_res.loc[ct] if ct in df_sub_res.index else None
        if row_info is not None and not np.isnan(row_info['p_val_lrt']):
            p_lrt = row_info['p_val_lrt']
            p_wald = row_info['p_val_wald']
            or_val = row_info['or']
            r_p = row_info['pearson_r']
            p_p = row_info['p_val_pearson']
            rho_s = row_info['spearman_rho']
            p_s = row_info['p_val_spearman']
            
            ast_lrt = p_value_to_asterisks(p_lrt)
            lrt_str = f"LRT p={format_p_val(p_lrt)} ({ast_lrt})"
            wald_str = f"Wald p={format_p_val(p_wald)} | OR={format_or_val(or_val)}"
            pearson_str = f"Pearson r={r_p:.2f} (p={format_p_val(p_p)})"
            spearman_str = f"Spearman rho={rho_s:.2f} (p={format_p_val(p_s)})"
            
            is_sig = (p_lrt < 0.05 or p_wald < 0.05 or p_p < 0.05 or p_s < 0.05)
            title_color = 'red' if is_sig else 'black'
            
            ax.set_title(
                f"{ct}\n{lrt_str}\n{wald_str}\n{pearson_str}\n{spearman_str}",
                fontsize=8, fontweight='bold', color=title_color
            )
        else:
            ax.set_title(ct, fontsize=9, fontweight='bold')
            
        ax.set_xlabel("Cell Type Fraction", fontsize=8)
        ax.set_ylabel("P(Responder)", fontsize=8)
        ax.set_xlim(-0.005, x_max)
        ax.set_ylim(-0.08, 1.08)
        ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(['NR (0)', '0.25', '0.50', '0.75', 'R (1)'])
        
    for idx in range(n_ct, len(axes_flat)):
        fig.delaxes(axes_flat[idx])
        
    n_r_total = sum(df_cohort['response'] == 'R')
    n_nr_total = sum(df_cohort['response'] == 'NR')
    
    fig.suptitle(
        f"Logistic Regression Curves: P(Responder) vs Cell Fraction ({combo_name}): {c_name}\n"
        f"Total Samples N={len(df_cohort)} (Responders R={n_r_total}, Non-Responders NR={n_nr_total}) | Fitted Sigmoid Curve with 95% CI Band",
        fontsize=15, fontweight='bold', y=1.01
    )
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    fig.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
    plt.close(fig)


def run_statistical_tests(X: np.ndarray, Y: np.ndarray):
    """
    Fits logistic regression and computes:
    1. Likelihood Ratio Test (LRT)
    2. Wald Test for Odds Ratio (z-test)
    3. Pearson Correlation (r, p-value)
    4. Spearman Rank Correlation (rho, p-value)
    """
    N = len(Y)
    if N < 10 or len(np.unique(Y)) < 2:
        return {
            'lrt_stat': np.nan, 'p_val_lrt': np.nan,
            'wald_z': np.nan, 'p_val_wald': np.nan,
            'pearson_r': np.nan, 'p_val_pearson': np.nan,
            'spearman_rho': np.nan, 'p_val_spearman': np.nan,
            'beta': np.nan, 'or': np.nan, 'or_ci_low': np.nan, 'or_ci_high': np.nan,
            'n_r': int(np.sum(Y == 1)), 'n_nr': int(np.sum(Y == 0))
        }
        
    try:
        r_p, p_p = pearsonr(X, Y)
    except Exception:
        r_p, p_p = np.nan, np.nan
        
    try:
        rho_s, p_s = spearmanr(X, Y)
    except Exception:
        rho_s, p_s = np.nan, np.nan
        
    try:
        X_null = np.ones((N, 1))
        null_mod = sm.Logit(Y, X_null)
        null_res = null_mod.fit(disp=False, maxiter=100)
        llf_null = null_res.llf
        
        X_full = sm.add_constant(X)
        full_mod = sm.Logit(Y, X_full)
        full_res = full_mod.fit(disp=False, maxiter=100)
        llf_full = full_res.llf
        
        lrt_stat = max(0.0, 2.0 * (llf_full - llf_null))
        p_val_lrt = float(chi2.sf(lrt_stat, df=1))
        
        beta = float(full_res.params[1])
        se = float(full_res.bse[1])
        
        wald_z = float(beta / se) if se > 0 else np.nan
        p_val_wald = float(2.0 * (1.0 - norm.cdf(abs(wald_z)))) if not np.isnan(wald_z) else np.nan
        
        odds_ratio = float(np.exp(beta))
        or_ci_low = float(np.exp(beta - 1.96 * se))
        or_ci_high = float(np.exp(beta + 1.96 * se))
        
        return {
            'lrt_stat': lrt_stat,
            'p_val_lrt': p_val_lrt,
            'wald_z': wald_z,
            'p_val_wald': p_val_wald,
            'pearson_r': float(r_p),
            'p_val_pearson': float(p_p),
            'spearman_rho': float(rho_s),
            'p_val_spearman': float(p_s),
            'beta': beta,
            'or': odds_ratio,
            'or_ci_low': or_ci_low,
            'or_ci_high': or_ci_high,
            'n_r': int(np.sum(Y == 1)),
            'n_nr': int(np.sum(Y == 0))
        }
    except Exception:
        return {
            'lrt_stat': np.nan, 'p_val_lrt': np.nan,
            'wald_z': np.nan, 'p_val_wald': np.nan,
            'pearson_r': float(r_p) if not np.isnan(r_p) else np.nan,
            'p_val_pearson': float(p_p) if not np.isnan(p_p) else np.nan,
            'spearman_rho': float(rho_s) if not np.isnan(rho_s) else np.nan,
            'p_val_spearman': float(p_s) if not np.isnan(p_s) else np.nan,
            'beta': np.nan, 'or': np.nan, 'or_ci_low': np.nan, 'or_ci_high': np.nan,
            'n_r': int(np.sum(Y == 1)), 'n_nr': int(np.sum(Y == 0))
        }


def format_or_val(val: float) -> str:
    if np.isnan(val) or np.isinf(val):
        return "N/A"
    elif val >= 1000.0 or (val <= 0.001 and val > 0):
        return f"{val:.2e}"
    else:
        return f"{val:.2f}"


def format_p_val(p: float) -> str:
    if np.isnan(p):
        return "N/A"
    elif p < 0.0001:
        return f"{p:.2e}"
    else:
        return f"{p:.4f}"


def compute_candle_limit(vals_r: np.ndarray, vals_nr: np.ndarray) -> float:
    max_candle = 0.05
    for vals in [vals_r, vals_nr]:
        if len(vals) == 0:
            continue
        q25, q75 = np.percentile(vals, [25, 75])
        iqr = q75 - q25
        upper_fence = q75 + 1.5 * iqr
        valid_whiskers = vals[vals <= upper_fence]
        candle = np.max(valid_whiskers) if len(valid_whiskers) > 0 else q75
        if candle > max_candle:
            max_candle = candle
            
    return min(1.02, float(max_candle * 1.25))


def p_value_to_asterisks(p: float) -> str:
    if np.isnan(p):
        return ""
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def main():
    lair_path = "/storage/halu/lair"
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    project_root = Path(__file__).resolve().parent.parent
    deconv_dir = project_root / "output" / "cohort-comparisons"
    output_dir = project_root / "output" / "outlier-filtered-response"
    output_dir.mkdir(parents=True, exist_ok=True)
    diag_dir = output_dir / "diagnostic_plots"
    diag_dir.mkdir(parents=True, exist_ok=True)
    
    lm22_src = Path("/storage/halu/manual-download/LM22/LM22.txt")
    lm22_df = pd.read_csv(lm22_src, sep="\t", index_col=0)
    lm22_df.index = lm22_df.index.str.upper()
    
    # ALL 10 iAtlas Datasets
    iatlas_cohorts = [
        "Hugo-iAtlas",
        "Riaz-iAtlas",
        "Liu-iAtlas",
        "Gide-iAtlas",
        "Rosenberg-iAtlas",
        "Padron-iAtlas",
        "Anders-iAtlas",
        "McDermott-iAtlas",
        "Choueiri-iAtlas",
        "Cloughesy-iAtlas"
    ]
    
    iatlas_cancer_map = {
        'Combined-Melanoma': ['Hugo-iAtlas', 'Riaz-iAtlas', 'Liu-iAtlas', 'Gide-iAtlas'],
        'Combined-BLCA': ['Rosenberg-iAtlas'],
        'Combined-PAAD': ['Padron-iAtlas'],
        'Combined-BRCA': ['Anders-iAtlas'],
        'Combined-RCC': ['McDermott-iAtlas', 'Choueiri-iAtlas'],
        'Combined-GBM': ['Cloughesy-iAtlas'],
        'Combined-All-iAtlas': iatlas_cohorts
    }
    
    # 1. Load Clinical Response and Fractions for LM22_no_IG
    print("\nLoading clinical response labels and fractions for ALL 10 iAtlas datasets...")
    clinical_map = {c: load_clinical_response(c, lair) for c in iatlas_cohorts if load_clinical_response(c, lair) is not None}
    fractions_map = {c: pd.read_csv(deconv_dir / f"deconv_{c}_LM22_no_IG.csv", index_col=0) for c in iatlas_cohorts if (deconv_dir / f"deconv_{c}_LM22_no_IG.csv").exists()}
    
    for c in iatlas_cohorts:
        n_c = len(clinical_map[c]) if c in clinical_map else 0
        n_f = len(fractions_map[c]) if c in fractions_map else 0
        print(f"  {c:<20}: Clinical Samples={n_c}, Fractions Samples={n_f}")
        
    # 2. Compute Outlier Metrics & Generate Diagnostic Plots per cohort
    print("\nComputing outlier metrics & generating diagnostic visualization plots...")
    outlier_map = {}
    summary_rows = []
    
    for cohort in iatlas_cohorts:
        if cohort in fractions_map:
            df_out, CLR_mat = compute_outlier_metrics(fractions_map[cohort], cohort, lair, lm22_df)
            outlier_map[cohort] = df_out
            df_out.to_csv(output_dir / f"outlier_diagnostics_{cohort}.csv")
            
            # Generate 4-panel diagnostic figure
            diag_fig_path = diag_dir / f"outlier_diagnostics_{cohort}.png"
            generate_cohort_outlier_diagnostic_plot(fractions_map[cohort], df_out, CLR_mat, cohort, diag_fig_path)
            
            counts = df_out['outlier_category'].value_counts()
            summary_rows.append({
                'Cohort': cohort,
                'Total_Samples': len(df_out),
                'Retained': counts.get('Retained', 0),
                'Recon_Only': counts.get('Recon_Only', 0),
                'CLR_Only': counts.get('CLR_Only', 0),
                'Purity_Only': counts.get('Purity_Only', 0),
                'Multiple_Flags': counts.get('Multiple_Flags', 0),
                'Total_Excluded': len(df_out) - counts.get('Retained', 0)
            })
            print(f"  {cohort:<20}: Retained={counts.get('Retained', 0)} / {len(df_out)} | Excluded={len(df_out) - counts.get('Retained', 0)}")
            
    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(output_dir / "outlier_breakdown_all_cohorts_summary.csv", index=False)
    
    # 3. Overview Outlier Breakdown Stacked Bar Chart across ALL cohorts
    plt.figure(figsize=(14, 7), dpi=300)
    df_plot = df_summary.set_index('Cohort')[['Retained', 'Recon_Only', 'CLR_Only', 'Purity_Only', 'Multiple_Flags']]
    df_plot.plot(kind='bar', stacked=True, color=[OUTLIER_COLORS[c] for c in df_plot.columns], figsize=(14, 7), ax=plt.gca())
    plt.title("Deconvolution Outlier Category Breakdown Across ALL 10 iAtlas Cohorts", fontsize=14, fontweight='bold', y=1.02)
    plt.xlabel("iAtlas Cohorts", fontweight='bold')
    plt.ylabel("Sample Count")
    plt.xticks(rotation=45, ha='right', fontweight='bold')
    plt.legend(title="Outlier Status", bbox_to_anchor=(1.02, 1), loc='upper left')
    
    bar_chart_path = output_dir / "outlier_breakdown_all_cohorts.png"
    plt.savefig(bar_chart_path, bbox_inches='tight')
    plt.savefig(output_dir / "outlier_breakdown_all_cohorts.svg", bbox_inches='tight')
    plt.close()
    print(f"\nSaved overview outlier breakdown bar chart to {bar_chart_path}")

    # 4. Define 8 Combinatorial Outlier Filtering Strategies
    filter_combos = {
        'Unfiltered': lambda df: df,
        'Filter_Reconstruction': lambda df: df[~df['flag_recon']],
        'Filter_CLR_Mahalanobis': lambda df: df[~df['flag_clr']],
        'Filter_Purity': lambda df: df[~df['flag_purity']],
        'Filter_Recon_and_CLR': lambda df: df[~(df['flag_recon'] | df['flag_clr'])],
        'Filter_Recon_and_Purity': lambda df: df[~(df['flag_recon'] | df['flag_purity'])],
        'Filter_CLR_and_Purity': lambda df: df[~(df['flag_clr'] | df['flag_purity'])],
        'Filter_Composite_All': lambda df: df[~(df['flag_recon'] | df['flag_clr'] | df['flag_purity'])]
    }

    for combo_name, filter_func in filter_combos.items():
        print(f"\n=======================================================")
        print(f"EVALUATING OUTLIER FILTERING STRATEGY: {combo_name}")
        print(f"=======================================================")
        
        combo_dir = output_dir / combo_name
        combo_dir.mkdir(parents=True, exist_ok=True)
        boxplots_dir = combo_dir / "cohort_boxplots"
        boxplots_dir.mkdir(parents=True, exist_ok=True)
        curves_dir = combo_dir / "logistic_curves"
        curves_dir.mkdir(parents=True, exist_ok=True)
        
        data_merged_map = {}
        for cohort in iatlas_cohorts:
            if cohort in fractions_map and cohort in clinical_map and cohort in outlier_map:
                df_f = fractions_map[cohort]
                df_c = clinical_map[cohort]
                df_o = outlier_map[cohort]
                
                valid_samples = filter_func(df_o).index
                common = df_f.index.intersection(df_c.index).intersection(valid_samples)
                
                if len(common) > 0:
                    df_m = df_f.loc[common].copy()
                    df_m['response'] = df_c.loc[common, 'response']
                    data_merged_map[cohort] = df_m
                    
        for comb_name, sub_cohorts in iatlas_cancer_map.items():
            sub_dfs = [data_merged_map[c] for c in sub_cohorts if c in data_merged_map]
            if sub_dfs:
                data_merged_map[comb_name] = pd.concat(sub_dfs, axis=0)
                
        if not data_merged_map:
            print(f"  Warning: No data remaining for strategy {combo_name}. Skipping.")
            continue
            
        ordered_eval_cohorts = [
            "Hugo-iAtlas", "Riaz-iAtlas", "Liu-iAtlas", "Gide-iAtlas", "Combined-Melanoma",
            "Rosenberg-iAtlas", "Combined-BLCA",
            "Padron-iAtlas", "Combined-PAAD",
            "Anders-iAtlas", "Combined-BRCA",
            "McDermott-iAtlas", "Choueiri-iAtlas", "Combined-RCC",
            "Cloughesy-iAtlas", "Combined-GBM",
            "Combined-All-iAtlas"
        ]
        eval_cohorts = [c for c in ordered_eval_cohorts if c in data_merged_map]
        
        sample_df = next(iter(data_merged_map.values()))
        cell_types = [c for c in sample_df.columns if c != 'response']
        
        results_list = []
        for c_name in eval_cohorts:
            df_cohort = data_merged_map[c_name]
            Y = (df_cohort['response'] == 'R').astype(int).values
            
            for ct in cell_types:
                X = df_cohort[ct].values
                res_dict = run_statistical_tests(X, Y)
                res_dict['Filter_Strategy'] = combo_name
                res_dict['Cohort'] = c_name
                res_dict['CellType'] = ct
                results_list.append(res_dict)
                
        df_results = pd.DataFrame(results_list)
        csv_path = combo_dir / f"lrt_summary_table_{combo_name}.csv"
        df_results.to_csv(csv_path, index=False)
        print(f"  Saved statistical summary table to {csv_path}")
        
        # OVERVIEW HEATMAP
        p_val_pivot = df_results.pivot(index='CellType', columns='Cohort', values='p_val_lrt')[eval_cohorts]
        or_pivot = df_results.pivot(index='CellType', columns='Cohort', values='or')[eval_cohorts]
        beta_pivot = df_results.pivot(index='CellType', columns='Cohort', values='beta')[eval_cohorts]
        
        annot_matrix = np.empty(p_val_pivot.shape, dtype=object)
        for i in range(p_val_pivot.shape[0]):
            for j in range(p_val_pivot.shape[1]):
                p_v = p_val_pivot.iloc[i, j]
                or_v = or_pivot.iloc[i, j]
                if np.isnan(p_v) or np.isnan(or_v):
                    annot_matrix[i, j] = "N/A"
                else:
                    ast = p_value_to_asterisks(p_v)
                    annot_matrix[i, j] = f"{format_or_val(or_v)}\n({ast})"
                    
        plt.figure(figsize=(19, 11), dpi=300)
        sns.heatmap(
            beta_pivot,
            cmap="vlag",
            center=0,
            annot=annot_matrix,
            fmt="",
            cbar_kws={'label': 'Log Odds Ratio (beta > 0 favors Responder R)'},
            linewidths=0.8,
            linecolor="white",
            annot_kws={"size": 7.5, "fontweight": "bold"}
        )
        plt.title(f"Likelihood Ratio Test (LRT) & Odds Ratios for Response: {combo_name}\nLogit(P(Responder)) = beta_0 + beta_1 * Cell_Fraction", fontsize=14, fontweight='bold', y=1.02)
        plt.xlabel("iAtlas & Combined Cohorts", fontweight='bold')
        plt.ylabel("LM22 Cell Types (No IG Genes)", fontweight='bold')
        plt.xticks(rotation=45, ha='right', fontweight='bold')
        
        heatmap_path = combo_dir / f"lrt_heatmap_{combo_name}.png"
        plt.savefig(heatmap_path, bbox_inches='tight')
        plt.savefig(combo_dir / f"lrt_heatmap_{combo_name}.svg", bbox_inches='tight')
        plt.close()
        
        # MULTI-PANEL RESPONDER BOXPLOTS & LOGISTIC REGRESSION CURVES
        for c_name in eval_cohorts:
            df_cohort = data_merged_map[c_name]
            df_sub_res = df_results[df_results['Cohort'] == c_name].set_index('CellType')
            
            n_ct = len(cell_types)
            n_cols = 4
            n_rows = int(np.ceil(n_ct / n_cols))
            
            # A. Multi-Panel Boxplots
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(24, 4.6 * n_rows), dpi=300)
            axes_flat = axes.flatten()
            
            for idx, ct in enumerate(cell_types):
                ax = axes_flat[idx]
                df_ct = df_cohort[[ct, 'response']].dropna()
                
                vals_r = df_ct[df_ct['response'] == 'R'][ct].values
                vals_nr = df_ct[df_ct['response'] == 'NR'][ct].values
                
                y_limit = compute_candle_limit(vals_r, vals_nr)
                
                sns.boxplot(
                    data=df_ct,
                    x='response',
                    y=ct,
                    hue='response',
                    order=['NR', 'R'],
                    palette=RESPONSE_COLORS,
                    fliersize=2,
                    linewidth=1.0,
                    ax=ax,
                    legend=False
                )
                
                row_info = df_sub_res.loc[ct] if ct in df_sub_res.index else None
                if row_info is not None and not np.isnan(row_info['p_val_lrt']):
                    p_lrt = row_info['p_val_lrt']
                    p_wald = row_info['p_val_wald']
                    or_val = row_info['or']
                    r_p = row_info['pearson_r']
                    p_p = row_info['p_val_pearson']
                    rho_s = row_info['spearman_rho']
                    p_s = row_info['p_val_spearman']
                    
                    ast_lrt = p_value_to_asterisks(p_lrt)
                    
                    lrt_str = f"LRT p={format_p_val(p_lrt)} ({ast_lrt})"
                    wald_str = f"Wald p={format_p_val(p_wald)} | OR={format_or_val(or_val)}"
                    pearson_str = f"Pearson r={r_p:.2f} (p={format_p_val(p_p)})"
                    spearman_str = f"Spearman rho={rho_s:.2f} (p={format_p_val(p_s)})"
                    
                    is_sig = (p_lrt < 0.05 or p_wald < 0.05 or p_p < 0.05 or p_s < 0.05)
                    title_color = 'red' if is_sig else 'black'
                    
                    ax.set_title(
                        f"{ct}\n{lrt_str}\n{wald_str}\n{pearson_str}\n{spearman_str}",
                        fontsize=8, fontweight='bold', color=title_color
                    )
                else:
                    ax.set_title(ct, fontsize=9, fontweight='bold')
                    
                ax.set_xlabel("")
                ax.set_ylabel("Fraction", fontsize=8)
                ax.set_ylim(-0.01, y_limit)
                
            for idx in range(n_ct, len(axes_flat)):
                fig.delaxes(axes_flat[idx])
                
            n_r_total = sum(df_cohort['response'] == 'R')
            n_nr_total = sum(df_cohort['response'] == 'NR')
            
            fig.suptitle(
                f"Cell Type Fraction Response Comparisons ({combo_name}): {c_name}\n"
                f"Total Samples N={len(df_cohort)} (Responders R={n_r_total}, Non-Responders NR={n_nr_total}) | Y-limits scaled by highest candle",
                fontsize=15, fontweight='bold', y=1.01
            )
            fig.tight_layout()
            
            cohort_clean = c_name.lower().replace("-", "_")
            out_box_path = boxplots_dir / f"response_boxplots_{cohort_clean}_{combo_name}.png"
            fig.savefig(out_box_path, bbox_inches='tight')
            fig.savefig(boxplots_dir / f"response_boxplots_{cohort_clean}_{combo_name}.svg", bbox_inches='tight')
            plt.close(fig)

            # B. Actual Logistic Regression Curves Grid
            out_curve_path = curves_dir / f"logistic_curves_{cohort_clean}_{combo_name}.png"
            generate_logistic_regression_curve_grid(df_cohort, df_sub_res, cell_types, c_name, combo_name, out_curve_path)

    print("\nOutlier detection, response sensitivity analysis, and logistic regression curve generation completed successfully for ALL 10 iAtlas datasets.")


if __name__ == "__main__":
    main()
