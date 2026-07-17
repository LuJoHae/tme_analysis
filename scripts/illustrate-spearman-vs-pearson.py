#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm

def compute_stats(vals, y):
    # Mann-Whitney U test
    vals_r = vals[y == 1]
    vals_nr = vals[y == 0]
    try:
        mwu_stat, mwu_p = stats.mannwhitneyu(vals_r, vals_nr, alternative='two-sided')
    except Exception:
        mwu_p = np.nan
        
    # Likelihood Ratio Test for Logistic Regression
    try:
        X_null = np.ones((len(y), 1))
        model_null = sm.Logit(y, X_null).fit(disp=False)
        X_full = sm.add_constant(vals)
        model_full = sm.Logit(y, X_full).fit(disp=False)
        lr_stat = 2 * (model_full.llf - model_null.llf)
        lrt_p = stats.chi2.sf(lr_stat, df=1)
    except Exception:
        lrt_p = np.nan
        
    return mwu_p, lrt_p

def kernel_regression(x_train, y_train, x_eval, h=None):
    if h is None:
        std = np.std(x_train)
        if std == 0:
            std = 1.0
        n = len(x_train)
        # Increased bandwidth (doubled Silverman's rule of thumb baseline for smoother curves)
        h = 2.0 * 1.06 * std * (n ** (-0.2))
        if h == 0:
            h = 0.1
    # Nadaraya-Watson Gaussian Kernel Regression
    dists = x_eval[:, np.newaxis] - x_train[np.newaxis, :]
    weights = np.exp(-0.5 * (dists / h) ** 2)
    sum_weights = np.sum(weights, axis=1)
    sum_weights = np.where(sum_weights == 0, 1e-10, sum_weights)
    y_pred = np.sum(weights * y_train, axis=1) / sum_weights
    return y_pred

def fit_custom_model(vals, y, lambda_0=None, lambda_1=None):
    from scipy.optimize import minimize
    
    # 1. Define S_low and S_high dynamically based on 20th and 80th percentiles
    p20 = np.percentile(vals, 20)
    p80 = np.percentile(vals, 80)
    
    if p20 == p80:
        s_low_mask = vals <= np.median(vals)
        s_high_mask = vals > np.median(vals)
    else:
        s_low_mask = vals <= p20
        s_high_mask = vals >= p80
        
    y_low = y[s_low_mask]
    y_high = y[s_high_mask]
    
    hat_p0 = np.mean(y_low) if len(y_low) > 0 else np.mean(y)
    hat_p1 = np.mean(y_high) if len(y_high) > 0 else np.mean(y)
    
    # Clip empirical estimates away from absolute 0/1 to keep initial guess well-defined
    hat_p0 = np.clip(hat_p0, 1e-4, 1 - 1e-4)
    hat_p1 = np.clip(hat_p1, 1e-4, 1 - 1e-4)
    
    # 2. Objective Function: Negative Log Likelihood (MLE)
    def objective(params):
        p0, p1, theta = params
        
        # pi(x) = p0 + (p1 - p0) * (x^theta)
        pi_x = p0 + (p1 - p0) * np.power(vals, theta)
        pi_x = np.clip(pi_x, 1e-15, 1 - 1e-15)
        
        # Log-likelihood
        log_lik = np.sum(y * np.log(pi_x) + (1 - y) * np.log(1 - pi_x))
        
        # Minimize negative log likelihood
        return -log_lik
        
    delta = 1e-5
    bounds = [
        (delta, 1 - delta),
        (delta, 1 - delta),
        (0.01, 100.0)
    ]
    
    # Start optimization at empirical estimates
    initial_guess = [hat_p0, hat_p1, 1.0]
    
    res = minimize(objective, initial_guess, method='L-BFGS-B', bounds=bounds)
    p0_opt, p1_opt, theta_opt = res.x
    return p0_opt, p1_opt, theta_opt, hat_p0, hat_p1

# Ensure local packages are in path
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))
import datalair
from single_cell_datasets import SingleCellDeconvolution

def load_clinical_response(cohort_name: str, lair):
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
        return df_clinical
    except Exception as e:
        print(f"Error loading clinical for {cohort_name}: {e}")
        return None

def main():
    lair_path = "/storage/halu/lair"
    lair = datalair.Lair(lair_path)
    ds_deconv = SingleCellDeconvolution()
    deconv_paths = lair.get_dataset_filepaths(ds_deconv)
    
    cohorts = ["Hugo-iAtlas", "Riaz-iAtlas", "Liu-iAtlas", "Gide-iAtlas"]
    resolution = "kmeans_subcluster_res_0.5"
    
    clinical_data = {}
    for cohort in cohorts:
        df_clin = load_clinical_response(cohort, lair)
        if df_clin is not None:
            clinical_data[cohort] = df_clin
            
    fracs_list = []
    y_list = []
    
    for c in cohorts:
        key = f"deconv_{c}_{resolution}.csv"
        path = deconv_paths.get(key)
        df_clin = clinical_data.get(c)
        if path and df_clin is not None:
            df_fracs = pd.read_csv(path, index_col=0)
            df_fracs.index = df_fracs.index.astype(str)
            common_samples = df_fracs.index.intersection(df_clin.index)
            if len(common_samples) >= 10:
                fracs_list.append(df_fracs.loc[common_samples])
                y_list.append(df_clin.loc[common_samples, 'response'].map({'R': 1, 'NR': 0}))
                
    df_fracs_comb = pd.concat(fracs_list, axis=0)
    y_series_comb = pd.concat(y_list, axis=0)
    
    comparison_rows = []
    for cell_type in df_fracs_comb.columns:
        vals = df_fracs_comb[cell_type].values
        y = y_series_comb.values
        if np.var(vals) == 0:
            continue
            
        spearman_r, spearman_p = stats.spearmanr(vals, y)
        pearson_r, pearson_p = stats.pearsonr(vals, y)
        
        # Measure skewness and outliers in fractions
        skewness = stats.skew(vals)
        # Outliers using IQR rule
        q75, q25 = np.percentile(vals, [75 ,25])
        iqr = q75 - q25
        upper_bound = q75 + 1.5 * iqr if iqr > 0 else 0
        n_outliers = np.sum(vals > upper_bound) if upper_bound > 0 else 0
        
        comparison_rows.append({
            'CellType': cell_type,
            'Spearman_R': spearman_r,
            'Spearman_P': spearman_p,
            'Pearson_R': pearson_r,
            'Pearson_P': pearson_p,
            'P_Ratio': pearson_p / spearman_p if spearman_p > 0 else 1.0,
            'Skewness': skewness,
            'Num_Outliers': n_outliers,
            'Mean_Fraction': np.mean(vals)
        })
        
    df_comp = pd.DataFrame(comparison_rows)
    output_dir = Path("output/single-cell-exploration")
    output_dir.mkdir(parents=True, exist_ok=True)
    df_comp.to_csv(output_dir / "spearman_vs_pearson_melanoma_comparison.csv", index=False)
    
    print("\n--- Top Features where Spearman is significantly MORE significant than Pearson ---")
    df_discrepant = df_comp[df_comp['Spearman_P'] < 0.05].sort_values(by='P_Ratio', ascending=False)
    for idx, row in df_discrepant.head(5).iterrows():
        print(f"CellType: {row['CellType']:20s} | Spearman P: {row['Spearman_P']:.6f} | Pearson P: {row['Pearson_P']:.6f} | Ratio: {row['P_Ratio']:.1f} | Skewness: {row['Skewness']:.2f} | Outliers: {row['Num_Outliers']}")
        
    # Filter features where mean fraction is > 0.03
    df_filtered = df_comp[(df_comp['Mean_Fraction'] > 0.03) & (df_comp['Spearman_P'] < 0.05)]
    df_filtered_discrepant = df_filtered.sort_values(by='P_Ratio', ascending=False)
    
    if not df_filtered_discrepant.empty:
        chosen_row = df_filtered_discrepant.iloc[0]
    else:
        chosen_row = df_discrepant.iloc[0]
        
    top_cell_type = chosen_row['CellType']
    vals = df_fracs_comb[top_cell_type].values
    y = y_series_comb.values
    y_labels = pd.Series(y).map({0: 'Non-Responder', 1: 'Responder'})
    
    fig, axes = plt.subplots(2, 4, figsize=(24, 10), dpi=300)
    
    mwu_p, lrt_p = compute_stats(vals, y)
    
    # Plot 1: Box/Violin plot showing the distribution of the raw fractions in R vs NR
    sns.boxplot(x=y_labels, y=vals, ax=axes[0, 0], palette={'Non-Responder': '#EC7063', 'Responder': '#5DADE2'}, hue=y_labels, legend=False)
    sns.stripplot(x=y_labels, y=vals, color='black', alpha=0.3, jitter=0.2, ax=axes[0, 0])
    axes[0, 0].set_title(f"Raw Fractions of {top_cell_type}\n(Highly Skewed: Skew={chosen_row['Skewness']:.2f})\n[MWU p={mwu_p:.6f}]", fontsize=11, fontweight='bold')
    axes[0, 0].set_xlabel("Response status")
    axes[0, 0].set_ylabel("Inferred Fraction")
    # Set y-axis limit to 95th percentile
    ymax = np.percentile(vals, 95)
    if ymax > 0:
        axes[0, 0].set_ylim(0, ymax)
    
    # Plot 2: Box plot of ranks (illustrating what Spearman sees)
    ranks = stats.rankdata(vals)
    sns.boxplot(x=y_labels, y=ranks, ax=axes[1, 0], palette={'Non-Responder': '#EC7063', 'Responder': '#5DADE2'}, hue=y_labels, legend=False)
    sns.stripplot(x=y_labels, y=ranks, color='black', alpha=0.3, jitter=0.2, ax=axes[1, 0])
    axes[1, 0].set_title(f"Rank-Transformed Fractions of {top_cell_type}\n(Spearman Correlation Basis)\n[MWU p={mwu_p:.6f}]", fontsize=11, fontweight='bold')
    axes[1, 0].set_xlabel("Response status")
    axes[1, 0].set_ylabel("Rank")
    
    # Plot 3: Scatter plot of raw fraction vs binary response to show the leverage of outliers
    sns.regplot(x=vals, y=y, logistic=True, ax=axes[0, 1], scatter_kws={'alpha': 0.4}, line_kws={'color': '#E67E22'})
    axes[0, 1].axhline(y=np.mean(y), color='gray', linestyle='--', alpha=0.7, label=f'Baseline Rate ({np.mean(y):.3f})')
    # Add vertical percentile lines
    pct_25 = np.percentile(vals, 25)
    pct_50 = np.percentile(vals, 50)
    pct_75 = np.percentile(vals, 75)
    pct_95 = np.percentile(vals, 95)
    axes[0, 1].axvline(x=pct_25, color='#8E44AD', linestyle=':', alpha=0.6, label=f'25th Pct ({pct_25:.4f})')
    axes[0, 1].axvline(x=pct_50, color='#8E44AD', linestyle='--', alpha=0.6, label=f'Median ({pct_50:.4f})')
    axes[0, 1].axvline(x=pct_75, color='#8E44AD', linestyle='-.', alpha=0.6, label=f'75th Pct ({pct_75:.4f})')
    axes[0, 1].axvline(x=pct_95, color='#C0392B', linestyle='--', alpha=0.6, label=f'95th Pct ({pct_95:.4f})')
    # Fit & plot Gaussian Kernel Regression in raw space
    x_grid_raw = np.linspace(np.min(vals), np.max(vals), 300)
    y_kernel_raw = kernel_regression(vals, y, x_grid_raw)
    axes[0, 1].plot(x_grid_raw, y_kernel_raw, color='#1ABC9C', lw=2, linestyle=':', label='Gaussian Kernel Reg')
    axes[0, 1].set_title(f"Linear/Logistic Regression Leverage\nPearson R={chosen_row['Pearson_R']:.3f} (p={chosen_row['Pearson_P']:.4f})\nSpearman R={chosen_row['Spearman_R']:.3f} (p={chosen_row['Spearman_P']:.4f})\n[LRT p={lrt_p:.6f}]", fontsize=11, fontweight='bold')
    axes[0, 1].set_xlabel("Raw Fraction")
    axes[0, 1].set_ylabel("Binary Response (0/1)")
    axes[0, 1].legend(loc='best', fontsize=8)
    
    # Plot 4: Custom model fit
    p0_opt, p1_opt, theta_opt, hat_p0, hat_p1 = fit_custom_model(vals, y)
    axes[0, 2].scatter(vals, y, alpha=0.4, color='#8E44AD', label='Data')
    x_grid = np.linspace(0, 1.0, 300)
    y_pred = p0_opt + (p1_opt - p0_opt) * np.power(x_grid, theta_opt)
    axes[0, 2].plot(x_grid, y_pred, color='#E74C3C', lw=2.5, label='Custom Model')
    # Overlay Gaussian Kernel Regression
    axes[0, 2].plot(x_grid_raw, y_kernel_raw, color='#1ABC9C', lw=2, linestyle=':', label='Gaussian Kernel Reg')
    axes[0, 2].axhline(y=p0_opt, color='#2980B9', linestyle='--', alpha=0.7, label=f'p0 Limit ({p0_opt:.3f})')
    axes[0, 2].axhline(y=p1_opt, color='#27AE60', linestyle='--', alpha=0.7, label=f'p1 Limit ({p1_opt:.3f})')
    axes[0, 2].scatter([0], [hat_p0], color='#2980B9', marker='x', s=80, lw=2, label=f'Empirical p0 ({hat_p0:.3f})', zorder=5)
    axes[0, 2].scatter([1], [hat_p1], color='#27AE60', marker='x', s=80, lw=2, label=f'Empirical p1 ({hat_p1:.3f})', zorder=5)
    axes[0, 2].axhline(y=np.mean(y), color='gray', linestyle=':', alpha=0.5, label=f'Baseline Rate ({np.mean(y):.3f})')
    axes[0, 2].set_title(f"Custom MLE Model Fit\np0={p0_opt:.3f}, p1={p1_opt:.3f}, theta={theta_opt:.2f}\n[Empirical: p0={hat_p0:.2f}, p1={hat_p1:.2f}]", fontsize=11, fontweight='bold')
    axes[0, 2].set_xlabel("Raw Fraction")
    axes[0, 2].set_ylabel("Response probability")
    axes[0, 2].set_xlim(axes[0, 1].get_xlim())
    axes[0, 2].set_ylim(-0.05, 1.05)
    axes[0, 2].legend(loc='best', fontsize=8)
    axes[0, 2].grid(True, linestyle='--', alpha=0.3)
    
    # Plot 5: Logistic Regression trained on raw, plotted on Ranks
    X_full = sm.add_constant(vals)
    model_raw = sm.Logit(y, X_full).fit(disp=False)
    probs_pred = model_raw.predict(X_full)
    sort_idx = np.argsort(ranks)
    ranks_sorted = ranks[sort_idx]
    probs_pred_sorted = probs_pred[sort_idx]
    axes[0, 3].scatter(ranks, y, alpha=0.4, color='#8E44AD', label='Data')
    axes[0, 3].plot(ranks_sorted, probs_pred_sorted, color='#E67E22', lw=2.5, label='Logistic (Trained on Raw)')
    # Overlay projected Gaussian Kernel Regression (Trained on Raw)
    y_kernel_pred = kernel_regression(vals, y, vals)
    y_kernel_pred_sorted = y_kernel_pred[sort_idx]
    axes[0, 3].plot(ranks_sorted, y_kernel_pred_sorted, color='#1ABC9C', lw=2, linestyle=':', label='Kernel Reg (Trained on Raw)')
    axes[0, 3].axhline(y=np.mean(y), color='gray', linestyle='--', alpha=0.7, label=f'Baseline Rate ({np.mean(y):.3f})')
    axes[0, 3].set_title(f"Logistic Reg (Fitted on Raw)\nPlotted on Ranks\n[LRT p={lrt_p:.6f}]", fontsize=11, fontweight='bold')
    axes[0, 3].set_xlabel("Rank")
    axes[0, 3].set_ylabel("Response probability")
    axes[0, 3].set_ylim(-0.05, 1.05)
    axes[0, 3].legend(loc='best', fontsize=8)
    axes[0, 3].grid(True, linestyle='--', alpha=0.3)
    
    # Plot 6: Logistic regression on ranks
    _, lrt_p_ranks = compute_stats(ranks, y)
    sns.regplot(x=ranks, y=y, logistic=True, ax=axes[1, 1], scatter_kws={'alpha': 0.4}, line_kws={'color': '#E67E22'})
    axes[1, 1].axhline(y=np.mean(y), color='gray', linestyle='--', alpha=0.7, label=f'Baseline Rate ({np.mean(y):.3f})')
    # Fit & plot Gaussian Kernel Regression in rank space
    x_grid_ranks = np.linspace(1, len(ranks), 300)
    y_kernel_ranks = kernel_regression(ranks, y, x_grid_ranks)
    axes[1, 1].plot(x_grid_ranks, y_kernel_ranks, color='#1ABC9C', lw=2, linestyle=':', label='Gaussian Kernel Reg')
    axes[1, 1].set_title(f"Logistic Regression on Ranks\nSpearman R={chosen_row['Spearman_R']:.3f} (p={chosen_row['Spearman_P']:.4f})\n[LRT p={lrt_p_ranks:.6f}]", fontsize=11, fontweight='bold')
    axes[1, 1].set_xlabel("Rank")
    axes[1, 1].set_ylabel("Binary Response (0/1)")
    axes[1, 1].legend(loc='best', fontsize=8)
    axes[1, 1].grid(True, linestyle='--', alpha=0.3)
    
    # Plot 7: Custom model on ranks
    ranks_norm = (ranks - 1) / (len(ranks) - 1)
    p0_r, p1_r, theta_r, hat_p0_r, hat_p1_r = fit_custom_model(ranks_norm, y)
    axes[1, 2].scatter(ranks, y, alpha=0.4, color='#8E44AD', label='Data')
    x_grid_norm = (x_grid_ranks - 1) / (len(ranks) - 1)
    y_pred_r = p0_r + (p1_r - p0_r) * np.power(x_grid_norm, theta_r)
    axes[1, 2].plot(x_grid_ranks, y_pred_r, color='#E74C3C', lw=2.5, label='Custom Model')
    # Overlay Gaussian Kernel Regression
    axes[1, 2].plot(x_grid_ranks, y_kernel_ranks, color='#1ABC9C', lw=2, linestyle=':', label='Gaussian Kernel Reg')
    axes[1, 2].axhline(y=p0_r, color='#2980B9', linestyle='--', alpha=0.7, label=f'p0 Limit ({p0_r:.3f})')
    axes[1, 2].axhline(y=p1_r, color='#27AE60', linestyle='--', alpha=0.7, label=f'p1 Limit ({p1_r:.3f})')
    axes[1, 2].scatter([1], [hat_p0_r], color='#2980B9', marker='x', s=80, lw=2, label=f'Empirical p0 ({hat_p0_r:.3f})', zorder=5)
    axes[1, 2].scatter([len(ranks)], [hat_p1_r], color='#27AE60', marker='x', s=80, lw=2, label=f'Empirical p1 ({hat_p1_r:.3f})', zorder=5)
    axes[1, 2].axhline(y=np.mean(y), color='gray', linestyle=':', alpha=0.5, label=f'Baseline Rate ({np.mean(y):.3f})')
    axes[1, 2].set_title(f"Custom Model on Ranks\np0={p0_r:.3f}, p1={p1_r:.3f}, theta={theta_r:.2f}\n[Empirical: p0={hat_p0_r:.2f}, p1={hat_p1_r:.2f}]", fontsize=11, fontweight='bold')
    axes[1, 2].set_xlabel("Rank")
    axes[1, 2].set_ylabel("Response probability")
    axes[1, 2].set_xlim(axes[1, 1].get_xlim())
    axes[1, 2].set_ylim(-0.05, 1.05)
    axes[1, 2].legend(loc='best', fontsize=8)
    axes[1, 2].grid(True, linestyle='--', alpha=0.3)
    
    # Plot 8: Custom MLE Model (Trained on Raw) Plotted on Ranks
    probs_custom = p0_opt + (p1_opt - p0_opt) * np.power(vals, theta_opt)
    probs_custom_sorted = probs_custom[sort_idx]
    axes[1, 3].scatter(ranks, y, alpha=0.4, color='#8E44AD', label='Data')
    axes[1, 3].plot(ranks_sorted, probs_custom_sorted, color='#E74C3C', lw=2.5, label='Custom (Trained on Raw)')
    # Overlay projected Gaussian Kernel Regression (Trained on Raw)
    axes[1, 3].plot(ranks_sorted, y_kernel_pred_sorted, color='#1ABC9C', lw=2, linestyle=':', label='Kernel Reg (Trained on Raw)')
    axes[1, 3].axhline(y=p0_opt, color='#2980B9', linestyle='--', alpha=0.7, label=f'p0 Limit ({p0_opt:.3f})')
    axes[1, 3].axhline(y=p1_opt, color='#27AE60', linestyle='--', alpha=0.7, label=f'p1 Limit ({p1_opt:.3f})')
    axes[1, 3].scatter([1], [hat_p0], color='#2980B9', marker='x', s=80, lw=2, label=f'Empirical p0 ({hat_p0:.3f})', zorder=5)
    axes[1, 3].scatter([len(ranks)], [hat_p1], color='#27AE60', marker='x', s=80, lw=2, label=f'Empirical p1 ({hat_p1:.3f})', zorder=5)
    axes[1, 3].axhline(y=np.mean(y), color='gray', linestyle=':', alpha=0.5, label=f'Baseline Rate ({np.mean(y):.3f})')
    axes[1, 3].set_title(f"Custom MLE (Fitted on Raw)\nPlotted on Ranks\np0={p0_opt:.3f}, p1={p1_opt:.3f}, theta={theta_opt:.2f}", fontsize=11, fontweight='bold')
    axes[1, 3].set_xlabel("Rank")
    axes[1, 3].set_ylabel("Response probability")
    axes[1, 3].set_xlim(axes[1, 1].get_xlim())
    axes[1, 3].set_ylim(-0.05, 1.05)
    axes[1, 3].legend(loc='best', fontsize=8)
    axes[1, 3].grid(True, linestyle='--', alpha=0.3)
    
    # Sync Plot 5's x-limits
    axes[0, 3].set_xlim(axes[1, 1].get_xlim())
    
    plt.suptitle(f"Spearman vs Pearson Discrepancy Analysis in Combined-Melanoma: {top_cell_type}", fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    plot_path = output_dir / "spearman_vs_pearson_illustration.svg"
    plt.savefig(plot_path, bbox_inches='tight')
    plt.close()
    print(f"\nSaved illustration plot to {plot_path}")

    # Generate scatter plot comparing all features' coefficients
    plot_path_all = output_dir / "spearman_vs_pearson_all_features.svg"
    plt.figure(figsize=(10, 8), dpi=300)
    
    sig_color = []
    for idx, row in df_comp.iterrows():
        sp_sig = row['Spearman_P'] < 0.05
        pe_sig = row['Pearson_P'] < 0.05
        if sp_sig and pe_sig:
            sig_color.append('Both Significant (p < 0.05)')
        elif sp_sig:
            sig_color.append('Spearman Only Significant (p < 0.05)')
        elif pe_sig:
            sig_color.append('Pearson Only Significant (p < 0.05)')
        else:
            sig_color.append('Neither Significant')
    df_comp['Significance'] = sig_color
    
    sns.scatterplot(
        data=df_comp,
        x='Pearson_R',
        y='Spearman_R',
        hue='Significance',
        size='Skewness',
        sizes=(40, 250),
        palette={
            'Both Significant (p < 0.05)': '#2ECC71',
            'Spearman Only Significant (p < 0.05)': '#3498DB',
            'Pearson Only Significant (p < 0.05)': '#E74C3C',
            'Neither Significant': '#BDC3C7'
        },
        alpha=0.8,
        edgecolor='black',
        linewidth=0.8
    )
    
    # Draw y = x diagonal line
    lims = [
        min(plt.xlim()[0], plt.ylim()[0]),
        max(plt.xlim()[1], plt.ylim()[1])
    ]
    plt.plot(lims, lims, color='#7F8C8D', linestyle='--', alpha=0.7, label='y = x (Perfect Agreement)')
    
    # Label top discrepant features
    for idx, row in df_discrepant.head(5).iterrows():
        plt.annotate(
            row['CellType'],
            xy=(row['Pearson_R'], row['Spearman_R']),
            xytext=(row['Pearson_R'] + 0.02, row['Spearman_R'] - 0.02),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=0.6),
            fontsize=9,
            fontweight='bold'
        )
        
    plt.title("Spearman vs. Pearson Correlation Coefficients (Combined-Melanoma)\nAcross All Inferred Cell Types / Subclusters", fontsize=12, fontweight='bold', pad=15)
    plt.xlabel("Pearson Correlation Coefficient (R vs NR)", fontsize=10, fontweight='bold')
    plt.ylabel("Spearman Correlation Coefficient (R vs NR)", fontsize=10, fontweight='bold')
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.legend(loc='lower right', frameon=True)
    sns.despine()
    plt.tight_layout()
    plt.savefig(plot_path_all, bbox_inches='tight')
    plt.close()
    print(f"Saved all features comparison scatter plot to {plot_path_all}")

    # Generate concordant case study plot (where both correlations are highly significant)
    concordant_cell_type = "18_sub_2"
    concordant_row = df_comp[df_comp['CellType'] == concordant_cell_type].iloc[0]
    
    vals_conc = df_fracs_comb[concordant_cell_type].values
    
    fig, axes = plt.subplots(2, 4, figsize=(24, 10), dpi=300)
    
    mwu_p_conc, lrt_p_conc = compute_stats(vals_conc, y)
    
    # Plot 1: Box plot of raw fractions
    sns.boxplot(x=y_labels, y=vals_conc, ax=axes[0, 0], palette={'Non-Responder': '#EC7063', 'Responder': '#5DADE2'}, hue=y_labels, legend=False)
    sns.stripplot(x=y_labels, y=vals_conc, color='black', alpha=0.3, jitter=0.2, ax=axes[0, 0])
    axes[0, 0].set_title(f"Raw Fractions of {concordant_cell_type}\n(CD8+ T cells - Highly Significant in Both)\n[MWU p={mwu_p_conc:.7f}]", fontsize=11, fontweight='bold')
    axes[0, 0].set_xlabel("Response status")
    axes[0, 0].set_ylabel("Inferred Fraction")
    ymax_c = np.percentile(vals_conc, 95)
    if ymax_c > 0:
        axes[0, 0].set_ylim(0, ymax_c)
        
    # Plot 2: Box plot of ranks
    ranks_c = stats.rankdata(vals_conc)
    sns.boxplot(x=y_labels, y=ranks_c, ax=axes[1, 0], palette={'Non-Responder': '#EC7063', 'Responder': '#5DADE2'}, hue=y_labels, legend=False)
    sns.stripplot(x=y_labels, y=ranks_c, color='black', alpha=0.3, jitter=0.2, ax=axes[1, 0])
    axes[1, 0].set_title(f"Rank-Transformed Fractions of {concordant_cell_type}\n(Spearman Correlation Basis)\n[MWU p={mwu_p_conc:.7f}]", fontsize=11, fontweight='bold')
    axes[1, 0].set_xlabel("Response status")
    axes[1, 0].set_ylabel("Rank")
    
    # Plot 3: Regplot of raw fractions vs binary response
    sns.regplot(x=vals_conc, y=y, logistic=True, ax=axes[0, 1], scatter_kws={'alpha': 0.4}, line_kws={'color': '#2ECC71'})
    axes[0, 1].axhline(y=np.mean(y), color='gray', linestyle='--', alpha=0.7, label=f'Baseline Rate ({np.mean(y):.3f})')
    # Add vertical percentile lines
    pct_25_c = np.percentile(vals_conc, 25)
    pct_50_c = np.percentile(vals_conc, 50)
    pct_75_c = np.percentile(vals_conc, 75)
    pct_95_c = np.percentile(vals_conc, 95)
    axes[0, 1].axvline(x=pct_25_c, color='#8E44AD', linestyle=':', alpha=0.6, label=f'25th Pct ({pct_25_c:.4f})')
    axes[0, 1].axvline(x=pct_50_c, color='#8E44AD', linestyle='--', alpha=0.6, label=f'Median ({pct_50_c:.4f})')
    axes[0, 1].axvline(x=pct_75_c, color='#8E44AD', linestyle='-.', alpha=0.6, label=f'75th Pct ({pct_75_c:.4f})')
    axes[0, 1].axvline(x=pct_95_c, color='#C0392B', linestyle='--', alpha=0.6, label=f'95th Pct ({pct_95_c:.4f})')
    # Fit & plot Gaussian Kernel Regression in raw space
    x_grid_raw = np.linspace(np.min(vals_conc), np.max(vals_conc), 300)
    y_kernel_raw = kernel_regression(vals_conc, y, x_grid_raw)
    axes[0, 1].plot(x_grid_raw, y_kernel_raw, color='#1ABC9C', lw=2, linestyle=':', label='Gaussian Kernel Reg')
    axes[0, 1].set_title(f"Linear/Logistic Regression Agreement\nPearson R={concordant_row['Pearson_R']:.3f} (p={concordant_row['Pearson_P']:.7f})\nSpearman R={concordant_row['Spearman_R']:.3f} (p={concordant_row['Spearman_P']:.7f})\n[LRT p={lrt_p_conc:.7f}]", fontsize=11, fontweight='bold')
    axes[0, 1].set_xlabel("Raw Fraction")
    axes[0, 1].set_ylabel("Binary Response (0/1)")
    axes[0, 1].legend(loc='best', fontsize=8)
    
    # Plot 4: Custom model fit
    p0_opt, p1_opt, theta_opt, hat_p0, hat_p1 = fit_custom_model(vals_conc, y)
    axes[0, 2].scatter(vals_conc, y, alpha=0.4, color='#8E44AD', label='Data')
    x_grid = np.linspace(0, 1.0, 300)
    y_pred = p0_opt + (p1_opt - p0_opt) * np.power(x_grid, theta_opt)
    axes[0, 2].plot(x_grid, y_pred, color='#E74C3C', lw=2.5, label='Custom Model')
    # Overlay Gaussian Kernel Regression
    axes[0, 2].plot(x_grid_raw, y_kernel_raw, color='#1ABC9C', lw=2, linestyle=':', label='Gaussian Kernel Reg')
    axes[0, 2].axhline(y=p0_opt, color='#2980B9', linestyle='--', alpha=0.7, label=f'p0 Limit ({p0_opt:.3f})')
    axes[0, 2].axhline(y=p1_opt, color='#27AE60', linestyle='--', alpha=0.7, label=f'p1 Limit ({p1_opt:.3f})')
    axes[0, 2].scatter([0], [hat_p0], color='#2980B9', marker='x', s=80, lw=2, label=f'Empirical p0 ({hat_p0:.3f})', zorder=5)
    axes[0, 2].scatter([1], [hat_p1], color='#27AE60', marker='x', s=80, lw=2, label=f'Empirical p1 ({hat_p1:.3f})', zorder=5)
    axes[0, 2].axhline(y=np.mean(y), color='gray', linestyle=':', alpha=0.5, label=f'Baseline Rate ({np.mean(y):.3f})')
    axes[0, 2].set_title(f"Custom MLE Model Fit\np0={p0_opt:.3f}, p1={p1_opt:.3f}, theta={theta_opt:.2f}\n[Empirical: p0={hat_p0:.2f}, p1={hat_p1:.2f}]", fontsize=11, fontweight='bold')
    axes[0, 2].set_xlabel("Raw Fraction")
    axes[0, 2].set_ylabel("Response probability")
    axes[0, 2].set_xlim(axes[0, 1].get_xlim())
    axes[0, 2].set_ylim(-0.05, 1.05)
    axes[0, 2].legend(loc='best', fontsize=8)
    axes[0, 2].grid(True, linestyle='--', alpha=0.3)
    
    # Plot 5: Logistic Regression trained on raw, plotted on Ranks
    X_full = sm.add_constant(vals_conc)
    model_raw = sm.Logit(y, X_full).fit(disp=False)
    probs_pred = model_raw.predict(X_full)
    sort_idx = np.argsort(ranks_c)
    ranks_sorted = ranks_c[sort_idx]
    probs_pred_sorted = probs_pred[sort_idx]
    axes[0, 3].scatter(ranks_c, y, alpha=0.4, color='#8E44AD', label='Data')
    axes[0, 3].plot(ranks_sorted, probs_pred_sorted, color='#2ECC71', lw=2.5, label='Logistic (Trained on Raw)')
    # Overlay projected Gaussian Kernel Regression (Trained on Raw)
    y_kernel_pred = kernel_regression(vals_conc, y, vals_conc)
    y_kernel_pred_sorted = y_kernel_pred[sort_idx]
    axes[0, 3].plot(ranks_sorted, y_kernel_pred_sorted, color='#1ABC9C', lw=2, linestyle=':', label='Kernel Reg (Trained on Raw)')
    axes[0, 3].axhline(y=np.mean(y), color='gray', linestyle='--', alpha=0.7, label=f'Baseline Rate ({np.mean(y):.3f})')
    axes[0, 3].set_title(f"Logistic Reg (Fitted on Raw)\nPlotted on Ranks\n[LRT p={lrt_p_conc:.7f}]", fontsize=11, fontweight='bold')
    axes[0, 3].set_xlabel("Rank")
    axes[0, 3].set_ylabel("Response probability")
    axes[0, 3].set_ylim(-0.05, 1.05)
    axes[0, 3].legend(loc='best', fontsize=8)
    axes[0, 3].grid(True, linestyle='--', alpha=0.3)
    
    # Plot 6: Logistic regression on ranks
    sns.regplot(x=ranks_c, y=y, logistic=True, ax=axes[1, 1], scatter_kws={'alpha': 0.4}, line_kws={'color': '#2ECC71'})
    axes[1, 1].axhline(y=np.mean(y), color='gray', linestyle='--', alpha=0.7, label=f'Baseline Rate ({np.mean(y):.3f})')
    # Fit & plot Gaussian Kernel Regression in rank space
    x_grid_ranks = np.linspace(1, len(ranks_c), 300)
    y_kernel_ranks = kernel_regression(ranks_c, y, x_grid_ranks)
    axes[1, 1].plot(x_grid_ranks, y_kernel_ranks, color='#1ABC9C', lw=2, linestyle=':', label='Gaussian Kernel Reg')
    axes[1, 1].set_title(f"Logistic Regression on Ranks\nSpearman R={concordant_row['Spearman_R']:.3f} (p={concordant_row['Spearman_P']:.7f})\n[LRT p={lrt_p_ranks:.7f}]", fontsize=11, fontweight='bold')
    axes[1, 1].set_xlabel("Rank")
    axes[1, 1].set_ylabel("Binary Response (0/1)")
    axes[1, 1].legend(loc='best', fontsize=8)
    axes[1, 1].grid(True, linestyle='--', alpha=0.3)
    
    # Plot 7: Custom model on ranks
    ranks_norm = (ranks_c - 1) / (len(ranks_c) - 1)
    p0_r, p1_r, theta_r, hat_p0_r, hat_p1_r = fit_custom_model(ranks_norm, y)
    axes[1, 2].scatter(ranks_c, y, alpha=0.4, color='#8E44AD', label='Data')
    x_grid_norm = (x_grid_ranks - 1) / (len(ranks_c) - 1)
    y_pred_r = p0_r + (p1_r - p0_r) * np.power(x_grid_norm, theta_r)
    axes[1, 2].plot(x_grid_ranks, y_pred_r, color='#E74C3C', lw=2.5, label='Custom Model')
    # Overlay Gaussian Kernel Regression
    axes[1, 2].plot(x_grid_ranks, y_kernel_ranks, color='#1ABC9C', lw=2, linestyle=':', label='Gaussian Kernel Reg')
    axes[1, 2].axhline(y=p0_r, color='#2980B9', linestyle='--', alpha=0.7, label=f'p0 Limit ({p0_r:.3f})')
    axes[1, 2].axhline(y=p1_r, color='#27AE60', linestyle='--', alpha=0.7, label=f'p1 Limit ({p1_r:.3f})')
    axes[1, 2].scatter([1], [hat_p0_r], color='#2980B9', marker='x', s=80, lw=2, label=f'Empirical p0 ({hat_p0_r:.3f})', zorder=5)
    axes[1, 2].scatter([len(ranks_c)], [hat_p1_r], color='#27AE60', marker='x', s=80, lw=2, label=f'Empirical p1 ({hat_p1_r:.3f})', zorder=5)
    axes[1, 2].axhline(y=np.mean(y), color='gray', linestyle=':', alpha=0.5, label=f'Baseline Rate ({np.mean(y):.3f})')
    axes[1, 2].set_title(f"Custom Model on Ranks\np0={p0_r:.3f}, p1={p1_r:.3f}, theta={theta_r:.2f}\n[Empirical: p0={hat_p0_r:.2f}, p1={hat_p1_r:.2f}]", fontsize=11, fontweight='bold')
    axes[1, 2].set_xlabel("Rank")
    axes[1, 2].set_ylabel("Response probability")
    axes[1, 2].set_xlim(axes[1, 1].get_xlim())
    axes[1, 2].set_ylim(-0.05, 1.05)
    axes[1, 2].legend(loc='best', fontsize=8)
    axes[1, 2].grid(True, linestyle='--', alpha=0.3)
    
    # Plot 8: Custom MLE Model (Trained on Raw) Plotted on Ranks
    probs_custom = p0_opt + (p1_opt - p0_opt) * np.power(vals_conc, theta_opt)
    probs_custom_sorted = probs_custom[sort_idx]
    axes[1, 3].scatter(ranks_c, y, alpha=0.4, color='#8E44AD', label='Data')
    axes[1, 3].plot(ranks_sorted, probs_custom_sorted, color='#E74C3C', lw=2.5, label='Custom (Trained on Raw)')
    # Overlay projected Gaussian Kernel Regression (Trained on Raw)
    axes[1, 3].plot(ranks_sorted, y_kernel_pred_sorted, color='#1ABC9C', lw=2, linestyle=':', label='Kernel Reg (Trained on Raw)')
    axes[1, 3].axhline(y=p0_opt, color='#2980B9', linestyle='--', alpha=0.7, label=f'p0 Limit ({p0_opt:.3f})')
    axes[1, 3].axhline(y=p1_opt, color='#27AE60', linestyle='--', alpha=0.7, label=f'p1 Limit ({p1_opt:.3f})')
    axes[1, 3].scatter([1], [hat_p0], color='#2980B9', marker='x', s=80, lw=2, label=f'Empirical p0 ({hat_p0:.3f})', zorder=5)
    axes[1, 3].scatter([len(ranks_c)], [hat_p1], color='#27AE60', marker='x', s=80, lw=2, label=f'Empirical p1 ({hat_p1:.3f})', zorder=5)
    axes[1, 3].axhline(y=np.mean(y), color='gray', linestyle=':', alpha=0.5, label=f'Baseline Rate ({np.mean(y):.3f})')
    axes[1, 3].set_title(f"Custom MLE (Fitted on Raw)\nPlotted on Ranks\np0={p0_opt:.3f}, p1={p1_opt:.3f}, theta={theta_opt:.2f}", fontsize=11, fontweight='bold')
    axes[1, 3].set_xlabel("Rank")
    axes[1, 3].set_ylabel("Response probability")
    axes[1, 3].set_xlim(axes[1, 1].get_xlim())
    axes[1, 3].set_ylim(-0.05, 1.05)
    axes[1, 3].legend(loc='best', fontsize=8)
    axes[1, 3].grid(True, linestyle='--', alpha=0.3)
    
    # Sync Plot 5's x-limits
    axes[0, 3].set_xlim(axes[1, 1].get_xlim())
    
    plt.suptitle(f"Spearman & Pearson Agreement in Combined-Melanoma: {concordant_cell_type}", fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    plot_path_conc = output_dir / "spearman_vs_pearson_concordant.svg"
    plt.savefig(plot_path_conc, bbox_inches='tight')
    plt.close()
    print(f"Saved concordant comparison plot to {plot_path_conc}")

    # Generate the 3-panel plot for all cell types
    generate_plots_for_all_cell_types(df_fracs_comb, y_series_comb, df_comp, output_dir)


def generate_plots_for_all_cell_types(df_fracs_comb, y_series_comb, df_comp, output_dir):
    plot_dir = output_dir / "spearman_vs_pearson_by_celltype"
    plot_dir.mkdir(parents=True, exist_ok=True)
    
    y = y_series_comb.values
    y_labels = pd.Series(y).map({0: 'Non-Responder', 1: 'Responder'})
    
    print(f"\nGenerating Spearman vs Pearson plots for all {len(df_comp)} cell types...")
    for idx, row in df_comp.iterrows():
        cell_type = row['CellType']
        vals = df_fracs_comb[cell_type].values
        
        fig, axes = plt.subplots(2, 4, figsize=(24, 10), dpi=150)
        
        mwu_p_all, lrt_p_all = compute_stats(vals, y)
        
        # Plot 1: Box plot of raw fractions
        sns.boxplot(x=y_labels, y=vals, ax=axes[0, 0], palette={'Non-Responder': '#EC7063', 'Responder': '#5DADE2'}, hue=y_labels, legend=False)
        sns.stripplot(x=y_labels, y=vals, color='black', alpha=0.3, jitter=0.2, ax=axes[0, 0])
        axes[0, 0].set_title(f"Raw Fractions of {cell_type}\n(Highly Skewed: Skew={row['Skewness']:.2f})\n[MWU p={mwu_p_all:.6f}]", fontsize=11, fontweight='bold')
        axes[0, 0].set_xlabel("Response status")
        axes[0, 0].set_ylabel("Inferred Fraction")
        ymax = np.percentile(vals, 95)
        if ymax > 0:
            axes[0, 0].set_ylim(0, ymax)
            
        # Plot 2: Box plot of ranks
        ranks = stats.rankdata(vals)
        sns.boxplot(x=y_labels, y=ranks, ax=axes[1, 0], palette={'Non-Responder': '#EC7063', 'Responder': '#5DADE2'}, hue=y_labels, legend=False)
        sns.stripplot(x=y_labels, y=ranks, color='black', alpha=0.3, jitter=0.2, ax=axes[1, 0])
        axes[1, 0].set_title(f"Rank-Transformed Fractions of {cell_type}\n(Spearman Correlation Basis)\n[MWU p={mwu_p_all:.6f}]", fontsize=11, fontweight='bold')
        axes[1, 0].set_xlabel("Response status")
        axes[1, 0].set_ylabel("Rank")
        
        # Plot 3: Regplot of raw fractions vs binary response
        sns.regplot(x=vals, y=y, logistic=True, ax=axes[0, 1], scatter_kws={'alpha': 0.4}, line_kws={'color': '#E67E22'})
        axes[0, 1].axhline(y=np.mean(y), color='gray', linestyle='--', alpha=0.7, label=f'Baseline Rate ({np.mean(y):.3f})')
        # Add vertical percentile lines
        pct_25_all = np.percentile(vals, 25)
        pct_50_all = np.percentile(vals, 50)
        pct_75_all = np.percentile(vals, 75)
        pct_95_all = np.percentile(vals, 95)
        axes[0, 1].axvline(x=pct_25_all, color='#8E44AD', linestyle=':', alpha=0.6, label=f'25th Pct ({pct_25_all:.4f})')
        axes[0, 1].axvline(x=pct_50_all, color='#8E44AD', linestyle='--', alpha=0.6, label=f'Median ({pct_50_all:.4f})')
        axes[0, 1].axvline(x=pct_75_all, color='#8E44AD', linestyle='-.', alpha=0.6, label=f'75th Pct ({pct_75_all:.4f})')
        axes[0, 1].axvline(x=pct_95_all, color='#C0392B', linestyle='--', alpha=0.6, label=f'95th Pct ({pct_95_all:.4f})')
        axes[0, 1].set_title(f"Linear/Logistic Regression Leverage\nPearson R={row['Pearson_R']:.3f} (p={row['Pearson_P']:.4f})\nSpearman R={row['Spearman_R']:.3f} (p={row['Spearman_P']:.4f})\n[LRT p={lrt_p_all:.6f}]", fontsize=11, fontweight='bold')
        axes[0, 1].set_xlabel("Raw Fraction")
        axes[0, 1].set_ylabel("Binary Response (0/1)")
        axes[0, 1].legend(loc='best', fontsize=8)
        
        # Plot 4: Custom model fit
        p0_opt, p1_opt, theta_opt, hat_p0, hat_p1 = fit_custom_model(vals, y)
        axes[0, 2].scatter(vals, y, alpha=0.4, color='#8E44AD', label='Data')
        x_grid = np.linspace(0, 1.0, 300)
        y_pred = p0_opt + (p1_opt - p0_opt) * np.power(x_grid, theta_opt)
        axes[0, 2].plot(x_grid, y_pred, color='#E74C3C', lw=2.5, label='Custom Model')
        axes[0, 2].axhline(y=p0_opt, color='#2980B9', linestyle='--', alpha=0.7, label=f'p0 Limit ({p0_opt:.3f})')
        axes[0, 2].axhline(y=p1_opt, color='#27AE60', linestyle='--', alpha=0.7, label=f'p1 Limit ({p1_opt:.3f})')
        axes[0, 2].scatter([0], [hat_p0], color='#2980B9', marker='x', s=80, lw=2, label=f'Empirical p0 ({hat_p0:.3f})', zorder=5)
        axes[0, 2].scatter([1], [hat_p1], color='#27AE60', marker='x', s=80, lw=2, label=f'Empirical p1 ({hat_p1:.3f})', zorder=5)
        axes[0, 2].axhline(y=np.mean(y), color='gray', linestyle=':', alpha=0.5, label=f'Baseline Rate ({np.mean(y):.3f})')
        axes[0, 2].set_title(f"Custom MLE Model Fit\np0={p0_opt:.3f}, p1={p1_opt:.3f}, theta={theta_opt:.2f}\n[Empirical: p0={hat_p0:.2f}, p1={hat_p1:.2f}]", fontsize=11, fontweight='bold')
        axes[0, 2].set_xlabel("Raw Fraction")
        axes[0, 2].set_ylabel("Response probability")
        axes[0, 2].set_xlim(axes[0, 1].get_xlim())
        axes[0, 2].set_ylim(-0.05, 1.05)
        axes[0, 2].legend(loc='best', fontsize=8)
        axes[0, 2].grid(True, linestyle='--', alpha=0.3)
        
        # Plot 5: Logistic Regression trained on raw, plotted on Ranks
        X_full = sm.add_constant(vals)
        model_raw = sm.Logit(y, X_full).fit(disp=False)
        probs_pred = model_raw.predict(X_full)
        sort_idx = np.argsort(ranks)
        ranks_sorted = ranks[sort_idx]
        probs_pred_sorted = probs_pred[sort_idx]
        axes[0, 3].scatter(ranks, y, alpha=0.4, color='#8E44AD', label='Data')
        axes[0, 3].plot(ranks_sorted, probs_pred_sorted, color='#E67E22', lw=2.5, label='Logistic (Trained on Raw)')
        axes[0, 3].axhline(y=np.mean(y), color='gray', linestyle='--', alpha=0.7, label=f'Baseline Rate ({np.mean(y):.3f})')
        axes[0, 3].set_title(f"Logistic Reg (Fitted on Raw)\nPlotted on Ranks\n[LRT p={lrt_p_all:.6f}]", fontsize=11, fontweight='bold')
        axes[0, 3].set_xlabel("Rank")
        axes[0, 3].set_ylabel("Response probability")
        axes[0, 3].set_ylim(-0.05, 1.05)
        axes[0, 3].legend(loc='best', fontsize=8)
        axes[0, 3].grid(True, linestyle='--', alpha=0.3)
        
        # Plot 6: Logistic regression on ranks
        _, lrt_p_ranks = compute_stats(ranks, y)
        sns.regplot(x=ranks, y=y, logistic=True, ax=axes[1, 1], scatter_kws={'alpha': 0.4}, line_kws={'color': '#E67E22'})
        axes[1, 1].axhline(y=np.mean(y), color='gray', linestyle='--', alpha=0.7, label=f'Baseline Rate ({np.mean(y):.3f})')
        axes[1, 1].set_title(f"Logistic Regression on Ranks\nSpearman R={row['Spearman_R']:.3f} (p={row['Spearman_P']:.4f})\n[LRT p={lrt_p_ranks:.6f}]", fontsize=11, fontweight='bold')
        axes[1, 1].set_xlabel("Rank")
        axes[1, 1].set_ylabel("Binary Response (0/1)")
        axes[1, 1].legend(loc='best', fontsize=8)
        axes[1, 1].grid(True, linestyle='--', alpha=0.3)
        
        # Plot 7: Custom model on ranks
        ranks_norm = (ranks - 1) / (len(ranks) - 1)
        p0_r, p1_r, theta_r, hat_p0_r, hat_p1_r = fit_custom_model(ranks_norm, y)
        axes[1, 2].scatter(ranks, y, alpha=0.4, color='#8E44AD', label='Data')
        x_grid_ranks = np.linspace(1, len(ranks), 300)
        x_grid_norm = (x_grid_ranks - 1) / (len(ranks) - 1)
        y_pred_r = p0_r + (p1_r - p0_r) * np.power(x_grid_norm, theta_r)
        axes[1, 2].plot(x_grid_ranks, y_pred_r, color='#E74C3C', lw=2.5, label='Custom Model')
        axes[1, 2].axhline(y=p0_r, color='#2980B9', linestyle='--', alpha=0.7, label=f'p0 Limit ({p0_r:.3f})')
        axes[1, 2].axhline(y=p1_r, color='#27AE60', linestyle='--', alpha=0.7, label=f'p1 Limit ({p1_r:.3f})')
        axes[1, 2].scatter([1], [hat_p0_r], color='#2980B9', marker='x', s=80, lw=2, label=f'Empirical p0 ({hat_p0_r:.3f})', zorder=5)
        axes[1, 2].scatter([len(ranks)], [hat_p1_r], color='#27AE60', marker='x', s=80, lw=2, label=f'Empirical p1 ({hat_p1_r:.3f})', zorder=5)
        axes[1, 2].axhline(y=np.mean(y), color='gray', linestyle=':', alpha=0.5, label=f'Baseline Rate ({np.mean(y):.3f})')
        axes[1, 2].set_title(f"Custom Model on Ranks\np0={p0_r:.3f}, p1={p1_r:.3f}, theta={theta_r:.2f}\n[Empirical: p0={hat_p0_r:.2f}, p1={hat_p1_r:.2f}]", fontsize=11, fontweight='bold')
        axes[1, 2].set_xlabel("Rank")
        axes[1, 2].set_ylabel("Response probability")
        axes[1, 2].set_xlim(axes[1, 1].get_xlim())
        axes[1, 2].set_ylim(-0.05, 1.05)
        axes[1, 2].legend(loc='best', fontsize=8)
        axes[1, 2].grid(True, linestyle='--', alpha=0.3)
        
        # Plot 8: Custom MLE Model (Trained on Raw) Plotted on Ranks
        probs_custom = p0_opt + (p1_opt - p0_opt) * np.power(vals, theta_opt)
        probs_custom_sorted = probs_custom[sort_idx]
        axes[1, 3].scatter(ranks, y, alpha=0.4, color='#8E44AD', label='Data')
        axes[1, 3].plot(ranks_sorted, probs_custom_sorted, color='#E74C3C', lw=2.5, label='Custom (Trained on Raw)')
        axes[1, 3].axhline(y=p0_opt, color='#2980B9', linestyle='--', alpha=0.7, label=f'p0 Limit ({p0_opt:.3f})')
        axes[1, 3].axhline(y=p1_opt, color='#27AE60', linestyle='--', alpha=0.7, label=f'p1 Limit ({p1_opt:.3f})')
        axes[1, 3].scatter([1], [hat_p0], color='#2980B9', marker='x', s=80, lw=2, label=f'Empirical p0 ({hat_p0:.3f})', zorder=5)
        axes[1, 3].scatter([len(ranks)], [hat_p1], color='#27AE60', marker='x', s=80, lw=2, label=f'Empirical p1 ({hat_p1:.3f})', zorder=5)
        axes[1, 3].axhline(y=np.mean(y), color='gray', linestyle=':', alpha=0.5, label=f'Baseline Rate ({np.mean(y):.3f})')
        axes[1, 3].set_title(f"Custom MLE (Fitted on Raw)\nPlotted on Ranks\np0={p0_opt:.3f}, p1={p1_opt:.3f}, theta={theta_opt:.2f}", fontsize=11, fontweight='bold')
        axes[1, 3].set_xlabel("Rank")
        axes[1, 3].set_ylabel("Response probability")
        axes[1, 3].set_xlim(axes[1, 1].get_xlim())
        axes[1, 3].set_ylim(-0.05, 1.05)
        axes[1, 3].legend(loc='best', fontsize=8)
        axes[1, 3].grid(True, linestyle='--', alpha=0.3)
        
        # Sync Plot 5's x-limits
        axes[0, 3].set_xlim(axes[1, 1].get_xlim())
        
        plt.suptitle(f"Spearman vs Pearson Analysis in Combined-Melanoma: {cell_type}", fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        plot_path = plot_dir / f"spearman_vs_pearson_{cell_type}.svg"
        plt.savefig(plot_path, bbox_inches='tight')
        plt.close()
        
    print(f"Completed generating all cell type plots in {plot_dir}")


if __name__ == "__main__":
    main()
