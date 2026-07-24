#!/usr/bin/env python3
"""
Logistic Regression, Likelihood Ratio Test (LRT), Wald Test, Pearson & Spearman Correlations
on LM22_no_IG Cell Type Fractions for Immunotherapy Response Prediction across iAtlas & Combined Cohorts.

Calculates and compares 4 statistical tests for every cell type:
1. Likelihood Ratio Test (LRT): full logit model vs null model
2. Wald Test (z-test): H0: Odds Ratio (OR) = 1 (beta_1 = 0)
3. Pearson Correlation: linear association (r, p-value)
4. Spearman Rank Correlation: monotonic rank association (rho, p-value)

Formats all metrics and p-values cleanly limiting decimal precision.
"""

import sys
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2, norm, pearsonr, spearmanr
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
        
    # Pearson and Spearman correlations
    try:
        r_p, p_p = pearsonr(X, Y)
    except Exception:
        r_p, p_p = np.nan, np.nan
        
    try:
        rho_s, p_s = spearmanr(X, Y)
    except Exception:
        rho_s, p_s = np.nan, np.nan
        
    # Logistic Regression, LRT & Wald Test
    try:
        X_null = np.ones((N, 1))
        null_mod = sm.Logit(Y, X_null)
        null_res = null_mod.fit(disp=False, maxiter=100)
        llf_null = null_res.llf
        
        X_full = sm.add_constant(X)
        full_mod = sm.Logit(Y, X_full)
        full_res = full_mod.fit(disp=False, maxiter=100)
        llf_full = full_res.llf
        
        # LRT Statistic & p-value
        lrt_stat = max(0.0, 2.0 * (llf_full - llf_null))
        p_val_lrt = float(chi2.sf(lrt_stat, df=1))
        
        beta = float(full_res.params[1])
        se = float(full_res.bse[1])
        
        # Wald Test (z-test for OR)
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
    """Formats Odds Ratio or CI value limiting decimal precision cleanly."""
    if np.isnan(val) or np.isinf(val):
        return "N/A"
    elif val >= 1000.0 or (val <= 0.001 and val > 0):
        return f"{val:.2e}"
    else:
        return f"{val:.2f}"


def format_p_val(p: float) -> str:
    """Formats p-value limiting decimal precision."""
    if np.isnan(p):
        return "N/A"
    elif p < 0.0001:
        return f"{p:.2e}"
    else:
        return f"{p:.4f}"


def compute_candle_limit(vals_r: np.ndarray, vals_nr: np.ndarray) -> float:
    """Computes upper Y-limit based on highest candle (upper fence: Q3 + 1.5 * IQR) between R and NR groups."""
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
    """Converts p-value into standard significance asterisks."""
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
    output_dir = project_root / "output" / "logistic-lrt-response"
    output_dir.mkdir(parents=True, exist_ok=True)
    boxplots_dir = output_dir / "cohort_boxplots"
    boxplots_dir.mkdir(parents=True, exist_ok=True)
    
    iatlas_cohorts = [
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
    
    iatlas_cancer_map = {
        'Combined-Melanoma': ['Hugo-iAtlas', 'Riaz-iAtlas', 'Liu-iAtlas', 'Gide-iAtlas'],
        'Combined-BLCA': ['Rosenberg-iAtlas'],
        'Combined-PAAD': ['Padron-iAtlas'],
        'Combined-BRCA': ['Anders-iAtlas'],
        'Combined-RCC': ['McDermott-iAtlas', 'Choueiri-iAtlas'],
        'Combined-All-iAtlas': iatlas_cohorts
    }
    
    # 1. Load Clinical Response Labels for all iAtlas cohorts
    print("\nLoading clinical response labels for iAtlas cohorts...")
    clinical_map = {}
    for cohort in iatlas_cohorts:
        df_clin = load_clinical_response(cohort, lair)
        if df_clin is not None:
            clinical_map[cohort] = df_clin
            print(f"  {cohort:<20}: N={len(df_clin)} samples (R={sum(df_clin['response']=='R')}, NR={sum(df_clin['response']=='NR')})")
            
    ref_variants = ["LM22_no_IG", "LM22_no_IG_iter10"]
    
    for ref_name in ref_variants:
        print(f"\n=======================================================")
        print(f"RUNNING MULTI-TEST STATISTICAL EVALUATION FOR: {ref_name}")
        print(f"=======================================================")
        
        fractions_map = {}
        for cohort in iatlas_cohorts:
            p = deconv_dir / f"deconv_{cohort}_{ref_name}.csv"
            if p.exists():
                df_f = pd.read_csv(p, index_col=0)
                df_f.index = df_f.index.astype(str)
                fractions_map[cohort] = df_f
                
        data_merged_map = {}
        for cohort in iatlas_cohorts:
            if cohort in fractions_map and cohort in clinical_map:
                df_f = fractions_map[cohort]
                df_c = clinical_map[cohort]
                common = df_f.index.intersection(df_c.index)
                if len(common) > 0:
                    df_merged = df_f.loc[common].copy()
                    df_merged['response'] = df_c.loc[common, 'response']
                    data_merged_map[cohort] = df_merged
                    
        for comb_name, sub_cohorts in iatlas_cancer_map.items():
            sub_dfs = [data_merged_map[c] for c in sub_cohorts if c in data_merged_map]
            if sub_dfs:
                data_merged_map[comb_name] = pd.concat(sub_dfs, axis=0)
                
        ordered_eval_cohorts = [
            "Hugo-iAtlas", "Riaz-iAtlas", "Liu-iAtlas", "Gide-iAtlas", "Combined-Melanoma",
            "Rosenberg-iAtlas", "Combined-BLCA",
            "Padron-iAtlas", "Combined-PAAD",
            "Anders-iAtlas", "Combined-BRCA",
            "McDermott-iAtlas", "Choueiri-iAtlas", "Combined-RCC",
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
                res_dict['Reference'] = ref_name
                res_dict['Cohort'] = c_name
                res_dict['CellType'] = ct
                results_list.append(res_dict)
                
        df_results = pd.DataFrame(results_list)
        csv_path = output_dir / f"lrt_summary_table_{ref_name}.csv"
        df_results.to_csv(csv_path, index=False)
        print(f"Saved comprehensive multi-test summary table to {csv_path}")
        
        # -------------------------------------------------------------
        # 2. OVERVIEW HEATMAP (LRT p-values & Odds Ratios)
        # -------------------------------------------------------------
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
                    
        plt.figure(figsize=(18, 11), dpi=300)
        sns.heatmap(
            beta_pivot,
            cmap="vlag",
            center=0,
            annot=annot_matrix,
            fmt="",
            cbar_kws={'label': 'Log Odds Ratio (beta > 0 favors Responder R)'},
            linewidths=0.8,
            linecolor="white",
            annot_kws={"size": 8, "fontweight": "bold"}
        )
        plt.title(f"Likelihood Ratio Test (LRT) & Odds Ratios for Immunotherapy Response: {ref_name}\nLogit(P(Responder)) = beta_0 + beta_1 * Cell_Fraction", fontsize=14, fontweight='bold', y=1.02)
        plt.xlabel("iAtlas & Combined Cohorts", fontweight='bold')
        plt.ylabel("LM22 Cell Types (No IG Genes)", fontweight='bold')
        plt.xticks(rotation=45, ha='right', fontweight='bold')
        
        heatmap_path = output_dir / f"lrt_heatmap_{ref_name}.png"
        plt.savefig(heatmap_path, bbox_inches='tight')
        plt.savefig(output_dir / f"lrt_heatmap_{ref_name}.svg", bbox_inches='tight')
        plt.close()
        print(f"Saved LRT heatmap to {heatmap_path}")
        
        # -------------------------------------------------------------
        # 3. MULTI-PANEL RESPONDER VS NON-RESPONDER BOXPLOTS WITH ALL 4 TESTS
        # -------------------------------------------------------------
        for c_name in eval_cohorts:
            df_cohort = data_merged_map[c_name]
            df_sub_res = df_results[df_results['Cohort'] == c_name].set_index('CellType')
            
            n_ct = len(cell_types)
            n_cols = 4
            n_rows = int(np.ceil(n_ct / n_cols))
            
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
                    ci_low = row_info['or_ci_low']
                    ci_high = row_info['or_ci_high']
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
                f"Cell Type Fraction Response Comparisons & Multi-Test Statistical Evaluation: {c_name} ({ref_name})\n"
                f"Total Samples N={len(df_cohort)} (Responders R={n_r_total}, Non-Responders NR={n_nr_total}) | Y-limits scaled by highest candle",
                fontsize=15, fontweight='bold', y=1.01
            )
            fig.tight_layout()
            
            cohort_clean = c_name.lower().replace("-", "_")
            out_box_path = boxplots_dir / f"response_boxplots_{cohort_clean}_{ref_name}.png"
            fig.savefig(out_box_path, bbox_inches='tight')
            fig.savefig(boxplots_dir / f"response_boxplots_{cohort_clean}_{ref_name}.svg", bbox_inches='tight')
            plt.close(fig)
            
        print(f"  Saved multi-test responder boxplots for all cohorts in {boxplots_dir}")

    print("\nMulti-test statistical evaluation completed successfully.")


if __name__ == "__main__":
    main()
