#!/usr/bin/env python3
"""
Calculate Spearman, Pearson, and Hilbert-Schmidt Independence Criterion (HSIC) 
correlations across multiple kernel bandwidths between inferred single-cell reference 
deconvolution fractions and clinical response in iAtlas cohorts.
Generates statistical summary tables, multi-bandwidth heatmaps, and cell cluster RKHS dependency curves.
"""

import os
import sys
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform

# Ensure local packages are in path
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))
import datalair
from single_cell_datasets._single_cell_reference import SingleCellDeconvolution

# Set publication style
sns.set_theme(style="whitegrid")
plt.rcParams.update({
    'font.size': 11,
    'figure.autolayout': True
})


def compute_rbf_kernel(X, bandwidth_factor=1.0, gamma=None):
    """Compute RBF kernel matrix using median heuristic for gamma scaled by bandwidth_factor."""
    X = np.atleast_2d(X)
    if X.shape[0] == 1:
        X = X.T
    dists_sq = squareform(pdist(X, metric="sqeuclidean"))
    if gamma is None:
        dists = np.sqrt(dists_sq)
        median_dist = np.median(dists[dists > 0])
        if median_dist == 0 or np.isnan(median_dist):
            median_dist = 1.0
        sigma = median_dist * bandwidth_factor
        gamma = 1.0 / (2.0 * sigma**2)
    K = np.exp(-gamma * dists_sq)
    return K, gamma


def compute_hsic(X, Y, bandwidth_factor=1.0, n_permutations=100, seed=42):
    """
    Compute Hilbert-Schmidt Independence Criterion (HSIC) statistic and permutation test p-value.
    Permutation resampling evaluates the empirical null distribution H0, preventing false positive inflation.
    """
    if seed is not None:
        np.random.seed(seed)
    X = np.ravel(X)[:, np.newaxis]
    Y = np.ravel(Y)[:, np.newaxis]
    N = len(X)
    if N < 4:
        return 0.0, 1.0

    K, _ = compute_rbf_kernel(X, bandwidth_factor=bandwidth_factor)
    L, _ = compute_rbf_kernel(Y, bandwidth_factor=bandwidth_factor)
    H = np.eye(N) - (1.0 / N) * np.ones((N, N))
    Kc = H @ K @ H
    Lc = H @ L @ H

    hsic_stat = float(np.sum(Kc * Lc) / ((N - 1) ** 2))

    if n_permutations > 0:
        perm_stats = np.zeros(n_permutations)
        for b in range(n_permutations):
            perm_idx = np.random.permutation(N)
            perm_stats[b] = np.sum(Kc * Lc[perm_idx, :][:, perm_idx]) / ((N - 1) ** 2)
        p_val = float((1.0 + np.sum(perm_stats >= hsic_stat)) / (n_permutations + 1.0))
    else:
        p_val = 1.0

    return hsic_stat, p_val


def compute_sample_influence(X, Y, bandwidth_factor=1.0):
    """Compute sample-wise HSIC influence vector h_i."""
    X = np.ravel(X)[:, np.newaxis]
    Y = np.ravel(Y)[:, np.newaxis]
    N = len(X)
    K, _ = compute_rbf_kernel(X, bandwidth_factor=bandwidth_factor)
    L, _ = compute_rbf_kernel(Y, bandwidth_factor=bandwidth_factor)
    H = np.eye(N) - (1.0 / N) * np.ones((N, N))
    Kc = H @ K @ H
    Lc = H @ L @ H
    C = Kc * Lc
    h_i = np.sum(C, axis=1) / (N - 1)
    return h_i


def compute_rkhs_conditional_moments(X, Y, n_grid=150, bandwidth_factor=0.25):
    """Compute RKHS conditional expectation E[Y | X = x] and conditional SD across grid."""
    X_vec = np.ravel(X)
    Y_vec = np.ravel(Y)
    x_grid = np.linspace(np.min(X_vec), np.max(X_vec), n_grid)
    
    dists_sq = squareform(pdist(X_vec[:, np.newaxis], metric="sqeuclidean"))
    dists = np.sqrt(dists_sq)
    med = np.median(dists[dists > 0])
    sigma = (med if med > 0 else 1.0) * bandwidth_factor
    gamma_x = 1.0 / (2.0 * sigma**2)

    cond_mean = np.zeros(n_grid)
    cond_sd = np.zeros(n_grid)

    for idx, x_val in enumerate(x_grid):
        w = np.exp(-gamma_x * (x_val - X_vec) ** 2)
        sum_w = np.sum(w)
        if sum_w == 0:
            sum_w = 1e-12
        w_norm = w / sum_w

        mu = np.sum(w_norm * Y_vec)
        var = np.sum(w_norm * (Y_vec - mu) ** 2)

        cond_mean[idx] = mu
        cond_sd[idx] = np.sqrt(var)

    return x_grid, cond_mean, cond_sd


def load_clinical_response(cohort_name: str, lair):
    """Helper to load clinical response from CBioPortal."""
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
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output" / "single-cell-exploration"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get deconvolution dataset filepaths
    ds_deconv = SingleCellDeconvolution()
    deconv_paths = lair.get_dataset_filepaths(ds_deconv)
    
    cohorts = [
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
    
    resolutions = ["leiden_res_0.5", "kmeans_subcluster_res_0.5"]
    
    # Pre-load clinical data
    clinical_data = {}
    for cohort in cohorts:
        df_clin = load_clinical_response(cohort, lair)
        if df_clin is not None:
            clinical_data[cohort] = df_clin
            print(f"Loaded clinical response for {cohort}: {len(df_clin)} samples ({sum(df_clin['response']=='R')} R, {sum(df_clin['response']=='NR')} NR)")
            
    correlation_rows = []
    combined_data_store = {}

    for resolution in resolutions:
        print(f"\n--- Analyzing correlations for {resolution} ---")
        
        for cohort in cohorts:
            key = f"deconv_{cohort}_{resolution}.csv"
            path = deconv_paths.get(key)
            df_clin = clinical_data.get(cohort)
            
            if not path or df_clin is None:
                continue
                
            df_fracs = pd.read_csv(path, index_col=0)
            df_fracs.index = df_fracs.index.astype(str)
            
            common_samples = df_fracs.index.intersection(df_clin.index)
            if len(common_samples) < 10:
                continue
                
            X_df = df_fracs.loc[common_samples]
            y_series = df_clin.loc[common_samples, 'response'].map({'R': 1, 'NR': 0})
            
            for cell_type in X_df.columns:
                vals = X_df[cell_type].values
                y = y_series.values
                
                if np.var(vals) == 0:
                    continue
                pearson_r, pearson_p = stats.pearsonr(vals, y)
                spearman_r, spearman_p = stats.spearmanr(vals, y)
                
                nr_vals = vals[y == 0]
                r_vals = vals[y == 1]
                if len(nr_vals) > 0 and len(r_vals) > 0:
                    _, mwu_p = stats.mannwhitneyu(nr_vals, r_vals, alternative='two-sided')
                else:
                    mwu_p = 1.0
                
                # HSIC at 4 bandwidth factors
                hsic10, hsic_p10 = compute_hsic(vals, y, bandwidth_factor=1.0)
                hsic05, hsic_p05 = compute_hsic(vals, y, bandwidth_factor=0.5)
                hsic025, hsic_p025 = compute_hsic(vals, y, bandwidth_factor=0.25)
                hsic0125, hsic_p0125 = compute_hsic(vals, y, bandwidth_factor=0.125)
                    
                correlation_rows.append({
                    'Resolution': resolution,
                    'Cohort': cohort,
                    'CellType': cell_type,
                    'Pearson_R': pearson_r,
                    'Pearson_P': pearson_p,
                    'Spearman_R': spearman_r,
                    'Spearman_P': spearman_p,
                    'MWU_P': mwu_p,
                    'HSIC': hsic10,
                    'HSIC_P': hsic_p10,
                    'HSIC_0.5': hsic05,
                    'HSIC_P_0.5': hsic_p05,
                    'HSIC_0.25': hsic025,
                    'HSIC_P_0.25': hsic_p025,
                    'HSIC_0.125': hsic0125,
                    'HSIC_P_0.125': hsic_p0125,
                    'IsSignificant': spearman_p < 0.05 or pearson_p < 0.05 or hsic_p10 < 0.05
                })

        # Combined Cohorts
        combined_specs = [
            ('Combined-Melanoma', ['Hugo-iAtlas', 'Riaz-iAtlas', 'Liu-iAtlas', 'Gide-iAtlas']),
            ('Combined-RCC', ['McDermott-iAtlas', 'Choueiri-iAtlas']),
            ('Combined-All', ['Hugo-iAtlas', 'Riaz-iAtlas', 'Liu-iAtlas', 'Gide-iAtlas', 'Rosenberg-iAtlas', 'Padron-iAtlas', 'Anders-iAtlas', 'McDermott-iAtlas', 'Choueiri-iAtlas'])
        ]
        
        for comb_name, sub_cohorts in combined_specs:
            fracs_list = []
            y_list = []
            for c in sub_cohorts:
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
            
            if not fracs_list:
                continue
                
            df_fracs_comb = pd.concat(fracs_list, axis=0)
            y_series_comb = pd.concat(y_list, axis=0)
            
            combined_data_store[(resolution, comb_name)] = (df_fracs_comb, y_series_comb)
            
            for cell_type in df_fracs_comb.columns:
                vals = df_fracs_comb[cell_type].values
                y = y_series_comb.values
                if np.var(vals) == 0:
                    continue
                pearson_r, pearson_p = stats.pearsonr(vals, y)
                spearman_r, spearman_p = stats.spearmanr(vals, y)
                
                nr_vals = vals[y == 0]
                r_vals = vals[y == 1]
                if len(nr_vals) > 0 and len(r_vals) > 0:
                    _, mwu_p = stats.mannwhitneyu(nr_vals, r_vals, alternative='two-sided')
                else:
                    mwu_p = 1.0
                
                hsic10, hsic_p10 = compute_hsic(vals, y, bandwidth_factor=1.0)
                hsic05, hsic_p05 = compute_hsic(vals, y, bandwidth_factor=0.5)
                hsic025, hsic_p025 = compute_hsic(vals, y, bandwidth_factor=0.25)
                hsic0125, hsic_p0125 = compute_hsic(vals, y, bandwidth_factor=0.125)

                correlation_rows.append({
                    'Resolution': resolution,
                    'Cohort': comb_name,
                    'CellType': cell_type,
                    'Pearson_R': pearson_r,
                    'Pearson_P': pearson_p,
                    'Spearman_R': spearman_r,
                    'Spearman_P': spearman_p,
                    'MWU_P': mwu_p,
                    'HSIC': hsic10,
                    'HSIC_P': hsic_p10,
                    'HSIC_0.5': hsic05,
                    'HSIC_P_0.5': hsic_p05,
                    'HSIC_0.25': hsic025,
                    'HSIC_P_0.25': hsic_p025,
                    'HSIC_0.125': hsic0125,
                    'HSIC_P_0.125': hsic_p0125,
                    'IsSignificant': spearman_p < 0.05 or pearson_p < 0.05 or hsic_p10 < 0.05
                })

    df_corrs = pd.DataFrame(correlation_rows)
    df_corrs.to_csv(output_dir / "self_ref_deconv_clinical_correlations.csv", index=False)
    print(f"Saved correlation summary to output/single-cell-exploration/self_ref_deconv_clinical_correlations.csv")

    CELL_TYPE_NAMES_LEIDEN_05 = {
        '0': 'CD4+ Helper T (c0)', '1': 'Epithelial Tumor c1', '2': 'M1 Macrophage c2',
        '3': 'CD8+ Cytotoxic T c3', '4': 'CD4+ Memory T c4', '5': 'B Cells c5',
        '6': 'Endothelial c6', '7': 'CAFs c7', '8': 'NK Cells c8', '9': 'Monocytes c9',
        '10': 'Plasma c10', '11': 'pDCs c11', '12': 'mDCs c12', '13': 'Tregs c13',
        '14': 'Mast Cells c14', '15': 'Neutrophils c15', '16': 'GammaDelta T c16',
        '17': 'Pericytes c17', '18': 'Myofibroblasts c18', '19': 'Erythrocytes c19',
        '20': 'Granulocytes c20', '21': 'Melanocytes c21'
    }

    # Generate Heatmaps for each resolution
    for resolution in resolutions:
        df_res = df_corrs[df_corrs['Resolution'] == resolution].copy()
        
        cohort_order = [
            'Hugo-iAtlas', 'Riaz-iAtlas', 'Liu-iAtlas', 'Gide-iAtlas', 
            'Rosenberg-iAtlas', 'Padron-iAtlas', 'Anders-iAtlas', 
            'McDermott-iAtlas', 'Choueiri-iAtlas',
            'Combined-Melanoma', 'Combined-RCC', 'Combined-All'
        ]
        existing_cohorts = [c for c in cohort_order if c in df_res['Cohort'].unique()]
        
        if resolution == "leiden_res_0.5":
            df_res['CellType_Name'] = df_res['CellType'].map(lambda x: CELL_TYPE_NAMES_LEIDEN_05.get(str(x), str(x)))
            new_index = [CELL_TYPE_NAMES_LEIDEN_05.get(str(i), str(i)) for i in range(22)]
        elif resolution == "kmeans_subcluster_res_0.5":
            df_res['CellType_Name'] = df_res['CellType']
            def parse_cluster_key(val):
                parts = str(val).split('_sub_')
                return (int(parts[0]), int(parts[1])) if len(parts) == 2 else (999, 999)
            unique_cats = df_res['CellType'].unique()
            sorted_cats = sorted(unique_cats, key=parse_cluster_key)
            
        # --- Spearman Heatmap ---
        pivot_r_sp = df_res.pivot(index='CellType', columns='Cohort', values='Spearman_R')[existing_cohorts]
        pivot_p_sp = df_res.pivot(index='CellType', columns='Cohort', values='Spearman_P')[existing_cohorts]
        if resolution == "leiden_res_0.5":
            pivot_r_sp.index = pd.Categorical(pivot_r_sp.index, categories=[str(i) for i in range(22)], ordered=True)
            pivot_r_sp = pivot_r_sp.sort_index(); pivot_r_sp.index = new_index
            pivot_p_sp.index = pd.Categorical(pivot_p_sp.index, categories=[str(i) for i in range(22)], ordered=True)
            pivot_p_sp = pivot_p_sp.sort_index(); pivot_p_sp.index = new_index
        elif resolution == "kmeans_subcluster_res_0.5":
            pivot_r_sp.index = pd.Categorical(pivot_r_sp.index, categories=sorted_cats, ordered=True); pivot_r_sp = pivot_r_sp.sort_index()
            pivot_p_sp.index = pd.Categorical(pivot_p_sp.index, categories=sorted_cats, ordered=True); pivot_p_sp = pivot_p_sp.sort_index()
            
        annot_sp = pd.DataFrame("", index=pivot_r_sp.index, columns=pivot_r_sp.columns)
        for r_idx in pivot_r_sp.index:
            for c_idx in pivot_r_sp.columns:
                p_val = pivot_p_sp.loc[r_idx, c_idx]
                if not pd.isna(p_val):
                    annot_sp.loc[r_idx, c_idx] = "***" if p_val < 0.001 else ("**" if p_val < 0.01 else ("*" if p_val < 0.05 else ""))
                    
        plt.figure(figsize=(14, 11), dpi=300)
        sns.heatmap(pivot_r_sp, annot=annot_sp, fmt="", cmap="vlag", center=0, cbar_kws={'label': "Spearman Correlation Coefficient (R vs NR)"}, linewidths=0.5, edgecolor='black')
        plt.title(f"Spearman Clinical Response Correlation: Single-Cell Reference Deconvolution ({resolution})\n(* p < 0.05, ** p < 0.01, *** p < 0.001)", fontsize=13, fontweight='bold', y=1.02)
        plt.xlabel("Cohort", fontweight='bold'); plt.ylabel("Inferred Cell Type / Cluster", fontweight='bold'); plt.tight_layout()
        plt.savefig(output_dir / f"self_ref_response_correlation_spearman_{resolution}.svg", bbox_inches='tight')
        plt.savefig(output_dir / f"self_ref_response_correlation_{resolution}.svg", bbox_inches='tight')
        plt.close()

        # --- Pearson Heatmap ---
        pivot_r_pe = df_res.pivot(index='CellType', columns='Cohort', values='Pearson_R')[existing_cohorts]
        pivot_p_pe = df_res.pivot(index='CellType', columns='Cohort', values='Pearson_P')[existing_cohorts]
        if resolution == "leiden_res_0.5":
            pivot_r_pe.index = pd.Categorical(pivot_r_pe.index, categories=[str(i) for i in range(22)], ordered=True); pivot_r_pe = pivot_r_pe.sort_index(); pivot_r_pe.index = new_index
            pivot_p_pe.index = pd.Categorical(pivot_p_pe.index, categories=[str(i) for i in range(22)], ordered=True); pivot_p_pe = pivot_p_pe.sort_index(); pivot_p_pe.index = new_index
        elif resolution == "kmeans_subcluster_res_0.5":
            pivot_r_pe.index = pd.Categorical(pivot_r_pe.index, categories=sorted_cats, ordered=True); pivot_r_pe = pivot_r_pe.sort_index()
            pivot_p_pe.index = pd.Categorical(pivot_p_pe.index, categories=sorted_cats, ordered=True); pivot_p_pe = pivot_p_pe.sort_index()
            
        annot_pe = pd.DataFrame("", index=pivot_r_pe.index, columns=pivot_r_pe.columns)
        for r_idx in pivot_r_pe.index:
            for c_idx in pivot_r_pe.columns:
                p_val = pivot_p_pe.loc[r_idx, c_idx]
                if not pd.isna(p_val):
                    annot_pe.loc[r_idx, c_idx] = "***" if p_val < 0.001 else ("**" if p_val < 0.01 else ("*" if p_val < 0.05 else ""))
                    
        plt.figure(figsize=(14, 11), dpi=300)
        sns.heatmap(pivot_r_pe, annot=annot_pe, fmt="", cmap="vlag", center=0, cbar_kws={'label': "Pearson Correlation Coefficient (R vs NR)"}, linewidths=0.5, edgecolor='black')
        plt.title(f"Pearson Clinical Response Correlation: Single-Cell Reference Deconvolution ({resolution})\n(* p < 0.05, ** p < 0.01, *** p < 0.001)", fontsize=13, fontweight='bold', y=1.02)
        plt.xlabel("Cohort", fontweight='bold'); plt.ylabel("Inferred Cell Type / Cluster", fontweight='bold'); plt.tight_layout()
        plt.savefig(output_dir / f"self_ref_response_correlation_pearson_{resolution}.svg", bbox_inches='tight')
        plt.close()

        # --- Multi-Bandwidth HSIC Heatmaps (2x2 Combined Grid Figure) ---
        bw_col_maps = [('HSIC', 'HSIC_P', 'Standard (1.0 × σ_med)'),
                       ('HSIC_0.5', 'HSIC_P_0.5', 'Intermediate (0.5 × σ_med)'),
                       ('HSIC_0.25', 'HSIC_P_0.25', 'Fine Local (0.25 × σ_med)'),
                       ('HSIC_0.125', 'HSIC_P_0.125', 'Ultra-Fine Local (1/8 × σ_med)')]

        fig, axes = plt.subplots(2, 2, figsize=(22, 18), dpi=300)

        for ax_idx, (r_col, p_col, title_lbl) in enumerate(bw_col_maps):
            row_i, col_j = divmod(ax_idx, 2)
            ax = axes[row_i, col_j]

            pivot_r_h = df_res.pivot(index='CellType', columns='Cohort', values=r_col)[existing_cohorts]
            pivot_p_h = df_res.pivot(index='CellType', columns='Cohort', values=p_col)[existing_cohorts]

            if resolution == "leiden_res_0.5":
                pivot_r_h.index = pd.Categorical(pivot_r_h.index, categories=[str(i) for i in range(22)], ordered=True); pivot_r_h = pivot_r_h.sort_index(); pivot_r_h.index = new_index
                pivot_p_h.index = pd.Categorical(pivot_p_h.index, categories=[str(i) for i in range(22)], ordered=True); pivot_p_h = pivot_p_h.sort_index(); pivot_p_h.index = new_index
            elif resolution == "kmeans_subcluster_res_0.5":
                pivot_r_h.index = pd.Categorical(pivot_r_h.index, categories=sorted_cats, ordered=True); pivot_r_h = pivot_r_h.sort_index()
                pivot_p_h.index = pd.Categorical(pivot_p_h.index, categories=sorted_cats, ordered=True); pivot_p_h = pivot_p_h.sort_index()

            annot_h = pd.DataFrame("", index=pivot_r_h.index, columns=pivot_r_h.columns)
            for r_idx in pivot_r_h.index:
                for c_idx in pivot_r_h.columns:
                    p_val = pivot_p_h.loc[r_idx, c_idx]
                    if not pd.isna(p_val):
                        annot_h.loc[r_idx, c_idx] = "***" if p_val < 0.001 else ("**" if p_val < 0.01 else ("*" if p_val < 0.05 else ""))

            sns.heatmap(pivot_r_h, annot=annot_h, fmt="", cmap="OrRd", cbar_kws={'label': "HSIC Statistic"}, linewidths=0.5, edgecolor='black', ax=ax)
            ax.set_title(f"HSIC Dependence ({title_lbl})\n(* p < 0.05, ** p < 0.01, *** p < 0.001)", fontsize=12, fontweight='bold')
            ax.set_xlabel("Cohort", fontweight='bold')
            ax.set_ylabel("Cell Type / Cluster", fontweight='bold')

            # Save individual bandwidth heatmaps
            plt_single, ax_s = plt.subplots(figsize=(14, 11), dpi=300)
            sns.heatmap(pivot_r_h, annot=annot_h, fmt="", cmap="OrRd", cbar_kws={'label': "HSIC Statistic"}, linewidths=0.5, edgecolor='black', ax=ax_s)
            ax_s.set_title(f"HSIC Response Dependence ({title_lbl}) - {resolution}\n(* p < 0.05, ** p < 0.01, *** p < 0.001)", fontsize=13, fontweight='bold')
            ax_s.set_xlabel("Cohort", fontweight='bold'); ax_s.set_ylabel("Inferred Cell Type / Cluster", fontweight='bold'); plt_single.tight_layout()
            plt_single.savefig(output_dir / f"self_ref_response_correlation_hsic_bw_{r_col}_{resolution}.svg", bbox_inches='tight')
            plt_single.savefig(output_dir / f"self_ref_response_correlation_hsic_bw_{r_col}_{resolution}.png", bbox_inches='tight')
            plt.close(plt_single)

        fig.suptitle(f"Multi-Bandwidth HSIC Clinical Response Dependence ({resolution})", fontsize=15, fontweight='bold')
        fig.tight_layout()
        multi_bw_path = output_dir / f"self_ref_hsic_heatmaps_multibandwidth_{resolution}.png"
        fig.savefig(multi_bw_path, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved multi-bandwidth HSIC heatmap grid to {multi_bw_path}")

        # Also save standard HSIC copy for backwards compatibility
        shutil.copy2(output_dir / f"self_ref_response_correlation_hsic_bw_HSIC_{resolution}.png", output_dir / f"self_ref_response_correlation_hsic_{resolution}.png")
        shutil.copy2(output_dir / f"self_ref_response_correlation_hsic_bw_HSIC_{resolution}.svg", output_dir / f"self_ref_response_correlation_hsic_{resolution}.svg")

    # --- Cell Cluster HSIC Dependency Visualizations (Scatter Influence + RKHS Response Curves) ---
    print("\nGenerating cell cluster HSIC dependency plots across bandwidths for top clusters...")
    
    # Pick top cell clusters for leiden_res_0.5 in Combined-All
    res_sel = "leiden_res_0.5"
    if (res_sel, "Combined-All") in combined_data_store:
        df_fracs_comb, y_series_comb = combined_data_store[(res_sel, "Combined-All")]
        X_mat = df_fracs_comb
        Y_vec = y_series_comb.values
        
        target_clusters = ['3', '2', '0', '1']
        bandwidth_factors = [1.0, 0.5, 0.25, 0.125]
        bw_names = ["1.0 × σ_med", "0.5 × σ_med", "0.25 × σ_med", "0.125 × σ_med (1/8)"]

        for c_id in target_clusters:
            if c_id not in X_mat.columns:
                continue
            c_name = CELL_TYPE_NAMES_LEIDEN_05.get(str(c_id), f"Cluster {c_id}")
            vals = X_mat[c_id].values
            
            fig, axes = plt.subplots(4, 2, figsize=(14, 15), sharex="col")
            
            for r_i, (bw, bw_lbl) in enumerate(zip(bandwidth_factors, bw_names)):
                # Left: Sample Influence Scatter
                h_inf = compute_sample_influence(vals, Y_vec, bandwidth_factor=bw)
                hsic_stat, hsic_p = compute_hsic(vals, Y_vec, bandwidth_factor=bw)
                
                sc = axes[r_i, 0].scatter(vals, Y_vec, c=h_inf, cmap="magma", s=25, alpha=0.85)
                cbar = fig.colorbar(sc, ax=axes[r_i, 0])
                cbar.set_label("hᵢ", fontweight="bold", fontsize=9)
                axes[r_i, 0].set_ylabel("Clinical Response (NR=0, R=1)", fontweight="bold", fontsize=9)
                axes[r_i, 0].set_title(f"Sample HSIC Influence (hᵢ) [{bw_lbl}]\nHSIC = {hsic_stat:.4f} (p = {hsic_p:.4f})", fontsize=10, fontweight="bold")
                
                # Right: RKHS Non-parametric Conditional Response Probability Curve E[Y|X]
                x_grid, mu_y, sd_y = compute_rkhs_conditional_moments(vals, Y_vec, bandwidth_factor=bw)
                axes[r_i, 1].scatter(vals, Y_vec, color="#2b5c8f", alpha=0.3, s=15, label="Samples")
                axes[r_i, 1].plot(x_grid, mu_y, "#d95f02", linewidth=2.5, label="RKHS P(Response | Fraction)")
                axes[r_i, 1].fill_between(x_grid, np.clip(mu_y - sd_y, 0, 1), np.clip(mu_y + sd_y, 0, 1), color="#d95f02", alpha=0.2, label="±1 SD Envelope")
                axes[r_i, 1].set_ylabel("P(Response | Fraction)", fontweight="bold", fontsize=9)
                axes[r_i, 1].set_title(f"RKHS Response Curve [{bw_lbl}]", fontsize=10, fontweight="bold")
                axes[r_i, 1].legend(loc="upper right", fontsize=8)
                
                if r_i == 3:
                    axes[r_i, 0].set_xlabel(f"Inferred {c_name} Fraction", fontweight="bold")
                    axes[r_i, 1].set_xlabel(f"Inferred {c_name} Fraction", fontweight="bold")

            fig.suptitle(f"HSIC Non-Linear Response Dependency: {c_name} (Combined-All Cohorts)\nComparison across 4 Kernel Bandwidths (1.0x, 0.5x, 0.25x, 0.125x)", fontsize=13, fontweight="bold")
            fig.tight_layout()
            
            clean_cname = c_name.lower().replace(" ", "_").replace("+", "").replace("(", "").replace(")", "").replace("/", "_")
            out_cluster_path = output_dir / f"self_ref_hsic_cluster_dependency_{clean_cname}_{res_sel}.png"
            fig.savefig(out_cluster_path, dpi=300)
            plt.close(fig)
            print(f"Saved cell cluster HSIC dependency plot to {out_cluster_path}")

    print("\nAll iAtlas self-reference response correlation & HSIC analyses completed successfully.")


if __name__ == "__main__":
    main()
