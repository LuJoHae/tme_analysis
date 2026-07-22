#!/usr/bin/env python3
"""
Calculate Spearman correlation and Mann-Whitney U test p-values between inferred cell type fractions
from our self-constructed single-cell reference deconvolution and immunotherapy response in iAtlas cohorts.
Generates statistical summary tables and heatmaps.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import seaborn as sns

# Ensure local packages are in path
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))
import datalair
from single_cell_datasets._single_cell_reference import SingleCellDeconvolution


def compute_rbf_kernel(X, gamma=None):
    """Compute RBF kernel matrix using median heuristic for gamma."""
    X = np.atleast_2d(X)
    if X.shape[0] == 1:
        X = X.T
    dists_sq = squareform(pdist(X, metric="sqeuclidean"))
    if gamma is None:
        dists = np.sqrt(dists_sq)
        median_dist = np.median(dists[dists > 0])
        if median_dist == 0 or np.isnan(median_dist):
            median_dist = 1.0
        sigma = median_dist
        gamma = 1.0 / (2.0 * sigma**2)
    K = np.exp(-gamma * dists_sq)
    return K, gamma


def compute_hsic(X, Y, n_permutations=200, seed=42):
    """
    Compute Hilbert-Schmidt Independence Criterion (HSIC) statistic and p-value.
    Uses asymptotic Gamma distribution approximation under H0 for fast p-value calculation.
    """
    if seed is not None:
        np.random.seed(seed)
    X = np.ravel(X)[:, np.newaxis]
    Y = np.ravel(Y)[:, np.newaxis]
    N = len(X)
    if N < 4:
        return 0.0, 1.0

    K, _ = compute_rbf_kernel(X)
    L, _ = compute_rbf_kernel(Y)
    H = np.eye(N) - (1.0 / N) * np.ones((N, N))
    Kc = H @ K @ H
    Lc = H @ L @ H

    hsic_stat = float(np.trace(Kc @ Lc) / ((N - 1) ** 2))

    if n_permutations > 0:
        perm_stats = np.zeros(n_permutations)
        for b in range(n_permutations):
            perm_idx = np.random.permutation(N)
            perm_stats[b] = np.sum(Kc * L[perm_idx, :][:, perm_idx]) / ((N - 1) ** 2)

        mu_null = float(np.mean(perm_stats))
        var_null = float(np.var(perm_stats))
        if var_null > 0 and hsic_stat > 0:
            shape_k = (mu_null ** 2) / var_null
            scale_theta = var_null / mu_null
            p_val = float(stats.gamma.sf(hsic_stat, a=shape_k, scale=scale_theta))
        else:
            p_val = float((1 + np.sum(perm_stats >= hsic_stat)) / (n_permutations + 1))
    else:
        p_val = 1.0

    return hsic_stat, p_val


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
            
    correlation_rows = []
    
    for resolution in resolutions:
        print(f"\n--- Analyzing correlations for {resolution} ---")
        
        for cohort in cohorts:
            key = f"deconv_{cohort}_{resolution}.csv"
            path = deconv_paths.get(key)
            df_clin = clinical_data.get(cohort)
            
            if not path or df_clin is None:
                continue
                
            # Load deconvolution fractions
            df_fracs = pd.read_csv(path, index_col=0)
            df_fracs.index = df_fracs.index.astype(str)
            
            # Align
            common_samples = df_fracs.index.intersection(df_clin.index)
            if len(common_samples) < 10:
                print(f"Skipping {cohort}: Too few samples ({len(common_samples)})")
                continue
                
            X_df = df_fracs.loc[common_samples]
            y_series = df_clin.loc[common_samples, 'response'].map({'R': 1, 'NR': 0})
            
            # Loop over each cell type (column)
            for cell_type in X_df.columns:
                vals = X_df[cell_type].values
                y = y_series.values
                
                # Check for zero variance
                if np.var(vals) == 0:
                    continue
                # Pearson correlation
                pearson_r, pearson_p = stats.pearsonr(vals, y)
                
                # Spearman correlation
                spearman_r, spearman_p = stats.spearmanr(vals, y)
                
                # Mann-Whitney U test p-value (NR vs R)
                nr_vals = vals[y == 0]
                r_vals = vals[y == 1]
                if len(nr_vals) > 0 and len(r_vals) > 0:
                    _, mwu_p = stats.mannwhitneyu(nr_vals, r_vals, alternative='two-sided')
                else:
                    mwu_p = 1.0
                
                # Hilbert-Schmidt Independence Criterion (HSIC)
                hsic_stat, hsic_p = compute_hsic(vals, y)
                    
                correlation_rows.append({
                    'Resolution': resolution,
                    'Cohort': cohort,
                    'CellType': cell_type,
                    'Pearson_R': pearson_r,
                    'Pearson_P': pearson_p,
                    'Spearman_R': spearman_r,
                    'Spearman_P': spearman_p,
                    'MWU_P': mwu_p,
                    'HSIC': hsic_stat,
                    'HSIC_P': hsic_p,
                    'IsSignificant': spearman_p < 0.05 or pearson_p < 0.05 or hsic_p < 0.05
                })
        # Calculate for combined cohorts
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
            
            for cell_type in df_fracs_comb.columns:
                vals = df_fracs_comb[cell_type].values
                y = y_series_comb.values
                if np.var(vals) == 0:
                    continue
                # Pearson correlation
                pearson_r, pearson_p = stats.pearsonr(vals, y)
                
                # Spearman correlation
                spearman_r, spearman_p = stats.spearmanr(vals, y)
                
                nr_vals = vals[y == 0]
                r_vals = vals[y == 1]
                if len(nr_vals) > 0 and len(r_vals) > 0:
                    _, mwu_p = stats.mannwhitneyu(nr_vals, r_vals, alternative='two-sided')
                else:
                    mwu_p = 1.0
                
                # Hilbert-Schmidt Independence Criterion (HSIC)
                hsic_stat, hsic_p = compute_hsic(vals, y)

                correlation_rows.append({
                    'Resolution': resolution,
                    'Cohort': comb_name,
                    'CellType': cell_type,
                    'Pearson_R': pearson_r,
                    'Pearson_P': pearson_p,
                    'Spearman_R': spearman_r,
                    'Spearman_P': spearman_p,
                    'MWU_P': mwu_p,
                    'HSIC': hsic_stat,
                    'HSIC_P': hsic_p,
                    'IsSignificant': spearman_p < 0.05 or pearson_p < 0.05 or hsic_p < 0.05
                })

    df_corrs = pd.DataFrame(correlation_rows)
    df_corrs.to_csv(output_dir / "self_ref_deconv_clinical_correlations.csv", index=False)
    print(f"Saved correlation summary to output/single-cell-exploration/self_ref_deconv_clinical_correlations.csv")
    
    # Let's print out the statistically significant correlations!
    sig_df = df_corrs[df_corrs['IsSignificant'] == True]
    print("\n--- STATISTICALLY SIGNIFICANT CORRELATIONS (p < 0.05) ---")
    if not sig_df.empty:
        for idx, row in sig_df.iterrows():
            print(f"[{row['Resolution']}] {row['Cohort']:20s} | {row['CellType']:25s} | Spearman R = {row['Spearman_R']:.3f} | p = {row['Spearman_P']:.4f} | MWU p = {row['MWU_P']:.4f}")
    else:
        print("None found.")
        
    # Unique name mapping for leiden_res_0.5 clusters based on marker gene analysis
    CELL_TYPE_NAMES_LEIDEN_05 = {
        '0': 'CD4+ Helper T (c0)',
        '1': 'Epithelial Tumor c1',
        '2': 'M1 Macrophage c2',
        '3': 'CD8+ Cytotoxic T c3',
        '4': 'M1 Macrophage c4',
        '5': 'Naive B (c5)',
        '6': 'CD8+ Cytotoxic T / NK (c6)',
        '7': 'CAF c7',
        '8': 'Melanoma Tumor c8',
        '9': 'Melanoma Tumor c9',
        '10': 'Normal Fibroblast c10',
        '11': 'Epithelial Tumor c11',
        '12': 'Plasma Cell c12',
        '13': 'Epithelial Tumor c13',
        '14': 'Blood Endothelial c14',
        '15': 'CD8+ Cytotoxic T c15',
        '16': 'Epithelial Tumor c16',
        '17': 'Epithelial Tumor c17',
        '18': 'Plasma Cell c18',
        '19': 'Mast Cell c19',
        '20': 'Melanoma Tumor c20',
        '21': 'Epithelial Tumor c21'
    }

    # Generate Heatmaps for each resolution
    for resolution in resolutions:
        df_res = df_corrs[df_corrs['Resolution'] == resolution]
        if df_res.empty:
            continue
            
        all_cohorts_ordered = cohorts + ['Combined-Melanoma', 'Combined-RCC', 'Combined-All']
        
        # --- Spearman Heatmap ---
        pivot_r_sp = df_res.pivot(index='CellType', columns='Cohort', values='Spearman_R')
        pivot_p_sp = df_res.pivot(index='CellType', columns='Cohort', values='Spearman_P')
        
        existing_cohorts = [c for c in all_cohorts_ordered if c in pivot_r_sp.columns]
        pivot_r_sp = pivot_r_sp[existing_cohorts]
        pivot_p_sp = pivot_p_sp[existing_cohorts]
        
        # Sort index numerically if possible
        if resolution == "leiden_res_0.5":
            pivot_r_sp.index = pd.Categorical(pivot_r_sp.index, categories=[str(i) for i in range(22)], ordered=True)
            pivot_r_sp = pivot_r_sp.sort_index()
            pivot_p_sp.index = pd.Categorical(pivot_p_sp.index, categories=[str(i) for i in range(22)], ordered=True)
            pivot_p_sp = pivot_p_sp.sort_index()
            
            new_index = [CELL_TYPE_NAMES_LEIDEN_05[str(i)] for i in pivot_r_sp.index]
            pivot_r_sp.index = new_index
            pivot_p_sp.index = new_index
        elif resolution == "kmeans_subcluster_res_0.5":
            # Sort subclusters cleanly
            sorted_cats = sorted(pivot_r_sp.index, key=lambda x: [int(s) if s.isdigit() else s for s in x.split('_')])
            pivot_r_sp.index = pd.Categorical(pivot_r_sp.index, categories=sorted_cats, ordered=True)
            pivot_r_sp = pivot_r_sp.sort_index()
            pivot_p_sp.index = pd.Categorical(pivot_p_sp.index, categories=sorted_cats, ordered=True)
            pivot_p_sp = pivot_p_sp.sort_index()
        
        annot_matrix_sp = pd.DataFrame("", index=pivot_r_sp.index, columns=pivot_r_sp.columns)
        for r_idx in pivot_r_sp.index:
            for c_idx in pivot_r_sp.columns:
                p_val = pivot_p_sp.loc[r_idx, c_idx]
                r_val = pivot_r_sp.loc[r_idx, c_idx]
                if pd.isna(p_val) or pd.isna(r_val):
                    continue
                if p_val < 0.001:
                    annot_matrix_sp.loc[r_idx, c_idx] = "***"
                elif p_val < 0.01:
                    annot_matrix_sp.loc[r_idx, c_idx] = "**"
                elif p_val < 0.05:
                    annot_matrix_sp.loc[r_idx, c_idx] = "*"
                    
        plt.figure(figsize=(14, 11), dpi=300)
        sns.heatmap(
            pivot_r_sp,
            annot=annot_matrix_sp,
            fmt="",
            cmap="vlag",
            center=0,
            cbar_kws={'label': "Spearman Correlation Coefficient (R vs NR)"},
            linewidths=0.5,
            edgecolor='black'
        )
        plt.title(f"Spearman Clinical Response Correlation: Single-Cell Reference Deconvolution ({resolution})\n(* p < 0.05, ** p < 0.01, *** p < 0.001)", fontsize=13, fontweight='bold', y=1.02)
        plt.xlabel("Cohort", fontweight='bold')
        plt.ylabel("Inferred Cell Type / Cluster", fontweight='bold')
        plt.tight_layout()
        
        plot_path_sp = output_dir / f"self_ref_response_correlation_spearman_{resolution}.svg"
        plt.savefig(plot_path_sp, bbox_inches='tight')
        
        # Save additional copy without "spearman_" in name as requested by the user
        plot_path_alt = output_dir / f"self_ref_response_correlation_{resolution}.svg"
        plt.savefig(plot_path_alt, bbox_inches='tight')
        plt.close()
        print(f"Saved Spearman response correlation heatmap to {plot_path_sp} and {plot_path_alt}")
        
        # --- Pearson Heatmap ---
        pivot_r_pe = df_res.pivot(index='CellType', columns='Cohort', values='Pearson_R')
        pivot_p_pe = df_res.pivot(index='CellType', columns='Cohort', values='Pearson_P')
        
        pivot_r_pe = pivot_r_pe[existing_cohorts]
        pivot_p_pe = pivot_p_pe[existing_cohorts]
        
        if resolution == "leiden_res_0.5":
            pivot_r_pe.index = pd.Categorical(pivot_r_pe.index, categories=[str(i) for i in range(22)], ordered=True)
            pivot_r_pe = pivot_r_pe.sort_index()
            pivot_p_pe.index = pd.Categorical(pivot_p_pe.index, categories=[str(i) for i in range(22)], ordered=True)
            pivot_p_pe = pivot_p_pe.sort_index()
            
            pivot_r_pe.index = new_index
            pivot_p_pe.index = new_index
        elif resolution == "kmeans_subcluster_res_0.5":
            pivot_r_pe.index = pd.Categorical(pivot_r_pe.index, categories=sorted_cats, ordered=True)
            pivot_r_pe = pivot_r_pe.sort_index()
            pivot_p_pe.index = pd.Categorical(pivot_p_pe.index, categories=sorted_cats, ordered=True)
            pivot_p_pe = pivot_p_pe.sort_index()
        
        annot_matrix_pe = pd.DataFrame("", index=pivot_r_pe.index, columns=pivot_r_pe.columns)
        for r_idx in pivot_r_pe.index:
            for c_idx in pivot_r_pe.columns:
                p_val = pivot_p_pe.loc[r_idx, c_idx]
                r_val = pivot_r_pe.loc[r_idx, c_idx]
                if pd.isna(p_val) or pd.isna(r_val):
                    continue
                if p_val < 0.001:
                    annot_matrix_pe.loc[r_idx, c_idx] = "***"
                elif p_val < 0.01:
                    annot_matrix_pe.loc[r_idx, c_idx] = "**"
                elif p_val < 0.05:
                    annot_matrix_pe.loc[r_idx, c_idx] = "*"
                    
        plt.figure(figsize=(14, 11), dpi=300)
        sns.heatmap(
            pivot_r_pe,
            annot=annot_matrix_pe,
            fmt="",
            cmap="vlag",
            center=0,
            cbar_kws={'label': "Pearson Correlation Coefficient (R vs NR)"},
            linewidths=0.5,
            edgecolor='black'
        )
        plt.title(f"Pearson Clinical Response Correlation: Single-Cell Reference Deconvolution ({resolution})\n(* p < 0.05, ** p < 0.01, *** p < 0.001)", fontsize=13, fontweight='bold', y=1.02)
        plt.xlabel("Cohort", fontweight='bold')
        plt.ylabel("Inferred Cell Type / Cluster", fontweight='bold')
        plt.tight_layout()
        
        plot_path_pe = output_dir / f"self_ref_response_correlation_pearson_{resolution}.svg"
        plt.savefig(plot_path_pe, bbox_inches='tight')
        
        # Save additional copy without "pearson_" in name for reference/safety
        plot_path_pe_alt = output_dir / f"self_ref_response_correlation_pearson_{resolution}_alt.svg"
        plt.savefig(plot_path_pe_alt, bbox_inches='tight')
        plt.close()
        print(f"Saved Pearson response correlation heatmap to {plot_path_pe}")
        
        # --- HSIC Heatmap ---
        pivot_r_hsic = df_res.pivot(index='CellType', columns='Cohort', values='HSIC')
        pivot_p_hsic = df_res.pivot(index='CellType', columns='Cohort', values='HSIC_P')
        
        pivot_r_hsic = pivot_r_hsic[existing_cohorts]
        pivot_p_hsic = pivot_p_hsic[existing_cohorts]
        
        if resolution == "leiden_res_0.5":
            pivot_r_hsic.index = pd.Categorical(pivot_r_hsic.index, categories=[str(i) for i in range(22)], ordered=True)
            pivot_r_hsic = pivot_r_hsic.sort_index()
            pivot_p_hsic.index = pd.Categorical(pivot_p_hsic.index, categories=[str(i) for i in range(22)], ordered=True)
            pivot_p_hsic = pivot_p_hsic.sort_index()
            
            pivot_r_hsic.index = new_index
            pivot_p_hsic.index = new_index
        elif resolution == "kmeans_subcluster_res_0.5":
            pivot_r_hsic.index = pd.Categorical(pivot_r_hsic.index, categories=sorted_cats, ordered=True)
            pivot_r_hsic = pivot_r_hsic.sort_index()
            pivot_p_hsic.index = pd.Categorical(pivot_p_hsic.index, categories=sorted_cats, ordered=True)
            pivot_p_hsic = pivot_p_hsic.sort_index()
        
        annot_matrix_hsic = pd.DataFrame("", index=pivot_r_hsic.index, columns=pivot_r_hsic.columns)
        for r_idx in pivot_r_hsic.index:
            for c_idx in pivot_r_hsic.columns:
                p_val = pivot_p_hsic.loc[r_idx, c_idx]
                r_val = pivot_r_hsic.loc[r_idx, c_idx]
                if pd.isna(p_val) or pd.isna(r_val):
                    continue
                if p_val < 0.001:
                    annot_matrix_hsic.loc[r_idx, c_idx] = "***"
                elif p_val < 0.01:
                    annot_matrix_hsic.loc[r_idx, c_idx] = "**"
                elif p_val < 0.05:
                    annot_matrix_hsic.loc[r_idx, c_idx] = "*"
                    
        plt.figure(figsize=(14, 11), dpi=300)
        sns.heatmap(
            pivot_r_hsic,
            annot=annot_matrix_hsic,
            fmt="",
            cmap="OrRd",
            cbar_kws={'label': "Hilbert-Schmidt Independence Criterion (HSIC)"},
            linewidths=0.5,
            edgecolor='black'
        )
        plt.title(f"HSIC Clinical Response Dependence: Single-Cell Reference Deconvolution ({resolution})\n(* p < 0.05, ** p < 0.01, *** p < 0.001)", fontsize=13, fontweight='bold', y=1.02)
        plt.xlabel("Cohort", fontweight='bold')
        plt.ylabel("Inferred Cell Type / Cluster", fontweight='bold')
        plt.tight_layout()
        
        plot_path_hsic = output_dir / f"self_ref_response_correlation_hsic_{resolution}.svg"
        plt.savefig(plot_path_hsic, bbox_inches='tight')
        
        # Save additional png for easy visual artifact rendering
        plot_path_hsic_png = output_dir / f"self_ref_response_correlation_hsic_{resolution}.png"
        plt.savefig(plot_path_hsic_png, bbox_inches='tight')
        plt.close()
        print(f"Saved HSIC response dependence heatmap to {plot_path_hsic} and {plot_path_hsic_png}")


if __name__ == "__main__":
    main()
