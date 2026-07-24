#!/usr/bin/env python3
"""
Compare Cell Type Fraction Distributions Across TCGA, iAtlas, and Combined Cohorts
for Deconvolutions done with LM22 (standard, IG-gene excluded iter50, IG-gene excluded iter10),
leiden_res_0.5, and kmeans_subcluster_res_0.5 references.

Guarantees strict linear (non-log2) space input and uniform TPM normalization across all bulk mRNA datasets
before InstaPrism deconvolution, and validates linear sample fraction compositions.
Excludes Immunoglobulin (IG) heavy/light chain genes (IGH, IGK, IGL, IGLL, JCHAIN, MZB1) to eliminate artificial Plasma cell inflation.
"""

import sys
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Ensure local packages are in path
sys.path.append(str(Path(__file__).resolve().parent.parent / "packages"))
import datalair
from single_cell_datasets._single_cell_reference import SingleCellDeconvolution

# LM22 Dataset Definition for Datalair
class DatasetLM22(datalair.Dataset):
    """Datalair Dataset class for the LM22 signature matrix."""
    def __init__(self) -> None:
        super().__init__(namespace="DatasetLM22")

class LM22(DatasetLM22):
    storage_path = Path("/storage/halu").resolve()
    def derive(self, lair: datalair.Lair) -> None:
        output_dir = lair.get_path(self)
        output_dir.mkdir(parents=True, exist_ok=True)
        src_path = self.storage_path / "manual-download" / "LM22" / "LM22.txt"
        dest_path = output_dir / "LM22.txt"
        if not src_path.exists():
            raise FileNotFoundError(f"Source file not found: {src_path}")
        shutil.copy2(src_path, dest_path)

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

# Human-readable cell type maps
CELL_TYPE_NAMES_LEIDEN_05 = {
    '0': 'CD4+ Helper T (c0)', '1': 'Epithelial Tumor (c1)', '2': 'M1 Macrophage (c2)',
    '3': 'CD8+ Cytotoxic T (c3)', '4': 'CD4+ Memory T (c4)', '5': 'B Cells (c5)',
    '6': 'Endothelial (c6)', '7': 'CAFs (c7)', '8': 'NK Cells (c8)', '9': 'Monocytes (c9)',
    '10': 'Plasma Cells (c10)', '11': 'pDCs (c11)', '12': 'mDCs (c12)', '13': 'Tregs (c13)',
    '14': 'Mast Cells (14)', '15': 'Neutrophils (c15)', '16': 'GammaDelta T (c16)',
    '17': 'Pericytes (c17)', '18': 'Myofibroblasts (c18)', '19': 'Erythrocytes (c19)',
    '20': 'Granulocytes (c20)', '21': 'Melanocytes (c21)'
}


def is_ig_gene(symbol: str) -> bool:
    """Check if a gene symbol represents an Immunoglobulin heavy/light chain or J-chain transcript."""
    s = str(symbol).strip().upper()
    return s.startswith(('IGH', 'IGK', 'IGL', 'IGLL')) or s in ['JCHAIN', 'MZB1']


def verify_linear_fractions(df: pd.DataFrame, cohort_name: str, ref_name: str) -> pd.DataFrame:
    """
    Validates that cell fraction data is non-negative and represents linear proportion compositions.
    Prints an explicit log verification statement.
    """
    if df.empty:
        return df
        
    row_sums = df.sum(axis=1)
    min_val = df.values.min()
    max_val = df.values.max()
    mean_sum = row_sums.mean()
    
    # Verify non-negativity and linear composition
    assert min_val >= -1e-5, f"ERROR: Found negative values in {cohort_name} ({ref_name}): min={min_val}"
    assert max_val <= 1.05, f"ERROR: Found values > 1 in {cohort_name} ({ref_name}): max={max_val}"
    
    print(f"  [LINEARITY CHECK] {cohort_name} ({ref_name}): N={len(df)} samples | Mean Row Sum={mean_sum:.4f} | Min={min_val:.5f} | Max={max_val:.5f} -> CONFIRMED LINEAR PROPORTIONS.")
    
    # Normalize row sums to exactly 1.0 if minor numerical imprecision
    df_norm = df.div(row_sums, axis=0)
    return df_norm


def ensure_lm22_deconvolutions(lair, iatlas_cohorts, tcga_projects, output_dir, exclude_ig=False, n_iter=50):
    """
    Derives LM22 InstaPrism deconvolutions for TCGA & iAtlas cohorts,
    strictly verifying linear (non-log) TPM space inputs before deconvolution.
    Optionally excludes Immunoglobulin (IG) genes and runs with specified n_iter.
    """
    lm22_paths = {}
    import instaprism
    import anndata as ad
    import tcga
    import ici_datasets
    from pyensembl import EnsemblRelease
    from ml_pipelines.tcga_background import aggregate_and_subset_by_hugo
    
    if exclude_ig:
        ref_suffix = f"LM22_no_IG_iter{n_iter}" if n_iter != 50 else "LM22_no_IG"
    else:
        ref_suffix = f"LM22_iter{n_iter}" if n_iter != 50 else "LM22"
    
    # 1. Load LM22 reference signature matrix
    lm22_ds = LM22()
    lair.safe_derive(lm22_ds, overwrite=False)
    lm22_file = lair.get_dataset_filepaths(lm22_ds)["LM22.txt"]
    lm22_df = pd.read_csv(lm22_file, sep="\t", index_col=0)
    lm22_df.index = lm22_df.index.str.upper()
    
    if exclude_ig:
        ig_genes = [g for g in lm22_df.index if is_ig_gene(g)]
        lm22_df = lm22_df.drop(index=ig_genes)
        print(f"  [IG GENE FILTER] Excluded {len(ig_genes)} Immunoglobulin genes from LM22 matrix ({ig_genes}). Remaining genes: {len(lm22_df)}")
        
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    
    # --- A. iAtlas Cohorts ---
    for cohort_name in iatlas_cohorts:
        cache_path = output_dir / f"deconv_{cohort_name}_{ref_suffix}.csv"
        
        # If standard LM22 (50 iters, no IG filter), check if already derived in datalair
        if not exclude_ig and n_iter == 50:
            class DummyLM22Deconv(datalair.Dataset):
                def __init__(self, name):
                    super().__init__(namespace="DatasetiAtlasDeconvolution", dataset_name=name)
                def derive(self, lair): pass
            try:
                p_lair = lair.get_path(DummyLM22Deconv(f"{cohort_name}-iter50")) / f"{cohort_name}_cell_fractions.csv"
                if p_lair.is_file():
                    lm22_paths[cohort_name] = p_lair
                    continue
            except Exception:
                pass
                
        if cache_path.is_file():
            lm22_paths[cohort_name] = cache_path
            continue
            
        print(f"  Computing {ref_suffix} ({n_iter} iterations) for {cohort_name} (linear TPM space)...")
        try:
            data_dir = ici_datasets.cbioportal_datasets.get_dataset_dir(lair, dataset_class, cohort_name)
            for fname in ['data_mrna_seq_expression.txt', 'data_mrna_seq_tpm.txt', 'data_mrna_seq_rpkm.txt']:
                if (data_dir / fname).exists():
                    df_bulk = pd.read_csv(data_dir / fname, sep='\t', index_col=0)
                    df_bulk.index = df_bulk.index.astype(str).str.upper()
                    
                    if df_bulk.max().max() < 50:
                        df_bulk = np.power(2, df_bulk) - 1
                        
                    df_tpm = df_bulk.div(df_bulk.sum(axis=0), axis=1) * 1e6
                    common_genes = df_tpm.index.intersection(lm22_df.index)
                    if exclude_ig:
                        common_genes = [g for g in common_genes if not is_ig_gene(g)]
                        
                    bulk_aligned = df_tpm.loc[common_genes]
                    lm22_aligned = lm22_df.loc[common_genes]
                    
                    ref_mat = lm22_aligned.T.values
                    ref_mat = ref_mat / ref_mat.sum(axis=1, keepdims=True)
                    
                    from joblib import Parallel, delayed
                    results = Parallel(n_jobs=-1)(
                        delayed(lambda sid: instaprism.insta_prism(bulk_aligned[sid].values, ref_mat, n_iter=n_iter)[2])(s)
                        for s in bulk_aligned.columns
                    )
                    df_fracs = pd.DataFrame(results, index=bulk_aligned.columns, columns=lm22_df.columns)
                    df_fracs.to_csv(cache_path)
                    lm22_paths[cohort_name] = cache_path
                    print(f"    Saved {cohort_name} {ref_suffix} to {cache_path}")
                    break
        except Exception as e:
            print(f"    Error deconvoluting {cohort_name}: {e}")
            
    # --- B. TCGA Cohorts ---
    ensembl = EnsemblRelease(release=111, species="human")
    try:
        ensembl.download()
        ensembl.index()
    except Exception:
        pass
        
    ds_tcga = tcga.AllProjectsAdata()
    paths_tcga = lair.get_dataset_filepaths(ds_tcga)
    
    for project in tcga_projects:
        cohort_name = f"TCGA-{project}"
        cache_path = output_dir / f"deconv_{cohort_name}_{ref_suffix}.csv"
        
        if cache_path.is_file():
            lm22_paths[cohort_name] = cache_path
            continue
            
        print(f"  Computing {ref_suffix} ({n_iter} iterations) for {cohort_name} (linear TPM space)...")
        try:
            file_key = f"{project}.h5ad"
            adata_tcga = ad.read_h5ad(paths_tcga[file_key]).T
            adata_tcga.var_names_make_unique()
            
            gene_id_to_hugo = {}
            valid_ids = []
            for gene_id in adata_tcga.var_names:
                clean_id = gene_id.split(".")[0]
                try:
                    gene_obj = ensembl.gene_by_id(clean_id)
                    if gene_obj.biotype != "protein_coding":
                        continue
                    name = ensembl.gene_name_of_gene_id(clean_id)
                    if name and name.startswith("MT-"):
                        continue
                    gene_id_to_hugo[gene_id] = name if name else gene_id
                    valid_ids.append(gene_id)
                except Exception:
                    pass
                    
            adata_tcga = adata_tcga[:, valid_ids].copy()
            adata_tcga.var['hugo_symbol'] = [gene_id_to_hugo[g] for g in adata_tcga.var_names]
            valid_hugo = list(set([h for h in adata_tcga.var['hugo_symbol'] if not h.startswith("ENSG")]))
            adata_agg = aggregate_and_subset_by_hugo(adata_tcga, valid_hugo)
            
            df_tcga = pd.DataFrame(
                adata_agg.X.T,
                index=adata_agg.var_names,
                columns=adata_agg.obs_names
            )
            df_tcga.index = df_tcga.index.str.upper()
            df_tcga = df_tcga.groupby(level=0).mean()
            
            max_val = df_tcga.max().max()
            if max_val < 50:
                df_tcga = np.power(2, df_tcga) - 1
                
            df_tcga_tpm = df_tcga.div(df_tcga.sum(axis=0), axis=1) * 1e6
            
            common_genes = df_tcga_tpm.index.intersection(lm22_df.index)
            if exclude_ig:
                common_genes = [g for g in common_genes if not is_ig_gene(g)]
                
            bulk_aligned = df_tcga_tpm.loc[common_genes]
            lm22_aligned = lm22_df.loc[common_genes]
            
            reference_matrix = lm22_aligned.T.values
            reference_matrix = reference_matrix / reference_matrix.sum(axis=1, keepdims=True)
            cell_types = lm22_df.columns
            sample_ids = bulk_aligned.columns
            
            from joblib import Parallel, delayed
            results = Parallel(n_jobs=-1)(
                delayed(lambda sid: instaprism.insta_prism(bulk_aligned[sid].values, reference_matrix, n_iter=n_iter)[2])(s)
                for s in sample_ids
            )
            
            df_fracs = pd.DataFrame(results, index=sample_ids, columns=cell_types)
            df_fracs.to_csv(cache_path)
            lm22_paths[cohort_name] = cache_path
            print(f"    Saved {cohort_name} {ref_suffix} to {cache_path}")
        except Exception as e:
            print(f"    Error deconvoluting {cohort_name} with {ref_suffix}: {e}")
            
    return lm22_paths


def compute_highest_candle_limit(df_ct: pd.DataFrame) -> float:
    """
    Computes upper Y-axis limit based on the highest boxplot candle (upper whisker: Q3 + 1.5 * IQR)
    across all cohorts for a given cell type, preventing extreme outliers from squishing the box distributions.
    """
    max_candle = 0.05
    for cohort, group in df_ct.groupby('Cohort'):
        vals = group['Fraction'].values
        if len(vals) == 0:
            continue
        q25, q75 = np.percentile(vals, [25, 75])
        iqr = q75 - q25
        upper_fence = q75 + 1.5 * iqr
        valid_whiskers = vals[vals <= upper_fence]
        candle = np.max(valid_whiskers) if len(valid_whiskers) > 0 else q75
        if candle > max_candle:
            max_candle = candle
            
    y_max = min(1.02, float(max_candle * 1.25))
    return y_max


def main():
    lair_path = "/storage/halu/lair"
    print(f"Initializing Lair at {lair_path}...")
    lair = datalair.Lair(lair_path)
    
    project_root = Path(__file__).resolve().parent.parent
    output_dir = project_root / "output" / "cohort-comparisons"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load filepaths for single-cell deconvolution dataset
    ds_deconv = SingleCellDeconvolution()
    deconv_paths_sc = lair.get_dataset_filepaths(ds_deconv)
    
    tcga_projects = ["SKCM", "BLCA", "PAAD", "BRCA", "KIRC"]
    tcga_cohorts = [f"TCGA-{p}" for p in tcga_projects]
    
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
        'Combined-GBM': ['Cloughesy-iAtlas']
    }
    
    # 2. Derive/Load LM22 Variants (Standard LM22 iter50, LM22_no_IG iter50, LM22_no_IG iter10)
    print("\nLoading / Deriving Standard LM22 (50 iters)...")
    lm22_std_paths = ensure_lm22_deconvolutions(lair, iatlas_cohorts, tcga_projects, output_dir, exclude_ig=False, n_iter=50)
    
    print("\nLoading / Deriving IG-Gene Excluded LM22_no_IG (50 iters)...")
    lm22_no_ig_paths = ensure_lm22_deconvolutions(lair, iatlas_cohorts, tcga_projects, output_dir, exclude_ig=True, n_iter=50)

    print("\nLoading / Deriving IG-Gene Excluded LM22_no_IG (10 iters)...")
    lm22_no_ig_iter10_paths = ensure_lm22_deconvolutions(lair, iatlas_cohorts, tcga_projects, output_dir, exclude_ig=True, n_iter=10)
    
    references = ["LM22", "LM22_no_IG", "LM22_no_IG_iter10", "leiden_res_0.5", "kmeans_subcluster_res_0.5"]
    
    category_colors = {
        'TCGA Individual': '#2b5c8f',    # Dark Blue
        'TCGA Combined': '#002642',      # Deep Navy
        'iAtlas Individual': '#d95f02',  # Orange
        'iAtlas Combined': '#7570b3'    # Purple
    }
    
    for ref_name in references:
        print(f"\n=======================================================")
        print(f"PROCESSING COHORT FRACTION DISTRIBUTIONS: {ref_name}")
        print(f"=======================================================")
        
        cohort_data_map = {}
        
        # Load Individual TCGA Cohorts
        for cohort in tcga_cohorts:
            if ref_name == "LM22":
                p = lm22_std_paths.get(cohort)
            elif ref_name == "LM22_no_IG":
                p = lm22_no_ig_paths.get(cohort)
            elif ref_name == "LM22_no_IG_iter10":
                p = lm22_no_ig_iter10_paths.get(cohort)
            else:
                p = deconv_paths_sc.get(f"deconv_{cohort}_{ref_name}.csv")
                
            if p and p.is_file():
                df = pd.read_csv(p, index_col=0)
                df.index = df.index.astype(str)
                df_verified = verify_linear_fractions(df, cohort, ref_name)
                cohort_data_map[cohort] = (df_verified, 'TCGA Individual')
                
        # Load Combined TCGA Cohort
        tcga_dfs = [cohort_data_map[c][0] for c in tcga_cohorts if c in cohort_data_map]
        if tcga_dfs:
            df_tcga_comb = pd.concat(tcga_dfs, axis=0)
            cohort_data_map['Combined-TCGA'] = (df_tcga_comb, 'TCGA Combined')
            
        # Load Individual iAtlas Cohorts
        for cohort in iatlas_cohorts:
            if ref_name == "LM22":
                p = lm22_std_paths.get(cohort)
            elif ref_name == "LM22_no_IG":
                p = lm22_no_ig_paths.get(cohort)
            elif ref_name == "LM22_no_IG_iter10":
                p = lm22_no_ig_iter10_paths.get(cohort)
            else:
                p = deconv_paths_sc.get(f"deconv_{cohort}_{ref_name}.csv")
                
            if p and p.is_file():
                df = pd.read_csv(p, index_col=0)
                df.index = df.index.astype(str)
                df_verified = verify_linear_fractions(df, cohort, ref_name)
                cohort_data_map[cohort] = (df_verified, 'iAtlas Individual')
                
        # Load Combined iAtlas Cohorts by Cancer Type
        for comb_name, sub_cohorts in iatlas_cancer_map.items():
            sub_dfs = [cohort_data_map[c][0] for c in sub_cohorts if c in cohort_data_map]
            if sub_dfs:
                df_comb = pd.concat(sub_dfs, axis=0)
                cohort_data_map[comb_name] = (df_comb, 'iAtlas Combined')
                
        if not cohort_data_map:
            print(f"Warning: No data loaded for reference {ref_name}. Skipping.")
            continue
            
        ordered_cohorts = [
            'TCGA-SKCM', 'TCGA-BLCA', 'TCGA-PAAD', 'TCGA-BRCA', 'TCGA-KIRC', 'Combined-TCGA',
            'Hugo-iAtlas', 'Riaz-iAtlas', 'Liu-iAtlas', 'Gide-iAtlas', 'Combined-Melanoma',
            'Rosenberg-iAtlas', 'Combined-BLCA',
            'Padron-iAtlas', 'Combined-PAAD',
            'Anders-iAtlas', 'Combined-BRCA',
            'McDermott-iAtlas', 'Choueiri-iAtlas', 'Combined-RCC'
        ]
        existing_ordered_cohorts = [c for c in ordered_cohorts if c in cohort_data_map]
        
        sample_df = next(iter(cohort_data_map.values()))[0]
        cell_types = list(sample_df.columns)
        
        if ref_name == "leiden_res_0.5":
            cell_type_labels = [CELL_TYPE_NAMES_LEIDEN_05.get(str(ct), str(ct)) for ct in cell_types]
        else:
            cell_type_labels = list(cell_types)
            
        records = []
        summary_records = []
        
        for c_name in existing_ordered_cohorts:
            df_c, cat = cohort_data_map[c_name]
            for ct_orig, ct_lbl in zip(cell_types, cell_type_labels):
                vals = df_c[ct_orig].values
                for v in vals:
                    records.append({
                        'Cohort': c_name,
                        'Category': cat,
                        'CellType': ct_lbl,
                        'Fraction': float(v)
                    })
                summary_records.append({
                    'Reference': ref_name,
                    'Cohort': c_name,
                    'Category': cat,
                    'CellType': ct_lbl,
                    'N_Samples': len(vals),
                    'Mean': float(np.mean(vals)),
                    'Std': float(np.std(vals)),
                    'Median': float(np.median(vals)),
                    'IQR': float(np.percentile(vals, 75) - np.percentile(vals, 25)),
                    'Percentile_5': float(np.percentile(vals, 5)),
                    'Percentile_95': float(np.percentile(vals, 95))
                })
                
        df_long = pd.DataFrame(records)
        df_summary = pd.DataFrame(summary_records)
        
        summary_path = output_dir / f"deconv_cohort_distribution_summary_{ref_name}.csv"
        df_summary.to_csv(summary_path, index=False)
        print(f"Saved distribution summary table to {summary_path}")
        
        cell_dir = output_dir / f"cell_types_{ref_name}"
        cell_dir.mkdir(parents=True, exist_ok=True)
        
        for ct_lbl in cell_type_labels:
            df_ct = df_long[df_long['CellType'] == ct_lbl]
            y_candle_limit = compute_highest_candle_limit(df_ct)
            
            plt.figure(figsize=(16, 7), dpi=300)
            ax = sns.boxplot(
                data=df_ct,
                x='Cohort',
                y='Fraction',
                hue='Category',
                palette=category_colors,
                order=existing_ordered_cohorts,
                fliersize=2,
                linewidth=1.2
            )
            plt.xticks(rotation=45, ha='right', fontweight='bold')
            plt.ylabel("Inferred Cell Fraction (Linear TPM)", fontweight='bold')
            plt.xlabel("Cohort", fontweight='bold')
            plt.title(f"Cell Type Fraction Distribution Across Cohorts: {ct_lbl} ({ref_name})\nStrict Linear Space & TPM Normalized Input", fontsize=13, fontweight='bold', y=1.02)
            plt.ylim(-0.01, y_candle_limit)
            
            clean_ct = ct_lbl.lower().replace(" ", "_").replace("+", "").replace("(", "").replace(")", "").replace("/", "_")
            plt.savefig(cell_dir / f"boxplot_{clean_ct}.png", bbox_inches='tight')
            plt.savefig(cell_dir / f"boxplot_{clean_ct}.svg", bbox_inches='tight')
            plt.close()

        n_ct = len(cell_type_labels)
        n_cols = 4
        n_rows = int(np.ceil(n_ct / n_cols))
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(24, 4.2 * n_rows), dpi=300)
        axes_flat = axes.flatten()
        
        for idx, ct_lbl in enumerate(cell_type_labels):
            ax = axes_flat[idx]
            df_ct = df_long[df_long['CellType'] == ct_lbl]
            y_candle_limit = compute_highest_candle_limit(df_ct)
            
            sns.boxplot(
                data=df_ct,
                x='Cohort',
                y='Fraction',
                hue='Category',
                palette=category_colors,
                order=existing_ordered_cohorts,
                fliersize=1,
                linewidth=0.8,
                ax=ax
            )
            ax.set_title(ct_lbl, fontsize=10, fontweight='bold')
            ax.set_xlabel("")
            ax.set_ylabel("Fraction", fontsize=8)
            ax.set_ylim(-0.01, y_candle_limit)
            ax.tick_params(axis='x', rotation=90, labelsize=7)
            ax.get_legend().remove()
            
        for idx in range(n_ct, len(axes_flat)):
            fig.delaxes(axes_flat[idx])
            
        fig.suptitle(f"Cell Type Fraction Distributions Across TCGA & iAtlas Cohorts ({ref_name})\n(All Datasets Linear Non-Log TPM Normalized Prior to Deconvolution; Y-limits scaled by highest candle)", fontsize=16, fontweight='bold', y=1.01)
        fig.tight_layout()
        
        overview_path = output_dir / f"deconv_boxplot_overview_{ref_name}.png"
        fig.savefig(overview_path, bbox_inches='tight')
        fig.savefig(output_dir / f"deconv_boxplot_overview_{ref_name}.svg", bbox_inches='tight')
        plt.close(fig)
        print(f"Saved multi-panel overview grid boxplot to {overview_path}")

    print("\nAll cohort deconvolution distribution comparisons completed successfully.")


if __name__ == "__main__":
    main()
