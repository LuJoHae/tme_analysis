import os
import shutil
import sys
from pathlib import Path
from typing import Callable
import numpy as np
import pandas as pd
import datalair
import instaprism
import ici_datasets
from pyensembl import EnsemblRelease
from joblib import Parallel, delayed
from scipy.stats import spearmanr, mannwhitneyu
from gene_utils import calculate_maf_tmb

# --- LM22 Dataset Definition ---

class DatasetLM22(datalair.Dataset):
    """Datalair Dataset class for the LM22 signature matrix."""

    def __init__(self) -> None:
        """Initialize the LM22 dataset namespace."""
        super().__init__(namespace="DatasetLM22")

class LM22(DatasetLM22):
    storage_path = Path("/storage/halu").resolve()

    def derive(self, lair: datalair.Lair) -> None:
        """Derive the LM22 dataset by copying the file from manual downloads.

        Parameters:
        lair: Lair instance handling the dataset repository.
        """
        output_dir = lair.get_path(self)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        src_path = self.storage_path / "manual-download" / "LM22" / "LM22.txt"
        dest_path = output_dir / "LM22.txt"
        
        if not src_path.exists():
            raise FileNotFoundError(f"Source file not found: {src_path}")
            
        shutil.copy2(src_path, dest_path)
        print(f"Successfully copied LM22.txt from {src_path} to {dest_path}")

# --- mRNA Loader & Helper Functions ---

def transform_expression_to_tpm(counts_df: pd.DataFrame, release: int = 110) -> pd.DataFrame:
    """Computes TPM from raw counts using pyensembl for local gene metrics lookup.

    Parameters:
    counts_df: DataFrame of raw counts (rows = Hugo symbols, cols = samples).
    release: Ensembl release version number.
    """
    data_source = EnsemblRelease(release)

    def extract_kb_span(symbol: str) -> float | None:
        """Maps a Hugo symbol to its genomic span in kilobases."""
        try:
            genes = data_source.genes_by_name(symbol)
            return (genes[0].end - genes[0].start + 1) / 1e3 if genes else None
        except ValueError:
            return None

    gene_lengths = pd.Series({
        symbol: length
        for symbol in counts_df.index
        if (length := extract_kb_span(symbol)) is not None
    })

    common_genes = counts_df.index.intersection(gene_lengths.index)
    counts = counts_df.loc[common_genes]
    lengths = gene_lengths.loc[common_genes]

    rpk = counts.div(lengths, axis=0)
    return rpk.div(rpk.sum(axis=0), axis=1) * 1e6

def transform_rpkm_to_tpm(rpkm_df: pd.DataFrame) -> pd.DataFrame:
    """Converts an RPKM/FPKM normalized DataFrame to TPM."""
    return rpkm_df.div(rpkm_df.sum(axis=0), axis=1) * 1e6

def transform_tpm_to_tpm(df: pd.DataFrame) -> pd.DataFrame:
    """Identity transformation for TPM normalized data."""
    return df

def load_and_transform_data_mrna(data_dir: str | Path) -> pd.DataFrame:
    """Loads the bulk RNA expression data and converts it to TPM."""
    p_dir = Path(data_dir)
    dispatch_map: dict[str, tuple[Callable[[pd.DataFrame], pd.DataFrame], str]] = {
        "data_mrna_seq_expression.txt": (transform_expression_to_tpm, "none"),
        "data_mrna_seq_tpm.txt": (transform_tpm_to_tpm, "tpm"),
        "data_mrna_seq_rpkm.txt": (transform_rpkm_to_tpm, "rpkm"),
    }

    existing_files = tuple(filter(lambda f: (p_dir / f).is_file(), dispatch_map.keys()))
    if len(existing_files) != 1:
        raise ValueError(
            f"Invariant violation: Expected exactly 1 file, but found {len(existing_files)} in {p_dir}."
        )
    target = existing_files[0]
    transform_fn, normalization = dispatch_map[target]

    data_mrna = transform_fn(pd.read_csv(p_dir / target, sep='\t', index_col=0))
    data_mrna.attrs["original_normalization"] = normalization
    if data_mrna.isna().any(axis=None):
        raise ValueError("Data integrity violation: DataFrame contains NA values.")
    return data_mrna

def get_dataset_dir(lair, dataset_class, ds_name):
    """Retrieves the directory of a CBioPortal dataset inside datalair."""
    ds = dataset_class(name=ds_name)
    filepaths = lair.get_dataset_filepaths(ds)
    unpacked_file_key = next(f for f in filepaths.keys() if not f.endswith('.tar.gz'))
    filepath = filepaths[unpacked_file_key]
    return filepath / filepath.name

# --- Deconvolution Dataset Classes ---

class DatasetiAtlasDeconvolution(datalair.Dataset):
    """Datalair Dataset class for deconvolution of iAtlas datasets."""

    def __init__(self, name: str, n_iter: int) -> None:
        """Initialize the dataset deconvolution.

        Parameters:
        name: Name of the iAtlas dataset (e.g. 'Gide-iAtlas').
        n_iter: Number of deconvolution iterations.
        """
        super().__init__(namespace="DatasetiAtlasDeconvolution", dataset_name=f"{name}-iter{n_iter}")
        self._cohort_name = name

class DeconvolutioniAtlas(DatasetiAtlasDeconvolution):
    def __init__(self, name: str, n_iter: int) -> None:
        super().__init__(name=name, n_iter=n_iter)
        self.n_iter = n_iter

    def derive(self, lair: datalair.Lair) -> None:
        """Runs instaprism deconvolution using LM22 as reference."""
        output_dir = lair.get_path(self)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 1. Derive LM22 reference
        lm22_ds = LM22()
        lair.safe_derive(lm22_ds)
        lm22_path = lair.get_dataset_filepaths(lm22_ds)["LM22.txt"]
        lm22_df = pd.read_csv(lm22_path, sep="\t", index_col=0)
        
        # 2. Safe-derive target CBioPortal dataset
        dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
        cb_ds = dataset_class(name=self._cohort_name)
        lair.safe_derive(cb_ds)
        
        # 3. Load bulk mRNA
        data_dir = get_dataset_dir(lair, dataset_class, self._cohort_name)
        bulk_df = load_and_transform_data_mrna(data_dir)
        
        # 4. Match genes (Hugo symbols)
        common_genes = bulk_df.index.intersection(lm22_df.index)
        if len(common_genes) == 0:
            raise ValueError(f"No common genes found between bulk expression and LM22 in {self._cohort_name}")
            
        print(f"[{self._cohort_name}-iter{self.n_iter}] Matching genes: {len(common_genes)} out of {len(lm22_df)} in LM22 reference.")
        
        # Align rows based on common genes
        bulk_aligned = bulk_df.loc[common_genes]
        lm22_aligned = lm22_df.loc[common_genes]
        
        # 5. Perform deconvolution
        # instaprism expects reference of shape (S, G), where S is cell types, G is genes
        reference_matrix = lm22_aligned.T.values
        cell_types = lm22_df.columns
        sample_ids = bulk_df.columns
        
        def deconvolute_sample(sample_id):
            bulk_sample = bulk_aligned[sample_id].values
            # insta_prism returns (probability_matrix, cell_state_gene_expression, cell_fractions, intermediate_cell_fractions)
            _, _, cell_fracs, _ = instaprism.insta_prism(
                bulk=bulk_sample,
                reference=reference_matrix,
                n_iter=self.n_iter
            )
            return cell_fracs

        print(f"[{self._cohort_name}-iter{self.n_iter}] Running deconvolution on {len(sample_ids)} samples...")
        results = Parallel(n_jobs=-1)(
            delayed(deconvolute_sample)(sid) for sid in sample_ids
        )
        
        # Create Output DataFrame and save to CSV
        deconv_df = pd.DataFrame(results, index=sample_ids, columns=cell_types)
        output_file = output_dir / f"{self._cohort_name}_cell_fractions.csv"
        deconv_df.to_csv(output_file)
        print(f"[{self._cohort_name}-iter{self.n_iter}] Successfully saved cell fractions to {output_file}")

# --- Immunotherapy Response Correlation Functions ---

def load_data_clinical(data_dir: Path) -> pd.DataFrame:
    """Loads response and TMB information from sample clinical and mutation files."""
    RESPONSE_COL_CANDIDATES = [
        "RESPONSE", "Response", "response",
        "BEST_RESPONSE", "Best_Response", "best_response",
        "RECIST", "recist",
    ]
    RESPONSE_VALUE_MAP = {
        "Complete Response": "R",
        "Partial Response": "R",
        "Stable Disease": "NR",
        "Progressive Disease": "NR",
        "CR": "R",
        "PR": "R",
        "SD": "NR",
        "PD": "NR",
    }

    def _resolve_column(df, candidates):
        normalized = {c.strip().lower(): c for c in df.columns}
        for cand in candidates:
            key = cand.strip().lower()
            if key in normalized:
                return normalized[key]
        raise KeyError(f"None of {candidates} found in clinical columns: {list(df.columns)}")

    filepath_clinical_sample = data_dir / "data_clinical_sample.txt"
    data_clinical_sample = pd.read_csv(filepath_clinical_sample, sep="\t", index_col=0, skiprows=4)
    data_clinical_sample = data_clinical_sample.set_index("SAMPLE_ID")
    
    response_col = _resolve_column(data_clinical_sample, RESPONSE_COL_CANDIDATES)
    data_clinical_sample["response"] = data_clinical_sample[response_col].map(RESPONSE_VALUE_MAP)
    df_out = data_clinical_sample[["response"]].copy()

    filepath_mutations = data_dir / "data_mutations.txt"
    if filepath_mutations.exists():
        data_mutations = pd.read_csv(filepath_mutations, sep="\t", low_memory=False)
        tmb_results = calculate_maf_tmb(data_mutations)
        tmb_results = tmb_results.set_index("Tumor_Sample_Barcode")
        df_out = pd.concat([tmb_results[["TMB_Score"]], df_out], join="inner", axis=1)
    else:
        df_out["TMB_Score"] = np.nan
        
    return df_out

def fdr_bh(pvals):
    """Computes Benjamini-Hochberg False Discovery Rate q-values."""
    pvals = np.asarray(pvals, dtype=np.float64)
    desc_idx = np.argsort(pvals)[::-1]
    idx = np.argsort(desc_idx)
    pvals = pvals[desc_idx]
    qvals = np.zeros_like(pvals)
    n = len(pvals)
    min_q = 1.0
    for i, p in enumerate(pvals):
        q = p * n / (n - i)
        min_q = min(min_q, q)
        qvals[i] = min_q
    return qvals[idx]

def compute_cohort_correlations(deconv_valid: pd.DataFrame, y_encoded: pd.Series, cohort_name: str) -> list[dict]:
    """Computes Spearman correlation and MWU p-values for a cohort subset."""
    cohort_records = []
    for cell_type in deconv_valid.columns:
        x = deconv_valid[cell_type]
        
        # Spearman correlation
        corr, p_corr = spearmanr(x, y_encoded)
        
        # Mann-Whitney U test
        x_r = x[y_encoded == 1]
        x_nr = x[y_encoded == 0]
        if len(x_r) > 0 and len(x_nr) > 0:
            stat, p_mwu = mannwhitneyu(x_r, x_nr, alternative="two-sided")
        else:
            corr, p_corr, p_mwu = np.nan, np.nan, np.nan
        
        cohort_records.append({
            "cohort": cohort_name,
            "cell_type": cell_type,
            "spearman_corr": corr if not np.isnan(corr) else 0.0,
            "spearman_p": p_corr,
            "mwu_p": p_mwu,
            "mean_R": x_r.mean() if len(x_r) > 0 else np.nan,
            "mean_NR": x_nr.mean() if len(x_nr) > 0 else np.nan,
        })
        
    # Adjust p-values using FDR Benjamini-Hochberg per cohort
    mwu_ps = [r["mwu_p"] for r in cohort_records]
    if not np.all(np.isnan(mwu_ps)):
        q_vals = fdr_bh([p if not np.isnan(p) else 1.0 for p in mwu_ps])
        for r, q in zip(cohort_records, q_vals):
            r["mwu_fdr_q"] = q if not np.isnan(r["mwu_p"]) else np.nan
    else:
        for r in cohort_records:
            r["mwu_fdr_q"] = np.nan
            
    return cohort_records

def plot_heatmap(df_results: pd.DataFrame, title: str, output_path: Path):
    """Generates and saves the correlation heatmap from results."""
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    corr_pivot = df_results.pivot(index="cell_type", columns="cohort", values="spearman_corr")
    p_pivot = df_results.pivot(index="cell_type", columns="cohort", values="mwu_p")
    
    corr_pivot = corr_pivot.fillna(0.0)
    p_pivot = p_pivot.fillna(1.0)
    
    # Significance annotation matrix
    annot_matrix = pd.DataFrame("", index=corr_pivot.index, columns=corr_pivot.columns)
    for r in corr_pivot.index:
        for c in corr_pivot.columns:
            pval = p_pivot.loc[r, c]
            if pval < 0.001:
                annot_matrix.loc[r, c] = "***"
            elif pval < 0.01:
                annot_matrix.loc[r, c] = "**"
            elif pval < 0.05:
                annot_matrix.loc[r, c] = "*"
                
    fig, ax = plt.subplots(figsize=(14, 11))
    sns.heatmap(
        corr_pivot,
        cmap="RdBu_r",
        center=0.0,
        annot=annot_matrix.values,
        fmt="",
        cbar_kws={'label': "Spearman Correlation Coefficient (R vs NR)"},
        linewidths=0.5,
        linecolor="gray",
        ax=ax
    )
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel("iAtlas Cohorts", fontsize=12)
    ax.set_ylabel("LM22 Cell Types", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"Successfully saved heatmap to {output_path}")

def run_correlation_analysis(lair: datalair.Lair, output_dir: Path, n_iter: int):
    """Computes correlation between cell type fractions and clinical response (overall & split by TMB)."""
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    iatlas_dataset_names = [s for s in dataset_class.get_dataset_name() if "iAtlas" in s]
    
    records_all = []
    records_tmb_low = []
    records_tmb_high = []
    
    for ds_name in iatlas_dataset_names:
        try:
            deconv_ds = DeconvolutioniAtlas(name=ds_name, n_iter=n_iter)
            deconv_path = lair.get_dataset_filepaths(deconv_ds)[f"{ds_name}_cell_fractions.csv"]
            deconv_df = pd.read_csv(deconv_path, index_col=0)
            
            data_dir = get_dataset_dir(lair, dataset_class, ds_name)
            clinical_df = load_data_clinical(data_dir)
        except Exception as e:
            print(f"Skipping {ds_name} due to missing data: {e}")
            continue
            
        common_samples = deconv_df.index.intersection(clinical_df.index)
        clinical_valid = clinical_df.loc[common_samples].dropna(subset=["response"])
        deconv_valid = deconv_df.loc[clinical_valid.index]
        
        if len(clinical_valid) < 5:
            print(f"Skipping {ds_name} - too few samples with response: {len(clinical_valid)}")
            continue
            
        y_encoded = (clinical_valid["response"] == "R").astype(int)
        if len(np.unique(y_encoded)) < 2:
            print(f"Skipping {ds_name} - single class present in response")
            continue
            
        # 1. Overall correlations
        records_all.extend(compute_cohort_correlations(deconv_valid, y_encoded, ds_name))
        
        # 2. TMB-split correlations
        if "TMB_Score" in clinical_valid.columns and not clinical_valid["TMB_Score"].isna().all():
            tmb_vals = clinical_valid["TMB_Score"].dropna()
            if len(tmb_vals) >= 4:
                median_tmb = tmb_vals.median()
                
                # Split indices
                low_idx = clinical_valid[clinical_valid["TMB_Score"] <= median_tmb].index
                high_idx = clinical_valid[clinical_valid["TMB_Score"] > median_tmb].index
                
                # Run Low TMB split
                if len(low_idx) >= 5:
                    y_low = (clinical_valid.loc[low_idx, "response"] == "R").astype(int)
                    if len(np.unique(y_low)) >= 2:
                        records_tmb_low.extend(compute_cohort_correlations(deconv_valid.loc[low_idx], y_low, ds_name))
                    else:
                        print(f"Skipping TMB-Low for {ds_name} - single response class present")
                else:
                    print(f"Skipping TMB-Low for {ds_name} - too few samples ({len(low_idx)})")
                    
                # Run High TMB split
                if len(high_idx) >= 5:
                    y_high = (clinical_valid.loc[high_idx, "response"] == "R").astype(int)
                    if len(np.unique(y_high)) >= 2:
                        records_tmb_high.extend(compute_cohort_correlations(deconv_valid.loc[high_idx], y_high, ds_name))
                    else:
                        print(f"Skipping TMB-High for {ds_name} - single response class present")
                else:
                    print(f"Skipping TMB-High for {ds_name} - too few samples ({len(high_idx)})")
            else:
                print(f"Skipping TMB-split for {ds_name} - too few TMB non-NaN values ({len(tmb_vals)})")
        else:
            print(f"Skipping TMB-split for {ds_name} - no TMB_Score column or all NaN")
            
    stats_dir = output_dir / "iAtlas-deconvolution"
    stats_dir.mkdir(parents=True, exist_ok=True)
    
    suffix = f"_iter{n_iter}" if n_iter != 1000 else ""
    
    # Process & Plot Overall
    if records_all:
        df_all = pd.DataFrame(records_all)
        df_all.to_csv(stats_dir / f"correlation_summary{suffix}.csv", index=False)
        plot_heatmap(
            df_all,
            f"Correlation of LM22 Inferred Cell Fractions & Immunotherapy Response (All Samples, {n_iter} Iters)\n(* p < 0.05, ** p < 0.01, *** p < 0.001)",
            stats_dir / f"cell_type_response_correlation{suffix}.svg"
        )
        
    # Process & Plot TMB Low
    if records_tmb_low:
        df_low = pd.DataFrame(records_tmb_low)
        df_low.to_csv(stats_dir / f"correlation_summary_tmb_low{suffix}.csv", index=False)
        plot_heatmap(
            df_low,
            f"Correlation of LM22 Inferred Cell Fractions & Immunotherapy Response (TMB-Low Samples, {n_iter} Iters)\n(* p < 0.05, ** p < 0.01, *** p < 0.001)",
            stats_dir / f"cell_type_response_correlation_tmb_low{suffix}.svg"
        )
        
    # Process & Plot TMB High
    if records_tmb_high:
        df_high = pd.DataFrame(records_tmb_high)
        df_high.to_csv(stats_dir / f"correlation_summary_tmb_high{suffix}.csv", index=False)
        plot_heatmap(
            df_high,
            f"Correlation of LM22 Inferred Cell Fractions & Immunotherapy Response (TMB-High Samples, {n_iter} Iters)\n(* p < 0.05, ** p < 0.01, *** p < 0.001)",
            stats_dir / f"cell_type_response_correlation_tmb_high{suffix}.svg"
        )

# --- Main Driver ---

def main():
    lair = datalair.Lair("/storage/halu/lair")
    lair.assert_ok_satus()
    
    # 1. Derive LM22 reference first
    print("Safe-deriving LM22 dataset...")
    lm22_ds = LM22()
    lair.safe_derive(lm22_ds)
    
    # 2. Loop over iterations
    iterations = [10, 20, 30, 40, 50, 1000]
    dataset_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
    iatlas_dataset_names = [s for s in dataset_class.get_dataset_name() if "iAtlas" in s]
    
    for n_iter in iterations:
        print(f"\n=========================================")
        print(f"RUNNING DECONVOLUTION: {n_iter} ITERATIONS")
        print(f"=========================================")
        for ds_name in iatlas_dataset_names:
            print(f"Deriving/loading deconvolution for {ds_name} ({n_iter} iterations)...")
            deconv_ds = DeconvolutioniAtlas(name=ds_name, n_iter=n_iter)
            lair.safe_derive(deconv_ds, overwrite=False) # Reuse cached results if already derived
            
        print(f"\n=========================================")
        print(f"Running response correlation analysis for {n_iter} iterations...")
        print(f"=========================================")
        output_dir = Path(".") / "output"
        run_correlation_analysis(lair, output_dir, n_iter=n_iter)
        
    print("\nAll deconvolution & correlation analyses completed successfully.")

if __name__ == "__main__":
    main()
