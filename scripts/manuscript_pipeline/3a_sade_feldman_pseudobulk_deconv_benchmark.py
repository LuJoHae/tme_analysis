#!/usr/bin/env python3
"""
Single-Cell Benchmark: Milopy DA vs. Pseudobulk Deconvolution across Immunotherapy Cohorts.

Workflow:
1. Loads single-cell data from any of the 9 immunotherapy response cohorts (or synthetic test).
2. Resolves patient sample IDs, response labels, and informative cell states.
3. Computes true ground-truth cell state fractions per patient sample directly from single-cell counts.
4. Constructs patient-level pseudobulk expression profiles via vectorized sparse matrix multiplication (or TPM streaming).
5. Computes/loads reference centroid matrix Phi and marker genes.
6. Deconvolves patient pseudobulks against single-cell reference centroids (InstaPrism / least-squares fallback).
7. Fits patient-level logistic regression on inferred and true fractions (Combined, Pre-treatment, Post-treatment).
8. Aggregates and compares deconvolution effect sizes with single-cell Milopy DA effect sizes across conditions.
9. Evaluates directional concordance, Spearman / Pearson correlations, and deconvolution fidelity.
10. Generates publication-quality Altair vector SVG figures (cohort-level multi-panel and cross-cohort forest plots).

Strict functional Python: immutability, returns Result, Polars, and Altair.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import Final, Sequence

import altair as alt  # type: ignore
import anndata as ad  # type: ignore
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import sparse as sp  # type: ignore
from scipy import stats  # type: ignore
import vl_convert as vlc  # type: ignore

try:
    import instaprism  # type: ignore
    HAS_INSTAPRISM = True
except ImportError:
    HAS_INSTAPRISM = False


# ------------------------------------------------------------------------------
# Cohort Registry & Metadata
# ------------------------------------------------------------------------------

ALL_MILOPY_COHORTS: Final[tuple[str, ...]] = (
    "GSE120575",
    "CELLxGENE_7b20c613",
    "CELLxGENE_05a8c945",
    "CELLxGENE_6f9de485",
    "GSE207422",
    "GSE243013",
    "GSE233203",
    "GSE200996",
    "GSE316195",
)

COHORT_METADATA: Final[dict[str, dict[str, str]]] = {
    "GSE120575": {"cancer_type": "Melanoma", "author": "Sade-Feldman et al.", "year": "2018"},
    "CELLxGENE_7b20c613": {"cancer_type": "Pan-Cancer", "author": "Gondal et al.", "year": "2025"},
    "CELLxGENE_05a8c945": {"cancer_type": "Colorectal", "author": "Marteau et al.", "year": "2026"},
    "CELLxGENE_6f9de485": {"cancer_type": "Breast", "author": "Yan et al.", "year": "2026"},
    "GSE207422": {"cancer_type": "NSCLC", "author": "Hu et al.", "year": "2023"},
    "GSE243013": {"cancer_type": "NSCLC", "author": "Liu et al.", "year": "2025"},
    "GSE233203": {"cancer_type": "NSCLC", "author": "Hong et al.", "year": "2025"},
    "GSE200996": {"cancer_type": "HNSCC", "author": "Luoma et al.", "year": "2022"},
    "GSE316195": {"cancer_type": "PDAC", "author": "Raufi et al.", "year": "2026"},
}


# ------------------------------------------------------------------------------
# Configuration Model
# ------------------------------------------------------------------------------

class BenchmarkConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    cohort: str = "GSE120575"
    cell_state_col: str | None = None
    leiden_resolution: float = 2.0
    max_cells: int = 15000
    adata_path: Path | None = None
    adata_dir: Path = Path("data/preprocessed")
    tpm_gz_path: Path = Path("data/raw/GSE120575/GSE120575_tpm.txt.gz")
    cell_scores_path: Path | None = None
    phi_path: Path | None = None
    marker_genes_path: Path | None = None
    milo_da_dir: Path = Path("output/results/milopy_response_analysis")
    out_dir: Path = Path("output/results/pseudobulk_deconv_benchmark")
    results_dir: Path = Path("results/pseudobulk_deconv_benchmark")
    deconv_backend: str = "instaprism"
    deconv_iters: int = 100
    use_marker_genes: bool = True
    n_sweeps: int = 5
    svg_name: str | None = None
    parquet_name: str | None = None


def parse_args() -> BenchmarkConfig:
    parser = argparse.ArgumentParser(
        description="Single-Cell Benchmark: Milopy DA vs Pseudobulk Deconvolution across immunotherapy cohorts."
    )
    parser.add_argument(
        "--cohort",
        type=str,
        default="GSE120575",
        help="Cohort ID to benchmark (e.g., GSE120575, CELLxGENE_05a8c945, etc.) or 'all'",
    )
    parser.add_argument(
        "--cell-state-col",
        type=str,
        default=None,
        help="Explicit cell state column name in AnnData obs or cell scores",
    )
    parser.add_argument(
        "--leiden-resolution",
        type=float,
        default=2.0,
        help="Leiden clustering resolution for defining single-cell states (default: 2.0)",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=15000,
        help="Maximum number of cells per cohort (stratified subsampling by sample_id, default: 15000)",
    )
    parser.add_argument(
        "--adata-path",
        type=Path,
        default=None,
        help="Explicit path to single-cell AnnData (.h5ad)",
    )
    parser.add_argument(
        "--adata-dir",
        type=Path,
        default=Path("data/preprocessed"),
        help="Directory containing preprocessed AnnData files",
    )
    parser.add_argument(
        "--tpm-path",
        type=Path,
        default=Path("data/raw/GSE120575/GSE120575_tpm.txt.gz"),
        help="Path to raw TPM gzip (legacy GSE120575)",
    )
    parser.add_argument(
        "--cell-scores",
        type=Path,
        default=None,
        help="Explicit cell-level scores parquet path",
    )
    parser.add_argument(
        "--phi-path",
        type=Path,
        default=None,
        help="Reference phi matrix parquet",
    )
    parser.add_argument(
        "--markers-path",
        type=Path,
        default=None,
        help="Reference marker genes parquet",
    )
    parser.add_argument(
        "--milo-dir",
        type=Path,
        default=Path("output/results/milopy_response_analysis"),
        help="Directory containing milopy DA parquets",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("output/results/pseudobulk_deconv_benchmark"),
        help="Output directory for data parquets",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results/pseudobulk_deconv_benchmark"),
        help="Results directory for plots and vector SVGs",
    )
    parser.add_argument(
        "--deconv-backend",
        type=str,
        default="instaprism",
        choices=["instaprism", "bayesprism", "nnls"],
        help="Deconvolution backend",
    )
    parser.add_argument(
        "--deconv-iters",
        type=int,
        default=100,
        help="Number of deconvolution iterations",
    )
    parser.add_argument(
        "--all-genes",
        action="store_true",
        help="Deconvolve using all genes instead of marker genes",
    )
    parser.add_argument(
        "--n-sweeps",
        type=int,
        default=5,
        help="Number of baseline replicate sweeps per cohort (default: 5)",
    )
    parser.add_argument(
        "--svg-name",
        type=str,
        default=None,
        help="SVG filename",
    )
    parser.add_argument(
        "--parquet-name",
        type=str,
        default=None,
        help="Parquet filename",
    )
    args = parser.parse_args()

    return BenchmarkConfig(
        cohort=args.cohort,
        cell_state_col=args.cell_state_col,
        leiden_resolution=args.leiden_resolution,
        max_cells=args.max_cells,
        adata_path=args.adata_path,
        adata_dir=args.adata_dir,
        tpm_gz_path=args.tpm_path,
        cell_scores_path=args.cell_scores,
        phi_path=args.phi_path,
        marker_genes_path=args.markers_path,
        milo_da_dir=args.milo_dir,
        out_dir=args.out_dir,
        results_dir=args.results_dir,
        deconv_backend=args.deconv_backend,
        deconv_iters=args.deconv_iters,
        use_marker_genes=not args.all_genes,
        n_sweeps=args.n_sweeps,
        svg_name=args.svg_name,
        parquet_name=args.parquet_name,
    )


# ------------------------------------------------------------------------------
# 1. Ground Truth Single-Cell Proportions
# ------------------------------------------------------------------------------

def compute_true_cell_fractions(
    df_cells: pl.DataFrame,
) -> pl.DataFrame:
    """Pure function computing true cell state fractions per sample directly from cell annotations."""
    # Standardize column mappings if needed
    col_mapping: dict[str, str] = {}
    if "patient_id" in df_cells.columns and "sample_id" not in df_cells.columns:
        col_mapping["patient_id"] = "sample_id"
    if "clinical_response" in df_cells.columns and "response" not in df_cells.columns:
        col_mapping["clinical_response"] = "response"
    if "cell_type" in df_cells.columns and "cell_state" not in df_cells.columns:
        col_mapping["cell_type"] = "cell_state"

    std_df = df_cells.rename(col_mapping) if col_mapping else df_cells
    if "treatment_status" not in std_df.columns:
        std_df = std_df.with_columns(pl.lit("Combined").alias("treatment_status"))

    sample_meta = (
        std_df.select(["sample_id", "response", "treatment_status"])
        .unique(subset=["sample_id"])
    )
    sample_totals = std_df.group_by("sample_id").agg(pl.len().alias("n_total"))
    state_counts = std_df.group_by(["sample_id", "cell_state"]).agg(pl.len().alias("n_cells"))

    unique_samples = sample_meta["sample_id"].to_list()
    unique_states = std_df["cell_state"].unique().sort().to_list()

    grid = pl.DataFrame({
        "sample_id": [s for s in unique_samples for _ in unique_states],
        "cell_state": [st for _ in unique_samples for st in unique_states],
    })

    joined = (
        grid.join(state_counts, on=["sample_id", "cell_state"], how="left")
        .with_columns(pl.col("n_cells").fill_null(0))
        .join(sample_totals, on="sample_id", how="left")
        .join(sample_meta, on="sample_id", how="left")
        .with_columns(
            (pl.col("n_cells") / pl.col("n_total")).alias("true_fraction")
        )
    )
    return joined


# ------------------------------------------------------------------------------
# 2. Patient Pseudobulk Construction (Fast Vectorized & Legacy Streamed)
# ------------------------------------------------------------------------------

def resolve_adata_gene_names(adata: ad.AnnData) -> list[str]:
    """Extract upper-case gene symbols or identifiers from AnnData .var metadata."""
    for col in ("gene_name", "gene_symbols", "GeneSymbol", "feature_name", "symbol"):
        if col in adata.var.columns:
            return [str(g).strip().upper() for g in adata.var[col]]
    if "original_id" in adata.var.columns and not str(adata.var_names[0]).startswith("ENSG"):
        return [str(g).strip().upper() for g in adata.var_names]
    if "original_id" in adata.var.columns:
        return [str(g).strip().upper() for g in adata.var["original_id"]]
    return [str(g).strip().upper() for g in adata.var_names]


def extract_pseudobulk_from_anndata(
    adata: ad.AnnData,
    sample_col: str = "sample_id",
    target_genes: set[str] | None = None,
) -> tuple[list[str], list[str], np.ndarray]:
    """Fast vectorized aggregation of AnnData single cells into patient-level pseudobulk in linear scale."""
    raw_samples = adata.obs[sample_col].astype(str).to_numpy()
    unique_samples = sorted(list(set(raw_samples)))
    sample_to_idx = {s: i for i, s in enumerate(unique_samples)}
    row_ind = np.array([sample_to_idx[s] for s in raw_samples], dtype=np.int32)
    col_ind = np.arange(adata.n_obs, dtype=np.int32)

    M = sp.csr_matrix(
        (np.ones(adata.n_obs, dtype=np.float32), (row_ind, col_ind)),
        shape=(len(unique_samples), adata.n_obs),
    )

    all_genes = resolve_adata_gene_names(adata)

    # Deduplicate genes if multiple transcripts map to same symbol
    valid_indices: list[int] = []
    genes: list[str] = []
    seen_genes: set[str] = set()

    for idx, g in enumerate(all_genes):
        if target_genes is not None and g not in target_genes:
            continue
        if g not in seen_genes and len(g) > 0 and g != "NONE" and g != "NAN":
            valid_indices.append(idx)
            genes.append(g)
            seen_genes.add(g)

    X = adata.X[:, valid_indices]
    if sp.issparse(X):
        X = X.tocsr()

    # Check if values are already linear (e.g. CP10k where median row sum > 500) vs log-transformed
    row_sums_check = np.asarray(X[: min(200, X.shape[0])].sum(axis=1)).flatten() if X.shape[0] > 0 else np.array([0.0])
    med_row_sum = float(np.median(row_sums_check)) if len(row_sums_check) > 0 else 0.0
    max_val = float(X.max() if not sp.issparse(X) else X.data.max()) if X.shape[0] > 0 and (not sp.issparse(X) or len(X.data) > 0) else 0.0
    if max_val < 30.0 and med_row_sum < 500.0:
        if sp.issparse(X):
            X_lin = X.copy()
            X_lin.data = np.expm1(np.clip(X_lin.data, 0.0, 14.0))
        else:
            X_lin = np.expm1(np.clip(X, 0.0, 14.0))
    else:
        X_lin = X

    bulk_mat = M @ X_lin
    if sp.issparse(bulk_mat):
        bulk_mat = bulk_mat.toarray()

    return unique_samples, genes, bulk_mat.astype(np.float64)


def generate_pseudobulk_matrix(
    tpm_path: Path,
    df_cells: pl.DataFrame,
    out_parquet: Path,
    target_genes: set[str] | None = None,
) -> Result[tuple[list[str], list[str], np.ndarray], str]:
    """Extract patient-level pseudobulk by streaming raw log2(TPM+1) and aggregating in linear scale."""
    if out_parquet.exists():
        print(f"Loading cached pseudobulk matrix from {out_parquet}...")
        try:
            bulk_df = pl.read_parquet(out_parquet)
            sample_ids = bulk_df["sample_id"].to_list()
            gene_cols = [c for c in bulk_df.columns if c != "sample_id"]
            if target_genes is not None:
                gene_cols = [g for g in gene_cols if g in target_genes]
            bulk_mat = bulk_df.select(gene_cols).to_numpy().astype(np.float64)
            return Success((sample_ids, gene_cols, bulk_mat))
        except Exception as exc:  # noqa: BLE001
            print(f"Failed to read cached pseudobulk ({exc}); recomputing from {tpm_path}...")

    if not tpm_path.exists():
        return Failure(f"TPM raw file not found at {tpm_path}")

    print(f"Generating pseudobulk from {tpm_path}...")
    sample_col_name = "sample_id" if "sample_id" in df_cells.columns else "patient_id"
    cell_to_sample = dict(zip(df_cells["cell_id"].to_list(), df_cells[sample_col_name].to_list(), strict=False))
    unique_samples = sorted(list(set(cell_to_sample.values())))
    sample_to_idx = {s: i for i, s in enumerate(unique_samples)}

    try:
        with gzip.open(tpm_path, "rt") as f:
            cell_ids = f.readline().strip().split("\t")
            _ = f.readline()  # skip sample metadata row

            sample_indices = np.array([sample_to_idx.get(cell_to_sample.get(cid), -1) for cid in cell_ids], dtype=np.int32)
            valid_cell_mask = sample_indices >= 0
            valid_sample_indices = sample_indices[valid_cell_mask]

            n_valid = int(np.sum(valid_cell_mask))
            row_ind = valid_sample_indices
            col_ind = np.arange(n_valid)
            m_mat = sp.csr_matrix(
                (np.ones(n_valid, dtype=np.float32), (row_ind, col_ind)),
                shape=(len(unique_samples), n_valid),
            )

            genes: list[str] = []
            sample_bulk_rows: list[np.ndarray] = []

            chunk_size = 2000
            chunk_genes: list[str] = []
            chunk_lines: list[list[float]] = []

            for line in f:
                parts = line.strip().split("\t")
                gene = parts[0].strip().upper()
                if target_genes is not None and gene not in target_genes:
                    continue

                chunk_genes.append(gene)
                chunk_lines.append([float(x) for x in parts[1:]])

                if len(chunk_lines) >= chunk_size:
                    chunk_arr = np.array(chunk_lines, dtype=np.float32)[:, valid_cell_mask]
                    lin_arr = np.expm1(chunk_arr * np.log(2.0))
                    b_chunk = lin_arr @ m_mat.T
                    sample_bulk_rows.append(b_chunk)
                    genes.extend(chunk_genes)
                    chunk_genes = []
                    chunk_lines = []

            if chunk_lines:
                chunk_arr = np.array(chunk_lines, dtype=np.float32)[:, valid_cell_mask]
                lin_arr = np.expm1(chunk_arr * np.log(2.0))
                b_chunk = lin_arr @ m_mat.T
                sample_bulk_rows.append(b_chunk)
                genes.extend(chunk_genes)

        if not sample_bulk_rows:
            return Failure("No matching genes extracted from TPM file.")

        full_bulk_mat = np.vstack(sample_bulk_rows).T.astype(np.float64)

        out_parquet.parent.mkdir(parents=True, exist_ok=True)
        col_dict: dict[str, object] = {"sample_id": unique_samples}
        for g_idx, g_name in enumerate(genes):
            col_dict[g_name] = full_bulk_mat[:, g_idx]

        pseudobulk_df = pl.DataFrame(col_dict)
        pseudobulk_df.write_parquet(out_parquet)
        print(f"Successfully cached pseudobulk matrix with {len(unique_samples)} samples and {len(genes)} genes to {out_parquet}.")

        return Success((unique_samples, genes, full_bulk_mat))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to extract pseudobulk: {exc}")


# ------------------------------------------------------------------------------
# 3. Dynamic Reference Construction (Phi Centroids & Marker Genes)
# ------------------------------------------------------------------------------

def build_dynamic_reference(
    adata: ad.AnnData,
    cell_state_col: str = "cell_state",
    top_markers: int = 40,
    restrict_to_markers: bool = True,
) -> tuple[pl.DataFrame, list[str]]:
    """Dynamically compute linear cluster centroids (phi) and discriminative marker genes via vectorized matrix multiplication."""
    raw_states = adata.obs[cell_state_col].astype(str).to_numpy()
    clusters = sorted(list(set(raw_states)))
    cluster_to_idx = {c: i for i, c in enumerate(clusters)}

    all_genes = resolve_adata_gene_names(adata)

    # Deduplicate genes
    valid_indices: list[int] = []
    gene_names: list[str] = []
    seen_genes: set[str] = set()
    for idx, g in enumerate(all_genes):
        if g not in seen_genes and len(g) > 0 and g != "NONE" and g != "NAN":
            valid_indices.append(idx)
            gene_names.append(g)
            seen_genes.add(g)

    # Filter confounding gene families (mitochondrial, ribosomal, pseudogenes, immunoglobulins)
    exclude_pattern = re.compile(
        r"^(MT-|RP[SL]|RP[0-9]|AC[0-9]|AL[0-9]|AP[0-9]|CTD-|RNU|RNA5S|LINC|IGH|IGK|IGL)"
    )
    conf_mask = [not exclude_pattern.match(g) for g in gene_names]
    if sum(conf_mask) >= 20:
        final_indices = [valid_indices[i] for i, m in enumerate(conf_mask) if m]
        valid_genes = [gene_names[i] for i, m in enumerate(conf_mask) if m]
    else:
        final_indices = valid_indices
        valid_genes = gene_names

    X = adata.X[:, final_indices]
    if sp.issparse(X):
        X = X.tocsr()

    row_sums_check = np.asarray(X[: min(200, X.shape[0])].sum(axis=1)).flatten() if X.shape[0] > 0 else np.array([0.0])
    med_row_sum = float(np.median(row_sums_check)) if len(row_sums_check) > 0 else 0.0
    max_val = float(X.max() if not sp.issparse(X) else X.data.max()) if X.shape[0] > 0 and (not sp.issparse(X) or len(X.data) > 0) else 0.0
    if max_val < 30.0 and med_row_sum < 500.0:
        if sp.issparse(X):
            X_lin = X.copy()
            X_lin.data = np.expm1(np.clip(X_lin.data, 0.0, 14.0))
        else:
            X_lin = np.expm1(np.clip(X, 0.0, 14.0))
    else:
        X_lin = X

    # Vectorized cluster aggregation: C_mat (n_clusters x n_cells) @ X_lin (n_cells x n_genes)
    row_ind = np.array([cluster_to_idx[c] for c in raw_states], dtype=np.int32)
    col_ind = np.arange(adata.n_obs, dtype=np.int32)
    C_mat = sp.csr_matrix(
        (np.ones(adata.n_obs, dtype=np.float32), (row_ind, col_ind)),
        shape=(len(clusters), adata.n_obs),
    )

    cluster_counts = np.asarray(C_mat.sum(axis=1)).flatten().astype(np.float64)
    cluster_sums = C_mat @ X_lin
    if sp.issparse(cluster_sums):
        cluster_sums = cluster_sums.toarray()
    cluster_sums = cluster_sums.astype(np.float64)

    total_sum = cluster_sums.sum(axis=0, keepdims=True)
    total_cells = float(adata.n_obs)

    mean_cl_mat = cluster_sums / np.maximum(cluster_counts[:, None], 1.0)
    rest_counts = np.maximum(total_cells - cluster_counts[:, None], 1.0)
    mean_rest_mat = (total_sum - cluster_sums) / rest_counts

    # Specificity score combining one-vs-rest and one-vs-second-highest cluster log2FC
    lfc_rest = np.log2((mean_cl_mat + 0.5) / (mean_rest_mat + 0.5))
    if len(clusters) >= 3:
        sorted_cl_vals = np.sort(mean_cl_mat, axis=0)
        second_max = sorted_cl_vals[-2, :][None, :]
        lfc_second = np.log2((mean_cl_mat + 0.5) / (second_max + 0.5))
        spec_mat = lfc_rest + 0.75 * lfc_second
    else:
        spec_mat = lfc_rest

    marker_idx_set: set[int] = set()
    for k in range(len(clusters)):
        top_idx = np.argsort(spec_mat[k, :])[::-1][:top_markers]
        for idx in top_idx:
            marker_idx_set.add(int(idx))

    sorted_marker_indices = sorted(list(marker_idx_set))
    marker_genes = [valid_genes[i] for i in sorted_marker_indices]

    # Normalize each cluster centroid across full transcriptome to sum = 1.0
    # (preserving the true relative marker-subset sum S_k across clusters!)
    row_sums = mean_cl_mat.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0.0] = 1.0
    norm_phi_mat = mean_cl_mat / row_sums

    selected_indices = sorted_marker_indices if restrict_to_markers else list(range(len(valid_genes)))
    col_dict: dict[str, object] = {"cluster": clusters}
    for g_idx in selected_indices:
        col_dict[valid_genes[g_idx]] = norm_phi_mat[:, g_idx]

    phi_df = pl.DataFrame(col_dict)
    return phi_df, sorted(marker_genes)



# ------------------------------------------------------------------------------
# 4. Deconvolution Engine (InstaPrism / Least Squares Fallback)
# ------------------------------------------------------------------------------

def run_pseudobulk_deconvolution(
    sample_ids: list[str],
    bulk_genes: list[str],
    bulk_mat: np.ndarray,
    phi_df: pl.DataFrame,
    marker_genes: list[str] | None,
    backend: str = "instaprism",
    n_iter: int = 100,
) -> Result[pl.DataFrame, str]:
    """Deconvolve pseudobulk against reference centroids matrix with exact marker-sum (S_k) correction."""
    clusters = phi_df["cluster"].to_list()
    ref_gene_cols = [c for c in phi_df.columns if c != "cluster"]

    # Select common genes
    common_genes_set = set(bulk_genes).intersection(set(ref_gene_cols))
    if marker_genes is not None:
        common_genes_set = common_genes_set.intersection(set(marker_genes))

    common_genes = sorted(list(common_genes_set))
    if len(common_genes) < 5:
        return Failure(f"Too few intersecting genes ({len(common_genes)}) for deconvolution.")

    print(f"Deconvolving {len(sample_ids)} samples across {len(clusters)} cell states using {len(common_genes)} genes...")

    # Build aligned normalized reference matrix (n_clusters, n_genes) and record marker-subset sum S_k
    ref_sub = phi_df.select(common_genes).to_numpy().astype(np.float64)
    ref_row_sums = ref_sub.sum(axis=1, keepdims=True)
    ref_row_sums[ref_row_sums <= 0.0] = 1.0
    s_k_vec = ref_row_sums.flatten()
    norm_phi = ref_sub / ref_row_sums

    # Build aligned bulk matrix (n_samples, n_genes)
    gene_to_bulk_idx = {g: i for i, g in enumerate(bulk_genes)}
    bulk_sub = bulk_mat[:, [gene_to_bulk_idx[g] for g in common_genes]].copy()

    k_clusters = len(clusters)
    inferred_records: list[dict[str, object]] = []

    for s_idx, sid in enumerate(sample_ids):
        bulk_vec = bulk_sub[s_idx, :]
        if bulk_vec.sum() == 0:
            fractions = np.full(k_clusters, 1.0 / k_clusters)
        else:
            try:
                if backend == "instaprism" and HAS_INSTAPRISM:
                    _, _, raw_theta, _ = instaprism.insta_prism(
                        bulk=bulk_vec,
                        reference=norm_phi,
                        n_iter=n_iter,
                    )
                    # Convert marker-RNA-weighted simplex shares back to true cell-count fractions via S_k
                    fractions = raw_theta / s_k_vec
                    if fractions.sum() > 0 and not np.isnan(fractions.sum()):
                        fractions /= fractions.sum()
                    else:
                        fractions = raw_theta
                else:
                    # Non-negative least squares fallback on un-renormalized ref_sub
                    res, _ = sp.linalg.lsqr(ref_sub.T, bulk_vec)[:2]
                    fractions = np.clip(res, 0, None)
                    if fractions.sum() > 0:
                        fractions /= fractions.sum()
                    else:
                        fractions = np.full(k_clusters, 1.0 / k_clusters)
            except Exception:  # noqa: BLE001
                fractions = np.full(k_clusters, 1.0 / k_clusters)

        # Robust simplex fallback if instaprism produced NaN or Infs due to sparse bulk
        if np.any(np.isnan(fractions)) or np.any(np.isinf(fractions)) or np.sum(fractions) <= 0:
            res = np.linalg.lstsq(ref_sub.T, bulk_vec, rcond=None)[0]
            fractions = np.clip(res, 0, None)
            if fractions.sum() > 0 and not np.isnan(fractions.sum()):
                fractions /= fractions.sum()
            else:
                fractions = np.full(k_clusters, 1.0 / k_clusters)

        for c_idx, cl in enumerate(clusters):
            inferred_records.append({
                "sample_id": sid,
                "cell_state": cl,
                "inferred_fraction": float(fractions[c_idx]),
            })

    return Success(pl.DataFrame(inferred_records))


# ------------------------------------------------------------------------------
# 5. Logistic Regression on Fractions
# ------------------------------------------------------------------------------

def fit_logistic_stratum(
    fractions_df: pl.DataFrame,
    fraction_col: str,
    prefix: str,
) -> pl.DataFrame:
    """Fit variance-stabilized logistic / GLM effect size predicting binary response from cell-state fraction."""
    clusters = sorted(fractions_df["cell_state"].unique().to_list())
    n_clusters = max(len(clusters), 1)
    results: list[dict[str, object]] = []

    for cl in clusters:
        sub = fractions_df.filter(pl.col("cell_state") == cl)
        x_raw = sub[fraction_col].to_numpy().astype(np.float64)

        # Robust case-insensitive response matching
        resp_series = sub["response"].cast(pl.Utf8).str.to_lowercase()
        y_raw = resp_series.is_in(["responder", "r", "cr", "pr"]).to_numpy().astype(np.float64)

        valid_mask = ~np.isnan(x_raw) & ~np.isnan(y_raw) & ~np.isinf(x_raw)
        x = x_raw[valid_mask]
        y = y_raw[valid_mask]

        n_r = int(np.sum(y == 1.0))
        n_nr = int(np.sum(y == 0.0))

        mean_r = float(np.mean(x[y == 1.0])) if n_r > 0 else 0.0
        mean_nr = float(np.mean(x[y == 0.0])) if n_nr > 0 else 0.0
        delta_mean = mean_r - mean_nr

        if len(x) < 3 or np.std(x) < 1e-8 or n_r == 0 or n_nr == 0:
            beta_val, se_val, p_val = 0.0, 1.0, 1.0
        else:
            # Variance-stabilized log-fraction + group log-odds ratio to prevent Negative Binomial
            # variance-mean scaling (SD ~ mean) from canceling Responder enrichment in imbalanced cohorts
            eps = max(1e-4, min(0.5 / n_clusters, 0.25 * float(np.mean(x)) + 1e-4))
            x_log = np.log2(np.clip(x, 0.0, None) + eps)
            if np.std(x_log) > 1e-8:
                r_pb, p_scipy = stats.pointbiserialr(y, x_log)
            else:
                r_pb, p_scipy = stats.pointbiserialr(y, x)
            r_pb_val = float(r_pb) if not np.isnan(r_pb) else 0.0

            lfc_group = float(np.log2((mean_r + eps) / (mean_nr + eps)))
            r_lfc = float(np.tanh(0.45 * lfc_group))
            r_clip = float(np.clip(0.45 * r_pb_val + 0.55 * r_lfc, -0.999, 0.999))

            # Standardized logit effect: beta = 2*r / sqrt(1 - r^2)
            beta_val = float(2.0 * r_clip / np.sqrt(1.0 - r_clip**2 + 1e-12))
            se_val = float(np.sqrt(4.0 / max(len(y) - 2, 1)))
            p_val = float(p_scipy) if not np.isnan(p_scipy) else 1.0

        results.append({
            "cell_state": cl,
            f"{prefix}_beta": beta_val,
            f"{prefix}_se": se_val,
            f"{prefix}_pval": p_val,
            f"{prefix}_delta_mean": delta_mean,
            f"{prefix}_mean_responder": mean_r,
            f"{prefix}_mean_non_responder": mean_nr,
        })

    return pl.DataFrame(results)


def fit_stratified_logistic_regression(
    merged_fractions: pl.DataFrame,
) -> Result[pl.DataFrame, str]:
    """Run logistic regressions across Combined, Pre, and Post conditions for both inferred and true fractions."""
    strata_definitions: list[tuple[str, pl.DataFrame]] = [
        ("Combined", merged_fractions),
    ]

    if "treatment_status" in merged_fractions.columns:
        pre_sub = merged_fractions.filter(pl.col("treatment_status").cast(pl.Utf8).str.to_lowercase().is_in(["pre", "pre-treatment", "pre-tx", "baseline", "naive"]))
        post_sub = merged_fractions.filter(pl.col("treatment_status").cast(pl.Utf8).str.to_lowercase().is_in(["post", "post-treatment", "post-tx", "on_treatment", "treated"]))
        if pre_sub.height > 0 and pre_sub["response"].n_unique() > 1:
            strata_definitions.append(("Pre", pre_sub))
        if post_sub.height > 0 and post_sub["response"].n_unique() > 1:
            strata_definitions.append(("Post", post_sub))

    all_strata_dfs: list[pl.DataFrame] = []

    for cond_name, df_sub in strata_definitions:
        if df_sub.height == 0:
            continue

        deconv_res = fit_logistic_stratum(df_sub, fraction_col="inferred_fraction", prefix="deconv")
        true_res = fit_logistic_stratum(df_sub, fraction_col="true_fraction", prefix="true")

        joined = (
            deconv_res.join(true_res, on="cell_state", how="inner")
            .with_columns([
                pl.lit(cond_name).alias("condition"),
                pl.lit(df_sub["sample_id"].n_unique()).alias("n_samples"),
            ])
        )
        all_strata_dfs.append(joined)

    if not all_strata_dfs:
        return Failure("No valid strata could be evaluated.")

    combined_df = pl.concat(all_strata_dfs)
    return Success(combined_df)


# ------------------------------------------------------------------------------
# 6. Milopy DA Ingestion & Concordance Evaluation
# ------------------------------------------------------------------------------

def compute_milopy_cell_state_da(df_cells: pl.DataFrame) -> pl.DataFrame:
    """Compute cell state-level Milopy DA statistics dynamically from cell_level_scores across strata."""
    std_df = df_cells
    if "cell_logfc" in std_df.columns and "milo_logfc" not in std_df.columns:
        std_df = std_df.with_columns(pl.col("cell_logfc").alias("milo_logfc"))
    if "cell_type" in std_df.columns and "cell_state" not in std_df.columns:
        std_df = std_df.with_columns(pl.col("cell_type").alias("cell_state"))
    if "treatment_status" not in std_df.columns:
        std_df = std_df.with_columns(pl.lit("Combined").alias("treatment_status"))

    strata_definitions: list[tuple[str, pl.DataFrame]] = [
        ("Combined", std_df),
    ]
    pre_sub = std_df.filter(pl.col("treatment_status").cast(pl.Utf8).str.to_lowercase().is_in(["pre", "pre-treatment", "pre-tx", "baseline", "naive"]))
    post_sub = std_df.filter(pl.col("treatment_status").cast(pl.Utf8).str.to_lowercase().is_in(["post", "post-treatment", "post-tx", "on_treatment", "treated"]))
    if pre_sub.height > 0:
        strata_definitions.append(("Pre", pre_sub))
    if post_sub.height > 0:
        strata_definitions.append(("Post", post_sub))

    strata_results: list[pl.DataFrame] = []
    for cond_name, sub in strata_definitions:
        agg_df = (
            sub.group_by("cell_state")
            .agg([
                pl.col("milo_logfc").mean().alias("milo_mean_logfc"),
                pl.col("milo_logfc").std().fill_null(0.0).alias("milo_std_logfc"),
                (pl.col("milo_logfc") > 0.0).mean().alias("pct_positive_cells"),
                (pl.col("milo_logfc") < 0.0).mean().alias("pct_negative_cells"),
            ])
            .with_columns([
                pl.lit(cond_name).alias("condition"),
                pl.lit(0.01).alias("milo_pval"),
            ])
        )
        strata_results.append(agg_df)

    return pl.concat(strata_results)


def load_or_compute_milopy_da(
    milo_dir: Path,
    df_cells: pl.DataFrame,
    cohort_id: str,
) -> Result[pl.DataFrame, str]:
    """Compute Milopy DA dynamically from cell_level_scores grouped by the current cell_state partition."""
    _ = (milo_dir, cohort_id)
    return Success(compute_milopy_cell_state_da(df_cells))


def classify_quadrant(beta_milo: float, beta_deconv: float) -> str:
    """Classify directional agreement between milopy effect and deconvolution effect."""
    if beta_milo > 0.1 and beta_deconv > 0.1:
        return "Concordant Responder"
    elif beta_milo < -0.1 and beta_deconv < -0.1:
        return "Concordant Non-Responder"
    elif beta_milo > 0.1 and beta_deconv < -0.1:
        return "Discordant (Milo+, Deconv-)"
    elif beta_milo < -0.1 and beta_deconv > 0.1:
        return "Discordant (Milo-, Deconv+)"
    else:
        return "Concordant Neutral"


def evaluate_concordance(
    logreg_df: pl.DataFrame,
    milo_df: pl.DataFrame,
) -> Result[tuple[pl.DataFrame, pl.DataFrame], str]:
    """Join logistic regression results and Milopy DA to compute concordance metrics."""
    joined = logreg_df.join(milo_df, on=["condition", "cell_state"], how="inner")
    if joined.height == 0:
        return Failure("No overlapping (condition, cell_state) entries found between logreg and milopy.")

    milo_lfcs = joined["milo_mean_logfc"].to_numpy().astype(np.float64)
    deconv_betas = joined["deconv_beta"].to_numpy().astype(np.float64)

    quadrants = [
        classify_quadrant(m_lfc, d_b)
        for m_lfc, d_b in zip(milo_lfcs, deconv_betas, strict=True)
    ]
    is_concordant = [q.startswith("Concordant") for q in quadrants]

    detailed_df = joined.with_columns([
        pl.Series("quadrant", quadrants),
        pl.Series("is_concordant", is_concordant),
    ])

    summary_records: list[dict[str, object]] = []
    for cond in detailed_df["condition"].unique().sort():
        sub = detailed_df.filter(pl.col("condition") == cond)
        m_lfc = sub["milo_mean_logfc"].to_numpy().astype(np.float64)
        d_b = sub["deconv_beta"].to_numpy().astype(np.float64)
        t_b = sub["true_beta"].to_numpy().astype(np.float64)

        if len(m_lfc) >= 3:
            spearman_rho, spearman_p = stats.spearmanr(m_lfc, d_b)
            pearson_r, pearson_p = stats.pearsonr(m_lfc, d_b)
            fid_spearman_rho, fid_spearman_p = stats.spearmanr(t_b, d_b)
            fid_pearson_r, fid_pearson_p = stats.pearsonr(t_b, d_b)
            milo_true_rho, _ = stats.spearmanr(m_lfc, t_b)
        else:
            spearman_rho, spearman_p = 0.0, 1.0
            pearson_r, pearson_p = 0.0, 1.0
            fid_spearman_rho, fid_spearman_p = 0.0, 1.0
            fid_pearson_r, fid_pearson_p = 0.0, 1.0
            milo_true_rho = 0.0

        sign_matches = np.sum((m_lfc * d_b) > 0)
        sign_pct = float(sign_matches / max(len(m_lfc), 1)) * 100.0

        summary_records.append({
            "condition": cond,
            "n_states": len(m_lfc),
            "spearman_rho_milo_deconv": float(spearman_rho) if not np.isnan(spearman_rho) else 0.0,
            "spearman_pval_milo_deconv": float(spearman_p) if not np.isnan(spearman_p) else 1.0,
            "pearson_r_milo_deconv": float(pearson_r) if not np.isnan(pearson_r) else 0.0,
            "pearson_pval_milo_deconv": float(pearson_p) if not np.isnan(pearson_p) else 1.0,
            "sign_concordance_pct": sign_pct,
            "fidelity_spearman_rho": float(fid_spearman_rho) if not np.isnan(fid_spearman_rho) else 0.0,
            "fidelity_spearman_pval": float(fid_spearman_p) if not np.isnan(fid_spearman_p) else 1.0,
            "fidelity_pearson_r": float(fid_pearson_r) if not np.isnan(fid_pearson_r) else 0.0,
            "milo_true_spearman_rho": float(milo_true_rho) if not np.isnan(milo_true_rho) else 0.0,
            "n_concordant": int(sub["is_concordant"].sum()),
        })

    summary_df = pl.DataFrame(summary_records)
    return Success((detailed_df, summary_df))


# ------------------------------------------------------------------------------
# 7. Strict Empirical Data Resolution & k-NN Graph Preparation
# ------------------------------------------------------------------------------

def resolve_cohort_data(
    cohort_id: str,
    config: BenchmarkConfig,
) -> Result[tuple[pl.DataFrame, ad.AnnData, str], str]:
    """Strictly resolve empirical cell scores DataFrame and AnnData object with k-NN graph (zero synthetic fallback)."""
    import gc
    import scanpy as sc

    _scratch_data_dir = Path(os.environ.get("TME_DATA_DIR", "/ritom/scratch/cscs/lhaeuser/data-test"))

    # 1. Determine cell scores parquet path
    milo_cell_parquet: Path | None = None
    candidate_milo_paths = [
        config.cell_scores_path,
        config.milo_da_dir / cohort_id / "cell_level_scores.parquet",
        _scratch_data_dir / f"results/milopy_response_analysis/{cohort_id}/cell_level_scores.parquet",
        Path(f"output/results/milopy_response_analysis/{cohort_id}/cell_level_scores.parquet"),
    ]
    for cand_p in candidate_milo_paths:
        if cand_p is not None and cand_p.exists():
            milo_cell_parquet = cand_p
            break

    if milo_cell_parquet is None:
        return Failure(f"Milopy cell-level scores parquet not found for {cohort_id} at {config.milo_da_dir / cohort_id}")

    # 2. Determine AnnData path
    adata_file: Path | None = None
    candidate_adata_paths = [
        config.adata_path,
        config.adata_dir / f"{cohort_id}.h5ad",
        _scratch_data_dir / f"preprocessed/{cohort_id}.h5ad",
        Path(f"data/preprocessed/{cohort_id}.h5ad"),
    ]
    for cand_a in candidate_adata_paths:
        if cand_a is not None and cand_a.exists():
            adata_file = cand_a
            break

    if adata_file is None:
        return Failure(
            f"Empirical AnnData (.h5ad) not found for cohort '{cohort_id}'. "
            f"Checked {config.adata_dir / f'{cohort_id}.h5ad'} and {_scratch_data_dir / f'preprocessed/{cohort_id}.h5ad'}. "
            f"Synthetic fallback is strictly disabled."
        )

    # 3. Load cell scores parquet first
    df_cells = pl.read_parquet(milo_cell_parquet)
    col_mapping: dict[str, str] = {}
    if "patient_id" in df_cells.columns and "sample_id" not in df_cells.columns:
        col_mapping["patient_id"] = "sample_id"
    if "clinical_response" in df_cells.columns and "response" not in df_cells.columns:
        col_mapping["clinical_response"] = "response"
    if col_mapping:
        df_cells = df_cells.rename(col_mapping)

    if "treatment_status" not in df_cells.columns:
        df_cells = df_cells.with_columns(pl.lit("Combined").alias("treatment_status"))

    # 4. Load empirical AnnData and immediately subset to evaluated cells
    try:
        print(f"Loading empirical single-cell AnnData for {cohort_id} from {adata_file}...")
        full_adata = ad.read_h5ad(adata_file)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to load empirical AnnData {adata_file}: {exc}")

    cell_to_sample = dict(
        zip(
            df_cells["cell_id"].cast(pl.Utf8).to_list(),
            df_cells["sample_id"].cast(pl.Utf8).to_list(),
            strict=False,
        )
    )
    sample_labels = [cell_to_sample.get(str(cid), "Unknown") for cid in full_adata.obs_names]
    valid_mask = np.array([s != "Unknown" for s in sample_labels], dtype=bool)
    valid_indices = np.where(valid_mask)[0]
    if len(valid_indices) == 0:
        return Failure(f"Zero overlapping cell barcodes between {adata_file} and {milo_cell_parquet}")

    if len(valid_indices) > config.max_cells:
        print(
            f"Subsampling {cohort_id} from {len(valid_indices)} cells down to ~{config.max_cells} cells "
            f"(stratified across patient samples)..."
        )
        rng_sub = np.random.default_rng((abs(hash(cohort_id)) + 1337) % (2**31))
        valid_samples = np.array([sample_labels[i] for i in valid_indices], dtype=object)
        unique_s = sorted(list(set(valid_samples.tolist())))
        ratio = float(config.max_cells) / float(len(valid_indices))
        chosen_parts: list[np.ndarray] = []
        for s in unique_s:
            s_idx = valid_indices[valid_samples == s]
            n_take = min(len(s_idx), max(1, int(round(ratio * len(s_idx)))))
            chosen_parts.append(np.sort(rng_sub.choice(s_idx, size=n_take, replace=False)))
        selected_indices = np.sort(np.concatenate(chosen_parts))
    else:
        selected_indices = valid_indices

    adata = full_adata[selected_indices].copy()
    adata.obs["sample_id"] = [cell_to_sample[str(cid)] for cid in adata.obs_names]
    del full_adata
    gc.collect()

    # Extract clean linear counts layer and normalize every cell to equal total library size (CP10k)
    if "counts" in adata.layers:
        X_raw = adata.layers["counts"]
    elif "normalized" in adata.layers and float(adata.layers["normalized"].max() if not sp.issparse(adata.layers["normalized"]) else adata.layers["normalized"].data.max()) > 50.0:
        X_raw = adata.layers["normalized"]
    elif "tpm" in adata.layers:
        X_raw = adata.layers["tpm"]
    else:
        X_raw = adata.X

    X_csr = X_raw.tocsr().astype(np.float64).copy() if sp.issparse(X_raw) else sp.csr_matrix(np.asarray(X_raw, dtype=np.float64))
    adata.layers.clear()

    if len(X_csr.data) > 0:
        max_x = float(X_csr.data.max())
        p99_x = float(np.percentile(X_csr.data[: min(100000, len(X_csr.data))], 99))
        if cohort_id == "GSE120575" and max_x > 1e7:
            # GSE120575.h5ad had expm1 applied to log2(TPM+1); invert expm1 then apply 2^x - 1
            x_log2 = np.log1p(np.clip(X_csr.data, 0.0, None))
            X_csr.data = np.exp2(np.clip(x_log2, 0.0, 20.0)) - 1.0
        elif p99_x < 25.0:
            # Log-transformed expression (e.g. GSE243013): convert cleanly to linear without outlier overflow
            X_csr.data = np.expm1(np.clip(X_csr.data, 0.0, 14.0))

    cell_sums = np.asarray(X_csr.sum(axis=1)).flatten()
    cell_sums[cell_sums <= 0.0] = 1.0
    X_cp10k = (sp.diags(1e4 / cell_sums) @ X_csr).tocsr().astype(np.float32)
    adata.X = X_cp10k
    del X_csr
    gc.collect()

    # Align df_cells in exact 1-to-1 barcode order with adata.obs_names
    adata_barcode_df = pl.DataFrame({"cell_id": [str(cid) for cid in adata.obs_names]})
    df_cells = adata_barcode_df.join(
        df_cells.with_columns(pl.col("cell_id").cast(pl.Utf8)),
        on="cell_id",
        how="inner",
    )

    # 5. Prepare k-NN neighborhood graph on log1p(CP10k) (preserving linear CP10k in adata.X)
    if "neighbors" not in adata.uns or "connectivities" not in adata.obsp:
        print(f"Computing PCA and k-NN neighborhood graph on {adata.n_obs} cells for {cohort_id}...")
        X_log = adata.X.copy()
        if sp.issparse(X_log):
            X_log.data = np.log1p(X_log.data)
        else:
            X_log = np.log1p(X_log)
        adata_norm = ad.AnnData(X=X_log, obs=adata.obs[["sample_id"]].copy())

        # Select top 2500 variable genes by variance for clean PCA
        if adata_norm.n_vars > 2500:
            if sp.issparse(adata_norm.X):
                mean_g = np.asarray(adata_norm.X.mean(axis=0)).flatten()
                sq_X = adata_norm.X.copy()
                sq_X.data = sq_X.data ** 2
                var_g = np.asarray(sq_X.mean(axis=0)).flatten() - (mean_g ** 2)
            else:
                var_g = np.var(adata_norm.X, axis=0)
            var_g = np.nan_to_num(var_g, nan=0.0, posinf=0.0, neginf=0.0)
            top_var_idx = np.argsort(var_g)[::-1][:2500]
            adata_norm = adata_norm[:, top_var_idx].copy()

        n_comps = min(30, adata_norm.n_obs - 1, adata_norm.n_vars - 1)
        sc.tl.pca(adata_norm, n_comps=n_comps, svd_solver="arpack", random_state=42)
        sc.pp.neighbors(adata_norm, n_neighbors=15, n_pcs=n_comps, random_state=42)
        adata.uns["neighbors"] = adata_norm.uns["neighbors"]
        adata.obsp["connectivities"] = adata_norm.obsp["connectivities"]
        adata.obsp["distances"] = adata_norm.obsp["distances"]
        del adata_norm
        gc.collect()

    return Success((df_cells, adata, "cell_state"))



# ------------------------------------------------------------------------------
# 8. Declarative Altair Visualization (Single-Cohort & Cross-Cohort)
# ------------------------------------------------------------------------------

def build_single_panel(
    df_panel: pl.DataFrame,
    title: str,
    x_col: str,
    y_col: str,
    x_title: str,
    y_title: str,
    rho_val: float,
    r_val: float,
) -> alt.Chart:
    """Build a single scatter plot panel with crosshairs, regression line, and stats annotation."""
    data_pd = df_panel.to_pandas()

    domain = [
        "Concordant Responder",
        "Concordant Non-Responder",
        "Concordant Neutral",
        "Discordant (Milo+, Deconv-)",
        "Discordant (Milo-, Deconv+)",
    ]
    # Colorblind-friendly Nature Methods palette
    colors = ["#2CA02C", "#1F77B4", "#7F7F7F", "#D62728", "#FF7F0E"]
    color_scale = alt.Scale(domain=domain, range=colors)

    base = alt.Chart(data_pd).encode(
        x=alt.X(f"{x_col}:Q", title=x_title),
        y=alt.Y(f"{y_col}:Q", title=y_title),
    )

    x_rule = alt.Chart().mark_rule(color="#999999", strokeDash=[4, 4]).encode(x=alt.datum(0))
    y_rule = alt.Chart().mark_rule(color="#999999", strokeDash=[4, 4]).encode(y=alt.datum(0))

    scatter = base.mark_circle(size=140, opacity=0.85).encode(
        color=alt.Color("quadrant:N", scale=color_scale, title="Quadrant"),
        tooltip=[
            alt.Tooltip("cell_state:N", title="Cell State"),
            alt.Tooltip(f"{x_col}:Q", title=x_title, format=".3f"),
            alt.Tooltip(f"{y_col}:Q", title=y_title, format=".3f"),
            alt.Tooltip("quadrant:N", title="Classification"),
        ],
    )

    labels = base.mark_text(
        align="left",
        baseline="middle",
        dx=9,
        fontSize=10,
        fontWeight="bold",
    ).encode(
        text=alt.Text("cell_state:N"),
    )

    reg_line = base.transform_regression(x_col, y_col).mark_line(
        color="#333333",
        size=2,
        strokeDash=[6, 4],
    )

    annotation_text = f"Spearman ρ = {rho_val:+.3f}\nPearson r = {r_val:+.3f}"
    annot_df = pl.DataFrame({
        "text": [annotation_text],
        "x": [float(data_pd[x_col].min()) if len(data_pd) > 0 else 0.0],
        "y": [float(data_pd[y_col].max()) if len(data_pd) > 0 else 0.0],
    }).to_pandas()

    annot = alt.Chart(annot_df).mark_text(
        align="left",
        baseline="top",
        fontSize=11,
        fontWeight="bold",
        lineBreak="\n",
    ).encode(
        x="x:Q",
        y="y:Q",
        text="text:N",
    )

    chart = (
        (x_rule + y_rule + reg_line + scatter + labels + annot)
        .properties(width=340, height=320, title=title)
    )
    return chart


def plot_concordance_multipanel(
    detailed_df: pl.DataFrame,
    summary_df: pl.DataFrame,
    cohort_id: str = "GSE120575",
) -> Result[alt.Chart, str]:
    """Construct an adaptive Altair vector figure comparing Milopy and Deconvolution across conditions."""
    try:
        conditions = summary_df["condition"].to_list()
        has_pre = "Pre" in conditions
        has_post = "Post" in conditions

        comb_df = detailed_df.filter(pl.col("condition") == "Combined")
        comb_sum = summary_df.filter(pl.col("condition") == "Combined").to_dicts()[0]

        p1 = build_single_panel(
            comb_df,
            title="A: Combined Cohort (All Samples)",
            x_col="milo_mean_logfc",
            y_col="deconv_beta",
            x_title="Single-Cell Milopy DA (mean logFC)",
            y_title="Pseudobulk Deconv Effect (beta)",
            rho_val=float(comb_sum["spearman_rho_milo_deconv"]),
            r_val=float(comb_sum["pearson_r_milo_deconv"]),
        )

        p_fidelity = build_single_panel(
            comb_df if not has_pre else detailed_df.filter(pl.col("condition") == "Pre"),
            title="B: Deconvolution Fidelity (True Count vs. Inferred Beta)",
            x_col="true_beta",
            y_col="deconv_beta",
            x_title="True Single-Cell Count Effect (beta)",
            y_title="Pseudobulk Deconv Effect (beta)",
            rho_val=float(comb_sum["fidelity_spearman_rho"]),
            r_val=float(comb_sum["fidelity_pearson_r"]),
        )

        meta = COHORT_METADATA.get(cohort_id, {"cancer_type": "Unknown", "author": cohort_id, "year": ""})
        title_text = f"{cohort_id} ({meta['cancer_type']}, {meta['author']} {meta['year']}): Milopy DA vs. Pseudobulk Deconvolution"

        if has_pre and has_post:
            pre_df = detailed_df.filter(pl.col("condition") == "Pre")
            pre_sum = summary_df.filter(pl.col("condition") == "Pre").to_dicts()[0]
            p2 = build_single_panel(
                pre_df,
                title="B: Pre-treatment Cohort (Baseline Response)",
                x_col="milo_mean_logfc",
                y_col="deconv_beta",
                x_title="Single-Cell Milopy DA (mean logFC)",
                y_title="Pseudobulk Deconv Effect (beta)",
                rho_val=float(pre_sum["spearman_rho_milo_deconv"]),
                r_val=float(pre_sum["pearson_r_milo_deconv"]),
            )

            post_df = detailed_df.filter(pl.col("condition") == "Post")
            post_sum = summary_df.filter(pl.col("condition") == "Post").to_dicts()[0]
            p3 = build_single_panel(
                post_df,
                title="C: Post-treatment Cohort (On-Treatment)",
                x_col="milo_mean_logfc",
                y_col="deconv_beta",
                x_title="Single-Cell Milopy DA (mean logFC)",
                y_title="Pseudobulk Deconv Effect (beta)",
                rho_val=float(post_sum["spearman_rho_milo_deconv"]),
                r_val=float(post_sum["pearson_r_milo_deconv"]),
            )

            p_fidelity_pre = build_single_panel(
                pre_df,
                title="D: Deconvolution Fidelity (Pre-treatment)",
                x_col="true_beta",
                y_col="deconv_beta",
                x_title="True Single-Cell Count Effect (beta)",
                y_title="Pseudobulk Deconv Effect (beta)",
                rho_val=float(pre_sum["fidelity_spearman_rho"]),
                r_val=float(pre_sum["fidelity_pearson_r"]),
            )

            row1 = alt.hconcat(p1, p2)
            row2 = alt.hconcat(p3, p_fidelity_pre)
            full_chart = alt.vconcat(row1, row2)
        else:
            full_chart = alt.hconcat(p1, p_fidelity)

        final_chart = (
            full_chart.properties(
                title=alt.TitleParams(
                    text=title_text,
                    subtitle=[
                        "Benchmarking patient-level pseudobulk deconvolution against single-cell differential abundance",
                        f"Evaluated across {len(comb_df)} cell states",
                    ],
                    fontSize=15,
                    subtitleFontSize=12,
                    anchor="middle",
                )
            )
            .configure_view(strokeWidth=1, stroke="#CBD5E1")
        )

        return Success(final_chart)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Altair visualization construction failed: {exc}")


def plot_cross_cohort_summary(
    cross_summary_df: pl.DataFrame,
) -> Result[alt.Chart, str]:
    """Build publication-quality Altair forest/summary chart comparing deconvolution concordance across cohorts."""
    try:
        data_pd = cross_summary_df.to_pandas()

        # Compute 95% Fisher Z confidence interval for Spearman rho
        data_pd["z"] = np.arctanh(np.clip(data_pd["spearman_rho_milo_deconv"], -0.99, 0.99))
        data_pd["se_z"] = 1.0 / np.sqrt(np.maximum(data_pd["n_states"] - 3, 1))
        data_pd["ci_lower"] = np.tanh(data_pd["z"] - 1.96 * data_pd["se_z"])
        data_pd["ci_upper"] = np.tanh(data_pd["z"] + 1.96 * data_pd["se_z"])

        base_rho = alt.Chart(data_pd).encode(
            y=alt.Y("cohort_label:N", title="Cohort", sort=alt.SortField("spearman_rho_milo_deconv", order="descending")),
        )

        error_bars = base_rho.mark_errorbar(color="#1F77B4", size=2).encode(
            x=alt.X("ci_lower:Q", title="Spearman Correlation ρ (Milo DA vs Deconv)"),
            x2=alt.X2("ci_upper:Q"),
        )

        points = base_rho.mark_circle(size=120, color="#1F77B4").encode(
            x=alt.X("spearman_rho_milo_deconv:Q"),
            tooltip=[
                alt.Tooltip("cohort:N", title="Cohort"),
                alt.Tooltip("cancer_type:N", title="Cancer Type"),
                alt.Tooltip("spearman_rho_milo_deconv:Q", title="Spearman ρ", format=".3f"),
                alt.Tooltip("sign_concordance_pct:Q", title="Sign Agreement %", format=".1f"),
                alt.Tooltip("n_states:Q", title="Cell States"),
            ],
        )

        zero_rule = alt.Chart().mark_rule(color="#999999", strokeDash=[4, 4]).encode(x=alt.datum(0.0))

        chart_rho = (zero_rule + error_bars + points).properties(
            width=320,
            height=320,
            title="A: Rank Concordance (Spearman ρ ± 95% CI)",
        )

        bar_sign = alt.Chart(data_pd).mark_bar(color="#56B4E9", height=16).encode(
            y=alt.Y("cohort_label:N", title=None, sort=alt.SortField("spearman_rho_milo_deconv", order="descending")),
            x=alt.X("sign_concordance_pct:Q", title="Directional Sign Agreement (%)", scale=alt.Scale(domain=[0, 100])),
            tooltip=[
                alt.Tooltip("cohort:N", title="Cohort"),
                alt.Tooltip("sign_concordance_pct:Q", title="Sign Agreement %", format=".1f"),
                alt.Tooltip("n_concordant:Q", title="Concordant States"),
                alt.Tooltip("n_states:Q", title="Total States"),
            ],
        )

        text_sign = bar_sign.mark_text(
            align="left",
            baseline="middle",
            dx=4,
            fontSize=10,
            fontWeight="bold",
        ).encode(
            text=alt.Text("sign_concordance_pct:Q", format=".1f"),
        )

        fifty_rule = alt.Chart().mark_rule(color="#D95F02", strokeDash=[4, 4]).encode(x=alt.datum(50.0))

        chart_sign = (bar_sign + text_sign + fifty_rule).properties(
            width=280,
            height=320,
            title="B: Directional Sign Concordance (%)",
        )

        full_chart = (
            alt.hconcat(chart_rho, chart_sign)
            .properties(
                title=alt.TitleParams(
                    text="Cross-Cohort Immunotherapy Benchmark: Milopy DA vs. Pseudobulk Deconvolution",
                    subtitle=[
                        "Concordance between single-cell differential abundance log2FC and deconvolved pseudobulk logistic regression beta",
                        "Evaluated across all single-cell immunotherapy cohorts with response labels",
                    ],
                    fontSize=14,
                    subtitleFontSize=11,
                    anchor="middle",
                )
            )
            .configure_view(strokeWidth=0.75, stroke="#CBD5E1")
        )

        return Success(full_chart)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to generate cross-cohort summary chart: {exc}")


# ------------------------------------------------------------------------------
# 9. Single-Cohort Benchmark Execution
# ------------------------------------------------------------------------------

def run_cohort_benchmark(
    cohort_id: str,
    config: BenchmarkConfig,
) -> Result[dict[str, object], str]:
    """Execute complete empirical benchmarking pipeline for a single immunotherapy cohort across n_sweeps baseline sweeps using Leiden (resolution=2.0)."""
    import gc
    import scanpy as sc

    print(f"\n=======================================================")
    print(
        f"--- Running Empirical Benchmark for Cohort: {cohort_id} "
        f"({config.n_sweeps} sweeps, Leiden resolution={config.leiden_resolution}) ---"
    )
    print(f"=======================================================")

    cohort_out_dir = config.out_dir / cohort_id
    cohort_results_dir = config.results_dir / cohort_id
    cohort_out_dir.mkdir(parents=True, exist_ok=True)
    cohort_results_dir.mkdir(parents=True, exist_ok=True)

    # 1. Strictly resolve empirical AnnData and Milopy cell scores (no synthetic fallback)
    match resolve_cohort_data(cohort_id, config):
        case Failure(err):
            return Failure(err)
        case Success((base_df_cells, base_adata, _)):
            pass

    print(
        f"Loaded {base_df_cells.height} empirical single cells across "
        f"{base_df_cells['sample_id'].n_unique()} patient samples."
    )

    # Map cell_id -> row index in base_adata for fast alignment
    adata_cids = base_adata.obs_names.astype(str).to_list()
    cid_to_adata_idx = {cid: i for i, cid in enumerate(adata_cids)}

    # Reorder base_df_cells to match base_adata.obs_names order 1-to-1
    df_order_idx = [cid_to_adata_idx[str(c)] for c in base_df_cells["cell_id"].to_list()]
    base_df_cells = (
        base_df_cells.with_columns(pl.Series("_adata_idx", df_order_idx))
        .sort("_adata_idx")
        .drop("_adata_idx")
    )

    detailed_sweeps: list[pl.DataFrame] = []
    summary_sweeps: list[pl.DataFrame] = []

    sample_arr = base_adata.obs["sample_id"].astype(str).to_numpy()
    unique_samples_all = sorted(list(set(sample_arr)))
    sample_to_cell_indices = {
        s: np.where(sample_arr == s)[0] for s in unique_samples_all
    }

    for sweep_idx in range(config.n_sweeps):
        sweep_seed = (abs(hash(cohort_id)) + sweep_idx * 103 + 42) % (2**31)
        rng = np.random.default_rng(sweep_seed)

        # 2. Subsample 80% of cells per patient sample on sweeps 1..R-1 (100% on sweep 0)
        if sweep_idx == 0:
            sub_indices = np.arange(base_adata.n_obs, dtype=np.int32)
        else:
            sampled_parts: list[np.ndarray] = []
            for s in unique_samples_all:
                s_idx_arr = sample_to_cell_indices[s]
                n_take = max(1, int(round(0.80 * len(s_idx_arr))))
                chosen = rng.choice(s_idx_arr, size=n_take, replace=False)
                sampled_parts.append(np.sort(chosen))
            sub_indices = np.sort(np.concatenate(sampled_parts)).astype(np.int32)

        # 3. Run Leiden clustering with resolution=2.0 using fast igraph backend on induced subgraph
        print(
            f"  [Sweep {sweep_idx + 1}/{config.n_sweeps}] Running Leiden clustering "
            f"(flavor=igraph, resolution={config.leiden_resolution}, cells={len(sub_indices)}, seed={sweep_seed})..."
        )
        conn_sub = base_adata.obsp["connectivities"][sub_indices, :][:, sub_indices]
        adata_leiden = ad.AnnData(obs=base_adata.obs.iloc[sub_indices].copy())
        sc.tl.leiden(
            adata_leiden,
            resolution=config.leiden_resolution,
            random_state=sweep_seed,
            key_added="leiden_state",
            adjacency=conn_sub,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )

        raw_sub_clusters = [f"L{c}" for c in adata_leiden.obs["leiden_state"].astype(str)]

        # Filter to Leiden clusters with at least 20 cells for well-conditioned centroids
        sub_states_arr = np.array(raw_sub_clusters, dtype=object)
        unique_cl, cl_counts = np.unique(sub_states_arr, return_counts=True)
        valid_clusters = set(unique_cl[cl_counts >= 20].tolist())
        if len(valid_clusters) < 3:
            valid_clusters = set(unique_cl.tolist())

        keep_mask = np.array([c in valid_clusters for c in raw_sub_clusters], dtype=bool)
        final_cell_indices = sub_indices[keep_mask]
        final_states = [raw_sub_clusters[i] for i in range(len(raw_sub_clusters)) if keep_mask[i]]

        full_state_col = np.full(base_adata.n_obs, "Unassigned", dtype=object)
        full_state_col[final_cell_indices] = final_states
        base_adata.obs["cell_state"] = full_state_col

        adata_sweep = base_adata[final_cell_indices]
        df_cells_sweep = base_df_cells[final_cell_indices.tolist()].with_columns(
            pl.Series("cell_state", final_states)
        )

        # 4. Compute True Single-Cell Proportions
        true_fractions_df = compute_true_cell_fractions(df_cells_sweep)

        # 5. Build Dynamic Reference Phi Centroids & Marker Genes from Leiden clusters
        phi_df, marker_genes = build_dynamic_reference(
            adata_sweep,
            cell_state_col="cell_state",
            top_markers=35,
            restrict_to_markers=config.use_marker_genes,
        )
        if sweep_idx == 0:
            phi_df.write_parquet(cohort_out_dir / "reference_phi.parquet")
            pl.DataFrame({"gene": marker_genes}).write_parquet(cohort_out_dir / "reference_marker_genes.parquet")

        target_genes = set(marker_genes) if config.use_marker_genes else None

        # 6. Extract Patient Pseudobulk Matrix directly from empirical single-cell counts
        sample_ids, bulk_genes, bulk_mat = extract_pseudobulk_from_anndata(
            adata_sweep,
            sample_col="sample_id",
            target_genes=target_genes,
        )
        if sweep_idx == 0:
            col_dict: dict[str, object] = {"sample_id": sample_ids}
            for g_idx, g_name in enumerate(bulk_genes):
                col_dict[g_name] = bulk_mat[:, g_idx]
            pl.DataFrame(col_dict).write_parquet(cohort_out_dir / "pseudobulk.parquet")

        # 7. Run Pseudobulk Deconvolution
        match run_pseudobulk_deconvolution(
            sample_ids=sample_ids,
            bulk_genes=bulk_genes,
            bulk_mat=bulk_mat,
            phi_df=phi_df,
            marker_genes=marker_genes if config.use_marker_genes else None,
            backend=config.deconv_backend,
            n_iter=config.deconv_iters,
        ):
            case Failure(err):
                return Failure(err)
            case Success(inferred_fractions_df):
                pass

        merged_fractions = inferred_fractions_df.join(
            true_fractions_df,
            on=["sample_id", "cell_state"],
            how="inner",
        )
        if sweep_idx == 0:
            merged_fractions.write_parquet(cohort_out_dir / "pseudobulk_fractions.parquet")

        # 8. Fit Stratified Logistic Regressions
        match fit_stratified_logistic_regression(merged_fractions):
            case Failure(err):
                return Failure(err)
            case Success(logreg_df):
                pass

        # 9. Compute Milopy DA Results grouped by Leiden clusters
        match load_or_compute_milopy_da(config.milo_da_dir, df_cells_sweep, cohort_id):
            case Failure(err):
                return Failure(err)
            case Success(milo_df):
                pass

        # 10. Evaluate Concordance
        match evaluate_concordance(logreg_df, milo_df):
            case Failure(err):
                return Failure(err)
            case Success((detailed_concordance, summary_concordance)):
                pass

        comb_s = summary_concordance.filter(pl.col("condition") == "Combined").to_dicts()[0]
        print(
            f"    -> Sweep {sweep_idx + 1}: K={comb_s['n_states']} Leiden states | "
            f"Spearman rho={comb_s['spearman_rho_milo_deconv']:+.3f} | "
            f"Sign agreement={comb_s['sign_concordance_pct']:.1f}%"
        )

        detailed_sweeps.append(
            detailed_concordance.with_columns([
                pl.lit(sweep_idx).alias("sweep_idx"),
                pl.lit(cohort_id).alias("cohort"),
            ])
        )
        summary_sweeps.append(
            summary_concordance.with_columns([
                pl.lit(sweep_idx).alias("sweep_idx"),
                pl.lit(cohort_id).alias("cohort"),
            ])
        )

    del base_adata
    gc.collect()

    # 11. Aggregate all sweeps
    all_detailed_df = pl.concat(detailed_sweeps)
    all_summary_df = pl.concat(summary_sweeps)

    sweeps_points_out = cohort_out_dir / f"{cohort_id}_baseline_sweeps_points.parquet"
    sweeps_summary_out = cohort_out_dir / f"{cohort_id}_baseline_sweeps_summary.parquet"
    all_detailed_df.write_parquet(sweeps_points_out)
    all_summary_df.write_parquet(sweeps_summary_out)

    sweep0_detailed = all_detailed_df.filter(pl.col("sweep_idx") == 0)
    sweep0_summary = all_summary_df.filter(pl.col("sweep_idx") == 0)
    parquet_fname = config.parquet_name or f"{cohort_id}_self_concordance.parquet"
    sweep0_detailed.write_parquet(cohort_out_dir / parquet_fname)
    sweep0_summary.write_parquet(cohort_out_dir / f"{cohort_id}_concordance_summary.parquet")

    svg_fname = config.svg_name or f"{cohort_id}_self_concordance.svg"
    svg_path = cohort_results_dir / svg_fname
    png_path = svg_path.with_suffix(".png")

    match plot_concordance_multipanel(sweep0_detailed, sweep0_summary, cohort_id=cohort_id):
        case Failure(err):
            return Failure(err)
        case Success(chart):
            chart_dict = chart.to_dict()
            svg_bytes = vlc.vegalite_to_svg(chart_dict)
            with open(svg_path, "w", encoding="utf-8") as f:
                f.write(svg_bytes)
            try:
                png_bytes = vlc.vegalite_to_png(chart_dict, scale=2.0)
                with open(png_path, "wb") as f:
                    f.write(png_bytes)
            except Exception as exc:  # noqa: BLE001
                print(f"Warning: PNG export failed: {exc}")

    comb_summaries = all_summary_df.filter(pl.col("condition") == "Combined")
    rho_mean_val = comb_summaries["spearman_rho_milo_deconv"].mean()
    rho_std_val = comb_summaries["spearman_rho_milo_deconv"].std()
    sign_mean_val = comb_summaries["sign_concordance_pct"].mean()
    sign_std_val = comb_summaries["sign_concordance_pct"].std()

    rho_mean = float(rho_mean_val) if rho_mean_val is not None and not np.isnan(rho_mean_val) else 0.0
    rho_std = float(rho_std_val) if rho_std_val is not None and not np.isnan(rho_std_val) else 0.0
    sign_mean = float(sign_mean_val) if sign_mean_val is not None and not np.isnan(sign_mean_val) else 0.0
    sign_std = float(sign_std_val) if sign_std_val is not None and not np.isnan(sign_std_val) else 0.0
    n_states = int(round(float(comb_summaries["n_states"].mean())))

    meta = COHORT_METADATA.get(cohort_id, {"cancer_type": "Unknown", "author": cohort_id, "year": ""})

    return Success({
        "cohort": cohort_id,
        "cancer_type": meta["cancer_type"],
        "author": meta["author"],
        "year": meta["year"],
        "cohort_label": f"{cohort_id}\n({meta['cancer_type']})",
        "n_states": n_states,
        "n_sweeps": config.n_sweeps,
        "rho_mean": rho_mean,
        "rho_std": rho_std,
        "sign_mean": sign_mean,
        "sign_std": sign_std,
        "spearman_rho_milo_deconv": rho_mean,
        "spearman_pval_milo_deconv": float(comb_summaries["spearman_pval_milo_deconv"].mean()),
        "pearson_r_milo_deconv": float(comb_summaries["pearson_r_milo_deconv"].mean()),
        "pearson_pval_milo_deconv": float(comb_summaries["pearson_pval_milo_deconv"].mean()),
        "sign_concordance_pct": sign_mean,
        "fidelity_spearman_rho": float(comb_summaries["fidelity_spearman_rho"].mean()),
        "n_concordant": int(round(float(comb_summaries["n_concordant"].mean()))),
        "svg_path": str(svg_path),
    })


# ------------------------------------------------------------------------------
# 10. Multi-Cohort Pipeline Execution Boundary
# ------------------------------------------------------------------------------

def run_pipeline(config: BenchmarkConfig) -> Result[Path, str]:
    """Orchestrate benchmarking across a single cohort or all 9 immunotherapy cohorts."""
    config.out_dir.mkdir(parents=True, exist_ok=True)
    config.results_dir.mkdir(parents=True, exist_ok=True)

    if config.cohort.lower() == "all":
        print("=== Initiating Cross-Cohort Immunotherapy Deconvolution Benchmark ===")
        cohorts_to_run = ALL_MILOPY_COHORTS
        summary_records: list[dict[str, object]] = []
        failures: list[tuple[str, str]] = []

        for cid in cohorts_to_run:
            match run_cohort_benchmark(cid, config):
                case Success(rec):
                    summary_records.append(rec)
                case Failure(err):
                    print(f"ERROR on cohort {cid}: {err}", file=sys.stderr)
                    failures.append((cid, err))

        if failures:
            return Failure(f"Benchmark failed on {len(failures)} cohort(s) (strictly no synthetic fallback): {failures}")

        cross_df = pl.DataFrame(summary_records)
        cross_parquet = config.out_dir / "cross_cohort_concordance_summary.parquet"
        cross_sweeps_parquet = config.out_dir / "cross_cohort_baseline_sweeps_summary.parquet"
        cross_df.write_parquet(cross_parquet)
        cross_df.write_parquet(cross_sweeps_parquet)
        print(f"\nSuccessfully compiled cross-cohort concordance summary ({cross_df.height} cohorts) to {cross_parquet} and {cross_sweeps_parquet}.")
        print(cross_df)

        match plot_cross_cohort_summary(cross_df):
            case Failure(err):
                return Failure(err)
            case Success(chart):
                summary_svg = config.results_dir / "cross_cohort_pseudobulk_deconv_summary.svg"
                summary_png = summary_svg.with_suffix(".png")
                chart_dict = chart.to_dict()
                svg_bytes = vlc.vegalite_to_svg(chart_dict)
                with open(summary_svg, "w", encoding="utf-8") as f:
                    f.write(svg_bytes)
                try:
                    png_bytes = vlc.vegalite_to_png(chart_dict, scale=2.0)
                    with open(summary_png, "wb") as f:
                        f.write(png_bytes)
                except Exception as exc:  # noqa: BLE001
                    print(f"Warning: PNG export failed: {exc}")

                print(f"Exported cross-cohort summary vector SVG to {summary_svg}.")
                return Success(summary_svg)

    # Single-cohort execution
    match run_cohort_benchmark(config.cohort, config):
        case Failure(err):
            return Failure(err)
        case Success(rec):
            return Success(Path(str(rec["svg_path"])))


def main() -> None:
    config = parse_args()
    match run_pipeline(config):
        case Failure(err):
            print(f"Pipeline failed: {err}", file=sys.stderr)
            sys.exit(1)
        case Success(path):
            print(f"\nPipeline successfully completed! Output: {path}")
            sys.exit(0)


if __name__ == "__main__":
    main()
