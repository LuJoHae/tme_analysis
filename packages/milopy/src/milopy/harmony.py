"""Harmony batch integration, graph construction, and neighborhood metadata annotation for Milopy.
"""

from __future__ import annotations

from typing import Sequence
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse as sp

try:
    import harmonypy as hm
    HAS_HARMONY = True
except ImportError:
    HAS_HARMONY = False

from .prevalence import ReplicatePrevalenceConfig, evaluate_replicate_prevalence
from returns.result import Failure, Success


def ensure_pca_and_graph(
    adata: ad.AnnData,
    patient_col: str,
    k: int = 30,
    d: int = 30,
    use_harmony: bool = True,
    seed: int = 42,
) -> str:
    """Ensures PCA coordinates, optional harmonypy batch correction, kNN graph connectivity, and fresh UMAP coordinates exist in AnnData.

    Returns the representation basis used for graph construction ('X_pca_harmony' or 'X_pca').
    """
    if "X_pca" not in adata.obsm:
        max_val = (
            adata.X.max()
            if not sp.issparse(adata.X)
            else adata.X.data.max()
            if adata.X.nnz > 0
            else 0
        )
        if max_val > 50:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)

        # Select highly variable genes if many genes exist to make PCA fast and memory-efficient
        if adata.n_vars > 2000:
            sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=False)
            has_hvg = "highly_variable" in adata.var.columns
        else:
            has_hvg = False

        n_comps = min(d, adata.n_vars - 1, adata.n_obs - 1)
        if has_hvg:
            sc.pp.pca(adata, n_comps=n_comps, mask_var="highly_variable", zero_center=False)
        else:
            sc.pp.pca(adata, n_comps=n_comps, zero_center=False)

    basis_rep = "X_pca"
    if use_harmony and HAS_HARMONY and patient_col in adata.obs.columns:
        n_pts = adata.obs[patient_col].nunique()
        if n_pts > 1:
            try:
                ho = hm.run_harmony(
                    adata.obsm["X_pca"],
                    adata.obs,
                    patient_col,
                    max_iter_harmony=10,
                    random_state=seed,
                    verbose=False,
                )
                adata.obsm["X_pca_raw"] = adata.obsm["X_pca"].copy()
                adata.obsm["X_pca_harmony"] = ho.Z_corr
                # milopy.core.make_nhoods uses adata.obsm['X_pca'] for refinement,
                # so setting X_pca to the harmonized embedding guarantees both graph construction
                # and neighborhood median refinement occur in harmonized latent space.
                adata.obsm["X_pca"] = ho.Z_corr
                basis_rep = "X_pca_harmony"
            except Exception:
                basis_rep = "X_pca"

    d_actual = min(d, adata.obsm["X_pca"].shape[1])
    k_actual = min(k, adata.n_obs - 1)
    if "connectivities" not in adata.obsp or "distances" not in adata.obsp or basis_rep == "X_pca_harmony":
        sc.pp.neighbors(adata, n_neighbors=k_actual, n_pcs=d_actual, use_rep="X_pca")

    # Compute fresh UMAP coordinates for the analyzed cell subset based on the graph
    sc.tl.umap(adata, min_dist=0.3, spread=1.0, random_state=seed)
    return basis_rep


def annotate_nhoods_with_metadata(
    adata: ad.AnnData,
    res_df: pd.DataFrame,
    patient_col: str,
    design_col: str,
    min_replicates: int = 2,
    min_prevalence_frac: float = 0.05,
    fdr_threshold: float = 0.1,
    dataset_col: str | None = None,
    permutation_pvalues: Sequence[float] | None = None,
    permutation_fdr: Sequence[float] | None = None,
) -> pd.DataFrame:
    """Identifies dominant cell type, purity, and evaluates biological donor prevalence and permutation FDR.

    Annotates total patients, responder patients, and non-responder patients present in
    each neighborhood. Evaluates replicate prevalence: significant hits (FDR < fdr_threshold)
    must be supported by >= min_replicates independent patients in the enriched arm and
    exceed minimum prevalence fraction, preventing single-patient private spikes from being reported.
    """
    res_df = res_df.copy()
    cell_type_col = None
    for cand in ("cell_type", "celltype", "CellType", "cell_type_annotation", "cluster"):
        if cand in adata.obs.columns:
            cell_type_col = cand
            break

    if "nhoods" not in adata.obsm:
        res_df["Nhood_CellType"] = "Unknown"
        res_df["Nhood_CellType_Purity"] = 1.0
        res_df["n_patients_total"] = 0
        res_df["n_patients_responder"] = 0
        res_df["n_patients_non_responder"] = 0
        res_df["neg_log10_fdr"] = -np.log10(res_df["FDR"].clip(lower=1e-300))
        res_df["is_significant"] = False
        res_df["status"] = "Not Significant"
        return res_df

    nhoods_mat = adata.obsm["nhoods"].tocsc()
    n_nhoods = nhoods_mat.shape[1]
    cell_types = (
        adata.obs[cell_type_col].astype(str).values
        if cell_type_col is not None
        else np.array(["Unknown"] * adata.n_obs)
    )

    majority_types: list[str] = []
    purities: list[float] = []

    for i in range(n_nhoods):
        idx = nhoods_mat[:, i].nonzero()[0]
        if len(idx) == 0:
            majority_types.append("Empty")
            purities.append(0.0)
            continue

        sample_types = cell_types[idx]
        vals, counts = np.unique(sample_types, return_counts=True)
        top_idx = int(np.argmax(counts))
        majority_types.append(str(vals[top_idx]))
        purities.append(float(counts[top_idx] / len(idx)))

    res_df["Nhood_CellType"] = majority_types
    res_df["Nhood_CellType_Purity"] = purities

    prev_res = evaluate_replicate_prevalence(
        adata=adata,
        res_df=res_df,
        patient_col=patient_col,
        design_col=design_col,
        dataset_col=dataset_col,
        permutation_pvalues=permutation_pvalues,
        permutation_fdr=permutation_fdr,
        config=ReplicatePrevalenceConfig(
            min_replicates=min_replicates,
            min_prevalence_frac=min_prevalence_frac,
            fdr_threshold=fdr_threshold,
            require_permutation_significance=(permutation_fdr is not None),
        ),
    )
    match prev_res:
        case Success(annotated_df):
            return annotated_df
        case Failure(err):
            return res_df
