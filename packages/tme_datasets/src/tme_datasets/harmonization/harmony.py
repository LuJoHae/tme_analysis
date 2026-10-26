"""Harmony batch integration for multi-cohort and multi-donor single-cell AnnData objects."""

from __future__ import annotations

from typing import Sequence
import anndata as ad
import numpy as np
from returns.result import Failure, Result, Success
import scanpy as sc
from scipy import sparse as sp

try:
    import harmonypy as hm

    HAS_HARMONY = True
except ImportError:
    HAS_HARMONY = False


def run_harmony_batch_integration(
    adata: ad.AnnData,
    batch_key: str | Sequence[str] = "dataset_id",
    n_comps: int = 30,
    max_iter: int = 15,
    seed: int = 42,
) -> Result[ad.AnnData, str]:
    """Applies Harmony latent space batch correction across datasets or biological donors.

    Args:
        adata: Single-cell AnnData object (must contain >= 3 cells).
        batch_key: Column name or sequence of column names in .obs defining batch partitions.
        n_comps: Number of principal components used for latent representation.
        max_iter: Maximum number of Harmony convergence iterations.
        seed: Random seed for reproducible clustering and alignment.

    Returns:
        Success(harmonized_adata) with .obsm['X_pca_raw'] and .obsm['X_pca_harmony'],
        or Failure(error_message).
    """
    if not HAS_HARMONY:
        return Failure("harmonypy is not installed. Please install harmonypy to use Harmony batch integration.")

    if adata.n_obs < 3:
        return Failure(f"Cannot run Harmony integration on fewer than 3 cells (received {adata.n_obs})")

    # Validate batch key(s) exist in obs
    batch_keys = [batch_key] if isinstance(batch_key, str) else list(batch_key)
    for k in batch_keys:
        if k not in adata.obs.columns:
            return Failure(f"Batch key '{k}' not found in adata.obs columns: {list(adata.obs.columns)}")

    # Check that there are at least 2 distinct batches
    total_batches = 1
    for k in batch_keys:
        total_batches = max(total_batches, int(adata.obs[k].nunique()))
    if total_batches < 2:
        new_adata = adata.copy()
        if "X_pca" in new_adata.obsm:
            new_adata.obsm["X_pca_raw"] = new_adata.obsm["X_pca"].copy()
            new_adata.obsm["X_pca_harmony"] = new_adata.obsm["X_pca"].copy()
        return Success(new_adata)

    try:
        new_adata = adata.copy()

        # 1. Ensure PCA coordinates exist
        if "X_pca" not in new_adata.obsm:
            # Check if counts need log transformation
            max_val = (
                new_adata.X.max()
                if not sp.issparse(new_adata.X)
                else new_adata.X.data.max()
                if new_adata.X.nnz > 0
                else 0
            )
            if max_val > 50.0:
                sc.pp.normalize_total(new_adata, target_sum=1e4)
                sc.pp.log1p(new_adata)

            n_pcs = min(n_comps, new_adata.n_vars - 1, new_adata.n_obs - 1)
            n_pcs = max(2, n_pcs)

            if new_adata.n_vars > 2000:
                sc.pp.highly_variable_genes(new_adata, n_top_genes=2000, subset=False)
                has_hvg = "highly_variable" in new_adata.var.columns
            else:
                has_hvg = False

            sc.pp.pca(
                new_adata,
                n_comps=n_pcs,
                mask_var="highly_variable" if has_hvg else None,
                zero_center=False,
                random_state=seed,
            )

        # 2. Run harmonypy
        ho = hm.run_harmony(
            new_adata.obsm["X_pca"],
            new_adata.obs,
            batch_key,
            max_iter_harmony=max_iter,
            random_state=seed,
            verbose=False,
        )

        # 3. Store harmonized coordinates
        new_adata.obsm["X_pca_raw"] = new_adata.obsm["X_pca"].copy()
        new_adata.obsm["X_pca_harmony"] = ho.Z_corr
        new_adata.obsm["X_pca"] = ho.Z_corr

        return Success(new_adata)
    except Exception as exc:
        return Failure(f"Harmony integration failed: {exc}")
