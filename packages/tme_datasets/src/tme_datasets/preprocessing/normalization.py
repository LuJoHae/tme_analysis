"""Pure functional library size normalization and linear transformations."""

from __future__ import annotations

import anndata as ad
import numpy as np
import scipy.sparse as sp
from returns.result import Failure, Result, Success


def normalize_total_counts(
    adata: ad.AnnData,
    target_sum: float = 1e6,
) -> Result[ad.AnnData, str]:
    """Normalize cell/sample library sizes to a fixed target sum (e.g. CPM = 10^6)."""
    try:
        new_adata = adata.copy()
        X = new_adata.X
        is_sparse = sp.issparse(X)

        counts_per_cell = np.asarray(X.sum(axis=1)).flatten()
        counts_per_cell[counts_per_cell == 0] = 1.0  # Avoid zero division
        scale_factors = target_sum / counts_per_cell

        if is_sparse:
            # Multiply each row by its scaling factor via diagonal matrix
            scaling_matrix = sp.diags(scale_factors)
            new_adata.X = scaling_matrix @ X
        else:
            new_adata.X = X * scale_factors[:, np.newaxis]

        return Success(new_adata)
    except Exception as exc:
        return Failure(f"Failed to normalize library sizes: {exc}")


def log1p_transform(adata: ad.AnnData) -> Result[ad.AnnData, str]:
    """Apply natural log(1 + x) transformation to AnnData expression values."""
    try:
        new_adata = adata.copy()
        X = new_adata.X
        if sp.issparse(X):
            new_adata.X = X.log1p()
        else:
            new_adata.X = np.log1p(X)
        return Success(new_adata)
    except Exception as exc:
        return Failure(f"Failed to apply log1p transform: {exc}")


def expm1_transform(adata: ad.AnnData) -> Result[ad.AnnData, str]:
    """Invert log1p transformation returning linear expression values: exp(x) - 1."""
    try:
        new_adata = adata.copy()
        X = new_adata.X
        if sp.issparse(X):
            new_adata.X = X.expm1()
        else:
            new_adata.X = np.expm1(X)
        return Success(new_adata)
    except Exception as exc:
        return Failure(f"Failed to apply expm1 transform: {exc}")
