"""In-silico targeted gene knockout and overexpression perturbations."""

from __future__ import annotations

import anndata as ad
import numpy as np
import scipy.sparse as sp
from returns.result import Failure, Result, Success


def in_silico_knockout(
    adata: ad.AnnData,
    genes: tuple[str, ...],
    efficiency: float = 1.0,
) -> Result[ad.AnnData, str]:
    """Simulate targeted gene inhibition or knockout by downscaling target genes by (1 - efficiency)."""
    if not 0.0 <= efficiency <= 1.0:
        return Failure(f"Knockout efficiency must be in [0, 1], got {efficiency}")

    try:
        new_adata = adata.copy()
        var_dict = {name: idx for idx, name in enumerate(adata.var_names)}
        target_indices = [var_dict[g] for g in genes if g in var_dict]

        if not target_indices:
            return Success(new_adata)

        X = new_adata.X.toarray() if sp.issparse(new_adata.X) else np.asarray(new_adata.X).copy()
        factor = 1.0 - efficiency
        X[:, target_indices] *= factor

        new_adata.X = sp.csr_matrix(X) if sp.issparse(adata.X) else X
        return Success(new_adata)
    except Exception as exc:
        return Failure(f"Failed to execute in-silico knockout: {exc}")


def in_silico_overexpression(
    adata: ad.AnnData,
    genes: tuple[str, ...],
    fold_change: float = 2.0,
) -> Result[ad.AnnData, str]:
    """Simulate targeted gene overexpression by multiplying target genes by fold_change."""
    if fold_change < 0.0:
        return Failure(f"Fold change must be non-negative, got {fold_change}")

    try:
        new_adata = adata.copy()
        var_dict = {name: idx for idx, name in enumerate(adata.var_names)}
        target_indices = [var_dict[g] for g in genes if g in var_dict]

        if not target_indices:
            return Success(new_adata)

        X = new_adata.X.toarray() if sp.issparse(new_adata.X) else np.asarray(new_adata.X).copy()
        X[:, target_indices] *= fold_change

        new_adata.X = sp.csr_matrix(X) if sp.issparse(adata.X) else X
        return Success(new_adata)
    except Exception as exc:
        return Failure(f"Failed to execute in-silico overexpression: {exc}")
