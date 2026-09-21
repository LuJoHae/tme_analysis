"""Pure reconciliation of gene features on AnnData objects."""

from __future__ import annotations

import anndata as ad
import numpy as np
import scipy.sparse as sp
from returns.maybe import Some
from returns.result import Failure, Result, Success

from ..models import GeneReconcileConfig
from .mapper import detect_gene_id_type, map_gene_identifier


def reconcile_genes(adata: ad.AnnData, config: GeneReconcileConfig) -> Result[ad.AnnData, str]:
    """Normalize and reconcile the feature names of an AnnData object."""
    try:
        current_names = tuple(str(g) for g in adata.var_names)
        source_type = detect_gene_id_type(current_names)

        if source_type == config.target_type:
            # Already matching target type, only strip version suffix if needed
            new_names = [
                g.split(".")[0] if config.strip_version_suffix and g.startswith("ENS") else g
                for g in current_names
            ]
        else:
            mapped = []
            for g in current_names:
                match map_gene_identifier(g, config.target_type, config.strip_version_suffix):
                    case Some(mapped_name):
                        mapped.append(mapped_name)
                    case _:
                        mapped.append(g)
            new_names = mapped

        # Handle potential duplicates in target names
        unique_names, inverse_indices, counts = np.unique(
            new_names, return_inverse=True, return_counts=True
        )

        if len(unique_names) == len(new_names):
            new_adata = adata.copy()
            new_adata.var_names = new_names
            return Success(new_adata)

        # Aggregate duplicates (default sum)
        X = adata.X
        is_sparse = sp.issparse(X)
        n_cells = adata.n_obs
        n_unique_genes = len(unique_names)

        # Construct projection aggregation matrix (n_genes x n_unique_genes)
        row_ind = np.arange(len(new_names))
        col_ind = inverse_indices
        data = np.ones(len(new_names), dtype=np.float32)
        proj = sp.csr_matrix((data, (row_ind, col_ind)), shape=(len(new_names), n_unique_genes))

        if is_sparse:
            new_X = X @ proj
        else:
            new_X = X @ proj.toarray()

        new_adata = ad.AnnData(
            X=new_X,
            obs=adata.obs.copy(),
            var=None,
        )
        new_adata.var_names = list(unique_names)
        return Success(new_adata)
    except Exception as exc:
        return Failure(f"Failed to reconcile genes on AnnData: {exc}")
