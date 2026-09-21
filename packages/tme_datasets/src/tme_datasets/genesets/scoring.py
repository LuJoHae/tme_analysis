"""Pure functional gene signature scoring algorithms for transcriptomic data."""

from __future__ import annotations

import anndata as ad
import numpy as np
import polars as pl
import scipy.sparse as sp
from returns.result import Failure, Result, Success

from .models import GeneSetCollection


def score_geneset_zscore(
    adata: ad.AnnData,
    collection: GeneSetCollection,
) -> Result[pl.DataFrame, str]:
    """Calculate mean standardized Z-scores for each gene set across cells/samples."""
    try:
        var_names = list(adata.var_names)
        var_dict = {name: idx for idx, name in enumerate(var_names)}
        X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)

        # Standardize genes across cells (Z-score)
        gene_means = np.mean(X, axis=0)
        gene_stds = np.std(X, axis=0)
        gene_stds[gene_stds == 0] = 1.0
        X_z = (X - gene_means) / gene_stds

        scores_dict: dict[str, list[float]] = {
            "sample_id": [str(x) for x in adata.obs_names]
        }

        for set_id, gs in collection.gene_sets.items():
            valid_indices = [var_dict[g] for g in gs.genes if g in var_dict]
            if not valid_indices:
                scores_dict[set_id] = [0.0] * adata.n_obs
            else:
                set_scores = np.mean(X_z[:, valid_indices], axis=1)
                scores_dict[set_id] = [float(v) for v in set_scores]

        return Success(pl.DataFrame(scores_dict))
    except Exception as exc:
        return Failure(f"Failed to compute Z-score signatures: {exc}")


def score_geneset_auc(
    adata: ad.AnnData,
    collection: GeneSetCollection,
    top_fraction: float = 0.20,
) -> Result[pl.DataFrame, str]:
    """Calculate rank-based AUCell enrichment scores for gene sets in each cell/sample."""
    try:
        var_names = list(adata.var_names)
        var_dict = {name: idx for idx, name in enumerate(var_names)}
        X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)

        n_cells, n_genes = X.shape
        k_threshold = max(1, int(n_genes * top_fraction))

        # Compute ranks per cell (descending: highest expression gets rank 1)
        ranks = np.argsort(np.argsort(-X, axis=1), axis=1) + 1

        scores_dict: dict[str, list[float]] = {
            "sample_id": [str(x) for x in adata.obs_names]
        }

        for set_id, gs in collection.gene_sets.items():
            valid_indices = [var_dict[g] for g in gs.genes if g in var_dict]
            if not valid_indices:
                scores_dict[set_id] = [0.0] * n_cells
                continue

            set_ranks = ranks[:, valid_indices]
            # AUC over ranks within top threshold
            recovered_in_top = np.sum(set_ranks <= k_threshold, axis=1) / len(valid_indices)
            scores_dict[set_id] = [float(v) for v in recovered_in_top]

        return Success(pl.DataFrame(scores_dict))
    except Exception as exc:
        return Failure(f"Failed to compute rank AUC signatures: {exc}")
