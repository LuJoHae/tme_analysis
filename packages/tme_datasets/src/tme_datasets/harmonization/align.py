"""Gene feature space alignment: intersection and union zero-filling."""

from __future__ import annotations

from typing import Sequence
import anndata as ad
import numpy as np
import scipy.sparse as sp
from returns.result import Failure, Result, Success

from ..models import HarmonizeConfig
from ..types import HarmonizeMode


def align_and_concatenate(
    adatas: Sequence[ad.AnnData],
    dataset_ids: Sequence[str],
    config: HarmonizeConfig,
) -> Result[ad.AnnData, str]:
    """Harmonize multiple AnnData objects across common or union feature spaces and concatenate."""
    if not adatas or len(adatas) != len(dataset_ids):
        return Failure("List of AnnData objects and dataset IDs must be non-empty and equal length")

    if len(adatas) == 1:
        single = adatas[0].copy()
        single.obs[config.batch_key] = dataset_ids[0]
        return Success(single)

    try:
        gene_sets = [set(a.var_names) for a in adatas]

        match config.mode:
            case HarmonizeMode.INTERSECTION:
                common_genes = sorted(list(set.intersection(*gene_sets)))
                if len(common_genes) < config.min_shared_genes:
                    return Failure(
                        f"Shared genes ({len(common_genes)}) below minimum threshold ({config.min_shared_genes})"
                    )

                subsets = []
                for a, ds_id in zip(adatas, dataset_ids):
                    sub = a[:, common_genes].copy()
                    sub.obs[config.batch_key] = ds_id
                    subsets.append(sub)

                combined = ad.concat(subsets, axis=0, join="inner", merge="first")
                combined.obs_names_make_unique()
                return Success(combined)

            case HarmonizeMode.UNION_ZERO_FILLED:
                all_genes = sorted(list(set.union(*gene_sets)))
                gene_to_idx = {g: i for i, g in enumerate(all_genes)}
                n_all_genes = len(all_genes)

                aligned_matrices = []
                all_obs = []

                for a, ds_id in zip(adatas, dataset_ids):
                    n_cells = a.n_obs
                    X_orig = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)

                    # Map columns from old index to new unified index
                    col_map = np.array([gene_to_idx[g] for g in a.var_names])
                    row_indices, col_indices = X_orig.nonzero()
                    data = np.asarray(X_orig[row_indices, col_indices]).flatten()
                    new_cols = col_map[col_indices]

                    new_X = sp.csr_matrix(
                        (data, (row_indices, new_cols)),
                        shape=(n_cells, n_all_genes),
                        dtype=np.float32,
                    )
                    aligned_matrices.append(new_X)

                    obs_copy = a.obs.copy()
                    obs_copy[config.batch_key] = ds_id
                    all_obs.append(obs_copy)

                combined_X = sp.vstack(aligned_matrices, format="csr")
                import pandas as pd
                combined_obs = pd.concat(all_obs, axis=0)
                combined_obs.index = [f"cell_union_{i:06d}" for i in range(len(combined_obs))]

                combined = ad.AnnData(
                    X=combined_X,
                    obs=combined_obs,
                    var=pd.DataFrame(index=all_genes),
                )
                return Success(combined)
    except Exception as exc:
        return Failure(f"Failed to harmonize and concatenate datasets: {exc}")
