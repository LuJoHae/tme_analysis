"""Pure functional subsampling and supersampling of transcriptomic datasets."""

from __future__ import annotations

import anndata as ad
import numpy as np
from returns.maybe import Maybe, Some
from returns.result import Failure, Result, Success

from ..models import SubsampleSpec


def subsample_cells(
    adata: ad.AnnData,
    spec: SubsampleSpec,
) -> Result[ad.AnnData, str]:
    """Subsample cells uniformly or stratified across clinical/biological annotations."""
    try:
        rng = np.random.default_rng(spec.seed.value_or(None) if isinstance(spec.seed, Some) else None)
        n_total = adata.n_obs

        target_n = (
            int(n_total * spec.n_or_fraction)
            if spec.n_or_fraction <= 1.0
            else min(n_total, int(spec.n_or_fraction))
        )

        match spec.stratify_by:
            case Some(strata_key) if strata_key in adata.obs.columns:
                strata = adata.obs[strata_key].to_numpy()
                unique_strata = np.unique(strata)
                chosen_indices = []

                if spec.balanced:
                    # Equal cells per category
                    n_per_stratum = max(1, target_n // len(unique_strata))
                    for val in unique_strata:
                        idx = np.where(strata == val)[0]
                        k = min(len(idx), n_per_stratum)
                        chosen = rng.choice(idx, size=k, replace=False)
                        chosen_indices.extend(chosen)
                else:
                    # Proportional sampling per category
                    for val in unique_strata:
                        idx = np.where(strata == val)[0]
                        k = max(1, int(round(len(idx) / n_total * target_n)))
                        k = min(len(idx), k)
                        chosen = rng.choice(idx, size=k, replace=False)
                        chosen_indices.extend(chosen)

                selected = np.sort(chosen_indices)
            case _:
                selected = np.sort(rng.choice(n_total, size=target_n, replace=False))

        return Success(adata[selected].copy())
    except Exception as exc:
        return Failure(f"Failed to subsample cells: {exc}")


def supersample_cells(
    adata: ad.AnnData,
    n_target: int,
    stratify_by: Maybe[str] = None, # type: ignore
    seed: Maybe[int] = None, # type: ignore
) -> Result[ad.AnnData, str]:
    """Supersample/bootstrap cells with replacement to achieve target cell count."""
    try:
        rng_seed = seed.value_or(None) if isinstance(seed, Some) else None
        rng = np.random.default_rng(rng_seed)
        n_total = adata.n_obs

        if n_target <= n_total:
            return Success(adata.copy())

        selected_indices = rng.choice(n_total, size=n_target, replace=True)
        new_adata = adata[selected_indices].copy()
        new_adata.obs_names = [f"cell_super_{i:06d}" for i in range(n_target)]
        return Success(new_adata)
    except Exception as exc:
        return Failure(f"Failed to supersample cells: {exc}")
