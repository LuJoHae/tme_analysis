"""Pure functional in-silico pseudobulk mixture simulation."""

from __future__ import annotations

import anndata as ad
import numpy as np
import polars as pl
import scipy.sparse as sp
from returns.maybe import Some
from returns.result import Failure, Result, Success

from ..models import PseudobulkConfig


def simulate_pseudobulk(
    adata: ad.AnnData,
    config: PseudobulkConfig,
    cell_type_key: str = "cell_type",
) -> Result[tuple[ad.AnnData, pl.DataFrame], str]:
    """Generate in-silico bulk mixtures with exact known ground-truth cell-type proportions.

    Returns:
        (bulk_adata, ground_truth_fractions_df)
    """
    if cell_type_key not in adata.obs.columns:
        return Failure(f"Cell type column '{cell_type_key}' not found in adata.obs")

    try:
        rng_seed = config.seed.value_or(None) if isinstance(config.seed, Some) else None
        rng = np.random.default_rng(rng_seed)

        cell_types = np.asarray(adata.obs[cell_type_key])
        unique_types = sorted(list(np.unique(cell_types)))
        n_types = len(unique_types)

        if n_types < 2:
            return Failure(f"At least 2 cell types required for mixtures, found {n_types}")

        # Cell indices partitioned by type
        type_indices = {ct: np.where(cell_types == ct)[0] for ct in unique_types}

        X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        n_samples = config.n_samples
        n_genes = adata.n_vars
        cells_per_sample = config.cells_per_sample

        bulk_matrix = np.zeros((n_samples, n_genes), dtype=np.float32)
        proportions_data: dict[str, list[float | str]] = {
            "sample_id": [f"Simulated_Bulk_{i:03d}" for i in range(n_samples)]
        }
        for ct in unique_types:
            proportions_data[ct] = []

        # Dirichlet alpha prior
        alpha_prior = np.ones(n_types, dtype=np.float32)

        for s_idx in range(n_samples):
            # Draw proportions from Dirichlet
            props = rng.dirichlet(alpha_prior)
            cell_counts = rng.multinomial(cells_per_sample, props)

            sample_cells = []
            for ct_idx, ct in enumerate(unique_types):
                n_draw = cell_counts[ct_idx]
                actual_fraction = float(n_draw / cells_per_sample)
                proportions_data[ct].append(actual_fraction)

                if n_draw > 0:
                    avail_idx = type_indices[ct]
                    # Sample with replacement if requested cells exceed available
                    replace = len(avail_idx) < n_draw
                    chosen = rng.choice(avail_idx, size=n_draw, replace=replace)
                    sample_cells.extend(chosen)

            # Sum expression of selected cells
            if sample_cells:
                bulk_matrix[s_idx, :] = np.sum(X[sample_cells, :], axis=0)

        # Optional Negative Binomial sequencing noise
        if isinstance(config.noise_dispersion, Some):
            disp = max(1e-5, config.noise_dispersion.value_or(0.1))
            shape = 1.0 / disp
            scale = disp * np.maximum(bulk_matrix, 0.0)
            lam = rng.gamma(shape=shape, scale=np.maximum(scale, 1e-8))
            bulk_matrix = rng.poisson(lam).astype(np.float32)

        sample_names = [f"Simulated_Bulk_{i:03d}" for i in range(n_samples)]
        bulk_obs = pl.DataFrame({"sample_id": sample_names}).to_pandas().set_index("sample_id")

        bulk_adata = ad.AnnData(
            X=bulk_matrix,
            obs=bulk_obs,
            var=adata.var.copy(),
        )

        truth_df = pl.DataFrame(proportions_data)
        return Success((bulk_adata, truth_df))
    except Exception as exc:
        return Failure(f"Failed to simulate pseudobulk mixtures: {exc}")
