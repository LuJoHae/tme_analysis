"""In silico pseudobulk benchmark validation for deconvolution and infiltrate estimators."""

from __future__ import annotations

import anndata as ad
import numpy as np
import polars as pl
from returns.result import Failure, Result, Success

from tme_datasets.simulation import simulate_pseudobulk


def benchmark_pseudobulk_deconvolution(
    adata_sc: ad.AnnData,
    sample_key: str = "sample_id",
    cell_type_key: str = "cell_type",
) -> Result[pl.DataFrame, str]:
    """Benchmark in silico deconvolution accuracy using simulated pseudobulk from single-cell data.

    Leverages tme_datasets.simulation.simulate_pseudobulk to generate mixtures with known cell fractions.
    """
    try:
        # Check metadata
        if sample_key not in adata_sc.obs.columns or cell_type_key not in adata_sc.obs.columns:
            return Failure(f"Single-cell AnnData missing '{sample_key}' or '{cell_type_key}' in .obs.")

        # Calculate true cell type fractions per sample
        obs_df = pl.from_pandas(adata_sc.obs)
        true_counts = (
            obs_df.group_by([sample_key, cell_type_key])
            .agg(pl.len().alias("count"))
        )
        sample_totals = (
            obs_df.group_by(sample_key)
            .agg(pl.len().alias("total"))
        )
        true_fractions = (
            true_counts.join(sample_totals, on=sample_key)
            .with_columns((pl.col("count") / pl.col("total")).alias("true_fraction"))
        )

        # Generate pseudobulk AnnData using tme_datasets simulation
        match simulate_pseudobulk(adata_sc, group_by=sample_key):
            case Failure(err):
                return Failure(f"Pseudobulk simulation failed: {err}")
            case Success(pb_adata):
                # Calculate CD8A / PRF1 surrogate score for CD8 T cell infiltrate
                var_names = list(pb_adata.var_names)
                var_map = {name: idx for idx, name in enumerate(var_names)}

                if "CD8A" not in var_map:
                    return Failure("CD8A not found in simulated pseudobulk.")

                X_dense = pb_adata.X.toarray() if hasattr(pb_adata.X, "toarray") else np.asarray(pb_adata.X)
                cd8_expr = np.log2(X_dense[:, var_map["CD8A"]] + 1.0)

                pred_df = pl.DataFrame({
                    sample_key: [str(x) for x in pb_adata.obs_names],
                    "est_cd8_score": [float(v) for v in cd8_expr],
                })

                # Join with true CD8+ T cell fractions
                cd8_true = (
                    true_fractions.filter(
                        pl.col(cell_type_key).str.to_lowercase().str.contains("cd8|cytotoxic|t_cell")
                    )
                    .group_by(sample_key)
                    .agg(pl.col("true_fraction").sum().alias("true_cd8_fraction"))
                )

                eval_df = pred_df.join(cd8_true, on=sample_key, how="inner")
                return Success(eval_df)
    except Exception as exc:
        return Failure(f"Pseudobulk benchmark failed: {exc}")
