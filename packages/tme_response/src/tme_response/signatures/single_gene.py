"""Litchfield et al. Single-gene transcriptomic predictors (CXCL9, CD8A, PDCD1, IFNG)."""

from __future__ import annotations

import polars as pl
from returns.result import Failure, Result, Success

from ..schemas import MultiOmicCohort, PredictionResult
from ..types import PredictorCategory


def compute_single_gene_score(
    cohort: MultiOmicCohort,
    gene_name: str = "CXCL9",
) -> Result[PredictionResult, str]:
    """Calculate single-gene log2(TPM + 1) predictor score (e.g. CXCL9, CD8A)."""
    try:
        expr = cohort.expression_tpm
        if gene_name not in expr.columns:
            return Failure(f"Gene '{gene_name}' not found in cohort {cohort.cohort_id}.")

        res_df = expr.select([
            pl.col("sample_id"),
            (pl.col(gene_name) + 1.0).log(base=2).alias("score"),
        ])

        return Success(
            PredictionResult(
                predictor_name=f"Gene_{gene_name}",
                category=PredictorCategory.SIGNATURE,
                predictions=res_df,
            )
        )
    except Exception as exc:
        return Failure(f"Failed to extract single gene {gene_name} for {cohort.cohort_id}: {exc}")
