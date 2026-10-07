"""Rooney et al. Cytolytic Activity (CYT) score calculation."""

from __future__ import annotations

import numpy as np
import polars as pl
from returns.result import Failure, Result, Success

from ..schemas import MultiOmicCohort, PredictionResult
from ..types import PredictorCategory


def compute_cyt_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Calculate the Rooney Cytolytic Activity score (geometric mean of GZMA and PRF1).

    Formula:
        CYT = (log2(GZMA + 1) + log2(PRF1 + 1)) / 2.0
    """
    try:
        expr = cohort.expression_tpm
        cols = set(expr.columns)

        if "GZMA" not in cols or "PRF1" not in cols:
            return Failure(f"Cohort {cohort.cohort_id} missing GZMA or PRF1 in expression data.")

        # Vectorized calculation in Polars
        res_df = expr.select([
            pl.col("sample_id"),
            (
                ((pl.col("GZMA") + 1.0).log(base=2) + (pl.col("PRF1") + 1.0).log(base=2))
                / 2.0
            ).alias("score")
        ])

        return Success(
            PredictionResult(
                predictor_name="CYT",
                category=PredictorCategory.SIGNATURE,
                predictions=res_df,
            )
        )
    except Exception as exc:
        return Failure(f"Failed to compute CYT score for {cohort.cohort_id}: {exc}")
