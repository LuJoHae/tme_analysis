"""Auslander et al. Immuno-Predictive Score (IMPRES) 15-pair non-parametric classifier."""

from __future__ import annotations

import polars as pl
from returns.result import Failure, Result, Success

from ..schemas import MultiOmicCohort, PredictionResult
from ..types import PredictorCategory

# 15 checkpoint/co-stimulatory gene pairs from Auslander et al. (Nature Medicine 2018)
IMPRES_PAIRS = (
    ("PDCD1", "TNFRSF4"),
    ("CD27", "CD40"),
    ("CD27", "CD80"),
    ("CD40LG", "CD80"),
    ("CD40LG", "CD86"),
    ("CD40LG", "TNFRSF9"),
    ("CD40LG", "ICOSLG"),
    ("CD28", "CD86"),
    ("CD80", "HAVCR2"),
    ("CD86", "HAVCR2"),
    ("CD274", "HAVCR2"),
    ("CTLA4", "HAVCR2"),
    ("CD27", "HAVCR2"),
    ("ICOS", "HAVCR2"),
    ("TNFRSF4", "HAVCR2"),
)


def compute_impres_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Calculate the 15-pair non-parametric IMPRES score.

    Score is the count of pairs (A, B) where Expr(A) > Expr(B), yielding an integer in [0, 15].
    """
    try:
        expr = cohort.expression_tpm
        cols = set(expr.columns)

        valid_pair_exprs = []
        for idx, (gene_a, gene_b) in enumerate(IMPRES_PAIRS):
            if gene_a in cols and gene_b in cols:
                # 1 if gene_a > gene_b else 0
                valid_pair_exprs.append(
                    (pl.col(gene_a) > pl.col(gene_b)).cast(pl.Int32)
                )

        if not valid_pair_exprs:
            return Failure(f"None of the 15 IMPRES gene pairs were found in cohort {cohort.cohort_id}.")

        # Sum valid pair comparisons
        sum_expr = valid_pair_exprs[0]
        for p_expr in valid_pair_exprs[1:]:
            sum_expr = sum_expr + p_expr

        res_df = expr.select([
            pl.col("sample_id"),
            sum_expr.cast(pl.Float64).alias("score"),
        ])

        return Success(
            PredictionResult(
                predictor_name="IMPRES",
                category=PredictorCategory.SIGNATURE,
                predictions=res_df,
            )
        )
    except Exception as exc:
        return Failure(f"Failed to compute IMPRES score for {cohort.cohort_id}: {exc}")
