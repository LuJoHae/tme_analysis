"""Multi-omic DNA + RNA synergy models integrating TMB with transcriptomics."""

from __future__ import annotations

import polars as pl
from returns.maybe import Nothing, Some
from returns.result import Failure, Result, Success

from ..schemas import MultiOmicCohort, PredictionResult
from ..types import PredictorCategory


def compute_dna_rna_composite(
    cohort: MultiOmicCohort,
    rna_prediction: PredictionResult,
    rna_weight: float = 1.0,
    tmb_weight: float = 0.5,
) -> Result[PredictionResult, str]:
    """Compute composite DNA-RNA response index.

    Formula:
        Composite = rna_weight * z(RNA_Score) + tmb_weight * z(log10(TMB + 1))
    """
    match cohort.tmb_scores:
        case Some(tmb_df):
            try:
                # Join RNA prediction with TMB
                joined = rna_prediction.predictions.join(
                    tmb_df.select([pl.col("sample_id"), pl.col("tmb_per_mb")]),
                    on="sample_id",
                    how="inner",
                )

                if joined.is_empty():
                    return Failure(f"No intersecting samples between RNA and TMB in {cohort.cohort_id}.")

                # Standardize RNA score
                rna_col = pl.col("score")
                z_rna = (rna_col - rna_col.mean()) / rna_col.std(ddof=0).fill_nan(1.0)

                # Standardize log10(TMB + 1)
                log_tmb = (pl.col("tmb_per_mb") + 1.0).log10()
                z_tmb = (log_tmb - log_tmb.mean()) / log_tmb.std(ddof=0).fill_nan(1.0)

                composite_expr = (rna_weight * z_rna) + (tmb_weight * z_tmb)

                res_df = joined.select([
                    pl.col("sample_id"),
                    composite_expr.alias("score"),
                    z_rna.alias("rna_zscore"),
                    z_tmb.alias("tmb_zscore"),
                ])

                pred_name = f"Composite_{rna_prediction.predictor_name}_TMB"
                return Success(
                    PredictionResult(
                        predictor_name=pred_name,
                        category=PredictorCategory.COMPOSITE_SYNERGY,
                        predictions=res_df,
                    )
                )
            except Exception as exc:
                return Failure(f"Failed to compute DNA-RNA composite score: {exc}")
        case _:
            return Failure(f"No TMB genomic scores available in cohort {cohort.cohort_id}.")
