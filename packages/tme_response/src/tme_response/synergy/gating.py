"""Antigen presentation and immune evasion mutation gating logic."""

from __future__ import annotations

import polars as pl
from returns.maybe import Nothing, Some
from returns.result import Failure, Result, Success

from ..schemas import MultiOmicCohort, PredictionResult
from ..types import PredictorCategory


def apply_antigen_presentation_gating(
    cohort: MultiOmicCohort,
    base_prediction: PredictionResult,
    gating_genes: tuple[str, ...] = ("B2M", "JAK1", "JAK2"),
) -> Result[PredictionResult, str]:
    """Gate response score by somatic mutations in core antigen presentation genes.

    If a sample harbors an inactivating somatic mutation in B2M, JAK1, or JAK2,
    its response score is gated to the minimum cohort score, reflecting structural resistance.
    """
    match cohort.driver_mutations:
        case Some(drv_df):
            try:
                # Find samples with inactivating mutations in gating genes
                inactivated_samples = (
                    drv_df.filter(
                        pl.col("gene").is_in(list(gating_genes))
                        & (pl.col("is_inactivating") == True)
                    )
                    .select("sample_id")
                    .unique()["sample_id"]
                    .to_list()
                )

                preds = base_prediction.predictions
                min_score = float(preds["score"].min() or 0.0)

                # Gate: if sample is in inactivated_samples, set to min_score - 1.0
                gated_df = preds.with_columns(
                    pl.when(pl.col("sample_id").is_in(inactivated_samples))
                    .then(pl.lit(min_score - 1.0))
                    .otherwise(pl.col("score"))
                    .alias("score"),
                    pl.col("sample_id").is_in(inactivated_samples).alias("is_ap_gated"),
                )

                pred_name = f"{base_prediction.predictor_name}_AP_Gated"
                return Success(
                    PredictionResult(
                        predictor_name=pred_name,
                        category=PredictorCategory.COMPOSITE_SYNERGY,
                        predictions=gated_df,
                    )
                )
            except Exception as exc:
                return Failure(f"Failed to apply antigen presentation gating: {exc}")
        case _:
            # If no driver mutations available, return base prediction unmodified
            return Success(base_prediction)

