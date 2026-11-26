"""Pure functional sample filtering and metadata harmonization."""

from __future__ import annotations

import polars as pl
from returns.result import Failure, Result, Success

from ..types import BiopsyTimepoint


def filter_baseline_samples(
    clinical_df: pl.DataFrame,
    allow_unknown_timepoint: bool = True,
) -> Result[pl.DataFrame, str]:
    """Filter clinical annotations to keep only pre-treatment samples with known response.

    Pure function returning a Result containing the filtered Polars DataFrame.
    """
    try:
        # Check required columns
        required_cols = {"sample_id", "response_binary"}
        missing = required_cols - set(clinical_df.columns)
        if missing:
            return Failure(f"Missing required columns in clinical annotations: {missing}")

        allowed_timepoints = [BiopsyTimepoint.PRE.value]
        if allow_unknown_timepoint:
            allowed_timepoints.append(BiopsyTimepoint.UNKNOWN.value)

        # Check if biopsy_timepoint exists, otherwise assume Pre/Unknown
        has_tp = "biopsy_timepoint" in clinical_df.columns
        filtered = (
            clinical_df.filter(
                pl.col("response_binary").is_not_null()
                & (
                    pl.col("biopsy_timepoint").is_in(allowed_timepoints)
                    if has_tp
                    else pl.lit(True)
                )
            )
        )

        if filtered.is_empty():
            return Failure("No baseline pre-treatment samples with valid response were found.")

        return Success(filtered)
    except Exception as exc:
        return Failure(f"Error filtering baseline samples: {exc}")
