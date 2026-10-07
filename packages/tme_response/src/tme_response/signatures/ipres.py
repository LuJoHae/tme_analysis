"""Hugo et al. Innate Anti-PD-1 Resistance (IPRES) signature scoring."""

from __future__ import annotations

import polars as pl
from returns.result import Failure, Result, Success

from ..schemas import MultiOmicCohort, PredictionResult
from ..types import PredictorCategory

# Core representative genes across IPRES mesenchymal, angiogenesis, and wound healing sets
IPRES_CORE_GENES = (
    "AXL", "ROR2", "WNT5A", "LOXL2", "TWIST2", "TAGLN", "FAP", "COL1A1",
    "COL3A1", "COL5A1", "VEGFA", "VEGFC", "ANGPT2", "PDGFRB", "FLT1", "KDR",
)


def compute_ipres_score(
    cohort: MultiOmicCohort,
    invert_for_response: bool = True,
) -> Result[PredictionResult, str]:
    """Calculate the Hugo IPRES resistance score.

    If invert_for_response is True, the score is multiplied by -1.0 so higher values
    consistently denote response / sensitivity across all benchmark metrics.
    """
    try:
        expr = cohort.expression_tpm
        cols = set(expr.columns)
        available = [g for g in IPRES_CORE_GENES if g in cols]

        if len(available) < 4:
            return Failure(f"Too few IPRES genes available in {cohort.cohort_id} ({len(available)}).")

        # Standardize genes
        z_exprs = []
        for g in available:
            log_col = (pl.col(g) + 1.0).log(base=2)
            z_col = (log_col - log_col.mean()) / log_col.std(ddof=0)
            z_exprs.append(z_col.fill_nan(0.0))

        sum_z = z_exprs[0]
        for z in z_exprs[1:]:
            sum_z = sum_z + z
        mean_z = sum_z / float(len(available))

        # Invert if response metric
        score_expr = (-1.0 * mean_z) if invert_for_response else mean_z

        res_df = expr.select([
            pl.col("sample_id"),
            score_expr.alias("score"),
        ])

        pred_name = "IPRES_Inverted" if invert_for_response else "IPRES_Resistance"
        return Success(
            PredictionResult(
                predictor_name=pred_name,
                category=PredictorCategory.SIGNATURE,
                predictions=res_df,
            )
        )
    except Exception as exc:
        return Failure(f"Failed to compute IPRES score for {cohort.cohort_id}: {exc}")
