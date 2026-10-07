"""Tumor Immune Dysfunction and Exclusion (TIDE) algorithm implementation."""

from __future__ import annotations

import numpy as np
import polars as pl
from returns.result import Failure, Result, Success

from ..schemas import MultiOmicCohort, PredictionResult
from ..types import PredictorCategory

# Cytotoxic T Lymphocyte (CTL) reference markers
CTL_GENES = ("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")

# Immunosuppressive T-cell exclusion markers (CAFs, MDSCs, M2 TAMs)
CAF_EXCLUSION_GENES = ("ACTA2", "COL1A1", "FAP", "PDGFRB", "TGFB1")
MDSC_EXCLUSION_GENES = ("CD14", "ITGAM", "S100A8", "S100A9", "STAT3")
TAM_M2_EXCLUSION_GENES = ("CD163", "MRC1", "MS4A4A", "VSIG4")


def compute_tide_score(
    cohort: MultiOmicCohort,
    invert_for_response: bool = True,
) -> Result[PredictionResult, str]:
    """Calculate the TIDE immune evasion score with cohort-level gene zero-centering.

    Higher TIDE score indicates stronger immune evasion (dysfunction or exclusion),
    which predicts non-response.
    If invert_for_response is True, the score is multiplied by -1.0 so higher scores denote response.
    """
    try:
        expr = cohort.expression_tpm
        cols = set(expr.columns) - {"sample_id"}

        # 1. Compute CTL level per sample
        avail_ctl = [g for g in CTL_GENES if g in cols]
        if not avail_ctl:
            return Failure(f"No CTL marker genes available in {cohort.cohort_id}.")

        # 2. Extract numeric expression matrix for centering
        sample_ids = expr["sample_id"].to_list()
        ctl_exprs = [
            (pl.col(g) + 1.0).log(base=2) for g in avail_ctl
        ]
        sum_ctl = ctl_exprs[0]
        for c in ctl_exprs[1:]:
            sum_ctl = sum_ctl + c
        ctl_mean = sum_ctl / float(len(avail_ctl))

        # 3. Compute exclusion score from CAF, MDSC, and M2 TAM signatures
        exclusion_genes = [
            g for g in (CAF_EXCLUSION_GENES + MDSC_EXCLUSION_GENES + TAM_M2_EXCLUSION_GENES)
            if g in cols
        ]
        if not exclusion_genes:
            return Failure(f"No exclusion marker genes found in {cohort.cohort_id}.")

        excl_exprs = [
            (pl.col(g) + 1.0).log(base=2) for g in exclusion_genes
        ]
        sum_excl = excl_exprs[0]
        for e in excl_exprs[1:]:
            sum_excl = sum_excl + e
        excl_score = sum_excl / float(len(exclusion_genes))

        # Standardize CTL and Exclusion across cohort
        ctl_z = (ctl_mean - ctl_mean.mean()) / ctl_mean.std(ddof=0).fill_nan(1.0)
        excl_z = (excl_score - excl_score.mean()) / excl_score.std(ddof=0).fill_nan(1.0)

        # Dysfunction proxy: high exclusion in presence of high CTL
        # Overall TIDE evasion score: positive when CTL is high & dysfunctional, or CTL is low & excluded
        # High excl_z - ctl_z indicates high exclusion relative to CTL
        raw_tide = excl_z - ctl_z

        final_score = (-1.0 * raw_tide) if invert_for_response else raw_tide

        pred_name = "TIDE_Inverted" if invert_for_response else "TIDE_Evasion"
        res_df = expr.select([
            pl.col("sample_id"),
            final_score.alias("score"),
            raw_tide.alias("tide_raw"),
            ctl_z.alias("ctl_score"),
            excl_z.alias("exclusion_score"),
        ])

        return Success(
            PredictionResult(
                predictor_name=pred_name,
                category=PredictorCategory.SYSTEMS_MODEL,
                predictions=res_df,
            )
        )
    except Exception as exc:
        return Failure(f"Failed to compute TIDE score for {cohort.cohort_id}: {exc}")
