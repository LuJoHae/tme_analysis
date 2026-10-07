"""Ayers et al. Expanded T-cell-inflamed Gene Expression Profile (GEP) score."""

from __future__ import annotations

import numpy as np
import polars as pl
from returns.result import Failure, Result, Success

from tme_datasets.genesets.collections import AYERS_T_CELL_INFLAMED_GEP

from ..schemas import MultiOmicCohort, PredictionResult
from ..types import PredictorCategory


def compute_gep_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Calculate the Ayers 18-gene expanded T-cell inflamed GEP score.

    Computes standardized mean z-score across available signature genes:
    CCL5, CD27, CD274, CD276, CD8A, CMKLR1, CXCL9, CXCR6, HLA-DQA1, HLA-DRB1,
    HLA-E, IDO1, LAG3, NKG7, PDCD1LG2, PSMB9, STAT1, TIGIT.
    """
    try:
        expr = cohort.expression_tpm
        cols = set(expr.columns)
        available_genes = [g for g in AYERS_T_CELL_INFLAMED_GEP.genes if g in cols]

        if len(available_genes) < 5:
            return Failure(
                f"Too few Ayers GEP genes available in {cohort.cohort_id} "
                f"({len(available_genes)}/{len(AYERS_T_CELL_INFLAMED_GEP.genes)})."
            )

        # Standardize each gene across samples: (log2(x + 1) - mean) / std
        z_exprs = []
        for g in available_genes:
            log_col = (pl.col(g) + 1.0).log(base=2)
            z_col = (log_col - log_col.mean()) / log_col.std(ddof=0)
            z_exprs.append(z_col.fill_nan(0.0))

        # Compute average z-score across genes
        sum_z = z_exprs[0]
        for z in z_exprs[1:]:
            sum_z = sum_z + z
        mean_z = sum_z / float(len(available_genes))

        res_df = expr.select([
            pl.col("sample_id"),
            mean_z.alias("score"),
        ])

        return Success(
            PredictionResult(
                predictor_name="Ayers_GEP",
                category=PredictorCategory.SIGNATURE,
                predictions=res_df,
            )
        )
    except Exception as exc:
        return Failure(f"Failed to compute Ayers GEP score for {cohort.cohort_id}: {exc}")
