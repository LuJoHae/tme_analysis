"""Microenvironment deconvolution wrappers and marker scoring."""

from __future__ import annotations

import polars as pl
from returns.result import Failure, Result, Success

from ..schemas import MultiOmicCohort, PredictionResult
from ..types import PredictorCategory

# Curated MCP-counter marker sets for primary cell lineages
MCP_MARKERS = {
    "CD8_T_cells": ("CD8A", "CD8B"),
    "Cytotoxic_lymphocytes": ("GZMB", "PRF1", "GNLY", "NKG7", "KLRD1"),
    "B_lineage": ("CD19", "MS4A1", "CD79A"),
    "Monocytic_lineage": ("CD14", "CD68", "CSF1R"),
    "Fibroblasts": ("COL1A1", "COL1A2", "DCN", "ACTA2", "FAP"),
    "Endothelial": ("PECAM1", "VWF", "CDH5"),
}


def compute_mcp_lineage_scores(cohort: MultiOmicCohort) -> Result[pl.DataFrame, str]:
    """Calculate geometric mean abundance scores for major immune and stromal populations."""
    try:
        expr = cohort.expression_tpm
        cols = set(expr.columns) - {"sample_id"}

        score_exprs = [pl.col("sample_id")]

        for pop_name, markers in MCP_MARKERS.items():
            avail = [m for m in markers if m in cols]
            if not avail:
                continue

            # Average log2(TPM + 1)
            pop_terms = [(pl.col(m) + 1.0).log(base=2) for m in avail]
            sum_terms = pop_terms[0]
            for t in pop_terms[1:]:
                sum_terms = sum_terms + t
            mean_term = sum_terms / float(len(avail))
            score_exprs.append(mean_term.alias(f"mcp_{pop_name}"))

        res_df = expr.select(score_exprs)
        return Success(res_df)
    except Exception as exc:
        return Failure(f"Failed to compute MCP-counter scores for {cohort.cohort_id}: {exc}")


def compute_cd8_infiltrate_predictor(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Use CD8+ T-cell infiltration estimate directly as an ICB response predictor."""
    match compute_mcp_lineage_scores(cohort):
        case Failure(err):
            return Failure(err)
        case Success(df):
            if "mcp_CD8_T_cells" not in df.columns:
                return Failure("CD8 T-cell markers not present in cohort.")

            pred_df = df.select([
                pl.col("sample_id"),
                pl.col("mcp_CD8_T_cells").alias("score"),
            ])
            return Success(
                PredictionResult(
                    predictor_name="MCP_CD8_T_cells",
                    category=PredictorCategory.CELL_DECONVOLUTION,
                    predictions=pred_df,
                )
            )
