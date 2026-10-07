"""Pure functional calculations for cohort predictability index and biomarker consensus."""

from __future__ import annotations

from typing import Sequence
import numpy as np
import polars as pl
from returns.result import Failure, Result, Success

from ..schemas import CohortPredictabilityResult

UNIVERSAL_RNA_PREDICTORS = (
    "CYT",
    "GEP",
    "IMPRES",
    "IPRES",
    "CXCL9",
    "CD8A",
    "TIDE",
    "MCP_CD8_T_cells",
)


def compute_cohort_predictability(
    benchmark_df: pl.DataFrame,
) -> Result[pl.DataFrame, str]:
    """Calculate the Cohort Predictability Index across predictors for each cohort and stratum.

    Separates the universal 8-biomarker RNA panel from all available predictors (including TMB/composites).
    Returns a Polars DataFrame sorted by mean RNA ROC-AUC in descending order.
    """
    if benchmark_df.is_empty():
        return Failure("Cannot compute cohort predictability from empty benchmark DataFrame.")

    required_cols = {
        "cohort_id",
        "cancer_type",
        "time_stratum",
        "response_stratum",
        "predictor_name",
        "roc_auc",
        "pr_auc",
        "delta_pr_auc",
        "baseline_prevalence",
        "n_samples",
        "n_responders",
    }
    missing = required_cols - set(benchmark_df.columns)
    if missing:
        return Failure(f"Missing required columns for predictability calculation: {sorted(missing)}")

    # Ensure pooling_strategy exists, defaulting to 'cohort' if absent
    df = (
        benchmark_df
        if "pooling_strategy" in benchmark_df.columns
        else benchmark_df.with_columns(pl.lit("cohort").alias("pooling_strategy"))
    )

    group_keys = ["cohort_id", "cancer_type", "time_stratum", "response_stratum", "pooling_strategy"]
    records: list[dict[str, object]] = []

    for group_vals, sub_df in df.group_by(group_keys):
        c_id, c_type, t_strat, r_strat, p_strat = group_vals

        n_samples = int(sub_df["n_samples"][0])
        n_resp = int(sub_df["n_responders"][0])
        prev = float(sub_df["baseline_prevalence"][0])

        # 1. Universal RNA predictors
        has_category = "category" in sub_df.columns
        pred_filter = (
            pl.col("predictor_name").is_in(list(UNIVERSAL_RNA_PREDICTORS))
            | pl.col("category").is_in(["signature", "cell_deconvolution", "systems_model"])
            if has_category
            else pl.col("predictor_name").is_in(list(UNIVERSAL_RNA_PREDICTORS))
        )
        rna_sub = sub_df.filter(pred_filter)
        if not rna_sub.is_empty():
            rna_aucs = rna_sub["roc_auc"].to_numpy()
            rna_prs = rna_sub["pr_auc"].to_numpy()
            rna_deltas = rna_sub["delta_pr_auc"].to_numpy()

            mean_auc_rna = float(np.mean(rna_aucs))
            median_auc_rna = float(np.median(rna_aucs))
            std_auc_rna = float(np.std(rna_aucs, ddof=1)) if len(rna_aucs) > 1 else 0.0
            mean_pr_rna = float(np.mean(rna_prs))
            mean_delta_rna = float(np.mean(rna_deltas))

            best_rna_idx = int(np.argmax(rna_aucs))
            best_pred_rna = str(rna_sub["predictor_name"][best_rna_idx])
            max_auc_rna = float(rna_aucs[best_rna_idx])
        else:
            mean_auc_rna = float(np.mean(sub_df["roc_auc"].to_numpy()))
            median_auc_rna = mean_auc_rna
            std_auc_rna = 0.0
            mean_pr_rna = float(np.mean(sub_df["pr_auc"].to_numpy()))
            mean_delta_rna = float(np.mean(sub_df["delta_pr_auc"].to_numpy()))
            best_pred_rna = "None"
            max_auc_rna = mean_auc_rna

        # 2. All available predictors (including TMB, composite, gating)
        all_aucs = sub_df["roc_auc"].to_numpy()
        mean_auc_all = float(np.mean(all_aucs))
        best_all_idx = int(np.argmax(all_aucs))
        best_pred_all = str(sub_df["predictor_name"][best_all_idx])
        max_auc_all = float(all_aucs[best_all_idx])
        n_eval = len(sub_df)

        records.append({
            "cohort_id": c_id,
            "cancer_type": c_type,
            "time_stratum": t_strat,
            "response_stratum": r_strat,
            "pooling_strategy": p_strat,
            "n_samples": n_samples,
            "n_responders": n_resp,
            "baseline_prevalence": prev,
            "mean_roc_auc_rna": round(mean_auc_rna, 4),
            "median_roc_auc_rna": round(median_auc_rna, 4),
            "std_roc_auc_rna": round(std_auc_rna, 4),
            "mean_pr_auc_rna": round(mean_pr_rna, 4),
            "mean_delta_pr_auc_rna": round(mean_delta_rna, 4),
            "best_predictor_rna": best_pred_rna,
            "max_roc_auc_rna": round(max_auc_rna, 4),
            "mean_roc_auc_all": round(mean_auc_all, 4),
            "best_predictor_all": best_pred_all,
            "max_roc_auc_all": round(max_auc_all, 4),
            "n_predictors_evaluated": n_eval,
        })

    predictability_df = pl.DataFrame(records).sort(
        by=["pooling_strategy", "mean_roc_auc_rna"], descending=[False, True]
    )

    return Success(predictability_df)
