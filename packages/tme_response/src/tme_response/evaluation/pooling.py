"""Pure functional cross-cohort score standardization and pooling algorithms."""

from __future__ import annotations

from typing import Mapping, Sequence
import numpy as np
import polars as pl
from scipy import stats
from returns.result import Failure, Result, Success

from ..schemas import PredictionResult
from ..types import PredictorCategory


def standardize_prediction_scores(prediction: PredictionResult) -> PredictionResult:
    """Standardize prediction scores within a cohort to zero mean and unit variance (z-score).

    Preserves monotonic intra-cohort ranking while aligning cross-study baselines.
    If standard deviation is zero (invariant score), values are centered to 0.0.
    """
    scores = prediction.predictions["score"].to_numpy()
    std = float(np.std(scores))
    mean = float(np.mean(scores))

    if std > 1e-9:
        z_scores = (scores - mean) / std
    else:
        z_scores = np.zeros_like(scores)

    updated_df = prediction.predictions.with_columns(
        pl.Series("score", z_scores)
    )

    return PredictionResult(
        predictor_name=prediction.predictor_name,
        category=prediction.category,
        predictions=updated_df,
    )


def pool_cohort_stratum_data(
    cohort_strata: Sequence[tuple[str, pl.DataFrame, Sequence[PredictionResult]]],
    group_id: str,
    cancer_type: str,
    standardize: bool = True,
) -> Result[tuple[pl.DataFrame, tuple[PredictionResult, ...]], str]:
    """Pool clinical response annotations and prediction scores across multiple cohorts.

    Args:
        cohort_strata: Sequence of (cohort_id, clinical_df, list_of_predictions)
        group_id: Identifier for the combined cohort (e.g. 'Melanoma-Combined')
        cancer_type: Cancer type descriptor (e.g. 'Melanoma' or 'Pan-Cancer')
        standardize: If True, applies intra-cohort z-scoring before concatenation.
                     If False, concatenates raw prediction scores directly.

    Returns:
        tuple of (pooled_clinical_df, tuple_of_pooled_predictions)
    """
    if not cohort_strata:
        return Failure("No cohort strata provided for pooling.")

    pooled_clin_dfs: list[pl.DataFrame] = []
    # Map predictor_name -> (category, list_of_prediction_dfs)
    pred_collector: dict[str, tuple[PredictorCategory, list[pl.DataFrame]]] = {}

    for cohort_id, clin_df, preds in cohort_strata:
        if clin_df.is_empty():
            continue

        # Select key clinical columns and prefix sample IDs to prevent barcode collisions
        clin_cols = [c for c in ["sample_id", "response_binary", "response_recist", "biopsy_timepoint"] if c in clin_df.columns]
        prefixed_clin = clin_df.select(clin_cols).with_columns(
            pl.concat_str([pl.lit(f"{cohort_id}::"), pl.col("sample_id")]).alias("sample_id")
        )
        pooled_clin_dfs.append(prefixed_clin)

        for p in preds:
            proc_p = standardize_prediction_scores(p) if standardize else p
            pred_cols = [c for c in ["sample_id", "score"] if c in proc_p.predictions.columns]
            prefixed_preds = proc_p.predictions.select(pred_cols).with_columns(
                pl.concat_str([pl.lit(f"{cohort_id}::"), pl.col("sample_id")]).alias("sample_id")
            )
            if p.predictor_name not in pred_collector:
                pred_collector[p.predictor_name] = (p.category, [])
            pred_collector[p.predictor_name][1].append(prefixed_preds)

    if not pooled_clin_dfs:
        return Failure(f"All cohort strata were empty for group '{group_id}'.")

    pooled_clinical = pl.concat(pooled_clin_dfs, how="diagonal_relaxed")

    pooled_predictions: list[PredictionResult] = []
    for pred_name, (cat, df_list) in pred_collector.items():
        if df_list:
            concatenated = pl.concat(df_list, how="diagonal_relaxed")
            pooled_predictions.append(
                PredictionResult(
                    predictor_name=pred_name,
                    category=cat,
                    predictions=concatenated,
                )
            )

    return Success((pooled_clinical, tuple(pooled_predictions)))


def compute_meta_analytic_auc(
    aucs: Sequence[float],
    ci_lowers: Sequence[float],
    ci_uppers: Sequence[float],
) -> dict[str, float]:
    """Compute fixed-effect and DerSimonian-Laird random-effects meta-analytic summary AUC.

    Approximates standard error from 95% bootstrap confidence interval width:
    SE ≈ (CI_upper - CI_lower) / (2 * 1.96).
    """
    arr_auc = np.asarray(aucs, dtype=float)
    arr_lower = np.asarray(ci_lowers, dtype=float)
    arr_upper = np.asarray(ci_uppers, dtype=float)

    k = len(arr_auc)
    if k == 0:
        return {
            "meta_auc_fixed": 0.5,
            "meta_auc_random": 0.5,
            "tau_squared": 0.0,
            "i_squared": 0.0,
            "p_heterogeneity": 1.0,
        }
    if k == 1:
        return {
            "meta_auc_fixed": float(arr_auc[0]),
            "meta_auc_random": float(arr_auc[0]),
            "tau_squared": 0.0,
            "i_squared": 0.0,
            "p_heterogeneity": 1.0,
        }

    # Estimate standard error
    se = (arr_upper - arr_lower) / (2.0 * 1.95996)
    se = np.maximum(se, 1e-4)
    w_fe = 1.0 / (se**2)

    # Fixed effects summary
    sum_w = float(np.sum(w_fe))
    theta_fe = float(np.sum(w_fe * arr_auc) / sum_w)

    # Cochran's Q
    q_stat = float(np.sum(w_fe * (arr_auc - theta_fe) ** 2))
    df = k - 1
    p_het = float(1.0 - stats.chi2.cdf(q_stat, df)) if df > 0 else 1.0

    # DerSimonian-Laird tau^2
    sum_w_sq = float(np.sum(w_fe**2))
    c_val = sum_w - (sum_w_sq / sum_w)
    tau_sq = float(max(0.0, (q_stat - df) / c_val)) if c_val > 0 else 0.0

    # Random effects summary
    w_re = 1.0 / (se**2 + tau_sq)
    theta_re = float(np.sum(w_re * arr_auc) / np.sum(w_re))

    # Higgins I^2 (%)
    i_sq = float(max(0.0, (q_stat - df) / q_stat) * 100.0) if q_stat > df and q_stat > 0 else 0.0

    return {
        "meta_auc_fixed": theta_fe,
        "meta_auc_random": theta_re,
        "tau_squared": tau_sq,
        "i_squared": i_sq,
        "p_heterogeneity": p_het,
    }
