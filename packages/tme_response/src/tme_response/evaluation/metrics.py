"""Pure functional evaluation metrics for response prediction."""

from __future__ import annotations

from typing import Sequence
import numpy as np
import polars as pl
from returns.result import Failure, Result, Success
from sklearn.metrics import auc, brier_score_loss, precision_recall_curve, roc_auc_score

from ..schemas import CohortBenchmarkResult, PredictionResult


def calculate_roc_auc(y_true: Sequence[float], y_score: Sequence[float]) -> float:
    """Calculate ROC-AUC handling degenerate cases gracefully."""
    arr_true = np.asarray(y_true)
    arr_score = np.asarray(y_score)
    # Check if more than one class present
    if len(np.unique(arr_true[~np.isnan(arr_true)])) < 2:
        return 0.5
    try:
        return float(roc_auc_score(arr_true, arr_score))
    except Exception:
        return 0.5


def calculate_discrimination_bootstrap(
    y_true: Sequence[float],
    y_score: Sequence[float],
    n_bootstrap: int = 1000,
    alpha: float = 0.05,
    seed: int = 42,
) -> dict[str, float]:
    """Calculate point-estimates, 95% bootstrap confidence intervals, and p-values for ROC-AUC and PR-AUC.

    Returns dict containing:
        - roc_auc, roc_auc_ci_lower, roc_auc_ci_upper, p_value_vs_half
        - pr_auc, pr_auc_ci_lower, pr_auc_ci_upper, p_value_prauc
        - baseline_prevalence, delta_pr_auc
    """
    arr_true = np.asarray(y_true)
    arr_score = np.asarray(y_score)
    n = len(arr_true)
    n_resp = int(np.sum(arr_true == 1.0))
    prevalence = float(n_resp / n) if n > 0 else 0.5

    point_roc = calculate_roc_auc(arr_true, arr_score)
    point_pr = calculate_pr_auc(arr_true, arr_score)
    delta_pr = point_pr - prevalence

    if n_bootstrap <= 0 or len(np.unique(arr_true)) < 2:
        return {
            "roc_auc": point_roc,
            "roc_auc_ci_lower": point_roc,
            "roc_auc_ci_upper": point_roc,
            "p_value_vs_half": 1.0,
            "pr_auc": point_pr,
            "pr_auc_ci_lower": point_pr,
            "pr_auc_ci_upper": point_pr,
            "p_value_prauc": 1.0,
            "baseline_prevalence": prevalence,
            "delta_pr_auc": delta_pr,
        }

    rng = np.random.default_rng(seed)
    boot_roc: list[float] = []
    boot_pr: list[float] = []

    for _ in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)
        b_true = arr_true[idx]
        b_score = arr_score[idx]
        if len(np.unique(b_true)) < 2:
            continue
        try:
            r_val = float(roc_auc_score(b_true, b_score))
            prec, rec, _ = precision_recall_curve(b_true, b_score)
            p_val = float(auc(rec, prec))
            boot_roc.append(r_val)
            boot_pr.append(p_val)
        except Exception:
            continue

    if not boot_roc or not boot_pr:
        return {
            "roc_auc": point_roc,
            "roc_auc_ci_lower": point_roc,
            "roc_auc_ci_upper": point_roc,
            "p_value_vs_half": 1.0,
            "pr_auc": point_pr,
            "pr_auc_ci_lower": point_pr,
            "pr_auc_ci_upper": point_pr,
            "p_value_prauc": 1.0,
            "baseline_prevalence": prevalence,
            "delta_pr_auc": delta_pr,
        }

    boot_roc_arr = np.asarray(boot_roc)
    boot_pr_arr = np.asarray(boot_pr)

    ci_lower_roc = float(np.percentile(boot_roc_arr, 100.0 * (alpha / 2.0)))
    ci_upper_roc = float(np.percentile(boot_roc_arr, 100.0 * (1.0 - alpha / 2.0)))

    ci_lower_pr = float(np.percentile(boot_pr_arr, 100.0 * (alpha / 2.0)))
    ci_upper_pr = float(np.percentile(boot_pr_arr, 100.0 * (1.0 - alpha / 2.0)))

    # Empirical two-sided p-value against chance (0.5 for ROC)
    p_le_roc = np.mean(boot_roc_arr <= 0.5)
    p_ge_roc = np.mean(boot_roc_arr >= 0.5)
    p_val_roc = float(min(1.0, 2.0 * min(p_le_roc, p_ge_roc)))

    # Empirical two-sided p-value against baseline prevalence for PR-AUC
    p_le_pr = np.mean(boot_pr_arr <= prevalence)
    p_ge_pr = np.mean(boot_pr_arr >= prevalence)
    p_val_pr = float(min(1.0, 2.0 * min(p_le_pr, p_ge_pr)))

    return {
        "roc_auc": point_roc,
        "roc_auc_ci_lower": ci_lower_roc,
        "roc_auc_ci_upper": ci_upper_roc,
        "p_value_vs_half": p_val_roc,
        "pr_auc": point_pr,
        "pr_auc_ci_lower": ci_lower_pr,
        "pr_auc_ci_upper": ci_upper_pr,
        "p_value_prauc": p_val_pr,
        "baseline_prevalence": prevalence,
        "delta_pr_auc": delta_pr,
    }


def calculate_roc_auc_bootstrap(
    y_true: Sequence[float],
    y_score: Sequence[float],
    n_bootstrap: int = 1000,
    alpha: float = 0.05,
    seed: int = 42,
) -> tuple[float, float, float, float]:
    """Calculate point-estimate ROC-AUC, 95% bootstrap confidence intervals, and p-value vs 0.5.

    Returns:
        tuple of (point_auc, ci_lower, ci_upper, p_value_vs_half)
    """
    res = calculate_discrimination_bootstrap(
        y_true, y_score, n_bootstrap=n_bootstrap, alpha=alpha, seed=seed
    )
    return (
        res["roc_auc"],
        res["roc_auc_ci_lower"],
        res["roc_auc_ci_upper"],
        res["p_value_vs_half"],
    )


def calculate_pr_auc(y_true: Sequence[float], y_score: Sequence[float]) -> float:
    """Calculate Precision-Recall AUC (Average Precision)."""
    arr_true = np.asarray(y_true)
    arr_score = np.asarray(y_score)
    if len(np.unique(arr_true[~np.isnan(arr_true)])) < 2:
        return 0.0
    try:
        prec, rec, _ = precision_recall_curve(arr_true, arr_score)
        return float(auc(rec, prec))
    except Exception:
        return 0.0


def calculate_brier_score(y_true: Sequence[float], y_score: Sequence[float]) -> float:
    """Calculate Brier score after min-max normalizing scores to [0, 1]."""
    arr_true = np.asarray(y_true)
    arr_score = np.asarray(y_score)
    min_v, max_v = np.min(arr_score), np.max(arr_score)
    if max_v > min_v:
        prob_est = (arr_score - min_v) / (max_v - min_v)
    else:
        prob_est = np.full_like(arr_score, 0.5)
    return float(brier_score_loss(arr_true, prob_est))


def evaluate_prediction(
    prediction: PredictionResult,
    clinical_df: pl.DataFrame,
    cohort_id: str,
    cancer_type: str,
    time_stratum: str = "Pre",
    response_stratum: str = "Standard",
    pooling_strategy: str = "cohort",
    n_bootstrap: int = 1000,
) -> Result[CohortBenchmarkResult, str]:
    """Join predictions with clinical binary response and calculate discrimination metrics."""
    try:
        joined = prediction.predictions.join(
            clinical_df.select([pl.col("sample_id"), pl.col("response_binary")]),
            on="sample_id",
            how="inner",
        ).filter(pl.col("response_binary").is_finite() & pl.col("score").is_finite())

        if joined.is_empty():
            return Failure(f"No valid annotated samples to evaluate {prediction.predictor_name} on {cohort_id}.")

        y_true = joined["response_binary"].to_list()
        y_score = joined["score"].to_list()

        n_resp = int(sum(y_true))
        n_total = len(y_true)

        disc_res = calculate_discrimination_bootstrap(
            y_true, y_score, n_bootstrap=n_bootstrap
        )
        brier = calculate_brier_score(y_true, y_score)

        return Success(
            CohortBenchmarkResult(
                cohort_id=cohort_id,
                cancer_type=cancer_type,
                predictor_name=prediction.predictor_name,
                category=prediction.category,
                n_samples=n_total,
                n_responders=n_resp,
                n_non_responders=n_total - n_resp,
                time_stratum=time_stratum,
                response_stratum=response_stratum,
                roc_auc=disc_res["roc_auc"],
                roc_auc_ci_lower=disc_res["roc_auc_ci_lower"],
                roc_auc_ci_upper=disc_res["roc_auc_ci_upper"],
                p_value_vs_half=disc_res["p_value_vs_half"],
                pr_auc=disc_res["pr_auc"],
                pr_auc_ci_lower=disc_res["pr_auc_ci_lower"],
                pr_auc_ci_upper=disc_res["pr_auc_ci_upper"],
                baseline_prevalence=disc_res["baseline_prevalence"],
                delta_pr_auc=disc_res["delta_pr_auc"],
                p_value_prauc=disc_res["p_value_prauc"],
                pooling_strategy=pooling_strategy,
                brier_score=brier,
            )
        )
    except Exception as exc:
        return Failure(f"Evaluation failed for {prediction.predictor_name} on {cohort_id}: {exc}")
