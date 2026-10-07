"""Decision Curve Analysis (DCA) and clinical utility benchmarking."""

from __future__ import annotations

from typing import Literal, Sequence
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from sklearn.linear_model import LogisticRegression


class DCAResultRecord(BaseModel):
    """Point evaluation for a strategy at a specific threshold probability."""
    model_config = ConfigDict(frozen=True)

    threshold: float
    strategy: str  # 'model', 'treat_all', 'treat_none'
    net_benefit: float
    true_positives: int
    false_positives: int
    interventions_avoided_per_100: float


def calibrate_probabilities(
    y_true: np.ndarray,
    y_score: np.ndarray,
    method: Literal["logistic", "minmax"] = "logistic",
) -> np.ndarray:
    """Calibrate continuous biomarker scores into well-behaved predicted probabilities in [0, 1]."""
    if method == "logistic":
        # Fit univariable logistic regression (Platt scaling)
        clf = LogisticRegression(solver="lbfgs")
        clf.fit(y_score.reshape(-1, 1), y_true)
        probs = clf.predict_proba(y_score.reshape(-1, 1))[:, 1]
        return np.asarray(probs, dtype=float)
    else:
        # Min-max normalization
        s_min, s_max = float(np.min(y_score)), float(np.max(y_score))
        if s_max > s_min:
            return (y_score - s_min) / (s_max - s_min)
        return np.full_like(y_score, 0.5)


def calculate_decision_curve(
    y_true: Sequence[float | int | bool],
    y_score: Sequence[float],
    predictor_name: str = "Biomarker",
    thresholds: Sequence[float] | None = None,
    calibration: Literal["logistic", "minmax"] = "logistic",
) -> Result[pl.DataFrame, str]:
    """Compute Decision Curve Analysis (DCA) comparing the biomarker model against Treat All and Treat None.

    Parameters
    ----------
    y_true : Sequence[float | int | bool]
        Binary ground truth clinical response (1 = responder, 0 = non-responder).
    y_score : Sequence[float]
        Continuous predictor score.
    predictor_name : str
        Display name of the biomarker predictor.
    thresholds : Sequence[float] | None
        Grid of decision threshold probabilities (defaults to 0.05 to 0.80 by step 0.025).
    calibration : 'logistic' | 'minmax'
        Method to map continuous scores to calibrated probabilities.

    Returns
    -------
    Result[pl.DataFrame, str]
        Polars DataFrame containing threshold, strategy, net_benefit, and interventions_avoided.
    """
    try:
        t_arr = np.asarray(y_true, dtype=int)
        s_arr = np.asarray(y_score, dtype=float)

        valid_mask = np.isfinite(t_arr) & np.isfinite(s_arr)
        t_clean = t_arr[valid_mask]
        s_clean = s_arr[valid_mask]

        n_samples = len(t_clean)
        n_pos = int(np.sum(t_clean))

        if n_samples < 5 or n_pos < 1 or n_pos == n_samples:
            return Failure(f"Invalid outcome distribution for DCA: N={n_samples}, Positives={n_pos}")

        if thresholds is None:
            threshold_grid = np.linspace(0.05, 0.80, 31)
        else:
            threshold_grid = np.asarray(thresholds, dtype=float)

        prob_pred = calibrate_probabilities(t_clean, s_clean, method=calibration)
        prevalence = n_pos / n_samples

        records: list[dict[str, object]] = []

        for p_t in threshold_grid:
            if p_t <= 0.0 or p_t >= 1.0:
                continue

            weight = p_t / (1.0 - p_t)

            # 1. Model Strategy: treat if prob_pred >= p_t
            pred_positive = prob_pred >= p_t
            tp_model = int(np.sum((pred_positive) & (t_clean == 1)))
            fp_model = int(np.sum((pred_positive) & (t_clean == 0)))
            nb_model = (tp_model / n_samples) - (fp_model / n_samples) * weight

            # 2. Treat All Strategy
            tp_all = n_pos
            fp_all = n_samples - n_pos
            nb_all = (tp_all / n_samples) - (fp_all / n_samples) * weight

            # 3. Treat None Strategy
            nb_none = 0.0

            # Interventions avoided per 100 patients compared to Treat All
            avoided = ((nb_model - nb_all) / weight) * 100.0 if weight > 0 else 0.0

            records.append({
                "threshold": float(p_t),
                "strategy": predictor_name,
                "net_benefit": float(nb_model),
                "true_positives": tp_model,
                "false_positives": fp_model,
                "interventions_avoided_per_100": float(avoided),
            })
            records.append({
                "threshold": float(p_t),
                "strategy": "Treat All",
                "net_benefit": float(nb_all),
                "true_positives": tp_all,
                "false_positives": fp_all,
                "interventions_avoided_per_100": 0.0,
            })
            records.append({
                "threshold": float(p_t),
                "strategy": "Treat None",
                "net_benefit": float(nb_none),
                "true_positives": 0,
                "false_positives": 0,
                "interventions_avoided_per_100": float(fp_all / n_samples * 100.0),
            })

        return Success(pl.DataFrame(records))
    except Exception as exc:
        return Failure(f"DCA calculation failed for {predictor_name}: {exc}")
