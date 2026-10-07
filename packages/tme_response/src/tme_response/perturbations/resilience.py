"""Pure functional calculation of the Perturbation Resilience Index (PRI) and noise decay metrics."""

from __future__ import annotations

import numpy as np
import polars as pl
from returns.result import Failure, Result, Success

from .models import PerturbationResilienceRecord


def compute_perturbation_resilience(
    perturbation_df: pl.DataFrame,
) -> Result[pl.DataFrame, str]:
    """Calculate the Perturbation Resilience Index (PRI) for each predictor across perturbation decay curves.

    PRI is defined as the normalized area under the AUC-retention curve:
    PRI = trapezoid_integral(AUC(x) dx) / ((x_max - x_min) * AUC_baseline)

    PRI = 1.0 represents perfect noise invariance (no degradation).
    PRI < 1.0 measures the rate and severity of performance collapse.
    """
    if perturbation_df.is_empty():
        return Failure("Cannot compute resilience from empty perturbation DataFrame.")

    required_cols = {
        "cohort_id",
        "predictor_name",
        "category",
        "perturbation_type",
        "intensity",
        "roc_auc",
    }
    missing = required_cols - set(perturbation_df.columns)
    if missing:
        return Failure(f"Missing required columns for resilience computation: {sorted(missing)}")

    group_keys = ["cohort_id", "predictor_name", "category", "perturbation_type"]
    records: list[dict[str, object]] = []

    for group_vals, sub_df in perturbation_df.group_by(group_keys):
        c_id, p_name, cat, p_type = group_vals

        sorted_sub = sub_df.sort("intensity")
        xs = sorted_sub["intensity"].to_numpy().astype(float)
        ys = sorted_sub["roc_auc"].to_numpy().astype(float)

        if len(xs) < 2:
            continue

        baseline_idx = int(np.argmin(xs)) if p_type != "dilution" else int(np.argmax(xs))
        baseline_auc = float(ys[baseline_idx])
        min_auc = float(np.min(ys))

        # Integrate trapezoidal area under the curve
        dx = float(np.max(xs) - np.min(xs))
        if dx > 1e-6 and baseline_auc > 0.05:
            area = float(np.trapezoid(ys, xs))
            ideal_area = baseline_auc * dx
            pri = float(np.clip(area / ideal_area, 0.0, 1.2))
        else:
            pri = 1.0

        worst_idx = int(np.argmax(xs)) if p_type != "dilution" else int(np.argmin(xs))
        retention = float(ys[worst_idx] / baseline_auc) if baseline_auc > 1e-4 else 1.0

        records.append({
            "cohort_id": c_id,
            "predictor_name": p_name,
            "category": cat,
            "perturbation_type": p_type,
            "baseline_roc_auc": round(baseline_auc, 4),
            "min_roc_auc": round(min_auc, 4),
            "pri_score": round(pri, 4),
            "relative_auc_retained": round(retention, 4),
        })

    if not records:
        return Failure("No valid multi-intensity perturbation curves found.")

    resilience_df = pl.DataFrame(records).sort(
        by=["cohort_id", "perturbation_type", "pri_score"],
        descending=[False, False, True],
    )

    return Success(resilience_df)
