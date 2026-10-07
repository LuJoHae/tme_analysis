"""Tests for survival analysis (C-index, Cox Proportional Hazards) and Decision Curve Analysis (DCA)."""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest
from returns.result import Failure, Success

from tme_response.evaluation.dca import (
    calculate_decision_curve,
    calibrate_probabilities,
)
from tme_response.evaluation.survival import (
    CIndexResult,
    CoxHazardResult,
    calculate_c_index_bootstrap,
    compute_c_index_raw,
    fit_univariable_cox,
)
from tme_response.visualization.survival import (
    create_c_index_forest_plot,
    create_cox_hr_forest_plot,
    create_dca_net_benefit_chart,
)


def test_compute_c_index_raw_perfect_and_inverse() -> None:
    # Patient 1: survived 5 months (event 1), score 1.0 (low score = short survival)
    # Patient 2: survived 15 months (event 1), score 3.0
    # Patient 3: survived 30 months (event 0, censored), score 5.0
    times = np.array([5.0, 15.0, 30.0])
    events = np.array([1, 1, 0])
    scores = np.array([1.0, 3.0, 5.0])

    c_idx, n_comp = compute_c_index_raw(times, events, scores, higher_is_better=True)
    assert c_idx == 1.0
    assert n_comp == 3

    # Inverse ordering
    scores_inv = np.array([5.0, 3.0, 1.0])
    c_idx_inv, n_comp_inv = compute_c_index_raw(times, events, scores_inv, higher_is_better=True)
    assert c_idx_inv == 0.0
    assert n_comp_inv == 3


def test_calculate_c_index_bootstrap_success() -> None:
    # Synthetic cohort where higher score strongly correlates with longer survival
    rng = np.random.default_rng(42)
    n = 50
    scores = rng.normal(0, 1, size=n)
    times = np.exp(2.0 + 0.8 * scores + rng.normal(0, 0.3, size=n))
    events = np.ones(n, dtype=int)
    # Censor 20%
    events[rng.random(n) < 0.2] = 0

    res = calculate_c_index_bootstrap(times, events, scores, n_bootstrap=100, seed=42)
    assert isinstance(res, Success)
    result = res.unwrap()
    assert isinstance(result, CIndexResult)
    assert result.c_index > 0.65
    assert 0.0 <= result.ci_lower <= result.c_index <= result.ci_upper <= 1.0
    assert result.p_value_vs_half < 0.05
    assert result.n_samples == n


def test_calculate_c_index_bootstrap_insufficient_data() -> None:
    res = calculate_c_index_bootstrap([5.0, 10.0], [1, 0], [1.0, 2.0])
    assert isinstance(res, Failure)


def test_fit_univariable_cox_protective() -> None:
    # Biomarker where higher score reduces hazard (longer survival, HR < 1.0)
    rng = np.random.default_rng(123)
    n = 60
    scores = rng.normal(0, 1, size=n)
    times = np.exp(3.0 + 0.7 * scores + rng.normal(0, 0.4, size=n))
    events = np.ones(n, dtype=int)
    events[:10] = 0  # 10 censored

    res = fit_univariable_cox(times, events, scores, standardize_score=True)
    assert isinstance(res, Success)
    cox_res = res.unwrap()
    assert isinstance(cox_res, CoxHazardResult)
    assert cox_res.hazard_ratio < 1.0
    assert cox_res.coefficient < 0.0
    assert cox_res.hr_ci_lower < cox_res.hazard_ratio < cox_res.hr_ci_upper
    assert cox_res.p_value < 0.05
    assert cox_res.n_samples == n


def test_fit_univariable_cox_deleterious() -> None:
    # Biomarker where higher score increases hazard (shorter survival, HR > 1.0)
    rng = np.random.default_rng(123)
    n = 60
    scores = rng.normal(0, 1, size=n)
    times = np.exp(3.0 - 0.7 * scores + rng.normal(0, 0.4, size=n))
    events = np.ones(n, dtype=int)

    res = fit_univariable_cox(times, events, scores, standardize_score=True)
    assert isinstance(res, Success)
    cox_res = res.unwrap()
    assert cox_res.hazard_ratio > 1.0
    assert cox_res.coefficient > 0.0


def test_fit_univariable_cox_zero_variance() -> None:
    times = [10.0] * 20
    events = [1] * 20
    scores = [2.0] * 20
    res = fit_univariable_cox(times, events, scores)
    assert isinstance(res, Failure)


def test_calibrate_probabilities() -> None:
    y_true = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    y_score = np.array([0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9])

    probs_log = calibrate_probabilities(y_true, y_score, method="logistic")
    assert probs_log.shape == (8,)
    assert (probs_log >= 0.0).all() and (probs_log <= 1.0).all()
    # Monotonic order should be preserved
    assert probs_log[0] < probs_log[-1]

    probs_mm = calibrate_probabilities(y_true, y_score, method="minmax")
    assert np.isclose(probs_mm[0], 0.0)
    assert np.isclose(probs_mm[-1], 1.0)


def test_calculate_decision_curve_structure() -> None:
    y_true = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
    y_score = [-1.5, -1.0, -0.5, 0.0, 0.2, 0.4, 0.8, 1.2, 1.5, 2.0]

    res = calculate_decision_curve(y_true, y_score, predictor_name="CYT", thresholds=[0.1, 0.2, 0.5])
    assert isinstance(res, Success)
    df = res.unwrap()
    assert isinstance(df, pl.DataFrame)

    required_cols = {"threshold", "strategy", "net_benefit", "interventions_avoided_per_100"}
    assert required_cols.issubset(set(df.columns))

    strategies = set(df["strategy"].unique().to_list())
    assert {"CYT", "Treat All", "Treat None"} == strategies

    # Net benefit of Treat None is always 0.0
    none_df = df.filter(pl.col("strategy") == "Treat None")
    assert (none_df["net_benefit"] == 0.0).all()


def test_survival_visualizations() -> None:
    # 1. C-index plot
    c_df = pl.DataFrame({
        "predictor_name": ["CYT", "GEP", "TIDE"],
        "category": ["Signature", "Signature", "Model"],
        "endpoint": ["OS", "OS", "OS"],
        "c_index": [0.65, 0.62, 0.58],
        "ci_lower": [0.55, 0.51, 0.47],
        "ci_upper": [0.74, 0.72, 0.69],
        "p_value": [0.01, 0.04, 0.12],
    })
    c_chart_res = create_c_index_forest_plot(c_df, endpoint="OS")
    assert isinstance(c_chart_res, Success)

    # 2. Cox HR plot
    cox_df = pl.DataFrame({
        "predictor_name": ["CYT", "GEP", "TIDE"],
        "category": ["Signature", "Signature", "Model"],
        "endpoint": ["OS", "OS", "OS"],
        "hazard_ratio": [0.68, 0.74, 0.95],
        "hr_ci_lower": [0.52, 0.58, 0.72],
        "hr_ci_upper": [0.89, 0.94, 1.25],
        "p_value": [0.005, 0.012, 0.70],
    })
    cox_chart_res = create_cox_hr_forest_plot(cox_df, endpoint="OS")
    assert isinstance(cox_chart_res, Success)

    # 3. DCA chart
    dca_df = pl.DataFrame({
        "threshold": [0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3],
        "strategy": ["CYT", "CYT", "CYT", "Treat All", "Treat All", "Treat All", "Treat None", "Treat None", "Treat None"],
        "net_benefit": [0.35, 0.28, 0.20, 0.30, 0.20, 0.05, 0.0, 0.0, 0.0],
        "interventions_avoided_per_100": [10.0, 15.0, 20.0, 0.0, 0.0, 0.0, 50.0, 50.0, 50.0],
    })
    dca_chart_res = create_dca_net_benefit_chart(dca_df)
    assert isinstance(dca_chart_res, Success)
