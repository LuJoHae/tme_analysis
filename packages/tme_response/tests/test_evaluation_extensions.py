"""Tests for cross-cohort pooling, meta-analysis, and cohort predictability calculations."""

from __future__ import annotations

import numpy as np
import polars as pl
from returns.result import Success

from tme_response.evaluation.pooling import (
    compute_meta_analytic_auc,
    pool_cohort_stratum_data,
    standardize_prediction_scores,
)
from tme_response.evaluation.predictability import (
    UNIVERSAL_RNA_PREDICTORS,
    compute_cohort_predictability,
)
from tme_response.schemas import PredictionResult
from tme_response.types import PredictorCategory


def test_standardize_prediction_scores() -> None:
    # 1. Normal variable
    scores = [10.0, 20.0, 30.0, 40.0, 50.0]
    df = pl.DataFrame({"sample_id": [f"S{i}" for i in range(5)], "score": scores})
    pred = PredictionResult(
        predictor_name="GEP",
        category=PredictorCategory.SIGNATURE,
        predictions=df,
    )
    std_pred = standardize_prediction_scores(pred)
    res_scores = std_pred.predictions["score"].to_numpy()

    assert np.isclose(np.mean(res_scores), 0.0, atol=1e-6)
    assert np.isclose(np.std(res_scores), 1.0, atol=1e-6)

    # 2. Constant variable
    df_const = pl.DataFrame({"sample_id": ["S1", "S2"], "score": [5.0, 5.0]})
    pred_const = PredictionResult(
        predictor_name="Const",
        category=PredictorCategory.SIGNATURE,
        predictions=df_const,
    )
    std_const = standardize_prediction_scores(pred_const)
    assert (std_const.predictions["score"].to_numpy() == 0.0).all()


def test_pool_cohort_stratum_data() -> None:
    # Cohort 1
    c1_clin = pl.DataFrame({"sample_id": ["A1", "A2"], "response_binary": [1.0, 0.0]})
    c1_pred = PredictionResult(
        predictor_name="CYT",
        category=PredictorCategory.SIGNATURE,
        predictions=pl.DataFrame({"sample_id": ["A1", "A2"], "score": [2.0, 1.0]}),
    )

    # Cohort 2 (study with higher baseline expression)
    c2_clin = pl.DataFrame({"sample_id": ["B1", "B2"], "response_binary": [1.0, 0.0]})
    c2_pred = PredictionResult(
        predictor_name="CYT",
        category=PredictorCategory.SIGNATURE,
        predictions=pl.DataFrame({"sample_id": ["B1", "B2"], "score": [12.0, 11.0]}),
    )

    # Pool with standardization
    match pool_cohort_stratum_data(
        [("C1", c1_clin, [c1_pred]), ("C2", c2_clin, [c2_pred])],
        group_id="Group1",
        cancer_type="Melanoma",
        standardize=True,
    ):
        case Success((pooled_clin, pooled_preds)):
            assert len(pooled_clin) == 4
            assert set(pooled_clin["sample_id"].to_list()) == {
                "C1::A1",
                "C1::A2",
                "C2::B1",
                "C2::B2",
            }
            assert len(pooled_preds) == 1
            cyt_preds = pooled_preds[0]
            assert cyt_preds.predictor_name == "CYT"
            # In both cohorts, responder score was higher than non-responder
            # Standardized: S(A1) > S(A2) and S(B1) > S(B2)
            c_scores = dict(
                zip(
                    cyt_preds.predictions["sample_id"].to_list(),
                    cyt_preds.predictions["score"].to_list(),
                )
            )
            assert c_scores["C1::A1"] > c_scores["C1::A2"]
            assert c_scores["C2::B1"] > c_scores["C2::B2"]
            # After z-scoring, responders in C1 and C2 have identical relative values
            assert np.isclose(c_scores["C1::A1"], c_scores["C2::B1"], atol=1e-5)
        case _ as fail:
            assert False, f"Pooling failed: {fail}"


def test_compute_meta_analytic_auc() -> None:
    # 3 cohorts with consistent AUC ~ 0.70
    aucs = [0.70, 0.72, 0.68]
    ci_lowers = [0.60, 0.62, 0.58]
    ci_uppers = [0.80, 0.82, 0.78]

    res = compute_meta_analytic_auc(aucs, ci_lowers, ci_uppers)
    assert 0.69 <= res["meta_auc_fixed"] <= 0.71
    assert 0.69 <= res["meta_auc_random"] <= 0.71
    assert res["i_squared"] < 10.0  # minimal heterogeneity


def test_compute_cohort_predictability() -> None:
    # Create mock benchmark df
    records = []
    for p in UNIVERSAL_RNA_PREDICTORS:
        records.append({
            "cohort_id": "CohortA",
            "cancer_type": "Melanoma",
            "time_stratum": "Pre",
            "response_stratum": "Standard",
            "pooling_strategy": "cohort",
            "predictor_name": p,
            "roc_auc": 0.80 if p == "GEP" else 0.60,
            "pr_auc": 0.55 if p == "GEP" else 0.40,
            "delta_pr_auc": 0.25 if p == "GEP" else 0.10,
            "baseline_prevalence": 0.30,
            "n_samples": 50,
            "n_responders": 15,
        })
    # Add non-universal predictor
    records.append({
        "cohort_id": "CohortA",
        "cancer_type": "Melanoma",
        "time_stratum": "Pre",
        "response_stratum": "Standard",
        "pooling_strategy": "cohort",
        "predictor_name": "TMB_Genomic",
        "roc_auc": 0.85,
        "pr_auc": 0.60,
        "delta_pr_auc": 0.30,
        "baseline_prevalence": 0.30,
        "n_samples": 50,
        "n_responders": 15,
    })

    df = pl.DataFrame(records)
    match compute_cohort_predictability(df):
        case Success(pred_df):
            assert len(pred_df) == 1
            row = pred_df.to_dicts()[0]
            assert row["cohort_id"] == "CohortA"
            assert row["best_predictor_rna"] == "GEP"
            assert row["max_roc_auc_rna"] == 0.80
            # 7 predictors at 0.60, 1 at 0.80 -> mean = (7*0.60 + 0.80)/8 = 5.0/8 = 0.625
            assert np.isclose(row["mean_roc_auc_rna"], 0.625)
            # best overall is TMB at 0.85
            assert row["best_predictor_all"] == "TMB_Genomic"
            assert row["max_roc_auc_all"] == 0.85
            assert row["n_predictors_evaluated"] == 9
        case _ as fail:
            assert False, f"Predictability computation failed: {fail}"
