"""Unit tests for single-cohort reference sampling and classifier HPO."""

import numpy as np
import pytest
from returns.result import Success

from tme_datasets.sampling.classifier_hpo import (
    ClassifierType,
    FeatureSelectorType,
    InnerClassifierConfig,
    TransformType,
    evaluate_single_cohort_stratified_cv,
    run_single_cohort_classifier_hpo,
)


def test_evaluate_single_cohort_stratified_cv() -> None:
    """Test that stratified cross-validation computes out-of-fold metrics across parallel folds."""
    rng = np.random.default_rng(42)
    n_samples = 40
    n_states = 6

    # Generate synthetic Dirichlet-distributed fractions (sum to 1)
    raw = rng.exponential(scale=1.0, size=(n_samples, n_states))
    frac_mat = raw / raw.sum(axis=1, keepdims=True)

    # Correlate response with state 0
    logits = 4.0 * frac_mat[:, 0] - 0.5 + rng.normal(scale=0.2, size=n_samples)
    y = (logits > np.median(logits)).astype(np.int64)

    state_labels = [f"State_{i}" for i in range(n_states)]

    cfg = InnerClassifierConfig(
        transform=TransformType.LOGIT,
        selector=FeatureSelectorType.NONE,
        classifier=ClassifierType.LOGISTIC_REGRESSION,
        c_regularization=1.0,
        n_jobs=2,
    )

    res = evaluate_single_cohort_stratified_cv(
        frac_mat=frac_mat,
        y=y,
        feature_names=state_labels,
        config=cfg,
        n_splits=5,
        n_repeats=2,
        seed=42,
    )

    assert isinstance(res, Success)
    mean_auc, mean_pr, sel_features = res.unwrap()
    assert 0.0 <= mean_auc <= 1.0
    assert 0.0 <= mean_pr <= 1.0
    assert len(sel_features) > 0


def test_run_single_cohort_classifier_hpo() -> None:
    """Test inner loop search over transforms and classifiers on single cohort."""
    rng = np.random.default_rng(123)
    n_samples = 30
    n_states = 5

    raw = rng.exponential(scale=1.0, size=(n_samples, n_states))
    frac_mat = raw / raw.sum(axis=1, keepdims=True)
    y = rng.choice([0, 1], size=n_samples, p=[0.5, 0.5])

    state_labels = [f"State_{i}" for i in range(n_states)]

    candidate_cfgs = (
        InnerClassifierConfig(
            transform=TransformType.RAW,
            selector=FeatureSelectorType.NONE,
            classifier=ClassifierType.LOGISTIC_REGRESSION,
            n_jobs=2,
        ),
        InnerClassifierConfig(
            transform=TransformType.LOGIT,
            selector=FeatureSelectorType.NONE,
            classifier=ClassifierType.LOGISTIC_REGRESSION,
            n_jobs=2,
        ),
    )

    res = run_single_cohort_classifier_hpo(
        frac_mat=frac_mat,
        y=y,
        feature_names=state_labels,
        candidate_configs=candidate_cfgs,
        n_splits=3,
        n_repeats=1,
        seed=42,
        n_jobs=2,
    )

    assert isinstance(res, Success)
    opt = res.unwrap()
    assert opt.best_loco_auc >= 0.0
    assert len(opt.selected_features) > 0
    assert opt.all_results_df.shape[0] == 2
