"""Unit tests for compositional transforms, feature selection, and downstream classifier HPO."""

from __future__ import annotations

import numpy as np
import polars as pl
import pytest
from returns.result import Success

from tme_datasets.sampling.classifier_hpo import (
    ClassifierType,
    FeatureSelectorType,
    InnerClassifierConfig,
    InnerOptimizationResult,
    TransformType,
    apply_compositional_transform,
    evaluate_inner_pipeline,
    run_inner_classifier_hpo,
)


def test_compositional_transforms() -> None:
    """Verify mathematical properties of simplex transforms."""
    # 4 samples x 5 cell states on probability simplex
    X = np.array([
        [0.4, 0.3, 0.1, 0.1, 0.1],
        [0.0, 0.5, 0.2, 0.2, 0.1],  # Contains exact 0.0
        [0.2, 0.2, 0.2, 0.2, 0.2],
        [0.7, 0.1, 0.1, 0.05, 0.05],
    ])

    # 1. Raw
    X_raw = apply_compositional_transform(X, TransformType.RAW)
    assert np.allclose(X_raw, X)

    # 2. CLR: Sum across states for each sample must be zero
    X_clr = apply_compositional_transform(X, TransformType.CLR, eps=1e-5)
    assert X_clr.shape == X.shape
    assert not np.any(np.isnan(X_clr))
    assert not np.any(np.isinf(X_clr))
    clr_sums = np.sum(X_clr, axis=1)
    assert np.allclose(clr_sums, 0.0, atol=1e-5)

    # 3. Logit
    X_logit = apply_compositional_transform(X, TransformType.LOGIT, eps=1e-5)
    assert X_logit.shape == X.shape
    assert not np.any(np.isnan(X_logit))
    assert not np.any(np.isinf(X_logit))


def test_inner_classifier_pipeline_and_hpo() -> None:
    """Verify inner Leave-One-Cohort-Out optimization across synthetic cell fractions."""
    rng = np.random.default_rng(999)
    n_states = 6
    state_names = tuple(f"State_{i}" for i in range(n_states))
    n_samples_per_cohort = 25

    cohort_fractions: dict[str, np.ndarray] = {}
    cohort_responses: dict[str, np.ndarray] = {}

    for cid in ("DiscoveryA", "DiscoveryB"):
        y = np.array([1] * 12 + [0] * 13)
        fracs = np.zeros((n_samples_per_cohort, n_states))

        for i in range(n_samples_per_cohort):
            if y[i] == 1:
                # Responders have elevated State 0 (e.g., CD8_Tex)
                w = np.array([0.5, 0.1, 0.1, 0.1, 0.1, 0.1])
            else:
                # Non-responders have elevated State 5 (e.g., M2 Macrophage)
                w = np.array([0.05, 0.15, 0.15, 0.15, 0.15, 0.35])

            noise = rng.uniform(0.01, 0.05, size=n_states)
            vec = w + noise
            fracs[i, :] = vec / vec.sum()

        cohort_fractions[cid] = fracs
        cohort_responses[cid] = y

    # Run single pipeline configuration
    cfg = InnerClassifierConfig(
        transform=TransformType.CLR,
        selector=FeatureSelectorType.NONE,
        classifier=ClassifierType.LOGISTIC_REGRESSION,
    )
    single_res = evaluate_inner_pipeline(
        cohort_fractions=cohort_fractions,
        cohort_responses=cohort_responses,
        feature_names=state_names,
        config=cfg,
    )
    assert isinstance(single_res, Success)
    mean_auc, mean_pr, c_aucs, sel_feats = single_res.unwrap()
    assert mean_auc >= 0.75
    assert len(c_aucs) == 2

    # Run full high-speed inner HPO across all candidate configurations
    hpo_res = run_inner_classifier_hpo(
        cohort_fractions=cohort_fractions,
        cohort_responses=cohort_responses,
        feature_names=state_names,
    )
    assert isinstance(hpo_res, Success)
    inner_summary = hpo_res.unwrap()

    assert inner_summary.best_loco_auc >= 0.75
    assert isinstance(inner_summary.all_results_df, pl.DataFrame)
    assert inner_summary.all_results_df.height >= 5
    assert inner_summary.elapsed_seconds >= 0.0
