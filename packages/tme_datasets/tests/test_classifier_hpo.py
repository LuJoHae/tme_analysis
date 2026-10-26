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


def test_pairwise_log_ratios_transform() -> None:
    """Verify pairwise log-ratios expansion and feature name alignment."""
    from tme_datasets.sampling.classifier_hpo import get_transformed_feature_names

    X = np.array([
        [0.4, 0.3, 0.3],
        [0.1, 0.2, 0.7],
    ])
    feature_names = ("StateA", "StateB", "StateC")

    X_ratios = apply_compositional_transform(X, TransformType.CLR_RATIOS, eps=1e-5)
    # Original 3 features + 3 pairwise differences = 6 features
    assert X_ratios.shape == (2, 6)

    transformed_names = get_transformed_feature_names(feature_names, TransformType.CLR_RATIOS)
    assert len(transformed_names) == 6
    assert transformed_names[:3] == feature_names
    assert "StateA_div_StateB" in transformed_names
    assert "StateA_div_StateC" in transformed_names
    assert "StateB_div_StateC" in transformed_names


def test_warm_start_richardson_lucy() -> None:
    """Verify warm-start Richardson-Lucy deconvolution is numerically identical to continuous iteration."""
    from tme_datasets.sampling.hpo_optimizer import batch_richardson_lucy_deconvolution

    rng = np.random.default_rng(42)
    n_samples, n_genes, n_states = 10, 50, 4

    bulk = rng.uniform(10, 500, size=(n_samples, n_genes))
    ref = rng.uniform(1, 100, size=(n_states, n_genes))
    ref = ref / ref.sum(axis=1, keepdims=True)

    # 1. 40 iterations continuous from scratch
    theta_continuous = batch_richardson_lucy_deconvolution(bulk, ref, n_iter=40)

    # 2. 20 iterations + 20 iterations warm start
    theta_20 = batch_richardson_lucy_deconvolution(bulk, ref, n_iter=20)
    theta_warm = batch_richardson_lucy_deconvolution(bulk, ref, n_iter=20, initial_theta=theta_20)

    # Must be numerically identical within floating point tolerance
    assert np.allclose(theta_continuous, theta_warm, atol=1e-12)


def test_single_cohort_stratified_cv_with_tree_ensembles() -> None:
    """Verify single-cohort pooled OOF CV works seamlessly with HistGradientBoosting and ExtraTrees."""
    from tme_datasets.sampling.classifier_hpo import evaluate_single_cohort_stratified_cv

    rng = np.random.default_rng(123)
    n_samples = 40
    n_states = 5
    feature_names = tuple(f"State_{i}" for i in range(n_states))

    y = np.array([1] * 20 + [0] * 20)
    frac_mat = np.zeros((n_samples, n_states))
    for i in range(n_samples):
        w = np.array([0.4, 0.2, 0.2, 0.1, 0.1]) if y[i] == 1 else np.array([0.1, 0.1, 0.2, 0.2, 0.4])
        noise = rng.uniform(0.01, 0.05, size=n_states)
        vec = w + noise
        frac_mat[i, :] = vec / vec.sum()

    for clf_type in (ClassifierType.HIST_GRADIENT_BOOSTING, ClassifierType.EXTRA_TREES):
        cfg = InnerClassifierConfig(
            transform=TransformType.CLR,
            selector=FeatureSelectorType.NONE,
            classifier=clf_type,
            n_jobs=2,
        )
        res = evaluate_single_cohort_stratified_cv(
            frac_mat=frac_mat,
            y=y,
            feature_names=feature_names,
            config=cfg,
            n_splits=5,
            n_repeats=1,
            seed=42,
        )
        assert isinstance(res, Success)
        auc_val, pr_val, sel_feats = res.unwrap()
        assert 0.0 <= auc_val <= 1.0
        assert 0.0 <= pr_val <= 1.0

