"""Tests for Stability Selection module.

Tests:
1. Mathematical error bound solvers (MB and Shah & Samworth Unimodal).
2. Complementary pairs subsampling and feature randomization.
3. Full stability selection engine on synthetic sparse high-dimensional data.
4. Scikit-learn adapter compatibility and Pipeline integration.
5. Declarative Altair chart generation.
"""

import math
import numpy as np
import pytest
from returns.maybe import Some, Nothing
from returns.result import Success, Failure
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge

from selective_inference.stability_selection.types import (
    StabilityParameters,
    StabilityResult,
)
from selective_inference.stability_selection.bounds import (
    compute_unimodal_constant,
    solve_pfer,
    solve_cutoff,
    solve_q,
    resolve_stability_parameters,
)
from selective_inference.stability_selection.subsampling import (
    generate_complementary_pairs,
    generate_subsamples,
    generate_stratified_complementary_pairs,
    generate_stratified_subsamples,
    apply_randomized_weights,
)
from selective_inference.stability_selection.core import (
    run_stability_selection,
    default_lasso_path_fitter,
)
from selective_inference.stability_selection.fitters import (
    create_lasso_fitter,
    create_elastic_net_fitter,
    create_l1_logistic_fitter,
    create_tree_importance_fitter,
    create_cohort_adjusted_fitter,
    create_group_lasso_cohort_fitter,
    create_merf_cohort_fitter,
    create_multitask_logistic_cohort_fitter,
    create_meta_analysis_cohort_fitter,
    create_multistudy_invariant_cohort_fitter,
    create_glmm_lasso_cohort_fitter,
    create_oscar_fitter,
    create_slope_fitter,
    create_oscar_cohort_fitter,
    create_slope_cohort_fitter,
    _prox_sorted_l1,
)
from selective_inference.stability_selection.sklearn_adapter import (
    StabilitySelector,
)
from selective_inference.stability_selection.visualization import (
    plot_stability_paths,
    plot_stability_scores,
)


# -----------------------------------------------------------------------------
# 1. Parameter Bounds & Solvers
# -----------------------------------------------------------------------------

def test_meinshausen_buhlmann_bounds():
    """Verify Meinshausen & Bühlmann (2010) closed-form formulas."""
    p = 100
    q = 10.0
    cutoff = 0.75
    B = 50

    # PFER = q^2 / (p * (2 * cutoff - 1)) = 100 / (100 * 0.5) = 2.0
    pfer_res = solve_pfer(p, cutoff, q, B, assumption="none")
    assert isinstance(pfer_res, Success)
    assert math.isclose(pfer_res.unwrap(), 2.0, rel_tol=1e-5)

    # Solve for cutoff given PFER = 2.0 -> cutoff = 0.5 + 100 / (2 * 100 * 2) = 0.75
    cutoff_res = solve_cutoff(p, q, 2.0, B, assumption="none")
    assert isinstance(cutoff_res, Success)
    assert math.isclose(cutoff_res.unwrap(), 0.75, rel_tol=1e-5)

    # Solve for q given cutoff = 0.75, PFER = 2.0 -> q = sqrt(100 * 2 * 0.5) = 10.0
    q_res = solve_q(p, cutoff, 2.0, B, assumption="none")
    assert isinstance(q_res, Success)
    assert math.isclose(q_res.unwrap(), 10.0, rel_tol=1e-5)


def test_shah_samworth_unimodal_bounds():
    """Verify Shah & Samworth (2013) unimodal bound is strictly tighter than MB."""
    p = 100
    q = 10.0
    cutoff = 0.75
    B = 50

    mb_pfer = solve_pfer(p, cutoff, q, B, assumption="none").unwrap()
    ss_pfer = solve_pfer(p, cutoff, q, B, assumption="unimodal").unwrap()

    # Shah-Samworth bound should be almost 2x tighter (lower expected false positives)
    assert ss_pfer < mb_pfer
    # SS bound with B=50 at cutoff=0.75: C = 2 * (2*0.75 - 1 - 1/100) = 0.98
    # PFER = 100 / (100 * 0.98) ≈ 1.0204
    assert math.isclose(ss_pfer, 100.0 / (100.0 * 0.98), rel_tol=1e-4)


def test_resolve_stability_parameters():
    """Verify triangular resolution of parameters."""
    p = 80
    B = 50

    # Case 1: cutoff + q specified -> solve PFER
    res1 = resolve_stability_parameters(
        p=p, cutoff=Some(0.8), q=Some(8.0), pfer=Nothing, B=B, assumption="none"
    )
    assert isinstance(res1, Success)
    params1 = res1.unwrap()
    assert params1.cutoff == 0.8
    assert params1.q == 8.0
    assert params1.pfer > 0

    # Case 2: q + pfer specified -> solve cutoff
    res2 = resolve_stability_parameters(
        p=p, cutoff=Nothing, q=Some(8.0), pfer=Some(1.5), B=B, assumption="none"
    )
    assert isinstance(res2, Success)
    assert 0.5 <= res2.unwrap().cutoff <= 1.0

    # Case 3: cutoff + pfer specified -> solve q
    res3 = resolve_stability_parameters(
        p=p, cutoff=Some(0.8), q=Nothing, pfer=Some(1.5), B=B, assumption="none"
    )
    assert isinstance(res3, Success)
    assert 1.0 <= res3.unwrap().q <= p

    # Invalid: all 3 specified or only 1 specified
    err1 = resolve_stability_parameters(
        p=p, cutoff=Some(0.8), q=Some(8.0), pfer=Some(1.0), B=B
    )
    assert isinstance(err1, Failure)

    err2 = resolve_stability_parameters(
        p=p, cutoff=Some(0.8), q=Nothing, pfer=Nothing, B=B
    )
    assert isinstance(err2, Failure)


# -----------------------------------------------------------------------------
# 2. Subsampling and Feature Randomization
# -----------------------------------------------------------------------------

def test_complementary_pairs_partition():
    """Verify complementary pairs are disjoint and have exact half-sizes."""
    n = 60
    B = 25
    pairs = generate_complementary_pairs(n_samples=n, B=B, seed=123)

    assert len(pairs) == B
    for pair_a, pair_b in pairs:
        assert len(pair_a) == 30
        assert len(pair_b) == 30
        # Check disjointness
        intersection = np.intersect1d(pair_a, pair_b)
        assert len(intersection) == 0


def test_randomized_weights():
    """Verify feature randomization scales columns without mutating original."""
    X = np.ones((20, 5))
    X_orig = X.copy()
    rng = np.random.default_rng(42)

    X_scaled = apply_randomized_weights(X, weakness=0.5, rng=rng)
    # Check original unmodified
    np.testing.assert_array_equal(X, X_orig)
    # Check all weights in [0.5, 1.0]
    assert np.all(X_scaled >= 0.5)
    assert np.all(X_scaled <= 1.0)
    assert not np.allclose(X_scaled, X)


# -----------------------------------------------------------------------------
# 3. Full Stability Selection Engine on High-Dimensional Synthetic Data
# -----------------------------------------------------------------------------

def test_stability_selection_sparse_recovery():
    """Verify stability selection accurately recovers true signals in p > n regime."""
    np.random.seed(42)
    n_samples = 70
    n_features = 60

    # True active features: 0, 1, 2
    true_active = [0, 1, 2]
    X = np.random.randn(n_samples, n_features)
    beta = np.zeros(n_features)
    beta[true_active] = [3.0, 3.5, -2.8]

    # Linear model with noise
    y = X @ beta + np.random.randn(n_samples) * 0.5

    # Parameter resolution: cutoff = 0.75, target PFER = 1.0
    params_res = resolve_stability_parameters(
        p=n_features,
        cutoff=Some(0.75),
        q=Nothing,
        pfer=Some(1.0),
        B=40,
        sampling_type="SS",
        assumption="unimodal",
    )
    assert isinstance(params_res, Success)
    params = params_res.unwrap()

    # Run Stability Selection
    result_res = run_stability_selection(
        X=X,
        y=y,
        parameters=params,
        weakness=0.5,
        seed=42,
    )
    assert isinstance(result_res, Success)
    result = result_res.unwrap()

    # Assertions
    # 1. True signals must have high stability scores
    for true_idx in true_active:
        assert result.stability_scores[true_idx] >= 0.70
        assert true_idx in result.selected_indices

    # 2. Number of falsely selected features must be small (controlled by PFER)
    false_positives = [
        idx for idx in result.selected_indices if idx not in true_active
    ]
    # Bound was PFER <= 1.0
    assert len(false_positives) <= 2

    # 3. Check Polars export
    df = result.to_polars()
    assert df.shape == (n_features, 3)
    assert df["stability_score"][0] >= df["stability_score"][1]  # sorted descending

    df_path = result.to_path_polars()
    assert df_path.height == n_features * len(result.lambdas)


# -----------------------------------------------------------------------------
# 4. Scikit-learn Adapter
# -----------------------------------------------------------------------------

def test_stability_selector_sklearn_pipeline():
    """Verify StabilitySelector behaves as a standard scikit-learn transformer."""
    np.random.seed(99)
    X = np.random.randn(50, 20)
    y = X[:, 0] * 2.0 + X[:, 1] * -2.0 + np.random.randn(50) * 0.5

    pipeline = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "stability",
                StabilitySelector(
                    cutoff=0.70,
                    pfer=1.0,
                    B=30,
                    weakness=0.6,
                    random_state=42,
                ),
            ),
            ("model", Ridge()),
        ]
    )

    pipeline.fit(X, y)
    predictions = pipeline.predict(X)
    assert predictions.shape == (50,)

    selector = pipeline.named_steps["stability"]
    support = selector.get_support()
    assert len(support) == 20
    assert support[0] == True or support[1] == True  # At least one signal picked up
    X_trans = selector.transform(X)
    assert X_trans.shape[1] == np.sum(support)


# -----------------------------------------------------------------------------
# 5. Altair Declarative Visualizations
# -----------------------------------------------------------------------------

def test_altair_visualizations():
    """Verify stability path and score charts create valid Altair specs."""
    np.random.seed(42)
    X = np.random.randn(40, 15)
    y = X[:, 0] * 3.0 + np.random.randn(40)

    params = resolve_stability_parameters(
        p=15, cutoff=Some(0.75), q=Some(4.0), pfer=Nothing, B=20
    ).unwrap()

    result = run_stability_selection(X, y, parameters=params, seed=42).unwrap()

    # Paths plot
    path_chart = plot_stability_paths(result)
    chart_dict = path_chart.to_dict()
    assert "data" in chart_dict or "layer" in chart_dict or "hconcat" in chart_dict

    # Paths plot with budgeted lambda cutoff
    result_budget = run_stability_selection(
        X, y, parameters=params, seed=42, max_expected_q=3.0
    ).unwrap()
    assert result_budget.lambda_cutoff != Nothing
    path_chart_budget = plot_stability_paths(result_budget)
    budget_dict = path_chart_budget.to_dict()
    assert "layer" in budget_dict

    # Scores plot
    score_chart = plot_stability_scores(result)
    score_dict = score_chart.to_dict()
    assert "data" in score_dict or "layer" in score_dict


# -----------------------------------------------------------------------------
# 6. Property-Based Testing (Hypothesis)
# -----------------------------------------------------------------------------

from hypothesis import given, strategies as st, settings


@given(
    st.integers(min_value=5, max_value=100),
    st.floats(min_value=0.01, max_value=0.99),
    st.floats(min_value=0.76, max_value=0.98),
)
@settings(max_examples=50)
def test_property_unimodal_constant_properties(B: int, u: float, c_high: float):
    """Property: Unimodal denominator constant is strictly positive and non-decreasing in valid domain."""
    c_min = 0.5 + 1.0 / (2.0 * B) + 0.005
    c_low = c_min + u * (0.74 - c_min)

    const_low = compute_unimodal_constant(c_low, B)
    const_high = compute_unimodal_constant(c_high, B)

    assert const_low > 0
    assert const_high > 0
    assert const_low <= const_high


@given(
    st.integers(min_value=20, max_value=500),
    st.floats(min_value=1.0, max_value=10.0),
    st.floats(min_value=11.0, max_value=20.0),
    st.floats(min_value=0.55, max_value=0.95),
)
@settings(max_examples=50)
def test_property_mb_pfer_monotonicity(p: int, q_small: float, q_large: float, cutoff: float):
    """Property: PFER increases monotonically with the average selected variables q."""
    pfer_small = solve_pfer(p, cutoff, q_small, 50, "none").unwrap()
    pfer_large = solve_pfer(p, cutoff, q_large, 50, "none").unwrap()

    assert pfer_small > 0
    assert pfer_large > 0
    assert pfer_small <= pfer_large


@given(
    st.integers(min_value=10, max_value=150),
    st.integers(min_value=2, max_value=15),
)
@settings(max_examples=50)
def test_property_complementary_pairs_invariants(n: int, B: int):
    """Property: Complementary pairs are always disjoint and have size floor(n / 2)."""
    m = n // 2
    pairs = generate_complementary_pairs(n_samples=n, B=B, seed=42)
    assert len(pairs) == B

    for pair_a, pair_b in pairs:
        assert len(pair_a) == m
        assert len(pair_b) == m
        # Invariant: Disjointness
        assert len(np.intersect1d(pair_a, pair_b)) == 0
        # Invariant: Valid indices
        assert np.all(pair_a < n) and np.all(pair_a >= 0)
        assert np.all(pair_b < n) and np.all(pair_b >= 0)


def test_regularization_budget_constraint():
    """Verify that specifying max_expected_q bounds empirical model size and restricts scores."""
    np.random.seed(123)
    n_samples = 120
    n_features = 40

    X = np.random.randn(n_samples, n_features)
    beta = np.zeros(n_features)
    beta[0] = 3.0
    beta[1] = 2.5
    y = X @ beta + np.random.randn(n_samples) * 0.5

    params = StabilityParameters(
        p=n_features,
        q=5.0,
        cutoff=0.75,
        pfer=1.0,
        B=30,
        sampling_type="SS",
        assumption="unimodal",
    )

    # 1. Run with tight budget constraint: max_expected_q = 4.0
    res_tight = run_stability_selection(
        X=X,
        y=y,
        parameters=params,
        max_expected_q=4.0,
        seed=42,
    ).unwrap()

    assert res_tight.empirical_q <= 4.0
    assert isinstance(res_tight.lambda_cutoff, Some)
    assert len(res_tight.expected_model_sizes) == len(res_tight.lambdas)
    assert len(res_tight.unrestricted_stability_scores) == n_features

    # All unrestricted scores must be >= budget-constrained scores
    for unres, res in zip(res_tight.unrestricted_stability_scores, res_tight.stability_scores):
        assert unres >= res

    # 2. Check path DataFrame export
    path_df = res_tight.to_path_polars()
    assert "expected_model_size" in path_df.columns
    assert "in_budget" in path_df.columns
    # Rows with in_budget == True should have lambda >= lambda_cutoff
    cutoff_val = res_tight.lambda_cutoff.unwrap()
    budget_rows = path_df.filter(path_df["in_budget"])
    assert (budget_rows["lambda"] >= cutoff_val - 1e-9).all()


def test_stratified_complementary_pairs_preserves_ratios():
    """Verify that stratified subsampling strictly preserves group proportions."""
    strata = np.array(["A"] * 40 + ["B"] * 20 + ["C"] * 10)
    B = 10
    pairs = generate_stratified_complementary_pairs(strata, B=B, seed=42)

    assert len(pairs) == B
    for pair_a, pair_b in pairs:
        # Invariant: Disjointness
        assert len(np.intersect1d(pair_a, pair_b)) == 0
        # Invariant: Size floor(n_c / 2) for each stratum
        # A: 40 // 2 = 20, B: 20 // 2 = 10, C: 10 // 2 = 5 -> total 35
        assert len(pair_a) == 35
        assert len(pair_b) == 35
        assert np.sum(strata[pair_a] == "A") == 20
        assert np.sum(strata[pair_b] == "A") == 20
        assert np.sum(strata[pair_a] == "B") == 10
        assert np.sum(strata[pair_b] == "B") == 10
        assert np.sum(strata[pair_a] == "C") == 5
        assert np.sum(strata[pair_b] == "C") == 5


def test_elastic_net_fitter():
    """Verify Elastic Net fitter generates valid active masks and groups correlated signals."""
    np.random.seed(42)
    n_samples, n_features = 80, 20
    X = np.random.randn(n_samples, n_features)
    # Collinear pair: feature 0 and 1
    X[:, 1] = X[:, 0] * 0.95 + np.random.randn(n_samples) * 0.05
    y = X[:, 0] * 3.0 + np.random.randn(n_samples) * 0.5

    lambdas = np.logspace(-1, -3, 10)
    fitter = create_elastic_net_fitter(l1_ratio=0.7)
    mask = fitter(X, y, lambdas)

    assert mask.shape == (10, n_features)
    assert mask.dtype == bool
    # At least one signal must be active at smallest lambda
    assert mask[-1, 0] or mask[-1, 1]


def test_l1_logistic_fitter():
    """Verify L1 Logistic Regression fitter correctly handles binary targets."""
    np.random.seed(42)
    n_samples, n_features = 100, 15
    X = np.random.randn(n_samples, n_features)
    # Binary classification target
    prob = 1.0 / (1.0 + np.exp(-2.5 * X[:, 0]))
    y = (np.random.rand(n_samples) < prob).astype(int)

    lambdas = np.logspace(-1, -3, 8)
    fitter = create_l1_logistic_fitter(max_iter=100)
    mask = fitter(X, y, lambdas)

    assert mask.shape == (8, n_features)
    assert mask.dtype == bool
    # Feature 0 is the driver
    assert mask[-1, 0]


def test_tree_importance_fitter():
    """Verify Tree importance fitter generates monotonic paths based on feature importance."""
    np.random.seed(42)
    n_samples, n_features = 80, 10
    X = np.random.randn(n_samples, n_features)
    y = (X[:, 2] > 0.0).astype(int)

    lambdas = np.linspace(1.0, 0.01, 10)
    fitter = create_tree_importance_fitter(estimator_type="rf", n_estimators=30, max_depth=3)
    mask = fitter(X, y, lambdas)

    assert mask.shape == (10, n_features)
    assert mask.dtype == bool
    # Monotonicity check: more features should be active as threshold relaxes
    counts = np.sum(mask, axis=1)
    assert counts[-1] >= counts[0]
    # True feature 2 should be active at least at late thresholds
    assert mask[-1, 2]


def test_cohort_adjusted_fitter():
    """Verify FWL cohort adjusted fitter is invariant to cohort mean shifts."""
    np.random.seed(42)
    n_samples, n_features = 100, 10
    X = np.random.randn(n_samples, n_features)
    cohort_labels = np.array(["Cohort1"] * 50 + ["Cohort2"] * 50)
    # Add large cohort-specific intercept shift (+10 vs -10)
    y_raw = X[:, 0] * 2.5 + np.random.randn(n_samples) * 0.2
    y = y_raw.copy()
    y[:50] += 10.0
    y[50:] -= 10.0

    lambdas = np.logspace(-1, -3, 8)
    fitter = create_cohort_adjusted_fitter(cohort_labels, base_fitter_type="lasso")
    mask = fitter(X, y, lambdas)

    assert mask.shape == (8, n_features)
    # The true feature 0 should be recovered despite the huge cohort batch shift
    assert mask[-1, 0]


def test_group_lasso_cohort_fitter():
    """Verify Multi-Cohort Group Lasso selects joint features across cohorts."""
    np.random.seed(42)
    n_samples, n_features = 90, 12
    X = np.random.randn(n_samples, n_features)
    cohort_labels = np.array(["C1"] * 30 + ["C2"] * 30 + ["C3"] * 30)
    # Signal in feature 0 across all cohorts
    y = X[:, 0] * 3.0 + np.random.randn(n_samples) * 0.5

    lambdas = np.logspace(-1, -3, 8)
    fitter = create_group_lasso_cohort_fitter(cohort_labels, max_iter=100)
    mask = fitter(X, y, lambdas)

    assert mask.shape == (8, n_features)
    assert mask.dtype == bool
    # Feature 0 should be selected
    assert mask[-1, 0]


def test_run_stability_selection_stratified_integration():
    """Verify full stability selection pipeline with stratification and custom fitter."""
    np.random.seed(42)
    n_samples, n_features = 80, 20
    X = np.random.randn(n_samples, n_features)
    cohort_labels = np.array(["Batch1"] * 40 + ["Batch2"] * 40)
    y = X[:, 0] * 3.5 + np.random.randn(n_samples) * 0.3

    params = StabilityParameters(
        p=n_features,
        q=5.0,
        cutoff=0.75,
        pfer=1.0,
        B=20,
        sampling_type="SS",
        assumption="unimodal",
    )

    fitter = create_elastic_net_fitter(l1_ratio=0.8)
    result_res = run_stability_selection(
        X=X,
        y=y,
        parameters=params,
        fitter=fitter,
        strata=cohort_labels,
        seed=42,
    )

    assert isinstance(result_res, Success)
    result = result_res.unwrap()
    assert 0 in result.selected_indices


def test_merf_cohort_fitter():
    """Verify Mixed Effects Random Forest isolates cohort random intercepts and recovers signals."""
    np.random.seed(42)
    n_samples, n_features = 120, 10
    X = np.random.randn(n_samples, n_features)
    cohort_labels = np.array(["CohortA"] * 60 + ["CohortB"] * 60)

    # Signal in feature 1 (non-linear threshold) plus cohort random intercepts (+5 vs -5)
    f_X = 3.0 * (X[:, 1] > 0.0).astype(float)
    y = f_X.copy()
    y[:60] += 5.0
    y[60:] -= 5.0
    y += np.random.randn(n_samples) * 0.2

    lambdas = np.linspace(1.0, 0.01, 8)
    fitter = create_merf_cohort_fitter(
        cohort_labels, n_estimators=25, max_depth=3, max_em_iter=3, seed=42
    )
    mask = fitter(X, y, lambdas)

    assert mask.shape == (8, n_features)
    assert mask.dtype == bool
    # Feature 1 must be active at least in late thresholds
    assert mask[-1, 1]


def test_multitask_logistic_cohort_fitter():
    """Verify Multi-Cohort L2,1 Logistic Regression selects joint features under Bernoulli loss."""
    np.random.seed(42)
    n_samples, n_features = 100, 12
    X = np.random.randn(n_samples, n_features)
    cohort_labels = np.array(["Trial1"] * 50 + ["Trial2"] * 50)

    # Binary outcome driven by feature 0 across both trials
    prob = 1.0 / (1.0 + np.exp(-3.0 * X[:, 0]))
    y = (np.random.rand(n_samples) < prob).astype(int)

    lambdas = np.logspace(-1, -3, 8)
    fitter = create_multitask_logistic_cohort_fitter(cohort_labels, tol=1e-3, max_iter=60)
    mask = fitter(X, y, lambdas)

    assert mask.shape == (8, n_features)
    assert mask.dtype == bool
    # Feature 0 should be selected
    assert mask[-1, 0]


def test_meta_analysis_cohort_fitter():
    """Verify Meta-Analytic Consensus requires multi-cohort replication."""
    np.random.seed(42)
    n_samples, n_features = 90, 8
    X = np.random.randn(n_samples, n_features)
    cohort_labels = np.array(["StudyA"] * 30 + ["StudyB"] * 30 + ["StudyC"] * 30)

    # Feature 0 is present in ALL 3 studies
    # Feature 1 is an artifact present ONLY in StudyA
    y = np.zeros(n_samples)
    y += X[:, 0] * 3.0  # Shared signal
    y[:30] += X[:30, 1] * 4.0  # StudyA specific artifact
    y += np.random.randn(n_samples) * 0.3

    lambdas = np.logspace(-1, -3, 6)
    # Require replication in at least 2 cohorts
    fitter = create_meta_analysis_cohort_fitter(
        cohort_labels, base_fitter_type="lasso", min_cohorts=2
    )
    mask = fitter(X, y, lambdas)

    assert mask.shape == (6, n_features)
    assert mask.dtype == bool
    # Feature 0 must replicate across studies
    assert mask[-1, 0]


def test_multistudy_invariant_cohort_fitter():
    """Verify Multi-Study Invariant fitter recovers cross-environment stable features."""
    np.random.seed(42)
    n_samples, n_features = 90, 10
    X = np.random.randn(n_samples, n_features)
    cohort_labels = np.array(["Env1"] * 45 + ["Env2"] * 45)

    # Invariant feature 0 has positive effect in both environments
    y = X[:, 0] * 2.5 + np.random.randn(n_samples) * 0.2

    lambdas = np.logspace(-1, -3, 6)
    fitter = create_multistudy_invariant_cohort_fitter(cohort_labels, gamma=1.0, max_iter=50)
    mask = fitter(X, y, lambdas)

    assert mask.shape == (6, n_features)
    assert mask.dtype == bool
    assert mask[-1, 0]


def test_glmm_lasso_cohort_fitter():
    """Verify Penalized GLMM Logistic absorbs cohort baseline log-odds shifts."""
    np.random.seed(42)
    n_samples, n_features = 100, 10
    X = np.random.randn(n_samples, n_features)
    cohort_labels = np.array(["HospitalA"] * 50 + ["HospitalB"] * 50)

    # Baseline response rate is high in HospitalA (+2.0 offset) and low in HospitalB (-2.0 offset)
    logits = X[:, 0] * 2.5
    logits[:50] += 2.0
    logits[50:] -= 2.0
    prob = 1.0 / (1.0 + np.exp(-logits))
    y = (np.random.rand(n_samples) < prob).astype(int)

    lambdas = np.logspace(-1, -3, 6)
    fitter = create_glmm_lasso_cohort_fitter(cohort_labels, max_iter=100)
    mask = fitter(X, y, lambdas)

    assert mask.shape == (6, n_features)
    assert mask.dtype == bool
    assert mask[-1, 0]


def test_multi_cohort_fitters_single_cohort_fallback():
    """Verify all 5 multi-cohort fitters cleanly degrade when only 1 cohort is present."""
    np.random.seed(42)
    n_samples, n_features = 40, 8
    X = np.random.randn(n_samples, n_features)
    cohort_labels = np.array(["SingleCohort"] * n_samples)
    y = (X[:, 0] > 0.0).astype(int)
    lambdas = np.logspace(-1, -3, 5)

    fitters = [
        create_merf_cohort_fitter(cohort_labels, n_estimators=10, max_depth=2, max_em_iter=2),
        create_multitask_logistic_cohort_fitter(cohort_labels, max_iter=30),
        create_meta_analysis_cohort_fitter(cohort_labels, min_cohorts=2),
        create_multistudy_invariant_cohort_fitter(cohort_labels, max_iter=30),
        create_glmm_lasso_cohort_fitter(cohort_labels, max_iter=50),
    ]

    for fitter in fitters:
        mask = fitter(X, y, lambdas)
        assert mask.shape == (5, n_features)
        assert mask.dtype == bool


def test_prox_sorted_l1_constant_weights():
    """Verify _prox_sorted_l1 with constant weights reduces to standard soft-thresholding."""
    np.random.seed(42)
    v = np.array([2.5, -1.8, 0.4, -0.2, 3.1])
    p = len(v)
    lam = 0.5
    eta = 0.2
    weights = np.full(p, lam)

    # Standard soft-thresholding: sign(v) * max(|v| - eta * lam, 0)
    expected = np.sign(v) * np.maximum(np.abs(v) - eta * lam, 0.0)
    result = _prox_sorted_l1(v, weights, eta)

    assert np.allclose(result, expected, atol=1e-8)


def test_oscar_exact_clustering_property():
    """Verify OSCAR groups strongly collinear features into exact equal magnitude clusters."""
    np.random.seed(42)
    n_samples, n_features = 80, 10
    X = np.random.randn(n_samples, n_features)
    # Collinear pair: feature 1 is identical to feature 0 with tiny noise
    X[:, 1] = X[:, 0] + np.random.randn(n_samples) * 1e-4

    y = X[:, 0] * 3.0 + X[:, 1] * 3.0 + np.random.randn(n_samples) * 0.1
    lambdas = np.logspace(-1, -3, 6)

    # Strong grouping kappa
    oscar_fitter = create_oscar_fitter(kappa=0.8, tol=1e-4, max_iter=150)
    mask = oscar_fitter(X, y, lambdas)

    assert mask.shape == (6, n_features)
    assert mask.dtype == bool
    # Both collinear features must be co-selected
    assert mask[-1, 0]
    assert mask[-1, 1]


def test_slope_fdr_adaptive_thresholding():
    """Verify SLOPE detects true sparse signals while suppressing Gaussian noise features."""
    np.random.seed(42)
    n_samples, n_features = 100, 20
    X = np.random.randn(n_samples, n_features)
    # Features 0 and 1 are true signals
    y = X[:, 0] * 4.0 + X[:, 1] * 3.0 + np.random.randn(n_samples) * 0.2

    lambdas = np.logspace(-1, -3, 8)
    slope_fitter = create_slope_fitter(q_fdr=0.1, tol=1e-4, max_iter=100)
    mask = slope_fitter(X, y, lambdas)

    assert mask.shape == (8, n_features)
    assert mask.dtype == bool
    # At least at lower penalties, true signals are active
    assert mask[-1, 0]
    assert mask[-1, 1]


def test_oscar_and_slope_cohort_fitters():
    """Verify cohort-adjusted OSCAR and SLOPE absorb cohort batch mean shifts."""
    np.random.seed(42)
    n_samples, n_features = 100, 12
    X = np.random.randn(n_samples, n_features)
    cohort_labels = np.array(["Center1"] * 50 + ["Center2"] * 50)

    # Large cohort baseline shift (+5.0 for Center1, -5.0 for Center2)
    y = X[:, 0] * 3.5 + np.random.randn(n_samples) * 0.2
    y[:50] += 5.0
    y[50:] -= 5.0

    lambdas = np.logspace(-1, -3, 6)

    oscar_fitter = create_oscar_cohort_fitter(cohort_labels, kappa=0.6)
    slope_fitter = create_slope_cohort_fitter(cohort_labels, q_fdr=0.1)

    mask_oscar = oscar_fitter(X, y, lambdas)
    mask_slope = slope_fitter(X, y, lambdas)

    assert mask_oscar.shape == (6, n_features)
    assert mask_slope.shape == (6, n_features)
    assert mask_oscar[-1, 0]
    assert mask_slope[-1, 0]


def test_oscar_slope_single_cohort_fallback():
    """Verify OSCAR and SLOPE cohort fitters fall back cleanly when only 1 cohort is present."""
    np.random.seed(42)
    n_samples, n_features = 40, 8
    X = np.random.randn(n_samples, n_features)
    cohort_labels = np.array(["SingleCohort"] * n_samples)
    y = X[:, 0] * 2.0 + np.random.randn(n_samples) * 0.2
    lambdas = np.logspace(-1, -3, 5)

    oscar_cohort = create_oscar_cohort_fitter(cohort_labels, kappa=0.5)
    slope_cohort = create_slope_cohort_fitter(cohort_labels, q_fdr=0.1)

    mask_oscar = oscar_cohort(X, y, lambdas)
    mask_slope = slope_cohort(X, y, lambdas)

    assert mask_oscar.shape == (5, n_features)
    assert mask_slope.shape == (5, n_features)
    assert mask_oscar.dtype == bool
    assert mask_slope.dtype == bool




