"""Pure functional Stability Selection computation engine.

Implements the complete algorithm:
1. Subsampling / Complementary Pairs
2. Column-wise randomized feature weighting (Randomized Lasso)
3. Regularization path fitting
4. Selection probability matrix and stability score aggregation
5. Thresholding with finite-sample error bound guarantees
"""

import warnings
from typing import Sequence, cast
import numpy as np
from joblib import Parallel, delayed  # type: ignore[import-untyped]
from sklearn.exceptions import ConvergenceWarning  # type: ignore[import-untyped]
from sklearn.linear_model import lasso_path  # type: ignore[import-untyped]
from returns.result import Result, Success, Failure
from returns.maybe import Maybe, Some, Nothing

from selective_inference.stability_selection.types import (
    PathFitter,
    StabilityParameters,
    StabilityResult,
)
from selective_inference.stability_selection.subsampling import (
    generate_complementary_pairs,
    generate_subsamples,
    generate_stratified_complementary_pairs,
    generate_stratified_subsamples,
    apply_randomized_weights,
)


def default_lasso_path_fitter(
    X: np.ndarray,
    y: np.ndarray,
    lambdas: Sequence[float] | np.ndarray,
    subsample_indices: np.ndarray | None = None,
) -> np.ndarray:
    """Default high-performance Lasso path fitter using coordinate descent.

    Parameters
    ----------
    X : np.ndarray of shape (n_samples, n_features)
    y : np.ndarray of shape (n_samples,)
    lambdas : Sequence of penalty values (sorted descending)

    Returns
    -------
    np.ndarray of shape (n_lambdas, n_features)
        Boolean mask of selected features (nonzero coefficients).
    """
    alphas = np.asarray(lambdas, dtype=np.float64)
    # lasso_path returns coefs of shape (n_features, n_alphas)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        _, coef_path, _ = lasso_path(X, y, alphas=alphas)
    # Transpose to (n_alphas, n_features) and check non-zero
    return cast(np.ndarray, (np.abs(coef_path.T) > 1e-6).astype(bool))


def generate_default_lambdas(
    X: np.ndarray, y: np.ndarray, n_lambdas: int = 25
) -> np.ndarray:
    """Generate a geometric grid of penalty parameters lambda.

    Starts at lambda_max = ||X^T y||_inf / n and descends logarithmically.
    """
    n = X.shape[0]
    # Center y and standardize columns to determine lambda_max accurately
    y_centered = y - np.mean(y)
    corrs = np.abs(X.T @ y_centered)
    lambda_max = float(np.max(corrs) / max(n, 1))
    if lambda_max <= 1e-8:
        lambda_max = 1.0
    lambda_min = lambda_max * 0.01
    return np.logspace(np.log10(lambda_max), np.log10(lambda_min), n_lambdas)


def run_stability_selection(
    X: np.ndarray,
    y: np.ndarray,
    parameters: StabilityParameters,
    lambdas: Sequence[float] | None = None,
    fitter: PathFitter = default_lasso_path_fitter,
    weakness: float = 0.5,
    feature_names: Sequence[str] | None = None,
    seed: int = 42,
    max_expected_q: float | None = None,
    strata: Sequence[object] | np.ndarray | None = None,
    n_jobs: int = 1,
) -> Result[StabilityResult, str]:
    """Execute Stability Selection on input data.

    Pure function: no mutation of inputs, deterministic execution given seed,
    returns immutable StabilityResult wrapped in a Result monad.

    Parameters
    ----------
    X : np.ndarray of shape (n_samples, n_features)
    y : np.ndarray of shape (n_samples,)
    parameters : StabilityParameters
        Resolved stability selection parameters including cutoff and B.
    lambdas : Sequence[float], optional
        Grid of regularization penalties. If omitted, generated automatically.
    fitter : PathFitter, default default_lasso_path_fitter
        Function to compute active sets along the penalty path.
    weakness : float, default 0.5
        Feature randomization factor alpha in (0, 1].
    feature_names : Sequence[str], optional
        Names of features. Defaults to ("x0", "x1", ...).
    seed : int, default 42
        Random seed for subsampling and weight generation.
    max_expected_q : float, optional
        Maximum allowed expected model size E[|S(lambda)|] along the path.
        If omitted, defaults to parameters.q. Penalties with larger expected
        model sizes are excluded from the stability score calculation to preserve
        finite-sample error bounds.
    strata : Sequence or np.ndarray, optional
        Stratification labels (e.g. cohort IDs) of length n_samples. When specified,
        subsampling partitions observations proportionally within each stratum.

    Returns
    -------
    Result[StabilityResult, str]
    """
    if X.ndim != 2:
        return Failure(f"X must be a 2D matrix, got shape {X.shape}")
    if y.ndim != 1:
        return Failure(f"y must be a 1D vector, got shape {y.shape}")

    n_samples, n_features = X.shape
    if n_samples != y.shape[0]:
        return Failure(
            f"Sample size mismatch: X has {n_samples}, y has {y.shape[0]}"
        )
    if n_features != parameters.p:
        return Failure(
            f"Feature count mismatch: X has {n_features}, parameters specify {parameters.p}"
        )

    # Resolve feature names
    names = (
        tuple(feature_names)
        if feature_names is not None
        else tuple(f"x_{i}" for i in range(n_features))
    )
    if len(names) != n_features:
        return Failure(
            f"Length of feature_names ({len(names)}) does not match n_features ({n_features})"
        )

    # Resolve regularization grid
    path_lambdas = (
        np.asarray(lambdas, dtype=np.float64)
        if lambdas is not None
        else generate_default_lambdas(X, y)
    )
    n_lambdas = len(path_lambdas)

    # Generate subsamples
    if strata is not None:
        strata_arr = np.asarray(strata)
        if len(strata_arr) != n_samples:
            return Failure(
                f"Length of strata ({len(strata_arr)}) must match n_samples ({n_samples})"
            )
        if parameters.sampling_type == "SS":
            strat_pairs = generate_stratified_complementary_pairs(strata_arr, parameters.B, seed=seed)
            subsample_indices: list[np.ndarray] = []
            for pair_a, pair_b in strat_pairs:
                subsample_indices.append(pair_a)
                subsample_indices.append(pair_b)
        else:
            strat_subsamples = generate_stratified_subsamples(strata_arr, parameters.B, seed=seed)
            subsample_indices = list(strat_subsamples)
    elif parameters.sampling_type == "SS":
        pairs = generate_complementary_pairs(n_samples, parameters.B, seed=seed)
        subsample_indices = []
        for pair_a, pair_b in pairs:
            subsample_indices.append(pair_a)
            subsample_indices.append(pair_b)
    else:
        subsamples = generate_subsamples(n_samples, parameters.B, seed=seed)
        subsample_indices = list(subsamples)

    total_subsamples = len(subsample_indices)

    def _fit_single_subsample(
        idx: int,
        sub_idx: np.ndarray,
    ) -> tuple[np.ndarray, float] | None:
        sub_rng = np.random.default_rng(seed + idx)
        X_sub = X[sub_idx]
        y_sub = y[sub_idx]
        X_perturbed = apply_randomized_weights(X_sub, weakness=weakness, rng=sub_rng)
        try:
            try:
                active_mask = fitter(X_perturbed, y_sub, path_lambdas, subsample_indices=sub_idx)
            except TypeError:
                active_mask = fitter(X_perturbed, y_sub, path_lambdas)
        except Exception:
            return None

        if active_mask.shape != (n_lambdas, n_features):
            return None

        return active_mask.astype(np.int32), float(np.mean(np.sum(active_mask, axis=1)))

    if n_jobs == 1 or total_subsamples <= 1:
        fits = [_fit_single_subsample(i, s_idx) for i, s_idx in enumerate(subsample_indices)]
    else:
        fits = Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(_fit_single_subsample)(i, s_idx)
            for i, s_idx in enumerate(subsample_indices)
        )

    # Accumulator for selections: shape (n_lambdas, n_features)
    selection_counts = np.zeros((n_lambdas, n_features), dtype=np.int32)
    selected_per_subsample: list[float] = []

    for fit_res in fits:
        if fit_res is None:
            continue
        mask, mean_k = fit_res
        selection_counts += mask
        selected_per_subsample.append(mean_k)

    if not selected_per_subsample:
        return Failure("All subsample fits failed in stability selection.")

    # Selection probability matrix Pi: shape (n_features, n_lambdas)
    pi_matrix = (selection_counts / total_subsamples).T  # Transpose to (n_features, n_lambdas)

    # Empirical expected model size at each penalty lambda_j: sum of selection probabilities across features
    expected_model_sizes = np.sum(pi_matrix, axis=0)  # shape: (n_lambdas,)

    # Determine valid lambda mask based on budget q
    target_q = float(max_expected_q) if max_expected_q is not None else float(parameters.q)
    valid_mask = expected_model_sizes <= target_q

    # Safeguard: if even lambda_max exceeds target_q, retain at least the first penalty
    if not np.any(valid_mask):
        valid_mask[0] = True

    # Identify the cutoff lambda (the smallest lambda satisfying the budget)
    valid_indices = np.where(valid_mask)[0]
    cutoff_idx = int(valid_indices[-1])
    lambda_cutoff: Maybe[float] = Some(float(path_lambdas[cutoff_idx]))

    # Compute budget-constrained stability scores and unrestricted full-path scores
    stability_scores = np.max(pi_matrix[:, valid_mask], axis=1)
    unrestricted_scores = np.max(pi_matrix, axis=1)

    # Selected features: stability_score >= cutoff
    selected_mask = stability_scores >= parameters.cutoff
    selected_indices = tuple(int(idx) for idx in np.where(selected_mask)[0])
    selected_features = tuple(names[idx] for idx in selected_indices)

    # Empirical q at the regularization cutoff
    empirical_q = float(expected_model_sizes[cutoff_idx])

    stability_matrix_tuple = tuple(
        tuple(float(val) for val in row) for row in pi_matrix
    )

    return Success(
        StabilityResult(
            feature_names=names,
            lambdas=tuple(float(l) for l in path_lambdas),
            stability_matrix=stability_matrix_tuple,
            stability_scores=tuple(float(score) for score in stability_scores),
            selected_features=selected_features,
            selected_indices=selected_indices,
            empirical_q=empirical_q,
            parameters=parameters,
            expected_model_sizes=tuple(float(sz) for sz in expected_model_sizes),
            lambda_cutoff=lambda_cutoff,
            unrestricted_stability_scores=tuple(float(score) for score in unrestricted_scores),
        )
    )
