"""Modular high-dimensional path fitters for Stability Selection.

All fitters implement the PathFitter Protocol:
(X: np.ndarray, y: np.ndarray, lambdas: Sequence[float] | np.ndarray) -> np.ndarray[bool]
returning an active feature selection mask of shape (n_lambdas, n_features).

Fitters included:
1. create_lasso_fitter: Coordinate descent L1 penalized linear regression.
2. create_elastic_net_fitter: L1 + L2 coordinate descent grouping correlated genes.
3. create_l1_logistic_fitter: Bernoulli log-likelihood for binary response.
4. create_tree_importance_fitter: Non-linear Random Forest feature importance paths.
5. create_cohort_adjusted_fitter: Unpenalized cohort fixed effects via FWL projection.
6. create_group_lasso_cohort_fitter: Multi-cohort L2,1 joint feature selection.
7. create_merf_cohort_fitter: Mixed Effects Random Forest (EM random intercepts + RF).
8. create_multitask_logistic_cohort_fitter: Multi-cohort L2,1 sparse logistic regression.
9. create_meta_analysis_cohort_fitter: Meta-analytic consensus path fitter over cohorts.
10. create_multistudy_invariant_cohort_fitter: Invariant multi-study shrinkage fitter.
11. create_glmm_lasso_cohort_fitter: High-dimensional penalized GLMM logistic regression.
12. create_oscar_fitter: Octagonal Shrinkage and Clustering Algorithm for Regression.
13. create_slope_fitter: Sorted L-One Penalized Estimation (BH adaptive FDR control).
14. create_oscar_cohort_fitter: FWL cohort-adjusted OSCAR with clustering.
15. create_slope_cohort_fitter: FWL cohort-adjusted SLOPE with FDR control.
"""

from typing import Callable, Literal, Sequence, cast
import numpy as np
from scipy.stats import norm  # type: ignore[import-untyped]
from sklearn.isotonic import isotonic_regression  # type: ignore[import-untyped]
from sklearn.linear_model import enet_path, lasso_path, LogisticRegression  # type: ignore[import-untyped]
from sklearn.ensemble import (  # type: ignore[import-untyped]
    ExtraTreesClassifier,
    RandomForestClassifier,
    RandomForestRegressor,
)

from selective_inference.stability_selection.types import PathFitter


def default_lasso_path_fitter(
    X: np.ndarray,
    y: np.ndarray,
    lambdas: Sequence[float] | np.ndarray,
    subsample_indices: np.ndarray | None = None,
) -> np.ndarray:
    """Default Lasso path fitter using coordinate descent."""
    alphas = np.asarray(lambdas, dtype=np.float64)
    _, coef_path, _ = lasso_path(X, y, alphas=alphas)
    return cast(np.ndarray, (np.abs(coef_path.T) > 1e-6).astype(bool))


def create_lasso_fitter() -> PathFitter:
    """Factory returning the standard coordinate descent Lasso path fitter."""
    return default_lasso_path_fitter


def create_elastic_net_fitter(l1_ratio: float = 0.7) -> PathFitter:
    """Create an Elastic Net path fitter (L1 + L2 regularization).

    Groups co-regulated and collinear features by adding an L2 quadratic penalty.

    Parameters
    ----------
    l1_ratio : float, default 0.7
        Mixing parameter between L1 and L2 penalty. 1.0 is pure Lasso,
        lower values increase the grouping effect.

    Returns
    -------
    PathFitter
    """
    ratio = float(np.clip(l1_ratio, 0.01, 1.0))

    def fitter(
        X: np.ndarray,
        y: np.ndarray,
        lambdas: Sequence[float] | np.ndarray,
        subsample_indices: np.ndarray | None = None,
    ) -> np.ndarray:
        alphas = np.asarray(lambdas, dtype=np.float64)
        _, coef_path, _ = enet_path(X, y, alphas=alphas, l1_ratio=ratio)
        return cast(np.ndarray, (np.abs(coef_path.T) > 1e-6).astype(bool))

    return fitter


def create_l1_logistic_fitter(
    tol: float = 1e-3, max_iter: int = 200
) -> PathFitter:
    """Create an L1-penalized Logistic Regression path fitter for binary outcomes.

    Optimizes the Bernoulli log-likelihood along the regularization path.
    Maps penalty lambda to inverse regularization C = 1.0 / (lambda * n_samples).

    Parameters
    ----------
    tol : float, default 1e-3
        Convergence tolerance for the SAGA coordinate descent solver.
    max_iter : int, default 200
        Maximum iterations per penalty.

    Returns
    -------
    PathFitter
    """
    def fitter(
        X: np.ndarray,
        y: np.ndarray,
        lambdas: Sequence[float] | np.ndarray,
        subsample_indices: np.ndarray | None = None,
    ) -> np.ndarray:
        n_samples, n_features = X.shape
        alphas = np.asarray(lambdas, dtype=np.float64)
        # SAGA with warm start along descending regularization (ascending C)
        Cs = 1.0 / np.maximum(alphas * max(n_samples, 1), 1e-8)

        # Sort ascending C for stable warm starts
        sort_order = np.argsort(Cs)
        rev_order = np.argsort(sort_order)
        sorted_Cs = Cs[sort_order]

        clf = LogisticRegression(
            solver="saga",
            l1_ratio=1.0,
            warm_start=True,
            tol=tol,
            max_iter=max_iter,
            random_state=42,
        )

        coef_list: list[np.ndarray] = []
        for c_val in sorted_Cs:
            clf.C = float(c_val)
            clf.fit(X, y)
            coef_list.append(np.squeeze(clf.coef_).copy())

        # Reorder back to original lambda order
        sorted_coefs = np.asarray(coef_list, dtype=np.float64)
        original_coefs = sorted_coefs[rev_order]
        return (np.abs(original_coefs) > 1e-5).astype(bool)

    return fitter


def create_tree_importance_fitter(
    estimator_type: Literal["rf", "extra_trees"] = "rf",
    n_estimators: int = 50,
    max_depth: int | None = 4,
    seed: int = 42,
) -> PathFitter:
    """Create a Tree-based feature importance path fitter.

    Fits an ensemble of decision trees and evaluates feature inclusion
    across a geometric sequence of importance thresholds.

    Parameters
    ----------
    estimator_type : {"rf", "extra_trees"}, default "rf"
        Ensemble type: RandomForestClassifier or ExtraTreesClassifier.
    n_estimators : int, default 50
        Number of trees in ensemble.
    max_depth : int or None, default 4
        Maximum depth to prevent overfitting on subsamples.
    seed : int, default 42
        Random state for reproducibility.

    Returns
    -------
    PathFitter
    """
    def fitter(
        X: np.ndarray,
        y: np.ndarray,
        lambdas: Sequence[float] | np.ndarray,
        subsample_indices: np.ndarray | None = None,
    ) -> np.ndarray:
        n_lambdas = len(lambdas)
        n_features = X.shape[1]

        cls = RandomForestClassifier if estimator_type == "rf" else ExtraTreesClassifier
        model = cls(
            n_estimators=n_estimators,
            max_depth=max_depth,
            random_state=seed,
            n_jobs=1,
        )
        model.fit(X, y)
        importances = np.asarray(model.feature_importances_, dtype=np.float64)
        max_imp = float(np.max(importances)) if len(importances) > 0 else 1.0
        if max_imp <= 1e-8:
            return np.zeros((n_lambdas, n_features), dtype=bool)

        rel_importances = importances / max_imp

        # Geometric threshold grid from 0.85 down to 0.02
        thresholds = np.logspace(np.log10(0.85), np.log10(0.02), n_lambdas)
        # Shape: (n_lambdas, n_features)
        active_mask = rel_importances[np.newaxis, :] >= thresholds[:, np.newaxis]
        return active_mask.astype(bool)

    return fitter


def create_cohort_adjusted_fitter(
    cohort_labels: np.ndarray | Sequence[object],
    base_fitter_type: Literal["lasso", "elastic_net"] = "lasso",
    l1_ratio: float = 0.7,
) -> PathFitter:
    """Create an unpenalized cohort fixed effects fitter using Frisch-Waugh-Lovell projection.

    Absorbs cohort-specific baseline response rates and batch mean shifts
    without consuming any gene sparsity budget.

    Parameters
    ----------
    cohort_labels : np.ndarray or Sequence
        Array of length n_samples containing cohort identifiers for each observation.
    base_fitter_type : {"lasso", "elastic_net"}, default "lasso"
        Base shrinkage model applied after cohort residualization.
    l1_ratio : float, default 0.7
        Elastic net ratio if base_fitter_type == "elastic_net".

    Returns
    -------
    PathFitter
    """
    labels_arr = np.asarray(cohort_labels)
    unique_cohorts = np.unique(labels_arr)

    def fitter(
        X: np.ndarray,
        y: np.ndarray,
        lambdas: Sequence[float] | np.ndarray,
        subsample_indices: np.ndarray | None = None,
    ) -> np.ndarray:
        n_samples = X.shape[0]
        if subsample_indices is not None and len(subsample_indices) == n_samples:
            sub_labels = labels_arr[subsample_indices]
        else:
            sub_labels = labels_arr[:n_samples] if len(labels_arr) >= n_samples else labels_arr

        # Center X and y within each cohort (FWL orthogonal projection)
        X_res = X.copy()
        y_res = y.astype(np.float64).copy()

        for c in unique_cohorts:
            c_mask = sub_labels == c
            if np.sum(c_mask) > 1:
                X_res[c_mask] = X_res[c_mask] - np.mean(X_res[c_mask], axis=0)
                y_res[c_mask] = y_res[c_mask] - np.mean(y_res[c_mask])

        if base_fitter_type == "elastic_net":
            return create_elastic_net_fitter(l1_ratio=l1_ratio)(
                X_res, y_res, lambdas, subsample_indices=subsample_indices
            )
        return default_lasso_path_fitter(
            X_res, y_res, lambdas, subsample_indices=subsample_indices
        )

    return fitter


def create_group_lasso_cohort_fitter(
    cohort_labels: np.ndarray | Sequence[object],
    tol: float = 1e-3,
    max_iter: int = 50,
) -> PathFitter:
    """Create a Multi-Cohort Group Lasso (L2,1) path fitter.

    Models each cohort as a distinct task while enforcing joint feature sparsity
    across all cohorts:
    min sum_t (1 / (2*n_t)) ||y_t - X_t beta_t||^2 + lambda sum_k ||beta_{k, :}||_2

    A feature k enters with cohort-specific slope adjustments, but is selected
    jointly across all studies using accelerated proximal gradient descent (FISTA).

    Parameters
    ----------
    cohort_labels : np.ndarray or Sequence
        Array of cohort identifiers for each observation.
    tol : float, default 1e-3
        Convergence tolerance for the FISTA solver.
    max_iter : int, default 50
        Maximum FISTA iterations per penalty value.

    Returns
    -------
    PathFitter
    """
    labels_arr = np.asarray(cohort_labels)
    unique_cohorts = tuple(np.unique(labels_arr))

    def fitter(
        X: np.ndarray,
        y: np.ndarray,
        lambdas: Sequence[float] | np.ndarray,
        subsample_indices: np.ndarray | None = None,
    ) -> np.ndarray:
        n_samples, n_features = X.shape
        n_lambdas = len(lambdas)
        if subsample_indices is not None and len(subsample_indices) == n_samples:
            sub_labels = labels_arr[subsample_indices]
        else:
            sub_labels = labels_arr[:n_samples] if len(labels_arr) >= n_samples else labels_arr

        # Group data by cohort
        cohort_masks = [sub_labels == c for c in unique_cohorts]
        # Filter out empty cohorts in this subsample
        active_cohort_indices = [t for t, m in enumerate(cohort_masks) if np.sum(m) > 1]
        T = len(active_cohort_indices)

        if T <= 1:
            # Fall back to standard Lasso if only 1 cohort is present in subsample
            return default_lasso_path_fitter(
                X, y, lambdas, subsample_indices=subsample_indices
            )

        X_t_list = [X[cohort_masks[t]] for t in active_cohort_indices]
        y_t_list = [y[cohort_masks[t]].astype(np.float64) for t in active_cohort_indices]
        n_t_list = [float(len(yt)) for yt in y_t_list]

        # Standardize X_t column norms: sum_i X_{t, i, k}^2 / n_t = 1
        X_norms = [
            np.sqrt(np.mean(Xt ** 2, axis=0) + 1e-8) for Xt in X_t_list
        ]
        X_scaled_list = [
            Xt / norm[np.newaxis, :] for Xt, norm in zip(X_t_list, X_norms)
        ]

        # Compute Lipschitz constant for step size eta
        L_list = [
            float(np.linalg.norm(Xt.T @ Xt, 2) / nt)
            for Xt, nt in zip(X_scaled_list, n_t_list)
        ]
        L_max = max(L_list) if L_list else 1.0
        eta = 1.0 / (L_max + 1e-4)

        # Coefficients: B has shape (n_features, T)
        B = np.zeros((n_features, T), dtype=np.float64)
        active_mask = np.zeros((n_lambdas, n_features), dtype=bool)
        path_lambdas = np.asarray(lambdas, dtype=np.float64)

        for l_idx, lam in enumerate(path_lambdas):
            Y = B.copy()
            for m in range(1, max_iter + 1):
                grad_Y = np.zeros((n_features, T), dtype=np.float64)
                for t in range(T):
                    r_t = y_t_list[t] - X_scaled_list[t] @ Y[:, t]
                    grad_Y[:, t] = - (X_scaled_list[t].T @ r_t) / n_t_list[t]

                tilde_B = Y - eta * grad_Y
                row_norms = np.linalg.norm(tilde_B, axis=1, keepdims=True)
                shrink = np.maximum(0.0, 1.0 - (eta * lam) / np.maximum(row_norms, 1e-12))
                B_new = tilde_B * shrink

                Y = B_new + ((m - 1) / (m + 2)) * (B_new - B)
                diff = float(np.max(np.abs(B_new - B)))
                B = B_new
                if diff < tol:
                    break

            # Feature k is active if its L2 norm across cohorts exceeds 1e-4
            active_mask[l_idx] = np.linalg.norm(B, axis=1) > 1e-4

        return active_mask

    return fitter


def _sigmoid(z: np.ndarray) -> np.ndarray:
    """Numerically stable sigmoid function."""
    z_clipped = np.clip(z, -35.0, 35.0)
    return 1.0 / (1.0 + np.exp(-z_clipped))


def create_merf_cohort_fitter(
    cohort_labels: np.ndarray | Sequence[object],
    n_estimators: int = 30,
    max_depth: int | None = 3,
    max_em_iter: int = 4,
    seed: int = 42,
) -> PathFitter:
    """Create a Mixed Effects Random Forest (MERF) path fitter.

    Isolates cohort-level random intercepts via an Expectation-Maximization (EM) loop:
        y_i = f(X_i) + b_{c_i} + e_i,  b_c ~ N(0, sigma_b^2), e_i ~ N(0, sigma_e^2)
    where f(X) is a non-linear Random Forest capturing fixed feature effects.

    Feature importance paths are extracted from the converged forest f(X)
    and evaluated across a geometric threshold sequence.

    Parameters
    ----------
    cohort_labels : np.ndarray or Sequence
        Array of cohort identifiers for each observation.
    n_estimators : int, default 30
        Number of trees in the Random Forest ensemble.
    max_depth : int or None, default 3
        Maximum tree depth to prevent overfitting on subsamples.
    max_em_iter : int, default 4
        Maximum EM iterations for random effect convergence.
    seed : int, default 42
        Random state for tree building.

    Returns
    -------
    PathFitter
    """
    labels_arr = np.asarray(cohort_labels)
    unique_cohorts = tuple(np.unique(labels_arr))

    def fitter(
        X: np.ndarray,
        y: np.ndarray,
        lambdas: Sequence[float] | np.ndarray,
        subsample_indices: np.ndarray | None = None,
    ) -> np.ndarray:
        n_samples, n_features = X.shape
        n_lambdas = len(lambdas)
        if subsample_indices is not None and len(subsample_indices) == n_samples:
            sub_labels = labels_arr[subsample_indices]
        else:
            sub_labels = labels_arr[:n_samples] if len(labels_arr) >= n_samples else labels_arr

        cohort_masks = {c: sub_labels == c for c in unique_cohorts}
        active_cohorts = tuple(c for c, m in cohort_masks.items() if np.sum(m) > 1)

        if len(active_cohorts) <= 1:
            return create_tree_importance_fitter(
                estimator_type="rf",
                n_estimators=n_estimators,
                max_depth=max_depth,
                seed=seed,
            )(X, y, lambdas, subsample_indices=subsample_indices)

        y_float = y.astype(np.float64)
        var_y = float(np.var(y_float))
        sigma_e2 = max(var_y, 1e-4)
        sigma_b2 = max(var_y * 0.5, 1e-4)
        b_dict = {c: 0.0 for c in active_cohorts}

        rf = RandomForestRegressor(
            n_estimators=n_estimators,
            max_depth=max_depth,
            random_state=seed,
            n_jobs=1,
        )

        for _ in range(max_em_iter):
            y_star = y_float.copy()
            for c in active_cohorts:
                y_star[cohort_masks[c]] -= b_dict[c]

            rf.fit(X, y_star)
            pred_fixed = rf.predict(X)
            residuals = y_float - pred_fixed

            b_sq_sum = 0.0
            for c in active_cohorts:
                c_mask = cohort_masks[c]
                n_c = float(np.sum(c_mask))
                r_bar_c = float(np.mean(residuals[c_mask]))
                shrinkage = (n_c * sigma_b2) / (n_c * sigma_b2 + sigma_e2 + 1e-8)
                b_c = shrinkage * r_bar_c
                b_dict[c] = b_c
                b_sq_sum += b_c ** 2

            e = y_float.copy()
            for c in active_cohorts:
                e[cohort_masks[c]] -= (pred_fixed[cohort_masks[c]] + b_dict[c])
            sigma_e2 = max(1e-4, float(np.mean(e ** 2)))
            sigma_b2 = max(1e-4, b_sq_sum / max(len(active_cohorts), 1))

        importances = np.asarray(rf.feature_importances_, dtype=np.float64)
        max_imp = float(np.max(importances)) if len(importances) > 0 else 1.0
        if max_imp <= 1e-8:
            return np.zeros((n_lambdas, n_features), dtype=bool)

        rel_importances = importances / max_imp
        thresholds = np.logspace(np.log10(0.85), np.log10(0.02), n_lambdas)
        active_mask = rel_importances[np.newaxis, :] >= thresholds[:, np.newaxis]
        return active_mask.astype(bool)

    return fitter


def create_multitask_logistic_cohort_fitter(
    cohort_labels: np.ndarray | Sequence[object],
    tol: float = 1e-3,
    max_iter: int = 60,
) -> PathFitter:
    """Create a Multi-Cohort L2,1 Sparse Logistic Regression path fitter.

    Optimizes exact Bernoulli log-likelihood across multiple cohorts with
    joint group-sparsity:
    min sum_t (1/n_t) sum_{i in c_t} log(1 + exp(- (2 y_i - 1) X_{t, i} beta_t)) + lambda sum_k ||beta_{k, :}||_2

    Solves with accelerated proximal gradient (FISTA) using row-wise group soft-thresholding.

    Parameters
    ----------
    cohort_labels : np.ndarray or Sequence
        Array of cohort identifiers for each observation.
    tol : float, default 1e-3
        Convergence tolerance for the FISTA solver.
    max_iter : int, default 60
        Maximum FISTA iterations per penalty value.

    Returns
    -------
    PathFitter
    """
    labels_arr = np.asarray(cohort_labels)
    unique_cohorts = tuple(np.unique(labels_arr))

    def fitter(
        X: np.ndarray,
        y: np.ndarray,
        lambdas: Sequence[float] | np.ndarray,
        subsample_indices: np.ndarray | None = None,
    ) -> np.ndarray:
        n_samples, n_features = X.shape
        n_lambdas = len(lambdas)
        if subsample_indices is not None and len(subsample_indices) == n_samples:
            sub_labels = labels_arr[subsample_indices]
        else:
            sub_labels = labels_arr[:n_samples] if len(labels_arr) >= n_samples else labels_arr

        cohort_masks = [sub_labels == c for c in unique_cohorts]
        active_cohort_indices = [t for t, m in enumerate(cohort_masks) if np.sum(m) > 1]
        T = len(active_cohort_indices)

        if T <= 1:
            return create_l1_logistic_fitter(tol=tol, max_iter=max_iter)(
                X, y, lambdas, subsample_indices=subsample_indices
            )

        X_t_list = [X[cohort_masks[t]] for t in active_cohort_indices]
        y_t_list = [y[cohort_masks[t]].astype(np.float64) for t in active_cohort_indices]
        n_t_list = [float(len(yt)) for yt in y_t_list]

        # Standardize X_t column norms: sum_i X_{t, i, k}^2 / n_t = 1
        X_norms = [np.sqrt(np.mean(Xt ** 2, axis=0) + 1e-8) for Xt in X_t_list]
        X_scaled_list = [Xt / norm[np.newaxis, :] for Xt, norm in zip(X_t_list, X_norms)]

        # Lipschitz constant for logistic loss: L_t <= 0.25 * ||X_t||_2^2 / n_t
        L_list = [
            float(0.25 * np.linalg.norm(Xt.T @ Xt, 2) / nt)
            for Xt, nt in zip(X_scaled_list, n_t_list)
        ]
        L_max = max(L_list) if L_list else 1.0
        eta = 1.0 / (L_max + 1e-4)

        # Coefficients: B has shape (n_features, T)
        B = np.zeros((n_features, T), dtype=np.float64)
        active_mask = np.zeros((n_lambdas, n_features), dtype=bool)
        path_lambdas = np.asarray(lambdas, dtype=np.float64)

        for l_idx, lam in enumerate(path_lambdas):
            Y = B.copy()
            for m in range(1, max_iter + 1):
                grad_Y = np.zeros((n_features, T), dtype=np.float64)
                for t in range(T):
                    p_t = _sigmoid(X_scaled_list[t] @ Y[:, t])
                    grad_Y[:, t] = (X_scaled_list[t].T @ (p_t - y_t_list[t])) / n_t_list[t]

                tilde_B = Y - eta * grad_Y
                row_norms = np.linalg.norm(tilde_B, axis=1, keepdims=True)
                shrink = np.maximum(0.0, 1.0 - (eta * lam) / np.maximum(row_norms, 1e-12))
                B_new = tilde_B * shrink

                Y = B_new + ((m - 1) / (m + 2)) * (B_new - B)
                diff = float(np.max(np.abs(B_new - B)))
                B = B_new
                if diff < tol:
                    break

            active_mask[l_idx] = np.linalg.norm(B, axis=1) > 1e-4

        return active_mask

    return fitter


def create_meta_analysis_cohort_fitter(
    cohort_labels: np.ndarray | Sequence[object],
    base_fitter_type: Literal["lasso", "logistic", "elastic_net"] = "logistic",
    min_cohorts: int = 2,
    l1_ratio: float = 0.7,
) -> PathFitter:
    """Create a Meta-Analytic Replicated Consensus Path Fitter.

    Fits independent high-dimensional regularization paths within each cohort
    and selects features that cross replication thresholds across multiple independent cohorts.

    A feature k is active at penalty lambda_j iff it is selected in at least
    min(min_cohorts, T_active) independent cohorts.

    Parameters
    ----------
    cohort_labels : np.ndarray or Sequence
        Array of cohort identifiers for each observation.
    base_fitter_type : {"lasso", "logistic", "elastic_net"}, default "logistic"
        Base shrinkage model applied to each cohort independently.
    min_cohorts : int, default 2
        Minimum number of independent cohorts in which a feature must be selected.
    l1_ratio : float, default 0.7
        Elastic Net ratio if base_fitter_type == "elastic_net".

    Returns
    -------
    PathFitter
    """
    labels_arr = np.asarray(cohort_labels)
    unique_cohorts = tuple(np.unique(labels_arr))

    def fitter(
        X: np.ndarray,
        y: np.ndarray,
        lambdas: Sequence[float] | np.ndarray,
        subsample_indices: np.ndarray | None = None,
    ) -> np.ndarray:
        n_samples, n_features = X.shape
        if subsample_indices is not None and len(subsample_indices) == n_samples:
            sub_labels = labels_arr[subsample_indices]
        else:
            sub_labels = labels_arr[:n_samples] if len(labels_arr) >= n_samples else labels_arr

        cohort_masks = [sub_labels == c for c in unique_cohorts]
        active_cohort_indices = [t for t, m in enumerate(cohort_masks) if np.sum(m) >= 4]
        T = len(active_cohort_indices)

        if T <= 1:
            if base_fitter_type == "logistic":
                return create_l1_logistic_fitter()(X, y, lambdas, subsample_indices=subsample_indices)
            elif base_fitter_type == "elastic_net":
                return create_elastic_net_fitter(l1_ratio=l1_ratio)(
                    X, y, lambdas, subsample_indices=subsample_indices
                )
            return default_lasso_path_fitter(X, y, lambdas, subsample_indices=subsample_indices)

        cohort_masks_active = [cohort_masks[t] for t in active_cohort_indices]
        cohort_masks_res = []

        for c_mask in cohort_masks_active:
            X_t = X[c_mask]
            y_t = y[c_mask]
            if base_fitter_type == "logistic" and len(np.unique(y_t)) >= 2:
                sub_fitter = create_l1_logistic_fitter()
            elif base_fitter_type == "elastic_net":
                sub_fitter = create_elastic_net_fitter(l1_ratio=l1_ratio)
            else:
                sub_fitter = default_lasso_path_fitter

            mask_t = sub_fitter(X_t, y_t, lambdas)
            cohort_masks_res.append(mask_t.astype(int))

        consensus_counts = np.sum(np.stack(cohort_masks_res, axis=0), axis=0)
        req_threshold = min(min_cohorts, T)
        return (consensus_counts >= req_threshold).astype(bool)

    return fitter


def create_multistudy_invariant_cohort_fitter(
    cohort_labels: np.ndarray | Sequence[object],
    gamma: float = 1.0,
    tol: float = 1e-3,
    max_iter: int = 50,
) -> PathFitter:
    """Create a Multi-Study Invariant Prediction Fitter.

    Objective:
    min sum_t (1 / (2 n_t)) ||y_t - X_t beta_t||^2 + lambda ||beta_bar||_1 + gamma sum_t ||beta_t - beta_bar||_2^2
    where beta_bar is the consensus invariant effect across environments and gamma penalizes cross-cohort heterogeneity.

    Parameters
    ----------
    cohort_labels : np.ndarray or Sequence
        Array of cohort identifiers for each observation.
    gamma : float, default 1.0
        Strength of penalty on cross-cohort parameter divergence.
    tol : float, default 1e-3
        Convergence tolerance.
    max_iter : int, default 50
        Maximum iterations per penalty value.

    Returns
    -------
    PathFitter
    """
    labels_arr = np.asarray(cohort_labels)
    unique_cohorts = tuple(np.unique(labels_arr))

    def fitter(
        X: np.ndarray,
        y: np.ndarray,
        lambdas: Sequence[float] | np.ndarray,
        subsample_indices: np.ndarray | None = None,
    ) -> np.ndarray:
        n_samples, n_features = X.shape
        n_lambdas = len(lambdas)
        if subsample_indices is not None and len(subsample_indices) == n_samples:
            sub_labels = labels_arr[subsample_indices]
        else:
            sub_labels = labels_arr[:n_samples] if len(labels_arr) >= n_samples else labels_arr

        cohort_masks = [sub_labels == c for c in unique_cohorts]
        active_cohort_indices = [t for t, m in enumerate(cohort_masks) if np.sum(m) > 1]
        T = len(active_cohort_indices)

        if T <= 1:
            return default_lasso_path_fitter(X, y, lambdas, subsample_indices=subsample_indices)

        X_t_list = [X[cohort_masks[t]] for t in active_cohort_indices]
        y_t_list = [y[cohort_masks[t]].astype(np.float64) for t in active_cohort_indices]
        n_t_list = [float(len(yt)) for yt in y_t_list]

        X_norms = [np.sqrt(np.mean(Xt ** 2, axis=0) + 1e-8) for Xt in X_t_list]
        X_scaled_list = [Xt / norm[np.newaxis, :] for Xt, norm in zip(X_t_list, X_norms)]

        L_list = [
            float(np.linalg.norm(Xt.T @ Xt, 2) / nt + 2.0 * gamma)
            for Xt, nt in zip(X_scaled_list, n_t_list)
        ]
        L_max = max(L_list) if L_list else 1.0
        eta = 1.0 / (L_max + 1e-4)

        B = np.zeros((n_features, T), dtype=np.float64)
        beta_bar = np.zeros(n_features, dtype=np.float64)
        active_mask = np.zeros((n_lambdas, n_features), dtype=bool)
        path_lambdas = np.asarray(lambdas, dtype=np.float64)

        for l_idx, lam in enumerate(path_lambdas):
            for _ in range(max_iter):
                for t in range(T):
                    r_t = y_t_list[t] - X_scaled_list[t] @ B[:, t]
                    grad_t = - (X_scaled_list[t].T @ r_t) / n_t_list[t] + 2.0 * gamma * (B[:, t] - beta_bar)
                    B[:, t] -= eta * grad_t

                grad_bar = 2.0 * gamma * (T * beta_bar - np.sum(B, axis=1))
                tilde_bar = beta_bar - eta * grad_bar
                beta_bar_new = np.sign(tilde_bar) * np.maximum(0.0, np.abs(tilde_bar) - eta * lam)
                diff = float(np.max(np.abs(beta_bar_new - beta_bar)))
                beta_bar = beta_bar_new
                if diff < tol:
                    break

            active_mask[l_idx] = np.abs(beta_bar) > 1e-4

        return active_mask

    return fitter


def create_glmm_lasso_cohort_fitter(
    cohort_labels: np.ndarray | Sequence[object],
    tol: float = 1e-3,
    max_iter: int = 200,
) -> PathFitter:
    """Create a Penalized GLMM Logistic Regression Fitter with cohort random intercepts.

    Absorbs cohort-level baseline response differences on the log-odds scale:
        logit(P(y_i = 1)) = X_i beta + a_{c_i}
    where a_c is the cohort log-odds offset and beta is L1-penalized.

    Parameters
    ----------
    cohort_labels : np.ndarray or Sequence
        Array of cohort identifiers for each observation.
    tol : float, default 1e-3
        Convergence tolerance.
    max_iter : int, default 200
        Maximum coordinate descent / FISTA iterations per penalty.

    Returns
    -------
    PathFitter
    """
    labels_arr = np.asarray(cohort_labels)
    unique_cohorts = tuple(np.unique(labels_arr))

    def fitter(
        X: np.ndarray,
        y: np.ndarray,
        lambdas: Sequence[float] | np.ndarray,
        subsample_indices: np.ndarray | None = None,
    ) -> np.ndarray:
        n_samples, n_features = X.shape
        n_lambdas = len(lambdas)
        if subsample_indices is not None and len(subsample_indices) == n_samples:
            sub_labels = labels_arr[subsample_indices]
        else:
            sub_labels = labels_arr[:n_samples] if len(labels_arr) >= n_samples else labels_arr

        y_bin = (y > np.median(y)).astype(np.float64) if len(np.unique(y)) > 2 else y.astype(np.float64)

        mean_p = (np.sum(y_bin) + 0.5) / (n_samples + 1.0)
        logit_mean = np.log(mean_p / (1.0 - mean_p))

        offsets = np.zeros(n_samples, dtype=np.float64)
        for c in unique_cohorts:
            c_mask = sub_labels == c
            if np.sum(c_mask) > 1:
                p_c = (np.sum(y_bin[c_mask]) + 0.5) / (np.sum(c_mask) + 1.0)
                logit_c = np.log(p_c / (1.0 - p_c))
                offsets[c_mask] = logit_c - logit_mean

        X_norms = np.sqrt(np.mean(X ** 2, axis=0) + 1e-8)
        X_scaled = X / X_norms[np.newaxis, :]

        L = float(0.25 * np.linalg.norm(X_scaled.T @ X_scaled, 2) / n_samples)
        eta = 1.0 / (L + 1e-4)

        beta = np.zeros(n_features, dtype=np.float64)
        active_mask = np.zeros((n_lambdas, n_features), dtype=bool)
        path_lambdas = np.asarray(lambdas, dtype=np.float64)

        for l_idx, lam in enumerate(path_lambdas):
            Y = beta.copy()
            for m in range(1, max_iter + 1):
                p_pred = _sigmoid(X_scaled @ Y + offsets)
                grad = (X_scaled.T @ (p_pred - y_bin)) / n_samples

                tilde_beta = Y - eta * grad
                beta_new = np.sign(tilde_beta) * np.maximum(0.0, np.abs(tilde_beta) - eta * lam)

                Y = beta_new + ((m - 1) / (m + 2)) * (beta_new - beta)
                diff = float(np.max(np.abs(beta_new - beta)))
                beta = beta_new
                if diff < tol:
                    break

            active_mask[l_idx] = np.abs(beta) > 1e-4

        return active_mask

    return fitter


def _prox_sorted_l1(
    v: np.ndarray,
    weights: np.ndarray,
    step_size: float,
) -> np.ndarray:
    """Proximal operator for Ordered Weighted L1 (OWL / SLOPE / OSCAR) norm.

    Solves:
        argmin_beta 0.5 * ||beta - v||_2^2 + step_size * sum_j w_j |beta|_{(j)}
    using the Pool Adjacent Violators Algorithm (PAVA) via isotonic regression.

    Parameters
    ----------
    v : np.ndarray
        Vector of coefficients to threshold, shape (p,).
    weights : np.ndarray
        Non-negative, non-increasing sequence of regularization weights, shape (p,).
    step_size : float
        Gradient descent step size (eta).

    Returns
    -------
    np.ndarray
        Thresholded and clustered coefficient vector, shape (p,).
    """
    if np.all(v == 0.0):
        return np.zeros_like(v)

    abs_v = np.abs(v)
    signs = np.sign(v)

    # Sort absolute values in descending order
    sort_order = np.argsort(-abs_v)
    rev_order = np.argsort(sort_order)
    sorted_abs_v = abs_v[sort_order]

    # Shift by step_size * weights
    shifted = sorted_abs_v - step_size * weights

    # Isotonic regression enforcing non-increasing order: z_1 >= z_2 >= ... >= z_p >= 0
    iso_fit = isotonic_regression(shifted, increasing=False)
    z = np.maximum(iso_fit, 0.0)

    # Restore original order and signs
    return signs * z[rev_order]


def _solve_sorted_l1_path(
    X: np.ndarray,
    y: np.ndarray,
    lambdas: Sequence[float] | np.ndarray,
    base_weights: np.ndarray,
    tol: float = 1e-4,
    max_iter: int = 100,
) -> np.ndarray:
    """Solve an Ordered Weighted L1 path using FISTA with warm restarts.

    Parameters
    ----------
    X : np.ndarray
        Feature matrix of shape (n_samples, n_features).
    y : np.ndarray
        Response vector of shape (n_samples,).
    lambdas : Sequence[float] | np.ndarray
        Sequence of regularization penalties.
    base_weights : np.ndarray
        Normalized penalty pattern of shape (n_features,). For a given penalty lambda,
        the weight vector applied is lambda * base_weights.
    tol : float, default 1e-4
        FISTA convergence tolerance.
    max_iter : int, default 100
        Maximum FISTA iterations per penalty value.

    Returns
    -------
    np.ndarray[bool]
        Active feature selection mask of shape (n_lambdas, n_features).
    """
    n_samples, n_features = X.shape
    n_lambdas = len(lambdas)

    if n_features == 0 or n_samples == 0:
        return np.zeros((n_lambdas, n_features), dtype=bool)

    path_lambdas = np.asarray(lambdas, dtype=np.float64)
    # Sort descending for stable warm starts
    sort_order = np.argsort(-path_lambdas)
    rev_order = np.argsort(sort_order)
    sorted_lambdas = path_lambdas[sort_order]

    # Lipschitz constant L = ||X||_2^2 / n_samples
    L = float(np.linalg.norm(X, 2) ** 2 / max(n_samples, 1))
    eta = 1.0 / (L + 1e-4)

    y_float = y.astype(np.float64)
    beta = np.zeros(n_features, dtype=np.float64)
    active_mask = np.zeros((n_lambdas, n_features), dtype=bool)

    for l_idx, lam in enumerate(sorted_lambdas):
        w = lam * base_weights
        Y = beta.copy()
        t0 = 1.0
        for m in range(1, max_iter + 1):
            grad = - (X.T @ (y_float - X @ Y)) / n_samples
            beta_new = _prox_sorted_l1(Y - eta * grad, w, eta)
            diff = float(np.max(np.abs(beta_new - beta)))

            t1 = (1.0 + np.sqrt(1.0 + 4.0 * t0**2)) / 2.0
            Y = beta_new + ((t0 - 1.0) / t1) * (beta_new - beta)
            t0 = t1
            beta = beta_new
            if diff < tol:
                break

        active_mask[l_idx] = np.abs(beta) > 1e-5

    return active_mask[rev_order]


def create_oscar_fitter(
    kappa: float = 0.5,
    tol: float = 1e-4,
    max_iter: int = 100,
) -> PathFitter:
    """Create an OSCAR (Octagonal Shrinkage and Clustering Algorithm) path fitter.

    OSCAR applies an Ordered Weighted L1 (OWL) penalty with linearly decreasing weights:
        Omega_OSCAR(beta) = lambda * sum_j [ (1 - kappa) + kappa * (p - j) / (p - 1) ] |beta|_{(j)}
    which induces octagonal sparsity polytopes that group collinear features into
    exact clusters with identical non-zero magnitudes.

    Parameters
    ----------
    kappa : float, default 0.5
        Clustering parameter in [0, 1]. kappa = 0 is pure Lasso; higher values
        strengthen pairwise grouping.
    tol : float, default 1e-4
        FISTA convergence tolerance.
    max_iter : int, default 100
        Maximum FISTA iterations per penalty.

    Returns
    -------
    PathFitter
    """
    clipped_kappa = float(np.clip(kappa, 0.0, 1.0))

    def fitter(
        X: np.ndarray,
        y: np.ndarray,
        lambdas: Sequence[float] | np.ndarray,
        subsample_indices: np.ndarray | None = None,
    ) -> np.ndarray:
        p = X.shape[1]
        if p <= 1:
            base_weights = np.ones(p, dtype=np.float64)
        else:
            base_weights = (1.0 - clipped_kappa) + clipped_kappa * np.linspace(1.0, 0.0, p)

        return _solve_sorted_l1_path(
            X=X,
            y=y,
            lambdas=lambdas,
            base_weights=base_weights,
            tol=tol,
            max_iter=max_iter,
        )

    return fitter


def create_slope_fitter(
    q_fdr: float = 0.1,
    tol: float = 1e-4,
    max_iter: int = 100,
) -> PathFitter:
    """Create a SLOPE (Sorted L-One Penalized Estimation) path fitter.

    SLOPE applies Benjamini-Hochberg inspired normal quantile weights:
        w_j = lambda * Phi^{-1}(1 - j * q_fdr / (2p))
    which provides adaptive shrinkage and controls the false discovery rate (FDR)
    of selected features under orthogonal designs.

    Parameters
    ----------
    q_fdr : float, default 0.1
        Target false discovery rate parameter in (0, 1).
    tol : float, default 1e-4
        FISTA convergence tolerance.
    max_iter : int, default 100
        Maximum FISTA iterations per penalty.

    Returns
    -------
    PathFitter
    """
    clipped_q = float(np.clip(q_fdr, 1e-4, 0.99))

    def fitter(
        X: np.ndarray,
        y: np.ndarray,
        lambdas: Sequence[float] | np.ndarray,
        subsample_indices: np.ndarray | None = None,
    ) -> np.ndarray:
        p = X.shape[1]
        ranks = np.arange(1, p + 1, dtype=np.float64)
        quantiles = 1.0 - (ranks * clipped_q) / (2.0 * max(p, 1))
        base_weights = np.maximum(norm.ppf(quantiles), 0.0)

        return _solve_sorted_l1_path(
            X=X,
            y=y,
            lambdas=lambdas,
            base_weights=base_weights,
            tol=tol,
            max_iter=max_iter,
        )

    return fitter


def create_oscar_cohort_fitter(
    cohort_labels: np.ndarray | Sequence[object],
    kappa: float = 0.5,
    tol: float = 1e-4,
    max_iter: int = 100,
) -> PathFitter:
    """Create a cohort-adjusted OSCAR path fitter using Frisch-Waugh-Lovell projection.

    Projects out cohort fixed effects and batch shifts prior to applying
    OSCAR octagonal shrinkage and collinear feature clustering.

    Parameters
    ----------
    cohort_labels : np.ndarray or Sequence
        Cohort identifiers for each observation.
    kappa : float, default 0.5
        Clustering parameter in [0, 1].
    tol : float, default 1e-4
        FISTA convergence tolerance.
    max_iter : int, default 100
        Maximum FISTA iterations per penalty.

    Returns
    -------
    PathFitter
    """
    labels_arr = np.asarray(cohort_labels)
    unique_cohorts = tuple(np.unique(labels_arr))
    base_fitter = create_oscar_fitter(kappa=kappa, tol=tol, max_iter=max_iter)

    def fitter(
        X: np.ndarray,
        y: np.ndarray,
        lambdas: Sequence[float] | np.ndarray,
        subsample_indices: np.ndarray | None = None,
    ) -> np.ndarray:
        n_samples = X.shape[0]
        if subsample_indices is not None and len(subsample_indices) == n_samples:
            sub_labels = labels_arr[subsample_indices]
        else:
            sub_labels = labels_arr[:n_samples] if len(labels_arr) >= n_samples else labels_arr

        cohort_masks = [sub_labels == c for c in unique_cohorts]
        active_cohort_indices = [t for t, m in enumerate(cohort_masks) if np.sum(m) > 1]
        if len(active_cohort_indices) <= 1:
            return base_fitter(X, y, lambdas, subsample_indices=subsample_indices)

        X_res = X.copy()
        y_res = y.astype(np.float64).copy()
        for c in unique_cohorts:
            c_mask = sub_labels == c
            if np.sum(c_mask) > 1:
                X_res[c_mask] = X_res[c_mask] - np.mean(X_res[c_mask], axis=0)
                y_res[c_mask] = y_res[c_mask] - np.mean(y_res[c_mask])

        return base_fitter(X_res, y_res, lambdas, subsample_indices=subsample_indices)

    return fitter


def create_slope_cohort_fitter(
    cohort_labels: np.ndarray | Sequence[object],
    q_fdr: float = 0.1,
    tol: float = 1e-4,
    max_iter: int = 100,
) -> PathFitter:
    """Create a cohort-adjusted SLOPE path fitter using Frisch-Waugh-Lovell projection.

    Projects out cohort fixed effects and batch shifts prior to applying
    SLOPE FDR-controlling sorted L1 penalization.

    Parameters
    ----------
    cohort_labels : np.ndarray or Sequence
        Cohort identifiers for each observation.
    q_fdr : float, default 0.1
        Target false discovery rate parameter in (0, 1).
    tol : float, default 1e-4
        FISTA convergence tolerance.
    max_iter : int, default 100
        Maximum FISTA iterations per penalty.

    Returns
    -------
    PathFitter
    """
    labels_arr = np.asarray(cohort_labels)
    unique_cohorts = tuple(np.unique(labels_arr))
    base_fitter = create_slope_fitter(q_fdr=q_fdr, tol=tol, max_iter=max_iter)

    def fitter(
        X: np.ndarray,
        y: np.ndarray,
        lambdas: Sequence[float] | np.ndarray,
        subsample_indices: np.ndarray | None = None,
    ) -> np.ndarray:
        n_samples = X.shape[0]
        if subsample_indices is not None and len(subsample_indices) == n_samples:
            sub_labels = labels_arr[subsample_indices]
        else:
            sub_labels = labels_arr[:n_samples] if len(labels_arr) >= n_samples else labels_arr

        cohort_masks = [sub_labels == c for c in unique_cohorts]
        active_cohort_indices = [t for t, m in enumerate(cohort_masks) if np.sum(m) > 1]
        if len(active_cohort_indices) <= 1:
            return base_fitter(X, y, lambdas, subsample_indices=subsample_indices)

        X_res = X.copy()
        y_res = y.astype(np.float64).copy()
        for c in unique_cohorts:
            c_mask = sub_labels == c
            if np.sum(c_mask) > 1:
                X_res[c_mask] = X_res[c_mask] - np.mean(X_res[c_mask], axis=0)
                y_res[c_mask] = y_res[c_mask] - np.mean(y_res[c_mask])

        return base_fitter(X_res, y_res, lambdas, subsample_indices=subsample_indices)

    return fitter


def build_fitter(
    fitter_name: str,
    cohort_labels: np.ndarray | None = None,
    l1_ratio: float = 0.7,
    kappa: float = 0.5,
    q_fdr: float = 0.1,
) -> PathFitter:
    """Instantiate a PathFitter by name with automatic multi-cohort fallback."""
    has_multi_cohort = cohort_labels is not None and len(np.unique(cohort_labels)) > 1

    match fitter_name:
        case "elastic_net":
            return create_elastic_net_fitter(l1_ratio=l1_ratio)
        case "logistic":
            return create_l1_logistic_fitter()
        case "rf":
            return create_tree_importance_fitter(estimator_type="rf")
        case "cohort_adjusted":
            if has_multi_cohort and cohort_labels is not None:
                return create_cohort_adjusted_fitter(cohort_labels, base_fitter_type="lasso")
            return create_lasso_fitter()
        case "group_lasso":
            if has_multi_cohort and cohort_labels is not None:
                return create_group_lasso_cohort_fitter(cohort_labels)
            return create_lasso_fitter()
        case "merf":
            if has_multi_cohort and cohort_labels is not None:
                return create_merf_cohort_fitter(cohort_labels)
            return create_tree_importance_fitter(estimator_type="rf")
        case "multitask_logistic":
            if has_multi_cohort and cohort_labels is not None:
                return create_multitask_logistic_cohort_fitter(cohort_labels)
            return create_l1_logistic_fitter()
        case "meta_analysis":
            if has_multi_cohort and cohort_labels is not None:
                return create_meta_analysis_cohort_fitter(cohort_labels, base_fitter_type="logistic")
            return create_l1_logistic_fitter()
        case "multistudy_invariant":
            if has_multi_cohort and cohort_labels is not None:
                return create_multistudy_invariant_cohort_fitter(cohort_labels)
            return create_lasso_fitter()
        case "glmm_lasso":
            if has_multi_cohort and cohort_labels is not None:
                return create_glmm_lasso_cohort_fitter(cohort_labels)
            return create_l1_logistic_fitter()
        case "oscar":
            if has_multi_cohort and cohort_labels is not None:
                return create_oscar_cohort_fitter(cohort_labels, kappa=kappa)
            return create_oscar_fitter(kappa=kappa)
        case "slope":
            if has_multi_cohort and cohort_labels is not None:
                return create_slope_cohort_fitter(cohort_labels, q_fdr=q_fdr)
            return create_slope_fitter(q_fdr=q_fdr)
        case "lasso" | _:
            return create_lasso_fitter()
