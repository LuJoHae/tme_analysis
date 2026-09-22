"""Parameter estimation for Negative Binomial (Gamma-Poisson) count models.

Provides three estimation strategies:
1. Method of Moments (MoM): Closed-form analytical estimator with depth scaling.
2. Maximum Likelihood Estimation (MLE): Newton-Raphson profile likelihood on digamma/trigamma.
3. Empirical Bayes / Trended MoM: Parametric dispersion-mean trend fitting with shrinkage.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import scipy.sparse as sp
from scipy.special import polygamma, psi
from returns.maybe import Maybe, Some
from returns.result import Failure, Result, Success

from ..logging import get_logger
from ..types import NBEstimationMethod

logger = get_logger("transforms.nb_inference")


def compute_size_factors(X: np.ndarray) -> np.ndarray:
    """Compute per-cell library size scaling factors normalized to median depth.

    s_i = library_size_i / median_library_size
    """
    depths = np.sum(X, axis=1)
    # Avoid divide-by-zero for empty cells
    depths = np.maximum(depths, 1.0)
    median_depth = np.median(depths)
    if median_depth <= 0:
        median_depth = 1.0
    return depths / median_depth


def fit_nb_moments(
    X: np.ndarray,
    size_factors: np.ndarray,
    min_dispersion: float = 1e-4,
    max_dispersion: float = 10.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit gene-level mean and dispersion using Method of Moments (MoM) in closed form.

    mu_g = sum_i(X_{i,g}) / sum_i(s_i)
    sigma2_g = 1 / (N - 1) sum_i (X_{i,g} / s_i - mu_g)^2
    alpha_g = (sigma2_g - mu_g) / (mu_g^2)
    """
    N, G = X.shape
    if N < 2:
        means = np.mean(X, axis=0)
        dispersions = np.full(G, min_dispersion, dtype=np.float64)
        return means, dispersions

    s = size_factors.reshape(-1, 1)
    sum_s = np.sum(size_factors)
    means = np.sum(X, axis=0) / max(sum_s, 1e-8)

    # Normalized expression
    norm_X = X / np.maximum(s, 1e-8)
    variances = np.var(norm_X, axis=0, ddof=1)

    # Moments estimator: alpha = (Var - mu) / (mu^2)
    mu_sq = np.maximum(means**2, 1e-12)
    raw_alpha = (variances - means) / mu_sq

    dispersions = np.clip(raw_alpha, min_dispersion, max_dispersion)
    return means.astype(np.float64), dispersions.astype(np.float64)


def fit_nb_mle(
    X: np.ndarray,
    size_factors: np.ndarray,
    max_iter: int = 25,
    tol: float = 1e-4,
    min_dispersion: float = 1e-4,
    max_dispersion: float = 10.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit gene-level mean and dispersion via Maximum Likelihood Estimation (MLE).

    Uses Method of Moments to initialize dispersion, then performs vectorized
    Newton-Raphson updates on the profile log-likelihood:
    f(alpha) = sum_i [ psi(y_i + 1/alpha) - psi(1/alpha) - log(1 + alpha * mu_i) - (y_i - mu_i)/(1/alpha + mu_i) ] = 0
    """
    N, G = X.shape
    means, mom_dispersions = fit_nb_moments(
        X, size_factors, min_dispersion=min_dispersion, max_dispersion=max_dispersion
    )

    s = size_factors.reshape(-1, 1)
    mu_mat = s * means.reshape(1, -1)  # shape (N, G)

    # Parameterize in log(alpha) for numerical stability and positivity guarantee
    log_alpha = np.log(np.clip(mom_dispersions, min_dispersion, max_dispersion))

    # Mask genes with zero mean expression
    active_genes = means > 1e-8

    for _ in range(max_iter):
        alpha = np.exp(log_alpha)
        inv_alpha = 1.0 / np.maximum(alpha, 1e-8)

        # Profile log-likelihood gradient w.r.t alpha:
        # dL/d(alpha) = -1/alpha^2 * sum_i [ psi(y_i + 1/alpha) - psi(1/alpha) ] + sum_i [ y_i/(alpha*(1 + alpha*mu)) - mu/(1 + alpha*mu) ]
        # More directly: dL/d(inv_alpha) = sum_i [ psi(y_i + inv_alpha) - psi(inv_alpha) + log(inv_alpha) - log(inv_alpha + mu_mat) ]
        # Let r = 1/alpha:
        r = inv_alpha.reshape(1, -1)
        # Score w.r.t log(alpha): dL/d(log(alpha)) = alpha * dL/d(alpha)
        # Using d/dr:
        psi_diff = psi(X + r) - psi(r)
        log_term = np.log(r) - np.log(np.maximum(r + mu_mat, 1e-8))
        score_r = np.sum(psi_diff + log_term, axis=0)  # dL/dr

        # Second derivative d2L/dr2:
        trigamma_diff = polygamma(1, X + r) - polygamma(1, r)
        rational_term = (1.0 / r) - (1.0 / np.maximum(r + mu_mat, 1e-8))
        hessian_r = np.sum(trigamma_diff + rational_term, axis=0)

        # Gradient w.r.t log(alpha): since r = exp(-log(alpha)), dr/d(log_alpha) = -r
        score_log_alpha = -r[0] * score_r
        hessian_log_alpha = (r[0] ** 2) * hessian_r + r[0] * score_r

        # Newton step: delta = -score / hessian
        # Ensure negative curvature (hessian < 0) for concave maximization
        safe_hessian = np.minimum(hessian_log_alpha, -1e-6)
        step = -score_log_alpha / safe_hessian

        step = np.clip(step, -2.0, 2.0)  # Step bounding
        step = np.where(active_genes, step, 0.0)

        log_alpha += step

        if np.max(np.abs(step)) < tol:
            break

    final_alpha = np.clip(np.exp(log_alpha), min_dispersion, max_dispersion)
    return means.astype(np.float64), final_alpha.astype(np.float64)


def fit_nb_empirical_bayes(
    X: np.ndarray,
    size_factors: np.ndarray,
    min_dispersion: float = 1e-4,
    max_dispersion: float = 10.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit gene-level dispersions with Empirical Bayes trend-shrinkage.

    1. Computes initial gene moments (means and raw dispersions).
    2. Fits parametric curve: alpha_trend(mu) = a0 + a1 / mu via robust regression.
    3. Shrinks gene dispersions toward the trend:
       alpha_shrunk = w * alpha_raw + (1 - w) * alpha_trend
       where w = N / (N + tau) is the Bayesian shrinkage factor.
    """
    N, G = X.shape
    means, raw_dispersions = fit_nb_moments(
        X, size_factors, min_dispersion=min_dispersion, max_dispersion=max_dispersion
    )

    # Filter informative genes for trend fitting
    valid_mask = (means > 0.01) & (raw_dispersions > min_dispersion)
    if np.sum(valid_mask) < 5:
        return means, raw_dispersions

    valid_means = means[valid_mask]
    valid_dispersions = raw_dispersions[valid_mask]

    # Model: alpha ~ a0 + a1 / mu  =>  alpha ~ [1, 1/mu] @ [a0, a1]
    A = np.column_stack([np.ones_like(valid_means), 1.0 / valid_means])
    try:
        coef, _, _, _ = np.linalg.lstsq(A, valid_dispersions, rcond=None)
        a0, a1 = coef[0], coef[1]
        a0 = max(min_dispersion, a0)
        a1 = max(0.0, a1)
    except Exception:
        a0, a1 = float(np.median(valid_dispersions)), 0.0

    # Trend value for all genes
    trend_dispersions = np.clip(
        a0 + a1 / np.maximum(means, 1e-6), min_dispersion, max_dispersion
    )

    # Shrinkage weight w = N / (N + prior_weight), default prior weight 20
    prior_weight = 20.0
    weight = float(N / (N + prior_weight))
    shrunk_dispersions = weight * raw_dispersions + (1.0 - weight) * trend_dispersions

    final_dispersions = np.clip(shrunk_dispersions, min_dispersion, max_dispersion)
    return means.astype(np.float64), final_dispersions.astype(np.float64)


def infer_dataset_nb_parameters(
    adata: ad.AnnData,
    method: NBEstimationMethod = NBEstimationMethod.MOMENTS,
    cluster_key: Maybe[str] = None,
    min_dispersion: float = 1e-4,
    max_dispersion: float = 10.0,
) -> Result[dict[str, tuple[np.ndarray, np.ndarray]], str]:
    """Pure functional wrapper to infer NB parameters across AnnData, optionally by cluster.

    Returns mapping: cluster_name -> (means, dispersions)
    """
    try:
        logger.info(
            "Inferring Negative Binomial parameters using %s (%d cells x %d genes)...",
            method.value,
            adata.n_obs,
            adata.n_vars,
        )
        X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X, dtype=np.float64)
        if method == NBEstimationMethod.SANITY:
            from ..preprocessing.sanity import run_sanity_normalization

            sanity_res = run_sanity_normalization(adata)
            if not isinstance(sanity_res, Success):
                return Failure(sanity_res.failure())
            s_adata = sanity_res.unwrap()
            # In Sanity, means correspond to baseline alphas * median depth, dispersions to variance v_g
            N_c = np.sum(X, axis=1)
            med_depth = float(np.median(N_c))
            means = s_adata.var["sanity_baseline_alpha"].values * med_depth
            dispersions = np.clip(
                s_adata.var["sanity_variance"].values, min_dispersion, max_dispersion
            )
            logger.info("Successfully inferred Sanity-based NB parameters (median dispersion: %.3f)", float(np.median(dispersions)))
            return Success({"global": (means, dispersions)})

        fit_fn = {
            NBEstimationMethod.MOMENTS: fit_nb_moments,
            NBEstimationMethod.MLE: fit_nb_mle,
            NBEstimationMethod.EMPIRICAL_BAYES: fit_nb_empirical_bayes,
        }[method]

        key = cluster_key.value_or(None) if isinstance(cluster_key, Some) else None

        if key is not None and key in adata.obs:
            clusters = adata.obs[key].astype(str).values
            unique_clusters = np.unique(clusters)
            logger.info("Inferring cluster-specific NB parameters for %d clusters (%s)...", len(unique_clusters), key)
            results: dict[str, tuple[np.ndarray, np.ndarray]] = {}

            for c in unique_clusters:
                idx = clusters == c
                sub_X = X[idx]
                sf = compute_size_factors(sub_X)
                mu, alpha = fit_fn(
                    sub_X,
                    sf,
                    min_dispersion=min_dispersion,
                    max_dispersion=max_dispersion,
                )
                results[str(c)] = (mu, alpha)

            logger.info("Successfully inferred NB parameters for %d clusters", len(results))
            return Success(results)

        # Global fit across all cells
        sf = compute_size_factors(X)
        mu, alpha = fit_fn(
            X,
            sf,
            min_dispersion=min_dispersion,
            max_dispersion=max_dispersion,
        )
        logger.info(
            "Successfully inferred global NB parameters (median mu: %.2f, median alpha: %.3f)",
            float(np.median(mu)),
            float(np.median(alpha)),
        )
        return Success({"global": (mu, alpha)})
    except Exception as exc:
        msg = f"Failed to infer Negative Binomial parameters: {exc}"
        logger.error(msg)
        return Failure(msg)
