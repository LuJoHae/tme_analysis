"""Sanity: Bayesian normalization and denoising from first principles.

Reference:
    Breda, J., Zavolan, M., & van Nimwegen, E. (2021).
    Bayesian inference of gene expression states from single-cell RNA-seq data.
    Nature Biotechnology, 39(8), 1008-1016.
    https://doi.org/10.1038/s41587-021-00875-x
"""

from __future__ import annotations

import time
import anndata as ad
import numpy as np
import scipy.sparse as sp
from scipy.optimize import root_scalar
from scipy.special import lambertw
from returns.result import Failure, Result, Success

from ..logging import get_logger
from ..models import SanityConfig

logger = get_logger("preprocessing.sanity")


def _solve_qg(y_gc: np.ndarray, v_g: float, n_g: float) -> float:
    """Solve the monotonic 1D consistency equation for q_g:

    sum_c f_{gc}(q_g) = 1
    where f_{gc} = W(exp(-q_g + y_{gc}) * v_g * (n_g + 1)) / (v_g * (n_g + 1)).
    """
    scale = v_g * (n_g + 1.0)
    # y_gc = log(N_c) + v_g * n_{gc}
    max_y = float(np.max(y_gc))

    def objective(q: float) -> float:
        # z = exp(-q + y_gc) * scale
        # Use log_z for numerical stability: log_z = -q + y_gc + log(scale)
        log_z = -q + y_gc + np.log(scale)
        # For large positive log_z, W(exp(log_z)) ~ log_z - log(log_z)
        # Scipy lambertw takes complex/float z
        # Bound z to avoid exp overflow
        clipped_log_z = np.clip(log_z, -50.0, 50.0)
        z = np.exp(clipped_log_z)
        w_val = np.real(lambertw(z))
        f_vals = w_val / scale
        return float(np.sum(f_vals) - 1.0)

    # Initial bracket based on y_gc
    # If q is very large, z -> 0, f -> 0, objective -> -1
    # If q is very small, z -> inf, objective > 0
    q_left = max_y - 20.0
    q_right = max_y + 20.0

    try:
        # Check signs to ensure valid bracket
        f_left = objective(q_left)
        f_right = objective(q_right)

        expand_count = 0
        while f_left <= 0 and expand_count < 10:
            q_left -= 10.0
            f_left = objective(q_left)
            expand_count += 1

        expand_count = 0
        while f_right >= 0 and expand_count < 10:
            q_right += 10.0
            f_right = objective(q_right)
            expand_count += 1

        sol = root_scalar(objective, bracket=[q_left, q_right], method="brentq")
        return float(sol.root)
    except Exception:
        # Fallback to mean approximation
        return float(np.mean(y_gc))


def run_sanity_single_gene(
    n_gc: np.ndarray,
    N_c: np.ndarray,
    v_grid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Run Sanity Bayesian inference for a single gene across all cells.

    Returns:
        delta_gc: posterior mean log-fold changes for each cell
        error_gc: posterior standard error bars for each cell
        v_g: inferred optimal true gene variance
    """
    C = len(n_gc)
    n_g = float(np.sum(n_gc))

    if n_g == 0:
        # Gene unexpressed across all cells: variance is minimal, delta is 0
        return np.zeros(C, dtype=np.float32), np.zeros(C, dtype=np.float32), float(v_grid[0])

    log_Nc = np.log(np.maximum(N_c, 1.0))
    n_v = len(v_grid)

    log_likelihoods = np.zeros(n_v, dtype=np.float64)
    delta_star_mat = np.zeros((n_v, C), dtype=np.float64)
    error_star_mat = np.zeros((n_v, C), dtype=np.float64)

    for i, v_g in enumerate(v_grid):
        inv_vg = 1.0 / v_g
        y_gc = log_Nc + v_g * n_gc
        q_g = _solve_qg(y_gc, v_g, n_g)

        scale = v_g * (n_g + 1.0)
        log_z = np.clip(-q_g + y_gc + np.log(scale), -50.0, 50.0)
        z = np.exp(log_z)
        w_val = np.real(lambertw(z))
        f_gc = w_val / scale
        f_gc = np.maximum(f_gc, 1e-12)

        delta_star = np.log(f_gc) - log_Nc + q_g
        delta_star_mat[i] = delta_star

        # Diagonal of M^g inverse (error bars):
        # M_{cc} = (n_g + 1) * f_gc + 1 / v_g
        denom = (n_g + 1.0) * f_gc + inv_vg
        diag_inv = 1.0 / np.maximum(denom, 1e-12)
        error_star_mat[i] = np.sqrt(diag_inv)

        # Log determinant of M^g using Matrix Determinant Lemma:
        # det(M) = [1 - sum_c ((n_g + 1)*f_gc^2 / ((n_g+1)*f_gc + 1/v_g))] * prod_c ((n_g+1)*f_gc + 1/v_g)
        factor = 1.0 - np.sum(((n_g + 1.0) * (f_gc**2)) / denom)
        factor = max(1e-8, factor)
        log_det_M = np.log(factor) + np.sum(np.log(denom))

        # Optimal log-likelihood L^*(v_g) from Eq. 20
        sum_exp = np.sum(N_c * np.exp(delta_star))
        l_star = (
            -0.5 * C * np.log(v_g)
            - 0.5 * inv_vg * np.sum(delta_star**2)
            + np.sum(n_gc * delta_star)
            - (n_g + 1.0) * np.log(max(1e-12, sum_exp))
        )

        # Laplace log marginal likelihood: log P(n_g | v_g) = L^*(v_g) - 0.5 * log_det(M)
        log_likelihoods[i] = l_star - 0.5 * log_det_M

    # Posterior over v_g with scale prior P(v_g) ~ 1 / v_g
    # log P(v_b | n_g) = log_likelihoods - log(v_grid)
    log_posterior = log_likelihoods - np.log(v_grid)
    # Normalize with logsumexp
    max_lp = np.max(log_posterior)
    weights = np.exp(log_posterior - max_lp)
    p_v = weights / np.sum(weights)

    # Posterior expectation values
    mean_delta = np.sum(p_v[:, None] * delta_star_mat, axis=0).astype(np.float32)
    mean_error = np.sum(p_v[:, None] * error_star_mat, axis=0).astype(np.float32)
    best_vg = float(np.sum(p_v * v_grid))

    return mean_delta, mean_error, best_vg


def run_sanity_normalization(
    adata: ad.AnnData,
    config: SanityConfig = SanityConfig(),
) -> Result[ad.AnnData, str]:
    """Perform Sanity Bayesian normalization and sampling-noise correction.

    Estimates Log-Transcription Quotients (LTQs) and error bars without tunable parameters.
    Results are added to AnnData:
      - `.layers["sanity_ltq"]`: Inferred mean log-fold changes <delta_{gc}>
      - `.layers["sanity_error"]`: Posterior standard errors (error bars) sigma_{gc}
      - `.var["sanity_variance"]`: Inferred gene true biological variance v_g
      - `.var["sanity_baseline_alpha"]`: Inferred dataset-wide mean quotient alpha_g
    """
    try:
        new_adata = adata.copy()
        X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X, dtype=np.float64)

        N_cells, G_genes = X.shape
        logger.info(
            "Starting Sanity Bayesian normalization for %d cells x %d genes (v_bins=%d, v_min=%.2e, v_max=%.2e)...",
            N_cells,
            G_genes,
            config.n_bins,
            config.v_min,
            config.v_max,
        )
        start_time = time.time()

        # Cell library sizes N_c
        N_c = np.sum(X, axis=1)
        N_c = np.maximum(N_c, 1.0)

        # Log-spaced grid for variance v_g
        v_grid = np.geomspace(config.v_min, config.v_max, num=config.n_bins)

        ltq_mat = np.zeros((N_cells, G_genes), dtype=np.float32)
        error_mat = np.zeros((N_cells, G_genes), dtype=np.float32)
        variances = np.zeros(G_genes, dtype=np.float64)

        # Baseline quotient alpha_g = sum_c(n_gc) / sum_c(N_c)
        total_depth = np.sum(N_c)
        baseline_alphas = (np.sum(X, axis=0) / max(total_depth, 1.0)).astype(np.float64)

        log_interval = max(500, G_genes // 5) if G_genes >= 1000 else max(100, G_genes // 2)

        for g in range(G_genes):
            if g > 0 and g % log_interval == 0:
                pct = (g / G_genes) * 100
                logger.info("Sanity Bayesian progress: %d / %d genes (%.1f%%)...", g, G_genes, pct)

            n_gc = X[:, g]
            delta_g, err_g, v_g = run_sanity_single_gene(n_gc, N_c, v_grid)
            ltq_mat[:, g] = delta_g
            error_mat[:, g] = err_g
            variances[g] = v_g

        # Save to AnnData
        new_adata.layers["sanity_ltq"] = ltq_mat
        new_adata.layers["sanity_error"] = error_mat
        new_adata.var["sanity_variance"] = variances
        new_adata.var["sanity_baseline_alpha"] = baseline_alphas

        elapsed = max(0.01, time.time() - start_time)
        logger.info(
            "Successfully completed Sanity normalization for %d genes across %d cells in %.2fs",
            G_genes,
            N_cells,
            elapsed,
        )
        return Success(new_adata)
    except Exception as exc:
        msg = f"Sanity normalization failed: {exc}"
        logger.error(msg)
        return Failure(msg)
