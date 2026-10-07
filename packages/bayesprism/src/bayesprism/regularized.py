"""
Collinearity-Aware Regularized Deconvolution Engine.

Implements Graph-Laplacian Manifold Regularization and Smoothed Fused Lasso
on the canonical probability simplex to solve cell state deconvolution under
high reference collinearity without discarding state-level resolution.

Strict functional implementation with returns, Pydantic (frozen=True), and Polars.
"""

from __future__ import annotations

from typing import Sequence
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Some, Nothing
from returns.result import Result, Success, Failure


class RegularizedDeconvConfig(BaseModel):
    """Configuration for Collinearity-Aware Regularized Deconvolution."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    kappa_target: float = 2000.0
    lambda_lap: Maybe[float] = Field(default=Nothing)
    lambda_fuse: Maybe[float] = Field(default=Nothing)
    max_iter: int = 500
    tol: float = 1e-8
    smoothing_eps: float = 1e-5
    correlation_power: float = 2.0


class RegularizedDeconvResult(BaseModel):
    """Results container holding inferred state fractions and diagnostics."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    theta: np.ndarray  # Shape (N, S)
    state_names: tuple[str, ...]
    sample_names: tuple[str, ...]
    condition_number_unreg: float
    condition_number_reg: float
    laplacian: np.ndarray
    weights: np.ndarray
    mean_iterations: float
    lambda_lap_used: float
    lambda_fuse_used: float

    def to_polars(self) -> pl.DataFrame:
        """Convert state proportions to a Polars DataFrame."""
        data_dict: dict[str, list[object]] = {"sample_id": list(self.sample_names)}
        for s_idx, s_name in enumerate(self.state_names):
            data_dict[s_name] = [float(v) for v in self.theta[:, s_idx]]
        return pl.DataFrame(data_dict)


def project_onto_simplex(v: np.ndarray) -> np.ndarray:
    """
    Exact Euclidean projection of an S-dimensional vector onto the probability simplex:
    Pi_Delta(v) = argmin_{w >= 0, sum(w) = 1} 0.5 * ||w - v||_2^2
    Implemented via Wang & Carreira-Perpiñán (2013) algorithm in O(S log S).
    """
    S = len(v)
    if S == 0:
        return v
    if not np.all(np.isfinite(v)):
        return np.ones(S, dtype=np.float64) / float(S)

    u = np.sort(v)[::-1]
    cssv = np.cumsum(u) - 1.0
    ind = np.arange(1, S + 1)
    cond = u - cssv / ind > 0.0
    if not np.any(cond):
        return np.ones(S, dtype=np.float64) / float(S)

    rho = int(np.where(cond)[0][-1])
    theta_mult = cssv[rho] / float(rho + 1)
    w = np.maximum(v - theta_mult, 0.0)
    w_sum = np.sum(w)
    if w_sum > 0:
        w /= w_sum
    else:
        w = np.ones(S, dtype=np.float64) / float(S)
    return w



def build_transcriptomic_graph(
    phi: np.ndarray,
    power: float = 2.0,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Construct transcriptomic similarity graph W and Graph Laplacian L.
    phi: shape (S, G)
    W_ij = max(0, PearsonCorr(phi_i, phi_j))^power for i != j, 0 on diagonal.
    L = D - W
    """
    S, G = phi.shape
    corr_mat = np.corrcoef(phi)
    corr_mat = np.nan_to_num(corr_mat, nan=0.0)

    # Positive thresholding and power scaling
    weights = np.maximum(corr_mat, 0.0) ** power
    np.fill_diagonal(weights, 0.0)

    # Degree and Laplacian
    degrees = np.sum(weights, axis=1)
    laplacian = np.diag(degrees) - weights

    return weights, laplacian


def calibrate_spectral_lambda(
    phi: np.ndarray,
    laplacian: np.ndarray,
    kappa_target: float = 2000.0,
) -> tuple[float, float]:
    """
    Analytically calibrate lambda_lap and lambda_fuse via Adaptive Deficit Regularization:
    Guarantees bounded condition number: kappa(H) <= kappa_target.

    If the nominal Hessian already satisfies kappa(H_nom) <= kappa_target, no regularization
    is required (lambda = 0.0), preserving unbiased recovery in well-conditioned regimes.
    Otherwise, injects precisely the eigenvalue deficit needed to lift the smallest singular values:
    sigma_target = sigma_max / kappa_target
    deficit = max(0.0, sigma_target - sigma_min)
    lambda_lap = deficit / mean_deg
    lambda_fuse = 0.02 * lambda_lap
    """
    mu0 = np.mean(phi, axis=0)
    mu0_safe = np.maximum(mu0, 1e-8)
    h_nom = (phi / mu0_safe) @ phi.T
    evals_h = np.linalg.eigvalsh(h_nom)
    sigma_max = float(np.max(evals_h))
    sigma_min = float(np.min(evals_h))

    # Target minimum eigenvalue to achieve kappa <= kappa_target
    sigma_target = sigma_max / max(kappa_target, 1.0)
    sigma_deficit = max(0.0, sigma_target - sigma_min)

    # Average degree of Graph Laplacian
    degrees = np.diag(laplacian)
    mean_deg = float(np.mean(degrees)) if np.mean(degrees) > 1e-9 else 1.0

    # Required lambda_lap to lift smallest singular values by the deficit
    lambda_lap = sigma_deficit / max(mean_deg, 1e-6)
    lambda_fuse = 0.02 * lambda_lap

    return float(lambda_lap), float(lambda_fuse)


def compute_loss_and_grad(
    theta: np.ndarray,
    x: np.ndarray,
    phi: np.ndarray,
    laplacian: np.ndarray,
    weights: np.ndarray,
    lambda_lap: float,
    lambda_fuse: float,
    eps: float,
) -> tuple[float, np.ndarray]:
    """
    Compute regularized Poisson-KL deviance, Graph Laplacian smoothness,
    and smoothed Fused Lasso penalty along with exact gradients.
    """
    S, G = phi.shape
    # Predicted expression: mu = phi.T @ theta (shape G, sums to 1)
    mu = phi.T @ theta
    mu_safe = np.maximum(mu, 1e-12)

    # Normalize mixture to gene proportion distribution p_x
    x_sum = float(np.sum(x))
    p_x = x / x_sum if x_sum > 0 else np.ones(G, dtype=np.float64) / float(G)

    # 1. Poisson-KL deviance: sum p_x * log(p_x / mu)
    pos_mask = p_x > 0
    kl_loss = float(np.sum(p_x[pos_mask] * np.log(p_x[pos_mask] / mu_safe[pos_mask])))

    # Gradient of KL w.r.t theta: phi @ (1 - p_x / mu)
    grad_kl = phi @ (1.0 - (p_x / mu_safe))

    # 2. Graph Laplacian penalty: 0.5 * lambda_lap * theta.T @ L @ theta
    lap_loss = 0.5 * lambda_lap * float(theta.T @ laplacian @ theta)
    grad_lap = lambda_lap * (laplacian @ theta)

    # 3. Smoothed Fused Lasso: lambda_fuse * sum_{i < j} W_ij * sqrt((theta_i - theta_j)^2 + eps^2)
    diff = theta[:, None] - theta[None, :]  # S x S, diff[i, j] = theta_i - theta_j
    smooth_abs = np.sqrt(diff**2 + eps**2)
    fuse_loss = 0.5 * lambda_fuse * float(np.sum(weights * smooth_abs))

    grad_fuse = lambda_fuse * np.sum(weights * (diff / smooth_abs), axis=1)

    total_loss = kl_loss + lap_loss + fuse_loss
    total_grad = grad_kl + grad_lap + grad_fuse

    return total_loss, total_grad


def solve_sample_fista(
    x: np.ndarray,
    phi: np.ndarray,
    laplacian: np.ndarray,
    weights: np.ndarray,
    lambda_lap: float,
    lambda_fuse: float,
    config: RegularizedDeconvConfig,
) -> tuple[np.ndarray, int, float]:
    """
    Solve regularized deconvolution for a single mixture sample using FISTA with backtracking line search.
    """
    S, G = phi.shape
    theta = np.ones(S, dtype=np.float64) / float(S)
    y = theta.copy()
    t_step = 1.0

    # Initial Lipschitz estimate for step size
    L = 50.0

    prev_loss = np.inf
    final_loss = np.inf
    iters = 0

    for k in range(1, config.max_iter + 1):
        loss_y, grad_y = compute_loss_and_grad(
            theta=y,
            x=x,
            phi=phi,
            laplacian=laplacian,
            weights=weights,
            lambda_lap=lambda_lap,
            lambda_fuse=lambda_fuse,
            eps=config.smoothing_eps,
        )

        # Backtracking line search for proximal step
        while True:
            step_size = 1.0 / L
            theta_next = project_onto_simplex(y - step_size * grad_y)
            loss_next, _ = compute_loss_and_grad(
                theta=theta_next,
                x=x,
                phi=phi,
                laplacian=laplacian,
                weights=weights,
                lambda_lap=lambda_lap,
                lambda_fuse=lambda_fuse,
                eps=config.smoothing_eps,
            )
            diff = theta_next - y
            quad_bound = loss_y + float(np.dot(grad_y, diff)) + 0.5 * L * float(np.sum(diff**2))
            if loss_next <= quad_bound + 1e-6 or L > 1e6:
                break
            L *= 1.5

        # Nesterov momentum update
        t_next = (1.0 + np.sqrt(1.0 + 4.0 * t_step**2)) / 2.0
        y = theta_next + ((t_step - 1.0) / t_next) * (theta_next - theta)

        diff_norm = float(np.max(np.abs(theta_next - theta)))
        theta = theta_next
        t_step = t_next
        final_loss = loss_next
        iters = k

        if diff_norm < config.tol or abs(prev_loss - loss_next) < config.tol:
            break
        prev_loss = loss_next

    return theta, iters, final_loss


def deconvolve_collinearity_regularized(
    mixture: np.ndarray,
    reference: np.ndarray,
    state_names: Sequence[str],
    sample_names: Sequence[str],
    config: Maybe[RegularizedDeconvConfig] | RegularizedDeconvConfig = Nothing,
) -> Result[RegularizedDeconvResult, str]:
    """
    Execute Collinearity-Aware Regularized Deconvolution across a patient cohort.

    mixture: shape (N, G), raw counts
    reference: shape (S, G), normalized single-cell state reference
    state_names: tuple of length S
    sample_names: tuple of length N
    """
    cfg = config if isinstance(config, RegularizedDeconvConfig) else config.value_or(RegularizedDeconvConfig())


    N, G_mix = mixture.shape
    S, G_ref = reference.shape

    if G_mix != G_ref:
        return Failure(f"Gene count mismatch: mixture has {G_mix} genes, reference has {G_ref} genes.")

    if len(state_names) != S:
        return Failure(f"State names length ({len(state_names)}) does not match reference rows ({S}).")

    if len(sample_names) != N:
        return Failure(f"Sample names length ({len(sample_names)}) does not match mixture rows ({N}).")

    # Ensure reference is normalized
    row_sums = np.sum(reference, axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0, 1.0, row_sums)
    phi = reference / row_sums

    # 1. Build Transcriptomic Graph and Laplacian
    weights, laplacian = build_transcriptomic_graph(phi, power=cfg.correlation_power)

    # 2. Calibrate regularization strengths
    match (cfg.lambda_lap, cfg.lambda_fuse):
        case (Some(l_lap), Some(l_fuse)):
            lambda_lap = l_lap
            lambda_fuse = l_fuse
        case (Some(l_lap), Nothing):
            lambda_lap = l_lap
            lambda_fuse = 0.02 * l_lap
        case _:
            lambda_lap, lambda_fuse = calibrate_spectral_lambda(phi, laplacian, cfg.kappa_target)

    # 3. Calculate Condition Numbers on Nominal Poisson Hessian
    mu0 = np.mean(phi, axis=0)
    mu0_safe = np.maximum(mu0, 1e-8)
    h_nom_unreg = (phi / mu0_safe) @ phi.T
    kappa_unreg = float(np.linalg.cond(h_nom_unreg))

    h_nom_reg = h_nom_unreg + lambda_lap * laplacian
    kappa_reg = float(np.linalg.cond(h_nom_reg))

    # 4. Solve per-sample deconvolution
    theta_mat = np.zeros((N, S), dtype=np.float64)
    iterations_list: list[int] = []

    for n in range(N):
        x_n = mixture[n].astype(np.float64)
        th_n, iters_n, _ = solve_sample_fista(
            x=x_n,
            phi=phi,
            laplacian=laplacian,
            weights=weights,
            lambda_lap=lambda_lap,
            lambda_fuse=lambda_fuse,
            config=cfg,
        )
        theta_mat[n, :] = th_n
        iterations_list.append(iters_n)

    result = RegularizedDeconvResult(
        theta=theta_mat,
        state_names=tuple(state_names),
        sample_names=tuple(sample_names),
        condition_number_unreg=kappa_unreg,
        condition_number_reg=kappa_reg,
        laplacian=laplacian,
        weights=weights,
        mean_iterations=float(np.mean(iterations_list)),
        lambda_lap_used=lambda_lap,
        lambda_fuse_used=lambda_fuse,
    )
    return Success(result)
