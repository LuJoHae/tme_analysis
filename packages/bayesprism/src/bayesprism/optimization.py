from typing import Union, Any
import numpy as np
import scipy.optimize as opt
import torch
from returns.result import Result, Success, Failure
from returns.maybe import Maybe, Some, Nothing
from bayesprism.models import (
    OptControl,
    RefPhi,
    RefTumor,
    Reference,
    JointPost,
)
from bayesprism.preprocessing import norm_to_one


def transform_phi_t(phi_t: np.ndarray, gamma_t: np.ndarray) -> np.ndarray:
    """
    Transform reference vector phi_t using log-fold-change gamma_t:
    psi_t = softmax(log(phi_t) + gamma_t)
    Stabilized by subtracting max(gamma_t).
    """
    stabilizing_constant = float(np.max(gamma_t))
    gamma_stab = gamma_t - stabilizing_constant
    psi_t = phi_t * np.exp(gamma_stab)
    sum_psi = float(np.sum(psi_t))
    if sum_psi == 0:
        sum_psi = 1.0
    return psi_t / sum_psi


def transform_phi(phi: np.ndarray, gamma: np.ndarray) -> np.ndarray:
    """Transform reference matrix phi row-wise using log-fold-change matrix gamma."""
    psi = np.zeros_like(phi)
    for t in range(phi.shape[0]):
        psi[t, :] = transform_phi_t(phi[t, :], gamma[t, :])
    return psi


def log_posterior_gamma(
    gamma_t: np.ndarray,
    phi_t: np.ndarray,
    phi_t_log: np.ndarray,
    Z_gt_t: np.ndarray,
    prior_num: float,
) -> float:
    """
    Compute negative log posterior over gamma_t for MAP optimization.
    """
    x = phi_t_log + gamma_t
    max_x = np.max(x)
    logsumexp_x = max_x + np.log(np.sum(np.exp(x - max_x)))
    psi_t_log = x - logsumexp_x

    log_likelihood = float(np.sum(Z_gt_t * psi_t_log))
    log_prior = float(np.sum(prior_num * (gamma_t**2)))
    log_posterior = log_likelihood + log_prior
    return -log_posterior


def log_posterior_gamma_grad(
    gamma_t: np.ndarray,
    phi_t: np.ndarray,
    phi_t_log: np.ndarray,
    Z_gt_t: np.ndarray,
    prior_num: float,
) -> np.ndarray:
    """
    Compute gradient of negative log posterior over gamma_t for MAP optimization.
    """
    psi_t = transform_phi_t(phi_t, gamma_t)
    Z_t_t = float(np.sum(Z_gt_t))

    log_likelihood_grad = Z_gt_t - (Z_t_t * psi_t)
    log_prior_grad = 2.0 * prior_num * gamma_t
    log_posterior_grad = log_likelihood_grad + log_prior_grad
    return -log_posterior_grad


def log_mle_gamma(
    gamma: np.ndarray,
    phi: np.ndarray,
    Z_gt: np.ndarray,
) -> float:
    """
    Compute negative log likelihood over a single gamma across cell types for MLE optimization.
    """
    phi_log = np.log(phi + 1e-12)
    x = phi_log + gamma[None, :]  # K x G
    max_x = np.max(x, axis=1, keepdims=True)
    logsumexp_x = max_x + np.log(np.sum(np.exp(x - max_x), axis=1, keepdims=True))
    psi_log = x - logsumexp_x
    log_likelihood = float(np.sum(Z_gt.T * psi_log))
    return -log_likelihood


def log_mle_gamma_grad(
    gamma: np.ndarray,
    phi: np.ndarray,
    Z_gt: np.ndarray,
) -> np.ndarray:
    """
    Compute gradient of negative log likelihood for MLE optimization.
    """
    psi = transform_phi(phi, np.tile(gamma, (phi.shape[0], 1)))
    Z_t = np.sum(Z_gt, axis=0, keepdims=True).T  # K x 1
    grad = np.sum(Z_gt.T - (Z_t * psi), axis=0)
    return -grad


def optimize_psi_map(
    phi: np.ndarray,
    Z_gt: np.ndarray,
    prior_num: float,
    opt_control: OptControl,
) -> tuple[np.ndarray, float]:
    """Optimize psi for each cell type independently using MAP."""
    K, G = phi.shape
    opt_gamma = np.zeros((K, G), dtype=np.float64)
    total_val = 0.0

    phi_log = np.log(np.maximum(phi, 1e-12))

    for t in range(K):
        res = opt.minimize(
            fun=log_posterior_gamma,
            x0=np.zeros(G, dtype=np.float64),
            args=(phi[t, :], phi_log[t, :], Z_gt[:, t], prior_num),
            method="L-BFGS-B",
            jac=log_posterior_gamma_grad,
            options={"maxiter": opt_control.maxit, "gtol": opt_control.eps},
        )
        gamma_t = res.x
        if np.max(np.abs(gamma_t)) > 20.0:
            gamma_t = np.zeros(G, dtype=np.float64)

        opt_gamma[t, :] = gamma_t
        total_val += float(res.fun)

    psi = transform_phi(phi, opt_gamma)
    return psi, total_val


def optimize_psi_mle(
    phi: np.ndarray,
    Z_gt: np.ndarray,
    opt_control: OptControl,
) -> tuple[np.ndarray, float]:
    """Optimize psi using a single gamma across cell types using MLE."""
    G = phi.shape[1]

    res = opt.minimize(
        fun=log_mle_gamma,
        x0=np.zeros(G, dtype=np.float64),
        args=(phi, Z_gt),
        method="L-BFGS-B",
        jac=log_mle_gamma_grad,
        options={"maxiter": opt_control.maxit, "gtol": opt_control.eps},
    )
    opt_gamma = res.x
    psi = transform_phi(phi, np.tile(opt_gamma, (phi.shape[0], 1)))
    return psi, float(res.fun)


def get_mle_psi_mal(Z_ng_mal: np.ndarray, pseudo_min: float) -> np.ndarray:
    """Compute MLE for malignant cell reference in each bulk sample."""
    row_sums = Z_ng_mal.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0, 1.0, row_sums)
    mle_psi = Z_ng_mal / row_sums

    # Apply norm_to_one
    mle_psi_tensor = norm_to_one(torch.from_numpy(mle_psi).float(), pseudo_min=pseudo_min)
    return mle_psi_tensor.numpy()


def update_reference(
    Z: torch.Tensor,
    phi_prime: RefPhi,
    state_to_type_map: dict[str, tuple[str, ...]],
    key: Maybe[str],
    opt_control: OptControl,
) -> Result[Union[RefPhi, RefTumor], str]:
    """
    Update reference matrix based on initial Gibbs sampling results.
    """
    Z_np = Z.detach().cpu().numpy()
    phi_prime_np = (
        phi_prime.phi.detach().cpu().numpy()
        if isinstance(phi_prime.phi, torch.Tensor)
        else np.asarray(phi_prime.phi)
    )

    sigma = opt_control.sigma
    prior_num = -1.0 / (2.0 * (sigma**2))

    if key == Nothing:
        # Non-malignant mode: optimize all cell types
        Z_gt = np.sum(Z_np, axis=0)  # G x K

        match opt_control.optimizer:
            case "MAP":
                psi, _ = optimize_psi_map(phi_prime_np, Z_gt, prior_num, opt_control)
            case "MLE":
                psi, _ = optimize_psi_mle(phi_prime_np, Z_gt, opt_control)
            case _:
                return Failure(f"Unknown optimizer: {opt_control.optimizer}")

        ref = RefPhi(
            phi=torch.from_numpy(psi).float(),
            cell_names=phi_prime.cell_names,
            gene_names=phi_prime.gene_names,
            pseudo_min=phi_prime.pseudo_min,
        )
        return Success(ref)

    else:
        key_name = key.unwrap()
        cell_types = phi_prime.cell_names
        if key_name not in cell_types:
            return Failure(f"Key {key_name} not found in cell types.")

        key_idx = cell_types.index(key_name)

        # Malignant expression profiles
        Z_ng_mal = Z_np[:, :, key_idx]
        psi_mal = get_mle_psi_mal(Z_ng_mal, pseudo_min=phi_prime.pseudo_min)

        # Environment expression profiles
        env_indices = [i for i, name in enumerate(cell_types) if name != key_name]
        env_cell_names = tuple(cell_types[i] for i in env_indices)

        Z_gt_env = np.sum(Z_np[:, :, env_indices], axis=0)
        phi_env = phi_prime_np[env_indices, :]

        match opt_control.optimizer:
            case "MAP":
                psi_env, _ = optimize_psi_map(phi_env, Z_gt_env, prior_num, opt_control)
            case "MLE":
                psi_env, _ = optimize_psi_mle(phi_env, Z_gt_env, opt_control)
            case _:
                return Failure(f"Unknown optimizer: {opt_control.optimizer}")

        bulk_names = tuple(f"sample_{n}" for n in range(Z_np.shape[0]))

        ref_tumor = RefTumor(
            psi_mal=torch.from_numpy(psi_mal).float(),
            psi_env=torch.from_numpy(psi_env).float(),
            key=key_name,
            bulk_names=bulk_names,
            env_cell_names=env_cell_names,
            gene_names=phi_prime.gene_names,
            pseudo_min=phi_prime.pseudo_min,
        )
        return Success(ref_tumor)
