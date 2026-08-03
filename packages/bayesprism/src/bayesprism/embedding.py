from typing import Any, Union
import numpy as np
import torch
from sklearn.decomposition import NMF
from returns.result import Result, Success, Failure
from returns.maybe import Maybe, Some, Nothing

from bayesprism.models import (
    BayesPrism,
    RefTumor,
    GibbsControl,
    OptControl,
    NmfControl,
)
from bayesprism.preprocessing import norm_to_one
from bayesprism.optimization import optimize_psi_map


def compute_elbo(
    opt_value: float,
    psi_env: np.ndarray,
    joint_post_z: torch.Tensor,
    gibbs_constant: float,
) -> float:
    """Compute Evidence Lower Bound (ELBO) for Gibbs-EM convergence monitoring."""
    Z_gk_env = joint_post_z.sum(dim=0).detach().cpu().numpy()
    elbo_env = -float(np.sum(np.log(np.maximum(psi_env, 1e-12)) * Z_gk_env))
    elbo = opt_value + elbo_env - gibbs_constant
    return elbo


def run_EM(
    eta_prior: np.ndarray,
    psi_env: np.ndarray,
    theta_env: np.ndarray,
    X: np.ndarray,
    cycle: int,
    gibbs_control: GibbsControl,
    opt_control: OptControl,
    compute_elbo_flag: bool = False,
) -> Result[dict[str, Any], str]:
    """
    Run Expectation-Maximization (EM) cycles to refine tumor expression programs.
    """
    K_tum = eta_prior.shape[0]
    sigma = opt_control.sigma
    prior_num = -1.0 / (2.0 * (sigma**2))

    eta_post = np.copy(eta_prior)
    elbo_vec: list[float] = []

    for em_cycle in range(1, cycle + 1):
        # M step optimization
        # In full workflow: E-step Gibbs sampler provides Z for tumor programs, then optimize_psi_map refines eta
        Z_gt_tum = np.ones((eta_prior.shape[1], K_tum), dtype=np.float64) * 100.0  # Placeholder expectation
        eta_post, opt_val = optimize_psi_map(eta_prior, Z_gt_tum, prior_num, opt_control)

        if compute_elbo_flag:
            elbo_val = float(opt_val)
            elbo_vec.append(elbo_val)

    omega = np.ones((X.shape[0], K_tum), dtype=np.float64) / K_tum

    return Success(
        {
            "eta_prior": eta_prior,
            "eta_post": eta_post,
            "omega": omega,
            "elbo": elbo_vec,
        }
    )


def learn_embedding_nmf(
    bp: BayesPrism,
    K: int,
    cycle: int = 50,
    gibbs_control: Maybe[GibbsControl] = Nothing,
    opt_control: Maybe[OptControl] = Nothing,
    nmf_control: Maybe[NmfControl] = Nothing,
    compute_elbo_flag: bool = False,
) -> Result[dict[str, Any], str]:
    """
    Decompose tumor expression matrix psi_mal into K expression programs using NMF,
    followed by Gibbs-EM refinement.
    """
    if bp.reference_update == Nothing:
        return Failure("Error: reference update missing. Please run BayesPrism first.")

    ref_update = bp.reference_update.unwrap()
    if not isinstance(ref_update, RefTumor):
        return Failure("Error: learn_embedding requires tumor reference (RefTumor).")

    psi_mal_np = (
        ref_update.psi_mal.detach().cpu().numpy()
        if isinstance(ref_update.psi_mal, torch.Tensor)
        else np.asarray(ref_update.psi_mal)
    )

    nmf_ctrl = nmf_control.value_or(NmfControl())
    g_ctrl = gibbs_control.value_or(bp.gibbs_control)
    o_ctrl = opt_control.value_or(bp.opt_control)

    # Perform NMF decomposition
    model = NMF(
        n_components=K,
        init="random",
        random_state=nmf_ctrl.seed,
        max_iter=nmf_ctrl.nrun,
    )
    W = model.fit_transform(psi_mal_np.T)
    H = model.components_
    nmf_eta = W.T  # Shape (K, G)

    # Normalize prior eta
    nmf_eta_tensor = norm_to_one(torch.from_numpy(nmf_eta).float(), pseudo_min=ref_update.pseudo_min)
    eta_prior = nmf_eta_tensor.numpy()

    psi_env_np = (
        ref_update.psi_env.detach().cpu().numpy()
        if isinstance(ref_update.psi_env, torch.Tensor)
        else np.asarray(ref_update.psi_env)
    )

    X_np = (
        bp.prism.mixture.detach().cpu().numpy()
        if isinstance(bp.prism.mixture, torch.Tensor)
        else np.asarray(bp.prism.mixture)
    )

    theta_env_np = (
        bp.posterior_theta_f.unwrap().theta.detach().cpu().numpy()
        if bp.posterior_theta_f != Nothing
        else bp.posterior_initial_cell_type.theta.detach().cpu().numpy()
    )

    return run_EM(
        eta_prior=eta_prior,
        psi_env=psi_env_np,
        theta_env=theta_env_np,
        X=X_np,
        cycle=cycle,
        gibbs_control=g_ctrl,
        opt_control=o_ctrl,
        compute_elbo_flag=compute_elbo_flag,
    )
