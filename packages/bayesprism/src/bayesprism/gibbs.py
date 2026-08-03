from typing import Union, Any
import numpy as np
import torch

from returns.result import Result, Success, Failure
from returns.maybe import Maybe, Some, Nothing
from bayesprism.models import (
    GibbsControl,
    GibbsSampler,
    RefPhi,
    RefTumor,
    JointPost,
    ThetaPost,
)


def get_gibbs_idx(gibbs_control: GibbsControl) -> torch.Tensor:
    """Compute iteration indices to retain after burn-in and thinning."""
    chain_length = gibbs_control.chain_length
    burn_in = gibbs_control.burn_in
    thinning = gibbs_control.thinning

    all_idx = torch.arange(1, chain_length + 1)
    burned_idx = all_idx[burn_in:]
    thinned_idx = burned_idx[::thinning]
    return thinned_idx


def rdirichlet(alpha: torch.Tensor, device: torch.device) -> torch.Tensor:
    """Draw a single sample from a Dirichlet distribution with parameter vector alpha."""
    gamma_samples = torch.distributions.Gamma(alpha, 1.0).sample()
    total = gamma_samples.sum()
    total = torch.where(total == 0, torch.ones_like(total), total)
    return gamma_samples / total


def sample_Z_theta_n(
    X_n: torch.Tensor,
    phi: torch.Tensor,
    alpha: float,
    gibbs_idx: torch.Tensor,
    device: torch.device,
) -> dict[str, Any]:
    """
    Joint MCMC sampling of latent profile Z_n and cell proportions theta_n for sample n.
    X_n: shape (G,)
    phi: shape (K, G)
    """
    K, G = phi.shape
    max_iter = int(gibbs_idx.max().item())
    gibbs_set = set(gibbs_idx.tolist())

    theta_n_i = torch.full((K,), 1.0 / K, dtype=torch.float32, device=device)
    Z_n_i = torch.zeros((G, K), dtype=torch.float32, device=device)

    Z_n_sum = torch.zeros((G, K), dtype=torch.float32, device=device)
    theta_n_sum = torch.zeros((K,), dtype=torch.float32, device=device)
    theta_n2_sum = torch.zeros((K,), dtype=torch.float32, device=device)
    multinom_coef = 0.0

    alpha_vec = torch.full((K,), float(alpha), dtype=torch.float32, device=device)

    for i in range(1, max_iter + 1):
        # prob matrix: shape (K, G) -> P_{k,g} = phi_{k,g} * theta_{n,k}
        prob_mat = phi * theta_n_i.unsqueeze(1)
        prob_sum = prob_mat.sum(dim=0, keepdim=True)
        prob_sum = torch.where(prob_sum == 0, torch.ones_like(prob_sum), prob_sum)
        prob_mat_norm = (prob_mat / prob_sum).t()  # shape (G, K)

        # Sample Z for each gene g
        for g in range(G):
            gene_count = int(X_n[g].item())
            if gene_count > 0:
                p = prob_mat_norm[g]
                Z_n_i[g, :] = torch.distributions.Multinomial(gene_count, probs=p).sample()
            else:
                Z_n_i[g, :] = 0.0

        # Sample theta for sample n
        Z_nk_i = Z_n_i.sum(dim=0)
        theta_n_i = rdirichlet(Z_nk_i + alpha_vec, device=device)

        if i in gibbs_set:
            Z_n_sum += Z_n_i
            theta_n_sum += theta_n_i
            theta_n2_sum += theta_n_i**2
            lfact_nk = torch.lgamma(Z_nk_i + 1.0).sum().item()
            lfact_ng = torch.lgamma(Z_n_i + 1.0).sum().item()
            multinom_coef += (lfact_nk - lfact_ng)

    sample_size = len(gibbs_set)
    Z_n = Z_n_sum / sample_size
    theta_n = theta_n_sum / sample_size
    var_theta = (theta_n2_sum / sample_size) - (theta_n**2)
    var_theta = torch.clamp(var_theta, min=0.0)
    theta_cv_n = torch.sqrt(var_theta) / torch.where(theta_n == 0, torch.ones_like(theta_n), theta_n)
    constant = multinom_coef / sample_size

    return {
        "Z_n": Z_n,
        "theta_n": theta_n,
        "theta_cv_n": theta_cv_n,
        "constant": constant,
    }


def sample_theta_n(
    X_n: torch.Tensor,
    phi: torch.Tensor,
    alpha: float,
    gibbs_idx: torch.Tensor,
    device: torch.device,
) -> dict[str, Any]:
    """
    Sampling for theta_n only (during reference refinement step).
    """
    K, G = phi.shape
    max_iter = int(gibbs_idx.max().item())
    gibbs_set = set(gibbs_idx.tolist())

    theta_n_i = torch.full((K,), 1.0 / K, dtype=torch.float32, device=device)
    Z_n_i = torch.zeros((G, K), dtype=torch.float32, device=device)

    theta_n_sum = torch.zeros((K,), dtype=torch.float32, device=device)
    theta_n2_sum = torch.zeros((K,), dtype=torch.float32, device=device)

    alpha_vec = torch.full((K,), float(alpha), dtype=torch.float32, device=device)

    for i in range(1, max_iter + 1):
        prob_mat = phi * theta_n_i.unsqueeze(1)
        prob_sum = prob_mat.sum(dim=0, keepdim=True)
        prob_sum = torch.where(prob_sum == 0, torch.ones_like(prob_sum), prob_sum)
        prob_mat_norm = (prob_mat / prob_sum).t()

        for g in range(G):
            gene_count = int(X_n[g].item())
            if gene_count > 0:
                p = prob_mat_norm[g]
                Z_n_i[g, :] = torch.distributions.Multinomial(gene_count, probs=p).sample()
            else:
                Z_n_i[g, :] = 0.0

        Z_nk_i = Z_n_i.sum(dim=0)
        theta_n_i = rdirichlet(Z_nk_i + alpha_vec, device=device)

        if i in gibbs_set:
            theta_n_sum += theta_n_i
            theta_n2_sum += theta_n_i**2

    sample_size = len(gibbs_set)
    theta_n = theta_n_sum / sample_size
    var_theta = (theta_n2_sum / sample_size) - (theta_n**2)
    var_theta = torch.clamp(var_theta, min=0.0)
    theta_cv_n = torch.sqrt(var_theta) / torch.where(theta_n == 0, torch.ones_like(theta_n), theta_n)

    return {
        "theta_n": theta_n,
        "theta_cv_n": theta_cv_n,
    }


def run_gibbs(
    sampler: GibbsSampler,
    final: bool = False,
) -> Result[Union[JointPost, ThetaPost], str]:
    """
    Run Gibbs sampling for all bulk samples in the gibbsSampler object.
    """
    control = sampler.gibbs_control
    device = torch.device(control.device)

    if isinstance(control.seed, Some):
        seed_val = control.seed.unwrap()
        torch.manual_seed(seed_val)
        np.random.seed(seed_val)

    torch.set_num_threads(control.num_threads)
    gibbs_idx = get_gibbs_idx(control).to(device=device)

    X_mat = sampler.X
    if not isinstance(X_mat, torch.Tensor):
        X_tensor = torch.from_numpy(np.asarray(X_mat)).to(dtype=torch.float32, device=device)
    else:
        X_tensor = X_mat.to(dtype=torch.float32, device=device)

    N, G = X_tensor.shape

    ref = sampler.reference

    if isinstance(ref, RefPhi):
        phi_mat = ref.phi
        if not isinstance(phi_mat, torch.Tensor):
            phi_tensor = torch.from_numpy(np.asarray(phi_mat)).to(dtype=torch.float32, device=device)
        else:
            phi_tensor = phi_mat.to(dtype=torch.float32, device=device)

        K = phi_tensor.shape[0]
        bulk_names = tuple(f"sample_{n}" for n in range(N))
        gene_names = ref.gene_names
        cell_names = ref.cell_names

        if not final:
            Z_array = torch.zeros((N, G, K), dtype=torch.float32, device=device)
            theta_mat = torch.zeros((N, K), dtype=torch.float32, device=device)
            theta_cv_mat = torch.zeros((N, K), dtype=torch.float32, device=device)
            total_constant = 0.0

            for n in range(N):
                res = sample_Z_theta_n(
                    X_n=X_tensor[n],
                    phi=phi_tensor,
                    alpha=control.alpha,
                    gibbs_idx=gibbs_idx,
                    device=device,
                )
                Z_array[n] = res["Z_n"]
                theta_mat[n] = res["theta_n"]
                theta_cv_mat[n] = res["theta_cv_n"]
                total_constant += res["constant"]

            joint_post = JointPost(
                Z=Z_array,
                theta=theta_mat,
                theta_cv=theta_cv_mat,
                constant=total_constant,
                bulk_names=bulk_names,
                gene_names=gene_names,
                cell_names=cell_names,
            )
            return Success(joint_post)
        else:
            theta_mat = torch.zeros((N, K), dtype=torch.float32, device=device)
            theta_cv_mat = torch.zeros((N, K), dtype=torch.float32, device=device)

            for n in range(N):
                res = sample_theta_n(
                    X_n=X_tensor[n],
                    phi=phi_tensor,
                    alpha=control.alpha,
                    gibbs_idx=gibbs_idx,
                    device=device,
                )
                theta_mat[n] = res["theta_n"]
                theta_cv_mat[n] = res["theta_cv_n"]

            theta_post = ThetaPost(
                theta=theta_mat,
                theta_cv=theta_cv_mat,
                bulk_names=bulk_names,
                cell_names=cell_names,
            )
            return Success(theta_post)

    elif isinstance(ref, RefTumor):
        psi_mal = ref.psi_mal
        psi_env = ref.psi_env
        if not isinstance(psi_mal, torch.Tensor):
            psi_mal_t = torch.from_numpy(np.asarray(psi_mal)).to(dtype=torch.float32, device=device)
        else:
            psi_mal_t = psi_mal.to(dtype=torch.float32, device=device)

        if not isinstance(psi_env, torch.Tensor):
            psi_env_t = torch.from_numpy(np.asarray(psi_env)).to(dtype=torch.float32, device=device)
        else:
            psi_env_t = psi_env.to(dtype=torch.float32, device=device)

        cell_names = (ref.key,) + ref.env_cell_names
        K = len(cell_names)
        bulk_names = ref.bulk_names

        theta_mat = torch.zeros((N, K), dtype=torch.float32, device=device)
        theta_cv_mat = torch.zeros((N, K), dtype=torch.float32, device=device)

        for n in range(N):
            phi_n = torch.cat([psi_mal_t[n : n + 1, :], psi_env_t], dim=0)
            nonzero_mask = phi_n.max(dim=0).values > 0
            res = sample_theta_n(
                X_n=X_tensor[n, nonzero_mask],
                phi=phi_n[:, nonzero_mask],
                alpha=control.alpha,
                gibbs_idx=gibbs_idx,
                device=device,
            )
            theta_mat[n] = res["theta_n"]
            theta_cv_mat[n] = res["theta_cv_n"]

        theta_post = ThetaPost(
            theta=theta_mat,
            theta_cv=theta_cv_mat,
            bulk_names=bulk_names,
            cell_names=cell_names,
        )
        return Success(theta_post)

    else:
        return Failure("Error: unknown reference object class.")
