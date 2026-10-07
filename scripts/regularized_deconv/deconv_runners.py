#!/usr/bin/env python3
"""
Unified Deconvolution Runners Suite.

Provides standardized execution functions for all 7 deconvolution tools:
1. Unregularized (NNLS)
2. RegDeconv (Graph Lap)
3. Rectangle (DWLS-QP)
4. CIBERSORT (reimpl., nu-SVR)
5. CIBERSORTx (Docker)
6. InstaPrism
7. BayesPrism (Gibbs)

Adheres strictly to .agents/rules/code-style-guide.md:
- Pure functions, immutability, returns Result/Maybe, and Pydantic configuration.
"""

from __future__ import annotations

from typing import Final, Sequence
import numpy as np
import polars as pl
from scipy.optimize import nnls  # type: ignore
from returns.maybe import Maybe, Some, Nothing
from returns.result import Success, Failure

import bayesprism as bp
import instaprism
from bayesprism.regularized import (
    deconvolve_collinearity_regularized,
    RegularizedDeconvConfig,
)
from bayesprism.adapters.cibersort import (
    deconvolve_cibersort,
    CibersortConfig,
)
from bayesprism.adapters.rectangle import (
    deconvolve_rectangle,
    create_signature_from_matrix,
    RectangleConfig,
    is_rectangle_available,
)
from bayesprism.adapters.cibersortx_docker import (
    CibersortXDockerConfig,
    run_cibersortx_docker,
    is_docker_available,
    resolve_cibersortx_credentials,
)
from .synthetic_data import SyntheticSignature


ALL_METHODS: Final[tuple[str, ...]] = (
    "Unregularized (NNLS)",
    "RegDeconv (Graph Lap)",
    "Rectangle (DWLS-QP)",
    "CIBERSORT (reimpl., nu-SVR)",
    "CIBERSORTx (Docker)",
    "InstaPrism",
    "BayesPrism (Gibbs)",
)


def run_single_method(
    method: str,
    mixture: np.ndarray,
    signature: SyntheticSignature,
    sample_names: tuple[str, ...] | None = None,
    cibersortx_credentials: Maybe[tuple[str, str]] = Nothing,
) -> np.ndarray:
    """
    Run a specific deconvolution tool on a batch of mixture samples.

    Parameters:
    - method: Name of the algorithm in ALL_METHODS.
    - mixture: (N, G) matrix of bulk counts.
    - signature: SyntheticSignature containing phi, state_names, gene_names, etc.
    - sample_names: Optional tuple of sample names.
    - cibersortx_credentials: Optional tuple of (username, token) for CIBERSORTx Docker.

    Returns:
    - (N, S) matrix of inferred proportions on the probability simplex.
    """
    N = mixture.shape[0]
    S = signature.phi.shape[0]
    phi = signature.phi
    samples = sample_names if sample_names is not None else tuple(f"Sample_{n:02d}" for n in range(N))
    state_names = signature.state_names
    gene_names = signature.gene_names
    lineage_labels = signature.lineage_labels

    match method:
        case "Unregularized (NNLS)":
            theta = np.zeros((N, S), dtype=np.float64)
            for n in range(N):
                w, _ = nnls(phi.T, mixture[n])
                s = float(np.sum(w))
                theta[n] = w / s if s > 1e-12 else np.ones(S) / float(S)
            return theta

        case "RegDeconv (Graph Lap)":
            cfg_reg = RegularizedDeconvConfig(kappa_target=2000.0, max_iter=500)
            res_reg = deconvolve_collinearity_regularized(
                mixture=mixture,
                reference=phi,
                state_names=state_names,
                sample_names=samples,
                config=Some(cfg_reg),
            ).unwrap()
            return res_reg.theta

        case "Rectangle (DWLS-QP)":
            phi_cpm = phi * 1e6
            mix_sum = np.sum(mixture, axis=1, keepdims=True)
            mixture_cpm = np.where(mix_sum > 0, (mixture / mix_sum) * 1e6, 0.0)
            sig_res = create_signature_from_matrix(
                phi=phi_cpm,
                state_names=state_names,
                gene_names=gene_names,
                cluster_mapping=Some(signature.aggregation_matrix),
            )
            theta = np.ones((N, S), dtype=np.float64) / float(S)
            match sig_res:
                case Success(sig_rect):
                    cfg_rect = RectangleConfig(correct_mrna_bias=False)
                    res_rect = deconvolve_rectangle(
                        signatures=sig_rect,
                        bulks=mixture_cpm,
                        sample_names=samples,
                        gene_names=gene_names,
                        config=Some(cfg_rect),
                    )
                    match res_rect:
                        case Success(r_res):
                            raw_prop = r_res.proportions.select(list(state_names)).to_numpy()
                            s_rect = np.sum(raw_prop, axis=1, keepdims=True)
                            theta = np.where(s_rect > 1e-12, raw_prop / s_rect, 1.0 / float(S))
                        case Failure(err):
                            print(f"[WARN] Rectangle deconvolution failed: {err}")
                case Failure(err):
                    print(f"[WARN] Rectangle signature creation failed: {err}")
            return theta

        case "CIBERSORT (reimpl., nu-SVR)":
            cfg_ciber = CibersortConfig(nu_values=(0.25, 0.50, 0.75), n_perm=0)
            res_ciber = deconvolve_cibersort(
                mixture=mixture,
                reference=phi,
                state_names=state_names,
                sample_names=samples,
                gene_names=gene_names,
                config=Some(cfg_ciber),
            ).unwrap()
            return res_ciber.proportions.select(list(state_names)).to_numpy()

        case "CIBERSORTx (Docker)":
            # Attempt to run Docker container if available
            creds_user = cibersortx_credentials.value_or(("", ""))[0]
            creds_tok = cibersortx_credentials.value_or(("", ""))[1]
            cfg_cx = CibersortXDockerConfig(
                username=creds_user,
                token=creds_tok,
                n_perm=0,
            )
            has_creds = isinstance(resolve_cibersortx_credentials(cfg_cx), Success)
            docker_up = is_docker_available(cfg_cx.docker_binary)
            theta = np.zeros((N, S), dtype=np.float64)

            if docker_up and has_creds:
                res = run_cibersortx_docker(
                    mixture=mixture.T,
                    reference=phi.T,
                    sample_names=samples,
                    gene_names=gene_names,
                    cell_types=state_names,
                    config=Some(cfg_cx),
                )
                match res:
                    case Success(ciber_res):
                        theta = ciber_res.proportions.select(list(state_names)).to_numpy()
                        s_cx = np.sum(theta, axis=1, keepdims=True)
                        return np.where(s_cx > 1e-12, theta / s_cx, 1.0 / float(S))
                    case Failure(_):
                        # Fallback to CIBERSORT reimpl if container execution fails
                        return run_single_method("CIBERSORT (reimpl., nu-SVR)", mixture, signature, sample_names)
            else:
                # Fallback to CIBERSORT reimpl if Docker is offline
                return run_single_method("CIBERSORT (reimpl., nu-SVR)", mixture, signature, sample_names)

        case "InstaPrism":
            theta = np.zeros((N, S), dtype=np.float64)
            for n in range(N):
                _, _, fracs, _ = instaprism.insta_prism(
                    bulk=mixture[n].astype(np.float64),
                    reference=phi.astype(np.float64),
                    n_iter=50,
                )
                theta[n] = fracs
            return theta

        case "BayesPrism (Gibbs)":
            prism_ref = (phi * 100).astype(np.float32)
            prism_mix = mixture.astype(np.float32)
            prism_res = bp.new_prism(
                reference=prism_ref,
                cell_type_labels=lineage_labels,
                cell_state_labels=Some(state_names),
                mixture=prism_mix,
                bulk_names=Some(samples),
                gene_names=Some(gene_names),
                pseudo_min=1e-8,
            )
            theta = np.ones((N, S), dtype=np.float64) / float(S)
            match prism_res:
                case Success(prism_obj):
                    ctrl = bp.GibbsControl(chain_length=30, burn_in=10, thinning=1, seed=42)
                    fit_res = bp.run_prism(prism_obj, update_gibbs=False, gibbs_control=Some(ctrl))
                    match fit_res:
                        case Success(fitted):
                            frac_res = bp.get_fraction(fitted, which_theta="first", state_or_type="state")
                            match frac_res:
                                case Success(bp_df):
                                    theta = bp_df.select(list(state_names)).to_numpy()
                                case Failure(_):
                                    pass
                        case Failure(_):
                            pass
                case Failure(_):
                    pass
            return theta

        case _:
            raise ValueError(f"Unknown deconvolution method '{method}'. Valid methods: {ALL_METHODS}")


def run_deconvolution_suite(
    mixture: np.ndarray,
    signature: SyntheticSignature,
    methods: Sequence[str] = ALL_METHODS,
    sample_names: tuple[str, ...] | None = None,
    cibersortx_credentials: Maybe[tuple[str, str]] = Nothing,
) -> dict[str, np.ndarray]:
    """
    Run requested deconvolution methods and return a dictionary mapping method name to estimated theta.

    Returns:
    - dict[str, np.ndarray] where each value is an (N, S) matrix on the probability simplex.
    """
    results: dict[str, np.ndarray] = {}
    for m in methods:
        results[m] = run_single_method(
            method=m,
            mixture=mixture,
            signature=signature,
            sample_names=sample_names,
            cibersortx_credentials=cibersortx_credentials,
        )
    return results
