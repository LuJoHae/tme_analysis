#!/usr/bin/env python3
"""
Pure Functional Synthetic Benchmark Data Generators.

Provides first-principles in silico data generation for deconvolution benchmarking:
- Synthetic reference signature matrices with controlled intra-lineage collinearity r.
- Ground truth cell state proportion vectors on the probability simplex.
- Discrete bulk mixture read count vectors via Multinomial sampling.

Adheres strictly to .agents/rules/code-style-guide.md:
- Pure functions, immutable Pydantic models / frozen dataclasses, no in-place mutation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final, Literal
import numpy as np


DEFAULT_HETERO_ALPHA: Final[tuple[float, ...]] = (
    1.2, 1.0, 0.8,  # Lineage 1
    1.0, 0.8, 0.6,  # Lineage 2
    1.2, 1.0, 0.8,  # Lineage 3
    0.8, 0.6, 0.4,  # Lineage 4
)


@dataclass(frozen=True)
class SyntheticSignature:
    """Immutable container holding reference expression signature and structural metadata."""

    phi: np.ndarray                     # (S, G) on probability simplex
    lineage_labels: tuple[str, ...]     # Length S (e.g. 'Lineage_1', ...)
    state_names: tuple[str, ...]        # Length S (e.g. 'L1_State_1', ...)
    gene_names: tuple[str, ...]         # Length G (e.g. 'Gene_000', ...)
    aggregation_matrix: np.ndarray      # (T, S) binary aggregation matrix M
    target_r: float                     # Requested target collinearity
    empirical_r: float                  # Actual mean intra-lineage correlation


def generate_synthetic_reference(
    n_lineages: int = 4,
    states_per_lineage: int = 3,
    n_genes: int = 400,
    target_r: float = 0.99,
    gamma_shape: float = 2.0,
    gamma_scale: float = 1.0,
    rng: np.random.Generator | None = None,
    seed: int = 42,
) -> SyntheticSignature:
    """
    Purely generate a synthetic reference matrix with controlled intra-lineage collinearity.

    Parameters:
    - n_lineages: Number of distinct cell lineages T (default 4).
    - states_per_lineage: Number of fine states per lineage (default 3, total S = 12).
    - n_genes: Number of genes G (default 400).
    - target_r: Target Pearson correlation r in [0.0, 0.99] between sibling states.
    - gamma_shape, gamma_scale: Parameters for the baseline Gamma distribution.
    - rng: Optional numpy Generator instance.

    Returns:
    - SyntheticSignature with simplex-normalized phi and metadata.
    """
    gen = rng if rng is not None else np.random.default_rng(seed)
    T = n_lineages
    S = T * states_per_lineage
    G = n_genes

    lineage_names = tuple(f"Lineage_{t+1}" for t in range(T))
    state_names = tuple(f"L{t+1}_State_{s+1}" for t in range(T) for s in range(states_per_lineage))
    lineage_labels = tuple(f"Lineage_{t+1}" for t in range(T) for _ in range(states_per_lineage))
    gene_names = tuple(f"Gene_{g:03d}" for g in range(G))

    # Binary aggregation matrix M mapping S fine states to T broad lineages
    M = np.zeros((T, S), dtype=np.float64)
    for t in range(T):
        for s in range(states_per_lineage):
            M[t, t * states_per_lineage + s] = 1.0

    # 1. Base profiles for each lineage
    base_profiles = [gen.gamma(gamma_shape, gamma_scale, G) for _ in range(T)]

    # 2. State profiles via linear superposition
    phi_raw = np.zeros((S, G), dtype=np.float64)
    for t in range(T):
        b_prof = base_profiles[t]
        for s in range(states_per_lineage):
            state_idx = t * states_per_lineage + s
            u_prof = gen.gamma(gamma_shape, gamma_scale, G)
            s_prof = np.sqrt(target_r) * b_prof + np.sqrt(1.0 - target_r) * u_prof
            phi_raw[state_idx] = s_prof

    # Simplex normalization: clip at 1e-6 and row-sum to 1.0
    phi_clipped = np.clip(phi_raw, 1e-6, None)
    phi = phi_clipped / np.sum(phi_clipped, axis=1, keepdims=True)

    # Compute empirical intra-lineage sibling correlation
    sibling_corrs: list[float] = []
    for t in range(T):
        for a in range(states_per_lineage):
            for b in range(a + 1, states_per_lineage):
                idx_a = t * states_per_lineage + a
                idx_b = t * states_per_lineage + b
                r_val = float(np.corrcoef(phi[idx_a], phi[idx_b])[0, 1])
                sibling_corrs.append(r_val)
    empirical_r = float(np.mean(sibling_corrs)) if sibling_corrs else target_r

    return SyntheticSignature(
        phi=phi,
        lineage_labels=lineage_labels,
        state_names=state_names,
        gene_names=gene_names,
        aggregation_matrix=M,
        target_r=target_r,
        empirical_r=empirical_r,
    )


def sample_true_proportions(
    n_samples: int,
    signature: SyntheticSignature,
    mode: Literal["skewed", "balanced", "state_dropout", "lineage_dropout", "rare_subpopulations"] = "skewed",
    alpha_within: float = 20.0,
    rng: np.random.Generator | None = None,
    seed: int = 42,
) -> np.ndarray:
    """
    Purely sample ground truth cell state proportions on the probability simplex.

    Modes:
    - 'skewed': Asymmetric Dirichlet proportions across 12 states (alpha = [1.2, 1.0, 0.8, ...]).
    - 'balanced': Lineages ~ Dir(1), sibling states ~ Dir(alpha_within) modeling continuous co-occurrence.
    - 'state_dropout': Sibling state L1_State_3 has exact true zero proportion (theta* = 0).
    - 'lineage_dropout': Entire Lineage_4 has exact true zero proportion (theta* = 0).
    - 'rare_subpopulations': Dominant lineage 1 is 75-85%, remaining fine states are 0.5% - 3.5%.

    Returns:
    - (N, S) matrix on the probability simplex (row-sums = 1.0).
    """
    gen = rng if rng is not None else np.random.default_rng(seed)
    S = signature.phi.shape[0]
    T = signature.aggregation_matrix.shape[0]
    states_per_lineage = S // T

    match mode:
        case "skewed":
            alpha = np.array(DEFAULT_HETERO_ALPHA[:S], dtype=np.float64)
            theta = gen.dirichlet(alpha, size=n_samples)
            return theta

        case "balanced":
            true_lineage = gen.dirichlet(np.ones(T, dtype=np.float64), size=n_samples)
            theta = np.zeros((n_samples, S), dtype=np.float64)
            for n in range(n_samples):
                for t in range(T):
                    alpha_w = np.full(states_per_lineage, alpha_within, dtype=np.float64)
                    w_sibling = gen.dirichlet(alpha_w)
                    theta[n, t * states_per_lineage : (t + 1) * states_per_lineage] = (
                        true_lineage[n, t] * w_sibling
                    )
            return theta

        case "state_dropout":
            # L1_State_3 (index 2) is strictly 0.0
            base_alpha = list(DEFAULT_HETERO_ALPHA[:S])
            base_alpha[2] = 0.0
            # Sample remaining 11 states
            active_indices = [i for i in range(S) if i != 2]
            active_alpha = np.array([base_alpha[i] for i in active_indices], dtype=np.float64)
            sub_theta = gen.dirichlet(active_alpha, size=n_samples)

            theta = np.zeros((n_samples, S), dtype=np.float64)
            for col_idx, active_i in enumerate(active_indices):
                theta[:, active_i] = sub_theta[:, col_idx]
            return theta

        case "lineage_dropout":
            # Lineage 4 (indices 9, 10, 11) is strictly 0.0
            active_indices = list(range(9))
            active_alpha = np.array([DEFAULT_HETERO_ALPHA[i] for i in active_indices], dtype=np.float64)
            sub_theta = gen.dirichlet(active_alpha, size=n_samples)

            theta = np.zeros((n_samples, S), dtype=np.float64)
            for col_idx, active_i in enumerate(active_indices):
                theta[:, active_i] = sub_theta[:, col_idx]
            return theta

        case "rare_subpopulations":
            # Lineage 1 is dominant (75% - 85% of total mixture)
            # Other 9 states share remaining 15% - 25% (individual states ~0.5% - 3.5%)
            theta = np.zeros((n_samples, S), dtype=np.float64)
            for n in range(n_samples):
                dom_frac = float(gen.uniform(0.75, 0.85))
                rare_frac = 1.0 - dom_frac
                w_dom = gen.dirichlet(np.array([2.0, 1.5, 1.0], dtype=np.float64))
                theta[n, :3] = dom_frac * w_dom
                w_rare = gen.dirichlet(np.ones(S - 3, dtype=np.float64))
                theta[n, 3:] = rare_frac * w_rare
            return theta


def generate_bulk_mixtures(
    true_theta: np.ndarray,
    phi: np.ndarray,
    n_total: int,
    drift_sd: float = 0.0,
    dispersion: float = 0.0,
    rng: np.random.Generator | None = None,
    seed: int = 42,
) -> np.ndarray:
    """
    Purely generate discrete bulk sequencing read counts via Multinomial or Negative Binomial sampling.

    Parameters:
    - true_theta: (N, S) matrix of true cell state proportions.
    - phi: (S, G) matrix of reference gene expression probabilities.
    - n_total: Multinomial total count parameter (number of trials / sample size).
    - drift_sd: Standard deviation of log-normal patient-specific expression drift (default 0.0).
    - dispersion: Negative binomial overdispersion alpha (default 0.0 for pure multinomial).
    - rng: Optional numpy Generator instance.

    Returns:
    - (N, G) matrix of non-negative integer read counts (float64).
    """
    gen = rng if rng is not None else np.random.default_rng(seed)
    N = true_theta.shape[0]
    S, G = phi.shape

    mixtures = np.zeros((N, G), dtype=np.float64)
    for n in range(N):
        if drift_sd > 1e-6:
            # Multiplicative lognormal patient-specific expression shift
            drift = gen.lognormal(mean=0.0, sigma=drift_sd, size=(S, G))
            phi_n_raw = phi * drift
            phi_n = phi_n_raw / np.sum(phi_n_raw, axis=1, keepdims=True)
            prob = true_theta[n] @ phi_n
            prob = prob / np.sum(prob)
        else:
            prob = true_theta[n] @ phi

        if dispersion > 1e-6:
            # Negative Binomial sampling via Gamma-Poisson mixture
            # Var(Y) = mu + dispersion * mu^2
            mu = np.clip(n_total * prob, 1e-8, None)
            k = 1.0 / dispersion
            scale = dispersion * mu
            lambdas = gen.gamma(shape=k, scale=scale)
            mixtures[n] = gen.poisson(lambdas).astype(np.float64)
        else:
            # Standard Multinomial counting noise
            mixtures[n] = gen.multinomial(n_total, prob).astype(np.float64)

    return mixtures
