#!/usr/bin/env python3
"""
Pure Functional Statistical Metrics for Deconvolution Benchmarking.

Provides exact bias-variance decomposition, state/lineage MSE, sibling cross-talk
correlation, spurious dropout rates, and false-positive phantom detection mass.

Adheres strictly to .agents/rules/code-style-guide.md:
- Pure functions, immutability, frozen dataclasses.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from .synthetic_data import SyntheticSignature


@dataclass(frozen=True)
class BiasVarianceMetrics:
    """Exact decomposition of estimator error into squared bias and variance."""

    method: str
    state_mse: float
    squared_bias: float
    variance: float
    variance_pct: float
    bias_pct: float
    lineage_mse: float


def compute_bias_variance_decomposition(
    method_name: str,
    estimated_replicates: np.ndarray,  # (B, S) across B sequencing draws
    true_theta: np.ndarray,            # (S,) or (1, S) fixed ground truth
    aggregation_matrix: np.ndarray,   # (T, S) binary lineage matrix
) -> BiasVarianceMetrics:
    """
    Compute exact bias-variance decomposition: MSE = Bias^2 + Variance.

    Parameters:
    - method_name: Identifier for the deconvolution algorithm.
    - estimated_replicates: (B, S) matrix of proportions across B technical replicate mixtures.
    - true_theta: True cell state proportion vector for this biological sample.
    - aggregation_matrix: (T, S) matrix mapping S states to T broad lineages.

    Returns:
    - BiasVarianceMetrics with exact mathematical components.
    """
    true_1d = np.squeeze(true_theta)
    # Expected prediction across technical replicate draws
    mean_pred = np.mean(estimated_replicates, axis=0)

    # 1. Squared Bias: || E[theta_hat] - theta* ||^2 / S
    sq_bias = float(np.mean((mean_pred - true_1d) ** 2))

    # 2. Estimator Variance: E[ || theta_hat - E[theta_hat] ||^2 ] / S
    # Notice: mean over states of sample variance
    var_per_state = np.var(estimated_replicates, axis=0, ddof=0)
    variance = float(np.mean(var_per_state))

    # Total state MSE = sq_bias + variance
    total_mse = sq_bias + variance
    var_pct = (variance / total_mse * 100.0) if total_mse > 1e-12 else 0.0
    bias_pct = (sq_bias / total_mse * 100.0) if total_mse > 1e-12 else 0.0

    # Lineage-level MSE across replicates
    true_lineage = true_1d @ aggregation_matrix.T
    pred_lineage = estimated_replicates @ aggregation_matrix.T
    lineage_mse = float(np.mean((pred_lineage - true_lineage) ** 2))

    return BiasVarianceMetrics(
        method=method_name,
        state_mse=total_mse,
        squared_bias=sq_bias,
        variance=variance,
        variance_pct=var_pct,
        bias_pct=bias_pct,
        lineage_mse=lineage_mse,
    )


def compute_population_metrics(
    estimated: np.ndarray,
    true_theta: np.ndarray,
    signature: SyntheticSignature,
    zero_threshold: float = 1e-4,
) -> dict[str, float]:
    """
    Compute population-level statistical metrics across N biological samples.

    Returns dictionary with:
    - 'state_mse': Mean Squared Error on fine state proportions
    - 'lineage_mse': Mean Squared Error on aggregated broad lineages
    - 'sibling_corr': Mean Pearson correlation between true sibling states
    - 'dropout_rate_pct': Percentage of true positive states collapsed to numerical zero
    """
    N, S = estimated.shape
    M = signature.aggregation_matrix
    T = M.shape[0]
    states_per_lineage = S // T

    # 1. State and Lineage MSE
    state_mse = float(np.mean((estimated - true_theta) ** 2))
    lineage_est = estimated @ M.T
    lineage_true = true_theta @ M.T
    lineage_mse = float(np.mean((lineage_est - lineage_true) ** 2))

    # 2. Sibling cross-talk correlation
    sibling_corrs: list[float] = []
    for t in range(T):
        for a in range(states_per_lineage):
            for b in range(a + 1, states_per_lineage):
                idx_a = t * states_per_lineage + a
                idx_b = t * states_per_lineage + b
                v_a = estimated[:, idx_a]
                v_b = estimated[:, idx_b]
                std_a = np.std(v_a)
                std_b = np.std(v_b)
                if std_a > 1e-9 and std_b > 1e-9:
                    r_val = float(np.corrcoef(v_a, v_b)[0, 1])
                    if not np.isnan(r_val):
                        sibling_corrs.append(r_val)
                else:
                    sibling_corrs.append(0.0)
    mean_sibling_corr = float(np.mean(sibling_corrs)) if sibling_corrs else 0.0

    # 3. Spurious dropout rate (% states collapsed to zero when ground truth > 0.01)
    active_mask = true_theta > 0.01
    if np.any(active_mask):
        collapsed_zeros = (estimated[active_mask] < zero_threshold)
        dropout_pct = float(np.mean(collapsed_zeros) * 100.0)
    else:
        dropout_pct = 0.0

    return {
        "state_mse": state_mse,
        "lineage_mse": lineage_mse,
        "sibling_corr": mean_sibling_corr,
        "dropout_rate_pct": dropout_pct,
    }


def compute_phantom_detection(
    estimated: np.ndarray,
    true_zero_mask: np.ndarray,
) -> float:
    """
    Compute total percentage mass erroneously assigned to true structural zeros.

    Parameters:
    - estimated: (N, S) matrix of inferred cell fractions.
    - true_zero_mask: Boolean matrix where True indicates theta* == 0.0.

    Returns:
    - Mean percentage mass allocated to false positive states.
    """
    if not np.any(true_zero_mask):
        return 0.0
    return float(np.mean(np.sum(estimated * true_zero_mask, axis=1)) * 100.0)
