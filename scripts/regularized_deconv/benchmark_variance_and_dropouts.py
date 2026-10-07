#!/usr/bin/env python3
"""
Diagnostic Benchmark: Estimator Variance & True Cell State/Type Dropouts.

1. Estimator Variance Stress-Test:
   - Evaluates deconvolution stability across B=30 technical sequencing replicates
     drawn from the SAME biological mixture at severe collinearity (r = 0.99).
   - Computes exact Bias-Variance Decomposition: MSE = Bias^2 + Variance.
   - Proves whether BayesPrism and InstaPrism trade fine-state correlation for
     near-zero estimator variance compared to the volatile fluctuations of NNLS.

2. True Cell State & Lineage Dropout Stress-Test:
   - Evaluates detection accuracy when states or entire lineages are TRULY ABSENT (θ* = 0.0).
   - Condition A (State Dropout): L1_State_3 is truly absent (θ* = 0.0).
   - Condition B (Lineage Dropout): Entire Lineage 4 is truly absent (θ* = 0.0).
   - Measures Phantom Detection Fraction (Dirichlet prior leakage on true zeros).

Evaluates all 7 deconvolution methods concurrently.
Adheres strictly to .agents/rules/code-style-guide.md: Polars, Pydantic, Returns.
"""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import polars as pl

try:
    from .synthetic_data import (
        generate_synthetic_reference,
        sample_true_proportions,
        generate_bulk_mixtures,
    )
    from .deconv_runners import ALL_METHODS, run_deconvolution_suite
    from .metrics import compute_bias_variance_decomposition
except ImportError:
    from regularized_deconv.synthetic_data import (
        generate_synthetic_reference,
        sample_true_proportions,
        generate_bulk_mixtures,
    )
    from regularized_deconv.deconv_runners import ALL_METHODS, run_deconvolution_suite
    from regularized_deconv.metrics import compute_bias_variance_decomposition


def run_variance_and_dropout_benchmark(
    n_replicates: int = 30,
    n_dropout_samples: int = 20,
    n_genes: int = 400,
    seed: int = 42,
    out_dir: Path = Path("output/concordance"),
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Run estimator variance test and true zero dropout stress-test across all 7 tools."""
    rng = np.random.default_rng(seed)
    out_dir.mkdir(parents=True, exist_ok=True)

    sig = generate_synthetic_reference(
        n_lineages=4,
        states_per_lineage=3,
        n_genes=n_genes,
        target_r=0.99,
        seed=seed,
    )
    S = sig.phi.shape[0]
    T = sig.aggregation_matrix.shape[0]

    # =========================================================================
    # PART 1: Estimator Variance across B=30 Technical Sequencing Replicates
    # =========================================================================
    print(f"[INFO] Running Part 1: Estimator Variance across B={n_replicates} replicates (r = 0.99)...")
    fixed_true_theta = sample_true_proportions(
        n_samples=1,
        signature=sig,
        mode="skewed",
        rng=rng,
    )[0]  # Shape (S,)

    # Generate B independent multinomial sequencing draws from the SAME tissue
    rep_theta = np.tile(fixed_true_theta, (n_replicates, 1))
    rep_mixtures = generate_bulk_mixtures(
        true_theta=rep_theta,
        phi=sig.phi,
        n_total=80_000,
        rng=rng,
    )
    rep_sample_names = tuple(f"Replicate_{b:02d}" for b in range(n_replicates))

    rep_results = run_deconvolution_suite(
        mixture=rep_mixtures,
        signature=sig,
        sample_names=rep_sample_names,
    )

    var_decomp_rows: list[dict[str, object]] = []
    for method_name in ALL_METHODS:
        est_mat = rep_results[method_name]  # Shape (B, S)
        decomp = compute_bias_variance_decomposition(
            method_name=method_name,
            estimated_replicates=est_mat,
            true_theta=fixed_true_theta,
            aggregation_matrix=sig.aggregation_matrix,
        )

        # Sibling cross-talk across replicates
        sib_corrs: list[float] = []
        for t in range(T):
            for a in range(3):
                for b in range(a + 1, 3):
                    ca = est_mat[:, t * 3 + a]
                    cb = est_mat[:, t * 3 + b]
                    if np.std(ca) > 1e-8 and np.std(cb) > 1e-8:
                        sib_corrs.append(float(np.corrcoef(ca, cb)[0, 1]))
                    else:
                        sib_corrs.append(0.0)
        rep_sib_corr = float(np.mean(sib_corrs)) if sib_corrs else 0.0

        var_decomp_rows.append({
            "method": method_name,
            "mean_variance": decomp.variance,
            "mean_bias_squared": decomp.squared_bias,
            "mean_mse": decomp.state_mse,
            "variance_pct_of_mse": decomp.variance_pct,
            "bias_pct_of_mse": decomp.bias_pct,
            "replicate_sibling_cross_talk": rep_sib_corr,
        })

    var_df = pl.DataFrame(var_decomp_rows)

    # =========================================================================
    # PART 2: True State & Lineage Dropout Stress-Test (Zero Proportions)
    # =========================================================================
    print(f"[INFO] Running Part 2: True Dropouts across N={n_dropout_samples} mixtures...")
    active_indices = [0, 1, 3, 4, 5, 6, 7, 8]  # 8 active states
    absent_state_idx = 2                        # L1_State_3 absent
    absent_lineage_indices = [9, 10, 11]       # Lineage_4 absent

    dropout_alpha = np.array([1.2, 1.0, 1.0, 0.8, 0.6, 1.2, 1.0, 0.8])
    dropout_true_theta = np.zeros((n_dropout_samples, S), dtype=np.float64)

    for i in range(n_dropout_samples):
        draw = rng.dirichlet(dropout_alpha)
        for sub_i, act_idx in enumerate(active_indices):
            dropout_true_theta[i, act_idx] = draw[sub_i]

    dropout_mixtures = generate_bulk_mixtures(
        true_theta=dropout_true_theta,
        phi=sig.phi,
        n_total=80_000,
        rng=rng,
    )
    dropout_sample_names = tuple(f"Dropout_Sample_{i:02d}" for i in range(n_dropout_samples))

    dropout_results = run_deconvolution_suite(
        mixture=dropout_mixtures,
        signature=sig,
        sample_names=dropout_sample_names,
    )

    dropout_metric_rows: list[dict[str, object]] = []
    for method_name in ALL_METHODS:
        est_mat = dropout_results[method_name]

        phantom_single_state = float(np.mean(est_mat[:, absent_state_idx]))
        phantom_lineage = float(np.mean(np.sum(est_mat[:, absent_lineage_indices], axis=1)))
        all_absent_indices = [absent_state_idx] + absent_lineage_indices
        total_phantom = float(np.mean(np.sum(est_mat[:, all_absent_indices], axis=1)))

        active_true = dropout_true_theta[:, active_indices]
        active_est = est_mat[:, active_indices]
        active_rmse = float(np.sqrt(np.mean((active_est - active_true) ** 2)))

        dropout_metric_rows.append({
            "method": method_name,
            "phantom_state_pct": phantom_single_state * 100.0,
            "phantom_lineage_pct": phantom_lineage * 100.0,
            "total_phantom_pct": total_phantom * 100.0,
            "active_states_rmse": active_rmse,
        })

    dropout_df = pl.DataFrame(dropout_metric_rows)

    var_df.write_parquet(out_dir / "synthetic_estimator_variance_summary.parquet")
    dropout_df.write_parquet(out_dir / "synthetic_dropout_stress_test_summary.parquet")

    print("\n" + "=" * 80)
    print(" PART 1: ESTIMATOR BIAS-VARIANCE DECOMPOSITION (B = 30 Replicates, r = 0.99)")
    print("=" * 80)
    print(var_df)

    print("\n" + "=" * 80)
    print(" PART 2: TRUE DROPOUT & PHANTOM DETECTION STRESS-TEST (N = 20 Samples)")
    print("=" * 80)
    print(dropout_df)

    return var_df, dropout_df


def main() -> None:
    run_variance_and_dropout_benchmark()


if __name__ == "__main__":
    main()
