#!/usr/bin/env python3
"""
Unified 3-Factorial Synthetic Deconvolution Benchmark Experiment.

Systematically decomposes total deconvolution error (MSE) into:
    MSE = Bias^2 + Variance
across a comprehensive factorial parameter grid:
1. Factor 1 (Collinearity r): [0.0, 0.6, 0.9, 0.99] (Condition numbers κ: ~33 to ~3,778)
2. Factor 2 (Count Parameter n_total): [2,500, 10,000, 80,000] (Poisson sampling shot noise)
3. Factor 3 (Biological State Architecture):
   - 'balanced': Intra-lineage phenotypic plasticity (α_within = 20)
   - 'skewed': Heterogeneous asymmetric proportions (α = [1.2, 1.0, 0.8, ...])
   - 'state_dropout': Structural biological zero in sibling state (θ* = 0)

Evaluates ALL 7 deconvolution tools concurrently:
1. Unregularized (NNLS)
2. RegDeconv (Graph Lap)
3. Rectangle (DWLS-QP)
4. CIBERSORT (reimpl., nu-SVR)
5. CIBERSORTx (Docker)
6. InstaPrism
7. BayesPrism (Gibbs)

Adheres strictly to .agents/rules/code-style-guide.md:
- Pure functional core, immutable models, Polars, and Returns.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import numpy as np
import polars as pl
from returns.maybe import Some, Nothing

from .synthetic_data import (
    SyntheticSignature,
    generate_synthetic_reference,
    sample_true_proportions,
    generate_bulk_mixtures,
)
from .deconv_runners import (
    ALL_METHODS,
    run_deconvolution_suite,
)
from .metrics import (
    compute_bias_variance_decomposition,
    compute_population_metrics,
    compute_phantom_detection,
)


def run_unified_experiment(
    collinearity_levels: tuple[float, ...] = (0.0, 0.6, 0.9, 0.99),
    count_levels: tuple[int, ...] = (2_500, 10_000, 80_000),
    architectures: tuple[str, ...] = ("balanced", "skewed", "state_dropout"),
    b_replicates: int = 30,
    n_samples: int = 30,
    n_genes: int = 400,
    out_dir: Path = Path("output/concordance"),
    cibersortx_username: str = "",
    cibersortx_token: str = "",
    seed: int = 42,
) -> int:
    """Execute the full 3-factorial synthetic deconvolution experiment."""
    rng = np.random.default_rng(seed)
    out_dir.mkdir(parents=True, exist_ok=True)

    creds = (
        Some((cibersortx_username, cibersortx_token))
        if cibersortx_username and cibersortx_token
        else Nothing
    )

    summary_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []

    print("=" * 80)
    print(" UNIFIED 3-FACTORIAL SYNTHETIC DECONVOLUTION BENCHMARK")
    print(f" Grid: {len(collinearity_levels)} Collinearity x {len(count_levels)} Count Levels x {len(architectures)} Architectures")
    print(f" Technical Replicates: B = {b_replicates} | Biological Samples: N = {n_samples}")
    print(f" Deconvolution Tools: {len(ALL_METHODS)} Methods (Including CIBERSORTx Docker)")
    print("=" * 80)

    for r_idx, target_r in enumerate(collinearity_levels):
        print(f"\n---> [Stage 1/4] Generating Reference Signature at Target Collinearity r = {target_r:.2f}...")
        signature = generate_synthetic_reference(
            n_lineages=4,
            states_per_lineage=3,
            n_genes=n_genes,
            target_r=target_r,
            rng=rng,
        )
        emp_r = signature.empirical_r
        S = signature.phi.shape[0]

        # Calculate nominal Hessian condition number kappa(H) at barycentric center
        theta_0 = np.ones(S) / float(S)
        mu_0 = theta_0 @ signature.phi
        H_nom = signature.phi @ np.diag(1.0 / np.clip(mu_0, 1e-12, None)) @ signature.phi.T
        eigs = np.linalg.eigvalsh(H_nom)
        kappa_val = float(eigs[-1] / max(eigs[0], 1e-12))

        for arch in architectures:
            print(f"  --> Architecture: {arch.upper()} | Condition Number kappa(H) = {kappa_val:,.1f}")

            # 1. Bias-Variance Decomposition (B replicates from fixed biological sample)
            # Sample single ground truth biological tissue
            fixed_theta = sample_true_proportions(1, signature, mode=arch, rng=rng)[0]  # (S,)

            for n_total in count_levels:
                print(f"    * Multinomial Count Parameter n_total = {n_total:,} counts | B = {b_replicates} Replicates...")
                # Generate B technical replicate mixtures from the same underlying tissue
                rep_thetas = np.tile(fixed_theta, (b_replicates, 1))
                rep_mixtures = generate_bulk_mixtures(rep_thetas, signature.phi, n_total=n_total, rng=rng)

                # Deconvolve all B replicates with all 7 tools
                rep_estimates = run_deconvolution_suite(
                    mixture=rep_mixtures,
                    signature=signature,
                    methods=ALL_METHODS,
                    cibersortx_credentials=creds,
                )

                # Compute exact Bias-Variance metrics for each tool
                bv_metrics = {
                    m: compute_bias_variance_decomposition(m, rep_estimates[m], fixed_theta, signature.aggregation_matrix)
                    for m in ALL_METHODS
                }

                # 2. Population Diversity & Ground Truth Correlation (N distinct samples)
                pop_thetas = sample_true_proportions(n_samples, signature, mode=arch, rng=rng)
                pop_mixtures = generate_bulk_mixtures(pop_thetas, signature.phi, n_total=n_total, rng=rng)

                pop_estimates = run_deconvolution_suite(
                    mixture=pop_mixtures,
                    signature=signature,
                    methods=ALL_METHODS,
                    cibersortx_credentials=creds,
                )

                for m in ALL_METHODS:
                    bv = bv_metrics[m]
                    pop_m = compute_population_metrics(pop_estimates[m], pop_thetas, signature)

                    # Compute phantom mass if architecture has structural zeros
                    phantom_mass = 0.0
                    if arch == "state_dropout":
                        true_zero_mask = (pop_thetas == 0.0)
                        phantom_mass = compute_phantom_detection(pop_estimates[m], true_zero_mask)

                    # Cross-sample Pearson correlation with ground truth
                    y_t = pop_thetas.ravel()
                    y_e = pop_estimates[m].ravel()
                    std_e = float(np.std(y_e))
                    pearson_r = float(np.corrcoef(y_t, y_e)[0, 1]) if std_e > 1e-8 else 0.0

                    summary_rows.append({
                        "target_r": target_r,
                        "collinearity_r": emp_r,
                        "condition_number": kappa_val,
                        "count_parameter_n": n_total,
                        "architecture": arch,
                        "method": m,
                        "state_mse": bv.state_mse,
                        "squared_bias": bv.squared_bias,
                        "variance": bv.variance,
                        "variance_pct": bv.variance_pct,
                        "bias_pct": bv.bias_pct,
                        "lineage_mse": bv.lineage_mse,
                        "population_mse": pop_m["state_mse"],
                        "sibling_cross_talk": pop_m["sibling_corr"],
                        "dropout_rate_pct": pop_m["dropout_rate_pct"],
                        "phantom_mass_pct": phantom_mass,
                        "pearson_r": pearson_r,
                    })

    # Persist summary Parquet
    summary_df = pl.DataFrame(summary_rows)
    out_path = out_dir / "unified_synthetic_bias_variance_summary.parquet"
    summary_df.write_parquet(out_path)
    print("\n" + "=" * 80)
    print(f"[OK] Experiment Complete! Summary saved to {out_path} ({summary_df.height} rows)")
    print("=" * 80)

    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description="Unified 3-Factorial Synthetic Deconvolution Benchmark")
    parser.add_argument("--collinearities", nargs="+", type=float, default=[0.0, 0.6, 0.9, 0.99])
    parser.add_argument("--counts", nargs="+", type=int, default=[2500, 10000, 80000])
    parser.add_argument("--architectures", nargs="+", type=str, default=["balanced", "skewed", "state_dropout"])
    parser.add_argument("--b-replicates", type=int, default=30)
    parser.add_argument("--n-samples", type=int, default=30)
    parser.add_argument("--cibersortx-username", type=str, default="")
    parser.add_argument("--cibersortx-token", type=str, default="")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    run_unified_experiment(
        collinearity_levels=tuple(args.collinearities),
        count_levels=tuple(args.counts),
        architectures=tuple(args.architectures),
        b_replicates=args.b_replicates,
        n_samples=args.n_samples,
        cibersortx_username=args.cibersortx_username,
        cibersortx_token=args.cibersortx_token,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
