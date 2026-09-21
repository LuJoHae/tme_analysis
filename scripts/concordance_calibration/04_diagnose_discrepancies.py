#!/usr/bin/env python3
"""
Step 4 (Calibration): Discrepancy Diagnostics and Resolution Trend Analysis.
Applies the 4-gate diagnostic decision tree to attribute discrepancies to:
- Deconvolution Collinear Leakage (Gate 1)
- mRNA Content / Cell Size Bias (Gate 2)
- Statistical Model Graph vs Cluster Discrepancy (Gate 3)
- Cross-Cohort Biological Heterogeneity (Gate 4)
And assesses concordance trends across resolutions from 0.5 up to 3.0.
Strict functional Python with returns, Pydantic, and Polars.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
import polars as pl
from scipy import stats  # type: ignore


class DiagnosticConfig(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    triangulation_path: Path
    fidelity_path: Path
    out_dir: Path


def parse_args() -> DiagnosticConfig:
    parser = argparse.ArgumentParser(
        description="Run discrepancy diagnostics across resolutions and modalities."
    )
    parser.add_argument(
        "--triangulation",
        type=Path,
        default=Path("output/concordance_calibration/triangulation_effect_estimates.parquet"),
        help="Path to triangulation_effect_estimates.parquet",
    )
    parser.add_argument(
        "--fidelity",
        type=Path,
        default=Path("output/concordance_calibration/deconv_fidelity_benchmark.parquet"),
        help="Path to deconv_fidelity_benchmark.parquet",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("output/concordance_calibration"),
        help="Output directory",
    )
    args = parser.parse_args()
    return DiagnosticConfig(
        triangulation_path=args.triangulation,
        fidelity_path=args.fidelity,
        out_dir=args.out_dir,
    )


def attribute_primary_cause(
    b_count: float,
    b_umi: float,
    b_deconv: float,
    b_milo: float,
    b_bulk: float,
    r_fidelity: float,
) -> str:
    """Classify the root cause of discrepancy for a cell state."""
    # Check if concordant across all 5 modalities
    signs = [np.sign(b) for b in [b_count, b_umi, b_deconv, b_milo, b_bulk] if abs(b) > 0.05]
    if len(signs) >= 4 and len(set(signs)) == 1:
        return "Robustly Concordant"

    # Gate 1: Deconvolution Identifiability Failure
    if r_fidelity < 0.5 or (np.sign(b_deconv) != np.sign(b_umi) and abs(b_deconv - b_umi) > 0.3):
        return "Deconvolution Collinear Leakage"

    # Gate 2: mRNA Mass / Cell Size Bias
    if np.sign(b_umi) != np.sign(b_count) and abs(b_umi - b_count) > 0.3:
        return "mRNA Mass / Cell Size Disparity"

    # Gate 3: Statistical Model (Milo k-NN vs Global Cluster)
    if np.sign(b_count) != np.sign(b_milo) and abs(b_count - b_milo) > 0.5:
        return "Milo Local Graph vs Cluster Discrepancy"

    # Gate 4: Cross-cohort Biological Divergence
    if np.sign(b_count) != np.sign(b_bulk) and abs(b_bulk) > 0.1:
        return "Cross-Cohort Biological Heterogeneity"

    return "Indeterminate / Low Effect Size"


def compute_resolution_trends(
    df: pl.DataFrame,
) -> pl.DataFrame:
    """Compute overall concordance and correlation metrics across resolutions."""
    trends: list[dict[str, object]] = []

    grouped = df.partition_by(["condition", "resolution"], as_dict=True)

    for (cond, res), sub_df in grouped.items():
        n_states = sub_df.height
        if n_states < 3:
            continue

        b_cnt = sub_df["beta_count"].to_numpy()
        b_umi = sub_df["beta_umi"].to_numpy()
        b_dec = sub_df["beta_deconv"].to_numpy()
        b_mil = sub_df["beta_milo"].to_numpy()
        b_blk = sub_df["beta_bulk"].to_numpy()

        # Pairwise correlations with variance guard
        r_dec_umi = float(stats.spearmanr(b_dec, b_umi)[0]) if np.std(b_dec) > 1e-6 and np.std(b_umi) > 1e-6 else 0.0
        r_umi_cnt = float(stats.spearmanr(b_umi, b_cnt)[0]) if np.std(b_umi) > 1e-6 and np.std(b_cnt) > 1e-6 else 0.0
        r_cnt_mil = float(stats.spearmanr(b_cnt, b_mil)[0]) if np.std(b_cnt) > 1e-6 and np.std(b_mil) > 1e-6 else 0.0
        r_dec_blk = float(stats.spearmanr(b_dec, b_blk)[0]) if np.std(b_dec) > 1e-6 and np.std(b_blk) > 1e-6 else 0.0

        # Directional concordance with bulk
        concordant_with_bulk = sum(
            np.sign(d) == np.sign(b) for d, b in zip(b_dec, b_blk, strict=False) if abs(d) > 0.05 and abs(b) > 0.05
        )
        total_eval = sum(1 for d, b in zip(b_dec, b_blk, strict=False) if abs(d) > 0.05 and abs(b) > 0.05)
        conc_rate = concordant_with_bulk / max(total_eval, 1)

        trends.append({
            "condition": cond,
            "resolution": float(res),
            "n_states": n_states,
            "spearman_deconv_vs_umi": float(r_dec_umi),
            "spearman_umi_vs_count": float(r_umi_cnt),
            "spearman_count_vs_milo": float(r_cnt_mil),
            "spearman_deconv_vs_bulk": float(r_dec_blk),
            "concordance_rate_bulk": float(conc_rate),
        })

    return pl.DataFrame(trends).sort(["condition", "resolution"])


def run_pipeline(config: DiagnosticConfig) -> Result[None, str]:
    """Execute Step 4 discrepancy diagnostics."""
    if not config.triangulation_path.exists():
        return Failure(f"Triangulation file not found: {config.triangulation_path}")

    triang_df = pl.read_parquet(config.triangulation_path)
    fidelity_df = pl.read_parquet(config.fidelity_path) if config.fidelity_path.exists() else pl.DataFrame()

    # Merge fidelity if available
    if not fidelity_df.is_empty():
        fid_sub = fidelity_df.select([
            pl.col("condition"),
            pl.col("resolution"),
            pl.col("cell_state"),
            pl.col("pearson_r_umi_vs_deconv").alias("r_fidelity"),
        ])
        merged = triang_df.join(fid_sub, on=["condition", "resolution", "cell_state"], how="left").fill_null(0.7)
    else:
        merged = triang_df.with_columns(pl.lit(0.7).alias("r_fidelity"))

    # Attribute primary cause per state
    attributions: list[str] = []
    for row in merged.iter_rows(named=True):
        cause = attribute_primary_cause(
            row["beta_count"],
            row["beta_umi"],
            row["beta_deconv"],
            row["beta_milo"],
            row["beta_bulk"],
            row["r_fidelity"],
        )
        attributions.append(cause)

    diag_df = merged.with_columns(pl.Series("primary_discrepancy_cause", attributions))

    # Resolution trend summary
    trends_df = compute_resolution_trends(merged)

    config.out_dir.mkdir(parents=True, exist_ok=True)
    diag_df.write_parquet(config.out_dir / "triangulation_diagnostics_summary.parquet")
    trends_df.write_parquet(config.out_dir / "resolution_trend_summary.parquet")

    print("==================================================================")
    print("RESOLUTION TRENDS & TRIANGULATION DIAGNOSTIC SUMMARY")
    print("==================================================================")
    for row in trends_df.filter(pl.col("condition") == "Combined").iter_rows(named=True):
        print(
            f"Res: {row['resolution']:<4.2f} | States: {row['n_states']:>2} | "
            f"Deconv vs UMI ρ: {row['spearman_deconv_vs_umi']:>5.2f} | "
            f"UMI vs Count ρ: {row['spearman_umi_vs_count']:>5.2f} | "
            f"Bulk Concordance: {row['concordance_rate_bulk']*100:>5.1f}%"
        )
    print("==================================================================")

    return Success(None)


def main() -> None:
    config = parse_args()
    match run_pipeline(config):
        case Success(_):
            sys.exit(0)
        case Failure(err):
            print(f"[ERROR] {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
