#!/usr/bin/env python3
"""
Step 5: Concordance Analysis between Bulk Deconvolution Logistic Regression and Single-Cell Milo DA.
Evaluates whether cell states correlating with immunotherapy response in bulk deconvolution
are concordant with single-cell differential abundance in Sade-Feldman dataset.
Outputs concordance_metrics.parquet and concordance_summary.parquet.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import stats  # type: ignore


class ConcordanceConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    logistic_path: Path
    milo_path: Path
    stratum: str = "Melanoma"
    out_dir: Path


def classify_quadrant(beta: float, milo_lfc: float) -> str:
    """Classify concordant vs discordant directionality."""
    if beta > 0 and milo_lfc > 0:
        return "Concordant Responder"
    elif beta < 0 and milo_lfc < 0:
        return "Concordant Non-Responder"
    elif beta > 0 and milo_lfc < 0:
        return "Discordant (Bulk+, SC-)"
    elif beta < 0 and milo_lfc > 0:
        return "Discordant (Bulk-, SC+)"
    else:
        return "Neutral"


def run_concordance_analysis(config: ConcordanceConfig) -> Result[Path, str]:
    """Pure analysis joining bulk and single-cell response associations."""
    if not config.logistic_path.exists():
        return Failure(f"Logistic results not found: {config.logistic_path}")
    if not config.milo_path.exists():
        return Failure(f"Milo results not found: {config.milo_path}")

    df_log = pl.read_parquet(config.logistic_path)
    df_milo = pl.read_parquet(config.milo_path)

    # Filter to target stratum
    df_log_sub = df_log.filter(pl.col("stratum") == config.stratum)
    if df_log_sub.height == 0:
        available = df_log["stratum"].unique().to_list()
        return Failure(
            f"Stratum '{config.stratum}' not found in logistic results. Available: {available}"
        )

    # Join on cell_state
    df_merged = df_log_sub.join(df_milo, on="cell_state", how="inner")
    if df_merged.height == 0:
        return Failure("No overlapping cell states between logistic results and Milo results.")

    print(f"Comparing {df_merged.height} cell states for stratum '{config.stratum}'...")

    betas = df_merged["beta"].to_numpy().astype(np.float64)
    milo_lfcs = df_merged["milo_mean_logfc"].to_numpy().astype(np.float64)
    milo_medians = df_merged["milo_median_logfc"].to_numpy().astype(np.float64)

    # Correlation statistics
    spearman_rho, spearman_p = stats.spearmanr(betas, milo_lfcs)
    pearson_r, pearson_p = stats.pearsonr(betas, milo_lfcs)
    spearman_rho_med, spearman_p_med = stats.spearmanr(betas, milo_medians)

    # Quadrant agreement
    quadrants = [classify_quadrant(b, m) for b, m in zip(betas, milo_lfcs)]
    n_concordant = sum(1 for q in quadrants if q.startswith("Concordant"))
    pct_agreement = float(n_concordant / len(quadrants)) * 100.0

    # Binomial test for direction concordance vs 50% chance
    binom_res = stats.binomtest(k=n_concordant, n=len(quadrants), p=0.5, alternative="greater")
    binom_p = float(binom_res.pvalue)

    # Add quadrant labels and signed agreement to dataframe
    df_enhanced = df_merged.with_columns(
        pl.Series("quadrant", quadrants),
        pl.Series(
            "is_concordant",
            [q.startswith("Concordant") for q in quadrants],
        ),
        (pl.col("beta") * pl.col("milo_mean_logfc") > 0).alias("sign_agreement"),
    )

    config.out_dir.mkdir(parents=True, exist_ok=True)
    out_table = config.out_dir / "concordance_metrics.parquet"
    df_enhanced.write_parquet(out_table)
    print(f"Saved concordance metrics table to: {out_table}")

    # Summary table
    df_summary = pl.DataFrame(
        {
            "stratum": [config.stratum],
            "n_cell_states": [df_merged.height],
            "spearman_rho": [float(spearman_rho)],
            "spearman_pvalue": [float(spearman_p)],
            "pearson_r": [float(pearson_r)],
            "pearson_pvalue": [float(pearson_p)],
            "spearman_rho_median": [float(spearman_rho_med)],
            "spearman_pvalue_median": [float(spearman_p_med)],
            "concordant_states_count": [n_concordant],
            "concordance_percentage": [pct_agreement],
            "binomial_pvalue": [binom_p],
        }
    )
    out_summary = config.out_dir / "concordance_summary.parquet"
    df_summary.write_parquet(out_summary)
    print(f"Saved concordance summary table to: {out_summary}")

    print("\n=======================================================")
    print("CONCORDANCE SUMMARY (Bulk Deconv vs. Single-Cell Milo)")
    print("=======================================================")
    print(f"Stratum:                  {config.stratum}")
    print(f"Cell states evaluated:    {df_merged.height}")
    print(f"Spearman Rank Correlation: rho = {spearman_rho:.3f} (p = {spearman_p:.4e})")
    print(f"Pearson Correlation:       r   = {pearson_r:.3f} (p = {pearson_p:.4e})")
    print(f"Directional Concordance:   {pct_agreement:.1f}% ({n_concordant}/{len(quadrants)} states, Binomial p = {binom_p:.4e})")
    print("=======================================================")

    return Success(out_table)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 5: Concordance analysis between bulk deconvolution and Milo DA."
    )
    parser.add_argument(
        "--logistic-results",
        type=str,
        default="output/sade_feldman_deconv_validation/logistic_regression_results.parquet",
        help="Path to logistic_regression_results.parquet from Step 3",
    )
    parser.add_argument(
        "--milo-results",
        type=str,
        default="output/sade_feldman_deconv_validation/milopy_cell_state_da.parquet",
        help="Path to milopy_cell_state_da.parquet from Step 4",
    )
    parser.add_argument(
        "--stratum",
        type=str,
        default="Melanoma",
        help="Stratum to evaluate (Melanoma or Pan-Cancer)",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="output/sade_feldman_deconv_validation",
        help="Directory to save output parquets",
    )
    args = parser.parse_args()

    config = ConcordanceConfig(
        logistic_path=Path(args.logistic_results),
        milo_path=Path(args.milo_results),
        stratum=args.stratum,
        out_dir=Path(args.out_dir),
    )

    match run_concordance_analysis(config):
        case Success(out_path):
            print(f"Step 5 completed successfully: {out_path}")
            sys.exit(0)
        case Failure(err):
            print(f"Step 5 failed with error: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
