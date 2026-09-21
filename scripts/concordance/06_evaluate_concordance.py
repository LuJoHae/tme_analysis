#!/usr/bin/env python3
"""
Step 6: Concordance Evaluation and Diagnostic Discrepancy Attribution.
Quantifies rank concordance, directional quadrant consistency, and statistical tests
between single-cell and bulk response associations, attributing discrepancies to
tumor purity, mRNA content disparity, or collinear deconvolution leakage.
Strict functional Python with returns, Pydantic, and Polars.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

import numpy as np
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
import polars as pl
from scipy import stats  # type: ignore


class EvalConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    sc_effects_path: Path
    bulk_effects_path: Path
    insilico_path: Path
    collinearity_path: Path
    out_dir: Path


def parse_args() -> EvalConfig:
    parser = argparse.ArgumentParser(
        description="Evaluate cross-cohort concordance and diagnose sources of discrepancy."
    )
    parser.add_argument("--sc-effects", type=Path, required=True, help="Single-cell effects parquet")
    parser.add_argument("--bulk-effects", type=Path, required=True, help="Bulk response effects parquet")
    parser.add_argument("--insilico", type=Path, required=True, help="In silico recovery comparison parquet")
    parser.add_argument("--collinearity", type=Path, required=True, help="Reference collinearity parquet")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    args = parser.parse_args()
    return EvalConfig(
        sc_effects_path=args.sc_effects,
        bulk_effects_path=args.bulk_effects,
        insilico_path=args.insilico,
        collinearity_path=args.collinearity,
        out_dir=args.out_dir,
    )


def classify_quadrant(beta_sc: float, beta_bulk: float) -> str:
    """Classify directional agreement between single-cell and bulk effect sizes."""
    if beta_sc > 0 and beta_bulk > 0:
        return "Concordant Responder"
    elif beta_sc < 0 and beta_bulk < 0:
        return "Concordant Non-Responder"
    elif beta_sc > 0 and beta_bulk < 0:
        return "Discordant (SC+, Bulk-)"
    elif beta_sc < 0 and beta_bulk > 0:
        return "Discordant (SC-, Bulk+)"
    else:
        return "Neutral"


def compute_concordance_for_fraction_type(
    sc_df: pl.DataFrame,
    bulk_df: pl.DataFrame,
    frac_type: str,
) -> tuple[dict[str, object], pl.DataFrame]:
    """Calculate concordance metrics and state-level quadrant labels for a fraction type."""
    sub_bulk = bulk_df.filter(pl.col("fraction_type") == frac_type)
    joined = sc_df.join(sub_bulk, on="cell_state", how="inner")

    n_states = joined.height
    if n_states < 3:
        summary: dict[str, object] = {
            "fraction_type": frac_type,
            "n_states": n_states,
            "spearman_rho": float("nan"),
            "spearman_pval": 1.0,
            "kendall_tau": float("nan"),
            "kendall_pval": 1.0,
            "pearson_r": float("nan"),
            "concordance_rate": 0.0,
            "binomial_pval": 1.0,
            "n_concordant": 0,
            "n_discordant": 0,
        }
        return summary, joined

    sc_betas = joined["beta_sc"].to_numpy().astype(np.float64)
    bulk_betas = joined["beta_bulk"].to_numpy().astype(np.float64)

    # Rank correlation
    rho, rho_p = stats.spearmanr(sc_betas, bulk_betas)
    tau, tau_p = stats.kendalltau(sc_betas, bulk_betas)
    r_val, _ = stats.pearsonr(sc_betas, bulk_betas)

    # Quadrants
    quadrants = [classify_quadrant(b_sc, b_bk) for b_sc, b_bk in zip(sc_betas, bulk_betas, strict=False)]
    concordant_mask = [q.startswith("Concordant") for q in quadrants]
    n_concordant = sum(concordant_mask)
    n_discordant = n_states - n_concordant

    # Binomial test against 50% random chance
    binom_res = stats.binomtest(n_concordant, n_states, p=0.5, alternative="greater")

    detailed_df = joined.with_columns([
        pl.Series("quadrant", quadrants),
        pl.Series("is_concordant", concordant_mask),
    ])

    summary_dict: dict[str, object] = {
        "fraction_type": frac_type,
        "n_states": n_states,
        "spearman_rho": float(rho),
        "spearman_pval": float(rho_p),
        "kendall_tau": float(tau),
        "kendall_pval": float(tau_p),
        "pearson_r": float(r_val),
        "concordance_rate": float(n_concordant / n_states),
        "binomial_pval": float(binom_res.pvalue),
        "n_concordant": int(n_concordant),
        "n_discordant": int(n_discordant),
    }

    return summary_dict, detailed_df


def diagnose_discordant_states(
    detailed_df: pl.DataFrame,
    insilico_df: pl.DataFrame,
    collinearity_df: pl.DataFrame,
) -> pl.DataFrame:
    """Attribute reasons for discordance: deconvolution identifiability vs collinearity."""
    discordant_states = (
        detailed_df.filter(pl.col("fraction_type") == "normalized")
        .filter(~pl.col("is_concordant"))
        ["cell_state"]
        .to_list()
    )

    insilico_map = {
        row["cell_state"]: (row["pearson_r_recovery"], row["is_identifiable"])
        for row in insilico_df.to_dicts()
    }

    diagnostics: list[dict[str, object]] = []

    for state in discordant_states:
        # Check in silico identifiability
        r_rec, is_id = insilico_map.get(state, (0.0, False))

        # Check collinearity with other states
        collin_partners = collinearity_df.filter(
            (pl.col("state_a") == state) | (pl.col("state_b") == state)
        ).filter(pl.col("is_collinear"))

        has_collinear = collin_partners.height > 0
        partners_str = (
            ",".join(
                [
                    row["state_b"] if row["state_a"] == state else row["state_a"]
                    for row in collin_partners.to_dicts()
                ]
            )
            if has_collinear
            else "None"
        )

        primary_cause = (
            "Collinear Deconvolution Leakage"
            if has_collinear
            else (
                "Unidentifiable scRNA-seq State"
                if not is_id
                else "Cohort Phenotypic Divergence / Expression Shift"
            )
        )

        diagnostics.append({
            "cell_state": state,
            "primary_discrepancy_cause": primary_cause,
            "insilico_recovery_r": float(r_rec),
            "is_insilico_identifiable": bool(is_id),
            "collinear_partners": partners_str,
        })

    return pl.DataFrame(diagnostics) if diagnostics else pl.DataFrame(schema={
        "cell_state": pl.String,
        "primary_discrepancy_cause": pl.String,
        "insilico_recovery_r": pl.Float64,
        "is_insilico_identifiable": pl.Boolean,
        "collinear_partners": pl.String,
    })


def run_pipeline(config: EvalConfig) -> Result[None, str]:
    """Execute Step 6 concordance evaluation and diagnostics."""
    for p in [config.sc_effects_path, config.bulk_effects_path, config.insilico_path, config.collinearity_path]:
        if not p.exists():
            return Failure(f"Input file does not exist: {p}")

    sc_df = pl.read_parquet(config.sc_effects_path)
    bulk_df = pl.read_parquet(config.bulk_effects_path)
    insilico_df = pl.read_parquet(config.insilico_path)
    collinearity_df = pl.read_parquet(config.collinearity_path)

    summaries: list[dict[str, object]] = []
    detailed_tables: list[pl.DataFrame] = []

    for f_type in ["raw", "normalized", "mrna_scaled"]:
        if f_type in bulk_df["fraction_type"].unique().to_list():
            summary, det_df = compute_concordance_for_fraction_type(sc_df, bulk_df, f_type)
            summaries.append(summary)
            detailed_tables.append(det_df)

    summary_df = pl.DataFrame(summaries)
    all_details_df = pl.concat(detailed_tables)

    # Diagnostics on normalized fractions
    diagnostics_df = diagnose_discordant_states(all_details_df, insilico_df, collinearity_df)

    config.out_dir.mkdir(parents=True, exist_ok=True)
    summary_df.write_parquet(config.out_dir / "concordance_metrics_summary.parquet")
    all_details_df.write_parquet(config.out_dir / "concordance_state_details.parquet")
    diagnostics_df.write_parquet(config.out_dir / "concordance_diagnostics.parquet")

    # Print summary
    print("==================================================================")
    print("CROSS-COHORT CONCORDANCE ANALYSIS SUMMARY")
    print("==================================================================")
    for row in summary_df.iter_rows(named=True):
        print(
            f"Fraction: {row['fraction_type']:<12} | "
            f"Spearman rho: {row['spearman_rho']:>6.3f} (p={row['spearman_pval']:.3e}) | "
            f"Concordance Rate: {row['concordance_rate']*100:>5.1f}% | "
            f"Binomial p: {row['binomial_pval']:.3e}"
        )
    print("==================================================================")

    return Success(None)


def main() -> None:
    config = parse_args()
    match run_pipeline(config):
        case Success(_):
            print(f"[SUCCESS] Concordance evaluation written to {config.out_dir}")
            sys.exit(0)
        case Failure(err):
            print(f"[ERROR] {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
