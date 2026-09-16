#!/usr/bin/env python3
"""
Step 5: Extended Concordance Analysis between Bulk Deconvolution and Single-Cell Milo DA.
Evaluates concordance across:
1. All 9 individual iAtlas ICI cohorts (Hugo, Riaz, Liu, Gide, Rosenberg, Padron, Anders, McDermott, Choueiri)
   and aggregate strata (Melanoma, Pan-Cancer).
2. Sade-Feldman treatment timepoints: Pre-treatment, Post-treatment, and Combined.
Outputs cohort_level_concordance_summary.parquet, timepoint_concordance_summary.parquet,
and concordance_metrics_full.parquet.
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


COHORT_CANCER_MAP: Final[dict[str, str]] = {
    "Hugo-iAtlas": "Melanoma",
    "Riaz-iAtlas": "Melanoma",
    "Liu-iAtlas": "Melanoma",
    "Gide-iAtlas": "Melanoma",
    "Rosenberg-iAtlas": "Bladder",
    "Padron-iAtlas": "Pancreatic",
    "Anders-iAtlas": "Breast",
    "McDermott-iAtlas": "Renal Cell",
    "Choueiri-iAtlas": "Renal Cell",
    "Melanoma": "Melanoma (Combined)",
    "Pan-Cancer": "Pan-Cancer (Combined)",
}


class ConcordanceConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    logistic_path: Path
    milo_dir: Path
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


def compute_pair_concordance(
    df_log_sub: pl.DataFrame,
    df_milo_sub: pl.DataFrame,
    stratum: str,
    condition: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Compute concordance statistics for a single (stratum, condition) pair."""
    # Clean cohort name
    clean_cohort = stratum.replace("Cohort_", "")
    cancer_type = COHORT_CANCER_MAP.get(clean_cohort, "Other")

    # Join on cell_state
    joined = df_log_sub.join(df_milo_sub, on="cell_state", how="inner")
    n_states = joined.height

    if n_states < 3:
        summary = {
            "stratum": stratum,
            "cohort": clean_cohort,
            "cancer_type": cancer_type,
            "condition": condition,
            "n_cell_states": n_states,
            "spearman_rho": np.nan,
            "spearman_pvalue": 1.0,
            "pearson_r": np.nan,
            "pearson_pvalue": 1.0,
            "concordance_percentage": 0.0,
            "binomial_pvalue": 1.0,
            "concordant_states_count": 0,
            "concordant_responders": 0,
            "concordant_non_responders": 0,
            "discordant_count": 0,
        }
        return summary, []

    betas = joined["beta"].to_numpy().astype(np.float64)
    milo_lfcs = joined["milo_mean_logfc"].to_numpy().astype(np.float64)

    # Correlations
    try:
        rho_val, rho_p = stats.spearmanr(betas, milo_lfcs)
        if np.isnan(rho_val):
            rho_val, rho_p = 0.0, 1.0
    except Exception:
        rho_val, rho_p = 0.0, 1.0

    try:
        r_val, r_p = stats.pearsonr(betas, milo_lfcs)
        if np.isnan(r_val):
            r_val, r_p = 0.0, 1.0
    except Exception:
        r_val, r_p = 0.0, 1.0

    # Quadrants
    quadrants = [classify_quadrant(b, m) for b, m in zip(betas, milo_lfcs)]
    n_concordant = sum(1 for q in quadrants if q.startswith("Concordant"))
    n_conc_r = sum(1 for q in quadrants if q == "Concordant Responder")
    n_conc_nr = sum(1 for q in quadrants if q == "Concordant Non-Responder")
    n_disc = sum(1 for q in quadrants if q.startswith("Discordant"))

    pct_agreement = float(n_concordant / n_states) * 100.0
    try:
        binom_p = float(stats.binomtest(k=n_concordant, n=n_states, p=0.5, alternative="greater").pvalue)
    except Exception:
        binom_p = 1.0

    summary = {
        "stratum": stratum,
        "cohort": clean_cohort,
        "cancer_type": cancer_type,
        "condition": condition,
        "n_cell_states": n_states,
        "spearman_rho": float(rho_val),
        "spearman_pvalue": float(rho_p),
        "pearson_r": float(r_val),
        "pearson_pvalue": float(r_p),
        "concordance_percentage": pct_agreement,
        "binomial_pvalue": binom_p,
        "concordant_states_count": n_concordant,
        "concordant_responders": n_conc_r,
        "concordant_non_responders": n_conc_nr,
        "discordant_count": n_disc,
    }

    metric_rows: list[dict[str, object]] = []
    for idx in range(n_states):
        row = joined.row(idx, named=True)
        metric_rows.append(
            {
                "stratum": stratum,
                "cohort": clean_cohort,
                "cancer_type": cancer_type,
                "condition": condition,
                "cell_state": row["cell_state"],
                "beta": float(row["beta"]),
                "or": float(row["or"]),
                "logistic_pvalue": float(row["p_value"]),
                "logistic_fdr": float(row.get("fdr", 1.0)),
                "auc": float(row.get("auc", 0.5)),
                "milo_mean_logfc": float(row["milo_mean_logfc"]),
                "milo_median_logfc": float(row.get("milo_median_logfc", 0.0)),
                "milo_wilcoxon_pval": float(row.get("milo_wilcoxon_pval", 1.0)),
                "quadrant": quadrants[idx],
                "is_concordant": quadrants[idx].startswith("Concordant"),
            }
        )

    return summary, metric_rows


def run_full_concordance(config: ConcordanceConfig) -> Result[Path, str]:
    """Orchestrates comprehensive concordance testing across cohorts and timepoints."""
    if not config.logistic_path.exists():
        return Failure(f"Logistic regression results missing: {config.logistic_path}")

    df_log = pl.read_parquet(config.logistic_path)
    all_strata = sorted(df_log["stratum"].unique().to_list())
    print(f"Loaded logistic results with {len(all_strata)} strata: {all_strata}")

    # Load Milo results for all available conditions
    conditions = ["Combined", "Pre", "Post"]
    milo_condition_dfs: dict[str, pl.DataFrame] = {}

    for cname in conditions:
        p_path = config.milo_dir / f"milopy_cell_state_da_{cname}.parquet"
        if not p_path.exists():
            # Fallback for Combined
            if cname == "Combined":
                p_path = config.milo_dir / "milopy_cell_state_da.parquet"

        if p_path.exists():
            milo_condition_dfs[cname] = pl.read_parquet(p_path)
            print(f"Loaded Milo cell state results for condition '{cname}' ({milo_condition_dfs[cname].height} rows)")

    if not milo_condition_dfs:
        return Failure(f"No Milo cell state result parquets found in {config.milo_dir}")

    summaries: list[dict[str, object]] = []
    all_metric_rows: list[dict[str, object]] = []

    for strat in all_strata:
        df_log_sub = df_log.filter(pl.col("stratum") == strat)

        for cname, df_milo_sub in milo_condition_dfs.items():
            sum_dict, m_rows = compute_pair_concordance(
                df_log_sub, df_milo_sub, strat, cname
            )
            summaries.append(sum_dict)
            all_metric_rows.extend(m_rows)

    df_summary_master = pl.DataFrame(summaries)
    df_metrics_master = pl.DataFrame(all_metric_rows)

    config.out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Full metrics table (all strata x conditions x cell states)
    out_metrics_full = config.out_dir / "concordance_metrics_full.parquet"
    df_metrics_master.write_parquet(out_metrics_full)
    print(f"Saved full concordance metrics ({df_metrics_master.height} rows) to: {out_metrics_full}")

    # 2. Master summary table
    out_summary_master = config.out_dir / "concordance_summary_all.parquet"
    df_summary_master.write_parquet(out_summary_master)

    # 3. Cohort-level summary table (individual cohorts vs Combined Milo DA)
    df_cohorts = df_summary_master.filter(
        pl.col("stratum").str.starts_with("Cohort_") & (pl.col("condition") == "Combined")
    )
    out_cohorts = config.out_dir / "cohort_level_concordance_summary.parquet"
    df_cohorts.write_parquet(out_cohorts)
    print(f"Saved cohort-level concordance summary ({df_cohorts.height} cohorts) to: {out_cohorts}")

    # 4. Timepoint summary table (Pre vs Post vs Combined for Melanoma and Pan-Cancer)
    df_timepoints = df_summary_master.filter(
        pl.col("stratum").is_in(["Melanoma", "Pan-Cancer"])
    )
    out_timepoints = config.out_dir / "timepoint_concordance_summary.parquet"
    df_timepoints.write_parquet(out_timepoints)
    print(f"Saved timepoint concordance summary to: {out_timepoints}")

    # 5. Default backward-compatible tables (Melanoma x Combined)
    df_default_metrics = df_metrics_master.filter(
        (pl.col("stratum") == "Melanoma") & (pl.col("condition") == "Combined")
    )
    df_default_metrics.write_parquet(config.out_dir / "concordance_metrics.parquet")

    df_default_summary = df_summary_master.filter(
        (pl.col("stratum") == "Melanoma") & (pl.col("condition") == "Combined")
    )
    df_default_summary.write_parquet(config.out_dir / "concordance_summary.parquet")

    # Print summary table of individual cohorts
    print("\n=========================================================================")
    print("INDIVIDUAL COHORT CONCORDANCE ANALYSIS (Bulk Deconv vs. Milo Combined)")
    print("=========================================================================")
    for row in df_cohorts.iter_rows(named=True):
        print(
            f"Cohort: {row['cohort']:<20} | Cancer: {row['cancer_type']:<12} | "
            f"Spearman rho = {row['spearman_rho']:>6.3f} (p={row['spearman_pvalue']:.3e}) | "
            f"Concordance = {row['concordance_percentage']:>5.1f}% ({row['concordant_states_count']}/{row['n_cell_states']})"
        )
    print("=========================================================================")

    return Success(out_metrics_full)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 5: Extended cross-modality concordance across cohorts and timepoints."
    )
    parser.add_argument(
        "--logistic-results",
        type=str,
        default="output/sade_feldman_deconv_validation/logistic_regression_results.parquet",
        help="Path to logistic_regression_results.parquet",
    )
    parser.add_argument(
        "--milo-dir",
        type=str,
        default="output/sade_feldman_deconv_validation",
        help="Directory containing milopy_cell_state_da parquets",
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
        milo_dir=Path(args.milo_dir),
        out_dir=Path(args.out_dir),
    )

    match run_full_concordance(config):
        case Success(out_path):
            print(f"Step 5 completed successfully: {out_path}")
            sys.exit(0)
        case Failure(err):
            print(f"Step 5 failed with error: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
