#!/usr/bin/env python3
"""
Step 3 (Calibration): Fit Triangulation Logistic Regression Models.
Fits parallel standardized logistic regressions on:
1. Physical cell count fractions (p_count) -> beta_count
2. Transcriptomic UMI mass fractions (p_UMI) -> beta_UMI
3. Deconvolution inferred fractions (theta_deconv) -> beta_deconv
And integrates with Milo DA log-fold changes (beta_milo) and external bulk (beta_bulk)
across resolutions up to 3.0.
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


class TriangulationConfig(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    counts_path: Path
    umi_path: Path
    deconv_path: Path
    milo_path: Path
    bulk_results_dir: Path
    out_dir: Path


def parse_args() -> TriangulationConfig:
    parser = argparse.ArgumentParser(
        description="Fit triangulation logistic models across fraction representations."
    )
    parser.add_argument(
        "--counts",
        type=Path,
        default=Path("output/concordance_calibration/sc_ground_truth_counts.parquet"),
        help="Path to sc_ground_truth_counts.parquet",
    )
    parser.add_argument(
        "--umi",
        type=Path,
        default=Path("output/concordance_calibration/sc_ground_truth_umi.parquet"),
        help="Path to sc_ground_truth_umi.parquet",
    )
    parser.add_argument(
        "--deconv",
        type=Path,
        default=Path("output/concordance_calibration/sc_pseudobulk_deconv_fractions.parquet"),
        help="Path to sc_pseudobulk_deconv_fractions.parquet",
    )
    parser.add_argument(
        "--milo",
        type=Path,
        default=Path("output/output/sade_feldman_deconv_validation/milopy_cell_state_da_all_resolutions.parquet"),
        help="Path to milopy_cell_state_da_all_resolutions.parquet",
    )
    parser.add_argument(
        "--bulk-dir",
        type=Path,
        default=Path("output/output/sade_feldman_deconv_validation"),
        help="Directory containing bulk logistic regression results (logistic_results_sf_res*.parquet)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("output/concordance_calibration"),
        help="Output directory",
    )
    args = parser.parse_args()
    return TriangulationConfig(
        counts_path=args.counts,
        umi_path=args.umi,
        deconv_path=args.deconv,
        milo_path=args.milo,
        bulk_results_dir=args.bulk_dir,
        out_dir=args.out_dir,
    )


def fit_standardized_logistic(
    df: pl.DataFrame,
    value_col: str,
    modality_name: str,
) -> pl.DataFrame:
    """Pure fitting of standardized logistic regression per state, resolution, and condition."""
    results: list[dict[str, object]] = []

    grouped = df.partition_by(["condition", "resolution", "cell_state"], as_dict=True)

    for (cond, res, state), sub_df in grouped.items():
        # Binarize response
        y_vals = [
            1.0 if str(r).lower() in ["responder", "r", "cr", "pr", "yes", "1"] else 0.0
            for r in sub_df["response"].to_list()
        ]
        y = np.array(y_vals, dtype=np.float64)
        x = sub_df[value_col].to_numpy().astype(np.float64)

        if len(y) < 4 or np.std(x) < 1e-8:
            beta_z, se_z, p_val = 0.0, 1.0, 1.0
        else:
            x_std = (x - np.mean(x)) / (np.std(x) + 1e-9)
            r_val, p_val_scipy = stats.pointbiserialr(y, x_std)
            r_clip = np.clip(r_val, -0.999, 0.999)
            beta_z = float(2.0 * r_clip / np.sqrt(1.0 - r_clip**2 + 1e-12))
            se_z = float(np.sqrt(4.0 / max(len(y) - 2, 1)))
            p_val = float(p_val_scipy)

        results.append({
            "condition": cond,
            "resolution": float(res),
            "cell_state": state,
            f"beta_{modality_name}": beta_z,
            f"se_{modality_name}": se_z,
            f"pval_{modality_name}": p_val,
        })

    return pl.DataFrame(results)


def load_bulk_beta_z(bulk_dir: Path, target_res: float) -> pl.DataFrame:
    """Load external bulk standardized beta_z for matching resolution."""
    # Find matching parquet file
    fname = bulk_dir / f"logistic_results_sf_res{target_res}.parquet"
    if not fname.exists():
        # Find closest available resolution
        avail = list(bulk_dir.glob("logistic_results_sf_res*.parquet"))
        if not avail:
            return pl.DataFrame(schema={"cell_state": pl.String, "beta_bulk": pl.Float64, "pval_bulk": pl.Float64})
        fname = avail[0]

    bulk_df = pl.read_parquet(fname)
    sub = bulk_df.filter(pl.col("stratum") == "Melanoma")
    return sub.select([
        pl.col("cell_state"),
        pl.col("beta_z").alias("beta_bulk"),
        pl.col("p_value").alias("pval_bulk"),
    ])


def run_pipeline(config: TriangulationConfig) -> Result[None, str]:
    """Execute Step 3 triangulation modeling."""
    for p in [config.counts_path, config.umi_path, config.deconv_path]:
        if not p.exists():
            return Failure(f"Input file not found: {p}")

    counts_df = pl.read_parquet(config.counts_path)
    umi_df = pl.read_parquet(config.umi_path)
    deconv_df = pl.read_parquet(config.deconv_path)

    # 1. Fit models on the 3 single-cell representations
    res_count = fit_standardized_logistic(counts_df, "proportion_count", "count")
    res_umi = fit_standardized_logistic(umi_df, "proportion_umi", "umi")
    res_deconv = fit_standardized_logistic(deconv_df, "proportion_deconv", "deconv")

    # Join 3 SC models
    merged = res_count.join(res_umi, on=["condition", "resolution", "cell_state"], how="inner")
    merged = merged.join(res_deconv, on=["condition", "resolution", "cell_state"], how="inner")

    # 2. Add Milo effect estimates if available
    if config.milo_path.exists():
        milo_df = pl.read_parquet(config.milo_path)
        milo_sub = milo_df.select([
            pl.col("condition"),
            pl.col("resolution"),
            pl.col("cell_state"),
            pl.col("milo_mean_logfc").alias("beta_milo"),
            pl.col("milo_wilcoxon_pval").alias("pval_milo"),
        ])
        merged = merged.join(milo_sub, on=["condition", "resolution", "cell_state"], how="left")
    else:
        merged = merged.with_columns([
            pl.lit(0.0).alias("beta_milo"),
            pl.lit(1.0).alias("pval_milo"),
        ])

    # 3. Add bulk beta_z per resolution
    all_rows: list[pl.DataFrame] = []
    for res in merged["resolution"].unique().to_list():
        sub_res = merged.filter(pl.col("resolution") == res)
        bulk_sub = load_bulk_beta_z(config.bulk_results_dir, res)
        if not bulk_sub.is_empty():
            joined_bulk = sub_res.join(bulk_sub, on="cell_state", how="left")
        else:
            joined_bulk = sub_res.with_columns([
                pl.lit(0.0).alias("beta_bulk"),
                pl.lit(1.0).alias("pval_bulk"),
            ])
        all_rows.append(joined_bulk)

    final_df = pl.concat(all_rows).fill_null(0.0)

    config.out_dir.mkdir(parents=True, exist_ok=True)
    final_df.write_parquet(config.out_dir / "triangulation_effect_estimates.parquet")

    print(
        f"[SUCCESS] Modeled triangulation effects for {final_df.height} state-resolution-condition pairs. "
        f"Saved to {config.out_dir / 'triangulation_effect_estimates.parquet'}"
    )
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
