#!/usr/bin/env python3
"""Step 5: Concordance Analysis and Unified Benchmark Metric Compilation.
Merges deconvolution response predictive metrics (AUCs, betas), single-cell Milo DA logFCs,
and reference condition numbers across all single-cell and combined reference strategies.
Outputs the master benchmark table: multi_resolution_benchmark_summary.parquet.

Follows functional programming standards: pure functions, immutable data structures,
and explicit error handling using Returns Result monads.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import stats  # type: ignore[import-untyped]


class ConcordanceConfig(BaseModel):
    """Immutable configuration for benchmark compilation."""

    model_config = ConfigDict(frozen=True)

    ref_dir: Path
    logistic_dir: Path
    milo_dir: Path
    out_dir: Path


def compile_benchmark_summary(config: ConcordanceConfig) -> Result[Path, str]:
    """Compiles multi-resolution benchmark summary table across all references."""
    config.out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load logistic regression results
    log_file = config.logistic_dir / "logistic_regression_results_all_resolutions.parquet"
    if not log_file.exists():
        return Failure(f"Logistic regression results not found at: {log_file}")
    df_log_all = pl.read_parquet(log_file)

    # 2. Gather all reference resolution metrics
    metric_files = sorted(config.ref_dir.glob("reference_resolution_metrics_*.parquet"))
    metrics_map: dict[str, dict[float, dict[str, Any]]] = {}

    for mf in metric_files:
        df_m = pl.read_parquet(mf)
        ref_stem = mf.stem.replace("reference_resolution_metrics_", "")
        for row in df_m.iter_rows(named=True):
            r_type = str(row.get("reference_type", ""))
            r_val = float(row.get("resolution", 0.5))
            if ref_stem not in metrics_map:
                metrics_map[ref_stem] = {}
            metrics_map[ref_stem][r_val] = row
            if r_type not in metrics_map:
                metrics_map[r_type] = {}
            metrics_map[r_type][r_val] = row

    # 3. Load Milo DA if available
    milo_file = config.milo_dir / "milopy_cell_state_da_all_resolutions.parquet"
    df_milo = pl.read_parquet(milo_file) if milo_file.exists() else None

    # Group logistic results by (reference_type, resolution)
    unique_configs = (
        df_log_all.select(["reference_type", "reference_category", "resolution"])
        .unique()
        .sort(["reference_type", "resolution"])
    )

    benchmark_records: list[dict[str, object]] = []

    for row in unique_configs.iter_rows(named=True):
        ref_id = str(row["reference_type"])
        ref_cat = str(row["reference_category"])
        res_val = float(row["resolution"])

        sub_df = df_log_all.filter(
            (pl.col("reference_type") == ref_id) & (pl.col("resolution") == res_val)
        )

        sub_mel = sub_df.filter(pl.col("stratum") == "Melanoma")
        sub_pan = sub_df.filter(pl.col("stratum") == "Pan-Cancer")
        cohort_rows = sub_df.filter(pl.col("stratum").str.starts_with("Cohort_"))

        mel_multi_auc = (
            float(sub_mel["multivariate_auc"][0])
            if sub_mel.height > 0 and "multivariate_auc" in sub_mel.columns
            else 0.5
        )
        pan_multi_auc = (
            float(sub_pan["multivariate_auc"][0])
            if sub_pan.height > 0 and "multivariate_auc" in sub_pan.columns
            else 0.5
        )
        mean_cohort_auc = (
            float(cohort_rows["multivariate_auc"].mean())
            if cohort_rows.height > 0 and "multivariate_auc" in cohort_rows.columns
            else 0.5
        )
        mean_univ_auc = float(sub_mel["auc"].mean()) if sub_mel.height > 0 else 0.5

        # Match metadata from metrics
        matched_meta = metrics_map.get(ref_id, {}).get(res_val, {})
        if not matched_meta:
            ref_clean = ref_id.lower().replace("-", "_")
            for k, res_dict in metrics_map.items():
                if k.lower().replace("-", "_") == ref_clean:
                    matched_meta = res_dict.get(res_val, res_dict.get(0.5, {}))
                    break

        # Accurately resolve reference category and canonical name
        ref_cat = str(matched_meta.get("reference_category", ""))
        lower = ref_id.lower()
        if not ref_cat or ref_cat == "Single Dataset":
            if "_sub" in lower or "subsample" in lower or "%" in lower:
                ref_cat = "Cell-Subsample"
            elif lower.startswith("random"):
                ref_cat = "Random-Combined"
            elif any(k in lower for k in ("duo", "trio", "quad", "quartet", "combined", "atlas")):
                ref_cat = "Criteria-Combined"
            else:
                ref_cat = "Single Dataset"

        cond_num = float(matched_meta.get("condition_number", np.nan))
        n_clusters = int(matched_meta.get("n_clusters", sub_mel.height if sub_mel.height > 0 else 10))
        n_sig_genes = int(matched_meta.get("n_signature_genes", 0))
        n_datasets = int(matched_meta.get("n_datasets", 1 if ref_cat == "Single Dataset" else 2))
        n_cells = int(matched_meta.get("n_cells", 20000))
        sub_frac = float(matched_meta.get("subsample_fraction", 1.0))
        canonical_name = str(matched_meta.get("reference_type", ref_id))

        # Concordance with Milo DA
        rho_val = float("nan")
        pct_conc = float("nan")
        if df_milo is not None:
            milo_res_sub = df_milo.filter(pl.col("resolution") == res_val)
            if milo_res_sub.height > 0:
                joined = sub_mel.join(milo_res_sub, on="cell_state", how="inner")
                if joined.height >= 3:
                    try:
                        r_val, _ = stats.spearmanr(
                            joined["beta"].to_numpy(), joined["milo_mean_logfc"].to_numpy()
                        )
                        rho_val = float(r_val) if not np.isnan(r_val) else 0.0
                    except Exception:
                        rho_val = 0.0
                    n_conc = sum(
                        1
                        for b, m in zip(joined["beta"], joined["milo_mean_logfc"])
                        if (b > 0 and m > 0) or (b < 0 and m < 0)
                    )
                    pct_conc = float(n_conc / joined.height) * 100.0

        cohort_dict: dict[str, float] = {}
        for c_row in cohort_rows.iter_rows(named=True):
            c_name = str(c_row["stratum"]).replace("Cohort_", "").replace("-iAtlas", "")
            c_auc = float(c_row.get("multivariate_auc", 0.5))
            cohort_dict[f"auc_{c_name.lower()}"] = c_auc

        worst_cohort_auc = float(cohort_rows["multivariate_auc"].min()) if cohort_rows.height > 0 else 0.5

        benchmark_records.append(
            {
                "reference_type": canonical_name,
                "reference_category": ref_cat,
                "resolution": res_val,
                "n_clusters": n_clusters,
                "n_signature_genes": n_sig_genes,
                "condition_number": cond_num,
                "spearman_rho": rho_val,
                "concordance_percentage": pct_conc,
                "melanoma_multivariate_auc": mel_multi_auc,
                "pancancer_multivariate_auc": pan_multi_auc,
                "mean_cohort_multivariate_auc": mean_cohort_auc,
                "worst_cohort_multivariate_auc": worst_cohort_auc,
                "mean_univariate_auc": mean_univ_auc,
                "n_datasets": n_datasets,
                "n_cells": n_cells,
                "subsample_fraction": sub_frac,
                **cohort_dict,
            }
        )

    if not benchmark_records:
        return Failure("No benchmark records compiled.")

    df_bench = pl.DataFrame(benchmark_records)

    # Compute Pareto optimality across mean_cohort_multivariate_auc (max) and condition_number (min)
    pareto_flags: list[bool] = []
    aucs = df_bench["mean_cohort_multivariate_auc"].to_numpy()
    conds = df_bench["condition_number"].to_numpy()
    n_pts = len(df_bench)

    for i in range(n_pts):
        auc_i = aucs[i]
        cond_i = conds[i]
        if np.isnan(auc_i) or np.isnan(cond_i):
            pareto_flags.append(False)
            continue
        is_dom = False
        for j in range(n_pts):
            if i == j:
                continue
            auc_j = aucs[j]
            cond_j = conds[j]
            if np.isnan(auc_j) or np.isnan(cond_j):
                continue
            if auc_j >= auc_i and cond_j <= cond_i and (auc_j > auc_i or cond_j < cond_i):
                is_dom = True
                break
        pareto_flags.append(not is_dom)

    df_bench = df_bench.with_columns(pl.Series("is_pareto_optimal", pareto_flags))
    out_master = config.out_dir / "multi_resolution_benchmark_summary.parquet"
    df_bench.write_parquet(out_master)
    print(f"\nSaved master multi-resolution benchmark summary ({df_bench.height} configurations) -> {out_master}")

    return Success(out_master)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 5: Compile concordance and predictive capacity benchmark summary."
    )
    parser.add_argument(
        "--ref-dir",
        type=str,
        default="output/sc9_milo_deconv_validation/01_references",
        help="Directory containing reference_resolution_metrics parquets",
    )
    parser.add_argument(
        "--logistic-dir",
        type=str,
        default="output/sc9_milo_deconv_validation/03_logistic_regression",
        help="Directory containing logistic regression parquet",
    )
    parser.add_argument(
        "--milo-dir",
        type=str,
        default="output/sc9_milo_deconv_validation/04_milopy",
        help="Directory containing Milo DA parquet",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="output/sc9_milo_deconv_validation/05_concordance",
        help="Directory to save compiled benchmark table",
    )
    args = parser.parse_args()

    config = ConcordanceConfig(
        ref_dir=Path(args.ref_dir),
        logistic_dir=Path(args.logistic_dir),
        milo_dir=Path(args.milo_dir),
        out_dir=Path(args.out_dir),
    )

    match compile_benchmark_summary(config):
        case Success(p):
            print(f"Benchmark compilation succeeded: {p}")
            sys.exit(0)
        case Failure(err):
            print(f"Benchmark compilation failed: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
