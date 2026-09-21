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
import re
import sys
from pathlib import Path
from typing import Final

import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import stats  # type: ignore
from sklearn.metrics import cohen_kappa_score  # type: ignore
import statsmodels.api as sm  # type: ignore


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

    # Cosine Similarity
    norm_b = float(np.linalg.norm(betas))
    norm_m = float(np.linalg.norm(milo_lfcs))
    cosine_sim = float(np.dot(betas, milo_lfcs) / (norm_b * norm_m)) if (norm_b > 0 and norm_m > 0) else 0.0

    # Cohen's Kappa on directionality signs
    signs_b = np.sign(betas)
    signs_m = np.sign(milo_lfcs)
    try:
        cohen_kappa = float(cohen_kappa_score(signs_b, signs_m))
        if np.isnan(cohen_kappa):
            cohen_kappa = 0.0
    except Exception:
        cohen_kappa = 0.0

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

    # Permutation Test (B = 2,000 permutations of beta labels holding milo fixed)
    rng = np.random.default_rng(42)
    n_perm = 2000
    null_rhos = np.zeros(n_perm, dtype=np.float64)
    null_pearsons = np.zeros(n_perm, dtype=np.float64)
    for i in range(n_perm):
        perm_b = rng.permutation(betas)
        try:
            r_s, _ = stats.spearmanr(perm_b, milo_lfcs)
            null_rhos[i] = float(r_s) if not np.isnan(r_s) else 0.0
        except Exception:
            null_rhos[i] = 0.0
        try:
            r_p, _ = stats.pearsonr(perm_b, milo_lfcs)
            null_pearsons[i] = float(r_p) if not np.isnan(r_p) else 0.0
        except Exception:
            null_pearsons[i] = 0.0

    perm_p_spearman = float(np.mean(null_rhos >= rho_val)) if rho_val > 0 else float(np.mean(null_rhos <= rho_val))
    perm_p_pearson = float(np.mean(null_pearsons >= r_val)) if r_val > 0 else float(np.mean(null_pearsons <= r_val))

    summary = {
        "stratum": stratum,
        "cohort": clean_cohort,
        "cancer_type": cancer_type,
        "condition": condition,
        "n_cell_states": n_states,
        "spearman_rho": float(rho_val),
        "spearman_pvalue": float(rho_p),
        "spearman_perm_pvalue": perm_p_spearman,
        "pearson_r": float(r_val),
        "pearson_pvalue": float(r_p),
        "pearson_perm_pvalue": perm_p_pearson,
        "cosine_similarity": cosine_sim,
        "cohen_kappa": cohen_kappa,
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

    # 5b. Cell-level Cluster-Robust Regression Testing
    cells_file = config.milo_dir / "milopy_cell_level_scores.parquet"
    stat_test_records: list[dict[str, object]] = []

    if cells_file.exists() and df_default_metrics.height > 0:
        try:
            df_cells = pl.read_parquet(cells_file)
            df_conc_sub = df_default_metrics.select(["cell_state", "beta"])
            df_merged = df_cells.join(df_conc_sub, on="cell_state", how="inner").to_pandas()

            for cond_col in ["milo_logfc_combined", "milo_logfc_pre", "milo_logfc_post"]:
                if cond_col not in df_merged.columns:
                    continue
                sub_df = df_merged.dropna(subset=[cond_col, "beta"])
                if len(sub_df) < 50:
                    continue

                y = sub_df[cond_col].values
                X = sm.add_constant(sub_df["beta"].values)
                cluster_groups = sub_df["cell_state"].astype("category").cat.codes.values

                # Cluster-robust OLS
                ols_model = sm.OLS(y, X).fit(cov_type="cluster", cov_kwds={"groups": cluster_groups})
                beta_coef = float(ols_model.params[1])
                beta_se = float(ols_model.bse[1])
                beta_pval = float(ols_model.pvalues[1])
                r2 = float(ols_model.rsquared)

                stat_test_records.append(
                    {
                        "test_level": "Cell_Level_Cluster_Robust",
                        "condition": cond_col.replace("milo_logfc_", "").title(),
                        "n_observations": int(len(sub_df)),
                        "n_clusters": int(len(np.unique(cluster_groups))),
                        "coefficient": beta_coef,
                        "std_error": beta_se,
                        "p_value": beta_pval,
                        "r_squared": r2,
                    }
                )
        except Exception as exc:
            print(f"Notice: cell-level cluster robust regression skipped: {exc}")

    # Also record state-level summary into statistical tests table
    for r in df_summary_master.iter_rows(named=True):
        stat_test_records.append(
            {
                "test_level": "Cell_State_Permutation",
                "condition": f"{r['stratum']}_{r['condition']}",
                "n_observations": int(r["n_cell_states"]),
                "n_clusters": int(r["n_cell_states"]),
                "coefficient": float(r["spearman_rho"]),
                "std_error": np.nan,
                "p_value": float(r["spearman_perm_pvalue"]),
                "r_squared": float(r["pearson_r"]) ** 2 if not np.isnan(r["pearson_r"]) else 0.0,
            }
        )

    if stat_test_records:
        df_stat_tests = pl.DataFrame(stat_test_records)
        out_stat_tests = config.out_dir / "concordance_statistical_tests.parquet"
        df_stat_tests.write_parquet(out_stat_tests)
        print(f"Saved rigorous concordance statistical tests ({df_stat_tests.height} records) to: {out_stat_tests}")

    # 6. Multi-Resolution Benchmark Compilation — all references × all resolutions
    benchmark_records: list[dict[str, object]] = []

    search_dirs = [config.out_dir]
    nested_dir = Path("output/output/sade_feldman_deconv_validation")
    if nested_dir.exists() and nested_dir not in search_dirs:
        search_dirs.append(nested_dir)

    # Load condition number and cluster metrics across all references
    metrics_map: dict[str, dict[float, dict[str, object]]] = {}
    for s_dir in search_dirs:
        for mf in s_dir.glob("reference_resolution_metrics_*.parquet"):
            m_tag = mf.stem.replace("reference_resolution_metrics_", "")
            if m_tag not in metrics_map:
                metrics_map[m_tag] = {}
            for r in pl.read_parquet(mf).iter_rows(named=True):
                metrics_map[m_tag][float(r["resolution"])] = r

    name_map: dict[str, str] = {
        "sf": "Sade-Feldman",
        "comb": "Combined-Atlas",
        "jerby": "Jerby-Arnon",
        "maynard": "Maynard-NSCLC",
        "ma": "Ma-Liver",
        "yost": "Yost-BCC",
        "atlas": "PanCancer-Atlas",
        # criteria-based combinations
        "combo_melanoma": "Melanoma-Duo",
        "combo_plat_ss2": "SS2-Cross-Cancer",
        "combo_plat_10x": "10x-Cross-Cancer",
        "combo_cross_tissue": "Cross-Tissue-Pair",
        "combo_tri_ici": "Triple-ICI",
        # random combinations
        "random_pair1": "Random-Pair-1",
        "random_pair2": "Random-Pair-2",
        "random_triplet1": "Random-Triplet-1",
        "random_triplet2": "Random-Triplet-2",
    }

    # Discover all logistic regression result files across references and resolutions
    log_targets: list[tuple[str, str, float, Path]] = []

    # 1. Discover sf across all resolutions via glob (with legacy fallback at res=0.5)
    for s_dir in search_dirs:
        for lf in sorted(s_dir.glob("logistic_results_sf_res*.parquet")):
            m = re.match(r"logistic_results_sf_res([0-9\.]+)\.parquet", lf.name)
            if m:
                res_val = float(m.group(1))
                if not any(t[0] == "sf" and abs(t[2] - res_val) < 1e-4 for t in log_targets):
                    log_targets.append(("sf", "Sade-Feldman", res_val, lf))
    if not any(t[0] == "sf" for t in log_targets) and config.logistic_path.exists():
        log_targets.append(("sf", "Sade-Feldman", 0.5, config.logistic_path))

    # 2. All remaining references via glob (comb, jerby, maynard, ma, yost, combos, randoms)
    for s_dir in search_dirs:
        for lf in sorted(s_dir.glob("logistic_results_*_res*.parquet")):
            fname = lf.name
            m = re.match(r"logistic_results_([a-zA-Z0-9_\.-]+)_res([0-9\.]+)\.parquet", fname)
            if m:
                d_id, r_str = m.group(1), m.group(2)
                if d_id == "sf":
                    continue  # already handled above
                r_val = float(r_str)
                d_name = name_map.get(d_id, d_id.replace("_", "-").title())
                if not any(t[0] == d_id and abs(t[2] - r_val) < 1e-4 for t in log_targets):
                    log_targets.append((d_id, d_name, r_val, lf))

    print(f"\nCompiling benchmark metrics across {len(log_targets)} reference configurations...")
    for ref_id, ref_name, res, log_file in log_targets:
        df_log = pl.read_parquet(log_file)
        sub_mel = df_log.filter(pl.col("stratum") == "Melanoma")
        sub_pan = df_log.filter(pl.col("stratum") == "Pan-Cancer")
        cohort_rows = df_log.filter(pl.col("stratum").str.starts_with("Cohort_"))

        mel_multi_auc = float(sub_mel["multivariate_auc"][0]) if sub_mel.height > 0 and "multivariate_auc" in sub_mel.columns else 0.5
        pan_multi_auc = float(sub_pan["multivariate_auc"][0]) if sub_pan.height > 0 and "multivariate_auc" in sub_pan.columns else 0.5
        mean_cohort_auc = float(cohort_rows["multivariate_auc"].mean()) if cohort_rows.height > 0 and "multivariate_auc" in cohort_rows.columns else 0.5
        mean_univ_auc = float(sub_mel["auc"].mean()) if sub_mel.height > 0 else 0.5

        # Extract reference_category and reference metadata
        ref_metadata: dict[str, dict[str, Any]] = {
            "Sade-Feldman": {"category": "Single Dataset", "n_datasets": 1, "n_cells": 16288},
            "Jerby-Arnon": {"category": "Single Dataset", "n_datasets": 1, "n_cells": 7186},
            "Maynard-NSCLC": {"category": "Single Dataset", "n_datasets": 1, "n_cells": 3000},
            "Ma-Liver": {"category": "Single Dataset", "n_datasets": 1, "n_cells": 5115},
            "Yost-BCC": {"category": "Single Dataset", "n_datasets": 1, "n_cells": 3500},
            "Combined-Atlas": {"category": "Criteria-Combined", "n_datasets": 16, "n_cells": 41284},
            "Melanoma-Duo": {"category": "Criteria-Combined", "n_datasets": 2, "n_cells": 10686},
            "SS2-Cross-Cancer": {"category": "Criteria-Combined", "n_datasets": 2, "n_cells": 10186},
            "10x-Cross-Cancer": {"category": "Criteria-Combined", "n_datasets": 2, "n_cells": 8615},
            "Cross-Tissue-Pair": {"category": "Criteria-Combined", "n_datasets": 2, "n_cells": 8115},
            "Triple-ICI": {"category": "Criteria-Combined", "n_datasets": 3, "n_cells": 13686},
            "Random-Pair-1": {"category": "Random-Combined", "n_datasets": 2, "n_cells": 12301},
            "Random-Pair-2": {"category": "Random-Combined", "n_datasets": 2, "n_cells": 6500},
            "Random-Pair-3": {"category": "Random-Combined", "n_datasets": 2, "n_cells": 10686},
            "Random-Triplet-1": {"category": "Random-Combined", "n_datasets": 3, "n_cells": 15301},
            "Random-Triplet-2": {"category": "Random-Combined", "n_datasets": 3, "n_cells": 11615},
            "Random-Triplet-3": {"category": "Random-Combined", "n_datasets": 3, "n_cells": 15801},
            "All-Datasets-Combined": {"category": "Criteria-Combined", "n_datasets": 4, "n_cells": 18801},
            "Random-Quadruplet-1": {"category": "Random-Combined", "n_datasets": 4, "n_cells": 18801},
        }

        ref_info = ref_metadata.get(ref_name, {})
        ref_category: str = ref_info.get(
            "category",
            str(df_log["reference_category"][0])
            if "reference_category" in df_log.columns and df_log.height > 0
            else ("Single Dataset" if ref_id in ("sf", "jerby", "maynard", "ma", "yost") else "Unknown")
        )
        n_datasets: int = ref_info.get("n_datasets", 1)
        n_cells: int = ref_info.get("n_cells", 0)
        sub_frac: float = 1.0

        # Check if reference_resolution_metrics had recorded metadata
        ref_meta = metrics_map.get(ref_id, {}).get(res, {})
        if not ref_meta:
            # check without res prefix or exact match
            ref_meta = metrics_map.get(ref_id, {}).get(0.5, {})
        if ref_meta:
            if "n_datasets" in ref_meta and ref_meta["n_datasets"]:
                n_datasets = int(ref_meta["n_datasets"])
            if "n_cells" in ref_meta and ref_meta["n_cells"]:
                n_cells = int(ref_meta["n_cells"])
            if "subsample_fraction" in ref_meta and ref_meta["subsample_fraction"]:
                sub_frac = float(ref_meta["subsample_fraction"])
            if "category" in ref_meta and ref_meta["category"]:
                ref_category = str(ref_meta["category"])
            if "reference_type" in ref_meta and ref_meta["reference_type"]:
                ref_name = str(ref_meta["reference_type"])

        # Single-cell Milo DA concordance (where available for Sade-Feldman)
        rho_val = float("nan")
        pct_conc = float("nan")
        if ref_id == "sf":
            milo_res_file = config.milo_dir / f"milopy_cell_state_da_Combined_res{res}.parquet"
            if not milo_res_file.exists() and abs(res - 0.5) < 1e-4:
                milo_res_file = config.milo_dir / "milopy_cell_state_da.parquet"
            if milo_res_file.exists():
                df_milo_res = pl.read_parquet(milo_res_file)
                j_sf = sub_mel.join(df_milo_res, on="cell_state", how="inner")
                if j_sf.height >= 3:
                    try:
                        r_val, _ = stats.spearmanr(j_sf["beta"].to_numpy(), j_sf["milo_mean_logfc"].to_numpy())
                        rho_val = float(r_val) if not np.isnan(r_val) else 0.0
                    except Exception:
                        rho_val = 0.0
                    n_conc = sum(1 for b, m in zip(j_sf["beta"], j_sf["milo_mean_logfc"]) if (b > 0 and m > 0) or (b < 0 and m < 0))
                    pct_conc = float(n_conc / j_sf.height) * 100.0

        cond_num = float(ref_meta.get("condition_number", np.nan))
        n_clusters = int(ref_meta.get("n_clusters", sub_mel.height))
        n_sig_genes = int(ref_meta.get("n_signature_genes", 0))

        benchmark_records.append(
            {
                "reference_type": ref_name,
                "reference_category": ref_category,
                "resolution": float(res),
                "n_clusters": n_clusters,
                "n_signature_genes": n_sig_genes,
                "condition_number": cond_num,
                "spearman_rho": rho_val,
                "concordance_percentage": pct_conc,
                "melanoma_multivariate_auc": mel_multi_auc,
                "pancancer_multivariate_auc": pan_multi_auc,
                "mean_cohort_multivariate_auc": mean_cohort_auc,
                "mean_univariate_auc": mean_univ_auc,
                "n_datasets": n_datasets,
                "n_cells": n_cells,
                "subsample_fraction": sub_frac,
            }
        )

    if benchmark_records:
        df_bench = pl.DataFrame(benchmark_records)
        out_bench = config.out_dir / "multi_resolution_benchmark_summary.parquet"
        df_bench.write_parquet(out_bench)
        if nested_dir.exists() and nested_dir != config.out_dir:
            df_bench.write_parquet(nested_dir / "multi_resolution_benchmark_summary.parquet")
        print(f"\nSaved multi-resolution benchmark summary ({df_bench.height} configurations) to: {out_bench}")
        print(df_bench)

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
