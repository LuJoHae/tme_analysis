#!/usr/bin/env python3
"""
Step 3: Logistic Regression of Inferred Cell State Fractions vs. Immunotherapy Response.
Computes univariate logistic regression effect sizes, Odds Ratios (95% CI), p-values, FDR, and ROC AUC.
Evaluates melanoma cohorts (matching tissue) and pan-cancer cohorts.
Outputs logistic_regression_results.parquet and sample_fractions_with_response.parquet.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import stats  # type: ignore
from sklearn.linear_model import LogisticRegression  # type: ignore
from sklearn.metrics import roc_auc_score  # type: ignore

# Local repo packages
sys.path.append(str(Path(__file__).resolve().parent.parent.parent / "packages"))
import datalair  # type: ignore
import ici_datasets  # type: ignore


RESPONSE_VALUE_MAP: Final[dict[str, int]] = {
    "complete response": 1,
    "partial response": 1,
    "cr": 1,
    "pr": 1,
    "r": 1,
    "responder": 1,
    "stable disease": 0,
    "progressive disease": 0,
    "sd": 0,
    "pd": 0,
    "nr": 0,
    "non-responder": 0,
    "nonresponder": 0,
}

RESPONSE_COL_CANDIDATES: Final[tuple[str, ...]] = (
    "response",
    "best_response",
    "recist",
    "clinical_response",
    "overall_response",
)


class LogRegConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    fractions_path: Path
    lair_dir: Path
    out_dir: Path
    output_name: str = "logistic_regression_results.parquet"
    benchmark_mode: bool = False


def load_cohort_clinical_response(cohort_name: str, lair: datalair.Lair) -> Result[pd.DataFrame, str]:
    """Extract clinical response from cBioPortal clinical metadata."""
    try:
        ds_class = ici_datasets.cbioportal_datasets.CBioPortalDataset
        data_dir = ici_datasets.cbioportal_datasets.get_dataset_dir(lair, ds_class, cohort_name)
        sample_file = data_dir / "data_clinical_sample.txt"
        if not sample_file.exists():
            return Failure(f"Clinical file not found: {sample_file}")

        df_clin = pd.read_csv(sample_file, sep="\t", skiprows=4)
        sample_col = next(
            (c for c in df_clin.columns if c.strip().lower() in ("sample_id", "sample id")),
            None,
        )
        if sample_col is None:
            return Failure(f"No sample_id column in clinical file for {cohort_name}")

        df_clin = df_clin.set_index(sample_col)
        df_clin.index = df_clin.index.astype(str)
        df_clin = df_clin[~df_clin.index.duplicated(keep="first")]

        # Locate response column
        resp_col = None
        for col in df_clin.columns:
            if col.strip().lower() in RESPONSE_COL_CANDIDATES:
                resp_col = col
                break

        if resp_col is None:
            return Failure(f"No response column recognized in {cohort_name}")

        # Map response to binary (1 = R, 0 = NR)
        clean_resp = (
            df_clin[resp_col].astype(str).str.strip().str.lower().map(RESPONSE_VALUE_MAP)
        )
        df_valid = pd.DataFrame({"response": clean_resp}).dropna()
        df_valid["response"] = df_valid["response"].astype(int)

        return Success(df_valid)
    except Exception as exc:
        return Failure(f"Failed to parse clinical for {cohort_name}: {exc}")


def fit_univariate_logistic(
    x_raw: np.ndarray,
    y: np.ndarray,
) -> dict[str, float]:
    """Fit univariate logistic regression and calculate Wald stats, OR, CI, and AUC."""
    n = len(y)
    n_r = int(np.sum(y == 1))
    n_nr = int(np.sum(y == 0))

    if n_r == 0 or n_nr == 0:
        return {
            "beta": 0.0,
            "se": np.nan,
            "or": 1.0,
            "or_ci_lower": 1.0,
            "or_ci_upper": 1.0,
            "p_value": 1.0,
            "beta_z": 0.0,
            "or_z": 1.0,
            "or_z_lower": 1.0,
            "or_z_upper": 1.0,
            "auc": 0.5,
            "n_total": float(n),
            "n_responders": float(n_r),
            "n_non_responders": float(n_nr),
        }

    # Standardize x for standardized effect size
    std_x = np.std(x_raw)
    x_z = (x_raw - np.mean(x_raw)) / (std_x if std_x > 0 else 1.0)

    # Scikit-learn LogisticRegression with no penalty for exact MLE
    try:
        clf = LogisticRegression(penalty=None, solver="lbfgs", max_iter=500)
        clf.fit(x_raw.reshape(-1, 1), y)
        beta = float(clf.coef_[0, 0])
        p_hat = clf.predict_proba(x_raw.reshape(-1, 1))[:, 1]
        p_hat = np.clip(p_hat, 1e-6, 1.0 - 1e-6)

        # Variance-covariance matrix of beta: inv(X^T W X)
        W = p_hat * (1.0 - p_hat)
        X_design = np.column_stack([np.ones_like(x_raw), x_raw])
        fisher_info = X_design.T @ (W[:, np.newaxis] * X_design)
        cov_mat = np.linalg.pinv(fisher_info)
        se_beta = float(np.sqrt(max(cov_mat[1, 1], 1e-12)))

        wald_stat = (beta / se_beta) ** 2
        p_val = float(1.0 - stats.chi2.cdf(wald_stat, df=1))

        # Standardized fit
        clf_z = LogisticRegression(penalty=None, solver="lbfgs", max_iter=500)
        clf_z.fit(x_z.reshape(-1, 1), y)
        beta_z = float(clf_z.coef_[0, 0])
        p_hat_z = np.clip(clf_z.predict_proba(x_z.reshape(-1, 1))[:, 1], 1e-6, 1.0 - 1e-6)
        W_z = p_hat_z * (1.0 - p_hat_z)
        X_design_z = np.column_stack([np.ones_like(x_z), x_z])
        cov_mat_z = np.linalg.pinv(X_design_z.T @ (W_z[:, np.newaxis] * X_design_z))
        se_beta_z = float(np.sqrt(max(cov_mat_z[1, 1], 1e-12)))

        or_val = float(np.exp(np.clip(beta, -15, 15)))
        or_ci_lower = float(np.exp(np.clip(beta - 1.96 * se_beta, -15, 15)))
        or_ci_upper = float(np.exp(np.clip(beta + 1.96 * se_beta, -15, 15)))

        or_z = float(np.exp(np.clip(beta_z, -15, 15)))
        or_z_lower = float(np.exp(np.clip(beta_z - 1.96 * se_beta_z, -15, 15)))
        or_z_upper = float(np.exp(np.clip(beta_z + 1.96 * se_beta_z, -15, 15)))

        auc = float(roc_auc_score(y, x_raw))
    except Exception:
        # Robust fallback using Mann-Whitney U test and simple rank stats
        u_stat, p_val = stats.mannwhitneyu(x_raw[y == 1], x_raw[y == 0], alternative="two-sided")
        auc = float(u_stat / (n_r * n_nr))
        beta = 0.0
        se_beta = np.nan
        or_val = 1.0
        or_ci_lower = 1.0
        or_ci_upper = 1.0
        beta_z = 0.0
        or_z = 1.0
        or_z_lower = 1.0
        or_z_upper = 1.0

    return {
        "beta": beta,
        "se": se_beta,
        "or": or_val,
        "or_ci_lower": or_ci_lower,
        "or_ci_upper": or_ci_upper,
        "p_value": p_val,
        "beta_z": beta_z,
        "or_z": or_z,
        "or_z_lower": or_z_lower,
        "or_z_upper": or_z_upper,
        "auc": auc,
        "n_total": float(n),
        "n_responders": float(n_r),
        "n_non_responders": float(n_nr),
    }


def compute_fdr(p_values: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR adjustment."""
    m = len(p_values)
    if m == 0:
        return []
    p_arr = np.asarray(p_values)
    sorted_idx = np.argsort(p_arr)
    sorted_p = p_arr[sorted_idx]
    fdr_sorted = np.zeros(m)
    for i in range(m):
        rank = i + 1
        fdr_sorted[i] = sorted_p[i] * m / rank
    # Monotonicity enforcement from back
    for i in range(m - 2, -1, -1):
        fdr_sorted[i] = min(fdr_sorted[i], fdr_sorted[i + 1])
    fdr_sorted = np.clip(fdr_sorted, 0.0, 1.0)
    orig_fdr = np.zeros(m)
    orig_fdr[sorted_idx] = fdr_sorted
    return orig_fdr.tolist()


def analyze_single_fractions_matrix(
    df_fracs: pl.DataFrame,
    df_clin_pl: pl.DataFrame,
    ref_type: str,
    resolution: float,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Fit univariate and multivariate models on a single deconvolution fractions table."""
    metadata_cols = {"sample_id", "cohort", "cancer_type"}
    cell_states = tuple(c for c in df_fracs.columns if c not in metadata_cols)
    cohorts = tuple(sorted(df_fracs["cohort"].unique().to_list()))

    df_merged = df_fracs.join(df_clin_pl, on="sample_id", how="inner")
    strata: list[tuple[str, pl.DataFrame]] = [
        ("Melanoma", df_merged.filter(pl.col("cancer_type") == "Melanoma")),
        ("Pan-Cancer", df_merged),
    ]
    for ch in cohorts:
        sub = df_merged.filter(pl.col("cohort") == ch)
        if sub.height >= 15 and sub.filter(pl.col("response") == 1).height >= 3:
            strata.append((f"Cohort_{ch}", sub))

    results_records: list[dict[str, object]] = []

    for stratum_name, stratum_df in strata:
        if stratum_df.height < 10:
            continue

        y_arr = stratum_df["response"].to_numpy().astype(int)
        stratum_pvals: list[float] = []
        stratum_temp_records: list[dict[str, object]] = []

        # Fit multivariate model using all cell states with L2 regularization
        if len(cell_states) > 1 and len(np.unique(y_arr)) == 2:
            try:
                X_mat = stratum_df.select(cell_states).to_numpy().astype(np.float64)
                clf_multi = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs", max_iter=500)
                clf_multi.fit(X_mat, y_arr)
                p_multi = clf_multi.predict_proba(X_mat)[:, 1]
                multi_auc = float(roc_auc_score(y_arr, p_multi))
            except Exception:
                multi_auc = 0.5
        else:
            multi_auc = 0.5

        for cs in cell_states:
            x_arr = stratum_df[cs].to_numpy().astype(np.float64)
            fit_res = fit_univariate_logistic(x_arr, y_arr)

            rec: dict[str, object] = {
                "reference_type": ref_type,
                "resolution": float(resolution),
                "stratum": stratum_name,
                "cell_state": cs,
                "multivariate_auc": multi_auc,
                **fit_res,
            }
            stratum_temp_records.append(rec)
            stratum_pvals.append(fit_res["p_value"])

        # Compute Benjamini-Hochberg FDR within this stratum
        fdrs = compute_fdr(stratum_pvals)
        for rec, fdr in zip(stratum_temp_records, fdrs):
            rec["fdr"] = fdr
            rec["log10_pval"] = float(-np.log10(max(rec["p_value"], 1e-12)))  # type: ignore
            rec["significant_fdr01"] = bool(fdr < 0.1)
            results_records.append(rec)

    return df_merged, pl.DataFrame(results_records)


def run_logistic_regression_pipeline(config: LogRegConfig) -> Result[Path, str]:
    """Run logistic regression across single or all benchmark fraction matrices."""
    config.out_dir.mkdir(parents=True, exist_ok=True)
    lair = datalair.Lair(str(config.lair_dir))

    # Identify cohorts from fractions
    sample_file = config.fractions_path if config.fractions_path.exists() else None
    if sample_file is None:
        sample_file = config.out_dir / "deconv_fractions.parquet"
        if not sample_file.exists():
            sample_candidates = list(config.out_dir.glob("deconv_fractions_*.parquet"))
            if sample_candidates:
                sample_file = sample_candidates[0]

    if sample_file is None or not sample_file.exists():
        return Failure(f"No deconvolution fractions files found in {config.out_dir}")

    df_sample = pl.read_parquet(sample_file)
    cohorts = tuple(sorted(df_sample["cohort"].unique().to_list()))

    clin_records: list[dict[str, object]] = []
    for cohort in cohorts:
        clin_res = load_cohort_clinical_response(cohort, lair)
        match clin_res:
            case Failure(err):
                print(f"Warning: {cohort} clinical response missing: {err}")
            case Success(df_clin):
                for sid, row in df_clin.iterrows():
                    clin_records.append({"sample_id": str(sid), "response": int(row["response"])})

    if not clin_records:
        return Failure("No clinical response records could be loaded from cBioPortal.")

    df_clin_pl = pl.DataFrame(clin_records)
    print(f"Loaded clinical response for {df_clin_pl.height} unique patient samples.")

    if config.benchmark_mode:
        print("\n=== RUNNING MULTI-RESOLUTION LOGISTIC REGRESSION BENCHMARK ===")
        targets: list[tuple[str, str, float, Path, str]] = []
        for res in (0.5, 1.0, 1.5, 2.0):
            sf_file = config.out_dir / f"deconv_fractions_sf_res{res}.parquet"
            if sf_file.exists():
                targets.append(("sf_res" + str(res), "Sade-Feldman", res, sf_file, f"logistic_results_sf_res{res}.parquet"))
            comb_file = config.out_dir / f"deconv_fractions_comb_res{res}.parquet"
            if comb_file.exists():
                targets.append(("comb_res" + str(res), "Combined-Atlas", res, comb_file, f"logistic_results_comb_res{res}.parquet"))

        if not targets:
            if config.fractions_path.exists():
                targets.append(("primary", "Sade-Feldman", 0.5, config.fractions_path, "logistic_regression_results.parquet"))
            else:
                return Failure("No fraction parquets found for benchmarking.")

        all_results_dfs: list[pl.DataFrame] = []
        last_out: Path = config.out_dir / "logistic_regression_results.parquet"

        for tag, ref_type, res, frac_path, out_name in targets:
            print(f"--- Fitting logistic models for {ref_type} (res={res}) ---")
            df_fracs = pl.read_parquet(frac_path)
            df_merged, df_res = analyze_single_fractions_matrix(df_fracs, df_clin_pl, ref_type, res)

            out_p = config.out_dir / out_name
            df_res.write_parquet(out_p)
            all_results_dfs.append(df_res)
            last_out = out_p

            if tag == "sf_res0.5":
                df_merged.write_parquet(config.out_dir / "sample_fractions_with_response.parquet")
                df_res.write_parquet(config.out_dir / "logistic_regression_results.parquet")
                print("Saved default logistic_regression_results.parquet (from sf_res0.5).")

        df_all_results = pl.concat(all_results_dfs, how="vertical")
        out_all_res = config.out_dir / "logistic_regression_results_all_resolutions.parquet"
        df_all_results.write_parquet(out_all_res)
        print(f"\nSaved combined multi-resolution logistic results ({df_all_results.height} rows) to: {out_all_res}")
        return Success(last_out)

    else:
        # Single execution
        if not config.fractions_path.exists():
            return Failure(f"Fractions file not found: {config.fractions_path}")

        df_fracs = pl.read_parquet(config.fractions_path)
        df_merged, df_results = analyze_single_fractions_matrix(
            df_fracs, df_clin_pl, "Sade-Feldman", 0.5
        )

        out_merged = config.out_dir / "sample_fractions_with_response.parquet"
        df_merged.write_parquet(out_merged)

        out_results = config.out_dir / config.output_name
        df_results.write_parquet(out_results)
        print(f"\nCompleted logistic regression analysis. Saved to: {out_results}")
        return Success(out_results)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 3: Univariate logistic regression of cell state fractions vs immunotherapy response."
    )
    parser.add_argument(
        "--fractions",
        type=str,
        default="output/sade_feldman_deconv_validation/deconv_fractions.parquet",
        help="Path to deconv_fractions.parquet from Step 2",
    )
    parser.add_argument(
        "--lair-dir",
        type=str,
        default="/storage/halu/lair",
        help="Path to datalair directory",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="output/sade_feldman_deconv_validation",
        help="Directory to save output parquets",
    )
    parser.add_argument(
        "--output-name",
        type=str,
        default="logistic_regression_results.parquet",
        help="Filename for single output results parquet",
    )
    parser.add_argument(
        "--benchmark-mode",
        action="store_true",
        help="Execute logistic regression across all 8 multi-resolution fractions matrices",
    )
    args = parser.parse_args()

    config = LogRegConfig(
        fractions_path=Path(args.fractions),
        lair_dir=Path(args.lair_dir),
        out_dir=Path(args.out_dir),
        output_name=args.output_name,
        benchmark_mode=args.benchmark_mode,
    )

    match run_logistic_regression_pipeline(config):
        case Success(out_file):
            print(f"Step 3 finished successfully: {out_file}")
            sys.exit(0)
        case Failure(err):
            print(f"Step 3 failed with error: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
