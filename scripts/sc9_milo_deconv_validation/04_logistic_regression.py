#!/usr/bin/env python3
"""Step 3: Logistic Regression & Response Prediction Performance Evaluation.
Evaluates the predictive capacity of cell state fractions estimated across all
multi-resolution references against ICB response across melanoma and pan-cancer cohorts.
Outputs logistic_regression_results_all_resolutions.parquet.

Follows functional programming standards: pure functions, immutable data structures,
and explicit error handling using Returns Result monads.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Final

import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import stats  # type: ignore[import-untyped]
from sklearn.linear_model import LogisticRegression  # type: ignore[import-untyped]
from sklearn.metrics import roc_auc_score  # type: ignore[import-untyped]

# Ensure repo packages are discoverable
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(REPO_ROOT / "packages" / "tme_datasets" / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "packages" / "tme_datasets" / "src"))

from tme_datasets import load_dataset

RESPONSE_VALUE_MAP: Final[dict[str, int]] = {
    "r": 1,
    "responder": 1,
    "cr": 1,
    "pr": 1,
    "nr": 0,
    "non-responder": 0,
    "sd": 0,
    "pd": 0,
    "nonresponder": 0,
}

RESPONSE_COL_CANDIDATES: Final[tuple[str, ...]] = (
    "response",
    "best_response",
    "recist",
    "clinical_response",
    "overall_response",
    "response_binary",
)


class LogRegConfig(BaseModel):
    """Immutable configuration for logistic regression evaluation."""

    model_config = ConfigDict(frozen=True)

    deconv_dir: Path
    out_dir: Path


def load_cohort_clinical_response(cohort_name: str) -> Result[pd.DataFrame, str]:
    """Extract binary clinical response (1 = R, 0 = NR) from tme_datasets."""
    match load_dataset(cohort_name, auto_download=False, normalize_ensembl=False):
        case Success(adata):
            try:
                obs = adata.obs.copy()
                obs.index = obs.index.astype(str)
                obs = obs[~obs.index.duplicated(keep="first")]
                if "response_binary" in obs.columns:
                    raw_bin = pd.to_numeric(obs["response_binary"], errors="coerce")
                    valid = raw_bin[raw_bin.isin([0.0, 1.0])].astype(int)
                    return Success(pd.DataFrame({"response": valid}))

                resp_col = next(
                    (c for c in obs.columns if c.strip().lower() in RESPONSE_COL_CANDIDATES),
                    None,
                )
                if resp_col is None:
                    return Failure(f"No response column recognized in {cohort_name}")
                clean_resp = obs[resp_col].astype(str).str.strip().str.lower().map(RESPONSE_VALUE_MAP)
                df_valid = pd.DataFrame({"response": clean_resp}).dropna()
                df_valid["response"] = df_valid["response"].astype(int)
                return Success(df_valid)
            except Exception as exc:
                return Failure(f"Failed to parse clinical response for {cohort_name}: {exc}")
        case Failure(tme_err):
            return Failure(f"Failed to load {cohort_name}: {tme_err}")


def fit_univariate_logistic(x_raw: np.ndarray, y: np.ndarray) -> dict[str, float]:
    """Fit univariate logistic regression with Wald test statistics and AUC."""
    n = len(y)
    n_r = int(np.sum(y == 1))
    n_nr = int(np.sum(y == 0))

    if n_r == 0 or n_nr == 0 or np.all(x_raw == x_raw[0]):
        return {
            "beta": 0.0,
            "se": np.nan,
            "or": 1.0,
            "or_ci_lower": 1.0,
            "or_ci_upper": 1.0,
            "p_value": 1.0,
            "auc": 0.5,
            "n_total": float(n),
            "n_responders": float(n_r),
            "n_non_responders": float(n_nr),
        }

    try:
        clf = LogisticRegression(penalty=None, solver="lbfgs", max_iter=500)
        clf.fit(x_raw.reshape(-1, 1), y)
        beta = float(clf.coef_[0, 0])
        p_hat = np.clip(clf.predict_proba(x_raw.reshape(-1, 1))[:, 1], 1e-6, 1.0 - 1e-6)

        W = p_hat * (1.0 - p_hat)
        X_design = np.column_stack([np.ones_like(x_raw), x_raw])
        fisher_info = X_design.T @ (W[:, np.newaxis] * X_design)
        cov_mat = np.linalg.pinv(fisher_info)
        se_beta = float(np.sqrt(max(cov_mat[1, 1], 1e-12)))

        wald_stat = (beta / se_beta) ** 2
        p_val = float(1.0 - stats.chi2.cdf(wald_stat, df=1))

        or_val = float(np.exp(np.clip(beta, -15, 15)))
        or_ci_lower = float(np.exp(np.clip(beta - 1.96 * se_beta, -15, 15)))
        or_ci_upper = float(np.exp(np.clip(beta + 1.96 * se_beta, -15, 15)))
        auc = float(roc_auc_score(y, x_raw))

    except Exception:
        u_stat, p_val = stats.mannwhitneyu(x_raw[y == 1], x_raw[y == 0], alternative="two-sided")
        auc = float(u_stat / (n_r * n_nr))
        beta = 0.0
        se_beta = np.nan
        or_val = 1.0
        or_ci_lower = 1.0
        or_ci_upper = 1.0

    return {
        "beta": beta,
        "se": se_beta,
        "or": or_val,
        "or_ci_lower": or_ci_lower,
        "or_ci_upper": or_ci_upper,
        "p_value": p_val,
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
    for i in range(m - 2, -1, -1):
        fdr_sorted[i] = min(fdr_sorted[i], fdr_sorted[i + 1])
    fdr_sorted = np.clip(fdr_sorted, 0.0, 1.0)
    orig_fdr = np.zeros(m)
    orig_fdr[sorted_idx] = fdr_sorted
    return orig_fdr.tolist()


def get_reference_category(name: str) -> str:
    """Categorizes reference strategies cleanly."""
    lower = name.lower()
    if "_sub" in lower or "subsample" in lower or "%" in lower:
        return "Cell-Subsample"
    if lower.startswith("random"):
        return "Random-Combined"
    if any(k in lower for k in ("duo", "trio", "quad", "quartet", "combined", "atlas")):
        return "Criteria-Combined"
    return "Single Dataset"


def run_logistic_regression(config: LogRegConfig) -> Result[Path, str]:
    """Discovers all deconvolution parquets and runs univariate + multivariate regressions."""
    config.out_dir.mkdir(parents=True, exist_ok=True)

    frac_files = sorted(config.deconv_dir.glob("deconv_fractions_*.parquet"))
    if not frac_files:
        return Failure(f"No deconv_fractions_*.parquet files found in {config.deconv_dir}")

    # 1. Discover unique cohorts and load clinical response labels
    sample_df = pl.read_parquet(frac_files[0])
    cohorts = sorted(sample_df["cohort"].unique().to_list())
    print(f"Loading clinical response for {len(cohorts)} cohorts...")

    clin_dfs: list[pd.DataFrame] = []
    for ch in cohorts:
        match load_cohort_clinical_response(ch):
            case Success(c_df):
                c_df["sample_id"] = c_df.index
                c_df["cohort"] = ch
                clin_dfs.append(c_df)
                print(f"  Loaded {ch}: {len(c_df)} annotated samples ({c_df['response'].sum()} R / {(1-c_df['response']).sum()} NR)")
            case Failure(err):
                print(f"  Warning: {err}")

    if not clin_dfs:
        return Failure("No clinical response data could be loaded.")

    master_clin = pd.concat(clin_dfs, ignore_index=True)
    df_clin_pl = pl.DataFrame(master_clin[["sample_id", "response"]])

    all_results: list[pl.DataFrame] = []
    print(f"\nRunning logistic regression across {len(frac_files)} deconvolution configurations...")

    for ff in frac_files:
        m = re.match(r"^deconv_fractions_(.+)_res([0-9\.]+)\.parquet$", ff.name)
        if not m:
            continue
        ref_id, res_str = m.group(1), m.group(2)
        res_val = float(res_str)

        out_single = config.out_dir / f"logistic_regression_results_{ref_id}_res{res_str}.parquet"
        if out_single.exists() and out_single.stat().st_size > 500:
            all_results.append(pl.read_parquet(out_single))
            continue

        df_fracs = pl.read_parquet(ff)
        metadata_cols = {"sample_id", "cohort", "cancer_type"}
        cell_states = tuple(c for c in df_fracs.columns if c not in metadata_cols)
        df_merged = df_fracs.join(df_clin_pl, on="sample_id", how="inner")

        strata: list[tuple[str, pl.DataFrame]] = [
            ("Melanoma", df_merged.filter(pl.col("cancer_type") == "Melanoma")),
            ("Pan-Cancer", df_merged),
        ]
        for ch in cohorts:
            sub = df_merged.filter(pl.col("cohort") == ch)
            if sub.height >= 12 and sub.filter(pl.col("response") == 1).height >= 2:
                strata.append((f"Cohort_{ch}", sub))

        records: list[dict[str, object]] = []

        for stratum_name, stratum_df in strata:
            if stratum_df.height < 8:
                continue

            y_arr = stratum_df["response"].to_numpy().astype(int)
            stratum_pvals: list[float] = []
            stratum_temp_records: list[dict[str, object]] = []

            # Multivariate L2 logistic regression
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
                    "reference_type": ref_id,
                    "reference_category": get_reference_category(ref_id),
                    "resolution": res_val,
                    "stratum": stratum_name,
                    "cell_state": cs,
                    "multivariate_auc": multi_auc,
                    **fit_res,
                }
                stratum_temp_records.append(rec)
                stratum_pvals.append(fit_res["p_value"])

            fdrs = compute_fdr(stratum_pvals)
            for rec, fdr in zip(stratum_temp_records, fdrs):
                rec["fdr"] = fdr
                rec["log10_pval"] = float(-np.log10(max(rec["p_value"], 1e-12)))  # type: ignore
                rec["significant_fdr01"] = bool(fdr < 0.1)
                records.append(rec)

        if records:
            df_res = pl.DataFrame(records)
            out_single = config.out_dir / f"logistic_regression_results_{ref_id}_res{res_str}.parquet"
            df_res.write_parquet(out_single)
            all_results.append(df_res)

    if not all_results:
        return Failure("No logistic regression results could be computed.")

    df_master = pl.concat(all_results, how="vertical")
    out_master = config.out_dir / "logistic_regression_results_all_resolutions.parquet"
    df_master.write_parquet(out_master)
    print(f"\nSaved master logistic regression results ({df_master.height} rows) -> {out_master}")

    return Success(out_master)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 3: Fit univariate and multivariate logistic regression predicting ICB response."
    )
    parser.add_argument(
        "--deconv-dir",
        type=str,
        default="output/sc9_milo_deconv_validation/02_deconvolution",
        help="Directory containing deconv_fractions_*.parquet files",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="output/sc9_milo_deconv_validation/03_logistic_regression",
        help="Directory to save logistic regression results",
    )
    args = parser.parse_args()

    config = LogRegConfig(
        deconv_dir=Path(args.deconv_dir),
        out_dir=Path(args.out_dir),
    )

    match run_logistic_regression(config):
        case Success(p):
            print(f"Logistic regression succeeded: {p}")
            sys.exit(0)
        case Failure(err):
            print(f"Logistic regression failed: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
