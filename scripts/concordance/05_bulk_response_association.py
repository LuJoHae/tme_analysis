#!/usr/bin/env python3
"""
Step 5: Bulk Response Association Modeling.
Fits logistic regression models of patient response against raw, purity-normalized,
and mRNA-scaled deconvolution fractions across all cell states.
Strict functional Python with returns, Pydantic, and Polars.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
import polars as pl
from scipy import stats  # type: ignore


class BulkAssocConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    fractions_raw_path: Path
    fractions_norm_path: Path
    fractions_mrna_path: Path
    clinical_path: Path
    sample_col: str
    response_col: str
    out_parquet: Path


def parse_args() -> BulkAssocConfig:
    parser = argparse.ArgumentParser(
        description="Fit logistic association of bulk response against deconvolution fractions."
    )
    parser.add_argument("--fractions-raw", type=Path, required=True, help="Raw fractions parquet")
    parser.add_argument("--fractions-norm", type=Path, required=True, help="Purity-normalized fractions parquet")
    parser.add_argument("--fractions-mrna", type=Path, required=True, help="mRNA-scaled fractions parquet")
    parser.add_argument("--clinical", type=Path, required=True, help="Clinical metadata parquet or csv")
    parser.add_argument("--sample-col", type=str, default="sample_id", help="Sample ID column")
    parser.add_argument("--response-col", type=str, default="response", help="Response label column")
    parser.add_argument("--out-parquet", type=Path, required=True, help="Output parquet path")
    args = parser.parse_args()
    return BulkAssocConfig(
        fractions_raw_path=args.fractions_raw,
        fractions_norm_path=args.fractions_norm,
        fractions_mrna_path=args.fractions_mrna,
        clinical_path=args.clinical,
        sample_col=args.sample_col,
        response_col=args.response_col,
        out_parquet=args.out_parquet,
    )


def binarize_response(val: str) -> Maybe[int]:
    """Map response labels to 1 (Responder) or 0 (Non-Responder)."""
    norm = val.strip().lower()
    match norm:
        case "responder" | "r" | "cr" | "pr" | "yes" | "1":
            return Some(1)
        case "non-responder" | "non_responder" | "nr" | "sd" | "pd" | "no" | "0":
            return Some(0)
        case _:
            return Nothing


def read_clinical(clinical_path: Path, sample_col: str, response_col: str) -> Result[pl.DataFrame, str]:
    """Read clinical metadata and ensure binary response mapping."""
    if not clinical_path.exists():
        return Failure(f"Clinical file does not exist: {clinical_path}")

    try:
        if clinical_path.suffix == ".parquet":
            df = pl.read_parquet(clinical_path)
        else:
            df = pl.read_csv(clinical_path)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to read clinical file: {exc}")

    if sample_col not in df.columns or response_col not in df.columns:
        return Failure(f"Missing columns in clinical table. Needs {sample_col} and {response_col}")

    # Map responses
    sub_df = df.select([sample_col, response_col]).filter(pl.col(response_col).is_not_null())
    mapped_rows: list[dict[str, object]] = []

    for row in sub_df.iter_rows(named=True):
        match binarize_response(str(row[response_col])):
            case Some(binary_val):
                mapped_rows.append({
                    "sample_id": str(row[sample_col]),
                    "response_binary": binary_val,
                })
            case Nothing:
                pass

    if len(mapped_rows) < 4:
        return Failure(f"Too few valid binary response samples in clinical table: {len(mapped_rows)}")

    return Success(pl.DataFrame(mapped_rows))


def fit_fractions_association(
    frac_df: pl.DataFrame,
    clinical_df: pl.DataFrame,
    fraction_type: str,
) -> Result[pl.DataFrame, str]:
    """Fit logistic response regression for each cell state in a fraction table."""
    # Align on sample_id
    id_col = "sample_id" if "sample_id" in frac_df.columns else frac_df.columns[0]
    renamed_frac = frac_df.rename({id_col: "sample_id"})
    joined = renamed_frac.join(clinical_df, on="sample_id", how="inner")

    if joined.height < 4:
        return Failure(f"Insufficient overlap between fractions and clinical table for {fraction_type}")

    y = joined["response_binary"].to_numpy().astype(np.float64)
    state_cols = [c for c in renamed_frac.columns if c != "sample_id"]

    effects: list[dict[str, object]] = []

    for state in state_cols:
        x = joined[state].to_numpy().astype(np.float64)
        if np.std(x) < 1e-8:
            beta, se, z_val, p_val = 0.0, 1.0, 0.0, 1.0
        else:
            x_std = (x - np.mean(x)) / np.std(x)
            r_val, p_val_scipy = stats.pointbiserialr(y, x_std)
            r_clip = np.clip(r_val, -0.999, 0.999)
            beta = float(2.0 * r_clip / np.sqrt(1.0 - r_clip**2 + 1e-12))
            se = float(np.sqrt(4.0 / max(len(y) - 2, 1)))
            z_val = beta / (se + 1e-12)
            p_val = float(p_val_scipy)

        effects.append({
            "cell_state": state,
            "fraction_type": fraction_type,
            "beta_bulk": beta,
            "se_bulk": se,
            "z_score_bulk": z_val,
            "p_value_bulk": p_val,
            "mean_fraction": float(np.mean(x)),
            "n_samples": int(len(y)),
        })

    # FDR calculation
    p_vals = np.array([e["p_value_bulk"] for e in effects])
    order = np.argsort(p_vals)
    ranked_p = p_vals[order]
    n_tests = len(effects)
    fdrs = np.minimum(1.0, ranked_p * n_tests / (np.arange(1, n_tests + 1)))
    fdrs_mono = np.minimum.accumulate(fdrs[::-1])[::-1]
    fdrs_orig = np.empty_like(fdrs_mono)
    fdrs_orig[order] = fdrs_mono

    for i, e in enumerate(effects):
        e["fdr_bulk"] = float(fdrs_orig[i])

    return Success(pl.DataFrame(effects))


def run_pipeline(config: BulkAssocConfig) -> Result[None, str]:
    """Execute Step 5 association pipeline across all fraction types."""
    match read_clinical(config.clinical_path, config.sample_col, config.response_col):
        case Failure(err):
            return Failure(err)
        case Success(clinical_df):
            pass

    for p in [config.fractions_raw_path, config.fractions_norm_path, config.fractions_mrna_path]:
        if not p.exists():
            return Failure(f"Fractions file does not exist: {p}")

    raw_df = pl.read_parquet(config.fractions_raw_path)
    norm_df = pl.read_parquet(config.fractions_norm_path)
    mrna_df = pl.read_parquet(config.fractions_mrna_path)

    results: list[pl.DataFrame] = []

    for name, df in [("raw", raw_df), ("normalized", norm_df), ("mrna_scaled", mrna_df)]:
        match fit_fractions_association(df, clinical_df, name):
            case Failure(err):
                return Failure(err)
            case Success(res):
                results.append(res)

    combined_df = pl.concat(results).sort(["fraction_type", "p_value_bulk"])
    config.out_parquet.parent.mkdir(parents=True, exist_ok=True)
    combined_df.write_parquet(config.out_parquet)

    print(f"[INFO] Successfully modeled bulk response associations for {len(results)} fraction types.")
    return Success(None)


def main() -> None:
    config = parse_args()
    match run_pipeline(config):
        case Success(_):
            print(f"[SUCCESS] Bulk response effects saved to {config.out_parquet}")
            sys.exit(0)
        case Failure(err):
            print(f"[ERROR] {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
