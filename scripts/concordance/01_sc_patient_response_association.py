#!/usr/bin/env python3
"""
Step 1: Patient-Level Single-Cell Response Association.
Computes patient-level cell proportions, applies Centered Log-Ratio (CLR) transformation,
and fits patient-level models against response labels to eliminate pseudoreplication.
Strict functional Python with returns and Polars.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence, assert_never

import anndata as ad  # type: ignore
import numpy as np
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
import polars as pl
from scipy import stats  # type: ignore


class AssociationConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    adata_path: Path
    patient_col: str
    response_col: str
    cluster_col: str
    out_parquet: Path
    pseudocount: float = 1e-5


def parse_args() -> AssociationConfig:
    parser = argparse.ArgumentParser(
        description="Compute patient-level response association for single-cell data."
    )
    parser.add_argument(
        "--adata",
        type=Path,
        required=True,
        help="Path to single-cell AnnData .h5ad",
    )
    parser.add_argument(
        "--patient-col",
        type=str,
        default="patient",
        help="Obs column name for patient identifier",
    )
    parser.add_argument(
        "--response-col",
        type=str,
        default="response",
        help="Obs column name for response label (e.g. Responder/Non-responder or R/NR)",
    )
    parser.add_argument(
        "--cluster-col",
        type=str,
        default="celltypist_leiden_0.5",
        help="Obs column name for cell type or state",
    )
    parser.add_argument(
        "--out-parquet",
        type=Path,
        required=True,
        help="Path to output effects parquet file",
    )
    parser.add_argument(
        "--pseudocount",
        type=float,
        default=1e-5,
        help="Pseudocount for CLR transformation",
    )
    args = parser.parse_args()
    return AssociationConfig(
        adata_path=args.adata,
        patient_col=args.patient_col,
        response_col=args.response_col,
        cluster_col=args.cluster_col,
        out_parquet=args.out_parquet,
        pseudocount=args.pseudocount,
    )


def extract_patient_contingency(
    adata: ad.AnnData,
    patient_col: str,
    response_col: str,
    cluster_col: str,
) -> Result[tuple[pl.DataFrame, pl.DataFrame], str]:
    """Pure function extracting patient cell-counts contingency and patient response metadata."""
    obs_df = pl.from_pandas(adata.obs.reset_index())

    missing_cols = [
        col for col in [patient_col, response_col, cluster_col] if col not in obs_df.columns
    ]
    if missing_cols:
        return Failure(f"Columns missing from AnnData.obs: {missing_cols}")

    # Patient metadata (id and binary response)
    patient_meta = (
        obs_df.select([patient_col, response_col])
        .unique()
        .filter(pl.col(response_col).is_not_null())
    )

    # Contingency: patient x cluster cell counts
    counts_df = (
        obs_df.group_by([patient_col, cluster_col])
        .agg(pl.len().alias("count"))
        .pivot(
            on=cluster_col,
            index=patient_col,
            values="count",
        )
        .fill_null(0)
    )

    return Success((counts_df, patient_meta))


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


def compute_clr_matrix(
    counts: np.ndarray,
    pseudocount: float,
) -> np.ndarray:
    """Pure computation of centered log-ratio (CLR) transformed proportions."""
    row_sums = counts.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0, 1.0, row_sums)
    props = (counts + pseudocount) / (row_sums + pseudocount * counts.shape[1])
    log_props = np.log(props)
    geometric_mean_log = log_props.mean(axis=1, keepdims=True)
    return log_props - geometric_mean_log


def fit_patient_association(
    counts_df: pl.DataFrame,
    patient_meta: pl.DataFrame,
    patient_col: str,
    response_col: str,
    pseudocount: float,
) -> Result[pl.DataFrame, str]:
    """Fit patient-level association for each cell state using CLR and logistic regression."""
    # Align counts and response
    joined = counts_df.join(patient_meta, on=patient_col, how="inner")
    if joined.is_empty():
        return Failure("No matching patients found between counts and response metadata.")

    # Filter patients with valid binary responses
    binary_resp = [binarize_response(str(v)) for v in joined[response_col].to_list()]
    valid_mask = [opt != Nothing for opt in binary_resp]

    if sum(valid_mask) < 4:
        return Failure(f"Insufficient number of patients with binary response labels: {sum(valid_mask)}")

    filtered_df = joined.filter(pl.Series(valid_mask))
    y = np.array([binarize_response(str(v)).unwrap() for v in filtered_df[response_col].to_list()], dtype=np.float64)

    cell_state_cols = [c for c in counts_df.columns if c != patient_col]
    counts_mat = filtered_df.select(cell_state_cols).to_numpy().astype(np.float64)

    # Compute CLR
    clr_mat = compute_clr_matrix(counts_mat, pseudocount)
    mean_props = counts_mat.sum(axis=0) / counts_mat.sum()

    effects: list[dict[str, object]] = []

    for idx, state_name in enumerate(cell_state_cols):
        x = clr_mat[:, idx]
        x_std = (x - np.mean(x)) / (np.std(x) + 1e-9)

        r_val, p_val = stats.pointbiserialr(y, x_std)
        r_clipped = np.clip(r_val, -0.999, 0.999)
        beta = float(2.0 * r_clipped / np.sqrt(1.0 - r_clipped**2 + 1e-12))
        se = float(np.sqrt(4.0 / max(len(y) - 2, 1)))
        z_score = beta / (se + 1e-12)

        effects.append({
            "cell_state": state_name,
            "beta_sc": beta,
            "se_sc": se,
            "z_score_sc": z_score,
            "p_value_sc": float(p_val),
            "mean_proportion": float(mean_props[idx]),
            "n_patients_eval": int(len(y)),
        })

    # Benjamini-Hochberg FDR
    p_vals = np.array([e["p_value_sc"] for e in effects])
    order = np.argsort(p_vals)
    ranked_p = p_vals[order]
    n_tests = len(effects)
    fdrs = np.minimum(1.0, ranked_p * n_tests / (np.arange(1, n_tests + 1)))
    fdrs_mono = np.minimum.accumulate(fdrs[::-1])[::-1]
    fdrs_orig = np.empty_like(fdrs_mono)
    fdrs_orig[order] = fdrs_mono

    for i, e in enumerate(effects):
        e["fdr_sc"] = float(fdrs_orig[i])

    res_df = pl.DataFrame(effects).sort("p_value_sc")
    return Success(res_df)


def run_pipeline(config: AssociationConfig) -> Result[pl.DataFrame, str]:
    """Execute Step 1 pure pipeline."""
    if not config.adata_path.exists():
        return Failure(f"AnnData file does not exist: {config.adata_path}")

    try:
        adata = ad.read_h5ad(config.adata_path, backed=False)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to read AnnData: {exc}")

    match extract_patient_contingency(adata, config.patient_col, config.response_col, config.cluster_col):
        case Failure(err):
            return Failure(err)
        case Success((counts_df, patient_meta)):
            pass

    match fit_patient_association(
        counts_df,
        patient_meta,
        config.patient_col,
        config.response_col,
        config.pseudocount,
    ):
        case Failure(err):
            return Failure(err)
        case Success(results_df):
            config.out_parquet.parent.mkdir(parents=True, exist_ok=True)
            results_df.write_parquet(config.out_parquet)
            return Success(results_df)


def main() -> None:
    config = parse_args()
    match run_pipeline(config):
        case Success(df):
            print(f"[SUCCESS] Saved single-cell patient response effects ({df.height} states) to {config.out_parquet}")
            sys.exit(0)
        case Failure(err):
            print(f"[ERROR] {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
