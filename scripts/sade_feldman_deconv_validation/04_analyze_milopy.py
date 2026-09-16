#!/usr/bin/env python3
"""
Step 4: Milopy Differential Abundance (DA) Analysis on Sade-Feldman Dataset (GSE120575).
Extracts neighborhood-level DA testing (~response: Responder vs Non-Responder),
projects neighborhood logFC onto individual cells, and aggregates metrics to reference cell states.
Outputs milopy_cell_state_da.parquet, milopy_cell_level_scores.parquet, and milopy_nhoods_results.parquet.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Final

import anndata as ad  # type: ignore
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import sparse as sp  # type: ignore
from scipy import stats  # type: ignore


class MiloConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    adata_path: Path
    cluster_col: str
    milo_dir: Path
    condition: str
    out_dir: Path


def load_or_run_milopy(
    adata_path: Path,
    milo_dir: Path,
    condition: str,
) -> Result[tuple[pd.DataFrame, sp.spmatrix], str]:
    """Load existing milopy results and nhoods matrix, or compute if absent."""
    csv_file = milo_dir / f"milopy_results_{condition}.csv"
    npz_file = milo_dir / f"milopy_nhoods_{condition}.npz"

    if csv_file.exists() and npz_file.exists():
        print(f"Loading existing Milo results from {csv_file} and {npz_file}...")
        try:
            df_milo = pd.read_csv(csv_file, index_col=0)
            mat_nhoods = sp.load_npz(npz_file)
            return Success((df_milo, mat_nhoods))
        except Exception as exc:
            print(f"Failed to load cached Milo files ({exc}), running fresh milopy...")

    # Run milopy
    print("Running milopy differential abundance analysis...")
    try:
        import milopy  # type: ignore

        adata = ad.read_h5ad(adata_path)
        if condition != "Combined":
            adata = adata[adata.obs["treatment_status"] == condition].copy()

        # Set response column
        if "response" not in adata.obs.columns:
            adata.obs["response"] = adata.obs["characteristics: response"]

        n_obs = adata.n_obs
        k_val = min(50, max(5, n_obs // 10))
        milopy.core.make_nhoods(adata, prop=0.2, k=k_val, d=40)
        milopy.core.count_cells(adata, sample_col="melanoma-sample")

        design_df = (
            adata.obs[["melanoma-sample", "response"]]
            .drop_duplicates()
            .set_index("melanoma-sample")
        )
        milopy.core.test_nhoods(adata, design="~response", design_df=design_df)

        df_res = adata.uns["nhood_test_results"]
        nhoods_mat = adata.obsm["nhoods"]

        # Cache to milo_dir
        milo_dir.mkdir(parents=True, exist_ok=True)
        df_res.to_csv(csv_file)
        sp.save_npz(npz_file, nhoods_mat)

        return Success((df_res, nhoods_mat))
    except Exception as exc:
        return Failure(f"Milopy analysis failed: {exc}")


def project_nhoods_to_cells(
    nhoods_mat: sp.spmatrix,
    df_milo: pd.DataFrame,
) -> np.ndarray:
    """Compute average neighborhood logFC for each cell."""
    n_cells, n_nhoods = nhoods_mat.shape
    logfc_vec = df_milo["logFC"].to_numpy().astype(np.float64)

    # nhoods_mat is (cells x nhoods)
    # Sum of logFC across all nhoods a cell belongs to
    cell_logfc_sum = nhoods_mat.dot(logfc_vec)
    # Degree of each cell (number of nhoods it belongs to)
    cell_degrees = np.asarray(nhoods_mat.sum(axis=1)).ravel()

    # Avoid div by zero
    valid_mask = cell_degrees > 0
    cell_mean_logfc = np.zeros(n_cells, dtype=np.float64)
    cell_mean_logfc[valid_mask] = cell_logfc_sum[valid_mask] / cell_degrees[valid_mask]

    return cell_mean_logfc


def run_milopy_aggregation(config: MiloConfig) -> Result[Path, str]:
    """Execute milopy projection and cell state level aggregation."""
    if not config.adata_path.exists():
        return Failure(f"AnnData not found: {config.adata_path}")

    # Load AnnData obs and UMAP
    adata = ad.read_h5ad(config.adata_path, backed="r")
    if config.cluster_col not in adata.obs.columns:
        return Failure(f"Cluster column '{config.cluster_col}' not in AnnData.")

    cell_clusters = adata.obs[config.cluster_col].astype(str).to_numpy()
    cell_ids = adata.obs_names.astype(str).to_numpy()
    umap_coords = (
        adata.obsm["X_umap"]
        if "X_umap" in adata.obsm
        else np.zeros((adata.n_obs, 2))
    )
    samples = adata.obs.get("melanoma-sample", pd.Series([""] * adata.n_obs)).astype(str).to_numpy()
    responses = adata.obs.get("characteristics: response", pd.Series([""] * adata.n_obs)).astype(str).to_numpy()
    treatment = adata.obs.get("treatment_status", pd.Series([""] * adata.n_obs)).astype(str).to_numpy()

    # Load or run Milo
    milo_res = load_or_run_milopy(config.adata_path, config.milo_dir, config.condition)
    match milo_res:
        case Failure(err):
            return Failure(err)
        case Success((df_milo, nhoods_mat)):
            pass

    print(f"Milopy results shape: {df_milo.shape}, nhoods matrix shape: {nhoods_mat.shape}")

    # Project to cells
    cell_logfc = project_nhoods_to_cells(nhoods_mat, df_milo)

    # Save cell-level parquet
    config.out_dir.mkdir(parents=True, exist_ok=True)
    df_cells = pl.DataFrame(
        {
            "cell_id": cell_ids,
            "umap_1": umap_coords[:, 0].astype(float),
            "umap_2": umap_coords[:, 1].astype(float),
            "cell_state": cell_clusters,
            "milo_logfc": cell_logfc,
            "sample_id": samples,
            "response": responses,
            "treatment_status": treatment,
        }
    )
    out_cells = config.out_dir / "milopy_cell_level_scores.parquet"
    df_cells.write_parquet(out_cells)
    print(f"Saved cell-level Milo scores to: {out_cells}")

    # Save neighborhood-level parquet
    nhood_recs: list[dict[str, object]] = []
    for idx, row in df_milo.iterrows():
        nhood_recs.append(
            {
                "nhood_id": int(idx) if isinstance(idx, (int, np.integer)) else str(idx),
                "logfc": float(row.get("logFC", 0.0)),
                "p_value": float(row.get("PValue", 1.0)),
                "fdr": float(row.get("FDR", 1.0)),
                "nhood_group": (
                    int(row["NhoodGroup"])
                    if "NhoodGroup" in row and pd.notna(row["NhoodGroup"])
                    else -1
                ),
            }
        )
    df_nhoods = pl.DataFrame(nhood_recs)
    out_nhoods = config.out_dir / "milopy_nhoods_results.parquet"
    df_nhoods.write_parquet(out_nhoods)
    print(f"Saved neighborhood-level Milo table to: {out_nhoods}")

    # Aggregate to cell states
    unique_clusters = sorted(list(set(cell_clusters)))
    state_records: list[dict[str, object]] = []

    for cl in unique_clusters:
        mask = cell_clusters == cl
        cl_logfc = cell_logfc[mask]
        n_cells = len(cl_logfc)

        mean_lfc = float(np.mean(cl_logfc))
        median_lfc = float(np.median(cl_logfc))
        std_lfc = float(np.std(cl_logfc))
        iqr_lfc = float(np.percentile(cl_logfc, 75) - np.percentile(cl_logfc, 25))

        # One-sample Wilcoxon test vs 0
        try:
            if np.all(cl_logfc == 0):
                w_pval = 1.0
            else:
                _, w_pval = stats.wilcoxon(cl_logfc)
                w_pval = float(w_pval)
        except Exception:
            w_pval = 1.0

        state_records.append(
            {
                "cell_state": cl,
                "n_cells": n_cells,
                "milo_mean_logfc": mean_lfc,
                "milo_median_logfc": median_lfc,
                "milo_std_logfc": std_lfc,
                "milo_iqr_logfc": iqr_lfc,
                "milo_wilcoxon_pval": w_pval,
                "pct_positive_cells": float(np.mean(cl_logfc > 0)),
                "pct_negative_cells": float(np.mean(cl_logfc < 0)),
            }
        )

    df_states = pl.DataFrame(state_records)
    out_states = config.out_dir / "milopy_cell_state_da.parquet"
    df_states.write_parquet(out_states)
    print(f"Saved cell-state aggregated Milo table to: {out_states}")

    return Success(out_states)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 4: Milopy differential abundance analysis and cell-state aggregation."
    )
    parser.add_argument(
        "--adata",
        type=str,
        default="/storage/halu/data/GSE120575/gse120575_processed.h5ad",
        help="Path to gse120575_processed.h5ad",
    )
    parser.add_argument(
        "--cluster-col",
        type=str,
        default="celltypist_leiden_0.5",
        help="Cluster column in adata.obs",
    )
    parser.add_argument(
        "--milo-dir",
        type=str,
        default="/storage/halu/data/output/whole_dataset",
        help="Directory containing or to store milopy results CSV/NPZ",
    )
    parser.add_argument(
        "--condition",
        type=str,
        default="Combined",
        choices=["Combined", "Pre", "Post"],
        help="Condition subset to evaluate",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="output/sade_feldman_deconv_validation",
        help="Directory to save output parquets",
    )
    args = parser.parse_args()

    config = MiloConfig(
        adata_path=Path(args.adata),
        cluster_col=args.cluster_col,
        milo_dir=Path(args.milo_dir),
        condition=args.condition,
        out_dir=Path(args.out_dir),
    )

    match run_milopy_aggregation(config):
        case Success(out_path):
            print(f"Step 4 finished successfully: {out_path}")
            sys.exit(0)
        case Failure(err):
            print(f"Step 4 failed with error: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
