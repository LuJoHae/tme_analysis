#!/usr/bin/env python3
"""
Step 4: Milopy Differential Abundance (DA) Analysis on Sade-Feldman Dataset (GSE120575).
Extracts neighborhood-level DA testing (~response: Responder vs Non-Responder)
stratified across Combined (full cohort), Pre-treatment, and Post-treatment samples.
Projects neighborhood logFC onto cells and aggregates metrics to reference cell states.
Outputs parquets for Combined, Pre, and Post conditions.
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
    out_dir: Path


def load_or_run_milopy_single(
    adata_full: ad.AnnData,
    milo_dir: Path,
    condition: str,
) -> Result[tuple[pd.DataFrame, sp.spmatrix, list[str]], str]:
    """Load existing milopy results and nhoods matrix for a condition, or compute."""
    csv_file = milo_dir / f"milopy_results_{condition}.csv"
    npz_file = milo_dir / f"milopy_nhoods_{condition}.npz"

    # Subset cell IDs
    if condition == "Combined":
        subset = adata_full
    else:
        subset = adata_full[adata_full.obs["treatment_status"] == condition]

    subset_cell_ids = list(subset.obs_names.astype(str))

    if csv_file.exists() and npz_file.exists():
        print(f"Loading existing Milo results for '{condition}' from {csv_file.name} and {npz_file.name}...")
        try:
            df_milo = pd.read_csv(csv_file, index_col=0)
            mat_nhoods = sp.load_npz(npz_file)
            return Success((df_milo, mat_nhoods, subset_cell_ids))
        except Exception as exc:
            print(f"Failed to load cached Milo files for {condition} ({exc}), running fresh...")

    print(f"Running milopy differential abundance analysis for condition: {condition}...")
    try:
        import milopy  # type: ignore

        sub_adata = subset.copy()
        if "response" not in sub_adata.obs.columns:
            sub_adata.obs["response"] = sub_adata.obs["characteristics: response"]

        n_obs = sub_adata.n_obs
        k_val = min(50, max(5, n_obs // 10))
        milopy.core.make_nhoods(sub_adata, prop=0.2, k=k_val, d=40)
        milopy.core.count_cells(sub_adata, sample_col="melanoma-sample")

        design_df = (
            sub_adata.obs[["melanoma-sample", "response"]]
            .drop_duplicates()
            .set_index("melanoma-sample")
        )
        milopy.core.test_nhoods(sub_adata, design="~response", design_df=design_df)

        df_res = sub_adata.uns["nhood_test_results"]
        nhoods_mat = sub_adata.obsm["nhoods"]

        milo_dir.mkdir(parents=True, exist_ok=True)
        df_res.to_csv(csv_file)
        sp.save_npz(npz_file, nhoods_mat)

        return Success((df_res, nhoods_mat, subset_cell_ids))
    except Exception as exc:
        return Failure(f"Milopy analysis failed for {condition}: {exc}")


def project_nhoods_to_cells(
    nhoods_mat: sp.spmatrix,
    df_milo: pd.DataFrame,
) -> np.ndarray:
    """Compute average neighborhood logFC for each cell in the subset."""
    n_cells, n_nhoods = nhoods_mat.shape
    logfc_vec = df_milo["logFC"].to_numpy().astype(np.float64)

    cell_logfc_sum = nhoods_mat.dot(logfc_vec)
    cell_degrees = np.asarray(nhoods_mat.sum(axis=1)).ravel()

    valid_mask = cell_degrees > 0
    cell_mean_logfc = np.zeros(n_cells, dtype=np.float64)
    cell_mean_logfc[valid_mask] = cell_logfc_sum[valid_mask] / cell_degrees[valid_mask]

    return cell_mean_logfc


def run_milopy_multi_condition_pipeline(config: MiloConfig) -> Result[Path, str]:
    """Run Milo DA across Combined, Pre, and Post cohorts and export aggregated summaries."""
    if not config.adata_path.exists():
        return Failure(f"AnnData not found: {config.adata_path}")

    adata = ad.read_h5ad(config.adata_path)
    if config.cluster_col not in adata.obs.columns:
        return Failure(f"Cluster column '{config.cluster_col}' not in AnnData.")

    all_cell_ids = list(adata.obs_names.astype(str))
    n_total_cells = len(all_cell_ids)
    cell_id_to_idx = {cid: idx for idx, cid in enumerate(all_cell_ids)}

    cell_clusters = adata.obs[config.cluster_col].astype(str).to_numpy()
    umap_coords = (
        adata.obsm["X_umap"]
        if "X_umap" in adata.obsm
        else np.zeros((n_total_cells, 2))
    )
    samples = adata.obs.get("melanoma-sample", pd.Series([""] * n_total_cells)).astype(str).to_numpy()
    responses = adata.obs.get("characteristics: response", pd.Series([""] * n_total_cells)).astype(str).to_numpy()
    treatment = adata.obs.get("treatment_status", pd.Series([""] * n_total_cells)).astype(str).to_numpy()

    config.out_dir.mkdir(parents=True, exist_ok=True)
    conditions = ["Combined", "Pre", "Post"]
    condition_cell_logfcs: dict[str, np.ndarray] = {}
    all_state_dfs: list[pl.DataFrame] = []

    for cname in conditions:
        milo_res = load_or_run_milopy_single(adata, config.milo_dir, cname)
        match milo_res:
            case Failure(err):
                print(f"Warning: Milo failed for {cname}: {err}")
                continue
            case Success((df_milo, nhoods_mat, sub_cell_ids)):
                pass

        print(f"[{cname}] Projecting {df_milo.shape[0]} nhoods onto {len(sub_cell_ids)} cells...")
        sub_logfc = project_nhoods_to_cells(nhoods_mat, df_milo)

        # Full array aligned to master cell IDs (fill NaN for cells outside condition)
        full_logfc = np.full(n_total_cells, np.nan, dtype=np.float64)
        for s_idx, cid in enumerate(sub_cell_ids):
            full_logfc[cell_id_to_idx[cid]] = sub_logfc[s_idx]
        condition_cell_logfcs[cname] = full_logfc

        # Save neighborhood-level parquet
        nhood_recs: list[dict[str, object]] = []
        for idx, row in df_milo.iterrows():
            nhood_recs.append(
                {
                    "condition": cname,
                    "nhood_id": int(idx) if isinstance(idx, (int, np.integer)) else str(idx),
                    "logfc": float(row.get("logFC", 0.0)),
                    "p_value": float(row.get("PValue", 1.0)),
                    "fdr": float(row.get("FDR", 1.0)),
                }
            )
        df_nhoods = pl.DataFrame(nhood_recs)
        out_nhoods = config.out_dir / f"milopy_nhoods_results_{cname}.parquet"
        df_nhoods.write_parquet(out_nhoods)
        if cname == "Combined":
            df_nhoods.write_parquet(config.out_dir / "milopy_nhoods_results.parquet")

        # Aggregate by cell state for this condition
        unique_clusters = sorted(list(set(cell_clusters)))
        state_recs: list[dict[str, object]] = []

        for cl in unique_clusters:
            mask = (cell_clusters == cl) & (~np.isnan(full_logfc))
            vals = full_logfc[mask]
            n_sub_cells = len(vals)

            if n_sub_cells > 0:
                m_mean = float(np.mean(vals))
                m_median = float(np.median(vals))
                m_std = float(np.std(vals))
                m_iqr = float(np.percentile(vals, 75) - np.percentile(vals, 25))
                try:
                    w_pval = float(stats.wilcoxon(vals)[1]) if not np.all(vals == 0) else 1.0
                except Exception:
                    w_pval = 1.0
                pct_pos = float(np.mean(vals > 0))
                pct_neg = float(np.mean(vals < 0))
            else:
                m_mean, m_median, m_std, m_iqr, w_pval, pct_pos, pct_neg = 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0

            state_recs.append(
                {
                    "condition": cname,
                    "cell_state": cl,
                    "n_cells": n_sub_cells,
                    "milo_mean_logfc": m_mean,
                    "milo_median_logfc": m_median,
                    "milo_std_logfc": m_std,
                    "milo_iqr_logfc": m_iqr,
                    "milo_wilcoxon_pval": w_pval,
                    "pct_positive_cells": pct_pos,
                    "pct_negative_cells": pct_neg,
                }
            )

        df_states = pl.DataFrame(state_recs)
        all_state_dfs.append(df_states)
        df_states.write_parquet(config.out_dir / f"milopy_cell_state_da_{cname}.parquet")
        if cname == "Combined":
            df_states.write_parquet(config.out_dir / "milopy_cell_state_da.parquet")

    # Combine all state aggregations into one master table
    df_all_states = pl.concat(all_state_dfs, how="vertical")
    out_all_states = config.out_dir / "milopy_cell_state_da_all.parquet"
    df_all_states.write_parquet(out_all_states)
    print(f"Saved master cell-state Milo table to: {out_all_states}")

    # Build master cell-level scores table
    df_cells = pl.DataFrame(
        {
            "cell_id": all_cell_ids,
            "umap_1": umap_coords[:, 0].astype(float),
            "umap_2": umap_coords[:, 1].astype(float),
            "cell_state": cell_clusters,
            "sample_id": samples,
            "response": responses,
            "treatment_status": treatment,
            "milo_logfc": condition_cell_logfcs.get("Combined", np.zeros(n_total_cells)),
            "milo_logfc_combined": condition_cell_logfcs.get("Combined", np.zeros(n_total_cells)),
            "milo_logfc_pre": condition_cell_logfcs.get("Pre", np.zeros(n_total_cells)),
            "milo_logfc_post": condition_cell_logfcs.get("Post", np.zeros(n_total_cells)),
        }
    )
    out_cells = config.out_dir / "milopy_cell_level_scores.parquet"
    df_cells.write_parquet(out_cells)
    print(f"Saved master cell-level scores table to: {out_cells}")

    return Success(out_all_states)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 4: Milopy DA analysis across Combined, Pre, and Post conditions."
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
        out_dir=Path(args.out_dir),
    )

    match run_milopy_multi_condition_pipeline(config):
        case Success(out_path):
            print(f"Step 4 multi-condition analysis completed successfully: {out_path}")
            sys.exit(0)
        case Failure(err):
            print(f"Step 4 failed with error: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
