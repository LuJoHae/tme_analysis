#!/usr/bin/env python3
"""
Step 1: Build Deconvolution Reference from Sade-Feldman Single-Cell Dataset (GSE120575).
Computes linear cluster mean expression profiles and selects top marker genes.
Outputs parquet files to the output directory.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Final

import anndata as ad  # type: ignore
import numpy as np
import numpy.typing as npt
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy.sparse import issparse, spmatrix  # type: ignore


class ReferenceConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    adata_path: Path
    cluster_col: str
    out_dir: Path
    top_markers: int = 10
    min_cluster_size: int = 10


def load_anndata(path: Path) -> Result[ad.AnnData, str]:
    """Pure boundary to safely load AnnData from disk."""
    try:
        if not path.exists():
            return Failure(f"AnnData file not found at: {path}")
        adata = ad.read_h5ad(path)
        return Success(adata)
    except Exception as exc:
        return Failure(f"Failed to read AnnData: {exc}")


def compute_cluster_means(
    expr_matrix: npt.NDArray[np.float64] | spmatrix,
    clusters: tuple[str, ...],
    unique_clusters: tuple[str, ...],
) -> npt.NDArray[np.float64]:
    """Pure function to compute cluster mean expression vectors in linear space."""
    clusters_arr = np.asarray(clusters)
    means_list: list[npt.NDArray[np.float64]] = []

    for cl in unique_clusters:
        mask = clusters_arr == cl
        sub_expr = expr_matrix[mask]
        mean_vec = (
            np.asarray(sub_expr.mean(axis=0)).ravel()
            if issparse(sub_expr)
            else np.mean(sub_expr, axis=0)
        )
        means_list.append(np.asarray(mean_vec, dtype=np.float64))

    return np.vstack(means_list)


def select_top_markers(
    mean_matrix: npt.NDArray[np.float64],
    genes: tuple[str, ...],
    unique_clusters: tuple[str, ...],
    top_n: int,
) -> pl.DataFrame:
    """Find top distinguishing marker genes per cluster via log2 fold-change over rest."""
    records: list[dict[str, object]] = []
    n_clusters = len(unique_clusters)
    eps: Final[float] = 1e-4

    for i, cl in enumerate(unique_clusters):
        cl_mean = mean_matrix[i, :]
        rest_indices = [j for j in range(n_clusters) if j != i]
        rest_mean = (
            np.mean(mean_matrix[rest_indices, :], axis=0)
            if rest_indices
            else np.ones_like(cl_mean) * eps
        )

        log2_fc = np.log2((cl_mean + eps) / (rest_mean + eps))
        # Top gene indices sorted descending
        top_idx = np.argsort(-log2_fc)[:top_n]

        for rank, g_idx in enumerate(top_idx, start=1):
            records.append(
                {
                    "cluster": cl,
                    "gene": str(genes[g_idx]),
                    "rank": rank,
                    "log2fc": float(log2_fc[g_idx]),
                    "cluster_mean": float(cl_mean[g_idx]),
                    "rest_mean": float(rest_mean[g_idx]),
                }
            )

    return pl.DataFrame(records)


def run_reference_pipeline(config: ReferenceConfig) -> Result[Path, str]:
    """Execute reference preparation steps functionally."""
    adata_res = load_anndata(config.adata_path)
    match adata_res:
        case Failure(err):
            return Failure(err)
        case Success(adata):
            pass

    if config.cluster_col not in adata.obs.columns:
        return Failure(
            f"Cluster column '{config.cluster_col}' not in adata.obs (columns: {list(adata.obs.columns)})"
        )

    # Convert gene names to uppercase strings for cross-platform harmonization
    gene_names: tuple[str, ...] = tuple(str(g).upper() for g in adata.var_names)
    cluster_series = adata.obs[config.cluster_col].astype(str)
    unique_clusters: tuple[str, ...] = tuple(sorted(cluster_series.unique().tolist()))
    clusters: tuple[str, ...] = tuple(cluster_series.tolist())

    print(f"Loaded AnnData with {adata.n_obs} cells and {adata.n_vars} genes.")
    print(f"Target clusters ({len(unique_clusters)}): {unique_clusters}")

    # Linearize expression (gse120575_processed.h5ad stores log1p normalized counts)
    print("Linearizing expression matrix (expm1)...")
    if issparse(adata.X):
        expr_linear = adata.X.expm1()
    else:
        expr_linear = np.expm1(adata.X)

    # Compute cluster linear means
    print("Computing cluster mean expression profiles...")
    mean_matrix = compute_cluster_means(expr_linear, clusters, unique_clusters)

    # Filter out genes with zero mean across all clusters
    total_mean = mean_matrix.sum(axis=0)
    valid_genes_mask = total_mean > 0
    filtered_genes: tuple[str, ...] = tuple(
        g for g, v in zip(gene_names, valid_genes_mask) if v
    )
    filtered_mean_matrix = mean_matrix[:, valid_genes_mask]

    print(
        f"Retained {len(filtered_genes)} genes with non-zero expression across clusters."
    )

    # Normalize each cluster profile to simplex (sum to 1) for deconvolution
    row_sums = filtered_mean_matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    phi_matrix = filtered_mean_matrix / row_sums

    # Select top markers
    print(f"Selecting top {config.top_markers} marker genes per cluster...")
    df_markers = select_top_markers(
        filtered_mean_matrix,
        filtered_genes,
        unique_clusters,
        config.top_markers,
    )

    # Prepare Polars DataFrames for export
    config.out_dir.mkdir(parents=True, exist_ok=True)

    # Reference Phi in tidy format: [cluster, gene, expression, weight]
    tidy_records: list[dict[str, object]] = []
    for c_idx, cl in enumerate(unique_clusters):
        for g_idx, g in enumerate(filtered_genes):
            tidy_records.append(
                {
                    "cluster": cl,
                    "gene": g,
                    "linear_mean": float(filtered_mean_matrix[c_idx, g_idx]),
                    "phi_weight": float(phi_matrix[c_idx, g_idx]),
                }
            )

    df_phi_tidy = pl.DataFrame(tidy_records)
    out_tidy = config.out_dir / "reference_phi_tidy.parquet"
    df_phi_tidy.write_parquet(out_tidy)
    print(f"Saved tidy reference Phi to: {out_tidy}")

    # Also save wide matrix: rows = clusters, columns = genes
    wide_dict: dict[str, object] = {"cluster": list(unique_clusters)}
    for g_idx, g in enumerate(filtered_genes):
        wide_dict[g] = [float(phi_matrix[c_idx, g_idx]) for c_idx in range(len(unique_clusters))]
    df_phi_wide = pl.DataFrame(wide_dict)
    out_wide = config.out_dir / "reference_phi.parquet"
    df_phi_wide.write_parquet(out_wide)
    print(f"Saved wide reference Phi to: {out_wide}")

    # Save markers
    out_markers = config.out_dir / "reference_marker_genes.parquet"
    df_markers.write_parquet(out_markers)
    print(f"Saved marker genes to: {out_markers}")

    return Success(out_wide)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 1: Build deconvolution reference from Sade-Feldman single cell dataset."
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
        "--out-dir",
        type=str,
        default="output/sade_feldman_deconv_validation",
        help="Directory to store output parquets",
    )
    parser.add_argument(
        "--top-markers",
        type=int,
        default=10,
        help="Number of top marker genes per cluster to export",
    )
    args = parser.parse_args()

    config = ReferenceConfig(
        adata_path=Path(args.adata),
        cluster_col=args.cluster_col,
        out_dir=Path(args.out_dir),
        top_markers=args.top_markers,
    )

    match run_reference_pipeline(config):
        case Success(out_path):
            print(f"Step 1 completed successfully. Reference saved at: {out_path}")
            sys.exit(0)
        case Failure(err):
            print(f"Step 1 failed with error: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
