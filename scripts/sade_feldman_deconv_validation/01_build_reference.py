#!/usr/bin/env python3
"""
Step 1: Build Deconvolution Reference from Sade-Feldman Single-Cell Dataset (GSE120575).
Computes linear cluster mean expression profiles from true linear TPM (2^x - 1),
filters confounding gene families (ribosomal, mitochondrial, immunoglobulins, pseudogenes, tumor markers),
selects canonical specific marker genes, and outputs a curated signature deconvolution matrix.
"""

from __future__ import annotations

import argparse
import re
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
    tpm_path: Path | None
    out_dir: Path
    resolutions: tuple[float, ...] = (0.5, 1.0, 1.5, 2.0)
    cluster_prefix: str = "celltypist_leiden_"
    top_markers: int = 35
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


def filter_confounding_genes(genes: list[str]) -> tuple[list[int], list[str]]:
    """Filters out confounding gene families:
    - Ribosomal proteins (RPS*, RPL*)
    - Mitochondrial genes (MT-*)
    - Immunoglobulins (IGH*, IGK*, IGL*)
    - Pseudogenes & non-coding RNAs (RP11*, AC00*, RNU*, RNA5S*, LINC*, CTD-*)
    - Melanoma parenchymal antigens (MLANA, PMEL, TYR, DCT, MITF, S100B, MAGEA*)
    """
    exclude_pattern = re.compile(
        r"^(MT-|RP[SL]|RP[0-9]|AC[0-9]|AL[0-9]|AP[0-9]|CTD-|RNU|RNA5S|LINC|IGH|IGK|IGL|MLANA$|PMEL$|TYR$|DCT$|MITF$|S100B$|MAGEA)"
    )
    valid_indices: list[int] = []
    valid_genes: list[str] = []

    for idx, g in enumerate(genes):
        clean_g = str(g).strip().upper()
        if len(clean_g) > 2 and not exclude_pattern.search(clean_g):
            valid_indices.append(idx)
            valid_genes.append(clean_g)

    return valid_indices, valid_genes


def compute_cluster_means(
    expr_matrix: npt.NDArray[np.float64],
    clusters: list[str],
    unique_clusters: list[str],
) -> npt.NDArray[np.float64]:
    """Pure function to compute cluster mean expression vectors in linear space.
    expr_matrix is shaped (n_genes, n_cells).
    Returns (n_clusters, n_genes).
    """
    clusters_arr = np.asarray(clusters)
    means_list: list[npt.NDArray[np.float64]] = []

    for cl in unique_clusters:
        mask = clusters_arr == cl
        if not np.any(mask):
            means_list.append(np.zeros(expr_matrix.shape[0], dtype=np.float64))
        else:
            sub_expr = expr_matrix[:, mask]
            means_list.append(np.asarray(sub_expr.mean(axis=1)).flatten())

    return np.vstack(means_list)


def select_top_markers(
    mean_matrix: npt.NDArray[np.float64],
    genes: list[str],
    unique_clusters: list[str],
    top_n: int,
) -> pl.DataFrame:
    """Select top marker genes per cluster based on one-vs-rest fold change in linear TPM."""
    n_clusters, n_genes = mean_matrix.shape
    records: list[dict[str, object]] = []

    for c_idx, cl in enumerate(unique_clusters):
        cl_mean = mean_matrix[c_idx, :]
        rest_indices = [i for i in range(n_clusters) if i != c_idx]
        rest_mean = (
            np.mean(mean_matrix[rest_indices, :], axis=0)
            if rest_indices
            else np.ones_like(cl_mean) * 1e-6
        )

        # Fold change: (cluster_mean + 1.0) / (rest_mean + 1.0)
        fold_change = (cl_mean + 1.0) / (rest_mean + 1.0)
        log2_fc = np.log2((cl_mean + 1e-3) / (rest_mean + 1e-3))

        # Filter to genes with cluster mean >= 2.0 TPM
        valid_candidate_mask = cl_mean >= 2.0
        candidate_indices = np.where(valid_candidate_mask)[0]

        if len(candidate_indices) > 0:
            sorted_candidates = candidate_indices[np.argsort(-fold_change[candidate_indices])]
            top_idx = sorted_candidates[:top_n]
        else:
            top_idx = np.argsort(-fold_change)[:top_n]

        for rank, g_idx in enumerate(top_idx, start=1):
            records.append(
                {
                    "cluster": cl,
                    "gene": str(genes[g_idx]),
                    "rank": rank,
                    "log2fc": float(log2_fc[g_idx]),
                    "fold_change": float(fold_change[g_idx]),
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

    cells: list[str] = adata.obs_names.tolist()
    print(f"Loaded AnnData with {adata.n_obs} cells.")

    # Determine if raw TPM parquet is provided / exists
    tpm_file = config.tpm_path if config.tpm_path and config.tpm_path.exists() else None
    if tpm_file is None:
        candidate_tpm = config.adata_path.parent / "gse120575_tpm.parquet"
        if candidate_tpm.exists():
            tpm_file = candidate_tpm

    if tpm_file is not None:
        print(f"Loading raw TPM parquet from {tpm_file}...")
        tpm_df = pl.read_parquet(tpm_file)
        raw_genes = [str(g).upper() for g in tpm_df["gene"].to_list()]

        # Filter confounding gene families
        valid_gene_idx, filtered_genes = filter_confounding_genes(raw_genes)
        print(f"Retained {len(filtered_genes)} protein-coding genes after confounding gene exclusion.")

        # Intersect cells
        common_cells = [c for c in cells if c in tpm_df.columns]
        cell_mask = [c in set(common_cells) for c in cells]

        print(f"Aligning {len(common_cells)} single cells with TPM data...")
        sub_tpm = tpm_df.select(common_cells).to_numpy()[valid_gene_idx, :]  # (genes, cells)

        # Invert log2(TPM + 1) to linear TPM: 2^x - 1
        print("Inverting log2(TPM + 1) to linear TPM scale (2^x - 1)...")
        lin_matrix = np.power(2.0, sub_tpm) - 1.0
    else:
        print("Warning: TPM parquet not found. Falling back to linearizing adata.X...")
        raw_genes = [str(g).upper() for g in adata.var_names]
        valid_gene_idx, filtered_genes = filter_confounding_genes(raw_genes)

        if issparse(adata.X):
            lin_matrix = adata.X[:, valid_gene_idx].expm1().toarray().T  # (genes, cells)
        else:
            lin_matrix = np.expm1(adata.X[:, valid_gene_idx]).T
        common_cells = cells
        cell_mask = [True] * len(cells)

    config.out_dir.mkdir(parents=True, exist_ok=True)
    resolution_records: list[dict[str, object]] = []
    primary_out_path = config.out_dir / "reference_phi.parquet"

    for res in config.resolutions:
        cluster_col = f"{config.cluster_prefix}{res}"
        if cluster_col not in adata.obs.columns:
            # Fallback to alternate naming if present
            alt_col = f"leiden_{res}"
            if alt_col in adata.obs.columns:
                cluster_col = alt_col
            else:
                print(f"Warning: Cluster column '{cluster_col}' not found in AnnData. Skipping res={res}...")
                continue

        cluster_series = adata.obs[cluster_col].astype(str)
        unique_clusters: list[str] = sorted(cluster_series.unique().tolist())
        clusters: list[str] = cluster_series.tolist()
        aligned_clusters = [cl for cl, m in zip(clusters, cell_mask) if m]

        print(f"\n--- Processing Sade-Feldman Resolution res={res} ({len(unique_clusters)} clusters) ---")
        mean_matrix = compute_cluster_means(lin_matrix, aligned_clusters, unique_clusters)

        # Filter out genes with zero mean across all clusters
        total_mean = mean_matrix.sum(axis=0)
        nonzero_mask = total_mean > 0
        final_genes = [g for g, v in zip(filtered_genes, nonzero_mask) if v]
        final_mean_matrix = mean_matrix[:, nonzero_mask]

        # Select top marker genes per cluster
        df_markers = select_top_markers(
            final_mean_matrix,
            final_genes,
            unique_clusters,
            config.top_markers,
        )

        signature_genes = sorted(list(set(df_markers["gene"].to_list())))
        sig_gene_indices = [final_genes.index(g) for g in signature_genes]
        sig_mean_matrix = final_mean_matrix[:, sig_gene_indices]

        # Normalize signature matrix to simplex
        row_sums = sig_mean_matrix.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        phi_signature = sig_mean_matrix / row_sums

        # Condition number of signature matrix (measure of collinearity)
        cond_num = float(np.linalg.cond(phi_signature))
        print(f"Resolution {res}: {len(unique_clusters)} clusters, {len(signature_genes)} signature genes, condition number = {cond_num:.2f}")

        resolution_records.append(
            {
                "reference_type": "Sade-Feldman",
                "resolution": float(res),
                "cluster_col": cluster_col,
                "n_clusters": len(unique_clusters),
                "n_signature_genes": len(signature_genes),
                "condition_number": cond_num,
            }
        )

        # 1. Save curated signature Phi matrix for this resolution
        wide_dict: dict[str, object] = {"cluster": list(unique_clusters)}
        for g_idx, g in enumerate(signature_genes):
            wide_dict[g] = [float(phi_signature[c_idx, g_idx]) for c_idx in range(len(unique_clusters))]
        df_phi_wide = pl.DataFrame(wide_dict)

        out_res_wide = config.out_dir / f"reference_phi_res{res}.parquet"
        df_phi_wide.write_parquet(out_res_wide)
        print(f"Saved signature Phi (res={res}) to: {out_res_wide}")

        # 2. Save markers table for this resolution
        out_res_markers = config.out_dir / f"reference_marker_genes_res{res}.parquet"
        df_markers.write_parquet(out_res_markers)

        # 3. Save tidy Phi for this resolution
        tidy_records: list[dict[str, object]] = []
        for c_idx, cl in enumerate(unique_clusters):
            for g_idx, g in enumerate(signature_genes):
                tidy_records.append(
                    {
                        "cluster": cl,
                        "gene": g,
                        "linear_mean": float(sig_mean_matrix[c_idx, g_idx]),
                        "phi_weight": float(phi_signature[c_idx, g_idx]),
                    }
                )
        df_phi_tidy = pl.DataFrame(tidy_records)
        df_phi_tidy.write_parquet(config.out_dir / f"reference_phi_tidy_res{res}.parquet")

        # Save default backward-compatible files for res == 0.5
        if abs(res - 0.5) < 1e-4:
            df_phi_wide.write_parquet(primary_out_path)
            df_markers.write_parquet(config.out_dir / "reference_marker_genes.parquet")
            df_phi_tidy.write_parquet(config.out_dir / "reference_phi_tidy.parquet")
            print("Saved backward-compatible primary reference_phi.parquet (res=0.5).")

    # Save resolution benchmarking metrics table
    df_metrics = pl.DataFrame(resolution_records)
    out_metrics = config.out_dir / "reference_resolution_metrics_sf.parquet"
    df_metrics.write_parquet(out_metrics)
    print(f"\nSaved Sade-Feldman multi-resolution metrics to: {out_metrics}")

    return Success(primary_out_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 1: Build multi-resolution deconvolution references from Sade-Feldman single cell dataset."
    )
    parser.add_argument(
        "--adata",
        type=str,
        default="/storage/halu/data/GSE120575/gse120575_processed.h5ad",
        help="Path to gse120575_processed.h5ad",
    )
    parser.add_argument(
        "--tpm",
        type=str,
        default="/storage/halu/data/GSE120575/gse120575_tpm.parquet",
        help="Path to gse120575_tpm.parquet",
    )
    parser.add_argument(
        "--resolutions",
        type=str,
        default="0.5,1.0,1.5,2.0",
        help="Comma-separated clustering resolutions to evaluate",
    )
    parser.add_argument(
        "--cluster-prefix",
        type=str,
        default="celltypist_leiden_",
        help="Prefix for cluster columns in adata.obs (e.g. celltypist_leiden_)",
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
        default=35,
        help="Number of top marker genes per cluster to include in signature",
    )
    args = parser.parse_args()

    res_list = tuple(float(r.strip()) for r in args.resolutions.split(",") if r.strip())

    config = ReferenceConfig(
        adata_path=Path(args.adata),
        tpm_path=Path(args.tpm) if args.tpm else None,
        resolutions=res_list,
        cluster_prefix=args.cluster_prefix,
        out_dir=Path(args.out_dir),
        top_markers=args.top_markers,
    )

    match run_reference_pipeline(config):
        case Success(out_path):
            print(f"Step 1 completed successfully. Multi-resolution references saved at: {config.out_dir}")
            sys.exit(0)
        case Failure(err):
            print(f"Step 1 failed with error: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
