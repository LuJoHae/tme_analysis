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
    cohort: str = "GSE120575"
    adata_path: Path | None = None
    tpm_path: Path | None = None
    out_dir: Path
    resolutions: tuple[float, ...] = (0.5, 1.0, 1.5, 2.0)
    cluster_prefix: str = "celltypist_leiden_"
    top_markers: int = 35
    min_cluster_size: int = 10


def load_reference_adata(config: ReferenceConfig) -> Result[ad.AnnData, str]:
    """Pure boundary to load AnnData strictly via unified tme_datasets.load_dataset()."""
    from tme_datasets import load_dataset

    return load_dataset(config.cohort, auto_download=False, normalize_ensembl=False)


def _extract_adata_gene_symbols(adata: ad.AnnData) -> list[str]:
    """Extract uppercase HGNC gene symbols from AnnData .var['gene_name'] or .var_names."""
    if "gene_name" in adata.var.columns:
        return [str(g).strip().upper() for g in adata.var["gene_name"].tolist()]
    return [str(g).strip().upper() for g in adata.var_names.tolist()]


def _ensure_resolutions_clustered(
    adata: ad.AnnData,
    resolutions: tuple[float, ...],
    cluster_prefix: str,
) -> ad.AnnData:
    """Compute Leiden clusters on the fly if requested resolution columns are absent from adata.obs."""
    missing = [
        res
        for res in resolutions
        if f"{cluster_prefix}{res}" not in adata.obs.columns and f"leiden_{res}" not in adata.obs.columns
    ]
    if not missing:
        return adata

    import scanpy as sc  # type: ignore
    from scipy.sparse import csr_matrix  # type: ignore

    print(f"Computing Scanpy PCA/neighbors and Leiden clusters for missing resolutions: {missing}...")
    work = ad.AnnData(
        X=adata.layers["log1p"].copy()
        if "log1p" in adata.layers
        else (adata.X.copy() if issparse(adata.X) else csr_matrix(adata.X)),
        obs=adata.obs.copy(),
    )
    max_v = float(work.X.data.max() if issparse(work.X) and work.X.nnz > 0 else np.max(work.X))
    if max_v > 50.0:
        sc.pp.normalize_total(work, target_sum=1e4)
        sc.pp.log1p(work)
    sc.pp.highly_variable_genes(work, n_top_genes=min(2000, work.n_vars), subset=False)
    sc.tl.pca(work, svd_solver="arpack", random_state=42)
    sc.pp.neighbors(work, n_neighbors=15, n_pcs=min(30, work.obsm["X_pca"].shape[1]), random_state=42)
    out = adata.copy()
    for res in missing:
        col = f"{cluster_prefix}{res}"
        sc.tl.leiden(work, resolution=res, key_added=col, flavor="igraph", directed=False, n_iterations=2)
        out.obs[col] = work.obs[col].astype(str).values
    return out


def filter_confounding_genes(genes: list[str]) -> tuple[list[int], list[str]]:
    """Filters out confounding gene families and deduplicates gene symbols:
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
    seen: set[str] = set()

    for idx, g in enumerate(genes):
        clean_g = str(g).strip().upper()
        if (
            len(clean_g) > 2
            and clean_g not in {"NAN", "NONE", "UNKNOWN"}
            and clean_g not in seen
            and not exclude_pattern.search(clean_g)
        ):
            seen.add(clean_g)
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
    adata_res = load_reference_adata(config)
    match adata_res:
        case Failure(err):
            return Failure(err)
        case Success(adata):
            pass

    adata = _ensure_resolutions_clustered(adata, config.resolutions, config.cluster_prefix)
    cells: list[str] = adata.obs_names.tolist()
    print(f"Loaded AnnData with {adata.n_obs} cells.")

    # Determine linear TPM matrix
    if "tpm" in adata.layers:
        print("Extracting linear TPM directly from AnnData.layers['tpm']...")
        raw_genes = _extract_adata_gene_symbols(adata)
        valid_gene_idx, filtered_genes = filter_confounding_genes(raw_genes)
        print(f"Retained {len(filtered_genes)} protein-coding genes after confounding gene exclusion.")
        tpm_mat = adata.layers["tpm"]
        if issparse(tpm_mat):
            lin_matrix = tpm_mat[:, valid_gene_idx].toarray().T
        else:
            lin_matrix = np.asarray(tpm_mat[:, valid_gene_idx]).T
        common_cells = cells
        cell_mask = [True] * len(cells)
    elif adata.uns.get("expression_type") == "tpm":
        print("Extracting linear TPM directly from AnnData.X (expression_type=tpm)...")
        raw_genes = _extract_adata_gene_symbols(adata)
        valid_gene_idx, filtered_genes = filter_confounding_genes(raw_genes)
        print(f"Retained {len(filtered_genes)} protein-coding genes after confounding gene exclusion.")
        if issparse(adata.X):
            lin_matrix = adata.X[:, valid_gene_idx].toarray().T
        else:
            lin_matrix = np.asarray(adata.X[:, valid_gene_idx]).T
        common_cells = cells
        cell_mask = [True] * len(cells)
    else:
        print("Extracting linear expression from adata.X via expm1...")
        raw_genes = _extract_adata_gene_symbols(adata)
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
        "--cohort",
        type=str,
        default="GSE120575",
        help="Cohort ID in tme_datasets registry (default: GSE120575)",
    )
    parser.add_argument(
        "--adata",
        type=str,
        default="",
        help="Optional explicit path to preprocessed AnnData file",
    )
    parser.add_argument(
        "--tpm",
        type=str,
        default="",
        help="Optional explicit path to legacy TPM parquet file",
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
        default="output/manuscript_pipeline/3f1_build_sf_reference",
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
        cohort=args.cohort,
        adata_path=Path(args.adata) if args.adata else None,
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
