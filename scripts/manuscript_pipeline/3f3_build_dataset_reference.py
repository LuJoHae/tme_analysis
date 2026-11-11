#!/usr/bin/env python3
"""
Step 1c: Build Multi-Resolution Deconvolution References from Additional Single-Cell Datasets.
Supports:
- Jerby-Arnon (GSE115978, Melanoma Smart-seq2)
- Yost et al. (GSE123813, BCC/SCC 10x Chromium)
- PanCancer-Atlas (10x Chromium)
Follows functional programming principles: pure functions, immutable data structures,
and explicit error handling with Result monads.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path
from typing import Final

import anndata as ad  # type: ignore
import numpy as np
import numpy.typing as npt
import pandas as pd
import polars as pl
import scanpy as sc  # type: ignore
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy.sparse import csr_matrix, issparse  # type: ignore


DEFAULT_RESOLUTIONS: Final[tuple[float, ...]] = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)


class DatasetReferenceConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    dataset_id: str
    dataset_name: str
    raw_dir: Path
    out_dir: Path
    resolutions: tuple[float, ...] = DEFAULT_RESOLUTIONS
    top_markers: int = 35
    n_subsample: int = 20000
    random_seed: int = 42


def filter_confounding_genes(genes: list[str]) -> tuple[list[int], list[str]]:
    """Filters out confounding gene families:
    - Ribosomal proteins (RPS*, RPL*)
    - Mitochondrial genes (MT-*)
    - Immunoglobulins (IGH*, IGK*, IGL*)
    - Pseudogenes & non-coding RNAs (RP11*, AC00*, RNU*, RNA5S*, LINC*, CTD-*)
    - Tumor antigens (MLANA, PMEL, TYR, DCT, MITF, S100B, MAGEA*)
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
            and clean_g not in {"NAN", "NONE", "UNKNOWN", "<NA>"}
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
    """Compute cluster mean expression vectors in linear space.
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
    """Select top marker genes per cluster based on one-vs-rest fold change in linear space."""
    n_clusters, _ = mean_matrix.shape
    records: list[dict[str, object]] = []

    for c_idx, cl in enumerate(unique_clusters):
        cl_mean = mean_matrix[c_idx, :]
        rest_indices = [i for i in range(n_clusters) if i != c_idx]
        rest_mean = (
            np.mean(mean_matrix[rest_indices, :], axis=0)
            if rest_indices
            else np.ones_like(cl_mean) * 1e-6
        )

        fold_change = (cl_mean + 1.0) / (rest_mean + 1.0)
        log2_fc = np.log2((cl_mean + 1e-3) / (rest_mean + 1e-3))

        valid_candidate_mask = cl_mean >= 1.0
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


def _load_preprocessed_sc_dataset(
    dataset_id: str,
    label: str,
    max_cells: int | None = None,
) -> Result[tuple[ad.AnnData, npt.NDArray[np.float64], list[str]], str]:
    """Load preprocessed single-cell AnnData via tme_datasets, extract linear TPM/CPM, and cluster."""
    try:
        from tme_datasets import load_dataset

        match load_dataset(dataset_id, auto_download=False, subsample_n=max_cells):
            case Failure(err):
                return Failure(err)
            case Success(raw_adata):
                adata_src = raw_adata

        raw_genes = (
            [str(g).strip().upper() for g in adata_src.var["gene_name"].tolist()]
            if "gene_name" in adata_src.var.columns
            else [str(g).strip().upper() for g in adata_src.var_names.tolist()]
        )
        valid_gene_idx, filtered_genes = filter_confounding_genes(raw_genes)
        if len(filtered_genes) < 200:
            return Failure(f"Insufficient filtered genes ({len(filtered_genes)}) in {dataset_id}")

        if "tpm" in adata_src.layers:
            sub = adata_src.layers["tpm"][:, valid_gene_idx]
            cell_by_gene = sub.toarray().astype(np.float64) if issparse(sub) else np.asarray(sub, dtype=np.float64)
        elif "counts" in adata_src.layers:
            sub = adata_src.layers["counts"][:, valid_gene_idx]
            arr = sub.toarray().astype(np.float64) if issparse(sub) else np.asarray(sub, dtype=np.float64)
            rsums = np.maximum(arr.sum(axis=1, keepdims=True), 1.0)
            cell_by_gene = (arr / rsums) * 1e6
        else:
            sub = adata_src.X[:, valid_gene_idx]
            arr = sub.toarray().astype(np.float64) if issparse(sub) else np.asarray(sub, dtype=np.float64)
            if float(np.max(arr)) < 35.0:
                arr = np.expm1(np.clip(arr, 0.0, 14.0))
            rsums = np.maximum(arr.sum(axis=1, keepdims=True), 1.0)
            cell_by_gene = (arr / rsums) * 1e6

        lin_mat: npt.NDArray[np.float64] = np.asarray(cell_by_gene.T, dtype=np.float64)
        adata_clean = ad.AnnData(
            X=csr_matrix(cell_by_gene, dtype=np.float32),
            obs=adata_src.obs.copy(),
            var=pd.DataFrame(index=filtered_genes),
        )
        print(f"Loaded {dataset_id} ({label}) via tme_datasets ({adata_clean.n_obs} cells, {len(filtered_genes)} genes). Clustering...")
        sc.pp.log1p(adata_clean)
        sc.pp.highly_variable_genes(adata_clean, n_top_genes=min(2000, len(filtered_genes)), subset=False)
        sc.tl.pca(adata_clean, svd_solver="arpack", random_state=42)
        sc.pp.neighbors(adata_clean, n_neighbors=15, n_pcs=min(30, adata_clean.obsm["X_pca"].shape[1]), random_state=42)
        for res in DEFAULT_RESOLUTIONS:
            sc.tl.leiden(adata_clean, resolution=res, key_added=f"leiden_{res}", flavor="igraph", directed=False, n_iterations=2)
            n_cl = adata_clean.obs[f"leiden_{res}"].nunique()
            print(f"  {label} res={res}: {n_cl} clusters")
        return Success((adata_clean, lin_mat, filtered_genes))
    except Exception as exc:
        return Failure(f"tme_datasets load failed for {dataset_id}: {exc}")


def load_jerby_arnon(raw_dir: Path) -> Result[tuple[ad.AnnData, npt.NDArray[np.float64], list[str]], str]:
    """Load Jerby-Arnon (GSE115978) single-cell data strictly via tme_datasets.load_dataset()."""
    _ = raw_dir
    return _load_preprocessed_sc_dataset("GSE115978", "Jerby-Arnon")


def load_maynard(repo_root: Path) -> Result[tuple[ad.AnnData, npt.NDArray[np.float64], list[str]], str]:
    """Load Maynard et al. 2020 NSCLC single-cell AnnData strictly via tme_datasets.load_dataset()."""
    _ = repo_root
    return _load_preprocessed_sc_dataset("Maynard_NSCLC", "Maynard NSCLC", max_cells=5000)


def load_ma_liver(raw_dir: Path) -> Result[tuple[ad.AnnData, npt.NDArray[np.float64], list[str]], str]:
    """Load Ma et al. (GSE125449) liver cancer dataset strictly via tme_datasets.load_dataset()."""
    _ = raw_dir
    return _load_preprocessed_sc_dataset("GSE125449", "Ma Liver")


def load_yost(raw_dir: Path, n_cells: int = 3500) -> Result[tuple[ad.AnnData, npt.NDArray[np.float64], list[str]], str]:
    """Load Yost et al. (GSE123813) BCC/SCC skin carcinoma dataset strictly via tme_datasets.load_dataset()."""
    _ = raw_dir
    return _load_preprocessed_sc_dataset("GSE123813", "Yost BCC", max_cells=n_cells)


def build_reference_signatures(
    adata: ad.AnnData,
    lin_matrix: npt.NDArray[np.float64],
    filtered_genes: list[str],
    config: DatasetReferenceConfig,
) -> Result[Path, str]:
    """Computes cluster mean signatures and exports multi-resolution parquet reference matrices."""
    config.out_dir.mkdir(parents=True, exist_ok=True)
    resolution_records: list[dict[str, object]] = []

    primary_phi_path = config.out_dir / f"reference_phi_{config.dataset_id}_res0.5.parquet"

    for res in config.resolutions:
        cluster_col = f"leiden_{res}"
        if cluster_col not in adata.obs.columns:
            continue

        cluster_series = adata.obs[cluster_col].astype(str)
        unique_clusters: list[str] = sorted(cluster_series.unique().tolist())
        clusters: list[str] = cluster_series.tolist()

        print(f"\n--- Processing {config.dataset_name} Resolution res={res} ({len(unique_clusters)} clusters) ---")
        mean_matrix = compute_cluster_means(lin_matrix, clusters, unique_clusters)

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

        cond_num = float(np.linalg.cond(phi_signature))
        print(f"[{config.dataset_name}] Res {res}: {len(unique_clusters)} clusters, {len(signature_genes)} signature genes, condition number = {cond_num:.2f}")

        resolution_records.append(
            {
                "reference_type": config.dataset_name,
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

        out_res_wide = config.out_dir / f"reference_phi_{config.dataset_id}_res{res}.parquet"
        df_phi_wide.write_parquet(out_res_wide)
        print(f"Saved signature Phi (res={res}) to: {out_res_wide}")

        # 2. Save markers table
        out_res_markers = config.out_dir / f"reference_marker_genes_{config.dataset_id}_res{res}.parquet"
        df_markers.write_parquet(out_res_markers)

        # 3. Save tidy Phi
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
        df_phi_tidy.write_parquet(config.out_dir / f"reference_phi_tidy_{config.dataset_id}_res{res}.parquet")

    # Save resolution benchmarking metrics table
    df_metrics = pl.DataFrame(resolution_records)
    out_metrics = config.out_dir / f"reference_resolution_metrics_{config.dataset_id}.parquet"
    df_metrics.write_parquet(out_metrics)
    print(f"\nSaved {config.dataset_name} multi-resolution metrics to: {out_metrics}")

    return Success(primary_phi_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 1c: Build multi-resolution deconvolution reference from arbitrary single-cell datasets."
    )
    parser.add_argument(
        "--dataset",
        type=str,
        default="jerby",
        choices=["jerby", "maynard", "ma", "yost"],
        help="Dataset identifier to process: 'jerby' (GSE115978), 'maynard' (NSCLC), 'ma' (GSE125449), or 'yost' (GSE123813)",
    )
    parser.add_argument(
        "--raw-dir",
        type=str,
        default="data/raw",
        help="Path to directory containing raw dataset folders",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="output/manuscript_pipeline/3f3_build_dataset_reference",
        help="Directory to save output reference parquets",
    )
    parser.add_argument(
        "--top-markers",
        type=int,
        default=35,
        help="Number of top marker genes per cluster",
    )
    parser.add_argument(
        "--n-subsample",
        type=int,
        default=20000,
        help="Subsample count for large datasets",
    )
    args = parser.parse_args()

    name_map = {
        "jerby": "Jerby-Arnon",
        "maynard": "Maynard-NSCLC",
        "ma": "Ma-Liver",
        "yost": "Yost-BCC",
    }

    config = DatasetReferenceConfig(
        dataset_id=args.dataset,
        dataset_name=name_map.get(args.dataset, args.dataset),
        raw_dir=Path(args.raw_dir),
        out_dir=Path(args.out_dir),
        top_markers=args.top_markers,
        n_subsample=args.n_subsample,
    )

    print(f"=== Building Reference for {config.dataset_name} ({config.dataset_id}) ===")
    repo_root = Path(__file__).resolve().parent.parent.parent

    if args.dataset == "jerby":
        load_res = load_jerby_arnon(config.raw_dir)
    elif args.dataset == "maynard":
        load_res = load_maynard(repo_root)
    elif args.dataset == "ma":
        load_res = load_ma_liver(config.raw_dir)
    elif args.dataset == "yost":
        load_res = load_yost(config.raw_dir)
    else:
        print(f"Unsupported dataset: {args.dataset}", file=sys.stderr)
        sys.exit(1)

    match load_res:
        case Failure(err):
            print(f"Error loading dataset: {err}", file=sys.stderr)
            sys.exit(1)
        case Success((adata, lin_mat, genes)):
            build_res = build_reference_signatures(adata, lin_mat, genes, config)
            match build_res:
                case Success(out_p):
                    print(f"Successfully generated reference signatures at: {out_p}")
                    sys.exit(0)
                case Failure(err):
                    print(f"Error building reference: {err}", file=sys.stderr)
                    sys.exit(1)


if __name__ == "__main__":
    main()
