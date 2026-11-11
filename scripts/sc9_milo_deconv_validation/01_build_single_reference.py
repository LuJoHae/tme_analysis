#!/usr/bin/env python3
"""Step 1a: Build Multi-Resolution Deconvolution References for Single-Cell Datasets.
Supports the 9 clinical single-cell response cohorts from the Milo cross-cohort study:
  - GSE120575 (Sade-Feldman et al. 2018, Melanoma)
  - CELLxGENE_7b20c613 (Gondal et al. 2025, Pan-Cancer / Melanoma)
  - CELLxGENE_05a8c945 (Marteau et al. 2026, Colorectal)
  - CELLxGENE_6f9de485 (Yan et al. 2026, Breast)
  - GSE207422 (Hu et al. 2023, NSCLC)
  - GSE243013 (Liu et al. 2025, NSCLC)
  - GSE233203 (Hong et al. 2025, NSCLC)
  - GSE200996 (Luoma et al. 2022, HNSCC)
  - GSE316195 (Raufi et al. 2026, PDAC)
along with legacy single-cell reference datasets (Jerby-Arnon, Maynard, Ma, Yost).

Follows functional programming standards: pure functions, immutable data structures,
and explicit error handling using Returns Result monads.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Final

import anndata as ad  # type: ignore[import-untyped]
import numpy as np
import numpy.typing as npt
import pandas as pd
import polars as pl
import scanpy as sc  # type: ignore[import-untyped]
from pydantic import BaseModel, ConfigDict, Field
from returns.result import Failure, Result, Success
from scipy.sparse import csr_matrix, issparse  # type: ignore[import-untyped]

# Ensure repo packages are discoverable
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(REPO_ROOT / "packages" / "tme_datasets" / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "packages" / "tme_datasets" / "src"))

DEFAULT_RESOLUTIONS: Final[tuple[float, ...]] = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)

COHORT_CANONICAL_NAMES: Final[dict[str, str]] = {
    "GSE120575": "Sade-Feldman",
    "CELLxGENE_7b20c613": "Gondal-PanCancer",
    "CELLxGENE_05a8c945": "Marteau-Colorectal",
    "CELLxGENE_6f9de485": "Yan-Breast",
    "GSE207422": "Hu-NSCLC",
    "GSE243013": "Liu-NSCLC",
    "GSE233203": "Hong-NSCLC",
    "GSE200996": "Luoma-HNSCC",
    "GSE316195": "Raufi-PDAC",
    # Legacy / Aux cohorts
    "GSE115978": "Jerby-Arnon",
    "Maynard_NSCLC": "Maynard-NSCLC",
    "GSE125449": "Ma-Liver",
    "GSE123813": "Yost-BCC",
}

COHORT_ALIASES: Final[dict[str, str]] = {
    "sf": "GSE120575",
    "sade_feldman": "GSE120575",
    "gondal": "CELLxGENE_7b20c613",
    "marteau": "CELLxGENE_05a8c945",
    "yan": "CELLxGENE_6f9de485",
    "hu": "GSE207422",
    "liu": "GSE243013",
    "hong": "GSE233203",
    "luoma": "GSE200996",
    "raufi": "GSE316195",
    "jerby": "GSE115978",
    "maynard": "Maynard_NSCLC",
    "ma": "GSE125449",
    "yost": "GSE123813",
}


class SingleReferenceConfig(BaseModel):
    """Immutable configuration for single-cell reference signature extraction."""

    model_config = ConfigDict(frozen=True)

    cohort_id: str
    display_name: str
    out_dir: Path
    adata_path: Path | None = None
    resolutions: tuple[float, ...] = DEFAULT_RESOLUTIONS
    top_markers: int = Field(default=35, ge=5)
    n_subsample: int = Field(default=20000, ge=100)
    random_seed: int = Field(default=42)


def extract_gene_symbols_from_var(var_df: pd.DataFrame) -> list[str]:
    """Extract canonical gene symbols from AnnData var DataFrame."""
    for col in ["gene_symbols", "feature_name", "gene_name", "GeneSymbol", "gene_symbol", "var_names"]:
        if col in var_df.columns:
            vals = [str(g).strip().upper() for g in var_df[col].tolist()]
            if any(not v.startswith("ENSG") for v in vals[:50]):
                return vals
    return [str(g).strip().upper() for g in var_df.index.tolist()]


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


def load_single_cell_adata(
    config: SingleReferenceConfig,
) -> Result[tuple[ad.AnnData, npt.NDArray[np.float64], list[str], int], str]:
    """Pure boundary to load and cluster single-cell dataset via tme_datasets or explicit path."""
    try:
        from tme_datasets import load_dataset

        if config.adata_path and config.adata_path.exists():
            print(f"Loading AnnData from explicit path (fast h5py mode): {config.adata_path}")
            import h5py
            import anndata._io.specs as specs

            with h5py.File(config.adata_path, "r") as f:
                if "X" not in f:
                    return Failure(f"No 'X' dataset found in {config.adata_path}")

                x_elem = f["X"]
                if isinstance(x_elem, h5py.Group) and "indptr" in x_elem:
                    indptr = x_elem["indptr"][:]
                    total_cells = len(indptr) - 1
                    n_genes = int(x_elem.attrs.get("shape", [total_cells, 0])[1])

                    if config.n_subsample and config.n_subsample < total_cells:
                        np.random.seed(config.random_seed)
                        keep_idx = np.sort(np.random.choice(total_cells, size=config.n_subsample, replace=False))
                    else:
                        keep_idx = np.arange(total_cells)

                    dset_data = x_elem["data"]
                    dset_indices = x_elem["indices"]
                    sub_data, sub_indices, sub_indptr = [], [], [0]
                    for row in keep_idx:
                        start, end = indptr[row], indptr[row + 1]
                        sub_data.append(dset_data[start:end])
                        sub_indices.append(dset_indices[start:end])
                        sub_indptr.append(sub_indptr[-1] + (end - start))

                    csr = csr_matrix(
                        (np.concatenate(sub_data), np.concatenate(sub_indices), np.array(sub_indptr)),
                        shape=(len(keep_idx), n_genes),
                        dtype=np.float32,
                    )
                else:
                    total_cells = x_elem.shape[0]
                    if config.n_subsample and config.n_subsample < total_cells:
                        np.random.seed(config.random_seed)
                        keep_idx = np.sort(np.random.choice(total_cells, size=config.n_subsample, replace=False))
                    else:
                        keep_idx = np.arange(total_cells)
                    csr = csr_matrix(x_elem[keep_idx, :], dtype=np.float32)

                # Extract gene annotations
                try:
                    var_df = specs.read_elem(f["var"])
                    raw_genes = extract_gene_symbols_from_var(var_df)
                except Exception:
                    raw_genes = [f"GENE_{i}" for i in range(csr.shape[1])]

                obs_df = pd.DataFrame(index=[f"cell_{i}" for i in range(csr.shape[0])])
                adata_src = ad.AnnData(X=csr, obs=obs_df, var=pd.DataFrame(index=raw_genes))
        else:
            print(f"Loading '{config.cohort_id}' via tme_datasets (subsample_n={config.n_subsample})...")
            match load_dataset(config.cohort_id, auto_download=False, subsample_n=config.n_subsample):
                case Failure(err):
                    return Failure(f"Failed to load dataset '{config.cohort_id}': {err}")
                case Success(loaded_adata):
                    adata_src = loaded_adata
                    total_cells = getattr(adata_src, "n_obs", len(adata_src))

        raw_genes = extract_gene_symbols_from_var(adata_src.var)
        valid_gene_idx, filtered_genes = filter_confounding_genes(raw_genes)
        if len(filtered_genes) < 200:
            return Failure(f"Insufficient filtered genes ({len(filtered_genes)}) in {config.cohort_id}")

        # Extract linear counts / normalized expression
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

        print(
            f"Preprocessed '{config.cohort_id}': {adata_clean.n_obs} cells x {len(filtered_genes)} genes. "
            "Running multi-resolution Leiden clustering..."
        )
        sc.pp.log1p(adata_clean)
        sc.pp.highly_variable_genes(adata_clean, n_top_genes=min(2000, len(filtered_genes)), subset=False)
        sc.tl.pca(adata_clean, svd_solver="arpack", random_state=config.random_seed)
        n_pcs = min(30, adata_clean.obsm["X_pca"].shape[1])
        sc.pp.neighbors(adata_clean, n_neighbors=15, n_pcs=n_pcs, random_state=config.random_seed)

        for res in config.resolutions:
            col_name = f"leiden_{res}"
            sc.tl.leiden(
                adata_clean,
                resolution=res,
                key_added=col_name,
                flavor="igraph",
                directed=False,
                n_iterations=2,
            )
            n_cl = adata_clean.obs[col_name].nunique()
            print(f"  [{config.cohort_id}] res={res}: {n_cl} clusters")

        return Success((adata_clean, lin_mat, filtered_genes, total_cells))

    except Exception as exc:
        return Failure(f"Exception during single cell load/clustering for {config.cohort_id}: {exc}")


def build_single_reference_signatures(
    adata: ad.AnnData,
    lin_matrix: npt.NDArray[np.float64],
    filtered_genes: list[str],
    total_cells: int,
    config: SingleReferenceConfig,
) -> Result[Path, str]:
    """Generates multi-resolution simplex Phi matrices and benchmark metrics."""
    config.out_dir.mkdir(parents=True, exist_ok=True)
    resolution_records: list[dict[str, object]] = []

    primary_phi_path = config.out_dir / f"reference_phi_{config.cohort_id}_res0.5.parquet"
    sub_frac = float(adata.n_obs / max(1, total_cells))

    for res in config.resolutions:
        cluster_col = f"leiden_{res}"
        if cluster_col not in adata.obs.columns:
            continue

        cluster_series = adata.obs[cluster_col].astype(str)
        unique_clusters: list[str] = sorted(cluster_series.unique().tolist())
        clusters: list[str] = cluster_series.tolist()

        mean_matrix = compute_cluster_means(lin_matrix, clusters, unique_clusters)

        # Remove genes with 0 mean across all clusters
        total_mean = mean_matrix.sum(axis=0)
        nonzero_mask = total_mean > 0
        final_genes = [g for g, v in zip(filtered_genes, nonzero_mask) if v]
        final_mean_matrix = mean_matrix[:, nonzero_mask]

        df_markers = select_top_markers(
            final_mean_matrix,
            final_genes,
            unique_clusters,
            config.top_markers,
        )

        signature_genes = sorted(list(set(df_markers["gene"].to_list())))
        sig_gene_indices = [final_genes.index(g) for g in signature_genes]
        sig_mean_matrix = final_mean_matrix[:, sig_gene_indices]

        # Simplex normalization (rows sum to 1)
        row_sums = sig_mean_matrix.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        phi_signature = sig_mean_matrix / row_sums

        cond_num = float(np.linalg.cond(phi_signature))
        print(
            f"  [{config.display_name}] res={res}: {len(unique_clusters)} clusters, "
            f"{len(signature_genes)} signature genes, cond={cond_num:.2f}"
        )

        resolution_records.append(
            {
                "reference_type": config.display_name,
                "reference_category": "Single Dataset",
                "resolution": float(res),
                "cluster_col": cluster_col,
                "n_clusters": len(unique_clusters),
                "n_signature_genes": len(signature_genes),
                "condition_number": cond_num,
                "n_datasets": 1,
                "n_cells": adata.n_obs,
                "subsample_fraction": sub_frac,
            }
        )

        # 1. Wide reference Phi
        wide_dict: dict[str, object] = {"cluster": list(unique_clusters)}
        for g_idx, g in enumerate(signature_genes):
            wide_dict[g] = [float(phi_signature[c_idx, g_idx]) for c_idx in range(len(unique_clusters))]
        df_phi_wide = pl.DataFrame(wide_dict)
        out_res_wide = config.out_dir / f"reference_phi_{config.cohort_id}_res{res}.parquet"
        df_phi_wide.write_parquet(out_res_wide)

        # 2. Markers table
        out_res_markers = config.out_dir / f"reference_marker_genes_{config.cohort_id}_res{res}.parquet"
        df_markers.write_parquet(out_res_markers)

        # 3. Tidy Phi
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
        df_phi_tidy.write_parquet(config.out_dir / f"reference_phi_tidy_{config.cohort_id}_res{res}.parquet")

    # Metrics summary table
    df_metrics = pl.DataFrame(resolution_records)
    out_metrics = config.out_dir / f"reference_resolution_metrics_{config.cohort_id}.parquet"
    df_metrics.write_parquet(out_metrics)
    print(f"Saved {config.display_name} multi-resolution metrics to: {out_metrics}")

    return Success(primary_phi_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 1a: Build multi-resolution deconvolution reference from a single-cell cohort."
    )
    parser.add_argument(
        "--cohort",
        type=str,
        required=True,
        help="Cohort accession or alias (e.g. GSE120575, CELLxGENE_7b20c613, sf, gondal, hu, etc.)",
    )
    parser.add_argument(
        "--adata",
        type=str,
        default="",
        help="Optional explicit path to preprocessed .h5ad file on disk",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="output/sc9_milo_deconv_validation/01_references",
        help="Directory to save generated reference parquets",
    )
    parser.add_argument(
        "--top-markers",
        type=int,
        default=35,
        help="Number of marker genes selected per cluster",
    )
    parser.add_argument(
        "--n-subsample",
        type=int,
        default=20000,
        help="Maximum cells to subsample for clustering and signature derivation",
    )
    parser.add_argument(
        "--resolutions",
        type=str,
        default="0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0",
        help="Comma-separated clustering resolutions",
    )
    args = parser.parse_args()

    # Resolve cohort ID and canonical display name
    raw_cohort = args.cohort.strip()
    cohort_id = COHORT_ALIASES.get(raw_cohort.lower(), raw_cohort)
    display_name = COHORT_CANONICAL_NAMES.get(cohort_id, cohort_id)

    res_tuple = tuple(float(r.strip()) for r in args.resolutions.split(",") if r.strip())

    config = SingleReferenceConfig(
        cohort_id=cohort_id,
        display_name=display_name,
        out_dir=Path(args.out_dir),
        adata_path=Path(args.adata) if args.adata else None,
        resolutions=res_tuple,
        top_markers=args.top_markers,
        n_subsample=args.n_subsample,
    )

    print(f"=== Building Reference Signatures for '{config.display_name}' ({config.cohort_id}) ===")

    match load_single_cell_adata(config):
        case Failure(err):
            print(f"Error loading single-cell dataset: {err}", file=sys.stderr)
            sys.exit(1)
        case Success((adata, lin_mat, genes, total_cells)):
            match build_single_reference_signatures(adata, lin_mat, genes, total_cells, config):
                case Success(phi_path):
                    print(f"Reference generation succeeded: {phi_path}")
                    sys.exit(0)
                case Failure(err):
                    print(f"Error creating reference signatures: {err}", file=sys.stderr)
                    sys.exit(1)


if __name__ == "__main__":
    main()
