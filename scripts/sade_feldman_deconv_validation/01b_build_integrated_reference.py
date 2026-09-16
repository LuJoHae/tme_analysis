#!/usr/bin/env python3
"""
Step 1b: Build Integrated Multi-Cohort Single-Cell Reference using Harmony Batch Correction.
Unifies Sade-Feldman (GSE120575, Smart-seq2) with landmark pan-cancer single-cell atlas data (10x Chromium),
performs Harmony batch integration to correct platform and cohort batch effects,
computes integrated clusters, and exports integrated_reference_phi.parquet.
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
import scanpy as sc  # type: ignore
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy.sparse import issparse, spmatrix, csr_matrix  # type: ignore
from sklearn.metrics import silhouette_score  # type: ignore

# Local repo packages
sys.path.append(str(Path(__file__).resolve().parent.parent.parent / "packages"))
import datalair  # type: ignore
from single_cell_datasets import SingleCellReference  # type: ignore


class IntegratedRefConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    adata_sf_path: Path
    lair_dir: Path
    out_dir: Path
    n_subsample_atlas: int = 25000
    resolution: float = 0.5
    top_markers: int = 10
    random_seed: int = 42


def load_sade_feldman(path: Path) -> Result[ad.AnnData, str]:
    """Load and prepare Sade-Feldman AnnData."""
    try:
        adata = ad.read_h5ad(path)
        adata.obs["dataset"] = "SadeFeldman_Melanoma"
        adata.obs["sequencing_tech"] = "Smart-seq2"
        adata.var_names = [str(g).upper() for g in adata.var_names]
        # Deduplicate gene names if any
        adata.var_names_make_unique()
        return Success(adata)
    except Exception as exc:
        return Failure(f"Failed to load Sade-Feldman from {path}: {exc}")


def load_and_subsample_atlas(lair_dir: Path, n_sample: int, seed: int) -> Result[ad.AnnData, str]:
    """Load SingleCellReference from datalair and perform stratified subsampling."""
    try:
        lair = datalair.Lair(str(lair_dir))
        ds_ref = SingleCellReference()
        paths = lair.get_dataset_filepaths(ds_ref)
        if "adata.h5ad" not in paths:
            return Failure(f"SingleCellReference adata.h5ad not found at {lair_dir}")

        h5_path = paths["adata.h5ad"]
        print(f"Reading SingleCellReference from {h5_path} (backed mode)...")
        adata = ad.read_h5ad(h5_path, backed="r")

        # Stratified sampling by cancer_code
        obs = adata.obs
        c_codes = obs["cancer_code"].unique()
        selected_idx: list[int] = []
        n_per_code = max(500, n_sample // len(c_codes))

        np.random.seed(seed)
        for code in c_codes:
            sub_indices = np.where(obs["cancer_code"] == code)[0]
            if len(sub_indices) > n_per_code:
                chosen = np.random.choice(sub_indices, size=n_per_code, replace=False)
            else:
                chosen = sub_indices
            selected_idx.extend(chosen.tolist())

        selected_idx = sorted(selected_idx[:n_sample])
        print(f"Subsampled {len(selected_idx)} representative cells from single cell atlas...")

        # Load subsample into memory
        sub_adata = adata[selected_idx].to_memory()
        sub_adata.obs["dataset"] = sub_adata.obs.get("dataset", "PanCancer_Atlas").astype(str)
        sub_adata.obs["sequencing_tech"] = "10x_Chromium"
        # Map Ensembl IDs to uppercase gene symbols using 'gene_name' column
        raw_gene_names = sub_adata.var["gene_name"].astype(str).str.upper().tolist()
        unique_symbols: list[str] = []
        seen: set[str] = set()
        for g in raw_gene_names:
            if not g or g in ("NONE", "NAN", "UNKNOWN"):
                g = "UNKNOWN"
            if g in seen:
                i = 1
                while f"{g}_{i}" in seen:
                    i += 1
                g = f"{g}_{i}"
            seen.add(g)
            unique_symbols.append(g)
        sub_adata.var_names = unique_symbols

        return Success(sub_adata)
    except Exception as exc:
        return Failure(f"Failed to load atlas data: {exc}")


def select_top_markers(
    mean_matrix: np.ndarray,
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


def run_integrated_pipeline(config: IntegratedRefConfig) -> Result[Path, str]:
    """Execute complete dataset harmonization, Harmony integration, and reference construction."""
    sf_res = load_sade_feldman(config.adata_sf_path)
    match sf_res:
        case Failure(err):
            return Failure(err)
        case Success(adata_sf):
            pass

    atlas_res = load_and_subsample_atlas(
        config.lair_dir,
        config.n_subsample_atlas,
        config.random_seed,
    )
    match atlas_res:
        case Failure(err):
            return Failure(err)
        case Success(adata_atlas):
            pass

    # Gene harmonization: intersect genes
    common_genes = sorted(list(set(adata_sf.var_names).intersection(set(adata_atlas.var_names))))
    if len(common_genes) < 5000:
        return Failure(f"Too few shared genes between Sade-Feldman and Atlas: {len(common_genes)}")

    print(f"Aligning {adata_sf.n_obs} Sade-Feldman cells and {adata_atlas.n_obs} Atlas cells on {len(common_genes)} shared genes...")

    sub_sf = adata_sf[:, common_genes].copy()
    sub_atlas = adata_atlas[:, common_genes].copy()

    # Linearize expression (Sade-Feldman log1p -> linear TPM; Atlas raw counts -> normalize to 1e4)
    if issparse(sub_sf.X):
        sf_linear = sub_sf.X.expm1()
    else:
        sf_linear = csr_matrix(np.expm1(sub_sf.X))

    if issparse(sub_atlas.X):
        atlas_linear = sub_atlas.X.copy()
    else:
        atlas_linear = csr_matrix(sub_atlas.X)

    # Normalize atlas library size to 1e4 to match Smart-seq2 scale
    atlas_sums = np.asarray(atlas_linear.sum(axis=1)).ravel()
    atlas_sums[atlas_sums == 0] = 1.0
    scale_factor = 1e4 / atlas_sums
    from scipy.sparse import diags
    atlas_linear = diags(scale_factor).dot(atlas_linear)

    # Stack linear expression matrix
    from scipy.sparse import vstack
    combined_linear = vstack([sf_linear, atlas_linear])

    # Build AnnData for integration
    obs_sf = pd.DataFrame(
        {
            "dataset": sub_sf.obs["dataset"].values,
            "sequencing_tech": sub_sf.obs["sequencing_tech"].values,
            "cell_type_original": sub_sf.obs.get("celltypist_leiden_0.5", pd.Series(["SadeFeldman"] * sub_sf.n_obs)).values,
        },
        index=[f"SF_{i}" for i in range(sub_sf.n_obs)],
    )

    obs_atlas = pd.DataFrame(
        {
            "dataset": sub_atlas.obs["dataset"].values,
            "sequencing_tech": sub_atlas.obs["sequencing_tech"].values,
            "cell_type_original": sub_atlas.obs.get("cancer_code", pd.Series(["Atlas"] * sub_atlas.n_obs)).values,
        },
        index=[f"Atlas_{i}" for i in range(sub_atlas.n_obs)],
    )

    combined_obs = pd.concat([obs_sf, obs_atlas], axis=0)
    var_df = pd.DataFrame(index=common_genes)

    adata_comb = ad.AnnData(
        X=combined_linear.copy(),
        obs=combined_obs,
        var=var_df,
    )
    # Store linear counts in layers for downstream cluster averaging
    adata_comb.layers["linear"] = combined_linear

    # Log1p transformation for PCA & Harmony
    print("Normalizing and computing PCA for integration...")
    sc.pp.normalize_total(adata_comb, target_sum=1e4)
    sc.pp.log1p(adata_comb)
    sc.pp.highly_variable_genes(adata_comb, n_top_genes=3000, inplace=True)
    sc.tl.pca(adata_comb, n_comps=40, mask_var="highly_variable")
    pca_uncorrected = adata_comb.obsm["X_pca"].copy()

    # Harmony Batch Integration
    print("Running Harmony batch correction across datasets and sequencing technologies...")
    import harmonypy  # type: ignore

    harmony_out = harmonypy.run_harmony(
        adata_comb.obsm["X_pca"].astype(np.float64),
        adata_comb.obs,
        "sequencing_tech",
        max_iter_harmony=20,
        random_state=0,
    )
    z_corr = harmony_out.Z_corr
    adata_comb.obsm["X_pca_harmony"] = z_corr if z_corr.shape[0] == adata_comb.n_obs else z_corr.T
    pca_harmony = adata_comb.obsm["X_pca_harmony"]

    # Compute neighbors and integrated Leiden clustering
    print(f"Computing integrated neighborhood graph and Leiden clustering (res={config.resolution})...")
    sc.pp.neighbors(adata_comb, use_rep="X_pca_harmony", n_neighbors=15)
    sc.tl.umap(adata_comb)
    cluster_key = f"integrated_leiden_{config.resolution}"
    sc.tl.leiden(adata_comb, resolution=config.resolution, key_added=cluster_key)

    # Integration evaluation metrics
    print("Evaluating integration quality metrics...")
    sub_sample_eval = min(5000, adata_comb.n_obs)
    np.random.seed(config.random_seed)
    eval_idx = np.random.choice(adata_comb.n_obs, size=sub_sample_eval, replace=False)

    sil_tech_uncorrected = float(silhouette_score(pca_uncorrected[eval_idx], adata_comb.obs["sequencing_tech"].iloc[eval_idx]))
    sil_tech_harmony = float(silhouette_score(pca_harmony[eval_idx], adata_comb.obs["sequencing_tech"].iloc[eval_idx]))
    sil_cluster_harmony = float(silhouette_score(pca_harmony[eval_idx], adata_comb.obs[cluster_key].iloc[eval_idx]))

    print("-------------------------------------------------------")
    print("HARMONY INTEGRATION QUALITY METRICS")
    print(f"Platform Silhouette (Uncorrected): {sil_tech_uncorrected:.4f} (high separation)")
    print(f"Platform Silhouette (Harmony):     {sil_tech_harmony:.4f} (low separation = well-mixed)")
    print(f"Cell Cluster Silhouette (Harmony): {sil_cluster_harmony:.4f} (distinct biological states)")
    print("-------------------------------------------------------")

    # Compute integrated cluster linear means
    print("Computing integrated linear cluster mean profiles...")
    cluster_labels = adata_comb.obs[cluster_key].astype(str).values
    unique_clusters = tuple(sorted(list(set(cluster_labels))))

    linear_mat = adata_comb.layers["linear"]
    cluster_means_list: list[np.ndarray] = []

    for cl in unique_clusters:
        cl_mask = cluster_labels == cl
        sub_lin = linear_mat[cl_mask]
        cl_mean = (
            np.asarray(sub_lin.mean(axis=0)).ravel()
            if issparse(sub_lin)
            else np.mean(sub_lin, axis=0)
        )
        cluster_means_list.append(np.asarray(cl_mean, dtype=np.float64))

    mean_matrix = np.vstack(cluster_means_list)

    # Normalize cluster profiles to sum to 1
    row_sums = mean_matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    phi_matrix = mean_matrix / row_sums

    # Top marker genes per integrated cluster
    df_markers = select_top_markers(
        mean_matrix,
        tuple(common_genes),
        unique_clusters,
        config.top_markers,
    )

    # Save outputs
    config.out_dir.mkdir(parents=True, exist_ok=True)

    # Wide reference Phi table: clusters x genes
    wide_dict: dict[str, object] = {"cluster": list(unique_clusters)}
    for g_idx, g in enumerate(common_genes):
        wide_dict[g] = [float(phi_matrix[c_idx, g_idx]) for c_idx in range(len(unique_clusters))]
    df_phi_wide = pl.DataFrame(wide_dict)
    out_phi = config.out_dir / "integrated_reference_phi.parquet"
    df_phi_wide.write_parquet(out_phi)
    print(f"Saved integrated reference Phi to: {out_phi}")

    # Markers table
    out_markers = config.out_dir / "integrated_marker_genes.parquet"
    df_markers.write_parquet(out_markers)
    print(f"Saved integrated marker genes to: {out_markers}")

    # Metadata and coordinates for plotting
    umap_coords = adata_comb.obsm["X_umap"]
    df_meta = pl.DataFrame(
        {
            "cell_id": list(adata_comb.obs_names),
            "dataset": adata_comb.obs["dataset"].values.tolist(),
            "sequencing_tech": adata_comb.obs["sequencing_tech"].values.tolist(),
            "cell_type_original": adata_comb.obs["cell_type_original"].values.tolist(),
            "integrated_cluster": cluster_labels.tolist(),
            "umap_1": umap_coords[:, 0].astype(float).tolist(),
            "umap_2": umap_coords[:, 1].astype(float).tolist(),
            "pca_uncorrected_1": pca_uncorrected[:, 0].astype(float).tolist(),
            "pca_uncorrected_2": pca_uncorrected[:, 1].astype(float).tolist(),
            "pca_harmony_1": pca_harmony[:, 0].astype(float).tolist(),
            "pca_harmony_2": pca_harmony[:, 1].astype(float).tolist(),
        }
    )
    out_meta = config.out_dir / "integrated_cell_metadata.parquet"
    df_meta.write_parquet(out_meta)
    print(f"Saved integrated cell metadata to: {out_meta}")

    # Metrics table
    df_metrics = pl.DataFrame(
        {
            "num_sade_feldman_cells": [sub_sf.n_obs],
            "num_atlas_cells": [sub_atlas.n_obs],
            "total_cells": [adata_comb.n_obs],
            "num_shared_genes": [len(common_genes)],
            "num_integrated_clusters": [len(unique_clusters)],
            "platform_silhouette_uncorrected": [sil_tech_uncorrected],
            "platform_silhouette_harmony": [sil_tech_harmony],
            "cluster_silhouette_harmony": [sil_cluster_harmony],
        }
    )
    out_metrics = config.out_dir / "integration_quality_metrics.parquet"
    df_metrics.write_parquet(out_metrics)
    print(f"Saved integration metrics to: {out_metrics}")

    return Success(out_phi)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 1b: Build integrated multi-cohort single-cell reference using Harmony."
    )
    parser.add_argument(
        "--adata-sf",
        type=str,
        default="/storage/halu/data/GSE120575/gse120575_processed.h5ad",
        help="Path to gse120575_processed.h5ad",
    )
    parser.add_argument(
        "--lair-dir",
        type=str,
        default="/storage/halu/lair",
        help="Path to datalair directory",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="output/sade_feldman_deconv_validation",
        help="Directory to store output parquets",
    )
    parser.add_argument(
        "--n-subsample-atlas",
        type=int,
        default=25000,
        help="Number of cells to subsample from atlas",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=0.5,
        help="Leiden resolution on Harmony latent space",
    )
    parser.add_argument(
        "--top-markers",
        type=int,
        default=10,
        help="Top marker genes per integrated cluster",
    )
    args = parser.parse_args()

    config = IntegratedRefConfig(
        adata_sf_path=Path(args.adata_sf),
        lair_dir=Path(args.lair_dir),
        out_dir=Path(args.out_dir),
        n_subsample_atlas=args.n_subsample_atlas,
        resolution=args.resolution,
        top_markers=args.top_markers,
    )

    match run_integrated_pipeline(config):
        case Success(out_path):
            print(f"Step 1b completed successfully. Integrated reference saved at: {out_path}")
            sys.exit(0)
        case Failure(err):
            print(f"Step 1b failed with error: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
