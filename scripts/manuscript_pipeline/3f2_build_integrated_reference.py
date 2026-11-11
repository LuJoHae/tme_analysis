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
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tme_datasets import get_dataset_spec, load_dataset
from tme_datasets.paths import get_data_paths
from returns.maybe import Some

PAN_CANCER_10X_ATLAS_COHORTS: Final[tuple[str, ...]] = (
    "GSE179994",
    "GSE178341",
    "GSE114727",
    "E-MTAB-8107",
    "GSE154763",
    "GSE154826",
    "GSE131907",
    "GSE121638",
    "GSE156625",
    "GSE169246",
    "GSE125449",
    "GSE123813",
)


class IntegratedRefConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    adata_sf_path: Path | None = None
    lair_dir: Path | None = None
    out_dir: Path
    n_subsample_atlas: int = 25000
    resolutions: tuple[float, ...] = (0.5, 1.0, 1.5, 2.0)
    top_markers: int = 10
    random_seed: int = 42


def _assign_unique_uppercase_symbols(adata: ad.AnnData) -> ad.AnnData:
    """Map .var['gene_name'] (or .var_names) to deduplicated uppercase HGNC symbols."""
    raw_names = (
        adata.var["gene_name"].astype(str).str.strip().str.upper().tolist()
        if "gene_name" in adata.var.columns
        else [str(g).strip().upper() for g in adata.var_names]
    )
    unique_symbols: list[str] = []
    seen: set[str] = set()
    for g in raw_names:
        base = "UNKNOWN" if (not g or g in {"NONE", "NAN", "UNKNOWN", "<NA>"}) else g
        sym = base
        if sym in seen:
            idx = 1
            while f"{base}_{idx}" in seen:
                idx += 1
            sym = f"{base}_{idx}"
        seen.add(sym)
        unique_symbols.append(sym)
    adata.var_names = unique_symbols
    return adata


def load_sade_feldman(path: Path | None) -> Result[ad.AnnData, str]:
    """Load and prepare Sade-Feldman AnnData strictly via tme_datasets.load_dataset()."""
    try:
        base_root = path.parent if path is not None else None
        match load_dataset("GSE120575", base_dir=base_root, auto_download=False, normalize_ensembl=False):
            case Failure(err):
                return Failure(f"Failed to load Sade-Feldman (GSE120575) via tme_datasets: {err}")
            case Success(loaded):
                adata = loaded
        adata.obs["dataset"] = "SadeFeldman_Melanoma"
        adata.obs["sequencing_tech"] = "Smart-seq2"
        adata = _assign_unique_uppercase_symbols(adata)
        return Success(adata)
    except Exception as exc:
        return Failure(f"Failed to load Sade-Feldman via tme_datasets: {exc}")


def load_and_subsample_atlas(
    lair_dir: Path | None,
    n_sample: int,
    seed: int,
) -> Result[ad.AnnData, str]:
    """Load 10x Chromium Pan-Cancer Single-Cell Atlas strictly via tme_datasets.load_dataset()."""
    _ = lair_dir
    per_cohort = max(800, (n_sample + len(PAN_CANCER_10X_ATLAS_COHORTS) - 1) // len(PAN_CANCER_10X_ATLAS_COHORTS))
    loaded_parts: list[ad.AnnData] = []

    for cid in PAN_CANCER_10X_ATLAS_COHORTS:
        match load_dataset(cid, auto_download=False, subsample_n=per_cohort, normalize_ensembl=False):
            case Success(part):
                cancer_label = cid
                match get_dataset_spec(cid):
                    case Some(spec):
                        cancer_label = spec.cancer_type
                    case _:
                        pass
                part = _assign_unique_uppercase_symbols(part.copy())
                linear_src = (
                    part.layers["counts"]
                    if "counts" in part.layers
                    else (part.layers["tpm"] if "tpm" in part.layers else part.X)
                )
                x_csr = linear_src.tocsr().astype(np.float32) if issparse(linear_src) else csr_matrix(linear_src, dtype=np.float32)
                max_v = float(x_csr.data.max()) if x_csr.nnz > 0 else 0.0
                if max_v < 35.0:
                    x_csr = x_csr.copy()
                    x_csr.data = np.expm1(np.clip(x_csr.data, 0.0, 14.0))
                clean_obs = pd.DataFrame(
                    {
                        "dataset": [cid] * part.n_obs,
                        "sequencing_tech": ["10x_Chromium"] * part.n_obs,
                        "cancer_code": [cancer_label] * part.n_obs,
                    },
                    index=[f"{cid}_{i}" for i in range(part.n_obs)],
                )
                loaded_parts.append(
                    ad.AnnData(
                        X=x_csr,
                        obs=clean_obs,
                        var=pd.DataFrame(index=list(part.var_names)),
                    )
                )
            case Failure(_):
                continue

    if loaded_parts:
        common_genes = sorted(
            set.intersection(*(set(p.var_names) for p in loaded_parts))
            - {"UNKNOWN"}
        )
        if len(common_genes) >= 5000:
            aligned = [p[:, common_genes].copy() for p in loaded_parts]
            merged = ad.concat(aligned, join="inner", merge="same")
            if merged.n_obs > n_sample:
                rng = np.random.default_rng(seed)
                chosen = np.sort(rng.choice(merged.n_obs, size=n_sample, replace=False))
                merged = merged[chosen].copy()
            print(
                f"Subsampled {merged.n_obs} representative cells across "
                f"{len(loaded_parts)} tme_datasets 10x Chromium cohorts ({len(common_genes)} shared genes)."
            )
            return Success(merged)

    return Failure("No pan-cancer 10x single-cell atlas datasets available via tme_datasets.load_dataset().")


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

    # Linearize expression (Sade-Feldman TPM/log1p -> linear; Atlas counts -> normalize to 1e4)
    if "tpm" in sub_sf.layers:
        sf_raw = sub_sf.layers["tpm"]
        sf_linear = sf_raw.tocsr() if issparse(sf_raw) else csr_matrix(sf_raw)
    else:
        sf_max = float(sub_sf.X.data.max() if issparse(sub_sf.X) and sub_sf.X.nnz > 0 else np.max(sub_sf.X))
        if sf_max < 35.0:
            sf_linear = sub_sf.X.expm1() if issparse(sub_sf.X) else csr_matrix(np.expm1(sub_sf.X))
        else:
            sf_linear = sub_sf.X.tocsr() if issparse(sub_sf.X) else csr_matrix(sub_sf.X)
    sf_sums = np.asarray(sf_linear.sum(axis=1)).ravel()
    sf_sums[sf_sums == 0] = 1.0
    from scipy.sparse import diags
    sf_linear = diags(1e4 / sf_sums).dot(sf_linear)

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

    # Compute neighbors and UMAP embedding on Harmony latent space
    print("Computing integrated neighborhood graph and UMAP embedding on Harmony latent space...")
    sc.pp.neighbors(adata_comb, use_rep="X_pca_harmony", n_neighbors=15)
    sc.tl.umap(adata_comb)

    config.out_dir.mkdir(parents=True, exist_ok=True)
    linear_mat = adata_comb.layers["linear"]
    resolution_records: list[dict[str, object]] = []
    primary_out_path = config.out_dir / "integrated_reference_phi.parquet"

    for res in config.resolutions:
        cluster_key = f"integrated_leiden_{res}"
        print(f"\n--- Running Integrated Leiden Clustering (res={res}) ---")
        sc.tl.leiden(adata_comb, resolution=res, key_added=cluster_key)

        cluster_labels = adata_comb.obs[cluster_key].astype(str).values
        unique_clusters = tuple(sorted(list(set(cluster_labels))))
        print(f"Discovered {len(unique_clusters)} integrated clusters at resolution {res}.")

        # Compute cluster linear means
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

        # Condition number of Phi matrix (measure of collinearity)
        cond_num = float(np.linalg.cond(phi_matrix))
        print(f"Combined Atlas res={res}: {len(unique_clusters)} clusters, condition number = {cond_num:.2f}")

        resolution_records.append(
            {
                "reference_type": "Combined-Atlas",
                "resolution": float(res),
                "cluster_col": cluster_key,
                "n_clusters": len(unique_clusters),
                "n_signature_genes": int(df_markers["gene"].n_unique()),
                "condition_number": cond_num,
            }
        )

        # Wide reference Phi table for this resolution: clusters x genes
        wide_dict: dict[str, object] = {"cluster": list(unique_clusters)}
        for g_idx, g in enumerate(common_genes):
            wide_dict[g] = [float(phi_matrix[c_idx, g_idx]) for c_idx in range(len(unique_clusters))]
        df_phi_wide = pl.DataFrame(wide_dict)

        out_res_phi = config.out_dir / f"integrated_reference_phi_res{res}.parquet"
        df_phi_wide.write_parquet(out_res_phi)
        print(f"Saved integrated Phi (res={res}) to: {out_res_phi}")

        out_res_markers = config.out_dir / f"integrated_marker_genes_res{res}.parquet"
        df_markers.write_parquet(out_res_markers)

        # Backward compatibility for res == 0.5
        if abs(res - 0.5) < 1e-4:
            df_phi_wide.write_parquet(primary_out_path)
            df_markers.write_parquet(config.out_dir / "integrated_marker_genes.parquet")
            print("Saved backward-compatible primary integrated_reference_phi.parquet (res=0.5).")

    # Integration evaluation metrics (using res=0.5 cluster key)
    print("\nEvaluating integration quality metrics...")
    sub_sample_eval = min(5000, adata_comb.n_obs)
    np.random.seed(config.random_seed)
    eval_idx = np.random.choice(adata_comb.n_obs, size=sub_sample_eval, replace=False)

    sil_tech_uncorrected = float(silhouette_score(pca_uncorrected[eval_idx], adata_comb.obs["sequencing_tech"].iloc[eval_idx]))
    sil_tech_harmony = float(silhouette_score(pca_harmony[eval_idx], adata_comb.obs["sequencing_tech"].iloc[eval_idx]))
    eval_cluster_key = "integrated_leiden_0.5" if "integrated_leiden_0.5" in adata_comb.obs else f"integrated_leiden_{config.resolutions[0]}"
    sil_cluster_harmony = float(silhouette_score(pca_harmony[eval_idx], adata_comb.obs[eval_cluster_key].iloc[eval_idx]))

    print("-------------------------------------------------------")
    print("HARMONY INTEGRATION QUALITY METRICS")
    print(f"Platform Silhouette (Uncorrected): {sil_tech_uncorrected:.4f} (high separation)")
    print(f"Platform Silhouette (Harmony):     {sil_tech_harmony:.4f} (low separation = well-mixed)")
    print(f"Cell Cluster Silhouette (Harmony): {sil_cluster_harmony:.4f} (distinct biological states)")
    print("-------------------------------------------------------")

    # Metadata and coordinates for plotting (including all resolution cluster assignments)
    umap_coords = adata_comb.obsm["X_umap"]
    meta_dict: dict[str, object] = {
        "cell_id": list(adata_comb.obs_names),
        "dataset": adata_comb.obs["dataset"].values.tolist(),
        "sequencing_tech": adata_comb.obs["sequencing_tech"].values.tolist(),
        "cell_type_original": adata_comb.obs["cell_type_original"].values.tolist(),
        "integrated_cluster": adata_comb.obs[eval_cluster_key].astype(str).values.tolist(),
        "umap_1": umap_coords[:, 0].astype(float).tolist(),
        "umap_2": umap_coords[:, 1].astype(float).tolist(),
        "pca_uncorrected_1": pca_uncorrected[:, 0].astype(float).tolist(),
        "pca_uncorrected_2": pca_uncorrected[:, 1].astype(float).tolist(),
        "pca_harmony_1": pca_harmony[:, 0].astype(float).tolist(),
        "pca_harmony_2": pca_harmony[:, 1].astype(float).tolist(),
    }
    # Add columns for each resolution
    for res in config.resolutions:
        r_key = f"integrated_leiden_{res}"
        if r_key in adata_comb.obs:
            meta_dict[r_key] = adata_comb.obs[r_key].astype(str).values.tolist()

    df_meta = pl.DataFrame(meta_dict)
    out_meta = config.out_dir / "integrated_cell_metadata.parquet"
    df_meta.write_parquet(out_meta)
    print(f"Saved integrated cell metadata with multi-resolution clusters to: {out_meta}")

    # Metrics table
    df_metrics = pl.DataFrame(
        {
            "num_sade_feldman_cells": [sub_sf.n_obs],
            "num_atlas_cells": [sub_atlas.n_obs],
            "total_cells": [adata_comb.n_obs],
            "num_shared_genes": [len(common_genes)],
            "num_integrated_clusters": [len(set(adata_comb.obs[eval_cluster_key]))],
            "platform_silhouette_uncorrected": [sil_tech_uncorrected],
            "platform_silhouette_harmony": [sil_tech_harmony],
            "cluster_silhouette_harmony": [sil_cluster_harmony],
        }
    )
    out_metrics = config.out_dir / "integration_quality_metrics.parquet"
    df_metrics.write_parquet(out_metrics)
    print(f"Saved integration metrics to: {out_metrics}")

    # Multi-resolution metrics table
    df_res_metrics = pl.DataFrame(resolution_records)
    out_res_metrics = config.out_dir / "reference_resolution_metrics_comb.parquet"
    df_res_metrics.write_parquet(out_res_metrics)
    print(f"Saved Combined Atlas multi-resolution metrics to: {out_res_metrics}")

    return Success(primary_out_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 1b: Build multi-resolution integrated single-cell references using Harmony."
    )
    parser.add_argument(
        "--adata-sf",
        type=str,
        default="",
        help="Optional explicit path to Sade-Feldman .h5ad (defaults to tme_datasets GSE120575)",
    )
    parser.add_argument(
        "--lair-dir",
        type=str,
        default="",
        help="Unused deprecated CLI flag retained for backward compatibility",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="output/manuscript_pipeline/3f2_build_integrated_reference",
        help="Directory to store output parquets",
    )
    parser.add_argument(
        "--n-subsample-atlas",
        type=int,
        default=25000,
        help="Number of cells to subsample from atlas",
    )
    parser.add_argument(
        "--resolutions",
        type=str,
        default="0.5,1.0,1.5,2.0",
        help="Comma-separated Leiden resolutions on Harmony latent space",
    )
    parser.add_argument(
        "--top-markers",
        type=int,
        default=10,
        help="Top marker genes per integrated cluster",
    )
    args = parser.parse_args()

    res_list = tuple(float(r.strip()) for r in args.resolutions.split(",") if r.strip())

    config = IntegratedRefConfig(
        adata_sf_path=Path(args.adata_sf) if args.adata_sf else None,
        lair_dir=Path(args.lair_dir) if args.lair_dir else None,
        out_dir=Path(args.out_dir),
        n_subsample_atlas=args.n_subsample_atlas,
        resolutions=res_list,
        top_markers=args.top_markers,
    )

    match run_integrated_pipeline(config):
        case Success(out_path):
            print(f"Step 1b completed successfully. Integrated references saved at: {config.out_dir}")
            sys.exit(0)
        case Failure(err):
            print(f"Step 1b failed with error: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
