#!/usr/bin/env python3
"""
Step 1d: Build Multi-Dataset Combined References (Criteria-Driven and Random Permutations)
Unifies single-cell datasets using Harmony batch integration across 8 Leiden clustering resolutions.
Supports criteria-driven combinations (tissue-matched, platform-matched, ICI-treated, domain-shift)
and random dataset combinations (permutation null models).
Outputs reference_phi_{combo_id}_res{res}.parquet and reference_resolution_metrics_{combo_id}.parquet.
"""

from __future__ import annotations

import argparse
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
from scipy.sparse import csr_matrix  # type: ignore

# Local imports
sys.path.append(str(Path(__file__).resolve().parent))
from importlib import import_module

build_mod = import_module("01c_build_dataset_reference")
load_jerby_arnon = build_mod.load_jerby_arnon
load_maynard = build_mod.load_maynard
load_ma_liver = build_mod.load_ma_liver
load_yost = build_mod.load_yost
select_top_markers = build_mod.select_top_markers
compute_cluster_means = build_mod.compute_cluster_means

DEFAULT_RESOLUTIONS: Final[tuple[float, ...]] = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)


class ComboConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    combo_id: str
    combo_name: str
    category: str  # "Criteria-Combined", "Random-Combined", "Cell-Subsample"
    datasets: tuple[str, ...]
    out_dir: Path
    resolutions: tuple[float, ...] = DEFAULT_RESOLUTIONS
    subsample_fraction: float = 1.0
    top_markers: int = 35
    random_seed: int = 42


CRITERIA_COMBOS: Final[dict[str, tuple[str, str, tuple[str, ...]]]] = {
    "combo_melanoma": ("Melanoma-Duo", "Criteria-Combined", ("jerby", "yost")),
    "combo_plat_ss2": ("Platform-SS2-Duo", "Criteria-Combined", ("jerby", "maynard")),
    "combo_plat_10x": ("Platform-10x-Duo", "Criteria-Combined", ("ma", "yost")),
    "combo_cross_tissue": ("Cross-Tissue-Contrast", "Criteria-Combined", ("maynard", "ma")),
    "combo_tri_ici": ("ICI-Trio", "Criteria-Combined", ("jerby", "maynard", "yost")),
    "all_4_combined": ("All-Datasets-Combined", "Criteria-Combined", ("jerby", "maynard", "ma", "yost")),
}

RANDOM_COMBOS: Final[dict[str, tuple[str, str, tuple[str, ...]]]] = {
    "random_pair1": ("Random-Pair-1", "Random-Combined", ("jerby", "ma")),
    "random_pair2": ("Random-Pair-2", "Random-Combined", ("maynard", "yost")),
    "random_pair3": ("Random-Pair-3", "Random-Combined", ("jerby", "yost")),
    "random_triplet1": ("Random-Triplet-1", "Random-Combined", ("jerby", "maynard", "ma")),
    "random_triplet2": ("Random-Triplet-2", "Random-Combined", ("maynard", "ma", "yost")),
    "random_triplet3": ("Random-Triplet-3", "Random-Combined", ("jerby", "ma", "yost")),
    "random_quad1": ("Random-Quadruplet-1", "Random-Combined", ("jerby", "maynard", "ma", "yost")),
}


def load_dataset_cached(
    d_id: str,
    raw_dir: Path,
    repo_root: Path,
    cache: dict[str, tuple[ad.AnnData, npt.NDArray[np.float64], list[str]]],
) -> Result[tuple[ad.AnnData, npt.NDArray[np.float64], list[str]], str]:
    """Load dataset or return from memory cache."""
    if d_id in cache:
        return Success(cache[d_id])

    match d_id:
        case "jerby":
            res = load_jerby_arnon(raw_dir)
        case "maynard":
            res = load_maynard(repo_root)
        case "ma":
            res = load_ma_liver(raw_dir)
        case "yost":
            res = load_yost(raw_dir)
        case _:
            return Failure(f"Unknown dataset ID: {d_id}")

    match res:
        case Success(data):
            cache[d_id] = data
            return Success(data)
        case Failure(err):
            return Failure(err)


def build_combination_reference(
    config: ComboConfig,
    raw_dir: Path,
    repo_root: Path,
    dataset_cache: dict[str, tuple[ad.AnnData, npt.NDArray[np.float64], list[str]]],
) -> Result[Path, str]:
    """Integrates datasets and constructs multi-resolution deconvolution reference signatures."""
    print(f"\n==========================================================================")
    print(f"Building Combined Reference: {config.combo_name} ({config.combo_id})")
    print(f"Category: {config.category} | Datasets: {config.datasets}")
    print(f"==========================================================================")

    loaded_datasets: list[tuple[str, ad.AnnData, npt.NDArray[np.float64], list[str]]] = []
    for d_id in config.datasets:
        d_res = load_dataset_cached(d_id, raw_dir, repo_root, dataset_cache)
        match d_res:
            case Failure(err):
                return Failure(f"Failed to load component {d_id}: {err}")
            case Success((adata, lin_mat, genes)):
                loaded_datasets.append((d_id, adata, lin_mat, genes))

    # Determine common genes
    common_genes_set = set(loaded_datasets[0][3])
    for _, _, _, g_list in loaded_datasets[1:]:
        common_genes_set = common_genes_set.intersection(set(g_list))

    common_genes = sorted(list(common_genes_set))
    print(f"Found {len(common_genes)} shared genes across {len(loaded_datasets)} datasets.")
    if len(common_genes) < 500:
        return Failure(f"Too few overlapping genes ({len(common_genes)}) across datasets.")

    # Align each dataset's linear matrix to common_genes
    mats_sub: list[npt.NDArray[np.float64]] = []
    cell_dataset_labels: list[str] = []
    cell_names_list: list[str] = []

    for d_id, adata, lin_mat, genes in loaded_datasets:
        gene_to_idx = {g: idx for idx, g in enumerate(genes)}
        sub_indices = [gene_to_idx[g] for g in common_genes]
        aligned_mat = lin_mat[sub_indices, :]  # (common_genes, n_cells)
        mats_sub.append(aligned_mat)
        n_c = aligned_mat.shape[1]
        cell_dataset_labels.extend([d_id] * n_c)
        cell_names_list.extend([f"{d_id}_{c}" for c in adata.obs_names])

    merged_lin_mat = np.hstack(mats_sub)  # (common_genes, total_cells)
    total_cells = merged_lin_mat.shape[1]
    print(f"Merged expression matrix shape: {merged_lin_mat.shape} ({total_cells} total cells).")

    # Cell-level subsampling if requested
    if config.subsample_fraction < 1.0:
        n_sample = max(500, int(round(total_cells * config.subsample_fraction)))
        np.random.seed(config.random_seed)
        chosen_idx = np.sort(np.random.choice(total_cells, size=n_sample, replace=False))
        merged_lin_mat = merged_lin_mat[:, chosen_idx]
        cell_dataset_labels = [cell_dataset_labels[i] for i in chosen_idx]
        cell_names_list = [cell_names_list[i] for i in chosen_idx]
        total_cells = merged_lin_mat.shape[1]
        print(f"Subsampled {total_cells} cells ({config.subsample_fraction*100:.0f}%) for {config.combo_name}.")

    # Construct combined AnnData for integration & clustering
    obs_df = pd.DataFrame(
        {"dataset": cell_dataset_labels},
        index=cell_names_list,
    )
    adata_comb = ad.AnnData(
        X=csr_matrix(merged_lin_mat.T, dtype=np.float32),
        obs=obs_df,
        var=pd.DataFrame(index=common_genes),
    )

    # Scanpy preprocessing & PCA
    print("Preprocessing merged reference...")
    sc.pp.log1p(adata_comb)
    sc.pp.highly_variable_genes(adata_comb, n_top_genes=2000, batch_key="dataset", subset=False)
    sc.tl.pca(adata_comb, svd_solver="arpack", random_state=config.random_seed)

    # Harmony batch integration
    print("Running Harmony batch correction across datasets...")
    try:
        import harmonypy
        harmony_out = harmonypy.run_harmony(
            adata_comb.obsm["X_pca"].astype(np.float64),
            adata_comb.obs,
            "dataset",
            max_iter_harmony=15,
            random_state=config.random_seed,
        )
        z_corr = harmony_out.Z_corr
        adata_comb.obsm["X_pca_harmony"] = z_corr if z_corr.shape[0] == adata_comb.n_obs else z_corr.T
        use_rep = "X_pca_harmony"
    except Exception as harm_err:
        print(f"Harmony integration fallback to PCA ({harm_err})")
        use_rep = "X_pca"

    sc.pp.neighbors(adata_comb, use_rep=use_rep, n_neighbors=15, n_pcs=30, random_state=config.random_seed)

    # Multi-resolution Leiden clustering & signature extraction
    resolution_records: list[dict[str, object]] = []
    primary_out = config.out_dir / f"reference_phi_{config.combo_id}_res0.5.parquet"

    for res in config.resolutions:
        sc.tl.leiden(
            adata_comb,
            resolution=res,
            key_added=f"leiden_{res}",
            flavor="igraph",
            directed=False,
            n_iterations=2,
        )
        cluster_series = adata_comb.obs[f"leiden_{res}"].astype(str)
        unique_clusters = sorted(cluster_series.unique().tolist())
        clusters = cluster_series.tolist()
        n_cl = len(unique_clusters)

        # Compute cluster mean vectors in linear space
        mean_matrix = compute_cluster_means(merged_lin_mat, clusters, unique_clusters)
        nonzero_mask = np.sum(mean_matrix, axis=0) > 0.0
        final_genes = [common_genes[i] for i, nz in enumerate(nonzero_mask) if nz]
        final_mean_matrix = mean_matrix[:, nonzero_mask]

        # Select top markers
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
        print(f"  [{config.combo_name}] Res {res}: {n_cl} clusters, {len(signature_genes)} sig genes, kappa = {cond_num:.2f}")

        # Save wide signature matrix
        wide_dict: dict[str, object] = {"cluster": list(unique_clusters)}
        for g_idx, g in enumerate(signature_genes):
            wide_dict[g] = [float(phi_signature[c_idx, g_idx]) for c_idx in range(n_cl)]
        df_phi_wide = pl.DataFrame(wide_dict)

        out_phi = config.out_dir / f"reference_phi_{config.combo_id}_res{res}.parquet"
        df_phi_wide.write_parquet(out_phi)

        out_markers = config.out_dir / f"reference_marker_genes_{config.combo_id}_res{res}.parquet"
        df_markers.write_parquet(out_markers)

        tidy_records = []
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
        df_phi_tidy.write_parquet(config.out_dir / f"reference_phi_tidy_{config.combo_id}_res{res}.parquet")

        resolution_records.append(
            {
                "reference_type": config.combo_name,
                "resolution": float(res),
                "n_clusters": n_cl,
                "n_signature_genes": len(signature_genes),
                "condition_number": cond_num,
                "category": config.category,
                "n_datasets": len(config.datasets),
                "n_cells": total_cells,
                "subsample_fraction": float(config.subsample_fraction),
            }
        )

    # Save resolution metrics
    df_metrics = pl.DataFrame(resolution_records)
    out_metrics = config.out_dir / f"reference_resolution_metrics_{config.combo_id}.parquet"
    df_metrics.write_parquet(out_metrics)
    print(f"Saved {config.combo_name} metrics to: {out_metrics}")

    return Success(primary_out)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 1d: Build multi-dataset combination references across 8 resolutions."
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
        default="output/output/sade_feldman_deconv_validation",
        help="Directory to save output reference parquets",
    )
    parser.add_argument(
        "--combos",
        type=str,
        default="all",
        help="Comma-separated combo IDs or 'all', 'criteria', 'random'",
    )
    parser.add_argument(
        "--subsample-combos",
        type=str,
        default="",
        help="Comma-separated combo IDs to evaluate 10%%-100%% subsamplings on (e.g. all_4_combined,random_triplet1)",
    )
    parser.add_argument(
        "--subsample-fractions",
        type=str,
        default="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0",
        help="Comma-separated cell subsampling fractions",
    )
    parser.add_argument(
        "--subsample-res",
        type=str,
        default="0.5",
        help="Comma-separated resolutions for cell subsamplings (default 0.5)",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    raw_dir = Path(args.raw_dir)
    repo_root = Path(__file__).resolve().parent.parent.parent

    target_combos: dict[str, tuple[str, str, tuple[str, ...]]] = {}
    if args.combos in ("all", "criteria"):
        target_combos.update(CRITERIA_COMBOS)
    if args.combos in ("all", "random"):
        target_combos.update(RANDOM_COMBOS)
    if args.combos not in ("all", "criteria", "random"):
        selected = [c.strip() for c in args.combos.split(",") if c.strip()]
        for c in selected:
            if c in CRITERIA_COMBOS:
                target_combos[c] = CRITERIA_COMBOS[c]
            elif c in RANDOM_COMBOS:
                target_combos[c] = RANDOM_COMBOS[c]
            else:
                print(f"Warning: combo '{c}' not found.")

    print(f"Queued {len(target_combos)} combined reference models to build.")

    dataset_cache: dict[str, tuple[ad.AnnData, npt.NDArray[np.float64], list[str]]] = {}
    successful = 0
    total_queued = len(target_combos)

    for combo_id, (combo_name, cat, ds_tuple) in target_combos.items():
        config = ComboConfig(
            combo_id=combo_id,
            combo_name=combo_name,
            category=cat,
            datasets=ds_tuple,
            out_dir=out_dir,
        )
        match build_combination_reference(config, raw_dir, repo_root, dataset_cache):
            case Success(_):
                successful += 1
            case Failure(err):
                print(f"Error building {combo_name}: {err}", file=sys.stderr)

    # Subsampling ladder execution
    if args.subsample_combos:
        sub_combos = [c.strip() for c in args.subsample_combos.split(",") if c.strip()]
        fractions = [float(f.strip()) for f in args.subsample_fractions.split(",") if f.strip()]
        sub_resolutions = tuple(float(r.strip()) for r in args.subsample_res.split(",") if r.strip())
        
        for c_id in sub_combos:
            combo_info = CRITERIA_COMBOS.get(c_id) or RANDOM_COMBOS.get(c_id)
            if not combo_info:
                print(f"Warning: subsample target combo '{c_id}' not found.")
                continue
            base_name, base_cat, ds_tuple = combo_info
            for frac in fractions:
                sub_id = f"{c_id}_sub{frac:.2f}"
                pct_str = f"{int(round(frac * 100))}%"
                sub_name = f"{base_name} ({pct_str} Cells)"
                total_queued += 1
                
                sub_config = ComboConfig(
                    combo_id=sub_id,
                    combo_name=sub_name,
                    category="Cell-Subsample",
                    datasets=ds_tuple,
                    out_dir=out_dir,
                    resolutions=sub_resolutions,
                    subsample_fraction=frac,
                )
                match build_combination_reference(sub_config, raw_dir, repo_root, dataset_cache):
                    case Success(_):
                        successful += 1
                    case Failure(err):
                        print(f"Error building {sub_name}: {err}", file=sys.stderr)

    print(f"\nCompleted {successful}/{total_queued} reference tasks successfully.")
    if successful == total_queued:
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
