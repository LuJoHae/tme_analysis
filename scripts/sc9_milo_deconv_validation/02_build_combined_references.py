#!/usr/bin/env python3
"""Step 1b: Build Multi-Dataset Combined References via Harmony Integration.
Combines single-cell cohorts across the 9 clinical Milo cohorts:
  - Criteria combinations: All-9-Milo-Combined, NSCLC-Trio, Melanoma-Duo, Platform-10x-Quad, Solid-Tumor-Quartet
  - Random combinations / permutation baselines: Random pairs, triplets, and quadruplets
  - Cell-subsample titration ladders (10% to 100% of cells)
Aligns gene spaces, applies Harmony latent space batch correction across datasets,
performs multi-resolution Leiden clustering, and computes simplex Phi matrices.

Follows functional programming standards: pure functions, immutable data structures,
and explicit error handling using Returns Result monads.
"""

from __future__ import annotations

import argparse
import os
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

sys.path.append(str(Path(__file__).resolve().parent))
from importlib import import_module

build_mod = import_module("01_build_single_reference")
DEFAULT_RESOLUTIONS = build_mod.DEFAULT_RESOLUTIONS
compute_cluster_means = build_mod.compute_cluster_means
filter_confounding_genes = build_mod.filter_confounding_genes
select_top_markers = build_mod.select_top_markers

CRITERIA_COMBINATIONS: Final[dict[str, tuple[str, tuple[str, ...]]]] = {
    "all_9_combined": ("All-9-Milo-Combined", (
        "GSE120575", "CELLxGENE_7b20c613", "CELLxGENE_05a8c945", "CELLxGENE_6f9de485",
        "GSE207422", "GSE243013", "GSE233203", "GSE200996", "GSE316195"
    )),
    "nsclc_trio": ("NSCLC-Trio", ("GSE207422", "GSE243013", "GSE233203")),
    "melanoma_duo": ("Melanoma-Duo", ("GSE120575", "CELLxGENE_7b20c613")),
    "platform_10x_quad": ("Platform-10x-Quad", ("CELLxGENE_7b20c613", "CELLxGENE_05a8c945", "CELLxGENE_6f9de485", "GSE233203")),
    "solid_tumor_quartet": ("Solid-Tumor-Quartet", ("GSE120575", "CELLxGENE_05a8c945", "CELLxGENE_6f9de485", "GSE200996")),
}

RANDOM_COMBINATIONS: Final[dict[str, tuple[str, tuple[str, ...]]]] = {
    # Size 2 (pairs)
    "random_pair1": ("Random-Pair-1", ("GSE120575", "CELLxGENE_05a8c945")),
    "random_pair2": ("Random-Pair-2", ("GSE207422", "CELLxGENE_6f9de485")),
    "random_pair3": ("Random-Pair-3", ("GSE243013", "GSE200996")),
    # Size 3 (triplets)
    "random_triplet1": ("Random-Triplet-1", ("GSE120575", "GSE207422", "GSE316195")),
    "random_triplet2": ("Random-Triplet-2", ("CELLxGENE_7b20c613", "CELLxGENE_05a8c945", "GSE243013")),
    "random_triplet3": ("Random-Triplet-3", ("CELLxGENE_6f9de485", "GSE233203", "GSE200996")),
    # Size 4 (quads)
    "random_quad1": ("Random-Quad-1", ("GSE120575", "CELLxGENE_7b20c613", "GSE207422", "GSE316195")),
    "random_quad2": ("Random-Quad-2", ("CELLxGENE_05a8c945", "CELLxGENE_6f9de485", "GSE243013", "GSE233203")),
    "random_quad3": ("Random-Quad-3", ("CELLxGENE_7b20c613", "GSE207422", "GSE200996", "GSE316195")),
    # Size 5 (quintets)
    "random_quintet1": ("Random-Quintet-1", ("GSE120575", "CELLxGENE_7b20c613", "GSE200996", "GSE243013", "GSE316195")),
    "random_quintet2": ("Random-Quintet-2", ("CELLxGENE_05a8c945", "CELLxGENE_6f9de485", "GSE207422", "GSE233203", "GSE316195")),
    "random_quintet3": ("Random-Quintet-3", ("GSE120575", "CELLxGENE_05a8c945", "GSE207422", "GSE233203", "GSE200996")),
    # Size 6 (sextets)
    "random_sextet1": ("Random-Sextet-1", ("CELLxGENE_05a8c945", "CELLxGENE_6f9de485", "GSE120575", "GSE200996", "GSE207422", "GSE243013")),
    "random_sextet2": ("Random-Sextet-2", ("CELLxGENE_7b20c613", "CELLxGENE_05a8c945", "GSE120575", "GSE207422", "GSE233203", "GSE316195")),
    "random_sextet3": ("Random-Sextet-3", ("CELLxGENE_6f9de485", "CELLxGENE_7b20c613", "GSE207422", "GSE233203", "GSE243013", "GSE200996")),
    # Size 7 (septets)
    "random_septet1": ("Random-Septet-1", ("CELLxGENE_05a8c945", "CELLxGENE_6f9de485", "GSE120575", "GSE200996", "GSE207422", "GSE243013", "GSE316195")),
    "random_septet2": ("Random-Septet-2", ("CELLxGENE_05a8c945", "CELLxGENE_7b20c613", "GSE120575", "GSE207422", "GSE233203", "GSE243013", "GSE316195")),
    "random_septet3": ("Random-Septet-3", ("CELLxGENE_6f9de485", "CELLxGENE_7b20c613", "GSE120575", "GSE207422", "GSE233203", "GSE200996", "GSE316195")),
    # Size 8 (octets - leave-one-dataset-out controls)
    "random_octet1": ("Random-Octet-1", ("CELLxGENE_7b20c613", "CELLxGENE_05a8c945", "CELLxGENE_6f9de485", "GSE207422", "GSE243013", "GSE233203", "GSE200996", "GSE316195")),
    "random_octet2": ("Random-Octet-2", ("GSE120575", "CELLxGENE_05a8c945", "CELLxGENE_6f9de485", "GSE207422", "GSE243013", "GSE233203", "GSE200996", "GSE316195")),
    "random_octet3": ("Random-Octet-3", ("GSE120575", "CELLxGENE_7b20c613", "CELLxGENE_05a8c945", "CELLxGENE_6f9de485", "GSE243013", "GSE233203", "GSE200996", "GSE316195")),
}


class CombinedReferenceConfig(BaseModel):
    """Immutable configuration for building multi-dataset combined references."""

    model_config = ConfigDict(frozen=True)

    combo_id: str
    combo_name: str
    category: str  # "Criteria-Combined", "Random-Combined", "Cell-Subsample"
    cohort_ids: tuple[str, ...]
    out_dir: Path
    resolutions: tuple[float, ...] = DEFAULT_RESOLUTIONS
    subsample_fraction: float = Field(default=1.0, ge=0.01, le=1.0)
    cells_per_dataset: int = Field(default=10000, ge=100)
    top_markers: int = Field(default=35, ge=5)
    random_seed: int = Field(default=42)


def extract_gene_symbols_from_var(var_df: pd.DataFrame) -> list[str]:
    """Extract canonical gene symbols from AnnData var DataFrame."""
    for col in ["gene_symbols", "feature_name", "gene_name", "GeneSymbol", "gene_symbol", "var_names"]:
        if col in var_df.columns:
            vals = [str(g).strip().upper() for g in var_df[col].tolist()]
            if any(not v.startswith("ENSG") for v in vals[:50]):
                return vals
    return [str(g).strip().upper() for g in var_df.index.tolist()]


def load_cohort_for_combination(
    cohort_id: str,
    max_cells: int,
    seed: int,
) -> Result[tuple[ad.AnnData, list[str]], str]:
    """Pure boundary to load and preprocess a cohort for multi-dataset combination."""
    try:
        from tme_datasets import load_dataset
        import h5py
        import anndata._io.specs as specs

        data_dir_env = os.environ.get("TME_DATA_DIR", "/ritom/scratch/cscs/lhaeuser/data-test")
        preproc_h5ad = Path(data_dir_env) / "preprocessed" / f"{cohort_id}.h5ad"

        if preproc_h5ad.exists() and preproc_h5ad.stat().st_size > 1024 * 1024 * 1024:
            print(f"Loading '{cohort_id}' via fast h5py CSR subsampling ({max_cells} cells)...")
            with h5py.File(preproc_h5ad, "r") as f:
                if "X" in f:
                    x_elem = f["X"]
                    if isinstance(x_elem, h5py.Group) and "indptr" in x_elem:
                        indptr = x_elem["indptr"][:]
                        total_cells = len(indptr) - 1
                        n_genes = int(x_elem.attrs.get("shape", [total_cells, 0])[1])

                        np.random.seed(seed)
                        n_pick = min(max_cells, total_cells)
                        keep_idx = np.sort(np.random.choice(total_cells, size=n_pick, replace=False))

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
                        n_pick = min(max_cells, total_cells)
                        np.random.seed(seed)
                        keep_idx = np.sort(np.random.choice(total_cells, size=n_pick, replace=False))
                        csr = csr_matrix(x_elem[keep_idx, :], dtype=np.float32)

                    try:
                        var_df = specs.read_elem(f["var"])
                        raw_genes = extract_gene_symbols_from_var(var_df)
                    except Exception:
                        raw_genes = [f"GENE_{i}" for i in range(csr.shape[1])]

                    adata_src = ad.AnnData(
                        X=csr,
                        obs=pd.DataFrame(index=[f"cell_{i}" for i in range(csr.shape[0])]),
                        var=pd.DataFrame(index=raw_genes),
                    )
                else:
                    return Failure(f"No 'X' dataset in {preproc_h5ad}")
        else:
            match load_dataset(cohort_id, auto_download=False, subsample_n=max_cells):
                case Failure(err):
                    return Failure(f"Failed to load '{cohort_id}': {err}")
                case Success(loaded_adata):
                    adata_src = loaded_adata

        raw_genes = extract_gene_symbols_from_var(adata_src.var)
        valid_gene_idx, filtered_genes = filter_confounding_genes(raw_genes)
        if len(filtered_genes) < 200:
            return Failure(f"Too few valid genes in {cohort_id}")

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

        sub_adata = ad.AnnData(
            X=csr_matrix(cell_by_gene, dtype=np.float32),
            obs=pd.DataFrame(
                {
                    "batch": cohort_id,
                    "dataset": cohort_id,
                },
                index=[f"{cohort_id}_{idx}" for idx in range(cell_by_gene.shape[0])],
            ),
            var=pd.DataFrame(index=filtered_genes),
        )
        return Success((sub_adata, filtered_genes))

    except Exception as exc:
        return Failure(f"Exception loading {cohort_id} for combination: {exc}")


def build_combined_reference(
    config: CombinedReferenceConfig,
    dataset_cache: dict[str, tuple[ad.AnnData, list[str]]],
) -> Result[Path, str]:
    """Combines single-cell cohorts, runs Harmony, clusters, and outputs Phi."""
    config.out_dir.mkdir(parents=True, exist_ok=True)
    metric_file = config.out_dir / f"reference_resolution_metrics_{config.combo_id}.parquet"
    primary_phi_path = config.out_dir / f"reference_phi_{config.combo_id}_res0.5.parquet"
    if metric_file.exists() and primary_phi_path.exists():
        print(f"[{config.combo_id}] Combined reference already exists. Skipping.")
        return Success(primary_phi_path)

    loaded_adatas: list[ad.AnnData] = []
    gene_sets: list[set[str]] = []

    for cid in config.cohort_ids:
        if cid not in dataset_cache:
            match load_cohort_for_combination(cid, config.cells_per_dataset, config.random_seed):
                case Failure(err):
                    print(f"Warning: Skipping {cid} in {config.combo_id}: {err}")
                    continue
                case Success(pair):
                    dataset_cache[cid] = pair

        ad_item, g_list = dataset_cache[cid]
        # Optional cell subsampling fraction titration
        if config.subsample_fraction < 0.999:
            n_keep = max(100, int(ad_item.n_obs * config.subsample_fraction))
            np.random.seed(config.random_seed)
            sub_idx = np.random.choice(ad_item.n_obs, size=n_keep, replace=False)
            ad_item = ad_item[sub_idx].copy()

        loaded_adatas.append(ad_item)
        gene_sets.append(set(g_list))

    if len(loaded_adatas) < 2:
        return Failure(f"Fewer than 2 cohorts available for combination {config.combo_id}")

    # Intersect shared genes
    common_genes = sorted(list(set.intersection(*gene_sets)))
    if len(common_genes) < 500:
        return Failure(f"Too few intersecting genes ({len(common_genes)}) across {config.combo_id}")

    # Harmonize gene matrices and concatenate
    subsetted: list[ad.AnnData] = []
    for ad_sub in loaded_adatas:
        ad_aligned = ad_sub[:, common_genes].copy()
        subsetted.append(ad_aligned)

    adata_comb = ad.concat(subsetted, join="inner")
    print(
        f"Combined '{config.combo_name}': {adata_comb.n_obs} cells across {len(loaded_adatas)} datasets, "
        f"{len(common_genes)} shared genes. Running Harmony batch correction..."
    )

    # Standard Scanpy HVG & PCA
    sc.pp.log1p(adata_comb)
    sc.pp.highly_variable_genes(adata_comb, n_top_genes=min(2500, len(common_genes)), subset=False)
    sc.tl.pca(adata_comb, n_comps=min(30, len(common_genes) - 1), random_state=config.random_seed)

    # Harmony batch correction
    try:
        import harmonypy  # type: ignore[import-untyped]

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
    except Exception as exc:
        print(f"Warning: Harmony integration encountered error ({exc}); falling back to standard PCA.")
        use_rep = "X_pca"

    sc.pp.neighbors(adata_comb, use_rep=use_rep, n_neighbors=15, random_state=config.random_seed)

    # Convert linear counts for cluster mean calculations
    lin_mat = adata_comb.X.toarray().T if issparse(adata_comb.X) else np.asarray(adata_comb.X.T, dtype=np.float64)
    if float(np.max(lin_mat)) < 35.0:
        lin_mat = np.expm1(np.clip(lin_mat, 0.0, 14.0))

    resolution_records: list[dict[str, object]] = []
    primary_phi_path = config.out_dir / f"reference_phi_{config.combo_id}_res0.5.parquet"

    for res in config.resolutions:
        col_name = f"leiden_{res}"
        sc.tl.leiden(
            adata_comb,
            resolution=res,
            key_added=col_name,
            flavor="igraph",
            directed=False,
            n_iterations=2,
        )

        cluster_series = adata_comb.obs[col_name].astype(str)
        unique_clusters: list[str] = sorted(cluster_series.unique().tolist())
        clusters: list[str] = cluster_series.tolist()

        mean_matrix = compute_cluster_means(lin_mat, clusters, unique_clusters)
        total_mean = mean_matrix.sum(axis=0)
        nonzero_mask = total_mean > 0
        final_genes = [g for g, v in zip(common_genes, nonzero_mask) if v]
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

        row_sums = sig_mean_matrix.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        phi_signature = sig_mean_matrix / row_sums

        cond_num = float(np.linalg.cond(phi_signature))
        print(
            f"  [{config.combo_name}] res={res}: {len(unique_clusters)} clusters, "
            f"{len(signature_genes)} signature genes, cond={cond_num:.2f}"
        )

        resolution_records.append(
            {
                "reference_type": config.combo_name,
                "reference_category": config.category,
                "resolution": float(res),
                "cluster_col": col_name,
                "n_clusters": len(unique_clusters),
                "n_signature_genes": len(signature_genes),
                "condition_number": cond_num,
                "n_datasets": len(loaded_adatas),
                "n_cells": adata_comb.n_obs,
                "subsample_fraction": config.subsample_fraction,
            }
        )

        # 1. Wide Phi
        wide_dict: dict[str, object] = {"cluster": list(unique_clusters)}
        for g_idx, g in enumerate(signature_genes):
            wide_dict[g] = [float(phi_signature[c_idx, g_idx]) for c_idx in range(len(unique_clusters))]
        df_phi_wide = pl.DataFrame(wide_dict)
        out_res_wide = config.out_dir / f"reference_phi_{config.combo_id}_res{res}.parquet"
        df_phi_wide.write_parquet(out_res_wide)

        # 2. Markers
        out_res_markers = config.out_dir / f"reference_marker_genes_{config.combo_id}_res{res}.parquet"
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
        df_phi_tidy.write_parquet(config.out_dir / f"reference_phi_tidy_{config.combo_id}_res{res}.parquet")

    # Metrics table
    df_metrics = pl.DataFrame(resolution_records)
    out_metrics = config.out_dir / f"reference_resolution_metrics_{config.combo_id}.parquet"
    df_metrics.write_parquet(out_metrics)
    print(f"Saved {config.combo_name} metrics to: {out_metrics}")

    return Success(primary_phi_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 1b: Build multi-dataset combined references via Harmony integration."
    )
    parser.add_argument(
        "--combo",
        type=str,
        default="all",
        help="Combination identifier (e.g. all_9_combined, nsclc_trio, melanoma_duo, or 'all')",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default="output/sc9_milo_deconv_validation/01_references",
        help="Directory to save output reference parquets",
    )
    parser.add_argument(
        "--cells-per-dataset",
        type=int,
        default=10000,
        help="Max cells to sample per cohort prior to integration",
    )
    parser.add_argument(
        "--resolutions",
        type=str,
        default="0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0",
        help="Comma-separated clustering resolutions",
    )
    args = parser.parse_args()

    res_tuple = tuple(float(r.strip()) for r in args.resolutions.split(",") if r.strip())
    out_path = Path(args.out_dir)

    all_combos: dict[str, tuple[str, str, tuple[str, ...]]] = {}
    for c_id, (name, cohort_list) in CRITERIA_COMBINATIONS.items():
        all_combos[c_id] = (name, "Criteria-Combined", cohort_list)
    for r_id, (name, cohort_list) in RANDOM_COMBINATIONS.items():
        all_combos[r_id] = (name, "Random-Combined", cohort_list)

    # Add cell titration series for the primary pan-cancer atlas
    for frac in (0.2, 0.4, 0.6, 0.8):
        frac_pct = int(frac * 100)
        c_id = f"all_9_combined_sub{frac_pct}"
        name = f"All-9-Milo-Combined ({frac_pct}% Cells)"
        all_combos[c_id] = (name, "Cell-Subsample", CRITERIA_COMBINATIONS["all_9_combined"][1])

    target_combos = (
        all_combos
        if args.combo == "all"
        else {args.combo: all_combos[args.combo]}
        if args.combo in all_combos
        else {}
    )

    if not target_combos:
        print(f"Unknown combination ID: {args.combo}. Choose from: {list(all_combos.keys())}", file=sys.stderr)
        sys.exit(1)

    dataset_cache: dict[str, tuple[ad.AnnData, list[str]]] = {}
    success_count = 0

    for combo_id, (combo_name, cat, cohorts) in target_combos.items():
        sub_frac = (
            float(int(combo_id.split("_sub")[1]) / 100.0)
            if "_sub" in combo_id
            else 1.0
        )
        config = CombinedReferenceConfig(
            combo_id=combo_id,
            combo_name=combo_name,
            category=cat,
            cohort_ids=cohorts,
            out_dir=out_path,
            resolutions=res_tuple,
            subsample_fraction=sub_frac,
            cells_per_dataset=args.cells_per_dataset,
        )
        print(f"\n=== Building Combined Reference: {config.combo_name} ({config.combo_id}) ===")
        match build_combined_reference(config, dataset_cache):
            case Success(p):
                success_count += 1
                print(f"  Successfully generated: {p}")
            case Failure(err):
                print(f"  Failed: {err}", file=sys.stderr)

    if success_count > 0:
        print(f"\nGenerated {success_count} combined reference configurations.")
        sys.exit(0)
    else:
        print("\nNo combined references were successfully generated.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
