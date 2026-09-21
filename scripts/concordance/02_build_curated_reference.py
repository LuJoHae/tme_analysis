#!/usr/bin/env python3
"""
Step 2: Build Curated Reference for BayesPrism Deconvolution.
Constructs unnormalized mean GEP profiles, flags collinear cell states,
calculates cell-type-specific mRNA scaling factors, and builds the cell hierarchy.
Strict functional Python with returns, Pydantic, and Polars.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

import anndata as ad  # type: ignore
import numpy as np
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
import polars as pl
from scipy import sparse  # type: ignore


class ReferenceBuildConfig(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    adata_path: Path
    cluster_col: str
    lineage_col: Maybe[str]
    malignant_label: str
    collinearity_threshold: float
    out_dir: Path


def parse_args() -> ReferenceBuildConfig:
    parser = argparse.ArgumentParser(
        description="Build curated single-cell reference and mRNA scaling factors for BayesPrism."
    )
    parser.add_argument(
        "--adata",
        type=Path,
        required=True,
        help="Path to single-cell AnnData .h5ad",
    )
    parser.add_argument(
        "--cluster-col",
        type=str,
        default="celltypist_leiden_0.5",
        help="Obs column name for cell state / cluster",
    )
    parser.add_argument(
        "--lineage-col",
        type=str,
        default=None,
        help="Obs column name for broad lineage (cell type). If not provided, inferred from cluster-col",
    )
    parser.add_argument(
        "--malignant-label",
        type=str,
        default="Malignant",
        help="Label denoting malignant/tumor cells for BayesPrism prior",
    )
    parser.add_argument(
        "--collinearity-threshold",
        type=float,
        default=0.85,
        help="Correlation threshold above which cell states are flagged for collinearity",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Output directory for reference parquets",
    )
    args = parser.parse_args()
    return ReferenceBuildConfig(
        adata_path=args.adata,
        cluster_col=args.cluster_col,
        lineage_col=Some(args.lineage_col) if args.lineage_col else Nothing,
        malignant_label=args.malignant_label,
        collinearity_threshold=args.collinearity_threshold,
        out_dir=args.out_dir,
    )


def filter_uninformative_genes(gene_names: Sequence[str]) -> list[bool]:
    """Flag uninformative/confounding genes: mitochondrial, ribosomal, and heat-shock."""
    keep_mask: list[bool] = []
    for g in gene_names:
        upper_g = g.upper()
        if (
            upper_g.startswith("MT-")
            or upper_g.startswith("RPS")
            or upper_g.startswith("RPL")
            or upper_g.startswith("HSP")
            or upper_g.startswith("DNAJ")
        ):
            keep_mask.append(False)
        else:
            keep_mask.append(True)
    return keep_mask


def compute_mrna_scaling(
    adata: ad.AnnData,
    cluster_col: str,
) -> pl.DataFrame:
    """Calculate mean UMI count (mRNA mass proxy) per cell state."""
    # Compute library size per cell
    counts = adata.raw.X if adata.raw is not None else adata.X
    lib_sizes = np.array(counts.sum(axis=1)).flatten()

    obs_df = pl.from_pandas(adata.obs.reset_index())
    obs_with_depth = obs_df.with_columns(pl.Series("total_umi", lib_sizes))

    scaling_df = (
        obs_with_depth.group_by(cluster_col)
        .agg([
            pl.col("total_umi").mean().alias("mean_umi_per_cell"),
            pl.col("total_umi").median().alias("median_umi_per_cell"),
            pl.len().alias("n_cells"),
        ])
        .rename({cluster_col: "cell_state"})
    )

    grand_mean = scaling_df["mean_umi_per_cell"].mean()
    return scaling_df.with_columns(
        (pl.col("mean_umi_per_cell") / grand_mean).alias("relative_rna_content")
    )


def compute_reference_gep(
    adata: ad.AnnData,
    cluster_col: str,
    gene_mask: list[bool],
) -> tuple[pl.DataFrame, np.ndarray, list[str], list[str]]:
    """Compute mean expression profile for each cell state across informative genes."""
    counts = adata.raw.X if adata.raw is not None else adata.X
    all_genes = list(adata.raw.var_names if adata.raw is not None else adata.var_names)

    filtered_genes = [g for g, keep in zip(all_genes, gene_mask, strict=True) if keep]
    gene_indices = [i for i, keep in enumerate(gene_mask) if keep]

    filtered_counts = counts[:, gene_indices]
    if sparse.issparse(filtered_counts):
        filtered_counts = filtered_counts.tocsr()

    obs_clusters = adata.obs[cluster_col].astype(str).to_numpy()
    unique_states = sorted(list(set(obs_clusters)))

    mean_profiles = np.zeros((len(unique_states), len(filtered_genes)), dtype=np.float64)

    for i, state in enumerate(unique_states):
        cell_mask = obs_clusters == state
        sub_counts = filtered_counts[cell_mask]
        mean_profiles[i, :] = np.array(sub_counts.mean(axis=0)).flatten()

    phi_dict: dict[str, list[float] | list[str]] = {"cell_state": unique_states}
    for j, g_name in enumerate(filtered_genes):
        phi_dict[g_name] = mean_profiles[:, j].tolist()

    phi_df = pl.DataFrame(phi_dict)
    return phi_df, mean_profiles, unique_states, filtered_genes


def evaluate_collinearity(
    mean_profiles: np.ndarray,
    states: list[str],
    threshold: float,
) -> pl.DataFrame:
    """Detect pairs of cell states with high cross-correlation that can cause signal leakage."""
    # Standardize rows
    means = mean_profiles.mean(axis=1, keepdims=True)
    stds = mean_profiles.std(axis=1, keepdims=True) + 1e-12
    norm_profiles = (mean_profiles - means) / stds

    corr_mat = np.dot(norm_profiles, norm_profiles.T) / mean_profiles.shape[1]

    pairs: list[dict[str, object]] = []
    for i in range(len(states)):
        for j in range(i + 1, len(states)):
            r_val = float(corr_mat[i, j])
            pairs.append({
                "state_a": states[i],
                "state_b": states[j],
                "correlation": r_val,
                "is_collinear": r_val >= threshold,
            })

    return pl.DataFrame(pairs).sort("correlation", descending=True)


def build_hierarchy_table(
    adata: ad.AnnData,
    cluster_col: str,
    lineage_col_opt: Maybe[str],
    malignant_label: str,
) -> pl.DataFrame:
    """Build two-tier hierarchy mapping table."""
    obs_df = pl.from_pandas(adata.obs.reset_index())

    lineage_col = lineage_col_opt.value_or(cluster_col)
    if lineage_col not in obs_df.columns:
        lineage_col = cluster_col

    mapping = (
        obs_df.select([cluster_col, lineage_col])
        .unique()
        .rename({cluster_col: "cell_state", lineage_col: "cell_type"})
    )

    # Flag malignant
    return mapping.with_columns(
        pl.when(
            pl.col("cell_state").str.to_lowercase().str.contains(malignant_label.lower())
            | pl.col("cell_type").str.to_lowercase().str.contains(malignant_label.lower())
        )
        .then(True)
        .otherwise(False)
        .alias("is_malignant")
    )


def run_pipeline(config: ReferenceBuildConfig) -> Result[None, str]:
    """Execute Step 2 pure reference preparation pipeline."""
    if not config.adata_path.exists():
        return Failure(f"AnnData file does not exist: {config.adata_path}")

    try:
        adata = ad.read_h5ad(config.adata_path, backed=False)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to read AnnData: {exc}")

    if config.cluster_col not in adata.obs.columns:
        return Failure(f"Cluster column '{config.cluster_col}' not found in AnnData.obs")

    # 1. Filter genes
    all_genes = list(adata.raw.var_names if adata.raw is not None else adata.var_names)
    gene_mask = filter_uninformative_genes(all_genes)
    retained_genes = sum(gene_mask)

    # 2. mRNA scaling factors
    mrna_df = compute_mrna_scaling(adata, config.cluster_col)

    # 3. Mean GEP profiles
    phi_df, mean_profiles, unique_states, _ = compute_reference_gep(adata, config.cluster_col, gene_mask)

    # 4. Collinearity analysis
    collinearity_df = evaluate_collinearity(mean_profiles, unique_states, config.collinearity_threshold)

    # 5. Cell hierarchy
    hierarchy_df = build_hierarchy_table(
        adata, config.cluster_col, config.lineage_col, config.malignant_label
    )

    # Save parquets
    config.out_dir.mkdir(parents=True, exist_ok=True)
    phi_df.write_parquet(config.out_dir / "curated_reference_phi.parquet")
    mrna_df.write_parquet(config.out_dir / "mrna_scaling_factors.parquet")
    collinearity_df.write_parquet(config.out_dir / "reference_collinearity.parquet")
    hierarchy_df.write_parquet(config.out_dir / "cell_hierarchy.parquet")

    high_collinear_count = collinearity_df.filter(pl.col("is_collinear")).height
    print(
        f"[INFO] Retained {retained_genes}/{len(all_genes)} genes. "
        f"Evaluated {len(unique_states)} cell states. "
        f"Flagged {high_collinear_count} collinear pairs (r >= {config.collinearity_threshold})."
    )

    return Success(None)


def main() -> None:
    config = parse_args()
    match run_pipeline(config):
        case Success(_):
            print(f"[SUCCESS] Curated reference files written to {config.out_dir}")
            sys.exit(0)
        case Failure(err):
            print(f"[ERROR] {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
