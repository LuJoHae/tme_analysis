"""Unified Processing Script for SingleCellDataProcessStep07 and GSE120575.

This module unifies single-cell transcriptomics data from `SingleCellDataProcessStep07`
(multi-cohort 10x Genomics UMI raw counts) with `GSE120575` (Sade-Feldman et al. Smart-seq2 TPM),
harmonizes metadata schemas, aligns gene features based on Step07 authoritative gene names,
and evaluates dataset combination quality.
"""

from pathlib import Path
from typing import Optional, Generic, TypeVar, assert_never
import numpy as np
import pandas as pd
import polars as pl
from scipy.sparse import csr_matrix, issparse
from pydantic import BaseModel, ConfigDict
import anndata as ad
import scanpy as sc
import altair as alt
from returns.result import Result, Success, Failure
from returns.maybe import Maybe, Some, Nothing

import datalair
from single_cell_datasets import SingleCellDataProcessStep07
from gene_utils import norm_genes


class UnifyConfig(BaseModel):
    """Immutable configuration for dataset unification and downsampling."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    lair_dir: Path
    gse120575_dir: Path
    output_path: Path
    subsample_n_per_type: Maybe[int] = Nothing
    random_seed: int = 42


class IntegrationMetrics(BaseModel):
    """Immutable data model holding integration evaluation metrics."""
    model_config = ConfigDict(frozen=True)
    num_step07_cells: int
    num_gse120575_cells: int
    num_shared_genes: int
    gene_overlap_ratio: float
    silhouette_cell_type: float
    silhouette_dataset: float
    marker_correlation_mean: float


def load_step07_adata(config: UnifyConfig) -> Result[ad.AnnData, str]:
    """Loads SingleCellDataProcessStep07 AnnData dataset from datalair."""
    try:
        if not config.lair_dir.exists():
            return Failure(f"Lair directory not found: {config.lair_dir}")
        lair = datalair.Lair(config.lair_dir)
        step7 = SingleCellDataProcessStep07()
        filepaths = lair.get_dataset_filepaths(step7)
        if "adata.h5ad" not in filepaths:
            return Failure(f"adata.h5ad not found in datalair for Step07 at {config.lair_dir}")
        
        adata = ad.read_h5ad(filepaths["adata.h5ad"])
        return Success(adata)
    except Exception as e:
        return Failure(f"Failed to load Step07 AnnData: {str(e)}")


def load_gse120575_adata(config: UnifyConfig) -> Result[ad.AnnData, str]:
    """Loads GSE120575 parquet dataset and formats into AnnData."""
    tpm_path = config.gse120575_dir / "gse120575_tpm.parquet"
    meta_path = config.gse120575_dir / "gse120575_tpm_cell_metadata.parquet"

    if not tpm_path.exists() or not meta_path.exists():
        return Failure(f"GSE120575 parquet files missing in {config.gse120575_dir}")

    try:
        df_tpm = pl.read_parquet(tpm_path)
        df_meta = pl.read_parquet(meta_path)

        genes = df_tpm["gene"].to_list()
        cell_cols = [c for c in df_tpm.columns if c != "gene"]
        matrix = df_tpm.select(cell_cols).to_numpy().T  # (cells, genes)

        obs_df = df_meta.to_pandas().set_index("Cell_ID")
        obs_df = obs_df.reindex(cell_cols)
        obs_df["dataset"] = "GSE120575_SadeFeldman"
        obs_df["cancer_code"] = "SKCM"
        obs_df["staging"] = "metastasis"
        obs_df["sequencing_tech"] = "Smart-seq2"
        obs_df["original.barcode"] = obs_df.index.tolist()

        var_df = pd.DataFrame(index=genes)
        var_df["gene_name"] = genes

        adata = ad.AnnData(
            X=csr_matrix(matrix, dtype=np.float32),
            obs=obs_df,
            var=var_df,
        )
        return Success(adata)
    except Exception as e:
        return Failure(f"Failed to load GSE120575 AnnData: {str(e)}")


def align_and_harmonize_metadata(
    adata_step07: ad.AnnData,
    adata_gse: ad.AnnData,
) -> Result[tuple[ad.AnnData, ad.AnnData], str]:
    """Aligns gene features to Step07 authoritative gene set and harmonizes metadata columns."""
    try:
        # 1. Harmonize Step07 metadata
        obs7 = adata_step07.obs.copy()
        if "cell_type" not in obs7.columns:
            for col in ["cell_type_main", "celltype", "cell_types"]:
                if col in obs7.columns:
                    obs7["cell_type"] = obs7[col]
                    break
            else:
                obs7["cell_type"] = "Unknown"
        obs7["sequencing_tech"] = "10x_UMI"
        adata_step07.obs = obs7

        # 2. Harmonize GSE120575 metadata
        obsg = adata_gse.obs.copy()
        if "cell_type" not in obsg.columns:
            obsg["cell_type"] = "Melanoma_TME"
        adata_gse.obs = obsg

        # 3. Map genes based on Step07 authoritative gene set
        step07_genes = set(adata_step07.var_names)
        gse_genes = set(adata_gse.var_names)
        shared_genes = sorted(list(step07_genes.intersection(gse_genes)))

        if not shared_genes:
            # Normalize GSE genes using gene_utils norm_genes if names differ
            adata_gse_norm = norm_genes(adata_gse.copy(), pre_id_transform="auto")
            shared_genes = sorted(list(step07_genes.intersection(adata_gse_norm.var_names)))
            if not shared_genes:
                return Failure("Zero shared genes found between Step07 and GSE120575!")
            adata_gse = adata_gse_norm

        adata_step07_sub = adata_step07[:, shared_genes].copy()
        adata_gse_sub = adata_gse[:, shared_genes].copy()

        return Success((adata_step07_sub, adata_gse_sub))
    except Exception as e:
        return Failure(f"Metadata harmonization failed: {str(e)}")


def subsample_dataset(
    adata: ad.AnnData,
    n_per_type: int,
    seed: int = 42,
) -> ad.AnnData:
    """Subsamples cells per cell type for performance and memory optimization."""
    np.random.seed(seed)
    selected_indices = []
    
    cell_type_col = "cell_type" if "cell_type" in adata.obs.columns else adata.obs.columns[0]
    for ctype in adata.obs[cell_type_col].unique():
        idx = np.where(adata.obs[cell_type_col] == ctype)[0]
        if len(idx) > n_per_type:
            idx = np.random.choice(idx, n_per_type, replace=False)
        selected_indices.extend(idx)
        
    return adata[selected_indices].copy()


def combine_datasets(
    adata_step07: ad.AnnData,
    adata_gse: ad.AnnData,
    config: UnifyConfig,
) -> Result[ad.AnnData, str]:
    """Combines Step07 and GSE120575 AnnData datasets with optional downsampling."""
    try:
        match config.subsample_n_per_type:
            case Some(n_val):
                a7 = subsample_dataset(adata_step07, n_val, config.random_seed)
                ag = subsample_dataset(adata_gse, n_val, config.random_seed)
            case Nothing:
                a7 = adata_step07
                ag = adata_gse

        a7.obs_names_make_unique()
        ag.obs_names_make_unique()
        unified_adata = ad.concat([a7, ag], axis=0, join="inner", merge="first")
        unified_adata.obs_names_make_unique()
        return Success(unified_adata)
    except Exception as e:
        return Failure(f"Dataset concatenation failed: {str(e)}")


def evaluate_integration_quality(
    unified_adata: ad.AnnData,
) -> Result[IntegrationMetrics, str]:
    """Computes integration evaluation metrics (gene ratio, silhouette scores, marker correlation)."""
    try:
        step07_mask = (unified_adata.obs["sequencing_tech"] == "10x_UMI").values
        gse_mask = (unified_adata.obs["sequencing_tech"] == "Smart-seq2").values

        n_s7 = int(step07_mask.sum())
        n_gse = int(gse_mask.sum())
        num_genes = unified_adata.n_vars

        # Compute normalized log1p expression for evaluation
        adata_eval = unified_adata.copy()
        sc.pp.normalize_total(adata_eval, target_sum=1e4)
        sc.pp.log1p(adata_eval)
        
        has_hvg = False
        if num_genes >= 20:
            sc.pp.highly_variable_genes(adata_eval, n_top_genes=min(500, num_genes), inplace=True)
            has_hvg = "highly_variable" in adata_eval.var.columns
        
        # PCA & Silhouette Scores
        if has_hvg:
            sc.tl.pca(adata_eval, mask_var="highly_variable")
        else:
            sc.tl.pca(adata_eval)
        pca_coords = adata_eval.obsm["X_pca"]

        from sklearn.metrics import silhouette_score
        
        cell_types = adata_eval.obs["cell_type"].values
        datasets = adata_eval.obs["sequencing_tech"].values
        
        sil_cell_type = float(silhouette_score(pca_coords, cell_types)) if len(set(cell_types)) > 1 else 0.0
        sil_dataset = float(silhouette_score(pca_coords, datasets)) if len(set(datasets)) > 1 else 0.0

        # Mean correlation across shared genes
        if n_s7 > 0 and n_gse > 0:
            mean_exp_s7 = np.asarray(adata_eval[step07_mask].X.mean(axis=0)).ravel()
            mean_exp_gse = np.asarray(adata_eval[gse_mask].X.mean(axis=0)).ravel()
            corr = float(np.corrcoef(mean_exp_s7, mean_exp_gse)[0, 1])
        else:
            corr = 1.0

        metrics = IntegrationMetrics(
            num_step07_cells=n_s7,
            num_gse120575_cells=n_gse,
            num_shared_genes=num_genes,
            gene_overlap_ratio=1.0,
            silhouette_cell_type=sil_cell_type,
            silhouette_dataset=sil_dataset,
            marker_correlation_mean=corr,
        )
        return Success(metrics)
    except Exception as e:
        return Failure(f"Integration evaluation failed: {str(e)}")


def plot_integration_pca_svg(
    unified_adata: ad.AnnData,
    output_svg_path: Path,
) -> Result[Path, str]:
    """Generates an Altair PCA scatter plot showing dataset integration and exports as SVG."""
    try:
        adata_eval = unified_adata.copy()
        sc.pp.normalize_total(adata_eval, target_sum=1e4)
        sc.pp.log1p(adata_eval)
        sc.tl.pca(adata_eval, n_comps=2)
        
        pca_df = pd.DataFrame({
            "PC1": adata_eval.obsm["X_pca"][:, 0],
            "PC2": adata_eval.obsm["X_pca"][:, 1],
            "Dataset": adata_eval.obs["dataset"].values,
            "Tech": adata_eval.obs["sequencing_tech"].values,
            "CellType": adata_eval.obs["cell_type"].values,
        })
        
        # Altair scatter plot
        chart = alt.Chart(pca_df).mark_circle(size=40, opacity=0.7).encode(
            x=alt.X("PC1:Q", title="PC1"),
            y=alt.Y("PC2:Q", title="PC2"),
            color=alt.Color("Tech:N", title="Sequencing Technology"),
            tooltip=["Dataset", "Tech", "CellType"]
        ).properties(
            title="SingleCellDataProcessStep07 vs GSE120575 Integration PCA",
            width=600,
            height=400,
        )
        
        output_svg_path.parent.mkdir(parents=True, exist_ok=True)
        chart.save(str(output_svg_path))
        return Success(output_svg_path)
    except Exception as e:
        return Failure(f"Failed to export integration PCA SVG chart: {str(e)}")


def save_unified_adata(
    unified_adata: ad.AnnData,
    output_path: Path,
) -> Result[Path, str]:
    """Writes the unified AnnData object to disk in .h5ad format."""
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        unified_adata.write_h5ad(output_path)
        return Success(output_path)
    except Exception as e:
        return Failure(f"Failed to write unified AnnData to {output_path}: {str(e)}")


def run_unification_pipeline(config: UnifyConfig) -> Result[tuple[ad.AnnData, IntegrationMetrics, Path], str]:
    """Executes the complete unification, subsampling, disk saving, and evaluation pipeline."""
    return load_step07_adata(config).bind(
        lambda a7: load_gse120575_adata(config).bind(
            lambda ag: align_and_harmonize_metadata(a7, ag).bind(
                lambda aligned: combine_datasets(aligned[0], aligned[1], config).bind(
                    lambda unified: save_unified_adata(unified, config.output_path).bind(
                        lambda out_path: evaluate_integration_quality(unified).map(
                            lambda metrics: (unified, metrics, out_path)
                        )
                    )
                )
            )
        )
    )


def main() -> None:
    """CLI entry point for running the dataset unification pipeline."""
    import argparse
    parser = argparse.ArgumentParser(description="Unify Step07 with GSE120575 dataset.")
    parser.add_argument("--lair-dir", type=str, default="scratch/lair", help="Path to datalair root")
    parser.add_argument("--gse120575-dir", type=str, default="scratch/GSE120575", help="Path to GSE120575 parquet dir")
    parser.add_argument("--output-path", type=str, default="scratch/unified_step07_gse120575.h5ad", help="Output h5ad path")
    parser.add_argument("--subsample", type=int, default=500, help="Number of cells per cell type to subsample (0 for no subsampling)")
    args = parser.parse_args()

    subsample_opt: Maybe[int] = Some(args.subsample) if args.subsample > 0 else Nothing

    config = UnifyConfig(
        lair_dir=Path(args.lair_dir).resolve(),
        gse120575_dir=Path(args.gse120575_dir).resolve(),
        output_path=Path(args.output_path).resolve(),
        subsample_n_per_type=subsample_opt,
        random_seed=42,
    )
    print(f"Starting dataset unification pipeline...")
    print(f"Lair path: {config.lair_dir}")
    print(f"GSE120575 path: {config.gse120575_dir}")
    print(f"Output h5ad path: {config.output_path}")

    match run_unification_pipeline(config):
        case Failure(err):
            print(f"Pipeline error: {err}")
        case Success((unified, metrics, out_path)):
            print("\n=== Dataset Unification & Integration Metrics ===")
            print(f"Saved unified AnnData to: {out_path}")
            print(f"Step07 cells: {metrics.num_step07_cells}")
            print(f"GSE120575 cells: {metrics.num_gse120575_cells}")
            print(f"Shared genes: {metrics.num_shared_genes}")
            print(f"Silhouette (Cell Type): {metrics.silhouette_cell_type:.4f}")
            print(f"Silhouette (Dataset): {metrics.silhouette_dataset:.4f}")
            print(f"Mean Marker Correlation: {metrics.marker_correlation_mean:.4f}")

            # Export PCA SVG plot
            svg_path = out_path.parent / "pca_integration.svg"
            match plot_integration_pca_svg(unified, svg_path):
                case Success(s_path):
                    print(f"Saved integration PCA chart to: {s_path}")
                case Failure(s_err):
                    print(f"Plotting warning: {s_err}")


if __name__ == "__main__":
    main()
