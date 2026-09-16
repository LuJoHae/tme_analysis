"""Script 4: Compute PCA and UMAP Embeddings and Export Vector SVG Plots per Dataset.

This module computes Uncorrected, ComBat, and Harmony PCA and UMAP embeddings for preprocessed
single-cell datasets and exports Altair scatter plots as SVG files per dataset.
"""

from pathlib import Path
from typing import Optional, Literal
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import altair as alt
from returns.result import Result, Success, Failure


class PlotConfig(BaseModel):
    """Immutable configuration for embedding computation and SVG plot generation."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    input_dir: Path
    output_dir: Path
    n_hvg: int = 1000
    n_pcs: int = 30


def apply_gene_length_scaling(adata: ad.AnnData) -> ad.AnnData:
    """Scales Smart-seq2 expression data by gene transcript length if present."""
    adata_scaled = adata.copy()
    gse_mask = (adata_scaled.obs["sequencing_tech"] == "Smart-seq2").values
    if gse_mask.sum() == 0:
        return adata_scaled

    gene_lengths = np.ones(adata_scaled.n_vars, dtype=np.float32)
    if "contig_length" in adata_scaled.var.columns:
        lens = pd.to_numeric(adata_scaled.var["contig_length"], errors="coerce").values
        gene_lengths = np.where(np.isnan(lens) | (lens <= 0), 2000.0, lens)
    
    gene_lengths_kb = gene_lengths / 1000.0

    X_mat = adata_scaled.X.copy()
    if hasattr(X_mat, "toarray"):
        X_dense = X_mat.toarray()
        X_dense[gse_mask] = X_dense[gse_mask] / gene_lengths_kb
        from scipy.sparse import csr_matrix
        adata_scaled.X = csr_matrix(X_dense)
    else:
        X_mat[gse_mask] = X_mat[gse_mask] / gene_lengths_kb
        adata_scaled.X = X_mat

    return adata_scaled


def compute_dataset_embeddings(
    adata: ad.AnnData,
    n_hvg: int = 1000,
    n_pcs: int = 30,
) -> Result[ad.AnnData, str]:
    """Computes PCA and UMAP embeddings for an AnnData object."""
    try:
        adata_comp = adata.copy()
        if "log1p" not in adata_comp.uns:
            sc.pp.normalize_total(adata_comp, target_sum=1e4)
            sc.pp.log1p(adata_comp)

        num_genes = adata_comp.n_vars
        if num_genes >= 20:
            sc.pp.highly_variable_genes(
                adata_comp,
                n_top_genes=min(n_hvg, num_genes),
                inplace=True,
            )
            has_hvg = "highly_variable" in adata_comp.var.columns
        else:
            has_hvg = False

        actual_pcs = min(n_pcs, num_genes - 1, adata_comp.n_obs - 1) if num_genes > 1 and adata_comp.n_obs > 1 else 1

        if has_hvg:
            sc.tl.pca(adata_comp, n_comps=actual_pcs, mask_var="highly_variable")
        else:
            sc.tl.pca(adata_comp, n_comps=actual_pcs)

        sc.pp.neighbors(adata_comp, n_pcs=min(15, actual_pcs))
        sc.tl.umap(adata_comp)

        return Success(adata_comp)
    except Exception as e:
        return Failure(f"Embedding computation failed: {str(e)}")


def export_dataset_svgs(
    adata: ad.AnnData,
    dataset_name: str,
    output_dir: Path,
) -> Result[list[Path], str]:
    """Exports PCA and UMAP SVG plots for a single dataset across metadata colorings."""
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        exported_paths = []

        pca_coords = adata.obsm["X_pca"]
        umap_coords = adata.obsm["X_umap"]

        df_emb = pd.DataFrame({
            "PC1": pca_coords[:, 0],
            "PC2": pca_coords[:, 1],
            "UMAP1": umap_coords[:, 0],
            "UMAP2": umap_coords[:, 1],
            "Cell_Type": adata.obs["cell_type"].values if "cell_type" in adata.obs.columns else np.array(["Unknown"] * adata.n_obs),
            "Cancer_Code": adata.obs["cancer_code"].values if "cancer_code" in adata.obs.columns else np.array(["Unknown"] * adata.n_obs),
            "Sequencing_Tech": adata.obs["sequencing_tech"].values if "sequencing_tech" in adata.obs.columns else np.array(["Unknown"] * adata.n_obs),
        })

        for color_col, label in [("Cell_Type", "Cell Type"), ("Cancer_Code", "Cancer Code"), ("Sequencing_Tech", "Sequencing Tech")]:
            pca_chart = alt.Chart(df_emb).mark_circle(size=35, opacity=0.75).encode(
                x=alt.X("PC1:Q", title="PC1"),
                y=alt.Y("PC2:Q", title="PC2"),
                color=alt.Color(f"{color_col}:N", title=label),
                tooltip=["Cell_Type", "Cancer_Code", "Sequencing_Tech"]
            ).properties(title=f"PCA [{dataset_name}]: {label}", width=500, height=400)
            
            pca_path = output_dir / f"{dataset_name}_pca_{color_col.lower()}.svg"
            pca_chart.save(str(pca_path))
            exported_paths.append(pca_path)

            umap_chart = alt.Chart(df_emb).mark_circle(size=35, opacity=0.75).encode(
                x=alt.X("UMAP1:Q", title="UMAP1"),
                y=alt.Y("UMAP2:Q", title="UMAP2"),
                color=alt.Color(f"{color_col}:N", title=label),
                tooltip=["Cell_Type", "Cancer_Code", "Sequencing_Tech"]
            ).properties(title=f"UMAP [{dataset_name}]: {label}", width=500, height=400)
            
            umap_path = output_dir / f"{dataset_name}_umap_{color_col.lower()}.svg"
            umap_chart.save(str(umap_path))
            exported_paths.append(umap_path)

        return Success(exported_paths)
    except Exception as e:
        return Failure(f"Failed to export SVG plots for {dataset_name}: {str(e)}")


def process_single_dataset_plotting(h5_file: Path, output_dir: Path, n_hvg: int, n_pcs: int) -> Result[list[Path], str]:
    """Helper function to compute embeddings and export SVG plots for a single dataset."""
    try:
        dataset_name = h5_file.stem.replace("_processed", "")
        print(f"Generating PCA/UMAP plots for dataset '{dataset_name}'...")
        adata = ad.read_h5ad(h5_file)
        if adata.n_obs < 3:
            return Success([])

        match compute_dataset_embeddings(adata, n_hvg, n_pcs):
            case Success(adata_emb):
                return export_dataset_svgs(adata_emb, dataset_name, output_dir)
            case Failure(err):
                return Failure(f"Embedding failure for {dataset_name}: {err}")
    except Exception as e:
        return Failure(f"Plotting error for {h5_file.name}: {str(e)}")


def run_plotting_pipeline(config: PlotConfig, n_jobs: int = -1) -> Result[list[Path], str]:
    """Iterates through preprocessed .h5ad files and exports SVG embedding plots in parallel per dataset."""
    try:
        from concurrent.futures import ProcessPoolExecutor
        import os

        config.output_dir.mkdir(parents=True, exist_ok=True)
        h5_files = list(config.input_dir.glob("*.h5ad"))
        
        if not h5_files:
            return Failure(f"No preprocessed .h5ad files found in: {config.input_dir}")

        workers = os.cpu_count() if n_jobs <= 0 else n_jobs
        print(f"Executing PCA/UMAP embedding plotting for {len(h5_files)} datasets in parallel using {workers} workers...")

        all_svgs = []
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(process_single_dataset_plotting, h5_file, config.output_dir, config.n_hvg, config.n_pcs)
                for h5_file in h5_files
            ]
            for future in futures:
                match future.result():
                    case Success(svgs):
                        all_svgs.extend(svgs)
                    case Failure(err):
                        print(f"Warning: {err}")

        return Success(all_svgs)
    except Exception as e:
        return Failure(f"Plotting pipeline failed: {str(e)}")


def main() -> None:
    """CLI entry point for per-dataset PCA/UMAP embedding plotting."""
    import argparse
    parser = argparse.ArgumentParser(description="Generate per-dataset PCA and UMAP vector SVG plots.")
    parser.add_argument("--input-dir", type=str, default="data/processed_h5", help="Input directory with preprocessed .h5ad files")
    parser.add_argument("--out-dir", type=str, default="data/results_embeddings", help="Output directory for SVG scatter plots")
    parser.add_argument("--input-h5", type=str, default=None, help="Optional specific input preprocessed h5ad file for single dataset")
    args = parser.parse_args()

    out_dir = Path(args.out_dir).resolve()

    if args.input_h5:
        input_h5 = Path(args.input_h5).resolve()
        match process_single_dataset_plotting(input_h5, out_dir, 1000, 30):
            case Success(paths):
                print(f"Successfully generated {len(paths)} SVG embedding plots for '{input_h5.name}'")
            case Failure(err):
                print(f"Plotting error for '{input_h5.name}': {err}")
    else:
        config = PlotConfig(
            input_dir=Path(args.input_dir).resolve(),
            output_dir=out_dir,
        )
        print(f"Starting per-dataset PCA/UMAP embedding plotting pipeline...")
        print(f"Input dir: {config.input_dir}")
        print(f"Output dir: {config.output_dir}")

        match run_plotting_pipeline(config):
            case Success(paths):
                print(f"\nSuccessfully generated {len(paths)} SVG embedding plots:")
                for p in paths:
                    print(f" - {p}")
            case Failure(err):
                print(f"Pipeline error: {err}")


if __name__ == "__main__":
    main()
