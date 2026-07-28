import sys
import argparse
from pathlib import Path
import polars as pl
import anndata as ad
import scanpy as sc
import altair as alt
import scipy.sparse as sp
import numpy as np
import math
from returns.result import Result, Success, Failure
from returns.methods import bind

# Allow large datasets in Altair
alt.data_transformers.disable_max_rows()

def load_data(tpm_path: Path, meta_path: Path) -> Result[tuple[pl.DataFrame, pl.DataFrame], str]:
    try:
        tpm_df = pl.read_parquet(tpm_path)
        meta_df = pl.read_parquet(meta_path)
        return Success((tpm_df, meta_df))
    except Exception as e:
        return Failure(f"Failed to load data: {e}")

def create_anndata(tpm_df: pl.DataFrame, meta_df: pl.DataFrame) -> Result[ad.AnnData, str]:
    try:
        genes = tpm_df["gene"].to_numpy()
        expr_df = tpm_df.drop("gene")
        cell_ids = expr_df.columns
        
        # Transpose expression matrix (cells x genes) and convert to sparse
        X_sparse = sp.csr_matrix(expr_df.to_numpy().T)
        
        # Prepare metadata dataframe
        obs_df = meta_df.to_pandas()
        obs_df.set_index("cell_id", inplace=True)
        
        # Create AnnData
        adata = ad.AnnData(
            X=X_sparse,
            obs=obs_df.reindex(cell_ids),
            var=pl.DataFrame({"gene": genes}).to_pandas().set_index("gene")
        )
        
        # Filter cells to only those that exist in our cleaned metadata
        valid_cells = ~adata.obs["melanoma-sample"].isna()
        adata = adata[valid_cells].copy()
        
        return Success(adata)
    except Exception as e:
        return Failure(f"Failed to create AnnData: {e}")

def process_scanpy(adata: ad.AnnData) -> Result[ad.AnnData, str]:
    try:
        # Filter exact duplicate cell indices (if any)
        _, unique_indices = np.unique(adata.obs_names, return_index=True)
        adata = adata[unique_indices].copy()
        
        # QC Filtering
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        
        # Log1p transformation (TPM is already size-factor normalized, just log it)
        sc.pp.log1p(adata)
        
        # PCA
        sc.tl.pca(adata, svd_solver='arpack')
        
        # Neighborhood graph
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        
        # UMAP
        sc.tl.umap(adata)
        
        return Success(adata)
    except Exception as e:
        return Failure(f"Failed in scanpy processing: {e}")

def create_umap_plots(adata: ad.AnnData) -> Result[alt.Chart, str]:
    try:
        umap_coords = adata.obsm['X_umap']
        
        plot_df = pl.DataFrame({
            "cell_id": adata.obs_names,
            "UMAP1": umap_coords[:, 0],
            "UMAP2": umap_coords[:, 1],
        })
        
        obs_df = pl.from_pandas(adata.obs.reset_index())
        if "index" in obs_df.columns:
            obs_df = obs_df.rename({"index": "cell_id"})
            
        plot_df = plot_df.join(obs_df, on="cell_id")
        
        exclude_cols = {"cell_id", "plate-row", "plate-col"}
        
        charts = []
        for col in adata.obs.columns:
            if col in exclude_cols or plot_df[col].null_count() == plot_df.height:
                continue
                
            chart = alt.Chart(plot_df).mark_circle(size=15, opacity=0.8).encode(
                x=alt.X("UMAP1:Q", title="UMAP 1"),
                y=alt.Y("UMAP2:Q", title="UMAP 2"),
                color=alt.Color(f"{col}:N", title=col),
                tooltip=["cell_id", col]
            ).properties(
                title=f"UMAP colored by {col}",
                width=350,
                height=350
            )
            charts.append(chart)
            
        if not charts:
            return Failure("No valid metadata columns found to plot.")
            
        # Layout in a grid (2 columns wide)
        cols = 2
        h_charts = []
        for i in range(0, len(charts), cols):
            row_charts = charts[i:i+cols]
            h_charts.append(alt.hconcat(*row_charts))
            
        final_chart = alt.vconcat(*h_charts).resolve_scale(color='independent')
        
        return Success(final_chart)
    except Exception as e:
        return Failure(f"Failed to create plots: {e}")

def save_plots(chart: alt.Chart, out_svg: Path) -> Result[bool, str]:
    try:
        out_svg.parent.mkdir(parents=True, exist_ok=True)
        chart.save(str(out_svg))
        return Success(True)
    except Exception as e:
        return Failure(f"Failed to save plots: {e}")

def run_pipeline(tpm_path: Path, meta_path: Path, out_svg: Path) -> Result[bool, str]:
    return load_data(tpm_path, meta_path).bind(
        lambda dfs: create_anndata(dfs[0], dfs[1]).bind(
            lambda adata: process_scanpy(adata).bind(
                lambda processed_adata: create_umap_plots(processed_adata).bind(
                    lambda chart: save_plots(chart, out_svg)
                )
            )
        )
    )

def main():
    parser = argparse.ArgumentParser(description="Analyze GSE120575 data using Scanpy")
    parser.add_argument("--tpm", required=True, help="Path to TPM parquet")
    parser.add_argument("--meta", required=True, help="Path to parsed metadata parquet")
    parser.add_argument("--out-svg", required=True, help="Output plot SVG")
    args = parser.parse_args()

    match run_pipeline(Path(args.tpm), Path(args.meta), Path(args.out_svg)):
        case Success(_):
            print("Successfully completed single-cell analysis and saved UMAP plots.")
            sys.exit(0)
        case Failure(err):
            print(f"Error processing data: {err}")
            sys.exit(1)

if __name__ == "__main__":
    main()
