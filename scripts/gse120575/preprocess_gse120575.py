import sys
import argparse
from pathlib import Path
import polars as pl  # type: ignore
import anndata as ad  # type: ignore
import scanpy as sc  # type: ignore
import scipy.sparse as sp  # type: ignore
import numpy as np  # type: ignore
from typing import Generator, Any
from returns.result import Result, Success, Failure  # type: ignore
from returns.decorators import do  # type: ignore

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
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Calculate highly variable genes
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata = adata[:, adata.var.highly_variable].copy()
        
        # PCA
        sc.tl.pca(adata, svd_solver='arpack')
        
        # Neighborhood graph
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        
        # UMAP
        sc.tl.umap(adata)
        
        return Success(adata)
    except Exception as e:
        return Failure(f"Failed in scanpy processing: {e}")

def save_anndata(adata: ad.AnnData, out_path: Path) -> Result[bool, str]:
    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(out_path)
        return Success(True)
    except Exception as e:
        return Failure(f"Failed to save AnnData: {e}")

@do(Result[bool, str])  # type: ignore
def run_pipeline(tpm_path: Path, meta_path: Path, out_path: Path) -> Generator[Any, Any, bool]:
    dfs = yield load_data(tpm_path, meta_path)
    adata = yield create_anndata(dfs[0], dfs[1])
    processed_adata = yield process_scanpy(adata)
    success = yield save_anndata(processed_adata, out_path)
    return bool(success)

def main() -> None:
    parser = argparse.ArgumentParser(description="Preprocess GSE120575 data using Scanpy")
    parser.add_argument("--tpm", required=True, help="Path to TPM parquet")
    parser.add_argument("--meta", required=True, help="Path to parsed metadata parquet")
    parser.add_argument("--out-h5ad", required=True, help="Output processed AnnData h5ad")
    args = parser.parse_args()

    match run_pipeline(Path(args.tpm), Path(args.meta), Path(args.out_h5ad)):
        case Success(_):
            print("Successfully completed preprocessing and saved AnnData.")
            sys.exit(0)
        case Failure(err):
            print(f"Error processing data: {err}")
            sys.exit(1)

if __name__ == "__main__":
    main()
