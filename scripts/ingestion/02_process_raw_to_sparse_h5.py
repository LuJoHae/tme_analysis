"""Script 2: Process Raw GEO Data into HDF5 (.h5ad) Files with Sparse CSR Encoding.

Inspired by the single_cell_datasets package (SingleCellDataProcessStep01 & Step04),
this script reads raw single-cell expression matrices, converts them to scipy.sparse.csr_matrix
for memory efficiency, normalizes gene IDs using gene_utils.norm_genes, and writes sparse .h5ad files.
"""

from pathlib import Path
from typing import Optional
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from pydantic import BaseModel, ConfigDict
import anndata as ad
from returns.result import Result, Success, Failure
from gene_utils import norm_genes


class RawProcessConfig(BaseModel):
    """Immutable configuration for raw expression matrix conversion to sparse HDF5."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    raw_dir: Path
    out_dir: Path
    min_cell_counts: int = 100


def convert_df_to_sparse_adata(
    df_counts: pd.DataFrame,
    dataset_name: str,
    min_cell_counts: int = 100,
    genes_as_rows: Optional[bool] = None,
) -> Result[ad.AnnData, str]:
    """Converts a pandas DataFrame expression matrix into a sparse CSR AnnData object."""
    try:
        # Determine orientation: single cell matrices typically have genes as rows (>1000 genes)
        if genes_as_rows is True or (genes_as_rows is None and df_counts.shape[0] > 5000 and df_counts.shape[0] > df_counts.shape[1]):
            df_mat = df_counts.T
        else:
            df_mat = df_counts

        cell_barcodes = df_mat.index.astype(str).values
        gene_names = df_mat.columns.astype(str).values

        # Convert to CSR sparse matrix
        matrix_vals = df_mat.values
        X_sparse = csr_matrix(matrix_vals, dtype=np.float32)

        obs = pd.DataFrame({
            "original.barcode": cell_barcodes,
            "dataset": dataset_name,
        }, index=cell_barcodes)

        var = pd.DataFrame(index=gene_names)

        adata = ad.AnnData(X=X_sparse, obs=obs, var=var)
        
        # Filter junk low-depth cells (< min_cell_counts total counts)
        cell_sum = np.asarray(adata.X.sum(axis=1)).ravel()
        adata = adata[cell_sum >= min_cell_counts, :].copy()

        # Harmonize gene names via gene_utils if available/network accessible
        try:
            adata_norm = norm_genes(adata, pre_id_transform="auto")
        except Exception:
            adata_norm = adata

        adata_norm.obs_names_make_unique()
        adata_norm.var_names_make_unique()

        return Success(adata_norm)
    except Exception as e:
        return Failure(f"Failed to convert matrix to sparse AnnData for {dataset_name}: {str(e)}")


def process_tsv_gz_to_sparse_h5(
    file_path: Path,
    dataset_name: str,
    output_path: Path,
    min_cell_counts: int = 100,
) -> Result[Path, str]:
    """Parses a TSV.gz count/TPM matrix or GEO file and exports as a sparse CSR HDF5 .h5ad file."""
    try:
        if not file_path.exists():
            return Failure(f"File not found: {file_path}")

        print(f"Loading raw matrix: {file_path}...")
        
        # Skip comment lines starting with '!' or '#' (GEO series matrix metadata)
        if "matrix" in file_path.name.lower() or "meta" in file_path.name.lower():
            df = pd.read_csv(file_path, sep="\t", comment="!", index_col=0, on_bad_lines="skip")
        elif "tpm" in file_path.name.lower():
            # GSE120575 TPM has cell IDs on row 1 and patient IDs on row 2
            df = pd.read_csv(file_path, sep="\t", skiprows=[1], index_col=0, on_bad_lines="skip")
        else:
            df = pd.read_csv(file_path, sep="\t", comment="#", index_col=0, on_bad_lines="skip")

        match convert_df_to_sparse_adata(df, dataset_name, min_cell_counts):
            case Success(adata):
                output_path.parent.mkdir(parents=True, exist_ok=True)
                ad.settings.allow_write_nullable_strings = True
                adata.write_h5ad(output_path, compression="gzip")
                print(f"Saved sparse AnnData ({adata.shape[0]} cells, {adata.shape[1]} genes) -> {output_path}")
                return Success(output_path)
            case Failure(err):
                return Failure(err)
    except Exception as e:
        return Failure(f"Error processing {file_path}: {str(e)}")


def run_raw_to_sparse_pipeline(config: RawProcessConfig) -> Result[list[Path], str]:
    """Finds all raw expression matrix files and converts them into sparse HDF5 (.h5ad) files."""
    try:
        config.out_dir.mkdir(parents=True, exist_ok=True)
        processed_paths = []

        # Find TSV/Parquet/CSV files in raw_dir
        for sub_dir in config.raw_dir.glob("*"):
            if not sub_dir.is_dir():
                continue
            dataset_name = sub_dir.name
            
            # Look for matrix files
            matrix_files = list(sub_dir.glob("*.txt.gz")) + list(sub_dir.glob("*.parquet")) + list(sub_dir.glob("*.tsv.gz"))
            for mfile in matrix_files:
                out_h5 = config.out_dir / f"{dataset_name}_{mfile.stem.replace('.txt', '')}.h5ad"
                match process_tsv_gz_to_sparse_h5(mfile, dataset_name, out_h5, config.min_cell_counts):
                    case Success(p):
                        processed_paths.append(p)
                    case Failure(err):
                        print(f"Warning on {mfile}: {err}")

        return Success(processed_paths)
    except Exception as e:
        return Failure(f"Raw to sparse pipeline failed: {str(e)}")


def main() -> None:
    """CLI entry point for converting raw GEO matrices into sparse HDF5 (.h5ad) files."""
    import argparse
    parser = argparse.ArgumentParser(description="Convert raw GEO matrices to sparse HDF5 (.h5ad) files.")
    parser.add_argument("--raw-dir", type=str, default="data/raw_geo", help="Input directory containing raw downloaded files")
    parser.add_argument("--out-dir", type=str, default="data/sparse_h5", help="Output directory for sparse .h5ad files")
    parser.add_argument("--min-counts", type=int, default=100, help="Minimum total cell counts threshold")
    args = parser.parse_args()

    config = RawProcessConfig(
        raw_dir=Path(args.raw_dir).resolve(),
        out_dir=Path(args.out_dir).resolve(),
        min_cell_counts=args.min_counts,
    )
    print(f"Starting raw to sparse HDF5 conversion...")
    print(f"Raw dir: {config.raw_dir}")
    print(f"Output dir: {config.out_dir}")

    match run_raw_to_sparse_pipeline(config):
        case Success(paths):
            print(f"\nSuccessfully generated {len(paths)} sparse HDF5 (.h5ad) files:")
            for p in paths:
                print(f" - {p}")
        case Failure(err):
            print(f"Pipeline error: {err}")


if __name__ == "__main__":
    main()
