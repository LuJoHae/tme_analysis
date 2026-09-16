"""Script 2: Convert Raw Single-Cell Matrices to Sparse HDF5 (.h5ad) Files with Per-Dataset Gene Normalization.

Inspired by single_cell_datasets (SingleCellDataProcessStep01 & Step04), this script:
1. Loads raw expression matrices (H5AD, Parquet, TSV.gz, CSV).
2. Converts expression data to scipy.sparse.csr_matrix (float32, gzip compressed) per dataset.
3. Normalizes all gene IDs for every dataset using gene_utils.norm_genes(pre_id_transform="auto").
4. Saves sparse .h5ad files to output directory.
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
    """Converts a pandas DataFrame expression matrix into a sparse CSR AnnData object with normalized gene IDs."""
    try:
        # Determine orientation: single cell matrices typically have genes as rows (>5000 genes)
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

        # Normalize all gene IDs per dataset using gene_utils
        try:
            adata_norm = norm_genes(adata, pre_id_transform="auto")
        except Exception:
            adata_norm = adata

        adata_norm.obs_names_make_unique()
        adata_norm.var_names_make_unique()

        return Success(adata_norm)
    except Exception as e:
        return Failure(f"Failed to convert matrix to sparse AnnData for {dataset_name}: {str(e)}")


def process_raw_file_to_sparse_h5(
    file_path: Path,
    dataset_name: str,
    output_path: Path,
    min_cell_counts: int = 100,
) -> Result[Path, str]:
    """Reads raw file (H5AD, Parquet, or TSV.gz) and exports a sparse CSR HDF5 .h5ad file."""
    try:
        if not file_path.exists():
            return Failure(f"File not found: {file_path}")

        print(f"Loading raw file: {file_path.name} for dataset '{dataset_name}'...")
        
        if file_path.suffix == ".h5ad":
            adata = ad.read_h5ad(file_path)
            adata.X = csr_matrix(adata.X, dtype=np.float32)
            adata.obs["dataset"] = dataset_name
            try:
                adata = norm_genes(adata, pre_id_transform="auto")
            except Exception:
                pass
        elif file_path.suffix == ".parquet":
            df = pd.read_parquet(file_path)
            match convert_df_to_sparse_adata(df, dataset_name, min_cell_counts):
                case Success(a):
                    adata = a
                case Failure(err):
                    return Failure(err)
        else:
            # Parse TSV / CSV / TXT.gz
            if "matrix" in file_path.name.lower() or "meta" in file_path.name.lower():
                df = pd.read_csv(file_path, sep="\t", comment="!", index_col=0, on_bad_lines="skip")
            elif "tpm" in file_path.name.lower():
                df = pd.read_csv(file_path, sep="\t", skiprows=[1], index_col=0, on_bad_lines="skip")
            else:
                df = pd.read_csv(file_path, sep="\t", comment="#", index_col=0, on_bad_lines="skip")

            match convert_df_to_sparse_adata(df, dataset_name, min_cell_counts):
                case Success(a):
                    adata = a
                case Failure(err):
                    return Failure(err)

        # Guarantee CSR sparse matrix format for HDF5 compatibility
        import scipy.sparse
        import numpy as np
        if not scipy.sparse.isspmatrix_csr(adata.X):
            if scipy.sparse.issparse(adata.X):
                adata.X = adata.X.tocsr()
            else:
                adata.X = scipy.sparse.csr_matrix(np.asarray(adata.X))

        output_path.parent.mkdir(parents=True, exist_ok=True)
        ad.settings.allow_write_nullable_strings = True
        adata.write_h5ad(output_path, compression="gzip")
        print(f"Saved sparse AnnData ({adata.n_obs} cells, {adata.n_vars} genes) -> {output_path}")
        return Success(output_path)
    except Exception as e:
        return Failure(f"Error processing {file_path}: {str(e)}")


def run_raw_to_sparse_pipeline(config: RawProcessConfig, n_jobs: int = -1) -> Result[list[Path], str]:
    """Finds all raw matrix files and converts each dataset into its sparse HDF5 (.h5ad) file in parallel."""
    try:
        from concurrent.futures import ProcessPoolExecutor
        import os

        config.out_dir.mkdir(parents=True, exist_ok=True)
        tasks = []

        # Find all raw matrix tasks
        for entry in config.raw_dir.glob("*"):
            if entry.is_dir():
                dataset_name = entry.name
                all_files = list(entry.glob("*.h5ad")) + list(entry.glob("*.parquet")) + list(entry.glob("*.txt.gz")) + list(entry.glob("*.tsv.gz"))
                mfile = find_expression_matrix_file(all_files)
                if mfile:
                    out_h5 = config.out_dir / f"{dataset_name}_sparse.h5ad"
                    tasks.append((mfile, dataset_name, out_h5, config.min_cell_counts))
            elif entry.suffix in [".h5ad", ".parquet", ".gz"]:
                if "meta" not in entry.name.lower() and "tcr" not in entry.name.lower():
                    dataset_name = entry.stem.split(".")[0]
                    out_h5 = config.out_dir / f"{dataset_name}_sparse.h5ad"
                    tasks.append((entry, dataset_name, out_h5, config.min_cell_counts))

        if not tasks:
            return Success([])

        workers = os.cpu_count() if n_jobs <= 0 else n_jobs
        print(f"Executing raw to sparse conversion for {len(tasks)} tasks in parallel using {workers} workers...")

        processed_paths = []
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(process_raw_file_to_sparse_h5, mfile, dname, out_h5, min_counts)
                for (mfile, dname, out_h5, min_counts) in tasks
            ]
            for future in futures:
                match future.result():
                    case Success(p):
                        processed_paths.append(p)
                    case Failure(err):
                        print(f"Warning: {err}")

        return Success(processed_paths)
    except Exception as e:
        return Failure(f"Raw to sparse pipeline failed: {str(e)}")


def find_expression_matrix_file(files: list[Path]) -> Optional[Path]:
    """Selects the actual expression matrix file from a list of candidate raw dataset files."""
    if not files:
        return None
    # 1. Filter out metadata, TCR, and non-expression files
    valid_candidates = [
        f for f in files 
        if "meta" not in f.name.lower() and "tcr" not in f.name.lower() and "series_matrix" not in f.name.lower()
    ]
    if not valid_candidates:
        valid_candidates = files

    # 2. Prioritize files containing expression keywords
    for kw in ["tpm", "count", "umitab", "matrix", "expr"]:
        for f in valid_candidates:
            if kw in f.name.lower():
                return f

    return valid_candidates[0]


def main() -> None:
    """CLI entry point for converting raw GEO matrices into sparse HDF5 (.h5ad) files."""
    import argparse
    parser = argparse.ArgumentParser(description="Convert raw GEO matrices to sparse HDF5 (.h5ad) files.")
    parser.add_argument("--raw-dir", type=str, default="data/raw_geo", help="Input directory containing raw downloaded files")
    parser.add_argument("--out-dir", type=str, default="data/sparse_h5", help="Output directory for sparse .h5ad files")
    parser.add_argument("--dataset", type=str, default=None, help="Optional single dataset name to process")
    parser.add_argument("--out-h5", type=str, default=None, help="Optional specific output h5ad path for single dataset")
    parser.add_argument("--min-counts", type=int, default=100, help="Minimum total cell counts threshold")
    args = parser.parse_args()

    raw_dir = Path(args.raw_dir).resolve()
    out_dir = Path(args.out_dir).resolve()

    if args.dataset and args.out_h5:
        # Single dataset execution mode for Snakemake wildcards
        dname = args.dataset
        out_h5 = Path(args.out_h5).resolve()
        
        # Find matching file in target_dir or raw_dir
        target_dir = raw_dir / dname
        raw_files = []
        if target_dir.exists():
            raw_files = (
                list(target_dir.glob("*.h5ad"))
                + list(target_dir.glob("*.parquet"))
                + list(target_dir.glob("*.txt.gz"))
                + list(target_dir.glob("*.tsv.gz"))
                + list(target_dir.glob("*.h5"))
            )
        
        if not raw_files:
            candidate_dirs = [
                raw_dir,
                Path("/storage/halu/dataset_datalair/SingleCellDataProcessStep01"),
                Path("/storage/halu/dataset_datalair/SingleCellDataProcessStep03"),
                Path("/storage/halu/dataset_datalair"),
                Path.home() / ".cache" / "datalair",
            ]
            for cdir in candidate_dirs:
                if cdir.exists():
                    matches = [
                        f for f in cdir.rglob("*")
                        if f.suffix in [".h5ad", ".parquet", ".gz", ".h5"]
                        and ("meta" not in f.name.lower() and "tcr" not in f.name.lower() and "series_matrix" not in f.name.lower())
                        and (dname.lower() in f.name.lower() or dname[:6].lower() in f.name.lower())
                    ]
                    if matches:
                        raw_files.extend(matches)

        expr_file = find_expression_matrix_file(raw_files)
        
        if not expr_file:
            print(f"Error: No raw expression file found for dataset '{dname}' in candidates")
            return

        print(f"Selected expression file for '{dname}': {expr_file.name}")
        match process_raw_file_to_sparse_h5(expr_file, dname, out_h5, args.min_counts):
            case Success(p):
                print(f"Successfully processed single dataset '{dname}' -> {p}")
            case Failure(err):
                print(f"Error processing single dataset '{dname}': {err}")
    else:
        config = RawProcessConfig(
            raw_dir=raw_dir,
            out_dir=out_dir,
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
