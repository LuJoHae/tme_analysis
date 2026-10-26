"""Script 3: Quality Control & Preprocessing of Single-Cell Datasets into HDF5 Files.

Inspired by single_cell_datasets (SingleCellDataProcessStep05 & Step07), this script:
1. Performs dataset-specific quality control (min/max genes, min/max counts, mitochondrial content filtering).
2. Removes non-canonical contig & mitochondrial junk genes.
3. Normalizes total counts to 1e4 and log1p transforms expression.
4. Saves one preprocessed .h5ad HDF5 file per dataset.
"""

from pathlib import Path
from typing import NamedTuple, Optional
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict
import anndata as ad
import scanpy as sc
from returns.result import Result, Success, Failure
from gene_utils import norm_genes


class DatasetQCSpec(BaseModel):
    """Quality control threshold specification for a single-cell dataset."""
    model_config = ConfigDict(frozen=True)
    min_genes: int = 300
    max_genes: int = 4500
    min_counts: int = 300
    max_counts: int = 15000
    max_mt_content: float = 20.0


class PreprocessConfig(BaseModel):
    """Immutable configuration for dataset preprocessing into HDF5 files."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    input_dir: Path
    output_dir: Path
    default_qc: DatasetQCSpec = DatasetQCSpec()


def check_and_normalize_transform_space(adata: ad.AnnData) -> tuple[ad.AnnData, dict[str, bool]]:
    """Detects if expression data has already been log-transformed or total-sum normalized.

    Returns updated AnnData and a metadata status dictionary indicating is_log1p and is_normalized.
    """
    adata_out = adata.copy()
    
    # 1. Check uns annotations
    has_uns_log1p = "log1p" in adata_out.uns or "log_transformed" in adata_out.uns
    has_uns_norm = "normalized" in adata_out.uns
    
    # 2. Inspect matrix statistics
    X_mat = adata_out.X
    if hasattr(X_mat, "toarray"):
        sample_vals = X_mat[:min(100, adata_out.n_obs)].toarray()
    else:
        sample_vals = X_mat[:min(100, adata_out.n_obs)]

    max_val = float(sample_vals.max()) if sample_vals.size > 0 else 0.0
    
    # Heuristic: max value <= 35 and presence of non-integer floats indicates log space (e.g. log2(TPM+1) or log1p)
    is_non_integer = not np.all(np.equal(np.mod(sample_vals, 1), 0))
    is_log_transformed = has_uns_log1p or (max_val <= 35.0 and is_non_integer and max_val > 0.0)

    # Check total sums per row
    row_sums = sample_vals.sum(axis=1)
    is_normalized = has_uns_norm or np.allclose(row_sums, 1e4, rtol=1e-2) or np.allclose(row_sums, 1e6, rtol=1e-2)

    status = {
        "is_log1p": is_log_transformed,
        "is_normalized": is_normalized,
        "max_value": max_val,
    }

    # Store detection status in uns metadata
    adata_out.uns["data_transform_state"] = status

    # If in log-space, convert back to linear space for raw count QC / BayesPrism reference preparation
    if is_log_transformed:
        adata_out.layers["log1p_original"] = adata_out.X.copy()
        if hasattr(adata_out.X, "data"):
            # Sparse matrix expm1
            adata_out.X.data = np.expm1(adata_out.X.data)
        else:
            adata_out.X = np.expm1(adata_out.X)

    return adata_out, status


def preprocess_single_dataset(
    h5ad_path: Path,
    qc_spec: DatasetQCSpec,
    output_path: Path,
) -> Result[Path, str]:
    """Applies QC filtering, junk gene removal, normalization, and saves as .h5ad HDF5 file."""
    try:
        if not h5ad_path.exists():
            return Failure(f"Input file not found: {h5ad_path}")

        print(f"Preprocessing dataset: {h5ad_path.name}...")
        adata = ad.read_h5ad(h5ad_path)

        # 1. Detect prior log-transform / normalization
        adata, status = check_and_normalize_transform_space(adata)
        print(f" -> Transform state for {h5ad_path.name}: log_transformed={status['is_log1p']}, normalized={status['is_normalized']}, max_val={status['max_value']:.2f}")

        # 2. Normalize gene names if possible
        try:
            adata = norm_genes(adata, pre_id_transform=None)
        except Exception:
            pass

        # 2. Calculate QC metrics & mitochondrial %
        if "contig" in adata.var.columns:
            adata.var["mt"] = (adata.var["contig"] == "MT").values
        else:
            adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")

        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=["mt"],
            percent_top=None,
            log1p=False,
            inplace=True,
        )

        # 3. Apply QC cell filtering
        sc.pp.filter_cells(adata, min_genes=qc_spec.min_genes)
        sc.pp.filter_cells(adata, min_counts=qc_spec.min_counts)
        sc.pp.filter_cells(adata, max_genes=qc_spec.max_genes)
        sc.pp.filter_cells(adata, max_counts=qc_spec.max_counts)

        if "pct_counts_mt" in adata.obs.columns:
            adata = adata[adata.obs["pct_counts_mt"] < qc_spec.max_mt_content, :].copy()

        # 4. Filter out non-canonical contig & mitochondrial junk genes
        if "contig" in adata.var.columns:
            junk_mask = (adata.var["contig"] == "MT") | (adata.var["contig"].astype(str).str.len() >= 3)
            adata = adata[:, ~junk_mask].copy()

        # 5. Save preprocessed AnnData object to disk as HDF5 (.h5ad)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        ad.settings.allow_write_nullable_strings = True
        adata.write_h5ad(output_path, compression="gzip")

        print(f"Saved preprocessed AnnData ({adata.n_obs} cells, {adata.n_vars} genes) -> {output_path}")
        return Success(output_path)
    except Exception as e:
        return Failure(f"Preprocessing failed for {h5ad_path.name}: {str(e)}")


def run_dataset_preprocessing_pipeline(config: PreprocessConfig) -> Result[list[Path], str]:
    """Iterates through all sparse .h5ad files and applies preprocessing to generate one HDF5 file per dataset."""
    try:
        config.output_dir.mkdir(parents=True, exist_ok=True)
        h5ad_files = list(config.input_dir.glob("*.h5ad"))
        
        if not h5ad_files:
            return Failure(f"No .h5ad files found in input directory: {config.input_dir}")

        processed_paths = []
        for h5_file in h5ad_files:
            out_file = config.output_dir / f"{h5_file.stem}_processed.h5ad"
            match preprocess_single_dataset(h5_file, config.default_qc, out_file):
                case Success(p):
                    processed_paths.append(p)
                case Failure(err):
                    print(f"Warning: {err}")

        return Success(processed_paths)
    except Exception as e:
        return Failure(f"Preprocessing pipeline failed: {str(e)}")


def main() -> None:
    """CLI entry point for dataset quality control and preprocessing."""
    import argparse
    parser = argparse.ArgumentParser(description="Preprocess single-cell datasets into HDF5 (.h5ad) files.")
    parser.add_argument("--input-dir", type=str, default="data/sparse_h5", help="Input directory with sparse .h5ad files")
    parser.add_argument("--out-dir", type=str, default="data/processed_h5", help="Output directory for preprocessed .h5ad files")
    args = parser.parse_args()

    config = PreprocessConfig(
        input_dir=Path(args.input_dir).resolve(),
        output_dir=Path(args.out_dir).resolve(),
    )
    print(f"Starting dataset preprocessing pipeline...")
    print(f"Input dir: {config.input_dir}")
    print(f"Output dir: {config.output_dir}")

    match run_dataset_preprocessing_pipeline(config):
        case Success(paths):
            print(f"\nSuccessfully generated {len(paths)} preprocessed HDF5 (.h5ad) files:")
            for p in paths:
                print(f" - {p}")
        case Failure(err):
            print(f"Pipeline error: {err}")


if __name__ == "__main__":
    main()
