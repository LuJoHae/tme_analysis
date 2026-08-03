"""Scanpy Preprocessing pipeline for GSE97168 single-cell AnnData dataset.

Detects normalization/log-transform status, applies total count normalization & log1p
if necessary, performs HVG identification (without discarding non-HVGs), PCA on HVGs,
kNN graph, Leiden clustering, and UMAP embedding.
"""

import argparse
import sys
from pathlib import Path
from typing import Final, assert_never

import anndata as ad  # type: ignore
import numpy as np  # type: ignore
import scanpy as sc  # type: ignore
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

DEFAULT_TARGET_SUM: Final[float] = 1e4


class PreprocessConfig(BaseModel):
    """Immutable configuration for GSE97168 Scanpy preprocessing."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    input_h5ad: Path
    out_h5ad: Path
    target_sum: float = DEFAULT_TARGET_SUM
    leiden_resolution: float = 0.5


def load_anndata(h5ad_path: Path) -> Result[ad.AnnData, str]:
    """Load AnnData from HDF5 file."""
    try:
        adata = ad.read_h5ad(h5ad_path)
        return Success(adata)
    except Exception as err:
        return Failure(f"Failed to load AnnData from {h5ad_path}: {err}")


def normalize_and_log_transform(adata: ad.AnnData, target_sum: float) -> Result[ad.AnnData, str]:
    """Detect matrix normalization & log-transform status, applying TPM/total count norm & log1p if needed."""
    try:
        X = adata.X
        if sp_is_sparse := hasattr(X, "data"):
            max_val = float(X.data.max()) if X.data.size > 0 else 0.0
            min_val = float(X.data.min()) if X.data.size > 0 else 0.0
            is_integer = bool(np.all(X.data == np.floor(X.data)))
        else:
            max_val = float(np.max(X))
            min_val = float(np.min(X))
            is_integer = bool(np.all(X == np.floor(X)))

        print(f"Matrix QC stats: min={min_val:.4f}, max={max_val:.4f}, integer_values={is_integer}")

        if max_val < 25.0 and not is_integer and min_val >= 0.0:
            print("Data detected as already log-transformed (max < 25.0, non-integer). Skipping log1p.")
        else:
            if is_integer or max_val >= 100.0:
                print(f"Raw count data detected. Normalizing total counts to {target_sum:.0f} per cell...")
                sc.pp.normalize_total(adata, target_sum=target_sum)
            else:
                print(f"Un-logged TPM/CPM detected. Rescaling total counts to {target_sum:.0f} per cell...")
                sc.pp.normalize_total(adata, target_sum=target_sum)

            print("Applying log1p transformation...")
            sc.pp.log1p(adata)

        return Success(adata)
    except Exception as err:
        return Failure(f"Failed during normalization/log-transform: {err}")


def compute_hvg_pca_umap(adata: ad.AnnData, leiden_res: float) -> Result[ad.AnnData, str]:
    """Compute HVGs (without discarding non-HVGs), PCA, KNN graph, Leiden clustering, and UMAP."""
    try:
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)

        print("Identifying highly variable genes (HVGs) without discarding non-HVGs...")
        sc.pp.highly_variable_genes(
            adata,
            min_mean=0.0125,
            max_mean=3.0,
            min_disp=0.5,
            inplace=True,
        )

        n_hvg = int(adata.var["highly_variable"].sum())
        print(f"Identified {n_hvg} highly variable genes out of {adata.n_vars} total genes.")

        print("Computing PCA based on HVGs...")
        sc.tl.pca(adata, svd_solver="arpack", mask_var="highly_variable")

        print("Computing nearest neighbor graph...")
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)

        print(f"Computing Leiden clustering at resolution {leiden_res}...")
        sc.tl.leiden(adata, resolution=leiden_res, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)

        print("Computing UMAP embedding...")
        sc.tl.umap(adata)

        return Success(adata)
    except Exception as err:
        return Failure(f"Failed computing HVG/PCA/KNN/Leiden/UMAP: {err}")


def save_anndata(adata: ad.AnnData, out_h5ad: Path) -> Result[Path, str]:
    """Save processed AnnData object to HDF5 file."""
    try:
        out_h5ad.parent.mkdir(parents=True, exist_ok=True)
        print(f"Saving processed AnnData ({adata.n_obs} cells x {adata.n_vars} genes) to {out_h5ad}...")

        # Sanitize metadata columns to prevent string encoding issues in HDF5
        for col in adata.obs.columns:
            adata.obs[col] = adata.obs[col].astype(str).replace("nan", "")

        adata.write_h5ad(str(out_h5ad), compression="gzip")
        return Success(out_h5ad)
    except Exception as err:
        return Failure(f"Failed to save AnnData to {out_h5ad}: {err}")


def run_pipeline(config: PreprocessConfig) -> Result[Path, str]:
    """Monadic processing pipeline for GSE97168 scanpy preprocessing."""
    match load_anndata(config.input_h5ad):
        case Failure(err):
            return Failure(err)
        case Success(adata):
            match normalize_and_log_transform(adata, config.target_sum):
                case Failure(err):
                    return Failure(err)
                case Success(norm_adata):
                    match compute_hvg_pca_umap(norm_adata, config.leiden_resolution):
                        case Failure(err):
                            return Failure(err)
                        case Success(final_adata):
                            return save_anndata(final_adata, config.out_h5ad)
                        case _ as unreachable:
                            assert_never(unreachable)
                case _ as unreachable:
                    assert_never(unreachable)
        case _ as unreachable:
            assert_never(unreachable)


def main() -> None:
    """CLI entry point for GSE97168 Scanpy preprocessing."""
    parser = argparse.ArgumentParser(description="Preprocess GSE97168 dataset (TPM/log1p, HVG, PCA, KNN, Leiden, UMAP).")
    _ = parser.add_argument("--input-h5ad", required=True, help="Input raw AnnData HDF5 file.")
    _ = parser.add_argument("--out-h5ad", required=True, help="Output preprocessed AnnData HDF5 file.")
    _ = parser.add_argument("--target-sum", type=float, default=DEFAULT_TARGET_SUM, help="Target sum per cell for normalization.")
    _ = parser.add_argument("--resolution", type=float, default=0.5, help="Leiden clustering resolution.")

    args = parser.parse_args()
    config = PreprocessConfig(
        input_h5ad=Path(args.input_h5ad),
        out_h5ad=Path(args.out_h5ad),
        target_sum=args.target_sum,
        leiden_resolution=args.resolution,
    )

    match run_pipeline(config):
        case Success(out_path):
            print(f"Successfully finished preprocessing GSE97168 dataset: {out_path}")
            sys.exit(0)
        case Failure(err):
            print(f"Preprocessing error: {err}", file=sys.stderr)
            sys.exit(1)
        case _ as unreachable:
            assert_never(unreachable)


if __name__ == "__main__":
    main()
