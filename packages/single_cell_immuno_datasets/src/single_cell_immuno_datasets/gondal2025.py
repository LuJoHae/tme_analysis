"""Dedicated module for Gondal et al. (2025) integrated ICB scRNA-seq resource."""

from pathlib import Path
import anndata as ad
from returns.result import Result, Success, Failure
from .config import DatasetSpec, DataDirectories, QualityControlSpec, TIER_1_DATASETS
from .downloader import download_single_file
from .preprocessor import preprocess_anndata


def fetch_gondal2025_spec() -> DatasetSpec:
    """Returns the DatasetSpec for Gondal et al. (2025)."""
    return [d for d in TIER_1_DATASETS if d.accession == "Gondal2025"][0]


def download_gondal2025(raw_dir: Path) -> Result[Path, str]:
    """Downloads raw Gondal et al. (2025) integrated h5ad from Zenodo."""
    spec = fetch_gondal2025_spec()
    dest_path = raw_dir / "Gondal2025" / "Gondal2025_scRNAseq_ICB_integrated.h5ad"
    return download_single_file(spec.download_urls["zenodo"], dest_path)


def preprocess_gondal2025(
    raw_h5ad: Path,
    preprocessed_dir: Path,
    qc_spec: QualityControlSpec = QualityControlSpec(),
) -> Result[Path, str]:
    """Preprocesses Gondal et al. (2025) integrated dataset and saves as h5ad."""
    try:
        spec = fetch_gondal2025_spec()
        print(f"Loading Gondal2025 integrated AnnData from {raw_h5ad}...")
        adata = ad.read_h5ad(raw_h5ad)
        
        out_path = preprocessed_dir / "Gondal2025_preprocessed.h5ad"
        return preprocess_anndata(adata, spec, qc_spec, out_path)
    except Exception as e:
        return Failure(f"Failed to process Gondal2025: {str(e)}")
