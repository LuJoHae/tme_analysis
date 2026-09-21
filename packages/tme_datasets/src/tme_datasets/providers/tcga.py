"""TCGA Pan-Cancer dataset provider and background reference cohorts."""

from __future__ import annotations

from pathlib import Path
import anndata as ad
from returns.result import Failure, Result, Success


def load_tcga_project(h5ad_path: Path) -> Result[ad.AnnData, str]:
    """Load a single TCGA cancer project AnnData file."""
    if not h5ad_path.exists():
        return Failure(f"TCGA file not found: {h5ad_path}")

    try:
        adata = ad.read_h5ad(h5ad_path)
        return Success(adata)
    except Exception as exc:
        return Failure(f"Failed to load TCGA project from {h5ad_path}: {exc}")
