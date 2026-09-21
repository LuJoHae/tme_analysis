"""Spatial transcriptomics dataset provider: 2D coordinates and neighborhood graphs."""

from __future__ import annotations

from pathlib import Path
import anndata as ad
import numpy as np
import scipy.sparse as sp
from returns.result import Failure, Result, Success
from sklearn.neighbors import NearestNeighbors


def load_spatial_dataset(
    h5ad_or_zarr_path: Path,
    coord_keys: tuple[str, str] = ("x", "y"),
) -> Result[ad.AnnData, str]:
    """Load spatial transcriptomics dataset ensuring 2D coordinates are in adata.obsm['spatial']."""
    if not h5ad_or_zarr_path.exists():
        return Failure(f"Spatial file not found: {h5ad_or_zarr_path}")

    try:
        adata = ad.read_h5ad(h5ad_or_zarr_path)
        if "spatial" not in adata.obsm:
            kx, ky = coord_keys
            if kx in adata.obs.columns and ky in adata.obs.columns:
                coords = np.column_stack([adata.obs[kx].to_numpy(), adata.obs[ky].to_numpy()])
                adata.obsm["spatial"] = coords

        return Success(adata)
    except Exception as exc:
        return Failure(f"Failed to load spatial dataset: {exc}")


def compute_spatial_graph(
    adata: ad.AnnData,
    k_neighbors: int = 6,
) -> Result[ad.AnnData, str]:
    """Calculate k-NN spatial proximity graph storing adjacency in adata.obsp['spatial_connectivities']."""
    if "spatial" not in adata.obsm:
        return Failure("Spatial coordinates not found in adata.obsm['spatial']")

    try:
        new_adata = adata.copy()
        coords = new_adata.obsm["spatial"]

        nn = NearestNeighbors(n_neighbors=k_neighbors + 1, metric="euclidean")
        nn.fit(coords)
        adj = nn.kneighbors_graph(coords, mode="connectivity")
        # Remove self-loops
        adj.setdiag(0)
        adj.eliminate_zeros()

        new_adata.obsp["spatial_connectivities"] = adj
        return Success(new_adata)
    except Exception as exc:
        return Failure(f"Failed to compute spatial neighborhood graph: {exc}")
