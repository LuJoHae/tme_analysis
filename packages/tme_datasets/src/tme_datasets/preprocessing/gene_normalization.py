"""Reusable functional normalization of dataset gene identifiers to Ensembl release."""

from __future__ import annotations

from pathlib import Path
import anndata as ad
from gene_utils import normalize_genes_to_ensembl as _norm_genes_engine
from returns.result import Failure, Result, Success

from ..logging import get_logger
from ..paths import get_data_paths, get_ensembl_dir

logger = get_logger("preprocessing.gene_normalization")


def normalize_dataset_to_ensembl(
    adata: ad.AnnData,
    release: int | None = None,
    species: str | None = None,
    ensembl_dir: Path | None = None,
    drop_unmapped: bool = True,
    aggregation: str = "sum",
) -> Result[ad.AnnData, str]:
    """Convert dataset gene identifiers to canonical Ensembl gene IDs and enrich .var attributes.

    Args:
        adata: Input AnnData container.
        release: Optional Ensembl release version (defaults to config default_ensembl_release, e.g. 111).
        species: Optional species name (defaults to 'human').
        ensembl_dir: Optional installation directory (defaults to config ensembl_dir 'data/ensembl').
        drop_unmapped: Whether to drop unmapped non-gene features (recording them in adata.uns).
        aggregation: Aggregation strategy for duplicate Ensembl IDs ('sum', 'mean', 'max').

    Returns:
        Success(normalized_adata) or Failure(error_message).
    """
    try:
        cfg = get_data_paths()
        target_release = release or cfg.default_ensembl_release
        target_species = species or cfg.default_species
        target_dir = ensembl_dir or cfg.ensembl_dir

        logger.info(
            "Starting Ensembl gene ID normalization for %d genes (Release %d, species=%s, ensembl_dir=%s)...",
            adata.n_vars,
            target_release,
            target_species,
            target_dir,
        )

        norm_adata = _norm_genes_engine(
            adata=adata,
            release=target_release,
            species=target_species,
            ensembl_dir=target_dir,
            drop_unmapped=drop_unmapped,
            aggregation=aggregation,
        )

        logger.info(
            "Ensembl gene ID normalization complete: %d obs x %d vars (all var_names are canonical Ensembl IDs)",
            norm_adata.n_obs,
            norm_adata.n_vars,
        )
        return Success(norm_adata)
    except Exception as exc:
        msg = f"Failed to normalize genes to Ensembl: {exc}"
        logger.error(msg)
        return Failure(msg)
