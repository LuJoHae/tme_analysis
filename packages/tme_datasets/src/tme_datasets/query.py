"""High-level declarative dataset querying and harmonization API."""

from __future__ import annotations

from pathlib import Path
import time
from typing import Sequence
import anndata as ad
from returns.maybe import Nothing, Some
from returns.result import Failure, Result, Success

from .genes.reconcile import reconcile_genes
from .genesets.collections import (
    get_bagaev_core_collection,
    get_tme_major_lineage_collection,
    get_tme_subtype_collection,
)
from .genesets.models import GeneSetCollection
from .harmonization.align import align_and_concatenate
from .logging import get_logger
from .models import GeneReconcileConfig, HarmonizeConfig
from .paths import (
    find_dataset_h5ad,
    get_preprocessed_h5ad_path,
    get_raw_dataset_dir,
    get_scratch_dataset_dir,
)
from .providers.bulk_iatlas import load_iatlas_cohort
from .providers.bulk_papers import load_genentech_egad, load_paper_h5ad
from .providers.single_cell import (
    load_jerby_arnon,
    load_ma_liver,
    load_maynard,
    load_sade_feldman,
)
from .registry import get_dataset_spec, list_registered_datasets

logger = get_logger("query")


def _dispatch_load(
    dataset_id: str,
    root: Path,
    auto_download: bool,
    force_download: bool = False,
) -> Result[ad.AnnData, str]:
    """Internal dispatch to dataset-specific loader using config-resolved paths."""
    raw_dir = get_raw_dataset_dir(dataset_id, repo_root=root)
    scratch_dir = get_scratch_dataset_dir(dataset_id, repo_root=root)

    match dataset_id:
        case "GSE120575":
            # Sade-Feldman
            res = load_sade_feldman(scratch_dir, auto_download=auto_download, force_download=force_download)
            if not isinstance(res, Success):
                res = load_sade_feldman(raw_dir, auto_download=auto_download, force_download=force_download)
            return res

        case "GSE115978":
            return load_jerby_arnon(raw_dir, auto_download=auto_download, force_download=force_download)

        case "GSE125449":
            return load_ma_liver(raw_dir)

        case "Maynard_NSCLC":
            return load_maynard(root)

        case _ if dataset_id.endswith("-iAtlas"):
            # cBioPortal iAtlas cohort
            return load_iatlas_cohort(raw_dir, cohort_name=dataset_id, auto_download=auto_download, force_download=force_download)

        case "EGAD00001006631":
            return load_genentech_egad(raw_dir)

        case _ if dataset_id in (
            "Auslander", "Chen-CTLA4", "Chen-PD1", "Freeman", "Gide", "Hugo",
            "Lauss", "Liu", "Prat", "Ravi", "Riaz", "Rose", "Snyder", "VanAllen"
        ):
            # Direct paper H5AD
            found = find_dataset_h5ad(dataset_id, repo_root=root)
            match found:
                case Some(h5ad_path):
                    return load_paper_h5ad(h5ad_path)
                case _:
                    return Failure(f"Paper H5AD for '{dataset_id}' not found in candidate paths")

        case _:
            return Failure(f"No loader implementation available for dataset '{dataset_id}'")


def load_dataset(
    dataset_id: str,
    base_dir: Path | None = None,
    auto_download: bool = True,
    force_recompute: bool = False,
    force_download: bool = False,
    cache_h5ad: bool = True,
    normalize_ensembl: bool = True,
    ensembl_release: int | None = None,
    drop_unmapped: bool = True,
) -> Result[ad.AnnData, str]:
    """Load an individual single-cell or bulk dataset by its registered identifier.

    Prioritizes loading directly from cached H5AD in <0.5s. If no H5AD exists (or
    force_recompute=True), auto-downloads raw files, parses them, normalizes gene IDs
    to canonical Ensembl IDs (Release 111 by default), enriches .var attributes,
    serializes to an H5AD cache file, and returns the AnnData object.

    Args:
        dataset_id: Registered dataset identifier (e.g. 'GSE120575', 'Hugo-iAtlas').
        base_dir: Optional root directory override.
        auto_download: Automatically download missing raw vendor files.
        force_recompute: Overwrite cached H5AD and re-parse from raw files.
        force_download: Force re-downloading raw vendor files from network.
        cache_h5ad: Write parsed AnnData to canonical H5AD cache upon completion.
        normalize_ensembl: Normalize gene IDs to canonical Ensembl identifiers and enrich .var.
        ensembl_release: Specific Ensembl release version (defaults to config release, e.g. 111).
        drop_unmapped: Drop unmapped non-gene features (saving them to adata.uns['unmapped_genes']).

    Returns:
        Success(adata) or Failure(error_message).
    """
    logger.info(
        "Loading dataset '%s' (auto_download=%s, force_recompute=%s, force_download=%s, normalize_ensembl=%s)...",
        dataset_id,
        auto_download,
        force_recompute,
        force_download,
        normalize_ensembl,
    )
    start_time = time.time()

    spec_maybe = get_dataset_spec(dataset_id)
    if not isinstance(spec_maybe, Some):
        msg = f"Dataset '{dataset_id}' is not recognized in the registry"
        logger.error(msg)
        return Failure(msg)

    root = base_dir or Path.cwd()

    # Step 1: Fast path - load from preprocessed H5AD if available
    if not force_recompute:
        h5ad_found = find_dataset_h5ad(dataset_id, repo_root=root)
        match h5ad_found:
            case Some(h5ad_path):
                size_mb = h5ad_path.stat().st_size / (1024 * 1024)
                logger.info(
                    "Found cached H5AD for '%s' at %s (%.1f MB). Loading directly in <0.5s...",
                    dataset_id,
                    h5ad_path.name,
                    size_mb,
                )
                try:
                    adata = ad.read_h5ad(h5ad_path)
                    elapsed = max(0.01, time.time() - start_time)
                    logger.info(
                        "Successfully loaded dataset '%s': %d obs x %d vars from H5AD (took %.2fs)",
                        dataset_id,
                        adata.n_obs,
                        adata.n_vars,
                        elapsed,
                    )
                    return Success(adata)
                except Exception as exc:
                    logger.warning(
                        "Failed to read cached H5AD at %s: %s. Falling back to raw ingestion.",
                        h5ad_path,
                        exc,
                    )

    # Step 2: Ingestion path - dispatch to provider loader
    res = _dispatch_load(dataset_id, root, auto_download=auto_download, force_download=force_download)

    # Step 3: Ensembl Normalization and H5AD Serialization
    match res:
        case Success(adata):
            if normalize_ensembl:
                from .preprocessing.gene_normalization import normalize_dataset_to_ensembl

                norm_res = normalize_dataset_to_ensembl(
                    adata,
                    release=ensembl_release,
                    drop_unmapped=drop_unmapped,
                )
                match norm_res:
                    case Success(norm_adata):
                        adata = norm_adata
                    case Failure(norm_err):
                        logger.warning(
                            "Ensembl gene normalization encountered an error: %s. Proceeding with raw identifiers.",
                            norm_err,
                        )

            if cache_h5ad:
                try:
                    target_h5ad = get_preprocessed_h5ad_path(dataset_id, repo_root=root)
                    target_h5ad.parent.mkdir(parents=True, exist_ok=True)

                    # Sanitize index names for robust h5py serialization
                    if adata.obs_names.name is None or not isinstance(adata.obs_names.name, str):
                        adata.obs_names.name = "sample_id" if "cell_id" not in adata.obs.columns else "cell_id"
                    if adata.var_names.name is None or not isinstance(adata.var_names.name, str):
                        adata.var_names.name = "gene_id"

                    logger.info("Writing parsed AnnData to H5AD cache at %s...", target_h5ad)
                    adata.write_h5ad(target_h5ad)
                    size_mb = target_h5ad.stat().st_size / (1024 * 1024)
                    logger.info(
                        "Successfully serialized dataset '%s' to H5AD: %s (%.1f MB). Future calls will load instantly.",
                        dataset_id,
                        target_h5ad.name,
                        size_mb,
                    )
                except Exception as cache_exc:
                    logger.warning("Failed to serialize H5AD cache for '%s': %s", dataset_id, cache_exc)

            elapsed = max(0.01, time.time() - start_time)
            logger.info(
                "Successfully loaded dataset '%s': %d obs x %d vars (took %.2fs)",
                dataset_id,
                adata.n_obs,
                adata.n_vars,
                elapsed,
            )
            return Success(adata)
        case Failure(err):
            elapsed = max(0.01, time.time() - start_time)
            logger.error("Failed to load dataset '%s': %s (took %.2fs)", dataset_id, err, elapsed)
            return Failure(err)

def query_datasets(
    dataset_ids: Sequence[str],
    config: HarmonizeConfig | None = None,
    base_dir: Path | None = None,
) -> Result[ad.AnnData, str]:
    """Query and harmonize multiple single-cell or bulk datasets into a unified AnnData object."""
    if not dataset_ids:
        msg = "At least one dataset ID must be provided"
        logger.error(msg)
        return Failure(msg)

    cfg = config or HarmonizeConfig()
    logger.info(
        "Querying and harmonizing %d datasets: %s (mode=%s, reconcile_genes=%s)",
        len(dataset_ids),
        list(dataset_ids),
        cfg.mode.value,
        cfg.reconcile_genes,
    )
    loaded_adatas = []

    for ds_id in dataset_ids:
        match load_dataset(ds_id, base_dir=base_dir):
            case Failure(err):
                return Failure(f"Failed to load dataset '{ds_id}': {err}")
            case Success(adata):
                if cfg.reconcile_genes:
                    logger.debug("Reconciling gene identifiers for '%s'...", ds_id)
                    reconcile_res = reconcile_genes(
                        adata,
                        GeneReconcileConfig(target_type=cfg.gene_target_type),
                    )
                    match reconcile_res:
                        case Success(rec_adata):
                            loaded_adatas.append(rec_adata)
                        case Failure(err):
                            return Failure(f"Gene reconciliation failed for '{ds_id}': {err}")
                else:
                    loaded_adatas.append(adata)

    logger.info("Aligning and concatenating %d loaded datasets...", len(loaded_adatas))
    res = align_and_concatenate(loaded_adatas, dataset_ids, cfg)
    match res:
        case Success(unified):
            logger.info(
                "Successfully harmonized %d datasets: %d obs x %d vars",
                len(loaded_adatas),
                unified.n_obs,
                unified.n_vars,
            )
        case Failure(err):
            logger.error("Failed to align and concatenate datasets: %s", err)
    return res


def load_geneset_collection(collection_id: str) -> Result[GeneSetCollection, str]:
    """Retrieve pre-registered gene set collection by ID."""
    match collection_id:
        case "tme_major_lineages":
            return Success(get_tme_major_lineage_collection())
        case "tme_subtypes":
            return Success(get_tme_subtype_collection())
        case "bagaev_core":
            return Success(get_bagaev_core_collection())
        case _:
            return Failure(f"Gene set collection '{collection_id}' not found")
