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
) -> Result[ad.AnnData, str]:
    """Internal dispatch to dataset-specific loader."""
    match dataset_id:
        case "GSE120575":
            # Sade-Feldman
            res = load_sade_feldman(root / "scratch/GSE120575", auto_download=auto_download)
            if not isinstance(res, Success):
                res = load_sade_feldman(root / "data/raw/GSE120575", auto_download=auto_download)
            return res

        case "GSE115978":
            return load_jerby_arnon(root / "data/raw/GSE115978")

        case "GSE125449":
            return load_ma_liver(root / "data/raw/GSE125449")

        case "Maynard_NSCLC":
            return load_maynard(root)

        case _ if dataset_id.endswith("-iAtlas"):
            # cBioPortal iAtlas cohort
            c_dir = root / f"scratch/lair/CBioPortalDataset-{dataset_id}"
            if not c_dir.exists():
                c_dir = root / f"output/CBioPortalDataset-{dataset_id}"
            return load_iatlas_cohort(c_dir, cohort_name=dataset_id, auto_download=auto_download)

        case "EGAD00001006631":
            return load_genentech_egad(root / "manual-download/EGAD00001006631-align")

        case _ if dataset_id in ("Auslander", "Chen-CTLA4", "Chen-PD1", "Freeman", "Gide", "Hugo", "Lauss", "Liu", "Prat", "Ravi", "Riaz", "Rose", "Snyder", "VanAllen"):
            # Direct paper H5AD
            candidate_paths = [
                root / f"dataset_papers/{dataset_id}.h5ad",
                root / f"data/preprocessed/{dataset_id}.h5ad",
                root / f"scratch/lair/ImmuneCheckpointTherapyResponseProcessedGeneNormalizedClinicalDataNormalized/{dataset_id}.h5ad",
                root / f"data/{dataset_id}.h5ad",
            ]
            found = next((p for p in candidate_paths if p.exists()), None)
            if found:
                return load_paper_h5ad(found)
            return Failure(
                f"Paper H5AD '{dataset_id}.h5ad' not found in candidate paths: "
                + ", ".join(str(p) for p in candidate_paths)
            )

        case _:
            return Failure(f"No loader implementation available for dataset '{dataset_id}'")


def load_dataset(
    dataset_id: str,
    base_dir: Path | None = None,
    auto_download: bool = True,
) -> Result[ad.AnnData, str]:
    """Load an individual single-cell or bulk dataset by its registered identifier."""
    logger.info("Loading dataset '%s' (auto_download=%s)...", dataset_id, auto_download)
    start_time = time.time()

    spec_maybe = get_dataset_spec(dataset_id)
    if not isinstance(spec_maybe, Some):
        msg = f"Dataset '{dataset_id}' is not recognized in the registry"
        logger.error(msg)
        return Failure(msg)

    root = base_dir or Path.cwd()
    res = _dispatch_load(dataset_id, root, auto_download)

    elapsed = max(0.01, time.time() - start_time)
    match res:
        case Success(adata):
            logger.info(
                "Successfully loaded dataset '%s': %d obs x %d vars (took %.2fs)",
                dataset_id,
                adata.n_obs,
                adata.n_vars,
                elapsed,
            )
        case Failure(err):
            logger.error("Failed to load dataset '%s': %s (took %.2fs)", dataset_id, err, elapsed)

    return res

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
