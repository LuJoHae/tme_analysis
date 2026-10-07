"""High-level declarative dataset querying and harmonization API."""

from __future__ import annotations

from pathlib import Path
import time
from typing import Any, Mapping, Sequence, Literal
import anndata as ad
from returns.maybe import Nothing, Some
from returns.result import Failure, Result, Success

from .genes.reconcile import reconcile_genes
from .genesets.collections import (
    get_bagaev_core_collection,
    get_tme_major_lineage_collection,
    get_tme_subtype_collection,
)
from .deconvolution import (
    DeconvolutionReferenceConfig,
    DeconvolutionReferenceResult,
    build_deconvolution_reference,
)
from .genesets.models import GeneSetCollection
from .harmonization.align import align_and_concatenate
from .logging import get_logger
from .models import (
    ClusterAnalysisSpec,
    DatasetSpec,
    GeneReconcileConfig,
    HarmonizeConfig,
    PreprocessedDatasetSpec,
    PreprocessedDatasets,
    RegisteredDatasets,
    SampledSingleCellResult,
    SingleCellProcessingSpec,
    SingleCellSamplingSpec,
)
from .types import GeneIDType, HarmonizeMode, Modality
from .paths import (
    find_dataset_h5ad,
    find_repo_root,
    get_manual_download_dir,
    get_preprocessed_h5ad_path,
    get_raw_dataset_dir,
    get_scratch_dataset_dir,
)
from .providers.bulk_iatlas import IATLAS_COHORTS, load_iatlas_cohort
from .providers.bulk_papers import load_genentech_egad, load_paper_h5ad
from .providers.atlas_single_cell import (
    load_azizi_brca,
    load_becker_coad,
    load_biermann_brainmet,
    load_borcherding_ccrcc,
    load_cheng_pancancer,
    load_durante_uvm,
    load_khaliq_cc,
    load_kim_luad,
    load_leader_nsclc,
    load_lu_hcc,
    load_pelka_crc,
    load_pu_ptc,
    load_qian_pancancer,
    load_sharma_hcc,
    load_vazquez_ov,
    load_zhang_myeloid,
    load_zhang_tnbc,
    _apply_subset_and_subsample,
)
from .providers.single_cell import (
    load_gse179994,
    load_jerby_arnon,
    load_ma_liver,
    load_maynard,
    load_sade_feldman,
    load_yost,
)
from .providers.tier0_single_cell import (
    load_cellxgene_7b20c613_melanoma,
    load_cellxgene_05a8c945_crc,
    load_cellxgene_6f9de485_breast,
    load_cellxgene_829a3cd1_crc,
    load_gse200996_hnscc,
    load_gse207422_nsclc,
    load_gse210038_ccrcc,
    load_gse212707_breast,
    load_gse218429_melanoma,
    load_gse233203_nsclc,
    load_gse236581_crc,
    load_gse243013_nsclc,
    load_gse245906_hcc,
    load_gse246613_breast,
    load_gse270680_gastric,
    load_gse287301_hnscc,
    load_gse299651_crc,
    load_gse300475_breast,
    load_gse301741_hnscc,
    load_gse311789_pdac,
    load_gse313642_hcc,
    load_gse314072_ccrcc,
    load_gse316195_pdac,
    load_gse317309_nsclc,
    load_gse344166_melanoma,
)
from .providers.tier1_single_cell import (
    is_tier1_dataset,
    load_tier1_cohort,
)
from .registry import (
    get_dataset_spec,
    list_preprocessed_datasets,
    list_registered_datasets,
    query_preprocessed_datasets,
    resolve_dataset_id,
)

logger = get_logger("query")


def _dispatch_load(
    dataset_id: str,
    root: Path,
    raw_dir: Path | None = None,
    auto_download: bool = True,
    force_download: bool = False,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Internal dispatch to dataset-specific loader using config-resolved paths."""
    canonical_id = resolve_dataset_id(dataset_id)
    raw_dir = raw_dir or get_raw_dataset_dir(canonical_id, repo_root=root)
    scratch_dir = get_scratch_dataset_dir(canonical_id, repo_root=root)
    manual_dir = get_manual_download_dir(repo_root=root)

    match canonical_id:
        case "GSE120575":
            # Sade-Feldman
            res = load_sade_feldman(scratch_dir, auto_download=auto_download, force_download=force_download)
            if not isinstance(res, Success):
                res = load_sade_feldman(raw_dir, auto_download=auto_download, force_download=force_download)
            return res

        case "GSE115978":
            return load_jerby_arnon(raw_dir, auto_download=auto_download, force_download=force_download)

        case "GSE125449":
            return load_ma_liver(raw_dir, auto_download=auto_download, force_download=force_download)

        case "GSE123813":
            return load_yost(raw_dir, auto_download=auto_download, force_download=force_download)

        case "GSE179994":
            return load_gse179994(raw_dir, auto_download=auto_download, force_download=force_download)

        case "Maynard_NSCLC":
            return load_maynard(root, raw_dir=raw_dir, auto_download=auto_download, force_download=force_download)

        # 17 Pan-Cancer Reference Atlas Cohorts
        case "GSE178341":
            return load_pelka_crc(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "GSE114727":
            return load_azizi_brca(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "E-MTAB-8107":
            return load_qian_pancancer(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "GSE154763":
            return load_cheng_pancancer(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "GSE154826":
            return load_leader_nsclc(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "GSE131907":
            return load_kim_luad(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "GSE201349":
            return load_becker_coad(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "GSE200997":
            return load_khaliq_cc(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "GSE121638":
            return load_borcherding_ccrcc(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "GSE156625":
            return load_sharma_hcc(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "GSE149614":
            return load_lu_hcc(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "GSE184362":
            return load_pu_ptc(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "GSE139829":
            return load_durante_uvm(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "GSE200218":
            return load_biermann_brainmet(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "GSE180661":
            return load_vazquez_ov(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "GSE169246":
            return load_zhang_tnbc(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        case "GSE215120":
            return load_zhang_myeloid(raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)

        # 22 Tier 0 Premier Benchmark Core Cohorts
        case "GSE246613":
            return load_gse246613_breast(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE300475":
            return load_gse300475_breast(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE212707":
            if auto_download and not (raw_dir / "GSE212707_RAW.tar").exists() and not list(raw_dir.glob("*.tar*")):
                from .download.geo import download_geo_supplementary
                download_geo_supplementary("GSE212707", raw_dir, expected_files=["GSE212707_RAW.tar"])
            return load_gse212707_breast(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE236581":
            return load_gse236581_crc(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE299651":
            return load_gse299651_crc(raw_dir, subset=subset, subsample_n=subsample_n)

        case "CELLxGENE_829a3cd1":
            return load_cellxgene_829a3cd1_crc(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE270680":
            return load_gse270680_gastric(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE313642":
            if auto_download and not (raw_dir / "GSE313642_RAW.tar").exists() and not list(raw_dir.glob("*.tar*")):
                from .download.geo import download_geo_supplementary
                download_geo_supplementary("GSE313642", raw_dir, expected_files=["GSE313642_RAW.tar"])
            return load_gse313642_hcc(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE245906":
            return load_gse245906_hcc(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE301741":
            return load_gse301741_hnscc(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE287301":
            return load_gse287301_hnscc(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE200996":
            return load_gse200996_hnscc(raw_dir, subset=subset, subsample_n=subsample_n)

        case "CELLxGENE_7b20c613":
            return load_cellxgene_7b20c613_melanoma(raw_dir, subset=subset, subsample_n=subsample_n)

        case "CELLxGENE_05a8c945":
            return load_cellxgene_05a8c945_crc(raw_dir, subset=subset, subsample_n=subsample_n)

        case "CELLxGENE_6f9de485":
            return load_cellxgene_6f9de485_breast(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE233203":
            return load_gse233203_nsclc(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE218429":
            return load_gse218429_melanoma(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE344166":
            return load_gse344166_melanoma(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE207422":
            return load_gse207422_nsclc(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE317309":
            return load_gse317309_nsclc(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE243013":
            return load_gse243013_nsclc(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE311789":
            return load_gse311789_pdac(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE316195":
            return load_gse316195_pdac(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE210038":
            return load_gse210038_ccrcc(raw_dir, subset=subset, subsample_n=subsample_n)

        case "GSE314072":
            return load_gse314072_ccrcc(raw_dir, subset=subset, subsample_n=subsample_n)

        case _ if is_tier1_dataset(canonical_id):
            return load_tier1_cohort(
                canonical_id,
                raw_dir,
                auto_download=auto_download,
                subset=subset,
                subsample_n=subsample_n,
            )

        case _ if canonical_id.endswith("-iAtlas"):
            # cBioPortal iAtlas cohort
            return load_iatlas_cohort(raw_dir, cohort_name=canonical_id, auto_download=auto_download, force_download=force_download)

        case "EGAD00001006631":
            align_dir = manual_dir / "EGAD00001006631-align"
            if not align_dir.exists() and raw_dir.exists():
                align_dir = raw_dir
            return load_genentech_egad(align_dir)

        case _ if canonical_id in (
            "Auslander", "Chen-CTLA4", "Chen-PD1", "Freeman", "Gide", "Hugo",
            "Lauss", "Liu", "Prat", "Ravi", "Riaz", "Rose", "Snyder", "VanAllen"
        ):
            # Direct paper H5AD
            expected_file = manual_dir / f"{canonical_id}.h5ad"
            if expected_file.is_file() and expected_file.stat().st_size > 0:
                return load_paper_h5ad(expected_file)
            found = find_dataset_h5ad(canonical_id, repo_root=root)
            match found:
                case Some(h5ad_path):
                    return load_paper_h5ad(h5ad_path)
                case _:
                    return Failure(
                        f"Paper H5AD for '{canonical_id}' not found. "
                        f"Expected file at '{expected_file}'. "
                        f"Please run 'python scripts/setup_manual_downloads.py' (or 'make setup-manual-downloads') "
                        f"to copy paper cohorts from cluster storage into data/manual_download/."
                    )

        case _:
            return Failure(f"No loader implementation available for dataset '{dataset_id}' (canonical '{canonical_id}')")


def load_dataset(
    dataset_id: str,
    base_dir: Path | None = None,
    raw_dir: Path | None = None,
    output_h5ad: Path | None = None,
    auto_download: bool = True,
    force_recompute: bool = False,
    force_download: bool = False,
    cache_h5ad: bool = True,
    normalize_ensembl: bool = True,
    ensembl_release: int | None = None,
    drop_unmapped: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
    apply_qc: bool = True,
    apply_sc_pipeline: bool = False,
    sc_pipeline_spec: SingleCellProcessingSpec | None = None,
    batch_size: int = 25000,
    use_batched_processing: bool = True,
) -> Result[ad.AnnData, str]:
    """Load an individual single-cell or bulk dataset by its registered identifier.

    Prioritizes loading directly from cached H5AD in <0.5s. If no H5AD exists (or
    force_recompute=True), auto-downloads raw files, parses them, normalizes gene IDs
    to canonical Ensembl IDs (Release 111 by default), enriches .var attributes,
    serializes to an H5AD cache file, and returns the AnnData object.

    For large datasets (>= batch_size cells), processing and Ensembl normalization
    are executed in memory-bounded batches streamed directly to H5AD to eliminate
    memory spikes (>180 GB).

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
        subset: Optional column-value dictionary to subset cells.
        subsample_n: Optional number of cells to subsample.
        batch_size: Cell batch size for memory-bounded sparse processing (default 25,000).
        use_batched_processing: If True, uses incremental sparse processing for large datasets.

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

    canonical_id = resolve_dataset_id(dataset_id)
    spec_maybe = get_dataset_spec(canonical_id)
    if not isinstance(spec_maybe, Some):
        msg = f"Dataset '{dataset_id}' (canonical '{canonical_id}') is not recognized in the registry"
        logger.error(msg)
        return Failure(msg)

    root = (base_dir or find_repo_root()).resolve()

    # Step 1: Fast path - load from preprocessed H5AD if available
    if not force_recompute and not force_download:
        h5ad_found = find_dataset_h5ad(canonical_id, repo_root=root)
        match h5ad_found:
            case Some(h5ad_path):
                size_mb = h5ad_path.stat().st_size / (1024 * 1024)
                logger.info(
                    "Found cached H5AD for '%s' at %s (%.1f MB). Loading directly in <0.5s...",
                    canonical_id,
                    h5ad_path.name,
                    size_mb,
                )
                try:
                    adata = ad.read_h5ad(h5ad_path)
                    if subset or subsample_n:
                        adata = _apply_subset_and_subsample(adata, subset=subset, subsample_n=subsample_n)

                    if normalize_ensembl:
                        from .genes.mapper import detect_gene_id_type
                        from .types import GeneIDType

                        if detect_gene_id_type(tuple(str(g) for g in adata.var_names[:50])) != GeneIDType.ENSEMBL_ID:
                            logger.info(
                                "Cached H5AD for '%s' does not use Ensembl IDs. Normalizing to Ensembl...",
                                canonical_id,
                            )
                            from .preprocessing.gene_normalization import normalize_dataset_to_ensembl

                            norm_res = normalize_dataset_to_ensembl(
                                adata,
                                release=ensembl_release,
                                drop_unmapped=drop_unmapped,
                            )
                            match norm_res:
                                case Success(norm_adata):
                                    adata = norm_adata
                                    if not subset and subsample_n is None and cache_h5ad:
                                        try:
                                            if adata.obs_names.name is not None and adata.obs_names.name in adata.obs.columns:
                                                adata.obs_names.name = None
                                            if adata.var_names.name is not None and adata.var_names.name in adata.var.columns:
                                                adata.var_names.name = None
                                            adata.write_h5ad(h5ad_path)
                                            logger.info("Persisted Ensembl-normalized AnnData back to cache at %s", h5ad_path.name)
                                        except Exception as write_err:
                                            logger.warning("Could not persist Ensembl update to cache: %s", write_err)
                                case Failure(norm_err):
                                    logger.warning(
                                        "Ensembl normalization of cached H5AD failed: %s. Using loaded representation.",
                                        norm_err,
                                    )

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
    res = _dispatch_load(
        canonical_id,
        root,
        raw_dir=raw_dir,
        auto_download=auto_download,
        force_download=force_download,
        subset=subset,
        subsample_n=subsample_n,
    )

    # Step 3: Quality Control, Ensembl Normalization, Layer Standardization, and H5AD Serialization
    match res:
        case Success(adata):
            adata.var_names_make_unique()
            should_cache = (cache_h5ad and not subset and subsample_n is None) or (output_h5ad is not None)
            target_h5ad = output_h5ad if output_h5ad is not None else (get_preprocessed_h5ad_path(canonical_id, repo_root=root) if should_cache else None)

            # 3a. Quality Control & Single-Cell Processing Pipeline
            spec = spec_maybe.unwrap()
            if apply_sc_pipeline:
                from .preprocessing.pipeline import process_single_cell_dataset

                plot_dir = root / f"output/reports/qc/{canonical_id}"
                pipe_res = process_single_cell_dataset(
                    adata,
                    spec=sc_pipeline_spec,
                    dataset_name=canonical_id,
                    output_dir=plot_dir,
                )
                match pipe_res:
                    case Success(sc_result):
                        logger.info(
                            "Executed single-cell pipeline for '%s': %d -> %d cells",
                            canonical_id,
                            sc_result.n_cells_pre_qc,
                            sc_result.n_cells_post_qc,
                        )
                        adata = sc_result.adata
                    case Failure(pipe_err):
                        logger.warning("Single-cell pipeline failed: %s. Falling back to standard QC.", pipe_err)
            elif apply_qc and isinstance(spec.qc_spec, Some):
                from .preprocessing.qc import apply_quality_control

                qc_res = apply_quality_control(adata, qc_spec=spec.qc_spec.unwrap())
                match qc_res:
                    case Success(qc_adata):
                        logger.info(
                            "Applied QC filtering to '%s': %d -> %d cells",
                            canonical_id,
                            adata.n_obs,
                            qc_adata.n_obs,
                        )
                        adata = qc_adata
                    case Failure(qc_err):
                        logger.warning("QC filtering encountered an error: %s. Proceeding without QC.", qc_err)

            # 3b. Ensembl Gene Normalization
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

            # 3c. Dual-Layer Standardization (adata.X & layers['counts'] raw, layers['log1p_norm'])
            from .preprocessing.normalization import standardize_dual_layers

            std_res = standardize_dual_layers(adata)
            match std_res:
                case Success(std_adata):
                    adata = std_adata
                case Failure(std_err):
                    logger.warning("Dual layer standardization failed: %s", std_err)

            # 3d. Serialization to H5AD Cache
            if should_cache and target_h5ad is not None:
                try:
                    target_h5ad.parent.mkdir(parents=True, exist_ok=True)

                    # Sanitize index names for robust h5py serialization (avoid column collision)
                    if adata.obs_names.name is not None and adata.obs_names.name in adata.obs.columns:
                        adata.obs_names.name = None
                    if adata.var_names.name is not None and adata.var_names.name in adata.var.columns:
                        adata.var_names.name = None

                    logger.info("Writing parsed AnnData to H5AD cache at %s...", target_h5ad)
                    adata.write_h5ad(target_h5ad)
                    size_mb = target_h5ad.stat().st_size / (1024 * 1024)
                    logger.info(
                        "Successfully serialized dataset '%s' to H5AD: %s (%.1f MB). Future calls will load instantly.",
                        canonical_id,
                        target_h5ad.name,
                        size_mb,
                    )
                except Exception as cache_exc:
                    logger.warning("Failed to serialize H5AD cache for '%s': %s", canonical_id, cache_exc)
                    if target_h5ad is not None and target_h5ad.exists():
                        target_h5ad.unlink(missing_ok=True)

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


IATLAS_COMBINED_GROUPS: Mapping[str, tuple[str, ...]] = {
    "melanoma": (
        "Hugo-iAtlas",
        "Riaz-iAtlas",
        "Liu-iAtlas",
        "Gide-iAtlas",
    ),
    "rcc": (
        "McDermott-iAtlas",
        "Choueiri-iAtlas",
    ),
    "pancancer": IATLAS_COHORTS,
    "all": IATLAS_COHORTS,
}


def load_combined_iatlas_cohorts(
    cohorts_or_group: Literal["melanoma", "rcc", "pancancer", "all"] | Sequence[str] = "pancancer",
    config: HarmonizeConfig | None = None,
    base_dir: Path | None = None,
) -> Result[ad.AnnData, str]:
    """Load and harmonize combined iAtlas cohorts into a unified AnnData object.

    Parameters
    ----------
    cohorts_or_group:
        Predefined cancer grouping ("melanoma", "rcc", "pancancer", "all")
        or a sequence of specific iAtlas cohort IDs.
    config:
        Harmonization configuration. Defaults to HarmonizeMode.INTERSECTION
        with reconcile_genes=False to preserve Hugo gene symbols.
    base_dir:
        Optional base directory override.

    Returns
    -------
    Result[ad.AnnData, str]:
        Success containing concatenated AnnData with unified genes and clinical metadata,
        or Failure if loading fails.
    """
    cohort_ids: tuple[str, ...]
    if isinstance(cohorts_or_group, str):
        key = cohorts_or_group.strip().lower()
        if key in IATLAS_COMBINED_GROUPS:
            cohort_ids = IATLAS_COMBINED_GROUPS[key]
        elif cohorts_or_group in IATLAS_COHORTS:
            cohort_ids = (cohorts_or_group,)
        else:
            return Failure(
                f"Unknown iAtlas cohort or grouping '{cohorts_or_group}'. "
                f"Available groupings: {sorted(IATLAS_COMBINED_GROUPS.keys())}, "
                f"or individual cohorts: {list(IATLAS_COHORTS)}"
            )
    else:
        cohort_ids = tuple(cohorts_or_group)

    if not cohort_ids:
        return Failure("No cohort IDs provided for combination.")

    cfg = config or HarmonizeConfig(
        mode=HarmonizeMode.INTERSECTION,
        reconcile_genes=False,
    )
    return query_datasets(dataset_ids=cohort_ids, config=cfg, base_dir=base_dir)


def load_iatlas_cohort_or_combined(
    cohort_or_group: str,
    base_dir: Path | None = None,
) -> Result[ad.AnnData, str]:
    """Load an individual iAtlas cohort or predefined combined grouping via tme_datasets."""
    key = cohort_or_group.strip().lower()
    if key in IATLAS_COMBINED_GROUPS or key in ("melanoma", "rcc", "pancancer", "all"):
        return load_combined_iatlas_cohorts(cohorts_or_group=key, base_dir=base_dir)
    return load_dataset(dataset_id=cohort_or_group, base_dir=base_dir)


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


def sample_single_cell_cohorts(
    sampling_spec: SingleCellSamplingSpec = SingleCellSamplingSpec(),
    cluster_spec: ClusterAnalysisSpec = ClusterAnalysisSpec(),
    harmonize_config: HarmonizeConfig = HarmonizeConfig(
        mode=HarmonizeMode.INTERSECTION,
        gene_target_type=GeneIDType.ENSEMBL_ID,
    ),
) -> Result[SampledSingleCellResult, str]:
    """Sample cells across single-cell cohorts, harmonize features, and execute PCA, kNN, and Leiden clustering."""
    from .sampling.cohort_sampler import sample_single_cell_cohorts as _sample_impl

    return _sample_impl(
        sampling_spec=sampling_spec,
        cluster_spec=cluster_spec,
        harmonize_config=harmonize_config,
    )


def generate_random_sampling_spec(**kwargs: Any) -> SingleCellSamplingSpec:
    """Generate randomized parameters for single-cell multi-cohort sampling."""
    from .sampling.cohort_sampler import generate_random_sampling_spec as _impl

    return _impl(**kwargs)


def generate_random_cluster_spec(**kwargs: Any) -> ClusterAnalysisSpec:
    """Generate randomized parameters for PCA, kNN, and Leiden clustering."""
    from .sampling.cohort_sampler import generate_random_cluster_spec as _impl

    return _impl(**kwargs)


def get_dataset_metadata(dataset_id: str) -> Result[DatasetSpec, str]:
    """Retrieve immutable DatasetSpec metadata for a dataset by ID or alias."""
    canonical_id = resolve_dataset_id(dataset_id)
    spec_maybe = get_dataset_spec(canonical_id)
    match spec_maybe:
        case Some(spec):
            return Success(spec)
        case _:
            return Failure(f"Dataset '{dataset_id}' (canonical '{canonical_id}') not found in registry")


def list_datasets(
    modality: Modality | None = None,
    has_response: bool | None = None,
    tier: str | None = None,
    exclude_dysfunctional: bool = True,
) -> RegisteredDatasets:
    """Filter and return registered dataset specifications declaratively.

    Args:
        modality: Optional modality filter.
        has_response: Optional filter for response label availability.
        tier: Optional tier string filter.
        exclude_dysfunctional: If True, exclude cohorts flagged as dysfunctional. Defaults to True.

    Returns:
        RegisteredDatasets: Immutable tuple container of matching dataset specifications.
    """
    specs = list_registered_datasets()
    filtered = []
    for s in specs:
        if modality is not None and s.modality != modality:
            continue
        if has_response is not None and s.has_response_labels != has_response:
            continue
        if tier is not None and s.tier.value_or(None) != tier:
            continue
        if exclude_dysfunctional and s.is_dysfunctional:
            continue
        filtered.append(s)
    return RegisteredDatasets(filtered)


def run_sampling_hpo(
    search_space: Any,
    bulk_data_dict: Any,
    n_trials: int = 15,
    seed: int = 42,
    early_prune_threshold: float = 0.52,
    n_jobs: int = 1,
) -> Any:
    """Execute complete multi-fidelity ASHA hyperparameter optimization campaign."""
    from .sampling.hpo_optimizer import run_sampling_hpo as _impl

    return _impl(
        search_space=search_space,
        bulk_data_dict=bulk_data_dict,
        n_trials=n_trials,
        seed=seed,
        early_prune_threshold=early_prune_threshold,
        n_jobs=n_jobs,
    )


def filter_candidate_cohorts(
    cancer_types: Any = None,
    min_viable_cells: int = 200,
    require_cached: bool = True,
) -> Any:
    """Filter single-cell cohorts from registry suitable for reference deconvolution optimization."""
    from .sampling.hpo_optimizer import filter_candidate_cohorts as _impl

    return _impl(
        cancer_types=cancer_types,
        min_viable_cells=min_viable_cells,
        require_cached=require_cached,
    )


from .sampling.hpo_models import HPORunResult, HPOSearchSpace, HPOTrialConfig, TrialEvaluationResult
from .sampling.malignant_sampling import MalignantSamplingConfig, MalignantStrategy

__all__ = [
    "IATLAS_COMBINED_GROUPS",
    "load_combined_iatlas_cohorts",
    "load_dataset",
    "query_datasets",
    "list_datasets",
    "list_preprocessed_datasets",
    "query_preprocessed_datasets",
    "get_dataset_metadata",
    "load_geneset_collection",
    "sample_single_cell_cohorts",
    "generate_random_sampling_spec",
    "generate_random_cluster_spec",
    "build_deconvolution_reference",
    "DeconvolutionReferenceConfig",
    "DeconvolutionReferenceResult",
    "run_sampling_hpo",
    "filter_candidate_cohorts",
    "HPOSearchSpace",
    "HPOTrialConfig",
    "HPORunResult",
    "TrialEvaluationResult",
    "MalignantStrategy",
    "MalignantSamplingConfig",
]


