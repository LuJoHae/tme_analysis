"""Multi-cohort single-cell sampling, gene alignment, joint PCA, and Leiden clustering."""

from __future__ import annotations

from typing import Mapping, Sequence
import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

from ..genes.reconcile import reconcile_genes
from ..harmonization.align import align_and_concatenate
from ..logging import get_logger
from ..models import (
    ClusterAnalysisSpec,
    GeneReconcileConfig,
    HarmonizeConfig,
    SampledSingleCellResult,
    SingleCellSamplingSpec,
    SubsampleSpec,
)
from ..paths import find_dataset_h5ad
from ..query import load_dataset
from ..registry import list_registered_datasets
from ..transforms.sampling import subsample_cells
from ..types import CohortSamplingMode, GeneIDType, HarmonizeMode, Modality

logger = get_logger("sampling.cohort_sampler")


def resolve_sampled_cohorts(
    spec: SingleCellSamplingSpec,
) -> Result[tuple[str, ...], str]:
    """Resolve and deterministically select single-cell cohort IDs according to sampling criteria.

    Args:
        spec: SingleCellSamplingSpec configuration.

    Returns:
        Success(tuple of dataset IDs) or Failure(error message).
    """
    try:
        registered = list_registered_datasets().filter(modality=Modality.SINGLE_CELL)

        # 1. If explicit cohort IDs provided, validate and filter
        match spec.cohort_ids:
            case Some(explicit_ids):
                valid_ids: list[str] = []
                for cid in explicit_ids:
                    if cid in registered:
                        valid_ids.append(cid)
                    else:
                        logger.warning("Requested cohort ID '%s' not registered as single-cell; skipping.", cid)
                candidates = valid_ids
            case _:
                candidates = [s.id for s in registered]

        # 2. Filter by cancer type if requested
        match spec.cancer_types:
            case Some(cancer_types):
                cancer_types_lower = {ct.lower() for ct in cancer_types}
                candidates = [
                    cid
                    for cid in candidates
                    if registered[cid].cancer_type.lower() in cancer_types_lower
                ]
            case _:
                pass

        # 3. Filter by response label presence if requested
        match spec.has_response_labels:
            case Some(has_resp):
                candidates = [
                    cid
                    for cid in candidates
                    if registered[cid].has_response_labels == has_resp
                ]
            case _:
                pass

        # 4. Filter by disk availability if require_cached_h5ad is True
        if spec.require_cached_h5ad:
            available_on_disk = [
                cid
                for cid in candidates
                if isinstance(find_dataset_h5ad(cid), Some)
            ]
            if not available_on_disk:
                return Failure(
                    f"No cached processed single-cell H5AD files found on disk for candidates: {candidates[:5]}..."
                )
            candidates = available_on_disk

        if not candidates:
            return Failure("No single-cell cohorts matched the specified sampling criteria.")

        # 5. Deterministic random cohort sampling if n_cohorts is specified
        match spec.n_cohorts:
            case Some(k):
                if k > len(candidates):
                    logger.warning(
                        "Requested %d cohorts, but only %d matched criteria. Using all %d.",
                        k,
                        len(candidates),
                        len(candidates),
                    )
                    selected_cohorts = tuple(candidates)
                else:
                    seed_val = spec.seed.value_or(None) if isinstance(spec.seed, Some) else None
                    rng = np.random.default_rng(seed_val)
                    chosen_idx = rng.choice(len(candidates), size=k, replace=False)
                    selected_cohorts = tuple(candidates[i] for i in sorted(chosen_idx))
            case _:
                selected_cohorts = tuple(candidates)

        logger.info(
            "Resolved %d single-cell cohorts for sampling: %s",
            len(selected_cohorts),
            list(selected_cohorts),
        )
        return Success(selected_cohorts)
    except Exception as exc:
        msg = f"Failed to resolve sampled cohorts: {exc}"
        logger.error(msg)
        return Failure(msg)


def compute_cohort_cell_allocations(
    cohort_ids: Sequence[str],
    cohort_total_cells: Mapping[str, int],
    spec: SingleCellSamplingSpec,
) -> Mapping[str, int]:
    """Compute the number of cells to sample from each selected cohort.

    Args:
        cohort_ids: Sequence of selected cohort IDs.
        cohort_total_cells: Mapping of cohort ID to total available cells in that cohort.
        spec: SingleCellSamplingSpec configuration.

    Returns:
        Mapping of cohort ID to target number of sampled cells.
    """
    allocations: dict[str, int] = {}

    match spec.mode:
        case CohortSamplingMode.FIXED_PER_COHORT:
            for cid in cohort_ids:
                total = cohort_total_cells.get(cid, spec.n_cells_per_cohort)
                allocations[cid] = max(1, min(total, spec.n_cells_per_cohort))

        case CohortSamplingMode.FRACTION_PER_COHORT:
            for cid in cohort_ids:
                total = cohort_total_cells.get(cid, 1000)
                alloc = int(round(total * spec.fraction_per_cohort))
                allocations[cid] = max(1, min(total, alloc))

        case CohortSamplingMode.EXPLICIT_COUNTS:
            for cid in cohort_ids:
                total = cohort_total_cells.get(cid, spec.n_cells_per_cohort)
                requested = spec.explicit_cell_counts.get(cid, spec.n_cells_per_cohort)
                allocations[cid] = max(1, min(total, requested))

        case CohortSamplingMode.EXPLICIT_FRACTIONS:
            for cid in cohort_ids:
                total = cohort_total_cells.get(cid, 1000)
                fraction = spec.explicit_cell_fractions.get(cid, spec.fraction_per_cohort)
                alloc = int(round(total * fraction))
                allocations[cid] = max(1, min(total, alloc))

        case CohortSamplingMode.GLOBAL_BUDGET:
            budget = spec.global_cell_budget.value_or(len(cohort_ids) * spec.n_cells_per_cohort)
            sum_total_cells = sum(cohort_total_cells.get(cid, 1000) for cid in cohort_ids)
            if sum_total_cells == 0:
                sum_total_cells = 1

            for cid in cohort_ids:
                total = cohort_total_cells.get(cid, 1000)
                # Proportional allocation of global budget
                prop = total / sum_total_cells
                alloc = max(1, min(total, int(round(budget * prop))))
                allocations[cid] = alloc

    return allocations


def sample_and_harmonize_cohorts(
    cohort_ids: Sequence[str],
    cell_allocations: Mapping[str, int],
    spec: SingleCellSamplingSpec,
    harmonize_config: HarmonizeConfig,
) -> Result[ad.AnnData, str]:
    """Load, subsample, and harmonize multiple single-cell cohorts into a unified AnnData object.

    Args:
        cohort_ids: Sequence of dataset IDs to load and sample.
        cell_allocations: Target number of cells for each cohort.
        spec: SingleCellSamplingSpec with stratification and seed parameters.
        harmonize_config: HarmonizeConfig defining intersection or union feature alignment.

    Returns:
        Success(harmonized_adata) or Failure(error message).
    """
    sampled_adatas: list[ad.AnnData] = []
    loaded_ids: list[str] = []

    seed_val = spec.seed.value_or(42) if isinstance(spec.seed, Some) else 42

    for idx, cid in enumerate(cohort_ids):
        target_n = cell_allocations.get(cid, spec.n_cells_per_cohort)
        logger.info("[%s] Loading dataset for sampling (%d cells requested)...", cid, target_n)

        load_res = load_dataset(cid, auto_download=False)
        match load_res:
            case Failure(err):
                logger.warning("[%s] Could not load cohort: %s. Skipping.", cid, err)
                continue
            case Success(raw_adata):
                cohort_adata = raw_adata

        # 1. Quality Control: Cell-level filtering
        if spec.only_qc_passing_cells:
            n_cells_pre = cohort_adata.n_obs
            if "is_retained_qc" in cohort_adata.obs.columns:
                cohort_adata = cohort_adata[cohort_adata.obs["is_retained_qc"].astype(bool)].copy()
                logger.info(
                    "[%s] Cell QC filtering via 'is_retained_qc': %d -> %d viable cells retained",
                    cid,
                    n_cells_pre,
                    cohort_adata.n_obs,
                )
            elif "pass_qc" in cohort_adata.obs.columns:
                cohort_adata = cohort_adata[cohort_adata.obs["pass_qc"].astype(bool)].copy()
                logger.info(
                    "[%s] Cell QC filtering via 'pass_qc': %d -> %d viable cells retained",
                    cid,
                    n_cells_pre,
                    cohort_adata.n_obs,
                )
            elif "sc_processing_summary" not in cohort_adata.uns:
                logger.info("[%s] Applying adaptive QC filtering to cohort before sampling...", cid)
                from ..preprocessing.qc import apply_quality_control

                qc_spec_to_use = spec.qc_spec.value_or(None) if isinstance(spec.qc_spec, Some) else None
                qc_res = apply_quality_control(
                    cohort_adata,
                    qc_spec=qc_spec_to_use,
                    use_adaptive_qc=(qc_spec_to_use is None),
                    min_cells_per_gene=spec.min_cells_per_gene,
                )
                match qc_res:
                    case Success(qc_adata):
                        cohort_adata = qc_adata
                    case Failure(err):
                        logger.warning("[%s] On-the-fly QC filtering failed (%s); proceeding.", cid, err)

        if cohort_adata.n_obs == 0:
            logger.warning("[%s] 0 cells passed QC filtering! Skipping cohort.", cid)
            continue

        # 2. Quality Control: Gene-level filtering
        if spec.only_qc_passing_genes:
            n_genes_pre = cohort_adata.n_vars
            if "is_retained_qc" in cohort_adata.var.columns:
                cohort_adata = cohort_adata[:, cohort_adata.var["is_retained_qc"].astype(bool)].copy()
                logger.info(
                    "[%s] Gene QC filtering via 'is_retained_qc': %d -> %d genes retained",
                    cid,
                    n_genes_pre,
                    cohort_adata.n_vars,
                )
            elif "pass_qc" in cohort_adata.var.columns:
                cohort_adata = cohort_adata[:, cohort_adata.var["pass_qc"].astype(bool)].copy()
                logger.info(
                    "[%s] Gene QC filtering via 'pass_qc': %d -> %d genes retained",
                    cid,
                    n_genes_pre,
                    cohort_adata.n_vars,
                )

            # Enforce minimum expression threshold across viable cells
            if spec.min_cells_per_gene > 0:
                X_mat = cohort_adata.X
                if sp.issparse(X_mat):
                    cells_per_gene = np.asarray((X_mat > 0).sum(axis=0)).ravel()
                else:
                    cells_per_gene = np.asarray((np.asarray(X_mat) > 0).sum(axis=0)).ravel()
                pass_gene_mask = cells_per_gene >= spec.min_cells_per_gene
                cohort_adata.var["cells_per_gene"] = cells_per_gene
                cohort_adata.var["is_retained_qc"] = pass_gene_mask
                cohort_adata = cohort_adata[:, pass_gene_mask].copy()
                logger.info(
                    "[%s] Filtered genes (expressed in >=%d cells): %d -> %d genes retained",
                    cid,
                    spec.min_cells_per_gene,
                    n_genes_pre,
                    cohort_adata.n_vars,
                )

        if cohort_adata.n_vars == 0:
            logger.warning("[%s] 0 genes passed QC filtering! Skipping cohort.", cid)
            continue

        # Subsample cells from this cohort
        subsample_spec = SubsampleSpec(
            n_or_fraction=float(target_n),
            stratify_by=spec.stratify_by,
            balanced=spec.balanced_strata,
            seed=Some(seed_val + idx),
        )
        sample_res = subsample_cells(cohort_adata, subsample_spec)
        match sample_res:
            case Failure(err):
                logger.warning("[%s] Failed to subsample cells: %s. Using first %d cells.", cid, err, target_n)
                sub_adata = cohort_adata[:target_n].copy()
            case Success(sampled):
                sub_adata = sampled

        # Annotate provenance metadata in .obs
        sub_adata.obs["dataset_id"] = cid
        sub_adata.obs_names = [f"{cid}_{name}" for name in sub_adata.obs_names]
        sub_adata.obs_names_make_unique()

        # Reconcile gene features to target nomenclature (e.g. HUGO symbols) if enabled
        if harmonize_config.reconcile_genes:
            reconcile_res = reconcile_genes(
                sub_adata,
                GeneReconcileConfig(target_type=harmonize_config.gene_target_type),
            )
            match reconcile_res:
                case Success(rec_adata):
                    sub_adata = rec_adata
                case Failure(err):
                    logger.warning("[%s] Gene reconciliation warning: %s", cid, err)

        sampled_adatas.append(sub_adata)
        loaded_ids.append(cid)

    if not sampled_adatas:
        return Failure("Failed to load and sample any of the requested cohorts.")

    logger.info(
        "Harmonizing %d sampled cohorts (%s) via %s...",
        len(sampled_adatas),
        loaded_ids,
        harmonize_config.mode.value,
    )
    return align_and_concatenate(sampled_adatas, loaded_ids, harmonize_config)


def run_pca_knn_leiden(
    combined_adata: ad.AnnData,
    cluster_spec: ClusterAnalysisSpec,
) -> Result[ad.AnnData, str]:
    """Execute Highly Variable Gene selection, joint PCA, kNN graph construction, and Leiden clustering.

    Args:
        combined_adata: Harmonized AnnData with sampled cells from multiple cohorts.
        cluster_spec: ClusterAnalysisSpec parameter configuration.

    Returns:
        Success(clustered_adata) or Failure(error message).
    """
    try:
        adata = combined_adata.copy()
        n_obs = adata.n_obs
        n_vars = adata.n_vars

        if n_obs < 3:
            return Failure(f"Cannot run PCA and clustering on fewer than 3 cells (received {n_obs})")

        # 1. Ensure log-transformed expression for PCA
        # If dataset already has log1p or normalized values in .X, use as is.
        # Check max value to verify log scale: if max > 50, likely raw counts; apply log1p.
        import scipy.sparse as sp
        max_val = adata.X.max() if not sp.issparse(adata.X) else adata.X.data.max() if len(adata.X.data) > 0 else 0
        if max_val > 50.0:
            logger.info("Counts in .X appear unlogged (max=%.1f); applying log1p prior to PCA...", max_val)
            sc.pp.log1p(adata)

        # 2. Highly Variable Gene (HVG) Selection
        hvg_target = min(cluster_spec.n_top_genes, n_vars)
        if hvg_target < n_vars:
            logger.info("Selecting top %d Highly Variable Genes across cohorts...", hvg_target)
            batch_key = (
                cluster_spec.batch_key
                if (
                    cluster_spec.batch_key in adata.obs.columns
                    and len(adata.obs[cluster_spec.batch_key].unique()) > 1
                )
                else None
            )
            try:
                sc.pp.highly_variable_genes(
                    adata,
                    n_top_genes=hvg_target,
                    flavor="seurat",
                    batch_key=batch_key,
                    subset=False,
                )
            except Exception as hvg_err:
                logger.warning("Batch-aware HVG selection failed (%s); falling back to global HVG.", hvg_err)
                sc.pp.highly_variable_genes(
                    adata,
                    n_top_genes=hvg_target,
                    flavor="seurat",
                    subset=False,
                )
            use_hvg = True
        else:
            use_hvg = False

        # 3. Joint PCA
        n_pcs = min(cluster_spec.n_pcs, n_obs - 1, n_vars - 1)
        n_pcs = max(2, n_pcs)
        logger.info("Computing joint PCA (%d principal components)...", n_pcs)
        sc.pp.pca(
            adata,
            n_comps=n_pcs,
            mask_var="highly_variable" if use_hvg else None,
            zero_center=True,
            svd_solver="arpack",
            random_state=cluster_spec.seed.value_or(42) if isinstance(cluster_spec.seed, Some) else 42,
        )

        # 4. k-Nearest Neighbors (kNN) Graph Construction
        n_neighbors = min(cluster_spec.n_neighbors, n_obs - 1)
        n_neighbors = max(2, n_neighbors)
        logger.info("Constructing kNN graph (k=%d, metric=%s)...", n_neighbors, cluster_spec.metric)
        sc.pp.neighbors(
            adata,
            n_neighbors=n_neighbors,
            n_pcs=min(n_pcs, adata.obsm["X_pca"].shape[1]),
            metric=cluster_spec.metric,
            use_rep="X_pca",
            random_state=cluster_spec.seed.value_or(42) if isinstance(cluster_spec.seed, Some) else 42,
        )

        # 5. Leiden Community Detection
        logger.info("Running Leiden clustering (resolution=%.2f)...", cluster_spec.leiden_resolution)
        try:
            sc.tl.leiden(
                adata,
                resolution=cluster_spec.leiden_resolution,
                key_added=cluster_spec.key_added,
                flavor="igraph",
                n_iterations=2,
                directed=False,
                random_state=cluster_spec.seed.value_or(42) if isinstance(cluster_spec.seed, Some) else 42,
            )
        except Exception:
            # Fallback for standard leidenalg
            sc.tl.leiden(
                adata,
                resolution=cluster_spec.leiden_resolution,
                key_added=cluster_spec.key_added,
                random_state=cluster_spec.seed.value_or(42) if isinstance(cluster_spec.seed, Some) else 42,
            )

        n_clusters = len(adata.obs[cluster_spec.key_added].unique())
        logger.info("Leiden clustering identified %d clusters.", n_clusters)
        return Success(adata)
    except Exception as exc:
        msg = f"Failed to compute PCA, kNN, and Leiden clustering: {exc}"
        logger.error(msg)
        return Failure(msg)


def sample_single_cell_cohorts(
    sampling_spec: SingleCellSamplingSpec = SingleCellSamplingSpec(),
    cluster_spec: ClusterAnalysisSpec = ClusterAnalysisSpec(),
    harmonize_config: HarmonizeConfig = HarmonizeConfig(
        mode=HarmonizeMode.INTERSECTION,
        gene_target_type=GeneIDType.ENSEMBL_ID,
    ),
) -> Result[SampledSingleCellResult, str]:
    """Sample cells across single-cell cohorts, harmonize features, and execute PCA, kNN, and Leiden.

    Pipeline Workflow:
        1. Resolve and filter candidate single-cell datasets (by cancer type, availability, random K).
        2. Inspect cohort sizes and calculate target cell allocations per cohort.
        3. Subsample cells from each cohort (uniform or stratified) and concatenate on shared genes.
        4. Detect Highly Variable Genes (HVGs) and compute joint PCA across all pooled cells.
        5. Build kNN graph in PCA representation and partition with Leiden graph clustering.
        6. Package result into an immutable SampledSingleCellResult model.

    Args:
        sampling_spec: SingleCellSamplingSpec with cohort selection, sampling mode, and cell counts.
        cluster_spec: ClusterAnalysisSpec with PCA, kNN, and Leiden clustering parameters.
        harmonize_config: HarmonizeConfig (default INTERSECTION across shared gene symbols).

    Returns:
        Success(SampledSingleCellResult) or Failure(error message).
    """
    logger.info("Starting multi-cohort single-cell sampling pipeline...")

    # Step 1: Resolve Cohorts
    cohorts_res = resolve_sampled_cohorts(sampling_spec)
    match cohorts_res:
        case Failure(err):
            return Failure(err)
        case Success(cohort_ids):
            pass

    # Step 2: Query Cohort Sizes and Compute Allocations
    registered = list_registered_datasets()
    cohort_sizes: dict[str, int] = {}
    for cid in cohort_ids:
        spec_item = registered[cid]
        size = spec_item.n_samples_or_cells.value_or(2000) if isinstance(spec_item.n_samples_or_cells, Some) else 2000
        cohort_sizes[cid] = size

    allocations = compute_cohort_cell_allocations(cohort_ids, cohort_sizes, sampling_spec)
    logger.info("Planned cell sampling allocations: %s", allocations)

    # Step 3: Subsample and Harmonize
    harmonized_res = sample_and_harmonize_cohorts(
        cohort_ids=cohort_ids,
        cell_allocations=allocations,
        spec=sampling_spec,
        harmonize_config=harmonize_config,
    )
    match harmonized_res:
        case Failure(err):
            return Failure(err)
        case Success(combined_adata):
            pass

    # Step 4: PCA, kNN, and Leiden Clustering
    cluster_res = run_pca_knn_leiden(combined_adata, cluster_spec)
    match cluster_res:
        case Failure(err):
            return Failure(err)
        case Success(clustered_adata):
            pass

    # Step 5: Packaging Result
    actual_cells_per_cohort = dict(clustered_adata.obs["dataset_id"].value_counts())
    n_clusters = len(clustered_adata.obs[cluster_spec.key_added].unique())
    total_cells = clustered_adata.n_obs
    n_shared_genes = clustered_adata.n_vars

    summary = {
        "sampled_cohort_ids": list(cohort_ids),
        "total_cells": total_cells,
        "cells_per_cohort": actual_cells_per_cohort,
        "n_shared_genes": n_shared_genes,
        "n_clusters": n_clusters,
        "n_pcs": clustered_adata.obsm["X_pca"].shape[1],
        "leiden_resolution": cluster_spec.leiden_resolution,
        "harmonize_mode": harmonize_config.mode.value,
    }
    clustered_adata.uns["sampling_summary"] = summary

    logger.info(
        "Successfully sampled %d cells across %d cohorts into %d Leiden clusters.",
        total_cells,
        len(cohort_ids),
        n_clusters,
    )

    return Success(
        SampledSingleCellResult(
            adata=clustered_adata,
            sampled_cohort_ids=cohort_ids,
            cells_per_cohort=actual_cells_per_cohort,
            total_cells=total_cells,
            n_shared_genes=n_shared_genes,
            n_clusters=n_clusters,
            summary=summary,
        )
    )


def generate_random_sampling_spec(
    available_cohorts: Sequence[str] | None = None,
    mode: CohortSamplingMode | str | None = None,
    min_cohorts: int = 2,
    max_cohorts: int = 5,
    min_cells_per_cohort: int = 50,
    max_cells_per_cohort: int = 1000,
    min_fraction: float = 0.02,
    max_fraction: float = 0.25,
    min_global_budget: int = 500,
    max_global_budget: int = 5000,
    stratify_options: Sequence[str | None] = (None, "cell_type", "response"),
    require_cached_h5ad: bool = True,
    seed: int | None = None,
) -> SingleCellSamplingSpec:
    """Generate randomized parameters for single-cell multi-cohort sampling.

    Generates valid, reproducible sampling specifications supporting:
        - Cohort selection: Random K cohorts from available cached or registered cohorts.
        - Sampling mode: Randomly selects from all 5 strategies if mode is None:
            * FIXED_PER_COHORT: Random N in [min_cells_per_cohort, max_cells_per_cohort].
            * FRACTION_PER_COHORT: Random fraction in [min_fraction, max_fraction].
            * EXPLICIT_COUNTS: Randomized count dictionary {cid: random_int} for each cohort.
            * EXPLICIT_FRACTIONS: Randomized fraction dictionary {cid: random_float} for each cohort.
            * GLOBAL_BUDGET: Random total budget in [min_global_budget, max_global_budget].
        - Stratification: Random biological stratification configuration.

    Args:
        available_cohorts: Optional pool of cohort IDs to choose from. Defaults to registered sc cohorts.
        mode: Optional forced mode (or None to randomly pick a mode).
        min_cohorts: Minimum number of cohorts to select.
        max_cohorts: Maximum number of cohorts to select.
        min_cells_per_cohort: Lower bound for cell count.
        max_cells_per_cohort: Upper bound for cell count.
        min_fraction: Lower bound for fractional sampling.
        max_fraction: Upper bound for fractional sampling.
        min_global_budget: Lower bound for total cell budget.
        max_global_budget: Upper bound for total cell budget.
        stratify_options: Candidate stratification columns to pick from.
        require_cached_h5ad: If True, only select from cohorts with cached H5AD on disk.
        seed: Random seed for reproducibility.

    Returns:
        An immutable SingleCellSamplingSpec instance.
    """
    rng = np.random.default_rng(seed)

    # 1. Resolve candidate cohorts pool
    if available_cohorts is not None:
        cohort_pool = list(available_cohorts)
    else:
        registered = list_registered_datasets().filter(modality=Modality.SINGLE_CELL)
        candidates = [s.id for s in registered]
        if require_cached_h5ad:
            cached = [c for c in candidates if isinstance(find_dataset_h5ad(c), Some)]
            cohort_pool = cached if cached else candidates
        else:
            cohort_pool = candidates

    if not cohort_pool:
        cohort_pool = ["GSE120575", "Maynard_NSCLC"]

    # 2. Select random number of cohorts K
    k = min(len(cohort_pool), int(rng.integers(min_cohorts, max_cohorts + 1)))
    chosen_idx = rng.choice(len(cohort_pool), size=k, replace=False)
    selected_cohorts = tuple(str(cohort_pool[i]) for i in chosen_idx)

    # 3. Determine sampling mode
    if mode is not None:
        chosen_mode = CohortSamplingMode(mode) if isinstance(mode, str) else mode
    else:
        modes = list(CohortSamplingMode)
        mode_idx = int(rng.integers(0, len(modes)))
        chosen_mode = modes[mode_idx]

    # 4. Generate mode-specific parameters
    n_cells_per_cohort = int(rng.integers(min_cells_per_cohort, max_cells_per_cohort + 1))
    fraction_per_cohort = float(np.round(rng.uniform(min_fraction, max_fraction), 3))

    explicit_cell_counts: dict[str, int] = {}
    explicit_cell_fractions: dict[str, float] = {}
    global_cell_budget: Maybe[int] = Nothing

    match chosen_mode:
        case CohortSamplingMode.FIXED_PER_COHORT:
            pass  # uses n_cells_per_cohort
        case CohortSamplingMode.FRACTION_PER_COHORT:
            pass  # uses fraction_per_cohort
        case CohortSamplingMode.EXPLICIT_COUNTS:
            for cid in selected_cohorts:
                explicit_cell_counts[cid] = int(rng.integers(min_cells_per_cohort, max_cells_per_cohort + 1))
        case CohortSamplingMode.EXPLICIT_FRACTIONS:
            for cid in selected_cohorts:
                explicit_cell_fractions[cid] = float(np.round(rng.uniform(min_fraction, max_fraction), 3))
        case CohortSamplingMode.GLOBAL_BUDGET:
            budget = int(rng.integers(min_global_budget, max_global_budget + 1))
            global_cell_budget = Some(budget)

    # 5. Stratification
    strat_idx = int(rng.integers(0, len(stratify_options)))
    strat_col = stratify_options[strat_idx]
    stratify_by = Some(str(strat_col)) if strat_col is not None else Nothing
    balanced_strata = bool(rng.integers(0, 2)) if stratify_by is not Nothing else False

    child_seed = int(rng.integers(1, 1000000)) if seed is not None else None

    return SingleCellSamplingSpec(
        cohort_ids=Some(selected_cohorts),
        mode=chosen_mode,
        n_cells_per_cohort=n_cells_per_cohort,
        fraction_per_cohort=fraction_per_cohort,
        explicit_cell_counts=explicit_cell_counts,
        explicit_cell_fractions=explicit_cell_fractions,
        global_cell_budget=global_cell_budget,
        stratify_by=stratify_by,
        balanced_strata=balanced_strata,
        require_cached_h5ad=require_cached_h5ad,
        seed=Some(child_seed) if child_seed is not None else Nothing,
    )


def generate_random_cluster_spec(
    min_pcs: int = 10,
    max_pcs: int = 50,
    min_neighbors: int = 5,
    max_neighbors: int = 30,
    min_resolution: float = 0.3,
    max_resolution: float = 1.5,
    seed: int | None = None,
) -> ClusterAnalysisSpec:
    """Generate randomized parameters for PCA, kNN, and Leiden clustering."""
    rng = np.random.default_rng(seed)
    n_pcs = int(rng.integers(min_pcs, max_pcs + 1))
    n_neighbors = int(rng.integers(min_neighbors, max_neighbors + 1))
    resolution = float(np.round(rng.uniform(min_resolution, max_resolution), 2))
    n_top_genes = int(rng.choice([500, 1000, 1500, 2000, 3000]))
    child_seed = int(rng.integers(1, 1000000)) if seed is not None else None

    return ClusterAnalysisSpec(
        n_top_genes=n_top_genes,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        leiden_resolution=resolution,
        seed=Some(child_seed) if child_seed is not None else Nothing,
    )

