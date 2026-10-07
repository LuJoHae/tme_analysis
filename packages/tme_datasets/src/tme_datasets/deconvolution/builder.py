"""Pure functional single-cell deconvolution reference construction engine."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Sequence
import anndata as ad
import numpy as np
import polars as pl
import scipy.sparse as sp
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

from ..logging import get_logger
from ..preprocessing.gene_filtering import filter_confounding_genes
from ..preprocessing.matrix_inspection import ExpressionType, inspect_expression_type
from .malignant import detect_malignant_cells
from .models import DeconvolutionReferenceConfig, DeconvolutionReferenceResult

logger = get_logger("deconvolution.builder")


def _compute_mrna_scaling(
    lib_sizes: np.ndarray,
    obs_states: np.ndarray,
    obs_types: np.ndarray,
    state_to_malignant: dict[str, bool],
) -> pl.DataFrame:
    """Compute cell-state library size and relative mRNA content scaling factors."""
    df_raw = pl.DataFrame({
        "cell_state": obs_states,
        "cell_type": obs_types,
        "total_umi": lib_sizes,
    })

    scaling_df = (
        df_raw.group_by(["cell_state", "cell_type"])
        .agg([
            pl.col("total_umi").mean().alias("mean_umi_per_cell"),
            pl.col("total_umi").median().alias("median_umi_per_cell"),
            pl.len().alias("n_cells"),
        ])
        .sort("cell_state")
    )

    grand_mean = float(scaling_df["mean_umi_per_cell"].mean())
    grand_mean_safe = grand_mean if grand_mean > 0 else 1.0

    return scaling_df.with_columns([
        (pl.col("mean_umi_per_cell") / grand_mean_safe).alias("relative_rna_content"),
        pl.col("cell_state")
        .map_elements(lambda s: state_to_malignant.get(str(s), False), return_dtype=pl.Boolean)
        .alias("is_malignant"),
    ])


def _compute_gep_matrix(
    X: sp.spmatrix | np.ndarray,
    group_labels: np.ndarray,
    unique_groups: list[str],
    gene_names: list[str],
    id_col: str,
    normalize: bool,
    pseudo_min: float,
) -> tuple[pl.DataFrame, np.ndarray]:
    """Compute average expression profiles and multinomial emission probabilities."""
    n_groups = len(unique_groups)
    n_genes = len(gene_names)
    profiles = np.zeros((n_groups, n_genes), dtype=np.float64)

    X_csr = X.tocsr() if sp.issparse(X) else X

    for i, grp in enumerate(unique_groups):
        mask = group_labels == grp
        sub = X_csr[mask]
        if sp.issparse(sub):
            profiles[i, :] = np.asarray(sub.mean(axis=0)).flatten()
        else:
            profiles[i, :] = np.mean(sub, axis=0)

    if normalize:
        row_sums = profiles.sum(axis=1, keepdims=True)
        row_sums = np.where(row_sums == 0, 1.0, row_sums)
        # Apply pseudo_min adjustment: phi = (ref / row_sums) * (1 - pseudo_min * G) + pseudo_min
        if pseudo_min > 0 and (1.0 - pseudo_min * n_genes) > 0:
            profiles = (profiles / row_sums) * (1.0 - pseudo_min * n_genes) + pseudo_min
        else:
            profiles = profiles / row_sums

    # Convert to Polars DataFrame
    col_dict: dict[str, list[object]] = {id_col: list(unique_groups)}
    for j, g in enumerate(gene_names):
        col_dict[g] = profiles[:, j].tolist()

    return pl.DataFrame(col_dict), profiles


def _evaluate_collinearity(
    profiles: np.ndarray,
    states: list[str],
    threshold: float,
) -> pl.DataFrame:
    """Evaluate pairwise Pearson correlation between cell states and flag collinear pairs."""
    means = profiles.mean(axis=1, keepdims=True)
    stds = profiles.std(axis=1, keepdims=True) + 1e-12
    normed = (profiles - means) / stds

    corr_mat = np.dot(normed, normed.T) / profiles.shape[1]

    pairs: list[dict[str, object]] = []
    n_states = len(states)
    for i in range(n_states):
        for j in range(i + 1, n_states):
            r_val = float(corr_mat[i, j])
            pairs.append({
                "state_a": states[i],
                "state_b": states[j],
                "correlation": r_val,
                "is_collinear": r_val >= threshold,
            })

    if not pairs:
        return pl.DataFrame(
            schema={
                "state_a": pl.String,
                "state_b": pl.String,
                "correlation": pl.Float64,
                "is_collinear": pl.Boolean,
            }
        )

    return pl.DataFrame(pairs).sort("correlation", descending=True)


def build_deconvolution_reference(
    data_source: ad.AnnData | Sequence[str],
    config: DeconvolutionReferenceConfig = DeconvolutionReferenceConfig(),
    repo_root: Path | None = None,
) -> Result[DeconvolutionReferenceResult, str]:
    """Build a publication-grade deconvolution reference from a single-cell dataset or cohort IDs.

    Args:
        data_source: Either an existing loaded AnnData object or a sequence of registered cohort IDs.
        config: DeconvolutionReferenceConfig specifying clustering keys, malignant handling, and thresholds.
        repo_root: Optional repository root path override.

    Returns:
        Success(DeconvolutionReferenceResult) with calibrated Mean GEP matrices and diagnostics,
        or Failure(error_message).
    """
    start_time = time.time()

    # 1. Dispatch data loading if cohort IDs were provided
    if isinstance(data_source, (list, tuple)):
        cohort_ids = list(data_source)
        if not cohort_ids:
            return Failure("Empty cohort IDs sequence provided to build_deconvolution_reference")

        logger.info("Harmonizing %d single-cell cohorts for reference construction: %s", len(cohort_ids), cohort_ids)
        from ..sampling.cohort_sampler import sample_single_cell_cohorts
        from ..models import SingleCellSamplingSpec, CohortSamplingMode

        sample_res = sample_single_cell_cohorts(
            sampling_spec=SingleCellSamplingSpec(
                cohort_ids=Some(tuple(cohort_ids)),
                mode=CohortSamplingMode.FRACTION_PER_COHORT,
                fraction_per_cohort=1.0,  # Retain all cells
                require_cached_h5ad=True,
            )
        )
        match sample_res:
            case Success(sampled_obj):
                working_adata = sampled_obj.adata
            case Failure(err):
                return Failure(f"Failed to load cohorts for reference building: {err}")
    elif isinstance(data_source, ad.AnnData):
        working_adata = data_source.copy()
    else:
        return Failure(f"Unsupported data_source type: {type(data_source)}. Expected AnnData or Sequence[str].")

    # 2. Validate cell_state key
    resolved_state_key = config.cell_state_key
    if resolved_state_key not in working_adata.obs.columns:
        if resolved_state_key == "cell_state":
            fallback_candidates = ["leiden", "cell_type", "cell_subtype"]
            found_fallback = next((k for k in fallback_candidates if k in working_adata.obs.columns), None)
            if found_fallback is not None:
                logger.info("Configured default cell_state_key 'cell_state' not in obs; falling back to '%s'", found_fallback)
                resolved_state_key = found_fallback
            else:
                return Failure(
                    f"Cell state column '{config.cell_state_key}' not found in AnnData.obs. "
                    f"Available columns: {list(working_adata.obs.columns)}"
                )
        else:
            return Failure(
                f"Cell state column '{config.cell_state_key}' not found in AnnData.obs. "
                f"Available columns: {list(working_adata.obs.columns)}"
            )

    # 3. Ensure linear expression counts
    inspection = inspect_expression_type(working_adata)
    match inspection.expression_type:
        case ExpressionType.LOG_NORMALIZED:
            logger.info("Detected log-normalized expression. Converting back to linear space...")
            X_arr = working_adata.X.toarray() if sp.issparse(working_adata.X) else np.asarray(working_adata.X)
            linear_X = np.expm1(X_arr)
            working_adata.X = sp.csr_matrix(linear_X) if sp.issparse(working_adata.X) else linear_X
        case _:
            pass

    # 4. Filter confounding / uninformative gene families
    if config.filter_confounding:
        filter_res = filter_confounding_genes(working_adata)
        match filter_res:
            case Success(filtered_adata):
                working_adata = filtered_adata
            case Failure(err):
                logger.warning("Confounding gene filtering returned failure: %s. Continuing with original genes.", err)

    # 5. Filter cell states below min_cells_per_state
    obs_states_raw = working_adata.obs[resolved_state_key].astype(str).to_numpy()
    unique_raw, counts_raw = np.unique(obs_states_raw, return_counts=True)
    valid_states = {s for s, count in zip(unique_raw, counts_raw) if count >= config.min_cells_per_state}

    if not valid_states:
        return Failure(
            f"No cell states met the minimum threshold of {config.min_cells_per_state} cells. "
            f"State counts: {dict(zip(unique_raw, counts_raw))}"
        )

    state_mask = np.isin(obs_states_raw, list(valid_states))
    if not np.all(state_mask):
        dropped_count = np.sum(~state_mask)
        logger.info("Filtered %d cells belonging to rare states (< %d cells)", dropped_count, config.min_cells_per_state)
        working_adata = working_adata[state_mask].copy()

    obs_states = working_adata.obs[resolved_state_key].astype(str).to_numpy()
    unique_states = sorted(list(set(obs_states)))

    # 6. Extract or map broad cell types (two-tier hierarchy)
    match config.cell_type_key:
        case config.cell_type_key if config.cell_type_key.value_or(None) is not None:
            type_col = config.cell_type_key.unwrap()
            if type_col in working_adata.obs.columns:
                obs_types = working_adata.obs[type_col].astype(str).to_numpy()
            else:
                logger.warning("Cell type key '%s' not found; defaulting to cell state labels.", type_col)
                obs_types = obs_states
        case _:
            obs_types = obs_states

    unique_types = sorted(list(set(obs_types)))

    # 7. Identify malignant cells and classify states
    is_mal_cells = detect_malignant_cells(working_adata, config)
    state_to_malignant: dict[str, bool] = {}
    for st in unique_states:
        st_mask = obs_states == st
        # If >= 50% of cells in cluster are malignant, state is considered malignant
        mal_fraction = float(np.mean(is_mal_cells[st_mask])) if np.any(st_mask) else 0.0
        state_to_malignant[st] = mal_fraction >= 0.5

    malignant_states = tuple(st for st, is_m in state_to_malignant.items() if is_m)

    # 8. Compute hierarchy table
    state_type_pairs = (
        pl.DataFrame({"cell_state": obs_states, "cell_type": obs_types})
        .group_by(["cell_state", "cell_type"])
        .len()
        .sort(["cell_type", "cell_state"])
    )
    hierarchy_table = state_type_pairs.with_columns(
        pl.col("cell_state")
        .map_elements(lambda s: state_to_malignant.get(str(s), False), return_dtype=pl.Boolean)
        .alias("is_malignant")
    ).rename({"len": "n_cells"})

    # 9. Compute library sizes & relative mRNA scaling factors
    counts_matrix = working_adata.X
    lib_sizes = np.asarray(counts_matrix.sum(axis=1)).flatten()
    mrna_scaling = _compute_mrna_scaling(lib_sizes, obs_states, obs_types, state_to_malignant)

    # 10. Compute Mean GEP Emission Matrices (Phi_state and Phi_type)
    gene_names = [str(g) for g in working_adata.var_names]
    phi_state, state_profiles = _compute_gep_matrix(
        X=counts_matrix,
        group_labels=obs_states,
        unique_groups=unique_states,
        gene_names=gene_names,
        id_col="cell_state",
        normalize=config.normalize_multinomial,
        pseudo_min=config.pseudo_min,
    )

    phi_type, _ = _compute_gep_matrix(
        X=counts_matrix,
        group_labels=obs_types,
        unique_groups=unique_types,
        gene_names=gene_names,
        id_col="cell_type",
        normalize=config.normalize_multinomial,
        pseudo_min=config.pseudo_min,
    )

    # 11. Collinearity analysis
    collinearity = _evaluate_collinearity(state_profiles, unique_states, config.collinearity_threshold)
    n_collinear = int(collinearity.filter(pl.col("is_collinear")).height) if collinearity.height > 0 else 0
    if n_collinear > 0:
        logger.warning("Identified %d collinear cell state pairs with r >= %.2f", n_collinear, config.collinearity_threshold)

    elapsed = time.time() - start_time
    logger.info(
        "Successfully constructed deconvolution reference (%d states, %d types, %d genes, %d malignant states) in %.2fs",
        len(unique_states),
        len(unique_types),
        len(gene_names),
        len(malignant_states),
        elapsed,
    )

    return Success(
        DeconvolutionReferenceResult(
            phi_state=phi_state,
            phi_type=phi_type,
            hierarchy_table=hierarchy_table,
            mrna_scaling=mrna_scaling,
            collinearity=collinearity,
            gene_names=tuple(gene_names),
            cell_states=tuple(unique_states),
            cell_types=tuple(unique_types),
            malignant_states=malignant_states,
            metadata={
                "n_cells": working_adata.n_obs,
                "n_genes": working_adata.n_vars,
                "n_states": len(unique_states),
                "n_types": len(unique_types),
                "n_malignant_states": len(malignant_states),
                "n_collinear_pairs": n_collinear,
                "elapsed_seconds": round(elapsed, 3),
            },
        )
    )
