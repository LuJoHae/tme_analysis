"""Shared domain models, data ingestion, feature filtering, and execution pipeline for stability selection.

Follows strict functional programming principles:
- Pure functions for matrix transformation, feature selection, and table assembly
- Monadic error handling with Result[T, str] and Maybe[T]
- Dataset loading strictly via tme_datasets API
- Columnar Parquet serialization with Polars
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path
from typing import Any, Final, Mapping, Sequence
import anndata as ad  # type: ignore[import-untyped]
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

from selective_inference.stability_selection.bounds import resolve_stability_parameters
from selective_inference.stability_selection.core import run_stability_selection
from selective_inference.stability_selection.fitters import build_fitter
from selective_inference.stability_selection.types import PathFitter, StabilityResult
from tme_datasets import (  # type: ignore[import-untyped]
    IMMUNE_CHECKPOINT_GENES,
    IMMUNOTHERAPY_GENE_PANEL,
    load_iatlas_cohort_or_combined,
)

DEFAULT_COHORTS: Final[tuple[str, ...]] = (
    "Hugo-iAtlas",
    "Riaz-iAtlas",
    "Liu-iAtlas",
    "Gide-iAtlas",
    "Rosenberg-iAtlas",
    "Padron-iAtlas",
    "Anders-iAtlas",
    "McDermott-iAtlas",
    "Choueiri-iAtlas",
    "melanoma",
    "rcc",
    "pancancer",
)

COMBINED_COHORTS: Final[tuple[str, ...]] = ("pancancer", "melanoma", "rcc")

FITTERS_FOR_COMBINED: Final[tuple[tuple[str, bool], ...]] = (
    ("lasso", False),
    ("elastic_net", False),
    ("logistic", False),
    ("rf", False),
    ("cohort_adjusted", True),
    ("group_lasso", True),
    ("merf", True),
    ("multitask_logistic", True),
    ("meta_analysis", True),
    ("multistudy_invariant", True),
    ("glmm_lasso", True),
    ("oscar", True),
    ("slope", True),
)

FITTERS_FOR_SINGLE: Final[tuple[tuple[str, bool], ...]] = (
    ("lasso", False),
    ("elastic_net", False),
    ("logistic", False),
    ("rf", False),
    ("oscar", False),
    ("slope", False),
)


class CohortSelectionOutput(BaseModel):
    """Immutable bundle of stability selection results for a single cohort."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    cohort: str
    n_samples: int
    n_responders: int
    n_non_responders: int
    n_features: int
    mb_result: StabilityResult
    ss_result: StabilityResult
    is_checkpoint_map: Mapping[str, bool]
    fitter: str = "lasso"
    is_stratified: bool = False


def load_cohort_adata(cohort_id: str) -> Result[ad.AnnData, str]:
    """Declarative loader routing strictly through tme_datasets API."""
    return load_iatlas_cohort_or_combined(cohort_id)


def resolve_gene_names(adata: ad.AnnData) -> tuple[str, ...]:
    """Extract clean, unique human-readable gene symbols from AnnData."""
    if "gene_name" in adata.var.columns:
        raw_names = adata.var["gene_name"].fillna("").astype(str).tolist()
        resolved = [name if name.strip() else str(idx) for name, idx in zip(raw_names, adata.var_names)]
    else:
        resolved = [str(idx) for idx in adata.var_names]
    return tuple(resolved)


def extract_clean_matrix(
    adata: ad.AnnData,
) -> Result[tuple[np.ndarray, np.ndarray, tuple[str, ...], np.ndarray | None], str]:
    """Extract dense expression matrix, binary response labels, and optional cohort labels."""
    if "response_binary" not in adata.obs.columns:
        return Failure("Dataset missing 'response_binary' column in obs metadata")

    response_series = adata.obs["response_binary"]
    valid_mask = response_series.notnull().values & (
        (response_series.values == 0.0) | (response_series.values == 1.0)
    )

    n_valid = int(np.sum(valid_mask))
    if n_valid < 10:
        return Failure(f"Insufficient valid binary response samples ({n_valid} < 10)")

    adata_clean = adata[valid_mask]
    y = adata_clean.obs["response_binary"].values.astype(np.float64)

    n_resp = int(np.sum(y == 1.0))
    n_non_resp = int(np.sum(y == 0.0))
    if n_resp < 2 or n_non_resp < 2:
        return Failure(
            f"Class imbalance too severe for stability selection: responders={n_resp}, non-responders={n_non_resp}"
        )

    raw_X: Any = adata_clean.X
    X_dense = raw_X.toarray() if hasattr(raw_X, "toarray") else np.asarray(raw_X)
    X = np.asarray(X_dense, dtype=np.float64)

    cohort_labels: np.ndarray | None = None
    if "dataset_id" in adata_clean.obs.columns:
        cohort_labels = np.asarray(adata_clean.obs["dataset_id"].values)
    elif "cohort" in adata_clean.obs.columns:
        cohort_labels = np.asarray(adata_clean.obs["cohort"].values)

    gene_names = resolve_gene_names(adata_clean)
    return Success((X, y, gene_names, cohort_labels))


def select_feature_matrix(
    X: np.ndarray,
    gene_names: tuple[str, ...],
    n_top_genes: int,
    checkpoint_genes: Sequence[str] = IMMUNE_CHECKPOINT_GENES,
) -> Result[tuple[np.ndarray, tuple[str, ...], tuple[bool, ...]], str]:
    """Filter to top highly variable genes unioned with immune checkpoint biomarkers."""
    n_samples, n_features = X.shape
    if n_features == 0 or n_samples == 0:
        return Failure("Empty feature matrix provided")

    variances = np.var(X, axis=0)
    finite_mask = np.isfinite(variances)
    variances_clean = np.where(finite_mask, variances, 0.0)

    sorted_indices = np.argsort(-variances_clean)
    target_hvg_count = min(n_top_genes, n_features)
    hvg_indices = set(sorted_indices[:target_hvg_count].tolist())

    gene_to_idx = {name: idx for idx, name in enumerate(gene_names)}
    checkpoint_indices = {gene_to_idx[g] for g in checkpoint_genes if g in gene_to_idx}

    combined_set = hvg_indices | checkpoint_indices
    selected_indices = sorted(combined_set)

    if not selected_indices:
        return Failure("No features selected after HVG and checkpoint filtering")

    X_sub = X[:, selected_indices]
    means = np.mean(X_sub, axis=0)
    stds = np.std(X_sub, axis=0)
    stds_safe = np.where(stds > 1e-8, stds, 1.0)
    X_scaled = (X_sub - means) / stds_safe

    selected_names = tuple(gene_names[i] for i in selected_indices)
    is_checkpoint = tuple(i in checkpoint_indices for i in selected_indices)

    return Success((X_scaled, selected_names, is_checkpoint))


def select_immunotherapy_feature_matrix(
    X: np.ndarray,
    gene_names: tuple[str, ...],
    gene_panel: Sequence[str] = IMMUNOTHERAPY_GENE_PANEL,
) -> Result[tuple[np.ndarray, tuple[str, ...], tuple[bool, ...]], str]:
    """Filter strictly to preselected immunotherapy genes and standardize columns."""
    n_samples, n_features = X.shape
    if n_features == 0 or n_samples == 0:
        return Failure("Empty feature matrix provided")

    gene_to_idx = {name: idx for idx, name in enumerate(gene_names)}
    matched_indices = [gene_to_idx[g] for g in gene_panel if g in gene_to_idx]
    if not matched_indices:
        return Failure("None of the preselected immunotherapy genes were found in the dataset")

    matched_indices_sorted = sorted(matched_indices)
    X_sub = X[:, matched_indices_sorted]

    means = np.mean(X_sub, axis=0)
    stds = np.std(X_sub, axis=0)
    stds_safe = np.where(stds > 1e-8, stds, 1.0)
    X_scaled = (X_sub - means) / stds_safe

    selected_names = tuple(gene_names[i] for i in matched_indices_sorted)
    is_checkpoint = tuple(gene_names[i] in IMMUNE_CHECKPOINT_GENES for i in matched_indices_sorted)

    return Success((X_scaled, selected_names, is_checkpoint))


def run_cohort_analysis(
    cohort_name: str,
    adata: ad.AnnData,
    feature_mode: str = "immunotherapy",
    n_top_genes: int = 500,
    pfer: float = 1.0,
    cutoff: float = 0.75,
    B: int = 50,
    seed: int = 42,
    max_expected_q: float | None = None,
    fitter_name: str = "lasso",
    l1_ratio: float = 0.7,
    stratified: bool = False,
    kappa: float = 0.5,
    q_fdr: float = 0.1,
) -> Result[CohortSelectionOutput, str]:
    """Execute both MB and SS-CPSS stability selection on a single cohort."""
    match extract_clean_matrix(adata):
        case Failure(err):
            return Failure(f"Data extraction failed for '{cohort_name}': {err}")
        case Success((X, y, gene_names, cohort_labels)):
            pass

    if feature_mode == "immunotherapy":
        feat_res = select_immunotherapy_feature_matrix(X, gene_names)
    else:
        feat_res = select_feature_matrix(X, gene_names, n_top_genes=n_top_genes)

    match feat_res:
        case Failure(err):
            return Failure(f"Feature selection failed for '{cohort_name}': {err}")
        case Success((X_scaled, features, is_checkpoint_tuple)):
            pass

    p = len(features)
    is_checkpoint_map = {feat: is_chk for feat, is_chk in zip(features, is_checkpoint_tuple)}

    fitter = build_fitter(
        fitter_name=fitter_name,
        cohort_labels=cohort_labels,
        l1_ratio=l1_ratio,
        kappa=kappa,
        q_fdr=q_fdr,
    )

    strata_for_run = cohort_labels if (stratified and cohort_labels is not None and len(np.unique(cohort_labels)) > 1) else None

    # Resolve parameters for MB (2010)
    mb_params_res = resolve_stability_parameters(
        p=p,
        cutoff=Some(cutoff),
        pfer=Some(pfer),
        B=B,
        sampling_type="MB",
        assumption="none",
    )
    match mb_params_res:
        case Failure(err):
            return Failure(f"Failed to resolve MB parameters for '{cohort_name}': {err}")
        case Success(mb_params):
            pass

    # Resolve parameters for SS-CPSS (2013)
    ss_params_res = resolve_stability_parameters(
        p=p,
        cutoff=Some(cutoff),
        pfer=Some(pfer),
        B=B,
        sampling_type="SS",
        assumption="unimodal",
    )
    match ss_params_res:
        case Failure(err):
            return Failure(f"Failed to resolve SS parameters for '{cohort_name}': {err}")
        case Success(ss_params):
            pass

    # Suppress coordinate descent convergence warnings during subsampling runs
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Run MB
        mb_run_res = run_stability_selection(
            X=X_scaled,
            y=y,
            parameters=mb_params,
            fitter=fitter,
            feature_names=features,
            seed=seed,
            max_expected_q=max_expected_q,
            strata=strata_for_run,
        )
        match mb_run_res:
            case Failure(err):
                return Failure(f"MB stability selection failed for '{cohort_name}': {err}")
            case Success(mb_result):
                pass

        # Run SS-CPSS
        ss_run_res = run_stability_selection(
            X=X_scaled,
            y=y,
            parameters=ss_params,
            fitter=fitter,
            feature_names=features,
            seed=seed,
            max_expected_q=max_expected_q,
            strata=strata_for_run,
        )
        match ss_run_res:
            case Failure(err):
                return Failure(f"SS stability selection failed for '{cohort_name}': {err}")
            case Success(ss_result):
                pass

    n_samples = int(X_scaled.shape[0])
    n_resp = int(np.sum(y == 1.0))
    n_non_resp = int(np.sum(y == 0.0))

    return Success(
        CohortSelectionOutput(
            cohort=cohort_name,
            n_samples=n_samples,
            n_responders=n_resp,
            n_non_responders=n_non_resp,
            n_features=p,
            mb_result=mb_result,
            ss_result=ss_result,
            is_checkpoint_map=is_checkpoint_map,
            fitter=fitter_name,
            is_stratified=bool(strata_for_run is not None),
        )
    )


def assemble_parquet_tables(
    outputs: Sequence[CohortSelectionOutput],
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Convert a sequence of cohort outputs into structured Polars DataFrames."""
    path_records = []
    score_records = []
    summary_records = []

    for out in outputs:
        mb_cutoff_val = out.mb_result.lambda_cutoff.value_or(0.0)
        ss_cutoff_val = out.ss_result.lambda_cutoff.value_or(0.0)

        summary_records.append({
            "cohort": out.cohort,
            "n_samples": out.n_samples,
            "n_responders": out.n_responders,
            "n_non_responders": out.n_non_responders,
            "n_features": out.n_features,
            "fitter": out.fitter,
            "is_stratified": out.is_stratified,
            "mb_q": out.mb_result.parameters.q,
            "ss_q": out.ss_result.parameters.q,
            "actual_mb_q": out.mb_result.empirical_q,
            "actual_ss_q": out.ss_result.empirical_q,
            "mb_selected_count": len(out.mb_result.selected_features),
            "ss_selected_count": len(out.ss_result.selected_features),
            "lambda_cutoff_mb": mb_cutoff_val,
            "lambda_cutoff_ss": ss_cutoff_val,
        })

        for method_name, res in (("MB", out.mb_result), ("SS-CPSS", out.ss_result)):
            cutoff_val = res.parameters.cutoff
            pfer_val = res.parameters.pfer
            q_val = res.parameters.q
            actual_q_val = res.empirical_q
            lam_cutoff_val = res.lambda_cutoff.value_or(0.0)

            n_feats = len(res.feature_names)
            for i in range(n_feats):
                feat = res.feature_names[i]
                score = float(res.stability_scores[i])
                unres_score = (
                    float(res.unrestricted_stability_scores[i])
                    if i < len(res.unrestricted_stability_scores)
                    else score
                )
                is_sel = bool(feat in res.selected_features)
                is_chk = bool(out.is_checkpoint_map.get(feat, False))
                score_records.append({
                    "cohort": out.cohort,
                    "method": method_name,
                    "fitter": out.fitter,
                    "is_stratified": out.is_stratified,
                    "feature": feat,
                    "stability_score": score,
                    "unrestricted_score": unres_score,
                    "selected": is_sel,
                    "cutoff": float(cutoff_val),
                    "pfer": float(pfer_val),
                    "q": float(q_val),
                    "actual_q": float(actual_q_val),
                    "lambda_cutoff": float(lam_cutoff_val),
                    "is_checkpoint": is_chk,
                })

            n_lambdas = len(res.lambdas)
            for i in range(n_feats):
                feat = res.feature_names[i]
                is_sel = bool(feat in res.selected_features)
                is_chk = bool(out.is_checkpoint_map.get(feat, False))
                for j in range(n_lambdas):
                    lam = float(res.lambdas[j])
                    prob = float(res.stability_matrix[i][j])
                    log_lam = float(np.log10(lam)) if lam > 0 else 0.0
                    exp_size = (
                        float(res.expected_model_sizes[j])
                        if j < len(res.expected_model_sizes)
                        else 0.0
                    )
                    in_budget = bool(lam >= lam_cutoff_val - 1e-9)
                    path_records.append({
                        "cohort": out.cohort,
                        "method": method_name,
                        "fitter": out.fitter,
                        "is_stratified": out.is_stratified,
                        "feature": feat,
                        "lambda": lam,
                        "log_lambda": log_lam,
                        "selection_probability": prob,
                        "expected_model_size": exp_size,
                        "in_budget": in_budget,
                        "lambda_cutoff": float(lam_cutoff_val),
                        "selected": is_sel,
                        "is_checkpoint": is_chk,
                        "cutoff": float(cutoff_val),
                    })

    paths_df = pl.DataFrame(path_records)
    scores_df = pl.DataFrame(score_records)
    summary_df = pl.DataFrame(summary_records)

    return paths_df, scores_df, summary_df
