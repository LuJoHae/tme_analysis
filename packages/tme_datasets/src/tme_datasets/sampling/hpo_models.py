"""Immutable data models and configuration for single-cell sampling HPO."""

from __future__ import annotations

from enum import Enum
from typing import Any, Mapping, Sequence
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Nothing, Some

from ..deconvolution.models import DeconvolutionReferenceConfig
from ..models import (
    ClusterAnalysisSpec,
    HarmonizeConfig,
    SingleCellSamplingSpec,
)
from ..types import CohortSamplingMode, GeneIDType, HarmonizeMode
from .classifier_hpo import InnerClassifierConfig, InnerOptimizationResult
from .malignant_sampling import MalignantSamplingConfig, MalignantStrategy


class FidelityRung(str, Enum):
    """Multi-fidelity evaluation rung for ASHA/Hyperband scheduling."""

    RUNG_0_SCREENING = "rung_0_screening"
    RUNG_1_REFINEMENT = "rung_1_refinement"
    RUNG_2_FULL = "rung_2_full"


class HPOObjectiveMetric(str, Enum):
    """Primary objective metric for downstream deconvolution evaluation."""

    LOCO_ROC_AUC = "loco_roc_auc"
    LOCO_PR_AUC = "loco_pr_auc"
    COMBINED_ROC_PR = "combined_roc_pr"


class HPOSearchSpace(BaseModel):
    """Search space specification for single-cell reference sampling and deconvolution."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    candidate_cohort_ids: tuple[str, ...]
    min_cohorts: int = 2
    max_cohorts: int = 5
    sampling_modes: tuple[CohortSamplingMode, ...] = (
        CohortSamplingMode.FIXED_PER_COHORT,
        CohortSamplingMode.GLOBAL_BUDGET,
        CohortSamplingMode.FRACTION_PER_COHORT,
    )
    min_cells_per_cohort: int = 150
    max_cells_per_cohort: int = 2500
    min_global_budget: int = 600
    max_global_budget: int = 6000
    stratify_options: tuple[str | None, ...] = (None, "cell_type", "patient_id")
    malignant_strategies: tuple[MalignantStrategy, ...] = (
        MalignantStrategy.PATIENT_STRATIFIED,
        MalignantStrategy.POOLED_GENERIC,
        MalignantStrategy.EXCLUDED_TME_ONLY,
    )
    malignant_fraction_range: tuple[float, float] = (0.05, 0.40)
    harmonize_modes: tuple[HarmonizeMode, ...] = (
        HarmonizeMode.INTERSECTION,
        HarmonizeMode.UNION_ZERO_FILLED,
    )
    gene_target_types: tuple[GeneIDType, ...] = (
        GeneIDType.ENSEMBL_ID,
        GeneIDType.HUGO_SYMBOL,
    )
    min_shared_genes_range: tuple[int, int] = (1000, 5000)
    leiden_resolution_range: tuple[float, float] = (0.3, 1.4)
    n_top_genes_range: tuple[int, int] = (1500, 4000)
    n_pcs_range: tuple[int, int] = (20, 60)
    collinearity_threshold_range: tuple[float, float] = (0.75, 0.92)
    min_cells_per_state_range: tuple[int, int] = (10, 40)
    require_cached_h5ad: bool = True


class HPOTrialConfig(BaseModel):
    """Immutable parameter configuration for an individual optimization trial."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    trial_id: int
    selected_cohort_ids: tuple[str, ...]
    sampling_spec: SingleCellSamplingSpec
    malignant_config: MalignantSamplingConfig
    cluster_spec: ClusterAnalysisSpec
    harmonize_config: HarmonizeConfig
    reference_config: DeconvolutionReferenceConfig


class TrialEvaluationResult(BaseModel):
    """Immutable evaluation metrics for a single optimization trial."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    trial_id: int
    rung: FidelityRung
    mean_loco_auc: float
    mean_loco_pr_auc: float
    cohort_aucs: Mapping[str, float]
    collinearity_max: float
    condition_number: float
    n_cell_states: int
    n_shared_genes: int
    total_cells_sampled: int
    elapsed_seconds: float
    inner_optimization: Maybe[InnerOptimizationResult] = Nothing
    is_pruned: bool = False
    prune_reason: Maybe[str] = Nothing


class HPORunResult(BaseModel):
    """Immutable summary and Pareto-optimal outputs of an entire HPO campaign."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    best_trial_id: int
    best_config: HPOTrialConfig
    best_loco_auc: float
    evaluations: tuple[TrialEvaluationResult, ...]
    evaluations_df: pl.DataFrame
    pareto_trials: tuple[int, ...]
    best_classifier_config: Maybe[InnerClassifierConfig] = Nothing
    total_trials: int
    pruned_trials: int
    total_elapsed_seconds: float
    metadata: Mapping[str, Any] = Field(default_factory=dict)
