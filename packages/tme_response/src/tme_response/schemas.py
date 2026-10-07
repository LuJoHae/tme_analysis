"""Immutable data models and schemas for multi-omic cohorts, predictors, and benchmark results."""

from __future__ import annotations

from typing import Mapping, Sequence
from pydantic import BaseModel, ConfigDict
import polars as pl
from returns.maybe import Maybe, Nothing

from .types import PredictorCategory


class GenomicVariant(BaseModel):
    """Immutable representation of a genomic somatic alteration."""
    model_config = ConfigDict(frozen=True)

    sample_id: str
    gene: str
    variant_classification: str
    is_frameshift: bool
    is_inactivating: bool


class MultiOmicCohort(BaseModel):
    """Synchronized multi-omic cohort containing bulk transcriptomics and genomics."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    cohort_id: str
    cancer_type: str
    sample_ids: tuple[str, ...]
    expression_tpm: pl.DataFrame          # [sample_id, gene_1, gene_2, ...]
    clinical_annotations: pl.DataFrame    # [sample_id, response_binary, response_recist, biopsy_timepoint]
    tmb_scores: Maybe[pl.DataFrame]       # [sample_id, tmb_per_mb, n_nonsynonymous]
    driver_mutations: Maybe[pl.DataFrame] # [sample_id, gene, variant_classification, is_inactivating]
    cna_scores: Maybe[pl.DataFrame]       # [sample_id, cdkn2a_loss, aneuploidy_score]


class PredictionResult(BaseModel):
    """Immutable output from a response predictor."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    predictor_name: str
    category: PredictorCategory
    predictions: pl.DataFrame  # [sample_id, score, Maybe[probability], Maybe[binary_call]]


class CohortBenchmarkResult(BaseModel):
    """Consolidated performance evaluation for a predictor on a specific cohort."""
    model_config = ConfigDict(frozen=True)

    cohort_id: str
    cancer_type: str
    predictor_name: str
    category: PredictorCategory
    n_samples: int
    n_responders: int
    n_non_responders: int
    time_stratum: str = "Pre"
    response_stratum: str = "Standard"
    roc_auc: float
    roc_auc_ci_lower: float = 0.5
    roc_auc_ci_upper: float = 0.5
    p_value_vs_half: float = 1.0
    pr_auc: float
    pr_auc_ci_lower: float = 0.0
    pr_auc_ci_upper: float = 1.0
    baseline_prevalence: float = 0.5
    delta_pr_auc: float = 0.0
    p_value_prauc: float = 1.0
    pooling_strategy: str = "cohort"
    brier_score: float


class CohortPredictabilityResult(BaseModel):
    """Consolidated intrinsic predictability metrics for a cohort across predictors."""
    model_config = ConfigDict(frozen=True)

    cohort_id: str
    cancer_type: str
    time_stratum: str = "Pre"
    response_stratum: str = "Standard"
    pooling_strategy: str = "cohort"
    n_samples: int
    n_responders: int
    baseline_prevalence: float
    mean_roc_auc_rna: float
    median_roc_auc_rna: float
    std_roc_auc_rna: float
    mean_pr_auc_rna: float
    mean_delta_pr_auc_rna: float
    best_predictor_rna: str
    max_roc_auc_rna: float
    mean_roc_auc_all: float
    best_predictor_all: str
    max_roc_auc_all: float
    n_predictors_evaluated: int


class SynergyConfig(BaseModel):
    """Configuration for composite multi-omic DNA-RNA models."""
    model_config = ConfigDict(frozen=True)

    rna_weight: float = 1.0
    tmb_weight: float = 0.5
    enable_antigen_presentation_gating: bool = True
    gating_genes: tuple[str, ...] = ("B2M", "JAK1", "JAK2")
