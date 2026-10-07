"""Immutable data models and configuration schemas for perturbation and robustness benchmarking."""

from __future__ import annotations

from typing import Sequence
from pydantic import BaseModel, ConfigDict


class JitterConfig(BaseModel):
    """Configuration for multiplicative log-normal expression jitter."""
    model_config = ConfigDict(frozen=True)

    sigma: float
    seed: int = 42


class DropoutConfig(BaseModel):
    """Configuration for technical gene dropout / unmeasured gene zero-masking."""
    model_config = ConfigDict(frozen=True)

    dropout_rate: float
    seed: int = 42


class DilutionConfig(BaseModel):
    """Configuration for non-immune stromal infiltration and tumor purity dilution."""
    model_config = ConfigDict(frozen=True)

    dilution_factor: float
    seed: int = 42


class LabelNoiseConfig(BaseModel):
    """Configuration for clinical response classification label noise."""
    model_config = ConfigDict(frozen=True)

    noise_rate: float
    seed: int = 42


class PerturbationSweepConfig(BaseModel):
    """Parametric sweep grid across 4 orthogonal perturbation modalities."""
    model_config = ConfigDict(frozen=True)

    jitter_sigmas: tuple[float, ...] = (0.0, 0.25, 0.50, 1.00, 1.50)
    dropout_rates: tuple[float, ...] = (0.0, 0.05, 0.10, 0.25, 0.50)
    dilution_factors: tuple[float, ...] = (1.0, 0.8, 0.6, 0.4, 0.2)
    label_noise_rates: tuple[float, ...] = (0.0, 0.05, 0.10, 0.20, 0.30)


class PerturbationBenchmarkRecord(BaseModel):
    """Point evaluation for a predictor under a specific perturbation intensity."""
    model_config = ConfigDict(frozen=True)

    cohort_id: str
    cancer_type: str
    predictor_name: str
    category: str
    perturbation_type: str  # 'jitter', 'dropout', 'dilution', 'label_noise'
    intensity: float
    roc_auc: float
    roc_auc_ci_lower: float
    roc_auc_ci_upper: float
    pr_auc: float
    pr_auc_ci_lower: float
    pr_auc_ci_upper: float
    delta_pr_auc: float
    baseline_prevalence: float
    n_samples: int


class PerturbationResilienceRecord(BaseModel):
    """Summary robustness metric (Perturbation Resilience Index, PRI) across an intensity curve."""
    model_config = ConfigDict(frozen=True)

    cohort_id: str
    predictor_name: str
    category: str
    perturbation_type: str
    baseline_roc_auc: float
    min_roc_auc: float
    pri_score: float  # Normalized Area Under Retention Curve (1.0 = perfect robustness)
    relative_auc_retained: float  # AUC(max_intensity) / AUC(baseline)
