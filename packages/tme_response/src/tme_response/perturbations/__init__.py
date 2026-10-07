"""Perturbation models and robustness benchmarking suite for tme_response."""

from __future__ import annotations

from .engine import (
    compute_standard_cohort_predictors,
    run_cohort_perturbation_sweep,
)
from .expression import (
    apply_expression_jitter,
    apply_gene_dropout,
    apply_immune_dilution,
)
from .labels import apply_label_noise
from .models import (
    DilutionConfig,
    DropoutConfig,
    JitterConfig,
    LabelNoiseConfig,
    PerturbationBenchmarkRecord,
    PerturbationResilienceRecord,
    PerturbationSweepConfig,
)
from .resilience import compute_perturbation_resilience

__all__ = [
    # Configs & Schemas
    "JitterConfig",
    "DropoutConfig",
    "DilutionConfig",
    "LabelNoiseConfig",
    "PerturbationSweepConfig",
    "PerturbationBenchmarkRecord",
    "PerturbationResilienceRecord",
    # Operators
    "apply_expression_jitter",
    "apply_gene_dropout",
    "apply_immune_dilution",
    "apply_label_noise",
    # Engine & Metrics
    "compute_perturbation_resilience",
    "compute_standard_cohort_predictors",
    "run_cohort_perturbation_sweep",
]
