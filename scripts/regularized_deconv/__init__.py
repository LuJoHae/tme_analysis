"""
Synthetic Deconvolution Benchmarking Core Package.
"""

from __future__ import annotations

from .synthetic_data import (
    DEFAULT_HETERO_ALPHA,
    SyntheticSignature,
    generate_synthetic_reference,
    sample_true_proportions,
    generate_bulk_mixtures,
)
from .deconv_runners import (
    ALL_METHODS,
    run_single_method,
    run_deconvolution_suite,
)
from .metrics import (
    BiasVarianceMetrics,
    compute_bias_variance_decomposition,
    compute_population_metrics,
    compute_phantom_detection,
)

__all__ = [
    "DEFAULT_HETERO_ALPHA",
    "SyntheticSignature",
    "generate_synthetic_reference",
    "sample_true_proportions",
    "generate_bulk_mixtures",
    "ALL_METHODS",
    "run_single_method",
    "run_deconvolution_suite",
    "BiasVarianceMetrics",
    "compute_bias_variance_decomposition",
    "compute_population_metrics",
    "compute_phantom_detection",
]
