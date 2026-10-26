"""Runtime bridge re-exporting symbols from `3b_pan_cohort_unified_deconv_benchmark.py`.

Python syntax forbids `from scripts.manuscript_pipeline.3b_... import ...` because
module identifiers cannot start with digits. This module loads stage `3b` via
`importlib.import_module` and is paired with `stage3b_bridge.pyi` so `mypy --strict`
validates all downstream callers with exact static types.
"""

from __future__ import annotations

import importlib

_mod = importlib.import_module(
    "scripts.manuscript_pipeline.3b_pan_cohort_unified_deconv_benchmark"
)

COHORT_ORDER = _mod.COHORT_ORDER
CohortEmpiricalBundle = _mod.CohortEmpiricalBundle
UnifiedBenchmarkConfig = _mod.UnifiedBenchmarkConfig
UnifiedReferenceAtlas = _mod.UnifiedReferenceAtlas
build_unified_reference_atlas = _mod.build_unified_reference_atlas
compute_permutation_pvalue = _mod.compute_permutation_pvalue
compute_within_cohort_specificity = _mod.compute_within_cohort_specificity
load_all_nine_cohorts = _mod.load_all_nine_cohorts

__all__ = [
    "COHORT_ORDER",
    "CohortEmpiricalBundle",
    "UnifiedBenchmarkConfig",
    "UnifiedReferenceAtlas",
    "build_unified_reference_atlas",
    "compute_permutation_pvalue",
    "compute_within_cohort_specificity",
    "load_all_nine_cohorts",
]
