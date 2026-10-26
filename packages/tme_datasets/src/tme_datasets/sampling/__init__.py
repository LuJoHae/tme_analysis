"""Multi-cohort sampling and graph clustering subpackage."""

from __future__ import annotations

from .classifier_hpo import (
    ClassifierType,
    FeatureSelectorType,
    InnerClassifierConfig,
    InnerOptimizationResult,
    TransformType,
    apply_compositional_transform,
    compute_single_cohort_oof_predictions,
    evaluate_inner_pipeline,
    evaluate_single_cohort_stratified_cv,
    run_inner_classifier_hpo,
    run_single_cohort_classifier_hpo,
)
from .single_cohort_hpo import (
    evaluate_single_cohort_trial_with_fidelity,
    run_single_cohort_sampling_hpo,
)
from .cohort_sampler import (
    compute_cohort_cell_allocations,
    generate_random_cluster_spec,
    generate_random_sampling_spec,
    resolve_sampled_cohorts,
    run_pca_knn_leiden,
    sample_and_harmonize_cohorts,
    sample_single_cell_cohorts,
)
from .hpo_models import (
    FidelityRung,
    HPOObjectiveMetric,
    HPORunResult,
    HPOSearchSpace,
    HPOTrialConfig,
    TrialEvaluationResult,
)
from .hpo_optimizer import (
    evaluate_trial_with_fidelity,
    fast_deconvolute_cohorts,
    fast_loco_cv_evaluation,
    filter_candidate_cohorts,
    generate_trial_config,
    run_sampling_hpo,
)
from .malignant_sampling import (
    MalignantSamplingConfig,
    MalignantStrategy,
    sample_malignant_and_tme_cells,
)

__all__ = [
    "compute_cohort_cell_allocations",
    "generate_random_cluster_spec",
    "generate_random_sampling_spec",
    "resolve_sampled_cohorts",
    "run_pca_knn_leiden",
    "sample_and_harmonize_cohorts",
    "sample_single_cell_cohorts",
    "FidelityRung",
    "HPOObjectiveMetric",
    "HPORunResult",
    "HPOSearchSpace",
    "HPOTrialConfig",
    "TrialEvaluationResult",
    "evaluate_trial_with_fidelity",
    "fast_deconvolute_cohorts",
    "fast_loco_cv_evaluation",
    "filter_candidate_cohorts",
    "generate_trial_config",
    "run_sampling_hpo",
    "run_single_cohort_sampling_hpo",
    "evaluate_single_cohort_trial_with_fidelity",
    "evaluate_single_cohort_stratified_cv",
    "compute_single_cohort_oof_predictions",
    "run_single_cohort_classifier_hpo",
    "ClassifierType",
    "FeatureSelectorType",
    "InnerClassifierConfig",
    "InnerOptimizationResult",
    "TransformType",
    "apply_compositional_transform",
    "evaluate_inner_pipeline",
    "run_inner_classifier_hpo",
    "MalignantSamplingConfig",
    "MalignantStrategy",
    "sample_malignant_and_tme_cells",
]
