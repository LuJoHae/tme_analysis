"""Stability Selection module with finite-sample error control.

Implements Meinshausen & Bühlmann (2010) and Shah & Samworth (2013)
Complementary Pairs Stability Selection (CPSS).
"""

from selective_inference.stability_selection.types import (
    PathFitter,
    StabilityParameters,
    StabilityResult,
    SamplingType,
    Assumption,
)
from selective_inference.stability_selection.bounds import (
    compute_unimodal_constant,
    solve_pfer,
    solve_cutoff,
    solve_q,
    resolve_stability_parameters,
)
from selective_inference.stability_selection.subsampling import (
    generate_complementary_pairs,
    generate_subsamples,
    generate_stratified_complementary_pairs,
    generate_stratified_subsamples,
    apply_randomized_weights,
)
from selective_inference.stability_selection.core import (
    run_stability_selection,
    default_lasso_path_fitter,
    generate_default_lambdas,
)
from selective_inference.stability_selection.fitters import (
    create_lasso_fitter,
    create_elastic_net_fitter,
    create_l1_logistic_fitter,
    create_tree_importance_fitter,
    create_cohort_adjusted_fitter,
    create_group_lasso_cohort_fitter,
    create_merf_cohort_fitter,
    create_multitask_logistic_cohort_fitter,
    create_meta_analysis_cohort_fitter,
    create_multistudy_invariant_cohort_fitter,
    create_glmm_lasso_cohort_fitter,
    create_oscar_fitter,
    create_slope_fitter,
    create_oscar_cohort_fitter,
    create_slope_cohort_fitter,
    build_fitter,
)
from selective_inference.stability_selection.sklearn_adapter import (
    StabilitySelector,
)
from selective_inference.stability_selection.visualization import (
    plot_stability_paths,
    plot_stability_scores,
)

__all__ = [
    "PathFitter",
    "StabilityParameters",
    "StabilityResult",
    "SamplingType",
    "Assumption",
    "compute_unimodal_constant",
    "solve_pfer",
    "solve_cutoff",
    "solve_q",
    "resolve_stability_parameters",
    "generate_complementary_pairs",
    "generate_subsamples",
    "generate_stratified_complementary_pairs",
    "generate_stratified_subsamples",
    "apply_randomized_weights",
    "run_stability_selection",
    "default_lasso_path_fitter",
    "generate_default_lambdas",
    "create_lasso_fitter",
    "create_elastic_net_fitter",
    "create_l1_logistic_fitter",
    "create_tree_importance_fitter",
    "create_cohort_adjusted_fitter",
    "create_group_lasso_cohort_fitter",
    "create_merf_cohort_fitter",
    "create_multitask_logistic_cohort_fitter",
    "create_meta_analysis_cohort_fitter",
    "create_multistudy_invariant_cohort_fitter",
    "create_glmm_lasso_cohort_fitter",
    "create_oscar_fitter",
    "create_slope_fitter",
    "create_oscar_cohort_fitter",
    "create_slope_cohort_fitter",
    "build_fitter",
    "StabilitySelector",
    "plot_stability_paths",
    "plot_stability_scores",
]
