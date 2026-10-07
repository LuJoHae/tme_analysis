"""tme_response: Multi-omic immune checkpoint therapy response prediction and systematic benchmarking framework."""

from __future__ import annotations

__version__ = "0.1.0"

# Core Models & Types
from .schemas import (
    CohortBenchmarkResult,
    CohortPredictabilityResult,
    GenomicVariant,
    MultiOmicCohort,
    PredictionResult,
    SynergyConfig,
)
from .types import (
    BiopsyTimepoint,
    ClinicalResponse,
    PredictorCategory,
    ValidationMetric,
)

# Data Ingestion
from .data import (
    filter_baseline_samples,
    load_multi_omic_cohort,
    parse_cna_file,
    parse_maf_mutations,
)

# Signatures
from .signatures import (
    compute_ayers_ifng6_score,
    compute_cd8_duo_score,
    compute_cristescu_gep_score,
    compute_cyt_score,
    compute_davoli_cis_score,
    compute_fehrenbacher_teff_score,
    compute_freeman_pgm_score,
    compute_genebio_target_score,
    compute_gep_score,
    compute_huang_nrs_score,
    compute_impres_score,
    compute_ipres_score,
    compute_jiang_ctls_score,
    compute_jiang_tams_score,
    compute_jiang_texh_score,
    compute_kong_netbio_score,
    compute_messina_cks_score,
    compute_nurmik_cafs_score,
    compute_roh_is_score,
    compute_single_gene_score,
    compute_standard_signatures,
    compute_wu_mias_score,
)

# Models
from .models import (
    compute_tide_score,
    export_cohort_for_easier,
    load_easier_predictions,
)

# Multi-Omic Synergy
from .synergy import (
    apply_antigen_presentation_gating,
    compute_dna_rna_composite,
)

# Deconvolution
from .deconvolution import (
    benchmark_pseudobulk_deconvolution,
    compute_cd8_infiltrate_predictor,
    compute_mcp_lineage_scores,
)

# Evaluation
from .evaluation import (
    CIndexResult,
    CoxHazardResult,
    DCAResultRecord,
    UNIVERSAL_RNA_PREDICTORS,
    calculate_brier_score,
    calculate_c_index_bootstrap,
    calculate_decision_curve,
    calculate_discrimination_bootstrap,
    calculate_pr_auc,
    calculate_roc_auc,
    calculate_roc_auc_bootstrap,
    calibrate_probabilities,
    compute_c_index_raw,
    compute_cohort_predictability,
    compute_meta_analytic_auc,
    evaluate_prediction,
    fit_univariable_cox,
    pool_cohort_stratum_data,
    run_cohort_benchmark,
    run_multi_cohort_benchmark,
    standardize_prediction_scores,
)

# Perturbations & Robustness
from .perturbations import (
    DilutionConfig,
    DropoutConfig,
    JitterConfig,
    LabelNoiseConfig,
    PerturbationBenchmarkRecord,
    PerturbationResilienceRecord,
    PerturbationSweepConfig,
    apply_expression_jitter,
    apply_gene_dropout,
    apply_immune_dilution,
    apply_label_noise,
    compute_perturbation_resilience,
    compute_standard_cohort_predictors,
    run_cohort_perturbation_sweep,
)

# Visualization
from .visualization import (
    UncertaintyMode,
    create_auc_heatmap,
    create_benchmark_bar_chart,
    create_c_index_forest_plot,
    create_cohort_predictability_forest_plot,
    create_cox_hr_forest_plot,
    create_cross_cohort_resilience_heatmap,
    create_dca_net_benefit_chart,
    create_meta_analytic_decay_chart,
    create_pan_cohort_faceted_decay_chart,
    create_perturbation_decay_chart,
    create_pooling_comparison_chart,
    create_resilience_ranking_chart,
    create_roc_chart,
    create_summary_forest_plot,
    export_chart_svg,
)


__all__ = [
    # Types & Models
    "PredictorCategory",
    "ValidationMetric",
    "BiopsyTimepoint",
    "ClinicalResponse",
    "GenomicVariant",
    "MultiOmicCohort",
    "PredictionResult",
    "CohortBenchmarkResult",
    "CohortPredictabilityResult",
    "SynergyConfig",
    # Data
    "load_multi_omic_cohort",
    "filter_baseline_samples",
    "parse_maf_mutations",
    "parse_cna_file",
    # Signatures
    "compute_cyt_score",
    "compute_impres_score",
    "compute_gep_score",
    "compute_single_gene_score",
    "compute_ipres_score",
    "compute_standard_signatures",
    "compute_genebio_target_score",
    "compute_cd8_duo_score",
    "compute_davoli_cis_score",
    "compute_fehrenbacher_teff_score",
    "compute_freeman_pgm_score",
    "compute_huang_nrs_score",
    "compute_ayers_ifng6_score",
    "compute_jiang_ctls_score",
    "compute_jiang_tams_score",
    "compute_jiang_texh_score",
    "compute_messina_cks_score",
    "compute_nurmik_cafs_score",
    "compute_roh_is_score",
    "compute_wu_mias_score",
    "compute_cristescu_gep_score",
    "compute_kong_netbio_score",
    # Models
    "compute_tide_score",
    "export_cohort_for_easier",
    "load_easier_predictions",
    # Synergy
    "compute_dna_rna_composite",
    "apply_antigen_presentation_gating",
    # Deconvolution
    "compute_mcp_lineage_scores",
    "compute_cd8_infiltrate_predictor",
    "benchmark_pseudobulk_deconvolution",
    # Evaluation
    "calculate_roc_auc",
    "calculate_roc_auc_bootstrap",
    "calculate_discrimination_bootstrap",
    "calculate_pr_auc",
    "calculate_brier_score",
    "evaluate_prediction",
    "run_cohort_benchmark",
    "run_multi_cohort_benchmark",
    "standardize_prediction_scores",
    "pool_cohort_stratum_data",
    "compute_meta_analytic_auc",
    "compute_cohort_predictability",
    "UNIVERSAL_RNA_PREDICTORS",
    "CIndexResult",
    "CoxHazardResult",
    "compute_c_index_raw",
    "calculate_c_index_bootstrap",
    "fit_univariable_cox",
    "DCAResultRecord",
    "calibrate_probabilities",
    "calculate_decision_curve",
    # Perturbations
    "UncertaintyMode",
    "JitterConfig",
    "DropoutConfig",
    "DilutionConfig",
    "LabelNoiseConfig",
    "PerturbationSweepConfig",
    "PerturbationBenchmarkRecord",
    "PerturbationResilienceRecord",
    "apply_expression_jitter",
    "apply_gene_dropout",
    "apply_immune_dilution",
    "apply_label_noise",
    "compute_perturbation_resilience",
    "compute_standard_cohort_predictors",
    "run_cohort_perturbation_sweep",
    # Visualization
    "create_roc_chart",
    "create_benchmark_bar_chart",
    "create_auc_heatmap",
    "create_summary_forest_plot",
    "create_cohort_predictability_forest_plot",
    "create_pooling_comparison_chart",
    "create_perturbation_decay_chart",
    "create_resilience_ranking_chart",
    "create_pan_cohort_faceted_decay_chart",
    "create_cross_cohort_resilience_heatmap",
    "create_meta_analytic_decay_chart",
    "create_c_index_forest_plot",
    "create_cox_hr_forest_plot",
    "create_dca_net_benefit_chart",
    "export_chart_svg",
]
