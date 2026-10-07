"""tme_datasets: Unified single-cell, spatial, bulk immunotherapy datasets, and perturbation framework."""

from __future__ import annotations

__version__ = "0.1.0"

# High-level query API
from .query import (
    IATLAS_COMBINED_GROUPS,
    build_deconvolution_reference,
    get_dataset_metadata,
    list_datasets,
    list_preprocessed_datasets,
    load_combined_iatlas_cohorts,
    load_iatlas_cohort_or_combined,
    load_dataset,
    load_geneset_collection,
    query_datasets,
    query_preprocessed_datasets,
    sample_single_cell_cohorts,
    run_sampling_hpo,
    filter_candidate_cohorts,
    HPOSearchSpace,
    HPOTrialConfig,
    HPORunResult,
    TrialEvaluationResult,
    MalignantStrategy,
    MalignantSamplingConfig,
)
from .genesets.collections import (
    IMMUNE_CHECKPOINT_GENES,
    IMMUNOTHERAPY_GENE_PANEL,
)

# Logging
from .logging import (
    configure_logging,
    get_logger,
    set_log_level,
)

# Registry
from .registry import (
    DATASET_ALIASES,
    DATASET_REGISTRY,
    get_dataset_spec,
    list_preprocessed_datasets,
    list_registered_datasets,
    query_preprocessed_datasets,
    resolve_dataset_id,
)

# Config-driven Paths API
from .paths import (
    DataPathsConfig,
    find_dataset_h5ad,
    find_repo_root,
    get_data_paths,
    get_ensembl_dir,
    get_manual_download_dir,
    get_preprocessed_h5ad_path,
    get_raw_dataset_dir,
    get_reference_h5ad_path,
    get_scratch_dataset_dir,
)

# Core types & models
from .types import (
    CohortSamplingMode,
    DatasetProvider,
    GeneIDType,
    HarmonizeMode,
    Modality,
    NBEstimationMethod,
    PerturbationTransform,
    SCTransformFlavor,
    StorageBackend,
)
from .models import (
    ChecksumSpec,
    ClusterAnalysisSpec,
    DatasetSpec,
    GeneReconcileConfig,
    HarmonizeConfig,
    IntegrationMetricsResult,
    NegativeBinomialConfig,
    PerturbationConfig,
    PreprocessedDatasetSpec,
    PreprocessedDatasets,
    PseudobulkConfig,
    QualityControlSpec,
    RegisteredDatasets,
    SampledSingleCellResult,
    SanityConfig,
    SCTransformConfig,
    SingleCellProcessingResult,
    SingleCellProcessingSpec,
    SingleCellSamplingSpec,
    SubsampleSpec,
)

# Preprocessing & Normalization
from .preprocessing import (
    ExpressionInspectionResult,
    ExpressionType,
    apply_quality_control,
    batch_normalize_to_sparse_h5ad,
    binarize_response,
    build_qc_dashboard,
    calculate_adaptive_thresholds,
    compute_qc_covariates,
    compute_tpm_matrix,
    display_qc_plots_inline,
    export_qc_plots,
    extract_qc_metrics_dataframe,
    filter_confounding_genes,
    harmonize_obs_metadata,
    inspect_expression_type,
    normalize_dataset_to_ensembl,
    normalize_sctransform,
    normalize_to_tpm,
    normalize_total_counts,
    process_single_cell_dataset,
    run_sanity_normalization,
    standardize_processed_layers,
    standardize_recist,
    standardize_timepoint,
    tag_expression_metadata,
)

# Deconvolution Reference Building
from .deconvolution import (
    DeconvolutionReferenceConfig,
    DeconvolutionReferenceResult,
    build_deconvolution_reference,
    calculate_cnv_proxy_scores,
    detect_malignant_cells,
    export_to_bayesprism,
    export_to_instaprism,
)

# Transforms & perturbations
from .transforms import (
    ComposeTransforms,
    add_expression_jitter,
    compute_size_factors,
    fit_nb_empirical_bayes,
    fit_nb_mle,
    fit_nb_moments,
    in_silico_knockout,
    in_silico_overexpression,
    infer_dataset_nb_parameters,
    randomize_negative_binomial,
    simulate_dropout,
    subsample_cells,
    supersample_cells,
)

# Gene sets
from .genesets import (
    AYERS_T_CELL_INFLAMED_GEP,
    TME_MAJOR_MARKERS,
    TME_SUBTYPE_MARKERS,
    GeneSet,
    GeneSetCollection,
    compute_geneset_overlap,
    export_gmt,
    get_bagaev_core_collection,
    get_tme_major_lineage_collection,
    get_tme_subtype_collection,
    parse_gmt,
    score_geneset_auc,
    score_geneset_zscore,
)

# Simulation
from .simulation import simulate_pseudobulk

# Harmonization & metrics
from .harmonization import (
    align_and_concatenate,
    evaluate_integration_metrics,
)

# PyTorch bridge
from .torch import (
    TmeTorchDataset,
    create_tme_dataloader,
)

# Out-of-core storage
from .storage import (
    H5ADSparseIncrementalWriter,
    convert_to_zarr,
    load_backed,
    slice_backed_dataset,
)

# Gene reconciliation
from .genes import (
    detect_gene_id_type,
    map_gene_identifier,
    reconcile_genes,
    strip_gene_version,
)

# Providers
from .providers import (
    download_and_build_maynard_full,
    load_maynard,
)

# Download & Verification
from .download import (
    compute_file_hash,
    download_gdrive_file,
    verify_checksum,
)

# Multi-Cohort Sampling & Clustering
from .sampling import (
    compute_cohort_cell_allocations,
    generate_random_cluster_spec,
    generate_random_sampling_spec,
    resolve_sampled_cohorts,
    run_pca_knn_leiden,
    sample_and_harmonize_cohorts,
    sample_single_cell_cohorts,
)

__all__ = [
    "__version__",
    # Query
    "load_dataset",
    "load_combined_iatlas_cohorts",
    "load_iatlas_cohort_or_combined",
    "IATLAS_COMBINED_GROUPS",
    "IMMUNE_CHECKPOINT_GENES",
    "IMMUNOTHERAPY_GENE_PANEL",
    "query_datasets",
    "list_preprocessed_datasets",
    "query_preprocessed_datasets",
    "load_geneset_collection",
    "sample_single_cell_cohorts",
    "generate_random_sampling_spec",
    "generate_random_cluster_spec",
    "build_deconvolution_reference",
    "run_sampling_hpo",
    "filter_candidate_cohorts",
    "HPOSearchSpace",
    "HPOTrialConfig",
    "HPORunResult",
    "TrialEvaluationResult",
    "MalignantStrategy",
    "MalignantSamplingConfig",
    # Registry
    "DATASET_REGISTRY",
    "get_dataset_spec",
    "list_registered_datasets",
    # Paths API
    "DataPathsConfig",
    "find_repo_root",
    "get_data_paths",
    "get_ensembl_dir",
    "get_manual_download_dir",
    "get_preprocessed_h5ad_path",
    "get_raw_dataset_dir",
    "get_scratch_dataset_dir",
    "get_reference_h5ad_path",
    "find_dataset_h5ad",
    # Types & Models
    "Modality",
    "HarmonizeMode",
    "CohortSamplingMode",
    "GeneIDType",
    "StorageBackend",
    "NBEstimationMethod",
    "SCTransformFlavor",
    "DatasetProvider",
    "PerturbationTransform",
    "ChecksumSpec",
    "DatasetSpec",
    "PreprocessedDatasetSpec",
    "RegisteredDatasets",
    "PreprocessedDatasets",
    "DeconvolutionReferenceConfig",
    "DeconvolutionReferenceResult",
    "QualityControlSpec",
    "SingleCellProcessingSpec",
    "SingleCellProcessingResult",
    "SingleCellSamplingSpec",
    "ClusterAnalysisSpec",
    "SampledSingleCellResult",
    "SubsampleSpec",
    "NegativeBinomialConfig",
    "SanityConfig",
    "SCTransformConfig",
    "PerturbationConfig",
    "GeneReconcileConfig",
    "PseudobulkConfig",
    "HarmonizeConfig",
    "IntegrationMetricsResult",
    # Sampling & Graph Clustering
    "compute_cohort_cell_allocations",
    "generate_random_cluster_spec",
    "generate_random_sampling_spec",
    "resolve_sampled_cohorts",
    "run_pca_knn_leiden",
    "sample_and_harmonize_cohorts",
    # Registry
    "DATASET_REGISTRY",
    "DATASET_ALIASES",
    "resolve_dataset_id",
    "get_dataset_spec",
    "list_registered_datasets",
    "list_preprocessed_datasets",
    "query_preprocessed_datasets",
    # Preprocessing
    "apply_quality_control",
    "calculate_adaptive_thresholds",
    "compute_qc_covariates",
    "extract_qc_metrics_dataframe",
    "build_qc_dashboard",
    "display_qc_plots_inline",
    "export_qc_plots",
    "process_single_cell_dataset",
    "normalize_to_tpm",
    "standardize_processed_layers",
    "normalize_total_counts",
    "normalize_dataset_to_ensembl",
    "batch_normalize_to_sparse_h5ad",
    "filter_confounding_genes",
    "binarize_response",
    "standardize_recist",
    "standardize_timepoint",
    "harmonize_obs_metadata",
    "run_sanity_normalization",
    "normalize_sctransform",
    "inspect_expression_type",
    "tag_expression_metadata",
    "ExpressionType",
    "ExpressionInspectionResult",
    # Deconvolution Reference Building
    "detect_malignant_cells",
    "calculate_cnv_proxy_scores",
    "export_to_bayesprism",
    "export_to_instaprism",
    # Transforms
    "subsample_cells",
    "supersample_cells",
    "randomize_negative_binomial",
    "simulate_dropout",
    "add_expression_jitter",
    "in_silico_knockout",
    "in_silico_overexpression",
    "ComposeTransforms",
    "compute_size_factors",
    "fit_nb_moments",
    "fit_nb_mle",
    "fit_nb_empirical_bayes",
    "infer_dataset_nb_parameters",
    # Gene sets
    "GeneSet",
    "GeneSetCollection",
    "AYERS_T_CELL_INFLAMED_GEP",
    "TME_MAJOR_MARKERS",
    "TME_SUBTYPE_MARKERS",
    "get_tme_major_lineage_collection",
    "get_tme_subtype_collection",
    "get_bagaev_core_collection",
    "parse_gmt",
    "export_gmt",
    "score_geneset_zscore",
    "score_geneset_auc",
    "compute_geneset_overlap",
    # Simulation
    "simulate_pseudobulk",
    # Harmonization
    "align_and_concatenate",
    "evaluate_integration_metrics",
    # PyTorch
    "TmeTorchDataset",
    "create_tme_dataloader",
    # Storage
    "load_backed",
    "slice_backed_dataset",
    "convert_to_zarr",
    "H5ADSparseIncrementalWriter",
    # Gene reconciliation
    "detect_gene_id_type",
    "map_gene_identifier",
    "reconcile_genes",
    "strip_gene_version",
    # Verification & Download
    "compute_file_hash",
    "verify_checksum",
    "download_gdrive_file",
    # Providers
    "load_maynard",
    "download_and_build_maynard_full",
    # Logging
    "configure_logging",
    "get_logger",
    "set_log_level",
]
