"""Preprocessing routines: gene filtering, library size normalization, and metadata harmonization."""

from .gene_filtering import CONFOUNDING_PATTERNS, filter_confounding_genes
from .metadata import (
    binarize_response,
    harmonize_obs_metadata,
    standardize_recist,
    standardize_timepoint,
)
from .gene_normalization import batch_normalize_to_sparse_h5ad, normalize_dataset_to_ensembl
from .normalization import (
    compute_log1p_norm_matrix,
    compute_tpm_matrix,
    expm1_transform,
    log1p_transform,
    normalize_to_tpm,
    normalize_total_counts,
    standardize_dual_layers,
    standardize_processed_layers,
)
from .pipeline import process_single_cell_dataset
from .qc import (
    apply_quality_control,
    calculate_adaptive_thresholds,
    compute_qc_covariates,
    extract_qc_metrics_dataframe,
)
from .qc_plots import build_qc_dashboard, display_qc_plots_inline, export_qc_plots
from .sanity import run_sanity_normalization
from .sctransform import normalize_sctransform

__all__ = [
    "apply_quality_control",
    "calculate_adaptive_thresholds",
    "compute_qc_covariates",
    "extract_qc_metrics_dataframe",
    "build_qc_dashboard",
    "display_qc_plots_inline",
    "export_qc_plots",
    "compute_log1p_norm_matrix",
    "compute_tpm_matrix",
    "normalize_to_tpm",
    "standardize_processed_layers",
    "process_single_cell_dataset",
    "CONFOUNDING_PATTERNS",
    "filter_confounding_genes",
    "normalize_dataset_to_ensembl",
    "batch_normalize_to_sparse_h5ad",
    "normalize_total_counts",
    "log1p_transform",
    "expm1_transform",
    "standardize_dual_layers",
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
]
from .matrix_inspection import (
    ExpressionInspectionResult,
    ExpressionType,
    inspect_expression_type,
    tag_expression_metadata,
)
