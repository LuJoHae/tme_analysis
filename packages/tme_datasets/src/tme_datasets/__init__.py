"""tme_datasets: Unified single-cell, spatial, bulk immunotherapy datasets, and perturbation framework."""

from __future__ import annotations

__version__ = "0.1.0"

# High-level query API
from .query import (
    load_dataset,
    load_geneset_collection,
    query_datasets,
)

# Registry
from .registry import (
    DATASET_REGISTRY,
    get_dataset_spec,
    list_registered_datasets,
)

# Core types & models
from .types import (
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
    DatasetSpec,
    GeneReconcileConfig,
    HarmonizeConfig,
    IntegrationMetricsResult,
    NegativeBinomialConfig,
    PerturbationConfig,
    PseudobulkConfig,
    QualityControlSpec,
    SanityConfig,
    SCTransformConfig,
    SubsampleSpec,
)

# Preprocessing & Normalization
from .preprocessing import (
    binarize_response,
    filter_confounding_genes,
    harmonize_obs_metadata,
    normalize_sctransform,
    normalize_total_counts,
    run_sanity_normalization,
    standardize_recist,
    standardize_timepoint,
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

# Cryptographic verification
from .download import (
    compute_file_hash,
    verify_checksum,
)

__all__ = [
    "__version__",
    # Query
    "load_dataset",
    "query_datasets",
    "load_geneset_collection",
    # Registry
    "DATASET_REGISTRY",
    "get_dataset_spec",
    "list_registered_datasets",
    # Types & Models
    "Modality",
    "HarmonizeMode",
    "GeneIDType",
    "StorageBackend",
    "NBEstimationMethod",
    "SCTransformFlavor",
    "DatasetProvider",
    "PerturbationTransform",
    "ChecksumSpec",
    "DatasetSpec",
    "QualityControlSpec",
    "SubsampleSpec",
    "NegativeBinomialConfig",
    "SanityConfig",
    "SCTransformConfig",
    "PerturbationConfig",
    "GeneReconcileConfig",
    "PseudobulkConfig",
    "HarmonizeConfig",
    "IntegrationMetricsResult",
    # Preprocessing
    "normalize_total_counts",
    "filter_confounding_genes",
    "binarize_response",
    "standardize_recist",
    "standardize_timepoint",
    "harmonize_obs_metadata",
    "run_sanity_normalization",
    "normalize_sctransform",
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
    # Gene reconciliation
    "detect_gene_id_type",
    "map_gene_identifier",
    "reconcile_genes",
    "strip_gene_version",
    # Verification
    "compute_file_hash",
    "verify_checksum",
]
