"""Immutable data models for tme_datasets."""

from __future__ import annotations

from pathlib import Path
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing

from .types import GeneIDType, HarmonizeMode, Modality, NBEstimationMethod


class ChecksumSpec(BaseModel):
    """Cryptographic hash specification for data integrity verification."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    algorithm: str = "sha256"
    expected_hash: str
    size_bytes: Maybe[int] = Nothing


class QualityControlSpec(BaseModel):
    """Quality control thresholds for single-cell or bulk filtering."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    min_genes_per_cell: int = 200
    max_genes_per_cell: int = 8000
    min_counts_per_cell: int = 500
    max_pct_mitochondrial: float = 20.0
    filter_confounding_genes: bool = True


class DatasetSpec(BaseModel):
    """Immutable metadata specification for a registered dataset."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    id: str
    title: str
    modality: Modality
    cancer_type: str
    platform: str
    organ: Maybe[str] = Nothing
    has_response_labels: bool = False
    n_samples_or_cells: Maybe[int] = Nothing
    raw_source_url: Maybe[str] = Nothing
    checksum: Maybe[ChecksumSpec] = Nothing
    qc_spec: Maybe[QualityControlSpec] = Nothing
    local_path: Maybe[Path] = Nothing


class SubsampleSpec(BaseModel):
    """Specification for subsampling or supersampling cells."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    n_or_fraction: float
    stratify_by: Maybe[str] = Nothing
    balanced: bool = False
    seed: Maybe[int] = Nothing


class NegativeBinomialConfig(BaseModel):
    """Configuration for Negative Binomial count randomization."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    dispersion: float = 0.1
    estimation_method: Maybe[NBEstimationMethod] = Nothing
    cluster_key: Maybe[str] = Nothing
    library_size_scaling: bool = True
    baseline_rate: float = 0.0
    min_dispersion: float = 1e-4
    max_dispersion: float = 10.0
    seed: Maybe[int] = Nothing


class SanityConfig(BaseModel):
    """Configuration for Sanity Bayesian Log-Normal Poisson normalization and denoising."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    v_min: float = 0.001
    v_max: float = 20.0
    n_bins: int = 40
    seed: Maybe[int] = Nothing


class PerturbationConfig(BaseModel):
    """Configuration for composite dataset perturbations."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    nb_config: Maybe[NegativeBinomialConfig] = Nothing
    dropout_rate: Maybe[float] = Nothing
    gaussian_jitter_sigma: Maybe[float] = Nothing
    seed: Maybe[int] = Nothing


class GeneReconcileConfig(BaseModel):
    """Configuration for gene symbol or Ensembl ID normalization."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    target_type: GeneIDType = GeneIDType.HUGO_SYMBOL
    ensembl_release: int = 110
    strip_version_suffix: bool = True
    handle_duplicates: str = "sum"


class PseudobulkConfig(BaseModel):
    """Configuration for in-silico pseudobulk mixture simulation."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    n_samples: int = 50
    cells_per_sample: int = 2000
    target_depth: Maybe[int] = Nothing
    noise_dispersion: Maybe[float] = Nothing
    seed: Maybe[int] = Nothing


class HarmonizeConfig(BaseModel):
    """Configuration for cross-dataset harmonization."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    mode: HarmonizeMode = HarmonizeMode.INTERSECTION
    batch_key: str = "dataset_id"
    reconcile_genes: bool = True
    gene_target_type: GeneIDType = GeneIDType.HUGO_SYMBOL
    min_shared_genes: int = 500


class IntegrationMetricsResult(BaseModel):
    """Quantitative evaluation metrics for dataset harmonization."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    mean_ilisi: float
    mean_clisi: float
    batch_silhouette: float
    cell_type_silhouette: float
    silhouette_ratio: float
    kbet_acceptance_rate: float
