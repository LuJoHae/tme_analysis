"""Immutable data models for tme_datasets."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence, overload
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some

from .types import (
    GeneIDType,
    HarmonizeMode,
    Modality,
    NBEstimationMethod,
    SCTransformFlavor,
)


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


class RegisteredDatasets(tuple[DatasetSpec, ...]):
    """Immutable collection of registered dataset specifications.

    Provides:
    - Dual indexing: by position (e.g. `specs[0]`) or by dataset ID (e.g. `specs["GSE120575"]`).
    - Transformation into a typed Polars DataFrame (`specs.to_polars()`).
    - List accessor methods for all attributes:
      - Non-optional: `ids()`, `titles()`, `modalities()`, `cancer_types()`, `platforms()`, `has_response_labels()`.
      - Optional: `organs(unwrapped=False)`, `n_samples_or_cells(...)`, `raw_source_urls(...)`,
        `checksums(...)`, `qc_specs(...)`, `local_paths(...)`.
    - Querying & filtering: `.get(id)`, `.filter(...)`.
    """

    def __new__(cls, specs: Iterable[DatasetSpec] = ()) -> RegisteredDatasets:
        return super().__new__(cls, tuple(specs))

    @overload
    def __getitem__(self, item: int) -> DatasetSpec: ...

    @overload
    def __getitem__(self, item: slice) -> RegisteredDatasets: ...

    @overload
    def __getitem__(self, item: str) -> DatasetSpec: ...

    def __getitem__(self, item: int | slice | str) -> DatasetSpec | RegisteredDatasets:
        if isinstance(item, str):
            for spec in self:
                if spec.id == item:
                    return spec
            raise KeyError(f"Dataset ID '{item}' not found in registered datasets")
        res = super().__getitem__(item)
        if isinstance(item, slice):
            return RegisteredDatasets(res)
        return res

    def __contains__(self, item: object) -> bool:
        if isinstance(item, str):
            return any(s.id == item for s in self)
        return super().__contains__(item)

    def get(self, dataset_id: str) -> Maybe[DatasetSpec]:
        """Look up a dataset specification by ID, returning Some(spec) or Nothing."""
        for spec in self:
            if spec.id == dataset_id:
                return Some(spec)
        return Nothing

    def filter(
        self,
        modality: Modality | None = None,
        cancer_type: str | None = None,
        has_response: bool | None = None,
    ) -> RegisteredDatasets:
        """Filter specifications by modality, cancer type, or response labels."""
        return RegisteredDatasets(
            spec
            for spec in self
            if (modality is None or spec.modality == modality)
            and (cancer_type is None or spec.cancer_type.lower() == cancer_type.lower())
            and (has_response is None or spec.has_response_labels == has_response)
        )

    # --- Attribute List Accessor Methods ---

    def ids(self) -> list[str]:
        """Return all registered dataset IDs as a list."""
        return [s.id for s in self]

    def titles(self) -> list[str]:
        """Return all registered dataset titles as a list."""
        return [s.title for s in self]

    def modalities(self) -> list[Modality]:
        """Return all dataset modalities as a list of Modality enums."""
        return [s.modality for s in self]

    def cancer_types(self) -> list[str]:
        """Return all cancer types as a list."""
        return [s.cancer_type for s in self]

    def platforms(self) -> list[str]:
        """Return all sequencing/assay platforms as a list."""
        return [s.platform for s in self]

    def has_response_labels(self) -> list[bool]:
        """Return boolean indicator of response labels for each dataset."""
        return [s.has_response_labels for s in self]

    def organs(self, unwrapped: bool = False) -> list[Maybe[str]] | list[str | None]:
        """Return organ annotations as list[Maybe[str]] (or list[str | None] if unwrapped=True)."""
        return [s.organ.value_or(None) if unwrapped else s.organ for s in self]

    def n_samples_or_cells(self, unwrapped: bool = False) -> list[Maybe[int]] | list[int | None]:
        """Return sample/cell counts as list[Maybe[int]] (or list[int | None] if unwrapped=True)."""
        return [s.n_samples_or_cells.value_or(None) if unwrapped else s.n_samples_or_cells for s in self]

    def raw_source_urls(self, unwrapped: bool = False) -> list[Maybe[str]] | list[str | None]:
        """Return source URLs as list[Maybe[str]] (or list[str | None] if unwrapped=True)."""
        return [s.raw_source_url.value_or(None) if unwrapped else s.raw_source_url for s in self]

    def checksums(self, unwrapped: bool = False) -> list[Maybe[ChecksumSpec]] | list[ChecksumSpec | None]:
        """Return checksum specifications as list[Maybe[ChecksumSpec]] (or unwrapped)."""
        return [s.checksum.value_or(None) if unwrapped else s.checksum for s in self]

    def qc_specs(self, unwrapped: bool = False) -> list[Maybe[QualityControlSpec]] | list[QualityControlSpec | None]:
        """Return quality control specs as list[Maybe[QualityControlSpec]] (or unwrapped)."""
        return [s.qc_spec.value_or(None) if unwrapped else s.qc_spec for s in self]

    def local_paths(self, unwrapped: bool = False) -> list[Maybe[Path]] | list[Path | None]:
        """Return local file paths as list[Maybe[Path]] (or list[Path | None] if unwrapped=True)."""
        return [s.local_path.value_or(None) if unwrapped else s.local_path for s in self]

    # --- Polars DataFrame Export ---

    def to_polars(self) -> pl.DataFrame:
        """Convert all registered specifications into a strongly typed Polars DataFrame."""
        return pl.DataFrame({
            "id": [s.id for s in self],
            "title": [s.title for s in self],
            "modality": [s.modality.value for s in self],
            "cancer_type": [s.cancer_type for s in self],
            "platform": [s.platform for s in self],
            "organ": [s.organ.value_or(None) for s in self],
            "has_response_labels": [s.has_response_labels for s in self],
            "n_samples_or_cells": [s.n_samples_or_cells.value_or(None) for s in self],
            "raw_source_url": [s.raw_source_url.value_or(None) for s in self],
            "checksum": [s.checksum.map(lambda c: c.expected_hash).value_or(None) for s in self],
            "has_qc_spec": [isinstance(s.qc_spec, Some) for s in self],
            "local_path": [s.local_path.map(str).value_or(None) for s in self],
        })

    def to_dataframe(self) -> pl.DataFrame:
        """Alias for to_polars()."""
        return self.to_polars()

    def __repr__(self) -> str:
        preview = ", ".join(s.id for s in self[:4])
        ellipsis = ", ..." if len(self) > 4 else ""
        return f"<RegisteredDatasets len={len(self)} [{preview}{ellipsis}]>"


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


class SCTransformConfig(BaseModel):
    """Configuration for SCTransform / Pearson residual normalization."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    flavor: SCTransformFlavor = SCTransformFlavor.ANALYTIC
    n_top_genes: Maybe[int] = Nothing
    theta: float = 100.0
    clip_residuals: bool = True
    max_residual: Maybe[float] = Nothing
    use_layer_as_x: bool = True


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
