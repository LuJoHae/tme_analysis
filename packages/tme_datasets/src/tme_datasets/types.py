"""Core protocol definitions, enums, and typeclasses for tme_datasets."""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Protocol, TypeVar, runtime_checkable
import anndata as ad
from returns.result import Result

T = TypeVar("T")


class Modality(str, Enum):
    """Data modality of the transcriptomic dataset."""

    SINGLE_CELL = "single_cell"
    BULK_RNA = "bulk_rna"
    SPATIAL = "spatial"
    GENE_SET = "gene_set"


class HarmonizeMode(str, Enum):
    """Method for harmonizing gene feature spaces across disparate datasets."""

    INTERSECTION = "intersection"
    UNION_ZERO_FILLED = "union_zero_filled"


class GeneIDType(str, Enum):
    """Target or source gene nomenclature type."""

    HUGO_SYMBOL = "hugo_symbol"
    ENSEMBL_ID = "ensembl_id"
    ENTREZ_ID = "entrez_id"
    AUTO_DETECT = "auto_detect"


class StorageBackend(str, Enum):
    """Storage backend for dataset persistence and streaming."""

    MEMORY = "memory"
    BACKED_H5AD = "backed_h5ad"
    ZARR = "zarr"


class NBEstimationMethod(str, Enum):
    """Estimation method for Negative Binomial parameters and Bayesian count modeling."""

    MOMENTS = "moments"
    MLE = "mle"
    EMPIRICAL_BAYES = "empirical_bayes"
    SANITY = "sanity"


class SCTransformFlavor(str, Enum):
    """Flavor of variance-stabilizing transformation."""

    ANALYTIC = "analytic"
    REGULARIZED_GLM = "regularized_glm"


@runtime_checkable
class DatasetProvider(Protocol):
    """Typeclass protocol for dataset loading and processing."""

    def download(self, target_dir: Path, force: bool = False) -> Result[Path, str]:
        """Download raw dataset assets into target_dir."""
        ...

    def preprocess(self, raw_path: Path, output_path: Path) -> Result[ad.AnnData, str]:
        """Preprocess raw assets into a clean AnnData object and save to output_path."""
        ...

    def load(self, base_dir: Path, auto_download: bool = True) -> Result[ad.AnnData, str]:
        """Load dataset, downloading and preprocessing if not present on disk."""
        ...


@runtime_checkable
class PerturbationTransform(Protocol):
    """Typeclass protocol for pure AnnData perturbation transforms."""

    def __call__(self, adata: ad.AnnData) -> Result[ad.AnnData, str]:
        """Transform AnnData returning a new modified copy without in-place mutation."""
        ...
