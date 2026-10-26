"""Strict type stub for `stage3b_bridge.py`."""

from __future__ import annotations

from pathlib import Path
from typing import Final, Sequence

import numpy as np
import numpy.typing as npt
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Result

COHORT_ORDER: Final[tuple[dict[str, str], ...]]


class UnifiedBenchmarkConfig(BaseModel):
    model_config = ConfigDict(frozen=True)

    benchmark_dir: Path = ...
    out_dir: Path = ...
    results_dir: Path = ...
    n_meta_states: int = ...
    softmax_temp: float = ...
    deconv_iters: int = ...
    n_permutations: int = ...
    seed: int = ...


class CohortEmpiricalBundle(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    cohort_id: str
    cancer_type: str
    author: str
    short_label: str
    year: str
    clusters: tuple[str, ...]
    phi_genes: tuple[str, ...]
    phi_mat: npt.NDArray[np.float64]
    sample_ids: tuple[str, ...]
    bulk_genes: tuple[str, ...]
    bulk_mat: npt.NDArray[np.float64]
    fractions_df: pl.DataFrame
    self_concordance_df: pl.DataFrame


class UnifiedReferenceAtlas(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    global_state_ids: tuple[str, ...]
    global_cohorts: tuple[str, ...]
    global_local_clusters: tuple[str, ...]
    union_genes: tuple[str, ...]
    similarity_matrix: npt.NDArray[np.float64]
    meta_state_labels: tuple[str, ...]
    phi_meta_df: pl.DataFrame


def load_all_nine_cohorts(
    config: UnifiedBenchmarkConfig,
) -> Result[tuple[CohortEmpiricalBundle, ...], str]: ...


def compute_within_cohort_specificity(
    phi_mat: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]: ...


def build_unified_reference_atlas(
    bundles: Sequence[CohortEmpiricalBundle],
    n_meta_states: int = ...,
) -> UnifiedReferenceAtlas: ...


def compute_permutation_pvalue(
    x_vec: npt.NDArray[np.float64],
    y_vec: npt.NDArray[np.float64],
    n_permutations: int = ...,
    seed: int = ...,
) -> tuple[float, float, float]: ...
