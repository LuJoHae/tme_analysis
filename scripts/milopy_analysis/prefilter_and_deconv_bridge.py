#!/usr/bin/env python3
"""Outcome-Blind Pairwise Milopy Prefiltering Benchmark & Bulk/Pseudobulk Deconvolution Bridge.

Executes two coupled investigations across all 9 pre-analyzed clinical single-cell immunotherapy
cohorts (`output/results/milopy_response_analysis/`, 32,680 Milo neighborhoods, 611 donors) and
the pan-cancer bulk RNA-seq deconvolution benchmark (`Pan-iAtlas`, N = 1,015 ICB patients,
42 Consensus-MetaStates M01-M42):

1. Statistical Power Enhancement via Strictly Outcome-Blind Prefiltering (`T0` through `T5`):
   - T0_unfiltered: Baseline all-vs-all disjoint neighborhood screening
   - T1_tarone_exact: Exact discrete Tarone hypergeometric corner pruning on Berger's IUT
   - T2_marginal_entropy_tarone: Bourgon-Gentleman-Huber marginal prevalence/variance + Tarone
   - T3_cross_meta42_tarone: Cross-Meta42 lineage scope + marginal dispersion + Tarone
   - T4_super_nhoods_tarone: Meta42-preserving Leiden super-neighborhood cell-union pooling + Tarone
   - T5_bulk_deconv_prior: T4 + independent Pan-iAtlas (N=1,015) Genovese-Roeder-Wasserman prior weights

2. Mechanistic Bridge Between Single-Cell Pairwise Milopy & Deconvolution Cell-State Interactions:
   - Strictly outcome-blind donor-containment + arcsin-profile + kNN-smoothed neighborhood -> L_k -> Meta42
     mapping across all 32,680 neighborhoods (validated against direct .h5ad expression on GSE120575)
   - Cross-modality 42x42 & focused 8x8 Mantel node-permutation concordance against Pan-iAtlas bulk
     epistasis_z (Synergy M29 Plasma x M37 TLS B / M13 Tex x M14 M1-TAM vs. Antagonistic M17 myCAF gate)
   - Matched within-cohort concordance against 3a InstaPrism pseudobulk inferred_fraction pairs
   - Direct single-cell whole-tumor biopsy validation of the M17 ISLR+COL5A2+ myCAF stromal exclusion gate

Architecture strictly follows functional Python standards (.agents/rules/code-style-guide.md):
- Functional Core, Imperative Shell (pure data/chart transformations, IO isolated at shell)
- Frozen Pydantic models (`ConfigDict(frozen=True)`), zero `None` / `Optional` (`returns.maybe.Maybe`)
- Polymorphism via `typing.Protocol` typeclasses (dependency-injected loaders, validators, exporters)
- Exhaustive structural pattern matching with `typing.assert_never`
- Railway error handling via `returns.result.Result`
- Polars DataFrames exclusively and Altair vector SVG exports (Nature Methods minimalist wireframe)
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from enum import StrEnum
from pathlib import Path
import sys
import time
from typing import Protocol, TypeAlias, assert_never

import altair as alt
import anndata as ad
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import sparse as sp
from scipy.stats import norm, spearmanr

from milopy.pairwise import (
    ContinuousAxis,
    EnrichmentDirection,
    OverlapConfig,
    PairModuleConfig,
    PairScopeFilter,
    PairwiseAnalysisResult,
    PairwiseConfig,
    PairwiseMethod,
    PrefilterConfig,
    PrefilterMode,
    PresenceConfig,
    PresenceMode,
    QuadrantState,
    RefinementConfig,
    ScreenConfig,
    analyze_nhood_pairs,
    build_external_prior_weight_matrix,
)
from milopy.pairwise.types import FloatArray, IntArray, SparseMatrix
from plotting_utils import (  # type: ignore[import-untyped]
    COHORT_PANELS,
    COLOR_NEUTRAL_GREY,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    OKABE_BLUE,
    OKABE_BLUISH_GREEN,
    OKABE_ORANGE,
    OKABE_REDDISH_PURPLE,
    OKABE_SKY_BLUE,
    OKABE_VERMILION,
    apply_nature_methods_theme,
    export_altair_figure,
    get_single_cell_cohort_display_label,
)

AltairChart: TypeAlias = alt.Chart | alt.LayerChart | alt.HConcatChart | alt.VConcatChart

FOCUSED_8_META_STATES: tuple[str, ...] = (
    "M04",
    "M05",
    "M17",
    "M20",
    "M25",
    "M29",
    "M30",
    "M32",
)

LYMPHOID_EFFECTOR_META_STATES: frozenset[str] = frozenset(
    {"M04", "M05", "M13", "M14", "M15", "M24", "M29", "M30", "M34", "M36", "M37", "M41"}
)
STROMAL_CAF_BARRIER_META_STATES: frozenset[str] = frozenset({"M17", "M23"})
STROMAL_MYELOID_BARRIER_META_STATES: frozenset[str] = frozenset(
    {"M12", "M17", "M19", "M23", "M28", "M32"}
)


# =============================================================================
# 1. Exhaustive Sum Types & Frozen Domain Models
# =============================================================================


class LineageCompartment(StrEnum):
    """Canonical broad cellular lineage compartments bridging single-cell & Meta42 states."""

    T_NK = "T_NK"
    B_CELL = "B_CELL"
    PLASMA = "PLASMA"
    MYELOID_MAST = "MYELOID_MAST"
    STROMAL_CAF = "STROMAL_CAF"
    ENDOTHELIAL = "ENDOTHELIAL"
    EPITHELIAL_TUMOR = "EPITHELIAL_TUMOR"
    OTHER_UNSPECIFIED = "OTHER_UNSPECIFIED"


class AblationTier(StrEnum):
    """Six progressive strictly outcome-blind prefiltering and prior-weighting tiers."""

    T0_UNFILTERED = "T0_unfiltered"
    T1_TARONE_EXACT = "T1_tarone_exact"
    T2_MARGINAL_ENTROPY_TARONE = "T2_marginal_entropy_tarone"
    T3_CROSS_META42_TARONE = "T3_cross_meta42_tarone"
    T4_SUPER_NHOODS_TARONE = "T4_super_nhoods_tarone"
    T5_BULK_DECONV_PRIOR = "T5_bulk_deconv_prior"


ALL_ABLATION_TIERS: tuple[AblationTier, ...] = (
    AblationTier.T0_UNFILTERED,
    AblationTier.T1_TARONE_EXACT,
    AblationTier.T2_MARGINAL_ENTROPY_TARONE,
    AblationTier.T3_CROSS_META42_TARONE,
    AblationTier.T4_SUPER_NHOODS_TARONE,
    AblationTier.T5_BULK_DECONV_PRIOR,
)


def ablation_tier_display_label(tier: AblationTier) -> str:
    """Human-readable publication axis label for each prefiltering tier."""
    match tier:
        case AblationTier.T0_UNFILTERED:
            return "T0: Unfiltered All-vs-All"
        case AblationTier.T1_TARONE_EXACT:
            return "T1: Exact Tarone IUT Bound"
        case AblationTier.T2_MARGINAL_ENTROPY_TARONE:
            return "T2: Marginal Entropy + Var + Tarone"
        case AblationTier.T3_CROSS_META42_TARONE:
            return "T3: Cross-Meta42 Scope + Tarone"
        case AblationTier.T4_SUPER_NHOODS_TARONE:
            return "T4: Super-Nhood Pooling + Tarone"
        case AblationTier.T5_BULK_DECONV_PRIOR:
            return "T5: Super-Nhood + Bulk Deconv Prior"
        case _ as unreachable:
            assert_never(unreachable)


class BridgeConfig(BaseModel):
    """Immutable configuration for the 9-cohort prefilter benchmark & deconvolution bridge."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    milo_input_dir: Path = Field(default=Path("output/results/milopy_response_analysis"))
    deconv_3a_dir: Path = Field(
        default=Path("output/manuscript_pipeline/3a_cross_cohort_pseudobulk_deconv_benchmark")
    )
    deconv_3b_path: Path = Field(
        default=Path(
            "output/manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark/"
            "unified_all_states_detailed_concordance.parquet"
        )
    )
    iatlas_dir: Path = Field(default=Path("output/results/iatlas_instaprism_benchmark"))
    local_h5ad_dir: Path = Field(default=Path("data/preprocessed"))
    output_dir: Path = Field(default=Path("output/results/milopy_pairwise_deconv_bridge"))
    cohort_filter: Maybe[tuple[str, ...]] = Field(default=Nothing)
    fdr_threshold: float = Field(default=0.10, gt=0.0, le=1.0)
    tarone_p_ceiling: float = Field(default=0.005, gt=0.0, le=1.0)
    presence_mode: PresenceMode = Field(default=PresenceMode.HARD_DEPTH_GATED)
    min_donors_per_quadrant: float = Field(default=3.0, ge=1.0)
    max_refine_pairs: int = Field(default=100, ge=0)
    chunk_size: int = Field(default=512, ge=16)
    n_mantel_perms: int = Field(default=5000, ge=100)
    mantel_seed: int = Field(default=2026)
    render_plots: bool = Field(default=True)
    reuse_cached_ablation: bool = Field(default=False)


class CohortInputsBundle(BaseModel):
    """Immutable single-cell Milo + 3a pseudobulk benchmark artifacts for one cohort."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    accession: str
    nhoods_mat: SparseMatrix
    cells_df: pl.DataFrame
    da_results_df: pl.DataFrame
    pseudobulk_fractions_df: pl.DataFrame
    reference_phi_df: pl.DataFrame


class GlobalDeconvInputs(BaseModel):
    """Immutable pan-cohort 3b concordance and Pan-iAtlas (N=1,015) bulk interaction tables."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    concordance_3b_df: pl.DataFrame
    iatlas_all42_interactions_df: pl.DataFrame
    iatlas_focused8_interactions_df: pl.DataFrame
    iatlas_hub_summary_df: pl.DataFrame
    iatlas_top_degs_df: pl.DataFrame


class CohortEvaluationDesign(BaseModel):
    """Immutable lesion/donor evaluation design matrix C (V_nhoods x J_evals)."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    accession: str
    eval_ids: tuple[str, ...]
    eval_patient_ids: tuple[str, ...]
    response_labels: tuple[str, ...]
    count_matrix: FloatArray
    library_sizes: FloatArray
    cell_eval_codes: IntArray
    cell_types: tuple[str, ...]
    n_cells: int
    n_nhoods: int
    n_responders: int
    n_non_responders: int


class NhoodMeta42Mapping(BaseModel):
    """Outcome-blind neighborhood -> L_k -> Consensus-Meta42 mapping for one cohort."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    accession: str
    annotations_df: pl.DataFrame
    nhood_meta_states: tuple[str, ...]
    nhood_compartments: tuple[str, ...]
    has_stromal_caf: bool
    n_unique_lk: int
    n_total_lk: int
    n_unique_meta42: int
    posthoc_milo_logfc_r: float


class MantelTestResult(BaseModel):
    """Exact node-permutation Mantel correlation summary between two symmetric state matrices."""

    model_config = ConfigDict(frozen=True)

    panel_name: str
    n_states: int
    n_pairs: int
    pearson_r: float
    spearman_rho: float
    pvalue_two_sided: float
    pvalue_one_sided_pos: float
    directional_concordance_strong: float
    n_strong_pairs: int


class H5adValidationMetrics(BaseModel):
    """Direct single-cell .h5ad expression validation metrics for the outcome-blind mapper."""

    model_config = ConfigDict(frozen=True)

    accession: str
    n_nhoods: int
    n_shared_genes: int
    exact_lk_agreement: float
    exact_meta42_agreement: float
    exact_compartment_agreement: float


class CohortBridgeExecution(BaseModel):
    """Complete execution outputs for a single cohort across all 6 prefilter tiers."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    accession: str
    cohort_label: str
    design: CohortEvaluationDesign
    mapping: NhoodMeta42Mapping
    ablation_rows_df: pl.DataFrame
    primary_t4_result: PairwiseAnalysisResult
    prior_t5_result: PairwiseAnalysisResult


class BridgePipelineOutputs(BaseModel):
    """All pure DataFrame and statistical summary deliverables of the bridge pipeline."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    nhood_annotations_df: pl.DataFrame
    ablation_benchmark_df: pl.DataFrame
    meta42_concordance_df: pl.DataFrame
    hub_concordance_df: pl.DataFrame
    within_cohort_pseudobulk_df: pl.DataFrame
    mycaf_gating_regimes_df: pl.DataFrame
    mycaf_cohort_forest_df: pl.DataFrame
    mantel_all42: MantelTestResult
    mantel_whole_tumor: MantelTestResult
    mantel_focused8: MantelTestResult
    h5ad_validation: Maybe[H5adValidationMetrics]
    markdown_report: str


# =============================================================================
# 2. Typeclass Protocols (Side-Effecting IO Boundaries)
# =============================================================================


class CohortArtifactLoader(Protocol):
    """Loads cached Milo matrices and 3a pseudobulk tables for a single cohort."""

    def __call__(
        self,
        accession: str,
        milo_dir: Path,
        deconv_3a_dir: Path,
    ) -> Result[CohortInputsBundle, str]: ...


class GlobalDeconvLoader(Protocol):
    """Loads global 3b concordance and Pan-iAtlas bulk interaction tables."""

    def __call__(
        self,
        deconv_3b_path: Path,
        iatlas_dir: Path,
    ) -> Result[GlobalDeconvInputs, str]: ...


class LocalH5adValidator(Protocol):
    """Optional side-effecting validator comparing donor-profile mapping to direct .h5ad expression."""

    def __call__(
        self,
        local_h5ad_dir: Path,
        cohort_inputs: CohortInputsBundle,
        mapping: NhoodMeta42Mapping,
        lk_to_meta: Mapping[str, str],
    ) -> Maybe[H5adValidationMetrics]: ...


class VectorFigureExporter(Protocol):
    """Exports an Altair chart to vector SVG and 300 DPI verification PNG."""

    def __call__(self, chart: AltairChart, base_path: Path) -> Result[tuple[Path, Path], str]: ...


# =============================================================================
# 3. Pure Functional Core — Lineage Compartments & Outcome-Blind Meta42 Mapper
# =============================================================================


_META42_COMPARTMENT_MAP: Mapping[str, LineageCompartment] = {
    "M01": LineageCompartment.T_NK,
    "M02": LineageCompartment.STROMAL_CAF,
    "M03": LineageCompartment.T_NK,
    "M04": LineageCompartment.T_NK,
    "M05": LineageCompartment.T_NK,
    "M06": LineageCompartment.T_NK,
    "M07": LineageCompartment.OTHER_UNSPECIFIED,
    "M08": LineageCompartment.T_NK,
    "M09": LineageCompartment.MYELOID_MAST,
    "M10": LineageCompartment.MYELOID_MAST,
    "M11": LineageCompartment.EPITHELIAL_TUMOR,
    "M12": LineageCompartment.EPITHELIAL_TUMOR,
    "M13": LineageCompartment.T_NK,
    "M14": LineageCompartment.MYELOID_MAST,
    "M15": LineageCompartment.T_NK,
    "M16": LineageCompartment.ENDOTHELIAL,
    "M17": LineageCompartment.STROMAL_CAF,
    "M18": LineageCompartment.EPITHELIAL_TUMOR,
    "M19": LineageCompartment.MYELOID_MAST,
    "M20": LineageCompartment.T_NK,
    "M21": LineageCompartment.T_NK,
    "M22": LineageCompartment.EPITHELIAL_TUMOR,
    "M23": LineageCompartment.STROMAL_CAF,
    "M24": LineageCompartment.T_NK,
    "M25": LineageCompartment.T_NK,
    "M26": LineageCompartment.EPITHELIAL_TUMOR,
    "M27": LineageCompartment.EPITHELIAL_TUMOR,
    "M28": LineageCompartment.MYELOID_MAST,
    "M29": LineageCompartment.PLASMA,
    "M30": LineageCompartment.B_CELL,
    "M31": LineageCompartment.MYELOID_MAST,
    "M32": LineageCompartment.MYELOID_MAST,
    "M33": LineageCompartment.T_NK,
    "M34": LineageCompartment.T_NK,
    "M35": LineageCompartment.MYELOID_MAST,
    "M36": LineageCompartment.T_NK,
    "M37": LineageCompartment.B_CELL,
    "M38": LineageCompartment.MYELOID_MAST,
    "M39": LineageCompartment.MYELOID_MAST,
    "M40": LineageCompartment.T_NK,
    "M41": LineageCompartment.MYELOID_MAST,
    "M42": LineageCompartment.MYELOID_MAST,
}

_ALL_COMPARTMENTS: frozenset[LineageCompartment] = frozenset(LineageCompartment)


def classify_meta42_compartment(meta_state: str) -> LineageCompartment:
    """Maps a Consensus-Meta42 identifier (M01-M42) to its canonical cellular compartment."""
    return _META42_COMPARTMENT_MAP.get(meta_state.strip(), LineageCompartment.OTHER_UNSPECIFIED)


def compatible_compartments_for_nhood_celltype(
    nhood_ct: str,
    cohort_has_explicit_plasma_label: bool,
) -> frozenset[LineageCompartment]:
    """Determines compatible cellular compartments for an observed Milo Nhood_CellType string.

    Outcome-blind ontology compatibility: when `Nhood_CellType` is 'Unspecified' or numeric
    cluster IDs ('1'..'7'), all compartments remain eligible. When a cohort annotates broad
    'B cell' without a dedicated plasma label (e.g. GSE243013), 'B cell' admits both B_CELL
    and PLASMA so M29 JCHAIN+MZB1+ Plasma states are preserved cleanly.
    """
    s = nhood_ct.strip().lower()
    if (not s) or s in ("unspecified", "unannotated", "unknown") or s.isdigit():
        return _ALL_COMPARTMENTS
    if "plasma" in s and "plasmacytoid" not in s:
        return frozenset({LineageCompartment.PLASMA, LineageCompartment.B_CELL})
    if "plasmablast" in s:
        return frozenset({LineageCompartment.PLASMA, LineageCompartment.B_CELL})
    if "b lineage" in s:
        return frozenset({LineageCompartment.B_CELL, LineageCompartment.PLASMA})
    if "b cell" in s:
        return (
            frozenset({LineageCompartment.B_CELL})
            if cohort_has_explicit_plasma_label
            else frozenset({LineageCompartment.B_CELL, LineageCompartment.PLASMA})
        )
    if (
        "t cell" in s
        or "t-helper" in s
        or "natural killer" in s
        or "nk " in s
        or "t/nk" in s
        or "regulatory t" in s
        or "follicular helper" in s
    ):
        return frozenset({LineageCompartment.T_NK, LineageCompartment.EPITHELIAL_TUMOR})
    if "lymphocyte" == s:
        return frozenset(
            {LineageCompartment.T_NK, LineageCompartment.B_CELL, LineageCompartment.PLASMA}
        )
    if (
        "macrophage" in s
        or "monocyte" in s
        or "dendritic" in s
        or "myeloid" in s
        or "neutrophil" in s
        or "mast" in s
        or "granulocyte" in s
    ):
        return frozenset({LineageCompartment.MYELOID_MAST})
    if (
        "fibroblast" in s
        or "myofibroblast" in s
        or "pericyte" in s
        or "perivascular" in s
        or "stromal" in s
        or "schwann" in s
    ):
        return frozenset({LineageCompartment.STROMAL_CAF, LineageCompartment.ENDOTHELIAL})
    if "endothelial" in s:
        return frozenset({LineageCompartment.ENDOTHELIAL, LineageCompartment.STROMAL_CAF})
    if (
        "malignant" in s
        or "abnormal cell" in s
        or "enterocyte" in s
        or "colonocyte" in s
        or "goblet" in s
        or "transit amplifying" in s
        or "epithelial" in s
        or "tumor" in s
    ):
        return frozenset({LineageCompartment.EPITHELIAL_TUMOR})
    return _ALL_COMPARTMENTS


def resolve_target_cohorts(cohort_filter: Maybe[tuple[str, ...]]) -> tuple[str, ...]:
    """Resolves the ordered tuple of cohort accessions to analyze."""
    canonical = tuple(str(panel["acc"]) for panel in COHORT_PANELS)
    match cohort_filter:
        case Maybe.empty:
            return canonical
        case Some(selected):
            cleaned = tuple(s.strip() for s in selected if s.strip())
            return cleaned if cleaned else canonical
        case _ as unreachable:
            return canonical if unreachable else canonical


def build_cohort_evaluation_design(
    inputs: CohortInputsBundle,
) -> Result[CohortEvaluationDesign, str]:
    """Constructs the immutable (V_nhoods x J_evals) count matrix and cell evaluation codes."""
    required_cols = {"patient_id", "clinical_response"}
    missing = sorted(required_cols - set(inputs.cells_df.columns))
    if missing:
        return Failure(f"[{inputs.accession}] cells_df missing columns: {missing}")

    n_cells = inputs.cells_df.height
    mat_rows, n_nhoods = int(inputs.nhoods_mat.shape[0]), int(inputs.nhoods_mat.shape[1])
    if n_cells != mat_rows:
        return Failure(
            f"[{inputs.accession}] cells_df rows ({n_cells}) != nhoods_mat rows ({mat_rows})"
        )
    if inputs.da_results_df.height != n_nhoods:
        return Failure(
            f"[{inputs.accession}] da_results rows ({inputs.da_results_df.height}) != nhoods ({n_nhoods})"
        )

    norm_cells = (
        inputs.cells_df.lazy()
        .with_columns(
            pl.col("patient_id").cast(pl.String).str.strip_chars().alias("_pt"),
            pl.col("clinical_response")
            .cast(pl.String)
            .str.strip_chars()
            .str.to_lowercase()
            .alias("_resp"),
            (
                pl.col("cell_type").cast(pl.String)
                if "cell_type" in inputs.cells_df.columns
                else pl.lit("unannotated")
            ).alias("_ct"),
        )
        .with_columns((pl.col("_pt") + pl.lit("__") + pl.col("_resp")).alias("_eval_id"))
        .collect()
    )

    invalid_resps = sorted(
        set(norm_cells["_resp"].unique().to_list()) - {"responder", "non-responder"}
    )
    if invalid_resps:
        return Failure(f"[{inputs.accession}] Invalid clinical_response labels: {invalid_resps}")

    eval_meta = (
        norm_cells.lazy()
        .group_by("_eval_id")
        .agg(
            pl.col("_pt").first().alias("patient_id"),
            pl.col("_resp").first().alias("clinical_response"),
            pl.len().cast(pl.Float64).alias("library_size"),
        )
        .sort("_eval_id")
        .collect()
    )

    eval_ids = tuple(str(x) for x in eval_meta["_eval_id"].to_list())
    eval_pts = tuple(str(x) for x in eval_meta["patient_id"].to_list())
    n_evals = len(eval_ids)
    if n_evals < 4:
        return Failure(f"[{inputs.accession}] Requires >= 4 evaluations, got {n_evals}")

    eval_to_idx = {eid: idx for idx, eid in enumerate(eval_ids)}
    cell_eval_codes = np.array(
        [eval_to_idx[str(eid)] for eid in norm_cells["_eval_id"].to_list()],
        dtype=np.int64,
    )
    sample_ind = sp.csr_matrix(
        (np.ones(n_cells, dtype=np.float64), (np.arange(n_cells, dtype=np.int64), cell_eval_codes)),
        shape=(n_cells, n_evals),
    )

    nhoods_csc = sp.csc_matrix(inputs.nhoods_mat)
    count_matrix = np.asarray(nhoods_csc.T.dot(sample_ind).toarray(), dtype=np.float64)
    library_sizes = eval_meta["library_size"].to_numpy().astype(np.float64)
    response_labels = tuple(str(r) for r in eval_meta["clinical_response"].to_list())
    n_resp = sum(1 for r in response_labels if r == "responder")
    cell_types = tuple(str(ct) for ct in norm_cells["_ct"].to_list())

    return Success(
        CohortEvaluationDesign(
            accession=inputs.accession,
            eval_ids=eval_ids,
            eval_patient_ids=eval_pts,
            response_labels=response_labels,
            count_matrix=count_matrix,
            library_sizes=library_sizes,
            cell_eval_codes=cell_eval_codes,
            cell_types=cell_types,
            n_cells=n_cells,
            n_nhoods=n_nhoods,
            n_responders=n_resp,
            n_non_responders=n_evals - n_resp,
        )
    )


def map_cohort_nhoods_to_meta42_outcome_blind(
    inputs: CohortInputsBundle,
    concordance_3b_df: pl.DataFrame,
    hub_summary_df: pl.DataFrame,
) -> Result[NhoodMeta42Mapping, str]:
    """Maps every Milo neighborhood to (L_k, Consensus-Meta42) strictly WITHOUT clinical outcomes.

    Uses only unlabelled donor cell-count containment `sum_s min(C_{v,s}, N_{k,s}) / sum_s C_{v,s}`,
    across-donor arcsin-sqrt fraction correlation, ontology compartment compatibility (when
    informative Nhood_CellType labels exist), and 1-hop kNN neighborhood overlap smoothing.
    """
    acc = inputs.accession
    sub_3b = concordance_3b_df.filter(pl.col("cohort") == acc)
    if sub_3b.is_empty():
        return Failure(f"[{acc}] No rows found in 3b concordance table")

    lk_to_meta = {
        str(r["cell_state"]): str(r["meta_state"])
        for r in sub_3b.select(["cell_state", "meta_state"]).to_dicts()
    }
    lk_ref_milo = {
        str(r["cell_state"]): float(r["milo_mean_logfc"])
        for r in sub_3b.select(["cell_state", "milo_mean_logfc"]).to_dicts()
    }
    meta_to_label = {
        str(r["cell_state"]): str(r["state_label"])
        for r in hub_summary_df.select(["cell_state", "state_label"]).to_dicts()
    }

    pb_df = inputs.pseudobulk_fractions_df
    sids = sorted(str(x).strip() for x in pb_df["sample_id"].unique().to_list())
    lks = sorted(str(x).strip() for x in pb_df["cell_state"].unique().to_list())
    missing_lks = [lk for lk in lks if lk not in lk_to_meta]
    if missing_lks:
        return Failure(f"[{acc}] States in pseudobulk_fractions missing from 3b: {missing_lks[:5]}")

    sid_to_idx = {s: idx for idx, s in enumerate(sids)}
    lk_to_idx = {lk: idx for idx, lk in enumerate(lks)}
    n_sids = len(sids)
    n_lks = len(lks)

    cell_pts = [str(p).strip() for p in inputs.cells_df["patient_id"].to_list()]
    unmatched_pts = sorted({p for p in cell_pts if p not in sid_to_idx})
    if unmatched_pts:
        return Failure(
            f"[{acc}] Patient IDs in cell_level_scores missing from pseudobulk_fractions: {unmatched_pts[:5]}"
        )

    pt_codes = np.array([sid_to_idx[p] for p in cell_pts], dtype=np.int64)
    n_cells = len(cell_pts)
    nhoods_csc = (sp.csc_matrix(inputs.nhoods_mat) > 0).astype(np.float64)
    n_nhoods = int(nhoods_csc.shape[1])

    donor_ind = sp.csr_matrix(
        (np.ones(n_cells, dtype=np.float64), (np.arange(n_cells, dtype=np.int64), pt_codes)),
        shape=(n_cells, n_sids),
    )
    c_vs = np.asarray(nhoods_csc.T.dot(donor_ind).toarray(), dtype=np.float64)
    tot_s = np.maximum(np.bincount(pt_codes, minlength=n_sids).astype(np.float64), 1.0)
    f_vs = np.arcsin(np.sqrt(np.clip(c_vs / np.expand_dims(tot_s, axis=0), 0.0, 1.0)))

    n_ks = np.zeros((n_lks, n_sids), dtype=np.float64)
    tf_ks = np.zeros((n_lks, n_sids), dtype=np.float64)
    for sid_raw, lk_raw, n_c_val, tf_val in pb_df.select(
        ["sample_id", "cell_state", "n_cells", "true_fraction"]
    ).iter_rows():
        ki = lk_to_idx[str(lk_raw).strip()]
        si = sid_to_idx[str(sid_raw).strip()]
        n_ks[ki, si] = float(n_c_val)
        tf_ks[ki, si] = float(tf_val)
    f_ks = np.arcsin(np.sqrt(np.clip(tf_ks, 0.0, 1.0)))

    # 1. Standardized donor-profile correlation across all S donors (zero clinical outcome usage!)
    fv_c = f_vs - np.mean(f_vs, axis=1, keepdims=True)
    fv_norm = np.linalg.norm(fv_c, axis=1, keepdims=True)
    fv_unit = np.where(fv_norm > 1e-12, fv_c / np.maximum(fv_norm, 1e-12), 0.0)

    fk_c = f_ks - np.mean(f_ks, axis=1, keepdims=True)
    fk_norm = np.linalg.norm(fk_c, axis=1, keepdims=True)
    fk_unit = np.where(fk_norm > 1e-12, fk_c / np.maximum(fk_norm, 1e-12), 0.0)
    corr_vk = fv_unit @ fk_unit.T

    # 2. Exact integer cell containment ratio: sum_s min(C_{v,s}, N_{k,s}) / sum_s C_{v,s}
    nh_tot_cells = np.maximum(np.sum(c_vs, axis=1, keepdims=True), 1.0)
    contain_vk = (
        np.sum(np.minimum(np.expand_dims(c_vs, axis=1), np.expand_dims(n_ks, axis=0)), axis=2)
        / nh_tot_cells
    )

    raw_affinity = 0.65 * corr_vk + 0.35 * contain_vk

    # 3. 1-hop graph smoothing over overlapping Milo kNN neighborhoods
    adj = nhoods_csc.T.dot(nhoods_csc).tocsr()
    deg = np.maximum(np.asarray(adj.sum(axis=1)).ravel(), 1.0)
    smoothed_affinity = 0.70 * raw_affinity + 0.30 * (
        adj.dot(raw_affinity) / np.expand_dims(deg, axis=1)
    )

    # 4. Lineage compartment compatibility constraint when informative Nhood_CellType labels exist
    nhood_cts = (
        tuple(str(x) for x in inputs.da_results_df["Nhood_CellType"].to_list())
        if "Nhood_CellType" in inputs.da_results_df.columns
        else tuple("Unspecified" for _ in range(n_nhoods))
    )
    has_explicit_plasma = any(
        ("plasma" in ct.lower() and "plasmacytoid" not in ct.lower()) for ct in nhood_cts
    )
    lk_compartments = tuple(classify_meta42_compartment(lk_to_meta[lk]) for lk in lks)

    compat_mask = np.ones((n_nhoods, n_lks), dtype=bool)
    for v_idx, nh_ct in enumerate(nhood_cts):
        allowed = compatible_compartments_for_nhood_celltype(nh_ct, has_explicit_plasma)
        if allowed != _ALL_COMPARTMENTS:
            row_ok = np.array([comp in allowed for comp in lk_compartments], dtype=bool)
            if bool(np.any(row_ok)):
                compat_mask[v_idx, :] = row_ok

    masked_affinity = np.where(compat_mask, smoothed_affinity, smoothed_affinity - 0.50)
    best_k_idx = np.argmax(masked_affinity, axis=1)

    v_arange = np.arange(n_nhoods, dtype=np.int64)
    chosen_lks = tuple(lks[int(k)] for k in best_k_idx.tolist())
    chosen_metas = tuple(lk_to_meta[lk] for lk in chosen_lks)
    chosen_labels = tuple(meta_to_label.get(ms, ms) for ms in chosen_metas)
    chosen_comps = tuple(classify_meta42_compartment(ms).value for ms in chosen_metas)
    chosen_corr = corr_vk[v_arange, best_k_idx]
    chosen_contain = contain_vk[v_arange, best_k_idx]
    chosen_aff = smoothed_affinity[v_arange, best_k_idx]

    # Post-hoc diagnostic ONLY (never feeds back into mapping!): correlation with 3b milo_mean_logfc
    posthoc_r = 0.0
    if "logFC" in inputs.da_results_df.columns:
        da_lfc = inputs.da_results_df["logFC"].cast(pl.Float64).to_numpy()
        emp_lfc: list[float] = []
        ref_lfc: list[float] = []
        for ki, lk in enumerate(lks):
            m_ki = best_k_idx == ki
            if int(np.sum(m_ki)) >= 3:
                emp_lfc.append(float(np.mean(da_lfc[m_ki])))
                ref_lfc.append(lk_ref_milo[lk])
        if len(emp_lfc) >= 3 and float(np.std(emp_lfc)) > 1e-12 and float(np.std(ref_lfc)) > 1e-12:
            posthoc_r = float(np.corrcoef(emp_lfc, ref_lfc)[0, 1])

    has_stromal_caf = any(
        ms in STROMAL_CAF_BARRIER_META_STATES for ms in chosen_metas
    )
    annotations_df = pl.DataFrame(
        {
            "accession": [acc] * n_nhoods,
            "nhood_idx": v_arange,
            "nhood_cell_type": list(nhood_cts),
            "mapped_cell_state": list(chosen_lks),
            "meta_state": list(chosen_metas),
            "meta_state_label": list(chosen_labels),
            "lineage_compartment": list(chosen_comps),
            "donor_profile_corr": chosen_corr,
            "containment_ratio": chosen_contain,
            "composite_affinity": chosen_aff,
        }
    )

    return Success(
        NhoodMeta42Mapping(
            accession=acc,
            annotations_df=annotations_df,
            nhood_meta_states=chosen_metas,
            nhood_compartments=chosen_comps,
            has_stromal_caf=has_stromal_caf,
            n_unique_lk=len(set(chosen_lks)),
            n_total_lk=n_lks,
            n_unique_meta42=len(set(chosen_metas)),
            posthoc_milo_logfc_r=posthoc_r,
        )
    )


# =============================================================================
# 4. Pure Functional Core — Six-Tier Outcome-Blind Prefiltering Benchmark
# =============================================================================


def build_bulk_deconv_state_pair_weights(
    iatlas_all42_df: pl.DataFrame,
    base_floor: float = 0.20,
    slope: float = 1.60,
    z_offset: float = 1.25,
    max_cap: float = 8.0,
) -> dict[tuple[str, str], float]:
    """Constructs external state-pair prior weights from independent Pan-iAtlas |epistasis_z|."""
    weights: dict[tuple[str, str], float] = {}
    for sa, sb, ez in iatlas_all42_df.select(["state_a", "state_b", "epistasis_z"]).iter_rows():
        sa_s = str(sa).strip()
        sb_s = str(sb).strip()
        if sa_s != sb_s:
            excess = max(abs(float(ez)) - z_offset, 0.0)
            w_val = float(np.clip(base_floor + slope * excess, base_floor, max_cap))
            key = (sa_s, sb_s) if sa_s < sb_s else (sb_s, sa_s)
            weights[key] = max(weights.get(key, base_floor), w_val)
    return weights


def build_tier_pairwise_config(
    tier: AblationTier,
    config: BridgeConfig,
    n_responders: int,
    n_non_responders: int,
) -> PairwiseConfig:
    """Resolves the exact immutable PairwiseConfig for a given prefiltering ablation tier."""
    min_arm = min(n_responders, n_non_responders)
    eff_min_donors = min(config.min_donors_per_quadrant, float(max(2, min_arm)))
    eff_min_total = min(6.0, float(max(4, n_responders + n_non_responders - 1)))

    base_presence = PresenceConfig(mode=config.presence_mode)
    base_refine = RefinementConfig(
        enabled=config.max_refine_pairs > 0,
        max_refine_pairs=config.max_refine_pairs,
    )
    base_modules = PairModuleConfig(enabled=True)

    match tier:
        case AblationTier.T0_UNFILTERED:
            overlap_cfg = OverlapConfig(max_jaccard_overlap=0.0, scope_filter=PairScopeFilter.ALL)
            prefilter_cfg = PrefilterConfig(mode=PrefilterMode.NONE, apply_tarone_iut_pruning=False)
            apply_tarone = False
        case AblationTier.T1_TARONE_EXACT:
            overlap_cfg = OverlapConfig(max_jaccard_overlap=0.0, scope_filter=PairScopeFilter.ALL)
            prefilter_cfg = PrefilterConfig(
                mode=PrefilterMode.TARONE_ONLY,
                apply_tarone_iut_pruning=True,
                tarone_p_ceiling=config.tarone_p_ceiling,
            )
            apply_tarone = True
        case AblationTier.T2_MARGINAL_ENTROPY_TARONE:
            overlap_cfg = OverlapConfig(max_jaccard_overlap=0.0, scope_filter=PairScopeFilter.ALL)
            prefilter_cfg = PrefilterConfig(
                mode=PrefilterMode.MARGINAL_ENTROPY_VARIANCE,
                min_marginal_prevalence=0.10,
                max_marginal_prevalence=0.90,
                min_variance_quantile=0.25,
                apply_tarone_iut_pruning=True,
                tarone_p_ceiling=config.tarone_p_ceiling,
            )
            apply_tarone = True
        case AblationTier.T3_CROSS_META42_TARONE:
            overlap_cfg = OverlapConfig(
                max_jaccard_overlap=0.0,
                scope_filter=PairScopeFilter.CROSS_CELL_TYPE_ONLY,
            )
            prefilter_cfg = PrefilterConfig(
                mode=PrefilterMode.MARGINAL_ENTROPY_VARIANCE,
                min_marginal_prevalence=0.10,
                max_marginal_prevalence=0.90,
                min_variance_quantile=0.20,
                apply_tarone_iut_pruning=True,
                tarone_p_ceiling=config.tarone_p_ceiling,
            )
            apply_tarone = True
        case AblationTier.T4_SUPER_NHOODS_TARONE | AblationTier.T5_BULK_DECONV_PRIOR:
            overlap_cfg = OverlapConfig(
                max_jaccard_overlap=0.0,
                scope_filter=PairScopeFilter.CROSS_CELL_TYPE_ONLY,
            )
            prefilter_cfg = PrefilterConfig(
                mode=PrefilterMode.COMPOSITE_SUPER_TARONE,
                min_marginal_prevalence=0.08,
                max_marginal_prevalence=0.92,
                min_variance_quantile=0.15,
                super_nhood_min_jaccard=0.10,
                super_nhood_resolution=1.0,
                apply_tarone_iut_pruning=True,
                tarone_p_ceiling=config.tarone_p_ceiling,
            )
            apply_tarone = True
        case _ as unreachable:
            assert_never(unreachable)

    screen_cfg = ScreenConfig(
        chunk_size=config.chunk_size,
        fdr_threshold=config.fdr_threshold,
        min_donors_per_quadrant=eff_min_donors,
        min_effective_samples=eff_min_total,
        max_retained_candidates=50000,
        apply_tarone_iut_pruning=apply_tarone,
        tarone_p_ceiling=config.tarone_p_ceiling,
    )

    return PairwiseConfig(
        method=PairwiseMethod.HYBRID_ALL,
        overlap=overlap_cfg,
        presence=base_presence,
        prefilter=prefilter_cfg,
        screen=screen_cfg,
        refinement=base_refine,
        modules=base_modules,
    )


def _summarize_tier_result(
    accession: str,
    cohort_label: str,
    tier: AblationTier,
    res: PairwiseAnalysisResult,
    baseline_pairs: int,
    runtime_s: float,
) -> dict[str, object]:
    """Extracts standardized power, hypothesis compression, and discovery metrics for one tier."""
    min_p_q = (
        float(np.min(res.pairs_df["pvalue_iut"].to_numpy()))
        if not res.pairs_df.is_empty()
        else 1.0
    )
    min_p_c = (
        float(np.min(res.continuous_pairs_df["pvalue_iut"].to_numpy()))
        if not res.continuous_pairs_df.is_empty()
        else 1.0
    )
    min_raw_p = max(min(min_p_q, min_p_c), 1e-300)

    min_fdr_q = (
        float(np.min(res.pairs_df["fdr_iut"].to_numpy()))
        if not res.pairs_df.is_empty()
        else 1.0
    )
    min_fdr_c = (
        float(np.min(res.continuous_pairs_df["fdr_iut"].to_numpy()))
        if not res.continuous_pairs_df.is_empty()
        else 1.0
    )
    min_fdr = max(min(min_fdr_q, min_fdr_c), 1e-12)

    sig_q = res.pairs_df.filter(pl.col("is_significant"))
    sig_c = res.continuous_pairs_df.filter(pl.col("is_significant"))

    n_pp = int(sig_q.filter(pl.col("quadrant") == QuadrantState.BOTH_PRESENT.value).height)
    n_gated = int(
        sig_q.filter(
            pl.col("quadrant").is_in(
                [QuadrantState.FIRST_ONLY.value, QuadrantState.SECOND_ONLY.value]
            )
        ).height
    )
    n_mm = int(sig_q.filter(pl.col("quadrant") == QuadrantState.BOTH_MISSING.value).height)
    n_bal = int(
        sig_c.filter(pl.col("axis") == ContinuousAxis.LOG_RATIO_BALANCE.value).height
    )
    n_prod = int(
        sig_c.filter(pl.col("axis") == ContinuousAxis.CO_ABUNDANCE_PRODUCT.value).height
    )

    denom_base = float(max(baseline_pairs, 1))
    pair_reduction_pct = (
        100.0 * (1.0 - float(res.n_evaluated_pairs) / denom_base)
        if baseline_pairs > 0
        else 0.0
    )

    return {
        "accession": accession,
        "cohort_label": cohort_label,
        "tier": tier.value,
        "tier_label": ablation_tier_display_label(tier),
        "n_total_nhoods": res.n_total_nhoods,
        "n_active_units": res.n_active_nhoods,
        "n_evaluated_pairs": res.n_evaluated_pairs,
        "pair_reduction_pct": float(np.clip(pair_reduction_pct, -100.0, 100.0)),
        "min_raw_pvalue_iut": min_raw_p,
        "max_neg_log10_p_iut": float(-np.log10(min_raw_p)),
        "min_fdr_iut": min_fdr,
        "max_neg_log10_fdr_iut": float(-np.log10(min_fdr)),
        "n_sig_quadrant": res.n_significant_quadrant_pairs,
        "n_sig_pp": n_pp,
        "n_sig_gated": n_gated,
        "n_sig_mm": n_mm,
        "n_sig_continuous": res.n_significant_continuous_pairs,
        "n_sig_log_ratio": n_bal,
        "n_sig_co_abundance": n_prod,
        "n_sig_total": res.n_significant_quadrant_pairs + res.n_significant_continuous_pairs,
        "n_modules": res.n_modules,
        "runtime_seconds": runtime_s,
    }


def execute_cohort_tiers_pure(
    nhoods_mat: SparseMatrix,
    design: CohortEvaluationDesign,
    mapping: NhoodMeta42Mapping,
    bulk_pair_weights: Mapping[tuple[str, str], float],
    config: BridgeConfig,
    cached_ablation_df: Maybe[pl.DataFrame] = Nothing,
) -> Result[CohortBridgeExecution, str]:
    """Pure execution of all 6 prefiltering tiers for one cohort using its exact sparse nhoods_mat."""
    acc = design.accession
    cohort_label = get_single_cell_cohort_display_label(acc, "cancer_author_year")
    lib_tuple = tuple(float(x) for x in design.library_sizes.tolist())

    prior_mat_res = build_external_prior_weight_matrix(
        nhood_state_labels=mapping.nhood_meta_states,
        state_pair_weights=bulk_pair_weights,
        base_floor=0.20,
        max_cap=8.0,
    )
    match prior_mat_res:
        case Failure(err):
            return Failure(f"[{acc}] Prior weight construction failed: {err}")
        case Success(prior_mat):
            pass

    early_tier_vals = [
        AblationTier.T0_UNFILTERED.value,
        AblationTier.T1_TARONE_EXACT.value,
        AblationTier.T2_MARGINAL_ENTROPY_TARONE.value,
        AblationTier.T3_CROSS_META42_TARONE.value,
    ]
    cached_early_maybe: Maybe[pl.DataFrame] = Nothing
    match cached_ablation_df:
        case Some(df_cache):
            if {"accession", "tier", "n_evaluated_pairs"}.issubset(set(df_cache.columns)):
                sub_early = df_cache.filter(
                    (pl.col("accession") == acc) & pl.col("tier").is_in(early_tier_vals)
                )
                if sub_early.height == 4:
                    cached_early_maybe = Some(sub_early)
        case Maybe.empty:
            pass
        case _ as unreachable_cache:
            _ = unreachable_cache

    tiers_to_run: tuple[AblationTier, ...] = (
        (AblationTier.T4_SUPER_NHOODS_TARONE, AblationTier.T5_BULK_DECONV_PRIOR)
        if isinstance(cached_early_maybe, Some)
        else ALL_ABLATION_TIERS
    )

    tier_results: dict[AblationTier, tuple[PairwiseAnalysisResult, float]] = {}
    for tier in tiers_to_run:
        t0 = time.perf_counter()
        pw_cfg = build_tier_pairwise_config(
            tier=tier,
            config=config,
            n_responders=design.n_responders,
            n_non_responders=design.n_non_responders,
        )
        use_meta_labels = tier in (
            AblationTier.T3_CROSS_META42_TARONE,
            AblationTier.T4_SUPER_NHOODS_TARONE,
            AblationTier.T5_BULK_DECONV_PRIOR,
        )
        group_arg: Maybe[Sequence[str]] = (
            Some(mapping.nhood_meta_states) if use_meta_labels else Nothing
        )
        prior_arg: Maybe[FloatArray] = (
            Some(prior_mat) if tier == AblationTier.T5_BULK_DECONV_PRIOR else Nothing
        )

        res_wrapped = analyze_nhood_pairs(
            nhoods_mat=nhoods_mat,
            count_matrix=design.count_matrix,
            response_labels=design.response_labels,
            sample_library_sizes=Some(lib_tuple),
            cell_types=Some(design.cell_types),
            nhood_group_labels=group_arg,
            cell_sample_codes=Some(design.cell_eval_codes),
            prior_pair_weights=prior_arg,
            config=pw_cfg,
        )
        match res_wrapped:
            case Failure(err):
                return Failure(f"[{acc}][{tier.value}] Pairwise analysis failed: {err}")
            case Success(pw_res):
                tier_results[tier] = (pw_res, time.perf_counter() - t0)

    match cached_early_maybe:
        case Some(sub_early_df):
            t0_row = sub_early_df.filter(pl.col("tier") == AblationTier.T0_UNFILTERED.value)
            baseline_pairs = int(t0_row["n_evaluated_pairs"][0])
            late_rows_df = pl.DataFrame(
                [
                    _summarize_tier_result(
                        accession=acc,
                        cohort_label=cohort_label,
                        tier=tier,
                        res=tier_results[tier][0],
                        baseline_pairs=baseline_pairs,
                        runtime_s=tier_results[tier][1],
                    )
                    for tier in tiers_to_run
                ]
            )
            aligned_early = sub_early_df.select(late_rows_df.columns)
            tier_rank = {t.value: idx for idx, t in enumerate(ALL_ABLATION_TIERS)}
            ablation_rows_df = (
                pl.concat([aligned_early, late_rows_df])
                .with_columns(
                    pl.col("tier").replace_strict(tier_rank, default=99).alias("_ord")
                )
                .sort("_ord")
                .drop("_ord")
            )
        case Maybe.empty:
            baseline_pairs = tier_results[AblationTier.T0_UNFILTERED][0].n_evaluated_pairs
            ablation_rows_df = pl.DataFrame(
                [
                    _summarize_tier_result(
                        accession=acc,
                        cohort_label=cohort_label,
                        tier=tier,
                        res=tier_results[tier][0],
                        baseline_pairs=baseline_pairs,
                        runtime_s=tier_results[tier][1],
                    )
                    for tier in ALL_ABLATION_TIERS
                ]
            )
        case _ as unreachable_ab:
            _ = unreachable_ab
            ablation_rows_df = pl.DataFrame()

    return Success(
        CohortBridgeExecution(
            accession=acc,
            cohort_label=cohort_label,
            design=design,
            mapping=mapping,
            ablation_rows_df=ablation_rows_df,
            primary_t4_result=tier_results[AblationTier.T4_SUPER_NHOODS_TARONE][0],
            prior_t5_result=tier_results[AblationTier.T5_BULK_DECONV_PRIOR][0],
        )
    )


# =============================================================================
# 5. Pure Functional Core — Deconvolution Concordance, Mantel & Stromal Gating
# =============================================================================


def _p_to_two_sided_z(pval: float) -> float:
    """Converts a two-sided p-value in (0, 1] to a non-negative standard normal Z magnitude."""
    p_safe = float(np.clip(pval, 1e-300, 1.0))
    return 0.0 if p_safe >= 1.0 else float(norm.isf(0.5 * p_safe))


def extract_cohort_meta42_pair_metrics(
    exec_bundle: CohortBridgeExecution,
) -> pl.DataFrame:
    """Summarizes T4 & T5 single-cell pairwise synergy-vs-gating Z-scores per unordered Meta42 pair.

    Computes one balanced net interaction Z per super-neighborhood pair `(nhood_u, nhood_v)` before
    aggregating onto `(state_a, state_b)`, matching bulk `consensus_metastate_all42_interactions.parquet`:
    - Positive (+): Co-activation synergy (`++` Responder, `--` Non-Responder immune desert,
      Responder `CO_ABUNDANCE_PRODUCT`)
    - Negative (-): Antagonistic gating or co-suppression (`+-`/`-+` discordant gating,
      `++` Non-Responder / `--` Responder co-suppression scaled by `0.65`, discordant `LOG_RATIO_BALANCE`)
    """
    t4 = exec_bundle.primary_t4_result
    t5 = exec_bundle.prior_t5_result
    n_evals = int(exec_bundle.design.n_responders + exec_bundle.design.n_non_responders)

    uv_meta: dict[tuple[int, int], tuple[str, str]] = {}
    uv_syn: dict[tuple[int, int], float] = {}
    uv_disc: dict[tuple[int, int], float] = {}
    uv_cosupp: dict[tuple[int, int], float] = {}
    uv_bal: dict[tuple[int, int], float] = {}
    pair_sig_t4: dict[tuple[str, str], int] = {}
    pair_sig_t5: dict[tuple[str, str], int] = {}

    if not t4.pairs_df.is_empty():
        for row in t4.pairs_df.select(
            [
                "nhood_u",
                "nhood_v",
                "cell_type_u",
                "cell_type_v",
                "quadrant",
                "direction",
                "pvalue_iut",
                "log_or_interaction",
                "is_significant",
            ]
        ).to_dicts():
            sa = str(row["cell_type_u"])
            sb = str(row["cell_type_v"])
            if sa == sb or not sa.startswith("M") or not sb.startswith("M"):
                continue
            key = (sa, sb) if sa < sb else (sb, sa)
            if bool(row["is_significant"]):
                pair_sig_t4[key] = pair_sig_t4.get(key, 0) + 1
            p_iut = float(row["pvalue_iut"])
            if p_iut >= 0.9999:
                continue
            u_id = int(row["nhood_u"])
            v_id = int(row["nhood_v"])
            uv = (u_id, v_id) if u_id < v_id else (v_id, u_id)
            uv_meta[uv] = key

            z_iut = _p_to_two_sided_z(p_iut)
            q_str = str(row["quadrant"])
            d_str = str(row["direction"])
            lor_int = float(row["log_or_interaction"])

            if q_str == QuadrantState.BOTH_PRESENT.value:
                if d_str == EnrichmentDirection.RESPONDER.value:
                    uv_syn[uv] = max(uv_syn.get(uv, 0.0), z_iut)
                else:
                    uv_cosupp[uv] = max(uv_cosupp.get(uv, 0.0), z_iut)
            elif q_str == QuadrantState.BOTH_MISSING.value:
                if d_str == EnrichmentDirection.NON_RESPONDER.value:
                    uv_syn[uv] = max(uv_syn.get(uv, 0.0), z_iut)
                else:
                    uv_cosupp[uv] = max(uv_cosupp.get(uv, 0.0), z_iut)
            else:
                disc_z = z_iut if lor_int <= 0.25 else 0.50 * z_iut
                uv_disc[uv] = max(uv_disc.get(uv, 0.0), disc_z)

    if not t4.continuous_pairs_df.is_empty():
        for row in t4.continuous_pairs_df.select(
            [
                "nhood_u",
                "nhood_v",
                "cell_type_u",
                "cell_type_v",
                "axis",
                "direction",
                "pvalue_iut",
                "is_significant",
            ]
        ).to_dicts():
            sa = str(row["cell_type_u"])
            sb = str(row["cell_type_v"])
            if sa == sb or not sa.startswith("M") or not sb.startswith("M"):
                continue
            key = (sa, sb) if sa < sb else (sb, sa)
            if bool(row["is_significant"]):
                pair_sig_t4[key] = pair_sig_t4.get(key, 0) + 1
            p_iut = float(row["pvalue_iut"])
            if p_iut >= 0.9999:
                continue
            u_id = int(row["nhood_u"])
            v_id = int(row["nhood_v"])
            uv = (u_id, v_id) if u_id < v_id else (v_id, u_id)
            uv_meta[uv] = key

            z_iut = _p_to_two_sided_z(p_iut)
            axis_s = str(row["axis"])
            d_str = str(row["direction"])
            if axis_s == ContinuousAxis.CO_ABUNDANCE_PRODUCT.value:
                if d_str == EnrichmentDirection.RESPONDER.value:
                    uv_syn[uv] = max(uv_syn.get(uv, 0.0), z_iut)
                else:
                    uv_cosupp[uv] = max(uv_cosupp.get(uv, 0.0), z_iut)
            elif axis_s == ContinuousAxis.LOG_RATIO_BALANCE.value:
                uv_bal[uv] = max(uv_bal.get(uv, 0.0), z_iut)

    pair_synergy: dict[tuple[str, str], list[float]] = {}
    pair_gating: dict[tuple[str, str], list[float]] = {}
    pair_epistasis: dict[tuple[str, str], list[float]] = {}
    for uv, key in uv_meta.items():
        syn_u = uv_syn.get(uv, 0.0)
        ant_u = max(
            uv_disc.get(uv, 0.0),
            0.65 * uv_cosupp.get(uv, 0.0),
            0.85 * uv_bal.get(uv, 0.0),
        )
        pair_synergy.setdefault(key, []).append(syn_u)
        pair_gating.setdefault(key, []).append(ant_u)
        pair_epistasis.setdefault(key, []).append(syn_u - ant_u)

    for df_t5 in (t5.pairs_df, t5.continuous_pairs_df):
        if not df_t5.is_empty():
            for row in df_t5.filter(pl.col("is_significant")).select(
                ["cell_type_u", "cell_type_v"]
            ).to_dicts():
                sa = str(row["cell_type_u"])
                sb = str(row["cell_type_v"])
                if sa != sb and sa.startswith("M") and sb.startswith("M"):
                    key = (sa, sb) if sa < sb else (sb, sa)
                    pair_sig_t5[key] = pair_sig_t5.get(key, 0) + 1

    all_keys = sorted(set(pair_epistasis.keys()) | set(pair_sig_t5.keys()))
    if not all_keys:
        return pl.DataFrame(
            schema={
                "accession": pl.String,
                "n_evals": pl.Int64,
                "has_stromal_caf": pl.Boolean,
                "state_a": pl.String,
                "state_b": pl.String,
                "cohort_sc_epistasis_z": pl.Float64,
                "cohort_max_synergy_z": pl.Float64,
                "cohort_max_gating_z": pl.Float64,
                "n_sig_t4": pl.Int64,
                "n_sig_t5": pl.Int64,
            }
        )

    out_rows: list[dict[str, object]] = []
    for key in all_keys:
        vals = pair_epistasis.get(key, [0.0])
        arr = np.asarray(vals, dtype=np.float64)
        top_order = np.argsort(-np.abs(arr), kind="mergesort")
        top_k = arr[top_order[: min(3, int(arr.shape[0]))]]
        cohort_z = 0.65 * float(np.mean(top_k)) + 0.35 * float(np.mean(arr))
        syn_vals = pair_synergy.get(key, [0.0])
        gat_vals = pair_gating.get(key, [0.0])
        out_rows.append(
            {
                "accession": exec_bundle.accession,
                "n_evals": n_evals,
                "has_stromal_caf": exec_bundle.mapping.has_stromal_caf,
                "state_a": key[0],
                "state_b": key[1],
                "cohort_sc_epistasis_z": cohort_z,
                "cohort_max_synergy_z": float(np.max(syn_vals)),
                "cohort_max_gating_z": float(np.max(gat_vals)),
                "n_sig_t4": pair_sig_t4.get(key, 0),
                "n_sig_t5": pair_sig_t5.get(key, 0),
            }
        )
    return pl.DataFrame(out_rows)


def run_mantel_node_permutation_test(
    states: Sequence[str],
    pair_to_x: Mapping[tuple[str, str], float],
    pair_to_y: Mapping[tuple[str, str], float],
    panel_name: str,
    n_perms: int = 5000,
    seed: int = 2026,
    strong_z_cutoff: float = 2.5,
) -> MantelTestResult:
    """Exact node-permutation Mantel test on symmetric off-diagonal interaction matrices.

    Simultaneously permutes row and column labels (`pi(a), pi(b)`) to preserve row/column hub
    degree structure under H_0.
    """
    state_list = sorted(set(states))
    m = len(state_list)
    if m < 3:
        return MantelTestResult(
            panel_name=panel_name,
            n_states=m,
            n_pairs=0,
            pearson_r=0.0,
            spearman_rho=0.0,
            pvalue_two_sided=1.0,
            pvalue_one_sided_pos=1.0,
            directional_concordance_strong=0.0,
            n_strong_pairs=0,
        )

    mat_x = np.zeros((m, m), dtype=np.float64)
    mat_y = np.zeros((m, m), dtype=np.float64)
    present_mask = np.zeros((m, m), dtype=bool)

    for i in range(m):
        for j in range(i + 1, m):
            key = (state_list[i], state_list[j])
            if key in pair_to_x and key in pair_to_y:
                vx = float(pair_to_x[key])
                vy = float(pair_to_y[key])
                mat_x[i, j] = mat_x[j, i] = vx
                mat_y[i, j] = mat_y[j, i] = vy
                present_mask[i, j] = present_mask[j, i] = True

    iu, ju = np.triu_indices(m, k=1)
    valid = present_mask[iu, ju]
    n_pairs = int(np.sum(valid))
    if n_pairs < 3:
        return MantelTestResult(
            panel_name=panel_name,
            n_states=m,
            n_pairs=n_pairs,
            pearson_r=0.0,
            spearman_rho=0.0,
            pvalue_two_sided=1.0,
            pvalue_one_sided_pos=1.0,
            directional_concordance_strong=0.0,
            n_strong_pairs=0,
        )

    x_obs = mat_x[iu[valid], ju[valid]]
    y_obs = mat_y[iu[valid], ju[valid]]
    sx = float(np.std(x_obs))
    sy = float(np.std(y_obs))
    if sx < 1e-12 or sy < 1e-12:
        return MantelTestResult(
            panel_name=panel_name,
            n_states=m,
            n_pairs=n_pairs,
            pearson_r=0.0,
            spearman_rho=0.0,
            pvalue_two_sided=1.0,
            pvalue_one_sided_pos=1.0,
            directional_concordance_strong=0.0,
            n_strong_pairs=0,
        )

    r_obs = float(np.corrcoef(x_obs, y_obs)[0, 1])
    rho_res = spearmanr(x_obs, y_obs)
    rho_obs = float(rho_res.statistic) if not np.isnan(float(rho_res.statistic)) else 0.0

    strong_mask = (np.abs(x_obs) >= strong_z_cutoff) & (np.abs(y_obs) > 1e-9)
    n_strong = int(np.sum(strong_mask))
    dir_conc = (
        float(np.mean(np.sign(x_obs[strong_mask]) == np.sign(y_obs[strong_mask])))
        if n_strong > 0
        else 0.0
    )

    rng = np.random.default_rng(seed)
    exceed_two = 0
    exceed_pos = 0
    ui_valid = iu[valid]
    uj_valid = ju[valid]

    for _ in range(n_perms):
        perm = rng.permutation(m)
        y_perm = mat_y[perm[ui_valid], perm[uj_valid]]
        m_perm = present_mask[perm[ui_valid], perm[uj_valid]]
        if int(np.sum(m_perm)) >= 3:
            xp = x_obs[m_perm]
            yp = y_perm[m_perm]
            sxp = float(np.std(xp))
            syp = float(np.std(yp))
            r_p = float(np.corrcoef(xp, yp)[0, 1]) if (sxp > 1e-12 and syp > 1e-12) else 0.0
        else:
            r_p = 0.0
        if abs(r_p) >= abs(r_obs) - 1e-12:
            exceed_two += 1
        if r_p >= r_obs - 1e-12:
            exceed_pos += 1

    p_two = float(1 + exceed_two) / float(1 + n_perms)
    p_pos = float(1 + exceed_pos) / float(1 + n_perms)

    return MantelTestResult(
        panel_name=panel_name,
        n_states=m,
        n_pairs=n_pairs,
        pearson_r=r_obs,
        spearman_rho=rho_obs,
        pvalue_two_sided=p_two,
        pvalue_one_sided_pos=p_pos,
        directional_concordance_strong=dir_conc,
        n_strong_pairs=n_strong,
    )


def aggregate_cross_modality_concordance(
    cohort_execs: Sequence[CohortBridgeExecution],
    global_deconv: GlobalDeconvInputs,
    n_mantel_perms: int,
    mantel_seed: int,
) -> tuple[
    pl.DataFrame,
    pl.DataFrame,
    MantelTestResult,
    MantelTestResult,
    MantelTestResult,
]:
    """Joins single-cell pairwise Milopy Meta42 scores with Pan-iAtlas (N=1,015) bulk interactions."""
    per_cohort_frames = [extract_cohort_meta42_pair_metrics(ex) for ex in cohort_execs]
    non_empty = [f for f in per_cohort_frames if not f.is_empty()]
    stacked_sc = (
        pl.concat(non_empty)
        if non_empty
        else pl.DataFrame(
            schema={
                "accession": pl.String,
                "n_evals": pl.Int64,
                "has_stromal_caf": pl.Boolean,
                "state_a": pl.String,
                "state_b": pl.String,
                "cohort_sc_epistasis_z": pl.Float64,
                "cohort_max_synergy_z": pl.Float64,
                "cohort_max_gating_z": pl.Float64,
                "n_sig_t4": pl.Int64,
                "n_sig_t5": pl.Int64,
            }
        )
    )

    sc_all9 = (
        stacked_sc.lazy()
        .group_by(["state_a", "state_b"])
        .agg(
            pl.len().cast(pl.Int64).alias("n_cohorts_observed"),
            (
                (pl.col("cohort_sc_epistasis_z") * pl.col("n_evals").cast(pl.Float64).sqrt()).sum()
                / pl.col("n_evals").cast(pl.Float64).sum().sqrt()
            ).alias("sc_epistasis_z_all9"),
            pl.col("cohort_max_synergy_z").max().alias("sc_max_synergy_z"),
            pl.col("cohort_max_gating_z").max().alias("sc_max_gating_z"),
            pl.col("n_sig_t4").sum().cast(pl.Int64).alias("total_sig_pairs_t4"),
            pl.col("n_sig_t5").sum().cast(pl.Int64).alias("total_sig_pairs_t5"),
        )
        .collect()
    )

    sc_wt = (
        stacked_sc.lazy()
        .filter(pl.col("has_stromal_caf"))
        .group_by(["state_a", "state_b"])
        .agg(
            pl.len().cast(pl.Int64).alias("n_whole_tumor_cohorts"),
            (
                (pl.col("cohort_sc_epistasis_z") * pl.col("n_evals").cast(pl.Float64).sqrt()).sum()
                / pl.col("n_evals").cast(pl.Float64).sum().sqrt()
            ).alias("sc_epistasis_z_whole_tumor"),
        )
        .collect()
    )

    bulk_upper = (
        global_deconv.iatlas_all42_interactions_df.lazy()
        .filter(pl.col("state_a") < pl.col("state_b"))
        .select(
            [
                "state_a",
                "state_b",
                "state_a_label",
                "state_b_label",
                "state_a_beta",
                "state_b_beta",
                "interaction_type",
                "friedman_h",
                pl.col("epistasis_z").alias("bulk_epistasis_z"),
                pl.col("signed_interaction").alias("bulk_signed_interaction"),
            ]
        )
        .collect()
    )

    focused_set = set(FOCUSED_8_META_STATES)
    joined_df = (
        bulk_upper.lazy()
        .join(sc_all9.lazy(), on=["state_a", "state_b"], how="left")
        .join(sc_wt.lazy(), on=["state_a", "state_b"], how="left")
        .with_columns(
            pl.col("n_cohorts_observed").fill_null(0),
            pl.col("n_whole_tumor_cohorts").fill_null(0),
            pl.col("sc_epistasis_z_all9").fill_null(0.0),
            pl.col("sc_epistasis_z_whole_tumor").fill_null(0.0),
            pl.col("sc_max_synergy_z").fill_null(0.0),
            pl.col("sc_max_gating_z").fill_null(0.0),
            pl.col("total_sig_pairs_t4").fill_null(0),
            pl.col("total_sig_pairs_t5").fill_null(0),
            (
                pl.col("state_a").is_in(list(focused_set))
                & pl.col("state_b").is_in(list(focused_set))
            ).alias("in_focused_8x8"),
        )
        .sort(["bulk_epistasis_z"], descending=True)
        .collect()
    )

    obs_rows = joined_df.filter(pl.col("n_cohorts_observed") > 0)
    pair_bulk: dict[tuple[str, str], float] = {}
    pair_sc_all9: dict[tuple[str, str], float] = {}
    pair_sc_wt: dict[tuple[str, str], float] = {}
    all_obs_states: set[str] = set()
    wt_obs_states: set[str] = set()

    for r in obs_rows.to_dicts():
        sa = str(r["state_a"])
        sb = str(r["state_b"])
        key = (sa, sb)
        pair_bulk[key] = float(r["bulk_epistasis_z"])
        pair_sc_all9[key] = float(r["sc_epistasis_z_all9"])
        all_obs_states.add(sa)
        all_obs_states.add(sb)
        if int(r["n_whole_tumor_cohorts"]) > 0:
            pair_sc_wt[key] = float(r["sc_epistasis_z_whole_tumor"])
            wt_obs_states.add(sa)
            wt_obs_states.add(sb)

    mantel_all42 = run_mantel_node_permutation_test(
        states=sorted(all_obs_states),
        pair_to_x=pair_bulk,
        pair_to_y=pair_sc_all9,
        panel_name="Consensus-Meta42 Observed Cross-Modality Matrix (All 9 Cohorts)",
        n_perms=n_mantel_perms,
        seed=mantel_seed,
    )
    mantel_wt = run_mantel_node_permutation_test(
        states=sorted(wt_obs_states),
        pair_to_x=pair_bulk,
        pair_to_y=pair_sc_wt,
        panel_name="Whole-Tumor Stroma-Complete Cohorts Matrix",
        n_perms=n_mantel_perms,
        seed=mantel_seed + 17,
    )
    mantel_f8 = run_mantel_node_permutation_test(
        states=FOCUSED_8_META_STATES,
        pair_to_x=pair_bulk,
        pair_to_y=pair_sc_wt if len(pair_sc_wt) >= 10 else pair_sc_all9,
        panel_name="Focused 8x8 Canonical Synergy/Gating Submatrix",
        n_perms=n_mantel_perms,
        seed=mantel_seed + 31,
        strong_z_cutoff=2.0,
    )

    # Node-level hub centrality concordance across Meta42 states
    end_a = obs_rows.select(
        pl.col("state_a").alias("cell_state"),
        pl.col("sc_epistasis_z_all9").abs().alias("abs_sc_z"),
        pl.col("sc_max_synergy_z").alias("syn_z"),
        pl.col("sc_max_gating_z").alias("gat_z"),
        pl.col("total_sig_pairs_t4").alias("sig_t4"),
        pl.col("total_sig_pairs_t5").alias("sig_t5"),
    )
    end_b = obs_rows.select(
        pl.col("state_b").alias("cell_state"),
        pl.col("sc_epistasis_z_all9").abs().alias("abs_sc_z"),
        pl.col("sc_max_synergy_z").alias("syn_z"),
        pl.col("sc_max_gating_z").alias("gat_z"),
        pl.col("total_sig_pairs_t4").alias("sig_t4"),
        pl.col("total_sig_pairs_t5").alias("sig_t5"),
    )
    sc_hub = (
        pl.concat([end_a, end_b])
        .lazy()
        .group_by("cell_state")
        .agg(
            pl.len().cast(pl.Int64).alias("sc_n_partner_states"),
            pl.col("abs_sc_z").mean().alias("sc_mean_abs_pair_z"),
            pl.col("syn_z").max().alias("sc_max_synergy_z"),
            pl.col("gat_z").max().alias("sc_max_gating_z"),
            pl.col("sig_t4").sum().cast(pl.Int64).alias("sc_total_sig_t4"),
            pl.col("sig_t5").sum().cast(pl.Int64).alias("sc_total_sig_t5"),
        )
        .collect()
    )

    hub_concordance_df = (
        global_deconv.iatlas_hub_summary_df.lazy()
        .join(sc_hub.lazy(), on="cell_state", how="left")
        .with_columns(
            pl.col("sc_n_partner_states").fill_null(0),
            pl.col("sc_mean_abs_pair_z").fill_null(0.0),
            pl.col("sc_max_synergy_z").fill_null(0.0),
            pl.col("sc_max_gating_z").fill_null(0.0),
            pl.col("sc_total_sig_t4").fill_null(0),
            pl.col("sc_total_sig_t5").fill_null(0),
        )
        .sort("mean_abs_epistasis_z", descending=True)
        .collect()
    )

    return joined_df, hub_concordance_df, mantel_all42, mantel_wt, mantel_f8


def compute_within_cohort_pseudobulk_concordance(
    cohort_execs: Sequence[CohortBridgeExecution],
    cohort_inputs_map: Mapping[str, CohortInputsBundle],
    concordance_3b_df: pl.DataFrame,
) -> pl.DataFrame:
    """Computes matched within-cohort concordance between single-cell pairwise Milo and 3a InstaPrism.

    Inside each cohort (J >= 10 donors), aggregates `inferred_fraction` from `pseudobulk_fractions.parquet`
    onto Meta42 states per donor, computes matched donor-level co-abundance vs. log-ratio balance
    association with clinical response, and evaluates correlation against single-cell T4 pair scores.
    """
    summary_rows: list[dict[str, object]] = []

    for ex in cohort_execs:
        acc = ex.accession
        sc_pairs_df = extract_cohort_meta42_pair_metrics(ex)
        if sc_pairs_df.height < 5:
            continue

        pb = cohort_inputs_map[acc].pseudobulk_fractions_df
        sub_3b = concordance_3b_df.filter(pl.col("cohort") == acc)
        lk_to_meta = {
            str(r["cell_state"]): str(r["meta_state"])
            for r in sub_3b.select(["cell_state", "meta_state"]).to_dicts()
        }

        pb_meta = (
            pb.lazy()
            .with_columns(
                pl.col("cell_state").replace_strict(lk_to_meta, default="OTHER").alias("meta_state"),
                pl.col("response").cast(pl.String).str.strip_chars().str.to_lowercase().alias("_resp"),
            )
            .filter(pl.col("_resp").is_in(["responder", "non-responder"]))
            .group_by(["sample_id", "meta_state"])
            .agg(
                pl.col("inferred_fraction").sum().alias("meta_inferred_frac"),
                pl.col("_resp").first().alias("response_norm"),
            )
            .collect()
        )

        sids = sorted(pb_meta["sample_id"].unique().to_list())
        n_donors = len(sids)
        if n_donors < 8:
            continue

        metas = sorted(pb_meta["meta_state"].unique().to_list())
        sid_to_i = {s: idx for idx, s in enumerate(sids)}
        meta_to_i = {m: idx for idx, m in enumerate(metas)}

        frac_mat = np.zeros((len(metas), n_donors), dtype=np.float64)
        y_vec = np.zeros(n_donors, dtype=np.float64)
        for sid, ms, val, resp_s in pb_meta.select(
            ["sample_id", "meta_state", "meta_inferred_frac", "response_norm"]
        ).iter_rows():
            mi = meta_to_i[str(ms)]
            si = sid_to_i[str(sid)]
            frac_mat[mi, si] = float(val)
            y_vec[si] = 1.0 if str(resp_s) == "responder" else 0.0

        n_r = int(np.sum(y_vec > 0.5))
        n_nr = n_donors - n_r
        if n_r < 2 or n_nr < 2:
            continue

        # Standardized arcsin-sqrt inferred fractions & log-fractions across the exact same donors
        asin_m = np.arcsin(np.sqrt(np.clip(frac_mat, 0.0, 1.0)))
        asin_c = asin_m - np.mean(asin_m, axis=1, keepdims=True)
        asin_sd = np.maximum(np.std(asin_c, axis=1, keepdims=True), 1e-12)
        z_mat = asin_c / asin_sd
        log_f = np.log2(frac_mat + 1e-4)

        y_c = y_vec - float(np.mean(y_vec))
        y_unit = y_c / max(float(np.linalg.norm(y_c)), 1e-12)
        sqrt_df = float(np.sqrt(max(n_donors - 2, 1)))

        sc_vals: list[float] = []
        pb_vals: list[float] = []
        for row in sc_pairs_df.to_dicts():
            sa = str(row["state_a"])
            sb = str(row["state_b"])
            if sa not in meta_to_i or sb not in meta_to_i:
                continue
            ia = meta_to_i[sa]
            ib = meta_to_i[sb]

            # Pseudobulk co-abundance product t-stat minus discordant log-ratio balance magnitude
            prod_v = z_mat[ia, :] * z_mat[ib, :]
            prod_c = prod_v - float(np.mean(prod_v))
            p_norm = float(np.linalg.norm(prod_c))
            r_prod = float(np.dot(prod_c / p_norm, y_unit)) if p_norm > 1e-12 else 0.0
            z_prod = r_prod * sqrt_df / float(np.sqrt(max(1.0 - r_prod * r_prod, 1e-6)))

            bal_v = log_f[ia, :] - log_f[ib, :]
            bal_c = bal_v - float(np.mean(bal_v))
            b_norm = float(np.linalg.norm(bal_c))
            r_bal = float(np.dot(bal_c / b_norm, y_unit)) if b_norm > 1e-12 else 0.0
            z_bal = abs(r_bal) * sqrt_df / float(np.sqrt(max(1.0 - r_bal * r_bal, 1e-6)))

            pb_score = 0.70 * z_prod - 0.30 * z_bal
            sc_vals.append(float(row["cohort_sc_epistasis_z"]))
            pb_vals.append(pb_score)

        if len(sc_vals) >= 5:
            x_arr = np.asarray(sc_vals, dtype=np.float64)
            y_arr = np.asarray(pb_vals, dtype=np.float64)
            r_p = (
                float(np.corrcoef(x_arr, y_arr)[0, 1])
                if (float(np.std(x_arr)) > 1e-12 and float(np.std(y_arr)) > 1e-12)
                else 0.0
            )
            sp_res = spearmanr(x_arr, y_arr)
            rho = float(sp_res.statistic) if not np.isnan(float(sp_res.statistic)) else 0.0
            pval_sp = float(sp_res.pvalue) if not np.isnan(float(sp_res.pvalue)) else 1.0
            strong = (np.abs(x_arr) >= 1.5) & (np.abs(y_arr) > 0.25)
            sign_agree = (
                float(np.mean(np.sign(x_arr[strong]) == np.sign(y_arr[strong])))
                if int(np.sum(strong)) > 0
                else float(np.mean(np.sign(x_arr) == np.sign(y_arr)))
            )
            summary_rows.append(
                {
                    "accession": acc,
                    "cohort_label": ex.cohort_label,
                    "n_donors": n_donors,
                    "n_shared_meta42_pairs": len(sc_vals),
                    "pearson_r": r_p,
                    "spearman_rho": rho,
                    "pvalue_spearman": pval_sp,
                    "sign_concordance": sign_agree,
                }
            )

    return (
        pl.DataFrame(summary_rows).sort("spearman_rho", descending=True)
        if summary_rows
        else pl.DataFrame(
            schema={
                "accession": pl.String,
                "cohort_label": pl.String,
                "n_donors": pl.Int64,
                "n_shared_meta42_pairs": pl.Int64,
                "pearson_r": pl.Float64,
                "spearman_rho": pl.Float64,
                "pvalue_spearman": pl.Float64,
                "sign_concordance": pl.Float64,
            }
        )
    )


def wilson_95_ci(n_success: int, n_total: int) -> tuple[float, float]:
    """Exact Wilson score 95% binomial confidence interval without continuity overshoot."""
    if n_total <= 0:
        return (0.0, 0.0)
    z = 1.959963984540054
    z2 = z * z
    n_f = float(n_total)
    p_hat = float(n_success) / n_f
    denom = 1.0 + z2 / n_f
    center = (p_hat + z2 / (2.0 * n_f)) / denom
    half = (z / denom) * float(np.sqrt((p_hat * (1.0 - p_hat)) / n_f + z2 / (4.0 * n_f * n_f)))
    return (float(np.clip(center - half, 0.0, 1.0)), float(np.clip(center + half, 0.0, 1.0)))


def compute_singlecell_mycaf_gating_regimes(
    cohort_execs: Sequence[CohortBridgeExecution],
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Tests the Pan-iAtlas M17 ISLR+COL5A2+ myCAF stromal exclusion gate directly in single-cell biopsies.

    In whole-tumor cohorts harboring stromal CAF neighborhoods (`has_stromal_caf == True`),
    computes each biopsy's Lymphoid/TLS effector fraction (`M04/M05/M13/M14/M15/M24/M29/M30/M34/M36/M37/M41`)
    and `M17/M23` myCAF/CAF stromal barrier fraction, stratifies within each cohort at the cohort
    median (outcome-blind marginal median!), and computes pooled empirical responder rates + Wilson
    95% CIs + per-cohort log-odds gating contrasts.
    """
    biopsy_records: list[dict[str, object]] = []
    forest_rows: list[dict[str, object]] = []

    for ex in cohort_execs:
        states = ex.mapping.nhood_meta_states
        eff_mask = np.array([ms in LYMPHOID_EFFECTOR_META_STATES for ms in states], dtype=bool)
        caf_mask = np.array([ms in STROMAL_CAF_BARRIER_META_STATES for ms in states], dtype=bool)
        bar_mask = np.array(
            [ms in STROMAL_MYELOID_BARRIER_META_STATES for ms in states],
            dtype=bool,
        )

        # Require unsorted / stroma-bearing cohort with both Lymphoid and M17/M23 CAF neighborhoods
        if not (bool(np.any(eff_mask)) and bool(np.any(caf_mask))):
            continue

        counts = ex.design.count_matrix
        libs = np.maximum(ex.design.library_sizes, 1.0)
        eff_score = np.sum(counts[eff_mask, :], axis=0) / libs
        caf_score = np.sum(counts[caf_mask, :], axis=0) / libs
        bar_score = np.sum(counts[bar_mask, :], axis=0) / libs

        # Outcome-blind marginal median split within cohort
        eff_med = float(np.median(eff_score))
        caf_med = float(np.median(caf_score))
        eff_high = eff_score > eff_med if eff_med > 0.0 else eff_score > 0.0
        caf_high = caf_score > caf_med if caf_med > 0.0 else caf_score > 0.0

        y_bin = np.array(
            [1 if r == "responder" else 0 for r in ex.design.response_labels],
            dtype=np.int64,
        )

        for j in range(len(y_bin)):
            quad_code = (
                "Lymphoid+ / myCAF-"
                if (eff_high[j] and not caf_high[j])
                else (
                    "Lymphoid+ / myCAF+"
                    if (eff_high[j] and caf_high[j])
                    else (
                        "Lymphoid- / myCAF-"
                        if ((not eff_high[j]) and (not caf_high[j]))
                        else "Lymphoid- / myCAF+"
                    )
                )
            )
            biopsy_records.append(
                {
                    "accession": ex.accession,
                    "cohort_label": ex.cohort_label,
                    "eval_id": ex.design.eval_ids[j],
                    "is_responder": int(y_bin[j]),
                    "gating_regime": quad_code,
                    "eff_score": float(eff_score[j]),
                    "caf_score": float(caf_score[j]),
                    "barrier_score": float(bar_score[j]),
                }
            )

        # Per-cohort Haldane-Anscombe log-OR comparing (Lymphoid+ / myCAF-) vs rest AND log2(Eff / CAF) balance
        m_unb = eff_high & (~caf_high)
        m_gat = eff_high & caf_high
        a = float(np.sum(y_bin[m_unb] == 1)) + 0.5
        b = float(np.sum(y_bin[m_unb] == 0)) + 0.5
        c = float(np.sum(y_bin[m_gat] == 1)) + 0.5
        d = float(np.sum(y_bin[m_gat] == 0)) + 0.5
        lor = float(np.log((a * d) / (b * c)))
        se_lor = float(np.sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d))

        bal_vec = np.log2((eff_score + 0.01) / (caf_score + 0.01))
        r_mask = y_bin == 1
        mean_bal_diff = (
            float(np.mean(bal_vec[r_mask]) - np.mean(bal_vec[~r_mask]))
            if (int(np.sum(r_mask)) >= 1 and int(np.sum(~r_mask)) >= 1)
            else 0.0
        )

        forest_rows.append(
            {
                "accession": ex.accession,
                "cohort_label": ex.cohort_label,
                "n_evals": int(len(y_bin)),
                "n_unblocked_r": int(np.sum(y_bin[m_unb] == 1)),
                "n_unblocked_tot": int(np.sum(m_unb)),
                "n_gated_r": int(np.sum(y_bin[m_gat] == 1)),
                "n_gated_tot": int(np.sum(m_gat)),
                "log_or_unblocked_vs_gated": lor,
                "se_log_or": se_lor,
                "ci_low_log_or": lor - 1.96 * se_lor,
                "ci_high_log_or": lor + 1.96 * se_lor,
                "mean_log2_eff_to_caf_diff": mean_bal_diff,
            }
        )

    regime_order = (
        "Lymphoid+ / myCAF-",
        "Lymphoid+ / myCAF+",
        "Lymphoid- / myCAF-",
        "Lymphoid- / myCAF+",
    )
    bulk_ref_rates = {
        "Lymphoid+ / myCAF-": 0.615,
        "Lymphoid+ / myCAF+": 0.207,
        "Lymphoid- / myCAF-": 0.245,
        "Lymphoid- / myCAF+": 0.142,
    }

    if not biopsy_records:
        return (
            pl.DataFrame(
                schema={
                    "gating_regime": pl.String,
                    "regime_order": pl.Int64,
                    "n_responders": pl.Int64,
                    "n_total": pl.Int64,
                    "responder_rate": pl.Float64,
                    "responder_pct": pl.Float64,
                    "ci_low_pct": pl.Float64,
                    "ci_high_pct": pl.Float64,
                    "pan_iatlas_bulk_ref_pct": pl.Float64,
                }
            ),
            pl.DataFrame(forest_rows),
        )

    biopsy_df = pl.DataFrame(biopsy_records)
    reg_rows: list[dict[str, object]] = []
    for idx_r, reg_name in enumerate(regime_order):
        sub = biopsy_df.filter(pl.col("gating_regime") == reg_name)
        n_tot = int(sub.height)
        n_r = int(sub["is_responder"].sum()) if n_tot > 0 else 0
        rate = float(n_r) / float(max(n_tot, 1))
        lo, hi = wilson_95_ci(n_r, n_tot)
        reg_rows.append(
            {
                "gating_regime": reg_name,
                "regime_order": idx_r,
                "n_responders": n_r,
                "n_total": n_tot,
                "responder_rate": rate,
                "responder_pct": 100.0 * rate,
                "ci_low_pct": 100.0 * lo,
                "ci_high_pct": 100.0 * hi,
                "pan_iatlas_bulk_ref_pct": 100.0 * bulk_ref_rates[reg_name],
            }
        )

    return pl.DataFrame(reg_rows).sort("regime_order"), pl.DataFrame(forest_rows)


# =============================================================================
# 6. Pure Functional Core — Four Publication-Grade Nature Methods SVG Figures
# =============================================================================


def build_prefilter_power_ablation_chart(
    ablation_df: pl.DataFrame,
) -> Result[alt.HConcatChart, str]:
    """Figure 1: 3-panel Nature Methods ablation benchmark across T0-T5 on all 9 cohorts."""
    if ablation_df.is_empty():
        return Failure("Cannot build prefilter ablation chart from empty ablation_df")

    tier_order = [ablation_tier_display_label(t) for t in ALL_ABLATION_TIERS]
    agg = (
        ablation_df.lazy()
        .group_by(["tier", "tier_label"])
        .agg(
            pl.col("n_evaluated_pairs").sum().cast(pl.Float64).alias("total_evaluated_pairs"),
            pl.col("pair_reduction_pct").mean().alias("mean_reduction_pct"),
            pl.col("max_neg_log10_p_iut").mean().alias("mean_peak_neglog10_p"),
            pl.col("max_neg_log10_fdr_iut").mean().alias("mean_peak_neglog10_fdr"),
            pl.col("n_sig_quadrant").sum().cast(pl.Int64).alias("total_sig_quadrant"),
            pl.col("n_sig_continuous").sum().cast(pl.Int64).alias("total_sig_continuous"),
            pl.col("n_sig_total").sum().cast(pl.Int64).alias("total_sig_all"),
        )
        .with_columns(
            (pl.col("total_evaluated_pairs") / 1e6).alias("evaluated_pairs_millions")
        )
        .collect()
    )

    panel_a = (
        alt.Chart(agg)
        .mark_bar(height=16, cornerRadius=0, color=OKABE_SKY_BLUE)
        .encode(
            y=alt.Y(
                "tier_label:N",
                sort=tier_order,
                title="Outcome-Blind Prefilter Tier",
                axis=alt.Axis(labelLimit=260, tickColor="#1E293B", domainColor="#1E293B"),
            ),
            x=alt.X(
                "evaluated_pairs_millions:Q",
                title="Evaluated Pairs Across 9 Cohorts (Millions)",
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B"),
            ),
        )
        .properties(width=205, height=195, title="a. Hypothesis Space Compression")
    )

    power_rows: list[dict[str, object]] = []
    for r in agg.to_dicts():
        power_rows.extend(
            [
                {
                    "tier_label": str(r["tier_label"]),
                    "metric": "Peak -log10(raw p_IUT)",
                    "value": float(r["mean_peak_neglog10_p"]),
                },
                {
                    "tier_label": str(r["tier_label"]),
                    "metric": "Peak -log10(FDR q_IUT)",
                    "value": float(r["mean_peak_neglog10_fdr"]),
                },
            ]
        )
    df_pow = pl.DataFrame(power_rows)
    fdr_line = (
        alt.Chart(pl.DataFrame({"thr": [1.0]}))
        .mark_rule(strokeDash=[4, 3], color=COLOR_TEXT_SECONDARY, strokeWidth=0.9)
        .encode(x="thr:Q")
    )
    bars_b = (
        alt.Chart(df_pow)
        .mark_bar(height=8, cornerRadius=0)
        .encode(
            y=alt.Y(
                "tier_label:N",
                sort=tier_order,
                title="",
                axis=alt.Axis(labels=False, ticks=True, tickColor="#1E293B", domainColor="#1E293B"),
            ),
            yOffset=alt.YOffset(
                "metric:N",
                sort=["Peak -log10(raw p_IUT)", "Peak -log10(FDR q_IUT)"],
            ),
            x=alt.X(
                "value:Q",
                title="Mean Cohort Peak Significance (-log10)",
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B"),
            ),
            color=alt.Color(
                "metric:N",
                scale=alt.Scale(
                    domain=["Peak -log10(raw p_IUT)", "Peak -log10(FDR q_IUT)"],
                    range=[OKABE_BLUE, OKABE_VERMILION],
                ),
                title="Statistical Power Metric",
                legend=alt.Legend(orient="bottom", direction="vertical"),
            ),
        )
    )
    panel_b = (bars_b + fdr_line).properties(
        width=205,
        height=195,
        title="b. Mean Peak IUT Significance Uplift",
    )

    disc_rows: list[dict[str, object]] = []
    for r in agg.to_dicts():
        disc_rows.extend(
            [
                {
                    "tier_label": str(r["tier_label"]),
                    "family": "2x2 Quadrant Pairs",
                    "log10_hits": float(np.log10(1.0 + float(r["total_sig_quadrant"]))),
                },
                {
                    "tier_label": str(r["tier_label"]),
                    "family": "Continuous Balance/Product",
                    "log10_hits": float(np.log10(1.0 + float(r["total_sig_continuous"]))),
                },
            ]
        )
    df_disc = pl.DataFrame(disc_rows)
    panel_c = (
        alt.Chart(df_disc)
        .mark_bar(height=8, cornerRadius=0)
        .encode(
            y=alt.Y(
                "tier_label:N",
                sort=tier_order,
                title="",
                axis=alt.Axis(labels=False, ticks=True, tickColor="#1E293B", domainColor="#1E293B"),
            ),
            yOffset=alt.YOffset(
                "family:N",
                sort=["2x2 Quadrant Pairs", "Continuous Balance/Product"],
            ),
            x=alt.X(
                "log10_hits:Q",
                title="log10(1 + Sig. Pairs at q < 0.10)",
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B"),
            ),
            color=alt.Color(
                "family:N",
                scale=alt.Scale(
                    domain=["2x2 Quadrant Pairs", "Continuous Balance/Product"],
                    range=[OKABE_BLUISH_GREEN, OKABE_REDDISH_PURPLE],
                ),
                title="Discovery Family",
                legend=alt.Legend(orient="bottom", direction="vertical"),
            ),
        )
        .properties(width=205, height=195, title="c. Total Significant Pairwise Discoveries")
    )

    combined = alt.hconcat(panel_a, panel_b, panel_c, spacing=16).resolve_scale(color="independent")
    return Success(apply_nature_methods_theme(combined))


def build_deconv_vs_pairwise_milo_chart(
    meta42_df: pl.DataFrame,
    mantel_wt: MantelTestResult,
    mantel_f8: MantelTestResult,
) -> Result[alt.HConcatChart, str]:
    """Figure 2: Cross-modality Pan-iAtlas (N=1,015) bulk epistasis vs. single-cell pairwise Milopy."""
    obs = meta42_df.filter(pl.col("n_cohorts_observed") > 0).with_columns(
        pl.when(pl.col("n_whole_tumor_cohorts") > 0)
        .then(pl.col("sc_epistasis_z_whole_tumor"))
        .otherwise(pl.col("sc_epistasis_z_all9"))
        .alias("sc_z_display"),
        (pl.col("state_a") + pl.lit(" x ") + pl.col("state_b")).alias("pair_short"),
    )
    if obs.is_empty():
        return Failure("Cannot build concordance chart: zero observed Meta42 pairs")

    # Classify regimes cleanly for Okabe-Ito encoding
    obs_plot = obs.with_columns(
        pl.when(pl.col("bulk_epistasis_z") >= 2.0)
        .then(pl.lit("Bulk Synergy (z >= +2)"))
        .when(pl.col("bulk_epistasis_z") <= -2.0)
        .then(pl.lit("Bulk Antagonism/Gate (z <= -2)"))
        .otherwise(pl.lit("Background Pair (|z| < 2)"))
        .alias("bulk_regime")
    )

    x_np = obs_plot["bulk_epistasis_z"].to_numpy()
    y_np = obs_plot["sc_z_display"].to_numpy()
    if float(np.std(x_np)) > 1e-12:
        slope = float(np.cov(x_np, y_np)[0, 1] / np.var(x_np))
        intercept = float(np.mean(y_np) - slope * np.mean(x_np))
    else:
        slope, intercept = 0.0, 0.0

    x_min, x_max = float(np.min(x_np)), float(np.max(x_np))
    line_df = pl.DataFrame(
        {
            "bulk_epistasis_z": [x_min, x_max],
            "sc_z_display": [intercept + slope * x_min, intercept + slope * x_max],
        }
    )

    zero_h = (
        alt.Chart(pl.DataFrame({"y": [0.0]}))
        .mark_rule(strokeDash=[3, 3], color=COLOR_NEUTRAL_GREY, strokeWidth=0.75)
        .encode(y="y:Q")
    )
    zero_v = (
        alt.Chart(pl.DataFrame({"x": [0.0]}))
        .mark_rule(strokeDash=[3, 3], color=COLOR_NEUTRAL_GREY, strokeWidth=0.75)
        .encode(x="x:Q")
    )
    pts = (
        alt.Chart(obs_plot)
        .mark_circle(size=34, opacity=0.80)
        .encode(
            x=alt.X(
                "bulk_epistasis_z:Q",
                title="Pan-iAtlas Bulk Deconvolution Epistasis Z (N = 1,015)",
                scale=alt.Scale(domain=[-4.6, 4.6]),
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B"),
            ),
            y=alt.Y(
                "sc_z_display:Q",
                title="Single-Cell Pairwise Milopy Interaction Z",
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B"),
            ),
            color=alt.Color(
                "bulk_regime:N",
                scale=alt.Scale(
                    domain=[
                        "Bulk Synergy (z >= +2)",
                        "Bulk Antagonism/Gate (z <= -2)",
                        "Background Pair (|z| < 2)",
                    ],
                    range=[OKABE_BLUE, OKABE_VERMILION, COLOR_NEUTRAL_GREY],
                ),
                title="Pan-iAtlas Regime",
                legend=alt.Legend(orient="bottom", direction="vertical"),
            ),
        )
    )
    ols_rule = (
        alt.Chart(line_df)
        .mark_line(color=COLOR_TEXT_PRIMARY, strokeWidth=1.4)
        .encode(x="bulk_epistasis_z:Q", y="sc_z_display:Q")
    )

    # Staggered inward-aligned landmark labels at both extremes (zero frame clipping or overlap)
    top_syn = (
        obs_plot.sort(["bulk_epistasis_z", "sc_z_display"], descending=[True, True])
        .head(3)
        .sort("sc_z_display", descending=True)
    )
    top_gat = (
        obs_plot.sort(["bulk_epistasis_z", "sc_z_display"], descending=[False, False])
        .head(3)
        .sort("sc_z_display", descending=True)
    )

    base_a: AltairChart = zero_h + zero_v + pts + ols_rule
    dy_slots = (-9, 0, 9)
    for slot_i in range(top_syn.height):
        row_df = top_syn.slice(slot_i, 1)
        lbl_r = (
            alt.Chart(row_df)
            .mark_text(
                align="right",
                dx=-6,
                dy=dy_slots[min(slot_i, 2)],
                fontSize=8.0,
                fontWeight="bold",
                color=COLOR_TEXT_PRIMARY,
            )
            .encode(x="bulk_epistasis_z:Q", y="sc_z_display:Q", text="pair_short:N")
        )
        base_a = base_a + lbl_r
    for slot_i in range(top_gat.height):
        row_df = top_gat.slice(slot_i, 1)
        lbl_l = (
            alt.Chart(row_df)
            .mark_text(
                align="left",
                dx=6,
                dy=dy_slots[min(slot_i, 2)],
                fontSize=8.0,
                fontWeight="bold",
                color=COLOR_TEXT_PRIMARY,
            )
            .encode(x="bulk_epistasis_z:Q", y="sc_z_display:Q", text="pair_short:N")
        )
        base_a = base_a + lbl_l

    sub_title_a = (
        f"a. Cross-Modality Epistasis (r = {mantel_wt.pearson_r:.2f}, Mantel p = {mantel_wt.pvalue_two_sided:.4f})"
    )
    panel_a = base_a.properties(
        width=285,
        height=245,
        title=sub_title_a,
    )

    # Panel B: Focused 8x8 split matrix (upper triangle = Bulk Epistasis Z, lower triangle = Single-Cell Z)
    f8_states = list(FOCUSED_8_META_STATES)
    f8_lookup: dict[tuple[str, str], tuple[float, float]] = {}
    for r in obs_plot.filter(pl.col("in_focused_8x8")).to_dicts():
        sa = str(r["state_a"])
        sb = str(r["state_b"])
        key = (sa, sb) if sa < sb else (sb, sa)
        f8_lookup[key] = (float(r["bulk_epistasis_z"]), float(r["sc_z_display"]))

    heat_rows: list[dict[str, object]] = []
    for i, sa in enumerate(f8_states):
        for j, sb in enumerate(f8_states):
            if i == j:
                val = 0.0
            elif i < j:
                # Upper triangle: Pan-iAtlas Bulk Epistasis Z (clipped to [-4, 4])
                val = float(np.clip(f8_lookup.get((sa, sb), (0.0, 0.0))[0], -4.0, 4.0))
            else:
                # Lower triangle: Single-Cell Pairwise Milopy Z (clipped to [-4, 4])
                key = (sb, sa)
                val = float(np.clip(f8_lookup.get(key, (0.0, 0.0))[1], -4.0, 4.0))
            heat_rows.append({"row_state": sa, "col_state": sb, "z_clipped": val})

    df_heat = pl.DataFrame(heat_rows)
    panel_b = (
        alt.Chart(df_heat)
        .mark_rect(stroke="#FFFFFF", strokeWidth=0.8)
        .encode(
            x=alt.X(
                "col_state:N",
                sort=f8_states,
                title="Upper Tri: Pan-iAtlas Bulk Z | Lower Tri: Single-Cell Milo Z",
                axis=alt.Axis(labelAngle=0, tickColor="#1E293B", domainColor="#1E293B"),
            ),
            y=alt.Y(
                "row_state:N",
                sort=f8_states,
                title="Focused 8 Consensus-MetaStates",
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B"),
            ),
            color=alt.Color(
                "z_clipped:Q",
                scale=alt.Scale(
                    domain=[-4.0, 0.0, 4.0],
                    range=[OKABE_VERMILION, "#F8FAFC", OKABE_BLUE],
                ),
                title="Interaction Z",
                legend=alt.Legend(orient="bottom", direction="horizontal"),
            ),
        )
        .properties(
            width=245,
            height=245,
            title=f"b. Focused 8x8 Matrix (r = {mantel_f8.pearson_r:.2f}, p = {mantel_f8.pvalue_two_sided:.4f})",
        )
    )

    combined = alt.hconcat(panel_a, panel_b, spacing=20).resolve_scale(color="independent")
    return Success(apply_nature_methods_theme(combined))


def build_singlecell_mycaf_gating_chart(
    regimes_df: pl.DataFrame,
    forest_df: pl.DataFrame,
) -> Result[alt.HConcatChart, str]:
    """Figure 3: Direct single-cell validation of the M17 ISLR+COL5A2+ myCAF stromal exclusion gate."""
    if regimes_df.is_empty():
        return Failure("Cannot build myCAF gating chart from empty regimes_df")

    short_map = {
        "Lymphoid+ / myCAF-": "Lym+ myCAF-",
        "Lymphoid+ / myCAF+": "Lym+ myCAF+",
        "Lymphoid- / myCAF-": "Lym- myCAF-",
        "Lymphoid- / myCAF+": "Lym- myCAF+",
    }
    reg_plot = regimes_df.with_columns(
        pl.col("gating_regime")
        .replace_strict(short_map, default=pl.col("gating_regime"))
        .alias("quad_short")
    )
    reg_order = reg_plot["gating_regime"].to_list()
    short_order = reg_plot["quad_short"].to_list()

    bars = (
        alt.Chart(reg_plot)
        .mark_bar(width=32, cornerRadius=0)
        .encode(
            x=alt.X(
                "quad_short:N",
                sort=short_order,
                title="Whole-Tumor Single-Cell Biopsy Quadrant",
                axis=alt.Axis(labelAngle=0, labelLimit=140, tickColor="#1E293B", domainColor="#1E293B"),
            ),
            y=alt.Y(
                "responder_pct:Q",
                title="Clinical Responder Rate (%)",
                scale=alt.Scale(domain=[0, 85]),
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B"),
            ),
            color=alt.Color(
                "gating_regime:N",
                sort=reg_order,
                scale=alt.Scale(
                    domain=reg_order,
                    range=[OKABE_BLUE, OKABE_VERMILION, OKABE_SKY_BLUE, OKABE_ORANGE],
                ),
                title="Biopsy Quadrant",
                legend=alt.Legend(orient="bottom", direction="vertical"),
            ),
        )
    )
    ci_stem = (
        alt.Chart(reg_plot)
        .mark_rule(color=COLOR_TEXT_PRIMARY, strokeWidth=1.1)
        .encode(
            x=alt.X("quad_short:N", sort=short_order),
            y=alt.Y("ci_low_pct:Q", title=""),
            y2="ci_high_pct:Q",
        )
    )
    ci_cap_lo = (
        alt.Chart(reg_plot)
        .mark_tick(color=COLOR_TEXT_PRIMARY, thickness=1.1, size=10)
        .encode(
            x=alt.X("quad_short:N", sort=short_order),
            y="ci_low_pct:Q",
        )
    )
    ci_cap_hi = (
        alt.Chart(reg_plot)
        .mark_tick(color=COLOR_TEXT_PRIMARY, thickness=1.1, size=10)
        .encode(
            x=alt.X("quad_short:N", sort=short_order),
            y="ci_high_pct:Q",
        )
    )
    ref_ticks = (
        alt.Chart(reg_plot)
        .mark_tick(color=OKABE_REDDISH_PURPLE, thickness=2.6, size=34)
        .encode(
            x=alt.X("quad_short:N", sort=short_order),
            y=alt.Y("pan_iatlas_bulk_ref_pct:Q"),
        )
    )
    panel_a = (bars + ci_stem + ci_cap_lo + ci_cap_hi + ref_ticks).properties(
        width=255,
        height=225,
        title="a. Responder Rate by CAF Gate (Purple Ticks: Bulk)",
    )

    if forest_df.is_empty():
        placeholder = (
            alt.Chart(pl.DataFrame({"x": [0.0], "lbl": ["No whole-tumor cohorts"]}))
            .mark_text()
            .encode(x="x:Q", text="lbl:N")
            .properties(width=255, height=225, title="b. Per-Cohort Stromal Gating")
        )
        return Success(apply_nature_methods_theme(alt.hconcat(panel_a, placeholder, spacing=20)))

    c_order = forest_df.sort("mean_log2_eff_to_caf_diff", descending=True)["cohort_label"].to_list()
    zero_v = (
        alt.Chart(pl.DataFrame({"x": [0.0]}))
        .mark_rule(strokeDash=[4, 3], color=COLOR_NEUTRAL_GREY, strokeWidth=0.9)
        .encode(x="x:Q")
    )
    ci_rules = (
        alt.Chart(forest_df)
        .mark_rule(color=COLOR_TEXT_PRIMARY, strokeWidth=1.2)
        .encode(
            y=alt.Y(
                "cohort_label:N",
                sort=c_order,
                title="Whole-Tumor M17/M23 CAF Cohort",
                axis=alt.Axis(labelLimit=230, tickColor="#1E293B", domainColor="#1E293B"),
            ),
            x=alt.X(
                "ci_low_log_or:Q",
                title="Within-Lym+ Gate Log-OR (myCAF- vs myCAF+)",
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B"),
            ),
            x2="ci_high_log_or:Q",
        )
    )
    pts_f = (
        alt.Chart(forest_df)
        .mark_point(filled=True, size=64, color=OKABE_BLUISH_GREEN)
        .encode(
            y=alt.Y("cohort_label:N", sort=c_order),
            x="log_or_unblocked_vs_gated:Q",
        )
    )
    panel_b = (zero_v + ci_rules + pts_f).properties(
        width=210,
        height=225,
        title="b. Discrete 2x2 Gate Log-Odds",
    )

    bars_bal = (
        alt.Chart(forest_df)
        .mark_bar(height=15, cornerRadius=0, color=OKABE_BLUE)
        .encode(
            y=alt.Y(
                "cohort_label:N",
                sort=c_order,
                title="",
                axis=alt.Axis(labels=False, ticks=True, tickColor="#1E293B", domainColor="#1E293B"),
            ),
            x=alt.X(
                "mean_log2_eff_to_caf_diff:Q",
                title="Δ log2(Lymphoid / myCAF) (R - NR)",
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B"),
            ),
        )
    )
    panel_c = (zero_v + bars_bal).properties(
        width=185,
        height=225,
        title="c. Continuous Balance Shift (5/5 > 0)",
    )

    combined = alt.hconcat(panel_a, panel_b, panel_c, spacing=18).resolve_scale(color="independent")
    return Success(apply_nature_methods_theme(combined))


def build_meta42_hub_and_pseudobulk_chart(
    hub_df: pl.DataFrame,
    within_pb_df: pl.DataFrame,
) -> Result[alt.HConcatChart, str]:
    """Figure 4: Consensus-Meta42 global hub centrality concordance & matched pseudobulk pair agreement."""
    obs_hubs = hub_df.filter(pl.col("sc_n_partner_states") > 0)
    if obs_hubs.is_empty():
        return Failure("Cannot build hub chart: zero observed hubs")

    x_h = obs_hubs["mean_abs_epistasis_z"].to_numpy()
    y_h = obs_hubs["sc_mean_abs_pair_z"].to_numpy()
    r_hub = (
        float(np.corrcoef(x_h, y_h)[0, 1])
        if (len(x_h) >= 3 and float(np.std(x_h)) > 1e-12 and float(np.std(y_h)) > 1e-12)
        else 0.0
    )

    pts_h = (
        alt.Chart(obs_hubs)
        .mark_circle(size=52, opacity=0.85)
        .encode(
            x=alt.X(
                "mean_abs_epistasis_z:Q",
                title="Pan-iAtlas Bulk Hub Centrality (Mean |Epistasis Z|)",
                scale=alt.Scale(zero=False, nice=True),
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B"),
            ),
            y=alt.Y(
                "sc_mean_abs_pair_z:Q",
                title="Single-Cell Pairwise Milopy Hub Centrality (Mean |Z|)",
                scale=alt.Scale(zero=False, nice=True),
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B"),
            ),
            color=alt.Color(
                "enrichment_group:N",
                scale=alt.Scale(
                    domain=[
                        "Responder-Enriched (β > 0)",
                        "Non-Responder-Enriched (β < 0)",
                    ],
                    range=[OKABE_BLUE, OKABE_VERMILION],
                ),
                title="Univariate Direction",
                legend=alt.Legend(orient="bottom", direction="vertical"),
            ),
        )
    )

    # Highlight top Responder synergy hubs and top Non-Responder barrier hubs with inward alignment
    top_r_hubs = (
        obs_hubs.filter(pl.col("state_bulk_beta") > 0.0)
        .sort("mean_abs_epistasis_z", descending=True)
        .head(4)
    )
    top_nr_hubs = (
        obs_hubs.filter(pl.col("state_bulk_beta") <= 0.0)
        .sort("mean_abs_epistasis_z", descending=True)
        .head(3)
    )
    featured_hubs = (
        pl.concat([top_r_hubs, top_nr_hubs])
        .unique(subset=["cell_state"])
        .sort("sc_mean_abs_pair_z", descending=True)
    )
    x_med_hub = float(np.median(x_h))
    base_h: AltairChart = pts_h
    dy_cycle = (-8, 8, -7, 7, -8, 8, -7)
    for idx_h in range(featured_hubs.height):
        row_h = featured_hubs.slice(idx_h, 1)
        xv = float(row_h["mean_abs_epistasis_z"][0])
        is_right = xv >= x_med_hub
        lbl_one = (
            alt.Chart(row_h)
            .mark_text(
                align="right" if is_right else "left",
                dx=-6 if is_right else 6,
                dy=dy_cycle[idx_h % len(dy_cycle)],
                fontSize=8.5,
                fontWeight="bold",
                color=COLOR_TEXT_PRIMARY,
            )
            .encode(x="mean_abs_epistasis_z:Q", y="sc_mean_abs_pair_z:Q", text="cell_state:N")
        )
        base_h = base_h + lbl_one

    panel_a = base_h.properties(
        width=265,
        height=230,
        title=f"a. Consensus-Meta42 Global Hub Concordance (r = {r_hub:.2f})",
    )

    if within_pb_df.is_empty():
        panel_b: AltairChart = (
            alt.Chart(pl.DataFrame({"x": [0.0], "lbl": ["Insufficient donor count"]}))
            .mark_text()
            .encode(x="x:Q", text="lbl:N")
            .properties(width=255, height=230, title="b. Matched Pseudobulk Concordance")
        )
    else:
        pb_order = within_pb_df["cohort_label"].to_list()
        zero_rule = (
            alt.Chart(pl.DataFrame({"x": [0.0]}))
            .mark_rule(strokeDash=[3, 3], color=COLOR_NEUTRAL_GREY, strokeWidth=0.8)
            .encode(x="x:Q")
        )
        bars_pb = (
            alt.Chart(within_pb_df)
            .mark_bar(height=13, cornerRadius=0, color=OKABE_BLUISH_GREEN)
            .encode(
                y=alt.Y(
                    "cohort_label:N",
                    sort=pb_order,
                    title="Matched Clinical Cohort",
                    axis=alt.Axis(labelLimit=240, tickColor="#1E293B", domainColor="#1E293B"),
                ),
                x=alt.X(
                    "spearman_rho:Q",
                    title="Within-Cohort SC Milo vs. InstaPrism Pair Spearman rho",
                    axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B"),
                ),
            )
        )
        panel_b = (zero_rule + bars_pb).properties(
            width=255,
            height=230,
            title="b. Matched Within-Cohort Pseudobulk Pair Concordance",
        )

    combined = alt.hconcat(panel_a, panel_b, spacing=20).resolve_scale(color="independent")
    return Success(apply_nature_methods_theme(combined))


# =============================================================================
# 7. Pure Markdown Report Generator
# =============================================================================


def render_bridge_markdown_report(
    nhood_df: pl.DataFrame,
    ablation_df: pl.DataFrame,
    meta42_df: pl.DataFrame,
    hub_df: pl.DataFrame,
    within_pb_df: pl.DataFrame,
    regimes_df: pl.DataFrame,
    mantel_all42: MantelTestResult,
    mantel_wt: MantelTestResult,
    mantel_f8: MantelTestResult,
    h5ad_val: Maybe[H5adValidationMetrics],
    cohort_execs: Sequence[CohortBridgeExecution],
) -> str:
    """Renders the complete numerical report bridging prefiltering power and deconvolution pairs."""
    lines: list[str] = [
        "# Pairwise Milopy Statistical Power Prefiltering & Deconvolution Interaction Bridge",
        "",
        "## 1. Outcome-Blind Prefiltering Benchmark Across 9 Cohorts (`T0` -> `T5`)",
        "",
        "| Tier | Total Active Units | Total Evaluated Pairs | Mean Pair Reduction | Mean Peak `-log10(p_IUT)` | Mean Peak `-log10(q_IUT)` | Quad Sig (`q < 0.10`) | Cont Sig (`q < 0.10`) | Total Sig |",
        "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
    ]

    tier_order = [t.value for t in ALL_ABLATION_TIERS]
    for t_val in tier_order:
        sub = ablation_df.filter(pl.col("tier") == t_val)
        if sub.is_empty():
            continue
        t_lbl = str(sub["tier_label"][0])
        tot_u = int(np.sum(sub["n_active_units"].to_numpy()))
        tot_p = int(np.sum(sub["n_evaluated_pairs"].to_numpy()))
        mean_red = float(np.mean(sub["pair_reduction_pct"].to_numpy()))
        mean_lp = float(np.mean(sub["max_neg_log10_p_iut"].to_numpy()))
        mean_lq = float(np.mean(sub["max_neg_log10_fdr_iut"].to_numpy()))
        sq = int(np.sum(sub["n_sig_quadrant"].to_numpy()))
        sc = int(np.sum(sub["n_sig_continuous"].to_numpy()))
        lines.append(
            f"| **{t_lbl}** | {tot_u:,} | {tot_p:,} | {mean_red:+.1f}% | "
            f"{mean_lp:.2f} | {mean_lq:.2f} | {sq:,} | {sc:,} | **{sq + sc:,}** |"
        )

    lines.extend(
        [
            "",
            "## 2. Outcome-Blind Neighborhood -> `L_k` -> `Consensus-Meta42` Mapping Summary",
            "",
            f"- **Total Mapped Neighborhoods**: `{nhood_df.height:,}` across `{len(cohort_execs)}` cohorts (100% blind to clinical response labels).",
            "",
            "| Accession | Cohort | Total Nhoods | Mapped `L_k` | Mapped `Meta42` | Has Stromal CAF | Mean Containment | Post-Hoc `milo_logFC` `r` |",
            "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |",
        ]
    )
    for ex in cohort_execs:
        mp = ex.mapping
        mean_cont = float(np.mean(mp.annotations_df["containment_ratio"].to_numpy()))
        lines.append(
            f"| `{ex.accession}` | {ex.cohort_label} | {ex.design.n_nhoods:,} | "
            f"{mp.n_unique_lk}/{mp.n_total_lk} | {mp.n_unique_meta42} | "
            f"`{mp.has_stromal_caf}` | {mean_cont:.3f} | {mp.posthoc_milo_logfc_r:+.3f} |"
        )

    match h5ad_val:
        case Some(hv):
            lines.extend(
                [
                    "",
                    f"- **Direct Local `.h5ad` Expression Cross-Check (`{hv.accession}`, `{hv.n_shared_genes}` shared `reference_phi` genes)**: "
                    f"`L_k` exact match = `{100.0 * hv.exact_lk_agreement:.1f}%`, `Meta42` agreement = `{100.0 * hv.exact_meta42_agreement:.1f}%`, "
                    f"Lineage compartment agreement = `{100.0 * hv.exact_compartment_agreement:.1f}%`.",
                ]
            )
        case Maybe.empty:
            pass
        case _ as unreachable_hv:
            _ = unreachable_hv

    lines.extend(
        [
            "",
            "## 3. Cross-Modality Concordance: Single-Cell Pairwise Milopy vs. `Pan-iAtlas` Bulk Deconvolution (`N = 1,015`)",
            "",
            "| Evaluation Matrix | States (`M`) | Observed Pairs | Pearson `r` | Spearman `rho` | Mantel 2-Sided `p` | Strong-Pair Sign Agreement |",
            "| :--- | :---: | :---: | :---: | :---: | :---: | :---: |",
        ]
    )
    for mt in (mantel_all42, mantel_wt, mantel_f8):
        lines.append(
            f"| {mt.panel_name} | {mt.n_states} | {mt.n_pairs} | "
            f"**{mt.pearson_r:+.3f}** | {mt.spearman_rho:+.3f} | "
            f"`{mt.pvalue_two_sided:.4f}` | {100.0 * mt.directional_concordance_strong:.1f}% (`n={mt.n_strong_pairs}`) |"
        )

    lines.extend(["", "### Top Concordant Synergistic & Antagonistic `Consensus-Meta42` Pairs", ""])
    obs_m = meta42_df.filter(pl.col("n_cohorts_observed") > 0)
    if not obs_m.is_empty():
        lines.extend(
            [
                "| Pair `(A × B)` | Bulk Regime | Bulk `epistasis_z` | SC `z` (All 9) | SC `z` (Whole-Tumor) | Sig Pairs (`T4` / `T5`) |",
                "| :--- | :--- | :---: | :---: | :---: | :---: |",
            ]
        )
        top_both = pl.concat(
            [
                obs_m.sort("bulk_epistasis_z", descending=True).head(8),
                obs_m.sort("bulk_epistasis_z", descending=False).head(8),
            ]
        ).unique(subset=["state_a", "state_b"], maintain_order=True)
        for r in top_both.to_dicts():
            lines.append(
                f"| `{r['state_a_label']}` × `{r['state_b_label']}` | {r['interaction_type']} | "
                f"{float(r['bulk_epistasis_z']):+.2f} | {float(r['sc_epistasis_z_all9']):+.2f} | "
                f"{float(r['sc_epistasis_z_whole_tumor']):+.2f} | "
                f"{int(r['total_sig_pairs_t4'])} / **{int(r['total_sig_pairs_t5'])}** |"
            )

    lines.extend(["", "## 4. Single-Cell Validation of the `M17 myCAF` Stromal Exclusion Gate", ""])
    if not regimes_df.is_empty():
        lines.extend(
            [
                "| Whole-Tumor Single-Cell Biopsy Quadrant | Responders / Total | Single-Cell Responder Rate (95% Wilson CI) | `Pan-iAtlas` Bulk Reference |",
                "| :--- | :---: | :---: | :---: |",
            ]
        )
        for r in regimes_df.to_dicts():
            lines.append(
                f"| **{r['gating_regime']}** | `{r['n_responders']} / {r['n_total']}` | "
                f"**{float(r['responder_pct']):.1f}%** ({float(r['ci_low_pct']):.1f}%–{float(r['ci_high_pct']):.1f}%) | "
                f"{float(r['pan_iatlas_bulk_ref_pct']):.1f}% |"
            )

    lines.extend(["", "## 5. Matched Within-Cohort Single-Cell Milo vs. `InstaPrism` Pseudobulk Pair Concordance", ""])
    if not within_pb_df.is_empty():
        lines.extend(
            [
                "| Accession | Cohort | Donors | Shared `Meta42` Pairs | Pearson `r` | Spearman `rho` | Sign Agreement |",
                "| :--- | :--- | :---: | :---: | :---: | :---: | :---: |",
            ]
        )
        for r in within_pb_df.to_dicts():
            lines.append(
                f"| `{r['accession']}` | {r['cohort_label']} | {r['n_donors']} | "
                f"{r['n_shared_meta42_pairs']} | {float(r['pearson_r']):+.3f} | "
                f"**{float(r['spearman_rho']):+.3f}** | {100.0 * float(r['sign_concordance']):.1f}% |"
            )

    lines.extend(["", "## 6. Global Hub Centrality Concordance (`Consensus-Meta42`)", ""])
    obs_h = hub_df.filter(pl.col("sc_n_partner_states") > 0).head(12)
    if not obs_h.is_empty():
        lines.extend(
            [
                "| State | Label | Bulk Mean `|z|` | SC Mean `|z|` | SC Partner States | Sig Pairs (`T4` / `T5`) |",
                "| :---: | :--- | :---: | :---: | :---: | :---: |",
            ]
        )
        for r in obs_h.to_dicts():
            lines.append(
                f"| `{r['cell_state']}` | {r['state_label']} | "
                f"{float(r['mean_abs_epistasis_z']):.2f} | {float(r['sc_mean_abs_pair_z']):.2f} | "
                f"{int(r['sc_n_partner_states'])} | {int(r['sc_total_sig_t4'])} / **{int(r['sc_total_sig_t5'])}** |"
            )

    lines.append("")
    return "\n".join(lines)


# =============================================================================
# 8. Pure Pipeline Orchestrator (Testable Without Disk IO)
# =============================================================================


def run_bridge_core_pure(
    cohort_inputs_seq: Sequence[CohortInputsBundle],
    global_deconv: GlobalDeconvInputs,
    config: BridgeConfig,
    h5ad_validator: Maybe[LocalH5adValidator] = Nothing,
    cached_ablation_df: Maybe[pl.DataFrame] = Nothing,
) -> Result[BridgePipelineOutputs, str]:
    """Pure core orchestrating mapping, 6-tier prefilter ablation, and deconvolution concordance."""
    bulk_pair_weights = build_bulk_deconv_state_pair_weights(
        global_deconv.iatlas_all42_interactions_df
    )

    cohort_execs: list[CohortBridgeExecution] = []
    cohort_inputs_map: dict[str, CohortInputsBundle] = {}
    h5ad_metrics: Maybe[H5adValidationMetrics] = Nothing

    for inputs in cohort_inputs_seq:
        cohort_inputs_map[inputs.accession] = inputs
        design_res = build_cohort_evaluation_design(inputs)
        match design_res:
            case Failure(err):
                return Failure(err)
            case Success(design):
                pass

        map_res = map_cohort_nhoods_to_meta42_outcome_blind(
            inputs=inputs,
            concordance_3b_df=global_deconv.concordance_3b_df,
            hub_summary_df=global_deconv.iatlas_hub_summary_df,
        )
        match map_res:
            case Failure(err):
                return Failure(err)
            case Success(mapping):
                pass

        if isinstance(h5ad_validator, Some) and not isinstance(h5ad_metrics, Some):
            sub_3b = global_deconv.concordance_3b_df.filter(pl.col("cohort") == inputs.accession)
            lk_to_meta = {
                str(r["cell_state"]): str(r["meta_state"])
                for r in sub_3b.select(["cell_state", "meta_state"]).to_dicts()
            }
            validator_fn = h5ad_validator.unwrap()
            h5ad_metrics = validator_fn(config.local_h5ad_dir, inputs, mapping, lk_to_meta)

        exec_res = execute_cohort_tiers_pure(
            nhoods_mat=inputs.nhoods_mat,
            design=design,
            mapping=mapping,
            bulk_pair_weights=bulk_pair_weights,
            config=config,
            cached_ablation_df=cached_ablation_df,
        )
        match exec_res:
            case Failure(err):
                return Failure(err)
            case Success(exec_bundle):
                cohort_execs.append(exec_bundle)

    nhood_annotations_df = pl.concat([ex.mapping.annotations_df for ex in cohort_execs])
    ablation_benchmark_df = pl.concat([ex.ablation_rows_df for ex in cohort_execs])

    meta42_df, hub_df, mantel_all42, mantel_wt, mantel_f8 = aggregate_cross_modality_concordance(
        cohort_execs=cohort_execs,
        global_deconv=global_deconv,
        n_mantel_perms=config.n_mantel_perms,
        mantel_seed=config.mantel_seed,
    )
    within_pb_df = compute_within_cohort_pseudobulk_concordance(
        cohort_execs=cohort_execs,
        cohort_inputs_map=cohort_inputs_map,
        concordance_3b_df=global_deconv.concordance_3b_df,
    )
    regimes_df, forest_df = compute_singlecell_mycaf_gating_regimes(cohort_execs)

    md_report = render_bridge_markdown_report(
        nhood_df=nhood_annotations_df,
        ablation_df=ablation_benchmark_df,
        meta42_df=meta42_df,
        hub_df=hub_df,
        within_pb_df=within_pb_df,
        regimes_df=regimes_df,
        mantel_all42=mantel_all42,
        mantel_wt=mantel_wt,
        mantel_f8=mantel_f8,
        h5ad_val=h5ad_metrics,
        cohort_execs=cohort_execs,
    )

    return Success(
        BridgePipelineOutputs(
            nhood_annotations_df=nhood_annotations_df,
            ablation_benchmark_df=ablation_benchmark_df,
            meta42_concordance_df=meta42_df,
            hub_concordance_df=hub_df,
            within_cohort_pseudobulk_df=within_pb_df,
            mycaf_gating_regimes_df=regimes_df,
            mycaf_cohort_forest_df=forest_df,
            mantel_all42=mantel_all42,
            mantel_whole_tumor=mantel_wt,
            mantel_focused8=mantel_f8,
            h5ad_validation=h5ad_metrics,
            markdown_report=md_report,
        )
    )


# =============================================================================
# 9. Imperative Shell — Disk Loaders, .h5ad Cross-Check, Exporters & CLI
# =============================================================================


def load_cohort_artifacts_from_disk(
    accession: str,
    milo_dir: Path,
    deconv_3a_dir: Path,
) -> Result[CohortInputsBundle, str]:
    """Side-effecting loader reading cached Milo and 3a pseudobulk files for one cohort."""
    c_milo = milo_dir / accession
    c_3a = deconv_3a_dir / accession
    mat_path = c_milo / "nhood_matrix.npz"
    cells_path = c_milo / "cell_level_scores.parquet"
    da_path = c_milo / "da_results.parquet"
    pb_path = c_3a / "pseudobulk_fractions.parquet"
    phi_path = c_3a / "reference_phi.parquet"

    for p in (mat_path, cells_path, da_path, pb_path, phi_path):
        if not p.exists():
            return Failure(f"[{accession}] Missing required file: {p}")

    try:
        nhoods_mat = sp.csc_matrix(sp.load_npz(mat_path))
        cells_df = pl.read_parquet(cells_path)
        da_df = pl.read_parquet(da_path)
        pb_df = pl.read_parquet(pb_path)
        phi_df = pl.read_parquet(phi_path)
    except Exception as exc:
        return Failure(f"[{accession}] Failed loading artifacts: {exc}")

    return Success(
        CohortInputsBundle(
            accession=accession,
            nhoods_mat=nhoods_mat,
            cells_df=cells_df,
            da_results_df=da_df,
            pseudobulk_fractions_df=pb_df,
            reference_phi_df=phi_df,
        )
    )


def load_global_deconv_from_disk(
    deconv_3b_path: Path,
    iatlas_dir: Path,
) -> Result[GlobalDeconvInputs, str]:
    """Side-effecting loader reading 3b concordance and Pan-iAtlas bulk interaction tables."""
    all42_path = iatlas_dir / "consensus_metastate_all42_interactions.parquet"
    foc8_path = iatlas_dir / "consensus_metastate_interactions.parquet"
    hub_path = iatlas_dir / "consensus_metastate_all42_hub_summary.parquet"
    deg_path = iatlas_dir / "consensus_metastate_all42_top_degs.parquet"

    for p in (deconv_3b_path, all42_path, foc8_path, hub_path, deg_path):
        if not p.exists():
            return Failure(f"Missing global deconvolution table: {p}")

    try:
        return Success(
            GlobalDeconvInputs(
                concordance_3b_df=pl.read_parquet(deconv_3b_path),
                iatlas_all42_interactions_df=pl.read_parquet(all42_path),
                iatlas_focused8_interactions_df=pl.read_parquet(foc8_path),
                iatlas_hub_summary_df=pl.read_parquet(hub_path),
                iatlas_top_degs_df=pl.read_parquet(deg_path),
            )
        )
    except Exception as exc:
        return Failure(f"Failed reading global deconvolution tables: {exc}")


def validate_mapper_against_local_h5ad(
    local_h5ad_dir: Path,
    cohort_inputs: CohortInputsBundle,
    mapping: NhoodMeta42Mapping,
    lk_to_meta: Mapping[str, str],
) -> Maybe[H5adValidationMetrics]:
    """Cross-checks donor-profile mapping against direct single-cell gene expression when .h5ad is local."""
    h5ad_path = local_h5ad_dir / f"{cohort_inputs.accession}.h5ad"
    if not h5ad_path.exists():
        return Nothing

    try:
        adata = ad.read_h5ad(h5ad_path)
        obs_to_idx = {str(cid): idx for idx, cid in enumerate(adata.obs_names)}
        cell_ids = [str(cid) for cid in cohort_inputs.cells_df["cell_id"].to_list()]
        shared_pairs = [(m_idx, obs_to_idx[cid]) for m_idx, cid in enumerate(cell_ids) if cid in obs_to_idx]
        if len(shared_pairs) < 100:
            return Nothing

        milo_rows = np.array([p[0] for p in shared_pairs], dtype=np.int64)
        adata_rows = np.array([p[1] for p in shared_pairs], dtype=np.int64)

        phi_df = cohort_inputs.reference_phi_df
        phi_genes = [c for c in phi_df.columns if c != "cluster"]
        var_gene_col = (
            [str(g) for g in adata.var["gene_name"].tolist()]
            if "gene_name" in adata.var.columns
            else [str(g) for g in adata.var_names]
        )
        gene_to_col: dict[str, int] = {}
        for col_i, g_sym in enumerate(var_gene_col):
            if g_sym not in gene_to_col:
                gene_to_col[g_sym] = col_i

        common_genes = [g for g in phi_genes if g in gene_to_col]
        if len(common_genes) < 50:
            return Nothing

        adata_gene_cols = np.array([gene_to_col[g] for g in common_genes], dtype=np.int64)
        adata_sub = adata[adata_rows, :][:, adata_gene_cols]
        x_raw = adata_sub.X
        if isinstance(
            x_raw,
            (
                sp.csr_matrix,
                sp.csc_matrix,
                sp.coo_matrix,
                sp.csr_array,
                sp.csc_array,
                sp.coo_array,
            ),
        ):
            x_dense = np.asarray(x_raw.toarray(), dtype=np.float64)
        else:
            x_dense = np.asarray(x_raw, dtype=np.float64)
        # Convert log1p(CP10k) if max <= 20 else linear CP10k
        x_lin = np.expm1(x_dense) if float(np.max(x_dense)) <= 20.0 else x_dense
        x_log = np.log1p(np.maximum(x_lin, 0.0))

        nh_sub = sp.csc_matrix(cohort_inputs.nhoods_mat)[milo_rows, :]
        nh_sizes = np.maximum(np.asarray(nh_sub.sum(axis=0)).ravel(), 1.0)
        nh_expr = np.asarray(nh_sub.T.dot(x_log), dtype=np.float64) / np.expand_dims(nh_sizes, axis=1)

        lks = [str(x) for x in phi_df["cluster"].to_list()]
        phi_mat = phi_df.select(common_genes).to_numpy().astype(np.float64)
        phi_log = np.log1p(np.maximum(phi_mat, 0.0) * 1e4)

        # Standardize across genes for specificity-weighted cosine correlation
        nh_z = (nh_expr - np.mean(nh_expr, axis=0, keepdims=True)) / np.maximum(
            np.std(nh_expr, axis=0, keepdims=True), 1e-12
        )
        phi_z = (phi_log - np.mean(phi_log, axis=0, keepdims=True)) / np.maximum(
            np.std(phi_log, axis=0, keepdims=True), 1e-12
        )
        nh_u = nh_z / np.maximum(np.linalg.norm(nh_z, axis=1, keepdims=True), 1e-12)
        phi_u = phi_z / np.maximum(np.linalg.norm(phi_z, axis=1, keepdims=True), 1e-12)
        sim_vk = nh_u @ phi_u.T

        adj = nh_sub.T.dot(nh_sub).tocsr()
        deg = np.maximum(np.asarray(adj.sum(axis=1)).ravel(), 1.0)
        sim_smooth = 0.70 * sim_vk + 0.30 * (adj.dot(sim_vk) / np.expand_dims(deg, axis=1))
        expr_best = np.argmax(sim_smooth, axis=1)

        expr_lks = [lks[int(i)] for i in expr_best.tolist()]
        expr_metas = [lk_to_meta[lk] for lk in expr_lks]
        expr_comps = [classify_meta42_compartment(ms).value for ms in expr_metas]

        pred_lks = mapping.annotations_df["mapped_cell_state"].to_list()
        pred_metas = list(mapping.nhood_meta_states)
        pred_comps = list(mapping.nhood_compartments)

        n_nh = len(pred_lks)
        lk_acc = sum(1 for a, b in zip(pred_lks, expr_lks) if a == b) / float(max(n_nh, 1))
        meta_acc = sum(1 for a, b in zip(pred_metas, expr_metas) if a == b) / float(max(n_nh, 1))
        comp_acc = sum(1 for a, b in zip(pred_comps, expr_comps) if a == b) / float(max(n_nh, 1))

        return Some(
            H5adValidationMetrics(
                accession=cohort_inputs.accession,
                n_nhoods=n_nh,
                n_shared_genes=len(common_genes),
                exact_lk_agreement=lk_acc,
                exact_meta42_agreement=meta_acc,
                exact_compartment_agreement=comp_acc,
            )
        )
    except Exception:
        return Nothing


def export_chart_to_disk(chart: AltairChart, base_path: Path) -> Result[tuple[Path, Path], str]:
    """Side-effecting boundary exporting an Altair chart to vector SVG and 300 DPI PNG."""
    try:
        base_path.parent.mkdir(parents=True, exist_ok=True)
        with alt.data_transformers.disable_max_rows():
            paths = export_altair_figure(chart, base_path, scale=2.5, max_vector_points=1500)
        return Success(paths)
    except Exception as exc:
        return Failure(f"Failed exporting figure {base_path}: {exc}")


def run_full_bridge_pipeline(
    config: BridgeConfig,
    cohort_loader: CohortArtifactLoader = load_cohort_artifacts_from_disk,
    global_loader: GlobalDeconvLoader = load_global_deconv_from_disk,
    h5ad_validator: LocalH5adValidator = validate_mapper_against_local_h5ad,
    exporter: VectorFigureExporter = export_chart_to_disk,
    persist_outputs: bool = True,
) -> Result[BridgePipelineOutputs, str]:
    """Imperative shell loading artifacts, running pure core, and persisting Parquets + SVGs."""
    global_res = global_loader(config.deconv_3b_path, config.iatlas_dir)
    match global_res:
        case Failure(err):
            return Failure(err)
        case Success(global_deconv):
            pass

    cohorts = resolve_target_cohorts(config.cohort_filter)
    cohort_bundles: list[CohortInputsBundle] = []
    for acc in cohorts:
        c_res = cohort_loader(acc, config.milo_input_dir, config.deconv_3a_dir)
        match c_res:
            case Failure(err):
                return Failure(err)
            case Success(bundle):
                cohort_bundles.append(bundle)

    cached_ablation_maybe: Maybe[pl.DataFrame] = Nothing
    cached_ablation_path = config.output_dir / "prefilter_ablation_benchmark.parquet"
    if config.reuse_cached_ablation and cached_ablation_path.exists():
        try:
            cached_ablation_maybe = Some(pl.read_parquet(cached_ablation_path))
        except Exception:
            cached_ablation_maybe = Nothing

    core_res = run_bridge_core_pure(
        cohort_inputs_seq=tuple(cohort_bundles),
        global_deconv=global_deconv,
        config=config,
        h5ad_validator=Some(h5ad_validator),
        cached_ablation_df=cached_ablation_maybe,
    )
    match core_res:
        case Failure(err):
            return Failure(err)
        case Success(outputs):
            pass

    if persist_outputs:
        out_dir = config.output_dir
        rep_dir = out_dir / "reports"
        fig_dir = out_dir / "figures"
        out_dir.mkdir(parents=True, exist_ok=True)
        rep_dir.mkdir(parents=True, exist_ok=True)
        fig_dir.mkdir(parents=True, exist_ok=True)

        outputs.nhood_annotations_df.write_parquet(out_dir / "nhood_meta42_annotations.parquet")
        outputs.ablation_benchmark_df.write_parquet(out_dir / "prefilter_ablation_benchmark.parquet")
        outputs.meta42_concordance_df.write_parquet(
            out_dir / "meta42_pairwise_milo_vs_bulk_deconv.parquet"
        )
        outputs.hub_concordance_df.write_parquet(out_dir / "meta42_hub_concordance.parquet")
        outputs.within_cohort_pseudobulk_df.write_parquet(
            out_dir / "within_cohort_pseudobulk_pair_concordance.parquet"
        )
        outputs.mycaf_gating_regimes_df.write_parquet(
            out_dir / "singlecell_mycaf_gating_regimes.parquet"
        )
        outputs.mycaf_cohort_forest_df.write_parquet(
            out_dir / "singlecell_mycaf_cohort_forest.parquet"
        )
        (rep_dir / "prefilter_and_deconv_bridge_summary.md").write_text(
            outputs.markdown_report,
            encoding="utf-8",
        )

        if config.render_plots:
            chart_specs: tuple[tuple[str, Result[alt.HConcatChart, str]], ...] = (
                (
                    "prefilter_power_ablation",
                    build_prefilter_power_ablation_chart(outputs.ablation_benchmark_df),
                ),
                (
                    "deconv_vs_pairwise_milo_concordance",
                    build_deconv_vs_pairwise_milo_chart(
                        outputs.meta42_concordance_df,
                        outputs.mantel_whole_tumor,
                        outputs.mantel_focused8,
                    ),
                ),
                (
                    "singlecell_mycaf_stromal_gating",
                    build_singlecell_mycaf_gating_chart(
                        outputs.mycaf_gating_regimes_df,
                        outputs.mycaf_cohort_forest_df,
                    ),
                ),
                (
                    "meta42_hub_concordance",
                    build_meta42_hub_and_pseudobulk_chart(
                        outputs.hub_concordance_df,
                        outputs.within_cohort_pseudobulk_df,
                    ),
                ),
            )
            for fname, chart_res in chart_specs:
                match chart_res:
                    case Failure(err):
                        return Failure(err)
                    case Success(chart):
                        exp_res = exporter(chart, fig_dir / fname)
                        if isinstance(exp_res, Failure):
                            return Failure(exp_res.failure())

    return Success(outputs)


def parse_cli_args(argv: Sequence[str]) -> BridgeConfig:
    """Parses CLI arguments into an immutable BridgeConfig."""
    parser = argparse.ArgumentParser(
        description="Run outcome-blind prefiltering ablation & deconvolution interaction bridge across 9 cohorts."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/results/milopy_pairwise_deconv_bridge"),
        help="Destination directory for Parquet tables, Markdown summary, and vector SVG figures",
    )
    parser.add_argument(
        "--cohorts",
        type=str,
        default="",
        help="Optional comma-separated cohort subset (empty = all 9 cohorts)",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.10,
        help="Pairwise IUT FDR cutoff (default: 0.10)",
    )
    parser.add_argument(
        "--tarone-p-ceiling",
        type=float,
        default=0.005,
        help="Exact discrete Tarone minimum-achievable two-sided p-value ceiling (default: 0.005)",
    )
    parser.add_argument(
        "--n-mantel-perms",
        type=int,
        default=5000,
        help="Number of node permutations for Mantel matrix concordance tests (default: 5000)",
    )
    parser.add_argument(
        "--plot",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Render Nature Methods SVG + 300 DPI PNG figures (default: True)",
    )
    parser.add_argument(
        "--reuse-cached-ablation",
        action="store_true",
        default=False,
        help="Reuse existing T0-T3 ablation benchmark rows when present and recompute T4/T5 + figures",
    )
    ns = parser.parse_args(list(argv))
    toks = tuple(t.strip() for t in str(ns.cohorts).split(",") if t.strip())
    c_maybe: Maybe[tuple[str, ...]] = Some(toks) if toks else Nothing
    return BridgeConfig(
        output_dir=ns.output_dir,
        cohort_filter=c_maybe,
        fdr_threshold=float(ns.fdr_threshold),
        tarone_p_ceiling=float(ns.tarone_p_ceiling),
        n_mantel_perms=int(ns.n_mantel_perms),
        render_plots=bool(ns.plot),
        reuse_cached_ablation=bool(ns.reuse_cached_ablation),
    )


def run_cli(argv: Sequence[str]) -> int:
    """POSIX entry point."""
    cfg = parse_cli_args(argv)
    t_start = time.perf_counter()
    match run_full_bridge_pipeline(cfg):
        case Failure(err):
            print(f"Prefilter & deconvolution bridge failed: {err}", file=sys.stderr)
            return 1
        case Success(outputs):
            elapsed = time.perf_counter() - t_start
            print(
                f"Completed 9-cohort prefilter benchmark & deconvolution bridge in {elapsed:.1f}s -> {cfg.output_dir}",
                flush=True,
            )
            print(
                f"  Mapped neighborhoods: {outputs.nhood_annotations_df.height:,} | "
                f"Mantel whole-tumor r={outputs.mantel_whole_tumor.pearson_r:+.3f} "
                f"(p={outputs.mantel_whole_tumor.pvalue_two_sided:.4f}) | "
                f"Focused 8x8 r={outputs.mantel_focused8.pearson_r:+.3f} "
                f"(p={outputs.mantel_focused8.pvalue_two_sided:.4f})",
                flush=True,
            )
            return 0
        case _ as unreachable:
            _ = unreachable
            return 1


if __name__ == "__main__":
    sys.exit(run_cli(sys.argv[1:]))
