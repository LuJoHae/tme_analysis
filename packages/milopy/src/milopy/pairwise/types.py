"""Immutable domain models, exhaustive sum types, and typeclass protocols for milopy.pairwise.

Adheres strictly to functional Python principles:
- Frozen Pydantic data models (`frozen=True`)
- Exhaustive pattern matching via `StrEnum` and `typing.assert_never`
- Zero `None` / `Optional` types (uses `returns.maybe.Maybe`)
- Polymorphism via `typing.Protocol` typeclasses (zero class inheritance for behavior)
"""

from __future__ import annotations

from enum import StrEnum
import re
from typing import Protocol, TypeAlias, assert_never

import numpy as np
from numpy.typing import NDArray
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Nothing
from returns.result import Result
from scipy import sparse as sp

FloatArray: TypeAlias = NDArray[np.float64]
IntArray: TypeAlias = NDArray[np.int64]
BoolArray: TypeAlias = NDArray[np.bool_]
SparseMatrix: TypeAlias = (
    sp.csr_matrix
    | sp.csc_matrix
    | sp.coo_matrix
    | sp.csr_array
    | sp.csc_array
    | sp.coo_array
)


class QuadrantState(StrEnum):
    """Four 2x2 neighborhood presence/absence states for an ordered pair (u, v)."""

    BOTH_PRESENT = "++"
    FIRST_ONLY = "+-"
    SECOND_ONLY = "-+"
    BOTH_MISSING = "--"


ALL_QUADRANT_STATES: tuple[QuadrantState, ...] = (
    QuadrantState.BOTH_PRESENT,
    QuadrantState.FIRST_ONLY,
    QuadrantState.SECOND_ONLY,
    QuadrantState.BOTH_MISSING,
)


def swap_quadrant(state: QuadrantState) -> QuadrantState:
    """Maps quadrant state of (u, v) to the equivalent quadrant state of (v, u)."""
    match state:
        case QuadrantState.BOTH_PRESENT:
            return QuadrantState.BOTH_PRESENT
        case QuadrantState.BOTH_MISSING:
            return QuadrantState.BOTH_MISSING
        case QuadrantState.FIRST_ONLY:
            return QuadrantState.SECOND_ONLY
        case QuadrantState.SECOND_ONLY:
            return QuadrantState.FIRST_ONLY
        case _ as unreachable:
            assert_never(unreachable)


def adjacent_flips(state: QuadrantState) -> tuple[QuadrantState, QuadrantState]:
    """Returns (flip_v, flip_u), the two Hamming-distance-1 neighbors on {-, +}^2.

    For state (a, b):
    - `flip_v` = (a, not b): tests necessity of neighborhood v given u = a.
    - `flip_u` = (not a, b): tests necessity of neighborhood u given v = b.
    """
    match state:
        case QuadrantState.BOTH_PRESENT:
            return (QuadrantState.FIRST_ONLY, QuadrantState.SECOND_ONLY)
        case QuadrantState.FIRST_ONLY:
            return (QuadrantState.BOTH_PRESENT, QuadrantState.BOTH_MISSING)
        case QuadrantState.SECOND_ONLY:
            return (QuadrantState.BOTH_MISSING, QuadrantState.BOTH_PRESENT)
        case QuadrantState.BOTH_MISSING:
            return (QuadrantState.SECOND_ONLY, QuadrantState.FIRST_ONLY)
        case _ as unreachable:
            assert_never(unreachable)


class PresenceMode(StrEnum):
    """Estimator mode for neighborhood presence/absence and soft gates across donors."""

    HARD_DEPTH_GATED = "hard_depth_gated"
    EMPIRICAL_BAYES_POSTERIOR = "empirical_bayes_posterior"
    QUANTILE_SOFT_GATE = "quantile_soft_gate"


class PairwiseMethod(StrEnum):
    """Pairwise testing family to execute."""

    QUADRANT_2X2 = "quadrant_2x2"
    LOG_RATIO_BALANCE = "log_ratio_balance"
    CONTINUOUS_INTERACTION = "continuous_interaction"
    HYBRID_ALL = "hybrid_all"


class TmeCompartment(StrEnum):
    """Canonical broad cellular compartments of the tumor microenvironment."""

    T_NK = "T_NK"
    B_LINEAGE = "B_LINEAGE"
    PLASMA = "PLASMA"
    MYELOID_MAST = "MYELOID_MAST"
    STROMAL_MESENCHYMAL = "STROMAL_MESENCHYMAL"
    ENDOTHELIAL_VASCULAR = "ENDOTHELIAL_VASCULAR"
    MALIGNANT_EPITHELIAL = "MALIGNANT_EPITHELIAL"
    OTHER = "OTHER"


_META42_COMPARTMENT_MAP: dict[str, TmeCompartment] = {
    "M01": TmeCompartment.T_NK,
    "M02": TmeCompartment.STROMAL_MESENCHYMAL,
    "M03": TmeCompartment.T_NK,
    "M04": TmeCompartment.T_NK,
    "M05": TmeCompartment.T_NK,
    "M06": TmeCompartment.T_NK,
    "M07": TmeCompartment.OTHER,
    "M08": TmeCompartment.T_NK,
    "M09": TmeCompartment.MYELOID_MAST,
    "M10": TmeCompartment.MYELOID_MAST,
    "M11": TmeCompartment.MALIGNANT_EPITHELIAL,
    "M12": TmeCompartment.MALIGNANT_EPITHELIAL,
    "M13": TmeCompartment.T_NK,
    "M14": TmeCompartment.MYELOID_MAST,
    "M15": TmeCompartment.T_NK,
    "M16": TmeCompartment.ENDOTHELIAL_VASCULAR,
    "M17": TmeCompartment.STROMAL_MESENCHYMAL,
    "M18": TmeCompartment.MALIGNANT_EPITHELIAL,
    "M19": TmeCompartment.MYELOID_MAST,
    "M20": TmeCompartment.T_NK,
    "M21": TmeCompartment.T_NK,
    "M22": TmeCompartment.MALIGNANT_EPITHELIAL,
    "M23": TmeCompartment.STROMAL_MESENCHYMAL,
    "M24": TmeCompartment.T_NK,
    "M25": TmeCompartment.T_NK,
    "M26": TmeCompartment.MALIGNANT_EPITHELIAL,
    "M27": TmeCompartment.MALIGNANT_EPITHELIAL,
    "M28": TmeCompartment.MYELOID_MAST,
    "M29": TmeCompartment.PLASMA,
    "M30": TmeCompartment.B_LINEAGE,
    "M31": TmeCompartment.MYELOID_MAST,
    "M32": TmeCompartment.MYELOID_MAST,
    "M33": TmeCompartment.T_NK,
    "M34": TmeCompartment.T_NK,
    "M35": TmeCompartment.MYELOID_MAST,
    "M36": TmeCompartment.T_NK,
    "M37": TmeCompartment.B_LINEAGE,
    "M38": TmeCompartment.MYELOID_MAST,
    "M39": TmeCompartment.MYELOID_MAST,
    "M40": TmeCompartment.T_NK,
    "M41": TmeCompartment.MYELOID_MAST,
    "M42": TmeCompartment.MYELOID_MAST,
}

_META42_CODE_RE: re.Pattern[str] = re.compile(r"\b(M(?:0[1-9]|[1-3][0-9]|4[0-2]))\b")


def classify_tme_compartment(label: str) -> TmeCompartment:
    """Deterministically normalizes Cell Ontology strings and Meta42 programs into broad TME compartments."""
    raw = str(label).strip()
    if not raw or raw.isdigit():
        return TmeCompartment.OTHER

    upper = raw.upper()
    m_match = _META42_CODE_RE.search(upper)
    if m_match:
        code = m_match.group(1)
        if code in _META42_COMPARTMENT_MAP:
            return _META42_COMPARTMENT_MAP[code]

    low = raw.lower()
    if low in ("t_nk",):
        return TmeCompartment.T_NK
    if low in ("b_lineage", "b_cell"):
        return TmeCompartment.B_LINEAGE
    if low in ("plasma",):
        return TmeCompartment.PLASMA
    if low in ("myeloid_mast",):
        return TmeCompartment.MYELOID_MAST
    if low in ("stromal_mesenchymal", "stromal_caf"):
        return TmeCompartment.STROMAL_MESENCHYMAL
    if low in ("endothelial_vascular", "endothelial"):
        return TmeCompartment.ENDOTHELIAL_VASCULAR
    if low in ("malignant_epithelial", "epithelial_tumor"):
        return TmeCompartment.MALIGNANT_EPITHELIAL
    if low in (
        "unknown",
        "unspecified",
        "unannotated",
        "unassigned",
        "empty",
        "other",
        "other_unspecified",
        "platelet",
    ) or any(tok in low for tok in ("erythro", "schwann", "platelet")):
        return TmeCompartment.OTHER

    # Check plasmacytoid DC before plasma cell so pDCs map strictly to MYELOID_MAST
    if "plasmacytoid" in low or re.search(r"\bpdc\b", low):
        return TmeCompartment.MYELOID_MAST
    if "plasma" in low or "plasmablast" in low:
        return TmeCompartment.PLASMA

    # Check B-cell lineage before generic 'lymphocyte' so 'lymphocyte of B lineage' -> B_LINEAGE
    if (
        "b lineage" in low
        or "tls b" in low
        or re.search(r"\bb[\s_-]?cells?\b", low)
    ):
        return TmeCompartment.B_LINEAGE

    if any(
        tok in low
        for tok in (
            "myeloid",
            "macrophage",
            "monocyte",
            "dendritic",
            "mregdc",
            "mast",
            "neutrophil",
            "granulocyte",
            "mdsc",
        )
    ) or re.search(r"\b(?:cdc|tam|tan)s?\b", low):
        return TmeCompartment.MYELOID_MAST

    if "endothelial" in low or "endo." in low:
        return TmeCompartment.ENDOTHELIAL_VASCULAR

    if any(
        tok in low
        for tok in (
            "fibroblast",
            "myofibroblast",
            "mycaf",
            "icaf",
            "stroma",
            "pericyte",
            "perivascular",
            "mesenchymal",
        )
    ) or re.search(r"\bcafs?\b", low):
        return TmeCompartment.STROMAL_MESENCHYMAL

    if any(
        tok in low
        for tok in (
            "malignant",
            "tumor",
            "cancer",
            "abnormal cell",
            "epithelial",
            "enterocyte",
            "colonocyte",
            "goblet",
            "tuft",
            "crypt",
            "transit amplifying",
            "enteroendocrine",
            "melanoc",
            "ductal",
            "neuroend",
        )
    ):
        return TmeCompartment.MALIGNANT_EPITHELIAL

    if any(
        tok in low
        for tok in (
            "t follicular",
            "t-helper",
            "regulatory t",
            "natural killer",
            "innate lymphoid",
            "lymphocyte",
            "cd8",
            "cd4",
            "treg",
            "tfh",
            "trm",
            "tex",
            "teff",
            "tem",
            "tcm",
        )
    ) or re.search(r"\b(?:t[\s_-]?cells?|nk|ilc|th1|th2|th17)\b", low):
        return TmeCompartment.T_NK

    return TmeCompartment.OTHER


class PairScope(StrEnum):
    """Lineage relationship between neighborhoods u and v."""

    CROSS_COMPARTMENT = "cross_compartment"
    CROSS_CELL_TYPE = "cross_cell_type"
    WITHIN_CELL_TYPE = "within_cell_type"
    UNANNOTATED = "unannotated"


class PairScopeFilter(StrEnum):
    """Filter restricting candidate neighborhood pairs by cell-type or compartment relationship."""

    ALL = "all"
    CROSS_CELL_TYPE_ONLY = "cross_cell_type_only"
    CROSS_COMPARTMENT_ONLY = "cross_compartment_only"
    WITHIN_CELL_TYPE_ONLY = "within_cell_type_only"


class EnrichmentDirection(StrEnum):
    """Clinical response arm enriched for a neighborhood pair regime."""

    RESPONDER = "responder"
    NON_RESPONDER = "non-responder"


class ContinuousAxis(StrEnum):
    """Continuous pairwise axis generalizing +-/-+ balance and ++/-- coupling."""

    LOG_RATIO_BALANCE = "log_ratio_balance"
    CO_ABUNDANCE_PRODUCT = "co_abundance_product"


class PrefilterMode(StrEnum):
    """Outcome-blind prefiltering and graph coarse-graining mode prior to pairwise testing."""

    NONE = "none"
    TARONE_ONLY = "tarone_only"
    MARGINAL_ENTROPY_VARIANCE = "marginal_entropy_variance"
    GRAPH_SUPER_NHOODS = "graph_super_nhoods"
    COMPOSITE_SUPER_TARONE = "composite_super_tarone"


class PrefilterConfig(BaseModel):
    """Immutable configuration for strictly outcome-blind neighborhood prefiltering & pooling."""

    model_config = ConfigDict(frozen=True)

    mode: PrefilterMode = Field(
        default=PrefilterMode.NONE,
        description="Outcome-blind prefiltering or super-neighborhood pooling mode",
    )
    min_marginal_prevalence: float = Field(
        default=0.10,
        ge=0.0,
        lt=0.5,
        description="Minimum unlabelled donor presence fraction pi_v to retain neighborhood v",
    )
    max_marginal_prevalence: float = Field(
        default=0.90,
        gt=0.5,
        le=1.0,
        description="Maximum unlabelled donor presence fraction pi_v to retain neighborhood v",
    )
    min_variance_quantile: float = Field(
        default=0.25,
        ge=0.0,
        lt=1.0,
        description="Minimum across-donor log2_cpm total variance quantile to retain neighborhood v",
    )
    max_active_nhoods_cap: int = Field(
        default=0,
        ge=0,
        description="Optional cap on top outcome-blind entropy/variance neighborhoods (0 = no cap)",
    )
    super_nhood_min_jaccard: float = Field(
        default=0.10,
        gt=0.0,
        le=1.0,
        description="Minimum within-lineage Jaccard overlap to link neighborhoods for Leiden pooling",
    )
    super_nhood_resolution: float = Field(
        default=1.0,
        gt=0.0,
        description="Leiden CPM resolution scaling factor for within-lineage super-neighborhood pooling",
    )
    apply_tarone_iut_pruning: bool = Field(
        default=False,
        description="Prune discrete pair-quadrants whose exact outcome-blind minimum IUT p-value > tarone_p_ceiling",
    )
    tarone_p_ceiling: float = Field(
        default=0.01,
        gt=0.0,
        le=1.0,
        description="Outcome-blind minimum attainable IUT p-value ceiling for Tarone testability pruning",
    )


class OverlapConfig(BaseModel):
    """Immutable configuration for kNN neighborhood overlap filtering and spatial weighting."""

    model_config = ConfigDict(frozen=True)

    max_jaccard_overlap: float = Field(
        default=0.0,
        ge=0.0,
        le=1.0,
        description="Maximum allowed cell-sharing Jaccard index J(u, v) for tested pairs (0.0 = strictly disjoint)",
    )
    scope_filter: PairScopeFilter = Field(
        default=PairScopeFilter.ALL,
        description="Whether to test all disjoint pairs, cross-cell-type only, or within-cell-type only",
    )


class PresenceConfig(BaseModel):
    """Immutable configuration for depth-calibrated neighborhood presence estimation."""

    model_config = ConfigDict(frozen=True)

    mode: PresenceMode = Field(
        default=PresenceMode.HARD_DEPTH_GATED,
        description="Presence representation mode across patients",
    )
    min_cell_count: int = Field(
        default=1,
        ge=1,
        description="Minimum cells from a donor in neighborhood v to declare hard presence",
    )
    min_fraction: float = Field(
        default=1e-4,
        ge=0.0,
        le=1.0,
        description="Minimum donor cell fraction C_{v,j} / L_j required for hard presence",
    )
    min_active_donors: int = Field(
        default=2,
        ge=1,
        description="Minimum donors with positive cells for neighborhood v to enter pairwise screening",
    )
    use_detection_power_weights: bool = Field(
        default=True,
        description="Down-weight unobserved zeros in shallow biopsies by Poisson detection power",
    )
    prior_cpm_pseudocount: float = Field(
        default=1.0,
        gt=0.0,
        description="Library-size-scaled pseudocount (in counts-per-thousand/CPM units) for log-ratio balances",
    )
    eb_background_rate: float = Field(
        default=1e-4,
        gt=0.0,
        le=0.5,
        description="Background ambient cell rate per library unit for empirical-Bayes hurdle model",
    )
    residualize_library_depth: bool = Field(
        default=False,
        description="Orthogonalize continuous log2_cpm and arcsin_sqrt_prop features against log10(library_sizes)",
    )
    mask_unobserved_compartments: bool = Field(
        default=False,
        description="Mask non-OTHER compartments absent from FACS-sorted cohorts (< min_cohort_compartment_share)",
    )
    min_cohort_compartment_share: float = Field(
        default=0.015,
        ge=0.0,
        le=0.25,
        description="Minimum cohort cell share in a TME compartment to treat that compartment as observable",
    )


class ScreenConfig(BaseModel):
    """Immutable configuration for Stage 1 chunked BLAS pairwise screening."""

    model_config = ConfigDict(frozen=True)

    chunk_size: int = Field(
        default=512,
        ge=16,
        description="Row block size for memory-bounded BLAS upper-triangular matrix screening",
    )
    fdr_threshold: float = Field(
        default=0.10,
        gt=0.0,
        le=1.0,
        description="Benjamini-Hochberg FDR cutoff for significance",
    )
    min_donors_per_quadrant: float = Field(
        default=3.0,
        ge=1.0,
        description="Minimum effective donors in the target quadrant within the enriched response arm",
    )
    min_effective_samples: float = Field(
        default=6.0,
        ge=2.0,
        description="Minimum total effective observable donors across both arms for a pair",
    )
    pseudocount: float = Field(
        default=0.5,
        gt=0.0,
        description="Haldane-Anscombe continuity pseudocount for 2x2 log-odds estimation",
    )
    use_spatial_fdr_weights: bool = Field(
        default=True,
        description="Weight Benjamini-Hochberg FDR by 2D inverse-overlap product weights w_u * w_v",
    )
    n_permutations: int = Field(
        default=0,
        ge=0,
        description="Number of stratified donor-label permutations for empirical null p-values (0 disables)",
    )
    permutation_seed: int = Field(
        default=42,
        description="Random seed for reproducible stratified donor permutations",
    )
    max_retained_candidates: int = Field(
        default=50000,
        ge=100,
        description="Maximum top pair-quadrant rows retained in memory after FDR calculation",
    )
    apply_tarone_iut_pruning: bool = Field(
        default=False,
        description="Prune discrete pair-quadrants whose exact outcome-blind minimum IUT p-value > tarone_p_ceiling",
    )
    tarone_p_ceiling: float = Field(
        default=0.01,
        gt=0.0,
        le=1.0,
        description="Outcome-blind minimum attainable IUT p-value ceiling for Tarone testability pruning",
    )
    dmpp_null_draws: int = Field(
        default=0,
        ge=0,
        le=500,
        description="Explicit decoupled partner-permutation count B for Knockoff+ FDR (0 uses adaptive default)",
    )


class RefinementConfig(BaseModel):
    """Immutable configuration for Stage 2 Firth penalized logistic LRT refinement."""

    model_config = ConfigDict(frozen=True)

    enabled: bool = Field(
        default=True,
        description="Refine top screened candidates with covariate-adjusted Firth logistic LRT",
    )
    max_refine_pairs: int = Field(
        default=250,
        ge=0,
        description="Maximum top IUT-ranked candidate rows refined via Firth logistic regression",
    )
    adjust_library_size: bool = Field(
        default=True,
        description="Include centered log(effective_library_size) as covariate in Firth GLM",
    )
    max_iter: int = Field(
        default=60,
        ge=5,
        description="Maximum Newton-Raphson iterations for Firth logistic regression",
    )


class PairModuleConfig(BaseModel):
    """Immutable configuration for grouping redundant neighborhood-pair hits into modules."""

    model_config = ConfigDict(frozen=True)

    enabled: bool = Field(
        default=True,
        description="Cluster overlapping neighborhood-pair hits into non-redundant ecosystem modules",
    )
    min_endpoint_jaccard: float = Field(
        default=0.05,
        gt=0.0,
        le=1.0,
        description="Minimum Jaccard overlap between u1-u2 and v1-v2 to link two pair hits into a module",
    )
    leiden_resolution: float = Field(
        default=1.0,
        gt=0.0,
        description="Resolution parameter for igraph Leiden clustering of dense pair-hit components",
    )


class MixingQCConfig(BaseModel):
    """Immutable outcome-blind donor and batch/sub-study mixing QC configuration."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    enabled: bool = Field(
        default=False,
        description="Gate active neighborhoods on outcome-blind donor-mixing LISI & dominance prior to pair screening",
    )
    min_distinct_donors: int = Field(
        default=3,
        ge=1,
        description="Minimum distinct donors contributing >= 1 cell to neighborhood v",
    )
    min_effective_donors_lisi: float = Field(
        default=2.5,
        ge=1.0,
        description="Minimum Inverse-Simpson effective donors LISI_pt = 1 / sum_j p_{v,j}^2",
    )
    max_single_donor_dominance: float = Field(
        default=0.50,
        gt=0.0,
        le=1.0,
        description="Maximum fraction of neighborhood cells derived from a single donor",
    )
    max_batch_fraction: Maybe[float] = Field(
        default=Nothing,
        description="Optional ceiling on single-batch/sub-study cell fraction when batch labels exist",
    )
    adaptive_small_cohort_scaling: bool = Field(
        default=True,
        description="Proportionally relax donor-mixing thresholds when cohort total donors J < 12",
    )


class PairwiseConfig(BaseModel):
    """Top-level immutable configuration for milopy pairwise neighborhood analysis."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    method: PairwiseMethod = Field(default=PairwiseMethod.HYBRID_ALL)
    mixing_qc: MixingQCConfig = Field(default=MixingQCConfig())
    prefilter: PrefilterConfig = Field(default=PrefilterConfig())
    overlap: OverlapConfig = Field(default=OverlapConfig())
    presence: PresenceConfig = Field(default=PresenceConfig())
    screen: ScreenConfig = Field(default=ScreenConfig())
    refinement: RefinementConfig = Field(default=RefinementConfig())
    modules: PairModuleConfig = Field(default=PairModuleConfig())


class CoarseGrainedNhoodData(BaseModel):
    """Immutable lineage-preserving super-neighborhood representation with exact traceability."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    nhoods_super: sp.csc_matrix
    counts_super: FloatArray
    super_nhood_labels: tuple[str, ...]
    member_nhood_indices: tuple[tuple[int, ...], ...]
    nhood_to_super_id: IntArray


class OverlapGraph(BaseModel):
    """Immutable kNN neighborhood overlap topology and spatial uniqueness weights."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    n_nhoods: int
    nhood_sizes: FloatArray
    jaccard_matrix: sp.csr_matrix
    spatial_weights_1d: FloatArray
    dominant_cell_types: tuple[str, ...]
    cell_type_purities: tuple[float, ...]
    has_cell_types: bool


class PresenceRepresentation(BaseModel):
    """Immutable depth-calibrated presence/absence and continuous abundance matrices."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    active_nhood_indices: IntArray
    presence_prob: FloatArray
    absence_prob: FloatArray
    observability: FloatArray
    raw_binary_presence: BoolArray
    arcsin_sqrt_prop: FloatArray
    log2_cpm: FloatArray
    effective_library_sizes: FloatArray


class PairwiseAnalysisResult(BaseModel):
    """Immutable result container for pairwise neighborhood differential association."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    pairs_df: pl.DataFrame
    continuous_pairs_df: pl.DataFrame
    modules_df: pl.DataFrame
    n_total_nhoods: int
    n_active_nhoods: int
    n_evaluated_pairs: int
    n_significant_quadrant_pairs: int
    n_significant_continuous_pairs: int
    n_modules: int


class PresenceEstimator(Protocol):
    """Typeclass protocol for mapping count matrices to depth-calibrated presence representations."""

    def __call__(
        self,
        count_matrix: FloatArray,
        library_sizes: Maybe[FloatArray],
        config: PresenceConfig,
    ) -> Result[PresenceRepresentation, str]: ...


class CandidatePairRefiner(Protocol):
    """Typeclass protocol for Stage 2 covariate-adjusted parametric refinement of candidate pairs."""

    def __call__(
        self,
        candidates_df: pl.DataFrame,
        presence: PresenceRepresentation,
        response_binary: FloatArray,
        strata_codes: Maybe[IntArray],
        config: RefinementConfig,
        fdr_threshold: float,
    ) -> Result[pl.DataFrame, str]: ...
