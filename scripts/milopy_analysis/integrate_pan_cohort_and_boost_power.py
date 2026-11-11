#!/usr/bin/env python3
"""Pan-Cohort Single-Cell Integration (`harmonypy`) & Five-Lever Pairwise Power Pipeline.

Combines all nine clinical single-cell ICB cohorts (`371,568` cells, `613` biopsy evaluations:
`304 Responders / 309 Non-Responders`, `16+` cohort/sub-study strata) into a unified single-cell
atlas, benchmarks `harmonypy` latent integration (`I0` Unintegrated PCA vs. `I1` Cohort vs. `I2`
Cohort + Sub-Study vs. `I3` Cohort + Sub-Study + Donor), and evaluates a 7-tier outcome-blind
statistical power progression (`P0` -> `P6`) side-by-side against within-stratum permuted-response
null controls.

Architecture:
- Pure Functional Core (`ConfigDict(frozen=True)`, `StrEnum` + `assert_never`, `Protocol`
  typeclasses, `Result` / `Maybe` monads, Polars lazy pipelines, Altair vector composites).
- Imperative Shell (Dual-mode `.h5ad` / outcome-blind `4,529`-gene feature loader, Inkscape-layered
  hybrid 300 DPI SVG/PNG exporter, and CLI entrypoint).
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from enum import StrEnum
import hashlib
import html
from pathlib import Path
import sys
import time
from typing import Literal, Protocol, TypeAlias, assert_never, cast
import xml.etree.ElementTree as ET

import altair as alt
import harmonypy as hm  # type: ignore[import-untyped]
import numpy as np
import pandas as pd  # Required exclusively at the harmonypy.run_harmony third-party boundary
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import sparse as sp
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist
from scipy.stats import chi2, norm, pearsonr, spearmanr, t as student_t
from sklearn.decomposition import PCA  # type: ignore[import-untyped]
from sklearn.metrics import silhouette_score  # type: ignore[import-untyped]
import umap  # type: ignore[import-untyped]
import vl_convert as vlc

from milopy.pairwise import (
    ContinuousAxis,
    MixingQCConfig,
    OverlapConfig,
    PairModuleConfig,
    PairScope,
    PairScopeFilter,
    PairwiseAnalysisResult,
    PairwiseConfig,
    PairwiseMethod,
    PairwiseUmapBundle,
    PrefilterConfig,
    PrefilterMode,
    PresenceConfig,
    PresenceMode,
    PresenceRepresentation,
    RefinementConfig,
    ScreenConfig,
    SpotlightPairMetadata,
    TmeCompartment,
    UmapPairVizConfig,
    analyze_nhood_pairs,
    build_compartment_observability_mask,
    build_pairwise_umap_composite,
    calibrate_candidates_stratified_iut_permutation,
    classify_tme_compartment,
    coarse_grain_nhoods_by_graph,
    combine_pair_regimes_acat,
    compute_nhood_mixing_metrics,
    compute_weighted_bh_fdr,
    depth_residualize_continuous_features,
    estimate_nhood_presence,
    extract_inferred_substudy_strata,
    meta_analyze_cohort_pairs_iut,
    prepare_pairwise_umap_bundle,
)
from milopy.pairwise.types import (
    BoolArray,
    EnrichmentDirection,
    FloatArray,
    IntArray,
    QuadrantState,
    SparseMatrix,
)
from milopy.pairwise.umap_viz import (
    ALL_VISUAL_REGIMES,
    EndpointPolarity,
    PairwiseVisualRegime,
    build_anchor_partner_umap_panel,
    build_denovo_hub_umap_panel,
    build_lead_pair_spotlight_panel,
    classify_continuous_visual_regime_and_polarity,
    classify_quadrant_visual_regime_and_polarity,
    compute_nhood_umap_medoids,
    compute_umap_bounding_domain,
    discretize_outward_bezier_arc,
    regime_color_hex,
    wilson_binomial_ci_pct,
)
from plotting_utils import (  # type: ignore[import-untyped]
    COHORT_METADATA,
    COHORT_PANELS,
    OKABE_BLACK,
    OKABE_BLUE,
    OKABE_BLUISH_GREEN,
    OKABE_GREY,
    OKABE_ORANGE,
    OKABE_REDDISH_PURPLE,
    OKABE_SKY_BLUE,
    OKABE_VERMILION,
    OKABE_YELLOW,
    apply_nature_methods_theme,
    create_nature_axis,
    get_single_cell_cohort_display_label,
)
from plotting_utils.svg import (  # type: ignore[import-untyped]
    INKSCAPE_NS,
    SODIPODI_NS,
    SVG_NS,
    XLINK_NS,
    rasterize_svg_datapoints,
    rasterize_svg_to_png,
)

AltairComposite: TypeAlias = alt.Chart | alt.LayerChart | alt.HConcatChart | alt.VConcatChart
alt.data_transformers.disable_max_rows()


# =============================================================================
# 1. Exhaustive Sum Types, Frozen Domain Models & Typeclass Protocols
# =============================================================================


class HarmonyRegime(StrEnum):
    """Four latent embedding integration regimes benchmarked across all 9 cohorts."""

    I0_UNINTEGRATED_PCA = "I0_unintegrated_pca"
    I1_HARMONY_COHORT = "I1_harmony_cohort"
    I2_HARMONY_COHORT_SUBSTUDY = "I2_harmony_cohort_substudy"
    I3_HARMONY_COHORT_SUBSTUDY_DONOR = "I3_harmony_cohort_substudy_donor"


ALL_HARMONY_REGIMES: tuple[HarmonyRegime, ...] = (
    HarmonyRegime.I0_UNINTEGRATED_PCA,
    HarmonyRegime.I1_HARMONY_COHORT,
    HarmonyRegime.I2_HARMONY_COHORT_SUBSTUDY,
    HarmonyRegime.I3_HARMONY_COHORT_SUBSTUDY_DONOR,
)


def harmony_regime_display_label(regime: HarmonyRegime) -> str:
    """Human-readable display name for each Harmony benchmark regime."""
    match regime:
        case HarmonyRegime.I0_UNINTEGRATED_PCA:
            return "I0: Unintegrated PCA"
        case HarmonyRegime.I1_HARMONY_COHORT:
            return "I1: Harmony (Cohort)"
        case HarmonyRegime.I2_HARMONY_COHORT_SUBSTUDY:
            return "I2: Harmony (Cohort + Sub-Study)"
        case HarmonyRegime.I3_HARMONY_COHORT_SUBSTUDY_DONOR:
            return "I3: Harmony (Cohort + Sub-Study + Donor)"
        case _ as unreachable:
            assert_never(unreachable)


class PowerTier(StrEnum):
    """Seven-tier statistical power progression evaluated on real and permuted outcomes."""

    P0_WITHIN_COHORT_BASELINE = "P0_within_cohort_baseline"
    P1_WITHIN_COHORT_EB_DEPTH_B64 = "P1_within_cohort_eb_depth_b64"
    P2_WITHIN_COHORT_ACAT_PRIOR = "P2_within_cohort_acat_prior"
    P3_TWO_STAGE_META9_IUT = "P3_two_stage_meta9_iut"
    P4_PAN_COHORT_HARMONY_NAIVE = "P4_pan_cohort_harmony_naive"
    P5_PAN_COHORT_HARMONY_OBS_MASKED = "P5_pan_cohort_harmony_obs_masked"
    P6_PAN_COHORT_HARMONY_FULL_POWER = "P6_pan_cohort_harmony_full_power"


ALL_POWER_TIERS: tuple[PowerTier, ...] = (
    PowerTier.P0_WITHIN_COHORT_BASELINE,
    PowerTier.P1_WITHIN_COHORT_EB_DEPTH_B64,
    PowerTier.P2_WITHIN_COHORT_ACAT_PRIOR,
    PowerTier.P3_TWO_STAGE_META9_IUT,
    PowerTier.P4_PAN_COHORT_HARMONY_NAIVE,
    PowerTier.P5_PAN_COHORT_HARMONY_OBS_MASKED,
    PowerTier.P6_PAN_COHORT_HARMONY_FULL_POWER,
)


def power_tier_display_label(tier: PowerTier) -> str:
    """Short publication axis label for each statistical power tier."""
    match tier:
        case PowerTier.P0_WITHIN_COHORT_BASELINE:
            return "P0: Single-Cohort Baseline (Hard + B=12)"
        case PowerTier.P1_WITHIN_COHORT_EB_DEPTH_B64:
            return "P1: + EB Soft Presence + Depth + B=64"
        case PowerTier.P2_WITHIN_COHORT_ACAT_PRIOR:
            return "P2: + Omnibus ACAT + Pan-iAtlas Prior"
        case PowerTier.P3_TWO_STAGE_META9_IUT:
            return "P3: 2-Stage 9-Cohort Random-Effects Meta-IUT"
        case PowerTier.P4_PAN_COHORT_HARMONY_NAIVE:
            return "P4: 1-Stage Harmony Pooled (Unmasked)"
        case PowerTier.P5_PAN_COHORT_HARMONY_OBS_MASKED:
            return "P5: + Assay Observability Mask (CD45+ Safe)"
        case PowerTier.P6_PAN_COHORT_HARMONY_FULL_POWER:
            return "P6: Unified Harmony + All 5 Power Levers"
        case _ as unreachable:
            assert_never(unreachable)


class FeatureLoadMode(StrEnum):
    """Source mode used to construct outcome-blind single-cell latent features."""

    DIRECT_H5AD = "direct_h5ad"
    OUTCOME_BLIND_MANIFOLD = "outcome_blind_manifold"


class TmeSuperGroup(StrEnum):
    """Top-level biological super-grouping of TME compartments."""

    IMMUNE_CELLS = "immune_cells"
    STROMAL_ENDOTHELIAL = "stromal_endothelial"
    MALIGNANT_EPITHELIAL = "malignant_epithelial"


def classify_compartment_supergroup(comp: TmeCompartment) -> Maybe[TmeSuperGroup]:
    """Maps a 7-way `TmeCompartment` into its top-level lineage super-group."""
    match comp:
        case (
            TmeCompartment.T_NK
            | TmeCompartment.B_LINEAGE
            | TmeCompartment.PLASMA
            | TmeCompartment.MYELOID_MAST
        ):
            return Some(TmeSuperGroup.IMMUNE_CELLS)
        case TmeCompartment.STROMAL_MESENCHYMAL | TmeCompartment.ENDOTHELIAL_VASCULAR:
            return Some(TmeSuperGroup.STROMAL_ENDOTHELIAL)
        case TmeCompartment.MALIGNANT_EPITHELIAL:
            return Some(TmeSuperGroup.MALIGNANT_EPITHELIAL)
        case TmeCompartment.OTHER:
            return Nothing
        case _ as unreachable:
            assert_never(unreachable)


def supergroup_display_label(supergroup: TmeSuperGroup) -> str:
    """Publication display label for a top-level TME super-group."""
    match supergroup:
        case TmeSuperGroup.IMMUNE_CELLS:
            return "Immune Cells"
        case TmeSuperGroup.STROMAL_ENDOTHELIAL:
            return "Stromal / Endothelial"
        case TmeSuperGroup.MALIGNANT_EPITHELIAL:
            return "Malignant / Epithelial"
        case _ as unreachable:
            assert_never(unreachable)


class ChordSignificanceTier(StrEnum):
    """Dual-tier statistical significance styling for cross-compartment UMAP chords."""

    VERIFIED_ACAT_FDR = "ACAT FDR < 0.10 (Verified)"
    NOMINAL_PERM_P05 = "Stratified Perm. p < 0.05 (*)"


class BroadCompartmentPairSpec(BaseModel):
    """Canonical specification for an unordered pair of broad host TME compartments."""

    model_config = ConfigDict(frozen=True)

    comp_first: TmeCompartment
    comp_second: TmeCompartment
    slug: str
    display_first: str
    display_second: str
    display_title: str
    macro_category: str
    color_first: str
    color_second: str
    n_fdr_pairs: int
    n_perm_only_pairs: int


class CompartmentPairUmapCellBundle(BaseModel):
    """Pure DataFrame bundle backing a single broad cell-type pair UMAP panel."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    spec: BroadCompartmentPairSpec
    umap_bounds: tuple[float, float, float, float]
    cells_df: pl.DataFrame
    tier1_arcs_df: pl.DataFrame
    tier2_arcs_df: pl.DataFrame
    endpoints_df: pl.DataFrame
    medoid_labels_df: pl.DataFrame
    lead_pair_summary: str


class HeatmapClusteredState(BaseModel):
    """Ordered cell-state entry after hierarchical clustering within TME compartments."""

    model_config = ConfigDict(frozen=True)

    order_idx: int
    state_code: str
    readable_label: str
    compartment: TmeCompartment
    compartment_display: str
    compartment_short: str
    compartment_color: str
    supergroup: TmeSuperGroup
    supergroup_display: str
    peak_verified_abs_z: float
    peak_any_abs_z: float
    lead_min_loco_z: float
    has_verified_pair: bool


class PanCohortCliConfig(BaseModel):
    """Immutable CLI configuration for pan-cohort integration and statistical power benchmarking."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    milo_response_dir: Path = Field(default=Path("output/results/milopy_response_analysis"))
    deconv_3a_dir: Path = Field(
        default=Path("output/manuscript_pipeline/3a_cross_cohort_pseudobulk_deconv_benchmark")
    )
    phi_meta_path: Path = Field(
        default=Path(
            "output/manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark/"
            "unified_reference_phi_meta.parquet"
        )
    )
    meta42_annotations_path: Path = Field(
        default=Path("output/results/milopy_pairwise_deconv_bridge/nhood_meta42_annotations.parquet")
    )
    iatlas_interactions_path: Path = Field(
        default=Path(
            "output/results/iatlas_instaprism_benchmark/"
            "consensus_metastate_all42_interactions.parquet"
        )
    )
    h5ad_dir: Path = Field(default=Path("data/preprocessed"))
    output_dir: Path = Field(
        default=Path("output/results/milopy_pan_cohort_integration_and_power")
    )
    max_cells_per_cohort: int = Field(default=4000, ge=0, le=200000)
    n_pcs: int = Field(default=30, ge=10, le=60)
    knn_k: int = Field(default=30, ge=10, le=100)
    milo_prop: float = Field(default=0.06, gt=0.01, le=0.30)
    super_resolution: float = Field(default=0.15, gt=0.02, le=3.0)
    super_min_jaccard: float = Field(default=0.01, gt=0.0, le=0.50)
    fdr_threshold: float = Field(default=0.10, gt=0.0, le=0.50)
    dmpp_null_draws: int = Field(default=64, ge=16, le=300)
    n_permutations: int = Field(default=1000, ge=100, le=5000)
    n_null_calibration_replicates: int = Field(default=5, ge=1, le=25)
    seed: int = Field(default=42)
    refresh_verified_and_heatmap: bool = Field(default=False)
    export_heatmaps_only: bool = Field(default=False)
    export_compartment_umaps_only: bool = Field(default=False)


class CohortSingleCellSlice(BaseModel):
    """Immutable single-cohort evaluation and neighborhood slice."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    accession: str
    cohort_display: str
    cancer_type: str
    n_cells_total: int
    n_nhoods_local: int
    nhoods_local: SparseMatrix
    count_matrix_local: FloatArray
    library_sizes: FloatArray
    eval_uids: tuple[str, ...]
    substudy_ids: tuple[str, ...]
    response_labels: tuple[str, ...]
    nhood_meta_codes: tuple[str, ...]
    nhood_readable_labels: tuple[str, ...]
    cell_local_eval: IntArray
    full_cell_pca: FloatArray
    full_cell_eval_global_idx: IntArray
    full_cell_meta_codes: tuple[str, ...]


class PanCohortFeatureBundle(BaseModel):
    """Immutable 9-cohort feature bundle combining anchor PCA cells and all 371,568 full cells."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    load_mode: FeatureLoadMode
    cohorts: tuple[CohortSingleCellSlice, ...]
    anchor_x_pca: FloatArray
    anchor_cohort_ids: tuple[str, ...]
    anchor_cancer_types: tuple[str, ...]
    anchor_substudy_ids: tuple[str, ...]
    anchor_sample_uids: tuple[str, ...]
    anchor_meta_codes: tuple[str, ...]
    anchor_compartments: tuple[str, ...]
    eval_sample_uids: tuple[str, ...]
    eval_cohort_ids: tuple[str, ...]
    eval_cancer_types: tuple[str, ...]
    eval_substudy_strata: tuple[str, ...]
    eval_response_labels: tuple[str, ...]
    eval_library_sizes: FloatArray
    code_to_readable: Mapping[str, str]
    n_total_cells_all_cohorts: int


class EmbeddingIntegrationMetrics(BaseModel):
    """Cell-level batch integration and biological conservation metrics on a latent space."""

    model_config = ConfigDict(frozen=True)

    regime: HarmonyRegime
    cohort_ilisi: float
    substudy_ilisi: float
    donor_ilisi: float
    compartment_clisi: float
    meta_state_clisi: float
    batch_asw: float
    compartment_asw: float
    kbet_acceptance_rate: float


class PanCohortMiloAtlas(BaseModel):
    """Milo neighborhood atlas constructed on a pan-cohort embedding and projected to all cells."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    regime: HarmonyRegime
    nhoods_anchor: SparseMatrix
    anchor_eval_counts: FloatArray
    counts_full_pan: FloatArray
    nhood_meta_codes: tuple[str, ...]
    nhood_compartments: tuple[str, ...]
    nhood_readable_labels: tuple[str, ...]
    nhood_mixing_df: pl.DataFrame
    median_cohorts_per_nhood: float
    median_cohort_ilisi_nhood: float
    median_cohort_entropy_bits: float
    median_cancer_entropy_bits: float
    median_donor_lisi: float
    passes_mixing_qc_pct: float


class CachedAuditTables(BaseModel):
    """Pre-computed integration and null-calibration tables injected for fast heatmap refresh."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    selected_regime: HarmonyRegime
    harmony_benchmark_df: pl.DataFrame
    pan_nhood_mixing_df: pl.DataFrame
    power_ablation_df: pl.DataFrame
    per_cohort_uplift_df: pl.DataFrame


class PanCohortPipelineOutputs(BaseModel):
    """All tabular, statistical, and figure deliverables of the pan-cohort pipeline."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    selected_regime: HarmonyRegime
    harmony_benchmark_df: pl.DataFrame
    pan_nhood_mixing_df: pl.DataFrame
    power_ablation_df: pl.DataFrame
    per_cohort_uplift_df: pl.DataFrame
    two_stage_meta9_df: pl.DataFrame
    concordance_1s_vs_2s_df: pl.DataFrame
    pan_verified_pairs_df: pl.DataFrame
    loco_forest_df: pl.DataFrame
    markdown_report: str


class PanCohortFeatureLoader(Protocol):
    """Side-effecting loader producing an outcome-blind `PanCohortFeatureBundle`."""

    def __call__(self, cfg: PanCohortCliConfig) -> Result[PanCohortFeatureBundle, str]: ...


class VectorFigureExporter(Protocol):
    """Side-effecting exporter writing Inkscape-layered `.svg` and 300 DPI `.png` files."""

    def __call__(
        self,
        chart: AltairComposite,
        base_path: Path,
        figure_title: str,
        layer_label: str,
    ) -> Result[tuple[Path, Path], str]: ...


# =============================================================================
# 2. Pure Functional Core — Harmony Integration & Multi-Scale Mixing Metrics
# =============================================================================


def _build_harmony_metadata_frame(
    cohort_ids: Sequence[str],
    substudy_ids: Sequence[str],
    sample_uids: Sequence[str],
) -> pd.DataFrame:
    """Isolated third-party adapter building the `pd.DataFrame` required by `harmonypy`."""
    return pd.DataFrame(
        {
            "cohort_id": [str(x) for x in cohort_ids],
            "substudy_id": [str(x) for x in substudy_ids],
            "sample_uid": [str(x) for x in sample_uids],
        }
    )


def run_harmony_embedding_regime(
    x_pca: FloatArray,
    cohort_ids: Sequence[str],
    substudy_ids: Sequence[str],
    sample_uids: Sequence[str],
    regime: HarmonyRegime,
    nclust: int = 36,
    sigma: float = 0.22,
    max_iter_harmony: int = 15,
    random_state: int = 42,
) -> Result[FloatArray, str]:
    """Executes `harmonypy` latent batch integration for `regime` without mutating inputs."""
    mat = np.asarray(x_pca, dtype=np.float64)
    if mat.ndim != 2 or int(mat.shape[0]) < 10:
        return Failure(f"x_pca must be 2D with >= 10 cells, got shape {mat.shape}")
    n_cells, n_dims = int(mat.shape[0]), int(mat.shape[1])
    if not (len(cohort_ids) == n_cells and len(substudy_ids) == n_cells and len(sample_uids) == n_cells):
        return Failure("Metadata sequence lengths must exactly match x_pca cell count")
    if not bool(np.all(np.isfinite(mat))):
        return Failure("x_pca must contain strictly finite values")

    match regime:
        case HarmonyRegime.I0_UNINTEGRATED_PCA:
            return Success(np.array(mat, dtype=np.float64, copy=True))
        case HarmonyRegime.I1_HARMONY_COHORT:
            vars_use: list[str] = ["cohort_id"]
            theta_vals: list[float] = [2.0]
        case HarmonyRegime.I2_HARMONY_COHORT_SUBSTUDY:
            vars_use = ["cohort_id", "substudy_id"]
            theta_vals = [2.0, 1.5]
        case HarmonyRegime.I3_HARMONY_COHORT_SUBSTUDY_DONOR:
            vars_use = ["cohort_id", "substudy_id", "sample_uid"]
            theta_vals = [2.0, 1.5, 1.0]
        case _ as unreachable:
            assert_never(unreachable)

    eff_nclust = min(max(int(nclust), 4), max(4, n_cells // 12))
    meta_pd = _build_harmony_metadata_frame(cohort_ids, substudy_ids, sample_uids)
    try:
        ho = hm.run_harmony(
            np.array(mat, dtype=np.float64, copy=True),
            meta_pd,
            vars_use,
            theta=theta_vals,
            sigma=float(sigma),
            nclust=eff_nclust,
            max_iter_harmony=int(max_iter_harmony),
            random_state=int(random_state),
            verbose=False,
        )
        z_raw = np.asarray(ho.Z_corr, dtype=np.float64)
    except Exception as exc:
        return Failure(f"harmonypy failed on {regime.value}: {exc}")

    z_aligned = (
        z_raw
        if z_raw.shape == (n_cells, n_dims)
        else (z_raw.T if z_raw.shape == (n_dims, n_cells) else np.empty((0, 0), dtype=np.float64))
    )
    if z_aligned.shape != (n_cells, n_dims):
        return Failure(
            f"Unexpected harmonypy Z_corr shape {z_raw.shape}, expected {(n_cells, n_dims)}"
        )
    return Success(np.array(z_aligned, dtype=np.float64, copy=True))


def _compute_knn_inverse_simpson(
    knn_indices: IntArray,
    knn_dists: FloatArray,
    labels: Sequence[str],
) -> FloatArray:
    """Computes Gaussian-kernel local inverse Simpson index (LISI; Korsunsky et al. 2019)."""
    n_eval, k_nn = int(knn_indices.shape[0]), int(knn_indices.shape[1])
    uniq_labels = sorted(set(str(x) for x in labels))
    if len(uniq_labels) <= 1 or n_eval == 0 or k_nn == 0:
        return np.ones(n_eval, dtype=np.float64)

    lbl_map = {name: idx for idx, name in enumerate(uniq_labels)}
    codes = np.array([lbl_map[str(x)] for x in labels], dtype=np.int64)
    nbr_codes = codes[knn_indices]

    scale = np.maximum(np.median(knn_dists, axis=1, keepdims=True), 1e-6)
    weights = np.exp(-0.5 * np.square(knn_dists / scale))
    w_sum = np.maximum(np.sum(weights, axis=1, keepdims=True), 1e-12)
    probs_nbr = weights / w_sum

    simpson = np.zeros(n_eval, dtype=np.float64)
    n_cats = len(uniq_labels)
    if n_cats <= 64:
        for c_id in range(n_cats):
            p_c = np.sum(np.where(nbr_codes == c_id, probs_nbr, 0.0), axis=1)
            simpson = simpson + np.square(p_c)
    else:
        for row_i in range(n_eval):
            r_codes = nbr_codes[row_i]
            r_probs = probs_nbr[row_i]
            uniq_r = np.unique(r_codes)
            s_val = sum(float(np.sum(r_probs[r_codes == u_c])) ** 2 for u_c in uniq_r.tolist())
            simpson[row_i] = s_val

    return np.asarray(1.0 / np.maximum(simpson, 1e-12), dtype=np.float64)


def _compute_compartment_aware_kbet_acceptance(
    knn_indices: IntArray,
    cohort_ids: Sequence[str],
    compartments: Sequence[str],
    alpha: float = 0.05,
) -> float:
    """Computes kBET chi-square acceptance rate relative to observable compartment cohort shares."""
    n_eval, k_nn = int(knn_indices.shape[0]), int(knn_indices.shape[1])
    if n_eval == 0 or k_nn < 5:
        return 0.0

    uniq_cohorts = sorted(set(str(c) for c in cohort_ids))
    if len(uniq_cohorts) <= 1:
        return 1.0

    coh_map = {c: idx for idx, c in enumerate(uniq_cohorts)}
    coh_codes = np.array([coh_map[str(c)] for c in cohort_ids], dtype=np.int64)
    comp_arr = np.array([str(cp) for cp in compartments], dtype=str)
    n_coh = len(uniq_cohorts)

    comp_expected: dict[str, FloatArray] = {}
    for cp in sorted(set(comp_arr.tolist())):
        cp_mask = comp_arr == cp
        counts_cp = np.bincount(coh_codes[cp_mask], minlength=n_coh).astype(np.float64)
        tot_cp = float(np.sum(counts_cp))
        comp_expected[cp] = counts_cp / max(tot_cp, 1.0)

    accepted = 0
    tested = 0
    for i in range(n_eval):
        cp_i = str(comp_arr[int(knn_indices[i, 0])])
        exp_props = comp_expected.get(cp_i, np.full(n_coh, 1.0 / float(n_coh)))
        active_c = exp_props > 0.005
        df_c = int(np.sum(active_c)) - 1
        if df_c < 1:
            continue
        obs_counts = np.bincount(coh_codes[knn_indices[i]], minlength=n_coh).astype(np.float64)[active_c]
        exp_counts = np.maximum(exp_props[active_c] * float(k_nn), 0.25)
        chi2_stat = float(np.sum(np.square(obs_counts - exp_counts) / exp_counts))
        p_val = float(chi2.sf(chi2_stat, df=df_c))
        tested += 1
        if p_val >= alpha:
            accepted += 1

    return float(accepted) / float(max(tested, 1))


def evaluate_embedding_integration_metrics(
    embedding: FloatArray,
    cohort_ids: Sequence[str],
    substudy_ids: Sequence[str],
    sample_uids: Sequence[str],
    compartments: Sequence[str],
    meta_states: Sequence[str],
    regime: HarmonyRegime = HarmonyRegime.I0_UNINTEGRATED_PCA,
    knn_k: int = 30,
    max_eval_cells: int = 8000,
    seed: int = 42,
) -> Result[EmbeddingIntegrationMetrics, str]:
    """Quantifies batch mixing (`iLISI`, `kBET`, `batch_ASW`) and bio-conservation (`cLISI`, `ASW`)."""
    emb = np.asarray(embedding, dtype=np.float64)
    n_cells = int(emb.shape[0])
    if emb.ndim != 2 or n_cells < 20:
        return Failure(f"embedding must be 2D with >= 20 cells, got {emb.shape}")
    if not (
        len(cohort_ids) == n_cells
        and len(substudy_ids) == n_cells
        and len(sample_uids) == n_cells
        and len(compartments) == n_cells
        and len(meta_states) == n_cells
    ):
        return Failure("Label lengths must equal embedding row count")

    eff_k = min(int(knn_k), n_cells - 1)
    rng = np.random.default_rng(seed)
    eval_idx = (
        np.sort(rng.choice(n_cells, size=max_eval_cells, replace=False))
        if n_cells > max_eval_cells
        else np.arange(n_cells, dtype=np.int64)
    )

    tree = cKDTree(emb)
    dists_raw, idx_raw = tree.query(emb[eval_idx], k=eff_k + 1)
    knn_dists = np.asarray(dists_raw, dtype=np.float64)[:, 1:]
    knn_idx = np.asarray(idx_raw, dtype=np.int64)[:, 1:]

    coh_ilisi = float(np.median(_compute_knn_inverse_simpson(knn_idx, knn_dists, cohort_ids)))
    sub_ilisi = float(np.median(_compute_knn_inverse_simpson(knn_idx, knn_dists, substudy_ids)))
    don_ilisi = float(np.median(_compute_knn_inverse_simpson(knn_idx, knn_dists, sample_uids)))
    comp_clisi = float(np.median(_compute_knn_inverse_simpson(knn_idx, knn_dists, compartments)))
    meta_clisi = float(np.median(_compute_knn_inverse_simpson(knn_idx, knn_dists, meta_states)))

    sil_n = min(int(eval_idx.shape[0]), 4000)
    sil_idx = eval_idx[:sil_n]
    coh_sub = [str(cohort_ids[int(i)]) for i in sil_idx.tolist()]
    comp_sub = [str(compartments[int(i)]) for i in sil_idx.tolist()]
    batch_asw = (
        float(silhouette_score(emb[sil_idx], coh_sub, metric="euclidean"))
        if len(set(coh_sub)) >= 2
        else 0.0
    )
    comp_asw = (
        float(silhouette_score(emb[sil_idx], comp_sub, metric="euclidean"))
        if len(set(comp_sub)) >= 2
        else 0.0
    )

    # Align knn_idx row 0 with self for compartment lookup in kBET
    knn_with_self = np.column_stack([eval_idx, knn_idx[:, : max(eff_k - 1, 1)]])
    kbet_acc = _compute_compartment_aware_kbet_acceptance(
        knn_indices=knn_with_self,
        cohort_ids=cohort_ids,
        compartments=compartments,
    )

    return Success(
        EmbeddingIntegrationMetrics(
            regime=regime,
            cohort_ilisi=coh_ilisi,
            substudy_ilisi=sub_ilisi,
            donor_ilisi=don_ilisi,
            compartment_clisi=comp_clisi,
            meta_state_clisi=meta_clisi,
            batch_asw=batch_asw,
            compartment_asw=comp_asw,
            kbet_acceptance_rate=kbet_acc,
        )
    )


def _shannon_entropy_bits_from_counts(counts_2d: FloatArray) -> FloatArray:
    """Row-wise Shannon entropy in bits (`log2`) across group columns."""
    totals = np.maximum(np.sum(counts_2d, axis=1, keepdims=True), 1e-12)
    probs = counts_2d / totals
    safe_log2 = np.where(probs > 1e-12, np.log2(np.maximum(probs, 1e-12)), 0.0)
    ent = -np.sum(probs * safe_log2, axis=1)
    return np.asarray(np.maximum(ent, 0.0), dtype=np.float64)


def _inverse_simpson_from_counts(counts_2d: FloatArray) -> FloatArray:
    """Row-wise inverse Simpson diversity (`1 / sum p_g^2`) across group columns."""
    totals = np.maximum(np.sum(counts_2d, axis=1, keepdims=True), 1e-12)
    probs = counts_2d / totals
    simp = np.sum(np.square(probs), axis=1)
    has_mass = np.ravel(np.sum(counts_2d, axis=1)) > 0.0
    return np.asarray(np.where(has_mass, 1.0 / np.maximum(simp, 1e-12), 0.0), dtype=np.float64)


def build_pan_cohort_milo_atlas(
    embedding: FloatArray,
    bundle: PanCohortFeatureBundle,
    regime: HarmonyRegime,
    knn_k: int = 30,
    milo_prop: float = 0.06,
    seed: int = 42,
) -> Result[PanCohortMiloAtlas, str]:
    """Constructs refined Milo kNN neighborhoods on `embedding` and projects all `371,568` cells."""
    emb = np.asarray(embedding, dtype=np.float64)
    n_anchor = int(emb.shape[0])
    eff_k = min(int(knn_k), n_anchor - 1)
    if eff_k < 5:
        return Failure("Insufficient anchor cells to construct pan-cohort Milo kNN graph")

    tree = cKDTree(emb)
    _, knn_idx = tree.query(emb, k=eff_k + 1)
    knn_all = np.asarray(knn_idx, dtype=np.int64)

    rng = np.random.default_rng(seed)
    n_init = max(24, int(np.round(float(n_anchor) * milo_prop)))
    init_verts = rng.choice(n_anchor, size=min(n_init, n_anchor), replace=False)

    refined_list: list[int] = []
    for v_init in init_verts.tolist():
        nbrs = knn_all[int(v_init)]
        med_prof = np.median(emb[nbrs], axis=0)
        dists = np.linalg.norm(emb[nbrs] - med_prof, axis=1)
        refined_list.append(int(nbrs[int(np.argmin(dists))]))

    unique_verts = np.unique(np.array(refined_list, dtype=np.int64))
    n_nhoods = int(unique_verts.shape[0])

    row_blocks: list[IntArray] = []
    col_blocks: list[IntArray] = []
    nhood_meta_codes: list[str] = []
    nhood_compartments: list[str] = []
    nhood_readable: list[str] = []

    immune_lineage_set = {
        TmeCompartment.T_NK.value,
        TmeCompartment.B_LINEAGE.value,
        TmeCompartment.PLASMA.value,
        TmeCompartment.MYELOID_MAST.value,
    }
    anchor_meta_arr = np.array(bundle.anchor_meta_codes, dtype=str)
    anchor_comp_arr = np.array(bundle.anchor_compartments, dtype=str)
    anchor_is_immune = np.isin(anchor_comp_arr, list(immune_lineage_set))
    for nh_col, v_idx in enumerate(unique_verts.tolist()):
        raw_members = knn_all[int(v_idx)]
        same_lineage = anchor_is_immune[raw_members] == bool(anchor_is_immune[int(v_idx)])
        members = raw_members[same_lineage] if int(np.sum(same_lineage)) >= 5 else raw_members
        row_blocks.append(members)
        col_blocks.append(np.full(members.shape[0], nh_col, dtype=np.int64))
        vals, cnts = np.unique(anchor_meta_arr[members], return_counts=True)
        dom_meta = str(vals[int(np.argmax(cnts))])
        nhood_meta_codes.append(dom_meta)
        nhood_compartments.append(classify_tme_compartment(dom_meta).value)
        nhood_readable.append(bundle.code_to_readable.get(dom_meta, dom_meta))

    rows_cat = np.concatenate(row_blocks)
    cols_cat = np.concatenate(col_blocks)
    nhoods_anchor = sp.csc_matrix(
        (np.ones(rows_cat.shape[0], dtype=np.float64), (rows_cat, cols_cat)),
        shape=(n_anchor, n_nhoods),
    )

    # Project all 371,568 cells per cohort onto anchor cells via within-cohort 1-NN in outcome-blind PCA
    n_evals = len(bundle.eval_sample_uids)
    anchor_eval_counts = np.zeros((n_anchor, n_evals), dtype=np.float64)
    anchor_coh_arr = np.array(bundle.anchor_cohort_ids, dtype=str)

    for coh_slice in bundle.cohorts:
        c_mask = np.flatnonzero(anchor_coh_arr == coh_slice.accession)
        if int(c_mask.shape[0]) == 0:
            continue
        coh_anchor_pca = bundle.anchor_x_pca[c_mask]
        coh_tree = cKDTree(coh_anchor_pca)
        _, nn_local = coh_tree.query(coh_slice.full_cell_pca, k=1)
        mapped_anchor_rows = c_mask[np.asarray(nn_local, dtype=np.int64)]
        np.add.at(anchor_eval_counts, (mapped_anchor_rows, coh_slice.full_cell_eval_global_idx), 1.0)

    counts_full = np.asarray(nhoods_anchor.T.dot(anchor_eval_counts), dtype=np.float64)

    # Audit post-Milo pan-cohort mixing across cohorts, cancer types, and donors
    uniq_cohorts = sorted(set(bundle.eval_cohort_ids))
    coh_to_col = {c: idx for idx, c in enumerate(uniq_cohorts)}
    coh_ind = np.zeros((n_evals, len(uniq_cohorts)), dtype=np.float64)
    for j_idx, c_name in enumerate(bundle.eval_cohort_ids):
        coh_ind[j_idx, coh_to_col[c_name]] = 1.0
    nhood_coh_counts = counts_full @ coh_ind

    uniq_cancers = sorted(set(bundle.eval_cancer_types))
    can_to_col = {c: idx for idx, c in enumerate(uniq_cancers)}
    can_ind = np.zeros((n_evals, len(uniq_cancers)), dtype=np.float64)
    for j_idx, c_type in enumerate(bundle.eval_cancer_types):
        can_ind[j_idx, can_to_col[c_type]] = 1.0
    nhood_can_counts = counts_full @ can_ind

    n_dist_cohorts = np.sum(nhood_coh_counts >= 1.0, axis=1).astype(np.int64)
    coh_ilisi_nh = _inverse_simpson_from_counts(nhood_coh_counts)
    coh_entropy = _shannon_entropy_bits_from_counts(nhood_coh_counts)
    n_dist_cancers = np.sum(nhood_can_counts >= 1.0, axis=1).astype(np.int64)
    can_entropy = _shannon_entropy_bits_from_counts(nhood_can_counts)

    mix_res = compute_nhood_mixing_metrics(
        nhoods_mat=nhoods_anchor,
        count_matrix=counts_full,
        sample_ids=bundle.eval_sample_uids,
        batch_labels=Some(bundle.eval_substudy_strata),
        nhood_labels=Some(nhood_meta_codes),
        config=MixingQCConfig(
            enabled=True,
            min_distinct_donors=5,
            min_effective_donors_lisi=3.0,
            max_single_donor_dominance=0.40,
        ),
    )
    match mix_res:
        case Failure(err):
            return Failure(err)
        case Success(base_mix_df):
            pass

    nhood_mixing_df = base_mix_df.with_columns(
        pl.lit(regime.value).alias("harmony_regime"),
        pl.Series("meta_state", nhood_meta_codes, dtype=pl.String),
        pl.Series("tme_compartment", nhood_compartments, dtype=pl.String),
        pl.Series("n_distinct_cohorts", n_dist_cohorts, dtype=pl.Int64),
        pl.Series("cohort_ilisi_nhood", coh_ilisi_nh, dtype=pl.Float64),
        pl.Series("cohort_shannon_entropy_bits", coh_entropy, dtype=pl.Float64),
        pl.Series("n_distinct_cancers", n_dist_cancers, dtype=pl.Int64),
        pl.Series("cancer_shannon_entropy_bits", can_entropy, dtype=pl.Float64),
    )

    qc_pct = 100.0 * float(np.mean(nhood_mixing_df["passes_mixing_qc"].to_numpy()))
    return Success(
        PanCohortMiloAtlas(
            regime=regime,
            nhoods_anchor=nhoods_anchor,
            anchor_eval_counts=anchor_eval_counts,
            counts_full_pan=counts_full,
            nhood_meta_codes=tuple(nhood_meta_codes),
            nhood_compartments=tuple(nhood_compartments),
            nhood_readable_labels=tuple(nhood_readable),
            nhood_mixing_df=nhood_mixing_df,
            median_cohorts_per_nhood=float(np.median(n_dist_cohorts)),
            median_cohort_ilisi_nhood=float(np.median(coh_ilisi_nh)),
            median_cohort_entropy_bits=float(np.median(coh_entropy)),
            median_cancer_entropy_bits=float(np.median(can_entropy)),
            median_donor_lisi=float(np.median(nhood_mixing_df["effective_donors_lisi"].to_numpy())),
            passes_mixing_qc_pct=qc_pct,
        )
    )


def select_pareto_optimal_harmony_regime(
    metrics_by_regime: Mapping[HarmonyRegime, EmbeddingIntegrationMetrics],
    atlases_by_regime: Mapping[HarmonyRegime, PanCohortMiloAtlas],
) -> HarmonyRegime:
    """Selects the Pareto-optimal Harmony regime maximizing batch mixing without bio-blurring.

    Guards against donor-level over-correction (`I3`) when `meta_state_clisi` inflates or
    compartment silhouette degrades relative to `I2_harmony_cohort_substudy`.
    """
    i0_comp_asw = metrics_by_regime[HarmonyRegime.I0_UNINTEGRATED_PCA].compartment_asw
    best_reg = HarmonyRegime.I2_HARMONY_COHORT_SUBSTUDY
    best_score = -1e9

    for reg in ALL_HARMONY_REGIMES:
        if reg == HarmonyRegime.I0_UNINTEGRATED_PCA:
            continue
        m = metrics_by_regime[reg]
        a = atlases_by_regime[reg]
        bio_ok = (m.compartment_clisi <= 1.15) and (m.meta_state_clisi <= 1.25)
        # Penalize excessive fine-state blurring (meta_state_clisi) and donor-level outcome confounding
        donor_overcorr_penalty = 1.60 if reg == HarmonyRegime.I3_HARMONY_COHORT_SUBSTUDY_DONOR else 0.0
        composite = (
            1.00 * m.cohort_ilisi
            + 0.45 * m.substudy_ilisi
            + 1.80 * m.kbet_acceptance_rate
            + 0.80 * a.median_cohort_entropy_bits
            + 0.015 * a.passes_mixing_qc_pct
            - 4.50 * max(m.compartment_clisi - 1.08, 0.0)
            - 2.50 * max(m.meta_state_clisi - 1.20, 0.0)
            - 0.20 * max(i0_comp_asw - m.compartment_asw, 0.0)
            - donor_overcorr_penalty
        )
        eff_score = composite if bio_ok else (composite - 10.0)
        if eff_score > best_score:
            best_score = eff_score
            best_reg = reg

    return best_reg


# =============================================================================
# 3. Pure Functional Core — Seven-Tier Power Engine, Null Control & LOCO Audit
# =============================================================================


def build_pan_iatlas_state_pair_weights(
    iatlas_df: pl.DataFrame,
    base_floor: float = 0.75,
    slope: float = 1.15,
    z_offset: float = 1.30,
    max_cap: float = 3.8,
) -> dict[tuple[str, str], float]:
    """Maps independent Pan-iAtlas (`N = 1,015` bulk RNA-seq) `|epistasis_z|` to prior weights."""
    if iatlas_df.is_empty() or not {"state_a", "state_b", "epistasis_z"}.issubset(
        set(iatlas_df.columns)
    ):
        return {}
    weights: dict[tuple[str, str], float] = {}
    for sa, sb, ez in iatlas_df.select(["state_a", "state_b", "epistasis_z"]).iter_rows():
        sa_s = str(sa).strip()
        sb_s = str(sb).strip()
        if sa_s and sb_s and sa_s != sb_s:
            excess = max(abs(float(ez)) - z_offset, 0.0)
            w_val = float(np.clip(base_floor + slope * excess, base_floor, max_cap))
            key = (sa_s, sb_s) if sa_s < sb_s else (sb_s, sa_s)
            weights[key] = max(weights.get(key, base_floor), w_val)
    return weights


def expand_state_weights_to_nhood_pairs(
    nhood_codes: Sequence[str],
    state_weights: Mapping[tuple[str, str], float],
) -> dict[tuple[int, int], float]:
    """Projects canonical `(state_a, state_b)` prior weights onto `(nhood_u, nhood_v)` index pairs."""
    if not state_weights:
        return {}
    out: dict[tuple[int, int], float] = {}
    n_nh = len(nhood_codes)
    codes = [str(c).strip() for c in nhood_codes]
    for u in range(n_nh):
        cu = codes[u]
        for v in range(u + 1, n_nh):
            cv = codes[v]
            if cu != cv:
                key = (cu, cv) if cu < cv else (cv, cu)
                if key in state_weights:
                    out[(u, v)] = float(state_weights[key])
    return out


def permute_outcomes_within_strata(
    response_labels: Sequence[str],
    strata_labels: Sequence[str],
    seed: int,
) -> tuple[str, ...]:
    """Permutes `responder` / `non-responder` strictly within each `cohort::substudy` stratum."""
    n = len(response_labels)
    out = [str(r) for r in response_labels]
    strata_arr = np.array([str(s) for s in strata_labels], dtype=str)
    rng = np.random.default_rng(seed)
    for st_val in sorted(set(strata_arr.tolist())):
        idx = np.flatnonzero(strata_arr == st_val)
        if int(idx.shape[0]) > 1:
            perm_vals = [out[int(i)] for i in rng.permutation(idx).tolist()]
            for k_pos, orig_i in enumerate(idx.tolist()):
                out[int(orig_i)] = perm_vals[k_pos]
    return tuple(out)


class PreparedSuperUnit(BaseModel):
    """Coarse-grained super-neighborhood representation for one cohort or the pan-cohort atlas."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    name: str
    nhoods_super: SparseMatrix
    counts_super: FloatArray
    library_sizes: FloatArray
    cohort_labels: tuple[str, ...]
    strata_labels: tuple[str, ...]
    super_meta_codes: tuple[str, ...]
    super_readable_labels: tuple[str, ...]


def prepare_super_nhood_unit(
    name: str,
    nhoods_mat: SparseMatrix,
    count_matrix: FloatArray,
    library_sizes: FloatArray,
    cohort_labels: Sequence[str],
    strata_labels: Sequence[str],
    nhood_meta_codes: Sequence[str],
    code_to_readable: Mapping[str, str],
    cell_sample_codes: Maybe[IntArray] = Nothing,
    anchor_eval_counts: Maybe[FloatArray] = Nothing,
    resolution: float = 0.25,
    min_jaccard: float = 0.02,
    seed: int = 42,
) -> Result[PreparedSuperUnit, str]:
    """Pools overlapping kNN spheres within each `Meta42` state into consolidated super-neighborhoods."""
    cg_res = coarse_grain_nhoods_by_graph(
        nhoods_mat=nhoods_mat,
        count_matrix=count_matrix,
        cell_sample_codes=cell_sample_codes,
        sample_library_sizes=Some(library_sizes),
        nhood_group_labels=Some(tuple(str(c) for c in nhood_meta_codes)),
        min_jaccard=min_jaccard,
        leiden_resolution=resolution,
    )
    match cg_res:
        case Failure(err):
            return Failure(f"[{name}] Super-nhood pooling failed: {err}")
        case Success(cg):
            pass

    # Outcome-blind consolidation of disconnected micro-islands within each Meta42 state
    # (retains up to 4 major spatial Leiden sub-domains per state while absorbing <3-sphere islands)
    state_clusters: dict[str, list[list[int]]] = {}
    for lbl, mem_tup in zip(cg.super_nhood_labels, cg.member_nhood_indices):
        state_clusters.setdefault(str(lbl), []).append(list(mem_tup))

    consolidated_groups: list[tuple[str, tuple[int, ...]]] = []
    for st_lbl in sorted(state_clusters.keys()):
        clusters = sorted(state_clusters[st_lbl], key=lambda m: (-len(m), m[0]))
        tot_st_nhoods = sum(len(m) for m in clusters)
        min_major = max(3, tot_st_nhoods // 8)
        majors = [list(m) for m in clusters[:4] if len(m) >= min_major]
        if not majors:
            majors = [list(clusters[0])]
            start_minor = 1
        else:
            start_minor = len(majors)
        for minor_m in clusters[start_minor:]:
            majors[0].extend(minor_m)
        for maj in majors:
            consolidated_groups.append((st_lbl, tuple(sorted(set(maj)))))

    n_nhoods_orig = int(nhoods_mat.shape[1])
    n_cells_orig = int(nhoods_mat.shape[0])
    n_samples = int(count_matrix.shape[1])
    n_super = len(consolidated_groups)

    map_rows: list[int] = []
    map_cols: list[int] = []
    super_codes_list: list[str] = []
    for g_idx, (st_lbl, mem_sorted) in enumerate(consolidated_groups):
        super_codes_list.append(st_lbl)
        for v_idx in mem_sorted:
            map_rows.append(int(v_idx))
            map_cols.append(g_idx)

    bin_nhoods = (sp.csc_matrix(nhoods_mat) > 0).astype(np.float64)
    group_map_csc = sp.csc_matrix(
        (np.ones(len(map_rows), dtype=np.float64), (np.asarray(map_rows), np.asarray(map_cols))),
        shape=(n_nhoods_orig, n_super),
    )
    nhoods_super = (bin_nhoods.dot(group_map_csc) > 0).astype(np.float64).tocsc()

    match anchor_eval_counts:
        case Some(anc_w):
            counts_super = np.asarray(nhoods_super.T.dot(anc_w), dtype=np.float64)
        case Maybe.empty:
            match cell_sample_codes:
                case Some(codes):
                    sample_ind = sp.csr_matrix(
                        (
                            np.ones(n_cells_orig, dtype=np.float64),
                            (np.arange(n_cells_orig, dtype=np.int64), np.asarray(codes, dtype=np.int64)),
                        ),
                        shape=(n_cells_orig, n_samples),
                    )
                    counts_super = np.asarray(nhoods_super.T.dot(sample_ind).toarray(), dtype=np.float64)
                case Maybe.empty:
                    counts_super = np.minimum(
                        np.asarray(group_map_csc.T.dot(count_matrix), dtype=np.float64),
                        np.expand_dims(np.asarray(library_sizes, dtype=np.float64), axis=0),
                    )
                case _ as unreachable_cs:
                    return Failure(f"Unexpected Maybe state: {unreachable_cs}")
        case _ as unreachable_anc:
            return Failure(f"Unexpected Maybe state: {unreachable_anc}")

    super_codes = tuple(super_codes_list)
    super_read = tuple(code_to_readable.get(code, code) for code in super_codes)

    return Success(
        PreparedSuperUnit(
            name=name,
            nhoods_super=nhoods_super,
            counts_super=counts_super,
            library_sizes=np.asarray(library_sizes, dtype=np.float64),
            cohort_labels=tuple(str(c) for c in cohort_labels),
            strata_labels=tuple(str(s) for s in strata_labels),
            super_meta_codes=super_codes,
            super_readable_labels=super_read,
        )
    )


def _extract_sig_state_pairs_from_pw(pw: PairwiseAnalysisResult) -> set[tuple[str, str]]:
    """Extracts canonical `(state_a, state_b)` cross-compartment discoveries from `PairwiseAnalysisResult`."""
    pairs: set[tuple[str, str]] = set()
    for df in (pw.pairs_df, pw.continuous_pairs_df):
        if df.is_empty() or "is_significant" not in df.columns:
            continue
        sig = df.filter(pl.col("is_significant"))
        for u_s, v_s in sig.select(["cell_type_u", "cell_type_v"]).iter_rows():
            su, sv = str(u_s).strip(), str(v_s).strip()
            cu, cv = classify_tme_compartment(su), classify_tme_compartment(sv)
            if (
                su
                and sv
                and su != sv
                and cu != cv
                and cu != TmeCompartment.OTHER
                and cv != TmeCompartment.OTHER
            ):
                pairs.add((su, sv) if su < sv else (sv, su))
    return pairs


def _extract_sig_state_pairs_from_acat(acat_df: pl.DataFrame) -> set[tuple[str, str]]:
    """Extracts canonical `(state_a, state_b)` cross-compartment discoveries from `combine_pair_regimes_acat`."""
    if acat_df.is_empty() or "is_significant" not in acat_df.columns:
        return set()
    pairs: set[tuple[str, str]] = set()
    for u_s, v_s in acat_df.filter(pl.col("is_significant")).select(["cell_type_u", "cell_type_v"]).iter_rows():
        su, sv = str(u_s).strip(), str(v_s).strip()
        cu, cv = classify_tme_compartment(su), classify_tme_compartment(sv)
        if (
            su
            and sv
            and su != sv
            and cu != cv
            and cu != TmeCompartment.OTHER
            and cv != TmeCompartment.OTHER
        ):
            pairs.add((su, sv) if su < sv else (sv, su))
    return pairs


def _run_unit_pairwise(
    unit: PreparedSuperUnit,
    response_labels: Sequence[str],
    presence_mode: PresenceMode,
    residualize_depth: bool,
    mask_compartments: bool,
    dmpp_draws: int,
    fdr_threshold: float,
    prior_weights: Maybe[Mapping[tuple[int, int], float]],
    is_pan_cohort: bool,
) -> Result[PairwiseAnalysisResult, str]:
    """Executes `analyze_nhood_pairs` on a prepared super-neighborhood unit."""
    n_r = sum(1 for r in response_labels if r == "responder")
    n_nr = len(response_labels) - n_r
    min_arm = min(n_r, n_nr)
    min_q_donors = (
        5.0
        if is_pan_cohort
        else float(min(3.0, max(2.0, float(min_arm) * 0.30)))
    )
    min_eff_samp = 14.0 if is_pan_cohort else float(max(6, min(10, len(response_labels) - 1)))

    uniq_strata = set(unit.strata_labels)
    strata_arg: Maybe[Sequence[str]] = (
        Some(unit.strata_labels) if len(uniq_strata) >= 2 else Nothing
    )
    cohort_arg: Maybe[Sequence[str]] = (
        Some(unit.cohort_labels) if mask_compartments else Nothing
    )
    mix_cfg = MixingQCConfig(
        enabled=len(response_labels) >= 12,
        min_distinct_donors=4 if is_pan_cohort else 3,
        min_effective_donors_lisi=2.8 if is_pan_cohort else 2.2,
        max_single_donor_dominance=0.45 if is_pan_cohort else 0.55,
    )

    pw_cfg = PairwiseConfig(
        method=PairwiseMethod.HYBRID_ALL,
        overlap=OverlapConfig(
            max_jaccard_overlap=0.05,
            scope_filter=PairScopeFilter.CROSS_COMPARTMENT_ONLY,
        ),
        presence=PresenceConfig(
            mode=presence_mode,
            residualize_library_depth=residualize_depth,
            mask_unobserved_compartments=mask_compartments,
            min_cohort_compartment_share=0.015,
        ),
        prefilter=PrefilterConfig(
            mode=PrefilterMode.MARGINAL_ENTROPY_VARIANCE,
            min_marginal_prevalence=0.04 if is_pan_cohort else 0.06,
            max_marginal_prevalence=1.0,
            min_variance_quantile=0.02,
            apply_tarone_iut_pruning=True,
            tarone_p_ceiling=0.02,
        ),
        mixing_qc=mix_cfg,
        screen=ScreenConfig(
            fdr_threshold=fdr_threshold,
            min_donors_per_quadrant=min_q_donors,
            min_effective_samples=min_eff_samp,
            max_retained_candidates=50000,
            apply_tarone_iut_pruning=True,
            tarone_p_ceiling=0.02,
            dmpp_null_draws=max(dmpp_draws, 96) if (is_pan_cohort and dmpp_draws > 16) else dmpp_draws,
        ),
        refinement=RefinementConfig(enabled=False),
        modules=PairModuleConfig(enabled=True),
    )

    n_v = int(unit.counts_super.shape[0])
    prior_mat_maybe: Maybe[FloatArray] = Nothing
    match prior_weights:
        case Some(pw_dict) if pw_dict:
            p_mat = np.ones((n_v, n_v), dtype=np.float64)
            for (u_i, v_i), w_val in pw_dict.items():
                if 0 <= u_i < n_v and 0 <= v_i < n_v and u_i != v_i:
                    p_mat[u_i, v_i] = float(w_val)
                    p_mat[v_i, u_i] = float(w_val)
            prior_mat_maybe = Some(p_mat)
        case _:
            pass

    return analyze_nhood_pairs(
        nhoods_mat=unit.nhoods_super,
        count_matrix=unit.counts_super,
        response_labels=response_labels,
        sample_library_sizes=Some(tuple(unit.library_sizes.tolist())),
        strata_labels=strata_arg,
        sample_ids=Some(tuple(f"S{j:04d}" for j in range(len(response_labels)))),
        cohort_labels=cohort_arg,
        nhood_group_labels=Some(unit.super_meta_codes),
        prior_pair_weights=prior_mat_maybe,
        config=pw_cfg,
    )


def _build_hybrid_prior_acat_df(
    pw_prior: PairwiseAnalysisResult,
    acat_df: pl.DataFrame,
    prior_weights: Maybe[Mapping[tuple[int, int], float]],
    fdr_threshold: float,
) -> pl.DataFrame:
    """Combines Roeder-Wasserman prior-weighted 3-leg DMPP FDR with multi-regime Cauchy ACAT FDR."""
    if acat_df.is_empty():
        return acat_df

    best_dmpp: dict[tuple[int, int], float] = {}
    if not pw_prior.pairs_df.is_empty():
        for u_q, v_q, fq0, fqi in pw_prior.pairs_df.select(
            ["nhood_u", "nhood_v", "fdr_quadrant", "fdr_iut"]
        ).iter_rows():
            kuv = (min(int(u_q), int(v_q)), max(int(u_q), int(v_q)))
            q3 = max(float(fq0), float(fqi))
            best_dmpp[kuv] = min(best_dmpp.get(kuv, 1.0), q3)

    if not pw_prior.continuous_pairs_df.is_empty():
        for u_c, v_c, fc0, fci in pw_prior.continuous_pairs_df.select(
            ["nhood_u", "nhood_v", "fdr_axis", "fdr_iut"]
        ).iter_rows():
            kuv = (min(int(u_c), int(v_c)), max(int(u_c), int(v_c)))
            c3 = max(float(fc0), float(fci))
            best_dmpp[kuv] = min(best_dmpp.get(kuv, 1.0), c3)

    u_arr = acat_df["nhood_u"].to_list()
    v_arr = acat_df["nhood_v"].to_list()
    fdr_ac = acat_df["fdr_acat_iut"].to_numpy()
    p_ac = acat_df["pvalue_acat_iut"].to_numpy()
    p_min_reg = acat_df["min_regime_pvalue_iut"].to_numpy()

    dmpp_arr = np.array(
        [best_dmpp.get((min(int(uu), int(vv)), max(int(uu), int(vv))), 1.0) for uu, vv in zip(u_arr, v_arr)],
        dtype=np.float64,
    )
    fdr_comb = np.minimum(fdr_ac, dmpp_arr)
    p_eff = np.where(dmpp_arr < fdr_ac, np.minimum(p_ac, p_min_reg), p_ac)

    match prior_weights:
        case Some(pmap):
            pw_col = np.array(
                [
                    float(pmap.get((min(int(uu), int(vv)), max(int(uu), int(vv))), 1.0))
                    for uu, vv in zip(u_arr, v_arr)
                ],
                dtype=np.float64,
            )
        case _:
            pw_col = acat_df["prior_weight"].to_numpy()

    return (
        acat_df.with_columns(
            pl.Series("prior_weight", pw_col, dtype=pl.Float64),
            pl.Series("pvalue_acat_iut", p_eff, dtype=pl.Float64),
            pl.Series("fdr_acat_iut", fdr_comb, dtype=pl.Float64),
            pl.Series("is_significant", fdr_comb < fdr_threshold, dtype=pl.Boolean),
        )
        .sort(
            ["is_significant", "fdr_acat_iut", "pvalue_acat_iut"],
            descending=[True, False, False],
        )
    )


class SinglePassPowerBundle(BaseModel):
    """Detailed single-pass outputs across `P0` -> `P6` for real or permuted clinical outcomes."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    tier_state_pair_counts: Mapping[PowerTier, int]
    tier_nhood_pair_counts: Mapping[PowerTier, int]
    per_cohort_p0_counts: Mapping[str, int]
    per_cohort_p1_counts: Mapping[str, int]
    per_cohort_p2_counts: Mapping[str, int]
    meta9_df: pl.DataFrame
    pan_p5_pw: PairwiseAnalysisResult
    pan_p6_pw: PairwiseAnalysisResult
    pan_p6_acat_df: pl.DataFrame


def evaluate_seven_power_tiers_single_pass(
    cohort_units: Sequence[PreparedSuperUnit],
    cohort_eval_indices: Mapping[str, IntArray],
    pan_unit: PreparedSuperUnit,
    global_response_labels: Sequence[str],
    state_prior_weights: Mapping[tuple[str, str], float],
    cfg: PanCohortCliConfig,
) -> Result[SinglePassPowerBundle, str]:
    """Evaluates all 7 power tiers (`P0` -> `P6`) on `global_response_labels`."""
    p0_union: set[tuple[str, str]] = set()
    p1_union: set[tuple[str, str]] = set()
    p2_union: set[tuple[str, str]] = set()
    p0_nhood_tot = 0
    p1_nhood_tot = 0
    p2_nhood_tot = 0

    coh_p0: dict[str, int] = {}
    coh_p1: dict[str, int] = {}
    coh_p2: dict[str, int] = {}
    p1_cohort_tables: dict[str, PairwiseAnalysisResult] = {}
    cohort_sizes: dict[str, int] = {}

    for c_unit in cohort_units:
        c_idx = cohort_eval_indices[c_unit.name]
        y_coh = tuple(str(global_response_labels[int(i)]) for i in c_idx.tolist())
        cohort_sizes[c_unit.name] = len(y_coh)
        c_prior = expand_state_weights_to_nhood_pairs(c_unit.super_meta_codes, state_prior_weights)
        c_prior_maybe: Maybe[Mapping[tuple[int, int], float]] = (
            Some(c_prior) if c_prior else Nothing
        )

        # P0: Current per-cohort baseline (HARD_DEPTH_GATED, no depth resid, default B=12 cap)
        res_p0 = _run_unit_pairwise(
            unit=c_unit,
            response_labels=y_coh,
            presence_mode=PresenceMode.HARD_DEPTH_GATED,
            residualize_depth=False,
            mask_compartments=False,
            dmpp_draws=12,
            fdr_threshold=cfg.fdr_threshold,
            prior_weights=Nothing,
            is_pan_cohort=False,
        )
        match res_p0:
            case Failure(err):
                return Failure(err)
            case Success(pw_p0):
                st_p0 = _extract_sig_state_pairs_from_pw(pw_p0)
                coh_p0[c_unit.name] = len(st_p0)
                p0_union |= st_p0
                p0_nhood_tot += pw_p0.n_significant_quadrant_pairs + pw_p0.n_significant_continuous_pairs

        # P1: + Empirical-Bayes Hurdle Soft Presence + Depth Residualization + B=64 DMPP
        res_p1 = _run_unit_pairwise(
            unit=c_unit,
            response_labels=y_coh,
            presence_mode=PresenceMode.EMPIRICAL_BAYES_POSTERIOR,
            residualize_depth=True,
            mask_compartments=False,
            dmpp_draws=cfg.dmpp_null_draws,
            fdr_threshold=cfg.fdr_threshold,
            prior_weights=Nothing,
            is_pan_cohort=False,
        )
        match res_p1:
            case Failure(err):
                return Failure(err)
            case Success(pw_p1):
                st_p1 = _extract_sig_state_pairs_from_pw(pw_p1)
                coh_p1[c_unit.name] = len(st_p1)
                p1_union |= st_p1
                p1_nhood_tot += pw_p1.n_significant_quadrant_pairs + pw_p1.n_significant_continuous_pairs
                p1_cohort_tables[c_unit.name] = pw_p1

        # P2: + Omnibus Cauchy ACAT Regime Combination + Independent Pan-iAtlas Prior
        res_p2_pw = _run_unit_pairwise(
            unit=c_unit,
            response_labels=y_coh,
            presence_mode=PresenceMode.EMPIRICAL_BAYES_POSTERIOR,
            residualize_depth=True,
            mask_compartments=False,
            dmpp_draws=cfg.dmpp_null_draws,
            fdr_threshold=cfg.fdr_threshold,
            prior_weights=c_prior_maybe,
            is_pan_cohort=False,
        )
        match res_p2_pw:
            case Failure(err):
                return Failure(err)
            case Success(pw_p2):
                pass

        acat_c_res = combine_pair_regimes_acat(
            quadrant_df=pw_p2.pairs_df,
            continuous_df=pw_p2.continuous_pairs_df,
            fdr_threshold=cfg.fdr_threshold,
            prior_weights=Nothing,
            n_evaluated_pairs=0,
        )
        match acat_c_res:
            case Failure(err):
                return Failure(err)
            case Success(raw_acat_c_df):
                acat_c_df = _build_hybrid_prior_acat_df(
                    pw_prior=pw_p2,
                    acat_df=raw_acat_c_df,
                    prior_weights=c_prior_maybe,
                    fdr_threshold=cfg.fdr_threshold,
                )
                st_p2 = _extract_sig_state_pairs_from_acat(acat_c_df)
                coh_p2[c_unit.name] = len(st_p2)
                p2_union |= st_p2
                p2_nhood_tot += int(acat_c_df.filter(pl.col("is_significant")).height)

    # P3: Two-Stage 9-Cohort DerSimonian-Laird Random-Effects 3-Leg Meta-IUT
    meta_res = meta_analyze_cohort_pairs_iut(
        cohort_pair_tables=p1_cohort_tables,
        cohort_sample_sizes=cohort_sizes,
        fdr_threshold=cfg.fdr_threshold,
        require_cross_compartment=True,
        min_cohorts=2,
    )
    match meta_res:
        case Failure(err):
            return Failure(err)
        case Success(meta9_df):
            pass

    p3_pairs: set[tuple[str, str]] = set()
    p3_nhood_tot = 0
    if not meta9_df.is_empty():
        sig_m = meta9_df.filter(pl.col("is_significant"))
        p3_nhood_tot = int(sig_m.height)
        for sa, sb in sig_m.select(["state_a", "state_b"]).iter_rows():
            sa_s, sb_s = str(sa), str(sb)
            p3_pairs.add((sa_s, sb_s) if sa_s < sb_s else (sb_s, sa_s))

    # P4: 1-Stage Harmony Pooled (Unmasked Observability, Hard Presence, B=12)
    res_p4 = _run_unit_pairwise(
        unit=pan_unit,
        response_labels=global_response_labels,
        presence_mode=PresenceMode.HARD_DEPTH_GATED,
        residualize_depth=False,
        mask_compartments=False,
        dmpp_draws=12,
        fdr_threshold=cfg.fdr_threshold,
        prior_weights=Nothing,
        is_pan_cohort=True,
    )
    match res_p4:
        case Failure(err):
            return Failure(err)
        case Success(pw_p4):
            p4_pairs = _extract_sig_state_pairs_from_pw(pw_p4)
            p4_nhood_tot = pw_p4.n_significant_quadrant_pairs + pw_p4.n_significant_continuous_pairs

    # P5: 1-Stage Harmony Pooled + Assay-Aware CD45+ Compartment Observability Masking
    res_p5 = _run_unit_pairwise(
        unit=pan_unit,
        response_labels=global_response_labels,
        presence_mode=PresenceMode.HARD_DEPTH_GATED,
        residualize_depth=False,
        mask_compartments=True,
        dmpp_draws=cfg.dmpp_null_draws,
        fdr_threshold=cfg.fdr_threshold,
        prior_weights=Nothing,
        is_pan_cohort=True,
    )
    match res_p5:
        case Failure(err):
            return Failure(err)
        case Success(pw_p5):
            p5_pairs = _extract_sig_state_pairs_from_pw(pw_p5)
            p5_nhood_tot = pw_p5.n_significant_quadrant_pairs + pw_p5.n_significant_continuous_pairs

    # P6: Full Unified Pan-Cohort Power Engine (Harmony + Obs Mask + EB Soft + Depth + B=96 + ACAT + Prior)
    pan_prior = expand_state_weights_to_nhood_pairs(pan_unit.super_meta_codes, state_prior_weights)
    pan_prior_maybe: Maybe[Mapping[tuple[int, int], float]] = (
        Some(pan_prior) if pan_prior else Nothing
    )
    res_p6_base = _run_unit_pairwise(
        unit=pan_unit,
        response_labels=global_response_labels,
        presence_mode=PresenceMode.EMPIRICAL_BAYES_POSTERIOR,
        residualize_depth=True,
        mask_compartments=True,
        dmpp_draws=cfg.dmpp_null_draws,
        fdr_threshold=cfg.fdr_threshold,
        prior_weights=pan_prior_maybe,
        is_pan_cohort=True,
    )
    match res_p6_base:
        case Failure(err):
            return Failure(err)
        case Success(pw_p6):
            pass

    acat_pan_res = combine_pair_regimes_acat(
        quadrant_df=pw_p6.pairs_df,
        continuous_df=pw_p6.continuous_pairs_df,
        fdr_threshold=cfg.fdr_threshold,
        prior_weights=Nothing,
        n_evaluated_pairs=0,
    )
    match acat_pan_res:
        case Failure(err):
            return Failure(err)
        case Success(raw_pan_acat_df):
            pan_acat_df = _build_hybrid_prior_acat_df(
                pw_prior=pw_p6,
                acat_df=raw_pan_acat_df,
                prior_weights=pan_prior_maybe,
                fdr_threshold=cfg.fdr_threshold,
            )
            p6_pairs = _extract_sig_state_pairs_from_acat(pan_acat_df)
            p6_nhood_tot = int(pan_acat_df.filter(pl.col("is_significant")).height)

    state_counts: dict[PowerTier, int] = {
        PowerTier.P0_WITHIN_COHORT_BASELINE: len(p0_union),
        PowerTier.P1_WITHIN_COHORT_EB_DEPTH_B64: len(p1_union),
        PowerTier.P2_WITHIN_COHORT_ACAT_PRIOR: len(p2_union),
        PowerTier.P3_TWO_STAGE_META9_IUT: len(p3_pairs),
        PowerTier.P4_PAN_COHORT_HARMONY_NAIVE: len(p4_pairs),
        PowerTier.P5_PAN_COHORT_HARMONY_OBS_MASKED: len(p5_pairs),
        PowerTier.P6_PAN_COHORT_HARMONY_FULL_POWER: len(p6_pairs),
    }
    nhood_counts: dict[PowerTier, int] = {
        PowerTier.P0_WITHIN_COHORT_BASELINE: p0_nhood_tot,
        PowerTier.P1_WITHIN_COHORT_EB_DEPTH_B64: p1_nhood_tot,
        PowerTier.P2_WITHIN_COHORT_ACAT_PRIOR: p2_nhood_tot,
        PowerTier.P3_TWO_STAGE_META9_IUT: p3_nhood_tot,
        PowerTier.P4_PAN_COHORT_HARMONY_NAIVE: p4_nhood_tot,
        PowerTier.P5_PAN_COHORT_HARMONY_OBS_MASKED: p5_nhood_tot,
        PowerTier.P6_PAN_COHORT_HARMONY_FULL_POWER: p6_nhood_tot,
    }

    return Success(
        SinglePassPowerBundle(
            tier_state_pair_counts=state_counts,
            tier_nhood_pair_counts=nhood_counts,
            per_cohort_p0_counts=coh_p0,
            per_cohort_p1_counts=coh_p1,
            per_cohort_p2_counts=coh_p2,
            meta9_df=meta9_df,
            pan_p5_pw=pw_p5,
            pan_p6_pw=pw_p6,
            pan_p6_acat_df=pan_acat_df,
        )
    )


def _compute_exact_stratified_perm_pvalues_for_acat_pairs(
    pan_unit: PreparedSuperUnit,
    dedup_acat: pl.DataFrame,
    quad_pairs_df: pl.DataFrame,
    global_response_labels: Sequence[str],
    n_permutations: int,
    fdr_threshold: float,
    seed: int,
) -> dict[tuple[int, int], float]:
    """Computes exact within-stratum (`cohort::substudy`) 3-leg IUT permutation p-values for all pairs."""
    if dedup_acat.is_empty():
        return {}

    pres_base = estimate_nhood_presence(
        count_matrix=pan_unit.counts_super,
        library_sizes=Some(pan_unit.library_sizes),
        config=PresenceConfig(
            mode=PresenceMode.EMPIRICAL_BAYES_POSTERIOR,
            residualize_library_depth=True,
        ),
    )
    match pres_base:
        case Failure(_):
            return {}
        case Success(pres_raw):
            pass

    if len(set(pan_unit.cohort_labels)) > 1 and int(pres_raw.active_nhood_indices.shape[0]) > 0:
        act_c = pan_unit.counts_super[pres_raw.active_nhood_indices, :]
        eff_libs = pres_raw.effective_library_sizes
        raw_p = act_c / np.maximum(np.expand_dims(eff_libs, axis=0), 1.0)
        raw_arc = np.arcsin(np.sqrt(np.clip(raw_p, 0.0, 1.0)))
        raw_cpm = np.log2(1e6 * raw_p + 1.0)
        g_cpm = depth_residualize_continuous_features(
            raw_cpm, eff_libs, preserve_row_means=True, group_labels=Some(pan_unit.cohort_labels)
        ).value_or(pres_raw.log2_cpm)
        g_arc = depth_residualize_continuous_features(
            raw_arc, eff_libs, preserve_row_means=True, group_labels=Some(pan_unit.cohort_labels)
        ).value_or(pres_raw.arcsin_sqrt_prop)
        pres_work = PresenceRepresentation(
            active_nhood_indices=pres_raw.active_nhood_indices,
            presence_prob=pres_raw.presence_prob,
            absence_prob=pres_raw.absence_prob,
            observability=pres_raw.observability,
            raw_binary_presence=pres_raw.raw_binary_presence,
            arcsin_sqrt_prop=np.where(act_c > 0.0, np.clip(g_arc, 1e-8, 0.5 * np.pi), 0.0),
            log2_cpm=g_cpm,
            effective_library_sizes=eff_libs,
        )
    else:
        pres_work = pres_raw

    masked_pres = build_compartment_observability_mask(
        presence=pres_work,
        nhood_labels=pan_unit.super_meta_codes,
        cohort_labels=pan_unit.cohort_labels,
        min_cohort_compartment_share=0.015,
    ).value_or(pres_work)

    global_to_active = {
        int(g_idx): a_idx for a_idx, g_idx in enumerate(masked_pres.active_nhood_indices.tolist())
    }
    uniq_st = sorted(set(pan_unit.strata_labels))
    st_to_id = {s: idx for idx, s in enumerate(uniq_st)}
    st_codes = np.array([st_to_id[s] for s in pan_unit.strata_labels], dtype=np.int64)
    y_bin = np.array(
        [1.0 if r == "responder" else 0.0 for r in global_response_labels],
        dtype=np.float64,
    )

    perm_map: dict[tuple[int, int], float] = {}

    # 1. Discrete quadrant pairs via exact stratified CMH 3-leg IUT permutations
    quad_targets = dedup_acat.filter(pl.col("dominant_regime").str.starts_with("quadrant:"))
    if not quad_targets.is_empty() and not quad_pairs_df.is_empty():
        wanted_quad_map: dict[tuple[int, int], str] = {}
        for u_r, v_r, cu_r, cv_r, reg_r in quad_targets.select(
            ["nhood_u", "nhood_v", "cell_type_u", "cell_type_v", "dominant_regime"]
        ).iter_rows():
            q_dom = str(reg_r).removeprefix("quadrant:")
            if str(cu_r) > str(cv_r):
                q_dom = "-+" if q_dom == "+-" else ("+-" if q_dom == "-+" else q_dom)
            wanted_quad_map[(int(u_r), int(v_r))] = q_dom

        q_mask: list[bool] = []
        for u_q, v_q, q_raw in quad_pairs_df.select(["nhood_u", "nhood_v", "quadrant"]).iter_rows():
            uu, vv = int(u_q), int(v_q)
            q_can = str(q_raw)
            if uu > vv:
                uu, vv = vv, uu
                q_can = "-+" if q_can == "+-" else ("+-" if q_can == "-+" else q_can)
            q_mask.append(wanted_quad_map.get((uu, vv), "") == q_can)

        sub_q_df = quad_pairs_df.filter(pl.Series(q_mask, dtype=pl.Boolean))
        if not sub_q_df.is_empty():
            q_perm_res = calibrate_candidates_stratified_iut_permutation(
                candidates_df=sub_q_df,
                presence=masked_pres,
                response_binary=y_bin,
                strata_codes=Some(st_codes),
                n_permutations=n_permutations,
                fdr_threshold=fdr_threshold,
                seed=seed,
            )
            match q_perm_res:
                case Success(perm_df) if not perm_df.is_empty():
                    for pu, pv, pp in perm_df.select(
                        ["nhood_u", "nhood_v", "pvalue_perm_iut"]
                    ).iter_rows():
                        key_uv = (min(int(pu), int(pv)), max(int(pu), int(pv)))
                        perm_map[key_uv] = min(perm_map.get(key_uv, 1.0), float(pp))
                case _:
                    pass

    # 2. Continuous `log_ratio_balance` & `co_abundance_product` exact stratified 3-leg IUT permutations
    all_inf_groups: list[IntArray] = []
    for s_id in range(len(uniq_st)):
        grp = np.flatnonzero(st_codes == s_id).astype(np.int64)
        n_g = int(grp.shape[0])
        n_r_g = int(np.sum(y_bin[grp] > 0.5))
        if n_g >= 2 and 0 < n_r_g < n_g:
            all_inf_groups.append(grp)
    if not all_inf_groups:
        return perm_map

    rng = np.random.default_rng(seed + 7919)
    y_all_draws = np.repeat(np.expand_dims(y_bin, axis=1), n_permutations + 1, axis=1)
    for b_idx in range(1, n_permutations + 1):
        for grp in all_inf_groups:
            y_all_draws[grp, b_idx] = rng.permutation(y_bin[grp])

    cont_targets = dedup_acat.filter(~pl.col("dominant_regime").str.starts_with("quadrant:"))
    domain_buckets: dict[tuple[bool, ...], list[tuple[int, int, int, int, str]]] = {}
    for u_raw, v_raw, reg_raw in cont_targets.select(
        ["nhood_u", "nhood_v", "dominant_regime"]
    ).iter_rows():
        u_id, v_id = int(u_raw), int(v_raw)
        if u_id not in global_to_active or v_id not in global_to_active:
            continue
        au, av = global_to_active[u_id], global_to_active[v_id]
        sig_dom = tuple(
            bool(
                np.all(masked_pres.observability[au, grp] > 0.5)
                and np.all(masked_pres.observability[av, grp] > 0.5)
            )
            for grp in all_inf_groups
        )
        if any(sig_dom):
            domain_buckets.setdefault(sig_dom, []).append((u_id, v_id, au, av, str(reg_raw)))

    v_act = int(masked_pres.active_nhood_indices.shape[0])
    for sig_dom, pair_entries in domain_buckets.items():
        dom_groups = [grp for grp, ok in zip(all_inf_groups, sig_dom) if ok]
        used_cols = np.concatenate(dom_groups)
        n_used = int(used_cols.shape[0])
        if n_used < 8:
            continue
        df_resid = max(n_used - len(dom_groups), 3)
        df_one = max(df_resid - 1, 1)
        df_two = max(df_resid - 2, 1)
        df_inter = max(df_resid - 3, 1)

        cpm_cent = np.zeros((v_act, n_used), dtype=np.float64)
        arc_cent = np.zeros((v_act, n_used), dtype=np.float64)
        y_cent = np.zeros((n_used, n_permutations + 1), dtype=np.float64)

        offset = 0
        for grp in dom_groups:
            g_len = int(grp.shape[0])
            sl = slice(offset, offset + g_len)
            cg = masked_pres.log2_cpm[:, grp]
            ag = masked_pres.arcsin_sqrt_prop[:, grp]
            yg = y_all_draws[grp, :]
            cpm_cent[:, sl] = cg - np.mean(cg, axis=1, keepdims=True)
            arc_cent[:, sl] = ag - np.mean(ag, axis=1, keepdims=True)
            y_cent[sl, :] = yg - np.mean(yg, axis=0, keepdims=True)
            offset += g_len

        cpm_norms = np.linalg.norm(cpm_cent, axis=1, keepdims=True)
        z_cpm = np.where(cpm_norms > 1e-12, cpm_cent / np.maximum(cpm_norms, 1e-12), 0.0)
        cpm_sq_norms = np.sum(cpm_cent * cpm_cent, axis=1)

        arc_norms = np.linalg.norm(arc_cent, axis=1, keepdims=True)
        z_arc = np.where(arc_norms > 1e-12, arc_cent / np.maximum(arc_norms, 1e-12), 0.0)

        y_norms = np.linalg.norm(y_cent, axis=0, keepdims=True)
        y_unit = np.where(y_norms > 1e-12, y_cent / np.maximum(y_norms, 1e-12), 0.0)

        r_cpm_draws = z_cpm @ y_unit
        cpm_dot_y_draws = cpm_cent @ y_unit
        r_arc_draws = z_arc @ y_unit
        sqrt_used = float(np.sqrt(float(n_used)))

        for u_id, v_id, au, av, reg_str in pair_entries:
            if reg_str == ContinuousAxis.LOG_RATIO_BALANCE.value:
                dot_uv = float(np.dot(cpm_cent[au], cpm_cent[av]))
                diff_sq = max(float(cpm_sq_norms[au] + cpm_sq_norms[av] - 2.0 * dot_uv), 1e-15)
                r_bal = np.clip(
                    (cpm_dot_y_draws[au] - cpm_dot_y_draws[av]) / np.sqrt(diff_sq),
                    -0.999999,
                    0.999999,
                )
                t_bal = r_bal * np.sqrt(float(df_one) / np.maximum(1.0 - r_bal * r_bal, 1e-12))
                p_bal = np.asarray(2.0 * student_t.sf(np.abs(t_bal), df=df_one), dtype=np.float64)

                rho_uv = float(np.clip(float(np.dot(z_cpm[au], z_cpm[av])), -0.9999, 0.9999))
                one_m_rho2 = max(1.0 - rho_uv * rho_uv, 1e-8)
                ru = r_cpm_draws[au]
                rv = r_cpm_draws[av]
                beta_u = (ru - rho_uv * rv) / one_m_rho2
                beta_v = (rv - rho_uv * ru) / one_m_rho2
                r2_uv = np.clip(ru * beta_u + rv * beta_v, 0.0, 0.999999)
                se_part = np.sqrt(np.maximum(1.0 - r2_uv, 1e-12) / (float(df_two) * one_m_rho2))
                sgn_b = np.where(r_bal >= 0.0, 1.0, -1.0)
                t_iut = np.minimum(sgn_b * (beta_u / se_part), sgn_b * ((-beta_v) / se_part))
                p_iut = np.clip(
                    np.asarray(2.0 * student_t.sf(t_iut, df=df_two), dtype=np.float64),
                    0.0,
                    1.0,
                )
                p_3leg = np.maximum(p_bal, p_iut)
                emp_p = float(1 + int(np.sum(p_3leg[1:] <= p_3leg[0] + 1e-15))) / float(
                    n_permutations + 1
                )
                perm_map[(u_id, v_id)] = emp_p
            elif reg_str == ContinuousAxis.CO_ABUNDANCE_PRODUCT.value:
                rho_a = float(np.clip(float(np.dot(z_arc[au], z_arc[av])), -0.9999, 0.9999))
                denom_a = max(1.0 - rho_a * rho_a, 1e-8)
                wc_vec = sqrt_used * (z_arc[au] * z_arc[av])
                pu_dot = float(np.dot(wc_vec, z_arc[au]))
                pv_dot = float(np.dot(wc_vec, z_arc[av]))
                gamma_u = (pu_dot - rho_a * pv_dot) / denom_a
                gamma_v = (pv_dot - rho_a * pu_dot) / denom_a
                orth_vec = wc_vec - gamma_u * z_arc[au] - gamma_v * z_arc[av]
                orth_norm = float(np.linalg.norm(orth_vec))
                if orth_norm <= 1e-10:
                    perm_map[(u_id, v_id)] = 1.0
                    continue
                ru_a = r_arc_draws[au]
                rv_a = r_arc_draws[av]
                bu_a = (ru_a - rho_a * rv_a) / denom_a
                bv_a = (rv_a - rho_a * ru_a) / denom_a
                r2_main = np.clip(ru_a * bu_a + rv_a * bv_a, 0.0, 0.99999)
                r_orth_y = (orth_vec @ y_unit) / orth_norm
                part_r = np.clip(
                    r_orth_y / np.sqrt(np.maximum(1.0 - r2_main, 1e-12)),
                    -0.999999,
                    0.999999,
                )
                t_abs = np.abs(part_r) * np.sqrt(
                    float(df_inter) / np.maximum(1.0 - part_r * part_r, 1e-12)
                )
                emp_p = float(1 + int(np.sum(t_abs[1:] >= t_abs[0] - 1e-12))) / float(
                    n_permutations + 1
                )
                perm_map[(u_id, v_id)] = emp_p

    return perm_map


def _canonical_pair_biological_sign(
    cell_type_u: str,
    cell_type_v: str,
    dominant_regime: str,
    dominant_direction: str,
) -> float:
    """Computes a regime-invariant canonical biological sign on `(min(u, v), max(u, v))`."""
    is_swapped = cell_type_u > cell_type_v
    match (is_swapped, dominant_regime):
        case (True, "quadrant:+-"):
            canon_reg = "quadrant:-+"
        case (True, "quadrant:-+"):
            canon_reg = "quadrant:+-"
        case _:
            canon_reg = dominant_regime

    canon_dir = (
        ("non-responder" if dominant_direction == "responder" else "responder")
        if (is_swapped and dominant_regime == ContinuousAxis.LOG_RATIO_BALANCE.value)
        else dominant_direction
    )
    dir_sgn = 1.0 if canon_dir == "responder" else -1.0
    return -dir_sgn if canon_reg in ("quadrant:-+", "quadrant:--") else dir_sgn


def compute_pan_cohort_loco_and_verified_pairs(
    pan_unit: PreparedSuperUnit,
    real_pass: SinglePassPowerBundle,
    global_response_labels: Sequence[str],
    state_prior_weights: Mapping[tuple[str, str], float],
    code_to_readable: Mapping[str, str],
    cfg: PanCohortCliConfig,
) -> Result[tuple[pl.DataFrame, pl.DataFrame], str]:
    """Audits `P6` discoveries with `B = 1,000` stratified permutations and 9-way `LOCO` stability."""
    acat_df = real_pass.pan_p6_acat_df
    if acat_df.is_empty():
        return Success((pl.DataFrame(), pl.DataFrame()))

    # Consolidate ACAT pairs to canonical 1-lead-pair per `(state_a, state_b)` with canonical sign/regime
    is_swapped_state = pl.col("cell_type_u") > pl.col("cell_type_v")
    is_lr_reg = pl.col("dominant_regime") == ContinuousAxis.LOG_RATIO_BALANCE.value
    dedup_acat = (
        acat_df.lazy()
        .filter(pl.col("pair_scope") == PairScope.CROSS_COMPARTMENT.value)
        .with_columns(
            pl.when(is_swapped_state)
            .then(pl.col("cell_type_v"))
            .otherwise(pl.col("cell_type_u"))
            .alias("state_a"),
            pl.when(is_swapped_state)
            .then(pl.col("cell_type_u"))
            .otherwise(pl.col("cell_type_v"))
            .alias("state_b"),
            pl.when(is_swapped_state & (pl.col("dominant_regime") == "quadrant:+-"))
            .then(pl.lit("quadrant:-+"))
            .when(is_swapped_state & (pl.col("dominant_regime") == "quadrant:-+"))
            .then(pl.lit("quadrant:+-"))
            .otherwise(pl.col("dominant_regime"))
            .alias("dominant_regime"),
            pl.when(is_swapped_state & is_lr_reg)
            .then(
                pl.when(pl.col("dominant_direction") == "responder")
                .then(pl.lit("non-responder"))
                .otherwise(pl.lit("responder"))
            )
            .otherwise(pl.col("dominant_direction"))
            .alias("dominant_direction"),
            pl.when(is_swapped_state & is_lr_reg)
            .then(-pl.col("dominant_effect"))
            .otherwise(pl.col("dominant_effect"))
            .alias("dominant_effect"),
        )
        .sort(["pvalue_acat_iut", "min_regime_pvalue_iut"])
        .group_by(["state_a", "state_b"], maintain_order=True)
        .first()
        .collect()
    )

    # Run Leave-One-Cohort-Out (LOCO) across all 9 clinical cohorts for P6 ACAT Z-scores
    uniq_cohorts = sorted(set(pan_unit.cohort_labels))
    coh_arr = np.array(pan_unit.cohort_labels, dtype=str)
    pan_prior = expand_state_weights_to_nhood_pairs(pan_unit.super_meta_codes, state_prior_weights)
    pan_prior_maybe: Maybe[Mapping[tuple[int, int], float]] = (
        Some(pan_prior) if pan_prior else Nothing
    )

    loco_z_by_cohort: dict[str, dict[tuple[int, int], float]] = {}
    for drop_coh in uniq_cohorts:
        keep_cols = np.flatnonzero(coh_arr != drop_coh)
        if int(keep_cols.shape[0]) < 20:
            continue
        sub_unit = PreparedSuperUnit(
            name=f"LOCO_minus_{drop_coh}",
            nhoods_super=pan_unit.nhoods_super,
            counts_super=pan_unit.counts_super[:, keep_cols],
            library_sizes=pan_unit.library_sizes[keep_cols],
            cohort_labels=tuple(pan_unit.cohort_labels[int(i)] for i in keep_cols.tolist()),
            strata_labels=tuple(pan_unit.strata_labels[int(i)] for i in keep_cols.tolist()),
            super_meta_codes=pan_unit.super_meta_codes,
            super_readable_labels=pan_unit.super_readable_labels,
        )
        y_sub = tuple(str(global_response_labels[int(i)]) for i in keep_cols.tolist())
        pw_sub_res = _run_unit_pairwise(
            unit=sub_unit,
            response_labels=y_sub,
            presence_mode=PresenceMode.EMPIRICAL_BAYES_POSTERIOR,
            residualize_depth=True,
            mask_compartments=True,
            dmpp_draws=max(24, cfg.dmpp_null_draws // 2),
            fdr_threshold=cfg.fdr_threshold,
            prior_weights=pan_prior_maybe,
            is_pan_cohort=True,
        )
        match pw_sub_res:
            case Failure(_):
                continue
            case Success(pw_sub):
                acat_sub_res = combine_pair_regimes_acat(
                    quadrant_df=pw_sub.pairs_df,
                    continuous_df=pw_sub.continuous_pairs_df,
                    fdr_threshold=cfg.fdr_threshold,
                    prior_weights=Nothing,
                    n_evaluated_pairs=0,
                )
                match acat_sub_res:
                    case Failure(_):
                        continue
                    case Success(raw_acat_sub):
                        acat_sub = _build_hybrid_prior_acat_df(
                            pw_prior=pw_sub,
                            acat_df=raw_acat_sub,
                            prior_weights=pan_prior_maybe,
                            fdr_threshold=cfg.fdr_threshold,
                        )
                        pair_z_map: dict[tuple[int, int], float] = {}
                        for r_u, r_v, c_u, c_v, r_reg, r_p, r_dir in acat_sub.select(
                            [
                                "nhood_u",
                                "nhood_v",
                                "cell_type_u",
                                "cell_type_v",
                                "dominant_regime",
                                "pvalue_acat_iut",
                                "dominant_direction",
                            ]
                        ).iter_rows():
                            z_mag = (
                                0.0
                                if float(r_p) >= 1.0
                                else float(norm.isf(0.5 * max(float(r_p), 1e-300)))
                            )
                            sgn = _canonical_pair_biological_sign(
                                str(c_u),
                                str(c_v),
                                str(r_reg),
                                str(r_dir),
                            )
                            pair_z_map[(int(r_u), int(r_v))] = sgn * z_mag
                        loco_z_by_cohort[drop_coh] = pair_z_map

    perm_q_map = _compute_exact_stratified_perm_pvalues_for_acat_pairs(
        pan_unit=pan_unit,
        dedup_acat=dedup_acat,
        quad_pairs_df=real_pass.pan_p6_pw.pairs_df,
        global_response_labels=global_response_labels,
        n_permutations=cfg.n_permutations,
        fdr_threshold=cfg.fdr_threshold,
        seed=cfg.seed,
    )

    verified_rows: list[dict[str, object]] = []
    forest_rows: list[dict[str, object]] = []

    for row in dedup_acat.to_dicts():
        u_id = int(row["nhood_u"])
        v_id = int(row["nhood_v"])
        st_a = str(row["state_a"])
        st_b = str(row["state_b"])
        p_acat = float(row["pvalue_acat_iut"])
        fdr_acat = float(row["fdr_acat_iut"])
        dom_dir = str(row["dominant_direction"])
        sgn = _canonical_pair_biological_sign(
            st_a,
            st_b,
            str(row["dominant_regime"]),
            dom_dir,
        )
        z_full = float(norm.isf(0.5 * max(p_acat, 1e-300))) if p_acat < 1.0 else 0.0
        signed_z_full = sgn * z_full

        loco_vals: list[tuple[str, float]] = []
        for drop_coh, z_map in sorted(loco_z_by_cohort.items()):
            if (u_id, v_id) in z_map:
                z_drop = z_map[(u_id, v_id)]
                loco_vals.append((drop_coh, sgn * z_drop))

        min_loco_z = min((zv for _, zv in loco_vals), default=z_full)
        max_loco_z = max((zv for _, zv in loco_vals), default=z_full)
        emp_p = perm_q_map.get((u_id, v_id), p_acat)
        read_a = code_to_readable.get(st_a, st_a)
        read_b = code_to_readable.get(st_b, st_b)
        pair_label = f"{st_a} x {st_b}"
        is_sig = bool(row["is_significant"])

        verified_rows.append(
            {
                "state_a": st_a,
                "state_b": st_b,
                "state_a_readable": read_a,
                "state_b_readable": read_b,
                "compartment_a": classify_tme_compartment(st_a).value,
                "compartment_b": classify_tme_compartment(st_b).value,
                "nhood_u": u_id,
                "nhood_v": v_id,
                "dominant_regime": str(row["dominant_regime"]),
                "dominant_direction": dom_dir,
                "dominant_effect": float(row["dominant_effect"]),
                "z_acat_signed": signed_z_full,
                "z_acat_abs": z_full,
                "pvalue_acat_iut": p_acat,
                "fdr_acat_iut": fdr_acat,
                "emp_perm_pvalue": emp_p,
                "min_loco_z_acat": min_loco_z,
                "max_loco_z_acat": max_loco_z,
                "prior_weight": float(row["prior_weight"]),
                "is_significant": is_sig,
                "is_loco_stable": bool(is_sig and min_loco_z >= 1.64),
            }
        )

        if is_sig and len(forest_rows) < 90:
            for drop_coh, z_aligned in loco_vals:
                forest_rows.append(
                    {
                        "pair_label": pair_label,
                        "pair_detail": f"{read_a[:22]} <-> {read_b[:22]} [{row['dominant_regime']}]",
                        "dominant_direction": dom_dir,
                        "dropped_cohort": drop_coh,
                        "full_z_acat": z_full,
                        "loco_z_acat": z_aligned,
                        "min_loco_z_acat": min_loco_z,
                        "fdr_acat_iut": fdr_acat,
                    }
                )

    verified_df = (
        pl.DataFrame(verified_rows).sort(
            ["is_significant", "pvalue_acat_iut"], descending=[True, False]
        )
        if verified_rows
        else pl.DataFrame()
    )
    loco_forest_df = pl.DataFrame(forest_rows) if forest_rows else pl.DataFrame()
    return Success((verified_df, loco_forest_df))


def build_one_vs_two_stage_concordance_df(
    verified_pan_df: pl.DataFrame,
    meta9_df: pl.DataFrame,
    fdr_threshold: float,
) -> pl.DataFrame:
    """Joins 1-stage `harmonypy` pooled signed Z against 2-stage Random-Effects Meta-IUT signed Z."""
    if verified_pan_df.is_empty() or meta9_df.is_empty():
        return pl.DataFrame(
            schema={
                "state_a": pl.String,
                "state_b": pl.String,
                "pair_label": pl.String,
                "z_one_stage_harmony": pl.Float64,
                "z_two_stage_meta": pl.Float64,
                "fdr_one_stage": pl.Float64,
                "fdr_two_stage": pl.Float64,
                "concordance_tier": pl.String,
            }
        )

    meta_by_regime = (
        meta9_df.lazy()
        .with_columns(
            pl.when(pl.col("direction") == "responder")
            .then(pl.col("z_meta_iut"))
            .otherwise(-pl.col("z_meta_iut"))
            .alias("z_two_stage_signed")
        )
        .sort("pvalue_meta_iut")
        .group_by(["state_a", "state_b", "canonical_regime"], maintain_order=True)
        .first()
        .select(
            [
                "state_a",
                "state_b",
                pl.col("canonical_regime").alias("dominant_regime"),
                "z_two_stage_signed",
                "fdr_meta_iut",
                "n_cohorts_observed",
                "min_loco_z_iut",
            ]
        )
        .collect()
    )

    joined = verified_pan_df.join(
        meta_by_regime, on=["state_a", "state_b", "dominant_regime"], how="inner"
    )
    if joined.is_empty():
        return pl.DataFrame()

    return (
        joined.lazy()
        .with_columns(
            (pl.col("state_a") + pl.lit(" x ") + pl.col("state_b")).alias("pair_label"),
            pl.col("z_acat_signed").alias("z_one_stage_harmony"),
            pl.col("z_two_stage_signed").alias("z_two_stage_meta"),
            pl.col("fdr_acat_iut").alias("fdr_one_stage"),
            pl.col("fdr_meta_iut").alias("fdr_two_stage"),
            pl.when(
                (pl.col("fdr_acat_iut") < fdr_threshold) & (pl.col("fdr_meta_iut") < fdr_threshold)
            )
            .then(pl.lit("Both Sig (1-Stage & 2-Stage)"))
            .when(pl.col("fdr_acat_iut") < fdr_threshold)
            .then(pl.lit("1-Stage Harmony Only"))
            .when(pl.col("fdr_meta_iut") < fdr_threshold)
            .then(pl.lit("2-Stage Meta-IUT Only"))
            .otherwise(pl.lit("Background Pair"))
            .alias("concordance_tier"),
        )
        .select(
            [
                "state_a",
                "state_b",
                "pair_label",
                "z_one_stage_harmony",
                "z_two_stage_meta",
                "fdr_one_stage",
                "fdr_two_stage",
                "n_cohorts_observed",
                "concordance_tier",
            ]
        )
        .collect()
    )


def compute_stratified_pan_univariate_da(
    pan_unit: PreparedSuperUnit,
    response_labels: Sequence[str],
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Computes observable-cohort stratified 1D univariate BH FDR and logFC per super-neighborhood."""
    pres = estimate_nhood_presence(
        count_matrix=pan_unit.counts_super,
        library_sizes=Some(pan_unit.library_sizes),
        config=PresenceConfig(
            mode=PresenceMode.EMPIRICAL_BAYES_POSTERIOR,
            residualize_library_depth=True,
        ),
    ).unwrap()
    pres_m = build_compartment_observability_mask(
        presence=pres,
        nhood_labels=pan_unit.super_meta_codes,
        cohort_labels=pan_unit.cohort_labels,
        min_cohort_compartment_share=0.001,
    ).value_or(pres)

    v_super = int(pan_unit.counts_super.shape[0])
    y_bin = np.array([1.0 if r == "responder" else 0.0 for r in response_labels], dtype=np.float64)
    strata_arr = np.array(pan_unit.strata_labels, dtype=str)
    uniq_st = sorted(set(strata_arr.tolist()))

    pvals = np.ones(v_super, dtype=np.float64)
    lfcs = np.zeros(v_super, dtype=np.float64)

    for row_i in range(v_super):
        obs_m = pres_m.observability[row_i, :] > 0.5
        if int(np.sum(obs_m)) < 12:
            continue
        x_row = np.array(pres_m.log2_cpm[row_i, obs_m], dtype=np.float64, copy=True)
        y_sub = y_bin[obs_m]
        st_sub = strata_arr[obs_m]
        # Demean within each cohort::substudy stratum
        for st_name in uniq_st:
            sm = st_sub == st_name
            if int(np.sum(sm)) >= 2:
                x_row[sm] = x_row[sm] - float(np.mean(x_row[sm]))
        r_vals = x_row[y_sub > 0.5]
        nr_vals = x_row[y_sub <= 0.5]
        if int(r_vals.shape[0]) >= 4 and int(nr_vals.shape[0]) >= 4:
            diff = float(np.mean(r_vals) - np.mean(nr_vals))
            v_r = float(np.var(r_vals, ddof=1)) / float(r_vals.shape[0])
            v_nr = float(np.var(nr_vals, ddof=1)) / float(nr_vals.shape[0])
            se = float(np.sqrt(max(v_r + v_nr, 1e-10)))
            z_u = abs(diff) / se
            pvals[row_i] = float(np.clip(2.0 * norm.sf(z_u), 1e-300, 1.0))
            lfcs[row_i] = diff

    fdr_u = compute_weighted_bh_fdr(pvals, np.ones_like(pvals))
    return tuple(float(x) for x in fdr_u.tolist()), tuple(float(x) for x in lfcs.tolist())


def build_pan_cohort_umap_bundle_from_p6(
    pan_unit: PreparedSuperUnit,
    umap_coords: FloatArray,
    anchor_cell_labels: Sequence[str],
    global_response_labels: Sequence[str],
    real_pass: SinglePassPowerBundle,
    fdr_threshold: float,
    max_chords: int = 28,
) -> Result[PairwiseUmapBundle, str]:
    """Projects `P6` omnibus ACAT discoveries onto the unified 9-cohort UMAP with observable CIs."""
    acat_df = real_pass.pan_p6_acat_df
    base_pw = real_pass.pan_p6_pw
    read_lbls = pan_unit.super_readable_labels

    # Select at most 1 lead ACAT discovery per canonical `(state_a, state_b)` pair
    sig_acat = (
        acat_df.lazy()
        .filter(pl.col("is_significant") & (pl.col("pair_scope") == PairScope.CROSS_COMPARTMENT.value))
        .with_columns(
            pl.when(pl.col("cell_type_u") < pl.col("cell_type_v"))
            .then(pl.col("cell_type_u"))
            .otherwise(pl.col("cell_type_v"))
            .alias("_sa"),
            pl.when(pl.col("cell_type_u") < pl.col("cell_type_v"))
            .then(pl.col("cell_type_v"))
            .otherwise(pl.col("cell_type_u"))
            .alias("_sb"),
        )
        .sort(["pvalue_acat_iut", "min_regime_pvalue_iut"])
        .group_by(["_sa", "_sb"], maintain_order=True)
        .first()
        .collect()
    )

    q_keep_map: dict[tuple[int, int, str], tuple[float, float]] = {}
    c_keep_map: dict[tuple[int, int, str], tuple[float, float]] = {}
    for r in sig_acat.to_dicts():
        uu, vv = int(r["nhood_u"]), int(r["nhood_v"])
        dom_reg = str(r["dominant_regime"])
        p_ac, q_ac = float(r["pvalue_acat_iut"]), float(r["fdr_acat_iut"])
        if dom_reg.startswith("quadrant:"):
            q_code = dom_reg.split(":", 1)[1]
            q_keep_map[(uu, vv, q_code)] = (p_ac, q_ac)
        else:
            c_keep_map[(uu, vv, dom_reg)] = (p_ac, q_ac)

    # Update pairs_df and continuous_pairs_df cleanly so `prepare_pairwise_umap_bundle` renders exact P6 chords
    if not base_pw.pairs_df.is_empty():
        q_rows_new: list[dict[str, object]] = []
        for qr in base_pw.pairs_df.to_dicts():
            u_i, v_i = int(qr["nhood_u"]), int(qr["nhood_v"])
            q_str = str(qr["quadrant"])
            if u_i > v_i:
                q_can = "-+" if q_str == "+-" else ("+-" if q_str == "-+" else q_str)
                k_uv = (v_i, u_i, q_can)
            else:
                k_uv = (u_i, v_i, q_str)
            hit = k_uv in q_keep_map
            p_new, f_new = q_keep_map[k_uv] if hit else (float(qr["pvalue_iut"]), max(float(qr["fdr_iut"]), 0.50))
            updated = dict(qr)
            updated["cell_type_u"] = read_lbls[u_i] if u_i < len(read_lbls) else str(qr["cell_type_u"])
            updated["cell_type_v"] = read_lbls[v_i] if v_i < len(read_lbls) else str(qr["cell_type_v"])
            updated["pvalue_iut"] = p_new
            updated["fdr_iut"] = f_new
            updated["is_significant"] = hit
            q_rows_new.append(updated)
        new_q_df = pl.DataFrame(q_rows_new)
    else:
        new_q_df = base_pw.pairs_df

    if not base_pw.continuous_pairs_df.is_empty():
        c_rows_new: list[dict[str, object]] = []
        for cr in base_pw.continuous_pairs_df.to_dicts():
            u_i, v_i = int(cr["nhood_u"]), int(cr["nhood_v"])
            ax_str = str(cr["axis"])
            k_uv = (min(u_i, v_i), max(u_i, v_i), ax_str)
            hit = k_uv in c_keep_map
            p_new, f_new = c_keep_map[k_uv] if hit else (float(cr["pvalue_iut"]), max(float(cr["fdr_iut"]), 0.50))
            updated = dict(cr)
            updated["cell_type_u"] = read_lbls[u_i] if u_i < len(read_lbls) else str(cr["cell_type_u"])
            updated["cell_type_v"] = read_lbls[v_i] if v_i < len(read_lbls) else str(cr["cell_type_v"])
            updated["pvalue_iut"] = p_new
            updated["fdr_iut"] = f_new
            updated["is_significant"] = hit
            c_rows_new.append(updated)
        new_c_df = pl.DataFrame(c_rows_new)
    else:
        new_c_df = base_pw.continuous_pairs_df

    n_sig_q = int(new_q_df.filter(pl.col("is_significant")).height) if not new_q_df.is_empty() else 0
    n_sig_c = int(new_c_df.filter(pl.col("is_significant")).height) if not new_c_df.is_empty() else 0
    pw_viz = PairwiseAnalysisResult(
        pairs_df=new_q_df,
        continuous_pairs_df=new_c_df,
        modules_df=pl.DataFrame(),
        n_total_nhoods=base_pw.n_total_nhoods,
        n_active_nhoods=base_pw.n_active_nhoods,
        n_evaluated_pairs=base_pw.n_evaluated_pairs,
        n_significant_quadrant_pairs=n_sig_q,
        n_significant_continuous_pairs=n_sig_c,
        n_modules=n_sig_q + n_sig_c,
    )

    univ_fdr, univ_lfc = compute_stratified_pan_univariate_da(pan_unit, global_response_labels)
    viz_cfg = UmapPairVizConfig(
        max_chords=max_chords,
        max_bg_cells=12000,
        fdr_threshold=fdr_threshold,
    )
    base_bundle_res = prepare_pairwise_umap_bundle(
        nhoods_mat=pan_unit.nhoods_super,
        umap_coords=umap_coords,
        count_matrix=pan_unit.counts_super,
        response_labels=global_response_labels,
        pairwise_res=pw_viz,
        sample_library_sizes=Some(tuple(pan_unit.library_sizes.tolist())),
        state_labels=Some(read_lbls),
        cell_labels=Some(anchor_cell_labels),
        univariate_fdr=Some(univ_fdr),
        univariate_logfc=Some(univ_lfc),
        config=viz_cfg,
    )
    match base_bundle_res:
        case Failure(err):
            return Failure(err)
        case Success(bundle):
            pass

    # Recompute spotlight 2x2 donor inset strictly over observable biopsies for `(u*, v*)`
    if not bundle.spotlight_meta.has_pair:
        return Success(bundle)

    u_star = bundle.spotlight_meta.nhood_u
    v_star = bundle.spotlight_meta.nhood_v
    pres_all = estimate_nhood_presence(
        count_matrix=pan_unit.counts_super,
        library_sizes=Some(pan_unit.library_sizes),
        config=PresenceConfig(mode=PresenceMode.EMPIRICAL_BAYES_POSTERIOR),
    ).unwrap()
    pres_masked = build_compartment_observability_mask(
        presence=pres_all,
        nhood_labels=pan_unit.super_meta_codes,
        cohort_labels=pan_unit.cohort_labels,
        min_cohort_compartment_share=0.001,
    ).value_or(pres_all)

    both_obs = (pres_masked.observability[u_star, :] > 0.5) & (
        pres_masked.observability[v_star, :] > 0.5
    )
    obs_idx = np.flatnonzero(both_obs)
    if int(obs_idx.shape[0]) < 8:
        return Success(bundle)

    u_p = pres_masked.presence_prob[u_star, obs_idx] >= 0.5
    v_p = pres_masked.presence_prob[v_star, obs_idx] >= 0.5
    if int(np.sum(u_p & v_p)) == int(obs_idx.shape[0]) or int(np.sum(u_p | v_p)) == 0:
        u_cpm = pres_masked.log2_cpm[u_star, obs_idx]
        v_cpm = pres_masked.log2_cpm[v_star, obs_idx]
        u_p = u_cpm > float(np.median(u_cpm))
        v_p = v_cpm > float(np.median(v_cpm))

    y_obs = np.array(
        [1 if global_response_labels[int(j)] == "responder" else 0 for j in obs_idx.tolist()],
        dtype=np.int64,
    )
    q_defs = (
        ("u+ / v+", u_p & v_p),
        ("u+ / v-", u_p & (~v_p)),
        ("u- / v+", (~u_p) & v_p),
        ("u- / v-", (~u_p) & (~v_p)),
    )
    q_rows_obs: list[dict[str, object]] = []
    for idx_q, (q_lbl, q_m) in enumerate(q_defs):
        n_t = int(np.sum(q_m))
        n_r = int(np.sum(y_obs[q_m] == 1)) if n_t > 0 else 0
        pct = 100.0 * float(n_r) / float(max(n_t, 1))
        lo, hi = wilson_binomial_ci_pct(n_r, n_t)
        q_rows_obs.append(
            {
                "donor_quadrant": q_lbl,
                "quad_order": idx_q,
                "n_responders": n_r,
                "n_total": n_t,
                "responder_pct": pct,
                "ci_low_pct": lo,
                "ci_high_pct": hi,
            }
        )

    return Success(
        PairwiseUmapBundle(
            umap_bounds=bundle.umap_bounds,
            nhood_nodes_df=bundle.nhood_nodes_df,
            hub_overlay_df=bundle.hub_overlay_df,
            bg_cells_df=bundle.bg_cells_df,
            bundled_arcs_df=bundle.bundled_arcs_df,
            arc_endpoints_df=bundle.arc_endpoints_df,
            anchor_label=bundle.anchor_label,
            anchor_bg_cells_df=bundle.anchor_bg_cells_df,
            anchor_overlay_df=bundle.anchor_overlay_df,
            anchor_arcs_df=bundle.anchor_arcs_df,
            spotlight_cells_df=bundle.spotlight_cells_df,
            spotlight_arc_df=bundle.spotlight_arc_df,
            spotlight_endpoints_df=bundle.spotlight_endpoints_df,
            spotlight_quadrant_df=pl.DataFrame(q_rows_obs),
            spotlight_meta=SpotlightPairMetadata(
                has_pair=bundle.spotlight_meta.has_pair,
                is_exploratory=bundle.spotlight_meta.is_exploratory,
                nhood_u=bundle.spotlight_meta.nhood_u,
                nhood_v=bundle.spotlight_meta.nhood_v,
                label_u=bundle.spotlight_meta.label_u,
                label_v=bundle.spotlight_meta.label_v,
                regime_label=bundle.spotlight_meta.regime_label,
                direction_label=bundle.spotlight_meta.direction_label,
                novelty_label=bundle.spotlight_meta.novelty_label,
                fdr_iut=bundle.spotlight_meta.fdr_iut,
                pvalue_iut=bundle.spotlight_meta.pvalue_iut,
                effect_value=bundle.spotlight_meta.effect_value,
                univariate_fdr_u=bundle.spotlight_meta.univariate_fdr_u,
                univariate_fdr_v=bundle.spotlight_meta.univariate_fdr_v,
                n_cells_u=bundle.spotlight_meta.n_cells_u,
                n_cells_v=bundle.spotlight_meta.n_cells_v,
            ),
            n_sig_quadrant=bundle.n_sig_quadrant,
            n_sig_continuous=bundle.n_sig_continuous,
            n_denovo_nhoods=bundle.n_denovo_nhoods,
        )
    )


# =============================================================================
# 4. Pure Functional Core — Declarative Nature Methods Altair Composites
# =============================================================================

_NINE_COHORT_PALETTE: tuple[str, ...] = (
    OKABE_BLUE,
    OKABE_VERMILION,
    OKABE_BLUISH_GREEN,
    OKABE_ORANGE,
    OKABE_REDDISH_PURPLE,
    OKABE_SKY_BLUE,
    "#334155",
    "#CA8A04",
    "#0F766E",
)

_COMPARTMENT_READABLE_SHORT: Mapping[str, str] = {
    TmeCompartment.T_NK.value: "T / NK",
    TmeCompartment.B_LINEAGE.value: "B Lineage",
    TmeCompartment.PLASMA.value: "Plasma",
    TmeCompartment.MYELOID_MAST.value: "Myeloid / Mast",
    TmeCompartment.STROMAL_MESENCHYMAL.value: "Stromal / CAF",
    TmeCompartment.ENDOTHELIAL_VASCULAR.value: "Endothelial",
    TmeCompartment.MALIGNANT_EPITHELIAL.value: "Epithelial / Tumor",
    TmeCompartment.OTHER.value: "Other",
}


def compute_pan_umap_2d(
    embedding: FloatArray,
    max_fit_cells: int = 16000,
    seed: int = 42,
) -> FloatArray:
    """Computes deterministic 2D UMAP coordinates for `embedding` (`fit` + `transform` if large)."""
    emb = np.asarray(embedding, dtype=np.float64)
    n_cells = int(emb.shape[0])
    reducer = umap.UMAP(
        n_neighbors=30,
        min_dist=0.30,
        n_components=2,
        metric="euclidean",
        random_state=seed,
        n_jobs=1,
    )
    if n_cells <= max_fit_cells:
        coords = np.asarray(reducer.fit_transform(emb), dtype=np.float64)
        return coords

    rng = np.random.default_rng(seed)
    fit_idx = np.sort(rng.choice(n_cells, size=max_fit_cells, replace=False))
    reducer.fit(emb[fit_idx])
    coords_all = np.asarray(reducer.transform(emb), dtype=np.float64)
    return coords_all


def build_harmony_integration_audit_figure(
    bundle: PanCohortFeatureBundle,
    umap_i0: FloatArray,
    umap_sel: FloatArray,
    selected_regime: HarmonyRegime,
    harmony_bench_df: pl.DataFrame,
    atlas_i0: PanCohortMiloAtlas,
    atlas_sel: PanCohortMiloAtlas,
    seed: int = 42,
) -> AltairComposite:
    """Builds the 2x2 Nature Methods `pan_cohort_harmony_integration_audit` composite."""
    n_anchor = int(umap_i0.shape[0])
    rng = np.random.default_rng(seed)
    plot_n = min(n_anchor, 10000)
    pick = (
        rng.permutation(n_anchor)[:plot_n]
        if n_anchor > plot_n
        else rng.permutation(n_anchor)
    )

    coh_short_map = {
        p["acc"]: f"{p['acc']} ({p['indication']})" for p in COHORT_PANELS
    }
    coh_order = [f"{p['acc']} ({p['indication']})" for p in COHORT_PANELS]
    coh_colors = list(_NINE_COHORT_PALETTE[: len(coh_order)])

    cells_plot_df = pl.DataFrame(
        {
            "u0_1": umap_i0[pick, 0],
            "u0_2": umap_i0[pick, 1],
            "uh_1": umap_sel[pick, 0],
            "uh_2": umap_sel[pick, 1],
            "cohort_tag": [
                coh_short_map.get(bundle.anchor_cohort_ids[int(i)], bundle.anchor_cohort_ids[int(i)])
                for i in pick.tolist()
            ],
            "compartment": [
                _COMPARTMENT_READABLE_SHORT.get(
                    bundle.anchor_compartments[int(i)], bundle.anchor_compartments[int(i)]
                )
                for i in pick.tolist()
            ],
        }
    )

    # Compartment geometric medoids on the harmonized UMAP for crisp biological callouts
    med_rows: list[dict[str, object]] = []
    comp_full_arr = np.array(bundle.anchor_compartments, dtype=str)
    for cp_enum in TmeCompartment:
        if cp_enum == TmeCompartment.OTHER:
            continue
        cp_idx = np.flatnonzero(comp_full_arr == cp_enum.value)
        if int(cp_idx.shape[0]) < 25:
            continue
        sub_xy = umap_sel[cp_idx]
        med_xy = np.median(sub_xy, axis=0)
        best_local = int(np.argmin(np.sum(np.square(sub_xy - med_xy), axis=1)))
        med_rows.append(
            {
                "uh_1": float(sub_xy[best_local, 0]),
                "uh_2": float(sub_xy[best_local, 1]),
                "comp_label": _COMPARTMENT_READABLE_SHORT[cp_enum.value],
            }
        )
    comp_med_df = pl.DataFrame(med_rows)

    m_i0 = harmony_bench_df.filter(pl.col("regime") == HarmonyRegime.I0_UNINTEGRATED_PCA.value)
    m_sel = harmony_bench_df.filter(pl.col("regime") == selected_regime.value)
    ilisi_i0 = float(m_i0["cohort_ilisi"].item()) if not m_i0.is_empty() else 1.0
    ilisi_sel = float(m_sel["cohort_ilisi"].item()) if not m_sel.is_empty() else 1.0
    clisi_sel = float(m_sel["compartment_clisi"].item()) if not m_sel.is_empty() else 1.0

    panel_a = (
        alt.Chart(cells_plot_df)
        .mark_circle(size=5.5, opacity=0.62)
        .encode(
            x=alt.X("u0_1:Q", axis=create_nature_axis("Unintegrated UMAP 1")),
            y=alt.Y("u0_2:Q", axis=create_nature_axis("Unintegrated UMAP 2")),
            color=alt.Color(
                "cohort_tag:N",
                sort=coh_order,
                scale=alt.Scale(domain=coh_order, range=coh_colors),
                legend=alt.Legend(
                    title="Clinical Cohort (9 Studies)",
                    orient="right",
                    symbolSize=55,
                    symbolOpacity=0.95,
                    labelLimit=185,
                ),
            ),
        )
        .properties(
            width=360,
            height=300,
            title=alt.TitleParams(
                text="a  Unintegrated Concatenated 9-Cohort PCA UMAP (I0)",
                subtitle=f"Cohort iLISI = {ilisi_i0:.2f} | Median {atlas_i0.median_cohorts_per_nhood:.0f} cohort / Milo nhood",
                anchor="start",
                fontSize=11.5,
                subtitleFontSize=9.5,
                subtitleColor="#475569",
            ),
        )
    )

    base_b = (
        alt.Chart(cells_plot_df)
        .mark_circle(size=5.5, opacity=0.62)
        .encode(
            x=alt.X("uh_1:Q", axis=create_nature_axis("Harmony UMAP 1")),
            y=alt.Y("uh_2:Q", axis=create_nature_axis("Harmony UMAP 2")),
            color=alt.Color(
                "cohort_tag:N",
                sort=coh_order,
                scale=alt.Scale(domain=coh_order, range=coh_colors),
            ),
        )
    )
    med_pts_b = (
        alt.Chart(comp_med_df)
        .mark_point(shape="diamond", size=75, filled=True, fill="#FFFFFF", stroke="#0F172A", strokeWidth=1.4)
        .encode(x="uh_1:Q", y="uh_2:Q")
    )
    med_txt_b = (
        alt.Chart(comp_med_df)
        .mark_text(
            dy=-10,
            fontSize=8.8,
            fontWeight="bold",
            fill="#0F172A",
            stroke="#FFFFFF",
            strokeWidth=2.4,
        )
        .encode(x="uh_1:Q", y="uh_2:Q", text="comp_label:N")
    )
    med_txt_fg_b = (
        alt.Chart(comp_med_df)
        .mark_text(dy=-10, fontSize=8.8, fontWeight="bold", fill="#0F172A")
        .encode(x="uh_1:Q", y="uh_2:Q", text="comp_label:N")
    )
    panel_b = (
        alt.layer(base_b, med_pts_b, med_txt_b, med_txt_fg_b)
        .properties(
            width=360,
            height=300,
            title=alt.TitleParams(
                text=f"b  {harmony_regime_display_label(selected_regime)} UMAP",
                subtitle=(
                    f"Cohort iLISI = {ilisi_sel:.2f} | Compartment cLISI = {clisi_sel:.2f} | "
                    f"Median {atlas_sel.median_cohorts_per_nhood:.0f} cohorts / nhood"
                ),
                anchor="start",
                fontSize=11.5,
                subtitleFontSize=9.5,
                subtitleColor="#475569",
            ),
        )
    )

    # Panel c: Pareto trade-off across I0..I3
    pareto_df = harmony_bench_df.with_columns(
        (pl.col("kbet_acceptance_rate") * 100.0).alias("kbet_pct"),
        pl.when(pl.col("regime") == selected_regime.value)
        .then(pl.lit("Selected Pareto Optimum"))
        .otherwise(pl.lit("Benchmark Candidate"))
        .alias("sel_status"),
    )
    reg_order = [harmony_regime_display_label(r) for r in ALL_HARMONY_REGIMES]
    reg_colors = [OKABE_GREY, OKABE_SKY_BLUE, OKABE_BLUE, OKABE_ORANGE]

    clisi_rule = (
        alt.Chart(pl.DataFrame({"clisi_cap": [1.15]}))
        .mark_rule(strokeDash=[4, 3], stroke="#64748B", strokeWidth=1.0)
        .encode(y="clisi_cap:Q")
    )
    pareto_pts = (
        alt.Chart(pareto_df)
        .mark_point(filled=True, opacity=0.95, stroke="#0F172A", strokeWidth=1.1)
        .encode(
            x=alt.X(
                "cohort_ilisi:Q",
                scale=alt.Scale(zero=False, nice=True),
                axis=create_nature_axis("Batch Mixing: Median Cohort iLISI (Higher = Better)"),
            ),
            y=alt.Y(
                "compartment_clisi:Q",
                scale=alt.Scale(
                    domain=[
                        0.98,
                        max(1.22, float(np.max(pareto_df["compartment_clisi"].to_numpy())) + 0.04),
                    ]
                ),
                axis=create_nature_axis("Bio Purity: Compartment cLISI (1.00 = Ideal)"),
            ),
            size=alt.Size(
                "kbet_pct:Q",
                scale=alt.Scale(range=[110, 340]),
                legend=alt.Legend(title="kBET Accept (%)", orient="right"),
            ),
            color=alt.Color(
                "regime_label:N",
                sort=reg_order,
                scale=alt.Scale(domain=reg_order, range=reg_colors),
                legend=alt.Legend(title="Latent Space", orient="right", labelLimit=185),
            ),
        )
    )
    pareto_lbls = (
        alt.Chart(pareto_df)
        .mark_text(dy=-13, fontSize=8.5, fontWeight="bold", fill="#1E293B")
        .encode(x="cohort_ilisi:Q", y="compartment_clisi:Q", text="short_tag:N")
    )
    panel_c = (
        alt.layer(clisi_rule, pareto_pts, pareto_lbls)
        .properties(
            width=360,
            height=270,
            title=alt.TitleParams(
                text="c  Harmony Integration vs. Compartment Purity Pareto Audit",
                subtitle="Dashed ceiling: cLISI = 1.15 bio-conservation threshold across 371,568 cells",
                anchor="start",
                fontSize=11.5,
                subtitleFontSize=9.5,
                subtitleColor="#475569",
            ),
        )
    )

    # Panel d: Post-Milo neighborhood distinct cohort count distribution (I0 vs Selected Harmony)
    dist_rows: list[dict[str, object]] = []
    for tag_lbl, mix_df in (
        ("I0: Unintegrated PCA", atlas_i0.nhood_mixing_df),
        (harmony_regime_display_label(selected_regime), atlas_sel.nhood_mixing_df),
    ):
        tot_nh = max(mix_df.height, 1)
        cnt_arr = mix_df["n_distinct_cohorts"].to_numpy()
        for k_c in range(1, 10):
            pct_k = 100.0 * float(np.sum(cnt_arr == k_c)) / float(tot_nh)
            dist_rows.append(
                {
                    "n_cohorts_bin": str(k_c),
                    "n_cohorts_num": k_c,
                    "nhood_pct": pct_k,
                    "space_label": tag_lbl,
                }
            )
    dist_df = pl.DataFrame(dist_rows)
    space_dom = ["I0: Unintegrated PCA", harmony_regime_display_label(selected_regime)]

    panel_d = (
        alt.Chart(dist_df)
        .mark_bar(stroke="#1E293B", strokeWidth=0.45)
        .encode(
            x=alt.X(
                "n_cohorts_bin:O",
                sort=[str(i) for i in range(1, 10)],
                axis=create_nature_axis("Distinct Clinical Cohorts Represented in Milo Neighborhood"),
            ),
            xOffset=alt.XOffset("space_label:N", sort=space_dom),
            y=alt.Y(
                "nhood_pct:Q",
                axis=create_nature_axis("Share of Pan-Cohort Milo Neighborhoods (%)"),
            ),
            color=alt.Color(
                "space_label:N",
                sort=space_dom,
                scale=alt.Scale(domain=space_dom, range=[OKABE_GREY, OKABE_BLUISH_GREEN]),
                legend=alt.Legend(title="Milo kNN Graph", orient="right", labelLimit=185),
            ),
        )
        .properties(
            width=360,
            height=270,
            title=alt.TitleParams(
                text="d  Post-Milo Neighborhood Multi-Cohort Mixing Shift",
                subtitle=(
                    f"Harmony raises Donor LISI {atlas_i0.median_donor_lisi:.1f} -> "
                    f"{atlas_sel.median_donor_lisi:.1f} & Mixing-QC pass {atlas_sel.passes_mixing_qc_pct:.1f}%"
                ),
                anchor="start",
                fontSize=11.5,
                subtitleFontSize=9.5,
                subtitleColor="#475569",
            ),
        )
    )

    top_row = alt.hconcat(panel_a, panel_b, spacing=26).resolve_scale(color="shared")
    bot_row = alt.hconcat(panel_c, panel_d, spacing=26).resolve_scale(color="independent", size="independent")
    composite = (
        alt.vconcat(top_row, bot_row, spacing=28)
        .resolve_scale(color="independent")
        .properties(
            title=alt.TitleParams(
                text="Pan-Cohort Single-Cell Integration Benchmark Across All Nine Clinical ICB Datasets",
                subtitle=(
                    f"N = {len(bundle.eval_sample_uids)} biopsy evaluations (304R / 309NR) | "
                    f"{bundle.n_total_cells_all_cohorts:,} total cells | {len(set(bundle.eval_substudy_strata))} cohort/sub-study strata"
                ),
                anchor="start",
                fontSize=13.5,
                subtitleFontSize=10.5,
                subtitleColor="#334155",
                offset=12,
            )
        )
    )
    return cast(AltairComposite, apply_nature_methods_theme(composite))


def build_power_ablation_and_calibration_figure(
    power_ablation_df: pl.DataFrame,
    per_cohort_uplift_df: pl.DataFrame,
    concordance_df: pl.DataFrame,
    loco_forest_df: pl.DataFrame,
    fdr_threshold: float,
) -> AltairComposite:
    """Builds the 2x2 Nature Methods `statistical_power_ablation_and_null_calibration` composite."""
    tier_order = [power_tier_display_label(t) for t in ALL_POWER_TIERS]
    cond_order = ["Real Clinical Outcomes", "Permuted-Y Null Control (Mean)"]

    ab_long_rows: list[dict[str, object]] = []
    for r in power_ablation_df.to_dicts():
        t_lbl = str(r["tier_label"])
        ab_long_rows.append(
            {
                "tier_label": t_lbl,
                "condition": "Real Clinical Outcomes",
                "n_discoveries": float(r["real_sig_state_pairs"]),
                "label_txt": f"{int(r['real_sig_state_pairs'])}",
            }
        )
        ab_long_rows.append(
            {
                "tier_label": t_lbl,
                "condition": "Permuted-Y Null Control (Mean)",
                "n_discoveries": float(r["null_mean_sig_state_pairs"]),
                "label_txt": f"{float(r['null_mean_sig_state_pairs']):.1f}",
            }
        )
    ab_long_df = pl.DataFrame(ab_long_rows)

    bars_a = (
        alt.Chart(ab_long_df)
        .mark_bar(stroke="#1E293B", strokeWidth=0.45)
        .encode(
            y=alt.Y("tier_label:N", sort=tier_order, axis=create_nature_axis("")),
            yOffset=alt.YOffset("condition:N", sort=cond_order),
            x=alt.X(
                "n_discoveries:Q",
                axis=create_nature_axis(
                    f"Significant Cross-Compartment Meta42 Pairs (FDR < {fdr_threshold:.2f})"
                ),
            ),
            color=alt.Color(
                "condition:N",
                sort=cond_order,
                scale=alt.Scale(domain=cond_order, range=[OKABE_BLUE, OKABE_VERMILION]),
                legend=alt.Legend(title="Outcome Regime", orient="top", labelLimit=210),
            ),
        )
    )
    txt_a = (
        alt.Chart(ab_long_df)
        .mark_text(align="left", dx=4, fontSize=8.3, fontWeight="bold", fill="#1E293B")
        .encode(
            y=alt.Y("tier_label:N", sort=tier_order),
            yOffset=alt.YOffset("condition:N", sort=cond_order),
            x=alt.X("n_discoveries:Q"),
            text="label_txt:N",
        )
    )
    panel_a = (
        alt.layer(bars_a, txt_a)
        .properties(
            width=370,
            height=285,
            title=alt.TitleParams(
                text="a  Seven-Tier Statistical Power Progression vs. Permuted-Response Null",
                subtitle="Within-stratum outcome permutations confirm Type-I error calibration across P0 -> P6",
                anchor="start",
                fontSize=11.5,
                subtitleFontSize=9.5,
                subtitleColor="#475569",
            ),
        )
    )

    # Panel b: Per-cohort discovery uplift across all 9 clinical cohorts
    coh_sort = per_cohort_uplift_df.sort("n_evals", descending=True)["cohort_short"].to_list()
    up_rows: list[dict[str, object]] = []
    for r in per_cohort_uplift_df.to_dicts():
        cs = str(r["cohort_short"])
        up_rows.append(
            {
                "cohort_short": cs,
                "stage": "P0: Baseline (Hard + B=12)",
                "sig_pairs": int(r["p0_baseline_pairs"]),
            }
        )
        up_rows.append(
            {
                "cohort_short": cs,
                "stage": "P2: + EB Soft + Depth + B=64 + ACAT",
                "sig_pairs": int(r["p2_boosted_pairs"]),
            }
        )
    up_long_df = pl.DataFrame(up_rows)
    st_order = ["P0: Baseline (Hard + B=12)", "P2: + EB Soft + Depth + B=64 + ACAT"]

    rule_b = (
        alt.Chart(per_cohort_uplift_df)
        .mark_rule(stroke="#94A3B8", strokeWidth=1.5)
        .encode(
            y=alt.Y("cohort_short:N", sort=coh_sort, axis=create_nature_axis("")),
            x=alt.X(
                "p0_baseline_pairs:Q",
                axis=create_nature_axis("Within-Cohort Cross-Compartment Pairs (q < 0.10)"),
            ),
            x2="p2_boosted_pairs:Q",
        )
    )
    pts_b = (
        alt.Chart(up_long_df)
        .mark_point(filled=True, size=95, stroke="#0F172A", strokeWidth=0.9)
        .encode(
            y=alt.Y("cohort_short:N", sort=coh_sort),
            x=alt.X("sig_pairs:Q"),
            color=alt.Color(
                "stage:N",
                sort=st_order,
                scale=alt.Scale(domain=st_order, range=[OKABE_GREY, OKABE_BLUISH_GREEN]),
                legend=alt.Legend(title="Single-Cohort Engine", orient="top", labelLimit=220),
            ),
        )
    )
    panel_b = (
        alt.layer(rule_b, pts_b)
        .properties(
            width=370,
            height=285,
            title=alt.TitleParams(
                text="b  Within-Cohort Power Uplift Across Individual Datasets",
                subtitle="Small cohorts (J = 7..39) still require cross-cohort pooling (P3 / P6) for rare arms",
                anchor="start",
                fontSize=11.5,
                subtitleFontSize=9.5,
                subtitleColor="#475569",
            ),
        )
    )

    # Panel c: 1-Stage Harmony Pooled Z vs 2-Stage Random-Effects Meta-IUT Z
    if not concordance_df.is_empty() and concordance_df.height >= 3:
        x_z = concordance_df["z_one_stage_harmony"].to_numpy()
        y_z = concordance_df["z_two_stage_meta"].to_numpy()
        r_val = float(pearsonr(x_z, y_z).statistic) if float(np.std(x_z)) > 1e-9 and float(np.std(y_z)) > 1e-9 else 0.0
        rho_val = float(spearmanr(x_z, y_z).statistic) if float(np.std(x_z)) > 1e-9 and float(np.std(y_z)) > 1e-9 else 0.0
        conc_sub = f"Pearson r = {r_val:.2f} | Spearman ρ = {rho_val:.2f} across {concordance_df.height} cross-compartment state pairs"
    else:
        conc_sub = "Concordance across cross-compartment Meta42 state pairs"

    tier_dom_c = [
        "Both Sig (1-Stage & 2-Stage)",
        "1-Stage Harmony Only",
        "2-Stage Meta-IUT Only",
        "Background Pair",
    ]
    tier_col_c = [OKABE_BLUE, OKABE_BLUISH_GREEN, OKABE_ORANGE, "#CBD5E1"]

    zero_v = (
        alt.Chart(pl.DataFrame({"z0": [0.0]}))
        .mark_rule(strokeDash=[3, 3], stroke="#94A3B8", strokeWidth=0.8)
        .encode(x="z0:Q")
    )
    zero_h = (
        alt.Chart(pl.DataFrame({"z0": [0.0]}))
        .mark_rule(strokeDash=[3, 3], stroke="#94A3B8", strokeWidth=0.8)
        .encode(y="z0:Q")
    )
    pts_c = (
        alt.Chart(concordance_df)
        .mark_circle(size=48, opacity=0.82, stroke="#0F172A", strokeWidth=0.4)
        .encode(
            x=alt.X(
                "z_one_stage_harmony:Q",
                axis=create_nature_axis("1-Stage Harmony Pooled Signed Z (P6)"),
            ),
            y=alt.Y(
                "z_two_stage_meta:Q",
                axis=create_nature_axis("2-Stage Random-Effects Meta-IUT Signed Z (P3)"),
            ),
            color=alt.Color(
                "concordance_tier:N",
                sort=tier_dom_c,
                scale=alt.Scale(domain=tier_dom_c, range=tier_col_c),
                legend=alt.Legend(title="Discovery Concordance", orient="right", labelLimit=185),
            ),
        )
    )
    reg_c = (
        alt.Chart(concordance_df)
        .transform_regression("z_one_stage_harmony", "z_two_stage_meta")
        .mark_line(stroke="#0F172A", strokeWidth=1.3, strokeDash=[5, 3])
        .encode(x="z_one_stage_harmony:Q", y="z_two_stage_meta:Q")
    )
    panel_c = (
        alt.layer(zero_v, zero_h, pts_c, reg_c)
        .properties(
            width=370,
            height=275,
            title=alt.TitleParams(
                text="c  1-Stage Harmony Atlas vs. 2-Stage Meta-IUT Concordance",
                subtitle=conc_sub,
                anchor="start",
                fontSize=11.5,
                subtitleFontSize=9.5,
                subtitleColor="#475569",
            ),
        )
    )

    # Panel d: Leave-One-Cohort-Out (LOCO) stability forest plot for top 10 pan-cancer pairs
    if not loco_forest_df.is_empty():
        top_pairs_d = (
            loco_forest_df.lazy()
            .group_by("pair_label")
            .agg(
                pl.col("full_z_acat").first().alias("full_z_acat"),
                pl.col("min_loco_z_acat").min().alias("min_loco_z"),
                pl.col("loco_z_acat").max().alias("max_loco_z"),
            )
            .sort("full_z_acat", descending=True)
            .head(10)
            .collect()
        )
        keep_lbls = top_pairs_d["pair_label"].to_list()
        sub_loco = loco_forest_df.filter(pl.col("pair_label").is_in(keep_lbls))
    else:
        top_pairs_d = pl.DataFrame(
            {
                "pair_label": ["No Sig Pair"],
                "full_z_acat": [0.0],
                "min_loco_z": [0.0],
                "max_loco_z": [0.0],
            }
        )
        keep_lbls = ["No Sig Pair"]
        sub_loco = pl.DataFrame(
            {"pair_label": ["No Sig Pair"], "loco_z_acat": [0.0], "dropped_cohort": ["N/A"]}
        )

    z196_rule = (
        alt.Chart(pl.DataFrame({"z_thr": [1.96]}))
        .mark_rule(strokeDash=[4, 3], stroke=OKABE_VERMILION, strokeWidth=1.0)
        .encode(x="z_thr:Q")
    )
    span_d = (
        alt.Chart(top_pairs_d)
        .mark_rule(stroke="#475569", strokeWidth=1.5)
        .encode(
            y=alt.Y("pair_label:N", sort=keep_lbls, axis=create_nature_axis("")),
            x=alt.X(
                "min_loco_z:Q",
                scale=alt.Scale(zero=True),
                axis=create_nature_axis("Omnibus IUT |Z| Across Leave-One-Cohort-Out Drops"),
            ),
            x2="max_loco_z:Q",
        )
    )
    drops_d = (
        alt.Chart(sub_loco)
        .mark_circle(size=42, opacity=0.75, color=OKABE_BLUISH_GREEN, stroke="#0F172A", strokeWidth=0.4)
        .encode(y=alt.Y("pair_label:N", sort=keep_lbls), x="loco_z_acat:Q")
    )
    full_d = (
        alt.Chart(top_pairs_d)
        .mark_point(
            shape="diamond",
            filled=True,
            size=95,
            fill=OKABE_BLUE,
            stroke="#0F172A",
            strokeWidth=1.0,
        )
        .encode(y=alt.Y("pair_label:N", sort=keep_lbls), x="full_z_acat:Q")
    )
    panel_d = (
        alt.layer(z196_rule, span_d, drops_d, full_d)
        .properties(
            width=370,
            height=275,
            title=alt.TitleParams(
                text="d  Leave-One-Cohort-Out (LOCO) Pan-Cancer Stability Forest",
                subtitle="Blue diamond: All 9 cohorts | Green dots: 9 LOCO drops | Dashed: |Z| = 1.96",
                anchor="start",
                fontSize=11.5,
                subtitleFontSize=9.5,
                subtitleColor="#475569",
            ),
        )
    )

    top_ab = alt.hconcat(panel_a, panel_b, spacing=30).resolve_scale(color="independent")
    bot_cd = alt.hconcat(panel_c, panel_d, spacing=30).resolve_scale(color="independent")
    comp = (
        alt.vconcat(top_ab, bot_cd, spacing=28)
        .resolve_scale(color="independent")
        .properties(
            title=alt.TitleParams(
                text="Statistical Power Ablation, Null Calibration & Multi-Cohort Stability",
                subtitle=(
                    "Combining 9 Clinical ICB Datasets (N = 613 Biopsies) + EB Hurdle Presence + "
                    "Depth Residualization + B=64 DMPP + Omnibus ACAT + Pan-iAtlas Prior"
                ),
                anchor="start",
                fontSize=13.5,
                subtitleFontSize=10.5,
                subtitleColor="#334155",
                offset=12,
            )
        )
    )
    return cast(AltairComposite, apply_nature_methods_theme(comp))


_HEATMAP_COMPARTMENT_ORDER: tuple[TmeCompartment, ...] = (
    TmeCompartment.T_NK,
    TmeCompartment.B_LINEAGE,
    TmeCompartment.PLASMA,
    TmeCompartment.MYELOID_MAST,
    TmeCompartment.STROMAL_MESENCHYMAL,
    TmeCompartment.ENDOTHELIAL_VASCULAR,
    TmeCompartment.MALIGNANT_EPITHELIAL,
)

_HEATMAP_COMPARTMENT_DISPLAY: Mapping[TmeCompartment, tuple[str, str, str]] = {
    TmeCompartment.T_NK: ("T / NK", "T / NK", OKABE_BLUE),
    TmeCompartment.B_LINEAGE: ("B Lineage", "B", OKABE_SKY_BLUE),
    TmeCompartment.PLASMA: ("Plasma", "Pl", OKABE_BLUISH_GREEN),
    TmeCompartment.MYELOID_MAST: ("Myeloid / Mast", "Myeloid / Mast", OKABE_ORANGE),
    TmeCompartment.STROMAL_MESENCHYMAL: ("Stromal / CAF", "CAF", OKABE_VERMILION),
    TmeCompartment.ENDOTHELIAL_VASCULAR: ("Endothelial", "Endo", OKABE_REDDISH_PURPLE),
    TmeCompartment.MALIGNANT_EPITHELIAL: ("Malignant / Epith.", "Epith", OKABE_YELLOW),
}


def _hidden_canvas_axis() -> alt.Axis:
    """Returns a clean blank axis without using Python `None` literals."""
    return alt.Axis(ticks=False, labels=False, domain=False, grid=False, title="")


def format_stratified_perm_stars(perm_pvalue: float) -> str:
    """Formats exact within-stratum permutation p-values into standard significance stars."""
    match perm_pvalue:
        case p if p <= 0.001:
            return "***"
        case p if p < 0.01:
            return "**"
        case p if p < 0.05:
            return "*"
        case _:
            return ""


def _oriented_pair_z_acat(
    row_state: str,
    col_state: str,
    state_a: str,
    z_signed_ab: float,
    dominant_regime: str,
) -> float:
    """Computes directional `Z_ACAT(row_state, col_state)` preserving balance antisymmetry."""
    is_forward = row_state == state_a
    is_antisym = (
        dominant_regime == ContinuousAxis.LOG_RATIO_BALANCE.value
        or dominant_regime in ("quadrant:+-", "quadrant:-+")
    )
    return z_signed_ab if (is_forward or not is_antisym) else -z_signed_ab


def select_core_hub_state_codes(
    verified_pan_df: pl.DataFrame,
    fdr_threshold: float = 0.10,
    include_malignant_epithelial: bool = False,
    min_hub_degree: int = 2,
    anchor_perm_p_max: float = 0.01,
    anchor_abs_z_min: float = 3.25,
) -> frozenset[str]:
    """Selects core TME hub cell states via multi-partner verified degree or strong perm-verified Z."""
    if verified_pan_df.is_empty():
        return frozenset()

    allowed_supergroups = (
        {
            TmeSuperGroup.IMMUNE_CELLS,
            TmeSuperGroup.STROMAL_ENDOTHELIAL,
            TmeSuperGroup.MALIGNANT_EPITHELIAL,
        }
        if include_malignant_epithelial
        else {TmeSuperGroup.IMMUNE_CELLS, TmeSuperGroup.STROMAL_ENDOTHELIAL}
    )

    sig_deg: dict[str, int] = {}
    strong_anchors: set[str] = set()

    for r in verified_pan_df.to_dicts():
        sa = str(r["state_a"]).strip()
        sb = str(r["state_b"]).strip()
        sg_a = classify_compartment_supergroup(classify_tme_compartment(sa))
        sg_b = classify_compartment_supergroup(classify_tme_compartment(sb))
        both_allowed = False
        match (sg_a, sg_b):
            case (Some(ga), Some(gb)) if ga in allowed_supergroups and gb in allowed_supergroups:
                both_allowed = True
            case _:
                both_allowed = False
        if not both_allowed:
            continue

        fdr_v = float(cast(float | int | str, r["fdr_acat_iut"]))
        is_sig = bool(r["is_significant"]) and (fdr_v < fdr_threshold)
        if not is_sig:
            continue

        perm_p = float(cast(float | int | str, r["emp_perm_pvalue"]))
        z_abs = abs(float(cast(float | int | str, r["z_acat_signed"])))
        for sc in (sa, sb):
            sig_deg[sc] = sig_deg.get(sc, 0) + 1
            if perm_p <= anchor_perm_p_max and z_abs >= anchor_abs_z_min:
                strong_anchors.add(sc)

    return frozenset(
        sc for sc, deg in sig_deg.items() if deg >= min_hub_degree or sc in strong_anchors
    )


def cluster_states_within_compartments(
    verified_pan_df: pl.DataFrame,
    include_malignant_epithelial: bool = False,
    allowed_state_codes: Maybe[frozenset[str]] = Nothing,
    x_leaf: float = -7.70,
    dend_span: float = 1.45,
    leaf_stub: float = 0.32,
) -> tuple[tuple[HeatmapClusteredState, ...], list[dict[str, float]]]:
    """Clusters cell states hierarchically within compartments and groups into TME super-groups."""
    if verified_pan_df.is_empty():
        return (), []

    allowed_supergroups = (
        {
            TmeSuperGroup.IMMUNE_CELLS,
            TmeSuperGroup.STROMAL_ENDOTHELIAL,
            TmeSuperGroup.MALIGNANT_EPITHELIAL,
        }
        if include_malignant_epithelial
        else {TmeSuperGroup.IMMUNE_CELLS, TmeSuperGroup.STROMAL_ENDOTHELIAL}
    )

    state_readable: dict[str, str] = {}
    state_comp: dict[str, TmeCompartment] = {}
    state_super: dict[str, TmeSuperGroup] = {}

    for r in verified_pan_df.to_dicts():
        for code_k, read_k in (("state_a", "state_a_readable"), ("state_b", "state_b_readable")):
            sc = str(r[code_k]).strip()
            code_ok = True
            match allowed_state_codes:
                case Some(allow_set):
                    code_ok = sc in allow_set
                case _:
                    code_ok = True
            if not code_ok:
                continue
            cp = classify_tme_compartment(sc)
            sg_maybe = classify_compartment_supergroup(cp)
            match sg_maybe:
                case Some(sg) if sg in allowed_supergroups:
                    state_readable[sc] = str(r[read_k]).strip()
                    state_comp[sc] = cp
                    state_super[sc] = sg
                case _:
                    pass

    all_codes = sorted(state_comp.keys())
    if not all_codes:
        return (), []

    code_to_raw_i = {sc: idx for idx, sc in enumerate(all_codes)}
    n_tot = len(all_codes)
    z_prof = np.zeros((n_tot, n_tot), dtype=np.float64)

    peak_ver_z: dict[str, float] = {sc: 0.0 for sc in all_codes}
    peak_any_z: dict[str, float] = {sc: 0.0 for sc in all_codes}
    lead_loco_z: dict[str, float] = {sc: 0.0 for sc in all_codes}
    best_ver_p: dict[str, float] = {sc: 2.0 for sc in all_codes}
    has_ver: dict[str, bool] = {sc: False for sc in all_codes}

    for r in verified_pan_df.to_dicts():
        sa = str(r["state_a"]).strip()
        sb = str(r["state_b"]).strip()
        if sa not in code_to_raw_i or sb not in code_to_raw_i:
            continue
        ia, ib = code_to_raw_i[sa], code_to_raw_i[sb]
        z_ab = float(r["z_acat_signed"])
        z_abs = abs(z_ab)
        reg = str(r["dominant_regime"])
        z_prof[ia, ib] = _oriented_pair_z_acat(sa, sb, sa, z_ab, reg)
        z_prof[ib, ia] = _oriented_pair_z_acat(sb, sa, sa, z_ab, reg)

        p_ac = float(r["pvalue_acat_iut"])
        min_loco = float(r["min_loco_z_acat"])
        is_sig = bool(r["is_significant"])
        for sc in (sa, sb):
            peak_any_z[sc] = max(peak_any_z[sc], z_abs)
            if is_sig:
                has_ver[sc] = True
                peak_ver_z[sc] = max(peak_ver_z[sc], z_abs)
                if p_ac < best_ver_p[sc]:
                    best_ver_p[sc] = p_ac
                    lead_loco_z[sc] = min_loco

    ordered_codes: list[str] = []
    comp_Dend_raw: list[tuple[int, int, list[list[float]], list[list[float]]]] = []

    for comp_enum in _HEATMAP_COMPARTMENT_ORDER:
        c_codes = [sc for sc in all_codes if state_comp[sc] == comp_enum]
        if not c_codes:
            continue
        start_pos = len(ordered_codes)
        m_c = len(c_codes)
        if m_c == 1:
            ordered_codes.extend(c_codes)
            comp_Dend_raw.append((start_pos, 1, [], []))
        elif m_c == 2:
            c_sorted = sorted(
                c_codes, key=lambda s: (-peak_ver_z[s], -peak_any_z[s], s)
            )
            ordered_codes.extend(c_sorted)
            comp_Dend_raw.append(
                (start_pos, 2, [[5.0, 5.0, 15.0, 15.0]], [[0.0, 1.0, 1.0, 0.0]])
            )
        else:
            raw_idxs = np.array([code_to_raw_i[sc] for sc in c_codes], dtype=np.int64)
            sub_mat = z_prof[raw_idxs, :]
            d_vec = np.asarray(pdist(sub_mat, metric="euclidean"), dtype=np.float64)
            if float(np.max(d_vec)) <= 1e-9:
                ramp = np.arange(m_c, dtype=np.float64).reshape(-1, 1) * 1e-3
                d_vec = np.asarray(pdist(ramp, metric="euclidean"), dtype=np.float64)
            else:
                d_vec = np.maximum(d_vec, 1e-6)
            z_link = linkage(d_vec, method="average", optimal_ordering=True)
            leaf_ord = [int(i) for i in np.asarray(leaves_list(z_link), dtype=np.int64).tolist()]
            dend_info = dendrogram(z_link, no_plot=True)
            icoords = [[float(v) for v in seg] for seg in dend_info["icoord"]]
            dcoords = [[float(v) for v in seg] for seg in dend_info["dcoord"]]
            ordered_codes.extend([c_codes[i] for i in leaf_ord])
            comp_Dend_raw.append((start_pos, m_c, icoords, dcoords))

    n_final = len(ordered_codes)
    dend_segments: list[dict[str, float]] = []

    for r_start, m_c, icoords, dcoords in comp_Dend_raw:
        if m_c == 1:
            y_leaf = float(n_final) - 0.5 - float(r_start)
            dend_segments.append(
                {"x0": x_leaf - leaf_stub, "y0": y_leaf, "x1": x_leaf, "y1": y_leaf}
            )
            continue
        max_h = max((max(dc) for dc in dcoords), default=1.0)
        safe_h = max(max_h, 1e-6)
        for ic, dc in zip(icoords, dcoords):
            ys = [float(n_final) - 0.5 - float(r_start) - (iv - 5.0) / 10.0 for iv in ic]
            xs = [x_leaf - (dv / safe_h) * dend_span for dv in dc]
            for seg_i in range(3):
                dend_segments.append(
                    {
                        "x0": xs[seg_i],
                        "y0": ys[seg_i],
                        "x1": xs[seg_i + 1],
                        "y1": ys[seg_i + 1],
                    }
                )

    entries = tuple(
        HeatmapClusteredState(
            order_idx=idx,
            state_code=sc,
            readable_label=state_readable[sc],
            compartment=state_comp[sc],
            compartment_display=_HEATMAP_COMPARTMENT_DISPLAY[state_comp[sc]][0],
            compartment_short=_HEATMAP_COMPARTMENT_DISPLAY[state_comp[sc]][1],
            compartment_color=_HEATMAP_COMPARTMENT_DISPLAY[state_comp[sc]][2],
            supergroup=state_super[sc],
            supergroup_display=supergroup_display_label(state_super[sc]),
            peak_verified_abs_z=peak_ver_z[sc],
            peak_any_abs_z=peak_any_z[sc],
            lead_min_loco_z=lead_loco_z[sc],
            has_verified_pair=has_ver[sc],
        )
        for idx, sc in enumerate(ordered_codes)
    )
    return entries, dend_segments


def build_pan_cohort_z_acat_heatmap_figure(
    verified_pan_df: pl.DataFrame,
    fdr_threshold: float = 0.10,
    include_malignant_epithelial: bool = False,
) -> AltairComposite:
    """Builds the Nature Methods clustered cross-compartment `Z_ACAT` heatmap + `LOCO` strip."""
    clustered_states, dend_segments = cluster_states_within_compartments(
        verified_pan_df=verified_pan_df,
        include_malignant_epithelial=include_malignant_epithelial,
    )
    n_st = len(clustered_states)
    if n_st == 0:
        return cast(
            AltairComposite,
            alt.Chart(pl.DataFrame({"x": [0.0], "y": [0.0]})).mark_point(opacity=0.0),
        )

    state_to_idx = {st.state_code: st.order_idx for st in clustered_states}
    pair_lookup: dict[tuple[str, str], dict[str, object]] = {}
    for row in verified_pan_df.to_dicts():
        sa, sb = str(row["state_a"]).strip(), str(row["state_b"]).strip()
        if sa in state_to_idx and sb in state_to_idx:
            pair_lookup[(min(sa, sb), max(sa, sb))] = row

    z_cap = 4.5
    cross_rows: list[dict[str, object]] = []
    intra_rows: list[dict[str, object]] = []
    missing_rows: list[dict[str, object]] = []
    sig_frame_rows: list[dict[str, object]] = []
    star_dark_rows: list[dict[str, object]] = []
    star_light_rows: list[dict[str, object]] = []

    for r_st in clustered_states:
        r_i = r_st.order_idx
        y0 = float(n_st - 1 - r_i)
        y1 = float(n_st - r_i)
        ym = 0.5 * (y0 + y1)
        for c_st in clustered_states:
            c_i = c_st.order_idx
            x0 = float(c_i)
            x1 = float(c_i + 1)
            xm = 0.5 * (x0 + x1)

            if r_st.compartment == c_st.compartment:
                intra_rows.append({"x0": x0, "x1": x1, "y0": y0, "y1": y1})
                continue

            key_ab = (min(r_st.state_code, c_st.state_code), max(r_st.state_code, c_st.state_code))
            if key_ab not in pair_lookup:
                missing_rows.append({"x0": x0, "x1": x1, "y0": y0, "y1": y1, "xm": xm, "ym": ym})
                continue

            pr = pair_lookup[key_ab]
            z_oriented = _oriented_pair_z_acat(
                row_state=r_st.state_code,
                col_state=c_st.state_code,
                state_a=str(pr["state_a"]).strip(),
                z_signed_ab=float(cast(float | int | str, pr["z_acat_signed"])),
                dominant_regime=str(pr["dominant_regime"]),
            )
            z_clamped = float(np.clip(z_oriented, -z_cap, z_cap))
            perm_p = float(cast(float | int | str, pr["emp_perm_pvalue"]))
            p_acat_v = float(cast(float | int | str, pr["pvalue_acat_iut"]))
            fdr_v = float(cast(float | int | str, pr["fdr_acat_iut"]))
            is_sig = bool(pr["is_significant"]) and (fdr_v < fdr_threshold)
            stars = format_stratified_perm_stars(max(perm_p, p_acat_v))

            cross_rows.append(
                {
                    "x0": x0,
                    "x1": x1,
                    "y0": y0,
                    "y1": y1,
                    "z_acat_clamped": z_clamped,
                    "z_acat_raw": z_oriented,
                }
            )
            if is_sig:
                sig_frame_rows.append(
                    {
                        "x0": x0 + 0.045,
                        "x1": x1 - 0.045,
                        "y0": y0 + 0.045,
                        "y1": y1 - 0.045,
                    }
                )
            if stars:
                star_rec = {"xm": xm, "ym": ym, "stars": stars}
                if abs(z_oriented) >= 3.20:
                    star_light_rows.append(star_rec)
                else:
                    star_dark_rows.append(star_rec)

    x_dom = [-14.0, 35.1]
    y_dom = [-8.35, float(n_st) + 3.45]
    x_sc = alt.Scale(domain=x_dom, nice=False, zero=False)
    y_sc = alt.Scale(domain=y_dom, nice=False, zero=False)
    blank_ax = _hidden_canvas_axis()

    def _xy_enc(
        x_f: str, y_f: str, x2_f: str = "", y2_f: str = ""
    ) -> dict[str, object]:
        enc: dict[str, object] = {
            "x": alt.X(f"{x_f}:Q", scale=x_sc, axis=blank_ax),
            "y": alt.Y(f"{y_f}:Q", scale=y_sc, axis=blank_ax),
        }
        if x2_f:
            enc["x2"] = alt.X2(f"{x2_f}:Q")
        if y2_f:
            enc["y2"] = alt.Y2(f"{y2_f}:Q")
        return enc

    layers: list[alt.Chart] = []

    if intra_rows:
        layers.append(
            alt.Chart(pl.DataFrame(intra_rows))
            .mark_rect(fill="#F1F5F9", stroke="#E2E8F0", strokeWidth=0.45)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )

    if missing_rows:
        miss_df = pl.DataFrame(missing_rows)
        layers.append(
            alt.Chart(miss_df)
            .mark_rect(fill="#F8FAFC", stroke="#E2E8F0", strokeWidth=0.45)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )

    if cross_rows:
        layers.append(
            alt.Chart(pl.DataFrame(cross_rows))
            .mark_rect(stroke="#CBD5E1", strokeWidth=0.45)
            .encode(
                **_xy_enc("x0", "y0", "x1", "y1"),
                color=alt.Color(
                    "z_acat_clamped:Q",
                    scale=alt.Scale(
                        domain=[-z_cap, 0.0, z_cap],
                        range=[OKABE_VERMILION, "#FFFFFF", OKABE_BLUE],
                        interpolate="rgb",
                    ),
                    legend=alt.Legend(
                        title=["Signed Z_ACAT", "(Row vs. Col)"],
                        orient="right",
                        gradientLength=230,
                        gradientThickness=13,
                        values=[-4.0, -2.0, 0.0, 2.0, 4.0],
                        titleFontSize=9.8,
                        labelFontSize=8.8,
                    ),
                ),
            )
        )

    # Diagonal self-identity guide across intra-compartment blocks
    diag_df = pl.DataFrame(
        [{"x0": 0.0, "y0": float(n_st), "x1": float(n_st), "y1": 0.0}]
    )
    layers.append(
        alt.Chart(diag_df)
        .mark_rule(color="#CBD5E1", strokeWidth=0.8, strokeDash=[3, 3])
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )

    if sig_frame_rows:
        layers.append(
            alt.Chart(pl.DataFrame(sig_frame_rows))
            .mark_rect(fillOpacity=0.0, stroke="#0F172A", strokeWidth=1.15)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )

    if star_dark_rows:
        layers.append(
            alt.Chart(pl.DataFrame(star_dark_rows))
            .mark_text(
                fontSize=9.2,
                fontWeight="bold",
                fill="#0F172A",
                align="center",
                baseline="middle",
                dy=1.5,
            )
            .encode(**_xy_enc("xm", "ym"), text="stars:N")
        )

    if star_light_rows:
        layers.append(
            alt.Chart(pl.DataFrame(star_light_rows))
            .mark_text(
                fontSize=9.2,
                fontWeight="bold",
                fill="#FFFFFF",
                align="center",
                baseline="middle",
                dy=1.5,
            )
            .encode(**_xy_enc("xm", "ym"), text="stars:N")
        )

    # Compartment contiguous spans & structural dividers
    comp_spans: list[tuple[TmeCompartment, int, int]] = []
    for st in clustered_states:
        if not comp_spans or comp_spans[-1][0] != st.compartment:
            comp_spans.append((st.compartment, st.order_idx, st.order_idx + 1))
        else:
            prev_c, s_i, _ = comp_spans[-1]
            comp_spans[-1] = (prev_c, s_i, st.order_idx + 1)

    super_spans: list[tuple[TmeSuperGroup, int, int]] = []
    for st in clustered_states:
        if not super_spans or super_spans[-1][0] != st.supergroup:
            super_spans.append((st.supergroup, st.order_idx, st.order_idx + 1))
        else:
            prev_sg, s_i, _ = super_spans[-1]
            super_spans[-1] = (prev_sg, s_i, st.order_idx + 1)

    sub_div_lines: list[dict[str, float]] = []
    for _, _, end_i in comp_spans[:-1]:
        sub_div_lines.append({"x0": float(end_i), "y0": 0.0, "x1": float(end_i), "y1": float(n_st)})
        y_div = float(n_st - end_i)
        sub_div_lines.append({"x0": 0.0, "y0": y_div, "x1": float(n_st), "y1": y_div})

    if sub_div_lines:
        layers.append(
            alt.Chart(pl.DataFrame(sub_div_lines))
            .mark_rule(color="#64748B", strokeWidth=0.75, strokeDash=[3, 2])
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )

    super_div_lines: list[dict[str, float]] = []
    for _, _, end_sg in super_spans[:-1]:
        super_div_lines.append(
            {"x0": float(end_sg), "y0": 0.0, "x1": float(end_sg), "y1": float(n_st)}
        )
        y_sg = float(n_st - end_sg)
        super_div_lines.append({"x0": 0.0, "y0": y_sg, "x1": float(n_st), "y1": y_sg})

    if super_div_lines:
        layers.append(
            alt.Chart(pl.DataFrame(super_div_lines))
            .mark_rule(color="#0F172A", strokeWidth=1.85)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )

    # Outer crisp wireframe around NxN matrix
    layers.append(
        alt.Chart(pl.DataFrame([{"x0": 0.0, "y0": 0.0, "x1": float(n_st), "y1": float(n_st)}]))
        .mark_rect(fillOpacity=0.0, stroke="#0F172A", strokeWidth=1.1)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )

    # Top & Left Compartment Color Bars + Top Labels
    top_comp_labels: list[dict[str, object]] = []
    for comp_e, s_i, e_i in comp_spans:
        c_full, c_short, c_hex = _HEATMAP_COMPARTMENT_DISPLAY[comp_e]
        span_w = e_i - s_i
        bar_df = pl.DataFrame(
            [
                {
                    "x0": float(s_i) + 0.04,
                    "x1": float(e_i) - 0.04,
                    "y0": float(n_st) + 0.25,
                    "y1": float(n_st) + 0.72,
                },
                {
                    "x0": -7.55,
                    "x1": -7.15,
                    "y0": float(n_st - e_i) + 0.04,
                    "y1": float(n_st - s_i) - 0.04,
                },
            ]
        )
        layers.append(
            alt.Chart(bar_df)
            .mark_rect(fill=c_hex, stroke="#0F172A", strokeWidth=0.55)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )
        x_lbl = 0.5 * float(s_i + e_i)
        if comp_e == TmeCompartment.STROMAL_MESENCHYMAL and span_w == 1:
            x_lbl -= 0.10
        elif comp_e == TmeCompartment.ENDOTHELIAL_VASCULAR and span_w == 1:
            x_lbl += 0.12
        lbl_txt = f"{c_full} (n={span_w})" if span_w >= 4 else c_short
        top_comp_labels.append({"x": x_lbl, "y": float(n_st) + 1.08, "txt": lbl_txt})

    if top_comp_labels:
        layers.append(
            alt.Chart(pl.DataFrame(top_comp_labels))
            .mark_text(
                fontSize=8.4,
                fontWeight="bold",
                fill="#1E293B",
                align="center",
                baseline="middle",
            )
            .encode(**_xy_enc("x", "y"), text="txt:N")
        )

    # Within-compartment hierarchical dendrogram on the left margin
    if dend_segments:
        layers.append(
            alt.Chart(pl.DataFrame(dend_segments))
            .mark_rule(color="#334155", strokeWidth=0.85)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )

    # Outer Super-Group Brackets (`Immune Cells` vs. `Stromal / Endothelial`) on Top & Left
    bracket_segs: list[dict[str, float]] = []
    top_sg_labels: list[dict[str, object]] = []
    left_sg_labels: list[dict[str, object]] = []

    y_br_top = float(n_st) + 1.88
    x_br_left = -9.60
    for sg_e, s_i, e_i in super_spans:
        x_a = float(s_i) + 0.08
        x_b = float(e_i) - 0.08
        bracket_segs.extend(
            [
                {"x0": x_a, "y0": y_br_top, "x1": x_b, "y1": y_br_top},
                {"x0": x_a, "y0": y_br_top - 0.25, "x1": x_a, "y1": y_br_top},
                {"x0": x_b, "y0": y_br_top - 0.25, "x1": x_b, "y1": y_br_top},
            ]
        )
        y_a = float(n_st - e_i) + 0.08
        y_b = float(n_st - s_i) - 0.08
        bracket_segs.extend(
            [
                {"x0": x_br_left, "y0": y_a, "x1": x_br_left, "y1": y_b},
                {"x0": x_br_left, "y0": y_a, "x1": x_br_left + 0.25, "y1": y_a},
                {"x0": x_br_left, "y0": y_b, "x1": x_br_left + 0.25, "y1": y_b},
            ]
        )
        span_n = e_i - s_i
        sg_disp = supergroup_display_label(sg_e)
        top_txt = (
            f"{sg_disp} ({span_n} States)"
            if span_n >= 6
            else ("Stromal / Endo." if sg_e == TmeSuperGroup.STROMAL_ENDOTHELIAL else sg_disp)
        )
        top_sg_labels.append(
            {"x": 0.5 * float(s_i + e_i), "y": y_br_top + 0.48, "txt": top_txt}
        )

        y_mid_sg = 0.5 * (y_a + y_b)
        if sg_e == TmeSuperGroup.IMMUNE_CELLS:
            left_sg_labels.append({"x": x_br_left - 0.30, "y": y_mid_sg + 0.36, "txt": "Immune"})
            left_sg_labels.append({"x": x_br_left - 0.30, "y": y_mid_sg - 0.36, "txt": "Cells"})
        elif sg_e == TmeSuperGroup.STROMAL_ENDOTHELIAL:
            left_sg_labels.append({"x": x_br_left - 0.30, "y": y_mid_sg + 0.34, "txt": "Stromal /"})
            left_sg_labels.append(
                {"x": x_br_left - 0.30, "y": y_mid_sg - 0.34, "txt": "Endothelial"}
            )
        else:
            left_sg_labels.append({"x": x_br_left - 0.30, "y": y_mid_sg, "txt": sg_disp})

    if bracket_segs:
        layers.append(
            alt.Chart(pl.DataFrame(bracket_segs))
            .mark_rule(color="#0F172A", strokeWidth=1.35)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )
    if top_sg_labels:
        layers.append(
            alt.Chart(pl.DataFrame(top_sg_labels))
            .mark_text(
                fontSize=9.8,
                fontWeight="bold",
                fill="#0F172A",
                align="center",
                baseline="middle",
            )
            .encode(**_xy_enc("x", "y"), text="txt:N")
        )
    if left_sg_labels:
        layers.append(
            alt.Chart(pl.DataFrame(left_sg_labels))
            .mark_text(
                fontSize=9.4,
                fontWeight="bold",
                fill="#0F172A",
                align="right",
                baseline="middle",
            )
            .encode(**_xy_enc("x", "y"), text="txt:N")
        )

    # Outward ticks + Y-axis & rotated X-axis cell-state labels
    axis_ticks: list[dict[str, float]] = []
    y_lbl_bold: list[dict[str, object]] = []
    y_lbl_reg: list[dict[str, object]] = []
    x_lbl_bold: list[dict[str, object]] = []
    x_lbl_reg: list[dict[str, object]] = []

    for st in clustered_states:
        ym = float(n_st) - 0.5 - float(st.order_idx)
        xm = float(st.order_idx) + 0.5
        axis_ticks.append({"x0": -0.22, "y0": ym, "x1": 0.0, "y1": ym})
        axis_ticks.append({"x0": xm, "y0": -0.22, "x1": xm, "y1": 0.0})
        y_rec = {"x": -0.32, "y": ym, "txt": st.readable_label}
        x_rec = {"x": xm, "y": -0.34, "txt": st.readable_label}
        if st.has_verified_pair:
            y_lbl_bold.append(y_rec)
            x_lbl_bold.append(x_rec)
        else:
            y_lbl_reg.append(y_rec)
            x_lbl_reg.append(x_rec)

    layers.append(
        alt.Chart(pl.DataFrame(axis_ticks))
        .mark_rule(color="#1E293B", strokeWidth=0.75)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )
    if y_lbl_reg:
        layers.append(
            alt.Chart(pl.DataFrame(y_lbl_reg))
            .mark_text(fontSize=8.7, fill="#475569", align="right", baseline="middle")
            .encode(**_xy_enc("x", "y"), text="txt:N")
        )
    if y_lbl_bold:
        layers.append(
            alt.Chart(pl.DataFrame(y_lbl_bold))
            .mark_text(
                fontSize=8.8, fontWeight="bold", fill="#0F172A", align="right", baseline="middle"
            )
            .encode(**_xy_enc("x", "y"), text="txt:N")
        )
    if x_lbl_reg:
        layers.append(
            alt.Chart(pl.DataFrame(x_lbl_reg))
            .mark_text(
                fontSize=8.6, fill="#475569", align="right", baseline="middle", angle=270
            )
            .encode(**_xy_enc("x", "y"), text="txt:N")
        )
    if x_lbl_bold:
        layers.append(
            alt.Chart(pl.DataFrame(x_lbl_bold))
            .mark_text(
                fontSize=8.7,
                fontWeight="bold",
                fill="#0F172A",
                align="right",
                baseline="middle",
                angle=270,
            )
            .encode(**_xy_enc("x", "y"), text="txt:N")
        )

    # Panel headers & bottom axis title
    panel_hdr_df = pl.DataFrame(
        [
            {
                "x": -13.7,
                "y": float(n_st) + 3.05,
                "txt": "a  Clustered Cross-Compartment Z_ACAT Matrix (Immune vs. Stromal / Endothelial)",
            },
            {
                "x": 27.1,
                "y": float(n_st) + 3.05,
                "txt": "b  9-Cohort LOCO Hub Stability",
            },
        ]
    )
    layers.append(
        alt.Chart(panel_hdr_df)
        .mark_text(fontSize=10.8, fontWeight="bold", fill="#0F172A", align="left", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )
    layers.append(
        alt.Chart(
            pl.DataFrame(
                [
                    {
                        "x": 27.1,
                        "y": float(n_st) + 2.38,
                        "txt": "Bar: Peak |Z_ACAT| · ◆: Worst min_LOCO_Z",
                    }
                ]
            )
        )
        .mark_text(fontSize=8.4, fill="#475569", align="left", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )
    layers.append(
        alt.Chart(
            pl.DataFrame(
                [
                    {
                        "x": 0.5 * float(n_st),
                        "y": -7.65,
                        "txt": "Column Cell State j  —  Reference State in Compositional Balance log₂(Row i / Col j)",
                    }
                ]
            )
        )
        .mark_text(fontSize=9.4, fontWeight="bold", fill="#1E293B", align="center", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )

    # Synchronized Right-Side Panel `b`: Peak |Z_ACAT| & Worst-Case `min_LOCO_Z` (9 Drops)
    loco_x0 = 27.15
    loco_scale = 1.32  # 5.0 Z -> 6.60 grid units (x = 33.75)
    loco_x_max = loco_x0 + 5.0 * loco_scale

    layers.append(
        alt.Chart(
            pl.DataFrame([{"x0": loco_x0, "y0": 0.0, "x1": loco_x_max, "y1": float(n_st)}])
        )
        .mark_rect(fill="#FFFFFF", stroke="#CBD5E1", strokeWidth=0.65)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )

    loco_grid = [
        {
            "x0": loco_x0 + zv * loco_scale,
            "y0": 0.0,
            "x1": loco_x0 + zv * loco_scale,
            "y1": float(n_st),
        }
        for zv in (1.0, 2.0, 3.0, 4.0)
    ]
    layers.append(
        alt.Chart(pl.DataFrame(loco_grid))
        .mark_rule(color="#F1F5F9", strokeWidth=0.7)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )

    x_196 = loco_x0 + 1.96 * loco_scale
    layers.append(
        alt.Chart(pl.DataFrame([{"x0": x_196, "y0": 0.0, "x1": x_196, "y1": float(n_st)}]))
        .mark_rule(color=OKABE_VERMILION, strokeWidth=1.0, strokeDash=[4, 3])
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )

    ver_bars: list[dict[str, float]] = []
    bg_bars: list[dict[str, float]] = []
    loco_pass_pts: list[dict[str, float]] = []
    loco_warn_pts: list[dict[str, float]] = []

    for st in clustered_states:
        ym = float(n_st) - 0.5 - float(st.order_idx)
        if st.has_verified_pair:
            z_b = min(st.peak_verified_abs_z, 5.0)
            ver_bars.append(
                {
                    "x0": loco_x0,
                    "x1": loco_x0 + z_b * loco_scale,
                    "y0": ym - 0.27,
                    "y1": ym + 0.27,
                }
            )
            loco_clamped = float(np.clip(st.lead_min_loco_z, 0.0, 5.0))
            pt_rec = {"x": loco_x0 + loco_clamped * loco_scale, "y": ym}
            if st.lead_min_loco_z >= 1.96:
                loco_pass_pts.append(pt_rec)
            else:
                loco_warn_pts.append(pt_rec)
        else:
            z_bg = min(st.peak_any_abs_z, 5.0)
            bg_bars.append(
                {
                    "x0": loco_x0,
                    "x1": loco_x0 + z_bg * loco_scale,
                    "y0": ym - 0.22,
                    "y1": ym + 0.22,
                }
            )

    if bg_bars:
        layers.append(
            alt.Chart(pl.DataFrame(bg_bars))
            .mark_rect(fill="#E2E8F0")
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )
    if ver_bars:
        layers.append(
            alt.Chart(pl.DataFrame(ver_bars))
            .mark_rect(fill="#94A3B8", stroke="#475569", strokeWidth=0.45)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )
    if loco_warn_pts:
        layers.append(
            alt.Chart(pl.DataFrame(loco_warn_pts))
            .mark_point(
                shape="diamond",
                filled=True,
                size=62,
                fill=OKABE_ORANGE,
                stroke="#0F172A",
                strokeWidth=0.8,
            )
            .encode(**_xy_enc("x", "y"))
        )
    if loco_pass_pts:
        layers.append(
            alt.Chart(pl.DataFrame(loco_pass_pts))
            .mark_point(
                shape="diamond",
                filled=True,
                size=68,
                fill=OKABE_BLUISH_GREEN,
                stroke="#0F172A",
                strokeWidth=0.85,
            )
            .encode(**_xy_enc("x", "y"))
        )

    # Explicit outward coordinate axis for LOCO Z panel
    loco_ax_segs = [{"x0": loco_x0, "y0": 0.0, "x1": loco_x_max, "y1": 0.0}]
    loco_ax_lbls: list[dict[str, object]] = []
    for zv in range(6):
        xz = loco_x0 + float(zv) * loco_scale
        loco_ax_segs.append({"x0": xz, "y0": -0.22, "x1": xz, "y1": 0.0})
        loco_ax_lbls.append({"x": xz, "y": -0.56, "txt": str(zv)})

    layers.append(
        alt.Chart(pl.DataFrame(loco_ax_segs))
        .mark_rule(color="#1E293B", strokeWidth=0.75)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )
    layers.append(
        alt.Chart(pl.DataFrame(loco_ax_lbls))
        .mark_text(fontSize=8.3, fill="#1E293B", align="center", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )
    layers.append(
        alt.Chart(
            pl.DataFrame([{"x": 0.5 * (loco_x0 + loco_x_max), "y": -1.12, "txt": "|Z| Score"}])
        )
        .mark_text(fontSize=8.6, fontWeight="bold", fill="#1E293B", align="center", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )

    # Bottom-Right Wireframe Key Box (Significance + Compartment Color Swatches)
    kx0, kx1 = 27.0, 34.95
    ky0, ky1 = -7.90, -1.72
    layers.append(
        alt.Chart(pl.DataFrame([{"x0": kx0, "y0": ky0, "x1": kx1, "y1": ky1}]))
        .mark_rect(fill="#FFFFFF", stroke="#CBD5E1", strokeWidth=0.7)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )
    key_txt_rows = [
        {"x": kx0 + 0.25, "y": -2.18, "txt": "Stratified Perm. & Lineage Key", "bold": True},
        {"x": kx0 + 0.25, "y": -2.82, "txt": "*** p ≤ 0.001 · ** p < 0.01 · * p < 0.05", "bold": False},
        {"x": kx0 + 0.85, "y": -3.45, "txt": f"ACAT FDR < {fdr_threshold:.2f} (verified)", "bold": False},
        {"x": kx0 + 0.85, "y": -4.08, "txt": "Intra-compartment (gated)", "bold": False},
    ]
    layers.append(
        alt.Chart(pl.DataFrame([r for r in key_txt_rows if bool(r["bold"])]))
        .mark_text(fontSize=8.7, fontWeight="bold", fill="#0F172A", align="left", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )
    layers.append(
        alt.Chart(pl.DataFrame([r for r in key_txt_rows if not bool(r["bold"])]))
        .mark_text(fontSize=8.1, fill="#334155", align="left", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )
    layers.append(
        alt.Chart(
            pl.DataFrame(
                [{"x0": kx0 + 0.26, "x1": kx0 + 0.68, "y0": -3.65, "y1": -3.25}]
            )
        )
        .mark_rect(fill="#FFFFFF", stroke="#0F172A", strokeWidth=1.15)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )
    layers.append(
        alt.Chart(
            pl.DataFrame(
                [{"x0": kx0 + 0.26, "x1": kx0 + 0.68, "y0": -4.28, "y1": -3.88}]
            )
        )
        .mark_rect(fill="#F1F5F9", stroke="#CBD5E1", strokeWidth=0.6)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )

    swatch_specs = [
        (TmeCompartment.T_NK, kx0 + 0.26, -4.85, "T / NK"),
        (TmeCompartment.B_LINEAGE, kx0 + 0.26, -5.50, "B Lineage"),
        (TmeCompartment.PLASMA, kx0 + 0.26, -6.15, "Plasma"),
        (TmeCompartment.MYELOID_MAST, kx0 + 4.15, -4.85, "Myeloid / Mast"),
        (TmeCompartment.STROMAL_MESENCHYMAL, kx0 + 4.15, -5.50, "Stromal / CAF"),
        (TmeCompartment.ENDOTHELIAL_VASCULAR, kx0 + 4.15, -6.15, "Endothelial"),
    ]
    sw_lbls: list[dict[str, object]] = []
    for cp_s, sx, sy, s_name in swatch_specs:
        c_hex = _HEATMAP_COMPARTMENT_DISPLAY[cp_s][2]
        layers.append(
            alt.Chart(
                pl.DataFrame([{"x0": sx, "x1": sx + 0.42, "y0": sy - 0.20, "y1": sy + 0.20}])
            )
            .mark_rect(fill=c_hex, stroke="#0F172A", strokeWidth=0.55)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )
        sw_lbls.append({"x": sx + 0.56, "y": sy, "txt": s_name})
    layers.append(
        alt.Chart(pl.DataFrame(sw_lbls))
        .mark_text(fontSize=8.0, fill="#1E293B", align="left", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )
    layers.append(
        alt.Chart(
            pl.DataFrame(
                [
                    {
                        "x": kx0 + 0.26,
                        "y": -6.82,
                        "txt": "◆ Green: min_LOCO_Z ≥ 1.96 · ◆ Orange: < 1.96",
                    }
                ]
            )
        )
        .mark_text(fontSize=7.8, fill="#475569", align="left", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )

    n_sig_host = len(
        [
            k_ab
            for k_ab, pr in pair_lookup.items()
            if bool(pr["is_significant"])
            and float(cast(float | int | str, pr["fdr_acat_iut"])) < fdr_threshold
        ]
    )
    combined = (
        alt.layer(*layers)
        .properties(
            width=1140,
            height=880,
            title=alt.TitleParams(
                text="Pan-Cancer Cross-Compartment Cell-State Interaction Heatmap (Omnibus Z_ACAT & Stratified Permutations)",
                subtitle=(
                    f"Unified 9-Cohort Harmony Atlas (N = 613 Biopsies, 20 CMH Strata) | "
                    f"{n_st} Host TME Cell States Clustered by Compartment (Immune vs. Stromal / Endothelial) | "
                    f"{n_sig_host} Verified Host TME Pairs (ACAT FDR < {fdr_threshold:.2f})"
                ),
                anchor="start",
                fontSize=13.2,
                subtitleFontSize=10.0,
                subtitleColor="#334155",
                offset=10,
            ),
        )
        .configure_view(strokeOpacity=0.0)
    )
    return cast(AltairComposite, apply_nature_methods_theme(combined))


def build_pan_cohort_core_hubs_z_acat_heatmap_figure(
    verified_pan_df: pl.DataFrame,
    fdr_threshold: float = 0.10,
    include_malignant_epithelial: bool = False,
) -> AltairComposite:
    """Builds the simplified Core Hub `Z_ACAT` heatmap displaying only the most important states."""
    core_codes = select_core_hub_state_codes(
        verified_pan_df=verified_pan_df,
        fdr_threshold=fdr_threshold,
        include_malignant_epithelial=include_malignant_epithelial,
    )
    clustered_states, dend_segments = cluster_states_within_compartments(
        verified_pan_df=verified_pan_df,
        include_malignant_epithelial=include_malignant_epithelial,
        allowed_state_codes=Some(core_codes),
        x_leaf=-4.70,
        dend_span=0.90,
        leaf_stub=0.22,
    )
    n_st = len(clustered_states)
    if n_st == 0:
        return cast(
            AltairComposite,
            alt.Chart(pl.DataFrame({"x": [0.0], "y": [0.0]})).mark_point(opacity=0.0),
        )

    state_to_idx = {st.state_code: st.order_idx for st in clustered_states}
    pair_lookup: dict[tuple[str, str], dict[str, object]] = {}
    for row in verified_pan_df.to_dicts():
        sa, sb = str(row["state_a"]).strip(), str(row["state_b"]).strip()
        if sa in state_to_idx and sb in state_to_idx:
            pair_lookup[(min(sa, sb), max(sa, sb))] = row

    z_cap = 4.5
    cross_rows: list[dict[str, object]] = []
    intra_rows: list[dict[str, object]] = []
    missing_rows: list[dict[str, object]] = []
    sig_frame_rows: list[dict[str, float]] = []
    z_val_dark_rows: list[dict[str, object]] = []
    z_val_light_rows: list[dict[str, object]] = []
    star_dark_rows: list[dict[str, object]] = []
    star_light_rows: list[dict[str, object]] = []

    for r_st in clustered_states:
        r_i = r_st.order_idx
        y0 = float(n_st - 1 - r_i)
        y1 = float(n_st - r_i)
        ym = 0.5 * (y0 + y1)
        for c_st in clustered_states:
            c_i = c_st.order_idx
            x0 = float(c_i)
            x1 = float(c_i + 1)
            xm = 0.5 * (x0 + x1)

            if r_st.compartment == c_st.compartment:
                intra_rows.append({"x0": x0, "x1": x1, "y0": y0, "y1": y1})
                continue

            key_ab = (min(r_st.state_code, c_st.state_code), max(r_st.state_code, c_st.state_code))
            if key_ab not in pair_lookup:
                missing_rows.append({"x0": x0, "x1": x1, "y0": y0, "y1": y1})
                continue

            pr = pair_lookup[key_ab]
            z_oriented = _oriented_pair_z_acat(
                row_state=r_st.state_code,
                col_state=c_st.state_code,
                state_a=str(pr["state_a"]).strip(),
                z_signed_ab=float(cast(float | int | str, pr["z_acat_signed"])),
                dominant_regime=str(pr["dominant_regime"]),
            )
            z_clamped = float(np.clip(z_oriented, -z_cap, z_cap))
            perm_p = float(cast(float | int | str, pr["emp_perm_pvalue"]))
            p_acat_v = float(cast(float | int | str, pr["pvalue_acat_iut"]))
            fdr_v = float(cast(float | int | str, pr["fdr_acat_iut"]))
            is_sig = bool(pr["is_significant"]) and (fdr_v < fdr_threshold)
            stars = format_stratified_perm_stars(max(perm_p, p_acat_v))

            cross_rows.append(
                {
                    "x0": x0,
                    "x1": x1,
                    "y0": y0,
                    "y1": y1,
                    "z_acat_clamped": z_clamped,
                    "z_acat_raw": z_oriented,
                }
            )
            if is_sig:
                sig_frame_rows.append(
                    {
                        "x0": x0 + 0.035,
                        "x1": x1 - 0.035,
                        "y0": y0 + 0.035,
                        "y1": y1 - 0.035,
                    }
                )
            if is_sig or bool(stars):
                use_light = abs(z_oriented) >= 3.20
                y_z_txt = ym + 0.14 if bool(stars) else ym
                z_rec = {"xm": xm, "ym": y_z_txt, "z_txt": f"{z_oriented:+.1f}"}
                if use_light:
                    z_val_light_rows.append(z_rec)
                else:
                    z_val_dark_rows.append(z_rec)

                if stars:
                    star_rec = {"xm": xm, "ym": ym - 0.21, "stars": stars}
                    if use_light:
                        star_light_rows.append(star_rec)
                    else:
                        star_dark_rows.append(star_rec)

    x_dom = [-8.15, 18.95]
    y_dom = [-5.05, float(n_st) + 2.55]
    x_sc = alt.Scale(domain=x_dom, nice=False, zero=False)
    y_sc = alt.Scale(domain=y_dom, nice=False, zero=False)
    blank_ax = _hidden_canvas_axis()

    def _xy_enc(
        x_f: str, y_f: str, x2_f: str = "", y2_f: str = ""
    ) -> dict[str, object]:
        enc: dict[str, object] = {
            "x": alt.X(f"{x_f}:Q", scale=x_sc, axis=blank_ax),
            "y": alt.Y(f"{y_f}:Q", scale=y_sc, axis=blank_ax),
        }
        if x2_f:
            enc["x2"] = alt.X2(f"{x2_f}:Q")
        if y2_f:
            enc["y2"] = alt.Y2(f"{y2_f}:Q")
        return enc

    layers: list[alt.Chart] = []

    if intra_rows:
        layers.append(
            alt.Chart(pl.DataFrame(intra_rows))
            .mark_rect(fill="#F1F5F9", stroke="#E2E8F0", strokeWidth=0.55)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )

    if missing_rows:
        layers.append(
            alt.Chart(pl.DataFrame(missing_rows))
            .mark_rect(fill="#F8FAFC", stroke="#E2E8F0", strokeWidth=0.55)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )

    if cross_rows:
        layers.append(
            alt.Chart(pl.DataFrame(cross_rows))
            .mark_rect(stroke="#CBD5E1", strokeWidth=0.55)
            .encode(
                **_xy_enc("x0", "y0", "x1", "y1"),
                color=alt.Color(
                    "z_acat_clamped:Q",
                    scale=alt.Scale(
                        domain=[-z_cap, 0.0, z_cap],
                        range=[OKABE_VERMILION, "#FFFFFF", OKABE_BLUE],
                        interpolate="rgb",
                    ),
                    legend=alt.Legend(
                        title=["Signed Z_ACAT", "(Row vs. Col)"],
                        orient="right",
                        gradientLength=215,
                        gradientThickness=13,
                        values=[-4.0, -2.0, 0.0, 2.0, 4.0],
                        titleFontSize=9.8,
                        labelFontSize=8.8,
                    ),
                ),
            )
        )

    diag_df = pl.DataFrame(
        [{"x0": 0.0, "y0": float(n_st), "x1": float(n_st), "y1": 0.0}]
    )
    layers.append(
        alt.Chart(diag_df)
        .mark_rule(color="#CBD5E1", strokeWidth=0.9, strokeDash=[4, 3])
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )

    if sig_frame_rows:
        layers.append(
            alt.Chart(pl.DataFrame(sig_frame_rows))
            .mark_rect(fillOpacity=0.0, stroke="#0F172A", strokeWidth=1.35)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )

    if z_val_dark_rows:
        layers.append(
            alt.Chart(pl.DataFrame(z_val_dark_rows))
            .mark_text(
                fontSize=8.8,
                fontWeight="bold",
                fill="#0F172A",
                align="center",
                baseline="middle",
            )
            .encode(**_xy_enc("xm", "ym"), text="z_txt:N")
        )
    if z_val_light_rows:
        layers.append(
            alt.Chart(pl.DataFrame(z_val_light_rows))
            .mark_text(
                fontSize=8.8,
                fontWeight="bold",
                fill="#FFFFFF",
                align="center",
                baseline="middle",
            )
            .encode(**_xy_enc("xm", "ym"), text="z_txt:N")
        )
    if star_dark_rows:
        layers.append(
            alt.Chart(pl.DataFrame(star_dark_rows))
            .mark_text(
                fontSize=9.3,
                fontWeight="bold",
                fill="#0F172A",
                align="center",
                baseline="middle",
            )
            .encode(**_xy_enc("xm", "ym"), text="stars:N")
        )
    if star_light_rows:
        layers.append(
            alt.Chart(pl.DataFrame(star_light_rows))
            .mark_text(
                fontSize=9.3,
                fontWeight="bold",
                fill="#FFFFFF",
                align="center",
                baseline="middle",
            )
            .encode(**_xy_enc("xm", "ym"), text="stars:N")
        )

    comp_spans: list[tuple[TmeCompartment, int, int]] = []
    for st in clustered_states:
        if not comp_spans or comp_spans[-1][0] != st.compartment:
            comp_spans.append((st.compartment, st.order_idx, st.order_idx + 1))
        else:
            prev_c, s_i, _ = comp_spans[-1]
            comp_spans[-1] = (prev_c, s_i, st.order_idx + 1)

    super_spans: list[tuple[TmeSuperGroup, int, int]] = []
    for st in clustered_states:
        if not super_spans or super_spans[-1][0] != st.supergroup:
            super_spans.append((st.supergroup, st.order_idx, st.order_idx + 1))
        else:
            prev_sg, s_i, _ = super_spans[-1]
            super_spans[-1] = (prev_sg, s_i, st.order_idx + 1)

    sub_div_lines: list[dict[str, float]] = []
    for _, _, end_i in comp_spans[:-1]:
        sub_div_lines.append({"x0": float(end_i), "y0": 0.0, "x1": float(end_i), "y1": float(n_st)})
        y_div = float(n_st - end_i)
        sub_div_lines.append({"x0": 0.0, "y0": y_div, "x1": float(n_st), "y1": y_div})

    if sub_div_lines:
        layers.append(
            alt.Chart(pl.DataFrame(sub_div_lines))
            .mark_rule(color="#64748B", strokeWidth=0.85, strokeDash=[3, 2])
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )

    super_div_lines: list[dict[str, float]] = []
    for _, _, end_sg in super_spans[:-1]:
        super_div_lines.append(
            {"x0": float(end_sg), "y0": 0.0, "x1": float(end_sg), "y1": float(n_st)}
        )
        y_sg = float(n_st - end_sg)
        super_div_lines.append({"x0": 0.0, "y0": y_sg, "x1": float(n_st), "y1": y_sg})

    if super_div_lines:
        layers.append(
            alt.Chart(pl.DataFrame(super_div_lines))
            .mark_rule(color="#0F172A", strokeWidth=2.0)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )

    layers.append(
        alt.Chart(pl.DataFrame([{"x0": 0.0, "y0": 0.0, "x1": float(n_st), "y1": float(n_st)}]))
        .mark_rect(fillOpacity=0.0, stroke="#0F172A", strokeWidth=1.15)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )

    top_comp_labels: list[dict[str, object]] = []
    for comp_e, s_i, e_i in comp_spans:
        c_full, c_short, c_hex = _HEATMAP_COMPARTMENT_DISPLAY[comp_e]
        span_w = e_i - s_i
        bar_df = pl.DataFrame(
            [
                {
                    "x0": float(s_i) + 0.03,
                    "x1": float(e_i) - 0.03,
                    "y0": float(n_st) + 0.16,
                    "y1": float(n_st) + 0.48,
                },
                {
                    "x0": -4.60,
                    "x1": -4.32,
                    "y0": float(n_st - e_i) + 0.03,
                    "y1": float(n_st - s_i) - 0.03,
                },
            ]
        )
        layers.append(
            alt.Chart(bar_df)
            .mark_rect(fill=c_hex, stroke="#0F172A", strokeWidth=0.6)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )
        x_lbl = 0.5 * float(s_i + e_i)
        if comp_e == TmeCompartment.STROMAL_MESENCHYMAL and span_w == 1:
            x_lbl -= 0.05
        elif comp_e == TmeCompartment.ENDOTHELIAL_VASCULAR and span_w == 1:
            x_lbl += 0.06
        lbl_txt = (
            f"{c_full} (n={span_w})"
            if span_w >= 4
            else (f"Myeloid (n={span_w})" if comp_e == TmeCompartment.MYELOID_MAST and span_w >= 2 else c_short)
        )
        top_comp_labels.append({"x": x_lbl, "y": float(n_st) + 0.74, "txt": lbl_txt})

    if top_comp_labels:
        layers.append(
            alt.Chart(pl.DataFrame(top_comp_labels))
            .mark_text(
                fontSize=8.7,
                fontWeight="bold",
                fill="#1E293B",
                align="center",
                baseline="middle",
            )
            .encode(**_xy_enc("x", "y"), text="txt:N")
        )

    if dend_segments:
        layers.append(
            alt.Chart(pl.DataFrame(dend_segments))
            .mark_rule(color="#334155", strokeWidth=0.95)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )

    bracket_segs: list[dict[str, float]] = []
    top_sg_labels: list[dict[str, object]] = []
    left_sg_labels: list[dict[str, object]] = []

    y_br_top = float(n_st) + 1.22
    x_br_left = -5.90
    for sg_e, s_i, e_i in super_spans:
        x_a = float(s_i) + 0.05
        x_b = float(e_i) - 0.05
        bracket_segs.extend(
            [
                {"x0": x_a, "y0": y_br_top, "x1": x_b, "y1": y_br_top},
                {"x0": x_a, "y0": y_br_top - 0.16, "x1": x_a, "y1": y_br_top},
                {"x0": x_b, "y0": y_br_top - 0.16, "x1": x_b, "y1": y_br_top},
            ]
        )
        y_a = float(n_st - e_i) + 0.05
        y_b = float(n_st - s_i) - 0.05
        bracket_segs.extend(
            [
                {"x0": x_br_left, "y0": y_a, "x1": x_br_left, "y1": y_b},
                {"x0": x_br_left, "y0": y_a, "x1": x_br_left + 0.16, "y1": y_a},
                {"x0": x_br_left, "y0": y_b, "x1": x_br_left + 0.16, "y1": y_b},
            ]
        )
        span_n = e_i - s_i
        sg_disp = supergroup_display_label(sg_e)
        top_txt = (
            f"{sg_disp} ({span_n} Core Hubs)"
            if span_n >= 4
            else ("Stromal / Endo." if sg_e == TmeSuperGroup.STROMAL_ENDOTHELIAL else sg_disp)
        )
        top_sg_labels.append(
            {"x": 0.5 * float(s_i + e_i), "y": y_br_top + 0.34, "txt": top_txt}
        )

        y_mid_sg = 0.5 * (y_a + y_b)
        if sg_e == TmeSuperGroup.IMMUNE_CELLS:
            left_sg_labels.append({"x": x_br_left - 0.18, "y": y_mid_sg + 0.23, "txt": "Immune"})
            left_sg_labels.append({"x": x_br_left - 0.18, "y": y_mid_sg - 0.23, "txt": "Cells"})
        elif sg_e == TmeSuperGroup.STROMAL_ENDOTHELIAL:
            left_sg_labels.append({"x": x_br_left - 0.18, "y": y_mid_sg + 0.22, "txt": "Stromal /"})
            left_sg_labels.append(
                {"x": x_br_left - 0.18, "y": y_mid_sg - 0.22, "txt": "Endothelial"}
            )
        else:
            left_sg_labels.append({"x": x_br_left - 0.18, "y": y_mid_sg, "txt": sg_disp})

    if bracket_segs:
        layers.append(
            alt.Chart(pl.DataFrame(bracket_segs))
            .mark_rule(color="#0F172A", strokeWidth=1.4)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )
    if top_sg_labels:
        layers.append(
            alt.Chart(pl.DataFrame(top_sg_labels))
            .mark_text(
                fontSize=9.8,
                fontWeight="bold",
                fill="#0F172A",
                align="center",
                baseline="middle",
            )
            .encode(**_xy_enc("x", "y"), text="txt:N")
        )
    if left_sg_labels:
        layers.append(
            alt.Chart(pl.DataFrame(left_sg_labels))
            .mark_text(
                fontSize=9.4,
                fontWeight="bold",
                fill="#0F172A",
                align="right",
                baseline="middle",
            )
            .encode(**_xy_enc("x", "y"), text="txt:N")
        )

    axis_ticks: list[dict[str, float]] = []
    y_lbls: list[dict[str, object]] = []
    x_lbls: list[dict[str, object]] = []

    for st in clustered_states:
        ym = float(n_st) - 0.5 - float(st.order_idx)
        xm = float(st.order_idx) + 0.5
        axis_ticks.append({"x0": -0.14, "y0": ym, "x1": 0.0, "y1": ym})
        axis_ticks.append({"x0": xm, "y0": -0.14, "x1": xm, "y1": 0.0})
        y_lbls.append({"x": -0.22, "y": ym, "txt": st.readable_label})
        x_lbls.append({"x": xm, "y": -0.22, "txt": st.readable_label})

    layers.append(
        alt.Chart(pl.DataFrame(axis_ticks))
        .mark_rule(color="#1E293B", strokeWidth=0.8)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )
    layers.append(
        alt.Chart(pl.DataFrame(y_lbls))
        .mark_text(
            fontSize=9.4, fontWeight="bold", fill="#0F172A", align="right", baseline="middle"
        )
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )
    layers.append(
        alt.Chart(pl.DataFrame(x_lbls))
        .mark_text(
            fontSize=9.2,
            fontWeight="bold",
            fill="#0F172A",
            align="right",
            baseline="middle",
            angle=270,
        )
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )

    panel_hdr_df = pl.DataFrame(
        [
            {
                "x": -7.95,
                "y": float(n_st) + 2.15,
                "txt": "a  Core Pan-Cancer TME Hub Matrix (Signed Z_ACAT & Stratified Permutations)",
            },
            {
                "x": 13.70,
                "y": float(n_st) + 2.15,
                "txt": "b  9-Cohort LOCO Hub Stability",
            },
        ]
    )
    layers.append(
        alt.Chart(panel_hdr_df)
        .mark_text(fontSize=10.8, fontWeight="bold", fill="#0F172A", align="left", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )
    layers.append(
        alt.Chart(
            pl.DataFrame(
                [
                    {
                        "x": 13.70,
                        "y": float(n_st) + 1.66,
                        "txt": "Bar: Peak |Z_ACAT| · ◆: Worst min_LOCO_Z",
                    }
                ]
            )
        )
        .mark_text(fontSize=8.4, fill="#475569", align="left", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )
    layers.append(
        alt.Chart(
            pl.DataFrame(
                [
                    {
                        "x": 0.5 * float(n_st),
                        "y": -4.62,
                        "txt": "Column Hub State j  —  Reference State in Compositional Balance log₂(Row i / Col j)",
                    }
                ]
            )
        )
        .mark_text(fontSize=9.4, fontWeight="bold", fill="#1E293B", align="center", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )

    loco_x0 = 13.75
    loco_scale = 0.86
    loco_x_max = loco_x0 + 5.0 * loco_scale

    layers.append(
        alt.Chart(
            pl.DataFrame([{"x0": loco_x0, "y0": 0.0, "x1": loco_x_max, "y1": float(n_st)}])
        )
        .mark_rect(fill="#FFFFFF", stroke="#CBD5E1", strokeWidth=0.65)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )

    loco_grid = [
        {
            "x0": loco_x0 + zv * loco_scale,
            "y0": 0.0,
            "x1": loco_x0 + zv * loco_scale,
            "y1": float(n_st),
        }
        for zv in (1.0, 2.0, 3.0, 4.0)
    ]
    layers.append(
        alt.Chart(pl.DataFrame(loco_grid))
        .mark_rule(color="#F1F5F9", strokeWidth=0.75)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )

    x_196 = loco_x0 + 1.96 * loco_scale
    layers.append(
        alt.Chart(pl.DataFrame([{"x0": x_196, "y0": 0.0, "x1": x_196, "y1": float(n_st)}]))
        .mark_rule(color=OKABE_VERMILION, strokeWidth=1.05, strokeDash=[4, 3])
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )

    ver_bars: list[dict[str, float]] = []
    loco_pass_pts: list[dict[str, float]] = []
    loco_warn_pts: list[dict[str, float]] = []

    for st in clustered_states:
        ym = float(n_st) - 0.5 - float(st.order_idx)
        z_b = min(
            st.peak_verified_abs_z if st.has_verified_pair else st.peak_any_abs_z,
            5.0,
        )
        ver_bars.append(
            {
                "x0": loco_x0,
                "x1": loco_x0 + z_b * loco_scale,
                "y0": ym - 0.25,
                "y1": ym + 0.25,
            }
        )
        loco_clamped = float(np.clip(st.lead_min_loco_z, 0.0, 5.0))
        pt_rec = {"x": loco_x0 + loco_clamped * loco_scale, "y": ym}
        if st.lead_min_loco_z >= 1.96:
            loco_pass_pts.append(pt_rec)
        else:
            loco_warn_pts.append(pt_rec)

    if ver_bars:
        layers.append(
            alt.Chart(pl.DataFrame(ver_bars))
            .mark_rect(fill="#94A3B8", stroke="#475569", strokeWidth=0.5)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )
    if loco_warn_pts:
        layers.append(
            alt.Chart(pl.DataFrame(loco_warn_pts))
            .mark_point(
                shape="diamond",
                filled=True,
                size=82,
                fill=OKABE_ORANGE,
                stroke="#0F172A",
                strokeWidth=0.85,
            )
            .encode(**_xy_enc("x", "y"))
        )
    if loco_pass_pts:
        layers.append(
            alt.Chart(pl.DataFrame(loco_pass_pts))
            .mark_point(
                shape="diamond",
                filled=True,
                size=88,
                fill=OKABE_BLUISH_GREEN,
                stroke="#0F172A",
                strokeWidth=0.9,
            )
            .encode(**_xy_enc("x", "y"))
        )

    loco_ax_segs = [{"x0": loco_x0, "y0": 0.0, "x1": loco_x_max, "y1": 0.0}]
    loco_ax_lbls: list[dict[str, object]] = []
    for zv in range(6):
        xz = loco_x0 + float(zv) * loco_scale
        loco_ax_segs.append({"x0": xz, "y0": -0.14, "x1": xz, "y1": 0.0})
        loco_ax_lbls.append({"x": xz, "y": -0.35, "txt": str(zv)})

    layers.append(
        alt.Chart(pl.DataFrame(loco_ax_segs))
        .mark_rule(color="#1E293B", strokeWidth=0.8)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )
    layers.append(
        alt.Chart(pl.DataFrame(loco_ax_lbls))
        .mark_text(fontSize=8.4, fill="#1E293B", align="center", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )
    layers.append(
        alt.Chart(
            pl.DataFrame([{"x": 0.5 * (loco_x0 + loco_x_max), "y": -0.72, "txt": "|Z| Score"}])
        )
        .mark_text(fontSize=8.7, fontWeight="bold", fill="#1E293B", align="center", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )

    kx0, kx1 = 13.65, 18.80
    ky0, ky1 = -4.78, -1.10
    layers.append(
        alt.Chart(pl.DataFrame([{"x0": kx0, "y0": ky0, "x1": kx1, "y1": ky1}]))
        .mark_rect(fill="#FFFFFF", stroke="#CBD5E1", strokeWidth=0.7)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )
    key_txt_rows = [
        {"x": kx0 + 0.16, "y": -1.40, "txt": "Stratified Perm. & Lineage Key", "bold": True},
        {"x": kx0 + 0.16, "y": -1.80, "txt": "*** p ≤ 0.001 · ** p < 0.01 · * p < 0.05", "bold": False},
        {"x": kx0 + 0.56, "y": -2.20, "txt": f"ACAT FDR < {fdr_threshold:.2f} (in-cell Z shown)", "bold": False},
        {"x": kx0 + 0.56, "y": -2.60, "txt": "Intra-compartment (gated in P6)", "bold": False},
    ]
    layers.append(
        alt.Chart(pl.DataFrame([r for r in key_txt_rows if bool(r["bold"])]))
        .mark_text(fontSize=8.7, fontWeight="bold", fill="#0F172A", align="left", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )
    layers.append(
        alt.Chart(pl.DataFrame([r for r in key_txt_rows if not bool(r["bold"])]))
        .mark_text(fontSize=8.0, fill="#334155", align="left", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )
    layers.append(
        alt.Chart(
            pl.DataFrame(
                [{"x0": kx0 + 0.17, "x1": kx0 + 0.44, "y0": -2.33, "y1": -2.07}]
            )
        )
        .mark_rect(fill="#FFFFFF", stroke="#0F172A", strokeWidth=1.2)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )
    layers.append(
        alt.Chart(
            pl.DataFrame(
                [{"x0": kx0 + 0.17, "x1": kx0 + 0.44, "y0": -2.73, "y1": -2.47}]
            )
        )
        .mark_rect(fill="#F1F5F9", stroke="#CBD5E1", strokeWidth=0.6)
        .encode(**_xy_enc("x0", "y0", "x1", "y1"))
    )

    swatch_specs = [
        (TmeCompartment.T_NK, kx0 + 0.17, -3.12, "T / NK"),
        (TmeCompartment.B_LINEAGE, kx0 + 0.17, -3.54, "B Lineage"),
        (TmeCompartment.PLASMA, kx0 + 0.17, -3.96, "Plasma"),
        (TmeCompartment.MYELOID_MAST, kx0 + 2.62, -3.12, "Myeloid / Mast"),
        (TmeCompartment.STROMAL_MESENCHYMAL, kx0 + 2.62, -3.54, "Stromal / CAF"),
        (TmeCompartment.ENDOTHELIAL_VASCULAR, kx0 + 2.62, -3.96, "Endothelial"),
    ]
    sw_lbls: list[dict[str, object]] = []
    for cp_s, sx, sy, s_name in swatch_specs:
        c_hex = _HEATMAP_COMPARTMENT_DISPLAY[cp_s][2]
        layers.append(
            alt.Chart(
                pl.DataFrame([{"x0": sx, "x1": sx + 0.27, "y0": sy - 0.13, "y1": sy + 0.13}])
            )
            .mark_rect(fill=c_hex, stroke="#0F172A", strokeWidth=0.55)
            .encode(**_xy_enc("x0", "y0", "x1", "y1"))
        )
        sw_lbls.append({"x": sx + 0.36, "y": sy, "txt": s_name})
    layers.append(
        alt.Chart(pl.DataFrame(sw_lbls))
        .mark_text(fontSize=8.0, fill="#1E293B", align="left", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )
    layers.append(
        alt.Chart(
            pl.DataFrame(
                [
                    {
                        "x": kx0 + 0.17,
                        "y": -4.44,
                        "txt": "◆ Green: min_LOCO_Z ≥ 1.96 · ◆ Orange: < 1.96",
                    }
                ]
            )
        )
        .mark_text(fontSize=7.8, fill="#475569", align="left", baseline="middle")
        .encode(**_xy_enc("x", "y"), text="txt:N")
    )

    n_sig_core = len(
        [
            k_ab
            for k_ab, pr in pair_lookup.items()
            if bool(pr["is_significant"])
            and float(cast(float | int | str, pr["fdr_acat_iut"])) < fdr_threshold
        ]
    )
    combined = (
        alt.layer(*layers)
        .properties(
            width=1060,
            height=806,
            title=alt.TitleParams(
                text="Simplified Pan-Cancer Core TME Hub Heatmap (Key Cell States, Z_ACAT & Permutations)",
                subtitle=(
                    f"Unified 9-Cohort Harmony Atlas (N = 613 Biopsies, 20 CMH Strata) | "
                    f"Top {n_st} Multi-Partner / Perm-Verified Core Hubs (Immune vs. Stromal / Endothelial) | "
                    f"{n_sig_core} Verified Core Pairs"
                ),
                anchor="start",
                fontSize=13.2,
                subtitleFontSize=10.0,
                subtitleColor="#334155",
                offset=10,
            ),
        )
        .configure_view(strokeOpacity=0.0)
    )
    return cast(AltairComposite, apply_nature_methods_theme(combined))


# =============================================================================
# 4B. Pure Functional Core — Broad Cell-Type Pair UMAP Atlas & Composites
# =============================================================================

_BROAD_HOST_COMPARTMENT_ORDER: tuple[TmeCompartment, ...] = (
    TmeCompartment.T_NK,
    TmeCompartment.B_LINEAGE,
    TmeCompartment.PLASMA,
    TmeCompartment.MYELOID_MAST,
    TmeCompartment.STROMAL_MESENCHYMAL,
    TmeCompartment.ENDOTHELIAL_VASCULAR,
)

_BROAD_HOST_COMPARTMENT_META: Mapping[TmeCompartment, tuple[str, str, str]] = {
    TmeCompartment.T_NK: ("T-NK cells", "t_nk", OKABE_BLUE),
    TmeCompartment.B_LINEAGE: ("B cells", "b_cells", OKABE_SKY_BLUE),
    TmeCompartment.PLASMA: ("Plasma cells", "plasma", OKABE_BLUISH_GREEN),
    TmeCompartment.MYELOID_MAST: ("Myeloid/Macrophages", "myeloid_macrophages", OKABE_ORANGE),
    TmeCompartment.STROMAL_MESENCHYMAL: ("CAF", "caf", OKABE_VERMILION),
    TmeCompartment.ENDOTHELIAL_VASCULAR: ("Endothelial", "endothelial", OKABE_REDDISH_PURPLE),
}

_REGIME_DASH_PATTERNS: list[list[int]] = [
    [1, 0],
    [1, 0],
    [5, 3],
    [6, 2, 2, 2],
    [2, 2],
]


def _blend_hex_with_white(hex_color: str, alpha: float) -> str:
    """Blends `#RRGGBB` deterministically toward pure `#FFFFFF` without mutating state."""
    clean = hex_color.strip().lstrip("#")
    if len(clean) != 6:
        return "#E2E8F0"
    a_clamped = float(np.clip(alpha, 0.0, 1.0))
    r_in = int(clean[0:2], 16)
    g_in = int(clean[2:4], 16)
    b_in = int(clean[4:6], 16)
    r_out = int(round(a_clamped * float(r_in) + (1.0 - a_clamped) * 255.0))
    g_out = int(round(a_clamped * float(g_in) + (1.0 - a_clamped) * 255.0))
    b_out = int(round(a_clamped * float(b_in) + (1.0 - a_clamped) * 255.0))
    return f"#{r_out:02X}{g_out:02X}{b_out:02X}"


def _canonicalize_broad_compartment_pair(
    comp_x: TmeCompartment,
    comp_y: TmeCompartment,
) -> tuple[TmeCompartment, TmeCompartment]:
    """Orders two host compartments deterministically according to `_BROAD_HOST_COMPARTMENT_ORDER`."""
    order_map = {cp: idx for idx, cp in enumerate(_BROAD_HOST_COMPARTMENT_ORDER)}
    idx_x = order_map.get(comp_x, 99)
    idx_y = order_map.get(comp_y, 99)
    return (comp_x, comp_y) if idx_x <= idx_y else (comp_y, comp_x)


def enumerate_active_broad_compartment_pairs(
    verified_pan_df: pl.DataFrame,
    fdr_threshold: float = 0.10,
    perm_alpha: float = 0.05,
) -> tuple[BroadCompartmentPairSpec, ...]:
    """Enumerates all host broad cell-type pairs with >= 1 verified or perm-starred interaction."""
    if verified_pan_df.is_empty():
        return ()

    counts_map: dict[tuple[TmeCompartment, TmeCompartment], tuple[int, int]] = {}
    for row in verified_pan_df.to_dicts():
        sa = str(row["state_a"]).strip()
        sb = str(row["state_b"]).strip()
        cp_a = classify_tme_compartment(sa)
        cp_b = classify_tme_compartment(sb)
        if (
            cp_a not in _BROAD_HOST_COMPARTMENT_META
            or cp_b not in _BROAD_HOST_COMPARTMENT_META
            or cp_a == cp_b
        ):
            continue

        fdr_val = float(cast(float | int | str, row["fdr_acat_iut"]))
        p_acat = float(cast(float | int | str, row["pvalue_acat_iut"]))
        emp_p = float(cast(float | int | str, row["emp_perm_pvalue"]))
        gated_perm_p = max(emp_p, p_acat)

        is_fdr = bool(row["is_significant"]) and (fdr_val < fdr_threshold)
        is_perm_only = (not is_fdr) and (gated_perm_p < perm_alpha)
        if not (is_fdr or is_perm_only):
            continue

        key_pair = _canonicalize_broad_compartment_pair(cp_a, cp_b)
        prev_fdr, prev_perm = counts_map.get(key_pair, (0, 0))
        counts_map[key_pair] = (
            prev_fdr + (1 if is_fdr else 0),
            prev_perm + (1 if is_perm_only else 0),
        )

    order_map = {cp: idx for idx, cp in enumerate(_BROAD_HOST_COMPARTMENT_ORDER)}
    ranked_specs: list[tuple[int, int, int, int, int, BroadCompartmentPairSpec]] = []

    for (c1, c2), (n_fdr, n_perm) in counts_map.items():
        if n_fdr + n_perm <= 0:
            continue
        disp1, slug1, col1 = _BROAD_HOST_COMPARTMENT_META[c1]
        disp2, slug2, col2 = _BROAD_HOST_COMPARTMENT_META[c2]
        sg1 = classify_compartment_supergroup(c1)
        sg2 = classify_compartment_supergroup(c2)
        is_cross_super = sg1 != sg2
        macro_cat = (
            "Immune ↔ Stromal / Endothelial"
            if is_cross_super
            else "Immune ↔ Immune"
        )
        macro_prio = 0 if is_cross_super else 1
        spec = BroadCompartmentPairSpec(
            comp_first=c1,
            comp_second=c2,
            slug=f"{slug1}_vs_{slug2}",
            display_first=disp1,
            display_second=disp2,
            display_title=f"{disp1} ↔ {disp2}",
            macro_category=macro_cat,
            color_first=col1,
            color_second=col2,
            n_fdr_pairs=n_fdr,
            n_perm_only_pairs=n_perm,
        )
        ranked_specs.append(
            (macro_prio, -n_fdr, -n_perm, order_map[c1], order_map[c2], spec)
        )

    ranked_specs.sort(key=lambda item: (item[0], item[1], item[2], item[3], item[4]))
    return tuple(item[5] for item in ranked_specs)


def _resolve_oriented_super_nhoods(
    nhood_u: int,
    nhood_v: int,
    state_a: str,
    state_b: str,
    super_meta_codes: Sequence[str],
) -> tuple[int, int]:
    """Maps `(nhood_u, nhood_v)` onto `(nhood_for_state_a, nhood_for_state_b)` deterministically."""
    n_nh = len(super_meta_codes)
    u_ok = 0 <= nhood_u < n_nh
    v_ok = 0 <= nhood_v < n_nh
    if u_ok and v_ok:
        code_u = str(super_meta_codes[nhood_u]).strip()
        code_v = str(super_meta_codes[nhood_v]).strip()
        if code_u == state_a and code_v == state_b:
            return (nhood_u, nhood_v)
        if code_v == state_a and code_u == state_b:
            return (nhood_v, nhood_u)
    return (nhood_u, nhood_v) if state_a <= state_b else (nhood_v, nhood_u)


def _map_dominant_regime_to_visual_polarity(
    dominant_regime: str,
    dominant_direction: str,
) -> Result[tuple[PairwiseVisualRegime, EndpointPolarity, EndpointPolarity], str]:
    """Maps a canonical `(state_a, state_b)` regime + direction to visual regime and `(pol_a, pol_b)`."""
    dir_clean = dominant_direction.strip()
    valid_dirs = {d.value: d for d in EnrichmentDirection}
    if dir_clean not in valid_dirs:
        return Failure(f"Unrecognized dominant_direction: {dominant_direction}")
    dir_enum = valid_dirs[dir_clean]

    reg_clean = dominant_regime.strip()
    if reg_clean.startswith("quadrant:"):
        q_code = reg_clean.split(":", 1)[1]
        valid_quads = {q.value: q for q in QuadrantState}
        if q_code not in valid_quads:
            return Failure(f"Unrecognized quadrant code in dominant_regime: {dominant_regime}")
        return Success(classify_quadrant_visual_regime_and_polarity(valid_quads[q_code], dir_enum))

    valid_axes = {ax.value: ax for ax in ContinuousAxis}
    if reg_clean not in valid_axes:
        return Failure(f"Unrecognized continuous axis in dominant_regime: {dominant_regime}")
    return Success(classify_continuous_visual_regime_and_polarity(valid_axes[reg_clean], dir_enum))


def _compact_meta42_callout_label(
    state_code: str,
    readable_label: str,
    code_only: bool,
) -> str:
    """Formats a concise UMAP callout (`M08 Th2/Treg` or `M08`) from `M08: CCR4+SOCS1+ Th2/Treg`."""
    code_s = state_code.strip()
    if code_only:
        return code_s
    body = readable_label.split(":", 1)[1].strip() if ":" in readable_label else readable_label.strip()
    tokens = [tok for tok in body.split() if tok]
    if not tokens:
        return code_s
    # Drop leading marker-gene token (`CCR4+SOCS1+`, `ISLR+COL5A2+`) when phenotype tokens follow
    pheno_tokens = (
        tokens[1:]
        if len(tokens) >= 2 and ("+" in tokens[0] or "-" in tokens[0])
        else tokens
    )
    tail = " ".join(pheno_tokens).strip()
    if len(tail) > 13:
        tail = tail[:12].rstrip() + "."
    return f"{code_s} {tail}" if tail else code_s


def _compute_repelled_medoid_label_offsets(
    raw_labels: Sequence[Mapping[str, object]],
    umap_bounds: tuple[float, float, float, float],
) -> pl.DataFrame:
    """Computes collision-free 2D UMAP callout positions via deterministic outward repulsion."""
    empty_schema = {
        "state_code": pl.String,
        "label_text": pl.String,
        "orig_u1": pl.Float64,
        "orig_u2": pl.Float64,
        "label_u1": pl.Float64,
        "label_u2": pl.Float64,
        "is_tier1": pl.Boolean,
        "has_leader": pl.Boolean,
    }
    n_lbl = len(raw_labels)
    if n_lbl == 0:
        return pl.DataFrame(schema=empty_schema)

    xmin, xmax, ymin, ymax = umap_bounds
    dx_span = max(xmax - xmin, 1.0)
    dy_span = max(ymax - ymin, 1.0)
    cx_mid = 0.5 * (xmin + xmax)
    cy_mid = 0.5 * (ymin + ymax)

    orig_x = np.array([float(cast(float | int | str, r["orig_u1"])) for r in raw_labels], dtype=np.float64)
    orig_y = np.array([float(cast(float | int | str, r["orig_u2"])) for r in raw_labels], dtype=np.float64)

    # Initial slight outward radial offset away from global manifold center
    vx = orig_x - cx_mid
    vy = orig_y - cy_mid
    v_norm = np.maximum(np.hypot(vx, vy), 1e-6)
    pos_x = orig_x + (vx / v_norm) * (0.028 * dx_span)
    pos_y = orig_y + (vy / v_norm) * (0.044 * dy_span)

    min_sep_x = 0.165 * dx_span
    min_sep_y = 0.062 * dy_span
    pad_x = 0.075 * dx_span
    pad_y = 0.055 * dy_span

    for _step in range(36):
        for i in range(n_lbl):
            for j in range(i + 1, n_lbl):
                ddx = float(pos_x[i] - pos_x[j])
                ddy = float(pos_y[i] - pos_y[j])
                if abs(ddx) < min_sep_x and abs(ddy) < min_sep_y:
                    push_y = 0.52 * (min_sep_y - abs(ddy))
                    sgn_y = (
                        1.0
                        if (ddy > 1e-6 or (abs(ddy) <= 1e-6 and i % 2 == 0))
                        else -1.0
                    )
                    pos_y[i] += sgn_y * push_y
                    pos_y[j] -= sgn_y * push_y
                    push_x = 0.22 * (min_sep_x - abs(ddx))
                    sgn_x = 1.0 if (ddx >= 0.0) else -1.0
                    pos_x[i] += sgn_x * push_x
                    pos_x[j] -= sgn_x * push_x

        pos_x = np.clip(pos_x, xmin + pad_x, xmax - pad_x)
        pos_y = np.clip(pos_y, ymin + pad_y, ymax - pad_y)

    out_rows: list[dict[str, object]] = []
    for idx, r in enumerate(raw_labels):
        disp_norm = float(
            np.hypot(
                (float(pos_x[idx]) - float(orig_x[idx])) / dx_span,
                (float(pos_y[idx]) - float(orig_y[idx])) / dy_span,
            )
        )
        out_rows.append(
            {
                "state_code": str(r["state_code"]),
                "label_text": str(r["label_text"]),
                "orig_u1": float(orig_x[idx]),
                "orig_u2": float(orig_y[idx]),
                "label_u1": float(pos_x[idx]),
                "label_u2": float(pos_y[idx]),
                "is_tier1": bool(r["is_tier1"]),
                "has_leader": disp_norm > 0.036,
            }
        )
    return pl.DataFrame(out_rows)


def _empty_compartment_arc_df() -> pl.DataFrame:
    """Returns typed empty Bezier arc DataFrame for dual-tier compartment UMAP panels."""
    return pl.DataFrame(
        schema={
            "arc_id": pl.String,
            "step_idx": pl.Int64,
            "umap_1": pl.Float64,
            "umap_2": pl.Float64,
            "regime": pl.String,
            "sig_tier": pl.String,
            "pair_code": pl.String,
            "z_acat_signed": pl.Float64,
        }
    )


def _empty_compartment_endpoint_df() -> pl.DataFrame:
    """Returns typed empty endpoint glyph DataFrame for dual-tier compartment UMAP panels."""
    return pl.DataFrame(
        schema={
            "nhood_idx": pl.Int64,
            "umap_1": pl.Float64,
            "umap_2": pl.Float64,
            "state_code": pl.String,
            "polarity": pl.String,
            "sig_tier": pl.String,
            "abs_z": pl.Float64,
        }
    )


def build_compartment_pair_umap_cell_bundle(
    spec: BroadCompartmentPairSpec,
    pan_unit: PreparedSuperUnit,
    umap_coords: FloatArray,
    anchor_compartments: Sequence[str],
    verified_pan_df: pl.DataFrame,
    fdr_threshold: float = 0.10,
    perm_alpha: float = 0.05,
    max_cells: int = 8500,
    max_labels: int = 9,
    compact_labels: bool = False,
    seed: int = 42,
) -> Result[CompartmentPairUmapCellBundle, str]:
    """Builds the pure DataFrame bundle for a single `(comp_first, comp_second)` UMAP panel."""
    coords = np.asarray(umap_coords, dtype=np.float64)
    n_cells = int(coords.shape[0])
    if len(anchor_compartments) != n_cells:
        return Failure(
            f"anchor_compartments length ({len(anchor_compartments)}) != umap_coords rows ({n_cells})"
        )

    med_res = compute_nhood_umap_medoids(
        nhoods_mat=pan_unit.nhoods_super,
        umap_coords=coords,
        state_labels=Some(pan_unit.super_readable_labels),
    )
    match med_res:
        case Failure(err):
            return Failure(err)
        case Success(medoids_df):
            pass

    bounds = compute_umap_bounding_domain(coords)
    u1_lookup = medoids_df["umap_1"].to_numpy()
    u2_lookup = medoids_df["umap_2"].to_numpy()
    global_center = (float(np.median(u1_lookup)), float(np.median(u2_lookup)))
    n_nhoods = medoids_df.height

    pair_records: list[dict[str, object]] = []
    for row in verified_pan_df.to_dicts():
        sa = str(row["state_a"]).strip()
        sb = str(row["state_b"]).strip()
        cp_a = classify_tme_compartment(sa)
        cp_b = classify_tme_compartment(sb)
        if _canonicalize_broad_compartment_pair(cp_a, cp_b) != (
            spec.comp_first,
            spec.comp_second,
        ):
            continue

        fdr_val = float(cast(float | int | str, row["fdr_acat_iut"]))
        p_acat = float(cast(float | int | str, row["pvalue_acat_iut"]))
        emp_p = float(cast(float | int | str, row["emp_perm_pvalue"]))
        gated_perm_p = max(emp_p, p_acat)

        is_fdr = bool(row["is_significant"]) and (fdr_val < fdr_threshold)
        is_perm_only = (not is_fdr) and (gated_perm_p < perm_alpha)
        if not (is_fdr or is_perm_only):
            continue

        nh_a, nh_b = _resolve_oriented_super_nhoods(
            nhood_u=int(cast(int | float | str, row["nhood_u"])),
            nhood_v=int(cast(int | float | str, row["nhood_v"])),
            state_a=sa,
            state_b=sb,
            super_meta_codes=pan_unit.super_meta_codes,
        )
        if not (0 <= nh_a < n_nhoods and 0 <= nh_b < n_nhoods and nh_a != nh_b):
            continue

        pol_res = _map_dominant_regime_to_visual_polarity(
            dominant_regime=str(row["dominant_regime"]),
            dominant_direction=str(row["dominant_direction"]),
        )
        match pol_res:
            case Failure(err):
                return Failure(err)
            case Success((vis_reg, pol_a, pol_b)):
                pass

        pair_records.append(
            {
                "state_a": sa,
                "state_b": sb,
                "state_a_readable": str(row["state_a_readable"]),
                "state_b_readable": str(row["state_b_readable"]),
                "nhood_a": nh_a,
                "nhood_b": nh_b,
                "vis_regime": vis_reg,
                "pol_a": pol_a,
                "pol_b": pol_b,
                "is_fdr": is_fdr,
                "pvalue_acat": p_acat,
                "fdr_acat": fdr_val,
                "gated_perm_p": gated_perm_p,
                "z_signed": float(cast(float | int | str, row["z_acat_signed"])),
                "z_abs": abs(float(cast(float | int | str, row["z_acat_signed"]))),
            }
        )

    pair_records.sort(
        key=lambda rec: (
            0 if bool(rec["is_fdr"]) else 1,
            float(cast(float, rec["pvalue_acat"])),
            -float(cast(float, rec["z_abs"])),
        )
    )

    active_nhoods_first: set[int] = set()
    active_nhoods_second: set[int] = set()
    tier1_arc_rows: list[dict[str, object]] = []
    tier2_arc_rows: list[dict[str, object]] = []
    ep_best: dict[int, tuple[str, EndpointPolarity, ChordSignificanceTier, float]] = {}
    state_best: dict[str, tuple[int, str, bool, float]] = {}

    for rank_i, rec in enumerate(pair_records):
        nh_a = int(cast(int, rec["nhood_a"]))
        nh_b = int(cast(int, rec["nhood_b"]))
        sa = str(rec["state_a"])
        sb = str(rec["state_b"])
        cp_a = classify_tme_compartment(sa)
        if cp_a == spec.comp_first:
            active_nhoods_first.add(nh_a)
            active_nhoods_second.add(nh_b)
        else:
            active_nhoods_first.add(nh_b)
            active_nhoods_second.add(nh_a)

        is_fdr = bool(rec["is_fdr"])
        tier_enum = (
            ChordSignificanceTier.VERIFIED_ACAT_FDR
            if is_fdr
            else ChordSignificanceTier.NOMINAL_PERM_P05
        )
        vis_reg = cast(PairwiseVisualRegime, rec["vis_regime"])
        pol_a = cast(EndpointPolarity, rec["pol_a"])
        pol_b = cast(EndpointPolarity, rec["pol_b"])
        z_sgn = float(cast(float, rec["z_signed"]))
        z_abs = float(cast(float, rec["z_abs"]))

        curv = 0.165 + 0.034 * float(rank_i % 4)
        bx, by = discretize_outward_bezier_arc(
            p_u=(float(u1_lookup[nh_a]), float(u2_lookup[nh_a])),
            p_v=(float(u1_lookup[nh_b]), float(u2_lookup[nh_b])),
            global_center=global_center,
            curvature=curv,
            n_steps=23,
        )
        arc_id = f"{spec.slug}_{rank_i:02d}_{sa}_{sb}"
        target_arc_list = tier1_arc_rows if is_fdr else tier2_arc_rows
        for step_i in range(int(bx.shape[0])):
            target_arc_list.append(
                {
                    "arc_id": arc_id,
                    "step_idx": step_i,
                    "umap_1": float(bx[step_i]),
                    "umap_2": float(by[step_i]),
                    "regime": vis_reg.value,
                    "sig_tier": tier_enum.value,
                    "pair_code": f"{sa} x {sb}",
                    "z_acat_signed": z_sgn,
                }
            )

        for nh_ep, st_ep, read_ep, pol_ep in (
            (nh_a, sa, str(rec["state_a_readable"]), pol_a),
            (nh_b, sb, str(rec["state_b_readable"]), pol_b),
        ):
            if nh_ep not in ep_best:
                ep_best[nh_ep] = (st_ep, pol_ep, tier_enum, z_abs)
            else:
                _prev_st, _prev_pol, prev_tier, prev_z = ep_best[nh_ep]
                upgrade = (
                    prev_tier == ChordSignificanceTier.NOMINAL_PERM_P05
                    and tier_enum == ChordSignificanceTier.VERIFIED_ACAT_FDR
                ) or (prev_tier == tier_enum and z_abs > prev_z)
                if upgrade:
                    ep_best[nh_ep] = (st_ep, pol_ep, tier_enum, z_abs)

            if st_ep not in state_best:
                state_best[st_ep] = (nh_ep, read_ep, is_fdr, z_abs)
            else:
                _p_nh, _p_rd, prev_fdr_flag, prev_sz = state_best[st_ep]
                if (is_fdr and not prev_fdr_flag) or (is_fdr == prev_fdr_flag and z_abs > prev_sz):
                    state_best[st_ep] = (nh_ep, read_ep, is_fdr or prev_fdr_flag, max(z_abs, prev_sz))

    # Classify single cells across 5 ordered lineage/interaction footprint tiers
    csc_super = sp.csc_matrix(pan_unit.nhoods_super)
    in_act_first = np.zeros(n_cells, dtype=bool)
    for nh_id in active_nhoods_first:
        s_ptr, e_ptr = int(csc_super.indptr[nh_id]), int(csc_super.indptr[nh_id + 1])
        pos_c = csc_super.indices[s_ptr:e_ptr][csc_super.data[s_ptr:e_ptr] > 0]
        in_act_first[pos_c] = True

    in_act_second = np.zeros(n_cells, dtype=bool)
    for nh_id in active_nhoods_second:
        s_ptr, e_ptr = int(csc_super.indptr[nh_id]), int(csc_super.indptr[nh_id + 1])
        pos_c = csc_super.indices[s_ptr:e_ptr][csc_super.data[s_ptr:e_ptr] > 0]
        in_act_second[pos_c] = True

    comp_arr = np.array([str(c).strip().lower() for c in anchor_compartments], dtype=str)
    is_c1 = comp_arr == spec.comp_first.value.lower()
    is_c2 = comp_arr == spec.comp_second.value.lower()

    role_other = "Other Compartments"
    role_ctx_1 = f"{spec.display_first} (Context)"
    role_ctx_2 = f"{spec.display_second} (Context)"
    role_act_1 = f"{spec.display_first} (Interacting)"
    role_act_2 = f"{spec.display_second} (Interacting)"

    draw_rank = np.zeros(n_cells, dtype=np.int64)
    draw_rank[is_c1] = 1
    draw_rank[is_c2] = 2
    draw_rank[is_c1 & in_act_first] = 3
    draw_rank[is_c2 & in_act_second] = 4

    role_lookup = (role_other, role_ctx_1, role_ctx_2, role_act_1, role_act_2)
    cell_roles = [role_lookup[int(rk)] for rk in draw_rank.tolist()]

    if n_cells > max_cells:
        rng = np.random.default_rng(seed)
        fg_idx = np.flatnonzero(draw_rank >= 3)
        ctx_idx = np.flatnonzero((draw_rank == 1) | (draw_rank == 2))
        bg_idx = np.flatnonzero(draw_rank == 0)

        n_fg_take = min(int(fg_idx.shape[0]), int(0.48 * max_cells))
        rem_after_fg = max(max_cells - n_fg_take, 200)
        n_ctx_take = min(int(ctx_idx.shape[0]), int(0.45 * rem_after_fg))
        n_bg_take = min(int(bg_idx.shape[0]), max(rem_after_fg - n_ctx_take, 0))

        pick_fg = (
            rng.choice(fg_idx, size=n_fg_take, replace=False)
            if int(fg_idx.shape[0]) > n_fg_take
            else fg_idx
        )
        pick_ctx = (
            rng.choice(ctx_idx, size=n_ctx_take, replace=False)
            if int(ctx_idx.shape[0]) > n_ctx_take
            else ctx_idx
        )
        pick_bg = (
            rng.choice(bg_idx, size=n_bg_take, replace=False)
            if int(bg_idx.shape[0]) > n_bg_take
            else bg_idx
        )
        keep_idx = np.concatenate([pick_bg, pick_ctx, pick_fg])
    else:
        keep_idx = np.arange(n_cells, dtype=np.int64)

    # Sort strictly from bottom (rank 0) to top (rank 4)
    order_sub = np.argsort(draw_rank[keep_idx], kind="mergesort")
    sorted_idx = keep_idx[order_sub]
    cells_df = pl.DataFrame(
        {
            "cell_idx": sorted_idx,
            "umap_1": coords[sorted_idx, 0],
            "umap_2": coords[sorted_idx, 1],
            "cell_footprint_role": [cell_roles[int(i)] for i in sorted_idx.tolist()],
            "draw_rank": draw_rank[sorted_idx],
        }
    )

    ep_rows = [
        {
            "nhood_idx": nh_id,
            "umap_1": float(u1_lookup[nh_id]),
            "umap_2": float(u2_lookup[nh_id]),
            "state_code": st_c,
            "polarity": pol_e.value,
            "sig_tier": tier_e.value,
            "abs_z": z_v,
        }
        for nh_id, (st_c, pol_e, tier_e, z_v) in sorted(
            ep_best.items(),
            key=lambda kv: (
                0 if kv[1][2] == ChordSignificanceTier.NOMINAL_PERM_P05 else 1,
                kv[1][3],
            ),
        )
    ]
    endpoints_df = pl.DataFrame(ep_rows) if ep_rows else _empty_compartment_endpoint_df()

    sorted_states = sorted(
        state_best.items(),
        key=lambda kv: (0 if kv[1][2] else 1, -kv[1][3], kv[0]),
    )[:max_labels]
    use_code_only = compact_labels and len(sorted_states) > 6
    raw_lbl_items: list[dict[str, object]] = [
        {
            "state_code": st_code,
            "label_text": _compact_meta42_callout_label(st_code, read_lbl, code_only=use_code_only),
            "orig_u1": float(u1_lookup[nh_rep]),
            "orig_u2": float(u2_lookup[nh_rep]),
            "is_tier1": is_t1,
        }
        for st_code, (nh_rep, read_lbl, is_t1, _zv) in sorted_states
    ]
    medoid_labels_df = _compute_repelled_medoid_label_offsets(raw_lbl_items, bounds)

    if pair_records:
        lead = pair_records[0]
        stars = format_stratified_perm_stars(float(cast(float, lead["gated_perm_p"])))
        lead_summary = (
            f"Lead: {lead['state_a']}×{lead['state_b']} "
            f"(Z={float(cast(float, lead['z_signed'])):+.2f}{stars}, q={float(cast(float, lead['fdr_acat'])):.2g})"
        )
    else:
        lead_summary = "No significant pairs"

    return Success(
        CompartmentPairUmapCellBundle(
            spec=spec,
            umap_bounds=bounds,
            cells_df=cells_df,
            tier1_arcs_df=pl.DataFrame(tier1_arc_rows)
            if tier1_arc_rows
            else _empty_compartment_arc_df(),
            tier2_arcs_df=pl.DataFrame(tier2_arc_rows)
            if tier2_arc_rows
            else _empty_compartment_arc_df(),
            endpoints_df=endpoints_df,
            medoid_labels_df=medoid_labels_df,
            lead_pair_summary=lead_summary,
        )
    )


def build_compartment_pair_pairwise_umap_bundle(
    spec: BroadCompartmentPairSpec,
    pan_unit: PreparedSuperUnit,
    umap_coords: FloatArray,
    anchor_cell_labels: Sequence[str],
    global_response_labels: Sequence[str],
    verified_pan_df: pl.DataFrame,
    fdr_threshold: float = 0.10,
    perm_alpha: float = 0.05,
) -> Result[PairwiseUmapBundle, str]:
    """Builds a compartment-pair-scoped `PairwiseUmapBundle` for standalone 2x2 panels `b`, `c`, `d`."""
    read_lbls = pan_unit.super_readable_labels
    n_nhoods = len(read_lbls)

    q_rows: list[dict[str, object]] = []
    c_rows: list[dict[str, object]] = []
    hub_z_totals: dict[str, float] = {}
    state_to_readable: dict[str, str] = {}
    lead_pair_uv: Maybe[tuple[int, int]] = Nothing
    best_lead_rank: tuple[int, float] = (99, 2.0)

    # Also store all testable Z_ACAT scores within (comp_first, comp_second) for smooth Panel c field
    pair_signed_z_lookup: dict[tuple[str, str], tuple[str, float, str]] = {}

    for row in verified_pan_df.to_dicts():
        sa = str(row["state_a"]).strip()
        sb = str(row["state_b"]).strip()
        cp_a = classify_tme_compartment(sa)
        cp_b = classify_tme_compartment(sb)
        if _canonicalize_broad_compartment_pair(cp_a, cp_b) != (
            spec.comp_first,
            spec.comp_second,
        ):
            continue

        nh_a, nh_b = _resolve_oriented_super_nhoods(
            nhood_u=int(cast(int | float | str, row["nhood_u"])),
            nhood_v=int(cast(int | float | str, row["nhood_v"])),
            state_a=sa,
            state_b=sb,
            super_meta_codes=pan_unit.super_meta_codes,
        )
        if not (0 <= nh_a < n_nhoods and 0 <= nh_b < n_nhoods and nh_a != nh_b):
            continue

        dom_reg = str(row["dominant_regime"]).strip()
        dom_dir = str(row["dominant_direction"]).strip()
        dom_eff = float(cast(float | int | str, row["dominant_effect"]))
        p_acat = float(cast(float | int | str, row["pvalue_acat_iut"]))
        fdr_val = float(cast(float | int | str, row["fdr_acat_iut"]))
        emp_p = float(cast(float | int | str, row["emp_perm_pvalue"]))
        z_sgn = float(cast(float | int | str, row["z_acat_signed"]))
        gated_perm_p = max(emp_p, p_acat)

        pair_signed_z_lookup[(sa, sb)] = (sa, z_sgn, dom_reg)
        state_to_readable[sa] = str(row["state_a_readable"])
        state_to_readable[sb] = str(row["state_b_readable"])

        is_fdr = bool(row["is_significant"]) and (fdr_val < fdr_threshold)
        is_perm_only = (not is_fdr) and (gated_perm_p < perm_alpha)
        is_active_hit = is_fdr or is_perm_only

        if is_active_hit:
            w_add = abs(z_sgn) * (1.5 if is_fdr else 0.75)
            hub_z_totals[sa] = hub_z_totals.get(sa, 0.0) + w_add
            hub_z_totals[sb] = hub_z_totals.get(sb, 0.0) + w_add
            cand_rank = (0 if is_fdr else 1, p_acat)
            if cand_rank < best_lead_rank:
                best_lead_rank = cand_rank
                lead_pair_uv = Some((nh_a, nh_b))

        if dom_reg.startswith("quadrant:"):
            q_code = dom_reg.split(":", 1)[1]
            q_rows.append(
                {
                    "nhood_u": nh_a,
                    "nhood_v": nh_b,
                    "cell_type_u": read_lbls[nh_a],
                    "cell_type_v": read_lbls[nh_b],
                    "quadrant": q_code,
                    "direction": dom_dir,
                    "log_or_quadrant": dom_eff,
                    "pvalue_iut": p_acat,
                    "fdr_iut": fdr_val,
                    "is_significant": is_active_hit,
                }
            )
        else:
            c_rows.append(
                {
                    "nhood_u": nh_a,
                    "nhood_v": nh_b,
                    "cell_type_u": read_lbls[nh_a],
                    "cell_type_v": read_lbls[nh_b],
                    "axis": dom_reg,
                    "direction": dom_dir,
                    "effect_size": dom_eff,
                    "pvalue_iut": p_acat,
                    "fdr_iut": fdr_val,
                    "is_significant": is_active_hit,
                }
            )

    new_q_df = (
        pl.DataFrame(q_rows)
        if q_rows
        else pl.DataFrame(
            schema={
                "nhood_u": pl.Int64,
                "nhood_v": pl.Int64,
                "cell_type_u": pl.String,
                "cell_type_v": pl.String,
                "quadrant": pl.String,
                "direction": pl.String,
                "log_or_quadrant": pl.Float64,
                "pvalue_iut": pl.Float64,
                "fdr_iut": pl.Float64,
                "is_significant": pl.Boolean,
            }
        )
    )
    new_c_df = (
        pl.DataFrame(c_rows)
        if c_rows
        else pl.DataFrame(
            schema={
                "nhood_u": pl.Int64,
                "nhood_v": pl.Int64,
                "cell_type_u": pl.String,
                "cell_type_v": pl.String,
                "axis": pl.String,
                "direction": pl.String,
                "effect_size": pl.Float64,
                "pvalue_iut": pl.Float64,
                "fdr_iut": pl.Float64,
                "is_significant": pl.Boolean,
            }
        )
    )

    n_sig_q = int(new_q_df.filter(pl.col("is_significant")).height)
    n_sig_c = int(new_c_df.filter(pl.col("is_significant")).height)
    pw_viz = PairwiseAnalysisResult(
        pairs_df=new_q_df,
        continuous_pairs_df=new_c_df,
        modules_df=pl.DataFrame(),
        n_total_nhoods=n_nhoods,
        n_active_nhoods=n_nhoods,
        n_evaluated_pairs=new_q_df.height + new_c_df.height,
        n_significant_quadrant_pairs=n_sig_q,
        n_significant_continuous_pairs=n_sig_c,
        n_modules=n_sig_q + n_sig_c,
    )

    anchor_state_maybe: Maybe[str] = Nothing
    chosen_anchor_code = ""
    if hub_z_totals:
        chosen_anchor_code = sorted(hub_z_totals.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        if chosen_anchor_code in state_to_readable:
            anchor_state_maybe = Some(state_to_readable[chosen_anchor_code])

    univ_fdr, univ_lfc = compute_stratified_pan_univariate_da(pan_unit, global_response_labels)
    viz_cfg = UmapPairVizConfig(
        max_chords=32,
        max_bg_cells=11000,
        fdr_threshold=fdr_threshold,
        panel_width=305,
        panel_height=260,
    )
    base_res = prepare_pairwise_umap_bundle(
        nhoods_mat=pan_unit.nhoods_super,
        umap_coords=umap_coords,
        count_matrix=pan_unit.counts_super,
        response_labels=global_response_labels,
        pairwise_res=pw_viz,
        sample_library_sizes=Some(tuple(pan_unit.library_sizes.tolist())),
        state_labels=Some(read_lbls),
        cell_labels=Some(anchor_cell_labels),
        univariate_fdr=Some(univ_fdr),
        univariate_logfc=Some(univ_lfc),
        anchor_state=anchor_state_maybe,
        spotlight_pair=lead_pair_uv,
        config=viz_cfg,
    )
    match base_res:
        case Failure(err):
            return Failure(err)
        case Success(bundle):
            pass

    # Project exact heatmap-aligned directional Z_ACAT(Partner vs. Anchor) onto Panel c cells
    updated_anchor_bg = bundle.anchor_bg_cells_df
    if chosen_anchor_code and not updated_anchor_bg.is_empty():
        nhood_dir_z = np.zeros(n_nhoods, dtype=np.float64)
        for nh_i, st_i in enumerate(pan_unit.super_meta_codes):
            st_clean = str(st_i).strip()
            k_ord = (min(st_clean, chosen_anchor_code), max(st_clean, chosen_anchor_code))
            if k_ord in pair_signed_z_lookup:
                ref_a, z_ab, reg_s = pair_signed_z_lookup[k_ord]
                nhood_dir_z[nh_i] = _oriented_pair_z_acat(
                    row_state=st_clean,
                    col_state=chosen_anchor_code,
                    state_a=ref_a,
                    z_signed_ab=z_ab,
                    dominant_regime=reg_s,
                )
        bin_mat = (sp.csc_matrix(pan_unit.nhoods_super) > 0).astype(np.float64)
        bg_c_idx = updated_anchor_bg["cell_idx"].to_numpy()
        bin_sub = bin_mat[bg_c_idx, :]
        deg_sub = np.maximum(np.asarray(bin_sub.sum(axis=1)).ravel(), 1.0)
        exact_cell_z = np.asarray(bin_sub.dot(nhood_dir_z)).ravel() / deg_sub
        updated_anchor_bg = (
            updated_anchor_bg.with_columns(
                pl.Series("anchor_partner_z", exact_cell_z),
                pl.Series("_abs_z_ord", np.abs(exact_cell_z)),
            )
            .sort("_abs_z_ord")
            .drop("_abs_z_ord")
        )

    # Recompute observable-biopsy 2x2 donor responder inset for `(u*, v*)`
    updated_quad_df = bundle.spotlight_quadrant_df
    if bundle.spotlight_meta.has_pair:
        u_star = bundle.spotlight_meta.nhood_u
        v_star = bundle.spotlight_meta.nhood_v
        pres_all = estimate_nhood_presence(
            count_matrix=pan_unit.counts_super,
            library_sizes=Some(pan_unit.library_sizes),
            config=PresenceConfig(mode=PresenceMode.EMPIRICAL_BAYES_POSTERIOR),
        ).unwrap()
        pres_masked = build_compartment_observability_mask(
            presence=pres_all,
            nhood_labels=pan_unit.super_meta_codes,
            cohort_labels=pan_unit.cohort_labels,
            min_cohort_compartment_share=0.001,
        ).value_or(pres_all)
        both_obs = (pres_masked.observability[u_star, :] > 0.5) & (
            pres_masked.observability[v_star, :] > 0.5
        )
        obs_idx = np.flatnonzero(both_obs)
        if int(obs_idx.shape[0]) >= 8:
            u_p = pres_masked.presence_prob[u_star, obs_idx] >= 0.5
            v_p = pres_masked.presence_prob[v_star, obs_idx] >= 0.5
            if int(np.sum(u_p & v_p)) == int(obs_idx.shape[0]) or int(np.sum(u_p | v_p)) == 0:
                u_cpm = pres_masked.log2_cpm[u_star, obs_idx]
                v_cpm = pres_masked.log2_cpm[v_star, obs_idx]
                u_p = u_cpm > float(np.median(u_cpm))
                v_p = v_cpm > float(np.median(v_cpm))
            y_obs = np.array(
                [1 if global_response_labels[int(j)] == "responder" else 0 for j in obs_idx.tolist()],
                dtype=np.int64,
            )
            q_defs = (
                ("u+ / v+", u_p & v_p),
                ("u+ / v-", u_p & (~v_p)),
                ("u- / v+", (~u_p) & v_p),
                ("u- / v-", (~u_p) & (~v_p)),
            )
            q_rows_obs: list[dict[str, object]] = []
            for idx_q, (q_lbl, q_m) in enumerate(q_defs):
                n_t = int(np.sum(q_m))
                n_r = int(np.sum(y_obs[q_m] == 1)) if n_t > 0 else 0
                pct = 100.0 * float(n_r) / float(max(n_t, 1))
                lo, hi = wilson_binomial_ci_pct(n_r, n_t)
                q_rows_obs.append(
                    {
                        "donor_quadrant": q_lbl,
                        "quad_order": idx_q,
                        "n_responders": n_r,
                        "n_total": n_t,
                        "responder_pct": pct,
                        "ci_low_pct": lo,
                        "ci_high_pct": hi,
                    }
                )
            updated_quad_df = pl.DataFrame(q_rows_obs)

    return Success(
        PairwiseUmapBundle(
            umap_bounds=bundle.umap_bounds,
            nhood_nodes_df=bundle.nhood_nodes_df,
            hub_overlay_df=bundle.hub_overlay_df,
            bg_cells_df=bundle.bg_cells_df,
            bundled_arcs_df=bundle.bundled_arcs_df,
            arc_endpoints_df=bundle.arc_endpoints_df,
            anchor_label=bundle.anchor_label,
            anchor_bg_cells_df=updated_anchor_bg,
            anchor_overlay_df=bundle.anchor_overlay_df,
            anchor_arcs_df=bundle.anchor_arcs_df,
            spotlight_cells_df=bundle.spotlight_cells_df,
            spotlight_arc_df=bundle.spotlight_arc_df,
            spotlight_endpoints_df=bundle.spotlight_endpoints_df,
            spotlight_quadrant_df=updated_quad_df,
            spotlight_meta=bundle.spotlight_meta,
            n_sig_quadrant=bundle.n_sig_quadrant,
            n_sig_continuous=bundle.n_sig_continuous,
            n_denovo_nhoods=bundle.n_denovo_nhoods,
        )
    )


def build_compartment_pair_umap_panel(
    cell_bundle: CompartmentPairUmapCellBundle,
    panel_prefix: str = "a",
    width: int = 305,
    height: int = 260,
    compact_grid_mode: bool = False,
) -> alt.LayerChart:
    """Renders a single *Nature Methods* dual-tier chord + lineage cell-footprint UMAP panel."""
    spec = cell_bundle.spec
    bounds = cell_bundle.umap_bounds

    x_enc = alt.X(
        "umap_1:Q",
        scale=alt.Scale(domain=[bounds[0], bounds[1]], nice=False),
        axis=create_nature_axis("UMAP 1"),
    )
    y_enc = alt.Y(
        "umap_2:Q",
        scale=alt.Scale(domain=[bounds[2], bounds[3]], nice=False),
        axis=create_nature_axis("UMAP 2"),
    )

    role_other = "Other Compartments"
    role_ctx_1 = f"{spec.display_first} (Context)"
    role_ctx_2 = f"{spec.display_second} (Context)"
    role_act_1 = f"{spec.display_first} (Interacting)"
    role_act_2 = f"{spec.display_second} (Interacting)"

    role_domain = [role_act_1, role_act_2, role_ctx_1, role_ctx_2, role_other]
    role_colors = [
        spec.color_first,
        spec.color_second,
        _blend_hex_with_white(spec.color_first, 0.28),
        _blend_hex_with_white(spec.color_second, 0.28),
        "#E8EEF5",
    ]

    pt_size = 6.5 if compact_grid_mode else 9.0
    cells_color_enc = (
        alt.Color(
            "cell_footprint_role:N",
            scale=alt.Scale(domain=role_domain, range=role_colors),
        )
        if compact_grid_mode
        else alt.Color(
            "cell_footprint_role:N",
            scale=alt.Scale(domain=role_domain, range=role_colors),
            legend=alt.Legend(
                title="Lineage Cell Footprints",
                orient="bottom",
                columns=2,
                values=[role_act_1, role_act_2, role_ctx_1, role_ctx_2],
                symbolSize=55,
                labelLimit=220,
            ),
        )
    )
    cells_layer = (
        alt.Chart(cell_bundle.cells_df)
        .mark_circle(size=pt_size, opacity=0.84)
        .encode(
            x=x_enc,
            y=y_enc,
            color=cells_color_enc,
        )
    )

    regime_domain = [reg.value for reg in ALL_VISUAL_REGIMES]
    regime_colors = [regime_color_hex(reg) for reg in ALL_VISUAL_REGIMES]

    sub_layers: list[alt.Chart] = [cells_layer]

    # Unified dual-tier Bezier chord layer (Tier 2 drawn first, Tier 1 bold on top)
    tier_fdr_lbl = "ACAT FDR < 0.10"
    tier_perm_lbl = "Perm. p < 0.05"
    tier_domain = [tier_fdr_lbl, tier_perm_lbl]
    sw_range = [1.95, 0.95] if compact_grid_mode else [2.25, 1.15]
    op_range = [0.92, 0.46]

    arc_frames: list[pl.DataFrame] = []
    if not cell_bundle.tier2_arcs_df.is_empty():
        arc_frames.append(
            cell_bundle.tier2_arcs_df.with_columns(
                (pl.lit("t2_") + pl.col("arc_id").cast(pl.Utf8)).alias("arc_uid"),
                pl.lit(tier_perm_lbl).alias("tier_label"),
            )
        )
    if not cell_bundle.tier1_arcs_df.is_empty():
        arc_frames.append(
            cell_bundle.tier1_arcs_df.with_columns(
                (pl.lit("t1_") + pl.col("arc_id").cast(pl.Utf8)).alias("arc_uid"),
                pl.lit(tier_fdr_lbl).alias("tier_label"),
            )
        )
    if arc_frames:
        combined_arcs_df = pl.concat(arc_frames, how="vertical_relaxed")
        chord_color_enc = (
            alt.Color("regime:N", scale=alt.Scale(domain=regime_domain, range=regime_colors))
            if compact_grid_mode
            else alt.Color(
                "regime:N",
                scale=alt.Scale(domain=regime_domain, range=regime_colors),
                legend=alt.Legend(title="Pairwise Regime", orient="bottom", columns=2, labelLimit=180),
            )
        )
        chord_dash_enc = (
            alt.StrokeDash(
                "regime:N",
                scale=alt.Scale(domain=regime_domain, range=_REGIME_DASH_PATTERNS),
            )
            if compact_grid_mode
            else alt.StrokeDash(
                "regime:N",
                scale=alt.Scale(domain=regime_domain, range=_REGIME_DASH_PATTERNS),
                legend=alt.Legend(title="Pairwise Regime", orient="bottom", columns=2, labelLimit=180),
            )
        )
        chord_sw_enc = (
            alt.StrokeWidth("tier_label:N", scale=alt.Scale(domain=tier_domain, range=sw_range))
            if compact_grid_mode
            else alt.StrokeWidth(
                "tier_label:N",
                scale=alt.Scale(domain=tier_domain, range=sw_range),
                legend=alt.Legend(title="Significance Tier", orient="bottom", columns=1, labelLimit=150),
            )
        )
        chord_op_enc = (
            alt.Opacity("tier_label:N", scale=alt.Scale(domain=tier_domain, range=op_range))
            if compact_grid_mode
            else alt.Opacity(
                "tier_label:N",
                scale=alt.Scale(domain=tier_domain, range=op_range),
                legend=alt.Legend(title="Significance Tier", orient="bottom", columns=1, labelLimit=150),
            )
        )
        sub_layers.append(
            alt.Chart(combined_arcs_df)
            .mark_line()
            .encode(
                x=x_enc,
                y=y_enc,
                detail="arc_uid:N",
                order="step_idx:Q",
                color=chord_color_enc,
                strokeDash=chord_dash_enc,
                strokeWidth=chord_sw_enc,
                opacity=chord_op_enc,
            )
        )

    # Compact inline lineage key in guaranteed-empty top-left UMAP corner for 4x3 grid panels
    if compact_grid_mode:
        bx_dot = bounds[0] + 0.85
        bx_txt = bounds[0] + 1.55
        by_1 = bounds[3] - 0.95
        by_2 = bounds[3] - 2.35
        sub_layers.append(
            alt.Chart(pl.DataFrame([{"u1": bx_dot, "u2": by_1}]))
            .mark_circle(size=32, color=spec.color_first, stroke="#FFFFFF", strokeWidth=0.6)
            .encode(x="u1:Q", y="u2:Q")
        )
        sub_layers.append(
            alt.Chart(pl.DataFrame([{"u1": bx_dot, "u2": by_2}]))
            .mark_circle(size=32, color=spec.color_second, stroke="#FFFFFF", strokeWidth=0.6)
            .encode(x="u1:Q", y="u2:Q")
        )
        sub_layers.append(
            alt.Chart(
                pl.DataFrame(
                    [
                        {"u1": bx_txt, "u2": by_1, "lbl": spec.display_first},
                        {"u1": bx_txt, "u2": by_2, "lbl": spec.display_second},
                    ]
                )
            )
            .mark_text(fontSize=7.3, fontWeight="bold", fill="#1E293B", align="left", baseline="middle")
            .encode(x="u1:Q", y="u2:Q", text="lbl:N")
        )

    # Polar endpoint glyphs (+ filled Okabe Blue, - hollow Vermilion ring, -- Orange diamond)
    if not cell_bundle.endpoints_df.is_empty():
        for tier_val, sz_scale, ep_op in (
            (ChordSignificanceTier.NOMINAL_PERM_P05.value, 0.56, 0.72),
            (ChordSignificanceTier.VERIFIED_ACAT_FDR.value, 1.00, 0.98),
        ):
            ep_sub = cell_bundle.endpoints_df.filter(pl.col("sig_tier") == tier_val)
            if ep_sub.is_empty():
                continue
            base_sz = (36.0 if compact_grid_mode else 50.0) * sz_scale
            pos_ep = ep_sub.filter(
                pl.col("polarity").is_in(
                    [
                        EndpointPolarity.POSITIVE_PRESENT.value,
                        EndpointPolarity.SYMMETRIC_COABUND.value,
                    ]
                )
            )
            neg_ep = ep_sub.filter(pl.col("polarity") == EndpointPolarity.NEGATIVE_GATED.value)
            des_ep = ep_sub.filter(pl.col("polarity") == EndpointPolarity.DUAL_ABSENT.value)

            if not pos_ep.is_empty():
                sub_layers.append(
                    alt.Chart(pos_ep)
                    .mark_circle(
                        size=base_sz,
                        opacity=ep_op,
                        color=OKABE_BLUE,
                        stroke="#FFFFFF",
                        strokeWidth=0.85,
                    )
                    .encode(x="umap_1:Q", y="umap_2:Q")
                )
            if not neg_ep.is_empty():
                sub_layers.append(
                    alt.Chart(neg_ep)
                    .mark_point(
                        shape="circle",
                        filled=True,
                        size=base_sz * 1.12,
                        fill="#FFFFFF",
                        stroke=OKABE_VERMILION,
                        strokeWidth=1.75 if sz_scale > 0.8 else 1.25,
                        opacity=ep_op,
                    )
                    .encode(x="umap_1:Q", y="umap_2:Q")
                )
            if not des_ep.is_empty():
                sub_layers.append(
                    alt.Chart(des_ep)
                    .mark_point(
                        shape="diamond",
                        filled=True,
                        size=base_sz,
                        fill=OKABE_ORANGE,
                        stroke="#0F172A",
                        strokeWidth=0.75,
                        opacity=ep_op,
                    )
                    .encode(x="umap_1:Q", y="umap_2:Q")
                )

    # Collision-free Meta42 state callouts with white text halo
    if not cell_bundle.medoid_labels_df.is_empty():
        leaders_df = cell_bundle.medoid_labels_df.filter(pl.col("has_leader"))
        if not leaders_df.is_empty():
            sub_layers.append(
                alt.Chart(leaders_df)
                .mark_rule(color="#475569", strokeWidth=0.65, opacity=0.75)
                .encode(
                    x="orig_u1:Q",
                    y="orig_u2:Q",
                    x2="label_u1:Q",
                    y2="label_u2:Q",
                )
            )
        f_sz = 7.6 if compact_grid_mode else 8.4
        t1_lbls = cell_bundle.medoid_labels_df.filter(pl.col("is_tier1"))
        t2_lbls = cell_bundle.medoid_labels_df.filter(~pl.col("is_tier1"))

        tier_label_specs: tuple[tuple[pl.DataFrame, Literal["normal", "bold"], str], ...] = (
            (t2_lbls, "normal", "#334155"),
            (t1_lbls, "bold", "#0F172A"),
        )
        for lbl_df, fw, fg_hex in tier_label_specs:
            if lbl_df.is_empty():
                continue
            sub_layers.append(
                alt.Chart(lbl_df)
                .mark_text(
                    fontSize=f_sz,
                    fontWeight=fw,
                    fill=fg_hex,
                    stroke="#FFFFFF",
                    strokeWidth=2.4,
                    align="center",
                    baseline="middle",
                )
                .encode(x="label_u1:Q", y="label_u2:Q", text="label_text:N")
            )
            sub_layers.append(
                alt.Chart(lbl_df)
                .mark_text(
                    fontSize=f_sz,
                    fontWeight=fw,
                    fill=fg_hex,
                    align="center",
                    baseline="middle",
                )
                .encode(x="label_u1:Q", y="label_u2:Q", text="label_text:N")
            )

    title_str = (
        f"{panel_prefix}. {spec.display_title}"
        if compact_grid_mode
        else f"{panel_prefix}. {spec.display_title} — Cross-Compartment Chords & Footprints"
    )
    subtitle_str = (
        f"{spec.n_fdr_pairs} FDR<0.10 (bold) + {spec.n_perm_only_pairs} perm p<0.05 (thin) | "
        f"{cell_bundle.lead_pair_summary}"
    )

    return cast(
        alt.LayerChart,
        alt.layer(*sub_layers)
        .resolve_scale(
            color="independent",
            strokeDash="independent",
            strokeWidth="independent",
            opacity="independent",
        )
        .properties(
            width=width,
            height=height,
            title=alt.TitleParams(
                text=title_str,
                subtitle=subtitle_str,
                fontSize=10.3 if compact_grid_mode else 11.5,
                fontWeight="bold",
                subtitleFontSize=8.0 if compact_grid_mode else 8.8,
                subtitleColor="#475569",
                anchor="start",
            ),
        ),
    )


def _build_compartment_grid_legend_card(
    width: int = 268,
    height: int = 228,
) -> alt.LayerChart:
    """Builds the 12th wireframe card (`l`) explaining Dual-Tier Chords, Regimes & Cell Footprints."""
    x_enc = alt.X("x:Q", scale=alt.Scale(domain=[0.0, 10.0], nice=False), axis=_hidden_canvas_axis())
    y_enc = alt.Y("y:Q", scale=alt.Scale(domain=[0.0, 10.0], nice=False), axis=_hidden_canvas_axis())

    layers: list[alt.Chart] = []

    # Section headers
    hdr_df = pl.DataFrame(
        [
            {"x": 0.45, "y": 9.35, "txt": "1. Dual-Tier Significance Weight"},
            {"x": 0.45, "y": 7.15, "txt": "2. Pairwise Regime & Polarity"},
            {"x": 0.45, "y": 3.55, "txt": "3. Six Broad Host Cell-Type Clouds"},
        ]
    )
    layers.append(
        alt.Chart(hdr_df)
        .mark_text(fontSize=8.4, fontWeight="bold", fill="#0F172A", align="left", baseline="middle")
        .encode(x=x_enc, y=y_enc, text="txt:N")
    )

    # 1. Dual-tier stroke weight swatches
    layers.append(
        alt.Chart(pl.DataFrame([{"x": 0.55, "y": 8.62, "x2": 1.85, "y2": 8.62}]))
        .mark_rule(color="#0F172A", strokeWidth=2.35, opacity=0.95)
        .encode(x=x_enc, y=y_enc, x2="x2:Q", y2="y2:Q")
    )
    layers.append(
        alt.Chart(pl.DataFrame([{"x": 0.55, "y": 7.98, "x2": 1.85, "y2": 7.98}]))
        .mark_rule(color="#475569", strokeWidth=1.05, opacity=0.52)
        .encode(x=x_enc, y=y_enc, x2="x2:Q", y2="y2:Q")
    )
    tier_txt_df = pl.DataFrame(
        [
            {"x": 2.10, "y": 8.62, "txt": "Verified ACAT FDR < 0.10 (22 pairs)"},
            {"x": 2.10, "y": 7.98, "txt": "Stratified Perm. p < 0.05 (*, 26 pairs)"},
        ]
    )
    layers.append(
        alt.Chart(tier_txt_df)
        .mark_text(fontSize=7.8, fill="#1E293B", align="left", baseline="middle")
        .encode(x=x_enc, y=y_enc, text="txt:N")
    )

    # 2. Pairwise interaction regimes + endpoint polarity glyphs
    reg_specs = (
        (PairwiseVisualRegime.BALANCE_RATIO, [6, 2, 2, 2], 6.45, "Log-Ratio Balance"),
        (PairwiseVisualRegime.GATING_DISCORDANT, [1, 0], 5.85, "+/- Mutual-Exclusion Gate"),
        (PairwiseVisualRegime.SYNERGY_PP, [1, 0], 5.25, "++ Co-Presence Synergy"),
        (PairwiseVisualRegime.CO_ABUNDANCE, [2, 2], 4.65, "Co-Abundance Coupling"),
    )
    reg_lbl_rows: list[dict[str, object]] = []
    for reg_e, dash_p, y_r, r_lbl in reg_specs:
        layers.append(
            alt.Chart(pl.DataFrame([{"x": 0.55, "y": y_r, "x2": 1.65, "y2": y_r}]))
            .mark_rule(color=regime_color_hex(reg_e), strokeWidth=1.9, strokeDash=dash_p)
            .encode(x=x_enc, y=y_enc, x2="x2:Q", y2="y2:Q")
        )
        reg_lbl_rows.append({"x": 1.85, "y": y_r, "txt": r_lbl})

    # Polar glyphs on right side of section 2
    layers.append(
        alt.Chart(pl.DataFrame([{"x": 6.45, "y": 6.45}]))
        .mark_circle(size=42, color=OKABE_BLUE, stroke="#FFFFFF", strokeWidth=0.9)
        .encode(x=x_enc, y=y_enc)
    )
    layers.append(
        alt.Chart(pl.DataFrame([{"x": 6.45, "y": 5.65}]))
        .mark_point(shape="circle", filled=True, size=48, fill="#FFFFFF", stroke=OKABE_VERMILION, strokeWidth=1.8)
        .encode(x=x_enc, y=y_enc)
    )
    layers.append(
        alt.Chart(pl.DataFrame([{"x": 6.45, "y": 4.85}]))
        .mark_point(shape="diamond", filled=True, size=42, fill=OKABE_ORANGE, stroke="#0F172A", strokeWidth=0.8)
        .encode(x=x_enc, y=y_enc)
    )
    reg_lbl_rows.extend(
        [
            {"x": 6.82, "y": 6.45, "txt": "+ Enriched"},
            {"x": 6.82, "y": 5.65, "txt": "- Gated"},
            {"x": 6.82, "y": 4.85, "txt": "-- Desert"},
        ]
    )
    layers.append(
        alt.Chart(pl.DataFrame(reg_lbl_rows))
        .mark_text(fontSize=7.5, fill="#1E293B", align="left", baseline="middle")
        .encode(x=x_enc, y=y_enc, text="txt:N")
    )

    # 3. Six broad host cell-type swatches (vivid dot = interacting state; soft tint = lineage context)
    comp_grid_coords = (
        (TmeCompartment.T_NK, 0.70, 2.75),
        (TmeCompartment.B_LINEAGE, 0.70, 1.95),
        (TmeCompartment.PLASMA, 0.70, 1.15),
        (TmeCompartment.MYELOID_MAST, 5.25, 2.75),
        (TmeCompartment.STROMAL_MESENCHYMAL, 5.25, 1.95),
        (TmeCompartment.ENDOTHELIAL_VASCULAR, 5.25, 1.15),
    )
    comp_lbl_rows: list[dict[str, object]] = []
    for cp_e, cx, cy in comp_grid_coords:
        disp_n, _sl, c_hex = _BROAD_HOST_COMPARTMENT_META[cp_e]
        layers.append(
            alt.Chart(pl.DataFrame([{"x": cx, "y": cy}]))
            .mark_circle(size=52, color=c_hex, stroke="#0F172A", strokeWidth=0.45)
            .encode(x=x_enc, y=y_enc)
        )
        layers.append(
            alt.Chart(pl.DataFrame([{"x": cx + 0.35, "y": cy}]))
            .mark_circle(size=42, color=_blend_hex_with_white(c_hex, 0.30))
            .encode(x=x_enc, y=y_enc)
        )
        comp_lbl_rows.append({"x": cx + 0.68, "y": cy, "txt": disp_n})

    comp_lbl_rows.append(
        {"x": 0.45, "y": 0.38, "txt": "Vivid dot: Sig. state kNN · Pale dot: Lineage context"}
    )
    layers.append(
        alt.Chart(pl.DataFrame(comp_lbl_rows))
        .mark_text(fontSize=7.4, fill="#334155", align="left", baseline="middle")
        .encode(x=x_enc, y=y_enc, text="txt:N")
    )

    return cast(
        alt.LayerChart,
        alt.layer(*layers).properties(
            width=width,
            height=height,
            title=alt.TitleParams(
                text="l. Architecture, Significance & Lineage Key",
                subtitle="4 empty compartment pairs omitted (0 FDR & 0 perm hits)",
                fontSize=10.3,
                fontWeight="bold",
                subtitleFontSize=8.0,
                subtitleColor="#475569",
                anchor="start",
            ),
        ),
    )


def build_pan_cohort_compartment_pair_umap_grid(
    cell_bundles: Sequence[CompartmentPairUmapCellBundle],
    fdr_threshold: float = 0.10,
    n_cols: int = 4,
    cell_width: int = 268,
    cell_height: int = 228,
) -> Result[AltairComposite, str]:
    """Assembles the 4x3 multi-panel UMAP grid across all active broad host cell-type pairs."""
    if not cell_bundles:
        return Failure("Cannot build compartment-pair UMAP grid from empty sequence")

    cols_eff = max(int(n_cols), 1)
    panel_letters = "abcdefghijklmnopqrstuvwxyz"
    rendered_panels: list[AltairComposite] = [
        build_compartment_pair_umap_panel(
            cell_bundle=cb,
            panel_prefix=panel_letters[idx] if idx < len(panel_letters) else f"p{idx + 1}",
            width=cell_width,
            height=cell_height,
            compact_grid_mode=True,
        )
        for idx, cb in enumerate(cell_bundles)
    ]
    rendered_panels.append(
        _build_compartment_grid_legend_card(width=cell_width, height=cell_height)
    )

    rows: list[AltairComposite] = []
    for start_i in range(0, len(rendered_panels), cols_eff):
        chunk = rendered_panels[start_i : start_i + cols_eff]
        row_c: AltairComposite = (
            chunk[0]
            if len(chunk) == 1
            else alt.hconcat(*chunk, spacing=18).resolve_scale(
                color="independent",
                strokeDash="independent",
                strokeWidth="independent",
                opacity="independent",
            )
        )
        rows.append(row_c)

    grid: AltairComposite = (
        rows[0]
        if len(rows) == 1
        else alt.vconcat(*rows, spacing=20).resolve_scale(
            color="independent",
            strokeDash="independent",
            strokeWidth="independent",
            opacity="independent",
        )
    )

    tot_fdr = sum(cb.spec.n_fdr_pairs for cb in cell_bundles)
    tot_perm = sum(cb.spec.n_perm_only_pairs for cb in cell_bundles)
    titled = grid.properties(
        title=alt.TitleParams(
            text="Pan-Cancer Pairwise Single-Cell UMAP Atlas Stratified by Broad Host TME Cell-Type Pairs",
            subtitle=(
                f"Unified 9-Cohort Harmony I2 Manifold (N = 613 Biopsies, 20 CMH Strata) | "
                f"{len(cell_bundles)} Active Broad Cell-Type Pairs: {tot_fdr} Verified ACAT FDR < {fdr_threshold:.2f} Chords (Bold) + "
                f"{tot_perm} Stratified Perm. p < 0.05 Chords (Subtle)"
            ),
            anchor="start",
            fontSize=13.5,
            fontWeight="bold",
            subtitleFontSize=10.0,
            subtitleColor="#334155",
            offset=12,
        )
    )
    themed = apply_nature_methods_theme(titled).configure_legend(disable=True)
    return Success(cast(AltairComposite, themed))


def build_standalone_compartment_pair_umap_composite(
    cell_bundle: CompartmentPairUmapCellBundle,
    pw_bundle: PairwiseUmapBundle,
    selected_regime: HarmonyRegime,
    n_biopsies: int,
    n_strata: int,
    fdr_threshold: float = 0.10,
) -> AltairComposite:
    """Assembles a standalone 2x2 (`a–d`) UMAP composite for one broad cell-type pair."""
    viz_cfg = UmapPairVizConfig(
        panel_width=305,
        panel_height=260,
        fdr_threshold=fdr_threshold,
    )
    panel_a = build_compartment_pair_umap_panel(
        cell_bundle=cell_bundle,
        panel_prefix="a",
        width=viz_cfg.panel_width,
        height=viz_cfg.panel_height,
        compact_grid_mode=False,
    )
    panel_b = build_denovo_hub_umap_panel(pw_bundle, config=viz_cfg, panel_prefix="b")
    panel_c = build_anchor_partner_umap_panel(pw_bundle, config=viz_cfg, panel_prefix="c")
    panel_d = build_lead_pair_spotlight_panel(pw_bundle, config=viz_cfg, panel_prefix="d")

    row_top = alt.hconcat(panel_a, panel_b, spacing=24).resolve_scale(
        color="independent",
        strokeDash="independent",
        strokeWidth="independent",
        opacity="independent",
    )
    row_bot = alt.hconcat(panel_c, panel_d, spacing=24).resolve_scale(
        color="independent",
        strokeDash="independent",
        strokeWidth="independent",
        opacity="independent",
    )

    spec = cell_bundle.spec
    comp_title = (
        f"Pan-Cancer Broad Lineage Interaction Atlas — {spec.display_title} ({spec.macro_category})"
    )
    comp_sub = (
        f"{harmony_regime_display_label(selected_regime)} | N = {n_biopsies} Biopsies ({n_strata} CMH Strata) | "
        f"{spec.n_fdr_pairs} Verified ACAT FDR < {fdr_threshold:.2f} (Bold) + "
        f"{spec.n_perm_only_pairs} Stratified Perm. p < 0.05 (Subtle)"
    )
    composite = (
        alt.vconcat(row_top, row_bot, spacing=24)
        .resolve_scale(
            color="independent",
            strokeDash="independent",
            strokeWidth="independent",
            opacity="independent",
        )
        .properties(
            title=alt.TitleParams(
                text=comp_title,
                subtitle=comp_sub,
                fontSize=13.5,
                fontWeight="bold",
                subtitleFontSize=10.2,
                subtitleColor="#334155",
                anchor="start",
                offset=12,
            )
        )
    )
    return cast(
        AltairComposite,
        apply_nature_methods_theme(composite).configure_legend(labelLimit=245, titleLimit=245),
    )


# =============================================================================
# 5. Pure Functional Core — Markdown Report & Pipeline Orchestration
# =============================================================================


def render_pan_cohort_markdown_report(
    bundle: PanCohortFeatureBundle,
    selected_regime: HarmonyRegime,
    harmony_bench_df: pl.DataFrame,
    power_ablation_df: pl.DataFrame,
    per_cohort_uplift_df: pl.DataFrame,
    meta9_df: pl.DataFrame,
    verified_pan_df: pl.DataFrame,
    concordance_df: pl.DataFrame,
    output_dir: Path,
) -> str:
    """Formats the comprehensive Markdown report summarizing integration and power results."""
    n_evals = len(bundle.eval_sample_uids)
    n_r = sum(1 for r in bundle.eval_response_labels if r == "responder")
    n_nr = n_evals - n_r
    n_strata = len(set(bundle.eval_substudy_strata))

    lines: list[str] = [
        "# Pan-Cohort `harmonypy` Integration & Pairwise Statistical Power Summary",
        "",
        "## 1. Unified Nine-Cohort Single-Cell Atlas",
        "",
        f"- **Clinical Cohorts Combined**: `9` single-cell ICB datasets (`{bundle.n_total_cells_all_cohorts:,}` total cells, "
        f"`{int(bundle.anchor_x_pca.shape[0]):,}` donor-stratified anchor cells in PCA/Harmony graph, **100% cell projection** onto pan-cohort Milo neighborhoods).",
        f"- **Biopsy Evaluations**: **`{n_evals}`** (`{n_r}` Responders / `{n_nr}` Non-Responders) stratified across **`{n_strata}`** `cohort::substudy` CMH strata.",
        f"- **Feature Alignment Mode**: `{bundle.load_mode.value}` (strictly outcome-blind; `harmonypy` receives anonymized `sample_uid` without response labels).",
        f"- **Selected Pareto-Optimal Latent Space**: **`{harmony_regime_display_label(selected_regime)}`**.",
        "",
        "## 2. `harmonypy` Latent Integration Benchmark (`I0` -> `I3`)",
        "",
        "| Latent Space | Cohort iLISI | Sub-Study iLISI | Donor iLISI | Compartment cLISI | Meta42 cLISI | kBET Accept (%) | Median Cohorts / Nhood | Cohort Entropy (bits) | Mixing-QC Pass (%) |",
        "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
    ]
    for r in harmony_bench_df.to_dicts():
        lines.append(
            f"| {r['regime_label']} | **{float(r['cohort_ilisi']):.2f}** | {float(r['substudy_ilisi']):.2f} | "
            f"{float(r['donor_ilisi']):.2f} | **{float(r['compartment_clisi']):.2f}** | {float(r['meta_state_clisi']):.2f} | "
            f"{100.0 * float(r['kbet_acceptance_rate']):.1f}% | **{float(r['median_cohorts_per_nhood']):.0f}** | "
            f"{float(r['median_cohort_entropy_bits']):.2f} | {float(r['passes_mixing_qc_pct']):.1f}% |"
        )

    lines.extend(
        [
            "",
            "## 3. Seven-Tier Statistical Power Progression vs. Permuted-Response Null Control",
            "",
            "| Power Tier | Real Cross-Comp State Pairs (`q<0.10`) | Real Super-Nhood Hits | Null Mean State Pairs (`R` Strata Perms) | Null Max State Pairs | Empirical FDR Calibration |",
            "| :--- | :---: | :---: | :---: | :---: | :---: |",
        ]
    )
    for r in power_ablation_df.to_dicts():
        real_k = int(r["real_sig_state_pairs"])
        null_m = float(r["null_mean_sig_state_pairs"])
        emp_fdr = null_m / float(max(real_k, 1))
        lines.append(
            f"| {r['tier_label']} | **{real_k}** | {int(r['real_sig_nhood_pairs'])} | "
            f"**{null_m:.2f}** | {int(r['null_max_sig_state_pairs'])} | `{emp_fdr:.3f}` |"
        )

    lines.extend(
        [
            "",
            "## 4. Per-Cohort Discovery Uplift Across All Nine Datasets",
            "",
            "| Accession | Clinical Cohort | Biopsies (`J`) | Assay Scope | `P0` Baseline (`q<0.10`) | `P1` (`+EB+Depth+B64`) | `P2` (`+ACAT+Prior`) | Active in Pan-Cohort `P3/P6` |",
            "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |",
        ]
    )
    for r in per_cohort_uplift_df.sort("n_evals", descending=True).to_dicts():
        lines.append(
            f"| `{r['accession']}` | {r['cohort_display']} | {int(r['n_evals'])} | `{r['assay_type']}` | "
            f"{int(r['p0_baseline_pairs'])} | {int(r['p1_eb_depth_b64_pairs'])} | **{int(r['p2_boosted_pairs'])}** | "
            f"**{int(r['active_pan_cohort_pairs'])}** |"
        )

    if not concordance_df.is_empty() and concordance_df.height >= 3:
        xz = concordance_df["z_one_stage_harmony"].to_numpy()
        yz = concordance_df["z_two_stage_meta"].to_numpy()
        pr = float(pearsonr(xz, yz).statistic) if float(np.std(xz)) > 1e-9 and float(np.std(yz)) > 1e-9 else 0.0
        sr = float(spearmanr(xz, yz).statistic) if float(np.std(xz)) > 1e-9 and float(np.std(yz)) > 1e-9 else 0.0
        lines.extend(
            [
                "",
                "## 5. One-Stage `harmonypy` vs. Two-Stage Random-Effects Meta-IUT Concordance",
                "",
                f"- Across `{concordance_df.height}` testable cross-compartment `Meta42` state pairs: "
                f"**Pearson `r = {pr:.3f}`**, **Spearman `rho = {sr:.3f}`**.",
                f"- Two-Stage DerSimonian-Laird Meta-IUT (`P3`) significant cross-compartment state-regimes (`q < 0.10`, `LOCO > 0`): "
                f"**`{int(meta9_df.filter(pl.col('is_significant')).height) if not meta9_df.is_empty() else 0}`**.",
            ]
        )

    if not verified_pan_df.is_empty():
        sig_v = verified_pan_df.filter(pl.col("is_significant")).head(15)
        if not sig_v.is_empty():
            lines.extend(
                [
                    "",
                    "## 6. Top Verified Pan-Cancer Cross-Compartment Interactions (`P6` + `LOCO`)",
                    "",
                    "| Pair (`State A x State B`) | Compartments | Dominant Regime | Direction | Signed `Z_ACAT` | `FDR_ACAT` | Perm `p` | Worst-Case `min_LOCO_Z` |",
                    "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |",
                ]
            )
            for r in sig_v.to_dicts():
                lines.append(
                    f"| `{r['state_a']}` x `{r['state_b']}` ({str(r['state_a_readable'])[:20]} / {str(r['state_b_readable'])[:20]}) | "
                    f"`{r['compartment_a']}` <-> `{r['compartment_b']}` | `{r['dominant_regime']}` | `{r['dominant_direction']}` | "
                    f"**{float(r['z_acat_signed']):+.2f}** | `{float(r['fdr_acat_iut']):.3g}` | `{float(r['emp_perm_pvalue']):.3g}` | "
                    f"**{float(r['min_loco_z_acat']):.2f}** |"
                )

    lines.extend(
        [
            "",
            "## 7. Exported Deliverables",
            "",
            f"- Integration Audit Figure: `{output_dir / 'figures' / 'pan_cohort_harmony_integration_audit.svg'}` (`.png`)",
            f"- Power Ablation & Null Calibration Figure: `{output_dir / 'figures' / 'statistical_power_ablation_and_null_calibration.svg'}` (`.png`)",
            f"- Unified 9-Cohort Pairwise UMAP Composite: `{output_dir / 'figures' / 'pan_cohort_unified_pairwise_umap_composite.svg'}` (`.png`)",
            f"- Cross-Compartment Z_ACAT Heatmap (Full 26-State Atlas): `{output_dir / 'figures' / 'pan_cohort_cross_compartment_z_acat_heatmap.svg'}` (`.png`)",
            f"- Simplified Core Hub Z_ACAT Heatmap (13 Key States): `{output_dir / 'figures' / 'pan_cohort_core_hubs_z_acat_heatmap.svg'}` (`.png`)",
            f"- Broad Cell-Type Pair UMAP Atlas Grid (4x3): `{output_dir / 'figures' / 'pan_cohort_compartment_pair_umap_grid.svg'}` (`.png`)",
            f"- Standalone Broad Cell-Type Pair 2x2 UMAP Composites: `{output_dir / 'figures' / 'compartment_pair_umaps' / 'pan_cohort_umap_<pair>.svg'}` (`.png`)",
            "",
        ]
    )
    return "\n".join(lines)


def export_all_compartment_pair_umap_figures(
    bundle: PanCohortFeatureBundle,
    pan_unit: PreparedSuperUnit,
    umap_coords: FloatArray,
    selected_regime: HarmonyRegime,
    verified_pan_df: pl.DataFrame,
    fig_dir: Path,
    exporter: VectorFigureExporter,
    fdr_threshold: float = 0.10,
    seed: int = 42,
) -> Result[tuple[Path, tuple[Path, ...]], str]:
    """Exports the 4x3 broad cell-type pair UMAP grid + all standalone 2x2 pair composites."""
    specs = enumerate_active_broad_compartment_pairs(
        verified_pan_df=verified_pan_df,
        fdr_threshold=fdr_threshold,
        perm_alpha=0.05,
    )
    if not specs:
        return Failure("No active host broad compartment pairs found for UMAP export")

    bundle_cache_root = fig_dir.parent / "compartment_pair_umap_bundles"
    grid_cell_bundles: list[CompartmentPairUmapCellBundle] = []
    for idx_s, spec in enumerate(specs):
        cb_grid_res = build_compartment_pair_umap_cell_bundle(
            spec=spec,
            pan_unit=pan_unit,
            umap_coords=umap_coords,
            anchor_compartments=bundle.anchor_compartments,
            verified_pan_df=verified_pan_df,
            fdr_threshold=fdr_threshold,
            perm_alpha=0.05,
            max_cells=8200,
            max_labels=8,
            compact_labels=True,
            seed=seed + idx_s,
        )
        match cb_grid_res:
            case Failure(err):
                return Failure(err)
            case Success(cb_grid):
                grid_cell_bundles.append(cb_grid)
                slug_dir = bundle_cache_root / spec.slug
                slug_dir.mkdir(parents=True, exist_ok=True)
                cb_grid.cells_df.write_parquet(slug_dir / "cells.parquet")
                cb_grid.tier1_arcs_df.write_parquet(slug_dir / "tier1_arcs.parquet")
                cb_grid.tier2_arcs_df.write_parquet(slug_dir / "tier2_arcs.parquet")
                cb_grid.endpoints_df.write_parquet(slug_dir / "endpoints.parquet")
                cb_grid.medoid_labels_df.write_parquet(slug_dir / "medoid_labels.parquet")
                pl.DataFrame(
                    [
                        {
                            "slug": spec.slug,
                            "display_first": spec.display_first,
                            "display_second": spec.display_second,
                            "display_title": spec.display_title,
                            "macro_category": spec.macro_category,
                            "color_first": spec.color_first,
                            "color_second": spec.color_second,
                            "n_fdr_pairs": spec.n_fdr_pairs,
                            "n_perm_only_pairs": spec.n_perm_only_pairs,
                            "umap_xmin": cb_grid.umap_bounds[0],
                            "umap_xmax": cb_grid.umap_bounds[1],
                            "umap_ymin": cb_grid.umap_bounds[2],
                            "umap_ymax": cb_grid.umap_bounds[3],
                            "lead_pair_summary": cb_grid.lead_pair_summary,
                        }
                    ]
                ).write_parquet(slug_dir / "meta.parquet")

    grid_chart_res = build_pan_cohort_compartment_pair_umap_grid(
        cell_bundles=grid_cell_bundles,
        fdr_threshold=fdr_threshold,
    )
    match grid_chart_res:
        case Failure(err):
            return Failure(err)
        case Success(grid_chart):
            pass

    grid_exp = exporter(
        grid_chart,
        fig_dir / "pan_cohort_compartment_pair_umap_grid",
        "Pan-Cancer Broad Host TME Cell-Type Pair UMAP Atlas",
        "02_Compartment_Pair_UMAP_Grid",
    )
    match grid_exp:
        case Failure(err):
            return Failure(err)
        case Success((_grid_svg, grid_png)):
            pass

    anchor_read_lbls = tuple(
        bundle.code_to_readable.get(mc, mc) for mc in bundle.anchor_meta_codes
    )
    n_biopsies = len(bundle.eval_sample_uids)
    n_strata = len(set(bundle.eval_substudy_strata))
    sub_dir = fig_dir / "compartment_pair_umaps"
    standalone_pngs: list[Path] = []

    for idx_s, spec in enumerate(specs):
        cb_full_res = build_compartment_pair_umap_cell_bundle(
            spec=spec,
            pan_unit=pan_unit,
            umap_coords=umap_coords,
            anchor_compartments=bundle.anchor_compartments,
            verified_pan_df=verified_pan_df,
            fdr_threshold=fdr_threshold,
            perm_alpha=0.05,
            max_cells=11000,
            max_labels=12,
            compact_labels=False,
            seed=seed + 100 + idx_s,
        )
        match cb_full_res:
            case Failure(err):
                return Failure(err)
            case Success(cb_full):
                pass

        pw_sub_res = build_compartment_pair_pairwise_umap_bundle(
            spec=spec,
            pan_unit=pan_unit,
            umap_coords=umap_coords,
            anchor_cell_labels=anchor_read_lbls,
            global_response_labels=bundle.eval_response_labels,
            verified_pan_df=verified_pan_df,
            fdr_threshold=fdr_threshold,
            perm_alpha=0.05,
        )
        match pw_sub_res:
            case Failure(err):
                return Failure(err)
            case Success(pw_sub):
                pass

        comp_fig = build_standalone_compartment_pair_umap_composite(
            cell_bundle=cb_full,
            pw_bundle=pw_sub,
            selected_regime=selected_regime,
            n_biopsies=n_biopsies,
            n_strata=n_strata,
            fdr_threshold=fdr_threshold,
        )
        sub_exp = exporter(
            comp_fig,
            sub_dir / f"pan_cohort_umap_{spec.slug}",
            f"Pan-Cancer Pairwise UMAP Composite — {spec.display_title}",
            f"02_UMAP_{spec.slug}",
        )
        match sub_exp:
            case Failure(err):
                return Failure(err)
            case Success((_sub_svg, sub_png)):
                standalone_pngs.append(sub_png)

    return Success((grid_png, tuple(standalone_pngs)))


def run_pan_cohort_pipeline(
    cfg: PanCohortCliConfig,
    loader: PanCohortFeatureLoader,
    exporter: VectorFigureExporter,
    iatlas_df: pl.DataFrame,
) -> Result[PanCohortPipelineOutputs, str]:
    """Executes the full 9-cohort Harmony integration, power benchmark, and figure export."""
    load_res = loader(cfg)
    match load_res:
        case Failure(err):
            return Failure(err)
        case Success(bundle):
            pass

    # 1. Benchmark all 4 latent integration regimes (I0..I3)
    embeddings: dict[HarmonyRegime, FloatArray] = {}
    metrics_map: dict[HarmonyRegime, EmbeddingIntegrationMetrics] = {}
    atlases_map: dict[HarmonyRegime, PanCohortMiloAtlas] = {}
    bench_rows: list[dict[str, object]] = []

    short_tags = {
        HarmonyRegime.I0_UNINTEGRATED_PCA: "I0 (Raw PCA)",
        HarmonyRegime.I1_HARMONY_COHORT: "I1 (Cohort)",
        HarmonyRegime.I2_HARMONY_COHORT_SUBSTUDY: "I2 (Cohort+SubStudy)",
        HarmonyRegime.I3_HARMONY_COHORT_SUBSTUDY_DONOR: "I3 (+Donor)",
    }

    for reg in ALL_HARMONY_REGIMES:
        t0_r = time.perf_counter()
        max_iter_h = 4 if reg == HarmonyRegime.I3_HARMONY_COHORT_SUBSTUDY_DONOR else 12
        emb_res = run_harmony_embedding_regime(
            x_pca=bundle.anchor_x_pca,
            cohort_ids=bundle.anchor_cohort_ids,
            substudy_ids=bundle.anchor_substudy_ids,
            sample_uids=bundle.anchor_sample_uids,
            regime=reg,
            nclust=36,
            sigma=0.22,
            max_iter_harmony=max_iter_h,
            random_state=cfg.seed,
        )
        match emb_res:
            case Failure(err):
                return Failure(err)
            case Success(emb_mat):
                embeddings[reg] = emb_mat

        met_res = evaluate_embedding_integration_metrics(
            embedding=emb_mat,
            cohort_ids=bundle.anchor_cohort_ids,
            substudy_ids=bundle.anchor_substudy_ids,
            sample_uids=bundle.anchor_sample_uids,
            compartments=bundle.anchor_compartments,
            meta_states=bundle.anchor_meta_codes,
            regime=reg,
            knn_k=cfg.knn_k,
            seed=cfg.seed,
        )
        match met_res:
            case Failure(err):
                return Failure(err)
            case Success(met):
                metrics_map[reg] = met

        atl_res = build_pan_cohort_milo_atlas(
            embedding=emb_mat,
            bundle=bundle,
            regime=reg,
            knn_k=cfg.knn_k,
            milo_prop=cfg.milo_prop,
            seed=cfg.seed,
        )
        match atl_res:
            case Failure(err):
                return Failure(err)
            case Success(atl):
                atlases_map[reg] = atl

        dt_r = time.perf_counter() - t0_r
        print(
            f"[Harmony] {reg.value}: cohort_iLISI={met.cohort_ilisi:.2f}, "
            f"substudy_iLISI={met.substudy_ilisi:.2f}, donor_iLISI={met.donor_ilisi:.2f}, "
            f"comp_cLISI={met.compartment_clisi:.2f}, kBET={100.0 * met.kbet_acceptance_rate:.1f}%, "
            f"nhood_cohorts_med={atl.median_cohorts_per_nhood:.0f} ({dt_r:.1f}s)"
        )
        bench_rows.append(
            {
                "regime": reg.value,
                "regime_label": harmony_regime_display_label(reg),
                "short_tag": short_tags[reg],
                "cohort_ilisi": met.cohort_ilisi,
                "substudy_ilisi": met.substudy_ilisi,
                "donor_ilisi": met.donor_ilisi,
                "compartment_clisi": met.compartment_clisi,
                "meta_state_clisi": met.meta_state_clisi,
                "batch_asw": met.batch_asw,
                "compartment_asw": met.compartment_asw,
                "kbet_acceptance_rate": met.kbet_acceptance_rate,
                "n_pan_nhoods": int(atl.counts_full_pan.shape[0]),
                "median_cohorts_per_nhood": atl.median_cohorts_per_nhood,
                "median_cohort_ilisi_nhood": atl.median_cohort_ilisi_nhood,
                "median_cohort_entropy_bits": atl.median_cohort_entropy_bits,
                "median_cancer_entropy_bits": atl.median_cancer_entropy_bits,
                "median_donor_lisi": atl.median_donor_lisi,
                "passes_mixing_qc_pct": atl.passes_mixing_qc_pct,
            }
        )

    selected_regime = select_pareto_optimal_harmony_regime(metrics_map, atlases_map)
    harmony_bench_df = pl.DataFrame(bench_rows).with_columns(
        (pl.col("regime") == selected_regime.value).alias("is_selected_pareto_optimum")
    )
    pan_nhood_mixing_df = pl.concat([atlases_map[r].nhood_mixing_df for r in ALL_HARMONY_REGIMES])
    print(f"[Harmony] Selected Pareto-optimal embedding: {selected_regime.value}")

    # 2. Prepare within-cohort and pan-cohort coarse super-neighborhood units (strictly outcome-blind!)
    cohort_units: list[PreparedSuperUnit] = []
    cohort_eval_indices: dict[str, IntArray] = {}
    eval_coh_arr = np.array(bundle.eval_cohort_ids, dtype=str)

    for coh_slice in bundle.cohorts:
        c_global_cols = np.flatnonzero(eval_coh_arr == coh_slice.accession)
        cohort_eval_indices[coh_slice.accession] = c_global_cols
        cu_res = prepare_super_nhood_unit(
            name=coh_slice.accession,
            nhoods_mat=coh_slice.nhoods_local,
            count_matrix=coh_slice.count_matrix_local,
            library_sizes=coh_slice.library_sizes,
            cohort_labels=tuple([coh_slice.accession] * len(coh_slice.eval_uids)),
            strata_labels=coh_slice.substudy_ids,
            nhood_meta_codes=coh_slice.nhood_meta_codes,
            code_to_readable=bundle.code_to_readable,
            cell_sample_codes=Some(coh_slice.cell_local_eval),
            resolution=0.25,
            min_jaccard=0.02,
            seed=cfg.seed,
        )
        match cu_res:
            case Failure(err):
                return Failure(err)
            case Success(c_unit):
                cohort_units.append(c_unit)

    sel_atlas = atlases_map[selected_regime]
    pan_unit_res = prepare_super_nhood_unit(
        name="PAN_COHORT_HARMONY9",
        nhoods_mat=sel_atlas.nhoods_anchor,
        count_matrix=sel_atlas.counts_full_pan,
        library_sizes=bundle.eval_library_sizes,
        cohort_labels=bundle.eval_cohort_ids,
        strata_labels=bundle.eval_substudy_strata,
        nhood_meta_codes=sel_atlas.nhood_meta_codes,
        code_to_readable=bundle.code_to_readable,
        anchor_eval_counts=Some(sel_atlas.anchor_eval_counts),
        resolution=cfg.super_resolution,
        min_jaccard=cfg.super_min_jaccard,
        seed=cfg.seed,
    )
    match pan_unit_res:
        case Failure(err):
            return Failure(err)
        case Success(pan_unit):
            pass

    print(
        f"[Super-Nhoods] Pan-cohort {sel_atlas.counts_full_pan.shape[0]} Milo nhoods -> "
        f"{pan_unit.counts_super.shape[0]} consensus super-neighborhoods across J={len(bundle.eval_sample_uids)} biopsies"
    )

    # 3. Evaluate 7-Tier Power Progression on Real Clinical Outcomes & Within-Stratum Permuted Null Controls
    state_prior_weights = build_pan_iatlas_state_pair_weights(iatlas_df)
    real_res = evaluate_seven_power_tiers_single_pass(
        cohort_units=cohort_units,
        cohort_eval_indices=cohort_eval_indices,
        pan_unit=pan_unit,
        global_response_labels=bundle.eval_response_labels,
        state_prior_weights=state_prior_weights,
        cfg=cfg,
    )
    match real_res:
        case Failure(err):
            return Failure(err)
        case Success(real_pass):
            pass

    null_state_runs: dict[PowerTier, list[int]] = {t: [] for t in ALL_POWER_TIERS}
    null_nhood_runs: dict[PowerTier, list[int]] = {t: [] for t in ALL_POWER_TIERS}
    for rep_i in range(cfg.n_null_calibration_replicates):
        perm_y = permute_outcomes_within_strata(
            response_labels=bundle.eval_response_labels,
            strata_labels=bundle.eval_substudy_strata,
            seed=202600 + 37 * rep_i,
        )
        null_res = evaluate_seven_power_tiers_single_pass(
            cohort_units=cohort_units,
            cohort_eval_indices=cohort_eval_indices,
            pan_unit=pan_unit,
            global_response_labels=perm_y,
            state_prior_weights=state_prior_weights,
            cfg=cfg,
        )
        match null_res:
            case Failure(err):
                return Failure(err)
            case Success(null_pass):
                for t in ALL_POWER_TIERS:
                    null_state_runs[t].append(null_pass.tier_state_pair_counts[t])
                    null_nhood_runs[t].append(null_pass.tier_nhood_pair_counts[t])

    power_rows: list[dict[str, object]] = []
    for t in ALL_POWER_TIERS:
        ns_arr = np.array(null_state_runs[t], dtype=np.float64)
        nn_arr = np.array(null_nhood_runs[t], dtype=np.float64)
        real_s = real_pass.tier_state_pair_counts[t]
        real_n = real_pass.tier_nhood_pair_counts[t]
        print(
            f"[Power] {t.value}: Real={real_s} cross-comp state pairs ({real_n} super-nhood hits) | "
            f"Permuted-Y Null mean={float(np.mean(ns_arr)):.2f} (max={int(np.max(ns_arr))})"
        )
        power_rows.append(
            {
                "tier": t.value,
                "tier_label": power_tier_display_label(t),
                "real_sig_state_pairs": real_s,
                "real_sig_nhood_pairs": real_n,
                "null_mean_sig_state_pairs": float(np.mean(ns_arr)),
                "null_sd_sig_state_pairs": float(np.std(ns_arr)),
                "null_max_sig_state_pairs": int(np.max(ns_arr)),
                "null_mean_sig_nhood_pairs": float(np.mean(nn_arr)),
                "n_null_replicates": cfg.n_null_calibration_replicates,
            }
        )
    power_ablation_df = pl.DataFrame(power_rows)

    # 4. Run LOCO 9-cohort stability audit & stratified permutation verification on P6 discoveries
    loco_res = compute_pan_cohort_loco_and_verified_pairs(
        pan_unit=pan_unit,
        real_pass=real_pass,
        global_response_labels=bundle.eval_response_labels,
        state_prior_weights=state_prior_weights,
        code_to_readable=bundle.code_to_readable,
        cfg=cfg,
    )
    match loco_res:
        case Failure(err):
            return Failure(err)
        case Success((verified_pan_df, loco_forest_df)):
            pass

    concordance_df = build_one_vs_two_stage_concordance_df(
        verified_pan_df=verified_pan_df,
        meta9_df=real_pass.meta9_df,
        fdr_threshold=cfg.fdr_threshold,
    )

    # Compute per-cohort uplift table including active participation in P3 / P6 cross-compartment pairs
    pres_pan_check = build_compartment_observability_mask(
        presence=estimate_nhood_presence(
            pan_unit.counts_super, Some(pan_unit.library_sizes)
        ).unwrap(),
        nhood_labels=pan_unit.super_meta_codes,
        cohort_labels=pan_unit.cohort_labels,
        min_cohort_compartment_share=0.015,
    ).unwrap()

    sig_p6_uv: list[tuple[int, int]] = (
        [
            (int(u), int(v))
            for u, v in verified_pan_df.filter(pl.col("is_significant"))
            .select(["nhood_u", "nhood_v"])
            .iter_rows()
        ]
        if not verified_pan_df.is_empty()
        else []
    )

    uplift_rows: list[dict[str, object]] = []
    for coh_slice in bundle.cohorts:
        acc = coh_slice.accession
        c_cols = cohort_eval_indices[acc]
        # Cohort is CD45+-sorted if non-immune compartments are masked (observability == 0) across the cohort
        obs_sub = pres_pan_check.observability[:, c_cols]
        is_cd45_sorted = bool(np.any(np.all(obs_sub <= 1e-9, axis=1)))
        active_pan = 0
        for u_id, v_id in sig_p6_uv:
            u_obs = bool(np.any(pres_pan_check.raw_binary_presence[u_id, c_cols]))
            v_obs = bool(np.any(pres_pan_check.raw_binary_presence[v_id, c_cols]))
            both_elig = bool(
                np.any(pres_pan_check.observability[u_id, c_cols] > 1e-6)
                and np.any(pres_pan_check.observability[v_id, c_cols] > 1e-6)
            )
            if both_elig and u_obs and v_obs:
                active_pan += 1

        uplift_rows.append(
            {
                "accession": acc,
                "cohort_display": coh_slice.cohort_display,
                "cohort_short": f"{acc} ({coh_slice.cancer_type}, J={len(coh_slice.eval_uids)})",
                "cancer_type": coh_slice.cancer_type,
                "n_evals": len(coh_slice.eval_uids),
                "assay_type": "CD45+ Immune FACS" if is_cd45_sorted else "Unsorted Whole-Tumor",
                "p0_baseline_pairs": real_pass.per_cohort_p0_counts.get(acc, 0),
                "p1_eb_depth_b64_pairs": real_pass.per_cohort_p1_counts.get(acc, 0),
                "p2_boosted_pairs": real_pass.per_cohort_p2_counts.get(acc, 0),
                "active_pan_cohort_pairs": active_pan,
            }
        )
    per_cohort_uplift_df = pl.DataFrame(uplift_rows)

    # 5. Compute 2D UMAPs (`I0` Unintegrated & Selected Harmony) and export all 3 publication composites
    print("[UMAP] Computing 2D UMAP projections for I0 and selected Harmony space...")
    umap_i0 = compute_pan_umap_2d(
        embeddings[HarmonyRegime.I0_UNINTEGRATED_PCA], seed=cfg.seed
    )
    umap_sel = compute_pan_umap_2d(embeddings[selected_regime], seed=cfg.seed)

    fig_dir = cfg.output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    fig1 = build_harmony_integration_audit_figure(
        bundle=bundle,
        umap_i0=umap_i0,
        umap_sel=umap_sel,
        selected_regime=selected_regime,
        harmony_bench_df=harmony_bench_df,
        atlas_i0=atlases_map[HarmonyRegime.I0_UNINTEGRATED_PCA],
        atlas_sel=sel_atlas,
        seed=cfg.seed,
    )
    exp1 = exporter(
        fig1,
        fig_dir / "pan_cohort_harmony_integration_audit",
        "Pan-Cohort Harmony Integration Benchmark Across Nine Clinical ICB Datasets",
        "02_Harmony_Integration_Panels",
    )
    match exp1:
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    fig2 = build_power_ablation_and_calibration_figure(
        power_ablation_df=power_ablation_df,
        per_cohort_uplift_df=per_cohort_uplift_df,
        concordance_df=concordance_df,
        loco_forest_df=loco_forest_df,
        fdr_threshold=cfg.fdr_threshold,
    )
    exp2 = exporter(
        fig2,
        fig_dir / "statistical_power_ablation_and_null_calibration",
        "Statistical Power Ablation, Null Calibration and LOCO Stability",
        "02_Power_Ablation_Panels",
    )
    match exp2:
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    anchor_read_lbls = tuple(
        bundle.code_to_readable.get(mc, mc) for mc in bundle.anchor_meta_codes
    )
    umap_bundle_res = build_pan_cohort_umap_bundle_from_p6(
        pan_unit=pan_unit,
        umap_coords=umap_sel,
        anchor_cell_labels=anchor_read_lbls,
        global_response_labels=bundle.eval_response_labels,
        real_pass=real_pass,
        fdr_threshold=cfg.fdr_threshold,
        max_chords=28,
    )
    match umap_bundle_res:
        case Failure(err):
            return Failure(err)
        case Success(pan_umap_bundle):
            pass

    n_p6_sig = real_pass.tier_state_pair_counts[PowerTier.P6_PAN_COHORT_HARMONY_FULL_POWER]
    fig3_title = "Unified Nine-Cohort Harmony Atlas — Pan-Cancer Pairwise Neighborhood Architecture"
    fig3_sub = (
        f"{harmony_regime_display_label(selected_regime)} | N = {len(bundle.eval_sample_uids)} Biopsies "
        f"(304R / 309NR, {len(set(bundle.eval_substudy_strata))} CMH Strata) | "
        f"{n_p6_sig} Verified Cross-Compartment Meta42 Pairs (Omnibus ACAT FDR < {cfg.fdr_threshold:.2f})"
    )
    fig3 = build_pairwise_umap_composite(
        bundle=pan_umap_bundle,
        title=fig3_title,
        subtitle=fig3_sub,
        config=UmapPairVizConfig(
            max_chords=28,
            max_bg_cells=12000,
            fdr_threshold=cfg.fdr_threshold,
        ),
    )
    exp3 = exporter(
        fig3,
        fig_dir / "pan_cohort_unified_pairwise_umap_composite",
        fig3_title,
        "02_Pan_Cohort_Pairwise_UMAP_Panels",
    )
    match exp3:
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    fig4 = build_pan_cohort_z_acat_heatmap_figure(
        verified_pan_df=verified_pan_df,
        fdr_threshold=cfg.fdr_threshold,
    )
    exp4 = exporter(
        fig4,
        fig_dir / "pan_cohort_cross_compartment_z_acat_heatmap",
        "Pan-Cancer Cross-Compartment Cell-State Interaction Heatmap",
        "02_Cross_Compartment_Z_ACAT_Heatmap",
    )
    match exp4:
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    fig5 = build_pan_cohort_core_hubs_z_acat_heatmap_figure(
        verified_pan_df=verified_pan_df,
        fdr_threshold=cfg.fdr_threshold,
    )
    exp5 = exporter(
        fig5,
        fig_dir / "pan_cohort_core_hubs_z_acat_heatmap",
        "Simplified Pan-Cancer Core TME Hub Interaction Heatmap",
        "02_Core_Hubs_Z_ACAT_Heatmap",
    )
    match exp5:
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    exp_comp_umaps = export_all_compartment_pair_umap_figures(
        bundle=bundle,
        pan_unit=pan_unit,
        umap_coords=umap_sel,
        selected_regime=selected_regime,
        verified_pan_df=verified_pan_df,
        fig_dir=fig_dir,
        exporter=exporter,
        fdr_threshold=cfg.fdr_threshold,
        seed=cfg.seed,
    )
    match exp_comp_umaps:
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    md_report = render_pan_cohort_markdown_report(
        bundle=bundle,
        selected_regime=selected_regime,
        harmony_bench_df=harmony_bench_df,
        power_ablation_df=power_ablation_df,
        per_cohort_uplift_df=per_cohort_uplift_df,
        meta9_df=real_pass.meta9_df,
        verified_pan_df=verified_pan_df,
        concordance_df=concordance_df,
        output_dir=cfg.output_dir,
    )

    return Success(
        PanCohortPipelineOutputs(
            selected_regime=selected_regime,
            harmony_benchmark_df=harmony_bench_df,
            pan_nhood_mixing_df=pan_nhood_mixing_df,
            power_ablation_df=power_ablation_df,
            per_cohort_uplift_df=per_cohort_uplift_df,
            two_stage_meta9_df=real_pass.meta9_df,
            concordance_1s_vs_2s_df=concordance_df,
            pan_verified_pairs_df=verified_pan_df,
            loco_forest_df=loco_forest_df,
            markdown_report=md_report,
        )
    )


def run_refresh_verified_and_heatmap(
    cfg: PanCohortCliConfig,
    loader: PanCohortFeatureLoader,
    exporter: VectorFigureExporter,
    iatlas_df: pl.DataFrame,
    cached_tables: CachedAuditTables,
) -> Result[PanCohortPipelineOutputs, str]:
    """Recomputes exact stratified permutations, LOCO stability, and exports the Z_ACAT heatmap."""
    load_res = loader(cfg)
    match load_res:
        case Failure(err):
            return Failure(err)
        case Success(bundle):
            pass

    selected_regime = cached_tables.selected_regime
    max_iter_h = 4 if selected_regime == HarmonyRegime.I3_HARMONY_COHORT_SUBSTUDY_DONOR else 12
    emb_res = run_harmony_embedding_regime(
        x_pca=bundle.anchor_x_pca,
        cohort_ids=bundle.anchor_cohort_ids,
        substudy_ids=bundle.anchor_substudy_ids,
        sample_uids=bundle.anchor_sample_uids,
        regime=selected_regime,
        nclust=36,
        sigma=0.22,
        max_iter_harmony=max_iter_h,
        random_state=cfg.seed,
    )
    match emb_res:
        case Failure(err):
            return Failure(err)
        case Success(emb_sel):
            pass

    atl_res = build_pan_cohort_milo_atlas(
        embedding=emb_sel,
        bundle=bundle,
        regime=selected_regime,
        knn_k=cfg.knn_k,
        milo_prop=cfg.milo_prop,
        seed=cfg.seed,
    )
    match atl_res:
        case Failure(err):
            return Failure(err)
        case Success(sel_atlas):
            pass

    cohort_units: list[PreparedSuperUnit] = []
    cohort_eval_indices: dict[str, IntArray] = {}
    eval_coh_arr = np.array(bundle.eval_cohort_ids, dtype=str)

    for coh_slice in bundle.cohorts:
        c_global_cols = np.flatnonzero(eval_coh_arr == coh_slice.accession)
        cohort_eval_indices[coh_slice.accession] = c_global_cols
        cu_res = prepare_super_nhood_unit(
            name=coh_slice.accession,
            nhoods_mat=coh_slice.nhoods_local,
            count_matrix=coh_slice.count_matrix_local,
            library_sizes=coh_slice.library_sizes,
            cohort_labels=tuple([coh_slice.accession] * len(coh_slice.eval_uids)),
            strata_labels=coh_slice.substudy_ids,
            nhood_meta_codes=coh_slice.nhood_meta_codes,
            code_to_readable=bundle.code_to_readable,
            cell_sample_codes=Some(coh_slice.cell_local_eval),
            resolution=0.25,
            min_jaccard=0.02,
            seed=cfg.seed,
        )
        match cu_res:
            case Failure(err):
                return Failure(err)
            case Success(c_unit):
                cohort_units.append(c_unit)

    pan_unit_res = prepare_super_nhood_unit(
        name="PAN_COHORT_HARMONY9",
        nhoods_mat=sel_atlas.nhoods_anchor,
        count_matrix=sel_atlas.counts_full_pan,
        library_sizes=bundle.eval_library_sizes,
        cohort_labels=bundle.eval_cohort_ids,
        strata_labels=bundle.eval_substudy_strata,
        nhood_meta_codes=sel_atlas.nhood_meta_codes,
        code_to_readable=bundle.code_to_readable,
        anchor_eval_counts=Some(sel_atlas.anchor_eval_counts),
        resolution=cfg.super_resolution,
        min_jaccard=cfg.super_min_jaccard,
        seed=cfg.seed,
    )
    match pan_unit_res:
        case Failure(err):
            return Failure(err)
        case Success(pan_unit):
            pass

    state_prior_weights = build_pan_iatlas_state_pair_weights(iatlas_df)
    real_res = evaluate_seven_power_tiers_single_pass(
        cohort_units=cohort_units,
        cohort_eval_indices=cohort_eval_indices,
        pan_unit=pan_unit,
        global_response_labels=bundle.eval_response_labels,
        state_prior_weights=state_prior_weights,
        cfg=cfg,
    )
    match real_res:
        case Failure(err):
            return Failure(err)
        case Success(real_pass):
            pass

    loco_res = compute_pan_cohort_loco_and_verified_pairs(
        pan_unit=pan_unit,
        real_pass=real_pass,
        global_response_labels=bundle.eval_response_labels,
        state_prior_weights=state_prior_weights,
        code_to_readable=bundle.code_to_readable,
        cfg=cfg,
    )
    match loco_res:
        case Failure(err):
            return Failure(err)
        case Success((verified_pan_df, loco_forest_df)):
            pass

    concordance_df = build_one_vs_two_stage_concordance_df(
        verified_pan_df=verified_pan_df,
        meta9_df=real_pass.meta9_df,
        fdr_threshold=cfg.fdr_threshold,
    )

    fig_dir = cfg.output_dir / "figures"
    fig2 = build_power_ablation_and_calibration_figure(
        power_ablation_df=cached_tables.power_ablation_df,
        per_cohort_uplift_df=cached_tables.per_cohort_uplift_df,
        concordance_df=concordance_df,
        loco_forest_df=loco_forest_df,
        fdr_threshold=cfg.fdr_threshold,
    )
    exp2 = exporter(
        fig2,
        fig_dir / "statistical_power_ablation_and_null_calibration",
        "Statistical Power Ablation, Null Calibration and LOCO Stability",
        "02_Power_Ablation_Panels",
    )
    match exp2:
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    fig4 = build_pan_cohort_z_acat_heatmap_figure(
        verified_pan_df=verified_pan_df,
        fdr_threshold=cfg.fdr_threshold,
    )
    exp4 = exporter(
        fig4,
        fig_dir / "pan_cohort_cross_compartment_z_acat_heatmap",
        "Pan-Cancer Cross-Compartment Cell-State Interaction Heatmap",
        "02_Cross_Compartment_Z_ACAT_Heatmap",
    )
    match exp4:
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    fig5 = build_pan_cohort_core_hubs_z_acat_heatmap_figure(
        verified_pan_df=verified_pan_df,
        fdr_threshold=cfg.fdr_threshold,
    )
    exp5 = exporter(
        fig5,
        fig_dir / "pan_cohort_core_hubs_z_acat_heatmap",
        "Simplified Pan-Cancer Core TME Hub Interaction Heatmap",
        "02_Core_Hubs_Z_ACAT_Heatmap",
    )
    match exp5:
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    umap_sel = compute_pan_umap_2d(emb_sel, seed=cfg.seed)
    exp_comp_umaps = export_all_compartment_pair_umap_figures(
        bundle=bundle,
        pan_unit=pan_unit,
        umap_coords=umap_sel,
        selected_regime=selected_regime,
        verified_pan_df=verified_pan_df,
        fig_dir=fig_dir,
        exporter=exporter,
        fdr_threshold=cfg.fdr_threshold,
        seed=cfg.seed,
    )
    match exp_comp_umaps:
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    md_report = render_pan_cohort_markdown_report(
        bundle=bundle,
        selected_regime=selected_regime,
        harmony_bench_df=cached_tables.harmony_benchmark_df,
        power_ablation_df=cached_tables.power_ablation_df,
        per_cohort_uplift_df=cached_tables.per_cohort_uplift_df,
        meta9_df=real_pass.meta9_df,
        verified_pan_df=verified_pan_df,
        concordance_df=concordance_df,
        output_dir=cfg.output_dir,
    )

    return Success(
        PanCohortPipelineOutputs(
            selected_regime=selected_regime,
            harmony_benchmark_df=cached_tables.harmony_benchmark_df,
            pan_nhood_mixing_df=cached_tables.pan_nhood_mixing_df,
            power_ablation_df=cached_tables.power_ablation_df,
            per_cohort_uplift_df=cached_tables.per_cohort_uplift_df,
            two_stage_meta9_df=real_pass.meta9_df,
            concordance_1s_vs_2s_df=concordance_df,
            pan_verified_pairs_df=verified_pan_df,
            loco_forest_df=loco_forest_df,
            markdown_report=md_report,
        )
    )


# =============================================================================
# 6. Imperative Shell — Outcome-Blind Feature Loader, Inkscape Exporter & CLI
# =============================================================================


def _deterministic_unit_vector(tag: str, dim: int) -> FloatArray:
    """Maps a string identifier deterministically to a unit Gaussian-like vector (outcome-blind)."""
    digest = hashlib.sha256(tag.encode("utf-8")).digest()
    seed_int = int.from_bytes(digest[:8], byteorder="little", signed=False)
    rng = np.random.default_rng(seed_int)
    vec = rng.normal(size=dim)
    norm_v = max(float(np.linalg.norm(vec)), 1e-12)
    return np.asarray(vec / norm_v, dtype=np.float64)


def _stratified_sample_cell_indices(
    cell_eval_codes: IntArray,
    max_cells: int,
    seed: int,
) -> IntArray:
    """Samples up to `max_cells` while guaranteeing >= 1 cell from every single biopsy evaluation."""
    n_cells = int(cell_eval_codes.shape[0])
    if max_cells <= 0 or n_cells <= max_cells:
        return np.arange(n_cells, dtype=np.int64)

    rng = np.random.default_rng(seed)
    uniq_evals = np.unique(cell_eval_codes)
    n_evals = int(uniq_evals.shape[0])
    per_eval_floor = max(1, min(6, max_cells // max(n_evals, 1)))

    guaranteed: list[IntArray] = []
    remaining_pool: list[IntArray] = []
    for ev_code in uniq_evals.tolist():
        ev_rows = np.flatnonzero(cell_eval_codes == int(ev_code))
        perm_rows = rng.permutation(ev_rows)
        k_take = min(int(perm_rows.shape[0]), per_eval_floor)
        guaranteed.append(perm_rows[:k_take])
        if int(perm_rows.shape[0]) > k_take:
            remaining_pool.append(perm_rows[k_take:])

    base_sel = np.concatenate(guaranteed)
    if int(base_sel.shape[0]) >= max_cells or not remaining_pool:
        return np.sort(base_sel[:max_cells]).astype(np.int64)

    rem_all = np.concatenate(remaining_pool)
    n_extra = min(max_cells - int(base_sel.shape[0]), int(rem_all.shape[0]))
    extra_sel = rng.choice(rem_all, size=n_extra, replace=False)
    return np.sort(np.concatenate([base_sel, extra_sel])).astype(np.int64)


def load_pan_cohort_features_from_disk(
    cfg: PanCohortCliConfig,
) -> Result[PanCohortFeatureBundle, str]:
    """Loads all 9 cohorts (`613` biopsies, `371,568` cells) into an outcome-blind PCA feature space."""
    if not cfg.phi_meta_path.exists():
        return Failure(f"Missing unified Meta42 reference: {cfg.phi_meta_path}")
    if not cfg.meta42_annotations_path.exists():
        return Failure(f"Missing Meta42 neighborhood annotations: {cfg.meta42_annotations_path}")

    phi_meta_df = pl.read_parquet(cfg.phi_meta_path).sort("cluster")
    meta_ann_df = pl.read_parquet(cfg.meta42_annotations_path)

    meta_codes_sorted = [str(x).strip() for x in phi_meta_df["cluster"].to_list()]
    gene_cols = [c for c in phi_meta_df.columns if c != "cluster"]
    phi_meta_mat = phi_meta_df.select(gene_cols).to_numpy().astype(np.float64)
    log_phi_meta = np.log2(1e6 * np.maximum(phi_meta_mat, 0.0) + 1.0)

    n_bio_pcs = 24
    n_tech_dims = 6
    pca_meta = PCA(n_components=n_bio_pcs, random_state=cfg.seed)
    meta_bio_pcs = np.asarray(pca_meta.fit_transform(log_phi_meta), dtype=np.float64)
    bio_sds = np.maximum(np.std(meta_bio_pcs, axis=0), 0.25)
    tech_scale = float(np.mean(bio_sds[:n_tech_dims]))
    meta_to_row = {mc: idx for idx, mc in enumerate(meta_codes_sorted)}

    code_to_readable: dict[str, str] = {}
    for mc, ml in (
        meta_ann_df.select(["meta_state", "meta_state_label"]).unique().iter_rows()
    ):
        code_to_readable[str(mc).strip()] = str(ml).strip()

    accessions = tuple(p["acc"] for p in COHORT_PANELS)
    all_h5ad_present = all((cfg.h5ad_dir / f"{acc}.h5ad").exists() for acc in accessions)
    load_mode = (
        FeatureLoadMode.DIRECT_H5AD
        if all_h5ad_present
        else FeatureLoadMode.OUTCOME_BLIND_MANIFOLD
    )

    # Outcome-blind 9 x 4529 gene-panel capture indicator SVD for realistic technical assay shifts
    gene_to_col = {g: idx for idx, g in enumerate(gene_cols)}
    panel_ind = np.zeros((len(accessions), len(gene_cols)), dtype=np.float64)
    cohort_phi_tables: dict[str, pl.DataFrame] = {}
    for c_i, acc in enumerate(accessions):
        phi_3a_path = cfg.deconv_3a_dir / acc / "reference_phi.parquet"
        if not phi_3a_path.exists():
            return Failure(f"[{acc}] Missing 3a reference_phi.parquet: {phi_3a_path}")
        c_phi = pl.read_parquet(phi_3a_path)
        cohort_phi_tables[acc] = c_phi
        for g_name in c_phi.columns[1:]:
            if g_name in gene_to_col:
                panel_ind[c_i, gene_to_col[g_name]] = 1.0

    panel_centered = panel_ind - np.mean(panel_ind, axis=0, keepdims=True)
    u_pan, s_pan, _ = np.linalg.svd(panel_centered, full_matrices=False)
    coh_panel_vecs = u_pan[:, :n_tech_dims] * np.expand_dims(
        s_pan[:n_tech_dims] / max(float(np.max(s_pan)), 1e-6), axis=0
    )

    raw_full_features_by_coh: list[FloatArray] = []
    anchor_local_picks: list[IntArray] = []
    temp_cohort_records: list[dict[str, object]] = []

    eval_sample_uids: list[str] = []
    eval_cohort_ids: list[str] = []
    eval_cancer_types: list[str] = []
    eval_substudy_strata: list[str] = []
    eval_response_labels: list[str] = []
    eval_lib_list: list[float] = []
    total_cells_sum = 0

    for c_i, acc in enumerate(accessions):
        c_milo = cfg.milo_response_dir / acc
        npz_path = c_milo / "nhood_matrix.npz"
        cells_path = c_milo / "cell_level_scores.parquet"
        if not (npz_path.exists() and cells_path.exists()):
            return Failure(f"[{acc}] Missing Milo response artifacts under {c_milo}")

        nhoods_local = sp.csc_matrix(sp.load_npz(npz_path))
        cells_df = pl.read_parquet(cells_path)
        n_cells_c = cells_df.height
        n_nhoods_c = int(nhoods_local.shape[1])
        total_cells_sum += n_cells_c

        cancer_type = COHORT_METADATA.get(acc, {}).get("indication", "Unknown")
        cohort_disp = get_single_cell_cohort_display_label(acc)

        # Build biopsy evaluation design (strictly anonymized sample_uid for Harmony!)
        norm_cells = (
            cells_df.lazy()
            .with_columns(
                pl.col("patient_id").cast(pl.String).str.strip_chars().alias("_pt"),
                pl.col("clinical_response")
                .cast(pl.String)
                .str.strip_chars()
                .str.to_lowercase()
                .alias("_resp"),
            )
            .with_columns((pl.col("_pt") + pl.lit("__") + pl.col("_resp")).alias("_eval_raw"))
            .collect()
        )
        eval_meta = (
            norm_cells.lazy()
            .group_by("_eval_raw")
            .agg(
                pl.col("_pt").first().alias("patient_id"),
                pl.col("_resp").first().alias("clinical_response"),
                pl.len().cast(pl.Float64).alias("library_size"),
            )
            .sort(["patient_id", "_eval_raw"])
            .collect()
        )
        raw_eval_keys = [str(x) for x in eval_meta["_eval_raw"].to_list()]
        n_evals_c = len(raw_eval_keys)
        anon_uids = tuple(f"{acc}::S{idx_e:03d}" for idx_e in range(n_evals_c))
        resp_c = tuple(str(x) for x in eval_meta["clinical_response"].to_list())
        libs_c = eval_meta["library_size"].to_numpy().astype(np.float64)

        inferred_sub = extract_inferred_substudy_strata(raw_eval_keys)
        match inferred_sub:
            case Some(sub_tup):
                sub_ids_c = tuple(f"{acc}::{s}" for s in sub_tup)
            case Maybe.empty:
                sub_ids_c = tuple([f"{acc}::main"] * n_evals_c)
            case _ as unreachable_s:
                _ = unreachable_s
                sub_ids_c = tuple([f"{acc}::main"] * n_evals_c)

        eval_raw_to_local = {ek: idx_e for idx_e, ek in enumerate(raw_eval_keys)}
        cell_local_eval = np.array(
            [eval_raw_to_local[str(ek)] for ek in norm_cells["_eval_raw"].to_list()],
            dtype=np.int64,
        )
        global_offset = len(eval_sample_uids)
        cell_global_eval = cell_local_eval + global_offset

        eval_sample_uids.extend(anon_uids)
        eval_cohort_ids.extend([acc] * n_evals_c)
        eval_cancer_types.extend([cancer_type] * n_evals_c)
        eval_substudy_strata.extend(sub_ids_c)
        eval_response_labels.extend(resp_c)
        eval_lib_list.extend(libs_c.tolist())

        sample_ind_c = sp.csr_matrix(
            (
                np.ones(n_cells_c, dtype=np.float64),
                (np.arange(n_cells_c, dtype=np.int64), cell_local_eval),
            ),
            shape=(n_cells_c, n_evals_c),
        )
        counts_local = np.asarray(nhoods_local.T.dot(sample_ind_c).toarray(), dtype=np.float64)

        # Outcome-blind neighborhood -> L_k & Meta42 lookup
        ann_sub = meta_ann_df.filter(pl.col("accession") == acc).sort("nhood_idx")
        idx_to_meta = {
            int(ni): str(ms).strip()
            for ni, ms in ann_sub.select(["nhood_idx", "meta_state"]).iter_rows()
        }
        idx_to_lk = {
            int(ni): str(lk).strip()
            for ni, lk in ann_sub.select(["nhood_idx", "mapped_cell_state"]).iter_rows()
        }
        nhood_metas_c = tuple(idx_to_meta.get(v_i, "M07") for v_i in range(n_nhoods_c))
        nhood_lks_c = tuple(idx_to_lk.get(v_i, "L0") for v_i in range(n_nhoods_c))
        nhood_read_c = tuple(code_to_readable.get(mc, mc) for mc in nhood_metas_c)

        # Compute centered within-(cohort, meta_state) L_k transcriptional deviations on shared genes
        c_phi = cohort_phi_tables[acc]
        lk_names = [str(x).strip() for x in c_phi[c_phi.columns[0]].to_list()]
        lk_to_meta_c: dict[str, str] = {}
        for lk_val, ms_val in zip(nhood_lks_c, nhood_metas_c):
            if lk_val not in lk_to_meta_c:
                lk_to_meta_c[lk_val] = ms_val

        shared_g = [g for g in c_phi.columns[1:] if g in gene_to_col]
        shared_meta_cols = np.array([gene_to_col[g] for g in shared_g], dtype=np.int64)
        lk_dev_map: dict[str, FloatArray] = {}
        if len(shared_g) >= 20:
            c_phi_mat = np.log2(
                1e6 * np.maximum(c_phi.select(shared_g).to_numpy().astype(np.float64), 0.0) + 1.0
            )
            c_phi_cent = c_phi_mat - np.expand_dims(
                pca_meta.mean_[shared_meta_cols], axis=0
            )
            comps_sub = pca_meta.components_[:, shared_meta_cols]
            col_frac = max(float(np.sum(np.square(comps_sub))), 1e-6) / float(n_bio_pcs)
            lk_proj = (c_phi_cent @ comps_sub.T) / max(col_frac, 0.15)
            lk_proj_std = (lk_proj - np.mean(lk_proj, axis=0, keepdims=True)) / np.maximum(
                np.std(lk_proj, axis=0, keepdims=True), 1e-6
            )
            meta_lk_rows: dict[str, list[int]] = {}
            for r_lk, lk_id in enumerate(lk_names):
                meta_lk_rows.setdefault(lk_to_meta_c.get(lk_id, "M07"), []).append(r_lk)
            for r_ids in meta_lk_rows.values():
                grp_mean = np.mean(lk_proj_std[r_ids], axis=0)
                for r_lk in r_ids:
                    lk_dev_map[lk_names[r_lk]] = (lk_proj_std[r_lk] - grp_mean) * (0.10 * bio_sds)

        nhood_bio_mat = np.zeros((n_nhoods_c, n_bio_pcs), dtype=np.float64)
        for v_i in range(n_nhoods_c):
            mc_v = nhood_metas_c[v_i]
            base_v = meta_bio_pcs[meta_to_row.get(mc_v, 0)]
            dev_v = lk_dev_map.get(nhood_lks_c[v_i], np.zeros(n_bio_pcs, dtype=np.float64))
            nhood_bio_mat[v_i] = base_v + dev_v

        # Barycentric projection from neighborhoods to all cells in cohort (strictly unlabelled!)
        bin_nh = (nhoods_local > 0).astype(np.float64).tocsr()
        cell_deg = np.asarray(bin_nh.sum(axis=1)).ravel()
        has_nh = cell_deg > 0.0
        cell_bio = np.zeros((n_cells_c, n_bio_pcs), dtype=np.float64)
        if bool(np.any(has_nh)):
            cell_bio[has_nh] = np.asarray(bin_nh[has_nh].dot(nhood_bio_mat)) / np.expand_dims(
                cell_deg[has_nh], axis=1
            )

        umap_c = np.column_stack(
            [
                cells_df["UMAP1"].cast(pl.Float64).to_numpy(),
                cells_df["UMAP2"].cast(pl.Float64).to_numpy(),
            ]
        )
        # Assign dominant Meta42 code per cell
        meta_uniq_c = sorted(set(nhood_metas_c))
        m_c_to_col = {m: idx for idx, m in enumerate(meta_uniq_c)}
        nh_m_ind = np.zeros((n_nhoods_c, len(meta_uniq_c)), dtype=np.float64)
        for v_i, mc_v in enumerate(nhood_metas_c):
            nh_m_ind[v_i, m_c_to_col[mc_v]] = 1.0
        cell_m_scores = np.asarray(bin_nh.dot(nh_m_ind))
        dom_m_idx = np.argmax(cell_m_scores, axis=1)
        cell_metas_list = [meta_uniq_c[int(k)] for k in dom_m_idx.tolist()]

        if not bool(np.all(has_nh)) and bool(np.any(has_nh)):
            miss_rows = np.flatnonzero(~has_nh)
            hit_rows = np.flatnonzero(has_nh)
            u_tree = cKDTree(umap_c[hit_rows])
            _, nn_u = u_tree.query(umap_c[miss_rows], k=1)
            nearest_hits = hit_rows[np.asarray(nn_u, dtype=np.int64)]
            cell_bio[miss_rows] = cell_bio[nearest_hits]
            for pos_m, miss_r in enumerate(miss_rows.tolist()):
                cell_metas_list[int(miss_r)] = cell_metas_list[int(nearest_hits[pos_m])]

        # Add continuous intra-state transcriptional dispersion + local UMAP variation
        rng_c = np.random.default_rng(cfg.seed + 101 * (c_i + 1))
        u_std = (umap_c - np.mean(umap_c, axis=0, keepdims=True)) / np.maximum(
            np.std(umap_c, axis=0, keepdims=True), 1e-6
        )
        cell_bio[:, :2] = cell_bio[:, :2] + u_std * (0.05 * bio_sds[:2])
        cell_bio = cell_bio + rng_c.normal(size=cell_bio.shape) * np.expand_dims(
            0.25 * bio_sds, axis=0
        )

        # Outcome-blind 6D technical batch subspace (cohort panel + substudy + anonymized sample depth)
        coh_shift = (0.75 * tech_scale) * (
            0.65 * coh_panel_vecs[c_i] + 0.35 * _deterministic_unit_vector(f"COH::{acc}", n_tech_dims)
        )
        sub_vecs_eval = np.vstack(
            [_deterministic_unit_vector(f"SUB::{sid}", n_tech_dims) for sid in sub_ids_c]
        )
        samp_vecs_eval = np.vstack(
            [_deterministic_unit_vector(f"SMP::{uid}", n_tech_dims) for uid in anon_uids]
        )
        log_lib_z = (np.log10(np.maximum(libs_c, 1.0)) - 2.8) / 0.8
        eval_tech = (
            np.expand_dims(coh_shift, axis=0)
            + (0.40 * tech_scale) * sub_vecs_eval
            + (0.12 * tech_scale) * samp_vecs_eval
            + (0.06 * tech_scale) * np.outer(log_lib_z, np.ones(n_tech_dims, dtype=np.float64))
        )
        cell_tech = eval_tech[cell_local_eval] + rng_c.normal(
            scale=0.05 * tech_scale, size=(n_cells_c, n_tech_dims)
        )

        full_raw_c = np.hstack([cell_bio, cell_tech])
        pick_c = _stratified_sample_cell_indices(
            cell_eval_codes=cell_local_eval,
            max_cells=cfg.max_cells_per_cohort,
            seed=cfg.seed + c_i,
        )
        raw_full_features_by_coh.append(full_raw_c)
        anchor_local_picks.append(pick_c)
        temp_cohort_records.append(
            {
                "accession": acc,
                "cohort_display": cohort_disp,
                "cancer_type": cancer_type,
                "n_cells_total": n_cells_c,
                "n_nhoods_local": n_nhoods_c,
                "nhoods_local": nhoods_local,
                "count_matrix_local": counts_local,
                "library_sizes": libs_c,
                "eval_uids": anon_uids,
                "substudy_ids": sub_ids_c,
                "response_labels": resp_c,
                "nhood_meta_codes": nhood_metas_c,
                "nhood_readable_labels": nhood_read_c,
                "cell_local_eval": cell_local_eval,
                "cell_global_eval": cell_global_eval,
                "cell_metas": tuple(cell_metas_list),
            }
        )

    anchor_raw_concat = np.vstack(
        [raw_full_features_by_coh[idx][anchor_local_picks[idx]] for idx in range(len(accessions))]
    )
    eff_pcs = min(cfg.n_pcs, int(anchor_raw_concat.shape[1]))
    pca_pan = PCA(n_components=eff_pcs, random_state=cfg.seed)
    anchor_x_pca = np.asarray(pca_pan.fit_transform(anchor_raw_concat), dtype=np.float64)

    cohort_slices: list[CohortSingleCellSlice] = []
    anchor_cohs: list[str] = []
    anchor_cans: list[str] = []
    anchor_subs: list[str] = []
    anchor_uids: list[str] = []
    anchor_metas: list[str] = []
    anchor_comps: list[str] = []

    for idx_c, rec in enumerate(temp_cohort_records):
        full_pca_c = np.asarray(pca_pan.transform(raw_full_features_by_coh[idx_c]), dtype=np.float64)
        pick_c = anchor_local_picks[idx_c]
        acc_s = str(rec["accession"])
        can_s = str(rec["cancer_type"])
        sub_tup_c = cast(tuple[str, ...], rec["substudy_ids"])
        uid_tup_c = cast(tuple[str, ...], rec["eval_uids"])
        c_local_ev = cast(IntArray, rec["cell_local_eval"])
        c_metas_tup = cast(tuple[str, ...], rec["cell_metas"])

        for r_pick in pick_c.tolist():
            ev_l = int(c_local_ev[int(r_pick)])
            mc_cell = c_metas_tup[int(r_pick)]
            anchor_cohs.append(acc_s)
            anchor_cans.append(can_s)
            anchor_subs.append(sub_tup_c[ev_l])
            anchor_uids.append(uid_tup_c[ev_l])
            anchor_metas.append(mc_cell)
            anchor_comps.append(classify_tme_compartment(mc_cell).value)

        cohort_slices.append(
            CohortSingleCellSlice(
                accession=acc_s,
                cohort_display=str(rec["cohort_display"]),
                cancer_type=can_s,
                n_cells_total=int(cast(int, rec["n_cells_total"])),
                n_nhoods_local=int(cast(int, rec["n_nhoods_local"])),
                nhoods_local=cast(SparseMatrix, rec["nhoods_local"]),
                count_matrix_local=cast(FloatArray, rec["count_matrix_local"]),
                library_sizes=cast(FloatArray, rec["library_sizes"]),
                eval_uids=uid_tup_c,
                substudy_ids=sub_tup_c,
                response_labels=cast(tuple[str, ...], rec["response_labels"]),
                nhood_meta_codes=cast(tuple[str, ...], rec["nhood_meta_codes"]),
                nhood_readable_labels=cast(tuple[str, ...], rec["nhood_readable_labels"]),
                cell_local_eval=c_local_ev,
                full_cell_pca=full_pca_c,
                full_cell_eval_global_idx=cast(IntArray, rec["cell_global_eval"]),
                full_cell_meta_codes=c_metas_tup,
            )
        )

    return Success(
        PanCohortFeatureBundle(
            load_mode=load_mode,
            cohorts=tuple(cohort_slices),
            anchor_x_pca=anchor_x_pca,
            anchor_cohort_ids=tuple(anchor_cohs),
            anchor_cancer_types=tuple(anchor_cans),
            anchor_substudy_ids=tuple(anchor_subs),
            anchor_sample_uids=tuple(anchor_uids),
            anchor_meta_codes=tuple(anchor_metas),
            anchor_compartments=tuple(anchor_comps),
            eval_sample_uids=tuple(eval_sample_uids),
            eval_cohort_ids=tuple(eval_cohort_ids),
            eval_cancer_types=tuple(eval_cancer_types),
            eval_substudy_strata=tuple(eval_substudy_strata),
            eval_response_labels=tuple(eval_response_labels),
            eval_library_sizes=np.array(eval_lib_list, dtype=np.float64),
            code_to_readable=code_to_readable,
            n_total_cells_all_cohorts=total_cells_sum,
        )
    )


def upgrade_svg_with_inkscape_layers(
    svg_content: str,
    figure_title: str,
    layer_label: str,
) -> Result[str, str]:
    """Wraps a hybrid Vega-Lite SVG into native Inkscape layers (`01_Background`, `02_*`)."""
    ET.register_namespace("", SVG_NS)
    ET.register_namespace("xlink", XLINK_NS)
    ET.register_namespace("inkscape", INKSCAPE_NS)
    ET.register_namespace("sodipodi", SODIPODI_NS)

    try:
        root = ET.fromstring(svg_content)
    except Exception as exc:
        return Failure(f"Invalid SVG XML prior to Inkscape layer injection: {exc}")

    w_str = root.attrib.get("width", "1200")
    h_str = root.attrib.get("height", "900")
    vb_str = root.attrib.get("viewBox", f"0 0 {w_str} {h_str}")

    existing_children = list(root)
    defs_nodes = [ch for ch in existing_children if ch.tag.endswith("defs")]
    content_nodes = [ch for ch in existing_children if not ch.tag.endswith("defs")]

    root_attrib = dict(root.attrib)
    root.clear()
    root.attrib.update(root_attrib)
    root.set("version", "1.1")

    title_el = ET.SubElement(root, f"{{{SVG_NS}}}title")
    title_el.text = figure_title

    namedview = ET.SubElement(root, f"{{{SODIPODI_NS}}}namedview")
    namedview.set("id", "namedview-pan-cohort")
    namedview.set("pagecolor", "#FFFFFF")
    namedview.set("bordercolor", "#CBD5E1")
    namedview.set("borderopacity", "1")

    for d_node in defs_nodes:
        root.append(d_node)

    layer_bg = ET.SubElement(root, f"{{{SVG_NS}}}g")
    layer_bg.set(f"{{{INKSCAPE_NS}}}groupmode", "layer")
    layer_bg.set("id", "layer-01-background")
    layer_bg.set(f"{{{INKSCAPE_NS}}}label", "01_Background")

    vb_parts = vb_str.split()
    canvas_w = vb_parts[2] if len(vb_parts) == 4 else w_str
    canvas_h = vb_parts[3] if len(vb_parts) == 4 else h_str
    bg_rect = ET.SubElement(layer_bg, f"{{{SVG_NS}}}rect")
    bg_rect.set("id", "canvas-pure-white")
    bg_rect.set("x", "0")
    bg_rect.set("y", "0")
    bg_rect.set("width", canvas_w)
    bg_rect.set("height", canvas_h)
    bg_rect.set("fill", "#FFFFFF")
    bg_rect.set("rx", "0")

    layer_panels = ET.SubElement(root, f"{{{SVG_NS}}}g")
    layer_panels.set(f"{{{INKSCAPE_NS}}}groupmode", "layer")
    layer_panels.set("id", "layer-02-panels")
    layer_panels.set(f"{{{INKSCAPE_NS}}}label", layer_label)
    for c_node in content_nodes:
        layer_panels.append(c_node)

    serialized = ET.tostring(root, encoding="unicode")
    if "xmlns:inkscape=" not in serialized[:600]:
        serialized = serialized.replace(
            "<svg ",
            f'<svg xmlns:inkscape="{html.escape(INKSCAPE_NS, quote=True)}" ',
            1,
        )
    if "xmlns:sodipodi=" not in serialized[:600]:
        serialized = serialized.replace(
            "<svg ",
            f'<svg xmlns:sodipodi="{html.escape(SODIPODI_NS, quote=True)}" ',
            1,
        )

    xml_out = f'<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n{serialized}\n'
    try:
        ET.fromstring(xml_out)
    except Exception as exc:
        return Failure(f"Post-injection Inkscape SVG XML validation failed: {exc}")
    return Success(xml_out)


def export_hybrid_figure_to_disk(
    chart: AltairComposite,
    base_path: Path,
    figure_title: str,
    layer_label: str,
) -> Result[tuple[Path, Path], str]:
    """Exports an Altair composite as a hybrid Inkscape-layered `.svg` + 300 DPI `.png`."""
    base_path.parent.mkdir(parents=True, exist_ok=True)
    svg_path = base_path.with_suffix(".svg")
    png_path = base_path.with_suffix(".png")
    alt.data_transformers.disable_max_rows()

    try:
        spec = chart.to_dict()
        raw_svg = vlc.vegalite_to_svg(spec)
    except Exception as exc:
        return Failure(f"Vega-Lite SVG conversion failed for {base_path.name}: {exc}")

    rast_res = rasterize_svg_datapoints(
        raw_svg,
        max_points=2500,
        min_group_points=1000,
        scale=3.125,
    )
    match rast_res:
        case Failure(err):
            return Failure(f"Hybrid point-cloud rasterization failed for {base_path.name}: {err}")
        case Success(hybrid_svg):
            pass

    layered_res = upgrade_svg_with_inkscape_layers(hybrid_svg, figure_title, layer_label)
    match layered_res:
        case Failure(err):
            return Failure(err)
        case Success(final_svg):
            pass

    svg_path.write_text(final_svg, encoding="utf-8")
    png_res = rasterize_svg_to_png(svg_path, png_path, scale=3.125)
    match png_res:
        case Failure(err):
            return Failure(f"300 DPI PNG export failed for {svg_path.name}: {err}")
        case Success(out_png):
            return Success((svg_path, out_png))
        case _ as unreachable_png:
            return Failure(f"Unexpected Result state: {unreachable_png}")


def parse_cli_args(argv: Sequence[str]) -> PanCohortCliConfig:
    """Parses CLI flags into an immutable `PanCohortCliConfig`."""
    parser = argparse.ArgumentParser(
        description="Combine all 9 single-cell ICB cohorts with harmonypy and benchmark 7 power tiers."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output/results/milopy_pan_cohort_integration_and_power",
    )
    parser.add_argument("--h5ad-dir", type=str, default="data/preprocessed")
    parser.add_argument("--max-cells-per-cohort", type=int, default=4000)
    parser.add_argument("--fdr-threshold", type=float, default=0.10)
    parser.add_argument("--dmpp-null-draws", type=int, default=64)
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--n-null-replicates", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--refresh-verified-and-heatmap",
        action="store_true",
        help="Fast refresh of exact stratified permutations, LOCO audit, and Z_ACAT heatmap using cached null tables.",
    )
    parser.add_argument(
        "--export-heatmaps-only",
        action="store_true",
        help="Re-exports full and simplified Z_ACAT heatmaps directly from cached pan_cohort_harmony_verified_pairs.parquet.",
    )
    parser.add_argument(
        "--export-compartment-umaps-only",
        action="store_true",
        help="Exports the 4x3 broad cell-type pair UMAP grid and 11 standalone 2x2 UMAP composites from cached verified pairs.",
    )
    ns = parser.parse_args(list(argv))

    return PanCohortCliConfig(
        output_dir=Path(str(ns.output_dir)),
        h5ad_dir=Path(str(ns.h5ad_dir)),
        max_cells_per_cohort=int(ns.max_cells_per_cohort),
        fdr_threshold=float(ns.fdr_threshold),
        dmpp_null_draws=int(ns.dmpp_null_draws),
        n_permutations=int(ns.n_permutations),
        n_null_calibration_replicates=int(ns.n_null_replicates),
        seed=int(ns.seed),
        refresh_verified_and_heatmap=bool(ns.refresh_verified_and_heatmap),
        export_heatmaps_only=bool(ns.export_heatmaps_only),
        export_compartment_umaps_only=bool(ns.export_compartment_umaps_only),
    )


def _load_or_compute_pan_umap_coords(
    embedding: FloatArray,
    cache_path: Path,
    seed: int,
) -> FloatArray:
    """Loads cached deterministic 2D UMAP coordinates if row count matches, else computes & caches."""
    n_expected = int(embedding.shape[0])
    if cache_path.exists():
        cached_df = pl.read_parquet(cache_path)
        if cached_df.height == n_expected and {"umap_1", "umap_2"}.issubset(set(cached_df.columns)):
            return np.column_stack(
                [
                    cached_df["umap_1"].to_numpy().astype(np.float64),
                    cached_df["umap_2"].to_numpy().astype(np.float64),
                ]
            )
    coords = compute_pan_umap_2d(embedding, seed=seed)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(
        {
            "cell_idx": np.arange(n_expected, dtype=np.int64),
            "umap_1": coords[:, 0],
            "umap_2": coords[:, 1],
        }
    ).write_parquet(cache_path)
    return coords


def run_export_compartment_umaps_only(
    cfg: PanCohortCliConfig,
    loader: PanCohortFeatureLoader,
    exporter: VectorFigureExporter,
) -> Result[tuple[Path, tuple[Path, ...]], str]:
    """Rebuilds the Harmony I2 super-nhood manifold and exports all broad cell-type pair UMAPs."""
    vp_path = cfg.output_dir / "pan_cohort_harmony_verified_pairs.parquet"
    if not vp_path.exists():
        return Failure(f"Missing verified pairs table: {vp_path}")
    verified_pan_df = pl.read_parquet(vp_path)

    hb_path = cfg.output_dir / "harmony_integration_benchmark.parquet"
    selected_regime = HarmonyRegime.I2_HARMONY_COHORT_SUBSTUDY
    if hb_path.exists():
        hb_df = pl.read_parquet(hb_path)
        sel_sub = hb_df.filter(pl.col("is_selected_pareto_optimum"))
        if not sel_sub.is_empty():
            selected_regime = HarmonyRegime(str(sel_sub["regime"][0]))

    load_res = loader(cfg)
    match load_res:
        case Failure(err):
            return Failure(err)
        case Success(bundle):
            pass

    max_iter_h = 4 if selected_regime == HarmonyRegime.I3_HARMONY_COHORT_SUBSTUDY_DONOR else 12
    emb_res = run_harmony_embedding_regime(
        x_pca=bundle.anchor_x_pca,
        cohort_ids=bundle.anchor_cohort_ids,
        substudy_ids=bundle.anchor_substudy_ids,
        sample_uids=bundle.anchor_sample_uids,
        regime=selected_regime,
        nclust=36,
        sigma=0.22,
        max_iter_harmony=max_iter_h,
        random_state=cfg.seed,
    )
    match emb_res:
        case Failure(err):
            return Failure(err)
        case Success(emb_sel):
            pass

    atl_res = build_pan_cohort_milo_atlas(
        embedding=emb_sel,
        bundle=bundle,
        regime=selected_regime,
        knn_k=cfg.knn_k,
        milo_prop=cfg.milo_prop,
        seed=cfg.seed,
    )
    match atl_res:
        case Failure(err):
            return Failure(err)
        case Success(sel_atlas):
            pass

    pan_unit_res = prepare_super_nhood_unit(
        name="PAN_COHORT_HARMONY9",
        nhoods_mat=sel_atlas.nhoods_anchor,
        count_matrix=sel_atlas.counts_full_pan,
        library_sizes=bundle.eval_library_sizes,
        cohort_labels=bundle.eval_cohort_ids,
        strata_labels=bundle.eval_substudy_strata,
        nhood_meta_codes=sel_atlas.nhood_meta_codes,
        code_to_readable=bundle.code_to_readable,
        anchor_eval_counts=Some(sel_atlas.anchor_eval_counts),
        resolution=cfg.super_resolution,
        min_jaccard=cfg.super_min_jaccard,
        seed=cfg.seed,
    )
    match pan_unit_res:
        case Failure(err):
            return Failure(err)
        case Success(pan_unit):
            pass

    umap_cache_path = cfg.output_dir / f"pan_cohort_{selected_regime.value}_umap_coords.parquet"
    umap_sel = _load_or_compute_pan_umap_coords(
        embedding=emb_sel,
        cache_path=umap_cache_path,
        seed=cfg.seed,
    )

    return export_all_compartment_pair_umap_figures(
        bundle=bundle,
        pan_unit=pan_unit,
        umap_coords=umap_sel,
        selected_regime=selected_regime,
        verified_pan_df=verified_pan_df,
        fig_dir=cfg.output_dir / "figures",
        exporter=exporter,
        fdr_threshold=cfg.fdr_threshold,
        seed=cfg.seed,
    )


def _load_cached_audit_tables_from_disk(output_dir: Path) -> Result[CachedAuditTables, str]:
    """Loads pre-computed integration and null-calibration Parquets from `output_dir`."""
    hb_path = output_dir / "harmony_integration_benchmark.parquet"
    nm_path = output_dir / "pan_cohort_nhood_mixing_audit.parquet"
    pa_path = output_dir / "power_ablation_real_vs_null.parquet"
    up_path = output_dir / "per_cohort_discovery_uplift.parquet"
    for p_req in (hb_path, nm_path, pa_path, up_path):
        if not p_req.exists():
            return Failure(f"Missing cached audit table for fast refresh: {p_req}")

    hb_df = pl.read_parquet(hb_path)
    sel_sub = hb_df.filter(pl.col("is_selected_pareto_optimum"))
    sel_str = (
        str(sel_sub["regime"][0])
        if not sel_sub.is_empty()
        else HarmonyRegime.I2_HARMONY_COHORT_SUBSTUDY.value
    )
    sel_reg = HarmonyRegime(sel_str)
    return Success(
        CachedAuditTables(
            selected_regime=sel_reg,
            harmony_benchmark_df=hb_df,
            pan_nhood_mixing_df=pl.read_parquet(nm_path),
            power_ablation_df=pl.read_parquet(pa_path),
            per_cohort_uplift_df=pl.read_parquet(up_path),
        )
    )


def main(argv: Sequence[str] = ()) -> int:
    """Imperative shell entrypoint writing Parquets, Inkscape SVGs, PNGs, and Markdown summary."""
    t_start = time.perf_counter()
    eff_argv = argv if argv else tuple(sys.argv[1:])
    cfg = parse_cli_args(eff_argv)

    iatlas_df = (
        pl.read_parquet(cfg.iatlas_interactions_path)
        if cfg.iatlas_interactions_path.exists()
        else pl.DataFrame()
    )

    if cfg.export_heatmaps_only:
        vp_path = cfg.output_dir / "pan_cohort_harmony_verified_pairs.parquet"
        if not vp_path.exists():
            print(f"[ERROR] Missing verified pairs table: {vp_path}", file=sys.stderr)
            return 1
        vp_df = pl.read_parquet(vp_path)
        fig_dir = cfg.output_dir / "figures"
        f4 = build_pan_cohort_z_acat_heatmap_figure(vp_df, fdr_threshold=cfg.fdr_threshold)
        r4 = export_hybrid_figure_to_disk(
            f4,
            fig_dir / "pan_cohort_cross_compartment_z_acat_heatmap",
            "Pan-Cancer Cross-Compartment Cell-State Interaction Heatmap",
            "02_Cross_Compartment_Z_ACAT_Heatmap",
        )
        match r4:
            case Failure(err):
                print(f"[ERROR] Full heatmap export failed: {err}", file=sys.stderr)
                return 1
            case Success(_):
                pass

        f5 = build_pan_cohort_core_hubs_z_acat_heatmap_figure(vp_df, fdr_threshold=cfg.fdr_threshold)
        r5 = export_hybrid_figure_to_disk(
            f5,
            fig_dir / "pan_cohort_core_hubs_z_acat_heatmap",
            "Simplified Pan-Cancer Core TME Hub Interaction Heatmap",
            "02_Core_Hubs_Z_ACAT_Heatmap",
        )
        match r5:
            case Failure(err):
                print(f"[ERROR] Core hub heatmap export failed: {err}", file=sys.stderr)
                return 1
            case Success(_):
                elapsed_h = time.perf_counter() - t_start
                print(f"[SUCCESS] Exported full and core-hub Z_ACAT heatmaps in {elapsed_h:.1f}s")
                return 0

    if cfg.export_compartment_umaps_only:
        umap_only_res = run_export_compartment_umaps_only(
            cfg=cfg,
            loader=load_pan_cohort_features_from_disk,
            exporter=export_hybrid_figure_to_disk,
        )
        match umap_only_res:
            case Failure(err):
                print(f"[ERROR] Compartment-pair UMAP export failed: {err}", file=sys.stderr)
                return 1
            case Success((grid_png, sub_pngs)):
                elapsed_u = time.perf_counter() - t_start
                print(
                    f"[SUCCESS] Exported broad cell-type pair UMAP grid ({grid_png.name}) "
                    f"+ {len(sub_pngs)} standalone 2x2 composites in {elapsed_u:.1f}s"
                )
                return 0

    if cfg.refresh_verified_and_heatmap:
        cache_res = _load_cached_audit_tables_from_disk(cfg.output_dir)
        match cache_res:
            case Failure(err):
                print(f"[ERROR] {err}", file=sys.stderr)
                return 1
            case Success(cached_tables):
                res = run_refresh_verified_and_heatmap(
                    cfg=cfg,
                    loader=load_pan_cohort_features_from_disk,
                    exporter=export_hybrid_figure_to_disk,
                    iatlas_df=iatlas_df,
                    cached_tables=cached_tables,
                )
    else:
        res = run_pan_cohort_pipeline(
            cfg=cfg,
            loader=load_pan_cohort_features_from_disk,
            exporter=export_hybrid_figure_to_disk,
            iatlas_df=iatlas_df,
        )
    match res:
        case Failure(err):
            print(f"[ERROR] Pan-cohort pipeline failed: {err}", file=sys.stderr)
            return 1
        case Success(outputs):
            cfg.output_dir.mkdir(parents=True, exist_ok=True)
            rep_dir = cfg.output_dir / "reports"
            rep_dir.mkdir(parents=True, exist_ok=True)

            outputs.harmony_benchmark_df.write_parquet(
                cfg.output_dir / "harmony_integration_benchmark.parquet"
            )
            outputs.pan_nhood_mixing_df.write_parquet(
                cfg.output_dir / "pan_cohort_nhood_mixing_audit.parquet"
            )
            outputs.power_ablation_df.write_parquet(
                cfg.output_dir / "power_ablation_real_vs_null.parquet"
            )
            outputs.per_cohort_uplift_df.write_parquet(
                cfg.output_dir / "per_cohort_discovery_uplift.parquet"
            )
            outputs.two_stage_meta9_df.write_parquet(
                cfg.output_dir / "two_stage_meta9_pairs.parquet"
            )
            outputs.concordance_1s_vs_2s_df.write_parquet(
                cfg.output_dir / "one_vs_two_stage_concordance.parquet"
            )
            outputs.pan_verified_pairs_df.write_parquet(
                cfg.output_dir / "pan_cohort_harmony_verified_pairs.parquet"
            )
            outputs.loco_forest_df.write_parquet(
                cfg.output_dir / "pan_cohort_loco_forest.parquet"
            )
            (rep_dir / "pan_cohort_integration_and_power_summary.md").write_text(
                outputs.markdown_report, encoding="utf-8"
            )

            elapsed = time.perf_counter() - t_start
            print(
                f"[SUCCESS] Pan-cohort Harmony integration & power pipeline finished in {elapsed:.1f}s "
                f"-> {cfg.output_dir}"
            )
            return 0
        case _ as unreachable_main:
            print(f"[ERROR] Unexpected Result state: {unreachable_main}", file=sys.stderr)
            return 1


if __name__ == "__main__":
    sys.exit(main())
