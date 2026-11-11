"""Pure-functional UMAP geometry, projection, and Nature Methods Altair charts for pairwise Milopy.

Solves the four core challenges of visualizing 2-body neighborhood interactions `(u, v)` on a 2D
single-cell UMAP embedding:
1. Manifold-faithful geometric medoids (`m_u = x_{i*}`, `i* in C_u`) preventing floating centroids
   in empty white space across UMAP tears.
2. Ecosystem-module bundled outward-bowing quadratic Bezier chords with polar endpoint glyphs
   (`+` filled Okabe Blue vs. `-` white-cored Okabe Vermilion ring) preventing hairball clutter.
3. De novo pairwise hubness projection isolating neighborhoods invisible to univariate Milo
   (`FDR_univ >= alpha`) that emerge strongly in pairwise testing (`FDR_IUT < alpha`).
4. Seed-conditioned 1-to-many partner scalar fields and single-pair kNN cell-cloud spotlights with
   empirical 2x2 donor responder rate Wilson 95% confidence interval insets.
"""

from __future__ import annotations

from collections.abc import Sequence
from enum import StrEnum
from typing import TypeAlias, assert_never, cast

import altair as alt
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import sparse as sp
from scipy.stats import norm

from .presence import estimate_nhood_presence
from .types import (
    ContinuousAxis,
    EnrichmentDirection,
    FloatArray,
    PairwiseAnalysisResult,
    PresenceConfig,
    PresenceMode,
    QuadrantState,
    SparseMatrix,
)

AltairComposite: TypeAlias = alt.Chart | alt.LayerChart | alt.HConcatChart | alt.VConcatChart

# Nature Methods minimalist wireframe & Okabe-Ito colorblind-safe tokens
_BG_WHITE: str = "#FFFFFF"
_HAIRLINE: str = "#CBD5E1"
_LIGHT_GRID: str = "#F1F5F9"
_CELL_INACTIVE: str = "#E2E8F0"
_TEXT_PRIMARY: str = "#0F172A"
_TEXT_SECONDARY: str = "#475569"

_OKABE_BLUE: str = "#0072B2"
_OKABE_VERMILION: str = "#D55E00"
_OKABE_ORANGE: str = "#E69F00"
_OKABE_BLUISH_GREEN: str = "#009E73"
_OKABE_REDDISH_PURPLE: str = "#CC79A7"
_OKABE_SKY_BLUE: str = "#56B4E9"


# =============================================================================
# 1. Exhaustive Sum Types & Frozen Domain Models
# =============================================================================


class NhoodCoordinateMethod(StrEnum):
    """Strategy for anchoring a kNN neighborhood onto 2D UMAP coordinates."""

    GEOMETRIC_MEDOID = "geometric_medoid"
    WEIGHTED_CENTROID = "weighted_centroid"


class PairwiseVisualRegime(StrEnum):
    """Five canonical pairwise visual interaction regimes on UMAP chords."""

    SYNERGY_PP = "++ Co-Presence Synergy"
    GATING_DISCORDANT = "+/- Mutual-Exclusion Gate"
    DESERT_MM = "-- Dual-Absence Desert"
    BALANCE_RATIO = "Log-Ratio Balance"
    CO_ABUNDANCE = "Co-Abundance Coupling"


ALL_VISUAL_REGIMES: tuple[PairwiseVisualRegime, ...] = (
    PairwiseVisualRegime.SYNERGY_PP,
    PairwiseVisualRegime.GATING_DISCORDANT,
    PairwiseVisualRegime.DESERT_MM,
    PairwiseVisualRegime.BALANCE_RATIO,
    PairwiseVisualRegime.CO_ABUNDANCE,
)


def regime_color_hex(regime: PairwiseVisualRegime) -> str:
    """Maps each pairwise interaction regime to its Okabe-Ito colorblind-safe hex code."""
    match regime:
        case PairwiseVisualRegime.SYNERGY_PP:
            return _OKABE_BLUE
        case PairwiseVisualRegime.GATING_DISCORDANT:
            return _OKABE_VERMILION
        case PairwiseVisualRegime.DESERT_MM:
            return _OKABE_ORANGE
        case PairwiseVisualRegime.BALANCE_RATIO:
            return _OKABE_BLUISH_GREEN
        case PairwiseVisualRegime.CO_ABUNDANCE:
            return _OKABE_REDDISH_PURPLE
        case _ as unreachable:
            assert_never(unreachable)


class EndpointPolarity(StrEnum):
    """Biological role of an individual neighborhood endpoint within a pairwise chord."""

    POSITIVE_PRESENT = "Present / Enriched (+)"
    NEGATIVE_GATED = "Excluded / Gated (-)"
    DUAL_ABSENT = "Dual-Absent (--)"
    SYMMETRIC_COABUND = "Co-Abundant"


class HubNoveltyStatus(StrEnum):
    """Comparison between 1D univariate Milo status and 2D pairwise Milo participation."""

    DE_NOVO_HUB = "De Novo Pairwise Hub (Univ. NS)"
    SHARED_HUB = "Shared Univ. + Pairwise Hub"
    UNIVARIATE_ONLY = "Univariate Only (No Sig. Pair)"
    INACTIVE = "Inactive Background"


class UmapPairVizConfig(BaseModel):
    """Immutable configuration controlling pairwise UMAP geometry and visual density."""

    model_config = ConfigDict(frozen=True)

    max_chords: int = Field(default=32, ge=1, le=200)
    max_anchor_chords: int = Field(default=14, ge=1, le=80)
    max_vector_hub_overlays: int = Field(default=160, ge=10, le=240)
    bezier_steps: int = Field(default=21, ge=5, le=101)
    base_curvature: float = Field(default=0.22, ge=0.0, le=0.80)
    min_spatial_chord_sep_frac: float = Field(default=0.035, ge=0.0, le=0.50)
    max_bg_cells: int = Field(default=12000, ge=200, le=100000)
    fdr_threshold: float = Field(default=0.10, gt=0.0, le=1.0)
    prefer_denovo_spotlight: bool = Field(default=True)
    allow_exploratory_fallback: bool = Field(default=False)
    coord_method: NhoodCoordinateMethod = Field(default=NhoodCoordinateMethod.GEOMETRIC_MEDOID)
    panel_width: int = Field(default=295, ge=120)
    panel_height: int = Field(default=255, ge=120)
    seed: int = Field(default=42)


class SpotlightPairMetadata(BaseModel):
    """Immutable summary of the featured #1 lead neighborhood pair in Panel d."""

    model_config = ConfigDict(frozen=True)

    has_pair: bool
    is_exploratory: bool
    nhood_u: int
    nhood_v: int
    label_u: str
    label_v: str
    regime_label: str
    direction_label: str
    novelty_label: str
    fdr_iut: float
    pvalue_iut: float
    effect_value: float
    univariate_fdr_u: float
    univariate_fdr_v: float
    n_cells_u: int
    n_cells_v: int


class PairwiseUmapBundle(BaseModel):
    """Complete pure-DataFrame bundle backing all four pairwise UMAP panels."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    bg_cells_df: pl.DataFrame
    nhood_nodes_df: pl.DataFrame
    hub_overlay_df: pl.DataFrame
    bundled_arcs_df: pl.DataFrame
    arc_endpoints_df: pl.DataFrame
    anchor_label: str
    anchor_bg_cells_df: pl.DataFrame
    anchor_overlay_df: pl.DataFrame
    anchor_arcs_df: pl.DataFrame
    spotlight_cells_df: pl.DataFrame
    spotlight_arc_df: pl.DataFrame
    spotlight_endpoints_df: pl.DataFrame
    spotlight_quadrant_df: pl.DataFrame
    spotlight_meta: SpotlightPairMetadata
    umap_bounds: tuple[float, float, float, float]
    n_sig_quadrant: int
    n_sig_continuous: int
    n_denovo_nhoods: int


# =============================================================================
# 2. Pure Mathematical Helpers (Medoids, Wilson CIs & Swap-Invariant Beziers)
# =============================================================================


def _p_to_z_mag(pval: float) -> float:
    """Maps a two-sided p-value in (0, 1] to a non-negative standard normal Z score."""
    p_safe = float(np.clip(pval, 1e-300, 1.0))
    return 0.0 if p_safe >= 1.0 else float(norm.isf(0.5 * p_safe))


def wilson_binomial_ci_pct(n_success: int, n_total: int) -> tuple[float, float]:
    """Exact Wilson score 95% confidence interval expressed on [0, 100] percentage scale."""
    if n_total <= 0:
        return (0.0, 0.0)
    z = 1.959963984540054
    z2 = z * z
    n_f = float(n_total)
    p_hat = float(n_success) / n_f
    denom = 1.0 + z2 / n_f
    center = (p_hat + z2 / (2.0 * n_f)) / denom
    half = (z / denom) * float(np.sqrt((p_hat * (1.0 - p_hat)) / n_f + z2 / (4.0 * n_f * n_f)))
    return (
        100.0 * float(np.clip(center - half, 0.0, 1.0)),
        100.0 * float(np.clip(center + half, 0.0, 1.0)),
    )


def compute_umap_bounding_domain(
    umap_coords: FloatArray,
    pad_fraction: float = 0.055,
) -> tuple[float, float, float, float]:
    """Computes shared padded `(xmin, xmax, ymin, ymax)` domain across all UMAP panels."""
    coords = np.asarray(umap_coords, dtype=np.float64)
    if coords.shape[0] == 0:
        return (-1.0, 1.0, -1.0, 1.0)
    x0, x1 = float(np.min(coords[:, 0])), float(np.max(coords[:, 0]))
    y0, y1 = float(np.min(coords[:, 1])), float(np.max(coords[:, 1]))
    dx = max((x1 - x0) * pad_fraction, 0.25)
    dy = max((y1 - y0) * pad_fraction, 0.25)
    return (x0 - dx, x1 + dx, y0 - dy, y1 + dy)


def compute_nhood_umap_medoids(
    nhoods_mat: SparseMatrix,
    umap_coords: FloatArray,
    state_labels: Maybe[Sequence[str]] = Nothing,
    method: NhoodCoordinateMethod = NhoodCoordinateMethod.GEOMETRIC_MEDOID,
) -> Result[pl.DataFrame, str]:
    """Computes manifold-faithful 2D UMAP coordinates for every neighborhood column in `nhoods_mat`.

    Under `GEOMETRIC_MEDOID`, each neighborhood `u` is anchored at the actual member cell
    `i* in C_u` minimizing squared Euclidean distance to the coordinate-wise median of `C_u`.
    Every medoid strictly coincides with an observed cell coordinate on the UMAP manifold.
    """
    coords = np.asarray(umap_coords, dtype=np.float64)
    if coords.ndim != 2 or int(coords.shape[1]) != 2:
        return Failure(f"umap_coords must have shape (n_cells, 2), got {coords.shape}")
    if not bool(np.all(np.isfinite(coords))):
        return Failure("umap_coords must contain strictly finite values")

    n_cells = int(nhoods_mat.shape[0])
    n_nhoods = int(nhoods_mat.shape[1])
    if int(coords.shape[0]) != n_cells:
        return Failure(
            f"umap_coords rows ({coords.shape[0]}) must match nhoods_mat rows ({n_cells})"
        )
    if n_cells == 0 or n_nhoods == 0:
        return Failure(f"Require positive cells and neighborhoods, got ({n_cells}, {n_nhoods})")

    match state_labels:
        case Some(lbls):
            if len(lbls) != n_nhoods:
                return Failure(
                    f"state_labels length ({len(lbls)}) must match n_nhoods ({n_nhoods})"
                )
            resolved_labels = tuple(str(x) for x in lbls)
        case Maybe.empty:
            resolved_labels = tuple(f"Nhood_{u}" for u in range(n_nhoods))
        case _ as unreachable_lbl:
            return Failure(f"Unexpected Maybe state: {unreachable_lbl}")

    csc = sp.csc_matrix(nhoods_mat)
    indptr = csc.indptr
    indices = csc.indices
    data = csc.data

    global_med_x = float(np.median(coords[:, 0]))
    global_med_y = float(np.median(coords[:, 1]))

    medoid_cells = np.full(n_nhoods, -1, dtype=np.int64)
    u1_arr = np.empty(n_nhoods, dtype=np.float64)
    u2_arr = np.empty(n_nhoods, dtype=np.float64)
    n_mem_arr = np.zeros(n_nhoods, dtype=np.int64)
    rad_arr = np.zeros(n_nhoods, dtype=np.float64)

    for u in range(n_nhoods):
        start, end = int(indptr[u]), int(indptr[u + 1])
        raw_rows = indices[start:end]
        raw_vals = data[start:end]
        pos_rows = raw_rows[raw_vals > 0]
        m_cnt = int(pos_rows.shape[0])
        n_mem_arr[u] = m_cnt
        if m_cnt == 0:
            u1_arr[u] = global_med_x
            u2_arr[u] = global_med_y
            continue

        pts = coords[pos_rows, :]
        match method:
            case NhoodCoordinateMethod.GEOMETRIC_MEDOID:
                med_xy = np.median(pts, axis=0)
                diff = pts - med_xy
                d2 = diff[:, 0] * diff[:, 0] + diff[:, 1] * diff[:, 1]
                best_loc = int(np.argmin(d2))
                chosen_cell = int(pos_rows[best_loc])
                medoid_cells[u] = chosen_cell
                u1_arr[u] = float(coords[chosen_cell, 0])
                u2_arr[u] = float(coords[chosen_cell, 1])
                d_chosen = pts - coords[chosen_cell, :]
                rad_arr[u] = float(
                    np.median(np.sqrt(d_chosen[:, 0] * d_chosen[:, 0] + d_chosen[:, 1] * d_chosen[:, 1]))
                )
            case NhoodCoordinateMethod.WEIGHTED_CENTROID:
                w = np.asarray(raw_vals[raw_vals > 0], dtype=np.float64)
                w_sum = max(float(np.sum(w)), 1e-12)
                cx = float(np.dot(w, pts[:, 0]) / w_sum)
                cy = float(np.dot(w, pts[:, 1]) / w_sum)
                diff = pts - np.array([cx, cy], dtype=np.float64)
                best_loc = int(np.argmin(diff[:, 0] * diff[:, 0] + diff[:, 1] * diff[:, 1]))
                medoid_cells[u] = int(pos_rows[best_loc])
                u1_arr[u] = cx
                u2_arr[u] = cy
                rad_arr[u] = float(
                    np.median(np.sqrt(diff[:, 0] * diff[:, 0] + diff[:, 1] * diff[:, 1]))
                )
            case _ as unreachable_m:
                assert_never(unreachable_m)

    return Success(
        pl.DataFrame(
            {
                "nhood_idx": np.arange(n_nhoods, dtype=np.int64),
                "medoid_cell_idx": medoid_cells,
                "umap_1": u1_arr,
                "umap_2": u2_arr,
                "n_member_cells": n_mem_arr,
                "radius_mad": rad_arr,
                "state_label": list(resolved_labels),
            }
        )
    )


def discretize_outward_bezier_arc(
    p_u: tuple[float, float],
    p_v: tuple[float, float],
    global_center: tuple[float, float],
    curvature: float,
    n_steps: int,
) -> tuple[FloatArray, FloatArray]:
    """Discretizes a swap-invariant outward-bowing quadratic Bezier curve from `p_u` to `p_v`.

    Canonicalizes the chord tangent lexicographically before choosing the outward normal away
    from `global_center`, guaranteeing `B_{vu}(1 - t) == B_{uv}(t)` to machine precision.
    """
    xu, yu = float(p_u[0]), float(p_u[1])
    xv, yv = float(p_v[0]), float(p_v[1])
    gx, gy = float(global_center[0]), float(global_center[1])

    steps = max(int(n_steps), 3)
    t_vec = np.linspace(0.0, 1.0, steps, dtype=np.float64)

    # Canonical lexicographic ordering ensures normal sign is 100% invariant to (u, v) <-> (v, u)
    if (xu, yu) <= (xv, yv):
        x0, y0, x1, y1 = xu, yu, xv, yv
    else:
        x0, y0, x1, y1 = xv, yv, xu, yu

    dx = x1 - x0
    dy = y1 - y0
    chord_len = float(np.hypot(dx, dy))
    if chord_len < 1e-12 or curvature <= 0.0:
        bx = (1.0 - t_vec) * xu + t_vec * xv
        by = (1.0 - t_vec) * yu + t_vec * yv
        return bx, by

    mx = 0.5 * (x0 + x1)
    my = 0.5 * (y0 + y1)
    nx = -dy / chord_len
    ny = dx / chord_len

    # Pick normal orientation pointing strictly further from global UMAP median
    d_pos = (mx + nx - gx) ** 2 + (my + ny - gy) ** 2
    d_neg = (mx - nx - gx) ** 2 + (my - ny - gy) ** 2
    s_out = 1.0 if d_pos >= d_neg else -1.0

    offset = curvature * chord_len
    cx = mx + s_out * offset * nx
    cy = my + s_out * offset * ny

    one_m_t = 1.0 - t_vec
    bx = (one_m_t * one_m_t) * xu + (2.0 * one_m_t * t_vec) * cx + (t_vec * t_vec) * xv
    by = (one_m_t * one_m_t) * yu + (2.0 * one_m_t * t_vec) * cy + (t_vec * t_vec) * yv
    return bx, by


def classify_quadrant_visual_regime_and_polarity(
    quadrant: QuadrantState,
    direction: EnrichmentDirection,
) -> tuple[PairwiseVisualRegime, EndpointPolarity, EndpointPolarity]:
    """Maps `(quadrant, direction)` to visual chord regime and `(polarity_u, polarity_v)`."""
    match quadrant:
        case QuadrantState.BOTH_PRESENT:
            return (
                (
                    PairwiseVisualRegime.SYNERGY_PP,
                    EndpointPolarity.POSITIVE_PRESENT,
                    EndpointPolarity.POSITIVE_PRESENT,
                )
                if direction == EnrichmentDirection.RESPONDER
                else (
                    PairwiseVisualRegime.GATING_DISCORDANT,
                    EndpointPolarity.NEGATIVE_GATED,
                    EndpointPolarity.NEGATIVE_GATED,
                )
            )
        case QuadrantState.BOTH_MISSING:
            return (
                PairwiseVisualRegime.DESERT_MM,
                EndpointPolarity.DUAL_ABSENT,
                EndpointPolarity.DUAL_ABSENT,
            )
        case QuadrantState.FIRST_ONLY:
            return (
                (
                    PairwiseVisualRegime.GATING_DISCORDANT,
                    EndpointPolarity.POSITIVE_PRESENT,
                    EndpointPolarity.NEGATIVE_GATED,
                )
                if direction == EnrichmentDirection.RESPONDER
                else (
                    PairwiseVisualRegime.GATING_DISCORDANT,
                    EndpointPolarity.NEGATIVE_GATED,
                    EndpointPolarity.POSITIVE_PRESENT,
                )
            )
        case QuadrantState.SECOND_ONLY:
            return (
                (
                    PairwiseVisualRegime.GATING_DISCORDANT,
                    EndpointPolarity.NEGATIVE_GATED,
                    EndpointPolarity.POSITIVE_PRESENT,
                )
                if direction == EnrichmentDirection.RESPONDER
                else (
                    PairwiseVisualRegime.GATING_DISCORDANT,
                    EndpointPolarity.POSITIVE_PRESENT,
                    EndpointPolarity.NEGATIVE_GATED,
                )
            )
        case _ as unreachable:
            assert_never(unreachable)


def classify_continuous_visual_regime_and_polarity(
    axis: ContinuousAxis,
    direction: EnrichmentDirection,
) -> tuple[PairwiseVisualRegime, EndpointPolarity, EndpointPolarity]:
    """Maps continuous `(axis, direction)` to visual chord regime and `(polarity_u, polarity_v)`."""
    match axis:
        case ContinuousAxis.LOG_RATIO_BALANCE:
            return (
                (
                    PairwiseVisualRegime.BALANCE_RATIO,
                    EndpointPolarity.POSITIVE_PRESENT,
                    EndpointPolarity.NEGATIVE_GATED,
                )
                if direction == EnrichmentDirection.RESPONDER
                else (
                    PairwiseVisualRegime.BALANCE_RATIO,
                    EndpointPolarity.NEGATIVE_GATED,
                    EndpointPolarity.POSITIVE_PRESENT,
                )
            )
        case ContinuousAxis.CO_ABUNDANCE_PRODUCT:
            return (
                (
                    PairwiseVisualRegime.CO_ABUNDANCE,
                    EndpointPolarity.SYMMETRIC_COABUND,
                    EndpointPolarity.SYMMETRIC_COABUND,
                )
                if direction == EnrichmentDirection.RESPONDER
                else (
                    PairwiseVisualRegime.CO_ABUNDANCE,
                    EndpointPolarity.NEGATIVE_GATED,
                    EndpointPolarity.NEGATIVE_GATED,
                )
            )
        case _ as unreachable:
            assert_never(unreachable)


# =============================================================================
# 3. Pure Hubness Projection, Chord Bundling, Anchor Field & Lead-Pair Spotlight
# =============================================================================


def summarize_nhood_pairwise_hubness(
    nhoods_mat: SparseMatrix,
    umap_coords: FloatArray,
    medoids_df: pl.DataFrame,
    pairwise_res: PairwiseAnalysisResult,
    univariate_fdr: Maybe[Sequence[float]] = Nothing,
    univariate_logfc: Maybe[Sequence[float]] = Nothing,
    cell_labels: Maybe[Sequence[str]] = Nothing,
    config: UmapPairVizConfig = UmapPairVizConfig(),
) -> Result[tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame], str]:
    """Projects pairwise significance degree & de novo emergence onto neighborhoods and cells."""
    n_cells = int(nhoods_mat.shape[0])
    n_nhoods = int(nhoods_mat.shape[1])
    if medoids_df.height != n_nhoods:
        return Failure(
            f"medoids_df height ({medoids_df.height}) must match n_nhoods ({n_nhoods})"
        )

    match univariate_fdr:
        case Some(fdrs):
            if len(fdrs) != n_nhoods:
                return Failure(f"univariate_fdr length ({len(fdrs)}) != n_nhoods ({n_nhoods})")
            univ_fdr_arr = np.asarray(fdrs, dtype=np.float64)
        case Maybe.empty:
            univ_fdr_arr = np.ones(n_nhoods, dtype=np.float64)
        case _ as unreachable_f:
            return Failure(f"Unexpected Maybe state: {unreachable_f}")

    match univariate_logfc:
        case Some(lfcs):
            if len(lfcs) != n_nhoods:
                return Failure(f"univariate_logfc length ({len(lfcs)}) != n_nhoods ({n_nhoods})")
            univ_lfc_arr = np.asarray(lfcs, dtype=np.float64)
        case Maybe.empty:
            univ_lfc_arr = np.zeros(n_nhoods, dtype=np.float64)
        case _ as unreachable_l:
            return Failure(f"Unexpected Maybe state: {unreachable_l}")

    match cell_labels:
        case Some(clbls):
            if len(clbls) != n_cells:
                return Failure(f"cell_labels length ({len(clbls)}) != n_cells ({n_cells})")
            resolved_cell_labels = list(str(c) for c in clbls)
        case Maybe.empty:
            resolved_cell_labels = ["Cell"] * n_cells
        case _ as unreachable_c:
            return Failure(f"Unexpected Maybe state: {unreachable_c}")

    partners: list[set[int]] = [set() for _ in range(n_nhoods)]
    quad_deg = np.zeros(n_nhoods, dtype=np.int64)
    cont_deg = np.zeros(n_nhoods, dtype=np.int64)
    denovo_deg = np.zeros(n_nhoods, dtype=np.int64)
    peak_neglog10_q = np.zeros(n_nhoods, dtype=np.float64)

    alpha = config.fdr_threshold
    univ_is_sig = univ_fdr_arr < alpha

    if not pairwise_res.pairs_df.is_empty():
        sig_q = pairwise_res.pairs_df.filter(pl.col("is_significant"))
        for r in sig_q.select(["nhood_u", "nhood_v", "fdr_iut"]).to_dicts():
            u = int(r["nhood_u"])
            v = int(r["nhood_v"])
            if 0 <= u < n_nhoods and 0 <= v < n_nhoods:
                partners[u].add(v)
                partners[v].add(u)
                quad_deg[u] += 1
                quad_deg[v] += 1
                nlq = float(-np.log10(max(float(r["fdr_iut"]), 1e-12)))
                peak_neglog10_q[u] = max(peak_neglog10_q[u], nlq)
                peak_neglog10_q[v] = max(peak_neglog10_q[v], nlq)
                if (not bool(univ_is_sig[u])) and (not bool(univ_is_sig[v])):
                    denovo_deg[u] += 1
                    denovo_deg[v] += 1

    if not pairwise_res.continuous_pairs_df.is_empty():
        sig_c = pairwise_res.continuous_pairs_df.filter(pl.col("is_significant"))
        for r in sig_c.select(["nhood_u", "nhood_v", "fdr_iut"]).to_dicts():
            u = int(r["nhood_u"])
            v = int(r["nhood_v"])
            if 0 <= u < n_nhoods and 0 <= v < n_nhoods:
                partners[u].add(v)
                partners[v].add(u)
                cont_deg[u] += 1
                cont_deg[v] += 1
                nlq = float(-np.log10(max(float(r["fdr_iut"]), 1e-12)))
                peak_neglog10_q[u] = max(peak_neglog10_q[u], nlq)
                peak_neglog10_q[v] = max(peak_neglog10_q[v], nlq)
                if (not bool(univ_is_sig[u])) and (not bool(univ_is_sig[v])):
                    denovo_deg[u] += 1
                    denovo_deg[v] += 1

    n_partners = np.array([len(s) for s in partners], dtype=np.int64)
    log10_deg = np.log10(1.0 + n_partners.astype(np.float64))

    hub_statuses: list[str] = []
    for u in range(n_nhoods):
        has_pw = int(n_partners[u]) > 0
        has_uv = bool(univ_is_sig[u])
        if has_pw and (not has_uv):
            hub_statuses.append(HubNoveltyStatus.DE_NOVO_HUB.value)
        elif has_pw and has_uv:
            hub_statuses.append(HubNoveltyStatus.SHARED_HUB.value)
        elif (not has_pw) and has_uv:
            hub_statuses.append(HubNoveltyStatus.UNIVARIATE_ONLY.value)
        else:
            hub_statuses.append(HubNoveltyStatus.INACTIVE.value)

    nhood_nodes_df = medoids_df.with_columns(
        pl.Series("univariate_fdr", univ_fdr_arr),
        pl.Series("univariate_logfc", univ_lfc_arr),
        pl.Series("is_univariate_sig", univ_is_sig),
        pl.Series("n_sig_partners", n_partners),
        pl.Series("log10_partner_degree", log10_deg),
        pl.Series("n_sig_quadrant_hits", quad_deg),
        pl.Series("n_sig_continuous_hits", cont_deg),
        pl.Series("n_denovo_hits", denovo_deg),
        pl.Series("peak_neg_log10_fdr_iut", peak_neglog10_q),
        pl.Series("hub_novelty_status", hub_statuses),
    )

    # Top active/de-novo vector overlay nodes (strictly capped < 250 so SVG marks remain vector)
    active_hubs = nhood_nodes_df.filter(pl.col("n_sig_partners") > 0)
    denovo_part = active_hubs.filter(
        pl.col("hub_novelty_status") == HubNoveltyStatus.DE_NOVO_HUB.value
    ).sort(["peak_neg_log10_fdr_iut", "n_sig_partners"], descending=[True, True])
    shared_part = active_hubs.filter(
        pl.col("hub_novelty_status") == HubNoveltyStatus.SHARED_HUB.value
    ).sort(["peak_neg_log10_fdr_iut", "n_sig_partners"], descending=[True, True])

    half_cap = max(config.max_vector_hub_overlays // 2, 5)
    n_dn_keep = min(denovo_part.height, max(half_cap, config.max_vector_hub_overlays - shared_part.height))
    n_sh_keep = min(shared_part.height, config.max_vector_hub_overlays - n_dn_keep)
    hub_overlay_df = pl.concat([denovo_part.head(n_dn_keep), shared_part.head(n_sh_keep)])

    # Project pairwise hubness & novelty status onto single cells for smooth 300 DPI manifold cloud
    bin_mat = (sp.csc_matrix(nhoods_mat) > 0).astype(np.float64)
    cell_nh_cnt = np.maximum(np.asarray(bin_mat.sum(axis=1)).ravel(), 1.0)
    cell_pw_peak = np.asarray(bin_mat.dot(peak_neglog10_q)).ravel() / cell_nh_cnt
    cell_dn_hits = np.asarray(
        bin_mat.dot((n_partners > 0).astype(np.float64) * (~univ_is_sig).astype(np.float64))
    ).ravel()
    cell_sh_hits = np.asarray(
        bin_mat.dot((n_partners > 0).astype(np.float64) * univ_is_sig.astype(np.float64))
    ).ravel()
    cell_uv_only = np.asarray(
        bin_mat.dot((n_partners == 0).astype(np.float64) * univ_is_sig.astype(np.float64))
    ).ravel()

    cell_hub_status: list[str] = []
    for i in range(n_cells):
        if cell_dn_hits[i] > 0.0:
            cell_hub_status.append(HubNoveltyStatus.DE_NOVO_HUB.value)
        elif cell_sh_hits[i] > 0.0:
            cell_hub_status.append(HubNoveltyStatus.SHARED_HUB.value)
        elif cell_uv_only[i] > 0.0:
            cell_hub_status.append(HubNoveltyStatus.UNIVARIATE_ONLY.value)
        else:
            cell_hub_status.append(HubNoveltyStatus.INACTIVE.value)

    coords = np.asarray(umap_coords, dtype=np.float64)
    full_cells_df = pl.DataFrame(
        {
            "cell_idx": np.arange(n_cells, dtype=np.int64),
            "umap_1": coords[:, 0],
            "umap_2": coords[:, 1],
            "cell_label": resolved_cell_labels,
            "cell_pairwise_peak_q": cell_pw_peak,
            "cell_hub_status": cell_hub_status,
        }
    )

    if n_cells > config.max_bg_cells:
        rng = np.random.default_rng(config.seed)
        active_c_mask = (cell_dn_hits + cell_sh_hits) > 0.0
        act_idx = np.flatnonzero(active_c_mask)
        bg_idx = np.flatnonzero(~active_c_mask)
        n_act_sample = min(int(act_idx.shape[0]), int(0.55 * config.max_bg_cells))
        n_bg_sample = min(int(bg_idx.shape[0]), config.max_bg_cells - n_act_sample)
        pick_act = (
            rng.choice(act_idx, size=n_act_sample, replace=False)
            if int(act_idx.shape[0]) > n_act_sample
            else act_idx
        )
        pick_bg = (
            rng.choice(bg_idx, size=n_bg_sample, replace=False)
            if int(bg_idx.shape[0]) > n_bg_sample
            else bg_idx
        )
        keep_idx = np.concatenate([np.sort(pick_bg), np.sort(pick_act)])
        bg_cells_df = full_cells_df[keep_idx.tolist()]
    else:
        bg_cells_df = full_cells_df

    # Sort so active/de-novo cells draw above inactive grey cells
    status_rank = {
        HubNoveltyStatus.INACTIVE.value: 0,
        HubNoveltyStatus.UNIVARIATE_ONLY.value: 1,
        HubNoveltyStatus.SHARED_HUB.value: 2,
        HubNoveltyStatus.DE_NOVO_HUB.value: 3,
    }
    bg_cells_df = (
        bg_cells_df.with_columns(
            pl.col("cell_hub_status").replace_strict(status_rank, default=0).alias("_draw_order")
        )
        .sort(["_draw_order", "cell_pairwise_peak_q"])
        .drop("_draw_order")
    )

    return Success((nhood_nodes_df, hub_overlay_df, bg_cells_df))


def _empty_arcs_and_endpoints() -> tuple[pl.DataFrame, pl.DataFrame]:
    """Returns typed empty arc-vertex and endpoint DataFrames when zero pairs are plotted."""
    empty_arcs = pl.DataFrame(
        schema={
            "arc_id": pl.String,
            "step_idx": pl.Int64,
            "umap_1": pl.Float64,
            "umap_2": pl.Float64,
            "nhood_u": pl.Int64,
            "nhood_v": pl.Int64,
            "regime": pl.String,
            "neg_log10_fdr_iut": pl.Float64,
            "stroke_width": pl.Float64,
            "is_exploratory": pl.Boolean,
        }
    )
    empty_ep = pl.DataFrame(
        schema={
            "nhood_idx": pl.Int64,
            "umap_1": pl.Float64,
            "umap_2": pl.Float64,
            "state_label": pl.String,
            "polarity": pl.String,
            "is_denovo": pl.Boolean,
            "neg_log10_fdr_iut": pl.Float64,
        }
    )
    return empty_arcs, empty_ep


def build_bundled_bezier_arcs(
    nhood_nodes_df: pl.DataFrame,
    pairwise_res: PairwiseAnalysisResult,
    umap_bounds: tuple[float, float, float, float],
    config: UmapPairVizConfig = UmapPairVizConfig(),
) -> Result[tuple[pl.DataFrame, pl.DataFrame], str]:
    """Builds regime-stratified, spatially deduplicated outward-bowing quadratic Bezier chords."""
    n_nhoods = nhood_nodes_df.height
    if n_nhoods == 0:
        return Success(_empty_arcs_and_endpoints())

    u1_lookup = nhood_nodes_df["umap_1"].to_numpy()
    u2_lookup = nhood_nodes_df["umap_2"].to_numpy()
    lbl_lookup = nhood_nodes_df["state_label"].to_list()
    univ_sig_lookup = (
        nhood_nodes_df["is_univariate_sig"].to_numpy()
        if "is_univariate_sig" in nhood_nodes_df.columns
        else np.zeros(n_nhoods, dtype=bool)
    )

    global_center = (float(np.median(u1_lookup)), float(np.median(u2_lookup)))
    span_diag = float(
        np.hypot(umap_bounds[1] - umap_bounds[0], umap_bounds[3] - umap_bounds[2])
    )
    min_sep = config.min_spatial_chord_sep_frac * max(span_diag, 1.0)

    # Collect candidate chords per visual regime: (nhood_u, nhood_v, regime, pol_u, pol_v, p_iut, fdr_iut, exploratory)
    candidates_by_regime: dict[
        PairwiseVisualRegime,
        list[tuple[int, int, PairwiseVisualRegime, EndpointPolarity, EndpointPolarity, float, float, bool]],
    ] = {reg: [] for reg in ALL_VISUAL_REGIMES}

    # 1. Quadrant ecosystem modules (already consolidated by Leiden CPM in group_nhood_pairs)
    used_quad_uv: set[tuple[int, int]] = set()
    if not pairwise_res.modules_df.is_empty():
        for r in pairwise_res.modules_df.sort("best_pvalue_iut", descending=False).to_dicts():
            u = int(r["lead_nhood_u"])
            v = int(r["lead_nhood_v"])
            if 0 <= u < n_nhoods and 0 <= v < n_nhoods and u != v:
                q_enum = QuadrantState(str(r["quadrant"]))
                d_enum = EnrichmentDirection(str(r["direction"]))
                reg, pu, pv = classify_quadrant_visual_regime_and_polarity(q_enum, d_enum)
                key = (min(u, v), max(u, v))
                used_quad_uv.add(key)
                candidates_by_regime[reg].append(
                    (
                        u,
                        v,
                        reg,
                        pu,
                        pv,
                        float(r["best_pvalue_iut"]),
                        float(r["best_fdr_iut"]),
                        False,
                    )
                )

    if not pairwise_res.pairs_df.is_empty():
        sig_q = pairwise_res.pairs_df.filter(pl.col("is_significant")).sort(
            "pvalue_iut", descending=False
        )
        for r in sig_q.to_dicts():
            u = int(r["nhood_u"])
            v = int(r["nhood_v"])
            key = (min(u, v), max(u, v))
            if 0 <= u < n_nhoods and 0 <= v < n_nhoods and u != v and key not in used_quad_uv:
                used_quad_uv.add(key)
                q_enum = QuadrantState(str(r["quadrant"]))
                d_enum = EnrichmentDirection(str(r["direction"]))
                reg, pu, pv = classify_quadrant_visual_regime_and_polarity(q_enum, d_enum)
                candidates_by_regime[reg].append(
                    (
                        u,
                        v,
                        reg,
                        pu,
                        pv,
                        float(r["pvalue_iut"]),
                        float(r["fdr_iut"]),
                        False,
                    )
                )

    # 2. Continuous balance & co-abundance hits
    if not pairwise_res.continuous_pairs_df.is_empty():
        sig_c = pairwise_res.continuous_pairs_df.filter(pl.col("is_significant")).sort(
            "pvalue_iut", descending=False
        )
        for r in sig_c.to_dicts():
            u = int(r["nhood_u"])
            v = int(r["nhood_v"])
            if 0 <= u < n_nhoods and 0 <= v < n_nhoods and u != v:
                ax_enum = ContinuousAxis(str(r["axis"]))
                d_enum = EnrichmentDirection(str(r["direction"]))
                reg, pu, pv = classify_continuous_visual_regime_and_polarity(ax_enum, d_enum)
                candidates_by_regime[reg].append(
                    (
                        u,
                        v,
                        reg,
                        pu,
                        pv,
                        float(r["pvalue_iut"]),
                        float(r["fdr_iut"]),
                        False,
                    )
                )

    total_sig_cands = sum(len(v) for v in candidates_by_regime.values())
    if total_sig_cands == 0 and config.allow_exploratory_fallback:
        if not pairwise_res.pairs_df.is_empty():
            for r in pairwise_res.pairs_df.sort("pvalue_iut", descending=False).head(12).to_dicts():
                u = int(r["nhood_u"])
                v = int(r["nhood_v"])
                if 0 <= u < n_nhoods and 0 <= v < n_nhoods and u != v:
                    q_enum = QuadrantState(str(r["quadrant"]))
                    d_enum = EnrichmentDirection(str(r["direction"]))
                    reg, pu, pv = classify_quadrant_visual_regime_and_polarity(q_enum, d_enum)
                    candidates_by_regime[reg].append(
                        (
                            u,
                            v,
                            reg,
                            pu,
                            pv,
                            float(r["pvalue_iut"]),
                            float(r["fdr_iut"]),
                            True,
                        )
                    )

    # Spatially deduplicate within each regime so parallel overlapping chords collapse cleanly
    deduped_by_regime: dict[
        PairwiseVisualRegime,
        list[tuple[int, int, PairwiseVisualRegime, EndpointPolarity, EndpointPolarity, float, float, bool]],
    ] = {}
    for reg in ALL_VISUAL_REGIMES:
        kept_reg: list[
            tuple[int, int, PairwiseVisualRegime, EndpointPolarity, EndpointPolarity, float, float, bool]
        ] = []
        kept_coords: list[tuple[float, float, float, float]] = []
        for cand in sorted(candidates_by_regime[reg], key=lambda item: (item[5], item[6])):
            u, v = cand[0], cand[1]
            xu, yu = float(u1_lookup[u]), float(u2_lookup[u])
            xv, yv = float(u1_lookup[v]), float(u2_lookup[v])
            too_close = False
            for oxu, oyu, oxv, oyv in kept_coords:
                d_dir = float(np.hypot(xu - oxu, yu - oyu) + np.hypot(xv - oxv, yv - oyv))
                d_swp = float(np.hypot(xu - oxv, yu - oyv) + np.hypot(xv - oxu, yv - oyu))
                if min(d_dir, d_swp) < min_sep:
                    too_close = True
                    break
            if not too_close:
                kept_reg.append(cand)
                kept_coords.append((xu, yu, xv, yv))
                if len(kept_reg) >= config.max_chords:
                    break
        if kept_reg:
            deduped_by_regime[reg] = kept_reg

    if not deduped_by_regime:
        return Success(_empty_arcs_and_endpoints())

    # Round-robin across active regimes so rare quadrant regimes are never drowned out by abundant axes
    selected: list[
        tuple[int, int, PairwiseVisualRegime, EndpointPolarity, EndpointPolarity, float, float, bool]
    ] = []
    seen_unordered: set[tuple[int, int]] = set()
    active_regs = [reg for reg in ALL_VISUAL_REGIMES if reg in deduped_by_regime]
    cursors = {reg: 0 for reg in active_regs}
    while len(selected) < config.max_chords:
        progressed = False
        for reg in active_regs:
            cur = cursors[reg]
            bucket = deduped_by_regime[reg]
            while cur < len(bucket):
                cand = bucket[cur]
                cur += 1
                uv_key = (min(cand[0], cand[1]), max(cand[0], cand[1]))
                if uv_key not in seen_unordered:
                    seen_unordered.add(uv_key)
                    selected.append(cand)
                    progressed = True
                    break
            cursors[reg] = cur
            if len(selected) >= config.max_chords:
                break
        if not progressed:
            break

    arc_rows: list[dict[str, object]] = []
    ep_best: dict[int, tuple[EndpointPolarity, bool, float]] = {}

    fan_offsets = (-0.045, 0.0, 0.045, -0.025, 0.025)
    for k_idx, (u, v, reg, pol_u, pol_v, _p_iut, fdr_iut, is_exp) in enumerate(selected):
        nlq = float(-np.log10(max(fdr_iut, 1e-12)))
        sw = float(np.clip(0.95 + 0.55 * nlq, 1.0, 3.0))
        curv = max(0.04, config.base_curvature + fan_offsets[k_idx % len(fan_offsets)])
        bx, by = discretize_outward_bezier_arc(
            p_u=(float(u1_lookup[u]), float(u2_lookup[u])),
            p_v=(float(u1_lookup[v]), float(u2_lookup[v])),
            global_center=global_center,
            curvature=curv,
            n_steps=config.bezier_steps,
        )
        arc_id = f"arc_{k_idx:03d}_{min(u, v)}_{max(u, v)}"
        for s_i in range(int(bx.shape[0])):
            arc_rows.append(
                {
                    "arc_id": arc_id,
                    "step_idx": s_i,
                    "umap_1": float(bx[s_i]),
                    "umap_2": float(by[s_i]),
                    "nhood_u": u,
                    "nhood_v": v,
                    "regime": reg.value,
                    "neg_log10_fdr_iut": nlq,
                    "stroke_width": sw,
                    "is_exploratory": is_exp,
                }
            )

        for ep_id, pol in ((u, pol_u), (v, pol_v)):
            is_dn = not bool(univ_sig_lookup[ep_id])
            if (ep_id not in ep_best) or (nlq > ep_best[ep_id][2]):
                ep_best[ep_id] = (pol, is_dn, nlq)

    ep_rows = [
        {
            "nhood_idx": ep_id,
            "umap_1": float(u1_lookup[ep_id]),
            "umap_2": float(u2_lookup[ep_id]),
            "state_label": str(lbl_lookup[ep_id]),
            "polarity": pol.value,
            "is_denovo": is_dn,
            "neg_log10_fdr_iut": nlq,
        }
        for ep_id, (pol, is_dn, nlq) in sorted(ep_best.items(), key=lambda kv: kv[1][2])
    ]

    return Success((pl.DataFrame(arc_rows), pl.DataFrame(ep_rows)))


def compute_anchor_partner_field(
    nhoods_mat: SparseMatrix,
    nhood_nodes_df: pl.DataFrame,
    bg_cells_df: pl.DataFrame,
    pairwise_res: PairwiseAnalysisResult,
    anchor_state: Maybe[str] = Nothing,
    anchor_nhood: Maybe[int] = Nothing,
    config: UmapPairVizConfig = UmapPairVizConfig(),
) -> Result[tuple[str, pl.DataFrame, pl.DataFrame, pl.DataFrame], str]:
    """Computes a seed-conditioned 1-to-many signed interaction field `Z_A(v)` across the UMAP."""
    n_nhoods = nhood_nodes_df.height
    if n_nhoods == 0:
        empty_arcs, _ = _empty_arcs_and_endpoints()
        return Success(("Unassigned", bg_cells_df, nhood_nodes_df.head(0), empty_arcs))

    state_list = [str(s) for s in nhood_nodes_df["state_label"].to_list()]
    deg_arr = nhood_nodes_df["n_sig_partners"].to_numpy()
    u1_lookup = nhood_nodes_df["umap_1"].to_numpy()
    u2_lookup = nhood_nodes_df["umap_2"].to_numpy()
    global_center = (float(np.median(u1_lookup)), float(np.median(u2_lookup)))

    # Resolve anchor set A
    anchor_indices: list[int] = []
    chosen_label = ""
    match anchor_nhood:
        case Some(nh_id):
            if 0 <= int(nh_id) < n_nhoods:
                anchor_indices = [int(nh_id)]
                chosen_label = f"{state_list[int(nh_id)]} (#{int(nh_id)})"
        case Maybe.empty:
            pass
        case _ as unreachable_an:
            _ = unreachable_an

    if not anchor_indices:
        match anchor_state:
            case Some(st_query):
                q_low = st_query.strip().lower()
                matched = [
                    idx
                    for idx, st in enumerate(state_list)
                    if q_low and (st.lower() == q_low or q_low in st.lower())
                ]
                if matched:
                    anchor_indices = matched
                    chosen_label = state_list[matched[0]]
            case Maybe.empty:
                pass
            case _ as unreachable_as:
                _ = unreachable_as

    if not anchor_indices:
        # Auto-select state with highest total significant partner degree (or peak degree if ties)
        state_deg_totals: dict[str, int] = {}
        for idx, st in enumerate(state_list):
            if st.lower() not in ("unspecified", "unannotated", "unknown"):
                state_deg_totals[st] = state_deg_totals.get(st, 0) + int(deg_arr[idx])
        if state_deg_totals and max(state_deg_totals.values()) > 0:
            best_st = sorted(
                state_deg_totals.items(),
                key=lambda kv: (-kv[1], kv[0]),
            )[0][0]
            anchor_indices = [idx for idx, st in enumerate(state_list) if st == best_st]
            chosen_label = best_st
        else:
            best_u = int(np.argmax(deg_arr)) if int(np.max(deg_arr)) > 0 else 0
            anchor_indices = [best_u]
            chosen_label = f"{state_list[best_u]} (#{best_u})"

    anchor_set = set(anchor_indices)
    partner_syn = np.zeros(n_nhoods, dtype=np.float64)
    partner_ant = np.zeros(n_nhoods, dtype=np.float64)
    best_edge_to_partner: dict[int, tuple[int, PairwiseVisualRegime, float, float]] = {}

    if not pairwise_res.pairs_df.is_empty():
        for r in pairwise_res.pairs_df.to_dicts():
            u = int(r["nhood_u"])
            v = int(r["nhood_v"])
            u_in, v_in = (u in anchor_set), (v in anchor_set)
            if u_in == v_in or not (0 <= u < n_nhoods and 0 <= v < n_nhoods):
                continue
            p_iut = float(r["pvalue_iut"])
            if p_iut >= 0.9999:
                continue
            anch_u = u if u_in else v
            part_v = v if u_in else u
            z_val = _p_to_z_mag(p_iut)
            q_enum = QuadrantState(str(r["quadrant"]))
            d_enum = EnrichmentDirection(str(r["direction"]))
            reg, _, _ = classify_quadrant_visual_regime_and_polarity(q_enum, d_enum)
            if reg == PairwiseVisualRegime.SYNERGY_PP or (
                q_enum == QuadrantState.BOTH_MISSING and d_enum == EnrichmentDirection.NON_RESPONDER
            ):
                partner_syn[part_v] = max(partner_syn[part_v], z_val)
            else:
                partner_ant[part_v] = max(partner_ant[part_v], z_val)
            if bool(r["is_significant"]):
                fdr_v = float(r["fdr_iut"])
                if (part_v not in best_edge_to_partner) or (
                    p_iut < best_edge_to_partner[part_v][2]
                ):
                    best_edge_to_partner[part_v] = (anch_u, reg, p_iut, fdr_v)

    if not pairwise_res.continuous_pairs_df.is_empty():
        for r in pairwise_res.continuous_pairs_df.to_dicts():
            u = int(r["nhood_u"])
            v = int(r["nhood_v"])
            u_in, v_in = (u in anchor_set), (v in anchor_set)
            if u_in == v_in or not (0 <= u < n_nhoods and 0 <= v < n_nhoods):
                continue
            p_iut = float(r["pvalue_iut"])
            if p_iut >= 0.9999:
                continue
            anch_u = u if u_in else v
            part_v = v if u_in else u
            z_val = _p_to_z_mag(p_iut)
            ax_enum = ContinuousAxis(str(r["axis"]))
            d_enum = EnrichmentDirection(str(r["direction"]))
            reg, _, _ = classify_continuous_visual_regime_and_polarity(ax_enum, d_enum)
            if ax_enum == ContinuousAxis.CO_ABUNDANCE_PRODUCT and d_enum == EnrichmentDirection.RESPONDER:
                partner_syn[part_v] = max(partner_syn[part_v], z_val)
            else:
                partner_ant[part_v] = max(partner_ant[part_v], 0.85 * z_val)
            if bool(r["is_significant"]):
                fdr_v = float(r["fdr_iut"])
                if (part_v not in best_edge_to_partner) or (
                    p_iut < best_edge_to_partner[part_v][2]
                ):
                    best_edge_to_partner[part_v] = (anch_u, reg, p_iut, fdr_v)

    partner_net_z = np.clip(partner_syn - partner_ant, -6.0, 6.0)
    is_anchor_arr = np.array([idx in anchor_set for idx in range(n_nhoods)], dtype=bool)

    # Project anchor field Z_A(v) onto the displayed background cells
    bin_mat = (sp.csc_matrix(nhoods_mat) > 0).astype(np.float64)
    bg_idx = bg_cells_df["cell_idx"].to_numpy()
    bin_sub = bin_mat[bg_idx, :]
    sub_cnt = np.maximum(np.asarray(bin_sub.sum(axis=1)).ravel(), 1.0)
    cell_anchor_z = np.asarray(bin_sub.dot(partner_net_z)).ravel() / sub_cnt
    cell_is_anchor = np.asarray(bin_sub.dot(is_anchor_arr.astype(np.float64))).ravel() > 0.0

    anchor_bg_df = (
        bg_cells_df.with_columns(
            pl.Series("anchor_partner_z", cell_anchor_z),
            pl.Series("is_anchor_cell", cell_is_anchor),
            pl.Series("abs_anchor_z", np.abs(cell_anchor_z)),
        )
        .sort("abs_anchor_z")
        .drop("abs_anchor_z")
    )

    # Build vector overlay of anchor medoids (up to 35) + top partner medoids
    anch_medoids = (
        nhood_nodes_df.with_columns(
            pl.Series("is_anchor", is_anchor_arr),
            pl.Series("anchor_partner_z", partner_net_z),
        )
        .filter(pl.col("is_anchor"))
        .sort("n_sig_partners", descending=True)
        .head(35)
    )
    top_partners_sorted = sorted(best_edge_to_partner.items(), key=lambda kv: (kv[1][2], kv[1][3]))[
        : config.max_anchor_chords
    ]
    part_ids = [kv[0] for kv in top_partners_sorted]
    part_medoids = (
        nhood_nodes_df.with_columns(
            pl.Series("is_anchor", is_anchor_arr),
            pl.Series("anchor_partner_z", partner_net_z),
        )
        .filter(pl.col("nhood_idx").is_in(part_ids))
    )
    anchor_overlay_df = pl.concat([part_medoids, anch_medoids])

    arc_rows: list[dict[str, object]] = []
    for k_i, (part_v, (anch_u, reg, _p_iut, fdr_v)) in enumerate(top_partners_sorted):
        nlq = float(-np.log10(max(fdr_v, 1e-12)))
        bx, by = discretize_outward_bezier_arc(
            p_u=(float(u1_lookup[anch_u]), float(u2_lookup[anch_u])),
            p_v=(float(u1_lookup[part_v]), float(u2_lookup[part_v])),
            global_center=global_center,
            curvature=config.base_curvature + 0.02 * ((k_i % 3) - 1),
            n_steps=config.bezier_steps,
        )
        arc_id = f"anch_arc_{k_i:02d}_{anch_u}_{part_v}"
        for s_i in range(int(bx.shape[0])):
            arc_rows.append(
                {
                    "arc_id": arc_id,
                    "step_idx": s_i,
                    "umap_1": float(bx[s_i]),
                    "umap_2": float(by[s_i]),
                    "nhood_u": anch_u,
                    "nhood_v": part_v,
                    "regime": reg.value,
                    "neg_log10_fdr_iut": nlq,
                    "stroke_width": float(np.clip(1.1 + 0.45 * nlq, 1.1, 2.8)),
                    "is_exploratory": False,
                }
            )

    anchor_arcs_df = pl.DataFrame(arc_rows) if arc_rows else _empty_arcs_and_endpoints()[0]
    return Success((chosen_label, anchor_bg_df, anchor_overlay_df, anchor_arcs_df))


def extract_lead_pair_spotlight(
    nhoods_mat: SparseMatrix,
    umap_coords: FloatArray,
    nhood_nodes_df: pl.DataFrame,
    bg_cells_df: pl.DataFrame,
    count_matrix: FloatArray,
    response_labels: Sequence[str],
    sample_library_sizes: Maybe[Sequence[float]],
    pairwise_res: PairwiseAnalysisResult,
    pair_override: Maybe[tuple[int, int]] = Nothing,
    config: UmapPairVizConfig = UmapPairVizConfig(),
) -> Result[tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame, SpotlightPairMetadata], str]:
    """Extracts the single-cell kNN clouds and 2x2 donor response inset for the #1 lead pair."""
    n_nhoods = nhood_nodes_df.height
    u1_lookup = nhood_nodes_df["umap_1"].to_numpy()
    u2_lookup = nhood_nodes_df["umap_2"].to_numpy()
    lbl_lookup = [str(s) for s in nhood_nodes_df["state_label"].to_list()]
    univ_fdr_lookup = nhood_nodes_df["univariate_fdr"].to_numpy()
    univ_sig_lookup = nhood_nodes_df["is_univariate_sig"].to_numpy()

    # Candidate selection: tuple of (priority_tier, p_iut, fdr_iut, u, v, reg, dir_s, eff_val, is_exp)
    chosen_tuple: Maybe[
        tuple[int, float, float, int, int, PairwiseVisualRegime, str, float, bool]
    ] = Nothing

    def _score_priority(u_id: int, v_id: int, is_quad: bool) -> int:
        both_ns = (not bool(univ_sig_lookup[u_id])) and (not bool(univ_sig_lookup[v_id]))
        cross_st = lbl_lookup[u_id] != lbl_lookup[v_id]
        if config.prefer_denovo_spotlight and both_ns and cross_st and is_quad:
            return 0
        if config.prefer_denovo_spotlight and both_ns and cross_st:
            return 1
        if cross_st and is_quad:
            return 2
        if is_quad:
            return 3
        if cross_st:
            return 4
        return 5

    cands: list[tuple[int, float, float, int, int, PairwiseVisualRegime, str, float, bool]] = []
    if not pairwise_res.pairs_df.is_empty():
        q_source = pairwise_res.pairs_df.filter(pl.col("is_significant"))
        is_exp_q = False
        if q_source.is_empty() and config.allow_exploratory_fallback:
            q_source = pairwise_res.pairs_df.sort("pvalue_iut").head(5)
            is_exp_q = True
        for r in q_source.to_dicts():
            u, v = int(r["nhood_u"]), int(r["nhood_v"])
            if 0 <= u < n_nhoods and 0 <= v < n_nhoods and u != v:
                if isinstance(pair_override, Some):
                    ov = pair_override.unwrap()
                    if (min(u, v), max(u, v)) != (min(ov[0], ov[1]), max(ov[0], ov[1])):
                        continue
                q_enum = QuadrantState(str(r["quadrant"]))
                d_enum = EnrichmentDirection(str(r["direction"]))
                reg, _, _ = classify_quadrant_visual_regime_and_polarity(q_enum, d_enum)
                cands.append(
                    (
                        _score_priority(u, v, True),
                        float(r["pvalue_iut"]),
                        float(r["fdr_iut"]),
                        u,
                        v,
                        reg,
                        f"{q_enum.value} ({d_enum.value})",
                        float(r["log_or_quadrant"]),
                        is_exp_q,
                    )
                )

    if not pairwise_res.continuous_pairs_df.is_empty():
        c_source = pairwise_res.continuous_pairs_df.filter(pl.col("is_significant"))
        is_exp_c = False
        if c_source.is_empty() and (not cands) and config.allow_exploratory_fallback:
            c_source = pairwise_res.continuous_pairs_df.sort("pvalue_iut").head(5)
            is_exp_c = True
        for r in c_source.to_dicts():
            u, v = int(r["nhood_u"]), int(r["nhood_v"])
            if 0 <= u < n_nhoods and 0 <= v < n_nhoods and u != v:
                if isinstance(pair_override, Some):
                    ov = pair_override.unwrap()
                    if (min(u, v), max(u, v)) != (min(ov[0], ov[1]), max(ov[0], ov[1])):
                        continue
                ax_enum = ContinuousAxis(str(r["axis"]))
                d_enum = EnrichmentDirection(str(r["direction"]))
                reg, _, _ = classify_continuous_visual_regime_and_polarity(ax_enum, d_enum)
                cands.append(
                    (
                        _score_priority(u, v, False),
                        float(r["pvalue_iut"]),
                        float(r["fdr_iut"]),
                        u,
                        v,
                        reg,
                        f"{ax_enum.value} ({d_enum.value})",
                        float(r["effect_size"]),
                        is_exp_c,
                    )
                )

    if cands:
        cands.sort(key=lambda item: (item[0], item[1], item[2]))
        chosen_tuple = Some(cands[0])

    empty_arcs, empty_ep = _empty_arcs_and_endpoints()
    if not isinstance(chosen_tuple, Some):
        spot_bg = bg_cells_df.select(["umap_1", "umap_2"]).with_columns(
            pl.lit("Background Manifold").alias("spotlight_membership")
        )
        empty_quad = pl.DataFrame(
            {
                "donor_quadrant": ["u+ / v+", "u+ / v-", "u- / v+", "u- / v-"],
                "quad_order": [0, 1, 2, 3],
                "n_responders": [0, 0, 0, 0],
                "n_total": [0, 0, 0, 0],
                "responder_pct": [0.0, 0.0, 0.0, 0.0],
                "ci_low_pct": [0.0, 0.0, 0.0, 0.0],
                "ci_high_pct": [0.0, 0.0, 0.0, 0.0],
            }
        )
        meta_empty = SpotlightPairMetadata(
            has_pair=False,
            is_exploratory=False,
            nhood_u=-1,
            nhood_v=-1,
            label_u="Unassigned",
            label_v="Unassigned",
            regime_label="No Significant Pairs",
            direction_label="FDR >= threshold",
            novelty_label="Unassigned",
            fdr_iut=1.0,
            pvalue_iut=1.0,
            effect_value=0.0,
            univariate_fdr_u=1.0,
            univariate_fdr_v=1.0,
            n_cells_u=0,
            n_cells_v=0,
        )
        return Success((spot_bg, empty_arcs, empty_ep, empty_quad, meta_empty))

    _, p_iut, fdr_iut, u_star, v_star, reg_star, dir_desc, eff_val, is_exp = chosen_tuple.unwrap()

    csc = sp.csc_matrix(nhoods_mat)
    rows_u = set(csc.indices[csc.indptr[u_star] : csc.indptr[u_star + 1]][
        csc.data[csc.indptr[u_star] : csc.indptr[u_star + 1]] > 0
    ].tolist())
    rows_v = set(csc.indices[csc.indptr[v_star] : csc.indptr[v_star + 1]][
        csc.data[csc.indptr[v_star] : csc.indptr[v_star + 1]] > 0
    ].tolist())

    def _compact_ontology(lbl: str) -> str:
        s = lbl.replace("CD8-positive, ", "CD8+ ").replace("CD4-positive, ", "CD4+ ")
        s = s.replace("CD28-negative, ", "CD28- ").replace("CD28-positive, ", "CD28+ ")
        s = s.replace("alpha-beta regulatory T cell", "ab Treg")
        s = s.replace("alpha-beta T cell", "ab T cell")
        s = s.replace("regulatory T cell", "Treg")
        return s

    lbl_u_short = f"u* (#{u_star}: {_compact_ontology(lbl_lookup[u_star])})"
    lbl_v_short = f"v* (#{v_star}: {_compact_ontology(lbl_lookup[v_star])})"

    coords = np.asarray(umap_coords, dtype=np.float64)
    bg_only_df = (
        bg_cells_df.filter(
            ~pl.col("cell_idx").is_in(sorted(rows_u | rows_v))
        )
        .select(["umap_1", "umap_2"])
        .with_columns(pl.lit("Background Manifold").alias("spotlight_membership"))
    )
    u_cells_df = pl.DataFrame(
        {
            "umap_1": coords[sorted(rows_u), 0] if rows_u else np.empty(0, dtype=np.float64),
            "umap_2": coords[sorted(rows_u), 1] if rows_u else np.empty(0, dtype=np.float64),
            "spotlight_membership": [lbl_u_short] * len(rows_u),
        }
    )
    v_cells_df = pl.DataFrame(
        {
            "umap_1": coords[sorted(rows_v), 0] if rows_v else np.empty(0, dtype=np.float64),
            "umap_2": coords[sorted(rows_v), 1] if rows_v else np.empty(0, dtype=np.float64),
            "spotlight_membership": [lbl_v_short] * len(rows_v),
        }
    )
    spotlight_cells_df = pl.concat([bg_only_df, u_cells_df, v_cells_df])

    global_center = (float(np.median(u1_lookup)), float(np.median(u2_lookup)))
    bx, by = discretize_outward_bezier_arc(
        p_u=(float(u1_lookup[u_star]), float(u2_lookup[u_star])),
        p_v=(float(u1_lookup[v_star]), float(u2_lookup[v_star])),
        global_center=global_center,
        curvature=config.base_curvature,
        n_steps=config.bezier_steps,
    )
    nlq = float(-np.log10(max(fdr_iut, 1e-12)))
    spotlight_arc_df = pl.DataFrame(
        {
            "arc_id": ["spotlight_arc"] * int(bx.shape[0]),
            "step_idx": list(range(int(bx.shape[0]))),
            "umap_1": bx,
            "umap_2": by,
            "nhood_u": [u_star] * int(bx.shape[0]),
            "nhood_v": [v_star] * int(bx.shape[0]),
            "regime": [reg_star.value] * int(bx.shape[0]),
            "neg_log10_fdr_iut": [nlq] * int(bx.shape[0]),
            "stroke_width": [2.4] * int(bx.shape[0]),
            "is_exploratory": [is_exp] * int(bx.shape[0]),
        }
    )
    spotlight_ep_df = pl.DataFrame(
        {
            "nhood_idx": [u_star, v_star],
            "umap_1": [float(u1_lookup[u_star]), float(u1_lookup[v_star])],
            "umap_2": [float(u2_lookup[u_star]), float(u2_lookup[v_star])],
            "state_label": [lbl_u_short, lbl_v_short],
            "polarity": [
                EndpointPolarity.POSITIVE_PRESENT.value,
                EndpointPolarity.NEGATIVE_GATED.value
                if reg_star in (PairwiseVisualRegime.GATING_DISCORDANT, PairwiseVisualRegime.BALANCE_RATIO)
                else EndpointPolarity.POSITIVE_PRESENT.value,
            ],
            "is_denovo": [not bool(univ_sig_lookup[u_star]), not bool(univ_sig_lookup[v_star])],
            "neg_log10_fdr_iut": [nlq, nlq],
        }
    )

    # 2x2 donor responder rate inset across (u+ / v+, u+ / v-, u- / v+, u- / v-)
    counts_sub = np.asarray(count_matrix[[u_star, v_star], :], dtype=np.float64)
    lib_maybe: Maybe[FloatArray] = (
        Some(np.asarray(sample_library_sizes.unwrap(), dtype=np.float64))
        if isinstance(sample_library_sizes, Some)
        else Nothing
    )
    pres_res = estimate_nhood_presence(
        count_matrix=counts_sub,
        library_sizes=lib_maybe,
        config=PresenceConfig(mode=PresenceMode.HARD_DEPTH_GATED),
    )
    match pres_res:
        case Success(pres_rep) if pres_rep.presence_prob.shape[0] == 2:
            u_pos = pres_rep.presence_prob[0, :] >= 0.5
            v_pos = pres_rep.presence_prob[1, :] >= 0.5
        case _:
            u_pos = counts_sub[0, :] > float(np.median(counts_sub[0, :]))
            v_pos = counts_sub[1, :] > float(np.median(counts_sub[1, :]))

    # If a continuous balance pair is ubiquitous (100% raw presence), split cleanly at marginal median CPM
    if int(np.sum(u_pos & v_pos)) == int(counts_sub.shape[1]) and reg_star in (
        PairwiseVisualRegime.BALANCE_RATIO,
        PairwiseVisualRegime.CO_ABUNDANCE,
    ):
        libs_f = (
            np.maximum(np.asarray(sample_library_sizes.unwrap(), dtype=np.float64), 1.0)
            if isinstance(sample_library_sizes, Some)
            else np.maximum(np.sum(count_matrix, axis=0), 1.0)
        )
        cpm_u = counts_sub[0, :] / libs_f
        cpm_v = counts_sub[1, :] / libs_f
        u_pos = cpm_u > float(np.median(cpm_u))
        v_pos = cpm_v > float(np.median(cpm_v))

    y_bin = np.array(
        [1 if str(r).strip().lower() in ("responder", "r", "1", "true", "cr", "pr") else 0 for r in response_labels],
        dtype=np.int64,
    )

    quad_defs = (
        ("u+ / v+", u_pos & v_pos),
        ("u+ / v-", u_pos & (~v_pos)),
        ("u- / v+", (~u_pos) & v_pos),
        ("u- / v-", (~u_pos) & (~v_pos)),
    )
    q_rows: list[dict[str, object]] = []
    for idx_q, (q_name, q_mask) in enumerate(quad_defs):
        n_tot = int(np.sum(q_mask))
        n_r = int(np.sum(y_bin[q_mask] == 1)) if n_tot > 0 else 0
        pct = 100.0 * float(n_r) / float(max(n_tot, 1))
        lo, hi = wilson_binomial_ci_pct(n_r, n_tot)
        q_rows.append(
            {
                "donor_quadrant": q_name,
                "quad_order": idx_q,
                "n_responders": n_r,
                "n_total": n_tot,
                "responder_pct": pct,
                "ci_low_pct": lo,
                "ci_high_pct": hi,
            }
        )

    both_ns = (not bool(univ_sig_lookup[u_star])) and (not bool(univ_sig_lookup[v_star]))
    either_ns = (not bool(univ_sig_lookup[u_star])) or (not bool(univ_sig_lookup[v_star]))
    nov_str = (
        "De Novo Pair (Both Univ. NS)"
        if both_ns
        else ("Semi-Novel (1 Univ. Sig)" if either_ns else "Both Univ. Sig")
    )

    meta = SpotlightPairMetadata(
        has_pair=True,
        is_exploratory=is_exp,
        nhood_u=u_star,
        nhood_v=v_star,
        label_u=lbl_u_short,
        label_v=lbl_v_short,
        regime_label=reg_star.value,
        direction_label=dir_desc,
        novelty_label=nov_str,
        fdr_iut=fdr_iut,
        pvalue_iut=p_iut,
        effect_value=eff_val,
        univariate_fdr_u=float(univ_fdr_lookup[u_star]),
        univariate_fdr_v=float(univ_fdr_lookup[v_star]),
        n_cells_u=len(rows_u),
        n_cells_v=len(rows_v),
    )
    return Success((spotlight_cells_df, spotlight_arc_df, spotlight_ep_df, pl.DataFrame(q_rows), meta))


def prepare_pairwise_umap_bundle(
    nhoods_mat: SparseMatrix,
    umap_coords: FloatArray,
    count_matrix: FloatArray,
    response_labels: Sequence[str],
    pairwise_res: PairwiseAnalysisResult,
    sample_library_sizes: Maybe[Sequence[float]] = Nothing,
    state_labels: Maybe[Sequence[str]] = Nothing,
    cell_labels: Maybe[Sequence[str]] = Nothing,
    univariate_fdr: Maybe[Sequence[float]] = Nothing,
    univariate_logfc: Maybe[Sequence[float]] = Nothing,
    anchor_state: Maybe[str] = Nothing,
    anchor_nhood: Maybe[int] = Nothing,
    spotlight_pair: Maybe[tuple[int, int]] = Nothing,
    config: UmapPairVizConfig = UmapPairVizConfig(),
) -> Result[PairwiseUmapBundle, str]:
    """Pure end-to-end orchestrator computing all DataFrames for the 4-panel pairwise UMAP composite."""
    med_res = compute_nhood_umap_medoids(
        nhoods_mat=nhoods_mat,
        umap_coords=umap_coords,
        state_labels=state_labels,
        method=config.coord_method,
    )
    match med_res:
        case Failure(err):
            return Failure(err)
        case Success(medoids_df):
            pass

    hub_res = summarize_nhood_pairwise_hubness(
        nhoods_mat=nhoods_mat,
        umap_coords=umap_coords,
        medoids_df=medoids_df,
        pairwise_res=pairwise_res,
        univariate_fdr=univariate_fdr,
        univariate_logfc=univariate_logfc,
        cell_labels=cell_labels,
        config=config,
    )
    match hub_res:
        case Failure(err):
            return Failure(err)
        case Success((nhood_nodes_df, hub_overlay_df, bg_cells_df)):
            pass

    bounds = compute_umap_bounding_domain(umap_coords)

    arcs_res = build_bundled_bezier_arcs(
        nhood_nodes_df=nhood_nodes_df,
        pairwise_res=pairwise_res,
        umap_bounds=bounds,
        config=config,
    )
    match arcs_res:
        case Failure(err):
            return Failure(err)
        case Success((bundled_arcs_df, arc_endpoints_df)):
            pass

    anch_res = compute_anchor_partner_field(
        nhoods_mat=nhoods_mat,
        nhood_nodes_df=nhood_nodes_df,
        bg_cells_df=bg_cells_df,
        pairwise_res=pairwise_res,
        anchor_state=anchor_state,
        anchor_nhood=anchor_nhood,
        config=config,
    )
    match anch_res:
        case Failure(err):
            return Failure(err)
        case Success((anchor_label, anchor_bg_df, anchor_overlay_df, anchor_arcs_df)):
            pass

    spot_res = extract_lead_pair_spotlight(
        nhoods_mat=nhoods_mat,
        umap_coords=umap_coords,
        nhood_nodes_df=nhood_nodes_df,
        bg_cells_df=bg_cells_df,
        count_matrix=count_matrix,
        response_labels=response_labels,
        sample_library_sizes=sample_library_sizes,
        pairwise_res=pairwise_res,
        pair_override=spotlight_pair,
        config=config,
    )
    match spot_res:
        case Failure(err):
            return Failure(err)
        case Success((spot_cells_df, spot_arc_df, spot_ep_df, spot_quad_df, spot_meta)):
            pass

    n_denovo = int(
        nhood_nodes_df.filter(
            pl.col("hub_novelty_status") == HubNoveltyStatus.DE_NOVO_HUB.value
        ).height
    )

    return Success(
        PairwiseUmapBundle(
            bg_cells_df=bg_cells_df,
            nhood_nodes_df=nhood_nodes_df,
            hub_overlay_df=hub_overlay_df,
            bundled_arcs_df=bundled_arcs_df,
            arc_endpoints_df=arc_endpoints_df,
            anchor_label=anchor_label,
            anchor_bg_cells_df=anchor_bg_df,
            anchor_overlay_df=anchor_overlay_df,
            anchor_arcs_df=anchor_arcs_df,
            spotlight_cells_df=spot_cells_df,
            spotlight_arc_df=spot_arc_df,
            spotlight_endpoints_df=spot_ep_df,
            spotlight_quadrant_df=spot_quad_df,
            spotlight_meta=spot_meta,
            umap_bounds=bounds,
            n_sig_quadrant=pairwise_res.n_significant_quadrant_pairs,
            n_sig_continuous=pairwise_res.n_significant_continuous_pairs,
            n_denovo_nhoods=n_denovo,
        )
    )


# =============================================================================
# 4. Pure Nature Methods Minimalist Wireframe Altair Chart Builders
# =============================================================================


def _umap_xy_encodings(
    bounds: tuple[float, float, float, float],
    show_x_title: bool = True,
    show_y_title: bool = True,
) -> tuple[alt.X, alt.Y]:
    """Shared padded UMAP coordinate axes adhering to Nature Methods wireframe styling."""
    x_enc = alt.X(
        "umap_1:Q",
        title="UMAP 1" if show_x_title else "",
        scale=alt.Scale(domain=[bounds[0], bounds[1]], nice=False),
        axis=alt.Axis(
            grid=False,
            tickCount=5,
            tickColor="#1E293B",
            domainColor=_HAIRLINE,
            labelColor=_TEXT_SECONDARY,
            titleColor=_TEXT_PRIMARY,
        ),
    )
    y_enc = alt.Y(
        "umap_2:Q",
        title="UMAP 2" if show_y_title else "",
        scale=alt.Scale(domain=[bounds[2], bounds[3]], nice=False),
        axis=alt.Axis(
            grid=False,
            tickCount=5,
            tickColor="#1E293B",
            domainColor=_HAIRLINE,
            labelColor=_TEXT_SECONDARY,
            titleColor=_TEXT_PRIMARY,
        ),
    )
    return x_enc, y_enc


def _apply_wireframe_theme(chart: AltairComposite) -> AltairComposite:
    """Applies strict Nature Methods pure-white wireframe theme (rx=0, #CBD5E1 border, no shadows)."""
    return cast(
        AltairComposite,
        chart.configure_view(
            strokeWidth=0.75,
            stroke=_HAIRLINE,
            fill=_BG_WHITE,
        )
        .configure_axis(
            labelFont="sans-serif",
            titleFont="sans-serif",
            labelFontSize=9.5,
            titleFontSize=10.5,
            labelColor=_TEXT_SECONDARY,
            titleColor=_TEXT_PRIMARY,
            gridColor=_LIGHT_GRID,
            domainColor=_HAIRLINE,
            tickColor="#1E293B",
        )
        .configure_legend(
            labelFont="sans-serif",
            titleFont="sans-serif",
            labelFontSize=9.0,
            titleFontSize=10.0,
            labelLimit=295,
            titleLimit=295,
            labelColor=_TEXT_SECONDARY,
            titleColor=_TEXT_PRIMARY,
        ),
    )


def build_chord_umap_panel(
    bundle: PairwiseUmapBundle,
    config: UmapPairVizConfig = UmapPairVizConfig(),
    panel_prefix: str = "a",
) -> alt.Chart | alt.LayerChart:
    """Panel a: Ecosystem-Bundled Outward-Bowing Bezier Chords & Polar Endpoint Glyphs."""
    x_enc, y_enc = _umap_xy_encodings(bundle.umap_bounds)
    bg_layer = (
        alt.Chart(bundle.bg_cells_df)
        .mark_circle(size=8, opacity=0.55, color=_CELL_INACTIVE)
        .encode(x=x_enc, y=y_enc)
    )

    regime_domain = [reg.value for reg in ALL_VISUAL_REGIMES]
    regime_colors = [regime_color_hex(reg) for reg in ALL_VISUAL_REGIMES]
    dash_domain = regime_domain
    dash_range: list[list[int]] = [
        [1, 0],
        [1, 0],
        [5, 3],
        [6, 2, 2, 2],
        [2, 2],
    ]

    layered: alt.Chart | alt.LayerChart = bg_layer
    n_arcs = (
        int(bundle.bundled_arcs_df["arc_id"].n_unique())
        if not bundle.bundled_arcs_df.is_empty()
        else 0
    )

    if n_arcs > 0:
        arc_layer = (
            alt.Chart(bundle.bundled_arcs_df)
            .mark_line(opacity=0.82, strokeWidth=1.65)
            .encode(
                x=x_enc,
                y=y_enc,
                detail="arc_id:N",
                order="step_idx:Q",
                color=alt.Color(
                    "regime:N",
                    scale=alt.Scale(domain=regime_domain, range=regime_colors),
                    title="Pairwise Regime",
                    legend=alt.Legend(orient="bottom", columns=2),
                ),
                strokeDash=alt.StrokeDash(
                    "regime:N",
                    scale=alt.Scale(domain=dash_domain, range=dash_range),
                    title="Pairwise Regime",
                    legend=alt.Legend(orient="bottom", columns=2),
                ),
            )
        )
        layered = layered + arc_layer

    if not bundle.arc_endpoints_df.is_empty():
        pos_df = bundle.arc_endpoints_df.filter(
            pl.col("polarity").is_in(
                [EndpointPolarity.POSITIVE_PRESENT.value, EndpointPolarity.SYMMETRIC_COABUND.value]
            )
        )
        neg_df = bundle.arc_endpoints_df.filter(
            pl.col("polarity") == EndpointPolarity.NEGATIVE_GATED.value
        )
        des_df = bundle.arc_endpoints_df.filter(
            pl.col("polarity") == EndpointPolarity.DUAL_ABSENT.value
        )
        if not pos_df.is_empty():
            ep_pos = (
                alt.Chart(pos_df)
                .mark_circle(size=46, opacity=0.95, color=_OKABE_BLUE, stroke=_BG_WHITE, strokeWidth=0.9)
                .encode(x=x_enc, y=y_enc)
            )
            layered = layered + ep_pos
        if not neg_df.is_empty():
            # Hollow white-core Vermilion target ring for excluded / gated (-) endpoints
            ep_neg = (
                alt.Chart(neg_df)
                .mark_point(
                    shape="circle",
                    filled=True,
                    size=52,
                    fill=_BG_WHITE,
                    stroke=_OKABE_VERMILION,
                    strokeWidth=1.9,
                    opacity=0.98,
                )
                .encode(x=x_enc, y=y_enc)
            )
            layered = layered + ep_neg
        if not des_df.is_empty():
            ep_des = (
                alt.Chart(des_df)
                .mark_point(
                    shape="diamond",
                    filled=True,
                    size=48,
                    fill=_OKABE_ORANGE,
                    stroke=_TEXT_PRIMARY,
                    strokeWidth=0.8,
                    opacity=0.95,
                )
                .encode(x=x_enc, y=y_enc)
            )
            layered = layered + ep_des

    sub = (
        f"Top {n_arcs} Non-Redundant Ecosystem Chords (Filled Blue: + | Hollow Red: - Gate)"
        if n_arcs > 0
        else "0 Significant Pairwise Chords at FDR < 0.10 (Strict Null Protection)"
    )
    return layered.properties(
        width=config.panel_width,
        height=config.panel_height,
        title=alt.TitleParams(
            text=f"{panel_prefix}. Bundled Pairwise Chord Architecture",
            subtitle=sub,
            fontSize=11.5,
            fontWeight="bold",
            subtitleFontSize=9.0,
            subtitleColor=_TEXT_SECONDARY,
            anchor="start",
        ),
    )


def build_denovo_hub_umap_panel(
    bundle: PairwiseUmapBundle,
    config: UmapPairVizConfig = UmapPairVizConfig(),
    panel_prefix: str = "b",
) -> alt.Chart | alt.LayerChart:
    """Panel b: Univariate vs. Emergent De Novo Pairwise Hubness Projection."""
    x_enc, y_enc = _umap_xy_encodings(bundle.umap_bounds)
    hub_domain = [
        HubNoveltyStatus.DE_NOVO_HUB.value,
        HubNoveltyStatus.SHARED_HUB.value,
        HubNoveltyStatus.UNIVARIATE_ONLY.value,
        HubNoveltyStatus.INACTIVE.value,
    ]
    hub_colors = [
        _OKABE_VERMILION,
        _OKABE_BLUE,
        _OKABE_SKY_BLUE,
        _CELL_INACTIVE,
    ]

    cell_layer = (
        alt.Chart(bundle.bg_cells_df)
        .mark_circle(size=9, opacity=0.78)
        .encode(
            x=x_enc,
            y=y_enc,
            color=alt.Color(
                "cell_hub_status:N",
                scale=alt.Scale(domain=hub_domain, range=hub_colors),
                title="Univariate vs. Pairwise Status",
                legend=alt.Legend(orient="bottom", columns=2),
            ),
        )
    )

    layered: alt.Chart | alt.LayerChart = cell_layer
    if not bundle.hub_overlay_df.is_empty():
        denovo_ov = bundle.hub_overlay_df.filter(
            pl.col("hub_novelty_status") == HubNoveltyStatus.DE_NOVO_HUB.value
        )
        if not denovo_ov.is_empty():
            outer_rings = (
                alt.Chart(denovo_ov)
                .mark_point(
                    shape="circle",
                    filled=False,
                    size=64,
                    stroke=_TEXT_PRIMARY,
                    strokeWidth=1.2,
                    opacity=0.90,
                )
                .encode(x=x_enc, y=y_enc)
            )
            inner_dots = (
                alt.Chart(denovo_ov)
                .mark_circle(size=24, color=_OKABE_VERMILION, opacity=0.95)
                .encode(x=x_enc, y=y_enc)
            )
            layered = layered + outer_rings + inner_dots

    sub = (
        f"{bundle.n_denovo_nhoods:,} De Novo Pairwise Hub Nhoods (Univariate Milo FDR >= {config.fdr_threshold:.2f})"
    )
    return layered.properties(
        width=config.panel_width,
        height=config.panel_height,
        title=alt.TitleParams(
            text=f"{panel_prefix}. Emergent De Novo Pairwise Hubs",
            subtitle=sub,
            fontSize=11.5,
            fontWeight="bold",
            subtitleFontSize=9.0,
            subtitleColor=_TEXT_SECONDARY,
            anchor="start",
        ),
    )


def build_anchor_partner_umap_panel(
    bundle: PairwiseUmapBundle,
    config: UmapPairVizConfig = UmapPairVizConfig(),
    panel_prefix: str = "c",
) -> alt.Chart | alt.LayerChart:
    """Panel c: Seed-Conditioned 1-to-Many Partner Interaction Field `Z_A(v)`."""
    x_enc, y_enc = _umap_xy_encodings(bundle.umap_bounds)
    z_vals = bundle.anchor_bg_cells_df["anchor_partner_z"].to_numpy()
    max_z = max(2.5, float(np.percentile(np.abs(z_vals), 99)) if z_vals.size > 0 else 2.5)
    cap_z = min(round(max_z, 1), 5.0)

    field_layer = (
        alt.Chart(bundle.anchor_bg_cells_df)
        .mark_circle(size=9, opacity=0.80)
        .encode(
            x=x_enc,
            y=y_enc,
            color=alt.Color(
                "anchor_partner_z:Q",
                scale=alt.Scale(
                    domain=[-cap_z, -0.5 * cap_z, 0.0, 0.5 * cap_z, cap_z],
                    range=[_OKABE_VERMILION, "#F4A582", _LIGHT_GRID, "#92C5DE", _OKABE_BLUE],
                    interpolate="rgb",
                    clamp=True,
                ),
                title="Partner Interaction Z (- Gate | + Synergy)",
                legend=alt.Legend(
                    orient="bottom",
                    direction="horizontal",
                    gradientLength=185,
                    titleLimit=295,
                ),
            ),
        )
    )

    layered: alt.Chart | alt.LayerChart = field_layer
    if not bundle.anchor_arcs_df.is_empty():
        rad_arcs = (
            alt.Chart(bundle.anchor_arcs_df)
            .mark_line(color=_TEXT_PRIMARY, opacity=0.58, strokeWidth=1.25, strokeDash=[4, 2])
            .encode(x=x_enc, y=y_enc, detail="arc_id:N", order="step_idx:Q")
        )
        layered = layered + rad_arcs

    if not bundle.anchor_overlay_df.is_empty():
        anch_only = bundle.anchor_overlay_df.filter(pl.col("is_anchor"))
        if not anch_only.is_empty():
            anch_marks = (
                alt.Chart(anch_only)
                .mark_point(
                    shape="diamond",
                    filled=True,
                    size=68,
                    fill=_OKABE_ORANGE,
                    stroke=_TEXT_PRIMARY,
                    strokeWidth=1.3,
                    opacity=0.98,
                )
                .encode(x=x_enc, y=y_enc)
            )
            layered = layered + anch_marks

    return layered.properties(
        width=config.panel_width,
        height=config.panel_height,
        title=alt.TitleParams(
            text=f"{panel_prefix}. Seed-Conditioned Partner Field",
            subtitle=f"Focal Anchor (Orange Diamonds): {bundle.anchor_label}",
            fontSize=11.5,
            fontWeight="bold",
            subtitleFontSize=9.0,
            subtitleColor=_TEXT_SECONDARY,
            anchor="start",
        ),
    )


def build_lead_pair_spotlight_panel(
    bundle: PairwiseUmapBundle,
    config: UmapPairVizConfig = UmapPairVizConfig(),
    panel_prefix: str = "d",
) -> alt.HConcatChart:
    """Panel d: Single-Cell kNN Footprint Spotlight + 2x2 Clinical Donor Quadrant Inset."""
    meta = bundle.spotlight_meta
    x_enc, y_enc = _umap_xy_encodings(bundle.umap_bounds)

    umap_w = max(int(config.panel_width * 0.60), 165)
    inset_w = max(config.panel_width - umap_w - 18, 95)

    if meta.has_pair:
        mem_domain = [meta.label_u, meta.label_v, "Background Manifold"]
        mem_colors = [_OKABE_BLUE, _OKABE_VERMILION, _CELL_INACTIVE]
    else:
        mem_domain = ["Background Manifold"]
        mem_colors = [_CELL_INACTIVE]

    spot_umap: alt.Chart | alt.LayerChart = (
        alt.Chart(bundle.spotlight_cells_df)
        .mark_circle(size=11, opacity=0.82)
        .encode(
            x=x_enc,
            y=y_enc,
            color=alt.Color(
                "spotlight_membership:N",
                scale=alt.Scale(domain=mem_domain, range=mem_colors),
                title="Lead Pair Cell Footprints",
                legend=alt.Legend(orient="bottom", direction="vertical"),
            ),
        )
    )

    if not bundle.spotlight_arc_df.is_empty():
        lead_arc = (
            alt.Chart(bundle.spotlight_arc_df)
            .mark_line(color=_TEXT_PRIMARY, strokeWidth=2.1, opacity=0.92)
            .encode(x=x_enc, y=y_enc, order="step_idx:Q")
        )
        spot_umap = spot_umap + lead_arc

    if not bundle.spotlight_endpoints_df.is_empty():
        ep_rings = (
            alt.Chart(bundle.spotlight_endpoints_df)
            .mark_point(
                shape="circle",
                filled=True,
                size=72,
                fill=_BG_WHITE,
                stroke=_TEXT_PRIMARY,
                strokeWidth=1.8,
            )
            .encode(x=x_enc, y=y_enc)
        )
        spot_umap = spot_umap + ep_rings

    sub_umap = (
        f"{meta.regime_label} | q_IUT = {meta.fdr_iut:.3g} ({meta.novelty_label})"
        if meta.has_pair
        else "No Significant Pair at Current Threshold"
    )
    left_view = spot_umap.properties(
        width=umap_w,
        height=config.panel_height,
        title=alt.TitleParams(
            text=f"{panel_prefix}. Lead Pair kNN Spotlight",
            subtitle=sub_umap,
            fontSize=11.5,
            fontWeight="bold",
            subtitleFontSize=8.8,
            subtitleColor=_TEXT_SECONDARY,
            anchor="start",
        ),
    )

    q_order = ["u+ / v+", "u+ / v-", "u- / v+", "u- / v-"]
    q_df = bundle.spotlight_quadrant_df
    non_empty_q = q_df.filter(pl.col("n_total") > 0)
    peak_q_name = (
        str(non_empty_q.sort(["responder_pct", "n_total"], descending=[True, True])["donor_quadrant"][0])
        if not non_empty_q.is_empty()
        else ""
    )
    q_x_enc = alt.X(
        "donor_quadrant:N",
        sort=q_order,
        title="Donor State",
        axis=alt.Axis(labelAngle=-28, tickColor="#1E293B", domainColor=_HAIRLINE),
    )
    q_y_enc = alt.Y(
        "responder_pct:Q",
        title="Responder %",
        scale=alt.Scale(domain=[0, 100]),
        axis=alt.Axis(tickCount=5, tickColor="#1E293B", domainColor=_HAIRLINE),
    )
    ref_bars = (
        alt.Chart(q_df.filter(pl.col("donor_quadrant") != peak_q_name))
        .mark_bar(width=18, cornerRadius=0, color="#94A3B8")
        .encode(x=q_x_enc, y=q_y_enc)
    )
    bars: alt.Chart | alt.LayerChart = ref_bars
    if peak_q_name:
        peak_bar = (
            alt.Chart(q_df.filter(pl.col("donor_quadrant") == peak_q_name))
            .mark_bar(width=18, cornerRadius=0, color=_OKABE_BLUE)
            .encode(x=q_x_enc, y=q_y_enc)
        )
        bars = bars + peak_bar
    ci_stem = (
        alt.Chart(bundle.spotlight_quadrant_df)
        .mark_rule(color=_TEXT_PRIMARY, strokeWidth=1.0)
        .encode(
            x=alt.X("donor_quadrant:N", sort=q_order),
            y="ci_low_pct:Q",
            y2="ci_high_pct:Q",
        )
    )
    ci_lo = (
        alt.Chart(bundle.spotlight_quadrant_df)
        .mark_tick(color=_TEXT_PRIMARY, thickness=1.0, size=8)
        .encode(x=alt.X("donor_quadrant:N", sort=q_order), y="ci_low_pct:Q")
    )
    ci_hi = (
        alt.Chart(bundle.spotlight_quadrant_df)
        .mark_tick(color=_TEXT_PRIMARY, thickness=1.0, size=8)
        .encode(x=alt.X("donor_quadrant:N", sort=q_order), y="ci_high_pct:Q")
    )
    right_inset = (bars + ci_stem + ci_lo + ci_hi).properties(
        width=inset_w,
        height=config.panel_height,
        title=alt.TitleParams(
            text="2x2 Biopsy Response",
            subtitle=f"u: q={meta.univariate_fdr_u:.2f}, v: q={meta.univariate_fdr_v:.2f}"
            if meta.has_pair
            else "Wilson 95% CI",
            fontSize=10.5,
            fontWeight="bold",
            subtitleFontSize=8.5,
            subtitleColor=_TEXT_SECONDARY,
            anchor="start",
        ),
    )

    return cast(
        alt.HConcatChart,
        alt.hconcat(left_view, right_inset, spacing=14).resolve_scale(color="independent"),
    )


def build_pairwise_umap_composite(
    bundle: PairwiseUmapBundle,
    title: str,
    subtitle: str,
    config: UmapPairVizConfig = UmapPairVizConfig(),
) -> AltairComposite:
    """Assembles the publication-grade 2x2 pairwise UMAP dashboard (`a`, `b`, `c`, `d`)."""
    panel_a = build_chord_umap_panel(bundle, config=config, panel_prefix="a")
    panel_b = build_denovo_hub_umap_panel(bundle, config=config, panel_prefix="b")
    panel_c = build_anchor_partner_umap_panel(bundle, config=config, panel_prefix="c")
    panel_d = build_lead_pair_spotlight_panel(bundle, config=config, panel_prefix="d")

    row_top = alt.hconcat(panel_a, panel_b, spacing=24).resolve_scale(
        color="independent",
        strokeDash="independent",
    )
    row_bot = alt.hconcat(panel_c, panel_d, spacing=24).resolve_scale(
        color="independent",
        strokeDash="independent",
    )

    composite = (
        alt.vconcat(row_top, row_bot, spacing=24)
        .resolve_scale(color="independent", strokeDash="independent")
        .properties(
            title=alt.TitleParams(
                text=title,
                subtitle=subtitle,
                fontSize=14.0,
                fontWeight="bold",
                subtitleFontSize=10.5,
                subtitleColor=_TEXT_SECONDARY,
                anchor="start",
                offset=12,
            )
        )
    )
    return _apply_wireframe_theme(composite)


def build_cohort_grid_umap_cell(
    bundle: PairwiseUmapBundle,
    cohort_title: str,
    cohort_subtitle: str,
    cell_width: int = 235,
    cell_height: int = 195,
) -> alt.Chart | alt.LayerChart:
    """Builds a single compact cohort UMAP cell for multi-cohort landscape comparison grids."""
    x_enc, y_enc = _umap_xy_encodings(bundle.umap_bounds)
    bg_layer = (
        alt.Chart(bundle.bg_cells_df)
        .mark_circle(size=5, opacity=0.52, color=_CELL_INACTIVE)
        .encode(x=x_enc, y=y_enc)
    )

    regime_domain = [reg.value for reg in ALL_VISUAL_REGIMES]
    regime_colors = [regime_color_hex(reg) for reg in ALL_VISUAL_REGIMES]
    dash_range: list[list[int]] = [
        [1, 0],
        [1, 0],
        [5, 3],
        [6, 2, 2, 2],
        [2, 2],
    ]

    layered: alt.Chart | alt.LayerChart = bg_layer
    if not bundle.hub_overlay_df.is_empty():
        dn_ov = bundle.hub_overlay_df.filter(
            pl.col("hub_novelty_status") == HubNoveltyStatus.DE_NOVO_HUB.value
        ).head(45)
        if not dn_ov.is_empty():
            dn_rings = (
                alt.Chart(dn_ov)
                .mark_point(
                    shape="circle",
                    filled=False,
                    size=36,
                    stroke="#64748B",
                    strokeWidth=0.85,
                    opacity=0.65,
                )
                .encode(x=x_enc, y=y_enc)
            )
            layered = layered + dn_rings

    # Always include arc_layer (even when empty) so shared legend resolves cleanly across 3x3 grid
    arc_layer = (
        alt.Chart(bundle.bundled_arcs_df)
        .mark_line(opacity=0.84, strokeWidth=1.45)
        .encode(
            x=x_enc,
            y=y_enc,
            detail="arc_id:N",
            order="step_idx:Q",
            color=alt.Color(
                "regime:N",
                scale=alt.Scale(domain=regime_domain, range=regime_colors),
                title="Pairwise Interaction Regime",
                legend=alt.Legend(orient="bottom", columns=3),
            ),
            strokeDash=alt.StrokeDash(
                "regime:N",
                scale=alt.Scale(domain=regime_domain, range=dash_range),
                title="Pairwise Interaction Regime",
                legend=alt.Legend(orient="bottom", columns=3),
            ),
        )
    )
    layered = layered + arc_layer

    if not bundle.arc_endpoints_df.is_empty():
        pos_df = bundle.arc_endpoints_df.filter(
            pl.col("polarity").is_in(
                [EndpointPolarity.POSITIVE_PRESENT.value, EndpointPolarity.SYMMETRIC_COABUND.value]
            )
        )
        neg_df = bundle.arc_endpoints_df.filter(
            pl.col("polarity") == EndpointPolarity.NEGATIVE_GATED.value
        )
        des_df = bundle.arc_endpoints_df.filter(
            pl.col("polarity") == EndpointPolarity.DUAL_ABSENT.value
        )
        if not pos_df.is_empty():
            layered = layered + (
                alt.Chart(pos_df)
                .mark_circle(size=32, opacity=0.95, color=_OKABE_BLUE, stroke=_BG_WHITE, strokeWidth=0.75)
                .encode(x=x_enc, y=y_enc)
            )
        if not neg_df.is_empty():
            layered = layered + (
                alt.Chart(neg_df)
                .mark_point(
                    shape="circle",
                    filled=True,
                    size=36,
                    fill=_BG_WHITE,
                    stroke=_OKABE_VERMILION,
                    strokeWidth=1.6,
                    opacity=0.98,
                )
                .encode(x=x_enc, y=y_enc)
            )
        if not des_df.is_empty():
            layered = layered + (
                alt.Chart(des_df)
                .mark_point(
                    shape="diamond",
                    filled=True,
                    size=34,
                    fill=_OKABE_ORANGE,
                    stroke=_TEXT_PRIMARY,
                    strokeWidth=0.7,
                    opacity=0.95,
                )
                .encode(x=x_enc, y=y_enc)
            )

    return layered.properties(
        width=cell_width,
        height=cell_height,
        title=alt.TitleParams(
            text=cohort_title,
            subtitle=cohort_subtitle,
            fontSize=10.5,
            fontWeight="bold",
            subtitleFontSize=8.5,
            subtitleColor=_TEXT_SECONDARY,
            anchor="start",
        ),
    )


def build_cross_cohort_umap_grid(
    cells: Sequence[tuple[PairwiseUmapBundle, str, str]],
    title: str,
    subtitle: str,
    n_cols: int = 3,
    cell_width: int = 235,
    cell_height: int = 195,
) -> Result[AltairComposite, str]:
    """Assembles an `R x C` landscape UMAP grid across cohorts with a single shared regime legend."""
    if not cells:
        return Failure("Cannot build cross-cohort UMAP grid from empty sequence")
    cols_eff = max(int(n_cols), 1)
    rendered_cells = [
        build_cohort_grid_umap_cell(
            bundle=b,
            cohort_title=ct,
            cohort_subtitle=cs,
            cell_width=cell_width,
            cell_height=cell_height,
        )
        for b, ct, cs in cells
    ]
    rows: list[AltairComposite] = []
    for start_i in range(0, len(rendered_cells), cols_eff):
        chunk = rendered_cells[start_i : start_i + cols_eff]
        row_chart: AltairComposite = (
            chunk[0] if len(chunk) == 1 else alt.hconcat(*chunk, spacing=18)
        )
        rows.append(row_chart)

    grid: AltairComposite = rows[0] if len(rows) == 1 else alt.vconcat(*rows, spacing=20)
    titled = grid.properties(
        title=alt.TitleParams(
            text=title,
            subtitle=subtitle,
            fontSize=13.5,
            fontWeight="bold",
            subtitleFontSize=10.0,
            subtitleColor=_TEXT_SECONDARY,
            anchor="start",
            offset=12,
        )
    )
    return Success(_apply_wireframe_theme(titled))
