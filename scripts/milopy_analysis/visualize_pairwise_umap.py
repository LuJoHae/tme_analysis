#!/usr/bin/env python3
"""Publication-Grade Pairwise Milopy UMAP Visualizer Across All 9 Clinical Immunotherapy Cohorts.

Renders four complementary Nature Methods vector UMAP panels (plus a 3x3 9-cohort landscape grid)
using the exact pre-computed single-cell `(UMAP1, UMAP2)` coordinates and Milo kNN neighborhood
graphs from `output/results/milopy_response_analysis/`:

1. Panel a — Ecosystem-Bundled Outward-Bowing Quadratic Bezier Chords & Polar Endpoint Glyphs:
   Top-K non-redundant module chords stratified across all five mechanistic regimes (`++` synergy,
   `+-`/`-+` mutual-exclusion gate, `--` immune desert, continuous log-ratio balance, continuous
   co-abundance coupling), anchored strictly at manifold-faithful geometric medoid cells.
2. Panel b — Emergent De Novo Pairwise Hubness Projection:
   Isolates neighborhoods invisible to 1D univariate Milo (`FDR >= 0.10`) that emerge strongly in
   2D pairwise testing (`FDR_IUT < 0.10`), marked with concentric `#0F172A` vector rings.
3. Panel c — Seed-Conditioned One-to-Many Partner Field `Z_A(v)`:
   Fixes the cohort's primary interaction hub (or user-pinned `--anchor-state` / `--anchor-nhood`)
   and projects its signed partner interaction field across the entire UMAP.
4. Panel d — Lead-Pair Single-Cell kNN Cloud Spotlight + 2x2 Clinical Biopsy Response Inset:
   Displays exact disjoint single-cell footprints `C_{u*}` and `C_{v*}` alongside donor Responder %
   (`+-` exact Wilson 95% CI) across `(u+ / v+, u+ / v-, u- / v+, u- / v-)`.

Hybrid Vector-Raster Export:
   Dense 12,000-cell UMAP manifold clouds (`>= 250` marks per group) are rasterized at 300 DPI
   inside native Inkscape SVG layers (`01_Background`, `02_Pairwise_UMAP_Panels`), while 100% of
   Bezier chord paths, polar endpoint glyphs, concentric de novo rings, axes, titles, and legends
   remain crisp editable vector SVG primitives.
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from enum import StrEnum
import html
from pathlib import Path
import sys
import time
from typing import Protocol, TypeAlias, assert_never
import xml.etree.ElementTree as ET

import altair as alt
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import sparse as sp
import vl_convert as vlc

from milopy.pairwise import (
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
    RefinementConfig,
    ScreenConfig,
    UmapPairVizConfig,
    analyze_nhood_pairs,
    build_anchor_partner_umap_panel,
    build_chord_umap_panel,
    build_cross_cohort_umap_grid,
    build_denovo_hub_umap_panel,
    build_lead_pair_spotlight_panel,
    build_pairwise_umap_composite,
    coarse_grain_nhoods_by_graph,
    deduplicate_continuous_pairs_to_lead_bundles,
    deduplicate_pairs_to_lead_bundles,
    extract_inferred_substudy_strata,
    prepare_pairwise_umap_bundle,
)
from milopy.pairwise.types import FloatArray, IntArray, SparseMatrix
from plotting_utils import (  # type: ignore[import-untyped]
    COHORT_PANELS,
    apply_nature_methods_theme,
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

DEFAULT_DISCOVERY_SPOTLIGHT_COHORTS: tuple[str, ...] = (
    "CELLxGENE_7b20c613",
    "GSE243013",
)


# =============================================================================
# 1. Exhaustive Sum Types, Frozen Domain Models & Typeclass Protocols
# =============================================================================


class PairRepresentationMode(StrEnum):
    """Pairwise hypothesis representation rendered on the UMAP."""

    BOTH = "both"
    SUPER_META42 = "super_meta42"
    RAW_MODULES = "raw_modules"


class ActivePairView(StrEnum):
    """Concrete pairwise view evaluated for a single composite figure."""

    SUPER_META42 = "super_meta42"
    RAW_MODULES = "raw_modules"


def expand_pair_views(mode: PairRepresentationMode) -> tuple[ActivePairView, ...]:
    """Expands CLI `--pair-mode` into concrete views via exhaustive pattern matching."""
    match mode:
        case PairRepresentationMode.BOTH:
            return (ActivePairView.SUPER_META42, ActivePairView.RAW_MODULES)
        case PairRepresentationMode.SUPER_META42:
            return (ActivePairView.SUPER_META42,)
        case PairRepresentationMode.RAW_MODULES:
            return (ActivePairView.RAW_MODULES,)
        case _ as unreachable:
            assert_never(unreachable)


class PairwiseUmapCliConfig(BaseModel):
    """Immutable CLI configuration for the 9-cohort pairwise UMAP pipeline."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    milo_response_dir: Path = Field(default=Path("output/results/milopy_response_analysis"))
    pairwise_t0_dir: Path = Field(default=Path("output/results/milopy_pairwise_analysis"))
    meta42_annotations_path: Path = Field(
        default=Path("output/results/milopy_pairwise_deconv_bridge/nhood_meta42_annotations.parquet")
    )
    output_dir: Path = Field(default=Path("output/results/milopy_pairwise_umap"))
    pair_mode: PairRepresentationMode = Field(default=PairRepresentationMode.BOTH)
    scope_filter: PairScopeFilter = Field(default=PairScopeFilter.CROSS_COMPARTMENT_ONLY)
    enable_mixing_qc: bool = Field(default=True)
    min_distinct_donors: int = Field(default=3, ge=1)
    min_effective_lisi: float = Field(default=2.5, ge=1.0)
    max_donor_dominance: float = Field(default=0.50, gt=0.0, le=1.0)
    use_substudy_strata: bool = Field(default=True)
    consolidate_state_pairs: bool = Field(default=True)
    cohort_filter: Maybe[tuple[str, ...]] = Field(default=Nothing)
    anchor_state: Maybe[str] = Field(default=Nothing)
    anchor_nhood: Maybe[int] = Field(default=Nothing)
    max_chords: int = Field(default=32, ge=4, le=160)
    max_bg_cells: int = Field(default=12000, ge=500, le=80000)
    fdr_threshold: float = Field(default=0.10, gt=0.0, le=1.0)
    tarone_p_ceiling: float = Field(default=0.005, gt=0.0, le=1.0)


def pair_view_display_title(
    view: ActivePairView,
    cli_cfg: PairwiseUmapCliConfig = PairwiseUmapCliConfig(),
) -> str:
    """Human-readable subtitle tag for each concrete pairwise representation."""
    match view:
        case ActivePairView.SUPER_META42:
            scope_tag = (
                "Cross-Compartment"
                if cli_cfg.scope_filter == PairScopeFilter.CROSS_COMPARTMENT_ONLY
                else f"Scope={cli_cfg.scope_filter.value}"
            )
            qc_tag = (
                f" + Donor-Mixing QC (LISI_pt>={cli_cfg.min_effective_lisi:.1f})"
                if cli_cfg.enable_mixing_qc
                else ""
            )
            cmh_tag = " + Sub-Study CMH" if cli_cfg.use_substudy_strata else ""
            dedup_tag = " + 1 Chord/State-Pair" if cli_cfg.consolidate_state_pairs else ""
            return f"Consensus-Meta42 Super-Nhoods ({scope_tag}{cmh_tag}{qc_tag}{dedup_tag})"
        case ActivePairView.RAW_MODULES:
            return "Unfiltered Disjoint Milo kNN Ecosystem Modules (T0)"
        case _ as unreachable:
            assert_never(unreachable)


class CohortRawVizArtifacts(BaseModel):
    """Immutable on-disk input bundle for a single clinical single-cell cohort."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    accession: str
    nhoods_mat: SparseMatrix
    cells_df: pl.DataFrame
    da_results_df: pl.DataFrame
    t0_quadrant_df: pl.DataFrame
    t0_continuous_df: pl.DataFrame
    t0_modules_df: pl.DataFrame
    meta42_sub_df: pl.DataFrame


class CohortDesignContext(BaseModel):
    """Immutable donor/lesion design matrices aligned to single cells and neighborhoods."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    accession: str
    cohort_label: str
    umap_coords: FloatArray
    cell_labels: tuple[str, ...]
    eval_ids: tuple[str, ...]
    inferred_strata: Maybe[tuple[str, ...]]
    response_labels: tuple[str, ...]
    count_matrix: FloatArray
    library_sizes: FloatArray
    cell_eval_codes: IntArray
    n_responders: int
    n_non_responders: int


class RenderedCohortViewSummary(BaseModel):
    """Summary statistics recorded per `(cohort, view)` for the Markdown report."""

    model_config = ConfigDict(frozen=True)

    accession: str
    cohort_label: str
    view: ActivePairView
    n_evals: int
    n_cells: int
    n_units: int
    n_sig_quadrant: int
    n_sig_continuous: int
    n_denovo_units: int
    n_displayed_chords: int
    anchor_label: str
    lead_pair_desc: str
    lead_fdr_iut: float
    exported_composite: bool


class CohortUmapArtifactLoader(Protocol):
    """Typeclass boundary for loading cached cohort UMAP, Milo, and pairwise Parquet artifacts."""

    def __call__(
        self,
        accession: str,
        config: PairwiseUmapCliConfig,
        all_meta42_df: pl.DataFrame,
    ) -> Result[CohortRawVizArtifacts, str]: ...


class HybridVectorFigureExporter(Protocol):
    """Typeclass boundary for exporting hybrid Inkscape SVG + 300 DPI PNG figures."""

    def __call__(
        self,
        chart: AltairComposite,
        base_path: Path,
        figure_title: str,
    ) -> Result[tuple[Path, Path], str]: ...


# =============================================================================
# 2. Pure Functional Core — Cohort Design, T0/T4 Pairwise Bundles & Inkscape XML
# =============================================================================


def resolve_cohorts_to_process(
    cohort_filter: Maybe[tuple[str, ...]],
) -> tuple[tuple[str, ...], frozenset[str]]:
    """Resolves `(all_grid_cohorts, full_composite_cohorts)` cleanly."""
    canonical = tuple(str(p["acc"]) for p in COHORT_PANELS)
    match cohort_filter:
        case Maybe.empty:
            return canonical, frozenset(DEFAULT_DISCOVERY_SPOTLIGHT_COHORTS)
        case Some(selected):
            cleaned = tuple(c.strip() for c in selected if c.strip())
            active = cleaned if cleaned else canonical
            return active, frozenset(active)
        case _ as unreachable:
            _ = unreachable
            return canonical, frozenset(DEFAULT_DISCOVERY_SPOTLIGHT_COHORTS)


def extract_cohort_design_context(
    raw: CohortRawVizArtifacts,
) -> Result[CohortDesignContext, str]:
    """Builds `(UMAP1, UMAP2)` arrays and donor/lesion count matrix `C` purely from `cells_df`."""
    req_cols = {"patient_id", "clinical_response", "UMAP1", "UMAP2"}
    missing = sorted(req_cols - set(raw.cells_df.columns))
    if missing:
        return Failure(f"[{raw.accession}] cell_level_scores missing columns: {missing}")

    n_cells = raw.cells_df.height
    mat_rows, n_nhoods = int(raw.nhoods_mat.shape[0]), int(raw.nhoods_mat.shape[1])
    if n_cells != mat_rows:
        return Failure(f"[{raw.accession}] cells ({n_cells}) != nhoods_mat rows ({mat_rows})")
    if raw.da_results_df.height != n_nhoods:
        return Failure(
            f"[{raw.accession}] da_results ({raw.da_results_df.height}) != nhoods ({n_nhoods})"
        )

    umap_coords = np.column_stack(
        [
            raw.cells_df["UMAP1"].cast(pl.Float64).to_numpy(),
            raw.cells_df["UMAP2"].cast(pl.Float64).to_numpy(),
        ]
    )

    norm_cells = (
        raw.cells_df.lazy()
        .with_columns(
            pl.col("patient_id").cast(pl.String).str.strip_chars().alias("_pt"),
            pl.col("clinical_response")
            .cast(pl.String)
            .str.strip_chars()
            .str.to_lowercase()
            .alias("_resp"),
            (
                pl.col("cell_type").cast(pl.String)
                if "cell_type" in raw.cells_df.columns
                else pl.lit("Cell")
            ).alias("_ct"),
        )
        .with_columns((pl.col("_pt") + pl.lit("__") + pl.col("_resp")).alias("_eval_id"))
        .collect()
    )

    eval_meta = (
        norm_cells.lazy()
        .group_by("_eval_id")
        .agg(
            pl.col("_resp").first().alias("clinical_response"),
            pl.len().cast(pl.Float64).alias("library_size"),
        )
        .sort("_eval_id")
        .collect()
    )

    eval_ids = tuple(str(e) for e in eval_meta["_eval_id"].to_list())
    n_evals = len(eval_ids)
    if n_evals < 4:
        return Failure(f"[{raw.accession}] Insufficient evaluations ({n_evals} < 4)")

    eval_to_idx = {eid: idx for idx, eid in enumerate(eval_ids)}
    cell_eval_codes = np.array(
        [eval_to_idx[str(e)] for e in norm_cells["_eval_id"].to_list()],
        dtype=np.int64,
    )
    sample_ind = sp.csr_matrix(
        (np.ones(n_cells, dtype=np.float64), (np.arange(n_cells, dtype=np.int64), cell_eval_codes)),
        shape=(n_cells, n_evals),
    )
    nhoods_csc = sp.csc_matrix(raw.nhoods_mat)
    count_matrix = np.asarray(nhoods_csc.T.dot(sample_ind).toarray(), dtype=np.float64)
    lib_sizes = eval_meta["library_size"].to_numpy().astype(np.float64)
    resp_labels = tuple(str(r) for r in eval_meta["clinical_response"].to_list())
    n_resp = sum(1 for r in resp_labels if r == "responder")
    cell_labels = tuple(str(c) for c in norm_cells["_ct"].to_list())
    inferred_strata = extract_inferred_substudy_strata(eval_ids)

    return Success(
        CohortDesignContext(
            accession=raw.accession,
            cohort_label=get_single_cell_cohort_display_label(raw.accession, "cancer_author_year"),
            umap_coords=umap_coords,
            cell_labels=cell_labels,
            eval_ids=eval_ids,
            inferred_strata=inferred_strata,
            response_labels=resp_labels,
            count_matrix=count_matrix,
            library_sizes=lib_sizes,
            cell_eval_codes=cell_eval_codes,
            n_responders=n_resp,
            n_non_responders=n_evals - n_resp,
        )
    )


def _resolve_raw_nhood_labels(
    da_results_df: pl.DataFrame,
    meta42_sub_df: pl.DataFrame,
) -> tuple[tuple[str, ...], tuple[str, ...], Mapping[str, str]]:
    """Returns `(raw_readable_labels, meta_state_codes, meta_code_to_readable_label)`."""
    n_nh = da_results_df.height
    raw_cts = (
        [str(x) for x in da_results_df["Nhood_CellType"].to_list()]
        if "Nhood_CellType" in da_results_df.columns
        else ["Unspecified"] * n_nh
    )

    meta_codes = ["Unassigned"] * n_nh
    meta_lbls = ["Unassigned"] * n_nh
    code_to_lbl: dict[str, str] = {}

    if not meta42_sub_df.is_empty():
        sorted_m = meta42_sub_df.sort("nhood_idx")
        for r in sorted_m.select(["nhood_idx", "meta_state", "meta_state_label"]).to_dicts():
            idx = int(r["nhood_idx"])
            ms = str(r["meta_state"])
            ml = str(r["meta_state_label"])
            code_to_lbl[ms] = ml
            if 0 <= idx < n_nh:
                meta_codes[idx] = ms
                meta_lbls[idx] = ml

    readable_raw: list[str] = []
    for idx in range(n_nh):
        ct = raw_cts[idx]
        if ct.strip().lower() not in ("unspecified", "unannotated", "unknown", ""):
            readable_raw.append(ct)
        elif meta_lbls[idx] != "Unassigned":
            readable_raw.append(meta_lbls[idx])
        else:
            readable_raw.append(f"Nhood_{idx}")

    fallback_codes = tuple(
        meta_codes[i] if meta_codes[i] != "Unassigned" else readable_raw[i]
        for i in range(n_nh)
    )
    return tuple(readable_raw), fallback_codes, code_to_lbl


def build_view_umap_bundle_pure(
    raw: CohortRawVizArtifacts,
    ctx: CohortDesignContext,
    view: ActivePairView,
    cli_cfg: PairwiseUmapCliConfig,
    viz_cfg: UmapPairVizConfig,
) -> Result[PairwiseUmapBundle, str]:
    """Constructs the complete `PairwiseUmapBundle` for either `SUPER_META42` (T4) or `RAW_MODULES` (T0)."""
    da_fdr = raw.da_results_df["FDR"].cast(pl.Float64).to_numpy()
    da_lfc = raw.da_results_df["logFC"].cast(pl.Float64).to_numpy()
    raw_readable, meta_codes, code_to_lbl = _resolve_raw_nhood_labels(
        raw.da_results_df, raw.meta42_sub_df
    )
    lib_seq = tuple(float(x) for x in ctx.library_sizes.tolist())

    match view:
        case ActivePairView.RAW_MODULES:
            n_sig_q = (
                int(raw.t0_quadrant_df.filter(pl.col("is_significant")).height)
                if (not raw.t0_quadrant_df.is_empty() and "is_significant" in raw.t0_quadrant_df.columns)
                else 0
            )
            n_sig_c = (
                int(raw.t0_continuous_df.filter(pl.col("is_significant")).height)
                if (not raw.t0_continuous_df.is_empty() and "is_significant" in raw.t0_continuous_df.columns)
                else 0
            )
            pw_res = PairwiseAnalysisResult(
                pairs_df=raw.t0_quadrant_df,
                continuous_pairs_df=raw.t0_continuous_df,
                modules_df=raw.t0_modules_df,
                n_total_nhoods=raw.da_results_df.height,
                n_active_nhoods=raw.da_results_df.height,
                n_evaluated_pairs=max(raw.t0_quadrant_df.height, raw.t0_continuous_df.height),
                n_significant_quadrant_pairs=n_sig_q,
                n_significant_continuous_pairs=n_sig_c,
                n_modules=raw.t0_modules_df.height,
            )
            return prepare_pairwise_umap_bundle(
                nhoods_mat=raw.nhoods_mat,
                umap_coords=ctx.umap_coords,
                count_matrix=ctx.count_matrix,
                response_labels=ctx.response_labels,
                pairwise_res=pw_res,
                sample_library_sizes=Some(lib_seq),
                state_labels=Some(raw_readable),
                cell_labels=Some(ctx.cell_labels),
                univariate_fdr=Some(tuple(float(x) for x in da_fdr.tolist())),
                univariate_logfc=Some(tuple(float(x) for x in da_lfc.tolist())),
                anchor_state=cli_cfg.anchor_state,
                anchor_nhood=cli_cfg.anchor_nhood,
                config=viz_cfg,
            )

        case ActivePairView.SUPER_META42:
            cg_res = coarse_grain_nhoods_by_graph(
                nhoods_mat=raw.nhoods_mat,
                count_matrix=ctx.count_matrix,
                cell_sample_codes=Some(ctx.cell_eval_codes),
                sample_library_sizes=Some(ctx.library_sizes),
                nhood_group_labels=Some(meta_codes),
                cell_types=Some(ctx.cell_labels),
                min_jaccard=0.10,
                leiden_resolution=1.0,
            )
            match cg_res:
                case Failure(err):
                    return Failure(f"[{raw.accession}] Super-nhood coarse-graining failed: {err}")
                case Success(cg_data):
                    pass

            # Pool univariate FDR (conservative min across constituent kNN neighborhoods: strict de novo criterion!)
            super_fdr: list[float] = []
            super_lfc: list[float] = []
            super_readable_list: list[str] = []
            for g_idx, mem_indices in enumerate(cg_data.member_nhood_indices):
                st_code = cg_data.super_nhood_labels[g_idx]
                super_readable_list.append(code_to_lbl.get(st_code, st_code))
                if mem_indices:
                    mem_arr = np.array(mem_indices, dtype=np.int64)
                    super_fdr.append(float(np.min(da_fdr[mem_arr])))
                    super_lfc.append(float(np.mean(da_lfc[mem_arr])))
                else:
                    super_fdr.append(1.0)
                    super_lfc.append(0.0)
            super_readable = tuple(super_readable_list)

            min_arm = min(ctx.n_responders, ctx.n_non_responders)
            eff_min_donors = min(3.0, float(max(2, min_arm)))
            eff_min_total = min(6.0, float(max(4, ctx.n_responders + ctx.n_non_responders - 1)))

            # Calibrate DMPP Knockoff+ step-up across CROSS_CELL_TYPE_ONLY when CROSS_COMPARTMENT_ONLY
            # is requested so Barber-Candes' (1 + #{null>=t}) / R(t) finite-sample quantization barrier
            # (R >= 1/alpha = 10) never censors true cross-compartment families with <10 bundles
            # (e.g., Gondal's 3 FWER<0.001 chords), then enforce strict CROSS_COMPARTMENT_ONLY post-screen.
            req_cross_comp = cli_cfg.scope_filter == PairScopeFilter.CROSS_COMPARTMENT_ONLY
            screen_scope = (
                PairScopeFilter.CROSS_CELL_TYPE_ONLY
                if req_cross_comp
                else cli_cfg.scope_filter
            )

            t4_post_cg_cfg = PairwiseConfig(
                method=PairwiseMethod.HYBRID_ALL,
                mixing_qc=MixingQCConfig(
                    enabled=cli_cfg.enable_mixing_qc,
                    min_distinct_donors=cli_cfg.min_distinct_donors,
                    min_effective_donors_lisi=cli_cfg.min_effective_lisi,
                    max_single_donor_dominance=cli_cfg.max_donor_dominance,
                    adaptive_small_cohort_scaling=True,
                ),
                overlap=OverlapConfig(
                    max_jaccard_overlap=0.0,
                    scope_filter=screen_scope,
                ),
                presence=PresenceConfig(mode=PresenceMode.HARD_DEPTH_GATED),
                prefilter=PrefilterConfig(
                    mode=PrefilterMode.MARGINAL_ENTROPY_VARIANCE,
                    min_marginal_prevalence=0.08,
                    max_marginal_prevalence=0.92,
                    min_variance_quantile=0.15,
                    apply_tarone_iut_pruning=True,
                    tarone_p_ceiling=cli_cfg.tarone_p_ceiling,
                ),
                screen=ScreenConfig(
                    chunk_size=512,
                    fdr_threshold=cli_cfg.fdr_threshold,
                    min_donors_per_quadrant=eff_min_donors,
                    min_effective_samples=eff_min_total,
                    max_retained_candidates=50000,
                    apply_tarone_iut_pruning=True,
                    tarone_p_ceiling=cli_cfg.tarone_p_ceiling,
                ),
                refinement=RefinementConfig(enabled=True, max_refine_pairs=250),
                modules=PairModuleConfig(enabled=not cli_cfg.consolidate_state_pairs),
            )
            pw_t4_res = analyze_nhood_pairs(
                nhoods_mat=cg_data.nhoods_super,
                count_matrix=cg_data.counts_super,
                response_labels=ctx.response_labels,
                sample_library_sizes=Some(lib_seq),
                strata_labels=ctx.inferred_strata if cli_cfg.use_substudy_strata else Nothing,
                cell_types=Some(ctx.cell_labels),
                nhood_group_labels=Some(super_readable),
                sample_ids=Some(ctx.eval_ids),
                cell_sample_codes=Some(ctx.cell_eval_codes),
                config=t4_post_cg_cfg,
            )
            match pw_t4_res:
                case Failure(err):
                    return Failure(f"[{raw.accession}] T4 pairwise analysis failed: {err}")
                case Success(pw_t4):
                    pass

            match cli_cfg.consolidate_state_pairs:
                case True:
                    dq_res = deduplicate_pairs_to_lead_bundles(
                        pw_t4.pairs_df,
                        require_cross_compartment=req_cross_comp,
                        one_per_state_pair=True,
                    )
                    dc_res = deduplicate_continuous_pairs_to_lead_bundles(
                        pw_t4.continuous_pairs_df,
                        require_cross_compartment=req_cross_comp,
                        one_per_state_pair=True,
                    )
                    match (dq_res, dc_res):
                        case (Failure(err_q), _):
                            return Failure(f"[{raw.accession}] Quadrant bundle consolidation failed: {err_q}")
                        case (_, Failure(err_c)):
                            return Failure(f"[{raw.accession}] Continuous bundle consolidation failed: {err_c}")
                        case (Success(dedup_q), Success(dedup_c)):
                            sig_q_cnt = (
                                int(dedup_q.filter(pl.col("is_significant")).height)
                                if (not dedup_q.is_empty() and "is_significant" in dedup_q.columns)
                                else 0
                            )
                            sig_c_cnt = (
                                int(dedup_c.filter(pl.col("is_significant")).height)
                                if (not dedup_c.is_empty() and "is_significant" in dedup_c.columns)
                                else 0
                            )
                            pw_for_viz = PairwiseAnalysisResult(
                                pairs_df=dedup_q,
                                continuous_pairs_df=dedup_c,
                                modules_df=pl.DataFrame(),
                                n_total_nhoods=pw_t4.n_total_nhoods,
                                n_active_nhoods=pw_t4.n_active_nhoods,
                                n_evaluated_pairs=pw_t4.n_evaluated_pairs,
                                n_significant_quadrant_pairs=sig_q_cnt,
                                n_significant_continuous_pairs=sig_c_cnt,
                                n_modules=sig_q_cnt + sig_c_cnt,
                            )
                case False:
                    filt_q = (
                        pw_t4.pairs_df.filter(pl.col("pair_scope") == PairScope.CROSS_COMPARTMENT.value)
                        if (req_cross_comp and not pw_t4.pairs_df.is_empty() and "pair_scope" in pw_t4.pairs_df.columns)
                        else pw_t4.pairs_df
                    )
                    filt_c = (
                        pw_t4.continuous_pairs_df.filter(pl.col("pair_scope") == PairScope.CROSS_COMPARTMENT.value)
                        if (
                            req_cross_comp
                            and not pw_t4.continuous_pairs_df.is_empty()
                            and "pair_scope" in pw_t4.continuous_pairs_df.columns
                        )
                        else pw_t4.continuous_pairs_df
                    )
                    sig_q_raw = (
                        int(filt_q.filter(pl.col("is_significant")).height)
                        if (not filt_q.is_empty() and "is_significant" in filt_q.columns)
                        else 0
                    )
                    sig_c_raw = (
                        int(filt_c.filter(pl.col("is_significant")).height)
                        if (not filt_c.is_empty() and "is_significant" in filt_c.columns)
                        else 0
                    )
                    pw_for_viz = (
                        PairwiseAnalysisResult(
                            pairs_df=filt_q,
                            continuous_pairs_df=filt_c,
                            modules_df=pl.DataFrame(),
                            n_total_nhoods=pw_t4.n_total_nhoods,
                            n_active_nhoods=pw_t4.n_active_nhoods,
                            n_evaluated_pairs=pw_t4.n_evaluated_pairs,
                            n_significant_quadrant_pairs=sig_q_raw,
                            n_significant_continuous_pairs=sig_c_raw,
                            n_modules=sig_q_raw + sig_c_raw,
                        )
                        if req_cross_comp
                        else pw_t4
                    )
                case _ as unreachable_cons:
                    assert_never(unreachable_cons)

            return prepare_pairwise_umap_bundle(
                nhoods_mat=cg_data.nhoods_super,
                umap_coords=ctx.umap_coords,
                count_matrix=cg_data.counts_super,
                response_labels=ctx.response_labels,
                pairwise_res=pw_for_viz,
                sample_library_sizes=Some(lib_seq),
                state_labels=Some(super_readable),
                cell_labels=Some(ctx.cell_labels),
                univariate_fdr=Some(tuple(super_fdr)),
                univariate_logfc=Some(tuple(super_lfc)),
                anchor_state=cli_cfg.anchor_state,
                anchor_nhood=cli_cfg.anchor_nhood,
                config=viz_cfg,
            )

        case _ as unreachable_v:
            assert_never(unreachable_v)


def upgrade_svg_with_inkscape_layers(
    svg_content: str,
    figure_title: str,
) -> Result[str, str]:
    """Wraps a hybrid Vega-Lite SVG into native Inkscape layers (`01_Background`, `02_Pairwise_UMAP_Panels`)."""
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
    namedview.set("id", "namedview-pairwise-umap")
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
    layer_panels.set(f"{{{INKSCAPE_NS}}}label", "02_Pairwise_UMAP_Panels")
    for c_node in content_nodes:
        layer_panels.append(c_node)

    serialized = ET.tostring(root, encoding="unicode")
    # Ensure explicit xmlns:inkscape and xmlns:sodipodi declarations on <svg> even when prefixes are only in attributes
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


def render_markdown_summary_report(
    summaries: Sequence[RenderedCohortViewSummary],
    output_dir: Path,
) -> str:
    """Purely formats the comprehensive Markdown report summarizing all pairwise UMAP figures."""
    lines: list[str] = [
        "# Pairwise Milopy Single-Cell UMAP Visualization Summary",
        "",
        "## 1. Architectural Overview",
        "",
        "- **Manifold-Faithful Geometric Medoids**: Every neighborhood or Consensus-Meta42 super-neighborhood "
        "`u` is anchored strictly at an observed member cell `i*(u) = argmin_{i in C_u} ||x_i - median(C_u)||_2^2`, "
        "eliminating floating arithmetic centroids across UMAP tears.",
        "- **Outward-Bowing Quadratic Bezier Chords**: Significant pairs (`FDR_IUT < 0.10`) are deduplicated "
        "into non-redundant ecosystem modules and routed along swap-invariant outward-bowing Bezier arcs "
        "with polar endpoint glyphs (**Filled Okabe Blue**: `+` Responder-enriched; **Hollow White-Core "
        "Okabe Vermilion Ring**: `-` excluded/gated partner; **Okabe Orange Diamond**: `--` dual-absence desert).",
        "- **Emergent De Novo Hub Isolation**: Neighborhoods that are non-significant in 1D univariate Milo "
        "(`FDR_univ >= 0.10`) yet participate in significant 2D pairwise interactions (`FDR_IUT < 0.10`) are "
        "highlighted with concentric `#0F172A` target rings.",
        "- **Cross-Compartment & Donor-Mixing QC Defaults (`SUPER_META42`)**: By default, `SUPER_META42` enforces "
        "`PairScopeFilter.CROSS_COMPARTMENT_ONLY` (eliminating intra-island T<->T and myeloid<->myeloid lineage "
        "self-pairs), pre-gates private single-patient neighborhoods via `MixingQCConfig(min_distinct_donors=3, "
        "min_effective_donors_lisi=2.5, max_single_donor_dominance=0.50)`, stratifies multi-study cohorts "
        "(e.g., `CELLxGENE_7b20c613` `<PMID>ggg<Patient>`) via Mantel-Haenszel (`CMH`) + within-stratum residualized "
        "Welch testing, and consolidates significant interactions to strictly 1 canonical lead chord per "
        "`(State_A, State_B)` pair across both discrete quadrant and continuous axes.",
        "- **Hybrid Vector-Raster Inkscape SVGs**: Dense 12,000-cell UMAP manifold clouds (`>= 250` marks) are "
        "embedded at 300 DPI inside native Inkscape layers while 100% of Bezier chords, endpoint glyphs, "
        "concentric rings, axes, titles, and legends remain crisp editable vector SVG primitives.",
        "",
        "## 2. Cohort-by-Cohort Pairwise UMAP Metrics",
        "",
        "| Accession | Clinical Cohort | Representation | Evals | Units | Sig Quad (`q<0.10`) | Sig Cont (`q<0.10`) | De Novo Hubs | Displayed Chords | Focal Anchor (Panel c) | #1 Lead Pair (Panel d) |",
        "| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- | :--- |",
    ]
    for s in summaries:
        lead_txt = (
            f"{s.lead_pair_desc} (`q={s.lead_fdr_iut:.3g}`)"
            if s.lead_fdr_iut < 1.0
            else "No sig pair (`q >= 0.10`)"
        )
        lines.append(
            f"| `{s.accession}` | {s.cohort_label} | `{s.view.value}` | {s.n_evals} | "
            f"{s.n_units:,} | **{s.n_sig_quadrant:,}** | **{s.n_sig_continuous:,}** | "
            f"{s.n_denovo_units:,} | {s.n_displayed_chords} | {s.anchor_label} | {lead_txt} |"
        )

    lines.extend(
        [
            "",
            "## 3. Exported Vector & 300 DPI Raster Figures",
            "",
            f"- Cross-Cohort 9-Panel Landscape Grid: `{output_dir / 'figures' / 'nine_cohorts_pairwise_umap_grid.svg'}` (`.png`)",
        ]
    )
    for s in summaries:
        if s.exported_composite:
            lines.append(
                f"- `{s.accession}` (`{s.view.value}` 2x2 Composite + Standalone Panels + Parquets): "
                f"`{output_dir / s.accession / f'{s.accession}_{s.view.value}_umap_composite.svg'}`"
            )
    lines.append("")
    return "\n".join(lines)


# =============================================================================
# 3. Imperative Shell — Disk Loaders, Hybrid Inkscape SVG Exporter & CLI Main
# =============================================================================


def load_cohort_viz_artifacts_from_disk(
    accession: str,
    config: PairwiseUmapCliConfig,
    all_meta42_df: pl.DataFrame,
) -> Result[CohortRawVizArtifacts, str]:
    """Side-effecting loader reading precomputed Milo UMAP, DA, T0 pairwise, and Meta42 tables."""
    c_milo = config.milo_response_dir / accession
    c_pw = config.pairwise_t0_dir / accession

    mat_path = c_milo / "nhood_matrix.npz"
    cells_path = c_milo / "cell_level_scores.parquet"
    da_path = c_milo / "da_results.parquet"
    q_path = c_pw / "quadrant_pairs.parquet"
    c_path = c_pw / "continuous_pairs.parquet"
    m_path = c_pw / "pair_modules.parquet"

    for req in (mat_path, cells_path, da_path):
        if not req.exists():
            return Failure(f"[{accession}] Missing Milo artifact: {req}")

    try:
        nhoods_mat = sp.csc_matrix(sp.load_npz(mat_path))
        cells_df = pl.read_parquet(cells_path)
        da_df = pl.read_parquet(da_path)
        q_df = pl.read_parquet(q_path) if q_path.exists() else pl.DataFrame()
        c_df = pl.read_parquet(c_path) if c_path.exists() else pl.DataFrame()
        m_df = pl.read_parquet(m_path) if m_path.exists() else pl.DataFrame()
        meta_sub = (
            all_meta42_df.filter(pl.col("accession") == accession)
            if not all_meta42_df.is_empty()
            else pl.DataFrame()
        )
    except Exception as exc:
        return Failure(f"[{accession}] Disk read failure: {exc}")

    return Success(
        CohortRawVizArtifacts(
            accession=accession,
            nhoods_mat=nhoods_mat,
            cells_df=cells_df,
            da_results_df=da_df,
            t0_quadrant_df=q_df,
            t0_continuous_df=c_df,
            t0_modules_df=m_df,
            meta42_sub_df=meta_sub,
        )
    )


def export_hybrid_inkscape_figure(
    chart: AltairComposite,
    base_path: Path,
    figure_title: str,
) -> Result[tuple[Path, Path], str]:
    """Compiles Altair chart -> hybrid SVG (`min_group_points=250`) -> Inkscape layers -> 300 DPI PNG."""
    svg_path = base_path.parent / f"{base_path.stem}.svg"
    png_path = base_path.parent / f"{base_path.stem}.png"
    try:
        svg_path.parent.mkdir(parents=True, exist_ok=True)
        with alt.data_transformers.disable_max_rows():
            spec = chart.to_dict()
        raw_svg = str(vlc.vegalite_to_svg(spec))
    except Exception as exc:
        return Failure(f"Altair-to-SVG compilation failed for {base_path}: {exc}")

    # Rasterize ONLY dense UMAP cell clouds (>= 250 marks); keep all chords/nodes/rings/labels vector!
    rast_res = rasterize_svg_datapoints(
        raw_svg,
        max_points=1000,
        min_group_points=250,
        scale=3.0,
    )
    match rast_res:
        case Failure(err):
            return Failure(f"Hybrid SVG datapoint rasterization failed for {base_path}: {err}")
        case Success(hybrid_svg):
            pass

    ink_res = upgrade_svg_with_inkscape_layers(hybrid_svg, figure_title=figure_title)
    match ink_res:
        case Failure(err):
            return Failure(err)
        case Success(inkscape_svg):
            pass

    try:
        svg_path.write_text(inkscape_svg, encoding="utf-8")
    except Exception as exc:
        return Failure(f"Failed writing SVG {svg_path}: {exc}")

    png_res = rasterize_svg_to_png(
        svg_path=svg_path,
        png_path=png_path,
        scale=2.5,
        background="white",
    )
    match png_res:
        case Failure(err):
            return Failure(err)
        case Success(written_png):
            return Success((svg_path, written_png))
        case _ as unreachable_p:
            return Failure(f"Unexpected Result state: {unreachable_p}")


def run_pairwise_umap_pipeline(
    cli_cfg: PairwiseUmapCliConfig,
    loader: CohortUmapArtifactLoader = load_cohort_viz_artifacts_from_disk,
    exporter: HybridVectorFigureExporter = export_hybrid_inkscape_figure,
) -> Result[tuple[RenderedCohortViewSummary, ...], str]:
    """Orchestrates 9-cohort UMAP bundle generation, composite/standalone SVG exports, and reporting."""
    t_start = time.perf_counter()
    all_meta42_df = (
        pl.read_parquet(cli_cfg.meta42_annotations_path)
        if cli_cfg.meta42_annotations_path.exists()
        else pl.DataFrame()
    )

    grid_cohorts, composite_cohorts = resolve_cohorts_to_process(cli_cfg.cohort_filter)
    active_views = expand_pair_views(cli_cfg.pair_mode)
    primary_grid_view = active_views[0]

    viz_cfg = UmapPairVizConfig(
        max_chords=cli_cfg.max_chords,
        max_bg_cells=cli_cfg.max_bg_cells,
        fdr_threshold=cli_cfg.fdr_threshold,
    )
    grid_viz_cfg = UmapPairVizConfig(
        max_chords=min(cli_cfg.max_chords, 24),
        max_bg_cells=min(cli_cfg.max_bg_cells, 8000),
        fdr_threshold=cli_cfg.fdr_threshold,
    )

    summaries: list[RenderedCohortViewSummary] = []
    grid_cells: list[tuple[PairwiseUmapBundle, str, str]] = []

    for acc in grid_cohorts:
        load_res = loader(acc, cli_cfg, all_meta42_df)
        match load_res:
            case Failure(err):
                return Failure(err)
            case Success(raw_art):
                pass

        ctx_res = extract_cohort_design_context(raw_art)
        match ctx_res:
            case Failure(err):
                return Failure(err)
            case Success(ctx):
                pass

        n_evals = ctx.n_responders + ctx.n_non_responders
        do_composite = acc in composite_cohorts
        views_for_acc = active_views if do_composite else (primary_grid_view,)

        for view in views_for_acc:
            cur_viz_cfg = viz_cfg if do_composite else grid_viz_cfg
            b_res = build_view_umap_bundle_pure(
                raw=raw_art,
                ctx=ctx,
                view=view,
                cli_cfg=cli_cfg,
                viz_cfg=cur_viz_cfg,
            )
            match b_res:
                case Failure(err):
                    return Failure(err)
                case Success(bundle):
                    pass

            n_arcs = (
                int(bundle.bundled_arcs_df["arc_id"].n_unique())
                if not bundle.bundled_arcs_df.is_empty()
                else 0
            )
            if view == primary_grid_view:
                cmh_suffix = (
                    " [CMH]"
                    if (
                        view == ActivePairView.SUPER_META42
                        and cli_cfg.use_substudy_strata
                        and isinstance(ctx.inferred_strata, Some)
                    )
                    else ""
                )
                grid_sub = (
                    f"N={n_evals} ({ctx.n_responders}R/{ctx.n_non_responders}NR){cmh_suffix} | "
                    f"{bundle.n_sig_quadrant:,} quad + {bundle.n_sig_continuous:,} cont sig (q<0.10)"
                )
                grid_cells.append((bundle, f"{acc} — {ctx.cohort_label}", grid_sub))

            if do_composite:
                cohort_out = cli_cfg.output_dir / acc
                cohort_out.mkdir(parents=True, exist_ok=True)
                bundle.nhood_nodes_df.write_parquet(
                    cohort_out / f"{acc}_{view.value}_nhood_umap_nodes.parquet"
                )
                bundle.bundled_arcs_df.write_parquet(
                    cohort_out / f"{acc}_{view.value}_bundled_arcs.parquet"
                )
                bundle.spotlight_quadrant_df.write_parquet(
                    cohort_out / f"{acc}_{view.value}_spotlight_quadrants.parquet"
                )

                comp_title = f"{acc} ({ctx.cohort_label}) — Pairwise Neighborhood UMAP Architecture"
                comp_sub = (
                    f"{pair_view_display_title(view, cli_cfg)} | N={n_evals} biopsies "
                    f"({ctx.n_responders}R / {ctx.n_non_responders}NR) | "
                    f"{bundle.n_sig_quadrant:,} quadrant & {bundle.n_sig_continuous:,} continuous hits (FDR < {cli_cfg.fdr_threshold:.2f})"
                )
                comp_chart = build_pairwise_umap_composite(
                    bundle=bundle,
                    title=comp_title,
                    subtitle=comp_sub,
                    config=viz_cfg,
                )
                exp_comp = exporter(
                    comp_chart,
                    cohort_out / f"{acc}_{view.value}_umap_composite",
                    comp_title,
                )
                match exp_comp:
                    case Failure(err):
                        return Failure(err)
                    case Success(_):
                        pass

                # Also export standalone high-resolution single panels for manuscript assembly
                standalone_panels: tuple[tuple[str, AltairComposite], ...] = (
                    (
                        "panel_a_bundled_chords",
                        apply_nature_methods_theme(
                            build_chord_umap_panel(bundle, config=viz_cfg, panel_prefix="a")
                        ),
                    ),
                    (
                        "panel_b_denovo_hubs",
                        apply_nature_methods_theme(
                            build_denovo_hub_umap_panel(bundle, config=viz_cfg, panel_prefix="b")
                        ),
                    ),
                    (
                        "panel_c_anchor_field",
                        apply_nature_methods_theme(
                            build_anchor_partner_umap_panel(bundle, config=viz_cfg, panel_prefix="c")
                        ),
                    ),
                    (
                        "panel_d_lead_spotlight",
                        apply_nature_methods_theme(
                            build_lead_pair_spotlight_panel(bundle, config=viz_cfg, panel_prefix="d")
                        ),
                    ),
                )
                for p_slug, p_chart in standalone_panels:
                    exp_p = exporter(
                        p_chart,
                        cohort_out / f"{acc}_{view.value}_{p_slug}",
                        f"{acc} {view.value} {p_slug}",
                    )
                    match exp_p:
                        case Failure(err):
                            return Failure(err)
                        case Success(_):
                            pass

            lead_meta = bundle.spotlight_meta
            lead_desc = (
                f"{lead_meta.label_u} x {lead_meta.label_v} [{lead_meta.regime_label}]"
                if lead_meta.has_pair
                else "No sig pair"
            )
            summaries.append(
                RenderedCohortViewSummary(
                    accession=acc,
                    cohort_label=ctx.cohort_label,
                    view=view,
                    n_evals=n_evals,
                    n_cells=raw_art.cells_df.height,
                    n_units=bundle.nhood_nodes_df.height,
                    n_sig_quadrant=bundle.n_sig_quadrant,
                    n_sig_continuous=bundle.n_sig_continuous,
                    n_denovo_units=bundle.n_denovo_nhoods,
                    n_displayed_chords=n_arcs,
                    anchor_label=bundle.anchor_label,
                    lead_pair_desc=lead_desc,
                    lead_fdr_iut=lead_meta.fdr_iut,
                    exported_composite=do_composite,
                )
            )
            print(
                f"[{acc}][{view.value}] units={bundle.nhood_nodes_df.height:,}, "
                f"sig_quad={bundle.n_sig_quadrant:,}, sig_cont={bundle.n_sig_continuous:,}, "
                f"denovo={bundle.n_denovo_nhoods:,}, chords={n_arcs}"
            )

    if grid_cells:
        grid_title = "Pan-Cohort Pairwise Neighborhood UMAP Atlas Across Nine Clinical ICB Datasets"
        grid_subtitle = (
            f"{pair_view_display_title(primary_grid_view, cli_cfg)} | Outward-bowing Bezier chords connect "
            f"cell-disjoint pairs at IUT FDR < {cli_cfg.fdr_threshold:.2f} "
            "(Filled Blue: + | Hollow Vermilion: - Gate | Grey Rings: De Novo Hubs)"
        )
        grid_res = build_cross_cohort_umap_grid(
            cells=grid_cells,
            title=grid_title,
            subtitle=grid_subtitle,
            n_cols=3,
        )
        match grid_res:
            case Failure(err):
                return Failure(err)
            case Success(grid_chart):
                exp_g = exporter(
                    grid_chart,
                    cli_cfg.output_dir / "figures" / "nine_cohorts_pairwise_umap_grid",
                    grid_title,
                )
                match exp_g:
                    case Failure(err):
                        return Failure(err)
                    case Success(_):
                        pass

    reports_dir = cli_cfg.output_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    report_md = render_markdown_summary_report(summaries, cli_cfg.output_dir)
    (reports_dir / "pairwise_umap_visualization_summary.md").write_text(report_md, encoding="utf-8")

    summary_rows = [s.model_dump(mode="json") for s in summaries]
    pl.DataFrame(summary_rows).write_parquet(
        cli_cfg.output_dir / "pairwise_umap_summary_metrics.parquet"
    )
    elapsed = time.perf_counter() - t_start
    print(f"Completed pairwise UMAP visualization pipeline in {elapsed:.1f}s -> {cli_cfg.output_dir}")
    return Success(tuple(summaries))


def parse_cli_args(argv: Sequence[str]) -> PairwiseUmapCliConfig:
    """Parses command-line flags into an immutable `PairwiseUmapCliConfig`."""
    parser = argparse.ArgumentParser(
        description="Render publication-grade pairwise Milopy UMAP composites and 9-cohort grids."
    )
    parser.add_argument(
        "--pair-mode",
        type=str,
        choices=[m.value for m in PairRepresentationMode],
        default=PairRepresentationMode.BOTH.value,
        help="Pairwise representation to visualize (default: both).",
    )
    parser.add_argument(
        "--cohorts",
        type=str,
        default="",
        help="Optional comma-separated cohort accessions (defaults to all 9 cohorts).",
    )
    parser.add_argument(
        "--max-chords",
        type=int,
        default=32,
        help="Maximum non-redundant ecosystem Bezier chords per UMAP panel (default: 32).",
    )
    parser.add_argument(
        "--max-bg-cells",
        type=int,
        default=12000,
        help="Maximum background manifold cells rendered at 300 DPI per panel (default: 12000).",
    )
    parser.add_argument(
        "--anchor-state",
        type=str,
        default="",
        help="Optional focal cell-state / Meta42 query for Panel c (auto-selects top hub if empty).",
    )
    parser.add_argument(
        "--anchor-nhood",
        type=int,
        default=-1,
        help="Optional focal neighborhood index >= 0 for Panel c.",
    )
    parser.add_argument(
        "--scope-filter",
        type=str,
        choices=[s.value for s in PairScopeFilter],
        default=PairScopeFilter.CROSS_COMPARTMENT_ONLY.value,
        help="Pairwise lineage scope filter for SUPER_META42 (default: cross_compartment_only).",
    )
    parser.add_argument(
        "--disable-mixing-qc",
        action="store_true",
        help="Disable per-neighborhood donor-mixing QC pre-gate (LISI_pt & max donor dominance).",
    )
    parser.add_argument(
        "--min-effective-lisi",
        type=float,
        default=2.5,
        help="Minimum effective donor LISI_pt threshold when mixing QC is enabled (default: 2.5).",
    )
    parser.add_argument(
        "--max-donor-dominance",
        type=float,
        default=0.50,
        help="Maximum single-donor cell share p_max when mixing QC is enabled (default: 0.50).",
    )
    parser.add_argument(
        "--disable-substudy-strata",
        action="store_true",
        help="Disable automatic <PMID>ggg<Patient> sub-study CMH/Welch stratification.",
    )
    parser.add_argument(
        "--no-consolidate-state-pairs",
        action="store_true",
        help="Disable 1-chord-per-(State_A, State_B) canonical bundle deduplication.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/results/milopy_pairwise_umap"),
        help="Output directory for hybrid Inkscape SVGs, 300 DPI PNGs, and Parquets.",
    )
    ns = parser.parse_args(list(argv))

    cohorts_clean = tuple(c.strip() for c in str(ns.cohorts).split(",") if c.strip())
    cohort_maybe: Maybe[tuple[str, ...]] = Some(cohorts_clean) if cohorts_clean else Nothing
    anchor_st_maybe: Maybe[str] = (
        Some(str(ns.anchor_state).strip()) if str(ns.anchor_state).strip() else Nothing
    )
    anchor_nh_maybe: Maybe[int] = Some(int(ns.anchor_nhood)) if int(ns.anchor_nhood) >= 0 else Nothing

    return PairwiseUmapCliConfig(
        output_dir=Path(ns.output_dir),
        pair_mode=PairRepresentationMode(str(ns.pair_mode)),
        scope_filter=PairScopeFilter(str(ns.scope_filter)),
        enable_mixing_qc=not bool(ns.disable_mixing_qc),
        min_effective_lisi=float(ns.min_effective_lisi),
        max_donor_dominance=float(ns.max_donor_dominance),
        use_substudy_strata=not bool(ns.disable_substudy_strata),
        consolidate_state_pairs=not bool(ns.no_consolidate_state_pairs),
        cohort_filter=cohort_maybe,
        anchor_state=anchor_st_maybe,
        anchor_nhood=anchor_nh_maybe,
        max_chords=int(ns.max_chords),
        max_bg_cells=int(ns.max_bg_cells),
    )


def main(argv: Sequence[str] = ()) -> int:
    """CLI entrypoint returning POSIX exit code."""
    eff_argv = argv if argv else tuple(sys.argv[1:])
    cfg = parse_cli_args(eff_argv)
    match run_pairwise_umap_pipeline(cfg):
        case Failure(err):
            sys.stderr.write(f"ERROR: {err}\n")
            return 1
        case Success(_):
            return 0
        case _ as unreachable_m:
            sys.stderr.write(f"Unexpected Result state: {unreachable_m}\n")
            return 1


if __name__ == "__main__":
    sys.exit(main())
