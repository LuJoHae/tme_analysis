#!/usr/bin/env python3
"""Pairwise Milopy Differential Neighborhood Co-Association Across Nine Clinical Cohorts.

Applies `milopy.pairwise` (4-quadrant presence/absence `++`, `+-`, `-+`, `--`, Aitchison
log-ratio balances, orthogonalized co-abundance products, Firth penalized LRT refinement,
and 2D product-graph ecosystem module clustering) to all 9 pre-analyzed single-cell
immunotherapy response cohorts in `output/results/milopy_response_analysis/`.

Architecture strictly follows functional Python standards:
- Functional Core, Imperative Shell (pure data/chart transformations, IO isolated at shell)
- Frozen Pydantic models (`ConfigDict(frozen=True)`), zero `None` / `Optional` (`returns.maybe.Maybe`)
- Polymorphism via `typing.Protocol` typeclasses (dependency-injected loaders and exporters)
- Exhaustive structural pattern matching with `typing.assert_never`
- Railway error handling via `returns.result.Result`
- Polars DataFrames exclusively and Altair vector SVG exports (Nature Methods wireframe)
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from enum import StrEnum
from pathlib import Path
import sys
import time
from typing import Protocol, TypeAlias, assert_never

import altair as alt
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import sparse as sp

from milopy.pairwise import (
    ContinuousAxis,
    EnrichmentDirection,
    OverlapConfig,
    PairModuleConfig,
    PairScope,
    PairwiseAnalysisResult,
    PairwiseConfig,
    PairwiseMethod,
    PresenceConfig,
    PresenceMode,
    QuadrantState,
    RefinementConfig,
    ScreenConfig,
    analyze_nhood_pairs,
)
from milopy.pairwise.types import FloatArray, SparseMatrix
from plotting_utils import (  # type: ignore[import-untyped]
    COHORT_METADATA,
    COHORT_PANELS,
    COLOR_BG_WHITE,
    COLOR_HAIRLINE,
    COLOR_LIGHT_GREY,
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


# =============================================================================
# 1. Exhaustive Sum Types & Immutable Domain Models
# =============================================================================


class PairNoveltyClass(StrEnum):
    """Relationship between pairwise significance and univariate Milo significance of (u, v)."""

    BOTH_UNIVARIATE_NS = "both_univariate_ns"
    ONE_UNIVARIATE_SIG = "one_univariate_sig"
    BOTH_UNIVARIATE_SIG = "both_univariate_sig"


class CohortSortOrder(StrEnum):
    """Ordering strategy for cross-cohort summary tables and charts."""

    PANEL_ORDER = "panel_order"
    DISCOVERY_COUNT = "discovery_count"
    ACCESSION = "accession"


class CohortPairwiseRunConfig(BaseModel):
    """Immutable configuration for running pairwise Milopy across the 9 clinical cohorts."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    input_dir: Path = Field(default=Path("output/results/milopy_response_analysis"))
    output_dir: Path = Field(default=Path("output/results/milopy_pairwise_analysis"))
    cohort_filter: Maybe[tuple[str, ...]] = Field(default=Nothing)
    fdr_threshold: float = Field(default=0.10, gt=0.0, le=1.0)
    max_jaccard_overlap: float = Field(default=0.0, ge=0.0, le=1.0)
    presence_mode: PresenceMode = Field(default=PresenceMode.HARD_DEPTH_GATED)
    pairwise_method: PairwiseMethod = Field(default=PairwiseMethod.HYBRID_ALL)
    min_donors_per_quadrant: float = Field(default=3.0, ge=1.0)
    adaptive_small_cohort_gate: bool = Field(default=True)
    max_refine_pairs: int = Field(default=250, ge=0)
    n_permutations: int = Field(default=0, ge=0)
    chunk_size: int = Field(default=512, ge=16)
    max_retained_candidates: int = Field(default=50000, ge=100)
    sort_order: CohortSortOrder = Field(default=CohortSortOrder.PANEL_ORDER)
    render_plots: bool = Field(default=True)


class CachedCohortInputs(BaseModel):
    """Immutable in-memory bundle of pre-analyzed single-cell Milo artifacts for one cohort."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    accession: str
    nhoods_mat: SparseMatrix
    cells_df: pl.DataFrame
    univariate_da_df: pl.DataFrame


class EvaluationDesign(BaseModel):
    """Immutable lesion/donor evaluation matrix C (V_nhoods x J_evals) and metadata."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    accession: str
    eval_ids: tuple[str, ...]
    response_labels: tuple[str, ...]
    count_matrix: FloatArray
    library_sizes: FloatArray
    cell_types: tuple[str, ...]
    n_cells: int
    n_nhoods: int
    n_responders: int
    n_non_responders: int


class CohortPairwiseSummary(BaseModel):
    """Immutable summary statistics comparing univariate and pairwise Milo for one cohort."""

    model_config = ConfigDict(frozen=True)

    accession: str
    cohort_label: str
    indication: str
    author: str
    year: str
    technology: str
    n_cells: int
    n_responders: int
    n_non_responders: int
    n_total_nhoods: int
    n_active_nhoods: int
    n_univariate_sig_r: int
    n_univariate_sig_nr: int
    n_evaluated_pairs: int
    n_sig_quadrant_total: int
    n_sig_pp: int
    n_sig_pm: int
    n_sig_mp: int
    n_sig_mm: int
    n_sig_quad_responder: int
    n_sig_quad_non_responder: int
    n_sig_quad_cross_celltype: int
    n_sig_quad_within_celltype: int
    n_sig_quad_de_novo: int
    n_sig_quad_one_univ: int
    n_sig_quad_both_univ: int
    n_sig_continuous_total: int
    n_sig_log_ratio: int
    n_sig_co_abundance: int
    n_sig_cont_de_novo: int
    n_modules: int
    runtime_seconds: float


class CohortPairwiseBundle(BaseModel):
    """Immutable complete result bundle for one cohort."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    summary: CohortPairwiseSummary
    quadrant_df: pl.DataFrame
    continuous_df: pl.DataFrame
    modules_df: pl.DataFrame


# =============================================================================
# 2. Typeclass Protocols (IO Boundaries)
# =============================================================================


class CohortDataLoader(Protocol):
    """Typeclass protocol for loading cached Milo matrices and tables for a cohort."""

    def __call__(self, accession: str, input_dir: Path) -> Result[CachedCohortInputs, str]: ...


class FigureExporter(Protocol):
    """Typeclass protocol for exporting Altair vector SVG + 300 DPI verification PNG."""

    def __call__(self, chart: AltairChart, base_path: Path) -> Result[tuple[Path, Path], str]: ...


# =============================================================================
# 3. Pure Functional Core — Evaluation Construction, Pairwise Execution & Joins
# =============================================================================


def classify_pair_novelty(sig_u: bool, sig_v: bool) -> PairNoveltyClass:
    """Exhaustively classifies whether a significant pair is de novo vs univariate-anchored."""
    match sig_u:
        case False:
            return (
                PairNoveltyClass.ONE_UNIVARIATE_SIG
                if sig_v
                else PairNoveltyClass.BOTH_UNIVARIATE_NS
            )
        case True:
            return (
                PairNoveltyClass.BOTH_UNIVARIATE_SIG
                if sig_v
                else PairNoveltyClass.ONE_UNIVARIATE_SIG
            )
        case _ as unreachable:
            assert_never(unreachable)


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


def build_evaluation_design(inputs: CachedCohortInputs) -> Result[EvaluationDesign, str]:
    """Builds the exact (V_nhoods x J_evals) count matrix and sample library sizes purely."""
    required_cell_cols = {"patient_id", "clinical_response"}
    missing_cols = sorted(required_cell_cols - set(inputs.cells_df.columns))
    if missing_cols:
        return Failure(
            f"[{inputs.accession}] cells_df missing required columns: {missing_cols}"
        )

    n_cells = inputs.cells_df.height
    mat_rows, n_nhoods = int(inputs.nhoods_mat.shape[0]), int(inputs.nhoods_mat.shape[1])
    if n_cells != mat_rows:
        return Failure(
            f"[{inputs.accession}] cells_df rows ({n_cells}) != nhoods_mat rows ({mat_rows})"
        )
    if inputs.univariate_da_df.height != n_nhoods:
        return Failure(
            f"[{inputs.accession}] univariate_da_df rows ({inputs.univariate_da_df.height}) "
            f"!= nhoods_mat cols ({n_nhoods})"
        )

    norm_cells = inputs.cells_df.lazy().with_columns(
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
    ).with_columns(
        (pl.col("_pt") + pl.lit("__") + pl.col("_resp")).alias("_eval_id")
    ).collect()

    valid_resps = {"responder", "non-responder"}
    observed_resps = set(norm_cells["_resp"].unique().to_list())
    invalid_resps = sorted(observed_resps - valid_resps)
    if invalid_resps:
        return Failure(
            f"[{inputs.accession}] Invalid clinical_response values: {invalid_resps}"
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

    eval_ids = tuple(str(x) for x in eval_meta["_eval_id"].to_list())
    n_evals = len(eval_ids)
    if n_evals < 4:
        return Failure(
            f"[{inputs.accession}] Pairwise analysis requires >= 4 evaluations, got {n_evals}"
        )

    eval_to_idx = {eid: idx for idx, eid in enumerate(eval_ids)}
    cell_eval_codes = np.array(
        [eval_to_idx[str(eid)] for eid in norm_cells["_eval_id"].to_list()],
        dtype=np.int64,
    )
    sample_indicator = sp.csr_matrix(
        (np.ones(n_cells, dtype=np.float64), (np.arange(n_cells, dtype=np.int64), cell_eval_codes)),
        shape=(n_cells, n_evals),
    )

    nhoods_csc = sp.csc_matrix(inputs.nhoods_mat)
    count_matrix = np.asarray(nhoods_csc.T.dot(sample_indicator).toarray(), dtype=np.float64)
    library_sizes = eval_meta["library_size"].to_numpy().astype(np.float64)
    response_labels = tuple(str(r) for r in eval_meta["clinical_response"].to_list())
    n_resp = sum(1 for r in response_labels if r == "responder")
    n_non_resp = n_evals - n_resp

    cell_types = tuple(str(ct) for ct in norm_cells["_ct"].to_list())

    return Success(
        EvaluationDesign(
            accession=inputs.accession,
            eval_ids=eval_ids,
            response_labels=response_labels,
            count_matrix=count_matrix,
            library_sizes=library_sizes,
            cell_types=cell_types,
            n_cells=n_cells,
            n_nhoods=n_nhoods,
            n_responders=n_resp,
            n_non_responders=n_non_resp,
        )
    )


def resolve_cohort_pairwise_config(
    run_cfg: CohortPairwiseRunConfig,
    n_responders: int,
    n_non_responders: int,
) -> PairwiseConfig:
    """Builds the immutable PairwiseConfig adapted cleanly to cohort sample arm sizes."""
    min_arm = min(n_responders, n_non_responders)
    eff_min_donors = (
        min(run_cfg.min_donors_per_quadrant, float(max(2, min_arm)))
        if run_cfg.adaptive_small_cohort_gate
        else run_cfg.min_donors_per_quadrant
    )
    eff_min_total = min(6.0, float(max(4, n_responders + n_non_responders - 1)))

    return PairwiseConfig(
        method=run_cfg.pairwise_method,
        overlap=OverlapConfig(max_jaccard_overlap=run_cfg.max_jaccard_overlap),
        presence=PresenceConfig(mode=run_cfg.presence_mode),
        screen=ScreenConfig(
            chunk_size=run_cfg.chunk_size,
            fdr_threshold=run_cfg.fdr_threshold,
            min_donors_per_quadrant=eff_min_donors,
            min_effective_samples=eff_min_total,
            n_permutations=run_cfg.n_permutations,
            max_retained_candidates=run_cfg.max_retained_candidates,
        ),
        refinement=RefinementConfig(
            enabled=run_cfg.max_refine_pairs > 0,
            max_refine_pairs=run_cfg.max_refine_pairs,
        ),
        modules=PairModuleConfig(enabled=True),
    )


def _prepare_univariate_lookup(
    univariate_da_df: pl.DataFrame,
    fdr_threshold: float,
) -> pl.DataFrame:
    """Extracts 0-based neighborhood index, univariate logFC, FDR, and significance flag."""
    base = (
        univariate_da_df.with_columns((pl.col("Nhood").cast(pl.Int64) - 1).alias("nhood_idx"))
        if (
            "Nhood" in univariate_da_df.columns
            and univariate_da_df.height > 0
            and int(univariate_da_df["Nhood"].cast(pl.Int64).to_numpy().min()) == 1
        )
        else (
            univariate_da_df.with_columns(pl.col("Nhood").cast(pl.Int64).alias("nhood_idx"))
            if "Nhood" in univariate_da_df.columns
            else univariate_da_df.with_row_index("nhood_idx").with_columns(
                pl.col("nhood_idx").cast(pl.Int64)
            )
        )
    )
    return (
        base.lazy()
        .select(
            pl.col("nhood_idx"),
            pl.col("logFC").cast(pl.Float64).alias("univariate_logfc"),
            pl.col("FDR").cast(pl.Float64).alias("univariate_fdr"),
            (pl.col("FDR").cast(pl.Float64) < fdr_threshold).alias("is_univariate_sig"),
        )
        .collect()
    )


def annotate_pairs_with_univariate_da(
    pairs_df: pl.DataFrame,
    univariate_da_df: pl.DataFrame,
    fdr_threshold: float,
) -> pl.DataFrame:
    """Purely joins pairwise results with univariate Milo DA statistics and novelty class."""
    if pairs_df.is_empty():
        return pairs_df.with_columns(
            pl.Series("univariate_logfc_u", [], dtype=pl.Float64),
            pl.Series("univariate_fdr_u", [], dtype=pl.Float64),
            pl.Series("is_univariate_sig_u", [], dtype=pl.Boolean),
            pl.Series("univariate_logfc_v", [], dtype=pl.Float64),
            pl.Series("univariate_fdr_v", [], dtype=pl.Float64),
            pl.Series("is_univariate_sig_v", [], dtype=pl.Boolean),
            pl.Series("novelty_class", [], dtype=pl.String),
        )

    lookup = _prepare_univariate_lookup(univariate_da_df, fdr_threshold)
    lookup_u = lookup.rename(
        {
            "nhood_idx": "nhood_u",
            "univariate_logfc": "univariate_logfc_u",
            "univariate_fdr": "univariate_fdr_u",
            "is_univariate_sig": "is_univariate_sig_u",
        }
    )
    lookup_v = lookup.rename(
        {
            "nhood_idx": "nhood_v",
            "univariate_logfc": "univariate_logfc_v",
            "univariate_fdr": "univariate_fdr_v",
            "is_univariate_sig": "is_univariate_sig_v",
        }
    )

    joined = (
        pairs_df.lazy()
        .join(lookup_u.lazy(), on="nhood_u", how="left")
        .join(lookup_v.lazy(), on="nhood_v", how="left")
        .with_columns(
            pl.col("is_univariate_sig_u").fill_null(False),
            pl.col("is_univariate_sig_v").fill_null(False),
        )
        .with_columns(
            pl.when((~pl.col("is_univariate_sig_u")) & (~pl.col("is_univariate_sig_v")))
            .then(pl.lit(PairNoveltyClass.BOTH_UNIVARIATE_NS.value))
            .when(pl.col("is_univariate_sig_u") & pl.col("is_univariate_sig_v"))
            .then(pl.lit(PairNoveltyClass.BOTH_UNIVARIATE_SIG.value))
            .otherwise(pl.lit(PairNoveltyClass.ONE_UNIVARIATE_SIG.value))
            .alias("novelty_class")
        )
        .collect()
    )
    return joined


def summarize_cohort_bundle(
    design: EvaluationDesign,
    univariate_da_df: pl.DataFrame,
    pairwise_res: PairwiseAnalysisResult,
    quad_df: pl.DataFrame,
    cont_df: pl.DataFrame,
    fdr_threshold: float,
    runtime_seconds: float,
) -> CohortPairwiseSummary:
    """Computes exhaustive univariate + pairwise summary metrics for a single cohort."""
    meta = COHORT_METADATA.get(
        design.accession,
        {
            "indication": "Unknown",
            "author": design.accession,
            "year": "",
            "tech": "scRNA-seq",
        },
    )
    label = get_single_cell_cohort_display_label(design.accession, "cancer_author_year")

    univ_sig = univariate_da_df.filter(pl.col("FDR").cast(pl.Float64) < fdr_threshold)
    n_univ_r = int(univ_sig.filter(pl.col("logFC").cast(pl.Float64) > 0.0).height)
    n_univ_nr = int(univ_sig.filter(pl.col("logFC").cast(pl.Float64) < 0.0).height)

    sig_q = quad_df.filter(pl.col("is_significant"))
    sig_c = cont_df.filter(pl.col("is_significant"))

    return CohortPairwiseSummary(
        accession=design.accession,
        cohort_label=label,
        indication=meta.get("indication", "Unknown"),
        author=meta.get("author", design.accession),
        year=meta.get("year", ""),
        technology=meta.get("tech", "scRNA-seq"),
        n_cells=design.n_cells,
        n_responders=design.n_responders,
        n_non_responders=design.n_non_responders,
        n_total_nhoods=pairwise_res.n_total_nhoods,
        n_active_nhoods=pairwise_res.n_active_nhoods,
        n_univariate_sig_r=n_univ_r,
        n_univariate_sig_nr=n_univ_nr,
        n_evaluated_pairs=pairwise_res.n_evaluated_pairs,
        n_sig_quadrant_total=int(sig_q.height),
        n_sig_pp=int(sig_q.filter(pl.col("quadrant") == QuadrantState.BOTH_PRESENT.value).height),
        n_sig_pm=int(sig_q.filter(pl.col("quadrant") == QuadrantState.FIRST_ONLY.value).height),
        n_sig_mp=int(sig_q.filter(pl.col("quadrant") == QuadrantState.SECOND_ONLY.value).height),
        n_sig_mm=int(sig_q.filter(pl.col("quadrant") == QuadrantState.BOTH_MISSING.value).height),
        n_sig_quad_responder=int(
            sig_q.filter(pl.col("direction") == EnrichmentDirection.RESPONDER.value).height
        ),
        n_sig_quad_non_responder=int(
            sig_q.filter(pl.col("direction") == EnrichmentDirection.NON_RESPONDER.value).height
        ),
        n_sig_quad_cross_celltype=int(
            sig_q.filter(pl.col("pair_scope") == PairScope.CROSS_CELL_TYPE.value).height
        ),
        n_sig_quad_within_celltype=int(
            sig_q.filter(pl.col("pair_scope") == PairScope.WITHIN_CELL_TYPE.value).height
        ),
        n_sig_quad_de_novo=int(
            sig_q.filter(
                pl.col("novelty_class") == PairNoveltyClass.BOTH_UNIVARIATE_NS.value
            ).height
        ),
        n_sig_quad_one_univ=int(
            sig_q.filter(
                pl.col("novelty_class") == PairNoveltyClass.ONE_UNIVARIATE_SIG.value
            ).height
        ),
        n_sig_quad_both_univ=int(
            sig_q.filter(
                pl.col("novelty_class") == PairNoveltyClass.BOTH_UNIVARIATE_SIG.value
            ).height
        ),
        n_sig_continuous_total=int(sig_c.height),
        n_sig_log_ratio=int(
            sig_c.filter(pl.col("axis") == ContinuousAxis.LOG_RATIO_BALANCE.value).height
        ),
        n_sig_co_abundance=int(
            sig_c.filter(pl.col("axis") == ContinuousAxis.CO_ABUNDANCE_PRODUCT.value).height
        ),
        n_sig_cont_de_novo=int(
            sig_c.filter(
                pl.col("novelty_class") == PairNoveltyClass.BOTH_UNIVARIATE_NS.value
            ).height
        ),
        n_modules=pairwise_res.n_modules,
        runtime_seconds=runtime_seconds,
    )


def analyze_single_cohort(
    inputs: CachedCohortInputs,
    run_cfg: CohortPairwiseRunConfig,
    runtime_seconds: float = 0.0,
) -> Result[CohortPairwiseBundle, str]:
    """Pure functional pipeline executing pairwise Milopy and univariate cross-annotation."""
    design_res = build_evaluation_design(inputs)
    match design_res:
        case Failure(err):
            return Failure(err)
        case Success(design):
            pass

    pw_cfg = resolve_cohort_pairwise_config(
        run_cfg=run_cfg,
        n_responders=design.n_responders,
        n_non_responders=design.n_non_responders,
    )

    pw_res = analyze_nhood_pairs(
        nhoods_mat=inputs.nhoods_mat,
        count_matrix=design.count_matrix,
        response_labels=design.response_labels,
        sample_library_sizes=Some(tuple(design.library_sizes.tolist())),
        cell_types=Some(design.cell_types),
        config=pw_cfg,
    )
    match pw_res:
        case Failure(err):
            return Failure(f"[{inputs.accession}] {err}")
        case Success(pw_out):
            pass

    quad_annotated = annotate_pairs_with_univariate_da(
        pairs_df=pw_out.pairs_df,
        univariate_da_df=inputs.univariate_da_df,
        fdr_threshold=run_cfg.fdr_threshold,
    ).with_columns(pl.lit(inputs.accession).alias("accession"))

    cont_annotated = annotate_pairs_with_univariate_da(
        pairs_df=pw_out.continuous_pairs_df,
        univariate_da_df=inputs.univariate_da_df,
        fdr_threshold=run_cfg.fdr_threshold,
    ).with_columns(pl.lit(inputs.accession).alias("accession"))

    modules_annotated = pw_out.modules_df.with_columns(
        pl.lit(inputs.accession).alias("accession")
    )

    summary = summarize_cohort_bundle(
        design=design,
        univariate_da_df=inputs.univariate_da_df,
        pairwise_res=pw_out,
        quad_df=quad_annotated,
        cont_df=cont_annotated,
        fdr_threshold=run_cfg.fdr_threshold,
        runtime_seconds=runtime_seconds,
    )

    return Success(
        CohortPairwiseBundle(
            summary=summary,
            quadrant_df=quad_annotated,
            continuous_df=cont_annotated,
            modules_df=modules_annotated,
        )
    )


def build_cohorts_summary_dataframe(
    bundles: Sequence[CohortPairwiseBundle],
    sort_order: CohortSortOrder = CohortSortOrder.PANEL_ORDER,
) -> pl.DataFrame:
    """Transforms CohortPairwiseBundles into an ordered Polars summary DataFrame."""
    panel_rank = {str(p["acc"]): idx for idx, p in enumerate(COHORT_PANELS)}
    summaries = [b.summary for b in bundles]

    match sort_order:
        case CohortSortOrder.PANEL_ORDER:
            ordered = sorted(summaries, key=lambda s: panel_rank.get(s.accession, 999))
        case CohortSortOrder.DISCOVERY_COUNT:
            ordered = sorted(
                summaries,
                key=lambda s: (s.n_sig_quadrant_total + s.n_sig_continuous_total, s.n_evaluated_pairs),
                reverse=True,
            )
        case CohortSortOrder.ACCESSION:
            ordered = sorted(summaries, key=lambda s: s.accession)
        case _ as unreachable:
            assert_never(unreachable)

    return pl.DataFrame([s.model_dump() for s in ordered])


def _canonicalize_quadrant_lineage_row(
    ct_u: str,
    ct_v: str,
    quadrant: str,
    direction: str,
) -> tuple[str, str, str, str]:
    """Maps index-ordered (ct_u, ct_v, quadrant) into index-invariant canonical lineage tuple."""
    q_enum = QuadrantState(quadrant)
    match q_enum:
        case QuadrantState.BOTH_PRESENT:
            c1, c2 = (ct_u, ct_v) if ct_u <= ct_v else (ct_v, ct_u)
            return (c1, c2, "++ (Both Present)", direction)
        case QuadrantState.BOTH_MISSING:
            c1, c2 = (ct_u, ct_v) if ct_u <= ct_v else (ct_v, ct_u)
            return (c1, c2, "-- (Both Missing)", direction)
        case QuadrantState.FIRST_ONLY:
            return (ct_u, ct_v, "+- (First Present / Second Missing)", direction)
        case QuadrantState.SECOND_ONLY:
            return (ct_v, ct_u, "+- (First Present / Second Missing)", direction)
        case _ as unreachable:
            assert_never(unreachable)


def aggregate_recurrent_celltype_pairs(bundles: Sequence[CohortPairwiseBundle]) -> pl.DataFrame:
    """Aggregates canonical lineage-pair regimes across all analyzed cohorts."""
    empty_schema = {
        "lineage_primary": pl.String,
        "lineage_secondary": pl.String,
        "regime": pl.String,
        "direction": pl.String,
        "pair_label": pl.String,
        "n_cohorts": pl.Int64,
        "cohorts": pl.String,
        "total_sig_pairs": pl.Int64,
        "total_de_novo_pairs": pl.Int64,
        "mean_effect_size": pl.Float64,
    }
    rows: list[dict[str, object]] = []

    for bundle in bundles:
        acc = bundle.summary.accession
        sig_q = bundle.quadrant_df.filter(pl.col("is_significant"))
        for r in sig_q.select(
            ["cell_type_u", "cell_type_v", "quadrant", "direction", "log_or_quadrant", "novelty_class"]
        ).to_dicts():
            c1, c2, reg, dir_s = _canonicalize_quadrant_lineage_row(
                str(r["cell_type_u"]),
                str(r["cell_type_v"]),
                str(r["quadrant"]),
                str(r["direction"]),
            )
            rows.append(
                {
                    "accession": acc,
                    "lineage_primary": c1,
                    "lineage_secondary": c2,
                    "regime": reg,
                    "direction": dir_s,
                    "effect": abs(float(r["log_or_quadrant"])),
                    "is_de_novo": 1 if str(r["novelty_class"]) == PairNoveltyClass.BOTH_UNIVARIATE_NS.value else 0,
                }
            )

        sig_c = bundle.continuous_df.filter(pl.col("is_significant"))
        for r in sig_c.select(
            ["cell_type_u", "cell_type_v", "axis", "direction", "effect_size", "novelty_class"]
        ).to_dicts():
            ax_str = str(r["axis"])
            eff_val = float(r["effect_size"])
            cu, cv = str(r["cell_type_u"]), str(r["cell_type_v"])
            dir_c = str(r["direction"])
            if ax_str == ContinuousAxis.LOG_RATIO_BALANCE.value:
                # log2(u/v) enriched in NR is algebraically identical to log2(v/u) enriched in R
                c1, c2 = (
                    (cu, cv)
                    if dir_c == EnrichmentDirection.RESPONDER.value
                    else (cv, cu)
                )
                dir_c = EnrichmentDirection.RESPONDER.value
                reg = "Log-Ratio (Primary > Secondary)"
            else:
                c1, c2 = (cu, cv) if cu <= cv else (cv, cu)
                reg = "Co-Abundance Coupling"
            rows.append(
                {
                    "accession": acc,
                    "lineage_primary": c1,
                    "lineage_secondary": c2,
                    "regime": reg,
                    "direction": dir_c,
                    "effect": abs(eff_val),
                    "is_de_novo": 1 if str(r["novelty_class"]) == PairNoveltyClass.BOTH_UNIVARIATE_NS.value else 0,
                }
            )

    if not rows:
        return pl.DataFrame(schema=empty_schema)

    raw_df = pl.DataFrame(rows)
    return (
        raw_df.lazy()
        .group_by(["lineage_primary", "lineage_secondary", "regime", "direction"])
        .agg(
            pl.col("accession").n_unique().cast(pl.Int64).alias("n_cohorts"),
            pl.col("accession").unique().sort().str.join(", ").alias("cohorts"),
            pl.len().cast(pl.Int64).alias("total_sig_pairs"),
            pl.col("is_de_novo").sum().cast(pl.Int64).alias("total_de_novo_pairs"),
            pl.col("effect").mean().cast(pl.Float64).alias("mean_effect_size"),
        )
        .with_columns(
            (
                pl.col("lineage_primary")
                + pl.lit(" × ")
                + pl.col("lineage_secondary")
                + pl.lit(" [")
                + pl.col("regime")
                + pl.lit("]")
            ).alias("pair_label")
        )
        .sort(["n_cohorts", "total_sig_pairs", "mean_effect_size"], descending=[True, True, True])
        .select(list(empty_schema.keys()))
        .collect()
    )


# =============================================================================
# 4. Pure Declarative Nature Methods Altair Visualizations (Zero Pandas / HTML)
# =============================================================================


def build_cohort_quadrant_volcano_chart(
    quadrant_df: pl.DataFrame,
    summary: CohortPairwiseSummary,
    fdr_threshold: float,
    max_points_display: int = 4000,
) -> Result[alt.LayerChart, str]:
    """Builds a Nature Methods 4-quadrant volcano plot (logOR vs -log10(p_IUT)) for one cohort."""
    if quadrant_df.is_empty():
        placeholder = pl.DataFrame(
            {
                "log_or_quadrant": [0.0],
                "neg_log10_p_iut": [0.0],
                "quadrant_group": ["Not Significant"],
            }
        )
        base_empty = (
            alt.Chart(placeholder)
            .mark_circle(size=20, opacity=0.0)
            .encode(
                x=alt.X("log_or_quadrant:Q", title="Quadrant Mantel–Haenszel log(OR) (R vs. NR)"),
                y=alt.Y("neg_log10_p_iut:Q", title="-log10(IUT Conditional Flip p-value)"),
            )
            .properties(width=340, height=260, title=f"{summary.cohort_label} — No Screened Pairs")
        )
        rule_zero = alt.Chart(pl.DataFrame({"x": [0.0]})).mark_rule(
            color=COLOR_NEUTRAL_GREY, strokeDash=[3, 3], strokeWidth=0.75
        ).encode(x="x:Q")
        return Success(apply_nature_methods_theme(base_empty + rule_zero))

    enriched = quadrant_df.lazy().with_columns(
        (-pl.col("pvalue_iut").clip(lower_bound=1e-300).log10()).alias("neg_log10_p_iut"),
        pl.when(~pl.col("is_significant"))
        .then(pl.lit("Not Significant"))
        .when(pl.col("quadrant") == QuadrantState.BOTH_PRESENT.value)
        .then(pl.lit("++ Both Present"))
        .when(pl.col("quadrant").is_in([QuadrantState.FIRST_ONLY.value, QuadrantState.SECOND_ONLY.value]))
        .then(pl.lit("+- / -+ Present-Missing Gate"))
        .otherwise(pl.lit("-- Both Missing"))
        .alias("quadrant_group"),
    ).collect()

    sig_part = enriched.filter(pl.col("is_significant"))
    ns_part = enriched.filter(~pl.col("is_significant"))
    n_sig = sig_part.height

    sampled = (
        pl.concat(
            [
                ns_part.head(max(500, max_points_display - min(n_sig, max_points_display - 500))),
                sig_part.head(max_points_display - 500),
            ]
        )
        if enriched.height > max_points_display
        else enriched
    ).select(["log_or_quadrant", "neg_log10_p_iut", "quadrant_group"])

    domain_groups = [
        "++ Both Present",
        "+- / -+ Present-Missing Gate",
        "-- Both Missing",
        "Not Significant",
    ]
    range_colors = [
        OKABE_BLUE,
        OKABE_BLUISH_GREEN,
        OKABE_VERMILION,
        "#CBD5E1",
    ]

    pts = (
        alt.Chart(sampled)
        .mark_circle(size=22, opacity=0.78)
        .encode(
            x=alt.X(
                "log_or_quadrant:Q",
                title="Quadrant Mantel–Haenszel log(OR) (R vs. NR)",
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B", gridColor=COLOR_LIGHT_GREY),
            ),
            y=alt.Y(
                "neg_log10_p_iut:Q",
                title="-log10(IUT Conditional Flip p-value)",
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B", gridColor=COLOR_LIGHT_GREY),
            ),
            color=alt.Color(
                "quadrant_group:N",
                scale=alt.Scale(domain=domain_groups, range=range_colors),
                title=f"Regime (FDR < {fdr_threshold:.2f})",
            ),
        )
        .properties(
            width=360,
            height=270,
            title={
                "text": f"{summary.cohort_label} [{summary.accession}]",
                "subtitle": (
                    f"{summary.n_responders} R / {summary.n_non_responders} NR | "
                    f"{summary.n_evaluated_pairs:,} disjoint pairs | "
                    f"{summary.n_sig_quadrant_total:,} sig. quadrant hits ({summary.n_modules} modules)"
                ),
                "fontSize": 11,
                "subtitleFontSize": 9,
                "color": COLOR_TEXT_PRIMARY,
                "subtitleColor": COLOR_TEXT_SECONDARY,
            },
        )
    )

    zero_line = (
        alt.Chart(pl.DataFrame({"x": [0.0]}))
        .mark_rule(color="#64748B", strokeWidth=0.75, strokeDash=[3, 3])
        .encode(x="x:Q")
    )

    return Success(apply_nature_methods_theme(zero_line + pts))


def build_cross_cohort_pairwise_comparison_chart(
    summary_df: pl.DataFrame,
) -> Result[alt.HConcatChart, str]:
    """Builds a 3-panel Nature Methods cross-cohort comparison SVG across all 9 cohorts."""
    if summary_df.is_empty():
        return Failure("Cannot build cross-cohort comparison chart from empty summary_df")

    cohort_order = summary_df["cohort_label"].to_list()
    rows_q: list[dict[str, object]] = []
    rows_c: list[dict[str, object]] = []
    rows_n: list[dict[str, object]] = []

    for r in summary_df.to_dicts():
        lbl = str(r["cohort_label"])
        rows_q.extend(
            [
                {
                    "cohort_label": lbl,
                    "regime": "++ Both Present",
                    "count": int(r["n_sig_pp"]),
                    "log10_count": float(np.log10(1.0 + float(r["n_sig_pp"]))),
                },
                {
                    "cohort_label": lbl,
                    "regime": "+- / -+ Gated",
                    "count": int(r["n_sig_pm"]) + int(r["n_sig_mp"]),
                    "log10_count": float(np.log10(1.0 + float(int(r["n_sig_pm"]) + int(r["n_sig_mp"])))),
                },
                {
                    "cohort_label": lbl,
                    "regime": "-- Both Missing",
                    "count": int(r["n_sig_mm"]),
                    "log10_count": float(np.log10(1.0 + float(r["n_sig_mm"]))),
                },
            ]
        )
        rows_c.extend(
            [
                {
                    "cohort_label": lbl,
                    "axis": "Log-Ratio Balance",
                    "count": int(r["n_sig_log_ratio"]),
                    "log10_count": float(np.log10(1.0 + float(r["n_sig_log_ratio"]))),
                },
                {
                    "cohort_label": lbl,
                    "axis": "Co-Abundance Product",
                    "count": int(r["n_sig_co_abundance"]),
                    "log10_count": float(np.log10(1.0 + float(r["n_sig_co_abundance"]))),
                },
            ]
        )
        de_novo_tot = int(r["n_sig_quad_de_novo"]) + int(r["n_sig_cont_de_novo"])
        anchored_tot = (
            int(r["n_sig_quad_one_univ"])
            + int(r["n_sig_quad_both_univ"])
            + (int(r["n_sig_continuous_total"]) - int(r["n_sig_cont_de_novo"]))
        )
        rows_n.extend(
            [
                {
                    "cohort_label": lbl,
                    "novelty": "De Novo (Both Univariate NS)",
                    "count": de_novo_tot,
                    "log10_count": float(np.log10(1.0 + float(de_novo_tot))),
                },
                {
                    "cohort_label": lbl,
                    "novelty": "Univariate-Anchored Synergy",
                    "count": max(anchored_tot, 0),
                    "log10_count": float(np.log10(1.0 + float(max(anchored_tot, 0)))),
                },
            ]
        )

    df_q = pl.DataFrame(rows_q)
    df_c = pl.DataFrame(rows_c)
    df_n = pl.DataFrame(rows_n)

    panel_a = (
        alt.Chart(df_q)
        .mark_bar(height=7, cornerRadius=0)
        .encode(
            y=alt.Y(
                "cohort_label:N",
                sort=cohort_order,
                title="Clinical Single-Cell Cohort",
                axis=alt.Axis(labelLimit=260, tickColor="#1E293B", domainColor="#1E293B"),
            ),
            yOffset=alt.YOffset("regime:N", sort=["++ Both Present", "+- / -+ Gated", "-- Both Missing"]),
            x=alt.X(
                "log10_count:Q",
                title="log10(1 + Sig. Quadrant Pairs)",
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B"),
            ),
            color=alt.Color(
                "regime:N",
                scale=alt.Scale(
                    domain=["++ Both Present", "+- / -+ Gated", "-- Both Missing"],
                    range=[OKABE_BLUE, OKABE_BLUISH_GREEN, OKABE_VERMILION],
                ),
                title="2x2 Quadrant Regime",
                legend=alt.Legend(orient="bottom", direction="vertical"),
            ),
        )
        .properties(width=200, height=290, title="a. Irreducible 2x2 Quadrant Pairs")
    )

    panel_b = (
        alt.Chart(df_c)
        .mark_bar(height=9, cornerRadius=0)
        .encode(
            y=alt.Y(
                "cohort_label:N",
                sort=cohort_order,
                title="",
                axis=alt.Axis(labels=False, ticks=True, tickColor="#1E293B", domainColor="#1E293B"),
            ),
            yOffset=alt.YOffset("axis:N", sort=["Log-Ratio Balance", "Co-Abundance Product"]),
            x=alt.X(
                "log10_count:Q",
                title="log10(1 + Sig. Continuous Pairs)",
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B"),
            ),
            color=alt.Color(
                "axis:N",
                scale=alt.Scale(
                    domain=["Log-Ratio Balance", "Co-Abundance Product"],
                    range=[OKABE_ORANGE, OKABE_REDDISH_PURPLE],
                ),
                title="Continuous Axis",
                legend=alt.Legend(orient="bottom", direction="vertical"),
            ),
        )
        .properties(width=200, height=290, title="b. Continuous Balance & Coupling")
    )

    panel_c = (
        alt.Chart(df_n)
        .mark_bar(height=9, cornerRadius=0)
        .encode(
            y=alt.Y(
                "cohort_label:N",
                sort=cohort_order,
                title="",
                axis=alt.Axis(labels=False, ticks=True, tickColor="#1E293B", domainColor="#1E293B"),
            ),
            yOffset=alt.YOffset(
                "novelty:N",
                sort=["De Novo (Both Univariate NS)", "Univariate-Anchored Synergy"],
            ),
            x=alt.X(
                "log10_count:Q",
                title="log10(1 + Sig. Pairs by Novelty)",
                axis=alt.Axis(tickColor="#1E293B", domainColor="#1E293B"),
            ),
            color=alt.Color(
                "novelty:N",
                scale=alt.Scale(
                    domain=["De Novo (Both Univariate NS)", "Univariate-Anchored Synergy"],
                    range=[OKABE_SKY_BLUE, "#475569"],
                ),
                title="Univariate Milo Overlap",
                legend=alt.Legend(orient="bottom", direction="vertical"),
            ),
        )
        .properties(width=200, height=290, title="c. Emergent De Novo vs. Anchored")
    )

    combined = alt.hconcat(panel_a, panel_b, panel_c, spacing=18).resolve_scale(color="independent")
    return Success(apply_nature_methods_theme(combined))


def build_recurrent_celltype_interactions_chart(
    bundles: Sequence[CohortPairwiseBundle],
    recurrent_df: pl.DataFrame,
    max_pairs: int = 18,
) -> Result[alt.Chart, str]:
    """Builds a Nature Methods dot-matrix chart of top lineage-pair interactions across cohorts."""
    top_labels = recurrent_df.head(max_pairs)["pair_label"].to_list()
    if not top_labels:
        empty_df = pl.DataFrame(
            {
                "cohort_label": ["No Significant Pairs"],
                "pair_label": ["No Significant Lineage Pairs"],
                "log10_pairs": [0.0],
                "direction": [EnrichmentDirection.RESPONDER.value],
            }
        )
        empty_chart = (
            alt.Chart(empty_df)
            .mark_circle(size=10, opacity=0.0)
            .encode(x="cohort_label:N", y="pair_label:N")
            .properties(width=320, height=200, title="Lineage-Pair Ecosystem Recurrence")
        )
        return Success(apply_nature_methods_theme(empty_chart))

    top_set = set(top_labels)
    per_cohort_rows: list[dict[str, object]] = []
    for bundle in bundles:
        c_label = bundle.summary.cohort_label
        sub_rec = aggregate_recurrent_celltype_pairs((bundle,))
        for r in sub_rec.to_dicts():
            p_lbl = str(r["pair_label"])
            if p_lbl in top_set:
                cnt = int(r["total_sig_pairs"])
                per_cohort_rows.append(
                    {
                        "cohort_label": c_label,
                        "pair_label": p_lbl,
                        "n_pairs": cnt,
                        "log10_pairs": float(np.log10(1.0 + float(cnt))),
                        "direction": str(r["direction"]),
                    }
                )

    plot_df = pl.DataFrame(per_cohort_rows)
    cohort_order = [b.summary.cohort_label for b in bundles]

    chart = (
        alt.Chart(plot_df)
        .mark_circle(opacity=0.88, stroke="#1E293B", strokeWidth=0.5)
        .encode(
            x=alt.X(
                "cohort_label:N",
                sort=cohort_order,
                title="Clinical Single-Cell Cohort",
                axis=alt.Axis(labelAngle=-35, labelLimit=200, tickColor="#1E293B", domainColor="#1E293B"),
            ),
            y=alt.Y(
                "pair_label:N",
                sort=top_labels,
                title="",
                axis=alt.Axis(labelLimit=380, tickColor="#1E293B", domainColor="#1E293B"),
            ),
            size=alt.Size(
                "log10_pairs:Q",
                scale=alt.Scale(range=[30, 260]),
                title="log10(1 + Pairs)",
            ),
            color=alt.Color(
                "direction:N",
                scale=alt.Scale(
                    domain=[EnrichmentDirection.RESPONDER.value, EnrichmentDirection.NON_RESPONDER.value],
                    range=[OKABE_BLUE, OKABE_VERMILION],
                ),
                title="Enriched Arm",
            ),
        )
        .properties(
            width=300,
            height=max(180, 18 * len(top_labels)),
            title="Cross-Cohort Lineage-Pair Ecosystem Recurrence",
        )
    )
    return Success(apply_nature_methods_theme(chart))


# =============================================================================
# 5. Pure Markdown Master Report Builder
# =============================================================================


def render_markdown_summary_report(
    summary_df: pl.DataFrame,
    recurrent_df: pl.DataFrame,
    bundles: Sequence[CohortPairwiseBundle],
    config: CohortPairwiseRunConfig,
) -> str:
    """Renders a comprehensive Markdown report comparing univariate vs pairwise Milo across cohorts."""
    total_cells = int(summary_df["n_cells"].to_numpy().sum()) if not summary_df.is_empty() else 0
    total_nhoods = int(summary_df["n_total_nhoods"].to_numpy().sum()) if not summary_df.is_empty() else 0
    total_eval_pairs = int(summary_df["n_evaluated_pairs"].to_numpy().sum()) if not summary_df.is_empty() else 0
    total_quad_sig = int(summary_df["n_sig_quadrant_total"].to_numpy().sum()) if not summary_df.is_empty() else 0
    total_cont_sig = int(summary_df["n_sig_continuous_total"].to_numpy().sum()) if not summary_df.is_empty() else 0
    total_de_novo = (
        int(
            summary_df["n_sig_quad_de_novo"].to_numpy().sum()
            + summary_df["n_sig_cont_de_novo"].to_numpy().sum()
        )
        if not summary_df.is_empty()
        else 0
    )
    total_modules = int(summary_df["n_modules"].to_numpy().sum()) if not summary_df.is_empty() else 0

    lines: list[str] = [
        "# Pairwise Milopy Across Nine Clinical Single-Cell Cohorts: Master Report",
        "",
        f"- **Pairwise Method**: `{config.pairwise_method.value}` (`max_jaccard_overlap={config.max_jaccard_overlap}`, `presence_mode={config.presence_mode.value}`)",
        f"- **FDR Threshold**: `{config.fdr_threshold}` (Spatial 2D Weighted BH + Berger's Intersection–Union Test)",
        f"- **Analyzed Totals**: `{total_cells:,}` cells | `{total_nhoods:,}` neighborhoods | `{total_eval_pairs:,}` disjoint neighborhood pairs",
        f"- **Irreducible Pairwise Discoveries**: `{total_quad_sig:,}` 2x2 quadrant hits (`{total_modules:,}` ecosystem modules) + `{total_cont_sig:,}` continuous balance/co-abundance hits (`{total_de_novo:,}` emergent *de novo* pairs with both endpoints univariate-NS)",
        "",
        "## 1. Cross-Cohort Univariate vs. Pairwise Execution Summary",
        "",
        "| Accession | Cohort | Evals (R / NR) | Nhoods (Act / Tot) | Univ. Sig (R / NR) | Disjoint Pairs | Quad Sig (`++` / `+-,-+` / `--`) | Cont Sig (Bal / Prod) | De Novo Pairs | Modules | Time (s) |",
        "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
    ]

    for r in summary_df.to_dicts():
        gated = int(r["n_sig_pm"]) + int(r["n_sig_mp"])
        de_novo = int(r["n_sig_quad_de_novo"]) + int(r["n_sig_cont_de_novo"])
        lines.append(
            f"| `{r['accession']}` | {r['cohort_label']} | "
            f"{r['n_responders']} / {r['n_non_responders']} | "
            f"{r['n_active_nhoods']:,} / {r['n_total_nhoods']:,} | "
            f"{r['n_univariate_sig_r']:,} / {r['n_univariate_sig_nr']:,} | "
            f"{r['n_evaluated_pairs']:,} | "
            f"**{r['n_sig_quadrant_total']:,}** ({r['n_sig_pp']:,} / {gated:,} / {r['n_sig_mm']:,}) | "
            f"**{r['n_sig_continuous_total']:,}** ({r['n_sig_log_ratio']:,} / {r['n_sig_co_abundance']:,}) | "
            f"{de_novo:,} | {r['n_modules']:,} | {float(r['runtime_seconds']):.1f} |"
        )

    lines.extend(["", "## 2. Top Recurrent Cell-Type Pair Regimes Across Cohorts", ""])
    if recurrent_df.is_empty():
        lines.append("_No significant lineage-pair interactions detected at the current threshold._")
    else:
        lines.extend(
            [
                "| Primary Lineage | Secondary Lineage | Pairwise Regime | Enriched Arm | Cohorts (`n`) | Total Sig. Pairs | De Novo Pairs | Mean Effect |",
                "| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: |",
            ]
        )
        for r in recurrent_df.head(20).to_dicts():
            lines.append(
                f"| {r['lineage_primary']} | {r['lineage_secondary']} | `{r['regime']}` | "
                f"`{r['direction']}` | **{r['n_cohorts']}** ({r['cohorts']}) | "
                f"{r['total_sig_pairs']:,} | {r['total_de_novo_pairs']:,} | "
                f"{float(r['mean_effect_size']):.2f} |"
            )

    lines.extend(["", "## 3. Top Emergent De Novo Pairwise Hits (`both_univariate_ns`)", ""])
    de_novo_q_frames = [
        b.quadrant_df.filter(
            pl.col("is_significant")
            & (pl.col("novelty_class") == PairNoveltyClass.BOTH_UNIVARIATE_NS.value)
        )
        for b in bundles
    ]
    non_empty_dn = [f for f in de_novo_q_frames if not f.is_empty()]
    if non_empty_dn:
        top_dn = (
            pl.concat(non_empty_dn)
            .sort(["fdr_iut", "pvalue_iut", "pvalue_quadrant"])
            .head(15)
        )
        lines.extend(
            [
                "| Accession | Pair `(u, v)` | Lineages `(u × v)` | Quadrant | Enriched Arm | `logOR` | `FDR_IUT` | Univ `FDR_u` | Univ `FDR_v` |",
                "| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |",
            ]
        )
        for r in top_dn.to_dicts():
            lines.append(
                f"| `{r['accession']}` | `({r['nhood_u']}, {r['nhood_v']})` | "
                f"{r['cell_type_u']} × {r['cell_type_v']} | `{r['quadrant']}` | "
                f"`{r['direction']}` | {float(r['log_or_quadrant']):+.2f} | "
                f"{float(r['fdr_iut']):.4f} | {float(r['univariate_fdr_u']):.3f} | "
                f"{float(r['univariate_fdr_v']):.3f} |"
            )
    else:
        lines.append("_No emergent quadrant de novo pairs detected; inspect `continuous_pairs.parquet` for compositional de novo balances._")

    lines.append("")
    return "\n".join(lines)


# =============================================================================
# 6. Imperative Shell — Filesystem Loaders, Exporters & CLI Orchestration
# =============================================================================


def load_cached_cohort_from_disk(accession: str, input_dir: Path) -> Result[CachedCohortInputs, str]:
    """Side-effecting boundary loader reading nhood_matrix.npz, cell_level_scores, and da_results."""
    cohort_dir = input_dir / accession
    mat_path = cohort_dir / "nhood_matrix.npz"
    cells_path = cohort_dir / "cell_level_scores.parquet"
    da_path = cohort_dir / "da_results.parquet"

    for p in (mat_path, cells_path, da_path):
        if not p.exists():
            return Failure(f"[{accession}] Missing cached file: {p}")

    try:
        nhoods_mat = sp.csc_matrix(sp.load_npz(mat_path))
        cells_df = pl.read_parquet(cells_path)
        da_df = pl.read_parquet(da_path)
    except Exception as exc:
        return Failure(f"[{accession}] Failed reading cached artifacts: {exc}")

    return Success(
        CachedCohortInputs(
            accession=accession,
            nhoods_mat=nhoods_mat,
            cells_df=cells_df,
            univariate_da_df=da_df,
        )
    )


def export_figure_to_disk(chart: AltairChart, base_path: Path) -> Result[tuple[Path, Path], str]:
    """Side-effecting boundary exporting an Altair chart to vector SVG and 300 DPI PNG."""
    try:
        with alt.data_transformers.disable_max_rows():
            paths = export_altair_figure(chart, base_path, scale=2.5, max_vector_points=1500)
        return Success(paths)
    except Exception as exc:
        return Failure(f"Failed exporting figure {base_path}: {exc}")


def run_all_cohorts_pipeline(
    config: CohortPairwiseRunConfig,
    loader: CohortDataLoader = load_cached_cohort_from_disk,
    exporter: FigureExporter = export_figure_to_disk,
    persist_outputs: bool = True,
) -> Result[tuple[tuple[CohortPairwiseBundle, ...], pl.DataFrame, pl.DataFrame], str]:
    """Executes pairwise Milopy across all resolved cohorts with injected IO boundaries."""
    cohorts = resolve_target_cohorts(config.cohort_filter)
    bundles: list[CohortPairwiseBundle] = []

    for acc in cohorts:
        t_start = time.perf_counter()
        loaded_res = loader(acc, config.input_dir)
        match loaded_res:
            case Failure(err):
                return Failure(err)
            case Success(inputs):
                pass

        bundle_res = analyze_single_cohort(inputs, config, runtime_seconds=0.0)
        match bundle_res:
            case Failure(err):
                return Failure(err)
            case Success(raw_bundle):
                elapsed = time.perf_counter() - t_start
                timed_summary = raw_bundle.summary.model_copy(update={"runtime_seconds": elapsed})
                bundle = CohortPairwiseBundle(
                    summary=timed_summary,
                    quadrant_df=raw_bundle.quadrant_df,
                    continuous_df=raw_bundle.continuous_df,
                    modules_df=raw_bundle.modules_df,
                )
                bundles.append(bundle)

        if persist_outputs:
            cohort_out = config.output_dir / acc
            cohort_out.mkdir(parents=True, exist_ok=True)
            bundle.quadrant_df.write_parquet(cohort_out / "quadrant_pairs.parquet")
            bundle.continuous_df.write_parquet(cohort_out / "continuous_pairs.parquet")
            bundle.modules_df.write_parquet(cohort_out / "pair_modules.parquet")

            if config.render_plots:
                vol_res = build_cohort_quadrant_volcano_chart(
                    quadrant_df=bundle.quadrant_df,
                    summary=bundle.summary,
                    fdr_threshold=config.fdr_threshold,
                )
                match vol_res:
                    case Failure(err):
                        return Failure(err)
                    case Success(vol_chart):
                        exp_res = exporter(vol_chart, cohort_out / "pairwise_quadrant_volcano")
                        if isinstance(exp_res, Failure):
                            return Failure(exp_res.failure())

        s = bundle.summary
        print(
            f"[{acc}] {s.n_responders}R/{s.n_non_responders}NR | "
            f"pairs={s.n_evaluated_pairs:,} | quad_sig={s.n_sig_quadrant_total:,} "
            f"(++:{s.n_sig_pp}, +-/-+:{s.n_sig_pm + s.n_sig_mp}, --:{s.n_sig_mm}) | "
            f"cont_sig={s.n_sig_continuous_total:,} | modules={s.n_modules} | {s.runtime_seconds:.2f}s",
            flush=True,
        )

    bundle_tuple = tuple(bundles)
    summary_df = build_cohorts_summary_dataframe(bundle_tuple, config.sort_order)
    recurrent_df = aggregate_recurrent_celltype_pairs(bundle_tuple)

    if persist_outputs:
        reports_dir = config.output_dir / "reports"
        reports_dir.mkdir(parents=True, exist_ok=True)
        summary_df.write_parquet(reports_dir / "pairwise_cohorts_summary.parquet")
        recurrent_df.write_parquet(reports_dir / "recurrent_celltype_pairs.parquet")

        md_text = render_markdown_summary_report(summary_df, recurrent_df, bundle_tuple, config)
        (reports_dir / "milopy_pairwise_cohorts_summary.md").write_text(md_text, encoding="utf-8")

        if config.render_plots:
            comp_res = build_cross_cohort_pairwise_comparison_chart(summary_df)
            match comp_res:
                case Failure(err):
                    return Failure(err)
                case Success(comp_chart):
                    exp_c = exporter(comp_chart, reports_dir / "pairwise_cohort_comparison")
                    if isinstance(exp_c, Failure):
                        return Failure(exp_c.failure())

            rec_res = build_recurrent_celltype_interactions_chart(bundle_tuple, recurrent_df)
            match rec_res:
                case Failure(err):
                    return Failure(err)
                case Success(rec_chart):
                    exp_r = exporter(rec_chart, reports_dir / "pairwise_celltype_interactions")
                    if isinstance(exp_r, Failure):
                        return Failure(exp_r.failure())

    return Success((bundle_tuple, summary_df, recurrent_df))


def parse_cli_args(argv: Sequence[str]) -> CohortPairwiseRunConfig:
    """Parses command-line flags into a frozen CohortPairwiseRunConfig."""
    parser = argparse.ArgumentParser(
        description="Apply milopy.pairwise across the nine pre-analyzed clinical single-cell cohorts."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("output/results/milopy_response_analysis"),
        help="Directory containing per-cohort cached Milo outputs",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/results/milopy_pairwise_analysis"),
        help="Destination directory for pairwise Parquet tables, reports, and SVG figures",
    )
    parser.add_argument(
        "--cohorts",
        type=str,
        default="",
        help="Optional comma-separated subset of cohort accessions (empty = all 9 cohorts)",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.10,
        help="Spatial weighted BH FDR cutoff (default: 0.10)",
    )
    parser.add_argument(
        "--max-jaccard-overlap",
        type=float,
        default=0.0,
        help="Maximum allowed cell-sharing Jaccard index between paired neighborhoods (default: 0.0)",
    )
    parser.add_argument(
        "--presence-mode",
        type=str,
        choices=[m.value for m in PresenceMode],
        default=PresenceMode.HARD_DEPTH_GATED.value,
        help="Depth-calibrated presence estimator mode",
    )
    parser.add_argument(
        "--method",
        type=str,
        choices=[m.value for m in PairwiseMethod],
        default=PairwiseMethod.HYBRID_ALL.value,
        help="Pairwise testing family (default: hybrid_all)",
    )
    parser.add_argument(
        "--min-donors-per-quadrant",
        type=float,
        default=3.0,
        help="Minimum donors in enriched arm for quadrant calls (default: 3.0)",
    )
    parser.add_argument(
        "--adaptive-small-cohort-gate",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Adapt min_donors_per_quadrant to 2.0 when min(n_R, n_NR) == 2 (default: True)",
    )
    parser.add_argument(
        "--max-refine-pairs",
        type=int,
        default=250,
        help="Top IUT-ranked candidate pairs refined via Firth logistic LRT per cohort",
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=0,
        help="Stratified donor-label permutations on retained candidates (default: 0)",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=512,
        help="Upper-triangular BLAS block size (default: 512)",
    )
    parser.add_argument(
        "--sort-order",
        type=str,
        choices=[s.value for s in CohortSortOrder],
        default=CohortSortOrder.PANEL_ORDER.value,
        help="Cohort ordering in summary outputs",
    )
    parser.add_argument(
        "--plot",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Export Nature Methods vector SVG figures (default: True)",
    )

    ns = parser.parse_args(list(argv))
    cohort_tokens = tuple(tok.strip() for tok in str(ns.cohorts).split(",") if tok.strip())
    cohort_maybe: Maybe[tuple[str, ...]] = Some(cohort_tokens) if cohort_tokens else Nothing

    return CohortPairwiseRunConfig(
        input_dir=ns.input_dir,
        output_dir=ns.output_dir,
        cohort_filter=cohort_maybe,
        fdr_threshold=float(ns.fdr_threshold),
        max_jaccard_overlap=float(ns.max_jaccard_overlap),
        presence_mode=PresenceMode(str(ns.presence_mode)),
        pairwise_method=PairwiseMethod(str(ns.method)),
        min_donors_per_quadrant=float(ns.min_donors_per_quadrant),
        adaptive_small_cohort_gate=bool(ns.adaptive_small_cohort_gate),
        max_refine_pairs=int(ns.max_refine_pairs),
        n_permutations=int(ns.n_permutations),
        chunk_size=int(ns.chunk_size),
        sort_order=CohortSortOrder(str(ns.sort_order)),
        render_plots=bool(ns.plot),
    )


def run_cli(argv: Sequence[str]) -> int:
    """Imperative entry point returning POSIX exit code."""
    cfg = parse_cli_args(argv)
    match run_all_cohorts_pipeline(cfg):
        case Failure(err):
            print(f"Pairwise cohort pipeline failed: {err}", file=sys.stderr)
            return 1
        case Success((_, summary_df, _)):
            print(
                f"\nCompleted pairwise Milopy across {summary_df.height} cohorts -> {cfg.output_dir}",
                flush=True,
            )
            return 0
        case _ as unreachable:
            _ = unreachable
            return 1


if __name__ == "__main__":
    sys.exit(run_cli(sys.argv[1:]))
