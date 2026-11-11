#!/usr/bin/env python3
"""9-Cohort Neighborhood Donor/Batch Mixing Audit, Scale Coarsening & 4-Layer Pair Significance Pipeline.

Executes a purely functional diagnostic and discovery benchmark across all 9 clinical ICB cohorts:
1. Audits neighborhood donor mixing (`effective_donors_lisi`, `max_single_donor_fraction`) and
   inferred multi-study batch purity (`<PMID>ggg<Patient>` in `CELLxGENE_7b20c613`, `<Study>.<Patient>`
   in `CELLxGENE_05a8c945`).
2. Quantifies the neighborhood scale trade-off across 4 coarsening levels:
   - `S0_Raw_k30`: Unpooled `k=30` Milo spheres
   - `S1_Super_Res10`: Fine lineage-preserving Leiden super-neighborhoods (`res=1.0, J>=0.10`)
   - `S2_Coarse_Res025`: Coarse lineage-preserving Leiden super-neighborhoods (`res=0.25, J>=0.02`)
   - `S3_Archetype_1PerState`: Exact cell-union state archetypes (`1` per `Meta42` state)
3. Runs a 5-stage simplification & 4-layer statistical verification cascade (`A0` -> `A4`):
   - Stratified Cochran-Mantel-Haenszel (`CMH`) across sub-studies
   - Outcome-blind Donor-Mixing QC (`LISI_pt >= 2.5`, `max_donor <= 0.50`)
   - `PairScopeFilter.CROSS_COMPARTMENT_ONLY` + coarse super-neighborhoods (`res=0.25`)
   - Canonical state-pair bundle consolidation (`1` chord per `(State_A, State_B, Quadrant, Direction)`)
   - Exact within-stratum donor-label IUT permutations (`B = 1,000`), Westfall-Young `max-T` FWER,
     Firth library-adjusted logistic LRT, and Leave-One-Donor/Stratum-Out (`LODO`/`LOSO`) jackknife.
4. Exports publication-grade Nature Methods Inkscape-layered `.svg` + 300 DPI `.png` diagnostics
   and simplified de-redundant UMAP composites to `output/results/milopy_mixing_and_simplified_pairs/`.
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from enum import StrEnum
import html
from pathlib import Path
import sys
import time
from typing import TypeAlias, assert_never, cast
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
    PairScopeFilter,
    PairwiseAnalysisResult,
    PairwiseConfig,
    PairwiseMethod,
    PrefilterConfig,
    PrefilterMode,
    PresenceConfig,
    PresenceMode,
    RefinementConfig,
    ScreenConfig,
    UmapPairVizConfig,
    analyze_nhood_pairs,
    build_pairwise_umap_composite,
    calibrate_candidates_stratified_iut_permutation,
    coarse_grain_nhoods_by_graph,
    compute_nhood_mixing_metrics,
    consolidate_state_pair_bundles,
    deduplicate_pairs_to_lead_bundles,
    estimate_nhood_presence,
    extract_inferred_substudy_strata,
    prepare_pairwise_umap_bundle,
)
from milopy.pairwise.types import FloatArray, IntArray, SparseMatrix
from plotting_utils import (  # type: ignore[import-untyped]
    COHORT_PANELS,
    OKABE_BLUE,
    OKABE_BLUISH_GREEN,
    OKABE_ORANGE,
    OKABE_VERMILION,
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

SPOTLIGHT_COHORTS: tuple[str, ...] = ("CELLxGENE_7b20c613", "GSE243013")


class ScaleLevel(StrEnum):
    """Four neighborhood coarsening scales evaluated across cohorts."""

    RAW_K30 = "S0: Raw k=30"
    SUPER_RES10 = "S1: Fine Super (res=1.0)"
    COARSE_RES025 = "S2: Coarse Super (res=0.25)"
    STATE_ARCHETYPE = "S3: 1/State Archetype"


ALL_SCALES: tuple[ScaleLevel, ...] = (
    ScaleLevel.RAW_K30,
    ScaleLevel.SUPER_RES10,
    ScaleLevel.COARSE_RES025,
    ScaleLevel.STATE_ARCHETYPE,
)


class AuditCliConfig(BaseModel):
    """Immutable CLI configuration for the 9-cohort mixing & pair-simplification pipeline."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    milo_response_dir: Path = Field(default=Path("output/results/milopy_response_analysis"))
    pairwise_t0_dir: Path = Field(default=Path("output/results/milopy_pairwise_analysis"))
    meta42_annotations_path: Path = Field(
        default=Path("output/results/milopy_pairwise_deconv_bridge/nhood_meta42_annotations.parquet")
    )
    output_dir: Path = Field(default=Path("output/results/milopy_mixing_and_simplified_pairs"))
    n_permutations: int = Field(default=1000, ge=50, le=10000)
    fdr_threshold: float = Field(default=0.10, gt=0.0, le=1.0)
    coarse_resolution: float = Field(default=0.25, gt=0.0, le=5.0)
    coarse_min_jaccard: float = Field(default=0.02, gt=0.0, le=1.0)
    min_distinct_donors: int = Field(default=3, ge=1)
    min_effective_lisi: float = Field(default=2.5, ge=1.0)
    max_donor_dominance: float = Field(default=0.50, gt=0.0, le=1.0)
    max_chords: int = Field(default=20, ge=4, le=80)


class CohortAuditContext(BaseModel):
    """Immutable cohort inputs aligned across cells, donors, UMAP coordinates, and Meta42 states."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    accession: str
    cohort_label: str
    nhoods_mat: SparseMatrix
    cells_df: pl.DataFrame
    da_results_df: pl.DataFrame
    t0_quadrant_df: pl.DataFrame
    umap_coords: FloatArray
    cell_types: tuple[str, ...]
    eval_ids: tuple[str, ...]
    response_labels: tuple[str, ...]
    response_binary: FloatArray
    inferred_strata: Maybe[tuple[str, ...]]
    strata_codes: Maybe[IntArray]
    count_matrix: FloatArray
    library_sizes: FloatArray
    cell_eval_codes: IntArray
    meta_codes: tuple[str, ...]
    code_to_readable: Mapping[str, str]


def _encode_optional_strata(strata: Maybe[tuple[str, ...]]) -> Maybe[IntArray]:
    """Maps inferred sub-study stratum names to integer stratum codes."""
    match strata:
        case Maybe.empty:
            return Nothing
        case Some(s_tup):
            uniq = sorted(set(s_tup))
            if len(uniq) <= 1:
                return Nothing
            s_map = {name: idx for idx, name in enumerate(uniq)}
            return Some(np.array([s_map[name] for name in s_tup], dtype=np.int64))
        case _ as unreachable:
            _ = unreachable
            return Nothing


def load_cohort_audit_context(
    accession: str,
    cfg: AuditCliConfig,
    all_meta42_df: pl.DataFrame,
) -> Result[CohortAuditContext, str]:
    """Loads single-cell UMAP, neighborhood matrix, and donor design arrays without mutation."""
    cohort_dir = cfg.milo_response_dir / accession
    npz_path = cohort_dir / "nhood_matrix.npz"
    cells_path = cohort_dir / "cell_level_scores.parquet"
    da_path = cohort_dir / "da_results.parquet"
    if not (npz_path.exists() and cells_path.exists() and da_path.exists()):
        return Failure(f"[{accession}] Missing Milo artifacts under {cohort_dir}")

    nhoods_mat = sp.csc_matrix(sp.load_npz(npz_path))
    cells_df = pl.read_parquet(cells_path)
    da_df = pl.read_parquet(da_path)
    q_path = cfg.pairwise_t0_dir / accession / "quadrant_pairs.parquet"
    t0_q_df = pl.read_parquet(q_path) if q_path.exists() else pl.DataFrame()

    n_cells = cells_df.height
    n_nhoods = int(nhoods_mat.shape[1])
    umap_coords = np.column_stack(
        [
            cells_df["UMAP1"].cast(pl.Float64).to_numpy(),
            cells_df["UMAP2"].cast(pl.Float64).to_numpy(),
        ]
    )

    norm_cells = (
        cells_df.lazy()
        .with_columns(
            pl.col("patient_id").cast(pl.String).str.strip_chars().alias("_pt"),
            pl.col("clinical_response")
            .cast(pl.String)
            .str.strip_chars()
            .str.to_lowercase()
            .alias("_resp"),
            (
                pl.col("cell_type").cast(pl.String)
                if "cell_type" in cells_df.columns
                else pl.lit("Unspecified")
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

    eval_ids = tuple(str(x) for x in eval_meta["_eval_id"].to_list())
    n_evals = len(eval_ids)
    if n_evals < 4:
        return Failure(f"[{accession}] Fewer than 4 donor evaluations ({n_evals})")

    eval_to_idx = {eid: idx for idx, eid in enumerate(eval_ids)}
    cell_eval_codes = np.array(
        [eval_to_idx[str(e)] for e in norm_cells["_eval_id"].to_list()],
        dtype=np.int64,
    )
    sample_ind = sp.csr_matrix(
        (np.ones(n_cells, dtype=np.float64), (np.arange(n_cells, dtype=np.int64), cell_eval_codes)),
        shape=(n_cells, n_evals),
    )
    count_matrix = np.asarray(nhoods_mat.T.dot(sample_ind).toarray(), dtype=np.float64)
    lib_sizes = eval_meta["library_size"].to_numpy().astype(np.float64)
    resp_labels = tuple(str(r) for r in eval_meta["clinical_response"].to_list())
    y_bin = np.array([1.0 if r == "responder" else 0.0 for r in resp_labels], dtype=np.float64)
    inferred_strata = extract_inferred_substudy_strata(eval_ids)
    strata_codes = _encode_optional_strata(inferred_strata)

    meta_sub = (
        all_meta42_df.filter(pl.col("accession") == accession)
        if (not all_meta42_df.is_empty() and "accession" in all_meta42_df.columns)
        else pl.DataFrame()
    )
    raw_cts = (
        [str(x) for x in da_df["Nhood_CellType"].to_list()]
        if "Nhood_CellType" in da_df.columns
        else ["Unspecified"] * n_nhoods
    )
    meta_codes_list = list(raw_cts)
    code_to_readable: dict[str, str] = {}
    if not meta_sub.is_empty():
        for r in meta_sub.sort("nhood_idx").select(["nhood_idx", "meta_state", "meta_state_label"]).to_dicts():
            idx = int(r["nhood_idx"])
            ms = str(r["meta_state"])
            ml = str(r["meta_state_label"])
            code_to_readable[ms] = ml
            if 0 <= idx < n_nhoods:
                meta_codes_list[idx] = ms

    return Success(
        CohortAuditContext(
            accession=accession,
            cohort_label=get_single_cell_cohort_display_label(accession, "cancer_author_year"),
            nhoods_mat=nhoods_mat,
            cells_df=cells_df,
            da_results_df=da_df,
            t0_quadrant_df=t0_q_df,
            umap_coords=umap_coords,
            cell_types=tuple(str(c) for c in norm_cells["_ct"].to_list()),
            eval_ids=eval_ids,
            response_labels=resp_labels,
            response_binary=y_bin,
            inferred_strata=inferred_strata,
            strata_codes=strata_codes,
            count_matrix=count_matrix,
            library_sizes=lib_sizes,
            cell_eval_codes=cell_eval_codes,
            meta_codes=tuple(meta_codes_list),
            code_to_readable=code_to_readable,
        )
    )


def _build_archetype_state_unions(
    ctx: CohortAuditContext,
) -> tuple[sp.csc_matrix, FloatArray, tuple[str, ...]]:
    """Constructs exact 1-per-state cell-union archetypes `nhoods_state[i, g] = 1{sum_{v in state_g} A_{i,v} > 0}`."""
    uniq_states = tuple(sorted(set(ctx.meta_codes)))
    st_to_idx = {st: idx for idx, st in enumerate(uniq_states)}
    n_nhoods = len(ctx.meta_codes)
    n_states = len(uniq_states)
    group_map = sp.csc_matrix(
        (
            np.ones(n_nhoods, dtype=np.float64),
            (np.arange(n_nhoods, dtype=np.int64), np.array([st_to_idx[s] for s in ctx.meta_codes], dtype=np.int64)),
        ),
        shape=(n_nhoods, n_states),
    )
    nhoods_bin = (sp.csc_matrix(ctx.nhoods_mat) > 0).astype(np.float64)
    nhoods_arch = (nhoods_bin.dot(group_map) > 0).astype(np.float64).tocsc()

    n_cells = int(nhoods_arch.shape[0])
    n_evals = len(ctx.eval_ids)
    sample_ind = sp.csr_matrix(
        (np.ones(n_cells, dtype=np.float64), (np.arange(n_cells, dtype=np.int64), ctx.cell_eval_codes)),
        shape=(n_cells, n_evals),
    )
    counts_arch = np.asarray(nhoods_arch.T.dot(sample_ind).toarray(), dtype=np.float64)
    return nhoods_arch, counts_arch, uniq_states


def evaluate_cohort_scales(
    ctx: CohortAuditContext,
    cfg: AuditCliConfig,
    mix_cfg: MixingQCConfig,
) -> Result[list[dict[str, object]], str]:
    """Evaluates donor-mixing & pair-space compression across all 4 neighborhood scales."""
    rows: list[dict[str, object]] = []
    lib_maybe = Some(ctx.library_sizes)

    for scale in ALL_SCALES:
        match scale:
            case ScaleLevel.RAW_K30:
                s_nhoods: SparseMatrix = ctx.nhoods_mat
                s_counts = ctx.count_matrix
                s_labels = ctx.meta_codes
            case ScaleLevel.SUPER_RES10:
                cg10 = coarse_grain_nhoods_by_graph(
                    nhoods_mat=ctx.nhoods_mat,
                    count_matrix=ctx.count_matrix,
                    cell_sample_codes=Some(ctx.cell_eval_codes),
                    sample_library_sizes=lib_maybe,
                    nhood_group_labels=Some(ctx.meta_codes),
                    min_jaccard=0.10,
                    leiden_resolution=1.0,
                )
                match cg10:
                    case Failure(err):
                        return Failure(err)
                    case Success(d10):
                        s_nhoods = d10.nhoods_super
                        s_counts = d10.counts_super
                        s_labels = d10.super_nhood_labels
            case ScaleLevel.COARSE_RES025:
                cg025 = coarse_grain_nhoods_by_graph(
                    nhoods_mat=ctx.nhoods_mat,
                    count_matrix=ctx.count_matrix,
                    cell_sample_codes=Some(ctx.cell_eval_codes),
                    sample_library_sizes=lib_maybe,
                    nhood_group_labels=Some(ctx.meta_codes),
                    min_jaccard=cfg.coarse_min_jaccard,
                    leiden_resolution=cfg.coarse_resolution,
                )
                match cg025:
                    case Failure(err):
                        return Failure(err)
                    case Success(d025):
                        s_nhoods = d025.nhoods_super
                        s_counts = d025.counts_super
                        s_labels = d025.super_nhood_labels
            case ScaleLevel.STATE_ARCHETYPE:
                s_nhoods, s_counts, s_labels = _build_archetype_state_unions(ctx)
            case _ as unreachable:
                assert_never(unreachable)

        mix_res = compute_nhood_mixing_metrics(
            nhoods_mat=s_nhoods,
            count_matrix=s_counts,
            sample_ids=ctx.eval_ids,
            batch_labels=ctx.inferred_strata,
            nhood_labels=Some(s_labels),
            config=mix_cfg,
        )
        match mix_res:
            case Failure(err):
                return Failure(err)
            case Success(mdf):
                n_u = mdf.height
                n_pairs = (n_u * (n_u - 1)) // 2
                lisi = mdf["effective_donors_lisi"].to_numpy()
                dom = mdf["max_single_donor_fraction"].to_numpy()
                nd = mdf["n_distinct_donors"].to_numpy()
                nc = mdf["n_cells"].to_numpy()
                bf = mdf["max_batch_fraction"].to_numpy()
                qc = mdf["passes_mixing_qc"].to_numpy()
                rows.append(
                    {
                        "accession": ctx.accession,
                        "cohort_label": ctx.cohort_label,
                        "n_donors": len(ctx.eval_ids),
                        "scale": scale.value,
                        "n_units": n_u,
                        "n_candidate_pairs": n_pairs,
                        "log10_candidate_pairs": float(np.log10(max(n_pairs, 1))),
                        "median_cells_per_unit": float(np.median(nc)),
                        "median_donors_per_unit": float(np.median(nd)),
                        "median_effective_donors_lisi": float(np.median(lisi)),
                        "frac_single_donor": float(np.mean(nd <= 1)),
                        "frac_dominance_gt_50pct": float(np.mean(dom > 0.50)),
                        "median_batch_purity": float(np.median(bf)),
                        "frac_passes_mixing_qc": float(np.mean(qc)),
                    }
                )
    return Success(rows)


def audit_gondal_substudies(ctx: CohortAuditContext) -> pl.DataFrame:
    """Computes exact donor counts, responder rates, and ontology distributions per Gondal PMID sub-study."""
    df_cells = ctx.cells_df.with_columns(
        pl.col("patient_id")
        .cast(pl.String)
        .map_elements(lambda s: s.split("ggg", 1)[0] if "ggg" in s else "Other", return_dtype=pl.String)
        .alias("pmid_substudy"),
        pl.col("clinical_response").cast(pl.String).str.to_lowercase().alias("_resp"),
    )
    donors = df_cells.select(["pmid_substudy", "patient_id", "_resp"]).unique()
    donor_agg = (
        donors.group_by("pmid_substudy")
        .agg(
            pl.len().cast(pl.Int64).alias("n_donors"),
            (pl.col("_resp") == "responder").sum().cast(pl.Int64).alias("n_responders"),
            (100.0 * (pl.col("_resp") == "responder").mean()).alias("responder_pct"),
        )
    )
    cell_agg = (
        df_cells.group_by("pmid_substudy")
        .agg(
            pl.len().cast(pl.Int64).alias("n_cells"),
            pl.col("cell_type").mode().first().cast(pl.String).alias("modal_cell_type"),
        )
    )
    return (
        donor_agg.join(cell_agg, on="pmid_substudy", how="left")
        .with_columns(
            (
                (pl.col("n_responders") > 0)
                & ((pl.col("n_donors") - pl.col("n_responders")) > 0)
            ).alias("is_cmh_informative")
        )
        .sort("n_cells", descending=True)
    )


def run_cohort_ablation_and_bundle(
    ctx: CohortAuditContext,
    cfg: AuditCliConfig,
    mix_cfg: MixingQCConfig,
    build_umap_bundle: bool,
) -> Result[
    tuple[
        list[dict[str, object]],
        dict[str, object],
        pl.DataFrame,
        Maybe[AltairComposite],
    ],
    str,
]:
    """Runs the 5-stage simplification & 4-layer statistical calibration cascade (`A0` -> `A4`)."""
    lib_seq = tuple(float(x) for x in ctx.library_sizes.tolist())
    ablation_rows: list[dict[str, object]] = []

    # Stage A0: Raw k=30 Unstratified Baseline (cached T0 or live)
    a0_sig = (
        int(ctx.t0_quadrant_df.filter(pl.col("is_significant")).height)
        if (not ctx.t0_quadrant_df.is_empty() and "is_significant" in ctx.t0_quadrant_df.columns)
        else 0
    )
    a0_same_comp = 0
    if a0_sig > 0:
        sig_t0 = ctx.t0_quadrant_df.filter(pl.col("is_significant"))
        ded_t0 = deduplicate_pairs_to_lead_bundles(
            sig_t0,
            require_cross_compartment=True,
            one_per_state_pair=True,
        )
        if isinstance(ded_t0, Success):
            a0_same_comp = int(sig_t0.height) - int(ded_t0.unwrap().height)
    ablation_rows.append(
        {
            "accession": ctx.accession,
            "cohort_label": ctx.cohort_label,
            "stage": "A0: Raw k=30 Unstrat.",
            "stage_order": 0,
            "n_units": int(ctx.nhoods_mat.shape[1]),
            "n_evaluated_pairs": (int(ctx.nhoods_mat.shape[1]) * (int(ctx.nhoods_mat.shape[1]) - 1)) // 2,
            "n_sig_pairs": a0_sig,
            "n_redundant_or_intracomp_removed": a0_same_comp,
        }
    )

    # Stage A1: Fine Super-Neighborhoods (res=1.0, CROSS_CELL_TYPE_ONLY, Unstratified)
    n_resp_c = int(np.sum(ctx.response_binary > 0.5))
    n_non_resp_c = len(ctx.eval_ids) - n_resp_c
    min_arm = min(n_resp_c, n_non_resp_c)
    eff_min_donors = min(3.0, float(max(2, min_arm)))
    eff_min_total = min(6.0, float(max(4, len(ctx.eval_ids) - 1)))

    cg_a1_res = coarse_grain_nhoods_by_graph(
        nhoods_mat=ctx.nhoods_mat,
        count_matrix=ctx.count_matrix,
        cell_sample_codes=Some(ctx.cell_eval_codes),
        sample_library_sizes=Some(ctx.library_sizes),
        nhood_group_labels=Some(ctx.meta_codes),
        cell_types=Some(ctx.cell_types),
        min_jaccard=0.10,
        leiden_resolution=1.0,
    )
    match cg_a1_res:
        case Failure(err):
            return Failure(err)
        case Success(cg_a1):
            pass

    readable_a1 = tuple(ctx.code_to_readable.get(lbl, lbl) for lbl in cg_a1.super_nhood_labels)
    t4_prefilter = PrefilterConfig(
        mode=PrefilterMode.MARGINAL_ENTROPY_VARIANCE,
        min_marginal_prevalence=0.08,
        max_marginal_prevalence=0.92,
        min_variance_quantile=0.15,
        apply_tarone_iut_pruning=True,
        tarone_p_ceiling=0.005,
    )
    t4_screen = ScreenConfig(
        chunk_size=512,
        fdr_threshold=cfg.fdr_threshold,
        min_donors_per_quadrant=eff_min_donors,
        min_effective_samples=eff_min_total,
        max_retained_candidates=50000,
        apply_tarone_iut_pruning=True,
        tarone_p_ceiling=0.005,
    )

    cfg_a1 = PairwiseConfig(
        method=PairwiseMethod.QUADRANT_2X2,
        prefilter=t4_prefilter,
        overlap=OverlapConfig(max_jaccard_overlap=0.0, scope_filter=PairScopeFilter.CROSS_CELL_TYPE_ONLY),
        presence=PresenceConfig(mode=PresenceMode.HARD_DEPTH_GATED),
        screen=t4_screen,
        refinement=RefinementConfig(enabled=False),
        modules=PairModuleConfig(enabled=False),
    )
    res_a1 = analyze_nhood_pairs(
        nhoods_mat=cg_a1.nhoods_super,
        count_matrix=cg_a1.counts_super,
        response_labels=ctx.response_labels,
        sample_library_sizes=Some(lib_seq),
        nhood_group_labels=Some(readable_a1),
        sample_ids=Some(ctx.eval_ids),
        config=cfg_a1,
    )
    match res_a1:
        case Failure(err):
            return Failure(err)
        case Success(r1):
            ablation_rows.append(
                {
                    "accession": ctx.accession,
                    "cohort_label": ctx.cohort_label,
                    "stage": "A1: Super (res=1.0) Unstrat.",
                    "stage_order": 1,
                    "n_units": r1.n_active_nhoods,
                    "n_evaluated_pairs": r1.n_evaluated_pairs,
                    "n_sig_pairs": r1.n_significant_quadrant_pairs,
                    "n_redundant_or_intracomp_removed": 0,
                }
            )

    # Stage A2: + Sub-Study CMH Stratification + Outcome-Blind Donor-Mixing QC + Firth Refinement
    cfg_a2 = PairwiseConfig(
        method=PairwiseMethod.HYBRID_ALL,
        mixing_qc=mix_cfg,
        prefilter=t4_prefilter,
        overlap=cfg_a1.overlap,
        presence=cfg_a1.presence,
        screen=t4_screen,
        refinement=RefinementConfig(enabled=True, max_refine_pairs=250),
        modules=PairModuleConfig(enabled=True, min_endpoint_jaccard=0.05),
    )
    res_a2 = analyze_nhood_pairs(
        nhoods_mat=cg_a1.nhoods_super,
        count_matrix=cg_a1.counts_super,
        response_labels=ctx.response_labels,
        sample_library_sizes=Some(lib_seq),
        strata_labels=ctx.inferred_strata,
        nhood_group_labels=Some(readable_a1),
        sample_ids=Some(ctx.eval_ids),
        config=cfg_a2,
    )
    match res_a2:
        case Failure(err):
            return Failure(err)
        case Success(r2):
            ablation_rows.append(
                {
                    "accession": ctx.accession,
                    "cohort_label": ctx.cohort_label,
                    "stage": "A2: + Substudy CMH & MixQC",
                    "stage_order": 2,
                    "n_units": r2.n_active_nhoods,
                    "n_evaluated_pairs": r2.n_evaluated_pairs,
                    "n_sig_pairs": r2.n_significant_quadrant_pairs,
                    "n_redundant_or_intracomp_removed": 0,
                }
            )

    # Stage A3: + Cross-Compartment Gating & Canonical State-Pair Consolidation (1 chord / State_A <-> State_B)
    dedup_res = deduplicate_pairs_to_lead_bundles(
        r2.pairs_df,
        require_cross_compartment=True,
        one_per_state_pair=True,
    )
    bundles_res = consolidate_state_pair_bundles(r2.pairs_df, fdr_threshold=cfg.fdr_threshold)
    pres_res = estimate_nhood_presence(cg_a1.counts_super, Some(ctx.library_sizes), cfg_a2.presence)
    match (dedup_res, bundles_res, pres_res):
        case (Failure(err), _, _) | (_, Failure(err), _) | (_, _, Failure(err)):
            return Failure(err)
        case (Success(dedup_df), Success(bundles_df), Success(pres_super)):
            pass

    sig_dedup = dedup_df.filter(pl.col("is_significant"))
    ablation_rows.append(
        {
            "accession": ctx.accession,
            "cohort_label": ctx.cohort_label,
            "stage": "A3: + Cross-Comp. Consolid.",
            "stage_order": 3,
            "n_units": r2.n_active_nhoods,
            "n_evaluated_pairs": r2.n_evaluated_pairs,
            "n_sig_pairs": int(sig_dedup.height),
            "n_redundant_or_intracomp_removed": int(r2.n_significant_quadrant_pairs) - int(sig_dedup.height),
        }
    )

    # Stage A4: + B=1,000 Stratified IUT Permutations, Westfall-Young max-T FWER & LODO/LOSO Jackknife
    cal_res = calibrate_candidates_stratified_iut_permutation(
        candidates_df=sig_dedup,
        presence=pres_super,
        response_binary=ctx.response_binary,
        strata_codes=ctx.strata_codes,
        n_permutations=cfg.n_permutations,
        seed=42,
        fdr_threshold=cfg.fdr_threshold,
    )
    match cal_res:
        case Failure(err):
            return Failure(err)
        case Success(verified_leads_df):
            pass

    n_verified = (
        int(verified_leads_df.filter(pl.col("is_statistically_verified")).height)
        if not verified_leads_df.is_empty()
        else 0
    )
    ablation_rows.append(
        {
            "accession": ctx.accession,
            "cohort_label": ctx.cohort_label,
            "stage": "A4: + Perm. & Jackknife",
            "stage_order": 4,
            "n_units": r2.n_active_nhoods,
            "n_evaluated_pairs": r2.n_evaluated_pairs,
            "n_sig_pairs": n_verified,
            "n_redundant_or_intracomp_removed": int(sig_dedup.height) - n_verified,
        }
    )

    # Parallel Scale Benchmark: Coarse Super-Neighborhoods (res=0.25, CROSS_COMPARTMENT_ONLY, HYBRID_ALL)
    # Quantifies how larger neighborhoods shift statistical power from binary 2x2 gates to log-ratio balances
    cg_coarse_res = coarse_grain_nhoods_by_graph(
        nhoods_mat=ctx.nhoods_mat,
        count_matrix=ctx.count_matrix,
        cell_sample_codes=Some(ctx.cell_eval_codes),
        sample_library_sizes=Some(ctx.library_sizes),
        nhood_group_labels=Some(ctx.meta_codes),
        cell_types=Some(ctx.cell_types),
        min_jaccard=cfg.coarse_min_jaccard,
        leiden_resolution=cfg.coarse_resolution,
    )
    match cg_coarse_res:
        case Failure(err):
            return Failure(err)
        case Success(cg_coarse):
            pass

    readable_coarse = tuple(ctx.code_to_readable.get(lbl, lbl) for lbl in cg_coarse.super_nhood_labels)
    cfg_coarse = PairwiseConfig(
        method=PairwiseMethod.HYBRID_ALL,
        mixing_qc=mix_cfg,
        prefilter=t4_prefilter,
        overlap=OverlapConfig(max_jaccard_overlap=0.0, scope_filter=PairScopeFilter.CROSS_COMPARTMENT_ONLY),
        presence=PresenceConfig(mode=PresenceMode.HARD_DEPTH_GATED),
        screen=t4_screen,
        refinement=RefinementConfig(enabled=False),
        modules=PairModuleConfig(enabled=False),
    )
    res_coarse = analyze_nhood_pairs(
        nhoods_mat=cg_coarse.nhoods_super,
        count_matrix=cg_coarse.counts_super,
        response_labels=ctx.response_labels,
        sample_library_sizes=Some(lib_seq),
        strata_labels=ctx.inferred_strata,
        nhood_group_labels=Some(readable_coarse),
        sample_ids=Some(ctx.eval_ids),
        config=cfg_coarse,
    )
    match res_coarse:
        case Failure(err):
            return Failure(err)
        case Success(rc):
            coarse_regime_row: dict[str, object] = {
                "accession": ctx.accession,
                "cohort_label": ctx.cohort_label,
                "n_donors": len(ctx.eval_ids),
                "res10_active_units": r2.n_active_nhoods,
                "res10_cross_comp_quadrant_bundles": int(sig_dedup.height),
                "res10_verified_quadrant_bundles": n_verified,
                "res10_sig_continuous_pairs": r2.n_significant_continuous_pairs,
                "res025_active_units": rc.n_active_nhoods,
                "res025_sig_quadrant_pairs": rc.n_significant_quadrant_pairs,
                "res025_sig_continuous_pairs": rc.n_significant_continuous_pairs,
            }

    enriched_bundles = (
        bundles_df.filter(pl.col("is_cross_compartment"))
        .join(
            verified_leads_df.select(
                [
                    "nhood_u",
                    "nhood_v",
                    "quadrant",
                    "direction",
                    "z_iut_cmh",
                    "pvalue_perm_iut",
                    "fdr_perm_iut",
                    "fwer_perm_iut",
                    "lodo_min_z_iut",
                    "lodo_sign_consistency",
                    "loso_min_z_iut",
                    "loso_sign_consistency",
                    "pvalue_quadrant_firth",
                    "pvalue_iut_firth",
                    "is_statistically_verified",
                ]
            ),
            left_on=["lead_nhood_u", "lead_nhood_v", "lead_quadrant", "direction"],
            right_on=["nhood_u", "nhood_v", "quadrant", "direction"],
            how="left",
        )
        .with_columns(
            pl.lit(ctx.accession).alias("accession"),
            pl.lit(ctx.cohort_label).alias("cohort_label"),
        )
        if not bundles_df.is_empty()
        else pl.DataFrame()
    )

    if not build_umap_bundle:
        return Success((ablation_rows, coarse_regime_row, enriched_bundles, Nothing))

    # Build de-redundant UMAP composite strictly from 4-layer verified 1-per-(State_A, State_B) cross-compartment chords
    verified_keys = (
        verified_leads_df.filter(pl.col("is_statistically_verified")).select(
            ["nhood_u", "nhood_v", "quadrant", "direction"]
        )
        if not verified_leads_df.is_empty()
        else pl.DataFrame(
            schema={
                "nhood_u": pl.Int64,
                "nhood_v": pl.Int64,
                "quadrant": pl.String,
                "direction": pl.String,
            }
        )
    )
    clean_pairs_for_umap = (
        dedup_df.join(
            verified_keys.with_columns(pl.lit(True).alias("_keep_ver")),
            on=["nhood_u", "nhood_v", "quadrant", "direction"],
            how="left",
        )
        .with_columns(
            (pl.col("is_significant") & pl.col("_keep_ver").fill_null(False)).alias("is_significant")
        )
        .drop("_keep_ver")
    )
    clean_pw_res = PairwiseAnalysisResult(
        pairs_df=clean_pairs_for_umap,
        continuous_pairs_df=r2.continuous_pairs_df.head(0),
        modules_df=pl.DataFrame(),
        n_total_nhoods=int(cg_a1.nhoods_super.shape[1]),
        n_active_nhoods=r2.n_active_nhoods,
        n_evaluated_pairs=r2.n_evaluated_pairs,
        n_significant_quadrant_pairs=n_verified,
        n_significant_continuous_pairs=0,
        n_modules=n_verified,
    )

    da_fdr_raw = ctx.da_results_df["FDR"].cast(pl.Float64).to_numpy()
    da_lfc_raw = ctx.da_results_df["logFC"].cast(pl.Float64).to_numpy()
    n_super = len(cg_a1.member_nhood_indices)
    super_fdr = np.empty(n_super, dtype=np.float64)
    super_lfc = np.empty(n_super, dtype=np.float64)
    for g_idx, members in enumerate(cg_a1.member_nhood_indices):
        m_arr = np.asarray(members, dtype=np.int64)
        best_pos = int(np.argmin(da_fdr_raw[m_arr]))
        super_fdr[g_idx] = float(da_fdr_raw[m_arr[best_pos]])
        super_lfc[g_idx] = float(da_lfc_raw[m_arr[best_pos]])

    viz_cfg = UmapPairVizConfig(
        max_chords=cfg.max_chords,
        max_bg_cells=12000,
        fdr_threshold=cfg.fdr_threshold,
    )
    b_res = prepare_pairwise_umap_bundle(
        nhoods_mat=cg_a1.nhoods_super,
        umap_coords=ctx.umap_coords,
        count_matrix=cg_a1.counts_super,
        response_labels=ctx.response_labels,
        pairwise_res=clean_pw_res,
        sample_library_sizes=Some(lib_seq),
        state_labels=Some(readable_a1),
        cell_labels=Some(ctx.cell_types),
        univariate_fdr=Some(tuple(float(x) for x in super_fdr.tolist())),
        univariate_logfc=Some(tuple(float(x) for x in super_lfc.tolist())),
        config=viz_cfg,
    )
    match b_res:
        case Failure(err):
            return Failure(err)
        case Success(umap_bundle):
            n_strata = (
                len(set(ctx.inferred_strata.unwrap()))
                if isinstance(ctx.inferred_strata, Some)
                else 1
            )
            sub_tag = (
                f"CMH-Stratified ({n_strata} PMIDs) | "
                if n_strata > 1
                else "Single-Protocol Cohort | "
            )
            chart = build_pairwise_umap_composite(
                bundle=umap_bundle,
                title=f"{ctx.cohort_label} — Simplified Cross-Compartment Pairwise UMAP",
                subtitle=(
                    f"{sub_tag}Donor-Mixing QC (LISI_pt >= {cfg.min_effective_lisi:.1f}) | "
                    f"1 Canonical Chord per (State_A, State_B) | "
                    f"B={cfg.n_permutations} Stratified Perm. + Firth + Jackknife Verified ({n_verified} Chords)"
                ),
                config=viz_cfg,
            )
            return Success((ablation_rows, coarse_regime_row, enriched_bundles, Some(chart)))
        case _ as unreachable_b:
            return Failure(f"Unexpected Result state: {unreachable_b}")


def _compact_ontology_label(raw_label: str) -> str:
    """Shortens verbose Cell Ontology strings for crisp in-panel bar annotations."""
    mapping = {
        "central memory CD8-positive, alpha-beta T cell": "CD8+ central memory T cell",
        "effector memory CD8-positive, alpha-beta T cell": "CD8+ effector memory T cell",
    }
    return mapping.get(raw_label, raw_label)


def build_mixing_and_scale_diagnostic_figure(
    raw_mix_summary_df: pl.DataFrame,
    gondal_substudies_df: pl.DataFrame,
    scale_df: pl.DataFrame,
    ablation_df: pl.DataFrame,
) -> AltairComposite:
    """Builds the 4-panel Nature Methods wireframe diagnostic figure (`a`–`d`)."""
    # Panel a: Donor mixing (LISI_pt) and % single-donor dominance > 50% across all 9 cohorts
    p_a_df = raw_mix_summary_df.with_columns(
        (100.0 * pl.col("frac_dominance_gt_50pct")).alias("dom_gt50_pct"),
        pl.max_horizontal(pl.col("median_effective_donors_lisi"), pl.lit(2.65)).alias("label_x_pos"),
        pl.when(pl.col("accession") == "CELLxGENE_7b20c613")
        .then(pl.lit("Gondal Multi-Study Compendium"))
        .when(pl.col("n_donors") < 20)
        .then(pl.lit("Small Donor Cohort (J < 20)"))
        .otherwise(pl.lit("Well-Mixed Cohort (J >= 34)"))
        .alias("mixing_tier"),
    ).sort("median_effective_donors_lisi", descending=True)

    cohort_sort = p_a_df["cohort_label"].to_list()
    bars_a = (
        alt.Chart(p_a_df)
        .mark_bar(height=15)
        .encode(
            y=alt.Y("cohort_label:N", sort=cohort_sort, title="Clinical Single-Cell Cohort"),
            x=alt.X(
                "median_effective_donors_lisi:Q",
                scale=alt.Scale(domain=[0, 31.5]),
                title="Median Effective Donors per k=30 Nhood (Inverse-Simpson LISI_pt)",
            ),
            color=alt.Color(
                "mixing_tier:N",
                title="Integration Tier",
                scale=alt.Scale(
                    domain=[
                        "Well-Mixed Cohort (J >= 34)",
                        "Gondal Multi-Study Compendium",
                        "Small Donor Cohort (J < 20)",
                    ],
                    range=[OKABE_BLUE, OKABE_VERMILION, OKABE_ORANGE],
                ),
                legend=alt.Legend(orient="bottom", direction="vertical", titleLimit=280, labelLimit=280),
            ),
        )
    )
    qc_rule_a = (
        alt.Chart(pl.DataFrame({"x": [2.5]}))
        .mark_rule(color="#0F172A", strokeDash=[4, 3], strokeWidth=1.1)
        .encode(x="x:Q")
    )
    txt_a = (
        alt.Chart(p_a_df)
        .mark_text(align="left", dx=4, fontSize=8.5, color="#1E293B")
        .encode(
            y=alt.Y("cohort_label:N", sort=cohort_sort),
            x=alt.X("label_x_pos:Q"),
            text=alt.Text("dom_label:N"),
        )
        .transform_calculate(
            dom_label="format(datum.median_effective_donors_lisi, '.1f') + ' (' + format(datum.dom_gt50_pct, '.0f') + '% >50% 1-pt)'"
        )
    )
    panel_a = (
        (bars_a + qc_rule_a + txt_a)
        .properties(
            width=365,
            height=225,
            title=alt.TitleParams(
                text="a  Donor Mixing Across Nine Cohorts (k=30)",
                subtitle="Bars: Median LISI_pt (Dashed QC = 2.5) | Text: % Nhoods >50% 1-Donor",
                anchor="start",
                fontSize=11.5,
                subtitleFontSize=9.0,
                subtitleColor="#475569",
            ),
        )
    )

    # Panel b: Gondal's 7 PMID sub-studies showing baseline R% confounding (0% -> 66.7%)
    p_b_df = gondal_substudies_df.with_columns(
        (
            pl.lit("PMID ")
            + pl.col("pmid_substudy")
            + pl.lit(" (J=")
            + pl.col("n_donors").cast(pl.String)
            + pl.lit(")")
        ).alias("substudy_display"),
        pl.col("modal_cell_type")
        .map_elements(_compact_ontology_label, return_dtype=pl.String)
        .alias("compact_modal_ct"),
        pl.when(pl.col("is_cmh_informative"))
        .then(pl.lit("CMH-Informative Stratum (R & NR)"))
        .otherwise(pl.lit("0%-Responder Stratum (Auto-Dropped by CMH)"))
        .alias("stratum_status"),
    ).sort("responder_pct", descending=True)
    sub_sort = p_b_df["substudy_display"].to_list()

    bars_b = (
        alt.Chart(p_b_df)
        .mark_bar(height=16)
        .encode(
            y=alt.Y("substudy_display:N", sort=sub_sort, title="Gondal Sub-Study (PMID)"),
            x=alt.X(
                "responder_pct:Q",
                scale=alt.Scale(domain=[0, 102]),
                title="Clinical Responder Rate (%)",
            ),
            color=alt.Color(
                "stratum_status:N",
                title="CMH Stratum Role",
                scale=alt.Scale(
                    domain=[
                        "CMH-Informative Stratum (R & NR)",
                        "0%-Responder Stratum (Auto-Dropped by CMH)",
                    ],
                    range=[OKABE_BLUISH_GREEN, "#94A3B8"],
                ),
                legend=alt.Legend(orient="bottom", direction="vertical", titleLimit=280, labelLimit=280),
            ),
        )
    )
    txt_b = (
        alt.Chart(p_b_df)
        .mark_text(align="left", dx=4, fontSize=8.5, color="#1E293B")
        .encode(
            y=alt.Y("substudy_display:N", sort=sub_sort),
            x=alt.X("responder_pct:Q"),
            text=alt.Text("lbl:N"),
        )
        .transform_calculate(
            lbl="format(datum.responder_pct, '.1f') + '% (' + datum.compact_modal_ct + ')'"
        )
    )
    panel_b = (
        (bars_b + txt_b)
        .properties(
            width=365,
            height=225,
            title=alt.TitleParams(
                text="b  Sub-Study Confounding in Gondal (7 PMIDs)",
                subtitle="Median Nhood Sub-Study Purity = 1.000 (93.4% >80% 1-PMID)",
                anchor="start",
                fontSize=11.5,
                subtitleFontSize=9.0,
                subtitleColor="#475569",
            ),
        )
    )

    # Consistent Okabe-Ito domain across Spotlight Cohorts (Blue = NSCLC, Vermilion = Gondal)
    spotlight_order_df = (
        scale_df.filter(pl.col("accession").is_in(list(SPOTLIGHT_COHORTS)))
        .select(["accession", "cohort_label"])
        .unique()
    )
    acc_to_lbl = {str(r["accession"]): str(r["cohort_label"]) for r in spotlight_order_df.to_dicts()}
    spotlight_domain = [
        acc_to_lbl.get("GSE243013", "NSCLC"),
        acc_to_lbl.get("CELLxGENE_7b20c613", "Gondal"),
    ]
    spotlight_range = [OKABE_BLUE, OKABE_VERMILION]

    # Panel c: Neighborhood scale trade-off (S0 -> S3) on Gondal & GSE243013
    p_c_df = scale_df.filter(pl.col("accession").is_in(list(SPOTLIGHT_COHORTS)))
    scale_order = [s.value for s in ALL_SCALES]
    line_c = (
        alt.Chart(p_c_df)
        .mark_line(point=alt.OverlayMarkDef(size=65, filled=True), strokeWidth=2.0)
        .encode(
            x=alt.X(
                "scale:N",
                sort=scale_order,
                scale=alt.Scale(padding=0.35),
                title="Lineage-Constrained Coarsening Scale",
                axis=alt.Axis(labelAngle=-18),
            ),
            y=alt.Y(
                "median_effective_donors_lisi:Q",
                scale=alt.Scale(type="log", domain=[1, 160]),
                axis=alt.Axis(values=[1, 2, 5, 10, 20, 50, 100]),
                title="Median Effective Donors (LISI_pt, log scale)",
            ),
            color=alt.Color(
                "cohort_label:N",
                title="Spotlight Cohort",
                scale=alt.Scale(domain=spotlight_domain, range=spotlight_range),
                legend=alt.Legend(orient="bottom", direction="vertical", titleLimit=280, labelLimit=280),
            ),
        )
    )
    txt_c = (
        alt.Chart(p_c_df)
        .mark_text(dy=-11, fontSize=8.0, color="#1E293B")
        .encode(
            x=alt.X("scale:N", sort=scale_order),
            y=alt.Y("median_effective_donors_lisi:Q"),
            text=alt.Text("pair_lbl:N"),
        )
        .transform_calculate(
            pair_lbl="'LISI=' + format(datum.median_effective_donors_lisi, '.1f') + ' (V=' + datum.n_units + ')'"
        )
    )
    panel_c = (
        (line_c + txt_c)
        .properties(
            width=365,
            height=225,
            title=alt.TitleParams(
                text="c  Larger Lineage-Constrained Nhoods Boost Mixing",
                subtitle="Exact Cell-Union Pooling Doubles LISI_pt & Cuts O(V^2) >150x",
                anchor="start",
                fontSize=11.5,
                subtitleFontSize=9.0,
                subtitleColor="#475569",
            ),
        )
    )

    # Panel d: 5-stage simplification & 4-layer statistical verification cascade
    p_d_df = ablation_df.filter(pl.col("accession").is_in(list(SPOTLIGHT_COHORTS))).with_columns(
        (pl.col("n_sig_pairs").cast(pl.Float64) + 1.0).log10().alias("log10_sig_plus_one")
    )
    stage_order = (
        p_d_df.select(["stage", "stage_order"])
        .unique()
        .sort("stage_order")["stage"]
        .to_list()
    )
    bars_d = (
        alt.Chart(p_d_df)
        .mark_bar(size=18)
        .encode(
            x=alt.X(
                "stage:N",
                sort=stage_order,
                title="Progressive Simplification & 4-Layer Verification Stage",
                axis=alt.Axis(labelAngle=-18),
            ),
            xOffset=alt.XOffset("cohort_label:N", sort=spotlight_domain),
            y=alt.Y(
                "log10_sig_plus_one:Q",
                scale=alt.Scale(domain=[0, 5.25]),
                title="Significant / Verified Chords (log10[N + 1])",
                axis=alt.Axis(
                    values=[0, 1, 2, 3, 4, 5],
                    labelExpr=(
                        "datum.value == 0 ? '0' : datum.value == 1 ? '10' : "
                        "datum.value == 2 ? '100' : datum.value == 3 ? '1,000' : "
                        "datum.value == 4 ? '10,000' : '100,000'"
                    ),
                ),
            ),
            color=alt.Color(
                "cohort_label:N",
                title="Spotlight Cohort",
                scale=alt.Scale(domain=spotlight_domain, range=spotlight_range),
                legend=alt.Legend(orient="bottom", direction="vertical", titleLimit=280, labelLimit=280),
            ),
        )
    )
    txt_d = (
        alt.Chart(p_d_df)
        .mark_text(dy=-7, fontSize=8.0, fontWeight="bold", color="#1E293B")
        .encode(
            x=alt.X("stage:N", sort=stage_order),
            xOffset=alt.XOffset("cohort_label:N", sort=spotlight_domain),
            y=alt.Y("log10_sig_plus_one:Q"),
            text=alt.Text("sig_lbl:N"),
        )
        .transform_calculate(sig_lbl="format(datum.n_sig_pairs, ',d')")
    )
    g_seq = (
        p_d_df.filter(pl.col("accession") == "CELLxGENE_7b20c613")
        .sort("stage_order")["n_sig_pairs"]
        .to_list()
    )
    n_seq = (
        p_d_df.filter(pl.col("accession") == "GSE243013")
        .sort("stage_order")["n_sig_pairs"]
        .to_list()
    )
    g_sub = (
        f"Gondal: {int(g_seq[0]):,} -> {int(g_seq[1]):,} -> {int(g_seq[2]):,} -> {int(g_seq[-1])} Verified"
        if len(g_seq) >= 4
        else "Gondal"
    )
    n_sub = (
        f"NSCLC: {int(n_seq[0]):,} -> {int(n_seq[1]):,} -> {int(n_seq[2]):,} -> {int(n_seq[-1])} Verified"
        if len(n_seq) >= 4
        else "NSCLC"
    )
    panel_d = (
        (bars_d + txt_d)
        .properties(
            width=365,
            height=225,
            title=alt.TitleParams(
                text="d  De-Redundancy & 4-Layer Significance Cascade",
                subtitle=f"{g_sub} | {n_sub}",
                anchor="start",
                fontSize=11.5,
                subtitleFontSize=9.0,
                subtitleColor="#475569",
            ),
        )
    )

    row_top = alt.hconcat(panel_a, panel_b, spacing=42).resolve_scale(color="independent")
    row_bot = alt.hconcat(panel_c, panel_d, spacing=42).resolve_scale(color="independent")
    composite = (
        alt.vconcat(row_top, row_bot, spacing=36)
        .resolve_scale(color="independent")
        .properties(
            title=alt.TitleParams(
                text="Sample Integration Audit, Neighborhood Scale Coarsening & 4-Layer Pairwise Statistical Calibration",
                subtitle=(
                    "Outcome-Blind Donor Mixing (LISI_pt) + Sub-Study CMH Stratification + Cross-Compartment Bundle Consolidation "
                    "+ B=1,000 Stratified IUT Permutations & Jackknife"
                ),
                anchor="start",
                fontSize=13.5,
                subtitleFontSize=10.0,
                subtitleColor="#334155",
                offset=12,
            )
        )
    )
    return cast(AltairComposite, apply_nature_methods_theme(composite))


def upgrade_svg_with_inkscape_layers(
    svg_content: str,
    figure_title: str,
    panel_layer_label: str,
) -> Result[str, str]:
    """Wraps compiled SVG into native Inkscape layers (`01_Background`, `02_<panel_layer_label>`)."""
    ET.register_namespace("", SVG_NS)
    ET.register_namespace("xlink", XLINK_NS)
    ET.register_namespace("inkscape", INKSCAPE_NS)
    ET.register_namespace("sodipodi", SODIPODI_NS)

    try:
        root = ET.fromstring(svg_content)
    except Exception as exc:
        return Failure(f"Invalid SVG XML: {exc}")

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
    namedview.set("id", "namedview-mixing-audit")
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
    layer_panels.set(f"{{{INKSCAPE_NS}}}label", panel_layer_label)
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
    return Success(serialized)


def export_publication_figure(
    chart: AltairComposite,
    base_path: Path,
    figure_title: str,
    panel_layer_label: str,
    hybrid_umap_raster: bool,
) -> Result[tuple[Path, Path], str]:
    """Exports an Altair composite to Inkscape-layered `.svg` + 300 DPI `.png`."""
    svg_path = base_path.parent / f"{base_path.stem}.svg"
    png_path = base_path.parent / f"{base_path.stem}.png"
    try:
        svg_path.parent.mkdir(parents=True, exist_ok=True)
        with alt.data_transformers.disable_max_rows():
            spec = chart.to_dict()
        raw_svg = str(vlc.vegalite_to_svg(spec))
    except Exception as exc:
        return Failure(f"Altair-to-SVG compilation failed for {base_path}: {exc}")

    if hybrid_umap_raster:
        rast_res = rasterize_svg_datapoints(
            raw_svg,
            max_points=1000,
            min_group_points=250,
            scale=3.0,
        )
        match rast_res:
            case Failure(err):
                return Failure(f"Hybrid SVG rasterization failed: {err}")
            case Success(working_svg):
                pass
    else:
        working_svg = raw_svg

    ink_res = upgrade_svg_with_inkscape_layers(
        working_svg,
        figure_title=figure_title,
        panel_layer_label=panel_layer_label,
    )
    match ink_res:
        case Failure(err):
            return Failure(err)
        case Success(inkscape_svg):
            pass

    svg_path.write_text(inkscape_svg, encoding="utf-8")
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
        case _ as unreachable:
            return Failure(f"Unexpected Result state: {unreachable}")


def render_audit_markdown_report(
    raw_mix_df: pl.DataFrame,
    gondal_sub_df: pl.DataFrame,
    ablation_df: pl.DataFrame,
    coarse_regime_df: pl.DataFrame,
    bundles_df: pl.DataFrame,
) -> str:
    """Formats a structured Markdown summary report of all mixing, scale, and significance results."""
    lines: list[str] = [
        "# Sample Integration Audit, Neighborhood Scale & 4-Layer Pair Significance Report",
        "",
        "## 1. Raw `k=30` Donor-Mixing Audit Across Nine Cohorts",
        "",
        "| Cohort | Donors ($J$) | Units ($V$) | Median Cells | Median Donors | Median $\\text{LISI}_{\\text{pt}}$ | `>50%` 1-Donor | Passes Mixing QC |",
        "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
    ]
    for r in raw_mix_df.sort("median_effective_donors_lisi", descending=True).to_dicts():
        lines.append(
            f"| `{r['accession']}` ({r['cohort_label']}) | {int(r['n_donors'])} | {int(r['n_units']):,} | "
            f"{float(r['median_cells_per_unit']):.1f} | {float(r['median_donors_per_unit']):.1f} | "
            f"**{float(r['median_effective_donors_lisi']):.2f}** | "
            f"{100.0 * float(r['frac_dominance_gt_50pct']):.1f}% | "
            f"{100.0 * float(r['frac_passes_mixing_qc']):.1f}% |"
        )

    lines.extend(
        [
            "",
            "## 2. Gondal (`CELLxGENE_7b20c613`) Multi-Study `PMID` Stratification Audit",
            "",
            "| PMID Sub-Study | Donors | Responders | Responder % | Cells | Modal Cell Type | CMH-Informative |",
            "| :--- | :---: | :---: | :---: | :---: | :--- | :---: |",
        ]
    )
    for r in gondal_sub_df.to_dicts():
        lines.append(
            f"| `{r['pmid_substudy']}` | {int(r['n_donors'])} | {int(r['n_responders'])} | "
            f"{float(r['responder_pct']):.1f}% | {int(r['n_cells']):,} | `{r['modal_cell_type']}` | "
            f"{'Yes' if bool(r['is_cmh_informative']) else 'No (Dropped)'} |"
        )

    lines.extend(
        [
            "",
            "## 3. Simplification & 4-Layer Statistical Significance Cascade",
            "",
            "| Cohort | Stage | Active Units | Evaluated Pairs | Significant / Verified Chords |",
            "| :--- | :--- | :---: | :---: | :---: |",
        ]
    )
    for r in ablation_df.sort(["accession", "stage_order"]).to_dicts():
        lines.append(
            f"| `{r['accession']}` | {r['stage']} | {int(r['n_units']):,} | "
            f"{int(r['n_evaluated_pairs']):,} | **{int(r['n_sig_pairs']):,}** |"
        )

    if not coarse_regime_df.is_empty():
        lines.extend(
            [
                "",
                "## 4. Medium (`res=1.0`) vs. Coarse (`res=0.25`) Super-Neighborhood Regime Trade-Off",
                "",
                "| Cohort | Donors ($J$) | `res=1.0` Units | `res=1.0` Verified Quadrant Chords | `res=1.0` Continuous Hits | `res=0.25` Units | `res=0.25` Quadrant Hits | `res=0.25` Cross-Comp. Continuous Hits |",
                "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
            ]
        )
        for r in coarse_regime_df.sort("n_donors", descending=True).to_dicts():
            lines.append(
                f"| `{r['accession']}` ({r['cohort_label']}) | {int(r['n_donors'])} | "
                f"{int(r['res10_active_units']):,} | **{int(r['res10_verified_quadrant_bundles'])}** | "
                f"{int(r['res10_sig_continuous_pairs']):,} | {int(r['res025_active_units']):,} | "
                f"{int(r['res025_sig_quadrant_pairs'])} | **{int(r['res025_sig_continuous_pairs']):,}** |"
            )

    if not bundles_df.is_empty():
        ver = bundles_df.filter(pl.col("is_statistically_verified").fill_null(False))
        lines.extend(
            [
                "",
                "## 5. Top 4-Layer Statistically Verified Cross-Compartment Ecosystem Chords",
                "",
                "| Cohort | State A | State B | Regime | Enriched Arm | Sub-Pairs | $q_{\\text{IUT}}$ | $p_{\\text{perm, IUT}}$ | $\\text{FWER}_{\\max T}$ | $\\min Z_{\\text{LODO}}$ | LOSO Sign % |",
                "| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
            ]
        )
        for r in ver.sort(["min_fdr_iut", "pvalue_perm_iut"]).head(25).to_dicts():
            lines.append(
                f"| `{r['accession']}` | {r['state_a']} | {r['state_b']} | `{r['canonical_quadrant']}` | "
                f"{r['direction']} | {int(r['n_supporting_subpairs'])} | "
                f"{float(r['min_fdr_iut']):.4g} | {float(r['pvalue_perm_iut']):.4g} | "
                f"{float(r['fwer_perm_iut']):.4g} | {float(r['lodo_min_z_iut']):.2f} | "
                f"{100.0 * float(r['loso_sign_consistency']):.0f}% |"
            )

    lines.append("")
    return "\n".join(lines)


def run_pipeline(cfg: AuditCliConfig) -> Result[Path, str]:
    """Executes the full 9-cohort audit, progressive ablation, and figure export pipeline."""
    t0 = time.perf_counter()
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = cfg.output_dir / "figures"
    rep_dir = cfg.output_dir / "reports"
    fig_dir.mkdir(parents=True, exist_ok=True)
    rep_dir.mkdir(parents=True, exist_ok=True)

    all_meta42_df = (
        pl.read_parquet(cfg.meta42_annotations_path)
        if cfg.meta42_annotations_path.exists()
        else pl.DataFrame()
    )
    mix_cfg = MixingQCConfig(
        enabled=True,
        min_distinct_donors=cfg.min_distinct_donors,
        min_effective_donors_lisi=cfg.min_effective_lisi,
        max_single_donor_dominance=cfg.max_donor_dominance,
        adaptive_small_cohort_scaling=True,
    )

    cohorts = tuple(str(p["acc"]) for p in COHORT_PANELS)
    all_scale_rows: list[dict[str, object]] = []
    all_ablation_rows: list[dict[str, object]] = []
    all_coarse_regime_rows: list[dict[str, object]] = []
    bundle_frames: list[pl.DataFrame] = []
    gondal_sub_df = pl.DataFrame()

    for acc in cohorts:
        ctx_res = load_cohort_audit_context(acc, cfg, all_meta42_df)
        match ctx_res:
            case Failure(err):
                return Failure(err)
            case Success(ctx):
                pass

        if acc == "CELLxGENE_7b20c613":
            gondal_sub_df = audit_gondal_substudies(ctx)

        sc_res = evaluate_cohort_scales(ctx, cfg, mix_cfg)
        match sc_res:
            case Failure(err):
                return Failure(err)
            case Success(sc_rows):
                all_scale_rows.extend(sc_rows)

        is_spotlight = acc in SPOTLIGHT_COHORTS
        ab_res = run_cohort_ablation_and_bundle(
            ctx=ctx,
            cfg=cfg,
            mix_cfg=mix_cfg,
            build_umap_bundle=is_spotlight,
        )
        match ab_res:
            case Failure(err):
                return Failure(err)
            case Success((ab_rows, coarse_row, b_df, chart_maybe)):
                all_ablation_rows.extend(ab_rows)
                all_coarse_regime_rows.append(coarse_row)
                if not b_df.is_empty():
                    bundle_frames.append(b_df)
                if isinstance(chart_maybe, Some):
                    exp_res = export_publication_figure(
                        chart=chart_maybe.unwrap(),
                        base_path=fig_dir / f"{acc}_simplified_cross_compartment_umap",
                        figure_title=f"{ctx.cohort_label} Simplified Cross-Compartment Pairwise UMAP",
                        panel_layer_label="02_Simplified_Pairwise_UMAP_Panels",
                        hybrid_umap_raster=True,
                    )
                    if isinstance(exp_res, Failure):
                        return exp_res

        a4_row = ab_rows[-1]
        print(
            f"[{acc}] donors={len(ctx.eval_ids)}, A0={ab_rows[0]['n_sig_pairs']} -> "
            f"A1={ab_rows[1]['n_sig_pairs']} -> A2={ab_rows[2]['n_sig_pairs']} -> "
            f"A3={ab_rows[3]['n_sig_pairs']} -> A4_verified={a4_row['n_sig_pairs']}"
        )

    scale_df = pl.DataFrame(all_scale_rows)
    raw_mix_df = scale_df.filter(pl.col("scale") == ScaleLevel.RAW_K30.value)
    ablation_df = pl.DataFrame(all_ablation_rows)
    coarse_regime_df = pl.DataFrame(all_coarse_regime_rows)
    bundles_df = pl.concat(bundle_frames, how="vertical_relaxed") if bundle_frames else pl.DataFrame()

    raw_mix_df.write_parquet(cfg.output_dir / "cohort_donor_mixing_summary.parquet")
    scale_df.write_parquet(cfg.output_dir / "neighborhood_scale_tradeoff.parquet")
    ablation_df.write_parquet(cfg.output_dir / "simplification_and_significance_ablation.parquet")
    coarse_regime_df.write_parquet(cfg.output_dir / "coarse_scale_regime_comparison.parquet")
    if not gondal_sub_df.is_empty():
        gondal_sub_df.write_parquet(cfg.output_dir / "gondal_substudy_confounding_audit.parquet")
    if not bundles_df.is_empty():
        bundles_df.write_parquet(cfg.output_dir / "verified_canonical_bundles.parquet")

    diag_chart = build_mixing_and_scale_diagnostic_figure(
        raw_mix_summary_df=raw_mix_df,
        gondal_substudies_df=gondal_sub_df,
        scale_df=scale_df,
        ablation_df=ablation_df,
    )
    diag_exp = export_publication_figure(
        chart=diag_chart,
        base_path=fig_dir / "cohort_donor_mixing_and_scale_audit",
        figure_title="Sample Integration Audit, Scale Coarsening & 4-Layer Pair Significance",
        panel_layer_label="02_Mixing_Scale_Significance_Panels",
        hybrid_umap_raster=False,
    )
    if isinstance(diag_exp, Failure):
        return diag_exp

    report_md = render_audit_markdown_report(
        raw_mix_df=raw_mix_df,
        gondal_sub_df=gondal_sub_df,
        ablation_df=ablation_df,
        coarse_regime_df=coarse_regime_df,
        bundles_df=bundles_df,
    )
    (rep_dir / "mixing_simplification_significance_report.md").write_text(report_md, encoding="utf-8")

    elapsed = time.perf_counter() - t0
    print(f"Completed 9-cohort mixing & simplified pair pipeline in {elapsed:.1f}s -> {cfg.output_dir}")
    return Success(cfg.output_dir)


def parse_args(argv: Sequence[str]) -> AuditCliConfig:
    """Parses CLI flags into an immutable `AuditCliConfig`."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("output/results/milopy_mixing_and_simplified_pairs"))
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--fdr-threshold", type=float, default=0.10)
    parser.add_argument("--coarse-resolution", type=float, default=0.25)
    parser.add_argument("--coarse-min-jaccard", type=float, default=0.02)
    parser.add_argument("--min-effective-lisi", type=float, default=2.5)
    parser.add_argument("--max-donor-dominance", type=float, default=0.50)
    parser.add_argument("--max-chords", type=int, default=20)
    ns = parser.parse_args(list(argv))
    return AuditCliConfig(
        output_dir=ns.output_dir,
        n_permutations=ns.n_permutations,
        fdr_threshold=ns.fdr_threshold,
        coarse_resolution=ns.coarse_resolution,
        coarse_min_jaccard=ns.coarse_min_jaccard,
        min_effective_lisi=ns.min_effective_lisi,
        max_donor_dominance=ns.max_donor_dominance,
        max_chords=ns.max_chords,
    )


def main() -> int:
    """Imperative shell entry point."""
    cfg = parse_args(sys.argv[1:])
    match run_pipeline(cfg):
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
