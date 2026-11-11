"""Pure functional monadic pipeline for pairwise neighborhood differential association.

Orchestrates overlap topology construction, depth-calibrated presence estimation,
vectorized 4-quadrant and continuous screening, Firth penalized LRT refinement,
and product-graph ecosystem module grouping without mutating input AnnData objects.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import assert_never

import anndata as ad
import numpy as np
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import sparse as sp

from .grouping import _empty_modules_df, group_nhood_pairs
from .mixing_qc import compute_nhood_mixing_metrics
from .overlap import compute_nhood_overlap_graph
from .power import build_compartment_observability_mask, depth_residualize_continuous_features
from .prefilter import (
    coarse_grain_nhoods_by_graph,
    select_informative_nhoods_outcome_blind,
    subset_presence_representation,
)
from .presence import estimate_nhood_presence
from .refine import refine_candidate_pairs_firth
from .screen import (
    _empty_continuous_df,
    _empty_quadrant_df,
    screen_continuous_pairs,
    screen_quadrant_pairs,
)
from .types import (
    CandidatePairRefiner,
    ContinuousAxis,
    FloatArray,
    IntArray,
    OverlapGraph,
    PairwiseAnalysisResult,
    PairwiseConfig,
    PairwiseMethod,
    PrefilterMode,
    PresenceEstimator,
    PresenceRepresentation,
    ScreenConfig,
    SparseMatrix,
)


def _encode_response_labels(response_labels: Sequence[str], n_samples: int) -> Result[FloatArray, str]:
    """Validates and maps binary clinical response labels to 1.0 (responder) / 0.0 (non-responder)."""
    if len(response_labels) != n_samples:
        return Failure(
            f"response_labels length ({len(response_labels)}) does not match sample count ({n_samples})"
        )
    normalized = tuple(str(lbl).strip().lower() for lbl in response_labels)
    invalid = sorted({lbl for lbl in normalized if lbl not in ("responder", "non-responder")})
    if invalid:
        return Failure(
            f"response_labels must contain only 'responder' and 'non-responder', found {invalid}"
        )

    y_bin = np.array([1.0 if lbl == "responder" else 0.0 for lbl in normalized], dtype=np.float64)
    n_resp = int(np.sum(y_bin > 0.5))
    n_non_resp = n_samples - n_resp
    if n_resp < 2 or n_non_resp < 2:
        return Failure(
            f"Pairwise analysis requires >= 2 responders and >= 2 non-responders, got {n_resp} R and {n_non_resp} NR"
        )
    return Success(y_bin)


def _encode_strata_labels(
    strata_labels: Maybe[Sequence[str]],
    n_samples: int,
) -> Result[Maybe[IntArray], str]:
    """Encodes optional cohort/study stratum strings into integer stratum codes."""
    match strata_labels:
        case Maybe.empty:
            return Success(Nothing)
        case Some(labels):
            if len(labels) != n_samples:
                return Failure(
                    f"strata_labels length ({len(labels)}) does not match sample count ({n_samples})"
                )
            norm_s = tuple(str(s) for s in labels)
            unique_s = sorted(set(norm_s))
            if len(unique_s) <= 1:
                return Success(Nothing)
            s_map = {name: idx for idx, name in enumerate(unique_s)}
            codes = np.array([s_map[name] for name in norm_s], dtype=np.int64)
            return Success(Some(codes))
        case _ as unreachable:
            return Failure(f"Unexpected Maybe state: {unreachable}")


def _pool_prior_weights_to_super(
    prior_weights: FloatArray,
    member_nhood_indices: tuple[tuple[int, ...], ...],
) -> FloatArray:
    """Symmetrically max-pools `(V, V)` prior weights onto `(G_super, G_super)` super-neighborhoods."""
    n_super = len(member_nhood_indices)
    out = np.zeros((n_super, n_super), dtype=np.float64)
    for g1 in range(n_super):
        idx1 = np.asarray(member_nhood_indices[g1], dtype=np.int64)
        sub_rows = prior_weights[idx1, :]
        for g2 in range(g1 + 1, n_super):
            idx2 = np.asarray(member_nhood_indices[g2], dtype=np.int64)
            val = float(np.max(sub_rows[:, idx2]))
            out[g1, g2] = val
            out[g2, g1] = val
    return out


def _resolve_effective_screen_config(config: PairwiseConfig) -> ScreenConfig:
    """Resolves whether exact outcome-blind Tarone IUT pruning is active."""
    use_tarone = (
        config.screen.apply_tarone_iut_pruning
        or config.prefilter.apply_tarone_iut_pruning
        or config.prefilter.mode in (PrefilterMode.TARONE_ONLY, PrefilterMode.COMPOSITE_SUPER_TARONE)
    )
    if not use_tarone:
        return config.screen
    p_ceil = (
        config.prefilter.tarone_p_ceiling
        if (
            config.prefilter.apply_tarone_iut_pruning
            or config.prefilter.mode in (PrefilterMode.TARONE_ONLY, PrefilterMode.COMPOSITE_SUPER_TARONE)
        )
        else config.screen.tarone_p_ceiling
    )
    return ScreenConfig(
        chunk_size=config.screen.chunk_size,
        fdr_threshold=config.screen.fdr_threshold,
        min_donors_per_quadrant=config.screen.min_donors_per_quadrant,
        min_effective_samples=config.screen.min_effective_samples,
        pseudocount=config.screen.pseudocount,
        use_spatial_fdr_weights=config.screen.use_spatial_fdr_weights,
        n_permutations=config.screen.n_permutations,
        permutation_seed=config.screen.permutation_seed,
        max_retained_candidates=config.screen.max_retained_candidates,
        apply_tarone_iut_pruning=True,
        tarone_p_ceiling=p_ceil,
        dmpp_null_draws=config.screen.dmpp_null_draws,
    )


def analyze_nhood_pairs(
    nhoods_mat: SparseMatrix,
    count_matrix: FloatArray,
    response_labels: Sequence[str],
    sample_library_sizes: Maybe[Sequence[float]] = Nothing,
    strata_labels: Maybe[Sequence[str]] = Nothing,
    cell_types: Maybe[Sequence[str]] = Nothing,
    nhood_group_labels: Maybe[Sequence[str]] = Nothing,
    cell_sample_codes: Maybe[IntArray] = Nothing,
    prior_pair_weights: Maybe[FloatArray] = Nothing,
    sample_ids: Maybe[Sequence[str]] = Nothing,
    cohort_labels: Maybe[Sequence[str]] = Nothing,
    config: PairwiseConfig = PairwiseConfig(),
    presence_estimator: Maybe[PresenceEstimator] = Nothing,
    refiner: Maybe[CandidatePairRefiner] = Nothing,
) -> Result[PairwiseAnalysisResult, str]:
    """Pure array/matrix pipeline discovering pairwise neighborhood associations with response."""
    counts_arr = np.asarray(count_matrix, dtype=np.float64)
    if counts_arr.ndim != 2:
        return Failure(f"count_matrix must be 2D, got shape {counts_arr.shape}")

    n_nhoods_orig, n_samples = int(counts_arr.shape[0]), int(counts_arr.shape[1])
    if int(nhoods_mat.shape[1]) != n_nhoods_orig:
        return Failure(
            f"nhoods_mat columns ({nhoods_mat.shape[1]}) must match count_matrix rows ({n_nhoods_orig})"
        )

    y_res = _encode_response_labels(response_labels, n_samples)
    match y_res:
        case Failure(err):
            return Failure(err)
        case Success(response_binary):
            pass

    strata_res = _encode_strata_labels(strata_labels, n_samples)
    match strata_res:
        case Failure(err):
            return Failure(err)
        case Success(strata_codes):
            pass

    lib_maybe: Maybe[FloatArray] = (
        Some(np.asarray(sample_library_sizes.unwrap(), dtype=np.float64))
        if isinstance(sample_library_sizes, Some)
        else Nothing
    )

    use_super = config.prefilter.mode in (
        PrefilterMode.GRAPH_SUPER_NHOODS,
        PrefilterMode.COMPOSITE_SUPER_TARONE,
    )
    use_marginal_filter = config.prefilter.mode in (
        PrefilterMode.MARGINAL_ENTROPY_VARIANCE,
        PrefilterMode.COMPOSITE_SUPER_TARONE,
    )

    work_nhoods: SparseMatrix = nhoods_mat
    work_counts: FloatArray = counts_arr
    work_nhood_labels: Maybe[Sequence[str]] = nhood_group_labels
    work_prior_weights: Maybe[FloatArray] = prior_pair_weights

    if use_super:
        cg_res = coarse_grain_nhoods_by_graph(
            nhoods_mat=nhoods_mat,
            count_matrix=counts_arr,
            cell_sample_codes=cell_sample_codes,
            sample_library_sizes=lib_maybe,
            nhood_group_labels=nhood_group_labels,
            cell_types=cell_types,
            min_jaccard=config.prefilter.super_nhood_min_jaccard,
            leiden_resolution=config.prefilter.super_nhood_resolution,
        )
        match cg_res:
            case Failure(err):
                return Failure(err)
            case Success(cg_data):
                work_nhoods = cg_data.nhoods_super
                work_counts = cg_data.counts_super
                work_nhood_labels = Some(cg_data.super_nhood_labels)
                if isinstance(prior_pair_weights, Some):
                    w_in = np.asarray(prior_pair_weights.unwrap(), dtype=np.float64)
                    if w_in.shape == (n_nhoods_orig, n_nhoods_orig):
                        work_prior_weights = Some(
                            _pool_prior_weights_to_super(w_in, cg_data.member_nhood_indices)
                        )

    overlap_res = compute_nhood_overlap_graph(
        nhoods_mat=work_nhoods,
        cell_types=cell_types,
        config=config.overlap,
    )
    match overlap_res:
        case Failure(err):
            return Failure(err)
        case Success(raw_overlap_graph):
            pass

    match work_nhood_labels:
        case Some(lbl_seq):
            if len(lbl_seq) != raw_overlap_graph.n_nhoods:
                return Failure(
                    f"nhood_group_labels length ({len(lbl_seq)}) does not match neighborhoods ({raw_overlap_graph.n_nhoods})"
                )
            resolved_lbls = tuple(str(x) for x in lbl_seq)
            has_lbls = any(lbl != "unannotated" for lbl in resolved_lbls)
            overlap_graph = OverlapGraph(
                n_nhoods=raw_overlap_graph.n_nhoods,
                nhood_sizes=raw_overlap_graph.nhood_sizes,
                jaccard_matrix=raw_overlap_graph.jaccard_matrix,
                spatial_weights_1d=raw_overlap_graph.spatial_weights_1d,
                dominant_cell_types=resolved_lbls,
                cell_type_purities=tuple(1.0 for _ in resolved_lbls),
                has_cell_types=has_lbls or raw_overlap_graph.has_cell_types,
            )
        case Maybe.empty:
            overlap_graph = raw_overlap_graph
        case _ as unreachable_lbl:
            return Failure(f"Unexpected Maybe state: {unreachable_lbl}")

    estimator_fn: PresenceEstimator = presence_estimator.value_or(estimate_nhood_presence)
    pres_res = estimator_fn(work_counts, lib_maybe, config.presence)
    match pres_res:
        case Failure(err):
            return Failure(err)
        case Success(raw_presence_rep):
            pass

    if config.mixing_qc.enabled:
        resolved_sids: Sequence[str] = (
            sample_ids.unwrap()
            if isinstance(sample_ids, Some)
            else tuple(f"sample_{j}" for j in range(n_samples))
        )
        mix_res = compute_nhood_mixing_metrics(
            nhoods_mat=work_nhoods,
            count_matrix=work_counts,
            sample_ids=resolved_sids,
            batch_labels=strata_labels,
            nhood_labels=Some(overlap_graph.dominant_cell_types),
            config=config.mixing_qc,
        )
        match mix_res:
            case Failure(err):
                return Failure(err)
            case Success(mix_df):
                qc_mask_global = mix_df["passes_mixing_qc"].to_numpy()
                keep_qc = qc_mask_global[raw_presence_rep.active_nhood_indices]
                qc_presence_rep = subset_presence_representation(raw_presence_rep, keep_qc)
    else:
        qc_presence_rep = raw_presence_rep

    effective_cohorts: Maybe[Sequence[str]] = (
        cohort_labels if isinstance(cohort_labels, Some) else strata_labels
    )
    if (
        config.presence.residualize_library_depth
        and isinstance(cohort_labels, Some)
        and len({str(x).strip() for x in cohort_labels.unwrap()}) > 1
        and int(qc_presence_rep.active_nhood_indices.shape[0]) > 0
    ):
        act_c = work_counts[qc_presence_rep.active_nhood_indices, :]
        eff_libs = qc_presence_rep.effective_library_sizes
        raw_p = act_c / np.maximum(np.expand_dims(eff_libs, axis=0), 1.0)
        raw_arc = np.arcsin(np.sqrt(np.clip(raw_p, 0.0, 1.0)))
        raw_cpm = np.log2(1e6 * raw_p + config.presence.prior_cpm_pseudocount)
        cpm_g_res = depth_residualize_continuous_features(
            raw_cpm, eff_libs, preserve_row_means=True, group_labels=cohort_labels
        )
        arc_g_res = depth_residualize_continuous_features(
            raw_arc, eff_libs, preserve_row_means=True, group_labels=cohort_labels
        )
        match (cpm_g_res, arc_g_res):
            case (Success(g_cpm), Success(g_arc)):
                qc_presence_rep = PresenceRepresentation(
                    active_nhood_indices=qc_presence_rep.active_nhood_indices,
                    presence_prob=qc_presence_rep.presence_prob,
                    absence_prob=qc_presence_rep.absence_prob,
                    observability=qc_presence_rep.observability,
                    raw_binary_presence=qc_presence_rep.raw_binary_presence,
                    arcsin_sqrt_prop=np.where(
                        act_c > 0.0, np.clip(g_arc, 1e-8, 0.5 * np.pi), 0.0
                    ),
                    log2_cpm=g_cpm,
                    effective_library_sizes=eff_libs,
                )
            case (Failure(err), _):
                return Failure(err)
            case (_, Failure(err)):
                return Failure(err)
            case _ as unreachable_grp:
                return Failure(f"Unexpected Result state: {unreachable_grp}")

    if (
        config.presence.mask_unobserved_compartments
        and isinstance(effective_cohorts, Some)
        and overlap_graph.has_cell_types
    ):
        obs_mask_res = build_compartment_observability_mask(
            presence=qc_presence_rep,
            nhood_labels=overlap_graph.dominant_cell_types,
            cohort_labels=effective_cohorts.unwrap(),
            min_cohort_compartment_share=config.presence.min_cohort_compartment_share,
        )
        match obs_mask_res:
            case Failure(err):
                return Failure(err)
            case Success(obs_presence_rep):
                pass
    else:
        obs_presence_rep = qc_presence_rep

    if use_marginal_filter:
        marg_res = select_informative_nhoods_outcome_blind(obs_presence_rep, config.prefilter)
        match marg_res:
            case Failure(err):
                return Failure(err)
            case Success(keep_act):
                presence_rep = subset_presence_representation(obs_presence_rep, keep_act)
    else:
        presence_rep = obs_presence_rep

    effective_screen_cfg = _resolve_effective_screen_config(config)

    run_quadrant = config.method in (PairwiseMethod.QUADRANT_2X2, PairwiseMethod.HYBRID_ALL)
    cont_axes: tuple[ContinuousAxis, ...] = ()
    match config.method:
        case PairwiseMethod.QUADRANT_2X2:
            cont_axes = ()
        case PairwiseMethod.LOG_RATIO_BALANCE:
            cont_axes = (ContinuousAxis.LOG_RATIO_BALANCE,)
        case PairwiseMethod.CONTINUOUS_INTERACTION:
            cont_axes = (ContinuousAxis.CO_ABUNDANCE_PRODUCT,)
        case PairwiseMethod.HYBRID_ALL:
            cont_axes = (ContinuousAxis.LOG_RATIO_BALANCE, ContinuousAxis.CO_ABUNDANCE_PRODUCT)
        case _ as unreachable:
            assert_never(unreachable)

    if run_quadrant:
        quad_res = screen_quadrant_pairs(
            presence=presence_rep,
            overlap=overlap_graph,
            response_binary=response_binary,
            strata_codes=strata_codes,
            overlap_config=config.overlap,
            screen_config=effective_screen_cfg,
            prior_pair_weights=work_prior_weights,
        )
        match quad_res:
            case Failure(err):
                return Failure(err)
            case Success((screened_df, n_eval_pairs)):
                pass

        refiner_fn: CandidatePairRefiner = refiner.value_or(refine_candidate_pairs_firth)
        ref_res = refiner_fn(
            screened_df,
            presence_rep,
            response_binary,
            strata_codes,
            config.refinement,
            effective_screen_cfg.fdr_threshold,
        )
        match ref_res:
            case Failure(err):
                return Failure(err)
            case Success(refined_df):
                pass

        mod_res = group_nhood_pairs(refined_df, overlap_graph, config.modules)
        match mod_res:
            case Failure(err):
                return Failure(err)
            case Success((final_pairs_df, modules_df)):
                pass
    else:
        final_pairs_df = _empty_quadrant_df().with_columns(
            pl.Series("log_or_quadrant_firth", [], dtype=pl.Float64),
            pl.Series("pvalue_quadrant_firth", [], dtype=pl.Float64),
            pl.Series("pvalue_iut_firth", [], dtype=pl.Float64),
            pl.Series("is_refined", [], dtype=pl.Boolean),
            pl.Series("pair_module_id", [], dtype=pl.Int64),
        )
        modules_df = _empty_modules_df()
        n_eval_pairs = 0

    if cont_axes:
        cont_res = screen_continuous_pairs(
            presence=presence_rep,
            overlap=overlap_graph,
            response_binary=response_binary,
            strata_codes=strata_codes,
            axes=cont_axes,
            overlap_config=config.overlap,
            screen_config=effective_screen_cfg,
            prior_pair_weights=work_prior_weights,
        )
        match cont_res:
            case Failure(err):
                return Failure(err)
            case Success(continuous_df):
                pass
        if n_eval_pairs == 0 and not continuous_df.is_empty():
            n_eval_pairs = int(
                continuous_df.select(pl.struct(["nhood_u", "nhood_v"]).n_unique()).item()
            )
    else:
        continuous_df = _empty_continuous_df()

    n_sig_quad = int(final_pairs_df.filter(pl.col("is_significant")).height)
    n_sig_cont = int(continuous_df.filter(pl.col("is_significant")).height)

    return Success(
        PairwiseAnalysisResult(
            pairs_df=final_pairs_df,
            continuous_pairs_df=continuous_df,
            modules_df=modules_df,
            n_total_nhoods=n_nhoods_orig,
            n_active_nhoods=int(presence_rep.active_nhood_indices.shape[0]),
            n_evaluated_pairs=n_eval_pairs,
            n_significant_quadrant_pairs=n_sig_quad,
            n_significant_continuous_pairs=n_sig_cont,
            n_modules=int(modules_df.height),
        )
    )


def test_nhood_pairs_adata(
    adata: ad.AnnData,
    sample_col: str,
    design_col: str,
    dataset_col: Maybe[str] = Nothing,
    cell_type_col: Maybe[str] = Nothing,
    config: PairwiseConfig = PairwiseConfig(),
    presence_estimator: Maybe[PresenceEstimator] = Nothing,
    refiner: Maybe[CandidatePairRefiner] = Nothing,
) -> Result[PairwiseAnalysisResult, str]:
    """Pure functional AnnData entry point for pairwise neighborhood differential association.

    Reads `adata.obsm['nhoods']` and `adata.obs` without mutating `adata`. Computes
    exact donor total cell library sizes from `adata.obs[sample_col]`.
    """
    if "nhoods" not in adata.obsm:
        return Failure("Neighborhood incidence matrix 'nhoods' not found in adata.obsm")
    if sample_col not in adata.obs.columns:
        return Failure(f"Sample column '{sample_col}' not found in adata.obs")
    if design_col not in adata.obs.columns:
        return Failure(f"Design column '{design_col}' not found in adata.obs")

    nhoods_mat = sp.csc_matrix(adata.obsm["nhoods"])
    n_cells = int(adata.n_obs)

    sample_vals = tuple(str(x) for x in adata.obs[sample_col].tolist())
    unique_samples = tuple(sorted(set(sample_vals)))
    n_samples = len(unique_samples)
    if n_samples < 4:
        return Failure(f"At least 4 biological samples are required, found {n_samples}")

    sample_to_idx = {s_name: idx for idx, s_name in enumerate(unique_samples)}
    sample_codes = np.array([sample_to_idx[s_name] for s_name in sample_vals], dtype=np.int64)

    # Pure sparse count aggregation C = A^T S of shape (V_nhoods, J_samples)
    sample_ind = sp.csr_matrix(
        (np.ones(n_cells, dtype=np.float64), (np.arange(n_cells, dtype=np.int64), sample_codes)),
        shape=(n_cells, n_samples),
    )
    count_matrix = np.asarray(nhoods_mat.T.dot(sample_ind).toarray(), dtype=np.float64)
    sample_lib_sizes = np.bincount(sample_codes, minlength=n_samples).astype(np.float64)

    # Extract sample-level clinical response and verify within-sample consistency
    resp_vals = tuple(str(x).strip().lower() for x in adata.obs[design_col].tolist())
    sample_resp_sets: list[set[str]] = [set() for _ in range(n_samples)]
    for c_idx, s_code in enumerate(sample_codes.tolist()):
        sample_resp_sets[s_code].add(resp_vals[c_idx])

    ambiguous = [
        unique_samples[idx]
        for idx, rset in enumerate(sample_resp_sets)
        if len(rset) != 1
    ]
    if ambiguous:
        return Failure(
            f"Samples have multiple conflicting '{design_col}' labels in adata.obs: {ambiguous[:5]}"
        )
    response_labels = tuple(next(iter(rset)) for rset in sample_resp_sets)

    strata_maybe: Maybe[Sequence[str]] = Nothing
    match dataset_col:
        case Maybe.empty:
            strata_maybe = Nothing
        case Some(ds_col):
            if ds_col not in adata.obs.columns:
                return Failure(f"Dataset column '{ds_col}' not found in adata.obs")
            ds_vals = tuple(str(x) for x in adata.obs[ds_col].tolist())
            sample_ds = [""] * n_samples
            for c_idx, s_code in enumerate(sample_codes.tolist()):
                sample_ds[s_code] = ds_vals[c_idx]
            strata_maybe = Some(tuple(sample_ds))
        case _ as unreachable:
            return Failure(f"Unexpected Maybe state: {unreachable}")

    ct_maybe: Maybe[Sequence[str]] = Nothing
    match cell_type_col:
        case Maybe.empty:
            ct_maybe = Nothing
        case Some(ct_col):
            if ct_col not in adata.obs.columns:
                return Failure(f"Cell type column '{ct_col}' not found in adata.obs")
            ct_maybe = Some(tuple(str(x) for x in adata.obs[ct_col].tolist()))
        case _ as unreachable:
            return Failure(f"Unexpected Maybe state: {unreachable}")

    return analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=count_matrix,
        response_labels=response_labels,
        sample_library_sizes=Some(tuple(float(x) for x in sample_lib_sizes.tolist())),
        strata_labels=strata_maybe,
        cell_types=ct_maybe,
        cell_sample_codes=Some(sample_codes),
        sample_ids=Some(unique_samples),
        config=config,
        presence_estimator=presence_estimator,
        refiner=refiner,
    )


setattr(test_nhood_pairs_adata, "__test__", False)


def project_pair_to_cells(
    nhoods_mat: SparseMatrix,
    nhood_u: int,
    nhood_v: int,
    cell_ids: Maybe[Sequence[str]] = Nothing,
) -> Result[pl.DataFrame, str]:
    """Pure projection of a neighborhood pair (nhood_u, nhood_v) onto single cells."""
    n_cells, n_nhoods = int(nhoods_mat.shape[0]), int(nhoods_mat.shape[1])
    if nhood_u < 0 or nhood_u >= n_nhoods or nhood_v < 0 or nhood_v >= n_nhoods:
        return Failure(
            f"Neighborhood indices ({nhood_u}, {nhood_v}) out of bounds for n_nhoods={n_nhoods}"
        )

    csc = sp.csc_matrix(nhoods_mat)
    in_u = np.asarray(csc[:, nhood_u].toarray(), dtype=np.float64).ravel() > 0.0
    in_v = np.asarray(csc[:, nhood_v].toarray(), dtype=np.float64).ravel() > 0.0

    membership = np.where(
        in_u & in_v,
        "shared",
        np.where(in_u, "nhood_u", np.where(in_v, "nhood_v", "background")),
    )

    match cell_ids:
        case Maybe.empty:
            ids = tuple(f"cell_{i}" for i in range(n_cells))
        case Some(seq):
            if len(seq) != n_cells:
                return Failure(f"cell_ids length ({len(seq)}) does not match n_cells ({n_cells})")
            ids = tuple(str(x) for x in seq)
        case _ as unreachable:
            return Failure(f"Unexpected Maybe state: {unreachable}")

    return Success(
        pl.DataFrame(
            {
                "cell_id": list(ids),
                "in_nhood_u": in_u,
                "in_nhood_v": in_v,
                "pair_membership": membership,
            }
        )
    )
