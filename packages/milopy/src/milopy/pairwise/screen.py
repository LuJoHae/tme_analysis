"""Memory-bounded vectorized BLAS screening for 2x2 neighborhood quadrants and continuous axes.

Evaluates all valid upper-triangular neighborhood pairs (u < v) in row blocks using:
- Stratified Cochran-Mantel-Haenszel (CMH) / Haldane-Anscombe log-odds & Robins-Breslow-Greenland SEs
- Berger's Intersection-Union Test (IUT) over Hamming-1 flip neighbors calibrated via Symmetrized
  Decoupled Partner-Permutation (DMPP) Knockoff+ step-up empirical FDP to block singleton piggybacking
- Outcome-blind marginal & co-positivity testability filtering
- Saturated 2x2 three-way epistatic interaction contrasts
- 2D product-graph weighted Benjamini-Hochberg FDR and optional stratified donor permutations
- Scale-invariant Aitchison log-ratio balance (+-/-+) and standardized co-abundance coupling (++/--)
"""

from __future__ import annotations

from typing import NamedTuple, TypeAlias, assert_never

import numpy as np
from numpy.typing import NDArray
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy.stats import norm, t as student_t

from .overlap import (
    build_block_pair_mask,
    compute_weighted_bh_fdr,
    subset_overlap_graph,
)
from .prefilter import (
    compute_tarone_max_iut_z_candidates_f32,
    compute_tarone_max_iut_z_strata,
    tarone_z_threshold_from_p_ceiling,
)
from .types import (
    ALL_QUADRANT_STATES,
    BoolArray,
    ContinuousAxis,
    EnrichmentDirection,
    FloatArray,
    IntArray,
    OverlapConfig,
    OverlapGraph,
    PairScopeFilter,
    PresenceRepresentation,
    QuadrantState,
    ScreenConfig,
    TmeCompartment,
    adjacent_flips,
    classify_tme_compartment,
)

Float32Array: TypeAlias = NDArray[np.float32]

_DMPP_GRID_MIN_Z: float = 1.5
_DMPP_GRID_MAX_Z: float = 12.0
_DMPP_GRID_BINS: int = 512


class _StratumSlice(NamedTuple):
    """Pre-sliced per-stratum matrices for O(V^2 J) stratified BLAS accumulation."""

    p_r: FloatArray
    q_r: FloatArray
    m_r: FloatArray
    p_nr: FloatArray
    q_nr: FloatArray
    m_nr: FloatArray
    p_all: FloatArray
    q_all: FloatArray


class _StratumSliceF32(NamedTuple):
    """Compact float32 stratum matrices for fast BLAS sgemm decoupled null evaluation."""

    p_r: Float32Array
    m_r: Float32Array
    p_nr: Float32Array
    m_nr: Float32Array


class _ChunkQuadrantArrays(NamedTuple):
    """Flat evaluated candidate arrays emitted by one row chunk."""

    u_local: IntArray
    v_local: IntArray
    quadrant_idx: IntArray
    jaccard: FloatArray
    scope: NDArray[np.str_]
    spatial_weight: FloatArray
    prior_multiplier: FloatArray
    n_r_q: FloatArray
    n_nr_q: FloatArray
    n_eff_total: FloatArray
    raw_donors_enriched: IntArray
    log_or_quadrant: FloatArray
    se_quadrant: FloatArray
    pvalue_quadrant: FloatArray
    min_cond_log_or: FloatArray
    z_iut_signed: FloatArray
    pvalue_iut: FloatArray
    quad_testable: BoolArray
    log_or_interaction: FloatArray
    se_interaction: FloatArray
    pvalue_interaction: FloatArray


def _build_strata_slices(
    presence: PresenceRepresentation,
    response_binary: FloatArray,
    strata_codes: Maybe[IntArray],
) -> tuple[list[_StratumSlice], list[IntArray]]:
    """Splits P, Q, M into responder/non-responder blocks per informative stratum."""
    p_mat = presence.presence_prob
    q_mat = presence.absence_prob
    m_mat = presence.observability
    n_samples = int(response_binary.shape[0])

    match strata_codes:
        case Some(codes):
            unique_codes = np.unique(codes)
            sample_groups = [np.flatnonzero(codes == c).astype(np.int64) for c in unique_codes]
        case Maybe.empty:
            sample_groups = [np.arange(n_samples, dtype=np.int64)]
        case _ as unreachable:
            sample_groups = [np.arange(n_samples, dtype=np.int64)]
            _ = unreachable

    slices: list[_StratumSlice] = []
    informative_groups: list[IntArray] = []
    for idx_grp in sample_groups:
        y_grp = response_binary[idx_grp]
        r_mask = y_grp > 0.5
        nr_mask = ~r_mask
        if bool(np.any(r_mask)) and bool(np.any(nr_mask)):
            r_idx = idx_grp[r_mask]
            nr_idx = idx_grp[nr_mask]
            slices.append(
                _StratumSlice(
                    p_r=p_mat[:, r_idx],
                    q_r=q_mat[:, r_idx],
                    m_r=m_mat[:, r_idx],
                    p_nr=p_mat[:, nr_idx],
                    q_nr=q_mat[:, nr_idx],
                    m_nr=m_mat[:, nr_idx],
                    p_all=p_mat[:, idx_grp],
                    q_all=q_mat[:, idx_grp],
                )
            )
            informative_groups.append(idx_grp)
    return slices, informative_groups


def _cmh_2x2_from_strata(
    a_list: list[FloatArray],
    b_list: list[FloatArray],
    c_list: list[FloatArray],
    d_list: list[FloatArray],
    pseudocount: float,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Computes Mantel-Haenszel log-OR, Robins-Breslow-Greenland SE, and exact CMH score Z across strata.

    Uses Haldane-Anscombe (+pseudocount/K) smoothing for finite log-OR and RBG SE estimation while
    computing exact Cochran-Mantel-Haenszel score Z-statistics on unsmoothed counts to avoid
    Hauck-Donner Wald collapse on zero-cell tables.
    """
    k_strata = len(a_list)
    delta = pseudocount / float(max(k_strata, 1))
    shape = a_list[0].shape

    r_mh = np.zeros(shape, dtype=np.float64)
    s_mh = np.zeros(shape, dtype=np.float64)
    t1 = np.zeros(shape, dtype=np.float64)
    t2 = np.zeros(shape, dtype=np.float64)
    t3 = np.zeros(shape, dtype=np.float64)
    cmh_num = np.zeros(shape, dtype=np.float64)
    cmh_var = np.zeros(shape, dtype=np.float64)

    for k in range(k_strata):
        a_raw, b_raw, c_raw, d_raw = a_list[k], b_list[k], c_list[k], d_list[k]
        n_raw = a_raw + b_raw + c_raw + d_raw
        safe_n = np.maximum(n_raw, 1e-12)
        r1 = a_raw + b_raw
        r0 = c_raw + d_raw
        c1 = a_raw + c_raw
        c0 = b_raw + d_raw
        exp_a = (r1 * c1) / safe_n
        denom_var = safe_n * safe_n * np.maximum(n_raw - 1.0, 1.0)
        var_a = np.where(
            n_raw > 1.0,
            (r1 * r0 * c1 * c0) / np.maximum(denom_var, 1e-12),
            0.0,
        )
        cmh_num = cmh_num + np.where(n_raw > 1.0, a_raw - exp_a, 0.0)
        cmh_var = cmh_var + np.maximum(var_a, 0.0)

        ak = a_raw + delta
        bk = b_raw + delta
        ck = c_raw + delta
        dk = d_raw + delta
        nk = np.maximum(ak + bk + ck + dk, 1e-12)

        ad_n = (ak * dk) / nk
        bc_n = (bk * ck) / nk
        pk = (ak + dk) / nk
        qk = (bk + ck) / nk

        r_mh = r_mh + ad_n
        s_mh = s_mh + bc_n
        t1 = t1 + pk * ad_n
        t2 = t2 + pk * bc_n + qk * ad_n
        t3 = t3 + qk * bc_n

    safe_r = np.maximum(r_mh, 1e-300)
    safe_s = np.maximum(s_mh, 1e-300)
    log_or = np.log(safe_r / safe_s)
    variance = (
        t1 / (2.0 * safe_r * safe_r)
        + t2 / (2.0 * safe_r * safe_s)
        + t3 / (2.0 * safe_s * safe_s)
    )
    se = np.sqrt(np.maximum(variance, 1e-15))
    z_score = np.where(
        cmh_var > 1e-12,
        cmh_num / np.sqrt(np.maximum(cmh_var, 1e-15)),
        0.0,
    )
    return log_or, se, z_score


def _quadrant_index(state: QuadrantState) -> int:
    match state:
        case QuadrantState.BOTH_PRESENT:
            return 0
        case QuadrantState.FIRST_ONLY:
            return 1
        case QuadrantState.SECOND_ONLY:
            return 2
        case QuadrantState.BOTH_MISSING:
            return 3
        case _ as unreachable:
            assert_never(unreachable)


def _stepup_knockoff_plus_fdr(
    stat_obs: FloatArray,
    weights_obs: FloatArray,
    null_hist_weights: FloatArray,
    grid_edges: FloatArray,
    n_null_draws: int,
) -> FloatArray:
    """Computes spatial-weighted multi-decoy Knockoff+ step-up empirical FDP q-values.

    Uses the finite-sample super-martingale estimator (Gimenez & Zou 2019; Emery & Keich 2020):
        FDP_+(z) = (1 + sum_{b=1}^B V_b(z)) / (B * max(R(z), 1))
    with conservative lower-bin lookup on `grid_edges` so null tail counts are bounded from below.
    """
    n_hyp = int(stat_obs.shape[0])
    if n_hyp == 0:
        return np.empty(0, dtype=np.float64)

    w_raw = np.maximum(np.asarray(weights_obs, dtype=np.float64), 0.0)
    w_sum = float(np.sum(w_raw))
    scale = (float(n_hyp) / w_sum) if w_sum > 0.0 else 1.0
    w_norm = (w_raw * scale) if w_sum > 0.0 else np.ones(n_hyp, dtype=np.float64)

    b_denom = float(max(n_null_draws, 1))
    null_tail_w = np.cumsum(null_hist_weights[::-1])[::-1] * scale
    mean_null_plus = (1.0 + null_tail_w) / b_denom

    pos_mask = stat_obs >= float(grid_edges[0])
    if not bool(np.any(pos_mask)):
        return np.ones(n_hyp, dtype=np.float64)

    s_pos = stat_obs[pos_mask]
    w_pos = w_norm[pos_mask]
    order_desc = np.argsort(-s_pos, kind="mergesort")
    s_sorted = s_pos[order_desc]
    cum_r = np.cumsum(w_pos[order_desc])

    bin_idx = np.clip(
        np.searchsorted(grid_edges[:-1], s_sorted, side="right") - 1,
        0,
        int(null_tail_w.shape[0]) - 1,
    )
    raw_fdp = mean_null_plus[bin_idx] / np.maximum(cum_r, 1.0)
    q_sorted = np.clip(np.minimum.accumulate(raw_fdp[::-1])[::-1], 0.0, 1.0)

    q_out = np.ones(n_hyp, dtype=np.float64)
    q_out[np.flatnonzero(pos_mask)[order_desc]] = q_sorted
    return q_out


def _build_symmetric_offdiag_mask(
    overlap_active: OverlapGraph,
    u_start: int,
    u_end: int,
    config: OverlapConfig,
) -> BoolArray:
    """Builds the symmetric off-diagonal validity mask (u != v) for decoupled null evaluation."""
    v_total = overlap_active.n_nhoods
    block_rows = u_end - u_start
    row_coords = np.expand_dims(np.arange(u_start, u_end, dtype=np.int64), axis=1)
    col_coords = np.expand_dims(np.arange(v_total, dtype=np.int64), axis=0)
    off_diag = col_coords != row_coords

    jaccard_block = np.asarray(
        overlap_active.jaccard_matrix[u_start:u_end, :].toarray(),
        dtype=np.float64,
    )
    disjoint_ok = jaccard_block <= (config.max_jaccard_overlap + 1e-12)

    if overlap_active.has_cell_types:
        ct_all = np.asarray(overlap_active.dominant_cell_types, dtype=str)
        is_same_ct = np.expand_dims(ct_all[u_start:u_end], axis=1) == np.expand_dims(ct_all, axis=0)
        match config.scope_filter:
            case PairScopeFilter.ALL:
                scope_ok = np.ones((block_rows, v_total), dtype=bool)
            case PairScopeFilter.CROSS_CELL_TYPE_ONLY:
                scope_ok = ~is_same_ct
            case PairScopeFilter.CROSS_COMPARTMENT_ONLY:
                comp_all = np.asarray(
                    [classify_tme_compartment(str(c)).value for c in overlap_active.dominant_cell_types],
                    dtype=str,
                )
                comp_u = np.expand_dims(comp_all[u_start:u_end], axis=1)
                comp_v = np.expand_dims(comp_all, axis=0)
                scope_ok = (
                    (comp_u != comp_v)
                    & (comp_u != TmeCompartment.OTHER.value)
                    & (comp_v != TmeCompartment.OTHER.value)
                )
            case PairScopeFilter.WITHIN_CELL_TYPE_ONLY:
                scope_ok = is_same_ct
            case _ as unreachable:
                assert_never(unreachable)
    else:
        match config.scope_filter:
            case PairScopeFilter.ALL:
                scope_ok = np.ones((block_rows, v_total), dtype=bool)
            case (
                PairScopeFilter.CROSS_CELL_TYPE_ONLY
                | PairScopeFilter.CROSS_COMPARTMENT_ONLY
                | PairScopeFilter.WITHIN_CELL_TYPE_ONLY
            ):
                scope_ok = np.zeros((block_rows, v_total), dtype=bool)
            case _ as unreachable:
                assert_never(unreachable)

    return np.asarray(off_diag & disjoint_ok & scope_ok, dtype=bool)


def _cmh_accum_step_f32(
    a: Float32Array,
    b: Float32Array,
    c: Float32Array,
    d: Float32Array,
    num_acc: Float32Array,
    var_acc: Float32Array,
) -> tuple[Float32Array, Float32Array]:
    """Accumulates one stratum's unsmoothed CMH numerator and hypergeometric variance in float32."""
    n = a + b + c + d
    r1 = a + b
    r0 = c + d
    c1 = a + c
    c0 = b + d
    sn = np.maximum(n, np.float32(1e-6))
    active = n > np.float32(1.0)
    num_k = np.where(active, (a * d - b * c) / sn, np.float32(0.0))
    den = sn * sn * np.maximum(n - np.float32(1.0), np.float32(1.0))
    var_k = np.where(
        active,
        (r1 * r0 * c1 * c0) / np.maximum(den, np.float32(1e-6)),
        np.float32(0.0),
    )
    return num_acc + num_k, var_acc + np.maximum(var_k, np.float32(0.0))


def _resolve_dmpp_draw_count(
    v_act: int,
    fdr_threshold: float,
    dmpp_null_draws: int = 0,
) -> int:
    """Determines decoupled partner-permutation count B so 1/B is strictly below fdr_threshold."""
    min_for_alpha = int(np.ceil(1.0 / max(fdr_threshold, 0.02))) + 2
    if dmpp_null_draws > 0:
        return max(int(dmpp_null_draws), min_for_alpha)
    return max(25, min_for_alpha) if v_act <= 250 else max(12, min_for_alpha)


def _apply_prior_weight_to_z(z_vals: FloatArray, prior_mult: FloatArray) -> FloatArray:
    """Monotonically maps normal Z-statistics via Genovese-Roeder-Wasserman p / w_prior scaling."""
    if int(z_vals.shape[0]) == 0:
        return z_vals
    pos = z_vals > 0.0
    if not bool(np.any(pos)):
        return z_vals
    out = np.array(z_vals, dtype=np.float64, copy=True)
    p_two = np.clip(2.0 * norm.sf(z_vals[pos]), 1e-300, 1.999999)
    p_adj = np.clip(p_two / np.maximum(prior_mult[pos], 1e-12), 1e-300, 1.999999)
    out[pos] = np.asarray(norm.isf(0.5 * p_adj), dtype=np.float64)
    return out


def _normalize_active_prior_matrix(
    prior_pair_weights: Maybe[FloatArray],
    active_nhood_indices: IntArray,
    overlap_active: OverlapGraph,
    overlap_config: OverlapConfig,
) -> Result[Maybe[FloatArray], str]:
    """Subsets and normalizes optional external prior weights over valid upper-triangular pairs."""
    match prior_pair_weights:
        case Maybe.empty:
            return Success(Nothing)
        case Some(w_full):
            w_arr = np.asarray(w_full, dtype=np.float64)
            if w_arr.ndim != 2 or w_arr.shape[0] != w_arr.shape[1]:
                return Failure(f"prior_pair_weights must be square 2D, got shape {w_arr.shape}")
            if int(active_nhood_indices.shape[0]) == 0:
                return Success(Some(np.empty((0, 0), dtype=np.float64)))
            if int(np.max(active_nhood_indices)) >= int(w_arr.shape[0]):
                return Failure(
                    f"prior_pair_weights dimension ({w_arr.shape[0]}) smaller than max active nhood index"
                )
            w_act = 0.5 * (
                w_arr[active_nhood_indices, :][:, active_nhood_indices]
                + w_arr[active_nhood_indices, :][:, active_nhood_indices].T
            )
            w_act = np.maximum(w_act, 0.0)
            np.fill_diagonal(w_act, 0.0)
            v_act = overlap_active.n_nhoods
            if v_act < 2:
                return Success(Some(w_act))
            valid_upper, _, _ = build_block_pair_mask(overlap_active, 0, v_act, overlap_config)
            pos_valid = valid_upper & (w_act > 0.0)
            if bool(np.any(pos_valid)):
                mean_pos = float(np.mean(w_act[pos_valid]))
                if mean_pos > 0.0:
                    w_act = np.where(w_act > 0.0, w_act / mean_pos, 0.0)
            return Success(Some(w_act))
        case _ as unreachable:
            return Failure(f"Unexpected Maybe state: {unreachable}")


def _compute_decoupled_quadrant_iut_fdr(
    z_iut_testable: FloatArray,
    sw_testable: FloatArray,
    prior_mult_testable: FloatArray,
    has_external_prior: bool,
    w_prior_active: Maybe[FloatArray],
    tarone_z_floor: float,
    presence: PresenceRepresentation,
    strata: list[_StratumSlice],
    informative_groups: list[IntArray],
    response_binary: FloatArray,
    overlap_active: OverlapGraph,
    overlap_config: OverlapConfig,
    screen_config: ScreenConfig,
) -> FloatArray:
    """Computes Symmetrized Decoupled Partner-Permutation (DMPP) Knockoff+ q-values for Berger's IUT.

    Evaluates (u_real, v_perm) across all symmetric off-diagonal non-overlapping pairs (u != v),
    where u < v contributes 0.5 weight for H_10 (u real, v permuted) and u > v contributes 0.5
    weight for H_01 (u permuted, v real) identically. Preserves exact 1D singleton signals,
    discrete hypergeometric support, kNN manifold correlation, and pair-exchange symmetry.
    """
    n_test = int(z_iut_testable.shape[0])
    if n_test == 0:
        return np.empty(0, dtype=np.float64)

    v_act = overlap_active.n_nhoods
    b_draws = _resolve_dmpp_draw_count(
        v_act,
        screen_config.fdr_threshold,
        screen_config.dmpp_null_draws,
    )
    grid_edges = np.linspace(_DMPP_GRID_MIN_Z, _DMPP_GRID_MAX_Z, _DMPP_GRID_BINS + 1, dtype=np.float64)
    null_hist = np.zeros(_DMPP_GRID_BINS, dtype=np.float64)

    real_f32 = [
        _StratumSliceF32(
            p_r=sl.p_r.astype(np.float32),
            m_r=sl.m_r.astype(np.float32),
            p_nr=sl.p_nr.astype(np.float32),
            m_nr=sl.m_nr.astype(np.float32),
        )
        for sl in strata
    ]
    k_strata = len(real_f32)

    # Precompute symmetric block masks and 1D spatial weights
    perm_block = max(screen_config.chunk_size, 1024)
    block_ranges = [
        (start, min(start + perm_block, v_act))
        for start in range(0, v_act, perm_block)
    ]
    w_prior_mat: FloatArray = (
        w_prior_active.unwrap()
        if isinstance(w_prior_active, Some)
        else np.empty((0, 0), dtype=np.float64)
    )
    sym_masks = [
        (
            _build_symmetric_offdiag_mask(overlap_active, s, e, overlap_config)
            & (w_prior_mat[s:e, :] > 0.0)
        )
        if has_external_prior
        else _build_symmetric_offdiag_mask(overlap_active, s, e, overlap_config)
        for s, e in block_ranges
    ]
    if not any(bool(np.any(m)) for m in sym_masks):
        return np.ones(n_test, dtype=np.float64)

    w1d = (
        overlap_active.spatial_weights_1d.astype(np.float64)
        if screen_config.use_spatial_fdr_weights
        else np.ones(v_act, dtype=np.float64)
    )

    n_samples = int(response_binary.shape[0])
    rng = np.random.default_rng(screen_config.permutation_seed)
    p_full_f32 = presence.presence_prob.astype(np.float32)
    m_full_f32 = presence.observability.astype(np.float32)
    stratum_r_nr = [
        (grp[response_binary[grp] > 0.5], grp[response_binary[grp] <= 0.5])
        for grp in informative_groups
    ]
    stratum_arm_counts = [(int(r_idx.shape[0]), int(nr_idx.shape[0])) for r_idx, nr_idx in stratum_r_nr]

    min_eff_f32 = np.float32(screen_config.min_effective_samples)
    min_don_f32 = np.float32(screen_config.min_donors_per_quadrant)
    min_flip_f32 = np.float32(min(1.5, 0.5 * screen_config.min_donors_per_quadrant))
    w_max = float(np.max(w_prior_mat)) if (has_external_prior and w_prior_mat.size > 0) else 1.0
    raw_min_z = (
        float(norm.isf(min(0.499, float(norm.sf(_DMPP_GRID_MIN_Z)) * max(w_max, 1.0))))
        if (has_external_prior and w_max > 1.0)
        else _DMPP_GRID_MIN_Z
    )
    z_floor_f32 = np.float32(max(0.25, raw_min_z - 1e-4))

    for _ in range(b_draws):
        perm_idx = np.arange(n_samples, dtype=np.int64)
        for grp in informative_groups:
            perm_idx[grp] = rng.permutation(grp)

        p_perm = p_full_f32[:, perm_idx]
        m_perm = m_full_f32[:, perm_idx]
        perm_f32 = [
            _StratumSliceF32(
                p_r=p_perm[:, r_idx],
                m_r=m_perm[:, r_idx],
                p_nr=p_perm[:, nr_idx],
                m_nr=m_perm[:, nr_idx],
            )
            for r_idx, nr_idx in stratum_r_nr
        ]

        for b_idx, (start, end) in enumerate(block_ranges):
            vs = sym_masks[b_idx]
            if not bool(np.any(vs)):
                continue

            b_rows = end - start
            shape_blk = (b_rows, v_act)
            num_u1 = np.zeros(shape_blk, dtype=np.float32)
            var_u1 = np.zeros(shape_blk, dtype=np.float32)
            num_u0 = np.zeros(shape_blk, dtype=np.float32)
            var_u0 = np.zeros(shape_blk, dtype=np.float32)
            num_v1 = np.zeros(shape_blk, dtype=np.float32)
            var_v1 = np.zeros(shape_blk, dtype=np.float32)
            num_v0 = np.zeros(shape_blk, dtype=np.float32)
            var_v0 = np.zeros(shape_blk, dtype=np.float32)

            tot_pp = np.zeros(shape_blk, dtype=np.float32)
            tot_pm = np.zeros(shape_blk, dtype=np.float32)
            tot_mp = np.zeros(shape_blk, dtype=np.float32)
            tot_mm = np.zeros(shape_blk, dtype=np.float32)
            eff_tot = np.zeros(shape_blk, dtype=np.float32)

            stratum_quads: list[
                tuple[
                    tuple[Float32Array, Float32Array, Float32Array, Float32Array],
                    tuple[Float32Array, Float32Array, Float32Array, Float32Array],
                    Float32Array,
                    Float32Array,
                ]
            ] = []

            for k in range(k_strata):
                rl = real_f32[k]
                pm = perm_f32[k]
                pr_b, mr_b = rl.p_r[start:end, :], rl.m_r[start:end, :]
                pn_b, mn_b = rl.p_nr[start:end, :], rl.m_nr[start:end, :]

                # Identity Q = M - P yields all 4 quadrants + total observability in 4 sgemm calls/arm
                r_pp = pr_b @ pm.p_r.T
                r_tu = pr_b @ pm.m_r.T
                r_tv = mr_b @ pm.p_r.T
                r_tot = mr_b @ pm.m_r.T
                r_pm = np.maximum(r_tu - r_pp, np.float32(0.0))
                r_mp = np.maximum(r_tv - r_pp, np.float32(0.0))
                r_mm = np.maximum(r_tot - r_tu - r_tv + r_pp, np.float32(0.0))

                n_pp = pn_b @ pm.p_nr.T
                n_tu = pn_b @ pm.m_nr.T
                n_tv = mn_b @ pm.p_nr.T
                n_tot = mn_b @ pm.m_nr.T
                n_pm = np.maximum(n_tu - n_pp, np.float32(0.0))
                n_mp = np.maximum(n_tv - n_pp, np.float32(0.0))
                n_mm = np.maximum(n_tot - n_tu - n_tv + n_pp, np.float32(0.0))

                num_u1, var_u1 = _cmh_accum_step_f32(r_pp, n_pp, r_pm, n_pm, num_u1, var_u1)
                num_u0, var_u0 = _cmh_accum_step_f32(r_mp, n_mp, r_mm, n_mm, num_u0, var_u0)
                num_v1, var_v1 = _cmh_accum_step_f32(r_pp, n_pp, r_mp, n_mp, num_v1, var_v1)
                num_v0, var_v0 = _cmh_accum_step_f32(r_pm, n_pm, r_mm, n_mm, num_v0, var_v0)

                tot_pp = tot_pp + r_pp + n_pp
                tot_pm = tot_pm + r_pm + n_pm
                tot_mp = tot_mp + r_mp + n_mp
                tot_mm = tot_mm + r_mm + n_mm
                eff_tot = eff_tot + r_tot + n_tot
                stratum_quads.append(
                    ((r_pp, r_pm, r_mp, r_mm), (n_pp, n_pm, n_mp, n_mm), r_tot, n_tot)
                )

            ze_u1 = np.where(var_u1 > 1e-8, num_u1 / np.sqrt(np.maximum(var_u1, np.float32(1e-12))), np.float32(0.0))
            ze_u0 = np.where(var_u0 > 1e-8, num_u0 / np.sqrt(np.maximum(var_u0, np.float32(1e-12))), np.float32(0.0))
            ze_v1 = np.where(var_v1 > 1e-8, num_v1 / np.sqrt(np.maximum(var_v1, np.float32(1e-12))), np.float32(0.0))
            ze_v0 = np.where(var_v0 > 1e-8, num_v0 / np.sqrt(np.maximum(var_v0, np.float32(1e-12))), np.float32(0.0))

            quad_specs = (
                (0, 1, 2, tot_pp, tot_pm, tot_mp, ze_u1, ze_v1),
                (1, 0, 3, tot_pm, tot_pp, tot_mm, -ze_u1, ze_v0),
                (2, 3, 0, tot_mp, tot_mm, tot_pp, ze_u0, -ze_v1),
                (3, 2, 1, tot_mm, tot_mp, tot_pm, -ze_u0, -ze_v0),
            )
            for qi, fv_i, fu_i, tq, tfv, tfu, zfv, zfu in quad_specs:
                cand = vs & ((np.minimum(zfv, zfu) >= z_floor_f32) | (np.maximum(zfv, zfu) <= -z_floor_f32))
                ii, jj = np.nonzero(cand)
                if int(ii.shape[0]) == 0:
                    continue
                et_c = eff_tot[ii, jj]
                tq_c = tq[ii, jj]
                ok = (
                    (et_c >= min_eff_f32)
                    & (tq_c >= min_don_f32)
                    & ((et_c - tq_c) >= min_don_f32)
                    & (tfv[ii, jj] >= min_flip_f32)
                    & (tfu[ii, jj] >= min_flip_f32)
                )
                if not bool(np.any(ok)):
                    continue
                ii_k = ii[ok]
                jj_k = jj[ok]
                if tarone_z_floor > 0.0:
                    z_max_cand = compute_tarone_max_iut_z_candidates_f32(
                        stratum_quads_f32=stratum_quads,
                        stratum_arm_counts=stratum_arm_counts,
                        ii=ii_k,
                        jj=jj_k,
                        qi=qi,
                        fv_i=fv_i,
                        fu_i=fu_i,
                    )
                    tar_ok = z_max_cand >= np.float32(tarone_z_floor - 1e-4)
                    if not bool(np.any(tar_ok)):
                        continue
                    ii_k = ii_k[tar_ok]
                    jj_k = jj_k[tar_ok]
                n_cand = int(ii_k.shape[0])
                num_q = np.zeros(n_cand, dtype=np.float32)
                var_q = np.zeros(n_cand, dtype=np.float32)
                for k in range(k_strata):
                    rq_tup, nq_tup, rt_k, nt_k = stratum_quads[k]
                    aq = rq_tup[qi][ii_k, jj_k]
                    bq = nq_tup[qi][ii_k, jj_k]
                    cq = np.maximum(rt_k[ii_k, jj_k] - aq, np.float32(0.0))
                    dq = np.maximum(nt_k[ii_k, jj_k] - bq, np.float32(0.0))
                    num_q, var_q = _cmh_accum_step_f32(aq, bq, cq, dq, num_q, var_q)

                zq = np.where(var_q > 1e-8, num_q / np.sqrt(np.maximum(var_q, np.float32(1e-12))), np.float32(0.0))
                sd = np.where(zq >= np.float32(0.0), np.float32(1.0), np.float32(-1.0))
                zsm_raw = np.minimum(sd * zfv[ii_k, jj_k], sd * zfu[ii_k, jj_k]).astype(np.float64) + 1e-5
                zsm = (
                    _apply_prior_weight_to_z(zsm_raw, w_prior_mat[ii_k + start, jj_k])
                    if has_external_prior
                    else zsm_raw
                )
                hi = zsm >= _DMPP_GRID_MIN_Z
                if bool(np.any(hi)):
                    bins = np.clip(
                        np.searchsorted(grid_edges[:-1], zsm[hi], side="right") - 1,
                        0,
                        _DMPP_GRID_BINS - 1,
                    )
                    w_cand = 0.5 * w1d[ii_k[hi] + start] * w1d[jj_k[hi]]
                    null_hist = null_hist + np.bincount(bins, weights=w_cand, minlength=_DMPP_GRID_BINS)

    obs_stat = (
        _apply_prior_weight_to_z(z_iut_testable, prior_mult_testable)
        if has_external_prior
        else z_iut_testable
    )
    return _stepup_knockoff_plus_fdr(
        stat_obs=obs_stat,
        weights_obs=sw_testable,
        null_hist_weights=null_hist,
        grid_edges=grid_edges,
        n_null_draws=b_draws,
    )


def _evaluate_quadrant_chunk(
    u_start: int,
    u_end: int,
    strata: list[_StratumSlice],
    stratum_arm_counts: list[tuple[int, int]],
    tarone_z_floor: float,
    w_prior_active: Maybe[FloatArray],
    raw_bin_r: FloatArray,
    not_bin_r: FloatArray,
    raw_bin_nr: FloatArray,
    not_bin_nr: FloatArray,
    overlap_active: OverlapGraph,
    overlap_config: OverlapConfig,
    screen_config: ScreenConfig,
) -> _ChunkQuadrantArrays:
    """Evaluates all 4 quadrants on row block [u_start, u_end) across all informative strata."""
    valid_mask, jaccard_block, scope_block = build_block_pair_mask(
        overlap_active=overlap_active,
        u_start=u_start,
        u_end=u_end,
        config=overlap_config,
    )
    row_rel, col_abs = np.nonzero(valid_mask)
    n_pairs = int(row_rel.shape[0])
    if n_pairs == 0:
        empty_f = np.empty(0, dtype=np.float64)
        empty_i = np.empty(0, dtype=np.int64)
        empty_b = np.empty(0, dtype=bool)
        return _ChunkQuadrantArrays(
            u_local=empty_i,
            v_local=empty_i,
            quadrant_idx=empty_i,
            jaccard=empty_f,
            scope=np.empty(0, dtype=str),
            spatial_weight=empty_f,
            prior_multiplier=empty_f,
            n_r_q=empty_f,
            n_nr_q=empty_f,
            n_eff_total=empty_f,
            raw_donors_enriched=empty_i,
            log_or_quadrant=empty_f,
            se_quadrant=empty_f,
            pvalue_quadrant=empty_f,
            min_cond_log_or=empty_f,
            z_iut_signed=empty_f,
            pvalue_iut=empty_f,
            quad_testable=empty_b,
            log_or_interaction=empty_f,
            se_interaction=empty_f,
            pvalue_interaction=empty_f,
        )

    # Extract per-stratum 4-quadrant vectors of length n_pairs via BLAS gemm
    # Index order: 0 -> ++, 1 -> +-, 2 -> -+, 3 -> --
    r_quad_strata: list[tuple[FloatArray, FloatArray, FloatArray, FloatArray]] = []
    nr_quad_strata: list[tuple[FloatArray, FloatArray, FloatArray, FloatArray]] = []
    r_tot_strata: list[FloatArray] = []
    nr_tot_strata: list[FloatArray] = []

    for sl in strata:
        pu_r, qu_r, mu_r = sl.p_r[u_start:u_end, :], sl.q_r[u_start:u_end, :], sl.m_r[u_start:u_end, :]
        pu_nr, qu_nr, mu_nr = sl.p_nr[u_start:u_end, :], sl.q_nr[u_start:u_end, :], sl.m_nr[u_start:u_end, :]

        r_pp = (pu_r @ sl.p_r.T)[row_rel, col_abs]
        r_pm = (pu_r @ sl.q_r.T)[row_rel, col_abs]
        r_mp = (qu_r @ sl.p_r.T)[row_rel, col_abs]
        r_mm = (qu_r @ sl.q_r.T)[row_rel, col_abs]
        r_tot = (mu_r @ sl.m_r.T)[row_rel, col_abs]

        nr_pp = (pu_nr @ sl.p_nr.T)[row_rel, col_abs]
        nr_pm = (pu_nr @ sl.q_nr.T)[row_rel, col_abs]
        nr_mp = (qu_nr @ sl.p_nr.T)[row_rel, col_abs]
        nr_mm = (qu_nr @ sl.q_nr.T)[row_rel, col_abs]
        nr_tot = (mu_nr @ sl.m_nr.T)[row_rel, col_abs]

        r_quad_strata.append((r_pp, r_pm, r_mp, r_mm))
        nr_quad_strata.append((nr_pp, nr_pm, nr_mp, nr_mm))
        r_tot_strata.append(r_tot)
        nr_tot_strata.append(nr_tot)

    # Filter to pairs with sufficient total effective observability
    eff_total: FloatArray = np.sum(np.stack(r_tot_strata, axis=0), axis=0) + np.sum(
        np.stack(nr_tot_strata, axis=0),
        axis=0,
    )
    keep_pair = eff_total >= screen_config.min_effective_samples
    if not bool(np.any(keep_pair)):
        empty_f = np.empty(0, dtype=np.float64)
        empty_i = np.empty(0, dtype=np.int64)
        empty_b = np.empty(0, dtype=bool)
        return _ChunkQuadrantArrays(
            u_local=empty_i,
            v_local=empty_i,
            quadrant_idx=empty_i,
            jaccard=empty_f,
            scope=np.empty(0, dtype=str),
            spatial_weight=empty_f,
            prior_multiplier=empty_f,
            n_r_q=empty_f,
            n_nr_q=empty_f,
            n_eff_total=empty_f,
            raw_donors_enriched=empty_i,
            log_or_quadrant=empty_f,
            se_quadrant=empty_f,
            pvalue_quadrant=empty_f,
            min_cond_log_or=empty_f,
            z_iut_signed=empty_f,
            pvalue_iut=empty_f,
            quad_testable=empty_b,
            log_or_interaction=empty_f,
            se_interaction=empty_f,
            pvalue_interaction=empty_f,
        )

    row_rel_k = row_rel[keep_pair]
    col_abs_k = col_abs[keep_pair]
    u_abs_k = row_rel_k + u_start
    eff_total_k = eff_total[keep_pair]
    jaccard_k = jaccard_block[row_rel_k, col_abs_k]
    scope_k = scope_block[row_rel_k, col_abs_k]

    w1d = overlap_active.spatial_weights_1d
    spatial_w_k = (
        (w1d[u_abs_k] * w1d[col_abs_k])
        if screen_config.use_spatial_fdr_weights
        else np.ones_like(eff_total_k)
    )
    prior_mult_k = (
        w_prior_active.unwrap()[u_abs_k, col_abs_k]
        if isinstance(w_prior_active, Some)
        else np.ones_like(eff_total_k)
    )

    # Raw integer donor support per arm via BLAS dgemm (avoids allocating O(n_pairs * J) slices)
    bu_r, nbu_r = raw_bin_r[u_start:u_end, :], not_bin_r[u_start:u_end, :]
    bu_nr, nbu_nr = raw_bin_nr[u_start:u_end, :], not_bin_nr[u_start:u_end, :]
    raw_r_counts = (
        np.rint((bu_r @ raw_bin_r.T)[row_rel_k, col_abs_k]).astype(np.int64),
        np.rint((bu_r @ not_bin_r.T)[row_rel_k, col_abs_k]).astype(np.int64),
        np.rint((nbu_r @ raw_bin_r.T)[row_rel_k, col_abs_k]).astype(np.int64),
        np.rint((nbu_r @ not_bin_r.T)[row_rel_k, col_abs_k]).astype(np.int64),
    )
    raw_nr_counts = (
        np.rint((bu_nr @ raw_bin_nr.T)[row_rel_k, col_abs_k]).astype(np.int64),
        np.rint((bu_nr @ not_bin_nr.T)[row_rel_k, col_abs_k]).astype(np.int64),
        np.rint((nbu_nr @ raw_bin_nr.T)[row_rel_k, col_abs_k]).astype(np.int64),
        np.rint((nbu_nr @ not_bin_nr.T)[row_rel_k, col_abs_k]).astype(np.int64),
    )

    # Slice per-stratum vectors to keep_pair
    r_q_k = [tuple(arr[keep_pair] for arr in quad) for quad in r_quad_strata]
    nr_q_k = [tuple(arr[keep_pair] for arr in quad) for quad in nr_quad_strata]
    r_tot_k = [arr[keep_pair] for arr in r_tot_strata]
    nr_tot_k = [arr[keep_pair] for arr in nr_tot_strata]

    # Outcome-blind quadrant marginal totals across all strata for testability filtering
    stratum_quad_totals = [
        (
            r_q_k[k][0] + nr_q_k[k][0],
            r_q_k[k][1] + nr_q_k[k][1],
            r_q_k[k][2] + nr_q_k[k][2],
            r_q_k[k][3] + nr_q_k[k][3],
        )
        for k in range(len(r_q_k))
    ]
    stratum_eff_totals = [r_tot_k[k] + nr_tot_k[k] for k in range(len(r_q_k))]
    tot_quad_all = tuple(
        np.sum(np.stack([stratum_quad_totals[k][qi] for k in range(len(r_q_k))], axis=0), axis=0)
        for qi in range(4)
    )
    min_flip_tot = min(1.5, 0.5 * screen_config.min_donors_per_quadrant)

    # Saturated 2x2 epistatic interaction across strata (inverse-variance weighted)
    k_strata = len(r_q_k)
    delta = screen_config.pseudocount / float(max(k_strata, 1))
    inter_num = np.zeros_like(eff_total_k)
    inter_wsum = np.zeros_like(eff_total_k)
    for k in range(k_strata):
        r0, r1, r2, r3 = (x + delta for x in r_q_k[k])
        n0, n1, n2, n3 = (x + delta for x in nr_q_k[k])
        beta_k = np.log((r0 * n1 * n2 * r3) / np.maximum(n0 * r1 * r2 * n3, 1e-300))
        var_k = (
            (1.0 / r0) + (1.0 / n0)
            + (1.0 / r1) + (1.0 / n1)
            + (1.0 / r2) + (1.0 / n2)
            + (1.0 / r3) + (1.0 / n3)
        )
        wk = 1.0 / np.maximum(var_k, 1e-12)
        inter_num = inter_num + wk * beta_k
        inter_wsum = inter_wsum + wk

    log_or_inter = inter_num / np.maximum(inter_wsum, 1e-12)
    se_inter = 1.0 / np.sqrt(np.maximum(inter_wsum, 1e-12))
    pval_inter = np.asarray(2.0 * norm.sf(np.abs(log_or_inter / se_inter)), dtype=np.float64)

    out_u: list[IntArray] = []
    out_v: list[IntArray] = []
    out_q: list[IntArray] = []
    out_jacc: list[FloatArray] = []
    out_scope: list[NDArray[np.str_]] = []
    out_sw: list[FloatArray] = []
    out_pm: list[FloatArray] = []
    out_nr_r: list[FloatArray] = []
    out_nr_nr: list[FloatArray] = []
    out_eff: list[FloatArray] = []
    out_raw_donors: list[IntArray] = []
    out_lor_q: list[FloatArray] = []
    out_se_q: list[FloatArray] = []
    out_p_q: list[FloatArray] = []
    out_min_cond: list[FloatArray] = []
    out_z_iut: list[FloatArray] = []
    out_p_iut: list[FloatArray] = []
    out_testable: list[BoolArray] = []
    out_lor_int: list[FloatArray] = []
    out_se_int: list[FloatArray] = []
    out_p_int: list[FloatArray] = []

    n_kept = int(u_abs_k.shape[0])
    for q_state in ALL_QUADRANT_STATES:
        qi = _quadrant_index(q_state)
        flip_v_state, flip_u_state = adjacent_flips(q_state)
        fv_i = _quadrant_index(flip_v_state)
        fu_i = _quadrant_index(flip_u_state)

        a_q = [r_q_k[k][qi] for k in range(k_strata)]
        b_q = [nr_q_k[k][qi] for k in range(k_strata)]
        c_not_q = [np.maximum(r_tot_k[k] - r_q_k[k][qi], 0.0) for k in range(k_strata)]
        d_not_q = [np.maximum(nr_tot_k[k] - nr_q_k[k][qi], 0.0) for k in range(k_strata)]

        lor_q, se_q, z_q = _cmh_2x2_from_strata(
            a_list=a_q,
            b_list=b_q,
            c_list=c_not_q,
            d_list=d_not_q,
            pseudocount=screen_config.pseudocount,
        )
        p_q = np.asarray(2.0 * norm.sf(np.abs(z_q)), dtype=np.float64)

        # Conditional flip 1: target q vs flip_v (necessity of v given u)
        lor_fv, _, z_fv = _cmh_2x2_from_strata(
            a_list=a_q,
            b_list=b_q,
            c_list=[r_q_k[k][fv_i] for k in range(k_strata)],
            d_list=[nr_q_k[k][fv_i] for k in range(k_strata)],
            pseudocount=screen_config.pseudocount,
        )
        # Conditional flip 2: target q vs flip_u (necessity of u given v)
        lor_fu, _, z_fu = _cmh_2x2_from_strata(
            a_list=a_q,
            b_list=b_q,
            c_list=[r_q_k[k][fu_i] for k in range(k_strata)],
            d_list=[nr_q_k[k][fu_i] for k in range(k_strata)],
            pseudocount=screen_config.pseudocount,
        )

        # Berger's Intersection-Union Test aligned with overall quadrant enrichment direction
        is_resp_dir = np.where(np.abs(z_q) > 1e-12, z_q >= 0.0, lor_q >= 0.0)
        sign_dir = np.where(is_resp_dir, 1.0, -1.0)
        z_iut_signed = np.minimum(sign_dir * z_fv, sign_dir * z_fu)
        p_iut = np.clip(
            np.asarray(2.0 * norm.sf(z_iut_signed), dtype=np.float64),
            0.0,
            1.0,
        )
        min_cond_lor = np.where(
            is_resp_dir,
            np.minimum(lor_fv, lor_fu),
            np.maximum(lor_fv, lor_fu),
        )

        # Outcome-blind Tarone-style marginal testability mask (strictly symmetric under u <-> v)
        tq_vec = tot_quad_all[qi]
        testable_q = (
            (tq_vec >= screen_config.min_donors_per_quadrant)
            & ((eff_total_k - tq_vec) >= screen_config.min_donors_per_quadrant)
            & (tot_quad_all[fv_i] >= min_flip_tot)
            & (tot_quad_all[fu_i] >= min_flip_tot)
            & (prior_mult_k > 0.0)
        )
        if tarone_z_floor > 0.0:
            z_max_iut = compute_tarone_max_iut_z_strata(
                stratum_quad_totals=stratum_quad_totals,
                stratum_eff_totals=stratum_eff_totals,
                stratum_arm_counts=stratum_arm_counts,
                qi=qi,
                fv_i=fv_i,
                fu_i=fu_i,
            )
            testable_q = testable_q & (z_max_iut >= (tarone_z_floor - 1e-6))

        sum_r_q: FloatArray = np.sum(np.stack(a_q, axis=0), axis=0)
        sum_nr_q: FloatArray = np.sum(np.stack(b_q, axis=0), axis=0)
        raw_don_enr = np.where(is_resp_dir, raw_r_counts[qi], raw_nr_counts[qi])

        q_parity = 1.0 if q_state in (QuadrantState.BOTH_PRESENT, QuadrantState.BOTH_MISSING) else -1.0

        out_u.append(u_abs_k)
        out_v.append(col_abs_k)
        out_q.append(np.full(n_kept, qi, dtype=np.int64))
        out_jacc.append(jaccard_k)
        out_scope.append(scope_k)
        out_sw.append(spatial_w_k)
        out_pm.append(prior_mult_k)
        out_nr_r.append(sum_r_q)
        out_nr_nr.append(sum_nr_q)
        out_eff.append(eff_total_k)
        out_raw_donors.append(raw_don_enr)
        out_lor_q.append(lor_q)
        out_se_q.append(se_q)
        out_p_q.append(p_q)
        out_min_cond.append(min_cond_lor)
        out_z_iut.append(z_iut_signed)
        out_p_iut.append(p_iut)
        out_testable.append(testable_q)
        out_lor_int.append(q_parity * log_or_inter)
        out_se_int.append(se_inter)
        out_p_int.append(pval_inter)

    return _ChunkQuadrantArrays(
        u_local=np.concatenate(out_u),
        v_local=np.concatenate(out_v),
        quadrant_idx=np.concatenate(out_q),
        jaccard=np.concatenate(out_jacc),
        scope=np.concatenate(out_scope),
        spatial_weight=np.concatenate(out_sw),
        prior_multiplier=np.concatenate(out_pm),
        n_r_q=np.concatenate(out_nr_r),
        n_nr_q=np.concatenate(out_nr_nr),
        n_eff_total=np.concatenate(out_eff),
        raw_donors_enriched=np.concatenate(out_raw_donors),
        log_or_quadrant=np.concatenate(out_lor_q),
        se_quadrant=np.concatenate(out_se_q),
        pvalue_quadrant=np.concatenate(out_p_q),
        min_cond_log_or=np.concatenate(out_min_cond),
        z_iut_signed=np.concatenate(out_z_iut),
        pvalue_iut=np.concatenate(out_p_iut),
        quad_testable=np.concatenate(out_testable),
        log_or_interaction=np.concatenate(out_lor_int),
        se_interaction=np.concatenate(out_se_int),
        pvalue_interaction=np.concatenate(out_p_int),
    )


def _compute_stratified_quadrant_permutations(
    u_local: IntArray,
    v_local: IntArray,
    q_idx: IntArray,
    presence: PresenceRepresentation,
    response_binary: FloatArray,
    informative_groups: list[IntArray],
    n_permutations: int,
    seed: int,
) -> FloatArray:
    """Computes exact stratified Rao score permutation p-values for retained candidates."""
    n_cand = int(u_local.shape[0])
    if n_cand == 0 or n_permutations <= 0:
        return np.ones(n_cand, dtype=np.float64)

    n_samples = int(response_binary.shape[0])
    p_mat = presence.presence_prob
    q_mat = presence.absence_prob

    u_is_pos = (q_idx == 0) | (q_idx == 1)
    v_is_pos = (q_idx == 0) | (q_idx == 2)
    s_u = np.where(np.expand_dims(u_is_pos, axis=1), p_mat[u_local, :], q_mat[u_local, :])
    s_v = np.where(np.expand_dims(v_is_pos, axis=1), p_mat[v_local, :], q_mat[v_local, :])
    s_cand = s_u * s_v

    rng = np.random.default_rng(seed)
    y_centered = np.zeros(n_samples, dtype=np.float64)
    y_perm_centered = np.zeros((n_samples, n_permutations), dtype=np.float64)

    for grp in informative_groups:
        y_g = response_binary[grp]
        mean_g = float(np.mean(y_g))
        y_centered[grp] = y_g - mean_g
        for b in range(n_permutations):
            y_perm_centered[grp, b] = rng.permutation(y_g) - mean_g

    t_obs = np.abs(s_cand @ y_centered)
    t_perm = np.abs(s_cand @ y_perm_centered)
    exceedances = np.sum(t_perm >= (np.expand_dims(t_obs, axis=1) - 1e-12), axis=1)
    return np.asarray(
        (1.0 + exceedances.astype(np.float64)) / (1.0 + float(n_permutations)),
        dtype=np.float64,
    )


def _empty_quadrant_df() -> pl.DataFrame:
    return pl.DataFrame(
        schema={
            "nhood_u": pl.Int64,
            "nhood_v": pl.Int64,
            "quadrant": pl.String,
            "direction": pl.String,
            "pair_scope": pl.String,
            "cell_type_u": pl.String,
            "cell_type_v": pl.String,
            "jaccard_overlap": pl.Float64,
            "spatial_weight": pl.Float64,
            "eff_donors_responder": pl.Float64,
            "eff_donors_non_responder": pl.Float64,
            "eff_donors_total": pl.Float64,
            "raw_donors_enriched_quadrant": pl.Int64,
            "log_or_quadrant": pl.Float64,
            "se_quadrant": pl.Float64,
            "pvalue_quadrant": pl.Float64,
            "fdr_quadrant": pl.Float64,
            "min_cond_log_or": pl.Float64,
            "pvalue_iut": pl.Float64,
            "fdr_iut": pl.Float64,
            "log_or_interaction": pl.Float64,
            "se_interaction": pl.Float64,
            "pvalue_interaction": pl.Float64,
            "fdr_interaction": pl.Float64,
            "pvalue_perm": pl.Float64,
            "fdr_perm": pl.Float64,
            "passes_donor_gate": pl.Boolean,
            "is_significant": pl.Boolean,
        }
    )


def screen_quadrant_pairs(
    presence: PresenceRepresentation,
    overlap: OverlapGraph,
    response_binary: FloatArray,
    strata_codes: Maybe[IntArray] = Nothing,
    overlap_config: OverlapConfig = OverlapConfig(),
    screen_config: ScreenConfig = ScreenConfig(),
    prior_pair_weights: Maybe[FloatArray] = Nothing,
) -> Result[tuple[pl.DataFrame, int], str]:
    """Screens all valid non-overlapping pairs across the 4 presence/absence quadrants."""
    strata, informative_groups = _build_strata_slices(presence, response_binary, strata_codes)
    if not strata:
        return Failure("No informative strata contain both responder and non-responder donors")

    overlap_active = subset_overlap_graph(overlap, presence.active_nhood_indices)
    v_act = overlap_active.n_nhoods
    if v_act < 2:
        return Success((_empty_quadrant_df(), 0))

    prior_res = _normalize_active_prior_matrix(
        prior_pair_weights=prior_pair_weights,
        active_nhood_indices=presence.active_nhood_indices,
        overlap_active=overlap_active,
        overlap_config=overlap_config,
    )
    match prior_res:
        case Failure(err):
            return Failure(err)
        case Success(w_prior_active):
            pass

    has_external_prior = isinstance(w_prior_active, Some)
    tarone_z_floor = (
        tarone_z_threshold_from_p_ceiling(screen_config.tarone_p_ceiling)
        if screen_config.apply_tarone_iut_pruning
        else 0.0
    )
    stratum_arm_counts = [(int(sl.p_r.shape[1]), int(sl.p_nr.shape[1])) for sl in strata]

    informative_idx = np.concatenate(informative_groups)
    r_inf_mask = response_binary[informative_idx] > 0.5
    r_inf_idx = informative_idx[r_inf_mask]
    nr_inf_idx = informative_idx[~r_inf_mask]
    raw_bin_r = presence.raw_binary_presence[:, r_inf_idx].astype(np.float64)
    not_bin_r = 1.0 - raw_bin_r
    raw_bin_nr = presence.raw_binary_presence[:, nr_inf_idx].astype(np.float64)
    not_bin_nr = 1.0 - raw_bin_nr

    chunk_results = [
        _evaluate_quadrant_chunk(
            u_start=start,
            u_end=min(start + screen_config.chunk_size, v_act),
            strata=strata,
            stratum_arm_counts=stratum_arm_counts,
            tarone_z_floor=tarone_z_floor,
            w_prior_active=w_prior_active,
            raw_bin_r=raw_bin_r,
            not_bin_r=not_bin_r,
            raw_bin_nr=raw_bin_nr,
            not_bin_nr=not_bin_nr,
            overlap_active=overlap_active,
            overlap_config=overlap_config,
            screen_config=screen_config,
        )
        for start in range(0, v_act, screen_config.chunk_size)
    ]

    non_empty = [c for c in chunk_results if int(c.u_local.shape[0]) > 0]
    if not non_empty:
        return Success((_empty_quadrant_df(), 0))

    u_loc = np.concatenate([c.u_local for c in non_empty])
    v_loc = np.concatenate([c.v_local for c in non_empty])
    q_idx = np.concatenate([c.quadrant_idx for c in non_empty])
    jacc = np.concatenate([c.jaccard for c in non_empty])
    scope = np.concatenate([c.scope for c in non_empty])
    sw = np.concatenate([c.spatial_weight for c in non_empty])
    pm = np.concatenate([c.prior_multiplier for c in non_empty])
    nr_r = np.concatenate([c.n_r_q for c in non_empty])
    nr_nr = np.concatenate([c.n_nr_q for c in non_empty])
    n_eff = np.concatenate([c.n_eff_total for c in non_empty])
    raw_don = np.concatenate([c.raw_donors_enriched for c in non_empty])
    lor_q = np.concatenate([c.log_or_quadrant for c in non_empty])
    se_q = np.concatenate([c.se_quadrant for c in non_empty])
    p_q = np.concatenate([c.pvalue_quadrant for c in non_empty])
    min_cond = np.concatenate([c.min_cond_log_or for c in non_empty])
    z_iut = np.concatenate([c.z_iut_signed for c in non_empty])
    p_iut = np.concatenate([c.pvalue_iut for c in non_empty])
    testable = np.concatenate([c.quad_testable for c in non_empty])
    lor_int = np.concatenate([c.log_or_interaction for c in non_empty])
    se_int = np.concatenate([c.se_interaction for c in non_empty])
    p_int = np.concatenate([c.pvalue_interaction for c in non_empty])

    n_evaluated_pairs = int(np.sum(pm > 0.0)) // 4
    sw_bh = sw * pm

    # Outcome-blind Tarone-filtered 2D spatial + prior-weighted BH for quadrant & epistatic contrasts
    fdr_q = np.ones_like(p_q)
    fdr_iut = np.ones_like(p_iut)
    fdr_int = np.ones_like(p_int)
    pos_prior = pm > 0.0
    if bool(np.any(pos_prior)):
        fdr_int[pos_prior] = compute_weighted_bh_fdr(p_int[pos_prior], sw_bh[pos_prior])

    if bool(np.any(testable)):
        fdr_q[testable] = compute_weighted_bh_fdr(p_q[testable], sw_bh[testable])
        fdr_iut[testable] = _compute_decoupled_quadrant_iut_fdr(
            z_iut_testable=z_iut[testable],
            sw_testable=sw[testable],
            prior_mult_testable=pm[testable],
            has_external_prior=has_external_prior,
            w_prior_active=w_prior_active,
            tarone_z_floor=tarone_z_floor,
            presence=presence,
            strata=strata,
            informative_groups=informative_groups,
            response_binary=response_binary,
            overlap_active=overlap_active,
            overlap_config=overlap_config,
            screen_config=screen_config,
        )

    is_resp = np.where(np.abs(lor_q) > 1e-12, lor_q >= 0.0, nr_r >= nr_nr)
    eff_enriched = np.where(is_resp, nr_r, nr_nr)
    passes_donor = (eff_enriched >= screen_config.min_donors_per_quadrant) & (
        raw_don.astype(np.float64) >= np.floor(screen_config.min_donors_per_quadrant)
    )
    is_sig = (
        testable
        & (fdr_q < screen_config.fdr_threshold)
        & (fdr_iut < screen_config.fdr_threshold)
        & passes_donor
    )

    # Retain all significant hits plus top candidates ranked by (testable, p_iut) up to max_retained_candidates
    if int(u_loc.shape[0]) > screen_config.max_retained_candidates:
        rank_score = np.where(testable, p_iut, 2.0 + p_iut)
        rank_order = np.argsort(rank_score, kind="mergesort")
        top_mask = np.zeros(u_loc.shape[0], dtype=bool)
        top_mask[rank_order[: screen_config.max_retained_candidates]] = True
        retain_mask = is_sig | top_mask
    else:
        retain_mask = np.ones(u_loc.shape[0], dtype=bool)

    u_ret = u_loc[retain_mask]
    v_ret = v_loc[retain_mask]
    q_ret = q_idx[retain_mask]
    sw_ret = sw[retain_mask]

    p_q_masked = np.where(testable, p_q, 1.0)
    p_iut_masked = np.where(testable, p_iut, 1.0)

    if screen_config.n_permutations > 0:
        p_perm = _compute_stratified_quadrant_permutations(
            u_local=u_ret,
            v_local=v_ret,
            q_idx=q_ret,
            presence=presence,
            response_binary=response_binary,
            informative_groups=informative_groups,
            n_permutations=screen_config.n_permutations,
            seed=screen_config.permutation_seed,
        )
        fdr_perm = compute_weighted_bh_fdr(p_perm, sw_ret)
    else:
        p_perm = p_q_masked[retain_mask]
        fdr_perm = fdr_q[retain_mask]

    act_idx = presence.active_nhood_indices
    ct_arr = np.asarray(overlap_active.dominant_cell_types, dtype=str)
    q_labels = np.asarray([s.value for s in ALL_QUADRANT_STATES], dtype=str)
    dir_labels = np.where(
        is_resp[retain_mask],
        EnrichmentDirection.RESPONDER.value,
        EnrichmentDirection.NON_RESPONDER.value,
    )

    df = (
        pl.DataFrame(
            {
                "nhood_u": act_idx[u_ret],
                "nhood_v": act_idx[v_ret],
                "quadrant": q_labels[q_ret],
                "direction": dir_labels,
                "pair_scope": scope[retain_mask],
                "cell_type_u": ct_arr[u_ret],
                "cell_type_v": ct_arr[v_ret],
                "jaccard_overlap": jacc[retain_mask],
                "spatial_weight": sw_ret,
                "eff_donors_responder": nr_r[retain_mask],
                "eff_donors_non_responder": nr_nr[retain_mask],
                "eff_donors_total": n_eff[retain_mask],
                "raw_donors_enriched_quadrant": raw_don[retain_mask],
                "log_or_quadrant": lor_q[retain_mask],
                "se_quadrant": se_q[retain_mask],
                "pvalue_quadrant": p_q_masked[retain_mask],
                "fdr_quadrant": fdr_q[retain_mask],
                "min_cond_log_or": min_cond[retain_mask],
                "pvalue_iut": p_iut_masked[retain_mask],
                "fdr_iut": fdr_iut[retain_mask],
                "log_or_interaction": lor_int[retain_mask],
                "se_interaction": se_int[retain_mask],
                "pvalue_interaction": p_int[retain_mask],
                "fdr_interaction": fdr_int[retain_mask],
                "pvalue_perm": p_perm,
                "fdr_perm": fdr_perm,
                "passes_donor_gate": passes_donor[retain_mask],
                "is_significant": is_sig[retain_mask],
            }
        )
        .lazy()
        .sort(["is_significant", "pvalue_iut", "pvalue_quadrant"], descending=[True, False, False])
        .collect()
    )
    return Success((df, n_evaluated_pairs))


def _stratum_center_and_standardize(
    features: FloatArray,
    response_binary: FloatArray,
    informative_groups: list[IntArray],
) -> tuple[FloatArray, FloatArray, int]:
    """Centers features and response within each informative stratum and scales to unit L2 norm."""
    used_idx = np.concatenate(informative_groups)
    n_used = int(used_idx.shape[0])
    n_strata = len(informative_groups)
    df_resid = max(n_used - n_strata, 3)

    x_centered = np.zeros((features.shape[0], n_used), dtype=np.float64)
    y_centered = np.zeros(n_used, dtype=np.float64)

    offset = 0
    for grp in informative_groups:
        g_len = int(grp.shape[0])
        sl = slice(offset, offset + g_len)
        xg = features[:, grp]
        yg = response_binary[grp]
        x_centered[:, sl] = xg - np.mean(xg, axis=1, keepdims=True)
        y_centered[sl] = yg - float(np.mean(yg))
        offset += g_len

    x_norms = np.linalg.norm(x_centered, axis=1, keepdims=True)
    x_unit = np.where(x_norms > 1e-12, x_centered / np.maximum(x_norms, 1e-12), 0.0)
    y_norm = float(np.linalg.norm(y_centered))
    y_unit = (y_centered / y_norm) if y_norm > 1e-12 else np.zeros_like(y_centered)
    return x_unit, y_unit, df_resid


def _apply_prior_weight_to_t(t_vals: FloatArray, prior_mult: FloatArray, df_t: int) -> FloatArray:
    """Monotonically maps Student-t statistics via Genovese-Roeder-Wasserman p / w_prior scaling."""
    if int(t_vals.shape[0]) == 0:
        return t_vals
    pos = t_vals > 0.0
    if not bool(np.any(pos)):
        return t_vals
    out = np.array(t_vals, dtype=np.float64, copy=True)
    p_two = np.clip(2.0 * student_t.sf(t_vals[pos], df=df_t), 1e-300, 1.999999)
    p_adj = np.clip(p_two / np.maximum(prior_mult[pos], 1e-12), 1e-300, 1.999999)
    out[pos] = np.asarray(student_t.isf(0.5 * p_adj, df=df_t), dtype=np.float64)
    return out


def _compute_decoupled_balance_iut_fdr(
    t_iut_testable: FloatArray,
    sw_testable: FloatArray,
    prior_mult_testable: FloatArray,
    has_external_prior: bool,
    w_prior_active: Maybe[FloatArray],
    z_cpm: FloatArray,
    cpm_centered: FloatArray,
    r_cpm: FloatArray,
    cpm_dot_y: FloatArray,
    y_unit: FloatArray,
    df_two: int,
    pos_support_ok: BoolArray,
    informative_groups: list[IntArray],
    overlap_active: OverlapGraph,
    overlap_config: OverlapConfig,
    screen_config: ScreenConfig,
) -> FloatArray:
    """Computes Symmetrized Decoupled Partner-Permutation (DMPP) Knockoff+ q-values for log-ratio IUT."""
    n_test = int(t_iut_testable.shape[0])
    if n_test == 0:
        return np.empty(0, dtype=np.float64)

    v_act = overlap_active.n_nhoods
    b_draws = _resolve_dmpp_draw_count(
        v_act,
        screen_config.fdr_threshold,
        screen_config.dmpp_null_draws,
    )
    grid_edges = np.linspace(1.5, 16.0, _DMPP_GRID_BINS + 1, dtype=np.float64)
    null_hist = np.zeros(_DMPP_GRID_BINS, dtype=np.float64)

    perm_block = max(screen_config.chunk_size, 1024)
    block_ranges = [
        (start, min(start + perm_block, v_act))
        for start in range(0, v_act, perm_block)
    ]
    w_prior_mat: FloatArray = (
        w_prior_active.unwrap()
        if isinstance(w_prior_active, Some)
        else np.empty((0, 0), dtype=np.float64)
    )
    support_2d_cols = np.expand_dims(pos_support_ok, axis=0)
    sym_masks = [
        (
            _build_symmetric_offdiag_mask(overlap_active, s, e, overlap_config)
            & np.expand_dims(pos_support_ok[s:e], axis=1)
            & support_2d_cols
            & (w_prior_mat[s:e, :] > 0.0)
        )
        if has_external_prior
        else (
            _build_symmetric_offdiag_mask(overlap_active, s, e, overlap_config)
            & np.expand_dims(pos_support_ok[s:e], axis=1)
            & support_2d_cols
        )
        for s, e in block_ranges
    ]
    if not any(bool(np.any(m)) for m in sym_masks):
        return np.ones(n_test, dtype=np.float64)

    w1d = (
        overlap_active.spatial_weights_1d.astype(np.float64)
        if screen_config.use_spatial_fdr_weights
        else np.ones(v_act, dtype=np.float64)
    )

    z_f32 = z_cpm.astype(np.float32)
    n_used = int(z_cpm.shape[1])
    group_lengths = [int(g.shape[0]) for g in informative_groups]
    rng = np.random.default_rng(screen_config.permutation_seed + 1009)
    df_two_f32 = np.float32(float(df_two))
    ru_all = r_cpm.astype(np.float32)
    cdy_all = cpm_dot_y.astype(np.float32)
    w_max = float(np.max(w_prior_mat)) if (has_external_prior and w_prior_mat.size > 0) else 1.0
    raw_min_t = (
        float(student_t.isf(min(0.499, float(student_t.sf(1.5, df=df_two)) * max(w_max, 1.0)), df=df_two))
        if (has_external_prior and w_max > 1.0)
        else 1.5
    )
    t_floor_f32 = np.float32(max(0.25, raw_min_t - 1e-4))

    for _ in range(b_draws):
        perm_used = np.arange(n_used, dtype=np.int64)
        offset = 0
        for g_len in group_lengths:
            local_slice = np.arange(offset, offset + g_len, dtype=np.int64)
            perm_used[offset : offset + g_len] = rng.permutation(local_slice)
            offset += g_len

        z_perm_f32 = z_f32[:, perm_used]
        r_perm = (z_cpm[:, perm_used] @ y_unit).astype(np.float32)
        cdy_perm = (cpm_centered[:, perm_used] @ y_unit).astype(np.float32)

        for b_idx, (start, end) in enumerate(block_ranges):
            vs = sym_masks[b_idx]
            if not bool(np.any(vs)):
                continue

            rho_blk = np.clip(z_f32[start:end, :] @ z_perm_f32.T, np.float32(-0.9999), np.float32(0.9999))
            one_m_rho2 = np.maximum(np.float32(1.0) - rho_blk * rho_blk, np.float32(1e-8))
            ru_b = np.expand_dims(ru_all[start:end], axis=1)
            rv_b = np.expand_dims(r_perm, axis=0)

            bu = (ru_b - rho_blk * rv_b) / one_m_rho2
            bv = (rv_b - rho_blk * ru_b) / one_m_rho2
            r2 = np.clip(ru_b * bu + rv_b * bv, np.float32(0.0), np.float32(0.999999))
            se = np.sqrt(np.maximum(np.float32(1.0) - r2, np.float32(1e-12)) / (df_two_f32 * one_m_rho2))

            sgn = np.where(
                np.expand_dims(cdy_all[start:end], axis=1) >= np.expand_dims(cdy_perm, axis=0),
                np.float32(1.0),
                np.float32(-1.0),
            )
            t_min_s = np.minimum(sgn * (bu / se), sgn * ((-bv) / se))
            hi_mask = vs & (t_min_s >= t_floor_f32)
            ii, jj = np.nonzero(hi_mask)
            if int(ii.shape[0]) > 0:
                t_vals_raw = t_min_s[ii, jj].astype(np.float64) + 1e-5
                t_vals = (
                    _apply_prior_weight_to_t(t_vals_raw, w_prior_mat[ii + start, jj], df_two)
                    if has_external_prior
                    else t_vals_raw
                )
                ok_v = t_vals >= 1.5
                if bool(np.any(ok_v)):
                    bins = np.clip(
                        np.searchsorted(grid_edges[:-1], t_vals[ok_v], side="right") - 1,
                        0,
                        _DMPP_GRID_BINS - 1,
                    )
                    w_cand = 0.5 * w1d[ii[ok_v] + start] * w1d[jj[ok_v]]
                    null_hist = null_hist + np.bincount(bins, weights=w_cand, minlength=_DMPP_GRID_BINS)

    obs_t = (
        _apply_prior_weight_to_t(t_iut_testable, prior_mult_testable, df_two)
        if has_external_prior
        else t_iut_testable
    )
    return _stepup_knockoff_plus_fdr(
        stat_obs=obs_t,
        weights_obs=sw_testable,
        null_hist_weights=null_hist,
        grid_edges=grid_edges,
        n_null_draws=b_draws,
    )


def _empty_continuous_df() -> pl.DataFrame:
    return pl.DataFrame(
        schema={
            "nhood_u": pl.Int64,
            "nhood_v": pl.Int64,
            "axis": pl.String,
            "direction": pl.String,
            "pair_scope": pl.String,
            "cell_type_u": pl.String,
            "cell_type_v": pl.String,
            "jaccard_overlap": pl.Float64,
            "spatial_weight": pl.Float64,
            "effect_size": pl.Float64,
            "pvalue_axis": pl.Float64,
            "fdr_axis": pl.Float64,
            "pvalue_iut": pl.Float64,
            "fdr_iut": pl.Float64,
            "is_significant": pl.Boolean,
        }
    )


def _screen_continuous_partitioned_by_observability(
    presence: PresenceRepresentation,
    overlap: OverlapGraph,
    response_binary: FloatArray,
    strata_codes: Maybe[IntArray],
    informative_groups: list[IntArray],
    stratum_obs: BoolArray,
    axes: tuple[ContinuousAxis, ...],
    overlap_config: OverlapConfig,
    screen_config: ScreenConfig,
    prior_pair_weights: Maybe[FloatArray],
) -> Result[pl.DataFrame, str]:
    """Evaluates continuous pairs strictly on strata where both endpoints `u` and `v` are observable."""
    v_act = int(presence.active_nhood_indices.shape[0])
    act_indices = presence.active_nhood_indices.tolist()
    nhood_to_sig: dict[int, tuple[bool, ...]] = {
        int(act_indices[i]): tuple(bool(x) for x in stratum_obs[i, :].tolist())
        for i in range(v_act)
    }
    unique_sigs = sorted(set(nhood_to_sig.values()))
    joint_sigs = sorted(
        {
            tuple(a and b for a, b in zip(sa, sb))
            for sa in unique_sigs
            for sb in unique_sigs
            if any(a and b for a, b in zip(sa, sb))
        }
    )
    if not joint_sigs:
        return Success(_empty_continuous_df())

    uncapped_cfg = ScreenConfig(
        chunk_size=screen_config.chunk_size,
        fdr_threshold=screen_config.fdr_threshold,
        min_donors_per_quadrant=screen_config.min_donors_per_quadrant,
        min_effective_samples=screen_config.min_effective_samples,
        pseudocount=screen_config.pseudocount,
        use_spatial_fdr_weights=screen_config.use_spatial_fdr_weights,
        n_permutations=screen_config.n_permutations,
        permutation_seed=screen_config.permutation_seed,
        max_retained_candidates=max(screen_config.max_retained_candidates, v_act * v_act + 100),
        apply_tarone_iut_pruning=screen_config.apply_tarone_iut_pruning,
        tarone_p_ceiling=screen_config.tarone_p_ceiling,
        dmpp_null_draws=screen_config.dmpp_null_draws,
    )

    domain_frames: list[pl.DataFrame] = []
    for j_sig in joint_sigs:
        keep_rows = np.array(
            [
                all((not req) or has_g for req, has_g in zip(j_sig, stratum_obs[i, :].tolist()))
                for i in range(v_act)
            ],
            dtype=bool,
        )
        n_sub_rows = int(np.sum(keep_rows))
        if n_sub_rows < 2:
            continue

        keep_cols = np.sort(
            np.concatenate([grp for grp, ok in zip(informative_groups, j_sig) if ok])
        )
        n_sub_cols = int(keep_cols.shape[0])
        if n_sub_cols < 6:
            continue

        sub_pres = PresenceRepresentation(
            active_nhood_indices=presence.active_nhood_indices[keep_rows],
            presence_prob=presence.presence_prob[keep_rows, :][:, keep_cols],
            absence_prob=presence.absence_prob[keep_rows, :][:, keep_cols],
            observability=np.ones((n_sub_rows, n_sub_cols), dtype=np.float64),
            raw_binary_presence=presence.raw_binary_presence[keep_rows, :][:, keep_cols],
            arcsin_sqrt_prop=presence.arcsin_sqrt_prop[keep_rows, :][:, keep_cols],
            log2_cpm=presence.log2_cpm[keep_rows, :][:, keep_cols],
            effective_library_sizes=presence.effective_library_sizes[keep_cols],
        )
        sub_strata: Maybe[IntArray] = (
            Some(strata_codes.unwrap()[keep_cols])
            if isinstance(strata_codes, Some)
            else Nothing
        )
        dom_res = screen_continuous_pairs(
            presence=sub_pres,
            overlap=overlap,
            response_binary=response_binary[keep_cols],
            strata_codes=sub_strata,
            axes=axes,
            overlap_config=overlap_config,
            screen_config=uncapped_cfg,
            prior_pair_weights=prior_pair_weights,
        )
        match dom_res:
            case Failure(err):
                return Failure(err)
            case Success(dom_df):
                if not dom_df.is_empty():
                    dom_ax_all = dom_df["axis"].to_numpy()
                    dom_p_all = dom_df["pvalue_iut"].to_numpy()
                    dom_q_all = dom_df["fdr_iut"].to_numpy()
                    dom_p_cal_all = np.array(dom_p_all, dtype=np.float64, copy=True)
                    for ax_e in axes:
                        ax_mask_all = dom_ax_all == ax_e.value
                        dom_family_m = max(int(np.sum(ax_mask_all)), 1)
                        m_d = ax_mask_all & (
                            (dom_p_all < 1.0 - 1e-15) | (dom_q_all < 1.0 - 1e-15)
                        )
                        dom_m_cnt = int(np.sum(m_d))
                        if dom_m_cnt > 0:
                            p_sub = dom_p_all[m_d]
                            q_sub = dom_q_all[m_d]
                            ord_d = np.lexsort((p_sub, q_sub))
                            ranks_d = np.empty(dom_m_cnt, dtype=np.float64)
                            ranks_d[ord_d] = np.arange(1, dom_m_cnt + 1, dtype=np.float64)
                            dom_p_cal_all[m_d] = np.minimum(
                                p_sub,
                                np.clip(q_sub * (ranks_d / float(dom_family_m)), 0.0, 1.0),
                            )

                    dom_with_cal = dom_df.with_columns(
                        pl.Series("_p_cal", dom_p_cal_all, dtype=pl.Float64)
                    )
                    u_arr = dom_with_cal["nhood_u"].to_list()
                    v_arr = dom_with_cal["nhood_v"].to_list()
                    exact_mask = [
                        tuple(
                            a and b
                            for a, b in zip(nhood_to_sig[int(uu)], nhood_to_sig[int(vv)])
                        )
                        == j_sig
                        for uu, vv in zip(u_arr, v_arr)
                    ]
                    filt_dom = dom_with_cal.filter(pl.Series(exact_mask, dtype=pl.Boolean))
                    if not filt_dom.is_empty():
                        domain_frames.append(filt_dom)

    if not domain_frames:
        return Success(_empty_continuous_df())

    merged = pl.concat(domain_frames, how="vertical")
    ax_col = merged["axis"].to_numpy()
    p_ax = merged["pvalue_axis"].to_numpy()
    p_iut = merged["pvalue_iut"].to_numpy()
    p_cal = merged["_p_cal"].to_numpy()
    dom_fdr_iut = merged["fdr_iut"].to_numpy()
    sw_col = merged["spatial_weight"].to_numpy()

    pm_col = np.ones_like(sw_col, dtype=np.float64)
    if isinstance(prior_pair_weights, Some):
        p_mat_raw = np.asarray(prior_pair_weights.unwrap(), dtype=np.float64)
        u_glob = merged["nhood_u"].to_numpy()
        v_glob = merged["nhood_v"].to_numpy()
        n_g = int(p_mat_raw.shape[0])
        in_b = (u_glob >= 0) & (u_glob < n_g) & (v_glob >= 0) & (v_glob < n_g)
        if bool(np.any(in_b)):
            pm_col[in_b] = np.maximum(p_mat_raw[u_glob[in_b], v_glob[in_b]], 0.0)
    sw_bh_col = sw_col * pm_col

    # Recompute global family-wide weighted BH across observability domains using DMPP-calibrated _p_cal
    fdr_ax_glob = np.ones_like(p_ax)
    fdr_iut_glob = np.ones_like(p_iut)
    for ax_enum in axes:
        m_ax = (ax_col == ax_enum.value) & (
            (dom_fdr_iut < 1.0 - 1e-15) | (p_iut < 1.0 - 1e-15)
        )
        if bool(np.any(m_ax)):
            fdr_ax_glob[m_ax] = compute_weighted_bh_fdr(p_ax[m_ax], sw_bh_col[m_ax])
            fdr_iut_glob[m_ax] = compute_weighted_bh_fdr(p_cal[m_ax], sw_bh_col[m_ax])

    sig_glob = (fdr_ax_glob < screen_config.fdr_threshold) & (
        fdr_iut_glob < screen_config.fdr_threshold
    )
    rebuilt = (
        merged.drop("_p_cal")
        .lazy()
        .with_columns(
            pl.Series("fdr_axis", fdr_ax_glob, dtype=pl.Float64),
            pl.Series("fdr_iut", fdr_iut_glob, dtype=pl.Float64),
            pl.Series("is_significant", sig_glob, dtype=pl.Boolean),
        )
        .sort(["is_significant", "pvalue_iut", "pvalue_axis"], descending=[True, False, False])
        .collect()
    )
    if rebuilt.height > screen_config.max_retained_candidates:
        sig_part = rebuilt.filter(pl.col("is_significant"))
        nonsig_part = rebuilt.filter(~pl.col("is_significant")).head(
            max(screen_config.max_retained_candidates - sig_part.height, 0)
        )
        return Success(pl.concat([sig_part, nonsig_part], how="vertical"))
    return Success(rebuilt)


def screen_continuous_pairs(
    presence: PresenceRepresentation,
    overlap: OverlapGraph,
    response_binary: FloatArray,
    strata_codes: Maybe[IntArray] = Nothing,
    axes: tuple[ContinuousAxis, ...] = (
        ContinuousAxis.LOG_RATIO_BALANCE,
        ContinuousAxis.CO_ABUNDANCE_PRODUCT,
    ),
    overlap_config: OverlapConfig = OverlapConfig(),
    screen_config: ScreenConfig = ScreenConfig(),
    prior_pair_weights: Maybe[FloatArray] = Nothing,
) -> Result[pl.DataFrame, str]:
    """Screens continuous Aitchison log-ratio balances (+-/-+) and co-abundance products (++/--)."""
    if not axes:
        return Success(_empty_continuous_df())

    strata, informative_groups = _build_strata_slices(presence, response_binary, strata_codes)
    if not strata:
        return Failure("No informative strata contain both responder and non-responder donors")

    if int(presence.active_nhood_indices.shape[0]) >= 2:
        all_inf_cols = np.sort(np.concatenate(informative_groups))
        if bool(np.any(presence.observability[:, all_inf_cols] <= 1e-6)):
            obs_blocks: list[IntArray] = []
            for grp in informative_groups:
                g_obs = presence.observability[:, grp] > 1e-6
                if bool(np.all(g_obs)):
                    obs_blocks.append(grp)
                else:
                    col_sigs = [tuple(bool(x) for x in g_obs[:, idx].tolist()) for idx in range(int(grp.shape[0]))]
                    for u_sig in sorted(set(col_sigs)):
                        sub_g = grp[np.array([cs == u_sig for cs in col_sigs], dtype=bool)]
                        if int(sub_g.shape[0]) > 0:
                            obs_blocks.append(sub_g)
            block_obs = np.column_stack(
                [np.all(presence.observability[:, blk] > 1e-6, axis=1) for blk in obs_blocks]
            )
            if not bool(np.all(block_obs)):
                return _screen_continuous_partitioned_by_observability(
                    presence=presence,
                    overlap=overlap,
                    response_binary=response_binary,
                    strata_codes=strata_codes,
                    informative_groups=obs_blocks,
                    stratum_obs=np.asarray(block_obs, dtype=bool),
                    axes=axes,
                    overlap_config=overlap_config,
                    screen_config=screen_config,
                    prior_pair_weights=prior_pair_weights,
                )

    used_samples = sum(int(g.shape[0]) for g in informative_groups)
    if used_samples < 6:
        return Success(_empty_continuous_df())

    n_r_used = sum(int(np.sum(response_binary[g] > 0.5)) for g in informative_groups)
    n_nr_used = used_samples - n_r_used
    arm_adequate = (used_samples >= 12) and (min(n_r_used, n_nr_used) >= 5)

    overlap_active = subset_overlap_graph(overlap, presence.active_nhood_indices)
    v_act = overlap_active.n_nhoods
    if v_act < 2:
        return Success(_empty_continuous_df())
    w1d = overlap_active.spatial_weights_1d

    prior_res = _normalize_active_prior_matrix(
        prior_pair_weights=prior_pair_weights,
        active_nhood_indices=presence.active_nhood_indices,
        overlap_active=overlap_active,
        overlap_config=overlap_config,
    )
    match prior_res:
        case Failure(err):
            return Failure(err)
        case Success(w_prior_active):
            pass

    has_external_prior = isinstance(w_prior_active, Some)

    # 1. Standardized within-stratum log2_cpm for scale-invariant LOG_RATIO_BALANCE
    z_cpm, y_unit, df_resid = _stratum_center_and_standardize(
        presence.log2_cpm,
        response_binary,
        informative_groups,
    )
    r_cpm = z_cpm @ y_unit

    # Raw stratum-centered log2_cpm for exact compositional balance B_{uv} = (c_u - c_v) / sqrt(2)
    used_idx = np.concatenate(informative_groups)
    cpm_centered = np.zeros((v_act, int(used_idx.shape[0])), dtype=np.float64)
    offset = 0
    for grp in informative_groups:
        g_len = int(grp.shape[0])
        cg = presence.log2_cpm[:, grp]
        cpm_centered[:, offset : offset + g_len] = cg - np.mean(cg, axis=1, keepdims=True)
        offset += g_len
    cpm_sq_norms = np.sum(cpm_centered * cpm_centered, axis=1)
    cpm_dot_y = cpm_centered @ y_unit

    # Outcome-blind & strictly scale-invariant positive-donor support (C_{v,j}/L_j > 0)
    pos_mat_f32 = (presence.arcsin_sqrt_prop[:, used_idx] > 0.0).astype(np.float32)
    pos_donors_1d = np.sum(pos_mat_f32, axis=1)
    min_cont_pos = float(max(6, int(np.ceil(0.08 * float(used_samples)))))
    min_copos = float(max(5, int(np.ceil(0.05 * float(used_samples)))))
    pos_support_ok: BoolArray = np.asarray(pos_donors_1d >= min_cont_pos, dtype=bool)

    # 2. Standardized within-stratum arcsin-sqrt proportions for CO_ABUNDANCE_PRODUCT
    z_arc, _, _ = _stratum_center_and_standardize(
        presence.arcsin_sqrt_prop,
        response_binary,
        informative_groups,
    )
    r_arc = z_arc @ y_unit
    z_arc_sq = z_arc * z_arc
    z_arc_row_sq = np.sum(z_arc_sq, axis=1)
    z_arc_y = z_arc * np.expand_dims(y_unit, axis=0)
    sqrt_used = float(np.sqrt(float(used_samples)))

    u_list: list[IntArray] = []
    v_list: list[IntArray] = []
    axis_list: list[NDArray[np.str_]] = []
    jacc_list: list[FloatArray] = []
    scope_list: list[NDArray[np.str_]] = []
    sw_list: list[FloatArray] = []
    pm_list: list[FloatArray] = []
    eff_list: list[FloatArray] = []
    p_axis_list: list[FloatArray] = []
    t_iut_list: list[FloatArray] = []
    p_iut_list: list[FloatArray] = []
    testable_list: list[BoolArray] = []

    df_two = max(df_resid - 2, 1)
    df_one = max(df_resid - 1, 1)
    df_inter = max(df_resid - 3, 1)

    for start in range(0, v_act, screen_config.chunk_size):
        end = min(start + screen_config.chunk_size, v_act)
        valid_mask, jaccard_block, scope_block = build_block_pair_mask(
            overlap_active=overlap_active,
            u_start=start,
            u_end=end,
            config=overlap_config,
        )
        row_rel, col_abs = np.nonzero(valid_mask)
        if int(row_rel.shape[0]) == 0:
            continue

        u_abs = row_rel + start
        n_p = int(u_abs.shape[0])
        jacc_p = jaccard_block[row_rel, col_abs]
        scope_p = scope_block[row_rel, col_abs]
        sw_p = (
            (w1d[u_abs] * w1d[col_abs])
            if screen_config.use_spatial_fdr_weights
            else np.ones(n_p, dtype=np.float64)
        )
        pm_p = (
            w_prior_active.unwrap()[u_abs, col_abs]
            if isinstance(w_prior_active, Some)
            else np.ones(n_p, dtype=np.float64)
        )
        endpt_ok = pos_support_ok[u_abs] & pos_support_ok[col_abs] & (pm_p > 0.0)

        if ContinuousAxis.LOG_RATIO_BALANCE in axes:
            dot_uv = (cpm_centered[start:end, :] @ cpm_centered.T)[row_rel, col_abs]
            diff_sq = np.maximum(
                cpm_sq_norms[u_abs] + cpm_sq_norms[col_abs] - 2.0 * dot_uv,
                1e-15,
            )
            r_bal = np.clip(
                (cpm_dot_y[u_abs] - cpm_dot_y[col_abs]) / np.sqrt(diff_sq),
                -0.999999,
                0.999999,
            )
            t_bal = r_bal * np.sqrt(float(df_one) / np.maximum(1.0 - r_bal * r_bal, 1e-12))
            p_bal = np.asarray(2.0 * student_t.sf(np.abs(t_bal), df=df_one), dtype=np.float64)

            rho_uv = np.clip(
                (z_cpm[start:end, :] @ z_cpm.T)[row_rel, col_abs],
                -0.9999,
                0.9999,
            )
            one_m_rho2 = np.maximum(1.0 - rho_uv * rho_uv, 1e-8)
            ru = r_cpm[u_abs]
            rv = r_cpm[col_abs]
            beta_u = (ru - rho_uv * rv) / one_m_rho2
            beta_v = (rv - rho_uv * ru) / one_m_rho2
            r2_uv = np.clip(ru * beta_u + rv * beta_v, 0.0, 0.999999)
            se_partial = np.sqrt(np.maximum(1.0 - r2_uv, 1e-12) / (float(df_two) * one_m_rho2))

            t_u_given_v = beta_u / se_partial
            t_neg_v_given_u = (-beta_v) / se_partial
            pos_dir = r_bal >= 0.0
            sign_b = np.where(pos_dir, 1.0, -1.0)
            t_iut_signed = np.minimum(sign_b * t_u_given_v, sign_b * t_neg_v_given_u)
            p_iut_b = np.clip(
                np.asarray(2.0 * student_t.sf(t_iut_signed, df=df_two), dtype=np.float64),
                0.0,
                1.0,
            )
            testable_lr = np.asarray(arm_adequate & endpt_ok, dtype=bool)

            u_list.append(u_abs)
            v_list.append(col_abs)
            axis_list.append(np.full(n_p, ContinuousAxis.LOG_RATIO_BALANCE.value))
            jacc_list.append(jacc_p)
            scope_list.append(scope_p)
            sw_list.append(sw_p)
            pm_list.append(pm_p)
            eff_list.append(r_bal)
            p_axis_list.append(p_bal)
            t_iut_list.append(t_iut_signed)
            p_iut_list.append(p_iut_b)
            testable_list.append(testable_lr)

        if ContinuousAxis.CO_ABUNDANCE_PRODUCT in axes:
            dot_uv_a = (z_arc[start:end, :] @ z_arc.T)[row_rel, col_abs]
            rho_a = np.clip(dot_uv_a, -0.9999, 0.9999)
            denom_a = np.maximum(1.0 - rho_a * rho_a, 1e-8)

            pu_dot = sqrt_used * ((z_arc_sq[start:end, :] @ z_arc.T)[row_rel, col_abs])
            pv_dot = sqrt_used * ((z_arc[start:end, :] @ z_arc_sq.T)[row_rel, col_abs])
            gamma_u = (pu_dot - rho_a * pv_dot) / denom_a
            gamma_v = (pv_dot - rho_a * pu_dot) / denom_a

            raw_sq_norm = float(used_samples) * ((z_arc_sq[start:end, :] @ z_arc_sq.T)[row_rel, col_abs])
            wc_sq_norm = np.maximum(raw_sq_norm - dot_uv_a * dot_uv_a, 0.0)
            norm_u_sq = z_arc_row_sq[u_abs]
            norm_v_sq = z_arc_row_sq[col_abs]
            orth_sq = np.maximum(
                wc_sq_norm
                - 2.0 * gamma_u * pu_dot
                - 2.0 * gamma_v * pv_dot
                + gamma_u * gamma_u * norm_u_sq
                + gamma_v * gamma_v * norm_v_sq
                + 2.0 * gamma_u * gamma_v * dot_uv_a,
                0.0,
            )
            orth_norm = np.sqrt(orth_sq)
            valid_orth = orth_norm > 1e-10

            ru_a = r_arc[u_abs]
            rv_a = r_arc[col_abs]
            bu_a = (ru_a - rho_a * rv_a) / denom_a
            bv_a = (rv_a - rho_a * ru_a) / denom_a
            r2_main = np.clip(ru_a * bu_a + rv_a * bv_a, 0.0, 0.99999)

            wc_dot_y = sqrt_used * ((z_arc_y[start:end, :] @ z_arc.T)[row_rel, col_abs])
            orth_dot_y = wc_dot_y - gamma_u * ru_a - gamma_v * rv_a
            r_orth_y = np.where(valid_orth, orth_dot_y / np.maximum(orth_norm, 1e-10), 0.0)

            partial_r = np.clip(
                r_orth_y / np.sqrt(np.maximum(1.0 - r2_main, 1e-12)),
                -0.999999,
                0.999999,
            )
            t_prod = partial_r * np.sqrt(float(df_inter) / np.maximum(1.0 - partial_r * partial_r, 1e-12))
            p_prod = np.where(
                valid_orth,
                np.asarray(2.0 * student_t.sf(np.abs(t_prod), df=df_inter), dtype=np.float64),
                1.0,
            )
            copos = (pos_mat_f32[start:end, :] @ pos_mat_f32.T)[row_rel, col_abs]
            testable_co = np.asarray(
                arm_adequate & endpt_ok & valid_orth & (copos >= (min_copos - 1e-4)),
                dtype=bool,
            )

            u_list.append(u_abs)
            v_list.append(col_abs)
            axis_list.append(np.full(n_p, ContinuousAxis.CO_ABUNDANCE_PRODUCT.value))
            jacc_list.append(jacc_p)
            scope_list.append(scope_p)
            sw_list.append(sw_p)
            pm_list.append(pm_p)
            eff_list.append(partial_r)
            p_axis_list.append(p_prod)
            t_iut_list.append(np.abs(t_prod))
            p_iut_list.append(p_prod)
            testable_list.append(testable_co)

    if not u_list:
        return Success(_empty_continuous_df())

    u_all = np.concatenate(u_list)
    v_all = np.concatenate(v_list)
    ax_all = np.concatenate(axis_list)
    jacc_all = np.concatenate(jacc_list)
    scope_all = np.concatenate(scope_list)
    sw_all = np.concatenate(sw_list)
    pm_all = np.concatenate(pm_list)
    eff_all = np.concatenate(eff_list)
    p_ax_all = np.concatenate(p_axis_list)
    t_iut_all = np.concatenate(t_iut_list)
    p_iut_all = np.concatenate(p_iut_list)
    testable_all = np.concatenate(testable_list)
    sw_bh_all = sw_all * pm_all

    # Compute weighted BH-FDR and DMPP Knockoff+ IUT FDR restricted to outcome-blind testable pairs
    fdr_ax_all = np.ones_like(p_ax_all)
    fdr_iut_all = np.ones_like(p_iut_all)
    for ax_enum in axes:
        mask_ax = (ax_all == ax_enum.value) & testable_all
        if bool(np.any(mask_ax)):
            fdr_ax_sub = compute_weighted_bh_fdr(p_ax_all[mask_ax], sw_bh_all[mask_ax])
            fdr_ax_all[mask_ax] = fdr_ax_sub
            if ax_enum == ContinuousAxis.LOG_RATIO_BALANCE:
                fdr_iut_all[mask_ax] = _compute_decoupled_balance_iut_fdr(
                    t_iut_testable=t_iut_all[mask_ax],
                    sw_testable=sw_all[mask_ax],
                    prior_mult_testable=pm_all[mask_ax],
                    has_external_prior=has_external_prior,
                    w_prior_active=w_prior_active,
                    z_cpm=z_cpm,
                    cpm_centered=cpm_centered,
                    r_cpm=r_cpm,
                    cpm_dot_y=cpm_dot_y,
                    y_unit=y_unit,
                    df_two=df_two,
                    pos_support_ok=pos_support_ok,
                    informative_groups=informative_groups,
                    overlap_active=overlap_active,
                    overlap_config=overlap_config,
                    screen_config=screen_config,
                )
            else:
                fdr_iut_all[mask_ax] = fdr_ax_sub

    is_sig = (
        testable_all
        & (fdr_ax_all < screen_config.fdr_threshold)
        & (fdr_iut_all < screen_config.fdr_threshold)
    )
    if int(u_all.shape[0]) > screen_config.max_retained_candidates:
        rank_score = np.where(testable_all, p_iut_all, 2.0 + p_iut_all)
        order = np.argsort(rank_score, kind="mergesort")
        top_m = np.zeros(u_all.shape[0], dtype=bool)
        top_m[order[: screen_config.max_retained_candidates]] = True
        keep = is_sig | top_m
    else:
        keep = np.ones(u_all.shape[0], dtype=bool)

    act_idx = presence.active_nhood_indices
    ct_arr = np.asarray(overlap_active.dominant_cell_types, dtype=str)
    dirs = np.where(
        eff_all[keep] >= 0.0,
        EnrichmentDirection.RESPONDER.value,
        EnrichmentDirection.NON_RESPONDER.value,
    )

    p_ax_masked = np.where(testable_all, p_ax_all, 1.0)
    p_iut_masked_c = np.where(testable_all, p_iut_all, 1.0)

    cont_df = (
        pl.DataFrame(
            {
                "nhood_u": act_idx[u_all[keep]],
                "nhood_v": act_idx[v_all[keep]],
                "axis": ax_all[keep],
                "direction": dirs,
                "pair_scope": scope_all[keep],
                "cell_type_u": ct_arr[u_all[keep]],
                "cell_type_v": ct_arr[v_all[keep]],
                "jaccard_overlap": jacc_all[keep],
                "spatial_weight": sw_all[keep],
                "effect_size": eff_all[keep],
                "pvalue_axis": p_ax_masked[keep],
                "fdr_axis": fdr_ax_all[keep],
                "pvalue_iut": p_iut_masked_c[keep],
                "fdr_iut": fdr_iut_all[keep],
                "is_significant": is_sig[keep],
            }
        )
        .lazy()
        .sort(["is_significant", "pvalue_iut", "pvalue_axis"], descending=[True, False, False])
        .collect()
    )
    return Success(cont_df)
