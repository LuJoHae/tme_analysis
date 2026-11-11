"""Strictly outcome-blind prefiltering and graph super-neighborhood pooling for milopy.pairwise.

All functions in this module are 100% blind to clinical response labels `y` (never inspecting
which specific donors are responders vs. non-responders), guaranteeing zero data snooping and
valid finite-sample / asymptotic FDR control:
1. Lineage/state-preserving Leiden super-neighborhood pooling (`coarse_grain_nhoods_by_graph`):
   Deduplicates shared kNN cells across locally overlapping spheres to reduce Poisson zero-dropout
   (sharpening raw `pvalue_quadrant` and `pvalue_iut`) while quadratically compressing O(V^2).
2. Bourgon-Gentleman-Huber independent marginal entropy & dispersion filtering
   (`select_informative_nhoods_outcome_blind`).
3. Exact discrete Tarone hypergeometric corner upper bounds on Berger's IUT CMH score
   (`compute_tarone_max_iut_z_strata`, `tarone_z_threshold_from_p_ceiling`).
4. External cross-modality prior pair weight construction (`build_external_prior_weight_matrix`)
   for Genovese-Roeder-Wasserman (2006) weighted FDR.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TypeAlias

import igraph as ig
import numpy as np
from numpy.typing import NDArray
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import sparse as sp
from scipy.stats import norm

from .overlap import compute_nhood_overlap_graph
from .types import (
    BoolArray,
    CoarseGrainedNhoodData,
    FloatArray,
    IntArray,
    OverlapConfig,
    PrefilterConfig,
    PresenceRepresentation,
    SparseMatrix,
)

Float32Array: TypeAlias = NDArray[np.float32]


def compute_binary_shannon_entropy(prevalence: FloatArray) -> FloatArray:
    """Computes binary Shannon entropy H_2(p) in [0, 1] bits without NaN at 0 or 1."""
    p_arr = np.clip(np.asarray(prevalence, dtype=np.float64), 0.0, 1.0)
    interior = (p_arr > 1e-15) & (p_arr < (1.0 - 1e-15))
    p_safe = np.where(interior, p_arr, 0.5)
    q_safe = 1.0 - p_safe
    h_raw = -(p_safe * np.log2(p_safe) + q_safe * np.log2(q_safe))
    return np.where(interior, np.clip(h_raw, 0.0, 1.0), 0.0)


def subset_presence_representation(
    presence: PresenceRepresentation,
    keep_active_mask: BoolArray,
) -> PresenceRepresentation:
    """Returns a new immutable PresenceRepresentation restricted to `keep_active_mask` rows."""
    mask = np.asarray(keep_active_mask, dtype=bool).ravel()
    return PresenceRepresentation(
        active_nhood_indices=presence.active_nhood_indices[mask],
        presence_prob=presence.presence_prob[mask, :],
        absence_prob=presence.absence_prob[mask, :],
        observability=presence.observability[mask, :],
        raw_binary_presence=presence.raw_binary_presence[mask, :],
        arcsin_sqrt_prop=presence.arcsin_sqrt_prop[mask, :],
        log2_cpm=presence.log2_cpm[mask, :],
        effective_library_sizes=presence.effective_library_sizes,
    )


def select_informative_nhoods_outcome_blind(
    presence: PresenceRepresentation,
    config: PrefilterConfig,
) -> Result[BoolArray, str]:
    """Selects informative active neighborhoods strictly from unlabelled across-donor marginals.

    Follows Bourgon, Gentleman & Huber (PNAS 2010): marginal donor prevalence `pi_v` and total
    pooled `log2_cpm` sample variance `s_v^2` across all donors combined are ancillary to two-group
    clinical response labels under exchangeable H_0.
    """
    v_act = int(presence.active_nhood_indices.shape[0])
    n_samples = int(presence.presence_prob.shape[1])
    if v_act == 0:
        return Success(np.empty(0, dtype=bool))
    if n_samples < 4:
        return Failure(f"At least 4 samples are required for marginal prefiltering, got {n_samples}")

    # Outcome-blind marginal presence rate & binary Shannon entropy across observable donors
    obs_donor_mask = presence.observability > 1e-6
    obs_donor_counts = np.maximum(np.sum(obs_donor_mask, axis=1, dtype=np.float64), 1.0)
    pi_v = np.sum(np.where(obs_donor_mask, presence.presence_prob, 0.0), axis=1) / obs_donor_counts
    entropy_v = compute_binary_shannon_entropy(pi_v)

    # Outcome-blind pooled sample variance of log2_cpm across all J donors
    total_var_v = np.var(presence.log2_cpm, axis=1, ddof=1)
    var_cutoff = (
        float(np.quantile(total_var_v, config.min_variance_quantile))
        if config.min_variance_quantile > 0.0
        else 0.0
    )

    band_ok = (
        (pi_v >= config.min_marginal_prevalence)
        & (pi_v <= config.max_marginal_prevalence)
        & (total_var_v >= (var_cutoff - 1e-15))
    )

    # Composite outcome-blind informativeness score in [0, 1] for tie-breaking / optional capping
    ranks = np.argsort(np.argsort(total_var_v, kind="mergesort"), kind="mergesort").astype(np.float64)
    var_pct = (ranks + 1.0) / float(max(v_act, 1))
    info_score = entropy_v * (0.5 + 0.5 * var_pct)

    # Guarantee >= 2 retained active neighborhoods on extreme/synthetic matrices using outcome-blind score
    if int(np.sum(band_ok)) < min(2, v_act):
        top_fallback = np.argsort(-info_score, kind="mergesort")[: min(max(4, 2), v_act)]
        fallback_mask = np.zeros(v_act, dtype=bool)
        fallback_mask[top_fallback] = True
        band_ok = fallback_mask

    if config.max_active_nhoods_cap > 0 and int(np.sum(band_ok)) > config.max_active_nhoods_cap:
        surviving_idx = np.flatnonzero(band_ok)
        order_surv = np.argsort(-info_score[surviving_idx], kind="mergesort")
        chosen = surviving_idx[order_surv[: config.max_active_nhoods_cap]]
        capped_mask = np.zeros(v_act, dtype=bool)
        capped_mask[chosen] = True
        return Success(capped_mask)

    return Success(np.asarray(band_ok, dtype=bool))


def coarse_grain_nhoods_by_graph(
    nhoods_mat: SparseMatrix,
    count_matrix: FloatArray,
    cell_sample_codes: Maybe[IntArray] = Nothing,
    sample_library_sizes: Maybe[FloatArray] = Nothing,
    nhood_group_labels: Maybe[Sequence[str]] = Nothing,
    cell_types: Maybe[Sequence[str]] = Nothing,
    min_jaccard: float = 0.10,
    leiden_resolution: float = 1.0,
) -> Result[CoarseGrainedNhoodData, str]:
    """Pools locally overlapping neighborhoods within identical lineage/state labels via Leiden CPM.

    Constructs the exact binary cell-union matrix `nhoods_super[i, g] = 1{sum_{v in C_g} A_{i, v} > 0}`
    so shared cells across overlapping k=30 spheres are never double-counted. When `cell_sample_codes`
    is provided, donor counts are computed exactly via `nhoods_super.T @ sample_indicator`; otherwise
    exact cell-overlap multiplicity weights `1 / k_i` are used.
    """
    counts_arr = np.asarray(count_matrix, dtype=np.float64)
    if counts_arr.ndim != 2:
        return Failure(f"count_matrix must be 2D, got shape {counts_arr.shape}")

    n_cells = int(nhoods_mat.shape[0])
    n_nhoods = int(nhoods_mat.shape[1])
    n_samples = int(counts_arr.shape[1])
    if int(counts_arr.shape[0]) != n_nhoods:
        return Failure(
            f"count_matrix rows ({counts_arr.shape[0]}) must match nhoods_mat cols ({n_nhoods})"
        )

    overlap_res = compute_nhood_overlap_graph(
        nhoods_mat=nhoods_mat,
        cell_types=cell_types,
        config=OverlapConfig(),
    )
    match overlap_res:
        case Failure(err):
            return Failure(err)
        case Success(overlap):
            pass

    match nhood_group_labels:
        case Some(lbls):
            if len(lbls) != n_nhoods:
                return Failure(
                    f"nhood_group_labels length ({len(lbls)}) must match n_nhoods ({n_nhoods})"
                )
            resolved_labels = tuple(str(x) for x in lbls)
        case Maybe.empty:
            resolved_labels = (
                overlap.dominant_cell_types
                if overlap.has_cell_types
                else tuple("unannotated" for _ in range(n_nhoods))
            )
        case _ as unreachable_lbl:
            return Failure(f"Unexpected Maybe state: {unreachable_lbl}")

    jacc_coo = overlap.jaccard_matrix.tocoo()
    edges: list[tuple[int, int]] = []
    weights: list[float] = []
    for r_i, c_i, val in zip(jacc_coo.row.tolist(), jacc_coo.col.tolist(), jacc_coo.data.tolist()):
        if r_i < c_i and val >= min_jaccard and resolved_labels[r_i] == resolved_labels[c_i]:
            edges.append((int(r_i), int(c_i)))
            weights.append(float(val))

    if edges:
        graph = ig.Graph(n=n_nhoods, edges=edges, directed=False)
        cpm_gamma = min(0.05 * leiden_resolution, 0.5 * min_jaccard)
        partition = graph.community_leiden(
            objective_function="CPM",
            weights=weights,
            resolution=cpm_gamma,
        )
        raw_membership = list(partition.membership)
    else:
        raw_membership = list(range(n_nhoods))

    # Deterministic canonical ordering of super-neighborhoods by (label, smallest original nhood index)
    buckets: dict[int, list[int]] = {}
    for v_idx, m_id in enumerate(raw_membership):
        buckets.setdefault(m_id, []).append(v_idx)

    ordered_clusters = sorted(
        buckets.values(),
        key=lambda members: (resolved_labels[members[0]], members[0]),
    )
    n_super = len(ordered_clusters)
    nhood_to_super = np.empty(n_nhoods, dtype=np.int64)
    member_tuples: list[tuple[int, ...]] = []
    super_labels: list[str] = []

    map_rows: list[int] = []
    map_cols: list[int] = []
    for g_idx, members in enumerate(ordered_clusters):
        mem_sorted = tuple(sorted(members))
        member_tuples.append(mem_sorted)
        super_labels.append(resolved_labels[mem_sorted[0]])
        for v_idx in mem_sorted:
            nhood_to_super[v_idx] = g_idx
            map_rows.append(v_idx)
            map_cols.append(g_idx)

    bin_nhoods = (sp.csc_matrix(nhoods_mat) > 0).astype(np.float64)
    group_map_csc = sp.csc_matrix(
        (np.ones(n_nhoods, dtype=np.float64), (np.asarray(map_rows), np.asarray(map_cols))),
        shape=(n_nhoods, n_super),
    )
    # multiplicities[i, g] = number of member neighborhoods in super-nhood g containing cell i
    multiplicities = bin_nhoods.dot(group_map_csc).tocsc()
    nhoods_super = (multiplicities > 0).astype(np.float64).tocsc()

    match cell_sample_codes:
        case Some(codes):
            codes_arr = np.asarray(codes, dtype=np.int64).ravel()
            if int(codes_arr.shape[0]) != n_cells:
                return Failure(
                    f"cell_sample_codes length ({codes_arr.shape[0]}) must match cells ({n_cells})"
                )
            if int(codes_arr.shape[0]) > 0 and (
                int(np.min(codes_arr)) < 0 or int(np.max(codes_arr)) >= n_samples
            ):
                return Failure(
                    f"cell_sample_codes must lie in [0, {n_samples - 1}]"
                )
            sample_ind = sp.csr_matrix(
                (np.ones(n_cells, dtype=np.float64), (np.arange(n_cells, dtype=np.int64), codes_arr)),
                shape=(n_cells, n_samples),
            )
            counts_super = np.asarray(nhoods_super.T.dot(sample_ind).toarray(), dtype=np.float64)
        case Maybe.empty:
            # Multiplicity-corrected aggregation: weight each cell by 1 / k_{i,g} so total cells == union size
            inv_mult = multiplicities.copy().tocoo()
            inv_data = np.where(inv_mult.data > 0.0, 1.0 / inv_mult.data, 0.0)
            inv_mult_csr = sp.csr_matrix(
                (inv_data, (inv_mult.row, inv_mult.col)),
                shape=(n_cells, n_super),
            )
            nhood_sizes = np.maximum(np.asarray(bin_nhoods.sum(axis=0)).ravel(), 1.0)
            dedup_weights = np.zeros(n_nhoods, dtype=np.float64)
            for g_idx, mem_tup in enumerate(member_tuples):
                col_vec = inv_mult_csr[:, g_idx]
                eff_contrib = np.asarray(bin_nhoods[:, list(mem_tup)].T.dot(col_vec).toarray()).ravel()
                for loc_i, v_idx in enumerate(mem_tup):
                    dedup_weights[v_idx] = float(eff_contrib[loc_i] / nhood_sizes[v_idx])

            weighted_map = sp.csr_matrix(
                (dedup_weights[np.asarray(map_rows)], (np.asarray(map_cols), np.asarray(map_rows))),
                shape=(n_super, n_nhoods),
            )
            counts_super = np.asarray(weighted_map.dot(counts_arr), dtype=np.float64)
            if isinstance(sample_library_sizes, Some):
                libs = np.asarray(sample_library_sizes.unwrap(), dtype=np.float64).ravel()
                if int(libs.shape[0]) == n_samples:
                    counts_super = np.minimum(counts_super, np.expand_dims(libs, axis=0))
        case _ as unreachable_codes:
            return Failure(f"Unexpected Maybe state: {unreachable_codes}")

    return Success(
        CoarseGrainedNhoodData(
            nhoods_super=nhoods_super,
            counts_super=counts_super,
            super_nhood_labels=tuple(super_labels),
            member_nhood_indices=tuple(member_tuples),
            nhood_to_super_id=nhood_to_super,
        )
    )


def tarone_z_threshold_from_p_ceiling(p_ceiling: float) -> float:
    """Converts a two-sided p-value ceiling into the corresponding non-negative normal Z floor."""
    p_clamped = float(np.clip(p_ceiling, 1e-15, 1.0))
    return 0.0 if p_clamped >= 1.0 else float(norm.isf(0.5 * p_clamped))


def _corner_cmh_pos_z_sq_f64(
    n_q: FloatArray,
    n_flip: FloatArray,
    pos_cap: FloatArray,
    neg_cap: FloatArray,
) -> FloatArray:
    """Exact single-stratum squared supremum CMH Z-score for enriching n_q in `pos_cap` arm."""
    a = np.minimum(n_q, pos_cap)
    b = np.maximum(n_q - a, 0.0)
    c = np.maximum(n_flip - neg_cap, 0.0)
    d = np.maximum(n_flip - c, 0.0)

    n = a + b + c + d
    r1 = a + b
    r0 = c + d
    c1 = a + c
    c0 = b + d
    active = (n > 1.0) & ((a * d) > (b * c))
    sn = np.maximum(n, 1e-12)
    num = np.where(active, (a * d - b * c) / sn, 0.0)
    den = sn * sn * np.maximum(n - 1.0, 1.0)
    var = np.where(active, (r1 * r0 * c1 * c0) / np.maximum(den, 1e-30), 0.0)
    return np.where((var > 1e-15) & (num > 0.0), (num * num) / np.maximum(var, 1e-30), 0.0)


def compute_tarone_max_iut_z_strata(
    stratum_quad_totals: Sequence[tuple[FloatArray, FloatArray, FloatArray, FloatArray]],
    stratum_eff_totals: Sequence[FloatArray],
    stratum_arm_counts: Sequence[tuple[int, int]],
    qi: int,
    fv_i: int,
    fu_i: int,
) -> FloatArray:
    """Computes exact outcome-blind upper bound Z_max^IUT across strata for quadrant `qi`.

    Depends strictly on unlabelled quadrant totals `r_k + nr_k`, total pair observability
    `eff_tot_k`, and fixed integer stratum arm sizes `(|R_k|, |NR_k|)` — 100% invariant to any
    within-stratum permutation of clinical response labels. By single-stratum quasiconcavity on
    the hypergeometric feasible polytope and Cauchy-Schwarz across strata, `z_iut_signed <= Z_max^IUT`
    holds universally.
    """
    k_strata = len(stratum_quad_totals)
    if k_strata == 0:
        return np.empty(0, dtype=np.float64)

    shape = stratum_eff_totals[0].shape
    z_sq_r_fv = np.zeros(shape, dtype=np.float64)
    z_sq_r_fu = np.zeros(shape, dtype=np.float64)
    z_sq_nr_fv = np.zeros(shape, dtype=np.float64)
    z_sq_nr_fu = np.zeros(shape, dtype=np.float64)

    for k in range(k_strata):
        qt_k = stratum_quad_totals[k]
        et_k = stratum_eff_totals[k]
        n_r_int, n_nr_int = stratum_arm_counts[k]
        r_cap = np.minimum(float(n_r_int), et_k)
        nr_cap = np.minimum(float(n_nr_int), et_k)

        nq_k = qt_k[qi]
        nfv_k = qt_k[fv_i]
        nfu_k = qt_k[fu_i]

        z_sq_r_fv = z_sq_r_fv + _corner_cmh_pos_z_sq_f64(nq_k, nfv_k, r_cap, nr_cap)
        z_sq_r_fu = z_sq_r_fu + _corner_cmh_pos_z_sq_f64(nq_k, nfu_k, r_cap, nr_cap)
        z_sq_nr_fv = z_sq_nr_fv + _corner_cmh_pos_z_sq_f64(nq_k, nfv_k, nr_cap, r_cap)
        z_sq_nr_fu = z_sq_nr_fu + _corner_cmh_pos_z_sq_f64(nq_k, nfu_k, nr_cap, r_cap)

    z_r_iut = np.minimum(np.sqrt(z_sq_r_fv), np.sqrt(z_sq_r_fu))
    z_nr_iut = np.minimum(np.sqrt(z_sq_nr_fv), np.sqrt(z_sq_nr_fu))
    return np.maximum(z_r_iut, z_nr_iut)


def _corner_cmh_pos_z_sq_f32(
    n_q: Float32Array,
    n_flip: Float32Array,
    pos_cap: Float32Array,
    neg_cap: Float32Array,
) -> Float32Array:
    """Float32 single-stratum squared supremum CMH Z-score for fast DMPP candidate pruning."""
    z0 = np.float32(0.0)
    a = np.minimum(n_q, pos_cap)
    b = np.maximum(n_q - a, z0)
    c = np.maximum(n_flip - neg_cap, z0)
    d = np.maximum(n_flip - c, z0)

    n = a + b + c + d
    r1 = a + b
    r0 = c + d
    c1 = a + c
    c0 = b + d
    active = (n > np.float32(1.0)) & ((a * d) > (b * c))
    sn = np.maximum(n, np.float32(1e-6))
    num = np.where(active, (a * d - b * c) / sn, z0)
    den = sn * sn * np.maximum(n - np.float32(1.0), np.float32(1.0))
    var = np.where(active, (r1 * r0 * c1 * c0) / np.maximum(den, np.float32(1e-12)), z0)
    return np.where(
        (var > np.float32(1e-8)) & (num > z0),
        (num * num) / np.maximum(var, np.float32(1e-12)),
        z0,
    )


def compute_tarone_max_iut_z_candidates_f32(
    stratum_quads_f32: Sequence[
        tuple[
            tuple[Float32Array, Float32Array, Float32Array, Float32Array],
            tuple[Float32Array, Float32Array, Float32Array, Float32Array],
            Float32Array,
            Float32Array,
        ]
    ],
    stratum_arm_counts: Sequence[tuple[int, int]],
    ii: IntArray,
    jj: IntArray,
    qi: int,
    fv_i: int,
    fu_i: int,
) -> Float32Array:
    """Float32 outcome-blind Tarone Z_max^IUT evaluated only on surviving DMPP candidate indices."""
    n_cand = int(ii.shape[0])
    if n_cand == 0:
        return np.empty(0, dtype=np.float32)

    z_sq_r_fv = np.zeros(n_cand, dtype=np.float32)
    z_sq_r_fu = np.zeros(n_cand, dtype=np.float32)
    z_sq_nr_fv = np.zeros(n_cand, dtype=np.float32)
    z_sq_nr_fu = np.zeros(n_cand, dtype=np.float32)

    for k, (rq_tup, nq_tup, rt_k, nt_k) in enumerate(stratum_quads_f32):
        et_k = rt_k[ii, jj] + nt_k[ii, jj]
        n_r_int, n_nr_int = stratum_arm_counts[k]
        r_cap = np.minimum(np.float32(float(n_r_int)), et_k)
        nr_cap = np.minimum(np.float32(float(n_nr_int)), et_k)

        nq_k = rq_tup[qi][ii, jj] + nq_tup[qi][ii, jj]
        nfv_k = rq_tup[fv_i][ii, jj] + nq_tup[fv_i][ii, jj]
        nfu_k = rq_tup[fu_i][ii, jj] + nq_tup[fu_i][ii, jj]

        z_sq_r_fv = z_sq_r_fv + _corner_cmh_pos_z_sq_f32(nq_k, nfv_k, r_cap, nr_cap)
        z_sq_r_fu = z_sq_r_fu + _corner_cmh_pos_z_sq_f32(nq_k, nfu_k, r_cap, nr_cap)
        z_sq_nr_fv = z_sq_nr_fv + _corner_cmh_pos_z_sq_f32(nq_k, nfv_k, nr_cap, r_cap)
        z_sq_nr_fu = z_sq_nr_fu + _corner_cmh_pos_z_sq_f32(nq_k, nfu_k, nr_cap, r_cap)

    z_r = np.minimum(np.sqrt(z_sq_r_fv), np.sqrt(z_sq_r_fu))
    z_nr = np.minimum(np.sqrt(z_sq_nr_fv), np.sqrt(z_sq_nr_fu))
    return np.maximum(z_r, z_nr)


def build_external_prior_weight_matrix(
    nhood_state_labels: Sequence[str],
    state_pair_weights: Mapping[tuple[str, str], float],
    base_floor: float = 0.20,
    max_cap: float = 8.0,
) -> Result[FloatArray, str]:
    """Constructs a symmetric non-negative `(V, V)` prior weight matrix from external state-pair priors.

    Designed for cross-modality prior weighting (Genovese, Roeder & Wasserman 2006) using
    statistically independent bulk deconvolution interaction strengths (e.g. Pan-iAtlas |epistasis_z|).
    """
    n_nhoods = len(nhood_state_labels)
    if n_nhoods < 2:
        return Failure(f"At least 2 neighborhoods are required, got {n_nhoods}")
    if base_floor < 0.0 or max_cap <= 0.0 or base_floor > max_cap:
        return Failure(f"Invalid prior bounds: base_floor={base_floor}, max_cap={max_cap}")

    labels = tuple(str(x) for x in nhood_state_labels)
    unique_states = sorted(set(labels))
    state_to_idx = {st: idx for idx, st in enumerate(unique_states)}
    k_states = len(unique_states)

    state_mat = np.full((k_states, k_states), base_floor, dtype=np.float64)
    for (sa, sb), raw_val in state_pair_weights.items():
        if sa in state_to_idx and sb in state_to_idx:
            ia = state_to_idx[sa]
            ib = state_to_idx[sb]
            clipped = float(np.clip(raw_val, base_floor, max_cap))
            best = max(float(state_mat[ia, ib]), float(state_mat[ib, ia]), clipped)
            state_mat[ia, ib] = best
            state_mat[ib, ia] = best

    nhood_codes = np.array([state_to_idx[lbl] for lbl in labels], dtype=np.int64)
    w_prior = state_mat[np.expand_dims(nhood_codes, axis=1), np.expand_dims(nhood_codes, axis=0)]
    np.fill_diagonal(w_prior, 0.0)
    return Success(np.asarray(w_prior, dtype=np.float64))
