"""Sparse kNN neighborhood overlap graph, pair scope masks, and 2D spatial weighted BH-FDR.

Computes exact cell-sharing Jaccard indices J(u, v) in sparse O(nnz(A^T A)) time and
generalizes Milo spatial FDR weighting to the 2D product space w_{uv} = w_u * w_v.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import assert_never

import numpy as np
from numpy.typing import NDArray
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import sparse as sp

from .types import (
    BoolArray,
    FloatArray,
    IntArray,
    OverlapConfig,
    OverlapGraph,
    PairScope,
    PairScopeFilter,
    SparseMatrix,
    TmeCompartment,
    classify_tme_compartment,
)


def compute_weighted_bh_fdr(p_values: FloatArray, weights: FloatArray) -> FloatArray:
    """Computes weighted Benjamini-Hochberg FDR q-values (Genovese, Roeder & Wasserman 2006).

    Normalizes non-negative prior weights so mean(w_norm) == 1, transforms
    p*_i = min(1, p_i / w_norm_i), and applies monotone step-up BH adjustment.
    Reduces identically to standard unweighted BH when weights are uniform.
    """
    n_hyp = int(p_values.shape[0])
    if n_hyp == 0:
        return np.empty(0, dtype=np.float64)

    p_clean = np.clip(np.asarray(p_values, dtype=np.float64), 0.0, 1.0)
    w_raw = np.maximum(np.asarray(weights, dtype=np.float64), 0.0)
    w_sum = float(np.sum(w_raw))
    w_norm = (w_raw * float(n_hyp) / w_sum) if w_sum > 0.0 else np.ones(n_hyp, dtype=np.float64)

    with np.errstate(divide="ignore", invalid="ignore"):
        p_weighted = np.where(
            w_norm > 0.0,
            np.minimum(1.0, p_clean / np.maximum(w_norm, 1e-300)),
            1.0,
        )

    order = np.argsort(p_weighted, kind="mergesort")
    sorted_p = p_weighted[order]
    ranks = np.arange(1, n_hyp + 1, dtype=np.float64)
    raw_q = (sorted_p * float(n_hyp)) / ranks
    monotone_q = np.clip(np.minimum.accumulate(raw_q[::-1])[::-1], 0.0, 1.0)

    inv_order = np.empty_like(order)
    inv_order[order] = np.arange(n_hyp)
    return monotone_q[inv_order]


def _compute_dominant_cell_types(
    nhoods_bin: sp.csr_matrix,
    cell_types: Sequence[str],
) -> Result[tuple[tuple[str, ...], tuple[float, ...]], str]:
    """Computes dominant cell type and purity per neighborhood via sparse B^T A."""
    n_cells, n_nhoods = nhoods_bin.shape
    if len(cell_types) != n_cells:
        return Failure(
            f"cell_types length ({len(cell_types)}) does not match nhoods_mat rows ({n_cells})"
        )

    labels = tuple(str(ct) for ct in cell_types)
    unique_classes = tuple(sorted(set(labels)))
    k_classes = len(unique_classes)
    if k_classes == 0:
        return Success((tuple("unannotated" for _ in range(n_nhoods)), tuple(1.0 for _ in range(n_nhoods))))

    class_to_idx = {cls_name: idx for idx, cls_name in enumerate(unique_classes)}
    col_indices = np.array([class_to_idx[val] for val in labels], dtype=np.int64)
    row_indices = np.arange(n_cells, dtype=np.int64)
    indicator = sp.csr_matrix(
        (np.ones(n_cells, dtype=np.float64), (row_indices, col_indices)),
        shape=(n_cells, k_classes),
    )

    counts_mat = np.asarray(indicator.T.dot(nhoods_bin).toarray(), dtype=np.float64)
    totals = np.sum(counts_mat, axis=0)
    safe_totals = np.where(totals > 0.0, totals, 1.0)
    proportions = counts_mat / np.expand_dims(safe_totals, axis=0)

    best_idx = np.argmax(proportions, axis=0)
    purities = np.max(proportions, axis=0)
    dom_names = tuple(
        unique_classes[int(best_idx[j])] if totals[j] > 0.0 else "empty"
        for j in range(n_nhoods)
    )
    dom_purities = tuple(float(purities[j]) if totals[j] > 0.0 else 0.0 for j in range(n_nhoods))
    return Success((dom_names, dom_purities))


def compute_nhood_overlap_graph(
    nhoods_mat: SparseMatrix,
    cell_types: Maybe[Sequence[str]],
    config: OverlapConfig = OverlapConfig(),
) -> Result[OverlapGraph, str]:
    """Constructs the sparse neighborhood Jaccard overlap graph and 1D spatial uniqueness weights."""
    _ = config
    if len(nhoods_mat.shape) != 2:
        return Failure(f"nhoods_mat must be 2D, got shape {nhoods_mat.shape}")

    n_cells, n_nhoods = int(nhoods_mat.shape[0]), int(nhoods_mat.shape[1])
    if n_cells == 0 or n_nhoods < 2:
        return Failure(
            f"Pairwise analysis requires n_cells >= 1 and n_nhoods >= 2, got ({n_cells}, {n_nhoods})"
        )

    csr_raw = sp.csr_matrix(nhoods_mat)
    nhoods_bin = sp.csr_matrix(
        (
            np.ones(csr_raw.nnz, dtype=np.float64),
            csr_raw.indices.copy(),
            csr_raw.indptr.copy(),
        ),
        shape=(n_cells, n_nhoods),
    )
    nhoods_bin.sum_duplicates()
    nhoods_bin.data[:] = 1.0

    overlap_coo = nhoods_bin.T.dot(nhoods_bin).tocoo()
    nhood_sizes = np.asarray(nhoods_bin.sum(axis=0), dtype=np.float64).ravel()

    denom = nhood_sizes[overlap_coo.row] + nhood_sizes[overlap_coo.col] - overlap_coo.data
    jaccard_vals = np.where(denom > 0.0, overlap_coo.data / np.maximum(denom, 1.0), 0.0)
    jaccard_csr = sp.csr_matrix(
        (jaccard_vals, (overlap_coo.row, overlap_coo.col)),
        shape=(n_nhoods, n_nhoods),
    )

    jaccard_row_sums = np.asarray(jaccard_csr.sum(axis=1), dtype=np.float64).ravel()
    spatial_weights_1d = 1.0 / np.maximum(jaccard_row_sums, 1.0)

    match cell_types:
        case Some(ct_seq):
            ct_res = _compute_dominant_cell_types(nhoods_bin, ct_seq)
            match ct_res:
                case Failure(err):
                    return Failure(err)
                case Success((dom_types, dom_purities)):
                    return Success(
                        OverlapGraph(
                            n_nhoods=n_nhoods,
                            nhood_sizes=nhood_sizes,
                            jaccard_matrix=jaccard_csr,
                            spatial_weights_1d=spatial_weights_1d,
                            dominant_cell_types=dom_types,
                            cell_type_purities=dom_purities,
                            has_cell_types=True,
                        )
                    )
                case _ as unreachable_res:
                    return Failure(f"Unexpected Result state: {unreachable_res}")
        case Maybe.empty:
            return Success(
                OverlapGraph(
                    n_nhoods=n_nhoods,
                    nhood_sizes=nhood_sizes,
                    jaccard_matrix=jaccard_csr,
                    spatial_weights_1d=spatial_weights_1d,
                    dominant_cell_types=tuple("unannotated" for _ in range(n_nhoods)),
                    cell_type_purities=tuple(1.0 for _ in range(n_nhoods)),
                    has_cell_types=False,
                )
            )
        case _ as unreachable:
            return Failure(f"Unexpected Maybe state: {unreachable}")


def subset_overlap_graph(overlap: OverlapGraph, active_indices: IntArray) -> OverlapGraph:
    """Slices an OverlapGraph to the subset of active neighborhoods while preserving 1D spatial weights."""
    sub_jaccard = overlap.jaccard_matrix[active_indices, :][:, active_indices].tocsr()
    return OverlapGraph(
        n_nhoods=int(active_indices.shape[0]),
        nhood_sizes=overlap.nhood_sizes[active_indices],
        jaccard_matrix=sub_jaccard,
        spatial_weights_1d=overlap.spatial_weights_1d[active_indices],
        dominant_cell_types=tuple(overlap.dominant_cell_types[int(i)] for i in active_indices),
        cell_type_purities=tuple(overlap.cell_type_purities[int(i)] for i in active_indices),
        has_cell_types=overlap.has_cell_types,
    )


def build_block_pair_mask(
    overlap_active: OverlapGraph,
    u_start: int,
    u_end: int,
    config: OverlapConfig,
) -> tuple[BoolArray, FloatArray, NDArray[np.str_]]:
    """Builds the boolean validity mask, Jaccard block, and scope label block for rows [u_start, u_end)."""
    v_total = overlap_active.n_nhoods
    block_rows = u_end - u_start
    row_coords = np.expand_dims(np.arange(u_start, u_end, dtype=np.int64), axis=1)
    col_coords = np.expand_dims(np.arange(v_total, dtype=np.int64), axis=0)
    upper_tri = col_coords > row_coords

    jaccard_block = np.asarray(
        overlap_active.jaccard_matrix[u_start:u_end, :].toarray(),
        dtype=np.float64,
    )
    disjoint_ok = jaccard_block <= (config.max_jaccard_overlap + 1e-12)

    if overlap_active.has_cell_types:
        ct_all = np.asarray(overlap_active.dominant_cell_types, dtype=str)
        ct_u = np.expand_dims(ct_all[u_start:u_end], axis=1)
        ct_v = np.expand_dims(ct_all, axis=0)
        is_same_ct = ct_u == ct_v
        base_scope_labels = np.where(
            is_same_ct,
            PairScope.WITHIN_CELL_TYPE.value,
            PairScope.CROSS_CELL_TYPE.value,
        )
        match config.scope_filter:
            case PairScopeFilter.ALL:
                scope_ok = np.ones((block_rows, v_total), dtype=bool)
                scope_labels = base_scope_labels
            case PairScopeFilter.CROSS_CELL_TYPE_ONLY:
                scope_ok = ~is_same_ct
                scope_labels = base_scope_labels
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
                scope_labels = np.where(
                    scope_ok,
                    PairScope.CROSS_COMPARTMENT.value,
                    base_scope_labels,
                )
            case PairScopeFilter.WITHIN_CELL_TYPE_ONLY:
                scope_ok = is_same_ct
                scope_labels = base_scope_labels
            case _ as unreachable:
                assert_never(unreachable)
    else:
        scope_labels = np.full((block_rows, v_total), PairScope.UNANNOTATED.value)
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

    valid_mask = upper_tri & disjoint_ok & scope_ok
    return valid_mask, jaccard_block, scope_labels
