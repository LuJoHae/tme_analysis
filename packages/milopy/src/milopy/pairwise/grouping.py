"""Graph-aware grouping of overlapping 1D neighborhoods and 2D neighborhood-pair hits.

Consolidates locally redundant kNN neighborhood discoveries into cohesive biological
modules using sparse Jaccard topology and Constant Potts Model (CPM) Leiden clustering.
"""

from __future__ import annotations

import igraph as ig
import numpy as np
import polars as pl
from returns.maybe import Nothing
from returns.result import Failure, Result, Success
from scipy import sparse as sp

from .overlap import compute_nhood_overlap_graph
from .types import (
    BoolArray,
    FloatArray,
    IntArray,
    OverlapConfig,
    OverlapGraph,
    PairModuleConfig,
    QuadrantState,
    SparseMatrix,
    swap_quadrant,
)


def group_nhoods(
    nhoods_mat: SparseMatrix,
    nhood_logfc: FloatArray,
    nhood_is_sig: BoolArray,
    min_jaccard: float = 0.10,
    leiden_resolution: float = 1.0,
) -> Result[IntArray, str]:
    """Groups significant concordant 1D Milo neighborhoods sharing kNN cells into modules.

    Non-significant neighborhoods receive group ID `-1`. Significant neighborhoods receive
    non-negative integer module IDs `0, 1, ...` sorted by descending module size.
    """
    overlap_res = compute_nhood_overlap_graph(
        nhoods_mat=nhoods_mat,
        cell_types=Nothing,
        config=OverlapConfig(),
    )
    match overlap_res:
        case Failure(err):
            return Failure(err)
        case Success(overlap):
            pass

    logfc_arr = np.asarray(nhood_logfc, dtype=np.float64).ravel()
    sig_arr = np.asarray(nhood_is_sig, dtype=bool).ravel()
    if int(logfc_arr.shape[0]) != overlap.n_nhoods or int(sig_arr.shape[0]) != overlap.n_nhoods:
        return Failure(
            f"Length mismatch: nhoods={overlap.n_nhoods}, logfc={logfc_arr.shape[0]}, sig={sig_arr.shape[0]}"
        )

    group_ids = np.full(overlap.n_nhoods, -1, dtype=np.int64)
    sig_indices = np.flatnonzero(sig_arr).astype(np.int64)
    n_sig = int(sig_indices.shape[0])
    if n_sig == 0:
        return Success(group_ids)
    if n_sig == 1:
        group_ids[sig_indices[0]] = 0
        return Success(group_ids)

    sub_jacc = overlap.jaccard_matrix[sig_indices, :][:, sig_indices].tocoo()
    signs = np.sign(logfc_arr[sig_indices])

    edges: list[tuple[int, int]] = []
    weights: list[float] = []
    for r_i, c_i, val in zip(sub_jacc.row.tolist(), sub_jacc.col.tolist(), sub_jacc.data.tolist()):
        if r_i < c_i and val >= min_jaccard and signs[r_i] == signs[c_i]:
            edges.append((int(r_i), int(c_i)))
            weights.append(float(val))

    graph = ig.Graph(n=n_sig, edges=edges, directed=False)
    if edges:
        cpm_gamma = min(0.05 * leiden_resolution, 0.5 * min_jaccard)
        partition = graph.community_leiden(
            objective_function="CPM",
            weights=weights,
            resolution=cpm_gamma,
        )
        raw_membership = list(partition.membership)
    else:
        raw_membership = list(range(n_sig))

    # Relabel modules deterministically by descending size then smallest member index
    unique_mods = sorted(
        set(raw_membership),
        key=lambda m: (-raw_membership.count(m), int(sig_indices[raw_membership.index(m)])),
    )
    mod_remap = {old_m: new_m for new_m, old_m in enumerate(unique_mods)}
    for local_i, orig_idx in enumerate(sig_indices.tolist()):
        group_ids[orig_idx] = mod_remap[raw_membership[local_i]]

    return Success(group_ids)


def _empty_modules_df() -> pl.DataFrame:
    return pl.DataFrame(
        schema={
            "pair_module_id": pl.Int64,
            "n_member_pairs": pl.Int64,
            "n_unique_nhoods_u": pl.Int64,
            "n_unique_nhoods_v": pl.Int64,
            "lead_nhood_u": pl.Int64,
            "lead_nhood_v": pl.Int64,
            "quadrant": pl.String,
            "direction": pl.String,
            "pair_scope": pl.String,
            "dominant_cell_type_u": pl.String,
            "dominant_cell_type_v": pl.String,
            "lead_log_or_quadrant": pl.Float64,
            "lead_min_cond_log_or": pl.Float64,
            "best_pvalue_iut": pl.Float64,
            "best_fdr_iut": pl.Float64,
        }
    )


def group_nhood_pairs(
    pairs_df: pl.DataFrame,
    overlap: OverlapGraph,
    config: PairModuleConfig = PairModuleConfig(),
) -> Result[tuple[pl.DataFrame, pl.DataFrame], str]:
    """Clusters significant neighborhood-pair hits into non-redundant ecosystem modules.

    Handles both direct alignment `(u1 ~ u2 and v1 ~ v2, q1 == q2)` and transposed
    index alignment `(u1 ~ v2 and v1 ~ u2, q1 == swap(q2))` so grouping is strictly
    invariant to neighborhood index ordering.
    """
    if pairs_df.is_empty() or (not config.enabled):
        tagged = pairs_df.with_columns(pl.lit(-1, dtype=pl.Int64).alias("pair_module_id"))
        return Success((tagged, _empty_modules_df()))

    sig_mask = pairs_df["is_significant"].to_numpy()
    sig_row_indices = np.flatnonzero(sig_mask).astype(np.int64)
    n_sig = int(sig_row_indices.shape[0])
    if n_sig == 0:
        tagged = pairs_df.with_columns(pl.lit(-1, dtype=pl.Int64).alias("pair_module_id"))
        return Success((tagged, _empty_modules_df()))

    sig_df = pairs_df.filter(pl.col("is_significant"))
    u_arr = sig_df["nhood_u"].to_numpy().astype(np.int64)
    v_arr = sig_df["nhood_v"].to_numpy().astype(np.int64)
    q_list = [QuadrantState(str(q)) for q in sig_df["quadrant"].to_list()]
    q_swapped = [swap_quadrant(q) for q in q_list]
    dir_list = [str(d) for d in sig_df["direction"].to_list()]

    # Extract dense Jaccard submatrix only over unique neighborhoods participating in significant hits
    unique_endpoints = np.unique(np.concatenate([u_arr, v_arr]))
    ep_to_loc = {int(ep): idx for idx, ep in enumerate(unique_endpoints.tolist())}
    jacc_sub = np.asarray(
        overlap.jaccard_matrix[unique_endpoints, :][:, unique_endpoints].toarray(),
        dtype=np.float64,
    )
    np.fill_diagonal(jacc_sub, 1.0)

    u_ep = np.array([ep_to_loc[int(x)] for x in u_arr], dtype=np.int64)
    v_ep = np.array([ep_to_loc[int(x)] for x in v_arr], dtype=np.int64)

    edges: list[tuple[int, int]] = []
    weights: list[float] = []
    thresh = config.min_endpoint_jaccard

    for i in range(n_sig):
        ui, vi, qi, di = int(u_ep[i]), int(v_ep[i]), q_list[i], dir_list[i]
        for j in range(i + 1, n_sig):
            if di != dir_list[j]:
                continue
            uj, vj = int(u_ep[j]), int(v_ep[j])
            direct_match = (
                (qi == q_list[j])
                and (jacc_sub[ui, uj] >= thresh)
                and (jacc_sub[vi, vj] >= thresh)
            )
            swap_match = (
                (qi == q_swapped[j])
                and (jacc_sub[ui, vj] >= thresh)
                and (jacc_sub[vi, uj] >= thresh)
            )
            if direct_match or swap_match:
                w_direct = 0.5 * (jacc_sub[ui, uj] + jacc_sub[vi, vj]) if direct_match else 0.0
                w_swap = 0.5 * (jacc_sub[ui, vj] + jacc_sub[vi, uj]) if swap_match else 0.0
                edges.append((i, j))
                weights.append(float(max(w_direct, w_swap)))

    if n_sig == 1:
        raw_membership = [0]
    elif edges:
        graph = ig.Graph(n=n_sig, edges=edges, directed=False)
        cpm_gamma = min(0.05 * config.leiden_resolution, 0.5 * thresh)
        partition = graph.community_leiden(
            objective_function="CPM",
            weights=weights,
            resolution=cpm_gamma,
        )
        raw_membership = list(partition.membership)
    else:
        raw_membership = list(range(n_sig))

    p_iut_sig = sig_df["pvalue_iut"].to_numpy()
    unique_mods = sorted(
        set(raw_membership),
        key=lambda m: (
            float(np.min(p_iut_sig[np.array(raw_membership) == m])),
            -raw_membership.count(m),
        ),
    )
    mod_remap = {old_m: new_m for new_m, old_m in enumerate(unique_mods)}
    sig_mod_ids = np.array([mod_remap[m] for m in raw_membership], dtype=np.int64)

    all_mod_ids = np.full(pairs_df.height, -1, dtype=np.int64)
    all_mod_ids[sig_row_indices] = sig_mod_ids
    tagged_pairs_df = pairs_df.with_columns(pl.Series("pair_module_id", all_mod_ids))

    # Build summarized modules_df sorted by module ID (best IUT significance first)
    sig_with_mod = sig_df.with_columns(pl.Series("pair_module_id", sig_mod_ids))
    modules_df = (
        sig_with_mod.lazy()
        .sort(["pair_module_id", "pvalue_iut", "pvalue_quadrant"], descending=[False, False, False])
        .group_by("pair_module_id", maintain_order=True)
        .agg(
            pl.len().cast(pl.Int64).alias("n_member_pairs"),
            pl.col("nhood_u").n_unique().cast(pl.Int64).alias("n_unique_nhoods_u"),
            pl.col("nhood_v").n_unique().cast(pl.Int64).alias("n_unique_nhoods_v"),
            pl.col("nhood_u").first().alias("lead_nhood_u"),
            pl.col("nhood_v").first().alias("lead_nhood_v"),
            pl.col("quadrant").first().alias("quadrant"),
            pl.col("direction").first().alias("direction"),
            pl.col("pair_scope").first().alias("pair_scope"),
            pl.col("cell_type_u").first().alias("dominant_cell_type_u"),
            pl.col("cell_type_v").first().alias("dominant_cell_type_v"),
            pl.col("log_or_quadrant").first().alias("lead_log_or_quadrant"),
            pl.col("min_cond_log_or").first().alias("lead_min_cond_log_or"),
            pl.col("pvalue_iut").min().alias("best_pvalue_iut"),
            pl.col("fdr_iut").min().alias("best_fdr_iut"),
        )
        .collect()
    )
    return Success((tagged_pairs_df, modules_df))
