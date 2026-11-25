"""Bagaev et al. (2021, Cancer Cell) Pan-Cancer Molecular Functional Portrait (MFP) TME Subtyping.

Implements:
1. Supervised K=35 Nearest-Neighbor classification against labeled TCGA reference portraits
   (indication-matched or pan-cancer) in 29 median/MAD-scaled Fges signature space.
2. Unsupervised 2-phase Blondel Louvain community detection on thresholded sample Pearson
   correlation graphs with Silhouette, Calinski-Harabasz, and Davies-Bouldin threshold selection
   and canonical quadrant matching onto the 4 conserved TME subtypes (IE/F, IE, F, D).
3. Statistical evaluation of subtype objective response rates (ORR, Wilson CIs, Fisher exact,
   chi-square, and cohort-stratified Cochran-Mantel-Haenszel common odds ratios).
4. Cross-modality integration with outcome-blind single-cell deconvolution fractions,
   patient-grouped StratifiedGroupKFold cross-validation, and 29x29 Fges conservation.
"""

from __future__ import annotations

from enum import Enum
from typing import Final, Literal, Sequence, assert_never
import numpy as np
import numpy.typing as npt
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import optimize, stats
from sklearn.cluster import AgglomerativeClustering  # type: ignore[import-untyped]
from sklearn.decomposition import PCA  # type: ignore[import-untyped]
from sklearn.linear_model import LogisticRegression  # type: ignore[import-untyped]
from sklearn.metrics import (  # type: ignore[import-untyped]
    adjusted_rand_score,
    average_precision_score,
    calinski_harabasz_score,
    cohen_kappa_score,
    davies_bouldin_score,
    normalized_mutual_info_score,
    roc_auc_score,
    silhouette_score,
)
from sklearn.model_selection import (  # type: ignore[import-untyped]
    RepeatedStratifiedKFold,
    StratifiedGroupKFold,
)
from sklearn.preprocessing import StandardScaler  # type: ignore[import-untyped]

from tme_datasets import (
    BAGAEV_MFP_FUNCTIONAL_GROUPS,
    BAGAEV_MFP_SIGNATURE_ORDER,
)

KnnMetric = Literal["correlation", "euclidean"]
DeconvTransformMethod = Literal["clr", "arcsin_sqrt", "raw"]
DeconvClusterAlgorithm = Literal["aitchison_ward", "clr_louvain"]


class BagaevSubtype(str, Enum):
    """Four conserved pan-cancer TME subtypes from Bagaev et al. (2021)."""

    IEF = "IE/F"
    IE = "IE"
    F = "F"
    D = "D"


CANONICAL_SUBTYPE_ORDER: Final[tuple[str, ...]] = ("IE", "IE/F", "F", "D")

IATLAS_TO_TCGA_PROJECTS: Final[dict[str, tuple[str, ...]]] = {
    "Hugo-iAtlas": ("SKCM",),
    "Riaz-iAtlas": ("SKCM",),
    "Liu-iAtlas": ("SKCM",),
    "Gide-iAtlas": ("SKCM",),
    "Rosenberg-iAtlas": ("BLCA",),
    "McDermott-iAtlas": ("KIRC",),
    "Choueiri-iAtlas": ("KIRC",),
    "Padron-iAtlas": ("PAAD",),
    "Anders-iAtlas": ("BRCA",),
}

INDICATION_TO_TCGA_PROJECTS: Final[dict[str, tuple[str, ...]]] = {
    "Melanoma": ("SKCM",),
    "Urothelial Bladder": ("BLCA",),
    "Bladder": ("BLCA",),
    "ccRCC": ("KIRC",),
    "RCC": ("KIRC",),
    "PDAC": ("PAAD",),
    "Pancreatic": ("PAAD",),
    "Breast": ("BRCA",),
    "TNBC": ("BRCA",),
    "NSCLC": ("LUAD", "LUSC"),
    "LUAD": ("LUAD",),
    "CRC": ("COAD", "READ"),
    "Colorectal": ("COAD", "READ"),
    "COAD": ("COAD", "READ"),
    "HNSCC": ("HNSC",),
    "HCC": ("LIHC",),
    "THCA": ("THCA",),
    "UVM": ("UVM",),
    "CHOL": ("CHOL",),
    "OV": ("OV",),
}


class BagaevKnnConfig(BaseModel):
    """Immutable configuration for Bagaev supervised K-NN TCGA classification."""

    model_config = ConfigDict(frozen=True)

    k_neighbors: int = 35
    metric: KnnMetric = "correlation"
    distance_weighted: bool = False


class BagaevLouvainConfig(BaseModel):
    """Immutable configuration for Bagaev correlation-graph Louvain clustering."""

    model_config = ConfigDict(frozen=True)

    thresholds: tuple[float, ...] = (0.30, 0.35, 0.40, 0.45, 0.50, 0.51, 0.55, 0.60)
    min_cluster_fraction: float = 0.05
    target_n_clusters: int = 4
    seed: int = 42


class BagaevLouvainResult(BaseModel):
    """Immutable result of Bagaev correlation-graph clustering and threshold sweep."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    assignments_df: pl.DataFrame
    threshold_sweep_df: pl.DataFrame
    selected_threshold: float
    n_raw_communities: int


def _compute_functional_block_means(
    mat: npt.NDArray[np.float64],
    sig_cols: Sequence[str],
) -> dict[str, npt.NDArray[np.float64]]:
    """Compute sample-wise mean scaled score for each of the 4 Bagaev functional blocks."""
    col_idx = {c: i for i, c in enumerate(sig_cols)}
    return {
        block_name: (
            np.mean(mat[:, idxs], axis=1)
            if (idxs := [col_idx[s] for s in block_sigs if s in col_idx])
            else np.zeros(mat.shape[0], dtype=np.float64)
        )
        for block_name, block_sigs in BAGAEV_MFP_FUNCTIONAL_GROUPS.items()
    }


def _pairwise_pearson_correlation(
    mat_a: npt.NDArray[np.float64],
    mat_b: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Compute exact row-wise Pearson correlation matrix between mat_a (N x D) and mat_b (M x D)."""
    a_centered = mat_a - np.mean(mat_a, axis=1, keepdims=True)
    b_centered = mat_b - np.mean(mat_b, axis=1, keepdims=True)
    a_norm = np.linalg.norm(a_centered, axis=1, keepdims=True)
    b_norm = np.linalg.norm(b_centered, axis=1, keepdims=True)
    a_safe = np.where(a_norm > 1e-12, a_centered / np.maximum(a_norm, 1e-12), 0.0)
    b_safe = np.where(b_norm > 1e-12, b_centered / np.maximum(b_norm, 1e-12), 0.0)
    corr = a_safe @ b_safe.T
    return np.asarray(np.clip(corr, -1.0, 1.0), dtype=np.float64)


def classify_bagaev_knn(
    query_scaled_df: pl.DataFrame,
    ref_scaled_df: pl.DataFrame,
    config: BagaevKnnConfig = BagaevKnnConfig(),
    signature_cols: Sequence[str] = BAGAEV_MFP_SIGNATURE_ORDER,
    ref_label_col: str = "MFP",
) -> Result[pl.DataFrame, str]:
    """Classify query samples into Bagaev TME subtypes (IE/F, IE, F, D) via TCGA K-NN."""
    if query_scaled_df.is_empty():
        return Failure("Query DataFrame is empty.")
    if ref_scaled_df.is_empty():
        return Failure("Reference TCGA DataFrame is empty.")
    if "sample_id" not in query_scaled_df.columns:
        return Failure("query_scaled_df must contain 'sample_id'.")
    if ref_label_col not in ref_scaled_df.columns:
        return Failure(f"ref_scaled_df must contain '{ref_label_col}'.")

    shared_sigs = tuple(
        c for c in signature_cols if c in query_scaled_df.columns and c in ref_scaled_df.columns
    )
    if len(shared_sigs) < 8:
        return Failure(f"Insufficient shared Bagaev signatures ({len(shared_sigs)} < 8).")

    valid_ref = ref_scaled_df.filter(pl.col(ref_label_col).is_in(list(CANONICAL_SUBTYPE_ORDER)))
    if valid_ref.height < 4:
        return Failure(f"Reference has fewer than 4 valid subtype samples ({valid_ref.height}).")

    k_eff = min(max(1, config.k_neighbors), valid_ref.height)
    q_mat = query_scaled_df.select(list(shared_sigs)).to_numpy().astype(np.float64)
    r_mat = valid_ref.select(list(shared_sigs)).to_numpy().astype(np.float64)
    r_labels = np.asarray(valid_ref[ref_label_col].to_list(), dtype=object)

    match config.metric:
        case "correlation":
            corr = _pairwise_pearson_correlation(q_mat, r_mat)
            dists = 1.0 - corr
        case "euclidean":
            diff = q_mat[:, np.newaxis, :] - r_mat[np.newaxis, :, :]
            dists = np.linalg.norm(diff, axis=2)
        case _ as unreachable:
            assert_never(unreachable)

    nn_idx = np.argpartition(dists, kth=k_eff - 1, axis=1)[:, :k_eff]
    n_query = q_mat.shape[0]
    n_classes = len(CANONICAL_SUBTYPE_ORDER)
    probs = np.zeros((n_query, n_classes), dtype=np.float64)
    mean_dists = np.full((n_query, n_classes), np.inf, dtype=np.float64)

    for i in range(n_query):
        idxs = nn_idx[i]
        d_i = dists[i, idxs]
        order = np.argsort(d_i)
        idxs = idxs[order]
        d_i = d_i[order]
        weights = (
            1.0 / np.maximum(d_i, 1e-6)
            if config.distance_weighted
            else np.ones(k_eff, dtype=np.float64)
        )
        w_sum = float(np.sum(weights))
        neigh_labels = r_labels[idxs]
        for s_idx, subtype in enumerate(CANONICAL_SUBTYPE_ORDER):
            st_mask = neigh_labels == subtype
            probs[i, s_idx] = float(np.sum(weights[st_mask])) / w_sum
            if np.any(st_mask):
                mean_dists[i, s_idx] = float(np.mean(d_i[st_mask]))

    # Break vote probability ties by closest mean neighbor correlation (smallest mean distance)
    best_idx = np.array(
        [
            min(
                range(n_classes),
                key=lambda s_idx: (-round(float(probs[row_i, s_idx]), 12), float(mean_dists[row_i, s_idx])),
            )
            for row_i in range(n_query)
        ],
        dtype=np.int64,
    )
    sorted_probs = np.sort(probs, axis=1)
    vote_margin = sorted_probs[:, -1] - sorted_probs[:, -2]
    safe_p = np.clip(probs, 1e-12, 1.0)
    entropy = -np.sum(probs * np.log2(safe_p), axis=1)

    block_means = _compute_functional_block_means(q_mat, shared_sigs)
    preds = [CANONICAL_SUBTYPE_ORDER[int(idx)] for idx in best_idx]

    return Success(
        pl.DataFrame(
            {
                "sample_id": query_scaled_df["sample_id"].cast(pl.String).to_list(),
                "bagaev_subtype": preds,
                "prob_IE": probs[:, 0].tolist(),
                "prob_IEF": probs[:, 1].tolist(),
                "prob_F": probs[:, 2].tolist(),
                "prob_D": probs[:, 3].tolist(),
                "vote_margin": vote_margin.tolist(),
                "vote_entropy": entropy.tolist(),
                "antitumor_immune_score": block_means["Antitumor_Immune"].tolist(),
                "protumor_immune_score": block_means["Protumor_Immune"].tolist(),
                "stromal_fibrotic_score": block_means["Angiogenesis_Fibroblasts"].tolist(),
                "malignant_properties_score": block_means["Malignant_Properties"].tolist(),
            }
        )
    )


def _louvain_modularity_partition(
    adj: npt.NDArray[np.float64],
    seed: int = 42,
    max_passes: int = 30,
    max_levels: int = 12,
) -> npt.NDArray[np.int64]:
    """Deterministic 2-phase Blondel Louvain modularity optimization with quotient coarsening."""
    n = adj.shape[0]
    if n == 0:
        return np.array([], dtype=np.int64)
    w_cur = np.maximum(adj, 0.0).copy()
    np.fill_diagonal(w_cur, 0.0)
    two_m = float(np.sum(w_cur))
    if two_m <= 1e-12:
        return np.arange(n, dtype=np.int64)

    rng = np.random.default_rng(seed)
    orig_comm = np.arange(n, dtype=np.int64)

    for _lvl in range(max_levels):
        m_nodes = w_cur.shape[0]
        deg = np.sum(w_cur, axis=1)
        comm = np.arange(m_nodes, dtype=np.int64)
        order = rng.permutation(m_nodes)
        tot_weight: dict[int, float] = {int(i): float(deg[i]) for i in range(m_nodes)}

        for _pass in range(max_passes):
            moved = False
            for u_raw in order:
                u = int(u_raw)
                c_old = int(comm[u])
                k_u = float(deg[u])
                if k_u <= 1e-12:
                    continue

                tot_weight[c_old] = tot_weight.get(c_old, 0.0) - k_u
                neigh_nodes = np.flatnonzero(w_cur[u] > 0.0)
                cand_weights: dict[int, float] = {c_old: 0.0}
                for v_raw in neigh_nodes:
                    v = int(v_raw)
                    if v != u:
                        c_v = int(comm[v])
                        cand_weights[c_v] = cand_weights.get(c_v, 0.0) + float(w_cur[u, v])

                best_c = c_old
                best_gain = 0.0
                k_u_old = cand_weights.get(c_old, 0.0)
                base_loss = k_u_old - (tot_weight.get(c_old, 0.0) * k_u) / two_m

                for c_cand, k_u_in in cand_weights.items():
                    gain = (k_u_in - (tot_weight.get(c_cand, 0.0) * k_u) / two_m) - base_loss
                    if gain > best_gain + 1e-10:
                        best_gain = gain
                        best_c = c_cand

                comm[u] = best_c
                tot_weight[best_c] = tot_weight.get(best_c, 0.0) + k_u
                if best_c != c_old:
                    moved = True
            if not moved:
                break

        _, relabeled = np.unique(comm, return_inverse=True)
        relabeled_arr = relabeled.astype(np.int64)
        k_new = int(np.max(relabeled_arr)) + 1 if m_nodes > 0 else 0
        orig_comm = relabeled_arr[orig_comm]
        if k_new >= m_nodes:
            break

        # Phase 2: coarsen communities into weighted super-nodes preserving self-loops
        proj = np.eye(k_new, dtype=np.float64)[relabeled_arr]
        w_cur = proj.T @ w_cur @ proj

    _, final_labels = np.unique(orig_comm, return_inverse=True)
    return final_labels.astype(np.int64)


def _consolidate_to_four_clusters(
    raw_labels: npt.NDArray[np.int64],
    mat: npt.NDArray[np.float64],
    corr: npt.NDArray[np.float64],
    min_cluster_fraction: float,
    target_k: int = 4,
) -> npt.NDArray[np.int64]:
    """Consolidate Louvain communities into exact `target_k` clusters."""
    n = mat.shape[0]
    k_goal = min(target_k, n)
    min_size = max(2, int(np.floor(n * min_cluster_fraction)))
    u_ids, counts = np.unique(raw_labels, return_counts=True)
    order = np.argsort(-counts)
    u_sorted = u_ids[order]
    c_sorted = counts[order]

    major_ids = [int(cid) for cid, cnt in zip(u_sorted, c_sorted, strict=True) if int(cnt) >= min_size]
    if len(major_ids) > k_goal:
        major_ids = major_ids[:k_goal]

    if len(major_ids) == k_goal:
        centroids = np.vstack([np.mean(mat[raw_labels == cid], axis=0) for cid in major_ids])
        sim_to_major = _pairwise_pearson_correlation(mat, centroids)
        major_lookup = {cid: idx for idx, cid in enumerate(major_ids)}
        assigned = np.array(
            [
                major_lookup[int(raw_labels[i])]
                if int(raw_labels[i]) in major_lookup
                else int(np.argmax(sim_to_major[i]))
                for i in range(n)
            ],
            dtype=np.int64,
        )
        if len(np.unique(assigned)) == k_goal:
            return assigned

    dist_mat = np.clip(1.0 - corr, 0.0, 2.0)
    np.fill_diagonal(dist_mat, 0.0)
    agg = AgglomerativeClustering(
        n_clusters=k_goal,
        metric="precomputed",
        linkage="average",
    )
    return np.asarray(agg.fit_predict(dist_mat), dtype=np.int64)


def assign_clusters_to_bagaev_subtypes(
    cluster_ids: npt.NDArray[np.int64],
    mat: npt.NDArray[np.float64],
    sig_cols: Sequence[str],
) -> tuple[list[str], dict[int, str]]:
    """Map unsupervised clusters bijectively onto {IE/F, IE, F, D} via Hungarian quadrant alignment."""
    blocks = _compute_functional_block_means(mat, sig_cols)
    imm = blocks["Antitumor_Immune"]
    fib = blocks["Angiogenesis_Fibroblasts"]

    unique_clusters = sorted(int(c) for c in np.unique(cluster_ids))
    raw_c_imm = np.array([float(np.mean(imm[cluster_ids == c])) for c in unique_clusters])
    raw_c_fib = np.array([float(np.mean(fib[cluster_ids == c])) for c in unique_clusters])
    std_imm = raw_c_imm - float(np.mean(raw_c_imm))
    std_fib = raw_c_fib - float(np.mean(raw_c_fib))

    targets = np.array(
        [
            [1.0, -1.0],   # IE
            [1.0, 1.0],    # IE/F
            [-1.0, 1.0],   # F
            [-1.0, -1.0],  # D
        ],
        dtype=np.float64,
    )
    centroids_2d = np.column_stack([std_imm, std_fib])

    if len(unique_clusters) <= len(CANONICAL_SUBTYPE_ORDER):
        cost = -(centroids_2d @ targets.T)
        row_ind, col_ind = optimize.linear_sum_assignment(cost)
        mapping = {
            unique_clusters[int(r)]: CANONICAL_SUBTYPE_ORDER[int(c)]
            for r, c in zip(row_ind, col_ind, strict=True)
        }
    else:
        mapping = {
            cid: CANONICAL_SUBTYPE_ORDER[int(np.argmax(centroids_2d[idx] @ targets.T))]
            for idx, cid in enumerate(unique_clusters)
        }

    labeled = [mapping[int(c)] for c in cluster_ids]
    return labeled, mapping


def cluster_bagaev_louvain(
    scaled_df: pl.DataFrame,
    config: BagaevLouvainConfig = BagaevLouvainConfig(),
    signature_cols: Sequence[str] = BAGAEV_MFP_SIGNATURE_ORDER,
) -> Result[BagaevLouvainResult, str]:
    """Perform Bagaev et al. correlation-graph Louvain clustering across threshold sweep."""
    if scaled_df.height < 8:
        return Failure(f"Need at least 8 samples for 4-subtype clustering; got {scaled_df.height}.")
    if "sample_id" not in scaled_df.columns:
        return Failure("scaled_df must contain 'sample_id'.")

    shared_sigs = tuple(c for c in signature_cols if c in scaled_df.columns)
    if len(shared_sigs) < 8:
        return Failure(f"Insufficient signature columns ({len(shared_sigs)} < 8).")

    mat = scaled_df.select(list(shared_sigs)).to_numpy().astype(np.float64)
    corr = _pairwise_pearson_correlation(mat, mat)
    corr_dist = np.clip(1.0 - corr, 0.0, 2.0)
    np.fill_diagonal(corr_dist, 0.0)

    sweep_rows: list[dict[str, float | int | bool]] = []
    candidate_partitions: dict[float, tuple[npt.NDArray[np.int64], int]] = {}

    for tau in config.thresholds:
        adj = np.where(corr >= tau, corr, 0.0)
        np.fill_diagonal(adj, 0.0)
        raw_comm = _louvain_modularity_partition(adj, seed=config.seed)
        n_raw = int(len(np.unique(raw_comm)))
        merged = _consolidate_to_four_clusters(
            raw_labels=raw_comm,
            mat=mat,
            corr=corr,
            min_cluster_fraction=config.min_cluster_fraction,
            target_k=config.target_n_clusters,
        )
        u_m, cnt_m = np.unique(merged, return_counts=True)
        n_eff = int(len(u_m))
        min_frac = float(np.min(cnt_m)) / float(mat.shape[0])
        is_valid = n_eff >= 2 and min_frac >= (config.min_cluster_fraction * 0.5)

        if 2 <= n_eff < mat.shape[0]:
            sil = float(silhouette_score(corr_dist, merged, metric="precomputed"))
            ch = float(calinski_harabasz_score(mat, merged))
            db = float(davies_bouldin_score(mat, merged))
        else:
            sil, ch, db = -1.0, 0.0, 99.0

        sweep_rows.append(
            {
                "threshold": float(tau),
                "n_raw_communities": n_raw,
                "n_clusters": n_eff,
                "min_cluster_fraction": min_frac,
                "silhouette_corr": sil,
                "calinski_harabasz": ch,
                "davies_bouldin": db,
                "is_valid_four_cluster": bool(is_valid and n_eff == config.target_n_clusters),
            }
        )
        candidate_partitions[float(tau)] = (merged, n_raw)

    sweep_df = pl.DataFrame(sweep_rows)
    valid_sub = sweep_df.filter(pl.col("is_valid_four_cluster"))
    eval_df = valid_sub if not valid_sub.is_empty() else sweep_df

    sil_r = stats.rankdata(eval_df["silhouette_corr"].to_numpy())
    ch_r = stats.rankdata(eval_df["calinski_harabasz"].to_numpy())
    db_r = stats.rankdata(-eval_df["davies_bouldin"].to_numpy())
    best_row_idx = int(np.argmax(sil_r + ch_r + db_r))
    best_tau = float(eval_df["threshold"][best_row_idx])

    best_clusters, best_n_raw = candidate_partitions[best_tau]
    subtypes, _ = assign_clusters_to_bagaev_subtypes(best_clusters, mat, shared_sigs)
    blocks = _compute_functional_block_means(mat, shared_sigs)

    assign_df = pl.DataFrame(
        {
            "sample_id": scaled_df["sample_id"].cast(pl.String).to_list(),
            "louvain_cluster_id": [int(c) for c in best_clusters],
            "louvain_subtype": subtypes,
            "antitumor_immune_score": blocks["Antitumor_Immune"].tolist(),
            "protumor_immune_score": blocks["Protumor_Immune"].tolist(),
            "stromal_fibrotic_score": blocks["Angiogenesis_Fibroblasts"].tolist(),
            "malignant_properties_score": blocks["Malignant_Properties"].tolist(),
        }
    )

    return Success(
        BagaevLouvainResult(
            assignments_df=assign_df,
            threshold_sweep_df=sweep_df,
            selected_threshold=best_tau,
            n_raw_communities=best_n_raw,
        )
    )


def _bh_fdr(raw_pvals: Sequence[float]) -> npt.NDArray[np.float64]:
    """Compute Benjamini-Hochberg FDR-adjusted q-values preserving input order."""
    m = len(raw_pvals)
    if m == 0:
        return np.array([], dtype=np.float64)
    p_arr = np.clip(np.asarray(raw_pvals, dtype=np.float64), 0.0, 1.0)
    order = np.argsort(p_arr)
    q_arr = np.empty(m, dtype=np.float64)
    prev_q = 1.0
    for rank_rev in range(m - 1, -1, -1):
        idx = int(order[rank_rev])
        q_val = min(prev_q, float(p_arr[idx]) * float(m) / float(rank_rev + 1))
        q_arr[idx] = min(1.0, q_val)
        prev_q = q_arr[idx]
    return q_arr


def _safe_corr_pair(
    x: npt.NDArray[np.float64],
    y: npt.NDArray[np.float64],
) -> tuple[float, float, float, float]:
    """Compute (spearman_rho, spearman_p, pearson_r, pearson_p) safely for degenerate inputs."""
    mask = np.isfinite(x) & np.isfinite(y)
    xv = x[mask]
    yv = y[mask]
    if xv.size < 4 or float(np.std(xv)) <= 1e-12 or float(np.std(yv)) <= 1e-12:
        return (0.0, 1.0, 0.0, 1.0)
    sp = stats.spearmanr(xv, yv)
    pr = stats.pearsonr(xv, yv)
    sp_r = float(sp.statistic) if np.isfinite(sp.statistic) else 0.0
    sp_p = float(sp.pvalue) if np.isfinite(sp.pvalue) else 1.0
    pr_r = float(pr.statistic) if np.isfinite(pr.statistic) else 0.0
    pr_p = float(pr.pvalue) if np.isfinite(pr.pvalue) else 1.0
    return (sp_r, sp_p, pr_r, pr_p)


def wilson_ci(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """Compute Wilson score 95% confidence interval for binomial proportion k / n."""
    if n <= 0:
        return (0.0, 0.0)
    p = float(k) / float(n)
    denom = 1.0 + (z * z) / float(n)
    center = (p + (z * z) / (2.0 * float(n))) / denom
    half = (z * np.sqrt((p * (1.0 - p) + (z * z) / (4.0 * float(n))) / float(n))) / denom
    return (float(max(0.0, center - half)), float(min(1.0, center + half)))


def compute_cmh_common_odds_ratio(
    df: pl.DataFrame,
    exposure_col: str,
    outcome_col: str = "response_binary",
    stratum_col: str = "cohort",
) -> tuple[float, float]:
    """Compute Cochran-Mantel-Haenszel common odds ratio and two-sided p-value across cohorts."""
    num_or = 0.0
    den_or = 0.0
    cmh_num = 0.0
    cmh_var = 0.0

    for grp_df in df.partition_by(stratum_col, maintain_order=True):
        n_k = float(grp_df.height)
        if n_k < 2:
            continue
        exp = grp_df[exposure_col].to_numpy().astype(bool)
        out = grp_df[outcome_col].to_numpy().astype(np.float64) == 1.0
        a = float(np.sum(exp & out))
        b = float(np.sum(exp & (~out)))
        c = float(np.sum((~exp) & out))
        d = float(np.sum((~exp) & (~out)))

        num_or += (a * d) / n_k
        den_or += (b * c) / n_k

        r1, r0, c1, c0 = a + b, c + d, a + c, b + d
        if n_k > 1.0 and r1 > 0 and r0 > 0 and c1 > 0 and c0 > 0:
            e_a = (r1 * c1) / n_k
            v_a = (r1 * r0 * c1 * c0) / (n_k * n_k * (n_k - 1.0))
            cmh_num += a - e_a
            cmh_var += v_a

    common_or = (num_or / den_or) if den_or > 1e-12 else (float("inf") if num_or > 0 else 1.0)
    if cmh_var <= 1e-12:
        return (float(common_or), 1.0)
    chi2_stat = (cmh_num * cmh_num) / cmh_var
    p_val = float(stats.chi2.sf(chi2_stat, df=1))
    return (float(common_or), p_val)


def compute_cmh_k_by_2_p_value(
    df: pl.DataFrame,
    cluster_col: str,
    outcome_col: str = "response_binary",
    stratum_col: str = "cohort",
) -> tuple[float, float]:
    """Compute generalized Cochran-Mantel-Haenszel Kx2 chi-square statistic and p-value."""
    valid = df.filter(
        pl.col(cluster_col).is_not_null()
        & pl.col(outcome_col).is_not_null()
        & pl.col(outcome_col).cast(pl.Float64).is_in([0.0, 1.0])
        & pl.col(stratum_col).is_not_null()
    )
    if valid.height < 4:
        return (0.0, 1.0)

    levels = sorted({str(v) for v in valid[cluster_col].cast(pl.String).to_list()})
    n_levels = len(levels)
    if n_levels < 2:
        return (0.0, 1.0)

    k_red = n_levels - 1
    d_sum = np.zeros(k_red, dtype=np.float64)
    v_sum = np.zeros((k_red, k_red), dtype=np.float64)

    for grp_df in valid.partition_by(stratum_col, maintain_order=True):
        n_h = grp_df.height
        if n_h < 2:
            continue
        y_h = grp_df[outcome_col].to_numpy().astype(np.float64) == 1.0
        tot_y = float(np.sum(y_h))
        if tot_y <= 0.0 or tot_y >= float(n_h):
            continue
        c_arr = np.asarray(grp_df[cluster_col].cast(pl.String).to_list(), dtype=object)
        n_vec = np.array([float(np.sum(c_arr == levels[k])) for k in range(k_red)], dtype=np.float64)
        y_vec = np.array([float(np.sum((c_arr == levels[k]) & y_h)) for k in range(k_red)], dtype=np.float64)
        e_vec = n_vec * (tot_y / float(n_h))
        coef = (tot_y * (float(n_h) - tot_y)) / (float(n_h) * (float(n_h) - 1.0))
        v_mat = coef * (np.diag(n_vec) - np.outer(n_vec, n_vec) / float(n_h))
        d_sum += y_vec - e_vec
        v_sum += v_mat

    s_vals = np.linalg.svd(v_sum, compute_uv=False)
    tol = max(1e-10, 1e-8 * float(np.max(s_vals))) if s_vals.size > 0 else 1e-10
    eff_df = int(np.sum(s_vals > tol))
    if eff_df < 1:
        return (0.0, 1.0)

    v_pinv = np.linalg.pinv(v_sum, rcond=1e-10)
    q_cmh = float(max(0.0, float(d_sum @ v_pinv @ d_sum)))
    p_cmh = float(stats.chi2.sf(q_cmh, df=eff_df))
    return (q_cmh, p_cmh)


def evaluate_subtype_response_associations(
    df: pl.DataFrame,
    subtype_col: str = "bagaev_subtype",
    response_col: str = "response_binary",
    cohort_col: Maybe[str] = Nothing,
) -> Result[pl.DataFrame, str]:
    """Compute per-subtype responder rates, Wilson CIs, Fisher ORs, and global/CMH p-values."""
    if subtype_col not in df.columns or response_col not in df.columns:
        return Failure(f"Required columns '{subtype_col}' and '{response_col}' missing.")

    valid = df.filter(
        pl.col(response_col).is_not_null()
        & pl.col(response_col).cast(pl.Float64).is_in([0.0, 1.0])
        & pl.col(subtype_col).is_not_null()
    )
    if valid.is_empty():
        return Failure("No evaluable binary response samples found.")

    y_all = valid[response_col].to_numpy().astype(np.float64) == 1.0
    s_all = np.asarray(valid[subtype_col].cast(pl.String).to_list(), dtype=object)

    present_subtypes = [s for s in CANONICAL_SUBTYPE_ORDER if np.any(s_all == s)]
    for extra in sorted(set(s_all.tolist()) - set(present_subtypes)):
        present_subtypes.append(str(extra))

    table_rows = [
        [int(np.sum((s_all == st_name) & y_all)), int(np.sum((s_all == st_name) & (~y_all)))]
        for st_name in present_subtypes
    ]
    valid_rows = [r for r in table_rows if (r[0] + r[1]) > 0]
    if len(valid_rows) >= 2 and sum(r[0] for r in valid_rows) > 0 and sum(r[1] for r in valid_rows) > 0:
        chi2_stat_raw, chi2_p_raw, _, _ = stats.chi2_contingency(np.asarray(valid_rows, dtype=np.float64))
        chi2_stat = float(chi2_stat_raw)
        chi2_p = float(chi2_p_raw)
    else:
        chi2_stat = 0.0
        chi2_p = 1.0

    fav_mask = np.isin(s_all, ["IE", "IE/F"])
    fav_or, fav_fisher_p = stats.fisher_exact(
        np.asarray(
            [
                [int(np.sum(fav_mask & y_all)), int(np.sum(fav_mask & (~y_all)))],
                [int(np.sum((~fav_mask) & y_all)), int(np.sum((~fav_mask) & (~y_all)))],
            ],
            dtype=np.int64,
        )
    )

    match cohort_col:
        case Some(c_col) if c_col in valid.columns:
            with_fav = valid.with_columns(
                pl.col(subtype_col).cast(pl.String).is_in(["IE", "IE/F"]).alias("_is_favorable")
            )
            cmh_or, cmh_p = compute_cmh_common_odds_ratio(
                with_fav, exposure_col="_is_favorable", outcome_col=response_col, stratum_col=c_col
            )
            cmh_omni_chi2, cmh_omni_p = compute_cmh_k_by_2_p_value(
                valid, cluster_col=subtype_col, outcome_col=response_col, stratum_col=c_col
            )
        case _:
            cmh_or, cmh_p = float(fav_or), float(fav_fisher_p)
            cmh_omni_chi2, cmh_omni_p = chi2_stat, chi2_p

    raw_rows: list[dict[str, str | int | float]] = []
    for st_name in present_subtypes:
        in_st = s_all == st_name
        n_st = int(np.sum(in_st))
        n_resp = int(np.sum(in_st & y_all))
        n_non = n_st - n_resp
        rest_resp = int(np.sum((~in_st) & y_all))
        rest_non = int(np.sum((~in_st) & (~y_all)))
        orr = float(n_resp) / float(n_st) if n_st > 0 else 0.0
        ci_lo, ci_hi = wilson_ci(n_resp, n_st)
        if n_st > 0 and (rest_resp + rest_non) > 0:
            st_or, st_p = stats.fisher_exact(np.asarray([[n_resp, n_non], [rest_resp, rest_non]], dtype=np.int64))
        else:
            st_or, st_p = 1.0, 1.0

        raw_rows.append(
            {
                "subtype": str(st_name),
                "n_samples": n_st,
                "n_responders": n_resp,
                "n_non_responders": n_non,
                "response_rate": orr,
                "ci95_low": ci_lo,
                "ci95_high": ci_hi,
                "vs_rest_odds_ratio": float(st_or),
                "vs_rest_fisher_p": float(st_p),
                "global_chi2_p": float(chi2_p),
                "favorable_IE_IEF_vs_FD_or": float(fav_or),
                "favorable_IE_IEF_vs_FD_fisher_p": float(fav_fisher_p),
                "cmh_favorable_or": float(cmh_or),
                "cmh_favorable_p": float(cmh_p),
            }
        )

    pop_orrs = [float(r["response_rate"]) for r in raw_rows if int(r["n_samples"]) >= 5]
    if len(pop_orrs) < 2:
        pop_orrs = [float(r["response_rate"]) for r in raw_rows if int(r["n_samples"]) >= 1]
    orr_spread = float(max(pop_orrs) - min(pop_orrs)) if len(pop_orrs) >= 2 else 0.0

    rows = [
        {
            **r,
            "orr_spread": orr_spread,
            "cmh_omnibus_chi2": float(cmh_omni_chi2),
            "cmh_omnibus_p": float(cmh_omni_p),
        }
        for r in raw_rows
    ]
    return Success(pl.DataFrame(rows))


def compare_partition_concordance(
    labels_a: Sequence[str],
    labels_b: Sequence[str],
) -> dict[str, float]:
    """Compute exact multi-class agreement, Cohen's kappa, ARI, and NMI between two labelings."""
    if len(labels_a) != len(labels_b) or len(labels_a) == 0:
        return {"n": 0.0, "accuracy": 0.0, "cohen_kappa": 0.0, "ari": 0.0, "nmi": 0.0}
    arr_a = np.asarray([str(x) for x in labels_a], dtype=object)
    arr_b = np.asarray([str(x) for x in labels_b], dtype=object)
    acc = float(np.mean(arr_a == arr_b))
    kappa = float(cohen_kappa_score(arr_a, arr_b)) if len(set(arr_a) | set(arr_b)) > 1 else 1.0
    ari = float(adjusted_rand_score(arr_a, arr_b))
    nmi = float(normalized_mutual_info_score(arr_a, arr_b))
    return {
        "n": float(len(arr_a)),
        "accuracy": acc,
        "cohen_kappa": kappa,
        "ari": ari,
        "nmi": nmi,
    }


def analyze_subtype_deconvolution_enrichment(
    joined_df: pl.DataFrame,
    state_cols: Sequence[str],
    subtype_col: str = "bagaev_subtype",
) -> Result[pl.DataFrame, str]:
    """Test enrichment of deconvolved cell states (e.g. M01..M42) across Bagaev TME subtypes."""
    if joined_df.is_empty() or subtype_col not in joined_df.columns:
        return Failure(f"Cannot compute deconvolution enrichment without '{subtype_col}'.")
    valid_states = [c for c in state_cols if c in joined_df.columns]
    if not valid_states:
        return Failure("None of the requested state_cols exist in joined_df.")

    subtypes = np.asarray(joined_df[subtype_col].cast(pl.String).to_list(), dtype=object)
    present = [s for s in CANONICAL_SUBTYPE_ORDER if np.sum(subtypes == s) >= 2]
    if len(present) < 2:
        return Failure("Need at least 2 populated Bagaev subtypes (>=2 samples) for Kruskal-Wallis.")

    rows: list[dict[str, str | float]] = []
    raw_pvals: list[float] = []

    for state in valid_states:
        vals = joined_df[state].to_numpy().astype(np.float64)
        groups = [vals[subtypes == s] for s in present]
        try:
            kw_stat, kw_p = stats.kruskal(*groups)
        except ValueError:
            kw_stat, kw_p = 0.0, 1.0

        n_tot = float(sum(len(g) for g in groups))
        k_grp = float(len(groups))
        eta_sq = max(0.0, (float(kw_stat) - k_grp + 1.0) / max(1.0, n_tot - k_grp))

        means_by_st = {
            s: float(np.mean(vals[subtypes == s])) if np.any(subtypes == s) else 0.0
            for s in CANONICAL_SUBTYPE_ORDER
        }
        top_st = max(present, key=lambda s: means_by_st[s])
        raw_pvals.append(float(kw_p))
        rows.append(
            {
                "cell_state": state,
                "dominant_subtype": top_st,
                "mean_IE": means_by_st["IE"],
                "mean_IEF": means_by_st["IE/F"],
                "mean_F": means_by_st["F"],
                "mean_D": means_by_st["D"],
                "kruskal_h": float(kw_stat),
                "kruskal_p": float(kw_p),
                "eta_squared": float(eta_sq),
            }
        )

    q_arr = _bh_fdr(raw_pvals)
    out_df = pl.DataFrame(rows).with_columns(pl.Series(name="fdr_q", values=q_arr)).sort("kruskal_p")
    return Success(out_df)


class DeconvClusterConfig(BaseModel):
    """Immutable configuration for unsupervised clustering of deconvoluted TME fraction matrices."""

    model_config = ConfigDict(frozen=True)

    k_values: tuple[int, ...] = (2, 3, 4, 5, 6, 7, 8)
    target_k: int = 4
    min_valid_k: int = 2
    transform: DeconvTransformMethod = "clr"
    algorithm: DeconvClusterAlgorithm = "aitchison_ward"
    center_by_cohort: bool = False
    standardize_features: bool = True
    min_cluster_fraction: float = 0.04
    louvain_thresholds: tuple[float, ...] = (0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
    seed: int = 42


class DeconvClusterResult(BaseModel):
    """Immutable result of deconvolution sample clustering across K sweep and K=4 centroids."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    assignments_df: pl.DataFrame
    k_sweep_df: pl.DataFrame
    optimal_k: int
    optimal_k_multi: int
    centroids_k4_df: pl.DataFrame


def transform_compositional_fractions(
    df: pl.DataFrame,
    state_cols: Sequence[str],
    method: DeconvTransformMethod = "clr",
    pseudocount: float = 1e-6,
    center_by_group: Maybe[str] = Nothing,
    standardize_cols: bool = False,
) -> Result[pl.DataFrame, str]:
    """Apply compositional simplex transformation (CLR / arcsin_sqrt / raw) and optional centering."""
    if df.is_empty():
        return Failure("Input DataFrame is empty.")
    valid_cols = [c for c in state_cols if c in df.columns]
    if len(valid_cols) < 2:
        return Failure(f"Need at least 2 state columns; found {len(valid_cols)}.")

    raw_mat = np.clip(df.select(valid_cols).to_numpy().astype(np.float64), 0.0, None)
    match method:
        case "clr":
            adj = raw_mat + float(pseudocount)
            row_sums = np.maximum(np.sum(adj, axis=1, keepdims=True), 1e-12)
            norm_mat = adj / row_sums
            log_mat = np.log(np.maximum(norm_mat, 1e-12))
            trans = log_mat - np.mean(log_mat, axis=1, keepdims=True)
        case "arcsin_sqrt":
            row_sums = np.sum(raw_mat, axis=1, keepdims=True)
            norm_mat = np.where(row_sums > 1e-12, raw_mat / np.maximum(row_sums, 1e-12), raw_mat)
            trans = np.arcsin(np.sqrt(np.clip(norm_mat, 0.0, 1.0)))
        case "raw":
            row_sums = np.sum(raw_mat, axis=1, keepdims=True)
            trans = np.where(row_sums > 1e-12, raw_mat / np.maximum(row_sums, 1e-12), raw_mat)
        case _ as unreachable:
            assert_never(unreachable)

    match center_by_group:
        case Some(g_col) if g_col in df.columns:
            groups = np.asarray(df[g_col].cast(pl.String).to_list(), dtype=object)
            centered = trans.copy()
            for g_val in np.unique(groups):
                mask = groups == g_val
                if int(np.sum(mask)) >= 1:
                    centered[mask] = centered[mask] - np.mean(centered[mask], axis=0, keepdims=True)
            trans = centered
        case _:
            pass

    if standardize_cols:
        trans = (trans - np.mean(trans, axis=0, keepdims=True)) / np.maximum(
            np.std(trans, axis=0, keepdims=True), 1e-8
        )

    out_df = df.with_columns(
        [pl.Series(name=col, values=trans[:, idx].tolist(), dtype=pl.Float64) for idx, col in enumerate(valid_cols)]
    )
    return Success(out_df)


def compute_deconv_subtype_correlations(
    joined_df: pl.DataFrame,
    state_cols: Sequence[str],
    subtype_col: str = "bagaev_subtype",
    louvain_col: Maybe[str] = Some("louvain_subtype"),
    prob_col_map: Maybe[dict[str, str]] = Nothing,
    cohort_col: Maybe[str] = Some("cohort"),
) -> Result[pl.DataFrame, str]:
    """Correlate deconvoluted cell-state fractions (raw & CLR) with Bagaev subtypes & meta-cohorts."""
    if joined_df.is_empty() or subtype_col not in joined_df.columns:
        return Failure(f"Missing required subtype column '{subtype_col}' or empty DataFrame.")
    valid_states = [c for c in state_cols if c in joined_df.columns]
    if not valid_states:
        return Failure("No valid deconvoluted state columns found.")

    clean_df = joined_df.filter(pl.col(subtype_col).is_not_null())
    clr_res = transform_compositional_fractions(clean_df, valid_states, method="clr", pseudocount=1e-6)
    if isinstance(clr_res, Failure):
        return Failure(clr_res.failure())
    clr_df = clr_res.unwrap()

    raw_mat = clean_df.select(valid_states).to_numpy().astype(np.float64)
    clr_mat = clr_df.select(valid_states).to_numpy().astype(np.float64)
    subtypes = np.asarray(clean_df[subtype_col].cast(pl.String).to_list(), dtype=object)

    louv_maybe: Maybe[npt.NDArray[np.object_]] = (
        Some(np.asarray(clean_df[louvain_col.unwrap()].cast(pl.String).to_list(), dtype=object))
        if isinstance(louvain_col, Some) and louvain_col.unwrap() in clean_df.columns
        else Nothing
    )

    default_prob_map = {"IE": "prob_IE", "IE/F": "prob_IEF", "F": "prob_F", "D": "prob_D"}
    active_prob_map = prob_col_map.value_or(default_prob_map)

    cohort_maybe: Maybe[tuple[npt.NDArray[np.object_], list[str]]] = Nothing
    match cohort_col:
        case Some(c_col) if c_col in clean_df.columns:
            c_arr = np.asarray(clean_df[c_col].cast(pl.String).to_list(), dtype=object)
            cohort_maybe = Some((c_arr, sorted(str(c) for c in np.unique(c_arr))))
        case _:
            pass

    rows: list[dict[str, str | int | float]] = []
    p_raw_list: list[float] = []
    p_clr_list: list[float] = []
    p_meta_list: list[float] = []

    for s_idx, st_name in enumerate(valid_states):
        x_raw = raw_mat[:, s_idx]
        x_clr = clr_mat[:, s_idx]
        for sub_name in CANONICAL_SUBTYPE_ORDER:
            ind = (subtypes == sub_name).astype(np.float64)
            sp_r, sp_p, pr_r, pr_p = _safe_corr_pair(x_raw, ind)
            sp_c, sp_c_p, pr_c, _ = _safe_corr_pair(x_clr, ind)

            p_col_val = active_prob_map.get(sub_name, "")
            if p_col_val and p_col_val in clean_df.columns:
                p_vec = clean_df[p_col_val].to_numpy().astype(np.float64)
                sp_prob, _, _, _ = _safe_corr_pair(x_raw, p_vec)
            else:
                sp_prob = float("nan")

            match louv_maybe:
                case Some(l_arr):
                    sp_louv, _, _, _ = _safe_corr_pair(x_raw, (l_arr == sub_name).astype(np.float64))
                case _:
                    sp_louv = float("nan")

            # Effective binary information Fisher-z meta-analysis across cohorts
            w_sum = 0.0
            wz_sum = 0.0
            n_coh_used = 0
            match cohort_maybe:
                case Some((cohort_arr, unique_cohorts)):
                    for coh_id in unique_cohorts:
                        c_mask = cohort_arr == coh_id
                        n_c = int(np.sum(c_mask))
                        ind_c = ind[c_mask]
                        n_pos = int(np.sum(ind_c == 1.0))
                        n_neg = n_c - n_pos
                        if n_c >= 10 and n_pos >= 3 and n_neg >= 3:
                            rho_c, _, _, _ = _safe_corr_pair(x_raw[c_mask], ind_c)
                            rho_clipped = float(np.clip(rho_c, -0.995, 0.995))
                            p_c = float(n_pos) / float(n_c)
                            w_c = float(n_c - 3) * (4.0 * p_c * (1.0 - p_c))
                            if w_c > 0.0:
                                w_sum += w_c
                                wz_sum += w_c * float(np.arctanh(rho_clipped))
                                n_coh_used += 1
                case _:
                    pass

            if w_sum > 0.0 and n_coh_used >= 1:
                z_meta = wz_sum / w_sum
                se_meta = 1.0 / float(np.sqrt(w_sum))
                meta_rho = float(np.tanh(z_meta))
                meta_p = float(2.0 * stats.norm.sf(abs(z_meta) / se_meta))
            else:
                meta_rho, meta_p = sp_r, sp_p

            p_raw_list.append(sp_p)
            p_clr_list.append(sp_c_p)
            p_meta_list.append(meta_p)
            rows.append(
                {
                    "cell_state": st_name,
                    "subtype": sub_name,
                    "spearman_rho_raw": sp_r,
                    "spearman_p_raw": sp_p,
                    "pearson_r_raw": pr_r,
                    "pearson_p_raw": pr_p,
                    "spearman_rho_clr": sp_c,
                    "spearman_p_clr": sp_c_p,
                    "pearson_r_clr": pr_c,
                    "spearman_rho_prob": sp_prob,
                    "spearman_rho_louvain": sp_louv,
                    "within_cohort_meta_rho": meta_rho,
                    "within_cohort_meta_p": meta_p,
                    "n_cohorts_meta": n_coh_used,
                }
            )

    raw_vec = np.array([float(r["spearman_rho_raw"]) for r in rows], dtype=np.float64)
    meta_vec = np.array([float(r["within_cohort_meta_rho"]) for r in rows], dtype=np.float64)
    pooled_vs_meta_r = (
        float(stats.pearsonr(raw_vec, meta_vec).statistic)
        if np.std(raw_vec) > 1e-12 and np.std(meta_vec) > 1e-12
        else 0.0
    )

    out_df = pl.DataFrame(rows).with_columns(
        [
            pl.Series("fdr_q_raw", _bh_fdr(p_raw_list)),
            pl.Series("fdr_q_clr", _bh_fdr(p_clr_list)),
            pl.Series("fdr_q_meta", _bh_fdr(p_meta_list)),
            pl.lit(pooled_vs_meta_r).alias("pooled_vs_meta_pearson_r"),
        ]
    )
    return Success(out_df)


def compute_state_vs_fges_correlations(
    joined_df: pl.DataFrame,
    state_cols: Sequence[str],
    signature_cols: Sequence[str] = BAGAEV_MFP_SIGNATURE_ORDER,
    use_clr: bool = False,
) -> Result[pl.DataFrame, str]:
    """Compute pairwise correlations between deconvoluted cell states and the 29 Bagaev Fges."""
    if joined_df.is_empty():
        return Failure("Empty DataFrame.")
    valid_states = [c for c in state_cols if c in joined_df.columns]
    valid_sigs = [c for c in signature_cols if c in joined_df.columns]
    if not valid_states or not valid_sigs:
        return Failure("Missing state_cols or signature_cols in joined_df.")

    src_df = joined_df
    if use_clr:
        clr_res = transform_compositional_fractions(joined_df, valid_states, method="clr", pseudocount=1e-6)
        if isinstance(clr_res, Failure):
            return Failure(clr_res.failure())
        src_df = clr_res.unwrap()

    s_mat = src_df.select(valid_states).to_numpy().astype(np.float64)
    g_mat = src_df.select(valid_sigs).to_numpy().astype(np.float64)

    rows: list[dict[str, str | float]] = []
    pvals: list[float] = []
    for i, st_name in enumerate(valid_states):
        x = s_mat[:, i]
        for j, sig_name in enumerate(valid_sigs):
            y = g_mat[:, j]
            sp_r, sp_p, pr_r, _ = _safe_corr_pair(x, y)
            pvals.append(sp_p)
            rows.append(
                {
                    "cell_state": st_name,
                    "fges_signature": sig_name,
                    "spearman_rho": sp_r,
                    "pearson_r": pr_r,
                    "p_value": sp_p,
                }
            )

    return Success(pl.DataFrame(rows).with_columns(pl.Series("fdr_q", _bh_fdr(pvals))))


def _order_clusters_deterministically(
    raw_labels: npt.NDArray[np.int64],
    pc1_scores: npt.NDArray[np.float64],
) -> npt.NDArray[np.int64]:
    """Relabel cluster IDs deterministically from high PC1 (DC1) to low PC1 (DCK)."""
    unique_ids = sorted(int(c) for c in np.unique(raw_labels))
    means = sorted(((float(np.mean(pc1_scores[raw_labels == cid])), cid) for cid in unique_ids), reverse=True)
    remap = {old_id: new_idx + 1 for new_idx, (_, old_id) in enumerate(means)}
    return np.array([remap[int(c)] for c in raw_labels], dtype=np.int64)


def _select_best_k(sweep_df: pl.DataFrame, valid_col: str, min_frac: float) -> int:
    """Select composite-validity optimal K cleanly without fallback bias."""
    valid_df = sweep_df.filter(pl.col(valid_col))
    eval_df = (
        valid_df
        if not valid_df.is_empty()
        else sweep_df.filter(pl.col("min_cluster_fraction") >= min_frac)
    )
    if eval_df.is_empty():
        eval_df = sweep_df
    opt_idx = int(np.argmax(eval_df["composite_validity_rank"].to_numpy()))
    return int(eval_df["k"][opt_idx])


def cluster_deconvoluted_samples(
    df: pl.DataFrame,
    state_cols: Sequence[str],
    config: DeconvClusterConfig = DeconvClusterConfig(),
    cohort_col: str = "cohort",
) -> Result[DeconvClusterResult, str]:
    """Unsupervised compositional clustering of deconvoluted samples (strictly outcome-blind)."""
    if df.height < 8:
        return Failure(f"Need at least 8 samples for deconvolution clustering; got {df.height}.")
    if "sample_id" not in df.columns:
        return Failure("df must contain 'sample_id'.")
    valid_states = [c for c in state_cols if c in df.columns]
    if len(valid_states) < 3:
        return Failure(f"Need at least 3 deconvoluted state columns; got {len(valid_states)}.")

    center_arg: Maybe[str] = Some(cohort_col) if (config.center_by_cohort and cohort_col in df.columns) else Nothing
    trans_res = transform_compositional_fractions(
        df,
        valid_states,
        method=config.transform,
        pseudocount=1e-6,
        center_by_group=center_arg,
        standardize_cols=config.standardize_features,
    )
    if isinstance(trans_res, Failure):
        return Failure(trans_res.failure())

    mat = trans_res.unwrap().select(valid_states).to_numpy().astype(np.float64)
    raw_mat = df.select(valid_states).to_numpy().astype(np.float64)

    n_samples = mat.shape[0]
    pca = PCA(n_components=min(2, mat.shape[1], n_samples), random_state=config.seed)
    pcs = pca.fit_transform(mat)
    pc1 = pcs[:, 0]
    pc2 = pcs[:, 1] if pcs.shape[1] > 1 else np.zeros(n_samples, dtype=np.float64)

    corr = _pairwise_pearson_correlation(mat, mat)
    corr_dist = np.clip(1.0 - corr, 0.0, 2.0)
    np.fill_diagonal(corr_dist, 0.0)

    k_list = sorted({k for k in config.k_values if 2 <= k < n_samples})
    if config.target_k not in k_list and 2 <= config.target_k < n_samples:
        k_list = sorted([*k_list, config.target_k])
    if not k_list:
        return Failure("No valid K values in [2, n_samples - 1].")

    partitions_by_k: dict[int, npt.NDArray[np.int64]] = {}
    sweep_rows: list[dict[str, float | int | bool]] = []

    for k_eff in k_list:
        match config.algorithm:
            case "aitchison_ward":
                raw_lbl = AgglomerativeClustering(
                    n_clusters=k_eff, metric="euclidean", linkage="ward"
                ).fit_predict(mat)
            case "clr_louvain":
                best_score = -1e9
                raw_lbl = AgglomerativeClustering(
                    n_clusters=k_eff, metric="euclidean", linkage="ward"
                ).fit_predict(mat)
                for tau in config.louvain_thresholds:
                    adj = np.where(corr >= tau, corr, 0.0)
                    np.fill_diagonal(adj, 0.0)
                    comm = _louvain_modularity_partition(adj, seed=config.seed)
                    merged = _consolidate_to_four_clusters(
                        raw_labels=comm,
                        mat=mat,
                        corr=corr,
                        min_cluster_fraction=config.min_cluster_fraction * 0.5,
                        target_k=k_eff,
                    )
                    if len(np.unique(merged)) == k_eff:
                        s_c = float(silhouette_score(corr_dist, merged, metric="precomputed"))
                        if s_c > best_score:
                            best_score = s_c
                            raw_lbl = merged
            case _ as unreachable:
                assert_never(unreachable)

        ordered = _order_clusters_deterministically(raw_lbl.astype(np.int64), pc1)
        partitions_by_k[k_eff] = ordered
        _, cnts = np.unique(ordered, return_counts=True)
        min_frac = float(np.min(cnts)) / float(n_samples)
        sil_euc = float(silhouette_score(mat, ordered, metric="euclidean"))
        sil_cor = float(silhouette_score(corr_dist, ordered, metric="precomputed"))
        ch = float(calinski_harabasz_score(mat, ordered))
        db = float(davies_bouldin_score(mat, ordered))

        sweep_rows.append(
            {
                "k": int(k_eff),
                "min_cluster_fraction": min_frac,
                "silhouette_euclidean": sil_euc,
                "silhouette_corr": sil_cor,
                "calinski_harabasz": ch,
                "davies_bouldin": db,
                "is_valid_k": bool(k_eff >= config.min_valid_k and min_frac >= config.min_cluster_fraction),
                "is_valid_multi_k": bool(k_eff >= 3 and min_frac >= config.min_cluster_fraction),
            }
        )

    sweep_base = pl.DataFrame(sweep_rows)
    sil_r = stats.rankdata(sweep_base["silhouette_euclidean"].to_numpy())
    ch_r = stats.rankdata(sweep_base["calinski_harabasz"].to_numpy())
    db_r = stats.rankdata(-sweep_base["davies_bouldin"].to_numpy())
    sweep_df = sweep_base.with_columns(
        pl.Series("composite_validity_rank", ((sil_r + ch_r + db_r) / 3.0).tolist())
    )

    optimal_k = _select_best_k(sweep_df, "is_valid_k", config.min_cluster_fraction)
    optimal_k_multi = _select_best_k(sweep_df, "is_valid_multi_k", config.min_cluster_fraction)

    k4_key = config.target_k if config.target_k in partitions_by_k else k_list[0]
    assign_cols: dict[str, object] = {
        "sample_id": df["sample_id"].cast(pl.String).to_list(),
        "cluster_k4": [f"DC{int(c)}" for c in partitions_by_k[k4_key]],
        "cluster_opt": [f"DC{int(c)}" for c in partitions_by_k[optimal_k]],
        "cluster_opt_multi": [f"DC{int(c)}" for c in partitions_by_k[optimal_k_multi]],
        "pca_1": pc1.tolist(),
        "pca_2": pc2.tolist(),
    }
    for k_eff, part in partitions_by_k.items():
        assign_cols[f"cluster_k{k_eff}"] = [f"DC{int(c)}" for c in part]
    assign_df = pl.DataFrame(assign_cols)

    k4_part = partitions_by_k[k4_key]
    centroid_rows = [
        {
            "cluster_id": f"DC{cid}",
            "cell_state": st_name,
            "mean_clr_z": float(np.mean(mat[k4_part == cid, s_idx])),
            "mean_raw_fraction": float(np.mean(raw_mat[k4_part == cid, s_idx])),
        }
        for cid in sorted(int(c) for c in np.unique(k4_part))
        for s_idx, st_name in enumerate(valid_states)
    ]

    return Success(
        DeconvClusterResult(
            assignments_df=assign_df,
            k_sweep_df=sweep_df,
            optimal_k=optimal_k,
            optimal_k_multi=optimal_k_multi,
            centroids_k4_df=pl.DataFrame(centroid_rows),
        )
    )


def evaluate_within_subtype_deconv_predictions(
    joined_df: pl.DataFrame,
    score_col: str = "pred_pct",
    subtype_col: str = "bagaev_subtype",
    response_col: str = "response_binary",
) -> Result[pl.DataFrame, str]:
    """Evaluate how well deconvolution response predictions stratify R vs NR within each Bagaev subtype."""
    for req in (score_col, subtype_col, response_col):
        if req not in joined_df.columns:
            return Failure(f"Missing column '{req}' in joined_df.")

    clean = joined_df.filter(
        pl.col(response_col).cast(pl.Float64).is_in([0.0, 1.0])
        & pl.col(score_col).is_not_null()
        & pl.col(score_col).cast(pl.Float64).is_finite()
        & pl.col(subtype_col).is_not_null()
    )
    if clean.is_empty():
        return Failure("No valid samples for within-subtype deconvolution evaluation.")

    rows: list[dict[str, str | int | float]] = []
    strata: list[tuple[str, pl.DataFrame]] = [("All_Combined", clean)] + [
        (st_name, clean.filter(pl.col(subtype_col).cast(pl.String) == st_name))
        for st_name in CANONICAL_SUBTYPE_ORDER
    ]

    for label, sub in strata:
        if sub.height < 4:
            continue
        y = sub[response_col].to_numpy().astype(np.int64)
        s = sub[score_col].to_numpy().astype(np.float64)
        n_pos = int(np.sum(y == 1))
        n_neg = int(np.sum(y == 0))
        mean_score = float(np.mean(s))
        if n_pos >= 1 and n_neg >= 1:
            auc = float(roc_auc_score(y, s))
            prauc = float(average_precision_score(y, s))
            _, mwu_p = stats.mannwhitneyu(s[y == 1], s[y == 0], alternative="two-sided")
            mean_r = float(np.mean(s[y == 1]))
            mean_nr = float(np.mean(s[y == 0]))
        else:
            auc, prauc, mwu_p, mean_r, mean_nr = float("nan"), float("nan"), 1.0, mean_score, mean_score

        rows.append(
            {
                "stratum": label,
                "n_samples": int(sub.height),
                "n_responders": n_pos,
                "observed_orr": float(n_pos) / float(sub.height),
                "mean_deconv_score": mean_score,
                "mean_score_responders": mean_r,
                "mean_score_non_responders": mean_nr,
                "roc_auc": auc,
                "pr_auc": prauc,
                "mannwhitney_p": float(mwu_p),
            }
        )

    return Success(pl.DataFrame(rows))


def compute_2d_feature_embedding(
    df: pl.DataFrame,
    feature_cols: Sequence[str],
    metric: KnnMetric = "correlation",
    seed: int = 42,
) -> tuple[list[float], list[float], list[float], list[float]]:
    """Compute deterministic 2D PCA and UMAP coordinates on numeric signature/CLR features."""
    valid_cols = [c for c in feature_cols if c in df.columns]
    if df.is_empty() or not valid_cols:
        n = df.height
        zeros = [0.0] * n
        return zeros, zeros, zeros, zeros

    mat = df.select(valid_cols).to_numpy().astype(np.float64)
    n_comp = min(10, mat.shape[1], mat.shape[0])
    if n_comp < 1:
        zeros = [0.0] * mat.shape[0]
        return zeros, zeros, zeros, zeros

    pcs = PCA(n_components=n_comp, random_state=seed).fit_transform(mat)
    pc1 = pcs[:, 0].tolist()
    pc2 = (pcs[:, 1] if pcs.shape[1] > 1 else np.zeros(mat.shape[0], dtype=np.float64)).tolist()

    if mat.shape[0] < 6:
        return pc1, pc2, pc1, pc2

    try:
        import anndata as ad  # type: ignore[import-untyped]
        import scanpy as sc  # type: ignore[import-untyped]

        tmp_ad = ad.AnnData(X=mat.astype(np.float32))
        sc.pp.neighbors(
            tmp_ad,
            n_neighbors=min(20, max(3, mat.shape[0] - 1)),
            use_rep="X",
            metric=metric,
            random_state=seed,
        )
        sc.tl.umap(tmp_ad, random_state=seed)
        coords = np.asarray(tmp_ad.obsm["X_umap"], dtype=np.float64)
        return pc1, pc2, coords[:, 0].tolist(), coords[:, 1].tolist()
    except ImportError:
        return pc1, pc2, pc1, pc2


def benchmark_bagaev_vs_deconv_cv(
    df: pl.DataFrame,
    state_cols: Sequence[str],
    signature_cols: Sequence[str] = BAGAEV_MFP_SIGNATURE_ORDER,
    subtype_col: str = "bagaev_subtype",
    response_col: str = "response_binary",
    group_col: Maybe[str] = Some("patient_id"),
    cohort_col: Maybe[str] = Some("cohort"),
    n_splits: int = 5,
    n_repeats: int = 5,
    seed: int = 42,
) -> Result[tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame], str]:
    """Run leakage-free patient-grouped 5x5 CV comparing Bagaev subtypes, 29-Fges & raw deconvolution."""
    req = ("sample_id", subtype_col, response_col)
    for col in req:
        if col not in df.columns:
            return Failure(f"Missing required column '{col}' for CV benchmark.")

    fges_cols = [c for c in signature_cols if c in df.columns]
    meta_cols = [c for c in state_cols if c in df.columns]
    if len(fges_cols) < 8 or len(meta_cols) < 4:
        return Failure("Insufficient Fges or deconvolution columns for CV benchmark.")

    clean = df.filter(
        pl.col(response_col).cast(pl.Float64).is_in([0.0, 1.0])
        & pl.col(subtype_col).is_not_null()
    )
    if clean.height < max(20, n_splits * 4):
        return Failure(f"Insufficient evaluable samples for {n_splits}-fold CV: {clean.height}.")

    y = clean[response_col].to_numpy().astype(np.int64)
    if int(np.sum(y == 1)) < n_splits or int(np.sum(y == 0)) < n_splits:
        return Failure("Insufficient class counts for stratified CV.")

    subtypes = clean[subtype_col].cast(pl.String).to_list()
    x_onehot = np.column_stack(
        [[1.0 if s == st_name else 0.0 for s in subtypes] for st_name in CANONICAL_SUBTYPE_ORDER]
    )
    x_fges = clean.select(fges_cols).to_numpy().astype(np.float64)
    x_meta = np.arcsin(np.sqrt(np.clip(clean.select(meta_cols).to_numpy().astype(np.float64), 0.0, 1.0)))
    x_comb = np.hstack([x_fges, x_meta])

    use_groups = False
    group_arr = np.arange(clean.height, dtype=object)
    match group_col:
        case Some(g_col) if g_col in clean.columns:
            coh_prefix = (
                clean[cohort_col.unwrap()].cast(pl.String)
                if isinstance(cohort_col, Some) and cohort_col.unwrap() in clean.columns
                else pl.lit("Coh")
            )
            group_keys = (
                clean.select([coh_prefix.alias("_c"), pl.col(g_col).cast(pl.String).alias("_g")])
                .select(pl.format("{}__{}", pl.col("_c"), pl.col("_g")).alias("_gk"))["_gk"]
                .to_numpy()
            )
            if len(np.unique(group_keys)) >= n_splits:
                group_arr = group_keys
                use_groups = True
        case _:
            pass

    fold_splits: list[tuple[npt.NDArray[np.int64], npt.NDArray[np.int64]]] = []
    if use_groups:
        for rep in range(n_repeats):
            sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed + rep)
            for tr_idx, te_idx in sgkf.split(x_fges, y, groups=group_arr):
                fold_splits.append((np.asarray(tr_idx, dtype=np.int64), np.asarray(te_idx, dtype=np.int64)))
    else:
        rskf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=seed)
        for tr_idx, te_idx in rskf.split(x_fges, y):
            fold_splits.append((np.asarray(tr_idx, dtype=np.int64), np.asarray(te_idx, dtype=np.int64)))

    model_specs: tuple[tuple[str, npt.NDArray[np.float64]], ...] = (
        ("Bagaev_4Subtype_Categorical", x_onehot),
        ("Bagaev_29Fges_L2_Logistic", x_fges),
        ("Deconv_RawMeta42_Fractions_L2", x_meta),
        ("Combined_Bagaev29_Plus_RawMeta42", x_comb),
    )

    pan_rows: list[dict[str, str | float | int]] = []
    oof_probs: dict[str, npt.NDArray[np.float64]] = {name: np.zeros(clean.height, dtype=np.float64) for name, _ in model_specs}
    oof_counts: dict[str, npt.NDArray[np.int64]] = {name: np.zeros(clean.height, dtype=np.int64) for name, _ in model_specs}

    for name, x_mat in model_specs:
        aucs: list[float] = []
        praucs: list[float] = []
        for tr_idx, te_idx in fold_splits:
            if len(np.unique(y[tr_idx])) < 2 or len(np.unique(y[te_idx])) < 2:
                continue
            scaler = StandardScaler()
            x_tr = scaler.fit_transform(x_mat[tr_idx])
            x_te = scaler.transform(x_mat[te_idx])
            clf = LogisticRegression(C=0.25, max_iter=1000, class_weight="balanced", random_state=seed)
            clf.fit(x_tr, y[tr_idx])
            prob = np.asarray(clf.predict_proba(x_te)[:, 1], dtype=np.float64)
            aucs.append(float(roc_auc_score(y[te_idx], prob)))
            praucs.append(float(average_precision_score(y[te_idx], prob)))
            oof_probs[name][te_idx] += prob
            oof_counts[name][te_idx] += 1

        if not aucs:
            return Failure(f"CV fold evaluation produced no valid folds for model '{name}'.")
        pan_rows.append(
            {
                "scope": "PanCancer_GroupedCV" if use_groups else "PanCancer_CV",
                "model": name,
                "n_samples": int(clean.height),
                "n_unique_patients": int(len(np.unique(group_arr))),
                "roc_auc_mean": float(np.mean(aucs)),
                "roc_auc_sd": float(np.std(aucs)),
                "pr_auc_mean": float(np.mean(praucs)),
                "pr_auc_sd": float(np.std(praucs)),
            }
        )

    def _to_oof_z(raw_sum: npt.NDArray[np.float64], cnts: npt.NDArray[np.int64]) -> npt.NDArray[np.float64]:
        p_mean = np.clip(raw_sum / np.maximum(cnts, 1), 1e-6, 1.0 - 1e-6)
        logits = np.log(p_mean / (1.0 - p_mean))
        sd = float(np.std(logits))
        return np.zeros_like(logits) if sd <= 1e-12 else (logits - float(np.mean(logits))) / sd

    oof_meta_z = _to_oof_z(oof_probs["Deconv_RawMeta42_Fractions_L2"], oof_counts["Deconv_RawMeta42_Fractions_L2"])
    oof_b29_z = _to_oof_z(oof_probs["Bagaev_29Fges_L2_Logistic"], oof_counts["Bagaev_29Fges_L2_Logistic"])
    oof_comb_z = _to_oof_z(oof_probs["Combined_Bagaev29_Plus_RawMeta42"], oof_counts["Combined_Bagaev29_Plus_RawMeta42"])

    oof_df = pl.DataFrame(
        {
            "sample_id": clean["sample_id"].cast(pl.String).to_list(),
            "oof_raw_meta42_prob": (
                oof_probs["Deconv_RawMeta42_Fractions_L2"] / np.maximum(oof_counts["Deconv_RawMeta42_Fractions_L2"], 1)
            ).tolist(),
            "oof_raw_meta42_z": oof_meta_z.tolist(),
            "oof_bagaev29_z": oof_b29_z.tolist(),
            "oof_combined_z": oof_comb_z.tolist(),
        }
    )

    per_cohort_rows: list[dict[str, str | float | int]] = []
    match cohort_col:
        case Some(c_col) if c_col in clean.columns:
            c_arr = clean[c_col].cast(pl.String).to_numpy()
            for cid in sorted(str(x) for x in np.unique(c_arr)):
                mask = c_arr == cid
                yc = y[mask]
                if int(np.sum(mask)) < 20 or int(np.sum(yc == 1)) < 5 or int(np.sum(yc == 0)) < 5:
                    continue
                for m_name in (
                    "Bagaev_4Subtype_Categorical",
                    "Bagaev_29Fges_L2_Logistic",
                    "Deconv_RawMeta42_Fractions_L2",
                    "Combined_Bagaev29_Plus_RawMeta42",
                ):
                    p_coh = oof_probs[m_name][mask] / np.maximum(oof_counts[m_name][mask], 1)
                    per_cohort_rows.append(
                        {
                            "cohort": cid,
                            "model": m_name,
                            "n_samples": int(np.sum(mask)),
                            "n_responders": int(np.sum(yc == 1)),
                            "oof_roc_auc": float(roc_auc_score(yc, p_coh)),
                            "oof_pr_auc": float(average_precision_score(yc, p_coh)),
                        }
                    )
        case _:
            pass

    per_cohort_df = (
        pl.DataFrame(per_cohort_rows)
        if per_cohort_rows
        else pl.DataFrame(
            schema={
                "cohort": pl.String,
                "model": pl.String,
                "n_samples": pl.Int64,
                "n_responders": pl.Int64,
                "oof_roc_auc": pl.Float64,
                "oof_pr_auc": pl.Float64,
            }
        )
    )
    return Success((pl.DataFrame(pan_rows), per_cohort_df, oof_df))


def evaluate_fges_correlation_conservation(
    query_scaled_df: pl.DataFrame,
    ref_scaled_df: pl.DataFrame,
    signature_cols: Sequence[str] = BAGAEV_MFP_SIGNATURE_ORDER,
) -> Result[dict[str, float], str]:
    """Quantify conservation of the 29x29 Fges inter-signature correlation matrix vs TCGA."""
    shared = [c for c in signature_cols if c in query_scaled_df.columns and c in ref_scaled_df.columns]
    if len(shared) < 8 or query_scaled_df.height < 4 or ref_scaled_df.height < 4:
        return Failure("Insufficient samples or shared signatures for 29x29 correlation conservation.")

    q_mat = query_scaled_df.select(shared).to_numpy().astype(np.float64).T
    r_mat = ref_scaled_df.select(shared).to_numpy().astype(np.float64).T
    c_q = _pairwise_pearson_correlation(q_mat, q_mat)
    c_r = _pairwise_pearson_correlation(r_mat, r_mat)

    tri_i, tri_j = np.triu_indices(len(shared), k=1)
    v_q = c_q[tri_i, tri_j]
    v_r = c_r[tri_i, tri_j]

    pearson_r = float(stats.pearsonr(v_q, v_r).statistic) if np.std(v_q) > 1e-12 and np.std(v_r) > 1e-12 else 0.0
    spearman_rho = float(stats.spearmanr(v_q, v_r).statistic) if np.std(v_q) > 1e-12 and np.std(v_r) > 1e-12 else 0.0
    frob_rmse = float(np.sqrt(np.mean((v_q - v_r) ** 2)))

    idx_map = {s: idx for idx, s in enumerate(shared)}
    within_vals: list[float] = []
    between_vals: list[float] = []
    sig_to_grp = {
        s: grp_name
        for grp_name, sigs in BAGAEV_MFP_FUNCTIONAL_GROUPS.items()
        for s in sigs
    }
    for i, j in zip(tri_i.tolist(), tri_j.tolist(), strict=True):
        s_i, s_j = shared[i], shared[j]
        if sig_to_grp.get(s_i) == sig_to_grp.get(s_j):
            within_vals.append(float(c_q[idx_map[s_i], idx_map[s_j]]))
        else:
            between_vals.append(float(c_q[idx_map[s_i], idx_map[s_j]]))

    mean_within = float(np.mean(within_vals)) if within_vals else 0.0
    mean_between = float(np.mean(between_vals)) if between_vals else 0.0

    return Success(
        {
            "n_signatures": float(len(shared)),
            "offdiag_pearson_r_vs_tcga": pearson_r,
            "offdiag_spearman_rho_vs_tcga": spearman_rho,
            "offdiag_rmse_vs_tcga": frob_rmse,
            "mean_within_block_corr": mean_within,
            "mean_between_block_corr": mean_between,
            "block_coherence_delta": mean_within - mean_between,
        }
    )


__all__ = [
    "BagaevKnnConfig",
    "BagaevLouvainConfig",
    "BagaevLouvainResult",
    "BagaevSubtype",
    "CANONICAL_SUBTYPE_ORDER",
    "DeconvClusterAlgorithm",
    "DeconvClusterConfig",
    "DeconvClusterResult",
    "DeconvTransformMethod",
    "IATLAS_TO_TCGA_PROJECTS",
    "INDICATION_TO_TCGA_PROJECTS",
    "KnnMetric",
    "analyze_subtype_deconvolution_enrichment",
    "assign_clusters_to_bagaev_subtypes",
    "benchmark_bagaev_vs_deconv_cv",
    "classify_bagaev_knn",
    "cluster_bagaev_louvain",
    "cluster_deconvoluted_samples",
    "compare_partition_concordance",
    "compute_2d_feature_embedding",
    "compute_cmh_common_odds_ratio",
    "compute_cmh_k_by_2_p_value",
    "compute_deconv_subtype_correlations",
    "compute_state_vs_fges_correlations",
    "evaluate_fges_correlation_conservation",
    "evaluate_subtype_response_associations",
    "evaluate_within_subtype_deconv_predictions",
    "transform_compositional_fractions",
    "wilson_ci",
]
