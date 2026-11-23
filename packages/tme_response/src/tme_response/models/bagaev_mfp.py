"""Bagaev et al. (2021, Cancer Cell) Pan-Cancer Molecular Functional Portrait (MFP) TME Subtyping.

Implements:
1. Supervised K=35 Nearest-Neighbor classification against labeled TCGA reference portraits
   (indication-matched or pan-cancer) in 29 median/MAD-scaled Fges signature space.
2. Unsupervised Louvain community detection on thresholded sample Pearson correlation graphs
   with Silhouette, Calinski-Harabasz, and Davies-Bouldin threshold selection and canonical
   quadrant matching onto the 4 conserved TME subtypes (IE/F, IE, F, D).
3. Statistical evaluation of subtype objective response rates (ORR, Wilson CIs, Fisher exact,
   chi-square, and cohort-stratified Cochran-Mantel-Haenszel common odds ratios).
4. Cross-modality integration with single-cell deconvolution fractions, out-of-fold response
   predictions (within-subtype ROC-AUC / PR-AUC), and 29x29 Fges correlation conservation.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Final, Sequence
import numpy as np
import numpy.typing as npt
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import optimize, stats
from sklearn.metrics import (
    adjusted_rand_score,
    average_precision_score,
    calinski_harabasz_score,
    cohen_kappa_score,
    davies_bouldin_score,
    normalized_mutual_info_score,
    roc_auc_score,
    silhouette_score,
)

from tme_datasets import (
    BAGAEV_MFP_FUNCTIONAL_GROUPS,
    BAGAEV_MFP_SIGNATURE_ORDER,
)


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
}


@dataclass(frozen=True)
class BagaevKnnConfig:
    """Configuration for Bagaev supervised K-NN TCGA classification."""

    k_neighbors: int = 35
    metric: str = "correlation"  # 'correlation' or 'euclidean'
    distance_weighted: bool = False


@dataclass(frozen=True)
class BagaevLouvainConfig:
    """Configuration for Bagaev correlation-graph Louvain clustering."""

    thresholds: tuple[float, ...] = (0.30, 0.35, 0.40, 0.45, 0.50, 0.51, 0.55, 0.60)
    min_cluster_fraction: float = 0.05
    target_n_clusters: int = 4
    seed: int = 42


@dataclass(frozen=True)
class BagaevLouvainResult:
    """Immutable result of Bagaev correlation-graph clustering and threshold sweep."""

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
    out: dict[str, npt.NDArray[np.float64]] = {}
    for block_name, block_sigs in BAGAEV_MFP_FUNCTIONAL_GROUPS.items():
        idxs = [col_idx[s] for s in block_sigs if s in col_idx]
        out[block_name] = (
            np.mean(mat[:, idxs], axis=1)
            if idxs
            else np.zeros(mat.shape[0], dtype=np.float64)
        )
    return out


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
    return np.clip(corr, -1.0, 1.0)


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

    if config.metric == "correlation":
        corr = _pairwise_pearson_correlation(q_mat, r_mat)
        dists = 1.0 - corr
    elif config.metric == "euclidean":
        diff = q_mat[:, np.newaxis, :] - r_mat[np.newaxis, :, :]
        dists = np.linalg.norm(diff, axis=2)
    else:
        return Failure(f"Unsupported distance metric '{config.metric}'.")

    nn_idx = np.argpartition(dists, kth=k_eff - 1, axis=1)[:, :k_eff]
    n_query = q_mat.shape[0]
    probs = np.zeros((n_query, len(CANONICAL_SUBTYPE_ORDER)), dtype=np.float64)

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
            probs[i, s_idx] = float(np.sum(weights[neigh_labels == subtype])) / w_sum

    best_idx = np.argmax(probs, axis=1)
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
) -> npt.NDArray[np.int64]:
    """Deterministic Blondel modularity-maximizing greedy Louvain community detection."""
    n = adj.shape[0]
    if n == 0:
        return np.array([], dtype=np.int64)
    w_clean = np.maximum(adj, 0.0).copy()
    np.fill_diagonal(w_clean, 0.0)
    two_m = float(np.sum(w_clean))
    if two_m <= 1e-12:
        return np.arange(n, dtype=np.int64)

    deg = np.sum(w_clean, axis=1)
    comm = np.arange(n, dtype=np.int64)
    rng = np.random.default_rng(seed)
    order = rng.permutation(n)

    tot_weight: dict[int, float] = {int(i): float(deg[i]) for i in range(n)}

    for _ in range(max_passes):
        moved = False
        for u_raw in order:
            u = int(u_raw)
            c_old = int(comm[u])
            k_u = float(deg[u])
            if k_u <= 1e-12:
                continue

            tot_weight[c_old] = tot_weight.get(c_old, 0.0) - k_u
            neigh_mask = w_clean[u] > 0.0
            neigh_nodes = np.flatnonzero(neigh_mask)
            cand_weights: dict[int, float] = {c_old: 0.0}
            for v_raw in neigh_nodes:
                v = int(v_raw)
                if v != u:
                    c_v = int(comm[v])
                    cand_weights[c_v] = cand_weights.get(c_v, 0.0) + float(w_clean[u, v])

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
    return relabeled.astype(np.int64)


def _consolidate_to_four_clusters(
    raw_labels: npt.NDArray[np.int64],
    mat: npt.NDArray[np.float64],
    corr: npt.NDArray[np.float64],
    min_cluster_fraction: float,
    target_k: int = 4,
) -> npt.NDArray[np.int64]:
    """Consolidate Louvain communities into `target_k` major clusters (>= min_cluster_fraction)."""
    n = mat.shape[0]
    min_size = max(2, int(np.floor(n * min_cluster_fraction)))
    u_ids, counts = np.unique(raw_labels, return_counts=True)
    order = np.argsort(-counts)
    u_sorted = u_ids[order]
    c_sorted = counts[order]

    major_ids = [int(cid) for cid, cnt in zip(u_sorted, c_sorted, strict=True) if int(cnt) >= min_size]
    if len(major_ids) > target_k:
        major_ids = major_ids[:target_k]

    if len(major_ids) < target_k:
        # Deterministic complete-linkage hierarchical clustering on correlation distance (1 - R)
        from sklearn.cluster import AgglomerativeClustering

        dist_mat = np.clip(1.0 - corr, 0.0, 2.0)
        np.fill_diagonal(dist_mat, 0.0)
        agg = AgglomerativeClustering(
            n_clusters=min(target_k, n),
            metric="precomputed",
            linkage="average",
        )
        return agg.fit_predict(dist_mat).astype(np.int64)

    # Reassign minor/singleton nodes to the major cluster with highest mean correlation
    centroids = np.vstack(
        [np.mean(mat[raw_labels == cid], axis=0) for cid in major_ids]
    )
    sim_to_major = _pairwise_pearson_correlation(mat, centroids)
    assigned = np.zeros(n, dtype=np.int64)
    major_lookup = {cid: idx for idx, cid in enumerate(major_ids)}
    for i in range(n):
        orig = int(raw_labels[i])
        if orig in major_lookup:
            assigned[i] = major_lookup[orig]
        else:
            assigned[i] = int(np.argmax(sim_to_major[i]))
    return assigned


def assign_clusters_to_bagaev_subtypes(
    cluster_ids: npt.NDArray[np.int64],
    mat: npt.NDArray[np.float64],
    sig_cols: Sequence[str],
) -> tuple[list[str], dict[int, str]]:
    """Map unsupervised clusters onto {IE/F, IE, F, D} via biological quadrant centroids."""
    blocks = _compute_functional_block_means(mat, sig_cols)
    imm = blocks["Antitumor_Immune"]
    fib = blocks["Angiogenesis_Fibroblasts"]

    unique_clusters = sorted(int(c) for c in np.unique(cluster_ids))
    # Standardize cluster centroids across clusters so relative high/low quadrants align cleanly
    raw_c_imm = np.array([float(np.mean(imm[cluster_ids == c])) for c in unique_clusters])
    raw_c_fib = np.array([float(np.mean(fib[cluster_ids == c])) for c in unique_clusters])
    std_imm = raw_c_imm - float(np.mean(raw_c_imm))
    std_fib = raw_c_fib - float(np.mean(raw_c_fib))

    # Canonical target vectors in (Immune, Fibrotic) plane:
    # IE: (+1, -1), IE/F: (+1, +1), F: (-1, +1), D: (-1, -1)
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

    if len(unique_clusters) == 4:
        # Maximize dot product with canonical quadrant directions via Hungarian algorithm
        cost = -(centroids_2d @ targets.T)
        row_ind, col_ind = optimize.linear_sum_assignment(cost)
        mapping = {
            unique_clusters[int(r)]: CANONICAL_SUBTYPE_ORDER[int(c)]
            for r, c in zip(row_ind, col_ind, strict=True)
        }
    else:
        mapping = {}
        for idx, cid in enumerate(unique_clusters):
            scores = centroids_2d[idx] @ targets.T
            mapping[cid] = CANONICAL_SUBTYPE_ORDER[int(np.argmax(scores))]

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

        if n_eff >= 2 and n_eff < mat.shape[0]:
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

    # Composite rank score across Silhouette (higher better), CH (higher better), DB (lower better)
    sil_r = stats.rankdata(eval_df["silhouette_corr"].to_numpy())
    ch_r = stats.rankdata(eval_df["calinski_harabasz"].to_numpy())
    db_r = stats.rankdata(-eval_df["davies_bouldin"].to_numpy())
    composite = sil_r + ch_r + db_r
    best_row_idx = int(np.argmax(composite))
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

        r1 = a + b
        r0 = c + d
        c1 = a + c
        c0 = b + d
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

    # 4-way chi-square test of independence
    present_subtypes = [s for s in CANONICAL_SUBTYPE_ORDER if np.any(s_all == s)]
    for extra in sorted(set(s_all.tolist()) - set(present_subtypes)):
        present_subtypes.append(str(extra))

    table_rows = [
        [int(np.sum((s_all == st_name) & y_all)), int(np.sum((s_all == st_name) & (~y_all)))]
        for st_name in present_subtypes
    ]
    valid_rows = [r for r in table_rows if (r[0] + r[1]) > 0]
    if len(valid_rows) >= 2 and sum(r[0] for r in valid_rows) > 0 and sum(r[1] for r in valid_rows) > 0:
        _, chi2_p, _, _ = stats.chi2_contingency(np.asarray(table_rows, dtype=np.float64))
    else:
        chi2_p = 1.0

    # Favorable (IE + IE/F) vs Unfavorable (F + D)
    fav_mask = np.isin(s_all, ["IE", "IE/F"])
    a_f = int(np.sum(fav_mask & y_all))
    b_f = int(np.sum(fav_mask & (~y_all)))
    c_f = int(np.sum((~fav_mask) & y_all))
    d_f = int(np.sum((~fav_mask) & (~y_all)))
    fav_or, fav_fisher_p = stats.fisher_exact([[a_f, b_f], [c_f, d_f]])

    match cohort_col:
        case Some(c_col) if c_col in valid.columns:
            with_fav = valid.with_columns(
                pl.col(subtype_col).cast(pl.String).is_in(["IE", "IE/F"]).alias("_is_favorable")
            )
            cmh_or, cmh_p = compute_cmh_common_odds_ratio(
                with_fav, exposure_col="_is_favorable", outcome_col=response_col, stratum_col=c_col
            )
        case _:
            cmh_or, cmh_p = float(fav_or), float(fav_fisher_p)

    rows: list[dict[str, str | int | float]] = []
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
            st_or, st_p = stats.fisher_exact([[n_resp, n_non], [rest_resp, rest_non]])
        else:
            st_or, st_p = 1.0, 1.0

        rows.append(
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

    # Benjamini-Hochberg FDR adjustment
    p_arr = np.asarray(raw_pvals, dtype=np.float64)
    m = len(p_arr)
    order = np.argsort(p_arr)
    q_arr = np.empty(m, dtype=np.float64)
    prev_q = 1.0
    for rank_rev in range(m - 1, -1, -1):
        idx = int(order[rank_rev])
        q_val = min(prev_q, float(p_arr[idx]) * float(m) / float(rank_rev + 1))
        q_arr[idx] = min(1.0, q_val)
        prev_q = q_arr[idx]

    out_df = pl.DataFrame(rows).with_columns(pl.Series(name="fdr_q", values=q_arr)).sort("kruskal_p")
    return Success(out_df)


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


def evaluate_fges_correlation_conservation(
    query_scaled_df: pl.DataFrame,
    ref_scaled_df: pl.DataFrame,
    signature_cols: Sequence[str] = BAGAEV_MFP_SIGNATURE_ORDER,
) -> Result[dict[str, float], str]:
    """Quantify conservation of the 29x29 Fges inter-signature correlation matrix vs TCGA."""
    shared = [c for c in signature_cols if c in query_scaled_df.columns and c in ref_scaled_df.columns]
    if len(shared) < 8 or query_scaled_df.height < 4 or ref_scaled_df.height < 4:
        return Failure("Insufficient samples or shared signatures for 29x29 correlation conservation.")

    q_mat = query_scaled_df.select(shared).to_numpy().astype(np.float64).T  # G x N_q
    r_mat = ref_scaled_df.select(shared).to_numpy().astype(np.float64).T    # G x N_r
    c_q = _pairwise_pearson_correlation(q_mat, q_mat)
    c_r = _pairwise_pearson_correlation(r_mat, r_mat)

    tri_i, tri_j = np.triu_indices(len(shared), k=1)
    v_q = c_q[tri_i, tri_j]
    v_r = c_r[tri_i, tri_j]

    pearson_r = float(stats.pearsonr(v_q, v_r).statistic) if np.std(v_q) > 1e-12 and np.std(v_r) > 1e-12 else 0.0
    spearman_rho = float(stats.spearmanr(v_q, v_r).statistic) if np.std(v_q) > 1e-12 and np.std(v_r) > 1e-12 else 0.0
    frob_rmse = float(np.sqrt(np.mean((v_q - v_r) ** 2)))

    # Within-block correlation coherence in query
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
    "IATLAS_TO_TCGA_PROJECTS",
    "INDICATION_TO_TCGA_PROJECTS",
    "analyze_subtype_deconvolution_enrichment",
    "assign_clusters_to_bagaev_subtypes",
    "classify_bagaev_knn",
    "cluster_bagaev_louvain",
    "compare_partition_concordance",
    "compute_cmh_common_odds_ratio",
    "evaluate_fges_correlation_conservation",
    "evaluate_subtype_response_associations",
    "evaluate_within_subtype_deconv_predictions",
    "wilson_ci",
]
