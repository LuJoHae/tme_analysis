"""Pure functional core for Meta42 <-> Coulton_TAM_Atlas substate matching and ICI outcome synthesis.

Strictly pure module (.agents/rules/code-style-guide.md):
- Zero file/network/database IO, zero mutation of inputs, zero `None` / `del` / `global` / `nonlocal`
- Exhaustive pattern matching with `typing.assert_never`
- Vectorized NumPy/SciPy sparse algebra and Polars tabular pipelines
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import math
from typing import Final, assert_never

import numpy as np
import numpy.typing as npt
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import sparse as sp
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist
from scipy.stats import norm, rankdata

from .models import (
    ALL_MYELOID_META42_IDS,
    AssaySource,
    BiclusterResult,
    BoolVector,
    BulkIcbCohortBundle,
    COULTON_LITERATURE_CATALOG,
    COULTON_ORDERED_CLUSTERS,
    CORE_TAM_META42_IDS,
    ComparisonConfig,
    ComparisonOutputs,
    CoultonCentroidAccumulator,
    CoultonCentroidBundle,
    CoultonClusterRole,
    DendrogramBracket,
    FloatMatrix,
    FloatVector,
    IciOutcomeDirection,
    IntVector,
    MACROPHAGE_CONTINUUM_META42_IDS,
    META42_BIO_LABELS,
    Meta42MacrophageTier,
    Meta42ReferenceBundle,
    OutcomeConcordanceClass,
    PAN_MACROPHAGE_CORE_GENES,
    TissueStratum,
    classify_meta42_compartment,
    classify_meta42_macrophage_tier,
)

_EXCLUDED_MARKER_PREFIXES: Final[tuple[str, ...]] = (
    "ENSG",
    "LINC",
    "MIR",
    "LOC",
    "RNU",
    "SNOR",
    "MT-",
    "RP11",
    "AC0",
    "AL0",
    "AP0",
    "RPS",
    "RPL",
)

_LYMPHOID_SPILLOVER_GENES: Final[frozenset[str]] = frozenset(
    {
        "CD3D",
        "CD3E",
        "CD3G",
        "CD8A",
        "CD8B",
        "TRAC",
        "TRBC1",
        "TRBC2",
        "NKG7",
        "GZMA",
        "GZMB",
        "GZMK",
        "PRF1",
        "CD79A",
        "CD79B",
        "MS4A1",
        "IGHG1",
        "IGHG2",
        "IGHG3",
        "IGHG4",
        "IGKC",
        "IGLC2",
        "JCHAIN",
        "MZB1",
    }
)


# =============================================================================
# 1. Pure Out-of-Core Sparse CSR Batch Centroid Reduction
# =============================================================================


def create_empty_coulton_accumulator(
    assay: AssaySource,
    cluster_labels: Sequence[str],
    gene_symbols: Sequence[str],
) -> CoultonCentroidAccumulator:
    """Create a zero-initialized immutable accumulator for out-of-core cluster aggregation."""
    k_clusters = len(cluster_labels)
    g_genes = len(gene_symbols)
    return CoultonCentroidAccumulator(
        assay=assay,
        cluster_labels=tuple(cluster_labels),
        gene_symbols=tuple(gene_symbols),
        all_cell_counts=np.zeros(k_clusters, dtype=np.int64),
        tumor_cell_counts=np.zeros(k_clusters, dtype=np.int64),
        all_linear_sum=np.zeros((k_clusters, g_genes), dtype=np.float64),
        tumor_linear_sum=np.zeros((k_clusters, g_genes), dtype=np.float64),
        all_log_sum=np.zeros((k_clusters, g_genes), dtype=np.float64),
        tumor_log_sum=np.zeros((k_clusters, g_genes), dtype=np.float64),
        all_detect_counts=np.zeros((k_clusters, g_genes), dtype=np.float64),
        tumor_detect_counts=np.zeros((k_clusters, g_genes), dtype=np.float64),
    )


def accumulate_coulton_batch(
    acc: CoultonCentroidAccumulator,
    batch_csr: sp.csr_matrix,
    cluster_codes: IntVector,
    is_tumor_mask: BoolVector,
) -> Result[CoultonCentroidAccumulator, str]:
    """Fold one sparse CSR batch of cells into cluster sums without densifying zeros."""
    n_batch, n_genes = batch_csr.shape
    k_clusters = len(acc.cluster_labels)
    if n_genes != len(acc.gene_symbols):
        return Failure(
            f"Batch gene count {n_genes} does not match accumulator {len(acc.gene_symbols)}"
        )
    if len(cluster_codes) != n_batch or len(is_tumor_mask) != n_batch:
        return Failure("Batch metadata vector lengths do not match CSR row count")

    valid_all = (cluster_codes >= 0) & (cluster_codes < k_clusters)
    valid_tum = valid_all & is_tumor_mask

    cell_idx_all = np.flatnonzero(valid_all)
    clust_idx_all = cluster_codes[cell_idx_all]
    ind_all = sp.csr_matrix(
        (
            np.ones(len(cell_idx_all), dtype=np.float64),
            (clust_idx_all, cell_idx_all),
        ),
        shape=(k_clusters, n_batch),
    )

    cell_idx_tum = np.flatnonzero(valid_tum)
    clust_idx_tum = cluster_codes[cell_idx_tum]
    ind_tum = sp.csr_matrix(
        (
            np.ones(len(cell_idx_tum), dtype=np.float64),
            (clust_idx_tum, cell_idx_tum),
        ),
        shape=(k_clusters, n_batch),
    )

    log_csr = batch_csr.tocsr().astype(np.float64, copy=False)
    lin_data = np.expm1(np.clip(log_csr.data, 0.0, 25.0))
    lin_csr = sp.csr_matrix(
        (lin_data, log_csr.indices, log_csr.indptr),
        shape=log_csr.shape,
    )
    det_data = (log_csr.data > 0.0).astype(np.float64)
    det_csr = sp.csr_matrix(
        (det_data, log_csr.indices, log_csr.indptr),
        shape=log_csr.shape,
    )

    add_all_counts = np.bincount(clust_idx_all, minlength=k_clusters).astype(np.int64)
    add_tum_counts = np.bincount(clust_idx_tum, minlength=k_clusters).astype(np.int64)

    add_all_lin = np.asarray((ind_all @ lin_csr).toarray(), dtype=np.float64)
    add_tum_lin = np.asarray((ind_tum @ lin_csr).toarray(), dtype=np.float64)
    add_all_log = np.asarray((ind_all @ log_csr).toarray(), dtype=np.float64)
    add_tum_log = np.asarray((ind_tum @ log_csr).toarray(), dtype=np.float64)
    add_all_det = np.asarray((ind_all @ det_csr).toarray(), dtype=np.float64)
    add_tum_det = np.asarray((ind_tum @ det_csr).toarray(), dtype=np.float64)

    return Success(
        CoultonCentroidAccumulator(
            assay=acc.assay,
            cluster_labels=acc.cluster_labels,
            gene_symbols=acc.gene_symbols,
            all_cell_counts=acc.all_cell_counts + add_all_counts,
            tumor_cell_counts=acc.tumor_cell_counts + add_tum_counts,
            all_linear_sum=acc.all_linear_sum + add_all_lin,
            tumor_linear_sum=acc.tumor_linear_sum + add_tum_lin,
            all_log_sum=acc.all_log_sum + add_all_log,
            tumor_log_sum=acc.tumor_log_sum + add_tum_log,
            all_detect_counts=np.asarray(acc.all_detect_counts, dtype=np.float64) + add_all_det,
            tumor_detect_counts=np.asarray(acc.tumor_detect_counts, dtype=np.float64) + add_tum_det,
        )
    )


def finalize_coulton_centroids(acc: CoultonCentroidAccumulator) -> CoultonCentroidBundle:
    """Convert accumulated sums into exact per-cell linear means, log1p means, and detection %."""
    denom_all = np.maximum(acc.all_cell_counts.reshape((-1, 1)).astype(np.float64), 1.0)
    denom_tum = np.maximum(acc.tumor_cell_counts.reshape((-1, 1)).astype(np.float64), 1.0)

    has_all = (acc.all_cell_counts > 0).reshape((-1, 1))
    has_tum = (acc.tumor_cell_counts > 0).reshape((-1, 1))

    all_det_mat = np.asarray(acc.all_detect_counts, dtype=np.float64)
    tum_det_mat = np.asarray(acc.tumor_detect_counts, dtype=np.float64)

    return CoultonCentroidBundle(
        assay=acc.assay,
        cluster_labels=acc.cluster_labels,
        gene_symbols=acc.gene_symbols,
        all_cell_counts=acc.all_cell_counts,
        tumor_cell_counts=acc.tumor_cell_counts,
        all_mean_linear=np.where(has_all, acc.all_linear_sum / denom_all, 0.0),
        tumor_mean_linear=np.where(has_tum, acc.tumor_linear_sum / denom_tum, 0.0),
        all_mean_log=np.where(has_all, acc.all_log_sum / denom_all, 0.0),
        tumor_mean_log=np.where(has_tum, acc.tumor_log_sum / denom_tum, 0.0),
        all_pct_expressing=np.where(has_all, (all_det_mat / denom_all) * 100.0, 0.0),
        tumor_pct_expressing=np.where(has_tum, (tum_det_mat / denom_tum) * 100.0, 0.0),
    )


def deduplicate_gene_columns(
    raw_symbols: Sequence[str],
    bundle: CoultonCentroidBundle,
) -> CoultonCentroidBundle:
    """Aggregate any duplicate HGNC gene symbols cleanly by linear addition."""
    unique_map: dict[str, list[int]] = {}
    ordered_syms: list[str] = []
    for idx, sym in enumerate(raw_symbols):
        cleaned = sym.strip()
        if not cleaned or cleaned.upper() in {"NAN", "NONE", "<NA>"}:
            continue
        if cleaned not in unique_map:
            unique_map[cleaned] = []
            ordered_syms.append(cleaned)
        unique_map[cleaned].append(idx)

    if len(ordered_syms) == len(raw_symbols) and all(len(v) == 1 for v in unique_map.values()):
        return CoultonCentroidBundle(
            assay=bundle.assay,
            cluster_labels=bundle.cluster_labels,
            gene_symbols=tuple(ordered_syms),
            all_cell_counts=bundle.all_cell_counts,
            tumor_cell_counts=bundle.tumor_cell_counts,
            all_mean_linear=bundle.all_mean_linear,
            tumor_mean_linear=bundle.tumor_mean_linear,
            all_mean_log=bundle.all_mean_log,
            tumor_mean_log=bundle.tumor_mean_log,
            all_pct_expressing=bundle.all_pct_expressing,
            tumor_pct_expressing=bundle.tumor_pct_expressing,
        )

    k_clusters = len(bundle.cluster_labels)
    g_out = len(ordered_syms)
    all_lin = np.zeros((k_clusters, g_out), dtype=np.float64)
    tum_lin = np.zeros((k_clusters, g_out), dtype=np.float64)
    all_log = np.zeros((k_clusters, g_out), dtype=np.float64)
    tum_log = np.zeros((k_clusters, g_out), dtype=np.float64)
    all_pct = np.zeros((k_clusters, g_out), dtype=np.float64)
    tum_pct = np.zeros((k_clusters, g_out), dtype=np.float64)

    for out_idx, sym in enumerate(ordered_syms):
        src_indices = unique_map[sym]
        if len(src_indices) == 1:
            s_i = src_indices[0]
            all_lin[:, out_idx] = bundle.all_mean_linear[:, s_i]
            tum_lin[:, out_idx] = bundle.tumor_mean_linear[:, s_i]
            all_log[:, out_idx] = bundle.all_mean_log[:, s_i]
            tum_log[:, out_idx] = bundle.tumor_mean_log[:, s_i]
            all_pct[:, out_idx] = bundle.all_pct_expressing[:, s_i]
            tum_pct[:, out_idx] = bundle.tumor_pct_expressing[:, s_i]
        else:
            all_lin[:, out_idx] = bundle.all_mean_linear[:, src_indices].sum(axis=1)
            tum_lin[:, out_idx] = bundle.tumor_mean_linear[:, src_indices].sum(axis=1)
            all_log[:, out_idx] = np.log1p(all_lin[:, out_idx])
            tum_log[:, out_idx] = np.log1p(tum_lin[:, out_idx])
            all_pct[:, out_idx] = bundle.all_pct_expressing[:, src_indices].max(axis=1)
            tum_pct[:, out_idx] = bundle.tumor_pct_expressing[:, src_indices].max(axis=1)

    return CoultonCentroidBundle(
        assay=bundle.assay,
        cluster_labels=bundle.cluster_labels,
        gene_symbols=tuple(ordered_syms),
        all_cell_counts=bundle.all_cell_counts,
        tumor_cell_counts=bundle.tumor_cell_counts,
        all_mean_linear=all_lin,
        tumor_mean_linear=tum_lin,
        all_mean_log=all_log,
        tumor_mean_log=tum_log,
        all_pct_expressing=all_pct,
        tumor_pct_expressing=tum_pct,
    )


def serialize_coulton_centroids_to_polars(
    bundles: Sequence[CoultonCentroidBundle],
) -> pl.DataFrame:
    """Serialize CoultonCentroidBundle sequence into a tidy, lossless Polars DataFrame."""
    frames: list[pl.DataFrame] = []
    for b in bundles:
        k_clusters = len(b.cluster_labels)
        g_genes = len(b.gene_symbols)
        cluster_idx_col = np.repeat(np.arange(k_clusters, dtype=np.int32), g_genes)
        cluster_lbl_col = np.repeat(np.asarray(b.cluster_labels, dtype=object), g_genes)
        all_cnt_col = np.repeat(b.all_cell_counts, g_genes)
        tum_cnt_col = np.repeat(b.tumor_cell_counts, g_genes)
        gene_idx_col = np.tile(np.arange(g_genes, dtype=np.int32), k_clusters)
        gene_sym_col = np.tile(np.asarray(b.gene_symbols, dtype=object), k_clusters)

        frames.append(
            pl.DataFrame(
                {
                    "assay": np.full(k_clusters * g_genes, b.assay.value, dtype=object),
                    "cluster_idx": cluster_idx_col,
                    "cluster_label": cluster_lbl_col,
                    "all_cell_count": all_cnt_col,
                    "tumor_cell_count": tum_cnt_col,
                    "gene_idx": gene_idx_col,
                    "gene_symbol": gene_sym_col,
                    "all_mean_linear": b.all_mean_linear.ravel().astype(np.float32),
                    "tumor_mean_linear": b.tumor_mean_linear.ravel().astype(np.float32),
                    "all_mean_log": b.all_mean_log.ravel().astype(np.float32),
                    "tumor_mean_log": b.tumor_mean_log.ravel().astype(np.float32),
                    "all_pct_expressing": b.all_pct_expressing.ravel().astype(np.float32),
                    "tumor_pct_expressing": b.tumor_pct_expressing.ravel().astype(np.float32),
                }
            )
        )
    return pl.concat(frames) if frames else pl.DataFrame()


def deserialize_coulton_centroids_from_polars(
    df: pl.DataFrame,
) -> Result[tuple[CoultonCentroidBundle, ...], str]:
    """Reconstruct CoultonCentroidBundle sequence from a cached Polars DataFrame."""
    if df.is_empty():
        return Failure("Cached centroid DataFrame is empty")

    bundles: list[CoultonCentroidBundle] = []
    for assay_val in df["assay"].unique(maintain_order=True).to_list():
        assay_enum = AssaySource.RNA_LOG1P if str(assay_val) == "RNA" else AssaySource.SCT_DATA
        sub = df.filter(pl.col("assay") == assay_val).sort(["cluster_idx", "gene_idx"])
        clust_df = sub.select(["cluster_idx", "cluster_label", "all_cell_count", "tumor_cell_count"]).unique(
            subset=["cluster_idx"], maintain_order=True
        ).sort("cluster_idx")
        gene_df = sub.select(["gene_idx", "gene_symbol"]).unique(
            subset=["gene_idx"], maintain_order=True
        ).sort("gene_idx")

        k_clusters = clust_df.height
        g_genes = gene_df.height
        if sub.height != k_clusters * g_genes:
            return Failure(
                f"Corrupted centroid cache for {assay_val}: expected {k_clusters * g_genes} rows, got {sub.height}"
            )

        bundles.append(
            CoultonCentroidBundle(
                assay=assay_enum,
                cluster_labels=tuple(clust_df["cluster_label"].cast(pl.Utf8).to_list()),
                gene_symbols=tuple(gene_df["gene_symbol"].cast(pl.Utf8).to_list()),
                all_cell_counts=clust_df["all_cell_count"].to_numpy().astype(np.int64),
                tumor_cell_counts=clust_df["tumor_cell_count"].to_numpy().astype(np.int64),
                all_mean_linear=sub["all_mean_linear"]
                .to_numpy()
                .astype(np.float64)
                .reshape((k_clusters, g_genes)),
                tumor_mean_linear=sub["tumor_mean_linear"]
                .to_numpy()
                .astype(np.float64)
                .reshape((k_clusters, g_genes)),
                all_mean_log=sub["all_mean_log"]
                .to_numpy()
                .astype(np.float64)
                .reshape((k_clusters, g_genes)),
                tumor_mean_log=sub["tumor_mean_log"]
                .to_numpy()
                .astype(np.float64)
                .reshape((k_clusters, g_genes)),
                all_pct_expressing=sub["all_pct_expressing"]
                .to_numpy()
                .astype(np.float64)
                .reshape((k_clusters, g_genes)),
                tumor_pct_expressing=sub["tumor_pct_expressing"]
                .to_numpy()
                .astype(np.float64)
                .reshape((k_clusters, g_genes)),
            )
        )
    return Success(tuple(bundles))


# =============================================================================
# 2. Substate Marker Extraction & Specificity Algebra
# =============================================================================


def is_informative_hgnc_symbol(gene_sym: str) -> bool:
    """Filter out mitochondrial, ribosomal, lncRNA, and clone symbols for clean marker rankings."""
    g = gene_sym.strip()
    if not g or g.startswith(_EXCLUDED_MARKER_PREFIXES):
        return False
    if g.startswith(("HLA-", "IG", "TR")):
        return not g.startswith(("IGKC", "IGHG", "IGLC", "IGHA", "IGHM"))
    return "-" not in g or g in {"HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-A", "HLA-B", "HLA-C"}


def _select_stratum_matrices(
    bundle: CoultonCentroidBundle,
    stratum: TissueStratum,
) -> tuple[FloatMatrix, FloatMatrix, IntVector]:
    """Select linear centroid, % expressing, and cell counts for the requested tissue stratum."""
    match stratum:
        case TissueStratum.TUMOR_ONLY:
            return bundle.tumor_mean_linear, bundle.tumor_pct_expressing, bundle.tumor_cell_counts
        case TissueStratum.ALL_TISSUES:
            return bundle.all_mean_linear, bundle.all_pct_expressing, bundle.all_cell_counts
        case _ as unreachable:
            assert_never(unreachable)


def compute_centered_log2_specificity(
    density_mat: FloatMatrix,
    reference_row_indices: Sequence[int],
    eps: float = 1e-6,
) -> FloatMatrix:
    """Compute gene-centered log2 relative density specificity relative to `reference_row_indices`.

    Each row of `density_mat` is first normalized to sum to 1.0 over the evaluated gene set,
    transformed via `log2(p + eps)`, and centered by subtracting the mean across `reference_row_indices`.
    """
    row_sums = np.maximum(density_mat.sum(axis=1, keepdims=True), 1e-12)
    normed = density_mat / row_sums
    log_mat = np.log2(normed + eps)
    ref_mean = log_mat[list(reference_row_indices), :].mean(axis=0, keepdims=True)
    return np.asarray(log_mat - ref_mean, dtype=np.float64)


def row_cosine_similarity(mat_a: FloatMatrix, mat_b: FloatMatrix) -> FloatMatrix:
    """Compute exact pairwise row cosine similarity in [-1.0, 1.0] between (M, G) and (K, G)."""
    norm_a = np.maximum(np.linalg.norm(mat_a, axis=1, keepdims=True), 1e-12)
    norm_b = np.maximum(np.linalg.norm(mat_b, axis=1, keepdims=True), 1e-12)
    cos = (mat_a / norm_a) @ (mat_b / norm_b).T
    return np.asarray(np.clip(cos, -1.0, 1.0), dtype=np.float64)


def row_pearson_correlation(mat_a: FloatMatrix, mat_b: FloatMatrix) -> FloatMatrix:
    """Compute exact pairwise row Pearson correlation in [-1.0, 1.0]."""
    cent_a = mat_a - mat_a.mean(axis=1, keepdims=True)
    cent_b = mat_b - mat_b.mean(axis=1, keepdims=True)
    return row_cosine_similarity(cent_a, cent_b)


def row_spearman_correlation(mat_a: FloatMatrix, mat_b: FloatMatrix) -> FloatMatrix:
    """Compute exact pairwise row Spearman rank correlation in [-1.0, 1.0]."""
    ranks_a = np.vstack([rankdata(row) for row in mat_a])
    ranks_b = np.vstack([rankdata(row) for row in mat_b])
    return row_pearson_correlation(ranks_a, ranks_b)


def extract_coulton_top_markers_df(
    bundle: CoultonCentroidBundle,
    stratum: TissueStratum = TissueStratum.TUMOR_ONLY,
    top_n: int = 25,
) -> pl.DataFrame:
    """Extract top positively enriched protein-coding markers per biological Coulton cluster."""
    lin_mat, pct_mat, cell_counts = _select_stratum_matrices(bundle, stratum)
    bio_indices = [
        idx
        for idx, lbl in enumerate(bundle.cluster_labels)
        if COULTON_LITERATURE_CATALOG.get(lbl, COULTON_LITERATURE_CATALOG["23_NA"]).role
        != CoultonClusterRole.UNANNOTATED_ARTIFACT
        and int(cell_counts[idx]) >= 10
    ]
    tam_bg_indices = [
        idx
        for idx in bio_indices
        if COULTON_LITERATURE_CATALOG[bundle.cluster_labels[idx]].role
        in {CoultonClusterRole.TRUE_TAM_SUBSTATE, CoultonClusterRole.CLASSICAL_MONOCYTE}
    ]

    row_sums = np.maximum(lin_mat.sum(axis=1, keepdims=True), 1e-12)
    dens_mat = lin_mat / row_sums
    eps = 1e-6

    informative_mask = np.array(
        [is_informative_hgnc_symbol(g) for g in bundle.gene_symbols],
        dtype=bool,
    )

    rows: list[dict[str, object]] = []
    for k_idx in bio_indices:
        c_label = bundle.cluster_labels[k_idx]
        lit = COULTON_LITERATURE_CATALOG[c_label]
        other_bg = [idx for idx in tam_bg_indices if idx != k_idx]
        bg_vec = dens_mat[other_bg, :].mean(axis=0)
        lfc_vec = np.log2((dens_mat[k_idx, :] + eps) / (bg_vec + eps))

        eligible_mask = informative_mask & (pct_mat[k_idx, :] >= 3.0)
        if lit.role != CoultonClusterRole.T_MYELOID_DOUBLET_CONTROL:
            non_lymph_mask = np.array(
                [g not in _LYMPHOID_SPILLOVER_GENES for g in bundle.gene_symbols],
                dtype=bool,
            )
            eligible_mask = eligible_mask & non_lymph_mask

        eligible_indices = np.flatnonzero(eligible_mask)
        ranked_indices = eligible_indices[np.argsort(-lfc_vec[eligible_indices])][:top_n]

        for rank_pos, g_idx in enumerate(ranked_indices, start=1):
            rows.append(
                {
                    "assay": bundle.assay.value,
                    "tissue_stratum": stratum.value,
                    "coulton_cluster": c_label,
                    "cluster_role": lit.role.value,
                    "cell_count": int(cell_counts[k_idx]),
                    "marker_rank": rank_pos,
                    "gene_symbol": bundle.gene_symbols[int(g_idx)],
                    "specificity_log2fc": float(lfc_vec[int(g_idx)]),
                    "pct_expressing": float(pct_mat[k_idx, int(g_idx)]),
                    "mean_linear_expr": float(lin_mat[k_idx, int(g_idx)]),
                }
            )
    return pl.DataFrame(rows)


# =============================================================================
# 3. Dual-Regime Meta42 <-> Coulton_TAM_Atlas Similarity & Nearest Substates
# =============================================================================


def compute_stratum_similarity_df(
    phi_meta_df: pl.DataFrame,
    coulton_bundle: CoultonCentroidBundle,
    stratum: TissueStratum,
    top_k_markers: int = 25,
    top_degs_df: Maybe[pl.DataFrame] = Nothing,
) -> Result[pl.DataFrame, str]:
    """Compute all pairwise similarities between 42 Meta42 states and 23 Coulton biological clusters.

    When `top_degs_df` (`consensus_metastate_all42_top_degs.parquet`, 168 canonical HGNC state-gene
    pairs from Stage 3c) is provided via `Some(...)`, computes Stage 3c's exact `alpha = 0.40`
    Marker-Anchored Hybrid similarity (`0.60` whole-transcriptome `Phi_meta` specificity cosine
    and HVG Spearman `rho` + `0.40 * tanh(canonical_deg_z / 2.0)`), ensuring full alignment with
    published `Consensus-Meta42` state definitions (`M14: CXCL9+STAT1+ M1-TAM`, `M32: MRC1+CD163+ M2-TAM`,
    `M42: SULT1B1+RGS2+ TAM`) while preserving standalone empirical fallback when `Nothing`.
    """
    meta_states = tuple(phi_meta_df["cluster"].cast(pl.Utf8).to_list())
    meta_gene_cols = tuple(c for c in phi_meta_df.columns if c != "cluster")
    meta_gene_map = {g: idx for idx, g in enumerate(meta_gene_cols)}
    coulton_gene_map = {g: idx for idx, g in enumerate(coulton_bundle.gene_symbols)}

    shared_genes = tuple(sorted(set(meta_gene_cols).intersection(set(coulton_bundle.gene_symbols))))
    if len(shared_genes) < 50:
        return Failure(
            f"Insufficient shared genes ({len(shared_genes)}) between Meta42 and Coulton_TAM_Atlas"
        )
    shared_gene_to_idx = {g: idx for idx, g in enumerate(shared_genes)}

    m_col_idx = [meta_gene_map[g] for g in shared_genes]
    c_col_idx = [coulton_gene_map[g] for g in shared_genes]

    raw_meta_mat = phi_meta_df.select(meta_gene_cols).to_numpy().astype(np.float64)[:, m_col_idx]
    lin_coulton_full, pct_coulton_full, cell_counts_full = _select_stratum_matrices(
        coulton_bundle, stratum
    )

    bio_coulton_indices = [
        idx
        for idx, lbl in enumerate(coulton_bundle.cluster_labels)
        if COULTON_LITERATURE_CATALOG.get(lbl, COULTON_LITERATURE_CATALOG["23_NA"]).role
        != CoultonClusterRole.UNANNOTATED_ARTIFACT
        and int(cell_counts_full[idx]) >= 5
    ]
    bio_coulton_labels = tuple(coulton_bundle.cluster_labels[idx] for idx in bio_coulton_indices)
    raw_coulton_mat = lin_coulton_full[np.ix_(bio_coulton_indices, c_col_idx)]
    pct_coulton_sub = pct_coulton_full[np.ix_(bio_coulton_indices, c_col_idx)]
    cell_counts_sub = cell_counts_full[bio_coulton_indices]

    # Normalize densities on identical shared-gene support
    meta_dens = raw_meta_mat / np.maximum(raw_meta_mat.sum(axis=1, keepdims=True), 1e-12)
    coulton_dens = raw_coulton_mat / np.maximum(raw_coulton_mat.sum(axis=1, keepdims=True), 1e-12)

    # Reference index sets for within-compartment and pan-lineage centering
    mac_meta_idx = [
        i for i, st in enumerate(meta_states) if st in MACROPHAGE_CONTINUUM_META42_IDS
    ]
    mye_meta_idx = [i for i, st in enumerate(meta_states) if st in ALL_MYELOID_META42_IDS]
    all_meta_idx = list(range(len(meta_states)))

    if not mac_meta_idx:
        mac_meta_idx = all_meta_idx
    if not mye_meta_idx:
        mye_meta_idx = all_meta_idx

    tam_mono_coulton_idx = [
        i
        for i, lbl in enumerate(bio_coulton_labels)
        if COULTON_LITERATURE_CATALOG[lbl].role
        in {CoultonClusterRole.TRUE_TAM_SUBSTATE, CoultonClusterRole.CLASSICAL_MONOCYTE}
    ]
    if not tam_mono_coulton_idx:
        tam_mono_coulton_idx = list(range(len(bio_coulton_labels)))

    # 1. Centered log2 relative-density specificities
    spec_meta_mac = compute_centered_log2_specificity(meta_dens, mac_meta_idx)
    spec_meta_mye = compute_centered_log2_specificity(meta_dens, mye_meta_idx)
    spec_meta_all = compute_centered_log2_specificity(meta_dens, all_meta_idx)
    spec_coulton_tam = compute_centered_log2_specificity(coulton_dens, tam_mono_coulton_idx)

    cos_within_mac = row_cosine_similarity(spec_meta_mac, spec_coulton_tam)
    cos_within_mye = row_cosine_similarity(spec_meta_mye, spec_coulton_tam)
    cos_pan_lineage = row_cosine_similarity(spec_meta_all, spec_coulton_tam)

    # 2. Uncentered pan-lineage log2 density Pearson correlation
    log_meta = np.log2(meta_dens + 1e-6)
    log_coulton = np.log2(coulton_dens + 1e-6)
    pearson_pan_log = row_pearson_correlation(log_meta, log_coulton)

    # 3. Spearman rank correlation across substate-variable informative genes (pan-centered to avoid subset distortion)
    informative_flags = np.array([is_informative_hgnc_symbol(g) for g in shared_genes], dtype=bool)
    tam_var = spec_coulton_tam[tam_mono_coulton_idx, :].std(axis=0)
    mye_var = spec_meta_mye[mye_meta_idx, :].std(axis=0)
    combined_var = np.where(informative_flags, tam_var + mye_var, -1.0)
    n_hvg = min(750, int(np.sum(informative_flags)))
    hvg_idx = np.argsort(-combined_var)[: max(n_hvg, 10)]
    spearman_within_mac = row_spearman_correlation(
        spec_meta_all[:, hvg_idx],
        spec_coulton_tam[:, hvg_idx],
    )

    # 4. Standardized gene-level specificity & Stage 3c Canonical 4-DEG Signature Enrichment
    z_meta_cols = (spec_meta_all - spec_meta_all.mean(axis=0, keepdims=True)) / np.maximum(
        spec_meta_all.std(axis=0, keepdims=True), 1e-6
    )
    z_coulton_cols = (
        spec_coulton_tam - spec_coulton_tam[tam_mono_coulton_idx, :].mean(axis=0, keepdims=True)
    ) / np.maximum(spec_coulton_tam[tam_mono_coulton_idx, :].std(axis=0, keepdims=True), 1e-6)

    canonical_deg_z_mat = np.zeros((len(meta_states), len(bio_coulton_labels)), dtype=np.float64)
    canonical_gene_sets: list[set[int]] = [set() for _ in meta_states]
    has_canonical_degs = False

    match top_degs_df:
        case Some(degs_table) if (
            not degs_table.is_empty()
            and {"cell_state", "gene", "specificity_log2fc"}.issubset(set(degs_table.columns))
        ):
            has_canonical_degs = True
            for m_i, m_st in enumerate(meta_states):
                st_rows = degs_table.filter(pl.col("cell_state") == m_st).sort("specificity_log2fc", descending=True)
                matched_g_idx: list[int] = []
                matched_w: list[float] = []
                for deg_row in st_rows.iter_rows(named=True):
                    g_sym = str(deg_row["gene"]).strip()
                    if g_sym in shared_gene_to_idx:
                        matched_g_idx.append(shared_gene_to_idx[g_sym])
                        matched_w.append(max(float(deg_row["specificity_log2fc"]), 0.25))
                if matched_g_idx:
                    canonical_gene_sets[m_i] = set(matched_g_idx)
                    w_arr = np.asarray(matched_w, dtype=np.float64)
                    canonical_deg_z_mat[m_i, :] = (
                        z_coulton_cols[:, matched_g_idx] @ w_arr
                    ) / max(float(w_arr.sum()), 1e-6)
        case _:
            has_canonical_degs = False

    inf_indices = np.flatnonzero(informative_flags)
    k_eff = min(top_k_markers, len(inf_indices))
    k_overlap = min(60, len(inf_indices))

    meta_top_k: list[npt.NDArray[np.intp]] = []
    meta_top_60_sets: list[set[int]] = []
    for m_i in range(len(meta_states)):
        ord_m = inf_indices[np.argsort(-spec_meta_all[m_i, inf_indices])]
        canon_first = sorted(canonical_gene_sets[m_i], key=lambda gi: -spec_meta_all[m_i, gi])
        merged_top = canon_first + [int(x) for x in ord_m if int(x) not in canonical_gene_sets[m_i]]
        meta_top_k.append(np.asarray(merged_top[:k_eff], dtype=np.intp))
        meta_top_60_sets.append(set(merged_top[:k_overlap]).union(canonical_gene_sets[m_i]))

    coulton_top_k: list[npt.NDArray[np.intp]] = []
    coulton_top_60_sets: list[set[int]] = []
    for c_i, c_lbl in enumerate(bio_coulton_labels):
        role = COULTON_LITERATURE_CATALOG[c_lbl].role
        cand_mask = informative_flags & (pct_coulton_sub[c_i, :] >= 2.0)
        if role != CoultonClusterRole.T_MYELOID_DOUBLET_CONTROL:
            non_lym = np.array([g not in _LYMPHOID_SPILLOVER_GENES for g in shared_genes], dtype=bool)
            cand_mask = cand_mask & non_lym
        cand_idx = np.flatnonzero(cand_mask)
        if len(cand_idx) < k_eff:
            cand_idx = inf_indices
        ord_c = cand_idx[np.argsort(-spec_coulton_tam[c_i, cand_idx])]
        coulton_top_k.append(ord_c[:k_eff])
        coulton_top_60_sets.append(set(int(x) for x in ord_c[:k_overlap]))

    bidir_z_mat = np.zeros((len(meta_states), len(bio_coulton_labels)), dtype=np.float64)
    shared_genes_str: list[list[str]] = [
        ["" for _ in bio_coulton_labels] for _ in meta_states
    ]
    shared_gene_counts = np.zeros((len(meta_states), len(bio_coulton_labels)), dtype=np.int32)

    for m_i in range(len(meta_states)):
        m_idx_vec = meta_top_k[m_i]
        c_set_m = canonical_gene_sets[m_i]
        for c_i in range(len(bio_coulton_labels)):
            c_idx_vec = coulton_top_k[c_i]
            z_m_in_c = float(z_coulton_cols[c_i, m_idx_vec].mean()) if len(m_idx_vec) > 0 else 0.0
            z_c_in_m = float(z_meta_cols[m_i, c_idx_vec].mean()) if len(c_idx_vec) > 0 else 0.0
            bidir_z_mat[m_i, c_i] = 0.5 * (z_m_in_c + z_c_in_m)

            ov_idx = sorted(
                meta_top_60_sets[m_i].intersection(coulton_top_60_sets[c_i]),
                key=lambda gi: -(
                    (4.0 if gi in c_set_m else 0.0)
                    + spec_meta_all[m_i, gi]
                    + spec_coulton_tam[c_i, gi]
                ),
            )
            shared_gene_counts[m_i, c_i] = len(ov_idx)
            shared_genes_str[m_i][c_i] = ", ".join(shared_genes[gi] for gi in ov_idx[:8])

    composite_mat = (
        (
            0.35 * cos_pan_lineage
            + 0.25 * spearman_within_mac
            + 0.40 * np.tanh(canonical_deg_z_mat / 2.0)
        )
        if has_canonical_degs
        else (
            0.50 * cos_within_mac
            + 0.30 * spearman_within_mac
            + 0.20 * np.tanh(bidir_z_mat / 3.0)
        )
    )
    hybrid_pan_affinity_mat = (
        0.60 * cos_pan_lineage + 0.40 * np.tanh(canonical_deg_z_mat / 2.5)
        if has_canonical_degs
        else cos_pan_lineage
    )

    true_tam_local_indices = [
        i
        for i, lbl in enumerate(bio_coulton_labels)
        if COULTON_LITERATURE_CATALOG[lbl].role == CoultonClusterRole.TRUE_TAM_SUBSTATE
    ]

    rows: list[dict[str, object]] = []
    for m_i, m_state in enumerate(meta_states):
        tam_scores = composite_mat[m_i, true_tam_local_indices]
        tam_rank_order = np.argsort(-tam_scores)
        tam_ranks: dict[int, int] = {
            true_tam_local_indices[int(pos)]: rank_num
            for rank_num, pos in enumerate(tam_rank_order, start=1)
        }

        all_rank_order = np.argsort(-composite_mat[m_i, :])
        all_ranks: dict[int, int] = {
            int(pos): rank_num for rank_num, pos in enumerate(all_rank_order, start=1)
        }

        m_label = META42_BIO_LABELS.get(m_state, m_state)
        m_comp = classify_meta42_compartment(m_state).value
        m_tier = classify_meta42_macrophage_tier(m_state).value

        for c_i, c_label in enumerate(bio_coulton_labels):
            lit = COULTON_LITERATURE_CATALOG[c_label]
            rows.append(
                {
                    "assay": coulton_bundle.assay.value,
                    "tissue_stratum": stratum.value,
                    "n_shared_genes": len(shared_genes),
                    "meta_state": m_state,
                    "meta_state_label": m_label,
                    "lineage_compartment": m_comp,
                    "macrophage_tier": m_tier,
                    "coulton_cluster": c_label,
                    "coulton_cluster_id": lit.cluster_id,
                    "coulton_role": lit.role.value,
                    "coulton_cell_count": int(cell_counts_sub[c_i]),
                    "within_mac_spec_cosine": float(cos_within_mac[m_i, c_i]),
                    "within_myeloid_spec_cosine": float(cos_within_mye[m_i, c_i]),
                    "within_mac_spec_spearman": float(spearman_within_mac[m_i, c_i]),
                    "bidirectional_marker_z": float(bidir_z_mat[m_i, c_i]),
                    "canonical_deg_z": float(canonical_deg_z_mat[m_i, c_i]),
                    "pan_lineage_log_pearson": float(pearson_pan_log[m_i, c_i]),
                    "pan_lineage_spec_cosine": float(cos_pan_lineage[m_i, c_i]),
                    "hybrid_pan_affinity": float(hybrid_pan_affinity_mat[m_i, c_i]),
                    "composite_substate_score": float(composite_mat[m_i, c_i]),
                    "rank_among_true_tams": tam_ranks.get(c_i, 99),
                    "rank_among_all_coulton": all_ranks[c_i],
                    "n_shared_top60_markers": int(shared_gene_counts[m_i, c_i]),
                    "shared_signature_genes": shared_genes_str[m_i][c_i],
                }
            )

    return Success(pl.DataFrame(rows))


# =============================================================================
# 4. Empirical Bulk ICB Logistic Meta-Analysis for Coulton Substates
# =============================================================================


def fit_penalized_logistic_newton(
    design_x: FloatMatrix,
    response_y: FloatVector,
    l2_penalty: float = 1e-3,
    max_iter: int = 35,
) -> tuple[FloatVector, FloatVector]:
    """Pure L2-stabilized Newton-Raphson logistic regression returning (beta, se).

    First column of `design_x` is assumed to be the unpenalized intercept `1.0`.
    """
    n_obs, n_feat = design_x.shape
    beta = np.zeros(n_feat, dtype=np.float64)
    pen_diag = np.full(n_feat, l2_penalty, dtype=np.float64)
    pen_diag[0] = 1e-8

    for _ in range(max_iter):
        eta = np.clip(design_x @ beta, -20.0, 20.0)
        prob = 1.0 / (1.0 + np.exp(-eta))
        weights = np.maximum(prob * (1.0 - prob), 1e-6)
        grad = design_x.T @ (response_y - prob) - pen_diag * beta
        hess = (design_x.T * weights) @ design_x + np.diag(pen_diag)
        step = np.linalg.solve(hess, grad)
        beta = beta + step
        if float(np.max(np.abs(step))) < 1e-7:
            break

    eta_f = np.clip(design_x @ beta, -20.0, 20.0)
    prob_f = 1.0 / (1.0 + np.exp(-eta_f))
    weights_f = np.maximum(prob_f * (1.0 - prob_f), 1e-6)
    fisher_info = (design_x.T * weights_f) @ design_x + np.diag(pen_diag)
    cov = np.linalg.inv(fisher_info)
    se = np.sqrt(np.maximum(np.diag(cov), 1e-12))
    return np.asarray(beta, dtype=np.float64), np.asarray(se, dtype=np.float64)


def compute_mann_whitney_auroc(scores: FloatVector, response_y: FloatVector) -> float:
    """Exact rank-based Mann-Whitney U / ROC-AUC without external sklearn state."""
    pos_mask = response_y == 1.0
    neg_mask = response_y == 0.0
    n_pos = int(np.sum(pos_mask))
    n_neg = int(np.sum(neg_mask))
    if n_pos == 0 or n_neg == 0:
        return 0.5
    ranks = rankdata(scores)
    u_stat = float(np.sum(ranks[pos_mask])) - (n_pos * (n_pos + 1)) / 2.0
    return float(np.clip(u_stat / float(n_pos * n_neg), 0.0, 1.0))


def benjamini_hochberg_fdr(pvalues: FloatVector) -> FloatVector:
    """Pure Benjamini-Hochberg FDR adjustment in [0.0, 1.0]."""
    m = len(pvalues)
    if m == 0:
        return np.array([], dtype=np.float64)
    order = np.argsort(pvalues)
    ranked_p = pvalues[order]
    ranks = np.arange(1, m + 1, dtype=np.float64)
    raw_q = (ranked_p * float(m)) / ranks
    monotone_q = np.minimum.accumulate(raw_q[::-1])[::-1]
    out = np.empty(m, dtype=np.float64)
    out[order] = np.clip(monotone_q, 0.0, 1.0)
    return out


def inverse_variance_meta_analysis(
    betas: FloatVector,
    ses: FloatVector,
) -> tuple[float, float, float]:
    """Fixed-effects inverse-variance meta-analysis returning (meta_beta, meta_se, two_sided_pval)."""
    valid = np.isfinite(betas) & np.isfinite(ses) & (ses > 1e-8)
    if int(np.sum(valid)) == 0:
        return 0.0, 1.0, 1.0
    w = 1.0 / np.square(ses[valid])
    w_sum = float(np.sum(w))
    meta_beta = float(np.sum(w * betas[valid]) / w_sum)
    meta_se = float(math.sqrt(1.0 / w_sum))
    z_val = meta_beta / max(meta_se, 1e-12)
    pval = float(2.0 * norm.sf(abs(z_val)))
    return meta_beta, meta_se, pval


def _standardize_vector(vec: FloatVector) -> FloatVector:
    std = float(np.std(vec))
    return (
        np.zeros_like(vec, dtype=np.float64)
        if std < 1e-8
        else np.asarray((vec - float(np.mean(vec))) / std, dtype=np.float64)
    )


def evaluate_coulton_signatures_on_bulk_icb(
    coulton_markers_df: pl.DataFrame,
    bulk_cohorts: Sequence[BulkIcbCohortBundle],
    n_sig_genes: int = 15,
) -> pl.DataFrame:
    """Evaluate all 23 biological Coulton substate signatures across the 11 clinical bulk ICB cohorts."""
    bio_clusters = [
        lbl
        for lbl in COULTON_ORDERED_CLUSTERS
        if COULTON_LITERATURE_CATALOG[lbl].role != CoultonClusterRole.UNANNOTATED_ARTIFACT
    ]

    sig_map: dict[str, tuple[str, ...]] = {}
    for c_lbl in bio_clusters:
        sub_m = (
            coulton_markers_df.filter(pl.col("coulton_cluster") == c_lbl)
            .sort("marker_rank")
            .head(n_sig_genes)
        )
        sig_map[c_lbl] = tuple(sub_m["gene_symbol"].cast(pl.Utf8).to_list())

    raw_rows: list[dict[str, object]] = []
    for c_lbl in bio_clusters:
        lit = COULTON_LITERATURE_CATALOG[c_lbl]
        sig_genes = sig_map[c_lbl]

        cohort_raw_betas: list[float] = []
        cohort_raw_ses: list[float] = []
        cohort_adj_betas: list[float] = []
        cohort_adj_ses: list[float] = []
        cohort_aurocs: list[float] = []
        cohort_weights: list[int] = []
        n_total_pts = 0

        for cb in bulk_cohorts:
            present_sig = [cb.gene_to_idx[g] for g in sig_genes if g in cb.gene_to_idx]
            present_pan = [
                cb.gene_to_idx[g] for g in PAN_MACROPHAGE_CORE_GENES if g in cb.gene_to_idx
            ]
            if len(present_sig) < 2 or int(np.sum(cb.response_binary == 1.0)) < 2:
                continue

            s_raw = _standardize_vector(cb.log2_z_mat[:, present_sig].mean(axis=1))
            intercept = np.ones(cb.n_samples, dtype=np.float64)
            x_uni = np.column_stack([intercept, s_raw])
            b_uni, se_uni = fit_penalized_logistic_newton(x_uni, cb.response_binary)
            cohort_raw_betas.append(float(b_uni[1]))
            cohort_raw_ses.append(float(se_uni[1]))

            if len(present_pan) >= 3:
                m_pan = _standardize_vector(cb.log2_z_mat[:, present_pan].mean(axis=1))
                x_adj = np.column_stack([intercept, s_raw, m_pan])
                b_adj, se_adj = fit_penalized_logistic_newton(x_adj, cb.response_binary)
                cohort_adj_betas.append(float(b_adj[1]))
                cohort_adj_ses.append(float(se_adj[1]))

            cohort_aurocs.append(compute_mann_whitney_auroc(s_raw, cb.response_binary))
            cohort_weights.append(cb.n_samples)
            n_total_pts += cb.n_samples

        raw_beta, raw_se, raw_p = inverse_variance_meta_analysis(
            np.asarray(cohort_raw_betas, dtype=np.float64),
            np.asarray(cohort_raw_ses, dtype=np.float64),
        )
        adj_beta, adj_se, adj_p = inverse_variance_meta_analysis(
            np.asarray(cohort_adj_betas, dtype=np.float64),
            np.asarray(cohort_adj_ses, dtype=np.float64),
        )
        weighted_auroc = (
            float(np.average(cohort_aurocs, weights=cohort_weights))
            if cohort_weights
            else 0.5
        )

        match lit.cpi1000_fgsea_qval:
            case Some(q_val):
                cpi_q = float(q_val)
                sign_mult = (
                    1.0
                    if lit.cpi1000_direction == IciOutcomeDirection.RESPONDER_ENRICHED
                    else -1.0
                )
                signed_neg_log10_q = sign_mult * (-math.log10(max(cpi_q, 1e-20)))
            case _:
                cpi_q = 1.0
                signed_neg_log10_q = 0.0

        raw_rows.append(
            {
                "coulton_cluster": c_lbl,
                "coulton_cluster_id": lit.cluster_id,
                "coulton_role": lit.role.value,
                "functional_program": lit.functional_program,
                "canonical_markers": ", ".join(lit.canonical_markers),
                "empirical_signature_genes": ", ".join(sig_genes[:8]),
                "is_ucell_gold_standard": lit.is_ucell_gold_standard,
                "cpi1000_published_direction": lit.cpi1000_direction.value,
                "cpi1000_fgsea_qval": cpi_q,
                "cpi1000_is_significant": cpi_q < 0.10,
                "cpi1000_signed_neg_log10_q": signed_neg_log10_q,
                "n_bulk_icb_patients": n_total_pts,
                "iatlas_raw_meta_beta": raw_beta,
                "iatlas_raw_meta_se": raw_se,
                "iatlas_raw_meta_pval": raw_p,
                "iatlas_mac_adj_meta_beta": adj_beta,
                "iatlas_mac_adj_meta_se": adj_se,
                "iatlas_mac_adj_meta_pval": adj_p,
                "iatlas_weighted_auroc": weighted_auroc,
                "mana_and_spatial_note": lit.mana_and_spatial_note,
            }
        )

    out_df = pl.DataFrame(raw_rows)
    raw_fdr = benjamini_hochberg_fdr(out_df["iatlas_raw_meta_pval"].to_numpy().astype(np.float64))
    adj_fdr = benjamini_hochberg_fdr(
        out_df["iatlas_mac_adj_meta_pval"].to_numpy().astype(np.float64)
    )
    return out_df.with_columns(
        [
            pl.Series("iatlas_raw_meta_fdr", raw_fdr),
            pl.Series("iatlas_mac_adj_meta_fdr", adj_fdr),
        ]
    )


# =============================================================================
# 5. Joint Nearest-Substate & Triangulated ICI Outcome Synthesis
# =============================================================================


def _classify_joint_concordance(
    meta_bulk_beta: float,
    meta_bulk_pval: float,
    coulton_cpi_dir: str,
    coulton_mac_adj_beta: float,
) -> OutcomeConcordanceClass:
    """Classify ICI therapy outcome concordance between a Meta42 state and its matched Coulton TAM."""
    meta_sign = 0 if (abs(meta_bulk_beta) < 0.20 and meta_bulk_pval > 0.15) else (1 if meta_bulk_beta > 0.0 else -1)

    if coulton_cpi_dir == IciOutcomeDirection.RESPONDER_ENRICHED.value:
        coulton_sign = 1
    elif coulton_cpi_dir == IciOutcomeDirection.NON_RESPONDER_ENRICHED.value:
        coulton_sign = -1
    elif abs(coulton_mac_adj_beta) >= 0.08:
        coulton_sign = 1 if coulton_mac_adj_beta > 0.0 else -1
    else:
        coulton_sign = 0

    if meta_sign == 0 or coulton_sign == 0:
        return OutcomeConcordanceClass.UNDETERMINED_NEUTRAL
    if meta_sign > 0 and coulton_sign > 0:
        return OutcomeConcordanceClass.CONCORDANT_PRO_RESPONSE
    if meta_sign < 0 and coulton_sign < 0:
        return OutcomeConcordanceClass.CONCORDANT_PRO_RESISTANCE
    return OutcomeConcordanceClass.POLARIZATION_SPLIT


def build_nearest_substates_and_ici_tables(
    full_sim_df: pl.DataFrame,
    coulton_ici_df: pl.DataFrame,
    meta42_bundle: Meta42ReferenceBundle,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Build ranked Top-3 Coulton TAM matches + ICI outcome summary + sensitivity + indication tables."""
    # Primary stratum: RNA + tumor_only
    primary_sim = full_sim_df.filter(
        (pl.col("assay") == AssaySource.RNA_LOG1P.value)
        & (pl.col("tissue_stratum") == TissueStratum.TUMOR_ONLY.value)
    )

    pan_eff = meta42_bundle.state_effects_df.filter(
        (pl.col("reference_id") == "Consensus-Meta42")
        & (pl.col("bulk_cohort") == "Pan-iAtlas-Meta")
    ).select(
        [
            pl.col("cell_state").alias("meta_state"),
            pl.col("milo_mean_logfc").alias("meta42_milo_logfc"),
            pl.col("milo_se").alias("meta42_milo_se"),
            pl.col("bulk_deconv_beta").alias("meta42_bulk_beta"),
            pl.col("bulk_deconv_se").alias("meta42_bulk_se"),
            pl.col("bulk_deconv_pval").alias("meta42_bulk_pval"),
            pl.col("bulk_multivariate_beta").alias("meta42_multivariate_beta"),
            pl.col("delta_clr").alias("meta42_delta_clr"),
            pl.col("quadrant").alias("meta42_concordance_quadrant"),
        ]
    )

    hub_sub = meta42_bundle.hub_summary_df.select(
        [
            pl.col("cell_state").alias("meta_state"),
            pl.col("enrichment_group").alias("meta42_enrichment_group"),
            pl.col("mean_abs_epistasis_z").alias("meta42_mean_abs_epistasis_z"),
            pl.col("max_synergy_z").alias("meta42_max_synergy_z"),
            pl.col("min_gating_z").alias("meta42_min_gating_z"),
        ]
    )

    best_overall = (
        primary_sim.filter(pl.col("rank_among_all_coulton") == 1)
        .select(
            [
                "meta_state",
                pl.col("coulton_cluster").alias("best_overall_coulton_cluster"),
                pl.col("coulton_role").alias("best_overall_coulton_role"),
                pl.col("composite_substate_score").alias("best_overall_score"),
            ]
        )
    )

    top3_true_tams = (
        primary_sim.filter(
            (pl.col("coulton_role") == CoultonClusterRole.TRUE_TAM_SUBSTATE.value)
            & (pl.col("rank_among_true_tams") <= 3)
        )
        .join(coulton_ici_df, on=["coulton_cluster", "coulton_cluster_id"], how="left")
        .join(pan_eff, on="meta_state", how="left")
        .join(hub_sub, on="meta_state", how="left")
        .join(best_overall, on="meta_state", how="left")
        .sort(["meta_state", "rank_among_true_tams"])
    )

    concordance_labels = [
        _classify_joint_concordance(
            meta_bulk_beta=float(r["meta42_bulk_beta"]),
            meta_bulk_pval=float(r["meta42_bulk_pval"]),
            coulton_cpi_dir=str(r["cpi1000_published_direction"]),
            coulton_mac_adj_beta=float(r["iatlas_mac_adj_meta_beta"]),
        ).value
        for r in top3_true_tams.iter_rows(named=True)
    ]
    meta_odds_ratios = np.exp(
        np.clip(top3_true_tams["meta42_bulk_beta"].to_numpy().astype(np.float64), -10.0, 10.0)
    )
    nearest_df = top3_true_tams.with_columns(
        [
            pl.Series("meta42_bulk_odds_ratio", meta_odds_ratios),
            pl.Series("ici_outcome_concordance", concordance_labels),
        ]
    )

    # Sensitivity table across (RNA vs SCT) x (tumor_only vs all_tissues) for Myeloid Meta42 states
    sens_top1 = (
        full_sim_df.filter(
            pl.col("meta_state").is_in(list(ALL_MYELOID_META42_IDS))
            & (pl.col("coulton_role") == CoultonClusterRole.TRUE_TAM_SUBSTATE.value)
            & (pl.col("rank_among_true_tams") == 1)
        )
        .select(
            [
                "meta_state",
                "meta_state_label",
                "macrophage_tier",
                "assay",
                "tissue_stratum",
                pl.col("coulton_cluster").alias("top1_tam_substate"),
                "pan_lineage_spec_cosine",
                "within_mac_spec_cosine",
                "within_mac_spec_spearman",
                "bidirectional_marker_z",
                "canonical_deg_z",
                "composite_substate_score",
                "shared_signature_genes",
            ]
        )
        .sort(["macrophage_tier", "meta_state", "assay", "tissue_stratum"])
    )

    # Per-indication bulk ICB effects for all Myeloid Meta42 states
    indication_cohorts = (
        "Pan-iAtlas-Meta",
        "PanCancer-Combined",
        "Melanoma-Combined",
        "RenalCell-Combined",
        "Rosenberg-iAtlas",
        "NSCLC-Combined",
        "Padron-iAtlas",
        "Anders-iAtlas",
    )
    indication_df = (
        meta42_bundle.state_effects_df.filter(
            (pl.col("reference_id") == "Consensus-Meta42")
            & pl.col("cell_state").is_in(list(ALL_MYELOID_META42_IDS))
            & pl.col("bulk_cohort").is_in(list(indication_cohorts))
        )
        .select(
            [
                pl.col("cell_state").alias("meta_state"),
                pl.col("bulk_cohort"),
                pl.col("bulk_short_name"),
                pl.col("bulk_deconv_beta"),
                pl.col("bulk_deconv_se"),
                pl.col("bulk_deconv_pval"),
                pl.col("milo_mean_logfc"),
                pl.col("quadrant"),
            ]
        )
        .with_columns(
            pl.col("meta_state")
            .replace_strict(dict(META42_BIO_LABELS), default="Unknown")
            .alias("meta_state_label")
        )
        .sort(["meta_state", "bulk_cohort"])
    )

    return nearest_df, sens_top1, indication_df


# =============================================================================
# 6. Pure Markdown Report Generator
# =============================================================================


def render_comparison_markdown_report(
    nearest_df: pl.DataFrame,
    sensitivity_df: pl.DataFrame,
    coulton_ici_df: pl.DataFrame,
    full_sim_df: pl.DataFrame,
    indication_df: pl.DataFrame,
) -> str:
    """Render a publication-grade Markdown synthesis report of Meta42 vs Coulton_TAM_Atlas."""
    n_shared_rna = int(
        full_sim_df.filter(pl.col("assay") == AssaySource.RNA_LOG1P.value)["n_shared_genes"][0]
    )
    lines: list[str] = [
        "# Consensus-Meta42 vs. Coulton et al. (2024) Pan-Cancer TAM Atlas: Substate Mapping & ICI Therapy Outcomes",
        "",
        "## 1. Study Scope & Alignment Summary",
        f"- **Shared Gene Space**: `{n_shared_rna:,} / 4,529` (`{100.0 * n_shared_rna / 4529.0:.2f}%`) HGNC symbols shared between `Consensus-Meta42` (`unified_reference_phi_meta.parquet` + `consensus_metastate_all42_top_degs.parquet`) and `Coulton_TAM_Atlas` (`363,315` total cells; `279,104` tumor-infiltrating cells across `32` studies and `22` cancer types).",
        "- **Primary Similarity Metric**: Stage `3c`-aligned Marker-Anchored Hybrid (`alpha = 0.40`) combining `0.60` whole-transcriptome `Phi_meta` specificity (`0.35` centered $\\log_2$ relative-density cosine + `0.25` HVG Spearman $\\rho$) with `0.40` canonical top-4 HGNC marker enrichment $\\tanh(z_{\\text{DEG}}/2)$ (`168` state-gene pairs), evaluated on tumor-infiltrating TAMs (`RNA@data` linear `expm1` centroids) and cross-checked across `All-Tissues` (`n = 363,315`) and `SCT@data`.",
        "- **Triangulated Clinical ICI Outcome Sources**:",
        "  1. **`Meta42` Pan-iAtlas Bulk Deconvolution (`N = 1,015`) & 9-Cohort scRNA-seq Milo DA (`N = 611` donors)**.",
        "  2. **Coulton et al. (*Nature Communications* 2024, `PMC11226649`) `CPI1000+` (`N = 1,446` ICI patients)** DESeq2 + `fgsea` enrichment FDR $q$-values and lung MANA neoantigen T-cell / NanoString CosMx spatial validation.",
        "  3. **Direct Empirical Signature Meta-Analysis across 11 `tme_datasets` Bulk ICB Cohorts (`N = 1,058` patients)** reporting both univariate logistic $\\beta_{\\text{ICB}}$ and pan-macrophage-adjusted polarization $\\beta_{\\text{ICB}}^{\\text{adj}}$.",
        "",
        "---",
        "",
        "## 2. Closest Coulton Macrophage Substates for `Meta42` Macrophage & Myeloid States",
        "",
        "### 2A. Tier 1 — Core Explicitly Annotated TAM States (`M14`, `M32`, `M42`)",
        "",
        "| Meta42 State | Meta42 Bulk $\\beta$ (`OR`, $p$) | Meta42 Milo $\\log_2\\text{FC}$ | Rank | Closest Coulton TAM Substate | Canonical DEG $z$ | Spec Cosine | Composite | Shared Top Markers | Coulton `CPI1000+` (`N=1,446`) | Coulton iAtlas $\\beta_{\\text{ICB}}^{\\text{adj}}$ ($p$) | Outcome Concordance |",
        "| :--- | :--- | :--- | :---: | :--- | :---: | :---: | :---: | :--- | :--- | :--- | :--- |",
    ]

    tier1_rows = nearest_df.filter(pl.col("meta_state").is_in(list(CORE_TAM_META42_IDS))).sort(
        ["meta_state", "rank_among_true_tams"]
    )
    for r in tier1_rows.iter_rows(named=True):
        cpi_str = (
            f"{r['cpi1000_published_direction']} (q={float(r['cpi1000_fgsea_qval']):.2e})"
            if bool(r["cpi1000_is_significant"])
            else "NS (q >= 0.10)"
        )
        lines.append(
            f"| **{r['meta_state_label']}** | `{float(r['meta42_bulk_beta']):+.3f}` (`{float(r['meta42_bulk_odds_ratio']):.2f}x`, `{float(r['meta42_bulk_pval']):.2e}`) | "
            f"`{float(r['meta42_milo_logfc']):+.2f}` | `#{int(r['rank_among_true_tams'])}` | **{r['coulton_cluster']}** | "
            f"`{float(r['canonical_deg_z']):+.2f}` | `{float(r['pan_lineage_spec_cosine']):+.3f}` | `{float(r['composite_substate_score']):+.3f}` | "
            f"`{r['shared_signature_genes']}` | {cpi_str} | `{float(r['iatlas_mac_adj_meta_beta']):+.3f}` (`{float(r['iatlas_mac_adj_meta_pval']):.2e}`) | {r['ici_outcome_concordance']} |"
        )

    lines.extend(
        [
            "",
            "### 2B. Tier 2 & Tier 3 — Monocyte / MDSC & Extended Myeloid–TAM Continuum (`M10`, `M19`, `M28`, `M35`, `M09`, `M31`, `M38`, `M39`, `M41`)",
            "",
            "| Meta42 State | Meta42 Bulk $\\beta$ ($p$) | Meta42 Milo $\\log_2\\text{FC}$ | Best Overall (incl. `19_ClassMono`) | Top-1 True TAM Substate | Top-2 True TAM Substate | Top-3 True TAM Substate | Top-1 Shared Markers | Top-1 Coulton `CPI1000+` / iAtlas $\\beta_{\\text{ICB}}^{\\text{adj}}$ |",
            "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |",
        ]
    )

    tier23_ids = ("M10", "M19", "M28", "M35", "M09", "M41", "M31", "M38", "M39")
    for m_id in tier23_ids:
        sub_m = nearest_df.filter(pl.col("meta_state") == m_id).sort("rank_among_true_tams")
        if sub_m.height < 3:
            continue
        r1 = sub_m.row(0, named=True)
        r2 = sub_m.row(1, named=True)
        r3 = sub_m.row(2, named=True)
        cpi1 = (
            f"{r1['cpi1000_published_direction']} (q={float(r1['cpi1000_fgsea_qval']):.2e})"
            if bool(r1["cpi1000_is_significant"])
            else "CPI1000+ NS"
        )
        lines.append(
            f"| **{r1['meta_state_label']}** | `{float(r1['meta42_bulk_beta']):+.3f}` (`{float(r1['meta42_bulk_pval']):.2e}`) | "
            f"`{float(r1['meta42_milo_logfc']):+.2f}` | `{r1['best_overall_coulton_cluster']}` (`{float(r1['best_overall_score']):+.3f}`) | "
            f"**{r1['coulton_cluster']}** (`{float(r1['composite_substate_score']):+.3f}`) | "
            f"{r2['coulton_cluster']} (`{float(r2['composite_substate_score']):+.3f}`) | "
            f"{r3['coulton_cluster']} (`{float(r3['composite_substate_score']):+.3f}`) | "
            f"`{r1['shared_signature_genes']}` | {cpi1}; adj $\\beta$=`{float(r1['iatlas_mac_adj_meta_beta']):+.3f}` |"
        )

    lines.extend(
        [
            "",
            "---",
            "",
            "## 3. Cross-Assay (`RNA` vs. `SCT`) & Tissue (`Tumor` vs. `All-Tissues`) Sensitivity Check",
            "",
            "| Meta42 State | Assay | Tissue Stratum | Top-1 True TAM Substate | Canonical DEG $z$ | Spec Cosine | Spearman $\\rho$ | Composite |",
            "| :--- | :---: | :--- | :--- | :---: | :---: | :---: | :---: |",
        ]
    )
    for r in sensitivity_df.filter(
        pl.col("meta_state").is_in(list(MACROPHAGE_CONTINUUM_META42_IDS) + ["M09"])
    ).iter_rows(named=True):
        lines.append(
            f"| {r['meta_state_label']} | `{r['assay']}` | `{r['tissue_stratum']}` | **{r['top1_tam_substate']}** | "
            f"`{float(r['canonical_deg_z']):+.2f}` | `{float(r['pan_lineage_spec_cosine']):+.3f}` | `{float(r['within_mac_spec_spearman']):+.3f}` | "
            f"`{float(r['composite_substate_score']):+.3f}` |"
        )

    lines.extend(
        [
            "",
            "---",
            "",
            "## 4. Complete Coulton et al. (2024) Substate ICI Outcome Catalog (`CPI1000+` & 11-Cohort `Pan-iAtlas`)",
            "",
            "| Coulton Substate | Role | Gold-Std Bulk Sig | Coulton `CPI1000+` (`N=1,446`) | `CPI1000+` FDR $q$ | Empirical iAtlas Raw $\\beta$ ($p$) | Pan-Mac-Adjusted $\\beta_{\\text{ICB}}^{\\text{adj}}$ ($p$) | Mean AUROC | Canonical Markers |",
            "| :--- | :--- | :---: | :--- | :---: | :--- | :--- | :---: | :--- |",
        ]
    )
    for r in coulton_ici_df.sort("coulton_cluster_id").iter_rows(named=True):
        q_txt = f"{float(r['cpi1000_fgsea_qval']):.2e}" if bool(r["cpi1000_is_significant"]) else "NS"
        gs_txt = "Yes" if bool(r["is_ucell_gold_standard"]) else "No"
        lines.append(
            f"| **{r['coulton_cluster']}** | {r['coulton_role']} | {gs_txt} | {r['cpi1000_published_direction']} | `{q_txt}` | "
            f"`{float(r['iatlas_raw_meta_beta']):+.3f}` (`{float(r['iatlas_raw_meta_pval']):.2e}`) | "
            f"`{float(r['iatlas_mac_adj_meta_beta']):+.3f}` (`{float(r['iatlas_mac_adj_meta_pval']):.2e}`) | "
            f"`{float(r['iatlas_weighted_auroc']):.3f}` | `{r['canonical_markers']}` |"
        )

    # Cross-lineage sanity checks (8_IFNGMac, 18_ECMMac, 14_ProliMac, 20_TDoub)
    pan_primary = full_sim_df.filter(
        (pl.col("assay") == AssaySource.RNA_LOG1P.value)
        & (pl.col("tissue_stratum") == TissueStratum.TUMOR_ONLY.value)
    )
    lines.extend(
        [
            "",
            "---",
            "",
            "## 5. Cross-Lineage Mechanistic Bridges Across All 42 `Meta42` States",
            "",
        ]
    )
    for check_cluster, desc in (
        ("8_IFNGMac", "IFN-gamma CXCL9+/CXCL10+ T-cell-recruiting TAM"),
        ("4_ICIMac2", "Lipid-associated TREM2+/APOE+/APOC1+ ICI-responder TAM"),
        ("18_ECMMac", "Collagen-remodeling fibroblast-adjacent ICI-resistance TAM"),
        ("14_ProliMac", "Cell-cycle MKI67+/TOP2A+ proliferating TAM"),
        ("20_TDoub", "T/Myeloid doublet positive control"),
    ):
        top_pan = (
            pan_primary.filter(pl.col("coulton_cluster") == check_cluster)
            .sort("hybrid_pan_affinity", descending=True)
            .head(3)
        )
        top_str = "; ".join(
            f"**{r['meta_state_label']}** (hybrid affinity = `{float(r['hybrid_pan_affinity']):+.3f}`, shared: `{r['shared_signature_genes']}`)"
            for r in top_pan.iter_rows(named=True)
        )
        lines.append(f"- **`{check_cluster}`** ({desc}) -> Top pan-lineage `Meta42` affinities: {top_str}.")

    lines.extend(
        [
            "",
            "---",
            "",
            "## 6. Indication-Stratified Bulk Deconvolution Effects for `Meta42` Macrophage Continuum",
            "",
            "| Meta42 State | Bulk ICB Stratum | Deconv $\\beta$ | SE | $p$-value | Quadrant |",
            "| :--- | :--- | :---: | :---: | :---: | :--- |",
        ]
    )
    for r in indication_df.filter(
        pl.col("meta_state").is_in(list(CORE_TAM_META42_IDS) + ["M09", "M19"])
    ).iter_rows(named=True):
        lines.append(
            f"| {r['meta_state_label']} | {r['bulk_short_name']} | `{float(r['bulk_deconv_beta']):+.3f}` | "
            f"`{float(r['bulk_deconv_se']):.3f}` | `{float(r['bulk_deconv_pval']):.2e}` | {r['quadrant']} |"
        )

    return "\n".join(lines) + "\n"


def run_comparison_core_pure(
    meta42_bundle: Meta42ReferenceBundle,
    coulton_bundles: Sequence[CoultonCentroidBundle],
    bulk_cohorts: Sequence[BulkIcbCohortBundle],
    config: ComparisonConfig,
) -> Result[ComparisonOutputs, str]:
    """Pure end-to-end functional pipeline orchestrating similarity, ICI meta-analysis, and reporting."""
    if not coulton_bundles:
        return Failure("No Coulton centroid bundles provided")

    rna_candidates = [b for b in coulton_bundles if b.assay == AssaySource.RNA_LOG1P]
    if not rna_candidates:
        return Failure("RNA CoultonCentroidBundle is required for primary substate matching")
    rna_bundle = rna_candidates[0]

    sim_frames: list[pl.DataFrame] = []
    for cb in coulton_bundles:
        for stratum in (TissueStratum.TUMOR_ONLY, TissueStratum.ALL_TISSUES):
            match compute_stratum_similarity_df(
                phi_meta_df=meta42_bundle.phi_meta_df,
                coulton_bundle=cb,
                stratum=stratum,
                top_k_markers=config.top_k_markers,
                top_degs_df=Some(meta42_bundle.top_degs_df),
            ):
                case Failure(err):
                    return Failure(err)
                case Success(sim_df):
                    sim_frames.append(sim_df)

    full_sim_df = pl.concat(sim_frames)
    coulton_markers_df = extract_coulton_top_markers_df(
        bundle=rna_bundle,
        stratum=TissueStratum.TUMOR_ONLY,
        top_n=max(config.top_k_markers, config.bulk_signature_genes),
    )
    coulton_ici_df = evaluate_coulton_signatures_on_bulk_icb(
        coulton_markers_df=coulton_markers_df,
        bulk_cohorts=bulk_cohorts,
        n_sig_genes=config.bulk_signature_genes,
    )
    nearest_df, sensitivity_df, indication_df = build_nearest_substates_and_ici_tables(
        full_sim_df=full_sim_df,
        coulton_ici_df=coulton_ici_df,
        meta42_bundle=meta42_bundle,
    )
    centroids_df = serialize_coulton_centroids_to_polars(coulton_bundles)
    report_md = render_comparison_markdown_report(
        nearest_df=nearest_df,
        sensitivity_df=sensitivity_df,
        coulton_ici_df=coulton_ici_df,
        full_sim_df=full_sim_df,
        indication_df=indication_df,
    )

    return Success(
        ComparisonOutputs(
            coulton_centroids_df=centroids_df,
            coulton_top_markers_df=coulton_markers_df,
            full_similarity_df=full_sim_df,
            nearest_substates_df=nearest_df,
            sensitivity_summary_df=sensitivity_df,
            coulton_ici_outcomes_df=coulton_ici_df,
            meta42_indication_effects_df=indication_df,
            markdown_report=report_md,
        )
    )


def strip_coulton_cluster_number(label: str) -> str:
    """Strip leading numeric Seurat cluster prefix `<digits>_` from a Coulton substate label."""
    parts = label.split("_", 1)
    return parts[1] if len(parts) == 2 and parts[0].isdigit() else label


def _build_brackets_from_linkage(
    linkage_mat: FloatMatrix,
    leaf_perm: Sequence[int],
) -> tuple[DendrogramBracket, ...]:
    """Construct normalized U-brackets in leaf-center coordinates `0.5, 1.5, ..., K - 0.5`."""
    n_leaves = len(leaf_perm)
    if n_leaves <= 1 or linkage_mat.shape[0] == 0:
        return ()

    max_dist = max(float(np.max(linkage_mat[:, 2])), 1e-9)
    node_u: dict[int, float] = {
        int(orig_idx): float(pos) + 0.5 for pos, orig_idx in enumerate(leaf_perm)
    }
    node_h: dict[int, float] = {int(orig_idx): 0.0 for orig_idx in leaf_perm}
    brackets: list[DendrogramBracket] = []

    for step in range(linkage_mat.shape[0]):
        c_left = int(linkage_mat[step, 0])
        c_right = int(linkage_mat[step, 1])
        h_merge = float(linkage_mat[step, 2]) / max_dist
        u_l = node_u[c_left]
        u_r = node_u[c_right]
        h_l = node_h[c_left]
        h_r = node_h[c_right]
        new_node = n_leaves + step
        node_u[new_node] = 0.5 * (u_l + u_r)
        node_h[new_node] = h_merge
        brackets.append(
            DendrogramBracket(
                bracket_id=step,
                u_left=u_l,
                h_left=h_l,
                u_right=u_r,
                h_right=h_r,
                h_merge=h_merge,
            )
        )
    return tuple(brackets)


def compute_ward_biclustering(
    matrix: FloatMatrix,
    row_labels: Sequence[str],
    col_labels: Sequence[str],
) -> Result[BiclusterResult, str]:
    """Hierarchically bicluster `(rows, cols)` via Ward linkage with optimal leaf ordering."""
    n_rows, n_cols = matrix.shape
    if n_rows == 0 or n_cols == 0:
        return Failure("Cannot bicluster an empty matrix")
    if len(row_labels) != n_rows or len(col_labels) != n_cols:
        return Failure(
            f"Label length mismatch in biclustering: matrix={matrix.shape}, "
            f"rows={len(row_labels)}, cols={len(col_labels)}"
        )

    if n_rows >= 2:
        row_d = np.asarray(pdist(matrix, metric="euclidean"), dtype=np.float64)
        row_z: FloatMatrix = np.asarray(
            linkage(row_d, method="ward", optimal_ordering=True),
            dtype=np.float64,
        )
        row_perm = [int(i) for i in leaves_list(row_z)]
    else:
        row_z = np.zeros((0, 4), dtype=np.float64)
        row_perm = [0]

    if n_cols >= 2:
        col_d = np.asarray(pdist(matrix.T, metric="euclidean"), dtype=np.float64)
        col_z: FloatMatrix = np.asarray(
            linkage(col_d, method="ward", optimal_ordering=True),
            dtype=np.float64,
        )
        col_perm = [int(j) for j in leaves_list(col_z)]
    else:
        col_z = np.zeros((0, 4), dtype=np.float64)
        col_perm = [0]

    # Orient column tree along main diagonal (top-left to bottom-right block affinity)
    ordered_mat = matrix[np.ix_(row_perm, col_perm)]
    r_centered = np.arange(n_rows, dtype=np.float64) - 0.5 * float(n_rows - 1)
    c_centered = np.arange(n_cols, dtype=np.float64) - 0.5 * float(n_cols - 1)
    diag_cov = float(np.sum(np.outer(r_centered, c_centered) * ordered_mat))
    final_col_perm = list(reversed(col_perm)) if diag_cov < 0.0 else col_perm

    row_brackets = _build_brackets_from_linkage(row_z, row_perm)
    col_brackets = _build_brackets_from_linkage(col_z, final_col_perm)
    return Success(
        BiclusterResult(
            ordered_rows=tuple(row_labels[i] for i in row_perm),
            ordered_cols=tuple(col_labels[j] for j in final_col_perm),
            row_brackets=row_brackets,
            col_brackets=col_brackets,
        )
    )
