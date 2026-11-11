"""Pure functional core for Pan-Cancer T-Cell, CAF & Myeloid Triad atlas comparison and bulk ICB gating.

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
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import pdist
from scipy.stats import norm, rankdata

from tme_datasets.providers.caf_atlases import (
    BUECHLER_CAF_CLUSTER_SPECS,
    CORDS_CAF_CLUSTER_SPECS,
    GAO_CAF_CLUSTER_SPECS,
)
from tme_datasets.providers.tcell_atlases import (
    CHU_TCELL_CLUSTER_SPECS,
    CHU_TCELL_HS_MARKERS,
    ZHENG_TCELL_CLUSTER_SPECS,
)

from .models import (
    AtlasCentroidAccumulator,
    AtlasCentroidBundle,
    BiclusterResult,
    BoolVector,
    BulkIcbCohortBundle,
    DendrogramBracket,
    ExternalAtlasId,
    FloatMatrix,
    FloatVector,
    IciDirection,
    IntVector,
    META42_BIO_LABELS,
    META42_CAF_AND_MESENCHYMAL_CONTEXT_IDS,
    META42_CORE_CAF_IDS,
    META42_MYELOID_IDS,
    META42_T_NK_IDS,
    Meta42ReferenceBundle,
    PAN_CAF_CORE_GENES,
    PAN_MACROPHAGE_CORE_GENES,
    PAN_TCELL_CORE_GENES,
    PanAtlasConfig,
    PanAtlasOutputs,
    ProfileOrigin,
    SignificanceTier,
    StromalGatingResult,
    TriadCompartment,
    atlas_citation_label,
    atlas_short_prefix,
    atlas_target_compartment,
    classify_meta42_triad_compartment,
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
    "IGHV",
    "IGKV",
    "IGLV",
    "TRAV",
    "TRBV",
)

_ALLOWED_TCR_RECEPTORS: Final[frozenset[str]] = frozenset(
    {
        "TRAC",
        "TRBC1",
        "TRBC2",
        "TRDC",
        "TRGC1",
        "TRGC2",
        "TRDV1",
        "TRDV2",
        "TRGV9",
        "TRAV1-2",
    }
)


# =============================================================================
# 1. Substate Metadata Resolution & Informative Gene Filtering
# =============================================================================


def is_informative_hgnc_symbol(gene_sym: str) -> bool:
    """Filter out mitochondrial, ribosomal, lncRNA, and hypervariable V-segment symbols."""
    g = gene_sym.strip()
    if not g:
        return False
    if g in _ALLOWED_TCR_RECEPTORS:
        return True
    if g.startswith(_EXCLUDED_MARKER_PREFIXES):
        return False
    if g.startswith(("IGKC", "IGHG", "IGLC", "IGHA", "IGHM")):
        return False
    return (
        "-" not in g
        or g.startswith("HLA-")
        or g in {"NKX2-1", "WIPI2-AS1"}
    )


def resolve_atlas_substate_display_metadata(
    atlas_id: ExternalAtlasId,
    cluster_label: str,
) -> tuple[str, str, str]:
    """Return `(short_display_label, lineage_branch, functional_state)` for any external substate."""
    lbl = cluster_label.strip()
    prefix = atlas_short_prefix(atlas_id)
    match atlas_id:
        case ExternalAtlasId.GAO_CAF:
            if lbl in GAO_CAF_CLUSTER_SPECS:
                sp_g = GAO_CAF_CLUSTER_SPECS[lbl]
                return (
                    f"{prefix} {sp_g.subtype_label}",
                    sp_g.lineage_branch,
                    sp_g.functional_state,
                )
            return (f"{prefix} {lbl}", "CAF / Fibroblast", "Cross-tissue fibroblast substate")
        case ExternalAtlasId.CORDS_CAF:
            if lbl in CORDS_CAF_CLUSTER_SPECS:
                sp_c = CORDS_CAF_CLUSTER_SPECS[lbl]
                return (
                    f"{prefix} {sp_c.subtype_label}",
                    sp_c.lineage_branch,
                    sp_c.functional_state,
                )
            return (f"{prefix} {lbl}", "CAF / Stromal", "Pan-cancer tumor CAF substate")
        case ExternalAtlasId.BUECHLER_CAF:
            if lbl in BUECHLER_CAF_CLUSTER_SPECS:
                sp_b = BUECHLER_CAF_CLUSTER_SPECS[lbl]
                return (
                    f"{prefix} {sp_b.subtype_label}",
                    sp_b.lineage_branch,
                    sp_b.functional_state,
                )
            return (f"{prefix} {lbl}", "Perturbed-State Fibroblast", "FibroXplorer human cluster")
        case ExternalAtlasId.ZHENG_TCELL:
            if lbl in ZHENG_TCELL_CLUSTER_SPECS:
                sp_z = ZHENG_TCELL_CLUSTER_SPECS[lbl]
                cd_tok = lbl.split(".")[0]
                return (
                    f"{prefix} {cd_tok}.{sp_z.short_name}",
                    sp_z.lineage_branch,
                    sp_z.functional_state,
                )
            return (f"{prefix} {lbl}", "Tumor-Infiltrating T", "Zheng pan-cancer T cluster")
        case ExternalAtlasId.CHU_TCELL:
            if lbl in CHU_TCELL_CLUSTER_SPECS:
                sp_ch = CHU_TCELL_CLUSTER_SPECS[lbl]
                return (
                    f"{prefix} {sp_ch.short_name}",
                    sp_ch.lineage_branch,
                    sp_ch.functional_state,
                )
            return (f"{prefix} {lbl}", "Tumor-Infiltrating T", "Chu pan-cancer T state")
        case ExternalAtlasId.SHI_TCELL:
            return (f"{prefix} {lbl}", "Pan-Cancer T / NK", f"Shi et al. fine T/NK state ({lbl})")
        case ExternalAtlasId.SHI_CAF:
            return (f"{prefix} {lbl}", "Pan-Cancer Stromal / CAF", f"Shi et al. fine CAF state ({lbl})")
        case ExternalAtlasId.COULTON_TAM:
            parts = lbl.split("_", 1)
            short_tam = parts[1] if len(parts) == 2 and parts[0].isdigit() else lbl
            return (f"{prefix} {short_tam}", "Myeloid / TAM", f"Coulton et al. TAM cluster ({lbl})")
        case _ as unreachable:
            assert_never(unreachable)


def format_significance_stars(pval: float) -> SignificanceTier:
    """Map a two-sided p-value or FDR q-value onto standard publication significance stars."""
    if not math.isfinite(pval) or pval >= 0.05:
        return SignificanceTier.NS
    if pval < 1e-4:
        return SignificanceTier.FOUR_STARS
    if pval < 1e-3:
        return SignificanceTier.THREE_STARS
    if pval < 1e-2:
        return SignificanceTier.TWO_STARS
    return SignificanceTier.ONE_STAR


# =============================================================================
# 2. Out-of-Core Sparse CSR Batch Centroid Reduction & Adapters
# =============================================================================


def create_empty_atlas_accumulator(
    atlas_id: ExternalAtlasId,
    origin: ProfileOrigin,
    cluster_labels: Sequence[str],
    gene_symbols: Sequence[str],
) -> AtlasCentroidAccumulator:
    """Create a zero-initialized immutable accumulator for out-of-core cluster aggregation."""
    k_clusters = len(cluster_labels)
    g_genes = len(gene_symbols)
    return AtlasCentroidAccumulator(
        atlas_id=atlas_id,
        origin=origin,
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


def accumulate_atlas_batch(
    acc: AtlasCentroidAccumulator,
    batch_csr: sp.csr_matrix,
    cluster_codes: IntVector,
    is_tumor_mask: BoolVector,
) -> Result[AtlasCentroidAccumulator, str]:
    """Fold one sparse CSR batch of `log1p`-normalized cells into cluster sums without densifying zeros."""
    n_batch, n_genes = batch_csr.shape
    k_clusters = len(acc.cluster_labels)
    if n_genes != len(acc.gene_symbols):
        return Failure(
            f"Batch gene count {n_genes} does not match accumulator {len(acc.gene_symbols)} for {acc.atlas_id.value}"
        )
    if len(cluster_codes) != n_batch or len(is_tumor_mask) != n_batch:
        return Failure(f"Batch metadata vector lengths do not match CSR row count ({n_batch})")

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
    lin_data = np.expm1(np.clip(log_csr.data, 0.0, 20.0))
    lin_csr = sp.csr_matrix((lin_data, log_csr.indices, log_csr.indptr), shape=log_csr.shape)
    det_data = (log_csr.data > 0.0).astype(np.float64)
    det_csr = sp.csr_matrix((det_data, log_csr.indices, log_csr.indptr), shape=log_csr.shape)

    add_all_counts = np.bincount(clust_idx_all, minlength=k_clusters).astype(np.int64)
    add_tum_counts = np.bincount(clust_idx_tum, minlength=k_clusters).astype(np.int64)

    add_all_lin = np.asarray((ind_all @ lin_csr).toarray(), dtype=np.float64)
    add_tum_lin = np.asarray((ind_tum @ lin_csr).toarray(), dtype=np.float64)
    add_all_log = np.asarray((ind_all @ log_csr).toarray(), dtype=np.float64)
    add_tum_log = np.asarray((ind_tum @ log_csr).toarray(), dtype=np.float64)
    add_all_det = np.asarray((ind_all @ det_csr).toarray(), dtype=np.float64)
    add_tum_det = np.asarray((ind_tum @ det_csr).toarray(), dtype=np.float64)

    return Success(
        AtlasCentroidAccumulator(
            atlas_id=acc.atlas_id,
            origin=acc.origin,
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


def finalize_atlas_centroids(acc: AtlasCentroidAccumulator) -> AtlasCentroidBundle:
    """Convert accumulated sums into exact per-cell linear means, log1p means, and detection %."""
    denom_all = np.maximum(acc.all_cell_counts.reshape((-1, 1)).astype(np.float64), 1.0)
    denom_tum = np.maximum(acc.tumor_cell_counts.reshape((-1, 1)).astype(np.float64), 1.0)
    has_all = (acc.all_cell_counts > 0).reshape((-1, 1))
    has_tum = (acc.tumor_cell_counts > 0).reshape((-1, 1))

    all_det = np.asarray(acc.all_detect_counts, dtype=np.float64)
    tum_det = np.asarray(acc.tumor_detect_counts, dtype=np.float64)

    return AtlasCentroidBundle(
        atlas_id=acc.atlas_id,
        origin=acc.origin,
        cluster_labels=acc.cluster_labels,
        gene_symbols=acc.gene_symbols,
        all_cell_counts=acc.all_cell_counts,
        tumor_cell_counts=acc.tumor_cell_counts,
        all_mean_linear=np.where(has_all, acc.all_linear_sum / denom_all, 0.0),
        tumor_mean_linear=np.where(has_tum, acc.tumor_linear_sum / denom_tum, 0.0),
        all_mean_log=np.where(has_all, acc.all_log_sum / denom_all, 0.0),
        tumor_mean_log=np.where(has_tum, acc.tumor_log_sum / denom_tum, 0.0),
        all_pct_expressing=np.where(has_all, (all_det / denom_all) * 100.0, 0.0),
        tumor_pct_expressing=np.where(has_tum, (tum_det / denom_tum) * 100.0, 0.0),
    )


def deduplicate_atlas_gene_columns(
    raw_symbols: Sequence[str],
    bundle: AtlasCentroidBundle,
) -> AtlasCentroidBundle:
    """Aggregate any duplicate HGNC gene symbols cleanly by linear addition."""
    unique_map: dict[str, list[int]] = {}
    ordered_syms: list[str] = []
    for idx, sym in enumerate(raw_symbols):
        cleaned = sym.strip()
        if not cleaned or cleaned.upper() in {"NAN", "NONE", "<NA>", "NA"}:
            continue
        if cleaned not in unique_map:
            unique_map[cleaned] = []
            ordered_syms.append(cleaned)
        unique_map[cleaned].append(idx)

    if len(ordered_syms) == len(raw_symbols) and all(len(v) == 1 for v in unique_map.values()):
        return AtlasCentroidBundle(
            atlas_id=bundle.atlas_id,
            origin=bundle.origin,
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

    return AtlasCentroidBundle(
        atlas_id=bundle.atlas_id,
        origin=bundle.origin,
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


def filter_atlas_min_cells(
    bundle: AtlasCentroidBundle,
    min_cells: int = 25,
) -> AtlasCentroidBundle:
    """Filter out sparse / unannotated clusters with fewer than `min_cells` total cells."""
    if bundle.origin == ProfileOrigin.CURATED_SIGNATURE_PRIOR:
        return bundle
    keep_idx = [
        idx
        for idx, cnt in enumerate(bundle.all_cell_counts)
        if int(cnt) >= min_cells
        and bundle.cluster_labels[idx].strip().lower() not in {"na", "unassigned", "unknown", "23_na"}
    ]
    if not keep_idx or len(keep_idx) == len(bundle.cluster_labels):
        return bundle
    return AtlasCentroidBundle(
        atlas_id=bundle.atlas_id,
        origin=bundle.origin,
        cluster_labels=tuple(bundle.cluster_labels[i] for i in keep_idx),
        gene_symbols=bundle.gene_symbols,
        all_cell_counts=bundle.all_cell_counts[keep_idx],
        tumor_cell_counts=bundle.tumor_cell_counts[keep_idx],
        all_mean_linear=bundle.all_mean_linear[keep_idx, :],
        tumor_mean_linear=bundle.tumor_mean_linear[keep_idx, :],
        all_mean_log=bundle.all_mean_log[keep_idx, :],
        tumor_mean_log=bundle.tumor_mean_log[keep_idx, :],
        all_pct_expressing=bundle.all_pct_expressing[keep_idx, :],
        tumor_pct_expressing=bundle.tumor_pct_expressing[keep_idx, :],
    )


def build_chu_curated_centroid_bundle(
    gene_universe: Sequence[str],
) -> AtlasCentroidBundle:
    """Construct `Chu_TCell_Atlas` curated marker-prior centroids over `gene_universe`.

    Used strictly when the MD Anderson portal Seurat `.rds` is not locally staged, transparently
    tagged with `ProfileOrigin.CURATED_SIGNATURE_PRIOR` so all 39 Chu et al. (*Nature Medicine* 2023)
    T-cell states—including heat-shock `CD8_c4_Tstr` and `CD4_c4_Tstr`—participate in cross-atlas
    ontology alignment and clinical bulk ICB meta-analysis.
    """
    genes = tuple(dict.fromkeys(g.strip() for g in gene_universe if g.strip()))
    g_map = {g: i for i, g in enumerate(genes)}
    clusters = tuple(CHU_TCELL_CLUSTER_SPECS.keys())
    k_c = len(clusters)
    g_n = len(genes)

    lin_mat = np.ones((k_c, g_n), dtype=np.float64)
    pct_mat = np.full((k_c, g_n), 1.5, dtype=np.float64)
    present_pan_idx = [g_map[pg] for pg in PAN_TCELL_CORE_GENES if pg in g_map]
    pan_idx_set = set(present_pan_idx)

    for row_i, cid in enumerate(clusters):
        spec = CHU_TCELL_CLUSTER_SPECS[cid]
        marker_weights: dict[int, float] = {}

        if cid.startswith("CD8"):
            for lg in ("CD8A", "CD8B", "NKG7", "CCL5"):
                if lg in g_map and g_map[lg] not in pan_idx_set:
                    marker_weights[g_map[lg]] = 3.0
                    pct_mat[row_i, g_map[lg]] = 70.0
        elif cid.startswith("CD4") or "Treg" in cid or "Th2" in cid:
            for lg in ("CD4", "IL7R", "CD40LG"):
                if lg in g_map and g_map[lg] not in pan_idx_set:
                    marker_weights[g_map[lg]] = 2.8
                    pct_mat[row_i, g_map[lg]] = 65.0

        if "Treg" in cid:
            for rg in ("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"):
                if rg in g_map and g_map[rg] not in pan_idx_set:
                    marker_weights[g_map[rg]] = 5.5
                    pct_mat[row_i, g_map[rg]] = 78.0

        if cid in ("CD8_c4_Tstr", "CD4_c4_Tstr"):
            for hsg in CHU_TCELL_HS_MARKERS:
                if hsg in g_map and g_map[hsg] not in pan_idx_set:
                    marker_weights[g_map[hsg]] = 6.5
                    pct_mat[row_i, g_map[hsg]] = 82.0

        for rank_idx, mk in enumerate(spec.canonical_markers):
            if mk in g_map and g_map[mk] not in pan_idx_set:
                boost = max(9.5 - 0.7 * float(rank_idx), 4.5)
                marker_weights[g_map[mk]] = max(marker_weights.get(g_map[mk], 0.0), boost)
                pct_mat[row_i, g_map[mk]] = 85.0

        # Allocate an identical total marker budget (300.0) across active subtype markers and an
        # identical pan-T core budget (80.0) across PAN_TCELL_CORE_GENES so every cluster has the
        # exact same row sum and untouched background genes have zero variance across clusters.
        w_sum = sum(marker_weights.values())
        if w_sum > 1e-9:
            for g_i, w_val in marker_weights.items():
                lin_mat[row_i, g_i] = 1.0 + 300.0 * (w_val / w_sum)
        if present_pan_idx:
            per_pan = 80.0 / float(len(present_pan_idx))
            for g_i in present_pan_idx:
                lin_mat[row_i, g_i] = 1.0 + per_pan
                pct_mat[row_i, g_i] = 75.0

    log_mat = np.log1p(lin_mat)
    pseudo_counts = np.full(k_c, 1_000, dtype=np.int64)
    return AtlasCentroidBundle(
        atlas_id=ExternalAtlasId.CHU_TCELL,
        origin=ProfileOrigin.CURATED_SIGNATURE_PRIOR,
        cluster_labels=clusters,
        gene_symbols=genes,
        all_cell_counts=pseudo_counts,
        tumor_cell_counts=pseudo_counts,
        all_mean_linear=lin_mat,
        tumor_mean_linear=lin_mat.copy(),
        all_mean_log=log_mat,
        tumor_mean_log=log_mat.copy(),
        all_pct_expressing=pct_mat,
        tumor_pct_expressing=pct_mat.copy(),
    )


def adapt_coulton_cached_parquet_to_bundle(
    coulton_df: pl.DataFrame,
) -> Result[AtlasCentroidBundle, str]:
    """Convert cached `coulton_substate_centroids.parquet` (`RNA` assay) into an `AtlasCentroidBundle`."""
    if coulton_df.is_empty():
        return Failure("Coulton centroid DataFrame is empty")
    rna_sub = (
        coulton_df.filter(pl.col("assay") == "RNA")
        if "assay" in coulton_df.columns
        else coulton_df
    )
    if rna_sub.is_empty():
        return Failure("No RNA assay rows found in Coulton centroid DataFrame")

    sub = rna_sub.sort(["cluster_idx", "gene_idx"])
    clust_df = (
        sub.select(["cluster_idx", "cluster_label", "all_cell_count", "tumor_cell_count"])
        .unique(subset=["cluster_idx"], maintain_order=True)
        .sort("cluster_idx")
    )
    gene_df = (
        sub.select(["gene_idx", "gene_symbol"])
        .unique(subset=["gene_idx"], maintain_order=True)
        .sort("gene_idx")
    )
    k_c = clust_df.height
    g_n = gene_df.height
    if sub.height != k_c * g_n:
        return Failure(
            f"Malformed Coulton centroid table: expected {k_c * g_n} rows, got {sub.height}"
        )

    raw_bundle = AtlasCentroidBundle(
        atlas_id=ExternalAtlasId.COULTON_TAM,
        origin=ProfileOrigin.EMPIRICAL_SINGLE_CELL_H5AD,
        cluster_labels=tuple(clust_df["cluster_label"].cast(pl.String).to_list()),
        gene_symbols=tuple(gene_df["gene_symbol"].cast(pl.String).to_list()),
        all_cell_counts=clust_df["all_cell_count"].to_numpy().astype(np.int64),
        tumor_cell_counts=clust_df["tumor_cell_count"].to_numpy().astype(np.int64),
        all_mean_linear=sub["all_mean_linear"].to_numpy().astype(np.float64).reshape((k_c, g_n)),
        tumor_mean_linear=sub["tumor_mean_linear"].to_numpy().astype(np.float64).reshape((k_c, g_n)),
        all_mean_log=sub["all_mean_log"].to_numpy().astype(np.float64).reshape((k_c, g_n)),
        tumor_mean_log=sub["tumor_mean_log"].to_numpy().astype(np.float64).reshape((k_c, g_n)),
        all_pct_expressing=sub["all_pct_expressing"].to_numpy().astype(np.float64).reshape((k_c, g_n)),
        tumor_pct_expressing=sub["tumor_pct_expressing"].to_numpy().astype(np.float64).reshape((k_c, g_n)),
    )
    return Success(filter_atlas_min_cells(raw_bundle, min_cells=10))


def serialize_atlas_centroids_to_polars(
    bundles: Sequence[AtlasCentroidBundle],
) -> pl.DataFrame:
    """Serialize `AtlasCentroidBundle` sequence into a tidy, lossless Polars DataFrame."""
    frames: list[pl.DataFrame] = []
    for b in bundles:
        k_clusters = len(b.cluster_labels)
        g_genes = len(b.gene_symbols)
        if k_clusters == 0 or g_genes == 0:
            continue
        frames.append(
            pl.DataFrame(
                {
                    "atlas_id": np.full(k_clusters * g_genes, b.atlas_id.value, dtype=object),
                    "origin": np.full(k_clusters * g_genes, b.origin.value, dtype=object),
                    "cluster_idx": np.repeat(np.arange(k_clusters, dtype=np.int32), g_genes),
                    "cluster_label": np.repeat(np.asarray(b.cluster_labels, dtype=object), g_genes),
                    "all_cell_count": np.repeat(b.all_cell_counts, g_genes),
                    "tumor_cell_count": np.repeat(b.tumor_cell_counts, g_genes),
                    "gene_idx": np.tile(np.arange(g_genes, dtype=np.int32), k_clusters),
                    "gene_symbol": np.tile(np.asarray(b.gene_symbols, dtype=object), k_clusters),
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


def deserialize_atlas_centroids_from_polars(
    df: pl.DataFrame,
) -> Result[tuple[AtlasCentroidBundle, ...], str]:
    """Reconstruct `AtlasCentroidBundle` sequence from a cached Polars DataFrame."""
    if df.is_empty():
        return Failure("Cached external atlas centroid DataFrame is empty")

    bundles: list[AtlasCentroidBundle] = []
    for raw_aid in df["atlas_id"].unique(maintain_order=True).to_list():
        try:
            aid_enum = ExternalAtlasId(str(raw_aid))
        except ValueError:
            return Failure(f"Unknown ExternalAtlasId in centroid cache: {raw_aid}")

        sub = df.filter(pl.col("atlas_id") == raw_aid).sort(["cluster_idx", "gene_idx"])
        raw_orig = str(sub["origin"][0]) if "origin" in sub.columns else ProfileOrigin.EMPIRICAL_SINGLE_CELL_H5AD.value
        orig_enum = (
            ProfileOrigin.CURATED_SIGNATURE_PRIOR
            if raw_orig == ProfileOrigin.CURATED_SIGNATURE_PRIOR.value
            else ProfileOrigin.EMPIRICAL_SINGLE_CELL_H5AD
        )
        clust_df = (
            sub.select(["cluster_idx", "cluster_label", "all_cell_count", "tumor_cell_count"])
            .unique(subset=["cluster_idx"], maintain_order=True)
            .sort("cluster_idx")
        )
        gene_df = (
            sub.select(["gene_idx", "gene_symbol"])
            .unique(subset=["gene_idx"], maintain_order=True)
            .sort("gene_idx")
        )
        k_c = clust_df.height
        g_n = gene_df.height
        if sub.height != k_c * g_n:
            return Failure(
                f"Corrupted centroid cache for {raw_aid}: expected {k_c * g_n} rows, got {sub.height}"
            )

        bundles.append(
            AtlasCentroidBundle(
                atlas_id=aid_enum,
                origin=orig_enum,
                cluster_labels=tuple(clust_df["cluster_label"].cast(pl.String).to_list()),
                gene_symbols=tuple(gene_df["gene_symbol"].cast(pl.String).to_list()),
                all_cell_counts=clust_df["all_cell_count"].to_numpy().astype(np.int64),
                tumor_cell_counts=clust_df["tumor_cell_count"].to_numpy().astype(np.int64),
                all_mean_linear=sub["all_mean_linear"].to_numpy().astype(np.float64).reshape((k_c, g_n)),
                tumor_mean_linear=sub["tumor_mean_linear"].to_numpy().astype(np.float64).reshape((k_c, g_n)),
                all_mean_log=sub["all_mean_log"].to_numpy().astype(np.float64).reshape((k_c, g_n)),
                tumor_mean_log=sub["tumor_mean_log"].to_numpy().astype(np.float64).reshape((k_c, g_n)),
                all_pct_expressing=sub["all_pct_expressing"].to_numpy().astype(np.float64).reshape((k_c, g_n)),
                tumor_pct_expressing=sub["tumor_pct_expressing"].to_numpy().astype(np.float64).reshape((k_c, g_n)),
            )
        )
    return Success(tuple(bundles))


def select_effective_bundle_matrices(
    bundle: AtlasCentroidBundle,
    min_tumor_cells: int = 15,
) -> tuple[FloatMatrix, FloatMatrix, IntVector]:
    """Select Tumor-Only expression per cluster when `>= min_tumor_cells`, falling back to All-Tissues.

    This preserves tumor-infiltrating specificity for tumor-enriched states (`LRRC15+ myCAF`,
    `CXCL13+ Tex`) while gracefully retaining normal/cross-tissue reference states (`PI16+`
    universal progenitor fibroblasts, naive blood/adjacent-normal T cells).
    """
    use_tum = (bundle.tumor_cell_counts >= min_tumor_cells).reshape((-1, 1))
    lin_eff = np.where(use_tum, bundle.tumor_mean_linear, bundle.all_mean_linear)
    pct_eff = np.where(use_tum, bundle.tumor_pct_expressing, bundle.all_pct_expressing)
    cnt_eff = np.where(bundle.tumor_cell_counts >= min_tumor_cells, bundle.tumor_cell_counts, bundle.all_cell_counts)
    return (
        np.asarray(lin_eff, dtype=np.float64),
        np.asarray(pct_eff, dtype=np.float64),
        np.asarray(cnt_eff, dtype=np.int64),
    )


# =============================================================================
# 3. Specificity Algebra, Stage 3c Hybrid Similarity & Cross-Atlas Equivalence
# =============================================================================


def compute_centered_log2_specificity(
    density_mat: FloatMatrix,
    reference_row_indices: Sequence[int],
    eps: float = 5e-5,
) -> FloatMatrix:
    """Compute gene-centered log2 relative density specificity relative to `reference_row_indices`."""
    row_sums = np.maximum(density_mat.sum(axis=1, keepdims=True), 1e-12)
    normed = density_mat / row_sums
    log_mat = np.log2(normed + eps)
    ref_rows = list(reference_row_indices) if reference_row_indices else list(range(density_mat.shape[0]))
    ref_mean = log_mat[ref_rows, :].mean(axis=0, keepdims=True)
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


def extract_atlas_top_markers_df(
    bundles: Sequence[AtlasCentroidBundle],
    top_n: int = 25,
) -> pl.DataFrame:
    """Extract top positively enriched protein-coding markers per substate across all external atlases."""
    rows: list[dict[str, object]] = []
    eps = 5e-5

    for bundle in bundles:
        k_clusters = len(bundle.cluster_labels)
        if k_clusters == 0 or len(bundle.gene_symbols) == 0:
            continue

        lin_mat, pct_mat, eff_counts = select_effective_bundle_matrices(bundle)
        row_sums = np.maximum(lin_mat.sum(axis=1, keepdims=True), 1e-12)
        dens_mat = lin_mat / row_sums
        inf_mask = np.array([is_informative_hgnc_symbol(g) for g in bundle.gene_symbols], dtype=bool)
        target_comp = atlas_target_compartment(bundle.atlas_id)

        for k_idx, c_lbl in enumerate(bundle.cluster_labels):
            other_bg = [i for i in range(k_clusters) if i != k_idx]
            bg_vec = dens_mat[other_bg, :].mean(axis=0) if other_bg else dens_mat[k_idx, :]
            lfc_vec = np.log2((dens_mat[k_idx, :] + eps) / (bg_vec + eps))

            min_pct = 2.0 if bundle.origin == ProfileOrigin.CURATED_SIGNATURE_PRIOR else 10.0
            eligible = np.flatnonzero(inf_mask & (pct_mat[k_idx, :] >= min_pct))
            if len(eligible) < top_n:
                eligible = np.flatnonzero(inf_mask & (pct_mat[k_idx, :] >= 3.0))
            if len(eligible) < top_n:
                eligible = np.flatnonzero(inf_mask)
            ranked = eligible[np.argsort(-lfc_vec[eligible])][:top_n]

            disp_lbl, branch_lbl, func_lbl = resolve_atlas_substate_display_metadata(
                bundle.atlas_id,
                c_lbl,
            )
            for rank_pos, g_idx in enumerate(ranked, start=1):
                rows.append(
                    {
                        "atlas_id": bundle.atlas_id.value,
                        "origin": bundle.origin.value,
                        "target_compartment": target_comp.value,
                        "cluster_label": c_lbl,
                        "substate_display": disp_lbl,
                        "lineage_branch": branch_lbl,
                        "functional_state": func_lbl,
                        "all_cell_count": int(bundle.all_cell_counts[k_idx]),
                        "tumor_cell_count": int(bundle.tumor_cell_counts[k_idx]),
                        "effective_cell_count": int(eff_counts[k_idx]),
                        "marker_rank": rank_pos,
                        "gene_symbol": bundle.gene_symbols[int(g_idx)],
                        "specificity_log2fc": float(lfc_vec[int(g_idx)]),
                        "pct_expressing": float(pct_mat[k_idx, int(g_idx)]),
                        "mean_linear_expr": float(lin_mat[k_idx, int(g_idx)]),
                    }
                )
    return pl.DataFrame(rows)


def _cognate_meta42_ids_for_compartment(comp: TriadCompartment) -> tuple[str, ...]:
    """Return the cognate Meta42 reference state IDs for within-compartment centering."""
    match comp:
        case TriadCompartment.T_NK:
            return META42_T_NK_IDS
        case TriadCompartment.STROMAL_CAF:
            return META42_CAF_AND_MESENCHYMAL_CONTEXT_IDS
        case TriadCompartment.MYELOID_TAM:
            return META42_MYELOID_IDS
        case (
            TriadCompartment.B_PLASMA
            | TriadCompartment.ENDOTHELIAL
            | TriadCompartment.EPITHELIAL_TUMOR
            | TriadCompartment.OTHER
        ):
            return tuple(META42_BIO_LABELS.keys())
        case _ as unreachable:
            assert_never(unreachable)


def compute_meta42_vs_atlas_similarity_df(
    phi_meta_df: pl.DataFrame,
    bundle: AtlasCentroidBundle,
    top_k_markers: int = 25,
    top_degs_df: Maybe[pl.DataFrame] = Nothing,
) -> Result[pl.DataFrame, str]:
    """Compute Stage 3c Marker-Anchored Hybrid similarity between 42 Meta42 states and an external atlas."""
    meta_states = tuple(phi_meta_df["cluster"].cast(pl.String).to_list())
    meta_gene_cols = tuple(c for c in phi_meta_df.columns if c != "cluster")
    meta_gene_map = {g: idx for idx, g in enumerate(meta_gene_cols)}
    ext_gene_map = {g: idx for idx, g in enumerate(bundle.gene_symbols)}

    shared_genes = tuple(sorted(set(meta_gene_cols).intersection(set(bundle.gene_symbols))))
    if len(shared_genes) < 30:
        return Failure(
            f"Insufficient shared genes ({len(shared_genes)}) between Consensus-Meta42 and {bundle.atlas_id.value}"
        )
    shared_to_idx = {g: idx for idx, g in enumerate(shared_genes)}

    m_col_idx = [meta_gene_map[g] for g in shared_genes]
    e_col_idx = [ext_gene_map[g] for g in shared_genes]

    raw_meta_mat = phi_meta_df.select(meta_gene_cols).to_numpy().astype(np.float64)[:, m_col_idx]
    lin_ext_full, pct_ext_full, eff_counts = select_effective_bundle_matrices(bundle)
    raw_ext_mat = lin_ext_full[:, e_col_idx]
    pct_ext_sub = pct_ext_full[:, e_col_idx]

    meta_dens = raw_meta_mat / np.maximum(raw_meta_mat.sum(axis=1, keepdims=True), 1e-12)
    ext_dens = raw_ext_mat / np.maximum(raw_ext_mat.sum(axis=1, keepdims=True), 1e-12)

    target_comp = atlas_target_compartment(bundle.atlas_id)
    cognate_ids = set(_cognate_meta42_ids_for_compartment(target_comp))
    comp_meta_idx = [i for i, st in enumerate(meta_states) if st in cognate_ids]
    all_meta_idx = list(range(len(meta_states)))
    if not comp_meta_idx:
        comp_meta_idx = all_meta_idx
    all_ext_idx = list(range(len(bundle.cluster_labels)))

    spec_meta_comp = compute_centered_log2_specificity(meta_dens, comp_meta_idx)
    spec_meta_all = compute_centered_log2_specificity(meta_dens, all_meta_idx)
    spec_ext = compute_centered_log2_specificity(ext_dens, all_ext_idx)

    # Pan-lineage reference frame: center external log2 density against the whole-tumor 42-lineage
    # background (rather than subtracting the external atlas's own single-lineage mean, which would
    # cancel pan-T, pan-CAF, or pan-macrophage identity genes!).
    meta_log_all = np.log2(meta_dens + 5e-5)
    ext_log_all = np.log2(ext_dens + 5e-5)
    pan42_log_bg = meta_log_all.mean(axis=0, keepdims=True)
    pan42_log_std = np.maximum(meta_log_all.std(axis=0, keepdims=True), 1e-4)
    spec_ext_pan = np.asarray(ext_log_all - pan42_log_bg, dtype=np.float64)
    z_ext_pan = spec_ext_pan / pan42_log_std

    cos_within_comp = row_cosine_similarity(spec_meta_comp, spec_ext)
    cos_pan_lineage = row_cosine_similarity(spec_meta_all, spec_ext_pan)
    pearson_pan_log = row_pearson_correlation(meta_log_all, ext_log_all)

    inf_flags = np.array([is_informative_hgnc_symbol(g) for g in shared_genes], dtype=bool)
    ext_var = spec_ext.std(axis=0)
    comp_var = spec_meta_comp[comp_meta_idx, :].std(axis=0)
    combined_var = np.where(inf_flags, ext_var + comp_var, -1.0)
    n_hvg = min(750, int(np.sum(inf_flags)))
    hvg_idx = np.argsort(-combined_var)[: max(n_hvg, 10)]
    spearman_hvg = row_spearman_correlation(spec_meta_comp[:, hvg_idx], spec_ext[:, hvg_idx])

    z_meta_cols = (spec_meta_comp - spec_meta_comp[comp_meta_idx, :].mean(axis=0, keepdims=True)) / np.maximum(
        spec_meta_comp[comp_meta_idx, :].std(axis=0, keepdims=True), 1e-6
    )
    z_ext_cols = (spec_ext - spec_ext.mean(axis=0, keepdims=True)) / np.maximum(
        spec_ext.std(axis=0, keepdims=True), 1e-6
    )

    canonical_deg_z = np.zeros((len(meta_states), len(bundle.cluster_labels)), dtype=np.float64)
    pan_deg_z = np.zeros((len(meta_states), len(bundle.cluster_labels)), dtype=np.float64)
    canonical_gene_sets: list[set[int]] = [set() for _ in meta_states]
    has_canonical_degs = False

    match top_degs_df:
        case Some(degs_table) if (
            not degs_table.is_empty()
            and {"cell_state", "gene", "specificity_log2fc"}.issubset(set(degs_table.columns))
        ):
            has_canonical_degs = True
            for m_i, m_st in enumerate(meta_states):
                st_rows = degs_table.filter(pl.col("cell_state") == m_st).sort(
                    "specificity_log2fc", descending=True
                )
                matched_g: list[int] = []
                matched_w: list[float] = []
                for deg_row in st_rows.iter_rows(named=True):
                    g_sym = str(deg_row["gene"]).strip()
                    if g_sym in shared_to_idx:
                        matched_g.append(shared_to_idx[g_sym])
                        matched_w.append(max(float(deg_row["specificity_log2fc"]), 0.25))
                if matched_g:
                    canonical_gene_sets[m_i] = set(matched_g)
                    w_arr = np.asarray(matched_w, dtype=np.float64)
                    w_denom = max(float(w_arr.sum()), 1e-6)
                    canonical_deg_z[m_i, :] = (z_ext_cols[:, matched_g] @ w_arr) / w_denom
                    pan_deg_z[m_i, :] = (z_ext_pan[:, matched_g] @ w_arr) / w_denom
        case _:
            has_canonical_degs = False

    inf_indices = np.flatnonzero(inf_flags)
    k_eff = min(top_k_markers, len(inf_indices))
    k_overlap = min(60, len(inf_indices))

    meta_top_k: list[npt.NDArray[np.intp]] = []
    meta_top_60_sets: list[set[int]] = []
    for m_i in range(len(meta_states)):
        ord_m = inf_indices[np.argsort(-spec_meta_comp[m_i, inf_indices])]
        canon_first = sorted(canonical_gene_sets[m_i], key=lambda gi: -spec_meta_comp[m_i, gi])
        merged_top = canon_first + [int(x) for x in ord_m if int(x) not in canonical_gene_sets[m_i]]
        meta_top_k.append(np.asarray(merged_top[:k_eff], dtype=np.intp))
        meta_top_60_sets.append(set(merged_top[:k_overlap]).union(canonical_gene_sets[m_i]))

    ext_top_k: list[npt.NDArray[np.intp]] = []
    ext_top_60_sets: list[set[int]] = []
    min_pct_cand = 2.0 if bundle.origin == ProfileOrigin.CURATED_SIGNATURE_PRIOR else 10.0
    for e_i in range(len(bundle.cluster_labels)):
        cand_idx = np.flatnonzero(inf_flags & (pct_ext_sub[e_i, :] >= min_pct_cand))
        if len(cand_idx) < k_eff:
            cand_idx = np.flatnonzero(inf_flags & (pct_ext_sub[e_i, :] >= 3.0))
        if len(cand_idx) < k_eff:
            cand_idx = inf_indices
        ord_e = cand_idx[np.argsort(-spec_ext[e_i, cand_idx])]
        ext_top_k.append(ord_e[:k_eff])
        ext_top_60_sets.append(set(int(x) for x in ord_e[:k_overlap]))

    bidir_z = np.zeros((len(meta_states), len(bundle.cluster_labels)), dtype=np.float64)
    shared_genes_str: list[list[str]] = [
        ["" for _ in bundle.cluster_labels] for _ in meta_states
    ]
    shared_counts = np.zeros((len(meta_states), len(bundle.cluster_labels)), dtype=np.int32)

    for m_i in range(len(meta_states)):
        m_vec = meta_top_k[m_i]
        c_set_m = canonical_gene_sets[m_i]
        for e_i in range(len(bundle.cluster_labels)):
            e_vec = ext_top_k[e_i]
            z_m_in_e = float(z_ext_cols[e_i, m_vec].mean()) if len(m_vec) > 0 else 0.0
            z_e_in_m = float(z_meta_cols[m_i, e_vec].mean()) if len(e_vec) > 0 else 0.0
            bidir_z[m_i, e_i] = 0.5 * (z_m_in_e + z_e_in_m)

            ov_idx = sorted(
                meta_top_60_sets[m_i].intersection(ext_top_60_sets[e_i]),
                key=lambda gi: -(
                    (4.0 if gi in c_set_m else 0.0) + spec_meta_comp[m_i, gi] + spec_ext[e_i, gi]
                ),
            )
            shared_counts[m_i, e_i] = len(ov_idx)
            shared_genes_str[m_i][e_i] = ", ".join(shared_genes[gi] for gi in ov_idx[:8])

    composite_mat = (
        (
            0.40 * cos_within_comp
            + 0.25 * spearman_hvg
            + 0.20 * np.tanh(bidir_z / 2.5)
            + 0.15 * np.tanh(canonical_deg_z / 2.0)
        )
        if has_canonical_degs
        else (
            0.50 * cos_within_comp
            + 0.30 * spearman_hvg
            + 0.20 * np.tanh(bidir_z / 3.0)
        )
    )
    hybrid_pan_affinity = (
        0.65 * cos_pan_lineage + 0.35 * np.tanh(pan_deg_z / 2.5)
        if has_canonical_degs
        else cos_pan_lineage
    )

    cite_lbl = atlas_citation_label(bundle.atlas_id)
    rows: list[dict[str, object]] = []
    for m_i, m_st in enumerate(meta_states):
        rank_order = np.argsort(-composite_mat[m_i, :])
        ranks = {int(pos): rk for rk, pos in enumerate(rank_order, start=1)}
        m_label = META42_BIO_LABELS.get(m_st, m_st)
        m_comp = classify_meta42_triad_compartment(m_st)

        for e_i, e_lbl in enumerate(bundle.cluster_labels):
            disp_lbl, branch_lbl, func_lbl = resolve_atlas_substate_display_metadata(
                bundle.atlas_id,
                e_lbl,
            )
            rows.append(
                {
                    "atlas_id": bundle.atlas_id.value,
                    "atlas_citation": cite_lbl,
                    "origin": bundle.origin.value,
                    "atlas_compartment": target_comp.value,
                    "n_shared_genes": len(shared_genes),
                    "meta_state": m_st,
                    "meta_state_label": m_label,
                    "meta_compartment": m_comp.value,
                    "is_cognate_compartment": m_comp == target_comp,
                    "external_cluster": e_lbl,
                    "external_display": disp_lbl,
                    "external_branch": branch_lbl,
                    "external_functional_state": func_lbl,
                    "external_cell_count": int(eff_counts[e_i]),
                    "within_compartment_spec_cosine": float(cos_within_comp[m_i, e_i]),
                    "within_compartment_spec_spearman": float(spearman_hvg[m_i, e_i]),
                    "bidirectional_marker_z": float(bidir_z[m_i, e_i]),
                    "canonical_deg_z": float(canonical_deg_z[m_i, e_i]),
                    "pan_lineage_log_pearson": float(pearson_pan_log[m_i, e_i]),
                    "pan_lineage_spec_cosine": float(cos_pan_lineage[m_i, e_i]),
                    "hybrid_pan_affinity": float(hybrid_pan_affinity[m_i, e_i]),
                    "composite_substate_score": float(composite_mat[m_i, e_i]),
                    "rank_within_atlas": ranks[e_i],
                    "n_shared_top60_markers": int(shared_counts[m_i, e_i]),
                    "shared_signature_genes": shared_genes_str[m_i][e_i],
                }
            )
    return Success(pl.DataFrame(rows))


def compute_cross_atlas_equivalence_df(
    bundles: Sequence[AtlasCentroidBundle],
    top_k_markers: int = 25,
) -> pl.DataFrame:
    """Compute top-1 cross-atlas substate equivalences within T/NK and Stromal/CAF compartments."""
    rows: list[dict[str, object]] = []
    for comp in (TriadCompartment.T_NK, TriadCompartment.STROMAL_CAF):
        comp_bundles = [b for b in bundles if atlas_target_compartment(b.atlas_id) == comp]
        for i_a, b_a in enumerate(comp_bundles):
            for i_b, b_b in enumerate(comp_bundles):
                if i_a >= i_b:
                    continue
                shared = tuple(sorted(set(b_a.gene_symbols).intersection(set(b_b.gene_symbols))))
                if len(shared) < 50:
                    continue
                map_a = {g: idx for idx, g in enumerate(b_a.gene_symbols)}
                map_b = {g: idx for idx, g in enumerate(b_b.gene_symbols)}
                idx_a = [map_a[g] for g in shared]
                idx_b = [map_b[g] for g in shared]

                lin_a, pct_a_full, _ = select_effective_bundle_matrices(b_a)
                lin_b, pct_b_full, _ = select_effective_bundle_matrices(b_b)
                sub_a = lin_a[:, idx_a]
                sub_b = lin_b[:, idx_b]
                pct_a = pct_a_full[:, idx_a]
                pct_b = pct_b_full[:, idx_b]
                dens_a = sub_a / np.maximum(sub_a.sum(axis=1, keepdims=True), 1e-12)
                dens_b = sub_b / np.maximum(sub_b.sum(axis=1, keepdims=True), 1e-12)

                spec_a = compute_centered_log2_specificity(dens_a, list(range(dens_a.shape[0])))
                spec_b = compute_centered_log2_specificity(dens_b, list(range(dens_b.shape[0])))

                inf_flags = np.array([is_informative_hgnc_symbol(g) for g in shared], dtype=bool)
                inf_idx = np.flatnonzero(inf_flags)
                hvg = inf_idx[np.argsort(-(spec_a[:, inf_idx].std(axis=0) + spec_b[:, inf_idx].std(axis=0)))][
                    : min(600, len(inf_idx))
                ]

                cos_ab = row_cosine_similarity(spec_a, spec_b)
                rho_ab = row_spearman_correlation(spec_a[:, hvg], spec_b[:, hvg])
                score_ab = 0.60 * cos_ab + 0.40 * rho_ab

                k_ov = min(top_k_markers * 2, len(inf_idx))
                min_p_a = 2.0 if b_a.origin == ProfileOrigin.CURATED_SIGNATURE_PRIOR else 10.0
                min_p_b = 2.0 if b_b.origin == ProfileOrigin.CURATED_SIGNATURE_PRIOR else 10.0
                top_sets_a: list[set[int]] = []
                for r in range(spec_a.shape[0]):
                    cand_a = np.flatnonzero(inf_flags & (pct_a[r, :] >= min_p_a))
                    if len(cand_a) < k_ov:
                        cand_a = inf_idx
                    top_sets_a.append(set(int(x) for x in cand_a[np.argsort(-spec_a[r, cand_a])][:k_ov]))
                top_sets_b: list[set[int]] = []
                for c in range(spec_b.shape[0]):
                    cand_b = np.flatnonzero(inf_flags & (pct_b[c, :] >= min_p_b))
                    if len(cand_b) < k_ov:
                        cand_b = inf_idx
                    top_sets_b.append(set(int(x) for x in cand_b[np.argsort(-spec_b[c, cand_b])][:k_ov]))

                for r_i, lbl_a in enumerate(b_a.cluster_labels):
                    best_c = int(np.argmax(score_ab[r_i, :]))
                    lbl_b = b_b.cluster_labels[best_c]
                    mutual_r = int(np.argmax(score_ab[:, best_c])) == r_i
                    ov = sorted(
                        top_sets_a[r_i].intersection(top_sets_b[best_c]),
                        key=lambda gi: -(spec_a[r_i, gi] + spec_b[best_c, gi]),
                    )
                    disp_a, branch_a, _ = resolve_atlas_substate_display_metadata(b_a.atlas_id, lbl_a)
                    disp_b, branch_b, _ = resolve_atlas_substate_display_metadata(b_b.atlas_id, lbl_b)
                    rows.append(
                        {
                            "compartment": comp.value,
                            "atlas_a": b_a.atlas_id.value,
                            "cluster_a": lbl_a,
                            "display_a": disp_a,
                            "branch_a": branch_a,
                            "atlas_b": b_b.atlas_id.value,
                            "best_match_cluster_b": lbl_b,
                            "display_b": disp_b,
                            "branch_b": branch_b,
                            "is_reciprocal_best_hit": mutual_r,
                            "spec_cosine": float(cos_ab[r_i, best_c]),
                            "spec_spearman": float(rho_ab[r_i, best_c]),
                            "equivalence_score": float(score_ab[r_i, best_c]),
                            "shared_top_markers": ", ".join(shared[gi] for gi in ov[:8]),
                        }
                    )
    return pl.DataFrame(rows)


# =============================================================================
# 4. Pan-Cohort Bulk ICB Logistic Meta-Analysis & Conditional Stromal Gating
# =============================================================================


def fit_penalized_logistic_newton(
    design_x: FloatMatrix,
    response_y: FloatVector,
    l2_penalty: float = 1e-3,
    max_iter: int = 35,
) -> tuple[FloatVector, FloatVector]:
    """Pure L2-stabilized Newton-Raphson logistic regression returning `(beta, se)`."""
    _, n_feat = design_x.shape
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
    """Exact rank-based Mann-Whitney U / ROC-AUC in [0.0, 1.0]."""
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
    """Fixed-effects inverse-variance meta-analysis returning `(fe_beta, fe_se, two_sided_pval)`."""
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


def dersimonian_laird_meta_analysis(
    betas: FloatVector,
    ses: FloatVector,
) -> tuple[float, float, float, float, float]:
    """DerSimonian-Laird random-effects meta-analysis returning `(dl_beta, dl_se, dl_pval, tau2, i2_pct)`."""
    valid = np.isfinite(betas) & np.isfinite(ses) & (ses > 1e-8)
    k = int(np.sum(valid))
    if k == 0:
        return 0.0, 1.0, 1.0, 0.0, 0.0
    b_v = betas[valid]
    v_v = np.square(ses[valid])
    w = 1.0 / v_v
    w_sum = float(np.sum(w))
    fe_beta = float(np.sum(w * b_v) / w_sum)
    if k == 1:
        se_1 = float(math.sqrt(1.0 / w_sum))
        p_1 = float(2.0 * norm.sf(abs(fe_beta / max(se_1, 1e-12))))
        return fe_beta, se_1, p_1, 0.0, 0.0

    q_stat = float(np.sum(w * np.square(b_v - fe_beta)))
    df_q = float(k - 1)
    c_denom = max(w_sum - float(np.sum(np.square(w))) / w_sum, 1e-12)
    tau2 = max(0.0, (q_stat - df_q) / c_denom)
    i2_pct = max(0.0, ((q_stat - df_q) / max(q_stat, 1e-12)) * 100.0) if q_stat > df_q else 0.0

    w_re = 1.0 / (v_v + tau2)
    w_re_sum = float(np.sum(w_re))
    dl_beta = float(np.sum(w_re * b_v) / w_re_sum)
    dl_se = float(math.sqrt(1.0 / w_re_sum))
    dl_p = float(2.0 * norm.sf(abs(dl_beta / max(dl_se, 1e-12))))
    return dl_beta, dl_se, dl_p, tau2, i2_pct


def _standardize_vec(vec: FloatVector) -> FloatVector:
    std = float(np.std(vec))
    return (
        np.zeros_like(vec, dtype=np.float64)
        if std < 1e-8
        else np.asarray((vec - float(np.mean(vec))) / std, dtype=np.float64)
    )


def _pan_lineage_control_genes(comp: TriadCompartment) -> tuple[str, ...]:
    match comp:
        case TriadCompartment.T_NK:
            return PAN_TCELL_CORE_GENES
        case TriadCompartment.STROMAL_CAF:
            return PAN_CAF_CORE_GENES
        case TriadCompartment.MYELOID_TAM:
            return PAN_MACROPHAGE_CORE_GENES
        case (
            TriadCompartment.B_PLASMA
            | TriadCompartment.ENDOTHELIAL
            | TriadCompartment.EPITHELIAL_TUMOR
            | TriadCompartment.OTHER
        ):
            return PAN_TCELL_CORE_GENES
        case _ as unreachable:
            assert_never(unreachable)


def _classify_ici_direction(beta_val: float, pval: float) -> IciDirection:
    if pval < 0.05 or abs(beta_val) >= 0.14:
        return IciDirection.PRO_RESPONSE if beta_val > 0.0 else IciDirection.PRO_RESISTANCE
    return IciDirection.NEUTRAL


def evaluate_atlas_signatures_on_bulk_icb(
    atlas_markers_df: pl.DataFrame,
    bulk_cohorts: Sequence[BulkIcbCohortBundle],
    n_sig_genes: int = 12,
) -> pl.DataFrame:
    """Evaluate all external atlas substate signatures across the 11 clinical bulk ICB cohorts (`N = 1,058`)."""
    if atlas_markers_df.is_empty():
        return pl.DataFrame()

    substate_keys = (
        atlas_markers_df.select(
            [
                "atlas_id",
                "origin",
                "target_compartment",
                "cluster_label",
                "substate_display",
                "lineage_branch",
                "functional_state",
                "effective_cell_count",
            ]
        )
        .unique(subset=["atlas_id", "cluster_label"], maintain_order=True)
    )

    raw_rows: list[dict[str, object]] = []
    for key_row in substate_keys.iter_rows(named=True):
        aid_str = str(key_row["atlas_id"])
        c_lbl = str(key_row["cluster_label"])
        comp_enum = TriadCompartment(str(key_row["target_compartment"]))
        pan_genes = _pan_lineage_control_genes(comp_enum)

        sig_genes = tuple(
            atlas_markers_df.filter(
                (pl.col("atlas_id") == aid_str) & (pl.col("cluster_label") == c_lbl)
            )
            .sort("marker_rank")
            .head(n_sig_genes)["gene_symbol"]
            .cast(pl.String)
            .to_list()
        )

        raw_betas: list[float] = []
        raw_ses: list[float] = []
        adj_betas: list[float] = []
        adj_ses: list[float] = []
        aurocs: list[float] = []
        weights: list[int] = []
        n_pts = 0

        for cb in bulk_cohorts:
            present_sig = [cb.gene_to_idx[g] for g in sig_genes if g in cb.gene_to_idx]
            present_pan = [cb.gene_to_idx[g] for g in pan_genes if g in cb.gene_to_idx]
            if len(present_sig) < 2 or int(np.sum(cb.response_binary == 1.0)) < 2:
                continue

            s_raw = _standardize_vec(cb.log2_z_mat[:, present_sig].mean(axis=1))
            intercept = np.ones(cb.n_samples, dtype=np.float64)
            b_u, se_u = fit_penalized_logistic_newton(
                np.column_stack([intercept, s_raw]),
                cb.response_binary,
            )
            raw_betas.append(float(b_u[1]))
            raw_ses.append(float(se_u[1]))

            if len(present_pan) >= 2:
                p_vec = _standardize_vec(cb.log2_z_mat[:, present_pan].mean(axis=1))
                b_a, se_a = fit_penalized_logistic_newton(
                    np.column_stack([intercept, s_raw, p_vec]),
                    cb.response_binary,
                )
                adj_betas.append(float(b_a[1]))
                adj_ses.append(float(se_a[1]))

            aurocs.append(compute_mann_whitney_auroc(s_raw, cb.response_binary))
            weights.append(cb.n_samples)
            n_pts += cb.n_samples

        b_raw_arr = np.asarray(raw_betas, dtype=np.float64)
        se_raw_arr = np.asarray(raw_ses, dtype=np.float64)
        b_adj_arr = np.asarray(adj_betas, dtype=np.float64)
        se_adj_arr = np.asarray(adj_ses, dtype=np.float64)

        fe_raw_b, fe_raw_se, fe_raw_p = inverse_variance_meta_analysis(b_raw_arr, se_raw_arr)
        dl_raw_b, dl_raw_se, dl_raw_p, _, raw_i2 = dersimonian_laird_meta_analysis(b_raw_arr, se_raw_arr)
        fe_adj_b, fe_adj_se, fe_adj_p = inverse_variance_meta_analysis(b_adj_arr, se_adj_arr)
        dl_adj_b, dl_adj_se, dl_adj_p, _, adj_i2 = dersimonian_laird_meta_analysis(b_adj_arr, se_adj_arr)
        w_auroc = float(np.average(aurocs, weights=weights)) if weights else 0.5

        raw_rows.append(
            {
                "atlas_id": aid_str,
                "origin": str(key_row["origin"]),
                "target_compartment": comp_enum.value,
                "cluster_label": c_lbl,
                "substate_display": str(key_row["substate_display"]),
                "lineage_branch": str(key_row["lineage_branch"]),
                "functional_state": str(key_row["functional_state"]),
                "effective_cell_count": int(key_row["effective_cell_count"]),
                "signature_genes": ", ".join(sig_genes[:8]),
                "n_bulk_cohorts": len(weights),
                "n_bulk_patients": n_pts,
                "raw_fe_beta": fe_raw_b,
                "raw_fe_se": fe_raw_se,
                "raw_fe_pval": fe_raw_p,
                "raw_dl_beta": dl_raw_b,
                "raw_dl_se": dl_raw_se,
                "raw_dl_pval": dl_raw_p,
                "raw_i2_pct": raw_i2,
                "lineage_adj_fe_beta": fe_adj_b,
                "lineage_adj_fe_se": fe_adj_se,
                "lineage_adj_fe_pval": fe_adj_p,
                "lineage_adj_dl_beta": dl_adj_b,
                "lineage_adj_dl_se": dl_adj_se,
                "lineage_adj_dl_pval": dl_adj_p,
                "lineage_adj_i2_pct": adj_i2,
                "weighted_auroc": w_auroc,
                "ici_direction": _classify_ici_direction(fe_adj_b, fe_adj_p).value,
            }
        )

    out_df = pl.DataFrame(raw_rows)
    if out_df.is_empty():
        return out_df

    # Compute Benjamini-Hochberg FDR within each external atlas
    fdr_raw = np.ones(out_df.height, dtype=np.float64)
    fdr_adj = np.ones(out_df.height, dtype=np.float64)
    atlas_col = out_df["atlas_id"].to_numpy()
    p_raw_all = out_df["raw_fe_pval"].to_numpy().astype(np.float64)
    p_adj_all = out_df["lineage_adj_fe_pval"].to_numpy().astype(np.float64)
    for aid_u in np.unique(atlas_col):
        mask = atlas_col == aid_u
        fdr_raw[mask] = benjamini_hochberg_fdr(p_raw_all[mask])
        fdr_adj[mask] = benjamini_hochberg_fdr(p_adj_all[mask])

    stars_adj = [format_significance_stars(float(q)).value for q in fdr_adj]
    return out_df.with_columns(
        [
            pl.Series("raw_fe_fdr", fdr_raw),
            pl.Series("lineage_adj_fe_fdr", fdr_adj),
            pl.Series("lineage_adj_sig_stars", stars_adj, dtype=pl.String),
        ]
    )


def evaluate_caf_tcell_stromal_gating_on_bulk_icb(
    atlas_markers_df: pl.DataFrame,
    bulk_cohorts: Sequence[BulkIcbCohortBundle],
    n_sig_genes: int = 12,
) -> pl.DataFrame:
    """Test conditional `TGF-beta myCAF` x `Tumor-Reactive CD8+ T` exclusion gating across `N = 1,058` patients."""
    def _lookup_sig(aid: ExternalAtlasId, clust: str, fallback: Sequence[str]) -> tuple[str, ...]:
        if not atlas_markers_df.is_empty():
            sub = (
                atlas_markers_df.filter(
                    (pl.col("atlas_id") == aid.value) & (pl.col("cluster_label") == clust)
                )
                .sort("marker_rank")
                .head(n_sig_genes)
            )
            if sub.height >= 4:
                return tuple(sub["gene_symbol"].cast(pl.String).to_list())
        return tuple(fallback)

    t_tex_genes = _lookup_sig(
        ExternalAtlasId.ZHENG_TCELL,
        "CD8.c12.Tex.CXCL13",
        ("CXCL13", "PDCD1", "HAVCR2", "LAG3", "ENTPD1", "LAYN", "TIGIT", "IFNG", "GZMB"),
    )
    t_tem_genes = _lookup_sig(
        ExternalAtlasId.ZHENG_TCELL,
        "CD8.c05.Tem.GZMK",
        ("GZMK", "CXCR4", "CXCR3", "EOMES", "NKG7", "CCL5", "GZMA", "CD8A"),
    )
    gao_mycaf_genes = _lookup_sig(
        ExternalAtlasId.GAO_CAF,
        "c04",
        ("LRRC15", "HOPX", "MMP11", "POSTN", "COL11A1", "CTHRC1", "COMP", "INHBA", "FAP"),
    )
    cords_mcaf_genes = _lookup_sig(
        ExternalAtlasId.CORDS_CAF,
        "mCAF",
        ("MMP11", "POSTN", "COL11A1", "LRRC15", "CTHRC1", "FAP", "FN1", "VCAN", "TGFB1"),
    )
    buechler_lrrc15_genes = _lookup_sig(
        ExternalAtlasId.BUECHLER_CAF,
        "LRRC15",
        ("LRRC15", "CTHRC1", "COL11A1", "MMP11", "POSTN", "FAP", "COMP", "TAGLN", "ACTA2"),
    )

    pairings: tuple[tuple[str, tuple[str, ...], str, tuple[str, ...]], ...] = (
        ("Zheng CD8.c12.Tex.CXCL13", t_tex_genes, "Gao c04 LRRC15+HOPX+ myCAF", gao_mycaf_genes),
        ("Zheng CD8.c12.Tex.CXCL13", t_tex_genes, "Cords mCAF (MMP11+POSTN+)", cords_mcaf_genes),
        ("Zheng CD8.c12.Tex.CXCL13", t_tex_genes, "Buechler LRRC15 Myofibroblast", buechler_lrrc15_genes),
        ("Zheng CD8.c05.Tem.GZMK", t_tem_genes, "Gao c04 LRRC15+HOPX+ myCAF", gao_mycaf_genes),
    )

    results: list[dict[str, object]] = []
    for t_lbl, t_genes, c_lbl, c_genes in pairings:
        b_t_list: list[float] = []
        se_t_list: list[float] = []
        b_c_list: list[float] = []
        se_c_list: list[float] = []
        b_i_list: list[float] = []
        se_i_list: list[float] = []
        b_t_low_caf: list[float] = []
        se_t_low_caf: list[float] = []
        b_t_high_caf: list[float] = []
        se_t_high_caf: list[float] = []

        y_ll: list[float] = []
        y_hl: list[float] = []
        y_lh: list[float] = []
        y_hh: list[float] = []
        n_used_cohorts = 0
        n_used_pts = 0

        for cb in bulk_cohorts:
            idx_t = [cb.gene_to_idx[g] for g in t_genes if g in cb.gene_to_idx]
            idx_c = [cb.gene_to_idx[g] for g in c_genes if g in cb.gene_to_idx]
            if len(idx_t) < 2 or len(idx_c) < 2 or int(np.sum(cb.response_binary == 1.0)) < 2:
                continue

            z_t = _standardize_vec(cb.log2_z_mat[:, idx_t].mean(axis=1))
            z_c = _standardize_vec(cb.log2_z_mat[:, idx_c].mean(axis=1))
            z_int = z_t * z_c
            ones = np.ones(cb.n_samples, dtype=np.float64)

            b_full, se_full = fit_penalized_logistic_newton(
                np.column_stack([ones, z_t, z_c, z_int]),
                cb.response_binary,
                l2_penalty=5e-3,
            )
            b_t_list.append(float(b_full[1]))
            se_t_list.append(float(se_full[1]))
            b_c_list.append(float(b_full[2]))
            se_c_list.append(float(se_full[2]))
            b_i_list.append(float(b_full[3]))
            se_i_list.append(float(se_full[3]))

            med_t = float(np.median(z_t))
            med_c = float(np.median(z_c))
            low_c_mask = z_c <= med_c
            high_c_mask = ~low_c_mask
            low_t_mask = z_t <= med_t
            high_t_mask = ~low_t_mask

            if int(np.sum(low_c_mask)) >= 6 and int(np.sum(cb.response_binary[low_c_mask] == 1.0)) >= 1:
                b_lc, se_lc = fit_penalized_logistic_newton(
                    np.column_stack([np.ones(int(np.sum(low_c_mask))), z_t[low_c_mask]]),
                    cb.response_binary[low_c_mask],
                    l2_penalty=1e-2,
                )
                b_t_low_caf.append(float(b_lc[1]))
                se_t_low_caf.append(float(se_lc[1]))

            if int(np.sum(high_c_mask)) >= 6 and int(np.sum(cb.response_binary[high_c_mask] == 1.0)) >= 1:
                b_hc, se_hc = fit_penalized_logistic_newton(
                    np.column_stack([np.ones(int(np.sum(high_c_mask))), z_t[high_c_mask]]),
                    cb.response_binary[high_c_mask],
                    l2_penalty=1e-2,
                )
                b_t_high_caf.append(float(b_hc[1]))
                se_t_high_caf.append(float(se_hc[1]))

            y_ll.extend(cb.response_binary[low_t_mask & low_c_mask].tolist())
            y_hl.extend(cb.response_binary[high_t_mask & low_c_mask].tolist())
            y_lh.extend(cb.response_binary[low_t_mask & high_c_mask].tolist())
            y_hh.extend(cb.response_binary[high_t_mask & high_c_mask].tolist())
            n_used_cohorts += 1
            n_used_pts += cb.n_samples

        bt_m, set_m, pt_m = inverse_variance_meta_analysis(
            np.asarray(b_t_list, dtype=np.float64),
            np.asarray(se_t_list, dtype=np.float64),
        )
        bc_m, sec_m, pc_m = inverse_variance_meta_analysis(
            np.asarray(b_c_list, dtype=np.float64),
            np.asarray(se_c_list, dtype=np.float64),
        )
        bi_m, sei_m, pi_m = inverse_variance_meta_analysis(
            np.asarray(b_i_list, dtype=np.float64),
            np.asarray(se_i_list, dtype=np.float64),
        )
        blc_m, _, plc_m = inverse_variance_meta_analysis(
            np.asarray(b_t_low_caf, dtype=np.float64),
            np.asarray(se_t_low_caf, dtype=np.float64),
        )
        bhc_m, _, phc_m = inverse_variance_meta_analysis(
            np.asarray(b_t_high_caf, dtype=np.float64),
            np.asarray(se_t_high_caf, dtype=np.float64),
        )

        gating = StromalGatingResult(
            tcell_axis_label=t_lbl,
            caf_axis_label=c_lbl,
            n_patients=n_used_pts,
            n_cohorts=n_used_cohorts,
            beta_tcell_main=bt_m,
            se_tcell_main=set_m,
            pval_tcell_main=pt_m,
            beta_caf_main=bc_m,
            se_caf_main=sec_m,
            pval_caf_main=pc_m,
            beta_interaction=bi_m,
            se_interaction=sei_m,
            pval_interaction=pi_m,
            or_tcell_low_caf=float(math.exp(np.clip(blc_m, -5.0, 5.0))),
            pval_tcell_low_caf=plc_m,
            or_tcell_high_caf=float(math.exp(np.clip(bhc_m, -5.0, 5.0))),
            pval_tcell_high_caf=phc_m,
            resp_rate_low_t_low_caf=float(np.mean(y_ll)) if y_ll else 0.0,
            n_low_t_low_caf=len(y_ll),
            resp_rate_high_t_low_caf=float(np.mean(y_hl)) if y_hl else 0.0,
            n_high_t_low_caf=len(y_hl),
            resp_rate_low_t_high_caf=float(np.mean(y_lh)) if y_lh else 0.0,
            n_low_t_high_caf=len(y_lh),
            resp_rate_high_t_high_caf=float(np.mean(y_hh)) if y_hh else 0.0,
            n_high_t_high_caf=len(y_hh),
        )
        results.append(gating.model_dump())

    return pl.DataFrame(results)


# =============================================================================
# 5. Nearest-Substate Synthesis, Triad Specificity, Ward Biclustering & Report
# =============================================================================


def _extract_meta42_ici_lookup(
    meta42_bundle: Meta42ReferenceBundle,
) -> tuple[dict[str, tuple[float, float]], dict[str, float]]:
    """Extract `Pan-iAtlas-Meta` bulk logistic `(beta, pval)` and Milo `log2FC` per Meta42 state."""
    bulk_map: dict[str, tuple[float, float]] = {}
    milo_map: dict[str, float] = {}

    se_df = meta42_bundle.state_effects_df
    if not se_df.is_empty():
        sub_se = se_df
        if "reference_id" in sub_se.columns:
            sub_se = sub_se.filter(pl.col("reference_id") == "Consensus-Meta42")
        if "bulk_cohort" in sub_se.columns:
            pan_rows = sub_se.filter(pl.col("bulk_cohort") == "Pan-iAtlas-Meta")
            sub_se = pan_rows if not pan_rows.is_empty() else sub_se
        for r in sub_se.iter_rows(named=True):
            st = str(r.get("cell_state", r.get("cluster", ""))).strip()
            if st and st not in bulk_map:
                raw_b = r.get("bulk_deconv_beta", r.get("meta_beta", 0.0))
                raw_p = r.get("bulk_deconv_pval", r.get("meta_pval", 1.0))
                b_val = float(raw_b) if isinstance(raw_b, (int, float)) else 0.0
                p_val = float(raw_p) if isinstance(raw_p, (int, float)) else 1.0
                bulk_map[st] = (b_val, p_val)
                raw_milo = r.get("milo_mean_logfc")
                if isinstance(raw_milo, (int, float)):
                    milo_map[st] = float(raw_milo)

    hub_df = meta42_bundle.hub_summary_df
    if not hub_df.is_empty():
        for r in hub_df.iter_rows(named=True):
            st = str(r.get("cell_state", r.get("cluster", ""))).strip()
            if st and st not in bulk_map and "bulk_meta_beta" in r:
                hb = r.get("bulk_meta_beta", 0.0)
                hp = r.get("bulk_meta_pval", 1.0)
                bulk_map[st] = (
                    float(hb) if isinstance(hb, (int, float)) else 0.0,
                    float(hp) if isinstance(hp, (int, float)) else 1.0,
                )
            hm = r.get("milo_mean_logfc")
            if st and st not in milo_map and isinstance(hm, (int, float)):
                milo_map[st] = float(hm)

    return bulk_map, milo_map


def assign_one_to_one_matches_df(full_sim_df: pl.DataFrame) -> pl.DataFrame:
    """Add `is_one_to_one_match`: optimal unique Meta42 <-> atlas substate pairing per atlas.

    Solves a linear sum assignment on `composite_substate_score` among the Meta42 states of the
    atlas's cognate compartment, so no Meta42 state or atlas substate is used more than once.
    """
    if full_sim_df.is_empty():
        return full_sim_df
    frames: list[pl.DataFrame] = []
    for key, grp in full_sim_df.group_by(["atlas_id", "atlas_compartment"], maintain_order=True):
        comp = TriadCompartment(str(key[1]))
        cognate = set(_cognate_meta42_ids_for_compartment(comp))
        sub = grp.filter(pl.col("meta_state").is_in(list(cognate)))
        sub = sub if not sub.is_empty() else grp
        metas = [str(x) for x in sub["meta_state"].unique(maintain_order=True).to_list()]
        exts = [str(x) for x in sub["external_cluster"].unique(maintain_order=True).to_list()]
        m_idx = {m: i for i, m in enumerate(metas)}
        e_idx = {e: i for i, e in enumerate(exts)}
        score = np.full((len(metas), len(exts)), -1e9, dtype=np.float64)
        for m, e, v in zip(
            sub["meta_state"].to_list(),
            sub["external_cluster"].to_list(),
            sub["composite_substate_score"].to_list(),
            strict=True,
        ):
            score[m_idx[str(m)], e_idx[str(e)]] = float(v)
        rows, cols = linear_sum_assignment(-score)
        pairs = {(metas[int(r)], exts[int(c)]) for r, c in zip(rows, cols, strict=True)}
        flags = [
            (str(m), str(e)) in pairs
            for m, e in zip(grp["meta_state"].to_list(), grp["external_cluster"].to_list(), strict=True)
        ]
        frames.append(grp.with_columns(pl.Series("is_one_to_one_match", flags, dtype=pl.Boolean)))
    return pl.concat(frames)


def build_nearest_compartment_matches_df(
    full_sim_df: pl.DataFrame,
    atlas_ici_df: pl.DataFrame,
    meta42_bundle: Meta42ReferenceBundle,
    target_compartment: TriadCompartment,
    meta_state_ids: Sequence[str],
    top_n_per_atlas: int = 3,
) -> pl.DataFrame:
    """Build top-N matched external substates per Meta42 state joined with bulk ICB outcomes."""
    if full_sim_df.is_empty():
        return pl.DataFrame()

    bulk_map, milo_map = _extract_meta42_ici_lookup(meta42_bundle)
    sub_sim = full_sim_df.filter(
        (pl.col("atlas_compartment") == target_compartment.value)
        & (pl.col("meta_state").is_in(list(meta_state_ids)))
        & ((pl.col("rank_within_atlas") <= top_n_per_atlas) | pl.col("is_one_to_one_match"))
    ).sort(["meta_state", "atlas_id", "rank_within_atlas"])

    if sub_sim.is_empty():
        return pl.DataFrame()

    ici_lookup: dict[tuple[str, str], dict[str, object]] = {}
    if not atlas_ici_df.is_empty():
        for ir in atlas_ici_df.iter_rows(named=True):
            ici_lookup[(str(ir["atlas_id"]), str(ir["cluster_label"]))] = dict(ir)

    enriched_rows: list[dict[str, object]] = []
    for sr in sub_sim.iter_rows(named=True):
        m_st = str(sr["meta_state"])
        aid = str(sr["atlas_id"])
        e_cl = str(sr["external_cluster"])
        m_beta, m_pval = bulk_map.get(m_st, (0.0, 1.0))
        m_milo = milo_map.get(m_st, 0.0)
        ici_rec = ici_lookup.get((aid, e_cl), {})
        raw_adj_b = ici_rec.get("lineage_adj_fe_beta", 0.0)
        raw_adj_p = ici_rec.get("lineage_adj_fe_pval", 1.0)
        raw_adj_q = ici_rec.get("lineage_adj_fe_fdr", 1.0)
        raw_fe_b = ici_rec.get("raw_fe_beta", 0.0)
        raw_auroc = ici_rec.get("weighted_auroc", 0.5)
        ext_adj_b = float(raw_adj_b) if isinstance(raw_adj_b, (int, float)) else 0.0
        ext_adj_p = float(raw_adj_p) if isinstance(raw_adj_p, (int, float)) else 1.0
        ext_adj_q = float(raw_adj_q) if isinstance(raw_adj_q, (int, float)) else 1.0
        ext_raw_b = float(raw_fe_b) if isinstance(raw_fe_b, (int, float)) else 0.0
        ext_auroc = float(raw_auroc) if isinstance(raw_auroc, (int, float)) else 0.5
        ext_dir = str(ici_rec.get("ici_direction", IciDirection.NEUTRAL.value))

        same_sign = (m_beta * ext_adj_b > 0.0) or (abs(m_beta) < 0.10 and abs(ext_adj_b) < 0.08)
        row_out = dict(sr)
        row_out.update(
            {
                "meta42_bulk_beta": m_beta,
                "meta42_bulk_odds_ratio": float(math.exp(np.clip(m_beta, -5.0, 5.0))),
                "meta42_bulk_pval": m_pval,
                "meta42_milo_logfc": m_milo,
                "ext_raw_fe_beta": ext_raw_b,
                "ext_lineage_adj_beta": ext_adj_b,
                "ext_lineage_adj_pval": ext_adj_p,
                "ext_lineage_adj_fdr": ext_adj_q,
                "ext_weighted_auroc": ext_auroc,
                "ext_ici_direction": ext_dir,
                "is_ici_direction_concordant": same_sign,
            }
        )
        enriched_rows.append(row_out)

    return pl.DataFrame(enriched_rows)


def compute_triad_block_specificity_summary_df(
    full_sim_df: pl.DataFrame,
) -> pl.DataFrame:
    """Summarize block-diagonal Triad specificity across all 42 `Consensus-Meta42` states."""
    if full_sim_df.is_empty():
        return pl.DataFrame()

    meta_states = sorted(full_sim_df["meta_state"].unique().to_list())
    triad_comps = (
        TriadCompartment.T_NK,
        TriadCompartment.STROMAL_CAF,
        TriadCompartment.MYELOID_TAM,
    )
    present_comps = [
        c for c in triad_comps if not full_sim_df.filter(pl.col("atlas_compartment") == c.value).is_empty()
    ]

    rows: list[dict[str, object]] = []
    for m_st in meta_states:
        m_df_raw = full_sim_df.filter(pl.col("meta_state") == m_st)
        m_df_emp = m_df_raw.filter(
            (pl.col("origin") == ProfileOrigin.EMPIRICAL_SINGLE_CELL_H5AD.value)
            & (pl.col("external_cluster") != "20_TDoub")
        )
        m_df = (
            m_df_emp
            if not m_df_emp.is_empty()
            else m_df_raw.filter(pl.col("external_cluster") != "20_TDoub")
        )
        if m_df.is_empty():
            m_df = m_df_raw
        m_comp = classify_meta42_triad_compartment(str(m_st))
        m_lbl = META42_BIO_LABELS.get(str(m_st), str(m_st))

        comp_max_affinity: dict[str, float] = {}
        comp_best_substate: dict[str, str] = {}
        for c_enum in present_comps:
            c_sub = m_df.filter(pl.col("atlas_compartment") == c_enum.value).sort(
                "hybrid_pan_affinity", descending=True
            )
            if not c_sub.is_empty():
                top_r = c_sub.row(0, named=True)
                comp_max_affinity[c_enum.value] = float(top_r["hybrid_pan_affinity"])
                comp_best_substate[c_enum.value] = str(top_r["external_display"])

        best_triad_comp = (
            max(comp_max_affinity.items(), key=lambda kv: kv[1])[0]
            if comp_max_affinity
            else "NONE"
        )
        best_overall_row = m_df.sort("hybrid_pan_affinity", descending=True).row(0, named=True)
        in_triad = m_comp in triad_comps and m_comp.value in comp_max_affinity
        cognate_aff = comp_max_affinity.get(m_comp.value, float("nan"))
        off_vals = [v for k, v in comp_max_affinity.items() if k != m_comp.value]
        max_off_aff = max(off_vals) if off_vals else 0.0
        margin = (cognate_aff - max_off_aff) if in_triad else float(best_overall_row["hybrid_pan_affinity"])

        rows.append(
            {
                "meta_state": str(m_st),
                "meta_state_label": m_lbl,
                "meta_compartment": m_comp.value,
                "is_in_triad_lineage": m_comp in triad_comps,
                "max_affinity_t_nk": comp_max_affinity.get(TriadCompartment.T_NK.value, 0.0),
                "max_affinity_stromal_caf": comp_max_affinity.get(TriadCompartment.STROMAL_CAF.value, 0.0),
                "max_affinity_myeloid_tam": comp_max_affinity.get(TriadCompartment.MYELOID_TAM.value, 0.0),
                "assigned_triad_compartment": best_triad_comp,
                "is_cognate_top1_compartment": (best_triad_comp == m_comp.value) if in_triad else False,
                "cognate_specificity_margin": margin,
                "best_external_substate": str(best_overall_row["external_display"]),
                "best_external_atlas": str(best_overall_row["atlas_id"]),
                "best_hybrid_pan_affinity": float(best_overall_row["hybrid_pan_affinity"]),
                "best_shared_markers": str(best_overall_row["shared_signature_genes"]),
            }
        )
    return pl.DataFrame(rows)


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

    ordered_mat = matrix[np.ix_(row_perm, col_perm)]
    r_centered = np.arange(n_rows, dtype=np.float64) - 0.5 * float(n_rows - 1)
    c_centered = np.arange(n_cols, dtype=np.float64) - 0.5 * float(n_cols - 1)
    diag_cov = float(np.sum(np.outer(r_centered, c_centered) * ordered_mat))
    final_col_perm = list(reversed(col_perm)) if diag_cov < 0.0 else col_perm

    return Success(
        BiclusterResult(
            ordered_rows=tuple(row_labels[i] for i in row_perm),
            ordered_cols=tuple(col_labels[j] for j in final_col_perm),
            row_brackets=_build_brackets_from_linkage(row_z, row_perm),
            col_brackets=_build_brackets_from_linkage(col_z, final_col_perm),
        )
    )


def render_pan_atlas_markdown_report(
    bundles: Sequence[AtlasCentroidBundle],
    nearest_tcell_df: pl.DataFrame,
    nearest_caf_df: pl.DataFrame,
    cross_equiv_df: pl.DataFrame,
    atlas_ici_df: pl.DataFrame,
    gating_df: pl.DataFrame,
    triad_summary_df: pl.DataFrame,
) -> str:
    """Render a publication-grade Markdown synthesis report across all external Triad atlases."""
    lines: list[str] = [
        "# Pan-Cancer External Reference Atlas Validation: T-Cell, CAF & Myeloid Triad",
        "",
        "## 1. Integrated External Reference Atlases",
        "",
        "| Atlas ID | Citation | Compartment | Provenance | Clusters | Total Cells | Tumor Cells |",
        "| :--- | :--- | :---: | :--- | :---: | :---: | :---: |",
    ]
    for b in bundles:
        lines.append(
            f"| `{b.atlas_id.value}` | {atlas_citation_label(b.atlas_id)} | "
            f"`{atlas_target_compartment(b.atlas_id).value}` | `{b.origin.value}` | "
            f"`{len(b.cluster_labels)}` | `{int(b.all_cell_counts.sum()):,}` | `{int(b.tumor_cell_counts.sum()):,}` |"
        )

    if not triad_summary_df.is_empty():
        triad_sub = triad_summary_df.filter(pl.col("is_in_triad_lineage"))
        n_triad = triad_sub.height
        n_correct = int(triad_sub.filter(pl.col("is_cognate_top1_compartment")).height)
        margin_arr = triad_sub["cognate_specificity_margin"].to_numpy().astype(np.float64)
        mean_margin = float(np.nanmean(margin_arr)) if margin_arr.size > 0 else 0.0
        lines.extend(
            [
                "",
                "## 2. Unified 42-State Triad Block-Diagonal Specificity",
                "",
                f"- **Cognate Triad Top-1 Assignment Accuracy**: `{n_correct} / {n_triad}` (`{(100.0 * n_correct / max(n_triad, 1)):.1f}%`)",
                f"- **Mean Cognate vs. Off-Target Hybrid Affinity Margin**: `{mean_margin:+.3f}`",
                "",
            ]
        )

    if not nearest_tcell_df.is_empty():
        lines.extend(
            [
                "## 3. One-to-One Matched External T-Cell Substates for `Consensus-Meta42` T/NK States",
                "",
                "| Meta42 State | Bulk ICB $\\beta$ ($p$) | Atlas | One-to-One Matched Substate | Composite | Shared Signature Markers | Lineage-Adj ICB $\\beta$ ($q$) |",
                "| :--- | :--- | :--- | :--- | :---: | :--- | :--- |",
            ]
        )
        for r in nearest_tcell_df.filter(pl.col("is_one_to_one_match")).iter_rows(named=True):
            lines.append(
                f"| **{r['meta_state_label']}** | `{float(r['meta42_bulk_beta']):+.3f}` (`{float(r['meta42_bulk_pval']):.2e}`) | "
                f"`{r['atlas_id']}` | **{r['external_display']}** | `{float(r['composite_substate_score']):+.3f}` | "
                f"`{r['shared_signature_genes']}` | `{float(r['ext_lineage_adj_beta']):+.3f}` (`{float(r['ext_lineage_adj_fdr']):.2e}`) |"
            )

    if not nearest_caf_df.is_empty():
        lines.extend(
            [
                "",
                "## 4. One-to-One Matched External CAF Substates for `Consensus-Meta42` Stromal & Mesenchymal States",
                "",
                "| Meta42 State | Bulk ICB $\\beta$ ($p$) | Atlas | One-to-One Matched Substate | Composite | Shared Signature Markers | Lineage-Adj ICB $\\beta$ ($q$) |",
                "| :--- | :--- | :--- | :--- | :---: | :--- | :--- |",
            ]
        )
        for r in nearest_caf_df.filter(pl.col("is_one_to_one_match")).iter_rows(named=True):
            lines.append(
                f"| **{r['meta_state_label']}** | `{float(r['meta42_bulk_beta']):+.3f}` (`{float(r['meta42_bulk_pval']):.2e}`) | "
                f"`{r['atlas_id']}` | **{r['external_display']}** | `{float(r['composite_substate_score']):+.3f}` | "
                f"`{r['shared_signature_genes']}` | `{float(r['ext_lineage_adj_beta']):+.3f}` (`{float(r['ext_lineage_adj_fdr']):.2e}`) |"
            )

    if not gating_df.is_empty():
        lines.extend(
            [
                "",
                "## 5. Conditional `TGF-beta myCAF` x `Tumor-Reactive CD8+ T` Stromal Gating (`N = 1,058` Bulk ICB)",
                "",
                "| T-Cell Axis | CAF Axis | $\\beta_{\\text{T}}$ ($p$) | $\\beta_{\\text{CAF}}$ ($p$) | $\\beta_{\\text{T}\\times\\text{CAF}}$ ($p$) | Resp Rate (High T / Low CAF) | Resp Rate (High T / High CAF) | Resp Rate (Low T / High CAF) |",
                "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |",
            ]
        )
        for gr in gating_df.iter_rows(named=True):
            lines.append(
                f"| **{gr['tcell_axis_label']}** | **{gr['caf_axis_label']}** | "
                f"`{float(gr['beta_tcell_main']):+.3f}` (`{float(gr['pval_tcell_main']):.2e}`) | "
                f"`{float(gr['beta_caf_main']):+.3f}` (`{float(gr['pval_caf_main']):.2e}`) | "
                f"`{float(gr['beta_interaction']):+.3f}` (`{float(gr['pval_interaction']):.2e}`) | "
                f"`{100.0 * float(gr['resp_rate_high_t_low_caf']):.1f}%` (`n={int(gr['n_high_t_low_caf'])}`) | "
                f"`{100.0 * float(gr['resp_rate_high_t_high_caf']):.1f}%` (`n={int(gr['n_high_t_high_caf'])}`) | "
                f"`{100.0 * float(gr['resp_rate_low_t_high_caf']):.1f}%` (`n={int(gr['n_low_t_high_caf'])}`) |"
            )

    if not cross_equiv_df.is_empty():
        rbh = cross_equiv_df.filter(pl.col("is_reciprocal_best_hit")).sort("equivalence_score", descending=True).head(20)
        lines.extend(
            [
                "",
                "## 6. Top Reciprocal Best-Hit Cross-Atlas Substate Equivalences",
                "",
                "| Compartment | Atlas A Substate | Atlas B Substate | Cosine | Spearman $\\rho$ | Equivalence | Shared Markers |",
                "| :---: | :--- | :--- | :---: | :---: | :---: | :--- |",
            ]
        )
        for er in rbh.iter_rows(named=True):
            lines.append(
                f"| `{er['compartment']}` | **{er['display_a']}** | **{er['display_b']}** | "
                f"`{float(er['spec_cosine']):+.3f}` | `{float(er['spec_spearman']):+.3f}` | "
                f"`{float(er['equivalence_score']):+.3f}` | `{er['shared_top_markers']}` |"
            )

    if not atlas_ici_df.is_empty():
        top_ici = atlas_ici_df.sort("lineage_adj_fe_pval").head(24)
        lines.extend(
            [
                "",
                "## 7. Strongest Lineage-Adjusted Clinical Bulk ICB Substate Associations (`11 Cohorts, N = 1,058`)",
                "",
                "| Atlas | Substate | Compartment | Raw FE $\\beta$ ($p$) | Lineage-Adj FE $\\beta$ ($p$, FDR $q$) | DerSimonian-Laird $\\beta$ | AUROC | Signature Genes |",
                "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- |",
            ]
        )
        for ir in top_ici.iter_rows(named=True):
            lines.append(
                f"| `{ir['atlas_id']}` | **{ir['substate_display']}** | `{ir['target_compartment']}` | "
                f"`{float(ir['raw_fe_beta']):+.3f}` (`{float(ir['raw_fe_pval']):.2e}`) | "
                f"`{float(ir['lineage_adj_fe_beta']):+.3f}` (`{float(ir['lineage_adj_fe_pval']):.2e}`, `{float(ir['lineage_adj_fe_fdr']):.2e}`) | "
                f"`{float(ir['lineage_adj_dl_beta']):+.3f}` | `{float(ir['weighted_auroc']):.3f}` | `{ir['signature_genes']}` |"
            )

    return "\n".join(lines) + "\n"


def run_pan_atlas_comparison_pure(
    meta42_bundle: Meta42ReferenceBundle,
    atlas_bundles: Sequence[AtlasCentroidBundle],
    bulk_cohorts: Sequence[BulkIcbCohortBundle],
    config: PanAtlasConfig,
) -> Result[PanAtlasOutputs, str]:
    """Pure end-to-end functional pipeline across all external T-Cell, CAF, and Myeloid Triad atlases."""
    if not atlas_bundles:
        return Failure("No external atlas centroid bundles provided")

    filtered_bundles = tuple(
        filter_atlas_min_cells(b, min_cells=config.min_cells_per_substate)
        for b in atlas_bundles
    )
    centroids_df = serialize_atlas_centroids_to_polars(filtered_bundles)
    markers_df = extract_atlas_top_markers_df(
        filtered_bundles,
        top_n=max(config.top_k_markers, config.bulk_signature_genes),
    )

    sim_frames: list[pl.DataFrame] = []
    for bundle in filtered_bundles:
        match compute_meta42_vs_atlas_similarity_df(
            phi_meta_df=meta42_bundle.phi_meta_df,
            bundle=bundle,
            top_k_markers=config.top_k_markers,
            top_degs_df=Some(meta42_bundle.top_degs_df),
        ):
            case Failure(err):
                return Failure(err)
            case Success(sdf):
                sim_frames.append(sdf)

    full_sim_df = assign_one_to_one_matches_df(pl.concat(sim_frames) if sim_frames else pl.DataFrame())
    cross_equiv_df = compute_cross_atlas_equivalence_df(
        filtered_bundles,
        top_k_markers=config.top_k_markers,
    )
    atlas_ici_df = evaluate_atlas_signatures_on_bulk_icb(
        atlas_markers_df=markers_df,
        bulk_cohorts=bulk_cohorts,
        n_sig_genes=config.bulk_signature_genes,
    )
    gating_df = evaluate_caf_tcell_stromal_gating_on_bulk_icb(
        atlas_markers_df=markers_df,
        bulk_cohorts=bulk_cohorts,
        n_sig_genes=config.bulk_signature_genes,
    )

    nearest_tcell_df = build_nearest_compartment_matches_df(
        full_sim_df=full_sim_df,
        atlas_ici_df=atlas_ici_df,
        meta42_bundle=meta42_bundle,
        target_compartment=TriadCompartment.T_NK,
        meta_state_ids=META42_T_NK_IDS,
        top_n_per_atlas=3,
    )
    nearest_caf_df = build_nearest_compartment_matches_df(
        full_sim_df=full_sim_df,
        atlas_ici_df=atlas_ici_df,
        meta42_bundle=meta42_bundle,
        target_compartment=TriadCompartment.STROMAL_CAF,
        meta_state_ids=META42_CAF_AND_MESENCHYMAL_CONTEXT_IDS,
        top_n_per_atlas=3,
    )
    triad_summary_df = compute_triad_block_specificity_summary_df(full_sim_df)
    report_md = render_pan_atlas_markdown_report(
        bundles=filtered_bundles,
        nearest_tcell_df=nearest_tcell_df,
        nearest_caf_df=nearest_caf_df,
        cross_equiv_df=cross_equiv_df,
        atlas_ici_df=atlas_ici_df,
        gating_df=gating_df,
        triad_summary_df=triad_summary_df,
    )

    return Success(
        PanAtlasOutputs(
            atlas_centroids_df=centroids_df,
            atlas_top_markers_df=markers_df,
            meta42_vs_atlases_similarity_df=full_sim_df,
            nearest_tcell_matches_df=nearest_tcell_df,
            nearest_caf_matches_df=nearest_caf_df,
            cross_atlas_equivalence_df=cross_equiv_df,
            atlas_ici_outcomes_df=atlas_ici_df,
            stromal_gating_df=gating_df,
            triad_block_summary_df=triad_summary_df,
            markdown_report=report_md,
        )
    )
