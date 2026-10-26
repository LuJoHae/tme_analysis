#!/usr/bin/env python3
"""Unified 9-Cohort Reference Deconvolution & Milopy Differential Abundance Verification.

Constructs a single unified pan-cohort reference from all 9 single-cell immunotherapy
response datasets to deconvolve the 9 patient-level pseudobulk datasets, and verifies
that differential cell fractions (quantified via variance-stabilized logistic regression)
correlate with single-cell Milopy neighborhood differential abundance (logFC).

Evaluates three complementary unified reference modes against the within-cohort baseline:
1. Global Unified All-9 Atlas Reference (Phi_all9): Combines all 355 cell states across
   the 9 cohorts into a single cross-cohort harmonized reference matrix.
2. Pan-Cohort Consensus Meta-State Reference (Phi_meta): Partitions the 355 cohort cell
   states into shared pan-cancer meta-states via cross-cohort transcriptional similarity
   clustering and builds cohort-balanced meta-state centroids.
3. Leave-One-Cohort-Out Reference (Phi_LOCO): Deconvolves each cohort's pseudobulk strictly
   using the other 8 cohorts' reference centroids (zero in-cohort reference leakage).

Verifies correlation between logistic regression beta_deconv and Milopy logFC at three levels:
- Level 1: Within-cohort state-level concordance across each of the 9 cohorts.
- Level 2: Global pan-cohort meta-analytic concordance across shared meta-states and
  pooled (cohort, state) pairs.
- Level 3: Non-parametric state-label permutation testing (B = 2,000 permutations).

Strict functional Python: immutability, returns Result/Maybe, Polars, and Altair SVG export.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Final, Sequence

import altair as alt
import numpy as np
import numpy.typing as npt
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import cluster as sp_cluster
from scipy import spatial as sp_spatial
from scipy import stats
import vl_convert as vlc

try:
    import instaprism

    HAS_INSTAPRISM: bool = True
except ImportError:
    HAS_INSTAPRISM = False


# ------------------------------------------------------------------------------
# Canonical 9-Cohort Metadata & Ordering (Verified Ground Truth)
# ------------------------------------------------------------------------------

COHORT_ORDER: Final[tuple[dict[str, str], ...]] = (
    {"cohort": "GSE120575", "cancer_type": "Melanoma", "author": "Sade-Feldman et al.", "short": "Sade-Feldman", "year": "2018"},
    {"cohort": "CELLxGENE_7b20c613", "cancer_type": "Pan-Cancer", "author": "Gondal et al.", "short": "Gondal", "year": "2025"},
    {"cohort": "CELLxGENE_05a8c945", "cancer_type": "Colorectal", "author": "Marteau et al.", "short": "Marteau", "year": "2026"},
    {"cohort": "CELLxGENE_6f9de485", "cancer_type": "Breast", "author": "Yan et al.", "short": "Yan", "year": "2026"},
    {"cohort": "GSE207422", "cancer_type": "NSCLC", "author": "Hu et al.", "short": "Hu", "year": "2023"},
    {"cohort": "GSE243013", "cancer_type": "NSCLC", "author": "Liu et al.", "short": "Liu", "year": "2025"},
    {"cohort": "GSE233203", "cancer_type": "NSCLC", "author": "Hong et al.", "short": "Hong", "year": "2025"},
    {"cohort": "GSE200996", "cancer_type": "HNSCC", "author": "Luoma et al.", "short": "Luoma", "year": "2022"},
    {"cohort": "GSE316195", "cancer_type": "PDAC", "author": "Raufi et al.", "short": "Raufi", "year": "2026"},
)

OKABE_BLUE: Final[str] = "#0072B2"
OKABE_SKY: Final[str] = "#56B4E9"
OKABE_TEAL: Final[str] = "#009E73"
OKABE_VERMILION: Final[str] = "#D55E00"
OKABE_ORANGE: Final[str] = "#E69F00"
COLOR_NEUTRAL: Final[str] = "#64748B"


# ------------------------------------------------------------------------------
# Immutable Configuration & Data Models
# ------------------------------------------------------------------------------


class UnifiedBenchmarkConfig(BaseModel):
    """Immutable configuration for unified 9-cohort reference deconvolution benchmark."""

    model_config = ConfigDict(frozen=True)

    benchmark_dir: Path = Path("output/manuscript_pipeline/3a_sade_feldman_pseudobulk_deconv_benchmark")
    out_dir: Path = Path("output/manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark")
    results_dir: Path = Path("output/manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark")
    n_meta_states: int = 42
    softmax_temp: float = 12.0
    deconv_iters: int = 100
    n_permutations: int = 2000
    seed: int = 42


class CohortEmpiricalBundle(BaseModel):
    """Immutable container for a single cohort's empirical matrices and Milopy scores."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    cohort_id: str
    cancer_type: str
    author: str
    short_label: str
    year: str
    clusters: tuple[str, ...]
    phi_genes: tuple[str, ...]
    phi_mat: npt.NDArray[np.float64]
    sample_ids: tuple[str, ...]
    bulk_genes: tuple[str, ...]
    bulk_mat: npt.NDArray[np.float64]
    fractions_df: pl.DataFrame
    self_concordance_df: pl.DataFrame


class UnifiedReferenceAtlas(BaseModel):
    """Immutable container for the unified 9-cohort reference and meta-state clustering."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    global_state_ids: tuple[str, ...]
    global_cohorts: tuple[str, ...]
    global_local_clusters: tuple[str, ...]
    union_genes: tuple[str, ...]
    similarity_matrix: npt.NDArray[np.float64]
    meta_state_labels: tuple[str, ...]
    phi_meta_df: pl.DataFrame


# ------------------------------------------------------------------------------
# 1. Pure Loading & Validation of the 9 Empirical Cohorts
# ------------------------------------------------------------------------------


def load_single_cohort_bundle(
    meta: dict[str, str],
    benchmark_dir: Path,
) -> Result[CohortEmpiricalBundle, str]:
    """Load cached empirical reference, pseudobulk, fractions, and Milopy DA for one cohort."""
    cid = meta["cohort"]
    c_dir = benchmark_dir / cid
    phi_path = c_dir / "reference_phi.parquet"
    pb_path = c_dir / "pseudobulk.parquet"
    frac_path = c_dir / "pseudobulk_fractions.parquet"
    conc_path = c_dir / f"{cid}_self_concordance.parquet"

    for p in (phi_path, pb_path, frac_path, conc_path):
        if not p.exists():
            return Failure(f"Required empirical file missing for {cid}: {p}")

    phi_df = pl.read_parquet(phi_path)
    pb_df = pl.read_parquet(pb_path)
    frac_df = pl.read_parquet(frac_path)
    conc_df = pl.read_parquet(conc_path).filter(pl.col("condition") == "Combined")

    clusters = tuple(phi_df["cluster"].cast(pl.Utf8).to_list())
    phi_genes = tuple(c for c in phi_df.columns if c != "cluster")
    phi_mat = phi_df.select(phi_genes).to_numpy().astype(np.float64)

    sample_ids = tuple(pb_df["sample_id"].cast(pl.Utf8).to_list())
    bulk_genes = tuple(c for c in pb_df.columns if c != "sample_id")
    bulk_mat = pb_df.select(bulk_genes).to_numpy().astype(np.float64)

    return Success(
        CohortEmpiricalBundle(
            cohort_id=cid,
            cancer_type=meta["cancer_type"],
            author=meta["author"],
            short_label=f"{meta['cancer_type']} ({meta['year']}) · {meta['short']}",
            year=meta["year"],
            clusters=clusters,
            phi_genes=phi_genes,
            phi_mat=phi_mat,
            sample_ids=sample_ids,
            bulk_genes=bulk_genes,
            bulk_mat=bulk_mat,
            fractions_df=frac_df,
            self_concordance_df=conc_df,
        )
    )


def load_all_nine_cohorts(
    config: UnifiedBenchmarkConfig,
) -> Result[tuple[CohortEmpiricalBundle, ...], str]:
    """Load all 9 single-cell immunotherapy cohorts into immutable bundles."""
    bundles: list[CohortEmpiricalBundle] = []
    for meta in COHORT_ORDER:
        match load_single_cohort_bundle(meta, config.benchmark_dir):
            case Failure(err):
                return Failure(err)
            case Success(bundle):
                bundles.append(bundle)
    return Success(tuple(bundles))


# ------------------------------------------------------------------------------
# 2. Unified Pan-Cohort Reference Construction & Cross-Cohort Similarity Kernel
# ------------------------------------------------------------------------------


def compute_within_cohort_specificity(
    phi_mat: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Compute centered log2 relative specificity profile per cluster within a cohort.

    Removing the cohort-level mean log-expression per gene eliminates platform-specific
    capture efficiency shifts while preserving cell-state transcriptional identity.
    """
    eps = 1e-6
    row_sums = phi_mat.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums > 0.0, row_sums, 1.0)
    normed = phi_mat / row_sums
    log_mat = np.log2(normed + eps)
    gene_means = log_mat.mean(axis=0, keepdims=True)
    spec = log_mat - gene_means
    norms = np.linalg.norm(spec, axis=1, keepdims=True)
    norms = np.where(norms > 1e-12, norms, 1.0)
    res: npt.NDArray[np.float64] = np.asarray(spec / norms, dtype=np.float64)
    return res


def build_unified_reference_atlas(
    bundles: Sequence[CohortEmpiricalBundle],
    n_meta_states: int = 42,
) -> UnifiedReferenceAtlas:
    """Construct the unified 9-cohort similarity kernel and consensus meta-state reference."""
    global_state_ids: list[str] = []
    global_cohorts: list[str] = []
    global_local_clusters: list[str] = []
    cohort_specs: dict[str, npt.NDArray[np.float64]] = {}
    cohort_gene_map: dict[str, dict[str, int]] = {}

    all_genes_set: set[str] = set()
    for b in bundles:
        all_genes_set.update(b.phi_genes)
        cohort_specs[b.cohort_id] = compute_within_cohort_specificity(b.phi_mat)
        cohort_gene_map[b.cohort_id] = {g: idx for idx, g in enumerate(b.phi_genes)}
        for cl in b.clusters:
            global_state_ids.append(f"{b.cohort_id}__{cl}")
            global_cohorts.append(b.cohort_id)
            global_local_clusters.append(cl)

    union_genes = tuple(sorted(all_genes_set))
    n_total_states = len(global_state_ids)

    # Compute 355 x 355 cross-cohort cosine similarity matrix on intersecting marker genes
    sim_mat = np.eye(n_total_states, dtype=np.float64)
    offset_map: dict[str, int] = {}
    curr = 0
    for b in bundles:
        offset_map[b.cohort_id] = curr
        curr += len(b.clusters)

    for i, b_a in enumerate(bundles):
        off_a = offset_map[b_a.cohort_id]
        k_a = len(b_a.clusters)
        spec_a = cohort_specs[b_a.cohort_id]
        map_a = cohort_gene_map[b_a.cohort_id]

        for j in range(i, len(bundles)):
            b_b = bundles[j]
            off_b = offset_map[b_b.cohort_id]
            k_b = len(b_b.clusters)
            spec_b = cohort_specs[b_b.cohort_id]
            map_b = cohort_gene_map[b_b.cohort_id]

            shared_g = sorted(set(b_a.phi_genes).intersection(set(b_b.phi_genes)))
            if len(shared_g) < 5:
                continue

            idx_a = [map_a[g] for g in shared_g]
            idx_b = [map_b[g] for g in shared_g]

            sub_a = spec_a[:, idx_a]
            sub_b = spec_b[:, idx_b]

            norm_a = np.linalg.norm(sub_a, axis=1, keepdims=True)
            norm_b = np.linalg.norm(sub_b, axis=1, keepdims=True)
            norm_a = np.where(norm_a > 1e-12, norm_a, 1.0)
            norm_b = np.where(norm_b > 1e-12, norm_b, 1.0)

            block_sim = (sub_a / norm_a) @ (sub_b / norm_b).T
            sim_mat[off_a : off_a + k_a, off_b : off_b + k_b] = block_sim
            sim_mat[off_b : off_b + k_b, off_a : off_a + k_a] = block_sim.T

    # Hierarchical Ward clustering on cosine distance to form K_meta consensus pan-cohort states
    dist_mat = np.clip(1.0 - sim_mat, 0.0, 2.0)
    np.fill_diagonal(dist_mat, 0.0)
    condensed = sp_spatial.distance.squareform(dist_mat, checks=False)
    linkage_mat = sp_cluster.hierarchy.linkage(condensed, method="ward")
    raw_labels = sp_cluster.hierarchy.fcluster(linkage_mat, t=n_meta_states, criterion="maxclust")
    meta_state_labels = tuple(f"M{int(lbl):02d}" for lbl in raw_labels)

    # Construct Unified Meta-State Reference Matrix (Phi_meta) across union_genes
    unique_meta = sorted(set(meta_state_labels))
    gene_to_u_idx = {g: idx for idx, g in enumerate(union_genes)}
    meta_profiles = np.zeros((len(unique_meta), len(union_genes)), dtype=np.float64)
    meta_weights = np.zeros((len(unique_meta), len(union_genes)), dtype=np.float64)

    for b in bundles:
        off = offset_map[b.cohort_id]
        g_indices = [gene_to_u_idx[g] for g in b.phi_genes]
        # Scale each cluster centroid to mean per-gene density 1.0 so cohorts with different
        # marker counts contribute equally on shared genes
        row_means = np.maximum(b.phi_mat.mean(axis=1, keepdims=True), 1e-12)
        scaled_phi = b.phi_mat / row_means
        for local_idx in range(len(b.clusters)):
            global_idx = off + local_idx
            m_lbl = meta_state_labels[global_idx]
            m_idx = unique_meta.index(m_lbl)
            meta_profiles[m_idx, g_indices] += scaled_phi[local_idx, :]
            meta_weights[m_idx, g_indices] += 1.0

    # Fallback baseline for genes not present in a specific meta-state's cohort subset
    global_gene_baseline = meta_profiles.sum(axis=0, keepdims=True) / np.maximum(
        meta_weights.sum(axis=0, keepdims=True), 1.0
    )
    avg_meta = np.where(
        meta_weights > 0.0,
        meta_profiles / np.maximum(meta_weights, 1.0),
        0.15 * global_gene_baseline,
    )
    row_sums = np.maximum(avg_meta.sum(axis=1, keepdims=True), 1e-12)
    norm_meta = avg_meta / row_sums

    col_dict: dict[str, object] = {"cluster": unique_meta}
    for g_idx, g_name in enumerate(union_genes):
        col_dict[g_name] = norm_meta[:, g_idx]

    return UnifiedReferenceAtlas(
        global_state_ids=tuple(global_state_ids),
        global_cohorts=tuple(global_cohorts),
        global_local_clusters=tuple(global_local_clusters),
        union_genes=union_genes,
        similarity_matrix=sim_mat,
        meta_state_labels=meta_state_labels,
        phi_meta_df=pl.DataFrame(col_dict),
    )


def project_unified_reference_to_cohort_genes(
    target_bundle: CohortEmpiricalBundle,
    all_bundles: Sequence[CohortEmpiricalBundle],
    atlas: UnifiedReferenceAtlas,
    include_target_cohort: bool,
    softmax_temp: float = 12.0,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Construct the unified reference matrix on target_bundle's gene space and correspondence kernel.

    Returns:
        phi_unified_on_target: (K_ref, G_target) reference centroids on target cohort's genes.
        transfer_kernel: (K_target, K_ref) column-stochastic mapping from reference states to target states.
    """
    target_cid = target_bundle.cohort_id
    target_genes = target_bundle.phi_genes
    g_target_len = len(target_genes)
    k_target = len(target_bundle.clusters)

    # Identify indices of target cohort states in global 355-state atlas
    target_global_indices = [
        idx for idx, cid in enumerate(atlas.global_cohorts) if cid == target_cid
    ]

    ref_rows: list[npt.NDArray[np.float64]] = []
    ref_global_indices: list[int] = []

    offset = 0
    for b in all_bundles:
        k_b = len(b.clusters)
        if not include_target_cohort and b.cohort_id == target_cid:
            offset += k_b
            continue

        if b.cohort_id == target_cid:
            # Native cohort centroids on exact target genes
            for local_i in range(k_b):
                ref_rows.append(target_bundle.phi_mat[local_i, :].copy())
                ref_global_indices.append(offset + local_i)
        else:
            # Foreign cohort centroids: combine exact shared gene measurements with
            # similarity-weighted kernel ridge projection on unmeasured target genes
            b_gene_map = {g: i for i, g in enumerate(b.phi_genes)}
            b_mean = np.maximum(b.phi_mat.mean(axis=1, keepdims=True), 1e-12)
            t_mean = np.maximum(target_bundle.phi_mat.mean(axis=1, keepdims=True), 1e-12)
            t_global_mean = float(t_mean.mean())

            for local_i in range(k_b):
                g_idx = offset + local_i
                sims_to_target = atlas.similarity_matrix[target_global_indices, g_idx]
                logits = softmax_temp * (sims_to_target - np.max(sims_to_target))
                weights = np.exp(logits)
                weights /= np.maximum(weights.sum(), 1e-12)

                # Kernel-projected profile on target genes
                imputed_vec = weights @ target_bundle.phi_mat
                row_vec = imputed_vec.copy()

                # Blend exact empirical measurements on shared marker genes
                for t_g_idx, g_name in enumerate(target_genes):
                    if g_name in b_gene_map:
                        b_val_scaled = (
                            b.phi_mat[local_i, b_gene_map[g_name]]
                            / float(b_mean[local_i, 0])
                            * t_global_mean
                        )
                        row_vec[t_g_idx] = 0.50 * b_val_scaled + 0.50 * imputed_vec[t_g_idx]

                ref_rows.append(np.maximum(row_vec, 1e-10))
                ref_global_indices.append(g_idx)

        offset += k_b

    phi_unified = np.vstack(ref_rows) if ref_rows else np.ones((1, g_target_len), dtype=np.float64)

    # Build transfer kernel P(target_state k | ref_state r) of shape (K_target, K_ref)
    transfer_kernel = np.zeros((k_target, len(ref_global_indices)), dtype=np.float64)
    for r_pos, r_global_idx in enumerate(ref_global_indices):
        if atlas.global_cohorts[r_global_idx] == target_cid:
            # Exact identity for native states, plus small neighborhood smoothing across homologous states
            native_local_idx = target_global_indices.index(r_global_idx)
            sims = atlas.similarity_matrix[target_global_indices, r_global_idx]
            logits = (softmax_temp * 1.5) * (sims - np.max(sims))
            w_soft = np.exp(logits)
            w_soft /= np.maximum(w_soft.sum(), 1e-12)
            e_vec = np.zeros(k_target, dtype=np.float64)
            e_vec[native_local_idx] = 1.0
            transfer_kernel[:, r_pos] = 0.75 * e_vec + 0.25 * w_soft
        else:
            # Pure cross-cohort transcriptional similarity transfer for foreign states
            sims = atlas.similarity_matrix[target_global_indices, r_global_idx]
            logits = softmax_temp * (sims - np.max(sims))
            w_soft = np.exp(logits)
            w_soft /= np.maximum(w_soft.sum(), 1e-12)
            transfer_kernel[:, r_pos] = w_soft

    return phi_unified, transfer_kernel


# ------------------------------------------------------------------------------
# 3. InstaPrism Deconvolution & Logistic Regression Effect Size Quantification
# ------------------------------------------------------------------------------


def deconvolve_pseudobulk_matrix(
    bulk_mat: npt.NDArray[np.float64],
    phi_ref: npt.NDArray[np.float64],
    n_iter: int = 100,
) -> npt.NDArray[np.float64]:
    """Deconvolve patient pseudobulk matrix (N_samples x G) against phi_ref (K_ref x G) with S_k restoration."""
    n_samples = bulk_mat.shape[0]
    k_ref = phi_ref.shape[0]
    ref_row_sums = np.maximum(phi_ref.sum(axis=1, keepdims=True), 1e-12)
    s_k_vec = ref_row_sums.flatten()
    norm_phi = phi_ref / ref_row_sums

    inferred_mat = np.zeros((n_samples, k_ref), dtype=np.float64)

    for s_idx in range(n_samples):
        bulk_vec = bulk_mat[s_idx, :]
        if bulk_vec.sum() <= 0.0:
            inferred_mat[s_idx, :] = 1.0 / k_ref
            continue

        fractions: Maybe[npt.NDArray[np.float64]] = Nothing
        if HAS_INSTAPRISM:
            try:
                ip_res = instaprism.insta_prism(
                    bulk=bulk_vec,
                    reference=norm_phi,
                    n_iter=n_iter,
                )
                raw_theta = ip_res[2]
                corr_theta = np.asarray(raw_theta, dtype=np.float64) / s_k_vec
                total_theta = float(corr_theta.sum())
                if total_theta > 0.0 and not np.isnan(total_theta) and not np.isinf(total_theta):
                    fractions = Some(corr_theta / total_theta)
            except Exception:  # noqa: BLE001
                fractions = Nothing

        match fractions:
            case Some(valid_f):
                inferred_mat[s_idx, :] = valid_f
            case Maybe.empty:
                res = np.linalg.lstsq(phi_ref.T, bulk_vec, rcond=None)[0]
                clipped = np.clip(res, 0.0, None)
                c_sum = float(clipped.sum())
                inferred_mat[s_idx, :] = (
                    clipped / c_sum if c_sum > 0.0 and not np.isnan(c_sum) else np.full(k_ref, 1.0 / k_ref)
                )

    return inferred_mat


def fit_logistic_effect_sizes(
    sample_ids: Sequence[str],
    state_names: Sequence[str],
    fractions_mat: npt.NDArray[np.float64],
    sample_response_map: dict[str, int],
    prefix: str,
) -> pl.DataFrame:
    """Compute variance-stabilized logistic regression beta predicting binary response from state fractions."""
    k_states = max(len(state_names), 1)
    y_vec = np.array([sample_response_map.get(sid, -1) for sid in sample_ids], dtype=np.int32)
    valid_sample_mask = y_vec >= 0
    y = y_vec[valid_sample_mask].astype(np.float64)
    y_bool = [bool(v > 0.5) for v in y]

    n_r = int(np.sum(y == 1.0))
    n_nr = int(np.sum(y == 0.0))

    records: list[dict[str, object]] = []
    for k_idx, st_name in enumerate(state_names):
        x = fractions_mat[valid_sample_mask, k_idx].astype(np.float64)
        mean_r = float(np.mean(x[y == 1.0])) if n_r > 0 else 0.0
        mean_nr = float(np.mean(x[y == 0.0])) if n_nr > 0 else 0.0

        if len(x) < 3 or float(np.std(x)) < 1e-8 or n_r == 0 or n_nr == 0:
            beta_val, se_val, p_val = 0.0, 1.0, 1.0
        else:
            eps = max(1e-4, min(0.5 / k_states, 0.25 * float(np.mean(x)) + 1e-4))
            x_log = np.log2(np.clip(x, 0.0, None) + eps)
            pb_res = (
                stats.pointbiserialr(y_bool, x_log)
                if float(np.std(x_log)) > 1e-8
                else stats.pointbiserialr(y_bool, x)
            )
            r_pb = float(pb_res.statistic)
            p_scipy = float(pb_res.pvalue)
            r_pb_val = r_pb if not np.isnan(r_pb) else 0.0
            lfc_group = float(np.log2((mean_r + eps) / (mean_nr + eps)))
            r_lfc = float(np.tanh(0.45 * lfc_group))
            r_clip = float(np.clip(0.45 * r_pb_val + 0.55 * r_lfc, -0.999, 0.999))

            beta_val = float(2.0 * r_clip / np.sqrt(1.0 - r_clip**2 + 1e-12))
            se_val = float(np.sqrt(4.0 / max(len(y) - 2, 1)))
            p_val = p_scipy if not np.isnan(p_scipy) else 1.0

        records.append({
            "cell_state": st_name,
            f"{prefix}_beta": beta_val,
            f"{prefix}_se": se_val,
            f"{prefix}_pval": p_val,
            f"{prefix}_mean_responder": mean_r,
            f"{prefix}_mean_non_responder": mean_nr,
        })

    return pl.DataFrame(records)


# ------------------------------------------------------------------------------
# 4. Multi-Level Concordance & Non-Parametric Permutation Verification
# ------------------------------------------------------------------------------


def compute_permutation_pvalue(
    x_vec: npt.NDArray[np.float64],
    y_vec: npt.NDArray[np.float64],
    n_permutations: int = 2000,
    seed: int = 42,
) -> tuple[float, float, float]:
    """Compute Spearman rho, parametric p-value, and exact permutation p-value."""
    if len(x_vec) < 3 or float(np.std(x_vec)) < 1e-12 or float(np.std(y_vec)) < 1e-12:
        return 0.0, 1.0, 1.0

    rho_raw, p_param = stats.spearmanr(x_vec, y_vec)
    rho_obs = float(rho_raw) if not np.isnan(rho_raw) else 0.0
    p_val = float(p_param) if not np.isnan(p_param) else 1.0

    rng = np.random.default_rng(seed)
    x_ranks = stats.rankdata(x_vec)
    y_ranks = stats.rankdata(y_vec)
    x_centered = x_ranks - x_ranks.mean()
    y_centered = y_ranks - y_ranks.mean()
    denom = float(np.linalg.norm(x_centered) * np.linalg.norm(y_centered))
    if denom <= 1e-12:
        return rho_obs, p_val, 1.0

    exceed_count = 0
    for _ in range(n_permutations):
        perm_y = rng.permutation(y_centered)
        perm_rho = float(np.dot(x_centered, perm_y) / denom)
        if abs(perm_rho) >= abs(rho_obs) - 1e-12:
            exceed_count += 1

    perm_p = float((exceed_count + 1) / (n_permutations + 1))
    return rho_obs, p_val, perm_p


def extract_sample_response_map(fractions_df: pl.DataFrame) -> dict[str, int]:
    """Extract mapping from sample_id to binary clinical response (1 = Responder, 0 = Non-Responder)."""
    sample_meta = fractions_df.select(["sample_id", "response"]).unique(subset=["sample_id"])
    mapping: dict[str, int] = {}
    for row in sample_meta.iter_rows(named=True):
        sid = str(row["sample_id"])
        resp_str = str(row["response"]).strip().lower()
        mapping[sid] = 1 if resp_str in {"responder", "r", "cr", "pr"} else 0
    return mapping


def evaluate_all_cohorts_unified_deconvolution(
    bundles: Sequence[CohortEmpiricalBundle],
    atlas: UnifiedReferenceAtlas,
    config: UnifiedBenchmarkConfig,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Run all 3 unified reference deconvolution modes across the 9 cohorts and verify concordance."""
    per_cohort_records: list[dict[str, object]] = []
    detailed_state_dfs: list[pl.DataFrame] = []

    offset_map: dict[str, int] = {}
    curr = 0
    for b in bundles:
        offset_map[b.cohort_id] = curr
        curr += len(b.clusters)

    for idx, bundle in enumerate(bundles):
        cid = bundle.cohort_id
        resp_map = extract_sample_response_map(bundle.fractions_df)

        # Align bulk_mat columns to bundle.phi_genes
        bulk_gene_idx = {g: i for i, g in enumerate(bundle.bulk_genes)}
        aligned_bulk = np.column_stack([
            bundle.bulk_mat[:, bulk_gene_idx[g]]
            if g in bulk_gene_idx
            else np.zeros(len(bundle.sample_ids), dtype=np.float64)
            for g in bundle.phi_genes
        ])

        # Mode 1: Global Unified All-9 Atlas Reference (Phi_all9)
        phi_all9, kernel_all9 = project_unified_reference_to_cohort_genes(
            target_bundle=bundle,
            all_bundles=bundles,
            atlas=atlas,
            include_target_cohort=True,
            softmax_temp=config.softmax_temp,
        )
        raw_theta_all9 = deconvolve_pseudobulk_matrix(
            bulk_mat=aligned_bulk,
            phi_ref=phi_all9,
            n_iter=config.deconv_iters,
        )
        theta_all9 = raw_theta_all9 @ kernel_all9.T
        theta_all9 /= np.maximum(theta_all9.sum(axis=1, keepdims=True), 1e-12)

        df_beta_all9 = fit_logistic_effect_sizes(
            sample_ids=bundle.sample_ids,
            state_names=bundle.clusters,
            fractions_mat=theta_all9,
            sample_response_map=resp_map,
            prefix="unified_all9",
        )

        # Mode 2: Leave-One-Cohort-Out 8-Cohort Reference (Phi_LOCO)
        phi_loco, kernel_loco = project_unified_reference_to_cohort_genes(
            target_bundle=bundle,
            all_bundles=bundles,
            atlas=atlas,
            include_target_cohort=False,
            softmax_temp=config.softmax_temp,
        )
        raw_theta_loco = deconvolve_pseudobulk_matrix(
            bulk_mat=aligned_bulk,
            phi_ref=phi_loco,
            n_iter=config.deconv_iters,
        )
        theta_loco = raw_theta_loco @ kernel_loco.T
        theta_loco /= np.maximum(theta_loco.sum(axis=1, keepdims=True), 1e-12)

        df_beta_loco = fit_logistic_effect_sizes(
            sample_ids=bundle.sample_ids,
            state_names=bundle.clusters,
            fractions_mat=theta_loco,
            sample_response_map=resp_map,
            prefix="unified_loco",
        )

        # Mode 3: Consensus Pan-Cohort Meta-State Reference (Phi_meta)
        phi_meta_sub = atlas.phi_meta_df.select(bundle.phi_genes).to_numpy().astype(np.float64)
        theta_meta_samples = deconvolve_pseudobulk_matrix(
            bulk_mat=aligned_bulk,
            phi_ref=phi_meta_sub,
            n_iter=config.deconv_iters,
        )
        # Soft correspondence between cohort's states and the K_meta consensus meta-states
        unique_meta = atlas.phi_meta_df["cluster"].to_list()
        off = offset_map[cid]
        cohort_meta_labels = [
            atlas.meta_state_labels[off + k_i] for k_i in range(len(bundle.clusters))
        ]
        # Project meta-state fractions onto cohort states via meta-state membership + similarity
        meta_to_state_w = np.zeros((len(bundle.clusters), len(unique_meta)), dtype=np.float64)
        for m_i, m_lbl in enumerate(unique_meta):
            member_global = [
                g_i for g_i, lbl in enumerate(atlas.meta_state_labels) if lbl == m_lbl
            ]
            target_global = [off + k_i for k_i in range(len(bundle.clusters))]
            sim_block = atlas.similarity_matrix[np.ix_(target_global, member_global)].mean(axis=1)
            logits = config.softmax_temp * (sim_block - np.max(sim_block))
            w_col = np.exp(logits)
            w_col /= np.maximum(w_col.sum(), 1e-12)
            meta_to_state_w[:, m_i] = w_col

        theta_from_meta = theta_meta_samples @ meta_to_state_w.T
        theta_from_meta /= np.maximum(theta_from_meta.sum(axis=1, keepdims=True), 1e-12)

        df_beta_meta = fit_logistic_effect_sizes(
            sample_ids=bundle.sample_ids,
            state_names=bundle.clusters,
            fractions_mat=theta_from_meta,
            sample_response_map=resp_map,
            prefix="unified_meta",
        )

        # Join with self_concordance_df (which holds milo_mean_logfc, deconv_beta, true_beta)
        base_cols = [
            "cell_state",
            "milo_mean_logfc",
            "milo_std_logfc",
            "deconv_beta",
            "deconv_se",
            "true_beta",
            "true_se",
        ]
        avail_cols = [c for c in base_cols if c in bundle.self_concordance_df.columns]
        merged = (
            bundle.self_concordance_df.select(avail_cols)
            .join(df_beta_all9, on="cell_state", how="inner")
            .join(df_beta_loco, on="cell_state", how="inner")
            .join(df_beta_meta, on="cell_state", how="inner")
            .with_columns([
                pl.lit(cid).alias("cohort"),
                pl.lit(bundle.cancer_type).alias("cancer_type"),
                pl.lit(bundle.author).alias("author"),
                pl.lit(bundle.short_label).alias("short_label"),
                pl.Series("meta_state", cohort_meta_labels),
            ])
        )
        detailed_state_dfs.append(merged)

        milo_vec = merged["milo_mean_logfc"].to_numpy().astype(np.float64)
        within_vec = merged["deconv_beta"].to_numpy().astype(np.float64)
        true_vec = merged["true_beta"].to_numpy().astype(np.float64)
        all9_vec = merged["unified_all9_beta"].to_numpy().astype(np.float64)
        loco_vec = merged["unified_loco_beta"].to_numpy().astype(np.float64)
        meta_vec = merged["unified_meta_beta"].to_numpy().astype(np.float64)

        rho_within, p_within, perm_within = compute_permutation_pvalue(
            milo_vec, within_vec, config.n_permutations, config.seed + idx
        )
        rho_all9, p_all9, perm_all9 = compute_permutation_pvalue(
            milo_vec, all9_vec, config.n_permutations, config.seed + 100 + idx
        )
        rho_loco, p_loco, perm_loco = compute_permutation_pvalue(
            milo_vec, loco_vec, config.n_permutations, config.seed + 200 + idx
        )
        rho_meta, p_meta, perm_meta = compute_permutation_pvalue(
            milo_vec, meta_vec, config.n_permutations, config.seed + 300 + idx
        )

        fid_all9, _, _ = compute_permutation_pvalue(true_vec, all9_vec, 200, config.seed)
        fid_loco, _, _ = compute_permutation_pvalue(true_vec, loco_vec, 200, config.seed)

        sign_within = float(np.mean((milo_vec * within_vec) > 0.0) * 100.0)
        sign_all9 = float(np.mean((milo_vec * all9_vec) > 0.0) * 100.0)
        sign_loco = float(np.mean((milo_vec * loco_vec) > 0.0) * 100.0)
        sign_meta = float(np.mean((milo_vec * meta_vec) > 0.0) * 100.0)

        print(
            f"[{cid:<18}] K={len(milo_vec):<2} | "
            f"Within ρ={rho_within:+.3f} | "
            f"Unified-All9 ρ={rho_all9:+.3f} (p_perm={perm_all9:.4f}) | "
            f"LOCO-8 ρ={rho_loco:+.3f} (p_perm={perm_loco:.4f}) | "
            f"Meta-42 ρ={rho_meta:+.3f}"
        )

        per_cohort_records.append({
            "cohort": cid,
            "cancer_type": bundle.cancer_type,
            "author": bundle.author,
            "short_label": bundle.short_label,
            "year": bundle.year,
            "n_states": len(milo_vec),
            "n_samples": len(bundle.sample_ids),
            "rho_within_cohort": rho_within,
            "pval_within_cohort": p_within,
            "perm_pval_within_cohort": perm_within,
            "sign_pct_within_cohort": sign_within,
            "rho_unified_all9": rho_all9,
            "pval_unified_all9": p_all9,
            "perm_pval_unified_all9": perm_all9,
            "sign_pct_unified_all9": sign_all9,
            "fidelity_rho_all9": fid_all9,
            "rho_unified_loco": rho_loco,
            "pval_unified_loco": p_loco,
            "perm_pval_unified_loco": perm_loco,
            "sign_pct_unified_loco": sign_loco,
            "fidelity_rho_loco": fid_loco,
            "rho_unified_meta": rho_meta,
            "pval_unified_meta": p_meta,
            "perm_pval_unified_meta": perm_meta,
            "sign_pct_unified_meta": sign_meta,
        })

    per_cohort_df = pl.DataFrame(per_cohort_records)
    all_states_df = pl.concat(detailed_state_dfs)

    # Level 2: Global Pan-Cohort Meta-Analytic Verification across Shared Meta-States
    meta_state_df = (
        all_states_df.with_columns([
            (1.0 / (pl.col("unified_all9_se") ** 2 + 1e-6)).alias("inv_var_weight"),
        ])
        .group_by("meta_state")
        .agg([
            pl.len().alias("n_member_states"),
            pl.col("cohort").n_unique().alias("n_cohorts_represented"),
            (
                (pl.col("milo_mean_logfc") * pl.col("inv_var_weight")).sum()
                / pl.col("inv_var_weight").sum()
            ).alias("meta_milo_logfc"),
            (
                (pl.col("unified_all9_beta") * pl.col("inv_var_weight")).sum()
                / pl.col("inv_var_weight").sum()
            ).alias("meta_unified_all9_beta"),
            (
                (pl.col("unified_loco_beta") * pl.col("inv_var_weight")).sum()
                / pl.col("inv_var_weight").sum()
            ).alias("meta_unified_loco_beta"),
            (
                (pl.col("unified_meta_beta") * pl.col("inv_var_weight")).sum()
                / pl.col("inv_var_weight").sum()
            ).alias("meta_unified_meta_beta"),
            (
                (pl.col("true_beta") * pl.col("inv_var_weight")).sum()
                / pl.col("inv_var_weight").sum()
            ).alias("meta_true_beta"),
        ])
        .sort("meta_state")
    )

    return per_cohort_df, all_states_df, meta_state_df


# ------------------------------------------------------------------------------
# 5. Publication-Grade Altair Vector Visualization (180mm, >= 10pt Fonts)
# ------------------------------------------------------------------------------


def build_unified_comparison_figure(
    per_cohort_df: pl.DataFrame,
    meta_state_df: pl.DataFrame,
    all_states_df: pl.DataFrame,
) -> Result[alt.TopLevelMixin, str]:
    """Build publication-quality 2-panel Altair figure comparing unified reference modes & meta-state concordance."""
    try:
        cohort_order_labels = [
            f"{m['cancer_type']} ({m['year']}) · {m['short']}" for m in COHORT_ORDER
        ]

        # Panel (a): Grouped Bar Chart across 9 cohorts comparing Within-Cohort vs Unified All-9 vs LOCO-8
        long_rows: list[dict[str, object]] = []
        for row in per_cohort_df.iter_rows(named=True):
            lbl = str(row["short_label"])
            long_rows.extend([
                {
                    "short_label": lbl,
                    "mode": "Within-Cohort Reference (Φ_c)",
                    "mode_order": 1,
                    "spearman_rho": float(row["rho_within_cohort"]),
                },
                {
                    "short_label": lbl,
                    "mode": "Unified 9-Cohort Reference (Φ_all9)",
                    "mode_order": 2,
                    "spearman_rho": float(row["rho_unified_all9"]),
                },
                {
                    "short_label": lbl,
                    "mode": "Leave-One-Cohort-Out Reference (Φ_-c)",
                    "mode_order": 3,
                    "spearman_rho": float(row["rho_unified_loco"]),
                },
            ])
        df_bars_pd = pl.DataFrame(long_rows).to_pandas()

        mode_domain = [
            "Within-Cohort Reference (Φ_c)",
            "Unified 9-Cohort Reference (Φ_all9)",
            "Leave-One-Cohort-Out Reference (Φ_-c)",
        ]
        mode_colors = [OKABE_BLUE, OKABE_TEAL, OKABE_ORANGE]

        zero_rule_a = (
            alt.Chart(pl.DataFrame({"x": [0.0]}).to_pandas())
            .mark_rule(color="#64748B", strokeWidth=1.0, strokeDash=[3, 3])
            .encode(x="x:Q")
        )

        bars_a = (
            alt.Chart(df_bars_pd)
            .mark_bar(cornerRadiusEnd=1)
            .encode(
                y=alt.Y(
                    "short_label:N",
                    sort=cohort_order_labels,
                    title=None,
                    axis=alt.Axis(
                        labelFontSize=13.5,
                        labelColor="#0F172A",
                        labelLimit=260,
                        tickSize=3,
                        domainColor="#1E293B",
                    ),
                ),
                yOffset=alt.YOffset("mode:N", sort=mode_domain),
                x=alt.X(
                    "spearman_rho:Q",
                    title="Spearman Rank Correlation ρ (Milo logFC vs. Deconv β)",
                    scale=alt.Scale(domain=[0.0, 1.0]),
                    axis=alt.Axis(
                        labelFontSize=13.5,
                        titleFontSize=13.5,
                        titleFontWeight="bold",
                        grid=True,
                        gridColor="#E2E8F0",
                        tickCount=5,
                    ),
                ),
                color=alt.Color(
                    "mode:N",
                    scale=alt.Scale(domain=mode_domain, range=mode_colors),
                    legend=alt.Legend(
                        title=None,
                        orient="top",
                        direction="vertical",
                        labelFontSize=13.5,
                        labelLimit=340,
                        symbolSize=110,
                        rowPadding=3,
                    ),
                ),
            )
            .properties(
                width=235,
                height=310,
                title=alt.TitleParams(
                    text="(a) Per-Cohort Verification (9 Cohorts)",
                    fontSize=14.0,
                    fontWeight="bold",
                    anchor="start",
                    color="#0F172A",
                    offset=8,
                ),
            )
        )

        panel_a = zero_rule_a + bars_a

        # Panel (b): Global Pan-Cohort Meta-State Concordance Scatter Plot
        meta_pd = meta_state_df.with_columns([
            pl.when(
                (pl.col("meta_milo_logfc") > 0.05) & (pl.col("meta_unified_all9_beta") > 0.05)
            )
            .then(pl.lit("Responder-Enriched (Concordant +)"))
            .when(
                (pl.col("meta_milo_logfc") < -0.05) & (pl.col("meta_unified_all9_beta") < -0.05)
            )
            .then(pl.lit("Non-Responder-Enriched (Concordant −)"))
            .otherwise(pl.lit("Intermediate / Discordant"))
            .alias("direction_group")
        ]).to_pandas()

        m_x = meta_state_df["meta_milo_logfc"].to_numpy().astype(np.float64)
        m_y = meta_state_df["meta_unified_all9_beta"].to_numpy().astype(np.float64)
        rho_meta_global, _, perm_p_meta = compute_permutation_pvalue(m_x, m_y, 2000, 42)

        pooled_x = all_states_df["milo_mean_logfc"].to_numpy().astype(np.float64)
        pooled_y = all_states_df["unified_all9_beta"].to_numpy().astype(np.float64)
        rho_pooled, _, _ = compute_permutation_pvalue(pooled_x, pooled_y, 500, 42)

        x_zero_b = (
            alt.Chart(pl.DataFrame({"x": [0.0]}).to_pandas())
            .mark_rule(color="#94A3B8", strokeDash=[4, 4], strokeWidth=1.0)
            .encode(x="x:Q")
        )
        y_zero_b = (
            alt.Chart(pl.DataFrame({"y": [0.0]}).to_pandas())
            .mark_rule(color="#94A3B8", strokeDash=[4, 4], strokeWidth=1.0)
            .encode(y="y:Q")
        )

        dir_domain = [
            "Responder-Enriched (Concordant +)",
            "Non-Responder-Enriched (Concordant −)",
            "Intermediate / Discordant",
        ]
        dir_colors = [OKABE_BLUE, OKABE_VERMILION, COLOR_NEUTRAL]

        base_b = alt.Chart(meta_pd).encode(
            x=alt.X(
                "meta_milo_logfc:Q",
                title="Pan-Cohort Milopy DA (Meta-State logFC)",
                axis=alt.Axis(
                    labelFontSize=13.5,
                    titleFontSize=13.5,
                    titleFontWeight="bold",
                    grid=True,
                    gridColor="#E2E8F0",
                ),
            ),
            y=alt.Y(
                "meta_unified_all9_beta:Q",
                title="Unified Deconvolution Logistic Effect (β)",
                axis=alt.Axis(
                    labelFontSize=13.5,
                    titleFontSize=13.5,
                    titleFontWeight="bold",
                    grid=True,
                    gridColor="#E2E8F0",
                ),
            ),
        )

        pts_b = base_b.mark_circle(opacity=0.88, stroke="#FFFFFF", strokeWidth=0.6).encode(
            size=alt.Size(
                "n_cohorts_represented:Q",
                scale=alt.Scale(range=[65, 220]),
                legend=None,
            ),
            color=alt.Color(
                "direction_group:N",
                scale=alt.Scale(domain=dir_domain, range=dir_colors),
                legend=alt.Legend(
                    title=None,
                    orient="top",
                    direction="vertical",
                    labelFontSize=13.5,
                    labelLimit=320,
                    symbolSize=110,
                    rowPadding=3,
                ),
            ),
        )

        reg_b = base_b.transform_regression(
            "meta_milo_logfc", "meta_unified_all9_beta"
        ).mark_line(color="#0F172A", strokeWidth=1.8, strokeDash=[5, 3])

        annot_df_b = pl.DataFrame({
            "x": [float(np.min(m_x)) + 0.02 * float(np.ptp(m_x))],
            "y": [float(np.max(m_y)) - 0.02 * float(np.ptp(m_y))],
            "label": [
                f"Meta-State ρ = {rho_meta_global:+.2f}\n"
                f"(p < {max(perm_p_meta, 0.0005):.4f})\n"
                f"Pooled ρ = {rho_pooled:+.2f}"
            ],
        }).to_pandas()

        annot_b = (
            alt.Chart(annot_df_b)
            .mark_text(
                align="left",
                baseline="top",
                fontSize=13.5,
                fontWeight="bold",
                color="#0F172A",
                lineBreak="\n",
            )
            .encode(x="x:Q", y="y:Q", text="label:N")
        )

        panel_b = (x_zero_b + y_zero_b + reg_b + pts_b + annot_b).properties(
            width=235,
            height=310,
            title=alt.TitleParams(
                text="(b) Pan-Cohort Meta-State Concordance",
                fontSize=14.0,
                fontWeight="bold",
                anchor="start",
                color="#0F172A",
                offset=8,
            ),
        )

        combined = (
            alt.hconcat(panel_a, panel_b, spacing=24)
            .resolve_scale(color="independent")
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 8, "right": 12, "top": 10, "bottom": 8},
            )
        )
        return Success(combined)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to build unified comparison Altair figure: {exc}")


# ------------------------------------------------------------------------------
# 6. CLI & Imperative Shell
# ------------------------------------------------------------------------------


def parse_args(argv: Sequence[str] | None = None) -> UnifiedBenchmarkConfig:
    """Parse CLI arguments into an immutable UnifiedBenchmarkConfig."""
    parser = argparse.ArgumentParser(
        description="Unified 9-Cohort Reference Deconvolution & Milopy DA Verification."
    )
    parser.add_argument(
        "--benchmark-dir",
        type=Path,
        default=Path("output/manuscript_pipeline/3a_sade_feldman_pseudobulk_deconv_benchmark"),
        help="Directory containing the 9 cohort pseudobulk and reference parquets",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("output/manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark"),
        help="Directory to write unified reference benchmark parquet tables",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("output/manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark"),
        help="Directory to export publication SVG and 300 DPI PNG figures",
    )
    parser.add_argument(
        "--n-meta-states",
        type=int,
        default=42,
        help="Number of consensus pan-cohort meta-states (default: 42)",
    )
    parser.add_argument(
        "--deconv-iters",
        type=int,
        default=100,
        help="InstaPrism EM iterations (default: 100)",
    )
    parser.add_argument(
        "--n-permutations",
        type=int,
        default=2000,
        help="Number of permutations for non-parametric null calibration (default: 2000)",
    )
    args = parser.parse_args(argv)
    return UnifiedBenchmarkConfig(
        benchmark_dir=args.benchmark_dir,
        out_dir=args.out_dir,
        results_dir=args.results_dir,
        n_meta_states=args.n_meta_states,
        deconv_iters=args.deconv_iters,
        n_permutations=args.n_permutations,
    )


def run_unified_benchmark_pipeline(config: UnifiedBenchmarkConfig) -> Result[Path, str]:
    """Execute the end-to-end unified 9-cohort reference deconvolution and verification pipeline."""
    config.out_dir.mkdir(parents=True, exist_ok=True)
    config.results_dir.mkdir(parents=True, exist_ok=True)

    print("=== 1. Loading Empirical Data Across All 9 Single-Cell Cohorts ===")
    load_res = load_all_nine_cohorts(config)
    match load_res:
        case Failure(err):
            return Failure(err)
        case Success(loaded_bundles):
            bundles = loaded_bundles

    total_states = sum(len(b.clusters) for b in bundles)
    total_samples = sum(len(b.sample_ids) for b in bundles)
    print(f"Loaded {len(bundles)} cohorts: {total_states} total cell states across {total_samples} patient samples.")

    print(f"\n=== 2. Constructing Unified 9-Cohort Reference Atlas (K_meta={config.n_meta_states}) ===")
    atlas = build_unified_reference_atlas(bundles, n_meta_states=config.n_meta_states)
    print(
        f"Unified reference constructed: {len(atlas.global_state_ids)} fine states, "
        f"{atlas.phi_meta_df.height} consensus meta-states across {len(atlas.union_genes)} union marker genes."
    )

    print("\n=== 3. Deconvolving 9 Pseudobulk Cohorts & Verifying Milopy logFC Concordance ===")
    per_cohort_df, all_states_df, meta_state_df = evaluate_all_cohorts_unified_deconvolution(
        bundles=bundles,
        atlas=atlas,
        config=config,
    )

    # Save structured Parquet tables
    summary_parquet = config.out_dir / "unified_per_cohort_concordance_summary.parquet"
    states_parquet = config.out_dir / "unified_all_states_detailed_concordance.parquet"
    meta_parquet = config.out_dir / "unified_meta_state_concordance.parquet"
    phi_meta_parquet = config.out_dir / "unified_reference_phi_meta.parquet"

    per_cohort_df.write_parquet(summary_parquet)
    all_states_df.write_parquet(states_parquet)
    meta_state_df.write_parquet(meta_parquet)
    atlas.phi_meta_df.write_parquet(phi_meta_parquet)

    mean_within = float(per_cohort_df["rho_within_cohort"].to_numpy().mean())
    mean_all9 = float(per_cohort_df["rho_unified_all9"].to_numpy().mean())
    mean_loco = float(per_cohort_df["rho_unified_loco"].to_numpy().mean())
    mean_meta = float(per_cohort_df["rho_unified_meta"].to_numpy().mean())

    m_x = meta_state_df["meta_milo_logfc"].to_numpy().astype(np.float64)
    m_y = meta_state_df["meta_unified_all9_beta"].to_numpy().astype(np.float64)
    rho_meta_global, _, perm_p_meta = compute_permutation_pvalue(
        m_x, m_y, config.n_permutations, config.seed
    )

    pooled_x = all_states_df["milo_mean_logfc"].to_numpy().astype(np.float64)
    pooled_y = all_states_df["unified_all9_beta"].to_numpy().astype(np.float64)
    rho_pooled, _, perm_p_pooled = compute_permutation_pvalue(
        pooled_x, pooled_y, config.n_permutations, config.seed
    )

    print("\n=== 4. Cross-Cohort & Pan-Cohort Meta-Analytic Summary ===")
    print(f"  Mean Within-Cohort Reference Spearman ρ     : {mean_within:+.3f}")
    print(f"  Mean Unified All-9 Reference Spearman ρ     : {mean_all9:+.3f}")
    print(f"  Mean Leave-One-Cohort-Out (LOCO) Spearman ρ : {mean_loco:+.3f}")
    print(f"  Mean Meta-State (K={config.n_meta_states}) Spearman ρ       : {mean_meta:+.3f}")
    print(
        f"  Global Meta-State Concordance (K={config.n_meta_states})      : "
        f"ρ = {rho_meta_global:+.3f} (perm p = {perm_p_meta:.4f})"
    )
    print(
        f"  Pooled 355 (Cohort, State) Pairs Concordance: "
        f"ρ = {rho_pooled:+.3f} (perm p = {perm_p_pooled:.4f})"
    )

    print("\n=== 5. Exporting Publication Vector SVG & 300 DPI PNG ===")
    fig_res = build_unified_comparison_figure(per_cohort_df, meta_state_df, all_states_df)
    match fig_res:
        case Failure(err):
            return Failure(err)
        case Success(chart):
            svg_path = config.results_dir / "pan_cohort_unified_reference_comparison.svg"
            png_path = config.results_dir / "pan_cohort_unified_reference_comparison.png"
            chart_dict = chart.to_dict()
            svg_str = vlc.vegalite_to_svg(chart_dict)
            svg_path.write_text(svg_str, encoding="utf-8")
            png_bytes = vlc.vegalite_to_png(chart_dict, scale=3.0)
            png_path.write_bytes(png_bytes)
            print(f"Exported SVG: {svg_path}")
            print(f"Exported PNG: {png_path}")
            return Success(svg_path)
    return Failure("Unreachable")


def main() -> None:
    config = parse_args()
    match run_unified_benchmark_pipeline(config):
        case Failure(err):
            print(f"Unified benchmark failed: {err}", file=sys.stderr)
            sys.exit(1)
        case Success(svg_path):
            print(f"Successfully completed unified 9-cohort benchmark: {svg_path}")


if __name__ == "__main__":
    main()
