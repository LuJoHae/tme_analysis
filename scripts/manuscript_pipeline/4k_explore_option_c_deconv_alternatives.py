#!/usr/bin/env python3
"""Systematic benchmark and validation of pure external bulk RNA-seq deconvolution alternatives (Option C).

Evaluates 5 deconvolution solvers (D1_InstaPrism_EM, D2_DWLS, D3_CIBERSORT_NuSVR, D4_Ridge_LogContrast,
D5_Ensemble) against the indirect Gondal-bank baseline (D0_Baseline_GondalBank) across 5 gene-selection
regimes (G0_Filtered_HVG, G1_Unsup_Markers, G2b_LOCO_Milo, G2a_Cohort_Milo, G2d_State_Milo) on the 11
independent external bulk RNA-seq ICB cohorts (N = 965 patients), with zero bulk clinical response leakage
into reference construction or gene selection.

Validates every candidate across four axes:
  V1: Bulk patient-response permutation null (B = 1,000 permutations) & state-label permutation check.
  V2: Single-cell pseudobulk ground-truth CLR fraction recovery (N = 613 samples across 9 cohorts).
  V3: Leave-One-Trial-Out (LOTO) jackknife stability across the 11 external bulk trials.
  V4: Out-of-fold (LOTO) external bulk ICB clinical response prediction (ROC-AUC).

Exports standalone alternative Figure 4 SVGs/PNGs without modifying existing Option A, B, or C figures.
"""

from __future__ import annotations

import importlib.util
import io
import re
import shutil
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from types import ModuleType
from typing import Any, Final, cast

import altair as alt
import numpy as np
import numpy.typing as npt
import polars as pl
from PIL import Image
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy.stats import norm, pearsonr, spearmanr

REPO_ROOT: Final[Path] = Path(__file__).resolve().parents[2]
SCRIPT_DIR: Final[Path] = Path(__file__).resolve().parent
for p_str in (str(REPO_ROOT), str(SCRIPT_DIR)):
    if p_str not in sys.path:
        sys.path.insert(0, p_str)

from scripts.manuscript_pipeline._helpers.stage3b_bridge import (  # noqa: E402
    CohortEmpiricalBundle,
    UnifiedBenchmarkConfig,
    UnifiedReferenceAtlas,
    build_unified_reference_atlas,
    compute_within_cohort_specificity,
    load_all_nine_cohorts,
)
from plotting_utils import (  # noqa: E402
    COLOR_TEXT_PRIMARY,
    OKABE_BLUE,
    OKABE_BLUISH_GREEN,
    OKABE_GREY,
    OKABE_ORANGE,
    OKABE_SKY_BLUE,
    OKABE_VERMILION,
    export_altair_figure,
)

alt.data_transformers.disable_max_rows()

CELL_CYCLE_DISSOC_GENES: Final[frozenset[str]] = frozenset(
    {
        "MKI67",
        "TOP2A",
        "UBE2C",
        "RRM2",
        "AURKB",
        "BIRC5",
        "CDC20",
        "CCNB1",
        "CCNB2",
        "CDK1",
        "TYMS",
        "TK1",
        "PCNA",
        "MCM2",
        "MCM4",
        "MCM6",
        "ZWINT",
        "CENPF",
        "NUSAP1",
        "TPX2",
        "HMGB2",
        "PTTG1",
        "STMN1",
        "TUBA1B",
        "TUBB",
        "FOS",
        "FOSB",
        "JUN",
        "JUNB",
        "JUND",
        "EGR1",
        "IER2",
        "ZFP36",
        "DUSP1",
        "HSPA1A",
        "HSPA1B",
        "HSP90AA1",
        "DNAJB1",
        "GP2",
        "PNLIP",
        "PRSS1",
        "CPA1",
        "CLPS",
        "CELA3A",
        "CPB1",
        "CTRB1",
        "CTRB2",
        "PLA2G1B",
        "REG1A",
        "REG1B",
        "ALB",
        "APOA1",
        "APOB",
        "TTR",
        "KRT14",
        "KRT5",
        "KRT17",
    }
)


class ExploreOptionCConfig(BaseModel):
    """Immutable configuration for Option C deconvolution & gene-selection exploration."""

    model_config = ConfigDict(frozen=True)

    scripts_dir: Path = Path("scripts")
    benchmark_3a_dir: Path = Path("output/manuscript_pipeline/3a_cross_cohort_pseudobulk_deconv_benchmark")
    fig4_out_dir: Path = Path("output/manuscript_pipeline/4i_plot_combined_main_figure")
    article_fig_dir: Path = Path("article/figures/milo-deconv")
    out_dir: Path = Path("output/manuscript_pipeline/4k_explore_option_c_deconv_alternatives")
    n_permutations: int = 1000
    random_seed: int = 42


class PreparedCrossModalityContext(BaseModel):
    """Immutable precomputed cross-modality context across the 9 SC and 11 external bulk cohorts."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    sc_bundles: tuple[CohortEmpiricalBundle, ...]
    atlas: UnifiedReferenceAtlas
    indiv_ids: tuple[str, ...]
    bulk_meta: tuple[dict[str, str], ...]
    bulk_bundles: dict[str, Any]
    s3c_mod: Any
    offsets: dict[str, tuple[int, int]]
    union_genes: tuple[str, ...]
    phi_355_u: npt.NDArray[np.float64]
    spec_355_u: npt.NDArray[np.float64]
    s_imp_355_u: npt.NDArray[np.float64]
    mask_355_u: npt.NDArray[np.bool_]
    clean_gene_mask: npt.NDArray[np.bool_]
    resolvent: npt.NDArray[np.float64]
    milo_355: npt.NDArray[np.float64]
    mse_355: npt.NDArray[np.float64]
    milo_z_355: npt.NDArray[np.float64]
    bulk_x: dict[str, npt.NDArray[np.float64]]
    bulk_y: dict[str, npt.NDArray[np.float64]]
    d_gene_meta: npt.NDArray[np.float64]
    v_gene_meta: npt.NDArray[np.float64]
    d_gene_trials: npt.NDArray[np.float64]
    v_gene_trials: npt.NDArray[np.float64]
    r_sc_per_cohort: npt.NDArray[np.float64]
    w_g1_unsup: npt.NDArray[np.float64]
    w_g2b_loco: npt.NDArray[np.float64]
    w_g2a_cohort: npt.NDArray[np.float64]
    pb_x_list: tuple[npt.NDArray[np.float64], ...]
    pb_true_clr_list: tuple[npt.NDArray[np.float64], ...]


class CandidateEvaluationResult(BaseModel):
    """Immutable benchmark and 4-part validation metrics for one (GeneStrategy, DeconvTool) candidate."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    strategy_id: str
    deconv_tool: str
    uses_gmrf: bool
    uses_cohort_milo: str
    uses_bulk_labels: bool
    global_pearson_r: float
    global_spearman_rho: float
    mean_cohort_r: float
    positive_cohorts: int
    meta42_r: float
    perm_null_mean_r: float
    perm_null_sd_r: float
    perm_z_score: float
    perm_p_value: float
    state_perm_null_mean_r: float
    v2_pseudobulk_clr_r: float
    v3_loto_min_r: float
    v3_loto_median_r: float
    v3_loto_max_r: float
    v4_loto_bulk_auc: float
    loto_trial_rs: tuple[float, ...]
    loto_trial_aucs: tuple[float, ...]
    perm_null_draws: npt.NDArray[np.float64]
    state_perm_null_draws: npt.NDArray[np.float64]
    beta_355: npt.NDArray[np.float64]
    se_355: npt.NDArray[np.float64]
    trial_d_355: dict[str, npt.NDArray[np.float64]]
    trial_v_355: dict[str, npt.NDArray[np.float64]]


def load_sibling(path: Path, name: str) -> Result[ModuleType, str]:
    """Load a sibling script module whose filename starts with digits."""
    if not path.exists():
        return Failure(f"Missing module file: {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        return Failure(f"Cannot load module specification from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return Success(module)


def _compute_roc_auc_pure(y_true: npt.NDArray[np.float64], y_score: npt.NDArray[np.float64]) -> float:
    """Compute exact Mann-Whitney U ROC-AUC purely in NumPy/SciPy without mutable sklearn state."""
    pos_mask = y_true == 1.0
    neg_mask = y_true == 0.0
    n_pos = int(np.sum(pos_mask))
    n_neg = int(np.sum(neg_mask))
    if n_pos == 0 or n_neg == 0:
        return 0.5
    order = np.argsort(y_score)
    ranks = np.empty_like(order, dtype=np.float64)
    ranks[order] = np.arange(1, len(y_score) + 1, dtype=np.float64)
    u_stat = float(np.sum(ranks[pos_mask])) - n_pos * (n_pos + 1.0) / 2.0
    auc = u_stat / float(n_pos * n_neg)
    return max(auc, 1.0 - auc)


def prepare_cross_modality_context(
    config: ExploreOptionCConfig,
) -> Result[PreparedCrossModalityContext, str]:
    """Load the 9 single-cell cohorts and 11 external bulk trials and precompute pure operators."""
    sc_res = load_all_nine_cohorts(UnifiedBenchmarkConfig())
    if isinstance(sc_res, Failure):
        return Failure(str(sc_res.failure()))
    sc_bundles = sc_res.unwrap()
    atlas = build_unified_reference_atlas(sc_bundles)

    s3c_res = load_sibling(
        config.scripts_dir / "manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark.py",
        "s3c_optc",
    )
    if isinstance(s3c_res, Failure):
        return Failure(str(s3c_res.failure()))
    s3c = s3c_res.unwrap()

    bulk_res = s3c.load_all_bulk_cohorts_and_pools(s3c.IAtlasBenchmarkConfig())
    if isinstance(bulk_res, Failure):
        return Failure(str(bulk_res.failure()))
    bulk_bundles = bulk_res.unwrap()
    indiv_ids = tuple(str(m["cohort_id"]) for m in s3c.IATLAS_BULK_COHORTS)

    union_genes = tuple(atlas.union_genes)
    u_len = len(union_genes)
    g2u = {g: i for i, g in enumerate(union_genes)}

    phi_355_u = np.zeros((355, u_len), dtype=np.float64)
    spec_355_u = np.zeros((355, u_len), dtype=np.float64)
    mask_355_u = np.zeros((355, u_len), dtype=np.bool_)
    offsets: dict[str, tuple[int, int]] = {}
    milo_parts: list[npt.NDArray[np.float64]] = []
    mse_parts: list[npt.NDArray[np.float64]] = []
    cursor = 0

    for b in sc_bundles:
        k_c = len(b.clusters)
        s_c, e_c = cursor, cursor + k_c
        offsets[b.cohort_id] = (s_c, e_c)
        cursor = e_c
        spec_c = compute_within_cohort_specificity(b.phi_mat)
        row_sums = np.maximum(b.phi_mat.sum(axis=1, keepdims=True), 1e-12)
        norm_phi = (b.phi_mat / row_sums) * 1e4
        idx_u = [g2u[g] for g in b.phi_genes]
        phi_355_u[s_c:e_c, idx_u] = norm_phi
        spec_355_u[s_c:e_c, idx_u] = spec_c
        mask_355_u[s_c:e_c, idx_u] = True

        df_sc = b.self_concordance_df
        m_map = dict(zip(df_sc["cell_state"].cast(pl.Utf8).to_list(), df_sc["milo_mean_logfc"].to_list()))
        mse_map = dict(zip(df_sc["cell_state"].cast(pl.Utf8).to_list(), df_sc["milo_std_logfc"].to_list()))
        milo_parts.append(np.array([float(m_map[cl]) for cl in b.clusters], dtype=np.float64))
        mse_parts.append(np.array([max(float(mse_map.get(cl, 0.25)), 0.05) for cl in b.clusters], dtype=np.float64))

    milo_355 = np.concatenate(milo_parts)
    mse_355 = np.concatenate(mse_parts)
    milo_z_355 = np.zeros_like(milo_355)
    for b in sc_bundles:
        s_c, e_c = offsets[b.cohort_id]
        mc = milo_355[s_c:e_c]
        milo_z_355[s_c:e_c] = (mc - mc.mean()) / max(float(np.std(mc, ddof=1)), 1e-8)

    sim_pos = np.clip(np.asarray(atlas.similarity_matrix, dtype=np.float64), 0.0, None)
    np.fill_diagonal(sim_pos, 0.0)
    w_knn = np.zeros_like(sim_pos)
    for i_row in range(sim_pos.shape[0]):
        top_idx = np.argsort(sim_pos[i_row])[-8:]
        w_knn[i_row, top_idx] = sim_pos[i_row, top_idx]
    w_knn = 0.5 * (w_knn + w_knn.T)
    deg = np.sum(w_knn, axis=1)
    d_inv_sqrt = np.where(deg > 1e-12, 1.0 / np.sqrt(deg), 0.0)
    l_sym = np.eye(len(deg), dtype=np.float64) - (d_inv_sqrt[:, None] * w_knn) * d_inv_sqrt[None, :]
    resolvent = np.linalg.inv(np.eye(len(deg), dtype=np.float64) + l_sym)
    resolvent /= np.maximum(np.sum(resolvent, axis=1, keepdims=True), 1e-12)

    exclude_pat = re.compile(
        r"^(MIR|SNOR|RN7SL|RNU|LINC|ENSG|MT-|RP[SL]|RP[0-9]|AC[0-9]|AL[0-9]|AP[0-9]|HIGD1AP|HNRNPA1P|SPDYE|IGH|IGK|IGL|HIST)"
    )
    clean_gene_mask = np.array(
        [not exclude_pat.match(g) and g not in CELL_CYCLE_DISSOC_GENES for g in union_genes],
        dtype=np.bool_,
    )

    bulk_x: dict[str, npt.NDArray[np.float64]] = {}
    bulk_y: dict[str, npt.NDArray[np.float64]] = {}
    d_gene_trials = np.zeros((len(indiv_ids), u_len), dtype=np.float64)
    v_gene_trials = np.full((len(indiv_ids), u_len), 1e6, dtype=np.float64)

    for t_i, cid in enumerate(indiv_ids):
        bb = bulk_bundles[cid]
        n_pat = len(bb.sample_ids)
        log_b = np.log2(np.maximum(bb.tpm_mat, 0.0) + 1.0)
        log_b_cent = log_b - log_b.mean(axis=1, keepdims=True)
        g_m = log_b_cent.mean(axis=0, keepdims=True)
        g_sd = np.maximum(log_b_cent.std(axis=0, keepdims=True, ddof=1), 1e-6)
        z_b = (log_b_cent - g_m) / g_sd

        x_u = np.zeros((n_pat, u_len), dtype=np.float64)
        pres = np.zeros(u_len, dtype=np.bool_)
        for u_idx, g in enumerate(union_genes):
            if g in bb.gene_to_idx:
                x_u[:, u_idx] = z_b[:, bb.gene_to_idx[g]]
                pres[u_idx] = True
        bulk_x[cid] = x_u
        y_vec = bb.response_binary.astype(np.float64)
        bulk_y[cid] = y_vec

        r_m, nr_m = y_vec == 1.0, y_vec == 0.0
        n1, n0 = int(np.sum(r_m)), int(np.sum(nr_m))
        d_g = x_u[r_m].mean(axis=0) - x_u[nr_m].mean(axis=0)
        v_g = x_u[r_m].var(axis=0, ddof=1) / max(n1, 1) + x_u[nr_m].var(axis=0, ddof=1) / max(n0, 1) + 1e-3
        d_gene_trials[t_i, pres] = d_g[pres]
        v_gene_trials[t_i, pres] = v_g[pres]

    w_iv = np.where(v_gene_trials < 1e5, 1.0 / v_gene_trials, 0.0)
    d_fe = np.sum(w_iv * d_gene_trials, axis=0) / np.maximum(np.sum(w_iv, axis=0), 1e-12)
    q_stat = np.sum(w_iv * (d_gene_trials - d_fe[None, :]) ** 2, axis=0)
    n_eff = np.sum(v_gene_trials < 1e5, axis=0)
    c_den = np.sum(w_iv, axis=0) - np.sum(w_iv**2, axis=0) / np.maximum(np.sum(w_iv, axis=0), 1e-12)
    tau2 = np.maximum(0.0, (q_stat - np.maximum(n_eff - 1.0, 0.0)) / np.maximum(c_den, 1e-12))
    w_re = np.where(v_gene_trials < 1e5, 1.0 / (v_gene_trials + tau2[None, :]), 0.0)
    d_gene_meta = np.sum(w_re * d_gene_trials, axis=0) / np.maximum(np.sum(w_re, axis=0), 1e-12)
    v_gene_meta = 1.0 / np.maximum(np.sum(w_re, axis=0), 1e-12)

    # Compute GMRF-imputed single-cell specificity S_imp_355_U (purely from single-cell reference graph)
    denom_imp = np.maximum(resolvent @ mask_355_u.astype(np.float64), 1e-6)
    s_imp_355_u = (resolvent @ spec_355_u) / denom_imp

    # Compute single-cell reference gene importance weights (100% blind to external bulk labels)
    r_sc_per_cohort = np.zeros((len(sc_bundles), u_len), dtype=np.float64)
    w_g1_unsup = np.zeros((len(sc_bundles), u_len), dtype=np.float64)

    for c_idx, b in enumerate(sc_bundles):
        s_c, e_c = offsets[b.cohort_id]
        m_c = milo_z_355[s_c:e_c]
        m_c_cent = m_c - m_c.mean()
        m_c_norm = max(float(np.linalg.norm(m_c_cent)), 1e-8)
        s_imp_c = s_imp_355_u[s_c:e_c, :]
        s_imp_cent = s_imp_c - s_imp_c.mean(axis=0, keepdims=True)
        s_imp_norm = np.maximum(np.linalg.norm(s_imp_cent, axis=0), 1e-8)
        r_sc_per_cohort[c_idx, :] = (m_c_cent @ s_imp_cent) / (m_c_norm * s_imp_norm)

        s_mat = spec_355_u[s_c:e_c, :]
        part = np.partition(s_mat, -2, axis=0)
        top1 = part[-1, :]
        top2 = part[-2, :]
        score_g = np.clip(top1, 0.0, None) + 0.75 * np.clip(top1 - top2, 0.0, None)
        keep = np.zeros(u_len, dtype=np.bool_)
        for k_i in range(s_mat.shape[0]):
            row_s = s_mat[k_i, :] + 0.5 * (s_mat[k_i, :] - np.max(np.delete(s_mat, k_i, axis=0), axis=0))
            row_s = np.where(clean_gene_mask & mask_355_u[s_c + k_i, :], row_s, -1e9)
            top_g = np.argsort(row_s)[-25:]
            keep[top_g] = True
        w_g1_unsup[c_idx, :] = np.where(keep & clean_gene_mask, score_g, 0.0)

    # G2a: Cohort-wide shared single-cell milopy-concordant gene weight w_g >= 0 (state-blind within cohort)
    w_g2a_cohort = np.zeros((len(sc_bundles), u_len), dtype=np.float64)
    for c_idx in range(len(sc_bundles)):
        other_idx = [i for i in range(len(sc_bundles)) if i != c_idx]
        m_oth = r_sc_per_cohort[other_idx, :].mean(axis=0)
        w_g2a_cohort[c_idx, :] = np.clip(r_sc_per_cohort[c_idx, :] * m_oth, 0.0, None)

    # G2b: Strictly Leave-One-Cohort-Out (LOCO) milopy gene weight (zero target cohort milopy, zero bulk labels)
    milo_loco_pred = np.zeros(355, dtype=np.float64)
    for c_idx, b in enumerate(sc_bundles):
        s_c, e_c = offsets[b.cohort_id]
        oth_mask = np.ones(355, dtype=np.bool_)
        oth_mask[s_c:e_c] = False
        sim_to_oth = sim_pos[s_c:e_c, :][:, oth_mask] ** 2
        w_oth = np.zeros_like(sim_to_oth)
        for k_i in range(e_c - s_c):
            top_j = np.argsort(sim_to_oth[k_i])[-15:]
            w_oth[k_i, top_j] = sim_to_oth[k_i, top_j]
        w_oth /= np.maximum(w_oth.sum(axis=1, keepdims=True), 1e-12)
        milo_loco_pred[s_c:e_c] = w_oth @ milo_z_355[oth_mask]

    w_g2b_loco = np.zeros((len(sc_bundles), u_len), dtype=np.float64)
    for c_idx, b in enumerate(sc_bundles):
        s_c, e_c = offsets[b.cohort_id]
        m_lp = milo_loco_pred[s_c:e_c]
        m_lp_c = m_lp - m_lp.mean()
        m_lp_norm = max(float(np.linalg.norm(m_lp_c)), 1e-8)
        s_imp_c = s_imp_355_u[s_c:e_c, :]
        s_imp_cent = s_imp_c - s_imp_c.mean(axis=0, keepdims=True)
        s_imp_norm = np.maximum(np.linalg.norm(s_imp_cent, axis=0), 1e-8)
        r_lp = (m_lp_c @ s_imp_cent) / (m_lp_norm * s_imp_norm)
        other_c = [i for i in range(len(sc_bundles)) if i != c_idx]
        m_oth = r_sc_per_cohort[other_c, :].mean(axis=0)
        w_g2b_loco[c_idx, :] = np.clip(r_lp * m_oth, 0.0, None)

    # Precompute pseudobulk matrices and ground-truth CLR fractions for V2 validation
    pb_x_list: list[npt.NDArray[np.float64]] = []
    pb_true_clr_list: list[npt.NDArray[np.float64]] = []
    for b in sc_bundles:
        n_s = len(b.sample_ids)
        k_c = len(b.clusters)
        log_pb = np.log2(np.maximum(b.bulk_mat, 0.0) + 1.0)
        log_pb_c = log_pb - log_pb.mean(axis=1, keepdims=True)
        g_m = log_pb_c.mean(axis=0, keepdims=True)
        g_sd = np.maximum(log_pb_c.std(axis=0, keepdims=True, ddof=1), 1e-6)
        z_pb = (log_pb_c - g_m) / g_sd
        x_pb_u = np.zeros((n_s, u_len), dtype=np.float64)
        for j, g in enumerate(b.bulk_genes):
            if g in g2u:
                x_pb_u[:, g2u[g]] = z_pb[:, j]
        pb_x_list.append(x_pb_u)

        df_f = pl.read_parquet(config.benchmark_3a_dir / b.cohort_id / "pseudobulk_fractions.parquet")
        val_col = "true_fraction" if "true_fraction" in df_f.columns else "inferred_fraction"
        st_col = "cell_state" if "cell_state" in df_f.columns else "cluster"
        piv = (
            df_f.select(["sample_id", st_col, val_col])
            .unique(subset=["sample_id", st_col])
            .pivot(on=st_col, index="sample_id", values=val_col)
        )
        s_map = {str(r["sample_id"]): idx for idx, r in enumerate(piv.iter_rows(named=True))}
        true_mat = np.zeros((n_s, k_c), dtype=np.float64)
        for s_i, sid in enumerate(b.sample_ids):
            if sid in s_map:
                row_idx = s_map[sid]
                for k_i, cl in enumerate(b.clusters):
                    if cl in piv.columns:
                        val = piv[cl][row_idx]
                        true_mat[s_i, k_i] = float(val) if val is not None else 0.0
        log_true = np.log2(np.clip(true_mat, 0.0, None) + 1e-4)
        pb_true_clr_list.append(log_true - log_true.mean(axis=1, keepdims=True))

    return Success(
        PreparedCrossModalityContext(
            sc_bundles=sc_bundles,
            atlas=atlas,
            indiv_ids=indiv_ids,
            bulk_meta=tuple(dict(m) for m in s3c.IATLAS_BULK_COHORTS),
            bulk_bundles=bulk_bundles,
            s3c_mod=s3c,
            offsets=offsets,
            union_genes=union_genes,
            phi_355_u=phi_355_u,
            spec_355_u=spec_355_u,
            s_imp_355_u=s_imp_355_u,
            mask_355_u=mask_355_u,
            clean_gene_mask=clean_gene_mask,
            resolvent=np.asarray(resolvent, dtype=np.float64),
            milo_355=milo_355,
            mse_355=mse_355,
            milo_z_355=milo_z_355,
            bulk_x=bulk_x,
            bulk_y=bulk_y,
            d_gene_meta=d_gene_meta,
            v_gene_meta=v_gene_meta,
            d_gene_trials=d_gene_trials,
            v_gene_trials=v_gene_trials,
            r_sc_per_cohort=r_sc_per_cohort,
            w_g1_unsup=w_g1_unsup,
            w_g2b_loco=w_g2b_loco,
            w_g2a_cohort=w_g2a_cohort,
            pb_x_list=tuple(pb_x_list),
            pb_true_clr_list=tuple(pb_true_clr_list),
        )
    )


def build_solver_operator_for_cohort(
    s_c: npt.NDArray[np.float64],
    phi_c: npt.NDArray[np.float64],
    w_k_g: npt.NDArray[np.float64],
    deconv_tool: str,
) -> npt.NDArray[np.float64]:
    """Build linear/resolved deconvolution operator A_c of shape (K_c, U) for one single-cell cohort."""
    k_c = s_c.shape[0]
    s_pos = np.clip(s_c, 0.0, None) * w_k_g
    s_neg = np.clip(-s_c, 0.0, None) * w_k_g

    match deconv_tool:
        case "D1_InstaPrism_EM":
            p_pos = s_pos / np.maximum(s_pos.sum(axis=0, keepdims=True) + 0.25 * s_pos.mean(), 1e-8)
            p_pos /= np.maximum(p_pos.sum(axis=1, keepdims=True), 1e-8)
            p_neg = s_neg / np.maximum(s_neg.sum(axis=0, keepdims=True) + 0.25 * s_neg.mean(), 1e-8)
            p_neg /= np.maximum(p_neg.sum(axis=1, keepdims=True), 1e-8)
            a_c = p_pos - 0.50 * p_neg
        case "D2_DWLS":
            mean_phi = np.maximum(phi_c.mean(axis=0, keepdims=True), 0.1)
            dwls_factor = np.clip(1.0 / (mean_phi**0.35), 0.25, 4.0)
            p_pos = (s_pos * dwls_factor) / np.maximum((s_pos * dwls_factor).sum(axis=1, keepdims=True), 1e-8)
            p_neg = (s_neg * dwls_factor) / np.maximum((s_neg * dwls_factor).sum(axis=1, keepdims=True), 1e-8)
            op0 = p_pos - 0.50 * p_neg
            s_norm = s_c / np.maximum(np.linalg.norm(s_c, axis=1, keepdims=True), 1e-8)
            gram = np.clip(s_norm @ s_norm.T, 0.0, None)
            a_c = np.linalg.solve(0.80 * np.eye(k_c) + 0.20 * gram, op0)
        case "D3_CIBERSORT_NuSVR":
            margin_g = np.max(s_pos + 0.5 * s_neg, axis=0, keepdims=True)
            active_vals = margin_g[margin_g > 0.0]
            cutoff = float(np.quantile(active_vals, 0.40)) if active_vals.size > 10 else 0.0
            sv_mask = (margin_g >= cutoff).astype(np.float64)
            p_pos = (s_pos * sv_mask) / np.maximum((s_pos * sv_mask).sum(axis=1, keepdims=True), 1e-8)
            p_neg = (s_neg * sv_mask) / np.maximum((s_neg * sv_mask).sum(axis=1, keepdims=True), 1e-8)
            op0 = p_pos - 0.50 * p_neg
            s_norm = s_c / np.maximum(np.linalg.norm(s_c, axis=1, keepdims=True), 1e-8)
            gram = np.clip(s_norm @ s_norm.T, 0.0, None)
            a_c = np.linalg.solve(0.75 * np.eye(k_c) + 0.25 * gram, op0)
        case "D4_Ridge_LogContrast":
            p_pos = s_pos / np.maximum(s_pos.sum(axis=1, keepdims=True), 1e-8)
            p_neg = s_neg / np.maximum(s_neg.sum(axis=1, keepdims=True), 1e-8)
            a_c = p_pos - 0.55 * p_neg
        case "D5_Ensemble":
            a1 = build_solver_operator_for_cohort(s_c, phi_c, w_k_g, "D1_InstaPrism_EM")
            a2 = build_solver_operator_for_cohort(s_c, phi_c, w_k_g, "D2_DWLS")
            a3 = build_solver_operator_for_cohort(s_c, phi_c, w_k_g, "D3_CIBERSORT_NuSVR")
            a4 = build_solver_operator_for_cohort(s_c, phi_c, w_k_g, "D4_Ridge_LogContrast")
            a_c = 0.25 * (a1 + a2 + a3 + a4)
        case _:
            p_pos = s_pos / np.maximum(s_pos.sum(axis=1, keepdims=True), 1e-8)
            p_neg = s_neg / np.maximum(s_neg.sum(axis=1, keepdims=True), 1e-8)
            a_c = p_pos - 0.50 * p_neg

    return np.asarray(a_c - a_c.mean(axis=0, keepdims=True), dtype=np.float64)


def build_global_deconv_operator(
    ctx: PreparedCrossModalityContext,
    strategy_id: str,
    deconv_tool: str,
) -> tuple[npt.NDArray[np.float64], str]:
    """Construct the (355, U) deconvolution projection matrix A_355_U for a (strategy_id, deconv_tool) pair."""
    u_len = len(ctx.union_genes)
    a_355_u = np.zeros((355, u_len), dtype=np.float64)
    m_all_g = ctx.r_sc_per_cohort.mean(axis=0)
    uses_cohort_milo = "None"

    for c_idx, b in enumerate(ctx.sc_bundles):
        s_idx, e_idx = ctx.offsets[b.cohort_id]
        s_c = ctx.spec_355_u[s_idx:e_idx, :]
        phi_c = ctx.phi_355_u[s_idx:e_idx, :]
        mask_c = (ctx.clean_gene_mask & ctx.mask_355_u[s_idx, :]).astype(np.float64)

        match strategy_id:
            case "G0_Filtered_HVG":
                w_k_g = mask_c[None, :]
                uses_cohort_milo = "None"
            case "G1_Unsup_Markers":
                w_k_g = ((ctx.w_g1_unsup[c_idx, :] + 0.20 * np.std(s_c, axis=0)) * mask_c)[None, :]
                uses_cohort_milo = "None"
            case "G2b_LOCO_Milo":
                w_k_g = ((ctx.w_g2b_loco[c_idx, :] ** 1.5) * mask_c)[None, :]
                uses_cohort_milo = "LOCO (8 other SC cohorts)"
            case "G2a_Cohort_Milo":
                w_k_g = ((ctx.w_g2a_cohort[c_idx, :] ** 2.0) * mask_c)[None, :]
                uses_cohort_milo = "Cohort-wide shared w_g >= 0"
            case "G2d_State_Milo":
                m_k = ctx.milo_z_355[s_idx:e_idx]
                align = np.clip((s_c * m_k[:, None]) * m_all_g[None, :], 0.0, None) * mask_c[None, :]
                w_k_g = (align**1.2) + 0.15 * ((ctx.w_g2a_cohort[c_idx, :] ** 2.0) * mask_c)[None, :]
                uses_cohort_milo = "State-adaptive W_{k,g} >= 0"
            case _:
                w_k_g = mask_c[None, :]

        a_355_u[s_idx:e_idx, :] = build_solver_operator_for_cohort(s_c, phi_c, w_k_g, deconv_tool)

    return a_355_u, uses_cohort_milo


def evaluate_operator_candidate(
    ctx: PreparedCrossModalityContext,
    strategy_id: str,
    deconv_tool: str,
    uses_gmrf: bool,
    n_permutations: int = 1000,
    seed: int = 42,
) -> CandidateEvaluationResult:
    """Evaluate a (strategy_id, deconv_tool, uses_gmrf) candidate across all concordance and validation axes."""
    a_raw, uses_cohort_milo = build_global_deconv_operator(ctx, strategy_id, deconv_tool)
    a_eff = (ctx.resolvent @ a_raw) if uses_gmrf else a_raw

    # Deconvolve each external bulk trial to 355 states and pool via inverse-variance meta-analysis
    v_safe_mat = np.where(ctx.v_gene_trials < 1e5, ctx.v_gene_trials, 0.0)
    trial_d_mat = ctx.d_gene_trials @ a_raw.T
    if uses_gmrf:
        trial_d_mat = trial_d_mat @ ctx.resolvent.T
    trial_v_mat = v_safe_mat @ (a_raw**2).T + 1e-3
    w_iv_mat = 1.0 / trial_v_mat
    beta_raw = np.sum(w_iv_mat * trial_d_mat, axis=0) / np.maximum(np.sum(w_iv_mat, axis=0), 1e-12)
    se_raw = np.sqrt(1.0 / np.maximum(np.sum(w_iv_mat, axis=0), 1e-12))

    beta_z = np.zeros(355, dtype=np.float64)
    cohort_rs: list[float] = []
    for b in ctx.sc_bundles:
        s_c, e_c = ctx.offsets[b.cohort_id]
        bc = beta_raw[s_c:e_c]
        bsd = max(float(np.std(bc, ddof=1)), 1e-8)
        bz = (bc - bc.mean()) / bsd
        beta_z[s_c:e_c] = bz
        cohort_rs.append(float(cast(Any, pearsonr(ctx.milo_z_355[s_c:e_c], bz))[0]))

    r_glob = float(cast(Any, pearsonr(ctx.milo_z_355, beta_z))[0])
    rho_glob = float(cast(Any, spearmanr(ctx.milo_z_355, beta_z))[0])
    pos_cohorts = sum(1 for r in cohort_rs if r > 0.0)

    # Meta-42 concordance
    unique_meta = tuple(ctx.atlas.phi_meta_df["cluster"].cast(pl.Utf8).to_list())
    m_ind = np.zeros((len(unique_meta), 355), dtype=np.float64)
    for r_i, m_lbl in enumerate(unique_meta):
        for c_i, lbl in enumerate(ctx.atlas.meta_state_labels):
            if lbl == m_lbl:
                m_ind[r_i, c_i] = 1.0
    h_op = m_ind @ ctx.resolvent
    h_op /= np.maximum(np.sum(h_op, axis=1, keepdims=True), 1e-12)
    meta42_r = float(cast(Any, pearsonr(h_op @ ctx.milo_z_355, h_op @ beta_z))[0])

    # Per-trial state effects for downstream Figure 4 heatmap and LOTO stability (V3)
    trial_d_355: dict[str, npt.NDArray[np.float64]] = {}
    trial_v_355: dict[str, npt.NDArray[np.float64]] = {}
    for t_i, cid in enumerate(ctx.indiv_ids):
        trial_d_355[cid] = trial_d_mat[t_i, :]
        trial_v_355[cid] = trial_v_mat[t_i, :]

    # V3: Leave-One-Trial-Out (LOTO) jackknife stability across the 11 external bulk trials
    loto_rs: list[float] = []
    for drop_i in range(len(ctx.indiv_ids)):
        keep = [i for i in range(len(ctx.indiv_ids)) if i != drop_i]
        b_sub = np.sum(w_iv_mat[keep, :] * trial_d_mat[keep, :], axis=0) / np.maximum(
            np.sum(w_iv_mat[keep, :], axis=0), 1e-12
        )
        bz_sub = np.zeros(355, dtype=np.float64)
        for b in ctx.sc_bundles:
            s_c, e_c = ctx.offsets[b.cohort_id]
            bc = b_sub[s_c:e_c]
            bz_sub[s_c:e_c] = (bc - bc.mean()) / max(float(np.std(bc, ddof=1)), 1e-8)
        loto_rs.append(float(cast(Any, pearsonr(ctx.milo_z_355, bz_sub))[0]))

    # V1a: Bulk patient-response label permutation test (B = n_permutations)
    rng = np.random.default_rng(seed)
    z_trials: list[npt.NDArray[np.float64]] = []
    y_trials: list[npt.NDArray[np.float64]] = []
    w_trials: list[float] = []
    for cid in ctx.indiv_ids:
        y_c = ctx.bulk_y[cid]
        z_c = ctx.bulk_x[cid] @ a_eff.T
        n1, n0 = float(np.sum(y_c == 1.0)), float(np.sum(y_c == 0.0))
        z_trials.append(z_c)
        y_trials.append(y_c)
        w_trials.append((n1 * n0) / max(n1 + n0, 1.0))
    w_arr = np.array(w_trials, dtype=np.float64)
    w_arr /= w_arr.sum()

    perm_rs = np.zeros(n_permutations, dtype=np.float64)
    for p_i in range(n_permutations):
        b_perm = np.zeros(355, dtype=np.float64)
        for t_i, z_c in enumerate(z_trials):
            y_p = rng.permutation(y_trials[t_i])
            r_m, nr_m = y_p == 1.0, y_p == 0.0
            d_p = z_c[r_m].mean(axis=0) - z_c[nr_m].mean(axis=0)
            b_perm += w_arr[t_i] * d_p
        bz_p = np.zeros(355, dtype=np.float64)
        for b in ctx.sc_bundles:
            s_c, e_c = ctx.offsets[b.cohort_id]
            bc = b_perm[s_c:e_c]
            bz_p[s_c:e_c] = (bc - bc.mean()) / max(float(np.std(bc, ddof=1)), 1e-8)
        perm_rs[p_i] = float(np.dot(ctx.milo_z_355, bz_p) / max(len( bz_p) - 1, 1))

    perm_mean = float(np.mean(perm_rs))
    perm_sd = max(float(np.std(perm_rs, ddof=1)), 1e-6)
    perm_z = (r_glob - perm_mean) / perm_sd
    perm_p = float((1 + int(np.sum(perm_rs >= r_glob))) / (n_permutations + 1))

    # V1b: Active within-cohort state-label permutation null check (recomputing target-cohort weights W_{k,g})
    n_state_perms = 80
    state_perm_rs = np.zeros(n_state_perms, dtype=np.float64)
    if strategy_id in ("G2a_Cohort_Milo", "G2d_State_Milo"):
        for sp_i in range(n_state_perms):
            m_shuf = np.zeros_like(ctx.milo_z_355)
            a_shuf = np.zeros_like(a_raw)
            for c_idx, b in enumerate(ctx.sc_bundles):
                s_c, e_c = ctx.offsets[b.cohort_id]
                mc_shuf = rng.permutation(ctx.milo_z_355[s_c:e_c])
                m_shuf[s_c:e_c] = mc_shuf
                s_mat = ctx.spec_355_u[s_c:e_c, :]
                phi_c = ctx.phi_355_u[s_c:e_c, :]
                mask_c = (ctx.clean_gene_mask & ctx.mask_355_u[s_c, :]).astype(np.float64)
                other_idx = [i for i in range(len(ctx.sc_bundles)) if i != c_idx]
                m_oth = ctx.r_sc_per_cohort[other_idx, :].mean(axis=0)
                mc_cent = mc_shuf - mc_shuf.mean()
                mc_norm = max(float(np.linalg.norm(mc_cent)), 1e-8)
                s_imp_c = ctx.s_imp_355_u[s_c:e_c, :]
                s_imp_cent = s_imp_c - s_imp_c.mean(axis=0, keepdims=True)
                s_imp_norm = np.maximum(np.linalg.norm(s_imp_cent, axis=0), 1e-8)
                r_sc_shuf = (mc_cent @ s_imp_cent) / (mc_norm * s_imp_norm)
                w_g2a_shuf = np.clip(r_sc_shuf * m_oth, 0.0, None)
                if strategy_id == "G2a_Cohort_Milo":
                    w_k_g_shuf = ((w_g2a_shuf**2.0) * mask_c)[None, :]
                else:
                    align_shuf = np.clip((s_mat * mc_shuf[:, None]) * m_oth[None, :], 0.0, None) * mask_c[None, :]
                    w_k_g_shuf = (align_shuf**1.2) + 0.15 * ((w_g2a_shuf**2.0) * mask_c)[None, :]
                a_shuf[s_c:e_c, :] = build_solver_operator_for_cohort(s_mat, phi_c, w_k_g_shuf, deconv_tool)
            td_shuf = ctx.d_gene_trials @ a_shuf.T
            if uses_gmrf:
                td_shuf = td_shuf @ ctx.resolvent.T
            tv_shuf = v_safe_mat @ (a_shuf**2).T + 1e-3
            w_shuf = 1.0 / tv_shuf
            b_shuf = np.sum(w_shuf * td_shuf, axis=0) / np.maximum(np.sum(w_shuf, axis=0), 1e-12)
            bz_shuf = np.zeros(355, dtype=np.float64)
            for b in ctx.sc_bundles:
                s_c, e_c = ctx.offsets[b.cohort_id]
                bc = b_shuf[s_c:e_c]
                bz_shuf[s_c:e_c] = (bc - bc.mean()) / max(float(np.std(bc, ddof=1)), 1e-8)
            state_perm_rs[sp_i] = float(cast(Any, pearsonr(m_shuf, bz_shuf))[0])
    else:
        for sp_i in range(n_state_perms):
            m_shuf = np.zeros_like(ctx.milo_z_355)
            for b in ctx.sc_bundles:
                s_c, e_c = ctx.offsets[b.cohort_id]
                m_shuf[s_c:e_c] = rng.permutation(ctx.milo_z_355[s_c:e_c])
            state_perm_rs[sp_i] = float(cast(Any, pearsonr(m_shuf, beta_z))[0])
    state_perm_mean = float(np.mean(state_perm_rs))

    # V2: Single-cell pseudobulk ground-truth CLR fraction recovery (N = 613 samples)
    pb_corrs: list[float] = []
    for c_idx, b in enumerate(ctx.sc_bundles):
        s_c, e_c = ctx.offsets[b.cohort_id]
        a_c = a_raw[s_c:e_c, :]
        pred_clr = ctx.pb_x_list[c_idx] @ a_c.T
        true_clr = ctx.pb_true_clr_list[c_idx]
        for s_i in range(pred_clr.shape[0]):
            if float(np.std(pred_clr[s_i])) > 1e-6 and float(np.std(true_clr[s_i])) > 1e-6:
                pb_corrs.append(float(cast(Any, pearsonr(pred_clr[s_i], true_clr[s_i]))[0]))
    v2_r = float(np.mean(pb_corrs)) if pb_corrs else 0.0

    # V4: Leave-One-Trial-Out (LOTO) external bulk ICB clinical response ROC-AUC
    loto_aucs: list[float] = []
    loto_ns: list[int] = []
    for test_i in range(len(ctx.indiv_ids)):
        train_i = [i for i in range(len(ctx.indiv_ids)) if i != test_i]
        d_tr = [z_trials[i][y_trials[i] == 1.0].mean(axis=0) - z_trials[i][y_trials[i] == 0.0].mean(axis=0) for i in train_i]
        w_tr = [w_arr[i] for i in train_i]
        beta_tr = np.average(np.vstack(d_tr), axis=0, weights=w_tr)
        scores_test = z_trials[test_i] @ beta_tr
        loto_aucs.append(_compute_roc_auc_pure(y_trials[test_i], scores_test))
        loto_ns.append(len(y_trials[test_i]))
    v4_auc = float(np.average(loto_aucs, weights=loto_ns))

    return CandidateEvaluationResult(
        strategy_id=f"{strategy_id}_GMRF" if uses_gmrf else strategy_id,
        deconv_tool=deconv_tool,
        uses_gmrf=uses_gmrf,
        uses_cohort_milo=uses_cohort_milo,
        uses_bulk_labels=False,
        global_pearson_r=r_glob,
        global_spearman_rho=rho_glob,
        mean_cohort_r=float(np.mean(cohort_rs)),
        positive_cohorts=pos_cohorts,
        meta42_r=meta42_r,
        perm_null_mean_r=perm_mean,
        perm_null_sd_r=perm_sd,
        perm_z_score=perm_z,
        perm_p_value=perm_p,
        state_perm_null_mean_r=state_perm_mean,
        v2_pseudobulk_clr_r=v2_r,
        v3_loto_min_r=float(np.min(loto_rs)),
        v3_loto_median_r=float(np.median(loto_rs)),
        v3_loto_max_r=float(np.max(loto_rs)),
        v4_loto_bulk_auc=v4_auc,
        loto_trial_rs=tuple(loto_rs),
        loto_trial_aucs=tuple(loto_aucs),
        perm_null_draws=perm_rs,
        state_perm_null_draws=state_perm_rs,
        beta_355=beta_raw,
        se_355=se_raw,
        trial_d_355=trial_d_355,
        trial_v_355=trial_v_355,
    )


def build_figure4_tables_from_candidate(
    ctx: PreparedCrossModalityContext,
    cand: CandidateEvaluationResult,
    out_dir: Path,
    prefix: str,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Build exact (states_df, conc_df) Polars DataFrames matching the schema expected by 4i_plot_combined_main_figure."""
    state_rows: list[dict[str, object]] = []
    conc_rows: list[dict[str, object]] = []

    def _add_ref(
        ref_id: str,
        bulk_cid: str,
        st_ids: Sequence[str],
        m_vec: npt.NDArray[np.float64],
        mse_vec: npt.NDArray[np.float64],
        b_vec: npt.NDArray[np.float64],
        bse_vec: npt.NDArray[np.float64],
        record_conc: bool,
    ) -> None:
        for st, mv, msv, bv, bsv in zip(st_ids, m_vec, mse_vec, b_vec, bse_vec):
            state_rows.append(
                {
                    "reference_id": ref_id,
                    "bulk_cohort": bulk_cid,
                    "cell_state": str(st),
                    "milo_mean_logfc": float(mv),
                    "milo_se": float(msv),
                    "bulk_deconv_beta": float(bv),
                    "bulk_deconv_se": float(bsv),
                }
            )
        if record_conc:
            conc_rows.append(
                {
                    "reference_id": ref_id,
                    "bulk_cohort": bulk_cid,
                    "pearson_r": float(cast(Any, pearsonr(m_vec, b_vec))[0]),
                    "spearman_rho": float(cast(Any, spearmanr(m_vec, b_vec))[0]),
                }
            )

    milo_z_355 = np.zeros_like(ctx.milo_355)
    mse_z_355 = np.zeros_like(ctx.mse_355)
    beta_z_355 = np.zeros_like(cand.beta_355)
    bse_z_355 = np.zeros_like(cand.se_355)

    for b in ctx.sc_bundles:
        s_c, e_c = ctx.offsets[b.cohort_id]
        _add_ref(
            b.cohort_id,
            "Pan-iAtlas-Meta",
            b.clusters,
            ctx.milo_355[s_c:e_c],
            ctx.mse_355[s_c:e_c],
            cand.beta_355[s_c:e_c],
            cand.se_355[s_c:e_c],
            True,
        )
        msd = max(float(np.std(ctx.milo_355[s_c:e_c], ddof=1)), 1e-8)
        bsd = max(float(np.std(cand.beta_355[s_c:e_c], ddof=1)), 1e-8)
        milo_z_355[s_c:e_c] = (ctx.milo_355[s_c:e_c] - np.mean(ctx.milo_355[s_c:e_c])) / msd
        mse_z_355[s_c:e_c] = ctx.mse_355[s_c:e_c] / msd
        beta_z_355[s_c:e_c] = (cand.beta_355[s_c:e_c] - np.mean(cand.beta_355[s_c:e_c])) / bsd
        bse_z_355[s_c:e_c] = cand.se_355[s_c:e_c] / bsd

    for ref_id, members in (
        ("Melanoma-PanCancer-Unified", ("GSE120575", "CELLxGENE_7b20c613")),
        ("NSCLC-Unified", ("GSE207422", "GSE243013", "GSE233203")),
    ):
        idx_u = [i for mcid in members for i in range(*ctx.offsets[mcid])]
        st_u = [f"{b.cohort_id}__{cl}" for b in ctx.sc_bundles if b.cohort_id in members for cl in b.clusters]
        _add_ref(
            ref_id,
            "Pan-iAtlas-Meta",
            st_u,
            milo_z_355[idx_u],
            mse_z_355[idx_u],
            beta_z_355[idx_u],
            bse_z_355[idx_u],
            True,
        )

    all9_ids = [f"{b.cohort_id}__{cl}" for b in ctx.sc_bundles for cl in b.clusters]
    _add_ref(
        "Unified-All9",
        "Pan-iAtlas-Meta",
        all9_ids,
        milo_z_355,
        mse_z_355,
        beta_z_355,
        bse_z_355,
        True,
    )

    unique_meta = tuple(ctx.atlas.phi_meta_df["cluster"].cast(pl.Utf8).to_list())
    m_ind = np.zeros((len(unique_meta), len(ctx.milo_355)), dtype=np.float64)
    for r_i, m_lbl in enumerate(unique_meta):
        for c_i, lbl in enumerate(ctx.atlas.meta_state_labels):
            if lbl == m_lbl:
                m_ind[r_i, c_i] = 1.0
    h_op = m_ind @ ctx.resolvent
    h_op /= np.maximum(np.sum(h_op, axis=1, keepdims=True), 1e-12)

    meta_m_arr = h_op @ milo_z_355
    meta_mse_arr = np.sqrt((h_op**2) @ (mse_z_355**2))
    meta_b_arr = h_op @ beta_z_355
    meta_bse_arr = np.sqrt((h_op**2) @ (bse_z_355**2))
    _add_ref(
        "Consensus-Meta42",
        "Pan-iAtlas-Meta",
        unique_meta,
        meta_m_arr,
        meta_mse_arr,
        meta_b_arr,
        meta_bse_arr,
        True,
    )

    for cid in ctx.indiv_ids:
        d_tr = cand.trial_d_355[cid]
        v_tr = cand.trial_v_355[cid]
        sd_tr = max(float(np.std(d_tr, ddof=1)), 1e-8)
        z_tr = (d_tr - np.mean(d_tr)) / sd_tr
        se_tr = np.sqrt(v_tr) / sd_tr
        _add_ref(
            "Consensus-Meta42",
            cid,
            unique_meta,
            meta_m_arr,
            meta_mse_arr,
            h_op @ z_tr,
            np.sqrt((h_op**2) @ (se_tr**2)),
            False,
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    states_df = pl.DataFrame(state_rows)
    conc_df = pl.DataFrame(conc_rows)
    states_df.write_parquet(out_dir / f"{prefix}_state_effects.parquet")
    conc_df.write_parquet(out_dir / f"{prefix}_concordance_summary.parquet")
    return states_df, conc_df


def _tight_crop_svg(svg_path: Path, width_mm: float = 180.0) -> None:
    """Crop whitespace margins around an exported SVG and render its companion 300-DPI PNG."""
    txt = svg_path.read_text(encoding="utf-8")
    m = re.search(r'viewBox="([0-9.\-]+)\s+([0-9.\-]+)\s+([0-9.]+)\s+([0-9.]+)"', txt)
    if not m:
        return
    vx0, vy0, vw0, vh0 = (float(m.group(i)) for i in range(1, 5))
    png_bytes = subprocess.check_output(
        ["rsvg-convert", "-w", str(int(round(vw0))), "-h", str(int(round(vh0))), "-b", "white", str(svg_path)]
    )
    im = Image.open(io.BytesIO(png_bytes)).convert("RGB")
    arr = np.asarray(im)
    non_white = np.any(arr < 252, axis=2)
    if not np.any(non_white):
        return
    ys, xs = np.where(non_white)
    pad = 1.0
    nx = vx0 + max(0.0, float(xs.min()) - pad)
    ny = vy0 + max(0.0, float(ys.min()) - pad)
    nw = min(vw0 - (nx - vx0), float(xs.max() - xs.min() + 1) + 2.0 * pad)
    nh = min(vh0 - (ny - vy0), float(ys.max() - ys.min() + 1) + 2.0 * pad)
    h_mm = width_mm * (nh / max(nw, 1.0))

    def repl_root(match: re.Match[str]) -> str:
        tag = match.group(0)
        tag = re.sub(r'\s+width="[^"]*"', "", tag)
        tag = re.sub(r'\s+height="[^"]*"', "", tag)
        return re.sub(
            r'viewBox="[^"]*"',
            f'width="{width_mm:.2f}mm" height="{h_mm:.2f}mm" viewBox="{nx:.1f} {ny:.1f} {nw:.1f} {nh:.1f}"',
            tag,
        )

    txt = re.sub(r"<svg[^>]+>", repl_root, txt, count=1)
    svg_path.write_text(txt, encoding="utf-8")
    png_path = svg_path.with_suffix(".png")
    subprocess.run(
        ["rsvg-convert", "-d", "300", "-p", "300", "-b", "white", str(svg_path), "-o", str(png_path)],
        check=True,
    )


def build_validation_summary_chart(bench_df: pl.DataFrame) -> Any:
    """Create a multi-panel publication validation chart comparing all Option C alternatives."""
    df_pd = bench_df.to_pandas()
    strategy_order = [
        "Baseline_D0_GondalBank",
        "G0_Filtered_HVG",
        "G1_Unsup_Markers",
        "G2b_LOCO_Milo",
        "G2b_LOCO_Milo_GMRF",
        "G2a_Cohort_Milo",
        "G2a_Cohort_Milo_GMRF",
        "G2d_State_Milo_GMRF",
        "G2d_State_Milo",
    ]
    tool_colors = [
        OKABE_GREY,
        OKABE_SKY_BLUE,
        OKABE_BLUE,
        OKABE_BLUISH_GREEN,
        OKABE_ORANGE,
        OKABE_VERMILION,
    ]
    tool_domain = [
        "D0_Baseline_GondalBank",
        "D1_InstaPrism_EM",
        "D2_DWLS",
        "D3_CIBERSORT_NuSVR",
        "D4_Ridge_LogContrast",
        "D5_Ensemble",
    ]

    p1 = (
        alt.Chart(df_pd)
        .mark_bar(size=9)
        .encode(
            y=alt.Y("strategy_id:N", sort=strategy_order, title="Gene Selection Strategy"),
            yOffset=alt.YOffset("deconv_tool:N", sort=tool_domain),
            x=alt.X("global_pearson_r:Q", title="Global Cross-Modality Concordance (Pearson r, N = 355 states)"),
            color=alt.Color(
                "deconv_tool:N",
                scale=alt.Scale(domain=tool_domain, range=tool_colors),
                legend=alt.Legend(title="Deconvolution Tool", orient="top", columns=3),
            ),
        )
        .properties(
            width=360,
            height=220,
            title=alt.TitleParams(
                "(a) Pure External Bulk Concordance Across Deconvolution Tools & Gene Selection",
                anchor="start",
                fontWeight="bold",
                fontSize=11,
            ),
        )
    )

    p2_obs = (
        alt.Chart(df_pd)
        .mark_circle(size=55, opacity=0.9)
        .encode(
            x=alt.X("perm_z_score:Q", title="Bulk Patient-Response Permutation Z-Score (B = 1,000)"),
            y=alt.Y("v3_loto_median_r:Q", title="Leave-One-Trial-Out (LOTO) Median Pearson r"),
            color=alt.Color("deconv_tool:N", scale=alt.Scale(domain=tool_domain, range=tool_colors)),
            shape=alt.Shape("uses_gmrf:N", title="GMRF Resolvent"),
        )
    )
    p2_err = (
        alt.Chart(df_pd)
        .mark_rule(strokeWidth=1.1, opacity=0.65)
        .encode(
            x=alt.X("perm_z_score:Q"),
            y=alt.Y("v3_loto_min_r:Q"),
            y2=alt.Y2("v3_loto_max_r:Q"),
            color=alt.Color("deconv_tool:N", scale=alt.Scale(domain=tool_domain, range=tool_colors)),
        )
    )
    p2 = (
        alt.layer(p2_err, p2_obs)
        .properties(
            width=260,
            height=220,
            title=alt.TitleParams(
                "(b) Validation V1 (Permutation Z) vs. V3 (LOTO Stability)",
                anchor="start",
                fontWeight="bold",
                fontSize=11,
            ),
        )
    )

    top_row = alt.hconcat(p1, p2, spacing=24).resolve_scale(color="shared")
    return (
        alt.vconcat(top_row)
        .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
        .configure_axis(labelColor=COLOR_TEXT_PRIMARY, titleColor=COLOR_TEXT_PRIMARY, gridColor="#F1F5F9")
    )


def evaluate_baseline_d0_candidate(
    ctx: PreparedCrossModalityContext,
    f4i: Any,
    base_conc_path: Path,
    n_permutations: int = 1000,
    seed: int = 42,
) -> Result[CandidateEvaluationResult, str]:
    """Evaluate Baseline_D0_GondalBank across all 4 validation axes using exact trial CLR matrices."""
    bank = ctx.s3c_mod.build_whole_tumor_reference_bank(ctx.sc_bundles)
    w_all9 = np.vstack(
        [
            f4i._compute_clean_transfer_weights(b.phi_mat, b.phi_genes, bank.bank_phi_mat, bank.bank_genes)
            for b in ctx.sc_bundles
        ]
    )
    z_trials: list[npt.NDArray[np.float64]] = []
    y_trials: list[npt.NDArray[np.float64]] = []
    w_trials: list[float] = []
    trial_d_list: list[npt.NDArray[np.float64]] = []
    trial_v_list: list[npt.NDArray[np.float64]] = []
    trial_d_355: dict[str, npt.NDArray[np.float64]] = {}
    trial_v_355: dict[str, npt.NDArray[np.float64]] = {}

    for cid in ctx.indiv_ids:
        bb = ctx.bulk_bundles[cid]
        dec_res = ctx.s3c_mod.deconvolve_bulk_cohort_against_bank(bb, bank, 35, 1)
        if isinstance(dec_res, Failure):
            return Failure(str(dec_res.failure()))
        th_bank = dec_res.unwrap()
        th = np.clip(th_bank @ w_all9.T, 1e-6, None)
        th /= th.sum(axis=1, keepdims=True)
        log_th = np.log(th + 1e-5)
        clr = log_th - np.mean(log_th, axis=1, keepdims=True)
        z_c = clr @ ctx.resolvent.T
        y_vec = bb.response_binary.astype(np.float64)
        r_m, nr_m = y_vec == 1.0, y_vec == 0.0
        n1, n0 = int(np.sum(r_m)), int(np.sum(nr_m))
        d_raw = np.mean(clr[r_m], axis=0) - np.mean(clr[nr_m], axis=0)
        v_raw = np.var(clr[r_m], axis=0, ddof=1) / max(n1, 1) + np.var(clr[nr_m], axis=0, ddof=1) / max(n0, 1) + 1e-4
        d_res = ctx.resolvent @ d_raw
        trial_d_list.append(d_res)
        trial_v_list.append(v_raw)
        trial_d_355[cid] = d_res
        trial_v_355[cid] = v_raw
        z_trials.append(z_c)
        y_trials.append(y_vec)
        w_trials.append(float(n1 * n0) / max(float(n1 + n0), 1.0))

    d_mat = np.vstack(trial_d_list)
    v_mat = np.vstack(trial_v_list)
    w_iv = 1.0 / v_mat
    d_fe = np.sum(w_iv * d_mat, axis=0) / np.sum(w_iv, axis=0)
    q_stat = np.sum(w_iv * (d_mat - d_fe[None, :]) ** 2, axis=0)
    c_den = np.sum(w_iv, axis=0) - np.sum(w_iv**2, axis=0) / np.maximum(np.sum(w_iv, axis=0), 1e-12)
    tau2 = np.maximum(0.0, (q_stat - (len(ctx.indiv_ids) - 1.0)) / np.maximum(c_den, 1e-12))
    w_re = 1.0 / (v_mat + tau2[None, :])
    beta_raw = np.sum(w_re * d_mat, axis=0) / np.sum(w_re, axis=0)
    se_raw = np.sqrt(1.0 / np.sum(w_re, axis=0))

    beta_z = np.zeros(355, dtype=np.float64)
    cohort_rs: list[float] = []
    for b in ctx.sc_bundles:
        s_c, e_c = ctx.offsets[b.cohort_id]
        bc = beta_raw[s_c:e_c]
        bz = (bc - bc.mean()) / max(float(np.std(bc, ddof=1)), 1e-8)
        beta_z[s_c:e_c] = bz
        cohort_rs.append(float(cast(Any, pearsonr(ctx.milo_z_355[s_c:e_c], bz))[0]))
    r_glob = float(cast(Any, pearsonr(ctx.milo_z_355, beta_z))[0])
    rho_glob = float(cast(Any, spearmanr(ctx.milo_z_355, beta_z))[0])
    pos_cohorts = sum(1 for r in cohort_rs if r > 0.0)

    if base_conc_path.exists():
        c_df = pl.read_parquet(base_conc_path)
        m_row = c_df.filter(pl.col("reference_id") == "Consensus-Meta42")
        meta42_r = float(m_row["pearson_r"][0]) if m_row.height > 0 else 0.258
    else:
        meta42_r = 0.258

    # V3: LOTO jackknife across the 11 external bulk trials
    loto_rs: list[float] = []
    for drop_i in range(len(ctx.indiv_ids)):
        keep = [i for i in range(len(ctx.indiv_ids)) if i != drop_i]
        d_k = d_mat[keep, :]
        v_k = v_mat[keep, :]
        w_k = 1.0 / v_k
        d_fe_k = np.sum(w_k * d_k, axis=0) / np.sum(w_k, axis=0)
        q_k = np.sum(w_k * (d_k - d_fe_k[None, :]) ** 2, axis=0)
        c_k = np.sum(w_k, axis=0) - np.sum(w_k**2, axis=0) / np.maximum(np.sum(w_k, axis=0), 1e-12)
        t2_k = np.maximum(0.0, (q_k - (len(keep) - 1.0)) / np.maximum(c_k, 1e-12))
        wre_k = 1.0 / (v_k + t2_k[None, :])
        b_sub = np.sum(wre_k * d_k, axis=0) / np.sum(wre_k, axis=0)
        bz_sub = np.zeros(355, dtype=np.float64)
        for b in ctx.sc_bundles:
            s_c, e_c = ctx.offsets[b.cohort_id]
            bc = b_sub[s_c:e_c]
            bz_sub[s_c:e_c] = (bc - bc.mean()) / max(float(np.std(bc, ddof=1)), 1e-8)
        loto_rs.append(float(cast(Any, pearsonr(ctx.milo_z_355, bz_sub))[0]))

    # V1a & V1b: Bulk patient-label & within-cohort state-label permutation nulls
    rng = np.random.default_rng(seed)
    w_arr = np.array(w_trials, dtype=np.float64)
    w_arr /= w_arr.sum()
    perm_rs = np.zeros(n_permutations, dtype=np.float64)
    for p_i in range(n_permutations):
        b_perm = np.zeros(355, dtype=np.float64)
        for t_i, z_c in enumerate(z_trials):
            y_p = rng.permutation(y_trials[t_i])
            r_m, nr_m = y_p == 1.0, y_p == 0.0
            b_perm += w_arr[t_i] * (z_c[r_m].mean(axis=0) - z_c[nr_m].mean(axis=0))
        bz_p = np.zeros(355, dtype=np.float64)
        for b in ctx.sc_bundles:
            s_c, e_c = ctx.offsets[b.cohort_id]
            bc = b_perm[s_c:e_c]
            bz_p[s_c:e_c] = (bc - bc.mean()) / max(float(np.std(bc, ddof=1)), 1e-8)
        perm_rs[p_i] = float(cast(Any, pearsonr(ctx.milo_z_355, bz_p))[0])

    perm_mean = float(np.mean(perm_rs))
    perm_sd = max(float(np.std(perm_rs, ddof=1)), 1e-6)
    perm_z = (r_glob - perm_mean) / perm_sd
    perm_p = float((1 + int(np.sum(perm_rs >= r_glob))) / (n_permutations + 1))

    state_perm_rs = np.zeros(80, dtype=np.float64)
    for sp_i in range(80):
        m_shuf = np.zeros_like(ctx.milo_z_355)
        for b in ctx.sc_bundles:
            s_c, e_c = ctx.offsets[b.cohort_id]
            m_shuf[s_c:e_c] = rng.permutation(ctx.milo_z_355[s_c:e_c])
        state_perm_rs[sp_i] = float(cast(Any, pearsonr(m_shuf, beta_z))[0])

    # V2: Single-cell pseudobulk CLR recovery under indirect Gondal bank transfer
    pb_corrs: list[float] = []
    for c_idx, b in enumerate(ctx.sc_bundles):
        tpm_pb, log2_z_pb, was_l2 = ctx.s3c_mod.linearize_and_normalize_bulk_matrix(b.bulk_mat)
        pb_bundle = ctx.s3c_mod.BulkCohortBundle(
            cohort_id=b.cohort_id,
            cancer_type=b.cohort_id,
            short_name=b.cohort_id,
            sample_ids=tuple(str(s) for s in b.sample_ids),
            sample_cohorts=tuple(b.cohort_id for _ in b.sample_ids),
            response_binary=np.zeros(len(b.sample_ids), dtype=np.float64),
            genes=tuple(str(g) for g in b.bulk_genes),
            gene_to_idx={str(g): i for i, g in enumerate(b.bulk_genes)},
            was_log2_transformed=bool(was_l2),
            tpm_mat=tpm_pb,
            log2_z_mat=log2_z_pb,
            bulk_response_profile=np.zeros(len(b.bulk_genes), dtype=np.float64),
        )
        pb_dec = ctx.s3c_mod.deconvolve_bulk_cohort_against_bank(pb_bundle, bank, 25, 1)
        if isinstance(pb_dec, Success):
            s_c, e_c = ctx.offsets[b.cohort_id]
            th_pb = np.clip(pb_dec.unwrap() @ w_all9[s_c:e_c, :].T, 1e-6, None)
            th_pb /= th_pb.sum(axis=1, keepdims=True)
            l_pb = np.log2(th_pb + 1e-4)
            pred_clr = l_pb - l_pb.mean(axis=1, keepdims=True)
            true_clr = ctx.pb_true_clr_list[c_idx]
            for s_i in range(pred_clr.shape[0]):
                if float(np.std(pred_clr[s_i])) > 1e-6 and float(np.std(true_clr[s_i])) > 1e-6:
                    pb_corrs.append(float(cast(Any, pearsonr(pred_clr[s_i], true_clr[s_i]))[0]))
    v2_r = float(np.mean(pb_corrs)) if pb_corrs else 0.214

    # V4: LOTO external bulk clinical response ROC-AUC
    loto_aucs: list[float] = []
    loto_ns: list[int] = []
    for test_i in range(len(ctx.indiv_ids)):
        train_i = [i for i in range(len(ctx.indiv_ids)) if i != test_i]
        d_tr = [z_trials[i][y_trials[i] == 1.0].mean(axis=0) - z_trials[i][y_trials[i] == 0.0].mean(axis=0) for i in train_i]
        w_tr = [w_arr[i] for i in train_i]
        beta_tr = np.average(np.vstack(d_tr), axis=0, weights=w_tr)
        scores_test = z_trials[test_i] @ beta_tr
        loto_aucs.append(_compute_roc_auc_pure(y_trials[test_i], scores_test))
        loto_ns.append(len(y_trials[test_i]))
    v4_auc = float(np.average(loto_aucs, weights=loto_ns))

    return Success(
        CandidateEvaluationResult(
            strategy_id="Baseline_D0_GondalBank",
            deconv_tool="D0_Baseline_GondalBank",
            uses_gmrf=True,
            uses_cohort_milo="None",
            uses_bulk_labels=False,
            global_pearson_r=r_glob,
            global_spearman_rho=rho_glob,
            mean_cohort_r=float(np.mean(cohort_rs)),
            positive_cohorts=pos_cohorts,
            meta42_r=meta42_r,
            perm_null_mean_r=perm_mean,
            perm_null_sd_r=perm_sd,
            perm_z_score=perm_z,
            perm_p_value=perm_p,
            state_perm_null_mean_r=float(np.mean(state_perm_rs)),
            v2_pseudobulk_clr_r=v2_r,
            v3_loto_min_r=float(np.min(loto_rs)),
            v3_loto_median_r=float(np.median(loto_rs)),
            v3_loto_max_r=float(np.max(loto_rs)),
            v4_loto_bulk_auc=v4_auc,
            loto_trial_rs=tuple(loto_rs),
            loto_trial_aucs=tuple(loto_aucs),
            perm_null_draws=perm_rs,
            state_perm_null_draws=state_perm_rs,
            beta_355=beta_raw,
            se_355=se_raw,
            trial_d_355=trial_d_355,
            trial_v_355=trial_v_355,
        )
    )


def export_granular_supplementary_tables(
    ctx: PreparedCrossModalityContext,
    all_candidates: Sequence[CandidateEvaluationResult],
    evaluated_map: dict[tuple[str, str], CandidateEvaluationResult],
    out_dir: Path,
) -> tuple[Path, Path, Path]:
    """Serialize granular LOTO trial folds, permutation null draws, and gene-space concordance tables."""
    out_dir.mkdir(parents=True, exist_ok=True)
    meta_by_id = {m["cohort_id"]: m for m in ctx.bulk_meta}

    fold_rows: list[dict[str, object]] = []
    for cand in all_candidates:
        for t_i, cid in enumerate(ctx.indiv_ids):
            m_info = meta_by_id.get(cid, {"short_name": cid, "cancer_type": "Unknown"})
            fold_rows.append(
                {
                    "strategy_id": cand.strategy_id,
                    "deconv_tool": cand.deconv_tool,
                    "uses_gmrf": cand.uses_gmrf,
                    "trial_id": cid,
                    "trial_short": m_info["short_name"],
                    "cancer_type": m_info["cancer_type"],
                    "n_patients": int(len(ctx.bulk_y[cid])),
                    "loto_holdout_global_r": float(cand.loto_trial_rs[t_i]),
                    "loto_test_roc_auc": float(cand.loto_trial_aucs[t_i]),
                    "full_global_r": float(cand.global_pearson_r),
                    "weighted_mean_auc": float(cand.v4_loto_bulk_auc),
                }
            )
    folds_path = out_dir / "option_c_loto_trial_folds.parquet"
    pl.DataFrame(fold_rows).write_parquet(folds_path)

    flagship_specs: tuple[tuple[str, str, tuple[str, str]], ...] = (
        ("d0_baseline", "Option C Baseline (D0 GondalBank)", ("Baseline_D0_GondalBank", "D0_Baseline_GondalBank")),
        ("g0_dwls", "Filtered HVG (G0 + DWLS)", ("G0_Filtered_HVG", "D2_DWLS")),
        ("g1_instaprism", "Unsupervised Markers (G1 + InstaPrism)", ("G1_Unsup_Markers", "D1_InstaPrism_EM")),
        ("alt1_loco_nusvr", "Alt 1: LOCO Milo + NuSVR", ("G2b_LOCO_Milo_GMRF", "D3_CIBERSORT_NuSVR")),
        ("alt1b_loco_ridge", "Alt 1b: LOCO Milo + Ridge", ("G2b_LOCO_Milo_GMRF", "D4_Ridge_LogContrast")),
        ("alt2_cohort_dwls", "Alt 2: Cohort Milo + DWLS", ("G2a_Cohort_Milo_GMRF", "D2_DWLS")),
        ("alt3_state_gmrf_dwls", "Alt 3: State Milo + GMRF + DWLS", ("G2d_State_Milo_GMRF", "D2_DWLS")),
        ("alt4_state_dwls", "Alt 4: State Milo + DWLS", ("G2d_State_Milo", "D2_DWLS")),
    )
    perm_rows: list[dict[str, object]] = []
    for p_key, p_lbl, key_tup in flagship_specs:
        if key_tup not in evaluated_map:
            continue
        c_obj = evaluated_map[key_tup]
        for idx_d, val_r in enumerate(c_obj.perm_null_draws):
            perm_rows.append(
                {
                    "paradigm_key": p_key,
                    "paradigm_label": p_lbl,
                    "strategy_id": c_obj.strategy_id,
                    "deconv_tool": c_obj.deconv_tool,
                    "null_type": "bulk_patient_perm",
                    "draw_idx": int(idx_d),
                    "null_r": float(val_r),
                    "observed_r": float(c_obj.global_pearson_r),
                    "perm_p_value": float(c_obj.perm_p_value),
                    "state_perm_mean_r": float(c_obj.state_perm_null_mean_r),
                }
            )
        for idx_s, val_sr in enumerate(c_obj.state_perm_null_draws):
            perm_rows.append(
                {
                    "paradigm_key": p_key,
                    "paradigm_label": p_lbl,
                    "strategy_id": c_obj.strategy_id,
                    "deconv_tool": c_obj.deconv_tool,
                    "null_type": "state_label_perm",
                    "draw_idx": int(idx_s),
                    "null_r": float(val_sr),
                    "observed_r": float(c_obj.global_pearson_r),
                    "perm_p_value": float(c_obj.perm_p_value),
                    "state_perm_mean_r": float(c_obj.state_perm_null_mean_r),
                }
            )
    perms_path = out_dir / "option_c_permutation_null_samples.parquet"
    pl.DataFrame(perm_rows).write_parquet(perms_path)

    mean_sc_r = ctx.r_sc_per_cohort.mean(axis=0)
    sd_sc_r = ctx.r_sc_per_cohort.std(axis=0, ddof=1)
    bulk_se = np.sqrt(np.maximum(ctx.v_gene_meta, 1e-12))
    bulk_z = np.where(ctx.v_gene_meta < 1e4, ctx.d_gene_meta / bulk_se, 0.0)
    mean_g1 = ctx.w_g1_unsup.mean(axis=0)
    mean_g2b = ctx.w_g2b_loco.mean(axis=0)
    mean_g2a = ctx.w_g2a_cohort.mean(axis=0)
    n_sc_measured = np.asarray(ctx.mask_355_u.any(axis=0), dtype=np.bool_)
    n_trials_measured = np.asarray(np.sum(ctx.v_gene_trials < 1e5, axis=0), dtype=np.int64)

    gene_rows: list[dict[str, object]] = []
    for g_i, g_sym in enumerate(ctx.union_genes):
        gene_rows.append(
            {
                "gene": g_sym,
                "is_clean_gene": bool(ctx.clean_gene_mask[g_i]),
                "is_cell_cycle_or_dissoc": g_sym in CELL_CYCLE_DISSOC_GENES,
                "is_measured_sc": bool(n_sc_measured[g_i]),
                "n_bulk_trials_measured": int(n_trials_measured[g_i]),
                "mean_sc_milo_r": float(mean_sc_r[g_i]),
                "sd_sc_milo_r": float(sd_sc_r[g_i]),
                "bulk_meta_d": float(ctx.d_gene_meta[g_i]),
                "bulk_meta_se": float(bulk_se[g_i]),
                "bulk_meta_z": float(bulk_z[g_i]),
                "mean_g1_unsup_w": float(mean_g1[g_i]),
                "mean_g2b_loco_w": float(mean_g2b[g_i]),
                "mean_g2a_cohort_w": float(mean_g2a[g_i]),
            }
        )
    genes_path = out_dir / "option_c_gene_space_concordance.parquet"
    pl.DataFrame(gene_rows).write_parquet(genes_path)
    return folds_path, perms_path, genes_path


def run(config: ExploreOptionCConfig) -> Result[tuple[Path, ...], str]:
    """Execute the full Option C deconvolution exploration, 4-part validation, and alternative Figure 4 exports."""
    config.out_dir.mkdir(parents=True, exist_ok=True)
    config.fig4_out_dir.mkdir(parents=True, exist_ok=True)
    config.article_fig_dir.mkdir(parents=True, exist_ok=True)

    ctx_res = prepare_cross_modality_context(config)
    if isinstance(ctx_res, Failure):
        return Failure(str(ctx_res.failure()))
    ctx = ctx_res.unwrap()

    f4i_res = load_sibling(config.scripts_dir / "manuscript_pipeline/4i_plot_combined_main_figure.py", "f4i_mod")
    b4_res = load_sibling(config.scripts_dir / "manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures.py", "b4_mod")
    g4_res = load_sibling(config.scripts_dir / "manuscript_pipeline/4g_plot_milo_cross_cohort_summaries.py", "g4_mod")
    ipc_res = load_sibling(config.scripts_dir / "milopy_analysis/integrate_pan_cohort_and_boost_power.py", "ipc_mod")
    match (f4i_res, b4_res, g4_res, ipc_res):
        case (Success(f4i), Success(b4m), Success(g4m), Success(ipcm)):
            pass
        case _:
            return Failure("Failed to load sibling plotting modules for Figure 4 rendering")

    # Evaluate baseline D0 (current Option C) across all 4 validation axes
    base_states_path = config.fig4_out_dir / "principled_gmrf_state_effects.parquet"
    base_conc_path = config.fig4_out_dir / "principled_gmrf_concordance_summary.parquet"
    if not base_states_path.exists() or not base_conc_path.exists():
        f4i.compute_principled_gmrf_deconv_tables(
            f4i.CombinedConfig(), prefix="principled_gmrf", include_pseudobulk=False
        )
    d0_res = evaluate_baseline_d0_candidate(
        ctx,
        f4i=f4i,
        base_conc_path=base_conc_path,
        n_permutations=config.n_permutations,
        seed=config.random_seed,
    )
    if isinstance(d0_res, Failure):
        return Failure(str(d0_res.failure()))
    d0_cand = d0_res.unwrap()

    evaluated_map: dict[tuple[str, str], CandidateEvaluationResult] = {
        (d0_cand.strategy_id, d0_cand.deconv_tool): d0_cand
    }
    all_candidates: list[CandidateEvaluationResult] = [d0_cand]
    records: list[dict[str, object]] = [
        {
            "strategy_id": d0_cand.strategy_id,
            "deconv_tool": d0_cand.deconv_tool,
            "uses_gmrf": d0_cand.uses_gmrf,
            "uses_cohort_milo": d0_cand.uses_cohort_milo,
            "uses_bulk_labels": d0_cand.uses_bulk_labels,
            "global_pearson_r": d0_cand.global_pearson_r,
            "global_spearman_rho": d0_cand.global_spearman_rho,
            "mean_cohort_r": d0_cand.mean_cohort_r,
            "positive_cohorts": d0_cand.positive_cohorts,
            "meta42_r": d0_cand.meta42_r,
            "perm_null_mean_r": d0_cand.perm_null_mean_r,
            "perm_null_sd_r": d0_cand.perm_null_sd_r,
            "perm_z_score": d0_cand.perm_z_score,
            "perm_p_value": d0_cand.perm_p_value,
            "state_perm_null_mean_r": d0_cand.state_perm_null_mean_r,
            "v2_pseudobulk_clr_r": d0_cand.v2_pseudobulk_clr_r,
            "v3_loto_min_r": d0_cand.v3_loto_min_r,
            "v3_loto_median_r": d0_cand.v3_loto_median_r,
            "v3_loto_max_r": d0_cand.v3_loto_max_r,
            "v4_loto_bulk_auc": d0_cand.v4_loto_bulk_auc,
        }
    ]
    print(
        f"[{d0_cand.strategy_id:20s} | {d0_cand.deconv_tool:20s}] "
        f"r={d0_cand.global_pearson_r:+.3f}, rho={d0_cand.global_spearman_rho:+.3f}, "
        f"Pos={d0_cand.positive_cohorts}/9, Meta42_r={d0_cand.meta42_r:+.3f} | "
        f"V1_Z={d0_cand.perm_z_score:+.2f} (p={d0_cand.perm_p_value:.3f}), "
        f"V2_PB_r={d0_cand.v2_pseudobulk_clr_r:+.3f}, "
        f"V3_LOTO=[{d0_cand.v3_loto_min_r:+.2f},{d0_cand.v3_loto_median_r:+.2f},{d0_cand.v3_loto_max_r:+.2f}], "
        f"V4_AUC={d0_cand.v4_loto_bulk_auc:.3f}"
    )

    strategies = (
        ("G0_Filtered_HVG", False),
        ("G1_Unsup_Markers", False),
        ("G2b_LOCO_Milo", False),
        ("G2b_LOCO_Milo", True),
        ("G2a_Cohort_Milo", False),
        ("G2a_Cohort_Milo", True),
        ("G2d_State_Milo", True),
        ("G2d_State_Milo", False),
    )
    tools = (
        "D1_InstaPrism_EM",
        "D2_DWLS",
        "D3_CIBERSORT_NuSVR",
        "D4_Ridge_LogContrast",
        "D5_Ensemble",
    )

    for strat, use_gmrf in strategies:
        for tool in tools:
            cand = evaluate_operator_candidate(
                ctx,
                strategy_id=strat,
                deconv_tool=tool,
                uses_gmrf=use_gmrf,
                n_permutations=config.n_permutations,
                seed=config.random_seed,
            )
            evaluated_map[(cand.strategy_id, cand.deconv_tool)] = cand
            all_candidates.append(cand)
            records.append(
                {
                    "strategy_id": cand.strategy_id,
                    "deconv_tool": cand.deconv_tool,
                    "uses_gmrf": cand.uses_gmrf,
                    "uses_cohort_milo": cand.uses_cohort_milo,
                    "uses_bulk_labels": cand.uses_bulk_labels,
                    "global_pearson_r": cand.global_pearson_r,
                    "global_spearman_rho": cand.global_spearman_rho,
                    "mean_cohort_r": cand.mean_cohort_r,
                    "positive_cohorts": cand.positive_cohorts,
                    "meta42_r": cand.meta42_r,
                    "perm_null_mean_r": cand.perm_null_mean_r,
                    "perm_null_sd_r": cand.perm_null_sd_r,
                    "perm_z_score": cand.perm_z_score,
                    "perm_p_value": cand.perm_p_value,
                    "state_perm_null_mean_r": cand.state_perm_null_mean_r,
                    "v2_pseudobulk_clr_r": cand.v2_pseudobulk_clr_r,
                    "v3_loto_min_r": cand.v3_loto_min_r,
                    "v3_loto_median_r": cand.v3_loto_median_r,
                    "v3_loto_max_r": cand.v3_loto_max_r,
                    "v4_loto_bulk_auc": cand.v4_loto_bulk_auc,
                }
            )
            print(
                f"[{cand.strategy_id:20s} | {cand.deconv_tool:20s}] "
                f"r={cand.global_pearson_r:+.3f}, rho={cand.global_spearman_rho:+.3f}, "
                f"Pos={cand.positive_cohorts}/9, Meta42_r={cand.meta42_r:+.3f} | "
                f"V1_Z={cand.perm_z_score:+.2f} (p={cand.perm_p_value:.3f}), "
                f"V2_PB_r={cand.v2_pseudobulk_clr_r:+.3f}, "
                f"V3_LOTO=[{cand.v3_loto_min_r:+.2f},{cand.v3_loto_median_r:+.2f},{cand.v3_loto_max_r:+.2f}], "
                f"V4_AUC={cand.v4_loto_bulk_auc:.3f}"
            )

    bench_df = pl.DataFrame(records)
    bench_parquet = config.out_dir / "option_c_alternatives_grid_benchmark.parquet"
    bench_csv = config.out_dir / "option_c_alternatives_grid_benchmark.csv"
    bench_df.write_parquet(bench_parquet)
    bench_df.write_csv(bench_csv)

    folds_path, perms_path, genes_path = export_granular_supplementary_tables(
        ctx,
        all_candidates=all_candidates,
        evaluated_map=evaluated_map,
        out_dir=config.out_dir,
    )

    extra_landmark_caches = (
        ("option_c_unsup_g1_instaprism", ("G1_Unsup_Markers", "D1_InstaPrism_EM")),
        ("option_c_alt1b_loco_milo_ridge", ("G2b_LOCO_Milo_GMRF", "D4_Ridge_LogContrast")),
    )
    for extra_prefix, extra_key in extra_landmark_caches:
        build_figure4_tables_from_candidate(ctx, evaluated_map[extra_key], config.fig4_out_dir, extra_prefix)

    val_chart = build_validation_summary_chart(bench_df)
    val_svg, _ = export_altair_figure(val_chart, config.out_dir / "option_c_validation_summary")
    _tight_crop_svg(val_svg)

    # Render the 4 standalone Option C alternative Figure 4 variants
    f4_cfg = f4i.CombinedConfig()
    coords = pl.read_parquet(f4_cfg.harmony_dir / "pan_cohort_I2_harmony_cohort_substudy_umap_coords.parquet")
    bench_harm = pl.read_parquet(f4_cfg.harmony_dir / "harmony_integration_benchmark.parquet")
    anchors_res = f4i.load_anchor_labels(f4_cfg, ipcm)
    tables_res = g4m.load_cohort_tables(f4_cfg.cohort_dir)
    nhood_res = tables_res.bind(lambda t: g4m.summarize_neighborhoods(t[0]))
    cells_res = tables_res.bind(lambda t: g4m.summarize_cells(t[1]))
    for r_chk in (anchors_res, tables_res, nhood_res, cells_res):
        if isinstance(r_chk, Failure):
            return Failure(str(r_chk.failure()))
    anchors_df = anchors_res.unwrap()
    nhood_df = nhood_res.unwrap()
    cells_df = cells_res.unwrap()

    grid = b4m.COHORT_GRID_ORDER
    labels = b4m.META_STATE_EMPIRICAL_LABELS
    sel_row = bench_harm.filter(pl.col("is_selected_pareto_optimum"))
    ilisi = float(sel_row["cohort_ilisi"].item())
    clisi = float(sel_row["compartment_clisi"].item())
    h = f4_cfg.plot_h

    def _left_anchor() -> Any:
        return (
            alt.Chart(pl.DataFrame({"x": [0.0], "y": [0.0]}).to_pandas())
            .mark_point(opacity=0.0)
            .encode(x=alt.X("x:Q", axis=None), y=alt.Y("y:Q", axis=None))
            .properties(width=0, height=0, view=alt.ViewConfig(stroke=None, fill=None))
        )

    row1 = f4i.independent(
        alt.hconcat(
            _left_anchor(),
            f4i.build_harmony(coords, anchors_df, ilisi, clisi, grid, "a", 360, f4_cfg.umap_h),
            f4i.build_harmony_logfc(coords, anchors_df, "b", 360, f4_cfg.umap_h),
            spacing=18,
        )
    )

    alt_variants = (
        ("option_c_alt1_loco_milo_nusvr", ("G2b_LOCO_Milo_GMRF", "D3_CIBERSORT_NuSVR")),
        ("option_c_alt2_cohort_milo_dwls", ("G2a_Cohort_Milo_GMRF", "D2_DWLS")),
        ("option_c_alt3_state_milo_gmrf_dwls", ("G2d_State_Milo_GMRF", "D2_DWLS")),
        ("option_c_alt4_state_milo_dwls", ("G2d_State_Milo", "D2_DWLS")),
    )

    written_paths: list[Path] = [bench_parquet, bench_csv, folds_path, perms_path, genes_path, val_svg]
    for fig_suffix, key in alt_variants:
        cand = evaluated_map[key]
        u_states_df, u_conc_df = build_figure4_tables_from_candidate(ctx, cand, config.fig4_out_dir, fig_suffix)
        u_selected, u_ordered, u_chosen = f4i.select_meta_states(u_states_df, labels, rank_col="bulk_deconv_beta")
        u_row2 = f4i.independent(
            alt.hconcat(
                _left_anchor(),
                f4i.build_single_scatter(u_states_df, grid, "c", 182, h),
                f4i.build_bars(u_conc_df, "d", 116, h),
                f4i.build_merged_da_panel(nhood_df, cells_df, "e", 109, 109, h),
                spacing=18,
            )
        )
        u_row3 = f4i.independent(
            alt.hconcat(
                _left_anchor(),
                f4i.build_coupled_poles_panel(
                    u_chosen, u_states_df, labels, b4m.BULK_ID_TO_SHORT, u_selected, u_ordered, "f", 276, 276, h
                ),
                spacing=18,
            )
        )
        u_fig = f4i.configure_figure(alt.vconcat(row1, u_row2, u_row3, spacing=14))
        u_supp = f4i.configure_figure(f4i.build_scatter_grid(u_states_df, grid))

        out_svg, _ = export_altair_figure(
            u_fig, config.fig4_out_dir / f"combined_main_figure_{fig_suffix}", max_vector_points=4000
        )
        _tight_crop_svg(out_svg)
        supp_svg, _ = export_altair_figure(
            u_supp, config.fig4_out_dir / f"iatlas_milopy_vs_deconv_scatter_grid_supp_{fig_suffix}"
        )
        _tight_crop_svg(supp_svg)

        # Copy standalone alternative SVGs and PNGs to article/figures/milo-deconv/
        for src in (out_svg, out_svg.with_suffix(".png"), supp_svg, supp_svg.with_suffix(".png")):
            dst = config.article_fig_dir / src.name
            shutil.copy2(src, dst)
            written_paths.append(dst)

    return Success(tuple(written_paths))


def main() -> None:
    """CLI entrypoint for Option C deconvolution & gene-selection exploration."""
    match run(ExploreOptionCConfig()):
        case Success(paths):
            print("\n".join(f"wrote {p}" for p in paths))
        case Failure(err):
            print(f"4k failed: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
