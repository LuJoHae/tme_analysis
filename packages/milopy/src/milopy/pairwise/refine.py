"""Stage 2 covariate-adjusted Firth penalized logistic LRT refinement for candidate pairs.

Fits bias-reduced Jeffreys-prior logistic regression (`firthmodels.FirthLogisticRegression`)
adjusting for donor log-library-size and cohort strata, resolving complete/quasi-complete
separation while vetoing capture-depth-confounded candidates.
"""

from __future__ import annotations

import warnings

from firthmodels import FirthLogisticRegression
import numpy as np
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy.special import chdtri
from scipy.stats import norm

from .types import (
    EnrichmentDirection,
    FloatArray,
    IntArray,
    PresenceRepresentation,
    QuadrantState,
    RefinementConfig,
    adjacent_flips,
)


def _quadrant_signal(
    state: QuadrantState,
    p_u: FloatArray,
    q_u: FloatArray,
    p_v: FloatArray,
    q_v: FloatArray,
) -> FloatArray:
    """Computes donor-level quadrant evidence vector S^{(q)}_{uv} in [0, 1]^J."""
    match state:
        case QuadrantState.BOTH_PRESENT:
            return p_u * p_v
        case QuadrantState.FIRST_ONLY:
            return p_u * q_v
        case QuadrantState.SECOND_ONLY:
            return q_u * p_v
        case QuadrantState.BOTH_MISSING:
            return q_u * q_v


def _build_nuisance_matrix(
    lib_sizes: FloatArray,
    strata_codes: Maybe[IntArray],
    adjust_library_size: bool,
) -> FloatArray:
    """Builds standardized log-library-size and stratum indicator nuisance columns."""
    n_samples = int(lib_sizes.shape[0])
    cols: list[FloatArray] = []

    if adjust_library_size:
        log_lib = np.log(np.maximum(lib_sizes, 1.0))
        sd_lib = float(np.std(log_lib))
        if sd_lib > 1e-8:
            cols.append((log_lib - float(np.mean(log_lib))) / sd_lib)

    match strata_codes:
        case Some(codes):
            unique_c = np.unique(codes)
            if int(unique_c.shape[0]) > 1:
                for c_val in unique_c[1:]:
                    dummy = (codes == c_val).astype(np.float64)
                    if float(np.std(dummy)) > 1e-8:
                        cols.append(dummy)
        case Maybe.empty:
            pass
        case _ as unreachable:
            _ = unreachable

    return (
        np.column_stack(cols)
        if cols
        else np.empty((n_samples, 0), dtype=np.float64)
    )


def _build_full_rank_design(
    primary_col: FloatArray,
    nuisance_cols: FloatArray,
) -> Maybe[FloatArray]:
    """Constructs a full-column-rank design matrix with primary_col at index 0."""
    p_centered = primary_col - float(np.mean(primary_col))
    p_norm = float(np.linalg.norm(p_centered))
    if p_norm <= 1e-7:
        return Nothing

    kept_cols: list[FloatArray] = [primary_col]
    ortho_basis: list[FloatArray] = [p_centered / p_norm]

    for j in range(int(nuisance_cols.shape[1])):
        cand = nuisance_cols[:, j]
        resid = cand - float(np.mean(cand))
        for b_vec in ortho_basis:
            resid = resid - float(np.dot(resid, b_vec)) * b_vec
        r_norm = float(np.linalg.norm(resid))
        if r_norm > 1e-5:
            kept_cols.append(cand)
            ortho_basis.append(resid / r_norm)

    return Some(np.column_stack(kept_cols))


def _fit_firth_feature0(
    primary_col: FloatArray,
    nuisance_cols: FloatArray,
    y_int: IntArray,
    max_iter: int,
) -> Maybe[tuple[float, float]]:
    """Fits Firth logistic regression and returns (coef_0, lrt_pvalue_0) for primary_col."""
    design_maybe = _build_full_rank_design(primary_col, nuisance_cols)
    match design_maybe:
        case Maybe.empty:
            return Nothing
        case Some(x_mat):
            try:
                model = FirthLogisticRegression(
                    backend="numpy",
                    max_iter=max_iter,
                    fit_intercept=True,
                )
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    model.fit(x_mat, y_int)
                    model.lrt(features=0, warm_start=True)
                beta_0 = float(model.coef_[0])
                pval_0 = float(np.clip(model.lrt_pvalues_[0], 1e-300, 1.0))
                if not (np.isfinite(beta_0) and np.isfinite(pval_0)):
                    return Nothing
                return Some((beta_0, pval_0))
            except Exception:
                return Nothing
        case _ as unreachable:
            _ = unreachable
            return Nothing


def _signed_z_from_lrt(beta: float, pval: float) -> float:
    """Converts a two-sided 1-df LRT p-value and effect sign into an equivalent signed Z-score."""
    p_clipped = float(np.clip(pval, 1e-300, 1.0 - 1e-15))
    chi2_val = float(max(chdtri(1.0, p_clipped), 0.0))
    sign = 1.0 if beta >= 0.0 else -1.0
    return sign * float(np.sqrt(chi2_val))


def refine_candidate_pairs_firth(
    candidates_df: pl.DataFrame,
    presence: PresenceRepresentation,
    response_binary: FloatArray,
    strata_codes: Maybe[IntArray] = Nothing,
    config: RefinementConfig = RefinementConfig(),
    fdr_threshold: float = 0.10,
) -> Result[pl.DataFrame, str]:
    """Refines top candidate pair-quadrant rows using covariate-adjusted Firth penalized LRT."""
    if candidates_df.is_empty():
        return Success(
            candidates_df.with_columns(
                pl.Series("log_or_quadrant_firth", [], dtype=pl.Float64),
                pl.Series("pvalue_quadrant_firth", [], dtype=pl.Float64),
                pl.Series("pvalue_iut_firth", [], dtype=pl.Float64),
                pl.Series("is_refined", [], dtype=pl.Boolean),
            )
        )

    n_rows = candidates_df.height
    lor_firth = candidates_df["log_or_quadrant"].to_numpy().astype(np.float64, copy=True)
    pq_firth = candidates_df["pvalue_quadrant"].to_numpy().astype(np.float64, copy=True)
    piut_firth = candidates_df["pvalue_iut"].to_numpy().astype(np.float64, copy=True)
    is_refined = np.zeros(n_rows, dtype=bool)

    if (not config.enabled) or config.max_refine_pairs <= 0:
        return Success(
            candidates_df.with_columns(
                pl.Series("log_or_quadrant_firth", lor_firth),
                pl.Series("pvalue_quadrant_firth", pq_firth),
                pl.Series("pvalue_iut_firth", piut_firth),
                pl.Series("is_refined", is_refined),
            )
        )

    y_int = (response_binary > 0.5).astype(np.int64)
    if int(np.unique(y_int).shape[0]) < 2:
        return Failure("response_binary must contain both responder and non-responder classes")

    nuisance_base = _build_nuisance_matrix(
        lib_sizes=presence.effective_library_sizes,
        strata_codes=strata_codes,
        adjust_library_size=config.adjust_library_size,
    )

    orig_to_local = {
        int(orig_idx): local_idx
        for local_idx, orig_idx in enumerate(presence.active_nhood_indices.tolist())
    }

    n_refine = min(n_rows, config.max_refine_pairs)
    u_col = candidates_df["nhood_u"].to_list()
    v_col = candidates_df["nhood_v"].to_list()
    q_col = candidates_df["quadrant"].to_list()

    p_mat = presence.presence_prob
    q_mat = presence.absence_prob
    m_mat = presence.observability

    for row_idx in range(n_refine):
        u_loc = orig_to_local.get(int(u_col[row_idx]), -1)
        v_loc = orig_to_local.get(int(v_col[row_idx]), -1)
        if u_loc < 0 or v_loc < 0:
            continue

        q_state = QuadrantState(str(q_col[row_idx]))
        flip_v_state, flip_u_state = adjacent_flips(q_state)

        pu, qu, mu = p_mat[u_loc, :], q_mat[u_loc, :], m_mat[u_loc, :]
        pv, qv, mv = p_mat[v_loc, :], q_mat[v_loc, :], m_mat[v_loc, :]
        m_pair = mu * mv

        s_target = _quadrant_signal(q_state, pu, qu, pv, qv)
        s_flip_v = _quadrant_signal(flip_v_state, pu, qu, pv, qv)
        s_flip_u = _quadrant_signal(flip_u_state, pu, qu, pv, qv)

        # 1. Overall quadrant Firth LRT adjusted for log(library_size) and strata
        quad_fit = _fit_firth_feature0(s_target, nuisance_base, y_int, config.max_iter)
        match quad_fit:
            case Maybe.empty:
                continue
            case Some((beta_q, pval_q)):
                lor_firth[row_idx] = beta_q
                pq_firth[row_idx] = pval_q
                is_refined[row_idx] = True

        # 2. Full-sample conditional contrast vs flip_v (absorbing complement m_pair - s_target - s_flip_v)
        comp_v = np.expand_dims(np.maximum(m_pair - s_target - s_flip_v, 0.0), axis=1)
        nuis_v = np.hstack([comp_v, nuisance_base])
        fit_fv = _fit_firth_feature0(s_target, nuis_v, y_int, config.max_iter)

        # 3. Full-sample conditional contrast vs flip_u
        comp_u = np.expand_dims(np.maximum(m_pair - s_target - s_flip_u, 0.0), axis=1)
        nuis_u = np.hstack([comp_u, nuisance_base])
        fit_fu = _fit_firth_feature0(s_target, nuis_u, y_int, config.max_iter)

        match (fit_fv, fit_fu):
            case (Some((b_fv, p_fv)), Some((b_fu, p_fu))):
                z_fv = _signed_z_from_lrt(b_fv, p_fv)
                z_fu = _signed_z_from_lrt(b_fu, p_fu)
                if lor_firth[row_idx] >= 0.0:
                    p_iut_ref = float(2.0 * norm.sf(min(z_fv, z_fu)))
                else:
                    p_iut_ref = float(2.0 * norm.cdf(max(z_fv, z_fu)))
                piut_firth[row_idx] = float(np.clip(p_iut_ref, 0.0, 1.0))
            case _:
                pass

    # Update direction and veto candidates whose signal collapses after covariate adjustment
    old_sig = candidates_df["is_significant"].to_numpy().copy()
    old_dir = candidates_df["direction"].to_numpy().astype(str, copy=True)

    firth_veto = is_refined & (
        (pq_firth > fdr_threshold)
        | (piut_firth > fdr_threshold)
        | (
            (lor_firth >= 0.0)
            != (candidates_df["log_or_quadrant"].to_numpy() >= 0.0)
        )
    )
    updated_sig = old_sig & (~firth_veto)
    updated_dir = np.where(
        is_refined,
        np.where(
            lor_firth >= 0.0,
            EnrichmentDirection.RESPONDER.value,
            EnrichmentDirection.NON_RESPONDER.value,
        ),
        old_dir,
    )

    refined_df = candidates_df.with_columns(
        pl.Series("direction", updated_dir),
        pl.Series("log_or_quadrant_firth", lor_firth),
        pl.Series("pvalue_quadrant_firth", pq_firth),
        pl.Series("pvalue_iut_firth", piut_firth),
        pl.Series("is_refined", is_refined),
        pl.Series("is_significant", updated_sig),
    )
    return Success(refined_df)
