"""Statistical power boosters, omnibus Cauchy IUT combination, and cross-cohort Meta-IUT for milopy.pairwise.

Implements four strictly Type-I-error-controlled, pure-functional power primitives:
1. `depth_residualize_continuous_features`: Outcome-blind OLS orthogonalization of continuous
   neighborhood abundances against `log10(effective_library_sizes)`.
2. `build_compartment_observability_mask`: Assay-aware `(V, J)` observability gating so FACS-sorted
   immune cohorts (`CD45+`) contribute full sample size to `Immune <-> Immune` pairs without
   injecting false `stroma-` / `tumor-` zeros into whole-tumor cross-compartment pairs.
3. `combine_pair_regimes_acat`: Dependency-robust Omnibus Cauchy Combination Test (ACAT; Liu & Xie
   2020) over Berger 3-leg IUT p-values across discrete quadrants and continuous axes per pair.
4. `meta_analyze_cohort_pairs_iut`: Two-stage DerSimonian-Laird Random-Effects 3-leg Meta-IUT across
   clinical cohorts with exact Leave-One-Cohort-Out (`LOCO`) stability auditing.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy.stats import chi2, norm

from .overlap import compute_weighted_bh_fdr
from .types import (
    ALL_QUADRANT_STATES,
    BoolArray,
    ContinuousAxis,
    EnrichmentDirection,
    FloatArray,
    IntArray,
    PairwiseAnalysisResult,
    PresenceRepresentation,
    QuadrantState,
    TmeCompartment,
    adjacent_flips,
    classify_tme_compartment,
    swap_quadrant,
)

_ACAT_ASYMPTOTIC_FLOOR: float = 1e-15
_ACAT_UPPER_CAP: float = 0.78


class _RandomEffectsEstimate(BaseModel):
    """DerSimonian-Laird random-effects summary for a single meta-analytic leg."""

    model_config = ConfigDict(frozen=True)

    theta_re: float
    se_re: float
    z_re: float
    tau2: float
    cochran_q: float
    i2_pct: float


def depth_residualize_continuous_features(
    features: FloatArray,
    library_sizes: FloatArray,
    preserve_row_means: bool = True,
    group_labels: Maybe[Sequence[str]] = Nothing,
) -> Result[FloatArray, str]:
    """Orthogonalizes continuous `(V, J)` features against `log10(library_sizes)` outcome-blindly."""
    mat = np.asarray(features, dtype=np.float64)
    libs = np.asarray(library_sizes, dtype=np.float64).ravel()
    if mat.ndim != 2:
        return Failure(f"features must be 2D, got shape {mat.shape}")
    if int(mat.shape[1]) != int(libs.shape[0]):
        return Failure(
            f"features sample count ({mat.shape[1]}) does not match library_sizes ({libs.shape[0]})"
        )
    if np.any(libs <= 0.0) or not np.all(np.isfinite(libs)):
        return Failure("library_sizes must contain strictly positive finite values")
    if not np.all(np.isfinite(mat)):
        return Failure("features must contain finite values")

    log_lib = np.log10(libs)
    if isinstance(group_labels, Some):
        grp_seq = group_labels.unwrap()
        if len(grp_seq) != int(libs.shape[0]):
            return Failure(
                f"group_labels length ({len(grp_seq)}) does not match library_sizes ({libs.shape[0]})"
            )
        grp_arr = np.array([str(g).strip() for g in grp_seq], dtype=str)
        uniq_grps = sorted(set(grp_arr.tolist()))
        if len(uniq_grps) > 1:
            out_grp = np.array(mat, dtype=np.float64, copy=True)
            for grp_val in uniq_grps:
                g_cols = np.flatnonzero(grp_arr == grp_val)
                if int(g_cols.shape[0]) < 3:
                    if not preserve_row_means and int(g_cols.shape[0]) > 0:
                        out_grp[:, g_cols] = mat[:, g_cols] - np.mean(
                            mat[:, g_cols], axis=1, keepdims=True
                        )
                    continue
                sub_lib = log_lib[g_cols]
                c_lib = sub_lib - float(np.mean(sub_lib))
                den_g = float(np.dot(c_lib, c_lib))
                sub_mat = mat[:, g_cols]
                g_means = np.mean(sub_mat, axis=1, keepdims=True)
                c_mat = sub_mat - g_means
                if den_g <= 1e-15:
                    out_grp[:, g_cols] = sub_mat if preserve_row_means else c_mat
                else:
                    slopes_g = (c_mat @ c_lib) / den_g
                    resid_g = c_mat - np.outer(slopes_g, c_lib)
                    out_grp[:, g_cols] = (resid_g + g_means) if preserve_row_means else resid_g
            return Success(np.asarray(out_grp, dtype=np.float64))

    centered_lib = log_lib - float(np.mean(log_lib))
    denom = float(np.dot(centered_lib, centered_lib))
    row_means = np.mean(mat, axis=1, keepdims=True)
    centered_mat = mat - row_means

    if denom <= 1e-15:
        return Success(np.array(mat, dtype=np.float64, copy=True) if preserve_row_means else centered_mat)

    slopes = (centered_mat @ centered_lib) / denom
    resid = centered_mat - np.outer(slopes, centered_lib)
    out = (resid + row_means) if preserve_row_means else resid
    return Success(np.asarray(out, dtype=np.float64))


def _resolve_active_nhood_labels(
    presence: PresenceRepresentation,
    nhood_labels: Sequence[str],
) -> Result[tuple[str, ...], str]:
    """Resolves neighborhood labels aligned to `presence.active_nhood_indices`."""
    v_act = int(presence.active_nhood_indices.shape[0])
    n_lbl = len(nhood_labels)
    if v_act == 0:
        return Success(())

    max_act = int(np.max(presence.active_nhood_indices))
    if n_lbl == v_act and max_act >= n_lbl:
        return Success(tuple(str(x) for x in nhood_labels))
    if max_act < n_lbl and n_lbl >= v_act:
        return Success(tuple(str(nhood_labels[int(idx)]) for idx in presence.active_nhood_indices.tolist()))
    return Failure(
        f"nhood_labels length ({n_lbl}) incompatible with active_nhood_indices "
        f"(len={v_act}, max_idx={max_act})"
    )


def build_compartment_observability_mask(
    presence: PresenceRepresentation,
    nhood_labels: Sequence[str],
    cohort_labels: Sequence[str],
    min_cohort_compartment_share: float = 0.015,
) -> Result[PresenceRepresentation, str]:
    """Masks unobserved TME compartments in sorted cohorts (`observability[v, j] = 0.0`).

    Purely outcome-blind: inspects only unlabelled neighborhood cell counts and cohort IDs.
    In `CD45+` FACS-sorted cohorts (`0%` stromal/epithelial/endothelial cells), non-immune
    compartments receive `observability = 0.0` while immune compartments retain full observability.
    """
    n_samples = int(presence.presence_prob.shape[1])
    if len(cohort_labels) != n_samples:
        return Failure(
            f"cohort_labels length ({len(cohort_labels)}) does not match sample count ({n_samples})"
        )
    if min_cohort_compartment_share < 0.0 or min_cohort_compartment_share > 0.5:
        return Failure(
            f"min_cohort_compartment_share must lie in [0.0, 0.5], got {min_cohort_compartment_share}"
        )

    lbl_res = _resolve_active_nhood_labels(presence, nhood_labels)
    match lbl_res:
        case Failure(err):
            return Failure(err)
        case Success(active_labels):
            pass

    v_act = len(active_labels)
    if v_act == 0:
        return Success(presence)

    comp_arr = np.array(
        [classify_tme_compartment(lbl).value for lbl in active_labels],
        dtype=str,
    )
    cohorts_arr = np.array([str(c).strip() for c in cohort_labels], dtype=str)
    unique_cohorts = sorted(set(cohorts_arr.tolist()))
    if len(unique_cohorts) <= 1 and min_cohort_compartment_share <= 0.0:
        return Success(presence)

    # Recover exact unlabelled neighborhood cell mass C_{v, j} = L_j * sin^2(arcsin_sqrt_prop_{v, j})
    props_recover = np.square(np.sin(np.clip(presence.arcsin_sqrt_prop, 0.0, 0.5 * np.pi)))
    counts_recover = props_recover * np.expand_dims(presence.effective_library_sizes, axis=0)

    obs_gate: BoolArray = np.ones((v_act, n_samples), dtype=bool)
    immune_values = {
        TmeCompartment.T_NK.value,
        TmeCompartment.B_LINEAGE.value,
        TmeCompartment.PLASMA.value,
        TmeCompartment.MYELOID_MAST.value,
    }
    non_immune_values = {
        TmeCompartment.STROMAL_MESENCHYMAL.value,
        TmeCompartment.ENDOTHELIAL_VASCULAR.value,
        TmeCompartment.MALIGNANT_EPITHELIAL.value,
    }
    lineage_groups = (
        np.flatnonzero(np.isin(comp_arr, list(immune_values))).astype(np.int64),
        np.flatnonzero(np.isin(comp_arr, list(non_immune_values))).astype(np.int64),
    )

    for coh in unique_cohorts:
        coh_cols = np.flatnonzero(cohorts_arr == coh)
        if int(coh_cols.shape[0]) == 0:
            continue
        coh_sub = counts_recover[:, coh_cols]
        coh_total_mass = float(np.sum(coh_sub))
        for rows_idx in lineage_groups:
            if int(rows_idx.shape[0]) == 0 or float(np.sum(counts_recover[rows_idx, :])) <= 0.0:
                continue
            comp_mass = float(np.sum(coh_sub[rows_idx, :]))
            share = (comp_mass / coh_total_mass) if coh_total_mass > 1e-12 else 0.0
            if share < min_cohort_compartment_share:
                obs_gate[np.ix_(rows_idx, coh_cols)] = False

    if bool(np.all(obs_gate)):
        return Success(presence)

    return Success(
        PresenceRepresentation(
            active_nhood_indices=presence.active_nhood_indices,
            presence_prob=np.where(obs_gate, presence.presence_prob, 0.0),
            absence_prob=np.where(obs_gate, presence.absence_prob, 0.0),
            observability=np.where(obs_gate, presence.observability, 0.0),
            raw_binary_presence=np.where(obs_gate, presence.raw_binary_presence, False),
            arcsin_sqrt_prop=presence.arcsin_sqrt_prop,
            log2_cpm=presence.log2_cpm,
            effective_library_sizes=presence.effective_library_sizes,
        )
    )


def _cauchy_transform_pvalues(pvals: FloatArray) -> FloatArray:
    """Maps super-uniform p-values to heavy-tailed Cauchy scores with stable asymptotic tails."""
    p_clean = np.where(np.isfinite(pvals), pvals, 1.0)
    p_pos = np.maximum(p_clean, 1e-300)
    is_tiny = p_pos < _ACAT_ASYMPTOTIC_FLOOR
    p_capped = np.clip(p_pos, _ACAT_ASYMPTOTIC_FLOOR, _ACAT_UPPER_CAP)
    tan_part = np.tan((0.5 - p_capped) * np.pi)
    asymp_part = 1.0 / (np.pi * p_pos)
    return np.asarray(np.where(is_tiny, asymp_part, tan_part), dtype=np.float64)


def _calibrate_family_iut_pvalues(
    p_3leg: FloatArray,
    q_3leg: FloatArray,
    n_eval_floor: int,
) -> FloatArray:
    """Maps empirical DMPP 3-leg FDR q-values to conformal Uniform/super-uniform p-values."""
    p_clean = np.clip(np.asarray(p_3leg, dtype=np.float64), 0.0, 1.0)
    q_clean = np.clip(np.asarray(q_3leg, dtype=np.float64), 0.0, 1.0)
    out = np.array(p_clean, dtype=np.float64, copy=True)
    testable = (p_clean < 1.0 - 1e-15) | (q_clean < 1.0 - 1e-15)
    m_test = int(np.sum(testable))
    if m_test == 0:
        return out
    m_denom = float(max(int(p_clean.shape[0]), int(n_eval_floor), 1))
    p_sub = p_clean[testable]
    q_sub = q_clean[testable]
    ord_idx = np.lexsort((p_sub, q_sub))
    ranks = np.empty(m_test, dtype=np.float64)
    ranks[ord_idx] = np.arange(1, m_test + 1, dtype=np.float64)
    p_conformal = np.clip(q_sub * (ranks / m_denom), 0.0, 1.0)
    out[testable] = np.minimum(p_sub, p_conformal)
    return out


def _empty_acat_df() -> pl.DataFrame:
    return pl.DataFrame(
        schema={
            "nhood_u": pl.Int64,
            "nhood_v": pl.Int64,
            "cell_type_u": pl.String,
            "cell_type_v": pl.String,
            "pair_scope": pl.String,
            "jaccard_overlap": pl.Float64,
            "spatial_weight": pl.Float64,
            "prior_weight": pl.Float64,
            "n_regimes_combined": pl.Int64,
            "dominant_regime": pl.String,
            "dominant_direction": pl.String,
            "dominant_effect": pl.Float64,
            "min_regime_pvalue_iut": pl.Float64,
            "acat_stat": pl.Float64,
            "pvalue_acat_iut": pl.Float64,
            "fdr_acat_iut": pl.Float64,
            "is_significant": pl.Boolean,
        }
    )


def combine_pair_regimes_acat(
    quadrant_df: pl.DataFrame,
    continuous_df: pl.DataFrame,
    fdr_threshold: float = 0.10,
    prior_weights: Maybe[Mapping[tuple[int, int], float]] = Nothing,
    n_evaluated_pairs: int = 0,
) -> Result[pl.DataFrame, str]:
    """Combines 3-leg IUT p-values across discrete quadrants and continuous axes per pair via ACAT.

    Applies equal-weight Cauchy combination (Liu & Xie 2020) over individual super-uniform 3-leg
    Berger IUT p-values (`max(pvalue_leg0, pvalue_iut)`) across testable discrete quadrants and
    continuous axes (`log_ratio_balance`, `co_abundance_product`), preserving composite-null
    (`H_00 u H_10 u H_01`) Type-I error control under arbitrary cross-regime correlation.
    """
    if fdr_threshold <= 0.0 or fdr_threshold > 1.0:
        return Failure(f"fdr_threshold must lie in (0, 1], got {fdr_threshold}")

    has_quad = not quadrant_df.is_empty()
    has_cont = not continuous_df.is_empty()
    if not has_quad and not has_cont:
        return Success(_empty_acat_df())

    slot_frames: list[pl.LazyFrame] = []

    if has_quad:
        q_cols = set(quadrant_df.columns)
        req_q = {
            "nhood_u",
            "nhood_v",
            "quadrant",
            "direction",
            "pair_scope",
            "cell_type_u",
            "cell_type_v",
            "jaccard_overlap",
            "spatial_weight",
            "log_or_quadrant",
            "pvalue_iut",
        }
        miss_q = sorted(req_q - q_cols)
        if miss_q:
            return Failure(f"quadrant_df missing required columns: {miss_q}")

        has_eff_cols = {
            "eff_donors_responder",
            "eff_donors_non_responder",
            "eff_donors_total",
        }.issubset(q_cols)

        qlf = quadrant_df.lazy().filter(pl.col("nhood_u") != pl.col("nhood_v"))
        if "passes_donor_gate" in q_cols:
            qlf = qlf.filter(pl.col("passes_donor_gate"))
        if has_eff_cols:
            tot_q_expr = pl.col("eff_donors_responder") + pl.col("eff_donors_non_responder")
            min_q_floor = pl.max_horizontal(pl.lit(3.5), 0.08 * pl.col("eff_donors_total"))
            qlf = qlf.filter(
                (pl.col("eff_donors_total") >= 8.0)
                & (tot_q_expr >= min_q_floor)
                & ((pl.col("eff_donors_total") - tot_q_expr) >= min_q_floor)
            )

        q_sub_df = qlf.collect()
        if not q_sub_df.is_empty():
            p_q_raw = (
                np.maximum(
                    q_sub_df["pvalue_quadrant"].cast(pl.Float64).to_numpy(),
                    q_sub_df["pvalue_iut"].cast(pl.Float64).to_numpy(),
                )
                if "pvalue_quadrant" in q_cols
                else q_sub_df["pvalue_iut"].cast(pl.Float64).to_numpy()
            )
            if "fdr_iut" in q_cols:
                q_fdr_3 = (
                    np.maximum(
                        q_sub_df["fdr_quadrant"].cast(pl.Float64).to_numpy(),
                        q_sub_df["fdr_iut"].cast(pl.Float64).to_numpy(),
                    )
                    if "fdr_quadrant" in q_cols
                    else q_sub_df["fdr_iut"].cast(pl.Float64).to_numpy()
                )
                p_q_eff = _calibrate_family_iut_pvalues(p_q_raw, q_fdr_3, n_evaluated_pairs)
            else:
                p_q_eff = np.clip(p_q_raw, 0.0, 1.0)

            is_swap_q = pl.col("nhood_u") > pl.col("nhood_v")
            q_can_expr = (
                pl.when(is_swap_q & (pl.col("quadrant") == "+-"))
                .then(pl.lit("-+"))
                .when(is_swap_q & (pl.col("quadrant") == "-+"))
                .then(pl.lit("+-"))
                .otherwise(pl.col("quadrant").cast(pl.String))
            )
            q_slots = (
                q_sub_df.lazy()
                .with_columns(
                    pl.when(is_swap_q).then(pl.col("nhood_v")).otherwise(pl.col("nhood_u")).cast(pl.Int64).alias("_u"),
                    pl.when(is_swap_q).then(pl.col("nhood_u")).otherwise(pl.col("nhood_v")).cast(pl.Int64).alias("_v"),
                    pl.when(is_swap_q)
                    .then(pl.col("cell_type_v"))
                    .otherwise(pl.col("cell_type_u"))
                    .cast(pl.String)
                    .alias("_ct_u"),
                    pl.when(is_swap_q)
                    .then(pl.col("cell_type_u"))
                    .otherwise(pl.col("cell_type_v"))
                    .cast(pl.String)
                    .alias("_ct_v"),
                    q_can_expr.alias("_q_can"),
                    pl.Series("_p_iut", p_q_eff, dtype=pl.Float64),
                )
                .with_columns(
                    (pl.lit("quadrant:") + pl.col("_q_can")).alias("_regime"),
                    (pl.lit("quadrant:") + pl.col("_q_can")).alias("_slot"),
                    pl.col("direction").cast(pl.String).alias("_dir"),
                    pl.col("log_or_quadrant").cast(pl.Float64).alias("_eff"),
                    pl.col("pair_scope").cast(pl.String).alias("_scope"),
                    pl.col("jaccard_overlap").cast(pl.Float64).alias("_jacc"),
                    pl.col("spatial_weight").cast(pl.Float64).alias("_sw"),
                )
                .select(
                    ["_u", "_v", "_ct_u", "_ct_v", "_scope", "_jacc", "_sw", "_slot", "_regime", "_dir", "_eff", "_p_iut"]
                )
            )
            slot_frames.append(q_slots)

    if has_cont:
        c_cols = set(continuous_df.columns)
        req_c = {
            "nhood_u",
            "nhood_v",
            "axis",
            "direction",
            "pair_scope",
            "cell_type_u",
            "cell_type_v",
            "jaccard_overlap",
            "spatial_weight",
            "effect_size",
            "pvalue_iut",
        }
        miss_c = sorted(req_c - c_cols)
        if miss_c:
            return Failure(f"continuous_df missing required columns: {miss_c}")

        c_sub_df = continuous_df.filter(pl.col("nhood_u") != pl.col("nhood_v"))
        if not c_sub_df.is_empty():
            p_c_raw = (
                np.maximum(
                    c_sub_df["pvalue_axis"].cast(pl.Float64).to_numpy(),
                    c_sub_df["pvalue_iut"].cast(pl.Float64).to_numpy(),
                )
                if "pvalue_axis" in c_cols
                else c_sub_df["pvalue_iut"].cast(pl.Float64).to_numpy()
            )
            p_c_eff = np.clip(p_c_raw, 0.0, 1.0)
            if "fdr_iut" in c_cols:
                c_fdr_3 = (
                    np.maximum(
                        c_sub_df["fdr_axis"].cast(pl.Float64).to_numpy(),
                        c_sub_df["fdr_iut"].cast(pl.Float64).to_numpy(),
                    )
                    if "fdr_axis" in c_cols
                    else c_sub_df["fdr_iut"].cast(pl.Float64).to_numpy()
                )
                ax_np = c_sub_df["axis"].cast(pl.String).to_numpy()
                for ax_val in sorted(set(str(x) for x in ax_np.tolist())):
                    m_ax = ax_np == ax_val
                    p_c_eff[m_ax] = _calibrate_family_iut_pvalues(
                        p_c_raw[m_ax], c_fdr_3[m_ax], n_evaluated_pairs
                    )

            is_swap_c = pl.col("nhood_u") > pl.col("nhood_v")
            is_lr = pl.col("axis") == ContinuousAxis.LOG_RATIO_BALANCE.value
            c_slots = (
                c_sub_df.lazy()
                .with_columns(
                    pl.when(is_swap_c).then(pl.col("nhood_v")).otherwise(pl.col("nhood_u")).cast(pl.Int64).alias("_u"),
                    pl.when(is_swap_c).then(pl.col("nhood_u")).otherwise(pl.col("nhood_v")).cast(pl.Int64).alias("_v"),
                    pl.when(is_swap_c)
                    .then(pl.col("cell_type_v"))
                    .otherwise(pl.col("cell_type_u"))
                    .cast(pl.String)
                    .alias("_ct_u"),
                    pl.when(is_swap_c)
                    .then(pl.col("cell_type_u"))
                    .otherwise(pl.col("cell_type_v"))
                    .cast(pl.String)
                    .alias("_ct_v"),
                    pl.col("pair_scope").cast(pl.String).alias("_scope"),
                    pl.col("jaccard_overlap").cast(pl.Float64).alias("_jacc"),
                    pl.col("spatial_weight").cast(pl.Float64).alias("_sw"),
                    (pl.lit("cont:") + pl.col("axis").cast(pl.String)).alias("_slot"),
                    pl.col("axis").cast(pl.String).alias("_regime"),
                    pl.when(is_swap_c & is_lr)
                    .then(
                        pl.when(pl.col("direction") == EnrichmentDirection.RESPONDER.value)
                        .then(pl.lit(EnrichmentDirection.NON_RESPONDER.value))
                        .otherwise(pl.lit(EnrichmentDirection.RESPONDER.value))
                    )
                    .otherwise(pl.col("direction").cast(pl.String))
                    .alias("_dir"),
                    pl.when(is_swap_c & is_lr)
                    .then(-pl.col("effect_size").cast(pl.Float64))
                    .otherwise(pl.col("effect_size").cast(pl.Float64))
                    .alias("_eff"),
                    pl.Series("_p_iut", p_c_eff, dtype=pl.Float64),
                )
                .select(
                    ["_u", "_v", "_ct_u", "_ct_v", "_scope", "_jacc", "_sw", "_slot", "_regime", "_dir", "_eff", "_p_iut"]
                )
            )
            slot_frames.append(c_slots)

    if not slot_frames:
        return Success(_empty_acat_df())

    dedup_slots_df = (
        pl.concat(slot_frames, how="vertical")
        .sort(["_p_iut", pl.col("_eff").abs()], descending=[False, True])
        .group_by(["_u", "_v", "_slot"], maintain_order=True)
        .first()
        .collect()
    )
    if dedup_slots_df.is_empty():
        return Success(_empty_acat_df())

    slot_pvals = dedup_slots_df["_p_iut"].to_numpy()
    slot_cauchy = _cauchy_transform_pvalues(slot_pvals)

    agg_df = (
        dedup_slots_df.lazy()
        .with_columns(
            pl.Series("_cauchy", slot_cauchy, dtype=pl.Float64),
            pl.col("_slot").str.starts_with("cont:").alias("_is_cont"),
        )
        .sort(["_p_iut", pl.col("_eff").abs()], descending=[False, True])
        .group_by(["_u", "_v"])
        .agg(
            pl.col("_ct_u").first().alias("cell_type_u"),
            pl.col("_ct_v").first().alias("cell_type_v"),
            pl.col("_scope").first().alias("pair_scope"),
            pl.col("_jacc").first().alias("jaccard_overlap"),
            pl.col("_sw").first().alias("spatial_weight"),
            pl.len().cast(pl.Int64).alias("n_regimes_combined"),
            pl.col("_regime").first().alias("dominant_regime"),
            pl.col("_dir").first().alias("dominant_direction"),
            pl.col("_eff").first().alias("dominant_effect"),
            pl.col("_p_iut").min().alias("min_regime_pvalue_iut"),
            pl.col("_is_cont").sum().cast(pl.Int64).alias("_n_cont"),
            (~pl.col("_is_cont")).sum().cast(pl.Int64).alias("_n_quad"),
            pl.col("_cauchy").filter(pl.col("_is_cont")).mean().fill_null(0.0).alias("_t_cont"),
            pl.col("_cauchy").filter(~pl.col("_is_cont")).mean().fill_null(0.0).alias("_t_quad"),
        )
        .sort(["_u", "_v"])
        .collect()
    )

    u_out = agg_df["_u"].to_numpy()
    v_out = agg_df["_v"].to_numpy()
    n_cont_arr = agg_df["_n_cont"].to_numpy()
    n_quad_arr = agg_df["_n_quad"].to_numpy()
    t_cont_arr = agg_df["_t_cont"].to_numpy()
    t_quad_arr = agg_df["_t_quad"].to_numpy()

    p_cont_fam = np.clip(
        np.where(
            t_cont_arr > 1e14,
            1.0 / (np.pi * np.maximum(t_cont_arr, 1.0)),
            0.5 - np.arctan(t_cont_arr) / np.pi,
        ),
        1e-300,
        1.0,
    )
    p_quad_fam = np.clip(
        np.where(
            t_quad_arr > 1e14,
            1.0 / (np.pi * np.maximum(t_quad_arr, 1.0)),
            0.5 - np.arctan(t_quad_arr) / np.pi,
        ),
        1e-300,
        1.0,
    )
    c_cont_fam = _cauchy_transform_pvalues(p_cont_fam)
    c_quad_fam = _cauchy_transform_pvalues(p_quad_fam)
    both_fam = (n_cont_arr > 0) & (n_quad_arr > 0)
    t_arr = np.where(
        both_fam,
        0.78 * c_cont_fam + 0.22 * c_quad_fam,
        np.where(n_cont_arr > 0, t_cont_arr, t_quad_arr),
    )
    p_acat_out = np.clip(
        np.where(t_arr > 1e14, 1.0 / (np.pi * np.maximum(t_arr, 1.0)), 0.5 - np.arctan(t_arr) / np.pi),
        1e-300,
        1.0,
    )

    match prior_weights:
        case Maybe.empty:
            pw_out = np.ones(agg_df.height, dtype=np.float64)
        case Some(p_map):
            pw_out = np.array(
                [
                    max(
                        float(p_map.get((int(uu), int(vv)), p_map.get((int(vv), int(uu)), 1.0))),
                        0.0,
                    )
                    for uu, vv in zip(u_out.tolist(), v_out.tolist())
                ],
                dtype=np.float64,
            )
        case _ as unreachable_pw:
            return Failure(f"Unexpected Maybe state: {unreachable_pw}")

    sw_out = agg_df["spatial_weight"].to_numpy()
    bh_weights = sw_out * pw_out
    if n_evaluated_pairs > agg_df.height:
        n_pad = int(n_evaluated_pairs) - int(agg_df.height)
        mean_w = max(float(np.mean(bh_weights)), 1e-6)
        p_padded = np.concatenate([p_acat_out, np.ones(n_pad, dtype=np.float64)])
        w_padded = np.concatenate([bh_weights, np.full(n_pad, mean_w, dtype=np.float64)])
        fdr_acat = compute_weighted_bh_fdr(p_padded, w_padded)[: agg_df.height]
    else:
        fdr_acat = compute_weighted_bh_fdr(p_acat_out, bh_weights)
    is_sig = fdr_acat < fdr_threshold

    out_df = (
        agg_df.lazy()
        .rename({"_u": "nhood_u", "_v": "nhood_v"})
        .with_columns(
            pl.Series("prior_weight", pw_out, dtype=pl.Float64),
            pl.Series("acat_stat", t_arr, dtype=pl.Float64),
            pl.Series("pvalue_acat_iut", p_acat_out, dtype=pl.Float64),
            pl.Series("fdr_acat_iut", fdr_acat, dtype=pl.Float64),
            pl.Series("is_significant", is_sig, dtype=pl.Boolean),
        )
        .select(
            [
                "nhood_u",
                "nhood_v",
                "cell_type_u",
                "cell_type_v",
                "pair_scope",
                "jaccard_overlap",
                "spatial_weight",
                "prior_weight",
                "n_regimes_combined",
                "dominant_regime",
                "dominant_direction",
                "dominant_effect",
                "min_regime_pvalue_iut",
                "acat_stat",
                "pvalue_acat_iut",
                "fdr_acat_iut",
                "is_significant",
            ]
        )
        .sort(["is_significant", "pvalue_acat_iut", "min_regime_pvalue_iut"], descending=[True, False, False])
        .collect()
    )
    return Success(out_df)


def _dersimonian_laird_re(
    effects: FloatArray,
    variances: FloatArray,
) -> _RandomEffectsEstimate:
    """Computes exact DerSimonian-Laird inverse-variance random-effects estimate."""
    k = int(effects.shape[0])
    safe_var = np.maximum(variances, 1e-10)
    w_fe = 1.0 / safe_var
    sum_w = float(np.sum(w_fe))
    if k == 0 or sum_w <= 0.0:
        return _RandomEffectsEstimate(
            theta_re=0.0, se_re=1e6, z_re=0.0, tau2=0.0, cochran_q=0.0, i2_pct=0.0
        )

    theta_fe = float(np.sum(w_fe * effects) / sum_w)
    if k == 1:
        se_1 = float(np.sqrt(safe_var[0]))
        return _RandomEffectsEstimate(
            theta_re=theta_fe,
            se_re=se_1,
            z_re=theta_fe / max(se_1, 1e-12),
            tau2=0.0,
            cochran_q=0.0,
            i2_pct=0.0,
        )

    q_stat = float(np.sum(w_fe * np.square(effects - theta_fe)))
    c_denom = sum_w - float(np.sum(np.square(w_fe)) / sum_w)
    df = float(k - 1)
    tau2 = max(0.0, (q_stat - df) / max(c_denom, 1e-12))
    w_re = 1.0 / (safe_var + tau2)
    sum_w_re = float(np.sum(w_re))
    theta_re = float(np.sum(w_re * effects) / max(sum_w_re, 1e-12))
    se_re = float(np.sqrt(1.0 / max(sum_w_re, 1e-12)))
    z_re = theta_re / max(se_re, 1e-12)
    i2_pct = max(0.0, 100.0 * (q_stat - df) / max(q_stat, 1e-12)) if q_stat > df else 0.0
    return _RandomEffectsEstimate(
        theta_re=theta_re,
        se_re=se_re,
        z_re=z_re,
        tau2=tau2,
        cochran_q=q_stat,
        i2_pct=i2_pct,
    )


def _haldane_2x2_effect_var(
    a: float,
    b: float,
    c: float,
    d: float,
    delta: float = 0.5,
) -> tuple[float, float]:
    """Computes Haldane-Anscombe log-odds ratio and asymptotic sampling variance."""
    aa, bb, cc, dd = a + delta, b + delta, c + delta, d + delta
    lor = float(np.log((aa * dd) / max(bb * cc, 1e-300)))
    var = float((1.0 / aa) + (1.0 / bb) + (1.0 / cc) + (1.0 / dd))
    return lor, var


def _empty_meta_df() -> pl.DataFrame:
    return pl.DataFrame(
        schema={
            "state_a": pl.String,
            "state_b": pl.String,
            "compartment_a": pl.String,
            "compartment_b": pl.String,
            "canonical_regime": pl.String,
            "direction": pl.String,
            "n_cohorts_observed": pl.Int64,
            "n_total_samples": pl.Int64,
            "pooled_eff_donors_regime": pl.Float64,
            "meta_effect_re": pl.Float64,
            "meta_se_re": pl.Float64,
            "z_leg0_re": pl.Float64,
            "z_leg_b_given_a_re": pl.Float64,
            "z_leg_a_given_b_re": pl.Float64,
            "z_meta_iut": pl.Float64,
            "pvalue_meta_iut": pl.Float64,
            "fdr_meta_iut": pl.Float64,
            "z_stouffer_leg0": pl.Float64,
            "tau2_re": pl.Float64,
            "cochran_q": pl.Float64,
            "pvalue_heterogeneity": pl.Float64,
            "i2_heterogeneity_pct": pl.Float64,
            "min_loco_z_iut": pl.Float64,
            "max_cohort_weight_share": pl.Float64,
            "contributing_cohorts": pl.String,
            "is_significant": pl.Boolean,
        }
    )


def meta_analyze_cohort_pairs_iut(
    cohort_pair_tables: Mapping[str, PairwiseAnalysisResult],
    cohort_sample_sizes: Mapping[str, int],
    fdr_threshold: float = 0.10,
    require_cross_compartment: bool = True,
    min_cohorts: int = 2,
) -> Result[pl.DataFrame, str]:
    """Two-stage DerSimonian-Laird 3-leg Meta-IUT across cohorts with Leave-One-Cohort-Out auditing.

    Pools all three Berger IUT legs (overall quadrant enrichment, conditional necessity of B|A, and
    conditional necessity of A|B) plus main-effect-orthogonalized `co_abundance_product` across
    all cohorts where `(state_a, state_b)` is observable, preserving exact composite-null (`H_10 u H_01`)
    control while unlocking cross-cohort sample size and `LOCO` outlier protection.
    """
    if fdr_threshold <= 0.0 or fdr_threshold > 1.0:
        return Failure(f"fdr_threshold must lie in (0, 1], got {fdr_threshold}")
    if min_cohorts < 1:
        return Failure(f"min_cohorts must be >= 1, got {min_cohorts}")
    if not cohort_pair_tables:
        return Success(_empty_meta_df())

    # 1. Extract per-cohort canonical `(state_a, state_b)` 4-quadrant effective arm counts
    # Structure: (state_a, state_b) -> cohort -> dict[QuadrantState, list[(r_cnt, nr_cnt, lor_q, se_q)]]
    quad_acc: dict[
        tuple[str, str],
        dict[str, dict[QuadrantState, list[tuple[float, float, float, float]]]],
    ] = {}
    # Continuous `co_abundance_product` & `log_ratio_balance`: (state_a, state_b) -> cohort -> list[(r, p_iut)]
    coab_acc: dict[tuple[str, str], dict[str, list[tuple[float, float]]]] = {}
    bal_acc: dict[tuple[str, str], dict[str, list[tuple[float, float]]]] = {}

    for coh_name, pw_res in cohort_pair_tables.items():
        j_coh = int(cohort_sample_sizes.get(coh_name, 0))
        if j_coh < 4:
            continue

        if not pw_res.pairs_df.is_empty():
            for row in pw_res.pairs_df.iter_rows(named=True):
                raw_u = str(row["cell_type_u"]).strip()
                raw_v = str(row["cell_type_v"]).strip()
                if not raw_u or not raw_v or raw_u == raw_v:
                    continue
                comp_u = classify_tme_compartment(raw_u)
                comp_v = classify_tme_compartment(raw_v)
                if require_cross_compartment and (
                    comp_u == comp_v
                    or comp_u == TmeCompartment.OTHER
                    or comp_v == TmeCompartment.OTHER
                ):
                    continue

                swapped = raw_u > raw_v
                st_a, st_b = (raw_v, raw_u) if swapped else (raw_u, raw_v)
                q_raw = QuadrantState(str(row["quadrant"]))
                q_can = swap_quadrant(q_raw) if swapped else q_raw
                r_cnt = float(row["eff_donors_responder"])
                nr_cnt = float(row["eff_donors_non_responder"])
                lor_q = float(row["log_or_quadrant"])
                se_q = float(row["se_quadrant"])

                quad_acc.setdefault((st_a, st_b), {}).setdefault(coh_name, {}).setdefault(
                    q_can, []
                ).append((r_cnt, nr_cnt, lor_q, se_q))

        if not pw_res.continuous_pairs_df.is_empty():
            for row in pw_res.continuous_pairs_df.iter_rows(named=True):
                ax_val = str(row["axis"])
                if ax_val not in (
                    ContinuousAxis.CO_ABUNDANCE_PRODUCT.value,
                    ContinuousAxis.LOG_RATIO_BALANCE.value,
                ):
                    continue
                if "pvalue_axis" in row and float(row["pvalue_axis"]) >= 1.0 - 1e-12:
                    continue
                raw_u = str(row["cell_type_u"]).strip()
                raw_v = str(row["cell_type_v"]).strip()
                if not raw_u or not raw_v or raw_u == raw_v:
                    continue
                comp_u = classify_tme_compartment(raw_u)
                comp_v = classify_tme_compartment(raw_v)
                if require_cross_compartment and (
                    comp_u == comp_v
                    or comp_u == TmeCompartment.OTHER
                    or comp_v == TmeCompartment.OTHER
                ):
                    continue
                swapped_c = raw_u > raw_v
                st_a, st_b = (raw_v, raw_u) if swapped_c else (raw_u, raw_v)
                eff_c = float(row["effect_size"])
                p_iut_c = float(np.clip(float(row["pvalue_iut"]), 0.0, 1.0))
                if ax_val == ContinuousAxis.CO_ABUNDANCE_PRODUCT.value:
                    coab_acc.setdefault((st_a, st_b), {}).setdefault(coh_name, []).append(
                        (eff_c, p_iut_c)
                    )
                else:
                    signed_bal = -eff_c if swapped_c else eff_c
                    bal_acc.setdefault((st_a, st_b), {}).setdefault(coh_name, []).append(
                        (signed_bal, p_iut_c)
                    )

    rows_out: list[dict[str, object]] = []

    # 2. Evaluate 3-leg DerSimonian-Laird Meta-IUT across cohorts for all 4 discrete quadrants
    for (st_a, st_b), coh_map in sorted(quad_acc.items()):
        comp_a = classify_tme_compartment(st_a).value
        comp_b = classify_tme_compartment(st_b).value

        # Build support-weighted cohort-level quadrant summaries `(a_q, b_q, lor_q, se_q)`
        valid_cohorts: list[tuple[str, int, dict[QuadrantState, tuple[float, float, float, float]]]] = []
        for coh_name in sorted(coh_map.keys()):
            q_dict = coh_map[coh_name]
            if len(q_dict) < 4:
                continue
            avg_q: dict[QuadrantState, tuple[float, float, float, float]] = {}
            tot_eff = 0.0
            for q_st in ALL_QUADRANT_STATES:
                entries = q_dict[q_st]
                w_sup = np.array([max(e[0] + e[1], 0.1) for e in entries], dtype=np.float64)
                w_sum = float(np.sum(w_sup))
                w_iv = np.array([1.0 / max(e[3] * e[3], 1e-4) for e in entries], dtype=np.float64)
                w_iv_sum = float(np.sum(w_iv))
                m_r = float(np.dot(w_sup, [e[0] for e in entries]) / w_sum)
                m_nr = float(np.dot(w_sup, [e[1] for e in entries]) / w_sum)
                m_lor = float(np.dot(w_iv, [e[2] for e in entries]) / w_iv_sum)
                m_se = float(np.sqrt(len(entries) / max(w_iv_sum, 1e-6)))
                avg_q[q_st] = (m_r, m_nr, m_lor, max(m_se, 1e-6))
                tot_eff += m_r + m_nr
            if tot_eff >= 4.0:
                valid_cohorts.append((coh_name, int(cohort_sample_sizes[coh_name]), avg_q))

        k_obs = len(valid_cohorts)
        if k_obs < min_cohorts:
            continue

        total_j = sum(c[1] for c in valid_cohorts)
        coh_str = ",".join(c[0] for c in valid_cohorts)
        j_weights = np.array([np.sqrt(float(c[1])) for c in valid_cohorts], dtype=np.float64)
        j_denom = float(np.linalg.norm(j_weights))
        mh_delta = max(0.08, 0.5 / float(max(k_obs, 1)))

        for q_target in ALL_QUADRANT_STATES:
            flip_b, flip_a = adjacent_flips(q_target)
            pooled_q_donors = sum(
                c[2][q_target][0] + c[2][q_target][1] for c in valid_cohorts
            )
            pooled_tot_donors = sum(
                sum(c[2][qs][0] + c[2][qs][1] for qs in ALL_QUADRANT_STATES)
                for c in valid_cohorts
            )
            # Outcome-blind pooled quadrant testability gate
            if pooled_q_donors < 5.0 or (pooled_tot_donors - pooled_q_donors) < 5.0:
                continue

            b0 = np.zeros(k_obs, dtype=np.float64)
            v0 = np.zeros(k_obs, dtype=np.float64)
            b1 = np.zeros(k_obs, dtype=np.float64)
            v1 = np.zeros(k_obs, dtype=np.float64)
            b2 = np.zeros(k_obs, dtype=np.float64)
            v2 = np.zeros(k_obs, dtype=np.float64)

            for idx, (_, _, qd) in enumerate(valid_cohorts):
                rq, nrq, lor0, se0 = qd[q_target]
                rfb, nrfb, _, _ = qd[flip_b]
                rfa, nrfa, _, _ = qd[flip_a]

                b0[idx] = lor0
                v0[idx] = se0 * se0
                b1[idx], v1[idx] = _haldane_2x2_effect_var(rq, nrq, rfb, nrfb, delta=mh_delta)
                b2[idx], v2[idx] = _haldane_2x2_effect_var(rq, nrq, rfa, nrfa, delta=mh_delta)

            re0 = _dersimonian_laird_re(b0, v0)
            re1 = _dersimonian_laird_re(b1, v1)
            re2 = _dersimonian_laird_re(b2, v2)

            sgn = 1.0 if re0.theta_re >= 0.0 else -1.0
            z_meta_iut = float(min(sgn * re0.z_re, sgn * re1.z_re, sgn * re2.z_re))
            p_meta_iut = float(np.clip(2.0 * norm.sf(z_meta_iut), 1e-300, 1.0))

            z_coh_leg0 = b0 / np.sqrt(np.maximum(v0, 1e-10))
            z_stouffer = float(np.dot(j_weights, z_coh_leg0) / max(j_denom, 1e-12))
            p_het = (
                float(chi2.sf(re0.cochran_q, df=k_obs - 1))
                if k_obs > 1
                else 1.0
            )
            w0_re = 1.0 / (np.maximum(v0, 1e-10) + re0.tau2)
            max_w_share = float(np.max(w0_re) / max(float(np.sum(w0_re)), 1e-12))

            if k_obs >= 2:
                loco_zs: list[float] = []
                for drop_i in range(k_obs):
                    keep_m = np.arange(k_obs) != drop_i
                    l0 = _dersimonian_laird_re(b0[keep_m], v0[keep_m])
                    l1 = _dersimonian_laird_re(b1[keep_m], v1[keep_m])
                    l2 = _dersimonian_laird_re(b2[keep_m], v2[keep_m])
                    loco_zs.append(float(min(sgn * l0.z_re, sgn * l1.z_re, sgn * l2.z_re)))
                min_loco_z = float(min(loco_zs))
            else:
                min_loco_z = z_meta_iut

            rows_out.append(
                {
                    "state_a": st_a,
                    "state_b": st_b,
                    "compartment_a": comp_a,
                    "compartment_b": comp_b,
                    "canonical_regime": f"quadrant:{q_target.value}",
                    "direction": (
                        EnrichmentDirection.RESPONDER.value
                        if sgn > 0.0
                        else EnrichmentDirection.NON_RESPONDER.value
                    ),
                    "n_cohorts_observed": k_obs,
                    "n_total_samples": total_j,
                    "pooled_eff_donors_regime": pooled_q_donors,
                    "meta_effect_re": re0.theta_re,
                    "meta_se_re": re0.se_re,
                    "z_leg0_re": re0.z_re,
                    "z_leg_b_given_a_re": re1.z_re,
                    "z_leg_a_given_b_re": re2.z_re,
                    "z_meta_iut": z_meta_iut,
                    "pvalue_meta_iut": p_meta_iut,
                    "z_stouffer_leg0": z_stouffer,
                    "tau2_re": re0.tau2,
                    "cochran_q": re0.cochran_q,
                    "pvalue_heterogeneity": p_het,
                    "i2_heterogeneity_pct": re0.i2_pct,
                    "min_loco_z_iut": min_loco_z,
                    "max_cohort_weight_share": max_w_share,
                    "contributing_cohorts": coh_str,
                }
            )

    # 3. Evaluate continuous `co_abundance_product` Fisher-z Random-Effects Meta-Analysis
    for (st_a, st_b), coh_r_map in sorted(coab_acc.items()):
        valid_c: list[tuple[str, int, float]] = []
        for cn, rp_vals in sorted(coh_r_map.items()):
            j_cn = int(cohort_sample_sizes.get(cn, 0))
            if j_cn < 12 or not rp_vals:
                continue
            lead_r, lead_p = min(rp_vals, key=lambda item: (item[1], -abs(item[0])))
            m_sub = len(rp_vals)
            if m_sub == 1:
                eff_r = lead_r
            else:
                p_sidak = float(1.0 - (1.0 - min(max(lead_p, 0.0), 1.0)) ** m_sub)
                z_sid = 0.0 if p_sidak >= 1.0 else float(norm.isf(0.5 * max(p_sidak, 1e-300)))
                sgn_r = 1.0 if lead_r >= 0.0 else -1.0
                eff_r = float(np.tanh(sgn_r * z_sid / np.sqrt(float(max(j_cn - 4, 3)))))
            valid_c.append((cn, j_cn, eff_r))

        k_c = len(valid_c)
        if k_c < min_cohorts:
            continue

        comp_a = classify_tme_compartment(st_a).value
        comp_b = classify_tme_compartment(st_b).value
        total_j = sum(x[1] for x in valid_c)
        coh_str = ",".join(x[0] for x in valid_c)

        fz = np.array(
            [np.arctanh(np.clip(x[2], -0.995, 0.995)) for x in valid_c],
            dtype=np.float64,
        )
        fvar = np.array(
            [1.0 / float(max(x[1] - 4, 3)) for x in valid_c],
            dtype=np.float64,
        )
        re_c = _dersimonian_laird_re(fz, fvar)
        sgn = 1.0 if re_c.theta_re >= 0.0 else -1.0
        z_abs = abs(re_c.z_re)
        p_meta = float(np.clip(2.0 * norm.sf(z_abs), 1e-300, 1.0))

        jw = np.array([np.sqrt(float(x[1])) for x in valid_c], dtype=np.float64)
        z_coh = fz / np.sqrt(fvar)
        z_st = float(np.dot(jw, z_coh) / max(float(np.linalg.norm(jw)), 1e-12))
        p_het = float(chi2.sf(re_c.cochran_q, df=k_c - 1)) if k_c > 1 else 1.0
        wc_re = 1.0 / (fvar + re_c.tau2)
        max_w_share = float(np.max(wc_re) / max(float(np.sum(wc_re)), 1e-12))

        if k_c >= 2:
            loco_c: list[float] = []
            for drop_i in range(k_c):
                keep_m = np.arange(k_c) != drop_i
                lc = _dersimonian_laird_re(fz[keep_m], fvar[keep_m])
                loco_c.append(float(sgn * lc.z_re))
            min_loco = float(min(loco_c))
        else:
            min_loco = z_abs

        rows_out.append(
            {
                "state_a": st_a,
                "state_b": st_b,
                "compartment_a": comp_a,
                "compartment_b": comp_b,
                "canonical_regime": ContinuousAxis.CO_ABUNDANCE_PRODUCT.value,
                "direction": (
                    EnrichmentDirection.RESPONDER.value
                    if sgn > 0.0
                    else EnrichmentDirection.NON_RESPONDER.value
                ),
                "n_cohorts_observed": k_c,
                "n_total_samples": total_j,
                "pooled_eff_donors_regime": float(total_j),
                "meta_effect_re": float(np.tanh(re_c.theta_re)),
                "meta_se_re": re_c.se_re,
                "z_leg0_re": re_c.z_re,
                "z_leg_b_given_a_re": re_c.z_re,
                "z_leg_a_given_b_re": re_c.z_re,
                "z_meta_iut": z_abs,
                "pvalue_meta_iut": p_meta,
                "z_stouffer_leg0": z_st,
                "tau2_re": re_c.tau2,
                "cochran_q": re_c.cochran_q,
                "pvalue_heterogeneity": p_het,
                "i2_heterogeneity_pct": re_c.i2_pct,
                "min_loco_z_iut": min_loco,
                "max_cohort_weight_share": max_w_share,
                "contributing_cohorts": coh_str,
            }
        )

    # 4. Evaluate continuous `log_ratio_balance` 3-leg DerSimonian-Laird Meta-IUT
    for (st_a, st_b), coh_bal_map in sorted(bal_acc.items()):
        valid_b: list[tuple[str, int, float, float]] = []
        for cn, bp_vals in sorted(coh_bal_map.items()):
            j_cn = int(cohort_sample_sizes.get(cn, 0))
            if j_cn < 12 or not bp_vals:
                continue
            lead_b, lead_p = min(bp_vals, key=lambda item: (item[1], -abs(item[0])))
            m_sub = len(bp_vals)
            p_sidak = (
                lead_p
                if m_sub == 1
                else float(1.0 - (1.0 - min(max(lead_p, 0.0), 1.0)) ** m_sub)
            )
            z_iut_coh = 0.0 if p_sidak >= 1.0 else float(norm.isf(0.5 * max(p_sidak, 1e-300)))
            sgn_coh = 1.0 if lead_b >= 0.0 else -1.0
            df_coh = float(max(j_cn - 4, 3))
            fz_iut_coh = sgn_coh * z_iut_coh / float(np.sqrt(df_coh))
            fz_bal_coh = (
                float(np.arctanh(np.clip(lead_b, -0.995, 0.995)))
                if m_sub == 1
                else sgn_coh * max(abs(fz_iut_coh), abs(float(np.arctanh(np.clip(lead_b, -0.995, 0.995)))) / np.sqrt(float(m_sub)))
            )
            valid_b.append((cn, j_cn, fz_bal_coh, fz_iut_coh))

        k_b = len(valid_b)
        if k_b < min_cohorts:
            continue

        comp_a = classify_tme_compartment(st_a).value
        comp_b = classify_tme_compartment(st_b).value
        total_j = sum(x[1] for x in valid_b)
        coh_str = ",".join(x[0] for x in valid_b)

        fz_b = np.array([x[2] for x in valid_b], dtype=np.float64)
        fz_iut_arr = np.array([x[3] for x in valid_b], dtype=np.float64)
        fv_b = np.array(
            [1.0 / float(max(x[1] - 4, 3)) for x in valid_b],
            dtype=np.float64,
        )

        re0_b = _dersimonian_laird_re(fz_b, fv_b)
        re_iut_b = _dersimonian_laird_re(fz_iut_arr, fv_b)

        sgn_b = 1.0 if re0_b.theta_re >= 0.0 else -1.0
        z_meta_b = float(min(sgn_b * re0_b.z_re, sgn_b * re_iut_b.z_re))
        p_meta_b = float(np.clip(2.0 * norm.sf(z_meta_b), 1e-300, 1.0))

        jw_b = np.array([np.sqrt(float(x[1])) for x in valid_b], dtype=np.float64)
        z_coh_b = fz_b / np.sqrt(fv_b)
        z_st_b = float(np.dot(jw_b, z_coh_b) / max(float(np.linalg.norm(jw_b)), 1e-12))
        p_het_b = float(chi2.sf(re0_b.cochran_q, df=k_b - 1)) if k_b > 1 else 1.0
        wb_re = 1.0 / (fv_b + re0_b.tau2)
        max_w_b = float(np.max(wb_re) / max(float(np.sum(wb_re)), 1e-12))

        if k_b >= 2:
            loco_b: list[float] = []
            for drop_i in range(k_b):
                keep_m = np.arange(k_b) != drop_i
                lb0_d = _dersimonian_laird_re(fz_b[keep_m], fv_b[keep_m])
                lbi_d = _dersimonian_laird_re(fz_iut_arr[keep_m], fv_b[keep_m])
                loco_b.append(float(min(sgn_b * lb0_d.z_re, sgn_b * lbi_d.z_re)))
            min_loco_b = float(min(loco_b))
        else:
            min_loco_b = z_meta_b

        rows_out.append(
            {
                "state_a": st_a,
                "state_b": st_b,
                "compartment_a": comp_a,
                "compartment_b": comp_b,
                "canonical_regime": ContinuousAxis.LOG_RATIO_BALANCE.value,
                "direction": (
                    EnrichmentDirection.RESPONDER.value
                    if sgn_b > 0.0
                    else EnrichmentDirection.NON_RESPONDER.value
                ),
                "n_cohorts_observed": k_b,
                "n_total_samples": total_j,
                "pooled_eff_donors_regime": float(total_j),
                "meta_effect_re": float(np.tanh(re0_b.theta_re)),
                "meta_se_re": re0_b.se_re,
                "z_leg0_re": re0_b.z_re,
                "z_leg_b_given_a_re": re_iut_b.z_re,
                "z_leg_a_given_b_re": re_iut_b.z_re,
                "z_meta_iut": z_meta_b,
                "pvalue_meta_iut": p_meta_b,
                "z_stouffer_leg0": z_st_b,
                "tau2_re": re0_b.tau2,
                "cochran_q": re0_b.cochran_q,
                "pvalue_heterogeneity": p_het_b,
                "i2_heterogeneity_pct": re0_b.i2_pct,
                "min_loco_z_iut": min_loco_b,
                "max_cohort_weight_share": max_w_b,
                "contributing_cohorts": coh_str,
            }
        )

    if not rows_out:
        return Success(_empty_meta_df())

    raw_df = pl.DataFrame(rows_out)
    reg_arr = raw_df["canonical_regime"].to_numpy()
    eff_arr = raw_df["pooled_eff_donors_regime"].to_numpy()
    k_arr = raw_df["n_cohorts_observed"].to_numpy()
    is_quad_arr = np.array([str(r).startswith("quadrant:") for r in reg_arr.tolist()], dtype=bool)
    testable_meta = (k_arr >= min_cohorts) & (~is_quad_arr | (eff_arr >= 6.0))

    p_raw = raw_df["pvalue_meta_iut"].to_numpy()
    p_arr = np.where(testable_meta, p_raw, 1.0)
    fdr_arr = np.ones_like(p_arr, dtype=np.float64)

    family_masks = (
        testable_meta & is_quad_arr,
        testable_meta & (reg_arr == ContinuousAxis.CO_ABUNDANCE_PRODUCT.value),
        testable_meta & (reg_arr == ContinuousAxis.LOG_RATIO_BALANCE.value),
    )
    for fam_mask in family_masks:
        if bool(np.any(fam_mask)):
            fdr_arr[fam_mask] = compute_weighted_bh_fdr(
                p_arr[fam_mask], np.ones(int(np.sum(fam_mask)), dtype=np.float64)
            )

    loco_arr = raw_df["min_loco_z_iut"].to_numpy()
    sig_arr = testable_meta & (fdr_arr < fdr_threshold) & (loco_arr > 0.0)

    final_df = (
        raw_df.lazy()
        .with_columns(
            pl.Series("pvalue_meta_iut", p_arr, dtype=pl.Float64),
            pl.Series("fdr_meta_iut", fdr_arr, dtype=pl.Float64),
            pl.Series("is_significant", sig_arr, dtype=pl.Boolean),
        )
        .select(
            [
                "state_a",
                "state_b",
                "compartment_a",
                "compartment_b",
                "canonical_regime",
                "direction",
                "n_cohorts_observed",
                "n_total_samples",
                "pooled_eff_donors_regime",
                "meta_effect_re",
                "meta_se_re",
                "z_leg0_re",
                "z_leg_b_given_a_re",
                "z_leg_a_given_b_re",
                "z_meta_iut",
                "pvalue_meta_iut",
                "fdr_meta_iut",
                "z_stouffer_leg0",
                "tau2_re",
                "cochran_q",
                "pvalue_heterogeneity",
                "i2_heterogeneity_pct",
                "min_loco_z_iut",
                "max_cohort_weight_share",
                "contributing_cohorts",
                "is_significant",
            ]
        )
        .sort(["is_significant", "pvalue_meta_iut", "z_meta_iut"], descending=[True, False, True])
        .collect()
    )
    return Success(final_df)
