"""Depth-calibrated neighborhood presence/absence estimation and compositional transformations.

Transforms raw neighborhood count matrices C (V x J) and donor library sizes L (J,) into:
1. Calibrated presence probabilities P in [0, 1] and detection-weighted absence probabilities Q in [0, 1]
   satisfying M = P + Q in [0, 1] and exact 2x2 outer-product factorizability.
2. Variance-stabilized arcsin(sqrt(C / L)) proportions and scale-invariant log2(CPM + eps) matrices.
"""

from __future__ import annotations

from typing import assert_never

import numpy as np
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy.special import expit, logit
from scipy.stats import rankdata

from .power import depth_residualize_continuous_features
from .types import (
    BoolArray,
    FloatArray,
    PresenceConfig,
    PresenceMode,
    PresenceRepresentation,
)


def _resolve_library_sizes(
    count_matrix: FloatArray,
    library_sizes: Maybe[FloatArray],
) -> Result[FloatArray, str]:
    """Resolves strictly positive effective donor library sizes."""
    n_samples = int(count_matrix.shape[1])
    match library_sizes:
        case Some(libs):
            lib_arr = np.asarray(libs, dtype=np.float64).ravel()
            if int(lib_arr.shape[0]) != n_samples:
                return Failure(
                    f"library_sizes length ({lib_arr.shape[0]}) does not match sample count ({n_samples})"
                )
            if np.any(lib_arr <= 0.0) or not np.all(np.isfinite(lib_arr)):
                return Failure("library_sizes must contain strictly positive finite values")
            return Success(lib_arr)
        case Maybe.empty:
            col_sums = np.sum(count_matrix, axis=0, dtype=np.float64)
            safe_libs = np.where(col_sums > 0.0, col_sums, 1.0)
            return Success(safe_libs)
        case _ as unreachable:
            return Failure(f"Unexpected Maybe state: {unreachable}")


def _compute_empirical_bayes_presence(
    counts_act: FloatArray,
    lib_sizes: FloatArray,
    raw_binary: BoolArray,
    eb_background_rate: float,
) -> FloatArray:
    """Computes zero-inflated Poisson vs. Exponential-Poisson (r=1) hurdle posterior P(Z=1 | C, L)."""
    n_samples = float(counts_act.shape[1])
    lib_row = np.expand_dims(lib_sizes, axis=0)

    # Empirical Bayes prior prevalence pi_v with Jeffreys Beta(0.5, 0.5) smoothing
    active_counts = np.sum(raw_binary, axis=1, dtype=np.float64)
    prior_pi = np.clip((active_counts + 0.5) / (n_samples + 1.0), 0.02, 0.98)
    prior_logit = np.expand_dims(logit(prior_pi), axis=1)

    # Active rate mu_v estimated across positive observations of each neighborhood
    pos_mask = counts_act > 0.0
    pos_counts_sum = np.sum(np.where(pos_mask, counts_act, 0.0), axis=1)
    pos_libs_sum = np.sum(np.where(pos_mask, lib_row, 0.0), axis=1)
    min_active_rate = 5.0 * eb_background_rate
    mu_rate = np.maximum(
        pos_counts_sum / np.maximum(pos_libs_sum, 1.0),
        min_active_rate,
    )

    lambda_0 = np.maximum(eb_background_rate * lib_row, 1e-8)
    mu_1 = np.maximum(np.expand_dims(mu_rate, axis=1) * lib_row, 1e-6)

    # Exact log Bayes factor: log NegBin(c; r=1, mean=mu_1) - log Poisson(c; lambda_0)
    log_bf = (
        counts_act * np.log(mu_1 / lambda_0)
        - (counts_act + 1.0) * np.log1p(mu_1)
        + lambda_0
    )
    posterior = np.asarray(expit(np.clip(prior_logit + log_bf, -40.0, 40.0)), dtype=np.float64)
    return np.clip(posterior, 0.0, 1.0)


def _compute_quantile_soft_gate(arcsin_props: FloatArray) -> FloatArray:
    """Maps variance-stabilized proportions per neighborhood to [0, 1] mid-rank quantiles."""
    n_nhoods, n_samples = arcsin_props.shape
    ranks = np.vstack(
        [
            rankdata(arcsin_props[idx, :], method="average")
            for idx in range(n_nhoods)
        ]
    )
    r_min = np.min(ranks, axis=1, keepdims=True)
    r_max = np.max(ranks, axis=1, keepdims=True)
    span = r_max - r_min
    scaled = np.where(span > 1e-12, (ranks - r_min) / np.maximum(span, 1e-12), 0.0)
    return np.clip(np.asarray(scaled, dtype=np.float64), 0.0, 1.0)


def estimate_nhood_presence(
    count_matrix: FloatArray,
    library_sizes: Maybe[FloatArray] = Nothing,
    config: PresenceConfig = PresenceConfig(),
) -> Result[PresenceRepresentation, str]:
    """Computes depth-calibrated presence/absence matrices and continuous compositional features."""
    arr = np.asarray(count_matrix, dtype=np.float64)
    if arr.ndim != 2:
        return Failure(f"count_matrix must be 2D, got shape {arr.shape}")

    n_nhoods, n_samples = int(arr.shape[0]), int(arr.shape[1])
    if n_nhoods < 2 or n_samples < 4:
        return Failure(
            f"Pairwise analysis requires >= 2 neighborhoods and >= 4 samples, got ({n_nhoods}, {n_samples})"
        )
    if np.any(arr < 0.0) or not np.all(np.isfinite(arr)):
        return Failure("count_matrix must contain non-negative finite values")

    lib_res = _resolve_library_sizes(arr, library_sizes)
    match lib_res:
        case Failure(err):
            return Failure(err)
        case Success(lib_sizes):
            pass

    lib_row = np.expand_dims(lib_sizes, axis=0)
    proportions_all = arr / lib_row
    raw_binary_all = (arr >= float(config.min_cell_count)) & (proportions_all >= config.min_fraction)

    n_pos_donors = np.sum(arr > 0.0, axis=1)
    active_indices = np.flatnonzero(n_pos_donors >= config.min_active_donors).astype(np.int64)
    if int(active_indices.shape[0]) < 2:
        return Failure(
            f"Fewer than 2 active neighborhoods pass min_active_donors={config.min_active_donors} "
            f"(found {active_indices.shape[0]} of {n_nhoods})"
        )

    counts_act = arr[active_indices, :]
    props_act = proportions_all[active_indices, :]
    raw_bin_act = raw_binary_all[active_indices, :]

    # Variance-stabilized arcsin-sqrt proportions and scale-invariant log2(CPM + eps)
    clipped_props = np.clip(props_act, 0.0, 1.0)
    raw_arcsin_sqrt = np.arcsin(np.sqrt(clipped_props))
    raw_log2_cpm = np.log2(1e6 * props_act + config.prior_cpm_pseudocount)

    if config.residualize_library_depth:
        cpm_res = depth_residualize_continuous_features(raw_log2_cpm, lib_sizes)
        arc_res = depth_residualize_continuous_features(raw_arcsin_sqrt, lib_sizes)
        match (cpm_res, arc_res):
            case (Success(res_cpm), Success(res_arc)):
                log2_cpm = res_cpm
                arcsin_sqrt = np.where(
                    counts_act > 0.0,
                    np.clip(res_arc, 1e-8, 0.5 * np.pi),
                    0.0,
                )
            case (Failure(err), _):
                return Failure(err)
            case (_, Failure(err)):
                return Failure(err)
            case _ as unreachable_res:
                return Failure(f"Unexpected Result state: {unreachable_res}")
    else:
        arcsin_sqrt = raw_arcsin_sqrt
        log2_cpm = raw_log2_cpm

    # Conditional active rate mu_{v | Z=1} across positive donors for Poisson zero-detection power
    pos_mask = counts_act > 0.0
    pos_counts_sum = np.sum(np.where(pos_mask, counts_act, 0.0), axis=1, keepdims=True)
    pos_libs_sum = np.sum(np.where(pos_mask, lib_row, 0.0), axis=1, keepdims=True)
    cond_active_rate = pos_counts_sum / np.maximum(pos_libs_sum, 1.0)
    detect_floor = max(config.min_fraction, 1e-5)
    detection_power = (
        (1.0 - np.exp(-lib_row * np.maximum(cond_active_rate, detect_floor)))
        if config.use_detection_power_weights
        else np.ones_like(counts_act, dtype=np.float64)
    )
    detection_power = np.clip(detection_power, 0.0, 1.0)

    match config.mode:
        case PresenceMode.HARD_DEPTH_GATED:
            presence_prob = raw_bin_act.astype(np.float64)
            absence_prob = detection_power * (1.0 - presence_prob)
        case PresenceMode.EMPIRICAL_BAYES_POSTERIOR:
            presence_prob = _compute_empirical_bayes_presence(
                counts_act=counts_act,
                lib_sizes=lib_sizes,
                raw_binary=raw_bin_act,
                eb_background_rate=config.eb_background_rate,
            )
            absence_prob = detection_power * (1.0 - presence_prob)
        case PresenceMode.QUANTILE_SOFT_GATE:
            presence_prob = _compute_quantile_soft_gate(arcsin_sqrt)
            # For ubiquitous neighborhoods where every donor is positive, absence is relative low density
            ubiquitous_mask = np.all(counts_act > 0.0, axis=1, keepdims=True)
            eff_detect = np.where(ubiquitous_mask, 1.0, detection_power)
            absence_prob = eff_detect * (1.0 - presence_prob)
        case _ as unreachable:
            assert_never(unreachable)

    observability = np.clip(presence_prob + absence_prob, 0.0, 1.0)

    return Success(
        PresenceRepresentation(
            active_nhood_indices=active_indices,
            presence_prob=presence_prob,
            absence_prob=absence_prob,
            observability=observability,
            raw_binary_presence=raw_bin_act,
            arcsin_sqrt_prop=arcsin_sqrt,
            log2_cpm=log2_cpm,
            effective_library_sizes=lib_sizes,
        )
    )
