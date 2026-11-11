"""Outcome-blind donor/batch mixing QC, cross-compartment simplification, and IUT permutation calibration.

Implements four strictly pure, functional capabilities for `milopy.pairwise`:
1. Neighborhood donor & multi-study batch mixing diagnostics (`compute_nhood_mixing_metrics`):
   Quantifies Inverse-Simpson effective donors (`LISI_pt`), single-donor dominance, normalized Shannon
   donor entropy, and optional sub-study/batch purity (`100%` blind to clinical response labels `y`).
2. Multi-study stratum inference (`extract_inferred_substudy_strata`):
   Detects compendium sub-study identifiers (`<PMID>ggg<Patient>` or `<Study>.<Patient>`) for
   Cochran-Mantel-Haenszel (`CMH`) stratification.
3. Canonical state-pair bundle consolidation (`consolidate_state_pair_bundles`,
   `deduplicate_pairs_to_lead_bundles`):
   Collapses parallel sub-neighborhood chords into unique `(state_a, state_b, quadrant, direction)`
   ecosystem interactions with Simes (1986) PRDS-valid combined p-values.
4. Stratified IUT permutation & LODO/LOSO jackknife verification
   (`calibrate_candidates_stratified_iut_permutation`):
   Vectorized within-stratum donor-label shuffles (`B` permutations) yielding empirical `pvalue_perm_iut`,
   Westfall-Young step-down/max-T `fwer_perm_iut`, and Leave-One-Donor-Out (`LODO`) +
   Leave-One-Stratum-Out (`LOSO`) stability metrics.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import NamedTuple

import numpy as np
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
from scipy import sparse as sp

from .overlap import compute_weighted_bh_fdr
from .types import (
    ALL_QUADRANT_STATES,
    BoolArray,
    ContinuousAxis,
    EnrichmentDirection,
    FloatArray,
    IntArray,
    MixingQCConfig,
    PresenceRepresentation,
    QuadrantState,
    SparseMatrix,
    TmeCompartment,
    adjacent_flips,
    classify_tme_compartment,
    swap_quadrant,
)


class _CandidateQuadVectors(NamedTuple):
    """Pre-extracted donor-level quadrant evidence arrays for C candidate rows across J samples."""

    x_q: FloatArray
    x_fv: FloatArray
    x_fu: FloatArray
    m_tot: FloatArray
    obs_sign: FloatArray


def _strip_eval_suffix(sample_id: str) -> str:
    """Strips optional `__responder` / `__non-responder` evaluation suffixes from sample IDs."""
    raw = str(sample_id).strip()
    if "__" in raw:
        head, tail = raw.rsplit("__", 1)
        if tail.strip().lower() in ("responder", "non-responder"):
            return head.strip()
    return raw


def _validate_strata_candidate(prefixes: tuple[str, ...]) -> Maybe[tuple[str, ...]]:
    """Accepts inferred strata iff >= 2 distinct groups exist and >= 2 groups contain >= 2 donors."""
    counts: dict[str, int] = {}
    for p in prefixes:
        counts[p] = counts.get(p, 0) + 1
    if len(counts) < 2:
        return Nothing
    multi_donor_groups = sum(1 for c in counts.values() if c >= 2)
    return Some(prefixes) if multi_donor_groups >= 2 else Nothing


def extract_inferred_substudy_strata(
    sample_ids: Sequence[str],
) -> Maybe[tuple[str, ...]]:
    """Infers sub-study / cohort batch strata from structured compendium sample identifiers.

    Supports `<PMID>ggg<PatientID>` (e.g. Gondal `CELLxGENE_7b20c613`) and
    `<Study_Year_Venue>.<PatientID>` (e.g. `CELLxGENE_05a8c945`), returning `Nothing`
    for homogeneous single-study cohorts (`P01`, `GSM...`, etc.).
    """
    n_samples = len(sample_ids)
    if n_samples < 4:
        return Nothing

    cleaned = tuple(_strip_eval_suffix(str(s)) for s in sample_ids)

    n_ggg = sum(1 for s in cleaned if "ggg" in s)
    if n_ggg >= int(np.ceil(0.80 * float(n_samples))):
        ggg_prefixes = tuple(
            s.split("ggg", 1)[0].strip() if "ggg" in s else "unassigned_stratum"
            for s in cleaned
        )
        res_ggg = _validate_strata_candidate(ggg_prefixes)
        if isinstance(res_ggg, Some):
            return res_ggg

    dot_candidates: list[str] = []
    valid_dot_count = 0
    for s in cleaned:
        if "." in s:
            head = s.split(".", 1)[0].strip()
            if len(head) >= 3 and any(ch.isalpha() for ch in head):
                dot_candidates.append(head)
                valid_dot_count += 1
            else:
                dot_candidates.append("unassigned_stratum")
        else:
            dot_candidates.append("unassigned_stratum")

    if valid_dot_count >= int(np.ceil(0.80 * float(n_samples))):
        return _validate_strata_candidate(tuple(dot_candidates))

    return Nothing


def resolve_effective_mixing_thresholds(
    n_total_donors: int,
    config: MixingQCConfig,
) -> tuple[int, float, float]:
    """Resolves `(min_donors, min_lisi, max_dominance)` with optional small-cohort scaling."""
    if (not config.adaptive_small_cohort_scaling) or n_total_donors >= 12:
        return (
            config.min_distinct_donors,
            config.min_effective_donors_lisi,
            config.max_single_donor_dominance,
        )

    j_safe = max(n_total_donors, 1)
    eff_min_don = min(config.min_distinct_donors, max(2, j_safe // 3))
    eff_min_lisi = min(config.min_effective_donors_lisi, max(1.30, 0.22 * float(j_safe)))
    eff_max_dom = max(config.max_single_donor_dominance, min(0.85, 3.5 / float(j_safe)))
    return eff_min_don, eff_min_lisi, eff_max_dom


def compute_nhood_mixing_metrics(
    nhoods_mat: SparseMatrix,
    count_matrix: FloatArray,
    sample_ids: Sequence[str],
    batch_labels: Maybe[Sequence[str]] = Nothing,
    nhood_labels: Maybe[Sequence[str]] = Nothing,
    config: MixingQCConfig = MixingQCConfig(),
) -> Result[pl.DataFrame, str]:
    """Computes strictly outcome-blind donor and batch mixing metrics for every neighborhood.

    Metrics per neighborhood `v`:
    - `n_cells`: total cell count across donors (`sum_j C_{v,j}`)
    - `n_distinct_donors`: donors contributing `C_{v,j} > 0`
    - `effective_donors_lisi`: Inverse-Simpson effective donors `1 / sum_j p_{v,j}^2` in `[1, J]`
    - `shannon_donor_entropy_norm`: normalized Shannon entropy in `[0, 1]`
    - `max_single_donor_fraction`: largest single-donor share `max_j p_{v,j}`
    - `n_distinct_batches`, `effective_batches_lisi`, `max_batch_fraction`: sub-study/batch mixing
    - `passes_mixing_qc`: outcome-blind QC boolean gate
    """
    counts_arr = np.asarray(count_matrix, dtype=np.float64)
    if counts_arr.ndim != 2:
        return Failure(f"count_matrix must be 2D, got shape {counts_arr.shape}")

    n_nhoods, n_samples = int(counts_arr.shape[0]), int(counts_arr.shape[1])
    if int(nhoods_mat.shape[1]) != n_nhoods:
        return Failure(
            f"nhoods_mat columns ({nhoods_mat.shape[1]}) must match count_matrix rows ({n_nhoods})"
        )
    if len(sample_ids) != n_samples:
        return Failure(
            f"sample_ids length ({len(sample_ids)}) must match count_matrix columns ({n_samples})"
        )
    if np.any(counts_arr < -1e-9):
        return Failure("count_matrix cannot contain negative counts")

    counts_pos = np.maximum(counts_arr, 0.0)
    totals = np.sum(counts_pos, axis=1)
    non_empty = totals > 0.0
    safe_totals = np.where(non_empty, totals, 1.0)
    donor_props = counts_pos / np.expand_dims(safe_totals, axis=1)

    n_distinct_donors = np.sum(counts_pos > 0.0, axis=1).astype(np.int64)
    simpson_donor = np.sum(donor_props * donor_props, axis=1)
    effective_lisi = np.where(
        non_empty,
        1.0 / np.maximum(simpson_donor, 1e-12),
        0.0,
    )
    max_donor_frac = np.where(
        non_empty,
        np.max(donor_props, axis=1),
        1.0,
    )

    p_pos = donor_props > 1e-15
    p_safe = np.where(p_pos, donor_props, 1.0)
    raw_entropy = -np.sum(np.where(p_pos, donor_props * np.log(p_safe), 0.0), axis=1)
    denom_log = float(np.log(max(n_samples, 2)))
    norm_entropy = np.where(non_empty, np.clip(raw_entropy / denom_log, 0.0, 1.0), 0.0)

    resolved_batches: Maybe[tuple[str, ...]] = Nothing
    match batch_labels:
        case Some(b_seq):
            if len(b_seq) != n_samples:
                return Failure(
                    f"batch_labels length ({len(b_seq)}) must match sample count ({n_samples})"
                )
            resolved_batches = Some(tuple(str(x) for x in b_seq))
        case Maybe.empty:
            resolved_batches = extract_inferred_substudy_strata(sample_ids)
        case _ as unreachable_b:
            return Failure(f"Unexpected Maybe state: {unreachable_b}")

    match resolved_batches:
        case Some(b_tup):
            uniq_b = tuple(sorted(set(b_tup)))
            n_batches = len(uniq_b)
            b_to_idx = {name: idx for idx, name in enumerate(uniq_b)}
            b_codes = np.array([b_to_idx[name] for name in b_tup], dtype=np.int64)
            b_ind = sp.csr_matrix(
                (np.ones(n_samples, dtype=np.float64), (np.arange(n_samples, dtype=np.int64), b_codes)),
                shape=(n_samples, n_batches),
            )
            batch_counts = np.asarray(counts_pos @ b_ind.toarray(), dtype=np.float64)
            batch_props = batch_counts / np.expand_dims(safe_totals, axis=1)
            n_distinct_batches = np.sum(batch_counts > 0.0, axis=1).astype(np.int64)
            batch_simpson = np.where(non_empty, np.sum(batch_props * batch_props, axis=1), 1.0)
            eff_batches_lisi = np.where(non_empty, 1.0 / np.maximum(batch_simpson, 1e-12), 0.0)
            max_batch_frac = np.where(non_empty, np.max(batch_props, axis=1), 1.0)
            has_multi_batch = n_batches >= 2
        case Maybe.empty:
            n_distinct_batches = np.where(non_empty, 1, 0).astype(np.int64)
            eff_batches_lisi = np.where(non_empty, 1.0, 0.0)
            max_batch_frac = np.ones(n_nhoods, dtype=np.float64)
            has_multi_batch = False
        case _ as unreachable_rb:
            return Failure(f"Unexpected Maybe state: {unreachable_rb}")

    match nhood_labels:
        case Some(lbl_seq):
            if len(lbl_seq) != n_nhoods:
                return Failure(
                    f"nhood_labels length ({len(lbl_seq)}) must match neighborhoods ({n_nhoods})"
                )
            resolved_lbls = tuple(str(x) for x in lbl_seq)
        case Maybe.empty:
            resolved_lbls = tuple(f"Nhood_{i}" for i in range(n_nhoods))
        case _ as unreachable_nl:
            return Failure(f"Unexpected Maybe state: {unreachable_nl}")

    compartments = tuple(classify_tme_compartment(lbl).value for lbl in resolved_lbls)
    eff_min_don, eff_min_lisi, eff_max_dom = resolve_effective_mixing_thresholds(n_samples, config)

    passes_donor_mix = (
        non_empty
        & (n_distinct_donors >= eff_min_don)
        & (effective_lisi >= (eff_min_lisi - 1e-12))
        & (max_donor_frac <= (eff_max_dom + 1e-12))
    )

    match config.max_batch_fraction:
        case Some(max_bf):
            passes_batch_mix = (
                (max_batch_frac <= (float(max_bf) + 1e-12))
                if has_multi_batch
                else np.ones(n_nhoods, dtype=bool)
            )
        case Maybe.empty:
            passes_batch_mix = np.ones(n_nhoods, dtype=bool)
        case _ as unreachable_mb:
            return Failure(f"Unexpected Maybe state: {unreachable_mb}")

    passes_qc = passes_donor_mix & passes_batch_mix
    nhood_sizes_mat = np.asarray((sp.csc_matrix(nhoods_mat) > 0).sum(axis=0), dtype=np.float64).ravel()

    return Success(
        pl.DataFrame(
            {
                "nhood_idx": np.arange(n_nhoods, dtype=np.int64),
                "nhood_label": list(resolved_lbls),
                "compartment": list(compartments),
                "nhood_matrix_cells": nhood_sizes_mat,
                "n_cells": totals,
                "n_distinct_donors": n_distinct_donors,
                "effective_donors_lisi": effective_lisi,
                "shannon_donor_entropy_norm": norm_entropy,
                "max_single_donor_fraction": max_donor_frac,
                "n_distinct_batches": n_distinct_batches,
                "effective_batches_lisi": eff_batches_lisi,
                "max_batch_fraction": max_batch_frac,
                "passes_mixing_qc": passes_qc,
            }
        )
    )


def _canonicalize_pair_rows(pairs_df: pl.DataFrame) -> Result[pl.DataFrame, str]:
    """Attaches canonical ordered `(state_a, state_b, canonical_quadrant, compartment_a, compartment_b)`."""
    req_cols = {
        "nhood_u",
        "nhood_v",
        "quadrant",
        "direction",
        "cell_type_u",
        "cell_type_v",
        "pvalue_iut",
        "fdr_iut",
        "log_or_quadrant",
        "is_significant",
    }
    if pairs_df.is_empty() and pairs_df.width == 0:
        return Success(pairs_df)

    missing = sorted(req_cols - set(pairs_df.columns))
    if missing:
        return Failure(f"pairs_df missing required columns for bundle consolidation: {missing}")

    if pairs_df.is_empty():
        return Success(
            pairs_df.with_columns(
                pl.Series("state_a", [], dtype=pl.String),
                pl.Series("state_b", [], dtype=pl.String),
                pl.Series("compartment_a", [], dtype=pl.String),
                pl.Series("compartment_b", [], dtype=pl.String),
                pl.Series("canonical_quadrant", [], dtype=pl.String),
            )
        )

    ct_u = [str(x) for x in pairs_df["cell_type_u"].to_list()]
    ct_v = [str(x) for x in pairs_df["cell_type_v"].to_list()]
    quads = [str(x) for x in pairs_df["quadrant"].to_list()]

    valid_q = {q.value for q in ALL_QUADRANT_STATES}
    state_a_list: list[str] = []
    state_b_list: list[str] = []
    comp_a_list: list[str] = []
    comp_b_list: list[str] = []
    canon_q_list: list[str] = []

    for u_lbl, v_lbl, q_str in zip(ct_u, ct_v, quads):
        comp_u = classify_tme_compartment(u_lbl).value
        comp_v = classify_tme_compartment(v_lbl).value
        if (u_lbl, comp_u) <= (v_lbl, comp_v):
            state_a_list.append(u_lbl)
            state_b_list.append(v_lbl)
            comp_a_list.append(comp_u)
            comp_b_list.append(comp_v)
            canon_q_list.append(q_str)
        else:
            state_a_list.append(v_lbl)
            state_b_list.append(u_lbl)
            comp_a_list.append(comp_v)
            comp_b_list.append(comp_u)
            swapped = swap_quadrant(QuadrantState(q_str)).value if q_str in valid_q else q_str
            canon_q_list.append(swapped)

    return Success(
        pairs_df.with_columns(
            pl.Series("state_a", state_a_list, dtype=pl.String),
            pl.Series("state_b", state_b_list, dtype=pl.String),
            pl.Series("compartment_a", comp_a_list, dtype=pl.String),
            pl.Series("compartment_b", comp_b_list, dtype=pl.String),
            pl.Series("canonical_quadrant", canon_q_list, dtype=pl.String),
        )
    )


def _simes_p_value(p_vals: FloatArray) -> float:
    """Computes Simes (1986) combined p-value `min_r (m * p_(r) / r)` valid under PRDS."""
    m = int(p_vals.shape[0])
    if m == 0:
        return 1.0
    p_sorted = np.sort(np.clip(p_vals, 0.0, 1.0))
    ranks = np.arange(1, m + 1, dtype=np.float64)
    return float(np.clip(np.min((float(m) * p_sorted) / ranks), 0.0, 1.0))


def consolidate_state_pair_bundles(
    pairs_df: pl.DataFrame,
    fdr_threshold: float = 0.10,
    significant_only: bool = True,
) -> Result[pl.DataFrame, str]:
    """Consolidates parallel sub-neighborhood pairs into canonical `(state_a, state_b, regime)` bundles."""
    canon_res = _canonicalize_pair_rows(pairs_df)
    match canon_res:
        case Failure(err):
            return Failure(err)
        case Success(canon_df):
            pass

    source_df = (
        canon_df.filter(pl.col("is_significant"))
        if significant_only
        else canon_df
    )
    if source_df.is_empty():
        return Success(
            pl.DataFrame(
                schema={
                    "bundle_id": pl.Int64,
                    "state_a": pl.String,
                    "state_b": pl.String,
                    "compartment_a": pl.String,
                    "compartment_b": pl.String,
                    "canonical_quadrant": pl.String,
                    "direction": pl.String,
                    "lead_nhood_u": pl.Int64,
                    "lead_nhood_v": pl.Int64,
                    "lead_quadrant": pl.String,
                    "n_supporting_subpairs": pl.Int64,
                    "n_significant_subpairs": pl.Int64,
                    "lead_log_or_quadrant": pl.Float64,
                    "median_log_or_quadrant": pl.Float64,
                    "min_pvalue_iut": pl.Float64,
                    "pvalue_simes_iut": pl.Float64,
                    "min_fdr_iut": pl.Float64,
                    "fdr_bundle_simes": pl.Float64,
                    "is_cross_compartment": pl.Boolean,
                    "is_bundle_significant": pl.Boolean,
                }
            )
        )

    ordered = (
        source_df.lazy()
        .sort(
            ["fdr_iut", "pvalue_iut", "nhood_u", "nhood_v"],
            descending=[False, False, False, False],
        )
        .collect()
    )
    group_cols = ["state_a", "state_b", "compartment_a", "compartment_b", "canonical_quadrant", "direction"]
    partitions = ordered.partition_by(group_cols, maintain_order=True)

    rows_out: list[dict[str, object]] = []
    for part in partitions:
        lead = part.row(0, named=True)
        p_arr = part["pvalue_iut"].cast(pl.Float64).to_numpy()
        lor_arr = part["log_or_quadrant"].cast(pl.Float64).to_numpy()
        fdr_arr = part["fdr_iut"].cast(pl.Float64).to_numpy()
        sig_arr = part["is_significant"].to_numpy()
        comp_a = str(lead["compartment_a"])
        comp_b = str(lead["compartment_b"])
        rows_out.append(
            {
                "state_a": str(lead["state_a"]),
                "state_b": str(lead["state_b"]),
                "compartment_a": comp_a,
                "compartment_b": comp_b,
                "canonical_quadrant": str(lead["canonical_quadrant"]),
                "direction": str(lead["direction"]),
                "lead_nhood_u": int(lead["nhood_u"]),
                "lead_nhood_v": int(lead["nhood_v"]),
                "lead_quadrant": str(lead["quadrant"]),
                "n_supporting_subpairs": int(part.height),
                "n_significant_subpairs": int(np.sum(sig_arr)),
                "lead_log_or_quadrant": float(lead["log_or_quadrant"]),
                "median_log_or_quadrant": float(np.median(lor_arr)),
                "min_pvalue_iut": float(np.min(p_arr)),
                "pvalue_simes_iut": _simes_p_value(p_arr),
                "min_fdr_iut": float(np.min(fdr_arr)),
                "is_cross_compartment": (
                    (comp_a != comp_b)
                    and (comp_a != TmeCompartment.OTHER.value)
                    and (comp_b != TmeCompartment.OTHER.value)
                ),
            }
        )

    raw_bundles = pl.DataFrame(rows_out)
    p_simes = raw_bundles["pvalue_simes_iut"].cast(pl.Float64).to_numpy()
    fdr_simes = compute_weighted_bh_fdr(p_simes, np.ones_like(p_simes))
    min_fdr = raw_bundles["min_fdr_iut"].cast(pl.Float64).to_numpy()
    is_b_sig = (min_fdr < fdr_threshold) & (fdr_simes < fdr_threshold)

    sorted_bundles = (
        raw_bundles.with_columns(
            pl.Series("fdr_bundle_simes", fdr_simes, dtype=pl.Float64),
            pl.Series("is_bundle_significant", is_b_sig, dtype=pl.Boolean),
        )
        .lazy()
        .sort(
            ["is_bundle_significant", "min_fdr_iut", "pvalue_simes_iut", "n_supporting_subpairs"],
            descending=[True, False, False, True],
        )
        .with_row_index("bundle_id")
        .with_columns(pl.col("bundle_id").cast(pl.Int64))
        .collect()
    )
    return Success(sorted_bundles)


def deduplicate_pairs_to_lead_bundles(
    pairs_df: pl.DataFrame,
    require_cross_compartment: bool = False,
    one_per_state_pair: bool = False,
) -> Result[pl.DataFrame, str]:
    """Filters `pairs_df` to strictly 1 lead `(nhood_u, nhood_v)` chord per state-pair bundle."""
    canon_res = _canonicalize_pair_rows(pairs_df)
    match canon_res:
        case Failure(err):
            return Failure(err)
        case Success(canon_df):
            pass

    if canon_df.is_empty():
        return Success(pairs_df)

    base_df = (
        canon_df.filter(
            (pl.col("compartment_a") != pl.col("compartment_b"))
            & (pl.col("compartment_a") != TmeCompartment.OTHER.value)
            & (pl.col("compartment_b") != TmeCompartment.OTHER.value)
        )
        if require_cross_compartment
        else canon_df
    )
    if base_df.is_empty():
        return Success(pairs_df.head(0))

    subset_cols = (
        ["state_a", "state_b"]
        if one_per_state_pair
        else ["state_a", "state_b", "canonical_quadrant", "direction"]
    )
    dedup = (
        base_df.lazy()
        .sort(
            ["is_significant", "fdr_iut", "pvalue_iut", "nhood_u", "nhood_v"],
            descending=[True, False, False, False, False],
        )
        .unique(
            subset=subset_cols,
            keep="first",
            maintain_order=True,
        )
        .drop(["state_a", "state_b", "compartment_a", "compartment_b", "canonical_quadrant"])
        .collect()
    )
    return Success(dedup)


def _flip_enrichment_direction_str(direction_str: str) -> str:
    """Flips `responder <-> non-responder` when swapping antisymmetric log-ratio endpoints."""
    if direction_str == EnrichmentDirection.RESPONDER.value:
        return EnrichmentDirection.NON_RESPONDER.value
    if direction_str == EnrichmentDirection.NON_RESPONDER.value:
        return EnrichmentDirection.RESPONDER.value
    return direction_str


def _canonicalize_continuous_pair_rows(
    continuous_pairs_df: pl.DataFrame,
) -> Result[pl.DataFrame, str]:
    """Attaches swap-invariant `(state_a, state_b, compartment_a, compartment_b, canonical_direction)` columns."""
    req_cols = {
        "nhood_u",
        "nhood_v",
        "axis",
        "direction",
        "cell_type_u",
        "cell_type_v",
        "pvalue_iut",
        "fdr_iut",
        "is_significant",
    }
    if continuous_pairs_df.is_empty() and continuous_pairs_df.width == 0:
        return Success(continuous_pairs_df)

    missing = sorted(req_cols - set(continuous_pairs_df.columns))
    if missing:
        return Failure(f"continuous_pairs_df missing required columns for consolidation: {missing}")

    if continuous_pairs_df.is_empty():
        return Success(
            continuous_pairs_df.with_columns(
                pl.Series("state_a", [], dtype=pl.String),
                pl.Series("state_b", [], dtype=pl.String),
                pl.Series("compartment_a", [], dtype=pl.String),
                pl.Series("compartment_b", [], dtype=pl.String),
                pl.Series("canonical_direction", [], dtype=pl.String),
            )
        )

    ct_u = [str(x) for x in continuous_pairs_df["cell_type_u"].to_list()]
    ct_v = [str(x) for x in continuous_pairs_df["cell_type_v"].to_list()]
    axes = [str(x) for x in continuous_pairs_df["axis"].to_list()]
    dirs = [str(x) for x in continuous_pairs_df["direction"].to_list()]

    state_a_list: list[str] = []
    state_b_list: list[str] = []
    comp_a_list: list[str] = []
    comp_b_list: list[str] = []
    canon_dir_list: list[str] = []

    for u_lbl, v_lbl, ax_str, dir_str in zip(ct_u, ct_v, axes, dirs):
        comp_u = classify_tme_compartment(u_lbl).value
        comp_v = classify_tme_compartment(v_lbl).value
        if (u_lbl, comp_u) <= (v_lbl, comp_v):
            state_a_list.append(u_lbl)
            state_b_list.append(v_lbl)
            comp_a_list.append(comp_u)
            comp_b_list.append(comp_v)
            canon_dir_list.append(dir_str)
        else:
            state_a_list.append(v_lbl)
            state_b_list.append(u_lbl)
            comp_a_list.append(comp_v)
            comp_b_list.append(comp_u)
            canon_dir = (
                _flip_enrichment_direction_str(dir_str)
                if ax_str == ContinuousAxis.LOG_RATIO_BALANCE.value
                else dir_str
            )
            canon_dir_list.append(canon_dir)

    return Success(
        continuous_pairs_df.with_columns(
            pl.Series("state_a", state_a_list, dtype=pl.String),
            pl.Series("state_b", state_b_list, dtype=pl.String),
            pl.Series("compartment_a", comp_a_list, dtype=pl.String),
            pl.Series("compartment_b", comp_b_list, dtype=pl.String),
            pl.Series("canonical_direction", canon_dir_list, dtype=pl.String),
        )
    )


def deduplicate_continuous_pairs_to_lead_bundles(
    continuous_pairs_df: pl.DataFrame,
    require_cross_compartment: bool = False,
    one_per_state_pair: bool = False,
) -> Result[pl.DataFrame, str]:
    """Filters `continuous_pairs_df` to 1 lead `(nhood_u, nhood_v)` chord per canonical state-pair bundle."""
    canon_res = _canonicalize_continuous_pair_rows(continuous_pairs_df)
    match canon_res:
        case Failure(err):
            return Failure(err)
        case Success(canon_df):
            pass

    if canon_df.is_empty():
        return Success(continuous_pairs_df)

    base_df = (
        canon_df.filter(
            (pl.col("compartment_a") != pl.col("compartment_b"))
            & (pl.col("compartment_a") != TmeCompartment.OTHER.value)
            & (pl.col("compartment_b") != TmeCompartment.OTHER.value)
        )
        if require_cross_compartment
        else canon_df
    )
    if base_df.is_empty():
        return Success(continuous_pairs_df.head(0))

    subset_cols = (
        ["state_a", "state_b"]
        if one_per_state_pair
        else ["state_a", "state_b", "axis", "canonical_direction"]
    )
    dedup = (
        base_df.lazy()
        .sort(
            ["is_significant", "fdr_iut", "pvalue_iut", "nhood_u", "nhood_v"],
            descending=[True, False, False, False, False],
        )
        .unique(
            subset=subset_cols,
            keep="first",
            maintain_order=True,
        )
        .drop(["state_a", "state_b", "compartment_a", "compartment_b", "canonical_direction"])
        .collect()
    )
    return Success(dedup)


def _quad_evidence_for_state(
    state: QuadrantState,
    p_u: FloatArray,
    q_u: FloatArray,
    p_v: FloatArray,
    q_v: FloatArray,
) -> FloatArray:
    """Computes `(C, J)` soft/hard quadrant donor evidence matrix for `state`."""
    match state:
        case QuadrantState.BOTH_PRESENT:
            return p_u * p_v
        case QuadrantState.FIRST_ONLY:
            return p_u * q_v
        case QuadrantState.SECOND_ONLY:
            return q_u * p_v
        case QuadrantState.BOTH_MISSING:
            return q_u * q_v


def _extract_candidate_vectors(
    candidates_df: pl.DataFrame,
    presence: PresenceRepresentation,
) -> Result[_CandidateQuadVectors, str]:
    """Maps global `(nhood_u, nhood_v, quadrant, direction)` rows into `(C, J)` evidence arrays."""
    act_map = {int(g_idx): loc for loc, g_idx in enumerate(presence.active_nhood_indices.tolist())}
    n_cand = candidates_df.height
    n_samples = int(presence.presence_prob.shape[1])

    u_locs = np.empty(n_cand, dtype=np.int64)
    v_locs = np.empty(n_cand, dtype=np.int64)
    q_states: list[QuadrantState] = []
    obs_sign = np.empty(n_cand, dtype=np.float64)

    for idx, row in enumerate(
        candidates_df.select(["nhood_u", "nhood_v", "quadrant", "direction"]).to_dicts()
    ):
        gu = int(row["nhood_u"])
        gv = int(row["nhood_v"])
        if gu not in act_map or gv not in act_map:
            return Failure(f"Candidate pair ({gu}, {gv}) not present in active_nhood_indices")
        u_locs[idx] = act_map[gu]
        v_locs[idx] = act_map[gv]
        q_states.append(QuadrantState(str(row["quadrant"])))
        obs_sign[idx] = (
            1.0 if str(row["direction"]) == EnrichmentDirection.RESPONDER.value else -1.0
        )

    p_u = presence.presence_prob[u_locs, :]
    q_u = presence.absence_prob[u_locs, :]
    p_v = presence.presence_prob[v_locs, :]
    q_v = presence.absence_prob[v_locs, :]
    m_tot = presence.observability[u_locs, :] * presence.observability[v_locs, :]

    x_q = np.empty((n_cand, n_samples), dtype=np.float64)
    x_fv = np.empty((n_cand, n_samples), dtype=np.float64)
    x_fu = np.empty((n_cand, n_samples), dtype=np.float64)

    for q_st in ALL_QUADRANT_STATES:
        mask = np.array([qs == q_st for qs in q_states], dtype=bool)
        if bool(np.any(mask)):
            fv_st, fu_st = adjacent_flips(q_st)
            x_q[mask, :] = _quad_evidence_for_state(q_st, p_u[mask, :], q_u[mask, :], p_v[mask, :], q_v[mask, :])
            x_fv[mask, :] = _quad_evidence_for_state(fv_st, p_u[mask, :], q_u[mask, :], p_v[mask, :], q_v[mask, :])
            x_fu[mask, :] = _quad_evidence_for_state(fu_st, p_u[mask, :], q_u[mask, :], p_v[mask, :], q_v[mask, :])

    return Success(_CandidateQuadVectors(x_q=x_q, x_fv=x_fv, x_fu=x_fu, m_tot=m_tot, obs_sign=obs_sign))


def _cmh_z_from_num_var(num: FloatArray, var: FloatArray) -> FloatArray:
    """Computes stabilized CMH Z-statistic from summed numerator and hypergeometric variance."""
    return np.where(var > 1e-12, num / np.sqrt(np.maximum(var, 1e-15)), 0.0)


def _stratum_cmh_components_2d(
    a_mat: FloatArray,
    r1_col: FloatArray,
    c_mat: FloatArray,
    r0_col: FloatArray,
) -> tuple[FloatArray, FloatArray]:
    """Computes `(C, B)` CMH numerator and variance for one stratum across `B` label vectors."""
    b_mat = np.maximum(r1_col - a_mat, 0.0)
    d_mat = np.maximum(r0_col - c_mat, 0.0)
    n_raw = a_mat + b_mat + c_mat + d_mat
    safe_n = np.maximum(n_raw, 1e-12)
    c1 = a_mat + c_mat
    c0 = b_mat + d_mat
    exp_a = (r1_col * c1) / safe_n
    denom = safe_n * safe_n * np.maximum(n_raw - 1.0, 1.0)
    active = n_raw > 1.0
    num = np.where(active, a_mat - exp_a, 0.0)
    var = np.where(active, (r1_col * r0_col * c1 * c0) / np.maximum(denom, 1e-12), 0.0)
    return num, np.maximum(var, 0.0)


def calibrate_candidates_stratified_iut_permutation(
    candidates_df: pl.DataFrame,
    presence: PresenceRepresentation,
    response_binary: FloatArray,
    strata_codes: Maybe[IntArray] = Nothing,
    n_permutations: int = 1000,
    seed: int = 42,
    fdr_threshold: float = 0.10,
    jackknife_z_floor: float = 1.6448536269514722,
) -> Result[pl.DataFrame, str]:
    """Verifies candidate pairs via stratified donor-label IUT permutations & LODO/LOSO jackknife.

    Appends:
    - `z_iut_cmh`: exact observed 3-way stratified Berger IUT score `min(s*Z_q, s*Z_{v|u}, s*Z_{u|v})`
    - `pvalue_perm_iut`: exact within-stratum donor-label permutation p-value (`(1 + exceed) / (1 + B)`)
    - `fdr_perm_iut`: Benjamini-Hochberg adjusted empirical IUT q-value
    - `fwer_perm_iut`: Westfall-Young max-statistic family-wise error rate across all evaluated candidates
    - `lodo_min_z_iut`, `lodo_sign_consistency`: Leave-One-Donor-Out minimum IUT Z and sign stability
    - `loso_min_z_iut`, `loso_sign_consistency`: Leave-One-Stratum-Out minimum IUT Z and sign stability
    - `is_statistically_verified`: joint pass across IUT FDR, stratified permutation FDR, and jackknife
    """
    if n_permutations < 10:
        return Failure(f"n_permutations must be >= 10, got {n_permutations}")

    if candidates_df.is_empty():
        return Success(
            candidates_df.with_columns(
                pl.Series("z_iut_cmh", [], dtype=pl.Float64),
                pl.Series("pvalue_perm_iut", [], dtype=pl.Float64),
                pl.Series("fdr_perm_iut", [], dtype=pl.Float64),
                pl.Series("fwer_perm_iut", [], dtype=pl.Float64),
                pl.Series("lodo_min_z_iut", [], dtype=pl.Float64),
                pl.Series("lodo_sign_consistency", [], dtype=pl.Float64),
                pl.Series("loso_min_z_iut", [], dtype=pl.Float64),
                pl.Series("loso_sign_consistency", [], dtype=pl.Float64),
                pl.Series("is_statistically_verified", [], dtype=pl.Boolean),
            )
        )

    y_arr = np.asarray(response_binary, dtype=np.float64).ravel()
    n_samples = int(y_arr.shape[0])
    match strata_codes:
        case Some(codes):
            c_arr = np.asarray(codes, dtype=np.int64).ravel()
            if int(c_arr.shape[0]) != n_samples:
                return Failure("strata_codes length mismatch with response_binary")
            raw_groups = [np.flatnonzero(c_arr == c_val).astype(np.int64) for c_val in np.unique(c_arr)]
        case Maybe.empty:
            raw_groups = [np.arange(n_samples, dtype=np.int64)]
        case _ as unreachable_sc:
            return Failure(f"Unexpected Maybe state: {unreachable_sc}")

    informative_groups = [
        grp
        for grp in raw_groups
        if bool(np.any(y_arr[grp] > 0.5)) and bool(np.any(y_arr[grp] <= 0.5))
    ]
    if not informative_groups:
        return Failure("No informative strata contain both responder and non-responder donors")

    vec_res = _extract_candidate_vectors(candidates_df, presence)
    match vec_res:
        case Failure(err):
            return Failure(err)
        case Success(vecs):
            pass

    n_cand = candidates_df.height
    k_strata = len(informative_groups)
    rng = np.random.default_rng(seed)

    # Accumulate observed (col 0) + B within-stratum permutations (cols 1..B) in one pass
    b_total = 1 + n_permutations
    num_q_all = np.zeros((n_cand, b_total), dtype=np.float64)
    var_q_all = np.zeros((n_cand, b_total), dtype=np.float64)
    num_fv_all = np.zeros((n_cand, b_total), dtype=np.float64)
    var_fv_all = np.zeros((n_cand, b_total), dtype=np.float64)
    num_fu_all = np.zeros((n_cand, b_total), dtype=np.float64)
    var_fu_all = np.zeros((n_cand, b_total), dtype=np.float64)

    # Store per-stratum observed CMH components for exact O(1) LODO & LOSO jackknife
    obs_num_q_s = np.zeros((k_strata, n_cand), dtype=np.float64)
    obs_var_q_s = np.zeros((k_strata, n_cand), dtype=np.float64)
    obs_num_fv_s = np.zeros((k_strata, n_cand), dtype=np.float64)
    obs_var_fv_s = np.zeros((k_strata, n_cand), dtype=np.float64)
    obs_num_fu_s = np.zeros((k_strata, n_cand), dtype=np.float64)
    obs_var_fu_s = np.zeros((k_strata, n_cand), dtype=np.float64)

    for s_idx, grp in enumerate(informative_groups):
        g_len = int(grp.shape[0])
        yg = (y_arr[grp] > 0.5).astype(np.float64)
        y_mat = np.empty((g_len, b_total), dtype=np.float64)
        y_mat[:, 0] = yg
        for b_i in range(1, b_total):
            y_mat[:, b_i] = rng.permutation(yg)

        xq_g = vecs.x_q[:, grp]
        xfv_g = vecs.x_fv[:, grp]
        xfu_g = vecs.x_fu[:, grp]
        mtot_g = vecs.m_tot[:, grp]

        r1_q = np.sum(xq_g, axis=1, keepdims=True)
        r0_notq = np.maximum(np.sum(mtot_g, axis=1, keepdims=True) - r1_q, 0.0)
        r0_fv = np.sum(xfv_g, axis=1, keepdims=True)
        r0_fu = np.sum(xfu_g, axis=1, keepdims=True)

        a_q = xq_g @ y_mat
        c_notq = np.maximum((mtot_g @ y_mat) - a_q, 0.0)
        c_fv = xfv_g @ y_mat
        c_fu = xfu_g @ y_mat

        nq, vq = _stratum_cmh_components_2d(a_q, r1_q, c_notq, r0_notq)
        nfv, vfv = _stratum_cmh_components_2d(a_q, r1_q, c_fv, r0_fv)
        nfu, vfu = _stratum_cmh_components_2d(a_q, r1_q, c_fu, r0_fu)

        num_q_all = num_q_all + nq
        var_q_all = var_q_all + vq
        num_fv_all = num_fv_all + nfv
        var_fv_all = var_fv_all + vfv
        num_fu_all = num_fu_all + nfu
        var_fu_all = var_fu_all + vfu

        obs_num_q_s[s_idx, :] = nq[:, 0]
        obs_var_q_s[s_idx, :] = vq[:, 0]
        obs_num_fv_s[s_idx, :] = nfv[:, 0]
        obs_var_fv_s[s_idx, :] = vfv[:, 0]
        obs_num_fu_s[s_idx, :] = nfu[:, 0]
        obs_var_fu_s[s_idx, :] = vfu[:, 0]

    zq_all = _cmh_z_from_num_var(num_q_all, var_q_all)
    zfv_all = _cmh_z_from_num_var(num_fv_all, var_fv_all)
    zfu_all = _cmh_z_from_num_var(num_fu_all, var_fu_all)

    # Observed directional 3-way Berger IUT score
    s_obs = vecs.obs_sign
    z_iut_obs = np.minimum(
        s_obs * zq_all[:, 0],
        np.minimum(s_obs * zfv_all[:, 0], s_obs * zfu_all[:, 0]),
    )
    t_obs_pos = np.maximum(z_iut_obs, 0.0)

    # Two-sided directional 3-way IUT score under B within-stratum permutations
    zq_perm = zq_all[:, 1:]
    zfv_perm = zfv_all[:, 1:]
    zfu_perm = zfu_all[:, 1:]
    t_perm_r = np.minimum(zq_perm, np.minimum(zfv_perm, zfu_perm))
    t_perm_nr = np.minimum(-zq_perm, np.minimum(-zfv_perm, -zfu_perm))
    t_perm = np.maximum(0.0, np.maximum(t_perm_r, t_perm_nr))

    exceed_pair = np.sum(t_perm >= (np.expand_dims(t_obs_pos, axis=1) - 1e-12), axis=1)
    p_perm_iut = (1.0 + exceed_pair.astype(np.float64)) / float(1 + n_permutations)

    sw_weights = (
        candidates_df["spatial_weight"].cast(pl.Float64).to_numpy()
        if "spatial_weight" in candidates_df.columns
        else np.ones(n_cand, dtype=np.float64)
    )
    fdr_perm_iut = compute_weighted_bh_fdr(p_perm_iut, sw_weights)

    max_t_per_perm = np.max(t_perm, axis=0)
    exceed_fwer = np.sum(
        np.expand_dims(max_t_per_perm, axis=0) >= (np.expand_dims(t_obs_pos, axis=1) - 1e-12),
        axis=1,
    )
    fwer_perm_iut = (1.0 + exceed_fwer.astype(np.float64)) / float(1 + n_permutations)

    # Exact O(1)-subtraction Leave-One-Donor-Out (LODO) jackknife across all informative donors
    tot_nq = num_q_all[:, 0]
    tot_vq = var_q_all[:, 0]
    tot_nfv = num_fv_all[:, 0]
    tot_vfv = var_fv_all[:, 0]
    tot_nfu = num_fu_all[:, 0]
    tot_vfu = var_fu_all[:, 0]

    lodo_z_list: list[FloatArray] = []
    for s_idx, grp in enumerate(informative_groups):
        yg = (y_arr[grp] > 0.5).astype(np.float64)
        xq_g = vecs.x_q[:, grp]
        xfv_g = vecs.x_fv[:, grp]
        xfu_g = vecs.x_fu[:, grp]
        mtot_g = vecs.m_tot[:, grp]

        sum_xq = np.sum(xq_g, axis=1, keepdims=True)
        sum_mtot = np.sum(mtot_g, axis=1, keepdims=True)
        sum_xfv = np.sum(xfv_g, axis=1, keepdims=True)
        sum_xfu = np.sum(xfu_g, axis=1, keepdims=True)

        sum_aq = np.sum(xq_g * yg, axis=1, keepdims=True)
        sum_rmtot = np.sum(mtot_g * yg, axis=1, keepdims=True)
        sum_cfv = np.sum(xfv_g * yg, axis=1, keepdims=True)
        sum_cfu = np.sum(xfu_g * yg, axis=1, keepdims=True)

        # Hold out each donor j in grp simultaneously as columns of shape (C, g_len)
        yg_row = np.expand_dims(yg, axis=0)
        r1_q_lo = np.maximum(sum_xq - xq_g, 0.0)
        mtot_lo = np.maximum(sum_mtot - mtot_g, 0.0)
        r0_notq_lo = np.maximum(mtot_lo - r1_q_lo, 0.0)
        r0_fv_lo = np.maximum(sum_xfv - xfv_g, 0.0)
        r0_fu_lo = np.maximum(sum_xfu - xfu_g, 0.0)

        aq_lo = np.maximum(sum_aq - xq_g * yg_row, 0.0)
        c_notq_lo = np.maximum((sum_rmtot - mtot_g * yg_row) - aq_lo, 0.0)
        cfv_lo = np.maximum(sum_cfv - xfv_g * yg_row, 0.0)
        cfu_lo = np.maximum(sum_cfu - xfu_g * yg_row, 0.0)

        # Only evaluate columns where leaving donor j out still leaves >= 1 R and >= 1 NR in stratum
        n_r_g = float(np.sum(yg))
        n_nr_g = float(grp.shape[0]) - n_r_g
        rem_r = n_r_g - yg
        rem_nr = n_nr_g - (1.0 - yg)
        stratum_still_inf = (rem_r >= 1.0) & (rem_nr >= 1.0)
        inf_row = np.expand_dims(stratum_still_inf, axis=0)

        nq_lo_raw, vq_lo_raw = _stratum_cmh_components_2d(aq_lo, r1_q_lo, c_notq_lo, r0_notq_lo)
        nfv_lo_raw, vfv_lo_raw = _stratum_cmh_components_2d(aq_lo, r1_q_lo, cfv_lo, r0_fv_lo)
        nfu_lo_raw, vfu_lo_raw = _stratum_cmh_components_2d(aq_lo, r1_q_lo, cfu_lo, r0_fu_lo)

        nq_lo = np.where(inf_row, nq_lo_raw, 0.0)
        vq_lo = np.where(inf_row, vq_lo_raw, 0.0)
        nfv_lo = np.where(inf_row, nfv_lo_raw, 0.0)
        vfv_lo = np.where(inf_row, vfv_lo_raw, 0.0)
        nfu_lo = np.where(inf_row, nfu_lo_raw, 0.0)
        vfu_lo = np.where(inf_row, vfu_lo_raw, 0.0)

        base_nq = np.expand_dims(tot_nq - obs_num_q_s[s_idx, :], axis=1) + nq_lo
        base_vq = np.expand_dims(tot_vq - obs_var_q_s[s_idx, :], axis=1) + vq_lo
        base_nfv = np.expand_dims(tot_nfv - obs_num_fv_s[s_idx, :], axis=1) + nfv_lo
        base_vfv = np.expand_dims(tot_vfv - obs_var_fv_s[s_idx, :], axis=1) + vfv_lo
        base_nfu = np.expand_dims(tot_nfu - obs_num_fu_s[s_idx, :], axis=1) + nfu_lo
        base_vfu = np.expand_dims(tot_vfu - obs_var_fu_s[s_idx, :], axis=1) + vfu_lo

        s_col = np.expand_dims(s_obs, axis=1)
        zq_j = _cmh_z_from_num_var(base_nq, base_vq)
        zfv_j = _cmh_z_from_num_var(base_nfv, base_vfv)
        zfu_j = _cmh_z_from_num_var(base_nfu, base_vfu)
        lodo_z_list.append(np.minimum(s_col * zq_j, np.minimum(s_col * zfv_j, s_col * zfu_j)))

    lodo_z_mat = np.concatenate(lodo_z_list, axis=1)
    lodo_min_z = np.min(lodo_z_mat, axis=1)
    lodo_sign_cons = np.mean(lodo_z_mat > 0.0, axis=1)

    # Leave-One-Stratum-Out (LOSO) jackknife when >= 2 informative strata exist
    if k_strata >= 2:
        loso_cols: list[FloatArray] = []
        for s_idx in range(k_strata):
            zq_s = _cmh_z_from_num_var(tot_nq - obs_num_q_s[s_idx, :], tot_vq - obs_var_q_s[s_idx, :])
            zfv_s = _cmh_z_from_num_var(tot_nfv - obs_num_fv_s[s_idx, :], tot_vfv - obs_var_fv_s[s_idx, :])
            zfu_s = _cmh_z_from_num_var(tot_nfu - obs_num_fu_s[s_idx, :], tot_vfu - obs_var_fu_s[s_idx, :])
            loso_cols.append(np.minimum(s_obs * zq_s, np.minimum(s_obs * zfv_s, s_obs * zfu_s)))
        loso_mat = np.column_stack(loso_cols)
        loso_min_z = np.min(loso_mat, axis=1)
        loso_sign_cons = np.mean(loso_mat > 0.0, axis=1)
    else:
        loso_min_z = lodo_min_z
        loso_sign_cons = lodo_sign_cons

    base_sig = (
        candidates_df["is_significant"].to_numpy()
        if "is_significant" in candidates_df.columns
        else (candidates_df["fdr_iut"].cast(pl.Float64).to_numpy() < fdr_threshold)
    )
    is_verified = (
        base_sig
        & (fdr_perm_iut < fdr_threshold)
        & (lodo_min_z >= jackknife_z_floor)
        & (loso_sign_cons >= 0.80 - 1e-12)
    )

    return Success(
        candidates_df.with_columns(
            pl.Series("z_iut_cmh", z_iut_obs, dtype=pl.Float64),
            pl.Series("pvalue_perm_iut", p_perm_iut, dtype=pl.Float64),
            pl.Series("fdr_perm_iut", fdr_perm_iut, dtype=pl.Float64),
            pl.Series("fwer_perm_iut", fwer_perm_iut, dtype=pl.Float64),
            pl.Series("lodo_min_z_iut", lodo_min_z, dtype=pl.Float64),
            pl.Series("lodo_sign_consistency", lodo_sign_cons, dtype=pl.Float64),
            pl.Series("loso_min_z_iut", loso_min_z, dtype=pl.Float64),
            pl.Series("loso_sign_consistency", loso_sign_cons, dtype=pl.Float64),
            pl.Series("is_statistically_verified", is_verified, dtype=pl.Boolean),
        )
    )
