"""Unit and Hypothesis property-based tests for milopy.pairwise.

Covers:
1. Algebraic & statistical invariants via Hypothesis property-based testing:
   - Quadrant swap involution and Hamming-1 flip equivariance
   - Outer-product probability conservation across all PresenceModes
   - Weighted Benjamini-Hochberg reduction to unweighted BH and bounds
   - Exact pair-exchange symmetry under neighborhood index reversal
   - Exact donor library-size scale invariance of compositional log-ratio balances
2. Canonical tumor-immunology synthetic benchmarks:
   - `++` synergistic co-presence (e.g., CXCL13+ CD8 T + LAMP3+ DC in responders)
   - `+-` / `-+` gated suppression where univariate Milo signal cancels to zero
   - `--` dual-absence immune desert (response OR-gate)
   - Adversarial singleton piggybacking null (blocked by Berger's IUT)
   - Adversarial capture-depth confounding null (blocked by detection weights + Firth log-L LRT)
   - Ubiquitous continuous log-ratio balance (CXCL9:SPP1 polarity)
   - Graph cell-sharing exclusion, 2D pair module deduplication, multi-cohort CMH stratification,
     permutation null calibration, dependency injection, and AnnData immutability.
"""

from __future__ import annotations

import anndata as ad
from hypothesis import given, settings, strategies as st
import numpy as np
import pandas as pd
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Result, Success
from scipy import sparse as sp
from scipy.stats import norm

import milopy
from milopy import (
    ALL_QUADRANT_STATES,
    ContinuousAxis,
    EnrichmentDirection,
    OverlapConfig,
    PairModuleConfig,
    PairScope,
    PairScopeFilter,
    PairwiseAnalysisResult,
    PairwiseConfig,
    PairwiseMethod,
    PrefilterConfig,
    PrefilterMode,
    PresenceConfig,
    PresenceMode,
    PresenceRepresentation,
    QuadrantState,
    RefinementConfig,
    ScreenConfig,
    adjacent_flips,
    analyze_nhood_pairs,
    build_external_prior_weight_matrix,
    coarse_grain_nhoods_by_graph,
    compute_benjamini_hochberg_fdr,
    compute_tarone_max_iut_z_strata,
    compute_weighted_bh_fdr,
    estimate_nhood_presence,
    group_nhood_pairs,
    group_nhoods,
    project_pair_to_cells,
    swap_quadrant,
    test_nhood_pairs_adata,
)


def _disjoint_nhoods_matrix(n_nhoods: int, cells_per_nhood: int = 5) -> sp.csc_matrix:
    """Builds a sparse incidence matrix with strictly disjoint neighborhoods."""
    n_cells = n_nhoods * cells_per_nhood
    rows = np.arange(n_cells, dtype=np.int64)
    cols = np.repeat(np.arange(n_nhoods, dtype=np.int64), cells_per_nhood)
    return sp.csc_matrix((np.ones(n_cells, dtype=np.float64), (rows, cols)), shape=(n_cells, n_nhoods))


# ---------------------------------------------------------------------------
# 1. Hypothesis Property-Based Invariant Tests
# ---------------------------------------------------------------------------


def test_quadrant_swap_involution_and_flip_equivariance() -> None:
    for q in ALL_QUADRANT_STATES:
        assert swap_quadrant(swap_quadrant(q)) == q
        flip_v, flip_u = adjacent_flips(q)
        swapped_fv, swapped_fu = adjacent_flips(swap_quadrant(q))
        assert swapped_fv == swap_quadrant(flip_u)
        assert swapped_fu == swap_quadrant(flip_v)


@settings(max_examples=35, deadline=None)
@given(
    n_nhoods=st.integers(min_value=3, max_value=8),
    n_samples=st.integers(min_value=6, max_value=16),
    mode=st.sampled_from(list(PresenceMode)),
    use_detect=st.booleans(),
    seed=st.integers(min_value=0, max_value=10000),
)
def test_prop_presence_outer_product_conservation(
    n_nhoods: int,
    n_samples: int,
    mode: PresenceMode,
    use_detect: bool,
    seed: int,
) -> None:
    rng = np.random.default_rng(seed)
    counts = rng.poisson(lam=2.5, size=(n_nhoods, n_samples)).astype(np.float64)
    # Ensure every neighborhood is active in >= 2 donors
    counts[:, :2] = np.maximum(counts[:, :2], 2.0)
    libs = rng.uniform(50.0, 1000.0, size=n_samples)

    res = estimate_nhood_presence(
        count_matrix=counts,
        library_sizes=Some(libs),
        config=PresenceConfig(
            mode=mode,
            min_active_donors=2,
            use_detection_power_weights=use_detect,
        ),
    )
    assert isinstance(res, Success)
    rep = res.unwrap()

    p, q, m = rep.presence_prob, rep.absence_prob, rep.observability
    assert np.all((p >= 0.0) & (p <= 1.0))
    assert np.all((q >= 0.0) & (q <= 1.0))
    assert np.all((m >= 0.0) & (m <= 1.0))

    # Exact outer-product additivity across all pairs and donors
    quad_sum = (
        p[0, :] * p[1, :]
        + p[0, :] * q[1, :]
        + q[0, :] * p[1, :]
        + q[0, :] * q[1, :]
    )
    np.testing.assert_allclose(quad_sum, m[0, :] * m[1, :], atol=1e-12)
    if not use_detect:
        np.testing.assert_allclose(quad_sum, np.ones(n_samples), atol=1e-12)


@settings(max_examples=40, deadline=None)
@given(
    p_list=st.lists(st.floats(min_value=0.0, max_value=1.0), min_size=1, max_size=40),
    w_scale=st.floats(min_value=0.1, max_value=25.0),
)
def test_prop_weighted_bh_matches_unweighted_on_uniform_weights(
    p_list: list[float],
    w_scale: float,
) -> None:
    p_arr = np.asarray(p_list, dtype=np.float64)
    uniform_w = np.full(p_arr.shape[0], w_scale, dtype=np.float64)

    q_weighted = compute_weighted_bh_fdr(p_arr, uniform_w)
    q_standard = compute_benjamini_hochberg_fdr(p_arr)
    assert np.all((q_weighted >= 0.0) & (q_weighted <= 1.0))
    np.testing.assert_allclose(q_weighted, q_standard, atol=1e-12)


@settings(max_examples=20, deadline=None)
@given(seed=st.integers(min_value=1, max_value=5000))
def test_prop_pair_exchange_symmetry(seed: int) -> None:
    rng = np.random.default_rng(seed)
    n_nhoods, n_samples = 5, 16
    counts = rng.integers(0, 8, size=(n_nhoods, n_samples)).astype(np.float64)
    counts[:, :3] = np.maximum(counts[:, :3], 1.0)
    libs = rng.uniform(100.0, 500.0, size=n_samples)
    labels = ["responder"] * 8 + ["non-responder"] * 8

    nhoods_fwd = _disjoint_nhoods_matrix(n_nhoods)
    cfg = PairwiseConfig(
        method=PairwiseMethod.QUADRANT_2X2,
        refinement=RefinementConfig(enabled=False),
        modules=PairModuleConfig(enabled=False),
        screen=ScreenConfig(min_effective_samples=2.0),
    )

    res_fwd = analyze_nhood_pairs(
        nhoods_mat=nhoods_fwd,
        count_matrix=counts,
        response_labels=labels,
        sample_library_sizes=Some(tuple(libs.tolist())),
        config=cfg,
    ).unwrap()

    rev_perm = np.arange(n_nhoods - 1, -1, -1, dtype=np.int64)
    res_rev = analyze_nhood_pairs(
        nhoods_mat=nhoods_fwd[:, rev_perm],
        count_matrix=counts[rev_perm, :],
        response_labels=labels,
        sample_library_sizes=Some(tuple(libs.tolist())),
        config=cfg,
    ).unwrap()

    fwd_map = {
        (int(r["nhood_u"]), int(r["nhood_v"]), str(r["quadrant"])): r
        for r in res_fwd.pairs_df.to_dicts()
    }
    rev_map = {
        (int(r["nhood_u"]), int(r["nhood_v"]), str(r["quadrant"])): r
        for r in res_rev.pairs_df.to_dicts()
    }

    assert len(fwd_map) == len(rev_map)
    for (u, v, q_str), row_f in fwd_map.items():
        u_rev = n_nhoods - 1 - v
        v_rev = n_nhoods - 1 - u
        q_rev = swap_quadrant(QuadrantState(q_str)).value
        row_r = rev_map[(u_rev, v_rev, q_rev)]
        np.testing.assert_allclose(row_f["log_or_quadrant"], row_r["log_or_quadrant"], atol=1e-11)
        np.testing.assert_allclose(row_f["pvalue_iut"], row_r["pvalue_iut"], atol=1e-11)
        np.testing.assert_allclose(row_f["fdr_quadrant"], row_r["fdr_quadrant"], atol=1e-11)
        np.testing.assert_allclose(row_f["fdr_iut"], row_r["fdr_iut"], atol=1e-11)
        np.testing.assert_allclose(row_f["log_or_interaction"], row_r["log_or_interaction"], atol=1e-11)


@settings(max_examples=20, deadline=None)
@given(seed=st.integers(min_value=1, max_value=5000))
def test_prop_compositional_log_ratio_scale_invariance(seed: int) -> None:
    rng = np.random.default_rng(seed)
    n_nhoods, n_samples = 4, 18
    counts = rng.uniform(1.0, 50.0, size=(n_nhoods, n_samples))
    libs = rng.uniform(200.0, 1500.0, size=n_samples)
    multipliers = rng.uniform(0.2, 20.0, size=n_samples)
    labels = ["responder"] * 9 + ["non-responder"] * 9
    nhoods_mat = _disjoint_nhoods_matrix(n_nhoods)

    cfg = PairwiseConfig(method=PairwiseMethod.LOG_RATIO_BALANCE)
    res_base = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=labels,
        sample_library_sizes=Some(tuple(libs.tolist())),
        config=cfg,
    ).unwrap()

    res_scaled = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts * np.expand_dims(multipliers, axis=0),
        response_labels=labels,
        sample_library_sizes=Some(tuple((libs * multipliers).tolist())),
        config=cfg,
    ).unwrap()

    df_b = res_base.continuous_pairs_df.sort(["nhood_u", "nhood_v"])
    df_s = res_scaled.continuous_pairs_df.sort(["nhood_u", "nhood_v"])
    np.testing.assert_allclose(df_b["effect_size"].to_numpy(), df_s["effect_size"].to_numpy(), atol=1e-10)
    np.testing.assert_allclose(df_b["pvalue_axis"].to_numpy(), df_s["pvalue_axis"].to_numpy(), atol=1e-10)
    np.testing.assert_allclose(df_b["pvalue_iut"].to_numpy(), df_s["pvalue_iut"].to_numpy(), atol=1e-10)


# ---------------------------------------------------------------------------
# 2. Canonical Immunological Regime Benchmarks
# ---------------------------------------------------------------------------


def test_regime_both_present_synergy_pp() -> None:
    """Regime 1 (++): Response occurs predominantly when BOTH nhood 0 and nhood 1 are present."""
    rng = np.random.default_rng(101)
    n_nhoods, n_samples = 10, 80
    nhoods_mat = _disjoint_nhoods_matrix(n_nhoods)
    libs = np.full(n_samples, 1000.0)

    # 20 donors in each quadrant of (0, 1): ++, +-, -+, --
    u_pres = np.array([1] * 20 + [1] * 20 + [0] * 20 + [0] * 20, dtype=np.float64)
    v_pres = np.array([1] * 20 + [0] * 20 + [1] * 20 + [0] * 20, dtype=np.float64)

    # ++ has 18/20 R (90%); +-, -+, -- each have 3/20 R (15%)
    y_resp = (
        ["responder"] * 18 + ["non-responder"] * 2
        + ["responder"] * 3 + ["non-responder"] * 17
        + ["responder"] * 3 + ["non-responder"] * 17
        + ["responder"] * 3 + ["non-responder"] * 17
    )

    counts = np.zeros((n_nhoods, n_samples), dtype=np.float64)
    counts[0, :] = u_pres * 25.0
    counts[1, :] = v_pres * 25.0
    for k in range(2, n_nhoods):
        counts[k, :] = (rng.random(n_samples) > 0.5).astype(np.float64) * 20.0

    res = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=y_resp,
        sample_library_sizes=Some(tuple(libs.tolist())),
        config=PairwiseConfig(method=PairwiseMethod.QUADRANT_2X2),
    ).unwrap()

    sig = res.pairs_df.filter(pl.col("is_significant"))
    target = sig.filter(
        (pl.col("nhood_u") == 0)
        & (pl.col("nhood_v") == 1)
        & (pl.col("quadrant") == QuadrantState.BOTH_PRESENT.value)
    )
    assert target.height == 1
    assert target["direction"][0] == EnrichmentDirection.RESPONDER.value
    assert float(target["log_or_interaction"][0]) > 1.5
    assert float(target["fdr_iut"][0]) < 0.05


def test_regime_present_missing_gating_pm_zero_marginal() -> None:
    """Regime 2 (+-): Nhood 0 (CD8) has ZERO marginal association with response because Treg (Nhood 1)
    cancels it in ++, yet (CD8+, Treg-) is strongly responder-enriched!
    """
    rng = np.random.default_rng(202)
    n_nhoods, n_samples = 8, 140
    nhoods_mat = _disjoint_nhoods_matrix(n_nhoods)
    libs = np.full(n_samples, 1000.0)

    # 35 donors per quadrant (++, +-, -+, --)
    u_pres = np.array([1] * 35 + [1] * 35 + [0] * 35 + [0] * 35, dtype=np.float64)
    v_pres = np.array([1] * 35 + [0] * 35 + [1] * 35 + [0] * 35, dtype=np.float64)

    # Marginal check on u: ++ has 2 R, +- has 33 R -> u+ has 35/70 R (50%!).
    # -+ has 17 R, -- has 18 R -> u- has 35/70 R (50%!).
    # Univariate u has EXACTLY 50% R vs 50% NR (zero marginal effect!).
    y_resp = (
        ["responder"] * 2 + ["non-responder"] * 33
        + ["responder"] * 33 + ["non-responder"] * 2
        + ["responder"] * 17 + ["non-responder"] * 18
        + ["responder"] * 18 + ["non-responder"] * 17
    )

    counts = np.zeros((n_nhoods, n_samples), dtype=np.float64)
    counts[0, :] = u_pres * 30.0
    counts[1, :] = v_pres * 30.0
    for k in range(2, n_nhoods):
        counts[k, :] = (rng.random(n_samples) > 0.5).astype(np.float64) * 20.0

    # Verify univariate u frequency in R vs NR is identical (35 vs 35)
    r_mask = np.array([lbl == "responder" for lbl in y_resp])
    assert int(np.sum(u_pres[r_mask])) == int(np.sum(u_pres[~r_mask])) == 35

    res = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=y_resp,
        sample_library_sizes=Some(tuple(libs.tolist())),
        config=PairwiseConfig(method=PairwiseMethod.QUADRANT_2X2),
    ).unwrap()

    sig = res.pairs_df.filter(pl.col("is_significant"))
    pm_hit = sig.filter(
        (pl.col("nhood_u") == 0)
        & (pl.col("nhood_v") == 1)
        & (pl.col("quadrant") == QuadrantState.FIRST_ONLY.value)
    )
    assert pm_hit.height == 1
    assert pm_hit["direction"][0] == EnrichmentDirection.RESPONDER.value
    assert float(pm_hit["fdr_iut"][0]) < 0.05


def test_regime_both_missing_immune_desert_mm() -> None:
    """Regime 3 (--): Having either nhood 0 or nhood 1 yields response (OR-gate), whereas
    missing BOTH (--) causes non-response (immune desert).
    """
    rng = np.random.default_rng(303)
    n_nhoods, n_samples = 8, 80
    nhoods_mat = _disjoint_nhoods_matrix(n_nhoods)
    libs = np.full(n_samples, 1000.0)

    u_pres = np.array([1] * 20 + [1] * 20 + [0] * 20 + [0] * 20, dtype=np.float64)
    v_pres = np.array([1] * 20 + [0] * 20 + [1] * 20 + [0] * 20, dtype=np.float64)

    # ++, +-, -+ each have 17/20 R (85%); -- has only 1/20 R (95% NR!)
    y_resp = (
        ["responder"] * 17 + ["non-responder"] * 3
        + ["responder"] * 17 + ["non-responder"] * 3
        + ["responder"] * 17 + ["non-responder"] * 3
        + ["responder"] * 1 + ["non-responder"] * 19
    )

    counts = np.zeros((n_nhoods, n_samples), dtype=np.float64)
    counts[0, :] = u_pres * 25.0
    counts[1, :] = v_pres * 25.0
    for k in range(2, n_nhoods):
        counts[k, :] = (rng.random(n_samples) > 0.5).astype(np.float64) * 20.0

    res = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=y_resp,
        sample_library_sizes=Some(tuple(libs.tolist())),
        config=PairwiseConfig(method=PairwiseMethod.QUADRANT_2X2),
    ).unwrap()

    mm_hit = res.pairs_df.filter(
        pl.col("is_significant")
        & (pl.col("nhood_u") == 0)
        & (pl.col("nhood_v") == 1)
        & (pl.col("quadrant") == QuadrantState.BOTH_MISSING.value)
    )
    assert mm_hit.height == 1
    assert mm_hit["direction"][0] == EnrichmentDirection.NON_RESPONDER.value
    assert float(mm_hit["fdr_iut"][0]) < 0.05


def test_adversarial_singleton_piggybacking_blocked_by_iut() -> None:
    """Regime 4: A single massive marginal responder neighborhood paired with 15 noise
    neighborhoods triggers raw quadrant enrichment, but Berger's IUT flags 0 false positives.
    """
    rng = np.random.default_rng(404)
    n_nhoods, n_samples = 16, 60
    nhoods_mat = _disjoint_nhoods_matrix(n_nhoods)
    libs = np.full(n_samples, 1000.0)
    y_resp = ["responder"] * 30 + ["non-responder"] * 30

    counts = np.zeros((n_nhoods, n_samples), dtype=np.float64)
    # Nhood 0: present in 29/30 R and 1/30 NR
    counts[0, :29] = 30.0
    counts[0, 30] = 30.0
    # Nhoods 1..15: pure independent coin flips (exactly balanced 15/30 in R and NR or random)
    for k in range(1, n_nhoods):
        counts[k, :] = (rng.random(n_samples) > 0.5).astype(np.float64) * 20.0

    res = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=y_resp,
        sample_library_sizes=Some(tuple(libs.tolist())),
        config=PairwiseConfig(
            method=PairwiseMethod.QUADRANT_2X2,
            screen=ScreenConfig(fdr_threshold=0.05),
        ),
    ).unwrap()

    # Raw quadrant FDR sees piggybacking on nhood 0, whereas IUT rejects piggybacking cleanly
    raw_quad_hits = res.pairs_df.filter(
        (pl.col("nhood_u") == 0) & (pl.col("fdr_quadrant") < 0.05)
    )
    assert raw_quad_hits.height > 0
    assert res.n_significant_quadrant_pairs == 0


def test_adversarial_capture_depth_confounding_null() -> None:
    """Regime 5: Non-responders have 12x shallower sequencing depth than responders under
    identical biological cell-state frequencies. Depth weighting + Firth log(L) adjustment
    must prevent false `--` and `++` discoveries.
    """
    rng = np.random.default_rng(505)
    n_nhoods, n_samples = 10, 60
    nhoods_mat = _disjoint_nhoods_matrix(n_nhoods)
    y_resp = ["responder"] * 30 + ["non-responder"] * 30

    # Responders: 500 cells/biopsy; Non-responders: 35 cells/biopsy
    libs = np.array([500.0] * 30 + [35.0] * 30, dtype=np.float64)
    true_rate = 0.008  # Expected count: 4.0 in R vs 0.28 in NR (75% technical zeros in NR!)
    counts = rng.poisson(lam=true_rate * np.expand_dims(libs, axis=0), size=(n_nhoods, n_samples)).astype(
        np.float64
    )

    res = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=y_resp,
        sample_library_sizes=Some(tuple(libs.tolist())),
        config=PairwiseConfig(
            method=PairwiseMethod.QUADRANT_2X2,
            presence=PresenceConfig(
                mode=PresenceMode.EMPIRICAL_BAYES_POSTERIOR,
                use_detection_power_weights=True,
            ),
            refinement=RefinementConfig(enabled=True, adjust_library_size=True),
        ),
    ).unwrap()

    assert res.n_significant_quadrant_pairs == 0


def test_ubiquitous_continuous_log_ratio_balance() -> None:
    """Regime 6: All neighborhoods have positive counts in 100% of donors; a collinear
    polarity ratio (nhood 0 / nhood 1) drives response while shared baseline infiltrate cancels.
    """
    rng = np.random.default_rng(606)
    n_nhoods, n_samples = 8, 60
    nhoods_mat = _disjoint_nhoods_matrix(n_nhoods)
    libs = rng.uniform(800.0, 2000.0, size=n_samples)
    y_bin = np.array([1.0] * 30 + [0.0] * 30)
    y_resp = ["responder"] * 30 + ["non-responder"] * 30

    shared_infiltrate = rng.normal(0.0, 1.2, size=n_samples)
    rates = np.zeros((n_nhoods, n_samples), dtype=np.float64)
    # Nhood 0 and Nhood 1 share heavy common variance, but log2(r0 / r1) tracks y_bin!
    rates[0, :] = 0.05 * np.exp2(shared_infiltrate + 0.9 * y_bin + rng.normal(0.0, 0.2, n_samples))
    rates[1, :] = 0.05 * np.exp2(shared_infiltrate - 0.9 * y_bin + rng.normal(0.0, 0.2, n_samples))
    for k in range(2, n_nhoods):
        rates[k, :] = 0.04 * np.exp2(rng.normal(0.0, 0.6, n_samples))

    counts = np.maximum(np.round(rates * np.expand_dims(libs, axis=0)), 5.0)
    assert np.all(counts > 0.0)

    res = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=y_resp,
        sample_library_sizes=Some(tuple(libs.tolist())),
        config=PairwiseConfig(method=PairwiseMethod.LOG_RATIO_BALANCE),
    ).unwrap()

    sig = res.continuous_pairs_df.filter(pl.col("is_significant"))
    bal_hit = sig.filter(
        (pl.col("nhood_u") == 0)
        & (pl.col("nhood_v") == 1)
        & (pl.col("axis") == ContinuousAxis.LOG_RATIO_BALANCE.value)
    )
    assert bal_hit.height == 1
    assert bal_hit["direction"][0] == EnrichmentDirection.RESPONDER.value


def test_overlap_exclusion_module_grouping_stratification_and_adata() -> None:
    """End-to-end AnnData integration test verifying:
    - Cell-sharing neighborhoods within a lineage are excluded at max_jaccard_overlap=0.0
    - Multiple overlapping neighborhoods across two lineages collapse into 1 Pair Module
    - Multi-cohort dataset_col stratification, permutations, projection, and 1D group_nhoods work
    - Input AnnData is strictly unmutated.
    """
    # 3 overlapping CD8 neighborhoods (cells 0..19) and 3 overlapping DC neighborhoods (cells 20..39)
    # plus 2 disjoint B-cell noise neighborhoods (cells 40..59) across 80 patients
    n_patients = 80
    cells_per_pt = 20
    n_cells = n_patients * cells_per_pt

    rng = np.random.default_rng(707)
    pt_ids = np.repeat([f"PT_{i:02d}" for i in range(n_patients)], cells_per_pt)
    datasets = np.repeat(["CohortA" if i % 2 == 0 else "CohortB" for i in range(n_patients)], cells_per_pt)

    # Quadrant assignment per patient across CD8 and DC lineages
    u_has = np.array([1] * 20 + [1] * 20 + [0] * 20 + [0] * 20)
    v_has = np.array([1] * 20 + [0] * 20 + [1] * 20 + [0] * 20)
    pt_resp = (
        ["responder"] * 18 + ["non-responder"] * 2
        + ["responder"] * 2 + ["non-responder"] * 18
        + ["responder"] * 2 + ["non-responder"] * 18
        + ["responder"] * 2 + ["non-responder"] * 18
    )
    cell_resp = np.repeat(pt_resp, cells_per_pt)

    cell_types: list[str] = []
    rows: list[int] = []
    cols: list[int] = []

    # Neighborhoods 0, 1, 2: overlapping CD8_T sub-spheres; 3, 4, 5: overlapping mregDC sub-spheres; 6, 7: B_cell
    for p_idx in range(n_patients):
        base = p_idx * cells_per_pt
        for offset in range(cells_per_pt):
            cid = base + offset
            if offset < 6 and u_has[p_idx] == 1:
                cell_types.append("CD8_T")
                if offset in (0, 1, 2, 3):
                    rows.append(cid)
                    cols.append(0)
                if offset in (1, 2, 3, 4):
                    rows.append(cid)
                    cols.append(1)
                if offset in (2, 3, 4, 5):
                    rows.append(cid)
                    cols.append(2)
            elif 6 <= offset < 12 and v_has[p_idx] == 1:
                cell_types.append("mregDC")
                if offset in (6, 7, 8, 9):
                    rows.append(cid)
                    cols.append(3)
                if offset in (7, 8, 9, 10):
                    rows.append(cid)
                    cols.append(4)
                if offset in (8, 9, 10, 11):
                    rows.append(cid)
                    cols.append(5)
            else:
                cell_types.append("B_cell")
                if offset == 12 and (rng.random() > 0.5):
                    rows.append(cid)
                    cols.append(6)
                elif offset == 13 and (rng.random() > 0.5):
                    rows.append(cid)
                    cols.append(7)

    nhoods_mat = sp.csc_matrix(
        (np.ones(len(rows), dtype=np.float64), (rows, cols)),
        shape=(n_cells, 8),
    )
    obs_df = pd.DataFrame(
        {
            "patient_id": pt_ids,
            "clinical_response": cell_resp,
            "dataset_id": datasets,
            "cell_type": cell_types,
        }
    )
    adata = ad.AnnData(X=sp.csr_matrix((n_cells, 4)), obs=obs_df)
    adata.obsm["nhoods"] = nhoods_mat

    uns_keys_before = set(adata.uns.keys())
    obs_cols_before = list(adata.obs.columns)

    cfg = PairwiseConfig(
        method=PairwiseMethod.QUADRANT_2X2,
        overlap=OverlapConfig(
            max_jaccard_overlap=0.0,
            scope_filter=PairScopeFilter.CROSS_CELL_TYPE_ONLY,
        ),
        screen=ScreenConfig(n_permutations=50, permutation_seed=42),
        modules=PairModuleConfig(enabled=True, min_endpoint_jaccard=0.20),
    )

    res_wrapped = test_nhood_pairs_adata(
        adata=adata,
        sample_col="patient_id",
        design_col="clinical_response",
        dataset_col=Some("dataset_id"),
        cell_type_col=Some("cell_type"),
        config=cfg,
    )
    assert isinstance(res_wrapped, Success)
    res: PairwiseAnalysisResult = res_wrapped.unwrap()

    # Verify strict immutability of input AnnData
    assert set(adata.uns.keys()) == uns_keys_before
    assert list(adata.obs.columns) == obs_cols_before

    # Overlapping within-lineage pairs (0,1), (1,2), (3,4), (4,5) must never appear
    evaluated_uv = set(
        zip(res.pairs_df["nhood_u"].to_list(), res.pairs_df["nhood_v"].to_list())
    )
    assert (0, 1) not in evaluated_uv
    assert (3, 4) not in evaluated_uv
    assert set(res.pairs_df["pair_scope"].to_list()) == {PairScope.CROSS_CELL_TYPE.value}

    # All 3 x 3 = 9 cross-lineage (CD8_T, mregDC) ++ pairs should be significant and collapse into 1 module!
    pp_mods = res.modules_df.filter(
        (pl.col("quadrant") == QuadrantState.BOTH_PRESENT.value)
        & (pl.col("direction") == EnrichmentDirection.RESPONDER.value)
    )
    assert pp_mods.height == 1
    assert int(pp_mods["n_member_pairs"][0]) == 9
    assert pp_mods["dominant_cell_type_u"][0] == "CD8_T"
    assert pp_mods["dominant_cell_type_v"][0] == "mregDC"

    # Also test 1D group_nhoods and cell-level pair projection
    grp_1d = group_nhoods(
        nhoods_mat=nhoods_mat,
        nhood_logfc=np.array([2.0, 2.1, 1.9, 1.8, 2.0, 1.7, 0.1, -0.1]),
        nhood_is_sig=np.array([True, True, True, True, True, True, False, False]),
        min_jaccard=0.20,
    ).unwrap()
    assert grp_1d[0] == grp_1d[1] == grp_1d[2] >= 0
    assert grp_1d[3] == grp_1d[4] == grp_1d[5] >= 0
    assert grp_1d[0] != grp_1d[3]
    assert grp_1d[6] == -1 and grp_1d[7] == -1

    proj = project_pair_to_cells(nhoods_mat, nhood_u=0, nhood_v=3).unwrap()
    assert proj.height == n_cells
    assert set(proj["pair_membership"].unique().to_list()).issubset(
        {"nhood_u", "nhood_v", "shared", "background"}
    )


def test_custom_refiner_dependency_injection_and_existing_exports() -> None:
    """Verifies Protocol dependency injection for CandidatePairRefiner and intact legacy exports."""
    for legacy_name in (
        "build_graph",
        "make_nhoods",
        "count_cells",
        "calc_nhood_distance",
        "test_nhoods",
        "test_nhoods_meta",
        "test_nhoods_mixed",
        "compute_score_permutation_null",
        "evaluate_replicate_prevalence",
        "compute_neighborhood_composition",
        "project_nhoods_to_cells",
    ):
        assert hasattr(milopy, legacy_name)

    called_flag: list[bool] = []

    def passthrough_refiner(
        candidates_df: pl.DataFrame,
        presence: PresenceRepresentation,
        response_binary: np.ndarray,
        strata_codes: Maybe[np.ndarray],
        config: RefinementConfig,
        fdr_threshold: float,
    ) -> Result[pl.DataFrame, str]:
        _ = (presence, response_binary, strata_codes, config, fdr_threshold)
        called_flag.append(True)
        return Success(
            candidates_df.with_columns(
                pl.col("log_or_quadrant").alias("log_or_quadrant_firth"),
                pl.col("pvalue_quadrant").alias("pvalue_quadrant_firth"),
                pl.col("pvalue_iut").alias("pvalue_iut_firth"),
                pl.lit(False).alias("is_refined"),
            )
        )

    nhoods_mat = _disjoint_nhoods_matrix(3)
    counts = np.array(
        [
            [10.0, 10.0, 10.0, 10.0, 0.0, 0.0, 0.0, 0.0],
            [10.0, 10.0, 10.0, 10.0, 0.0, 0.0, 0.0, 0.0],
            [5.0, 0.0, 5.0, 0.0, 5.0, 0.0, 5.0, 0.0],
        ]
    )
    labels = ["responder"] * 4 + ["non-responder"] * 4
    res = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=labels,
        refiner=Some(passthrough_refiner),
    )
    assert isinstance(res, Success)
    assert called_flag == [True]


@settings(max_examples=25, deadline=None)
@given(
    n_nhoods=st.integers(min_value=4, max_value=9),
    n_per_arm_stratum=st.integers(min_value=4, max_value=10),
    n_strata=st.integers(min_value=1, max_value=3),
    seed=st.integers(min_value=1, max_value=10000),
)
def test_prop_co_abundance_gram_matches_explicit_projection(
    n_nhoods: int,
    n_per_arm_stratum: int,
    n_strata: int,
    seed: int,
) -> None:
    """Proves exact numerical equivalence between BLAS Gram co-abundance and explicit (n_pairs, J) residuals."""
    rng = np.random.default_rng(seed)
    n_samples = 2 * n_per_arm_stratum * n_strata
    counts = rng.poisson(lam=4.0, size=(n_nhoods, n_samples)).astype(np.float64)
    counts[:, :2] = np.maximum(counts[:, :2], 1.0)
    libs = rng.uniform(150.0, 1200.0, size=n_samples)

    labels: list[str] = []
    strata_names: list[str] = []
    for k in range(n_strata):
        labels.extend(["responder"] * n_per_arm_stratum + ["non-responder"] * n_per_arm_stratum)
        strata_names.extend([f"S{k}"] * (2 * n_per_arm_stratum))

    nhoods_mat = _disjoint_nhoods_matrix(n_nhoods)
    strata_arg: Maybe[tuple[str, ...]] = Some(tuple(strata_names)) if n_strata > 1 else Nothing

    res = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=labels,
        sample_library_sizes=Some(tuple(libs.tolist())),
        strata_labels=strata_arg,
        config=PairwiseConfig(method=PairwiseMethod.CONTINUOUS_INTERACTION),
    ).unwrap()

    # Explicit reference calculation per pair (u, v)
    pres = estimate_nhood_presence(counts, Some(libs)).unwrap()
    y_bin = np.array([1.0 if s == "responder" else 0.0 for s in labels], dtype=np.float64)
    groups = [
        np.arange(k * 2 * n_per_arm_stratum, (k + 1) * 2 * n_per_arm_stratum, dtype=np.int64)
        for k in range(n_strata)
    ]
    x_c = np.zeros_like(pres.arcsin_sqrt_prop)
    y_c = np.zeros(n_samples, dtype=np.float64)
    offset = 0
    for grp in groups:
        glen = int(grp.shape[0])
        x_c[:, offset : offset + glen] = pres.arcsin_sqrt_prop[:, grp] - np.mean(
            pres.arcsin_sqrt_prop[:, grp], axis=1, keepdims=True
        )
        y_c[offset : offset + glen] = y_bin[grp] - float(np.mean(y_bin[grp]))
        offset += glen

    x_norms = np.linalg.norm(x_c, axis=1, keepdims=True)
    z_ref = np.where(x_norms > 1e-12, x_c / np.maximum(x_norms, 1e-12), 0.0)
    y_unit = y_c / float(np.linalg.norm(y_c))

    df_sorted = res.continuous_pairs_df.sort(["nhood_u", "nhood_v"])
    for row in df_sorted.to_dicts():
        u = int(row["nhood_u"])
        v = int(row["nhood_v"])
        xu = z_ref[u, :]
        xv = z_ref[v, :]
        w_raw = xu * xv * float(np.sqrt(float(n_samples)))
        w_c = w_raw - float(np.mean(w_raw))
        rho = float(np.clip(np.dot(xu, xv), -0.9999, 0.9999))
        denom = max(1.0 - rho * rho, 1e-8)
        pu = float(np.dot(w_c, xu))
        pv = float(np.dot(w_c, xv))
        gu = (pu - rho * pv) / denom
        gv = (pv - rho * pu) / denom
        w_orth = w_c - gu * xu - gv * xv
        onorm = float(np.linalg.norm(w_orth))
        ru = float(np.dot(xu, y_unit))
        rv = float(np.dot(xv, y_unit))
        r2_m = float(np.clip(ru * ((ru - rho * rv) / denom) + rv * ((rv - rho * ru) / denom), 0.0, 0.99999))
        ref_eff = (
            float(np.clip((float(np.dot(w_orth / onorm, y_unit))) / np.sqrt(max(1.0 - r2_m, 1e-12)), -0.999999, 0.999999))
            if onorm > 1e-10
            else 0.0
        )
        np.testing.assert_allclose(float(row["effect_size"]), ref_eff, atol=1e-10)


def test_cmh_score_recovers_zero_cell_separation_under_high_multiplicity() -> None:
    """Verifies CMH score + Fisher-ancillary IUT recovers zero-cell separated pairs across >10,000 tests."""
    rng = np.random.default_rng(808)
    n_nhoods, n_samples = 80, 76  # 3,160 disjoint pairs = 12,640 quadrant hypotheses
    nhoods_mat = _disjoint_nhoods_matrix(n_nhoods, cells_per_nhood=3)
    libs = np.full(n_samples, 500.0)

    # Zero-cell separated ++ synergy: 22/24 R in ++ vs 0/24 R in +-, -+, and 2/28 R in --
    u_pres = np.array([1] * 24 + [1] * 12 + [0] * 12 + [0] * 28, dtype=np.float64)
    v_pres = np.array([1] * 24 + [0] * 12 + [1] * 12 + [0] * 28, dtype=np.float64)
    y_resp = (
        ["responder"] * 22 + ["non-responder"] * 2
        + ["non-responder"] * 12
        + ["non-responder"] * 12
        + ["responder"] * 2 + ["non-responder"] * 26
    )

    counts = np.zeros((n_nhoods, n_samples), dtype=np.float64)
    counts[0, :] = u_pres * 15.0
    counts[1, :] = v_pres * 15.0
    for k in range(2, n_nhoods):
        counts[k, :] = (rng.random(n_samples) < 0.35).astype(np.float64) * 12.0

    res = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=y_resp,
        sample_library_sizes=Some(tuple(libs.tolist())),
        config=PairwiseConfig(method=PairwiseMethod.QUADRANT_2X2),
    ).unwrap()

    target = res.pairs_df.filter(
        pl.col("is_significant")
        & (pl.col("nhood_u") == 0)
        & (pl.col("nhood_v") == 1)
        & (pl.col("quadrant") == QuadrantState.BOTH_PRESENT.value)
    )
    assert target.height == 1
    assert target["direction"][0] == EnrichmentDirection.RESPONDER.value
    assert float(target["fdr_iut"][0]) < 0.05


def test_multi_singleton_piggybacking_ancillary_iut_control() -> None:
    """Verifies Fisher-ancillary Stage-2 IUT blocks false positives when multiple singletons pair with noise."""
    rng = np.random.default_rng(910)
    n_nhoods, n_samples = 40, 80
    nhoods_mat = _disjoint_nhoods_matrix(n_nhoods, cells_per_nhood=3)
    libs = np.full(n_samples, 600.0)
    y_resp = ["responder"] * 40 + ["non-responder"] * 40

    counts = np.zeros((n_nhoods, n_samples), dtype=np.float64)
    # Nhood 0 is a strong R singleton; Nhood 1 is a strong NR singleton; 2..39 are pure coin-flip noise
    counts[0, :36] = 20.0
    counts[0, 40:43] = 20.0
    counts[1, :3] = 20.0
    counts[1, 40:76] = 20.0
    for k in range(2, n_nhoods):
        counts[k, :] = (rng.random(n_samples) < 0.40).astype(np.float64) * 15.0

    res = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=y_resp,
        sample_library_sizes=Some(tuple(libs.tolist())),
        config=PairwiseConfig(
            method=PairwiseMethod.HYBRID_ALL,
            screen=ScreenConfig(fdr_threshold=0.05),
        ),
    ).unwrap()

    # Neither singleton 0 nor singleton 1 may piggyback onto noise partners 2..39
    piggy_quad = res.pairs_df.filter(
        pl.col("is_significant") & (pl.col("nhood_u").is_in([0, 1]) & (pl.col("nhood_v") >= 2))
    )
    piggy_cont = res.continuous_pairs_df.filter(
        pl.col("is_significant") & (pl.col("nhood_u").is_in([0, 1]) & (pl.col("nhood_v") >= 2))
    )
    assert piggy_quad.height == 0
    assert piggy_cont.height == 0
    # Whereas the true opposite-arm pair (0, 1) itself IS genuinely +- / -+ and log-ratio enriched!
    true_pair = res.pairs_df.filter(
        pl.col("is_significant") & (pl.col("nhood_u") == 0) & (pl.col("nhood_v") == 1)
    )
    assert true_pair.height >= 1


def test_small_imbalanced_cohort_blocks_sparse_continuous_leverage_spikes() -> None:
    """Verifies small/imbalanced cohorts (e.g. 15 R / 2 NR) block sparse 2-donor Student-t leverage spikes."""
    rng = np.random.default_rng(1011)
    n_nhoods, n_samples = 30, 17
    nhoods_mat = _disjoint_nhoods_matrix(n_nhoods, cells_per_nhood=3)
    libs = np.full(n_samples, 400.0)
    y_resp = ["responder"] * 15 + ["non-responder"] * 2

    counts = np.zeros((n_nhoods, n_samples), dtype=np.float64)
    # Sparse 2-donor spikes across the 2 NR donors that would trigger r -> 1.0 without arm/support gating
    for k in range(n_nhoods):
        chosen = rng.choice(15, size=2, replace=False) if k >= 6 else np.array([15, 16])
        counts[k, chosen] = 15.0

    res = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=y_resp,
        sample_library_sizes=Some(tuple(libs.tolist())),
        config=PairwiseConfig(method=PairwiseMethod.HYBRID_ALL),
    ).unwrap()

    assert res.n_significant_continuous_pairs == 0
    assert res.n_significant_quadrant_pairs == 0


# ---------------------------------------------------------------------------
# 3. Outcome-Blind Prefiltering, Super-Neighborhood Pooling & External Priors
# ---------------------------------------------------------------------------


@settings(max_examples=30, deadline=None)
@given(
    n_nhoods=st.integers(min_value=4, max_value=8),
    n_per_arm_stratum=st.integers(min_value=4, max_value=8),
    n_strata=st.integers(min_value=1, max_value=3),
    mode=st.sampled_from(list(PresenceMode)),
    seed=st.integers(min_value=1, max_value=20000),
)
def test_prop_tarone_max_iut_z_universal_upper_bound_and_label_invariance(
    n_nhoods: int,
    n_per_arm_stratum: int,
    n_strata: int,
    mode: PresenceMode,
    seed: int,
) -> None:
    """Proves Z_max^IUT is (a) a universal upper bound on observed IUT Z and (b) 100% invariant to label permutations."""
    rng = np.random.default_rng(seed)
    n_samples = 2 * n_per_arm_stratum * n_strata
    counts = rng.poisson(lam=2.2, size=(n_nhoods, n_samples)).astype(np.float64)
    counts[:, :2] = np.maximum(counts[:, :2], 2.0)
    libs = rng.uniform(80.0, 900.0, size=n_samples)

    labels_orig: list[str] = []
    labels_perm: list[str] = []
    strata_names: list[str] = []
    for k in range(n_strata):
        arm_labels = ["responder"] * n_per_arm_stratum + ["non-responder"] * n_per_arm_stratum
        perm_labels = list(arm_labels)
        rng.shuffle(perm_labels)
        labels_orig.extend(arm_labels)
        labels_perm.extend(perm_labels)
        strata_names.extend([f"S{k}"] * (2 * n_per_arm_stratum))

    pres = estimate_nhood_presence(
        count_matrix=counts,
        library_sizes=Some(libs),
        config=PresenceConfig(mode=mode, min_active_donors=2),
    ).unwrap()
    v_act = int(pres.active_nhood_indices.shape[0])
    assert v_act >= 2

    p = pres.presence_prob
    q = pres.absence_prob
    m = pres.observability

    def _build_tarone_inputs(
        lbls: list[str],
    ) -> tuple[
        list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
        list[np.ndarray],
        list[tuple[int, int]],
    ]:
        quad_totals: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []
        eff_totals: list[np.ndarray] = []
        arm_counts: list[tuple[int, int]] = []
        for k in range(n_strata):
            sl = slice(k * 2 * n_per_arm_stratum, (k + 1) * 2 * n_per_arm_stratum)
            y_k = np.array([lbl == "responder" for lbl in lbls[sl]], dtype=bool)
            # Notice: unlabelled totals depend ONLY on sl, never on y_k, except integer arm counts!
            pk, qk, mk = p[:, sl], q[:, sl], m[:, sl]
            t_pp = pk @ pk.T
            t_pm = pk @ qk.T
            t_mp = qk @ pk.T
            t_mm = qk @ qk.T
            quad_totals.append((t_pp, t_pm, t_mp, t_mm))
            eff_totals.append(mk @ mk.T)
            arm_counts.append((int(np.sum(y_k)), int(np.sum(~y_k))))
        return quad_totals, eff_totals, arm_counts

    qt_orig, et_orig, ac_orig = _build_tarone_inputs(labels_orig)
    qt_perm, et_perm, ac_perm = _build_tarone_inputs(labels_perm)

    quad_to_idx = {qs: i for i, qs in enumerate(ALL_QUADRANT_STATES)}
    z_max_by_q: dict[str, np.ndarray] = {}
    for qi, qs in enumerate(ALL_QUADRANT_STATES):
        fv, fu = adjacent_flips(qs)
        z_orig = compute_tarone_max_iut_z_strata(
            qt_orig, et_orig, ac_orig, qi, quad_to_idx[fv], quad_to_idx[fu]
        )
        z_perm = compute_tarone_max_iut_z_strata(
            qt_perm, et_perm, ac_perm, qi, quad_to_idx[fv], quad_to_idx[fu]
        )
        # Theorem 1b: Exact invariance under within-stratum clinical response permutations
        np.testing.assert_allclose(z_orig, z_perm, atol=1e-12)
        z_max_by_q[qs.value] = z_orig

    # Theorem 1a: Observed IUT Z never exceeds Z_max^IUT on any evaluated candidate
    nhoods_mat = _disjoint_nhoods_matrix(n_nhoods)
    strata_arg: Maybe[tuple[str, ...]] = Some(tuple(strata_names)) if n_strata > 1 else Nothing
    res = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=labels_orig,
        sample_library_sizes=Some(tuple(libs.tolist())),
        strata_labels=strata_arg,
        config=PairwiseConfig(
            method=PairwiseMethod.QUADRANT_2X2,
            presence=PresenceConfig(mode=mode, min_active_donors=2),
            screen=ScreenConfig(min_donors_per_quadrant=1.0, min_effective_samples=2.0),
            refinement=RefinementConfig(enabled=False),
            modules=PairModuleConfig(enabled=False),
        ),
    ).unwrap()

    act_to_local = {int(orig_i): loc_i for loc_i, orig_i in enumerate(pres.active_nhood_indices.tolist())}
    for row in res.pairs_df.to_dicts():
        u_loc = act_to_local[int(row["nhood_u"])]
        v_loc = act_to_local[int(row["nhood_v"])]
        q_str = str(row["quadrant"])
        p_iut = float(row["pvalue_iut"])
        z_obs = 0.0 if p_iut >= 1.0 else float(norm.isf(0.5 * max(p_iut, 1e-300)))
        z_bound = float(z_max_by_q[q_str][u_loc, v_loc])
        assert z_obs <= z_bound + 1e-5


def test_super_nhood_pooling_preserves_lineages_deduplicates_cells_and_boosts_power() -> None:
    """Verifies Leiden super-neighborhood pooling deduplicates shared cells, preserves lineages,
    and rescues split-sphere synergies from Poisson dropout + O(V^2) multiplicity.
    """
    rng = np.random.default_rng(2026)
    n_samples = 72
    cells_per_sample = 30
    n_cells = n_samples * cells_per_sample

    # 18 donors per quadrant of (CD8_T, mregDC): ++ has 16/18 R; +-, -+, -- each have 3/18 R
    u_true = np.array([1] * 18 + [1] * 18 + [0] * 18 + [0] * 18)
    v_true = np.array([1] * 18 + [0] * 18 + [1] * 18 + [0] * 18)
    y_resp = (
        ["responder"] * 16 + ["non-responder"] * 2
        + ["responder"] * 3 + ["non-responder"] * 15
        + ["responder"] * 3 + ["non-responder"] * 15
        + ["responder"] * 3 + ["non-responder"] * 15
    )

    # 4 overlapping CD8_T sub-spheres (0..3), 4 overlapping mregDC sub-spheres (4..7),
    # and 52 disjoint noise neighborhoods (8..59) -> V = 60 neighborhoods total
    n_nhoods = 60
    nhood_labels = ["M04_EffT"] * 4 + ["M41_mregDC"] * 4 + [f"M_Noise_{k % 6}" for k in range(52)]
    rows: list[int] = []
    cols: list[int] = []
    sample_codes = np.repeat(np.arange(n_samples, dtype=np.int64), cells_per_sample)

    # Each positive donor gets a small number of lineage cells dispersed across overlapping sub-spheres
    # so individual sub-spheres suffer stochastic dropout (~45% dropout per individual sphere!),
    # while the pooled super-neighborhood union detects the lineage reliably!
    for j in range(n_samples):
        base = j * cells_per_sample
        if u_true[j] == 1:
            # Activate 2 overlapping sub-spheres out of 4 in donor j, sharing anchor cells across donors
            anchor = int(j % 3)
            for cid_off, nh_pair in ((0, (anchor, anchor + 1)), (1, (anchor, anchor + 1))):
                for nh_i in nh_pair:
                    rows.append(base + cid_off)
                    cols.append(nh_i)
        if v_true[j] == 1:
            anchor_v = 4 + int((j + 1) % 3)
            for cid_off, nh_pair in ((4, (anchor_v, anchor_v + 1)), (5, (anchor_v, anchor_v + 1))):
                for nh_i in nh_pair:
                    rows.append(base + cid_off)
                    cols.append(nh_i)
        # Noise neighborhoods 8..59 occupy offsets 8..29
        for k_noise in range(52):
            if rng.random() < 0.35:
                cid_noise = base + 8 + (k_noise % 22)
                rows.append(cid_noise)
                cols.append(8 + k_noise)

    nhoods_mat = sp.csc_matrix(
        (np.ones(len(rows), dtype=np.float64), (rows, cols)),
        shape=(n_cells, n_nhoods),
    )
    sample_ind = sp.csr_matrix(
        (np.ones(n_cells, dtype=np.float64), (np.arange(n_cells, dtype=np.int64), sample_codes)),
        shape=(n_cells, n_samples),
    )
    counts = np.asarray(nhoods_mat.T.dot(sample_ind).toarray(), dtype=np.float64)
    libs = np.full(n_samples, float(cells_per_sample), dtype=np.float64)

    # 1. Verify exact deduplication & strict lineage preservation in coarse_grain_nhoods_by_graph
    cg_exact = coarse_grain_nhoods_by_graph(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        cell_sample_codes=Some(sample_codes),
        sample_library_sizes=Some(libs),
        nhood_group_labels=Some(nhood_labels),
        min_jaccard=0.10,
    ).unwrap()
    cg_mult = coarse_grain_nhoods_by_graph(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        cell_sample_codes=Nothing,
        sample_library_sizes=Some(libs),
        nhood_group_labels=Some(nhood_labels),
        min_jaccard=0.10,
    ).unwrap()

    assert len(cg_exact.super_nhood_labels) < n_nhoods
    for members in cg_exact.member_nhood_indices:
        member_lbls = {nhood_labels[m] for m in members}
        assert len(member_lbls) == 1
    assert np.all(np.sum(cg_exact.counts_super, axis=0) <= libs + 1e-9)
    assert np.all(cg_mult.counts_super <= np.expand_dims(libs, axis=0) + 1e-9)

    # 2. Compare Unfiltered (PrefilterMode.NONE) vs Composite Super-Nhood + Tarone
    res_unfiltered = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=y_resp,
        sample_library_sizes=Some(tuple(libs.tolist())),
        nhood_group_labels=Some(nhood_labels),
        cell_sample_codes=Some(sample_codes),
        config=PairwiseConfig(
            method=PairwiseMethod.QUADRANT_2X2,
            overlap=OverlapConfig(scope_filter=PairScopeFilter.CROSS_CELL_TYPE_ONLY),
            presence=PresenceConfig(mode=PresenceMode.HARD_DEPTH_GATED, min_cell_count=1.0),
            prefilter=PrefilterConfig(mode=PrefilterMode.NONE),
            refinement=RefinementConfig(enabled=False),
        ),
    ).unwrap()

    res_super = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=y_resp,
        sample_library_sizes=Some(tuple(libs.tolist())),
        nhood_group_labels=Some(nhood_labels),
        cell_sample_codes=Some(sample_codes),
        config=PairwiseConfig(
            method=PairwiseMethod.QUADRANT_2X2,
            overlap=OverlapConfig(scope_filter=PairScopeFilter.CROSS_CELL_TYPE_ONLY),
            presence=PresenceConfig(mode=PresenceMode.HARD_DEPTH_GATED, min_cell_count=1.0),
            prefilter=PrefilterConfig(
                mode=PrefilterMode.COMPOSITE_SUPER_TARONE,
                super_nhood_min_jaccard=0.10,
                min_variance_quantile=0.0,
            ),
            refinement=RefinementConfig(enabled=False),
        ),
    ).unwrap()

    best_p_unfiltered = float(res_unfiltered.pairs_df["pvalue_iut"].min())
    best_p_super = float(res_super.pairs_df["pvalue_iut"].min())
    # Super-neighborhood pooling both reduces evaluated pairs AND sharpens raw pvalue_iut (>10x lower p)!
    assert res_super.n_evaluated_pairs < res_unfiltered.n_evaluated_pairs
    assert best_p_super < best_p_unfiltered * 0.10

    sig_super = res_super.pairs_df.filter(
        pl.col("is_significant")
        & (pl.col("quadrant") == QuadrantState.BOTH_PRESENT.value)
        & (pl.col("cell_type_u") == "M04_EffT")
        & (pl.col("cell_type_v") == "M41_mregDC")
    )
    assert sig_super.height == 1
    assert float(sig_super["fdr_iut"][0]) < 0.05


def test_outcome_blind_prefilters_and_external_priors_control_h10_piggybacking() -> None:
    """Verifies zero false positives under adversarial strong-singleton H_10 piggybacking across
    all outcome-blind prefiltering modes and external prior weightings.
    """
    rng = np.random.default_rng(31415)
    n_nhoods, n_samples = 36, 76
    nhoods_mat = _disjoint_nhoods_matrix(n_nhoods, cells_per_nhood=4)
    libs = np.full(n_samples, 800.0)
    y_resp = ["responder"] * 38 + ["non-responder"] * 38

    counts = np.zeros((n_nhoods, n_samples), dtype=np.float64)
    # Nhood 0 is a massive univariate R singleton (35/38 R vs 2/38 NR); 1..35 are pure noise
    counts[0, :35] = 25.0
    counts[0, 38:40] = 25.0
    for k in range(1, n_nhoods):
        counts[k, :] = (rng.random(n_samples) < 0.42).astype(np.float64) * 18.0

    labels = ["M13_Tex"] + [f"M_Noise_{k % 5}" for k in range(1, n_nhoods)]
    w_prior = build_external_prior_weight_matrix(
        nhood_state_labels=labels,
        state_pair_weights={("M13_Tex", "M_Noise_1"): 6.0, ("M13_Tex", "M_Noise_2"): 4.0},
        base_floor=0.25,
        max_cap=8.0,
    ).unwrap()

    for mode in (
        PrefilterMode.TARONE_ONLY,
        PrefilterMode.MARGINAL_ENTROPY_VARIANCE,
        PrefilterMode.COMPOSITE_SUPER_TARONE,
    ):
        res = analyze_nhood_pairs(
            nhoods_mat=nhoods_mat,
            count_matrix=counts,
            response_labels=y_resp,
            sample_library_sizes=Some(tuple(libs.tolist())),
            nhood_group_labels=Some(labels),
            prior_pair_weights=Some(w_prior),
            config=PairwiseConfig(
                method=PairwiseMethod.HYBRID_ALL,
                prefilter=PrefilterConfig(mode=mode, apply_tarone_iut_pruning=True),
                screen=ScreenConfig(fdr_threshold=0.05),
            ),
        ).unwrap()
        assert res.n_significant_quadrant_pairs == 0
        assert res.n_significant_continuous_pairs == 0


def test_external_prior_weights_unlock_targeted_true_pair() -> None:
    """Verifies independent cross-modality prior weights lower weighted FDR on true prioritized pairs
    buried under heavy background multiplicity without inflating background noise.
    """
    rng = np.random.default_rng(27182)
    n_nhoods, n_samples = 80, 68
    nhoods_mat = _disjoint_nhoods_matrix(n_nhoods, cells_per_nhood=3)
    libs = np.full(n_samples, 600.0)

    # Moderate ++ synergy on (0, 1): 13/17 R in ++ vs 4/17 R in +-, 4/17 R in -+, 4/17 R in --
    u_pres = np.array([1] * 17 + [1] * 17 + [0] * 17 + [0] * 17, dtype=np.float64)
    v_pres = np.array([1] * 17 + [0] * 17 + [1] * 17 + [0] * 17, dtype=np.float64)
    y_resp = (
        ["responder"] * 13 + ["non-responder"] * 4
        + ["responder"] * 4 + ["non-responder"] * 13
        + ["responder"] * 4 + ["non-responder"] * 13
        + ["responder"] * 4 + ["non-responder"] * 13
    )

    counts = np.zeros((n_nhoods, n_samples), dtype=np.float64)
    counts[0, :] = u_pres * 16.0
    counts[1, :] = v_pres * 16.0
    for k in range(2, n_nhoods):
        counts[k, :] = (rng.random(n_samples) < 0.45).astype(np.float64) * 14.0

    state_labels = ["M29", "M37"] + [f"M{(k % 12) + 1:02d}" for k in range(2, n_nhoods)]
    w_prior = build_external_prior_weight_matrix(
        nhood_state_labels=state_labels,
        state_pair_weights={("M29", "M37"): 8.0},
        base_floor=0.15,
        max_cap=8.0,
    ).unwrap()

    cfg = PairwiseConfig(
        method=PairwiseMethod.QUADRANT_2X2,
        screen=ScreenConfig(fdr_threshold=0.10),
        refinement=RefinementConfig(enabled=False),
        modules=PairModuleConfig(enabled=False),
    )
    res_unweighted = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=y_resp,
        sample_library_sizes=Some(tuple(libs.tolist())),
        nhood_group_labels=Some(state_labels),
        config=cfg,
    ).unwrap()
    res_prior = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=y_resp,
        sample_library_sizes=Some(tuple(libs.tolist())),
        nhood_group_labels=Some(state_labels),
        prior_pair_weights=Some(w_prior),
        config=cfg,
    ).unwrap()

    row_unw = res_unweighted.pairs_df.filter(
        (pl.col("nhood_u") == 0)
        & (pl.col("nhood_v") == 1)
        & (pl.col("quadrant") == QuadrantState.BOTH_PRESENT.value)
    )
    row_pri = res_prior.pairs_df.filter(
        (pl.col("nhood_u") == 0)
        & (pl.col("nhood_v") == 1)
        & (pl.col("quadrant") == QuadrantState.BOTH_PRESENT.value)
    )
    assert row_unw.height == 1 and row_pri.height == 1
    assert float(row_pri["fdr_quadrant"][0]) < float(row_unw["fdr_quadrant"][0])
    assert float(row_pri["fdr_iut"][0]) < float(row_unw["fdr_iut"][0])
    assert float(row_pri["fdr_iut"][0]) < 0.10
    noise_sig = res_prior.pairs_df.filter(
        pl.col("is_significant") & ~((pl.col("nhood_u") == 0) & (pl.col("nhood_v") == 1))
    )
    assert noise_sig.height == 0


