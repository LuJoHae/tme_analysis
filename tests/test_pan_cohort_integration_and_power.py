"""Property-based and statistical calibration tests for pan-cohort Harmony integration & pairwise power."""

from __future__ import annotations

from hypothesis import given, settings, strategies as st
import numpy as np
import polars as pl
from returns.maybe import Nothing, Some
from returns.result import Success
from scipy import sparse as sp
from scipy.stats import norm

from milopy.pairwise import (
    OverlapConfig,
    PairScope,
    PairScopeFilter,
    PairwiseConfig,
    PairwiseMethod,
    PresenceConfig,
    PresenceMode,
    RefinementConfig,
    ScreenConfig,
    TmeCompartment,
    analyze_nhood_pairs,
    build_compartment_observability_mask,
    combine_pair_regimes_acat,
    depth_residualize_continuous_features,
    estimate_nhood_presence,
    meta_analyze_cohort_pairs_iut,
)
from scripts.milopy_analysis.integrate_pan_cohort_and_boost_power import (
    HarmonyRegime,
    PreparedSuperUnit,
    TmeSuperGroup,
    _map_dominant_regime_to_visual_polarity,
    _oriented_pair_z_acat,
    _resolve_oriented_super_nhoods,
    build_compartment_pair_pairwise_umap_bundle,
    build_compartment_pair_umap_cell_bundle,
    build_pan_cohort_compartment_pair_umap_grid,
    build_pan_cohort_core_hubs_z_acat_heatmap_figure,
    build_pan_cohort_z_acat_heatmap_figure,
    build_standalone_compartment_pair_umap_composite,
    classify_compartment_supergroup,
    cluster_states_within_compartments,
    enumerate_active_broad_compartment_pairs,
    evaluate_embedding_integration_metrics,
    format_stratified_perm_stars,
    run_harmony_embedding_regime,
    select_core_hub_state_codes,
)


def _make_disjoint_incidence(n_nhoods: int, cells_per_nhood: int = 8) -> sp.csc_matrix:
    """Constructs strictly disjoint single-cell neighborhood incidence columns."""
    n_cells = n_nhoods * cells_per_nhood
    rows = np.arange(n_cells, dtype=np.int64)
    cols = np.repeat(np.arange(n_nhoods, dtype=np.int64), cells_per_nhood)
    data = np.ones(n_cells, dtype=np.float64)
    return sp.csc_matrix((data, (rows, cols)), shape=(n_cells, n_nhoods))


@settings(max_examples=15, deadline=5000)
@given(
    corr=st.floats(min_value=-0.80, max_value=0.80, allow_nan=False, allow_infinity=False),
    seed=st.integers(min_value=1, max_value=5000),
)
def test_prop_acat_omnibus_composite_null_control_and_power(corr: float, seed: int) -> None:
    """Verifies ACAT controls composite-null H_00 / H_10 Type-I error and boosts dual-regime power."""
    rng = np.random.default_rng(seed)
    n_null_pairs = 500

    # Generate H_10 singleton-piggybacking pairs (u -> +inf, v ~ correlated Gaussian noise)
    zv1 = rng.normal(size=n_null_pairs)
    zv0 = rng.normal(size=n_null_pairs)
    z_lr = corr * zv1 + np.sqrt(max(1.0 - corr * corr, 1e-6)) * rng.normal(size=n_null_pairs)
    z_co = corr * zv0 + np.sqrt(max(1.0 - corr * corr, 1e-6)) * rng.normal(size=n_null_pairs)

    p_pp = np.clip(2.0 * norm.sf(zv1), 0.0, 1.0)
    p_pm = np.clip(2.0 * norm.sf(-zv1), 0.0, 1.0)
    p_mp = np.clip(2.0 * norm.sf(zv0), 0.0, 1.0)
    p_mm = np.clip(2.0 * norm.sf(-zv0), 0.0, 1.0)
    p_lr = np.clip(2.0 * norm.sf(z_lr), 0.0, 1.0)
    p_co = np.clip(2.0 * norm.sf(np.abs(z_co)), 0.0, 1.0)

    q_rows: list[dict[str, object]] = []
    c_rows: list[dict[str, object]] = []
    for idx in range(n_null_pairs):
        u_id = 2 * idx
        v_id = 2 * idx + 1
        for q_str, p_val in (("++", p_pp[idx]), ("+-", p_pm[idx]), ("-+", p_mp[idx]), ("--", p_mm[idx])):
            q_rows.append(
                {
                    "nhood_u": u_id,
                    "nhood_v": v_id,
                    "quadrant": q_str,
                    "direction": "responder",
                    "pair_scope": PairScope.CROSS_COMPARTMENT.value,
                    "cell_type_u": "M04",
                    "cell_type_v": "M17",
                    "jaccard_overlap": 0.0,
                    "spatial_weight": 1.0,
                    "eff_donors_responder": 12.0,
                    "eff_donors_non_responder": 12.0,
                    "eff_donors_total": 48.0,
                    "log_or_quadrant": 0.2,
                    "pvalue_iut": float(p_val),
                }
            )
        for ax_str, p_val in (("log_ratio_balance", p_lr[idx]), ("co_abundance_product", p_co[idx])):
            c_rows.append(
                {
                    "nhood_u": u_id,
                    "nhood_v": v_id,
                    "axis": ax_str,
                    "direction": "responder",
                    "pair_scope": PairScope.CROSS_COMPARTMENT.value,
                    "cell_type_u": "M04",
                    "cell_type_v": "M17",
                    "jaccard_overlap": 0.0,
                    "spatial_weight": 1.0,
                    "effect_size": 0.1,
                    "pvalue_iut": float(p_val),
                }
            )

    null_res = combine_pair_regimes_acat(pl.DataFrame(q_rows), pl.DataFrame(c_rows), fdr_threshold=0.10)
    assert isinstance(null_res, Success)
    null_df = null_res.unwrap()
    assert null_df.height == n_null_pairs
    emp_alpha_05 = float(np.mean(null_df["pvalue_acat_iut"].to_numpy() <= 0.05))
    assert emp_alpha_05 <= 0.085

    # Inject 10 dual discrete + continuous cross-compartment signal pairs
    for s_idx in range(10):
        u_id = 10_000 + 2 * s_idx
        v_id = 10_000 + 2 * s_idx + 1
        for q_str, p_val in (("++", 1.0), ("+-", 2.5e-4), ("-+", 1.0), ("--", 1.0)):
            q_rows.append(
                {
                    "nhood_u": u_id,
                    "nhood_v": v_id,
                    "quadrant": q_str,
                    "direction": "responder",
                    "pair_scope": PairScope.CROSS_COMPARTMENT.value,
                    "cell_type_u": "M04",
                    "cell_type_v": "M17",
                    "jaccard_overlap": 0.0,
                    "spatial_weight": 1.0,
                    "eff_donors_responder": 14.0,
                    "eff_donors_non_responder": 10.0,
                    "eff_donors_total": 48.0,
                    "log_or_quadrant": 2.4,
                    "pvalue_iut": p_val,
                }
            )
        c_rows.append(
            {
                "nhood_u": u_id,
                "nhood_v": v_id,
                "axis": "log_ratio_balance",
                "direction": "responder",
                "pair_scope": PairScope.CROSS_COMPARTMENT.value,
                "cell_type_u": "M04",
                "cell_type_v": "M17",
                "jaccard_overlap": 0.0,
                "spatial_weight": 1.0,
                "effect_size": 0.55,
                "pvalue_iut": 2.0e-4,
            }
        )

    alt_res = combine_pair_regimes_acat(pl.DataFrame(q_rows), pl.DataFrame(c_rows), fdr_threshold=0.10)
    assert isinstance(alt_res, Success)
    alt_df = alt_res.unwrap()
    sig_df = alt_df.filter(pl.col("is_significant"))
    tp_hits = int(sig_df.filter(pl.col("nhood_u") >= 10_000).height)
    assert tp_hits == 10


def test_depth_residualize_continuous_features_removes_library_confounding() -> None:
    """Verifies outcome-blind depth residualization zeroes log10(library_sizes) correlation."""
    rng = np.random.default_rng(123)
    libs = 10.0 ** rng.uniform(2.0, 4.2, size=80)
    log_libs = np.log10(libs)
    noise = rng.normal(scale=0.35, size=(6, 80))
    features = 4.0 + np.outer(np.array([1.5, -1.2, 0.8, 0.0, 2.1, -0.6]), log_libs) + noise

    res = depth_residualize_continuous_features(features, libs, preserve_row_means=True)
    assert isinstance(res, Success)
    cleaned = res.unwrap()
    assert cleaned.shape == features.shape
    np.testing.assert_allclose(np.mean(cleaned, axis=1), np.mean(features, axis=1), atol=1e-10)
    centered_l = log_libs - np.mean(log_libs)
    dots = np.abs(cleaned @ centered_l)
    np.testing.assert_allclose(dots, np.zeros(6), atol=1e-10)


def test_compartment_observability_mask_protects_cd45_sorted_cohorts() -> None:
    """Verifies CD45+ FACS-sorted cohort masking preserves Immune<->Immune J and protects Stromal gates."""
    rng = np.random.default_rng(77)
    j_wt = 48
    j_cd45 = 64
    j_tot = j_wt + j_cd45

    y_wt = ["responder"] * (j_wt // 2) + ["non-responder"] * (j_wt // 2)
    y_cd45 = ["responder"] * (j_cd45 // 2) + ["non-responder"] * (j_cd45 // 2)
    y_all = tuple(y_wt + y_cd45)
    cohorts = tuple(["WHOLE_TUMOR"] * j_wt + ["CD45_SORTED"] * j_cd45)
    libs = np.full(j_tot, 1000.0, dtype=np.float64)

    # Row 0: M04 (T_NK), Row 1: M30 (B_LINEAGE), Row 2: M17 (STROMAL_MESENCHYMAL), Rows 3..5: background
    counts = np.zeros((6, j_tot), dtype=np.float64)
    nhood_labels = ("M04", "M30", "M17", "M09", "M10", "M14")

    # Whole-tumor cohort: M04+ & M17- strongly enriched in Responders; M04+ & M30+ synergy in Responders
    for j in range(j_wt):
        is_r = j < (j_wt // 2)
        if is_r:
            has_m04 = j < 20
            has_m30 = j < 19
            has_m17 = j >= 21
        else:
            has_m04 = j < (j_wt // 2 + 10)
            has_m30 = j >= (j_wt // 2 + 12)
            has_m17 = j < (j_wt // 2 + 19)
        counts[0, j] = 25.0 if has_m04 else 0.0
        counts[1, j] = 20.0 if has_m30 else 0.0
        counts[2, j] = 30.0 if has_m17 else 0.0

    # CD45+ cohort: M17 is 0 across ALL 64 donors (FACS stripped); M04+ & M30+ synergy holds!
    for j in range(j_wt, j_tot):
        local_j = j - j_wt
        is_r = local_j < (j_cd45 // 2)
        if is_r:
            has_m04 = local_j < 27
            has_m30 = local_j < 26
        else:
            has_m04 = local_j < 12
            has_m30 = local_j >= 16
        counts[0, j] = 25.0 if has_m04 else 0.0
        counts[1, j] = 20.0 if has_m30 else 0.0
        counts[2, j] = 0.0  # FACS sorted CD45+: zero stromal cells!

    for bg_row in (3, 4, 5):
        counts[bg_row, :] = rng.binomial(1, 0.45, size=j_tot) * 15.0

    pres_raw = estimate_nhood_presence(counts, Some(libs)).unwrap()
    masked_res = build_compartment_observability_mask(
        presence=pres_raw,
        nhood_labels=nhood_labels,
        cohort_labels=cohorts,
        min_cohort_compartment_share=0.005,
    )
    assert isinstance(masked_res, Success)
    pres_masked = masked_res.unwrap()

    # M04 (row 0) and M30 (row 1) remain observable across all 112 donors; M17 (row 2) is 0 in CD45_SORTED
    assert np.all(pres_masked.observability[0, :] > 0.9)
    assert np.all(pres_masked.observability[1, :] > 0.9)
    assert np.all(pres_masked.observability[2, :j_wt] > 0.9)
    np.testing.assert_allclose(pres_masked.observability[2, j_wt:], 0.0, atol=1e-12)

    nhoods_mat = _make_disjoint_incidence(6)
    cfg_masked = PairwiseConfig(
        method=PairwiseMethod.HYBRID_ALL,
        overlap=OverlapConfig(scope_filter=PairScopeFilter.CROSS_COMPARTMENT_ONLY),
        presence=PresenceConfig(
            mode=PresenceMode.HARD_DEPTH_GATED,
            mask_unobserved_compartments=True,
            min_cohort_compartment_share=0.005,
        ),
        screen=ScreenConfig(fdr_threshold=0.10, dmpp_null_draws=40),
        refinement=RefinementConfig(enabled=False),
    )
    res_pooled = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=y_all,
        sample_library_sizes=Some(tuple(libs.tolist())),
        strata_labels=Some(cohorts),
        cohort_labels=Some(cohorts),
        nhood_group_labels=Some(nhood_labels),
        config=cfg_masked,
    )
    assert isinstance(res_pooled, Success)
    pw = res_pooled.unwrap()

    # Check M04 <-> M30 uses full J=112 while M04 <-> M17 uses strictly J=48 whole-tumor donors
    m04_m30 = pw.pairs_df.filter(
        (pl.col("nhood_u") == 0) & (pl.col("nhood_v") == 1) & (pl.col("quadrant") == "++")
    )
    m04_m17 = pw.pairs_df.filter(
        (pl.col("nhood_u") == 0) & (pl.col("nhood_v") == 2) & (pl.col("quadrant") == "+-")
    )
    assert m04_m30.height == 1
    assert m04_m17.height == 1
    assert float(m04_m30["eff_donors_total"].item()) > 105.0
    assert float(m04_m17["eff_donors_total"].item()) <= float(j_wt) + 0.01
    assert bool(m04_m30["is_significant"].item())
    assert bool(m04_m17["is_significant"].item())


def test_two_stage_meta_analysis_recovers_subthreshold_pan_cohort_pairs() -> None:
    """Verifies 2-stage DerSimonian-Laird Meta-IUT unlocks 5 sub-threshold cohorts & rejects outliers."""
    rng = np.random.default_rng(2026)
    nhoods_mat = _make_disjoint_incidence(6)
    labels = ("M04", "M17", "M09", "M11", "M30", "M16")
    cfg = PairwiseConfig(
        method=PairwiseMethod.QUADRANT_2X2,
        overlap=OverlapConfig(scope_filter=PairScopeFilter.CROSS_COMPARTMENT_ONLY),
        screen=ScreenConfig(fdr_threshold=0.10, min_donors_per_quadrant=2.0),
        refinement=RefinementConfig(enabled=False),
    )

    cohort_tables = {}
    cohort_sizes = {}
    for c_idx in range(5):
        cname = f"COHORT_{c_idx + 1}"
        j_c = 18
        y_c = ("responder",) * 9 + ("non-responder",) * 9
        libs = tuple([1000.0] * j_c)
        counts = np.zeros((6, j_c), dtype=np.float64)

        # Consistent moderate M04+ / M17- gate across all 5 small cohorts (subthreshold individually)
        # Responders: 5 (+-), 1 (++), 2 (-+), 1 (--)
        # Non-Responders: 1 (+-), 3 (++), 2 (-+), 3 (--)
        r_states = [("+-", 5), ("++", 1), ("-+", 2), ("--", 1)]
        nr_states = [("+-", 1), ("++", 3), ("-+", 2), ("--", 3)]
        col = 0
        for state_list in (r_states, nr_states):
            for q_code, reps in state_list:
                for _ in range(reps):
                    counts[0, col] = 20.0 if q_code[0] == "+" else 0.0
                    counts[1, col] = 20.0 if q_code[1] == "+" else 0.0
                    col += 1

        # Single-cohort outlier on M09 (row 2) <-> M11 (row 3) ONLY in COHORT_1; pure opposite/null in 2..5
        if c_idx == 0:
            counts[2, :8] = 25.0
            counts[3, :8] = 25.0
            counts[2, 9:13] = 25.0
            counts[3, 13:18] = 25.0
        else:
            counts[2, :] = rng.binomial(1, 0.5, size=j_c) * 20.0
            counts[3, :] = rng.binomial(1, 0.5, size=j_c) * 20.0

        counts[4, :] = rng.binomial(1, 0.5, size=j_c) * 20.0
        counts[5, :] = rng.binomial(1, 0.5, size=j_c) * 20.0

        pw_c = analyze_nhood_pairs(
            nhoods_mat=nhoods_mat,
            count_matrix=counts,
            response_labels=y_c,
            sample_library_sizes=Some(libs),
            nhood_group_labels=Some(labels),
            config=cfg,
        ).unwrap()
        # Each individual J=16 cohort has 0 significant hits at FDR < 0.10 for M04 <-> M17
        sub_m04_m17 = pw_c.pairs_df.filter(
            (pl.col("cell_type_u") == "M04")
            & (pl.col("cell_type_v") == "M17")
            & (pl.col("quadrant") == "+-")
        )
        assert not bool(sub_m04_m17["is_significant"].item())
        cohort_tables[cname] = pw_c
        cohort_sizes[cname] = j_c

    meta_res = meta_analyze_cohort_pairs_iut(
        cohort_pair_tables=cohort_tables,
        cohort_sample_sizes=cohort_sizes,
        fdr_threshold=0.10,
        require_cross_compartment=True,
    )
    assert isinstance(meta_res, Success)
    meta_df = meta_res.unwrap()

    target = meta_df.filter(
        (pl.col("state_a") == "M04")
        & (pl.col("state_b") == "M17")
        & (pl.col("canonical_regime") == "quadrant:+-")
    )
    assert target.height == 1
    assert bool(target["is_significant"].item())
    assert float(target["fdr_meta_iut"].item()) < 0.02
    assert float(target["min_loco_z_iut"].item()) > 2.0

    outlier = meta_df.filter(
        (pl.col("state_a") == "M09")
        & (pl.col("state_b") == "M11")
        & (pl.col("canonical_regime") == "quadrant:++")
    )
    if outlier.height > 0:
        assert not bool(outlier["is_significant"].item())



def test_harmonypy_multi_cohort_integration_improves_ilisi_and_preserves_clisi() -> None:
    """Verifies harmonypy increases cohort_iLISI > 2x while keeping compartment_cLISI < 1.10."""
    rng = np.random.default_rng(4242)
    n_cohorts = 3
    n_compartments = 4
    cells_per_block = 75
    n_dims = 16

    comp_names = ("t_nk", "b_lineage", "myeloid_mast", "stromal_mesenchymal")
    meta_names = ("M04", "M30", "M10", "M17")

    # Well-separated biological compartment centroids in dims 0..7
    comp_centroids = np.zeros((n_compartments, n_dims), dtype=np.float64)
    for k in range(n_compartments):
        comp_centroids[k, 2 * k : 2 * k + 2] = 6.5

    # Orthogonal technical cohort batch shifts in dims 8..13
    cohort_shifts = np.zeros((n_cohorts, n_dims), dtype=np.float64)
    for c in range(n_cohorts):
        cohort_shifts[c, 8 + 2 * c : 10 + 2 * c] = 3.2

    blocks: list[np.ndarray] = []
    cohort_ids: list[str] = []
    substudy_ids: list[str] = []
    sample_uids: list[str] = []
    compartments: list[str] = []
    meta_states: list[str] = []

    for c_idx in range(n_cohorts):
        c_tag = f"COHORT_{c_idx + 1}"
        for k_idx in range(n_compartments):
            pts = (
                comp_centroids[k_idx]
                + cohort_shifts[c_idx]
                + rng.normal(scale=0.55, size=(cells_per_block, n_dims))
            )
            blocks.append(pts)
            cohort_ids.extend([c_tag] * cells_per_block)
            substudy_ids.extend([f"{c_tag}::sub{(i % 2) + 1}" for i in range(cells_per_block)])
            sample_uids.extend([f"{c_tag}::S{(i % 6):02d}" for i in range(cells_per_block)])
            compartments.extend([comp_names[k_idx]] * cells_per_block)
            meta_states.extend([meta_names[k_idx]] * cells_per_block)

    x_pca = np.vstack(blocks)

    m_i0_res = evaluate_embedding_integration_metrics(
        embedding=x_pca,
        cohort_ids=cohort_ids,
        substudy_ids=substudy_ids,
        sample_uids=sample_uids,
        compartments=compartments,
        meta_states=meta_states,
        regime=HarmonyRegime.I0_UNINTEGRATED_PCA,
        knn_k=25,
        seed=42,
    )
    assert isinstance(m_i0_res, Success)
    m_i0 = m_i0_res.unwrap()

    emb_i2_res = run_harmony_embedding_regime(
        x_pca=x_pca,
        cohort_ids=cohort_ids,
        substudy_ids=substudy_ids,
        sample_uids=sample_uids,
        regime=HarmonyRegime.I2_HARMONY_COHORT_SUBSTUDY,
        nclust=8,
        sigma=0.25,
        max_iter_harmony=15,
        random_state=42,
    )
    assert isinstance(emb_i2_res, Success)
    emb_i2 = emb_i2_res.unwrap()
    assert emb_i2.shape == x_pca.shape

    m_i2_res = evaluate_embedding_integration_metrics(
        embedding=emb_i2,
        cohort_ids=cohort_ids,
        substudy_ids=substudy_ids,
        sample_uids=sample_uids,
        compartments=compartments,
        meta_states=meta_states,
        regime=HarmonyRegime.I2_HARMONY_COHORT_SUBSTUDY,
        knn_k=25,
        seed=42,
    )
    assert isinstance(m_i2_res, Success)
    m_i2 = m_i2_res.unwrap()

    assert m_i2.cohort_ilisi > 2.0 * m_i0.cohort_ilisi
    assert m_i2.cohort_ilisi > 2.2
    assert m_i2.compartment_clisi < 1.10


def test_z_acat_heatmap_clustering_supergroups_and_perm_stars() -> None:
    """Verifies Immune vs. Stromal/Endo super-grouping, within-compartment clustering, and Z_ACAT antisymmetry."""
    assert classify_compartment_supergroup(TmeCompartment.T_NK) == Some(TmeSuperGroup.IMMUNE_CELLS)
    assert classify_compartment_supergroup(TmeCompartment.B_LINEAGE) == Some(
        TmeSuperGroup.IMMUNE_CELLS
    )
    assert classify_compartment_supergroup(TmeCompartment.PLASMA) == Some(
        TmeSuperGroup.IMMUNE_CELLS
    )
    assert classify_compartment_supergroup(TmeCompartment.MYELOID_MAST) == Some(
        TmeSuperGroup.IMMUNE_CELLS
    )
    assert classify_compartment_supergroup(TmeCompartment.STROMAL_MESENCHYMAL) == Some(
        TmeSuperGroup.STROMAL_ENDOTHELIAL
    )
    assert classify_compartment_supergroup(TmeCompartment.ENDOTHELIAL_VASCULAR) == Some(
        TmeSuperGroup.STROMAL_ENDOTHELIAL
    )
    assert classify_compartment_supergroup(TmeCompartment.MALIGNANT_EPITHELIAL) == Some(
        TmeSuperGroup.MALIGNANT_EPITHELIAL
    )
    assert classify_compartment_supergroup(TmeCompartment.OTHER) == Nothing

    assert format_stratified_perm_stars(0.000999) == "***"
    assert format_stratified_perm_stars(0.004) == "**"
    assert format_stratified_perm_stars(0.025) == "*"
    assert format_stratified_perm_stars(0.12) == ""

    assert (
        _oriented_pair_z_acat("M17", "M33", "M17", -3.36, "log_ratio_balance") == -3.36
    )
    assert (
        _oriented_pair_z_acat("M33", "M17", "M17", -3.36, "log_ratio_balance") == 3.36
    )
    assert (
        _oriented_pair_z_acat("M30", "M04", "M04", 3.49, "co_abundance_product") == 3.49
    )

    synth_rows = [
        {
            "state_a": "M04",
            "state_b": "M17",
            "state_a_readable": "M04: CXCL13+LAG3+ Tex",
            "state_b_readable": "M17: ISLR+COL5A2+ myCAF",
            "compartment_a": "t_nk",
            "compartment_b": "stromal_mesenchymal",
            "dominant_regime": "log_ratio_balance",
            "dominant_direction": "responder",
            "z_acat_signed": 4.05,
            "pvalue_acat_iut": 5.0e-5,
            "fdr_acat_iut": 0.004,
            "emp_perm_pvalue": 0.000999,
            "min_loco_z_acat": 3.42,
            "is_significant": True,
        },
        {
            "state_a": "M08",
            "state_b": "M17",
            "state_a_readable": "M08: CCR4+SOCS1+ Th2/Treg",
            "state_b_readable": "M17: ISLR+COL5A2+ myCAF",
            "compartment_a": "t_nk",
            "compartment_b": "stromal_mesenchymal",
            "dominant_regime": "log_ratio_balance",
            "dominant_direction": "responder",
            "z_acat_signed": 4.14,
            "pvalue_acat_iut": 3.5e-5,
            "fdr_acat_iut": 0.004,
            "emp_perm_pvalue": 0.000999,
            "min_loco_z_acat": 3.41,
            "is_significant": True,
        },
        {
            "state_a": "M01",
            "state_b": "M17",
            "state_a_readable": "M01: TRBV18+GCNT4+ Mem T",
            "state_b_readable": "M17: ISLR+COL5A2+ myCAF",
            "compartment_a": "t_nk",
            "compartment_b": "stromal_mesenchymal",
            "dominant_regime": "log_ratio_balance",
            "dominant_direction": "non-responder",
            "z_acat_signed": -0.35,
            "pvalue_acat_iut": 0.72,
            "fdr_acat_iut": 0.91,
            "emp_perm_pvalue": 0.68,
            "min_loco_z_acat": -0.10,
            "is_significant": False,
        },
        {
            "state_a": "M04",
            "state_b": "M30",
            "state_a_readable": "M04: CXCL13+LAG3+ Tex",
            "state_b_readable": "M30: MS4A1+CD79A+ B Cell",
            "compartment_a": "t_nk",
            "compartment_b": "b_lineage",
            "dominant_regime": "co_abundance_product",
            "dominant_direction": "responder",
            "z_acat_signed": 3.49,
            "pvalue_acat_iut": 4.8e-4,
            "fdr_acat_iut": 0.015,
            "emp_perm_pvalue": 0.002997,
            "min_loco_z_acat": 2.91,
            "is_significant": True,
        },
        {
            "state_a": "M10",
            "state_b": "M17",
            "state_a_readable": "M10: FCN1+CD1A+ Mono/cDC",
            "state_b_readable": "M17: ISLR+COL5A2+ myCAF",
            "compartment_a": "myeloid_mast",
            "compartment_b": "stromal_mesenchymal",
            "dominant_regime": "log_ratio_balance",
            "dominant_direction": "responder",
            "z_acat_signed": 3.85,
            "pvalue_acat_iut": 1.2e-4,
            "fdr_acat_iut": 0.007,
            "emp_perm_pvalue": 0.000999,
            "min_loco_z_acat": 2.86,
            "is_significant": True,
        },
        {
            "state_a": "M14",
            "state_b": "M16",
            "state_a_readable": "M14: CXCL9+STAT1+ M1-TAM",
            "state_b_readable": "M16: VEGFA+KDR+ Endo.",
            "compartment_a": "myeloid_mast",
            "compartment_b": "endothelial_vascular",
            "dominant_regime": "log_ratio_balance",
            "dominant_direction": "responder",
            "z_acat_signed": 2.15,
            "pvalue_acat_iut": 0.031,
            "fdr_acat_iut": 0.14,
            "emp_perm_pvalue": 0.038,
            "min_loco_z_acat": 1.62,
            "is_significant": False,
        },
        {
            "state_a": "M17",
            "state_b": "M18",
            "state_a_readable": "M17: ISLR+COL5A2+ myCAF",
            "state_b_readable": "M18: S100B+TYR+ Melanoc.",
            "compartment_a": "stromal_mesenchymal",
            "compartment_b": "malignant_epithelial",
            "dominant_regime": "log_ratio_balance",
            "dominant_direction": "non-responder",
            "z_acat_signed": -2.40,
            "pvalue_acat_iut": 0.016,
            "fdr_acat_iut": 0.09,
            "emp_perm_pvalue": 0.019,
            "min_loco_z_acat": 0.12,
            "is_significant": True,
        },
    ]
    synth_df = pl.DataFrame(synth_rows)
    clustered, dend_segs = cluster_states_within_compartments(
        synth_df, include_malignant_epithelial=False
    )
    codes = [st.state_code for st in clustered]
    assert "M18" not in codes
    assert len(codes) == 8
    assert len(dend_segs) > 0

    # All Immune Cells precede all Stromal / Endothelial states
    supergroups = [st.supergroup for st in clustered]
    assert supergroups[:6] == [TmeSuperGroup.IMMUNE_CELLS] * 6
    assert supergroups[6:] == [TmeSuperGroup.STROMAL_ENDOTHELIAL] * 2

    # Correlated T_NK states M04 and M08 cluster next to each other
    assert abs(codes.index("M04") - codes.index("M08")) == 1

    chart = build_pan_cohort_z_acat_heatmap_figure(synth_df, fdr_threshold=0.10)
    spec = chart.to_dict()
    assert "layer" in spec
    assert len(spec["layer"]) >= 12

    core_codes = select_core_hub_state_codes(synth_df, fdr_threshold=0.10)
    # M04 (deg 2), M17 (deg 3), M08 (perm<=0.01 & |Z|>=3.25), M30 (perm<=0.01 & |Z|>=3.25), M10 (perm<=0.01 & |Z|>=3.25)
    assert core_codes == frozenset({"M04", "M08", "M10", "M17", "M30"})
    assert "M01" not in core_codes
    assert "M14" not in core_codes
    assert "M18" not in core_codes

    core_chart = build_pan_cohort_core_hubs_z_acat_heatmap_figure(synth_df, fdr_threshold=0.10)
    core_spec = core_chart.to_dict()
    assert "layer" in core_spec
    assert len(core_spec["layer"]) >= 14


@settings(max_examples=20, deadline=3000)
@given(
    swap_nhoods=st.booleans(),
    regime=st.sampled_from(
        [
            "log_ratio_balance",
            "co_abundance_product",
            "quadrant:++",
            "quadrant:+-",
            "quadrant:-+",
            "quadrant:--",
        ]
    ),
    direction=st.sampled_from(["responder", "non-responder"]),
)
def test_prop_oriented_super_nhoods_and_regime_polarity_invariance(
    swap_nhoods: bool,
    regime: str,
    direction: str,
) -> None:
    """Verifies endpoint resolution & visual polarity stay locked to `(state_a, state_b)` under nhood swaps."""
    super_codes = ("M30", "M17", "M04", "M42")
    raw_u, raw_v = (0, 2) if swap_nhoods else (2, 0)
    nh_a, nh_b = _resolve_oriented_super_nhoods(
        nhood_u=raw_u,
        nhood_v=raw_v,
        state_a="M04",
        state_b="M30",
        super_meta_codes=super_codes,
    )
    assert super_codes[nh_a] == "M04"
    assert super_codes[nh_b] == "M30"

    pol_res = _map_dominant_regime_to_visual_polarity(regime, direction)
    assert isinstance(pol_res, Success)
    vis_reg, pol_a, pol_b = pol_res.unwrap()
    assert len(vis_reg.value) > 0
    if regime == "log_ratio_balance" and direction == "responder":
        assert pol_a.value != pol_b.value


def test_broad_compartment_pair_umap_grid_and_standalone_composites() -> None:
    """Verifies active broad compartment-pair stratification, dual-tier chords, and Vega-Lite specs."""
    rng = np.random.default_rng(314)
    super_codes = ("M04", "M08", "M17", "M30", "M16", "M18")
    read_lbls = (
        "M04: GNLY+PRF1+ Eff. T",
        "M08: CCR4+SOCS1+ Th2/Treg",
        "M17: ISLR+COL5A2+ myCAF",
        "M30: CR2+PAX5+ TLS B",
        "M16: VEGFA+KDR+ Endo.",
        "M18: S100B+TYR+ Melanoc.",
    )
    comp_per_nh = (
        "t_nk",
        "t_nk",
        "stromal_mesenchymal",
        "b_lineage",
        "endothelial_vascular",
        "malignant_epithelial",
    )
    n_nh = len(super_codes)
    cells_per_nh = 20
    nhoods_mat = _make_disjoint_incidence(n_nh, cells_per_nhood=cells_per_nh)
    n_cells = int(nhoods_mat.shape[0])

    island_centers = np.array(
        [
            [2.0, -5.0],
            [2.8, -4.4],
            [-2.0, 8.8],
            [3.2, 5.2],
            [-3.8, 7.0],
            [3.0, 11.5],
        ],
        dtype=np.float64,
    )
    umap_coords = np.vstack(
        [
            island_centers[k] + rng.normal(scale=0.35, size=(cells_per_nh, 2))
            for k in range(n_nh)
        ]
    )
    anchor_comps = tuple(comp_per_nh[k] for k in range(n_nh) for _ in range(cells_per_nh))
    anchor_lbls = tuple(read_lbls[k] for k in range(n_nh) for _ in range(cells_per_nh))

    n_evals = 24
    counts_super = rng.poisson(lam=12.0, size=(n_nh, n_evals)).astype(np.float64)
    libs = np.full(n_evals, 1000.0, dtype=np.float64)
    cohorts = tuple(["COH_A"] * 12 + ["COH_B"] * 12)
    strata = tuple(["COH_A::main"] * 12 + ["COH_B::main"] * 12)
    resp = tuple(["responder", "non-responder"] * (n_evals // 2))

    pan_unit = PreparedSuperUnit(
        name="SYNTH_PAN",
        nhoods_super=nhoods_mat,
        counts_super=counts_super,
        library_sizes=libs,
        cohort_labels=cohorts,
        strata_labels=strata,
        super_meta_codes=super_codes,
        super_readable_labels=read_lbls,
    )

    verified_df = pl.DataFrame(
        [
            # T_NK <-> CAF: 1 verified FDR < 0.10 + 1 perm-only p < 0.05
            {
                "state_a": "M08",
                "state_b": "M17",
                "state_a_readable": read_lbls[1],
                "state_b_readable": read_lbls[2],
                "compartment_a": "T_NK",
                "compartment_b": "STROMAL_MESENCHYMAL",
                "nhood_u": 2,  # Intentionally swapped nhood_u / nhood_v!
                "nhood_v": 1,
                "dominant_regime": "log_ratio_balance",
                "dominant_direction": "responder",
                "dominant_effect": 0.62,
                "z_acat_signed": 4.14,
                "z_acat_abs": 4.14,
                "pvalue_acat_iut": 3.5e-5,
                "fdr_acat_iut": 0.007,
                "emp_perm_pvalue": 0.001,
                "min_loco_z_acat": 3.41,
                "max_loco_z_acat": 4.20,
                "prior_weight": 2.0,
                "is_significant": True,
                "is_loco_stable": True,
            },
            {
                "state_a": "M04",
                "state_b": "M17",
                "state_a_readable": read_lbls[0],
                "state_b_readable": read_lbls[2],
                "compartment_a": "T_NK",
                "compartment_b": "STROMAL_MESENCHYMAL",
                "nhood_u": 0,
                "nhood_v": 2,
                "dominant_regime": "quadrant:+-",
                "dominant_direction": "responder",
                "dominant_effect": 1.45,
                "z_acat_signed": 2.28,
                "z_acat_abs": 2.28,
                "pvalue_acat_iut": 0.022,
                "fdr_acat_iut": 0.24,
                "emp_perm_pvalue": 0.026,
                "min_loco_z_acat": 1.75,
                "max_loco_z_acat": 2.40,
                "prior_weight": 1.2,
                "is_significant": False,
                "is_loco_stable": False,
            },
            # B <-> Endothelial: 0 FDR, 1 perm-only p < 0.05
            {
                "state_a": "M16",
                "state_b": "M30",
                "state_a_readable": read_lbls[4],
                "state_b_readable": read_lbls[3],
                "compartment_a": "ENDOTHELIAL_VASCULAR",
                "compartment_b": "B_LINEAGE",
                "nhood_u": 3,
                "nhood_v": 4,
                "dominant_regime": "log_ratio_balance",
                "dominant_direction": "non-responder",
                "dominant_effect": -0.48,
                "z_acat_signed": -2.31,
                "z_acat_abs": 2.31,
                "pvalue_acat_iut": 0.020,
                "fdr_acat_iut": 0.24,
                "emp_perm_pvalue": 0.0207,
                "min_loco_z_acat": 1.80,
                "max_loco_z_acat": 2.45,
                "prior_weight": 1.0,
                "is_significant": False,
                "is_loco_stable": False,
            },
            # Malignant epithelial hit must be excluded from host broad cell-type UMAPs
            {
                "state_a": "M17",
                "state_b": "M18",
                "state_a_readable": read_lbls[2],
                "state_b_readable": read_lbls[5],
                "compartment_a": "STROMAL_MESENCHYMAL",
                "compartment_b": "MALIGNANT_EPITHELIAL",
                "nhood_u": 2,
                "nhood_v": 5,
                "dominant_regime": "log_ratio_balance",
                "dominant_direction": "non-responder",
                "dominant_effect": -0.55,
                "z_acat_signed": -3.10,
                "z_acat_abs": 3.10,
                "pvalue_acat_iut": 0.002,
                "fdr_acat_iut": 0.04,
                "emp_perm_pvalue": 0.005,
                "min_loco_z_acat": 2.10,
                "max_loco_z_acat": 3.20,
                "prior_weight": 1.0,
                "is_significant": True,
                "is_loco_stable": True,
            },
        ]
    )

    specs = enumerate_active_broad_compartment_pairs(verified_df, fdr_threshold=0.10, perm_alpha=0.05)
    assert len(specs) == 2
    assert specs[0].slug == "t_nk_vs_caf"
    assert specs[0].n_fdr_pairs == 1
    assert specs[0].n_perm_only_pairs == 1
    assert specs[1].slug == "b_cells_vs_endothelial"
    assert specs[1].n_fdr_pairs == 0
    assert specs[1].n_perm_only_pairs == 1

    cell_bundles = []
    for sp_item in specs:
        cb_res = build_compartment_pair_umap_cell_bundle(
            spec=sp_item,
            pan_unit=pan_unit,
            umap_coords=umap_coords,
            anchor_compartments=anchor_comps,
            verified_pan_df=verified_df,
            fdr_threshold=0.10,
            perm_alpha=0.05,
        )
        assert isinstance(cb_res, Success)
        cb = cb_res.unwrap()
        n_t1 = cb.tier1_arcs_df["arc_id"].n_unique() if not cb.tier1_arcs_df.is_empty() else 0
        n_t2 = cb.tier2_arcs_df["arc_id"].n_unique() if not cb.tier2_arcs_df.is_empty() else 0
        assert n_t1 == sp_item.n_fdr_pairs
        assert n_t2 == sp_item.n_perm_only_pairs
        assert cb.cells_df.height == n_cells
        cell_bundles.append(cb)

    grid_res = build_pan_cohort_compartment_pair_umap_grid(cell_bundles, fdr_threshold=0.10)
    assert isinstance(grid_res, Success)
    grid_dict = grid_res.unwrap().to_dict()
    assert grid_dict["config"]["legend"]["disable"] is True
    assert "feDropShadow" not in str(grid_dict)

    pw_res = build_compartment_pair_pairwise_umap_bundle(
        spec=specs[1],
        pan_unit=pan_unit,
        umap_coords=umap_coords,
        anchor_cell_labels=anchor_lbls,
        global_response_labels=resp,
        verified_pan_df=verified_df,
        fdr_threshold=0.10,
        perm_alpha=0.05,
    )
    assert isinstance(pw_res, Success)
    pw_b = pw_res.unwrap()
    assert pw_b.spotlight_meta.has_pair is True

    comp_chart = build_standalone_compartment_pair_umap_composite(
        cell_bundle=cell_bundles[1],
        pw_bundle=pw_b,
        selected_regime=HarmonyRegime.I2_HARMONY_COHORT_SUBSTUDY,
        n_biopsies=n_evals,
        n_strata=2,
        fdr_threshold=0.10,
    )
    comp_dict = comp_chart.to_dict()
    assert "vconcat" in comp_dict
    assert "feDropShadow" not in str(comp_dict)

