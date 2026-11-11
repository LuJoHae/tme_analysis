"""Property-based and statistical calibration tests for donor mixing QC, cross-compartment
pair simplification, and stratified IUT permutation + LODO/LOSO jackknife verification.
"""

from __future__ import annotations

from hypothesis import given, settings, strategies as st
import numpy as np
import polars as pl
from returns.maybe import Nothing, Some
from returns.result import Success
from scipy import sparse as sp

from milopy.pairwise import (
    ContinuousAxis,
    EnrichmentDirection,
    MixingQCConfig,
    OverlapConfig,
    PairScope,
    PairScopeFilter,
    PairwiseConfig,
    PairwiseMethod,
    PresenceConfig,
    PresenceMode,
    QuadrantState,
    RefinementConfig,
    ScreenConfig,
    TmeCompartment,
    analyze_nhood_pairs,
    calibrate_candidates_stratified_iut_permutation,
    classify_tme_compartment,
    compute_nhood_mixing_metrics,
    consolidate_state_pair_bundles,
    deduplicate_continuous_pairs_to_lead_bundles,
    deduplicate_pairs_to_lead_bundles,
    estimate_nhood_presence,
    extract_inferred_substudy_strata,
)


def test_classify_tme_compartment_resolves_gondal_ontology_synonyms() -> None:
    """Verifies ontology synonyms and Meta42 programs map cleanly without cross-lineage confusion."""
    for t_label in (
        "T cell",
        "lymphocyte",
        "central memory CD8-positive, alpha-beta T cell",
        "T follicular helper cell",
        "natural killer cell",
        "M03",
        "M03: CD200+ELMO1+ Tfh T",
        "M05",
        "M13",
    ):
        assert classify_tme_compartment(t_label) == TmeCompartment.T_NK

    assert classify_tme_compartment("lymphocyte of B lineage") == TmeCompartment.B_LINEAGE
    assert classify_tme_compartment("M30: CR2+PAX5+ TLS B") == TmeCompartment.B_LINEAGE
    assert classify_tme_compartment("IgG plasma cell") == TmeCompartment.PLASMA
    assert classify_tme_compartment("M29") == TmeCompartment.PLASMA
    assert classify_tme_compartment("plasmacytoid dendritic cell, human") == TmeCompartment.MYELOID_MAST
    assert classify_tme_compartment("M41: LAMP3+CCR7+ mregDC") == TmeCompartment.MYELOID_MAST
    assert classify_tme_compartment("myofibroblast cell") == TmeCompartment.STROMAL_MESENCHYMAL
    assert classify_tme_compartment("M17: ISLR+COL5A2+ myCAF") == TmeCompartment.STROMAL_MESENCHYMAL
    assert classify_tme_compartment("vein endothelial cell") == TmeCompartment.ENDOTHELIAL_VASCULAR
    assert classify_tme_compartment("malignant cell") == TmeCompartment.MALIGNANT_EPITHELIAL
    assert classify_tme_compartment("Unspecified") == TmeCompartment.OTHER
    assert classify_tme_compartment("3") == TmeCompartment.OTHER


def test_extract_inferred_substudy_strata_compendium_vs_single_study() -> None:
    """Detects `<PMID>ggg<Patient>` and `<Study>.<Patient>` strata while ignoring single-study IDs."""
    gondal_ids = (
        "27124452gggMel58__non-responder",
        "27124452gggMel60__non-responder",
        "33958794gggP01__responder",
        "33958794gggP02__non-responder",
        "31359002gggPt4__responder",
        "31359002gggPt9__non-responder",
    )
    res_g = extract_inferred_substudy_strata(gondal_ids)
    assert isinstance(res_g, Some)
    assert res_g.unwrap() == (
        "27124452",
        "27124452",
        "33958794",
        "33958794",
        "31359002",
        "31359002",
    )

    crc_ids = (
        "Che_2021_Cell_Discov.COL15",
        "Che_2021_Cell_Discov.COL17",
        "Chen_2024_Cancer_Cell.P01",
        "Chen_2024_Cancer_Cell.P02",
    )
    res_c = extract_inferred_substudy_strata(crc_ids)
    assert isinstance(res_c, Some)
    assert res_c.unwrap() == (
        "Che_2021_Cell_Discov",
        "Che_2021_Cell_Discov",
        "Chen_2024_Cancer_Cell",
        "Chen_2024_Cancer_Cell",
    )

    single_ids = ("P01", "P02", "P03", "P04", "P05", "P06")
    assert extract_inferred_substudy_strata(single_ids) == Nothing


@settings(max_examples=40, deadline=None)
@given(
    n_nhoods=st.integers(min_value=2, max_value=14),
    n_donors=st.integers(min_value=4, max_value=20),
    seed=st.integers(min_value=0, max_value=10000),
)
def test_nhood_mixing_metrics_mathematical_invariants(
    n_nhoods: int,
    n_donors: int,
    seed: int,
) -> None:
    """Property test: 1 <= LISI_pt <= n_distinct_donors <= J, uniform == J, permutation-invariant."""
    rng = np.random.default_rng(seed)
    counts = rng.poisson(lam=2.5, size=(n_nhoods, n_donors)).astype(np.float64)
    # Force row 0 to be strictly uniform across all J donors, row 1 strictly single-donor private
    counts[0, :] = 5.0
    counts[1, :] = 0.0
    counts[1, 0] = 25.0

    nhoods_mat = sp.eye(n_nhoods, format="csc", dtype=np.float64)
    sids = tuple(f"P{j:02d}" for j in range(n_donors))
    cfg = MixingQCConfig(
        enabled=True,
        min_distinct_donors=3,
        min_effective_donors_lisi=2.5,
        max_single_donor_dominance=0.50,
        adaptive_small_cohort_scaling=False,
    )

    res = compute_nhood_mixing_metrics(nhoods_mat, counts, sids, config=cfg)
    assert isinstance(res, Success)
    df = res.unwrap()

    lisi = df["effective_donors_lisi"].to_numpy()
    nd = df["n_distinct_donors"].to_numpy()
    dom = df["max_single_donor_fraction"].to_numpy()
    ent = df["shannon_donor_entropy_norm"].to_numpy()
    nc = df["n_cells"].to_numpy()

    pos = nc > 0.0
    assert np.all(lisi[pos] >= 1.0 - 1e-9)
    assert np.all(lisi[pos] <= nd[pos].astype(np.float64) + 1e-9)
    assert np.all((ent >= -1e-9) & (ent <= 1.0 + 1e-9))
    assert np.all(dom[pos] * nd[pos].astype(np.float64) >= 1.0 - 1e-9)

    # Row 0 (uniform across all J donors) -> LISI == J, entropy == 1, dominance == 1/J
    assert np.isclose(float(lisi[0]), float(n_donors), atol=1e-9)
    assert np.isclose(float(ent[0]), 1.0, atol=1e-9)
    assert np.isclose(float(dom[0]), 1.0 / float(n_donors), atol=1e-9)
    assert bool(df["passes_mixing_qc"][0]) is True

    # Row 1 (single-donor private clone) -> LISI == 1, dominance == 1, fails QC
    assert np.isclose(float(lisi[1]), 1.0, atol=1e-9)
    assert np.isclose(float(dom[1]), 1.0, atol=1e-9)
    assert bool(df["passes_mixing_qc"][1]) is False

    # Permuting donor columns preserves all neighborhood-level mixing metrics identically
    perm = rng.permutation(n_donors)
    perm_sids = tuple(sids[int(p)] for p in perm)
    res_perm = compute_nhood_mixing_metrics(nhoods_mat, counts[:, perm], perm_sids, config=cfg)
    assert isinstance(res_perm, Success)
    df_perm = res_perm.unwrap()
    assert np.allclose(lisi, df_perm["effective_donors_lisi"].to_numpy(), atol=1e-12)
    assert np.array_equal(
        df["passes_mixing_qc"].to_numpy(),
        df_perm["passes_mixing_qc"].to_numpy(),
    )


def test_state_pair_bundle_consolidation_swap_symmetry_and_idempotence() -> None:
    """Verifies (M03, M17, +- ) and (M17, M03, -+) collapse into 1 canonical bundle idempotently."""
    pairs_df = pl.DataFrame(
        {
            "nhood_u": [1, 4, 8],
            "nhood_v": [10, 2, 12],
            "quadrant": ["+-", "-+", "++"],
            "direction": ["responder", "responder", "responder"],
            "cell_type_u": ["M03", "M17", "M03"],
            "cell_type_v": ["M17", "M03", "M05"],
            "log_or_quadrant": [2.4, 2.1, 1.9],
            "pvalue_iut": [1e-4, 4e-4, 2e-4],
            "fdr_iut": [0.01, 0.03, 0.02],
            "is_significant": [True, True, True],
        }
    )

    b_res = consolidate_state_pair_bundles(pairs_df, fdr_threshold=0.10)
    assert isinstance(b_res, Success)
    bundles = b_res.unwrap()
    # Two unique bundles: (M03, M17, +-) with 2 sub-pairs, and intra-T (M03, M05, ++) with 1 sub-pair
    assert bundles.height == 2
    m03_m17 = bundles.filter((pl.col("state_a") == "M03") & (pl.col("state_b") == "M17"))
    assert m03_m17.height == 1
    assert int(m03_m17["n_supporting_subpairs"][0]) == 2
    assert str(m03_m17["canonical_quadrant"][0]) == QuadrantState.FIRST_ONLY.value
    assert bool(m03_m17["is_cross_compartment"][0]) is True

    # Deduplicating with require_cross_compartment=True drops M03<->M05 and keeps 1 M03<->M17 lead chord
    d1_res = deduplicate_pairs_to_lead_bundles(pairs_df, require_cross_compartment=True)
    assert isinstance(d1_res, Success)
    d1 = d1_res.unwrap()
    assert d1.height == 1
    assert int(d1["nhood_u"][0]) == 1
    assert int(d1["nhood_v"][0]) == 10

    # Strict idempotence
    d2_res = deduplicate_pairs_to_lead_bundles(d1, require_cross_compartment=True)
    assert isinstance(d2_res, Success)
    assert d2_res.unwrap().equals(d1)

    # Multi-regime state pair collapses to strictly 1 chord under one_per_state_pair=True
    multi_regime_df = pl.concat(
        [
            pairs_df,
            pl.DataFrame(
                {
                    "nhood_u": [3],
                    "nhood_v": [11],
                    "quadrant": ["--"],
                    "direction": ["non-responder"],
                    "cell_type_u": ["M03"],
                    "cell_type_v": ["M17"],
                    "log_or_quadrant": [2.0],
                    "pvalue_iut": [3e-4],
                    "fdr_iut": [0.025],
                    "is_significant": [True],
                }
            ),
        ]
    )
    assert deduplicate_pairs_to_lead_bundles(
        multi_regime_df, require_cross_compartment=True, one_per_state_pair=False
    ).unwrap().height == 2
    d_one = deduplicate_pairs_to_lead_bundles(
        multi_regime_df, require_cross_compartment=True, one_per_state_pair=True
    ).unwrap()
    assert d_one.height == 1
    assert int(d_one["nhood_u"][0]) == 1


@settings(max_examples=40, deadline=None)
@given(
    p_best=st.floats(min_value=1e-6, max_value=1e-3, allow_nan=False, allow_infinity=False),
    delta_p=st.floats(min_value=1e-4, max_value=2e-2, allow_nan=False, allow_infinity=False),
    swap_lead_endpoints=st.booleans(),
)
def test_prop_continuous_bundle_deduplication_swap_sign_and_idempotence(
    p_best: float,
    delta_p: float,
    swap_lead_endpoints: bool,
) -> None:
    """Property test: antisymmetric log-ratio sign flip, symmetric co-abundance, cross-compartment gate & idempotence."""
    p_second = min(0.08, p_best + delta_p)
    lead_u_lbl = "M17: ISLR+COL5A2+ myCAF" if swap_lead_endpoints else "M03: CD200+ELMO1+ Tfh T"
    lead_v_lbl = "M03: CD200+ELMO1+ Tfh T" if swap_lead_endpoints else "M17: ISLR+COL5A2+ myCAF"
    # Antisymmetric log-ratio: log(u/v) enriched in R <=> log(v/u) enriched in NR!
    lead_dir = (
        EnrichmentDirection.NON_RESPONDER.value
        if swap_lead_endpoints
        else EnrichmentDirection.RESPONDER.value
    )
    second_u_lbl = "M03: CD200+ELMO1+ Tfh T" if swap_lead_endpoints else "M17: ISLR+COL5A2+ myCAF"
    second_v_lbl = "M17: ISLR+COL5A2+ myCAF" if swap_lead_endpoints else "M03: CD200+ELMO1+ Tfh T"
    second_dir = (
        EnrichmentDirection.RESPONDER.value
        if swap_lead_endpoints
        else EnrichmentDirection.NON_RESPONDER.value
    )

    cont_df = pl.DataFrame(
        {
            "nhood_u": [2, 9, 4, 11, 6],
            "nhood_v": [8, 3, 10, 5, 7],
            "axis": [
                ContinuousAxis.LOG_RATIO_BALANCE.value,
                ContinuousAxis.LOG_RATIO_BALANCE.value,
                ContinuousAxis.CO_ABUNDANCE_PRODUCT.value,
                ContinuousAxis.CO_ABUNDANCE_PRODUCT.value,
                ContinuousAxis.LOG_RATIO_BALANCE.value,
            ],
            "direction": [
                lead_dir,
                second_dir,
                EnrichmentDirection.RESPONDER.value,
                EnrichmentDirection.RESPONDER.value,
                EnrichmentDirection.RESPONDER.value,
            ],
            "cell_type_u": [
                lead_u_lbl,
                second_u_lbl,
                "M03: CD200+ELMO1+ Tfh T",
                "M17: ISLR+COL5A2+ myCAF",
                "M03: CD200+ELMO1+ Tfh T",
            ],
            "cell_type_v": [
                lead_v_lbl,
                second_v_lbl,
                "M17: ISLR+COL5A2+ myCAF",
                "M03: CD200+ELMO1+ Tfh T",
                "M05: CD96+TRAF3+ Mem. T",
            ],
            "effect_size": [-2.8 if swap_lead_endpoints else 2.8, 2.3, 2.5, 2.1, 3.0],
            "pvalue_iut": [p_best, p_second, p_second * 1.1, p_second * 1.4, p_best * 0.5],
            "fdr_iut": [0.01, 0.04, 0.045, 0.06, 0.005],
            "is_significant": [True, True, True, True, True],
        }
    )

    # 1. With one_per_state_pair=False & require_cross_compartment=True:
    #    - Intra-T_NK (M03 <-> M05) row is dropped
    #    - Antisymmetric log_ratio_balance pair collapses (2 rows -> 1 row, lead nhood_u=2)
    #    - Symmetric co_abundance_product pair collapses (2 rows -> 1 row, lead nhood_u=4)
    res_by_axis = deduplicate_continuous_pairs_to_lead_bundles(
        cont_df,
        require_cross_compartment=True,
        one_per_state_pair=False,
    )
    assert isinstance(res_by_axis, Success)
    df_by_axis = res_by_axis.unwrap()
    assert df_by_axis.height == 2
    assert set(df_by_axis["nhood_u"].to_list()) == {2, 4}

    # 2. With one_per_state_pair=True & require_cross_compartment=True:
    #    Both axes of (M03, M17) collapse to strictly the single strongest lead chord (nhood_u=2)
    res_one = deduplicate_continuous_pairs_to_lead_bundles(
        cont_df,
        require_cross_compartment=True,
        one_per_state_pair=True,
    )
    assert isinstance(res_one, Success)
    df_one = res_one.unwrap()
    assert df_one.height == 1
    assert int(df_one["nhood_u"][0]) == 2
    assert int(df_one["nhood_v"][0]) == 8

    # 3. Strict idempotence under repeated application
    res_idem = deduplicate_continuous_pairs_to_lead_bundles(
        df_one,
        require_cross_compartment=True,
        one_per_state_pair=True,
    )
    assert isinstance(res_idem, Success)
    assert res_idem.unwrap().equals(df_one)


def test_confounded_substudy_simpsons_paradox_blocked_by_cmh_and_mixing_qc() -> None:
    """Synthetic multi-study compendium where responder rate varies 0% -> 80% across sub-studies.

    Neighborhoods are purely study-associated (zero within-study R vs NR effect), plus one
    0%-Responder sub-study (`27124452`). Unstratified screening yields false positives, whereas
    CMH sub-study stratification + donor mixing QC produces 0 false discoveries.
    """
    rng = np.random.default_rng(123)
    # 4 sub-studies: S0 (20 NR, 0 R), S1 (4 R, 41 NR), S2 (22 R, 23 NR), S3 (41 R, 4 NR) -> J = 155
    substudies = (
        ["27124452"] * 20
        + ["38181739"] * 45
        + ["33958794"] * 45
        + ["33711272"] * 45
    )
    responses = (
        ["non-responder"] * 20
        + ["responder"] * 4
        + ["non-responder"] * 41
        + ["responder"] * 22
        + ["non-responder"] * 23
        + ["responder"] * 41
        + ["non-responder"] * 4
    )
    sample_ids = tuple(f"{substudies[j]}gggPt{j:03d}" for j in range(len(substudies)))
    n_samples = len(sample_ids)

    # 12 disjoint neighborhoods across compartments:
    # Nhoods (0, 2, 6, 8) are protocol-enriched in high-responder study 33711272 (90% regardless of R/NR!)
    # and absent/rare in low-responder studies (5% regardless of R/NR!) -> pure Simpson's paradox!
    n_nhoods = 12
    counts = np.zeros((n_nhoods, n_samples), dtype=np.float64)
    for j in range(n_samples):
        st_id = substudies[j]
        high_prob = 0.92 if st_id == "33711272" else (0.35 if st_id == "33958794" else 0.04)
        for v in (0, 2, 6, 8):
            if rng.random() < high_prob:
                counts[v, j] = float(rng.integers(5, 14))
        for v in (3, 4, 5, 7, 9, 10, 11):
            if rng.random() < 0.45:
                counts[v, j] = float(rng.integers(3, 9))

    # Make nhood 1 a single-patient clonal spike in one responder of 33711272
    counts[1, :] = 0.0
    counts[1, 80] = 60.0

    nhoods_mat = sp.eye(n_nhoods, format="csc", dtype=np.float64)
    nhood_labels = ("M03", "M03", "M04", "M05", "M13", "M15", "M17", "M23", "M10", "M14", "M30", "M29")
    lib_sizes = tuple([200.0] * n_samples)

    unstrat_cfg = PairwiseConfig(
        method=PairwiseMethod.QUADRANT_2X2,
        overlap=OverlapConfig(scope_filter=PairScopeFilter.CROSS_COMPARTMENT_ONLY),
        presence=PresenceConfig(mode=PresenceMode.HARD_DEPTH_GATED, min_cell_count=1, min_fraction=0.0),
        refinement=RefinementConfig(enabled=False),
    )
    res_unstrat = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=responses,
        sample_library_sizes=Some(lib_sizes),
        nhood_group_labels=Some(nhood_labels),
        sample_ids=Some(sample_ids),
        config=unstrat_cfg,
    )
    assert isinstance(res_unstrat, Success)
    # Unstratified analysis is confounded by Simpson's paradox across studies (0% vs 80% R)
    assert res_unstrat.unwrap().n_significant_quadrant_pairs >= 1

    # Now apply inferred PMID strata + MixingQCConfig -> 0 false discoveries!
    inferred_strata = extract_inferred_substudy_strata(sample_ids)
    assert isinstance(inferred_strata, Some)
    strat_cfg = PairwiseConfig(
        method=PairwiseMethod.QUADRANT_2X2,
        mixing_qc=MixingQCConfig(
            enabled=True,
            min_distinct_donors=3,
            min_effective_donors_lisi=2.5,
            max_single_donor_dominance=0.50,
        ),
        overlap=OverlapConfig(scope_filter=PairScopeFilter.CROSS_COMPARTMENT_ONLY),
        presence=PresenceConfig(mode=PresenceMode.HARD_DEPTH_GATED, min_cell_count=1, min_fraction=0.0),
        refinement=RefinementConfig(enabled=False),
    )
    res_strat = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=responses,
        sample_library_sizes=Some(lib_sizes),
        strata_labels=Some(inferred_strata.unwrap()),
        nhood_group_labels=Some(nhood_labels),
        sample_ids=Some(sample_ids),
        config=strat_cfg,
    )
    assert isinstance(res_strat, Success)
    strat_out = res_strat.unwrap()
    assert strat_out.n_significant_quadrant_pairs == 0
    # Clonal single-donor nhood 1 was removed by MixingQCConfig before pair formation
    evaluated_nodes = set(strat_out.pairs_df["nhood_u"].to_list()) | set(
        strat_out.pairs_df["nhood_v"].to_list()
    )
    assert 1 not in evaluated_nodes


def test_four_layer_pair_significance_recovers_true_gate_and_blocks_singleton_piggybacking() -> None:
    """Verifies true cross-compartment +- gate passes all 4 layers while H_10 singleton piggybacking fails."""
    rng = np.random.default_rng(2026)
    # 2 balanced strata x 40 donors = 80 donors (40 R, 40 NR)
    n_samples = 80
    strata = tuple(["S1"] * 40 + ["S2"] * 40)
    strata_codes = np.array([0] * 40 + [1] * 40, dtype=np.int64)
    responses = tuple((["responder"] * 20 + ["non-responder"] * 20) * 2)
    y_bin = np.array([1.0 if r == "responder" else 0.0 for r in responses], dtype=np.float64)

    # 8 disjoint cross-compartment neighborhoods
    # Pair (0: M03 T_NK, 1: M17 STROMAL_MESENCHYMAL): genuine Responder +- gate!
    #   Responders: 85% are (u+, v-); Non-responders: split evenly across (u+, v+), (u-, v+), (u-, v-)
    #Singleton 2 (M13 T_NK): strongly R-enriched alone, paired with pure noise 3..7 (H_10 piggybacking!)
    n_nhoods = 8
    counts = np.zeros((n_nhoods, n_samples), dtype=np.float64)
    for j in range(n_samples):
        is_r = y_bin[j] > 0.5
        if is_r:
            if rng.random() < 0.88:
                counts[0, j] = 8.0
                counts[1, j] = 0.0
            else:
                counts[0, j] = 8.0 if rng.random() < 0.5 else 0.0
                counts[1, j] = 8.0
            if rng.random() < 0.85:
                counts[2, j] = 7.0
        else:
            roll = rng.random()
            if roll < 0.42:
                counts[0, j] = 8.0
                counts[1, j] = 8.0
            elif roll < 0.78:
                counts[0, j] = 0.0
                counts[1, j] = 8.0
            elif roll < 0.86:
                counts[0, j] = 8.0
                counts[1, j] = 0.0
            if rng.random() < 0.15:
                counts[2, j] = 7.0

        for v in range(3, n_nhoods):
            if rng.random() < 0.50:
                counts[v, j] = 6.0

    nhoods_mat = sp.eye(n_nhoods, format="csc", dtype=np.float64)
    nhood_labels = ("M03", "M17", "M13", "M30", "M10", "M14", "M16", "M29")
    lib_sizes = np.full(n_samples, 200.0, dtype=np.float64)

    cfg = PairwiseConfig(
        method=PairwiseMethod.QUADRANT_2X2,
        mixing_qc=MixingQCConfig(enabled=True, min_distinct_donors=3, min_effective_donors_lisi=2.5),
        overlap=OverlapConfig(scope_filter=PairScopeFilter.CROSS_COMPARTMENT_ONLY),
        presence=PresenceConfig(mode=PresenceMode.HARD_DEPTH_GATED, min_cell_count=1, min_fraction=0.0),
        screen=ScreenConfig(fdr_threshold=0.10),
        refinement=RefinementConfig(enabled=True, max_refine_pairs=50),
    )
    res_wrapped = analyze_nhood_pairs(
        nhoods_mat=nhoods_mat,
        count_matrix=counts,
        response_labels=responses,
        sample_library_sizes=Some(tuple(lib_sizes.tolist())),
        strata_labels=Some(strata),
        nhood_group_labels=Some(nhood_labels),
        config=cfg,
    )
    assert isinstance(res_wrapped, Success)
    res = res_wrapped.unwrap()
    assert set(res.pairs_df["pair_scope"].to_list()) == {PairScope.CROSS_COMPARTMENT.value}

    # True (0, 1, +-) gate must be significant
    true_hit = res.pairs_df.filter(
        (pl.col("nhood_u") == 0)
        & (pl.col("nhood_v") == 1)
        & (pl.col("quadrant") == QuadrantState.FIRST_ONLY.value)
    )
    assert true_hit.height == 1
    assert bool(true_hit["is_significant"][0]) is True

    # Singleton piggybacking H_10 pairs: true singleton nhood 2 paired with noise partners 3..7 must NOT pass IUT!
    piggy_sig = res.pairs_df.filter(
        (pl.col("nhood_u") == 2)
        & pl.col("nhood_v").is_in([3, 4, 5, 6, 7])
        & pl.col("is_significant")
    )
    assert piggy_sig.height == 0

    # Run Layer 3 & Layer 4 stratified IUT permutation + LODO/LOSO jackknife verification
    pres = estimate_nhood_presence(counts, Some(lib_sizes), cfg.presence).unwrap()
    cal_res = calibrate_candidates_stratified_iut_permutation(
        candidates_df=res.pairs_df.head(20),
        presence=pres,
        response_binary=y_bin,
        strata_codes=Some(strata_codes),
        n_permutations=500,
        seed=42,
        fdr_threshold=0.10,
    )
    assert isinstance(cal_res, Success)
    cal_df = cal_res.unwrap()
    ver_hit = cal_df.filter(
        (pl.col("nhood_u") == 0)
        & (pl.col("nhood_v") == 1)
        & (pl.col("quadrant") == QuadrantState.FIRST_ONLY.value)
    )
    assert ver_hit.height == 1
    assert float(ver_hit["pvalue_perm_iut"][0]) < 0.01
    assert float(ver_hit["fwer_perm_iut"][0]) < 0.05
    assert float(ver_hit["lodo_min_z_iut"][0]) > 2.0
    assert float(ver_hit["loso_sign_consistency"][0]) == 1.0
    assert bool(ver_hit["is_statistically_verified"][0]) is True
