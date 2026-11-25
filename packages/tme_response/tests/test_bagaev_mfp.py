"""Unit and property-based tests for Bagaev et al. (2021) MFP TME subtyping and response analysis."""

from __future__ import annotations

from hypothesis import given, settings, strategies as st
import numpy as np
import polars as pl
from returns.maybe import Some
from returns.result import Success

from tme_datasets import BAGAEV_MFP_FUNCTIONAL_GROUPS, BAGAEV_MFP_SIGNATURE_ORDER
from tme_response.models.bagaev_mfp import (
    BagaevKnnConfig,
    BagaevLouvainConfig,
    DeconvClusterConfig,
    _louvain_modularity_partition,
    analyze_subtype_deconvolution_enrichment,
    assign_clusters_to_bagaev_subtypes,
    benchmark_bagaev_vs_deconv_cv,
    classify_bagaev_knn,
    cluster_bagaev_louvain,
    cluster_deconvoluted_samples,
    compare_partition_concordance,
    compute_cmh_k_by_2_p_value,
    compute_deconv_subtype_correlations,
    compute_state_vs_fges_correlations,
    evaluate_fges_correlation_conservation,
    evaluate_subtype_response_associations,
    evaluate_within_subtype_deconv_predictions,
    transform_compositional_fractions,
    wilson_ci,
)


def _synthesize_bagaev_quadrant_df(n_per_subtype: int = 15, seed: int = 42) -> pl.DataFrame:
    """Synthesize 4 canonical Bagaev TME subtypes (IE, IE/F, F, D) across all 29 Fges."""
    rng = np.random.default_rng(seed)
    imm_sigs = set(BAGAEV_MFP_FUNCTIONAL_GROUPS["Antitumor_Immune"])
    fib_sigs = set(BAGAEV_MFP_FUNCTIONAL_GROUPS["Angiogenesis_Fibroblasts"])
    pro_sigs = set(BAGAEV_MFP_FUNCTIONAL_GROUPS["Protumor_Immune"])

    archetypes = (
        ("IE", 2.0, -1.8, 0.3, -0.5, 0.70),
        ("IE/F", 1.8, 1.9, 1.0, -0.6, 0.45),
        ("F", -1.7, 2.1, 0.6, 0.2, 0.12),
        ("D", -1.9, -1.7, -1.2, 1.8, 0.08),
    )

    rows: list[dict[str, str | int | float]] = []
    idx = 0
    for subtype, mu_imm, mu_fib, mu_pro, mu_mal, p_resp in archetypes:
        for k in range(n_per_subtype):
            row: dict[str, str | int | float] = {
                "sample_id": f"S_{idx:03d}",
                "MFP": subtype,
                "cohort": "C1" if k % 2 == 0 else "C2",
                "response_binary": int(rng.random() < p_resp),
            }
            for sig in BAGAEV_MFP_SIGNATURE_ORDER:
                center = (
                    mu_imm
                    if sig in imm_sigs
                    else mu_fib
                    if sig in fib_sigs
                    else mu_pro
                    if sig in pro_sigs
                    else mu_mal
                )
                row[sig] = float(center + rng.normal(0.0, 0.30))
            rows.append(row)
            idx += 1

    return pl.DataFrame(rows)


def test_classify_bagaev_knn_recovers_four_subtypes() -> None:
    """Verify K=35 nearest-neighbor classification accurately recovers IE, IE/F, F, and D."""
    ref_df = _synthesize_bagaev_quadrant_df(n_per_subtype=25, seed=10)
    query_df = _synthesize_bagaev_quadrant_df(n_per_subtype=12, seed=99)

    res = classify_bagaev_knn(query_df, ref_df, config=BagaevKnnConfig(k_neighbors=35))
    assert isinstance(res, Success)
    pred_df = res.unwrap()
    assert pred_df.height == query_df.height

    metrics = compare_partition_concordance(
        query_df["MFP"].to_list(),
        pred_df["bagaev_subtype"].to_list(),
    )
    assert metrics["accuracy"] >= 0.95
    assert metrics["cohen_kappa"] >= 0.92


def test_cluster_bagaev_louvain_discovers_and_labels_four_subtypes() -> None:
    """Verify unsupervised correlation-graph Louvain clustering recovers the 4 biological subtypes."""
    df = _synthesize_bagaev_quadrant_df(n_per_subtype=20, seed=7)
    res = cluster_bagaev_louvain(df, config=BagaevLouvainConfig(seed=7))
    assert isinstance(res, Success)
    out = res.unwrap()
    assert out.assignments_df.height == df.height
    assert set(out.assignments_df["louvain_subtype"].to_list()) == {"IE", "IE/F", "F", "D"}

    metrics = compare_partition_concordance(
        df["MFP"].to_list(),
        out.assignments_df["louvain_subtype"].to_list(),
    )
    assert metrics["ari"] >= 0.85
    assert metrics["accuracy"] >= 0.90


def test_evaluate_subtype_response_and_deconvolution_integration() -> None:
    """Verify ORR, CMH odds ratios, state enrichment, and within-subtype AUC evaluation."""
    df = _synthesize_bagaev_quadrant_df(n_per_subtype=30, seed=21)
    rng = np.random.default_rng(21)
    y = df["response_binary"].to_numpy().astype(np.float64)
    st_arr = df["MFP"].to_numpy()
    enriched = df.with_columns(
        [
            pl.col("MFP").alias("bagaev_subtype"),
            pl.Series("M01", np.where(st_arr == "IE", 0.35, 0.05) + rng.uniform(0.0, 0.03, df.height)),
            pl.Series("M02", np.where(st_arr == "F", 0.40, 0.04) + rng.uniform(0.0, 0.03, df.height)),
            pl.Series("pred_pct", 25.0 + 40.0 * y + rng.normal(0.0, 10.0, df.height)),
        ]
    )

    resp_res = evaluate_subtype_response_associations(
        enriched, subtype_col="bagaev_subtype", cohort_col=Some("cohort")
    )
    assert isinstance(resp_res, Success)
    resp_df = resp_res.unwrap()
    assert resp_df.height == 4
    ie_orr = float(resp_df.filter(pl.col("subtype") == "IE")["response_rate"][0])
    d_orr = float(resp_df.filter(pl.col("subtype") == "D")["response_rate"][0])
    assert ie_orr > d_orr
    assert float(resp_df["cmh_favorable_or"][0]) > 1.5

    enr_res = analyze_subtype_deconvolution_enrichment(enriched, ("M01", "M02"))
    assert isinstance(enr_res, Success)
    enr_df = enr_res.unwrap()
    assert enr_df.filter(pl.col("cell_state") == "M01")["dominant_subtype"][0] == "IE"
    assert enr_df.filter(pl.col("cell_state") == "M02")["dominant_subtype"][0] == "F"

    within_res = evaluate_within_subtype_deconv_predictions(enriched)
    assert isinstance(within_res, Success)
    within_df = within_res.unwrap()
    assert float(within_df.filter(pl.col("stratum") == "IE")["roc_auc"][0]) > 0.75

    cons_res = evaluate_fges_correlation_conservation(enriched, enriched)
    assert isinstance(cons_res, Success)
    assert np.isclose(cons_res.unwrap()["offdiag_pearson_r_vs_tcga"], 1.0)


@given(
    k=st.integers(min_value=0, max_value=100),
    n=st.integers(min_value=1, max_value=100),
)
@settings(max_examples=40)
def test_wilson_ci_bounds_property(k: int, n: int) -> None:
    """Property: Wilson 95% CI always brackets k/n inside [0, 1]."""
    k_clamped = min(k, n)
    lo, hi = wilson_ci(k_clamped, n)
    p = float(k_clamped) / float(n)
    assert 0.0 <= lo <= p + 1e-12
    assert p - 1e-12 <= hi <= 1.0


@given(
    scale=st.floats(min_value=0.05, max_value=50.0, allow_nan=False, allow_infinity=False),
)
@settings(max_examples=30)
def test_transform_compositional_fractions_clr_scale_invariance_property(scale: float) -> None:
    """Property: CLR rows sum to zero and are invariant to positive multiplicative scaling."""
    cols = ("S1", "S2", "S3", "S4")
    base_vals = np.array([[0.10, 0.25, 0.45, 0.20], [0.55, 0.15, 0.10, 0.20]])
    df1 = pl.DataFrame({"sample_id": ["A", "B"], **{c: base_vals[:, i].tolist() for i, c in enumerate(cols)}})
    df2 = pl.DataFrame(
        {"sample_id": ["A", "B"], **{c: (base_vals[:, i] * scale).tolist() for i, c in enumerate(cols)}}
    )

    r1 = transform_compositional_fractions(df1, cols, method="clr", pseudocount=0.0)
    r2 = transform_compositional_fractions(df2, cols, method="clr", pseudocount=0.0)
    assert isinstance(r1, Success) and isinstance(r2, Success)

    m1 = r1.unwrap().select(list(cols)).to_numpy()
    m2 = r2.unwrap().select(list(cols)).to_numpy()
    assert np.allclose(np.sum(m1, axis=1), 0.0, atol=1e-10)
    assert np.allclose(m1, m2, atol=1e-9)


def test_compute_deconv_subtype_correlations_and_fges() -> None:
    """Verify state-vs-subtype Spearman rho, within-cohort Fisher-z meta rho, and Fges links."""
    df = _synthesize_bagaev_quadrant_df(n_per_subtype=25, seed=77)
    rng = np.random.default_rng(77)
    st_arr = df["MFP"].to_numpy()
    # Build closed simplex fractions aligned strongly with subtypes
    s_ie = np.where(st_arr == "IE", 0.55, 0.08) + rng.uniform(0.0, 0.04, df.height)
    s_ief = np.where(st_arr == "IE/F", 0.55, 0.08) + rng.uniform(0.0, 0.04, df.height)
    s_f = np.where(st_arr == "F", 0.55, 0.08) + rng.uniform(0.0, 0.04, df.height)
    s_d = np.where(st_arr == "D", 0.55, 0.08) + rng.uniform(0.0, 0.04, df.height)
    raw = np.column_stack([s_ie, s_ief, s_f, s_d])
    raw = raw / np.sum(raw, axis=1, keepdims=True)

    joined = df.with_columns(
        [
            pl.col("MFP").alias("bagaev_subtype"),
            pl.Series("M_IE", raw[:, 0]),
            pl.Series("M_IEF", raw[:, 1]),
            pl.Series("M_F", raw[:, 2]),
            pl.Series("M_D", raw[:, 3]),
        ]
    )
    corr_res = compute_deconv_subtype_correlations(joined, ("M_IE", "M_IEF", "M_F", "M_D"))
    assert isinstance(corr_res, Success)
    c_df = corr_res.unwrap()
    assert c_df.height == 16

    row_ie = c_df.filter((pl.col("cell_state") == "M_IE") & (pl.col("subtype") == "IE"))
    assert float(row_ie["spearman_rho_raw"][0]) > 0.60
    assert float(row_ie["within_cohort_meta_rho"][0]) > 0.60
    assert float(row_ie["fdr_q_raw"][0]) < 1e-4

    fges_res = compute_state_vs_fges_correlations(joined, ("M_IE", "M_F"))
    assert isinstance(fges_res, Success)
    f_df = fges_res.unwrap()
    assert f_df.height == 2 * len(BAGAEV_MFP_SIGNATURE_ORDER)
    caf_link = f_df.filter((pl.col("cell_state") == "M_F") & (pl.col("fges_signature") == "CAF"))
    ie_caf_link = f_df.filter((pl.col("cell_state") == "M_IE") & (pl.col("fges_signature") == "CAF"))
    assert float(caf_link["spearman_rho"][0]) > 0.30
    assert float(caf_link["fdr_q"][0]) < 0.01
    assert float(ie_caf_link["spearman_rho"][0]) < -0.25


def test_cluster_deconvoluted_samples_ward_and_louvain_and_cmh_k_by_2() -> None:
    """Verify planted 4-cluster simplex recovery (Ward & CLR-Louvain) and Kx2 CMH stratification."""
    rng = np.random.default_rng(123)
    n_per = 25
    state_cols = ("S1", "S2", "S3", "S4", "S5", "S6")
    # 4 planted cellular compositions
    proto = np.array(
        [
            [0.45, 0.35, 0.05, 0.05, 0.05, 0.05],
            [0.25, 0.05, 0.45, 0.15, 0.05, 0.05],
            [0.05, 0.05, 0.15, 0.45, 0.25, 0.05],
            [0.05, 0.05, 0.05, 0.05, 0.35, 0.45],
        ]
    )
    rows: list[dict[str, str | int | float]] = []
    idx = 0
    for c_idx in range(4):
        p_true = (0.75, 0.55, 0.25, 0.15)[c_idx]
        for k in range(n_per):
            coh = "TrialA" if k % 2 == 0 else "TrialB"
            vec = proto[c_idx] + rng.uniform(0.0, 0.04, len(state_cols))
            vec = vec / np.sum(vec)
            row: dict[str, str | int | float] = {
                "sample_id": f"D_{idx:03d}",
                "true_cluster": f"T{c_idx + 1}",
                "cohort": coh,
                "response_binary": int(rng.random() < p_true),
            }
            for s_i, sc in enumerate(state_cols):
                row[sc] = float(vec[s_i])
            rows.append(row)
            idx += 1

    df = pl.DataFrame(rows)
    w_res = cluster_deconvoluted_samples(
        df, state_cols, config=DeconvClusterConfig(k_values=(2, 3, 4, 5), algorithm="aitchison_ward")
    )
    assert isinstance(w_res, Success)
    w_out = w_res.unwrap()
    assert w_out.assignments_df.height == df.height
    assert compare_partition_concordance(df["true_cluster"].to_list(), w_out.assignments_df["cluster_k4"].to_list())["ari"] >= 0.90

    l_res = cluster_deconvoluted_samples(
        df, state_cols, config=DeconvClusterConfig(k_values=(3, 4), algorithm="clr_louvain")
    )
    assert isinstance(l_res, Success)
    assert compare_partition_concordance(df["true_cluster"].to_list(), l_res.unwrap().assignments_df["cluster_k4"].to_list())["ari"] >= 0.85

    joined_eval = df.join(w_out.assignments_df.select(["sample_id", "cluster_k4"]), on="sample_id")
    q_sig, p_sig = compute_cmh_k_by_2_p_value(joined_eval, cluster_col="cluster_k4")
    assert q_sig > 10.0
    assert p_sig < 0.01

    # Pure cohort confounding (Simpson's paradox): within each cohort cluster has identical ORR,
    # while cluster prevalence & baseline ORR differ strongly between TrialA (60%) and TrialB (10%).
    conf_rows: list[dict[str, str | int]] = []
    # TrialA: 40 A (24 R, 16 NR = 60%), 10 B (6 R, 4 NR = 60%)
    for idx_c in range(40):
        conf_rows.append({"sample_id": f"TA_A_{idx_c}", "cohort": "TrialA", "cl": "A", "response_binary": 1 if idx_c < 24 else 0})
    for idx_c in range(10):
        conf_rows.append({"sample_id": f"TA_B_{idx_c}", "cohort": "TrialA", "cl": "B", "response_binary": 1 if idx_c < 6 else 0})
    # TrialB: 10 A (1 R, 9 NR = 10%), 40 B (4 R, 36 NR = 10%)
    for idx_c in range(10):
        conf_rows.append({"sample_id": f"TB_A_{idx_c}", "cohort": "TrialB", "cl": "A", "response_binary": 1 if idx_c < 1 else 0})
    for idx_c in range(40):
        conf_rows.append({"sample_id": f"TB_B_{idx_c}", "cohort": "TrialB", "cl": "B", "response_binary": 1 if idx_c < 4 else 0})

    conf_df = pl.DataFrame(conf_rows)
    q_conf, p_conf = compute_cmh_k_by_2_p_value(conf_df, cluster_col="cl")
    assert q_conf < 1e-12
    assert p_conf > 0.90

    # Exact closed-form 2x2 hypergeometric check (N=80: A=25/40 R vs B=10/40 R)
    # Q_CMH = (N - 1) * (25*30 - 15*10)^2 / (40 * 40 * 35 * 45)
    exact_rows = (
        [{"sample_id": f"EA_{i}", "cohort": "One", "cl": "A", "response_binary": 1 if i < 25 else 0} for i in range(40)]
        + [{"sample_id": f"EB_{i}", "cohort": "One", "cl": "B", "response_binary": 1 if i < 10 else 0} for i in range(40)]
    )
    q_exact, _ = compute_cmh_k_by_2_p_value(pl.DataFrame(exact_rows), cluster_col="cl")
    expected_q = (79.0 * (25.0 * 30.0 - 15.0 * 10.0) ** 2) / (40.0 * 40.0 * 35.0 * 45.0)
    assert abs(q_exact - expected_q) < 1e-9


def test_knn_tie_break_coarsening_bijection_and_k2_validity() -> None:
    """Verify unbiased k-NN tie breaking, 2-phase Louvain coarsening, bijection & K=2/K>=3 selection."""
    # 1. Exact 2-vs-2 vote tie between IE (distant) and F (near-identical) must select F
    sigs = BAGAEV_MFP_SIGNATURE_ORDER
    query = pl.DataFrame(
        {"sample_id": ["Q1"], **{s: [1.0 if s in BAGAEV_MFP_FUNCTIONAL_GROUPS["Angiogenesis_Fibroblasts"] else -1.0] for s in sigs}}
    )
    ref = pl.DataFrame(
        [
            {"sample_id": "R1", "MFP": "F", **{s: (1.0 if s in BAGAEV_MFP_FUNCTIONAL_GROUPS["Angiogenesis_Fibroblasts"] else -1.0) for s in sigs}},
            {"sample_id": "R2", "MFP": "F", **{s: (0.95 if s in BAGAEV_MFP_FUNCTIONAL_GROUPS["Angiogenesis_Fibroblasts"] else -0.95) for s in sigs}},
            {"sample_id": "R3", "MFP": "IE", **{s: (-1.0 if s in BAGAEV_MFP_FUNCTIONAL_GROUPS["Antitumor_Immune"] else -0.9) for s in sigs}},
            {"sample_id": "R4", "MFP": "IE", **{s: (-1.0 if s in BAGAEV_MFP_FUNCTIONAL_GROUPS["Antitumor_Immune"] else -0.8) for s in sigs}},
        ]
    )
    tie_res = classify_bagaev_knn(query, ref, config=BagaevKnnConfig(k_neighbors=4, metric="correlation"))
    assert isinstance(tie_res, Success)
    assert tie_res.unwrap()["bagaev_subtype"][0] == "F"

    # 2. Bijective subtype mapping when M=3 clusters (must produce 3 distinct subtypes)
    mat3 = np.vstack(
        [
            np.array([1.5 if s in BAGAEV_MFP_FUNCTIONAL_GROUPS["Antitumor_Immune"] else -1.0 for s in sigs]),
            np.array([1.6 if s in BAGAEV_MFP_FUNCTIONAL_GROUPS["Antitumor_Immune"] else 0.9 if s in BAGAEV_MFP_FUNCTIONAL_GROUPS["Angiogenesis_Fibroblasts"] else -1.0 for s in sigs]),
            np.array([-1.5 for _ in sigs]),
        ]
    )
    lbls3, mapping3 = assign_clusters_to_bagaev_subtypes(np.array([0, 1, 2]), mat3, sigs)
    assert len(set(lbls3)) == 3
    assert len(set(mapping3.values())) == 3

    # 3. 2-phase Louvain recovers 2 strong cliques cleanly
    adj = np.block([[np.ones((6, 6)), 0.05 * np.ones((6, 6))], [0.05 * np.ones((6, 6)), np.ones((6, 6))]])
    comm = _louvain_modularity_partition(adj, seed=42)
    assert len(np.unique(comm)) == 2
    assert len(set(comm[:6])) == 1 and len(set(comm[6:])) == 1 and comm[0] != comm[6]

    # 4. Bimodal simplex gives optimal_k == 2 vs optimal_k_multi == 3 & clr_louvain respects exact k_eff
    rng = np.random.default_rng(7)
    st_cols = ("S1", "S2", "S3", "S4")
    bi_rows = []
    for i in range(40):
        vec = ([0.45, 0.45, 0.05, 0.05] if i < 20 else [0.05, 0.05, 0.45, 0.45]) + rng.uniform(0.0, 0.01, 4)
        vec = vec / np.sum(vec)
        bi_rows.append({"sample_id": f"B_{i}", **{st_cols[j]: float(vec[j]) for j in range(4)}})
    bi_res = cluster_deconvoluted_samples(
        pl.DataFrame(bi_rows), st_cols, config=DeconvClusterConfig(k_values=(2, 3, 4), algorithm="clr_louvain")
    )
    assert isinstance(bi_res, Success)
    bi_out = bi_res.unwrap()
    assert bi_out.optimal_k == 2
    assert bi_out.optimal_k_multi >= 3
    assert bi_out.assignments_df["cluster_k3"].n_unique() == 3
    assert bi_out.assignments_df["cluster_k4"].n_unique() == 4


def test_benchmark_bagaev_vs_deconv_cv_patient_group_isolation() -> None:
    """Verify 5x5 CV groups repeated biopsies by patient_id and outputs OOF scores."""
    df = _synthesize_bagaev_quadrant_df(n_per_subtype=20, seed=55)
    # Assign shared cohort & patient IDs in pairs so repeated biopsies belong to the same patient
    pids = [f"P_{i // 2:02d}" for i in range(df.height)]
    cohs = ["C1" if (i // 2) % 2 == 0 else "C2" for i in range(df.height)]
    mat_df = df.with_columns(
        [
            pl.col("MFP").alias("bagaev_subtype"),
            pl.Series("patient_id", pids),
            pl.Series("cohort", cohs),
            pl.when(pl.col("MFP").is_in(["IE", "IE/F"])).then(0.45).otherwise(0.15).alias("M01"),
            pl.when(pl.col("MFP") == "IE").then(0.35).otherwise(0.10).alias("M02"),
            pl.when(pl.col("MFP") == "F").then(0.45).otherwise(0.15).alias("M03"),
            pl.when(pl.col("MFP") == "D").then(0.45).otherwise(0.15).alias("M04"),
        ]
    )
    cv_res = benchmark_bagaev_vs_deconv_cv(
        mat_df,
        state_cols=("M01", "M02", "M03", "M04"),
        group_col=Some("patient_id"),
        n_splits=4,
        n_repeats=2,
    )
    assert isinstance(cv_res, Success)
    pan_df, coh_df, oof_df = cv_res.unwrap()
    assert pan_df.height == 4
    assert oof_df.height == mat_df.height
    assert "oof_raw_meta42_z" in oof_df.columns
    assert int(pan_df["n_unique_patients"][0]) < mat_df.height
    assert coh_df.height >= 4
