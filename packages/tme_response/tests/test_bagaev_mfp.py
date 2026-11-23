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
    analyze_subtype_deconvolution_enrichment,
    classify_bagaev_knn,
    cluster_bagaev_louvain,
    compare_partition_concordance,
    evaluate_fges_correlation_conservation,
    evaluate_subtype_response_associations,
    evaluate_within_subtype_deconv_predictions,
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
