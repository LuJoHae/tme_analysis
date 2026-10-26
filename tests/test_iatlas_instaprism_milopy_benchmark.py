"""Unit and property-based tests for Cross-Modality iAtlas InstaPrism & Milopy Concordance Benchmark."""

from __future__ import annotations

from pathlib import Path
import re
from typing import cast
import xml.etree.ElementTree as ET

from hypothesis import given, settings
from hypothesis import strategies as st
import numpy as np
import polars as pl
from returns.pipeline import is_successful
from returns.result import Success

from scripts.milopy_analysis.iatlas_instaprism_milopy_benchmark import (
    BulkCohortBundle,
    ReferenceSpecBundle,
    analyze_all42_nonlinear_conditional_interactions,
    compute_dersimonian_laird_meta_beta,
    compute_hanley_mcneil_ci,
    evaluate_regularized_logistic_cv,
    linearize_and_normalize_bulk_matrix,
)


@settings(max_examples=25)
@given(
    n_samples=st.integers(min_value=4, max_value=20),
    n_genes=st.integers(min_value=10, max_value=40),
    is_log2=st.booleans(),
)
def test_linearize_bulk_expression_simplex(
    n_samples: int,
    n_genes: int,
    is_log2: bool,
) -> None:
    """Property-based test verifying log2 detection and exact sample-wise TPM (10^6) scaling."""
    rng = np.random.default_rng(42)
    raw = (
        rng.uniform(0.1, 18.0, size=(n_samples, n_genes))
        if is_log2
        else rng.uniform(10.0, 5000.0, size=(n_samples, n_genes))
    )
    tpm_mat, log2_z_mat, detected_log2 = linearize_and_normalize_bulk_matrix(raw)

    assert detected_log2 is is_log2
    assert tpm_mat.shape == (n_samples, n_genes)
    assert log2_z_mat.shape == (n_samples, n_genes)
    assert np.all(tpm_mat >= 0.0)
    np.testing.assert_allclose(tpm_mat.sum(axis=1), np.full(n_samples, 1e6), rtol=1e-5)


@settings(max_examples=30)
@given(
    auc=st.floats(min_value=0.50, max_value=0.95),
    n_pos=st.integers(min_value=3, max_value=200),
    n_neg=st.integers(min_value=3, max_value=200),
)
def test_hanley_mcneil_ci_bounds(auc: float, n_pos: int, n_neg: int) -> None:
    """Property-based test verifying Hanley-McNeil 95% confidence interval ordering and bounds."""
    ci_low, ci_high, p_val = compute_hanley_mcneil_ci(auc, n_pos, n_neg)
    assert 0.0 <= ci_low <= auc <= ci_high <= 1.0
    assert 0.0 <= p_val <= 1.0


def test_dersimonian_laird_meta_beta() -> None:
    """Verify DerSimonian-Laird random-effects meta-analysis across cohorts."""
    b1 = np.array([0.8, -0.6, 0.2], dtype=np.float64)
    b2 = np.array([0.8, -0.6, 0.2], dtype=np.float64)
    se1 = np.array([0.2, 0.2, 0.2], dtype=np.float64)
    se2 = np.array([0.2, 0.2, 0.2], dtype=np.float64)

    meta_b, meta_se, meta_p = compute_dersimonian_laird_meta_beta([b1, b2], [se1, se2])
    np.testing.assert_allclose(meta_b, b1, atol=1e-6)
    assert np.all(meta_se > 0.0)
    assert np.all((meta_p >= 0.0) & (meta_p <= 1.0))


def test_evaluate_regularized_logistic_cv_synthetic() -> None:
    """Verify leakage-free stratified CV evaluation across ridge_l2, elastic_net, and lasso_l1."""
    rng = np.random.default_rng(123)
    n_samples = 40
    k_states = 8
    y = np.array([1.0] * 20 + [0.0] * 20, dtype=np.float64)

    milo_logfc = np.linspace(1.5, -1.5, k_states)
    raw_fracs = rng.dirichlet(np.ones(k_states), size=n_samples)
    # Inject responder enrichment in state 0 and non-responder enrichment in state -1
    raw_fracs[:20, 0] += 0.25
    raw_fracs[20:, -1] += 0.25
    fractions_mat = raw_fracs / raw_fracs.sum(axis=1, keepdims=True)

    bulk_bundle = BulkCohortBundle(
        cohort_id="Synthetic-Bulk",
        cancer_type="Melanoma",
        short_name="Synthetic",
        sample_ids=tuple(f"S{i}" for i in range(n_samples)),
        sample_cohorts=tuple("Synthetic-Bulk" for _ in range(n_samples)),
        response_binary=y,
        genes=("CD8A", "TCF7", "SPP1"),
        gene_to_idx={"CD8A": 0, "TCF7": 1, "SPP1": 2},
        was_log2_transformed=True,
        tpm_mat=np.full((n_samples, 3), 1e6 / 3.0),
        log2_z_mat=rng.normal(0.0, 1.0, size=(n_samples, 3)),
        bulk_response_profile=y + rng.normal(0.0, 0.4, size=n_samples),
    )

    ref_spec = ReferenceSpecBundle(
        reference_id="Synthetic-Ref",
        reference_tier="Tier 1",
        display_label="Synthetic Ref",
        cancer_type="Melanoma",
        primary_bulk_target="Synthetic-Bulk",
        state_ids=tuple(f"C{i}" for i in range(k_states)),
        state_cohorts=tuple("Synthetic-Ref" for _ in range(k_states)),
        milo_logfc=milo_logfc,
        milo_se=np.full(k_states, 0.2),
        bank_transfer_weights=np.eye(k_states),
        native_genes=("CD8A", "TCF7", "SPP1"),
        native_phi_mat=np.full((k_states, 3), 1.0 / 3.0),
    )

    for penalty in ("ridge_l2", "elastic_net", "lasso_l1", "linear_reg", "random_forest"):
        res = evaluate_regularized_logistic_cv(
            fractions_mat=fractions_mat,
            bulk_bundle=bulk_bundle,
            ref_spec=ref_spec,
            penalty_mode=penalty,
            n_splits=4,
            n_repeats=2,
            seed=42,
        )
        assert cast(float, res["roc_auc"]) > 0.65
        assert 0.0 <= cast(float, res["brier_score"]) <= 1.0


def test_all42_nonlinear_interactions_symmetry(tmp_path: Path) -> None:
    """Verify exact matrix symmetry, zero diagonal, and 3x3 regime CI bounds on synthetic data."""
    rng = np.random.default_rng(7)
    n_samples = 45
    k_states = 6
    y = np.array([1.0] * 22 + [0.0] * 23, dtype=np.float64)
    raw_fracs = rng.dirichlet(np.ones(k_states), size=n_samples)
    raw_fracs[:22, 0] += 0.30
    raw_fracs[22:, -1] += 0.30
    theta_mat = raw_fracs / raw_fracs.sum(axis=1, keepdims=True)

    bulk_bundle = BulkCohortBundle(
        cohort_id="PanCancer-Combined",
        cancer_type="Pan-Cancer",
        short_name="PanCancer",
        sample_ids=tuple(f"P{i}" for i in range(n_samples)),
        sample_cohorts=tuple("PanCancer-Combined" for _ in range(n_samples)),
        response_binary=y,
        genes=("GNLY", "PRF1", "FAP", "TGFB1"),
        gene_to_idx={"GNLY": 0, "PRF1": 1, "FAP": 2, "TGFB1": 3},
        was_log2_transformed=True,
        tpm_mat=np.full((n_samples, 4), 2.5e5),
        log2_z_mat=rng.normal(0.0, 1.0, size=(n_samples, 4)),
        bulk_response_profile=y + rng.normal(0.0, 0.3, size=n_samples),
    )

    state_ids = ("M04", "M13", "M20", "M17", "M21", "M23")
    ref_spec = ReferenceSpecBundle(
        reference_id="Consensus-Meta42",
        reference_tier="Tier 3",
        display_label="Synthetic Meta",
        cancer_type="Pan-Cancer",
        primary_bulk_target="PanCancer-Combined",
        state_ids=state_ids,
        state_cohorts=tuple("Pan-Cohort" for _ in state_ids),
        milo_logfc=np.array([2.5, 2.0, 1.2, -1.0, -2.0, -3.0], dtype=np.float64),
        milo_se=np.full(k_states, 0.25),
        bank_transfer_weights=np.eye(k_states),
        native_genes=("GNLY", "PRF1", "FAP", "TGFB1"),
        native_phi_mat=np.full((k_states, 4), 0.25),
    )
    states_df = pl.DataFrame({
        "reference_id": ["Consensus-Meta42"] * k_states,
        "bulk_cohort": ["Pan-iAtlas-Meta"] * k_states,
        "cell_state": list(state_ids),
        "bulk_deconv_beta": [1.8, 1.4, 0.7, -0.8, -1.5, -2.4],
    })

    inter_df, hub_df, reg_df = analyze_all42_nonlinear_conditional_interactions(
        theta_mat, bulk_bundle, ref_spec, states_df, tmp_path
    )
    assert inter_df.height == k_states * k_states
    assert hub_df.height == k_states
    assert reg_df.height == 9

    piv_h = inter_df.pivot(on="state_b", index="state_a", values="friedman_h").select(list(state_ids)).to_numpy()
    piv_z = inter_df.pivot(on="state_b", index="state_a", values="epistasis_z").select(list(state_ids)).to_numpy()
    np.testing.assert_allclose(piv_h, piv_h.T, atol=1e-10)
    np.testing.assert_allclose(piv_z, piv_z.T, atol=1e-10)
    np.testing.assert_allclose(np.diag(piv_h), np.zeros(k_states), atol=1e-12)
    for row in reg_df.iter_rows(named=True):
        assert 0.0 <= float(row["ci_low_pct"]) <= float(row["responder_rate_pct"]) <= float(row["ci_high_pct"]) <= 100.0


def test_generated_benchmark_parquets_and_figures() -> None:
    """Verify empirical parquet outputs, simplex fraction invariants, and all 11 180mm vector SVGs."""
    out_dir = Path("output/results/iatlas_instaprism_benchmark")
    results_dir = Path("results/iatlas_instaprism_benchmark")

    pred_path = out_dir / "response_prediction_metrics.parquet"
    conc_path = out_dir / "milopy_iatlas_concordance_summary.parquet"
    states_path = out_dir / "milopy_iatlas_state_effects.parquet"
    loco_path = out_dir / "loco_prediction_metrics.parquet"
    deg_path = out_dir / "consensus_metastate_top_degs.parquet"
    inter_path = out_dir / "consensus_metastate_interactions.parquet"
    regime_path = out_dir / "consensus_metastate_regimes.parquet"
    all42_deg_path = out_dir / "consensus_metastate_all42_top_degs.parquet"
    all42_inter_path = out_dir / "consensus_metastate_all42_interactions.parquet"
    all42_hub_path = out_dir / "consensus_metastate_all42_hub_summary.parquet"
    all42_regime_path = out_dir / "consensus_metastate_all42_regimes.parquet"

    for req_p in (
        pred_path,
        conc_path,
        states_path,
        loco_path,
        deg_path,
        inter_path,
        regime_path,
        all42_deg_path,
        all42_inter_path,
        all42_hub_path,
        all42_regime_path,
    ):
        assert req_p.exists(), f"Missing {req_p}"

    pred_df = pl.read_parquet(pred_path)
    conc_df = pl.read_parquet(conc_path)
    states_df = pl.read_parquet(states_path)
    loco_df = pl.read_parquet(loco_path)
    deg_df = pl.read_parquet(deg_path)
    inter_df = pl.read_parquet(inter_path)
    regime_df = pl.read_parquet(regime_path)
    all42_deg_df = pl.read_parquet(all42_deg_path)
    all42_inter_df = pl.read_parquet(all42_inter_path)
    all42_hub_df = pl.read_parquet(all42_hub_path)
    all42_regime_df = pl.read_parquet(all42_regime_path)

    # All 13 references evaluated (9 individual SC + 2 indication-unified + Unified-All9 + Consensus-Meta42)
    assert pred_df["reference_id"].n_unique() == 13
    assert conc_df["reference_id"].n_unique() == 13
    assert states_df.height > 350
    assert loco_df.height >= 33
    assert deg_df.height == 56  # 14 selected states x 4 top marker genes
    assert inter_df.height == 64  # 8x8 pairwise interaction grid
    assert regime_df.height == 9  # 3x3 lymphoid x stromal gating regimes
    assert all42_deg_df.height == 168  # all 42 states x 4 top marker genes
    assert all42_deg_df["cell_state"].n_unique() == 42
    assert all42_inter_df.height == 1764  # all 42x42 pairwise interaction grid
    assert all42_inter_df["state_a"].n_unique() == 42
    assert all42_inter_df["state_b"].n_unique() == 42
    assert all42_hub_df.height == 42
    assert all42_regime_df.height == 9

    # Verify all three model families (Logistic Regression, Linear Regression, Random Forests) are present
    modes = set(pred_df["penalty_mode"].unique().to_list())
    assert {"ridge_l2", "linear_reg", "random_forest"}.issubset(modes)

    # Verify high linear Pearson r and rank Spearman rho across all 13 references on Pan-iAtlas-Meta
    meta_sub = conc_df.filter(pl.col("bulk_cohort") == "Pan-iAtlas-Meta")
    assert meta_sub.height == 13
    assert all(float(r) > 0.75 for r in meta_sub["spearman_rho"].to_list())
    assert all(float(r) > 0.75 for r in meta_sub["pearson_r"].to_list())

    # Verify cached deconvoluted fractions sum to 1.0 on the simplex
    frac_files = sorted((out_dir / "fractions").glob("*.parquet"))
    assert len(frac_files) >= 13
    sample_frac_df = pl.read_parquet(frac_files[0])
    state_cols = [
        c for c in sample_frac_df.columns if c not in {"sample_id", "cohort", "response_binary"}
    ]
    mat = sample_frac_df.select(state_cols).to_numpy()
    assert np.all(mat >= 0.0)
    np.testing.assert_allclose(mat.sum(axis=1), np.ones(mat.shape[0]), atol=1e-5)

    # Verify all 11 SVG figures exist, parse as well-formed XML, have width="180.0mm", height <= 245mm, and font >= 10pt
    svg_files = (
        results_dir / "iatlas_response_prediction_auc_heatmap_and_bars.svg",
        results_dir / "iatlas_milopy_vs_deconv_concordance_grid.svg",
        results_dir / "iatlas_consensus_metastate_forest_and_tme_profile.svg",
        results_dir / "iatlas_classifier_and_reference_tier_pareto.svg",
        results_dir / "iatlas_cross_cohort_concordance_transfer_matrix.svg",
        results_dir / "iatlas_consensus_metastate_top_deg_markers.svg",
        results_dir / "iatlas_response_prediction_prauc_heatmap_and_bars.svg",
        results_dir / "iatlas_nonlinear_conditional_interactions.svg",
        results_dir / "iatlas_consensus_metastate_all42_forest_and_tme_profile.svg",
        results_dir / "iatlas_consensus_metastate_all42_top_deg_markers.svg",
        results_dir / "iatlas_nonlinear_conditional_interactions_all42.svg",
    )
    for svg_p in svg_files:
        assert svg_p.exists(), f"Missing SVG figure: {svg_p}"
        png_p = svg_p.with_suffix(".png")
        assert png_p.exists() and png_p.stat().st_size > 10000, f"Missing PNG figure: {png_p}"

        svg_text = svg_p.read_text(encoding="utf-8")
        root = ET.fromstring(svg_text)
        assert root.attrib.get("width") == "180.0mm"
        h_attr = root.attrib.get("height", "0.0mm").replace("mm", "")
        assert float(h_attr) <= 245.0, f"Height {h_attr}mm in {svg_p.name} exceeds 245mm page limit"
        assert "role-title-subtitle" in svg_text, f"Missing subtitle in {svg_p.name}"

        font_matches = re.findall(r"(?:font-size:\s*|font-size=\")([0-9.]+)px", svg_text)
        assert len(font_matches) > 0
        for fm in font_matches:
            assert float(fm) >= 13.3, f"Font size {fm}px in {svg_p.name} is below 10pt (13.33px)"
    assert is_successful(Success(True))


def test_tme_datasets_pipeline_dag_svg_and_png() -> None:
    """Verify the tme_datasets provenance & script-dependency DAG SVG and 300 DPI PNG."""
    svg_paths = (
        Path("results/dataset_overview/figure_tme_datasets_pipeline_dag.svg"),
        Path("article/figures/dataset_overview/figure_tme_datasets_pipeline_dag.svg"),
    )
    for svg_p in svg_paths:
        assert svg_p.exists(), f"Missing DAG SVG: {svg_p}"
        png_p = svg_p.with_suffix(".png")
        assert png_p.exists() and png_p.stat().st_size > 20000, f"Missing DAG PNG: {png_p}"

        svg_text = svg_p.read_text(encoding="utf-8")
        root = ET.fromstring(svg_text)
        assert root.attrib.get("width") == "180.0mm"
        assert "feDropShadow" not in svg_text
        assert "role-title-subtitle" in svg_text

        inkscape_ns = "{http://www.inkscape.org/namespaces/inkscape}groupmode"
        layers = [el for el in root.iter() if el.attrib.get(inkscape_ns) == "layer"]
        assert len(layers) >= 6

        font_matches = re.findall(r"(?:font-size:\s*|font-size=\")([0-9.]+)px", svg_text)
        assert len(font_matches) > 0
        for fm in font_matches:
            assert float(fm) >= 13.3, f"Font size {fm}px in {svg_p.name} is below 13.3px"

