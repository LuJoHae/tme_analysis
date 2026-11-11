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

import importlib

from scipy.special import expit

_stage3c = importlib.import_module(
    "scripts.manuscript_pipeline.3c_iatlas_instaprism_milopy_benchmark"
)
BulkCohortBundle = _stage3c.BulkCohortBundle
ReferenceSpecBundle = _stage3c.ReferenceSpecBundle
_compute_pan42_calibrated_epistasis_pvalues = (
    _stage3c._compute_pan42_calibrated_epistasis_pvalues
)
analyze_all42_pure_clr_interactions = (
    _stage3c.analyze_all42_pure_clr_interactions
)
benchmark_epistasis_discovery_strategies = (
    _stage3c.benchmark_epistasis_discovery_strategies
)
benchmark_zero_leakage_gating_strategies = (
    _stage3c.benchmark_zero_leakage_gating_strategies
)
compute_dersimonian_laird_meta_beta = _stage3c.compute_dersimonian_laird_meta_beta
compute_hanley_mcneil_ci = _stage3c.compute_hanley_mcneil_ci
evaluate_regularized_logistic_cv = _stage3c.evaluate_regularized_logistic_cv
linearize_and_normalize_bulk_matrix = _stage3c.linearize_and_normalize_bulk_matrix


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

    inter_df = analyze_all42_pure_clr_interactions(
        theta_mat, bulk_bundle, ref_spec, states_df, tmp_path
    )
    assert inter_df.height == k_states * k_states

    piv_h = inter_df.pivot(on="state_b", index="state_a", values="friedman_h").select(list(state_ids)).to_numpy()
    piv_z = inter_df.pivot(on="state_b", index="state_a", values="epistasis_z").select(list(state_ids)).to_numpy()
    np.testing.assert_allclose(piv_h, piv_h.T, atol=1e-10)
    np.testing.assert_allclose(piv_z, piv_z.T, atol=1e-10)
    np.testing.assert_allclose(np.diag(piv_h), np.zeros(k_states), atol=1e-12)


def test_pan42_quad_logit_spa_calibration_vs_naive_2state() -> None:
    """Verify Pan42-Quad-Logit-SPA controls compositional CLR false positives vs. naive 2-state perm."""
    rng = np.random.default_rng(2026)
    n_samples = 420
    k_states = 14
    cohorts_arr = np.asarray(
        (["Trial1"] * 140) + (["Trial2"] * 140) + (["Trial3"] * 140),
        dtype=np.str_,
    )
    # Compartment-structured Dirichlet simplex with tumor-purity dilution (omitted third-state coskewness)
    # plus co-infiltrating sibling states (0 & 1) with quadratic saturation and ZERO pairwise epistasis
    tumor_frac = rng.beta(2.2, 2.2, size=n_samples)
    immune_raw = rng.dirichlet(np.ones(k_states - 1) * 1.2, size=n_samples)
    sibling_drive = rng.exponential(0.8, size=n_samples)
    immune_raw[:, 0] += 1.6 * sibling_drive
    immune_raw[:, 1] += 1.6 * sibling_drive
    immune_raw = immune_raw / np.sum(immune_raw, axis=1, keepdims=True)
    theta_syn = np.column_stack([tumor_frac, (1.0 - tumor_frac)[:, None] * immune_raw])
    log_theta = np.log2(theta_syn + 1e-4)
    clr_syn = log_theta - np.mean(log_theta, axis=1, keepdims=True)
    z_syn = (clr_syn - np.mean(clr_syn, axis=0, keepdims=True)) / np.maximum(
        np.std(clr_syn, axis=0, keepdims=True), 1e-6
    )

    beta_true = np.array([-1.6, 1.4, -1.3, 0.9, -0.8, 0.5, -0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    delta_quad = np.array([0.0, -0.45, -0.40, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    eta_null = -0.5 + z_syn @ beta_true + ((z_syn**2) - 1.0) @ delta_quad
    p_null = np.asarray(expit(eta_null), dtype=np.float64)
    y_null = rng.binomial(1, p_null).astype(np.float64)

    cal = _compute_pan42_calibrated_epistasis_pvalues(
        z_ordered=z_syn,
        y_obs=y_null,
        cohorts_arr=cohorts_arr,
        n_boot=80,
        n_perm_diag=80,
        seed=2026,
    )
    np.testing.assert_allclose(cal.z_spa_mat, cal.z_spa_mat.T, atol=1e-12)
    np.testing.assert_allclose(cal.p_spa_mat, cal.p_spa_mat.T, atol=1e-12)
    np.testing.assert_allclose(np.diag(cal.z_spa_mat), np.zeros(k_states), atol=1e-12)
    np.testing.assert_allclose(np.diag(cal.p_spa_mat), np.ones(k_states), atol=1e-12)

    tri_i, tri_j = np.triu_indices(k_states, k=1)
    fpr_spa = float(np.mean(cal.p_spa_mat[tri_i, tri_j] <= 0.05))
    fpr_2s = float(np.mean(cal.p_perm_label_2state_mat[tri_i, tri_j] <= 0.05))
    assert fpr_spa < 0.10
    assert fpr_spa < fpr_2s
    assert cal.lambda_gc_spa < cal.lambda_gc_2state
    assert 0.0 < cal.omnibus_boot_p <= 1.0


def test_continuous_fwer_and_fdr_no_discretization_floor() -> None:
    """Verify continuous-|Z| Westfall-Young FWER & Storey/Ge FDR detect strong epistasis (<0.05, no 0.58 floor)."""
    rng = np.random.default_rng(777)
    n_samples = 480
    k_states = 16
    cohorts_arr = np.asarray(
        (["Trial1"] * 160) + (["Trial2"] * 160) + (["Trial3"] * 160),
        dtype=np.str_,
    )
    z_syn = rng.normal(0.0, 1.0, size=(n_samples, k_states))
    z_syn = (z_syn - np.mean(z_syn, axis=0, keepdims=True)) / np.std(
        z_syn, axis=0, keepdims=True
    )
    eta = -0.3 + 0.45 * z_syn[:, 0] - 0.40 * z_syn[:, 1] - 1.35 * (z_syn[:, 0] * z_syn[:, 1])
    y_syn = rng.binomial(1, np.asarray(expit(eta), dtype=np.float64)).astype(np.float64)

    cal = _compute_pan42_calibrated_epistasis_pvalues(
        z_ordered=z_syn,
        y_obs=y_syn,
        cohorts_arr=cohorts_arr,
        n_boot=250,
        n_perm_diag=60,
        seed=777,
    )
    assert float(cal.p_spa_mat[0, 1]) < 1e-4
    assert float(cal.q_bh_spa_mat[0, 1]) < 0.01
    assert float(cal.p_fwer_boot_mat[0, 1]) < 0.05
    assert float(cal.q_emp_boot_mat[0, 1]) < 0.05
    assert cal.omnibus_boot_p < 0.05


def test_epistasis_discovery_strategies_benchmark_synthetic(tmp_path: Path) -> None:
    """Verify multi-strategy zero-leakage epistasis discovery benchmark runs cleanly on synthetic cohorts."""
    rng = np.random.default_rng(314)
    n_samples = 120
    k_states = 10
    sample_cohorts = tuple(
        (["Hugo-iAtlas"] * 40) + (["Rosenberg-iAtlas"] * 40) + (["McDermott-iAtlas"] * 40)
    )
    y = rng.binomial(1, 0.35, size=n_samples).astype(np.float64)
    raw_fracs = rng.dirichlet(np.ones(k_states), size=n_samples)
    theta_raw_mat = raw_fracs / raw_fracs.sum(axis=1, keepdims=True)
    gene_names = ("GNLY", "PRF1", "JCHAIN", "MZB1", "VIM", "ZEB1", "CXCL8", "ITGAM", "FCGR3B", "FPR2")
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    bulk_bundle = BulkCohortBundle(
        cohort_id="PanCancer-Combined",
        cancer_type="Pan-Cancer",
        short_name="PanCancer",
        sample_ids=tuple(f"E{i}" for i in range(n_samples)),
        sample_cohorts=sample_cohorts,
        response_binary=y,
        genes=gene_names,
        gene_to_idx=gene_to_idx,
        was_log2_transformed=True,
        tpm_mat=np.full((n_samples, len(gene_names)), 1e6 / len(gene_names)),
        log2_z_mat=rng.normal(0.0, 1.0, size=(n_samples, len(gene_names))),
        bulk_response_profile=y + rng.normal(0.0, 0.25, size=n_samples),
    )
    state_ids = tuple(f"M{idx:02d}" for idx in range(1, k_states + 1))
    ref_spec = ReferenceSpecBundle(
        reference_id="Consensus-Meta42",
        reference_tier="Tier 3",
        display_label="Synthetic Meta42",
        cancer_type="Pan-Cancer",
        primary_bulk_target="PanCancer-Combined",
        state_ids=state_ids,
        state_cohorts=tuple("Pan-Cohort" for _ in state_ids),
        milo_logfc=np.linspace(2.2, -2.2, k_states),
        milo_se=np.full(k_states, 0.25),
        bank_transfer_weights=np.eye(k_states),
        native_genes=gene_names,
        native_phi_mat=np.full((k_states, len(gene_names)), 1.0 / len(gene_names)),
    )
    states_df = pl.DataFrame(
        {
            "reference_id": ["Consensus-Meta42"] * k_states,
            "bulk_cohort": ["Pan-iAtlas-Meta"] * k_states,
            "cell_state": list(state_ids),
            "bulk_deconv_beta": np.linspace(0.45, -0.45, k_states).tolist(),
        }
    )

    strat_df, sig_df = benchmark_epistasis_discovery_strategies(
        theta_raw_meta42=theta_raw_mat,
        bulk_bundle=bulk_bundle,
        meta42_spec=ref_spec,
        states_df=states_df,
        out_dir=tmp_path,
        n_boot_emp=40,
        n_null_reps=1,
        n_boot_null=30,
        seed=314,
    )
    assert strat_df.height == 6
    assert sig_df.filter(pl.col("record_type") == "tertile_2x2_cell").height == 4
    assert sig_df.filter(pl.col("record_type") == "indication_forest").height == 5
    assert sig_df.filter(pl.col("record_type") == "all861_pair_volcano").height == (
        k_states * (k_states - 1) // 2
    )
    assert sig_df.filter(pl.col("record_type") == "m29_dose_curve").height == 122
    assert sig_df.filter(pl.col("record_type") == "m29_quintile_bin").height == 10
    assert (tmp_path / "epistasis_discovery_strategies_summary.parquet").exists()
    assert (tmp_path / "epistasis_significant_pairs_calibrated.parquet").exists()


def test_zero_leakage_gating_synthetic_null(tmp_path: Path) -> None:
    """Verify permutation-null calibration separates leaky y-residualization from strict zero-leakage methods."""
    rng = np.random.default_rng(19)
    n_samples = 90
    k_states = 42
    sample_cohorts = tuple(["TrialA"] * 30 + ["TrialB"] * 30 + ["TrialC"] * 30)
    y = np.array(([1.0, 0.0] * 45), dtype=np.float64)
    rng.shuffle(y)

    raw_fracs = rng.dirichlet(np.ones(k_states), size=n_samples)
    theta_raw_mat = raw_fracs / raw_fracs.sum(axis=1, keepdims=True)

    bulk_bundle = BulkCohortBundle(
        cohort_id="PanCancer-Combined",
        cancer_type="Pan-Cancer",
        short_name="PanCancer",
        sample_ids=tuple(f"Z{i}" for i in range(n_samples)),
        sample_cohorts=sample_cohorts,
        response_binary=y,
        genes=("GNLY", "PRF1", "FAP", "TGFB1"),
        gene_to_idx={"GNLY": 0, "PRF1": 1, "FAP": 2, "TGFB1": 3},
        was_log2_transformed=True,
        tpm_mat=np.full((n_samples, 4), 2.5e5),
        log2_z_mat=rng.normal(0.0, 1.0, size=(n_samples, 4)),
        bulk_response_profile=y + rng.normal(0.0, 0.25, size=n_samples),
        response_recist=tuple("CR" if v > 0.5 else "PD" for v in y),
        os_months=np.linspace(3.0, 36.0, n_samples),
        os_event=np.array(([1.0, 0.0, 1.0] * 30), dtype=np.float64),
    )

    state_ids = tuple(f"M{idx:02d}" for idx in range(1, k_states + 1))
    ref_spec = ReferenceSpecBundle(
        reference_id="Consensus-Meta42",
        reference_tier="Tier 3",
        display_label="Synthetic Meta42",
        cancer_type="Pan-Cancer",
        primary_bulk_target="PanCancer-Combined",
        state_ids=state_ids,
        state_cohorts=tuple("Pan-Cohort" for _ in state_ids),
        milo_logfc=np.linspace(1.8, -1.8, k_states),
        milo_se=np.full(k_states, 0.25),
        bank_transfer_weights=np.eye(k_states),
        native_genes=("GNLY", "PRF1", "FAP", "TGFB1"),
        native_phi_mat=np.full((k_states, 4), 0.25),
    )

    summary_df, cells_df = benchmark_zero_leakage_gating_strategies(
        theta_raw_meta42=theta_raw_mat,
        bulk_bundle=bulk_bundle,
        meta42_spec=ref_spec,
        out_dir=tmp_path,
        n_perm=40,
        seed=19,
    )
    assert summary_df.height == 7
    assert cells_df.height == 63
    assert (tmp_path / "oof_2x2_extreme_regimes_cells.parquet").exists()
    assert (tmp_path / "oof_2x2_extreme_regimes_contrasts.parquet").exists()
    syn_oof_cells = pl.read_parquet(tmp_path / "oof_2x2_extreme_regimes_cells.parquet")
    syn_oof_cons = pl.read_parquet(tmp_path / "oof_2x2_extreme_regimes_contrasts.parquet")
    assert syn_oof_cells.height == 4
    assert syn_oof_cons.height == 5
    assert (tmp_path / "oof_9a_9b_patient_landscape_and_curves.parquet").exists()
    assert (tmp_path / "oof_9a_9b_clinical_and_histology_breakdown.parquet").exists()
    assert (tmp_path / "oof_9a_9b_fold_state_drivers.parquet").exists()
    syn_land = pl.read_parquet(tmp_path / "oof_9a_9b_patient_landscape_and_curves.parquet")
    syn_clin = pl.read_parquet(tmp_path / "oof_9a_9b_clinical_and_histology_breakdown.parquet")
    syn_drv = pl.read_parquet(tmp_path / "oof_9a_9b_fold_state_drivers.parquet")
    assert syn_land.filter(pl.col("record_type") == "patient").height == n_samples
    assert syn_clin.filter(pl.col("record_type") == "km_logrank").height == 3
    assert syn_drv.height == k_states

    leg_row = summary_df.filter(pl.col("method_id") == "legacy_explicit_y").row(0, named=True)
    assert float(leg_row["perm_null_mean_delta_t3_pct"]) > 15.0

    zl_rows = summary_df.filter(pl.col("is_zero_leakage"))
    assert zl_rows.height == 5
    for r in zl_rows.iter_rows(named=True):
        assert abs(float(r["perm_null_mean_delta_t3_pct"])) < 8.0
        assert 0.0 < float(r["perm_pvalue"]) <= 1.0

    for c_row in cells_df.iter_rows(named=True):
        assert (
            0.0
            <= float(c_row["ci_low_pct"])
            <= float(c_row["responder_rate_pct"])
            <= float(c_row["ci_high_pct"])
            <= 100.0
        )


def test_generated_benchmark_parquets_and_figures() -> None:
    """Verify empirical unblended parquet outputs, simplex invariants, Typst captions, and all 14 180mm SVGs."""
    out_dir = Path("output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark")
    results_dir = Path("output/manuscript_pipeline/4b_plot_iatlas_instaprism_milopy_figures")
    article_fig_dir = Path("article/figures/milo-deconv")
    captions_dir = Path("article/captions")

    pred_path = out_dir / "response_prediction_metrics.parquet"
    conc_path = out_dir / "milopy_iatlas_concordance_summary.parquet"
    states_path = out_dir / "milopy_iatlas_state_effects.parquet"
    loco_path = out_dir / "loco_prediction_metrics.parquet"
    deg_path = out_dir / "consensus_metastate_top_degs.parquet"
    all42_deg_path = out_dir / "consensus_metastate_all42_top_degs.parquet"
    all42_pure_inter_path = out_dir / "consensus_metastate_all42_interactions_pure_clr.parquet"
    zl_summary_path = out_dir / "zero_leakage_gating_methods_summary.parquet"
    zl_cells_path = out_dir / "zero_leakage_gating_regimes_cells.parquet"
    oof_2x2_cells_path = out_dir / "oof_2x2_extreme_regimes_cells.parquet"
    oof_2x2_contrasts_path = out_dir / "oof_2x2_extreme_regimes_contrasts.parquet"
    epi_strat_path = out_dir / "epistasis_discovery_strategies_summary.parquet"
    epi_sig_path = out_dir / "epistasis_significant_pairs_calibrated.parquet"
    oof_9a9b_land_path = out_dir / "oof_9a_9b_patient_landscape_and_curves.parquet"
    oof_9a9b_clin_path = out_dir / "oof_9a_9b_clinical_and_histology_breakdown.parquet"
    oof_9a9b_drv_path = out_dir / "oof_9a_9b_fold_state_drivers.parquet"
    raw_frac_path = out_dir / "fractions" / "raw_instaprism_Consensus-Meta42_on_PanCancer-Combined.parquet"

    for req_p in (
        pred_path,
        conc_path,
        states_path,
        loco_path,
        deg_path,
        all42_deg_path,
        all42_pure_inter_path,
        zl_summary_path,
        zl_cells_path,
        oof_2x2_cells_path,
        oof_2x2_contrasts_path,
        epi_strat_path,
        epi_sig_path,
        oof_9a9b_land_path,
        oof_9a9b_clin_path,
        oof_9a9b_drv_path,
        raw_frac_path,
    ):
        assert req_p.exists(), f"Missing {req_p}"

    # Verify all obsolete Hybrid α=0.40 parquets are removed
    for obs_name in (
        "consensus_metastate_interactions.parquet",
        "consensus_metastate_regimes.parquet",
        "consensus_metastate_all42_interactions.parquet",
        "consensus_metastate_all42_hub_summary.parquet",
        "consensus_metastate_all42_regimes.parquet",
        "consensus_metastate_all42_hub_summary_pure_clr.parquet",
        "consensus_metastate_all42_regimes_pure_clr.parquet",
        "consensus_metastate_mixture_sweep.parquet",
        "consensus_metastate_mixture_sibling_collinearity.parquet",
        "consensus_metastate_mixture_regimes_comparison.parquet",
    ):
        assert not (out_dir / obs_name).exists(), f"Obsolete hybrid file still exists: {obs_name}"

    pred_df = pl.read_parquet(pred_path)
    conc_df = pl.read_parquet(conc_path)
    states_df = pl.read_parquet(states_path)
    loco_df = pl.read_parquet(loco_path)
    deg_df = pl.read_parquet(deg_path)
    all42_deg_df = pl.read_parquet(all42_deg_path)
    all42_pure_inter_df = pl.read_parquet(all42_pure_inter_path)
    zl_summary_df = pl.read_parquet(zl_summary_path)
    zl_cells_df = pl.read_parquet(zl_cells_path)
    oof_2x2_cells_df = pl.read_parquet(oof_2x2_cells_path)
    oof_2x2_contrasts_df = pl.read_parquet(oof_2x2_contrasts_path)
    epi_strat_df = pl.read_parquet(epi_strat_path)
    epi_sig_df = pl.read_parquet(epi_sig_path)
    oof_9a9b_land_df = pl.read_parquet(oof_9a9b_land_path)
    oof_9a9b_clin_df = pl.read_parquet(oof_9a9b_clin_path)
    oof_9a9b_drv_df = pl.read_parquet(oof_9a9b_drv_path)
    raw_frac_df = pl.read_parquet(raw_frac_path)

    # All 13 references evaluated (9 individual SC + 2 indication-unified + Unified-All9 + Consensus-Meta42)
    assert pred_df["reference_id"].n_unique() == 13
    assert conc_df["reference_id"].n_unique() == 13
    assert states_df.height > 350
    assert loco_df.height >= 33
    assert deg_df.height == 56  # 14 selected states x 4 top marker genes
    assert all42_deg_df.height == 168  # all 42 states x 4 top marker genes
    assert all42_deg_df["cell_state"].n_unique() == 42

    # Verify Pure Deconvolved CLR (α = 0.00) 42x42 interaction table
    assert all42_pure_inter_df.height == 1764
    assert all42_pure_inter_df["state_a"].n_unique() == 42
    assert all42_pure_inter_df["state_b"].n_unique() == 42
    assert all(
        "Pure Deconvolved CLR" in str(v)
        for v in all42_pure_inter_df["feature_representation"].to_list()
    )
    req_cal_cols = {
        "z_spa",
        "p_spa",
        "q_bh_spa",
        "p_param_boot",
        "q_emp_boot",
        "p_fwer_boot",
        "wald_hc3_2state_t",
        "p_perm_label_2state",
        "lambda_gc_spa",
        "lambda_gc_2state",
        "omnibus_boot_p",
    }
    assert req_cal_cols.issubset(set(all42_pure_inter_df.columns))
    diag_pure = all42_pure_inter_df.filter(pl.col("state_a") == pl.col("state_b"))
    off_pure = all42_pure_inter_df.filter(pl.col("state_a") != pl.col("state_b"))
    assert diag_pure.height == 42
    assert off_pure.height == 1722
    assert all(float(v) == 1.0 for v in diag_pure["p_spa"].to_list())
    assert all(0.0 < float(v) <= 1.0 for v in off_pure["p_spa"].to_list())
    assert all(0.0 < float(v) <= 1.0 for v in off_pure["p_param_boot"].to_list())
    lam_spa = float(off_pure["lambda_gc_spa"][0])
    lam_2s = float(off_pure["lambda_gc_2state"][0])
    assert 0.5 < lam_spa < 1.35
    assert lam_spa < lam_2s

    # Verify Zero-Leakage 3x3 Stromal Gating Benchmark & Within-Trial Permutation Null invariants (unblended)
    assert zl_summary_df.height == 7
    assert zl_cells_df.height == 63
    assert raw_frac_df.height == 1015
    raw_state_cols = [
        c for c in raw_frac_df.columns if c not in {"sample_id", "cohort", "response_binary"}
    ]
    assert len(raw_state_cols) == 42
    raw_mat = raw_frac_df.select(raw_state_cols).to_numpy()
    assert np.all(raw_mat >= 0.0)
    np.testing.assert_allclose(raw_mat.sum(axis=1), np.ones(raw_mat.shape[0]), atol=1e-5)

    leg_emp = zl_summary_df.filter(pl.col("method_id") == "legacy_explicit_y").row(0, named=True)
    ins_emp = zl_summary_df.filter(pl.col("method_id") == "insample_bulk_beta").row(0, named=True)
    assert float(leg_emp["perm_null_mean_delta_t3_pct"]) > 25.0
    assert float(ins_emp["perm_null_mean_delta_t3_pct"]) > 4.0

    zl_emp_df = zl_summary_df.filter(pl.col("is_zero_leakage"))
    assert zl_emp_df.height == 5
    for zl_r in zl_emp_df.iter_rows(named=True):
        assert abs(float(zl_r["perm_null_mean_delta_t3_pct"])) < 4.0
        assert float(zl_r["observed_delta_t3_pct"]) > 5.0

    for strong_id in ("external_sc_milopy", "honest_5fold_oof", "nested_9trial_loco"):
        s_row = zl_summary_df.filter(pl.col("method_id") == strong_id).row(0, named=True)
        assert float(s_row["observed_delta_t3_pct"]) > 15.0
        assert float(s_row["perm_pvalue"]) <= 0.01

    # Verify Honest 5-Fold OOF 2x2 Extreme-Regime Cells (4 corners) & Contrasts (5 tests)
    assert oof_2x2_cells_df.height == 4
    assert oof_2x2_contrasts_df.height == 5
    for c2_row in oof_2x2_contrasts_df.iter_rows(named=True):
        assert abs(float(c2_row["perm_null_mean_pct"])) < 2.0
    c_perm_row = oof_2x2_contrasts_df.filter(
        pl.col("contrast_id") == "permissive_high_vs_low_r"
    ).row(0, named=True)
    c_high_row = oof_2x2_contrasts_df.filter(
        pl.col("contrast_id") == "high_r_permissive_vs_excluded"
    ).row(0, named=True)
    c_did_row = oof_2x2_contrasts_df.filter(
        pl.col("contrast_id") == "did_2x2_gating_interaction"
    ).row(0, named=True)
    c_excl_row = oof_2x2_contrasts_df.filter(
        pl.col("contrast_id") == "excluded_high_vs_low_r"
    ).row(0, named=True)
    c_low_row = oof_2x2_contrasts_df.filter(
        pl.col("contrast_id") == "low_r_permissive_vs_excluded"
    ).row(0, named=True)
    assert float(c_perm_row["fisher_two_sided_p"]) < 0.005
    assert float(c_high_row["fisher_two_sided_p"]) < 0.005
    assert float(c_did_row["wald_two_sided_p"]) < 0.05
    assert float(c_excl_row["fisher_two_sided_p"]) > 0.50
    assert float(c_low_row["fisher_two_sided_p"]) > 0.50

    # Verify Non-Hybrid Multi-Strategy Calibrated Epistasis Discovery Benchmark (6 strategies, 1002 sig/volcano/curve rows)
    assert epi_strat_df.height == 6
    assert epi_sig_df.height == 1002
    assert epi_sig_df.filter(pl.col("record_type") == "tertile_2x2_cell").height == 4
    assert epi_sig_df.filter(pl.col("record_type") == "indication_forest").height == 5
    assert epi_sig_df.filter(pl.col("record_type") == "all861_pair_volcano").height == 861
    assert epi_sig_df.filter(pl.col("record_type") == "m29_dose_curve").height == 122
    assert epi_sig_df.filter(pl.col("record_type") == "m29_quintile_bin").height == 10
    volc_top = (
        epi_sig_df.filter(pl.col("record_type") == "all861_pair_volcano")
        .sort("p_spa")
        .row(0, named=True)
    )
    assert str(volc_top["pair_id"]) == "M29_x_M26"
    volc_tnk = epi_sig_df.filter(
        (pl.col("record_type") == "all861_pair_volcano") & (pl.col("stratum_id") == "t_nk_pair")
    )
    assert all(float(qv) >= 0.35 for qv in volc_tnk["q_bh_all861"].to_list())
    m29_lo_row = epi_sig_df.filter(
        (pl.col("record_type") == "m29_dose_curve") & (pl.col("stratum_id") == "low_emt_m26_t1")
    ).row(0, named=True)
    m29_hi_row = epi_sig_df.filter(
        (pl.col("record_type") == "m29_dose_curve") & (pl.col("stratum_id") == "high_emt_m26_t3")
    ).row(0, named=True)
    assert float(m29_lo_row["beta_inter"]) > 1.50
    assert float(m29_lo_row["p_spa"]) < 1e-4
    assert float(m29_hi_row["p_spa"]) > 0.50
    inflated_rows = epi_strat_df.filter(~pl.col("is_zero_leakage_calibrated"))
    calibrated_rows = epi_strat_df.filter(pl.col("is_zero_leakage_calibrated"))
    assert inflated_rows.height == 2
    assert calibrated_rows.height == 4
    naive_row = inflated_rows.filter(pl.col("strategy_id") == "naive_2state_perm_clr").row(
        0, named=True
    )
    peek_row = inflated_rows.filter(pl.col("strategy_id") == "outcome_peeking_top12_screen").row(
        0, named=True
    )
    assert float(naive_row["null_b_fpr_p05"]) > 0.35
    assert float(peek_row["null_b_fpr_p05"]) > 0.12
    for cal_r in calibrated_rows.iter_rows(named=True):
        assert float(cal_r["perm_null_fpr_p05"]) < 0.08
    assert calibrated_rows.filter(pl.col("emp_min_q_bh") < 0.05).height == 3
    sc_deg_row = calibrated_rows.filter(pl.col("strategy_id") == "sc_deg_100_all861").row(
        0, named=True
    )
    bip_row = calibrated_rows.filter(pl.col("strategy_id") == "sc_deg_100_milo_bipartite").row(
        0, named=True
    )
    oof_gate_row = calibrated_rows.filter(pl.col("strategy_id") == "oof_2x2_composite_gate").row(
        0, named=True
    )
    assert float(sc_deg_row["emp_min_q_bh"]) < 0.05
    assert float(sc_deg_row["emp_min_fwer_p"]) < 0.05
    assert float(bip_row["emp_min_q_bh"]) < 0.01
    assert float(oof_gate_row["emp_min_q_bh"]) < 0.005

    # Verify OOF Hits 9a & 9b continuous dose-response, 4/4 histology replication, OS & 5-fold drivers
    assert oof_9a9b_land_df.filter(pl.col("record_type") == "patient").height == 1015
    slope_9a = oof_9a9b_land_df.filter(
        (pl.col("record_type") == "dose_curve_9a") & (pl.col("stratum_id") == "permissive_nr_t1")
    ).row(0, named=True)
    slope_9b = oof_9a9b_land_df.filter(
        (pl.col("record_type") == "dose_curve_9b") & (pl.col("stratum_id") == "high_r_t3")
    ).row(0, named=True)
    assert float(slope_9a["wald_p"]) < 0.001
    assert float(slope_9a["or_per_sd"]) > 1.25
    assert float(slope_9b["wald_p"]) < 0.005
    assert float(slope_9b["or_per_sd"]) < 0.80

    hist_9a = oof_9a9b_clin_df.filter(
        (pl.col("record_type") == "histology_contrast")
        & (pl.col("hit_id") == "hit_9a")
        & (pl.col("group_id") != "pancancer")
    )
    hist_9b = oof_9a9b_clin_df.filter(
        (pl.col("record_type") == "histology_contrast")
        & (pl.col("hit_id") == "hit_9b")
        & (pl.col("group_id") != "pancancer")
    )
    assert hist_9a.height == 4
    assert hist_9b.height == 4
    assert all(float(d) > 0.0 for d in hist_9a["delta_pct"].to_list())
    assert all(float(d) > 0.0 for d in hist_9b["delta_pct"].to_list())

    lr_9a = oof_9a9b_clin_df.filter(
        (pl.col("record_type") == "km_logrank") & (pl.col("hit_id") == "hit_9a")
    ).row(0, named=True)
    lr_9b = oof_9a9b_clin_df.filter(
        (pl.col("record_type") == "km_logrank") & (pl.col("hit_id") == "hit_9b")
    ).row(0, named=True)
    assert float(lr_9a["p_value"]) < 0.005
    assert float(lr_9b["p_value"]) < 0.005
    assert oof_9a9b_drv_df.height == 42

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

    # Verify all 14 SVG figures exist, parse as well-formed XML, width="180.0mm", height <= 245mm, font >= 10pt
    legacy_titled_stems = (
        "iatlas_response_prediction_auc_heatmap_and_bars",
        "iatlas_milopy_vs_deconv_concordance_grid",
        "iatlas_consensus_metastate_forest_and_tme_profile",
        "iatlas_classifier_and_reference_tier_pareto",
        "iatlas_cross_cohort_concordance_transfer_matrix",
        "iatlas_consensus_metastate_top_deg_markers",
        "iatlas_response_prediction_prauc_heatmap_and_bars",
        "iatlas_consensus_metastate_all42_forest_and_tme_profile",
        "iatlas_consensus_metastate_all42_top_deg_markers",
    )
    heading_free_stems = (
        "iatlas_epistasis_and_stromal_gating_main",
        "iatlas_nonlinear_conditional_interactions_all42_pure_clr",
        "iatlas_epistasis_discovery_power_benchmark",
        "iatlas_zero_leakage_gating_comparison",
        "iatlas_oof_supplementary_diagnostics",
    )
    for stem in legacy_titled_stems + heading_free_stems:
        svg_p = results_dir / f"{stem}.svg"
        png_p = results_dir / f"{stem}.png"
        art_svg_p = article_fig_dir / f"{stem}.svg"
        art_png_p = article_fig_dir / f"{stem}.png"
        assert svg_p.exists() and art_svg_p.exists(), f"Missing SVG figure: {stem}"
        assert png_p.exists() and png_p.stat().st_size > 10000, f"Missing PNG figure: {png_p}"
        assert art_png_p.exists() and art_png_p.stat().st_size > 10000, f"Missing article PNG: {art_png_p}"

        svg_text = svg_p.read_text(encoding="utf-8")
        root = ET.fromstring(svg_text)
        assert root.attrib.get("width") == "180.0mm"
        h_attr = root.attrib.get("height", "0.0mm").replace("mm", "")
        assert float(h_attr) <= 245.0, f"Height {h_attr}mm in {svg_p.name} exceeds 245mm page limit"
        if stem in heading_free_stems:
            assert "role-title" not in svg_text, f"Heading/title still present in {svg_p.name}"
        else:
            assert "role-title-subtitle" in svg_text, f"Missing subtitle in {svg_p.name}"

        font_matches = re.findall(r"(?:font-size:\s*|font-size=\")([0-9.]+)px", svg_text)
        assert len(font_matches) > 0
        for fm in font_matches:
            assert float(fm) >= 13.3, f"Font size {fm}px in {svg_p.name} is below 10pt (13.33px)"

    # Verify standalone importable Typst caption files exist in article/captions/ and results_dir
    main_cap = captions_dir / "iatlas_epistasis_and_stromal_gating_main_caption.typ"
    supp_cap = captions_dir / "iatlas_epistasis_supplementary_captions.typ"
    assert main_cap.exists() and (results_dir / main_cap.name).exists()
    assert supp_cap.exists() and (results_dir / supp_cap.name).exists()
    assert "#let caption_epistasis_stromal_gating_main" in main_cap.read_text(encoding="utf-8")
    supp_txt = supp_cap.read_text(encoding="utf-8")
    for sym in (
        "#let caption_supp_pure_clr_all42",
        "#let caption_supp_epistasis_power_benchmark",
        "#let caption_supp_zero_leakage_gating",
        "#let caption_supp_oof_diagnostics",
    ):
        assert sym in supp_txt
    assert is_successful(Success(True))


def test_tme_datasets_pipeline_dag_svg_and_png() -> None:
    """Verify the tme_datasets provenance & script-dependency DAG SVG and 300 DPI PNG."""
    svg_paths = (
        Path(
            "output/manuscript_pipeline/4e_generate_tme_datasets_pipeline_dag/figure_tme_datasets_pipeline_dag.svg"
        ),
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

