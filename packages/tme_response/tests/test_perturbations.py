"""Unit and property tests for transcriptomic, compositional, and clinical perturbations."""

from __future__ import annotations

import numpy as np
import polars as pl
from returns.maybe import Nothing
from returns.result import Success

from tme_response.perturbations import (
    DilutionConfig,
    DropoutConfig,
    JitterConfig,
    LabelNoiseConfig,
    apply_expression_jitter,
    apply_gene_dropout,
    apply_immune_dilution,
    apply_label_noise,
    compute_perturbation_resilience,
)
from tme_response.schemas import MultiOmicCohort


def make_mock_cohort() -> MultiOmicCohort:
    """Construct a clean minimal MultiOmicCohort for perturbation testing."""
    sample_ids = ("S1", "S2", "S3", "S4")
    expr_df = pl.DataFrame({
        "sample_id": list(sample_ids),
        "GZMA": [10.0, 20.0, 30.0, 40.0],
        "PRF1": [5.0, 15.0, 25.0, 35.0],
        "CXCL9": [2.0, 4.0, 8.0, 16.0],
        "HOUSEKEEPING": [100.0, 100.0, 100.0, 100.0],
    })
    clin_df = pl.DataFrame({
        "sample_id": list(sample_ids),
        "response_binary": [0.0, 0.0, 1.0, 1.0],
        "response_recist": ["PD", "PD", "PR", "CR"],
        "biopsy_timepoint": ["Pre", "Pre", "Pre", "Pre"],
    })
    return MultiOmicCohort(
        cohort_id="TestCohort",
        cancer_type="TestCancer",
        sample_ids=sample_ids,
        expression_tpm=expr_df,
        clinical_annotations=clin_df,
        tmb_scores=Nothing,
        driver_mutations=Nothing,
        cna_scores=Nothing,
    )


def test_apply_expression_jitter() -> None:
    cohort = make_mock_cohort()

    # 1. Zero jitter = identity
    p0 = apply_expression_jitter(cohort, JitterConfig(sigma=0.0, seed=42))
    assert (p0.expression_tpm.select(["GZMA", "PRF1"]).to_numpy() == cohort.expression_tpm.select(["GZMA", "PRF1"]).to_numpy()).all()

    # 2. Non-zero jitter preserves non-negativity and perturbs values
    p1 = apply_expression_jitter(cohort, JitterConfig(sigma=0.5, seed=42))
    mat1 = p1.expression_tpm.select(["GZMA", "PRF1"]).to_numpy()
    assert (mat1 >= 0.0).all()
    assert not np.allclose(mat1, cohort.expression_tpm.select(["GZMA", "PRF1"]).to_numpy())

    # 3. Determinism with same seed
    p2 = apply_expression_jitter(cohort, JitterConfig(sigma=0.5, seed=42))
    mat2 = p2.expression_tpm.select(["GZMA", "PRF1"]).to_numpy()
    assert np.allclose(mat1, mat2)


def test_apply_gene_dropout() -> None:
    cohort = make_mock_cohort()

    # 1. Rate 0.0 = identity
    p0 = apply_gene_dropout(cohort, DropoutConfig(dropout_rate=0.0, seed=42))
    assert (p0.expression_tpm.select(["GZMA"]).to_numpy() == cohort.expression_tpm.select(["GZMA"]).to_numpy()).all()

    # 2. Rate 1.0 = all zeros
    p1 = apply_gene_dropout(cohort, DropoutConfig(dropout_rate=1.0, seed=42))
    mat1 = p1.expression_tpm.select(["GZMA", "PRF1", "CXCL9", "HOUSEKEEPING"]).to_numpy()
    assert (mat1 == 0.0).all()


def test_apply_immune_dilution() -> None:
    cohort = make_mock_cohort()

    # Dilute GZMA, PRF1, CXCL9 by 0.5; HOUSEKEEPING should remain 100.0
    p = apply_immune_dilution(cohort, DilutionConfig(dilution_factor=0.5, seed=42), immune_genes=("GZMA", "PRF1", "CXCL9"))
    gzma_orig = cohort.expression_tpm["GZMA"].to_numpy()
    gzma_diluted = p.expression_tpm["GZMA"].to_numpy()
    assert np.allclose(gzma_diluted, gzma_orig * 0.5)

    hk_orig = cohort.expression_tpm["HOUSEKEEPING"].to_numpy()
    hk_after = p.expression_tpm["HOUSEKEEPING"].to_numpy()
    assert np.allclose(hk_orig, hk_after)


def test_apply_label_noise() -> None:
    clin_df = pl.DataFrame({
        "sample_id": [f"S{i}" for i in range(100)],
        "response_binary": [1.0] * 50 + [0.0] * 50,
    })

    # Rate 0.0 = unchanged
    c0 = apply_label_noise(clin_df, LabelNoiseConfig(noise_rate=0.0, seed=42))
    assert (c0["response_binary"].to_numpy() == clin_df["response_binary"].to_numpy()).all()

    # Rate 1.0 = full flip
    c1 = apply_label_noise(clin_df, LabelNoiseConfig(noise_rate=1.0, seed=42))
    assert np.allclose(c1["response_binary"].to_numpy(), 1.0 - clin_df["response_binary"].to_numpy())

    # Rate 0.2 = ~20 flips
    c_noisy = apply_label_noise(clin_df, LabelNoiseConfig(noise_rate=0.2, seed=42))
    flips = np.sum(c_noisy["response_binary"].to_numpy() != clin_df["response_binary"].to_numpy())
    assert 10 <= flips <= 30


def test_compute_perturbation_resilience() -> None:
    # Predictor 1: Constant robust performance
    records = []
    for s in [0.0, 0.5, 1.0]:
        records.append({
            "cohort_id": "C1",
            "predictor_name": "RobustPred",
            "category": "Signature",
            "perturbation_type": "jitter",
            "intensity": s,
            "roc_auc": 0.80,
        })
    # Predictor 2: Decays linearly from 0.80 to 0.50
    for s, auc_val in zip([0.0, 0.5, 1.0], [0.80, 0.65, 0.50]):
        records.append({
            "cohort_id": "C1",
            "predictor_name": "FragilePred",
            "category": "SingleGene",
            "perturbation_type": "jitter",
            "intensity": s,
            "roc_auc": auc_val,
        })

    df = pl.DataFrame(records)
    match compute_perturbation_resilience(df):
        case Success(res_df):
            assert len(res_df) == 2
            d = dict(zip(res_df["predictor_name"].to_list(), res_df["pri_score"].to_list()))
            assert np.isclose(d["RobustPred"], 1.0, atol=1e-3)
            # Area under [0.80, 0.65, 0.50] from 0 to 1 = 0.65. Ideal = 0.80. PRI = 0.65 / 0.80 = 0.8125
            assert np.isclose(d["FragilePred"], 0.8125, atol=1e-3)
            assert d["RobustPred"] > d["FragilePred"]
        case _ as fail:
            assert False, f"Resilience calculation failed: {fail}"


def test_pan_cohort_visualizations() -> None:
    from tme_response.visualization import (
        create_cross_cohort_resilience_heatmap,
        create_meta_analytic_decay_chart,
        create_pan_cohort_faceted_decay_chart,
    )

    records = []
    for cid in ["CohortA", "CohortB"]:
        for s in [0.0, 0.5, 1.0]:
            records.append({
                "cohort_id": cid,
                "predictor_name": "GEP",
                "category": "Signature",
                "perturbation_type": "jitter",
                "intensity": s,
                "roc_auc": 0.80 - s * 0.1,
                "roc_auc_ci_lower": 0.70 - s * 0.1,
                "roc_auc_ci_upper": 0.90 - s * 0.1,
                "pr_auc": 0.60,
                "pr_auc_ci_lower": 0.50,
                "pr_auc_ci_upper": 0.70,
            })
    df = pl.DataFrame(records)

    # 1. Faceted chart (test all 3 uncertainty modes)
    faceted_shaded = create_pan_cohort_faceted_decay_chart(df, "jitter", uncertainty="shaded")
    assert isinstance(faceted_shaded, Success)
    faceted_bars = create_pan_cohort_faceted_decay_chart(df, "jitter", uncertainty="errorbars")
    assert isinstance(faceted_bars, Success)
    faceted_none = create_pan_cohort_faceted_decay_chart(df, "jitter", uncertainty="none")
    assert isinstance(faceted_none, Success)

    # 2. Meta-analytic chart (test all 3 uncertainty modes)
    meta_shaded = create_meta_analytic_decay_chart(df, "jitter", uncertainty="shaded")
    assert isinstance(meta_shaded, Success)
    meta_bars = create_meta_analytic_decay_chart(df, "jitter", uncertainty="errorbars")
    assert isinstance(meta_bars, Success)
    meta_none = create_meta_analytic_decay_chart(df, "jitter", uncertainty="none")
    assert isinstance(meta_none, Success)

    # 3. Single-cohort decay chart (test all 3 uncertainty modes)
    from tme_response.visualization import create_perturbation_decay_chart
    decay_shaded = create_perturbation_decay_chart(df, "jitter", uncertainty="shaded")
    assert isinstance(decay_shaded, Success)
    decay_bars = create_perturbation_decay_chart(df, "jitter", uncertainty="errorbars")
    assert isinstance(decay_bars, Success)
    decay_none = create_perturbation_decay_chart(df, "jitter", uncertainty="none")
    assert isinstance(decay_none, Success)

    # 4. Cross-cohort heatmap
    res_df = compute_perturbation_resilience(df).unwrap()
    heatmap = create_cross_cohort_resilience_heatmap(res_df)
    assert isinstance(heatmap, Success)

