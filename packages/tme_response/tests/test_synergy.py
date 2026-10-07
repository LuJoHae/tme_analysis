"""Unit tests for DNA-RNA composite synergy and biological gating logic."""

import polars as pl
from returns.maybe import Some
from returns.result import Success

from tme_response.schemas import MultiOmicCohort, PredictionResult
from tme_response.signatures.single_gene import compute_single_gene_score
from tme_response.synergy.dna_rna import compute_dna_rna_composite
from tme_response.synergy.gating import apply_antigen_presentation_gating
from tme_response.types import PredictorCategory


def test_dna_rna_composite_calculation() -> None:
    """Verify DNA-RNA linear combination properly standardizes and weights features."""
    sample_ids = ("S1", "S2", "S3")
    expr_df = pl.DataFrame({
        "sample_id": list(sample_ids),
        "CXCL9": [10.0, 5.0, 1.0],
    })
    clin_df = pl.DataFrame({
        "sample_id": list(sample_ids),
        "response_binary": [1.0, 1.0, 0.0],
    })
    tmb_df = pl.DataFrame({
        "sample_id": list(sample_ids),
        "tmb_per_mb": [20.0, 5.0, 1.0],
    })

    cohort = MultiOmicCohort(
        cohort_id="TestSynergy",
        cancer_type="Melanoma",
        sample_ids=sample_ids,
        expression_tpm=expr_df,
        clinical_annotations=clin_df,
        tmb_scores=Some(tmb_df),
        driver_mutations=Some(pl.DataFrame({
            "sample_id": ["S1"],
            "gene": ["BRAF"],
            "variant_classification": ["Missense_Mutation"],
            "is_inactivating": [False],
        })),
        cna_scores=Some(pl.DataFrame({
            "sample_id": list(sample_ids),
            "cdkn2a_loss": [False, False, True],
            "aneuploidy_score": [0.1, 0.2, 0.8],
        })),
    )

    rna_res = compute_single_gene_score(cohort, "CXCL9")
    assert isinstance(rna_res, Success)

    comp_res = compute_dna_rna_composite(cohort, rna_res.unwrap(), rna_weight=1.0, tmb_weight=0.5)
    assert isinstance(comp_res, Success)
    preds = comp_res.unwrap().predictions
    assert "score" in preds.columns
    assert "rna_zscore" in preds.columns
    assert "tmb_zscore" in preds.columns
    # S1 has highest RNA and highest TMB, so it must have highest composite score
    assert preds["score"][0] > preds["score"][1] > preds["score"][2]


def test_antigen_presentation_gating() -> None:
    """Verify inactivating B2M mutation gates score down to minimum."""
    sample_ids = ("S1", "S2")
    # S1 has high score, but has inactivating B2M mutation
    base_preds = PredictionResult(
        predictor_name="BaseScore",
        category=PredictorCategory.SIGNATURE,
        predictions=pl.DataFrame({
            "sample_id": ["S1", "S2"],
            "score": [10.0, 2.0],
        }),
    )

    driver_df = pl.DataFrame({
        "sample_id": ["S1"],
        "gene": ["B2M"],
        "variant_classification": ["Frame_Shift_Del"],
        "is_inactivating": [True],
    })

    cohort = MultiOmicCohort(
        cohort_id="GatingTest",
        cancer_type="Melanoma",
        sample_ids=sample_ids,
        expression_tpm=pl.DataFrame({"sample_id": list(sample_ids)}),
        clinical_annotations=pl.DataFrame({"sample_id": list(sample_ids)}),
        tmb_scores=Some(pl.DataFrame({"sample_id": list(sample_ids), "tmb_per_mb": [1.0, 1.0]})),
        driver_mutations=Some(driver_df),
        cna_scores=Some(pl.DataFrame({"sample_id": list(sample_ids)})),
    )

    gated_res = apply_antigen_presentation_gating(cohort, base_preds)
    assert isinstance(gated_res, Success)
    gated_df = gated_res.unwrap().predictions

    # S1 was gated down below S2 because of B2M inactivating deletion
    s1_score = gated_df.filter(pl.col("sample_id") == "S1")["score"][0]
    s2_score = gated_df.filter(pl.col("sample_id") == "S2")["score"][0]
    assert s1_score < s2_score
