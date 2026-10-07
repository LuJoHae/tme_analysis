"""Unit tests for pure functional transcriptomic signatures and COMPASS Table S2 baselines."""

import numpy as np
import polars as pl
import pytest
from returns.maybe import Nothing, Some
from returns.result import Success

from tme_response.schemas import MultiOmicCohort
from tme_response.signatures.compass_baselines import (
    compute_ayers_ifng6_score,
    compute_cd8_duo_score,
    compute_cristescu_gep_score,
    compute_davoli_cis_score,
    compute_fehrenbacher_teff_score,
    compute_freeman_pgm_score,
    compute_genebio_target_score,
    compute_huang_nrs_score,
    compute_jiang_ctls_score,
    compute_jiang_tams_score,
    compute_jiang_texh_score,
    compute_kong_netbio_score,
    compute_messina_cks_score,
    compute_nurmik_cafs_score,
    compute_roh_is_score,
    compute_wu_mias_score,
)
from tme_response.signatures.cyt import compute_cyt_score
from tme_response.signatures.gep import compute_gep_score
from tme_response.signatures.impres import compute_impres_score
from tme_response.signatures.ipres import compute_ipres_score
from tme_response.signatures.scoring import compute_standard_signatures
from tme_response.signatures.single_gene import compute_single_gene_score
from tme_response.models.tide import compute_tide_score


@pytest.fixture
def mock_cohort() -> MultiOmicCohort:
    """Construct an analytical toy cohort with known expression values.

    S1: Inflamed / High Immune / Low Exclusion
    S2: Intermediate
    S3: Cold / Low Immune / High Exclusion
    """
    sample_ids = ("S1", "S2", "S3")
    expr_dict: dict[str, list[object]] = {
        "sample_id": list(sample_ids),
        # CYT: GZMA & PRF1
        # S1: GZMA=3, PRF1=7 -> log2(4)=2, log2(8)=3 -> CYT = 2.5
        # S2: GZMA=1, PRF1=1 -> log2(2)=1, log2(2)=1 -> CYT = 1.0
        # S3: GZMA=0, PRF1=0 -> log2(1)=0, log2(1)=0 -> CYT = 0.0
        "GZMA": [3.0, 1.0, 0.0],
        "PRF1": [7.0, 1.0, 0.0],
        # CXCL9 & CD8
        "CXCL9": [15.0, 3.0, 0.0],
        "CD8A": [7.0, 3.0, 1.0],
        "CD8B": [7.0, 3.0, 1.0],
        "GZMB": [7.0, 3.0, 1.0],
        # Target ICI genes
        "PDCD1": [10.0, 1.0, 0.0],
        "CD274": [5.0, 0.0, 0.0],
        "CTLA4": [5.0, 0.0, 0.0],
        # Freeman PGM
        # S1: log2(16)=4 - log2(1)=0 -> 4.0
        # S2: log2(4)=2 - log2(4)=2 -> 0.0
        # S3: log2(1)=0 - log2(16)=4 -> -4.0
        "MAP4K1": [15.0, 3.0, 0.0],
        "TBX3": [0.0, 3.0, 15.0],
        # IMPRES pairs
        "TNFRSF4": [2.0, 5.0, 0.0],
        "CD27": [10.0, 0.0, 0.0],
        "CD40": [1.0, 5.0, 0.0],
        "CD80": [1.0, 1.0, 1.0],
        "CD40LG": [5.0, 0.0, 0.0],
        "CD86": [1.0, 1.0, 1.0],
        "TNFRSF9": [1.0, 1.0, 1.0],
        "ICOSLG": [1.0, 1.0, 1.0],
        "CD28": [5.0, 0.0, 0.0],
        "HAVCR2": [1.0, 5.0, 5.0],
        "ICOS": [5.0, 0.0, 0.0],
        # Ayers GEP / IFNG
        "CCL5": [5.0, 2.0, 1.0],
        "CD276": [2.0, 2.0, 2.0],
        "CMKLR1": [1.0, 1.0, 1.0],
        "CXCR6": [3.0, 1.0, 0.0],
        "HLA-DQA1": [10.0, 5.0, 1.0],
        "HLA-DRB1": [10.0, 5.0, 1.0],
        "HLA-E": [5.0, 5.0, 5.0],
        "IDO1": [2.0, 0.0, 0.0],
        "LAG3": [2.0, 0.0, 0.0],
        "NKG7": [5.0, 2.0, 0.0],
        "PDCD1LG2": [1.0, 1.0, 1.0],
        "PSMB9": [5.0, 2.0, 1.0],
        "PSMB10": [5.0, 2.0, 1.0],
        "STAT1": [5.0, 2.0, 1.0],
        "TIGIT": [2.0, 1.0, 0.0],
        "IFNG": [5.0, 2.0, 0.0],
        "CXCL10": [5.0, 2.0, 0.0],
        "HLA-DRA": [5.0, 2.0, 0.0],
        # Davoli CIS & Teff
        "CD247": [5.0, 2.0, 0.0],
        "CD2": [5.0, 2.0, 0.0],
        "CD3E": [5.0, 2.0, 0.0],
        "GZMH": [5.0, 2.0, 0.0],
        "GZMK": [5.0, 2.0, 0.0],
        "EOMES": [5.0, 2.0, 0.0],
        "TBX21": [5.0, 2.0, 0.0],
        # Messina CKS
        "CCL2": [5.0, 2.0, 0.0],
        "CCL3": [5.0, 2.0, 0.0],
        "CCL4": [5.0, 2.0, 0.0],
        "CCL8": [5.0, 2.0, 0.0],
        "CCL18": [5.0, 2.0, 0.0],
        "CCL19": [5.0, 2.0, 0.0],
        "CCL21": [5.0, 2.0, 0.0],
        "CXCL11": [5.0, 2.0, 0.0],
        "CXCL13": [5.0, 2.0, 0.0],
        # Exclusion / CAF markers
        "ACTA2": [0.0, 5.0, 10.0],
        "COL1A1": [0.0, 5.0, 10.0],
        "FAP": [0.0, 5.0, 10.0],
        "PDGFRB": [0.0, 5.0, 10.0],
        "TGFB1": [0.0, 5.0, 10.0],
        "MFAP5": [0.0, 5.0, 10.0],
        "COL11A1": [0.0, 5.0, 10.0],
        "TNC": [0.0, 5.0, 10.0],
        # MDSC & TAM markers
        "CD14": [0.0, 5.0, 10.0],
        "ITGAM": [0.0, 5.0, 10.0],
        "S100A8": [0.0, 5.0, 10.0],
        "S100A9": [0.0, 5.0, 10.0],
        "STAT3": [0.0, 5.0, 10.0],
        "CD163": [0.0, 5.0, 10.0],
        "MRC1": [0.0, 5.0, 10.0],
        "MS4A4A": [0.0, 5.0, 10.0],
        "VSIG4": [0.0, 5.0, 10.0],
        "F13A1": [0.0, 5.0, 10.0],
        "FCER1A": [0.0, 5.0, 10.0],
        "CCL17": [0.0, 5.0, 10.0],
        "FOXQ1": [0.0, 5.0, 10.0],
        "ESPNL": [0.0, 5.0, 10.0],
        "CD1A": [0.0, 5.0, 10.0],
        "GATM": [0.0, 5.0, 10.0],
        "CCL13": [0.0, 5.0, 10.0],
        "PALLD": [0.0, 5.0, 10.0],
        "GALNT18": [0.0, 5.0, 10.0],
        "MAOA": [0.0, 5.0, 10.0],
        # T-cell exhaustion markers
        "HSPA1B": [0.0, 5.0, 10.0],
        "HSPA1A": [0.0, 5.0, 10.0],
        "NR4A2": [0.0, 5.0, 10.0],
        "RGS1": [0.0, 5.0, 10.0],
        "TNFAIP3": [0.0, 5.0, 10.0],
        "CAMK2N1": [0.0, 5.0, 10.0],
        "DUSP1": [0.0, 5.0, 10.0],
        "NR4A3": [0.0, 5.0, 10.0],
        "IFIT3": [0.0, 5.0, 10.0],
        "IFIT1B": [0.0, 5.0, 10.0],
        "FOSB": [0.0, 5.0, 10.0],
        # Roh IS / Wu MIAS / NRS
        "GNLY": [5.0, 2.0, 0.0],
        "HLA-A": [5.0, 2.0, 0.0],
        "HLA-B": [5.0, 2.0, 0.0],
        "HLA-C": [5.0, 2.0, 0.0],
        "PTPN6": [5.0, 2.0, 0.0],
        "LCK": [5.0, 2.0, 0.0],
        "CD3D": [5.0, 2.0, 0.0],
        "CD3G": [5.0, 2.0, 0.0],
        "CD4": [5.0, 2.0, 0.0],
        "HLA-DPB1": [5.0, 2.0, 0.0],
        "HLA-DPA1": [5.0, 2.0, 0.0],
        "IRF8": [5.0, 2.0, 0.0],
        "CLEC5A": [5.0, 2.0, 0.0],
        "TNFSF8": [5.0, 2.0, 0.0],
        "LILRB2": [5.0, 2.0, 0.0],
        "FCGR2A": [5.0, 2.0, 0.0],
        "CLEC7A": [5.0, 2.0, 0.0],
        "LAIR1": [5.0, 2.0, 0.0],
        # IPRES markers
        "AXL": [0.0, 2.0, 8.0],
        "ROR2": [0.0, 2.0, 8.0],
        "WNT5A": [0.0, 2.0, 8.0],
        "LOXL2": [0.0, 2.0, 8.0],
    }

    clinical_dict = {
        "sample_id": list(sample_ids),
        "response_binary": [1.0, 0.0, 0.0],
        "response_recist": ["PR", "PD", "PD"],
        "biopsy_timepoint": ["Pre", "Pre", "Pre"],
    }

    return MultiOmicCohort(
        cohort_id="ToyCohort",
        cancer_type="TestCancer",
        sample_ids=sample_ids,
        expression_tpm=pl.DataFrame(expr_dict),
        clinical_annotations=pl.DataFrame(clinical_dict),
        tmb_scores=Nothing,
        driver_mutations=Nothing,
        cna_scores=Nothing,
    )


# ==============================================================================
# Existing Core Signature Tests
# ==============================================================================

def test_cyt_analytical_values(mock_cohort: MultiOmicCohort) -> None:
    """Verify CYT calculations match exact mathematical definitions."""
    res = compute_cyt_score(mock_cohort)
    assert isinstance(res, Success)
    df = res.unwrap().predictions
    scores = df["score"].to_list()
    assert np.isclose(scores[0], 2.5)
    assert np.isclose(scores[1], 1.0)
    assert np.isclose(scores[2], 0.0)


def test_impres_score_bounds(mock_cohort: MultiOmicCohort) -> None:
    """Verify IMPRES returns values in [0, 15] and ranks hot tumor highest."""
    res = compute_impres_score(mock_cohort)
    assert isinstance(res, Success)
    df = res.unwrap().predictions
    scores = df["score"].to_list()
    for s in scores:
        assert 0.0 <= s <= 15.0
    assert scores[0] > scores[2]


def test_gep_score_ranking(mock_cohort: MultiOmicCohort) -> None:
    """Verify Ayers GEP ranks inflamed S1 > intermediate S2 > cold S3."""
    res = compute_gep_score(mock_cohort)
    assert isinstance(res, Success)
    scores = res.unwrap().predictions["score"].to_list()
    assert scores[0] > scores[1] > scores[2]


def test_single_gene_cxcl9(mock_cohort: MultiOmicCohort) -> None:
    """Verify single gene log2 extraction."""
    res = compute_single_gene_score(mock_cohort, "CXCL9")
    assert isinstance(res, Success)
    scores = res.unwrap().predictions["score"].to_list()
    assert np.isclose(scores[0], 4.0)
    assert np.isclose(scores[1], 2.0)
    assert np.isclose(scores[2], 0.0)


def test_tide_score_inversion(mock_cohort: MultiOmicCohort) -> None:
    """Verify TIDE inverted score ranks inflamed/responsive tumor highest."""
    res = compute_tide_score(mock_cohort, invert_for_response=True)
    assert isinstance(res, Success)
    scores = res.unwrap().predictions["score"].to_list()
    assert scores[0] > scores[2]


# ==============================================================================
# COMPASS Table S2 Baseline Tests
# ==============================================================================

def test_genebio_target_composite(mock_cohort: MultiOmicCohort) -> None:
    """Verify GeneBio target composite computes mean of PDCD1, CD274, CTLA4."""
    res = compute_genebio_target_score(mock_cohort)
    assert isinstance(res, Success)
    pred = res.unwrap()
    assert pred.predictor_name == "GeneBio_TargetComposite"
    scores = pred.predictions["score"].to_list()
    # S1: log2(11)=3.459, log2(6)=2.585, log2(6)=2.585
    assert scores[0] > scores[1] > scores[2]


def test_cd8_duo(mock_cohort: MultiOmicCohort) -> None:
    """Verify CD8 duo computes mean log2(TPM+1) of CD8A and CD8B."""
    res = compute_cd8_duo_score(mock_cohort)
    assert isinstance(res, Success)
    scores = res.unwrap().predictions["score"].to_list()
    # S1: log2(8)=3.0, S2: log2(4)=2.0, S3: log2(2)=1.0
    assert np.isclose(scores[0], 3.0)
    assert np.isclose(scores[1], 2.0)
    assert np.isclose(scores[2], 1.0)


def test_freeman_pgm(mock_cohort: MultiOmicCohort) -> None:
    """Verify Freeman PGM computes log2(MAP4K1+1) - log2(TBX3+1)."""
    res = compute_freeman_pgm_score(mock_cohort)
    assert isinstance(res, Success)
    scores = res.unwrap().predictions["score"].to_list()
    assert np.isclose(scores[0], 4.0)
    assert np.isclose(scores[1], 0.0)
    assert np.isclose(scores[2], -4.0)


def test_davoli_cis_and_fehrenbacher_teff(mock_cohort: MultiOmicCohort) -> None:
    """Verify Davoli CIS and Fehrenbacher Teff rank hot S1 > intermediate S2 > cold S3."""
    res_cis = compute_davoli_cis_score(mock_cohort)
    assert isinstance(res_cis, Success)
    scores_cis = res_cis.unwrap().predictions["score"].to_list()
    assert scores_cis[0] > scores_cis[1] > scores_cis[2]

    res_teff = compute_fehrenbacher_teff_score(mock_cohort)
    assert isinstance(res_teff, Success)
    scores_teff = res_teff.unwrap().predictions["score"].to_list()
    assert scores_teff[0] > scores_teff[1] > scores_teff[2]


def test_huang_nrs_and_ayers_ifng6(mock_cohort: MultiOmicCohort) -> None:
    """Verify Huang NRS and Ayers IFNG6 compute cleanly and rank correctly."""
    res_nrs = compute_huang_nrs_score(mock_cohort)
    assert isinstance(res_nrs, Success)
    scores_nrs = res_nrs.unwrap().predictions["score"].to_list()
    assert scores_nrs[0] > scores_nrs[1] > scores_nrs[2]

    res_ifng6 = compute_ayers_ifng6_score(mock_cohort)
    assert isinstance(res_ifng6, Success)
    scores_ifng6 = res_ifng6.unwrap().predictions["score"].to_list()
    assert scores_ifng6[0] > scores_ifng6[1] > scores_ifng6[2]


def test_jiang_ctl_and_resistance_inversion(mock_cohort: MultiOmicCohort) -> None:
    """Verify Jiang CTLs positive ranking, and TAMs/Texh/CAFs inverted ranking."""
    res_ctl = compute_jiang_ctls_score(mock_cohort)
    assert isinstance(res_ctl, Success)
    scores_ctl = res_ctl.unwrap().predictions["score"].to_list()
    assert scores_ctl[0] > scores_ctl[1] > scores_ctl[2]

    # Inverted signatures: high exclusion in S3 -> more negative score -> S1 > S2 > S3
    res_tam = compute_jiang_tams_score(mock_cohort, invert_for_response=True)
    assert isinstance(res_tam, Success)
    scores_tam = res_tam.unwrap().predictions["score"].to_list()
    assert scores_tam[0] > scores_tam[1] > scores_tam[2]

    res_texh = compute_jiang_texh_score(mock_cohort, invert_for_response=True)
    assert isinstance(res_texh, Success)
    scores_texh = res_texh.unwrap().predictions["score"].to_list()
    assert scores_texh[0] > scores_texh[1] > scores_texh[2]

    res_caf = compute_nurmik_cafs_score(mock_cohort, invert_for_response=True)
    assert isinstance(res_caf, Success)
    scores_caf = res_caf.unwrap().predictions["score"].to_list()
    assert scores_caf[0] > scores_caf[1] > scores_caf[2]


def test_messina_cks_pca(mock_cohort: MultiOmicCohort) -> None:
    """Verify Messina CKS computes PCA PC1 aligned with chemokine presence."""
    res = compute_messina_cks_score(mock_cohort)
    assert isinstance(res, Success)
    scores = res.unwrap().predictions["score"].to_list()
    assert scores[0] > scores[2]


def test_roh_is_and_wu_mias_and_cristescu_gep(mock_cohort: MultiOmicCohort) -> None:
    """Verify Roh IS, Wu MIAS, and Cristescu GEP rank hot tumor highest."""
    res_is = compute_roh_is_score(mock_cohort)
    assert isinstance(res_is, Success)
    scores_is = res_is.unwrap().predictions["score"].to_list()
    assert scores_is[0] > scores_is[1] > scores_is[2]

    res_mias = compute_wu_mias_score(mock_cohort)
    assert isinstance(res_mias, Success)
    scores_mias = res_mias.unwrap().predictions["score"].to_list()
    assert scores_mias[0] > scores_mias[1] > scores_mias[2]

    res_cgep = compute_cristescu_gep_score(mock_cohort)
    assert isinstance(res_cgep, Success)
    scores_cgep = res_cgep.unwrap().predictions["score"].to_list()
    assert scores_cgep[0] > scores_cgep[1] > scores_cgep[2]


def test_kong_netbio(mock_cohort: MultiOmicCohort) -> None:
    """Verify Kong NetBio calculates mean across available target-proximal genes."""
    res = compute_kong_netbio_score(mock_cohort)
    assert isinstance(res, Success)
    scores = res.unwrap().predictions["score"].to_list()
    assert scores[0] > scores[1] > scores[2]


def test_consolidate_signatures_all_baselines(mock_cohort: MultiOmicCohort) -> None:
    """Verify compute_standard_signatures returns all signatures and Table S2 baselines."""
    res = compute_standard_signatures(mock_cohort)
    assert isinstance(res, Success)
    df = res.unwrap()
    assert "score_CYT" in df.columns
    assert "score_IMPRES" in df.columns
    assert "score_Ayers_GEP" in df.columns
    assert "score_CXCL9" in df.columns
    assert "score_CD8A" in df.columns
    assert "score_PDCD1" in df.columns
    assert "score_GeneBio" in df.columns
    assert "score_CD8_Duo" in df.columns
    assert "score_Davoli_CIS" in df.columns
    assert "score_Fehrenbacher_Teff" in df.columns
    assert "score_Freeman_PGM" in df.columns
    assert "score_Huang_NRS" in df.columns
    assert "score_Ayers_IFNG_6" in df.columns
    assert "score_Jiang_CTLs" in df.columns
    assert "score_Jiang_TAMs_Inverted" in df.columns
    assert "score_Jiang_Texh_Inverted" in df.columns
    assert "score_Messina_CKS" in df.columns
    assert "score_Nurmik_CAFs_Inverted" in df.columns
    assert "score_Roh_IS" in df.columns
    assert "score_Wu_MIAS" in df.columns
    assert "score_Cristescu_GEP" in df.columns
    assert "score_Kong_NetBio_PD1" in df.columns
    assert len(df) == 3
