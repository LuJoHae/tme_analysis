"""Curated baseline immuno-oncology response predictors from COMPASS Table S2 (Shen et al.).

Implements 22 baseline predictors as pure functions returning Result[PredictionResult, str]:
- Single gene targets: PD1 (PDCD1), PDL1 (CD274), CTLA4
- Oligo markers: CD8 duo (CD8A, CD8B), GeneBio target composite
- Multi-gene activation signatures: Davoli CIS, Fehrenbacher Teff, Huang NRS, Ayers IFNG6,
  Jiang CTLs, Roh IS, Wu MIAS, Kong NetBio
- Inverted resistance/exclusion signatures: Jiang TAMs, Jiang Texh, Nurmik CAFs
- Pairwise differential models: Freeman PGM (MAP4K1 - TBX3)
- Dimensionality reduction: Messina CKS (12-chemokine PCA PC1)
- Weighted profile: Cristescu GEP
"""

from __future__ import annotations

import numpy as np
import polars as pl
from returns.result import Failure, Result, Success
from sklearn.decomposition import PCA

from ..schemas import MultiOmicCohort, PredictionResult
from ..types import PredictorCategory

# ==============================================================================
# Master Curated Gene Sets from COMPASS Table S2
# ==============================================================================

GENEBIO_TARGET_GENES: tuple[str, ...] = ("PDCD1", "CD274", "CTLA4")

CD8_DUO_GENES: tuple[str, ...] = ("CD8A", "CD8B")

DAVOLI_CIS_GENES: tuple[str, ...] = (
    "CD247", "CD2", "CD3E", "GZMH", "NKG7", "PRF1", "GZMK",
)

FEHRENBACHER_TEFF_GENES: tuple[str, ...] = (
    "CD8A", "GZMA", "GZMB", "IFNG", "EOMES", "CXCL9", "CXCL10", "TBX21",
)

FREEMAN_PGM_GENES: tuple[str, ...] = ("MAP4K1", "TBX3")

HUANG_NRS_GENES: tuple[str, ...] = (
    "CLEC5A", "TNFSF8", "LILRB2", "FCGR2A", "CLEC7A", "CD86", "CD14", "LAIR1",
    "ITGAM", "FCGR3A", "LILRB4", "CD33", "CCR1", "SIGLEC9", "MARCO", "CD163",
    "EMILIN2", "CCL8", "CCL2", "HOPX", "EFEMP1", "FN1", "SERPING1", "SULF1",
    "KLRC2", "GZMB", "IL2RA", "TNFSF10", "SERPINA1", "PSTPIP2", "GZMH", "SLAMF8",
    "HCLS1", "PDCD1LG2", "TNFAIP8L2", "GZMA", "CCL5", "NKG7", "CCR5", "CLEC12A",
    "CST7", "CXCL11", "CXCL10", "GBP1", "LILRB1", "HCST", "PRF1", "LAG3",
    "IFNG", "TNFSF14", "CXCL13", "SIGLEC10", "IL32", "CXCL9", "VCAM1", "CASP10",
    "CD38", "CCR2", "THY1", "FAP", "CDH11", "RARRES2", "COL1A1", "COL3A1",
    "COL1A2", "MMP1", "IL27", "CLEC4E", "TNFSF18",
)

AYERS_IFNG6_GENES: tuple[str, ...] = (
    "IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG",
)

JIANG_CTLS_GENES: tuple[str, ...] = ("CD8A", "CD8B", "GZMA", "GZMB", "PRF1")

JIANG_TAMS_GENES: tuple[str, ...] = (
    "F13A1", "FCER1A", "CCL17", "FOXQ1", "ESPNL", "CD1A", "GATM", "CCL13",
    "PALLD", "GALNT18", "MAOA", "RAMP1", "STAB1", "CCL26", "CCL23", "PARM1",
    "CD1E", "ITM2C", "CALCRL", "CRH", "RGS18", "MS4A6A", "DHRS2", "PON2",
    "ALOX15", "RAB33A", "MOCOS", "CCL18", "IL17RB", "FABP4", "CMTM8", "QPRT",
    "CDR2L", "DUOXA1", "ABCC4", "SYT17", "PPP1R14A", "PDGFC", "GPT", "FAM189A2",
    "RASAL1", "IPCEF1", "ZNF366", "MAP4K1", "RAB30", "PCED1B", "TMIGD3", "SH3BP4",
    "RRS1", "RNASE1",
)

JIANG_TEXH_GENES: tuple[str, ...] = (
    "HSPA1B", "HSPA1A", "NR4A2", "RGS1", "TNFAIP3", "CAMK2N1", "DUSP1", "NR4A3",
    "IFIT3", "IFIT1B", "FOSB", "CCL5", "SAMD3", "CCL3L3", "MACC1", "VPS37B",
    "TF", "DSEL", "SPATA20", "ADSSL1", "DHX58", "SPP1", "TRIM15", "DUSP26",
    "ABI3", "SFN", "C1QC", "RTP4", "JUN", "BCKDHB", "CCL4", "POU6F1",
    "CYSLTR2", "SLC14A1", "RGS16", "CUEDC1", "DEDD2", "CPEB1", "CTSS", "TMEM88",
    "FAM196B", "SKI", "BTG1", "SYCP2L", "RGS2", "ABCB1", "ADRB2", "CD69",
    "RNF166", "TNN",
)

MESSINA_CKS_GENES: tuple[str, ...] = (
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21",
    "CXCL9", "CXCL10", "CXCL11", "CXCL13",
)

NURMIK_CAFS_GENES: tuple[str, ...] = ("FAP", "ACTA2", "MFAP5", "COL11A1", "TNC")

ROH_IS_GENES: tuple[str, ...] = (
    "GZMA", "GZMB", "PRF1", "GNLY", "HLA-A", "HLA-B", "HLA-C", "HLA-E",
    "HLA-F", "HLA-G", "HLA-H", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
    "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1",
    "IFNG", "IFNGR1", "IFNGR2", "IRF1", "STAT1", "PSMB9", "CCR5", "CCL3",
    "CCL4", "CCL5", "CXCL9", "CXCL10", "CXCL11", "ICAM1", "ICAM2", "ICAM3",
    "ICAM4", "ICAM5", "VCAM1",
)

WU_MIAS_GENES: tuple[str, ...] = (
    "PTPN6", "CD3E", "CD247", "LCK", "CD3D", "CD3G", "PDCD1", "HLA-DRA",
    "HLA-DPB1", "CD4", "TBX21", "HLA-DPA1", "IRF8",
)

CRISTESCU_GEP_WEIGHTS: dict[str, float] = {
    "CCL5": 0.008346,
    "CD27": 0.072293,
    "CD274": 0.042853,
    "CD276": -0.023900,
    "CD8A": 0.031021,
    "CMKLR1": 0.151253,
    "CXCL9": 0.074135,
    "CXCR6": 0.004313,
    "HLA-DQA1": 0.020091,
    "HLA-DRB1": 0.058806,
    "HLA-E": 0.071750,
    "IDO1": 0.060679,
    "LAG3": 0.123895,
    "NKG7": 0.075524,
    "PDCD1LG2": 0.003734,
    "PSMB10": 0.032999,
    "STAT1": 0.250229,
    "TIGIT": 0.084767,
}

KONG_NETBIO_PD1_GENES: tuple[str, ...] = (
    "PDCD1", "CTLA4", "CD274", "HLA-DRB1", "HLA-DRA", "PTPN11", "LAG3", "FOXP3",
    "LCK", "HLA-DPA1", "HLA-DRB5", "HLA-DQA1", "HLA-DQA2", "PTPN6", "HLA-DPB1",
    "HLA-DQB1", "HLA-DQB2", "CD28", "CD80", "CD3G", "CD86", "CD3D", "CD4",
    "PDCD1LG2", "CD3E", "HAVCR2", "BTLA", "CSK", "TBX21", "TNFRSF4", "IDO1",
    "ICOS", "IL2", "TNF", "IL4", "PIK3CA", "B2M", "IL10", "HLA-A", "PIK3R1",
    "IL6", "CD40LG", "PTPRC", "HLA-E", "STAT3", "HLA-G", "CD44", "FYN", "HLA-C",
    "HLA-B", "ITGAM", "HLA-F", "CD40", "CD160", "TNFRSF9", "PIK3R2", "CSF2",
    "AKT1", "EGFR", "ZAP70", "NCAM1", "IRF1", "IRF4", "ICAM1", "ITGAX", "UBA52",
    "STAT1", "FCGR1A", "PLCG1", "RPS27A", "IFNG", "IRF3", "PIK3CB", "SRC",
    "CBL", "IRF7", "UBC", "UBB", "PIK3R3", "CD247", "RAC1", "PTAFR", "GRAP2",
    "VCAM1", "VAV1", "SYK", "CDC42", "IRF9", "IL7R", "IL17A", "HRAS", "IRF2",
    "LCP2", "EGF", "IRF5", "TIGIT", "LYN", "SHC1", "DYNLL1", "MAPK1", "OAS2",
    "TNFRSF18", "OASL", "OAS3", "OAS1", "GZMB", "IL2RB", "IRF6", "GNB1", "TLR4",
    "PRKCQ", "IL13", "IRF8", "PML", "DYNC1H1", "INS", "AP2A1", "AP2M1", "TP53",
    "CIITA", "PIK3CD", "IL7", "AP2A2", "GBP2", "IL2RA", "AP2S1", "JAK2", "CLTC",
    "SP100", "CD19", "CD27", "FCGR1B", "KRAS", "CLTA", "CD69", "GBP1", "MT2A",
    "STAT5A", "GRB2", "RHOA", "IFI30", "CCR7", "AP2B1", "SH3GL2", "DNM1",
    "DNM2", "ARF1", "LAT", "MAPK3", "DYNLL2", "GBP5", "DYNC1I2", "SELL",
    "GBP3", "JAK1", "IL15", "JUN", "VTCN1", "TNFRSF14", "GBP7", "GBP4",
    "GBP6", "NRAS", "YES1", "ITGB1", "SOS1", "DCTN2", "ITK", "STAT5B", "KIF2C",
    "GNG2", "TRIM25", "CXCL10", "DYNC1LI2", "CD8A", "DYNC1LI1", "TRIM21",
    "DNM3", "DYNC1I1", "PTPN22", "VAMP8", "CD2", "CENPE", "TNFSF4", "IL5",
    "KIF11", "ACTR1A", "ACTR10", "PAG1", "DCTN3", "ITGB2", "KIF18A", "DCTN1",
    "TNFSF14", "KIF4A", "SEC13", "KIF2A", "SOCS3", "TRIM45", "CD74",
)


# ==============================================================================
# Pure Functional Scoring Core Helpers
# ==============================================================================

def _compute_mean_log2_score(
    cohort: MultiOmicCohort,
    genes: tuple[str, ...],
    predictor_name: str,
    min_fraction: float = 0.4,
    invert: bool = False,
) -> Result[PredictionResult, str]:
    """Pure helper to calculate the mean log2(TPM + 1) across available signature genes."""
    try:
        expr = cohort.expression_tpm
        cols = set(expr.columns)
        available = [g for g in genes if g in cols]
        min_required = max(1, int(len(genes) * min_fraction))

        if len(available) < min_required:
            return Failure(
                f"Too few genes available for {predictor_name} in {cohort.cohort_id}: "
                f"{len(available)}/{len(genes)} (min required: {min_required})."
            )

        log_cols = [(pl.col(g) + 1.0).log(base=2) for g in available]
        mean_expr = pl.mean_horizontal(log_cols)
        score_expr = (-1.0 * mean_expr) if invert else mean_expr

        res_df = expr.select([
            pl.col("sample_id"),
            score_expr.alias("score"),
        ])

        return Success(
            PredictionResult(
                predictor_name=predictor_name,
                category=PredictorCategory.SIGNATURE,
                predictions=res_df,
            )
        )
    except Exception as exc:
        return Failure(f"Failed to compute {predictor_name} for {cohort.cohort_id}: {exc}")


def _compute_standardized_z_score(
    cohort: MultiOmicCohort,
    genes: tuple[str, ...],
    predictor_name: str,
    weights: dict[str, float] | None = None,
    min_fraction: float = 0.4,
) -> Result[PredictionResult, str]:
    """Pure helper to calculate standardized z-score mean or weighted sum across genes."""
    try:
        expr = cohort.expression_tpm
        cols = set(expr.columns)
        available = [g for g in genes if g in cols]
        min_required = max(1, int(len(genes) * min_fraction))

        if len(available) < min_required:
            return Failure(
                f"Too few genes available for {predictor_name} in {cohort.cohort_id}: "
                f"{len(available)}/{len(genes)} (min required: {min_required})."
            )

        # Standardize each gene column across cohort samples
        z_exprs = []
        weight_vals = []
        for g in available:
            log_col = (pl.col(g) + 1.0).log(base=2)
            z_col = (log_col - log_col.mean()) / log_col.std(ddof=0)
            z_filled = z_col.fill_nan(0.0)
            w = weights.get(g, 1.0) if weights is not None else 1.0
            z_exprs.append(z_filled * w)
            weight_vals.append(w)

        sum_z = z_exprs[0]
        for z in z_exprs[1:]:
            sum_z = sum_z + z

        total_weight = sum(weight_vals) if weights is None else 1.0
        final_score = sum_z / float(total_weight) if total_weight != 0.0 else sum_z

        res_df = expr.select([
            pl.col("sample_id"),
            final_score.alias("score"),
        ])

        return Success(
            PredictionResult(
                predictor_name=predictor_name,
                category=PredictorCategory.SIGNATURE,
                predictions=res_df,
            )
        )
    except Exception as exc:
        return Failure(f"Failed to compute {predictor_name} for {cohort.cohort_id}: {exc}")


# ==============================================================================
# Public Pure Scoring Functions for Table S2 Baselines
# ==============================================================================

def compute_genebio_target_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Kong et al. 2022 GeneBio composite target score (PDCD1 + CD274 + CTLA4)."""
    return _compute_mean_log2_score(
        cohort,
        GENEBIO_TARGET_GENES,
        predictor_name="GeneBio_TargetComposite",
        min_fraction=0.33,
    )


def compute_cd8_duo_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Chen et al. 2016 / Kong et al. 2022 CD8 duo score (CD8A + CD8B)."""
    return _compute_mean_log2_score(
        cohort,
        CD8_DUO_GENES,
        predictor_name="Gene_CD8_Duo",
        min_fraction=0.5,
    )


def compute_davoli_cis_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Davoli et al. 2017 Cytotoxic Immune Signature (CIS) 7-gene mean."""
    return _compute_mean_log2_score(
        cohort,
        DAVOLI_CIS_GENES,
        predictor_name="Davoli_CIS",
        min_fraction=0.5,
    )


def compute_fehrenbacher_teff_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Fehrenbacher et al. 2016 T-effector IFN-gamma signature 8-gene mean."""
    return _compute_mean_log2_score(
        cohort,
        FEHRENBACHER_TEFF_GENES,
        predictor_name="Fehrenbacher_Teff",
        min_fraction=0.5,
    )


def compute_freeman_pgm_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Freeman et al. 2022 Paired Gene Marker score: log2(MAP4K1+1) - log2(TBX3+1)."""
    try:
        expr = cohort.expression_tpm
        cols = set(expr.columns)
        if "MAP4K1" not in cols or "TBX3" not in cols:
            return Failure(
                f"Cohort {cohort.cohort_id} missing MAP4K1 or TBX3 for Freeman PGM."
            )

        score_expr = (
            (pl.col("MAP4K1") + 1.0).log(base=2) - (pl.col("TBX3") + 1.0).log(base=2)
        )

        res_df = expr.select([
            pl.col("sample_id"),
            score_expr.alias("score"),
        ])

        return Success(
            PredictionResult(
                predictor_name="Freeman_PGM",
                category=PredictorCategory.SIGNATURE,
                predictions=res_df,
            )
        )
    except Exception as exc:
        return Failure(f"Failed to compute Freeman PGM for {cohort.cohort_id}: {exc}")


def compute_huang_nrs_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Huang et al. 2019 Neoadjuvant Response Signature (NRS) 69-gene mean."""
    return _compute_mean_log2_score(
        cohort,
        HUANG_NRS_GENES,
        predictor_name="Huang_NRS",
        min_fraction=0.35,
    )


def compute_ayers_ifng6_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Ayers et al. 2017 IFN-gamma 6-gene core signature mean."""
    return _compute_mean_log2_score(
        cohort,
        AYERS_IFNG6_GENES,
        predictor_name="Ayers_IFNG_6",
        min_fraction=0.5,
    )


def compute_jiang_ctls_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Jiang et al. 2018 Cytotoxic T Lymphocyte (CTL) 5-gene core mean."""
    return _compute_mean_log2_score(
        cohort,
        JIANG_CTLS_GENES,
        predictor_name="Jiang_CTLs",
        min_fraction=0.6,
    )


def compute_jiang_tams_score(
    cohort: MultiOmicCohort,
    invert_for_response: bool = True,
) -> Result[PredictionResult, str]:
    """Joyce / Jiang 2018 Tumor-Associated Macrophages (TAM) 50-gene mean (inverted for response)."""
    p_name = "Jiang_TAMs_Inverted" if invert_for_response else "Jiang_TAMs"
    return _compute_mean_log2_score(
        cohort,
        JIANG_TAMS_GENES,
        predictor_name=p_name,
        min_fraction=0.20,
        invert=invert_for_response,
    )


def compute_jiang_texh_score(
    cohort: MultiOmicCohort,
    invert_for_response: bool = True,
) -> Result[PredictionResult, str]:
    """Giordano / Jiang 2018 T-cell exhaustion (Texh) 50-gene mean (inverted for response)."""
    p_name = "Jiang_Texh_Inverted" if invert_for_response else "Jiang_Texh"
    return _compute_mean_log2_score(
        cohort,
        JIANG_TEXH_GENES,
        predictor_name=p_name,
        min_fraction=0.20,
        invert=invert_for_response,
    )


def compute_nurmik_cafs_score(
    cohort: MultiOmicCohort,
    invert_for_response: bool = True,
) -> Result[PredictionResult, str]:
    """Nurmik et al. 2020 Cancer-Associated Fibroblast (CAF) signature (inverted for response)."""
    # Accept TNC or alias TN-C
    expr = cohort.expression_tpm
    cols = set(expr.columns)
    caf_genes = list(NURMIK_CAFS_GENES)
    if "TNC" not in cols and "TN-C" in cols:
        caf_genes = [g if g != "TNC" else "TN-C" for g in caf_genes]

    p_name = "Nurmik_CAFs_Inverted" if invert_for_response else "Nurmik_CAFs"
    return _compute_mean_log2_score(
        cohort,
        tuple(caf_genes),
        predictor_name=p_name,
        min_fraction=0.4,
        invert=invert_for_response,
    )


def compute_roh_is_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Roh et al. 2017 Immune Score (IS) 40-gene mean."""
    return _compute_mean_log2_score(
        cohort,
        ROH_IS_GENES,
        predictor_name="Roh_IS",
        min_fraction=0.25,
    )


def compute_wu_mias_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Wu et al. 2022 MHC-I Association Immunoscore (MIAS) standardized z-score mean."""
    return _compute_standardized_z_score(
        cohort,
        WU_MIAS_GENES,
        predictor_name="Wu_MIAS",
        min_fraction=0.4,
    )


def compute_cristescu_gep_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Cristescu et al. 2018 T-cell inflamed GEP weighted standardized sum."""
    expr = cohort.expression_tpm
    cols = set(expr.columns)
    weights = dict(CRISTESCU_GEP_WEIGHTS)
    # Handle PSMB9 alias for PSMB10 if PSMB10 is absent
    if "PSMB10" not in cols and "PSMB9" in cols:
        weights["PSMB9"] = weights.pop("PSMB10")

    return _compute_standardized_z_score(
        cohort,
        tuple(weights.keys()),
        predictor_name="Cristescu_GEP",
        weights=weights,
        min_fraction=0.4,
    )


def compute_kong_netbio_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Kong et al. 2022 NetBio top-200 target-proximal network diffusion signature mean."""
    return _compute_mean_log2_score(
        cohort,
        KONG_NETBIO_PD1_GENES,
        predictor_name="Kong_NetBio_PD1",
        min_fraction=0.20,
    )


def compute_messina_cks_score(cohort: MultiOmicCohort) -> Result[PredictionResult, str]:
    """Messina et al. 2012 12-Chemokine Signature (CKS) via PCA Component 1.

    Sign is oriented to positively correlate with CXCL9/CXCL10 so higher score denotes response.
    """
    try:
        expr = cohort.expression_tpm
        cols = set(expr.columns)
        available = [g for g in MESSINA_CKS_GENES if g in cols]

        if len(available) < 4:
            return Failure(
                f"Too few Messina CKS chemokines available in {cohort.cohort_id}: "
                f"{len(available)}/{len(MESSINA_CKS_GENES)}."
            )

        # Standardize log2 expression matrix for PCA
        log_mat = np.log2(expr.select(available).to_numpy() + 1.0)
        mean_vec = np.nanmean(log_mat, axis=0, keepdims=True)
        std_vec = np.nanstd(log_mat, axis=0, ddof=0, keepdims=True)
        std_vec[std_vec == 0.0] = 1.0
        z_mat = np.nan_to_num((log_mat - mean_vec) / std_vec, nan=0.0)

        pca = PCA(n_components=1)
        pc1 = pca.fit_transform(z_mat).flatten()

        # Sign orientation: correlate with CXCL9 or CXCL10 if present, else mean of chemokines
        ref_gene = "CXCL9" if "CXCL9" in cols else ("CXCL10" if "CXCL10" in cols else available[0])
        ref_vec = np.log2(expr.select(pl.col(ref_gene)).to_series().to_numpy() + 1.0)

        if np.std(pc1) > 0 and np.std(ref_vec) > 0:
            corr = np.corrcoef(pc1, ref_vec)[0, 1]
            if corr < 0:
                pc1 = -pc1

        res_df = pl.DataFrame({
            "sample_id": expr["sample_id"].to_list(),
            "score": pc1.tolist(),
        })

        return Success(
            PredictionResult(
                predictor_name="Messina_CKS",
                category=PredictorCategory.SIGNATURE,
                predictions=res_df,
            )
        )
    except Exception as exc:
        return Failure(f"Failed to compute Messina CKS for {cohort.cohort_id}: {exc}")
