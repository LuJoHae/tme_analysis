"""Pure functional loaders for Tier 1 scRNA-seq and spatial transcriptomics cohorts.

Covers all 84 distinct Tier 1 cohorts across 9 human solid tumor indications:
- Breast: GSE274141, GSE222859, GSE274139, GSE262288, GSE303346, GSE254991, GSE199219, GSE332708, GSE331487, GSE302453, GSE299267
- CRC: GSE205506, GSE146771, GSE164522, GSE235917, GSE278406, GSE274321, GSE309346, CELLxGENE_2554a654, CELLxGENE_5ee552f5, CELLxGENE_4b5afdf9, GSE188711, GSE336564, GSE216534, CELLxGENE_387acac5, CELLxGENE_ef0d813e, CELLxGENE_2e95d453
- Gastric: GSE239676
- HCC: GSE319709, GSE272347, GSE318418, GSE299340, GSE282343, GSE255830, GSE281110, GSE233405, GSE278324, GSE265770, GSE318420, GSE272348, GSE224411, GSE215428
- HNSCC: GSE296954, GSE301720, GSE247582, GSE296867, GSE327189
- Melanoma: GSE286410, GSE198265, GSE242477, GSE303948, GSE294273, GSE320040, GSE210963, GSE256291, GSE244983, GSE300446, GSE270464, GSE211068
- NSCLC: GSE241934, GSE270148, GSE303762, GSE205049, GSE205354, GSE233203, GSE253718, GSE223779, GSE307811, GSE285888, GSE276139
- PDAC: GSE279781, GSE212966, GSE283206, GSE311788, GSE211644, GSE348275, GSE318413, GSE156405, GSE335452
- ccRCC: GSE285701, GSE304466, GSE223808, GSE220313, GSE254498

All loaders:
- Return Result[ad.AnnData, str] (Railway Oriented Programming).
- Produce cells x genes CSR sparse matrix in adata.X with raw non-negative integer counts.
- Enforce unique gene and cell identifiers.
- Standardize clinical observation metadata in adata.obs (patient_id, sample_id, treatment_status, clinical_response, indication, technology, cell_selection_strategy).
- Tag expression metadata via tag_expression_metadata().
"""

from __future__ import annotations

from enum import Enum
import gzip
from pathlib import Path
import re
import shutil
import tarfile
from typing import Callable, Mapping, Sequence
import urllib.request

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
import scanpy as sc
import scipy.io as sio
import scipy.sparse as sp

from ..download.fetcher import download_single_file
from ..download.geo import download_geo_supplementary
from ..logging import get_logger
from ..paths import find_repo_root
from ..preprocessing.matrix_inspection import tag_expression_metadata
from .tier0_single_cell import (
    _apply_subset_and_subsample,
    _extract_tar_safely,
    _harmonize_single_response,
    _load_10x_triplets_or_h5,
    _parse_geo_series_matrix,
    _restore_raw_anndata,
    _to_csr,
)

logger = get_logger("providers.tier1_single_cell")


# =========================================================================
# Tier 1 Data Models & Cohort Specification
# =========================================================================
class Archetype(str, Enum):
    """Categorical classification of raw vendor single-cell package formats."""
    GEO_10X_TAR = "geo_10x_tar"
    GEO_10X_H5_DIR = "geo_10x_h5_dir"
    CELLXGENE_H5AD = "cellxgene_h5ad"
    FLAT_COUNT_TABLE = "flat_count_table"
    SPATIAL_ASSAY = "spatial_assay"


class Tier1CohortConfig(BaseModel):
    """Immutable metadata and layout configuration for a Tier 1 scRNA-seq cohort."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    accession: str
    title: str
    indication: str
    technology: str
    cell_selection: str
    archetype: Archetype
    tar_filename: Maybe[str] = Nothing
    matrix_filename: Maybe[str] = Nothing
    patient_col: Maybe[str] = Nothing
    sample_col: Maybe[str] = Nothing
    response_col: Maybe[str] = Nothing
    response_map: Maybe[Mapping[str, str]] = Nothing
    treatment_col: Maybe[str] = Nothing
    cell_type_col: Maybe[str] = Nothing
    min_counts: int = 100
    min_genes: int = 50
    gex_only: bool = True
    has_response_labels: bool = True


# Auto-generated Tier 1 Cohort Specifications
TIER1_COHORTS: Mapping[str, Tier1CohortConfig] = {
    "GSE274141": Tier1CohortConfig(
        accession="GSE274141",
        title="Single-Cell RNA Sequencing Identifies Molecular Biomarkers Predicting Late Progression to CDK4/6 Inhibition in Patients with HR+/HER2- Metastatic Breast Cancer [FFPE]",
        indication="Breast",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE274141_read-counts-n54-new.csv.gz"),
        has_response_labels=True,
    ),
    "GSE222859": Tier1CohortConfig(
        accession="GSE222859",
        title="Gene expression profile at single cell level of lymphocytes cells for Pan-T cell analysis",
        indication="Breast",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE222859_matrix.mtx.gz"),
        has_response_labels=True,
    ),
    "GSE274139": Tier1CohortConfig(
        accession="GSE274139",
        title="Single-Cell RNA Sequencing Identifies Molecular Biomarkers Predicting Late Progression to CDK4/6 Inhibition in Patients with HR+/HER2- Metastatic Breast Cancer [tissue]",
        indication="Breast",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE274139_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE262288": Tier1CohortConfig(
        accession="GSE262288",
        title="Single-Cell RNA Sequencing Identifies Molecular Biomarkers Predicting Response to CDK4/6 Inhibition in Metastatic HR+/HER2- Breast Cancer",
        indication="Breast",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE262288_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE303346": Tier1CohortConfig(
        accession="GSE303346",
        title="Single-Cell RNA-Seq Reveals Immunosuppressive Effects of Triple-Negative Breast Cancer–Derived Exosomes on NK Cell Responses",
        indication="Breast",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_H5_DIR,
        tar_filename=Nothing,
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE254991": Tier1CohortConfig(
        accession="GSE254991",
        title="Interleukin-16 Establishes a Th1-dominant Tumor Immune Microenvironment And Potentiates Immune Checkpoint Therapies [scRNA-Seq]",
        indication="Breast",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE254991_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE199219": Tier1CohortConfig(
        accession="GSE199219",
        title="Proteogenomic integration of single-cell RNA and protein analysis identifies novel tumour-infiltrating lymphocyte phenotypes in breast cancer",
        indication="Breast",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE199219_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE332708": Tier1CohortConfig(
        accession="GSE332708",
        title="Hybrid In Vivo Breast Cancer Model Reveals Transcriptomic Insights into Cancer Progression with Age",
        indication="Breast",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_H5_DIR,
        tar_filename=Nothing,
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE331487": Tier1CohortConfig(
        accession="GSE331487",
        title="Immunosuppressive myeloid cells induce mesenchymal-like breast cancer stem cells by a membrane-bound TGF-β1-dependent mechanism [scRNA-seq]",
        indication="Breast",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE331487_counts_raw_integrated_samples.tsv.gz"),
        has_response_labels=True,
    ),
    "GSE302453": Tier1CohortConfig(
        accession="GSE302453",
        title="A single-cell map of intratumoral heterogeneity during combination treatment of anti-PD1 and chemotherapy in triple-negative breast cancer",
        indication="Breast",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE302453_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE299267": Tier1CohortConfig(
        accession="GSE299267",
        title="Single cell gene expression profiling of triple negative breast cancer organoids",
        indication="Breast",
        technology="Single-Nucleus RNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE299267_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE205506": Tier1CohortConfig(
        accession="GSE205506",
        title="Remodeling of the Immune and Stromal Cell Compartment by PD-1 Blockade in Mismatch Repair-Deficient Colorectal Cancer",
        indication="CRC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE205506_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE146771": Tier1CohortConfig(
        accession="GSE146771",
        title="Single-Cell Analyses Inform Mechanisms of Myeloid-Targeted therapies in colon cancer",
        indication="CRC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE146771_CRC.Leukocyte.10x.Metadata.txt.gz"),
        has_response_labels=True,
    ),
    "GSE164522": Tier1CohortConfig(
        accession="GSE164522",
        title="Single-cell analyses reveal phenotypic linkage between colorectal cancer and liver metastasis",
        indication="CRC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD45+ Immune-enriched)",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE164522_CRLM_LN_expression.csv.gz"),
        has_response_labels=True,
    ),
    "GSE235917": Tier1CohortConfig(
        accession="GSE235917",
        title="First-line durvalumab and tremelimumab with chemotherapy in RAS-mutated metastatic colorectal cancer: a phase 1b/2 trial [scRNA-Seq]",
        indication="CRC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE235917_5prim_matrix.mtx.gz"),
        has_response_labels=True,
    ),
    "GSE278406": Tier1CohortConfig(
        accession="GSE278406",
        title="Phenotypic plasticity and increased tissue infiltration of TREM1+ mono-macrophages following radiotherapy in rectal cancer. [scRNA-Seq]",
        indication="CRC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE278406_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE274321": Tier1CohortConfig(
        accession="GSE274321",
        title="IKKa modulates colorectal cancer metastasis by preventing tight junction stabilization and collective cell migration [scRNA-seq]",
        indication="CRC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE274321_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE309346": Tier1CohortConfig(
        accession="GSE309346",
        title="PI3K and MAPK signaling nodes as divergent drivers of phenotypic plasticity in cancer-associated fibroblasts in colorectal cancer [scRNA-Seq]",
        indication="CRC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE309346_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "CELLxGENE_2554a654": Tier1CohortConfig(
        accession="CELLxGENE_2554a654",
        title="progressive_plasticity_during_crc_metastasis_tumor",
        indication="CRC",
        technology="10x Chromium 3' v3/v3.1",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.CELLXGENE_H5AD,
        tar_filename=Nothing,
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "CELLxGENE_5ee552f5": Tier1CohortConfig(
        accession="CELLxGENE_5ee552f5",
        title="progressive_plasticity_during_crc_metastasis_non-tumor_epithelial",
        indication="CRC",
        technology="10x Chromium 3' v3/v3.1",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.CELLXGENE_H5AD,
        tar_filename=Nothing,
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "CELLxGENE_4b5afdf9": Tier1CohortConfig(
        accession="CELLxGENE_4b5afdf9",
        title="progressive_plasticity_during_crc_metastasis_untreated_epithelial",
        indication="CRC",
        technology="10x Chromium 3' v3/v3.1",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.CELLXGENE_H5AD,
        tar_filename=Nothing,
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE188711": Tier1CohortConfig(
        accession="GSE188711",
        title="Resolving the Difference Between Left-sided and Right-sided Colorectal Cancer by Single-cell Sequencing",
        indication="CRC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE188711_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE336564": Tier1CohortConfig(
        accession="GSE336564",
        title="ZFP36L2 orchestrates stress-adaptive plasticity in intestinal regeneration and colorectal cancer metastasis (PDO scRNAseq)",
        indication="CRC",
        technology="10x Chromium (3' unspecified)",
        cell_selection="FACS-sorted (Viability DAPI-/7-AAD- only)",
        archetype=Archetype.CELLXGENE_H5AD,
        tar_filename=Nothing,
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE216534": Tier1CohortConfig(
        accession="GSE216534",
        title="γδ T cells are effectors of immunotherapy in cancers with HLA class I defects",
        indication="CRC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_H5_DIR,
        tar_filename=Nothing,
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "CELLxGENE_387acac5": Tier1CohortConfig(
        accession="CELLxGENE_387acac5",
        title="progressive_plasticity_during_crc_metastasis_kg150_tumor",
        indication="CRC",
        technology="10x Chromium 3' v3/v3.1",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.CELLXGENE_H5AD,
        tar_filename=Nothing,
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "CELLxGENE_ef0d813e": Tier1CohortConfig(
        accession="CELLxGENE_ef0d813e",
        title="progressive_plasticity_during_crc_metastasis_kg183_tumor",
        indication="CRC",
        technology="10x Chromium 3' v3/v3.1",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.CELLXGENE_H5AD,
        tar_filename=Nothing,
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "CELLxGENE_2e95d453": Tier1CohortConfig(
        accession="CELLxGENE_2e95d453",
        title="progressive_plasticity_during_crc_metastasis_kg146_tumor",
        indication="CRC",
        technology="10x Chromium 3' v3/v3.1",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.CELLXGENE_H5AD,
        tar_filename=Nothing,
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE239676": Tier1CohortConfig(
        accession="GSE239676",
        title="Atlas of metastatic gastric cancer links ferroptosis to disease progression and immunotherapy response",
        indication="Gastric",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE239676_count_matrix.mtx.gz"),
        has_response_labels=True,
    ),
    "GSE319709": Tier1CohortConfig(
        accession="GSE319709",
        title="Single-cell transcriptomics reveals etiology-specific T-cell heterogeneity in hepatocellular carcinoma and implicates regulatory",
        indication="HCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE319709_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE272347": Tier1CohortConfig(
        accession="GSE272347",
        title="Late-stage tertiary lymphoid structures in hepatocellular carcinoma treated with neoadjuvant immune checkpoint blockade [scRNA-seq]",
        indication="HCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE272347_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE318418": Tier1CohortConfig(
        accession="GSE318418",
        title="Immunosuppressive monocytes are enriched in hepatocellular carcinoma patients with liver dysfunction in a phase II trial of combination sorafenib and nivolumab [scRNA-Seq]",
        indication="HCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE318418_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE299340": Tier1CohortConfig(
        accession="GSE299340",
        title="Single-cell RNA sequencing reveals B cell-related immunosuppressive landscape and a potential suppressor in hepatocellular carcinoma",
        indication="HCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE299340_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE282343": Tier1CohortConfig(
        accession="GSE282343",
        title="Viral-Track integrated single-cell RNA-sequencing reveals HBV lymphotropism and immunosuppressive microenvironment in HBV-associated hepatocellular carcinoma",
        indication="HCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE282343_ScRNA.counts.txt.gz"),
        has_response_labels=True,
    ),
    "GSE255830": Tier1CohortConfig(
        accession="GSE255830",
        title="Gene expression and T cell repertoire profile at single cell level of peripheral blood T cells after treatment with a personalized neoantigen vaccine (GNOS-PV02) and Pembrolizumab for advanced hepatocellular carcinoma.",
        indication="HCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE255830_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE281110": Tier1CohortConfig(
        accession="GSE281110",
        title="Molecular landscape of tumor-associated tissue-resident memory T cells in tumor microenvironment of hepatocellular carcinoma [HCC_scRNA]",
        indication="HCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE281110_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE233405": Tier1CohortConfig(
        accession="GSE233405",
        title="Immunohistochemical scoring of LAG-3 in conjunction with CD8 in  the tumor microenvironment predicts response to immunotherapy in  hepatocellular carcinoma",
        indication="HCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE233405_processed_scRNAseq_data.csv.gz"),
        has_response_labels=True,
    ),
    "GSE278324": Tier1CohortConfig(
        accession="GSE278324",
        title="Gene regulatory network analysis on snRNAseq revealed key regulators for hepatocellular carcinoma progression",
        indication="HCC",
        technology="Single-Nucleus RNA-seq",
        cell_selection="Nuclei Isolation (snRNA-seq)",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE278324_HCC.combined.count.matrix.txt.gz"),
        has_response_labels=True,
    ),
    "GSE265770": Tier1CohortConfig(
        accession="GSE265770",
        title="Gene expression profile at single cell level of CD56+ natural killer cells and CD8+ T cells from blood, spleen and HCC-PDX in humanized mice",
        indication="HCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE265770_NKsubset_integrated_metadata.csv.gz"),
        has_response_labels=True,
    ),
    "GSE318420": Tier1CohortConfig(
        accession="GSE318420",
        title="Immunosuppressive monocytes are enriched in hepatocellular carcinoma patients with liver dysfunction in a phase II trial of combination sorafenib and nivolumab",
        indication="HCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE318420_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE272348": Tier1CohortConfig(
        accession="GSE272348",
        title="Late-stage tertiary lymphoid structures in hepatocellular carcinoma treated with neoadjuvant immune checkpoint blockade",
        indication="HCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE272348_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE224411": Tier1CohortConfig(
        accession="GSE224411",
        title="Uncovering the spatial landscape of molecular interactions within the tumor microenvironment through latent spaces",
        indication="HCC",
        technology="10x Visium Spatial Transcriptomics",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.SPATIAL_ASSAY,
        tar_filename=Some("GSE224411_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE215428": Tier1CohortConfig(
        accession="GSE215428",
        title="Single-cell RNA sequencing of immune landscape in hepatocellular carcinoma treated with sintilimab and sorafenib",
        indication="HCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE215428_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE296954": Tier1CohortConfig(
        accession="GSE296954",
        title="Differentiation of tumor-infiltrating GZMK+ effector memory T cells associates with response to neoadjuvant immunotherapy in head and neck cancer",
        indication="HNSCC",
        technology="10x Chromium 5' (Immune Profiling)",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE296954_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE301720": Tier1CohortConfig(
        accession="GSE301720",
        title="Integrated single-cell and spatial analysis identifies context-dependent myeloid-T cell interactions in head and neck cancer immune checkpoint blockade response",
        indication="HNSCC",
        technology="10x Visium Spatial Transcriptomics",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.SPATIAL_ASSAY,
        tar_filename=Some("GSE301720_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE247582": Tier1CohortConfig(
        accession="GSE247582",
        title="Single-cell analysis of CD4+ cytotoxic T lymphocytes in human oral squamous cell carcinoma",
        indication="HNSCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE247582_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE296867": Tier1CohortConfig(
        accession="GSE296867",
        title="Post-treatment peripheral-blood T cells from HNSCC patients undergoing neoadjuvant immunotherapy treatment.",
        indication="HNSCC",
        technology="10x Chromium 5' (Immune Profiling)",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE296867_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE327189": Tier1CohortConfig(
        accession="GSE327189",
        title="Viral-based individualized neoantigen vaccine as adjuvant treatment in resected head and neck squamous cell carcinoma: immunogenicity and efficacy from a randomized Phase I trial",
        indication="HNSCC",
        technology="10x Chromium 5' (Immune Profiling)",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE327189_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE286410": Tier1CohortConfig(
        accession="GSE286410",
        title="Multi-modal Omics Analysis of a Paediatric Melanoma Highlights Mechanisms Underlying Treatment Resistance [Seq]",
        indication="Melanoma",
        technology="Single-Nucleus RNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE286410_Annotation.txt.gz"),
        has_response_labels=True,
    ),
    "GSE198265": Tier1CohortConfig(
        accession="GSE198265",
        title="Neoantigen specific CD4+ T cells in human melanoma have diverse differentiation states and correlate with CD8+ T cell, macrophage, and B cell function",
        indication="Melanoma",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE198265_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE242477": Tier1CohortConfig(
        accession="GSE242477",
        title="Single-cell profiling of acral melanoma infiltrating lymphocytes reveals a suppressive tumor microenvironment",
        indication="Melanoma",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE242477_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE303948": Tier1CohortConfig(
        accession="GSE303948",
        title="Genome accessibility profile at single cell level for conventional dendritic cells in human metastatic melanoma samples [snATAC-Seq]",
        indication="Melanoma",
        technology="Single-Nucleus RNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE303948_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE294273": Tier1CohortConfig(
        accession="GSE294273",
        title="Tumor-resident T cells and dendritic cells form an in situ archetype for immunotherapy response in melanoma [scRNA-seq]",
        indication="Melanoma",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE294273_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE320040": Tier1CohortConfig(
        accession="GSE320040",
        title="High-resolution and noninvasive profiling of the tumor microenvironment with spatial ecotypes [scRNA-seq]",
        indication="Melanoma",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE320040_CRC_tumor_scrna_counts.mtx.gz"),
        has_response_labels=True,
    ),
    "GSE210963": Tier1CohortConfig(
        accession="GSE210963",
        title="Single-Cell RNA-Seq Analysis of Patient Myeloid-Derived Suppressor Cells and the Response to Inhibition of Bruton’s Tyrosine Kinase",
        indication="Melanoma",
        technology="10x Chromium (3' unspecified)",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE210963_counts.csv.gz"),
        has_response_labels=True,
    ),
    "GSE256291": Tier1CohortConfig(
        accession="GSE256291",
        title="Comparing transcriptional profiles of CD14+  monocytes from melanoma patients and healthy donors by scRNA-seq.",
        indication="Melanoma",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE256291_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE244983": Tier1CohortConfig(
        accession="GSE244983",
        title="Molecular patterns of resistance to immune checkpoint blockade in melanoma [scRNA-Seq]",
        indication="Melanoma",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE244983_NormalizedData_scRNAseq.txt.gz"),
        has_response_labels=True,
    ),
    "GSE300446": Tier1CohortConfig(
        accession="GSE300446",
        title="Spatial tumour-immune ecosystems shape the efficacy of anti-PD1 immunotherapy in primary cutaneous melanoma [snRNAseq]",
        indication="Melanoma",
        technology="Single-Nucleus RNA-seq",
        cell_selection="Nuclei Isolation (snRNA-seq)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE300446_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE270464": Tier1CohortConfig(
        accession="GSE270464",
        title="Specific oncogene activation of the cell of origin in mucosal melanoma",
        indication="Melanoma",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE270464_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE211068": Tier1CohortConfig(
        accession="GSE211068",
        title="Integrated multiomics profiling identifies the differentiation program of regulatory T cells in human tumors [scRNA-seq]",
        indication="Melanoma",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE211068_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE241934": Tier1CohortConfig(
        accession="GSE241934",
        title="Neoadjuvant sintilimab plus chemotherapy in early-stage EGFR-mutant NSCLC: phase 2 trial interim results (NEOTIDE/CTONG2104)",
        indication="NSCLC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE241934_IIT_Matrix.mtx.gz"),
        has_response_labels=True,
    ),
    "GSE270148": Tier1CohortConfig(
        accession="GSE270148",
        title="Myeloid progenitor dysregulation fuels immunosuppressive macrophages in tumours",
        indication="NSCLC",
        technology="High-Throughput scRNA-seq",
        cell_selection="MACS Bead-selected",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE270148_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE303762": Tier1CohortConfig(
        accession="GSE303762",
        title="Benchmarking long-read RNA-sequencing technologies with LongBench: a cross-platform reference dataset profiling cancer cell lines with bulk and single-cell approaches",
        indication="NSCLC",
        technology="Single-Nucleus RNA-seq",
        cell_selection="Nuclei Isolation (snRNA-seq)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE303762_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE205049": Tier1CohortConfig(
        accession="GSE205049",
        title="Spatially Resolved Multi-Omics Single-Cell Analyses Inform Mechanisms of Immune Dysfunction in Pancreatic Cancer",
        indication="NSCLC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE205049_scRNA-seq-annotations.csv.gz"),
        has_response_labels=True,
    ),
    "GSE205354": Tier1CohortConfig(
        accession="GSE205354",
        title="Spatially resolved multi-omics single-cell analyses inform mechanisms of immune-dysfunction in pancreatic cancer",
        indication="NSCLC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE205354_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE233203": Tier1CohortConfig(
        accession="GSE233203",
        title="The single-cell level molecular characteristics according to combination immunotherapy response of non-small cell lung cancer patients.",
        indication="NSCLC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE233203_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE253718": Tier1CohortConfig(
        accession="GSE253718",
        title="Single-cell transcriptome reveals drug-resistance signature and immunosuppressive microenvironment in lung adenocarcinoma harboring EGFR mutation",
        indication="NSCLC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE253718_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE223779": Tier1CohortConfig(
        accession="GSE223779",
        title="Single-cell transcriptomic analysis uncovers intratumoral heterogeneity and drug-tolerant persister in ALK-rearranged lung adenocarcinoma",
        indication="NSCLC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE223779_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE307811": Tier1CohortConfig(
        accession="GSE307811",
        title="Development of antibody-drug conjugates targeting L1CAM to treat metastatic cancer",
        indication="NSCLC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE307811_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE285888": Tier1CohortConfig(
        accession="GSE285888",
        title="Single-Cell RNA Sequencing of Baseline PBMCs Predicts ICI efficacy and irAE Severity in NSCLC Patients",
        indication="NSCLC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE285888_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE276139": Tier1CohortConfig(
        accession="GSE276139",
        title="Gene expression profile at single cell level of cerebrospinal fluid (CSF) cells from lung adenocarcinoma leptomeningeal metastases patients (LUAD LM)",
        indication="NSCLC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE276139_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE279781": Tier1CohortConfig(
        accession="GSE279781",
        title="CD137 agonism enhances anti-PD1 induced activation of clonally expanded CD8+ T cells in a neoadjuvant pancreatic cancer clinical trial",
        indication="PDAC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE279781_matrix.mtx.gz"),
        has_response_labels=True,
    ),
    "GSE212966": Tier1CohortConfig(
        accession="GSE212966",
        title="Single-cell RNA-seq reveals immune landscape of pancreatic cancer",
        indication="PDAC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE212966_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE283206": Tier1CohortConfig(
        accession="GSE283206",
        title="Combined Flt3L and CD40 agonism restores dendritic cell driven T cell immunity in mouse models and patients with pancreatic cancer [human]",
        indication="PDAC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE283206_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE311788": Tier1CohortConfig(
        accession="GSE311788",
        title="DeCAF redefines fibroblast states uncovering multidimensional tumor-stroma relationships driving clinical tumor progression and immunotherapy response [scRNA-Seq]",
        indication="PDAC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE311788_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE211644": Tier1CohortConfig(
        accession="GSE211644",
        title="Single cell transcriptomic and T cell repertoire analysis reveals trajectory of tumor - infiltrating lymphocyte states in pancreatic cancer",
        indication="PDAC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE211644_fresh_matrix.mtx.gz"),
        has_response_labels=True,
    ),
    "GSE348275": Tier1CohortConfig(
        accession="GSE348275",
        title="Clonal lineage tracing and parallel multiomics profiling reveal transcriptional heterogeneity induced by ARID1A deficiency",
        indication="PDAC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE348275_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE318413": Tier1CohortConfig(
        accession="GSE318413",
        title="Patient-derived orthotopic xenograft models recapitulate the peritoneal dissemination of pancreatic cancer and delineate its transcriptional and regulatory programs [scRNA-seq]",
        indication="PDAC",
        technology="Single-Nucleus RNA-seq",
        cell_selection="Nuclei Isolation (snRNA-seq)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE318413_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE156405": Tier1CohortConfig(
        accession="GSE156405",
        title="Elucidation of tumor-stromal heterogeneity and the ligand-receptor interactome by single cell transcriptomics in real-world pancreatic cancer biopsies",
        indication="PDAC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE156405_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE335452": Tier1CohortConfig(
        accession="GSE335452",
        title="Single-cell RNA-seq and spatial transcriptomics characterize CD8+ exhausted T cells in pancreatic ductal adenocarcinoma",
        indication="PDAC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE335452_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE285701": Tier1CohortConfig(
        accession="GSE285701",
        title="The paradoxical significance of CD39+CD8+ T cells in clear cell renal cell carcinoma [scRNA-seq]",
        indication="ccRCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Nothing,
        matrix_filename=Some("GSE285701_single_cell_adt_data.txt.gz"),
        has_response_labels=True,
    ),
    "GSE304466": Tier1CohortConfig(
        accession="GSE304466",
        title="Single-cell transcriptome combined with spatial transcriptome to investigate the molecular mechanism associated with autophagy of clear renal cell carcinoma",
        indication="ccRCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE304466_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE223808": Tier1CohortConfig(
        accession="GSE223808",
        title="Exhausted intratumoral Vδ2- γδ T cells in human kidney cancer retain effector function [scRNA-seq]",
        indication="ccRCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD3+/CD8+ T-cell enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE223808_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE220313": Tier1CohortConfig(
        accession="GSE220313",
        title="Gene expression profile and TCR sequencing data of CD45+ cells sorted from single cell suspensions of tumors after in vitro anti-CD3 stimulation and treatments",
        indication="ccRCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="FACS-sorted (CD45+ Immune-enriched)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE220313_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE254498": Tier1CohortConfig(
        accession="GSE254498",
        title="Integrating whole-exome sequencing and scRNA-seq reveal the characteristic in one clear cell renal cell carcinoma sample arising in the setting of VHL disease",
        indication="ccRCC",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE254498_RAW.tar"),
        matrix_filename=Nothing,
        has_response_labels=True,
    ),
    "GSE222315": Tier1CohortConfig(
        accession="GSE222315",
        title="Single-cell transcriptome profiling of human bladder cancer and normal adjacent tissues",
        indication="Bladder",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.FLAT_COUNT_TABLE,
        tar_filename=Some("GSE222315_RAW.tar"),
        has_response_labels=False,
    ),
    "GSE302781": Tier1CohortConfig(
        accession="GSE302781",
        title="snRNASeq of Metastatic Urothelial carcinoma samples",
        indication="Bladder",
        technology="Single-Nucleus RNA-seq",
        cell_selection="Nuclei Isolation (snRNA-seq)",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE302781_RAW.tar"),
        patient_col=Some("cell line"),
        sample_col=Some("Sample_title"),
        has_response_labels=False,
    ),
    "GSE326225": Tier1CohortConfig(
        accession="GSE326225",
        title="A Single-Cell Atlas of Muscle-Invasive Bladder Cancer Reveals Lineage-Specific Vulnerabilities",
        indication="Bladder",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE326225_RAW.tar"),
        patient_col=Some("Sample_title"),
        sample_col=Some("Sample_title"),
        has_response_labels=False,
    ),
    "GSE301651": Tier1CohortConfig(
        accession="GSE301651",
        title="Single cell sequencing of bladder cancer patients",
        indication="Bladder",
        technology="High-Throughput scRNA-seq",
        cell_selection="Unselected / Total Single-Cell Suspension",
        archetype=Archetype.GEO_10X_TAR,
        tar_filename=Some("GSE301651_RAW.tar"),
        has_response_labels=False,
    ),
}


# =========================================================================
# Clinical Schema Harmonization
# =========================================================================
def _standardize_tier1_obs(
    adata: ad.AnnData,
    config: Tier1CohortConfig,
    meta_df: pd.DataFrame | None = None,
) -> ad.AnnData:
    """Standardize .obs clinical and technical metadata purely."""
    obs = adata.obs.copy()
    obs["cell_id"] = adata.obs_names.astype(str)
    obs["dataset_id"] = config.accession
    obs["indication"] = config.indication
    obs["cancer_type"] = config.indication
    obs["technology"] = config.technology
    obs["cell_selection_strategy"] = config.cell_selection

    # Merge metadata from GEO series matrix if provided
    if meta_df is not None and not meta_df.empty:
        sample_keys = obs["sample_id"].astype(str) if "sample_id" in obs.columns else obs.index.astype(str)
        gsm_matches: dict[str, str] = {}
        for s in set(sample_keys):
            if s in meta_df.index:
                gsm_matches[s] = s
            else:
                if "Sample_title" in meta_df.columns:
                    title_match = meta_df[meta_df["Sample_title"].astype(str) == s]
                    if not title_match.empty:
                        gsm_matches[s] = str(title_match.index[0])
                        continue
                for gsm_idx in meta_df.index:
                    if str(gsm_idx) in s:
                        gsm_matches[s] = str(gsm_idx)
                        break

        if gsm_matches:
            matched_gsm_series = sample_keys.map(gsm_matches)
            for col in meta_df.columns:
                if col not in obs.columns:
                    val_series = matched_gsm_series.map(meta_df[col].to_dict())
                    if val_series.notna().any():
                        obs[col] = val_series.fillna("NA")

    # Map patient/donor
    patient_col = config.patient_col.value_or(None)
    if patient_col and patient_col in obs.columns:
        obs["patient_id"] = obs[patient_col].astype(str)
    elif "patient_id" not in obs.columns:
        for c in ("PMID_donor_id", "donor_id", "patient", "patientID", "donor", "subject", "sample"):
            if c in obs.columns:
                obs["patient_id"] = obs[c].astype(str)
                break
        else:
            obs["patient_id"] = f"{config.accession}_donor_unspecified"

    # Map sample
    sample_col = config.sample_col.value_or(None)
    if sample_col and sample_col in obs.columns:
        obs["sample_id"] = obs[sample_col].astype(str)
    elif "sample_id" not in obs.columns:
        obs["sample_id"] = obs["patient_id"]

    if config.accession == "GSE301651" and "Sample_title" in obs.columns:
        # Sample_title is e.g. "YCK,bladder cancer patient" or "BC-TWM,bladder cancer patient"
        sample_codes = obs["Sample_title"].astype(str).str.split(",").str[0].str.strip()
        obs["sample_id"] = sample_codes
        patient_map = {
            "BC-TWM": "TWM",
            "BC-LNM-TWM": "TWM",
            "PB-TWM": "TWM",
        }
        obs["patient_id"] = [patient_map.get(sc, sc) for sc in sample_codes]
        if "tissue" in obs.columns:
            obs["specimen_type"] = obs["tissue"].astype(str)

    # Map response
    response_col = config.response_col.value_or(None)
    raw_response_series = None
    if response_col and response_col in obs.columns:
        raw_response_series = obs[response_col]
    elif "clinical_response_raw" in obs.columns:
        raw_response_series = obs["clinical_response_raw"]
    elif "clinical_response" in obs.columns:
        raw_response_series = obs["clinical_response"]
    else:
        for c in (
            "response", "Response", "RECIST", "benefit", "responder", "outcome",
            "Combined_outcome", "pCR_status", "Pathologic_Response",
            "pathological_response", "Path_response",
        ):
            if c in obs.columns:
                raw_response_series = obs[c]
                break

    if raw_response_series is not None:
        obs["clinical_response_raw"] = raw_response_series.astype(str)
        response_map = config.response_map.value_or(None)
        if response_map:
            obs["clinical_response"] = obs["clinical_response_raw"].map(response_map).fillna("not-evaluable")
        else:
            harmonized = [_harmonize_single_response(v) for v in raw_response_series]
            obs["clinical_response"] = [h[0] for h in harmonized]
            obs["clinical_response_raw"] = [h[1] for h in harmonized]
    else:
        obs["clinical_response"] = "not-evaluable"
        obs["clinical_response_raw"] = "Unspecified"

    # Map treatment
    treatment_col = config.treatment_col.value_or(None)
    if treatment_col and treatment_col in obs.columns:
        obs["treatment_status"] = obs[treatment_col].astype(str)
    elif "treatment_status" not in obs.columns:
        for c in ("treatment", "treatment_status", "timepoint", "Timepoint", "cycle"):
            if c in obs.columns:
                obs["treatment_status"] = obs[c].astype(str)
                break
        else:
            obs["treatment_status"] = "Pre-treatment" if config.has_response_labels else "Baseline"

    # Map cell type
    cell_type_col = config.cell_type_col.value_or(None)
    if cell_type_col and cell_type_col in obs.columns:
        obs["cell_type"] = obs[cell_type_col].astype(str)
    elif "cell_type" not in obs.columns:
        for c in ("cell_type", "celltype", "CellType", "major_cell_type", "cluster"):
            if c in obs.columns:
                obs["cell_type"] = obs[c].astype(str)
                break
        else:
            obs["cell_type"] = "Unspecified"

    # Clean object columns to prevent HDF5 serialization crashes
    for c in obs.columns:
        if obs[c].dtype == "object":
            obs[c] = obs[c].fillna("NA").astype(str)

    new_adata = adata.copy()
    new_adata.obs = obs
    new_adata.obs_names = obs["cell_id"]
    new_adata.obs_names.name = "cell_id"
    new_adata.var_names_make_unique()
    return new_adata


# =========================================================================
# Archetype Ingestion Handlers
# =========================================================================
def _ingest_geo_10x_tar(config: Tier1CohortConfig, raw_dir: Path) -> Result[ad.AnnData, str]:
    """Ingest GEO RAW tarball containing 10x MTX triplets or sample .h5 files."""
    tar_stem = config.tar_filename.value_or(f"{config.accession}_RAW.tar")
    tar_path = raw_dir / tar_stem
    if not tar_path.exists():
        tars = list(raw_dir.glob("*.tar")) + list(raw_dir.glob("*.tar.gz"))
        if tars:
            tar_path = tars[0]
        else:
            return Failure(f"Tar archive for {config.accession} not found at {tar_path}")

    extract_dir = raw_dir / f"unpacked_{config.accession.lower()}"
    ext_res = _extract_tar_safely(tar_path, extract_dir)
    match ext_res:
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    # Unpack any sub-tarballs (e.g. GSM*.tar.gz) inside extract_dir
    for st in list(extract_dir.glob("*.tar.gz")) + list(extract_dir.glob("*.tar")):
        sub_dir = extract_dir / st.name.replace(".tar.gz", "").replace(".tar", "")
        if not sub_dir.exists():
            _extract_tar_safely(st, sub_dir)

    load_res = _load_10x_triplets_or_h5(
        extract_dir,
        min_counts=config.min_counts,
        min_genes=config.min_genes,
    )
    match load_res:
        case Failure(_):
            return _ingest_flat_count_table(config, extract_dir)
        case Success(adatas):
            combined = adatas[0] if len(adatas) == 1 else ad.concat(adatas, join="outer")
            combined.obs_names_make_unique()

    meta_df = None
    meta_res = _parse_geo_series_matrix(raw_dir, config.accession)
    if isinstance(meta_res, Success):
        meta_df = meta_res.unwrap()

    combined = _standardize_tier1_obs(combined, config, meta_df=meta_df)
    combined = tag_expression_metadata(combined)
    return Success(combined)


def _ingest_geo_10x_h5_dir(config: Tier1CohortConfig, raw_dir: Path) -> Result[ad.AnnData, str]:
    """Ingest directory of 10x HDF5 (.h5 / .h5.gz) feature matrices."""
    for hg in raw_dir.glob("*.h5.gz"):
        decomp = hg.with_suffix("")
        if not decomp.exists():
            with gzip.open(hg, "rb") as f_in, open(decomp, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

    h5_files = sorted(list(raw_dir.glob("*.h5"))) or sorted(list(raw_dir.glob("**/*.h5")))
    if not h5_files:
        return Failure(f"No 10x .h5 files found in {raw_dir}")

    adatas: list[ad.AnnData] = []
    for hf in h5_files:
        sample_name = hf.name.replace(".h5", "").replace("filtered_feature_bc_matrix_", "")
        try:
            sub_a = sc.read_10x_h5(hf)
            sub_a.obs_names = [f"{sample_name}_{b}" for b in sub_a.obs_names]
            sub_a.var_names_make_unique()
            sub_a.obs["sample_id"] = sample_name
            sub_a.obs["patient_id"] = sample_name.split("_")[0]

            if sp.issparse(sub_a.X):
                counts = np.asarray(sub_a.X.sum(axis=1)).ravel()
                genes = np.asarray((sub_a.X > 0).sum(axis=1)).ravel()
            else:
                counts = np.asarray(np.sum(sub_a.X, axis=1)).ravel()
                genes = np.asarray(np.sum(sub_a.X > 0, axis=1)).ravel()

            valid = (counts >= config.min_counts) & (genes >= config.min_genes)
            if np.any(valid):
                adatas.append(sub_a[valid].copy())
        except Exception as exc:
            logger.warning("Failed to parse h5 file %s: %s", hf.name, exc)

    if not adatas:
        return Failure(f"Failed to load any valid cells from h5 files in {raw_dir}")

    combined = adatas[0] if len(adatas) == 1 else ad.concat(adatas, join="outer")
    combined.obs_names_make_unique()

    meta_df = None
    meta_res = _parse_geo_series_matrix(raw_dir, config.accession)
    if isinstance(meta_res, Success):
        meta_df = meta_res.unwrap()

    combined = _standardize_tier1_obs(combined, config, meta_df=meta_df)
    combined = tag_expression_metadata(combined)
    return Success(combined)


def _ingest_cellxgene_h5ad(config: Tier1CohortConfig, raw_dir: Path) -> Result[ad.AnnData, str]:
    """Ingest CZ CELLxGENE or curated AnnData H5AD file restoring raw counts."""
    h5ad_files = list(raw_dir.glob("*.h5ad"))
    if not h5ad_files:
        return Failure(f"H5AD file not found for {config.accession} in {raw_dir}")

    try:
        adata = ad.read_h5ad(h5ad_files[0])
        adata = _restore_raw_anndata(adata)
        adata = _standardize_tier1_obs(adata, config)
        adata = tag_expression_metadata(adata)
        return Success(adata)
    except Exception as exc:
        return Failure(f"Failed to load H5AD {h5ad_files[0].name}: {exc}")


def _ingest_flat_count_table(config: Tier1CohortConfig, raw_dir: Path) -> Result[ad.AnnData, str]:
    """Ingest flat count tables (CSV, TSV, TXT, or standalone MTX)."""
    # 0. Unpack any tar archives if present and not yet extracted
    extract_dir = raw_dir / f"unpacked_{config.accession.lower()}"
    for tf in list(raw_dir.glob("*.tar")) + list(raw_dir.glob("*.tar.gz")):
        if not extract_dir.exists():
            _extract_tar_safely(tf, extract_dir)

    search_dirs = [extract_dir, raw_dir] if extract_dir.exists() else [raw_dir]

    # 1. First, check if 10x triplets exist in search dirs
    for d in search_dirs:
        load_10x_res = _load_10x_triplets_or_h5(d, min_counts=config.min_counts, min_genes=config.min_genes)
        if isinstance(load_10x_res, Success):
            adatas = load_10x_res.unwrap()
            combined = adatas[0] if len(adatas) == 1 else ad.concat(adatas, join="outer")
            combined.obs_names_make_unique()
            meta_res = _parse_geo_series_matrix(raw_dir, config.accession)
            meta_df = meta_res.unwrap() if isinstance(meta_res, Success) else None
            combined = _standardize_tier1_obs(combined, config, meta_df=meta_df)
            combined = tag_expression_metadata(combined)
            return Success(combined)

    # 2. Check for standalone MTX files across search_dirs
    mtx_candidates: list[Path] = []
    for d in search_dirs:
        mtx_candidates.extend(d.glob("*.mtx.gz"))
        mtx_candidates.extend(d.glob("*.mtx"))
        mtx_candidates.extend(d.glob("**/*.mtx.gz"))
        mtx_candidates.extend(d.glob("**/*.mtx"))
    mtx_files = sorted(list(dict.fromkeys(mtx_candidates)))

    if mtx_files:
        mtx_path = mtx_files[0]
        try:
            mat = sio.mmread(mtx_path)
            mat_csr = _to_csr(mat)

            # Discover barcodes and features in directory
            bc_candidates: list[Path] = []
            feat_candidates: list[Path] = []
            for d in search_dirs:
                bc_candidates.extend(d.glob("*barcode*"))
                bc_candidates.extend(d.glob("*cell*"))
                feat_candidates.extend(d.glob("*feature*"))
                feat_candidates.extend(d.glob("*gene*"))

            bcs: list[str] = []
            gene_names: list[str] = []
            if bc_candidates:
                try:
                    bcs = pd.read_csv(bc_candidates[0], header=None, sep=None, engine="python")[0].astype(str).tolist()
                except Exception:
                    pass
            if feat_candidates:
                try:
                    f_df = pd.read_csv(feat_candidates[0], header=None, sep=None, engine="python")
                    gene_names = (f_df[1] if len(f_df.columns) > 1 else f_df[0]).astype(str).tolist()
                except Exception:
                    pass

            # Handle shape orientation (if 10x MTX, shape is genes x cells)
            if bcs and len(bcs) == mat_csr.shape[1] and mat_csr.shape[0] != len(bcs):
                mat_csr = mat_csr.T

            n_cells = mat_csr.shape[0]
            n_genes = mat_csr.shape[1]

            if not bcs or len(bcs) != n_cells:
                bcs = [f"{config.accession}_cell_{i}" for i in range(n_cells)]
            else:
                bcs = [f"{config.accession}_{b}" for b in bcs]

            if not gene_names or len(gene_names) != n_genes:
                gene_names = [f"gene_{j}" for j in range(n_genes)]

            adata = ad.AnnData(
                X=mat_csr,
                obs=pd.DataFrame(index=bcs),
                var=pd.DataFrame(index=gene_names),
            )
            adata.var_names_make_unique()
            adata.obs_names_make_unique()
            meta_res = _parse_geo_series_matrix(raw_dir, config.accession)
            meta_df = meta_res.unwrap() if isinstance(meta_res, Success) else None
            adata = _standardize_tier1_obs(adata, config, meta_df=meta_df)
            adata = tag_expression_metadata(adata)
            return Success(adata)
        except Exception as exc:
            logger.warning("MTX parse failed for %s: %s", mtx_path.name, exc)

    # 3. Check for CSV / TSV / TXT count tables across search dirs
    candidates: list[Path] = []
    target_pattern = config.matrix_filename.value_or(None)
    for d in search_dirs:
        if target_pattern:
            candidates.extend(d.glob(f"*{target_pattern}*"))
            candidates.extend(d.glob(f"**/*{target_pattern}*"))
        for ext in ("*.csv.gz", "*.tsv.gz", "*.txt.gz", "*.csv", "*.tsv", "*.txt"):
            candidates.extend(d.glob(ext))
            candidates.extend(d.glob(f"**/{ext}"))

    candidates = list(dict.fromkeys(candidates))
    if target_pattern:
        target_matches = [c for c in candidates if target_pattern in c.name]
        if target_matches:
            candidates = target_matches

    # If not explicitly matched via target_pattern, prioritize expression tables over metadata
    if not (target_pattern and any(target_pattern in c.name for c in candidates)):
        expr_candidates = [
            c for c in candidates
            if not any(x in c.name.lower() for x in (
                "meta", "filelist", "series_matrix", "pheno", "barcode", "feature", "gene",
                "sample", "readme", "annotation"
            ))
            and not c.name.endswith(".mtx")
            and not c.name.endswith(".mtx.gz")
        ]
        candidates = expr_candidates if expr_candidates else [c for c in candidates if not c.name.endswith((".mtx", ".mtx.gz"))]

    if not candidates:
        return Failure(f"No count table candidates found for {config.accession} in {raw_dir}")

    mat_file = candidates[0]
    try:
        sep = "," if ".csv" in mat_file.name else ("\t" if ".tsv" in mat_file.name else None)
        df = pd.read_csv(mat_file, sep=sep, index_col=0, comment="#", engine="c" if sep else "python")

        try:
            if df.shape[0] > df.shape[1]:
                X = sp.csr_matrix(df.values.T, dtype=np.float32)
                obs_names = [f"{config.accession}_{c}" for c in df.columns.astype(str)]
                var_names = df.index.astype(str).tolist()
                obs_df = pd.DataFrame(index=obs_names)
            else:
                X = sp.csr_matrix(df.values, dtype=np.float32)
                obs_names = [f"{config.accession}_{c}" for c in df.index.astype(str)]
                var_names = df.columns.astype(str).tolist()
                obs_df = pd.DataFrame(index=obs_names)
        except (ValueError, TypeError):
            # Non-numeric table (e.g. metadata/annotation table uploaded instead of count matrix)
            logger.info("Count table %s is non-numeric; preserving as cell metadata.", mat_file.name)
            X = sp.csr_matrix((df.shape[0], 0), dtype=np.float32)
            obs_names = [f"{config.accession}_{c}" for c in df.index.astype(str)]
            var_names = []
            obs_df = df.copy()
            obs_df.index = obs_names

        adata = ad.AnnData(
            X=X,
            obs=obs_df,
            var=pd.DataFrame(index=var_names),
        )
        adata.var_names_make_unique()
        adata.obs_names_make_unique()

        meta_res = _parse_geo_series_matrix(raw_dir, config.accession)
        meta_df = meta_res.unwrap() if isinstance(meta_res, Success) else None
        adata = _standardize_tier1_obs(adata, config, meta_df=meta_df)
        adata = tag_expression_metadata(adata)
        return Success(adata)
    except Exception as exc:
        return Failure(f"Failed to read count table {mat_file.name} for {config.accession}: {exc}")


def _ingest_spatial_assay(config: Tier1CohortConfig, raw_dir: Path) -> Result[ad.AnnData, str]:
    """Ingest 10x Visium or CosMx spatial profiling data."""
    tars = list(raw_dir.glob("*.tar")) + list(raw_dir.glob("*.tar.gz"))
    extract_dir = raw_dir / f"unpacked_{config.accession.lower()}"
    if tars and not extract_dir.exists():
        ext_res = _extract_tar_safely(tars[0], extract_dir)
        match ext_res:
            case Failure(err):
                return Failure(err)
            case Success(_):
                pass

    search_dir = extract_dir if extract_dir.exists() else raw_dir

    load_res = _load_10x_triplets_or_h5(search_dir, min_counts=config.min_counts, min_genes=config.min_genes)
    match load_res:
        case Failure(_):
            flat_res = _ingest_flat_count_table(config, search_dir)
            match flat_res:
                case Failure(err):
                    return Failure(err)
                case Success(ad_flat):
                    combined = ad_flat
        case Success(adatas):
            combined = adatas[0] if len(adatas) == 1 else ad.concat(adatas, join="outer")
            combined.obs_names_make_unique()

    for pos_file in list(search_dir.glob("**/tissue_positions*.csv")):
        try:
            pos_df = pd.read_csv(pos_file, header=None)
            if len(pos_df.columns) >= 6:
                pos_df.set_index(0, inplace=True)
                coords = np.zeros((combined.n_obs, 2), dtype=np.float32)
                for i, bc in enumerate(combined.obs.index):
                    raw_bc = bc.split("_")[-1]
                    if raw_bc in pos_df.index:
                        coords[i] = [float(pos_df.loc[raw_bc, 4]), float(pos_df.loc[raw_bc, 5])]
                combined.obsm["spatial"] = coords
                break
        except Exception as exc:
            logger.warning("Could not parse spatial coordinates from %s: %s", pos_file.name, exc)

    meta_res = _parse_geo_series_matrix(raw_dir, config.accession)
    meta_df = meta_res.unwrap() if isinstance(meta_res, Success) else None
    combined = _standardize_tier1_obs(combined, config, meta_df=meta_df)
    combined = tag_expression_metadata(combined)
    return Success(combined)


def _ingest_gse222315_expression_tables(config: Tier1CohortConfig, raw_dir: Path) -> Result[ad.AnnData, str]:
    """Ingest GSE222315 paired MTX matrices and expression table headers across 13 samples."""
    tar_stem = config.tar_filename.value_or(f"{config.accession}_RAW.tar")
    tar_path = raw_dir / tar_stem
    if not tar_path.exists():
        tars = list(raw_dir.glob("*.tar")) + list(raw_dir.glob("*.tar.gz"))
        if tars:
            tar_path = tars[0]
        else:
            return Failure(f"Tar archive for {config.accession} not found at {tar_path}")

    extract_dir = raw_dir / f"unpacked_{config.accession.lower()}"
    ext_res = _extract_tar_safely(tar_path, extract_dir)
    match ext_res:
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    # Read reference gene list once from any sample expression file
    expr_files = sorted(list(extract_dir.glob("*_expression.*.gz")))
    if not expr_files:
        return Failure(f"No expression files (*_expression.*.gz) found in {extract_dir}")

    ref_exp = expr_files[0]
    genes: list[str] = []
    gene_names: list[str] = []
    try:
        with gzip.open(ref_exp, "rt", encoding="utf-8", errors="replace") as gz:
            gz.readline()  # Skip header
            for line in gz:
                parts = line.strip().split("\t", 2)
                if len(parts) >= 2:
                    genes.append(parts[0])
                    gene_names.append(parts[1])
    except Exception as exc:
        return Failure(f"Failed to read gene definitions from {ref_exp.name}: {exc}")

    if not genes:
        return Failure(f"No genes parsed from reference file {ref_exp.name}")

    var_df = pd.DataFrame(index=genes)
    var_df["gene_name"] = gene_names

    # Ingest each of the 13 sample pairs
    adatas: list[ad.AnnData] = []
    mtx_files = sorted(list(extract_dir.glob("*_matrix.mtx.gz")))
    if not mtx_files:
        return Failure(f"No matrix files (*_matrix.mtx.gz) found in {extract_dir}")

    for mtx_p in mtx_files:
        prefix = mtx_p.name.replace("_matrix.mtx.gz", "")
        matching_exp = [f for f in expr_files if prefix in f.name]
        if not matching_exp:
            continue
        exp_p = matching_exp[0]

        try:
            with gzip.open(exp_p, "rt", encoding="utf-8", errors="replace") as gz:
                header_line = gz.readline().strip().split("\t")
                barcodes = header_line[2:]

            with gzip.open(mtx_p, "rb") as gz:
                mat = sio.mmread(gz)

            mat_csr = _to_csr(mat)
            if mat_csr.shape[0] != len(barcodes) and mat_csr.shape[1] == len(barcodes):
                mat_csr = mat_csr.T

            if mat_csr.shape[0] != len(barcodes) or mat_csr.shape[1] != len(genes):
                logger.warning(
                    "Dimension mismatch in sample %s: mat %s vs (%d cells, %d genes)",
                    prefix,
                    mat_csr.shape,
                    len(barcodes),
                    len(genes),
                )
                continue

            sample_id = prefix
            sub_a = ad.AnnData(
                X=mat_csr.astype(np.float32),
                obs=pd.DataFrame(index=[f"{sample_id}_{b}" for b in barcodes]),
                var=var_df.copy(),
            )
            sub_a.obs["sample_id"] = sample_id
            tokens = prefix.split("_")
            patient = tokens[1].upper() if len(tokens) > 1 else prefix
            sub_a.obs["patient_id"] = patient
            sub_a.obs["tissue_type"] = "tumor" if "_BCa" in prefix else "adjacent normal"
            adatas.append(sub_a)
        except Exception as exc:
            logger.warning("Failed to parse sample %s: %s", prefix, exc)

    if not adatas:
        return Failure(f"Failed to load any valid samples from {extract_dir}")

    combined = adatas[0] if len(adatas) == 1 else ad.concat(adatas, join="outer")
    combined.obs_names_make_unique()
    combined.var_names_make_unique()

    meta_df = None
    meta_res = _parse_geo_series_matrix(raw_dir, config.accession)
    if isinstance(meta_res, Success):
        meta_df = meta_res.unwrap()

    combined = _standardize_tier1_obs(combined, config, meta_df=meta_df)
    combined = tag_expression_metadata(combined)
    return Success(combined)


def _auto_download_tier1_assets(config: Tier1CohortConfig, raw_dir: Path) -> Result[None, str]:
    """Ensure raw vendor data files for a Tier 1 cohort are present on disk, downloading if missing."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    match config.archetype:
        case Archetype.GEO_10X_TAR:
            tar_stem = config.tar_filename.value_or(f"{config.accession}_RAW.tar")
            tar_path = raw_dir / tar_stem
            existing_tars = list(raw_dir.glob("*.tar")) + list(raw_dir.glob("*.tar.gz"))
            if not tar_path.exists() and not existing_tars:
                logger.info("Auto-downloading RAW tarball for %s to %s...", config.accession, raw_dir)
                match download_geo_supplementary(config.accession, raw_dir, expected_files=[tar_stem]):
                    case Failure(_):
                        match download_geo_supplementary(config.accession, raw_dir):
                            case Failure(err2):
                                return Failure(f"Failed to auto-download {config.accession} tarball: {err2}")
                            case Success(_):
                                pass
                    case Success(_):
                        pass

        case Archetype.GEO_10X_H5_DIR:
            h5_files = list(raw_dir.glob("*.h5")) + list(raw_dir.glob("*.h5.gz"))
            if not h5_files:
                logger.info("Auto-downloading 10x .h5 files for %s to %s...", config.accession, raw_dir)
                match download_geo_supplementary(config.accession, raw_dir):
                    case Failure(err):
                        return Failure(f"Failed to auto-download {config.accession} h5 files: {err}")
                    case Success(_):
                        pass

        case Archetype.FLAT_COUNT_TABLE:
            mat_files = (
                list(raw_dir.glob("*matrix.mtx*"))
                + list(raw_dir.glob("*.csv*"))
                + list(raw_dir.glob("*.tsv*"))
                + list(raw_dir.glob("*.txt*"))
                + list(raw_dir.glob("*.tar*"))
            )
            if not mat_files:
                logger.info("Auto-downloading count table for %s to %s...", config.accession, raw_dir)
                expected = [config.matrix_filename.unwrap()] if isinstance(config.matrix_filename, Some) else None
                match download_geo_supplementary(config.accession, raw_dir, expected_files=expected):
                    case Failure(_):
                        match download_geo_supplementary(config.accession, raw_dir):
                            case Failure(err2):
                                return Failure(f"Failed to auto-download {config.accession} count table: {err2}")
                            case Success(_):
                                pass
                    case Success(_):
                        pass

            # Extract tar if present and unpacked dir missing
            extract_dir = raw_dir / f"unpacked_{config.accession.lower()}"
            for tf in list(raw_dir.glob("*.tar")) + list(raw_dir.glob("*.tar.gz")):
                if not extract_dir.exists():
                    _extract_tar_safely(tf, extract_dir)

        case Archetype.SPATIAL_ASSAY:
            existing = list(raw_dir.glob("*.tar*")) + list(raw_dir.glob("*.h5*")) + list(raw_dir.glob("*.mtx*"))
            if not existing:
                logger.info("Auto-downloading spatial assay for %s to %s...", config.accession, raw_dir)
                match download_geo_supplementary(config.accession, raw_dir):
                    case Failure(err):
                        return Failure(f"Failed to auto-download spatial data for {config.accession}: {err}")
                    case Success(_):
                        pass

        case Archetype.CELLXGENE_H5AD:
            h5ad_files = list(raw_dir.glob("*.h5ad"))
            if not h5ad_files:
                logger.info("Auto-downloading CELLxGENE .h5ad for %s...", config.accession)
                repo_root = find_repo_root()
                reg_path = repo_root / "data/registry/discovered_solid_tumor_sc_datasets.parquet"
                if reg_path.exists():
                    df = pl.read_parquet(reg_path)
                    matched = df.filter(pl.col("accession") == config.accession)
                    if not matched.is_empty():
                        urls_str = matched["download_urls"][0]
                        urls = [u.strip() for u in urls_str.split(";") if u.strip()]
                        if urls:
                            dest_file = raw_dir / f"{config.accession}.h5ad"
                            match download_single_file(urls[0], dest_file):
                                case Failure(err):
                                    return Failure(f"Failed to download CELLxGENE {config.accession}: {err}")
                                case Success(_):
                                    pass
                if not list(raw_dir.glob("*.h5ad")):
                    return Failure(f"No .h5ad file found or downloaded for {config.accession}")

    return Success(None)


# =========================================================================
# Unified Dispatcher & Public API
# =========================================================================
def load_tier1_cohort(
    accession: str,
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Unified declarative loader for any registered Tier 1 single-cell cohort."""
    if accession not in TIER1_COHORTS:
        return Failure(f"Dataset '{accession}' is not a registered Tier 1 cohort.")

    config = TIER1_COHORTS[accession]
    logger.info("Loading Tier 1 cohort %s (%s, Archetype=%s)...", config.accession, config.indication, config.archetype.value)

    if auto_download:
        dl_res = _auto_download_tier1_assets(config, raw_dir)
        match dl_res:
            case Failure(err):
                logger.warning("Auto-download warning for %s: %s", config.accession, err)
            case Success(_):
                pass

    match config.archetype:
        case Archetype.GEO_10X_TAR:
            res = _ingest_geo_10x_tar(config, raw_dir)
        case Archetype.GEO_10X_H5_DIR:
            res = _ingest_geo_10x_h5_dir(config, raw_dir)
        case Archetype.CELLXGENE_H5AD:
            res = _ingest_cellxgene_h5ad(config, raw_dir)
        case Archetype.FLAT_COUNT_TABLE:
            if config.accession == "GSE222315":
                res = _ingest_gse222315_expression_tables(config, raw_dir)
            else:
                res = _ingest_flat_count_table(config, raw_dir)
        case Archetype.SPATIAL_ASSAY:
            res = _ingest_spatial_assay(config, raw_dir)

    match res:
        case Failure(err):
            return Failure(err)
        case Success(adata):
            filtered = _apply_subset_and_subsample(adata, subset, subsample_n)
            return Success(filtered)


def is_tier1_dataset(accession: str) -> bool:
    """Check if accession is a registered Tier 1 cohort."""
    return accession in TIER1_COHORTS


# =========================================================================
# Exported Named Loaders for All 84 Tier 1 Cohorts
# =========================================================================
def load_gse274141_breast(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE274141 (Breast — Single-Cell RNA Sequencing Identifies Molecular Biomarkers Predicting Late Progression to CDK4/6 Inhibition in Patients with HR+/HER2- Metastatic Breast Cancer [FFPE])."""
    return load_tier1_cohort("GSE274141", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse222859_breast(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE222859 (Breast — Gene expression profile at single cell level of lymphocytes cells for Pan-T cell analysis)."""
    return load_tier1_cohort("GSE222859", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse274139_breast(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE274139 (Breast — Single-Cell RNA Sequencing Identifies Molecular Biomarkers Predicting Late Progression to CDK4/6 Inhibition in Patients with HR+/HER2- Metastatic Breast Cancer [tissue])."""
    return load_tier1_cohort("GSE274139", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse262288_breast(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE262288 (Breast — Single-Cell RNA Sequencing Identifies Molecular Biomarkers Predicting Response to CDK4/6 Inhibition in Metastatic HR+/HER2- Breast Cancer)."""
    return load_tier1_cohort("GSE262288", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse303346_breast(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE303346 (Breast — Single-Cell RNA-Seq Reveals Immunosuppressive Effects of Triple-Negative Breast Cancer–Derived Exosomes on NK Cell Responses)."""
    return load_tier1_cohort("GSE303346", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse254991_breast(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE254991 (Breast — Interleukin-16 Establishes a Th1-dominant Tumor Immune Microenvironment And Potentiates Immune Checkpoint Therapies [scRNA-Seq])."""
    return load_tier1_cohort("GSE254991", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse199219_breast(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE199219 (Breast — Proteogenomic integration of single-cell RNA and protein analysis identifies novel tumour-infiltrating lymphocyte phenotypes in breast cancer)."""
    return load_tier1_cohort("GSE199219", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse332708_breast(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE332708 (Breast — Hybrid In Vivo Breast Cancer Model Reveals Transcriptomic Insights into Cancer Progression with Age)."""
    return load_tier1_cohort("GSE332708", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse331487_breast(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE331487 (Breast — Immunosuppressive myeloid cells induce mesenchymal-like breast cancer stem cells by a membrane-bound TGF-β1-dependent mechanism [scRNA-seq])."""
    return load_tier1_cohort("GSE331487", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse302453_breast(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE302453 (Breast — A single-cell map of intratumoral heterogeneity during combination treatment of anti-PD1 and chemotherapy in triple-negative breast cancer)."""
    return load_tier1_cohort("GSE302453", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse299267_breast(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE299267 (Breast — Single cell gene expression profiling of triple negative breast cancer organoids)."""
    return load_tier1_cohort("GSE299267", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse205506_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE205506 (CRC — Remodeling of the Immune and Stromal Cell Compartment by PD-1 Blockade in Mismatch Repair-Deficient Colorectal Cancer)."""
    return load_tier1_cohort("GSE205506", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse146771_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE146771 (CRC — Single-Cell Analyses Inform Mechanisms of Myeloid-Targeted therapies in colon cancer)."""
    return load_tier1_cohort("GSE146771", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse164522_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE164522 (CRC — Single-cell analyses reveal phenotypic linkage between colorectal cancer and liver metastasis)."""
    return load_tier1_cohort("GSE164522", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse235917_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE235917 (CRC — First-line durvalumab and tremelimumab with chemotherapy in RAS-mutated metastatic colorectal cancer: a phase 1b/2 trial [scRNA-Seq])."""
    return load_tier1_cohort("GSE235917", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse278406_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE278406 (CRC — Phenotypic plasticity and increased tissue infiltration of TREM1+ mono-macrophages following radiotherapy in rectal cancer. [scRNA-Seq])."""
    return load_tier1_cohort("GSE278406", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse274321_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE274321 (CRC — IKKa modulates colorectal cancer metastasis by preventing tight junction stabilization and collective cell migration [scRNA-seq])."""
    return load_tier1_cohort("GSE274321", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse309346_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE309346 (CRC — PI3K and MAPK signaling nodes as divergent drivers of phenotypic plasticity in cancer-associated fibroblasts in colorectal cancer [scRNA-Seq])."""
    return load_tier1_cohort("GSE309346", raw_dir, subset=subset, subsample_n=subsample_n)

def load_cellxgene_2554a654_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort CELLxGENE_2554a654 (CRC — progressive_plasticity_during_crc_metastasis_tumor)."""
    return load_tier1_cohort("CELLxGENE_2554a654", raw_dir, subset=subset, subsample_n=subsample_n)

def load_cellxgene_5ee552f5_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort CELLxGENE_5ee552f5 (CRC — progressive_plasticity_during_crc_metastasis_non-tumor_epithelial)."""
    return load_tier1_cohort("CELLxGENE_5ee552f5", raw_dir, subset=subset, subsample_n=subsample_n)

def load_cellxgene_4b5afdf9_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort CELLxGENE_4b5afdf9 (CRC — progressive_plasticity_during_crc_metastasis_untreated_epithelial)."""
    return load_tier1_cohort("CELLxGENE_4b5afdf9", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse188711_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE188711 (CRC — Resolving the Difference Between Left-sided and Right-sided Colorectal Cancer by Single-cell Sequencing)."""
    return load_tier1_cohort("GSE188711", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse336564_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE336564 (CRC — ZFP36L2 orchestrates stress-adaptive plasticity in intestinal regeneration and colorectal cancer metastasis (PDO scRNAseq))."""
    return load_tier1_cohort("GSE336564", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse216534_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE216534 (CRC — γδ T cells are effectors of immunotherapy in cancers with HLA class I defects)."""
    return load_tier1_cohort("GSE216534", raw_dir, subset=subset, subsample_n=subsample_n)

def load_cellxgene_387acac5_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort CELLxGENE_387acac5 (CRC — progressive_plasticity_during_crc_metastasis_kg150_tumor)."""
    return load_tier1_cohort("CELLxGENE_387acac5", raw_dir, subset=subset, subsample_n=subsample_n)

def load_cellxgene_ef0d813e_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort CELLxGENE_ef0d813e (CRC — progressive_plasticity_during_crc_metastasis_kg183_tumor)."""
    return load_tier1_cohort("CELLxGENE_ef0d813e", raw_dir, subset=subset, subsample_n=subsample_n)

def load_cellxgene_2e95d453_crc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort CELLxGENE_2e95d453 (CRC — progressive_plasticity_during_crc_metastasis_kg146_tumor)."""
    return load_tier1_cohort("CELLxGENE_2e95d453", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse239676_gastric(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE239676 (Gastric — Atlas of metastatic gastric cancer links ferroptosis to disease progression and immunotherapy response)."""
    return load_tier1_cohort("GSE239676", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse319709_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE319709 (HCC — Single-cell transcriptomics reveals etiology-specific T-cell heterogeneity in hepatocellular carcinoma and implicates regulatory)."""
    return load_tier1_cohort("GSE319709", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse272347_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE272347 (HCC — Late-stage tertiary lymphoid structures in hepatocellular carcinoma treated with neoadjuvant immune checkpoint blockade [scRNA-seq])."""
    return load_tier1_cohort("GSE272347", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse318418_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE318418 (HCC — Immunosuppressive monocytes are enriched in hepatocellular carcinoma patients with liver dysfunction in a phase II trial of combination sorafenib and nivolumab [scRNA-Seq])."""
    return load_tier1_cohort("GSE318418", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse299340_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE299340 (HCC — Single-cell RNA sequencing reveals B cell-related immunosuppressive landscape and a potential suppressor in hepatocellular carcinoma)."""
    return load_tier1_cohort("GSE299340", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse282343_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE282343 (HCC — Viral-Track integrated single-cell RNA-sequencing reveals HBV lymphotropism and immunosuppressive microenvironment in HBV-associated hepatocellular carcinoma)."""
    return load_tier1_cohort("GSE282343", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse255830_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE255830 (HCC — Gene expression and T cell repertoire profile at single cell level of peripheral blood T cells after treatment with a personalized neoantigen vaccine (GNOS-PV02) and Pembrolizumab for advanced hepatocellular carcinoma.)."""
    return load_tier1_cohort("GSE255830", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse281110_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE281110 (HCC — Molecular landscape of tumor-associated tissue-resident memory T cells in tumor microenvironment of hepatocellular carcinoma [HCC_scRNA])."""
    return load_tier1_cohort("GSE281110", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse233405_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE233405 (HCC — Immunohistochemical scoring of LAG-3 in conjunction with CD8 in  the tumor microenvironment predicts response to immunotherapy in  hepatocellular carcinoma)."""
    return load_tier1_cohort("GSE233405", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse278324_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE278324 (HCC — Gene regulatory network analysis on snRNAseq revealed key regulators for hepatocellular carcinoma progression)."""
    return load_tier1_cohort("GSE278324", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse265770_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE265770 (HCC — Gene expression profile at single cell level of CD56+ natural killer cells and CD8+ T cells from blood, spleen and HCC-PDX in humanized mice)."""
    return load_tier1_cohort("GSE265770", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse318420_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE318420 (HCC — Immunosuppressive monocytes are enriched in hepatocellular carcinoma patients with liver dysfunction in a phase II trial of combination sorafenib and nivolumab)."""
    return load_tier1_cohort("GSE318420", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse272348_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE272348 (HCC — Late-stage tertiary lymphoid structures in hepatocellular carcinoma treated with neoadjuvant immune checkpoint blockade)."""
    return load_tier1_cohort("GSE272348", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse224411_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE224411 (HCC — Uncovering the spatial landscape of molecular interactions within the tumor microenvironment through latent spaces)."""
    return load_tier1_cohort("GSE224411", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse215428_hcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE215428 (HCC — Single-cell RNA sequencing of immune landscape in hepatocellular carcinoma treated with sintilimab and sorafenib)."""
    return load_tier1_cohort("GSE215428", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse296954_hnscc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE296954 (HNSCC — Differentiation of tumor-infiltrating GZMK+ effector memory T cells associates with response to neoadjuvant immunotherapy in head and neck cancer)."""
    return load_tier1_cohort("GSE296954", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse301720_hnscc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE301720 (HNSCC — Integrated single-cell and spatial analysis identifies context-dependent myeloid-T cell interactions in head and neck cancer immune checkpoint blockade response)."""
    return load_tier1_cohort("GSE301720", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse247582_hnscc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE247582 (HNSCC — Single-cell analysis of CD4+ cytotoxic T lymphocytes in human oral squamous cell carcinoma)."""
    return load_tier1_cohort("GSE247582", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse296867_hnscc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE296867 (HNSCC — Post-treatment peripheral-blood T cells from HNSCC patients undergoing neoadjuvant immunotherapy treatment.)."""
    return load_tier1_cohort("GSE296867", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse327189_hnscc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE327189 (HNSCC — Viral-based individualized neoantigen vaccine as adjuvant treatment in resected head and neck squamous cell carcinoma: immunogenicity and efficacy from a randomized Phase I trial)."""
    return load_tier1_cohort("GSE327189", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse286410_melanoma(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE286410 (Melanoma — Multi-modal Omics Analysis of a Paediatric Melanoma Highlights Mechanisms Underlying Treatment Resistance [Seq])."""
    return load_tier1_cohort("GSE286410", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse198265_melanoma(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE198265 (Melanoma — Neoantigen specific CD4+ T cells in human melanoma have diverse differentiation states and correlate with CD8+ T cell, macrophage, and B cell function)."""
    return load_tier1_cohort("GSE198265", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse242477_melanoma(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE242477 (Melanoma — Single-cell profiling of acral melanoma infiltrating lymphocytes reveals a suppressive tumor microenvironment)."""
    return load_tier1_cohort("GSE242477", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse303948_melanoma(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE303948 (Melanoma — Genome accessibility profile at single cell level for conventional dendritic cells in human metastatic melanoma samples [snATAC-Seq])."""
    return load_tier1_cohort("GSE303948", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse294273_melanoma(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE294273 (Melanoma — Tumor-resident T cells and dendritic cells form an in situ archetype for immunotherapy response in melanoma [scRNA-seq])."""
    return load_tier1_cohort("GSE294273", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse320040_melanoma(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE320040 (Melanoma — High-resolution and noninvasive profiling of the tumor microenvironment with spatial ecotypes [scRNA-seq])."""
    return load_tier1_cohort("GSE320040", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse210963_melanoma(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE210963 (Melanoma — Single-Cell RNA-Seq Analysis of Patient Myeloid-Derived Suppressor Cells and the Response to Inhibition of Bruton’s Tyrosine Kinase)."""
    return load_tier1_cohort("GSE210963", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse256291_melanoma(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE256291 (Melanoma — Comparing transcriptional profiles of CD14+  monocytes from melanoma patients and healthy donors by scRNA-seq.)."""
    return load_tier1_cohort("GSE256291", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse244983_melanoma(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE244983 (Melanoma — Molecular patterns of resistance to immune checkpoint blockade in melanoma [scRNA-Seq])."""
    return load_tier1_cohort("GSE244983", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse300446_melanoma(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE300446 (Melanoma — Spatial tumour-immune ecosystems shape the efficacy of anti-PD1 immunotherapy in primary cutaneous melanoma [snRNAseq])."""
    return load_tier1_cohort("GSE300446", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse270464_melanoma(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE270464 (Melanoma — Specific oncogene activation of the cell of origin in mucosal melanoma)."""
    return load_tier1_cohort("GSE270464", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse211068_melanoma(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE211068 (Melanoma — Integrated multiomics profiling identifies the differentiation program of regulatory T cells in human tumors [scRNA-seq])."""
    return load_tier1_cohort("GSE211068", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse241934_nsclc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE241934 (NSCLC — Neoadjuvant sintilimab plus chemotherapy in early-stage EGFR-mutant NSCLC: phase 2 trial interim results (NEOTIDE/CTONG2104))."""
    return load_tier1_cohort("GSE241934", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse270148_nsclc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE270148 (NSCLC — Myeloid progenitor dysregulation fuels immunosuppressive macrophages in tumours)."""
    return load_tier1_cohort("GSE270148", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse303762_nsclc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE303762 (NSCLC — Benchmarking long-read RNA-sequencing technologies with LongBench: a cross-platform reference dataset profiling cancer cell lines with bulk and single-cell approaches)."""
    return load_tier1_cohort("GSE303762", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse205049_nsclc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE205049 (NSCLC — Spatially Resolved Multi-Omics Single-Cell Analyses Inform Mechanisms of Immune Dysfunction in Pancreatic Cancer)."""
    return load_tier1_cohort("GSE205049", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse205354_nsclc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE205354 (NSCLC — Spatially resolved multi-omics single-cell analyses inform mechanisms of immune-dysfunction in pancreatic cancer)."""
    return load_tier1_cohort("GSE205354", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse233203_nsclc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE233203 (NSCLC — The single-cell level molecular characteristics according to combination immunotherapy response of non-small cell lung cancer patients.)."""
    return load_tier1_cohort("GSE233203", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse253718_nsclc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE253718 (NSCLC — Single-cell transcriptome reveals drug-resistance signature and immunosuppressive microenvironment in lung adenocarcinoma harboring EGFR mutation)."""
    return load_tier1_cohort("GSE253718", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse223779_nsclc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE223779 (NSCLC — Single-cell transcriptomic analysis uncovers intratumoral heterogeneity and drug-tolerant persister in ALK-rearranged lung adenocarcinoma)."""
    return load_tier1_cohort("GSE223779", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse307811_nsclc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE307811 (NSCLC — Development of antibody-drug conjugates targeting L1CAM to treat metastatic cancer)."""
    return load_tier1_cohort("GSE307811", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse285888_nsclc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE285888 (NSCLC — Single-Cell RNA Sequencing of Baseline PBMCs Predicts ICI efficacy and irAE Severity in NSCLC Patients)."""
    return load_tier1_cohort("GSE285888", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse276139_nsclc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE276139 (NSCLC — Gene expression profile at single cell level of cerebrospinal fluid (CSF) cells from lung adenocarcinoma leptomeningeal metastases patients (LUAD LM))."""
    return load_tier1_cohort("GSE276139", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse279781_pdac(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE279781 (PDAC — CD137 agonism enhances anti-PD1 induced activation of clonally expanded CD8+ T cells in a neoadjuvant pancreatic cancer clinical trial)."""
    return load_tier1_cohort("GSE279781", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse212966_pdac(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE212966 (PDAC — Single-cell RNA-seq reveals immune landscape of pancreatic cancer)."""
    return load_tier1_cohort("GSE212966", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse283206_pdac(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE283206 (PDAC — Combined Flt3L and CD40 agonism restores dendritic cell driven T cell immunity in mouse models and patients with pancreatic cancer [human])."""
    return load_tier1_cohort("GSE283206", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse311788_pdac(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE311788 (PDAC — DeCAF redefines fibroblast states uncovering multidimensional tumor-stroma relationships driving clinical tumor progression and immunotherapy response [scRNA-Seq])."""
    return load_tier1_cohort("GSE311788", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse211644_pdac(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE211644 (PDAC — Single cell transcriptomic and T cell repertoire analysis reveals trajectory of tumor - infiltrating lymphocyte states in pancreatic cancer)."""
    return load_tier1_cohort("GSE211644", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse348275_pdac(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE348275 (PDAC — Clonal lineage tracing and parallel multiomics profiling reveal transcriptional heterogeneity induced by ARID1A deficiency)."""
    return load_tier1_cohort("GSE348275", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse318413_pdac(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE318413 (PDAC — Patient-derived orthotopic xenograft models recapitulate the peritoneal dissemination of pancreatic cancer and delineate its transcriptional and regulatory programs [scRNA-seq])."""
    return load_tier1_cohort("GSE318413", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse156405_pdac(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE156405 (PDAC — Elucidation of tumor-stromal heterogeneity and the ligand-receptor interactome by single cell transcriptomics in real-world pancreatic cancer biopsies)."""
    return load_tier1_cohort("GSE156405", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse335452_pdac(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE335452 (PDAC — Single-cell RNA-seq and spatial transcriptomics characterize CD8+ exhausted T cells in pancreatic ductal adenocarcinoma)."""
    return load_tier1_cohort("GSE335452", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse285701_ccrcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE285701 (ccRCC — The paradoxical significance of CD39+CD8+ T cells in clear cell renal cell carcinoma [scRNA-seq])."""
    return load_tier1_cohort("GSE285701", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse304466_ccrcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE304466 (ccRCC — Single-cell transcriptome combined with spatial transcriptome to investigate the molecular mechanism associated with autophagy of clear renal cell carcinoma)."""
    return load_tier1_cohort("GSE304466", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse223808_ccrcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE223808 (ccRCC — Exhausted intratumoral Vδ2- γδ T cells in human kidney cancer retain effector function [scRNA-seq])."""
    return load_tier1_cohort("GSE223808", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse220313_ccrcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE220313 (ccRCC — Gene expression profile and TCR sequencing data of CD45+ cells sorted from single cell suspensions of tumors after in vitro anti-CD3 stimulation and treatments)."""
    return load_tier1_cohort("GSE220313", raw_dir, subset=subset, subsample_n=subsample_n)

def load_gse254498_ccrcc(
    raw_dir: Path,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load Tier 1 cohort GSE254498 (ccRCC — Integrating whole-exome sequencing and scRNA-seq reveal the characteristic in one clear cell renal cell carcinoma sample arising in the setting of VHL disease)."""
    return load_tier1_cohort("GSE254498", raw_dir, subset=subset, subsample_n=subsample_n)


def load_gse222315_bladder(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load bladder single-cell cohort GSE222315 (Chen et al. 2024, matched tumor vs normal adjacent)."""
    return load_tier1_cohort("GSE222315", raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)


def load_gse302781_bladder(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load metastatic urothelial carcinoma rapid autopsy snRNA-seq cohort GSE302781."""
    return load_tier1_cohort("GSE302781", raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)


def load_gse326225_bladder(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load muscle-invasive bladder cancer scRNA-seq cohort GSE326225."""
    return load_tier1_cohort("GSE326225", raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)


def load_gse301651_bladder(
    raw_dir: Path,
    auto_download: bool = True,
    subset: Mapping[str, str] | None = None,
    subsample_n: int | None = None,
) -> Result[ad.AnnData, str]:
    """Load bladder single-cell cohort GSE301651 (13 samples across primary tumor, metastasis, and PBMC)."""
    return load_tier1_cohort("GSE301651", raw_dir, auto_download=auto_download, subset=subset, subsample_n=subsample_n)


