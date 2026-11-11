"""Pan-Cancer Cancer-Associated Fibroblast (CAF) Single-Cell Atlas providers.

Integrates three peer-reviewed pan-cancer / cross-tissue human fibroblast & CAF atlases:
1. `Gao_CAF_Atlas` — Gao et al., "Cross-tissue human fibroblast atlas reveals myofibroblast
   subtypes with distinct roles in immune modulation", Cancer Cell 42(10), 1764-1783 (2024).
   DOI: 10.1016/j.ccell.2024.08.020 | PMID: 39303725 | Zenodo: 13950610 (269,899 cells)
2. `Cords_CAF_Atlas` — Cords et al., "Cancer-associated fibroblast classification in single-cell
   and spatial proteomics data", Nature Communications 14, 4294 (2023) & Cancer Cell 42(3),
   396-412 (2024). DOI: 10.1016/j.ccell.2023.12.021 | PMID: 38242124 | Zenodo: 7540604 (29,741 cells)
3. `Buechler_CAF_Atlas` — Buechler et al., "Cross-tissue organization of the fibroblast lineage",
   Nature 593, 575-579 (2021). DOI: 10.1038/s41586-021-03549-5 | PMID: 33981032 (10,355 cells)

Strictly follows functional Python standards (.agents/rules/code-style-guide.md):
- Pure Polars `.obs` harmonization preserving 100% of original author metadata columns
- Frozen Pydantic models (`ConfigDict(frozen=True)`), exhaustive `match...case` with `assert_never`
- Zero `None` / `Optional` annotations (`returns.maybe.Maybe` & `returns.result.Result`)
"""

from __future__ import annotations

from contextlib import ExitStack
from enum import StrEnum
import gc
from pathlib import Path
import shutil
from typing import Any, Final, Mapping, Sequence, assert_never
import zipfile

import anndata as ad  # type: ignore[import-untyped]
import anndata.io  # type: ignore[import-untyped]
import h5py  # type: ignore[import-untyped]
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
import scipy.sparse as sp

from ..download.fetcher import download_single_file
from ..logging import get_logger
from ..models import ChecksumSpec
from ..paths import (
    find_dataset_h5ad,
    find_repo_root,
    get_ensembl_dir,
    get_manual_download_dir,
    get_preprocessed_h5ad_path,
    get_raw_dataset_dir,
    get_scratch_dataset_dir,
    is_cscs_environment,
)
from ..preprocessing.gene_normalization import (
    _build_ensembl_projection_operator,
    _load_or_update_mapping_cache,
    _project_sparse_chunk_to_ensembl,
    ensure_ensembl_release_installed,
)
from ..storage.incremental_writer import (
    append_csr_batch_to_group,
    init_compressed_csr_group,
    write_anndata_dataframe_elem,
)
from ..storage.rds_stream import (
    SeuratStagingManifest,
    StagedAssaySpec,
    read_r_dataframe,
    stream_r_object_to_staging_h5,
)
from .atlas_single_cell import _apply_subset_and_subsample

logger = get_logger("providers.caf_atlases")


# =============================================================================
# 1. Immutable Atlas Constants, Checksums & Rich Cluster Ontologies
# =============================================================================


class CafSubstateMetadata(BaseModel):
    """Immutable biological ontology entry for a fibroblast / mural substate."""

    model_config = ConfigDict(frozen=True)

    cluster_id: str
    subtype_label: str
    broad_cell_type: str
    lineage_compartment: str
    lineage_branch: str
    functional_state: str
    canonical_markers: tuple[str, ...]


# --- 1A. Gao et al. (2024, Cancer Cell) ---
GAO_CAF_DATASET_ID: Final[str] = "Gao_CAF_Atlas"
GAO_CAF_ZENODO_RECORD: Final[str] = "13950610"
GAO_CAF_ZENODO_DOI: Final[str] = "10.5281/zenodo.13950610"
GAO_CAF_PUBLICATION_DOI: Final[str] = "10.1016/j.ccell.2024.08.020"
GAO_CAF_PUBLICATION_PMID: Final[str] = "39303725"

GAO_CAF_EXPR_FILENAME: Final[str] = "Fib_Atlas_expression.rda"
GAO_CAF_META_FILENAME: Final[str] = "Fib_Atlas_metadata.rda"
GAO_CAF_EXPR_URL: Final[str] = (
    "https://zenodo.org/records/13950610/files/Fib_Atlas_expression.rda?download=1"
)
GAO_CAF_META_URL: Final[str] = (
    "https://zenodo.org/records/13950610/files/Fib_Atlas_metadata.rda?download=1"
)
GAO_CAF_EXPR_CHECKSUM: Final[ChecksumSpec] = ChecksumSpec(
    algorithm="md5",
    expected_hash="61bf329bb40f0c883412e6b694c11d5f",
    size_bytes=Some(1_190_536_415),
)
GAO_CAF_META_CHECKSUM: Final[ChecksumSpec] = ChecksumSpec(
    algorithm="md5",
    expected_hash="1995ce2cee5a8a8494b350bda036155d",
    size_bytes=Some(7_071_242),
)

GAO_CAF_CLUSTER_SPECS: Final[Mapping[str, CafSubstateMetadata]] = {
    "c01": CafSubstateMetadata(
        cluster_id="c01",
        subtype_label="c01_MYH11_ACTA2_vSMC",
        broad_cell_type="Smooth Muscle Cell",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Mural / Pericyte / SMC",
        functional_state="Vascular smooth muscle contractile state (ACTA2+MYH11+TAGLN+)",
        canonical_markers=("ACTA2", "TAGLN", "MYH11", "MYL9", "ADIRF"),
    ),
    "c02": CafSubstateMetadata(
        cluster_id="c02",
        subtype_label="c02_ADAMDEC1_CXCL14_lpFib",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Inflammatory CAF (iCAF)",
        functional_state="Mucosal / lamina propria BMP-signaling fibroblast (ADAMDEC1+CXCL14+APOE+)",
        canonical_markers=("ADAMDEC1", "CXCL14", "APOE", "CCL11", "FABP4"),
    ),
    "c03": CafSubstateMetadata(
        cluster_id="c03",
        subtype_label="c03_COL15A1_APOD_progFib",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Universal Progenitor Fibroblast",
        functional_state="COL15A1+/APOD+ basement-membrane universal progenitor fibroblast",
        canonical_markers=("COL15A1", "APOD", "CFD", "GSN", "CXCL12"),
    ),
    "c04": CafSubstateMetadata(
        cluster_id="c04",
        subtype_label="c04_LRRC15_CTHRC1_myCAF",
        broad_cell_type="Myofibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Myofibroblastic CAF (myCAF)",
        functional_state="Immunosuppressive TGF-beta myofibroblastic CAF (LRRC15+CTHRC1+MMP11+FAP+)",
        canonical_markers=("LRRC15", "CTHRC1", "MMP11", "COL1A1", "FAP", "POSTN"),
    ),
    "c05": CafSubstateMetadata(
        cluster_id="c05",
        subtype_label="c05_PI16_MFAP5_advFib",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Universal Progenitor Fibroblast",
        functional_state="PI16+ universal adventitial stem-like progenitor fibroblast",
        canonical_markers=("PI16", "MFAP5", "IGFBP6", "CD34", "CFD", "GSN"),
    ),
    "c06": CafSubstateMetadata(
        cluster_id="c06",
        subtype_label="c06_MGP_SPARCL1_homeoFib",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Homeostatic Matrix Fibroblast",
        functional_state="Homeostatic interstitial matrix fibroblast (MGP+DCN+FBLN1+SPARCL1+)",
        canonical_markers=("MGP", "SPARCL1", "DCN", "FBLN1", "LUM"),
    ),
    "c07": CafSubstateMetadata(
        cluster_id="c07",
        subtype_label="c07_MMP3_CXCL8_iCAF",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Inflammatory CAF (iCAF)",
        functional_state="Chemokine/MMP inflammatory CAF recruiting myeloid cells (MMP1+MMP3+CXCL8+IL6+)",
        canonical_markers=("CXCL8", "CXCL1", "CXCL2", "MMP3", "MMP1", "IL6", "TNFAIP6"),
    ),
    "c08": CafSubstateMetadata(
        cluster_id="c08",
        subtype_label="c08_RGS5_NOTCH3_Pericyte",
        broad_cell_type="Pericyte",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Mural / Pericyte / SMC",
        functional_state="Angiogenic microvascular pericyte (RGS5+NOTCH3+PDGFRB+COL4A1+)",
        canonical_markers=("RGS5", "NOTCH3", "NDUFA4L2", "PDGFRB", "COL4A1", "MCAM"),
    ),
    "c09": CafSubstateMetadata(
        cluster_id="c09",
        subtype_label="c09_COL1A1_MDK_fetalFib",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Developmental Stromal Fibroblast",
        functional_state="Fetal developmental stromal fibroblast (COL1A1+MDK+HBG2+)",
        canonical_markers=("COL1A1", "COL3A1", "MDK", "HBG2", "PTN"),
    ),
    "c10": CafSubstateMetadata(
        cluster_id="c10",
        subtype_label="c10_PRG4_CRTAC1_liningFib",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Tissue-Specialized Fibroblast",
        functional_state="Synovial lining lubricin fibroblast (PRG4+CRTAC1+CLU+)",
        canonical_markers=("PRG4", "CRTAC1", "PLA2G2A", "CLU", "HTRA1"),
    ),
    "c11": CafSubstateMetadata(
        cluster_id="c11",
        subtype_label="c11_TAGLN_HOPX_myoFib",
        broad_cell_type="Myofibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Myofibroblastic CAF (myCAF)",
        functional_state="Contractile tissue-remodeling myofibroblast (TAGLN+HOPX+BGN+ACTA2+)",
        canonical_markers=("TAGLN", "HOPX", "BGN", "ACTA2", "TPM2"),
    ),
    "c12": CafSubstateMetadata(
        cluster_id="c12",
        subtype_label="c12_C3_KRT18_mesoFib",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Antigen-Presenting / Mesothelial CAF",
        functional_state="Mesothelial-like complement/antigen-presenting fibroblast (C3+KRT18+SLPI+CD74+)",
        canonical_markers=("C3", "RARRES1", "KRT18", "SLPI", "TM4SF1", "CD74"),
    ),
    "c13": CafSubstateMetadata(
        cluster_id="c13",
        subtype_label="c13_C7_SFRP1_compFib",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Inflammatory CAF (iCAF)",
        functional_state="Complement-secreting immunomodulatory fibroblast (C7+SELENOP+A2M+SFRP1+)",
        canonical_markers=("C7", "SELENOP", "A2M", "SFRP1", "C3"),
    ),
    "c14": CafSubstateMetadata(
        cluster_id="c14",
        subtype_label="c14_HSPA1A_DNAJB1_hspCAF",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Stress-Response CAF (hspCAF)",
        functional_state="Heat-shock stress-response fibroblast (HSPA1A+HSPA1B+DNAJB1+HSP90AA1+)",
        canonical_markers=("HSPA1A", "HSPA1B", "DNAJB1", "HSP90AA1", "FOS"),
    ),
    "c15": CafSubstateMetadata(
        cluster_id="c15",
        subtype_label="c15_POSTN_BMP4_teloFib",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Tissue-Specialized Fibroblast",
        functional_state="Subepithelial BMP4+/FOXL1+ telocyte niche fibroblast (POSTN+CXCL14+F3+BMP4+)",
        canonical_markers=("POSTN", "BMP4", "F3", "CXCL14", "PLAT", "SOX6"),
    ),
    "c16": CafSubstateMetadata(
        cluster_id="c16",
        subtype_label="c16_SFRP2_SFRP4_matFib",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Myofibroblastic CAF (myCAF)",
        functional_state="SFRP2+/SFRP4+ Wnt-modulating matrix fibroblast (COL1A1+SFRP2+SFRP4+CCDC80+)",
        canonical_markers=("SFRP2", "SFRP4", "CCDC80", "COL1A1", "COMP"),
    ),
    "c17": CafSubstateMetadata(
        cluster_id="c17",
        subtype_label="c17_STMN1_MDK_prolifCAF",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Proliferating CAF (pCAF)",
        functional_state="Cell-cycle proliferating / developmental fibroblast (STMN1+MDK+MKI67+TOP2A+)",
        canonical_markers=("STMN1", "MDK", "MARCKSL1", "RBP1", "MKI67", "TOP2A"),
    ),
    "c18": CafSubstateMetadata(
        cluster_id="c18",
        subtype_label="c18_ACTG2_MYLK_viscSMC",
        broad_cell_type="Smooth Muscle Cell",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Mural / Pericyte / SMC",
        functional_state="Visceral smooth muscle contractile cell (ACTG2+MYLK+DES+ACTA2+)",
        canonical_markers=("ACTG2", "MYLK", "ACTA2", "TAGLN", "DES"),
    ),
    "c19": CafSubstateMetadata(
        cluster_id="c19",
        subtype_label="c19_CTSK_CTHRC1_invCAF",
        broad_cell_type="Myofibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Myofibroblastic CAF (myCAF)",
        functional_state="Cathepsin-K invasive ECM-degrading myCAF (CTSK+CTHRC1+FAP+COL6A3+)",
        canonical_markers=("CTSK", "CTHRC1", "COL1A1", "COL6A3", "FAP", "TGFB1"),
    ),
    "c20": CafSubstateMetadata(
        cluster_id="c20",
        subtype_label="c20_CXCL12_CCL19_immunoCAF",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Lymphoid-Organizing / Immunomodulatory CAF",
        functional_state="CXCL12+/CCL19+ lymphoid-recruiting immunomodulatory fibroblast (CXCL12+PTGDS+CCL19+)",
        canonical_markers=("CXCL12", "PTGDS", "PLA2G2A", "CCL19", "APOE"),
    ),
}


# --- 1B. Cords et al. (2024, Cancer Cell & 2023, Nat. Commun.) ---
CORDS_CAF_DATASET_ID: Final[str] = "Cords_CAF_Atlas"
CORDS_CAF_ZENODO_RECORD: Final[str] = "7540604"
CORDS_CAF_ZENODO_DOI: Final[str] = "10.5281/zenodo.7540604"
CORDS_CAF_PUBLICATION_DOI: Final[str] = "10.1016/j.ccell.2023.12.021"
CORDS_CAF_PUBLICATION_PMID: Final[str] = "38242124"
CORDS_CAF_ZIP_FILENAME: Final[str] = "scRNA-seq_dataobjects.zip"
CORDS_CAF_ZIP_URL: Final[str] = (
    "https://zenodo.org/records/7540604/files/scRNA-seq_dataobjects.zip?download=1"
)
CORDS_CAF_ZIP_CHECKSUM: Final[ChecksumSpec] = ChecksumSpec(
    algorithm="md5",
    expected_hash="af5c346c93944c643d2bad97950e0f52",
    size_bytes=Some(2_172_690_224),
)

CORDS_COHORT_RDS_MEMBERS: Final[tuple[tuple[str, str, str, str], ...]] = (
    ("scRNA-seq/BREAST_fibro_tumour.rds", "Cords_Breast", "BC", "Breast"),
    ("scRNA-seq/PDAC_gse21_fibro_tumour.rds", "Cords_PDAC_GSE212966", "PDAC", "Pancreas"),
    ("scRNA-seq/PDAC_gse15_fibro_tumour.rds", "Cords_PDAC_GSE155698", "PDAC", "Pancreas"),
    ("scRNA-seq/COLON_fibro_tumour.rds", "Cords_Colon", "CRC", "Colon / Rectum"),
    ("scRNA-seq/LUNG_fibro_tumour.rds", "Cords_Lung", "NSCLC", "Lung"),
    ("scRNA-seq/HNSCC_fibro_tumour.rds", "Cords_HNSCC", "HNSCC", "Head and Neck"),
)

CORDS_CAF_CLUSTER_SPECS: Final[Mapping[str, CafSubstateMetadata]] = {
    "mCAF": CafSubstateMetadata(
        cluster_id="mCAF",
        subtype_label="mCAF_MMP11_POSTN_matrix",
        broad_cell_type="Myofibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Myofibroblastic CAF (myCAF)",
        functional_state="Matrix-remodeling desmoplastic myofibroblastic CAF (MMP11+COL11A1+POSTN+FAP+LRRC15+)",
        canonical_markers=("MMP11", "COL11A1", "POSTN", "FAP", "LRRC15", "COL1A1"),
    ),
    "iCAF": CafSubstateMetadata(
        cluster_id="iCAF",
        subtype_label="iCAF_CFD_CXCL12_inflammatory",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Inflammatory CAF (iCAF)",
        functional_state="Inflammatory cytokine/complement CAF (CFD+PLA2G2A+IL6+CXCL12+CD34+)",
        canonical_markers=("CFD", "PLA2G2A", "CXCL12", "IL6", "CD34", "DPT"),
    ),
    "vCAF": CafSubstateMetadata(
        cluster_id="vCAF",
        subtype_label="vCAF_NOTCH3_MCAM_vascular",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Vascular CAF / Mural (vCAF)",
        functional_state="Vascular / perivascular angiogenic CAF (NOTCH3+COL18A1+MCAM+ACTA2+)",
        canonical_markers=("NOTCH3", "MCAM", "COL18A1", "ACTA2", "STEAP4"),
    ),
    "tpCAF": CafSubstateMetadata(
        cluster_id="tpCAF",
        subtype_label="tpCAF_MME_NDRG1_tumorlike",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Glycolytic / Proliferating CAF (tCAF/pCAF)",
        functional_state="Tumor-like hypoxic/glycolytic & proliferating CAF (MME+TMEM158+NDRG1+SLC2A1+MKI67+)",
        canonical_markers=("MME", "TMEM158", "NDRG1", "ENO1", "PGK1", "MKI67"),
    ),
    "hsp_tpCAF": CafSubstateMetadata(
        cluster_id="hsp_tpCAF",
        subtype_label="hsp_tpCAF_HSPA1A_stress",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Stress-Response CAF (hspCAF)",
        functional_state="Heat-shock stress-activated tumor-like CAF (HSPA1A+HSPA1B+HSP90AA1+DNAJB1+)",
        canonical_markers=("HSPA1A", "HSPA1B", "DNAJB1", "HSP90AA1", "MME"),
    ),
    "IDO_CAF": CafSubstateMetadata(
        cluster_id="IDO_CAF",
        subtype_label="ifnCAF_IDO1_CXCL9_interferon",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Interferon-Response CAF (ifnCAF)",
        functional_state="Interferon-stimulated IDO1+CXCL9+CXCL10+ immunomodulatory CAF (ifnCAF)",
        canonical_markers=("IDO1", "CXCL9", "CXCL10", "CXCL11", "STAT1", "ISG15"),
    ),
    "apCAF": CafSubstateMetadata(
        cluster_id="apCAF",
        subtype_label="apCAF_CD74_HLADRA_antigen",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Antigen-Presenting / Mesothelial CAF",
        functional_state="MHC-II antigen-presenting CAF (CD74+HLA-DRA+HLA-DRB1+SLPI+)",
        canonical_markers=("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "SLPI"),
    ),
    "rCAF": CafSubstateMetadata(
        cluster_id="rCAF",
        subtype_label="rCAF_CCL21_CCL19_reticular",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Lymphoid-Organizing / Immunomodulatory CAF",
        functional_state="TLS-organizing reticular fibroblast (CCL21+CCL19+FDCSP+CR2+)",
        canonical_markers=("CCL21", "CCL19", "CLU", "CXCL13"),
    ),
    "dCAF": CafSubstateMetadata(
        cluster_id="dCAF",
        subtype_label="dCAF_MKI67_TOP2A_dividing",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Proliferating CAF (pCAF)",
        functional_state="Dividing / mitotic CAF (MKI67+TOP2A+BIRC5+STMN1+)",
        canonical_markers=("MKI67", "TOP2A", "BIRC5", "UBE2C", "STMN1"),
    ),
    "Pericyte": CafSubstateMetadata(
        cluster_id="Pericyte",
        subtype_label="Pericyte_RGS5_PDGFRB_mural",
        broad_cell_type="Pericyte",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Mural / Pericyte / SMC",
        functional_state="Microvascular mural pericyte (RGS5+CSP4+PDGFRB+NDUFA4L2+)",
        canonical_markers=("RGS5", "PDGFRB", "NDUFA4L2", "HIGD1B", "COX4I2"),
    ),
    "other": CafSubstateMetadata(
        cluster_id="other",
        subtype_label="other_CAF_unassigned",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="CAF / Fibroblast",
        functional_state="Minor unassigned stromal fibroblast cluster",
        canonical_markers=("DCN", "LUM", "VIM", "COL1A2"),
    ),
}


# --- 1C. Buechler et al. (2021, Nature) ---
BUECHLER_CAF_DATASET_ID: Final[str] = "Buechler_CAF_Atlas"
BUECHLER_CAF_PUBLICATION_DOI: Final[str] = "10.1038/s41586-021-03549-5"
BUECHLER_CAF_PUBLICATION_PMID: Final[str] = "33981032"
BUECHLER_CAF_RDS_FILENAME: Final[str] = "Human_PS_Fibro.RDS"
BUECHLER_CAF_GDRIVE_URL: Final[str] = (
    "https://drive.usercontent.google.com/download?id=1n-XLSYacOIXSlFQbe0XZP3BhhW0CBl7Z"
    "&export=download&confirm=t"
)
BUECHLER_CAF_CHECKSUM: Final[ChecksumSpec] = ChecksumSpec(
    algorithm="md5",
    expected_hash="f357cf64d2239b7ec7d2988ecb1334dd",
    size_bytes=Some(174_346_125),
)

BUECHLER_CAF_CLUSTER_SPECS: Final[Mapping[str, CafSubstateMetadata]] = {
    "COL3A1": CafSubstateMetadata(
        cluster_id="COL3A1",
        subtype_label="0_COL3A1_panCAF",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Homeostatic Matrix Fibroblast",
        functional_state="COL3A1+/COL15A1+ universal parenchymal matrix fibroblast",
        canonical_markers=("COL3A1", "COL15A1", "COL1A2", "SPARC", "DCN"),
    ),
    "LRRC15": CafSubstateMetadata(
        cluster_id="LRRC15",
        subtype_label="1_LRRC15_myCAF",
        broad_cell_type="Myofibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Myofibroblastic CAF (myCAF)",
        functional_state="TGF-beta-driven LRRC15+CTHRC1+FAP+ myofibroblastic CAF",
        canonical_markers=("LRRC15", "CTHRC1", "COL11A1", "MMP11", "FAP", "POSTN"),
    ),
    "PI16": CafSubstateMetadata(
        cluster_id="PI16",
        subtype_label="2_PI16_advProgenitor",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Universal Progenitor Fibroblast",
        functional_state="PI16+CD34+DPT+ universal adventitial stem-like fibroblast progenitor",
        canonical_markers=("PI16", "CD34", "MFAP5", "PCOLCE2", "IGFBP6", "CFD"),
    ),
    "ADAMDEC1": CafSubstateMetadata(
        cluster_id="ADAMDEC1",
        subtype_label="3_ADAMDEC1_laminaPropria",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Inflammatory CAF (iCAF)",
        functional_state="ADAMDEC1+CXCL14+ mucosal lamina propria inflammatory fibroblast",
        canonical_markers=("ADAMDEC1", "CXCL14", "CCL11", "BMP4", "APOE"),
    ),
    "NPNT": CafSubstateMetadata(
        cluster_id="NPNT",
        subtype_label="4_NPNT_alveolar",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Tissue-Specialized Fibroblast",
        functional_state="NPNT+CES1+INMT+ specialized lung alveolar fibroblast",
        canonical_markers=("NPNT", "CES1", "INMT", "MACF1", "TCF21"),
    ),
    "CCL19": CafSubstateMetadata(
        cluster_id="CCL19",
        subtype_label="5_CCL19_FRC_TLS",
        broad_cell_type="Fibroblast",
        lineage_compartment="STROMAL_CAF",
        lineage_branch="Lymphoid-Organizing / Immunomodulatory CAF",
        functional_state="CCL19+CCL21+CXCL13+ lymphoid-organizing fibroblastic reticular cell (TLS niche)",
        canonical_markers=("CCL19", "CCL21", "CXCL12", "CXCL13", "CLU"),
    ),
}

_BUECHLER_NUM_TO_NAME: Final[Mapping[str, str]] = {
    "0": "COL3A1",
    "1": "LRRC15",
    "2": "PI16",
    "3": "ADAMDEC1",
    "4": "NPNT",
    "5": "CCL19",
}


class CafAtlasChoice(StrEnum):
    """Supported pan-cancer CAF reference atlases."""

    GAO = "Gao_CAF_Atlas"
    CORDS = "Cords_CAF_Atlas"
    BUECHLER = "Buechler_CAF_Atlas"


# =============================================================================
# 2. Pure Functional Core: Polars `.obs` Harmonizers (100% Column Preservation)
# =============================================================================


def _pandas_to_polars_preserving_all(
    raw_obs: pd.DataFrame,
    cell_names_fallback: Sequence[str] = (),
) -> pl.DataFrame:
    """Convert any Pandas DataFrame into a typed Polars DataFrame without dropping any column."""
    n_rows = len(raw_obs) if len(raw_obs) > 0 else len(cell_names_fallback)
    if n_rows == 0:
        return pl.DataFrame({"cell_id": pl.Series("cell_id", [], dtype=pl.String)})

    barcodes = (
        [str(idx) for idx in raw_obs.index]
        if len(raw_obs) == n_rows and len(raw_obs.index) == n_rows
        else [str(c) for c in cell_names_fallback[:n_rows]]
    )
    records: dict[str, Any] = {"cell_id": barcodes}
    for col in raw_obs.columns:
        col_str = str(col).strip()
        if not col_str or col_str == "cell_id":
            continue
        series = raw_obs[col]
        if pd.api.types.is_bool_dtype(series):
            records[col_str] = series.fillna(False).to_numpy(dtype=bool)
        elif pd.api.types.is_numeric_dtype(series):
            records[col_str] = pd.to_numeric(series, errors="coerce").to_numpy(dtype=np.float64)
        else:
            records[col_str] = [
                "NA" if pd.isna(v) or str(v).strip() == "" else str(v).strip()
                for v in series.tolist()
            ]
    return pl.DataFrame(records)


def _finalize_polars_to_obs_pandas(df: pl.DataFrame) -> pd.DataFrame:
    """Convert enriched Polars `.obs` table back to an AnnData-safe Pandas DataFrame."""
    if df.is_empty():
        return pd.DataFrame(index=pd.Index([], name="cell_id_index"))
    out_pdf = df.to_pandas()
    out_pdf = out_pdf.set_index("cell_id", drop=False)
    out_pdf.index.name = "cell_id_index"
    return out_pdf


def _map_gao_tissue_and_phenotype_to_cancer(tissue: str, phenotype: str) -> tuple[str, str, str]:
    """Map `(Tissue, Phenotype)` in Gao et al. (2024) to `(cancer_code, cancer_type, organ)`."""
    organ_map: Mapping[str, str] = {
        "Breast": "Breast",
        "Colon": "Colon / Rectum",
        "Kidney": "Kidney",
        "Synovium": "Synovium",
        "Lung": "Lung",
        "Pancreas": "Pancreas",
        "Ovary": "Ovary",
        "Skin": "Skin",
        "Stomach": "Stomach",
        "Liver": "Liver",
        "HNSC": "Head and Neck",
    }
    tumor_code_map: Mapping[str, tuple[str, str]] = {
        "Breast": ("BC", "Breast Cancer"),
        "Colon": ("CRC", "Colorectal Cancer"),
        "Kidney": ("RCC", "Renal Cell Carcinoma"),
        "Lung": ("NSCLC", "Non-Small Cell Lung"),
        "Pancreas": ("PDAC", "Pancreatic Ductal Adenocarcinoma"),
        "Ovary": ("OV", "Ovarian Cancer"),
        "Skin": ("MEL", "Skin Cancer / Melanoma"),
        "Stomach": ("STAD", "Gastric Adenocarcinoma"),
        "Liver": ("LIHC", "Liver Cancer"),
        "HNSC": ("HNSC", "Head and Neck Squamous Cell Carcinoma"),
        "Synovium": ("NON_MALIGNANT", "Synovial Inflammation / Non-Malignant"),
    }
    organ = organ_map.get(tissue, tissue if tissue else "Multi-organ")
    pheno_low = phenotype.strip().lower()
    if pheno_low in {"tumor", "metastasis"}:
        c_code, c_name = tumor_code_map.get(tissue, (tissue.upper() or "PAN_CANCER", f"{tissue} Cancer"))
        return c_code, c_name, organ
    if pheno_low in {"inflamed", "non-inflamed"}:
        return "INFLAMMATORY_DISEASE", f"{organ} ({phenotype})", organ
    return "NON_TUMOR_CONTROL", f"{organ} ({phenotype or 'Healthy/Normal'})", organ


def _normalize_gao_sample_type(phenotype: str) -> str:
    """Normalize Gao `Phenotype` (`Tumor`, `Metastasis`, `Normal`, `Healthy`, `Inflamed`, etc.)."""
    low = phenotype.strip().lower()
    match low:
        case "tumor":
            return "tumor"
        case "metastasis":
            return "metastasis"
        case "normal":
            return "adjacent_normal"
        case "healthy":
            return "healthy"
        case "inflamed":
            return "inflamed"
        case "non-inflamed":
            return "non_inflamed"
        case "fetal":
            return "fetal"
        case _:
            return "unknown"


def harmonize_gao_caf_obs_polars(
    raw_meta_df: pd.DataFrame,
    expr_cell_names: Sequence[str] = (),
) -> tuple[pd.DataFrame, Maybe[np.ndarray[Any, Any]]]:
    """Harmonize Gao et al. (2024) fibroblast metadata while preserving 100% of author columns.

    Also extracts `(UMAP_1, UMAP_2)` into a `(n_cells, 2)` `float32` array wrapped in `Maybe`
    so `.obsm["X_umap"]` is populated alongside `.obs["UMAP_1"]` and `.obs["UMAP_2"]`.
    """
    pl_meta = _pandas_to_polars_preserving_all(raw_meta_df)
    if pl_meta.is_empty():
        return _finalize_polars_to_obs_pandas(pl_meta), Nothing

    if "CellName" in pl_meta.columns:
        pl_meta = pl_meta.with_columns(pl.col("CellName").cast(pl.String).alias("cell_id"))

    if len(expr_cell_names) > 0:
        order_df = pl.DataFrame({"cell_id": [str(c) for c in expr_cell_names]})
        pl_meta = order_df.join(pl_meta, on="cell_id", how="left")

    cluster_col = (
        pl.col("Cluster").cast(pl.String).fill_null("NA")
        if "Cluster" in pl_meta.columns
        else pl.lit("NA")
    )
    tissue_col = (
        pl.col("Tissue").cast(pl.String).fill_null("Multi-organ")
        if "Tissue" in pl_meta.columns
        else pl.lit("Multi-organ")
    )
    pheno_col = (
        pl.col("Phenotype").cast(pl.String).fill_null("Tumor")
        if "Phenotype" in pl_meta.columns
        else pl.lit("Tumor")
    )
    study_col = (
        pl.col("Dataset").cast(pl.String).fill_null("Gao2024")
        if "Dataset" in pl_meta.columns
        else pl.lit("Gao2024")
    )
    patient_col = (
        pl.col("Patient").cast(pl.String).fill_null("NA")
        if "Patient" in pl_meta.columns
        else pl.lit("NA")
    )
    sample_col = (
        pl.col("Sample").cast(pl.String).fill_null("NA")
        if "Sample" in pl_meta.columns
        else pl.lit("NA")
    )

    subtype_map = {k: v.subtype_label for k, v in GAO_CAF_CLUSTER_SPECS.items()}
    broad_map = {k: v.broad_cell_type for k, v in GAO_CAF_CLUSTER_SPECS.items()}
    comp_map = {k: v.lineage_compartment for k, v in GAO_CAF_CLUSTER_SPECS.items()}
    branch_map = {k: v.lineage_branch for k, v in GAO_CAF_CLUSTER_SPECS.items()}
    func_map = {k: v.functional_state for k, v in GAO_CAF_CLUSTER_SPECS.items()}

    base_df = pl_meta.lazy().with_columns(
        cluster_col.alias("author_cluster_id"),
        tissue_col.alias("_tissue_clean"),
        pheno_col.alias("_pheno_clean"),
        study_col.alias("source_study"),
        pl.when(patient_col.is_in(["", "NA", "nan"]))
        .then(sample_col)
        .otherwise(patient_col)
        .alias("_patient_clean"),
        sample_col.alias("_sample_clean"),
    ).collect()

    tissues_list = base_df["_tissue_clean"].to_list()
    phenos_list = base_df["_pheno_clean"].to_list()
    triples = [
        _map_gao_tissue_and_phenotype_to_cancer(str(t), str(p))
        for t, p in zip(tissues_list, phenos_list, strict=False)
    ]
    c_codes = [tr[0] for tr in triples]
    c_types = [tr[1] for tr in triples]
    organs = [tr[2] for tr in triples]
    sample_types = [_normalize_gao_sample_type(str(p)) for p in phenos_list]

    enriched = (
        base_df.lazy()
        .with_columns(
            pl.Series("cancer_code", c_codes, dtype=pl.String),
            pl.Series("cancer_type", c_types, dtype=pl.String),
            pl.Series("organ", organs, dtype=pl.String),
            pl.Series("sample_type", sample_types, dtype=pl.String),
            pl.col("_pheno_clean").alias("tissue"),
            pl.col("_patient_clean").alias("patient"),
            pl.col("_patient_clean").alias("patient_id"),
            pl.when(pl.col("_sample_clean").is_in(["", "NA", "nan"]))
            .then(pl.col("_patient_clean"))
            .otherwise(pl.col("_sample_clean"))
            .alias("sample"),
            pl.col("author_cluster_id")
            .replace_strict(broad_map, default=pl.lit("Fibroblast"))
            .alias("cell_type"),
            pl.col("author_cluster_id")
            .replace_strict(subtype_map, default=pl.col("author_cluster_id"))
            .alias("cell_subtype"),
            pl.col("author_cluster_id")
            .replace_strict(comp_map, default=pl.lit("STROMAL_CAF"))
            .alias("lineage_compartment"),
            pl.col("author_cluster_id")
            .replace_strict(branch_map, default=pl.lit("CAF / Fibroblast"))
            .alias("lineage_branch"),
            pl.col("author_cluster_id")
            .replace_strict(func_map, default=pl.lit("Cross-tissue fibroblast state"))
            .alias("functional_state"),
        )
        .with_columns(
            pl.col("sample").alias("sample_id"),
            pl.col("cell_subtype").alias("cell_state"),
            pl.col("author_cluster_id").alias("cell_type_author"),
            pl.lit(GAO_CAF_DATASET_ID).alias("dataset"),
            pl.lit(GAO_CAF_DATASET_ID).alias("dataset_id"),
        )
        .drop(["_tissue_clean", "_pheno_clean", "_patient_clean", "_sample_clean"])
        .collect()
    )

    umap_maybe: Maybe[np.ndarray[Any, Any]] = Nothing
    if "UMAP_1" in enriched.columns and "UMAP_2" in enriched.columns:
        u1 = enriched["UMAP_1"].cast(pl.Float32, strict=False).fill_null(float("nan")).to_numpy()
        u2 = enriched["UMAP_2"].cast(pl.Float32, strict=False).fill_null(float("nan")).to_numpy()
        xy = np.column_stack([u1, u2]).astype(np.float32, copy=False)
        if bool(np.isfinite(xy).any()):
            umap_maybe = Some(xy)

    return _finalize_polars_to_obs_pandas(enriched), umap_maybe


def _canonicalize_cords_caftype(raw_label: str) -> str:
    """Normalize Cords et al. cohort-level `CAFtype` / `cluster_ft` tokens onto the 11 canonical CAF phenotypes."""
    cleaned = raw_label.strip()
    low = cleaned.lower()
    match low:
        case "mcaf":
            return "mCAF"
        case "icaf":
            return "iCAF"
        case "vcaf" | "vcaf+rcaf":
            return "vCAF"
        case "tpcaf" | "tcaf" | "pcaf":
            return "tpCAF"
        case "hsp_tpcaf" | "hspcaf":
            return "hsp_tpCAF"
        case "ido_caf" | "ifncaf":
            return "IDO_CAF"
        case "apcaf":
            return "apCAF"
        case "rcaf" | "rcaf+apcaf":
            return "rCAF"
        case "dcaf":
            return "dCAF"
        case "pericyte" | "pericytes" | "pericyte type 1" | "pericyte type 2":
            return "Pericyte"
        case _:
            return "other"


def harmonize_cords_caf_obs_polars(
    raw_obs: pd.DataFrame,
    cohort_id: str,
    cancer_code: str,
    organ: str,
    cell_names_fallback: Sequence[str] = (),
) -> pd.DataFrame:
    """Harmonize a Cords et al. (2024) tumor CAF cohort `.obs` table while retaining curated columns."""
    cancer_full_map: Mapping[str, str] = {
        "BC": "Breast Cancer",
        "PDAC": "Pancreatic Ductal Adenocarcinoma",
        "CRC": "Colorectal Cancer",
        "NSCLC": "Non-Small Cell Lung",
        "HNSCC": "Head and Neck Squamous Cell Carcinoma",
    }

    # Drop verbose multi-resolution SNN sweep & GEO boilerplate columns only when >30 cols exist?
    # Wait: user requested preserving rich metadata like cell types! Let's keep ALL author phenotype,
    # cell-type, QC, patient, sample, activation, and cluster columns, and only drop multi-kilobyte
    # GEO contact/protocol prose strings (`treatment_protocol_ch1`, `extract_protocol_ch1*`, `data_processing*`, `contact_*`)
    # that repeat identical 1,500-char text across rows, OR keep all informative columns cleanly!
    # Even better: keep every non-prose column AND keep `geo_accession`, `platform_id`, `instrument_model`, `tumor.site.ch1`!
    prose_boilerplate = {
        "treatment_protocol_ch1",
        "extract_protocol_ch1",
        "extract_protocol_ch1.1",
        "data_processing",
        "data_processing.1",
        "data_processing.2",
        "data_processing.3",
        "contact_name",
        "contact_email",
        "contact_institute",
        "contact_address",
        "contact_city",
        "contact_state",
        "contact_zip.postal_code",
        "contact_country",
        "supplementary_file_1",
        "relation",
    }
    retained_cols = [
        c
        for c in raw_obs.columns
        if str(c) not in prose_boilerplate
        and not (
            any(
                str(c).startswith(pref)
                for pref in (
                    "SCT_snn_res.",
                    "all.colon_res.",
                    "all.subs_res.",
                    "all.cells_res.",
                    "all.HNSCC_res.",
                    "all.fibros_res.",
                    "fibro.pdac_res.",
                )
            )
        )
    ]
    trimmed_obs = raw_obs[retained_cols].copy()
    pl_df = _pandas_to_polars_preserving_all(trimmed_obs, cell_names_fallback=cell_names_fallback)
    if pl_df.is_empty():
        return _finalize_polars_to_obs_pandas(pl_df)

    # Prefix cell_id with cohort_id to guarantee global uniqueness across the 6 tumor cohorts
    pl_df = pl_df.with_columns(
        pl.col("cell_id").alias("original_barcode"),
        (pl.lit(f"{cohort_id}__") + pl.col("cell_id")).alias("cell_id"),
    )

    caf_candidates = ("CAFtype", "cluster_ft", "phenotype_final", "Cell_subtype", "seurat_clusters")
    avail_caf = [c for c in caf_candidates if c in pl_df.columns]
    raw_caf_expr = (
        pl.coalesce(
            [
                pl.when(
                    pl.col(c).cast(pl.String).is_not_null()
                    & ~pl.col(c).cast(pl.String).is_in(["", "NA", "nan"])
                )
                .then(pl.col(c).cast(pl.String))
                .otherwise(pl.lit(""))
                for c in avail_caf
            ]
        )
        if avail_caf
        else pl.lit("other")
    )
    # Note: pl.coalesce skips nulls, not empty strings; build nested pl.when cleanly:
    nested_caf_expr: pl.Expr = pl.lit("other")
    for c_name in reversed(avail_caf):
        c_str = pl.col(c_name).cast(pl.String).fill_null("NA")
        nested_caf_expr = (
            pl.when(~c_str.is_in(["", "NA", "nan"]))
            .then(c_str)
            .otherwise(nested_caf_expr)
        )
    _ = raw_caf_expr

    sample_candidates = ("Sample", "sample", "orig.ident")
    avail_samp = [c for c in sample_candidates if c in pl_df.columns]
    nested_samp_expr: pl.Expr = pl.lit(cohort_id)
    for s_name in reversed(avail_samp):
        s_str = pl.col(s_name).cast(pl.String).fill_null("NA")
        nested_samp_expr = (
            pl.when(~s_str.is_in(["", "NA", "nan", "SeuratProject"]))
            .then(s_str)
            .otherwise(nested_samp_expr)
        )

    patient_candidates = ("Patient", "Sample", "sample", "orig.ident")
    avail_pat = [c for c in patient_candidates if c in pl_df.columns]
    nested_pat_expr: pl.Expr = pl.lit(cohort_id)
    for p_name in reversed(avail_pat):
        p_str = pl.col(p_name).cast(pl.String).fill_null("NA")
        nested_pat_expr = (
            pl.when(~p_str.is_in(["", "NA", "nan", "SeuratProject"]))
            .then(p_str)
            .otherwise(nested_pat_expr)
        )

    act_expr = (
        pl.col("activation_state").cast(pl.String).fill_null("NA")
        if "activation_state" in pl_df.columns
        else pl.lit("NA")
    )

    subtype_map = {k: v.subtype_label for k, v in CORDS_CAF_CLUSTER_SPECS.items()}
    broad_map = {k: v.broad_cell_type for k, v in CORDS_CAF_CLUSTER_SPECS.items()}
    comp_map = {k: v.lineage_compartment for k, v in CORDS_CAF_CLUSTER_SPECS.items()}
    branch_map = {k: v.lineage_branch for k, v in CORDS_CAF_CLUSTER_SPECS.items()}
    func_map = {k: v.functional_state for k, v in CORDS_CAF_CLUSTER_SPECS.items()}

    enriched = (
        pl_df.lazy()
        .with_columns(
            nested_caf_expr.alias("cell_type_author"),
            nested_caf_expr.map_elements(_canonicalize_cords_caftype, return_dtype=pl.String).alias(
                "author_cluster_id"
            ),
            act_expr.alias("caf_activation_state"),
            pl.lit(cohort_id).alias("source_study"),
            nested_pat_expr.alias("patient"),
            nested_pat_expr.alias("patient_id"),
            nested_samp_expr.alias("sample"),
            nested_samp_expr.alias("sample_id"),
            pl.lit(organ).alias("organ"),
            pl.lit(cancer_code).alias("cancer_code"),
            pl.lit(cancer_full_map.get(cancer_code, cancer_code)).alias("cancer_type"),
            pl.lit("tumor").alias("sample_type"),
            pl.lit("Tumor").alias("tissue"),
        )
        .with_columns(
            pl.col("author_cluster_id")
            .replace_strict(broad_map, default=pl.lit("Fibroblast"))
            .alias("cell_type"),
            pl.col("author_cluster_id")
            .replace_strict(subtype_map, default=pl.col("author_cluster_id"))
            .alias("cell_subtype"),
            pl.col("author_cluster_id")
            .replace_strict(comp_map, default=pl.lit("STROMAL_CAF"))
            .alias("lineage_compartment"),
            pl.col("author_cluster_id")
            .replace_strict(branch_map, default=pl.lit("CAF / Fibroblast"))
            .alias("lineage_branch"),
            pl.col("author_cluster_id")
            .replace_strict(func_map, default=pl.lit("Human tumor CAF state"))
            .alias("functional_state"),
        )
        .with_columns(
            pl.col("cell_subtype").alias("cell_state"),
            pl.lit(CORDS_CAF_DATASET_ID).alias("dataset"),
            pl.lit(CORDS_CAF_DATASET_ID).alias("dataset_id"),
        )
        .collect()
    )
    return _finalize_polars_to_obs_pandas(enriched)


def harmonize_buechler_caf_obs_polars(
    raw_obs: pd.DataFrame,
    active_ident: Sequence[str] = (),
    cell_names_fallback: Sequence[str] = (),
) -> pd.DataFrame:
    """Harmonize Buechler et al. (2021, Nature) Human Perturbed-State Fibroblast `.obs` via Polars."""
    pl_df = _pandas_to_polars_preserving_all(raw_obs, cell_names_fallback=cell_names_fallback)
    if pl_df.is_empty():
        return _finalize_polars_to_obs_pandas(pl_df)

    n_rows = len(pl_df)
    ident_vals = [str(v) for v in active_ident] if len(active_ident) == n_rows else ["NA"] * n_rows
    pl_df = pl_df.with_columns(pl.Series("_active_ident", ident_vals, dtype=pl.String))

    clust_name_expr = (
        pl.col("ClustName").cast(pl.String).fill_null("NA")
        if "ClustName" in pl_df.columns
        else pl.lit("NA")
    )
    clust_num_expr = (
        pl.col("Cluster").cast(pl.String).fill_null("NA")
        if "Cluster" in pl_df.columns
        else pl.col("_active_ident")
    )
    tissue_expr = (
        pl.col("Tissue").cast(pl.String).fill_null("Multi-organ")
        if "Tissue" in pl_df.columns
        else pl.lit("Multi-organ")
    )

    subtype_map = {k: v.subtype_label for k, v in BUECHLER_CAF_CLUSTER_SPECS.items()}
    broad_map = {k: v.broad_cell_type for k, v in BUECHLER_CAF_CLUSTER_SPECS.items()}
    comp_map = {k: v.lineage_compartment for k, v in BUECHLER_CAF_CLUSTER_SPECS.items()}
    branch_map = {k: v.lineage_branch for k, v in BUECHLER_CAF_CLUSTER_SPECS.items()}
    func_map = {k: v.functional_state for k, v in BUECHLER_CAF_CLUSTER_SPECS.items()}

    tissue_to_organ: Mapping[str, str] = {
        "NSCLC": "Lung",
        "PDAC": "Pancreas",
        "IPF": "Lung",
        "COVID": "Lung",
        "Ulcerative Colitis": "Colon / Rectum",
    }
    tissue_to_cancer_code: Mapping[str, str] = {
        "NSCLC": "NSCLC",
        "PDAC": "PDAC",
        "IPF": "FIBROSIS_IPF",
        "COVID": "INFLAMMATORY_COVID",
        "Ulcerative Colitis": "INFLAMMATORY_UC",
    }
    tissue_to_cancer_type: Mapping[str, str] = {
        "NSCLC": "Non-Small Cell Lung",
        "PDAC": "Pancreatic Ductal Adenocarcinoma",
        "IPF": "Idiopathic Pulmonary Fibrosis",
        "COVID": "COVID-19 Lung Injury",
        "Ulcerative Colitis": "Ulcerative Colitis",
    }
    tissue_to_sample_type: Mapping[str, str] = {
        "NSCLC": "tumor",
        "PDAC": "tumor",
        "IPF": "fibrotic",
        "COVID": "inflamed",
        "Ulcerative Colitis": "inflamed",
    }

    enriched = (
        pl_df.lazy()
        .with_columns(
            pl.when(~clust_name_expr.is_in(["", "NA", "nan"]))
            .then(clust_name_expr)
            .otherwise(clust_num_expr.replace_strict(_BUECHLER_NUM_TO_NAME, default=clust_num_expr))
            .alias("author_cluster_id"),
            tissue_expr.alias("disease_condition"),
            tissue_expr.replace_strict(tissue_to_organ, default=pl.lit("Multi-organ")).alias("organ"),
            tissue_expr.replace_strict(tissue_to_cancer_code, default=tissue_expr).alias("cancer_code"),
            tissue_expr.replace_strict(tissue_to_cancer_type, default=tissue_expr).alias("cancer_type"),
            tissue_expr.replace_strict(tissue_to_sample_type, default=pl.lit("perturbed")).alias(
                "sample_type"
            ),
            pl.when(tissue_expr.is_in(["NSCLC", "PDAC"]))
            .then(pl.lit("Tumor"))
            .otherwise(tissue_expr)
            .alias("tissue"),
            (pl.lit("Buechler2021_") + tissue_expr.str.replace_all(" ", "_")).alias("source_study"),
            # Extract sample prefix from cell_id when separated by `_`
            pl.col("cell_id").str.split("_").list.first().fill_null("Buechler_Sample").alias("_bar_prefix"),
        )
        .with_columns(
            (pl.col("disease_condition").str.replace_all(" ", "_") + pl.lit("_") + pl.col("_bar_prefix")).alias(
                "sample"
            ),
            pl.col("author_cluster_id")
            .replace_strict(broad_map, default=pl.lit("Fibroblast"))
            .alias("cell_type"),
            pl.col("author_cluster_id")
            .replace_strict(subtype_map, default=pl.col("author_cluster_id"))
            .alias("cell_subtype"),
            pl.col("author_cluster_id")
            .replace_strict(comp_map, default=pl.lit("STROMAL_CAF"))
            .alias("lineage_compartment"),
            pl.col("author_cluster_id")
            .replace_strict(branch_map, default=pl.lit("CAF / Fibroblast"))
            .alias("lineage_branch"),
            pl.col("author_cluster_id")
            .replace_strict(func_map, default=pl.lit("Human perturbed-state fibroblast"))
            .alias("functional_state"),
        )
        .with_columns(
            pl.col("sample").alias("sample_id"),
            pl.col("sample").alias("patient"),
            pl.col("sample").alias("patient_id"),
            pl.col("cell_subtype").alias("cell_state"),
            pl.col("author_cluster_id").alias("cell_type_author"),
            pl.lit(BUECHLER_CAF_DATASET_ID).alias("dataset"),
            pl.lit(BUECHLER_CAF_DATASET_ID).alias("dataset_id"),
        )
        .drop(["_active_ident", "_bar_prefix"])
        .collect()
    )
    return _finalize_polars_to_obs_pandas(enriched)


# =============================================================================
# 3. Imperative Shell: Streaming Helpers, Downloaders, Preprocessors & Loaders
# =============================================================================


def _slice_staged_csr_batch(
    assay_grp: h5py.Group,
    start_row: int,
    end_row: int,
    n_cols: int,
) -> sp.csr_matrix:
    """Slice `[start_row:end_row]` from a staged CSR HDF5 group in O(batch_nnz) RAM."""
    indptr_slice = np.asarray(assay_grp["indptr"][start_row : end_row + 1], dtype=np.int64)
    p0 = int(indptr_slice[0])
    p1 = int(indptr_slice[-1])
    if p1 <= p0:
        return sp.csr_matrix(
            (np.empty(0, dtype=np.float32), (np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32))),
            shape=(end_row - start_row, n_cols),
            dtype=np.float32,
        )
    data_slice = np.asarray(assay_grp["data"][p0:p1], dtype=np.float32)
    indices_slice = np.asarray(assay_grp["indices"][p0:p1], dtype=np.int32)
    local_indptr = indptr_slice - p0
    return sp.csr_matrix(
        (data_slice, indices_slice, local_indptr),
        shape=(end_row - start_row, n_cols),
        dtype=np.float32,
    )


def _download_verified_file_with_sidecar(
    url: str,
    dest_path: Path,
    checksum: ChecksumSpec,
    force_download: bool,
    verify_md5: bool,
) -> Result[Path, str]:
    """Download a remote file with CSCS cluster guard, HTTP 206 resume, and fast MD5 sidecar cache."""
    if is_cscs_environment():
        return Failure(
            f"Network downloads are disabled on CSCS Alps Eiger (cannot fetch {url} -> {dest_path}). "
            f"Pre-stage raw or preprocessed files on the login/transfer node."
        )
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    marker = dest_path.parent / f".{dest_path.name}.md5.ok"

    if force_download:
        dest_path.unlink(missing_ok=True)
        marker.unlink(missing_ok=True)
        dest_path.with_name(dest_path.name + ".part").unlink(missing_ok=True)

    if verify_md5 and dest_path.is_file() and marker.is_file():
        st = dest_path.stat()
        expected_size_ok = (
            st.st_size == checksum.size_bytes.unwrap()
            if isinstance(checksum.size_bytes, Some)
            else st.st_size > 0
        )
        stamp = f"{st.st_size}:{st.st_mtime_ns}:{checksum.expected_hash}"
        try:
            if expected_size_ok and marker.read_text(encoding="utf-8").strip() == stamp:
                return Success(dest_path)
        except OSError:
            pass

    res = download_single_file(
        url=url,
        dest_path=dest_path,
        timeout_sec=180.0,
        max_retries=10,
        resume=True,
        checksum=Some(checksum) if verify_md5 else Nothing,
    )
    match res:
        case Success(verified_path):
            if verify_md5 and verified_path.is_file():
                try:
                    st = verified_path.stat()
                    marker.write_text(
                        f"{st.st_size}:{st.st_mtime_ns}:{checksum.expected_hash}\n",
                        encoding="utf-8",
                    )
                except OSError:
                    pass
            return Success(verified_path)
        case Failure(err):
            return Failure(err)
        case _:
            return Failure("Unexpected download result state")


def _apply_maybe_subset_and_subsample(
    adata: ad.AnnData,
    subset: Maybe[Mapping[str, str]] = Nothing,
    subsample_n: Maybe[int] = Nothing,
) -> ad.AnnData:
    """Apply `Maybe`-typed `subset` and `subsample_n` filters cleanly without `None` signatures."""
    match (subset, subsample_n):
        case (Some(sub_map), Some(n_sub)):
            return _apply_subset_and_subsample(adata, subset=sub_map, subsample_n=n_sub)
        case (Some(sub_map), _):
            return _apply_subset_and_subsample(adata, subset=sub_map)
        case (_, Some(n_sub)):
            return _apply_subset_and_subsample(adata, subsample_n=n_sub)
        case _:
            return adata


# -----------------------------------------------------------------------------
# 3A. Gao et al. (2024, Cancer Cell) — `Gao_CAF_Atlas`
# -----------------------------------------------------------------------------


def download_gao_caf_atlas(
    raw_dir: Maybe[Path] = Nothing,
    force_download: bool = False,
    verify_md5: bool = True,
) -> Result[tuple[Path, Path], str]:
    """Download `Fib_Atlas_expression.rda` and `Fib_Atlas_metadata.rda` from Zenodo 13950610."""
    target_dir = (
        raw_dir.unwrap()
        if isinstance(raw_dir, Some)
        else get_raw_dataset_dir(GAO_CAF_DATASET_ID)
    )
    meta_res = _download_verified_file_with_sidecar(
        url=GAO_CAF_META_URL,
        dest_path=target_dir / GAO_CAF_META_FILENAME,
        checksum=GAO_CAF_META_CHECKSUM,
        force_download=force_download,
        verify_md5=verify_md5,
    )
    match meta_res:
        case Failure(err):
            return Failure(err)
        case Success(meta_path):
            pass
        case _:
            return Failure("Unexpected metadata download result")

    expr_res = _download_verified_file_with_sidecar(
        url=GAO_CAF_EXPR_URL,
        dest_path=target_dir / GAO_CAF_EXPR_FILENAME,
        checksum=GAO_CAF_EXPR_CHECKSUM,
        force_download=force_download,
        verify_md5=verify_md5,
    )
    match expr_res:
        case Failure(err):
            return Failure(err)
        case Success(expr_path):
            return Success((expr_path, meta_path))
        case _:
            return Failure("Unexpected expression download result")


def preprocess_gao_caf_atlas(
    raw_dir: Path,
    output_h5ad: Maybe[Path] = Nothing,
    normalize_ensembl: bool = True,
    ensembl_release: int = 111,
    drop_unmapped: bool = True,
    batch_size: int = 25_000,
) -> Result[Path, str]:
    """Stream `Fib_Atlas_expression.rda` + `Fib_Atlas_metadata.rda` to Ensembl-aligned `.h5ad` (< 1.5 GB RAM)."""
    expr_path = raw_dir / GAO_CAF_EXPR_FILENAME if raw_dir.is_dir() else raw_dir
    meta_path = expr_path.parent / GAO_CAF_META_FILENAME
    if not expr_path.is_file():
        return Failure(f"Gao CAF expression .rda not found at {expr_path}")
    if not meta_path.is_file():
        return Failure(f"Gao CAF metadata .rda not found at {meta_path}")

    target_h5ad = (
        output_h5ad.unwrap().resolve()
        if isinstance(output_h5ad, Some)
        else get_preprocessed_h5ad_path(GAO_CAF_DATASET_ID).resolve()
    )
    target_h5ad.parent.mkdir(parents=True, exist_ok=True)
    tmp_h5ad = target_h5ad.with_name(target_h5ad.name + ".tmp")

    scratch_dir = get_scratch_dataset_dir(GAO_CAF_DATASET_ID)
    scratch_dir.mkdir(parents=True, exist_ok=True)
    staging_h5 = scratch_dir / "gao_caf_staging.h5"

    try:
        meta_df_res = read_r_dataframe(meta_path)
        match meta_df_res:
            case Failure(m_err):
                return Failure(m_err)
            case Success(raw_meta_df):
                pass
            case _:
                return Failure("Unexpected metadata read result")

        stage_res = stream_r_object_to_staging_h5(
            r_path=expr_path,
            staging_h5_path=staging_h5,
            skip_assay_counts=False,
        )
        match stage_res:
            case Failure(s_err):
                return Failure(s_err)
            case Success(manifest):
                pass
            case _:
                return Failure("Unexpected staging result")

        assay_key = manifest.active_assay if manifest.active_assay in manifest.assays else next(
            iter(manifest.assays.keys()), ""
        )
        if not assay_key:
            return Failure(f"No expression matrix found in {expr_path.name}")

        rna_spec: StagedAssaySpec = manifest.assays[assay_key]
        n_cells = rna_spec.n_cells
        raw_genes = (
            list(rna_spec.gene_names)
            if rna_spec.gene_names
            else [f"GENE_{i}" for i in range(rna_spec.n_genes)]
        )

        obs_clean, umap_maybe = harmonize_gao_caf_obs_polars(
            raw_meta_df=raw_meta_df,
            expr_cell_names=rna_spec.cell_names,
        )
        if len(obs_clean) != n_cells:
            return Failure(
                f"Gao CAF cell count mismatch: expression={n_cells}, metadata={len(obs_clean)}"
            )

        raw_var = pd.DataFrame(index=pd.Index(raw_genes, name="original_id"))
        unmapped_genes: list[str] = []
        if normalize_ensembl:
            ensembl_dir = get_ensembl_dir()
            var_final, rna_keep_mask, rna_M, unmapped_genes = _build_ensembl_projection_operator(
                raw_var_names=raw_genes,
                raw_var=raw_var,
                release=ensembl_release,
                species="human",
                ensembl_dir=ensembl_dir,
                drop_unmapped=drop_unmapped,
                aggregation="sum",
            )
            rna_agg = "log1p_sum"
        else:
            unique_idx = ad.utils.make_index_unique(pd.Index(raw_genes, name="gene_id"))
            var_final = pd.DataFrame(
                {"gene_name": raw_genes, "original_id": raw_genes},
                index=unique_idx,
            )
            rna_keep_mask = np.ones(len(raw_genes), dtype=bool)
            rna_M = sp.eye(len(raw_genes), dtype=np.float32, format="csr")
            rna_agg = "sum"

        n_final_vars = len(var_final)
        tmp_h5ad.unlink(missing_ok=True)

        with ExitStack() as stack:
            src_h5 = stack.enter_context(h5py.File(staging_h5, "r"))
            dst_h5 = stack.enter_context(h5py.File(tmp_h5ad, "w"))
            dst_h5.attrs["encoding-type"] = "anndata"
            dst_h5.attrs["encoding-version"] = "0.1.0"

            write_anndata_dataframe_elem(dst_h5, "obs", obs_clean, default_index_name="cell_id")
            write_anndata_dataframe_elem(dst_h5, "var", var_final, default_index_name="gene_id")

            layers_grp = dst_h5.create_group("layers")
            layers_grp.attrs["encoding-type"] = "dict"
            layers_grp.attrs["encoding-version"] = "0.1.0"

            grp_x = init_compressed_csr_group(dst_h5, "X", n_obs=n_cells, n_vars=n_final_vars)
            grp_log1p = init_compressed_csr_group(
                layers_grp,
                "log1p_norm",
                n_obs=n_cells,
                n_vars=n_final_vars,
            )

            src_grp = src_h5[rna_spec.h5_group]
            for start in range(0, n_cells, batch_size):
                end = min(start + batch_size, n_cells)
                raw_chunk = (
                    _slice_staged_csr_batch(src_grp, start, end, rna_spec.n_genes)
                    if rna_spec.is_sparse
                    else sp.csr_matrix(src_grp["dense"][start:end, :], dtype=np.float32)
                )
                proj_chunk = _project_sparse_chunk_to_ensembl(
                    raw_chunk,
                    rna_keep_mask,
                    rna_M,
                    rna_agg,
                )
                append_csr_batch_to_group(grp_x, proj_chunk)
                append_csr_batch_to_group(grp_log1p, proj_chunk)
                gc.collect()

            obsm_grp = dst_h5.create_group("obsm")
            obsm_grp.attrs["encoding-type"] = "dict"
            obsm_grp.attrs["encoding-version"] = "0.1.0"
            match umap_maybe:
                case Some(umap_arr):
                    anndata.io.write_elem(obsm_grp, "X_umap", umap_arr)
                case _:
                    pass

            uns_payload: dict[str, Any] = {
                "dataset_id": GAO_CAF_DATASET_ID,
                "title": "Gao et al. 2024 Cross-Tissue & Pan-Cancer Human Fibroblast Atlas",
                "publication_doi": GAO_CAF_PUBLICATION_DOI,
                "publication_pmid": GAO_CAF_PUBLICATION_PMID,
                "zenodo_record": GAO_CAF_ZENODO_RECORD,
                "zenodo_doi": GAO_CAF_ZENODO_DOI,
                "is_raw_counts": False,
                "has_raw_counts": False,
                "expression_type": "log_normalized",
                "n_unmapped_genes": len(unmapped_genes),
                "atlas_metadata": {
                    "journal": "Cancer Cell",
                    "year": 2024,
                    "original_obs_columns": [str(c) for c in raw_meta_df.columns],
                    "cluster_subtypes": {
                        cid: spec.subtype_label for cid, spec in GAO_CAF_CLUSTER_SPECS.items()
                    },
                    "cluster_functional_states": {
                        cid: spec.functional_state for cid, spec in GAO_CAF_CLUSTER_SPECS.items()
                    },
                },
            }
            anndata.io.write_elem(dst_h5, "uns", uns_payload)

        if target_h5ad.exists():
            target_h5ad.unlink()
        tmp_h5ad.rename(target_h5ad)
        logger.info(
            "Successfully wrote Gao_CAF_Atlas H5AD to %s (%d cells x %d genes)",
            target_h5ad,
            n_cells,
            n_final_vars,
        )
        return Success(target_h5ad)
    except Exception as exc:  # noqa: BLE001
        tmp_h5ad.unlink(missing_ok=True)
        return Failure(f"Failed to preprocess Gao_CAF_Atlas: {exc}")
    finally:
        staging_h5.unlink(missing_ok=True)


def _read_selective_h5ad(h5ad_path: Path, include_all_layers: bool = False) -> ad.AnnData:
    """Read `.h5ad` selectively (`obs`, `var`, `X`, `obsm`, `uns`) to keep RAM minimal."""
    if include_all_layers:
        return ad.read_h5ad(h5ad_path)
    with h5py.File(h5ad_path, "r") as f:
        obs = anndata.io.read_elem(f["obs"])
        var = anndata.io.read_elem(f["var"])
        x_mat = anndata.io.read_elem(f["X"])
        uns = anndata.io.read_elem(f["uns"]) if "uns" in f else {}
        obsm = (
            {k: anndata.io.read_elem(f["obsm"][k]) for k in f["obsm"].keys()}
            if "obsm" in f
            else {}
        )
        return ad.AnnData(X=x_mat, obs=obs, var=var, uns=uns, obsm=obsm)


def load_gao_caf_atlas(
    raw_dir: Maybe[Path] = Nothing,
    output_h5ad: Maybe[Path] = Nothing,
    auto_download: bool = True,
    force_download: bool = False,
    force_recompute: bool = False,
    normalize_ensembl: bool = True,
    ensembl_release: int = 111,
    include_all_layers: bool = False,
    subset: Maybe[Mapping[str, str]] = Nothing,
    subsample_n: Maybe[int] = Nothing,
    repo_root: Maybe[Path] = Nothing,
) -> Result[ad.AnnData, str]:
    """Load the Gao et al. (2024, Cancer Cell) Pan-Cancer / Cross-Tissue Fibroblast Atlas (`Gao_CAF_Atlas`)."""
    root = repo_root.unwrap().resolve() if isinstance(repo_root, Some) else find_repo_root()
    resolved_raw_dir = (
        raw_dir.unwrap().resolve()
        if isinstance(raw_dir, Some)
        else get_raw_dataset_dir(GAO_CAF_DATASET_ID, repo_root=root)
    )
    target_h5ad = (
        output_h5ad.unwrap().resolve()
        if isinstance(output_h5ad, Some)
        else get_preprocessed_h5ad_path(GAO_CAF_DATASET_ID, repo_root=root)
    )

    cached_h5ad: Maybe[Path] = Nothing
    if not force_recompute and not force_download:
        cached_h5ad = (
            Some(target_h5ad)
            if target_h5ad.is_file() and target_h5ad.stat().st_size > 0
            else find_dataset_h5ad(GAO_CAF_DATASET_ID, repo_root=root)
        )

    match cached_h5ad:
        case Some(existing_path):
            try:
                adata = _read_selective_h5ad(existing_path, include_all_layers=include_all_layers)
                return Success(
                    _apply_maybe_subset_and_subsample(adata, subset=subset, subsample_n=subsample_n)
                )
            except Exception as exc:  # noqa: BLE001
                return Failure(f"Failed to read cached {existing_path}: {exc}")
        case _:
            pass

    expr_candidate = resolved_raw_dir / GAO_CAF_EXPR_FILENAME
    meta_candidate = resolved_raw_dir / GAO_CAF_META_FILENAME
    if force_download or not (expr_candidate.is_file() and meta_candidate.is_file()):
        if not auto_download:
            return Failure(
                f"Gao_CAF_Atlas raw files not found in {resolved_raw_dir} and auto_download=False. "
                f"Download from Zenodo {GAO_CAF_ZENODO_RECORD} or set auto_download=True."
            )
        match download_gao_caf_atlas(
            raw_dir=Some(resolved_raw_dir),
            force_download=force_download,
            verify_md5=True,
        ):
            case Failure(dl_err):
                return Failure(dl_err)
            case Success(_):
                pass
            case _:
                return Failure("Unexpected Gao download state")

    match preprocess_gao_caf_atlas(
        raw_dir=resolved_raw_dir,
        output_h5ad=Some(target_h5ad),
        normalize_ensembl=normalize_ensembl,
        ensembl_release=ensembl_release,
    ):
        case Failure(prep_err):
            return Failure(prep_err)
        case Success(written_path):
            try:
                adata = _read_selective_h5ad(written_path, include_all_layers=include_all_layers)
                return Success(
                    _apply_maybe_subset_and_subsample(adata, subset=subset, subsample_n=subsample_n)
                )
            except Exception as exc:  # noqa: BLE001
                return Failure(f"Preprocessed {written_path} succeeded, but load failed: {exc}")
        case _:
            return Failure("Unexpected Gao preprocess state")


# -----------------------------------------------------------------------------
# 3B. Cords et al. (2024, Cancer Cell) — `Cords_CAF_Atlas`
# -----------------------------------------------------------------------------


def download_cords_caf_atlas(
    raw_dir: Maybe[Path] = Nothing,
    force_download: bool = False,
    verify_md5: bool = True,
) -> Result[Path, str]:
    """Download `scRNA-seq_dataobjects.zip` (2.17 GB) from Zenodo 7540604."""
    target_dir = (
        raw_dir.unwrap()
        if isinstance(raw_dir, Some)
        else get_raw_dataset_dir(CORDS_CAF_DATASET_ID)
    )
    return _download_verified_file_with_sidecar(
        url=CORDS_CAF_ZIP_URL,
        dest_path=target_dir / CORDS_CAF_ZIP_FILENAME,
        checksum=CORDS_CAF_ZIP_CHECKSUM,
        force_download=force_download,
        verify_md5=verify_md5,
    )


def _align_cohort_genes_to_canonical_var(
    cohort_genes: Sequence[str],
    canonical_ids: Sequence[str],
    symbol_to_canonical: Mapping[str, str],
) -> tuple[np.ndarray[Any, Any], sp.csr_matrix]:
    """Build sparse projection matrix mapping a Cords cohort's raw gene list onto `canonical_ids`."""
    id_to_col = {gid: idx for idx, gid in enumerate(canonical_ids)}
    n_target = len(canonical_ids)
    keep_flags: list[bool] = []
    target_cols: list[int] = []

    for sym in cohort_genes:
        mapped = symbol_to_canonical.get(sym, sym)
        if mapped in id_to_col:
            keep_flags.append(True)
            target_cols.append(id_to_col[mapped])
        else:
            keep_flags.append(False)

    keep_mask = np.asarray(keep_flags, dtype=bool)
    n_kept = int(keep_mask.sum())
    if n_kept == 0:
        empty_proj = sp.csr_matrix(
            (np.empty(0, dtype=np.float32), (np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32))),
            shape=(0, n_target),
            dtype=np.float32,
        )
        return keep_mask, empty_proj

    proj = sp.csr_matrix(
        (
            np.ones(n_kept, dtype=np.float32),
            (np.arange(n_kept, dtype=np.int32), np.asarray(target_cols, dtype=np.int32)),
        ),
        shape=(n_kept, n_target),
        dtype=np.float32,
    )
    return keep_mask, proj


def preprocess_cords_caf_atlas(
    zip_path_or_dir: Path,
    output_h5ad: Maybe[Path] = Nothing,
    normalize_ensembl: bool = True,
    ensembl_release: int = 111,
    drop_unmapped: bool = True,
    batch_size: int = 10_000,
) -> Result[Path, str]:
    """Stream all 6 Cords et al. (2024) tumor CAF Seurat `.rds` objects into a unified `.h5ad`."""
    zip_path = (
        zip_path_or_dir / CORDS_CAF_ZIP_FILENAME
        if zip_path_or_dir.is_dir()
        else zip_path_or_dir
    )
    if not zip_path.is_file():
        return Failure(f"Cords CAF zip archive not found at {zip_path}")

    target_h5ad = (
        output_h5ad.unwrap().resolve()
        if isinstance(output_h5ad, Some)
        else get_preprocessed_h5ad_path(CORDS_CAF_DATASET_ID).resolve()
    )
    target_h5ad.parent.mkdir(parents=True, exist_ok=True)
    tmp_h5ad = target_h5ad.with_name(target_h5ad.name + ".tmp")

    scratch_dir = get_scratch_dataset_dir(CORDS_CAF_DATASET_ID)
    scratch_dir.mkdir(parents=True, exist_ok=True)

    staged_items: list[tuple[str, str, str, Path, SeuratStagingManifest]] = []
    extracted_rds_paths: list[Path] = []
    staged_h5_paths: list[Path] = []

    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            available_members = set(zf.namelist())
            for member_path, cohort_id, cancer_code, organ in CORDS_COHORT_RDS_MEMBERS:
                if member_path not in available_members:
                    continue
                rds_out = scratch_dir / Path(member_path).name
                h5_out = scratch_dir / f"{cohort_id}_staging.h5"
                extracted_rds_paths.append(rds_out)
                staged_h5_paths.append(h5_out)

                with zf.open(member_path) as src_s, rds_out.open("wb") as dst_s:
                    shutil.copyfileobj(src_s, dst_s)

                match stream_r_object_to_staging_h5(
                    r_path=rds_out,
                    staging_h5_path=h5_out,
                    skip_assay_counts=True,
                ):
                    case Failure(err):
                        return Failure(f"Failed staging {member_path}: {err}")
                    case Success(manifest):
                        staged_items.append((cohort_id, cancer_code, organ, h5_out, manifest))
                    case _:
                        return Failure(f"Unexpected staging state for {member_path}")

                rds_out.unlink(missing_ok=True)

        if not staged_items:
            return Failure(f"No tumor CAF .rds cohorts found inside {zip_path}")

        # Build union gene list across all staged cohorts
        union_genes_ordered: list[str] = []
        seen_raw_genes: set[str] = set()
        for _, _, _, _, man in staged_items:
            if "RNA" not in man.assays:
                return Failure(f"Cohort {man.project_name} lacks 'RNA' assay")
            for g in man.assays["RNA"].gene_names:
                if g not in seen_raw_genes:
                    seen_raw_genes.add(g)
                    union_genes_ordered.append(g)

        raw_union_var = pd.DataFrame(index=pd.Index(union_genes_ordered, name="original_id"))
        symbol_to_canonical: dict[str, str] = {}
        unmapped_genes: list[str] = []

        if normalize_ensembl:
            ensembl_dir = get_ensembl_dir()
            ensembl = ensure_ensembl_release_installed(
                release=ensembl_release,
                species="human",
                ensembl_dir=ensembl_dir,
            )
            cache_path = Path(ensembl_dir) / f"gene_mapping_cache_release_{ensembl_release}.parquet"
            mapping_df = _load_or_update_mapping_cache(union_genes_ordered, ensembl, cache_path)
            for q_sym, gid_val, m_status in zip(
                mapping_df["query"].tolist(),
                mapping_df["gene_id"].tolist(),
                mapping_df["mapping_status"].tolist(),
                strict=False,
            ):
                if gid_val and str(m_status) != "unmapped":
                    symbol_to_canonical[str(q_sym)] = str(gid_val)
                elif not drop_unmapped:
                    symbol_to_canonical[str(q_sym)] = str(q_sym)

            var_final, _, _, unmapped_genes = _build_ensembl_projection_operator(
                raw_var_names=union_genes_ordered,
                raw_var=raw_union_var,
                release=ensembl_release,
                species="human",
                ensembl_dir=ensembl_dir,
                drop_unmapped=drop_unmapped,
                aggregation="sum",
            )
            rna_agg = "log1p_sum"
        else:
            unique_idx = ad.utils.make_index_unique(pd.Index(union_genes_ordered, name="gene_id"))
            var_final = pd.DataFrame(
                {"gene_name": union_genes_ordered, "original_id": union_genes_ordered},
                index=unique_idx,
            )
            symbol_to_canonical = {
                raw_g: str(uid)
                for raw_g, uid in zip(union_genes_ordered, unique_idx.tolist(), strict=False)
            }
            rna_agg = "sum"

        canonical_ids = [str(v) for v in var_final.index.tolist()]
        n_final_vars = len(canonical_ids)

        # Harmonize `.obs` and collect per-cohort UMAP coordinates across all 6 cohorts
        obs_parts: list[pd.DataFrame] = []
        umap_parts: list[np.ndarray[Any, Any]] = []
        for cohort_id, cancer_code, organ, h5_p, man in staged_items:
            rna_spec = man.assays["RNA"]
            cohort_obs = harmonize_cords_caf_obs_polars(
                raw_obs=man.obs_df,
                cohort_id=cohort_id,
                cancer_code=cancer_code,
                organ=organ,
                cell_names_fallback=rna_spec.cell_names,
            )
            with h5py.File(h5_p, "r") as src_f:
                if "umap" in man.reductions and man.reductions["umap"].n_cells == len(cohort_obs):
                    emb_u = np.asarray(src_f[man.reductions["umap"].h5_dataset][:], dtype=np.float32)
                    u_2d = (
                        emb_u[:, :2]
                        if emb_u.shape[1] >= 2
                        else np.pad(emb_u, ((0, 0), (0, 2 - emb_u.shape[1])))
                    )
                else:
                    u_2d = np.zeros((len(cohort_obs), 2), dtype=np.float32)
            cohort_obs["cohort_umap_1"] = u_2d[:, 0]
            cohort_obs["cohort_umap_2"] = u_2d[:, 1]
            obs_parts.append(cohort_obs)
            umap_parts.append(u_2d)

        # Align union of `.obs` columns cleanly via Polars diagonal_relaxed concat
        pl_obs_parts = [_pandas_to_polars_preserving_all(part.drop(columns=["cell_id"], errors="ignore"), cell_names_fallback=part["cell_id"].tolist()) for part in obs_parts]
        merged_pl = pl.concat(pl_obs_parts, how="diagonal_relaxed")
        filled_exprs: list[pl.Expr] = []
        for col_nm, dtype in merged_pl.schema.items():
            if dtype == pl.String:
                filled_exprs.append(pl.col(col_nm).fill_null("NA").alias(col_nm))
            elif dtype == pl.Boolean:
                filled_exprs.append(pl.col(col_nm).fill_null(False).alias(col_nm))
        if filled_exprs:
            merged_pl = merged_pl.with_columns(filled_exprs)

        combined_obs = _finalize_polars_to_obs_pandas(merged_pl)
        total_cells = len(combined_obs)
        combined_umap = np.vstack(umap_parts).astype(np.float32, copy=False)

        tmp_h5ad.unlink(missing_ok=True)
        with h5py.File(tmp_h5ad, "w") as dst_h5:
            dst_h5.attrs["encoding-type"] = "anndata"
            dst_h5.attrs["encoding-version"] = "0.1.0"

            write_anndata_dataframe_elem(dst_h5, "obs", combined_obs, default_index_name="cell_id")
            write_anndata_dataframe_elem(dst_h5, "var", var_final, default_index_name="gene_id")

            layers_grp = dst_h5.create_group("layers")
            layers_grp.attrs["encoding-type"] = "dict"
            layers_grp.attrs["encoding-version"] = "0.1.0"

            grp_x = init_compressed_csr_group(dst_h5, "X", n_obs=total_cells, n_vars=n_final_vars)
            grp_log1p = init_compressed_csr_group(
                layers_grp,
                "log1p_norm",
                n_obs=total_cells,
                n_vars=n_final_vars,
            )

            for _, _, _, h5_p, man in staged_items:
                rna_spec = man.assays["RNA"]
                keep_mask, proj_m = _align_cohort_genes_to_canonical_var(
                    cohort_genes=rna_spec.gene_names,
                    canonical_ids=canonical_ids,
                    symbol_to_canonical=symbol_to_canonical,
                )
                with h5py.File(h5_p, "r") as src_h5:
                    src_grp = src_h5[rna_spec.h5_group]
                    for start in range(0, rna_spec.n_cells, batch_size):
                        end = min(start + batch_size, rna_spec.n_cells)
                        raw_chunk = (
                            _slice_staged_csr_batch(src_grp, start, end, rna_spec.n_genes)
                            if rna_spec.is_sparse
                            else sp.csr_matrix(src_grp["dense"][start:end, :], dtype=np.float32)
                        )
                        proj_chunk = _project_sparse_chunk_to_ensembl(
                            raw_chunk,
                            keep_mask,
                            proj_m,
                            rna_agg,
                        )
                        append_csr_batch_to_group(grp_x, proj_chunk)
                        append_csr_batch_to_group(grp_log1p, proj_chunk)
                        gc.collect()

            obsm_grp = dst_h5.create_group("obsm")
            obsm_grp.attrs["encoding-type"] = "dict"
            obsm_grp.attrs["encoding-version"] = "0.1.0"
            anndata.io.write_elem(obsm_grp, "X_umap", combined_umap)

            uns_payload: dict[str, Any] = {
                "dataset_id": CORDS_CAF_DATASET_ID,
                "title": "Cords et al. 2024 Multi-Cancer Human CAF Classification Atlas",
                "publication_doi": CORDS_CAF_PUBLICATION_DOI,
                "publication_pmid": CORDS_CAF_PUBLICATION_PMID,
                "zenodo_record": CORDS_CAF_ZENODO_RECORD,
                "zenodo_doi": CORDS_CAF_ZENODO_DOI,
                "is_raw_counts": False,
                "has_raw_counts": False,
                "expression_type": "log_normalized",
                "n_unmapped_genes": len(unmapped_genes),
                "atlas_metadata": {
                    "journal": "Cancer Cell / Nature Communications",
                    "year": 2024,
                    "cohorts_integrated": [cid for cid, _, _, _, _ in staged_items],
                    "cluster_subtypes": {
                        cid: spec.subtype_label for cid, spec in CORDS_CAF_CLUSTER_SPECS.items()
                    },
                    "cluster_functional_states": {
                        cid: spec.functional_state for cid, spec in CORDS_CAF_CLUSTER_SPECS.items()
                    },
                },
            }
            anndata.io.write_elem(dst_h5, "uns", uns_payload)

        if target_h5ad.exists():
            target_h5ad.unlink()
        tmp_h5ad.rename(target_h5ad)
        logger.info(
            "Successfully wrote Cords_CAF_Atlas H5AD to %s (%d cells x %d genes)",
            target_h5ad,
            total_cells,
            n_final_vars,
        )
        return Success(target_h5ad)
    except Exception as exc:  # noqa: BLE001
        tmp_h5ad.unlink(missing_ok=True)
        return Failure(f"Failed to preprocess Cords_CAF_Atlas: {exc}")
    finally:
        for p in [*extracted_rds_paths, *staged_h5_paths]:
            p.unlink(missing_ok=True)


def load_cords_caf_atlas(
    raw_dir: Maybe[Path] = Nothing,
    output_h5ad: Maybe[Path] = Nothing,
    auto_download: bool = True,
    force_download: bool = False,
    force_recompute: bool = False,
    normalize_ensembl: bool = True,
    ensembl_release: int = 111,
    include_all_layers: bool = False,
    subset: Maybe[Mapping[str, str]] = Nothing,
    subsample_n: Maybe[int] = Nothing,
    repo_root: Maybe[Path] = Nothing,
) -> Result[ad.AnnData, str]:
    """Load the Cords et al. (2024, Cancer Cell) Multi-Cancer Human CAF Atlas (`Cords_CAF_Atlas`)."""
    root = repo_root.unwrap().resolve() if isinstance(repo_root, Some) else find_repo_root()
    resolved_raw_dir = (
        raw_dir.unwrap().resolve()
        if isinstance(raw_dir, Some)
        else get_raw_dataset_dir(CORDS_CAF_DATASET_ID, repo_root=root)
    )
    target_h5ad = (
        output_h5ad.unwrap().resolve()
        if isinstance(output_h5ad, Some)
        else get_preprocessed_h5ad_path(CORDS_CAF_DATASET_ID, repo_root=root)
    )

    cached_h5ad: Maybe[Path] = Nothing
    if not force_recompute and not force_download:
        cached_h5ad = (
            Some(target_h5ad)
            if target_h5ad.is_file() and target_h5ad.stat().st_size > 0
            else find_dataset_h5ad(CORDS_CAF_DATASET_ID, repo_root=root)
        )

    match cached_h5ad:
        case Some(existing_path):
            try:
                adata = _read_selective_h5ad(existing_path, include_all_layers=include_all_layers)
                return Success(
                    _apply_maybe_subset_and_subsample(adata, subset=subset, subsample_n=subsample_n)
                )
            except Exception as exc:  # noqa: BLE001
                return Failure(f"Failed to read cached {existing_path}: {exc}")
        case _:
            pass

    zip_candidate = resolved_raw_dir / CORDS_CAF_ZIP_FILENAME
    if force_download or not zip_candidate.is_file():
        if not auto_download:
            return Failure(
                f"Cords_CAF_Atlas raw zip not found at {zip_candidate} and auto_download=False. "
                f"Download from Zenodo {CORDS_CAF_ZENODO_RECORD} or set auto_download=True."
            )
        match download_cords_caf_atlas(
            raw_dir=Some(resolved_raw_dir),
            force_download=force_download,
            verify_md5=True,
        ):
            case Failure(dl_err):
                return Failure(dl_err)
            case Success(downloaded_zip):
                zip_candidate = downloaded_zip
            case _:
                return Failure("Unexpected Cords download state")

    match preprocess_cords_caf_atlas(
        zip_path_or_dir=zip_candidate,
        output_h5ad=Some(target_h5ad),
        normalize_ensembl=normalize_ensembl,
        ensembl_release=ensembl_release,
    ):
        case Failure(prep_err):
            return Failure(prep_err)
        case Success(written_path):
            try:
                adata = _read_selective_h5ad(written_path, include_all_layers=include_all_layers)
                return Success(
                    _apply_maybe_subset_and_subsample(adata, subset=subset, subsample_n=subsample_n)
                )
            except Exception as exc:  # noqa: BLE001
                return Failure(f"Preprocessed {written_path} succeeded, but load failed: {exc}")
        case _:
            return Failure("Unexpected Cords preprocess state")


# -----------------------------------------------------------------------------
# 3C. Buechler et al. (2021, Nature) — `Buechler_CAF_Atlas`
# -----------------------------------------------------------------------------


def download_buechler_caf_atlas(
    raw_dir: Maybe[Path] = Nothing,
    force_download: bool = False,
    verify_md5: bool = True,
    repo_root: Maybe[Path] = Nothing,
) -> Result[Path, str]:
    """Download `Human_PS_Fibro.RDS` from FibroXplorer or resolve from `data/manual_download/Buechler_CAF_Atlas/`."""
    root = repo_root.unwrap().resolve() if isinstance(repo_root, Some) else find_repo_root()
    target_dir = (
        raw_dir.unwrap()
        if isinstance(raw_dir, Some)
        else get_raw_dataset_dir(BUECHLER_CAF_DATASET_ID, repo_root=root)
    )
    dest_rds = target_dir / BUECHLER_CAF_RDS_FILENAME

    # Check manual_download staging before hitting FibroXplorer portal
    manual_candidate = (
        get_manual_download_dir(repo_root=root)
        / BUECHLER_CAF_DATASET_ID
        / BUECHLER_CAF_RDS_FILENAME
    )
    if not force_download and not dest_rds.is_file() and manual_candidate.is_file():
        return Success(manual_candidate)

    dl_res = _download_verified_file_with_sidecar(
        url=BUECHLER_CAF_GDRIVE_URL,
        dest_path=dest_rds,
        checksum=BUECHLER_CAF_CHECKSUM,
        force_download=force_download,
        verify_md5=verify_md5,
    )
    match dl_res:
        case Success(verified_path):
            return Success(verified_path)
        case Failure(err):
            return Failure(
                f"Failed to download Buechler_CAF_Atlas from FibroXplorer ({err}). "
                f"If portal-throttled, download '{BUECHLER_CAF_RDS_FILENAME}' from "
                f"https://www.fibroxplorer.com/download and place it at '{manual_candidate}'."
            )
        case _:
            return Failure("Unexpected Buechler download state")


def preprocess_buechler_caf_atlas(
    rds_path_or_dir: Path,
    output_h5ad: Maybe[Path] = Nothing,
    normalize_ensembl: bool = True,
    ensembl_release: int = 111,
    drop_unmapped: bool = True,
    batch_size: int = 10_000,
) -> Result[Path, str]:
    """Stream `Human_PS_Fibro.RDS` into an Ensembl-normalized `.h5ad` preserving PCA, Harmony, and UMAP."""
    rds_path = (
        rds_path_or_dir / BUECHLER_CAF_RDS_FILENAME
        if rds_path_or_dir.is_dir()
        else rds_path_or_dir
    )
    if not rds_path.is_file():
        return Failure(f"Buechler Human_PS_Fibro.RDS not found at {rds_path}")

    target_h5ad = (
        output_h5ad.unwrap().resolve()
        if isinstance(output_h5ad, Some)
        else get_preprocessed_h5ad_path(BUECHLER_CAF_DATASET_ID).resolve()
    )
    target_h5ad.parent.mkdir(parents=True, exist_ok=True)
    tmp_h5ad = target_h5ad.with_name(target_h5ad.name + ".tmp")

    scratch_dir = get_scratch_dataset_dir(BUECHLER_CAF_DATASET_ID)
    scratch_dir.mkdir(parents=True, exist_ok=True)
    staging_h5 = scratch_dir / "buechler_caf_staging.h5"

    try:
        match stream_r_object_to_staging_h5(
            r_path=rds_path,
            staging_h5_path=staging_h5,
            skip_assay_counts=True,
        ):
            case Failure(err):
                return Failure(err)
            case Success(manifest):
                pass
            case _:
                return Failure("Unexpected Buechler staging state")

        assay_key = "RNA" if "RNA" in manifest.assays else next(iter(manifest.assays.keys()), "")
        if not assay_key:
            return Failure(f"No RNA assay found in {rds_path.name}")

        rna_spec = manifest.assays[assay_key]
        n_cells = rna_spec.n_cells
        raw_genes = (
            list(rna_spec.gene_names)
            if rna_spec.gene_names
            else [f"GENE_{i}" for i in range(rna_spec.n_genes)]
        )

        obs_clean = harmonize_buechler_caf_obs_polars(
            raw_obs=manifest.obs_df,
            active_ident=manifest.active_ident,
            cell_names_fallback=rna_spec.cell_names,
        )
        if len(obs_clean) != n_cells:
            return Failure(
                f"Buechler cell count mismatch: RNA={n_cells}, meta.data={len(obs_clean)}"
            )

        raw_var = pd.DataFrame(index=pd.Index(raw_genes, name="original_id"))
        unmapped_genes: list[str] = []
        if normalize_ensembl:
            ensembl_dir = get_ensembl_dir()
            var_final, rna_keep_mask, rna_M, unmapped_genes = _build_ensembl_projection_operator(
                raw_var_names=raw_genes,
                raw_var=raw_var,
                release=ensembl_release,
                species="human",
                ensembl_dir=ensembl_dir,
                drop_unmapped=drop_unmapped,
                aggregation="sum",
            )
            rna_agg = "log1p_sum"
        else:
            unique_idx = ad.utils.make_index_unique(pd.Index(raw_genes, name="gene_id"))
            var_final = pd.DataFrame(
                {"gene_name": raw_genes, "original_id": raw_genes},
                index=unique_idx,
            )
            rna_keep_mask = np.ones(len(raw_genes), dtype=bool)
            rna_M = sp.eye(len(raw_genes), dtype=np.float32, format="csr")
            rna_agg = "sum"

        n_final_vars = len(var_final)
        tmp_h5ad.unlink(missing_ok=True)

        with ExitStack() as stack:
            src_h5 = stack.enter_context(h5py.File(staging_h5, "r"))
            dst_h5 = stack.enter_context(h5py.File(tmp_h5ad, "w"))
            dst_h5.attrs["encoding-type"] = "anndata"
            dst_h5.attrs["encoding-version"] = "0.1.0"

            write_anndata_dataframe_elem(dst_h5, "obs", obs_clean, default_index_name="cell_id")
            write_anndata_dataframe_elem(dst_h5, "var", var_final, default_index_name="gene_id")

            layers_grp = dst_h5.create_group("layers")
            layers_grp.attrs["encoding-type"] = "dict"
            layers_grp.attrs["encoding-version"] = "0.1.0"

            grp_x = init_compressed_csr_group(dst_h5, "X", n_obs=n_cells, n_vars=n_final_vars)
            grp_log1p = init_compressed_csr_group(
                layers_grp,
                "log1p_norm",
                n_obs=n_cells,
                n_vars=n_final_vars,
            )

            src_grp = src_h5[rna_spec.h5_group]
            for start in range(0, n_cells, batch_size):
                end = min(start + batch_size, n_cells)
                raw_chunk = (
                    _slice_staged_csr_batch(src_grp, start, end, rna_spec.n_genes)
                    if rna_spec.is_sparse
                    else sp.csr_matrix(src_grp["dense"][start:end, :], dtype=np.float32)
                )
                proj_chunk = _project_sparse_chunk_to_ensembl(
                    raw_chunk,
                    rna_keep_mask,
                    rna_M,
                    rna_agg,
                )
                append_csr_batch_to_group(grp_x, proj_chunk)
                append_csr_batch_to_group(grp_log1p, proj_chunk)
                gc.collect()

            obsm_grp = dst_h5.create_group("obsm")
            obsm_grp.attrs["encoding-type"] = "dict"
            obsm_grp.attrs["encoding-version"] = "0.1.0"
            for red_key, red_spec in manifest.reductions.items():
                if red_spec.n_cells == n_cells:
                    obsm_name = red_key if red_key.startswith("X_") else f"X_{red_key.lower()}"
                    emb_arr = np.asarray(src_h5[red_spec.h5_dataset][:], dtype=np.float32)
                    anndata.io.write_elem(obsm_grp, obsm_name, emb_arr)

            uns_payload: dict[str, Any] = {
                "dataset_id": BUECHLER_CAF_DATASET_ID,
                "title": "Buechler et al. 2021 Human Perturbed-State Fibroblast Lineage Atlas",
                "publication_doi": BUECHLER_CAF_PUBLICATION_DOI,
                "publication_pmid": BUECHLER_CAF_PUBLICATION_PMID,
                "source_file": BUECHLER_CAF_RDS_FILENAME,
                "is_raw_counts": False,
                "has_raw_counts": False,
                "expression_type": "log_normalized",
                "n_unmapped_genes": len(unmapped_genes),
                "atlas_metadata": {
                    "journal": "Nature",
                    "year": 2021,
                    "original_obs_columns": [str(c) for c in manifest.obs_df.columns],
                    "cluster_subtypes": {
                        cid: spec.subtype_label for cid, spec in BUECHLER_CAF_CLUSTER_SPECS.items()
                    },
                    "cluster_functional_states": {
                        cid: spec.functional_state
                        for cid, spec in BUECHLER_CAF_CLUSTER_SPECS.items()
                    },
                },
            }
            anndata.io.write_elem(dst_h5, "uns", uns_payload)

        if target_h5ad.exists():
            target_h5ad.unlink()
        tmp_h5ad.rename(target_h5ad)
        logger.info(
            "Successfully wrote Buechler_CAF_Atlas H5AD to %s (%d cells x %d genes)",
            target_h5ad,
            n_cells,
            n_final_vars,
        )
        return Success(target_h5ad)
    except Exception as exc:  # noqa: BLE001
        tmp_h5ad.unlink(missing_ok=True)
        return Failure(f"Failed to preprocess Buechler_CAF_Atlas: {exc}")
    finally:
        staging_h5.unlink(missing_ok=True)


def load_buechler_caf_atlas(
    raw_dir: Maybe[Path] = Nothing,
    output_h5ad: Maybe[Path] = Nothing,
    auto_download: bool = True,
    force_download: bool = False,
    force_recompute: bool = False,
    normalize_ensembl: bool = True,
    ensembl_release: int = 111,
    include_all_layers: bool = False,
    subset: Maybe[Mapping[str, str]] = Nothing,
    subsample_n: Maybe[int] = Nothing,
    repo_root: Maybe[Path] = Nothing,
) -> Result[ad.AnnData, str]:
    """Load the Buechler et al. (2021, Nature) Human Perturbed-State Fibroblast Atlas (`Buechler_CAF_Atlas`)."""
    root = repo_root.unwrap().resolve() if isinstance(repo_root, Some) else find_repo_root()
    resolved_raw_dir = (
        raw_dir.unwrap().resolve()
        if isinstance(raw_dir, Some)
        else get_raw_dataset_dir(BUECHLER_CAF_DATASET_ID, repo_root=root)
    )
    target_h5ad = (
        output_h5ad.unwrap().resolve()
        if isinstance(output_h5ad, Some)
        else get_preprocessed_h5ad_path(BUECHLER_CAF_DATASET_ID, repo_root=root)
    )

    cached_h5ad: Maybe[Path] = Nothing
    if not force_recompute and not force_download:
        cached_h5ad = (
            Some(target_h5ad)
            if target_h5ad.is_file() and target_h5ad.stat().st_size > 0
            else find_dataset_h5ad(BUECHLER_CAF_DATASET_ID, repo_root=root)
        )

    match cached_h5ad:
        case Some(existing_path):
            try:
                adata = _read_selective_h5ad(existing_path, include_all_layers=include_all_layers)
                return Success(
                    _apply_maybe_subset_and_subsample(adata, subset=subset, subsample_n=subsample_n)
                )
            except Exception as exc:  # noqa: BLE001
                return Failure(f"Failed to read cached {existing_path}: {exc}")
        case _:
            pass

    manual_rds = (
        get_manual_download_dir(repo_root=root)
        / BUECHLER_CAF_DATASET_ID
        / BUECHLER_CAF_RDS_FILENAME
    )
    raw_rds = resolved_raw_dir / BUECHLER_CAF_RDS_FILENAME
    rds_candidate = raw_rds if raw_rds.is_file() else manual_rds

    if force_download or not rds_candidate.is_file():
        if not auto_download:
            return Failure(
                f"Buechler_CAF_Atlas raw file not found at '{raw_rds}' or '{manual_rds}' "
                f"and auto_download=False. Place '{BUECHLER_CAF_RDS_FILENAME}' in "
                f"'{manual_rds.parent}' or set auto_download=True."
            )
        match download_buechler_caf_atlas(
            raw_dir=Some(resolved_raw_dir),
            force_download=force_download,
            verify_md5=True,
            repo_root=Some(root),
        ):
            case Failure(dl_err):
                return Failure(dl_err)
            case Success(downloaded_rds):
                rds_candidate = downloaded_rds
            case _:
                return Failure("Unexpected Buechler download state")

    match preprocess_buechler_caf_atlas(
        rds_path_or_dir=rds_candidate,
        output_h5ad=Some(target_h5ad),
        normalize_ensembl=normalize_ensembl,
        ensembl_release=ensembl_release,
    ):
        case Failure(prep_err):
            return Failure(prep_err)
        case Success(written_path):
            try:
                adata = _read_selective_h5ad(written_path, include_all_layers=include_all_layers)
                return Success(
                    _apply_maybe_subset_and_subsample(adata, subset=subset, subsample_n=subsample_n)
                )
            except Exception as exc:  # noqa: BLE001
                return Failure(f"Preprocessed {written_path} succeeded, but load failed: {exc}")
        case _:
            return Failure("Unexpected Buechler preprocess state")


def dispatch_caf_atlas_load(
    atlas: CafAtlasChoice,
    raw_dir: Maybe[Path] = Nothing,
    auto_download: bool = True,
    force_download: bool = False,
    subset: Maybe[Mapping[str, str]] = Nothing,
    subsample_n: Maybe[int] = Nothing,
    repo_root: Maybe[Path] = Nothing,
) -> Result[ad.AnnData, str]:
    """Exhaustive `match...case` router across all three CAF atlases."""
    match atlas:
        case CafAtlasChoice.GAO:
            return load_gao_caf_atlas(
                raw_dir=raw_dir,
                auto_download=auto_download,
                force_download=force_download,
                subset=subset,
                subsample_n=subsample_n,
                repo_root=repo_root,
            )
        case CafAtlasChoice.CORDS:
            return load_cords_caf_atlas(
                raw_dir=raw_dir,
                auto_download=auto_download,
                force_download=force_download,
                subset=subset,
                subsample_n=subsample_n,
                repo_root=repo_root,
            )
        case CafAtlasChoice.BUECHLER:
            return load_buechler_caf_atlas(
                raw_dir=raw_dir,
                auto_download=auto_download,
                force_download=force_download,
                subset=subset,
                subsample_n=subsample_n,
                repo_root=repo_root,
            )
        case _ as unreachable:
            assert_never(unreachable)


__all__ = [
    "BUECHLER_CAF_CHECKSUM",
    "BUECHLER_CAF_CLUSTER_SPECS",
    "BUECHLER_CAF_DATASET_ID",
    "BUECHLER_CAF_GDRIVE_URL",
    "BUECHLER_CAF_PUBLICATION_DOI",
    "BUECHLER_CAF_PUBLICATION_PMID",
    "BUECHLER_CAF_RDS_FILENAME",
    "CORDS_CAF_CLUSTER_SPECS",
    "CORDS_CAF_DATASET_ID",
    "CORDS_CAF_PUBLICATION_DOI",
    "CORDS_CAF_PUBLICATION_PMID",
    "CORDS_CAF_ZENODO_DOI",
    "CORDS_CAF_ZENODO_RECORD",
    "CORDS_CAF_ZIP_CHECKSUM",
    "CORDS_CAF_ZIP_FILENAME",
    "CORDS_CAF_ZIP_URL",
    "CafAtlasChoice",
    "CafSubstateMetadata",
    "GAO_CAF_CLUSTER_SPECS",
    "GAO_CAF_DATASET_ID",
    "GAO_CAF_EXPR_CHECKSUM",
    "GAO_CAF_EXPR_FILENAME",
    "GAO_CAF_EXPR_URL",
    "GAO_CAF_META_CHECKSUM",
    "GAO_CAF_META_FILENAME",
    "GAO_CAF_META_URL",
    "GAO_CAF_PUBLICATION_DOI",
    "GAO_CAF_PUBLICATION_PMID",
    "GAO_CAF_ZENODO_DOI",
    "GAO_CAF_ZENODO_RECORD",
    "dispatch_caf_atlas_load",
    "download_buechler_caf_atlas",
    "download_cords_caf_atlas",
    "download_gao_caf_atlas",
    "harmonize_buechler_caf_obs_polars",
    "harmonize_cords_caf_obs_polars",
    "harmonize_gao_caf_obs_polars",
    "load_buechler_caf_atlas",
    "load_cords_caf_atlas",
    "load_gao_caf_atlas",
    "preprocess_buechler_caf_atlas",
    "preprocess_cords_caf_atlas",
    "preprocess_gao_caf_atlas",
]
