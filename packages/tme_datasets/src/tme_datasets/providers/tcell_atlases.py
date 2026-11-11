"""Pan-Cancer Tumor-Infiltrating T-Cell & Multicellular Reference Atlas providers.

Integrates three peer-reviewed pan-cancer T-cell and cross-compartment atlases:
1. `Zheng_TCell_Atlas` — Zheng et al., "Pan-cancer single-cell landscape of tumor-infiltrating T cells",
   Science 374(6574), abe6474 (2021). DOI: 10.1126/science.abe6474 | PMID: 34914499 |
   Zenodo: 5461803 & GEO: GSE156728 (397,810 T cells: 234,550 CD8+ & 163,260 CD4+, 41 meta-clusters).
2. `Chu_TCell_Atlas` — Chu et al., "Pan-cancer T cell atlas links a cellular stress response state
   to immunotherapy resistance", Nature Medicine 29(6), 1550-1562 (2023).
   DOI: 10.1038/s41591-023-02371-y | PMID: 37248301 (308,048 T cells, 39 fine clusters incl. Tstr).
3. `Shi_PanCancer_Atlas` — Shi et al., "Cross-tissue multicellular coordination and its rewiring in cancer",
   Nature 643(8071), 241-251 (2025). DOI: 10.1038/s41586-025-09053-4 | PMID: 40437094 |
   Zenodo: 15554080 (`PanCancer_igt_s9_fine_counts.h5ad`, streaming T/NK & Stromal/CAF slices).

Strictly follows functional Python standards (.agents/rules/code-style-guide.md):
- Pure Polars `.obs` harmonization preserving 100% of original author metadata columns
- Frozen Pydantic models (`ConfigDict(frozen=True)`), exhaustive `match...case` with `assert_never`
- Zero `None` / `Optional` annotations (`returns.maybe.Maybe` & `returns.result.Result`)
"""

from __future__ import annotations

from contextlib import ExitStack
from enum import StrEnum
import gc
import gzip
import io
from pathlib import Path
import shutil
import tarfile
from typing import Any, Final, Mapping, Sequence, assert_never

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
    read_r_dataframe,
    stream_r_object_to_staging_h5,
)
from .atlas_single_cell import _apply_subset_and_subsample

logger = get_logger("providers.tcell_atlases")


# =============================================================================
# 1. Immutable Atlas Constants, Checksums & Rich T-Cell Cluster Ontologies
# =============================================================================


class TCellSubstateMetadata(BaseModel):
    """Immutable biological ontology entry for a pan-cancer T-cell substate."""

    model_config = ConfigDict(frozen=True)

    cluster_id: str
    short_name: str
    broad_cell_type: str
    lineage_compartment: str
    lineage_branch: str
    trajectory_group: str
    functional_state: str
    canonical_markers: tuple[str, ...]


class TCellAtlasChoice(StrEnum):
    """Exhaustive enum of T-cell / multicellular pan-cancer atlases."""

    ZHENG = "Zheng_TCell_Atlas"
    CHU = "Chu_TCell_Atlas"
    SHI = "Shi_PanCancer_Atlas"


class ShiCompartmentFilter(StrEnum):
    """Exhaustive compartment selector for out-of-core streaming of `Shi_PanCancer_Atlas`."""

    T_AND_CAF = "t_and_caf"
    T_NK_ONLY = "t_nk_only"
    STROMAL_CAF_ONLY = "stromal_caf_only"
    ALL_CELLS = "all_cells"


# --- 1A. Zheng et al. (2021, Science) ---
ZHENG_TCELL_DATASET_ID: Final[str] = "Zheng_TCell_Atlas"
ZHENG_TCELL_ZENODO_RECORD: Final[str] = "5461803"
ZHENG_TCELL_ZENODO_DOI: Final[str] = "10.5281/zenodo.5461803"
ZHENG_TCELL_GEO_ACCESSION: Final[str] = "GSE156728"
ZHENG_TCELL_PUBLICATION_DOI: Final[str] = "10.1126/science.abe6474"
ZHENG_TCELL_PUBLICATION_PMID: Final[str] = "34914499"

ZHENG_TCELL_EXPR_TAR_FILENAME: Final[str] = "data.expression.tar.gz"
ZHENG_TCELL_META_TAR_FILENAME: Final[str] = "data.metaInfo.tar.gz"
ZHENG_TCELL_TCR_TAR_FILENAME: Final[str] = "data.tcr.tar.gz"

ZHENG_TCELL_EXPR_URL: Final[str] = (
    "https://zenodo.org/records/5461803/files/data.expression.tar.gz?download=1"
)
ZHENG_TCELL_META_URL: Final[str] = (
    "https://zenodo.org/records/5461803/files/data.metaInfo.tar.gz?download=1"
)
ZHENG_TCELL_TCR_URL: Final[str] = (
    "https://zenodo.org/records/5461803/files/data.tcr.tar.gz?download=1"
)

ZHENG_TCELL_EXPR_CHECKSUM: Final[ChecksumSpec] = ChecksumSpec(
    algorithm="md5",
    expected_hash="e4abd31e6a8dfe8b9b33b6f22b264b3f",
    size_bytes=Some(14_004_090_334),
)
ZHENG_TCELL_META_CHECKSUM: Final[ChecksumSpec] = ChecksumSpec(
    algorithm="md5",
    expected_hash="be18bdf0f2a6ac7cd54d8d7ce64dc2a8",
    size_bytes=Some(60_547_827),
)
ZHENG_TCELL_TCR_CHECKSUM: Final[ChecksumSpec] = ChecksumSpec(
    algorithm="md5",
    expected_hash="63b1633d6d7558ef0ec04503f9a39060",
    size_bytes=Some(16_783_417),
)

# All 41 official meta-clusters (17 CD8+ and 24 CD4+) from Zheng et al. Science 2021 (`name.conversion.txt`)
ZHENG_TCELL_CLUSTER_SPECS: Final[Mapping[str, TCellSubstateMetadata]] = {
    "CD8.c01.Tn.MAL": TCellSubstateMetadata(
        cluster_id="CD8.c01.Tn.MAL",
        short_name="Tn",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Central Memory CD8",
        trajectory_group="Start",
        functional_state="Naive CD8+ T cell (MAL+CCR7+LEF1+TCF7+SELL+)",
        canonical_markers=("MAL", "CCR7", "LEF1", "TCF7", "SELL"),
    ),
    "CD8.c02.Tm.IL7R": TCellSubstateMetadata(
        cluster_id="CD8.c02.Tm.IL7R",
        short_name="IL7R+ Tm",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Central Memory CD8",
        trajectory_group="PathTrm",
        functional_state="IL7R+ memory CD8+ T cell (IL7R+GPR183+FOS+)",
        canonical_markers=("IL7R", "GPR183", "CD69", "ZFP36", "FOS"),
    ),
    "CD8.c03.Tm.RPS12": TCellSubstateMetadata(
        cluster_id="CD8.c03.Tm.RPS12",
        short_name="uncharacterized",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Central Memory CD8",
        trajectory_group="PathO",
        functional_state="Ribosomal-enriched resting CD8+ memory T cell (RPS12+RPL13+LTB+)",
        canonical_markers=("RPS12", "RPL13", "RPL32", "LTB", "IL7R"),
    ),
    "CD8.c04.Tm.CD52": TCellSubstateMetadata(
        cluster_id="CD8.c04.Tm.CD52",
        short_name="ZNF683+CXCR6- Tm",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Tissue-Resident Memory CD8 (ZNF683+)",
        trajectory_group="PathO",
        functional_state="Early HOBIT+ pre-resident memory CD8+ T cell (CD52+ZNF683+CXCR6-)",
        canonical_markers=("CD52", "ZNF683", "S100A4", "TMSB10", "GZMM"),
    ),
    "CD8.c05.Tem.CXCR5": TCellSubstateMetadata(
        cluster_id="CD8.c05.Tem.CXCR5",
        short_name="GZMK+ early Tem",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Effector Memory CD8 (GZMK+)",
        trajectory_group="PathTem",
        functional_state="CXCR5+/GZMK+ early transitional effector-memory CD8+ T cell",
        canonical_markers=("CXCR5", "GZMK", "EOMES", "CMC1", "CD74"),
    ),
    "CD8.c06.Tem.GZMK": TCellSubstateMetadata(
        cluster_id="CD8.c06.Tem.GZMK",
        short_name="GZMK+ Tem",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Effector Memory CD8 (GZMK+)",
        trajectory_group="PathTem",
        functional_state="GZMK+ effector-memory CD8+ T cell (GZMK+CXCR4+CCL4+CCL5+)",
        canonical_markers=("GZMK", "CCL5", "CCL4", "CXCR4", "CMC1", "NKG7"),
    ),
    "CD8.c07.Temra.CX3CR1": TCellSubstateMetadata(
        cluster_id="CD8.c07.Temra.CX3CR1",
        short_name="Temra",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Cytotoxic Temra / NK-like CD8",
        trajectory_group="PathO",
        functional_state="Terminally differentiated cytotoxic Temra CD8+ T cell (CX3CR1+FGFBP2+FCGR3A+PRF1+)",
        canonical_markers=("CX3CR1", "FGFBP2", "FCGR3A", "PRF1", "GZMB", "GNLY", "KLRG1"),
    ),
    "CD8.c08.Tk.TYROBP": TCellSubstateMetadata(
        cluster_id="CD8.c08.Tk.TYROBP",
        short_name="KIR+EOMES+ NK-like",
        broad_cell_type="NK-like CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Cytotoxic Temra / NK-like CD8",
        trajectory_group="PathO",
        functional_state="Innate-like KIR+EOMES+TYROBP+ NK-like CD8+ T cell",
        canonical_markers=("TYROBP", "FCER1G", "KLRC2", "KLRD1", "EOMES", "GNLY"),
    ),
    "CD8.c09.Tk.KIR2DL4": TCellSubstateMetadata(
        cluster_id="CD8.c09.Tk.KIR2DL4",
        short_name="KIR+TXK+ NK-like",
        broad_cell_type="NK-like CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Cytotoxic Temra / NK-like CD8",
        trajectory_group="PathO",
        functional_state="KIR2DL4+TXK+ innate-like cytotoxic CD8+ T cell",
        canonical_markers=("KIR2DL4", "TXK", "KLRC1", "KLRF1", "GZMB", "XCL1"),
    ),
    "CD8.c10.Trm.ZNF683": TCellSubstateMetadata(
        cluster_id="CD8.c10.Trm.ZNF683",
        short_name="ZNF683+CXCR6+ Trm",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Tissue-Resident Memory CD8 (ZNF683+)",
        trajectory_group="PathTrm",
        functional_state="Tissue-resident memory CD8+ T cell (ZNF683+ITGAE+CXCR6+XCL1+)",
        canonical_markers=("ZNF683", "ITGAE", "CXCR6", "XCL1", "XCL2", "CD69"),
    ),
    "CD8.c11.Tex.PDCD1": TCellSubstateMetadata(
        cluster_id="CD8.c11.Tex.PDCD1",
        short_name="GZMK+ Tex",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Progenitor Exhausted CD8 (TCF7+)",
        trajectory_group="PathTem",
        functional_state="GZMK+PDCD1+ pre-exhausted transitional CD8+ T cell along Tem->Tex path",
        canonical_markers=("PDCD1", "GZMK", "LAG3", "TIGIT", "DUSP4", "CD27"),
    ),
    "CD8.c12.Tex.CXCL13": TCellSubstateMetadata(
        cluster_id="CD8.c12.Tex.CXCL13",
        short_name="terminal Tex",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Terminal Exhausted CD8 (CXCL13+)",
        trajectory_group="PathTex",
        functional_state="Tumor-reactive CXCL13+HAVCR2+ENTPD1+LAYN+ terminally exhausted CD8+ T cell",
        canonical_markers=("CXCL13", "HAVCR2", "ENTPD1", "LAYN", "PDCD1", "CTLA4", "TOX"),
    ),
    "CD8.c13.Tex.myl12a": TCellSubstateMetadata(
        cluster_id="CD8.c13.Tex.myl12a",
        short_name="OXPHOS- Tex",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Terminal Exhausted CD8 (CXCL13+)",
        trajectory_group="PathTex",
        functional_state="MYL12A+ metabolically hypo-oxidative exhausted CD8+ T cell",
        canonical_markers=("MYL12A", "PDCD1", "HAVCR2", "TIGIT", "LAG3", "TOX"),
    ),
    "CD8.c14.Tex.TCF7": TCellSubstateMetadata(
        cluster_id="CD8.c14.Tex.TCF7",
        short_name="TCF7+ Tex",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Progenitor Exhausted CD8 (TCF7+)",
        trajectory_group="PathTex",
        functional_state="Stem-like TCF7+PDCD1+XCL1+SLAMF6+ progenitor exhausted CD8+ T cell",
        canonical_markers=("TCF7", "PDCD1", "SLAMF6", "XCL1", "CCR7", "CD200"),
    ),
    "CD8.c15.ISG.IFIT1": TCellSubstateMetadata(
        cluster_id="CD8.c15.ISG.IFIT1",
        short_name="ISG+ CD8+ T",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Interferon-Stimulated T (ISG+)",
        trajectory_group="PathISG",
        functional_state="Type-I interferon-stimulated CD8+ T cell (IFIT1+IFIT3+ISG15+MX1+OAS1+)",
        canonical_markers=("IFIT1", "IFIT3", "ISG15", "MX1", "IFI44L", "STAT1"),
    ),
    "CD8.c16.MAIT.SLC4A10": TCellSubstateMetadata(
        cluster_id="CD8.c16.MAIT.SLC4A10",
        short_name="Tc17",
        broad_cell_type="MAIT / Tc17 Cell",
        lineage_compartment="T_CELL",
        lineage_branch="MAIT / Mucosal Tc17",
        trajectory_group="PathO",
        functional_state="Mucosal-associated invariant / RORC+IL17A+ Tc17 CD8+ T cell (SLC4A10+KLRB1+RORC+)",
        canonical_markers=("SLC4A10", "KLRB1", "RORC", "ZBTB16", "CCR6", "IL23R"),
    ),
    "CD8.c17.Tm.NME1": TCellSubstateMetadata(
        cluster_id="CD8.c17.Tm.NME1",
        short_name="NME1+ T",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Metabolic / Cycling NME1+ T",
        trajectory_group="PathO",
        functional_state="Biosynthetically active / cycling NME1+MKI67+ CD8+ T cell",
        canonical_markers=("NME1", "NME2", "MKI67", "TOP2A", "FABP5", "MYC"),
    ),
    "CD4.c01.Tn.TCF7": TCellSubstateMetadata(
        cluster_id="CD4.c01.Tn.TCF7",
        short_name="Tn",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="PathO",
        functional_state="Naive CCR7+TCF7+LEF1+SELL+ CD4+ T cell",
        canonical_markers=("TCF7", "CCR7", "LEF1", "SELL", "MAL"),
    ),
    "CD4.c02.Tn.PASK": TCellSubstateMetadata(
        cluster_id="CD4.c02.Tn.PASK",
        short_name="CXCR5+ pre-Tfh",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Follicular / Th1 Helper CD4 (CXCL13+/CXCR5+)",
        trajectory_group="PathTfh",
        functional_state="PASK+CXCR5+ early lymphoid pre-Tfh CD4+ progenitor",
        canonical_markers=("PASK", "CXCR5", "CCR7", "TCF7", "KLF2"),
    ),
    "CD4.c03.Tn.ADSL": TCellSubstateMetadata(
        cluster_id="CD4.c03.Tn.ADSL",
        short_name="ADSL+ Tn",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="PathO",
        functional_state="Purine-biosynthetic ADSL+ naive CD4+ T cell",
        canonical_markers=("ADSL", "CCR7", "TCF7", "LEF1", "NPM1"),
    ),
    "CD4.c04.Tn.il7r": TCellSubstateMetadata(
        cluster_id="CD4.c04.Tn.il7r",
        short_name="IL7R- Tn",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="PathO",
        functional_state="IL7R-low naive-like CD4+ T cell",
        canonical_markers=("SOX4", "CCR7", "LEF1", "TCF7", "SATB1"),
    ),
    "CD4.c05.Tm.TNF": TCellSubstateMetadata(
        cluster_id="CD4.c05.Tm.TNF",
        short_name="TNF+ T",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="PathO",
        functional_state="Immediate-early cytokine-producing TNF+IFNG+CD69+ memory CD4+ T cell",
        canonical_markers=("TNF", "CD69", "IFNG", "FOS", "JUN", "CD40LG"),
    ),
    "CD4.c06.Tm.ANXA1": TCellSubstateMetadata(
        cluster_id="CD4.c06.Tm.ANXA1",
        short_name="AREG+ Tm",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="PathO",
        functional_state="ANXA1+AREG+IL7R+ tissue-repair central-memory CD4+ T cell",
        canonical_markers=("ANXA1", "AREG", "IL7R", "GPR183", "LMNA"),
    ),
    "CD4.c07.Tm.ANXA2": TCellSubstateMetadata(
        cluster_id="CD4.c07.Tm.ANXA2",
        short_name="TIMP1+ Tm",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="PathO",
        functional_state="ANXA2+TIMP1+S100A4+ effector-memory CD4+ T cell",
        canonical_markers=("ANXA2", "TIMP1", "S100A4", "S100A11", "ITGB1"),
    ),
    "CD4.c08.Tm.CREM": TCellSubstateMetadata(
        cluster_id="CD4.c08.Tm.CREM",
        short_name="CREM+ Tm",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="PathO",
        functional_state="cAMP-responsive CREM+NR4A2+ memory CD4+ T cell",
        canonical_markers=("CREM", "NR4A2", "YPEL5", "ZNF331", "CXCR4"),
    ),
    "CD4.c09.Tm.CCL5": TCellSubstateMetadata(
        cluster_id="CD4.c09.Tm.CCL5",
        short_name="CCL5+ Tm",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="PathO",
        functional_state="CCL5+GZMA+ cytotoxic-memory CD4+ T cell",
        canonical_markers=("CCL5", "GZMA", "CXCR6", "ID2", "BHLHE40"),
    ),
    "CD4.c10.Tm.CAPG": TCellSubstateMetadata(
        cluster_id="CD4.c10.Tm.CAPG",
        short_name="CAPG+ Tm",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="PathO",
        functional_state="CAPG+SH3BGRL3+ cytoskeletal-remodeling memory CD4+ T cell",
        canonical_markers=("CAPG", "SH3BGRL3", "S100A4", "VIM", "COTL1"),
    ),
    "CD4.c11.Tm.GZMA": TCellSubstateMetadata(
        cluster_id="CD4.c11.Tm.GZMA",
        short_name="CAPG+CREM- Tm",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="PathO",
        functional_state="GZMA+CAPG+CREM- effector-memory CD4+ T cell",
        canonical_markers=("GZMA", "CAPG", "CCL5", "KLRB1", "CD40LG"),
    ),
    "CD4.c12.Tem.GZMK": TCellSubstateMetadata(
        cluster_id="CD4.c12.Tem.GZMK",
        short_name="GZMK+ Tem",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="PathO",
        functional_state="GZMK+EOMES+ effector-memory CD4+ T cell",
        canonical_markers=("GZMK", "EOMES", "CCL4", "CCL5", "GZMA"),
    ),
    "CD4.c13.Temra.CX3CR1": TCellSubstateMetadata(
        cluster_id="CD4.c13.Temra.CX3CR1",
        short_name="Temra",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Cytotoxic CD4 (Temra)",
        trajectory_group="PathO",
        functional_state="Cytotoxic CX3CR1+PRF1+GZMB+GNLY+ CD4+ CTL / Temra cell",
        canonical_markers=("CX3CR1", "PRF1", "GZMB", "GNLY", "FGFBP2", "NKG7"),
    ),
    "CD4.c14.Th17.SLC4A10": TCellSubstateMetadata(
        cluster_id="CD4.c14.Th17.SLC4A10",
        short_name="CCR6+ Th17",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Th17 Helper CD4",
        trajectory_group="PathO",
        functional_state="CCR6+KLRB1+RORC+ early Th17 helper CD4+ T cell",
        canonical_markers=("SLC4A10", "CCR6", "KLRB1", "RORC", "CTSSH", "IL7R"),
    ),
    "CD4.c15.Th17.IL23R": TCellSubstateMetadata(
        cluster_id="CD4.c15.Th17.IL23R",
        short_name="IL26+ Th17",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Th17 Helper CD4",
        trajectory_group="PathO",
        functional_state="Pathogenic IL23R+IL17A+IL26+RORC+ Th17 helper CD4+ T cell",
        canonical_markers=("IL23R", "IL17A", "IL17F", "IL26", "RORC", "CCL20"),
    ),
    "CD4.c16.Tfh.CXCR5": TCellSubstateMetadata(
        cluster_id="CD4.c16.Tfh.CXCR5",
        short_name="IL21+ Tfh",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Follicular / Th1 Helper CD4 (CXCL13+/CXCR5+)",
        trajectory_group="PathTfh",
        functional_state="TLS-organizing CXCR5+BCL6+IL21+CXCL13+PDCD1+ follicular helper CD4+ T cell",
        canonical_markers=("CXCR5", "BCL6", "IL21", "CXCL13", "PDCD1", "ICOS", "TOX2"),
    ),
    "CD4.c17.TfhTh1.CXCL13": TCellSubstateMetadata(
        cluster_id="CD4.c17.TfhTh1.CXCL13",
        short_name="IFNG+ Tfh/Th1",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Follicular / Th1 Helper CD4 (CXCL13+/CXCR5+)",
        trajectory_group="PathTfh",
        functional_state="Tumor-reactive CXCL13+IFNG+PDCD1+BHLHE40+ Tfh/Th1 dual helper CD4+ T cell",
        canonical_markers=("CXCL13", "IFNG", "PDCD1", "BHLHE40", "TIGIT", "CTLA4", "TOX"),
    ),
    "CD4.c18.Treg.RTKN2": TCellSubstateMetadata(
        cluster_id="CD4.c18.Treg.RTKN2",
        short_name="TNFRSF9- Treg",
        broad_cell_type="Regulatory T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Resting / S1PR1+ Treg",
        trajectory_group="PathTreg",
        functional_state="RTKN2+FOXP3+IL2RA+TNFRSF9- resting Treg",
        canonical_markers=("RTKN2", "FOXP3", "IL2RA", "IKZF2", "CTLA4"),
    ),
    "CD4.c19.Treg.S1PR1": TCellSubstateMetadata(
        cluster_id="CD4.c19.Treg.S1PR1",
        short_name="S1PR1+ Treg",
        broad_cell_type="Regulatory T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Resting / S1PR1+ Treg",
        trajectory_group="PathTreg",
        functional_state="Recirculating S1PR1+CCR7+KLF2+FOXP3+ naive/resting Treg",
        canonical_markers=("S1PR1", "FOXP3", "CCR7", "KLF2", "SELL", "LEF1"),
    ),
    "CD4.c20.Treg.TNFRSF9": TCellSubstateMetadata(
        cluster_id="CD4.c20.Treg.TNFRSF9",
        short_name="TNFRSF9+ Treg",
        broad_cell_type="Regulatory T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Suppressive Treg (FOXP3+/TNFRSF9+)",
        trajectory_group="PathTreg",
        functional_state="Intratumoral activated immunosuppressive TNFRSF9+CCR8+LAYN+CTLA4+FOXP3+ effector Treg",
        canonical_markers=("TNFRSF9", "FOXP3", "CTLA4", "CCR8", "LAYN", "IL2RA", "BATF"),
    ),
    "CD4.c21.Treg.OAS1": TCellSubstateMetadata(
        cluster_id="CD4.c21.Treg.OAS1",
        short_name="ISG+ Treg",
        broad_cell_type="Regulatory T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Suppressive Treg (FOXP3+/TNFRSF9+)",
        trajectory_group="PathTreg",
        functional_state="Interferon-stimulated OAS1+ISG15+IFIT1+FOXP3+ Treg",
        canonical_markers=("OAS1", "ISG15", "IFIT1", "IFIT3", "FOXP3", "IL2RA"),
    ),
    "CD4.c22.ISG.IFIT1": TCellSubstateMetadata(
        cluster_id="CD4.c22.ISG.IFIT1",
        short_name="ISG+ Th",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Interferon-Stimulated T (ISG+)",
        trajectory_group="PathO",
        functional_state="Type-I interferon-stimulated IFIT1+ISG15+MX1+ helper CD4+ T cell",
        canonical_markers=("IFIT1", "IFIT3", "ISG15", "MX1", "OAS1", "STAT1"),
    ),
    "CD4.c23.Mix.NME1": TCellSubstateMetadata(
        cluster_id="CD4.c23.Mix.NME1",
        short_name="NME1+CCR4- T",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Metabolic / Cycling NME1+ T",
        trajectory_group="PathO",
        functional_state="Metabolically active / proliferating NME1+MKI67+CCR4- CD4+ T cell",
        canonical_markers=("NME1", "NME2", "MKI67", "TOP2A", "TYMS"),
    ),
    "CD4.c24.Mix.NME2": TCellSubstateMetadata(
        cluster_id="CD4.c24.Mix.NME2",
        short_name="NME1+CCR4+ T",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Metabolic / Cycling NME1+ T",
        trajectory_group="PathO",
        functional_state="Biosynthetic NME1+NME2+CCR4+ activated CD4+ T cell",
        canonical_markers=("NME2", "NME1", "CCR4", "FABP5", "NPM1"),
    ),
}


# --- 1B. Chu et al. (2023, Nature Medicine — Pan-Cancer T-Cell Stress Response Atlas) ---
CHU_TCELL_DATASET_ID: Final[str] = "Chu_TCell_Atlas"
CHU_TCELL_PUBLICATION_DOI: Final[str] = "10.1038/s41591-023-02371-y"
CHU_TCELL_PUBLICATION_PMID: Final[str] = "37248301"
CHU_TCELL_PORTAL_URL: Final[str] = "https://singlecell.mdanderson.org/TCM/"
CHU_TCELL_CANDIDATE_URLS: Final[tuple[str, ...]] = (
    "https://singlecell.mdanderson.org/TCM/data/panC_T_all.rds",
    "https://singlecell.mdanderson.org/TCM/download/panC_T_all.rds",
)
CHU_TCELL_CANDIDATE_FILENAMES: Final[tuple[str, ...]] = (
    "panC_T_all.rds",
    "Chu_TCell_Atlas.rds",
    "T_all.rds",
    "CD8.rds",
    "CD4.rds",
)

# Official 15 heat-shock stress-response genes (`HSMarkers` in `Coolgenome/TCM`)
CHU_TCELL_HS_MARKERS: Final[tuple[str, ...]] = (
    "HSPA6",
    "HSPA1A",
    "HSPA1B",
    "DNAJB1",
    "HSPH1",
    "HSP90AA1",
    "HSPE1",
    "HSPB1",
    "BAG3",
    "HSPD1",
    "DNAJA1",
    "HSP90AB1",
    "HSPA8",
    "DNAJB4",
    "DNAJA4",
)

# All 39 fine T-cell states from Chu et al. 2023 (`Coolgenome/TCM`: 14 CD8 + 12 CD4 + 13 prolif/unconventional)
CHU_TCELL_CLUSTER_SPECS: Final[Mapping[str, TCellSubstateMetadata]] = {
    "CD8_c0_t-Teff": TCellSubstateMetadata(
        cluster_id="CD8_c0_t-Teff",
        short_name="CD8_t-Teff",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Effector Memory CD8 (GZMK+)",
        trajectory_group="CD8_Effector",
        functional_state="Transitional GZMK+NKG7+CST7+ effector CD8+ T cell",
        canonical_markers=("GZMK", "NKG7", "CST7", "CCL5", "GZMA", "CD8A"),
    ),
    "CD8_c1_Tex": TCellSubstateMetadata(
        cluster_id="CD8_c1_Tex",
        short_name="CD8_Tex",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Terminal Exhausted CD8 (CXCL13+)",
        trajectory_group="CD8_Exhausted",
        functional_state="Terminally exhausted CXCL13+HAVCR2+LAG3+PDCD1+ENTPD1+TOX+ CD8+ T cell",
        canonical_markers=("CXCL13", "HAVCR2", "PDCD1", "LAG3", "TIGIT", "ENTPD1", "TOX", "LAYN"),
    ),
    "CD8_c2_Teff": TCellSubstateMetadata(
        cluster_id="CD8_c2_Teff",
        short_name="CD8_Teff",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Cytotoxic Temra / NK-like CD8",
        trajectory_group="CD8_Effector",
        functional_state="Cytotoxic CX3CR1+FGFBP2+FCGR3A+PRF1+GZMB+ effector CD8+ T cell",
        canonical_markers=("CX3CR1", "FGFBP2", "FCGR3A", "PRF1", "GZMB", "GNLY", "KLRG1"),
    ),
    "CD8_c3_Tn": TCellSubstateMetadata(
        cluster_id="CD8_c3_Tn",
        short_name="CD8_Tn",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Central Memory CD8",
        trajectory_group="CD8_Naive",
        functional_state="Naive CCR7+LEF1+SELL+TCF7+ CD8+ T cell",
        canonical_markers=("CCR7", "LEF1", "SELL", "TCF7", "MAL"),
    ),
    "CD8_c4_Tstr": TCellSubstateMetadata(
        cluster_id="CD8_c4_Tstr",
        short_name="CD8_Tstr",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Stress-Response T (Tstr)",
        trajectory_group="CD8_Stress",
        functional_state="Heat-shock stress-response CD8+ Tstr state linked to ICB non-response (HSPA1A+HSPA1B+DNAJB1+)",
        canonical_markers=("HSPA1A", "HSPA1B", "HSPA6", "DNAJB1", "HSP90AA1", "HSPH1", "BAG3", "FOS"),
    ),
    "CD8_c5_Tisg": TCellSubstateMetadata(
        cluster_id="CD8_c5_Tisg",
        short_name="CD8_Tisg",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Interferon-Stimulated T (ISG+)",
        trajectory_group="CD8_ISG",
        functional_state="Interferon-stimulated ISG15+IFIT1+IFIT3+MX1+ CD8+ T cell",
        canonical_markers=("ISG15", "IFIT1", "IFIT3", "MX1", "IFI6", "OAS1", "STAT1"),
    ),
    "CD8_c6_Tcm": TCellSubstateMetadata(
        cluster_id="CD8_c6_Tcm",
        short_name="CD8_Tcm",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Central Memory CD8",
        trajectory_group="CD8_Memory",
        functional_state="Central-memory IL7R+GPR183+TCF7+CD28+ CD8+ T cell",
        canonical_markers=("IL7R", "GPR183", "CD28", "TCF7", "LTB", "CCR7"),
    ),
    "CD8_c7_p-Tex": TCellSubstateMetadata(
        cluster_id="CD8_c7_p-Tex",
        short_name="CD8_p-Tex",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Progenitor Exhausted CD8 (TCF7+)",
        trajectory_group="CD8_Exhausted",
        functional_state="Progenitor exhausted TCF7+SLAMF6+PDCD1+XCL1+ CD8+ T cell",
        canonical_markers=("TCF7", "SLAMF6", "PDCD1", "XCL1", "CD200", "TOX"),
    ),
    "CD8_c8_Teff_KLRG1": TCellSubstateMetadata(
        cluster_id="CD8_c8_Teff_KLRG1",
        short_name="CD8_Teff_KLRG1",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Cytotoxic Temra / NK-like CD8",
        trajectory_group="CD8_Effector",
        functional_state="Short-lived KLRG1+GZMH+ZEB2+ effector CD8+ T cell",
        canonical_markers=("KLRG1", "GZMH", "ZEB2", "NKG7", "CCL5", "PRF1"),
    ),
    "CD8_c9_Trm": TCellSubstateMetadata(
        cluster_id="CD8_c9_Trm",
        short_name="CD8_Trm",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Tissue-Resident Memory CD8 (ZNF683+)",
        trajectory_group="CD8_Resident",
        functional_state="Tissue-resident memory ZNF683+ITGAE+CXCR6+CD69+ CD8+ T cell",
        canonical_markers=("ZNF683", "ITGAE", "CXCR6", "CD69", "KLRC1", "XCL1"),
    ),
    "CD8_c10_Teff_CD244": TCellSubstateMetadata(
        cluster_id="CD8_c10_Teff_CD244",
        short_name="CD8_Teff_CD244",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Cytotoxic Temra / NK-like CD8",
        trajectory_group="CD8_Effector",
        functional_state="CD244+EOMES+GZMB+ senescent-like cytotoxic CD8+ T cell",
        canonical_markers=("CD244", "EOMES", "GZMB", "PRF1", "KLRD1", "CCL4"),
    ),
    "CD8_c11_Teff_SEMA4A": TCellSubstateMetadata(
        cluster_id="CD8_c11_Teff_SEMA4A",
        short_name="CD8_Teff_SEMA4A",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Effector Memory CD8 (GZMK+)",
        trajectory_group="CD8_Effector",
        functional_state="SEMA4A+IFNG+TNF+ cytokine-effector CD8+ T cell",
        canonical_markers=("SEMA4A", "IFNG", "TNF", "CCL4", "NR4A2", "CD69"),
    ),
    "CD8_c12_Trm_KLRC1": TCellSubstateMetadata(
        cluster_id="CD8_c12_Trm_KLRC1",
        short_name="CD8_Trm_KLRC1",
        broad_cell_type="NK-like CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Tissue-Resident Memory CD8 (ZNF683+)",
        trajectory_group="CD8_Resident",
        functional_state="NKG2A+ (KLRC1+ITGAE+ZNF683+) intraepithelial resident CD8+ T cell",
        canonical_markers=("KLRC1", "ITGAE", "ZNF683", "KLRD1", "GNLY", "XCL2"),
    ),
    "CD8_c13_Tgd": TCellSubstateMetadata(
        cluster_id="CD8_c13_Tgd",
        short_name="CD8_Tgd",
        broad_cell_type="Gamma-Delta T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Innate / Gamma-Delta T",
        trajectory_group="Unconventional",
        functional_state="TRDC+TRGC1+TRGC2+ cytotoxic gamma-delta T cell",
        canonical_markers=("TRDC", "TRGC1", "TRGC2", "KLRC2", "GNLY", "FCER1G"),
    ),
    "CD4_c0_Tcm": TCellSubstateMetadata(
        cluster_id="CD4_c0_Tcm",
        short_name="CD4_Tcm",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="CD4_Memory",
        functional_state="IL7R+ANXA1+GPR183+ central-memory CD4+ T cell",
        canonical_markers=("IL7R", "ANXA1", "GPR183", "CD40LG", "FOS", "LTB"),
    ),
    "CD4_c1_Treg": TCellSubstateMetadata(
        cluster_id="CD4_c1_Treg",
        short_name="CD4_eTreg",
        broad_cell_type="Regulatory T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Suppressive Treg (FOXP3+/TNFRSF9+)",
        trajectory_group="CD4_Treg",
        functional_state="Activated intratumoral FOXP3+IL2RA+TNFRSF9+CTLA4+CCR8+ effector Treg",
        canonical_markers=("FOXP3", "IL2RA", "TNFRSF9", "CTLA4", "CCR8", "TIGIT", "BATF"),
    ),
    "CD4_c2_Tn": TCellSubstateMetadata(
        cluster_id="CD4_c2_Tn",
        short_name="CD4_Tn",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="CD4_Naive",
        functional_state="Naive CCR7+TCF7+LEF1+SELL+ CD4+ T cell",
        canonical_markers=("CCR7", "TCF7", "LEF1", "SELL", "MAL"),
    ),
    "CD4_c3_Tfh": TCellSubstateMetadata(
        cluster_id="CD4_c3_Tfh",
        short_name="CD4_Tfh",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Follicular / Th1 Helper CD4 (CXCL13+/CXCR5+)",
        trajectory_group="CD4_Tfh",
        functional_state="CXCL13+CXCR5+PDCD1+BCL6+ICOS+IL21+ follicular helper CD4+ T cell",
        canonical_markers=("CXCL13", "CXCR5", "PDCD1", "BCL6", "ICOS", "IL21", "TOX2"),
    ),
    "CD4_c4_Tstr": TCellSubstateMetadata(
        cluster_id="CD4_c4_Tstr",
        short_name="CD4_Tstr",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Stress-Response T (Tstr)",
        trajectory_group="CD4_Stress",
        functional_state="Heat-shock stress-response CD4+ Tstr state linked to ICB resistance (HSPA1A+HSPA1B+DNAJB1+)",
        canonical_markers=("HSPA1A", "HSPA1B", "HSPA6", "DNAJB1", "HSP90AA1", "HSPH1", "BAG3", "HSPD1"),
    ),
    "CD4_c5_Tctl": TCellSubstateMetadata(
        cluster_id="CD4_c5_Tctl",
        short_name="CD4_Tctl",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Cytotoxic CD4 (Temra)",
        trajectory_group="CD4_Effector",
        functional_state="Cytotoxic GZMB+PRF1+GNLY+NKG7+CX3CR1+ CD4+ T lymphocyte",
        canonical_markers=("GZMB", "PRF1", "GNLY", "NKG7", "CX3CR1", "FGFBP2"),
    ),
    "CD4_c6_nr-Treg": TCellSubstateMetadata(
        cluster_id="CD4_c6_nr-Treg",
        short_name="CD4_nr-Treg",
        broad_cell_type="Regulatory T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Resting / S1PR1+ Treg",
        trajectory_group="CD4_Treg",
        functional_state="Naive / resting FOXP3+CCR7+S1PR1+RTKN2+TNFRSF9- regulatory T cell",
        canonical_markers=("FOXP3", "IL2RA", "RTKN2", "CCR7", "S1PR1", "IKZF2"),
    ),
    "CD4_c7_Tem": TCellSubstateMetadata(
        cluster_id="CD4_c7_Tem",
        short_name="CD4_Tem",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="CD4_Memory",
        functional_state="GZMK+CCL5+GZMA+CXCR6+ effector-memory CD4+ T cell",
        canonical_markers=("GZMK", "GZMA", "CCL5", "CXCR6", "KLRB1"),
    ),
    "CD4_c8_Th17": TCellSubstateMetadata(
        cluster_id="CD4_c8_Th17",
        short_name="CD4_Th17",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Th17 Helper CD4",
        trajectory_group="CD4_Helper",
        functional_state="RORC+IL17A+IL23R+CCR6+CCL20+ Th17 helper CD4+ T cell",
        canonical_markers=("IL17A", "IL17F", "RORC", "CCR6", "IL23R", "CCL20", "KLRB1"),
    ),
    "CD4_c9_Tm_CREM": TCellSubstateMetadata(
        cluster_id="CD4_c9_Tm_CREM",
        short_name="CD4_Tm_CREM",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="CD4_Memory",
        functional_state="CREM+NR4A2+AREG+ cAMP-modulated tissue memory CD4+ T cell",
        canonical_markers=("CREM", "NR4A2", "AREG", "TNFAIP3", "CXCR4"),
    ),
    "CD4_c10_Th1_CXCL13": TCellSubstateMetadata(
        cluster_id="CD4_c10_Th1_CXCL13",
        short_name="CD4_Th1_CXCL13",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Follicular / Th1 Helper CD4 (CXCL13+/CXCR5+)",
        trajectory_group="CD4_Tfh",
        functional_state="CXCL13+IFNG+BHLHE40+PDCD1+ Th1-polarized tumor-reactive CD4+ T cell",
        canonical_markers=("CXCL13", "IFNG", "BHLHE40", "PDCD1", "CTLA4", "TIGIT"),
    ),
    "CD4_c11_Tisg": TCellSubstateMetadata(
        cluster_id="CD4_c11_Tisg",
        short_name="CD4_Tisg",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Interferon-Stimulated T (ISG+)",
        trajectory_group="CD4_ISG",
        functional_state="Type-I interferon-stimulated ISG15+IFIT1+IFIT3+MX1+ CD4+ T cell",
        canonical_markers=("ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "IFI44L", "STAT1"),
    ),
    "T_c0_Prolif_MKI67": TCellSubstateMetadata(
        cluster_id="T_c0_Prolif_MKI67",
        short_name="Tprolif_MKI67",
        broad_cell_type="Proliferating T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Metabolic / Cycling NME1+ T",
        trajectory_group="Cycling",
        functional_state="G2/M-phase proliferating MKI67+TOP2A+BIRC5+CDK1+ T cell",
        canonical_markers=("MKI67", "TOP2A", "BIRC5", "UBE2C", "CENPF"),
    ),
    "T_c1_Prolif_STMN1": TCellSubstateMetadata(
        cluster_id="T_c1_Prolif_STMN1",
        short_name="Tprolif_STMN1",
        broad_cell_type="Proliferating T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Metabolic / Cycling NME1+ T",
        trajectory_group="Cycling",
        functional_state="S-phase DNA-replicating STMN1+TYMS+PCNA+MCM5+ T cell",
        canonical_markers=("STMN1", "TYMS", "PCNA", "MCM5", "DUT", "TUBA1B"),
    ),
    "T_c2_Prolif_Tex": TCellSubstateMetadata(
        cluster_id="T_c2_Prolif_Tex",
        short_name="Tprolif_Tex",
        broad_cell_type="Proliferating T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Terminal Exhausted CD8 (CXCL13+)",
        trajectory_group="Cycling",
        functional_state="Cycling exhausted CXCL13+HAVCR2+MKI67+TOP2A+ CD8+ T cell",
        canonical_markers=("MKI67", "TOP2A", "CXCL13", "HAVCR2", "PDCD1", "LAG3"),
    ),
    "T_c3_Prolif_Treg": TCellSubstateMetadata(
        cluster_id="T_c3_Prolif_Treg",
        short_name="Tprolif_Treg",
        broad_cell_type="Regulatory T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Suppressive Treg (FOXP3+/TNFRSF9+)",
        trajectory_group="Cycling",
        functional_state="Proliferating intratumoral FOXP3+TNFRSF9+MKI67+TOP2A+ Treg",
        canonical_markers=("FOXP3", "IL2RA", "TNFRSF9", "MKI67", "TOP2A", "TYMS"),
    ),
    "T_c4_MAIT": TCellSubstateMetadata(
        cluster_id="T_c4_MAIT",
        short_name="MAIT",
        broad_cell_type="MAIT / Tc17 Cell",
        lineage_compartment="T_CELL",
        lineage_branch="MAIT / Mucosal Tc17",
        trajectory_group="Unconventional",
        functional_state="Mucosal-associated invariant SLC4A10+KLRB1+ZBTB16+TRAV1-2+ T cell",
        canonical_markers=("SLC4A10", "KLRB1", "ZBTB16", "RORC", "DPP4", "NCR3"),
    ),
    "T_c5_NKT_FCGR3A": TCellSubstateMetadata(
        cluster_id="T_c5_NKT_FCGR3A",
        short_name="NKT_FCGR3A",
        broad_cell_type="NK-like CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Cytotoxic Temra / NK-like CD8",
        trajectory_group="Unconventional",
        functional_state="FCGR3A+TYROBP+NCAM1+PRF1+ natural killer T cell",
        canonical_markers=("FCGR3A", "TYROBP", "KLRD1", "GNLY", "PRF1", "GZMB"),
    ),
    "T_c6_NKT_KIR": TCellSubstateMetadata(
        cluster_id="T_c6_NKT_KIR",
        short_name="NKT_KIR",
        broad_cell_type="NK-like CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Cytotoxic Temra / NK-like CD8",
        trajectory_group="Unconventional",
        functional_state="KIR2DL3+KIR3DL1+KLRC2+ innate cytotoxic NKT cell",
        canonical_markers=("KIR2DL3", "KIR2DL4", "KLRC2", "KLRF1", "TYROBP", "GZMB"),
    ),
    "T_c7_Tgd_TRDV1": TCellSubstateMetadata(
        cluster_id="T_c7_Tgd_TRDV1",
        short_name="Tgd_TRDV1",
        broad_cell_type="Gamma-Delta T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Innate / Gamma-Delta T",
        trajectory_group="Unconventional",
        functional_state="Tissue-resident Vdelta1+ (TRDV1+TRDC+KLRC1+) gamma-delta T cell",
        canonical_markers=("TRDV1", "TRDC", "TRGC2", "KLRC1", "ITGAE", "XCL1"),
    ),
    "T_c8_Tgd_TRDV2": TCellSubstateMetadata(
        cluster_id="T_c8_Tgd_TRDV2",
        short_name="Tgd_TRDV2",
        broad_cell_type="Gamma-Delta T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Innate / Gamma-Delta T",
        trajectory_group="Unconventional",
        functional_state="Circulating cytotoxic Vgamma9Vdelta2+ (TRDV2+TRGV9+GZMB+) gamma-delta T cell",
        canonical_markers=("TRDV2", "TRGV9", "TRDC", "GZMB", "PRF1", "NKG7"),
    ),
    "T_c9_DN_SOX4": TCellSubstateMetadata(
        cluster_id="T_c9_DN_SOX4",
        short_name="DN_SOX4",
        broad_cell_type="Unconventional T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Innate / Gamma-Delta T",
        trajectory_group="Unconventional",
        functional_state="CD4-CD8- double-negative SOX4+IKZF2+ unconventional T cell",
        canonical_markers=("SOX4", "IKZF2", "TRAC", "CD3D", "LAT"),
    ),
    "T_c10_NME1": TCellSubstateMetadata(
        cluster_id="T_c10_NME1",
        short_name="T_NME1",
        broad_cell_type="CD8+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Metabolic / Cycling NME1+ T",
        trajectory_group="Cycling",
        functional_state="Nucleotide-biosynthetic NME1+NME2+FABP5+HSP90AB1+ activated T cell",
        canonical_markers=("NME1", "NME2", "FABP5", "MYC", "C1QBP"),
    ),
    "T_c11_Th2_GATA3": TCellSubstateMetadata(
        cluster_id="T_c11_Th2_GATA3",
        short_name="CD4_Th2",
        broad_cell_type="CD4+ T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Naive / Memory CD4",
        trajectory_group="CD4_Helper",
        functional_state="GATA3+CCR4+IL4I1+PTGDR2+ Th2 helper CD4+ T cell",
        canonical_markers=("GATA3", "CCR4", "IL13", "IL4R", "AREG"),
    ),
    "T_c12_Treg_ISG": TCellSubstateMetadata(
        cluster_id="T_c12_Treg_ISG",
        short_name="CD4_Treg_ISG",
        broad_cell_type="Regulatory T Cell",
        lineage_compartment="T_CELL",
        lineage_branch="Suppressive Treg (FOXP3+/TNFRSF9+)",
        trajectory_group="CD4_Treg",
        functional_state="Interferon-stimulated FOXP3+ISG15+IFIT1+OAS1+ regulatory T cell",
        canonical_markers=("FOXP3", "IL2RA", "ISG15", "IFIT1", "IFIT3", "OAS1"),
    ),
}


# --- 1C. Shi et al. (2025, Nature — Cross-Tissue Multicellular Coordination) ---
SHI_PANCANCER_DATASET_ID: Final[str] = "Shi_PanCancer_Atlas"
SHI_PANCANCER_ZENODO_RECORD: Final[str] = "15554080"
SHI_PANCANCER_ZENODO_DOI: Final[str] = "10.5281/zenodo.15554080"
SHI_PANCANCER_PUBLICATION_DOI: Final[str] = "10.1038/s41586-025-09053-4"
SHI_PANCANCER_PUBLICATION_PMID: Final[str] = "40437094"

SHI_PANCANCER_H5AD_FILENAME: Final[str] = "PanCancer_igt_s9_fine_counts.h5ad"
SHI_PANCANCER_ZENODO_URL: Final[str] = (
    "https://zenodo.org/records/15554080/files/PanCancer_igt_s9_fine_counts.h5ad?download=1"
)
SHI_PANCANCER_CHECKSUM: Final[ChecksumSpec] = ChecksumSpec(
    algorithm="md5",
    expected_hash="95547a82eebf08e1c1d93c837693f23d",
    size_bytes=Some(18_204_274_013),
)

_ZHENG_CANCER_CODE_TO_FULL: Final[Mapping[str, tuple[str, str]]] = {
    "BC": ("Breast Cancer", "Breast"),
    "BRCA": ("Breast Cancer", "Breast"),
    "BCL": ("B-Cell Lymphoma", "Lymph Node"),
    "CHOL": ("Cholangiocarcinoma", "Liver / Bile Duct"),
    "CRC": ("Colorectal Cancer", "Colon / Rectum"),
    "ESCA": ("Esophageal Carcinoma", "Esophagus"),
    "FTC": ("Follicular Thyroid Carcinoma", "Thyroid"),
    "GIST": ("Gastrointestinal Stromal Tumor", "Stomach"),
    "HCC": ("Hepatocellular Carcinoma", "Liver"),
    "HNSCC": ("Head and Neck Squamous Cell Carcinoma", "Head and Neck"),
    "LC": ("Lung Cancer", "Lung"),
    "NSCLC": ("Non-Small Cell Lung", "Lung"),
    "LUAD": ("Lung Adenocarcinoma", "Lung"),
    "MELA": ("Melanoma", "Skin"),
    "MEL": ("Melanoma", "Skin"),
    "MM": ("Multiple Myeloma", "Bone Marrow"),
    "NPC": ("Nasopharyngeal Carcinoma", "Nasopharynx"),
    "OV": ("Ovarian Cancer", "Ovary"),
    "PACA": ("Pancreatic Adenocarcinoma", "Pancreas"),
    "PDAC": ("Pancreatic Ductal Adenocarcinoma", "Pancreas"),
    "RC": ("Renal Cell Carcinoma", "Kidney"),
    "RCC": ("Renal Cell Carcinoma", "Kidney"),
    "STAD": ("Gastric Adenocarcinoma", "Stomach"),
    "THCA": ("Thyroid Carcinoma", "Thyroid"),
    "UCEC": ("Uterine Corpus Endometrial Carcinoma", "Uterus"),
    "UM": ("Uveal Melanoma", "Eye"),
    "Normal": ("Healthy / Normal Donor", "Blood / Normal"),
}


# =============================================================================
# 2. Pure Functional Core: Polars Metadata Harmonizers & Sparse CSR Operators
# =============================================================================


def _pandas_to_polars_preserving_all_columns(
    raw_obs: pd.DataFrame,
    cell_names_fallback: Sequence[str] = (),
) -> pl.DataFrame:
    """Convert any pandas `.obs` DataFrame into a typed Polars DataFrame preserving 100% of columns."""
    n_rows = len(raw_obs) if len(raw_obs) > 0 else len(cell_names_fallback)
    if n_rows == 0:
        return pl.DataFrame({"cell_id": pl.Series("cell_id", [], dtype=pl.String)})

    barcodes = (
        [str(idx) for idx in raw_obs.index]
        if len(raw_obs) == n_rows and len(raw_obs.index) == n_rows
        else [str(c) for c in cell_names_fallback[:n_rows]]
    )

    series_list: list[pl.Series] = [pl.Series("cell_id", barcodes, dtype=pl.String)]
    for col in raw_obs.columns:
        col_str = str(col).strip()
        if not col_str or col_str == "cell_id":
            continue
        s = raw_obs[col]
        if pd.api.types.is_bool_dtype(s):
            series_list.append(pl.Series(col_str, s.fillna(False).to_numpy(dtype=bool), dtype=pl.Boolean))
        elif pd.api.types.is_integer_dtype(s):
            series_list.append(pl.Series(col_str, s.to_numpy(dtype=np.int64), dtype=pl.Int64))
        elif pd.api.types.is_float_dtype(s):
            series_list.append(pl.Series(col_str, s.to_numpy(dtype=np.float64), dtype=pl.Float64))
        else:
            cleaned = [
                "NA" if pd.isna(v) or str(v).strip() in ("", "nan", "None", "<NA>") else str(v).strip()
                for v in s.tolist()
            ]
            series_list.append(pl.Series(col_str, cleaned, dtype=pl.String))

    return pl.DataFrame(series_list)


def _finalize_polars_to_obs_pandas(df_pl: pl.DataFrame) -> pd.DataFrame:
    """Convert harmonized Polars `.obs` back to an AnnData-ready pandas DataFrame indexed by `cell_id`."""
    if df_pl.is_empty():
        return pd.DataFrame(index=pd.Index([], name="cell_id"))
    pdf = df_pl.to_pandas()
    unique_ids = ad.utils.make_index_unique(pd.Index([str(x) for x in pdf["cell_id"].tolist()], name="cell_id"))
    pdf["cell_id"] = [str(x) for x in unique_ids.tolist()]
    pdf.index = unique_ids
    return pdf


def _pick_first_existing_str_expr(
    schema_cols: Sequence[str],
    candidates: Sequence[str],
    fallback: str,
) -> pl.Expr:
    """Return a Polars string expression picking the first non-NA value among candidate columns."""
    valid_exprs = [
        pl.when(
            pl.col(c).cast(pl.String).is_not_null()
            & ~pl.col(c).cast(pl.String).str.strip_chars().is_in(["", "NA", "nan", "None", "unknown"])
        )
        .then(pl.col(c).cast(pl.String).str.strip_chars())
        .otherwise(pl.lit(""))
        for c in candidates
        if c in schema_cols
    ]
    if not valid_exprs:
        return pl.lit(fallback)
    chained = pl.lit(fallback)
    for expr in reversed(valid_exprs):
        chained = pl.when(expr != "").then(expr).otherwise(chained)
    return chained


def _normalize_zheng_sample_loc(loc_tok: str) -> str:
    """Map Zheng et al. `loc` codes (`T`, `N`, `P`, `L`) to standardized `sample_type`."""
    match loc_tok.strip().upper():
        case "T" | "TUMOR" | "TUMOUR":
            return "tumor"
        case "N" | "NORMAL" | "ADJACENT":
            return "adjacent_normal"
        case "P" | "PBMC" | "BLOOD":
            return "blood"
        case "L" | "LN" | "LYMPH_NODE":
            return "lymph_node"
        case _:
            return "unknown"


def _resolve_zheng_cluster_key(raw_tok: str) -> str:
    """Resolve a raw Zheng cluster token (`CD8.c12.Tex.CXCL13`, `CD8.c12`, `CD8_c12`, etc.) to canonical key."""
    cleaned = raw_tok.strip()
    if cleaned in ZHENG_TCELL_CLUSTER_SPECS:
        return cleaned
    norm_dots = cleaned.replace("_", ".")
    if norm_dots in ZHENG_TCELL_CLUSTER_SPECS:
        return norm_dots
    prefix = ".".join(norm_dots.split(".")[:2])
    for full_key in ZHENG_TCELL_CLUSTER_SPECS:
        if full_key.startswith(prefix + "."):
            return full_key
    return cleaned


def harmonize_zheng_tcell_obs_polars(
    raw_obs: pd.DataFrame,
    cell_names_fallback: Sequence[str] = (),
    clinical_df: Maybe[pl.DataFrame] = Nothing,
    tcr_df: Maybe[pl.DataFrame] = Nothing,
    pseudotime_df: Maybe[pl.DataFrame] = Nothing,
) -> pd.DataFrame:
    """Harmonize Zheng et al. (2021, Science) T-cell metadata via pure Polars transforms.

    Preserves 100% of original author metadata columns (`meta.cluster`, `meta.cluster.coarse`,
    `miniCluster`, `ClusterID.harmony`, `DIG.Score1`, `S.Score`, `G2M.Score`, `Phase`, etc.),
    enriches donor clinical Stage/Age/Gender (`TableS1.clinicalInfo.P316.txt`), TCR clonality &
    pTRT tumor-reactivity scores (`data.tcr.tar.gz`), and trajectory pseudotime when provided.
    """
    df = _pandas_to_polars_preserving_all_columns(raw_obs, cell_names_fallback=cell_names_fallback)
    if df.is_empty():
        return _finalize_polars_to_obs_pandas(df)

    cols = list(df.columns)
    raw_cluster_expr = _pick_first_existing_str_expr(
        cols,
        ("meta.cluster", "majorCluster", "meta_cluster", "cluster.name", "ClusterID"),
        "Unassigned",
    )
    cancer_code_expr = _pick_first_existing_str_expr(
        cols,
        ("cancerType", "cancer_type", "cancer"),
        "Pan-Cancer",
    )
    loc_expr = _pick_first_existing_str_expr(cols, ("loc", "tissue", "sample_type"), "T")
    patient_expr = _pick_first_existing_str_expr(cols, ("patient", "patient_id", "donor"), "NA")
    sample_expr = _pick_first_existing_str_expr(cols, ("libraryID", "batchV", "sample", "patient"), "NA")
    study_expr = _pick_first_existing_str_expr(cols, ("dataset", "dataset.tech", "study"), ZHENG_TCELL_DATASET_ID)
    raw_cell_expr = _pick_first_existing_str_expr(cols, ("cellID", "Cell_Name", "cell_id"), "NA")

    df = df.with_columns(
        pl.lit(ZHENG_TCELL_DATASET_ID).alias("dataset_id"),
        study_expr.alias("source_study"),
        patient_expr.alias("patient_id"),
        sample_expr.alias("sample_id"),
        cancer_code_expr.alias("cancer_code"),
        loc_expr.map_elements(_normalize_zheng_sample_loc, return_dtype=pl.String).alias("sample_type"),
        raw_cluster_expr.map_elements(_resolve_zheng_cluster_key, return_dtype=pl.String).alias("author_cluster_id"),
        raw_cell_expr.alias("_join_cell_name"),
    )

    cluster_ids = df["author_cluster_id"].to_list()
    cancer_codes = df["cancer_code"].to_list()

    cell_types: list[str] = []
    subtypes: list[str] = []
    states: list[str] = []
    compartments: list[str] = []
    branches: list[str] = []
    traj_groups: list[str] = []
    func_states: list[str] = []
    for cid in cluster_ids:
        if cid in ZHENG_TCELL_CLUSTER_SPECS:
            spec = ZHENG_TCELL_CLUSTER_SPECS[cid]
            cell_types.append(spec.broad_cell_type)
            subtypes.append(f"{cid.split('.')[0]}_{spec.short_name}")
            states.append(cid)
            compartments.append(spec.lineage_compartment)
            branches.append(spec.lineage_branch)
            traj_groups.append(spec.trajectory_group)
            func_states.append(spec.functional_state)
        else:
            inferred_bt = "CD4+ T Cell" if cid.startswith("CD4") else "CD8+ T Cell"
            cell_types.append(inferred_bt)
            subtypes.append(cid)
            states.append(cid)
            compartments.append("T_CELL")
            branches.append("Unclassified T Cell")
            traj_groups.append("NA")
            func_states.append(cid)

    cancer_names: list[str] = []
    organs: list[str] = []
    for cc in cancer_codes:
        full_nm, org_nm = _ZHENG_CANCER_CODE_TO_FULL.get(cc, (cc, "Multi-organ"))
        cancer_names.append(full_nm)
        organs.append(org_nm)

    df = df.with_columns(
        pl.Series("cancer_type", cancer_names, dtype=pl.String),
        pl.Series("organ", organs, dtype=pl.String),
        pl.Series("cell_type", cell_types, dtype=pl.String),
        pl.Series("cell_subtype", subtypes, dtype=pl.String),
        pl.Series("cell_state", states, dtype=pl.String),
        pl.Series("lineage_compartment", compartments, dtype=pl.String),
        pl.Series("lineage_branch", branches, dtype=pl.String),
        pl.Series("trajectory_group", traj_groups, dtype=pl.String),
        pl.Series("functional_state", func_states, dtype=pl.String),
        (pl.col("source_study") + pl.lit(".") + pl.col("patient_id")).alias("_donor_join_key"),
    )

    match clinical_df:
        case Some(clin) if not clin.is_empty() and "_donor_join_key" in clin.columns:
            clin_dedup = clin.unique(subset=["_donor_join_key"], keep="first")
            df = df.join(clin_dedup, on="_donor_join_key", how="left")
        case _:
            pass

    match tcr_df:
        case Some(tcr) if not tcr.is_empty() and "_join_cell_name" in tcr.columns:
            join_keys = (
                ["source_study", "_join_cell_name"]
                if "source_study" in tcr.columns
                else ["_join_cell_name"]
            )
            tcr_dedup = tcr.unique(subset=join_keys, keep="first")
            df = df.join(tcr_dedup, on=join_keys, how="left")
        case _:
            pass

    match pseudotime_df:
        case Some(pt) if not pt.is_empty() and "miniCluster" in pt.columns and "miniCluster" in df.columns:
            pt_dedup = pt.unique(subset=["miniCluster"], keep="first")
            df = df.join(pt_dedup, on="miniCluster", how="left")
        case _:
            pass

    drop_helpers = [c for c in ("_join_cell_name", "_donor_join_key") if c in df.columns]
    if drop_helpers:
        df = df.drop(drop_helpers)

    fill_exprs: list[pl.Expr] = [
        pl.col(c).fill_null("NA").alias(c)
        for c, dt in df.schema.items()
        if dt == pl.String
    ]
    if fill_exprs:
        df = df.with_columns(fill_exprs)

    return _finalize_polars_to_obs_pandas(df)


def _resolve_chu_cluster_key(raw_tok: str) -> str:
    """Resolve raw Chu et al. (2023) cluster labels (`CD8_c4_Tstr`, `c4_Tstr`, `Tstr`, etc.) to canonical key."""
    cleaned = raw_tok.strip()
    if cleaned in CHU_TCELL_CLUSTER_SPECS:
        return cleaned
    norm = cleaned.replace(".", "_").replace(" ", "_")
    if norm in CHU_TCELL_CLUSTER_SPECS:
        return norm
    low = norm.lower()
    for key, spec in CHU_TCELL_CLUSTER_SPECS.items():
        if key.lower() == low or spec.short_name.lower() == low:
            return key
    if "tstr" in low or "hsp" in low or "stress" in low:
        return "CD4_c4_Tstr" if "cd4" in low else "CD8_c4_Tstr"
    return cleaned


def harmonize_chu_tcell_obs_polars(
    raw_obs: pd.DataFrame,
    active_ident: Sequence[str] = (),
    cell_names_fallback: Sequence[str] = (),
) -> pd.DataFrame:
    """Harmonize Chu et al. (2023, Nature Medicine) Pan-Cancer T-cell metadata via pure Polars transforms."""
    df = _pandas_to_polars_preserving_all_columns(raw_obs, cell_names_fallback=cell_names_fallback)
    if df.is_empty():
        return _finalize_polars_to_obs_pandas(df)

    if len(active_ident) == len(df) and "active_ident" not in df.columns:
        df = df.with_columns(
            pl.Series("active_ident", [str(x) for x in active_ident], dtype=pl.String)
        )

    cols = list(df.columns)
    cluster_expr = _pick_first_existing_str_expr(
        cols,
        ("cell.type", "cell_type", "Cluster", "cluster", "subcluster", "meta.cluster", "active_ident", "seurat_clusters"),
        "Unassigned",
    )
    cancer_expr = _pick_first_existing_str_expr(
        cols,
        ("CancerType", "cancer_type", "cancer", "TumorType", "cohort"),
        "Pan-Cancer",
    )
    organ_expr = _pick_first_existing_str_expr(
        cols,
        ("Tissue", "tissue", "organ", "Organ"),
        "Multi-organ",
    )
    sample_type_expr = _pick_first_existing_str_expr(
        cols,
        ("TissueType", "sample_type", "Source", "loc"),
        "tumor",
    )
    patient_expr = _pick_first_existing_str_expr(
        cols,
        ("Patient", "patient", "patient_id", "Donor", "donor_id"),
        "NA",
    )
    sample_expr = _pick_first_existing_str_expr(
        cols,
        ("Sample", "sample", "sample_id", "orig.ident"),
        "NA",
    )
    study_expr = _pick_first_existing_str_expr(
        cols,
        ("Dataset", "dataset", "Study", "study"),
        CHU_TCELL_DATASET_ID,
    )

    df = df.with_columns(
        pl.lit(CHU_TCELL_DATASET_ID).alias("dataset_id"),
        study_expr.alias("source_study"),
        patient_expr.alias("patient_id"),
        sample_expr.alias("sample_id"),
        cancer_expr.alias("cancer_type"),
        organ_expr.alias("organ"),
        sample_type_expr.alias("sample_type"),
        cluster_expr.map_elements(_resolve_chu_cluster_key, return_dtype=pl.String).alias("author_cluster_id"),
    )

    cids = df["author_cluster_id"].to_list()
    cell_types: list[str] = []
    subtypes: list[str] = []
    branches: list[str] = []
    traj_groups: list[str] = []
    func_states: list[str] = []
    is_tstr: list[bool] = []

    for cid in cids:
        if cid in CHU_TCELL_CLUSTER_SPECS:
            spec = CHU_TCELL_CLUSTER_SPECS[cid]
            cell_types.append(spec.broad_cell_type)
            subtypes.append(spec.short_name)
            branches.append(spec.lineage_branch)
            traj_groups.append(spec.trajectory_group)
            func_states.append(spec.functional_state)
            is_tstr.append(cid in ("CD8_c4_Tstr", "CD4_c4_Tstr"))
        else:
            bt = (
                "Regulatory T Cell"
                if "treg" in cid.lower()
                else ("CD4+ T Cell" if "cd4" in cid.lower() else "CD8+ T Cell")
            )
            cell_types.append(bt)
            subtypes.append(cid)
            branches.append("Unclassified T Cell")
            traj_groups.append("NA")
            func_states.append(cid)
            is_tstr.append("tstr" in cid.lower())

    df = df.with_columns(
        pl.Series("cell_type", cell_types, dtype=pl.String),
        pl.Series("cell_subtype", subtypes, dtype=pl.String),
        pl.col("author_cluster_id").alias("cell_state"),
        pl.lit("T_CELL").alias("lineage_compartment"),
        pl.Series("lineage_branch", branches, dtype=pl.String),
        pl.Series("trajectory_group", traj_groups, dtype=pl.String),
        pl.Series("functional_state", func_states, dtype=pl.String),
        pl.Series("is_stress_tstr", is_tstr, dtype=pl.Boolean),
    )

    return _finalize_polars_to_obs_pandas(df)


def _classify_shi_substate_ontology(
    fine_label: str,
    major_label: str,
) -> tuple[str, str, str, str]:
    """Classify a Shi et al. (2025, Nature) cell annotation into `(compartment, branch, broad_type, functional_state)`."""
    combined = f"{major_label} {fine_label}".lower()

    # 0. Exclude specialized non-CAF parenchymal/neural-crest clusters grouped under `Stromal` in Shi et al.
    if any(
        non_caf_tok in combined
        for non_caf_tok in ("melanocyte", "glia_", "cmc_", "skmc_", "satellite_")
    ):
        return (
            "OTHER",
            "Specialized Non-CAF Parenchymal",
            fine_label or "Other",
            f"Shi et al. non-CAF parenchymal/neural-crest state ({fine_label})",
        )

    # 1. Stromal / CAF / Fibroblast / Mural compartment (check before endothelial/other)
    if any(
        tok in combined
        for tok in (
            "fibro",
            "caf",
            "myofib",
            "icaf",
            "mcaf",
            "mycaf",
            "apcaf",
            "pericyte",
            "smooth muscle",
            "vsmc",
            "smc",
            "mural",
            "mesenchym",
            "stromal",
            "stellate",
        )
    ) and "endothel" not in combined:
        if any(t in combined for t in ("pericyte", "rgs5", "mural")):
            return ("STROMAL_CAF", "Mural / Pericyte / SMC", "Pericyte", f"Shi et al. mural/pericyte state ({fine_label})")
        if any(t in combined for t in ("smooth muscle", "vsmc", "smc", "myh11", "actg2", "rergl", "fbln5")):
            return ("STROMAL_CAF", "Mural / Pericyte / SMC", "Smooth Muscle Cell", f"Shi et al. vascular/visceral SMC state ({fine_label})")
        if any(t in combined for t in ("mycaf", "mcaf", "myofib", "lrrc15", "fap", "acta2", "postn", "mmp11", "cthrc1")):
            return ("STROMAL_CAF", "Myofibroblastic CAF (myCAF)", "Myofibroblast", f"Shi et al. myofibroblastic CAF state ({fine_label})")
        if any(t in combined for t in ("icaf", "inflam", "il6", "cxcl", "ccl", "cfd", "complement", "c3", "pla2g2a")):
            return ("STROMAL_CAF", "Inflammatory CAF (iCAF)", "Fibroblast", f"Shi et al. inflammatory CAF state ({fine_label})")
        if any(t in combined for t in ("apcaf", "cd74", "hla", "antigen", "mesothel")):
            return ("STROMAL_CAF", "Antigen-Presenting / Mesothelial CAF", "Fibroblast", f"Shi et al. antigen-presenting CAF state ({fine_label})")
        if any(t in combined for t in ("pi16", "cd34", "adventit", "progenitor", "universal")):
            return ("STROMAL_CAF", "Universal Progenitor Fibroblast", "Fibroblast", f"Shi et al. universal progenitor fibroblast ({fine_label})")
        return ("STROMAL_CAF", "Homeostatic Matrix Fibroblast", "Fibroblast", f"Shi et al. stromal fibroblast state ({fine_label})")

    # 2. NK / Innate Lymphoid Cell compartment (exclude gamma-delta T cells such as `I09_gdT_ITGA1`)
    if any(tok in combined for tok in ("nk cell", "natural killer", "ilc", "nk_", "cd16hink", "cd16lonk", "cd56", "klrf1")) and not any(
        t_tok in combined for t_tok in ("nkt", "nk-like", "cd8", "cd4", "t cell", "t_cell", "gdt", "gd_t")
    ):
        return ("NK_ILC", "Cytotoxic NK / ILC", "NK Cell", f"Shi et al. innate lymphoid / NK state ({fine_label})")

    # 3. T-Cell compartment
    if any(
        tok in combined
        for tok in (
            "t cell",
            "t_cell",
            "t-cell",
            "cd8",
            "cd4",
            "treg",
            "tfh",
            "th1",
            "th17",
            "th2",
            "mait",
            "gd_t",
            "gdt",
            "gamma-delta",
            "gamma delta",
            "nkt",
            "tex",
            "temra",
            "trm",
            "lymphocyte",
        )
    ) and not any(b_tok in combined for b_tok in ("b cell", "b_cell", "b-cell", "plasma", "plasmablast")):
        if any(t in combined for t in ("treg", "foxp3", "regulatory")):
            return ("T_CELL", "Suppressive Treg (FOXP3+/TNFRSF9+)", "Regulatory T Cell", f"Shi et al. regulatory T-cell state ({fine_label})")
        if any(t in combined for t in ("tex", "exhaust", "cxcl13", "havcr2", "pdcd1")) and "cd4" not in combined:
            if any(p in combined for p in ("tcf7", "prog", "pre", "stem", "gzmk")):
                return ("T_CELL", "Progenitor Exhausted CD8 (TCF7+)", "CD8+ T Cell", f"Shi et al. progenitor exhausted CD8 state ({fine_label})")
            return ("T_CELL", "Terminal Exhausted CD8 (CXCL13+)", "CD8+ T Cell", f"Shi et al. exhausted CD8 state ({fine_label})")
        if any(t in combined for t in ("tfh", "cxcr5", "th1", "bcl6")) or ("cxcl13" in combined and "cd4" in combined):
            return ("T_CELL", "Follicular / Th1 Helper CD4 (CXCL13+/CXCR5+)", "CD4+ T Cell", f"Shi et al. Tfh/Th1 helper CD4 state ({fine_label})")
        if any(t in combined for t in ("th17", "il17", "rorc")) and "cd8" not in combined:
            return ("T_CELL", "Th17 Helper CD4", "CD4+ T Cell", f"Shi et al. Th17 CD4 state ({fine_label})")
        if any(t in combined for t in ("isg", "ifit", "interferon", "stat1", "mx1")):
            bt_isg = "CD4+ T Cell" if "cd4" in combined else "CD8+ T Cell"
            return ("T_CELL", "Interferon-Stimulated T (ISG+)", bt_isg, f"Shi et al. interferon-stimulated T state ({fine_label})")
        if any(t in combined for t in ("hsp", "stress", "tstr", "dnajb1")):
            bt_str = "CD4+ T Cell" if "cd4" in combined else "CD8+ T Cell"
            return ("T_CELL", "Stress-Response T (Tstr)", bt_str, f"Shi et al. stress-response T state ({fine_label})")
        if any(t in combined for t in ("prolif", "cycling", "mki67", "top2a", "stmn1")):
            return ("T_CELL", "Metabolic / Cycling NME1+ T", "Proliferating T Cell", f"Shi et al. proliferating T state ({fine_label})")
        if any(t in combined for t in ("gd", "gamma", "delta", "trdc", "trdv")) and "cd8" not in combined and "cd4" not in combined:
            return ("T_CELL", "Innate / Gamma-Delta T", "Gamma-Delta T Cell", f"Shi et al. gamma-delta T state ({fine_label})")
        if any(t in combined for t in ("trm", "znf683", "itgae", "resident")):
            return ("T_CELL", "Tissue-Resident Memory CD8 (ZNF683+)", "CD8+ T Cell", f"Shi et al. tissue-resident memory T state ({fine_label})")
        if any(t in combined for t in ("temra", "cx3cr1", "fgfbp2", "eff", "ctl", "cytotox", "nkt")):
            bt_eff = "CD4+ T Cell" if "cd4" in combined else "CD8+ T Cell"
            return ("T_CELL", "Cytotoxic Temra / NK-like CD8", bt_eff, f"Shi et al. cytotoxic effector T state ({fine_label})")
        if any(t in combined for t in ("gzmk", "tem")) and "cd4" not in combined:
            return ("T_CELL", "Effector Memory CD8 (GZMK+)", "CD8+ T Cell", f"Shi et al. GZMK+ effector-memory CD8 state ({fine_label})")
        if any(t in combined for t in ("gzmk", "tem")) and "cd4" in combined:
            return ("T_CELL", "Cytotoxic / Effector Memory CD4", "CD4+ T Cell", f"Shi et al. GZMK+ effector-memory CD4 state ({fine_label})")
        if any(t in combined for t in ("mait", "slc4a10", "tc17")):
            return ("T_CELL", "MAIT / Mucosal Tc17", "MAIT / Tc17 Cell", f"Shi et al. MAIT/Tc17 state ({fine_label})")
        if "cd4" in combined:
            return ("T_CELL", "Naive / Memory CD4", "CD4+ T Cell", f"Shi et al. CD4+ T state ({fine_label})")
        return ("T_CELL", "Naive / Central Memory CD8", "CD8+ T Cell", f"Shi et al. CD8+/T state ({fine_label})")

    if any(tok in combined for tok in ("macro", "mono", "dc", "dendritic", "myeloid", "neutrophil", "mast", "microglia")):
        return ("MYELOID", "Myeloid / Mononuclear Phagocyte", "Myeloid Cell", f"Shi et al. myeloid state ({fine_label})")
    if any(tok in combined for tok in ("b cell", "b_cell", "b-cell", "plasma", "plasmablast")):
        return ("B_PLASMA", "B / Plasma Lineage", "B / Plasma Cell", f"Shi et al. B/plasma state ({fine_label})")
    if any(tok in combined for tok in ("endothel", "vascular", "lymphatic", "lec", "vec")):
        return ("ENDOTHELIAL", "Endothelial Lineage", "Endothelial Cell", f"Shi et al. endothelial state ({fine_label})")
    if any(tok in combined for tok in ("epithel", "malignant", "tumor", "cancer", "carcinoma", "melanocyte", "hepatocyte")):
        return ("EPITHELIAL_MALIGNANT", "Epithelial / Malignant", "Epithelial / Malignant Cell", f"Shi et al. epithelial/malignant state ({fine_label})")

    return ("OTHER", "Other Lineage", fine_label or "Other", f"Shi et al. cell state ({fine_label})")


def harmonize_shi_pancancer_obs_polars(
    raw_obs: pd.DataFrame,
    cell_names_fallback: Sequence[str] = (),
) -> pd.DataFrame:
    """Harmonize Shi et al. (2025, Nature) `PanCancer_igt_s9_fine_counts.h5ad` `.obs` via Polars."""
    df = _pandas_to_polars_preserving_all_columns(raw_obs, cell_names_fallback=cell_names_fallback)
    if df.is_empty():
        return _finalize_polars_to_obs_pandas(df)

    cols = list(df.columns)
    fine_expr = _pick_first_existing_str_expr(
        cols,
        (
            "subCluster",
            "Prediction",
            "fine_cell_type",
            "cell_type_fine",
            "sub_cell_type",
            "subtype",
            "cell_subtype",
            "minor_celltype",
            "celltype_fine",
            "s9_fine",
            "igt_s9_fine",
            "cluster",
            "annotation",
            "cellType2",
            "cell_type",
            "celltype",
        ),
        "Unassigned",
    )
    major_expr = _pick_first_existing_str_expr(
        cols,
        (
            "majorCluster",
            "ann1",
            "major_cell_type",
            "cell_type_major",
            "major_celltype",
            "compartment",
            "lineage",
            "broad_cell_type",
            "cellType1",
            "cell_type",
            "celltype",
        ),
        "NA",
    )
    cancer_expr = _pick_first_existing_str_expr(
        cols,
        (
            "cancerType",
            "cancer_type",
            "CancerType",
            "sub_cancerType",
            "tumor_type",
            "disease",
            "cancer",
            "cohort_cancer",
        ),
        "Pan-Cancer",
    )
    organ_expr = _pick_first_existing_str_expr(
        cols,
        ("cancerType", "tissue", "Tissue", "organ", "Organ", "site", "tissue_general"),
        "Multi-organ",
    )
    raw_site_expr = _pick_first_existing_str_expr(
        cols,
        ("sample_type", "tissue", "tissue_type", "origin", "condition", "tumor_stage"),
        "tumor",
    )
    sample_type_expr = (
        pl.when(raw_site_expr.str.to_lowercase().is_in(["tumor", "tumour", "primary"]))
        .then(pl.lit("tumor"))
        .when(raw_site_expr.str.to_lowercase() == "metastasis")
        .then(pl.lit("metastasis"))
        .when(raw_site_expr.str.to_lowercase().is_in(["adjacent", "normal_adjacent", "normal"]))
        .then(pl.lit("normal_adjacent"))
        .when(raw_site_expr.str.to_lowercase().is_in(["blood", "pbmc"]))
        .then(pl.lit("blood"))
        .when(raw_site_expr.str.to_lowercase().is_in(["prelesion", "pre_lesion"]))
        .then(pl.lit("pre_lesion"))
        .when(raw_site_expr.str.to_lowercase().is_in(["pleuralfluids", "ascites", "effusion"]))
        .then(pl.lit("effusion"))
        .otherwise(raw_site_expr.str.to_lowercase())
    )
    patient_expr = _pick_first_existing_str_expr(
        cols,
        ("donorID", "donor_id", "patient_id", "patient", "Patient", "subject_id"),
        "NA",
    )
    sample_expr = _pick_first_existing_str_expr(
        cols,
        ("sampleID", "sample_id", "sample", "Sample", "library_id", "orig.ident"),
        "NA",
    )
    study_expr = _pick_first_existing_str_expr(
        cols,
        ("cohortName", "datasetID", "cohortID", "dataset", "study", "Dataset", "Study", "batch", "cohort"),
        SHI_PANCANCER_DATASET_ID,
    )

    df = df.with_columns(
        pl.lit(SHI_PANCANCER_DATASET_ID).alias("dataset_id"),
        study_expr.alias("source_study"),
        patient_expr.alias("patient_id"),
        sample_expr.alias("sample_id"),
        cancer_expr.alias("cancer_type"),
        organ_expr.alias("organ"),
        sample_type_expr.alias("sample_type"),
        fine_expr.alias("author_cluster_id"),
        major_expr.alias("author_major_lineage"),
    )

    fine_list = df["author_cluster_id"].to_list()
    major_list = df["author_major_lineage"].to_list()
    unique_pairs = list(dict.fromkeys(zip(fine_list, major_list)))
    pair_lookup = {pair: _classify_shi_substate_ontology(pair[0], pair[1]) for pair in unique_pairs}

    compartments: list[str] = []
    branches: list[str] = []
    broad_types: list[str] = []
    func_states: list[str] = []
    for f_lbl, m_lbl in zip(fine_list, major_list):
        comp, br, bt, fs = pair_lookup[(f_lbl, m_lbl)]
        compartments.append(comp)
        branches.append(br)
        broad_types.append(bt)
        func_states.append(fs)

    df = df.with_columns(
        pl.Series("cell_type", broad_types, dtype=pl.String),
        pl.col("author_cluster_id").alias("cell_subtype"),
        pl.col("author_cluster_id").alias("cell_state"),
        pl.Series("lineage_compartment", compartments, dtype=pl.String),
        pl.Series("lineage_branch", branches, dtype=pl.String),
        pl.Series("functional_state", func_states, dtype=pl.String),
    )

    return _finalize_polars_to_obs_pandas(df)


def select_shi_compartment_mask(
    obs_df: pd.DataFrame,
    compartment_filter: ShiCompartmentFilter,
) -> np.ndarray[Any, Any]:
    """Exhaustive `match...case` selector returning a boolean row mask for `Shi_PanCancer_Atlas`."""
    comp_arr = obs_df["lineage_compartment"].astype(str).to_numpy()
    match compartment_filter:
        case ShiCompartmentFilter.T_AND_CAF:
            return np.isin(comp_arr, ["T_CELL", "NK_ILC", "STROMAL_CAF"])
        case ShiCompartmentFilter.T_NK_ONLY:
            return np.isin(comp_arr, ["T_CELL", "NK_ILC"])
        case ShiCompartmentFilter.STROMAL_CAF_ONLY:
            return np.asarray(comp_arr == "STROMAL_CAF", dtype=bool)
        case ShiCompartmentFilter.ALL_CELLS:
            return np.ones(len(obs_df), dtype=bool)
        case _ as unreachable:
            assert_never(unreachable)


def _build_canonical_ensembl_var_and_lookup(
    ensembl_release: int = 111,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Build deterministic canonical Ensembl GRCh38 `.var` table (`1..22, X, Y, MT`) and fast index lookup."""
    ensembl_dir = get_ensembl_dir()
    ensembl = ensure_ensembl_release_installed(
        release=ensembl_release,
        species="human",
        ensembl_dir=ensembl_dir,
    )
    canonical_contigs = {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M"}
    genes = sorted(
        (g for g in ensembl.genes() if str(g.contig).replace("chr", "") in canonical_contigs),
        key=lambda g: g.gene_id,
    )

    gene_ids = [g.gene_id for g in genes]
    gene_names = [g.gene_name if g.gene_name else g.gene_id for g in genes]
    contigs = [str(g.contig) for g in genes]
    starts = [int(g.start) for g in genes]
    ends = [int(g.end) for g in genes]
    strands = [str(g.strand) for g in genes]
    biotypes = [str(g.biotype) for g in genes]

    var_df = pd.DataFrame(
        {
            "gene_name": gene_names,
            "original_id": gene_ids,
            "contig": contigs,
            "start": starts,
            "end": ends,
            "strand": strands,
            "biotype": biotypes,
        },
        index=pd.Index(gene_ids, name="gene_id"),
    )

    id_to_col: dict[str, int] = {}
    for col_idx, (gid, gname) in enumerate(zip(gene_ids, gene_names)):
        id_to_col[gid] = col_idx
        if gname and gname not in id_to_col:
            id_to_col[gname] = col_idx

    cache_path = ensembl_dir / f"mapping_cache_human_rel{ensembl_release}.parquet"
    if cache_path.is_file():
        try:
            pl_cache = pl.read_parquet(cache_path)
            if {"query", "gene_id"}.issubset(set(pl_cache.columns)):
                for q_raw, gid_raw in zip(
                    pl_cache["query"].to_list(),
                    pl_cache["gene_id"].to_list(),
                    strict=False,
                ):
                    q_key = str(q_raw).strip()
                    gid_key = str(gid_raw).strip()
                    if not q_key or q_key in id_to_col:
                        continue
                    id_to_col[q_key] = id_to_col.get(gid_key, -1)
        except Exception:  # noqa: BLE001
            pass

    return var_df, id_to_col


def _build_cohort_projection_to_canonical(
    cohort_genes: Sequence[str],
    n_canonical_vars: int,
    id_to_col: Mapping[str, int],
    ensembl_release: int = 111,
) -> tuple[np.ndarray[Any, Any], sp.csr_matrix, dict[str, int]]:
    """Map any cohort's gene identifiers (`ENSG...`, `ENSG.ver`, or HGNC symbol) onto `var_df`."""
    resolved_cols = np.full(len(cohort_genes), -1, dtype=np.int32)
    unresolved_indices: list[int] = []
    unresolved_queries: list[str] = []
    next_lookup: dict[str, int] = dict(id_to_col)

    for idx, raw_g in enumerate(cohort_genes):
        clean_g = str(raw_g).strip()
        stripped = clean_g.split(".")[0] if clean_g.startswith("ENSG") else clean_g
        if stripped in next_lookup:
            resolved_cols[idx] = next_lookup[stripped]
        elif clean_g in next_lookup:
            resolved_cols[idx] = next_lookup[clean_g]
        elif clean_g:
            unresolved_indices.append(idx)
            unresolved_queries.append(clean_g)

    if unresolved_queries:
        unique_queries = list(dict.fromkeys(unresolved_queries))
        ensembl_dir = get_ensembl_dir()
        ensembl = ensure_ensembl_release_installed(
            release=ensembl_release,
            species="human",
            ensembl_dir=ensembl_dir,
        )
        cache_path = ensembl_dir / f"mapping_cache_human_rel{ensembl_release}.parquet"
        map_df = _load_or_update_mapping_cache(unique_queries, ensembl, cache_path)
        q_col = "query" if "query" in map_df.columns else "original_id"
        e_col = "gene_id" if "gene_id" in map_df.columns else "ensembl_id"
        query_to_gid: dict[str, str] = {}
        for q_val, eid_val in zip(map_df[q_col].tolist(), map_df[e_col].tolist(), strict=False):
            if pd.notna(eid_val) and str(eid_val).startswith("ENSG"):
                query_to_gid[str(q_val).strip()] = str(eid_val).strip()
        for q_u in unique_queries:
            mapped_gid = query_to_gid.get(q_u, "")
            next_lookup[q_u] = next_lookup.get(mapped_gid, -1) if mapped_gid else -1
        for raw_idx, q_str in zip(unresolved_indices, unresolved_queries, strict=False):
            resolved_cols[raw_idx] = next_lookup.get(q_str, -1)

    keep_mask = resolved_cols >= 0
    kept_targets = resolved_cols[keep_mask]
    n_kept = int(keep_mask.sum())
    proj_matrix = sp.csr_matrix(
        (
            np.ones(n_kept, dtype=np.float32),
            (np.arange(n_kept, dtype=np.int32), kept_targets),
        ),
        shape=(n_kept, n_canonical_vars),
        dtype=np.float32,
    )
    return keep_mask, proj_matrix, next_lookup


def _normalize_counts_csr_to_log1p_cp10k(
    counts_csr: sp.csr_matrix,
    target_sum: float = 10_000.0,
) -> sp.csr_matrix:
    """Pure CSR row-wise `ln(1 + counts / total_counts * 10_000)` transform (`float32`)."""
    if counts_csr.nnz == 0:
        return counts_csr.astype(np.float32)
    row_sums = np.asarray(counts_csr.sum(axis=1)).ravel().astype(np.float64)
    inv_scale = np.where(row_sums > 0.0, target_sum / np.maximum(row_sums, 1e-12), 0.0)
    row_nnz = np.diff(counts_csr.indptr)
    per_nz_scale = np.repeat(inv_scale, row_nnz)
    log_data = np.log1p(counts_csr.data.astype(np.float64) * per_nz_scale).astype(np.float32)
    return sp.csr_matrix(
        (log_data, counts_csr.indices.copy(), counts_csr.indptr.copy()),
        shape=counts_csr.shape,
        dtype=np.float32,
    )


def _slice_staged_csr_batch(
    assay_grp: h5py.Group,
    start_row: int,
    end_row: int,
    n_cols: int,
) -> sp.csr_matrix:
    """Slice rows `[start_row:end_row]` from a staged CSR HDF5 group with bounded RAM."""
    indptr_slice = np.asarray(assay_grp["indptr"][start_row : end_row + 1], dtype=np.int64)
    p_start = int(indptr_slice[0])
    p_end = int(indptr_slice[-1])
    data_slice = np.asarray(assay_grp["data"][p_start:p_end], dtype=np.float32)
    indices_slice = np.asarray(assay_grp["indices"][p_start:p_end], dtype=np.int32)
    shifted_indptr = indptr_slice - p_start
    return sp.csr_matrix(
        (data_slice, indices_slice, shifted_indptr),
        shape=(end_row - start_row, n_cols),
        dtype=np.float32,
    )


def _slice_h5_csr_selected_rows(
    csr_grp: h5py.Group,
    window_start: int,
    window_end: int,
    selected_global_rows: np.ndarray[Any, Any],
    n_cols: int,
) -> sp.csr_matrix:
    """Slice a contiguous window `[window_start:window_end]` from HDF5 CSR and gather `selected_global_rows`."""
    full_window = _slice_staged_csr_batch(csr_grp, window_start, window_end, n_cols)
    local_rows = selected_global_rows - window_start
    return full_window[local_rows, :].tocsr()


def _apply_maybe_subset_and_subsample(
    adata: ad.AnnData,
    subset: Maybe[Mapping[str, str]],
    subsample_n: Maybe[int],
) -> ad.AnnData:
    """Apply `Maybe`-typed subsetting and subsampling without `Optional`/`None` annotations."""
    has_sub = isinstance(subset, Some)
    has_n = isinstance(subsample_n, Some)
    if not has_sub and not has_n:
        return adata
    if has_sub and has_n:
        return _apply_subset_and_subsample(
            adata,
            subset=subset.unwrap(),
            subsample_n=subsample_n.unwrap(),
        )
    if has_sub:
        return _apply_subset_and_subsample(adata, subset=subset.unwrap())
    return _apply_subset_and_subsample(adata, subsample_n=subsample_n.unwrap())


def _read_selective_h5ad(h5ad_path: Path, include_all_layers: bool) -> ad.AnnData:
    """Load an AnnData `.h5ad` selectively so memory stays bounded unless `include_all_layers=True`."""
    if include_all_layers:
        return ad.read_h5ad(h5ad_path)

    with h5py.File(h5ad_path, "r") as f:
        obs = anndata.io.read_elem(f["obs"])
        var = anndata.io.read_elem(f["var"])
        x_mat = anndata.io.read_elem(f["X"])
        uns = anndata.io.read_elem(f["uns"]) if "uns" in f else {}
        obsm: dict[str, Any] = {}
        if "obsm" in f:
            for k in f["obsm"].keys():
                obsm[k] = anndata.io.read_elem(f["obsm"][k])

    return ad.AnnData(
        X=x_mat,
        obs=obs,
        var=var,
        uns=uns,
        obsm=obsm,
        layers={"log1p_norm": x_mat},
    )


# =============================================================================
# 3. Imperative Shell: Downloaders, Out-of-Core Preprocessors & Loaders
# =============================================================================


def _write_md5_sidecar(file_path: Path, expected_md5: str) -> int:
    """Write a deterministic `.md5` sidecar file and return bytes written."""
    sidecar = file_path.with_name(file_path.name + ".md5")
    text = f"{expected_md5}  {file_path.name}\n"
    sidecar.write_text(text, encoding="utf-8")
    return len(text)


# -----------------------------------------------------------------------------
# 3A. Zheng et al. (2021, Science) — `Zheng_TCell_Atlas`
# -----------------------------------------------------------------------------


def _extract_zheng_metainfo_artifacts(
    meta_tar_path: Path,
    scratch_dir: Path,
) -> tuple[
    Maybe[pl.DataFrame],
    Maybe[pl.DataFrame],
    dict[str, np.ndarray[Any, Any]],
    dict[str, list[str]],
]:
    """Extract donor clinical table, pseudotime, metacell embeddings, and signature genes from `data.metaInfo.tar.gz`."""
    clinical_maybe: Maybe[pl.DataFrame] = Nothing
    pseudotime_maybe: Maybe[pl.DataFrame] = Nothing
    metacell_umap: dict[str, np.ndarray[Any, Any]] = {}
    sig_genes: dict[str, list[str]] = {}

    if not meta_tar_path.is_file():
        return (clinical_maybe, pseudotime_maybe, metacell_umap, sig_genes)

    scratch_dir.mkdir(parents=True, exist_ok=True)
    pt_records: dict[str, dict[str, float]] = {}

    with tarfile.open(meta_tar_path, "r:gz") as tf:
        for member in tf.getmembers():
            if not member.isfile():
                continue
            bname = Path(member.name).name

            if bname == "TableS1.clinicalInfo.P316.txt":
                f_obj = tf.extractfile(member)
                if f_obj:
                    raw_txt = f_obj.read().decode("utf-8", errors="replace")
                    pdf = pd.read_csv(io.StringIO(raw_txt), sep="\t", dtype=str).fillna("NA")
                    if "Patient ID" in pdf.columns:
                        clin_pl = pl.DataFrame(
                            {
                                "_donor_join_key": [str(x).strip() for x in pdf["Patient ID"].tolist()],
                                "donor_age": [
                                    str(x).strip() for x in pdf.get("Age", pd.Series(["NA"] * len(pdf))).tolist()
                                ],
                                "donor_gender": [
                                    str(x).strip() for x in pdf.get("Gender", pd.Series(["NA"] * len(pdf))).tolist()
                                ],
                                "donor_ptnm": [
                                    str(x).strip() for x in pdf.get("pTNM", pd.Series(["NA"] * len(pdf))).tolist()
                                ],
                                "donor_clinical_stage": [
                                    str(x).strip()
                                    for x in pdf.get("Clinical Stage", pd.Series(["NA"] * len(pdf))).tolist()
                                ],
                            }
                        )
                        clinical_maybe = Some(clin_pl)

            elif bname == "CD4.umap_harmony.txt.gz":
                f_obj = tf.extractfile(member)
                if f_obj:
                    dec = gzip.decompress(f_obj.read()).decode("utf-8", errors="replace")
                    for ln in dec.splitlines()[1:]:
                        parts = ln.split("\t")
                        if len(parts) >= 3 and parts[0].strip():
                            try:
                                metacell_umap[parts[0].strip()] = np.array(
                                    [float(parts[1]), float(parts[2])],
                                    dtype=np.float32,
                                )
                            except ValueError:
                                pass

            elif bname in (
                "CD8_Tex.dpt_pseudotime.txt.gz",
                "CD4_Treg.dpt_pseudotime.txt.gz",
                "CD4_Tfh.dpt_pseudotime.txt.gz",
                "CD4_Treg.monocle3_pseudotime.txt.gz",
                "CD4_Tfh.monocle3_pseudotime.txt.gz",
            ):
                f_obj = tf.extractfile(member)
                if f_obj:
                    col_key = "monocle3_pseudotime" if "monocle3" in bname else "dpt_pseudotime"
                    dec = gzip.decompress(f_obj.read()).decode("utf-8", errors="replace")
                    for ln in dec.splitlines()[1:]:
                        parts = ln.split("\t")
                        if len(parts) >= 2 and parts[0].strip():
                            try:
                                mc_id = parts[0].strip()
                                val = float(parts[-1])
                                pt_records.setdefault(mc_id, {})[col_key] = val
                            except ValueError:
                                pass

            elif bname == "signature_genes.txt.gz":
                f_obj = tf.extractfile(member)
                if f_obj:
                    dec = gzip.decompress(f_obj.read()).decode("utf-8", errors="replace")
                    for ln in dec.splitlines():
                        parts = ln.strip().split("\t")
                        if len(parts) >= 2 and parts[0] and parts[1]:
                            sig_genes.setdefault(parts[0], []).append(parts[1])

    if pt_records:
        mc_keys = sorted(pt_records.keys())
        pseudotime_maybe = Some(
            pl.DataFrame(
                {
                    "miniCluster": mc_keys,
                    "dpt_pseudotime": [
                        pt_records[k].get("dpt_pseudotime", float("nan")) for k in mc_keys
                    ],
                    "monocle3_pseudotime": [
                        pt_records[k].get("monocle3_pseudotime", float("nan")) for k in mc_keys
                    ],
                }
            )
        )

    return (clinical_maybe, pseudotime_maybe, metacell_umap, sig_genes)


def _extract_zheng_tcr_polars(
    tcr_tar_path: Path,
    scratch_dir: Path,
) -> Maybe[pl.DataFrame]:
    """Extract cell-level TCR clonotype & pTRT tumor-reactivity scores from `data.tcr.tar.gz`."""
    if not tcr_tar_path.is_file():
        return Nothing

    scratch_dir.mkdir(parents=True, exist_ok=True)
    tmp_comb = scratch_dir / "_zheng_tcr_comb.rds"
    tmp_ptrt = scratch_dir / "_zheng_tcr_ptrt.rds"
    frames: list[pl.DataFrame] = []

    try:
        with tarfile.open(tcr_tar_path, "r:gz") as tf:
            names = {m.name: m for m in tf.getmembers() if m.isfile()}
            if "data/tcr/pTRT/TRT.signaling.tcr.ext.sc.tb.rds" in names:
                f_ptrt = tf.extractfile(names["data/tcr/pTRT/TRT.signaling.tcr.ext.sc.tb.rds"])
                if f_ptrt:
                    tmp_ptrt.write_bytes(f_ptrt.read())
                    match read_r_dataframe(tmp_ptrt):
                        case Success(pdf):
                            frames.append(
                                pl.DataFrame(
                                    {
                                        "source_study": [str(x).strip() for x in pdf["dataset"].tolist()],
                                        "_join_cell_name": [str(x).strip() for x in pdf["Cell_Name"].tolist()],
                                        "tcr_clone_id": [str(x).strip() for x in pdf["cloneID"].tolist()],
                                        "tcr_clone_status": [str(x).strip() for x in pdf["cloneStatus"].tolist()],
                                        "tcr_clone_size": pd.to_numeric(pdf["cloneSize"], errors="coerce").fillna(0).to_numpy(dtype=np.int64),
                                        "tcr_cdr3_alpha": [str(x).strip() for x in pdf["CDR3.A1"].tolist()],
                                        "tcr_cdr3_beta": [str(x).strip() for x in pdf["CDR3.B1"].tolist()],
                                        "ptrt_nes": pd.to_numeric(pdf["NES"], errors="coerce").to_numpy(dtype=np.float64),
                                        "ptrt_pvalue": pd.to_numeric(pdf["pvalue"], errors="coerce").to_numpy(dtype=np.float64),
                                    }
                                )
                            )
                        case _:
                            pass

            if "data/tcr/byCell/tcr.zhangLab.comb.flt.rds" in names:
                f_comb = tf.extractfile(names["data/tcr/byCell/tcr.zhangLab.comb.flt.rds"])
                if f_comb:
                    tmp_comb.write_bytes(f_comb.read())
                    match read_r_dataframe(tmp_comb):
                        case Success(pdf2):
                            frames.append(
                                pl.DataFrame(
                                    {
                                        "source_study": [str(x).strip() for x in pdf2["dataset"].tolist()],
                                        "_join_cell_name": [str(x).strip() for x in pdf2["Cell_Name"].tolist()],
                                        "tcr_clone_id": [str(x).strip() for x in pdf2["cloneID"].tolist()],
                                        "tcr_clone_status": [str(x).strip() for x in pdf2["cloneStatus"].tolist()],
                                        "tcr_clone_size": pd.to_numeric(pdf2["cloneSize"], errors="coerce").fillna(0).to_numpy(dtype=np.int64),
                                        "tcr_cdr3_alpha": [str(x).strip() for x in pdf2["CDR3.A1"].tolist()],
                                        "tcr_cdr3_beta": [str(x).strip() for x in pdf2["CDR3.B1"].tolist()],
                                        "ptrt_nes": np.full(len(pdf2), np.nan, dtype=np.float64),
                                        "ptrt_pvalue": np.full(len(pdf2), np.nan, dtype=np.float64),
                                    }
                                )
                            )
                        case _:
                            pass
    finally:
        tmp_comb.unlink(missing_ok=True)
        tmp_ptrt.unlink(missing_ok=True)

    if not frames:
        return Nothing
    merged = pl.concat(frames, how="vertical").unique(
        subset=["source_study", "_join_cell_name"],
        keep="first",
    )
    return Some(merged)


def download_zheng_tcell_atlas(
    raw_dir: Maybe[Path] = Nothing,
    force_download: bool = False,
    verify_md5: bool = True,
    repo_root: Maybe[Path] = Nothing,
) -> Result[Path, str]:
    """Download Zheng et al. (2021, Science) Zenodo `5461803` archives (`metaInfo`, `tcr`, `expression`)."""
    if is_cscs_environment():
        return Failure(
            "Internet dataset downloads are strictly disabled on CSCS Eiger. "
            "Download Zheng_TCell_Atlas on 'olm' first and stage via Globus."
        )

    root = repo_root.unwrap().resolve() if isinstance(repo_root, Some) else find_repo_root()
    target_dir = (
        raw_dir.unwrap().resolve()
        if isinstance(raw_dir, Some)
        else get_raw_dataset_dir(ZHENG_TCELL_DATASET_ID, repo_root=root)
    )
    target_dir.mkdir(parents=True, exist_ok=True)

    jobs: tuple[tuple[str, str, ChecksumSpec], ...] = (
        (ZHENG_TCELL_META_URL, ZHENG_TCELL_META_TAR_FILENAME, ZHENG_TCELL_META_CHECKSUM),
        (ZHENG_TCELL_TCR_URL, ZHENG_TCELL_TCR_TAR_FILENAME, ZHENG_TCELL_TCR_CHECKSUM),
        (ZHENG_TCELL_EXPR_URL, ZHENG_TCELL_EXPR_TAR_FILENAME, ZHENG_TCELL_EXPR_CHECKSUM),
    )

    for url, fname, spec in jobs:
        dest = target_dir / fname
        if force_download and dest.exists():
            dest.unlink(missing_ok=True)
        csum: Maybe[ChecksumSpec] = (
            Some(spec)
            if verify_md5
            else Some(ChecksumSpec(algorithm="size_only", expected_hash="", size_bytes=spec.size_bytes))
        )
        match download_single_file(
            url=url,
            dest_path=dest,
            timeout_sec=180.0,
            max_retries=15,
            resume=True,
            checksum=csum,
        ):
            case Failure(err):
                return Failure(f"Zheng_TCell_Atlas download failed for {fname}: {err}")
            case Success(p):
                if verify_md5:
                    _ = _write_md5_sidecar(p, spec.expected_hash)
            case _:
                return Failure(f"Unexpected download state for {fname}")

    return Success(target_dir)


def _preprocess_zheng_from_geo_counts(
    raw_dir: Path,
    count_files: Sequence[Path],
    target_h5ad: Path,
    normalize_ensembl: bool,
    ensembl_release: int,
    drop_unmapped: bool,
) -> Result[Path, str]:
    """Fallback preprocessor for GEO `GSE156728` `.counts.txt.gz` + `GSE156728_metadata.txt.gz` files."""
    meta_gz = raw_dir / "GSE156728_metadata.txt.gz"
    meta_lookup: dict[str, dict[str, Any]] = {}
    if meta_gz.is_file():
        with gzip.open(meta_gz, "rt", encoding="utf-8", errors="replace") as mf:
            meta_pdf = pd.read_csv(mf, sep="\t", dtype=str).fillna("NA")
        key_col = "cellID.uniq" if "cellID.uniq" in meta_pdf.columns else meta_pdf.columns[0]
        for rec in meta_pdf.to_dict(orient="records"):
            rec_str: dict[str, Any] = {str(k): v for k, v in rec.items()}
            meta_lookup[str(rec_str.get(key_col, ""))] = rec_str
            if "cellID" in rec_str:
                meta_lookup[str(rec_str["cellID"])] = rec_str

    cohort_adatas: list[ad.AnnData] = []
    for cf in sorted(count_files):
        with gzip.open(cf, "rt", encoding="utf-8", errors="replace") as gf:
            mat_df = pd.read_csv(gf, sep="\t", index_col=0)
        genes = [str(g) for g in mat_df.index.tolist()]
        cells = [str(c) for c in mat_df.columns.tolist()]
        counts_csr = sp.csr_matrix(mat_df.to_numpy(dtype=np.float32).T)
        rows_meta = [meta_lookup.get(cid, {"cellID": cid}) for cid in cells]
        obs_raw = pd.DataFrame(rows_meta, index=pd.Index(cells, name="cell_id"))
        obs_clean = harmonize_zheng_tcell_obs_polars(obs_raw, cell_names_fallback=cells)
        cohort_adatas.append(
            ad.AnnData(
                X=_normalize_counts_csr_to_log1p_cp10k(counts_csr),
                obs=obs_clean,
                var=pd.DataFrame(index=pd.Index(genes, name="original_id")),
                layers={"counts": counts_csr},
            )
        )

    merged = ad.concat(cohort_adatas, axis=0, join="outer", fill_value=0.0)
    merged.obs_names_make_unique()
    raw_genes = [str(g) for g in merged.var_names.tolist()]
    if normalize_ensembl:
        var_final, keep_mask, proj_m, unmapped = _build_ensembl_projection_operator(
            raw_var_names=raw_genes,
            raw_var=pd.DataFrame(index=pd.Index(raw_genes, name="original_id")),
            release=ensembl_release,
            species="human",
            ensembl_dir=get_ensembl_dir(),
            drop_unmapped=drop_unmapped,
            aggregation="sum",
        )
        proj_counts = _project_sparse_chunk_to_ensembl(
            sp.csr_matrix(merged.layers["counts"]),
            keep_mask,
            proj_m,
            "sum",
        )
        proj_log = _normalize_counts_csr_to_log1p_cp10k(proj_counts)
    else:
        unmapped = []
        var_final = pd.DataFrame(
            {"gene_name": raw_genes, "original_id": raw_genes},
            index=ad.utils.make_index_unique(pd.Index(raw_genes, name="gene_id")),
        )
        proj_counts = sp.csr_matrix(merged.layers["counts"], dtype=np.float32)
        proj_log = _normalize_counts_csr_to_log1p_cp10k(proj_counts)

    out_ad = ad.AnnData(
        X=proj_log,
        obs=merged.obs,
        var=var_final,
        layers={"counts": proj_counts, "log1p_norm": proj_log},
        uns={
            "dataset_id": ZHENG_TCELL_DATASET_ID,
            "title": "Zheng et al. 2021 Pan-Cancer Tumor-Infiltrating T-Cell Atlas",
            "publication_doi": ZHENG_TCELL_PUBLICATION_DOI,
            "publication_pmid": ZHENG_TCELL_PUBLICATION_PMID,
            "is_raw_counts": False,
            "has_raw_counts": True,
            "expression_type": "log_normalized",
            "n_unmapped_genes": len(unmapped),
            "atlas_metadata": {
                "journal": "Science",
                "year": 2021,
                "source_mode": "GEO_GSE156728",
                "cluster_subtypes": {k: v.short_name for k, v in ZHENG_TCELL_CLUSTER_SPECS.items()},
            },
        },
    )
    target_h5ad.parent.mkdir(parents=True, exist_ok=True)
    out_ad.write_h5ad(target_h5ad, compression="gzip")
    return Success(target_h5ad)


def preprocess_zheng_tcell_atlas(
    raw_dir_or_tar: Path,
    output_h5ad: Maybe[Path] = Nothing,
    normalize_ensembl: bool = True,
    ensembl_release: int = 111,
    drop_unmapped: bool = True,
    batch_size: int = 10_000,
) -> Result[Path, str]:
    """Stream Zheng et al. (2021, Science) Zenodo `5461803` (397,810 cells) or GEO `GSE156728` into `.h5ad`.

    Executes a single sequential streaming pass over `data.expression.tar.gz` (`CD8/byDataset/*.sce.rds`
    and `CD4/byDataset/*.sce.rds`), projecting every cell's full transcriptome onto canonical Ensembl v111
    with < 500 MB temporary disk overhead and < 1.5 GB RAM, while enriching `.obs` and `.obsm` with
    clinical Stage/Age/Gender, TCR clonality & pTRT scores, Harmony UMAP, cohort UMAP, and pseudotime.
    """
    raw_dir = raw_dir_or_tar if raw_dir_or_tar.is_dir() else raw_dir_or_tar.parent
    expr_tar = (
        raw_dir_or_tar
        if raw_dir_or_tar.is_file() and raw_dir_or_tar.name.endswith(".tar.gz")
        else raw_dir / ZHENG_TCELL_EXPR_TAR_FILENAME
    )
    meta_tar = raw_dir / ZHENG_TCELL_META_TAR_FILENAME
    tcr_tar = raw_dir / ZHENG_TCELL_TCR_TAR_FILENAME

    target_h5ad = (
        output_h5ad.unwrap().resolve()
        if isinstance(output_h5ad, Some)
        else get_preprocessed_h5ad_path(ZHENG_TCELL_DATASET_ID).resolve()
    )
    target_h5ad.parent.mkdir(parents=True, exist_ok=True)

    if not expr_tar.is_file():
        geo_counts = sorted(raw_dir.glob("*counts.txt.gz"))
        if geo_counts:
            return _preprocess_zheng_from_geo_counts(
                raw_dir=raw_dir,
                count_files=geo_counts,
                target_h5ad=target_h5ad,
                normalize_ensembl=normalize_ensembl,
                ensembl_release=ensembl_release,
                drop_unmapped=drop_unmapped,
            )
        return Failure(
            f"Zheng_TCell_Atlas expression archive not found at {expr_tar} "
            f"(nor any *counts.txt.gz in {raw_dir})."
        )

    scratch_dir = get_scratch_dataset_dir(ZHENG_TCELL_DATASET_ID)
    scratch_dir.mkdir(parents=True, exist_ok=True)
    tmp_h5ad = target_h5ad.with_name(target_h5ad.name + ".tmp")
    tmp_rds = scratch_dir / "_zheng_cohort_stream.sce.rds"
    tmp_stage_h5 = scratch_dir / "_zheng_cohort_stream.h5"

    try:
        clinical_maybe, pseudotime_maybe, metacell_umap, sig_genes = _extract_zheng_metainfo_artifacts(
            meta_tar_path=meta_tar,
            scratch_dir=scratch_dir,
        )
        tcr_maybe = _extract_zheng_tcr_polars(tcr_tar_path=tcr_tar, scratch_dir=scratch_dir)

        var_final, id_to_col = _build_canonical_ensembl_var_and_lookup(ensembl_release=ensembl_release)
        n_final_vars = len(var_final)

        obs_frames: list[pd.DataFrame] = []
        cohort_umap_parts: list[np.ndarray[Any, Any]] = []
        streamed_cohorts: list[str] = []
        has_any_counts = False
        total_cells = 0

        tmp_h5ad.unlink(missing_ok=True)
        with h5py.File(tmp_h5ad, "w") as dst_h5:
            dst_h5.attrs["encoding-type"] = "anndata"
            dst_h5.attrs["encoding-version"] = "0.1.0"
            write_anndata_dataframe_elem(dst_h5, "var", var_final, default_index_name="gene_id")

            layers_grp = dst_h5.create_group("layers")
            layers_grp.attrs["encoding-type"] = "dict"
            layers_grp.attrs["encoding-version"] = "0.1.0"

            grp_x = init_compressed_csr_group(dst_h5, "X", n_obs=400_000, n_vars=n_final_vars)
            grp_log1p = init_compressed_csr_group(
                layers_grp,
                "log1p_norm",
                n_obs=400_000,
                n_vars=n_final_vars,
            )
            grp_counts = init_compressed_csr_group(
                layers_grp,
                "counts",
                n_obs=400_000,
                n_vars=n_final_vars,
            )

            with tarfile.open(expr_tar, mode="r|gz") as tf:
                for member in tf:
                    if not member.isfile():
                        continue
                    mname = member.name

                    # Capture CD8/CD4 integrated metacell harmony.umap if encountered in stream
                    if mname.endswith(("int.CD8.S35.sce.merged.rds", "int.CD4.S35.sce.merged.rds")):
                        f_int = tf.extractfile(member)
                        if f_int:
                            with open(tmp_rds, "wb") as out_f:
                                shutil.copyfileobj(f_int, out_f, length=8 * 1024 * 1024)
                            match stream_r_object_to_staging_h5(
                                r_path=tmp_rds,
                                staging_h5_path=tmp_stage_h5,
                                skip_assay_counts=True,
                            ):
                                case Success(int_man):
                                    red_key = (
                                        "harmony.umap"
                                        if "harmony.umap" in int_man.reductions
                                        else ("umap" if "umap" in int_man.reductions else "")
                                    )
                                    if red_key and not int_man.obs_df.empty:
                                        with h5py.File(tmp_stage_h5, "r") as int_h5:
                                            coords = np.asarray(
                                                int_h5[int_man.reductions[red_key].h5_dataset][:, :2],
                                                dtype=np.float32,
                                            )
                                        for mc_id, row_xy in zip(int_man.obs_df.index.tolist(), coords):
                                            metacell_umap[str(mc_id).strip()] = row_xy
                                case _:
                                    pass
                            tmp_rds.unlink(missing_ok=True)
                            tmp_stage_h5.unlink(missing_ok=True)
                        continue

                    is_by_dataset = (
                        ("/CD8/byDataset/" in mname or "/CD4/byDataset/" in mname)
                        and mname.endswith(".sce.rds")
                    )
                    if not is_by_dataset:
                        continue

                    f_sce = tf.extractfile(member)
                    if not f_sce:
                        continue
                    with open(tmp_rds, "wb") as out_f:
                        shutil.copyfileobj(f_sce, out_f, length=8 * 1024 * 1024)

                    match stream_r_object_to_staging_h5(
                        r_path=tmp_rds,
                        staging_h5_path=tmp_stage_h5,
                        skip_assay_counts=False,
                    ):
                        case Failure(st_err):
                            logger.warning("Skipping %s due to R stream error: %s", mname, st_err)
                            tmp_rds.unlink(missing_ok=True)
                            tmp_stage_h5.unlink(missing_ok=True)
                            continue
                        case Success(man):
                            pass
                        case _:
                            continue

                    if not man.assays:
                        tmp_rds.unlink(missing_ok=True)
                        tmp_stage_h5.unlink(missing_ok=True)
                        continue

                    has_counts_assay = "counts" in man.assays
                    primary_key = (
                        "counts"
                        if has_counts_assay
                        else (
                            "norm_exprs"
                            if "norm_exprs" in man.assays
                            else next(iter(man.assays.keys()))
                        )
                    )
                    a_spec = man.assays[primary_key]
                    n_c = a_spec.n_cells
                    if n_c == 0:
                        tmp_rds.unlink(missing_ok=True)
                        tmp_stage_h5.unlink(missing_ok=True)
                        continue

                    keep_mask, proj_m, id_to_col = _build_cohort_projection_to_canonical(
                        cohort_genes=a_spec.gene_names,
                        n_canonical_vars=n_final_vars,
                        id_to_col=id_to_col,
                        ensembl_release=ensembl_release,
                    )

                    obs_part = harmonize_zheng_tcell_obs_polars(
                        raw_obs=man.obs_df,
                        cell_names_fallback=a_spec.cell_names,
                        clinical_df=clinical_maybe,
                        tcr_df=tcr_maybe,
                        pseudotime_df=pseudotime_maybe,
                    )
                    if len(obs_part) != n_c:
                        obs_part = obs_part.iloc[:n_c].copy()
                    obs_frames.append(obs_part)

                    with h5py.File(tmp_stage_h5, "r") as src_h5:
                        umap_cand = next(
                            (rk for rk in ("seurat.umap", "umap", "harmony.umap") if rk in man.reductions),
                            "",
                        )
                        if umap_cand and man.reductions[umap_cand].n_cells == n_c:
                            cohort_umap_parts.append(
                                np.asarray(src_h5[man.reductions[umap_cand].h5_dataset][:, :2], dtype=np.float32)
                            )
                        else:
                            cohort_umap_parts.append(np.full((n_c, 2), np.nan, dtype=np.float32))

                        src_grp = src_h5[a_spec.h5_group]
                        for start in range(0, n_c, batch_size):
                            end = min(start + batch_size, n_c)
                            raw_chunk = (
                                _slice_staged_csr_batch(src_grp, start, end, a_spec.n_genes)
                                if a_spec.is_sparse
                                else sp.csr_matrix(src_grp["dense"][start:end, :], dtype=np.float32)
                            )
                            if has_counts_assay:
                                has_any_counts = True
                                proj_counts = _project_sparse_chunk_to_ensembl(
                                    raw_chunk,
                                    keep_mask,
                                    proj_m,
                                    "sum",
                                )
                                proj_log = _normalize_counts_csr_to_log1p_cp10k(proj_counts)
                                append_csr_batch_to_group(grp_counts, proj_counts)
                                append_csr_batch_to_group(grp_x, proj_log)
                                append_csr_batch_to_group(grp_log1p, proj_log)
                            else:
                                proj_log = _project_sparse_chunk_to_ensembl(
                                    raw_chunk,
                                    keep_mask,
                                    proj_m,
                                    "log1p_sum",
                                )
                                empty_counts = sp.csr_matrix(
                                    (
                                        np.empty(0, dtype=np.float32),
                                        (np.empty(0, dtype=np.int32), np.empty(0, dtype=np.int32)),
                                    ),
                                    shape=(end - start, n_final_vars),
                                    dtype=np.float32,
                                )
                                append_csr_batch_to_group(grp_counts, empty_counts)
                                append_csr_batch_to_group(grp_x, proj_log)
                                append_csr_batch_to_group(grp_log1p, proj_log)

                    total_cells += n_c
                    streamed_cohorts.append(Path(mname).name.replace(".sce.rds", ""))
                    tmp_rds.unlink(missing_ok=True)
                    tmp_stage_h5.unlink(missing_ok=True)
                    gc.collect()

            if total_cells == 0 or not obs_frames:
                return Failure(f"No valid SingleCellExperiment cohorts streamed from {expr_tar.name}")

            grp_x.attrs["shape"] = np.array([total_cells, n_final_vars], dtype=np.int64)
            grp_log1p.attrs["shape"] = np.array([total_cells, n_final_vars], dtype=np.int64)
            grp_counts.attrs["shape"] = np.array([total_cells, n_final_vars], dtype=np.int64)
            if not has_any_counts:
                _ = layers_grp.pop("counts", False)

            combined_obs = pd.concat(obs_frames, axis=0)
            unique_idx = ad.utils.make_index_unique(pd.Index([str(x) for x in combined_obs.index], name="cell_id"))
            combined_obs["cell_id"] = [str(x) for x in unique_idx]
            combined_obs.index = unique_idx
            write_anndata_dataframe_elem(dst_h5, "obs", combined_obs, default_index_name="cell_id")

            obsm_grp = dst_h5.create_group("obsm")
            obsm_grp.attrs["encoding-type"] = "dict"
            obsm_grp.attrs["encoding-version"] = "0.1.0"

            cohort_umap = np.vstack(cohort_umap_parts).astype(np.float32, copy=False)
            anndata.io.write_elem(obsm_grp, "X_umap_cohort", cohort_umap)

            if metacell_umap and "miniCluster" in combined_obs.columns:
                harmony_xy = np.full((total_cells, 2), np.nan, dtype=np.float32)
                for row_i, mc_tok in enumerate(combined_obs["miniCluster"].astype(str).tolist()):
                    if mc_tok in metacell_umap:
                        harmony_xy[row_i] = metacell_umap[mc_tok]
                anndata.io.write_elem(obsm_grp, "X_umap", harmony_xy)
                anndata.io.write_elem(obsm_grp, "X_umap_harmony", harmony_xy)
            else:
                anndata.io.write_elem(obsm_grp, "X_umap", cohort_umap)

            uns_payload: dict[str, Any] = {
                "dataset_id": ZHENG_TCELL_DATASET_ID,
                "title": "Zheng et al. 2021 Pan-Cancer Tumor-Infiltrating T-Cell Atlas",
                "publication_doi": ZHENG_TCELL_PUBLICATION_DOI,
                "publication_pmid": ZHENG_TCELL_PUBLICATION_PMID,
                "zenodo_record": ZHENG_TCELL_ZENODO_RECORD,
                "zenodo_doi": ZHENG_TCELL_ZENODO_DOI,
                "geo_accession": ZHENG_TCELL_GEO_ACCESSION,
                "is_raw_counts": False,
                "has_raw_counts": has_any_counts,
                "expression_type": "log_normalized",
                "atlas_metadata": {
                    "journal": "Science",
                    "year": 2021,
                    "n_cohorts_streamed": len(streamed_cohorts),
                    "cohorts_streamed": streamed_cohorts,
                    "cluster_subtypes": {
                        cid: spec.short_name for cid, spec in ZHENG_TCELL_CLUSTER_SPECS.items()
                    },
                    "cluster_functional_states": {
                        cid: spec.functional_state for cid, spec in ZHENG_TCELL_CLUSTER_SPECS.items()
                    },
                    "author_signature_genes": {
                        cid: genes_list[:50] for cid, genes_list in sig_genes.items()
                    },
                },
            }
            anndata.io.write_elem(dst_h5, "uns", uns_payload)

        if target_h5ad.exists():
            target_h5ad.unlink()
        tmp_h5ad.rename(target_h5ad)
        logger.info(
            "Successfully wrote Zheng_TCell_Atlas H5AD to %s (%d cells x %d genes across %d cohorts)",
            target_h5ad,
            total_cells,
            n_final_vars,
            len(streamed_cohorts),
        )
        return Success(target_h5ad)
    except Exception as exc:  # noqa: BLE001
        tmp_h5ad.unlink(missing_ok=True)
        return Failure(f"Failed to preprocess Zheng_TCell_Atlas: {exc}")
    finally:
        tmp_rds.unlink(missing_ok=True)
        tmp_stage_h5.unlink(missing_ok=True)


def load_zheng_tcell_atlas(
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
    """Load the Zheng et al. (2021, Science) Pan-Cancer Tumor-Infiltrating T-Cell Atlas (`Zheng_TCell_Atlas`)."""
    root = repo_root.unwrap().resolve() if isinstance(repo_root, Some) else find_repo_root()
    resolved_raw_dir = (
        raw_dir.unwrap().resolve()
        if isinstance(raw_dir, Some)
        else get_raw_dataset_dir(ZHENG_TCELL_DATASET_ID, repo_root=root)
    )
    target_h5ad = (
        output_h5ad.unwrap().resolve()
        if isinstance(output_h5ad, Some)
        else get_preprocessed_h5ad_path(ZHENG_TCELL_DATASET_ID, repo_root=root)
    )

    cached_h5ad: Maybe[Path] = Nothing
    if not force_recompute and not force_download:
        cached_h5ad = (
            Some(target_h5ad)
            if target_h5ad.is_file() and target_h5ad.stat().st_size > 0
            else find_dataset_h5ad(ZHENG_TCELL_DATASET_ID, repo_root=root)
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

    expr_tar = resolved_raw_dir / ZHENG_TCELL_EXPR_TAR_FILENAME
    has_geo_fallback = bool(list(resolved_raw_dir.glob("*counts.txt.gz")))

    if force_download or (not expr_tar.is_file() and not has_geo_fallback):
        if not auto_download:
            return Failure(
                f"Zheng_TCell_Atlas raw archive not found at {expr_tar} and auto_download=False. "
                f"Download from Zenodo {ZHENG_TCELL_ZENODO_RECORD} or set auto_download=True."
            )
        match download_zheng_tcell_atlas(
            raw_dir=Some(resolved_raw_dir),
            force_download=force_download,
            verify_md5=True,
            repo_root=Some(root),
        ):
            case Failure(dl_err):
                return Failure(dl_err)
            case Success(_):
                pass
            case _:
                return Failure("Unexpected Zheng_TCell_Atlas download state")

    match preprocess_zheng_tcell_atlas(
        raw_dir_or_tar=resolved_raw_dir,
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
            return Failure("Unexpected Zheng_TCell_Atlas preprocess state")


# -----------------------------------------------------------------------------
# 3B. Chu et al. (2023, Nature Medicine) — `Chu_TCell_Atlas`
# -----------------------------------------------------------------------------


def _discover_chu_staged_rds_files(
    raw_dir: Path,
    manual_dir: Path,
) -> list[Path]:
    """Locate any staged Chu et al. (2023) `.rds`/`.rda` files in `raw_dir` or `data/manual_download/Chu_TCell_Atlas`."""
    found: list[Path] = []
    for search_root in (raw_dir, manual_dir / CHU_TCELL_DATASET_ID):
        if not search_root.is_dir():
            continue
        for cand in CHU_TCELL_CANDIDATE_FILENAMES:
            p = search_root / cand
            if p.is_file() and p.stat().st_size > 0:
                found.append(p)
        if not found:
            for p_any in sorted(search_root.glob("*.rds")) + sorted(search_root.glob("*.RDS")):
                if p_any.is_file() and p_any.stat().st_size > 0:
                    found.append(p_any)
        if found:
            return found
    return found


def download_chu_tcell_atlas(
    raw_dir: Maybe[Path] = Nothing,
    force_download: bool = False,
    repo_root: Maybe[Path] = Nothing,
) -> Result[Path, str]:
    """Attempt portal download of Chu et al. (2023, Nature Medicine) T-cell atlas or resolve manual staging."""
    if is_cscs_environment():
        return Failure(
            "Internet dataset downloads are strictly disabled on CSCS Eiger. "
            "Stage Chu_TCell_Atlas (`panC_T_all.rds` or `CD8.rds` + `CD4.rds`) via Globus."
        )

    root = repo_root.unwrap().resolve() if isinstance(repo_root, Some) else find_repo_root()
    target_dir = (
        raw_dir.unwrap().resolve()
        if isinstance(raw_dir, Some)
        else get_raw_dataset_dir(CHU_TCELL_DATASET_ID, repo_root=root)
    )
    manual_root = get_manual_download_dir(repo_root=root)

    if not force_download:
        staged = _discover_chu_staged_rds_files(target_dir, manual_root)
        if staged:
            return Success(staged[0])

    target_dir.mkdir(parents=True, exist_ok=True)
    dest_rds = target_dir / "panC_T_all.rds"
    last_err = "Portal connection unavailable"
    for cand_url in CHU_TCELL_CANDIDATE_URLS:
        match download_single_file(
            url=cand_url,
            dest_path=dest_rds,
            timeout_sec=15.0,
            max_retries=1,
            resume=True,
            checksum=Nothing,
        ):
            case Success(downloaded_path) if downloaded_path.is_file() and downloaded_path.stat().st_size > 1_048_576:
                return Success(downloaded_path)
            case Failure(err):
                last_err = err
            case _:
                dest_rds.unlink(missing_ok=True)

    dest_rds.unlink(missing_ok=True)
    manual_target = manual_root / CHU_TCELL_DATASET_ID
    return Failure(
        f"Chu_TCell_Atlas automated fetch from MD Anderson TCM portal ({CHU_TCELL_PORTAL_URL}) failed "
        f"({last_err}). Download 'panC_T_all.rds' (or 'CD8.rds' and 'CD4.rds') from "
        f"{CHU_TCELL_PORTAL_URL} (Chu et al., Nature Medicine 2023, DOI: {CHU_TCELL_PUBLICATION_DOI}) "
        f"and stage into '{manual_target}/'."
    )


def preprocess_chu_tcell_atlas(
    rds_path_or_dir: Path,
    output_h5ad: Maybe[Path] = Nothing,
    normalize_ensembl: bool = True,
    ensembl_release: int = 111,
    drop_unmapped: bool = True,
    batch_size: int = 10_000,
) -> Result[Path, str]:
    """Stream Chu et al. (2023, Nature Medicine) Seurat `.rds` file(s) into an Ensembl-normalized `.h5ad`."""
    rds_files: list[Path] = []
    if rds_path_or_dir.is_file():
        rds_files = [rds_path_or_dir]
    elif rds_path_or_dir.is_dir():
        rds_files = _discover_chu_staged_rds_files(rds_path_or_dir, rds_path_or_dir.parent)

    if not rds_files:
        return Failure(
            f"No Chu_TCell_Atlas .rds files found at {rds_path_or_dir}. "
            f"Stage 'panC_T_all.rds' (or 'CD8.rds' + 'CD4.rds') from {CHU_TCELL_PORTAL_URL}."
        )

    target_h5ad = (
        output_h5ad.unwrap().resolve()
        if isinstance(output_h5ad, Some)
        else get_preprocessed_h5ad_path(CHU_TCELL_DATASET_ID).resolve()
    )
    target_h5ad.parent.mkdir(parents=True, exist_ok=True)
    tmp_h5ad = target_h5ad.with_name(target_h5ad.name + ".tmp")

    scratch_dir = get_scratch_dataset_dir(CHU_TCELL_DATASET_ID)
    scratch_dir.mkdir(parents=True, exist_ok=True)
    staging_h5 = scratch_dir / "_chu_tcell_staging.h5"

    try:
        if normalize_ensembl:
            var_final, id_to_col = _build_canonical_ensembl_var_and_lookup(ensembl_release=ensembl_release)
        else:
            # First-file fallback when `normalize_ensembl=False` (synthetic unit tests)
            match stream_r_object_to_staging_h5(rds_files[0], staging_h5, skip_assay_counts=True):
                case Failure(err0):
                    return Failure(err0)
                case Success(man0):
                    ak0 = "RNA" if "RNA" in man0.assays else next(iter(man0.assays.keys()), "")
                    if not ak0:
                        return Failure(f"No assay in {rds_files[0].name}")
                    g0 = list(man0.assays[ak0].gene_names)
                    var_final = pd.DataFrame(
                        {"gene_name": g0, "original_id": g0},
                        index=ad.utils.make_index_unique(pd.Index(g0, name="gene_id")),
                    )
                    id_to_col = {g: i for i, g in enumerate(g0)}
                case _:
                    return Failure("Unexpected Chu staging state")
            staging_h5.unlink(missing_ok=True)

        n_final_vars = len(var_final)
        obs_frames: list[pd.DataFrame] = []
        umap_parts: list[np.ndarray[Any, Any]] = []
        total_cells = 0

        tmp_h5ad.unlink(missing_ok=True)
        with h5py.File(tmp_h5ad, "w") as dst_h5:
            dst_h5.attrs["encoding-type"] = "anndata"
            dst_h5.attrs["encoding-version"] = "0.1.0"
            write_anndata_dataframe_elem(dst_h5, "var", var_final, default_index_name="gene_id")

            layers_grp = dst_h5.create_group("layers")
            layers_grp.attrs["encoding-type"] = "dict"
            layers_grp.attrs["encoding-version"] = "0.1.0"

            grp_x = init_compressed_csr_group(dst_h5, "X", n_obs=320_000, n_vars=n_final_vars)
            grp_log1p = init_compressed_csr_group(layers_grp, "log1p_norm", n_obs=320_000, n_vars=n_final_vars)

            for rds_file in rds_files:
                match stream_r_object_to_staging_h5(rds_file, staging_h5, skip_assay_counts=True):
                    case Failure(st_err):
                        return Failure(st_err)
                    case Success(man):
                        pass
                    case _:
                        return Failure("Unexpected Chu staging state")

                akey = "RNA" if "RNA" in man.assays else next(iter(man.assays.keys()), "")
                if not akey:
                    staging_h5.unlink(missing_ok=True)
                    return Failure(f"No RNA assay found in {rds_file.name}")

                aspec = man.assays[akey]
                n_c = aspec.n_cells
                obs_part = harmonize_chu_tcell_obs_polars(
                    raw_obs=man.obs_df,
                    active_ident=man.active_ident,
                    cell_names_fallback=aspec.cell_names,
                )
                obs_frames.append(obs_part)

                if normalize_ensembl:
                    keep_mask, proj_m, id_to_col = _build_cohort_projection_to_canonical(
                        cohort_genes=aspec.gene_names,
                        n_canonical_vars=n_final_vars,
                        id_to_col=id_to_col,
                        ensembl_release=ensembl_release,
                    )
                    agg_mode = "log1p_sum"
                else:
                    keep_mask = np.ones(aspec.n_genes, dtype=bool)
                    proj_m = sp.eye(aspec.n_genes, n=n_final_vars, dtype=np.float32, format="csr")
                    agg_mode = "sum"

                with h5py.File(staging_h5, "r") as src_h5:
                    um_key = next((rk for rk in ("umap", "harmony", "tsne") if rk in man.reductions), "")
                    if um_key and man.reductions[um_key].n_cells == n_c:
                        umap_parts.append(
                            np.asarray(src_h5[man.reductions[um_key].h5_dataset][:, :2], dtype=np.float32)
                        )
                    else:
                        umap_parts.append(np.full((n_c, 2), np.nan, dtype=np.float32))

                    src_grp = src_h5[aspec.h5_group]
                    for start in range(0, n_c, batch_size):
                        end = min(start + batch_size, n_c)
                        raw_chunk = (
                            _slice_staged_csr_batch(src_grp, start, end, aspec.n_genes)
                            if aspec.is_sparse
                            else sp.csr_matrix(src_grp["dense"][start:end, :], dtype=np.float32)
                        )
                        proj_chunk = _project_sparse_chunk_to_ensembl(
                            raw_chunk,
                            keep_mask,
                            proj_m,
                            agg_mode,
                        )
                        append_csr_batch_to_group(grp_x, proj_chunk)
                        append_csr_batch_to_group(grp_log1p, proj_chunk)

                total_cells += n_c
                staging_h5.unlink(missing_ok=True)
                gc.collect()

            grp_x.attrs["shape"] = np.array([total_cells, n_final_vars], dtype=np.int64)
            grp_log1p.attrs["shape"] = np.array([total_cells, n_final_vars], dtype=np.int64)

            combined_obs = pd.concat(obs_frames, axis=0)
            unique_idx = ad.utils.make_index_unique(pd.Index([str(x) for x in combined_obs.index], name="cell_id"))
            combined_obs["cell_id"] = [str(x) for x in unique_idx]
            combined_obs.index = unique_idx
            write_anndata_dataframe_elem(dst_h5, "obs", combined_obs, default_index_name="cell_id")

            obsm_grp = dst_h5.create_group("obsm")
            obsm_grp.attrs["encoding-type"] = "dict"
            obsm_grp.attrs["encoding-version"] = "0.1.0"
            anndata.io.write_elem(obsm_grp, "X_umap", np.vstack(umap_parts).astype(np.float32, copy=False))

            uns_payload: dict[str, Any] = {
                "dataset_id": CHU_TCELL_DATASET_ID,
                "title": "Chu et al. 2023 Pan-Cancer T-Cell Stress Response Atlas",
                "publication_doi": CHU_TCELL_PUBLICATION_DOI,
                "publication_pmid": CHU_TCELL_PUBLICATION_PMID,
                "portal_url": CHU_TCELL_PORTAL_URL,
                "is_raw_counts": False,
                "has_raw_counts": False,
                "expression_type": "log_normalized",
                "atlas_metadata": {
                    "journal": "Nature Medicine",
                    "year": 2023,
                    "hs_stress_markers": list(CHU_TCELL_HS_MARKERS),
                    "cluster_subtypes": {
                        cid: spec.short_name for cid, spec in CHU_TCELL_CLUSTER_SPECS.items()
                    },
                    "cluster_functional_states": {
                        cid: spec.functional_state for cid, spec in CHU_TCELL_CLUSTER_SPECS.items()
                    },
                },
            }
            anndata.io.write_elem(dst_h5, "uns", uns_payload)

        if target_h5ad.exists():
            target_h5ad.unlink()
        tmp_h5ad.rename(target_h5ad)
        return Success(target_h5ad)
    except Exception as exc:  # noqa: BLE001
        tmp_h5ad.unlink(missing_ok=True)
        return Failure(f"Failed to preprocess Chu_TCell_Atlas: {exc}")
    finally:
        staging_h5.unlink(missing_ok=True)


def load_chu_tcell_atlas(
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
    """Load the Chu et al. (2023, Nature Medicine) Pan-Cancer T-Cell Atlas (`Chu_TCell_Atlas`)."""
    root = repo_root.unwrap().resolve() if isinstance(repo_root, Some) else find_repo_root()
    resolved_raw_dir = (
        raw_dir.unwrap().resolve()
        if isinstance(raw_dir, Some)
        else get_raw_dataset_dir(CHU_TCELL_DATASET_ID, repo_root=root)
    )
    manual_root = get_manual_download_dir(repo_root=root)
    target_h5ad = (
        output_h5ad.unwrap().resolve()
        if isinstance(output_h5ad, Some)
        else get_preprocessed_h5ad_path(CHU_TCELL_DATASET_ID, repo_root=root)
    )

    cached_h5ad: Maybe[Path] = Nothing
    if not force_recompute and not force_download:
        cached_h5ad = (
            Some(target_h5ad)
            if target_h5ad.is_file() and target_h5ad.stat().st_size > 0
            else find_dataset_h5ad(CHU_TCELL_DATASET_ID, repo_root=root)
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

    staged_files = _discover_chu_staged_rds_files(resolved_raw_dir, manual_root)
    if force_download or not staged_files:
        if not auto_download:
            manual_target = manual_root / CHU_TCELL_DATASET_ID
            return Failure(
                f"Chu_TCell_Atlas raw Seurat .rds not found in '{resolved_raw_dir}' or '{manual_target}' "
                f"and auto_download=False. Download 'panC_T_all.rds' (or 'CD8.rds' + 'CD4.rds') from "
                f"{CHU_TCELL_PORTAL_URL} (Chu et al., Nature Medicine 2023, DOI: {CHU_TCELL_PUBLICATION_DOI}) "
                f"and stage into '{manual_target}/'."
            )
        match download_chu_tcell_atlas(
            raw_dir=Some(resolved_raw_dir),
            force_download=force_download,
            repo_root=Some(root),
        ):
            case Failure(dl_err):
                return Failure(dl_err)
            case Success(downloaded_rds):
                staged_files = [downloaded_rds]
            case _:
                return Failure("Unexpected Chu_TCell_Atlas download state")

    rds_input = staged_files[0].parent if len(staged_files) > 1 else staged_files[0]
    match preprocess_chu_tcell_atlas(
        rds_path_or_dir=rds_input,
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
            return Failure("Unexpected Chu_TCell_Atlas preprocess state")


# -----------------------------------------------------------------------------
# 3C. Shi et al. (2025, Nature) — `Shi_PanCancer_Atlas`
# -----------------------------------------------------------------------------


def download_shi_pancancer_atlas(
    raw_dir: Maybe[Path] = Nothing,
    force_download: bool = False,
    verify_md5: bool = True,
    repo_root: Maybe[Path] = Nothing,
) -> Result[Path, str]:
    """Download Shi et al. (2025, Nature) `PanCancer_igt_s9_fine_counts.h5ad` from Zenodo `15554080`."""
    if is_cscs_environment():
        return Failure(
            "Internet dataset downloads are strictly disabled on CSCS Eiger. "
            "Download Shi_PanCancer_Atlas on 'olm' first and stage via Globus."
        )

    root = repo_root.unwrap().resolve() if isinstance(repo_root, Some) else find_repo_root()
    target_dir = (
        raw_dir.unwrap().resolve()
        if isinstance(raw_dir, Some)
        else get_raw_dataset_dir(SHI_PANCANCER_DATASET_ID, repo_root=root)
    )
    target_dir.mkdir(parents=True, exist_ok=True)
    dest = target_dir / SHI_PANCANCER_H5AD_FILENAME
    if force_download and dest.exists():
        dest.unlink(missing_ok=True)

    csum: Maybe[ChecksumSpec] = (
        Some(SHI_PANCANCER_CHECKSUM)
        if verify_md5
        else Some(
            ChecksumSpec(
                algorithm="size_only",
                expected_hash="",
                size_bytes=SHI_PANCANCER_CHECKSUM.size_bytes,
            )
        )
    )
    match download_single_file(
        url=SHI_PANCANCER_ZENODO_URL,
        dest_path=dest,
        timeout_sec=180.0,
        max_retries=20,
        resume=True,
        checksum=csum,
    ):
        case Failure(err):
            return Failure(f"Shi_PanCancer_Atlas download failed: {err}")
        case Success(p):
            if verify_md5:
                _ = _write_md5_sidecar(p, SHI_PANCANCER_CHECKSUM.expected_hash)
            return Success(p)
        case _:
            return Failure("Unexpected Shi_PanCancer_Atlas download state")


def preprocess_shi_pancancer_atlas(
    h5ad_path_or_dir: Path,
    output_h5ad: Maybe[Path] = Nothing,
    compartment_filter: ShiCompartmentFilter = ShiCompartmentFilter.T_AND_CAF,
    max_cells_per_substate: Maybe[int] = Nothing,
    normalize_ensembl: bool = True,
    ensembl_release: int = 111,
    drop_unmapped: bool = True,
    batch_size: int = 25_000,
    random_seed: int = 42,
) -> Result[Path, str]:
    """Stream `PanCancer_igt_s9_fine_counts.h5ad` out-of-core via `h5py` CSR row-window slicing.

    Never calls `ad.read_h5ad()` on the 18.2 GB raw file: reads only `obs` and `var` metadata,
    harmonizes 100% of original `.obs` columns via Polars, filters rows by `ShiCompartmentFilter`,
    streams matching CSR rows in `batch_size` windows, projects UMI counts to Ensembl v111, and
    materializes both `.layers["counts"]` and `.X`/`.layers["log1p_norm"]` (`log1p(CP10K)`).
    """
    src_path = (
        h5ad_path_or_dir / SHI_PANCANCER_H5AD_FILENAME
        if h5ad_path_or_dir.is_dir()
        else h5ad_path_or_dir
    )
    if not src_path.is_file():
        return Failure(f"Shi_PanCancer_Atlas raw H5AD not found at {src_path}")

    target_h5ad = (
        output_h5ad.unwrap().resolve()
        if isinstance(output_h5ad, Some)
        else get_preprocessed_h5ad_path(SHI_PANCANCER_DATASET_ID).resolve()
    )
    target_h5ad.parent.mkdir(parents=True, exist_ok=True)
    tmp_h5ad = target_h5ad.with_name(target_h5ad.name + ".tmp")

    try:
        with h5py.File(src_path, "r") as src_h5:
            raw_obs = anndata.io.read_elem(src_h5["obs"])
            raw_var = anndata.io.read_elem(src_h5["var"])

            obs_full = harmonize_shi_pancancer_obs_polars(raw_obs)
            mask = select_shi_compartment_mask(obs_full, compartment_filter=compartment_filter)
            selected_indices = np.flatnonzero(mask).astype(np.int64)

            if len(selected_indices) == 0:
                return Failure(
                    f"Shi_PanCancer_Atlas filter '{compartment_filter.value}' matched 0 cells out of {len(obs_full)}"
                )

            match max_cells_per_substate:
                case Some(cap) if cap > 0:
                    rng = np.random.default_rng(random_seed)
                    sub_labels = obs_full["author_cluster_id"].astype(str).to_numpy()[selected_indices]
                    kept_chunks: list[np.ndarray[Any, Any]] = []
                    for u_lbl in sorted(set(sub_labels.tolist())):
                        lbl_idx = selected_indices[sub_labels == u_lbl]
                        if len(lbl_idx) > cap:
                            chosen = np.sort(rng.choice(lbl_idx, size=cap, replace=False))
                            kept_chunks.append(chosen)
                        else:
                            kept_chunks.append(lbl_idx)
                    selected_indices = np.sort(np.concatenate(kept_chunks))
                case _:
                    pass

            obs_selected = obs_full.iloc[selected_indices].copy()
            n_selected = len(obs_selected)

            # Prefer explicit gene symbol/Ensembl ID columns in `raw_var` if index is integer
            raw_idx_strs = [str(x) for x in raw_var.index.tolist()]
            if raw_idx_strs and raw_idx_strs[0].isdigit():
                sym_col = next(
                    (
                        c
                        for c in ("ensembl_id", "gene_ids", "gene_id", "feature_id", "gene_name", "feature_name", "symbol")
                        if c in raw_var.columns
                    ),
                    "",
                )
                raw_genes = (
                    [str(x) for x in raw_var[sym_col].tolist()]
                    if sym_col
                    else raw_idx_strs
                )
            else:
                raw_genes = raw_idx_strs

            n_raw_genes = len(raw_genes)
            unmapped_genes: list[str] = []
            if normalize_ensembl:
                var_final, keep_mask, proj_m, unmapped_genes = _build_ensembl_projection_operator(
                    raw_var_names=raw_genes,
                    raw_var=raw_var,
                    release=ensembl_release,
                    species="human",
                    ensembl_dir=get_ensembl_dir(),
                    drop_unmapped=drop_unmapped,
                    aggregation="sum",
                )
            else:
                var_final = pd.DataFrame(
                    {"gene_name": raw_genes, "original_id": raw_genes},
                    index=ad.utils.make_index_unique(pd.Index(raw_genes, name="gene_id")),
                )
                keep_mask = np.ones(n_raw_genes, dtype=bool)
                proj_m = sp.eye(n_raw_genes, dtype=np.float32, format="csr")

            n_final_vars = len(var_final)
            total_raw_obs = len(obs_full)
            x_node = src_h5["X"]
            x_enc = str(x_node.attrs.get("encoding-type", "csr_matrix"))

            tmp_h5ad.unlink(missing_ok=True)
            with h5py.File(tmp_h5ad, "w") as dst_h5:
                dst_h5.attrs["encoding-type"] = "anndata"
                dst_h5.attrs["encoding-version"] = "0.1.0"

                write_anndata_dataframe_elem(dst_h5, "obs", obs_selected, default_index_name="cell_id")
                write_anndata_dataframe_elem(dst_h5, "var", var_final, default_index_name="gene_id")

                layers_grp = dst_h5.create_group("layers")
                layers_grp.attrs["encoding-type"] = "dict"
                layers_grp.attrs["encoding-version"] = "0.1.0"

                grp_x = init_compressed_csr_group(dst_h5, "X", n_obs=n_selected, n_vars=n_final_vars)
                grp_log1p = init_compressed_csr_group(layers_grp, "log1p_norm", n_obs=n_selected, n_vars=n_final_vars)
                grp_counts = init_compressed_csr_group(layers_grp, "counts", n_obs=n_selected, n_vars=n_final_vars)

                if isinstance(x_node, h5py.Group) and "csr" in x_enc:
                    ptr_pos = 0
                    n_sel_total = len(selected_indices)
                    for w_start in range(0, total_raw_obs, batch_size):
                        w_end = min(total_raw_obs, w_start + batch_size)
                        next_pos = int(np.searchsorted(selected_indices, w_end, side="left"))
                        if next_pos <= ptr_pos:
                            continue
                        win_sel = selected_indices[ptr_pos:next_pos]
                        ptr_pos = next_pos

                        raw_sub_csr = _slice_h5_csr_selected_rows(
                            csr_grp=x_node,
                            window_start=w_start,
                            window_end=w_end,
                            selected_global_rows=win_sel,
                            n_cols=n_raw_genes,
                        )
                        proj_counts = _project_sparse_chunk_to_ensembl(
                            raw_sub_csr,
                            keep_mask,
                            proj_m,
                            "sum",
                        )
                        proj_log = _normalize_counts_csr_to_log1p_cp10k(proj_counts)
                        append_csr_batch_to_group(grp_counts, proj_counts)
                        append_csr_batch_to_group(grp_x, proj_log)
                        append_csr_batch_to_group(grp_log1p, proj_log)
                        gc.collect()
                else:
                    # Fallback for dense or CSC X (`anndata.io.read_elem` only in bounded slices if Dataset)
                    if isinstance(x_node, h5py.Dataset):
                        for pos_s in range(0, len(selected_indices), batch_size):
                            win_sel = selected_indices[pos_s : pos_s + batch_size]
                            raw_sub_csr = sp.csr_matrix(x_node[win_sel, :], dtype=np.float32)
                            proj_counts = _project_sparse_chunk_to_ensembl(raw_sub_csr, keep_mask, proj_m, "sum")
                            proj_log = _normalize_counts_csr_to_log1p_cp10k(proj_counts)
                            append_csr_batch_to_group(grp_counts, proj_counts)
                            append_csr_batch_to_group(grp_x, proj_log)
                            append_csr_batch_to_group(grp_log1p, proj_log)
                    else:
                        full_mat = sp.csr_matrix(anndata.io.read_elem(x_node))[selected_indices, :]
                        for pos_s in range(0, full_mat.shape[0], batch_size):
                            raw_sub_csr = full_mat[pos_s : pos_s + batch_size, :].tocsr()
                            proj_counts = _project_sparse_chunk_to_ensembl(raw_sub_csr, keep_mask, proj_m, "sum")
                            proj_log = _normalize_counts_csr_to_log1p_cp10k(proj_counts)
                            append_csr_batch_to_group(grp_counts, proj_counts)
                            append_csr_batch_to_group(grp_x, proj_log)
                            append_csr_batch_to_group(grp_log1p, proj_log)

                grp_x.attrs["shape"] = np.array([n_selected, n_final_vars], dtype=np.int64)
                grp_log1p.attrs["shape"] = np.array([n_selected, n_final_vars], dtype=np.int64)
                grp_counts.attrs["shape"] = np.array([n_selected, n_final_vars], dtype=np.int64)

                obsm_grp = dst_h5.create_group("obsm")
                obsm_grp.attrs["encoding-type"] = "dict"
                obsm_grp.attrs["encoding-version"] = "0.1.0"
                if "obsm" in src_h5:
                    for obsm_k in src_h5["obsm"].keys():
                        emb_elem = anndata.io.read_elem(src_h5["obsm"][obsm_k])
                        if isinstance(emb_elem, np.ndarray) and emb_elem.shape[0] == total_raw_obs:
                            anndata.io.write_elem(
                                obsm_grp,
                                obsm_k,
                                np.asarray(emb_elem[selected_indices], dtype=np.float32),
                            )
                        elif isinstance(emb_elem, pd.DataFrame) and len(emb_elem) == total_raw_obs:
                            anndata.io.write_elem(
                                obsm_grp,
                                obsm_k,
                                emb_elem.iloc[selected_indices].to_numpy(dtype=np.float32),
                            )

                comp_counts = obs_selected["lineage_compartment"].value_counts().to_dict()
                uns_payload: dict[str, Any] = {
                    "dataset_id": SHI_PANCANCER_DATASET_ID,
                    "title": "Shi et al. 2025 Cross-Tissue Multicellular Coordination Cancer Atlas",
                    "publication_doi": SHI_PANCANCER_PUBLICATION_DOI,
                    "publication_pmid": SHI_PANCANCER_PUBLICATION_PMID,
                    "zenodo_record": SHI_PANCANCER_ZENODO_RECORD,
                    "zenodo_doi": SHI_PANCANCER_ZENODO_DOI,
                    "source_file": SHI_PANCANCER_H5AD_FILENAME,
                    "compartment_filter": compartment_filter.value,
                    "is_raw_counts": False,
                    "has_raw_counts": True,
                    "expression_type": "log_normalized",
                    "n_unmapped_genes": len(unmapped_genes),
                    "atlas_metadata": {
                        "journal": "Nature",
                        "year": 2025,
                        "total_raw_atlas_cells": int(total_raw_obs),
                        "selected_cells": int(n_selected),
                        "compartment_counts": {str(k): int(v) for k, v in comp_counts.items()},
                        "original_obs_columns": [str(c) for c in raw_obs.columns],
                    },
                }
                anndata.io.write_elem(dst_h5, "uns", uns_payload)

        if target_h5ad.exists():
            target_h5ad.unlink()
        tmp_h5ad.rename(target_h5ad)
        logger.info(
            "Successfully wrote Shi_PanCancer_Atlas H5AD to %s (%d/%d cells x %d genes, filter=%s)",
            target_h5ad,
            n_selected,
            total_raw_obs,
            n_final_vars,
            compartment_filter.value,
        )
        return Success(target_h5ad)
    except Exception as exc:  # noqa: BLE001
        tmp_h5ad.unlink(missing_ok=True)
        return Failure(f"Failed to preprocess Shi_PanCancer_Atlas: {exc}")


def load_shi_pancancer_atlas(
    raw_dir: Maybe[Path] = Nothing,
    output_h5ad: Maybe[Path] = Nothing,
    compartment_filter: ShiCompartmentFilter = ShiCompartmentFilter.T_AND_CAF,
    max_cells_per_substate: Maybe[int] = Nothing,
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
    """Load the Shi et al. (2025, Nature) Pan-Cancer Multicellular Atlas (`Shi_PanCancer_Atlas`)."""
    root = repo_root.unwrap().resolve() if isinstance(repo_root, Some) else find_repo_root()
    resolved_raw_dir = (
        raw_dir.unwrap().resolve()
        if isinstance(raw_dir, Some)
        else get_raw_dataset_dir(SHI_PANCANCER_DATASET_ID, repo_root=root)
    )
    target_h5ad = (
        output_h5ad.unwrap().resolve()
        if isinstance(output_h5ad, Some)
        else get_preprocessed_h5ad_path(SHI_PANCANCER_DATASET_ID, repo_root=root)
    )

    cached_h5ad: Maybe[Path] = Nothing
    if not force_recompute and not force_download:
        cached_h5ad = (
            Some(target_h5ad)
            if target_h5ad.is_file() and target_h5ad.stat().st_size > 0
            else find_dataset_h5ad(SHI_PANCANCER_DATASET_ID, repo_root=root)
        )

    match cached_h5ad:
        case Some(existing_path):
            try:
                adata = _read_selective_h5ad(existing_path, include_all_layers=include_all_layers)
                if compartment_filter != ShiCompartmentFilter.T_AND_CAF and "lineage_compartment" in adata.obs.columns:
                    sub_mask = select_shi_compartment_mask(adata.obs, compartment_filter)
                    adata = adata[sub_mask].copy()
                return Success(
                    _apply_maybe_subset_and_subsample(adata, subset=subset, subsample_n=subsample_n)
                )
            except Exception as exc:  # noqa: BLE001
                return Failure(f"Failed to read cached {existing_path}: {exc}")
        case _:
            pass

    raw_h5ad = resolved_raw_dir / SHI_PANCANCER_H5AD_FILENAME
    if force_download or not raw_h5ad.is_file():
        if not auto_download:
            return Failure(
                f"Shi_PanCancer_Atlas raw file not found at {raw_h5ad} and auto_download=False. "
                f"Download from Zenodo {SHI_PANCANCER_ZENODO_RECORD} or set auto_download=True."
            )
        match download_shi_pancancer_atlas(
            raw_dir=Some(resolved_raw_dir),
            force_download=force_download,
            verify_md5=True,
            repo_root=Some(root),
        ):
            case Failure(dl_err):
                return Failure(dl_err)
            case Success(downloaded_h5ad):
                raw_h5ad = downloaded_h5ad
            case _:
                return Failure("Unexpected Shi_PanCancer_Atlas download state")

    match preprocess_shi_pancancer_atlas(
        h5ad_path_or_dir=raw_h5ad,
        output_h5ad=Some(target_h5ad),
        compartment_filter=compartment_filter,
        max_cells_per_substate=max_cells_per_substate,
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
            return Failure("Unexpected Shi_PanCancer_Atlas preprocess state")


def dispatch_tcell_atlas_load(
    atlas: TCellAtlasChoice,
    raw_dir: Maybe[Path] = Nothing,
    auto_download: bool = True,
    force_download: bool = False,
    subset: Maybe[Mapping[str, str]] = Nothing,
    subsample_n: Maybe[int] = Nothing,
    repo_root: Maybe[Path] = Nothing,
) -> Result[ad.AnnData, str]:
    """Exhaustive `match...case` router across all three T-cell / multicellular atlases."""
    match atlas:
        case TCellAtlasChoice.ZHENG:
            return load_zheng_tcell_atlas(
                raw_dir=raw_dir,
                auto_download=auto_download,
                force_download=force_download,
                subset=subset,
                subsample_n=subsample_n,
                repo_root=repo_root,
            )
        case TCellAtlasChoice.CHU:
            return load_chu_tcell_atlas(
                raw_dir=raw_dir,
                auto_download=auto_download,
                force_download=force_download,
                subset=subset,
                subsample_n=subsample_n,
                repo_root=repo_root,
            )
        case TCellAtlasChoice.SHI:
            return load_shi_pancancer_atlas(
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
    "CHU_TCELL_CANDIDATE_FILENAMES",
    "CHU_TCELL_CLUSTER_SPECS",
    "CHU_TCELL_DATASET_ID",
    "CHU_TCELL_HS_MARKERS",
    "CHU_TCELL_PORTAL_URL",
    "CHU_TCELL_PUBLICATION_DOI",
    "CHU_TCELL_PUBLICATION_PMID",
    "SHI_PANCANCER_CHECKSUM",
    "SHI_PANCANCER_DATASET_ID",
    "SHI_PANCANCER_H5AD_FILENAME",
    "SHI_PANCANCER_PUBLICATION_DOI",
    "SHI_PANCANCER_PUBLICATION_PMID",
    "SHI_PANCANCER_ZENODO_DOI",
    "SHI_PANCANCER_ZENODO_RECORD",
    "SHI_PANCANCER_ZENODO_URL",
    "ShiCompartmentFilter",
    "TCellAtlasChoice",
    "TCellSubstateMetadata",
    "ZHENG_TCELL_CLUSTER_SPECS",
    "ZHENG_TCELL_DATASET_ID",
    "ZHENG_TCELL_EXPR_CHECKSUM",
    "ZHENG_TCELL_EXPR_TAR_FILENAME",
    "ZHENG_TCELL_EXPR_URL",
    "ZHENG_TCELL_GEO_ACCESSION",
    "ZHENG_TCELL_META_CHECKSUM",
    "ZHENG_TCELL_META_TAR_FILENAME",
    "ZHENG_TCELL_META_URL",
    "ZHENG_TCELL_PUBLICATION_DOI",
    "ZHENG_TCELL_PUBLICATION_PMID",
    "ZHENG_TCELL_TCR_CHECKSUM",
    "ZHENG_TCELL_TCR_TAR_FILENAME",
    "ZHENG_TCELL_TCR_URL",
    "ZHENG_TCELL_ZENODO_DOI",
    "ZHENG_TCELL_ZENODO_RECORD",
    "dispatch_tcell_atlas_load",
    "download_chu_tcell_atlas",
    "download_shi_pancancer_atlas",
    "download_zheng_tcell_atlas",
    "harmonize_chu_tcell_obs_polars",
    "harmonize_shi_pancancer_obs_polars",
    "harmonize_zheng_tcell_obs_polars",
    "load_chu_tcell_atlas",
    "load_shi_pancancer_atlas",
    "load_zheng_tcell_atlas",
    "preprocess_chu_tcell_atlas",
    "preprocess_shi_pancancer_atlas",
    "preprocess_zheng_tcell_atlas",
    "select_shi_compartment_mask",
]
