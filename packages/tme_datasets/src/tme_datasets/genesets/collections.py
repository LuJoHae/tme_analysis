"""Pre-defined curated gene set collections for immuno-oncology and TME analysis."""

from __future__ import annotations

from .models import GeneSet, GeneSetCollection

# Ayers et al. Expanded Interferon-gamma / T-cell-inflamed GEP (18 genes)
AYERS_T_CELL_INFLAMED_GEP = GeneSet(
    id="Ayers_Tcell_Inflamed_GEP",
    name="Ayers 18-gene T-cell Inflamed GEP",
    description="Clinical 18-gene signature predictive of response to pembrolizumab",
    genes=(
        "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
        "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
        "PSMB9", "STAT1", "TIGIT",
    ),
)

# TME Major Lineage Markers
TME_MAJOR_MARKERS = {
    "T_NK": GeneSet(
        id="T_NK",
        name="T & NK Cells",
        description="Lineage markers for T and Natural Killer cells",
        genes=("CD3D", "CD3E", "CD8A", "CD4", "NKG7", "NCAM1", "TRAC"),
    ),
    "B_Plasma": GeneSet(
        id="B_Plasma",
        name="B & Plasma Cells",
        description="Lineage markers for B cells and plasma cells",
        genes=("MS4A1", "CD19", "MZB1", "SDC1", "CD79A", "IGHG1"),
    ),
    "Myeloid": GeneSet(
        id="Myeloid",
        name="Myeloid Lineage",
        description="Lineage markers for monocytes, macrophages, and dendritic cells",
        genes=("CD14", "CD68", "LYZ", "CLEC9A", "CD1C", "LILRA4", "FCGR3A"),
    ),
    "Endothelial": GeneSet(
        id="Endothelial",
        name="Endothelial Cells",
        description="Vascular and lymphatic endothelial markers",
        genes=("PECAM1", "VWF", "PLVAP", "CDH5", "KDR"),
    ),
    "Fibroblasts": GeneSet(
        id="Fibroblasts",
        name="Cancer-Associated Fibroblasts",
        description="Stroma and fibroblast activation markers",
        genes=("COL1A1", "COL3A1", "DCN", "LUM", "ACTA2", "FAP"),
    ),
    "Tumor_Epithelial": GeneSet(
        id="Tumor_Epithelial",
        name="Tumor & Epithelial Cells",
        description="Malignant and epithelial markers",
        genes=("EPCAM", "KRT18", "KRT8", "KRT19", "PMEL", "MLANA"),
    ),
}

# TME Subtype Granular Markers
TME_SUBTYPE_MARKERS = {
    "CD8_Cytotoxic_T": GeneSet(
        id="CD8_Cytotoxic_T",
        name="CD8 Cytotoxic T Cells",
        description="Effector cytotoxic T cell markers",
        genes=("CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "IFNG"),
    ),
    "CD8_Exhausted_T": GeneSet(
        id="CD8_Exhausted_T",
        name="CD8 Exhausted T Cells",
        description="T-cell exhaustion and checkpoint inhibitory receptors",
        genes=("PDCD1", "LAG3", "HAVCR2", "TOX", "TIGIT", "CTLA4", "ENTPD1"),
    ),
    "Regulatory_T": GeneSet(
        id="Regulatory_T",
        name="Regulatory T Cells (Treg)",
        description="Immunosuppressive regulatory T cell markers",
        genes=("FOXP3", "IL2RA", "BATF", "IKZF2", "CTLA4"),
    ),
    "M1_Macrophage": GeneSet(
        id="M1_Macrophage",
        name="M1 Pro-inflammatory Macrophages",
        description="Pro-inflammatory anti-tumor macrophage markers",
        genes=("CD68", "NOS2", "TNF", "IL1B", "CXCL10", "IL6"),
    ),
    "M2_Macrophage": GeneSet(
        id="M2_Macrophage",
        name="M2 Immunosuppressive Macrophages",
        description="Pro-tumor wound-healing macrophage markers",
        genes=("CD68", "CD163", "MRC1", "ARG1", "TGFB1", "MSR1"),
    ),
}

# Bagaev Multifunctional Portrait (MFP) core signatures subset
BAGAEV_CORE_SIGNATURES = {
    "T_cells": GeneSet(id="T_cells", name="T Cells", genes=("CD3D", "CD3E", "CD3G", "CD2")),
    "CD8_T_cells": GeneSet(id="CD8_T_cells", name="CD8 T Cells", genes=("CD8A", "CD8B")),
    "Cytotoxic_cells": GeneSet(
        id="Cytotoxic_cells", name="Cytotoxic Cells", genes=("GZMA", "GZMB", "PRF1", "NKG7")
    ),
    "Checkpoint_inhibition": GeneSet(
        id="Checkpoint_inhibition",
        name="Checkpoint Inhibition",
        genes=("PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT"),
    ),
    "Endothelium": GeneSet(
        id="Endothelium", name="Endothelium", genes=("PECAM1", "VWF", "CDH5", "KDR")
    ),
    "Fibroblasts": GeneSet(
        id="Fibroblasts", name="Fibroblasts", genes=("COL1A1", "COL1A2", "ACTA2", "DCN")
    ),
    "Angiogenesis": GeneSet(
        id="Angiogenesis", name="Angiogenesis", genes=("VEGFA", "KDR", "FLT1", "ANGPT2")
    ),
    "MHC_class_I": GeneSet(
        id="MHC_class_I", name="MHC Class I", genes=("HLA-A", "HLA-B", "HLA-C", "B2M")
    ),
    "MHC_class_II": GeneSet(
        id="MHC_class_II",
        name="MHC Class II",
        genes=("HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1"),
    ),
}


def get_tme_major_lineage_collection() -> GeneSetCollection:
    """Return collection of primary TME lineage marker gene sets."""
    return GeneSetCollection(
        id="tme_major_lineages",
        name="TME Major Lineages",
        description="Major cellular lineages across the tumor microenvironment",
        gene_sets=TME_MAJOR_MARKERS,
    )


def get_tme_subtype_collection() -> GeneSetCollection:
    """Return collection of granular TME cellular subtype signatures."""
    return GeneSetCollection(
        id="tme_subtypes",
        name="TME Cell Subtypes",
        description="Granular immune and stromal subtype markers",
        gene_sets=TME_SUBTYPE_MARKERS,
    )


def get_bagaev_core_collection() -> GeneSetCollection:
    """Return collection of core Bagaev MFP functional TME signatures."""
    return GeneSetCollection(
        id="bagaev_core",
        name="Bagaev MFP Core Signatures",
        description="Functional gene expression signatures defining the 4 TME subtypes",
        gene_sets=BAGAEV_CORE_SIGNATURES,
    )
