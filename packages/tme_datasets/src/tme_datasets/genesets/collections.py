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


# Canonical ordering of all 29 Bagaev et al. (2021, Cancer Cell) Fges signatures
# matching BostonGene/MFP `gene_signatures_order.tsv`
BAGAEV_MFP_SIGNATURE_ORDER: tuple[str, ...] = (
    "Angiogenesis",
    "Endothelium",
    "CAF",
    "Matrix",
    "Matrix_remodeling",
    "Protumor_cytokines",
    "Neutrophil_signature",
    "Granulocyte_traffic",
    "Macrophages",
    "Macrophage_DC_traffic",
    "MDSC_traffic",
    "MDSC",
    "Th2_signature",
    "T_reg_traffic",
    "Treg",
    "M1_signatures",
    "MHCII",
    "Antitumor_cytokines",
    "Coactivation_molecules",
    "B_cells",
    "NK_cells",
    "Checkpoint_inhibition",
    "Effector_cells",
    "T_cells",
    "Th1_signature",
    "T_cell_traffic",
    "MHCI",
    "EMT_signature",
    "Proliferation_rate",
)

BAGAEV_MFP_FUNCTIONAL_GROUPS: dict[str, tuple[str, ...]] = {
    "Angiogenesis_Fibroblasts": (
        "Angiogenesis",
        "Endothelium",
        "CAF",
        "Matrix",
        "Matrix_remodeling",
        "Protumor_cytokines",
    ),
    "Protumor_Immune": (
        "Neutrophil_signature",
        "Granulocyte_traffic",
        "Macrophages",
        "Macrophage_DC_traffic",
        "MDSC_traffic",
        "MDSC",
        "Th2_signature",
        "T_reg_traffic",
        "Treg",
    ),
    "Antitumor_Immune": (
        "M1_signatures",
        "MHCII",
        "Antitumor_cytokines",
        "Coactivation_molecules",
        "B_cells",
        "NK_cells",
        "Checkpoint_inhibition",
        "Effector_cells",
        "T_cells",
        "Th1_signature",
        "T_cell_traffic",
        "MHCI",
    ),
    "Malignant_Properties": (
        "EMT_signature",
        "Proliferation_rate",
    ),
}

# Complete 29 Functional Gene Expression Signatures (Fges) from Bagaev et al. (2021)
BAGAEV_MFP_29_SIGNATURES: dict[str, GeneSet] = {
    "MHCI": GeneSet(
        id="MHCI",
        name="MHCI",
        genes=("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5", "TAPBP"),
    ),
    "MHCII": GeneSet(
        id="MHCII",
        name="MHCII",
        genes=(
            "HLA-DRA", "HLA-DRB1", "HLA-DMA", "HLA-DPA1", "HLA-DPB1",
            "HLA-DMB", "HLA-DQB1", "HLA-DQA1", "CIITA",
        ),
    ),
    "Coactivation_molecules": GeneSet(
        id="Coactivation_molecules",
        name="Co-stimulatory molecules",
        genes=(
            "CD28", "CD40", "TNFRSF4", "ICOS", "TNFRSF9", "CD27", "CD80",
            "CD86", "CD40LG", "CD83", "TNFSF4", "ICOSLG", "TNFSF9", "CD70",
        ),
    ),
    "Effector_cells": GeneSet(
        id="Effector_cells",
        name="Effector cells",
        genes=(
            "IFNG", "GZMA", "GZMB", "PRF1", "GZMK", "ZAP70", "GNLY",
            "FASLG", "TBX21", "EOMES", "CD8A", "CD8B",
        ),
    ),
    "T_cell_traffic": GeneSet(
        id="T_cell_traffic",
        name="Effector cell traffic",
        genes=("CXCL9", "CXCL10", "CXCL11", "CX3CL1", "CCL3", "CCL4", "CX3CR1", "CCL5", "CXCR3"),
    ),
    "NK_cells": GeneSet(
        id="NK_cells",
        name="NK cells",
        genes=(
            "NKG7", "CD160", "CD244", "NCR1", "KLRC2", "KLRK1", "CD226", "GZMH",
            "GNLY", "IFNG", "KIR2DL4", "EOMES", "GZMB", "FGFBP2", "KLRF1", "SH2D1B", "NCR3",
        ),
    ),
    "T_cells": GeneSet(
        id="T_cells",
        name="T cells",
        genes=("TBX21", "ITK", "CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "CD28", "CD5", "TRAT1"),
    ),
    "B_cells": GeneSet(
        id="B_cells",
        name="B cells",
        genes=(
            "CD19", "MS4A1", "TNFRSF13C", "CR2", "TNFRSF17", "TNFRSF13B",
            "CD22", "CD79A", "CD79B", "BLK", "FCRL5", "PAX5", "STAP1",
        ),
    ),
    "M1_signatures": GeneSet(
        id="M1_signatures",
        name="M1 signature",
        genes=("NOS2", "TNF", "IL1B", "SOCS3", "CMKLR1", "IRF5", "IL12A", "IL12B", "IL23A"),
    ),
    "Th1_signature": GeneSet(
        id="Th1_signature",
        name="Th1 signature",
        genes=("IFNG", "IL2", "CD40LG", "IL21", "TBX21", "STAT4", "IL12RB2"),
    ),
    "Antitumor_cytokines": GeneSet(
        id="Antitumor_cytokines",
        name="Antitumor cytokines",
        genes=("TNF", "IFNB1", "IFNA2", "CCL3", "TNFSF10", "IL21"),
    ),
    "Checkpoint_inhibition": GeneSet(
        id="Checkpoint_inhibition",
        name="Checkpoint molecules",
        genes=("PDCD1", "CD274", "CTLA4", "LAG3", "PDCD1LG2", "BTLA", "HAVCR2", "TIGIT", "VSIR"),
    ),
    "Treg": GeneSet(
        id="Treg",
        name="Treg",
        genes=("FOXP3", "CTLA4", "IL10", "TNFRSF18", "CCR8", "IKZF4", "IKZF2"),
    ),
    "T_reg_traffic": GeneSet(
        id="T_reg_traffic",
        name="Treg and Th2 traffic",
        genes=("CCL17", "CCL22", "CCL1", "CCL28", "CCR4", "CCR8", "CCR10"),
    ),
    "Neutrophil_signature": GeneSet(
        id="Neutrophil_signature",
        name="Neutrophil signature",
        genes=("MPO", "ELANE", "PRTN3", "CTSG", "CXCR1", "CXCR2", "FCGR3B", "CD177", "FFAR2", "PGLYRP1"),
    ),
    "Granulocyte_traffic": GeneSet(
        id="Granulocyte_traffic",
        name="Granulocyte traffic",
        genes=("CXCL8", "CXCL1", "CXCL2", "CXCL5", "CCL11", "KITLG", "CXCR1", "CXCR2", "CCR3"),
    ),
    "MDSC": GeneSet(
        id="MDSC",
        name="Immune Suppression by Myeloid Cells",
        genes=("IDO1", "ARG1", "IL10", "CYBB", "PTGS2", "IL4I1", "IL6"),
    ),
    "MDSC_traffic": GeneSet(
        id="MDSC_traffic",
        name="Myeloid cells traffic",
        genes=(
            "CSF2", "CSF3", "CXCL12", "CCL26", "IL6", "CXCL8", "CXCL5",
            "CSF1R", "CSF2RA", "CSF3R", "CXCR4", "IL6R", "CXCR2", "CCL15", "CSF1",
        ),
    ),
    "Macrophages": GeneSet(
        id="Macrophages",
        name="Tumor-associated Macrophages",
        genes=("IL10", "MRC1", "MSR1", "CD163", "CSF1R", "IL4I1", "SIGLEC1", "CD68"),
    ),
    "Macrophage_DC_traffic": GeneSet(
        id="Macrophage_DC_traffic",
        name="Macrophage and DC traffic",
        genes=("CCL2", "CCL7", "CCL8", "XCL1", "CCR2", "XCR1", "CSF1R", "CSF1"),
    ),
    "Th2_signature": GeneSet(
        id="Th2_signature",
        name="Th2 signature",
        genes=("IL4", "IL5", "IL13", "IL10", "CCR4"),
    ),
    "Protumor_cytokines": GeneSet(
        id="Protumor_cytokines",
        name="Protumor cytokines",
        genes=("IL10", "TGFB1", "TGFB2", "TGFB3", "IL22", "MIF", "IL6"),
    ),
    "CAF": GeneSet(
        id="CAF",
        name="Fibroblasts",
        genes=(
            "COL1A1", "COL1A2", "COL5A1", "ACTA2", "FGF2", "FAP", "LRP1",
            "CD248", "COL6A1", "COL6A2", "COL6A3", "CXCL12", "FBLN1", "LUM",
            "MFAP5", "MMP3", "MMP2", "PDGFRB", "PDGFRA",
        ),
    ),
    "Matrix": GeneSet(
        id="Matrix",
        name="Matrix",
        genes=(
            "FN1", "COL1A1", "COL1A2", "COL4A1", "COL3A1", "VTN", "LGALS7",
            "LGALS9", "LAMA3", "LAMB3", "LAMC2", "TNC", "ELN", "COL5A1", "COL11A1",
        ),
    ),
    "Matrix_remodeling": GeneSet(
        id="Matrix_remodeling",
        name="Matrix remodeling",
        genes=(
            "CA9", "MMP9", "MMP2", "MMP1", "MMP3", "MMP12",
            "MMP7", "MMP11", "PLOD2", "ADAMTS4", "ADAMTS5", "LOX",
        ),
    ),
    "Angiogenesis": GeneSet(
        id="Angiogenesis",
        name="Angiogenesis",
        genes=(
            "VEGFA", "VEGFB", "VEGFC", "PDGFC", "CXCL8", "CXCR2", "FLT1",
            "PGF", "CXCL5", "KDR", "ANGPT1", "ANGPT2", "TEK", "VWF", "CDH5",
        ),
    ),
    "Endothelium": GeneSet(
        id="Endothelium",
        name="Endothelium",
        genes=("NOS3", "KDR", "FLT1", "VCAM1", "VWF", "CDH5", "MMRN1", "ENG", "CLEC14A", "MMRN2"),
    ),
    "Proliferation_rate": GeneSet(
        id="Proliferation_rate",
        name="Tumor proliferation rate",
        genes=(
            "MKI67", "ESCO2", "CETN3", "CDK2", "CCND1", "CCNE1", "AURKA",
            "AURKB", "E2F1", "MYBL2", "BUB1", "PLK1", "CCNB1", "MCM2", "MCM6",
        ),
    ),
    "EMT_signature": GeneSet(
        id="EMT_signature",
        name="EMT signature",
        genes=("SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2", "CDH2"),
    ),
}


def get_bagaev_mfp_collection() -> GeneSetCollection:
    """Return the complete 29-signature Bagaev et al. (2021) Molecular Functional Portrait collection."""
    ordered_sets = {
        sig_id: BAGAEV_MFP_29_SIGNATURES[sig_id]
        for sig_id in BAGAEV_MFP_SIGNATURE_ORDER
    }
    return GeneSetCollection(
        id="bagaev_mfp_29",
        name="Bagaev MFP 29 Functional Gene Expression Signatures",
        description="Complete 29 functional gene expression signatures defining conserved pan-cancer TME subtypes (IE/F, IE, F, D)",
        gene_sets=ordered_sets,
    )

# Immune Checkpoint Hallmark Markers (13 genes)
IMMUNE_CHECKPOINT_GENES: tuple[str, ...] = (
    "CD274",   # PD-L1
    "PDCD1",   # PD-1
    "CTLA4",   # CTLA-4
    "LAG3",    # LAG-3
    "TIGIT",   # TIGIT
    "HAVCR2",  # TIM-3
    "IFNG",    # Interferon gamma
    "CXCL9",   # CXCL9
    "CXCL10",  # CXCL10
    "IDO1",    # IDO1
    "FOXP3",   # FOXP3
    "GZMB",    # Granzyme B
    "PRF1",    # Perforin 1
)

# Comprehensive Curated Immunotherapy Gene Panel (68 genes)
IMMUNOTHERAPY_GENE_PANEL: tuple[str, ...] = (
    # Checkpoints & Co-inhibitory Receptors / Ligands
    "CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "BTLA",
    "CD160", "KLRG1", "LAIR1", "CD276", "VTCN1", "LGALS9", "PVR", "CD80", "CD86",
    "HHLA2", "VSIR",
    # Co-stimulatory Receptors & Ligands
    "CD28", "ICOS", "ICOSLG", "CD27", "CD70", "CD40", "CD40LG", "TNFRSF4",
    "TNFSF4", "TNFRSF9", "TNFSF9", "TNFRSF18", "TNFSF18", "TNFRSF14",
    # Cytolytic Effectors & T Cell Activation
    "IFNG", "GZMA", "GZMB", "GZMK", "GZMH", "PRF1", "NKG7", "GNLY", "FASLG",
    "TNF", "IL2", "STAT1",
    # Chemokines & Homing Receptors
    "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CCL5", "CCL4", "CCL3", "CXCR3",
    "CXCR6", "CCR5", "CCR7",
    # Antigen Presentation Machinery & MHC
    "B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9",
    "NLRC5", "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-E", "CIITA",
    # Immunosuppression, Metabolic Enzymes, & Stroma
    "IDO1", "ARG1", "NOS2", "ENTPD1", "NT5E", "PTGS2", "TGFB1", "IL10", "VEGFA",
    "ANGPT2", "IL6", "AXL", "FAP", "COL1A1", "COL3A1",
    # Exhaustion & Differentiation Regulators
    "FOXP3", "TOX", "PRDM1", "EOMES", "TBX21", "BATF", "IKZF2", "TCF7",
    # Lineage Co-receptors
    "CD8A", "CD8B", "CD4", "CD3D", "CD3E", "CD3G",
)
