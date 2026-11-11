"""Immutable domain models, exhaustive sum types, Protocols, and Coulton (2024) clinical metadata.

Strictly follows functional Python standards (.agents/rules/code-style-guide.md):
- Frozen Pydantic models (`ConfigDict(frozen=True)`), zero `None` / `Optional` (`returns.maybe.Maybe`)
- Exhaustive `StrEnum` sum types with `typing.assert_never`
- Side-effecting boundaries typed via `typing.Protocol` typeclasses
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from enum import StrEnum
from pathlib import Path
from typing import Final, Protocol, TypeAlias, assert_never

import altair as alt
import numpy as np
import numpy.typing as npt
import polars as pl
from pydantic import BaseModel, ConfigDict, Field
from returns.maybe import Maybe, Nothing, Some
from returns.result import Result

FloatMatrix: TypeAlias = npt.NDArray[np.float64]
FloatVector: TypeAlias = npt.NDArray[np.float64]
IntVector: TypeAlias = npt.NDArray[np.int64]
BoolVector: TypeAlias = npt.NDArray[np.bool_]
AltairChart: TypeAlias = alt.Chart | alt.LayerChart | alt.HConcatChart | alt.VConcatChart


# =============================================================================
# 1. Exhaustive Sum Types (StrEnum)
# =============================================================================


class LineageCompartment(StrEnum):
    """Canonical broad cellular lineage compartments across Consensus-Meta42 (M01-M42)."""

    T_NK = "T_NK"
    B_CELL = "B_CELL"
    PLASMA = "PLASMA"
    MYELOID_MAST = "MYELOID_MAST"
    STROMAL_CAF = "STROMAL_CAF"
    ENDOTHELIAL = "ENDOTHELIAL"
    EPITHELIAL_TUMOR = "EPITHELIAL_TUMOR"
    OTHER_UNSPECIFIED = "OTHER_UNSPECIFIED"


class Meta42MacrophageTier(StrEnum):
    """Hierarchical macrophage/myeloid proximity tier for Consensus-Meta42 states."""

    CORE_TAM = "Tier 1: Core TAM States (M14, M32, M42)"
    MONO_MDSC_CONTINUUM = "Tier 2: Monocyte / MDSC Continuum (M10, M19, M28, M35)"
    OTHER_MYELOID = "Tier 3: Other Myeloid / DC / Mast (M09, M31, M38, M39, M41)"
    NON_MYELOID_CONTEXT = "Tier 4: Non-Myeloid TME Context (30 States)"


class CoultonClusterRole(StrEnum):
    """Biological ontology role of each cluster in Coulton et al. (2024) (`0_AlvMac`..`23_NA`)."""

    TRUE_TAM_SUBSTATE = "True Macrophage / TAM Substate"
    CLASSICAL_MONOCYTE = "Classical Monocyte Precursor"
    T_MYELOID_DOUBLET_CONTROL = "T/Myeloid Doublet Control"
    UNANNOTATED_ARTIFACT = "Unannotated Residual (n=2)"


class TissueStratum(StrEnum):
    """Tissue sampling scope within Coulton_TAM_Atlas."""

    TUMOR_ONLY = "tumor_only"
    ALL_TISSUES = "all_tissues"


class AssaySource(StrEnum):
    """Processed Seurat assay slot in Coulton_TAM_Atlas."""

    RNA_LOG1P = "RNA"
    SCT_DATA = "SCT"


class IciOutcomeDirection(StrEnum):
    """Clinical ICI checkpoint therapy association direction."""

    RESPONDER_ENRICHED = "Responder-Enriched"
    NON_RESPONDER_ENRICHED = "Non-Responder-Enriched"
    NEUTRAL_OR_NS = "Neutral / Non-Significant"


class OutcomeConcordanceClass(StrEnum):
    """Joint ICI outcome agreement between a Meta42 state and its closest Coulton TAM substate."""

    CONCORDANT_PRO_RESPONSE = "Concordant Pro-Response (+/+)"
    CONCORDANT_PRO_RESISTANCE = "Concordant Pro-Resistance (-/-)"
    POLARIZATION_SPLIT = "Polarization-Split (+/- or -/+)"
    UNDETERMINED_NEUTRAL = "Neutral / Sub-Threshold"


class SignificanceTier(StrEnum):
    """Standard four-tier statistical significance scale (*, **, ***, ****)."""

    NS = "ns"
    ONE_STAR = "*"
    TWO_STARS = "**"
    THREE_STARS = "***"
    FOUR_STARS = "****"


def tissue_stratum_display_label(stratum: TissueStratum) -> str:
    """Human-readable label for each Coulton tissue stratum."""
    match stratum:
        case TissueStratum.TUMOR_ONLY:
            return "Tumor Tissue Only (n = 279,104 cells)"
        case TissueStratum.ALL_TISSUES:
            return "All Atlas Tissues (n = 363,315 cells)"
        case _ as unreachable:
            assert_never(unreachable)


def assay_source_display_label(assay: AssaySource) -> str:
    """Human-readable label for each Coulton assay slot."""
    match assay:
        case AssaySource.RNA_LOG1P:
            return "RNA@data (log1p_norm, linear expm1 centroid)"
        case AssaySource.SCT_DATA:
            return "SCT@data (SCTransform log-corrected centroid)"
        case _ as unreachable:
            assert_never(unreachable)


# =============================================================================
# 2. Frozen Domain Models
# =============================================================================


class ComparisonConfig(BaseModel):
    """Immutable configuration for the Meta42 <-> Coulton_TAM_Atlas comparison pipeline."""

    model_config = ConfigDict(frozen=True)

    phi_meta_path: Path = Field(
        default=Path(
            "output/manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark/"
            "unified_reference_phi_meta.parquet"
        )
    )
    concordance_3b_path: Path = Field(
        default=Path(
            "output/manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark/"
            "unified_all_states_detailed_concordance.parquet"
        )
    )
    hub_summary_path: Path = Field(
        default=Path(
            "output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark/"
            "consensus_metastate_all42_hub_summary.parquet"
        )
    )
    state_effects_path: Path = Field(
        default=Path(
            "output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark/"
            "milopy_iatlas_state_effects.parquet"
        )
    )
    top_degs_path: Path = Field(
        default=Path(
            "output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark/"
            "consensus_metastate_all42_top_degs.parquet"
        )
    )
    output_dir: Path = Field(default=Path("output/results/coulton_tam_comparison"))
    batch_size: int = Field(default=25_000, ge=1_000)
    top_k_markers: int = Field(default=25, ge=5, le=100)
    bulk_signature_genes: int = Field(default=15, ge=5, le=50)
    include_sct_sensitivity: bool = Field(default=True)
    reuse_cached_centroids: bool = Field(default=True)
    render_plots: bool = Field(default=True)


class CoultonLiteratureSpec(BaseModel):
    """Published phenotype & CPI1000+ (N=1,446) ICI outcome metadata from Coulton et al. (2024)."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    cluster_id: int
    short_label: str
    role: CoultonClusterRole
    functional_program: str
    canonical_markers: tuple[str, ...]
    is_ucell_gold_standard: bool
    cpi1000_direction: IciOutcomeDirection
    cpi1000_fgsea_qval: Maybe[float]
    mana_and_spatial_note: str


class CoultonCentroidAccumulator(BaseModel):
    """Immutable running accumulator for out-of-core sparse CSR cluster centroid reduction."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    assay: AssaySource
    cluster_labels: tuple[str, ...]
    gene_symbols: tuple[str, ...]
    all_cell_counts: IntVector
    tumor_cell_counts: IntVector
    all_linear_sum: FloatMatrix
    tumor_linear_sum: FloatMatrix
    all_log_sum: FloatMatrix
    tumor_log_sum: FloatMatrix
    all_detect_counts: IntVector | FloatMatrix
    tumor_detect_counts: IntVector | FloatMatrix


class CoultonCentroidBundle(BaseModel):
    """Finalized cluster expression profiles for one assay across Tumor-Only and All-Tissues."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    assay: AssaySource
    cluster_labels: tuple[str, ...]
    gene_symbols: tuple[str, ...]
    all_cell_counts: IntVector
    tumor_cell_counts: IntVector
    all_mean_linear: FloatMatrix
    tumor_mean_linear: FloatMatrix
    all_mean_log: FloatMatrix
    tumor_mean_log: FloatMatrix
    all_pct_expressing: FloatMatrix
    tumor_pct_expressing: FloatMatrix


class Meta42ReferenceBundle(BaseModel):
    """Immutable bundle of Consensus-Meta42 profiles and pan-cohort ICI outcome tables."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    phi_meta_df: pl.DataFrame
    concordance_3b_df: pl.DataFrame
    hub_summary_df: pl.DataFrame
    state_effects_df: pl.DataFrame
    top_degs_df: pl.DataFrame


class BulkIcbCohortBundle(BaseModel):
    """Immutable standardized bulk RNA-seq clinical ICB cohort for empirical TAM scoring."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    cohort_id: str
    cancer_type: str
    short_name: str
    n_samples: int
    n_responders: int
    response_binary: FloatVector
    genes: tuple[str, ...]
    gene_to_idx: Mapping[str, int]
    log2_z_mat: FloatMatrix


class ComparisonOutputs(BaseModel):
    """Complete pure deliverables of the Meta42 <-> Coulton_TAM_Atlas comparison pipeline."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    coulton_centroids_df: pl.DataFrame
    coulton_top_markers_df: pl.DataFrame
    full_similarity_df: pl.DataFrame
    nearest_substates_df: pl.DataFrame
    sensitivity_summary_df: pl.DataFrame
    coulton_ici_outcomes_df: pl.DataFrame
    meta42_indication_effects_df: pl.DataFrame
    markdown_report: str


class DendrogramBracket(BaseModel):
    """Four-point orthogonal U-bracket `(u_left, h_left) -> (u_left, h_merge) -> (u_right, h_merge) -> (u_right, h_right)`."""

    model_config = ConfigDict(frozen=True)

    bracket_id: int
    u_left: float
    h_left: float
    u_right: float
    h_right: float
    h_merge: float


class BiclusterResult(BaseModel):
    """Immutable hierarchical Ward biclustering order and normalized marginal dendrogram brackets."""

    model_config = ConfigDict(frozen=True)

    ordered_rows: tuple[str, ...]
    ordered_cols: tuple[str, ...]
    row_brackets: tuple[DendrogramBracket, ...]
    col_brackets: tuple[DendrogramBracket, ...]


# =============================================================================
# 3. Typeclass Protocols (Side-Effecting IO Boundaries)
# =============================================================================


class CoultonAtlasStreamer(Protocol):
    """Streams cluster centroids from Coulton_TAM_Atlas strictly via tme_datasets."""

    def __call__(
        self,
        assay: AssaySource,
        batch_size: int,
    ) -> Result[CoultonCentroidBundle, str]: ...


class Meta42ArtifactLoader(Protocol):
    """Loads Consensus-Meta42 reference profiles and ICI effect parquets."""

    def __call__(self, config: ComparisonConfig, /) -> Result[Meta42ReferenceBundle, str]: ...


class BulkIcbSuiteLoader(Protocol):
    """Loads clinical bulk RNA-seq ICB cohorts strictly via tme_datasets.load_dataset."""

    def __call__(self) -> Result[tuple[BulkIcbCohortBundle, ...], str]: ...


class VectorChartExporter(Protocol):
    """Exports a declarative Altair chart to vector SVG and 300 DPI PNG."""

    def __call__(self, chart: AltairChart, base_path: Path, /) -> Result[tuple[Path, Path], str]: ...


# =============================================================================
# 4. Canonical Metadata Maps (`Meta42` & Coulton et al. 2024 `PMC11226649`)
# =============================================================================


META42_BIO_LABELS: Final[Mapping[str, str]] = {
    "M01": "M01: TRBV18+GCNT4+ Mem T",
    "M02": "M02: GATA4+CADM2+ Stroma",
    "M03": "M03: CD200+ELMO1+ Tfh T",
    "M04": "M04: GNLY+PRF1+ Eff. T",
    "M05": "M05: CD96+TRAF3+ Mem. T",
    "M06": "M06: CX3CR1+FGFBP2+ NK/T",
    "M07": "M07: ANK1+TRIM58+ Eryth.",
    "M08": "M08: CCR4+SOCS1+ Th2/Treg",
    "M09": "M09: SIGLEC8+CR1+ Granul.",
    "M10": "M10: FCN1+CD1A+ Mono/cDC",
    "M11": "M11: AQP4+CTSE+ Ductal",
    "M12": "M12: CA9+SLC2A1+ Hypox.",
    "M13": "M13: CXCL13+PDCD1+ Tex",
    "M14": "M14: CXCL9+STAT1+ M1-TAM",
    "M15": "M15: TCF7+IL7R+ Stem CD8",
    "M16": "M16: VEGFA+KDR+ Endo.",
    "M17": "M17: ISLR+COL5A2+ myCAF",
    "M18": "M18: S100B+TYR+ Melanoc.",
    "M19": "M19: S100A8+CXCL8+ MDSC",
    "M20": "M20: IL7R+CD40LG+ Trm T",
    "M21": "M21: FOXP3+CTLA4+ Treg",
    "M22": "M22: MKI67+TOP2A+ Prolif",
    "M23": "M23: FAP+TGFB1+ Mat. CAF",
    "M24": "M24: CXCL13+NMB+ Tfh/Trm",
    "M25": "M25: CTLA4+TIGIT+ Tex",
    "M26": "M26: VIM+ZEB1+ EMT Tumor",
    "M27": "M27: ABCC8+PPY+ Neuroend",
    "M28": "M28: CXCL8+ITGAM+ MDSC",
    "M29": "M29: JCHAIN+MZB1+ Plasma",
    "M30": "M30: CR2+PAX5+ TLS B",
    "M31": "M31: CLEC4C+LILRA4+ pDC",
    "M32": "M32: MRC1+CD163+ M2-TAM",
    "M33": "M33: CCL5+GNLY+ CD8 Teff",
    "M34": "M34: GZMK+NKG7+ Tem CD8+",
    "M35": "M35: MPO+CD79A+ Myeloid",
    "M36": "M36: FCGR3A+GNLY+ NK",
    "M37": "M37: MS4A1+CD79A+ TLS B",
    "M38": "M38: FCGR3B+FPR2+ TAN",
    "M39": "M39: MS4A2+TPSAB1+ Mast",
    "M40": "M40: FOXP3+CCR8+ Treg",
    "M41": "M41: LAMP3+CCR7+ mregDC",
    "M42": "M42: SULT1B1+RGS2+ TAM",
}

META42_COMPARTMENT_MAP: Final[Mapping[str, LineageCompartment]] = {
    "M01": LineageCompartment.T_NK,
    "M02": LineageCompartment.STROMAL_CAF,
    "M03": LineageCompartment.T_NK,
    "M04": LineageCompartment.T_NK,
    "M05": LineageCompartment.T_NK,
    "M06": LineageCompartment.T_NK,
    "M07": LineageCompartment.OTHER_UNSPECIFIED,
    "M08": LineageCompartment.T_NK,
    "M09": LineageCompartment.MYELOID_MAST,
    "M10": LineageCompartment.MYELOID_MAST,
    "M11": LineageCompartment.EPITHELIAL_TUMOR,
    "M12": LineageCompartment.EPITHELIAL_TUMOR,
    "M13": LineageCompartment.T_NK,
    "M14": LineageCompartment.MYELOID_MAST,
    "M15": LineageCompartment.T_NK,
    "M16": LineageCompartment.ENDOTHELIAL,
    "M17": LineageCompartment.STROMAL_CAF,
    "M18": LineageCompartment.EPITHELIAL_TUMOR,
    "M19": LineageCompartment.MYELOID_MAST,
    "M20": LineageCompartment.T_NK,
    "M21": LineageCompartment.T_NK,
    "M22": LineageCompartment.EPITHELIAL_TUMOR,
    "M23": LineageCompartment.STROMAL_CAF,
    "M24": LineageCompartment.T_NK,
    "M25": LineageCompartment.T_NK,
    "M26": LineageCompartment.EPITHELIAL_TUMOR,
    "M27": LineageCompartment.EPITHELIAL_TUMOR,
    "M28": LineageCompartment.MYELOID_MAST,
    "M29": LineageCompartment.PLASMA,
    "M30": LineageCompartment.B_CELL,
    "M31": LineageCompartment.MYELOID_MAST,
    "M32": LineageCompartment.MYELOID_MAST,
    "M33": LineageCompartment.T_NK,
    "M34": LineageCompartment.T_NK,
    "M35": LineageCompartment.MYELOID_MAST,
    "M36": LineageCompartment.T_NK,
    "M37": LineageCompartment.B_CELL,
    "M38": LineageCompartment.MYELOID_MAST,
    "M39": LineageCompartment.MYELOID_MAST,
    "M40": LineageCompartment.T_NK,
    "M41": LineageCompartment.MYELOID_MAST,
    "M42": LineageCompartment.MYELOID_MAST,
}

CORE_TAM_META42_IDS: Final[tuple[str, ...]] = ("M14", "M32", "M42")
MONO_MDSC_META42_IDS: Final[tuple[str, ...]] = ("M10", "M19", "M28", "M35")
OTHER_MYELOID_META42_IDS: Final[tuple[str, ...]] = ("M09", "M31", "M38", "M39", "M41")
MACROPHAGE_CONTINUUM_META42_IDS: Final[tuple[str, ...]] = (
    CORE_TAM_META42_IDS + MONO_MDSC_META42_IDS
)
ALL_MYELOID_META42_IDS: Final[tuple[str, ...]] = (
    CORE_TAM_META42_IDS + MONO_MDSC_META42_IDS + OTHER_MYELOID_META42_IDS
)

PAN_MACROPHAGE_CORE_GENES: Final[tuple[str, ...]] = (
    "CD68",
    "CD14",
    "CSF1R",
    "TYROBP",
    "FCER1G",
    "AIF1",
    "C1QA",
    "C1QB",
    "C1QC",
    "MSR1",
)

BULK_ICB_COHORT_SPECS: Final[tuple[tuple[str, str, str], ...]] = (
    ("Hugo-iAtlas", "Melanoma", "Hugo (Melanoma)"),
    ("Riaz-iAtlas", "Melanoma", "Riaz (Melanoma)"),
    ("Liu-iAtlas", "Melanoma", "Liu (Melanoma)"),
    ("Gide-iAtlas", "Melanoma", "Gide (Melanoma)"),
    ("Rosenberg-iAtlas", "Urothelial Bladder", "Rosenberg (Bladder)"),
    ("Padron-iAtlas", "PDAC", "Padron (PDAC)"),
    ("Anders-iAtlas", "Breast", "Anders (TNBC)"),
    ("McDermott-iAtlas", "ccRCC", "McDermott (ccRCC)"),
    ("Choueiri-iAtlas", "ccRCC", "Choueiri (ccRCC)"),
    ("Jung", "NSCLC", "Jung (NSCLC)"),
    ("Cho", "NSCLC", "Cho (NSCLC)"),
)


def classify_meta42_macrophage_tier(meta_state: str) -> Meta42MacrophageTier:
    """Classify a Consensus-Meta42 state into its macrophage/myeloid proximity tier."""
    st = meta_state.strip()
    if st in CORE_TAM_META42_IDS:
        return Meta42MacrophageTier.CORE_TAM
    if st in MONO_MDSC_META42_IDS:
        return Meta42MacrophageTier.MONO_MDSC_CONTINUUM
    if st in OTHER_MYELOID_META42_IDS:
        return Meta42MacrophageTier.OTHER_MYELOID
    return Meta42MacrophageTier.NON_MYELOID_CONTEXT


def classify_meta42_compartment(meta_state: str) -> LineageCompartment:
    """Map a Consensus-Meta42 state ID to its canonical cellular compartment."""
    return META42_COMPARTMENT_MAP.get(meta_state.strip(), LineageCompartment.OTHER_UNSPECIFIED)


COULTON_ORDERED_CLUSTERS: Final[tuple[str, ...]] = (
    "0_AlvMac",
    "1_MetM2Mac",
    "2_C3Mac",
    "3_ICIMac1",
    "4_ICIMac2",
    "5_StressMac",
    "6_SPP1AREGMac",
    "7_IFNMac",
    "8_IFNGMac",
    "9_AngioMac",
    "10_InflamMac",
    "11_MetalloMac",
    "12_MBMMac",
    "13_CalciumMac",
    "14_ProliMac",
    "15_LYZMac",
    "16_ECMHomeoMac",
    "17_IFNMac3",
    "18_ECMMac",
    "19_ClassMono",
    "20_TDoub",
    "21_HemeMac",
    "22_IFNMac4",
    "23_NA",
)

# Exact published findings from Coulton et al., Nature Communications 15, 5665 (2024), PMC11226649
COULTON_LITERATURE_CATALOG: Final[Mapping[str, CoultonLiteratureSpec]] = {
    "0_AlvMac": CoultonLiteratureSpec(
        cluster_id=0,
        short_label="0_AlvMac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="Lung tissue-resident alveolar macrophages (lipid/surfactant clearance)",
        canonical_markers=("FABP4", "MCEMP1", "CD52", "MARCO", "PPARG"),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="Predominantly lung-resident macrophage pool",
    ),
    "1_MetM2Mac": CoultonLiteratureSpec(
        cluster_id=1,
        short_label="1_MetM2Mac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="Immunoregulatory iron/selenium-metabolic M2-polarized TAMs",
        canonical_markers=("SELENOP", "SLC40A1", "FOLR2", "LYVE1", "F13A1", "PLTP"),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="Perivascular/tissue-resident M2 regulatory phenotype (FOLR2+/LYVE1+)",
    ),
    "2_C3Mac": CoultonLiteratureSpec(
        cluster_id=2,
        short_label="2_C3Mac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="Complement C3+ & MHC-II antigen-presenting phagocytic TAMs",
        canonical_markers=("C3", "PLD4", "C1QB", "HLA-DPA1", "HLA-DPB1"),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="Maps to C1QB+ TAMs; enriched in primary melanoma/GBM vs brain metastases",
    ),
    "3_ICIMac1": CoultonLiteratureSpec(
        cluster_id=3,
        short_label="3_ICIMac1",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="TREM2+/SPP1+/RNASE1+ lipid-remodeling ICI-associated TAMs",
        canonical_markers=("SPP1", "RNASE1", "NUPR1", "TREM2"),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.RESPONDER_ENRICHED,
        cpi1000_fgsea_qval=Some(0.027749715367563134),
        mana_and_spatial_note="CPI1000+ responder-enriched signature (q = 0.0277); TREM2 blockade target",
    ),
    "4_ICIMac2": CoultonLiteratureSpec(
        cluster_id=4,
        short_label="4_ICIMac2",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="TREM2+/APOE+/APOC1+ lipid-metabolic ICI-associated TAMs",
        canonical_markers=("TREM2", "APOE", "APOC1"),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.RESPONDER_ENRICHED,
        cpi1000_fgsea_qval=Some(0.002756081126333197),
        mana_and_spatial_note="CPI1000+ responder-enriched signature (q = 0.00276)",
    ),
    "5_StressMac": CoultonLiteratureSpec(
        cluster_id=5,
        short_label="5_StressMac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="Heat-shock / unfolded-protein stress-inducible TAMs",
        canonical_markers=("HSPA6", "HSPA1B", "HSPA1A", "DNAJB1", "HSPB1", "BAG3"),
        is_ucell_gold_standard=True,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="Gold-standard bulk-specific TAM signature; NS in CPI1000+",
    ),
    "6_SPP1AREGMac": CoultonLiteratureSpec(
        cluster_id=6,
        short_label="6_SPP1AREGMac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="Chemokine/EGFR-ligand migratory & pro-metastatic TAMs",
        canonical_markers=("SPP1", "AREG", "CCL20", "CXCL3", "IL1B", "CXCL2", "CXCL8"),
        is_ucell_gold_standard=True,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="Gold-standard signature; maps to CXCL8+ TAMs; enriched in LUSC & liver mets",
    ),
    "7_IFNMac": CoultonLiteratureSpec(
        cluster_id=7,
        short_label="7_IFNMac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="CC-chemokine & interferon-primed inflammatory TAMs",
        canonical_markers=("CCL2", "CCL8", "CCL4L2", "CCL3L3", "SPP1"),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="Monocyte-recruiting CCL2+/CCL8+ inflammatory TAM state",
    ),
    "8_IFNGMac": CoultonLiteratureSpec(
        cluster_id=8,
        short_label="8_IFNGMac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="IFN-gamma-stimulated CXCL9+/CXCL10+ T-cell-recruiting TAMs",
        canonical_markers=("CXCL9", "CXCL10", "MMP9", "VAMP5", "STAT1", "GBP1"),
        is_ucell_gold_standard=True,
        cpi1000_direction=IciOutcomeDirection.RESPONDER_ENRICHED,
        cpi1000_fgsea_qval=Some(1.3715289e-11),
        mana_and_spatial_note=(
            "Top gold-standard responder TAM in CPI1000+ (q = 1.37e-11); enriched in MANA-high "
            "neoantigen T-cell tumors (q = 0.0596) and spatially adjacent to CD4/CD8 memory T cells"
        ),
    ),
    "9_AngioMac": CoultonLiteratureSpec(
        cluster_id=9,
        short_label="9_AngioMac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="Pro-angiogenic VEGFA+/VCAN+/THBS1+ monocyte-derived TAMs",
        canonical_markers=("VEGFA", "VCAN", "THBS1", "AREG", "EREG", "IL1B"),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="Ontogenetically adjacent to 19_ClassMono; drives tumor angiogenesis",
    ),
    "10_InflamMac": CoultonLiteratureSpec(
        cluster_id=10,
        short_label="10_InflamMac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="TNF+/IL1B+/CXCL8+ acute pro-inflammatory cytokine TAMs",
        canonical_markers=("CCL3L3", "CCL4L2", "CXCL8", "IL1B", "TNF", "CCL4", "CCL3"),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="NF-kappaB inflammatory chemokine hub",
    ),
    "11_MetalloMac": CoultonLiteratureSpec(
        cluster_id=11,
        short_label="11_MetalloMac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="Metallothionein zinc-homeostasis & redox-buffering TAMs",
        canonical_markers=("MT1G", "MT1X", "MT2A", "MT1E", "MT1H", "MT1F", "MT1M"),
        is_ucell_gold_standard=True,
        cpi1000_direction=IciOutcomeDirection.RESPONDER_ENRICHED,
        cpi1000_fgsea_qval=Some(0.036466593541472025),
        mana_and_spatial_note="Gold-standard bulk-specific signature; responder-enriched in CPI1000+ (q = 0.0365)",
    ),
    "12_MBMMac": CoultonLiteratureSpec(
        cluster_id=12,
        short_label="12_MBMMac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="LRMDA+ melanoma brain metastasis niche-specialized TAMs",
        canonical_markers=("LRMDA",),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="Strongly enriched in melanoma brain metastases (q = 1.02e-9)",
    ),
    "13_CalciumMac": CoultonLiteratureSpec(
        cluster_id=13,
        short_label="13_CalciumMac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="FN1+ & S100 calcium-binding matrix/vesicle TAMs",
        canonical_markers=("FN1", "S100A6", "S100A10", "S100A4"),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="Extracellular-vesicle FN1 & S100A4/A6/A10 signaling state",
    ),
    "14_ProliMac": CoultonLiteratureSpec(
        cluster_id=14,
        short_label="14_ProliMac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="Cell-cycle proliferating MKI67+/TOP2A+ self-renewing TAMs",
        canonical_markers=("MKI67", "TOP2A", "CDK1", "CDKN3", "CENPF", "STMN1"),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.RESPONDER_ENRICHED,
        cpi1000_fgsea_qval=Some(6.084325e-12),
        mana_and_spatial_note="Responder-enriched in CPI1000+ (q = 6.08e-12; reflects proliferative immune-hot TME)",
    ),
    "15_LYZMac": CoultonLiteratureSpec(
        cluster_id=15,
        short_label="15_LYZMac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="Lysozyme-high antibacterial / homeostatic myeloid TAMs",
        canonical_markers=("LYZ",),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="Enriched in LUAD vs LUSC (q = 0.080)",
    ),
    "16_ECMHomeoMac": CoultonLiteratureSpec(
        cluster_id=16,
        short_label="16_ECMHomeoMac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="MMP9+/TIMP1+ balanced matrix-turnover TAMs",
        canonical_markers=("MMP9", "TIMP1", "SPP1"),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="Maps to SPP1+ TAMs in oral cancer; enriched in LUSC (q = 0.0317)",
    ),
    "17_IFNMac3": CoultonLiteratureSpec(
        cluster_id=17,
        short_label="17_IFNMac3",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="ISG15+/CXCL10+/CCL8+ interferon-stimulated antiviral TAMs",
        canonical_markers=("ISG15", "CXCL10", "CCL8", "IFIT1", "IFIT3"),
        is_ucell_gold_standard=True,
        cpi1000_direction=IciOutcomeDirection.RESPONDER_ENRICHED,
        cpi1000_fgsea_qval=Some(1.7930358164e-8),
        mana_and_spatial_note="Gold-standard responder TAM in CPI1000+ (q = 1.79e-8)",
    ),
    "18_ECMMac": CoultonLiteratureSpec(
        cluster_id=18,
        short_label="18_ECMMac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="Collagen-depositing COL1A1+/COL1A2+/SPARC+ ECM-remodeling TAMs",
        canonical_markers=("COL1A2", "COL1A1", "COL3A1", "SPARC", "COL6A1", "COL6A2", "COL6A3"),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.NON_RESPONDER_ENRICHED,
        cpi1000_fgsea_qval=Some(3.8213695118505e-5),
        mana_and_spatial_note=(
            "Primary ICI resistance TAM in CPI1000+ (q = 3.82e-5); enriched in MANA-low T-cell "
            "dysfunction tumors (q = 0.0780) and spatially co-localized with fibroblasts (CosMx)"
        ),
    ),
    "19_ClassMono": CoultonLiteratureSpec(
        cluster_id=19,
        short_label="19_ClassMono",
        role=CoultonClusterRole.CLASSICAL_MONOCYTE,
        functional_program="S100A8+/S100A9+/S100A12+/VCAN+ classical infiltrating monocytes",
        canonical_markers=("S100A8", "S100A9", "S100A12", "VCAN", "LYZ", "FCN1"),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="Maps to CD14+ classical monocytes; precursor to infiltrating TAMs",
    ),
    "20_TDoub": CoultonLiteratureSpec(
        cluster_id=20,
        short_label="20_TDoub",
        role=CoultonClusterRole.T_MYELOID_DOUBLET_CONTROL,
        functional_program="NKG7+/TRBC2+/CD3D+ T-cell / macrophage doublet control cluster",
        canonical_markers=("NKG7", "TRBC2", "CD3D", "IGKC"),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.RESPONDER_ENRICHED,
        cpi1000_fgsea_qval=Some(0.001668273617609862),
        mana_and_spatial_note="Retained by Coulton et al. as T-cell doublet positive control (q = 0.00167)",
    ),
    "21_HemeMac": CoultonLiteratureSpec(
        cluster_id=21,
        short_label="21_HemeMac",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="CD163+/HMOX1+/SLC40A1+ heme-clearance iron-scavenging TAMs",
        canonical_markers=("CD163", "HMOX1", "SLC40A1", "SELENOP"),
        is_ucell_gold_standard=True,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="Gold-standard bulk-specific heme-clearance TAM signature",
    ),
    "22_IFNMac4": CoultonLiteratureSpec(
        cluster_id=22,
        short_label="22_IFNMac4",
        role=CoultonClusterRole.TRUE_TAM_SUBSTATE,
        functional_program="IFITM2+/LST1+ interferon-exposed patrolling/inflammatory TAMs",
        canonical_markers=("IFITM2", "LST1", "FCGR3A", "AIF1"),
        is_ucell_gold_standard=True,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="Gold-standard bulk-specific interferon-exposed TAM signature",
    ),
    "23_NA": CoultonLiteratureSpec(
        cluster_id=23,
        short_label="23_NA",
        role=CoultonClusterRole.UNANNOTATED_ARTIFACT,
        functional_program="Residual 2-cell unannotated cluster excluded from substate rankings",
        canonical_markers=(),
        is_ucell_gold_standard=False,
        cpi1000_direction=IciOutcomeDirection.NEUTRAL_OR_NS,
        cpi1000_fgsea_qval=Nothing,
        mana_and_spatial_note="Tiny n=2 artifact cluster",
    ),
}
