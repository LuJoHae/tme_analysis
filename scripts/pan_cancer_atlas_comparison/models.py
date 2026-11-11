"""Immutable domain models, exhaustive sum types, Protocols, and Triad ontology constants.

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
from returns.maybe import Maybe, Nothing
from returns.result import Result

FloatMatrix: TypeAlias = npt.NDArray[np.float64]
FloatVector: TypeAlias = npt.NDArray[np.float64]
IntVector: TypeAlias = npt.NDArray[np.int64]
BoolVector: TypeAlias = npt.NDArray[np.bool_]
AltairChart: TypeAlias = alt.Chart | alt.LayerChart | alt.HConcatChart | alt.VConcatChart


# =============================================================================
# 1. Exhaustive Sum Types (StrEnum)
# =============================================================================


class TriadCompartment(StrEnum):
    """Canonical broad cellular lineage compartments across Consensus-Meta42 and the Triad atlases."""

    T_NK = "T_NK"
    STROMAL_CAF = "STROMAL_CAF"
    MYELOID_TAM = "MYELOID_TAM"
    B_PLASMA = "B_PLASMA"
    ENDOTHELIAL = "ENDOTHELIAL"
    EPITHELIAL_TUMOR = "EPITHELIAL_TUMOR"
    OTHER = "OTHER"


class ExternalAtlasId(StrEnum):
    """Peer-reviewed pan-cancer external single-cell reference atlases across the TME Triad."""

    ZHENG_TCELL = "Zheng_TCell_Atlas"
    CHU_TCELL = "Chu_TCell_Atlas"
    SHI_TCELL = "Shi_PanCancer_Atlas_TNK"
    GAO_CAF = "Gao_CAF_Atlas"
    CORDS_CAF = "Cords_CAF_Atlas"
    BUECHLER_CAF = "Buechler_CAF_Atlas"
    SHI_CAF = "Shi_PanCancer_Atlas_CAF"
    COULTON_TAM = "Coulton_TAM_Atlas"


class ProfileOrigin(StrEnum):
    """Provenance of a substate expression profile."""

    EMPIRICAL_SINGLE_CELL_H5AD = "empirical_single_cell_h5ad"
    CURATED_SIGNATURE_PRIOR = "curated_signature_prior"


class IciDirection(StrEnum):
    """Clinical immune checkpoint blockade association direction."""

    PRO_RESPONSE = "Responder-Enriched"
    PRO_RESISTANCE = "Non-Responder-Enriched"
    NEUTRAL = "Neutral / Sub-Threshold"


class SignificanceTier(StrEnum):
    """Standard four-tier statistical significance scale (*, **, ***, ****)."""

    NS = "ns"
    ONE_STAR = "*"
    TWO_STARS = "**"
    THREE_STARS = "***"
    FOUR_STARS = "****"


def atlas_target_compartment(atlas_id: ExternalAtlasId) -> TriadCompartment:
    """Exhaustively map each external atlas identifier to its primary TME Triad compartment."""
    match atlas_id:
        case ExternalAtlasId.ZHENG_TCELL | ExternalAtlasId.CHU_TCELL | ExternalAtlasId.SHI_TCELL:
            return TriadCompartment.T_NK
        case (
            ExternalAtlasId.GAO_CAF
            | ExternalAtlasId.CORDS_CAF
            | ExternalAtlasId.BUECHLER_CAF
            | ExternalAtlasId.SHI_CAF
        ):
            return TriadCompartment.STROMAL_CAF
        case ExternalAtlasId.COULTON_TAM:
            return TriadCompartment.MYELOID_TAM
        case _ as unreachable:
            assert_never(unreachable)


def atlas_citation_label(atlas_id: ExternalAtlasId) -> str:
    """Human-readable journal citation label for each external reference atlas."""
    match atlas_id:
        case ExternalAtlasId.ZHENG_TCELL:
            return "Zheng et al. (Science 2021)"
        case ExternalAtlasId.CHU_TCELL:
            return "Chu et al. (Nat. Med. 2023)"
        case ExternalAtlasId.SHI_TCELL:
            return "Shi et al. (Nature 2025, T/NK)"
        case ExternalAtlasId.GAO_CAF:
            return "Gao et al. (Cancer Cell 2024)"
        case ExternalAtlasId.CORDS_CAF:
            return "Cords et al. (Cancer Cell 2024)"
        case ExternalAtlasId.BUECHLER_CAF:
            return "Buechler et al. (Nature 2021)"
        case ExternalAtlasId.SHI_CAF:
            return "Shi et al. (Nature 2025, CAF)"
        case ExternalAtlasId.COULTON_TAM:
            return "Coulton et al. (Nat. Commun. 2024)"
        case _ as unreachable:
            assert_never(unreachable)


def atlas_short_prefix(atlas_id: ExternalAtlasId) -> str:
    """Concise prefix for combined cross-atlas axis labels."""
    match atlas_id:
        case ExternalAtlasId.ZHENG_TCELL:
            return "Zheng"
        case ExternalAtlasId.CHU_TCELL:
            return "Chu"
        case ExternalAtlasId.SHI_TCELL:
            return "Shi-T"
        case ExternalAtlasId.GAO_CAF:
            return "Gao"
        case ExternalAtlasId.CORDS_CAF:
            return "Cords"
        case ExternalAtlasId.BUECHLER_CAF:
            return "Buechler"
        case ExternalAtlasId.SHI_CAF:
            return "Shi-CAF"
        case ExternalAtlasId.COULTON_TAM:
            return "Coulton"
        case _ as unreachable:
            assert_never(unreachable)


# =============================================================================
# 2. Frozen Domain Models
# =============================================================================


class PanAtlasConfig(BaseModel):
    """Immutable configuration for the Pan-Cancer T-Cell, CAF & Myeloid Triad comparison pipeline."""

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
    coulton_centroids_path: Path = Field(
        default=Path("output/results/coulton_tam_comparison/coulton_substate_centroids.parquet")
    )
    output_dir: Path = Field(default=Path("output/results/pan_cancer_atlas_comparison"))
    figures_dir: Path = Field(default=Path("output/figures/pan_cancer_atlas_comparison"))
    batch_size: int = Field(default=25_000, ge=1_000)
    top_k_markers: int = Field(default=25, ge=5, le=100)
    bulk_signature_genes: int = Field(default=12, ge=4, le=40)
    min_cells_per_substate: int = Field(default=25, ge=5)
    reuse_cached_centroids: bool = Field(default=True)
    render_plots: bool = Field(default=True)


class AtlasCentroidAccumulator(BaseModel):
    """Immutable running accumulator for out-of-core sparse CSR cluster centroid reduction."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    atlas_id: ExternalAtlasId
    origin: ProfileOrigin
    cluster_labels: tuple[str, ...]
    gene_symbols: tuple[str, ...]
    all_cell_counts: IntVector
    tumor_cell_counts: IntVector
    all_linear_sum: FloatMatrix
    tumor_linear_sum: FloatMatrix
    all_log_sum: FloatMatrix
    tumor_log_sum: FloatMatrix
    all_detect_counts: FloatMatrix
    tumor_detect_counts: FloatMatrix


class AtlasCentroidBundle(BaseModel):
    """Finalized cluster expression profiles for one external atlas across Tumor-Only and All-Tissues."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    atlas_id: ExternalAtlasId
    origin: ProfileOrigin
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
    """Immutable standardized bulk RNA-seq clinical ICB cohort (`N = 1,058` across 11 cohorts)."""

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


class StromalGatingResult(BaseModel):
    """Empirical conditional `TGF-beta myCAF` x `Tumor-Reactive CD8 T` gating across `N = 1,058` ICB patients."""

    model_config = ConfigDict(frozen=True)

    tcell_axis_label: str
    caf_axis_label: str
    n_patients: int
    n_cohorts: int
    beta_tcell_main: float
    se_tcell_main: float
    pval_tcell_main: float
    beta_caf_main: float
    se_caf_main: float
    pval_caf_main: float
    beta_interaction: float
    se_interaction: float
    pval_interaction: float
    or_tcell_low_caf: float
    pval_tcell_low_caf: float
    or_tcell_high_caf: float
    pval_tcell_high_caf: float
    resp_rate_low_t_low_caf: float
    n_low_t_low_caf: int
    resp_rate_high_t_low_caf: float
    n_high_t_low_caf: int
    resp_rate_low_t_high_caf: float
    n_low_t_high_caf: int
    resp_rate_high_t_high_caf: float
    n_high_t_high_caf: int


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


class PanAtlasOutputs(BaseModel):
    """Complete pure deliverables of the Pan-Cancer T-Cell, CAF & Myeloid Triad comparison pipeline."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    atlas_centroids_df: pl.DataFrame
    atlas_top_markers_df: pl.DataFrame
    meta42_vs_atlases_similarity_df: pl.DataFrame
    nearest_tcell_matches_df: pl.DataFrame
    nearest_caf_matches_df: pl.DataFrame
    cross_atlas_equivalence_df: pl.DataFrame
    atlas_ici_outcomes_df: pl.DataFrame
    stromal_gating_df: pl.DataFrame
    triad_block_summary_df: pl.DataFrame
    markdown_report: str


# =============================================================================
# 3. Typeclass Protocols (Side-Effecting IO Boundaries)
# =============================================================================


class ExternalAtlasStreamer(Protocol):
    """Streams cluster centroids from an external atlas strictly via `tme_datasets`."""

    def __call__(
        self,
        atlas_id: ExternalAtlasId,
        batch_size: int,
        min_cells: int,
    ) -> Result[AtlasCentroidBundle, str]: ...


class Meta42ArtifactLoader(Protocol):
    """Loads Consensus-Meta42 reference profiles and ICI effect parquets."""

    def __call__(self, config: PanAtlasConfig, /) -> Result[Meta42ReferenceBundle, str]: ...


class BulkIcbSuiteLoader(Protocol):
    """Loads clinical bulk RNA-seq ICB cohorts strictly via `tme_datasets.load_dataset`."""

    def __call__(self) -> Result[tuple[BulkIcbCohortBundle, ...], str]: ...


class VectorChartExporter(Protocol):
    """Exports a declarative Altair chart to Inkscape-layered `.svg` and 300 DPI `.png`."""

    def __call__(self, chart: AltairChart, base_path: Path, /) -> Result[tuple[Path, Path], str]: ...


# =============================================================================
# 4. Canonical `Consensus-Meta42` & Clinical Bulk ICB Constants
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

META42_TRIAD_MAP: Final[Mapping[str, TriadCompartment]] = {
    "M01": TriadCompartment.T_NK,
    "M02": TriadCompartment.STROMAL_CAF,
    "M03": TriadCompartment.T_NK,
    "M04": TriadCompartment.T_NK,
    "M05": TriadCompartment.T_NK,
    "M06": TriadCompartment.T_NK,
    "M07": TriadCompartment.OTHER,
    "M08": TriadCompartment.T_NK,
    "M09": TriadCompartment.MYELOID_TAM,
    "M10": TriadCompartment.MYELOID_TAM,
    "M11": TriadCompartment.EPITHELIAL_TUMOR,
    "M12": TriadCompartment.EPITHELIAL_TUMOR,
    "M13": TriadCompartment.T_NK,
    "M14": TriadCompartment.MYELOID_TAM,
    "M15": TriadCompartment.T_NK,
    "M16": TriadCompartment.ENDOTHELIAL,
    "M17": TriadCompartment.STROMAL_CAF,
    "M18": TriadCompartment.EPITHELIAL_TUMOR,
    "M19": TriadCompartment.MYELOID_TAM,
    "M20": TriadCompartment.T_NK,
    "M21": TriadCompartment.T_NK,
    "M22": TriadCompartment.EPITHELIAL_TUMOR,
    "M23": TriadCompartment.STROMAL_CAF,
    "M24": TriadCompartment.T_NK,
    "M25": TriadCompartment.T_NK,
    "M26": TriadCompartment.EPITHELIAL_TUMOR,
    "M27": TriadCompartment.EPITHELIAL_TUMOR,
    "M28": TriadCompartment.MYELOID_TAM,
    "M29": TriadCompartment.B_PLASMA,
    "M30": TriadCompartment.B_PLASMA,
    "M31": TriadCompartment.MYELOID_TAM,
    "M32": TriadCompartment.MYELOID_TAM,
    "M33": TriadCompartment.T_NK,
    "M34": TriadCompartment.T_NK,
    "M35": TriadCompartment.MYELOID_TAM,
    "M36": TriadCompartment.T_NK,
    "M37": TriadCompartment.B_PLASMA,
    "M38": TriadCompartment.MYELOID_TAM,
    "M39": TriadCompartment.MYELOID_TAM,
    "M40": TriadCompartment.T_NK,
    "M41": TriadCompartment.MYELOID_TAM,
    "M42": TriadCompartment.MYELOID_TAM,
}

META42_T_NK_IDS: Final[tuple[str, ...]] = (
    "M13",
    "M25",
    "M15",
    "M34",
    "M33",
    "M04",
    "M06",
    "M36",
    "M20",
    "M01",
    "M05",
    "M24",
    "M03",
    "M21",
    "M40",
    "M08",
)

META42_CORE_CAF_IDS: Final[tuple[str, ...]] = ("M17", "M23", "M02")
META42_CAF_AND_MESENCHYMAL_CONTEXT_IDS: Final[tuple[str, ...]] = (
    "M17",
    "M23",
    "M02",
    "M26",
    "M16",
    "M12",
)
META42_MYELOID_IDS: Final[tuple[str, ...]] = (
    "M14",
    "M42",
    "M32",
    "M10",
    "M28",
    "M35",
    "M19",
    "M41",
    "M09",
    "M31",
    "M39",
    "M38",
)

PAN_TCELL_CORE_GENES: Final[tuple[str, ...]] = (
    "CD3D",
    "CD3E",
    "CD3G",
    "TRAC",
    "TRBC1",
    "TRBC2",
    "CD2",
    "LCK",
)

PAN_CAF_CORE_GENES: Final[tuple[str, ...]] = (
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "SPARC",
    "COL6A1",
    "VIM",
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


def classify_meta42_triad_compartment(meta_state: str) -> TriadCompartment:
    """Classify a Consensus-Meta42 state into its canonical Triad compartment."""
    return META42_TRIAD_MAP.get(meta_state.strip(), TriadCompartment.OTHER)
