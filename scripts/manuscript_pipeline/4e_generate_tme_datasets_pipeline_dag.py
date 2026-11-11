#!/usr/bin/env python3
"""Generate Publication Vector Diagram: tme_datasets Package & 37-Script Provenance DAG.

Constructs a 180.0mm double-column vector diagram adhering strictly to Nature Methods
minimalist wireframe and Inkscape multi-layer standards:
- Pure white canvas (#FFFFFF), sharp square corners (rx="0"), 0.75pt hairline borders (#CBD5E1),
  and zero drop shadows or filter glows.
- Okabe-Ito colorblind-safe accents strictly applied to functional top-bars, node ports, and arrows.
- Native Inkscape layers (inkscape:groupmode="layer") and semantic groups (<g id="...">).
- Minimum font size >= 13.5px across all labels.
- Tier I maps all 6 domains in `packages/tme_datasets` (`SC-9`, `BULK-15`, `SC-REF`, `TCGA-5`,
  `IO-SETS`, and `PIPE` = `tme_datasets.pipeline`).
- Zero arrows are drawn for `tme_datasets` API loads; instead, each script card explicitly
  states which `tme_datasets` API function(s) it calls, while directed arrows strictly denote
  script-to-script `.parquet` and `.svg` artifact handoffs across the 37 scripts (`2a1`–`4e`).
"""

from __future__ import annotations

import html
from pathlib import Path
import subprocess
import sys
from typing import Final
import xml.etree.ElementTree as ET

from pydantic import BaseModel, ConfigDict
from returns.pipeline import is_successful
from returns.result import Failure, Result, Success
import vl_convert as vlc

# Ensure repository root is importable for shared Nature Methods style tokens
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.figures.nature_style_config import (  # noqa: E402
    COLOR_BORDER_HAIRLINE,
    COLOR_CANVAS_BG,
    COLOR_DIVIDER_RULE,
    COLOR_PANEL_BG,
    COLOR_SUBTLE_FILL,
    COLOR_TEXT_MUTED,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    FONT_SANS,
    OKABE_BLUE,
    OKABE_BLUISH_GREEN,
    OKABE_ORANGE,
    OKABE_REDDISH_PURPLE,
    OKABE_SKY_BLUE,
    OKABE_VERMILION,
)

VIEWBOX_W: Final[int] = 1400
VIEWBOX_H: Final[int] = 1664
FIG_WIDTH_MM: Final[float] = 180.0
FIG_HEIGHT_MM: Final[float] = round(FIG_WIDTH_MM * (VIEWBOX_H / VIEWBOX_W), 2)
MIN_FONT_PX: Final[float] = 13.5
API_PREFIX: Final[str] = "tme_datasets API: "


class DagCardSpec(BaseModel):
    """Immutable specification for a rectangular dataset or script card in the provenance DAG."""

    model_config = ConfigDict(frozen=True)

    card_id: str
    x: int
    y: int
    width: int
    height: int
    accent_color: str
    badge_tag: str
    title: str
    body_lines: tuple[str, ...]
    footer_label: str
    footer_color: str


class TierBandSpec(BaseModel):
    """Immutable specification for a horizontal tier container band."""

    model_config = ConfigDict(frozen=True)

    band_id: str
    letter: str
    heading: str
    right_summary: str
    y: int
    height: int


class ConnectorSpec(BaseModel):
    """Immutable specification for a gutter-routed orthogonal connector arrow."""

    model_config = ConfigDict(frozen=True)

    edge_id: str
    points: tuple[tuple[int, int], ...]
    stroke_color: str
    marker_id: str
    dashed: bool = False


class PillLabelSpec(BaseModel):
    """Immutable specification for a compact data-artifact pill label placed in a gutter."""

    model_config = ConfigDict(frozen=True)

    pill_id: str
    cx: int
    cy: int
    width: int
    text: str
    border_color: str
    text_color: str


class DagExportConfig(BaseModel):
    """Immutable export paths for the provenance DAG figure."""

    model_config = ConfigDict(frozen=True)

    results_svg_path: Path = Path(
        "output/manuscript_pipeline/4e_generate_tme_datasets_pipeline_dag/figure_tme_datasets_pipeline_dag.svg"
    )
    results_png_path: Path = Path(
        "output/manuscript_pipeline/4e_generate_tme_datasets_pipeline_dag/figure_tme_datasets_pipeline_dag.png"
    )
    article_svg_path: Path = Path("article/figures/dataset_overview/figure_tme_datasets_pipeline_dag.svg")
    article_png_path: Path = Path("article/figures/dataset_overview/figure_tme_datasets_pipeline_dag.png")


TIER_BANDS: Final[tuple[TierBandSpec, ...]] = (
    TierBandSpec(
        band_id="band-tier-1",
        letter="a",
        heading="Tier I — packages/tme_datasets: Registry, Asset Pools & Ingestion/QC Pipeline",
        right_summary="Ensembl v111 · RECIST v1.1 · API-Only Preprocessed Access",
        y=82,
        height=344,
    ),
    TierBandSpec(
        band_id="band-tier-2",
        letter="b",
        heading="Tier II: Primary Cohort Analysis Scripts (3 scripts)",
        right_summary="scripts/manuscript_pipeline/2a1, 2a2, 2b",
        y=446,
        height=198,
    ),
    TierBandSpec(
        band_id="band-tier-3",
        letter="c",
        heading="Tier III: Chained Compute Engines",
        right_summary="",
        y=702,
        height=530,
    ),
    TierBandSpec(
        band_id="band-tier-4",
        letter="d",
        heading="Tier IV: Article Figure & Table Scripts",
        right_summary="",
        y=1286,
        height=362,
    ),
)

TIER1_CARDS: Final[tuple[DagCardSpec, ...]] = (
    DagCardSpec(
        card_id="card-t1-sc9",
        x=34,
        y=122,
        width=440,
        height=138,
        accent_color=OKABE_BLUE,
        badge_tag="SC-9",
        title="9 Single-Cell ICB Cohorts (3.60M Cells)",
        body_lines=(
            "GSE120575 (Melanoma), CELLxGENE_7b20c613 (Pan-Ca),",
            "05a8c945 (CRC), 6f9de485 (TNBC), GSE207422 / 243013 /",
            "233203 (NSCLC), GSE200996 (HNSCC), GSE316195 (PDAC)",
        ),
        footer_label="API: load_dataset() · load_empirical_pseudobulk() · load_backed()",
        footer_color=OKABE_BLUE,
    ),
    DagCardSpec(
        card_id="card-t1-bulk15",
        x=492,
        y=122,
        width=424,
        height=138,
        accent_color=OKABE_BLUISH_GREEN,
        badge_tag="BULK-15",
        title="11 Bulk ICB Trials + 4 Pools (N = 1,058)",
        body_lines=(
            "9 iAtlas (N = 1,015): Hugo, Riaz, Liu, Gide, Rosenberg,",
            "Padron, Anders, McDermott, Choueiri (providers.bulk_iatlas)",
            "2 NSCLC (N = 43): Jung, Cho + 4 Indication/Pan-Ca Pools",
        ),
        footer_label="API: load_dataset() · load_iatlas_cohort_or_combined()",
        footer_color=OKABE_BLUISH_GREEN,
    ),
    DagCardSpec(
        card_id="card-t1-scref",
        x=934,
        y=122,
        width=432,
        height=138,
        accent_color=OKABE_SKY_BLUE,
        badge_tag="SC-REF",
        title="17+ Single-Cell Reference Atlases",
        body_lines=(
            "GSE179994 (NSCLC T), Maynard_NSCLC, GSE115978 (Mel.),",
            "GSE125449 (HCC/iCCA), GSE123813 (BCC/SCC) + 12 Solid",
            "Tumor Atlases (providers.single_cell & cellxgene)",
        ),
        footer_label="API: load_dataset() · sample_single_cell_cohorts()",
        footer_color=OKABE_BLUE,
    ),
    DagCardSpec(
        card_id="card-t1-tcga5",
        x=34,
        y=274,
        width=440,
        height=138,
        accent_color=OKABE_VERMILION,
        badge_tag="TCGA-5",
        title="5 TCGA Primary Baselines (N = 2,932)",
        body_lines=(
            "Untreated primary tumor RNA-seq baseline cohorts:",
            "TCGA-SKCM, TCGA-BLCA, TCGA-KIRC, TCGA-BRCA, TCGA-PAAD",
            "Harmonized TPM & tumor-purity background (providers.tcga)",
        ),
        footer_label='API: load_dataset("TCGA_*") · list_registered_datasets()',
        footer_color=OKABE_VERMILION,
    ),
    DagCardSpec(
        card_id="card-t1-iosets",
        x=492,
        y=274,
        width=424,
        height=138,
        accent_color=OKABE_REDDISH_PURPLE,
        badge_tag="IO-SETS",
        title="IO Gene Panels & Biomarker Signatures",
        body_lines=(
            "IMMUNOTHERAPY_GENE_PANEL (101 curated IO symbols)",
            "IMMUNE_CHECKPOINT_GENES · Bagaev 29 TME Signatures",
            "TIDE · Ayers IFN-γ · Jiang Texh · 32-Biomarker Suite",
        ),
        footer_label="API: load_geneset_collection() · IO Gene Panels",
        footer_color=OKABE_REDDISH_PURPLE,
    ),
    DagCardSpec(
        card_id="card-t1-pipe",
        x=934,
        y=274,
        width=432,
        height=138,
        accent_color=OKABE_ORANGE,
        badge_tag="PIPE",
        title="tme_datasets.pipeline (Ingest & QC)",
        body_lines=(
            "download_preprocess · batch_preprocess · qc_orchestrator",
            "Ensembl v111 symbol reconciliation + Luecken–Theis MAD",
            "outlier/doublet QC -> data/{preprocessed,pseudobulk}/",
        ),
        footer_label="CLI: tme-datasets-{ingest,batch-preprocess,qc}",
        footer_color="#B45309",
    ),
)

TIER2_CARDS: Final[tuple[DagCardSpec, ...]] = (
    DagCardSpec(
        card_id="card-t2-milopy",
        x=34,
        y=488,
        width=654,
        height=144,
        accent_color=OKABE_BLUE,
        badge_tag="2A",
        title="2a1_sample_and_analyze_cohorts.py + 2a2_combine_cohorts_and_analyze.py",
        body_lines=(
            "tme_datasets API: resolve_milopy_cohorts(), load_dataset(), align_and_concatenate()",
            "Fits patient-aware Harmony + milopy kNN Negative Binomial GLM across SC-9 & 4 pools",
            "Projects kNN logFC onto single cells (nhood_matrix @ da_logFC -> cell_level_scores)",
        ),
        footer_label="Out: output/manuscript_pipeline/2a1_sample_and_analyze_cohorts/*.parquet & 2a2_*/",
        footer_color=OKABE_BLUE,
    ),
    DagCardSpec(
        card_id="card-t2-stability",
        x=712,
        y=488,
        width=654,
        height=144,
        accent_color=OKABE_REDDISH_PURPLE,
        badge_tag="2B",
        title="2b_run_stability_selection_iatlas.py",
        body_lines=(
            "tme_datasets API: load_iatlas_cohort_or_combined(), IMMUNOTHERAPY_GENE_PANEL",
            "Evaluates 12 bulk cohorts/pools × 13 fitters (Lasso, LassoFWL, MERF, OWL, SLOPE)",
            "Complementary-pairs subsampling (B = 50) under finite-sample PFER q_budget",
        ),
        footer_label="Out: output/manuscript_pipeline/2b_run_stability_selection_iatlas/{scores,paths}.parquet",
        footer_color=OKABE_REDDISH_PURPLE,
    ),
)

TIER3_CARDS: Final[tuple[DagCardSpec, ...]] = (
    DagCardSpec(
        card_id="card-t3-pseudobulk",
        x=34,
        y=742,
        width=440,
        height=148,
        accent_color=OKABE_BLUE,
        badge_tag="3A",
        title="3a_cross_cohort_pseudobulk_deconv_benchmark.py",
        body_lines=(
            "tme_datasets API: load_dataset(), load_empirical_pseudobulk()",
            "Pairs SC-9 with cell_level_scores.parquet from [2A1]",
            "Builds 355 Leiden r=2.0 reference centroids (Phi_1..Phi_9)",
            "Evaluates 5-sweep Milo vs. pseudobulk β concordance",
        ),
        footer_label="Out: .../3a_cross_cohort_pseudobulk_deconv_benchmark/",
        footer_color=OKABE_BLUE,
    ),
    DagCardSpec(
        card_id="card-t3-unified",
        x=34,
        y=906,
        width=440,
        height=148,
        accent_color=OKABE_SKY_BLUE,
        badge_tag="3B",
        title="3b_pan_cohort_unified_deconv_benchmark.py",
        body_lines=(
            "tme_datasets API: load_empirical_pseudobulk() [SC-9 cache]",
            "Reads 355 centroids & 5-sweep Milo betas from [3A]",
            "Ward-clusters 355 SC states -> 42 Consensus-Meta42",
            "Verifies Unified-All9, Meta42 & LOCO concordance",
        ),
        footer_label="Out: .../3b_pan_cohort_unified_deconv_benchmark/",
        footer_color=OKABE_BLUE,
    ),
    DagCardSpec(
        card_id="card-t3-multires",
        x=34,
        y=1070,
        width=440,
        height=148,
        accent_color=OKABE_SKY_BLUE,
        badge_tag="3F",
        title="3f1–3f8 Multi-Res Reference Scaling Suite",
        body_lines=(
            "tme_datasets API: load_dataset() [SC-9, SC-REF, BULK-15]",
            "3f1..3f4: SF / Harmony / {2,3,4,16}-dataset references",
            "3f5..3f6: Multi-res InstaPrism deconv & L2 ROC-AUC",
            "3f7..3f8: Multi-res Milo DA & 18-strategy saturation",
        ),
        footer_label="Out: .../3f1_*/..3f8_*/ (multi_resolution_benchmark_*)",
        footer_color=OKABE_BLUE,
    ),
    DagCardSpec(
        card_id="card-t3-iatlas-bench",
        x=492,
        y=742,
        width=424,
        height=476,
        accent_color=OKABE_BLUISH_GREEN,
        badge_tag="3C",
        title="3c_iatlas_instaprism_milopy_benchmark.py",
        body_lines=(
            "tme_datasets API: load_dataset() [11 Bulk Trials]",
            "Upstream: 355 SC [3A] & Meta-42 [3B] .parquet refs",
            "1. Detects log2 bulk scale, linearizes (2^X - 1 -> TPM)",
            "2. InstaPrism EM across 13 Refs (Tiers 1–3) × 15 Cohorts",
            "3. Aitchison CLR + arcsin-sqrt compositional features",
            "4. Evaluates 5 ML families (Ridge, ENet, Lasso, OLS, RF)",
            "5. Zero-shot Leave-One-Cohort-Out (LOCO) across 11",
            "6. DerSimonian–Laird Pan-iAtlas-Meta β vs. Milo logFC",
            "7. Extracts 168 top-4 HGNC DEGs (14-state & All-42)",
            "8. Unblended CLR/SC-Markers, OOF Gating & Epistasis",
        ),
        footer_label="Out: .../3c_iatlas_instaprism_milopy_benchmark/ (16)",
        footer_color=OKABE_BLUISH_GREEN,
    ),
    DagCardSpec(
        card_id="card-t3-dossier",
        x=934,
        y=742,
        width=432,
        height=192,
        accent_color=OKABE_REDDISH_PURPLE,
        badge_tag="3D",
        title="3d_curate_expert_gene_list.py",
        body_lines=(
            "tme_datasets API: load_iatlas_cohort_or_combined()",
            "Reads stability {scores,paths}.parquet from [2B]",
            "Computes pan-cancer bulk fold-changes & Mann–Whitney p",
            "Synthesizes composite Evidence Index E(g) across 91 runs",
            "Exports curated biomarker dossier & Typst tables",
        ),
        footer_label="Out: .../3d_curate_expert_gene_list/ & table_*.typ",
        footer_color=OKABE_REDDISH_PURPLE,
    ),
    DagCardSpec(
        card_id="card-t3-predictors",
        x=934,
        y=950,
        width=432,
        height=268,
        accent_color=OKABE_ORANGE,
        badge_tag="3E",
        title="3e1–3e5 32-Biomarker, Survival & HPO Suite",
        body_lines=(
            "tme_datasets API: load_dataset() [BULK-15 & SC]",
            "tme_datasets API: sample_single_cell_cohorts() [3e2]",
            "3e1_benchmark_iatlas_icb_predictors (27 RNA / 31 all)",
            "3e2_run_hpo_pipeline · 3e3_benchmark_survival_and_dca",
            "3e4_update_texh_benchmark (Jiang_Texh orientation)",
            "3e5_integrate_deconvolution_our (reads [3C] AUC ->",
            "  32nd biomarker, 28-RNA PR-AUC & forest SVGs/PNGs)",
        ),
        footer_label="Out: .../3e1_*/..3e5_*/ (*.parquet, *.svg, *.png)",
        footer_color="#B45309",
    ),
)

TIER4_CARDS: Final[tuple[DagCardSpec, ...]] = (
    DagCardSpec(
        card_id="card-t4-sc-figs",
        x=34,
        y=1324,
        width=440,
        height=152,
        accent_color=OKABE_BLUE,
        badge_tag="4A",
        title="4a1–4a7 Single-Cell Volcano & 2×2 Figures",
        body_lines=(
            "tme_datasets API: load_backed() [4a1 UMAPs only]",
            "4a2: 9-cohort Milo volcano & logFC grids (reads [2A])",
            "4a3..4a6: UMAP, DA bar & scatter/ρ grids ([2A], [3A])",
            "4a7: Assembles 180mm 2×2 immunotherapy master figure",
        ),
        footer_label="Out: .../4a1_*/..4a7_*/ -> figures/single-cell/*",
        footer_color=OKABE_BLUE,
    ),
    DagCardSpec(
        card_id="card-t4-iatlas-figs",
        x=492,
        y=1324,
        width=424,
        height=152,
        accent_color=OKABE_VERMILION,
        badge_tag="4B",
        title="4b_plot_iatlas_instaprism_milopy_figures.py",
        body_lines=(
            "tme_datasets API: none (reads 16 [3C] .parquet tables)",
            "Renders 14 cross-modality & epistasis figures + .typ:",
            "Focused & All-42 Concordance/DEGs (9) + Heading-Free",
            "  Main Epistasis/Gating (1) & Supp S19–S22 (4)",
        ),
        footer_label="Out: .../4b_plot_iatlas_*/ -> figures/ & captions/",
        footer_color=OKABE_VERMILION,
    ),
    DagCardSpec(
        card_id="card-t4-stab-figs",
        x=934,
        y=1324,
        width=432,
        height=152,
        accent_color=OKABE_REDDISH_PURPLE,
        badge_tag="4C",
        title="4c1–4c6 Stability, HPO & 32-Biomarker Forest",
        body_lines=(
            "tme_datasets API: none (reads [2B]/[3D]/[3E] parquets)",
            "4c1..4c4: Stability paths, cohort & fitter overlaps",
            "4c5: Optuna HPO sampling & benchmark figures ([3E2])",
            "4c6: Combined 32-biomarker summary forest ([3E5])",
        ),
        footer_label="Out: .../4c1_*/..4c6_*/ (*.svg & 300 DPI *.png)",
        footer_color=OKABE_REDDISH_PURPLE,
    ),
    DagCardSpec(
        card_id="card-t4-multires-figs",
        x=34,
        y=1494,
        width=654,
        height=140,
        accent_color=OKABE_SKY_BLUE,
        badge_tag="4D",
        title="4d1–4d2 Multi-Res Validation & Reference Selection Figures",
        body_lines=(
            "tme_datasets API: none (reads [3F] multi-res .parquet tables)",
            "4d1_plot_sf_multires_validation_figures: 18-panel (a)–(r) saturation [3F]",
            "4d2_plot_cohort_auc_reference_selection: {1, 2, 3, 16}-dataset gain [3F8]",
        ),
        footer_label="Out: .../4d1_*/ & .../4d2_*/ (fig1_predictive_capacity_* & fig1_selected_*)",
        footer_color=OKABE_BLUE,
    ),
    DagCardSpec(
        card_id="card-t4-dag",
        x=712,
        y=1494,
        width=654,
        height=140,
        accent_color=OKABE_ORANGE,
        badge_tag="4E",
        title="4e_generate_tme_datasets_pipeline_dag.py, Snakefile & Makefile",
        body_lines=(
            "tme_datasets API: none (renders 37-script multi-layer provenance DAG)",
            "Snakefile + workflow/rules/manuscript_pipeline.smk (isolated stage outputs)",
            "article/Makefile: rsvg-convert 300 DPI rasterization -> Typst article/main.pdf",
        ),
        footer_label="Out: .../4e_generate_tme_datasets_pipeline_dag/ -> article/main.pdf",
        footer_color="#B45309",
    ),
)

# Orthogonal connectors strictly inside gutters and clear of all tier header text (x > 365).
# ZERO outgoing arrows from Tier I (packages/tme_datasets): only script-to-script handoffs!
CONNECTORS: Final[tuple[ConnectorSpec, ...]] = (
    # Tier II -> Tier III artifact handoffs
    ConnectorSpec(
        edge_id="edge-milopy-to-pseudobulk",
        points=((230, 632), (230, 692), (385, 692), (385, 740)),
        stroke_color=OKABE_BLUE,
        marker_id="arr-blue",
    ),
    ConnectorSpec(
        edge_id="edge-stability-to-dossier",
        points=((1150, 632), (1150, 740)),
        stroke_color=OKABE_REDDISH_PURPLE,
        marker_id="arr-purple",
    ),
    # Intra-Tier III script-to-script handoffs
    ConnectorSpec(
        edge_id="edge-pseudobulk-to-unified",
        points=((254, 890), (254, 904)),
        stroke_color=OKABE_BLUE,
        marker_id="arr-blue",
    ),
    ConnectorSpec(
        edge_id="edge-pseudobulk-to-iatlas-bench",
        points=((474, 816), (490, 816)),
        stroke_color=OKABE_BLUE,
        marker_id="arr-blue",
    ),
    ConnectorSpec(
        edge_id="edge-unified-to-iatlas-bench",
        points=((474, 980), (490, 980)),
        stroke_color=OKABE_BLUE,
        marker_id="arr-blue",
    ),
    ConnectorSpec(
        edge_id="edge-iatlas-bench-to-predictors",
        points=((916, 1084), (932, 1084)),
        stroke_color=OKABE_BLUISH_GREEN,
        marker_id="arr-teal",
        dashed=True,
    ),
    # Left outer gutter from [2A] & [3A] down to [4A] (enters at x >= 375, clear of Tier IV title)
    ConnectorSpec(
        edge_id="edge-milopy-to-sc-figs",
        points=((34, 560), (22, 560), (22, 1278), (375, 1278), (375, 1322)),
        stroke_color=OKABE_BLUE,
        marker_id="arr-blue",
        dashed=True,
    ),
    ConnectorSpec(
        edge_id="edge-pseudobulk-to-sc-figs",
        points=((34, 816), (28, 816), (28, 1272), (415, 1272), (415, 1322)),
        stroke_color=OKABE_BLUE,
        marker_id="arr-blue",
        dashed=True,
    ),
    # Col 1-2 vertical gutter (x=483) from [3F] straight onto [4D] top edge (y=1494)
    ConnectorSpec(
        edge_id="edge-multires-to-4d",
        points=((455, 1218), (455, 1275), (483, 1275), (483, 1492)),
        stroke_color=OKABE_BLUE,
        marker_id="arr-blue",
    ),
    # Tier III -> Tier IV figure & manuscript handoffs
    ConnectorSpec(
        edge_id="edge-iatlas-bench-to-figs",
        points=((704, 1218), (704, 1322)),
        stroke_color=OKABE_VERMILION,
        marker_id="arr-vermilion",
    ),
    ConnectorSpec(
        edge_id="edge-predictors-to-stab-figs",
        points=((1150, 1218), (1150, 1322)),
        stroke_color=OKABE_REDDISH_PURPLE,
        marker_id="arr-purple",
    ),
    # Right outer gutter from [2B] & [3D] down to [4C]
    ConnectorSpec(
        edge_id="edge-dossier-to-stab-figs",
        points=((1366, 838), (1372, 838), (1372, 1272), (1295, 1272), (1295, 1322)),
        stroke_color=OKABE_REDDISH_PURPLE,
        marker_id="arr-purple",
        dashed=True,
    ),
    ConnectorSpec(
        edge_id="edge-stability-to-stab-figs",
        points=((1366, 560), (1378, 560), (1378, 1278), (1335, 1278), (1335, 1322)),
        stroke_color=OKABE_REDDISH_PURPLE,
        marker_id="arr-purple",
        dashed=True,
    ),
    # Intra-Tier IV synthesis handoffs into 4E (Makefile + Typst main.pdf)
    ConnectorSpec(
        edge_id="edge-4d-to-4e",
        points=((688, 1564), (710, 1564)),
        stroke_color=OKABE_ORANGE,
        marker_id="arr-orange",
        dashed=True,
    ),
    ConnectorSpec(
        edge_id="edge-4b-to-4e",
        points=((800, 1476), (800, 1492)),
        stroke_color=OKABE_ORANGE,
        marker_id="arr-orange",
        dashed=True,
    ),
    ConnectorSpec(
        edge_id="edge-4c-to-4e",
        points=((1150, 1476), (1150, 1492)),
        stroke_color=OKABE_ORANGE,
        marker_id="arr-orange",
        dashed=True,
    ),
)

GUTTER_PILLS: Final[tuple[PillLabelSpec, ...]] = (
    PillLabelSpec(
        pill_id="pill-milo-parquet",
        cx=230,
        cy=666,
        width=338,
        text="cell_level_scores.parquet (kNN->cell logFC)",
        border_color=OKABE_BLUE,
        text_color=OKABE_BLUE,
    ),
    PillLabelSpec(
        pill_id="pill-stab-parquet",
        cx=1150,
        cy=666,
        width=248,
        text="stability {scores,paths}.parquet",
        border_color=OKABE_REDDISH_PURPLE,
        text_color=OKABE_REDDISH_PURPLE,
    ),
    PillLabelSpec(
        pill_id="pill-16-parquets",
        cx=704,
        cy=1256,
        width=248,
        text="16 Cross-Modality .parquet Tables",
        border_color=OKABE_VERMILION,
        text_color=OKABE_VERMILION,
    ),
    PillLabelSpec(
        pill_id="pill-32-biomarkers",
        cx=1150,
        cy=1256,
        width=242,
        text="32-Biomarker & HPO .parquet",
        border_color=OKABE_REDDISH_PURPLE,
        text_color=OKABE_REDDISH_PURPLE,
    ),
)


def _esc(raw: str) -> str:
    """XML-escape text content safely."""
    return html.escape(raw, quote=True)


def _render_body_line(line: str, y_pos: int, accent_color: str) -> str:
    """Render a body line, highlighting tme_datasets API declarations cleanly."""
    if line.startswith(API_PREFIX):
        rest = line[len(API_PREFIX) :]
        return (
            f'      <text x="12" y="{y_pos}" font-size="13.5px" fill="{COLOR_TEXT_PRIMARY}">'
            f'<tspan font-weight="700">tme_datasets API:</tspan>'
            f'<tspan dx="4" font-weight="600" fill="{accent_color}">{_esc(rest)}</tspan></text>'
        )
    return (
        f'      <text x="12" y="{y_pos}" font-size="13.5px" '
        f'fill="{COLOR_TEXT_SECONDARY}">{_esc(line)}</text>'
    )


def _render_card(spec: DagCardSpec) -> str:
    """Render a single sharp-cornered Nature Methods card using local (0, 0) coordinates."""
    footer_h = 28
    footer_y = spec.height - footer_h
    line_y_start = 46
    line_step = 22 if len(spec.body_lines) <= 4 else (21 if len(spec.body_lines) <= 8 else 38)
    if len(spec.body_lines) > 8:
        line_step = 38

    body_svg = "\n".join(
        _render_body_line(
            line,
            line_y_start + idx * line_step,
            spec.footer_color,
        )
        for idx, line in enumerate(spec.body_lines)
    )

    tag_w = max(36, int(len(spec.badge_tag) * 8.6) + 12)
    title_x = 12 + tag_w + 8
    title_font_px = "13.5px" if len(spec.title) > 40 or len(spec.badge_tag) > 3 else "14.5px"

    return f"""    <g id="{spec.card_id}" inkscape:label="{spec.card_id}" transform="translate({spec.x}, {spec.y})">
      <rect x="0" y="0" width="{spec.width}" height="{spec.height}" rx="0"
            fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      <rect x="0" y="0" width="{spec.width}" height="3.5" rx="0" fill="{spec.accent_color}" />
      <rect x="12" y="9" width="{tag_w}" height="20" rx="0"
            fill="{COLOR_SUBTLE_FILL}" stroke="{spec.accent_color}" stroke-width="0.75" />
      <text x="{12 + tag_w / 2:.1f}" y="23.5" font-size="13.5px" font-weight="700"
            fill="{spec.accent_color}" text-anchor="middle">{_esc(spec.badge_tag)}</text>
      <text x="{title_x}" y="24" font-size="{title_font_px}" font-weight="700"
            fill="{COLOR_TEXT_PRIMARY}">{_esc(spec.title)}</text>
{body_svg}
      <line x1="0" y1="{footer_y}" x2="{spec.width}" y2="{footer_y}"
            stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
      <rect x="1" y="{footer_y + 1}" width="{spec.width - 2}" height="{footer_h - 2}" rx="0"
            fill="{COLOR_SUBTLE_FILL}" />
      <text x="12" y="{footer_y + 19}" font-size="13.5px" font-weight="700"
            fill="{spec.footer_color}">{_esc(spec.footer_label)}</text>
    </g>"""


def _render_tier_band(band: TierBandSpec) -> str:
    """Render an outer horizontal tier frame and its header bar."""
    right_svg = (
        f'\n      <text x="{VIEWBOX_W - 28}" y="{band.y + 21}" font-size="13.5px" font-weight="600" '
        f'fill="{COLOR_TEXT_MUTED}" text-anchor="end">{_esc(band.right_summary)}</text>'
        if band.right_summary
        else ""
    )
    return f"""    <g id="{band.band_id}" inkscape:label="{band.band_id}">
      <rect x="16" y="{band.y}" width="{VIEWBOX_W - 32}" height="{band.height}" rx="0"
            fill="{COLOR_CANVAS_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      <line x1="16" y1="{band.y + 30}" x2="{VIEWBOX_W - 16}" y2="{band.y + 30}"
            stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.75" />
      <text x="28" y="{band.y + 21}" font-size="16.0px" font-weight="700"
            fill="{COLOR_TEXT_PRIMARY}">{_esc(band.letter)}</text>
      <text x="48" y="{band.y + 21}" font-size="14.5px" font-weight="700"
            fill="{COLOR_TEXT_PRIMARY}">{_esc(band.heading)}</text>{right_svg}
    </g>"""


def _render_connector(conn: ConnectorSpec) -> str:
    """Render an orthogonal SVG path with optional dash pattern and arrowhead."""
    path_cmds = " ".join(
        f"{'M' if idx == 0 else 'L'} {pt[0]} {pt[1]}" for idx, pt in enumerate(conn.points)
    )
    dash_attr = ' stroke-dasharray="5,3"' if conn.dashed else ""
    return (
        f'    <path id="{conn.edge_id}" d="{path_cmds}" fill="none" '
        f'stroke="{conn.stroke_color}" stroke-width="1.6"{dash_attr} '
        f'marker-end="url(#{conn.marker_id})" />'
    )


def _render_pill(pill: PillLabelSpec) -> str:
    """Render a white gutter badge over a connector channel."""
    half_w = pill.width / 2.0
    return f"""    <g id="{pill.pill_id}" transform="translate({pill.cx - half_w:.1f}, {pill.cy - 11})">
      <rect x="0" y="0" width="{pill.width}" height="22" rx="0"
            fill="{COLOR_CANVAS_BG}" stroke="{pill.border_color}" stroke-width="0.75" />
      <text x="{half_w:.1f}" y="15.5" font-size="13.5px" font-weight="700"
            fill="{pill.text_color}" text-anchor="middle">{_esc(pill.text)}</text>
    </g>"""


def build_pipeline_dag_svg() -> str:
    """Construct the complete multi-layer Inkscape SVG string for the provenance DAG."""
    bands_svg = "\n".join(_render_tier_band(b) for b in TIER_BANDS)
    t1_svg = "\n".join(_render_card(c) for c in TIER1_CARDS)
    t2_svg = "\n".join(_render_card(c) for c in TIER2_CARDS)
    t3_svg = "\n".join(_render_card(c) for c in TIER3_CARDS)
    t4_svg = "\n".join(_render_card(c) for c in TIER4_CARDS)
    conns_svg = "\n".join(_render_connector(e) for e in CONNECTORS)
    pills_svg = "\n".join(_render_pill(p) for p in GUTTER_PILLS)

    return f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{FIG_WIDTH_MM:.1f}mm" height="{FIG_HEIGHT_MM:.2f}mm"
     viewBox="0 0 {VIEWBOX_W} {VIEWBOX_H}"
     version="1.1" id="figure-tme-datasets-pipeline-dag"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">

  <defs>
    <marker id="arr-blue" viewBox="0 0 8 8" refX="6.5" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{OKABE_BLUE}" />
    </marker>
    <marker id="arr-teal" viewBox="0 0 8 8" refX="6.5" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{OKABE_BLUISH_GREEN}" />
    </marker>
    <marker id="arr-orange" viewBox="0 0 8 8" refX="6.5" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{OKABE_ORANGE}" />
    </marker>
    <marker id="arr-purple" viewBox="0 0 8 8" refX="6.5" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{OKABE_REDDISH_PURPLE}" />
    </marker>
    <marker id="arr-vermilion" viewBox="0 0 8 8" refX="6.5" refY="4" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 1 1 L 7 4 L 1 7 z" fill="{OKABE_VERMILION}" />
    </marker>
  </defs>

  <!-- LAYER 00: CANVAS BACKGROUND & TIER FRAMES -->
  <g inkscape:groupmode="layer" id="layer-00-background" inkscape:label="00_Canvas_Background">
    <rect id="bg-canvas" x="0" y="0" width="{VIEWBOX_W}" height="{VIEWBOX_H}" rx="0" fill="{COLOR_CANVAS_BG}" />
{bands_svg}
  </g>

  <!-- LAYER 01: HEADER & COLOR LEGEND -->
  <g inkscape:groupmode="layer" id="layer-01-header" inkscape:label="01_Figure_Header_And_Legend">
    <text id="txt-main-title" x="16" y="28" font-size="18.0px" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">
      Unified tme_datasets Package &amp; 37-Script Computational Provenance DAG
    </text>
    <g class="role-title-subtitle" id="grp-subtitle">
      <text x="16" y="50" font-size="13.5px" fill="{COLOR_TEXT_MUTED}">
        Ingestion &amp; QC live in packages/tme_datasets; 37 analytical/figure scripts reside in scripts/manuscript_pipeline/&lt;tier&gt;_*.py
      </text>
      <text x="16" y="68" font-size="13.5px" fill="{COLOR_TEXT_MUTED}">
        Cards annotate exact tme_datasets API loaders; arrows strictly denote script-to-script .parquet/.svg handoffs
      </text>
    </g>

    <g id="grp-stream-legend" transform="translate(918, 12)">
      <rect x="0" y="0" width="466" height="58" rx="0"
            fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      <rect x="12" y="11" width="12" height="12" rx="0" fill="{OKABE_BLUE}" />
      <text x="30" y="22" font-size="13.5px" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">SC Milopy &amp; Multi-Res</text>
      <rect x="234" y="11" width="12" height="12" rx="0" fill="{OKABE_BLUISH_GREEN}" />
      <text x="252" y="22" font-size="13.5px" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Cross-Modality Bulk Deconv.</text>
      <rect x="12" y="35" width="12" height="12" rx="0" fill="{OKABE_ORANGE}" />
      <text x="30" y="46" font-size="13.5px" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Pipeline QC &amp; Biomarkers</text>
      <rect x="234" y="35" width="12" height="12" rx="0" fill="{OKABE_REDDISH_PURPLE}" />
      <text x="252" y="46" font-size="13.5px" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Stability &amp; Gene Dossier</text>
    </g>
  </g>

  <!-- LAYER 02: TIER I — TME_DATASETS PACKAGE -->
  <g inkscape:groupmode="layer" id="layer-02-tier1" inkscape:label="02_Tier1_TME_Datasets_Package">
{t1_svg}
  </g>

  <!-- LAYER 03: TIER II — PRIMARY COHORT ANALYSIS SCRIPTS -->
  <g inkscape:groupmode="layer" id="layer-03-tier2" inkscape:label="03_Tier2_Primary_Analysis_Scripts">
{t2_svg}
  </g>

  <!-- LAYER 04: TIER III — DOWNSTREAM COMPUTE SCRIPTS -->
  <g inkscape:groupmode="layer" id="layer-04-tier3" inkscape:label="04_Tier3_Downstream_Compute_Scripts">
{t3_svg}
  </g>

  <!-- LAYER 05: TIER IV — PUBLICATION FIGURE & TABLE SCRIPTS -->
  <g inkscape:groupmode="layer" id="layer-05-tier4" inkscape:label="05_Tier4_Article_Figure_Scripts">
{t4_svg}
  </g>

  <!-- LAYER 06: GUTTER-CONSTRAINED CONNECTORS & ARTIFACT PILLS -->
  <g inkscape:groupmode="layer" id="layer-06-connectors" inkscape:label="06_Gutter_Connectors_And_Pills">
{conns_svg}
{pills_svg}
  </g>
</svg>
"""


def validate_dag_svg(svg_text: str) -> Result[int, str]:
    """Validate XML well-formedness, Inkscape layer count, and zero prohibited filter elements."""
    try:
        root = ET.fromstring(svg_text)
    except ET.ParseError as exc:
        return Failure(f"Malformed SVG XML: {exc}")

    if "feDropShadow" in svg_text:
        return Failure("Prohibited feDropShadow found in SVG.")

    inkscape_ns = "{http://www.inkscape.org/namespaces/inkscape}groupmode"
    layers = [elem for elem in root.iter() if elem.attrib.get(inkscape_ns) == "layer"]
    return (
        Failure("Expected at least 6 Inkscape layers.")
        if len(layers) < 6
        else Success(len(layers))
    )


def export_pipeline_dag(config: DagExportConfig) -> Result[tuple[Path, Path], str]:
    """Render, validate, and export the provenance DAG to vector SVG and 300 DPI PNG."""
    svg_text = build_pipeline_dag_svg()
    match validate_dag_svg(svg_text):
        case Failure(err):
            return Failure(err)
        case Success(_):
            pass

    try:
        for path in (config.results_svg_path, config.article_svg_path):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(svg_text, encoding="utf-8")

        target_w_px = int(round(FIG_WIDTH_MM * 300.0 / 25.4))
        target_h_px = int(round(FIG_HEIGHT_MM * 300.0 / 25.4))

        for png_p, svg_p in (
            (config.results_png_path, config.results_svg_path),
            (config.article_png_path, config.article_svg_path),
        ):
            try:
                subprocess.run(
                    [
                        "rsvg-convert",
                        "-w",
                        str(target_w_px),
                        "-h",
                        str(target_h_px),
                        "-b",
                        "white",
                        str(svg_p),
                        "-o",
                        str(png_p),
                    ],
                    check=True,
                )
            except Exception:  # noqa: BLE001
                png_bytes = vlc.svg_to_png(svg_text, scale=2.0)
                png_p.write_bytes(png_bytes)

        return Success((config.article_svg_path, config.article_png_path))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to export pipeline DAG SVG/PNG: {exc}")


def main() -> int:
    """Imperative shell CLI entrypoint."""
    outcome = export_pipeline_dag(DagExportConfig())
    if is_successful(outcome):
        svg_p, png_p = outcome.unwrap()
        sys.stdout.write(
            f"Exported tme_datasets provenance DAG: {svg_p} & {png_p} "
            f"({FIG_WIDTH_MM:.1f}mm x {FIG_HEIGHT_MM:.2f}mm)\n"
        )
        return 0
    sys.stderr.write(f"ERROR: {outcome.failure()}\n")
    return 1


if __name__ == "__main__":
    sys.exit(main())
