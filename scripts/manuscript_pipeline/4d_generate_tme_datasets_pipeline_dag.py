#!/usr/bin/env python3
"""Generate Publication Vector Diagram: tme_datasets Registry Queries & Script-to-Script Provenance DAG.

Constructs a 180.0mm double-column vector diagram adhering strictly to Nature Methods
minimalist wireframe and Inkscape multi-layer standards:
- Pure white canvas (#FFFFFF), sharp square corners (rx="0"), 0.75pt hairline borders (#CBD5E1),
  and zero drop shadows or filter glows.
- Okabe-Ito colorblind-safe accents strictly applied to functional top-bars, node ports, and arrows.
- Native Inkscape layers (inkscape:groupmode="layer") and semantic groups (<g id="...">).
- Minimum font size >= 13.5px across all labels.
- Maps all datasets and gene panels queried from `packages/tme_datasets` by Stage 1/2 scripts,
  and all intermediate `.h5ad` / `.parquet` handoffs consumed by Stage 2, 3, and 4 scripts.
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
VIEWBOX_H: Final[int] = 1188
FIG_WIDTH_MM: Final[float] = 180.0
FIG_HEIGHT_MM: Final[float] = round(FIG_WIDTH_MM * (VIEWBOX_H / VIEWBOX_W), 2)
MIN_FONT_PX: Final[float] = 13.5


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
        "output/manuscript_pipeline/4d_generate_tme_datasets_pipeline_dag/figure_tme_datasets_pipeline_dag.svg"
    )
    results_png_path: Path = Path(
        "output/manuscript_pipeline/4d_generate_tme_datasets_pipeline_dag/figure_tme_datasets_pipeline_dag.png"
    )
    article_svg_path: Path = Path("article/figures/dataset_overview/figure_tme_datasets_pipeline_dag.svg")
    article_png_path: Path = Path("article/figures/dataset_overview/figure_tme_datasets_pipeline_dag.png")


TIER_BANDS: Final[tuple[TierBandSpec, ...]] = (
    TierBandSpec(
        band_id="band-tier-1",
        letter="a",
        heading="Tier I — packages/tme_datasets Unified Registry, Providers & Gene Panels",
        right_summary="Typed DatasetSpec · Ensembl v111 · RECIST v1.1",
        y=82,
        height=194,
    ),
    TierBandSpec(
        band_id="band-tier-2",
        letter="b",
        heading="Tier II: Primary Query & .h5ad Scripts",
        right_summary="",
        y=328,
        height=190,
    ),
    TierBandSpec(
        band_id="band-tier-3",
        letter="c",
        heading="Tier III: Chained Compute Engines",
        right_summary="",
        y=570,
        height=352,
    ),
    TierBandSpec(
        band_id="band-tier-4",
        letter="d",
        heading="Tier IV: Article Figure & Table Scripts",
        right_summary="",
        y=968,
        height=202,
    ),
)

TIER1_CARDS: Final[tuple[DagCardSpec, ...]] = (
    DagCardSpec(
        card_id="card-t1-sc9",
        x=34,
        y=120,
        width=434,
        height=144,
        accent_color=OKABE_BLUE,
        badge_tag="SC-9",
        title="9 Single-Cell ICB Cohorts (3.60M Cells)",
        body_lines=(
            "GSE120575 (Melanoma), CELLxGENE_7b20c613 (Pan-Ca),",
            "05a8c945 (CRC), 6f9de485 (TNBC), GSE207422 / 243013 /",
            "233203 (NSCLC), GSE200996 (HNSCC), GSE316195 (PDAC)",
        ),
        footer_label="API: load_dataset() · resolve_milopy_cohorts() · load_backed()",
        footer_color=OKABE_BLUE,
    ),
    DagCardSpec(
        card_id="card-t1-bulk11",
        x=488,
        y=120,
        width=428,
        height=144,
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
        card_id="card-t1-panels",
        x=936,
        y=120,
        width=430,
        height=144,
        accent_color=OKABE_REDDISH_PURPLE,
        badge_tag="PANELS",
        title="SC Reference Pool, TCGA & IO Gene Sets",
        body_lines=(
            "17+ SC Reference Atlases (GSE179994, Maynard_NSCLC, ...)",
            "5 TCGA Primary Baselines (SKCM, BLCA, KIRC, BRCA, PAAD)",
            "IMMUNOTHERAPY_GENE_PANEL (101) · Checkpoint · Bagaev",
        ),
        footer_label="API: load_geneset_collection() · list_registered_datasets()",
        footer_color=OKABE_REDDISH_PURPLE,
    ),
)

TIER2_CARDS: Final[tuple[DagCardSpec, ...]] = (
    DagCardSpec(
        card_id="card-t2-milopy",
        x=34,
        y=364,
        width=434,
        height=142,
        accent_color=OKABE_BLUE,
        badge_tag="2A",
        title="2a1_sample_and_analyze_cohorts.py",
        body_lines=(
            "+ 2a2_combine_cohorts_and_analyze.py (manuscript_pipeline/)",
            "Queries SC-9 via resolve_milopy_cohorts() & Harmony",
            "Fits milopy kNN Negative Binomial GLM response DA",
        ),
        footer_label="Out: .../2a1_sample_and_analyze_cohorts/*.parquet",
        footer_color=OKABE_BLUE,
    ),
    DagCardSpec(
        card_id="card-t2-ingest",
        x=488,
        y=364,
        width=428,
        height=142,
        accent_color=OKABE_ORANGE,
        badge_tag="1A–1C",
        title="1a_download_and_preprocess_all.py",
        body_lines=(
            "+ 1b_batch_preprocess_cohorts.py & 1c_run_qc_*.py",
            "Queries SC-9 & BULK-15 via load_dataset() & Ensembl v111",
            "Writes atomic QC-validated .h5ad + summary parquets",
        ),
        footer_label="Out: data/preprocessed/*.h5ad & .../1a_*/ .../1c_*/",
        footer_color="#B45309",
    ),
    DagCardSpec(
        card_id="card-t2-stability",
        x=936,
        y=364,
        width=430,
        height=142,
        accent_color=OKABE_REDDISH_PURPLE,
        badge_tag="2B",
        title="2b_run_stability_selection_iatlas.py",
        body_lines=(
            "Queries load_iatlas_cohort_or_combined() + IO Panels",
            "12 cohorts × 13 fitters (Lasso, LassoFWL, MERF, SLOPE)",
            "Complementary-pairs subsampling under PFER q_budget",
        ),
        footer_label="Out: .../2b_run_stability_selection_iatlas/*.parquet",
        footer_color=OKABE_REDDISH_PURPLE,
    ),
)

TIER3_CARDS: Final[tuple[DagCardSpec, ...]] = (
    DagCardSpec(
        card_id="card-t3-pseudobulk",
        x=34,
        y=606,
        width=434,
        height=142,
        accent_color=OKABE_BLUE,
        badge_tag="3A",
        title="3a_sade_feldman_pseudobulk_deconv_benchmark.py",
        body_lines=(
            "Reads 9 SC .h5ad [1A] + cell_level_scores [2A1]",
            "Builds 355 Leiden r=2.0 reference centroids (Phi_1..Phi_9)",
            "Evaluates 5-sweep Milo vs. pseudobulk β concordance",
        ),
        footer_label="Out: .../3a_sade_feldman_pseudobulk_deconv_benchmark/",
        footer_color=OKABE_BLUE,
    ),
    DagCardSpec(
        card_id="card-t3-unified",
        x=34,
        y=768,
        width=434,
        height=142,
        accent_color=OKABE_SKY_BLUE,
        badge_tag="3B",
        title="3b_pan_cohort_unified_deconv_benchmark.py",
        body_lines=(
            "Reads 355 cohort centroids & sweeps from [3A]",
            "Ward-clusters 355 SC states -> 42 Consensus-Meta42",
            "Verifies Unified-All9, Meta42 & LOCO concordance",
        ),
        footer_label="Out: .../3b_pan_cohort_unified_deconv_benchmark/",
        footer_color=OKABE_BLUE,
    ),
    DagCardSpec(
        card_id="card-t3-iatlas-bench",
        x=488,
        y=606,
        width=428,
        height=304,
        accent_color=OKABE_BLUISH_GREEN,
        badge_tag="3C",
        title="3c_iatlas_instaprism_milopy_benchmark.py",
        body_lines=(
            "Inputs: 11 Bulk .h5ad [1A] + 355 SC & Meta-42",
            "  reference/milopy parquets (from [3A] & [3B])",
            "1. Detects log2 bulk scale, linearizes (2^X - 1 -> 10^6 TPM)",
            "2. Runs InstaPrism EM across 13 Refs (Tiers 1–3) × 15 Cohorts",
            "3. Evaluates 5 ML families (Ridge, ENet, Lasso, OLS, RF) + LOCO",
            "4. DerSimonian–Laird Pan-iAtlas-Meta β vs. milopy logFC",
            "5. Extracts 168 top-4 HGNC DEGs (14-state & All-42 atlases)",
            "6. Computes 8×8 & 42×42 epistasis z, RF 2D H & 3×3 regimes",
        ),
        footer_label="Out: .../3c_iatlas_instaprism_milopy_benchmark/ (11)",
        footer_color=OKABE_BLUISH_GREEN,
    ),
    DagCardSpec(
        card_id="card-t3-dossier",
        x=936,
        y=606,
        width=430,
        height=142,
        accent_color=OKABE_REDDISH_PURPLE,
        badge_tag="3D",
        title="3d_curate_expert_gene_list.py",
        body_lines=(
            "Reads 2b_run_stability_selection_iatlas/*.parquet [2B]",
            "+ queries pan-cancer bulk fold-changes via tme_datasets",
            "Computes composite Evidence Index E(g) across 91 runs",
        ),
        footer_label="Out: .../3d_curate_expert_gene_list/ & table_*.typ",
        footer_color=OKABE_REDDISH_PURPLE,
    ),
    DagCardSpec(
        card_id="card-t3-predictors",
        x=936,
        y=768,
        width=430,
        height=142,
        accent_color=OKABE_ORANGE,
        badge_tag="3E",
        title="3e1_benchmark_iatlas_icb_predictors.py",
        body_lines=(
            "+ 3e2_run_hpo_pipeline.py (Optuna TPE + LOCO)",
            "Queries 9 iAtlas & SC cohorts via tme_datasets.load_dataset()",
            "Benchmarks 27 RNA/multi-omic signatures & HPO sampling",
        ),
        footer_label="Out: .../3e1_benchmark_*/ & .../3e2_run_hpo_*/",
        footer_color="#B45309",
    ),
)

TIER4_CARDS: Final[tuple[DagCardSpec, ...]] = (
    DagCardSpec(
        card_id="card-t4-sc-figs",
        x=34,
        y=1004,
        width=434,
        height=152,
        accent_color=OKABE_BLUE,
        badge_tag="4A",
        title="4a1–4a7 Single-Cell & 2×2 Master Figures",
        body_lines=(
            "4a1_visualize_cohort_umaps · 4a2_generate_combined_*",
            "4a3..4a6_plot_cross_cohort_{umaps,da,deconv_*}",
            "4a7_assemble_immunotherapy_2x2_figure (reads 2A1 & 3A)",
        ),
        footer_label="Out: .../4a1_*/..4a7_*/ -> article/figures/*",
        footer_color=OKABE_BLUE,
    ),
    DagCardSpec(
        card_id="card-t4-iatlas-figs",
        x=488,
        y=1004,
        width=428,
        height=152,
        accent_color=OKABE_VERMILION,
        badge_tag="4B",
        title="4b_plot_iatlas_instaprism_milopy_figures.py",
        body_lines=(
            "Consumes all 11 parquets from [3C] -> renders 11 Figures:",
            "Focused (8): ROC/PR-AUC, Milo Grid, 14-State Forest/DEG, 8×8",
            "Full Atlas (3): All-42 Forest, 168-Gene DEG, 42×42 Epistasis",
        ),
        footer_label="Out: .../4b_plot_iatlas_*/ -> figures/milo-deconv/",
        footer_color=OKABE_VERMILION,
    ),
    DagCardSpec(
        card_id="card-t4-stab-figs",
        x=936,
        y=1004,
        width=430,
        height=152,
        accent_color=OKABE_REDDISH_PURPLE,
        badge_tag="4C–4D",
        title="4c1–4c5 Stability/HPO & 4d Provenance DAG",
        body_lines=(
            "4c1..4c4_*{stability_paths,cohort_*,fitter_*}_overlap",
            "4c5_generate_hpo_summary_figures · 4d_generate_*_dag",
            "article/Makefile: rsvg-convert (300 DPI) + typst compile",
        ),
        footer_label="Out: .../4c1_*/..4d_*/ -> article/main.pdf",
        footer_color=OKABE_REDDISH_PURPLE,
    ),
)

# Orthogonal connectors strictly inside gutters and clear of all tier header text (x > 345).
CONNECTORS: Final[tuple[ConnectorSpec, ...]] = (
    # Tier I -> Tier II queries
    ConnectorSpec(
        edge_id="edge-sc9-to-milopy",
        points=((210, 264), (210, 316), (365, 316), (365, 362)),
        stroke_color=OKABE_BLUE,
        marker_id="arr-blue",
    ),
    ConnectorSpec(
        edge_id="edge-sc9-to-ingest",
        points=((442, 264), (442, 316), (575, 316), (575, 362)),
        stroke_color=OKABE_ORANGE,
        marker_id="arr-orange",
    ),
    ConnectorSpec(
        edge_id="edge-bulk11-to-ingest",
        points=((702, 264), (702, 362)),
        stroke_color=OKABE_ORANGE,
        marker_id="arr-orange",
    ),
    ConnectorSpec(
        edge_id="edge-bulk11-to-stability",
        points=((855, 264), (855, 316), (1045, 316), (1045, 362)),
        stroke_color=OKABE_REDDISH_PURPLE,
        marker_id="arr-purple",
    ),
    ConnectorSpec(
        edge_id="edge-panels-to-stability",
        points=((1210, 264), (1210, 362)),
        stroke_color=OKABE_REDDISH_PURPLE,
        marker_id="arr-purple",
    ),
    # Tier II -> Tier III artifact handoffs
    ConnectorSpec(
        edge_id="edge-milopy-to-pseudobulk",
        points=((210, 506), (210, 558), (365, 558), (365, 604)),
        stroke_color=OKABE_BLUE,
        marker_id="arr-blue",
    ),
    ConnectorSpec(
        edge_id="edge-ingest-sc-to-pseudobulk",
        points=((565, 506), (565, 534), (442, 534), (442, 604)),
        stroke_color=OKABE_ORANGE,
        marker_id="arr-orange",
    ),
    ConnectorSpec(
        edge_id="edge-ingest-bulk-to-iatlas-bench",
        points=((702, 506), (702, 604)),
        stroke_color=OKABE_BLUISH_GREEN,
        marker_id="arr-teal",
    ),
    ConnectorSpec(
        edge_id="edge-ingest-to-predictors",
        points=((855, 506), (855, 534), (926, 534), (926, 838), (934, 838)),
        stroke_color=OKABE_ORANGE,
        marker_id="arr-orange",
        dashed=True,
    ),
    ConnectorSpec(
        edge_id="edge-stability-to-dossier",
        points=((1160, 506), (1160, 604)),
        stroke_color=OKABE_REDDISH_PURPLE,
        marker_id="arr-purple",
    ),
    # Intra-Tier III script-to-script handoffs
    ConnectorSpec(
        edge_id="edge-pseudobulk-to-unified",
        points=((250, 748), (250, 766)),
        stroke_color=OKABE_BLUE,
        marker_id="arr-blue",
    ),
    ConnectorSpec(
        edge_id="edge-pseudobulk-to-iatlas-bench",
        points=((468, 678), (486, 678)),
        stroke_color=OKABE_BLUE,
        marker_id="arr-blue",
    ),
    ConnectorSpec(
        edge_id="edge-unified-to-iatlas-bench",
        points=((468, 838), (486, 838)),
        stroke_color=OKABE_BLUE,
        marker_id="arr-blue",
    ),
    # Left outer gutter from 2A/3A down to 4A (enters at x=365 > 335, clear of Tier IV title)
    ConnectorSpec(
        edge_id="edge-pseudobulk-to-sc-figs",
        points=((34, 678), (25, 678), (25, 944), (365, 944), (365, 1002)),
        stroke_color=OKABE_BLUE,
        marker_id="arr-blue",
        dashed=True,
    ),
    # Tier III -> Tier IV figure & manuscript handoffs
    ConnectorSpec(
        edge_id="edge-iatlas-bench-to-figs",
        points=((702, 910), (702, 1002)),
        stroke_color=OKABE_VERMILION,
        marker_id="arr-vermilion",
    ),
    ConnectorSpec(
        edge_id="edge-predictors-to-stab-figs",
        points=((1150, 910), (1150, 1002)),
        stroke_color=OKABE_REDDISH_PURPLE,
        marker_id="arr-purple",
    ),
    ConnectorSpec(
        edge_id="edge-dossier-to-stab-figs",
        points=((1366, 678), (1375, 678), (1375, 944), (1270, 944), (1270, 1002)),
        stroke_color=OKABE_REDDISH_PURPLE,
        marker_id="arr-purple",
        dashed=True,
    ),
)

GUTTER_PILLS: Final[tuple[PillLabelSpec, ...]] = (
    PillLabelSpec(
        pill_id="pill-q-sc",
        cx=210,
        cy=294,
        width=194,
        text="resolve_milopy_cohorts()",
        border_color=OKABE_BLUE,
        text_color=OKABE_BLUE,
    ),
    PillLabelSpec(
        pill_id="pill-q-h5ad",
        cx=702,
        cy=294,
        width=188,
        text="load_dataset() [SC + Bulk]",
        border_color=OKABE_ORANGE,
        text_color="#B45309",
    ),
    PillLabelSpec(
        pill_id="pill-q-stab",
        cx=1210,
        cy=294,
        width=220,
        text="iAtlas Pools + 101-Gene Panel",
        border_color=OKABE_REDDISH_PURPLE,
        text_color=OKABE_REDDISH_PURPLE,
    ),
    PillLabelSpec(
        pill_id="pill-milo-parquet",
        cx=210,
        cy=534,
        width=206,
        text="cell_level_scores.parquet",
        border_color=OKABE_BLUE,
        text_color=OKABE_BLUE,
    ),
    PillLabelSpec(
        pill_id="pill-sc-h5ad",
        cx=504,
        cy=534,
        width=106,
        text="9 SC .h5ad",
        border_color=OKABE_ORANGE,
        text_color="#B45309",
    ),
    PillLabelSpec(
        pill_id="pill-bulk-h5ad",
        cx=702,
        cy=534,
        width=136,
        text="11 Bulk .h5ad",
        border_color=OKABE_BLUISH_GREEN,
        text_color=OKABE_BLUISH_GREEN,
    ),
    PillLabelSpec(
        pill_id="pill-stab-parquet",
        cx=1160,
        cy=534,
        width=232,
        text="stability {scores,paths}.parquet",
        border_color=OKABE_REDDISH_PURPLE,
        text_color=OKABE_REDDISH_PURPLE,
    ),
    PillLabelSpec(
        pill_id="pill-11-parquets",
        cx=702,
        cy=942,
        width=244,
        text="11 Cross-Modality .parquet Tables",
        border_color=OKABE_VERMILION,
        text_color=OKABE_VERMILION,
    ),
)


def _esc(raw: str) -> str:
    """XML-escape text content safely."""
    return html.escape(raw, quote=True)


def _render_card(spec: DagCardSpec) -> str:
    """Render a single sharp-cornered Nature Methods card using local (0, 0) coordinates."""
    footer_h = 28
    footer_y = spec.height - footer_h
    line_y_start = 46
    line_step = 22 if len(spec.body_lines) <= 4 else 21

    body_svg = "\n".join(
        f'      <text x="12" y="{line_y_start + idx * line_step}" '
        f'font-size="13.5px" fill="{COLOR_TEXT_SECONDARY}">{_esc(line)}</text>'
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
      Unified tme_datasets Registry Queries &amp; Multi-Stage Script Provenance DAG
    </text>
    <g class="role-title-subtitle" id="grp-subtitle">
      <text x="16" y="50" font-size="13.5px" fill="{COLOR_TEXT_MUTED}">
        All 26 pipeline scripts reside in scripts/manuscript_pipeline/&lt;tier&gt;_*.py and query packages/tme_datasets
      </text>
      <text x="16" y="68" font-size="13.5px" fill="{COLOR_TEXT_MUTED}">
        Each script writes to output/manuscript_pipeline/&lt;script_stem&gt;/ (solid = primary .h5ad/.parquet; dashed = aux)
      </text>
    </g>

    <g id="grp-stream-legend" transform="translate(918, 12)">
      <rect x="0" y="0" width="466" height="58" rx="0"
            fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.75" />
      <rect x="12" y="11" width="12" height="12" rx="0" fill="{OKABE_BLUE}" />
      <text x="30" y="22" font-size="13.5px" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">SC Milopy &amp; Pseudobulk</text>
      <rect x="234" y="11" width="12" height="12" rx="0" fill="{OKABE_BLUISH_GREEN}" />
      <text x="252" y="22" font-size="13.5px" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Cross-Modality Bulk Deconv.</text>
      <rect x="12" y="35" width="12" height="12" rx="0" fill="{OKABE_ORANGE}" />
      <text x="30" y="46" font-size="13.5px" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">.h5ad Cache &amp; Signatures</text>
      <rect x="234" y="35" width="12" height="12" rx="0" fill="{OKABE_REDDISH_PURPLE}" />
      <text x="252" y="46" font-size="13.5px" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Stability Selection &amp; Dossier</text>
    </g>
  </g>

  <!-- LAYER 02: TIER I — TME_DATASETS REGISTRY -->
  <g inkscape:groupmode="layer" id="layer-02-tier1" inkscape:label="02_Tier1_TME_Datasets_Registry">
{t1_svg}
  </g>

  <!-- LAYER 03: TIER II — INGESTION & PRIMARY QUERY SCRIPTS -->
  <g inkscape:groupmode="layer" id="layer-03-tier2" inkscape:label="03_Tier2_Primary_Query_Scripts">
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
