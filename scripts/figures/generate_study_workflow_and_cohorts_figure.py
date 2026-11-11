#!/usr/bin/env python3
"""Generate the 5-panel Study Workflow & Clinical/Single-Cell Cohorts Figure (SVG + 300 DPI PNG).

Constructs a publication-grade 3-row, 5-panel vector figure (180mm x 230mm, viewBox="0 0 1080 1380")
adhering strictly to Nature Methods minimalist wireframe standards, Okabe-Ito colorblind-safe
scientific palette, and Inkscape multi-layer architecture:
  - Panel (a): Tumor Microenvironment Ecosystem (embeds assets/tme_tumor_microenvironment_ecosystem.svg
               and cell icons from assets/cell_*.svg across 4 annotated compartments).
  - Panel (b): Specimen Processing & High-Throughput Sequencing Workflow (Resection -> Pathological
               Grossing & Spatial Sampling -> Tissue Preservation -> RNA Extraction & Library Prep ->
               High-Throughput Sequencing).
  - Panel (c): Single-Cell & Bulk Deconvolution Classification Pipeline (Single Cell Transcriptomic
               Datasets -> Preprocessing -> Cell State Reference Construction [using assets/cell_*.svg]
               + Patient Cohort -> Bulk RNA-sequencing data -> Deconvolution -> Driver Cell States
               Identification -> Response Classification).
  - Panel (d): Bulk RNA-seq ICB Clinical Cohorts (n = 1,097 patients across 9 trials, RECIST v1.1 R/NR
               stacked bars, ORR%, Overall summary row, and RECIST response legend).
  - Panel (e): Single-Cell ICB Response Cohorts: Milo Compendium (9 Cohorts, n = 613 evaluations,
               3,597,092 cells with updated manuscript study citations and Milo Total summary row).
"""

from __future__ import annotations

import html
from pathlib import Path
import re
import sys
from typing import Final, Literal, Sequence, assert_never
import xml.etree.ElementTree as ET

from pydantic import BaseModel, ConfigDict
from returns.pipeline import is_successful
from returns.result import Failure, Result, Success

# Ensure repository root is importable
REPO_ROOT: Final[Path] = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

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

# Physical dimensions and coordinate grid (6.0 px per mm -> 1080 x 1216)
SVG_WIDTH_MM: Final[str] = "180mm"
SVG_HEIGHT_MM: Final[str] = "203mm"
VIEWBOX_W: Final[int] = 1080
VIEWBOX_H: Final[int] = 1216
MIN_FONT_SIZE_PX: Final[float] = 9.5

COLOR_SUMMARY_FILL: Final[str] = "#F1F5F9"
COLOR_AXIS_LINE: Final[str] = "#1E293B"


# ==============================================================================
# IMMUTABLE DATA MODELS (FROZEN PYDANTIC, NO NONE)
# ==============================================================================
class TmeCompartmentSpec(BaseModel):
    """Immutable specification for a TME compartment callout card in Panel (a)."""

    model_config = ConfigDict(frozen=True)

    badge_num: str
    title: str
    subtitle: str
    accent_color: str
    card_y: float
    pin_x: float
    pin_y: float
    asset_key: str
    asset_scale: float


WetLabStepId = Literal["resection", "grossing", "preservation", "rna_prep", "sequencing"]


class WetLabStepSpec(BaseModel):
    """Immutable specification for a clinical specimen processing step in Panel (b)."""

    model_config = ConfigDict(frozen=True)

    step_id: WetLabStepId
    lines: tuple[str, ...]
    sub_line: str
    center_y: float


class BulkTrialCohortSpec(BaseModel):
    """Immutable model for a bulk RNA-seq ICB clinical trial cohort in Panel (d)."""

    model_config = ConfigDict(frozen=True)

    study: str
    cancer_type: str
    cancer_color: str
    regimen: str
    n_total: int
    n_resp: int
    n_non_resp: int
    orr_pct: float


class SingleCellMiloCohortSpec(BaseModel):
    """Immutable model for a single-cell ICB cohort analyzed with Milo in Panel (e)."""

    model_config = ConfigDict(frozen=True)

    publication: str
    cancer_type: str
    cancer_color: str
    regimen: str
    cells_raw: int
    n_resp: int
    n_non_resp: int
    orr_pct: float


class AssetBundle(BaseModel):
    """Immutable bundle of extracted SVG defs and body strings from assets/."""

    model_config = ConfigDict(frozen=True)

    defs_svg: str
    tme_body: str
    cell_cancer: str
    cell_t_lymph: str
    cell_macro_m2: str
    cell_caf: str


# ==============================================================================
# DECLARATIVE PANEL SPECIFICATIONS
# ==============================================================================
TME_COMPARTMENTS: Final[tuple[TmeCompartmentSpec, ...]] = (
    TmeCompartmentSpec(
        badge_num="1",
        title="1. Malignant Cells",
        subtitle="Proliferating & hypoxic tumor",
        accent_color="#C2410C",
        card_y=62.0,
        pin_x=126.9,
        pin_y=145.6,
        asset_key="cell_cancer",
        asset_scale=0.115,
    ),
    TmeCompartmentSpec(
        badge_num="2",
        title="2. Effector Immune",
        subtitle="CD8+ T, NK & B lymphocytes",
        accent_color=OKABE_BLUE,
        card_y=129.0,
        pin_x=170.3,
        pin_y=187.9,
        asset_key="cell_t_lymph",
        asset_scale=0.13,
    ),
    TmeCompartmentSpec(
        badge_num="3",
        title="3. Suppressive Immune",
        subtitle="Tregs, M2 TAMs & MDSCs",
        accent_color="#6C32A6",
        card_y=196.0,
        pin_x=149.8,
        pin_y=237.2,
        asset_key="cell_macro_m2",
        asset_scale=0.115,
    ),
    TmeCompartmentSpec(
        badge_num="4",
        title="4. Stroma & Vessels",
        subtitle="CAFs, collagen ECM & vessels",
        accent_color="#8C6239",
        card_y=263.0,
        pin_x=110.9,
        pin_y=302.2,
        asset_key="cell_caf",
        asset_scale=0.098,
    ),
)

WETLAB_STEPS: Final[tuple[WetLabStepSpec, ...]] = (
    WetLabStepSpec(
        step_id="resection",
        lines=("Resection",),
        sub_line="",
        center_y=72.0,
    ),
    WetLabStepSpec(
        step_id="grossing",
        lines=("Pathological Grossing", "& Spatial Sampling"),
        sub_line="",
        center_y=136.0,
    ),
    WetLabStepSpec(
        step_id="preservation",
        lines=("Tissue Preservation",),
        sub_line="(Snap-Frozen)",
        center_y=200.0,
    ),
    WetLabStepSpec(
        step_id="rna_prep",
        lines=("RNA Extraction &", "Library Preparation"),
        sub_line="",
        center_y=264.0,
    ),
    WetLabStepSpec(
        step_id="sequencing",
        lines=("High-Throughput", "Sequencing"),
        sub_line="",
        center_y=328.0,
    ),
)

BULK_TRIAL_COHORTS: Final[tuple[BulkTrialCohortSpec, ...]] = (
    BulkTrialCohortSpec(
        study="Hugo (2016)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        regimen="Pembrolizumab (PD-1)",
        n_total=27,
        n_resp=14,
        n_non_resp=13,
        orr_pct=51.9,
    ),
    BulkTrialCohortSpec(
        study="Riaz (2017)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        regimen="Nivolumab (PD-1)",
        n_total=107,
        n_resp=34,
        n_non_resp=73,
        orr_pct=31.8,
    ),
    BulkTrialCohortSpec(
        study="Liu (2019)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        regimen="Nivolumab / Pembro (PD-1)",
        n_total=122,
        n_resp=48,
        n_non_resp=74,
        orr_pct=39.3,
    ),
    BulkTrialCohortSpec(
        study="Gide (2019)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        regimen="Anti-PD-1 ± Ipi (Combo)",
        n_total=91,
        n_resp=63,
        n_non_resp=28,
        orr_pct=69.2,
    ),
    BulkTrialCohortSpec(
        study="Rosenberg (2016)",
        cancer_type="Bladder",
        cancer_color=OKABE_SKY_BLUE,
        regimen="Atezolizumab (PD-L1)",
        n_total=347,
        n_resp=68,
        n_non_resp=279,
        orr_pct=19.6,
    ),
    BulkTrialCohortSpec(
        study="McDermott (2018)",
        cancer_type="ccRCC",
        cancer_color=OKABE_BLUISH_GREEN,
        regimen="Atezo ± Bevacizumab",
        n_total=263,
        n_resp=97,
        n_non_resp=166,
        orr_pct=36.9,
    ),
    BulkTrialCohortSpec(
        study="Choueiri (2016)",
        cancer_type="ccRCC",
        cancer_color=OKABE_BLUISH_GREEN,
        regimen="Nivolumab (PD-1)",
        n_total=16,
        n_resp=6,
        n_non_resp=10,
        orr_pct=37.5,
    ),
    BulkTrialCohortSpec(
        study="Padron (2022)",
        cancer_type="PDAC",
        cancer_color=OKABE_BLUE,
        regimen="Nivo + Sotigalimab",
        n_total=93,
        n_resp=16,
        n_non_resp=77,
        orr_pct=17.2,
    ),
    BulkTrialCohortSpec(
        study="Anders (2021)",
        cancer_type="TNBC",
        cancer_color=OKABE_REDDISH_PURPLE,
        regimen="Atezo + Nab-paclitaxel",
        n_total=31,
        n_resp=11,
        n_non_resp=20,
        orr_pct=35.5,
    ),
)

SC_MILO_COHORTS: Final[tuple[SingleCellMiloCohortSpec, ...]] = (
    SingleCellMiloCohortSpec(
        publication="Sade-Feldman (2018)",
        cancer_type="Melanoma",
        cancer_color=OKABE_ORANGE,
        regimen="Anti-PD-1 ± CTLA-4",
        cells_raw=16291,
        n_resp=14,
        n_non_resp=25,
        orr_pct=35.9,
    ),
    SingleCellMiloCohortSpec(
        publication="Gondal (2025)",
        cancer_type="Pan-Cancer",
        cancer_color=OKABE_ORANGE,
        regimen="Nivolumab / Pembro",
        cells_raw=332972,
        n_resp=41,
        n_non_resp=82,
        orr_pct=33.3,
    ),
    SingleCellMiloCohortSpec(
        publication="Marteau (2026)",
        cancer_type="Colorectal",
        cancer_color=OKABE_BLUISH_GREEN,
        regimen="Pembro / Nivolumab",
        cells_raw=1107923,
        n_resp=48,
        n_non_resp=29,
        orr_pct=62.3,
    ),
    SingleCellMiloCohortSpec(
        publication="Yan (2026)",
        cancer_type="Breast",
        cancer_color=OKABE_REDDISH_PURPLE,
        regimen="Neoadjuvant AC Chemo",
        cells_raw=424114,
        n_resp=45,
        n_non_resp=37,
        orr_pct=54.9,
    ),
    SingleCellMiloCohortSpec(
        publication="Hu (2023)",
        cancer_type="NSCLC",
        cancer_color=OKABE_SKY_BLUE,
        regimen="Anti-PD-1 + Chemo",
        cells_raw=91997,
        n_resp=4,
        n_non_resp=10,
        orr_pct=28.6,
    ),
    SingleCellMiloCohortSpec(
        publication="Liu (2025)",
        cancer_type="NSCLC",
        cancer_color=OKABE_SKY_BLUE,
        regimen="PD-(L)1 + Chemo",
        cells_raw=1254749,
        n_resp=130,
        n_non_resp=112,
        orr_pct=53.7,
    ),
    SingleCellMiloCohortSpec(
        publication="Hong (2025)",
        cancer_type="NSCLC",
        cancer_color=OKABE_SKY_BLUE,
        regimen="Atezo + Bev + Chemo",
        cells_raw=40142,
        n_resp=3,
        n_non_resp=4,
        orr_pct=42.9,
    ),
    SingleCellMiloCohortSpec(
        publication="Luoma (2022)",
        cancer_type="HNSCC",
        cancer_color=OKABE_VERMILION,
        regimen="Nivolumab ± Ipi",
        cells_raw=245253,
        n_resp=4,
        n_non_resp=8,
        orr_pct=33.3,
    ),
    SingleCellMiloCohortSpec(
        publication="Raufi (2026)",
        cancer_type="Pancreatic",
        cancer_color=OKABE_BLUE,
        regimen="Cemiplimab + CXCR4i",
        cells_raw=83651,
        n_resp=15,
        n_non_resp=2,
        orr_pct=88.2,
    ),
)


# ==============================================================================
# PURE SVG ASSET EXTRACTION & NAMESPACING
# ==============================================================================
def extract_svg_asset(
    asset_path: Path,
    prefix: str,
    strip_body_ids: bool = False,
) -> Result[tuple[str, str], str]:
    """Extract and namespace <defs> children and artwork <g> elements from an SVG asset."""
    if not asset_path.exists():
        return Failure(f"SVG asset not found: {asset_path}")

    raw_svg = asset_path.read_text(encoding="utf-8")
    ids_found = set(re.findall(r'\bid="([^"]+)"', raw_svg))
    namespaced_svg = raw_svg
    for old_id in sorted(ids_found, key=len, reverse=True):
        new_id = f"{prefix}{old_id}"
        namespaced_svg = namespaced_svg.replace(f'id="{old_id}"', f'id="{new_id}"')
        namespaced_svg = namespaced_svg.replace(f"url(#{old_id})", f"url(#{new_id})")
        namespaced_svg = namespaced_svg.replace(f'href="#{old_id}"', f'href="#{new_id}"')

    defs_match = re.search(r"<defs[^>]*>(.*?)</defs>", namespaced_svg, flags=re.DOTALL)
    defs_inner = defs_match.group(1).strip() if defs_match else ""

    body_match = re.search(r"</defs>(.*)</svg>", namespaced_svg, flags=re.DOTALL)
    if not body_match:
        body_match = re.search(r"<svg[^>]*>(.*)</svg>", namespaced_svg, flags=re.DOTALL)
    if not body_match:
        return Failure(f"Could not extract SVG body from {asset_path}")

    body_content = body_match.group(1)
    body_clean = re.sub(
        r"<!-- 00_Canvas_Background -->.*?<g[^>]*id=\"[^\"]*layer-00-background\"[^>]*>.*?</g>",
        "",
        body_content,
        flags=re.DOTALL,
    )
    body_clean = body_clean.replace('inkscape:groupmode="layer"', "")
    if strip_body_ids:
        body_clean = re.sub(r'\s+id="[^"]+"', "", body_clean)
    return Success((defs_inner, body_clean.strip()))


def load_asset_bundle(assets_dir: Path) -> Result[AssetBundle, str]:
    """Load and namespace all required SVG assets from the assets/ directory."""
    asset_specs: tuple[tuple[str, str, str, bool], ...] = (
        ("tme", "tme_tumor_microenvironment_ecosystem.svg", "ast_tme_", False),
        ("cancer", "cell_cancer_cell.svg", "ast_canc_", True),
        ("tcell", "cell_t_lymphocyte.svg", "ast_tcell_", True),
        ("macro", "cell_macrophage_m2.svg", "ast_mac_", True),
        ("caf", "cell_fibroblast_caf.svg", "ast_caf_", True),
    )

    defs_parts: list[str] = []
    bodies: dict[str, str] = {}

    for key, filename, prefix, strip_ids in asset_specs:
        res = extract_svg_asset(assets_dir / filename, prefix, strip_body_ids=strip_ids)
        if not is_successful(res):
            return Failure(res.failure())
        d_str, b_str = res.unwrap()
        if d_str:
            defs_parts.append(d_str)
        bodies[key] = b_str

    return Success(
        AssetBundle(
            defs_svg="\n".join(defs_parts),
            tme_body=bodies["tme"],
            cell_cancer=bodies["cancer"],
            cell_t_lymph=bodies["tcell"],
            cell_macro_m2=bodies["macro"],
            cell_caf=bodies["caf"],
        )
    )


# ==============================================================================
# PANEL (A): TUMOR MICROENVIRONMENT ECOSYSTEM RENDERER
# ==============================================================================
def render_panel_a_tme_ecosystem(bundle: AssetBundle) -> str:
    """Render Panel (a): Heterogeneous TME illustration and 4 compartment callout cards."""
    asset_lookup: dict[str, str] = {
        "cell_cancer": bundle.cell_cancer,
        "cell_t_lymph": bundle.cell_t_lymph,
        "cell_macro_m2": bundle.cell_macro_m2,
        "cell_caf": bundle.cell_caf,
    }

    cards_svg: list[str] = []
    pins_svg: list[str] = []

    for spec in TME_COMPARTMENTS:
        cell_art = asset_lookup[spec.asset_key]
        card_mid_y = spec.card_y + 29.5
        cards_svg.append(f"""
      <g id="pa-compartment-card-{spec.badge_num}" transform="translate(254, {spec.card_y})">
        <rect x="0" y="0" width="198" height="59.0" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />
        <rect x="0" y="0" width="4.5" height="59.0" fill="{spec.accent_color}" />
        <g id="pa-swatch-{spec.badge_num}" transform="translate(8.0, 14.0) scale({spec.asset_scale})">
          {cell_art}
        </g>
        <text x="46.0" y="24.5" font-size="12.2" font-weight="700" fill="{spec.accent_color}">{html.escape(spec.title)}</text>
        <text x="46.0" y="43.5" font-size="10.4" fill="{COLOR_TEXT_SECONDARY}">{html.escape(spec.subtitle)}</text>
      </g>""")

        pins_svg.append(f"""
      <g id="pa-callout-pin-{spec.badge_num}">
        <path d="M {spec.pin_x} {spec.pin_y} L 243 {card_mid_y} L 254 {card_mid_y}"
              stroke="{spec.accent_color}" stroke-width="1.2" stroke-dasharray="3.0,2.0" fill="none" />
        <circle cx="{spec.pin_x}" cy="{spec.pin_y}" r="9.5" fill="{spec.accent_color}" stroke="#FFFFFF" stroke-width="1.2" />
        <text x="{spec.pin_x}" y="{spec.pin_y + 3.8}" font-size="10.8" font-weight="700" fill="#FFFFFF" text-anchor="middle">{spec.badge_num}</text>
      </g>""")

    return f"""
  <!-- ====================================================================== -->
  <!-- LAYER 02: PANEL (A) — TUMOR MICROENVIRONMENT ECOSYSTEM                 -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-02-panel-a" inkscape:label="02_Panel_A_TME_Ecosystem">
    <!-- Left Frame: Embedded TME Ecosystem Asset (flush top=62, bottom=322 with Cards 1-4) -->
    <rect id="pa-tme-frame" x="16" y="62.0" width="224" height="260.0" fill="#FFFFFF" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />
    <g id="pa-tme-artwork-clipped" clip-path="url(#clip-pa-tme)">
      <g id="pa-tme-artwork" inkscape:label="Embedded_TME_Ecosystem_Asset" transform="translate(16.0, 100.0) scale(0.28)">
        {bundle.tme_body}
      </g>
    </g>
    <rect id="pa-tme-badge-bg" x="22.0" y="68.0" width="144.0" height="22.0" fill="#FFFFFF" opacity="0.94" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.8" />
    <text id="pa-tme-badge-txt" x="94.0" y="83.0" font-size="11.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="middle">Heterogeneous TME</text>

    {"".join(pins_svg)}
    {"".join(cards_svg)}
  </g>"""


# ==============================================================================
# PANEL (B): SPECIMEN PROCESSING & SEQUENCING WORKFLOW RENDERER
# ==============================================================================
def render_wetlab_icon(step_id: WetLabStepId) -> str:
    """Render crisp vector illustration for each of the 5 wet-lab specimen processing steps."""
    match step_id:
        case "resection":
            return """
        <g id="icon-wetlab-resection" transform="translate(154, -34)">
          <!-- Trachea & Bronchi -->
          <rect x="33" y="2" width="8" height="22" rx="2" fill="#E2E8F0" stroke="#475569" stroke-width="1.1" />
          <line x1="33" y1="7" x2="41" y2="7" stroke="#64748B" stroke-width="0.9" />
          <line x1="33" y1="12" x2="41" y2="12" stroke="#64748B" stroke-width="0.9" />
          <line x1="33" y1="17" x2="41" y2="17" stroke="#64748B" stroke-width="0.9" />
          <!-- Right Lung (Viewer Left) with Resected Wedge Cutout -->
          <path d="M 32 22 C 24 18, 12 32, 10 58 C 12 64, 32 64, 33 56 Z" fill="#F4A2A2" stroke="#7F1D1D" stroke-width="1.2" />
          <!-- Resected Tumor Wedge (Exploded Top-Left) -->
          <path d="M 4 24 L 22 12 L 26 28 Z" fill="#F4A2A2" stroke="#DC2626" stroke-width="1.1" stroke-dasharray="2.5,1.5" />
          <ellipse cx="17" cy="21" rx="4.5" ry="3.0" fill="#FEF08A" stroke="#CA8A04" stroke-width="0.8" />
          <!-- Left Lung (Viewer Right) -->
          <path d="M 42 22 C 52 14, 64 28, 66 58 C 62 64, 42 64, 41 56 Z" fill="#F4A2A2" stroke="#7F1D1D" stroke-width="1.2" />
        </g>"""
        case "grossing":
            return """
        <g id="icon-wetlab-grossing" transform="translate(154, -32)">
          <!-- Petri Dish -->
          <ellipse cx="37" cy="38" rx="32" ry="20" fill="#F1F5F9" stroke="#64748B" stroke-width="1.1" />
          <ellipse cx="37" cy="35" rx="30" ry="18" fill="#F8FAFC" stroke="#94A3B8" stroke-width="0.8" />
          <!-- Tissue Wedge on Dish -->
          <polygon points="18,40 36,24 54,36" fill="#F4A2A2" stroke="#7F1D1D" stroke-width="1.1" />
          <ellipse cx="35" cy="33" rx="6.5" ry="3.8" fill="#FEF08A" stroke="#CA8A04" stroke-width="0.8" />
          <!-- Biopsy Punch Cores -->
          <circle cx="38" cy="31" r="2.2" fill="#FDE047" stroke="#713F12" stroke-width="0.8" />
          <circle cx="47" cy="35" r="2.2" fill="#FCA5A5" stroke="#7F1D1D" stroke-width="0.8" />
          <!-- Scalpel & Punch Tools -->
          <line x1="42" y1="20" x2="58" y2="4" stroke="#475569" stroke-width="2.2" />
          <polygon points="36,26 42,18 45,21" fill="#CBD5E1" stroke="#334155" stroke-width="0.8" />
          <line x1="48" y1="25" x2="66" y2="9" stroke="#475569" stroke-width="2.4" />
        </g>"""
        case "preservation":
            return """
        <g id="icon-wetlab-preservation" transform="translate(154, -32)">
          <!-- Cryovial (Snap-Frozen) -->
          <rect x="12" y="14" width="14" height="44" rx="4" fill="#E0F2FE" stroke="#0284C7" stroke-width="1.1" />
          <rect x="10" y="8" width="18" height="7" rx="1.5" fill="#93C5FD" stroke="#0369A1" stroke-width="1.1" />
          <ellipse cx="19" cy="40" rx="4.2" ry="5.5" fill="#F4A2A2" stroke="#991B1B" stroke-width="0.9" />
          <ellipse cx="19" cy="40" rx="2.2" ry="3.0" fill="#FEF08A" />
          <!-- Histology Cassette / Slide Block -->
          <polygon points="34,42 56,22 70,29 48,49" fill="#E2E8F0" stroke="#475569" stroke-width="1.1" />
          <ellipse cx="49" cy="37" rx="6.5" ry="2.8" transform="rotate(-35 49 37)" fill="#F4A2A2" stroke="#991B1B" stroke-width="0.9" />
        </g>"""
        case "rna_prep":
            return """
        <g id="icon-wetlab-rnaprep" transform="translate(154, -32)">
          <!-- Open Microcentrifuge Tube -->
          <path d="M 10 22 C 6 14, 14 10, 20 18" fill="#E2E8F0" stroke="#475569" stroke-width="1.1" />
          <polygon points="18,22 32,22 29,48 25,58 21,48" fill="#F1F5F9" stroke="#475569" stroke-width="1.1" />
          <polygon points="19.5,34 30.5,34 28.5,48 25,56 21.5,48" fill="#FCA5A5" />
          <!-- RNA Hairpins & Adapter-Ligated Fragments -->
          <path d="M 30 18 Q 36 8 42 18 T 54 18" fill="none" stroke="#451A03" stroke-width="1.3" />
          <circle cx="38" cy="7" r="3.5" fill="none" stroke="#451A03" stroke-width="1.2" />
          <!-- Adapter-ligated cDNA duplexes -->
          <rect x="46" y="14" width="5" height="3.2" fill="#FCA5A5" stroke="#7F1D1D" stroke-width="0.6" />
          <rect x="51" y="14" width="12" height="3.2" fill="#E2E8F0" stroke="#334155" stroke-width="0.6" />
          <rect x="63" y="14" width="5" height="3.2" fill="#FDE047" stroke="#854D0E" stroke-width="0.6" />
          <rect x="44" y="21" width="5" height="3.2" fill="#FCA5A5" stroke="#7F1D1D" stroke-width="0.6" />
          <rect x="49" y="21" width="12" height="3.2" fill="#E2E8F0" stroke="#334155" stroke-width="0.6" />
          <rect x="61" y="21" width="5" height="3.2" fill="#FDE047" stroke="#854D0E" stroke-width="0.6" />
        </g>"""
        case "sequencing":
            return """
        <g id="icon-wetlab-sequencing" transform="translate(154, -34)">
          <!-- Benchtop High-Throughput Sequencer -->
          <rect x="6" y="8" width="46" height="38" rx="4" fill="#CBD5E1" stroke="#334155" stroke-width="1.2" />
          <rect x="11" y="13" width="26" height="18" rx="1.5" fill="#64748B" stroke="#1E293B" stroke-width="1.0" />
          <line x1="15" y1="29" x2="26" y2="15" stroke="#94A3B8" stroke-width="1.8" />
          <!-- Flow Cell Cartridge -->
          <rect x="30" y="34" width="12" height="26" rx="2.5" fill="#E0F2FE" stroke="#1E293B" stroke-width="1.1" />
          <line x1="34" y1="38" x2="34" y2="56" stroke="#0284C7" stroke-width="0.8" />
          <line x1="38" y1="38" x2="38" y2="56" stroke="#0284C7" stroke-width="0.8" />
          <!-- Read Pileup Alignment Bars -->
          <rect x="45" y="38" width="10" height="2.6" fill="#F87171" />
          <rect x="56" y="38" width="10" height="2.6" fill="#FACC15" />
          <rect x="46" y="43" width="10" height="2.6" fill="#F87171" />
          <rect x="57" y="43" width="10" height="2.6" fill="#60A5FA" />
          <rect x="44" y="48" width="10" height="2.6" fill="#60A5FA" />
          <rect x="55" y="48" width="11" height="2.6" fill="#F87171" />
          <rect x="47" y="53" width="10" height="2.6" fill="#60A5FA" />
          <rect x="58" y="53" width="10" height="2.6" fill="#FACC15" />
        </g>"""
        case _ as unreachable:
            assert_never(unreachable)


def render_panel_b_wetlab_workflow() -> str:
    """Render Panel (b): 5-step clinical specimen processing and sequencing workflow."""
    steps_svg: list[str] = []
    arrows_svg: list[str] = []

    for idx, spec in enumerate(WETLAB_STEPS):
        text_lines_svg: list[str] = []
        n_lines = len(spec.lines) + (1 if spec.sub_line else 0)
        start_y = -5.0 * (n_lines - 1)
        for l_idx, line_str in enumerate(spec.lines):
            ly = start_y + l_idx * 13.5
            text_lines_svg.append(
                f'<text x="134" y="{ly}" font-size="11.2" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">{html.escape(line_str)}</text>'
            )
        if spec.sub_line:
            sub_y = start_y + len(spec.lines) * 13.5
            text_lines_svg.append(
                f'<text x="134" y="{sub_y}" font-size="10.4" font-weight="500" fill="{COLOR_TEXT_MUTED}" text-anchor="end">{html.escape(spec.sub_line)}</text>'
            )

        steps_svg.append(f"""
      <g id="pb-step-{idx + 1}-{spec.step_id}" transform="translate(476, {spec.center_y})">
        {"".join(text_lines_svg)}
        <g transform="translate(10, 2) scale(0.84)">
          {render_wetlab_icon(spec.step_id)}
        </g>
      </g>""")

        if idx < len(WETLAB_STEPS) - 1:
            next_y = WETLAB_STEPS[idx + 1].center_y
            mid_y = (spec.center_y + next_y) / 2.0
            arrows_svg.append(f"""
      <g id="pb-arrow-{idx + 1}" transform="translate(640, {mid_y})">
        <line x1="0" y1="-8.5" x2="0" y2="2.5" stroke="#0F172A" stroke-width="2.6" />
        <polygon points="-5.2,1.5 5.2,1.5 0,9.0" fill="#0F172A" />
      </g>""")

    return f"""
  <!-- ====================================================================== -->
  <!-- LAYER 03: PANEL (B) — SPECIMEN PROCESSING & SEQUENCING WORKFLOW        -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-03-panel-b" inkscape:label="03_Panel_B_Specimen_Sequencing">
    {"".join(steps_svg)}
    {"".join(arrows_svg)}
  </g>"""


# ==============================================================================
# PANEL (C): DECONVOLUTION & CLASSIFICATION PIPELINE RENDERER
# ==============================================================================
def render_panel_c_deconv_pipeline(bundle: AssetBundle) -> str:
    """Render Panel (c): Single-Cell & Bulk Deconvolution Classification Pipeline mimicking Image 1."""
    # Y-centers for the 6 left-stream steps inside y = 44..364
    y_steps: tuple[float, ...] = (66.0, 120.0, 174.0, 228.0, 282.0, 336.0)

    down_arrows: list[str] = []
    for idx in range(len(y_steps) - 1):
        mid_y = (y_steps[idx] + y_steps[idx + 1]) / 2.0
        down_arrows.append(f"""
      <g id="pc-down-arrow-{idx + 1}" transform="translate(748, {mid_y})">
        <line x1="0" y1="-5.0" x2="0" y2="1.5" stroke="#0F172A" stroke-width="2.2" />
        <polygon points="-4.2,0.5 4.2,0.5 0,6.2" fill="#0F172A" />
      </g>""")

    return f"""
  <!-- ====================================================================== -->
  <!-- LAYER 04: PANEL (C) — DECONVOLUTION & CLASSIFICATION PIPELINE          -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-04-panel-c" inkscape:label="04_Panel_C_Deconvolution_Pipeline">

    <!-- Step 1: Single Cell Transcriptomic Datasets -->
    <g id="pc-step-1-scrna" transform="translate(722, {y_steps[0]})">
      <g id="pc-icon-scrna" transform="translate(0, -20) scale(0.76)">
        <!-- Axes -->
        <line x1="4" y1="50" x2="4" y2="6" stroke="#334155" stroke-width="1.2" />
        <polygon points="1.5,8 6.5,8 4,3" fill="#334155" />
        <line x1="4" y1="50" x2="54" y2="50" stroke="#334155" stroke-width="1.2" />
        <polygon points="52,47.5 52,52.5 57,50" fill="#334155" />
        <!-- 4 Dashed UMAP Clusters -->
        <ellipse cx="18" cy="20" rx="9" ry="7" fill="none" stroke="#475569" stroke-width="0.9" stroke-dasharray="2.0,1.5" />
        <circle cx="15" cy="18" r="2.1" fill="#93C5FD" stroke="#1E3A8A" stroke-width="0.5" />
        <circle cx="20" cy="17" r="2.1" fill="#93C5FD" stroke="#1E3A8A" stroke-width="0.5" />
        <circle cx="18" cy="22" r="2.1" fill="#93C5FD" stroke="#1E3A8A" stroke-width="0.5" />

        <ellipse cx="40" cy="16" rx="11" ry="8" fill="none" stroke="#475569" stroke-width="0.9" stroke-dasharray="2.0,1.5" />
        <circle cx="36" cy="14" r="2.1" fill="#FDE047" stroke="#713F12" stroke-width="0.5" />
        <circle cx="42" cy="13" r="2.1" fill="#FDE047" stroke="#713F12" stroke-width="0.5" />
        <circle cx="39" cy="18" r="2.1" fill="#FDE047" stroke="#713F12" stroke-width="0.5" />
        <circle cx="45" cy="17" r="2.1" fill="#FDE047" stroke="#713F12" stroke-width="0.5" />

        <ellipse cx="22" cy="38" rx="10" ry="7" fill="none" stroke="#475569" stroke-width="0.9" stroke-dasharray="2.0,1.5" />
        <circle cx="18" cy="38" r="2.1" fill="#86EFAC" stroke="#14532D" stroke-width="0.5" />
        <circle cx="23" cy="36" r="2.1" fill="#86EFAC" stroke="#14532D" stroke-width="0.5" />
        <circle cx="25" cy="40" r="2.1" fill="#86EFAC" stroke="#14532D" stroke-width="0.5" />

        <ellipse cx="44" cy="36" rx="9" ry="7" fill="none" stroke="#475569" stroke-width="0.9" stroke-dasharray="2.0,1.5" />
        <circle cx="41" cy="35" r="2.1" fill="#FCA5A5" stroke="#7F1D1D" stroke-width="0.5" />
        <circle cx="46" cy="34" r="2.1" fill="#FCA5A5" stroke="#7F1D1D" stroke-width="0.5" />
        <circle cx="44" cy="39" r="2.1" fill="#FCA5A5" stroke="#7F1D1D" stroke-width="0.5" />
      </g>
      <text x="56" y="-8" font-size="10.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Single Cell</text>
      <text x="56" y="4" font-size="10.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Transcriptomic</text>
      <text x="56" y="16" font-size="10.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Datasets</text>
    </g>

    <!-- Step 2: Preprocessing (QC, HVG, PCA Harmonization) -->
    <g id="pc-step-2-preproc" transform="translate(722, {y_steps[1]})">
      <g id="pc-icon-preproc" transform="translate(0, -20) scale(0.76)">
        <!-- Mini Expression Matrix -->
        <rect x="4" y="6" width="24" height="42" fill="#F8FAFC" stroke="#334155" stroke-width="0.9" />
        <rect x="6" y="8" width="6" height="6" fill="#FCA5A5" />
        <rect x="13" y="8" width="6" height="6" fill="#FDE047" />
        <rect x="20" y="8" width="6" height="6" fill="#93C5FD" />
        <rect x="6" y="16" width="6" height="6" fill="#FCA5A5" />
        <rect x="13" y="16" width="6" height="6" fill="#93C5FD" />
        <rect x="20" y="16" width="6" height="6" fill="#FDE047" />
        <rect x="6" y="24" width="6" height="6" fill="#FDE047" />
        <rect x="13" y="24" width="6" height="6" fill="#FCA5A5" />
        <rect x="20" y="24" width="6" height="6" fill="#93C5FD" />
        <rect x="6" y="32" width="6" height="6" fill="#93C5FD" />
        <rect x="13" y="32" width="6" height="6" fill="#FDE047" />
        <rect x="20" y="32" width="6" height="6" fill="#FCA5A5" />
        <rect x="6" y="40" width="6" height="6" fill="#CBD5E1" />
        <rect x="13" y="40" width="6" height="6" fill="#93C5FD" />
        <rect x="20" y="40" width="6" height="6" fill="#FCA5A5" />
        <!-- kNN Graph + Curve -->
        <line x1="36" y1="14" x2="46" y2="10" stroke="#334155" stroke-width="0.9" />
        <line x1="46" y1="10" x2="54" y2="16" stroke="#334155" stroke-width="0.9" />
        <line x1="36" y1="14" x2="44" y2="22" stroke="#334155" stroke-width="0.9" />
        <line x1="46" y1="10" x2="44" y2="22" stroke="#334155" stroke-width="0.9" />
        <circle cx="36" cy="14" r="3.0" fill="#FCA5A5" stroke="#1E293B" stroke-width="0.7" />
        <circle cx="46" cy="10" r="3.0" fill="#93C5FD" stroke="#1E293B" stroke-width="0.7" />
        <circle cx="54" cy="16" r="3.0" fill="#FDE047" stroke="#1E293B" stroke-width="0.7" />
        <circle cx="44" cy="22" r="3.0" fill="#FCA5A5" stroke="#1E293B" stroke-width="0.7" />
        <polyline points="32,28 32,48 56,48" fill="none" stroke="#334155" stroke-width="1.0" />
        <path d="M 34 46 C 42 45, 46 32, 55 30" fill="none" stroke="#1E293B" stroke-width="1.2" />
      </g>
      <text x="56" y="-8" font-size="10.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}">Preprocessing</text>
      <text x="56" y="4" font-size="10.0" font-weight="500" fill="{COLOR_TEXT_MUTED}">(QC, HVG, PCA</text>
      <text x="56" y="16" font-size="10.0" font-weight="500" fill="{COLOR_TEXT_MUTED}">Harmonization)</text>
    </g>

    <!-- Step 3: Cell State Reference Construction (Embeds 4 cell icons from assets/) -->
    <g id="pc-step-3-ref" transform="translate(722, {y_steps[2]})">
      <g id="pc-icon-ref" transform="translate(0, -20) scale(0.76)">
        <!-- 4 Cell Assets on Left -->
        <g transform="translate(2, 2) scale(0.045)">{bundle.cell_cancer}</g>
        <g transform="translate(2, 15) scale(0.05)">{bundle.cell_t_lymph}</g>
        <g transform="translate(2, 28) scale(0.048)">{bundle.cell_macro_m2}</g>
        <g transform="translate(1, 41) scale(0.042)">{bundle.cell_caf}</g>
        <!-- Dotted Mapping Lines -->
        <line x1="17" y1="8" x2="25" y2="8" stroke="#334155" stroke-width="0.9" stroke-dasharray="1.5,1.5" />
        <line x1="17" y1="21" x2="25" y2="21" stroke="#334155" stroke-width="0.9" stroke-dasharray="1.5,1.5" />
        <line x1="17" y1="34" x2="25" y2="34" stroke="#334155" stroke-width="0.9" stroke-dasharray="1.5,1.5" />
        <line x1="17" y1="46" x2="25" y2="46" stroke="#334155" stroke-width="0.9" stroke-dasharray="1.5,1.5" />
        <!-- 4-Column Reference Signature Matrix -->
        <g transform="translate(26, 4)">
          <rect x="0" y="0" width="7" height="44" fill="#FDE047" stroke="#1E293B" stroke-width="0.7" />
          <rect x="8" y="0" width="7" height="44" fill="#93C5FD" stroke="#1E293B" stroke-width="0.7" />
          <rect x="16" y="0" width="7" height="44" fill="#FCA5A5" stroke="#1E293B" stroke-width="0.7" />
          <rect x="24" y="0" width="7" height="44" fill="#86EFAC" stroke="#1E293B" stroke-width="0.7" />
          <line x1="0" y1="9" x2="31" y2="9" stroke="#1E293B" stroke-width="0.6" />
          <line x1="0" y1="18" x2="31" y2="18" stroke="#1E293B" stroke-width="0.6" />
          <line x1="0" y1="27" x2="31" y2="27" stroke="#1E293B" stroke-width="0.6" />
          <line x1="0" y1="36" x2="31" y2="36" stroke="#1E293B" stroke-width="0.6" />
        </g>
      </g>
      <text x="56" y="-8" font-size="10.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Cell State</text>
      <text x="56" y="4" font-size="10.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Reference</text>
      <text x="56" y="16" font-size="10.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Construction</text>
    </g>

    <!-- Step 4: Deconvolution (Prism + 4 Posterior Density Curves) -->
    <g id="pc-step-4-deconv" transform="translate(722, {y_steps[3]})">
      <g id="pc-icon-deconv" transform="translate(0, -20) scale(0.76)">
        <!-- Braided Incoming Bulk Mixture Beam -->
        <path d="M 0 26 Q 5 22 10 26 T 18 26" fill="none" stroke="#FCA5A5" stroke-width="2.0" />
        <path d="M 0 28 Q 5 32 10 28 T 18 28" fill="none" stroke="#93C5FD" stroke-width="2.0" />
        <path d="M 0 27 L 18 27" stroke="#FDE047" stroke-width="1.4" />
        <!-- Optical Prism -->
        <polygon points="24,8 10,42 38,42" fill="#E0F2FE" stroke="#1E293B" stroke-width="1.1" />
        <!-- Dispersed Rays -->
        <polygon points="29,19 46,6 46,12 31,23" fill="#FCA5A5" opacity="0.85" />
        <polygon points="31,23 46,18 46,24 33,27" fill="#FDE047" opacity="0.85" />
        <polygon points="33,27 46,30 46,36 35,31" fill="#93C5FD" opacity="0.85" />
        <polygon points="35,31 46,42 46,48 37,36" fill="#86EFAC" opacity="0.85" />
        <!-- 4 Mini Posterior Density Curves -->
        <path d="M 47 11 Q 52 2 57 11" fill="#FCA5A5" stroke="#1E293B" stroke-width="0.7" />
        <path d="M 47 23 Q 52 14 57 23" fill="#FDE047" stroke="#1E293B" stroke-width="0.7" />
        <path d="M 47 35 Q 52 26 57 35" fill="#93C5FD" stroke="#1E293B" stroke-width="0.7" />
        <path d="M 47 47 Q 52 38 57 47" fill="#86EFAC" stroke="#1E293B" stroke-width="0.7" />
      </g>
      <text x="56" y="-2" font-size="11.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Deconvolution</text>
      <text x="56" y="10" font-size="10.4" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">(Our Approach)</text>
    </g>

    <!-- Step 5: Driver Cell States Identification (Funnel + Logistic Sigmoid Curve) -->
    <g id="pc-step-5-driver" transform="translate(722, {y_steps[4]})">
      <g id="pc-icon-driver" transform="translate(0, -20) scale(0.76)">
        <!-- Axes -->
        <line x1="12" y1="48" x2="12" y2="6" stroke="#334155" stroke-width="1.0" />
        <line x1="12" y1="48" x2="56" y2="48" stroke="#334155" stroke-width="1.0" />
        <line x1="12" y1="27" x2="56" y2="27" stroke="#64748B" stroke-width="0.8" stroke-dasharray="2.0,1.5" />
        <!-- Filter Funnel on Y-axis -->
        <polygon points="4,20 20,20 14,29 14,35 10,37 10,29" fill="#E0F2FE" stroke="#1E293B" stroke-width="0.9" />
        <!-- Sigmoid Curve & Dots -->
        <path d="M 15 43 C 28 43, 30 11, 52 11" fill="none" stroke="#1E293B" stroke-width="1.3" />
        <circle cx="32" cy="12" r="2.8" fill="#93C5FD" stroke="#1E293B" stroke-width="0.7" />
        <circle cx="38" cy="14" r="2.8" fill="#FCA5A5" stroke="#1E293B" stroke-width="0.7" />
        <circle cx="44" cy="10" r="2.8" fill="#FDE047" stroke="#1E293B" stroke-width="0.7" />
        <circle cx="22" cy="44" r="2.8" fill="#CBD5E1" stroke="#1E293B" stroke-width="0.7" />
        <circle cx="29" cy="43" r="2.8" fill="#CBD5E1" stroke="#1E293B" stroke-width="0.7" />
      </g>
      <text x="56" y="-2" font-size="10.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Driver Cell States</text>
      <text x="56" y="10" font-size="10.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Identification</text>
    </g>

    <!-- Step 6: Response Classification (Random Forest Trees -> Consensus Output) -->
    <g id="pc-step-6-rf" transform="translate(722, {y_steps[5]})">
      <g id="pc-icon-rf" transform="translate(0, -20) scale(0.76)">
        <!-- Tree 1 -->
        <line x1="11" y1="8" x2="6" y2="18" stroke="#1E293B" stroke-width="0.9" />
        <line x1="11" y1="8" x2="16" y2="18" stroke="#1E293B" stroke-width="0.9" />
        <circle cx="11" cy="8" r="2.5" fill="#FFFFFF" stroke="#1E293B" stroke-width="0.9" />
        <rect x="3" y="18" width="5" height="6" rx="1" fill="#93C5FD" stroke="#1E293B" stroke-width="0.6" />
        <rect x="9" y="18" width="5" height="6" rx="1" fill="#FCA5A5" stroke="#1E293B" stroke-width="0.6" />
        <rect x="15" y="18" width="5" height="6" rx="1" fill="#86EFAC" stroke="#1E293B" stroke-width="0.6" />
        <!-- Tree 2 -->
        <line x1="30" y1="8" x2="25" y2="18" stroke="#1E293B" stroke-width="0.9" />
        <line x1="30" y1="8" x2="35" y2="18" stroke="#1E293B" stroke-width="0.9" />
        <circle cx="30" cy="8" r="2.5" fill="#FFFFFF" stroke="#1E293B" stroke-width="0.9" />
        <rect x="22" y="18" width="5" height="6" rx="1" fill="#FCA5A5" stroke="#1E293B" stroke-width="0.6" />
        <rect x="28" y="18" width="5" height="6" rx="1" fill="#FDE047" stroke="#1E293B" stroke-width="0.6" />
        <rect x="34" y="18" width="5" height="6" rx="1" fill="#93C5FD" stroke="#1E293B" stroke-width="0.6" />
        <!-- Tree 3 -->
        <line x1="49" y1="8" x2="44" y2="18" stroke="#1E293B" stroke-width="0.9" />
        <line x1="49" y1="8" x2="54" y2="18" stroke="#1E293B" stroke-width="0.9" />
        <circle cx="49" cy="8" r="2.5" fill="#FFFFFF" stroke="#1E293B" stroke-width="0.9" />
        <rect x="41" y="18" width="5" height="6" rx="1" fill="#FDE047" stroke="#1E293B" stroke-width="0.6" />
        <rect x="47" y="18" width="5" height="6" rx="1" fill="#FCA5A5" stroke="#1E293B" stroke-width="0.6" />
        <rect x="53" y="18" width="5" height="6" rx="1" fill="#86EFAC" stroke="#1E293B" stroke-width="0.6" />
        <!-- Dotted Aggregation Arrows & Output Pill -->
        <path d="M 11 26 L 11 33 L 26 40" fill="none" stroke="#334155" stroke-width="0.9" stroke-dasharray="1.5,1.5" />
        <line x1="30" y1="26" x2="30" y2="39" stroke="#334155" stroke-width="0.9" stroke-dasharray="1.5,1.5" />
        <path d="M 49 26 L 49 33 L 34 40" fill="none" stroke="#334155" stroke-width="0.9" stroke-dasharray="1.5,1.5" />
        <rect x="21" y="41" width="18" height="7" rx="3" fill="#FCA5A5" stroke="#1E293B" stroke-width="1.0" />
      </g>
      <text x="56" y="-2" font-size="10.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Response</text>
      <text x="56" y="10" font-size="10.8" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Classification</text>
    </g>

    <!-- Downward Arrows for Left Stream -->
    {"".join(down_arrows)}

    <!-- Right Stream Step 1: Patient Cohort (5 Silhouettes) -->
    <g id="pc-right-cohort" transform="translate(918, 68)">
      <text x="44" y="-3" font-size="10.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">Patient</text>
      <text x="44" y="10" font-size="10.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">Cohort</text>
      <!-- 5 Overlapping Patient Bust Silhouettes -->
      <g id="pc-icon-patients" transform="translate(50, -20) scale(0.84)">
        <!-- Back Row Silhouettes -->
        <circle cx="18" cy="12" r="5.5" fill="#F5D0B5" stroke="#451A03" stroke-width="0.9" />
        <path d="M 9 36 C 9 23, 27 23, 27 36 Z" fill="#F5D0B5" stroke="#451A03" stroke-width="0.9" />
        <circle cx="42" cy="12" r="5.5" fill="#F5D0B5" stroke="#451A03" stroke-width="0.9" />
        <path d="M 33 36 C 33 23, 51 23, 51 36 Z" fill="#F5D0B5" stroke="#451A03" stroke-width="0.9" />
        <!-- Front Row Silhouettes -->
        <circle cx="8" cy="17" r="5.5" fill="#F5D0B5" stroke="#451A03" stroke-width="1.0" />
        <path d="M -1 44 C -1 28, 17 28, 17 44 Z" fill="#F5D0B5" stroke="#451A03" stroke-width="1.0" />
        <circle cx="52" cy="17" r="5.5" fill="#F5D0B5" stroke="#451A03" stroke-width="1.0" />
        <path d="M 43 44 C 43 28, 61 28, 61 44 Z" fill="#F5D0B5" stroke="#451A03" stroke-width="1.0" />
        <circle cx="30" cy="16" r="6.0" fill="#F5D0B5" stroke="#451A03" stroke-width="1.1" />
        <path d="M 19 44 C 19 27, 41 27, 41 44 Z" fill="#F5D0B5" stroke="#451A03" stroke-width="1.1" />
      </g>
    </g>

    <!-- Arrow: Patient Cohort -> Bulk RNA-sequencing data -->
    <g id="pc-right-arrow-1" transform="translate(996, 108)">
      <line x1="0" y1="-6" x2="0" y2="2" stroke="#0F172A" stroke-width="2.4" />
      <polygon points="-4.8,1.0 4.8,1.0 0,7.5" fill="#0F172A" />
    </g>

    <!-- Right Stream Step 2: Bulk RNA-sequencing data -->
    <g id="pc-right-bulk" transform="translate(918, 152)">
      <text x="44" y="-8" font-size="10.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">Bulk RNA-</text>
      <text x="44" y="4" font-size="10.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">sequencing</text>
      <text x="44" y="16" font-size="10.8" font-weight="600" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">data</text>
      <!-- RNA Hairpins + Bulk Expression Matrix -->
      <g id="pc-icon-bulk-matrix" transform="translate(50, -22) scale(0.84)">
        <!-- Left RNA Hairpin Loops & Fragments -->
        <path d="M 4 14 Q 12 6 18 14 T 28 14" fill="none" stroke="#451A03" stroke-width="1.1" />
        <circle cx="12" cy="6" r="3.0" fill="none" stroke="#451A03" stroke-width="1.0" />
        <path d="M 2 28 Q 10 20 18 28 T 28 28" fill="none" stroke="#451A03" stroke-width="1.1" />
        <circle cx="8" cy="21" r="2.8" fill="none" stroke="#451A03" stroke-width="1.0" />
        <path d="M 4 42 Q 12 34 20 42 T 28 42" fill="none" stroke="#451A03" stroke-width="1.1" />
        <rect x="12" y="16" width="12" height="3.0" fill="#FCA5A5" stroke="#451A03" stroke-width="0.6" />
        <rect x="10" y="32" width="12" height="3.0" fill="#FDE047" stroke="#451A03" stroke-width="0.6" />
        <rect x="14" y="46" width="12" height="3.0" fill="#93C5FD" stroke="#451A03" stroke-width="0.6" />
        <!-- Right 5x7 Colored Bulk Expression Matrix -->
        <g transform="translate(30, 4)">
          <rect x="0" y="0" width="32" height="48" fill="#F8FAFC" stroke="#1E293B" stroke-width="0.9" />
          <rect x="1.5" y="1.5" width="5.5" height="6.0" fill="#FCA5A5" />
          <rect x="7.5" y="1.5" width="5.5" height="6.0" fill="#FCA5A5" />
          <rect x="13.5" y="1.5" width="5.5" height="6.0" fill="#FDE047" />
          <rect x="19.5" y="1.5" width="5.5" height="6.0" fill="#93C5FD" />
          <rect x="25.5" y="1.5" width="5.0" height="6.0" fill="#93C5FD" />

          <rect x="1.5" y="8.2" width="5.5" height="6.0" fill="#FCA5A5" />
          <rect x="7.5" y="8.2" width="5.5" height="6.0" fill="#FDE047" />
          <rect x="13.5" y="8.2" width="5.5" height="6.0" fill="#FDE047" />
          <rect x="19.5" y="8.2" width="5.5" height="6.0" fill="#93C5FD" />
          <rect x="25.5" y="8.2" width="5.0" height="6.0" fill="#93C5FD" />

          <rect x="1.5" y="14.9" width="5.5" height="6.0" fill="#FCA5A5" />
          <rect x="7.5" y="14.9" width="5.5" height="6.0" fill="#FCA5A5" />
          <rect x="13.5" y="14.9" width="5.5" height="6.0" fill="#FDE047" />
          <rect x="19.5" y="14.9" width="5.5" height="6.0" fill="#93C5FD" />
          <rect x="25.5" y="14.9" width="5.0" height="6.0" fill="#86EFAC" />

          <rect x="1.5" y="21.6" width="5.5" height="6.0" fill="#86EFAC" />
          <rect x="7.5" y="21.6" width="5.5" height="6.0" fill="#93C5FD" />
          <rect x="13.5" y="21.6" width="5.5" height="6.0" fill="#FCA5A5" />
          <rect x="19.5" y="21.6" width="5.5" height="6.0" fill="#86EFAC" />
          <rect x="25.5" y="21.6" width="5.0" height="6.0" fill="#86EFAC" />

          <rect x="1.5" y="28.3" width="5.5" height="6.0" fill="#93C5FD" />
          <rect x="7.5" y="28.3" width="5.5" height="6.0" fill="#86EFAC" />
          <rect x="13.5" y="28.3" width="5.5" height="6.0" fill="#FDE047" />
          <rect x="19.5" y="28.3" width="5.5" height="6.0" fill="#FCA5A5" />
          <rect x="25.5" y="28.3" width="5.0" height="6.0" fill="#FCA5A5" />

          <rect x="1.5" y="35.0" width="5.5" height="6.0" fill="#86EFAC" />
          <rect x="7.5" y="35.0" width="5.5" height="6.0" fill="#86EFAC" />
          <rect x="13.5" y="35.0" width="5.5" height="6.0" fill="#93C5FD" />
          <rect x="19.5" y="35.0" width="5.5" height="6.0" fill="#FDE047" />
          <rect x="25.5" y="35.0" width="5.0" height="6.0" fill="#FCA5A5" />

          <rect x="1.5" y="41.5" width="5.5" height="5.5" fill="#86EFAC" />
          <rect x="7.5" y="41.5" width="5.5" height="5.5" fill="#93C5FD" />
          <rect x="13.5" y="41.5" width="5.5" height="5.5" fill="#FDE047" />
          <rect x="19.5" y="41.5" width="5.5" height="5.5" fill="#FCA5A5" />
          <rect x="25.5" y="41.5" width="5.0" height="5.5" fill="#FCA5A5" />
        </g>
      </g>
    </g>

    <!-- L-Shaped Connector: Bulk RNA-seq data -> Deconvolution -->
    <g id="pc-l-connector-bulk-to-deconv">
      <polyline points="996,184 996,228 866,228" fill="none" stroke="#0F172A" stroke-width="1.5" />
      <polygon points="868,224.0 868,232.0 860,228" fill="#0F172A" />
    </g>
  </g>"""


# ==============================================================================
# PANEL (D): BULK RNA-SEQ ICB CLINICAL COHORTS TABLE RENDERER
# ==============================================================================
def render_bulk_cohort_row(idx: int, item: BulkTrialCohortSpec) -> str:
    """Render a single bulk RNA-seq trial row in Panel (d) with study names inside the table like Panel (e)."""
    row_h = 28
    y_offset = 28 + idx * row_h
    table_w = 1044
    bar_w = 240
    raw_w_r = int(round((item.n_resp / item.n_total) * bar_w))
    w_r = max(58, min(bar_w - 68, raw_w_r))
    w_nr = bar_w - w_r
    row_bg = COLOR_PANEL_BG if idx % 2 == 0 else COLOR_SUBTLE_FILL
    r_center = w_r // 2
    nr_center = w_r + (w_nr // 2)

    return f"""
      <g id="pd-bulk-row-{idx + 1:02d}" transform="translate(0, {y_offset})">
        <rect x="0" y="0" width="{table_w}" height="{row_h}" fill="{row_bg}" />
        <line x1="0" y1="{row_h}" x2="{table_w}" y2="{row_h}" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.6" />
        <!-- Primary Study / Trial Cohort inside table (matching Panel e) -->
        <text x="10" y="18.5" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">{html.escape(item.study)}</text>
        <!-- Cancer Type with Okabe-Ito indicator dot -->
        <circle cx="226" cy="14.0" r="4.5" fill="{item.cancer_color}" />
        <text x="238" y="18.5" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">{html.escape(item.cancer_type)}</text>
        <!-- Regimen & Target -->
        <text x="382" y="18.5" font-size="12.0" font-weight="500" fill="{COLOR_TEXT_SECONDARY}">{html.escape(item.regimen)}</text>
        <!-- Sample Size N -->
        <text x="662" y="18.5" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">{item.n_total}</text>
        <!-- Stacked RECIST Bar (R / NR) -->
        <g transform="translate(684, 4.0)">
          <rect x="0" y="0" width="{w_r}" height="20" fill="{OKABE_BLUISH_GREEN}" />
          <rect x="{w_r}" y="0" width="{w_nr}" height="20" fill="{OKABE_VERMILION}" />
          <text x="{r_center}" y="14.5" font-size="11.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">{item.n_resp}R</text>
          <text x="{nr_center}" y="14.5" font-size="11.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">{item.n_non_resp}NR</text>
        </g>
        <!-- ORR % -->
        <text x="1032" y="18.5" font-size="12.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="end">{item.orr_pct:.1f}%</text>
      </g>"""


def render_panel_d_bulk_cohorts() -> str:
    """Render Panel (d): Bulk RNA-seq ICB Clinical Cohorts (n = 1,097) table & legend."""
    total_n = sum(c.n_total for c in BULK_TRIAL_COHORTS)
    total_r = sum(c.n_resp for c in BULK_TRIAL_COHORTS)
    total_nr = sum(c.n_non_resp for c in BULK_TRIAL_COHORTS)
    total_orr = (total_r / total_n) * 100.0

    row_h = 28
    table_w = 1044
    bar_w = 240
    total_w_r = max(68, int(round((total_r / total_n) * bar_w)))
    total_w_nr = bar_w - total_w_r
    summary_y = 28 + len(BULK_TRIAL_COHORTS) * row_h

    rows_svg = "".join(render_bulk_cohort_row(idx, item) for idx, item in enumerate(BULK_TRIAL_COHORTS))

    return f"""
  <!-- ====================================================================== -->
  <!-- LAYER 05: PANEL (D) — BULK RNA-SEQ ICB CLINICAL COHORTS (N = 1,097)    -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-05-panel-d" inkscape:label="05_Panel_D_Bulk_ICB_Cohorts">
    <g id="pd-table-group" transform="translate(18, 424)">
      <!-- Outer Table Border -->
      <rect x="0" y="0" width="{table_w}" height="{summary_y + row_h}" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />

      <!-- Header Bar -->
      <rect x="0" y="0" width="{table_w}" height="{row_h}" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.8" />
      <text x="10" y="19.0" font-size="13.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Primary Study</text>
      <text x="238" y="19.0" font-size="13.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Cancer</text>
      <text x="382" y="19.0" font-size="13.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Regimen &amp; Target</text>
      <text x="662" y="19.0" font-size="13.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">N</text>
      <text x="684" y="19.0" font-size="13.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">RECIST (R / NR)</text>
      <text x="1032" y="19.0" font-size="13.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">ORR</text>

      {rows_svg}

      <!-- Summary Row: Overall (9 Cohorts) -->
      <g id="pd-bulk-row-summary" transform="translate(0, {summary_y})">
        <rect x="0" y="0" width="{table_w}" height="{row_h}" fill="{COLOR_SUMMARY_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.8" />
        <text x="10" y="18.5" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Overall (9 Cohorts)</text>
        <text x="238" y="18.5" font-size="12.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">5 Solid Types</text>
        <text x="382" y="18.5" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Anti-PD-(L)1 ± Combo</text>
        <text x="662" y="18.5" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">{total_n:,}</text>
        <g transform="translate(684, 4.0)">
          <rect x="0" y="0" width="{total_w_r}" height="20" fill="{OKABE_BLUISH_GREEN}" />
          <rect x="{total_w_r}" y="0" width="{total_w_nr}" height="20" fill="{OKABE_VERMILION}" />
          <text x="{total_w_r // 2}" y="14.5" font-size="11.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">{total_r}R</text>
          <text x="{total_w_r + (total_w_nr // 2)}" y="14.5" font-size="11.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">{total_nr}NR</text>
        </g>
        <text x="1032" y="18.5" font-size="12.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="end">{total_orr:.1f}%</text>
      </g>

      <!-- Bottom Legend Box: RECIST v1.1 Response -->
      <g id="pd-recist-legend" transform="translate(0, {summary_y + row_h + 10})">
        <rect x="0" y="0" width="{table_w}" height="32" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.8" />
        <text x="14" y="21" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">RECIST v1.1 Response:</text>
        <rect x="200" y="9" width="14" height="14" fill="{OKABE_BLUISH_GREEN}" />
        <text x="220" y="21" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">R: Responder (357, 32.5%)</text>
        <rect x="460" y="9" width="14" height="14" fill="{OKABE_VERMILION}" />
        <text x="480" y="21" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">NR: Non-Responder (740, 67.5%)</text>
      </g>
    </g>
  </g>"""


# ==============================================================================
# PANEL (E): SINGLE-CELL MILO COMPENDIUM TABLE RENDERER
# ==============================================================================
def render_sc_milo_cohort_row(idx: int, item: SingleCellMiloCohortSpec) -> str:
    """Render a single-cell ICB cohort row in Panel (e) mimicking Image 5."""
    row_h = 28
    y_offset = 28 + idx * row_h
    table_w = 1044
    bar_w = 240
    n_eval = item.n_resp + item.n_non_resp
    raw_w_r = int(round((item.n_resp / n_eval) * bar_w))
    w_r = max(56, min(bar_w - 64, raw_w_r))
    w_nr = bar_w - w_r
    row_bg = COLOR_PANEL_BG if idx % 2 == 0 else COLOR_SUBTLE_FILL
    r_center = w_r // 2
    nr_center = w_r + (w_nr // 2)

    return f"""
      <g id="pe-sc-row-{idx + 1:02d}" transform="translate(0, {y_offset})">
        <rect x="0" y="0" width="{table_w}" height="{row_h}" fill="{row_bg}" />
        <line x1="0" y1="{row_h}" x2="{table_w}" y2="{row_h}" stroke="{COLOR_DIVIDER_RULE}" stroke-width="0.6" />
        <text x="10" y="18.5" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">{html.escape(item.publication)}</text>
        <circle cx="226" cy="14.0" r="4.5" fill="{item.cancer_color}" />
        <text x="238" y="18.5" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">{html.escape(item.cancer_type)}</text>
        <text x="382" y="18.5" font-size="12.0" font-weight="500" fill="{COLOR_TEXT_SECONDARY}">{html.escape(item.regimen)}</text>
        <text x="662" y="18.5" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">{item.cells_raw:,}</text>
        <g transform="translate(684, 4.0)">
          <rect x="0" y="0" width="{w_r}" height="20" fill="{OKABE_BLUISH_GREEN}" />
          <rect x="{w_r}" y="0" width="{w_nr}" height="20" fill="{OKABE_VERMILION}" />
          <text x="{r_center}" y="14.5" font-size="11.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">{item.n_resp}R</text>
          <text x="{nr_center}" y="14.5" font-size="11.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">{item.n_non_resp}NR</text>
        </g>
        <text x="1032" y="18.5" font-size="12.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="end">{item.orr_pct:.1f}%</text>
      </g>"""


def render_panel_e_sc_milo_cohorts() -> str:
    """Render Panel (e): Single-Cell ICB Response Cohorts (Milo Compendium, 9 Cohorts)."""
    sc_total_cells = sum(c.cells_raw for c in SC_MILO_COHORTS)
    sc_total_r = sum(c.n_resp for c in SC_MILO_COHORTS)
    sc_total_nr = sum(c.n_non_resp for c in SC_MILO_COHORTS)
    sc_total_eval = sc_total_r + sc_total_nr
    sc_orr = (sc_total_r / sc_total_eval) * 100.0

    row_h = 28
    table_w = 1044
    bar_w = 240
    sc_w_r = int(round((sc_total_r / sc_total_eval) * bar_w))
    sc_w_nr = bar_w - sc_w_r
    summary_y = 28 + len(SC_MILO_COHORTS) * row_h

    rows_svg = "".join(render_sc_milo_cohort_row(idx, item) for idx, item in enumerate(SC_MILO_COHORTS))

    return f"""
  <!-- ====================================================================== -->
  <!-- LAYER 06: PANEL (E) — SINGLE-CELL MILO COMPENDIUM (9 COHORTS)          -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-06-panel-e" inkscape:label="06_Panel_E_Single_Cell_Milo_Cohorts">
    <g id="pe-milo-table" transform="translate(18, 858)">
      <!-- Outer Table Border -->
      <rect x="0" y="0" width="{table_w}" height="{summary_y + row_h}" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.9" />

      <!-- Header Bar -->
      <rect x="0" y="0" width="{table_w}" height="{row_h}" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.8" />
      <text x="10" y="19.0" font-size="13.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Primary Study</text>
      <text x="238" y="19.0" font-size="13.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Cancer</text>
      <text x="382" y="19.0" font-size="13.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Regimen</text>
      <text x="662" y="19.0" font-size="13.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">Cell Number</text>
      <text x="684" y="19.0" font-size="13.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Response (R / NR)</text>
      <text x="1032" y="19.0" font-size="13.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">ORR</text>

      {rows_svg}

      <!-- Summary Row: Milo Total (9 Cohorts) -->
      <g id="pe-sc-row-summary" transform="translate(0, {summary_y})">
        <rect x="0" y="0" width="{table_w}" height="{row_h}" fill="{COLOR_SUMMARY_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.8" />
        <text x="10" y="18.5" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">Milo Total (9 Cohorts)</text>
        <text x="238" y="18.5" font-size="12.5" font-weight="600" fill="{COLOR_TEXT_SECONDARY}">6 Solid Types</text>
        <text x="382" y="18.5" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">PD-(L)1 ± Combo</text>
        <text x="662" y="18.5" font-size="12.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}" text-anchor="end">{sc_total_cells:,}</text>
        <g transform="translate(684, 4.0)">
          <rect x="0" y="0" width="{sc_w_r}" height="20" fill="{OKABE_BLUISH_GREEN}" />
          <rect x="{sc_w_r}" y="0" width="{sc_w_nr}" height="20" fill="{OKABE_VERMILION}" />
          <text x="{sc_w_r // 2}" y="14.5" font-size="11.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">{sc_total_r}R</text>
          <text x="{sc_w_r + (sc_w_nr // 2)}" y="14.5" font-size="11.5" font-weight="700" fill="#FFFFFF" text-anchor="middle">{sc_total_nr}NR</text>
        </g>
        <text x="1032" y="18.5" font-size="12.5" font-weight="700" fill="{OKABE_BLUISH_GREEN}" text-anchor="end">{sc_orr:.1f}%</text>
      </g>
    </g>
  </g>"""


# ==============================================================================
# MASTER SVG ASSEMBLY & VALIDATION
# ==============================================================================
def render_complete_svg(bundle: AssetBundle) -> str:
    """Assemble the complete 180mm x 230mm 5-panel SVG document with native Inkscape layers."""
    return f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
     width="{SVG_WIDTH_MM}" height="{SVG_HEIGHT_MM}"
     viewBox="0 0 {VIEWBOX_W} {VIEWBOX_H}"
     version="1.1" id="figure-study-workflow-and-cohorts"
     style="background-color: {COLOR_CANVAS_BG}; font-family: {FONT_SANS};">
  <sodipodi:namedview
     id="namedview-study-workflow"
     pagecolor="#ffffff"
     bordercolor="#666666"
     borderopacity="1.0"
     inkscape:document-units="mm" />
  <defs id="defs-main">
    <clipPath id="clip-pa-tme">
      <rect x="17.0" y="63.0" width="222.0" height="258.0" />
    </clipPath>
    {bundle.defs_svg}
  </defs>

  <!-- ====================================================================== -->
  <!-- LAYER 00: CANVAS BACKGROUND & OUTER PANEL WIREFRAME CONTAINERS         -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-00-background" inkscape:label="00_Canvas_Background">
    <rect id="bg-canvas" x="0" y="0" width="{VIEWBOX_W}" height="{VIEWBOX_H}" fill="{COLOR_CANVAS_BG}" />
    <!-- Row 1: Panels (a), (b), (c) -->
    <rect id="card-panel-a" x="8" y="8" width="452" height="360" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1.0" />
    <rect id="card-panel-b" x="468" y="8" width="240" height="360" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1.0" />
    <rect id="card-panel-c" x="716" y="8" width="356" height="360" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1.0" />
    <!-- Row 2: Panel (d) Bulk Clinical Cohorts -->
    <rect id="card-panel-d" x="8" y="378" width="1064" height="422" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1.0" />
    <!-- Row 3: Panel (e) Single-Cell Milo Cohorts -->
    <rect id="card-panel-e" x="8" y="810" width="1064" height="398" fill="{COLOR_PANEL_BG}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="1.0" />
  </g>

  <!-- ====================================================================== -->
  <!-- LAYER 01: PANEL SECTION HEADERS & BADGES                               -->
  <!-- ====================================================================== -->
  <g inkscape:groupmode="layer" id="layer-01-headers" inkscape:label="01_Panel_Headers">
    <!-- Panel (a) Banner -->
    <g id="hdr-panel-a">
      <rect x="8" y="8" width="452" height="32" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.8" />
      <text x="18" y="29" font-size="14.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">(a) TUMOR MICROENVIRONMENT ECOSYSTEM</text>
    </g>
    <!-- Panel (b) Banner -->
    <g id="hdr-panel-b">
      <rect x="468" y="8" width="240" height="32" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.8" />
      <text x="478" y="29" font-size="14.0" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">(b) SPECIMEN &amp; SEQUENCING</text>
    </g>
    <!-- Panel (c) Banner -->
    <g id="hdr-panel-c">
      <rect x="716" y="8" width="356" height="32" fill="{COLOR_SUBTLE_FILL}" stroke="{COLOR_BORDER_HAIRLINE}" stroke-width="0.8" />
      <text x="724" y="29" font-size="13.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">(c) DECONVOLUTION &amp; CLASSIFICATION PIPELINE</text>
    </g>
    <!-- Panel (d) Title -->
    <g id="hdr-panel-d">
      <text x="18" y="408" font-size="15.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">(d) Bulk RNA-seq ICB Clinical Cohorts: Regimens &amp; RECIST Response Rates (9 Cohorts, n = 1,097)</text>
    </g>
    <!-- Panel (e) Title -->
    <g id="hdr-panel-e">
      <text x="18" y="840" font-size="15.5" font-weight="700" fill="{COLOR_TEXT_PRIMARY}">(e) Single-Cell ICB Response Cohorts: Milo Compendium (9 Cohorts, n = 613, 3,597,092 Cells)</text>
    </g>
  </g>
  {render_panel_a_tme_ecosystem(bundle)}
  {render_panel_b_wetlab_workflow()}
  {render_panel_c_deconv_pipeline(bundle)}
  {render_panel_d_bulk_cohorts()}
  {render_panel_e_sc_milo_cohorts()}
</svg>
"""


def validate_svg(svg_text: str) -> Result[int, str]:
    """Validate XML well-formedness, Inkscape layers, unique IDs, and minimum font size."""
    try:
        root = ET.fromstring(svg_text)
    except ET.ParseError as exc:
        return Failure(f"Invalid SVG XML: {exc}")

    ns_inkscape = "{http://www.inkscape.org/namespaces/inkscape}groupmode"
    layers = [el for el in root.iter() if el.attrib.get(ns_inkscape) == "layer"]
    if len(layers) < 7:
        return Failure(f"Expected 7 Inkscape layers, found {len(layers)}")

    seen_ids: set[str] = set()
    dup_ids: list[str] = []
    for el in root.iter():
        eid = el.attrib.get("id", "")
        if eid:
            if eid in seen_ids:
                dup_ids.append(eid)
            seen_ids.add(eid)
        fs_raw = el.attrib.get("font-size", "")
        if fs_raw:
            fs_val = float(re.sub(r"[a-zA-Z%]+$", "", fs_raw))
            if fs_val < MIN_FONT_SIZE_PX:
                return Failure(f"Font size {fs_raw} on <{el.tag}> is below minimum {MIN_FONT_SIZE_PX}")

    return Failure(f"Duplicate XML IDs: {dup_ids[:10]}") if dup_ids else Success(len(layers))


def export_figure(svg_text: str, svg_targets: Sequence[Path], png_targets: Sequence[Path]) -> Result[int, str]:
    """Write validated SVG to disk and rasterize 300 DPI PNGs via plotting_utils."""
    from plotting_utils import rasterize_svg_to_png

    for svg_p in svg_targets:
        svg_p.parent.mkdir(parents=True, exist_ok=True)
        svg_p.write_text(svg_text, encoding="utf-8")

    for png_p in png_targets:
        png_p.parent.mkdir(parents=True, exist_ok=True)
        res = rasterize_svg_to_png(
            svg_path=svg_targets[0],
            png_path=png_p,
            width_px=2160,
            height_px=2432,
            background="white",
        )
        if not is_successful(res):
            return Failure(res.failure())

    return Success(len(svg_targets) + len(png_targets))


def main() -> int:
    """Imperative shell entry point."""
    assets_dir = REPO_ROOT / "assets"
    svg_targets: tuple[Path, ...] = (
        REPO_ROOT / "article" / "figures" / "dataset_overview" / "figure_study_workflow_and_cohorts.svg",
        REPO_ROOT / "assets" / "figure_study_workflow_and_cohorts.svg",
    )
    png_targets: tuple[Path, ...] = (
        REPO_ROOT / "article" / "figures" / "dataset_overview" / "figure_study_workflow_and_cohorts.png",
        REPO_ROOT / "assets" / "rendered_png" / "figure_study_workflow_and_cohorts.png",
    )

    bundle_res = load_asset_bundle(assets_dir)
    if not is_successful(bundle_res):
        sys.stderr.write(f"[ASSET ERROR] {bundle_res.failure()}\n")
        return 1

    svg_content = render_complete_svg(bundle_res.unwrap())
    val_res = validate_svg(svg_content)
    if not is_successful(val_res):
        sys.stderr.write(f"[VALIDATION ERROR] {val_res.failure()}\n")
        return 1

    n_layers = val_res.unwrap()
    sys.stdout.write(
        f"[OK] Validated SVG ({SVG_WIDTH_MM} x {SVG_HEIGHT_MM}) with {n_layers} Inkscape layers.\n"
    )

    exp_res = export_figure(svg_content, svg_targets, png_targets)
    if not is_successful(exp_res):
        sys.stderr.write(f"[EXPORT ERROR] {exp_res.failure()}\n")
        return 1

    for p in (*svg_targets, *png_targets):
        sys.stdout.write(f"[SAVED] {p.relative_to(REPO_ROOT)} ({p.stat().st_size:,} bytes)\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
