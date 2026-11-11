#!/usr/bin/env python3
"""Standalone figure exporter for the Meta42 <-> Coulton_TAM_Atlas comparison suite.

Reads cached Polars `.parquet` outputs from `output/results/coulton_tam_comparison/`
and exports publication-grade Altair vector `.svg` (with native Inkscape layer annotations)
and 300 DPI verification `.png` figures via `plotting_utils`.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path
import sys
from typing import Final
import xml.etree.ElementTree as ET

import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

_REPO_ROOT: Final[Path] = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from plotting_utils import export_altair_figure

from scripts.coulton_tam_comparison.models import AltairChart, VectorChartExporter
from scripts.coulton_tam_comparison.plotting import (
    build_all42_pan_lineage_coulton_chart,
    build_myeloid_coulton_similarity_and_ici_chart,
    build_nearest_substate_outcome_bridge_chart,
    build_simplified_bicluster_typst_caption,
    build_simplified_biclustered_myeloid_coulton_chart,
)

_SVG_NS: Final[str] = "http://www.w3.org/2000/svg"
_XLINK_NS: Final[str] = "http://www.w3.org/1999/xlink"
_INKSCAPE_NS: Final[str] = "http://www.inkscape.org/namespaces/inkscape"
_SODIPODI_NS: Final[str] = "http://sodipodi.sourceforge.org/DTD/sodipodi-0.dtd"


def inject_inkscape_layers_pure(svg_text: str, figure_id: str) -> Result[str, str]:
    """Purely enrich a compiled Vega-Lite SVG with explicit Inkscape XML namespaces and layers."""
    ET.register_namespace("", _SVG_NS)
    ET.register_namespace("xlink", _XLINK_NS)
    ET.register_namespace("inkscape", _INKSCAPE_NS)
    ET.register_namespace("sodipodi", _SODIPODI_NS)

    try:
        root = ET.fromstring(svg_text)
    except ET.ParseError as exc:
        return Failure(f"Invalid SVG XML generated for {figure_id}: {exc}")

    root.set("id", figure_id)

    children = list(root)
    for ch in children:
        root.remove(ch)

    defs_elements = [c for c in children if c.tag.endswith("defs") or c.tag.endswith("style")]
    body_elements = [c for c in children if c not in defs_elements]

    for d_el in defs_elements:
        root.append(d_el)

    bg_layer = ET.Element(
        f"{{{_SVG_NS}}}g",
        {
            f"{{{_INKSCAPE_NS}}}groupmode": "layer",
            "id": "layer-00-background",
            f"{{{_INKSCAPE_NS}}}label": "00_Canvas_Background",
        },
    )
    bg_rect_attrs: dict[str, str] = {
        "id": f"bg-{figure_id}",
        "x": "0",
        "y": "0",
        "width": root.attrib.get("width", "100%"),
        "height": root.attrib.get("height", "100%"),
        "rx": "0",
        "fill": "#FFFFFF",
    }
    bg_layer.append(ET.Element(f"{{{_SVG_NS}}}rect", bg_rect_attrs))
    root.append(bg_layer)

    panel_layer = ET.Element(
        f"{{{_SVG_NS}}}g",
        {
            f"{{{_INKSCAPE_NS}}}groupmode": "layer",
            "id": "layer-01-panels",
            f"{{{_INKSCAPE_NS}}}label": "01_Data_Panels_And_Axes",
        },
    )
    for b_el in body_elements:
        panel_layer.append(b_el)
    root.append(panel_layer)

    serialized = ET.tostring(root, encoding="unicode")
    # Ensure explicit namespace attributes remain visible on root even when unused in attributes
    if "xmlns:sodipodi=" not in serialized:
        serialized = serialized.replace(
            "<svg ",
            f'<svg xmlns:sodipodi="{_SODIPODI_NS}" ',
            1,
        )
    if "xmlns:inkscape=" not in serialized:
        serialized = serialized.replace(
            "<svg ",
            f'<svg xmlns:inkscape="{_INKSCAPE_NS}" ',
            1,
        )
    return Success(serialized)


def default_vector_exporter(chart: AltairChart, base_path: Path) -> Result[tuple[Path, Path], str]:
    """Imperative exporter compiling Altair to Inkscape-layered SVG and 300 DPI PNG."""
    try:
        svg_path, png_path = export_altair_figure(
            chart=chart,
            base_path=base_path,
            scale=2.5,
            max_vector_points=2500,
        )
        raw_svg = svg_path.read_text(encoding="utf-8")
        match inject_inkscape_layers_pure(raw_svg, figure_id=base_path.stem):
            case Success(layered_svg):
                svg_path.write_text(layered_svg, encoding="utf-8")
            case Failure(err):
                return Failure(err)
        return Success((svg_path, png_path))
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to export figure {base_path.name}: {exc}")


def export_all_comparison_figures(
    full_sim_df: pl.DataFrame,
    nearest_df: pl.DataFrame,
    coulton_ici_df: pl.DataFrame,
    figures_dir: Path,
    exporter: VectorChartExporter = default_vector_exporter,
) -> Result[tuple[tuple[Path, Path], ...], str]:
    """Build and export all four publication vector figures."""
    figures_dir.mkdir(parents=True, exist_ok=True)
    specs: tuple[tuple[str, Result[AltairChart, str]], ...] = (
        (
            "fig1_meta42_myeloid_coulton_similarity_and_ici",
            build_myeloid_coulton_similarity_and_ici_chart(
                full_sim_df=full_sim_df,
                nearest_df=nearest_df,
                coulton_ici_df=coulton_ici_df,
            ),
        ),
        (
            "fig1b_meta42_myeloid_coulton_simplified_bicluster",
            build_simplified_biclustered_myeloid_coulton_chart(
                full_sim_df=full_sim_df,
                nearest_df=nearest_df,
                coulton_ici_df=coulton_ici_df,
            ),
        ),
        (
            "fig2_macrophage_nearest_substate_ici_forest",
            build_nearest_substate_outcome_bridge_chart(nearest_df=nearest_df),
        ),
        (
            "fig3_all42_pan_lineage_coulton_specificity",
            build_all42_pan_lineage_coulton_chart(full_sim_df=full_sim_df),
        ),
    )

    exported: list[tuple[Path, Path]] = []
    for stem, chart_res in specs:
        match chart_res:
            case Failure(c_err):
                return Failure(c_err)
            case Success(chart):
                match exporter(chart, figures_dir / stem):
                    case Failure(e_err):
                        return Failure(e_err)
                    case Success(paths):
                        exported.append(paths)
    return Success(tuple(exported))


def export_main_paper_typst_caption(
    full_sim_df: pl.DataFrame,
    nearest_df: pl.DataFrame,
    coulton_ici_df: pl.DataFrame,
    results_dir: Path,
) -> Result[Path, str]:
    """Export the standalone importable Typst caption module for the heading-free bicluster figure."""
    match build_simplified_bicluster_typst_caption(
        full_sim_df=full_sim_df,
        nearest_df=nearest_df,
        coulton_ici_df=coulton_ici_df,
    ):
        case Failure(c_err):
            return Failure(c_err)
        case Success(typst_text):
            try:
                results_dir.mkdir(parents=True, exist_ok=True)
                cap_path = results_dir / "coulton_tam_simplified_bicluster_main_caption.typ"
                cap_path.write_text(typst_text, encoding="utf-8")
                return Success(cap_path)
            except Exception as exc:  # noqa: BLE001
                return Failure(f"Failed to write Typst caption file: {exc}")
        case _:
            return Failure("Unexpected Result state in build_simplified_bicluster_typst_caption")


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    """Parse CLI arguments for standalone figure rendering."""
    parser = argparse.ArgumentParser(
        description="Render publication SVG/PNG figures from cached Meta42 <-> Coulton_TAM_Atlas parquets."
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("output/results/coulton_tam_comparison"),
        help="Directory containing generated comparison parquet tables.",
    )
    return parser.parse_args(list(argv))


def main(argv: Maybe[Sequence[str]] = Nothing) -> int:
    """Imperative CLI entrypoint for standalone figure regeneration."""
    raw_argv: Sequence[str] = argv.unwrap() if isinstance(argv, Some) else sys.argv[1:]
    args = parse_args(raw_argv)
    res_dir: Path = args.results_dir

    sim_path = res_dir / "meta42_vs_coulton_full_similarity.parquet"
    near_path = res_dir / "meta42_macrophage_nearest_coulton_substates.parquet"
    ici_path = res_dir / "coulton_substate_ici_outcomes.parquet"

    for p in (sim_path, near_path, ici_path):
        if not p.is_file():
            print(f"[ERROR] Missing required parquet file: {p}", file=sys.stderr)
            return 1

    full_sim_df = pl.read_parquet(sim_path)
    nearest_df = pl.read_parquet(near_path)
    coulton_ici_df = pl.read_parquet(ici_path)

    export_res = export_all_comparison_figures(
        full_sim_df=full_sim_df,
        nearest_df=nearest_df,
        coulton_ici_df=coulton_ici_df,
        figures_dir=res_dir / "figures",
    )
    match export_res:
        case Failure(err):
            print(f"[ERROR] Figure export failed: {err}", file=sys.stderr)
            return 1
        case Success(paths):
            for svg_p, png_p in paths:
                print(f"[OK] Exported {svg_p} & {png_p}")
        case _:
            return 1

    match export_main_paper_typst_caption(
        full_sim_df=full_sim_df,
        nearest_df=nearest_df,
        coulton_ici_df=coulton_ici_df,
        results_dir=res_dir,
    ):
        case Failure(cap_err):
            print(f"[ERROR] Typst caption export failed: {cap_err}", file=sys.stderr)
            return 1
        case Success(cap_p):
            print(f"[OK] Exported {cap_p}")
            return 0
        case _:
            return 1


if __name__ == "__main__":
    raise SystemExit(main())
