#!/usr/bin/env python3
"""Standalone figure exporter for Pan-Cancer T-Cell, CAF & Myeloid Triad comparison.

Reads cached Polars `.parquet` outputs from `output/results/pan_cancer_atlas_comparison/`
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

from plotting_utils import export_altair_figure  # noqa: E402

from scripts.pan_cancer_atlas_comparison.models import (  # noqa: E402
    AltairChart,
    TriadCompartment,
    VectorChartExporter,
)
from scripts.pan_cancer_atlas_comparison.plotting import (  # noqa: E402
    build_atlas_similarity_and_ici_chart,
    build_triad_all42_specificity_chart,
)

_SVG_NS: Final[str] = "http://www.w3.org/2000/svg"
_XLINK_NS: Final[str] = "http://www.w3.org/1999/xlink"
_INKSCAPE_NS: Final[str] = "http://www.inkscape.org/namespaces/inkscape"
_SODIPODI_NS: Final[str] = "http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"


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
            max_vector_points=3500,
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


def export_all_pan_atlas_figures(
    full_sim_df: pl.DataFrame,
    nearest_tcell_df: pl.DataFrame,
    nearest_caf_df: pl.DataFrame,
    atlas_ici_df: pl.DataFrame,
    triad_summary_df: pl.DataFrame,
    figures_dir: Path,
    exporter: VectorChartExporter = default_vector_exporter,
) -> Result[tuple[tuple[Path, Path], ...], str]:
    """Build and export one figure per T-cell/CAF atlas plus the Triad figure (`.svg` + 300 DPI `.png`)."""
    figures_dir.mkdir(parents=True, exist_ok=True)
    meta_beta_df = pl.concat(
        [
            d.select(["meta_state_label", "meta42_bulk_beta"])
            for d in (nearest_tcell_df, nearest_caf_df)
            if not d.is_empty()
        ]
    ) if (not nearest_tcell_df.is_empty() or not nearest_caf_df.is_empty()) else pl.DataFrame()
    lineage_atlas_ids = (
        full_sim_df.filter(
            pl.col("atlas_compartment").is_in(
                [TriadCompartment.T_NK.value, TriadCompartment.STROMAL_CAF.value]
            )
        )["atlas_id"]
        .unique(maintain_order=True)
        .to_list()
    )
    specs: tuple[tuple[str, Result[AltairChart, str]], ...] = tuple(
        (
            f"meta42_{aid}_similarity_and_ici",
            build_atlas_similarity_and_ici_chart(
                atlas_id=str(aid),
                full_sim_df=full_sim_df,
                atlas_ici_df=atlas_ici_df,
                meta_beta_df=meta_beta_df,
            ),
        )
        for aid in lineage_atlas_ids
    ) + (
        (
            "meta42_triad_all42_external_specificity",
            build_triad_all42_specificity_chart(
                full_sim_df=full_sim_df,
                triad_summary_df=triad_summary_df,
            ),
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


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    """Parse CLI arguments for standalone figure rendering."""
    parser = argparse.ArgumentParser(
        description="Render publication SVG/PNG figures from cached Pan-Cancer Triad atlas parquets."
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("output/results/pan_cancer_atlas_comparison"),
        help="Directory containing generated comparison parquet tables.",
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=Path("output/figures/pan_cancer_atlas_comparison"),
        help="Destination directory for publication SVG and 300 DPI PNG figures.",
    )
    return parser.parse_args(list(argv))


def main(argv: Maybe[Sequence[str]] = Nothing) -> int:
    """Imperative CLI entrypoint for standalone figure regeneration."""
    raw_argv: Sequence[str] = argv.unwrap() if isinstance(argv, Some) else sys.argv[1:]
    args = parse_args(raw_argv)
    res_dir: Path = args.results_dir if args.results_dir.is_absolute() else (_REPO_ROOT / args.results_dir)
    fig_dir: Path = args.figures_dir if args.figures_dir.is_absolute() else (_REPO_ROOT / args.figures_dir)

    sim_p = res_dir / "meta42_vs_external_atlases_similarity.parquet"
    nt_p = res_dir / "meta42_tcell_nearest_substates.parquet"
    nc_p = res_dir / "meta42_caf_nearest_substates.parquet"
    ici_p = res_dir / "external_atlas_substate_ici_outcomes.parquet"
    tri_p = res_dir / "triad_block_specificity_summary.parquet"

    for p in (sim_p, nt_p, nc_p, ici_p, tri_p):
        if not p.is_file():
            print(f"[ERROR] Missing required parquet file: {p}", file=sys.stderr)
            return 1

    match export_all_pan_atlas_figures(
        full_sim_df=pl.read_parquet(sim_p),
        nearest_tcell_df=pl.read_parquet(nt_p),
        nearest_caf_df=pl.read_parquet(nc_p),
        atlas_ici_df=pl.read_parquet(ici_p),
        triad_summary_df=pl.read_parquet(tri_p),
        figures_dir=fig_dir,
    ):
        case Failure(err):
            print(f"[ERROR] Figure export failed: {err}", file=sys.stderr)
            return 1
        case Success(paths):
            for svg_p, png_p in paths:
                print(f"[OK] Exported {svg_p} & {png_p}")
            return 0
        case _:
            return 1


if __name__ == "__main__":
    raise SystemExit(main())
