#!/usr/bin/env python3
"""Cross-cohort Milo summary figures and single-panel vector crops for the main article.

Produces four manuscript figure assets (all SVG, Okabe-Ito, Nature Methods theme):

1. ``milo_cross_cohort_volcano_summary``: neighborhood-level summary of the nine Milo
   volcano panels (% significant neighborhoods per direction; effect sizes of the
   significant neighborhoods).
2. ``milo_cross_cohort_umap_summary``: cell-level summary of the nine logFC UMAP panels
   (cell fractions per DA status; cell-level log2FC distribution).
3. ``meta42_observed_recist_radar_panel_b``: panel (b) of the Consensus-Meta42 radar figure.
4. ``pan_cohort_harmony_i2_umap_panel_b``: panel (b) of the Harmony integration audit.

Functional core (pure functions returning ``Result``), imperative shell in ``main``.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Final, cast

import altair as alt
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

from plotting_utils import (
    COHORT_PANELS,
    COLOR_HAIRLINE,
    COLOR_TEXT_PRIMARY,
    OKABE_BLUE,
    OKABE_GREY,
    OKABE_VERMILION,
    apply_nature_methods_theme,
    export_altair_figure,
)

alt.data_transformers.disable_max_rows()

R_LABEL: Final[str] = "Enriched in responders"
NR_LABEL: Final[str] = "Enriched in non-responders"
NS_LABEL: Final[str] = "Not significant"
STATUS_DOMAIN: Final[tuple[str, ...]] = (R_LABEL, NS_LABEL, NR_LABEL)
STATUS_COLORS: Final[tuple[str, ...]] = (OKABE_VERMILION, OKABE_GREY, OKABE_BLUE)


class CropBox(BaseModel):
    """Crop rectangle in source SVG user units."""

    model_config = ConfigDict(frozen=True)
    x: float
    y: float
    w: float
    h: float


class SummaryConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    results_dir: Path = Path("output/manuscript_pipeline/2a1_sample_and_analyze_cohorts")
    out_dir: Path = Path("output/manuscript_pipeline/4g_plot_milo_cross_cohort_summaries")
    radar_svg: Path = Path("article/figures/milo-deconv/iatlas_consensus_metastate_responder_radar.svg")
    harmony_svg: Path = Path(
        "output/results/milopy_pan_cohort_integration_and_power/figures/pan_cohort_harmony_integration_audit.svg"
    )
    radar_crop: CropBox = CropBox(x=790.0, y=0.0, w=700.0, h=425.0)
    harmony_crop: CropBox = CropBox(x=425.0, y=30.0, w=635.0, h=395.0)
    fdr_threshold: float = 0.10
    volcano_letters: tuple[str, str] = ("e", "f")
    umap_letters: tuple[str, str] = ("g", "h")
    tag_cohorts: bool = False
    show_title: bool = False


# ---------------------------------------------------------------------------
# Pure core
# ---------------------------------------------------------------------------

def cohort_label(meta: dict[str, str], tagged: bool) -> str:
    """Display label; optionally prefixed with the volcano / UMAP grid panel tag (a-i)."""
    base = f"{meta['author'].replace(' et al.', '')} ({meta['indication']})"
    return f"{meta['tag']}  {base}" if tagged else base


def labels_for(tagged: bool) -> dict[str, str]:
    return {m["acc"]: cohort_label(m, tagged) for m in COHORT_PANELS}


def order_for(tagged: bool) -> tuple[str, ...]:
    labels = labels_for(tagged)
    return tuple(labels[m["acc"]] for m in COHORT_PANELS)


def _is_responder_enriched(status: pl.Expr) -> pl.Expr:
    return status.str.starts_with("DA+") | status.str.starts_with("Enriched in Responders")


def _is_nonresponder_enriched(status: pl.Expr) -> pl.Expr:
    return status.str.starts_with("DA-") | status.str.starts_with("Enriched in Non-Responders")


def harmonize_status(status_col: str) -> pl.Expr:
    """Map the heterogeneous DA status labels to the three canonical categories."""
    status = pl.col(status_col)
    return (
        pl.when(_is_responder_enriched(status))
        .then(pl.lit(R_LABEL))
        .when(_is_nonresponder_enriched(status))
        .then(pl.lit(NR_LABEL))
        .otherwise(pl.lit(NS_LABEL))
        .alias("direction")
    )


def summarize_neighborhoods(da: pl.DataFrame) -> Result[pl.DataFrame, str]:
    """Per-cohort share of significant neighborhoods and effect sizes of the significant ones."""
    return (
        Failure("No neighborhood-level DA rows available.")
        if da.is_empty()
        else Success(
            da.lazy()
            .with_columns(harmonize_status("status"))
            .group_by("cohort")
            .agg(
                pl.len().alias("n_nhoods"),
                (pl.col("direction") == R_LABEL).sum().alias("n_r"),
                (pl.col("direction") == NR_LABEL).sum().alias("n_nr"),
                pl.col("logFC").filter(pl.col("direction") == R_LABEL).median().alias("med_r"),
                pl.col("logFC").filter(pl.col("direction") == R_LABEL).quantile(0.25).alias("q25_r"),
                pl.col("logFC").filter(pl.col("direction") == R_LABEL).quantile(0.75).alias("q75_r"),
                pl.col("logFC").filter(pl.col("direction") == NR_LABEL).median().alias("med_nr"),
                pl.col("logFC").filter(pl.col("direction") == NR_LABEL).quantile(0.25).alias("q25_nr"),
                pl.col("logFC").filter(pl.col("direction") == NR_LABEL).quantile(0.75).alias("q75_nr"),
            )
            .with_columns(
                (100.0 * pl.col("n_r") / pl.col("n_nhoods")).alias("pct_r"),
                (-100.0 * pl.col("n_nr") / pl.col("n_nhoods")).alias("pct_nr"),
            )
            .collect()
        )
    )


def summarize_cells(cells: pl.DataFrame) -> Result[pl.DataFrame, str]:
    """Per-cohort cell fractions by DA direction and cell-level log2FC quantiles."""
    return (
        Failure("No cell-level rows available.")
        if cells.is_empty()
        else Success(
            cells.lazy()
            .with_columns(harmonize_status("da_status"))
            .group_by("cohort")
            .agg(
                pl.len().alias("n_cells"),
                (pl.col("direction") == R_LABEL).sum().alias("n_r"),
                (pl.col("direction") == NR_LABEL).sum().alias("n_nr"),
                pl.col("cell_logfc").quantile(0.05).alias("q05"),
                pl.col("cell_logfc").quantile(0.25).alias("q25"),
                pl.col("cell_logfc").median().alias("med"),
                pl.col("cell_logfc").quantile(0.75).alias("q75"),
                pl.col("cell_logfc").quantile(0.95).alias("q95"),
            )
            .with_columns(
                (100.0 * pl.col("n_r") / pl.col("n_cells")).alias("pct_r"),
                (100.0 * pl.col("n_nr") / pl.col("n_cells")).alias("pct_nr"),
                (100.0 * (pl.col("n_cells") - pl.col("n_r") - pl.col("n_nr")) / pl.col("n_cells")).alias(
                    "pct_ns"
                ),
            )
            .collect()
        )
    )


def with_labels(df: pl.DataFrame, tagged: bool) -> pl.DataFrame:
    return df.with_columns(
        pl.col("cohort").replace_strict(labels_for(tagged), default=pl.col("cohort")).alias("label")
    )


def _y_axis(show_labels: bool, tagged: bool) -> alt.Y:
    return alt.Y(
        "label:N",
        sort=list(order_for(tagged)),
        title=None,
        axis=alt.Axis(labels=show_labels, ticks=show_labels, domain=show_labels, labelLimit=260),
    )


def _legend(columns: int) -> alt.Legend:
    return (
        alt.Legend(title=None, orient="bottom", columns=columns, labelLimit=170)
        if columns > 0
        else alt.Legend(title=None, orient="bottom")
    )


def _status_scale() -> alt.Scale:
    return alt.Scale(domain=list(STATUS_DOMAIN), range=list(STATUS_COLORS))


class PanelSize(BaseModel):
    """Plot-area size (px) of the left and right panel of a summary figure."""

    model_config = ConfigDict(frozen=True)
    left_w: int = 300
    right_w: int = 260
    height: int = 230
    legend_columns: int = 0
    show_left_y: bool = True


def _alt_data(records: list[dict[str, Any]]) -> Any:
    data_fn: Any = alt.Data
    return data_fn(values=records)


def build_neighborhood_summary_panels(
    summary: pl.DataFrame, letters: tuple[str, str], tagged: bool, size: PanelSize = PanelSize()
) -> tuple[alt.LayerChart, alt.LayerChart]:
    data = _alt_data(with_labels(summary, tagged).to_dicts())
    bars = (
        alt.Chart(data)
        .transform_fold(["pct_r", "pct_nr"], as_=["key", "pct"])
        .transform_calculate(
            direction="datum.key === 'pct_r' ? '" + R_LABEL + "' : '" + NR_LABEL + "'"
        )
        .mark_bar(height=15, stroke=COLOR_TEXT_PRIMARY, strokeWidth=0.5)
        .encode(
            y=_y_axis(size.show_left_y, tagged),
            x=alt.X("pct:Q", title="Significant neighborhoods (%)"),
            color=alt.Color("direction:N", scale=_status_scale(), legend=_legend(size.legend_columns)),
        )
        .properties(width=size.left_w, height=size.height, title=f"({letters[0]})  Significant neighborhoods")
    )
    zero = alt.Chart(_alt_data([{"z": 0}])).mark_rule(strokeDash=[3, 3], color=COLOR_TEXT_PRIMARY).encode(x="z:Q")
    labelled_rows = with_labels(summary, tagged).to_dicts()
    effect_rows = [
        row
        for row in labelled_rows
        if row["n_r"] > 0 or row["n_nr"] > 0
    ]
    long_rows = tuple(
        {"label": r["label"], "direction": R_LABEL, "med": r["med_r"], "q25": r["q25_r"], "q75": r["q75_r"], "n": r["n_r"]}
        for r in effect_rows
        if r["n_r"] > 0
    ) + tuple(
        {"label": r["label"], "direction": NR_LABEL, "med": r["med_nr"], "q25": r["q25_nr"], "q75": r["q75_nr"], "n": r["n_nr"]}
        for r in effect_rows
        if r["n_nr"] > 0
    )
    eff = alt.Chart(_alt_data(list(long_rows)))
    ranges = eff.mark_rule(strokeWidth=2).encode(
        y=_y_axis(False, tagged),
        x=alt.X("q25:Q", title="Neighborhood log\u2082FC (med., IQR)"),
        x2="q75:Q",
        color=alt.Color("direction:N", scale=_status_scale(), legend=None),
    )
    dots = eff.mark_point(filled=True, size=65, stroke=COLOR_TEXT_PRIMARY, strokeWidth=0.5).encode(
        y=_y_axis(False, tagged),
        x="med:Q",
        color=alt.Color("direction:N", scale=_status_scale(), legend=None),
    )
    panel_b = (ranges + dots + zero).properties(
        width=size.right_w, height=size.height, title=f"({letters[1]})  Effect sizes (FDR < 0.10)"
    )
    return bars + zero, panel_b


def build_neighborhood_summary_chart(
    summary: pl.DataFrame, letters: tuple[str, str], tagged: bool
) -> alt.HConcatChart:
    left, right = build_neighborhood_summary_panels(summary, letters, tagged)
    return alt.hconcat(left, right, spacing=24)


def build_cell_summary_panels(
    summary: pl.DataFrame, letters: tuple[str, str], tagged: bool, size: PanelSize = PanelSize()
) -> tuple[alt.Chart, alt.LayerChart]:
    labelled = with_labels(summary, tagged)
    stack_rows = tuple(
        {"label": r["label"], "direction": d, "pct": r[k]}
        for r in labelled.to_dicts()
        for d, k in ((R_LABEL, "pct_r"), (NS_LABEL, "pct_ns"), (NR_LABEL, "pct_nr"))
    )
    stacked = (
        alt.Chart(_alt_data(list(stack_rows)))
        .mark_bar(height=15, stroke=COLOR_TEXT_PRIMARY, strokeWidth=0.5)
        .encode(
            y=_y_axis(size.show_left_y, tagged),
            x=alt.X("pct:Q", title="Cells by DA status (%)", scale=alt.Scale(domain=[0, 100])),
            color=alt.Color("direction:N", scale=_status_scale(), legend=_legend(size.legend_columns)),
            order=alt.Order("direction:N", sort="ascending"),
        )
        .properties(width=size.left_w, height=size.height, title=f"({letters[0]})  Cells per DA status")
    )
    base = alt.Chart(_alt_data(labelled.to_dicts()))
    whisker = base.mark_rule(color=COLOR_TEXT_PRIMARY).encode(
        y=_y_axis(False, tagged), x=alt.X("q05:Q", title="Cell log\u2082FC (med., IQR, 5\u201395%)"), x2="q95:Q"
    )
    box = base.mark_bar(height=10, color=OKABE_GREY, stroke=COLOR_TEXT_PRIMARY, strokeWidth=0.5).encode(
        y=_y_axis(False, tagged), x="q25:Q", x2="q75:Q"
    )
    median = base.mark_tick(color=COLOR_TEXT_PRIMARY, thickness=2, size=14).encode(y=_y_axis(False, tagged), x="med:Q")
    zero = alt.Chart(_alt_data([{"z": 0}])).mark_rule(strokeDash=[3, 3], color=COLOR_TEXT_PRIMARY).encode(x="z:Q")
    panel_b = (whisker + box + median + zero).properties(
        width=size.right_w, height=size.height, title=f"({letters[1]})  Cell-level log\u2082FC"
    )
    return stacked, panel_b


def build_cell_summary_chart(
    summary: pl.DataFrame, letters: tuple[str, str], tagged: bool
) -> alt.HConcatChart:
    left, right = build_cell_summary_panels(summary, letters, tagged)
    return alt.hconcat(left, right, spacing=24)


def finalize_chart(chart: alt.HConcatChart, title: str, subtitle: str, show_title: bool) -> alt.HConcatChart:
    titled = (
        chart.properties(
            title=alt.TitleParams(text=title, subtitle=subtitle, anchor="start", fontSize=14, subtitleFontSize=10)
        )
        if show_title
        else chart
    )
    return cast(alt.HConcatChart, apply_nature_methods_theme(titled.configure_view(stroke=COLOR_HAIRLINE, strokeWidth=0.75)))


_ROOT_TAG: Final[re.Pattern[str]] = re.compile(r"<svg\b[^>]*>", re.S)
_ATTR: Final[str] = r'{name}="([-0-9.]+)([a-z%]*)"'


def _attr(tag: str, name: str) -> Result[tuple[float, str], str]:
    found = re.search(_ATTR.format(name=name), tag)
    return Failure(f"Missing root attribute {name}") if found is None else Success((float(found[1]), found[2]))


def crop_svg_viewbox(svg_text: str, box: CropBox) -> Result[str, str]:
    """Return a copy of ``svg_text`` whose visible area is restricted to ``box`` (vector preserved)."""
    root = _ROOT_TAG.search(svg_text)
    if root is None:
        return Failure("No <svg> root element found.")
    tag = root.group(0)
    view = re.search(r'viewBox="([-0-9.]+) ([-0-9.]+) ([-0-9.]+) ([-0-9.]+)"', tag)
    width = _attr(tag, "width")
    height = _attr(tag, "height")
    match (view, width, height):
        case (re.Match() as vb, Success((w_val, w_unit)), Success((h_val, h_unit))):
            scale_x = w_val / float(vb[3])
            scale_y = h_val / float(vb[4])
            new_tag = re.sub(r'viewBox="[^"]*"', f'viewBox="{box.x} {box.y} {box.w} {box.h}"', tag)
            new_tag = re.sub(r'width="[^"]*"', f'width="{round(box.w * scale_x, 3)}{w_unit}"', new_tag, count=1)
            new_tag = re.sub(r'height="[^"]*"', f'height="{round(box.h * scale_y, 3)}{h_unit}"', new_tag, count=1)
            return Success(svg_text[: root.start()] + new_tag + svg_text[root.end():])
        case _:
            return Failure("Root <svg> lacks viewBox/width/height; cannot crop.")


# ---------------------------------------------------------------------------
# Imperative shell
# ---------------------------------------------------------------------------

def load_cohort_tables(results_dir: Path) -> Result[tuple[pl.DataFrame, pl.DataFrame], str]:
    da_paths = tuple((m["acc"], results_dir / m["acc"] / "da_results.parquet") for m in COHORT_PANELS)
    cell_paths = tuple((m["acc"], results_dir / m["acc"] / "cell_level_scores.parquet") for m in COHORT_PANELS)
    missing = tuple(str(p) for _, p in da_paths + cell_paths if not p.exists())
    return (
        Failure("Missing inputs: " + ", ".join(missing))
        if missing
        else Success(
            (
                pl.concat(
                    (pl.read_parquet(p).with_columns(pl.lit(a).alias("cohort")) for a, p in da_paths),
                    how="diagonal_relaxed",
                ),
                pl.concat(
                    (
                        pl.read_parquet(p, columns=["cell_logfc", "da_status"]).with_columns(pl.lit(a).alias("cohort"))
                        for a, p in cell_paths
                    ),
                    how="diagonal_relaxed",
                ),
            )
        )
    )


def _write_text(dst: Path, text: str) -> Path:
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(text, encoding="utf-8")
    return dst


def write_crop(src: Path, dst: Path, box: CropBox) -> Result[Path, str]:
    return (
        Failure(f"Missing source SVG: {src}")
        if not src.exists()
        else crop_svg_viewbox(src.read_text(encoding="utf-8"), box).map(lambda text: _write_text(dst, text))
    )


def run(config: SummaryConfig) -> Result[Sequence[Path], str]:
    config.out_dir.mkdir(parents=True, exist_ok=True)
    tables = load_cohort_tables(config.results_dir)
    nhood = tables.bind(lambda t: summarize_neighborhoods(t[0]))
    cells = tables.bind(lambda t: summarize_cells(t[1]))
    outputs = (
        nhood.map(
            lambda s: export_altair_figure(
                finalize_chart(
                    build_neighborhood_summary_chart(s, config.volcano_letters, config.tag_cohorts),
                    "Single-cell differential abundance summary across 9 ICB cohorts",
                    "Milo kNN neighborhoods, quasi-likelihood NB GLM, FDR < 0.10; responder vs. non-responder",
                    config.show_title,
                ),
                config.out_dir / "milo_cross_cohort_volcano_summary",
            )[0]
        ),
        cells.map(
            lambda s: export_altair_figure(
                finalize_chart(
                    build_cell_summary_chart(s, config.umap_letters, config.tag_cohorts),
                    "Single-cell landscape summary across 9 ICB cohorts",
                    "Cell-level Milo log2FC (responder vs. non-responder) and neighborhood DA status",
                    config.show_title,
                ),
                config.out_dir / "milo_cross_cohort_umap_summary",
            )[0]
        ),
        write_crop(config.radar_svg, config.out_dir / "meta42_observed_recist_radar_panel_b.svg", config.radar_crop),
        write_crop(config.harmony_svg, config.out_dir / "pan_cohort_harmony_i2_umap_panel_b.svg", config.harmony_crop),
    )
    failures = tuple(o.failure() for o in outputs if isinstance(o, Failure))
    return (
        Failure("; ".join(failures))
        if failures
        else Success(tuple(o.unwrap() for o in outputs))
    )


def parse_args(argv: Sequence[str]) -> SummaryConfig:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=SummaryConfig().results_dir)
    parser.add_argument("--out-dir", type=Path, default=SummaryConfig().out_dir)
    ns = parser.parse_args(argv)
    return SummaryConfig(results_dir=ns.results_dir, out_dir=ns.out_dir)


def main() -> None:
    match run(parse_args(sys.argv[1:])):
        case Success(paths):
            print("\n".join(f"wrote {p}" for p in paths))
        case Failure(err):
            print(f"4g failed: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
