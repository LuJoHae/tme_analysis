"""Declarative UMAP plotting script for GSE97168 using Altair."""

import argparse
import sys
from pathlib import Path
from typing import Final, assert_never

import altair as alt  # type: ignore
import anndata as ad  # type: ignore
import polars as pl  # type: ignore
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

alt.data_transformers.disable_max_rows()

EXCLUDE_FIELDS: Final[set[str]] = {"cell_id", "n_genes", "n_counts", "plate-row", "plate-col"}


class PlotConfig(BaseModel):
    """Immutable configuration for GSE97168 UMAP plotting."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    input_h5ad: Path
    out_svg: Path
    out_png: Path


def load_anndata(h5ad_path: Path) -> Result[ad.AnnData, str]:
    """Load AnnData from HDF5 file."""
    try:
        adata = ad.read_h5ad(h5ad_path)
        return Success(adata)
    except Exception as err:
        return Failure(f"Failed to load AnnData from {h5ad_path}: {err}")


def build_umap_chart(adata: ad.AnnData) -> Result[alt.Chart, str]:
    """Construct multi-panel UMAP scatter plot using Altair with Polars data processing."""
    try:
        if "X_umap" not in adata.obsm:
            return Failure("AnnData object does not contain UMAP coordinates ('X_umap'). Run preprocessing first.")

        umap_coords = adata.obsm["X_umap"]

        coords_df = pl.DataFrame(
            {
                "cell_id": [str(c) for c in adata.obs_names],
                "UMAP 1": umap_coords[:, 0],
                "UMAP 2": umap_coords[:, 1],
            }
        )

        obs_pd = adata.obs.reset_index()
        obs_df = pl.from_pandas(obs_pd)
        id_col = obs_df.columns[0]
        obs_df = obs_df.rename({id_col: "cell_id"})

        plot_df = coords_df.join(obs_df, on="cell_id", how="left")

        if plot_df.height > 30000:
            print(f"Downsampling {plot_df.height} cells to 30,000 for vector visualization rendering...")
            plot_df = plot_df.sample(n=30000, seed=42)

        # Select all informative metadata columns for plotting
        target_fields: list[str] = []
        priority_keys = ["leiden", "patient", "patient_id", "tissue", "condition", "sample", "cell_type", "cluster"]

        for key in priority_keys:
            matched_cols = [c for c in plot_df.columns if key in c.lower() and c not in EXCLUDE_FIELDS]
            for col in matched_cols:
                if col not in target_fields and plot_df[col].n_unique() > 1:
                    target_fields.append(col)

        # Fallback to general obs columns if priority keys not matched
        if not target_fields:
            for col in adata.obs.columns:
                if col not in EXCLUDE_FIELDS and plot_df[col].n_unique() > 1:
                    target_fields.append(col)

        # Cap max panels to 6 for readable grid
        target_fields = target_fields[:6]

        charts: list[alt.Chart] = []
        for field in target_fields:
            field_df = plot_df.with_columns(
                pl.when(pl.col(field).is_null() | (pl.col(field) == ""))
                .then(pl.lit("Unannotated"))
                .otherwise(pl.col(field))
                .alias(field)
            )

            chart = (
                alt.Chart(field_df)
                .mark_circle(size=6, opacity=0.75)
                .encode(
                    x=alt.X("UMAP 1:Q", title="UMAP 1", scale=alt.Scale(zero=False)),
                    y=alt.Y("UMAP 2:Q", title="UMAP 2", scale=alt.Scale(zero=False)),
                    color=alt.Color(
                        f"{field}:N",
                        title=field.replace("_", " ").title(),
                        legend=alt.Legend(symbolLimit=25),
                    ),
                    tooltip=[
                        "cell_id",
                        alt.Tooltip(f"{field}:N", title=field.replace("_", " ").title()),
                    ],
                )
                .properties(
                    title=f"GSE97168 UMAP: {field.replace('_', ' ').title()}",
                    width=380,
                    height=380,
                )
            )
            charts.append(chart)

        if not charts:
            return Failure("No valid metadata fields found for plotting UMAPs.")

        concat_chart = alt.concat(*charts, columns=2).resolve_scale(color="independent")
        return Success(concat_chart)
    except Exception as err:
        return Failure(f"Failed to build UMAP Altair chart: {err}")


def export_chart(chart: alt.Chart, out_svg: Path, out_png: Path) -> Result[tuple[Path, Path], str]:
    """Export Altair chart as SVG and PNG files."""
    try:
        out_svg.parent.mkdir(parents=True, exist_ok=True)
        out_png.parent.mkdir(parents=True, exist_ok=True)

        print(f"Saving SVG plot to {out_svg}...")
        chart.save(str(out_svg))

        print(f"Saving PNG plot to {out_png}...")
        chart.save(str(out_png), ppi=300)

        return Success((out_svg, out_png))
    except Exception as err:
        return Failure(f"Failed exporting plots: {err}")


def run_pipeline(config: PlotConfig) -> Result[tuple[Path, Path], str]:
    """Monadic workflow for loading AnnData, building UMAP charts, and saving plots."""
    match load_anndata(config.input_h5ad):
        case Failure(err):
            return Failure(err)
        case Success(adata):
            match build_umap_chart(adata):
                case Failure(err):
                    return Failure(err)
                case Success(chart):
                    return export_chart(chart, config.out_svg, config.out_png)
                case _ as unreachable:
                    assert_never(unreachable)
        case _ as unreachable:
            assert_never(unreachable)


def main() -> None:
    """CLI entry point for GSE97168 UMAP plotting."""
    parser = argparse.ArgumentParser(description="Plot GSE97168 UMAP embeddings with multiple metadata colorings.")
    _ = parser.add_argument("--adata", required=True, help="Input preprocessed AnnData HDF5 file.")
    _ = parser.add_argument("--out-svg", required=True, help="Output filepath for SVG plot.")
    _ = parser.add_argument("--out-png", required=True, help="Output filepath for PNG plot.")

    args = parser.parse_args()
    config = PlotConfig(
        input_h5ad=Path(args.adata),
        out_svg=Path(args.out_svg),
        out_png=Path(args.out_png),
    )

    match run_pipeline(config):
        case Success((svg_path, png_path)):
            print(f"Successfully generated GSE97168 UMAP plots:\n  - {svg_path}\n  - {png_path}")
            sys.exit(0)
        case Failure(err):
            print(f"Plotting error: {err}", file=sys.stderr)
            sys.exit(1)
        case _ as unreachable:
            assert_never(unreachable)


if __name__ == "__main__":
    main()
