"""Declarative plotting script for GSE123139 UMAP embeddings using Altair."""

import argparse
import sys
from pathlib import Path
from typing import Final, assert_never

import altair as alt  # type: ignore
import anndata as ad  # type: ignore
import polars as pl  # type: ignore
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

# Allow large single-cell datasets in Altair
alt.data_transformers.disable_max_rows()

COLOR_FIELDS: Final[tuple[str, ...]] = (
    "sample_source",
    "leiden",
    "patient_id",
    "seq_batch",
    "facs_gate",
)


class PlotConfig(BaseModel):
    """Immutable configuration for GSE123139 UMAP plotting."""

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

        # Build coordinates dataframe in Polars
        coords_df = pl.DataFrame(
            {
                "cell_id": [str(c) for c in adata.obs_names],
                "UMAP 1": umap_coords[:, 0],
                "UMAP 2": umap_coords[:, 1],
            }
        )

        # Build metadata dataframe in Polars
        obs_pd = adata.obs.reset_index()
        obs_df = pl.from_pandas(obs_pd)

        # Standardize cell_id column for join
        id_col = obs_df.columns[0]
        obs_df = obs_df.rename({id_col: "cell_id"})

        # Join UMAP coordinates with metadata
        plot_df = coords_df.join(obs_df, on="cell_id", how="left")

        # Downsample for responsive vector rendering if cell count > 30,000
        if plot_df.height > 30000:
            print(f"Downsampling {plot_df.height} cells to 30,000 for vector visualization rendering...")
            plot_df = plot_df.sample(n=30000, seed=42)

        charts: list[alt.Chart] = []

        for field in COLOR_FIELDS:
            if field not in plot_df.columns:
                print(f"Field '{field}' not found in AnnData obs metadata. Skipping.")
                continue

            # Replace empty strings with 'Unannotated' for clean color mapping
            field_df = plot_df.with_columns(
                pl.when(pl.col(field).is_null() | (pl.col(field) == ""))
                .then(pl.lit("Unannotated"))
                .otherwise(pl.col(field))
                .alias(field)
            )

            if field_df[field].n_unique() <= 1:
                continue

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
                        alt.Tooltip("sample_source:N", title="Sample Source"),
                    ],
                )
                .properties(
                    title=f"GSE123139 Melanoma UMAP: {field.replace('_', ' ').title()}",
                    width=380,
                    height=380,
                )
            )
            charts.append(chart)

        if not charts:
            return Failure("No valid metadata fields were available for plotting UMAPs.")

        # Grid arrangement: 2 columns with independent color legends
        concat_chart = alt.concat(*charts, columns=2).resolve_scale(color="independent")
        return Success(concat_chart)
    except Exception as err:
        return Failure(f"Failed to build UMAP Altair chart: {err}")


def export_chart(chart: alt.Chart, out_svg: Path, out_png: Path) -> Result[tuple[Path, Path], str]:
    """Export Altair chart as SVG and PNG images."""
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
    """CLI entry point for GSE123139 UMAP plotting."""
    parser = argparse.ArgumentParser(description="Plot GSE123139 UMAP embeddings with multiple metadata colorings.")
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
            print(f"Successfully generated GSE123139 UMAP plots:\n  - {svg_path}\n  - {png_path}")
            sys.exit(0)
        case Failure(err):
            print(f"Plotting error: {err}", file=sys.stderr)
            sys.exit(1)
        case _ as unreachable:
            assert_never(unreachable)


if __name__ == "__main__":
    main()
