import sys
from pathlib import Path
import polars as pl
import altair as alt
from returns.result import Result, Success, Failure

import argparse

# Enable Altair to handle somewhat larger datasets, though we will downsample
alt.data_transformers.disable_max_rows()

def read_data(tpm_path: Path, meta_path: Path) -> Result[tuple[pl.DataFrame, pl.DataFrame], str]:
    try:
        df_tpm = pl.read_parquet(tpm_path)
        df_meta = pl.read_parquet(meta_path)
        return Success((df_tpm, df_meta))
    except Exception as e:
        return Failure(f"Failed to read parquet data: {str(e)}")

def calculate_summary(df_meta: pl.DataFrame) -> Result[pl.DataFrame, str]:
    """Calculates summary statistics from metadata by finding low-cardinality group columns."""
    try:
        cols = df_meta.columns
        if len(cols) > 1:
            # Group by string columns with low cardinality (likely Patient/Response)
            group_cols = [c for c in cols if df_meta[c].n_unique() < 50 and df_meta[c].dtype in (pl.String, pl.Categorical)]
            if not group_cols:
                group_cols = [cols[1]]
                
            summary = df_meta.group_by(group_cols).agg(pl.len().alias("cell_count")).sort("cell_count", descending=True)
            return Success(summary)
        else:
            return Success(df_meta.describe())
    except Exception as e:
        return Failure(f"Failed to calculate summary: {str(e)}")

def create_plots(df_meta: pl.DataFrame) -> Result[alt.Chart, str]:
    """Creates downsampled plots using Altair."""
    try:
        n_samples = min(5000, df_meta.height)
        df_plot = df_meta.head(n_samples)
        
        cols = df_plot.columns
        if len(cols) > 1:
            candidate_cols = [c for c in cols if df_plot[c].n_unique() < 50 and df_plot[c].dtype in (pl.String, pl.Categorical)]
            x_col = candidate_cols[0] if candidate_cols else cols[1]
            color_col = candidate_cols[1] if len(candidate_cols) > 1 else x_col
            
            chart = alt.Chart(df_plot.to_pandas()).mark_bar().encode(
                x=alt.X(f"{x_col}:N", title=x_col, sort='-y'),
                y=alt.Y("count():Q", title="Number of Cells"),
                color=alt.Color(f"{color_col}:N", title=color_col)
            ).properties(
                title=f"Cell Counts by {x_col} (Sampled)"
            )
            return Success(chart)
        else:
            return Failure("Not enough columns to plot meaningful data.")
    except Exception as e:
        return Failure(f"Failed to create plots: {str(e)}")

def save_outputs(summary: pl.DataFrame, chart: alt.Chart, out_csv: Path, out_html: Path, out_svg: Path) -> Result[bool, str]:
    try:
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        summary.write_csv(out_csv)
        
        out_html.parent.mkdir(parents=True, exist_ok=True)
        chart.save(str(out_html))
        
        out_svg.parent.mkdir(parents=True, exist_ok=True)
        chart.save(str(out_svg))
        
        return Success(True)
    except Exception as e:
        return Failure(f"Failed to save outputs: {str(e)}")

def run_pipeline(tpm_path: Path, meta_path: Path, out_csv: Path, out_html: Path, out_svg: Path) -> Result[bool, str]:
    return read_data(tpm_path, meta_path).bind(
        lambda dfs: calculate_summary(dfs[1]).bind(
            lambda summary: create_plots(dfs[1]).bind(
                lambda chart: save_outputs(summary, chart, out_csv, out_html, out_svg)
            )
        )
    )

def main():
    parser = argparse.ArgumentParser(description="Analyze GSE120575 data")
    parser.add_argument("--tpm", required=True, help="Path to TPM parquet")
    parser.add_argument("--meta", required=True, help="Path to Meta parquet")
    parser.add_argument("--out-csv", required=True, help="Output summary CSV")
    parser.add_argument("--out-html", required=True, help="Output plot HTML")
    parser.add_argument("--out-svg", required=True, help="Output plot SVG")
    args = parser.parse_args()

    result = run_pipeline(
        Path(args.tpm),
        Path(args.meta),
        Path(args.out_csv),
        Path(args.out_html),
        Path(args.out_svg)
    )

    match result:
        case Success(_):
            print("Successfully processed GSE120575 data.")
            sys.exit(0)
        case Failure(err):
            print(f"Error processing data: {err}")
            sys.exit(1)

if __name__ == "__main__":
    main()
