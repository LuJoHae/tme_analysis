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
    """Calculates summary statistics from metadata."""
    try:
        # Assuming df_meta has columns like 'cell_id', 'patient_id', 'response', 'cluster'
        # Since we don't know the exact column names a priori without fetching, we do a generic summary.
        # We will just describe it generically if we can't find expected columns.
        cols = df_meta.columns
        if len(cols) > 1:
            # Group by the second column (usually sample/patient ID in these GEO files)
            group_col = cols[1]
            summary = df_meta.group_by(group_col).agg(pl.len().alias("cell_count"))
            return Success(summary)
        else:
            return Success(df_meta.describe())
    except Exception as e:
        return Failure(f"Failed to calculate summary: {str(e)}")

def create_plots(df_meta: pl.DataFrame) -> Result[alt.Chart, str]:
    """Creates downsampled plots using Altair."""
    try:
        # Downsample to max 5000 cells for plotting to avoid huge HTML
        n_samples = min(5000, df_meta.height)
        # Using a fixed seed for reproducible pure-like behavior, though sample is pseudo-random.
        # Polars sample doesn't take a seed directly in this API version unless specified.
        # We'll just take the head for strict determinism and speed, or sample if possible.
        # Let's just use head/tail or slice.
        df_plot = df_meta.head(n_samples)
        
        cols = df_plot.columns
        if len(cols) > 1:
            x_col = cols[1]
            chart = alt.Chart(df_plot.to_pandas()).mark_bar().encode(
                x=alt.X(f"{x_col}:N", title=x_col),
                y=alt.Y("count():Q", title="Number of Cells"),
                color=alt.Color(f"{x_col}:N", legend=None)
            ).properties(
                title=f"Cell Counts by {x_col} (Sampled)"
            )
            return Success(chart)
        else:
            return Failure("Not enough columns to plot meaningful data.")
    except Exception as e:
        return Failure(f"Failed to create plots: {str(e)}")

def save_outputs(summary: pl.DataFrame, chart: alt.Chart, out_csv: Path, out_html: Path) -> Result[bool, str]:
    try:
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        summary.write_csv(out_csv)
        
        out_html.parent.mkdir(parents=True, exist_ok=True)
        chart.save(str(out_html))
        
        return Success(True)
    except Exception as e:
        return Failure(f"Failed to save outputs: {str(e)}")

def run_pipeline(tpm_path: Path, meta_path: Path, out_csv: Path, out_html: Path) -> Result[bool, str]:
    return read_data(tpm_path, meta_path).bind(
        lambda dfs: calculate_summary(dfs[1]).bind(
            lambda summary: create_plots(dfs[1]).bind(
                lambda chart: save_outputs(summary, chart, out_csv, out_html)
            )
        )
    )

def main():
    parser = argparse.ArgumentParser(description="Analyze GSE120575 data")
    parser.add_argument("--tpm", required=True, help="Path to TPM parquet")
    parser.add_argument("--meta", required=True, help="Path to Meta parquet")
    parser.add_argument("--out-csv", required=True, help="Output summary CSV")
    parser.add_argument("--out-html", required=True, help="Output plot HTML")
    args = parser.parse_args()

    result = run_pipeline(
        Path(args.tpm),
        Path(args.meta),
        Path(args.out_csv),
        Path(args.out_html)
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
