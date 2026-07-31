import argparse
import sys
from pathlib import Path
from typing import Callable
import polars as pl  # type: ignore
import altair as alt  # type: ignore
from returns.result import Result, Success, Failure  # type: ignore
from returns.pipeline import flow  # type: ignore
from returns.pointfree import bind  # type: ignore

# Ensure altair can handle large datasets
alt.data_transformers.disable_max_rows()

def load_results(csv_path: Path) -> Result[pl.DataFrame, str]:
    try:
        df = pl.read_csv(csv_path)
        # Rename the index column if it exists
        if "" in df.columns:
            df = df.rename({"": "Nhood"})
        elif "Unnamed: 0" in df.columns:
            df = df.rename({"Unnamed: 0": "Nhood"})
        return Success(df)
    except Exception as e:
        return Failure(f"Failed to load {csv_path}: {e}")

def determine_significance(df: pl.DataFrame) -> pl.DataFrame:
    # Add a column indicating if FDR < 0.1
    return df.with_columns(
        (pl.col("FDR") < 0.1).alias("Significant"),
        (-pl.col("FDR").log10()).alias("log10_FDR_inv")
    )

def plot_volcano(df: pl.DataFrame, condition_name: str, out_dir: Path) -> Result[bool, str]:
    try:
        chart = alt.Chart(df).mark_point(filled=True, opacity=0.7).encode(
            x=alt.X("logFC:Q", title="Log2 Fold Change"),
            y=alt.Y("log10_FDR_inv:Q", title="-log10(FDR)"),
            color=alt.Color("Significant:N", 
                            title="FDR < 0.1", 
                            scale=alt.Scale(domain=[True, False], range=['#d62728', '#aec7e8'])),
            tooltip=["Nhood", "logFC", "FDR", "PValue"]
        ).properties(
            title=f"Differential Abundance - {condition_name}",
            width=400,
            height=400
        ).interactive()
        
        plot_path = out_dir / f"milopy_volcano_{condition_name}.svg"
        chart.save(str(plot_path))
        return Success(True)
    except Exception as e:
        return Failure(f"Volcano plotting failed: {e}")

def plot_pval_hist(df: pl.DataFrame, condition_name: str, out_dir: Path) -> Result[bool, str]:
    try:
        chart = alt.Chart(df).mark_bar(opacity=0.8, color="#4c78a8").encode(
            x=alt.X("PValue:Q", bin=alt.Bin(maxbins=50), title="P-Value"),
            y=alt.Y("count():Q", title="Count")
        ).properties(
            title=f"P-Value Distribution - {condition_name}",
            width=400,
            height=300
        )
        
        plot_path = out_dir / f"milopy_pval_hist_{condition_name}.svg"
        chart.save(str(plot_path))
        return Success(True)
    except Exception as e:
        return Failure(f"P-Value histogram plotting failed: {e}")

def process_file(csv_path: Path, out_dir: Path) -> Result[bool, str]:
    # Extract condition from milopy_results_Combined.csv
    condition_name = csv_path.stem.replace("milopy_results_", "")
    
    def run_plots(df: pl.DataFrame) -> Result[bool, str]:
        df_sig = determine_significance(df)
        res1 = plot_volcano(df_sig, condition_name, out_dir)
        if not isinstance(res1, Success): return res1
        
        res2 = plot_pval_hist(df_sig, condition_name, out_dir)
        if not isinstance(res2, Success): return res2
        
        return Success(True)

    res = flow(
        load_results(csv_path),
        bind(run_plots)
    )
    return res

def run_plotting(data_dir: Path, out_dir: Path) -> Result[bool, str]:
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        files = list(data_dir.glob("milopy_results_*.csv"))
        
        if not files:
            return Failure("No milopy result CSV files found.")
            
        for f in files:
            print(f"Plotting {f.name}...")
            res = process_file(f, out_dir)
            if not isinstance(res, Success):
                print(f"Warning: {res.failure()}")
        return Success(True)
    except Exception as e:
        return Failure(str(e))

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, help="Directory with milopy result CSVs")
    parser.add_argument("--out-dir", required=True, help="Directory to save plots")
    args = parser.parse_args()
    
    match run_plotting(Path(args.data_dir), Path(args.out_dir)):
        case Success(_):
            print("Successfully completed plotting.")
            sys.exit(0)
        case Failure(err):
            print(f"Error: {err}")
            sys.exit(1)

if __name__ == "__main__":
    main()
