import argparse
import sys
from pathlib import Path
import pandas as pd # type: ignore
import altair as alt # type: ignore
from returns.result import Result, Success, Failure # type: ignore
from returns.pipeline import flow
from returns.pointfree import bind

def plot_single_overlap(csv_path: Path, out_dir: Path) -> Result[bool, str]:
    try:
        df = pd.read_csv(csv_path, index_col=0)
        if df.empty:
            return Failure(f"Empty dataframe for {csv_path.name}")
            
        k1 = df.index.name
        k2 = df.columns.name
        
        # Normalize by row to see proportion of k1 that falls into k2
        df_norm = df.div(df.sum(axis=1), axis=0).reset_index()
        
        long_df = df_norm.melt(id_vars=[df_norm.columns[0]], var_name=k2, value_name="Overlap Proportion")
        long_df.rename(columns={df_norm.columns[0]: k1}, inplace=True)
        long_df[k1] = long_df[k1].astype(str)
        long_df[k2] = long_df[k2].astype(str)
        
        chart = alt.Chart(long_df).mark_rect().encode(
            x=alt.X(f"{k2}:N", title=k2),
            y=alt.Y(f"{k1}:N", title=k1),
            color=alt.Color("Overlap Proportion:Q", scale=alt.Scale(scheme="viridis")),
            tooltip=[k1, k2, "Overlap Proportion"]
        ).properties(
            title=f"Overlap: {k1} vs {k2}",
            width=500,
            height=500
        )
        
        plot_path = out_dir / csv_path.name.replace(".csv", ".svg")
        chart.save(str(plot_path))
        
        return Success(True)
    except Exception as e:
        return Failure(f"Failed to plot {csv_path.name}: {e}")

def run_plotting(data_dir: Path, out_dir: Path) -> Result[bool, str]:
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        overlap_files = list(data_dir.glob("cluster_overlap_*.csv"))
        
        if not overlap_files:
            return Failure("No overlap CSVs found.")
            
        for file in overlap_files:
            print(f"Plotting {file.name}...")
            res = plot_single_overlap(file, out_dir)
            if not isinstance(res, Success):
                print(f"Warning: {res.failure()}")
                
        return Success(True)
    except Exception as e:
        return Failure(str(e))

def main() -> None:
    parser = argparse.ArgumentParser(description="Plot pairwise cluster overlaps")
    parser.add_argument("--data-dir", required=True, help="Directory containing overlap CSVs")
    parser.add_argument("--out-dir", required=True, help="Output directory for plots")
    args = parser.parse_args()
    
    match run_plotting(Path(args.data_dir), Path(args.out_dir)):
        case Success(_):
            print("Successfully completed overlap plotting.")
            sys.exit(0)
        case Failure(err):
            print(f"Error: {err}")
            sys.exit(1)

if __name__ == "__main__":
    main()
