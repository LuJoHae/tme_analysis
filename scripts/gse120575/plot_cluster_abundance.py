import argparse
import sys
from pathlib import Path
import pandas as pd  # type: ignore
import altair as alt  # type: ignore
from returns.result import Result, Success, Failure  # type: ignore
from returns.pipeline import flow
from returns.pointfree import bind
from returns.decorators import do  # type: ignore

# Ensure altair can handle large datasets if necessary
alt.data_transformers.disable_max_rows()

def load_proportions(csv_path: Path) -> Result[pd.DataFrame, str]:
    try:
        df = pd.read_csv(csv_path)
        return Success(df)
    except Exception as e:
        return Failure(f"Failed to load {csv_path}: {e}")

def plot_proportions(df: pd.DataFrame, cluster_key: str, condition_name: str, out_dir: Path) -> Result[bool, str]:
    try:
        # We need to unpivot from wide format to long format
        id_vars = ["melanoma-sample", "response"]
        value_vars = [c for c in df.columns if c not in id_vars]
        
        long_df = df.melt(id_vars=id_vars, value_vars=value_vars, var_name="Cluster", value_name="Proportion")
        
        # Filter out samples that have 0 sum across all clusters (dropped samples)
        sample_sums = long_df.groupby("melanoma-sample")["Proportion"].sum().reset_index()
        valid_samples = sample_sums[sample_sums["Proportion"] > 0]["melanoma-sample"].tolist()
        
        long_df = long_df[long_df["melanoma-sample"].isin(valid_samples)]
        
        if long_df.empty:
            return Failure(f"No valid data to plot for {condition_name} {cluster_key}")
            
        chart = alt.Chart(long_df).mark_bar().encode(
            x=alt.X("melanoma-sample:N", 
                    sort=alt.EncodingSortField(field="response"),
                    title="Sample", 
                    axis=alt.Axis(labelAngle=-90)),
            y=alt.Y("Proportion:Q", title="Proportion"),
            color=alt.Color("Cluster:N", scale=alt.Scale(scheme="category20")),
            column=alt.Column("response:N", title="Response", header=alt.Header(labelOrient="bottom"))
        ).properties(
            title=f"Cluster Proportions - {condition_name} ({cluster_key})",
            width=300,
            height=400
        ).resolve_scale(x='independent')
        
        plot_path = out_dir / f"sccoda_abundance_{condition_name}_{cluster_key}.png"
        
        # Since vl-convert-python is used to save altair charts
        chart.save(str(plot_path), dpi=300)
        return Success(True)
    except Exception as e:
        return Failure(f"Plotting failed: {e}")

@do(Result[bool, str])
def process_file(csv_path: Path, out_dir: Path) -> bool:
    name_parts = csv_path.stem.replace("sccoda_proportions_", "").split("_")
    condition_name = name_parts[0]
    cluster_key = "_".join(name_parts[1:])
    
    df = yield load_proportions(csv_path)
    yield plot_proportions(df, cluster_key, condition_name, out_dir)
    return True

def run_plotting(data_dir: Path, out_dir: Path) -> Result[bool, str]:
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        files = list(data_dir.glob("sccoda_proportions_*.csv"))
        
        if not files:
            return Failure("No proportion CSV files found.")
            
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
    parser.add_argument("--data-dir", required=True, help="Directory with proportion CSVs")
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
