import argparse
import sys
from pathlib import Path
import pandas as pd # type: ignore
import altair as alt # type: ignore
import anndata as ad # type: ignore
from returns.result import Result, Success, Failure # type: ignore

def plot_single_jaccard(k1: str, k2: str, df: pd.DataFrame, out_dir: Path) -> Result[bool, str]:
    try:
        if df.empty:
            return Failure(f"Empty dataframe for {k1} vs {k2}")
            
        df = df.copy()
        df.index.name = k1
        df.columns.name = k2
            
        long_df = df.reset_index().melt(id_vars=[k1], var_name=k2, value_name="Jaccard Similarity")
        
        long_df[k1] = long_df[k1].astype(str)
        long_df[k2] = long_df[k2].astype(str)
        
        chart = alt.Chart(long_df).mark_rect().encode(
            x=alt.X(f"{k2}:N", title=k2),
            y=alt.Y(f"{k1}:N", title=k1),
            color=alt.Color("Jaccard Similarity:Q", scale=alt.Scale(scheme="viridis")),
            tooltip=[k1, k2, "Jaccard Similarity"]
        ).properties(
            title=f"Jaccard Similarity: {k1} vs {k2}",
            width=500,
            height=500
        )
        
        plot_path = out_dir / f"jaccard_overlap_{k1}_vs_{k2}.svg"
        chart.save(str(plot_path))
        
        return Success(True)
    except Exception as e:
        return Failure(f"Failed to plot {k1} vs {k2}: {e}")

def run_plotting(adata_path: Path, out_dir: Path) -> Result[bool, str]:
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"Loading AnnData from {adata_path}...")
        adata = ad.read_h5ad(adata_path)
        
        jaccard_keys = [k for k in adata.uns.keys() if k.startswith("jaccard_")]
        if not jaccard_keys:
            return Failure("No Jaccard matrices found in adata.uns.")
            
        for key in jaccard_keys:
            print(f"Plotting {key}...")
            # Key format: jaccard_k1_vs_k2
            parts = key.replace("jaccard_", "").split("_vs_")
            k1 = parts[0]
            k2 = parts[1]
            
            df = adata.uns[key]
            res = plot_single_jaccard(k1, k2, df, out_dir)
            if not isinstance(res, Success):
                print(f"Warning: {res.failure()}")
                
        return Success(True)
    except Exception as e:
        return Failure(str(e))

def main() -> None:
    parser = argparse.ArgumentParser(description="Plot pairwise cluster Jaccard overlaps")
    parser.add_argument("--adata", required=True, help="Path to processed h5ad file")
    parser.add_argument("--out-dir", required=True, help="Output directory for plots")
    args = parser.parse_args()
    
    match run_plotting(Path(args.adata), Path(args.out_dir)):
        case Success(_):
            print("Successfully completed Jaccard overlap plotting.")
            sys.exit(0)
        case Failure(err):
            print(f"Error: {err}")
            sys.exit(1)

if __name__ == "__main__":
    main()
