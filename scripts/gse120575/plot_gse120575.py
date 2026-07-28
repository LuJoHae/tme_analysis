import sys
import argparse
from pathlib import Path
import polars as pl
import anndata as ad
import altair as alt
from returns.result import Result, Success, Failure

# Allow large datasets in Altair
alt.data_transformers.disable_max_rows()

def load_anndata(h5ad_path: Path) -> Result[ad.AnnData, str]:
    try:
        adata = ad.read_h5ad(h5ad_path)
        return Success(adata)
    except Exception as e:
        return Failure(f"Failed to load AnnData: {e}")

def create_umap_plots(adata: ad.AnnData) -> Result[alt.Chart, str]:
    try:
        umap_coords = adata.obsm['X_umap']
        
        plot_df = pl.DataFrame({
            "cell_id": adata.obs_names,
            "UMAP1": umap_coords[:, 0],
            "UMAP2": umap_coords[:, 1],
        })
        
        obs_df = pl.from_pandas(adata.obs.reset_index())
        if "index" in obs_df.columns:
            obs_df = obs_df.rename({"index": "cell_id"})
            
        plot_df = plot_df.join(obs_df, on="cell_id")
        
        exclude_cols = {"cell_id", "plate-row", "plate-col", "n_genes", "n_counts"}
        
        charts = []
        for col in adata.obs.columns:
            if col in exclude_cols or plot_df[col].null_count() == plot_df.height:
                continue
                
            # Skip columns that have only a single unique value across all cells (e.g. organism=Homo sapiens)
            if plot_df[col].n_unique() <= 1:
                continue
                
            chart = alt.Chart(plot_df).mark_circle(size=5, opacity=0.8).encode(
                x=alt.X("UMAP1:Q", title="UMAP 1"),
                y=alt.Y("UMAP2:Q", title="UMAP 2"),
                color=alt.Color(field=col, type="nominal", title=col),
                tooltip=["cell_id", alt.Tooltip(field=col, type="nominal")]
            ).properties(
                title=f"UMAP colored by {col}",
                width=350,
                height=350
            )
            charts.append(chart)
            
        if not charts:
            return Failure("No valid metadata columns found to plot.")
            
        # Layout in a grid (2 columns wide) using alt.concat.
        # resolve_scale ensures each subplot gets its own legend instead of sharing one large legend area.
        final_chart = alt.concat(*charts, columns=2).resolve_scale(color='independent')
        
        return Success(final_chart)
    except Exception as e:
        return Failure(f"Failed to create plots: {e}")

def save_plots(chart: alt.Chart, out_svg: Path, out_png: Path) -> Result[bool, str]:
    try:
        out_svg.parent.mkdir(parents=True, exist_ok=True)
        chart.save(str(out_svg))
        
        out_png.parent.mkdir(parents=True, exist_ok=True)
        chart.save(str(out_png), ppi=300)
        return Success(True)
    except Exception as e:
        return Failure(f"Failed to save plots: {e}")

def run_pipeline(h5ad_path: Path, out_svg: Path, out_png: Path) -> Result[bool, str]:
    return load_anndata(h5ad_path).bind(
        lambda adata: create_umap_plots(adata).bind(
            lambda chart: save_plots(chart, out_svg, out_png)
        )
    )

def main():
    parser = argparse.ArgumentParser(description="Plot GSE120575 UMAPs")
    parser.add_argument("--adata", required=True, help="Path to processed AnnData h5ad")
    parser.add_argument("--out-svg", required=True, help="Output plot SVG")
    parser.add_argument("--out-png", required=True, help="Output plot PNG")
    args = parser.parse_args()

    match run_pipeline(Path(args.adata), Path(args.out_svg), Path(args.out_png)):
        case Success(_):
            print("Successfully completed plotting and saved UMAP plots.")
            sys.exit(0)
        case Failure(err):
            print(f"Error plotting data: {err}")
            sys.exit(1)

if __name__ == "__main__":
    main()
