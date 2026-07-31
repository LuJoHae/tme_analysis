import argparse
import sys
from pathlib import Path
import polars as pl  # type: ignore
import pandas as pd  # type: ignore
import numpy as np  # type: ignore
import scipy.sparse as sp  # type: ignore
import anndata as ad  # type: ignore
import altair as alt  # type: ignore
from returns.result import Result, Success, Failure  # type: ignore

alt.data_transformers.disable_max_rows()

def load_results(csv_path: Path) -> Result[pl.DataFrame, str]:
    try:
        df = pl.read_csv(csv_path)
        if "" in df.columns:
            df = df.rename({"": "Nhood"})
        elif "Unnamed: 0" in df.columns:
            df = df.rename({"Unnamed: 0": "Nhood"})
        return Success(df)
    except Exception as e:
        return Failure(f"Failed to load {csv_path}: {e}")

def get_significant_nhoods(df: pl.DataFrame, fdr_thresh: float = 0.1) -> pl.DataFrame:
    # Returns DataFrame with Nhood and logFC
    return df.filter(pl.col("FDR") < fdr_thresh).select(["Nhood", "logFC"])

def plot_composition(long_df: pd.DataFrame, condition: str, col: str, out_dir: Path) -> None:
    # Top chart: logFC (effect size and direction)
    effect_chart = alt.Chart(long_df).mark_bar().encode(
        x=alt.X("Nhood:O", title=None, axis=alt.Axis(labels=False, ticks=False), sort=alt.EncodingSortField(field="logFC", order="descending")),
        y=alt.Y("logFC:Q", title="log2(FC)"),
        color=alt.condition(
            alt.datum.logFC > 0,
            alt.value("#d62728"),  # Red for positive logFC (Enriched)
            alt.value("#1f77b4")   # Blue for negative logFC (Depleted)
        ),
        tooltip=["Nhood", "logFC"]
    ).properties(
        width=800,
        height=100
    )

    # Bottom chart: Composition
    comp_chart = alt.Chart(long_df).mark_bar().encode(
        x=alt.X("Nhood:O", title="Significant Neighborhood ID", sort=alt.EncodingSortField(field="logFC", order="descending")),
        y=alt.Y("Proportion:Q", title="Fraction of Cells"),
        color=alt.Color("Cluster:N", scale=alt.Scale(scheme="category20")),
        tooltip=["Nhood", "Cluster", "Proportion", "logFC"]
    ).properties(
        width=800,
        height=300
    )

    chart = alt.vconcat(effect_chart, comp_chart).resolve_scale(x='shared').properties(
        title=f"{condition} - {col} composition of significant neighborhoods"
    )
    
    chart.save(str(out_dir / f"milopy_nhood_composition_{condition}_{col}.svg"))

def process_condition(condition: str, data_dir: Path, out_dir: Path, adata: ad.AnnData) -> Result[bool, str]:
    print(f"Processing condition {condition}...")
    
    csv_path = data_dir / f"milopy_results_{condition}.csv"
    if not csv_path.exists():
        return Failure(f"Results CSV not found: {csv_path}")
        
    df_res = load_results(csv_path)
    if not isinstance(df_res, Success):
        return df_res
        
    sig_df = get_significant_nhoods(df_res.unwrap(), fdr_thresh=0.1)
    if sig_df.height == 0:
        print(f"No significant neighborhoods found for {condition} at FDR < 0.1. Skipping.")
        return Success(True)
        
    sig_nhoods = sig_df.get_column("Nhood").to_list()
    sig_logfc = sig_df.get_column("logFC").to_list()
    logfc_dict = dict(zip(sig_nhoods, sig_logfc))
    
    if condition != "Combined":
        obs = adata.obs[adata.obs["treatment_status"] == condition].copy()
    else:
        obs = adata.obs.copy()
        
    npz_path = data_dir / f"milopy_nhoods_{condition}.npz"
    if not npz_path.exists():
        return Failure(f"Sparse matrix not found: {npz_path}")
        
    nhoods = sp.load_npz(npz_path)
    if nhoods.shape[0] != len(obs):
        return Failure(f"Shape mismatch: {nhoods.shape[0]} cells in npz vs {len(obs)} in metadata for {condition}")
        
    clustering_cols = ["celltypist_leiden_0.5", "celltypist_leiden_1.0", "celltypist_leiden_1.5", "celltypist_leiden_2.0"]
    
    for col in clustering_cols:
        if col not in obs.columns:
            print(f"Warning: {col} not found in adata.obs. Skipping.")
            continue
            
        clusters = obs[col].astype("category")
        cluster_names = clusters.cat.categories
        
        one_hot = sp.csr_matrix(pd.get_dummies(clusters).values)
        comp = (one_hot.T @ nhoods).T
        comp_sig = comp[sig_nhoods, :].toarray()
        
        row_sums = comp_sig.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1 
        props = comp_sig / row_sums
        
        plot_df = pd.DataFrame(props, columns=cluster_names, index=sig_nhoods)
        plot_df.index.name = "Nhood"
        long_df = plot_df.reset_index().melt(id_vars="Nhood", var_name="Cluster", value_name="Proportion")
        
        long_df = long_df[long_df["Proportion"] > 0].copy()
        long_df["logFC"] = long_df["Nhood"].map(logfc_dict)
        
        plot_composition(long_df, condition, col, out_dir)
        
    return Success(True)

def run_all(adata_path: Path, data_dir: Path, out_dir: Path) -> Result[bool, str]:
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        print("Loading AnnData...")
        adata = ad.read_h5ad(adata_path, backed="r")
        
        for cond in ["Pre", "Post", "Combined"]:
            res = process_condition(cond, data_dir, out_dir, adata)
            if not isinstance(res, Success):
                print(f"Error processing {cond}: {res.failure()}")
                
        return Success(True)
    except Exception as e:
        return Failure(str(e))

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--adata", required=True, help="Path to processed h5ad file")
    parser.add_argument("--data-dir", required=True, help="Directory with milopy results and npz files")
    parser.add_argument("--out-dir", required=True, help="Directory to save plots")
    args = parser.parse_args()
    
    match run_all(Path(args.adata), Path(args.data_dir), Path(args.out_dir)):
        case Success(_):
            print("Successfully completed composition plotting.")
            sys.exit(0)
        case Failure(err):
            print(f"Error: {err}")
            sys.exit(1)

if __name__ == "__main__":
    main()
