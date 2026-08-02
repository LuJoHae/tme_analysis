import argparse
import sys
from pathlib import Path
import polars as pl  # type: ignore
import pandas as pd  # type: ignore
import numpy as np  # type: ignore
import scipy.sparse as sp  # type: ignore
import anndata as ad  # type: ignore
import altair as alt  # type: ignore
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scanpy as sc
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
    # Returns DataFrame with Nhood, logFC, and NhoodGroup
    if "NhoodGroup" not in df.columns:
        df = df.with_columns(pl.lit(None).alias("NhoodGroup"))
    return df.filter(pl.col("FDR") < fdr_thresh).filter(pl.col("NhoodGroup").is_not_null()).select(["Nhood", "logFC", "NhoodGroup"])

def plot_composition(long_df: pd.DataFrame, condition: str, col: str, out_dir: Path, x_col: str, prefix: str, title_suffix: str) -> None:
    long_df = long_df.copy()
    if x_col == "NhoodGroup":
        long_df[x_col] = long_df[x_col].astype(int).astype(str)
    else:
        long_df[x_col] = long_df[x_col].astype(str)
        
    effect_df = long_df[[x_col, "logFC"]].drop_duplicates()
    if "stats_text" in long_df.columns:
        effect_df = effect_df.merge(long_df[[x_col, "stats_text"]].drop_duplicates(), on=x_col)
        
    # Extract the exact desired order from pandas to guarantee bulletproof sorting
    sort_order = effect_df.sort_values("logFC", ascending=False)[x_col].tolist()
        
    effect_chart = alt.Chart(effect_df).mark_bar().encode(
        x=alt.X(f"{x_col}:N", title=None, axis=alt.Axis(labels=False, ticks=False), sort=sort_order),
        y=alt.Y("logFC:Q", title="log2(FC)"),
        color=alt.condition(
            alt.datum.logFC > 0,
            alt.value("#d62728"),
            alt.value("#1f77b4")
        ),
        tooltip=[x_col, "logFC"]
    ).properties(
        width=800,
        height=100
    )

    comp_chart = alt.Chart(long_df).mark_bar().encode(
        x=alt.X(f"{x_col}:N", title=title_suffix, sort=sort_order),
        y=alt.Y("Proportion:Q", title="Fraction of Cells"),
        color=alt.Color("Cluster:N", scale=alt.Scale(scheme="category20")),
        tooltip=[x_col, "Cluster", "Proportion", "logFC"]
    ).properties(
        width=800,
        height=300
    )

    if "stats_text" in effect_df.columns:
        text_chart = alt.Chart(effect_df).mark_text(
            align='left',
            baseline='middle',
            fontSize=10,
            angle=270,
            lineBreak='|'
        ).encode(
            x=alt.X(f"{x_col}:N", title=None, axis=alt.Axis(labels=False, ticks=False), sort=sort_order),
            text="stats_text:N",
            tooltip=[x_col, "stats_text"]
        ).properties(
            width=800,
            height=60
        )
        
        chart = alt.vconcat(text_chart, effect_chart, comp_chart).resolve_scale(x='shared').properties(
            title=f"{condition} - {col} composition of {title_suffix.lower()}"
        )
    else:
        chart = alt.vconcat(effect_chart, comp_chart).resolve_scale(x='shared').properties(
            title=f"{condition} - {col} composition of {title_suffix.lower()}"
        )
    
    chart.save(str(out_dir / f"milopy_nhood_composition_{condition}_{col}_{prefix}.svg"))

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
    nhood_groups = sig_df.get_column("NhoodGroup").to_list()
    group_map = dict(zip(sig_nhoods, nhood_groups))
    
    # Calculate mean logFC per module
    group_logfc = sig_df.group_by("NhoodGroup").agg(pl.col("logFC").mean()).to_dict(as_series=False)
    mean_logfc_dict = dict(zip(group_logfc["NhoodGroup"], group_logfc["logFC"]))
    
    if condition != "Combined":
        adata_sub = adata[adata.obs["treatment_status"] == condition].to_memory()
    else:
        adata_sub = adata.to_memory()
        
    obs = adata_sub.obs
        
    npz_path = data_dir / f"milopy_nhoods_{condition}.npz"
    if not npz_path.exists():
        return Failure(f"Sparse matrix not found: {npz_path}")
        
    nhoods = sp.load_npz(npz_path)
    if nhoods.shape[0] != len(obs):
        return Failure(f"Shape mismatch: {nhoods.shape[0]} cells in npz vs {len(obs)} in metadata for {condition}")
        
    # --- Milo UMAP Plotting ---
    logfc_array = np.array(sig_df.get_column("logFC").to_list())
    sig_nhoods_idx = [int(i) for i in sig_nhoods]
    nhoods_sig = nhoods[:, sig_nhoods_idx]
    
    pos_mask = logfc_array > 0
    neg_mask = logfc_array < 0
    
    n_pos = nhoods_sig[:, pos_mask].sum(axis=1).A1
    n_neg = nhoods_sig[:, neg_mask].sum(axis=1).A1
    
    total = n_pos + n_neg
    f_pos = np.full(total.shape, np.nan)
    mask = total > 0
    f_pos[mask] = n_pos[mask] / total[mask]
    
    adata_sub.obs["milo_gradient"] = f_pos
    cmap = mcolors.LinearSegmentedColormap.from_list("milo_cmap", ["blue", "purple", "red"])
    
    sc.settings.set_figure_params(dpi=150, frameon=False, figsize=(6, 6))
    sc.pl.umap(adata_sub, color="milo_gradient", cmap=cmap, na_color="lightgray", title=f"Milo Neighborhoods - {condition}", show=False, s=15, alpha=0.8)
    plt.savefig(out_dir / f"milopy_umap_gradient_{condition}.png", dpi=300, bbox_inches="tight")
    plt.close()
    # ---------------------------
        
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
        
        # 1. Plot individual neighborhoods
        long_df_nhoods = plot_df.reset_index().melt(id_vars="Nhood", var_name="Cluster", value_name="Proportion")
        long_df_nhoods = long_df_nhoods[long_df_nhoods["Proportion"] > 0].copy()
        
        # We need a logfc mapping for nhoods
        nhood_logfc_dict = dict(zip(sig_nhoods, sig_df.get_column("logFC").to_list()))
        long_df_nhoods["logFC"] = long_df_nhoods["Nhood"].map(nhood_logfc_dict)
        
        plot_composition(long_df_nhoods, condition, col, out_dir, x_col="Nhood", prefix="nhoods", title_suffix="Significant Neighborhoods")
        
        # 2. Plot aggregated modules
        plot_df["NhoodGroup"] = plot_df.index.map(group_map)
        
        group_df = plot_df.groupby("NhoodGroup").mean().reset_index()
        long_df_modules = group_df.melt(id_vars="NhoodGroup", var_name="Cluster", value_name="Proportion")
        
        long_df_modules = long_df_modules[long_df_modules["Proportion"] > 0].copy()
        long_df_modules["logFC"] = long_df_modules["NhoodGroup"].map(mean_logfc_dict)
        
        # Calculate stats text
        group_stats = []
        for g in group_df["NhoodGroup"]:
            g_nhoods = [n for n, grp in group_map.items() if grp == g]
            n_nhoods = len(g_nhoods)
            g_cells = np.array(nhoods[:, g_nhoods].sum(axis=1) > 0).sum()
            group_stats.append({"NhoodGroup": g, "stats_text": f"{g_cells} cells|{n_nhoods} nhoods"})
            
        stats_df = pd.DataFrame(group_stats)
        long_df_modules = long_df_modules.merge(stats_df, on="NhoodGroup")
        
        plot_composition(long_df_modules, condition, col, out_dir, x_col="NhoodGroup", prefix="modules", title_suffix="Milo Modules")
        
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
