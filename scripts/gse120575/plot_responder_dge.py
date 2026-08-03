import argparse
import sys
from pathlib import Path
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

def plot_volcano(csv_path: Path, out_dir: Path, prefix: str, title: str) -> None:
    df = pd.read_csv(csv_path)
    df["-log10(FDR)"] = -np.log10(df["pvals_adj"] + 1e-300)
    
    plt.figure(figsize=(8, 6))
    
    df["sig"] = "Not Sig"
    df.loc[(df["pvals_adj"] < 0.05) & (df["logfoldchanges"] > 0.5), "sig"] = "Up in Responder"
    df.loc[(df["pvals_adj"] < 0.05) & (df["logfoldchanges"] < -0.5), "sig"] = "Down in Responder"
    
    sns.scatterplot(
        data=df, x="logfoldchanges", y="-log10(FDR)", hue="sig", 
        palette={"Not Sig": "grey", "Up in Responder": "red", "Down in Responder": "blue"},
        s=20, alpha=0.8
    )
    
    # Annotate top 5 up and down
    top_up = df[df["sig"] == "Up in Responder"].sort_values("pvals_adj").head(10)
    top_down = df[df["sig"] == "Down in Responder"].sort_values("pvals_adj").head(10)
    
    for _, row in pd.concat([top_up, top_down]).iterrows():
        plt.text(row["logfoldchanges"] + 0.05, row["-log10(FDR)"], row["names"], fontsize=9)
        
    plt.axvline(x=0.5, color="black", linestyle="--", linewidth=0.5)
    plt.axvline(x=-0.5, color="black", linestyle="--", linewidth=0.5)
    plt.axhline(y=-np.log10(0.05), color="black", linestyle="--", linewidth=0.5)
    
    plt.title(f"Responder vs Non-Responder DGE ({title})")
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10(FDR)")
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(out_dir / f"{prefix}_responder_dge_volcano.svg", format="svg")
    plt.savefig(out_dir / f"{prefix}_responder_dge_volcano.png", format="png", dpi=300)
    plt.close()

def plot_dotplot(adata_path: Path, csv_path: Path, out_dir: Path, prefix: str) -> None:
    adata = ad.read_h5ad(adata_path)
    adata.obs["response"] = adata.obs["characteristics: response"]
    adata = adata[~adata.obs["response"].isna()].copy()
    
    df = pd.read_csv(csv_path)
    top_up = df[(df["pvals_adj"] < 0.05) & (df["logfoldchanges"] > 0)].sort_values("pvals_adj").head(15)["names"].tolist()
    top_down = df[(df["pvals_adj"] < 0.05) & (df["logfoldchanges"] < 0)].sort_values("pvals_adj").head(15)["names"].tolist()
    
    genes = top_up + top_down
    if genes:
        sc.settings.set_figure_params(dpi=150, frameon=False, figsize=(8, 4))
        # Need to recompute highly_variable or just force it so dotplot can see all values,
        # but the adata is already subsetted to 2000 genes so it's fine.
        sc.pl.dotplot(adata, var_names=genes, groupby="response", show=False, dendrogram=False, standard_scale='var')
        plt.savefig(out_dir / f"{prefix}_responder_dge_dotplot.svg", bbox_inches='tight')
        plt.savefig(out_dir / f"{prefix}_responder_dge_dotplot.png", format="png", dpi=300, bbox_inches='tight')
        plt.close()

def main() -> None:
    parser = argparse.ArgumentParser(description="Plot DEGs between Responder and Non-responder cells")
    parser.add_argument("--adata", required=True, help="Input cell subset h5ad path")
    parser.add_argument("--csv", required=True, help="Input CSV path for DGE results")
    parser.add_argument("--out-dir", required=True, help="Output directory for plots")
    parser.add_argument("--prefix", required=True, help="Prefix for output files (e.g. cd8 or cd4)")
    parser.add_argument("--title", required=True, help="Plot title suffix (e.g. 'CD8+ T cells')")
    args = parser.parse_args()
    
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("Plotting Volcano...")
    plot_volcano(Path(args.csv), out_dir, args.prefix, args.title)
    
    print("Plotting DotPlot...")
    plot_dotplot(Path(args.adata), Path(args.csv), out_dir, args.prefix)
    
    print("Done!")

if __name__ == "__main__":
    main()
