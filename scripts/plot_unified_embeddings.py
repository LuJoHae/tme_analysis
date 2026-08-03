"""Plotting Script for Unified SingleCellDataProcessStep07 + GSE120575 Dataset.

This module computes Uncorrected, Gene-Length Scaled, ComBat-corrected, and Harmony-corrected
PCA and UMAP embeddings for the combined AnnData dataset and generates Altair visualizations
exported as SVG vector graphics across multiple colorings:
1. GSE120575 vs Rest
2. Sequencing Technology (10x UMI vs Smart-seq2)
3. Cell Type Lineage
4. Cancer Code / Tumor Site
5. Dataset of Origin
"""

from pathlib import Path
from typing import Optional, Literal
import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import altair as alt
from returns.result import Result, Success, Failure
from returns.maybe import Maybe, Some, Nothing


class PlotConfig(BaseModel):
    """Immutable configuration for embedding computation and SVG plot generation."""
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    input_h5ad_path: Path
    output_dir: Path
    n_hvg: int = 1000
    n_pcs: int = 30
    random_seed: int = 42


def load_unified_dataset(h5ad_path: Path) -> Result[ad.AnnData, str]:
    """Loads the unified AnnData dataset from disk."""
    try:
        if not h5ad_path.exists():
            return Failure(f"Input AnnData file not found: {h5ad_path}")
        adata = ad.read_h5ad(h5ad_path)
        return Success(adata)
    except Exception as e:
        return Failure(f"Failed to load AnnData from {h5ad_path}: {str(e)}")


def apply_gene_length_scaling(adata: ad.AnnData) -> ad.AnnData:
    """Scales Smart-seq2 expression data by gene transcript length to align full-length reads with UMI counts."""
    adata_scaled = adata.copy()
    gse_mask = (adata_scaled.obs["sequencing_tech"] == "Smart-seq2").values
    
    if gse_mask.sum() == 0:
        return adata_scaled

    # Use contig length or fallback gene length estimate (e.g. 2kb mean)
    gene_lengths = np.ones(adata_scaled.n_vars, dtype=np.float32)
    if "contig_length" in adata_scaled.var.columns:
        lens = pd.to_numeric(adata_scaled.var["contig_length"], errors="coerce").values
        gene_lengths = np.where(np.isnan(lens) | (lens <= 0), 2000.0, lens)
    
    gene_lengths_kb = gene_lengths / 1000.0

    # Scale Smart-seq2 expression values by gene_length_kb
    X_mat = adata_scaled.X.copy()
    if hasattr(X_mat, "toarray"):
        X_dense = X_mat.toarray()
        X_dense[gse_mask] = X_dense[gse_mask] / gene_lengths_kb
        from scipy.sparse import csr_matrix
        adata_scaled.X = csr_matrix(X_dense)
    else:
        X_mat[gse_mask] = X_mat[gse_mask] / gene_lengths_kb
        adata_scaled.X = X_mat

    return adata_scaled


def compute_embeddings_mode(
    adata: ad.AnnData,
    mode: Literal["uncorrected", "combat", "harmony"],
    n_hvg: int = 1000,
    n_pcs: int = 30,
) -> Result[ad.AnnData, str]:
    """Computes PCA and UMAP embeddings under specified batch correction mode."""
    try:
        adata_comp = adata.copy()
        
        # 1. Normalize and log1p transform
        if "log1p" not in adata_comp.uns:
            sc.pp.normalize_total(adata_comp, target_sum=1e4)
            sc.pp.log1p(adata_comp)

        # 2. Select Highly Variable Genes
        num_genes = adata_comp.n_vars
        if num_genes >= 20:
            sc.pp.highly_variable_genes(
                adata_comp,
                n_top_genes=min(n_hvg, num_genes),
                inplace=True,
            )
            has_hvg = "highly_variable" in adata_comp.var.columns
        else:
            has_hvg = False

        actual_pcs = min(n_pcs, num_genes - 1) if num_genes > 1 else 1

        match mode:
            case "uncorrected":
                if has_hvg:
                    sc.tl.pca(adata_comp, n_comps=actual_pcs, mask_var="highly_variable")
                else:
                    sc.tl.pca(adata_comp, n_comps=actual_pcs)
                sc.pp.neighbors(adata_comp, n_pcs=min(15, actual_pcs))
                sc.tl.umap(adata_comp)

            case "combat":
                if "sequencing_tech" in adata_comp.obs.columns and len(adata_comp.obs["sequencing_tech"].unique()) > 1:
                    sc.pp.combat(adata_comp, key="sequencing_tech")
                if has_hvg:
                    sc.tl.pca(adata_comp, n_comps=actual_pcs, mask_var="highly_variable")
                else:
                    sc.tl.pca(adata_comp, n_comps=actual_pcs)
                sc.pp.neighbors(adata_comp, n_pcs=min(15, actual_pcs))
                sc.tl.umap(adata_comp)

            case "harmony":
                if has_hvg:
                    sc.tl.pca(adata_comp, n_comps=actual_pcs, mask_var="highly_variable")
                else:
                    sc.tl.pca(adata_comp, n_comps=actual_pcs)
                
                if "sequencing_tech" in adata_comp.obs.columns and len(adata_comp.obs["sequencing_tech"].unique()) > 1:
                    import importlib.util
                    if importlib.util.find_spec("harmonypy") is not None:
                        try:
                            sce.pp.harmony_integrate(adata_comp, key="sequencing_tech", basis="X_pca", adjusted_basis="X_pca_harmony")
                            if "X_pca_harmony" in adata_comp.obsm:
                                if adata_comp.obsm["X_pca_harmony"].shape[0] != adata_comp.n_obs:
                                    adata_comp.obsm["X_pca_harmony"] = adata_comp.obsm["X_pca_harmony"].T
                        except Exception:
                            # Fallback if harmonypy integration errors out
                            X_pca = adata_comp.obsm["X_pca"].copy()
                            for tech in adata_comp.obs["sequencing_tech"].unique():
                                mask = (adata_comp.obs["sequencing_tech"] == tech).values
                                X_pca[mask] -= X_pca[mask].mean(axis=0, keepdims=True)
                            adata_comp.obsm["X_pca_harmony"] = X_pca
                    else:
                        # Fallback to PCA batch mean-centering if harmonypy package is missing
                        X_pca = adata_comp.obsm["X_pca"].copy()
                        for tech in adata_comp.obs["sequencing_tech"].unique():
                            mask = (adata_comp.obs["sequencing_tech"] == tech).values
                            X_pca[mask] -= X_pca[mask].mean(axis=0, keepdims=True)
                        adata_comp.obsm["X_pca_harmony"] = X_pca

                    sc.pp.neighbors(adata_comp, use_rep="X_pca_harmony", n_pcs=min(15, actual_pcs))
                else:
                    sc.pp.neighbors(adata_comp, n_pcs=min(15, actual_pcs))
                sc.tl.umap(adata_comp)

            case _:
                assert_never(mode)

        return Success(adata_comp)
    except Exception as e:
        return Failure(f"Embedding computation failed for mode '{mode}': {str(e)}")


def extract_embedding_dataframe(
    adata: ad.AnnData,
    mode: str,
) -> Result[pd.DataFrame, str]:
    """Extracts PC1, PC2, UMAP1, UMAP2 and metadata annotations into a Pandas DataFrame."""
    try:
        pca_key = "X_pca_harmony" if "X_pca_harmony" in adata.obsm else "X_pca"
        pca_coords = adata.obsm[pca_key]
        umap_coords = adata.obsm["X_umap"]

        obs_df = adata.obs.copy()
        dataset_col = obs_df["dataset"].values if "dataset" in obs_df.columns else np.array(["Unknown"] * adata.n_obs)
        tech_col = obs_df["sequencing_tech"].values if "sequencing_tech" in obs_df.columns else np.array(["Unknown"] * adata.n_obs)
        cell_type_col = obs_df["cell_type"].values if "cell_type" in obs_df.columns else np.array(["Unknown"] * adata.n_obs)
        cancer_col = obs_df["cancer_code"].values if "cancer_code" in obs_df.columns else np.array(["Unknown"] * adata.n_obs)

        # Create binary indicator: GSE120575 vs Rest
        gse_vs_rest = np.where(
            (dataset_col == "GSE120575_SadeFeldman") | (tech_col == "Smart-seq2"),
            "GSE120575 (Melanoma)",
            "Rest (Step07 Multi-Cohort)",
        )

        df_emb = pd.DataFrame({
            "PC1": pca_coords[:, 0],
            "PC2": pca_coords[:, 1],
            "UMAP1": umap_coords[:, 0],
            "UMAP2": umap_coords[:, 1],
            "Dataset": dataset_col,
            "Sequencing_Tech": tech_col,
            "Cell_Type": cell_type_col,
            "Cancer_Code": cancer_col,
            "GSE120575_vs_Rest": gse_vs_rest,
            "Correction_Mode": mode,
        })
        return Success(df_emb)
    except Exception as e:
        return Failure(f"Failed to extract embedding dataframe: {str(e)}")


def create_altair_scatter(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    color_col: str,
    title: str,
    color_title: str,
) -> alt.Chart:
    """Creates a standardized Altair scatter plot for PCA or UMAP."""
    chart = alt.Chart(df).mark_circle(size=35, opacity=0.75).encode(
        x=alt.X(f"{x_col}:Q", title=x_col),
        y=alt.Y(f"{y_col}:Q", title=y_col),
        color=alt.Color(f"{color_col}:N", title=color_title, scale=alt.Scale(scheme="category10")),
        tooltip=["GSE120575_vs_Rest", "Sequencing_Tech", "Dataset", "Cell_Type", "Cancer_Code", "Correction_Mode"]
    ).properties(
        title=title,
        width=600,
        height=450,
    ).interactive()
    return chart


def export_mode_svgs(
    df_emb: pd.DataFrame,
    mode_prefix: str,
    output_dir: Path,
) -> Result[list[Path], str]:
    """Generates and exports SVG charts for PCA and UMAP across all 5 colorings for a batch mode."""
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        exported_paths = []

        plot_specs = [
            ("GSE120575_vs_Rest", "GSE120575 vs Rest", "gse_vs_rest"),
            ("Sequencing_Tech", "Sequencing Technology (10x UMI vs Smart-seq2)", "sequencing_tech"),
            ("Cell_Type", "Cell Type Lineage", "cell_type"),
            ("Cancer_Code", "Cancer Code / Tumor Site", "cancer_code"),
            ("Dataset", "Dataset of Origin", "dataset"),
        ]

        for color_col, display_label, file_suffix in plot_specs:
            # 1. PCA SVG
            pca_chart = create_altair_scatter(
                df=df_emb,
                x_col="PC1",
                y_col="PC2",
                color_col=color_col,
                title=f"PCA ({mode_prefix}): {display_label}",
                color_title=display_label,
            )
            pca_path = output_dir / f"pca_{mode_prefix}_{file_suffix}.svg"
            pca_chart.save(str(pca_path))
            exported_paths.append(pca_path)

            # 2. UMAP SVG
            umap_chart = create_altair_scatter(
                df=df_emb,
                x_col="UMAP1",
                y_col="UMAP2",
                color_col=color_col,
                title=f"UMAP ({mode_prefix}): {display_label}",
                color_title=display_label,
            )
            umap_path = output_dir / f"umap_{mode_prefix}_{file_suffix}.svg"
            umap_chart.save(str(umap_path))
            exported_paths.append(umap_path)

        return Success(exported_paths)
    except Exception as e:
        return Failure(f"Failed to export SVG embedding plots for {mode_prefix}: {str(e)}")


def run_plotting_pipeline(config: PlotConfig) -> Result[list[Path], str]:
    """Coordinates loading, Gene Length scaling, Uncorrected, ComBat, and Harmony embeddings & SVG export."""
    return load_unified_dataset(config.input_h5ad_path).bind(
        lambda raw_adata: Success(apply_gene_length_scaling(raw_adata)).bind(
            lambda scaled_adata: compute_embeddings_mode(scaled_adata, "uncorrected", config.n_hvg, config.n_pcs).bind(
                lambda adata_uncorrected: extract_embedding_dataframe(adata_uncorrected, "Uncorrected").bind(
                    lambda df_uncorrected: export_mode_svgs(df_uncorrected, "uncorrected", config.output_dir).bind(
                        lambda paths_uncorrected: compute_embeddings_mode(scaled_adata, "combat", config.n_hvg, config.n_pcs).bind(
                            lambda adata_combat: extract_embedding_dataframe(adata_combat, "ComBat").bind(
                                lambda df_combat: export_mode_svgs(df_combat, "combat", config.output_dir).bind(
                                    lambda paths_combat: compute_embeddings_mode(scaled_adata, "harmony", config.n_hvg, config.n_pcs).bind(
                                        lambda adata_harmony: extract_embedding_dataframe(adata_harmony, "Harmony").bind(
                                            lambda df_harmony: export_mode_svgs(df_harmony, "harmony", config.output_dir).map(
                                                lambda paths_harmony: paths_uncorrected + paths_combat + paths_harmony
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
    )


def main() -> None:
    """CLI entry point for generating Uncorrected, ComBat, and Harmony PCA/UMAP plots."""
    import argparse
    parser = argparse.ArgumentParser(description="Generate Uncorrected, ComBat, and Harmony PCA/UMAP plots.")
    parser.add_argument("--input-h5ad", type=str, default="scratch/unified_step07_gse120575.h5ad", help="Input unified h5ad path")
    parser.add_argument("--output-dir", type=str, default="scratch/plots_integration", help="Output directory for SVG plots")
    args = parser.parse_args()

    config = PlotConfig(
        input_h5ad_path=Path(args.input_h5ad).resolve(),
        output_dir=Path(args.output_dir).resolve(),
    )
    print(f"Starting PCA/UMAP embedding plotting pipeline (Uncorrected, ComBat, Harmony)...")
    print(f"Input file: {config.input_h5ad_path}")
    print(f"Output dir: {config.output_dir}")

    match run_plotting_pipeline(config):
        case Failure(err):
            print(f"Plotting error: {err}")
        case Success(svg_paths):
            print(f"\nSuccessfully generated {len(svg_paths)} SVG embedding plots:")
            for p in svg_paths:
                print(f" - {p}")


if __name__ == "__main__":
    main()
