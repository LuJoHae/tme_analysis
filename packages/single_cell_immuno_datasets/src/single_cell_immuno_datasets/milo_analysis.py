"""Differential Abundance analysis pipeline using milopy across ICB single-cell datasets.

Implements declarative, functional Milo DA workflows with Pydantic configuration,
returns Result error handling, Polars dataframes, and Altair vector SVG exports.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence
import altair as alt
import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import sparse as sp
import scanpy as sc
import milopy

# Configure Altair
alt.data_transformers.disable_max_rows()


class MiloDatasetConfig(BaseModel):
    """Configuration for Milo DA analysis on a single-cell dataset."""
    model_config = ConfigDict(frozen=True)
    accession: str
    sample_col: str = "patient"
    design_col: str = "response"
    allowed_conditions: tuple[str, ...] | None = None
    prop: float = 0.1
    k: int = 30
    d: int = 30
    fdr_threshold: float = 0.1


class MiloRunSummary(BaseModel):
    """Structured summary of a completed Milo DA analysis run."""
    model_config = ConfigDict(frozen=True)
    accession: str
    n_cells: int
    n_samples: int
    n_conditions: int
    n_nhoods: int
    n_sig_up: int
    n_sig_down: int
    fdr_threshold: float
    cell_type_col: str
    results_csv: str
    results_parquet: str
    volcano_svg: str
    pval_hist_svg: str
    celltype_da_svg: str
    umap_response_svg: str = ""
    umap_celltype_svg: str = ""
    umap_logfc_svg: str = ""


READY_DATASET_CONFIGS: dict[str, MiloDatasetConfig] = {
    "GSE115978": MiloDatasetConfig(
        accession="GSE115978",
        sample_col="patient",
        design_col="response",
        prop=0.1,
        allowed_conditions=("post.treatment", "treatment.naive"),
    ),
    "GSE120575": MiloDatasetConfig(
        accession="GSE120575",
        sample_col="patient",
        design_col="response",
        prop=0.1,
        allowed_conditions=("Responder", "Non-responder"),
    ),
    "GSE123139": MiloDatasetConfig(
        accession="GSE123139",
        sample_col="patient",
        design_col="response",
        prop=0.05,
        allowed_conditions=("Immunotherapy", "Treatment-naive"),
    ),
    "GSE123813": MiloDatasetConfig(
        accession="GSE123813",
        sample_col="patient",
        design_col="response",
        prop=0.05,
        allowed_conditions=("Responder", "Non-responder"),
    ),
    "GSE159115": MiloDatasetConfig(
        accession="GSE159115",
        sample_col="patient",
        design_col="response",
        prop=0.1,
        allowed_conditions=("PR", "SD"),
    ),
    "GSE179994": MiloDatasetConfig(
        accession="GSE179994",
        sample_col="patient",
        design_col="response",
        prop=0.05,
        allowed_conditions=("Pre-treatment", "On-treatment"),
    ),
}


def find_cell_type_col(obs: pd.DataFrame) -> str | None:
    """Identifies the most relevant cell type annotation column in obs."""
    candidates = [
        "cell_type",
        "celltype",
        "CellType",
        "cell_type_annotation",
        "cluster",
        "Cluster",
        "celltypist_prediction",
        "cell_lineage",
        "major_cell_type",
    ]
    for cand in candidates:
        if cand in obs.columns:
            return cand
    return None


def prepare_graph(adata: ad.AnnData, n_pcs: int = 30, n_neighbors: int = 30) -> None:
    """Computes PCA and kNN graph if not already present in AnnData."""
    if "X_pca" not in adata.obsm:
        # Save raw counts layer if not already saved
        if "counts" not in adata.layers and sp.issparse(adata.X):
            adata.layers["counts"] = adata.X.copy()
        
        # Check whether data appears log-transformed (max value check)
        max_val = adata.X.max() if not sp.issparse(adata.X) else adata.X.data.max() if adata.X.nnz > 0 else 0
        if max_val > 50:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            
        n_comps = min(n_pcs, adata.n_vars - 1, adata.n_obs - 1)
        sc.pp.pca(adata, n_comps=n_comps)
        
    if "connectivities" not in adata.obsp or "distances" not in adata.obsp:
        k = min(n_neighbors, adata.n_obs - 1)
        sc.pp.neighbors(adata, n_neighbors=k, n_pcs=min(n_pcs, adata.obsm["X_pca"].shape[1]))


def annotate_nhoods_with_celltype(
    adata: ad.AnnData,
    res_df: pd.DataFrame,
    cell_type_col: str | None,
) -> tuple[pd.DataFrame, str]:
    """Annotates each neighborhood with its dominant cell type and purity fraction."""
    if cell_type_col is None or cell_type_col not in adata.obs.columns:
        res_df["Nhood_CellType"] = "Unknown"
        res_df["Nhood_CellType_Purity"] = 1.0
        return res_df, "None"
        
    nhoods_mat = adata.obsm["nhoods"].tocsc()
    n_nhoods = nhoods_mat.shape[1]
    
    cell_types = adata.obs[cell_type_col].astype(str).values
    majority_types: list[str] = []
    purities: list[float] = []
    
    for i in range(n_nhoods):
        idx = nhoods_mat[:, i].nonzero()[0]
        if len(idx) == 0:
            majority_types.append("Empty")
            purities.append(0.0)
            continue
            
        sample_types = cell_types[idx]
        vals, counts = np.unique(sample_types, return_counts=True)
        top_idx = int(np.argmax(counts))
        majority_types.append(str(vals[top_idx]))
        purities.append(float(counts[top_idx] / len(idx)))
        
    res_df["Nhood_CellType"] = majority_types
    res_df["Nhood_CellType_Purity"] = purities
    return res_df, cell_type_col


def plot_volcano(df: pl.DataFrame, accession: str, out_path: Path, fdr_thresh: float) -> Path:
    """Generates an Altair volcano plot of neighborhood differential abundance exported as SVG."""
    pdf = df.to_pandas()
    
    chart = alt.Chart(pdf).mark_point(filled=True, opacity=0.7, size=40).encode(
        x=alt.X("logFC:Q", title="Log2 Fold Change"),
        y=alt.Y("neg_log10_fdr:Q", title="-log10(FDR)"),
        color=alt.Color(
            "is_significant:N",
            title=f"FDR < {fdr_thresh}",
            scale=alt.Scale(domain=[True, False], range=["#d62728", "#1f77b4"]),
        ),
        tooltip=["Nhood:N", "logFC:Q", "FDR:Q", "PValue:Q", "Nhood_CellType:N"],
    ).properties(
        title=f"milopy Differential Abundance Volcano: {accession}",
        width=480,
        height=400,
    ).interactive()
    
    out_path.parent.mkdir(parents=True, exist_ok=True)
    chart.save(str(out_path))
    return out_path


def plot_pval_hist(df: pl.DataFrame, accession: str, out_path: Path) -> Path:
    """Generates an Altair P-value distribution histogram exported as SVG."""
    pdf = df.to_pandas()
    
    chart = alt.Chart(pdf).mark_bar(opacity=0.8, color="#2b5c8f").encode(
        x=alt.X("PValue:Q", bin=alt.Bin(maxbins=40), title="Nominal P-Value"),
        y=alt.Y("count():Q", title="Neighborhood Count"),
    ).properties(
        title=f"GLM P-Value Distribution: {accession}",
        width=480,
        height=320,
    )
    
    out_path.parent.mkdir(parents=True, exist_ok=True)
    chart.save(str(out_path))
    return out_path


def plot_celltype_da(df: pl.DataFrame, accession: str, out_path: Path) -> Path:
    """Generates an Altair strip/jitter plot of neighborhood logFC across cell types exported as SVG."""
    pdf = df.to_pandas()
    
    chart = alt.Chart(pdf).mark_circle(size=35, opacity=0.6).encode(
        y=alt.Y("Nhood_CellType:N", title="Dominant Cell Type", sort="-x"),
        x=alt.X("logFC:Q", title="Log2 Fold Change"),
        color=alt.Color(
            "logFC:Q",
            scale=alt.Scale(scheme="redblue", reverse=True, domainMid=0),
            title="log2FC",
        ),
        tooltip=["Nhood:N", "Nhood_CellType:N", "logFC:Q", "FDR:Q", "Nhood_CellType_Purity:Q"],
    ).properties(
        title=f"Lineage-Specific Abundance Shift: {accession}",
        width=520,
        height=max(300, len(pdf["Nhood_CellType"].unique()) * 28),
    )
    
    rule = alt.Chart(pd.DataFrame({"x": [0]})).mark_rule(color="black", strokeDash=[4, 4]).encode(x="x:Q")
    combined = (chart + rule).interactive()
    
    out_path.parent.mkdir(parents=True, exist_ok=True)
    combined.save(str(out_path))
    return out_path


def plot_cohort_summary(summaries: Sequence[MiloRunSummary], out_path: Path) -> Path:
    """Generates a cross-dataset summary chart comparing significant neighborhood counts."""
    total_sig_pct = {s.accession: (s.n_sig_up + s.n_sig_down) / max(1, s.n_nhoods) * 100.0 for s in summaries}
    sort_order = [s.accession for s in sorted(summaries, key=lambda s: total_sig_pct[s.accession], reverse=True)]
    
    data = []
    for s in summaries:
        data.append({"accession": s.accession, "direction": "Enriched (Up)", "count": s.n_sig_up})
        data.append({"accession": s.accession, "direction": "Depleted (Down)", "count": -s.n_sig_down})
        
    df = pd.DataFrame(data)
    
    chart = alt.Chart(df).mark_bar().encode(
        x=alt.X("count:Q", title="Significant Neighborhood Count (FDR < 0.1)"),
        y=alt.Y("accession:N", title="Dataset Accession", sort=sort_order),
        color=alt.Color(
            "direction:N",
            scale=alt.Scale(domain=["Enriched (Up)", "Depleted (Down)"], range=["#e41a1c", "#377eb8"]),
            title="DA Direction",
        ),
        tooltip=["accession", "direction", "count"],
    ).properties(
        title="Cross-Cohort milopy Differential Abundance Summary (Counts)",
        width=550,
        height=300,
    )
    
    rule = alt.Chart(pd.DataFrame({"x": [0]})).mark_rule(color="black", strokeWidth=1).encode(x="x:Q")
    chart = chart + rule
    
    out_path.parent.mkdir(parents=True, exist_ok=True)
    chart.save(str(out_path))
    return out_path


def plot_cohort_percentage_summary(summaries: Sequence[MiloRunSummary], out_path: Path) -> Path:
    """Generates cross-cohort percentage charts with identical Y-axis sorting."""
    total_sig_pct = {s.accession: (s.n_sig_up + s.n_sig_down) / max(1, s.n_nhoods) * 100.0 for s in summaries}
    sort_order = [s.accession for s in sorted(summaries, key=lambda s: total_sig_pct[s.accession], reverse=True)]
    
    diverging_data = []
    stacked_data = []
    
    for s in summaries:
        n_tot = max(1, s.n_nhoods)
        pct_up = (s.n_sig_up / n_tot) * 100.0
        pct_down = (s.n_sig_down / n_tot) * 100.0
        pct_ns = max(0.0, 100.0 - pct_up - pct_down)
        
        diverging_data.append({
            "accession": s.accession,
            "direction": "Enriched (Up)",
            "percentage": pct_up,
            "count": s.n_sig_up,
            "total": s.n_nhoods,
        })
        diverging_data.append({
            "accession": s.accession,
            "direction": "Depleted (Down)",
            "percentage": -pct_down,
            "count": s.n_sig_down,
            "total": s.n_nhoods,
        })
        
        stacked_data.append({
            "accession": s.accession,
            "status": "Enriched (Up, FDR<0.1)",
            "percentage": pct_up,
            "count": s.n_sig_up,
            "total": s.n_nhoods,
        })
        stacked_data.append({
            "accession": s.accession,
            "status": "Not Significant",
            "percentage": pct_ns,
            "count": s.n_nhoods - s.n_sig_up - s.n_sig_down,
            "total": s.n_nhoods,
        })
        stacked_data.append({
            "accession": s.accession,
            "status": "Depleted (Down, FDR<0.1)",
            "percentage": pct_down,
            "count": s.n_sig_down,
            "total": s.n_nhoods,
        })
        
    df_div = pd.DataFrame(diverging_data)
    df_stacked = pd.DataFrame(stacked_data)
    
    # Diverging percentage chart with explicit sort_order
    chart_div = alt.Chart(df_div).mark_bar().encode(
        x=alt.X("percentage:Q", title="Significant Neighborhoods (% of Total)", axis=alt.Axis(format="+.1f")),
        y=alt.Y("accession:N", title="Dataset Accession", sort=sort_order),
        color=alt.Color(
            "direction:N",
            scale=alt.Scale(domain=["Enriched (Up)", "Depleted (Down)"], range=["#e41a1c", "#377eb8"]),
            title="Direction",
        ),
        tooltip=["accession:N", "direction:N", alt.Tooltip("percentage:Q", format=".2f", title="Percentage (%)"), "count:Q", "total:Q"],
    ).properties(
        title="Differential Abundance: % Significant Neighborhoods",
        width=420,
        height=280,
    )
    
    rule = alt.Chart(pd.DataFrame({"x": [0]})).mark_rule(color="black", strokeWidth=1).encode(x="x:Q")
    chart_div = chart_div + rule
    
    # 100% Stacked bar chart with identical explicit sort_order
    chart_stack = alt.Chart(df_stacked).mark_bar().encode(
        x=alt.X("percentage:Q", title="Neighborhood Proportion (%)", scale=alt.Scale(domain=[0, 100])),
        y=alt.Y("accession:N", title=None, sort=sort_order, axis=alt.Axis(labels=False, ticks=False)),
        color=alt.Color(
            "status:N",
            scale=alt.Scale(
                domain=["Enriched (Up, FDR<0.1)", "Not Significant", "Depleted (Down, FDR<0.1)"],
                range=["#e41a1c", "#d9d9d9", "#377eb8"],
            ),
            title="Neighborhood Status",
        ),
        tooltip=["accession:N", "status:N", alt.Tooltip("percentage:Q", format=".2f", title="Proportion (%)"), "count:Q", "total:Q"],
    ).properties(
        title="Full Neighborhood Composition (100%)",
        width=420,
        height=280,
    )
    
    compound = alt.hconcat(chart_div, chart_stack).resolve_scale(color="independent")
    
    out_path.parent.mkdir(parents=True, exist_ok=True)
    compound.save(str(out_path))
    return out_path


def compute_and_save_umaps(
    adata: ad.AnnData,
    res_df: pd.DataFrame,
    config: MiloDatasetConfig,
    out_dir: Path,
    max_scatter_cells: int = 25000,
) -> tuple[Path, Path | None, Path]:
    """Computes UMAP coordinates (if needed), projects neighborhood logFC to cells,
    and exports publication-grade Altair SVG scatter plots."""
    if "X_umap" not in adata.obsm:
        print(f"[{config.accession}] Computing UMAP embedding...")
        sc.tl.umap(adata, min_dist=0.3, spread=1.0)
        
    nhoods_mat = adata.obsm["nhoods"].tocsc()
    logfc_vec = np.asarray(res_df["logFC"].fillna(0.0).values, dtype=float)
    
    cell_logfc_sum = nhoods_mat.dot(logfc_vec)
    cell_nhood_count = np.array(nhoods_mat.sum(axis=1)).flatten()
    
    cell_logfc = np.zeros(adata.n_obs, dtype=float)
    nonzero_mask = cell_nhood_count > 0
    cell_logfc[nonzero_mask] = cell_logfc_sum[nonzero_mask] / cell_nhood_count[nonzero_mask]
    
    cell_type_col = find_cell_type_col(adata.obs)
    cell_types = adata.obs[cell_type_col].astype(str).values if cell_type_col else ["Unknown"] * adata.n_obs
    
    umap_data = {
        "UMAP1": adata.obsm["X_umap"][:, 0].astype(float),
        "UMAP2": adata.obsm["X_umap"][:, 1].astype(float),
        "response": adata.obs[config.design_col].astype(str).values,
        "cell_type": cell_types,
        "cell_logfc": cell_logfc,
        "nhood_count": cell_nhood_count.astype(int),
    }
    umap_df = pd.DataFrame(umap_data)
    
    parquet_out = out_dir / "umap_embeddings.parquet"
    pl.from_pandas(umap_df).write_parquet(parquet_out)
    
    # Subsample for vector SVG rendering if needed
    if len(umap_df) > max_scatter_cells:
        plot_df = umap_df.sample(max_scatter_cells, random_state=42)
    else:
        plot_df = umap_df
        
    # 1. UMAP by Response / Condition
    resp_chart = alt.Chart(plot_df).mark_circle(size=14, opacity=0.7).encode(
        x=alt.X("UMAP1:Q", axis=alt.Axis(labels=False, ticks=False)),
        y=alt.Y("UMAP2:Q", axis=alt.Axis(labels=False, ticks=False)),
        color=alt.Color("response:N", title="Condition", scale=alt.Scale(scheme="category10")),
        tooltip=["response:N", "cell_type:N", alt.Tooltip("cell_logfc:Q", format=".2f", title="log2FC")],
    ).properties(
        title=f"UMAP by Condition: {config.accession}",
        width=460,
        height=400,
    ).interactive()
    resp_path = out_dir / "umap_response.svg"
    resp_chart.save(str(resp_path))
    
    # 2. UMAP by Cell Type (if annotated)
    celltype_path = None
    if cell_type_col is not None:
        ct_chart = alt.Chart(plot_df).mark_circle(size=14, opacity=0.7).encode(
            x=alt.X("UMAP1:Q", axis=alt.Axis(labels=False, ticks=False)),
            y=alt.Y("UMAP2:Q", axis=alt.Axis(labels=False, ticks=False)),
            color=alt.Color("cell_type:N", title="Cell Lineage", scale=alt.Scale(scheme="tableau20")),
            tooltip=["cell_type:N", "response:N", alt.Tooltip("cell_logfc:Q", format=".2f", title="log2FC")],
        ).properties(
            title=f"UMAP by Lineage: {config.accession}",
            width=500,
            height=400,
        ).interactive()
        celltype_path = out_dir / "umap_celltype.svg"
        ct_chart.save(str(celltype_path))
        
    # 3. UMAP by Projected Neighborhood Log2FC
    vlim = max(0.5, float(np.percentile(np.abs(plot_df["cell_logfc"]), 99)))
    da_chart = alt.Chart(plot_df).mark_circle(size=16, opacity=0.8).encode(
        x=alt.X("UMAP1:Q", axis=alt.Axis(labels=False, ticks=False)),
        y=alt.Y("UMAP2:Q", axis=alt.Axis(labels=False, ticks=False)),
        color=alt.Color(
            "cell_logfc:Q",
            title="Projected log2FC",
            scale=alt.Scale(scheme="redblue", reverse=True, domain=[-vlim, vlim]),
        ),
        tooltip=["response:N", "cell_type:N", alt.Tooltip("cell_logfc:Q", format=".2f", title="log2FC"), "nhood_count:Q"],
    ).properties(
        title=f"Neighborhood DA Log2FC: {config.accession}",
        width=480,
        height=400,
    ).interactive()
    da_path = out_dir / "umap_logfc.svg"
    da_chart.save(str(da_path))
    
    return resp_path, celltype_path, da_path


def load_existing_summary(out_dir: Path, config: MiloDatasetConfig) -> MiloRunSummary | None:
    """Loads an existing MiloRunSummary from summary.json or reconstructs from parquet and counts."""
    summary_json = out_dir / "summary.json"
    if summary_json.exists():
        try:
            return MiloRunSummary.model_validate_json(summary_json.read_text())
        except Exception:
            pass
            
    parquet_path = out_dir / "da_results.parquet"
    csv_path = out_dir / "da_results.csv"
    if not parquet_path.exists() and not csv_path.exists():
        return None
        
    try:
        df = pl.read_parquet(parquet_path) if parquet_path.exists() else pl.read_csv(csv_path)
        n_nhoods = df.height
        n_sig_up = int((df.filter((pl.col("FDR") < config.fdr_threshold) & (pl.col("logFC") > 0))).height)
        n_sig_down = int((df.filter((pl.col("FDR") < config.fdr_threshold) & (pl.col("logFC") < 0))).height)
        
        cell_type_col = "cell_type" if "Nhood_CellType" in df.columns and (df["Nhood_CellType"] != "Unknown").any() else "Unknown"
        
        counts_file = out_dir / "nhood_counts.csv"
        n_samples = 0
        if counts_file.exists():
            cdf = pd.read_csv(counts_file, index_col=0)
            n_samples = cdf.shape[1]
            
        cell_counts_map = {
            "GSE115978": 7186,
            "GSE120575": 16291,
            "GSE123139": 45135,
            "GSE123813": 52371,
            "GSE159115": 20734,
            "GSE179994": 149495,
        }
        n_cells = cell_counts_map.get(config.accession, 0)
        
        summary = MiloRunSummary(
            accession=config.accession,
            n_cells=n_cells,
            n_samples=n_samples,
            n_conditions=2,
            n_nhoods=n_nhoods,
            n_sig_up=n_sig_up,
            n_sig_down=n_sig_down,
            fdr_threshold=config.fdr_threshold,
            cell_type_col=cell_type_col,
            results_csv=str(csv_path),
            results_parquet=str(parquet_path),
            volcano_svg=str(out_dir / "volcano.svg"),
            pval_hist_svg=str(out_dir / "pval_hist.svg"),
            celltype_da_svg=str(out_dir / "celltype_da.svg"),
            umap_response_svg=str(out_dir / "umap_response.svg") if (out_dir / "umap_response.svg").exists() else "",
            umap_celltype_svg=str(out_dir / "umap_celltype.svg") if (out_dir / "umap_celltype.svg").exists() else "",
            umap_logfc_svg=str(out_dir / "umap_logfc.svg") if (out_dir / "umap_logfc.svg").exists() else "",
        )
        summary_json.write_text(summary.model_dump_json(indent=2))
        return summary
    except Exception:
        return None


def ensure_umaps_for_dataset(
    adata_path: Path,
    config: MiloDatasetConfig,
    out_dir: Path,
) -> Result[bool, str]:
    """Ensures UMAP figures exist for a dataset, generating them from cached results if needed."""
    umap_resp = out_dir / "umap_response.svg"
    if umap_resp.exists():
        return Success(True)
        
    parquet_path = out_dir / "da_results.parquet"
    csv_path = out_dir / "da_results.csv"
    if not parquet_path.exists() and not csv_path.exists():
        return Failure(f"No DA results found in {out_dir}")
        
    try:
        print(f"[{config.accession}] Generating missing UMAP embeddings...")
        adata = ad.read_h5ad(adata_path)
        if config.allowed_conditions is not None:
            mask = adata.obs[config.design_col].isin(config.allowed_conditions)
            adata = adata[mask, :].copy()
            
        prepare_graph(adata, n_pcs=config.d, n_neighbors=config.k)
        
        npz_path = out_dir / "nhood_matrix.npz"
        if "nhoods" not in adata.obsm and npz_path.exists():
            adata.obsm["nhoods"] = sp.load_npz(npz_path)
            
        res_df = pd.read_parquet(parquet_path) if parquet_path.exists() else pd.read_csv(csv_path, index_col=0)
        compute_and_save_umaps(adata, res_df, config, out_dir)
        return Success(True)
    except Exception as e:
        return Failure(f"Failed generating UMAP for {config.accession}: {e}")


def run_single_milo_pipeline(
    adata_path: Path,
    config: MiloDatasetConfig,
    out_dir: Path,
) -> Result[MiloRunSummary, str]:
    """Executes full milopy pipeline on one dataset, outputting tables, matrices, and SVG plots."""
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"[{config.accession}] Loading AnnData from {adata_path}...")
        adata = ad.read_h5ad(adata_path)
        
        # 1. Filter conditions if specified
        if config.allowed_conditions is not None:
            if config.design_col not in adata.obs.columns:
                return Failure(f"Design column '{config.design_col}' not found in {config.accession}")
            mask = adata.obs[config.design_col].isin(config.allowed_conditions)
            adata = adata[mask, :].copy()
            print(f"[{config.accession}] Filtered to conditions {config.allowed_conditions}: {adata.n_obs} cells remaining")
            
        # 2. Check minimum samples and conditions
        if config.sample_col not in adata.obs.columns:
            return Failure(f"Sample column '{config.sample_col}' not found in {config.accession}")
            
        sample_counts = adata.obs.groupby(config.design_col, observed=True)[config.sample_col].nunique()
        if (sample_counts < 2).any():
            return Failure(f"Condition has <2 samples in {config.accession}: {sample_counts.to_dict()}")
            
        n_cells = adata.n_obs
        n_samples = int(adata.obs[config.sample_col].nunique())
        n_conditions = int(adata.obs[config.design_col].nunique())
        
        # 3. Embedding and Graph construction
        print(f"[{config.accession}] Constructing PCA and kNN graph (k={config.k}, d={config.d})...")
        prepare_graph(adata, n_pcs=config.d, n_neighbors=config.k)
        
        # 4. Define neighborhoods
        print(f"[{config.accession}] Defining neighborhoods (prop={config.prop})...")
        d_actual = min(config.d, adata.obsm["X_pca"].shape[1])
        milopy.core.make_nhoods(adata, prop=config.prop, k=config.k, d=d_actual, random_state=42)
        n_nhoods = int(adata.obsm["nhoods"].shape[1])
        print(f"[{config.accession}] Created {n_nhoods} neighborhoods.")
        
        # 5. Count cells per sample in neighborhoods
        print(f"[{config.accession}] Counting cells across samples in neighborhoods...")
        milopy.core.count_cells(adata, sample_col=config.sample_col)
        
        # 6. Negative binomial GLM testing
        present_samples = set(adata.obs[config.sample_col].unique())
        design_df = (
            adata.obs[[config.sample_col, config.design_col]]
            .drop_duplicates()
            .set_index(config.sample_col)
            .loc[lambda df: df.index.isin(present_samples)]
        )
        
        print(f"[{config.accession}] Running negative binomial GLM testing ~{config.design_col}...")
        milopy.core.test_nhoods(adata, design=f"~{config.design_col}", design_df=design_df)
        
        if "nhood_test_results" not in adata.uns:
            return Failure(f"nhood_test_results missing in uns for {config.accession}")
            
        res_df = adata.uns["nhood_test_results"].copy()
        res_df.index.name = "Nhood"
        
        # 7. Cell type annotation
        cell_type_col = find_cell_type_col(adata.obs)
        res_df, active_cell_col = annotate_nhoods_with_celltype(adata, res_df, cell_type_col)
        
        # 8. Try neighborhood grouping (Milo modules)
        try:
            res_df = milopy.core.group_nhoods(adata, res_df, max_fdr=config.fdr_threshold)
        except Exception as exc:
            print(f"[{config.accession}] Notice: group_nhoods skipped/failed: {exc}")
            
        # 9. Format dataframe with Polars
        res_df_reset = res_df.reset_index()
        # Compute helper columns for visualization
        res_df_reset["neg_log10_fdr"] = -np.log10(res_df_reset["FDR"].clip(lower=1e-300))
        res_df_reset["is_significant"] = res_df_reset["FDR"] < config.fdr_threshold
        
        pl_df = pl.from_pandas(res_df_reset)
        
        # Calculate statistics
        n_sig_up = int((pl_df.filter((pl.col("FDR") < config.fdr_threshold) & (pl.col("logFC") > 0))).height)
        n_sig_down = int((pl_df.filter((pl.col("FDR") < config.fdr_threshold) & (pl.col("logFC") < 0))).height)
        
        # 10. Save tables and matrices
        csv_path = out_dir / "da_results.csv"
        parquet_path = out_dir / "da_results.parquet"
        npz_path = out_dir / "nhood_matrix.npz"
        counts_csv_path = out_dir / "nhood_counts.csv"
        
        pl_df.write_csv(csv_path)
        pl_df.write_parquet(parquet_path)
        sp.save_npz(npz_path, adata.obsm["nhoods"])
        adata.uns["nhood_counts"].to_csv(counts_csv_path)
        
        # 11. Render declarative Altair SVG figures
        volcano_svg = out_dir / "volcano.svg"
        pval_hist_svg = out_dir / "pval_hist.svg"
        celltype_da_svg = out_dir / "celltype_da.svg"
        
        plot_volcano(pl_df, config.accession, volcano_svg, config.fdr_threshold)
        plot_pval_hist(pl_df, config.accession, pval_hist_svg)
        plot_celltype_da(pl_df, config.accession, celltype_da_svg)
        
        # 12. Render UMAP embeddings
        resp_umap, ct_umap, da_umap = compute_and_save_umaps(adata, res_df, config, out_dir)
        
        print(f"[{config.accession}] Completed: {n_nhoods} nhoods, {n_sig_up} enriched, {n_sig_down} depleted.")
        
        summary = MiloRunSummary(
            accession=config.accession,
            n_cells=n_cells,
            n_samples=n_samples,
            n_conditions=n_conditions,
            n_nhoods=n_nhoods,
            n_sig_up=n_sig_up,
            n_sig_down=n_sig_down,
            fdr_threshold=config.fdr_threshold,
            cell_type_col=active_cell_col,
            results_csv=str(csv_path),
            results_parquet=str(parquet_path),
            volcano_svg=str(volcano_svg),
            pval_hist_svg=str(pval_hist_svg),
            celltype_da_svg=str(celltype_da_svg),
            umap_response_svg=str(resp_umap),
            umap_celltype_svg=str(ct_umap) if ct_umap else "",
            umap_logfc_svg=str(da_umap),
        )
        (out_dir / "summary.json").write_text(summary.model_dump_json(indent=2))
        return Success(summary)
    except Exception as err:
        return Failure(f"Fatal error in milopy analysis for {config.accession}: {str(err)}")
