#!/usr/bin/env python3
"""Synthetic Single-Cell UMAP Visualizer: Milopy Differential Abundance.

Simulates single-cell RNA-seq counts across K ground-truth clusters, assigns cells
to 30 clinical patients (15 Responders vs 15 Non-Responders), executes Milopy DA testing,
projects neighborhood effect sizes to single cells, computes 2D UMAP embeddings, and
exports publication-grade Nature Methods compliant Altair vector SVGs and Parquet tables.

Adheres strictly to functional Python, immutability, Result monads, Polars, and Altair.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

import altair as alt  # type: ignore
import anndata as ad  # type: ignore
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
from scipy import sparse as sp  # type: ignore
import scanpy as sc  # type: ignore
import vl_convert as vlc  # type: ignore

try:
    import milopy  # type: ignore
    HAS_MILOPY = True
except ImportError:
    HAS_MILOPY = False

# Configure Altair for publication-grade exports
alt.data_transformers.disable_max_rows()

from plotting_utils import (
    COLOR_HAIRLINE,
    COLOR_LIGHT_GREY,
    COLOR_NEUTRAL_GREY,
    COLOR_TEXT_PRIMARY,
    OKABE_BLUE,
    OKABE_BLUISH_GREEN,
    OKABE_ORANGE,
    OKABE_REDDISH_PURPLE,
    OKABE_SKY_BLUE,
    OKABE_VERMILION,
    OKABE_YELLOW,
)

CLUSTER_PALETTE: tuple[str, ...] = (
    OKABE_BLUE,
    OKABE_ORANGE,
    OKABE_BLUISH_GREEN,
    OKABE_YELLOW,
    OKABE_SKY_BLUE,
    OKABE_REDDISH_PURPLE,
)


class UmapVisualizerConfig(BaseModel):
    """Immutable configuration for synthetic single-cell UMAP generation and visualization."""

    model_config = ConfigDict(frozen=True)

    n_cells: int = 3000
    n_genes: int = 300
    n_clusters: int = 6
    n_patients: int = 30
    cluster_width: float = 0.85  # Drastically increases cluster dispersion / width
    marker_overlap: float = 0.90  # Kernel bandwidth for cross-cluster marker program overlap
    bg_lambda: float = 0.80  # Baseline ambient expression
    umap_min_dist: float = 0.40
    umap_spread: float = 1.20
    match_prob: float = 0.85
    variable_ratios: bool = True
    min_prob: float = 0.10
    max_prob: float = 0.90
    leiden_resolution: float = 0.5
    fdr_threshold: float = 0.10
    milo_prop: float = 0.15
    seed: int = 42
    out_dir: Path = Path("output/synthetic_benchmark")
    results_dir: Path = Path("results/synthetic_benchmark")
    input_parquet: Path | None = None


def parse_args() -> UmapVisualizerConfig:
    """Parses command-line arguments into an immutable configuration model."""
    parser = argparse.ArgumentParser(
        description="Synthetic scRNA-seq UMAP visualizer for Milopy differential abundance results."
    )
    parser.add_argument("--n-cells", type=int, default=3000, help="Total synthetic cell count")
    parser.add_argument("--n-genes", type=int, default=300, help="Total gene count")
    parser.add_argument("--n-clusters", type=int, default=6, help="Number of ground-truth clusters")
    parser.add_argument("--n-patients", type=int, default=30, help="Total number of patients (half R, half NR)")
    parser.add_argument("--cluster-width", type=float, default=0.85, help="Cluster spatial width / dispersion (controls overlap)")
    parser.add_argument("--marker-overlap", type=float, default=0.90, help="Bandwidth for cross-cluster marker expression overlap")
    parser.add_argument("--bg-lambda", type=float, default=0.80, help="Ambient baseline expression rate")
    parser.add_argument("--umap-min-dist", type=float, default=0.40, help="UMAP minimum distance parameter")
    parser.add_argument("--umap-spread", type=float, default=1.20, help="UMAP spread parameter")
    parser.add_argument("--fdr-threshold", type=float, default=0.10, help="FDR threshold for DA significance")
    parser.add_argument("--resolution", type=float, default=0.5, help="Leiden clustering resolution")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument("--out-dir", type=Path, default=Path("output/synthetic_benchmark"), help="Output directory for parquets")
    parser.add_argument("--results-dir", type=Path, default=Path("results/synthetic_benchmark"), help="Results directory for SVGs")
    parser.add_argument("--input-parquet", type=Path, default=None, help="Optional pre-computed cell-level parquet to visualize directly")
    args = parser.parse_args()

    return UmapVisualizerConfig(
        n_cells=args.n_cells,
        n_genes=args.n_genes,
        n_clusters=args.n_clusters,
        n_patients=args.n_patients,
        cluster_width=args.cluster_width,
        marker_overlap=args.marker_overlap,
        bg_lambda=args.bg_lambda,
        umap_min_dist=args.umap_min_dist,
        umap_spread=args.umap_spread,
        fdr_threshold=args.fdr_threshold,
        leiden_resolution=args.resolution,
        seed=args.seed,
        out_dir=args.out_dir,
        results_dir=args.results_dir,
        input_parquet=args.input_parquet,
    )


# =============================================================================
# 1. Pure Simulation Engine
# =============================================================================

def generate_synthetic_scrna(
    config: UmapVisualizerConfig,
    rng: np.random.Generator,
) -> ad.AnnData:
    """Pure function generating synthetic single-cell expression matrix with wide, overlapping clusters."""
    n_cells = config.n_cells
    n_genes = config.n_genes
    k = config.n_clusters

    cells_per_cluster = n_cells // k
    cluster_assignments = np.concatenate([
        np.full(cells_per_cluster, fill_value=i, dtype=np.int32)
        for i in range(k)
    ])
    remainder = n_cells - len(cluster_assignments)
    if remainder > 0:
        cluster_assignments = np.concatenate([
            cluster_assignments,
            rng.choice(k, size=remainder).astype(np.int32),
        ])

    # Cluster positions on a circular manifold
    cluster_angles = np.linspace(0, 2 * np.pi, k, endpoint=False)
    cluster_centers = np.column_stack([np.cos(cluster_angles), np.sin(cluster_angles)])

    # Drastically disperse cell latent positions (large cluster_width -> wide, overlapping clouds)
    cell_latent = np.zeros((n_cells, 2), dtype=np.float32)
    for i in range(k):
        mask = (cluster_assignments == i)
        cell_latent[mask] = cluster_centers[i] + rng.normal(0, config.cluster_width, size=(int(mask.sum()), 2))

    # Continuous marker activation kernel: cells near cluster boundaries express shared markers
    dists = np.linalg.norm(cell_latent[:, None, :] - cluster_centers[None, :, :], axis=2)
    weights = np.exp(-0.5 * (dists / config.marker_overlap) ** 2)
    weights /= weights.sum(axis=1, keepdims=True)

    # Base background: Poisson expression
    raw_counts = rng.poisson(lam=config.bg_lambda, size=(n_cells, n_genes)).astype(np.float32)

    # Expression of marker programs with continuous cross-cluster gradients
    genes_per_marker = max(10, n_genes // k)
    for j in range(k):
        start_g = j * genes_per_marker
        end_g = min(n_genes, (j + 1) * genes_per_marker)
        if end_g > start_g:
            marker_strength = weights[:, j][:, None] * 7.5
            noise = rng.negative_binomial(n=4, p=0.35, size=(n_cells, end_g - start_g))
            raw_counts[:, start_g:end_g] += (marker_strength * (noise + 1.0)).astype(np.float32)

    x_sparse = sp.csr_matrix(raw_counts)
    obs_names = [f"cell_{i:05d}" for i in range(n_cells)]
    var_names = [f"gene_{g:04d}" for g in range(n_genes)]

    adata = ad.AnnData(
        X=x_sparse,
        obs=pd.DataFrame(
            {"synthetic_ground_truth": cluster_assignments},
            index=obs_names,
        ),
        var=pd.DataFrame({"gene_id": var_names}, index=var_names),
    )
    return adata


def cluster_and_project_umap(
    adata: ad.AnnData,
    config: UmapVisualizerConfig,
    rng: np.random.Generator,
) -> Result[ad.AnnData, str]:
    """Computes normalization, PCA, kNN graph, Leiden clustering, and 2D UMAP coordinates."""
    try:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        n_comps = min(25, adata.n_vars - 1)
        sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=False, random_state=42, zero_center=False)
        sc.pp.neighbors(adata, n_neighbors=20, n_pcs=min(20, n_comps), random_state=42)
        sc.tl.leiden(adata, key_added="leiden", resolution=config.leiden_resolution, random_state=42)
        sc.tl.umap(adata, min_dist=config.umap_min_dist, spread=config.umap_spread, random_state=42)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Scanpy embedding pipeline failed: {exc}")

    clusters = sorted(adata.obs["leiden"].unique().tolist())
    n_cl = len(clusters)
    if n_cl < 2:
        return Failure(f"Clustering yielded only {n_cl} cluster(s); requires at least 2.")

    if config.variable_ratios:
        raw_probs = np.linspace(config.max_prob, config.min_prob, n_cl)
        rng.shuffle(raw_probs)
        cluster_probs = {cl: round(float(p), 3) for cl, p in zip(clusters, raw_probs, strict=True)}
    else:
        half = n_cl // 2
        binary_choices = [config.max_prob] * half + [config.min_prob] * (n_cl - half)
        rng.shuffle(binary_choices)
        cluster_probs = {cl: round(float(p), 3) for cl, p in zip(clusters, binary_choices, strict=True)}

    cell_probs = [cluster_probs[str(cl)] for cl in adata.obs["leiden"]]
    adata.obs["cluster_true_prob"] = cell_probs
    adata.obs["cluster_label"] = [1 if p >= 0.5 else 0 for p in cell_probs]
    return Success(adata)


def assign_patients(
    adata: ad.AnnData,
    n_patients: int,
    rng: np.random.Generator,
) -> Result[ad.AnnData, str]:
    """Assigns cells to responder and non-responder patients based on cluster probabilities."""
    if n_patients < 4 or n_patients % 2 != 0:
        return Failure("n_patients must be an even integer >= 4.")

    half_p = n_patients // 2
    true_patients = [f"patient_{j:02d}" for j in range(half_p)]
    false_patients = [f"patient_{j:02d}" for j in range(half_p, n_patients)]

    patient_assignments: list[str] = []
    patient_responses: list[str] = []
    cell_probs = adata.obs["cluster_true_prob"].to_numpy()

    for p_val in cell_probs:
        is_true = rng.random() < p_val
        if is_true:
            chosen_p = rng.choice(true_patients)
            resp = "Responder"
        else:
            chosen_p = rng.choice(false_patients)
            resp = "Non-Responder"

        patient_assignments.append(str(chosen_p))
        patient_responses.append(resp)

    adata.obs["patient"] = patient_assignments
    adata.obs["response"] = patient_responses
    return Success(adata)


def run_milopy_and_project_cells(
    adata: ad.AnnData,
    prop: float,
    fdr_threshold: float,
) -> Result[pl.DataFrame, str]:
    """Runs Milopy DA testing and projects neighborhood effect sizes and significance to single cells."""
    if not HAS_MILOPY:
        return Failure("milopy is not installed in the environment.")

    try:
        milopy.core.make_nhoods(adata, prop=prop, k=15, d=20, refined=True, random_state=42)
        milopy.core.count_cells(adata, sample_col="patient")

        design_df = (
            adata.obs[["patient", "response"]]
            .drop_duplicates()
            .set_index("patient")
        )
        milopy.core.test_nhoods(adata, design="~response", design_df=design_df)

        if "nhood_test_results" not in adata.uns:
            return Failure("nhood_test_results not found after test_nhoods.")

        nhood_res = adata.uns["nhood_test_results"]
        nhoods_mat = adata.obsm["nhoods"].tocsc()
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Milopy DA execution error: {exc}")

    logfc_col = "logFC" if "logFC" in nhood_res.columns else nhood_res.columns[0]
    nhood_lfc = nhood_res[logfc_col].fillna(0.0).to_numpy().astype(np.float64)
    fdr_col = "SpatialFDR" if "SpatialFDR" in nhood_res.columns else ("FDR" if "FDR" in nhood_res.columns else nhood_res.columns[1])
    nhood_fdr = nhood_res[fdr_col].fillna(1.0).to_numpy().astype(np.float64)

    # Cell-level weighted projection
    nhoods_sum = np.array(nhoods_mat.sum(axis=1)).flatten()
    cell_mask = nhoods_sum > 0

    cell_lfc = np.zeros(adata.n_obs, dtype=np.float64)
    cell_lfc[cell_mask] = np.array(nhoods_mat[cell_mask] @ nhood_lfc).flatten() / nhoods_sum[cell_mask]

    # Cell-level minimum FDR across neighborhoods
    # Convert CSC to COO for fast row-wise min
    nhoods_coo = nhoods_mat.tocoo()
    cell_min_fdr = np.ones(adata.n_obs, dtype=np.float64)
    for row_idx, col_idx in zip(nhoods_coo.row, nhoods_coo.col, strict=True):
        fdr_val = nhood_fdr[col_idx]
        if fdr_val < cell_min_fdr[row_idx]:
            cell_min_fdr[row_idx] = fdr_val

    # Status classification
    status_list: list[str] = []
    for lfc_val, fdr_val in zip(cell_lfc, cell_min_fdr, strict=True):
        if fdr_val < fdr_threshold and lfc_val > 0.05:
            status_list.append("DA+ (Responder Enriched)")
        elif fdr_val < fdr_threshold and lfc_val < -0.05:
            status_list.append("DA- (Non-Responder Enriched)")
        else:
            status_list.append("Not Significant")

    umap_coords = adata.obsm["X_umap"]
    cell_table = pl.DataFrame({
        "cell_id": adata.obs_names.tolist(),
        "UMAP_1": umap_coords[:, 0].tolist(),
        "UMAP_2": umap_coords[:, 1].tolist(),
        "cluster": [f"Cluster {c}" for c in adata.obs["leiden"]],
        "ground_truth": [f"State {g}" for g in adata.obs["synthetic_ground_truth"]],
        "patient": adata.obs["patient"].tolist(),
        "response": adata.obs["response"].tolist(),
        "cell_logFC": cell_lfc.tolist(),
        "spatial_fdr": cell_min_fdr.tolist(),
        "status": status_list,
    })

    return Success(cell_table)


# =============================================================================
# 2. Altair Vector SVG Visualizations (Nature Methods Wireframe Standards)
# =============================================================================

def build_cluster_chart(df: pl.DataFrame, width: int = 380, height: int = 320) -> alt.Chart:
    """Panel A: UMAP colored by discrete cell cluster / ground truth state."""
    pdf = df.to_pandas()
    clusters = sorted(pdf["cluster"].unique().tolist())
    palette = list(CLUSTER_PALETTE[: len(clusters)])

    chart = (
        alt.Chart(pdf)
        .mark_circle(size=28, opacity=0.85, stroke="#334155", strokeWidth=0.25)
        .encode(
            x=alt.X("UMAP_1:Q", title="UMAP 1", axis=alt.Axis(grid=False, tickCount=5)),
            y=alt.Y("UMAP_2:Q", title="UMAP 2", axis=alt.Axis(grid=False, tickCount=5)),
            color=alt.Color(
                "cluster:N",
                title="Cell State / Cluster",
                scale=alt.Scale(domain=clusters, range=palette),
                legend=alt.Legend(orient="bottom", columns=3, symbolSize=40, titleFontSize=11, labelFontSize=10),
            ),
            tooltip=[
                "cell_id:N",
                "cluster:N",
                "ground_truth:N",
                "patient:N",
                "response:N",
            ],
        )
        .properties(
            title=alt.TitleParams(
                text="a   Ground-Truth Cell Clusters",
                subtitle="kNN Graph Leiden Partitioning (Resolution = 0.5)",
                fontSize=13,
                fontWeight="bold",
                subtitleFontSize=10,
                anchor="start",
            ),
            width=width,
            height=height,
        )
    )
    return chart


def build_response_chart(df: pl.DataFrame, width: int = 380, height: int = 320) -> alt.Chart:
    """Panel B: UMAP colored by clinical response phenotype (Responder vs Non-Responder)."""
    pdf = df.to_pandas()
    chart = (
        alt.Chart(pdf)
        .mark_circle(size=28, opacity=0.85, stroke="#334155", strokeWidth=0.25)
        .encode(
            x=alt.X("UMAP_1:Q", title="UMAP 1", axis=alt.Axis(grid=False, tickCount=5)),
            y=alt.Y("UMAP_2:Q", title="UMAP 2", axis=alt.Axis(grid=False, tickCount=5)),
            color=alt.Color(
                "response:N",
                title="Clinical Response",
                scale=alt.Scale(
                    domain=["Responder", "Non-Responder"],
                    range=[OKABE_BLUE, OKABE_VERMILION],
                ),
                legend=alt.Legend(orient="bottom", symbolSize=40, titleFontSize=11, labelFontSize=10),
            ),
            tooltip=[
                "cell_id:N",
                "patient:N",
                "response:N",
                "cluster:N",
            ],
        )
        .properties(
            title=alt.TitleParams(
                text="b   Patient Clinical Outcome",
                subtitle="Single-Cell Distribution across Phenotypic Response Arms",
                fontSize=13,
                fontWeight="bold",
                subtitleFontSize=10,
                anchor="start",
            ),
            width=width,
            height=height,
        )
    )
    return chart


def build_milopy_logfc_chart(df: pl.DataFrame, width: int = 380, height: int = 320) -> alt.Chart:
    """Panel C: UMAP colored by continuous projected Milo DA log2FC (Diverging)."""
    pdf = df.to_pandas()
    max_abs = float(np.percentile(np.abs(pdf["cell_logFC"]), 98))
    max_val = max(1.5, round(max_abs, 1))

    chart = (
        alt.Chart(pdf)
        .mark_circle(size=28, opacity=0.85, stroke="#334155", strokeWidth=0.25)
        .encode(
            x=alt.X("UMAP_1:Q", title="UMAP 1", axis=alt.Axis(grid=False, tickCount=5)),
            y=alt.Y("UMAP_2:Q", title="UMAP 2", axis=alt.Axis(grid=False, tickCount=5)),
            color=alt.Color(
                "cell_logFC:Q",
                title="Milo log2FC (R vs NR)",
                scale=alt.Scale(
                    domain=[-max_val, 0, max_val],
                    range=[OKABE_VERMILION, COLOR_LIGHT_GREY, OKABE_BLUE],
                    clamp=True,
                ),
                legend=alt.Legend(
                    orient="bottom",
                    titleFontSize=11,
                    labelFontSize=10,
                    gradientLength=160,
                ),
            ),
            tooltip=[
                "cell_id:N",
                alt.Tooltip("cell_logFC:Q", format=".2f", title="log2FC"),
                alt.Tooltip("spatial_fdr:Q", format=".2e", title="SpatialFDR"),
                "status:N",
            ],
        )
        .properties(
            title=alt.TitleParams(
                text="c   Milopy Differential Abundance",
                subtitle="Projected Neighborhood log2 Fold-Change (~response)",
                fontSize=13,
                fontWeight="bold",
                subtitleFontSize=10,
                anchor="start",
            ),
            width=width,
            height=height,
        )
    )
    return chart


def build_significance_chart(df: pl.DataFrame, width: int = 380, height: int = 320) -> alt.Chart:
    """Panel D: UMAP colored by categorical DA significance."""
    pdf = df.to_pandas()
    domain = [
        "DA+ (Responder Enriched)",
        "DA- (Non-Responder Enriched)",
        "Not Significant",
    ]
    range_colors = [OKABE_BLUE, OKABE_VERMILION, COLOR_NEUTRAL_GREY]

    chart = (
        alt.Chart(pdf)
        .mark_circle(size=28, opacity=0.85, stroke="#334155", strokeWidth=0.25)
        .encode(
            x=alt.X("UMAP_1:Q", title="UMAP 1", axis=alt.Axis(grid=False, tickCount=5)),
            y=alt.Y("UMAP_2:Q", title="UMAP 2", axis=alt.Axis(grid=False, tickCount=5)),
            color=alt.Color(
                "status:N",
                title="DA Status (SpatialFDR < 0.1)",
                scale=alt.Scale(domain=domain, range=range_colors),
                legend=alt.Legend(orient="bottom", columns=1, symbolSize=40, titleFontSize=11, labelFontSize=10),
            ),
            tooltip=[
                "cell_id:N",
                "status:N",
                alt.Tooltip("cell_logFC:Q", format=".2f", title="log2FC"),
                alt.Tooltip("spatial_fdr:Q", format=".2e", title="SpatialFDR"),
            ],
        )
        .properties(
            title=alt.TitleParams(
                text="d   Statistically Significant Neighborhoods",
                subtitle="SpatialFDR < 0.10 Threshold Envelope",
                fontSize=13,
                fontWeight="bold",
                subtitleFontSize=10,
                anchor="start",
            ),
            width=width,
            height=height,
        )
    )
    return chart


def build_composite_figure(
    chart_a: alt.Chart,
    chart_b: alt.Chart,
    chart_c: alt.Chart,
    chart_d: alt.Chart,
) -> alt.Chart:
    """Builds unified 2x2 publication figure adhering to Nature Methods layout standards."""
    top_row = alt.hconcat(chart_a, chart_b, spacing=35)
    bottom_row = alt.hconcat(chart_c, chart_d, spacing=35)
    composite = (
        alt.vconcat(top_row, bottom_row, spacing=35)
        .properties(
            title=alt.TitleParams(
                text="Differential Abundance Single-Cell Resolution Architecture",
                subtitle="Milopy Graph Neighborhood Testing on Benchmark scRNA-seq Data",
                fontSize=16,
                fontWeight="bold",
                subtitleFontSize=12,
                anchor="start",
                offset=20,
            )
        )
        .configure_view(stroke=COLOR_HAIRLINE, strokeWidth=0.75)
    )
    return composite


def export_svg(chart: alt.Chart, out_path: Path) -> Result[Path, str]:
    """Purely converts Altair chart to SVG vector format and saves to out_path."""
    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        chart_dict = chart.to_dict()
        svg_content = vlc.vegalite_to_svg(chart_dict)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(svg_content)
        return Success(out_path)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"SVG export failed for {out_path}: {exc}")


# =============================================================================
# 3. Main Pipeline Coordination
# =============================================================================

def run_pipeline(config: UmapVisualizerConfig) -> Result[tuple[pl.DataFrame, list[Path]], str]:
    """Executes the complete pipeline: simulation/load, DA, UMAP, and SVG generation."""
    print("=" * 70)
    print("STARTING SYNTHETIC SINGLE-CELL MILOPY UMAP VISUALIZATION PIPELINE")
    print("=" * 70)

    # Step 1: Ingest or Simulate Data
    if config.input_parquet and config.input_parquet.exists():
        print(f"[1/4] Loading pre-computed cell-level metrics from {config.input_parquet}...")
        try:
            cell_table = pl.read_parquet(config.input_parquet)
        except Exception as exc:  # noqa: BLE001
            return Failure(f"Failed to read input parquet {config.input_parquet}: {exc}")
    else:
        print("[1/4] Simulating synthetic single-cell expression with marker programs...")
        rng = np.random.default_rng(config.seed)
        adata = generate_synthetic_scrna(config, rng)
        print(f"      Generated {adata.n_obs:,} cells x {adata.n_vars:,} genes across {config.n_clusters} clusters.")

        print("[2/4] Embedding cells (PCA, kNN Graph, Leiden, UMAP)...")
        match cluster_and_project_umap(adata, config=config, rng=rng):
            case Failure(err):
                return Failure(err)
            case Success(clustered_adata):
                adata = clustered_adata

        print(f"[3/4] Assigning cells to {config.n_patients} patients conditioned on cluster response...")
        match assign_patients(adata, config.n_patients, rng):
            case Failure(err):
                return Failure(err)
            case Success(assigned_adata):
                adata = assigned_adata

        print("[4/4] Executing Milopy differential abundance testing (~response)...")
        match run_milopy_and_project_cells(adata, prop=config.milo_prop, fdr_threshold=config.fdr_threshold):
            case Failure(err):
                return Failure(err)
            case Success(table):
                cell_table = table

        # Save Parquet
        config.out_dir.mkdir(parents=True, exist_ok=True)
        parquet_path = config.out_dir / "synthetic_cells_umap.parquet"
        cell_table.write_parquet(parquet_path)
        print(f"      Saved cell-level UMAP coordinates to {parquet_path}")

    # Generate Altair Charts
    print("\n[PLOTTING] Rendering Nature Methods minimalist wireframe UMAP SVGs...")
    chart_a = build_cluster_chart(cell_table)
    chart_b = build_response_chart(cell_table)
    chart_c = build_milopy_logfc_chart(cell_table)
    chart_d = build_significance_chart(cell_table)
    composite = build_composite_figure(chart_a, chart_b, chart_c, chart_d)

    # Export SVGs
    config.results_dir.mkdir(parents=True, exist_ok=True)
    svg_paths: list[Path] = []

    exports = [
        (chart_a, config.results_dir / "synthetic_umap_clusters.svg"),
        (chart_b, config.results_dir / "synthetic_umap_response.svg"),
        (chart_c, config.results_dir / "synthetic_umap_milopy_logfc.svg"),
        (chart_d, config.results_dir / "synthetic_umap_significance.svg"),
        (composite, config.results_dir / "synthetic_milo_umaps_composite.svg"),
    ]

    for ch, path in exports:
        match export_svg(ch, path):
            case Failure(err):
                print(f"Warning: Failed to export {path.name}: {err}", file=sys.stderr)
            case Success(p):
                print(f"      [EXPORTED] {p}")
                svg_paths.append(p)

    print("=" * 70)
    print(f"PIPELINE COMPLETED SUCCESSFULLY: {len(svg_paths)} SVG figures generated.")
    print("=" * 70)
    return Success((cell_table, svg_paths))


def main() -> None:
    """CLI entrypoint."""
    config = parse_args()
    match run_pipeline(config):
        case Success(_):
            sys.exit(0)
        case Failure(err):
            print(f"[FATAL ERROR] {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
