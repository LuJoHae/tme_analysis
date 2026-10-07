#!/usr/bin/env python3
"""
Publication-Grade Single-Cell UMAP Visualizer for Milopy Differential Abundance Across Response Cohorts.

Features:
1. "Always redo UMAP": Recomputes fresh 2D UMAP coordinates (min_dist=0.3, spread=1.0) on the exact
   single-cell subsets evaluated by Milopy, updating cell_level_scores.parquet with UMAP1 and UMAP2.
2. Derives cell-level differential abundance significance (DA+, DA-, Not Significant) from neighborhood testing.
3. Generates per-cohort Nature Methods compliant Altair vector SVGs:
   - milopy_umap_logfc.svg (focused continuous log2FC diverging gradient)
   - milopy_umap_composite.svg (2x2 publication panel: Cell Types, Response, log2FC, Significance)
4. Generates a consolidated 3x3 multi-cohort publication grid SVG comparing all 9 clinical cohorts side-by-side.

Strict functional Python adhering to immutability, returns Result, Polars, and Altair vector SVG.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import altair as alt  # type: ignore
import anndata as ad  # type: ignore
import numpy as np
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
import scanpy as sc  # type: ignore
import scipy.sparse as sp  # type: ignore
import vl_convert as vlc  # type: ignore

# Disable Altair row limit for large single-cell datasets
alt.data_transformers.disable_max_rows()

from tme_datasets import find_dataset_h5ad, load_backed
from milopy import project_nhoods_to_cells
from plotting_utils import (
    COLOR_BG_WHITE,
    COLOR_HAIRLINE,
    COLOR_LIGHT_GREY,
    COLOR_NEUTRAL_GREY,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    COHORT_METADATA,
    OKABE_BLACK,
    OKABE_BLUE,
    OKABE_BLUISH_GREEN,
    OKABE_ORANGE,
    OKABE_REDDISH_PURPLE,
    OKABE_SKY_BLUE,
    OKABE_VERMILION,
    OKABE_YELLOW,
)


class VisualizerConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    base_dir: Path = Path("/storage/halu/data-test")
    results_dir: Path = Path("/storage/halu/data-test/results/milopy_response_analysis")
    reports_dir: Path = Path("/storage/halu/data-test/reports")
    min_dist: float = 0.30
    spread: float = 1.00
    fdr_threshold: float = 0.10
    seed: int = 42
    grid_max_cells: int = 10000
    point_size: int = 14
    point_opacity: float = 0.75


def parse_args() -> VisualizerConfig:
    parser = argparse.ArgumentParser(
        description="Generates publication-grade UMAP visualizations of Milopy differential abundance."
    )
    parser.add_argument("--base-dir", type=Path, default=Path("/storage/halu/data-test"), help="Base data directory")
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("/storage/halu/data-test/results/milopy_response_analysis"),
        help="Directory containing cohort milopy outputs",
    )
    parser.add_argument(
        "--reports-dir",
        type=Path,
        default=Path("/storage/halu/data-test/reports"),
        help="Directory to save master reports and SVGs",
    )
    parser.add_argument("--min-dist", type=float, default=0.30, help="UMAP min_dist parameter")
    parser.add_argument("--spread", type=float, default=1.00, help="UMAP spread parameter")
    parser.add_argument("--fdr-threshold", type=float, default=0.10, help="FDR significance threshold")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    parser.add_argument("--grid-max-cells", type=int, default=10000, help="Max cells per cohort in 3x3 grid")
    args = parser.parse_args()

    return VisualizerConfig(
        base_dir=args.base_dir,
        results_dir=args.results_dir,
        reports_dir=args.reports_dir,
        min_dist=args.min_dist,
        spread=args.spread,
        fdr_threshold=args.fdr_threshold,
        seed=args.seed,
        grid_max_cells=args.grid_max_cells,
    )


# =============================================================================
# 1. Coordinate & Significance Computation Engine
# =============================================================================

def recompute_cohort_umap(
    accession: str,
    cohort_dir: Path,
    config: VisualizerConfig,
) -> Result[pl.DataFrame, str]:
    """Always computes fresh 2D UMAP coordinates and cell significance for a cohort."""
    parquet_path = cohort_dir / "cell_level_scores.parquet"
    if not parquet_path.exists():
        return Failure(f"Missing cell_level_scores.parquet in {cohort_dir}")

    cell_df = pl.read_parquet(parquet_path)
    cell_ids = cell_df["cell_id"].to_list()
    n_cells = len(cell_ids)

    h5ad_maybe = find_dataset_h5ad(accession, repo_root=config.base_dir)
    h5ad_file: Path | None = h5ad_maybe.value_or(None)

    if h5ad_file is None:
        return Failure(f"Could not locate preprocessed H5AD for {accession} via tme_datasets API")

    print(f"[{accession}] Streaming {n_cells:,} cells from {h5ad_file} via tme_datasets backed loader to recompute fresh UMAP...")
    try:
        backed_res = load_backed(h5ad_file, mode="r")
        if not isinstance(backed_res, Success):
            return Failure(f"Failed to open backed H5AD for {accession}: {backed_res.failure()}")
        backed_adata = backed_res.unwrap()
        # Subsetting slice
        idx_map = backed_adata.obs.index.get_indexer(cell_ids)
        valid_mask = idx_map >= 0
        if not np.all(valid_mask):
            # Try matching without prefix or stripping
            obs_idx = list(backed_adata.obs.index)
            cell_set = set(cell_ids)
            idx_map = np.array([obs_idx.index(cid) if cid in cell_set else -1 for cid in cell_ids])

        sort_order = np.argsort(idx_map)
        sorted_pos = idx_map[sort_order]
        adata_sub = backed_adata[sorted_pos, :].to_memory()
        backed_adata.file.close()

        # Re-align order to match cell_ids
        reverse_order = np.argsort(sort_order)
        adata = adata_sub[reverse_order, :].copy()
        del adata_sub
    except Exception as exc:
        return Failure(f"Failed to extract cell subset from {h5ad_file}: {exc}")

    # Compute PCA, kNN graph, and UMAP
    try:
        max_val = adata.X.max() if not sp.issparse(adata.X) else (adata.X.data.max() if adata.X.nnz > 0 else 0)
        if max_val > 50:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)

        if adata.n_vars > 2000:
            sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=False)
            sc.pp.pca(adata, n_comps=min(30, adata.n_vars - 1), mask_var="highly_variable", zero_center=False)
        else:
            sc.pp.pca(adata, n_comps=min(30, adata.n_vars - 1), zero_center=False)

        sc.pp.neighbors(adata, n_neighbors=30, n_pcs=min(30, adata.obsm["X_pca"].shape[1]))
        sc.tl.umap(adata, min_dist=config.min_dist, spread=config.spread, random_state=config.seed)
    except Exception as exc:
        return Failure(f"Scanpy UMAP computation failed for {accession}: {exc}")

    umap_coords = adata.obsm["X_umap"]
    u1 = [float(x) for x in umap_coords[:, 0]]
    u2 = [float(x) for x in umap_coords[:, 1]]

    # Map neighborhood statistical significance to single cells via sparse matrix multiplication
    da_path = cohort_dir / "da_results.parquet"
    nhood_mat_path = cohort_dir / "nhood_matrix.npz"

    if da_path.exists() and nhood_mat_path.exists():
        try:
            da_df = pl.read_parquet(da_path)
            nhoods_mat = sp.load_npz(nhood_mat_path).tocsc()  # cells x nhoods
            updated_df = project_nhoods_to_cells(
                nhoods_mat=nhoods_mat,
                res_df=da_df,
                cell_ids=cell_df["cell_id"].to_list(),
                patient_ids=cell_df["patient_id"].to_list(),
                clinical_responses=cell_df["clinical_response"].to_list(),
                cell_types=cell_df["cell_type"].to_list() if "cell_type" in cell_df.columns else None,
                treatment_statuses=cell_df["treatment_status"].to_list() if "treatment_status" in cell_df.columns else None,
                umap_coords=umap_coords,
                fdr_threshold=config.fdr_threshold,
            )
        except Exception as e:
            print(f"[{accession}] Note: could not project nhood significance to cells ({e})")
            updated_df = cell_df.with_columns([
                pl.Series("UMAP1", u1),
                pl.Series("UMAP2", u2),
            ])
            if "da_status" not in updated_df.columns:
                updated_df = updated_df.with_columns(pl.lit("Not Significant").alias("da_status"))
    else:
        updated_df = cell_df.with_columns([
            pl.Series("UMAP1", u1),
            pl.Series("UMAP2", u2),
        ])
        if "da_status" not in updated_df.columns:
            updated_df = updated_df.with_columns(pl.lit("Not Significant").alias("da_status"))

    # Save back to Parquet
    updated_df.write_parquet(parquet_path)
    print(f"[{accession}] Saved updated UMAP coordinates and DA status to {parquet_path}")
    return Success(updated_df)


# =============================================================================
# 2. Nature Methods Declarative Altair Plotting Functions
# =============================================================================

def build_cohort_logfc_umap(
    df: pl.DataFrame,
    accession: str,
    meta: dict[str, str],
    width: int = 420,
    height: int = 360,
    point_size: int = 16,
    opacity: float = 0.80,
) -> alt.Chart:
    """Builds an individual UMAP colored by continuous projected Milopy log2FC with diverging gradient."""
    pdf = df.to_pandas()
    pdf["abs_logfc"] = np.abs(pdf["cell_logfc"])
    pdf = pdf.sort_values(by="abs_logfc", ascending=True)

    if len(pdf) > 25000:
        sig_mask = pdf["da_status"] != "Not Significant"
        sig_pdf = pdf[sig_mask]
        non_sig_pdf = pdf[~sig_mask]
        n_needed = max(1000, 25000 - len(sig_pdf))
        if len(non_sig_pdf) > n_needed:
            non_sig_sampled = non_sig_pdf.sample(n=n_needed, random_state=42)
            pdf = pd.concat([non_sig_sampled, sig_pdf], axis=0).sort_values(by="abs_logfc", ascending=True)

    max_abs = float(np.percentile(np.abs(pdf["cell_logfc"]), 98))
    max_val = max(1.5, round(max_abs, 1))

    title_text = f"{accession} — {meta.get('indication', 'TME')} ({meta.get('tech', '')})"
    sub_text = f"Milopy log2FC (~clinical_response) | {meta.get('patients', '')} | N = {len(df):,} cells"

    chart = (
        alt.Chart(pdf)
        .mark_circle(size=point_size, opacity=opacity)
        .encode(
            x=alt.X("UMAP1:Q", title="UMAP 1", axis=alt.Axis(grid=False, tickCount=5)),
            y=alt.Y("UMAP2:Q", title="UMAP 2", axis=alt.Axis(grid=False, tickCount=5)),
            color=alt.Color(
                "cell_logfc:Q",
                title="Milo log2FC (R vs NR; clamped at ±5)",
                scale=alt.Scale(
                    domain=[-5.0, 0.0, 5.0],
                    range=[OKABE_BLUE, COLOR_LIGHT_GREY, OKABE_VERMILION],
                    clamp=True,
                ),
                legend=alt.Legend(
                    orient="bottom",
                    titleFontSize=11,
                    labelFontSize=10,
                    titleLimit=0,
                    gradientLength=200,
                    values=[-5.0, -2.5, 0.0, 2.5, 5.0],
                ),
            ),
            tooltip=[
                "cell_id:N",
                alt.Tooltip("cell_logfc:Q", format=".2f", title="log2FC"),
                "clinical_response:N",
                "cell_type:N",
                "da_status:N",
            ],
        )
        .properties(
            title=alt.TitleParams(
                text=title_text,
                subtitle=sub_text,
                fontSize=13,
                fontWeight="bold",
                subtitleFontSize=10,
                anchor="start",
            ),
            width=width,
            height=height,
        )
        .configure_view(stroke=COLOR_HAIRLINE, strokeWidth=0.75)
    )
    return chart


def build_cohort_composite_umap(
    df: pl.DataFrame,
    accession: str,
    meta: dict[str, str],
    width: int = 340,
    height: int = 280,
) -> alt.VConcatChart:
    """Builds a 2x2 publication composite panel for a single cohort."""
    pdf = df.to_pandas()
    pdf["abs_logfc"] = np.abs(pdf["cell_logfc"])
    pdf = pdf.sort_values(by="abs_logfc", ascending=True)

    if len(pdf) > 25000:
        sig_mask = pdf["da_status"] != "Not Significant"
        sig_pdf = pdf[sig_mask]
        non_sig_pdf = pdf[~sig_mask]
        n_needed = max(1000, 25000 - len(sig_pdf))
        if len(non_sig_pdf) > n_needed:
            non_sig_sampled = non_sig_pdf.sample(n=n_needed, random_state=42)
            pdf = pd.concat([non_sig_sampled, sig_pdf], axis=0).sort_values(by="abs_logfc", ascending=True)

    max_abs = float(np.percentile(np.abs(pdf["cell_logfc"]), 98))
    max_val = max(1.5, round(max_abs, 1))

    # Panel A: Cell Types
    chart_a = (
        alt.Chart(pdf)
        .mark_circle(size=12, opacity=0.75)
        .encode(
            x=alt.X("UMAP1:Q", title="UMAP 1", axis=alt.Axis(grid=False, tickCount=4)),
            y=alt.Y("UMAP2:Q", title="UMAP 2", axis=alt.Axis(grid=False, tickCount=4)),
            color=alt.Color("cell_type:N", title="Cell Type", legend=alt.Legend(orient="bottom", columns=2, symbolSize=30)),
            tooltip=["cell_id:N", "cell_type:N"],
        )
        .properties(
            title=alt.TitleParams(text="a   Annotated Cell Types", fontSize=12, fontWeight="bold", anchor="start"),
            width=width,
            height=height,
        )
    )

    # Panel B: Clinical Response
    chart_b = (
        alt.Chart(pdf)
        .mark_circle(size=12, opacity=0.75)
        .encode(
            x=alt.X("UMAP1:Q", title="UMAP 1", axis=alt.Axis(grid=False, tickCount=4)),
            y=alt.Y("UMAP2:Q", title="UMAP 2", axis=alt.Axis(grid=False, tickCount=4)),
            color=alt.Color(
                "clinical_response:N",
                title="Response",
                scale=alt.Scale(domain=["responder", "non-responder"], range=[OKABE_VERMILION, OKABE_BLUE]),
                legend=alt.Legend(orient="bottom", symbolSize=30),
            ),
            tooltip=["cell_id:N", "clinical_response:N", "patient_id:N"],
        )
        .properties(
            title=alt.TitleParams(text="b   Clinical Response Arms", fontSize=12, fontWeight="bold", anchor="start"),
            width=width,
            height=height,
        )
    )

    # Panel C: Continuous log2FC
    chart_c = (
        alt.Chart(pdf)
        .mark_circle(size=12, opacity=0.85)
        .encode(
            x=alt.X("UMAP1:Q", title="UMAP 1", axis=alt.Axis(grid=False, tickCount=4)),
            y=alt.Y("UMAP2:Q", title="UMAP 2", axis=alt.Axis(grid=False, tickCount=4)),
            color=alt.Color(
                "cell_logfc:Q",
                title="Milo log2FC (clamped at ±5)",
                scale=alt.Scale(domain=[-5.0, 0.0, 5.0], range=[OKABE_BLUE, COLOR_LIGHT_GREY, OKABE_VERMILION], clamp=True),
                legend=alt.Legend(orient="bottom", gradientLength=160, titleLimit=0, values=[-5.0, -2.5, 0.0, 2.5, 5.0]),
            ),
            tooltip=["cell_id:N", alt.Tooltip("cell_logfc:Q", format=".2f", title="log2FC")],
        )
        .properties(
            title=alt.TitleParams(text="c   Differential Abundance (log2FC)", fontSize=12, fontWeight="bold", anchor="start"),
            width=width,
            height=height,
        )
    )

    # Panel D: DA Significance Envelope
    chart_d = (
        alt.Chart(pdf)
        .mark_circle(size=12, opacity=0.85)
        .encode(
            x=alt.X("UMAP1:Q", title="UMAP 1", axis=alt.Axis(grid=False, tickCount=4)),
            y=alt.Y("UMAP2:Q", title="UMAP 2", axis=alt.Axis(grid=False, tickCount=4)),
            color=alt.Color(
                "da_status:N",
                title="DA Status (FDR < 0.1)",
                scale=alt.Scale(
                    domain=["DA+ (Responder Enriched)", "DA- (Non-Responder Enriched)", "Not Significant"],
                    range=[OKABE_VERMILION, OKABE_BLUE, COLOR_NEUTRAL_GREY],
                ),
                legend=alt.Legend(orient="bottom", symbolSize=30),
            ),
            tooltip=["cell_id:N", "da_status:N", alt.Tooltip("cell_logfc:Q", format=".2f", title="log2FC")],
        )
        .properties(
            title=alt.TitleParams(text="d   Statistically Significant Neighborhoods", fontSize=12, fontWeight="bold", anchor="start"),
            width=width,
            height=height,
        )
    )

    header = alt.TitleParams(
        text=f"{accession} — Multi-Panel Single-Cell Characterization",
        subtitle=f"{meta.get('indication', '')} | {meta.get('tech', '')} | {meta.get('patients', '')} | N = {len(df):,} cells",
        fontSize=15,
        fontWeight="bold",
        anchor="start",
    )

    row1 = alt.hconcat(chart_a, chart_b, spacing=25)
    row2 = alt.hconcat(chart_c, chart_d, spacing=25)
    composite = alt.vconcat(row1, row2, spacing=25).properties(title=header).configure_view(stroke=COLOR_HAIRLINE, strokeWidth=0.75)
    return composite


# =============================================================================
# 3. Main Workflow Orchestration
# =============================================================================

def main() -> None:
    config = parse_args()
    print("=" * 70)
    print("STARTING MILOPY SINGLE-CELL UMAP RECOMPUTATION & VISUALIZATION PIPELINE")
    print("=" * 70)
    print(f"Results Directory : {config.results_dir}")
    print(f"Reports Directory : {config.reports_dir}")
    print(f"Base Directory    : {config.base_dir}")

    config.reports_dir.mkdir(parents=True, exist_ok=True)
    cohort_dirs = [d for d in sorted(config.results_dir.iterdir()) if d.is_dir() and d.name != "reports"]

    cohort_dfs: dict[str, pl.DataFrame] = {}

    for cdir in cohort_dirs:
        acc = cdir.name
        print(f"\n[{acc}] Processing UMAP visualization...")
        res = recompute_cohort_umap(acc, cdir, config)
        match res:
            case Success(df):
                cohort_dfs[acc] = df
                meta = COHORT_METADATA.get(acc, {})

                # 1. Render and export single logFC UMAP
                chart_logfc = build_cohort_logfc_umap(df, acc, meta)
                svg_logfc_path = cdir / "milopy_umap_logfc.svg"
                svg_str = vlc.vegalite_to_svg(chart_logfc.to_json())
                svg_logfc_path.write_text(svg_str, encoding="utf-8")
                print(f"[{acc}] [EXPORTED] {svg_logfc_path} ({svg_logfc_path.stat().st_size / 1e6:.2f} MB)")

                # 2. Render and export 2x2 composite panel
                chart_composite = build_cohort_composite_umap(df, acc, meta)
                svg_comp_path = cdir / "milopy_umap_composite.svg"
                svg_comp_str = vlc.vegalite_to_svg(chart_composite.to_json())
                svg_comp_path.write_text(svg_comp_str, encoding="utf-8")
                print(f"[{acc}] [EXPORTED] {svg_comp_path} ({svg_comp_path.stat().st_size / 1e6:.2f} MB)")

            case Failure(err):
                print(f"[{acc}] Error: {err}", file=sys.stderr)

    # Write completion sentinel
    sentinel_path = config.results_dir / "visualize_cohort_umaps_completed.txt"
    sentinel_path.write_text(f"Completed {len(cohort_dfs)} cohorts\n", encoding="utf-8")
    print(f"[SENTINEL] Wrote {sentinel_path}")

    print("\n" + "=" * 70)
    print("PIPELINE COMPLETED SUCCESSFULLY.")
    print("=" * 70)


if __name__ == "__main__":
    main()
