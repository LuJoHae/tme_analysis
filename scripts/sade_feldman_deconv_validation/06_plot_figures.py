#!/usr/bin/env python3
"""
Step 6: Pure Altair Plotting Engine for Sade-Feldman Deconvolution & Validation Pipeline.
Generates publication-quality resolution-independent vector SVGs (and PNGs) for all analytical steps.
Outputs to results/sade_feldman_deconv_validation/.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Final

import altair as alt  # type: ignore
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success

# Disable row limit for single-cell UMAP coordinates
alt.data_transformers.disable_max_rows()


class PlotConfig(BaseModel):
    model_config = ConfigDict(frozen=True)
    data_dir: Path
    results_dir: Path
    stratum: str = "Melanoma"


def plot_step1_reference_heatmap(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot marker gene expression heatmap across reference cell states."""
    tidy_path = data_dir / "reference_phi_tidy.parquet"
    markers_path = data_dir / "reference_marker_genes.parquet"

    if not tidy_path.exists() or not markers_path.exists():
        return Failure(f"Step 1 files missing in {data_dir}")

    try:
        df_markers = pl.read_parquet(markers_path)
        top_genes = df_markers.filter(pl.col("rank") <= 3)["gene"].unique().to_list()

        df_tidy = pl.read_parquet(tidy_path)
        df_sub = df_tidy.filter(pl.col("gene").is_in(top_genes))

        # Standardize gene expression across clusters for visual contrast
        df_plot = (
            df_sub.with_columns(
                (
                    (pl.col("linear_mean") - pl.col("linear_mean").mean().over("gene"))
                    / (pl.col("linear_mean").std().over("gene") + 1e-6)
                ).alias("z_score")
            )
            .to_pandas()
        )

        chart = (
            alt.Chart(df_plot)
            .mark_rect()
            .encode(
                x=alt.X("gene:N", title="Marker Gene", sort=top_genes),
                y=alt.Y("cluster:N", title="Cell State", sort="ascending"),
                color=alt.Color(
                    "z_score:Q",
                    title="Relative Expression (Z-score)",
                    scale=alt.Scale(scheme="viridis"),
                ),
                tooltip=["cluster", "gene", "linear_mean", "z_score"],
            )
            .properties(
                title="Sade-Feldman Deconvolution Reference: Top Marker Genes",
                width=650,
                height=350,
            )
        )

        out_file = results_dir / "step01_reference_marker_heatmap.svg"
        chart.save(str(out_file))
        return Success(out_file)
    except Exception as exc:
        return Failure(f"Failed to plot step 1 heatmap: {exc}")


def plot_step2_fractions_distribution(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot cell state fraction distributions across cohorts."""
    fracs_path = data_dir / "deconv_fractions.parquet"
    if not fracs_path.exists():
        return Failure(f"Fractions file missing: {fracs_path}")

    try:
        df_fracs = pl.read_parquet(fracs_path)
        meta_cols = ["sample_id", "cohort", "cancer_type"]
        cell_states = [c for c in df_fracs.columns if c not in meta_cols]

        # Melt to long format for Altair
        df_long = (
            df_fracs.unpivot(
                index=meta_cols,
                on=cell_states,
                variable_name="cell_state",
                value_name="fraction",
            )
            .to_pandas()
        )

        chart = (
            alt.Chart(df_long)
            .mark_boxplot(extent="min-max", size=15)
            .encode(
                x=alt.X("cell_state:N", title="Cell State", axis=alt.Axis(labelAngle=-45)),
                y=alt.Y("fraction:Q", title="Inferred Cell Fraction", scale=alt.Scale(zero=True)),
                color=alt.Color("cohort:N", title="Cohort"),
            )
            .properties(
                title="Deconvoluted Cell State Fractions across iAtlas Cohorts",
                width=750,
                height=350,
            )
        )

        out_file = results_dir / "step02_deconv_fractions_distribution.svg"
        chart.save(str(out_file))
        return Success(out_file)
    except Exception as exc:
        return Failure(f"Failed to plot step 2 fractions: {exc}")


def plot_step3_logistic_regression(
    data_dir: Path,
    results_dir: Path,
    stratum: str,
) -> Result[tuple[Path, Path], str]:
    """Plot volcano and forest plot for logistic regression response association."""
    res_path = data_dir / "logistic_regression_results.parquet"
    if not res_path.exists():
        return Failure(f"Logistic results missing: {res_path}")

    try:
        df_res = pl.read_parquet(res_path).filter(pl.col("stratum") == stratum).to_pandas()
        if len(df_res) == 0:
            return Failure(f"No results found for stratum '{stratum}'")

        # Volcano plot: log(OR) vs -log10(p-value)
        volcano = (
            alt.Chart(df_res)
            .mark_circle(size=90, opacity=0.85)
            .encode(
                x=alt.X("beta:Q", title="Log Odds Ratio (Effect Size β)"),
                y=alt.Y("log10_pval:Q", title="-log10(P-value)"),
                color=alt.Color(
                    "significant_fdr01:N",
                    title="FDR < 0.1",
                    scale=alt.Scale(domain=[True, False], range=["#e41a1c", "#377eb8"]),
                ),
                tooltip=["cell_state", "beta", "or", "p_value", "fdr", "auc"],
            )
            .properties(
                title=f"Immunotherapy Response Association: {stratum} Cohorts",
                width=450,
                height=380,
            )
        )
        # Add horizontal significance threshold line (p = 0.05)
        hline = (
            alt.Chart(pd.DataFrame({"y": [-np.log10(0.05)]}))
            .mark_rule(strokeDash=[4, 4], color="gray")
            .encode(y="y:Q")
        )
        volcano_final = volcano + hline
        out_volcano = results_dir / "step03_logistic_regression_volcano.svg"
        volcano_final.save(str(out_volcano))

        # Forest plot: Odds Ratios with 95% Confidence Intervals
        points = (
            alt.Chart(df_res)
            .mark_point(filled=True, size=60)
            .encode(
                x=alt.X("or:Q", title="Odds Ratio (95% CI)", scale=alt.Scale(type="log")),
                y=alt.Y("cell_state:N", title="Cell State", sort=alt.EncodingSortField(field="or", order="descending")),
                color=alt.Color("beta:Q", title="Beta", scale=alt.Scale(scheme="redblue", reverse=True, domain=[-3, 3])),
                tooltip=["cell_state", "or", "or_ci_lower", "or_ci_upper", "p_value"],
            )
        )
        error_bars = (
            alt.Chart(df_res)
            .mark_rule()
            .encode(
                x="or_ci_lower:Q",
                x2="or_ci_upper:Q",
                y=alt.Y("cell_state:N", sort=alt.EncodingSortField(field="or", order="descending")),
                color=alt.value("#555555"),
            )
        )
        vline = (
            alt.Chart(pd.DataFrame({"x": [1.0]}))
            .mark_rule(strokeDash=[3, 3], color="black")
            .encode(x="x:Q")
        )
        forest_final = (error_bars + points + vline).properties(
            title=f"Forest Plot: Odds of Immunotherapy Response ({stratum})",
            width=500,
            height=380,
        )
        out_forest = results_dir / "step03_logistic_regression_forest.svg"
        forest_final.save(str(out_forest))

        return Success((out_volcano, out_forest))
    except Exception as exc:
        return Failure(f"Failed to plot step 3 logistic regression: {exc}")


def plot_step4_milopy(data_dir: Path, results_dir: Path) -> Result[tuple[Path, Path], str]:
    """Plot milopy neighborhood volcano and cell-state differential abundance."""
    nhoods_path = data_dir / "milopy_nhoods_results.parquet"
    states_path = data_dir / "milopy_cell_state_da.parquet"

    if not nhoods_path.exists() or not states_path.exists():
        return Failure(f"Milo files missing in {data_dir}")

    try:
        df_nhoods = pl.read_parquet(nhoods_path).with_columns(
            (-pl.col("p_value").log10()).alias("log10_pval"),
            (pl.col("fdr") < 0.1).alias("significant"),
        ).to_pandas()

        # Nhood volcano
        volcano = (
            alt.Chart(df_nhoods)
            .mark_circle(size=40, opacity=0.7)
            .encode(
                x=alt.X("logfc:Q", title="Log2 Fold Change (Responder vs Non-Responder)"),
                y=alt.Y("log10_pval:Q", title="-log10(P-value)"),
                color=alt.Color(
                    "significant:N",
                    title="FDR < 0.1",
                    scale=alt.Scale(domain=[True, False], range=["#e41a1c", "#9ecae1"]),
                ),
                tooltip=["nhood_id", "logfc", "p_value", "fdr"],
            )
            .properties(
                title="Milopy Single-Cell Differential Abundance (GSE120575)",
                width=450,
                height=380,
            )
        )
        out_volcano = results_dir / "step04_milopy_nhood_volcano.svg"
        volcano.save(str(out_volcano))

        # Cell state DA bar chart with Wilcoxon significance
        df_states = pl.read_parquet(states_path).to_pandas()
        bar = (
            alt.Chart(df_states)
            .mark_bar()
            .encode(
                x=alt.X("milo_mean_logfc:Q", title="Mean Milo Log2FC"),
                y=alt.Y("cell_state:N", title="Cell State", sort=alt.EncodingSortField(field="milo_mean_logfc", order="descending")),
                color=alt.Color(
                    "milo_mean_logfc:Q",
                    title="Milo Log2FC",
                    scale=alt.Scale(scheme="redblue", domain=[-1.5, 1.5], reverse=True),
                ),
                tooltip=["cell_state", "n_cells", "milo_mean_logfc", "milo_median_logfc", "milo_wilcoxon_pval"],
            )
            .properties(
                title="Single-Cell Differential Abundance by Reference Cell State",
                width=500,
                height=380,
            )
        )
        out_bar = results_dir / "step04_milopy_cell_state_da.svg"
        bar.save(str(out_bar))

        return Success((out_volcano, out_bar))
    except Exception as exc:
        return Failure(f"Failed to plot step 4 milopy: {exc}")


def plot_step5_concordance_scatter(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot concordance scatter plot between bulk deconv Beta and Milo logFC."""
    conc_path = data_dir / "concordance_metrics.parquet"
    summary_path = data_dir / "concordance_summary.parquet"

    if not conc_path.exists():
        return Failure(f"Concordance file missing: {conc_path}")

    try:
        df_conc = pl.read_parquet(conc_path).to_pandas()
        rho_val, p_val = 0.0, 1.0
        if summary_path.exists():
            df_sum = pl.read_parquet(summary_path)
            rho_val = float(df_sum["spearman_rho"][0])
            p_val = float(df_sum["spearman_pvalue"][0])

        # Quadrant lines
        hline = alt.Chart(pd.DataFrame({"y": [0.0]})).mark_rule(strokeDash=[3, 3], color="gray").encode(y="y:Q")
        vline = alt.Chart(pd.DataFrame({"x": [0.0]})).mark_rule(strokeDash=[3, 3], color="gray").encode(x="x:Q")

        scatter = (
            alt.Chart(df_conc)
            .mark_circle(size=120, opacity=0.9)
            .encode(
                x=alt.X("beta:Q", title="Bulk Deconvolution Logistic Regression Beta (Log Odds Ratio)"),
                y=alt.Y("milo_mean_logfc:Q", title="Single-Cell Milo Mean Log2FC (~Response)"),
                color=alt.Color(
                    "quadrant:N",
                    title="Concordance Category",
                    scale=alt.Scale(
                        domain=[
                            "Concordant Responder",
                            "Concordant Non-Responder",
                            "Discordant (Bulk+, SC-)",
                            "Discordant (Bulk-, SC+)",
                        ],
                        range=["#e41a1c", "#377eb8", "#ff7f00", "#984ea3"],
                    ),
                ),
                tooltip=["cell_state", "beta", "or", "milo_mean_logfc", "quadrant"],
            )
        )

        labels = (
            alt.Chart(df_conc)
            .mark_text(align="left", baseline="middle", dx=7, fontSize=10)
            .encode(
                x="beta:Q",
                y="milo_mean_logfc:Q",
                text="cell_state:N",
                color=alt.value("#333333"),
            )
        )

        trend = (
            alt.Chart(df_conc)
            .transform_regression("beta", "milo_mean_logfc")
            .mark_line(color="black", strokeDash=[4, 4])
            .encode(x="beta:Q", y="milo_mean_logfc:Q")
        )

        chart = (hline + vline + trend + scatter + labels).properties(
            title=f"Cross-Modality Concordance: Bulk Deconv Beta vs. Single-Cell Milo DA (Spearman rho = {rho_val:.2f}, p = {p_val:.3e})",
            width=580,
            height=460,
        )

        out_file = results_dir / "step05_concordance_scatter.svg"
        chart.save(str(out_file))
        return Success(out_file)
    except Exception as exc:
        return Failure(f"Failed to plot step 5 concordance scatter: {exc}")


def plot_step6_dual_umap(data_dir: Path, results_dir: Path) -> Result[Path, str]:
    """Plot multi-panel UMAP: Reference Cell States, Single-Cell Milo DA, and Bulk Deconvolution Beta."""
    cells_path = data_dir / "milopy_cell_level_scores.parquet"
    conc_path = data_dir / "concordance_metrics.parquet"

    if not cells_path.exists():
        return Failure(f"Cell level scores missing: {cells_path}")

    try:
        df_cells = pl.read_parquet(cells_path)
        # Map deconvolution beta to cells
        if conc_path.exists():
            df_conc = pl.read_parquet(conc_path).select(["cell_state", "beta"])
            df_cells = df_cells.join(df_conc, on="cell_state", how="left")
        else:
            df_cells = df_cells.with_columns(pl.lit(0.0).alias("beta"))

        # Convert to pandas
        df_plot = df_cells.fill_null(0.0).to_pandas()

        # Panel A: Reference Cell States
        panel_a = (
            alt.Chart(df_plot)
            .mark_circle(size=12, opacity=0.8)
            .encode(
                x=alt.X("umap_1:Q", title="UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
                y=alt.Y("umap_2:Q", title="UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
                color=alt.Color("cell_state:N", title="Reference Cell State", scale=alt.Scale(scheme="tableau20")),
                tooltip=["cell_state", "sample_id", "treatment_status"],
            )
            .properties(title="A. Reference Cell States (Sade-Feldman)", width=320, height=320)
        )

        # Panel B: Milopy Single-Cell DA (logFC)
        panel_b = (
            alt.Chart(df_plot)
            .mark_circle(size=12, opacity=0.85)
            .encode(
                x=alt.X("umap_1:Q", title="UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
                y=alt.Y("umap_2:Q", title="UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
                color=alt.Color(
                    "milo_logfc:Q",
                    title="Milo Log2FC (~Response)",
                    scale=alt.Scale(scheme="redblue", reverse=True, domain=[-2.0, 2.0]),
                ),
                tooltip=["cell_state", "milo_logfc"],
            )
            .properties(title="B. Single-Cell Milo DA (~Response)", width=320, height=320)
        )

        # Panel C: Bulk Deconvolution Logistic Regression Beta
        panel_c = (
            alt.Chart(df_plot)
            .mark_circle(size=12, opacity=0.85)
            .encode(
                x=alt.X("umap_1:Q", title="UMAP 1", axis=alt.Axis(labels=False, ticks=False)),
                y=alt.Y("umap_2:Q", title="UMAP 2", axis=alt.Axis(labels=False, ticks=False)),
                color=alt.Color(
                    "beta:Q",
                    title="Bulk Deconv Beta",
                    scale=alt.Scale(scheme="redblue", reverse=True, domain=[-2.0, 2.0]),
                ),
                tooltip=["cell_state", "beta"],
            )
            .properties(title="C. Bulk Deconv Beta (Mapped to Manifold)", width=320, height=320)
        )

        # Combine into side-by-side composite
        composite = alt.hconcat(panel_a, panel_b, panel_c).properties(
            title="Comparison of Single-Cell DA vs. Bulk Deconvolution Response Predictors on Sade-Feldman Manifold"
        )

        out_svg = results_dir / "step06_dual_umap_validation.svg"
        composite.save(str(out_svg))

        # Also save PNG for immediate visual inspection
        try:
            out_png = results_dir / "step06_dual_umap_validation.png"
            composite.save(str(out_png))
        except Exception:
            pass

        return Success(out_svg)
    except Exception as exc:
        return Failure(f"Failed to plot step 6 dual UMAP: {exc}")


def run_all_plots(config: PlotConfig) -> Result[list[Path], str]:
    """Generate all figures across all steps."""
    config.results_dir.mkdir(parents=True, exist_ok=True)
    generated: list[Path] = []

    print("Generating Step 1 Reference Marker Heatmap...")
    match plot_step1_reference_heatmap(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 2 Deconvolution Fractions Distribution...")
    match plot_step2_fractions_distribution(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 3 Logistic Regression Volcano & Forest Plots...")
    match plot_step3_logistic_regression(config.data_dir, config.results_dir, config.stratum):
        case Success((v_path, f_path)):
            generated.extend([v_path, f_path])
            print(f"  Saved: {v_path.name}, {f_path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 4 Milopy Differential Abundance Plots...")
    match plot_step4_milopy(config.data_dir, config.results_dir):
        case Success((v_path, b_path)):
            generated.extend([v_path, b_path])
            print(f"  Saved: {v_path.name}, {b_path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 5 Cross-Modality Concordance Scatter Plot...")
    match plot_step5_concordance_scatter(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    print("Generating Step 6 Multi-Panel UMAP Manifold Visualization...")
    match plot_step6_dual_umap(config.data_dir, config.results_dir):
        case Success(path):
            generated.append(path)
            print(f"  Saved: {path.name}")
        case Failure(err):
            print(f"  Warning: {err}")

    if not generated:
        return Failure("No figures could be generated.")

    return Success(generated)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Step 6: Pure Altair SVG Plotting Engine for Sade-Feldman Pipeline."
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="output/sade_feldman_deconv_validation",
        help="Directory containing output parquet files",
    )
    parser.add_argument(
        "--results-dir",
        type=str,
        default="results/sade_feldman_deconv_validation",
        help="Directory to save vector SVG figures",
    )
    parser.add_argument(
        "--stratum",
        type=str,
        default="Melanoma",
        help="Stratum to highlight in logistic regression plots",
    )
    args = parser.parse_args()

    config = PlotConfig(
        data_dir=Path(args.data_dir),
        results_dir=Path(args.results_dir),
        stratum=args.stratum,
    )

    match run_all_plots(config):
        case Success(paths):
            print(f"\nStep 6 completed successfully. Generated {len(paths)} SVG figures in: {config.results_dir}")
            sys.exit(0)
        case Failure(err):
            print(f"Step 6 failed with error: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
