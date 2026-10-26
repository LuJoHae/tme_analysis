#!/usr/bin/env python3
"""Plot Comparative Stability Paths (MB 2010 vs SS-CPSS 2013) for iAtlas Cohorts.

Follows strict functional programming principles:
- Pure functions mapping Polars DataFrames to Altair Chart specifications
- Nature Methods minimalist wireframe standard via plotting_utils
- Okabe-Ito colorblind-safe palette for active biomarkers
- Vector SVG exports via vl-convert-python
- Monadic error handling with Result[T, str]
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Final, Sequence
import altair as alt
import numpy as np
import polars as pl
from returns.result import Failure, Result, Success
import vl_convert as vlc

from plotting_utils import (
    COLOR_BORDER_HAIRLINE,
    COLOR_NEUTRAL_GREY,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    OKABE_PALETTE,
    THRESHOLD_COLOR,
    apply_nature_methods_theme,
    create_nature_axis,
)
from plotting_utils.svg import verify_and_rasterize_svg

# Disable 5,000 row limit for full regularization paths
alt.data_transformers.disable_max_rows()


def build_cohort_paths_chart(
    df_cohort: pl.DataFrame,
    cohort_name: str,
    summary_row: pl.DataFrame | None = None,
) -> Result[alt.HConcatChart, str]:
    """Pure function generating a comparative two-panel stability paths chart (MB vs SS-CPSS)."""
    if df_cohort.is_empty():
        return Failure(f"No data available for cohort '{cohort_name}'")

    # Add negative log lambda for standard left-to-right penalty relaxation
    df_prep = df_cohort.with_columns(
        (-pl.col("log_lambda")).alias("neg_log_lambda")
    )

    cutoff_val = float(df_prep["cutoff"][0]) if "cutoff" in df_prep.columns else 0.75

    stable_features = (
        df_prep.filter(pl.col("selected"))
        .select("feature")
        .unique()
        .to_series()
        .to_list()
    )

    methods = ("MB", "SS-CPSS")
    method_titles = {
        "MB": "Meinshausen & Bühlmann (2010)",
        "SS-CPSS": "Shah & Samworth (2013) CPSS",
    }

    panel_charts: list[alt.LayerChart] = []

    for method_key in methods:
        df_method = df_prep.filter(pl.col("method") == method_key)
        if df_method.is_empty():
            continue

        n_selected = (
            df_method.filter(pl.col("selected"))
            .select("feature")
            .n_unique()
        )
        lam_cutoff_val = 0.0
        if "lambda_cutoff" in df_method.columns and not df_method["lambda_cutoff"].is_null().all():
            first_cut = df_method["lambda_cutoff"][0]
            lam_cutoff_val = float(first_cut) if isinstance(first_cut, (int, float)) else 0.0

        if summary_row is not None and not summary_row.is_empty():
            cut_col = "lambda_cutoff_mb" if method_key == "MB" else "lambda_cutoff_ss"
            if cut_col in summary_row.columns:
                lam_cutoff_val = float(summary_row[cut_col][0])

        cut_str = f" | λ_cut = {lam_cutoff_val:.4g}" if lam_cutoff_val > 0.0 else ""
        subtitle_info = f"Threshold π_thr = {cutoff_val:.2f}{cut_str} | {n_selected} stable features"

        if summary_row is not None and not summary_row.is_empty():
            q_col = "mb_q" if method_key == "MB" else "ss_q"
            act_q_col = "actual_mb_q" if method_key == "MB" else "actual_ss_q"
            if q_col in summary_row.columns:
                q_val = float(summary_row[q_col][0])
                act_str = ""
                if act_q_col in summary_row.columns:
                    act_q = float(summary_row[act_q_col][0])
                    act_str = f" (actual {act_q:.1f})"
                subtitle_info = f"Budget q = {q_val:.1f}{act_str}{cut_str} | π_thr = {cutoff_val:.2f} | {n_selected} selected"

        df_bg = df_method.filter(~pl.col("selected"))
        df_fg = df_method.filter(pl.col("selected"))

        min_elem = df_prep["neg_log_lambda"].min()
        max_elem = df_prep["neg_log_lambda"].max()
        x_min = float(min_elem) if isinstance(min_elem, (int, float)) else 0.0
        x_max = float(max_elem) if isinstance(max_elem, (int, float)) else 4.0

        x_axis = create_nature_axis(title="-log₁₀(λ)  [Penalty Relaxation →]")
        y_axis = create_nature_axis(
            title="Selection Probability  Π̂(λ)",
            values=[0.0, 0.25, 0.5, 0.75, 1.0],
        )

        scale_x = alt.Scale(domain=[x_min, x_max])
        scale_y = alt.Scale(domain=[0.0, 1.05])

        # 1. Background paths
        bg_chart = (
            alt.Chart(df_bg if not df_bg.is_empty() else df_method)
            .mark_line(strokeWidth=0.75, opacity=0.25, color=COLOR_NEUTRAL_GREY)
            .encode(
                x=alt.X("neg_log_lambda:Q", axis=x_axis, scale=scale_x),
                y=alt.Y("selection_probability:Q", axis=y_axis, scale=scale_y),
                detail="feature:N",
            )
        )

        # 2. Foreground paths (stable biomarkers)
        fg_chart = None
        if not df_fg.is_empty():
            fg_chart = (
                alt.Chart(df_fg)
                .mark_line(strokeWidth=2.2, opacity=0.9)
                .encode(
                    x=alt.X("neg_log_lambda:Q", scale=scale_x),
                    y=alt.Y("selection_probability:Q", scale=scale_y),
                    color=alt.Color(
                        "feature:N",
                        scale=alt.Scale(scheme="tableau10"),
                        legend=alt.Legend(
                            title="Stable Biomarkers",
                            titleFontSize=10,
                            labelFontSize=9,
                            symbolStrokeWidth=2,
                            orient="right",
                            columns=1 if len(stable_features) <= 15 else 2,
                        ),
                    ),
                    detail="feature:N",
                    tooltip=[
                        alt.Tooltip("feature:N", title="Gene"),
                        alt.Tooltip("neg_log_lambda:Q", title="-log10(lambda)", format=".2f"),
                        alt.Tooltip("selection_probability:Q", title="Probability", format=".3f"),
                        alt.Tooltip("is_checkpoint:N", title="Checkpoint"),
                    ],
                )
            )

        # 3. Horizontal threshold rule
        rule_threshold = (
            alt.Chart(pl.DataFrame({"y": [cutoff_val]}))
            .mark_rule(
                strokeDash=[5, 3],
                strokeWidth=1.2,
                color=THRESHOLD_COLOR,
                opacity=0.8,
            )
            .encode(y="y:Q")
        )

        # 4. Vertical regularization budget cutoff rule
        budget_rule = None
        budget_text = None
        if lam_cutoff_val > 0.0:
            neg_log_cut = -float(np.log10(lam_cutoff_val))
            budget_rule = (
                alt.Chart(pl.DataFrame({"x": [neg_log_cut]}))
                .mark_rule(
                    strokeDash=[4, 4],
                    strokeWidth=1.5,
                    color="#0284C7",
                    opacity=0.85,
                )
                .encode(x="x:Q")
            )
            budget_text = (
                alt.Chart(pl.DataFrame({"x": [neg_log_cut], "y": [0.03], "text": ["Budget Limit"]}))
                .mark_text(
                    align="right",
                    dx=-5,
                    fontSize=8.5,
                    fontWeight=600,
                    color="#0284C7",
                    font="sans-serif",
                )
                .encode(x="x:Q", y="y:Q", text="text:N")
            )

        layers = [bg_chart]
        if fg_chart is not None:
            layers.append(fg_chart)
        layers.append(rule_threshold)
        if budget_rule is not None:
            layers.append(budget_rule)
        if budget_text is not None:
            layers.append(budget_text)

        panel = (
            alt.layer(*layers)
            .properties(
                title=alt.TitleParams(
                    text=method_titles.get(method_key, method_key),
                    subtitle=subtitle_info,
                    fontSize=13,
                    subtitleFontSize=10,
                    fontWeight="bold",
                    color=COLOR_TEXT_PRIMARY,
                    subtitleColor=COLOR_TEXT_SECONDARY,
                    anchor="start",
                    offset=10,
                ),
                width=360,
                height=260,
            )
        )
        panel_charts.append(panel)

    if not panel_charts:
        return Failure(f"No valid panel charts generated for '{cohort_name}'")

    combined_chart = (
        alt.hconcat(*panel_charts)
        .properties(
            title=alt.TitleParams(
                text=f"{cohort_name}: Finite-Sample Stability Paths",
                subtitle="Selection probability paths along the regularization trajectory",
                fontSize=15,
                subtitleFontSize=11,
                fontWeight="bold",
                color=COLOR_TEXT_PRIMARY,
                subtitleColor=COLOR_TEXT_SECONDARY,
                anchor="start",
                offset=14,
            ),
            spacing=30,
        )
    )

    themed_chart = apply_nature_methods_theme(combined_chart)
    return Success(themed_chart)


def export_chart_to_svg(chart: alt.Chart | alt.HConcatChart | alt.VConcatChart, target_path: Path) -> Result[Path, str]:
    """Pure function serializing an Altair specification to an SVG file via vl-convert."""
    try:
        target_path.parent.mkdir(parents=True, exist_ok=True)
        chart_dict = chart.to_dict()
        svg_content = vlc.vegalite_to_svg(chart_dict)
        return verify_and_rasterize_svg(svg_content, target_path, min_layers=0)
    except Exception as exc:
        return Failure(f"Failed to export SVG to {target_path}: {exc}")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot comparative stability selection paths (MB vs SS-CPSS) for iAtlas datasets."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help="Directory containing stability_paths.parquet and cohort_summary.parquet (default: auto-detected)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory to export SVG vector figures (default: <input-dir>/figures)",
    )
    parser.add_argument(
        "--cohorts",
        nargs="+",
        default=None,
        help="Optional subset of cohorts to plot (default: all available in Parquet)",
    )
    parser.add_argument(
        "--fitter",
        type=str,
        default="lasso",
        help="Base fitter to visualize (default: lasso)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    input_dir: Path
    if args.input_dir is not None:
        input_dir = args.input_dir
    elif Path("output/stability_selection_fitters_benchmark").exists():
        input_dir = Path("output/stability_selection_fitters_benchmark")
    elif Path("output/stability_selection_iatlas_immunotherapy").exists():
        input_dir = Path("output/stability_selection_iatlas_immunotherapy")
    else:
        input_dir = Path("output/stability_selection_iatlas")

    output_dir: Path = args.output_dir if args.output_dir is not None else (input_dir / "figures")

    paths_file = input_dir / "stability_paths.parquet"
    summary_file = (
        input_dir / "cohort_fitter_summary.parquet"
        if (input_dir / "cohort_fitter_summary.parquet").exists()
        else input_dir / "cohort_summary.parquet"
    )

    if not paths_file.exists():
        print(f"[-] Error: Missing {paths_file}. Run run_stability_selection_iatlas.py first.")
        return 1

    print("==================================================================")
    print(f" Generating Publication SVG Stability Paths Figures (Fitter: {args.fitter.upper()})")
    print(f" Input:  {paths_file}")
    print(f" Output: {output_dir}")
    print("==================================================================")

    paths_df = pl.read_parquet(paths_file)
    if "fitter" in paths_df.columns:
        paths_df = paths_df.filter(pl.col("fitter") == args.fitter)

    summary_df = None
    if summary_file.exists():
        summary_df = pl.read_parquet(summary_file)
        if "fitter" in summary_df.columns:
            summary_df = summary_df.filter(pl.col("fitter") == args.fitter)

    available_cohorts = sorted(paths_df["cohort"].unique().to_list())
    target_cohorts = args.cohorts if args.cohorts else available_cohorts

    print(f"[+] Found {len(available_cohorts)} cohorts: {available_cohorts}")
    print(f"[+] Generating figures for: {target_cohorts}")

    successful_count = 0
    for cohort_id in target_cohorts:
        print(f"\n[+] Rendering stability paths for '{cohort_id}'...")
        cohort_df = paths_df.filter(pl.col("cohort") == cohort_id)
        if cohort_df.is_empty():
            print(f"    [-] Cohort '{cohort_id}' not found in paths table. Skipping.")
            continue

        cohort_summary = (
            summary_df.filter(pl.col("cohort") == cohort_id)
            if summary_df is not None
            else None
        )

        match build_cohort_paths_chart(cohort_df, cohort_id, cohort_summary):
            case Failure(err):
                print(f"    [-] Failed to build chart for '{cohort_id}': {err}")
            case Success(chart):
                clean_name = cohort_id.lower().replace("-", "_")
                svg_target = output_dir / f"paths_comparison_{clean_name}.svg"
                match export_chart_to_svg(chart, svg_target):
                    case Failure(err):
                        print(f"    [-] SVG export failed: {err}")
                    case Success(path):
                        print(f"    [x] Exported SVG: {path} ({path.stat().st_size:,} bytes)")
                        successful_count += 1

    print("\n==================================================================")
    print(f" Successfully exported {successful_count} publication SVG figures to:")
    print(f" {output_dir}")
    print("==================================================================")
    return 0 if successful_count > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
