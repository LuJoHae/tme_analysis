"""Publication-grade Altair vector SVG diagnostics for reference sampling HPO."""

from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence
import altair as alt  # type: ignore[import-untyped]
import polars as pl
from returns.result import Failure, Result, Success

# Okabe-Ito colorblind-safe palette
OKABE_ITO_PALETTE = [
    "#0072B2",  # Blue
    "#E69F00",  # Orange
    "#009E73",  # Bluish green
    "#D55E00",  # Vermillion
    "#CC79A7",  # Reddish purple
    "#56B4E9",  # Sky blue
]


def plot_pareto_frontier(
    eval_df: pl.DataFrame,
    pareto_trial_ids: Sequence[int],
    out_file: Path,
) -> Result[Path, str]:
    """Generate publication-standard Altair SVG scatter plot of Pareto frontier."""
    try:
        out_file.parent.mkdir(parents=True, exist_ok=True)
        pareto_set = set(pareto_trial_ids)

        plot_df = eval_df.with_columns([
            pl.col("trial_id").is_in(list(pareto_set)).alias("is_pareto"),
            pl.when(pl.col("is_pruned"))
            .then(pl.lit("Pruned (Early Stop)"))
            .otherwise(pl.col("rung"))
            .alias("status"),
        ])

        # Convert to pandas for Altair
        pdf = plot_df.to_pandas()

        base = (
            alt.Chart(pdf)
            .encode(
                x=alt.X(
                    "collinearity_max:Q",
                    title="Max Reference State Collinearity (r)",
                    scale=alt.Scale(zero=False, padding=10),
                ),
                y=alt.Y(
                    "mean_loco_auc:Q",
                    title="Leave-One-Cohort-Out ROC-AUC",
                    scale=alt.Scale(zero=False, padding=10),
                ),
            )
        )

        all_points = base.mark_circle().encode(
            color=alt.Color(
                "status:N",
                scale=alt.Scale(
                    domain=["rung_0_screening", "rung_1_refinement", "rung_2_full", "Pruned (Early Stop)"],
                    range=["#56B4E9", "#E69F00", "#0072B2", "#94A3B8"],
                ),
                title="ASHA Rung",
            ),
            size=alt.Size("total_cells_sampled:Q", title="Sampled Cells", scale=alt.Scale(range=[40, 200])),
            tooltip=[
                "trial_id:Q",
                "status:N",
                "mean_loco_auc:Q",
                "collinearity_max:Q",
                "n_cell_states:Q",
                "total_cells_sampled:Q",
            ],
        )

        # Highlight Pareto non-dominated trials
        pareto_df = pdf[pdf["is_pareto"]]
        pareto_points = (
            alt.Chart(pareto_df)
            .mark_point(shape="diamond", size=180, stroke="#D55E00", strokeWidth=2.0, fill="none")
            .encode(
                x="collinearity_max:Q",
                y="mean_loco_auc:Q",
                tooltip=["trial_id:Q", "mean_loco_auc:Q", "collinearity_max:Q"],
            )
        )

        chart = (
            (all_points + pareto_points)
            .properties(
                title=alt.TitleParams(
                    text="Reference Sampling HPO: Multi-Objective Pareto Frontier",
                    subtitle="Leave-One-Cohort-Out AUC vs. Collinearity across Multi-Fidelity ASHA Rungs (Diamonds = Pareto Front)",
                    fontSize=13,
                    fontWeight="bold",
                ),
                width=520,
                height=380,
            )
            .configure_view(strokeWidth=0.75, stroke="#CBD5E1")
            .configure_axis(
                labelFontSize=10,
                titleFontSize=11,
                gridColor="#F1F5F9",
                tickColor="#1E293B",
            )
        )

        chart.save(str(out_file))
        return Success(out_file)
    except Exception as exc:
        return Failure(f"Failed to generate Pareto frontier SVG: {exc}")


def plot_held_out_validation(
    held_out_aucs: Mapping[str, float],
    discovery_mean_auc: float,
    out_file: Path,
) -> Result[Path, str]:
    """Generate publication-standard Altair SVG bar chart for held-out validation cohorts."""
    try:
        out_file.parent.mkdir(parents=True, exist_ok=True)

        rows = [
            {"cohort": cid, "roc_auc": float(auc_val), "type": "Held-Out Test"}
            for cid, auc_val in held_out_aucs.items()
        ]
        rows.append({"cohort": "Discovery Mean", "roc_auc": float(discovery_mean_auc), "type": "Discovery Benchmark"})
        df_plot = pl.DataFrame(rows).to_pandas()

        bars = (
            alt.Chart(df_plot)
            .mark_bar(cornerRadiusTopLeft=0, cornerRadiusTopRight=0)
            .encode(
                x=alt.X("cohort:N", title="Clinical Cohort", sort=None),
                y=alt.Y("roc_auc:Q", title="Out-of-Cohort ROC-AUC", scale=alt.Scale(domain=[0.4, 0.9])),
                color=alt.Color(
                    "type:N",
                    scale=alt.Scale(
                        domain=["Discovery Benchmark", "Held-Out Test"],
                        range=["#0072B2", "#009E73"],
                    ),
                    title="Evaluation Set",
                ),
                tooltip=["cohort:N", "roc_auc:Q", "type:N"],
            )
        )

        # Baseline reference rule at 0.50 (chance)
        rule = (
            alt.Chart(pl.DataFrame([{"y": 0.50}]).to_pandas())
            .mark_rule(strokeDash=[4, 4], stroke="#94A3B8", strokeWidth=1.0)
            .encode(y="y:Q")
        )

        chart = (
            (bars + rule)
            .properties(
                title=alt.TitleParams(
                    text="Independent Generalization: Held-Out iAtlas Cohorts",
                    subtitle="Validation of Best Deconvolution Reference & Classifier on Locked Test Cohorts",
                    fontSize=13,
                    fontWeight="bold",
                ),
                width=360,
                height=300,
            )
            .configure_view(strokeWidth=0.75, stroke="#CBD5E1")
            .configure_axis(
                labelFontSize=10,
                titleFontSize=11,
                gridColor="#F1F5F9",
            )
        )

        chart.save(str(out_file))
        return Success(out_file)
    except Exception as exc:
        return Failure(f"Failed to generate held-out validation SVG: {exc}")
