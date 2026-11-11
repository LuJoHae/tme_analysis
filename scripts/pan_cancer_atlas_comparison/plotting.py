"""Declarative Altair vector figures for Pan-Cancer T-Cell, CAF & Myeloid Triad validation.

Adheres strictly to Nature Methods minimalist wireframe standards (.agents/rules/code-style-guide.md):
- Pure `#FFFFFF` canvas, hairline `#CBD5E1` borders, `rx="0"`, zero drop shadows
- Okabe-Ito colorblind-safe palette strictly on data marks
- Standalone Unicode math symbols (`β`, `ρ`, `×`) without combining accents
"""

from __future__ import annotations

import altair as alt
import polars as pl
from returns.result import Failure, Result, Success

from plotting_utils import (
    COLOR_HAIRLINE,
    COLOR_NEUTRAL_GREY,
    COLOR_TEXT_PRIMARY,
    OKABE_BLUE,
    OKABE_BLUISH_GREEN,
    OKABE_ORANGE,
    OKABE_REDDISH_PURPLE,
    OKABE_SKY_BLUE,
    OKABE_VERMILION,
    apply_nature_methods_theme,
)

from .core import format_significance_stars
from .models import (
    AltairChart,
    IciDirection,
    META42_CAF_AND_MESENCHYMAL_CONTEXT_IDS,
    META42_T_NK_IDS,
    ProfileOrigin,
    SignificanceTier,
    TriadCompartment,
)

_TRIAD_ROW_COMPARTMENT_ORDER: dict[str, int] = {
    TriadCompartment.T_NK.value: 0,
    TriadCompartment.STROMAL_CAF.value: 1,
    TriadCompartment.MYELOID_TAM.value: 2,
    TriadCompartment.B_PLASMA.value: 3,
    TriadCompartment.ENDOTHELIAL.value: 4,
    TriadCompartment.EPITHELIAL_TUMOR.value: 5,
    TriadCompartment.OTHER.value: 6,
}

_TRIAD_COMPARTMENT_COLORS: tuple[tuple[str, str], ...] = (
    (TriadCompartment.T_NK.value, OKABE_BLUE),
    (TriadCompartment.STROMAL_CAF.value, OKABE_ORANGE),
    (TriadCompartment.MYELOID_TAM.value, OKABE_BLUISH_GREEN),
    (TriadCompartment.B_PLASMA.value, OKABE_SKY_BLUE),
    (TriadCompartment.ENDOTHELIAL.value, OKABE_REDDISH_PURPLE),
    (TriadCompartment.EPITHELIAL_TUMOR.value, COLOR_NEUTRAL_GREY),
    (TriadCompartment.OTHER.value, COLOR_NEUTRAL_GREY),
)


def _stars_text(pval: float) -> str:
    tier = format_significance_stars(pval)
    return "" if tier == SignificanceTier.NS else tier.value


def _select_informative_external_columns(
    sim_sub_df: pl.DataFrame,
    max_per_atlas: int = 10,
) -> list[str]:
    """Select all top-1 matched external substates plus highest-variance substates per atlas."""
    if sim_sub_df.is_empty():
        return []

    chosen: list[str] = []
    for aid in sim_sub_df["atlas_id"].unique(maintain_order=True).to_list():
        a_df = sim_sub_df.filter(pl.col("atlas_id") == aid)
        top1_disp = (
            a_df.filter(pl.col("rank_within_atlas") == 1)["external_display"]
            .unique(maintain_order=True)
            .to_list()
        )
        peak_df = (
            a_df.group_by("external_display")
            .agg(pl.col("composite_substate_score").max().alias("peak_score"))
            .sort("peak_score", descending=True)
        )
        merged_a: list[str] = list(top1_disp)
        for cand in peak_df["external_display"].to_list():
            if len(merged_a) >= max_per_atlas:
                break
            if str(cand) not in merged_a:
                merged_a.append(str(cand))
        for item in merged_a[:max_per_atlas]:
            if str(item) not in chosen:
                chosen.append(str(item))
    return chosen


def _compartment_meta_ids(comp: str) -> tuple[str, ...]:
    """Return the ordered Meta42 reference state IDs shown for an atlas compartment."""
    return (
        META42_T_NK_IDS
        if comp == TriadCompartment.T_NK.value
        else META42_CAF_AND_MESENCHYMAL_CONTEXT_IDS
    )


def build_atlas_similarity_and_ici_chart(
    atlas_id: str,
    full_sim_df: pl.DataFrame,
    atlas_ici_df: pl.DataFrame,
    meta_beta_df: pl.DataFrame,
) -> Result[AltairChart, str]:
    """Single-atlas figure: ICB β bars (a), one-to-one Meta42 alignment heatmap (b), Meta42 bulk β (c)."""
    a_sim_all = full_sim_df.filter(pl.col("atlas_id") == atlas_id)
    if a_sim_all.is_empty():
        return Failure(f"No similarity rows for atlas {atlas_id}")
    comp = str(a_sim_all["atlas_compartment"][0])
    meta_ids = _compartment_meta_ids(comp)
    a_sim = a_sim_all.filter(pl.col("meta_state").is_in(list(meta_ids)))
    a_ici = atlas_ici_df.filter(pl.col("atlas_id") == atlas_id).unique(
        subset=["substate_display"], maintain_order=True
    )
    if a_sim.is_empty() or a_ici.is_empty():
        return Failure(f"Empty similarity or ICB table for atlas {atlas_id}")

    cite = str(a_sim["atlas_citation"][0])
    lineage_word = "T-Cell" if comp == TriadCompartment.T_NK.value else "CAF"
    pan_word = "Pan-T" if comp == TriadCompartment.T_NK.value else "Pan-CAF"

    row_ord = {st: i for i, st in enumerate(meta_ids)}
    matches = a_sim.filter(pl.col("is_one_to_one_match")).with_columns(
        pl.col("meta_state").replace_strict(row_ord, default=99).alias("row_ord")
    )
    matched_cols = matches.sort("row_ord")["external_display"].to_list()
    other_cols = sorted(
        c for c in a_sim["external_display"].unique().to_list() if c not in set(matched_cols)
    )
    col_order: list[str] = [str(c) for c in matched_cols] + [str(c) for c in other_cols]

    plot_df = a_sim.with_columns(
        [
            pl.col("meta_state").replace_strict(row_ord, default=99).alias("row_ord"),
            pl.when(pl.col("is_one_to_one_match"))
            .then(pl.lit("●"))
            .otherwise(pl.lit(""))
            .alias("match_badge"),
        ]
    ).sort("row_ord")
    y_domain = (
        plot_df.select(["row_ord", "meta_state_label"])
        .unique(subset=["meta_state_label"], maintain_order=True)
        .sort("row_ord")["meta_state_label"]
        .to_list()
    )

    ici_plot = a_ici.with_columns(
        pl.Series(
            "star_txt",
            [_stars_text(float(q)) for q in a_ici["lineage_adj_fe_fdr"].to_list()],
            dtype=pl.String,
        )
    )
    w_px = max(320, 24 * len(col_order))
    h_px = max(120, 24 * len(y_domain))

    zero_h = (
        alt.Chart(pl.DataFrame({"y0": [0.0]}))
        .mark_rule(strokeDash=[3, 3], color=COLOR_NEUTRAL_GREY, strokeWidth=0.75)
        .encode(y="y0:Q")
    )
    top_bars = (
        alt.Chart(ici_plot)
        .mark_bar(size=14)
        .encode(
            x=alt.X("substate_display:N", sort=col_order, axis=alt.Axis(labels=False, ticks=False, title="")),
            y=alt.Y(
                "lineage_adj_fe_beta:Q",
                title=f"{pan_word} Adj. β_ICB",
                axis=alt.Axis(tickCount=5, titleFontSize=9.5, labelFontSize=8.5),
            ),
            color=alt.Color(
                "ici_direction:N",
                title="11-Cohort Pan-iAtlas ICB (N = 1,058)",
                scale=alt.Scale(
                    domain=[
                        IciDirection.PRO_RESPONSE.value,
                        IciDirection.PRO_RESISTANCE.value,
                        IciDirection.NEUTRAL.value,
                    ],
                    range=[OKABE_BLUISH_GREEN, OKABE_VERMILION, COLOR_NEUTRAL_GREY],
                ),
                legend=alt.Legend(orient="top", direction="horizontal", titleLimit=420, labelLimit=260),
            ),
        )
    )
    top_stars = (
        alt.Chart(ici_plot.filter(pl.col("star_txt") != ""))
        .mark_text(fontSize=8.5, fontWeight="bold", dy=-5, color=COLOR_TEXT_PRIMARY)
        .encode(
            x=alt.X("substate_display:N", sort=col_order),
            y=alt.Y("lineage_adj_fe_beta:Q"),
            text="star_txt:N",
        )
    )
    top_panel = (zero_h + top_bars + top_stars).properties(
        width=w_px,
        height=95,
        title=f"a. {cite} {lineage_word} Substate ICB Association ({pan_word} Load Adjusted; BH FDR * <0.05, ** <0.01, *** <1e-3, **** <1e-4)",
    )

    heat_rects = (
        alt.Chart(plot_df)
        .mark_rect(stroke=COLOR_HAIRLINE, strokeWidth=0.35)
        .encode(
            x=alt.X(
                "external_display:N",
                sort=col_order,
                title=f"{cite} Substates (matched pairs first)",
                axis=alt.Axis(labelAngle=-42, labelFontSize=8.5, titleFontSize=10, labelLimit=185),
            ),
            y=alt.Y(
                "meta_state_label:N",
                sort=y_domain,
                title=f"Consensus-Meta42 {lineage_word} States ({len(y_domain)})",
                axis=alt.Axis(labelFontSize=9, titleFontSize=10, labelLimit=215),
            ),
            color=alt.Color(
                "composite_substate_score:Q",
                title="Hybrid Proximity",
                scale=alt.Scale(
                    domain=[-0.55, 0.0, 0.55],
                    range=[OKABE_VERMILION, "#FFFFFF", OKABE_BLUE],
                    clamp=True,
                ),
                legend=alt.Legend(orient="bottom", gradientLength=170),
            ),
        )
    )
    heat_badges = (
        alt.Chart(plot_df.filter(pl.col("match_badge") != ""))
        .mark_text(fontSize=11, fontWeight="bold")
        .encode(
            x=alt.X("external_display:N", sort=col_order),
            y=alt.Y("meta_state_label:N", sort=y_domain),
            text="match_badge:N",
            color=alt.condition(
                alt.datum.composite_substate_score > 0.30,
                alt.value("#FFFFFF"),
                alt.value(COLOR_TEXT_PRIMARY),
            ),
        )
    )
    center_panel = (heat_rects + heat_badges).properties(
        width=w_px,
        height=h_px,
        title="b. One-to-One Meta42 ↔ Atlas Substate Assignment (● = optimal unique pair)",
    )

    side_df = (
        meta_beta_df.filter(pl.col("meta_state_label").is_in(y_domain))
        .unique(subset=["meta_state_label"], maintain_order=True)
        .with_columns(
            pl.when(pl.col("meta42_bulk_beta") >= 0.0)
            .then(pl.lit("Responder (β ≥ 0)"))
            .otherwise(pl.lit("Non-Responder (β < 0)"))
            .alias("meta_dir")
        )
    )
    zero_v = (
        alt.Chart(pl.DataFrame({"x0": [0.0]}))
        .mark_rule(strokeDash=[3, 3], color=COLOR_NEUTRAL_GREY, strokeWidth=0.75)
        .encode(x="x0:Q")
    )
    side_bars = (
        alt.Chart(side_df)
        .mark_bar(size=12)
        .encode(
            y=alt.Y("meta_state_label:N", sort=y_domain, axis=alt.Axis(labels=False, ticks=False, title="")),
            x=alt.X("meta42_bulk_beta:Q", title="Meta42 Bulk β", axis=alt.Axis(tickCount=4, labelFontSize=8.5)),
            color=alt.Color(
                "meta_dir:N",
                title="Meta42 Outcome",
                scale=alt.Scale(
                    domain=["Responder (β ≥ 0)", "Non-Responder (β < 0)"],
                    range=[OKABE_BLUISH_GREEN, OKABE_VERMILION],
                ),
                legend=alt.Legend(orient="bottom", direction="vertical"),
            ),
        )
    )
    right_panel = (zero_v + side_bars).properties(width=130, height=h_px, title="c. Meta42 ICB β")

    combined = alt.vconcat(
        top_panel,
        alt.hconcat(center_panel, right_panel, spacing=8).resolve_scale(color="independent"),
        spacing=8,
    ).resolve_scale(color="independent")
    return Success(apply_nature_methods_theme(combined))



def build_triad_all42_specificity_chart(
    full_sim_df: pl.DataFrame,
    triad_summary_df: pl.DataFrame,
) -> Result[AltairChart, str]:
    """Figure S-Triad: Unified 42-State `Consensus-Meta42` x Myeloid/TAM + T/NK + Stromal/CAF Triad Specificity."""
    if full_sim_df.is_empty() or triad_summary_df.is_empty():
        return Failure("Cannot build unified Triad chart: empty similarity or summary DataFrame")

    col_order: list[str] = []
    for comp in (TriadCompartment.T_NK, TriadCompartment.STROMAL_CAF, TriadCompartment.MYELOID_TAM):
        sub_c_raw = full_sim_df.filter(
            (pl.col("atlas_compartment") == comp.value)
            & (pl.col("external_cluster") != "20_TDoub")
        )
        sub_c_emp = sub_c_raw.filter(pl.col("origin") == ProfileOrigin.EMPIRICAL_SINGLE_CELL_H5AD.value)
        sub_c = sub_c_emp if not sub_c_emp.is_empty() else sub_c_raw
        for col_lbl in _select_informative_external_columns(sub_c, max_per_atlas=5):
            if col_lbl not in col_order:
                col_order.append(col_lbl)

    if not col_order:
        return Failure("No Triad columns selected for unified 42-state matrix")

    row_ordered_summary = triad_summary_df.with_columns(
        pl.col("meta_compartment").replace_strict(_TRIAD_ROW_COMPARTMENT_ORDER, default=99).alias("comp_ord")
    ).sort(["comp_ord", "cognate_specificity_margin"], descending=[False, True])
    y_domain = row_ordered_summary["meta_state_label"].to_list()

    triad_plot = full_sim_df.filter(pl.col("external_display").is_in(col_order))
    w_px = max(580, 18 * len(col_order))

    heat_all42 = (
        alt.Chart(triad_plot)
        .mark_rect(stroke=COLOR_HAIRLINE, strokeWidth=0.25)
        .encode(
            x=alt.X(
                "external_display:N",
                sort=col_order,
                title="Representative Pan-Cancer Triad Substates (T/NK · Stromal/CAF · Myeloid/TAM)",
                axis=alt.Axis(labelAngle=-45, labelFontSize=8, titleFontSize=10, labelLimit=180),
            ),
            y=alt.Y(
                "meta_state_label:N",
                sort=y_domain,
                title="All 42 Consensus-Meta42 Deconvolution States (Ordered by Lineage Compartment)",
                axis=alt.Axis(labelFontSize=8, titleFontSize=10, labelLimit=215),
            ),
            color=alt.Color(
                "hybrid_pan_affinity:Q",
                title="Pan-Lineage Affinity",
                scale=alt.Scale(
                    domain=[-0.50, 0.0, 0.55],
                    range=[OKABE_VERMILION, "#FFFFFF", OKABE_BLUE],
                    clamp=True,
                ),
                legend=alt.Legend(orient="bottom", gradientLength=180),
            ),
        )
        .properties(
            width=w_px,
            height=630,
            title="a. Unified 42-State Block-Diagonal Specificity Across T/NK, Stromal/CAF, and Myeloid/TAM Reference Atlases",
        )
    )

    comp_domain = [k for k, _ in _TRIAD_COMPARTMENT_COLORS]
    comp_range = [v for _, v in _TRIAD_COMPARTMENT_COLORS]
    margin_bars = (
        alt.Chart(row_ordered_summary)
        .mark_bar(size=10)
        .encode(
            y=alt.Y("meta_state_label:N", sort=y_domain, axis=alt.Axis(labels=False, ticks=False, title="")),
            x=alt.X(
                "best_hybrid_pan_affinity:Q",
                title="Peak Triad Affinity",
                axis=alt.Axis(tickCount=4, labelFontSize=8.5),
            ),
            color=alt.Color(
                "meta_compartment:N",
                title="Meta42 Lineage",
                scale=alt.Scale(domain=comp_domain, range=comp_range),
                legend=alt.Legend(orient="right"),
            ),
        )
        .properties(
            width=145,
            height=630,
            title="b. Peak Affinity",
        )
    )

    combined_triad = alt.hconcat(heat_all42, margin_bars, spacing=8).resolve_scale(color="independent")
    return Success(apply_nature_methods_theme(combined_triad))
