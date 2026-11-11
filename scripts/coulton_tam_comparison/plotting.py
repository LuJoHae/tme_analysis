"""Declarative Altair vector visualizations for Meta42 <-> Coulton_TAM_Atlas comparison.

Adheres strictly to Nature Methods minimalist wireframe standards (.agents/rules/code-style-guide.md):
- Pure `#FFFFFF` canvas, hairline `#CBD5E1` borders, zero drop shadows
- Okabe-Ito colorblind-safe palette strictly on data marks
- Vector `.svg` + 300 DPI `.png` export via `plotting_utils`
"""

from __future__ import annotations

from collections.abc import Sequence
import math
from typing import assert_never, cast

import altair as alt
import numpy as np
import polars as pl
from returns.result import Failure, Result, Success

from plotting_utils import (
    COLOR_HAIRLINE,
    COLOR_NEUTRAL_GREY,
    COLOR_TEXT_PRIMARY,
    COLOR_TEXT_SECONDARY,
    OKABE_BLUE,
    OKABE_BLUISH_GREEN,
    OKABE_ORANGE,
    OKABE_SKY_BLUE,
    OKABE_VERMILION,
    apply_nature_methods_theme,
)

from .core import compute_ward_biclustering, strip_coulton_cluster_number
from .models import (
    ALL_MYELOID_META42_IDS,
    AltairChart,
    AssaySource,
    COULTON_LITERATURE_CATALOG,
    COULTON_ORDERED_CLUSTERS,
    CoultonClusterRole,
    DendrogramBracket,
    IciOutcomeDirection,
    LineageCompartment,
    MACROPHAGE_CONTINUUM_META42_IDS,
    SignificanceTier,
    TissueStratum,
)

_MYELOID_DISPLAY_ORDER: tuple[str, ...] = (
    "M14",
    "M42",
    "M32",
    "M10",
    "M28",
    "M35",
    "M19",
    "M41",
    "M09",
    "M31",
    "M39",
    "M38",
)

_COMPARTMENT_SORT_ORDER: dict[str, int] = {
    LineageCompartment.MYELOID_MAST.value: 0,
    LineageCompartment.T_NK.value: 1,
    LineageCompartment.B_CELL.value: 2,
    LineageCompartment.PLASMA.value: 3,
    LineageCompartment.STROMAL_CAF.value: 4,
    LineageCompartment.ENDOTHELIAL.value: 5,
    LineageCompartment.EPITHELIAL_TUMOR.value: 6,
    LineageCompartment.OTHER_UNSPECIFIED.value: 7,
}

_COMPARTMENT_DISPLAY_LABELS: dict[str, str] = {
    LineageCompartment.MYELOID_MAST.value: "Myeloid & TAM",
    LineageCompartment.T_NK.value: "T & NK Cells",
    LineageCompartment.B_CELL.value: "B Cells",
    LineageCompartment.PLASMA.value: "Plasma Cells",
    LineageCompartment.STROMAL_CAF.value: "Stromal / CAF",
    LineageCompartment.ENDOTHELIAL.value: "Endothelial",
    LineageCompartment.EPITHELIAL_TUMOR.value: "Epithelial / Tumor",
    LineageCompartment.OTHER_UNSPECIFIED.value: "Other / Erythroid",
}


def classify_significance_tier(pval: float) -> SignificanceTier:
    """Classify a p-value or FDR q-value into the standard 4-tier significance scale."""
    match pval:
        case p if math.isnan(p) or p >= 0.05:
            return SignificanceTier.NS
        case p if p < 0.0001:
            return SignificanceTier.FOUR_STARS
        case p if p < 0.001:
            return SignificanceTier.THREE_STARS
        case p if p < 0.01:
            return SignificanceTier.TWO_STARS
        case _:
            return SignificanceTier.ONE_STAR


def format_significance_stars(pval: float) -> str:
    """Format a p-value or FDR q-value onto '', '*', '**', '***', or '****'."""
    tier = classify_significance_tier(pval)
    match tier:
        case SignificanceTier.NS:
            return ""
        case SignificanceTier.ONE_STAR:
            return "*"
        case SignificanceTier.TWO_STARS:
            return "**"
        case SignificanceTier.THREE_STARS:
            return "***"
        case SignificanceTier.FOUR_STARS:
            return "****"
        case _ as unreachable:
            assert_never(unreachable)


def build_myeloid_coulton_similarity_and_ici_chart(
    full_sim_df: pl.DataFrame,
    nearest_df: pl.DataFrame,
    coulton_ici_df: pl.DataFrame,
) -> Result[AltairChart, str]:
    """Figure 1: Coupled Myeloid Meta42 x Coulton TAM substate similarity heatmap + ICI outcome tracks."""
    sub_sim = full_sim_df.filter(
        (pl.col("assay") == AssaySource.RNA_LOG1P.value)
        & (pl.col("tissue_stratum") == TissueStratum.TUMOR_ONLY.value)
        & pl.col("meta_state").is_in(list(ALL_MYELOID_META42_IDS))
    )
    if sub_sim.is_empty() or coulton_ici_df.is_empty():
        return Failure("Cannot build Figure 1: empty similarity or ICI DataFrame")

    coulton_order = [
        lbl
        for lbl in COULTON_ORDERED_CLUSTERS
        if COULTON_LITERATURE_CATALOG[lbl].role != CoultonClusterRole.UNANNOTATED_ARTIFACT
    ]
    order_map = {st: idx for idx, st in enumerate(_MYELOID_DISPLAY_ORDER)}
    sub_sim_ord = sub_sim.with_columns(
        [
            pl.col("meta_state").replace_strict(order_map, default=99).alias("row_ord"),
            pl.when(
                (pl.col("coulton_role") == CoultonClusterRole.TRUE_TAM_SUBSTATE.value)
                & (pl.col("rank_among_true_tams") <= 3)
            )
            .then(pl.lit("#") + pl.col("rank_among_true_tams").cast(pl.Utf8))
            .otherwise(pl.lit(""))
            .alias("rank_badge"),
        ]
    ).sort(["row_ord", "coulton_cluster_id"])

    y_domain = (
        sub_sim_ord.select(["row_ord", "meta_state_label"])
        .unique(subset=["meta_state_label"], maintain_order=True)
        .sort("row_ord")["meta_state_label"]
        .to_list()
    )

    # Top panel: Coulton substate ICI outcomes (CPI1000+ color & 4-tier FDR stars + empirical Pan-iAtlas adj beta)
    top_base = coulton_ici_df.filter(pl.col("coulton_cluster").is_in(coulton_order))
    cpi_stars = [
        format_significance_stars(float(q)) if bool(sig) else ""
        for q, sig in zip(
            top_base["cpi1000_fgsea_qval"].to_list(),
            top_base["cpi1000_is_significant"].to_list(),
            strict=True,
        )
    ]
    top_df = top_base.with_columns(pl.Series("cpi_star", cpi_stars))
    zero_h = (
        alt.Chart(pl.DataFrame({"y0": [0.0]}))
        .mark_rule(strokeDash=[3, 3], color=COLOR_NEUTRAL_GREY, strokeWidth=0.75)
        .encode(y="y0:Q")
    )
    coulton_bars = (
        alt.Chart(top_df)
        .mark_bar(size=16)
        .encode(
            x=alt.X(
                "coulton_cluster:N",
                sort=coulton_order,
                axis=alt.Axis(labels=False, ticks=False, title=""),
            ),
            y=alt.Y(
                "iatlas_mac_adj_meta_beta:Q",
                title="Coulton Adj. β_ICB",
                scale=alt.Scale(domain=[-0.46, 0.68]),
                axis=alt.Axis(tickCount=5, titleFontSize=10, labelFontSize=9),
            ),
            color=alt.Color(
                "cpi1000_published_direction:N",
                title="Coulton CPI1000+ (N=1,446)",
                scale=alt.Scale(
                    domain=[
                        IciOutcomeDirection.RESPONDER_ENRICHED.value,
                        IciOutcomeDirection.NON_RESPONDER_ENRICHED.value,
                        IciOutcomeDirection.NEUTRAL_OR_NS.value,
                    ],
                    range=[OKABE_BLUISH_GREEN, OKABE_VERMILION, COLOR_NEUTRAL_GREY],
                ),
                legend=alt.Legend(orient="top", direction="horizontal"),
            ),
        )
    )
    coulton_stars_pos = (
        alt.Chart(
            top_df.filter(
                (pl.col("cpi_star") != "") & (pl.col("iatlas_mac_adj_meta_beta") >= 0.0)
            )
        )
        .mark_text(fontSize=9.5, fontWeight="bold", dy=-6, color=COLOR_TEXT_PRIMARY)
        .encode(
            x=alt.X("coulton_cluster:N", sort=coulton_order),
            y=alt.Y("iatlas_mac_adj_meta_beta:Q"),
            text="cpi_star:N",
        )
    )
    coulton_stars_neg = (
        alt.Chart(
            top_df.filter(
                (pl.col("cpi_star") != "") & (pl.col("iatlas_mac_adj_meta_beta") < 0.0)
            )
        )
        .mark_text(fontSize=9.5, fontWeight="bold", dy=10, color=COLOR_TEXT_PRIMARY)
        .encode(
            x=alt.X("coulton_cluster:N", sort=coulton_order),
            y=alt.Y("iatlas_mac_adj_meta_beta:Q"),
            text="cpi_star:N",
        )
    )
    top_panel = (zero_h + coulton_bars + coulton_stars_pos + coulton_stars_neg).properties(
        width=598,
        height=105,
        title="a. Coulton TAM Substate ICI Therapy Effect (Pan-iAtlas Adj. β; CPI1000+ FDR q: * <0.05, ** <0.01, *** <1e-3, **** <1e-4)",
    )

    # Main central heatmap
    heat_rects = (
        alt.Chart(sub_sim_ord)
        .mark_rect(stroke=COLOR_HAIRLINE, strokeWidth=0.35)
        .encode(
            x=alt.X(
                "coulton_cluster:N",
                sort=coulton_order,
                title="Coulton et al. (2024) Pan-Cancer TAM & Monocyte Clusters (Tumor Tissue)",
                axis=alt.Axis(labelAngle=-42, labelFontSize=9.5, titleFontSize=10.5),
            ),
            y=alt.Y(
                "meta_state_label:N",
                sort=y_domain,
                title="Consensus-Meta42 Myeloid & Macrophage States",
                axis=alt.Axis(labelFontSize=9.5, titleFontSize=10.5, labelLimit=220),
            ),
            color=alt.Color(
                "composite_substate_score:Q",
                title="Substate Proximity",
                scale=alt.Scale(
                    domain=[-0.55, 0.0, 0.55],
                    range=[OKABE_VERMILION, "#FFFFFF", OKABE_BLUE],
                    clamp=True,
                ),
                legend=alt.Legend(orient="bottom", gradientLength=180),
            ),
        )
    )
    badge_df = sub_sim_ord.filter(pl.col("rank_badge") != "")
    heat_badges = (
        alt.Chart(badge_df)
        .mark_text(fontSize=8.5, fontWeight="bold")
        .encode(
            x=alt.X("coulton_cluster:N", sort=coulton_order),
            y=alt.Y("meta_state_label:N", sort=y_domain),
            text="rank_badge:N",
            color=alt.condition(
                alt.datum.composite_substate_score > 0.30,
                alt.value("#FFFFFF"),
                alt.value(COLOR_TEXT_PRIMARY),
            ),
        )
    )
    center_panel = (heat_rects + heat_badges).properties(
        width=598,
        height=288,
        title="b. Within-Compartment Transcriptomic Proximity & Top-3 Closest True TAM Substates (#1, #2, #3)",
    )

    # Right marginal track: Meta42 Pan-iAtlas Bulk Deconvolution Beta (95% CI + 4-tier p-value stars)
    meta_side_base = nearest_df.filter(
        pl.col("meta_state").is_in(list(ALL_MYELOID_META42_IDS))
        & (pl.col("rank_among_true_tams") == 1)
    )
    meta_stars = [
        format_significance_stars(float(p))
        for p in meta_side_base["meta42_bulk_pval"].to_list()
    ]
    meta_side_df = meta_side_base.with_columns(
        [
            pl.Series("meta_star", meta_stars),
            (pl.col("meta42_bulk_beta") - 1.96 * pl.col("meta42_bulk_se")).alias("ci_low"),
            (pl.col("meta42_bulk_beta") + 1.96 * pl.col("meta42_bulk_se")).alias("ci_high"),
            pl.when(pl.col("meta42_bulk_beta") >= 0.0)
            .then(pl.lit("Responder (β > 0)"))
            .otherwise(pl.lit("Non-Responder (β < 0)"))
            .alias("meta_dir"),
        ]
    )
    zero_v = (
        alt.Chart(pl.DataFrame({"x0": [0.0]}))
        .mark_rule(strokeDash=[3, 3], color=COLOR_NEUTRAL_GREY, strokeWidth=0.75)
        .encode(x="x0:Q")
    )
    side_ci = (
        alt.Chart(meta_side_df)
        .mark_rule(strokeWidth=1.4)
        .encode(
            y=alt.Y("meta_state_label:N", sort=y_domain, axis=alt.Axis(labels=False, ticks=False, title="")),
            x=alt.X(
                "ci_low:Q",
                title="Meta42 Bulk β (95% CI)",
                scale=alt.Scale(domain=[-2.45, 2.95]),
                axis=alt.Axis(tickCount=5, labelFontSize=9),
            ),
            x2="ci_high:Q",
            color=alt.Color(
                "meta_dir:N",
                title="Meta42 Outcome",
                scale=alt.Scale(
                    domain=["Responder (β > 0)", "Non-Responder (β < 0)"],
                    range=[OKABE_BLUISH_GREEN, OKABE_VERMILION],
                ),
                legend=alt.Legend(orient="bottom", direction="vertical"),
            ),
        )
    )
    side_pts = (
        alt.Chart(meta_side_df)
        .mark_circle(size=54, opacity=1.0)
        .encode(
            y=alt.Y("meta_state_label:N", sort=y_domain),
            x="meta42_bulk_beta:Q",
            color=alt.Color(
                "meta_dir:N",
                scale=alt.Scale(
                    domain=["Responder (β > 0)", "Non-Responder (β < 0)"],
                    range=[OKABE_BLUISH_GREEN, OKABE_VERMILION],
                ),
            ),
        )
    )
    side_stars_pos = (
        alt.Chart(
            meta_side_df.filter(
                (pl.col("meta_star") != "") & (pl.col("meta42_bulk_beta") >= 0.0)
            )
        )
        .mark_text(align="left", dx=4, dy=1, fontSize=9.5, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(
            y=alt.Y("meta_state_label:N", sort=y_domain),
            x=alt.X("ci_high:Q"),
            text="meta_star:N",
        )
    )
    side_stars_neg = (
        alt.Chart(
            meta_side_df.filter(
                (pl.col("meta_star") != "") & (pl.col("meta42_bulk_beta") < 0.0)
            )
        )
        .mark_text(align="right", dx=-4, dy=1, fontSize=9.5, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(
            y=alt.Y("meta_state_label:N", sort=y_domain),
            x=alt.X("ci_low:Q"),
            text="meta_star:N",
        )
    )
    right_panel = (zero_v + side_ci + side_pts + side_stars_pos + side_stars_neg).properties(
        width=205,
        height=288,
        title="c. Meta42 ICB Effect",
    )

    bottom_row = alt.hconcat(center_panel, right_panel, spacing=12).resolve_scale(color="independent")
    full_fig = alt.vconcat(top_panel, bottom_row, spacing=10).resolve_scale(color="independent")
    return Success(apply_nature_methods_theme(full_fig))


def build_nearest_substate_outcome_bridge_chart(
    nearest_df: pl.DataFrame,
) -> Result[AltairChart, str]:
    """Figure 2: Side-by-side Top-3 Coulton TAM matches & triangulated ICI forest plot for Macrophage Continuum."""
    focus_order = ("M14", "M42", "M32", "M09", "M10", "M28", "M35", "M19")
    focus_ids = list(MACROPHAGE_CONTINUUM_META42_IDS) + ["M09"]
    ord_map = {st: idx for idx, st in enumerate(focus_order)}
    raw_sub = nearest_df.filter(pl.col("meta_state").is_in(focus_ids))
    if raw_sub.is_empty():
        return Failure("Cannot build Figure 2: nearest_df has zero macrophage continuum rows")

    short_genes = [
        ", ".join(p.strip() for p in str(g).split(",")[:5] if p.strip())
        for g in raw_sub["shared_signature_genes"].to_list()
    ]
    coulton_adj_stars = [
        format_significance_stars(float(p))
        for p in raw_sub["iatlas_mac_adj_meta_pval"].to_list()
    ]
    meta_bulk_stars = [
        format_significance_stars(float(p))
        for p in raw_sub["meta42_bulk_pval"].to_list()
    ]
    sub = (
        raw_sub.with_columns(
            [
                pl.Series("short_shared_genes", short_genes),
                pl.Series("coulton_adj_star", coulton_adj_stars),
                pl.Series("meta_bulk_star", meta_bulk_stars),
                pl.col("meta_state").replace_strict(ord_map, default=99).alias("m_ord"),
                (
                    pl.col("meta_state")
                    + pl.lit(" #")
                    + pl.col("rank_among_true_tams").cast(pl.Utf8)
                    + pl.lit(": ")
                    + pl.col("coulton_cluster")
                ).alias("pair_row_label"),
                (pl.col("iatlas_mac_adj_meta_beta") - 1.96 * pl.col("iatlas_mac_adj_meta_se")).alias(
                    "coulton_ci_low"
                ),
                (pl.col("iatlas_mac_adj_meta_beta") + 1.96 * pl.col("iatlas_mac_adj_meta_se")).alias(
                    "coulton_ci_high"
                ),
                (pl.col("meta42_bulk_beta") - 1.96 * pl.col("meta42_bulk_se")).alias("meta_ci_low"),
                (pl.col("meta42_bulk_beta") + 1.96 * pl.col("meta42_bulk_se")).alias("meta_ci_high"),
            ]
        )
        .sort(["m_ord", "rank_among_true_tams"])
    )

    pair_order = sub["pair_row_label"].to_list()

    # Panel a: Proximity bar + shared markers text
    prox_bars = (
        alt.Chart(sub)
        .mark_bar(height=12)
        .encode(
            y=alt.Y(
                "pair_row_label:N",
                sort=pair_order,
                title="Meta42 State -> Top-3 Closest Coulton TAM Substates",
                axis=alt.Axis(labelFontSize=9.5, titleFontSize=10.5, labelLimit=240),
            ),
            x=alt.X(
                "composite_substate_score:Q",
                title="Marker-Anchored Substate Proximity Score",
                scale=alt.Scale(domain=[0.0, 0.85]),
                axis=alt.Axis(tickCount=5, labelFontSize=9),
            ),
            color=alt.Color(
                "macrophage_tier:N",
                title="Meta42 Myeloid Tier",
                scale=alt.Scale(
                    domain=[
                        "Tier 1: Core TAM States (M14, M32, M42)",
                        "Tier 2: Monocyte / MDSC Continuum (M10, M19, M28, M35)",
                        "Tier 3: Other Myeloid / DC / Mast (M09, M31, M38, M39, M41)",
                    ],
                    range=[OKABE_BLUE, OKABE_ORANGE, "#CC79A7"],
                ),
                legend=alt.Legend(orient="bottom", direction="vertical", labelLimit=360),
            ),
        )
    )
    marker_txt = (
        alt.Chart(sub)
        .mark_text(align="left", dx=4, fontSize=8, color=COLOR_TEXT_SECONDARY)
        .encode(
            y=alt.Y("pair_row_label:N", sort=pair_order),
            x=alt.X("composite_substate_score:Q"),
            text="short_shared_genes:N",
        )
    )
    panel_a = (prox_bars + marker_txt).properties(
        width=365,
        height=485,
        title="a. Nearest Coulton TAM Substates & Shared Signature Genes",
    )

    # Panel b: Coulton TAM Substate ICI Outcome Forest (Empirical Pan-iAtlas adj beta + CPI1000+ color + 4-tier stars)
    zero_rule_b = (
        alt.Chart(pl.DataFrame({"x0": [0.0]}))
        .mark_rule(strokeDash=[3, 3], color=COLOR_NEUTRAL_GREY, strokeWidth=0.75)
        .encode(x="x0:Q")
    )
    coulton_ci = (
        alt.Chart(sub)
        .mark_rule(strokeWidth=1.4)
        .encode(
            y=alt.Y("pair_row_label:N", sort=pair_order, axis=alt.Axis(labels=False, ticks=False, title="")),
            x=alt.X(
                "coulton_ci_low:Q",
                title="Matched Coulton Substate Adj. β_ICB (95% CI)",
                scale=alt.Scale(domain=[-0.68, 0.85]),
                axis=alt.Axis(tickCount=5, labelFontSize=9),
            ),
            x2="coulton_ci_high:Q",
            color=alt.Color(
                "cpi1000_published_direction:N",
                title="Coulton CPI1000+ Outcome",
                scale=alt.Scale(
                    domain=[
                        IciOutcomeDirection.RESPONDER_ENRICHED.value,
                        IciOutcomeDirection.NON_RESPONDER_ENRICHED.value,
                        IciOutcomeDirection.NEUTRAL_OR_NS.value,
                    ],
                    range=[OKABE_BLUISH_GREEN, OKABE_VERMILION, COLOR_NEUTRAL_GREY],
                ),
                legend=alt.Legend(orient="bottom", direction="vertical", labelLimit=280),
            ),
        )
    )
    coulton_pts = (
        alt.Chart(sub)
        .mark_circle(size=52, opacity=1.0)
        .encode(
            y=alt.Y("pair_row_label:N", sort=pair_order),
            x="iatlas_mac_adj_meta_beta:Q",
            color=alt.Color(
                "cpi1000_published_direction:N",
                scale=alt.Scale(
                    domain=[
                        IciOutcomeDirection.RESPONDER_ENRICHED.value,
                        IciOutcomeDirection.NON_RESPONDER_ENRICHED.value,
                        IciOutcomeDirection.NEUTRAL_OR_NS.value,
                    ],
                    range=[OKABE_BLUISH_GREEN, OKABE_VERMILION, COLOR_NEUTRAL_GREY],
                ),
            ),
        )
    )
    coulton_b_stars_pos = (
        alt.Chart(
            sub.filter(
                (pl.col("coulton_adj_star") != "") & (pl.col("iatlas_mac_adj_meta_beta") >= 0.0)
            )
        )
        .mark_text(align="left", dx=4, dy=1, fontSize=8.5, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(
            y=alt.Y("pair_row_label:N", sort=pair_order),
            x=alt.X("coulton_ci_high:Q"),
            text="coulton_adj_star:N",
        )
    )
    coulton_b_stars_neg = (
        alt.Chart(
            sub.filter(
                (pl.col("coulton_adj_star") != "") & (pl.col("iatlas_mac_adj_meta_beta") < 0.0)
            )
        )
        .mark_text(align="right", dx=-4, dy=1, fontSize=8.5, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(
            y=alt.Y("pair_row_label:N", sort=pair_order),
            x=alt.X("coulton_ci_low:Q"),
            text="coulton_adj_star:N",
        )
    )
    panel_b = (
        zero_rule_b + coulton_ci + coulton_pts + coulton_b_stars_pos + coulton_b_stars_neg
    ).properties(
        width=225,
        height=485,
        title="b. Coulton Substate ICB Outcome",
    )

    # Panel c: Parent Meta42 Bulk Deconv Beta (95% CI + 4-tier p-value stars)
    zero_rule_c = (
        alt.Chart(pl.DataFrame({"x0": [0.0]}))
        .mark_rule(strokeDash=[3, 3], color=COLOR_NEUTRAL_GREY, strokeWidth=0.75)
        .encode(x="x0:Q")
    )
    meta_ci = (
        alt.Chart(sub)
        .mark_rule(strokeWidth=1.4)
        .encode(
            y=alt.Y("pair_row_label:N", sort=pair_order, axis=alt.Axis(labels=False, ticks=False, title="")),
            x=alt.X(
                "meta_ci_low:Q",
                title="Parent Meta42 Bulk β (95% CI)",
                scale=alt.Scale(domain=[-2.45, 2.95]),
                axis=alt.Axis(tickCount=5, labelFontSize=9),
            ),
            x2="meta_ci_high:Q",
            color=alt.Color(
                "ici_outcome_concordance:N",
                title="Joint Outcome Agreement",
                scale=alt.Scale(
                    domain=[
                        "Concordant Pro-Response (+/+)",
                        "Concordant Pro-Resistance (-/-)",
                        "Polarization-Split (+/- or -/+)",
                        "Neutral / Sub-Threshold",
                    ],
                    range=[OKABE_BLUISH_GREEN, OKABE_VERMILION, OKABE_ORANGE, OKABE_SKY_BLUE],
                ),
                legend=alt.Legend(orient="bottom", direction="vertical", labelLimit=280),
            ),
        )
    )
    meta_pts = (
        alt.Chart(sub)
        .mark_square(size=48, opacity=1.0)
        .encode(
            y=alt.Y("pair_row_label:N", sort=pair_order),
            x="meta42_bulk_beta:Q",
            color=alt.Color(
                "ici_outcome_concordance:N",
                scale=alt.Scale(
                    domain=[
                        "Concordant Pro-Response (+/+)",
                        "Concordant Pro-Resistance (-/-)",
                        "Polarization-Split (+/- or -/+)",
                        "Neutral / Sub-Threshold",
                    ],
                    range=[OKABE_BLUISH_GREEN, OKABE_VERMILION, OKABE_ORANGE, OKABE_SKY_BLUE],
                ),
            ),
        )
    )
    meta_c_stars_pos = (
        alt.Chart(
            sub.filter(
                (pl.col("meta_bulk_star") != "") & (pl.col("meta42_bulk_beta") >= 0.0)
            )
        )
        .mark_text(align="left", dx=4, dy=1, fontSize=8.5, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(
            y=alt.Y("pair_row_label:N", sort=pair_order),
            x=alt.X("meta_ci_high:Q"),
            text="meta_bulk_star:N",
        )
    )
    meta_c_stars_neg = (
        alt.Chart(
            sub.filter(
                (pl.col("meta_bulk_star") != "") & (pl.col("meta42_bulk_beta") < 0.0)
            )
        )
        .mark_text(align="right", dx=-4, dy=1, fontSize=8.5, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(
            y=alt.Y("pair_row_label:N", sort=pair_order),
            x=alt.X("meta_ci_low:Q"),
            text="meta_bulk_star:N",
        )
    )
    panel_c = (zero_rule_c + meta_ci + meta_pts + meta_c_stars_pos + meta_c_stars_neg).properties(
        width=215,
        height=485,
        title="c. Parent Meta42 Bulk ICB Effect",
    )

    combined = alt.hconcat(panel_a, panel_b, panel_c, spacing=14).resolve_scale(color="independent")
    return Success(apply_nature_methods_theme(combined))


def build_all42_pan_lineage_coulton_chart(
    full_sim_df: pl.DataFrame,
) -> Result[AltairChart, str]:
    """Figure 3: Complete 42 x 23 pan-lineage Meta42 vs Coulton_TAM_Atlas specificity & lineage validation."""
    pan_df = full_sim_df.filter(
        (pl.col("assay") == AssaySource.RNA_LOG1P.value)
        & (pl.col("tissue_stratum") == TissueStratum.TUMOR_ONLY.value)
    )
    if pan_df.is_empty():
        return Failure("Cannot build Figure 3: empty similarity DataFrame")

    metric_col = (
        "hybrid_pan_affinity"
        if "hybrid_pan_affinity" in pan_df.columns
        else "pan_lineage_spec_cosine"
    )
    coulton_order = [
        lbl
        for lbl in COULTON_ORDERED_CLUSTERS
        if COULTON_LITERATURE_CATALOG[lbl].role != CoultonClusterRole.UNANNOTATED_ARTIFACT
    ]
    pan_ord = pan_df.with_columns(
        [
            pl.col("lineage_compartment")
            .replace_strict(_COMPARTMENT_SORT_ORDER, default=99)
            .alias("comp_ord"),
            pl.col("lineage_compartment")
            .replace_strict(_COMPARTMENT_DISPLAY_LABELS, default="Other / Erythroid")
            .alias("compartment_display"),
        ]
    ).sort(["comp_ord", "meta_state", "coulton_cluster_id"])

    y_all_order = (
        pan_ord.select(["comp_ord", "meta_state", "meta_state_label"])
        .unique(subset=["meta_state"], maintain_order=True)
        .sort(["comp_ord", "meta_state"])["meta_state_label"]
        .to_list()
    )

    heat_all = (
        alt.Chart(pan_ord)
        .mark_rect(stroke=COLOR_HAIRLINE, strokeWidth=0.25)
        .encode(
            x=alt.X(
                "coulton_cluster:N",
                sort=coulton_order,
                title="Coulton et al. (2024) Clusters (21 TAM Substates + 19_ClassMono + 20_TDoub Control)",
                axis=alt.Axis(labelAngle=-42, labelFontSize=9, titleFontSize=10.5),
            ),
            y=alt.Y(
                "meta_state_label:N",
                sort=y_all_order,
                title="All 42 Consensus-Meta42 States (Grouped by Lineage Compartment)",
                axis=alt.Axis(labelFontSize=8.5, titleFontSize=10.5, labelLimit=210),
            ),
            color=alt.Color(
                f"{metric_col}:Q",
                title="Pan-Lineage TAM Affinity",
                scale=alt.Scale(
                    domain=[-0.45, 0.0, 0.55],
                    range=[OKABE_VERMILION, "#FFFFFF", OKABE_BLUE],
                    clamp=True,
                ),
                legend=alt.Legend(orient="bottom", gradientLength=200),
            ),
        )
        .properties(
            width=585,
            height=640,
            title="a. Full 42 x 23 Pan-Lineage Specificity Matrix (Validating Myeloid Enrichment & Cross-Lineage Bridges)",
        )
    )

    max_tam_corr = (
        pan_ord.filter(pl.col("coulton_role") == CoultonClusterRole.TRUE_TAM_SUBSTATE.value)
        .group_by(["comp_ord", "meta_state", "meta_state_label", "compartment_display"])
        .agg(pl.col(metric_col).max().alias("max_tam_affinity"))
        .sort(["comp_ord", "meta_state"])
    )
    zero_v = (
        alt.Chart(pl.DataFrame({"x0": [0.0]}))
        .mark_rule(strokeDash=[3, 3], color=COLOR_NEUTRAL_GREY, strokeWidth=0.75)
        .encode(x="x0:Q")
    )
    max_bars = (
        alt.Chart(max_tam_corr)
        .mark_bar(height=10)
        .encode(
            y=alt.Y("meta_state_label:N", sort=y_all_order, axis=alt.Axis(labels=False, ticks=False, title="")),
            x=alt.X(
                "max_tam_affinity:Q",
                title="Peak True-TAM Affinity",
                axis=alt.Axis(tickCount=4, labelFontSize=9),
            ),
            color=alt.Color(
                "compartment_display:N",
                title="Meta42 Compartment",
                scale=alt.Scale(
                    domain=[
                        "Myeloid & TAM",
                        "T & NK Cells",
                        "B Cells",
                        "Plasma Cells",
                        "Stromal / CAF",
                        "Endothelial",
                        "Epithelial / Tumor",
                        "Other / Erythroid",
                    ],
                    range=[
                        OKABE_ORANGE,
                        OKABE_BLUE,
                        OKABE_BLUISH_GREEN,
                        OKABE_SKY_BLUE,
                        OKABE_VERMILION,
                        "#CC79A7",
                        "#64748B",
                        COLOR_NEUTRAL_GREY,
                    ],
                ),
                legend=alt.Legend(orient="bottom", direction="vertical", columns=2, labelLimit=180),
            ),
        )
    )
    right_bars = (zero_v + max_bars).properties(
        width=190,
        height=640,
        title="b. Peak True-TAM Affinity",
    )

    combined = alt.hconcat(heat_all, right_bars, spacing=12).resolve_scale(color="independent")
    return Success(apply_nature_methods_theme(combined))


def _brackets_to_polyline_df(brackets: Sequence[DendrogramBracket]) -> pl.DataFrame:
    """Convert orthogonal dendrogram U-brackets into a 4-vertex polyline Polars DataFrame."""
    if not brackets:
        return pl.DataFrame(
            {"bracket_id": [0], "pt_ord": [0], "u": [0.5], "h": [0.0]}
        )
    rows: list[dict[str, int | float]] = []
    for br in brackets:
        rows.extend(
            (
                {"bracket_id": br.bracket_id, "pt_ord": 0, "u": br.u_left, "h": br.h_left},
                {"bracket_id": br.bracket_id, "pt_ord": 1, "u": br.u_left, "h": br.h_merge},
                {"bracket_id": br.bracket_id, "pt_ord": 2, "u": br.u_right, "h": br.h_merge},
                {"bracket_id": br.bracket_id, "pt_ord": 3, "u": br.u_right, "h": br.h_right},
            )
        )
    return pl.DataFrame(rows)


def _panel_letter_layer(letter: str, x_px: float, y_px: float) -> alt.Chart:
    """Render a Nature Methods bold lowercase panel letter without Vega-Lite role-title."""
    return cast(
        alt.Chart,
        alt.Chart(pl.DataFrame({"_panel_lbl": [letter]}))
        .mark_text(
            align="left",
            baseline="bottom",
            fontSize=15.0,
            fontWeight="bold",
            color=COLOR_TEXT_PRIMARY,
            clip=False,
        )
        .encode(
            x=alt.value(x_px),
            y=alt.value(y_px),
            text="_panel_lbl:N",
        ),
    )


def _extract_simplified_bicluster_tables(
    full_sim_df: pl.DataFrame,
    nearest_df: pl.DataFrame,
    coulton_ici_df: pl.DataFrame,
) -> Result[tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, list[str], list[str], tuple[DendrogramBracket, ...], tuple[DendrogramBracket, ...]], str]:
    """Filter (>= ** myeloid rows x >= ** Top-3 true TAM columns), mark only #1, strip numbers, and Ward-bicluster."""
    meta_side_all = nearest_df.filter(
        pl.col("meta_state").is_in(list(ALL_MYELOID_META42_IDS))
        & (pl.col("rank_among_true_tams") == 1)
    )
    if meta_side_all.is_empty() or full_sim_df.is_empty() or coulton_ici_df.is_empty():
        return Failure("Cannot build Figure 1b: empty similarity, nearest, or ICI DataFrame")

    all_stars = [
        format_significance_stars(float(p))
        for p in meta_side_all["meta42_bulk_pval"].to_list()
    ]
    meta_side_filtered = meta_side_all.with_columns(
        pl.Series("meta_star", all_stars)
    ).filter(pl.col("meta_star").str.len_chars() >= 2)
    if meta_side_filtered.is_empty():
        return Failure("Cannot build Figure 1b: zero myeloid states reach >= 2 significance asterisks")

    all_cpi_stars = [
        format_significance_stars(float(q)) if bool(sig) else ""
        for q, sig in zip(
            coulton_ici_df["cpi1000_fgsea_qval"].to_list(),
            coulton_ici_df["cpi1000_is_significant"].to_list(),
            strict=True,
        )
    ]
    coulton_ici_with_stars = coulton_ici_df.with_columns(
        pl.Series("cpi_star", all_cpi_stars)
    )
    sig_coulton_clusters = (
        coulton_ici_with_stars.filter(pl.col("cpi_star").str.len_chars() >= 2)["coulton_cluster"]
        .to_list()
    )

    retained_meta_ids = meta_side_filtered["meta_state"].to_list()
    sub_rows = full_sim_df.filter(
        (pl.col("assay") == AssaySource.RNA_LOG1P.value)
        & (pl.col("tissue_stratum") == TissueStratum.TUMOR_ONLY.value)
        & pl.col("meta_state").is_in(retained_meta_ids)
    )
    retained_coulton_raw = (
        sub_rows.filter(
            (pl.col("coulton_role") == CoultonClusterRole.TRUE_TAM_SUBSTATE.value)
            & (pl.col("rank_among_true_tams") <= 3)
            & pl.col("coulton_cluster").is_in(sig_coulton_clusters)
        )["coulton_cluster"]
        .unique(maintain_order=True)
        .to_list()
    )
    if not retained_coulton_raw:
        return Failure("Cannot build Figure 1b: zero >= ** Top-3 Coulton TAM columns retained")

    sub_clean = sub_rows.filter(pl.col("coulton_cluster").is_in(retained_coulton_raw)).with_columns(
        [
            pl.Series(
                "coulton_clean_label",
                [
                    strip_coulton_cluster_number(str(c))
                    for c in sub_rows.filter(pl.col("coulton_cluster").is_in(retained_coulton_raw))[
                        "coulton_cluster"
                    ].to_list()
                ],
            ),
            pl.when(
                (pl.col("coulton_role") == CoultonClusterRole.TRUE_TAM_SUBSTATE.value)
                & (pl.col("rank_among_true_tams") == 1)
            )
            .then(pl.lit("#1"))
            .otherwise(pl.lit(""))
            .alias("rank_badge"),
        ]
    )

    init_rows = sorted(sub_clean["meta_state_label"].unique().to_list())
    init_cols = sorted(sub_clean["coulton_clean_label"].unique().to_list())
    row_to_idx = {r: idx for idx, r in enumerate(init_rows)}
    col_to_idx = {c: idx for idx, c in enumerate(init_cols)}
    score_mat = np.zeros((len(init_rows), len(init_cols)), dtype=np.float64)
    for r_lbl, c_lbl, val in zip(
        sub_clean["meta_state_label"].to_list(),
        sub_clean["coulton_clean_label"].to_list(),
        sub_clean["composite_substate_score"].to_list(),
        strict=True,
    ):
        score_mat[row_to_idx[str(r_lbl)], col_to_idx[str(c_lbl)]] = float(val)

    match compute_ward_biclustering(score_mat, init_rows, init_cols):
        case Failure(b_err):
            return Failure(b_err)
        case Success(bic):
            bicluster = bic
        case _:
            return Failure("Unexpected Result state in compute_ward_biclustering")

    top_raw = coulton_ici_with_stars.filter(pl.col("coulton_cluster").is_in(retained_coulton_raw))
    top_df = top_raw.with_columns(
        pl.Series(
            "coulton_clean_label",
            [strip_coulton_cluster_number(str(c)) for c in top_raw["coulton_cluster"].to_list()],
        )
    )

    return Success(
        (
            sub_clean,
            meta_side_filtered,
            top_df,
            list(bicluster.ordered_rows),
            list(bicluster.ordered_cols),
            bicluster.row_brackets,
            bicluster.col_brackets,
        )
    )


def build_simplified_biclustered_myeloid_coulton_chart(
    full_sim_df: pl.DataFrame,
    nearest_df: pl.DataFrame,
    coulton_ici_df: pl.DataFrame,
    show_panel_letters: bool = True,
) -> Result[AltairChart, str]:
    """Figure 1b: Heading-free paper-ready Ward-biclustered Myeloid (>= **) x Coulton TAM (>= **) figure."""
    match _extract_simplified_bicluster_tables(full_sim_df, nearest_df, coulton_ici_df):
        case Failure(err):
            return Failure(err)
        case Success(extracted):
            sub_clean, meta_side_filtered, top_df, y_domain, x_domain, row_brackets, col_brackets = extracted
        case _:
            return Failure("Unexpected Result state in _extract_simplified_bicluster_tables")

    n_rows = len(y_domain)
    n_cols = len(x_domain)
    heat_w = max(272, 68 * n_cols)
    heat_h = max(144, 42 * n_rows)

    # Top marginal dendrogram for Coulton TAM columns (heading-free; zero role-title)
    col_dendro_df = _brackets_to_polyline_df(col_brackets)
    col_dendro_lines = (
        alt.Chart(col_dendro_df)
        .mark_line(color=COLOR_TEXT_PRIMARY, strokeWidth=1.05)
        .encode(
            x=alt.X(
                "u:Q",
                scale=alt.Scale(domain=[0.0, float(n_cols)], nice=False, zero=False),
                axis=alt.Axis(labels=False, ticks=False, domain=False, grid=False, title=""),
            ),
            y=alt.Y(
                "h:Q",
                scale=alt.Scale(domain=[0.0, 1.08], nice=False, zero=False),
                axis=alt.Axis(labels=False, ticks=False, domain=False, grid=False, title=""),
            ),
            detail="bracket_id:N",
            order="pt_ord:Q",
        )
    )
    col_dendro = alt.layer(
        *(
            [col_dendro_lines, _panel_letter_layer("a", x_px=-58.0, y_px=15.0)]
            if show_panel_letters
            else [col_dendro_lines]
        ),
        view=alt.ViewConfig(strokeOpacity=0),
    ).properties(width=heat_w, height=44)

    # Panel a: Coulton TAM Substate ICI outcomes over retained & biclustered (>= **) columns
    zero_h = (
        alt.Chart(pl.DataFrame({"y0": [0.0]}))
        .mark_rule(strokeDash=[3, 3], color=COLOR_NEUTRAL_GREY, strokeWidth=0.75)
        .encode(y="y0:Q")
    )
    coulton_bars = (
        alt.Chart(top_df)
        .mark_bar(size=32)
        .encode(
            x=alt.X(
                "coulton_clean_label:N",
                sort=x_domain,
                scale=alt.Scale(paddingInner=0.2, paddingOuter=0.1),
                axis=alt.Axis(labels=False, ticks=False, title=""),
            ),
            y=alt.Y(
                "iatlas_mac_adj_meta_beta:Q",
                title="Coulton Adj. β_ICB",
                scale=alt.Scale(domain=[-0.48, 0.70]),
                axis=alt.Axis(tickCount=5, titleFontSize=11, labelFontSize=10),
            ),
            color=alt.Color(
                "cpi1000_published_direction:N",
                title="Coulton CPI1000+ (N=1,446)",
                scale=alt.Scale(
                    domain=[
                        IciOutcomeDirection.RESPONDER_ENRICHED.value,
                        IciOutcomeDirection.NON_RESPONDER_ENRICHED.value,
                    ],
                    range=[OKABE_BLUISH_GREEN, OKABE_VERMILION],
                ),
                legend=alt.Legend(
                    orient="right",
                    direction="vertical",
                    labelLimit=240,
                    titleFontSize=10.5,
                    labelFontSize=10,
                ),
            ),
        )
    )
    coulton_stars_pos = (
        alt.Chart(
            top_df.filter(
                (pl.col("cpi_star") != "") & (pl.col("iatlas_mac_adj_meta_beta") >= 0.0)
            )
        )
        .mark_text(fontSize=10.5, fontWeight="bold", dy=-7, color=COLOR_TEXT_PRIMARY)
        .encode(
            x=alt.X("coulton_clean_label:N", sort=x_domain),
            y=alt.Y("iatlas_mac_adj_meta_beta:Q"),
            text="cpi_star:N",
        )
    )
    coulton_stars_neg = (
        alt.Chart(
            top_df.filter(
                (pl.col("cpi_star") != "") & (pl.col("iatlas_mac_adj_meta_beta") < 0.0)
            )
        )
        .mark_text(fontSize=10.5, fontWeight="bold", dy=11, color=COLOR_TEXT_PRIMARY)
        .encode(
            x=alt.X("coulton_clean_label:N", sort=x_domain),
            y=alt.Y("iatlas_mac_adj_meta_beta:Q"),
            text="cpi_star:N",
        )
    )
    top_panel = (zero_h + coulton_bars + coulton_stars_pos + coulton_stars_neg).properties(
        width=heat_w,
        height=108,
    )

    # Panel b: Biclustered central heatmap (heading-free; zero role-title; only #1 marked)
    heat_rects = (
        alt.Chart(sub_clean)
        .mark_rect(stroke=COLOR_HAIRLINE, strokeWidth=0.5)
        .encode(
            x=alt.X(
                "coulton_clean_label:N",
                sort=x_domain,
                scale=alt.Scale(paddingInner=0.0, paddingOuter=0.0),
                title="Coulton et al. (2024) Significant TAM Substates (q < 0.01)",
                axis=alt.Axis(labelAngle=-35, labelFontSize=10.5, titleFontSize=11),
            ),
            y=alt.Y(
                "meta_state_label:N",
                sort=y_domain,
                scale=alt.Scale(paddingInner=0.0, paddingOuter=0.0),
                title="Significant Myeloid States (p < 0.01)",
                axis=alt.Axis(labelFontSize=10.5, titleFontSize=11, labelLimit=240),
            ),
            color=alt.Color(
                "composite_substate_score:Q",
                title="Substate Proximity",
                scale=alt.Scale(
                    domain=[-0.55, 0.0, 0.55],
                    range=[OKABE_VERMILION, "#FFFFFF", OKABE_BLUE],
                    clamp=True,
                ),
                legend=alt.Legend(
                    orient="bottom",
                    gradientLength=185,
                    titleFontSize=10.5,
                    labelFontSize=10,
                ),
            ),
        )
    )
    badge_df = sub_clean.filter(pl.col("rank_badge") != "")
    heat_badges = (
        alt.Chart(badge_df)
        .mark_text(fontSize=10.5, fontWeight="bold")
        .encode(
            x=alt.X("coulton_clean_label:N", sort=x_domain),
            y=alt.Y("meta_state_label:N", sort=y_domain),
            text="rank_badge:N",
            color=alt.condition(
                alt.datum.composite_substate_score > 0.30,
                alt.value("#FFFFFF"),
                alt.value(COLOR_TEXT_PRIMARY),
            ),
        )
    )
    center_layers = (
        heat_rects + heat_badges + _panel_letter_layer("b", x_px=-58.0, y_px=-8.0)
        if show_panel_letters
        else heat_rects + heat_badges
    )
    center_panel = center_layers.properties(
        width=heat_w,
        height=heat_h,
    )

    # Panel c: Meta42 ICB Effect (>= ** states, heading-free; zero role-title)
    meta_side_df = meta_side_filtered.with_columns(
        [
            (pl.col("meta42_bulk_beta") - 1.96 * pl.col("meta42_bulk_se")).alias("ci_low"),
            (pl.col("meta42_bulk_beta") + 1.96 * pl.col("meta42_bulk_se")).alias("ci_high"),
            pl.when(pl.col("meta42_bulk_beta") >= 0.0)
            .then(pl.lit("Responder (β > 0)"))
            .otherwise(pl.lit("Non-Responder (β < 0)"))
            .alias("meta_dir"),
        ]
    )
    zero_v = (
        alt.Chart(pl.DataFrame({"x0": [0.0]}))
        .mark_rule(strokeDash=[3, 3], color=COLOR_NEUTRAL_GREY, strokeWidth=0.75)
        .encode(x="x0:Q")
    )
    side_ci = (
        alt.Chart(meta_side_df)
        .mark_rule(strokeWidth=1.6)
        .encode(
            y=alt.Y(
                "meta_state_label:N",
                sort=y_domain,
                scale=alt.Scale(paddingInner=0.2, paddingOuter=0.1),
                axis=alt.Axis(labels=False, ticks=False, title=""),
            ),
            x=alt.X(
                "ci_low:Q",
                title="Meta42 Bulk β (95% CI)",
                scale=alt.Scale(domain=[-2.50, 3.00]),
                axis=alt.Axis(tickCount=5, labelFontSize=10, titleFontSize=11),
            ),
            x2="ci_high:Q",
            color=alt.Color(
                "meta_dir:N",
                title="Meta42 Outcome",
                scale=alt.Scale(
                    domain=["Responder (β > 0)", "Non-Responder (β < 0)"],
                    range=[OKABE_BLUISH_GREEN, OKABE_VERMILION],
                ),
                legend=alt.Legend(
                    orient="bottom",
                    direction="vertical",
                    titleFontSize=10.5,
                    labelFontSize=10,
                ),
            ),
        )
    )
    side_pts = (
        alt.Chart(meta_side_df)
        .mark_circle(size=66, opacity=1.0)
        .encode(
            y=alt.Y("meta_state_label:N", sort=y_domain),
            x="meta42_bulk_beta:Q",
            color=alt.Color(
                "meta_dir:N",
                scale=alt.Scale(
                    domain=["Responder (β > 0)", "Non-Responder (β < 0)"],
                    range=[OKABE_BLUISH_GREEN, OKABE_VERMILION],
                ),
            ),
        )
    )
    side_stars_pos = (
        alt.Chart(
            meta_side_df.filter(
                (pl.col("meta_star") != "") & (pl.col("meta42_bulk_beta") >= 0.0)
            )
        )
        .mark_text(align="left", dx=5, dy=1, fontSize=10.5, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(
            y=alt.Y("meta_state_label:N", sort=y_domain),
            x=alt.X("ci_high:Q"),
            text="meta_star:N",
        )
    )
    side_stars_neg = (
        alt.Chart(
            meta_side_df.filter(
                (pl.col("meta_star") != "") & (pl.col("meta42_bulk_beta") < 0.0)
            )
        )
        .mark_text(align="right", dx=-5, dy=1, fontSize=10.5, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(
            y=alt.Y("meta_state_label:N", sort=y_domain),
            x=alt.X("ci_low:Q"),
            text="meta_star:N",
        )
    )
    right_base = zero_v + side_ci + side_pts + side_stars_pos + side_stars_neg
    right_layers = (
        right_base + _panel_letter_layer("c", x_px=-2.0, y_px=-8.0)
        if show_panel_letters
        else right_base
    )
    right_panel = right_layers.properties(
        width=210,
        height=heat_h,
    )

    # Right marginal dendrogram for Myeloid Meta42 rows
    row_dendro_df = _brackets_to_polyline_df(row_brackets)
    row_dendro = (
        alt.Chart(row_dendro_df, view=alt.ViewConfig(strokeOpacity=0))
        .mark_line(color=COLOR_TEXT_PRIMARY, strokeWidth=1.05)
        .encode(
            x=alt.X(
                "h:Q",
                scale=alt.Scale(domain=[0.0, 1.08], nice=False, zero=False),
                axis=alt.Axis(labels=False, ticks=False, domain=False, grid=False, title=""),
            ),
            y=alt.Y(
                "u:Q",
                scale=alt.Scale(domain=[float(n_rows), 0.0], nice=False, zero=False),
                axis=alt.Axis(labels=False, ticks=False, domain=False, grid=False, title=""),
            ),
            detail="bracket_id:N",
            order="pt_ord:Q",
        )
        .properties(width=46, height=heat_h)
    )

    upper_stack = alt.vconcat(col_dendro, top_panel, spacing=2).resolve_scale(color="independent")
    bottom_row = alt.hconcat(center_panel, right_panel, row_dendro, spacing=10).resolve_scale(
        color="independent"
    )
    full_fig = alt.vconcat(upper_stack, bottom_row, spacing=18).resolve_scale(color="independent")
    return Success(apply_nature_methods_theme(full_fig))


def build_simplified_bicluster_typst_caption(
    full_sim_df: pl.DataFrame,
    nearest_df: pl.DataFrame,
    coulton_ici_df: pl.DataFrame,
) -> Result[str, str]:
    """Generate a standalone data-synchronized Typst caption module for the heading-free bicluster figure."""
    match _extract_simplified_bicluster_tables(full_sim_df, nearest_df, coulton_ici_df):
        case Failure(err):
            return Failure(err)
        case Success(extracted):
            sub_clean, meta_side_filtered, _top_df, y_domain, x_domain, _rb, _cb = extracted
        case _:
            return Failure("Unexpected Result state in _extract_simplified_bicluster_tables")

    n_rows = len(y_domain)
    n_cols = len(x_domain)

    state_meta_map: dict[str, tuple[str, float, str]] = {}
    for row in meta_side_filtered.iter_rows(named=True):
        state_meta_map[str(row["meta_state_label"])] = (
            str(row["meta_state"]),
            float(row["meta42_bulk_beta"]),
            str(row["meta_star"]),
        )

    row_summaries: list[str] = []
    for r_lbl in y_domain:
        m_id, m_beta, m_star = state_meta_map[r_lbl]
        top1_sub = sub_clean.filter(
            (pl.col("meta_state_label") == r_lbl) & (pl.col("rank_among_true_tams") == 1)
        )
        top_tags = ", ".join(
            f"`{r['coulton_clean_label']}` `#1`"
            for r in top1_sub.iter_rows(named=True)
        )
        match_clause = f" $<->$ {top_tags}" if top_tags else ""
        row_summaries.append(
            f"`{m_id}` ($hat(beta)_(\"meta\") = {m_beta:+.2f}$, `{m_star}`{match_clause})"
        )

    cols_formatted = ", ".join(f"`{c}`" for c in x_domain)
    rows_formatted = "; ".join(row_summaries)

    typst_text = (
        "// Auto-generated by scripts/coulton_tam_comparison/render_figures.py\n"
        "// Import in Typst via:\n"
        "// #import \"../captions/coulton_tam_simplified_bicluster_main_caption.typ\": caption_coulton_simplified_bicluster_main\n\n"
        "#let caption_coulton_simplified_bicluster_main = [\n"
        "  *Ward-biclustered cross-atlas alignment and triangulated checkpoint blockade outcomes between significant `Consensus-Meta42` myeloid states and Coulton et al. pan-cancer TAM substates.*\n"
        f"  *(a)* Optimal-leaf-ordered Ward column dendrogram and pan-macrophage-adjusted meta-analytic logistic log-odds slopes $hat(beta)_(\"ICB\")^(\"adj\")$ across $11$ independent bulk RNA-seq checkpoint blockade trials ($N = 1,058$ patients; adjusted for pan-macrophage lineage load `CD68`, `CD14`, `CSF1R`, `TYROBP`, `C1QA/B/C`) across the ${n_cols}$ true tumour-associated macrophage (`TAM`) substates ({cols_formatted}; numeric Seurat cluster prefixes omitted) from the Coulton et al. pan-cancer atlas ($n = 279,104$ tumor-infiltrating cells across $32$ studies and $22$ solid cancer types) @coultonUsingPancancerAtlas2024 that achieve $q < 0.01$ (`**` or stronger) in `CPI1000+` ($N = 1,446$) `DESeq2`/`fgsea` responder versus non-responder enrichment (`*` $q < 0.05$, `**` $q < 0.01$, `***` $q < 10^(-3)$, `****` $q < 10^(-4)$).\n"
        f"  *(b)* Simultaneously biclustered ${n_rows} times {n_cols}$ within-compartment Marker-Anchored Hybrid proximity matrix ($0.60$ whole-transcriptome `Phi_meta` specificity cosine/Spearman plus $0.40$ canonical top-4 HGNC marker enrichment $tanh(z_(\"DEG\") / 2)$) across all ${n_rows}$ `Consensus-Meta42` myeloid states achieving $p < 0.01$ (`**` or stronger) in bulk deconvolution, with `#1` badges marking each state's closest Coulton `TAM` substate: {rows_formatted}.\n"
        "  *(c)* Paired `Consensus-Meta42` `Pan-iAtlas-Meta` bulk deconvolution log-odds $hat(beta)_(\"meta\")$ ($plus.minus 95%$ CI, $N = 1,015$ patients; `*` $p < 0.05$, `**` $p < 0.01$, `***` $p < 10^(-3)$, `****` $p < 10^(-4)$) aligned with the right-hand Ward row dendrogram.\n"
        "]\n"
    )
    return Success(typst_text)
