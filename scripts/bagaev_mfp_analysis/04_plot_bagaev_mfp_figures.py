#!/usr/bin/env python3
"""Step 04: Pure Plotting Script for Bagaev MFP & Zero-Leakage Deconvolution Vector Figures.

Strictly loads pre-computed Parquet artifacts from Steps 01–03 (zero inline statistics or scipy calls),
applies `plotting_utils.apply_nature_methods_theme` + Okabe-Ito palettes, and exports Inkscape-layered
SVGs + 300 DPI PNGs via `export_altair_inkscape_svg`.
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Final
import altair as alt
import polars as pl
from returns.result import Failure

from plotting_utils import (
    COLOR_BG_WHITE,
    COLOR_DIVIDER_RULE,
    COLOR_HAIRLINE,
    COLOR_SUBTLE_FILL,
    COLOR_TEXT_MUTED,
    COLOR_TEXT_PRIMARY,
    OKABE_BLUE,
    OKABE_BLUISH_GREEN,
    OKABE_ORANGE,
    OKABE_REDDISH_PURPLE,
    OKABE_SKY_BLUE,
    OKABE_VERMILION,
    apply_nature_methods_theme,
    export_altair_inkscape_svg,
)
from tme_datasets import BAGAEV_MFP_SIGNATURE_ORDER, get_data_paths
from tme_response import CANONICAL_SUBTYPE_ORDER

SUBTYPE_COLORS: Final[dict[str, str]] = {
    "IE": OKABE_BLUISH_GREEN,
    "IE/F": OKABE_BLUE,
    "F": OKABE_ORANGE,
    "D": OKABE_VERMILION,
}
DECONV_CLUSTER_PALETTE: Final[list[str]] = [
    OKABE_BLUISH_GREEN,
    OKABE_BLUE,
    OKABE_ORANGE,
    OKABE_VERMILION,
    OKABE_SKY_BLUE,
    OKABE_REDDISH_PURPLE,
]
SUBTYPE_DOMAIN: Final[list[str]] = list(CANONICAL_SUBTYPE_ORDER)
SUBTYPE_RANGE: Final[list[str]] = [SUBTYPE_COLORS[k] for k in SUBTYPE_DOMAIN]


def _clean_records(df: pl.DataFrame) -> list[dict[str, object]]:
    """Replace non-finite floats (NaN / Inf) with JSON-compliant nulls before Vega-Lite export."""
    float_cols = [c for c, dt in df.schema.items() if dt in (pl.Float32, pl.Float64)]
    cleaned = (
        df.with_columns([pl.when(pl.col(c).is_finite()).then(pl.col(c)).otherwise(None).alias(c) for c in float_cols])
        if float_cols
        else df
    )
    return cleaned.to_dicts()


def _alt_data(records: Sequence[Mapping[str, object]]) -> alt.Data:
    """Wrap inline records in Altair Data safely under strict mypy."""
    return alt.Data(values=records)  # type: ignore[no-untyped-call]


def _export_layered_figure(chart: alt.TopLevelMixin, out_svg: Path) -> None:
    """Apply Nature Methods theme and export Inkscape-layered SVG + 300 DPI PNG."""
    themed = apply_nature_methods_theme(chart)
    res = export_altair_inkscape_svg(themed, out_svg, scale=2.5, max_vector_points=1500)
    if isinstance(res, Failure):
        raise RuntimeError(f"Failed exporting {out_svg.name}: {res.failure()}")
    print(f"Exported Inkscape-layered {out_svg.name} (+ .png)")


def plot_fig1_iatlas_portrait_and_clustering(data_dir: Path, results_dir: Path) -> None:
    """Fig 1: 29-Fges Subtype Mean Portrait Heatmap + Immune-vs-Fibrotic Quadrant & UMAP."""
    samples_df = pl.read_parquet(data_dir / "iatlas_bagaev_sample_subtypes.parquet")
    means_df = pl.read_parquet(data_dir / "iatlas_subtype_29fges_means.parquet")
    sig_cols = [c for c in BAGAEV_MFP_SIGNATURE_ORDER if c in samples_df.columns]
    n_total = samples_df.height

    panel_a = (
        alt.Chart(_alt_data(_clean_records(means_df)))
        .mark_rect(stroke=COLOR_BG_WHITE, strokeWidth=0.6)
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN, title="Bagaev TME Subtype (TCGA k-NN)"),
            y=alt.Y("signature:N", sort=sig_cols, title="29 Functional Gene Expression Signatures (Fges)"),
            color=alt.Color(
                "mean_z:Q",
                scale=alt.Scale(domain=[-1.6, 0.0, 1.6], range=[OKABE_BLUE, COLOR_SUBTLE_FILL, OKABE_VERMILION], clamp=True),
                title="Mean Scaled Fges",
            ),
        )
        .properties(width=200, height=430, title=f"A. iAtlas 29-Fges Subtype Profiles (N = {n_total:,})")
    )

    quad_records = _clean_records(
        samples_df.select(
            [
                "sample_id",
                "cohort",
                "bagaev_subtype",
                "louvain_subtype",
                "antitumor_immune_score",
                "stromal_fibrotic_score",
                "umap_1",
                "umap_2",
            ]
        )
    )
    zero_h = alt.Chart(_alt_data([{"y": 0.0}])).mark_rule(stroke=COLOR_TEXT_MUTED, strokeDash=[4, 4], strokeWidth=0.9).encode(y="y:Q")
    zero_v = alt.Chart(_alt_data([{"x": 0.0}])).mark_rule(stroke=COLOR_TEXT_MUTED, strokeDash=[4, 4], strokeWidth=0.9).encode(x="x:Q")

    pts_quad = (
        alt.Chart(_alt_data(quad_records))
        .mark_circle(size=26, opacity=0.72)
        .encode(
            x=alt.X("antitumor_immune_score:Q", title="Anti-Tumor Immune Fges Mean (MAD-Scaled)"),
            y=alt.Y("stromal_fibrotic_score:Q", title="Angiogenesis & Fibroblast Fges Mean (MAD-Scaled)"),
            color=alt.Color("bagaev_subtype:N", scale=alt.Scale(domain=SUBTYPE_DOMAIN, range=SUBTYPE_RANGE), title="TME Subtype"),
        )
    )
    panel_b = (pts_quad + zero_h + zero_v).properties(
        width=280, height=200, title="B. Immune vs. Fibrotic Architecture (TCGA k-NN)"
    )

    panel_c = (
        alt.Chart(_alt_data(quad_records))
        .mark_circle(size=26, opacity=0.72)
        .encode(
            x=alt.X("umap_1:Q", title="29-Fges UMAP 1"),
            y=alt.Y("umap_2:Q", title="29-Fges UMAP 2"),
            color=alt.Color("louvain_subtype:N", scale=alt.Scale(domain=SUBTYPE_DOMAIN, range=SUBTYPE_RANGE), title="TME Subtype"),
        )
        .properties(width=280, height=200, title="C. Unsupervised De Novo Louvain Clusters")
    )

    _export_layered_figure(
        alt.hconcat(panel_a, alt.vconcat(panel_b, panel_c, spacing=24), spacing=28),
        results_dir / "fig1_iatlas_29fges_portrait_and_clustering.svg",
    )


def plot_fig2_iatlas_response_rates(data_dir: Path, results_dir: Path) -> None:
    """Fig 2: Patient-Deduplicated Pan-Cancer & Per-Cohort ORR Across Bagaev Subtypes."""
    resp_df = pl.read_parquet(data_dir / "iatlas_subtype_response_rates.parquet").filter(
        pl.col("sample_scope") == "Baseline_Pre"
    )
    pan_sub_df = resp_df.filter(
        (pl.col("cohort") == "PanCancer-iAtlas")
        & (pl.col("method") == "Supervised_TCGA_kNN_Matched")
    )
    n_pan = int(pan_sub_df["n_samples"].sum()) if not pan_sub_df.is_empty() else 0

    pan_sub = _clean_records(
        resp_df.filter(
            (pl.col("cohort") == "PanCancer-iAtlas")
            & pl.col("method").is_in(["Supervised_TCGA_kNN_Matched", "Unsupervised_DeNovo_Louvain"])
        ).with_columns(
            (pl.col("response_rate") * 100.0).alias("orr_pct"),
            (pl.col("ci95_low") * 100.0).alias("ci_lo_pct"),
            (pl.col("ci95_high") * 100.0).alias("ci_hi_pct"),
            pl.when(pl.col("method") == "Supervised_TCGA_kNN_Matched")
            .then(pl.lit("Supervised TCGA k-NN (K=35)"))
            .otherwise(pl.lit("De Novo Louvain Graph"))
            .alias("method_display"),
        )
    )

    bars_a = (
        alt.Chart(_alt_data(pan_sub))
        .mark_bar(size=24)
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN, title="TME Subtype"),
            xOffset=alt.XOffset("method_display:N"),
            y=alt.Y("orr_pct:Q", title="Objective Response Rate (% CR/PR)"),
            color=alt.Color(
                "method_display:N",
                scale=alt.Scale(domain=["Supervised TCGA k-NN (K=35)", "De Novo Louvain Graph"], range=[OKABE_BLUE, OKABE_SKY_BLUE]),
                title="Classification Method",
            ),
        )
    )
    err_a = (
        alt.Chart(_alt_data(pan_sub))
        .mark_errorbar(color=COLOR_TEXT_PRIMARY, thickness=1.2)
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN, title="TME Subtype"),
            xOffset=alt.XOffset("method_display:N"),
            y=alt.Y("ci_lo_pct:Q", title="Objective Response Rate (% CR/PR)"),
            y2=alt.Y2("ci_hi_pct:Q"),
        )
    )
    panel_a = (bars_a + err_a).properties(
        width=280, height=280, title=f"A. Pan-iAtlas Unique Baseline Patients ORR (N = {n_pan:,})"
    )

    cohort_sub = _clean_records(
        resp_df.filter(
            (pl.col("cohort") != "PanCancer-iAtlas")
            & (pl.col("method") == "Supervised_TCGA_kNN_Matched")
            & (pl.col("n_samples") >= 2)
        ).with_columns(
            (pl.col("response_rate") * 100.0).alias("orr_pct"),
            pl.format(
                "{}% ({}/{})",
                (pl.col("response_rate") * 100.0).round(0).cast(pl.Int64),
                pl.col("n_responders"),
                pl.col("n_samples"),
            ).alias("cell_label"),
        )
    )

    hm_b = (
        alt.Chart(_alt_data(cohort_sub))
        .mark_rect(stroke=COLOR_BG_WHITE, strokeWidth=1.0)
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN, title="Bagaev TME Subtype (TCGA k-NN)"),
            y=alt.Y("cohort:N", title="iAtlas Clinical Trial Cohort"),
            color=alt.Color(
                "orr_pct:Q",
                scale=alt.Scale(domain=[0.0, 65.0], range=[COLOR_SUBTLE_FILL, OKABE_BLUISH_GREEN], clamp=True),
                title="ORR (%)",
            ),
        )
    )
    txt_b = (
        alt.Chart(_alt_data(cohort_sub))
        .mark_text(fontSize=9.5, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN), y=alt.Y("cohort:N"), text=alt.Text("cell_label:N"))
    )
    panel_b = (hm_b + txt_b).properties(width=290, height=280, title="B. Per-Trial Patient ORR % (Responders / Total)")

    _export_layered_figure(alt.hconcat(panel_a, panel_b, spacing=32), results_dir / "fig2_iatlas_subtype_response_rates.svg")


def plot_fig3_iatlas_deconvolution_integration(data_dir: Path, results_dir: Path) -> None:
    """Fig 3: Outcome-Blind Meta42 Enrichment Across Subtypes + Zero-Leakage OOF Discrimination."""
    hm_df = pl.read_parquet(data_dir / "iatlas_meta42_enrichment_heatmap_z.parquet").sort("order_idx")
    joined_df = pl.read_parquet(data_dir / "iatlas_bagaev_joined_meta42_oof.parquet")
    within_df = pl.read_parquet(data_dir / "iatlas_within_subtype_deconv_auc.parquet")

    panel_a = (
        alt.Chart(_alt_data(_clean_records(hm_df)))
        .mark_rect(stroke=COLOR_BG_WHITE, strokeWidth=0.8)
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN, title="Bagaev TME Subtype"),
            y=alt.Y("state_label:N", sort=hm_df["state_label"].unique(maintain_order=True).to_list(), title="Top Subtype-Differentiated Meta42 States"),
            color=alt.Color(
                "rel_z:Q",
                scale=alt.Scale(domain=[-1.5, 0.0, 1.5], range=[OKABE_BLUE, COLOR_SUBTLE_FILL, OKABE_VERMILION], clamp=True),
                title="Relative Enrichment (Z)",
            ),
        )
        .properties(width=200, height=320, title=f"A. Raw InstaPrism Meta42 Enrichment (N = {joined_df.height:,})")
    )

    oof_sub = joined_df.filter(
        pl.col("is_primary_baseline_patient")
        & pl.col("response_binary").cast(pl.Float64).is_in([0.0, 1.0])
        & pl.col("oof_raw_meta42_z").is_not_null()
        & pl.col("oof_raw_meta42_z").cast(pl.Float64).is_finite()
    )
    box_recs = _clean_records(
        oof_sub.with_columns(
            pl.when(pl.col("response_binary") == 1.0)
            .then(pl.lit("Responder (CR/PR)"))
            .otherwise(pl.lit("Non-Responder (SD/PD)"))
            .alias("clinical_outcome")
        ).select(["bagaev_subtype", "clinical_outcome", "oof_raw_meta42_z"])
    )

    panel_b = (
        alt.Chart(_alt_data(box_recs))
        .mark_boxplot(extent=1.5, size=20)
        .encode(
            x=alt.X("bagaev_subtype:N", sort=SUBTYPE_DOMAIN, title="Bagaev TME Subtype"),
            xOffset=alt.XOffset("clinical_outcome:N"),
            y=alt.Y("oof_raw_meta42_z:Q", title="Zero-Leakage OOF Meta42 Score (Z)"),
            color=alt.Color(
                "clinical_outcome:N",
                scale=alt.Scale(domain=["Responder (CR/PR)", "Non-Responder (SD/PD)"], range=[OKABE_BLUISH_GREEN, OKABE_VERMILION]),
                title="Clinical Outcome",
            ),
        )
        .properties(width=270, height=150, title=f"B. Patient-Grouped OOF Meta42 (N = {oof_sub.height:,})")
    )

    within_recs = _clean_records(within_df.filter(pl.col("roc_auc").is_finite()))
    rule_50 = alt.Chart(_alt_data([{"x": 0.5}])).mark_rule(stroke=COLOR_TEXT_MUTED, strokeDash=[4, 4]).encode(x="x:Q")
    bars_c = (
        alt.Chart(_alt_data(within_recs))
        .mark_bar(height=18, color=OKABE_BLUE)
        .encode(
            y=alt.Y("stratum:N", sort=["All_Combined", *SUBTYPE_DOMAIN], title="Evaluated Patient Stratum"),
            x=alt.X("roc_auc:Q", scale=alt.Scale(domain=[0.0, 0.80]), title="Grouped-CV Raw Meta42 ROC-AUC"),
        )
    )
    txt_c = (
        alt.Chart(_alt_data(within_recs))
        .mark_text(align="left", dx=4, fontSize=9.5, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(y=alt.Y("stratum:N", sort=["All_Combined", *SUBTYPE_DOMAIN]), x=alt.X("roc_auc:Q"), text=alt.Text("roc_auc:Q", format=".3f"))
    )
    panel_c = (bars_c + rule_50 + txt_c).properties(
        width=270, height=135, title="C. Within-Subtype Zero-Leakage ROC-AUC"
    )

    _export_layered_figure(
        alt.hconcat(panel_a, alt.vconcat(panel_b, panel_c, spacing=26).resolve_scale(color="independent"), spacing=30).resolve_scale(color="independent"),
        results_dir / "fig3_iatlas_deconvolution_integration.svg",
    )


def plot_fig4_sc23_pseudobulk_conservation(data_dir: Path, results_dir: Path) -> None:
    """Fig 4: 29x29 Correlation Conservation vs TCGA & Subtype Proportions Across All 23 SC Cohorts."""
    sum_df = pl.read_parquet(data_dir / "sc23_pseudobulk_cohort_summary.parquet").sort(
        ["sorting_regime", "offdiag_pearson_r_vs_tcga"], descending=[True, True], nulls_last=True
    )
    cohort_order = sum_df["display_name"].to_list()

    bars_a = (
        alt.Chart(
            _alt_data(
                _clean_records(
                    sum_df.filter(pl.col("offdiag_pearson_r_vs_tcga").is_finite()).select(
                        ["display_name", "offdiag_pearson_r_vs_tcga", "sorting_regime"]
                    )
                )
            )
        )
        .mark_bar(height=13)
        .encode(
            y=alt.Y("display_name:N", sort=cohort_order, title="Tumor-Only Single-Cell Pseudobulk Cohort (9 SC9 + 14 Atlases)"),
            x=alt.X("offdiag_pearson_r_vs_tcga:Q", title="29x29 Fges Correlation Conservation vs. TCGA (Pearson r)"),
            color=alt.Color(
                "sorting_regime:N",
                scale=alt.Scale(domain=["Whole_Tumor", "FACS_Immune_Sorted"], range=[OKABE_BLUE, OKABE_ORANGE]),
                title="Cell Selection Strategy",
            ),
        )
    )
    small_n_df = sum_df.filter(~pl.col("offdiag_pearson_r_vs_tcga").is_finite()).with_columns(
        pl.lit(0.03).alias("x_pos"),
        pl.format("n={} (<4 min)", pl.col("n_samples")).alias("note"),
    )
    panel_a_layers = bars_a
    if not small_n_df.is_empty():
        txt_small = (
            alt.Chart(_alt_data(_clean_records(small_n_df.select(["display_name", "x_pos", "note"]))))
            .mark_text(align="left", fontSize=8.5, fontStyle="italic", color=COLOR_TEXT_MUTED)
            .encode(y=alt.Y("display_name:N", sort=cohort_order), x=alt.X("x_pos:Q"), text="note:N")
        )
        panel_a_layers = bars_a + txt_small
    panel_a = panel_a_layers.properties(width=260, height=390, title="A. Conservation of Bagaev 29-Fges vs. TCGA")

    prop_rows = [
        {
            "display_name": str(r["display_name"]),
            "subtype": st_key,
            "pct": 100.0 * float(r[col_k]) / float(max(1, int(r["n_samples"]))),
        }
        for r in sum_df.iter_rows(named=True)
        for st_key, col_k in (("IE", "n_IE"), ("IE/F", "n_IEF"), ("F", "n_F"), ("D", "n_D"))
    ]
    panel_b = (
        alt.Chart(_alt_data(prop_rows))
        .mark_bar(height=13)
        .encode(
            y=alt.Y("display_name:N", sort=cohort_order, axis=alt.Axis(labels=False, title=None)),
            x=alt.X("pct:Q", stack="normalize", axis=alt.Axis(format="%"), title="Tumor Pseudobulk Subtype Proportion"),
            color=alt.Color("subtype:N", scale=alt.Scale(domain=SUBTYPE_DOMAIN, range=SUBTYPE_RANGE), title="TME Subtype"),
        )
        .properties(width=210, height=390, title="B. Tumor Pseudobulk Subtype Distribution")
    )

    _export_layered_figure(
        alt.hconcat(panel_a, panel_b, spacing=18).resolve_scale(color="independent"),
        results_dir / "fig4_sc23_pseudobulk_clustering_conservation.svg",
    )


def plot_fig5_sc9_milopy_response(data_dir: Path, results_dir: Path) -> None:
    """Fig 5: Patient-Deduplicated Milopy SC9 Pseudobulk Subtype Responder Rates & AUCs."""
    resp_path = data_dir / "sc9_milopy_subtype_response_rates.parquet"
    auc_path = data_dir / "sc9_milopy_cohort_bagaev_aucs.parquet"
    if not resp_path.is_file() or not auc_path.is_file():
        return

    resp_raw = pl.read_parquet(resp_path)
    if "sample_scope" in resp_raw.columns:
        resp_raw = resp_raw.filter(pl.col("sample_scope") == "Unique_Patients_PreferPre")
    resp_df = resp_raw.filter(
        pl.col("cohort").is_in(["SC9_WholeTumor_7_Pooled", "SC9_FACS_ImmuneSorted_2_Pooled"])
        & (pl.col("method") == "TCGA_kNN")
    )
    recs_a = _clean_records(
        resp_df.with_columns(
            (pl.col("response_rate") * 100.0).alias("orr_pct"),
            (pl.col("ci95_low") * 100.0).alias("ci_lo_pct"),
            (pl.col("ci95_high") * 100.0).alias("ci_hi_pct"),
            pl.when(pl.col("cohort") == "SC9_WholeTumor_7_Pooled")
            .then(pl.lit("Whole-Tumor SC9 (7 Cohorts)"))
            .otherwise(pl.lit("FACS Immune-Sorted SC9 (2 Cohorts)"))
            .alias("stratum_label"),
        ).select(["subtype", "stratum_label", "orr_pct", "ci_lo_pct", "ci_hi_pct"])
    )
    n_sc9_pat = int(resp_df["n_samples"].sum()) if not resp_df.is_empty() else 0

    bars_a = (
        alt.Chart(_alt_data(recs_a))
        .mark_bar(size=22)
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN, title="Bagaev TME Subtype (Tumor Pseudobulk k-NN)"),
            xOffset=alt.XOffset("stratum_label:N"),
            y=alt.Y("orr_pct:Q", title="ICB Responder Rate (% R)"),
            color=alt.Color(
                "stratum_label:N",
                scale=alt.Scale(domain=["Whole-Tumor SC9 (7 Cohorts)", "FACS Immune-Sorted SC9 (2 Cohorts)"], range=[OKABE_BLUE, OKABE_ORANGE]),
                title="SC9 Compartment Coverage",
            ),
        )
    )
    err_a = (
        alt.Chart(_alt_data(recs_a))
        .mark_errorbar(color=COLOR_TEXT_PRIMARY, thickness=1.1)
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN, title="Bagaev TME Subtype (Tumor Pseudobulk k-NN)"),
            xOffset=alt.XOffset("stratum_label:N"),
            y=alt.Y("ci_lo_pct:Q", title="ICB Responder Rate (% R)"),
            y2=alt.Y2("ci_hi_pct:Q"),
        )
    )
    panel_a = (bars_a + err_a).properties(
        width=280, height=260, title=f"A. Milopy SC9 Patient ORR by Subtype (N = {n_sc9_pat:,})"
    )

    auc_df = pl.read_parquet(auc_path).sort("auc_bagaev_immune_minus_fibrotic", descending=True)
    c_order = auc_df["display_name"].to_list()
    rule_50 = alt.Chart(_alt_data([{"x": 0.5}])).mark_rule(stroke=COLOR_TEXT_MUTED, strokeDash=[4, 4]).encode(x="x:Q")
    pts_b = (
        alt.Chart(_alt_data(_clean_records(auc_df)))
        .mark_bar(height=18)
        .encode(
            y=alt.Y("display_name:N", sort=c_order, title="Milopy Single-Cell ICB Cohort (Unique Patients)"),
            x=alt.X("auc_bagaev_immune_minus_fibrotic:Q", scale=alt.Scale(domain=[0.0, 1.0]), title="Pseudobulk Bagaev Immune - Fibrotic ROC-AUC"),
            color=alt.Color(
                "sorting_regime:N",
                scale=alt.Scale(domain=["Whole_Tumor", "FACS_Immune_Sorted"], range=[OKABE_BLUE, OKABE_ORANGE]),
                title="Sorting Regime",
            ),
        )
    )
    panel_b = (pts_b + rule_50).properties(width=280, height=260, title="B. Per-Cohort Patient ROC-AUC (SC9)")

    _export_layered_figure(
        alt.hconcat(panel_a, panel_b, spacing=30).resolve_scale(color="independent"),
        results_dir / "fig5_sc9_milopy_response_correlation.svg",
    )


def plot_fig6_deconv_subtype_correlations(data_dir: Path, results_dir: Path) -> None:
    """Fig 6: Outcome-Blind Deconvolved Cell State Correlations with Bagaev Subtypes & 29 Fges."""
    corr_path = data_dir / "iatlas_deconv_subtype_correlations.parquet"
    fges_path = data_dir / "iatlas_meta42_vs_29fges_correlations.parquet"
    if not corr_path.is_file():
        return

    corr_df = pl.read_parquet(corr_path).with_columns(
        pl.coalesce([pl.col("state_label"), pl.col("cell_state")]).alias("display_state")
    )
    n_eval = int(corr_df["n_evaluated_samples"][0]) if "n_evaluated_samples" in corr_df.columns else 1015
    meta_r_val = float(corr_df["pooled_vs_meta_pearson_r"][0]) if "pooled_vs_meta_pearson_r" in corr_df.columns else 0.0

    piv_ie = corr_df.filter(pl.col("subtype") == "IE").select(["cell_state", "display_state", pl.col("spearman_rho_raw").alias("rho_ie")])
    piv_f = corr_df.filter(pl.col("subtype") == "F").select(["cell_state", pl.col("spearman_rho_raw").alias("rho_f")])
    state_order_df = (
        piv_ie.join(piv_f, on="cell_state", how="left")
        .with_columns((pl.col("rho_ie") - pl.col("rho_f").fill_null(0.0)).alias("ie_minus_f"))
        .sort("ie_minus_f", descending=True)
    )
    state_order = state_order_df["display_state"].to_list()
    top_bottom_ids = state_order_df["cell_state"].head(9).to_list() + state_order_df["cell_state"].tail(9).to_list()

    hm_a_df = corr_df.with_columns(
        pl.when(pl.col("fdr_q_raw") < 0.001)
        .then(pl.lit("***"))
        .when(pl.col("fdr_q_raw") < 0.01)
        .then(pl.lit("**"))
        .when(pl.col("fdr_q_raw") < 0.05)
        .then(pl.lit("*"))
        .otherwise(pl.lit(""))
        .alias("fdr_stars"),
        pl.col("spearman_rho_raw").abs().alias("abs_rho"),
    ).with_columns(
        pl.format("{}{}", pl.col("spearman_rho_raw").round(2), pl.col("fdr_stars")).alias("cell_txt"),
        (pl.col("abs_rho") >= 0.34).alias("use_white_txt"),
    )
    recs_a = _clean_records(hm_a_df.select(["display_state", "subtype", "spearman_rho_raw", "cell_txt", "use_white_txt"]))

    rect_a = (
        alt.Chart(_alt_data(recs_a))
        .mark_rect(stroke=COLOR_BG_WHITE, strokeWidth=0.6)
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN, title="Bagaev TME Subtype (TCGA k-NN)"),
            y=alt.Y("display_state:N", sort=state_order, title="Consensus-Meta42 Deconvolved Cell States (42)"),
            color=alt.Color(
                "spearman_rho_raw:Q",
                scale=alt.Scale(domain=[-0.55, 0.0, 0.55], range=[OKABE_BLUE, COLOR_SUBTLE_FILL, OKABE_VERMILION], clamp=True),
                title="Spearman ρ (Fraction vs Subtype)",
            ),
        )
    )
    txt_a_dark = (
        alt.Chart(_alt_data([r for r in recs_a if not r["use_white_txt"]]))
        .mark_text(fontSize=8.5, color=COLOR_TEXT_PRIMARY)
        .encode(x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN), y=alt.Y("display_state:N", sort=state_order), text="cell_txt:N")
    )
    txt_a_light = (
        alt.Chart(_alt_data([r for r in recs_a if r["use_white_txt"]]))
        .mark_text(fontSize=8.5, fontWeight="bold", color=COLOR_BG_WHITE)
        .encode(x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN), y=alt.Y("display_state:N", sort=state_order), text="cell_txt:N")
    )
    panel_a = (rect_a + txt_a_dark + txt_a_light).properties(
        width=220, height=580, title=f"A. Raw InstaPrism Meta42 vs. Subtype Correlations (N = {n_eval:,})"
    )

    b_raw = corr_df.select(["cell_state", "display_state", "subtype", "spearman_rho_raw", "within_cohort_meta_rho"])
    picked_keys: set[tuple[str, str]] = set()
    min_row = b_raw.sort("spearman_rho_raw", descending=False).row(0, named=True)
    picked_keys.add((str(min_row["cell_state"]), str(min_row["subtype"])))
    for st_name in SUBTYPE_DOMAIN:
        sub_st = b_raw.filter(pl.col("subtype") == st_name).sort("spearman_rho_raw", descending=True)
        if not sub_st.is_empty():
            r_top = sub_st.row(0, named=True)
            picked_keys.add((str(r_top["cell_state"]), str(r_top["subtype"])))

    ordered_landmarks = sorted(
        (r for r in b_raw.iter_rows(named=True) if (str(r["cell_state"]), str(r["subtype"])) in picked_keys),
        key=lambda r: float(r["spearman_rho_raw"]),
    )
    kept_landmarks: dict[tuple[str, str], str] = {}
    last_rho = -99.0
    for idx_l, row_l in enumerate(ordered_landmarks):
        rho_l = float(row_l["spearman_rho_raw"])
        if rho_l - last_rho >= 0.08 or rho_l > 0.55 or rho_l < -0.55:
            key_l = (str(row_l["cell_state"]), str(row_l["subtype"]))
            side = "right_of_dot" if (rho_l < -0.50 or (0.40 <= rho_l <= 0.54) or (idx_l % 2 == 0 and rho_l < 0.54)) else "left_of_dot"
            kept_landmarks[key_l] = side
            last_rho = rho_l

    recs_b: list[dict[str, object]] = []
    for rec in _clean_records(b_raw):
        k_pair = (str(rec["cell_state"]), str(rec["subtype"]))
        rec_copy = dict(rec)
        rec_copy["lbl"] = f"{k_pair[0]}:{k_pair[1]}" if k_pair in kept_landmarks else ""
        rec_copy["lbl_side"] = kept_landmarks.get(k_pair, "")
        recs_b.append(rec_copy)

    diag_rule = alt.Chart(_alt_data([{"x": -0.65, "y": -0.65}, {"x": 0.65, "y": 0.65}])).mark_line(stroke=COLOR_TEXT_MUTED, strokeDash=[4, 4], strokeWidth=1.0).encode(x="x:Q", y="y:Q")
    zero_x = alt.Chart(_alt_data([{"x": 0.0}])).mark_rule(stroke=COLOR_HAIRLINE, strokeWidth=0.8).encode(x="x:Q")
    zero_y = alt.Chart(_alt_data([{"y": 0.0}])).mark_rule(stroke=COLOR_HAIRLINE, strokeWidth=0.8).encode(y="y:Q")
    pts_b = (
        alt.Chart(_alt_data(recs_b))
        .mark_circle(size=42, opacity=0.82)
        .encode(
            x=alt.X("spearman_rho_raw:Q", scale=alt.Scale(domain=[-0.68, 0.68]), title="Pooled Pan-iAtlas Spearman ρ"),
            y=alt.Y("within_cohort_meta_rho:Q", scale=alt.Scale(domain=[-0.68, 0.68]), title="Within-Cohort Weighted Fisher-z Meta ρ"),
            color=alt.Color("subtype:N", scale=alt.Scale(domain=SUBTYPE_DOMAIN, range=SUBTYPE_RANGE), title="Bagaev Subtype"),
        )
    )
    txt_b_left = (
        alt.Chart(_alt_data([r for r in recs_b if r["lbl_side"] == "left_of_dot"]))
        .mark_text(align="right", dx=-6, dy=-3, fontSize=8.5, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(x="spearman_rho_raw:Q", y="within_cohort_meta_rho:Q", text="lbl:N")
    )
    txt_b_right = (
        alt.Chart(_alt_data([r for r in recs_b if r["lbl_side"] == "right_of_dot"]))
        .mark_text(align="left", dx=6, dy=5, fontSize=8.5, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
        .encode(x="spearman_rho_raw:Q", y="within_cohort_meta_rho:Q", text="lbl:N")
    )
    panel_b = (diag_rule + zero_x + zero_y + pts_b + txt_b_left + txt_b_right).properties(
        width=280, height=240, title=f"B. Within-Cohort Conservation (r = {meta_r_val:.3f}, 168 pairs)"
    )

    if fges_path.is_file():
        core_fges = [
            "Effector_cells", "T_cells", "Th1_signature", "Checkpoint_inhibition",
            "B_cells", "NK_cells", "M1_signatures", "MHCII",
            "Treg", "Macrophages", "MDSC",
            "Angiogenesis", "Endothelium", "CAF", "Matrix_remodeling", "Proliferation_rate",
        ]
        fges_df = (
            pl.read_parquet(fges_path)
            .filter(pl.col("cell_state").is_in(top_bottom_ids) & pl.col("fges_signature").is_in(core_fges))
            .with_columns(pl.coalesce([pl.col("state_label"), pl.col("cell_state")]).alias("display_state"))
        )
        sub_state_order = [s for s in state_order if any(s == ds for ds in fges_df["display_state"].unique().to_list())]
        recs_c = _clean_records(fges_df.select(["display_state", "fges_signature", "spearman_rho"]))
        panel_c = (
            alt.Chart(_alt_data(recs_c))
            .mark_rect(stroke=COLOR_BG_WHITE, strokeWidth=0.6)
            .encode(
                x=alt.X("fges_signature:N", sort=core_fges, axis=alt.Axis(labelAngle=-45), title="Core Bagaev Fges"),
                y=alt.Y("display_state:N", sort=sub_state_order, title="Top Subtype-Associated States"),
                color=alt.Color(
                    "spearman_rho:Q",
                    scale=alt.Scale(domain=[-0.70, 0.0, 0.70], range=[OKABE_BLUE, COLOR_SUBTLE_FILL, OKABE_VERMILION], clamp=True),
                    title="Spearman ρ (State vs Fges)",
                ),
            )
            .properties(width=280, height=290, title="C. Cell State vs. Bagaev 29-Fges Mapping")
        )
        right_col = alt.vconcat(panel_b, panel_c, spacing=30).resolve_scale(color="independent")
    else:
        right_col = panel_b

    _export_layered_figure(
        alt.hconcat(panel_a, right_col, spacing=34).resolve_scale(color="independent"),
        results_dir / "fig6_deconv_subtype_correlations.svg",
    )


def plot_fig7_deconv_sample_clusters_and_response(data_dir: Path, results_dir: Path) -> None:
    """Fig 7: Unsupervised Zero-Leakage Deconvolution Clusters, Bagaev Overlap & Patient ORR."""
    samp_path = data_dir / "iatlas_deconv_cluster_samples.parquet"
    cent_path = data_dir / "iatlas_deconv_cluster_centroids_k4.parquet"
    resp_path = data_dir / "iatlas_deconv_cluster_response_rates.parquet"
    cont_path = data_dir / "iatlas_deconv_clusters_vs_bagaev_contingency.parquet"
    bench_path = data_dir / "iatlas_multires_deconv_cluster_benchmark.parquet"
    if not (samp_path.is_file() and cent_path.is_file() and resp_path.is_file()):
        return

    samp_df = pl.read_parquet(samp_path)
    cent_df = pl.read_parquet(cent_path)
    resp_df = pl.read_parquet(resp_path)

    cluster_ids = sorted(samp_df["deconv_cluster_pan_ward_k4"].unique().to_list())
    label_order = [
        str(samp_df.filter(pl.col("deconv_cluster_pan_ward_k4") == cid)["deconv_cluster_pan_ward_k4_label"][0])
        for cid in cluster_ids
    ]
    cl_colors = DECONV_CLUSTER_PALETTE[: len(cluster_ids)]

    recs_a = _clean_records(samp_df.select(["clr_umap_1", "clr_umap_2", "deconv_cluster_pan_ward_k4", "deconv_cluster_pan_ward_k4_label"]))
    panel_a = (
        alt.Chart(_alt_data(recs_a))
        .mark_circle(size=24, opacity=0.76)
        .encode(
            x=alt.X("clr_umap_1:Q", title="Raw Meta42 CLR UMAP 1"),
            y=alt.Y("clr_umap_2:Q", title="Raw Meta42 CLR UMAP 2"),
            color=alt.Color("deconv_cluster_pan_ward_k4_label:N", scale=alt.Scale(domain=label_order, range=cl_colors), title="Deconv Cluster (CLR-Ward K=4)"),
        )
        .properties(width=250, height=230, title=f"A. Zero-Leakage Tumor Ecosystem Clusters (N = {samp_df.height:,})")
    )

    st_spread = (
        cent_df.group_by("cell_state")
        .agg([(pl.col("mean_clr_z").max() - pl.col("mean_clr_z").min()).alias("z_range")])
        .sort("z_range", descending=True)
        .head(16)
    )
    cent_sub = cent_df.filter(pl.col("cell_state").is_in(st_spread["cell_state"].to_list()))
    y_cent_order = (
        cent_sub.filter(pl.col("cluster_id") == "DC1")
        .select(["state_label", "mean_clr_z"])
        .sort("mean_clr_z", descending=True)["state_label"]
        .to_list()
    )
    panel_b = (
        alt.Chart(_alt_data(_clean_records(cent_sub)))
        .mark_rect(stroke=COLOR_BG_WHITE, strokeWidth=0.6)
        .encode(
            x=alt.X("cluster_id:N", sort=cluster_ids, title="Deconv Cluster"),
            y=alt.Y("state_label:N", sort=y_cent_order, title="Top Cluster-Defining Meta42 States"),
            color=alt.Color(
                "mean_clr_z:Q",
                scale=alt.Scale(domain=[-1.2, 0.0, 1.2], range=[OKABE_BLUE, COLOR_SUBTLE_FILL, OKABE_VERMILION], clamp=True),
                title="Mean CLR Z",
            ),
        )
        .properties(width=185, height=230, title="B. Cluster Cellular Archetypes")
    )

    cont_df = pl.read_parquet(cont_path)
    panel_c = (
        alt.Chart(_alt_data(_clean_records(cont_df)))
        .mark_bar(size=28)
        .encode(
            x=alt.X("deconv_cluster_pan_ward_k4:N", sort=cluster_ids, title="Deconv Cluster (CLR-Ward K=4)"),
            y=alt.Y("n_samples:Q", stack="normalize", axis=alt.Axis(format="%"), title="Bagaev Subtype Composition"),
            color=alt.Color("bagaev_subtype:N", scale=alt.Scale(domain=SUBTYPE_DOMAIN, range=SUBTYPE_RANGE), title="Bagaev Subtype"),
        )
        .properties(width=185, height=230, title="C. Bagaev TME Subtype Overlap")
    )

    d_sub = (
        resp_df.filter(
            (pl.col("cohort") == "PanCancer-iAtlas")
            & (pl.col("sample_scope") == "Baseline_Pre")
            & pl.col("method").is_in(["Meta42_PanCancer_CLR_Ward_K4", "Meta42_CohortCentered_CLR_Ward_K4"])
        )
        .with_columns(
            (pl.col("response_rate") * 100.0).alias("orr_pct"),
            (pl.col("ci95_low") * 100.0).alias("ci_lo_pct"),
            (pl.col("ci95_high") * 100.0).alias("ci_hi_pct"),
            pl.when(pl.col("method") == "Meta42_PanCancer_CLR_Ward_K4")
            .then(pl.lit("Pan-Cancer CLR Ward (K=4)"))
            .otherwise(pl.lit("Cohort-Centered CLR Ward (K=4)"))
            .alias("clustering_regime"),
        )
        .with_columns(
            pl.format("{} ({}/{})", pl.col("orr_pct").round(0).cast(pl.Int64), pl.col("n_responders"), pl.col("n_samples")).alias("orr_lbl")
        )
    )
    recs_d = _clean_records(d_sub)
    bars_d = (
        alt.Chart(_alt_data(recs_d))
        .mark_bar(size=22)
        .encode(
            x=alt.X("subtype:N", sort=cluster_ids, title="Deconvolution Sample Cluster"),
            xOffset=alt.XOffset("clustering_regime:N"),
            y=alt.Y("orr_pct:Q", scale=alt.Scale(domain=[0, 62]), title="Baseline Patient ORR (% CR/PR)"),
            color=alt.Color(
                "clustering_regime:N",
                scale=alt.Scale(domain=["Pan-Cancer CLR Ward (K=4)", "Cohort-Centered CLR Ward (K=4)"], range=[OKABE_BLUE, OKABE_SKY_BLUE]),
                title="CLR Clustering Regime",
            ),
        )
    )
    err_d = (
        alt.Chart(_alt_data(recs_d))
        .mark_errorbar(color=COLOR_TEXT_PRIMARY, thickness=1.2)
        .encode(
            x=alt.X("subtype:N", sort=cluster_ids, title="Deconvolution Sample Cluster"),
            xOffset=alt.XOffset("clustering_regime:N"),
            y=alt.Y("ci_lo_pct:Q", title="Baseline Patient ORR (% CR/PR)"),
            y2=alt.Y2("ci_hi_pct:Q"),
        )
    )
    pan_cmh_p = float(d_sub.filter(pl.col("method") == "Meta42_PanCancer_CLR_Ward_K4")["cmh_omnibus_p"][0])
    coh_cmh_p = float(d_sub.filter(pl.col("method") == "Meta42_CohortCentered_CLR_Ward_K4")["cmh_omnibus_p"][0])
    pan_p_str = f"{pan_cmh_p:.2e}" if pan_cmh_p < 1e-3 else f"{pan_cmh_p:.3f}"
    coh_p_str = f"{coh_cmh_p:.2e}" if coh_cmh_p < 1e-3 else f"{coh_cmh_p:.3f}"
    panel_d = (bars_d + err_d).properties(
        width=310, height=230, title=f"D. Patient Baseline ORR (CMH p = {pan_p_str} / {coh_p_str})"
    )

    if bench_path.is_file():
        b_raw = pl.read_parquet(bench_path).filter(
            (pl.col("clustering_regime") == "Bulk_29Fges_kNN")
            | pl.col("reference_id").is_in(
                [
                    "Consensus-Meta42_RawInstaPrism",
                    "SC9_all_9_combined_res0.5",
                    "SC9_all_9_combined_res1.0",
                    "SC9_solid_tumor_quartet_res0.5",
                    "SC9_solid_tumor_quartet_res1.0",
                ]
            )
        )
        b_df = (
            b_raw.with_columns(
                pl.when(pl.col("clustering_regime") == "Bulk_29Fges_kNN")
                .then(pl.col("display_name"))
                .when(pl.col("clustering_regime") == "PanCancer_CLR_Ward")
                .then(pl.col("display_name") + pl.lit(" [Pan-CLR]"))
                .otherwise(pl.col("display_name") + pl.lit(" [Coh-CLR]"))
                .alias("bench_label"),
                (pl.col("baseline_orr_spread") * 100.0).alias("spread_pct"),
                pl.when(pl.col("baseline_cmh_omnibus_p") < 0.001)
                .then(pl.lit("CMH p < 0.001"))
                .when(pl.col("baseline_cmh_omnibus_p") < 0.05)
                .then(pl.lit("CMH p < 0.05"))
                .otherwise(pl.lit("CMH p >= 0.05"))
                .alias("sig_tier"),
            )
            .sort("spread_pct", descending=True)
        )
        b_order = b_df["bench_label"].to_list()
        bg_spread = b_df.filter(pl.col("clustering_regime") == "Bulk_29Fges_kNN")["spread_pct"].to_list()
        bg_val = float(bg_spread[0]) if bg_spread else 15.0
        recs_e = _clean_records(b_df)
        bars_e = (
            alt.Chart(_alt_data(recs_e))
            .mark_bar(height=15)
            .encode(
                y=alt.Y("bench_label:N", sort=b_order, title="Zero-Leakage Deconv Reference (K = 4)"),
                x=alt.X("spread_pct:Q", title="Baseline Patient ORR Spread (Max - Min Cluster, %)"),
                color=alt.Color(
                    "sig_tier:N",
                    scale=alt.Scale(domain=["CMH p < 0.001", "CMH p < 0.05", "CMH p >= 0.05"], range=[OKABE_BLUISH_GREEN, OKABE_BLUE, COLOR_DIVIDER_RULE]),
                    title="Cohort-Stratified CMH p",
                ),
            )
        )
        txt_e = (
            alt.Chart(_alt_data(recs_e))
            .mark_text(align="left", dx=4, fontSize=8.5, fontWeight="bold", color=COLOR_TEXT_PRIMARY)
            .encode(y=alt.Y("bench_label:N", sort=b_order), x=alt.X("spread_pct:Q"), text=alt.Text("spread_pct:Q", format=".1f"))
        )
        rule_bg = alt.Chart(_alt_data([{"x": bg_val}])).mark_rule(stroke=OKABE_VERMILION, strokeDash=[4, 4], strokeWidth=1.2).encode(x="x:Q")
        panel_e = (bars_e + rule_bg + txt_e).properties(
            width=310, height=230, title="E. Cluster ORR Spread Across Outcome-Blind References"
        )
        bot_row: alt.LayerChart | alt.HConcatChart = alt.hconcat(panel_d, panel_e, spacing=36).resolve_scale(color="independent")
    else:
        bot_row = panel_d

    top_row = alt.hconcat(panel_a, panel_b, panel_c, spacing=28).resolve_scale(color="independent")
    _export_layered_figure(
        alt.vconcat(top_row, bot_row, spacing=32).resolve_scale(color="independent"),
        results_dir / "fig7_deconv_sample_clusters_and_response.svg",
    )


def run_step04(data_dir: Path, results_dir: Path) -> None:
    """Render all 7 Bagaev MFP & Deconvolution analysis Inkscape-layered SVGs."""
    results_dir.mkdir(parents=True, exist_ok=True)
    plot_fig1_iatlas_portrait_and_clustering(data_dir, results_dir)
    plot_fig2_iatlas_response_rates(data_dir, results_dir)
    plot_fig3_iatlas_deconvolution_integration(data_dir, results_dir)
    plot_fig4_sc23_pseudobulk_conservation(data_dir, results_dir)
    plot_fig5_sc9_milopy_response(data_dir, results_dir)
    plot_fig6_deconv_subtype_correlations(data_dir, results_dir)
    plot_fig7_deconv_sample_clusters_and_response(data_dir, results_dir)
    print(f"Step 04 complete. Exported Inkscape SVGs & PNGs to {results_dir}")


def main() -> None:
    paths = get_data_paths()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=paths.output_dir / "bagaev_mfp_analysis")
    parser.add_argument("--results-dir", type=Path, default=paths.results_dir / "bagaev_mfp_analysis")
    args = parser.parse_args()
    run_step04(data_dir=args.data_dir, results_dir=args.results_dir)


if __name__ == "__main__":
    main()
