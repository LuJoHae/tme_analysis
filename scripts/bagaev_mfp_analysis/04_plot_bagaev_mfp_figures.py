#!/usr/bin/env python3
"""Step 04: Publication-Grade Altair Vector Figures (SVG) for Bagaev MFP TME Subtyping Analysis.

Loads pre-computed Parquet artifacts from `output/bagaev_mfp_analysis/` and exports crisp
Nature Methods minimalist vector SVG figures (with Okabe-Ito colorblind-safe palette) to
`results/bagaev_mfp_analysis/`.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Final
import altair as alt
import polars as pl
import vl_convert as vlc

from tme_datasets import BAGAEV_MFP_SIGNATURE_ORDER, get_data_paths
from tme_response import CANONICAL_SUBTYPE_ORDER

# Okabe-Ito Colorblind-Safe Subtype & Status Palette
SUBTYPE_COLORS: Final[dict[str, str]] = {
    "IE": "#009E73",    # Bluish Green (Immune-Enriched, Non-Fibrotic)
    "IE/F": "#0072B2",  # Blue (Immune-Enriched, Fibrotic)
    "F": "#E69F00",     # Orange (Fibrotic / Stromal-Excluded)
    "D": "#D55E00",     # Vermilion (Immune-Depleted / Desert)
}
SUBTYPE_DOMAIN: Final[list[str]] = list(CANONICAL_SUBTYPE_ORDER)
SUBTYPE_RANGE: Final[list[str]] = [SUBTYPE_COLORS[k] for k in SUBTYPE_DOMAIN]


def _clean_records(df: pl.DataFrame) -> list[dict[str, object]]:
    """Replace non-finite floats (NaN / Inf) with JSON-compliant nulls prior to Vega-Lite serialization."""
    float_cols = [c for c, dt in df.schema.items() if dt in (pl.Float32, pl.Float64)]
    cleaned = (
        df.with_columns([pl.when(pl.col(c).is_finite()).then(pl.col(c)).otherwise(None).alias(c) for c in float_cols])
        if float_cols
        else df
    )
    return cleaned.to_dicts()


def _save_svg_and_png(chart: alt.TopLevelMixin, out_svg: Path) -> None:
    """Export Altair chart to publication vector SVG and 300-DPI PNG via vl-convert."""
    spec_json = chart.to_json()
    svg_str = vlc.vegalite_to_svg(spec_json)
    out_svg.write_text(svg_str, encoding="utf-8")
    png_bytes = vlc.vegalite_to_png(spec_json, scale=2.5)
    out_svg.with_suffix(".png").write_bytes(png_bytes)
    print(f"Exported {out_svg.name} (+ .png)")


def _configure_nature_chart(chart: alt.TopLevelMixin) -> alt.TopLevelMixin:
    """Apply Nature Methods minimalist wireframe styling (white canvas, clean sans-serif, hairline axes)."""
    return (
        chart.configure(background="#FFFFFF", font="Helvetica Neue, Arial, sans-serif")
        .configure_axis(
            labelColor="#1E293B",
            titleColor="#0F172A",
            labelFontSize=10,
            titleFontSize=11,
            titleFontWeight="bold",
            domainColor="#1E293B",
            domainWidth=0.8,
            tickColor="#1E293B",
            tickWidth=0.8,
            gridColor="#E2E8F0",
            gridWidth=0.5,
        )
        .configure_header(
            labelColor="#0F172A",
            titleColor="#0F172A",
            labelFontSize=11,
            labelFontWeight="bold",
        )
        .configure_legend(
            labelColor="#1E293B",
            titleColor="#0F172A",
            labelFontSize=10,
            titleFontSize=10.5,
            titleFontWeight="bold",
            labelLimit=260,
            titleLimit=260,
        )
        .configure_title(
            color="#0F172A",
            fontSize=12.5,
            fontWeight="bold",
            anchor="start",
        )
        .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
    )


def plot_fig1_iatlas_portrait_and_clustering(data_dir: Path, results_dir: Path) -> None:
    """Fig 1: 29-Fges Subtype Mean Portrait Heatmap + Immune-vs-Fibrotic Quadrant & UMAP."""
    samples_df = pl.read_parquet(data_dir / "iatlas_bagaev_sample_subtypes.parquet")
    sig_cols = [c for c in BAGAEV_MFP_SIGNATURE_ORDER if c in samples_df.columns]

    # Panel A: Mean scaled 29-Fges profile per Bagaev subtype (4 columns x 29 signatures)
    mean_rows: list[dict[str, str | float | int]] = []
    for s_idx, sig in enumerate(sig_cols):
        for st_name in SUBTYPE_DOMAIN:
            sub = samples_df.filter(pl.col("bagaev_subtype") == st_name)
            val = float(sub[sig].mean()) if not sub.is_empty() else 0.0
            mean_rows.append(
                {
                    "signature": sig,
                    "sig_order": s_idx,
                    "subtype": st_name,
                    "mean_z": val,
                }
            )
    hm_records = mean_rows

    panel_a = (
        alt.Chart(alt.Data(values=hm_records))
        .mark_rect(stroke="#FFFFFF", strokeWidth=0.6)
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN, title="Bagaev TME Subtype (TCGA k-NN)"),
            y=alt.Y("signature:N", sort=sig_cols, title="29 Functional Gene Expression Signatures (Fges)"),
            color=alt.Color(
                "mean_z:Q",
                scale=alt.Scale(domain=[-1.6, 0.0, 1.6], range=["#0072B2", "#F8FAFC", "#D55E00"], clamp=True),
                title="Mean Scaled Fges",
            ),
        )
        .properties(width=200, height=430, title="A. iAtlas 29-Fges Subtype Profiles (N = 1,097)")
    )

    # Panel B: Anti-tumor Immune vs Stromal/Fibrotic Quadrant Colored by Supervised Subtype
    quad_records = samples_df.select(
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
    ).to_dicts()

    zero_h = alt.Chart(alt.Data(values=[{"y": 0.0}])).mark_rule(stroke="#64748B", strokeDash=[4, 4], strokeWidth=0.9).encode(y="y:Q")
    zero_v = alt.Chart(alt.Data(values=[{"x": 0.0}])).mark_rule(stroke="#64748B", strokeDash=[4, 4], strokeWidth=0.9).encode(x="x:Q")

    pts_quad = (
        alt.Chart(alt.Data(values=quad_records))
        .mark_circle(size=26, opacity=0.72)
        .encode(
            x=alt.X("antitumor_immune_score:Q", title="Anti-Tumor Immune Fges Mean (MAD-Scaled)"),
            y=alt.Y("stromal_fibrotic_score:Q", title="Angiogenesis & Fibroblast Fges Mean (MAD-Scaled)"),
            color=alt.Color(
                "bagaev_subtype:N",
                scale=alt.Scale(domain=SUBTYPE_DOMAIN, range=SUBTYPE_RANGE),
                title="TME Subtype",
            ),
        )
    )
    panel_b = (pts_quad + zero_h + zero_v).properties(
        width=280, height=200, title="B. Immune vs. Fibrotic Architecture (TCGA k-NN)"
    )

    panel_c = (
        alt.Chart(alt.Data(values=quad_records))
        .mark_circle(size=26, opacity=0.72)
        .encode(
            x=alt.X("umap_1:Q", title="29-Fges UMAP 1"),
            y=alt.Y("umap_2:Q", title="29-Fges UMAP 2"),
            color=alt.Color(
                "louvain_subtype:N",
                scale=alt.Scale(domain=SUBTYPE_DOMAIN, range=SUBTYPE_RANGE),
                title="TME Subtype",
            ),
        )
        .properties(width=280, height=200, title="C. Unsupervised De Novo Louvain Clusters")
    )

    combined = _configure_nature_chart(alt.hconcat(panel_a, alt.vconcat(panel_b, panel_c, spacing=24), spacing=28))
    _save_svg_and_png(combined, results_dir / "fig1_iatlas_29fges_portrait_and_clustering.svg")


def plot_fig2_iatlas_response_rates(data_dir: Path, results_dir: Path) -> None:
    """Fig 2: Pan-Cancer & Per-Cohort Objective Response Rates Across Bagaev Subtypes."""
    resp_df = pl.read_parquet(data_dir / "iatlas_subtype_response_rates.parquet").filter(
        pl.col("sample_scope") == "Baseline_Pre"
    )

    pan_sub = (
        resp_df.filter(
            (pl.col("cohort") == "PanCancer-iAtlas")
            & pl.col("method").is_in(["Supervised_TCGA_kNN_Matched", "Unsupervised_DeNovo_Louvain"])
        )
        .with_columns(
            (pl.col("response_rate") * 100.0).alias("orr_pct"),
            (pl.col("ci95_low") * 100.0).alias("ci_lo_pct"),
            (pl.col("ci95_high") * 100.0).alias("ci_hi_pct"),
            pl.when(pl.col("method") == "Supervised_TCGA_kNN_Matched")
            .then(pl.lit("Supervised TCGA k-NN (K=35)"))
            .otherwise(pl.lit("De Novo Louvain Graph"))
            .alias("method_display"),
        )
        .to_dicts()
    )

    bars_a = (
        alt.Chart(alt.Data(values=pan_sub))
        .mark_bar(size=24)
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN, title="TME Subtype"),
            xOffset=alt.XOffset("method_display:N"),
            y=alt.Y("orr_pct:Q", title="Objective Response Rate (% CR/PR)"),
            color=alt.Color(
                "method_display:N",
                scale=alt.Scale(
                    domain=["Supervised TCGA k-NN (K=35)", "De Novo Louvain Graph"],
                    range=["#0072B2", "#56B4E9"],
                ),
                title="Classification Method",
            ),
        )
    )
    err_a = (
        alt.Chart(alt.Data(values=pan_sub))
        .mark_errorbar(color="#0F172A", thickness=1.2)
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN, title="TME Subtype"),
            xOffset=alt.XOffset("method_display:N"),
            y=alt.Y("ci_lo_pct:Q", title="Objective Response Rate (% CR/PR)"),
            y2=alt.Y2("ci_hi_pct:Q"),
        )
    )
    panel_a = (bars_a + err_a).properties(
        width=280, height=280, title="A. Pan-iAtlas Baseline Response Rate by Subtype"
    )

    cohort_sub = (
        resp_df.filter(
            (pl.col("cohort") != "PanCancer-iAtlas")
            & (pl.col("method") == "Supervised_TCGA_kNN_Matched")
            & (pl.col("n_samples") >= 2)
        )
        .with_columns(
            (pl.col("response_rate") * 100.0).alias("orr_pct"),
            pl.format(
                "{}% ({}/{})",
                (pl.col("response_rate") * 100.0).round(0).cast(pl.Int64),
                pl.col("n_responders"),
                pl.col("n_samples"),
            ).alias("cell_label"),
        )
        .to_dicts()
    )

    hm_b = (
        alt.Chart(alt.Data(values=cohort_sub))
        .mark_rect(stroke="#FFFFFF", strokeWidth=1.0)
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN, title="Bagaev TME Subtype (TCGA k-NN)"),
            y=alt.Y("cohort:N", title="iAtlas Clinical Trial Cohort"),
            color=alt.Color(
                "orr_pct:Q",
                scale=alt.Scale(domain=[0.0, 65.0], range=["#F8FAFC", "#009E73"], clamp=True),
                title="ORR (%)",
            ),
        )
    )
    txt_b = (
        alt.Chart(alt.Data(values=cohort_sub))
        .mark_text(fontSize=9.5, fontWeight="bold", color="#0F172A")
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN),
            y=alt.Y("cohort:N"),
            text=alt.Text("cell_label:N"),
        )
    )
    panel_b = (hm_b + txt_b).properties(
        width=290, height=280, title="B. Per-Trial ORR % (Responders / Total)"
    )

    combined = _configure_nature_chart(alt.hconcat(panel_a, panel_b, spacing=32))
    _save_svg_and_png(combined, results_dir / "fig2_iatlas_subtype_response_rates.svg")


def plot_fig3_iatlas_deconvolution_integration(data_dir: Path, results_dir: Path) -> None:
    """Fig 3: Consensus-Meta42 Enrichment Across Subtypes + Within-Subtype OOF Prediction Discrimination."""
    enr_df = pl.read_parquet(data_dir / "iatlas_meta42_enrichment_by_bagaev_subtype.parquet")
    joined_df = pl.read_parquet(data_dir / "iatlas_bagaev_joined_meta42_oof.parquet")
    within_df = pl.read_parquet(data_dir / "iatlas_within_subtype_deconv_auc.parquet")
    cv_df = pl.read_parquet(data_dir / "iatlas_bagaev_vs_deconv_cv_benchmark.parquet")

    # Panel A: Top 16 most subtype-differentiated Consensus-Meta42 states (row z-scored across subtypes)
    top_enr = enr_df.head(16)
    hm_rows: list[dict[str, str | float]] = []
    for r in top_enr.iter_rows(named=True):
        vals = [float(r["mean_IE"]), float(r["mean_IEF"]), float(r["mean_F"]), float(r["mean_D"])]
        mu = sum(vals) / 4.0
        sd = max(1e-8, (sum((v - mu) ** 2 for v in vals) / 4.0) ** 0.5)
        markers = str(r.get("top_deg_markers") or "")
        short_m = ", ".join(markers.split(", ")[:2]) if markers else ""
        st_label = f"{r['cell_state']} ({short_m})" if short_m else str(r["cell_state"])
        for st_name, raw_v in zip(SUBTYPE_DOMAIN, vals, strict=True):
            hm_rows.append(
                {
                    "state_label": st_label,
                    "subtype": st_name,
                    "rel_z": (raw_v - mu) / sd,
                }
            )

    panel_a = (
        alt.Chart(alt.Data(values=hm_rows))
        .mark_rect(stroke="#FFFFFF", strokeWidth=0.8)
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN, title="Bagaev TME Subtype"),
            y=alt.Y("state_label:N", title="Top Subtype-Differentiated Meta42 States"),
            color=alt.Color(
                "rel_z:Q",
                scale=alt.Scale(domain=[-1.5, 0.0, 1.5], range=["#0072B2", "#F8FAFC", "#D55E00"], clamp=True),
                title="Relative Enrichment (Z)",
            ),
        )
        .properties(width=200, height=320, title="A. Deconvolved Meta42 State Enrichment")
    )

    # Panel B: OOF Deconv Net Score (Rz - Nz) by Bagaev Subtype & True Clinical Response
    box_recs = _clean_records(
        joined_df.filter(pl.col("response_binary").cast(pl.Float64).is_in([0.0, 1.0]) & pl.col("oof_deconv_net_z").is_not_null())
        .with_columns(
            pl.when(pl.col("response_binary") == 1.0)
            .then(pl.lit("Responder (CR/PR)"))
            .otherwise(pl.lit("Non-Responder (SD/PD)"))
            .alias("clinical_outcome")
        )
        .select(["bagaev_subtype", "clinical_outcome", "oof_deconv_net_z"])
    )

    panel_b = (
        alt.Chart(alt.Data(values=box_recs))
        .mark_boxplot(extent=1.5, size=20)
        .encode(
            x=alt.X("bagaev_subtype:N", sort=SUBTYPE_DOMAIN, title="Bagaev TME Subtype"),
            xOffset=alt.XOffset("clinical_outcome:N"),
            y=alt.Y("oof_deconv_net_z:Q", title="Deconv OOF Net Score (Rz - Nz)"),
            color=alt.Color(
                "clinical_outcome:N",
                scale=alt.Scale(
                    domain=["Responder (CR/PR)", "Non-Responder (SD/PD)"],
                    range=["#009E73", "#D55E00"],
                ),
                title="Clinical Outcome",
            ),
        )
        .properties(width=270, height=150, title="B. Deconv OOF Score Stratifies R vs. NR Inside Subtypes")
    )

    # Panel C: Within-Subtype Deconvolution ROC-AUC
    within_recs = _clean_records(within_df.filter(pl.col("roc_auc").is_finite()))
    rule_50 = alt.Chart(alt.Data(values=[{"x": 0.5}])).mark_rule(stroke="#64748B", strokeDash=[4, 4]).encode(x="x:Q")
    bars_c = (
        alt.Chart(alt.Data(values=within_recs))
        .mark_bar(height=18, color="#0072B2")
        .encode(
            y=alt.Y("stratum:N", sort=["All_Combined", *SUBTYPE_DOMAIN], title="Evaluated Patient Stratum"),
            x=alt.X("roc_auc:Q", scale=alt.Scale(domain=[0.0, 0.80]), title="Deconvolution OOF ROC-AUC"),
        )
    )
    txt_c = (
        alt.Chart(alt.Data(values=within_recs))
        .mark_text(align="left", dx=4, fontSize=9.5, fontWeight="bold", color="#0F172A")
        .encode(
            y=alt.Y("stratum:N", sort=["All_Combined", *SUBTYPE_DOMAIN]),
            x=alt.X("roc_auc:Q"),
            text=alt.Text("roc_auc:Q", format=".3f"),
        )
    )
    panel_c = (bars_c + rule_50 + txt_c).properties(
        width=270, height=135, title="C. Within-Subtype Deconv Discrimination (ROC-AUC)"
    )

    _ = cv_df  # Verified via TSV/Parquet summary
    combined = _configure_nature_chart(
        alt.hconcat(panel_a, alt.vconcat(panel_b, panel_c, spacing=26).resolve_scale(color="independent"), spacing=30)
        .resolve_scale(color="independent")
    )
    _save_svg_and_png(combined, results_dir / "fig3_iatlas_deconvolution_integration.svg")


def plot_fig4_sc23_pseudobulk_conservation(data_dir: Path, results_dir: Path) -> None:
    """Fig 4: 29x29 Correlation Conservation vs TCGA & Subtype Proportions Across All 23 SC Cohorts."""
    sum_df = pl.read_parquet(data_dir / "sc23_pseudobulk_cohort_summary.parquet").sort(
        ["sorting_regime", "offdiag_pearson_r_vs_tcga"], descending=[True, True]
    )
    cohort_order = sum_df["display_name"].to_list()

    panel_a = (
        alt.Chart(alt.Data(values=_clean_records(sum_df.select(["display_name", "offdiag_pearson_r_vs_tcga", "sorting_regime"]))))
        .mark_bar(height=13)
        .encode(
            y=alt.Y("display_name:N", sort=cohort_order, title="Single-Cell Pseudobulk Cohort (9 SC9 + 14 Atlases)"),
            x=alt.X("offdiag_pearson_r_vs_tcga:Q", title="29x29 Fges Correlation Conservation vs. TCGA (Pearson r)"),
            color=alt.Color(
                "sorting_regime:N",
                scale=alt.Scale(
                    domain=["Whole_Tumor", "FACS_Immune_Sorted"],
                    range=["#0072B2", "#E69F00"],
                ),
                title="Cell Selection Strategy",
            ),
        )
        .properties(width=260, height=390, title="A. Conservation of Bagaev 29-Fges Structure vs. TCGA")
    )

    # Stacked subtype proportions across all 23 SC cohorts
    prop_rows: list[dict[str, str | float]] = []
    for r in sum_df.iter_rows(named=True):
        n_tot = max(1, int(r["n_samples"]))
        for st_key, col_k in (("IE", "n_IE"), ("IE/F", "n_IEF"), ("F", "n_F"), ("D", "n_D")):
            prop_rows.append(
                {
                    "display_name": str(r["display_name"]),
                    "subtype": st_key,
                    "pct": 100.0 * float(r[col_k]) / float(n_tot),
                }
            )

    panel_b = (
        alt.Chart(alt.Data(values=prop_rows))
        .mark_bar(height=13)
        .encode(
            y=alt.Y("display_name:N", sort=cohort_order, axis=alt.Axis(labels=False, title=None)),
            x=alt.X("pct:Q", stack="normalize", axis=alt.Axis(format="%"), title="Sample Subtype Proportion"),
            color=alt.Color(
                "subtype:N",
                scale=alt.Scale(domain=SUBTYPE_DOMAIN, range=SUBTYPE_RANGE),
                title="TME Subtype",
            ),
        )
        .properties(width=210, height=390, title="B. Pseudobulk Subtype Distribution")
    )

    combined = _configure_nature_chart(alt.hconcat(panel_a, panel_b, spacing=18).resolve_scale(color="independent"))
    _save_svg_and_png(combined, results_dir / "fig4_sc23_pseudobulk_clustering_conservation.svg")


def plot_fig5_sc9_milopy_response(data_dir: Path, results_dir: Path) -> None:
    """Fig 5: Milopy SC9 Pseudobulk Subtype Responder Rates & Per-Cohort Discrimination."""
    resp_path = data_dir / "sc9_milopy_subtype_response_rates.parquet"
    auc_path = data_dir / "sc9_milopy_cohort_bagaev_aucs.parquet"
    if not resp_path.is_file() or not auc_path.is_file():
        return

    resp_df = pl.read_parquet(resp_path).filter(
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

    bars_a = (
        alt.Chart(alt.Data(values=recs_a))
        .mark_bar(size=22)
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN, title="Bagaev TME Subtype (Pseudobulk TCGA k-NN)"),
            xOffset=alt.XOffset("stratum_label:N"),
            y=alt.Y("orr_pct:Q", title="ICB Responder Rate (% R)"),
            color=alt.Color(
                "stratum_label:N",
                scale=alt.Scale(
                    domain=["Whole-Tumor SC9 (7 Cohorts)", "FACS Immune-Sorted SC9 (2 Cohorts)"],
                    range=["#0072B2", "#E69F00"],
                ),
                title="SC9 Compartment Coverage",
            ),
        )
    )
    err_a = (
        alt.Chart(alt.Data(values=recs_a))
        .mark_errorbar(color="#0F172A", thickness=1.1)
        .encode(
            x=alt.X("subtype:N", sort=SUBTYPE_DOMAIN, title="Bagaev TME Subtype (Pseudobulk TCGA k-NN)"),
            xOffset=alt.XOffset("stratum_label:N"),
            y=alt.Y("ci_lo_pct:Q", title="ICB Responder Rate (% R)"),
            y2=alt.Y2("ci_hi_pct:Q"),
        )
    )
    panel_a = (bars_a + err_a).properties(
        width=280, height=260, title="A. Milopy SC9 Responder Rate by Pseudobulk Subtype"
    )

    auc_df = pl.read_parquet(auc_path).sort("auc_bagaev_immune_minus_fibrotic", descending=True)
    c_order = auc_df["display_name"].to_list()
    rule_50 = alt.Chart(alt.Data(values=[{"x": 0.5}])).mark_rule(stroke="#64748B", strokeDash=[4, 4]).encode(x="x:Q")
    pts_b = (
        alt.Chart(alt.Data(values=_clean_records(auc_df)))
        .mark_bar(height=18)
        .encode(
            y=alt.Y("display_name:N", sort=c_order, title="Milopy Single-Cell ICB Cohort"),
            x=alt.X(
                "auc_bagaev_immune_minus_fibrotic:Q",
                scale=alt.Scale(domain=[0.0, 1.0]),
                title="Pseudobulk Bagaev Immune - Fibrotic ROC-AUC",
            ),
            color=alt.Color(
                "sorting_regime:N",
                scale=alt.Scale(
                    domain=["Whole_Tumor", "FACS_Immune_Sorted"],
                    range=["#0072B2", "#E69F00"],
                ),
                title="Sorting Regime",
            ),
        )
    )
    panel_b = (pts_b + rule_50).properties(
        width=280, height=260, title="B. Per-Cohort Response Discrimination (SC9)"
    )

    combined = _configure_nature_chart(alt.hconcat(panel_a, panel_b, spacing=30).resolve_scale(color="independent"))
    _save_svg_and_png(combined, results_dir / "fig5_sc9_milopy_response_correlation.svg")


def run_step04(data_dir: Path, results_dir: Path) -> None:
    """Render all 5 Bagaev MFP analysis vector figures."""
    results_dir.mkdir(parents=True, exist_ok=True)
    plot_fig1_iatlas_portrait_and_clustering(data_dir, results_dir)
    plot_fig2_iatlas_response_rates(data_dir, results_dir)
    plot_fig3_iatlas_deconvolution_integration(data_dir, results_dir)
    plot_fig4_sc23_pseudobulk_conservation(data_dir, results_dir)
    plot_fig5_sc9_milopy_response(data_dir, results_dir)
    print(f"Step 04 complete. Exported figures to {results_dir}")


def main() -> None:
    paths = get_data_paths()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=paths.output_dir / "bagaev_mfp_analysis")
    parser.add_argument("--results-dir", type=Path, default=paths.results_dir / "bagaev_mfp_analysis")
    args = parser.parse_args()
    run_step04(data_dir=args.data_dir, results_dir=args.results_dir)


if __name__ == "__main__":
    main()
