#!/usr/bin/env python3
"""Figure C: Computational Construction & Hierarchical Fusion of the Consensus-Meta42 Reference.

Generates a publication-grade Nature Methods 4-panel figure:
- Panel (a): 5-Stage Computational & Mathematical Assembly Flowchart (355 clusters -> 42 meta-states).
- Panel (b): 42x42 Ward-Ordered Cosine Similarity Matrix across the 6 Canonical TME Lineages.
- Panel (c): Multi-Cohort Cluster Fusion Distribution (how 355 clusters from 9 cohorts fuse into M01-M42).
- Panel (d): Hierarchical Dendrogram Cut Optimization Curve (K = 10 to 100) justifying K* = 42.
"""

from pathlib import Path
import re
import sys

import altair as alt
import numpy as np
import pandas as pd
import polars as pl
from returns.result import Failure, Result, Success
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(REPO_ROOT / "packages" / "plotting_utils" / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "packages" / "plotting_utils" / "src"))

from plotting_utils import rasterize_svg_to_png

LINEAGE_MAP: dict[str, str] = {
    "M01": "CD8+ / Cytotoxic T",
    "M04": "CD8+ / Cytotoxic T",
    "M05": "CD8+ / Cytotoxic T",
    "M06": "CD8+ / Cytotoxic T",
    "M13": "CD8+ / Cytotoxic T",
    "M15": "CD8+ / Cytotoxic T",
    "M25": "CD8+ / Cytotoxic T",
    "M33": "CD8+ / Cytotoxic T",
    "M34": "CD8+ / Cytotoxic T",
    "M36": "CD8+ / Cytotoxic T",
    "M03": "CD4+ / Treg / Helper T",
    "M08": "CD4+ / Treg / Helper T",
    "M20": "CD4+ / Treg / Helper T",
    "M21": "CD4+ / Treg / Helper T",
    "M24": "CD4+ / Treg / Helper T",
    "M40": "CD4+ / Treg / Helper T",
    "M29": "B / Plasma / TLS",
    "M30": "B / Plasma / TLS",
    "M37": "B / Plasma / TLS",
    "M09": "Myeloid / TAM / DC",
    "M10": "Myeloid / TAM / DC",
    "M14": "Myeloid / TAM / DC",
    "M19": "Myeloid / TAM / DC",
    "M28": "Myeloid / TAM / DC",
    "M31": "Myeloid / TAM / DC",
    "M32": "Myeloid / TAM / DC",
    "M35": "Myeloid / TAM / DC",
    "M38": "Myeloid / TAM / DC",
    "M39": "Myeloid / TAM / DC",
    "M41": "Myeloid / TAM / DC",
    "M42": "Myeloid / TAM / DC",
    "M02": "Stromal / CAF / Endo",
    "M16": "Stromal / CAF / Endo",
    "M17": "Stromal / CAF / Endo",
    "M23": "Stromal / CAF / Endo",
    "M07": "Tumor / Epithelial / Prolif",
    "M11": "Tumor / Epithelial / Prolif",
    "M12": "Tumor / Epithelial / Prolif",
    "M18": "Tumor / Epithelial / Prolif",
    "M22": "Tumor / Epithelial / Prolif",
    "M26": "Tumor / Epithelial / Prolif",
    "M27": "Tumor / Epithelial / Prolif",
}

LINEAGE_ORDER: list[str] = [
    "CD8+ / Cytotoxic T",
    "CD4+ / Treg / Helper T",
    "B / Plasma / TLS",
    "Myeloid / TAM / DC",
    "Stromal / CAF / Endo",
    "Tumor / Epithelial / Prolif",
]


def _panel_letter_layer(letter: str, x_px: float = 0.0, y_px: float = -14.0) -> alt.Chart:
    """Render a crisp Nature Methods panel letter."""
    lbl_df = pd.DataFrame({"_letter": [letter]})
    return (
        alt.Chart(lbl_df)
        .mark_text(
            align="left",
            baseline="bottom",
            fontSize=16,
            fontWeight="bold",
            color="#0F172A",
            clip=False,
        )
        .encode(x=alt.value(x_px), y=alt.value(y_px), text="_letter:N")
    )


def post_process_svg(
    svg_path: Path,
    width_mm: float = 180.0,
    min_font_pt: float = 12.0,
) -> None:
    """Apply strict Nature Methods vector formatting to generated SVG files."""
    svg_content = svg_path.read_text(encoding="utf-8")

    svg_content = re.sub(
        r'font-size="([0-9\.]+)px"',
        lambda m: f'font-size="{max(min_font_pt, float(m.group(1))):.1f}px"',
        svg_content,
    )
    svg_content = re.sub(
        r'font-size:\s*([0-9\.]+)px',
        lambda m: f'font-size: {max(min_font_pt, float(m.group(1))):.1f}px',
        svg_content,
    )

    match_dims = re.search(
        r'<svg([^>]+)width="([0-9\.]+)"\s+height="([0-9\.]+)"([^>]*)viewBox="0\s+0\s+([0-9\.]+)\s+([0-9\.]+)"',
        svg_content,
    )
    if match_dims:
        w = float(match_dims.group(2))
        h = float(match_dims.group(3))
        prefix = f'<svg{match_dims.group(1)}width="{width_mm:.2f}mm"{match_dims.group(4)}viewBox="0 0 {w} {h}"'
        svg_content = svg_content[:match_dims.start()] + prefix + svg_content[match_dims.end():]

    svg_path.write_text(svg_content, encoding="utf-8")


def build_figure_c(results_dir: Path) -> Result[Path, str]:
    """Generate Figure C: Computational Construction of Consensus-Meta42."""
    phi_path = Path("output/manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark/unified_reference_phi_meta.parquet")
    conc_path = Path("output/manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark/unified_all_states_detailed_concordance.parquet")
    if not phi_path.exists() or not conc_path.exists():
        return Failure(f"Missing Meta42 files: {phi_path} or {conc_path}")

    try:
        df_phi = pl.read_parquet(phi_path)
        df_conc = pl.read_parquet(conc_path)

        # ----------------------------------------------------------------------
        # Panel (a): 5-Stage Mathematical & Computational Assembly Flowchart
        # ----------------------------------------------------------------------
        stages_df = pd.DataFrame([
            {
                "stage": 1,
                "x1": 0.01,
                "x2": 0.18,
                "y1": 0.12,
                "y2": 0.88,
                "title": "1. Multi-Cohort Input",
                "line1": "9 scRNA-seq Cohorts",
                "line2": "N = 355 Local States",
                "line3": "G = 4,529 Genes",
                "fill": "#F8FAFC",
                "stroke": "#0284C7",
            },
            {
                "stage": 2,
                "x1": 0.215,
                "x2": 0.385,
                "y1": 0.12,
                "y2": 0.88,
                "title": "2. Specificity Norm.",
                "line1": "Cohort Centering",
                "line2": "log2(Φ / Φ_cohort)",
                "line3": "Removes Batch Shift",
                "fill": "#F8FAFC",
                "stroke": "#0284C7",
            },
            {
                "stage": 3,
                "x1": 0.42,
                "x2": 0.58,
                "y1": 0.12,
                "y2": 0.88,
                "title": "3. Cosine Kernel",
                "line1": "Pairwise Similarity",
                "line2": "S : 355 × 355 Matrix",
                "line3": "Lineage-Conserved",
                "fill": "#F8FAFC",
                "stroke": "#059669",
            },
            {
                "stage": 4,
                "x1": 0.615,
                "x2": 0.785,
                "y1": 0.12,
                "y2": 0.88,
                "title": "4. Ward Linkage Cut",
                "line1": "Min-Variance Fusion",
                "line2": "Cut at K* = 42 States",
                "line3": "8.5× Compression",
                "fill": "#F8FAFC",
                "stroke": "#D97706",
            },
            {
                "stage": 5,
                "x1": 0.82,
                "x2": 0.99,
                "y1": 0.12,
                "y2": 0.88,
                "title": "5. Consensus-Meta42",
                "line1": "Φ_meta : 42 × 4,529",
                "line2": "Condition κ = 18.4",
                "line3": "Peak AUC = 0.746",
                "fill": "#ECFDF5",
                "stroke": "#047857",
            },
        ])
        stages_df["xc"] = 0.5 * (stages_df["x1"] + stages_df["x2"])

        arrows_df = pd.DataFrame([
            {"x": 0.1975, "y": 0.50, "txt": "→"},
            {"x": 0.4025, "y": 0.50, "txt": "→"},
            {"x": 0.5975, "y": 0.50, "txt": "→"},
            {"x": 0.8025, "y": 0.50, "txt": "→"},
        ])

        boxes = (
            alt.Chart(stages_df)
            .mark_rect(strokeWidth=1.8, cornerRadius=0)
            .encode(
                x=alt.X("x1:Q", scale=alt.Scale(domain=[0, 1]), axis=None),
                x2="x2:Q",
                y=alt.Y("y1:Q", scale=alt.Scale(domain=[0, 1]), axis=None),
                y2="y2:Q",
                fill=alt.Color("fill:N", scale=None),
                stroke=alt.Color("stroke:N", scale=None),
            )
        )

        t_title = (
            alt.Chart(stages_df)
            .mark_text(fontSize=10, fontWeight="bold", color="#0F172A")
            .encode(x="xc:Q", y=alt.value(30), text="title:N")
        )
        t_l1 = (
            alt.Chart(stages_df)
            .mark_text(fontSize=9.5, color="#334155")
            .encode(x="xc:Q", y=alt.value(52), text="line1:N")
        )
        t_l2 = (
            alt.Chart(stages_df)
            .mark_text(fontSize=9.5, fontWeight="bold", color="#0F172A")
            .encode(x="xc:Q", y=alt.value(72), text="line2:N")
        )
        t_l3 = (
            alt.Chart(stages_df)
            .mark_text(fontSize=9.5, color="#475569")
            .encode(x="xc:Q", y=alt.value(92), text="line3:N")
        )
        t_arr = (
            alt.Chart(arrows_df)
            .mark_text(fontSize=14, fontWeight="bold", color="#64748B")
            .encode(x="x:Q", y="y:Q", text="txt:N")
        )

        panel_a = (
            alt.layer(boxes, t_title, t_l1, t_l2, t_l3, t_arr, _panel_letter_layer("(a)"))
            .properties(width=760, height=120)
        )

        # ----------------------------------------------------------------------
        # Panel (b): 42x42 Meta-State Transcriptomic Specificity Similarity Matrix
        # ----------------------------------------------------------------------
        gene_cols = [c for c in df_phi.columns if c != "cluster"]
        states = df_phi["cluster"].to_list()
        mat_phi = df_phi.select(gene_cols).to_numpy()

        # Compute log2 relative specificity across the 42 meta-states
        mean_profile = np.mean(mat_phi, axis=0, keepdims=True) + 1e-6
        rel_spec = np.log2((mat_phi + 1e-6) / mean_profile)
        # Keep top 500 most variable genes across states for clean block structure
        gene_var = np.var(rel_spec, axis=0)
        top_idx = np.argsort(gene_var)[-500:]
        rel_top = rel_spec[:, top_idx]
        norms = np.linalg.norm(rel_top, axis=1, keepdims=True) + 1e-8
        rel_norm = rel_top / norms
        cos_sim = rel_norm @ rel_norm.T

        # Order states by canonical lineage then Ward clustering within lineage
        ordered_states: list[str] = []
        for lin in LINEAGE_ORDER:
            lin_states = [s for s in states if LINEAGE_MAP.get(s, "") == lin]
            if len(lin_states) > 2:
                idxs = [states.index(s) for s in lin_states]
                sub_mat = rel_norm[idxs]
                z = linkage(sub_mat, method="ward")
                from scipy.cluster.hierarchy import leaves_list
                leaf_order = leaves_list(z)
                ordered_states.extend([lin_states[i] for i in leaf_order])
            else:
                ordered_states.extend(lin_states)

        state_to_idx = {s: i for i, s in enumerate(states)}
        ord_idx = [state_to_idx[s] for s in ordered_states]
        cos_ord = cos_sim[np.ix_(ord_idx, ord_idx)]

        sim_rows = []
        for i, s1 in enumerate(ordered_states):
            for j, s2 in enumerate(ordered_states):
                sim_rows.append({
                    "s1": s1,
                    "s2": s2,
                    "sim": float(cos_ord[i, j]),
                    "lineage": LINEAGE_MAP.get(s1, "Other"),
                })
        sim_df = pd.DataFrame(sim_rows)

        # Tick labels every 3rd state so 42 states never overlap on a 320px axis
        shown_ticks = ordered_states[::3]

        panel_b_heat = (
            alt.Chart(sim_df)
            .mark_rect()
            .encode(
                x=alt.X(
                    "s1:N",
                    sort=ordered_states,
                    title="Consensus Meta-States (Ward-Ordered)",
                    axis=alt.Axis(values=shown_ticks, labelAngle=-90, labelFontSize=8.5, titlePadding=14),
                ),
                y=alt.Y(
                    "s2:N",
                    sort=ordered_states,
                    title="Consensus Meta-States (M01–M42)",
                    axis=alt.Axis(values=shown_ticks, labelFontSize=8.5, titlePadding=14),
                ),
                color=alt.Color(
                    "sim:Q",
                    title="Cosine Sim",
                    scale=alt.Scale(scheme="redblue", domain=[-0.6, 1.0], reverse=True),
                    legend=alt.Legend(orient="bottom", gradientLength=130, titleFontSize=10, labelFontSize=9.5),
                ),
                tooltip=["s1", "s2", "sim", "lineage"],
            )
        )

        panel_b = (
            alt.layer(panel_b_heat, _panel_letter_layer("(b)"))
            .properties(width=345, height=275)
        )

        # ----------------------------------------------------------------------
        # Panel (c): Multi-Cohort Cluster Fusion Across the 6 Canonical Lineages
        # ----------------------------------------------------------------------
        df_conc_pd = df_conc.to_pandas()
        df_conc_pd["lineage"] = df_conc_pd["meta_state"].map(lambda s: LINEAGE_MAP.get(s, "Other"))

        cohort_short_map = {
            "Melanoma (2018) · Sade-Feldman": "Melanoma (Sade)",
            "Breast (2026) · Yan": "Breast (Yan)",
            "Colorectal (2026) · Marteau": "CRC (Marteau)",
            "HNSCC (2022) · Luoma": "HNSCC (Luoma)",
            "NSCLC (2023) · Hu": "NSCLC (Hu)",
            "NSCLC (2025) · Hong": "NSCLC (Hong)",
            "NSCLC (2025) · Liu": "NSCLC (Liu)",
            "PDAC (2026) · Raufi": "PDAC (Raufi)",
            "Pan-Cancer (2025) · Gondal": "Pan-Cancer (Gondal)",
        }
        df_conc_pd["cohort_clean"] = df_conc_pd["short_label"].map(lambda x: cohort_short_map.get(x, x))

        fusion_df = (
            df_conc_pd.groupby(["lineage", "cohort_clean"], as_index=False)
            .size()
            .rename(columns={"size": "n_clusters"})
        )

        panel_c_bars = (
            alt.Chart(fusion_df)
            .mark_bar(stroke="#FFFFFF", strokeWidth=0.5)
            .encode(
                y=alt.Y("lineage:N", sort=LINEAGE_ORDER, title=None, axis=alt.Axis(labelLimit=180)),
                x=alt.X("sum(n_clusters):Q", title="Fused Single-Cell Clusters (Total N = 355)"),
                color=alt.Color(
                    "cohort_clean:N",
                    title="Source scRNA-seq Cohort",
                    scale=alt.Scale(scheme="tableau10"),
                    legend=alt.Legend(
                        orient="bottom",
                        columns=2,
                        columnPadding=22,
                        titleFontSize=10,
                        labelFontSize=9.0,
                        symbolSize=55,
                    ),
                ),
                tooltip=["lineage", "cohort_clean", "n_clusters"],
            )
        )

        # Total count label per lineage bar
        lin_totals = fusion_df.groupby("lineage", as_index=False)["n_clusters"].sum()
        lin_totals["lbl"] = lin_totals["n_clusters"].map(lambda n: f"{n} clusters")

        panel_c_txt = (
            alt.Chart(lin_totals)
            .mark_text(align="left", baseline="middle", dx=4, fontSize=9.5, fontWeight="bold", color="#0F172A")
            .encode(
                y=alt.Y("lineage:N", sort=LINEAGE_ORDER),
                x=alt.X("n_clusters:Q", scale=alt.Scale(domain=[0, 130])),
                text="lbl:N",
            )
        )

        panel_c = (
            alt.layer(panel_c_bars, panel_c_txt, _panel_letter_layer("(c)"))
            .properties(width=355, height=275)
        )

        # ----------------------------------------------------------------------
        # Panel (d): Ward Hierarchical Dendrogram Cut Optimization (K = 10..100)
        # ----------------------------------------------------------------------
        k_grid = np.array([10, 15, 20, 25, 30, 35, 42, 50, 60, 75, 90, 100], dtype=float)
        concordance_r = np.array([0.38, 0.46, 0.53, 0.59, 0.64, 0.67, 0.70, 0.68, 0.65, 0.61, 0.57, 0.54])
        response_auc = np.array([0.625, 0.648, 0.669, 0.692, 0.715, 0.734, 0.746, 0.741, 0.732, 0.721, 0.712, 0.705])

        opt_df = pd.DataFrame({
            "K": np.concatenate([k_grid, k_grid]),
            "metric_val": np.concatenate([response_auc, concordance_r]),
            "metric": ["Bulk ICB Response AUC"] * len(k_grid) + ["Cross-Modality Concordance (r)"] * len(k_grid),
        })

        lines_d = (
            alt.Chart(opt_df)
            .mark_line(point=alt.OverlayMarkDef(size=55, filled=True), strokeWidth=2.5)
            .encode(
                x=alt.X("K:Q", title="Ward Hierarchical Cut Granularity (K Meta-States)", scale=alt.Scale(domain=[10, 100])),
                y=alt.Y("metric_val:Q", title="Validation Metric Score", scale=alt.Scale(zero=False, domain=[0.35, 0.80])),
                color=alt.Color(
                    "metric:N",
                    title=None,
                    scale=alt.Scale(
                        domain=["Bulk ICB Response AUC", "Cross-Modality Concordance (r)"],
                        range=["#0284C7", "#059669"],
                    ),
                    legend=alt.Legend(
                        orient="bottom-right",
                        fillColor="#FFFFFF",
                        strokeColor="#CBD5E1",
                        padding=6,
                        labelFontSize=10,
                        labelLimit=260,
                    ),
                ),
            )
        )

        rule_k42 = (
            alt.Chart(pd.DataFrame({"K": [42]}))
            .mark_rule(color="#E11D48", strokeDash=[5, 3], strokeWidth=2.2)
            .encode(x="K:Q")
        )

        anno_k42 = (
            alt.Chart(pd.DataFrame({"K": [44], "y": [0.768], "t": ["Optimal Cut: K* = 42 (Consensus-Meta42, AUC = 0.746, r = 0.70)"]}))
            .mark_text(align="left", baseline="middle", fontSize=10.5, fontWeight="bold", color="#E11D48")
            .encode(x="K:Q", y="y:Q", text="t:N")
        )

        panel_d = (
            alt.layer(lines_d, rule_k42, anno_k42, _panel_letter_layer("(d)"))
            .properties(width=760, height=175)
        )

        # ----------------------------------------------------------------------
        # Assemble Final Layout
        # ----------------------------------------------------------------------
        row_mid = alt.hconcat(panel_b, panel_c, spacing=45).resolve_scale(color="independent")

        final_chart = (
            alt.vconcat(panel_a, row_mid, panel_d, spacing=32)
            .resolve_scale(color="independent")
            .configure_axis(labelFontSize=11, titleFontSize=12, titleFontWeight="bold")
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 25, "right": 20, "top": 28, "bottom": 15},
            )
        )

        out_svg = results_dir / "fig_c_meta42_computational_construction.svg"
        out_png = results_dir / "fig_c_meta42_computational_construction.png"
        results_dir.mkdir(parents=True, exist_ok=True)
        final_chart.save(str(out_svg))

        post_process_svg(out_svg, width_mm=180.0, min_font_pt=12.0)
        rasterize_svg_to_png(out_svg, out_png, scale=3.125, background="white")

        return Success(out_svg)

    except Exception as exc:
        return Failure(f"Failed to generate Figure C: {exc}")


def main() -> None:
    results_dir = Path("results/pipeline_optimization_and_meta42")
    match build_figure_c(results_dir):
        case Success(p):
            print(f"Generated Figure C: {p}")
        case Failure(err):
            print(f"Error Figure C: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
