#!/usr/bin/env python3
"""Figure D: The Consensus-Meta42 Reference Atlas (Both Polar Wheel & Cartesian Signature Matrix).

Generates a publication-grade Nature Methods 4-panel figure:
- Panel (a): The Radial "Meta42 Wheel" (Circular Polar Plot of all 42 states grouped by 6 lineages,
             showing radial bulk responder beta_meta spokes, Milo logFC node size, and top marker genes).
- Panel (b): Cartesian Block-Diagonal Reference Signature Heatmap (Phi_meta across 42 states x top markers).
- Panel (c): Cross-State Transcriptomic Co-Regulation Network (Effector-TLS Hub vs. CAF-Myeloid Barrier Hub).
- Panel (d): Pan-Cancer Cohort Ubiquity Matrix (42 meta-states x 9 scRNA-seq cohorts).
"""

from pathlib import Path
import re
import sys

import altair as alt
import numpy as np
import pandas as pd
import polars as pl
from returns.result import Failure, Result, Success

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(REPO_ROOT / "packages" / "plotting_utils" / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "packages" / "plotting_utils" / "src"))

from plotting_utils import rasterize_svg_to_png

LINEAGE_MAP: dict[str, str] = {
    "M01": "CD8+ / NK T",
    "M04": "CD8+ / NK T",
    "M05": "CD8+ / NK T",
    "M06": "CD8+ / NK T",
    "M13": "CD8+ / NK T",
    "M15": "CD8+ / NK T",
    "M25": "CD8+ / NK T",
    "M33": "CD8+ / NK T",
    "M34": "CD8+ / NK T",
    "M36": "CD8+ / NK T",
    "M03": "CD4+ / Treg",
    "M08": "CD4+ / Treg",
    "M20": "CD4+ / Treg",
    "M21": "CD4+ / Treg",
    "M24": "CD4+ / Treg",
    "M40": "CD4+ / Treg",
    "M29": "B / Plasma / TLS",
    "M30": "B / Plasma / TLS",
    "M37": "B / Plasma / TLS",
    "M09": "Myeloid / TAM",
    "M10": "Myeloid / TAM",
    "M14": "Myeloid / TAM",
    "M19": "Myeloid / TAM",
    "M28": "Myeloid / TAM",
    "M31": "Myeloid / TAM",
    "M32": "Myeloid / TAM",
    "M35": "Myeloid / TAM",
    "M38": "Myeloid / TAM",
    "M39": "Myeloid / TAM",
    "M41": "Myeloid / TAM",
    "M42": "Myeloid / TAM",
    "M02": "Stromal / CAF",
    "M16": "Stromal / CAF",
    "M17": "Stromal / CAF",
    "M23": "Stromal / CAF",
    "M07": "Tumor / Epith.",
    "M11": "Tumor / Epith.",
    "M12": "Tumor / Epith.",
    "M18": "Tumor / Epith.",
    "M22": "Tumor / Epith.",
    "M26": "Tumor / Epith.",
    "M27": "Tumor / Epith.",
}

LINEAGE_ORDER: list[str] = [
    "CD8+ / NK T",
    "CD4+ / Treg",
    "B / Plasma / TLS",
    "Myeloid / TAM",
    "Stromal / CAF",
    "Tumor / Epith.",
]

LINEAGE_COLORS: dict[str, str] = {
    "CD8+ / NK T": "#0284C7",
    "CD4+ / Treg": "#7C3AED",
    "B / Plasma / TLS": "#059669",
    "Myeloid / TAM": "#D97706",
    "Stromal / CAF": "#E11D48",
    "Tumor / Epith.": "#475569",
}


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


def build_figure_d(results_dir: Path) -> Result[Path, str]:
    """Generate Figure D: Consensus-Meta42 Reference Atlas (Polar Wheel + Cartesian Signatures)."""
    phi_path = Path("output/manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark/unified_reference_phi_meta.parquet")
    conc_path = Path("output/manuscript_pipeline/3b_pan_cohort_unified_deconv_benchmark/unified_all_states_detailed_concordance.parquet")
    deg_path = Path("output/manuscript_pipeline/3c_iatlas_instaprism_milopy_benchmark/consensus_metastate_all42_top_degs.parquet")

    if not phi_path.exists() or not conc_path.exists() or not deg_path.exists():
        return Failure("Missing required Meta42 input parquet files.")

    try:
        df_phi = pl.read_parquet(phi_path)
        df_conc = pl.read_parquet(conc_path)
        df_deg = pl.read_parquet(deg_path)

        # Extract per-state metadata from top_degs (rank == 1)
        state_meta = (
            df_deg.filter(pl.col("gene_rank") == 1)
            .select(["cell_state", "state_label", "gene", "state_bulk_beta", "state_milo_logfc"])
            .to_pandas()
        )
        state_meta["lineage"] = state_meta["cell_state"].map(lambda s: LINEAGE_MAP.get(s, "Other"))
        state_meta["lin_order"] = state_meta["lineage"].map(lambda l: LINEAGE_ORDER.index(l) if l in LINEAGE_ORDER else 99)
        state_meta = state_meta.sort_values(["lin_order", "state_bulk_beta"], ascending=[True, False]).reset_index(drop=True)

        ordered_states = state_meta["cell_state"].tolist()

        # ----------------------------------------------------------------------
        # Panel (a): Radial "Meta42 Wheel" (Circular Polar Spoke & Marker Map)
        # ----------------------------------------------------------------------
        n_states = len(state_meta)
        angles = np.linspace(0, 2 * np.pi, n_states, endpoint=False)
        # Start at top-right (75 deg) and go clockwise so top/bottom poles have only 1 label
        thetas = (78.0 * np.pi / 180.0) - angles

        r_base = 0.48
        # Map beta in [-3.2, +2.8] to radius in [0.27, 0.72] so inward spokes never hit center text
        betas = state_meta["state_bulk_beta"].to_numpy()
        r_tip = np.clip(r_base + 0.08 * betas, 0.27, 0.72)

        wheel_rows = []
        for idx, row in state_meta.iterrows():
            th = float(thetas[idx])
            rt = float(r_tip[idx])
            r_lbl = 0.86
            x0, y0 = r_base * np.cos(th), r_base * np.sin(th)
            x1, y1 = rt * np.cos(th), rt * np.sin(th)
            xl, yl = r_lbl * np.cos(th), r_lbl * np.sin(th)
            wheel_rows.append({
                "cell_state": row["cell_state"],
                "gene": row["gene"],
                "short_tag": f"{row['cell_state']}:{row['gene']}",
                "lineage": row["lineage"],
                "beta": float(row["state_bulk_beta"]),
                "abs_milo": float(abs(row["state_milo_logfc"])),
                "x0": float(x0),
                "y0": float(y0),
                "x1": float(x1),
                "y1": float(y1),
                "xl": float(xl),
                "yl": float(yl),
                "side": "right" if np.cos(th) >= 0 else "left",
            })
        wheel_df = pd.DataFrame(wheel_rows)

        # Concentric reference rings (inner = -2, baseline = 0, outer = +2)
        ring_angles = np.linspace(0, 2 * np.pi, 160)
        rings_list = []
        for r_val, r_name in [(0.32, "β = -2"), (r_base, "β = 0 (Baseline)"), (0.64, "β = +2")]:
            for step_i, a in enumerate(ring_angles):
                rings_list.append({
                    "ring": r_name,
                    "step": step_i,
                    "x": float(r_val * np.cos(a)),
                    "y": float(r_val * np.sin(a)),
                })
        rings_df = pd.DataFrame(rings_list)

        ring_chart = (
            alt.Chart(rings_df)
            .mark_line(strokeDash=[3, 3], strokeWidth=1.0, color="#94A3B8")
            .encode(
                x=alt.X("x:Q", scale=alt.Scale(domain=[-1.35, 1.35]), axis=None),
                y=alt.Y("y:Q", scale=alt.Scale(domain=[-1.12, 1.12]), axis=None),
                detail="ring:N",
                order="step:Q",
            )
        )

        # Highlight baseline ring (beta = 0)
        base_ring_df = rings_df[rings_df["ring"] == "β = 0 (Baseline)"]
        base_ring_chart = (
            alt.Chart(base_ring_df)
            .mark_line(strokeWidth=1.5, color="#475569")
            .encode(x="x:Q", y="y:Q", order="step:Q")
        )

        # Radial spokes from r_base to r_tip
        spokes_chart = (
            alt.Chart(wheel_df)
            .mark_rule(strokeWidth=2.4)
            .encode(
                x="x0:Q",
                y="y0:Q",
                x2="x1:Q",
                y2="y1:Q",
                color=alt.Color(
                    "lineage:N",
                    title="Canonical Lineage",
                    scale=alt.Scale(domain=LINEAGE_ORDER, range=[LINEAGE_COLORS[l] for l in LINEAGE_ORDER]),
                    legend=alt.Legend(
                        orient="bottom",
                        columns=2,
                        columnPadding=28,
                        titleFontSize=10,
                        labelFontSize=9.0,
                        symbolSize=55,
                    ),
                ),
            )
        )

        # Tip circles sized by |Milo logFC|
        tips_chart = (
            alt.Chart(wheel_df)
            .mark_circle(opacity=0.92, stroke="#FFFFFF", strokeWidth=0.6)
            .encode(
                x="x1:Q",
                y="y1:Q",
                color=alt.Color("lineage:N", scale=alt.Scale(domain=LINEAGE_ORDER, range=[LINEAGE_COLORS[l] for l in LINEAGE_ORDER])),
                size=alt.Size("abs_milo:Q", scale=alt.Scale(range=[30, 150]), legend=None),
                tooltip=["cell_state", "gene", "lineage", "beta", "abs_milo"],
            )
        )

        # Select 14 well-spaced states around the wheel (every 3rd state = 25.7 deg spacing)
        # Split into right-hemisphere (align='left') and left-hemisphere (align='right') for zero collision!
        lbl_subset = wheel_df.iloc[::3].copy()
        lbl_right = lbl_subset[lbl_subset["side"] == "right"]
        lbl_left = lbl_subset[lbl_subset["side"] == "left"]

        lbl_chart_r = (
            alt.Chart(lbl_right)
            .mark_text(align="left", baseline="middle", fontSize=8.5, fontWeight="bold", color="#0F172A")
            .encode(x="xl:Q", y="yl:Q", text="short_tag:N")
        )
        lbl_chart_l = (
            alt.Chart(lbl_left)
            .mark_text(align="right", baseline="middle", fontSize=8.5, fontWeight="bold", color="#0F172A")
            .encode(x="xl:Q", y="yl:Q", text="short_tag:N")
        )

        center_lbl = (
            alt.Chart(pd.DataFrame([
                {"x": 0.0, "y": 0.06, "t": "Meta42"},
                {"x": 0.0, "y": -0.06, "t": "(K = 42)"},
            ]))
            .mark_text(align="center", baseline="middle", fontSize=9.5, fontWeight="bold", color="#0F172A")
            .encode(x="x:Q", y="y:Q", text="t:N")
        )

        panel_a = (
            alt.layer(
                ring_chart,
                base_ring_chart,
                spokes_chart,
                tips_chart,
                lbl_chart_r,
                lbl_chart_l,
                center_lbl,
                _panel_letter_layer("(a)"),
            )
            .properties(width=355, height=315)
        )

        # ----------------------------------------------------------------------
        # Panel (b): Cartesian Block-Diagonal Reference Signature Heatmap (Phi_meta)
        # ----------------------------------------------------------------------
        # Select 21 unique top marker genes ordered along the 42 states so the block-diagonal shines
        rep_genes: list[str] = []
        for idx, row in state_meta.iterrows():
            g = row["gene"]
            if g in df_phi.columns and g not in rep_genes:
                rep_genes.append(g)
        rep_genes = rep_genes[:21]

        phi_pd = df_phi.select(["cluster"] + rep_genes).to_pandas().set_index("cluster")
        phi_pd = phi_pd.loc[ordered_states]
        # Log2 specificity relative to column median, scaled so each marker peaks cleanly along the diagonal
        col_mean = phi_pd.mean(axis=0) + 1e-6
        log2_spec = np.log2((phi_pd + 1e-6) / col_mean)
        z_mat = (log2_spec - log2_spec.mean(axis=0)) / (log2_spec.std(axis=0) + 1e-8)
        z_mat = z_mat.clip(-1.0, 2.8)

        heat_rows = []
        for st in ordered_states:
            for g in rep_genes:
                heat_rows.append({
                    "state": st,
                    "gene": g,
                    "z_expr": float(z_mat.loc[st, g]),
                    "lineage": LINEAGE_MAP.get(st, "Other"),
                })
        heat_df = pd.DataFrame(heat_rows)

        shown_y_states = ordered_states[::2]

        panel_b_heat = (
            alt.Chart(heat_df)
            .mark_rect()
            .encode(
                x=alt.X(
                    "gene:N",
                    sort=rep_genes,
                    title="Canonical Lineage & State Marker Genes (HGNC)",
                    axis=alt.Axis(labelAngle=-60, labelFontSize=8.5, titlePadding=12),
                ),
                y=alt.Y(
                    "state:N",
                    sort=ordered_states,
                    title="Consensus Meta-States (M01–M42)",
                    axis=alt.Axis(values=shown_y_states, labelFontSize=8.5, titlePadding=12),
                ),
                color=alt.Color(
                    "z_expr:Q",
                    title="Specificity Z",
                    scale=alt.Scale(scheme="viridis", domain=[-1.0, 2.8]),
                    legend=alt.Legend(orient="bottom", gradientLength=130, titleFontSize=10, labelFontSize=9.5),
                ),
                tooltip=["state", "gene", "z_expr", "lineage"],
            )
        )

        panel_b = (
            alt.layer(panel_b_heat, _panel_letter_layer("(b)"))
            .properties(width=355, height=315)
        )

        # ----------------------------------------------------------------------
        # Panel (c): Cross-State Co-Regulation & Antagonism Landscape (Beta vs. Milo logFC)
        # ----------------------------------------------------------------------
        node_df = state_meta.copy()
        cluster_counts = df_conc.group_by("meta_state").len().to_pandas().rename(columns={"meta_state": "cell_state", "len": "n_clusters"})
        node_df = node_df.merge(cluster_counts, on="cell_state", how="left")
        node_df["n_clusters"] = node_df["n_clusters"].fillna(5)

        corr_mat = np.corrcoef(z_mat.to_numpy())
        edge_rows = []
        for i in range(len(ordered_states)):
            for j in range(i + 1, len(ordered_states)):
                if corr_mat[i, j] > 0.48:
                    r_i = node_df[node_df["cell_state"] == ordered_states[i]].iloc[0]
                    r_j = node_df[node_df["cell_state"] == ordered_states[j]].iloc[0]
                    edge_rows.append({
                        "x1": float(r_i["state_milo_logfc"]),
                        "y1": float(r_i["state_bulk_beta"]),
                        "x2": float(r_j["state_milo_logfc"]),
                        "y2": float(r_j["state_bulk_beta"]),
                        "corr": float(corr_mat[i, j]),
                    })
        edge_df = pd.DataFrame(edge_rows)

        zero_x = alt.Chart(pd.DataFrame({"x": [0.0]})).mark_rule(color="#CBD5E1", strokeDash=[3, 3]).encode(x="x:Q")
        zero_y = alt.Chart(pd.DataFrame({"y": [0.0]})).mark_rule(color="#CBD5E1", strokeDash=[3, 3]).encode(y="y:Q")

        edges_c = (
            alt.Chart(edge_df)
            .mark_rule(color="#94A3B8", opacity=0.35, strokeWidth=1.2)
            .encode(x="x1:Q", y="y1:Q", x2="x2:Q", y2="y2:Q")
        )

        nodes_c = (
            alt.Chart(node_df)
            .mark_circle(opacity=0.9, stroke="#0F172A", strokeWidth=0.6)
            .encode(
                x=alt.X("state_milo_logfc:Q", title="Single-Cell Milo Neighborhood Enrichment (logFC)", scale=alt.Scale(domain=[-9.5, 9.5])),
                y=alt.Y("state_bulk_beta:Q", title="Bulk ICB Responder Effect (β_meta)", scale=alt.Scale(domain=[-3.6, 3.3])),
                color=alt.Color(
                    "state_bulk_beta:Q",
                    title="Responder β",
                    scale=alt.Scale(scheme="redblue", domain=[-3.0, 3.0], reverse=True),
                    legend=alt.Legend(orient="bottom", gradientLength=130, titleFontSize=10, labelFontSize=9.5),
                ),
                size=alt.Size("n_clusters:Q", title="Clusters", scale=alt.Scale(range=[45, 220]), legend=None),
                tooltip=["cell_state", "state_label", "state_milo_logfc", "state_bulk_beta", "n_clusters"],
            )
        )

        # Explicit non-overlapping label coordinates for key hub states
        hub_annos = pd.DataFrame([
            {"x": 6.96, "y": 2.92, "lbl": "M13 (CXCL13)", "align": "center"},
            {"x": 4.40, "y": 2.28, "lbl": "M14 (CXCL9)", "align": "right"},
            {"x": 5.50, "y": 1.85, "lbl": "M15 (TCF7)", "align": "left"},
            {"x": 2.50, "y": 1.65, "lbl": "M37 (MS4A1)", "align": "right"},
            {"x": 5.35, "y": 1.30, "lbl": "M34 (GZMK)", "align": "left"},
            {"x": -2.40, "y": -1.35, "lbl": "M19 (S100A8)", "align": "left"},
            {"x": -3.10, "y": -1.85, "lbl": "M21 (FOXP3)", "align": "left"},
            {"x": -7.60, "y": -2.85, "lbl": "M23 (FAP)", "align": "left"},
        ])

        lbl_c_left = (
            alt.Chart(hub_annos[hub_annos["align"] == "left"])
            .mark_text(align="left", baseline="middle", fontSize=8.5, fontWeight="bold", color="#0F172A")
            .encode(x="x:Q", y="y:Q", text="lbl:N")
        )
        lbl_c_right = (
            alt.Chart(hub_annos[hub_annos["align"] == "right"])
            .mark_text(align="right", baseline="middle", fontSize=8.5, fontWeight="bold", color="#0F172A")
            .encode(x="x:Q", y="y:Q", text="lbl:N")
        )
        lbl_c_center = (
            alt.Chart(hub_annos[hub_annos["align"] == "center"])
            .mark_text(align="center", baseline="bottom", fontSize=8.5, fontWeight="bold", color="#0F172A")
            .encode(x="x:Q", y="y:Q", text="lbl:N")
        )

        panel_c = (
            alt.layer(zero_x, zero_y, edges_c, nodes_c, lbl_c_left, lbl_c_right, lbl_c_center, _panel_letter_layer("(c)"))
            .properties(width=355, height=255)
        )

        # ----------------------------------------------------------------------
        # Panel (d): Pan-Cancer Cohort Ubiquity Matrix (42 States x 9 Cohorts)
        # ----------------------------------------------------------------------
        df_conc_pd = df_conc.to_pandas()
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

        ub_df = (
            df_conc_pd.groupby(["meta_state", "cohort_clean"], as_index=False)
            .size()
            .rename(columns={"size": "n_clusters"})
        )
        ub_df["lineage"] = ub_df["meta_state"].map(lambda s: LINEAGE_MAP.get(s, "Other"))

        panel_d_dots = (
            alt.Chart(ub_df)
            .mark_rect(stroke="#FFFFFF", strokeWidth=0.5)
            .encode(
                x=alt.X(
                    "meta_state:N",
                    sort=ordered_states,
                    title="Consensus Meta-States (M01–M42, Ordered by Lineage)",
                    axis=alt.Axis(values=ordered_states[::3], labelAngle=-90, labelFontSize=8.5, titlePadding=12),
                ),
                y=alt.Y(
                    "cohort_clean:N",
                    title=None,
                    axis=alt.Axis(labelFontSize=9.5, labelLimit=180),
                ),
                color=alt.Color(
                    "n_clusters:Q",
                    title="Fused Clusters",
                    scale=alt.Scale(scheme="blues", domain=[0, 5]),
                    legend=alt.Legend(orient="bottom", gradientLength=130, titleFontSize=10, labelFontSize=9.5),
                ),
                tooltip=["meta_state", "cohort_clean", "n_clusters", "lineage"],
            )
        )

        panel_d = (
            alt.layer(panel_d_dots, _panel_letter_layer("(d)"))
            .properties(width=355, height=255)
        )

        # ----------------------------------------------------------------------
        # Assemble Final Layout
        # ----------------------------------------------------------------------
        row1 = alt.hconcat(panel_a, panel_b, spacing=42).resolve_scale(color="independent", size="independent")
        row2 = alt.hconcat(panel_c, panel_d, spacing=42).resolve_scale(color="independent", size="independent")

        final_chart = (
            alt.vconcat(row1, row2, spacing=34)
            .resolve_scale(color="independent", size="independent")
            .configure_axis(labelFontSize=10.5, titleFontSize=11.5, titleFontWeight="bold")
            .configure_view(stroke="#CBD5E1", strokeWidth=0.75)
            .configure(
                font="system-ui, -apple-system, sans-serif",
                background="#FFFFFF",
                padding={"left": 25, "right": 20, "top": 28, "bottom": 15},
            )
        )

        out_svg = results_dir / "fig_d_meta42_atlas_polar_and_signatures.svg"
        out_png = results_dir / "fig_d_meta42_atlas_polar_and_signatures.png"
        results_dir.mkdir(parents=True, exist_ok=True)
        final_chart.save(str(out_svg))

        post_process_svg(out_svg, width_mm=180.0, min_font_pt=12.0)
        rasterize_svg_to_png(out_svg, out_png, scale=3.125, background="white")

        return Success(out_svg)

    except Exception as exc:
        return Failure(f"Failed to generate Figure D: {exc}")


def main() -> None:
    results_dir = Path("results/pipeline_optimization_and_meta42")
    match build_figure_d(results_dir):
        case Success(p):
            print(f"Generated Figure D: {p}")
        case Failure(err):
            print(f"Error Figure D: {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
