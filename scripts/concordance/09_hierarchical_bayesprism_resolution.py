#!/usr/bin/env python3
"""
Step 9: Resolve Sibling State Discrepancy via Hierarchical BayesPrism
and Data-Driven Collinearity Consolidation.
Implements:
1. Reference collinearity pruning and cluster consolidation (r >= 0.85)
2. Condition number comparison (kappa_original vs kappa_consolidated)
3. Hierarchical BayesPrism (cell_type lineages vs cell_state sub-clusters)
4. Consolidated BayesPrism deconvolution
5. Two-tier standardized logistic regression vs patient response
6. Single-cell concordances and publication vector SVG export
Strict functional Python with returns, Pydantic, Polars, and Altair.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import altair as alt  # type: ignore
import bayesprism as bp  # type: ignore
import numpy as np
import pandas as pd  # type: ignore
from pydantic import BaseModel, ConfigDict
from returns.maybe import Nothing, Some
from returns.result import Failure, Result, Success
import polars as pl
from scipy import stats  # type: ignore
import vl_convert as vlc  # type: ignore


class ResolutionConfig(BaseModel):
    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)
    ref_path: Path
    hugo_tpm_path: Path
    hugo_clin_path: Path
    sc_effects_path: Path
    out_dir: Path
    results_dir: Path
    chain_length: int = 40
    burn_in: int = 20
    collinear_threshold: float = 0.85


def parse_args() -> ResolutionConfig:
    parser = argparse.ArgumentParser(
        description="Resolve sibling state discrepancy via hierarchical BayesPrism."
    )
    parser.add_argument(
        "--ref",
        type=Path,
        default=Path("output/output/sade_feldman_deconv_validation/reference_phi_res0.5.parquet"),
        help="Path to reference phi parquet",
    )
    parser.add_argument(
        "--hugo-tpm",
        type=Path,
        default=Path(
            "scratch/lair/CBioPortalDataset-Hugo-iAtlas/mel_iatlas_hugo_ucla_2016/mel_iatlas_hugo_ucla_2016/data_mrna_seq_tpm.txt"
        ),
        help="Path to Hugo bulk TPM file",
    )
    parser.add_argument(
        "--hugo-clin",
        type=Path,
        default=Path(
            "scratch/lair/CBioPortalDataset-Hugo-iAtlas/mel_iatlas_hugo_ucla_2016/mel_iatlas_hugo_ucla_2016/data_clinical_sample.txt"
        ),
        help="Path to Hugo clinical sample file",
    )
    parser.add_argument(
        "--sc-effects",
        type=Path,
        default=Path("output/concordance/sc_patient_response_effects.parquet"),
        help="Path to single-cell response effects parquet",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("output/concordance"),
        help="Output directory for parquets",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results/concordance"),
        help="Output directory for figures",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.85,
        help="Collinearity correlation threshold for cluster consolidation (e.g. 0.85)",
    )
    args = parser.parse_args()
    return ResolutionConfig(
        ref_path=args.ref,
        hugo_tpm_path=args.hugo_tpm,
        hugo_clin_path=args.hugo_clin,
        sc_effects_path=args.sc_effects,
        out_dir=args.out_dir,
        results_dir=args.results_dir,
        collinear_threshold=args.threshold,
    )


def map_coarse_lineage(cluster_name: str) -> str:
    """Pure mapping of fine cell state to major parental lineage."""
    name_lower = cluster_name.lower()
    if "cytotoxic" in name_lower or "temra" in name_lower or "tem/trm" in name_lower:
        return "Cytotoxic_T"
    if "regulatory" in name_lower:
        return "Treg"
    if "b cell" in name_lower or "plasma" in name_lower:
        return "B_lineage"
    if "macrophage" in name_lower:
        return "Macrophage"
    if "pdc" in name_lower:
        return "pDC"
    return "Other"


def consolidate_collinear_reference(
    ref_mat: np.ndarray,
    clusters: list[str],
    threshold: float,
) -> tuple[np.ndarray, list[str], dict[str, list[str]]]:
    """Merge sibling states with Pearson r >= threshold within the same parental lineage."""
    n_states = len(clusters)
    corr_mat = np.corrcoef(ref_mat)

    # Track merged groups
    merged_groups: list[list[int]] = []
    visited = set()

    for i in range(n_states):
        if i in visited:
            continue
        group = [i]
        visited.add(i)
        lineage_i = map_coarse_lineage(clusters[i])

        for j in range(i + 1, n_states):
            if j not in visited and map_coarse_lineage(clusters[j]) == lineage_i:
                if corr_mat[i, j] >= threshold:
                    group.append(j)
                    visited.add(j)
        merged_groups.append(group)

    consolidated_mat_rows: list[np.ndarray] = []
    consolidated_names: list[str] = []
    cluster_mapping: dict[str, list[str]] = {}

    for grp in merged_groups:
        grp_names = [clusters[idx] for idx in grp]
        if len(grp) == 1:
            name = grp_names[0]
            profile = ref_mat[grp[0]]
        else:
            # Create consolidated name
            lin = map_coarse_lineage(grp_names[0])
            name = f"{lin}_Consolidated_{grp_names[0][:2]}_{grp_names[-1][:2]}"
            profile = np.mean(ref_mat[grp, :], axis=0)

        consolidated_names.append(name)
        consolidated_mat_rows.append(profile)
        cluster_mapping[name] = grp_names

    consolidated_mat = np.vstack(consolidated_mat_rows)
    return consolidated_mat, consolidated_names, cluster_mapping


def run_deconvolution_pipeline(
    config: ResolutionConfig,
) -> Result[
    tuple[
        pl.DataFrame,
        pl.DataFrame,
        pl.DataFrame,
        float,
        float,
        list[str],
        list[str],
        dict[str, list[str]],
    ],
    str,
]:
    """Execute both Hierarchical BayesPrism and Consolidated BayesPrism."""
    # 1. Load inputs
    ref_df = pl.read_parquet(config.ref_path)
    bulk_df = pd.read_csv(config.hugo_tpm_path, sep="\t", index_col=0)
    bulk_df.index = bulk_df.index.astype(str).str.upper()
    bulk_df = bulk_df.groupby(level=0).mean()

    clin_df = pd.read_csv(config.hugo_clin_path, sep="\t", skiprows=4).set_index("SAMPLE_ID")

    clusters = ref_df["cluster"].to_list()
    ref_genes = [c for c in ref_df.columns if c != "cluster"]
    ref_mat = ref_df.select(ref_genes).to_numpy().astype(np.float64)

    common_genes = sorted(list(set(bulk_df.index).intersection(set(ref_genes))))
    ref_gene_indices = [ref_genes.index(g) for g in common_genes]
    sub_ref = ref_mat[:, ref_gene_indices]

    sub_bulk = bulk_df.loc[common_genes].T  # samples x genes
    sample_ids = [s for s in list(sub_bulk.index) if s in clin_df.index]
    bulk_counts = np.round(sub_bulk.loc[sample_ids].to_numpy().astype(np.float64)).astype(np.int64)

    # Calculate Condition Numbers
    kappa_original = float(np.linalg.cond(sub_ref.T))

    # Consolidate collinear clusters
    sub_ref_cons, cons_names, cluster_map = consolidate_collinear_reference(
        sub_ref, clusters, config.collinear_threshold
    )
    kappa_consolidated = float(np.linalg.cond(sub_ref_cons.T))

    print("==================================================================")
    print("REFERENCE CONDITION NUMBER ANALYSIS")
    print(f"Original Reference (12 States) Condition Number kappa: {kappa_original:8.1f}")
    print(f"Consolidated Reference ({len(cons_names)} States) Condition Number kappa: {kappa_consolidated:8.1f}")
    print(f"Condition Number Reduction: {(1.0 - kappa_consolidated / kappa_original) * 100:.1f}%")
    print("==================================================================")

    # 2. Setup Hierarchical BayesPrism
    lineages = [map_coarse_lineage(c) for c in clusters]
    ctrl = bp.GibbsControl(chain_length=config.chain_length, burn_in=config.burn_in, thinning=2)

    print("\n[INFO] Running Hierarchical BayesPrism (Coarse Lineages + Fine States)...")
    prism_hier = bp.new_prism(
        reference=sub_ref,
        cell_type_labels=lineages,      # Coarse 5 lineages
        cell_state_labels=Some(clusters),  # Fine 12 states
        mixture=bulk_counts,
        gene_names=Some(common_genes),
        bulk_names=Some(sample_ids),
    )
    match prism_hier:
        case Failure(err):
            return Failure(f"Hierarchical Prism init failed: {err}")
        case Success(prism_h_obj):
            fit_h = bp.run_prism(prism_h_obj, update_gibbs=False, gibbs_control=Some(ctrl))
            match fit_h:
                case Failure(err):
                    return Failure(f"Hierarchical Prism run failed: {err}")
                case Success(fitted_h):
                    # Tier 1 Lineage fractions
                    frac_lineage_res = bp.get_fraction(fitted_h, which_theta="first", state_or_type="type")
                    # Tier 2 Sub-state fractions
                    frac_state_res = bp.get_fraction(fitted_h, which_theta="first", state_or_type="state")
                    match (frac_lineage_res, frac_state_res):
                        case (Success(df_lin), Success(df_state)):
                            frac_hier_lineage = df_lin
                            frac_hier_state = df_state
                        case _:
                            return Failure("Failed to extract hierarchical fractions")

    # 3. Setup Consolidated BayesPrism
    print("[INFO] Running Consolidated BayesPrism (Pruned Collinear Clusters)...")
    prism_cons = bp.new_prism(
        reference=sub_ref_cons,
        cell_type_labels=cons_names,
        cell_state_labels=Nothing,
        mixture=bulk_counts,
        gene_names=Some(common_genes),
        bulk_names=Some(sample_ids),
    )
    match prism_cons:
        case Failure(err):
            return Failure(f"Consolidated Prism init failed: {err}")
        case Success(prism_c_obj):
            fit_c = bp.run_prism(prism_c_obj, update_gibbs=False, gibbs_control=Some(ctrl))
            match fit_c:
                case Failure(err):
                    return Failure(f"Consolidated Prism run failed: {err}")
                case Success(fitted_c):
                    frac_cons_res = bp.get_fraction(fitted_c, which_theta="first", state_or_type="type")
                    match frac_cons_res:
                        case Failure(err):
                            return Failure(f"Consolidated fraction extraction failed: {err}")
                        case Success(df_cons):
                            frac_consolidated = df_cons

    return Success((
        frac_hier_lineage,
        frac_hier_state,
        frac_consolidated,
        kappa_original,
        kappa_consolidated,
        sample_ids,
        clusters,
        cluster_map,
    ))


def fit_models_and_evaluate_concordance(
    frac_lin: pl.DataFrame,
    frac_state: pl.DataFrame,
    frac_cons: pl.DataFrame,
    clin_df: pd.DataFrame,
    sc_df: pl.DataFrame,
    cluster_map: dict[str, list[str]],
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Fit standardized response models and compute concordance across tiers."""
    valid_samples = frac_lin["bulk_id"].to_list()
    responses = np.array([
        1.0 if str(clin_df.loc[s, "RESPONSE"]).strip().lower() in ["complete response", "partial response", "responder"] else 0.0
        for s in valid_samples
    ], dtype=np.float64)

    effects_rows: list[dict[str, object]] = []

    # Map SC effects to lineage level
    sc_with_lineage = sc_df.with_columns(
        pl.col("cell_state").map_elements(map_coarse_lineage, return_dtype=pl.String).alias("lineage")
    )
    sc_lineage_df = sc_with_lineage.group_by("lineage").agg([
        pl.col("beta_sc").mean().alias("beta_sc_lineage"),
    ])

    # 1. Tier 1: Hierarchical Lineage Models
    lin_cols = [c for c in frac_lin.columns if c != "bulk_id"]
    for lin in lin_cols:
        x = frac_lin[lin].to_numpy().astype(np.float64)
        x_std = (x - np.mean(x)) / (np.std(x) + 1e-9)
        r_val, p_val = stats.pointbiserialr(responses, x_std)
        r_clip = np.clip(r_val, -0.999, 0.999)
        beta_z = float(2.0 * r_clip / np.sqrt(1.0 - r_clip**2 + 1e-12))

        sc_val = sc_lineage_df.filter(pl.col("lineage") == lin)["beta_sc_lineage"][0] if lin in sc_lineage_df["lineage"] else 0.0

        effects_rows.append({
            "evaluation_tier": "Tier 1: Major Lineage (Hierarchical)",
            "entity": lin,
            "beta_bulk": beta_z,
            "p_value_bulk": float(p_val),
            "beta_sc": float(sc_val),
        })

    # 2. Tier 2: Consolidated Clusters Models
    cons_cols = [c for c in frac_cons.columns if c != "bulk_id"]
    for cons_name in cons_cols:
        x = frac_cons[cons_name].to_numpy().astype(np.float64)
        x_std = (x - np.mean(x)) / (np.std(x) + 1e-9)
        r_val, p_val = stats.pointbiserialr(responses, x_std)
        r_clip = np.clip(r_val, -0.999, 0.999)
        beta_z = float(2.0 * r_clip / np.sqrt(1.0 - r_clip**2 + 1e-12))

        # Average SC beta for the constituent clusters
        constituent = cluster_map.get(cons_name, [cons_name])
        sc_sub = sc_df.filter(pl.col("cell_state").is_in(constituent))
        sc_val = float(sc_sub["beta_sc"].mean()) if not sc_sub.is_empty() else 0.0

        effects_rows.append({
            "evaluation_tier": "Tier 2: Consolidated Clusters (Pruned r >= 0.85)",
            "entity": cons_name,
            "beta_bulk": beta_z,
            "p_value_bulk": float(p_val),
            "beta_sc": sc_val,
        })

    effects_df = pl.DataFrame(effects_rows)

    # Compute concordance metrics
    metrics_rows: list[dict[str, object]] = []
    for tier in effects_df["evaluation_tier"].unique().to_list():
        sub = effects_df.filter(pl.col("evaluation_tier") == tier)
        b_sc = sub["beta_sc"].to_numpy()
        b_blk = sub["beta_bulk"].to_numpy()

        rho, p_rho = stats.spearmanr(b_sc, b_blk)
        r, p_r = stats.pearsonr(b_sc, b_blk)

        concordant = sum(
            np.sign(s) == np.sign(b) for s, b in zip(b_sc, b_blk, strict=False) if abs(s) > 0.05 and abs(b) > 0.05
        )
        total = sum(1 for s, b in zip(b_sc, b_blk, strict=False) if abs(s) > 0.05 and abs(b) > 0.05)
        rate = concordant / max(total, 1)

        metrics_rows.append({
            "evaluation_tier": tier,
            "n_entities": sub.height,
            "spearman_rho": float(rho),
            "spearman_pval": float(p_rho),
            "pearson_r": float(r),
            "concordance_rate": float(rate),
        })

    metrics_df = pl.DataFrame(metrics_rows)
    return effects_df, metrics_df


def plot_hierarchical_comparison(
    effects_df: pl.DataFrame,
    metrics_df: pl.DataFrame,
    out_path: Path,
) -> Result[None, str]:
    """Plot publication SVG comparing Tier 1 vs Tier 2 concordance."""
    scatter = (
        alt.Chart(effects_df)
        .mark_circle(size=120, opacity=0.9)
        .encode(
            x=alt.X("beta_sc:Q", title="Single-Cell Effect (Milo DA β_sc)"),
            y=alt.Y("beta_bulk:Q", title="Bulk Deconvolution Effect (Standardized β_bulk)"),
            color=alt.Color("evaluation_tier:N", legend=alt.Legend(title="Resolution Tier")),
            tooltip=["entity:N", "evaluation_tier:N", "beta_sc:Q", "beta_bulk:Q"],
        )
    )

    rule_zero = alt.Chart(pl.DataFrame({"x": [0]})).mark_rule(strokeDash=[4, 4], color="#888888").encode(x="x:Q")
    rule_y = alt.Chart(pl.DataFrame({"y": [0]})).mark_rule(strokeDash=[4, 4], color="#888888").encode(y="y:Q")

    chart = (
        (scatter + rule_zero + rule_y)
        .properties(
            title="Resolving Sibling Discrepancy: Tier 1 Lineages & Consolidated Clusters",
            width=520,
            height=400,
        )
        .configure_axis(grid=True, gridColor="#f0f0f0")
    )

    try:
        svg_str = vlc.vegalite_to_svg(chart.to_dict())
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(svg_str, encoding="utf-8")
        return Success(None)
    except Exception as exc:  # noqa: BLE001
        return Failure(f"Failed to save SVG: {exc}")


def run_pipeline(config: ResolutionConfig) -> Result[None, str]:
    """Execute Step 9 resolution pipeline."""
    match run_deconvolution_pipeline(config):
        case Failure(err):
            return Failure(err)
        case Success((frac_lin, frac_state, frac_cons, k_orig, k_cons, sample_ids, clusters, cluster_map)):
            pass

    clin_df = pd.read_csv(config.hugo_clin_path, sep="\t", skiprows=4).set_index("SAMPLE_ID")
    sc_df = pl.read_parquet(config.sc_effects_path)

    effects_df, metrics_df = fit_models_and_evaluate_concordance(
        frac_lin, frac_state, frac_cons, clin_df, sc_df, cluster_map
    )

    # Save parquets
    config.out_dir.mkdir(parents=True, exist_ok=True)
    effects_df.write_parquet(config.out_dir / "hierarchical_bayesprism_effects.parquet")
    metrics_df.write_parquet(config.out_dir / "hierarchical_concordance_metrics.parquet")

    # Save SVG
    config.results_dir.mkdir(parents=True, exist_ok=True)
    svg_path = config.results_dir / "fig6_hierarchical_resolution_comparison.svg"
    plot_hierarchical_comparison(effects_df, metrics_df, svg_path)

    print("\n==================================================================")
    print("CONCORDANCE AFTER SIBLING RESOLUTION (HIERARCHICAL & CONSOLIDATED)")
    print("==================================================================")
    for row in metrics_df.iter_rows(named=True):
        print(
            f"Tier: {row['evaluation_tier']:<45} | "
            f"Entities: {row['n_entities']:>2} | "
            f"Spearman ρ: {row['spearman_rho']:>5.3f} | "
            f"Concordance: {row['concordance_rate']*100:>5.1f}%"
        )
    print("==================================================================")

    print("\nDetailed Consolidated Entities:")
    for row in effects_df.filter(pl.col("evaluation_tier").str.contains("Consolidated")).iter_rows(named=True):
        print(f"  {row['entity']:<38} -> Bulk beta_z={row['beta_bulk']:>6.2f} | SC beta={row['beta_sc']:>6.2f}")
    print("==================================================================")

    return Success(None)


def main() -> None:
    config = parse_args()
    match run_pipeline(config):
        case Success(_):
            sys.exit(0)
        case Failure(err):
            print(f"[ERROR] {err}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
