#!/usr/bin/env python3
"""Step 02: Integrate Bagaev TME Subtypes with Outcome-Blind Consensus-Meta42 & SC9 Deconvolution.

Uses strictly outcome-blind `raw_instaprism_Consensus-Meta42_on_PanCancer-Combined.parquet`
and pure `SC9` `deconv_fractions_*.parquet` matrices (zero outcome-projected fraction leakage):
1. Enrichment of Consensus-Meta42 states (M01..M42) across Bagaev IE, IE/F, F, and D subtypes.
2. Within-subtype responder stratification via strictly patient-grouped out-of-fold Meta42 scores.
3. Head-to-head 5x5 patient-grouped StratifiedGroupKFold CV: Bagaev 4-Subtype vs 29-Fges vs Raw Meta42 vs Combined.
4. Direct cell-state vs Bagaev subtype correlations (raw, CLR & effective-information Fisher-z meta).
5. Outcome-blind compositional CLR clustering across K=2..8 (reporting K=2, K=3 & K=4) with
   external single-cell Milo logFC cluster ranking and cohort-stratified Cox PH evaluation.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Success

from tme_datasets import (
    BAGAEV_MFP_FUNCTIONAL_GROUPS,
    BAGAEV_MFP_SIGNATURE_ORDER,
    get_data_paths,
)
from tme_response import (
    CANONICAL_SUBTYPE_ORDER,
    DeconvClusterConfig,
    analyze_subtype_deconvolution_enrichment,
    benchmark_bagaev_vs_deconv_cv,
    cluster_deconvoluted_samples,
    compare_partition_concordance,
    compute_2d_feature_embedding,
    compute_deconv_subtype_correlations,
    compute_state_vs_fges_correlations,
    evaluate_subtype_response_associations,
    evaluate_within_subtype_deconv_predictions,
    fit_univariable_cox,
    transform_compositional_fractions,
)


def _load_and_join_frac_table(
    frac_path: Path,
    subtypes_df: pl.DataFrame,
) -> tuple[pl.DataFrame, list[str]]:
    """Load fraction Parquet and join 1-to-1 with `subtypes_df` regardless of ID prefixing."""
    frac_df = pl.read_parquet(frac_path)
    exclude_cols = {
        "sample_id", "cohort", "cohort_id", "cancer_type",
        "response", "response_binary", "dataset_id",
    }
    state_cols = [
        c for c in frac_df.columns
        if c not in exclude_cols and frac_df[c].dtype in (pl.Float32, pl.Float64)
    ]
    sample_preview = str(frac_df["sample_id"][0]) if frac_df.height > 0 else ""
    frac_norm = (
        frac_df.with_columns(pl.format("{}__{}", pl.col("cohort"), pl.col("sample_id")).alias("sample_id"))
        if "__" not in sample_preview and "cohort" in frac_df.columns
        else frac_df
    )

    joined = subtypes_df.join(frac_norm.select(["sample_id", *state_cols]), on="sample_id", how="inner")
    if joined.height < 50 and "raw_sample_id" in subtypes_df.columns:
        joined = subtypes_df.join(
            frac_df.select([pl.col("sample_id").alias("raw_sample_id"), *state_cols]),
            on="raw_sample_id",
            how="inner",
        )
    return joined, state_cols


def _short_bio_tag(state_id: str, state_label: str) -> str:
    """Condense 'M13: CXCL13+PDCD1+ Tex' into 'CXCL13+ Tex' for clean plot legends."""
    raw = state_label.split(":", 1)[1].strip() if ":" in state_label else state_label.strip()
    parts = raw.split()
    if len(parts) >= 2 and "+" in parts[0]:
        return f"{parts[0].split('+')[0]}+ {' '.join(parts[1:])}"
    return raw if raw else state_id


def _build_cluster_bio_labels(
    centroids_df: pl.DataFrame,
    state_label_map: dict[str, str],
) -> dict[str, str]:
    """Derive concise biological labels for clusters DC1..DCK from top CLR-enriched states."""
    return {
        str(cid): (
            f"{cid}: "
            + " / ".join(
                _short_bio_tag(str(st), state_label_map.get(str(st), str(st)))
                for st in (
                    centroids_df.filter(pl.col("cluster_id") == cid)
                    .sort("mean_clr_z", descending=True)
                    .head(2)["cell_state"]
                    .to_list()
                )
            )
        )
        for cid in sorted(centroids_df["cluster_id"].unique().to_list())
    }


def _run_cluster_eval_across_k_and_variants(
    joined: pl.DataFrame,
    meta_cols: list[str],
    state_label_map: dict[str, str],
    state_milo_map: dict[str, float],
    bagaev_dir: Path,
) -> None:
    """Cluster outcome-blind Meta42 fractions (K=2..8) and evaluate ORR, stratified Cox & concordance."""
    cfg_pan_ward = DeconvClusterConfig(
        k_values=(2, 3, 4, 5, 6, 7, 8),
        target_k=4,
        min_valid_k=2,
        transform="clr",
        algorithm="aitchison_ward",
        center_by_cohort=False,
        seed=42,
    )
    cfg_pan_louv = DeconvClusterConfig(
        k_values=(2, 3, 4, 5, 6, 7, 8),
        target_k=4,
        min_valid_k=2,
        transform="clr",
        algorithm="clr_louvain",
        center_by_cohort=False,
        seed=42,
    )
    cfg_coh_ward = DeconvClusterConfig(
        k_values=(2, 3, 4, 5, 6, 7, 8),
        target_k=4,
        min_valid_k=2,
        transform="clr",
        algorithm="aitchison_ward",
        center_by_cohort=True,
        seed=42,
    )
    cfg_arcsin_ward = DeconvClusterConfig(
        k_values=(4,),
        target_k=4,
        min_valid_k=2,
        transform="arcsin_sqrt",
        algorithm="aitchison_ward",
        center_by_cohort=False,
        seed=42,
    )

    res_pan_ward = cluster_deconvoluted_samples(joined, meta_cols, config=cfg_pan_ward).unwrap()
    res_pan_louv = cluster_deconvoluted_samples(joined, meta_cols, config=cfg_pan_louv).unwrap()
    res_coh_ward = cluster_deconvoluted_samples(joined, meta_cols, config=cfg_coh_ward).unwrap()
    res_arcsin = cluster_deconvoluted_samples(joined, meta_cols, config=cfg_arcsin_ward).unwrap()

    clr_df = transform_compositional_fractions(
        joined, meta_cols, method="clr", pseudocount=1e-6, standardize_cols=True
    ).unwrap()
    _, _, umap_x, umap_y = compute_2d_feature_embedding(clr_df, meta_cols, metric="correlation", seed=42)

    bio_labels_k4 = _build_cluster_bio_labels(res_pan_ward.centroids_k4_df, state_label_map)
    bio_labels_coh = _build_cluster_bio_labels(res_coh_ward.centroids_k4_df, state_label_map)

    centroids_out = (
        res_pan_ward.centroids_k4_df.with_columns(
            [
                pl.col("cluster_id").replace_strict(bio_labels_k4, default=pl.col("cluster_id")).alias("cluster_label"),
                pl.col("cell_state").replace_strict(state_label_map, default=pl.col("cell_state")).alias("state_label"),
                pl.col("cell_state").replace_strict(state_milo_map, default=0.0).alias("state_milo_logfc"),
            ]
        )
        .sort(["cluster_id", "mean_clr_z"], descending=[False, True])
    )
    centroids_out.write_parquet(bagaev_dir / "iatlas_deconv_cluster_centroids_k4.parquet")
    centroids_out.write_csv(bagaev_dir / "iatlas_deconv_cluster_centroids_k4.tsv", separator="\t")

    pan_rename = {
        "cluster_k2": "deconv_cluster_pan_ward_k2",
        "cluster_k3": "deconv_cluster_pan_ward_k3",
        "cluster_k4": "deconv_cluster_pan_ward_k4",
        "cluster_opt": "deconv_cluster_pan_ward_opt",
        "cluster_opt_multi": "deconv_cluster_pan_ward_opt_multi",
        "pca_1": "clr_pca_1",
        "pca_2": "clr_pca_2",
    }
    pan_assign = res_pan_ward.assignments_df.rename(
        {k: v for k, v in pan_rename.items() if k in res_pan_ward.assignments_df.columns}
    )
    louv_assign = res_pan_louv.assignments_df.select(
        [pl.col("sample_id"), pl.col("cluster_k4").alias("deconv_cluster_pan_louvain_k4")]
    )
    coh_assign = res_coh_ward.assignments_df.select(
        [
            pl.col("sample_id"),
            pl.col("cluster_k4").alias("deconv_cluster_coh_ward_k4"),
            pl.col("cluster_opt").alias("deconv_cluster_coh_ward_opt"),
            pl.col("cluster_opt_multi").alias("deconv_cluster_coh_ward_opt_multi"),
        ]
    )
    arc_assign = res_arcsin.assignments_df.select(
        [pl.col("sample_id"), pl.col("cluster_k4").alias("deconv_cluster_arcsin_ward_k4")]
    )

    sample_clusters_df = (
        joined.select(
            [
                "sample_id", "patient_id", "cohort", "cancer_type",
                "is_baseline_pre", "is_primary_baseline_patient",
                "response_binary", "response_recist", "os_months", "os_event",
                "pfs_months", "pfs_event", "bagaev_subtype", "louvain_subtype",
                "prob_IE", "prob_IEF", "prob_F", "prob_D",
                "antitumor_immune_score", "stromal_fibrotic_score", "oof_raw_meta42_z",
            ]
        )
        .join(pan_assign, on="sample_id", how="left")
        .join(louv_assign, on="sample_id", how="left")
        .join(coh_assign, on="sample_id", how="left")
        .join(arc_assign, on="sample_id", how="left")
        .with_columns(
            [
                pl.Series("clr_umap_1", umap_x),
                pl.Series("clr_umap_2", umap_y),
                pl.col("deconv_cluster_pan_ward_k4")
                .replace_strict(bio_labels_k4, default=pl.col("deconv_cluster_pan_ward_k4"))
                .alias("deconv_cluster_pan_ward_k4_label"),
                pl.col("deconv_cluster_coh_ward_k4")
                .replace_strict(bio_labels_coh, default=pl.col("deconv_cluster_coh_ward_k4"))
                .alias("deconv_cluster_coh_ward_k4_label"),
            ]
        )
    )
    sample_clusters_df.write_parquet(bagaev_dir / "iatlas_deconv_cluster_samples.parquet")
    sample_clusters_df.write_csv(bagaev_dir / "iatlas_deconv_cluster_samples.tsv", separator="\t")

    # K=2..8 sweep combining internal validity & baseline patient ORR separation
    sweep_tables: list[pl.DataFrame] = []
    for cfg_label, res_obj in (
        ("PanCancer_CLR_Ward", res_pan_ward),
        ("PanCancer_CLR_Louvain", res_pan_louv),
        ("CohortCentered_CLR_Ward", res_coh_ward),
    ):
        eval_base = joined.select(["sample_id", "cohort", "is_primary_baseline_patient", "response_binary"]).join(
            res_obj.assignments_df, on="sample_id", how="inner"
        )
        extra_rows: list[dict[str, float | int | str]] = []
        for r_k in res_obj.k_sweep_df.iter_rows(named=True):
            k_val = int(r_k["k"])
            k_col = f"cluster_k{k_val}"
            rec: dict[str, float | int | str] = {
                "method": cfg_label,
                "k": k_val,
                "optimal_k": int(res_obj.optimal_k),
                "optimal_k_multi": int(res_obj.optimal_k_multi),
            }
            for scope_name, sc_df in (
                ("baseline", eval_base.filter(pl.col("is_primary_baseline_patient"))),
                ("all", eval_base),
            ):
                rr_res = evaluate_subtype_response_associations(
                    sc_df, subtype_col=k_col, cohort_col=Some("cohort")
                )
                if isinstance(rr_res, Success):
                    rr_df = rr_res.unwrap()
                    orrs = rr_df["response_rate"].to_numpy()
                    rec[f"{scope_name}_min_orr"] = float(np.min(orrs))
                    rec[f"{scope_name}_max_orr"] = float(np.max(orrs))
                    rec[f"{scope_name}_orr_spread"] = float(rr_df["orr_spread"][0])
                    rec[f"{scope_name}_global_chi2_p"] = float(rr_df["global_chi2_p"][0])
                    rec[f"{scope_name}_cmh_omnibus_p"] = float(rr_df["cmh_omnibus_p"][0])
                else:
                    rec[f"{scope_name}_min_orr"] = float("nan")
                    rec[f"{scope_name}_max_orr"] = float("nan")
                    rec[f"{scope_name}_orr_spread"] = float("nan")
                    rec[f"{scope_name}_global_chi2_p"] = 1.0
                    rec[f"{scope_name}_cmh_omnibus_p"] = 1.0
            extra_rows.append(rec)

        sweep_tables.append(res_obj.k_sweep_df.join(pl.DataFrame(extra_rows), on="k", how="inner"))

    all_sweep_df = pl.concat(sweep_tables, how="diagonal_relaxed")
    all_sweep_df.write_parquet(bagaev_dir / "iatlas_deconv_cluster_k_sweep.parquet")
    all_sweep_df.write_csv(bagaev_dir / "iatlas_deconv_cluster_k_sweep.tsv", separator="\t")

    # Per-cluster response rates transparently covering K=2 (global opt), K=3 (multi opt) & K=4
    method_specs: list[tuple[str, str, dict[str, str]]] = [
        ("deconv_cluster_pan_ward_k2", "Meta42_PanCancer_CLR_Ward_K2", {}),
        ("deconv_cluster_pan_ward_k3", "Meta42_PanCancer_CLR_Ward_K3", {}),
        ("deconv_cluster_pan_ward_k4", "Meta42_PanCancer_CLR_Ward_K4", bio_labels_k4),
        ("deconv_cluster_pan_ward_opt", f"Meta42_PanCancer_CLR_Ward_OptK{res_pan_ward.optimal_k}", {}),
        ("deconv_cluster_pan_ward_opt_multi", f"Meta42_PanCancer_CLR_Ward_MultiOptK{res_pan_ward.optimal_k_multi}", {}),
        ("deconv_cluster_pan_louvain_k4", "Meta42_PanCancer_CLR_Louvain_K4", {}),
        ("deconv_cluster_coh_ward_k4", "Meta42_CohortCentered_CLR_Ward_K4", bio_labels_coh),
        ("deconv_cluster_coh_ward_opt", f"Meta42_CohortCentered_CLR_Ward_OptK{res_coh_ward.optimal_k}", {}),
        ("deconv_cluster_arcsin_ward_k4", "Meta42_ArcsinSqrt_Ward_K4", {}),
        ("bagaev_subtype", "Bagaev_4Subtype_TCGA_kNN_Ref", {}),
    ]

    resp_tables: list[pl.DataFrame] = []
    cohorts_unique = sorted(sample_clusters_df["cohort"].unique().to_list())
    for scope_name, sc_df in (
        ("Baseline_Pre", sample_clusters_df.filter(pl.col("is_primary_baseline_patient"))),
        ("All_Samples", sample_clusters_df),
    ):
        for col_name, m_label, lbl_map in method_specs:
            if col_name not in sc_df.columns:
                continue
            label_expr = (
                pl.col("subtype").replace_strict(lbl_map, default=pl.col("subtype"))
                if lbl_map
                else pl.col("subtype")
            ).alias("cluster_label")
            pan_res = evaluate_subtype_response_associations(sc_df, subtype_col=col_name, cohort_col=Some("cohort"))
            if isinstance(pan_res, Success):
                resp_tables.append(
                    pan_res.unwrap().with_columns(
                        [
                            pl.lit("PanCancer-iAtlas").alias("cohort"),
                            pl.lit(scope_name).alias("sample_scope"),
                            pl.lit(m_label).alias("method"),
                            label_expr,
                        ]
                    )
                )
            if col_name in ("deconv_cluster_pan_ward_k4", "deconv_cluster_coh_ward_k4"):
                for cid in cohorts_unique:
                    c_df = sc_df.filter(pl.col("cohort") == cid)
                    c_res = evaluate_subtype_response_associations(c_df, subtype_col=col_name, cohort_col=Nothing)
                    if isinstance(c_res, Success):
                        resp_tables.append(
                            c_res.unwrap().with_columns(
                                [
                                    pl.lit(cid).alias("cohort"),
                                    pl.lit(scope_name).alias("sample_scope"),
                                    pl.lit(m_label).alias("method"),
                                    label_expr,
                                ]
                            )
                        )

    all_resp_df = pl.concat(resp_tables, how="diagonal_relaxed")
    all_resp_df.write_parquet(bagaev_dir / "iatlas_deconv_cluster_response_rates.parquet")
    all_resp_df.write_csv(bagaev_dir / "iatlas_deconv_cluster_response_rates.tsv", separator="\t")

    cont_df = (
        sample_clusters_df.group_by(["deconv_cluster_pan_ward_k4", "deconv_cluster_pan_ward_k4_label", "bagaev_subtype"])
        .agg(
            [
                pl.len().alias("n_samples"),
                pl.col("response_binary").filter(pl.col("response_binary").cast(pl.Float64).is_in([0.0, 1.0])).mean().alias("observed_orr"),
            ]
        )
        .with_columns(
            (100.0 * pl.col("n_samples") / pl.col("n_samples").sum().over("deconv_cluster_pan_ward_k4")).alias("pct_within_deconv_cluster")
        )
        .sort(["deconv_cluster_pan_ward_k4", "bagaev_subtype"])
    )
    cont_df.write_parquet(bagaev_dir / "iatlas_deconv_clusters_vs_bagaev_contingency.parquet")
    cont_df.write_csv(bagaev_dir / "iatlas_deconv_clusters_vs_bagaev_contingency.tsv", separator="\t")

    conc_rows: list[dict[str, str | float]] = []
    for col_name, m_label, _ in method_specs:
        if col_name == "bagaev_subtype" or col_name not in sample_clusters_df.columns:
            continue
        for ref_col, ref_name in (("bagaev_subtype", "Bagaev_TCGA_kNN"), ("louvain_subtype", "Bagaev_DeNovo_Louvain")):
            c_met = compare_partition_concordance(
                sample_clusters_df[col_name].cast(pl.String).to_list(),
                sample_clusters_df[ref_col].cast(pl.String).to_list(),
            )
            conc_rows.append({"deconv_method": m_label, "bagaev_reference": ref_name, **c_met})
    conc_df = pl.DataFrame(conc_rows)
    conc_df.write_parquet(bagaev_dir / "iatlas_deconv_vs_bagaev_concordance.parquet")
    conc_df.write_csv(bagaev_dir / "iatlas_deconv_vs_bagaev_concordance.tsv", separator="\t")

    # External single-cell Milo logFC cluster ranking (zero bulk outcome leakage)
    cluster_imm_rank = (
        centroids_out.with_columns((pl.col("mean_clr_z") * pl.col("state_milo_logfc")).alias("imm_score_contrib"))
        .group_by("cluster_id")
        .agg(pl.col("imm_score_contrib").sum().alias("immune_axis_score"))
        .sort("immune_axis_score", descending=True)
    )
    top_bio_cluster = str(cluster_imm_rank["cluster_id"][0])
    top_two_bio_clusters = [str(c) for c in cluster_imm_rank["cluster_id"].head(2).to_list()]

    surv_rows: list[dict[str, str | float | int]] = []
    pre_s = sample_clusters_df.filter(pl.col("is_primary_baseline_patient"))
    for ep, t_col, e_col in (("OS", "os_months", "os_event"), ("PFS", "pfs_months", "pfs_event")):
        for sc_name, sc_sub in [("PanCancer-iAtlas", pre_s)] + [
            (cid, pre_s.filter(pl.col("cohort") == cid)) for cid in cohorts_unique
        ]:
            v_s = sc_sub.filter(
                pl.col(t_col).is_not_null()
                & pl.col(t_col).is_finite()
                & (pl.col(t_col) > 0.0)
                & pl.col(e_col).cast(pl.Float64).is_in([0.0, 1.0])
            )
            if v_s.height < 15:
                continue
            times = v_s[t_col].to_numpy().astype(np.float64)
            events = v_s[e_col].to_numpy().astype(np.int64)
            if int(np.sum(events)) < 5:
                continue
            strata_arg = Some(v_s["cohort"].cast(pl.String).to_list()) if sc_name == "PanCancer-iAtlas" else Nothing
            for contrast_name, ind in (
                (
                    f"TopImmune_{top_bio_cluster}_vs_Rest",
                    (v_s["deconv_cluster_pan_ward_k4"] == top_bio_cluster).to_numpy().astype(np.float64),
                ),
                (
                    "ImmuneHighHalf_vs_LowHalf",
                    v_s["deconv_cluster_pan_ward_k4"].is_in(top_two_bio_clusters).to_numpy().astype(np.float64),
                ),
            ):
                cx_r = fit_univariable_cox(
                    times, events, ind, standardize_score=False, strata=strata_arg
                )
                if isinstance(cx_r, Success):
                    cx = cx_r.unwrap()
                    surv_rows.append(
                        {
                            "cohort": sc_name,
                            "endpoint": ep,
                            "contrast": contrast_name,
                            "n_samples": int(cx.n_samples),
                            "n_events": int(cx.n_events),
                            "hazard_ratio": float(cx.hazard_ratio),
                            "hr_ci95_low": float(cx.hr_ci_lower),
                            "hr_ci95_high": float(cx.hr_ci_upper),
                            "wald_p": float(cx.p_value),
                        }
                    )
    if surv_rows:
        surv_df = pl.DataFrame(surv_rows)
        surv_df.write_parquet(bagaev_dir / "iatlas_deconv_cluster_survival_cox.parquet")
        surv_df.write_csv(bagaev_dir / "iatlas_deconv_cluster_survival_cox.tsv", separator="\t")


def _run_multires_deconv_benchmark(
    subtypes_df: pl.DataFrame,
    bench_3c_dir: Path,
    sc9_deconv_dir: Path,
    bagaev_dir: Path,
) -> None:
    """Benchmark strictly outcome-blind deconvolution references across SC9 res0.25..2.0 & Meta42 Raw."""
    candidates: list[tuple[str, str, Path]] = [
        (
            "Consensus-Meta42_RawInstaPrism",
            "Meta42 (Raw InstaPrism)",
            bench_3c_dir / "fractions" / "raw_instaprism_Consensus-Meta42_on_PanCancer-Combined.parquet",
        ),
    ]
    for fam, fam_disp in (("all_9_combined", "SC9 All-9"), ("solid_tumor_quartet", "SC9 Quartet")):
        for res_str in ("0.25", "0.5", "0.75", "1.0", "1.25", "1.5", "1.75", "2.0"):
            candidates.append(
                (
                    f"SC9_{fam}_res{res_str}",
                    f"{fam_disp} r={res_str}",
                    sc9_deconv_dir / f"deconv_fractions_{fam}_res{res_str}.parquet",
                )
            )

    rows: list[dict[str, str | int | float]] = []
    bag_base = evaluate_subtype_response_associations(
        subtypes_df.filter(pl.col("is_primary_baseline_patient")),
        subtype_col="bagaev_subtype",
        cohort_col=Some("cohort"),
    )
    if isinstance(bag_base, Success):
        bdf = bag_base.unwrap()
        rows.append(
            {
                "reference_id": "Bagaev_4Subtype_TCGA_kNN",
                "display_name": "Bagaev 4-Subtype (Ref)",
                "clustering_regime": "Bulk_29Fges_kNN",
                "n_states": 29,
                "k_clusters": 4,
                "optimal_k": 4,
                "optimal_k_multi": 4,
                "n_evaluable_baseline": int(bdf["n_samples"].sum()),
                "silhouette_euclidean": float("nan"),
                "ari_vs_bagaev_knn": 1.0,
                "baseline_orr_min": float(bdf["response_rate"].min()),
                "baseline_orr_max": float(bdf["response_rate"].max()),
                "baseline_orr_spread": float(bdf["orr_spread"][0]),
                "baseline_global_chi2_p": float(bdf["global_chi2_p"][0]),
                "baseline_cmh_omnibus_p": float(bdf["cmh_omnibus_p"][0]),
            }
        )

    for ref_id, disp_name, fpath in candidates:
        if not fpath.is_file():
            continue
        j_df, st_cols = _load_and_join_frac_table(fpath, subtypes_df)
        if j_df.height < 100 or len(st_cols) < 4:
            continue
        for reg_name, center_flag in (("PanCancer_CLR_Ward", False), ("CohortCentered_CLR_Ward", True)):
            cfg = DeconvClusterConfig(
                k_values=(2, 3, 4),
                target_k=4,
                min_valid_k=2,
                transform="clr",
                algorithm="aitchison_ward",
                center_by_cohort=center_flag,
                seed=42,
            )
            cl_res = cluster_deconvoluted_samples(j_df, st_cols, config=cfg)
            if isinstance(cl_res, Failure):
                continue
            cl_out = cl_res.unwrap()
            sil_euc = float(cl_out.k_sweep_df.filter(pl.col("k") == 4)["silhouette_euclidean"][0])
            ev_df = j_df.select(
                ["sample_id", "cohort", "is_primary_baseline_patient", "response_binary", "bagaev_subtype"]
            ).join(cl_out.assignments_df.select(["sample_id", "cluster_k4"]), on="sample_id", how="inner")
            ari_bg = compare_partition_concordance(
                ev_df["cluster_k4"].cast(pl.String).to_list(),
                ev_df["bagaev_subtype"].cast(pl.String).to_list(),
            )["ari"]
            base_sub = ev_df.filter(
                pl.col("is_primary_baseline_patient") & pl.col("response_binary").cast(pl.Float64).is_in([0.0, 1.0])
            )
            rr_res = evaluate_subtype_response_associations(base_sub, subtype_col="cluster_k4", cohort_col=Some("cohort"))
            if isinstance(rr_res, Success):
                rdf = rr_res.unwrap()
                rows.append(
                    {
                        "reference_id": ref_id,
                        "display_name": disp_name,
                        "clustering_regime": reg_name,
                        "n_states": len(st_cols),
                        "k_clusters": 4,
                        "optimal_k": int(cl_out.optimal_k),
                        "optimal_k_multi": int(cl_out.optimal_k_multi),
                        "n_evaluable_baseline": int(rdf["n_samples"].sum()),
                        "silhouette_euclidean": sil_euc,
                        "ari_vs_bagaev_knn": float(ari_bg),
                        "baseline_orr_min": float(rdf["response_rate"].min()),
                        "baseline_orr_max": float(rdf["response_rate"].max()),
                        "baseline_orr_spread": float(rdf["orr_spread"][0]),
                        "baseline_global_chi2_p": float(rdf["global_chi2_p"][0]),
                        "baseline_cmh_omnibus_p": float(rdf["cmh_omnibus_p"][0]),
                    }
                )

    bench_df = pl.DataFrame(rows)
    bench_df.write_parquet(bagaev_dir / "iatlas_multires_deconv_cluster_benchmark.parquet")
    bench_df.write_csv(bagaev_dir / "iatlas_multires_deconv_cluster_benchmark.tsv", separator="\t")


def run_step02(
    bagaev_dir: Path,
    bench_3c_dir: Path,
    sc9_deconv_dir: Path,
) -> None:
    """Execute Step 02 integration of iAtlas Bagaev subtypes with outcome-blind deconvolution."""
    subtypes_path = bagaev_dir / "iatlas_bagaev_sample_subtypes.parquet"
    if not subtypes_path.is_file():
        raise FileNotFoundError(f"Step 01 output not found at {subtypes_path}")
    subtypes_df = pl.read_parquet(subtypes_path)

    # Strictly use zero-leakage raw InstaPrism Consensus-Meta42 fractions
    raw_frac_path = bench_3c_dir / "fractions" / "raw_instaprism_Consensus-Meta42_on_PanCancer-Combined.parquet"
    hub_path = bench_3c_dir / "consensus_metastate_all42_hub_summary.parquet"
    deg_path = bench_3c_dir / "consensus_metastate_all42_top_degs.parquet"

    if not raw_frac_path.is_file():
        raise FileNotFoundError(f"Missing zero-leakage InstaPrism fractions at {raw_frac_path}")

    raw_frac_df = pl.read_parquet(raw_frac_path)
    meta_cols = [f"M{i:02d}" for i in range(1, 43) if f"M{i:02d}" in raw_frac_df.columns]
    joined_base = subtypes_df.join(raw_frac_df.select(["sample_id", *meta_cols]), on="sample_id", how="inner")

    # Zero-leakage patient-grouped 5x5 CV (returns pan-cancer, per-cohort, and OOF sample scores)
    cv_input = joined_base.filter(pl.col("is_primary_baseline_patient"))
    cv_res = benchmark_bagaev_vs_deconv_cv(
        cv_input,
        state_cols=meta_cols,
        signature_cols=BAGAEV_MFP_SIGNATURE_ORDER,
        subtype_col="bagaev_subtype",
        response_col="response_binary",
        group_col=Some("patient_id"),
        cohort_col=Some("cohort"),
    )
    if isinstance(cv_res, Failure):
        raise RuntimeError(f"CV benchmark failed: {cv_res.failure()}")
    pan_cv_df, per_cohort_cv_df, oof_preds_df = cv_res.unwrap()

    pan_cv_df.write_parquet(bagaev_dir / "iatlas_bagaev_vs_deconv_cv_benchmark.parquet")
    pan_cv_df.write_csv(bagaev_dir / "iatlas_bagaev_vs_deconv_cv_benchmark.tsv", separator="\t")
    per_cohort_cv_df.write_parquet(bagaev_dir / "iatlas_bagaev_vs_deconv_per_cohort_cv.parquet")
    per_cohort_cv_df.write_csv(bagaev_dir / "iatlas_bagaev_vs_deconv_per_cohort_cv.tsv", separator="\t")

    joined = joined_base.join(oof_preds_df, on="sample_id", how="left")
    joined.write_parquet(bagaev_dir / "iatlas_bagaev_joined_meta42_oof.parquet")
    print(f"Joined {joined.height} iAtlas samples ({cv_input.height} unique baseline patients) with zero-leakage InstaPrism Meta42.")

    # Annotate Meta42 states with external single-cell Milo logFC & DEG markers
    state_label_map: dict[str, str] = {c: c for c in meta_cols}
    state_milo_map: dict[str, float] = {c: 0.0 for c in meta_cols}
    hub_df = pl.read_parquet(hub_path) if hub_path.is_file() else pl.DataFrame()
    if not hub_df.is_empty():
        for r in hub_df.iter_rows(named=True):
            cs = str(r.get("cell_state") or "")
            if cs:
                state_label_map[cs] = str(r.get("state_label") or cs)
                state_milo_map[cs] = float(r.get("state_milo_logfc") or 0.0)

    top_markers_maybe: Maybe[pl.DataFrame] = Nothing
    if deg_path.is_file():
        deg_df = pl.read_parquet(deg_path)
        st_c = next((c for c in ("meta_state", "state_id", "cell_state") if c in deg_df.columns), "")
        gn_c = next((c for c in ("gene_symbol", "gene_name", "gene") if c in deg_df.columns), "")
        if st_c and gn_c:
            top_markers_maybe = Some(
                deg_df.group_by(st_c)
                .agg(pl.col(gn_c).head(5).str.join(", ").alias("top_deg_markers"))
                .rename({st_c: "cell_state"})
            )

    enr_res = analyze_subtype_deconvolution_enrichment(joined, meta_cols, subtype_col="bagaev_subtype")
    if isinstance(enr_res, Failure):
        raise RuntimeError(f"Meta42 enrichment failed: {enr_res.failure()}")
    enr_df = enr_res.unwrap()

    if not hub_df.is_empty():
        id_col = next((c for c in ("cell_state", "meta_state", "state_id") if c in hub_df.columns), "")
        if id_col:
            keep_cols = [c for c in ("state_label", "enrichment_group", "state_milo_logfc") if c in hub_df.columns and c not in enr_df.columns]
            if keep_cols:
                enr_df = enr_df.join(
                    hub_df.select([pl.col(id_col).alias("cell_state"), *keep_cols]),
                    on="cell_state",
                    how="left",
                )
    if isinstance(top_markers_maybe, Some):
        enr_df = enr_df.join(top_markers_maybe.unwrap(), on="cell_state", how="left")

    enr_df.write_parquet(bagaev_dir / "iatlas_meta42_enrichment_by_bagaev_subtype.parquet")
    enr_df.write_csv(bagaev_dir / "iatlas_meta42_enrichment_by_bagaev_subtype.tsv", separator="\t")

    # Pre-compute top-16 subtype-differentiated heatmap row Z-scores for Step 04 Fig 3 Panel A
    top_enr = enr_df.head(16)
    hm_rows: list[dict[str, str | float | int]] = []
    for row_idx, r in enumerate(top_enr.iter_rows(named=True)):
        vals = [float(r["mean_IE"]), float(r["mean_IEF"]), float(r["mean_F"]), float(r["mean_D"])]
        mu = sum(vals) / 4.0
        sd = max(1e-8, (sum((v - mu) ** 2 for v in vals) / 4.0) ** 0.5)
        markers = str(r.get("top_deg_markers") or "")
        short_m = ", ".join(markers.split(", ")[:2]) if markers else ""
        st_label = f"{r['cell_state']} ({short_m})" if short_m else str(r["cell_state"])
        for st_name, raw_v in zip(CANONICAL_SUBTYPE_ORDER, vals, strict=True):
            hm_rows.append(
                {
                    "order_idx": row_idx,
                    "cell_state": str(r["cell_state"]),
                    "state_label": st_label,
                    "subtype": st_name,
                    "rel_z": (raw_v - mu) / sd,
                }
            )
    pl.DataFrame(hm_rows).write_parquet(bagaev_dir / "iatlas_meta42_enrichment_heatmap_z.parquet")

    # Direct correlations of Consensus-Meta42 states with Bagaev subtypes & 29 Fges
    sub_corr_res = compute_deconv_subtype_correlations(
        joined,
        meta_cols,
        subtype_col="bagaev_subtype",
        louvain_col=Some("louvain_subtype"),
        cohort_col=Some("cohort"),
    )
    if isinstance(sub_corr_res, Failure):
        raise RuntimeError(f"Subtype correlation failed: {sub_corr_res.failure()}")
    sub_corr_df = sub_corr_res.unwrap().with_columns(
        pl.lit(joined.height).alias("n_evaluated_samples")
    )
    if not hub_df.is_empty() and "cell_state" in hub_df.columns:
        hub_cols = [c for c in ("state_label", "enrichment_group", "state_milo_logfc") if c in hub_df.columns]
        if hub_cols:
            sub_corr_df = sub_corr_df.join(hub_df.select(["cell_state", *hub_cols]), on="cell_state", how="left")
    if isinstance(top_markers_maybe, Some):
        sub_corr_df = sub_corr_df.join(top_markers_maybe.unwrap(), on="cell_state", how="left")

    sub_corr_df.write_parquet(bagaev_dir / "iatlas_deconv_subtype_correlations.parquet")
    sub_corr_df.write_csv(bagaev_dir / "iatlas_deconv_subtype_correlations.tsv", separator="\t")

    fges_raw_res = compute_state_vs_fges_correlations(joined, meta_cols, BAGAEV_MFP_SIGNATURE_ORDER, use_clr=False)
    fges_clr_res = compute_state_vs_fges_correlations(joined, meta_cols, BAGAEV_MFP_SIGNATURE_ORDER, use_clr=True)
    if isinstance(fges_raw_res, Success):
        sig_to_mod = {sig: grp for grp, sigs in BAGAEV_MFP_FUNCTIONAL_GROUPS.items() for sig in sigs}
        fges_df = fges_raw_res.unwrap()
        if isinstance(fges_clr_res, Success):
            fges_df = fges_df.join(
                fges_clr_res.unwrap().select(
                    ["cell_state", "fges_signature", pl.col("spearman_rho").alias("spearman_rho_clr")]
                ),
                on=["cell_state", "fges_signature"],
                how="left",
            )
        fges_df = fges_df.with_columns(
            [
                pl.col("fges_signature").replace_strict(sig_to_mod, default="Other").alias("fges_module"),
                pl.col("cell_state").replace_strict(state_label_map, default=pl.col("cell_state")).alias("state_label"),
            ]
        )
        fges_df.write_parquet(bagaev_dir / "iatlas_meta42_vs_29fges_correlations.parquet")
        fges_df.write_csv(bagaev_dir / "iatlas_meta42_vs_29fges_correlations.tsv", separator="\t")

    # Within-subtype out-of-fold discrimination using strictly zero-leakage Meta42 CV predictions
    within_res = evaluate_within_subtype_deconv_predictions(
        joined.filter(pl.col("is_primary_baseline_patient")),
        score_col="oof_raw_meta42_z",
        subtype_col="bagaev_subtype",
        response_col="response_binary",
    )
    if isinstance(within_res, Failure):
        raise RuntimeError(f"Within-subtype evaluation failed: {within_res.failure()}")
    within_df = within_res.unwrap()
    within_df.write_parquet(bagaev_dir / "iatlas_within_subtype_deconv_auc.parquet")
    within_df.write_csv(bagaev_dir / "iatlas_within_subtype_deconv_auc.tsv", separator="\t")

    _run_cluster_eval_across_k_and_variants(
        joined=joined,
        meta_cols=meta_cols,
        state_label_map=state_label_map,
        state_milo_map=state_milo_map,
        bagaev_dir=bagaev_dir,
    )
    _run_multires_deconv_benchmark(
        subtypes_df=subtypes_df,
        bench_3c_dir=bench_3c_dir,
        sc9_deconv_dir=sc9_deconv_dir,
        bagaev_dir=bagaev_dir,
    )

    sc9_frac_file = sc9_deconv_dir / "deconv_fractions_all_9_combined_res0.5.parquet"
    if sc9_frac_file.is_file():
        sc9_join, sc9_state_cols = _load_and_join_frac_table(sc9_frac_file, subtypes_df)
        if sc9_join.height >= 20 and sc9_state_cols:
            sc9_enr = analyze_subtype_deconvolution_enrichment(sc9_join, sc9_state_cols)
            if isinstance(sc9_enr, Success):
                sc9_enr.unwrap().write_parquet(bagaev_dir / "iatlas_sc9_res05_enrichment_by_bagaev_subtype.parquet")

    print(f"Step 02 complete. Saved zero-leakage deconvolution integration & clustering tables to {bagaev_dir}")


def main() -> None:
    paths = get_data_paths()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bagaev-dir", type=Path, default=paths.output_dir / "bagaev_mfp_analysis")
    parser.add_argument(
        "--bench-3c-dir",
        type=Path,
        default=paths.output_dir / "manuscript_pipeline" / "3c_iatlas_instaprism_milopy_benchmark",
    )
    parser.add_argument(
        "--sc9-deconv-dir",
        type=Path,
        default=paths.output_dir / "sc9_milo_deconv_validation" / "02_deconvolution",
    )
    args = parser.parse_args()
    run_step02(bagaev_dir=args.bagaev_dir, bench_3c_dir=args.bench_3c_dir, sc9_deconv_dir=args.sc9_deconv_dir)


if __name__ == "__main__":
    main()
