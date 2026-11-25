#!/usr/bin/env python3
"""Step 01: Bagaev et al. (2021) 29-Fges Subtyping & De Novo Clustering Across All 9 iAtlas Cohorts.

Orchestrates:
1. Exact 29-Fges ssGSEA scores across all 9 preprocessed iAtlas bulk RNA-seq cohorts.
2. Cohort-wise and indication-pooled median/MAD signature normalization.
3. Supervised K=35 Nearest-Neighbor classification against labeled TCGA reference portraits.
4. Unsupervised de novo 2-phase Louvain clustering + 2D PCA/UMAP embedding.
5. Timepoint-preserving validation against Bagaev et al.'s published Pan-Melanoma MFP labels.
6. Patient-deduplicated ORR, RECIST 1.1 distributions, CMH odds ratios, and cohort-stratified
   unscaled binary OS/PFS Cox hazard ratios.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Final
import numpy as np
import polars as pl
from returns.maybe import Nothing, Some
from returns.result import Failure, Success

from tme_datasets import (
    BAGAEV_MFP_SIGNATURE_ORDER,
    deduplicate_patient_profiles,
    get_bagaev_mfp_collection,
    get_data_paths,
    load_dataset,
    load_panmelanoma_bagaev_labels,
    load_tcga_bagaev_reference,
    match_iatlas_to_panmelanoma_labels,
    median_mad_scale_signatures,
    score_geneset_ssgsea,
)
from tme_datasets.providers.bulk_iatlas import IATLAS_COHORTS
from tme_response import (
    CANONICAL_SUBTYPE_ORDER,
    IATLAS_TO_TCGA_PROJECTS,
    BagaevKnnConfig,
    BagaevLouvainConfig,
    classify_bagaev_knn,
    cluster_bagaev_louvain,
    compare_partition_concordance,
    compute_2d_feature_embedding,
    evaluate_fges_correlation_conservation,
    evaluate_subtype_response_associations,
    fit_univariable_cox,
)

COHORT_INDICATION_MAP: Final[dict[str, str]] = {
    "Hugo-iAtlas": "Melanoma",
    "Riaz-iAtlas": "Melanoma",
    "Liu-iAtlas": "Melanoma",
    "Gide-iAtlas": "Melanoma",
    "Rosenberg-iAtlas": "Urothelial Bladder",
    "Padron-iAtlas": "PDAC",
    "Anders-iAtlas": "Breast",
    "McDermott-iAtlas": "ccRCC",
    "Choueiri-iAtlas": "ccRCC",
}


def _extract_numeric_obs_col(obs: object, col_name: str, n_obs: int) -> list[float]:
    """Safely extract numeric column from AnnData `.obs` as float list with NaN fallback."""
    import pandas as pd

    if not isinstance(obs, pd.DataFrame) or col_name not in obs.columns:
        return [float("nan")] * n_obs
    return [
        float(v) if v is not None and not np.isnan(float(v)) else float("nan")
        for v in obs[col_name].tolist()
    ]


def _extract_str_obs_col(obs: object, col_name: str, n_obs: int, default_val: str) -> list[str]:
    """Safely extract string column from AnnData `.obs` with fallback."""
    import pandas as pd

    if not isinstance(obs, pd.DataFrame) or col_name not in obs.columns:
        return [default_val] * n_obs
    return [
        str(v).strip() if v is not None and str(v).strip() not in ("", "nan", "None") else default_val
        for v in obs[col_name].tolist()
    ]


def run_step01(out_dir: Path, lair_dir: Path) -> None:
    """Execute Step 01 iAtlas Bagaev subtyping, clustering, validation, and clinical evaluation."""
    out_dir.mkdir(parents=True, exist_ok=True)
    collection = get_bagaev_mfp_collection()
    tcga_res = load_tcga_bagaev_reference(lair_dir)
    if isinstance(tcga_res, Failure):
        raise RuntimeError(tcga_res.failure())
    tcga_ref = tcga_res.unwrap()
    print(f"Loaded Bagaev TCGA reference: {tcga_ref.height} samples across {tcga_ref['TCGA_project'].n_unique()} projects.")
    tcga_ref.write_parquet(out_dir / "tcga_bagaev_reference_29fges.parquet")

    raw_cohort_dfs: list[pl.DataFrame] = []
    for cid in IATLAS_COHORTS:
        res = load_dataset(cid, auto_download=False)
        if isinstance(res, Failure):
            raise RuntimeError(f"Failed to load {cid}: {res.failure()}")
        adata = res.unwrap()
        ssg_res = score_geneset_ssgsea(adata, collection)
        if isinstance(ssg_res, Failure):
            raise RuntimeError(f"ssGSEA failed on {cid}: {ssg_res.failure()}")
        ssg_df = ssg_res.unwrap()

        obs = adata.obs
        raw_sids = [str(x).strip() for x in adata.obs_names]
        prefixed_sids = [f"{cid}__{s}" for s in raw_sids]

        patient_ids = (
            _extract_str_obs_col(obs, "PATIENT_ID", adata.n_obs, "")
            if hasattr(obs, "columns") and "PATIENT_ID" in obs.columns
            else raw_sids
        )
        patient_ids = [p if p else sid for p, sid in zip(patient_ids, raw_sids, strict=True)]
        sample_trt = _extract_str_obs_col(obs, "SAMPLE_TREATMENT", adata.n_obs, "Pre")
        biopsy_tp = [
            "Pre" if t.strip().lower() in ("pre", "baseline", "nan", "none", "") else "On-Treatment"
            for t in sample_trt
        ]

        meta_df = pl.DataFrame(
            {
                "sample_id": prefixed_sids,
                "raw_sample_id": raw_sids,
                "patient_id": patient_ids,
                "cohort": [cid] * adata.n_obs,
                "cancer_type": [COHORT_INDICATION_MAP[cid]] * adata.n_obs,
                "sample_treatment": sample_trt,
                "biopsy_timepoint": biopsy_tp,
                "is_baseline_pre": [t == "Pre" for t in biopsy_tp],
                "response_binary": _extract_numeric_obs_col(obs, "response_binary", adata.n_obs),
                "response_recist": _extract_str_obs_col(obs, "response_recist", adata.n_obs, "Unknown"),
                "os_months": _extract_numeric_obs_col(obs, "os_months", adata.n_obs),
                "os_event": _extract_numeric_obs_col(obs, "os_event", adata.n_obs),
                "pfs_months": _extract_numeric_obs_col(obs, "pfs_months", adata.n_obs),
                "pfs_event": _extract_numeric_obs_col(obs, "pfs_event", adata.n_obs),
            }
        )
        merged = meta_df.hstack(ssg_df.drop("sample_id"))
        raw_cohort_dfs.append(merged)
        print(f"Scored 29 Bagaev Fges on {cid:18s}: {adata.n_obs:3d} biopsies ({len(set(patient_ids)):3d} patients)")

    all_raw_df = pl.concat(raw_cohort_dfs, how="vertical")
    all_raw_df.write_parquet(out_dir / "iatlas_29fges_raw_ssgsea.parquet")

    # Identify strictly independent 1-per-patient baseline tumor profiles
    dedup_base_df = deduplicate_patient_profiles(
        all_raw_df,
        mode="strict_pre",
        patient_col="patient_id",
        cohort_col="cohort",
        timepoint_col="biopsy_timepoint",
    )
    primary_baseline_ids = set(dedup_base_df["sample_id"].to_list())

    scaled_cohort_res = median_mad_scale_signatures(
        all_raw_df, BAGAEV_MFP_SIGNATURE_ORDER, group_col=Some("cohort"), clip=Some(4.0)
    )
    scaled_ind_res = median_mad_scale_signatures(
        all_raw_df, BAGAEV_MFP_SIGNATURE_ORDER, group_col=Some("cancer_type"), clip=Some(4.0)
    )
    if isinstance(scaled_cohort_res, Failure) or isinstance(scaled_ind_res, Failure):
        raise RuntimeError("Median/MAD scaling failed on iAtlas.")
    scaled_df = scaled_cohort_res.unwrap().with_columns(
        pl.col("sample_id").is_in(list(primary_baseline_ids)).alias("is_primary_baseline_patient")
    )
    scaled_ind_df = scaled_ind_res.unwrap()

    # Supervised K=35 TCGA k-NN (indication-matched & Pan-Cancer TCGA)
    knn_cfg = BagaevKnnConfig(k_neighbors=35, metric="correlation")
    ind_preds_parts: list[pl.DataFrame] = []
    ind_pool_knn_parts: list[pl.DataFrame] = []
    for cid in IATLAS_COHORTS:
        sub_q = scaled_df.filter(pl.col("cohort") == cid)
        sub_ref = tcga_ref.filter(pl.col("TCGA_project").is_in(list(IATLAS_TO_TCGA_PROJECTS[cid])))
        pred_res = classify_bagaev_knn(sub_q, sub_ref, config=knn_cfg)
        if isinstance(pred_res, Failure):
            raise RuntimeError(f"Indication k-NN failed for {cid}: {pred_res.failure()}")
        ind_preds_parts.append(pred_res.unwrap())

        r_ind = classify_bagaev_knn(scaled_ind_df.filter(pl.col("cohort") == cid), sub_ref, config=knn_cfg)
        if isinstance(r_ind, Success):
            ind_pool_knn_parts.append(
                r_ind.unwrap().select(
                    [pl.col("sample_id"), pl.col("bagaev_subtype").alias("subtype_knn_indication_pooled_scale")]
                )
            )

    ind_knn_df = pl.concat(ind_preds_parts, how="vertical")
    pan_knn_res = classify_bagaev_knn(scaled_df, tcga_ref, config=knn_cfg)
    if isinstance(pan_knn_res, Failure):
        raise RuntimeError(f"Pan-Cancer k-NN failed: {pan_knn_res.failure()}")
    pan_knn_df = pan_knn_res.unwrap().select(
        [pl.col("sample_id"), pl.col("bagaev_subtype").alias("subtype_knn_pancancer_tcga")]
    )

    # Unsupervised 2-phase Louvain clustering + deterministic 2D embedding
    louv_res = cluster_bagaev_louvain(scaled_df, config=BagaevLouvainConfig(seed=42))
    if isinstance(louv_res, Failure):
        raise RuntimeError(f"Louvain clustering failed: {louv_res.failure()}")
    louv_out = louv_res.unwrap()
    louv_out.threshold_sweep_df.write_parquet(out_dir / "iatlas_louvain_threshold_sweep.parquet")

    per_ind_louv_parts: list[pl.DataFrame] = []
    for ind_name in sorted(set(COHORT_INDICATION_MAP.values())):
        sub_ind = scaled_df.filter(pl.col("cancer_type") == ind_name)
        if sub_ind.height >= 12:
            sub_louv = cluster_bagaev_louvain(sub_ind, config=BagaevLouvainConfig(seed=42))
            if isinstance(sub_louv, Success):
                per_ind_louv_parts.append(
                    sub_louv.unwrap().assignments_df.select(
                        [pl.col("sample_id"), pl.col("louvain_subtype").alias("subtype_louvain_within_indication")]
                    )
                )

    pc1, pc2, umap1, umap2 = compute_2d_feature_embedding(
        scaled_df, BAGAEV_MFP_SIGNATURE_ORDER, metric="correlation", seed=42
    )

    master_df = (
        scaled_df.join(ind_knn_df, on="sample_id", how="left")
        .join(pan_knn_df, on="sample_id", how="left")
        .join(
            louv_out.assignments_df.select(["sample_id", "louvain_cluster_id", "louvain_subtype"]),
            on="sample_id",
            how="left",
        )
        .with_columns(
            [
                pl.Series("pca_1", pc1),
                pl.Series("pca_2", pc2),
                pl.Series("umap_1", umap1),
                pl.Series("umap_2", umap2),
            ]
        )
    )
    if ind_pool_knn_parts:
        master_df = master_df.join(pl.concat(ind_pool_knn_parts), on="sample_id", how="left")
    if per_ind_louv_parts:
        master_df = master_df.join(pl.concat(per_ind_louv_parts), on="sample_id", how="left")

    # Pre-compute subtype mean 29-Fges portrait table for Step 04 Fig 1 Panel A
    subtype_mean_rows: list[dict[str, str | int | float]] = []
    for s_idx, sig in enumerate(BAGAEV_MFP_SIGNATURE_ORDER):
        for st_name in CANONICAL_SUBTYPE_ORDER:
            sub = master_df.filter(pl.col("bagaev_subtype") == st_name)
            raw_m = sub[sig].mean() if not sub.is_empty() else 0.0
            mean_z = float(raw_m) if isinstance(raw_m, (int, float)) else 0.0
            subtype_mean_rows.append(
                {
                    "signature": sig,
                    "sig_order": s_idx,
                    "subtype": st_name,
                    "mean_z": mean_z,
                    "n_total_samples": int(master_df.height),
                }
            )
    pl.DataFrame(subtype_mean_rows).write_parquet(out_dir / "iatlas_subtype_29fges_means.parquet")

    # Timepoint-preserving Pan-Melanoma published label validation
    pan_mel_res = load_panmelanoma_bagaev_labels(lair_dir)
    if isinstance(pan_mel_res, Success):
        master_df = match_iatlas_to_panmelanoma_labels(master_df, pan_mel_res.unwrap())
    else:
        master_df = master_df.with_columns(pl.lit("").alias("published_bagaev_mfp"))

    master_df.write_parquet(out_dir / "iatlas_bagaev_sample_subtypes.parquet")
    master_df.select(
        [
            "sample_id", "patient_id", "cohort", "cancer_type",
            "is_baseline_pre", "is_primary_baseline_patient",
            "response_binary", "response_recist", "bagaev_subtype",
            "published_bagaev_mfp", "subtype_knn_pancancer_tcga", "louvain_subtype",
            "prob_IE", "prob_IEF", "prob_F", "prob_D",
            "antitumor_immune_score", "stromal_fibrotic_score",
        ]
    ).write_csv(out_dir / "iatlas_bagaev_sample_subtypes.tsv", separator="\t")

    concordance_rows: list[dict[str, str | float]] = [
        {
            "comparison": "Indication_TCGA_kNN_vs_PanCancer_TCGA_kNN",
            **compare_partition_concordance(
                master_df["bagaev_subtype"].to_list(),
                master_df["subtype_knn_pancancer_tcga"].to_list(),
            ),
        },
        {
            "comparison": "Indication_TCGA_kNN_vs_DeNovo_Louvain",
            **compare_partition_concordance(
                master_df["bagaev_subtype"].to_list(),
                master_df["louvain_subtype"].to_list(),
            ),
        },
    ]
    if "subtype_knn_indication_pooled_scale" in master_df.columns:
        concordance_rows.append(
            {
                "comparison": "Cohort_MAD_Scale_vs_Indication_Pooled_Scale",
                **compare_partition_concordance(
                    master_df["bagaev_subtype"].to_list(),
                    master_df["subtype_knn_indication_pooled_scale"].to_list(),
                ),
            }
        )

    matched_pub_df = master_df.filter(pl.col("published_bagaev_mfp") != "")
    if not matched_pub_df.is_empty():
        for iat_cid in ("Hugo-iAtlas", "Riaz-iAtlas", "Liu-iAtlas", "Gide-iAtlas"):
            sub_m = matched_pub_df.filter(pl.col("cohort") == iat_cid)
            if not sub_m.is_empty():
                concordance_rows.append(
                    {
                        "comparison": f"Published_Bagaev_vs_Our_kNN_{iat_cid}",
                        **compare_partition_concordance(
                            sub_m["published_bagaev_mfp"].to_list(),
                            sub_m["bagaev_subtype"].to_list(),
                        ),
                    }
                )
        concordance_rows.append(
            {
                "comparison": "Published_Bagaev_vs_Our_kNN_All_Melanoma",
                **compare_partition_concordance(
                    matched_pub_df["published_bagaev_mfp"].to_list(),
                    matched_pub_df["bagaev_subtype"].to_list(),
                ),
            }
        )

    cons_iatlas = evaluate_fges_correlation_conservation(scaled_df, tcga_ref)
    if isinstance(cons_iatlas, Success):
        pl.DataFrame([{"dataset": "Pan_iAtlas_All9", **cons_iatlas.unwrap()}]).write_parquet(
            out_dir / "iatlas_vs_tcga_29fges_conservation.parquet"
        )

    concordance_df = pl.DataFrame(concordance_rows)
    concordance_df.write_parquet(out_dir / "iatlas_classification_concordance.parquet")
    concordance_df.write_csv(out_dir / "iatlas_classification_concordance.tsv", separator="\t")

    # Response rate associations: primary Baseline_Pre = deduplicated unique patients
    resp_tables: list[pl.DataFrame] = []
    for sample_filter_name, df_sub in (
        ("Baseline_Pre", master_df.filter(pl.col("is_primary_baseline_patient"))),
        ("Baseline_AllPreBiopsies", master_df.filter(pl.col("is_baseline_pre"))),
        ("All_Samples", master_df),
    ):
        for method_col, method_label in (
            ("bagaev_subtype", "Supervised_TCGA_kNN_Matched"),
            ("subtype_knn_pancancer_tcga", "Supervised_TCGA_kNN_PanCancer"),
            ("louvain_subtype", "Unsupervised_DeNovo_Louvain"),
        ):
            pan_res = evaluate_subtype_response_associations(
                df_sub, subtype_col=method_col, cohort_col=Some("cohort")
            )
            if isinstance(pan_res, Success):
                resp_tables.append(
                    pan_res.unwrap().with_columns(
                        [
                            pl.lit("PanCancer-iAtlas").alias("cohort"),
                            pl.lit("Pan-Cancer").alias("cancer_type"),
                            pl.lit(sample_filter_name).alias("sample_scope"),
                            pl.lit(method_label).alias("method"),
                        ]
                    )
                )
            for cid in IATLAS_COHORTS:
                c_sub = df_sub.filter(pl.col("cohort") == cid)
                c_res = evaluate_subtype_response_associations(c_sub, subtype_col=method_col, cohort_col=Nothing)
                if isinstance(c_res, Success):
                    resp_tables.append(
                        c_res.unwrap().with_columns(
                            [
                                pl.lit(cid).alias("cohort"),
                                pl.lit(COHORT_INDICATION_MAP[cid]).alias("cancer_type"),
                                pl.lit(sample_filter_name).alias("sample_scope"),
                                pl.lit(method_label).alias("method"),
                            ]
                        )
                    )

    all_resp_df = pl.concat(resp_tables, how="vertical")
    all_resp_df.write_parquet(out_dir / "iatlas_subtype_response_rates.parquet")
    all_resp_df.write_csv(out_dir / "iatlas_subtype_response_rates.tsv", separator="\t")

    recist_df = (
        master_df.filter(pl.col("is_primary_baseline_patient") & (pl.col("response_recist") != "Unknown"))
        .group_by(["cohort", "bagaev_subtype", "response_recist"])
        .len()
        .rename({"len": "n_samples"})
        .sort(["cohort", "bagaev_subtype", "response_recist"])
    )
    recist_df.write_parquet(out_dir / "iatlas_subtype_recist_breakdown.parquet")

    # Survival (OS & PFS) unscaled binary Cox PH (cohort-stratified for PanCancer-iAtlas)
    surv_rows: list[dict[str, str | float | int]] = []
    pre_df = master_df.filter(pl.col("is_primary_baseline_patient"))
    for endpoint, t_col, e_col in (("OS", "os_months", "os_event"), ("PFS", "pfs_months", "pfs_event")):
        for scope_id, sub_s in [("PanCancer-iAtlas", pre_df)] + [
            (cid, pre_df.filter(pl.col("cohort") == cid)) for cid in IATLAS_COHORTS
        ]:
            valid_s = sub_s.filter(
                pl.col(t_col).is_not_null()
                & pl.col(t_col).is_finite()
                & (pl.col(t_col) > 0.0)
                & pl.col(e_col).is_in([0.0, 1.0])
            )
            if valid_s.height < 15:
                continue
            times = valid_s[t_col].to_numpy().astype(np.float64)
            events = valid_s[e_col].to_numpy().astype(np.int64)
            if int(np.sum(events)) < 5:
                continue
            strata_arg = Some(valid_s["cohort"].cast(pl.String).to_list()) if scope_id == "PanCancer-iAtlas" else Nothing

            for contrast_name, indicator in (
                ("Favorable_IE_IEF_vs_FD", valid_s["bagaev_subtype"].is_in(["IE", "IE/F"]).to_numpy().astype(np.float64)),
                ("IE_vs_Rest", (valid_s["bagaev_subtype"] == "IE").to_numpy().astype(np.float64)),
            ):
                cox_res = fit_univariable_cox(
                    times.tolist(), events.tolist(), indicator.tolist(), standardize_score=False, strata=strata_arg
                )
                if isinstance(cox_res, Success):
                    cx = cox_res.unwrap()
                    surv_rows.append(
                        {
                            "cohort": scope_id,
                            "endpoint": endpoint,
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
        surv_df.write_parquet(out_dir / "iatlas_subtype_survival_cox.parquet")
        surv_df.write_csv(out_dir / "iatlas_subtype_survival_cox.tsv", separator="\t")

    print(f"Step 01 complete ({master_df.height} samples, {len(primary_baseline_ids)} unique baseline patients). Saved outputs to {out_dir}")


def main() -> None:
    paths = get_data_paths()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=paths.output_dir / "bagaev_mfp_analysis",
        help="Output directory for Step 01 Parquet/TSV artifacts.",
    )
    parser.add_argument(
        "--lair-dir",
        type=Path,
        default=paths.lair_dir,
        help="Directory containing DatasetBagaev-TCGA and DatasetBagaev-PanMelanoma.",
    )
    args = parser.parse_args()
    run_step01(out_dir=args.out_dir, lair_dir=args.lair_dir)


if __name__ == "__main__":
    main()
