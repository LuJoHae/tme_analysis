#!/usr/bin/env python3
"""Step 01: Bagaev et al. (2021) 29-Fges Subtyping & De Novo Clustering Across All 9 iAtlas Cohorts."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Final
import numpy as np
import pandas as pd
import polars as pl
from returns.maybe import Some
from returns.result import Success

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.bagaev_mfp_analysis.shared import (
    assert_inside_data_root,
    eval_tagged_response,
    evaluate_cox_contrasts,
    unwrap_or_raise,
    write_table,
)
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


def _obs_num(obs: object, col: str, n: int) -> list[float]:
    if not isinstance(obs, pd.DataFrame) or col not in obs.columns:
        return [float("nan")] * n
    return [float(v) if pd.notna(v) and not np.isnan(float(v)) else float("nan") for v in obs[col].tolist()]


def _obs_str(obs: object, col: str, n: int, default: str) -> list[str]:
    if not isinstance(obs, pd.DataFrame) or col not in obs.columns:
        return [default] * n
    return [str(v).strip() if pd.notna(v) and str(v).strip() not in ("", "nan", "None") else default for v in obs[col].tolist()]


def _score_cohort(cid: str, collection: object) -> pl.DataFrame:
    adata = unwrap_or_raise(load_dataset(cid, auto_download=False), f"Failed to load {cid}")
    ssg_df = unwrap_or_raise(score_geneset_ssgsea(adata, collection), f"ssGSEA failed on {cid}")  # type: ignore[arg-type]
    obs, n = adata.obs, adata.n_obs
    raw_sids = [str(x).strip() for x in adata.obs_names]
    raw_pids = _obs_str(obs, "PATIENT_ID", n, "")
    pids = [p if p else s for p, s in zip(raw_pids, raw_sids, strict=True)]
    trt = _obs_str(obs, "SAMPLE_TREATMENT", n, "Pre")
    tp = ["Pre" if t.lower() in ("pre", "baseline", "nan", "none", "") else "On-Treatment" for t in trt]
    print(f"Scored 29 Bagaev Fges on {cid:18s}: {n:3d} biopsies ({len(set(pids)):3d} patients)")
    return pl.DataFrame(
        {
            "sample_id": [f"{cid}__{s}" for s in raw_sids],
            "raw_sample_id": raw_sids,
            "patient_id": pids,
            "cohort": [cid] * n,
            "cancer_type": [COHORT_INDICATION_MAP[cid]] * n,
            "sample_treatment": trt,
            "biopsy_timepoint": tp,
            "is_baseline_pre": [t == "Pre" for t in tp],
            "response_binary": _obs_num(obs, "response_binary", n),
            "response_recist": _obs_str(obs, "response_recist", n, "Unknown"),
            "os_months": _obs_num(obs, "os_months", n),
            "os_event": _obs_num(obs, "os_event", n),
            "pfs_months": _obs_num(obs, "pfs_months", n),
            "pfs_event": _obs_num(obs, "pfs_event", n),
        }
    ).hstack(ssg_df.drop("sample_id"))


def run_step01(out_dir: Path) -> None:
    """Execute Step 01 iAtlas Bagaev subtyping, clustering, validation, and clinical evaluation."""
    paths = get_data_paths()
    assert_inside_data_root(paths.raw_dir, paths.data_root)
    assert_inside_data_root(paths.preprocessed_dir, paths.data_root)
    out_dir.mkdir(parents=True, exist_ok=True)
    collection = get_bagaev_mfp_collection()
    tcga_ref = unwrap_or_raise(load_tcga_bagaev_reference())
    print(f"Loaded Bagaev TCGA reference: {tcga_ref.height} samples across {tcga_ref['TCGA_project'].n_unique()} projects.")

    all_raw_df = pl.concat([_score_cohort(cid, collection) for cid in IATLAS_COHORTS], how="vertical")
    write_table(all_raw_df, out_dir, "iatlas_29fges_raw_ssgsea", write_tsv=False)

    primary_baseline_ids = set(
        deduplicate_patient_profiles(
            all_raw_df, mode="strict_pre", patient_col="patient_id", cohort_col="cohort", timepoint_col="biopsy_timepoint"
        )["sample_id"].to_list()
    )
    scaled_df = unwrap_or_raise(
        median_mad_scale_signatures(all_raw_df, BAGAEV_MFP_SIGNATURE_ORDER, group_col=Some("cohort"), clip=Some(4.0)),
        "Median/MAD cohort scaling failed",
    ).with_columns(pl.col("sample_id").is_in(list(primary_baseline_ids)).alias("is_primary_baseline_patient"))
    scaled_ind_df = unwrap_or_raise(
        median_mad_scale_signatures(all_raw_df, BAGAEV_MFP_SIGNATURE_ORDER, group_col=Some("cancer_type"), clip=Some(4.0)),
        "Median/MAD indication scaling failed",
    )

    knn_cfg = BagaevKnnConfig(k_neighbors=35, metric="correlation")
    ind_preds_parts: list[pl.DataFrame] = []
    ind_pool_knn_parts: list[pl.DataFrame] = []
    for cid in IATLAS_COHORTS:
        sub_ref = tcga_ref.filter(pl.col("TCGA_project").is_in(list(IATLAS_TO_TCGA_PROJECTS[cid])))
        ind_preds_parts.append(
            unwrap_or_raise(classify_bagaev_knn(scaled_df.filter(pl.col("cohort") == cid), sub_ref, config=knn_cfg), f"Indication k-NN failed for {cid}")
        )
        r_ind = classify_bagaev_knn(scaled_ind_df.filter(pl.col("cohort") == cid), sub_ref, config=knn_cfg)
        if isinstance(r_ind, Success):
            ind_pool_knn_parts.append(
                r_ind.unwrap().select(["sample_id", pl.col("bagaev_subtype").alias("subtype_knn_indication_pooled_scale")])
            )

    pan_knn_df = unwrap_or_raise(classify_bagaev_knn(scaled_df, tcga_ref, config=knn_cfg), "Pan-Cancer k-NN failed").select(
        ["sample_id", pl.col("bagaev_subtype").alias("subtype_knn_pancancer_tcga")]
    )
    louv_out = unwrap_or_raise(cluster_bagaev_louvain(scaled_df, config=BagaevLouvainConfig(seed=42)), "Louvain clustering failed")
    write_table(louv_out.threshold_sweep_df, out_dir, "iatlas_louvain_threshold_sweep", write_tsv=False)

    per_ind_louv_parts = [
        res.unwrap().assignments_df.select(["sample_id", pl.col("louvain_subtype").alias("subtype_louvain_within_indication")])
        for ind_name in sorted(set(COHORT_INDICATION_MAP.values()))
        if (sub_ind := scaled_df.filter(pl.col("cancer_type") == ind_name)).height >= 12
        for res in (cluster_bagaev_louvain(sub_ind, config=BagaevLouvainConfig(seed=42)),)
        if isinstance(res, Success)
    ]

    pc1, pc2, umap1, umap2 = compute_2d_feature_embedding(scaled_df, BAGAEV_MFP_SIGNATURE_ORDER, metric="correlation", seed=42)
    master_df = (
        scaled_df.join(pl.concat(ind_preds_parts, how="vertical"), on="sample_id", how="left")
        .join(pan_knn_df, on="sample_id", how="left")
        .join(louv_out.assignments_df.select(["sample_id", "louvain_cluster_id", "louvain_subtype"]), on="sample_id", how="left")
        .with_columns([pl.Series("pca_1", pc1), pl.Series("pca_2", pc2), pl.Series("umap_1", umap1), pl.Series("umap_2", umap2)])
    )
    for opt_parts in (ind_pool_knn_parts, per_ind_louv_parts):
        if opt_parts:
            master_df = master_df.join(pl.concat(opt_parts), on="sample_id", how="left")

    means_by_st = {
        str(r["bagaev_subtype"]): r
        for r in master_df.group_by("bagaev_subtype").agg([pl.col(s).mean() for s in BAGAEV_MFP_SIGNATURE_ORDER]).iter_rows(named=True)
    }
    write_table(
        pl.DataFrame(
            [
                {
                    "signature": sig,
                    "sig_order": s_idx,
                    "subtype": st_name,
                    "mean_z": float(means_by_st.get(st_name, {}).get(sig) or 0.0),
                    "n_total_samples": int(master_df.height),
                }
                for s_idx, sig in enumerate(BAGAEV_MFP_SIGNATURE_ORDER)
                for st_name in CANONICAL_SUBTYPE_ORDER
            ]
        ),
        out_dir,
        "iatlas_subtype_29fges_means",
        write_tsv=False,
    )

    pan_mel_res = load_panmelanoma_bagaev_labels()
    master_df = (
        match_iatlas_to_panmelanoma_labels(master_df, pan_mel_res.unwrap())
        if isinstance(pan_mel_res, Success)
        else master_df.with_columns(pl.lit("").alias("published_bagaev_mfp"))
    )
    write_table(master_df, out_dir, "iatlas_bagaev_sample_subtypes", write_tsv=False)
    master_df.select(
        [
            "sample_id", "patient_id", "cohort", "cancer_type", "is_baseline_pre", "is_primary_baseline_patient",
            "response_binary", "response_recist", "bagaev_subtype", "published_bagaev_mfp", "subtype_knn_pancancer_tcga",
            "louvain_subtype", "prob_IE", "prob_IEF", "prob_F", "prob_D", "antitumor_immune_score", "stromal_fibrotic_score",
        ]
    ).write_csv(out_dir / "iatlas_bagaev_sample_subtypes.tsv", separator="\t")

    conc_pairs: list[tuple[str, pl.DataFrame, str, str]] = [
        ("Indication_TCGA_kNN_vs_PanCancer_TCGA_kNN", master_df, "bagaev_subtype", "subtype_knn_pancancer_tcga"),
        ("Indication_TCGA_kNN_vs_DeNovo_Louvain", master_df, "bagaev_subtype", "louvain_subtype"),
    ]
    if "subtype_knn_indication_pooled_scale" in master_df.columns:
        conc_pairs.append(("Cohort_MAD_Scale_vs_Indication_Pooled_Scale", master_df, "bagaev_subtype", "subtype_knn_indication_pooled_scale"))
    matched_pub_df = master_df.filter(pl.col("published_bagaev_mfp") != "")
    if not matched_pub_df.is_empty():
        for iat_cid in ("Hugo-iAtlas", "Riaz-iAtlas", "Liu-iAtlas", "Gide-iAtlas"):
            sub_m = matched_pub_df.filter(pl.col("cohort") == iat_cid)
            if not sub_m.is_empty():
                conc_pairs.append((f"Published_Bagaev_vs_Our_kNN_{iat_cid}", sub_m, "published_bagaev_mfp", "bagaev_subtype"))
        conc_pairs.append(("Published_Bagaev_vs_Our_kNN_All_Melanoma", matched_pub_df, "published_bagaev_mfp", "bagaev_subtype"))

    write_table(
        pl.DataFrame(
            [{"comparison": name, **compare_partition_concordance(df_c[c1].to_list(), df_c[c2].to_list())} for name, df_c, c1, c2 in conc_pairs]
        ),
        out_dir,
        "iatlas_classification_concordance",
    )

    cons_iatlas = evaluate_fges_correlation_conservation(scaled_df, tcga_ref)
    if isinstance(cons_iatlas, Success):
        write_table(pl.DataFrame([{"dataset": "Pan_iAtlas_All9", **cons_iatlas.unwrap()}]), out_dir, "iatlas_vs_tcga_29fges_conservation", write_tsv=False)

    scopes = (
        ("Baseline_Pre", master_df.filter(pl.col("is_primary_baseline_patient"))),
        ("Baseline_AllPreBiopsies", master_df.filter(pl.col("is_baseline_pre"))),
        ("All_Samples", master_df),
    )
    methods = (
        ("bagaev_subtype", "Supervised_TCGA_kNN_Matched"),
        ("subtype_knn_pancancer_tcga", "Supervised_TCGA_kNN_PanCancer"),
        ("louvain_subtype", "Unsupervised_DeNovo_Louvain"),
    )
    resp_tables = [
        tbl
        for scope_name, df_sub in scopes
        for m_col, m_label in methods
        for cid, ctype, sub, strat in [
            ("PanCancer-iAtlas", "Pan-Cancer", df_sub, True),
            *((c, COHORT_INDICATION_MAP[c], df_sub.filter(pl.col("cohort") == c), False) for c in IATLAS_COHORTS),
        ]
        for tbl in eval_tagged_response(
            sub,
            subtype_col=m_col,
            cohort_name=cid,
            sample_scope=scope_name,
            method=m_label,
            stratify_by_cohort=strat,
            extra_exprs=[pl.lit(ctype).alias("cancer_type")],
        )
    ]
    write_table(pl.concat(resp_tables, how="vertical"), out_dir, "iatlas_subtype_response_rates")

    recist_df = (
        master_df.filter(pl.col("is_primary_baseline_patient") & (pl.col("response_recist") != "Unknown"))
        .group_by(["cohort", "bagaev_subtype", "response_recist"])
        .len()
        .rename({"len": "n_samples"})
        .sort(["cohort", "bagaev_subtype", "response_recist"])
    )
    write_table(recist_df, out_dir, "iatlas_subtype_recist_breakdown", write_tsv=False)

    surv_rows = evaluate_cox_contrasts(
        master_df.filter(pl.col("is_primary_baseline_patient")),
        IATLAS_COHORTS,
        (
            ("Favorable_IE_IEF_vs_FD", pl.col("bagaev_subtype").is_in(["IE", "IE/F"])),
            ("IE_vs_Rest", pl.col("bagaev_subtype") == "IE"),
        ),
    )
    if surv_rows:
        write_table(pl.DataFrame(surv_rows), out_dir, "iatlas_subtype_survival_cox")

    print(f"Step 01 complete ({master_df.height} samples, {len(primary_baseline_ids)} unique baseline patients). Saved outputs to {out_dir}")


def main() -> None:
    paths = get_data_paths()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=paths.output_dir / "bagaev_mfp_analysis")
    args = parser.parse_args()
    run_step01(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
