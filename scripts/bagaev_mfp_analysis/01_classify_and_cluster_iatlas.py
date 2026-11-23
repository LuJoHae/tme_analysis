#!/usr/bin/env python3
"""Step 01: Bagaev et al. (2021) 29-Fges Subtyping & De Novo Clustering Across All 9 iAtlas Cohorts.

Computes:
1. Exact 29-Fges ssGSEA scores across all 9 preprocessed iAtlas bulk RNA-seq cohorts.
2. Cohort-wise and indication-pooled median/MAD signature normalization.
3. Supervised K=35 Nearest-Neighbor classification against labeled TCGA reference portraits
   (both indication-matched TCGA projects and Pan-Cancer TCGA, N=7,979).
4. Unsupervised de novo Louvain community detection on thresholded sample Pearson correlation
   graphs (sweeping tau in [0.30, 0.60]) + 2D PCA/UMAP embedding coordinates.
5. Patient-level validation against Bagaev et al.'s published Pan-Melanoma ground-truth MFP
   calls (Hugo, Riaz, Liu, Gide).
6. Objective response rates (ORR), RECIST 1.1 distributions, Fisher/chi2/CMH odds ratios, and
   univariable OS/PFS Cox hazard ratios across subtypes (Pre-treatment baseline & All samples).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Final
import numpy as np
import polars as pl
from returns.maybe import Nothing, Some
from returns.result import Failure, Success
from sklearn.decomposition import PCA

from tme_datasets import (
    BAGAEV_MFP_SIGNATURE_ORDER,
    get_bagaev_mfp_collection,
    get_data_paths,
    load_dataset,
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


def _load_tcga_bagaev_reference(lair_dir: Path) -> pl.DataFrame:
    """Load BostonGene/Bagaev labeled TCGA reference signatures and annotations into Polars."""
    tcga_dir = lair_dir / "DatasetBagaev-TCGA"
    annot_path = tcga_dir / "annotations.tsv"
    sig_path = tcga_dir / "signature.tsv"
    if not annot_path.is_file() or not sig_path.is_file():
        raise FileNotFoundError(
            f"Missing Bagaev TCGA reference files in {tcga_dir}. Expected annotations.tsv and signature.tsv."
        )

    # Read annotations TSV (first column unnamed sample ID)
    annot_lines = annot_path.read_text(encoding="utf-8", errors="replace").splitlines()
    header = ["sample_id", *annot_lines[0].split("\t")[1:]]
    annot_rows: list[dict[str, str]] = []
    for line in annot_lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        row = {header[i]: (parts[i].strip() if i < len(parts) else "") for i in range(len(header))}
        annot_rows.append(row)
    annot_df = pl.DataFrame(annot_rows).select(["sample_id", "TCGA_project", "HISTOLOGICAL_SUBTYPE", "MFP"])

    # Read signature TSV (rows = 29 signatures, columns = TCGA sample IDs)
    sig_lines = sig_path.read_text(encoding="utf-8", errors="replace").splitlines()
    sample_cols = [s.strip() for s in sig_lines[0].split("\t")[1:]]
    sig_matrix: dict[str, list[float]] = {}
    for line in sig_lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        sig_name = parts[0].strip()
        vals = [float(x) if x.strip() not in ("", "NA", "nan") else 0.0 for x in parts[1 : 1 + len(sample_cols)]]
        sig_matrix[sig_name] = vals

    sig_df = pl.DataFrame({"sample_id": sample_cols, **sig_matrix})
    joined = annot_df.join(sig_df, on="sample_id", how="inner").filter(
        pl.col("MFP").is_in(list(CANONICAL_SUBTYPE_ORDER))
    )
    return joined


def _load_panmelanoma_published_labels(lair_dir: Path) -> pl.DataFrame:
    """Load Bagaev et al. published Pan-Melanoma clinical trial MFP labels if available."""
    pan_path = lair_dir / "DatasetBagaev-PanMelanoma" / "annotations.tsv"
    if not pan_path.is_file():
        return pl.DataFrame(schema={"Sample": pl.String, "Cohort": pl.String, "MFP": pl.String, "Patient": pl.String})
    return pl.read_csv(pan_path, separator="\t", infer_schema_length=5000, ignore_errors=True).select(
        [c for c in ("Sample", "Cohort", "MFP", "Patient", "Response", "Time") if c in pl.read_csv(pan_path, separator="\t", n_rows=2).columns]
    )


def _compute_2d_embedding(scaled_df: pl.DataFrame) -> tuple[list[float], list[float], list[float], list[float]]:
    """Compute deterministic 2D PCA and UMAP coordinates on the 29 scaled Fges matrix."""
    mat = scaled_df.select(list(BAGAEV_MFP_SIGNATURE_ORDER)).to_numpy().astype(np.float64)
    pca = PCA(n_components=min(10, mat.shape[1], mat.shape[0]), random_state=42)
    pcs = pca.fit_transform(mat)
    pc1 = pcs[:, 0].tolist()
    pc2 = (pcs[:, 1] if pcs.shape[1] > 1 else np.zeros(mat.shape[0])).tolist()

    try:
        import anndata as ad
        import scanpy as sc

        tmp_ad = ad.AnnData(X=mat.astype(np.float32))
        sc.pp.neighbors(tmp_ad, n_neighbors=min(20, max(3, mat.shape[0] - 1)), use_rep="X", metric="correlation", random_state=42)
        sc.tl.umap(tmp_ad, random_state=42)
        umap_coords = np.asarray(tmp_ad.obsm["X_umap"], dtype=np.float64)
        return pc1, pc2, umap_coords[:, 0].tolist(), umap_coords[:, 1].tolist()
    except Exception:
        return pc1, pc2, pc1, pc2


def _normalize_sample_key(raw_id: str) -> str:
    """Normalize sample/patient identifiers to facilitate matching against published PanMelanoma IDs."""
    s = raw_id.strip()
    if "__" in s:
        s = s.split("__", 1)[1]
    return s.upper().replace("-", "_").replace(".", "_")


def run_step01(out_dir: Path, lair_dir: Path) -> None:
    """Execute Step 01 iAtlas Bagaev subtyping, clustering, validation, and clinical evaluation."""
    out_dir.mkdir(parents=True, exist_ok=True)
    collection = get_bagaev_mfp_collection()
    tcga_ref = _load_tcga_bagaev_reference(lair_dir)
    print(f"Loaded Bagaev TCGA reference: {tcga_ref.height} samples across {tcga_ref['TCGA_project'].n_unique()} projects.")

    # Cache TCGA reference Parquet for downstream steps (Step 03 & Step 04)
    tcga_ref.write_parquet(out_dir / "tcga_bagaev_reference_29fges.parquet")

    raw_cohort_dfs: list[pl.DataFrame] = []
    for cid in IATLAS_COHORTS:
        res = load_dataset(cid, auto_download=False)
        match res:
            case Failure(err):
                raise RuntimeError(f"Failed to load {cid}: {err}")
            case Success(adata):
                ssg_res = score_geneset_ssgsea(adata, collection)
                if isinstance(ssg_res, Failure):
                    raise RuntimeError(f"ssGSEA failed on {cid}: {ssg_res.failure()}")
                ssg_df = ssg_res.unwrap()

                obs = adata.obs
                raw_sids = [str(x) for x in adata.obs_names]
                prefixed_sids = [f"{cid}__{s}" for s in raw_sids]

                resp_bin = (
                    [float(v) if v is not None and not np.isnan(float(v)) else float("nan") for v in obs["response_binary"].tolist()]
                    if "response_binary" in obs.columns
                    else [float("nan")] * adata.n_obs
                )
                resp_rec = (
                    [str(v) for v in obs["response_recist"].tolist()]
                    if "response_recist" in obs.columns
                    else ["Unknown"] * adata.n_obs
                )
                sample_trt = (
                    [str(v) for v in obs["SAMPLE_TREATMENT"].tolist()]
                    if "SAMPLE_TREATMENT" in obs.columns
                    else ["Pre"] * adata.n_obs
                )
                os_m = (
                    [float(v) if v is not None and not np.isnan(float(v)) else float("nan") for v in obs["os_months"].tolist()]
                    if "os_months" in obs.columns
                    else [float("nan")] * adata.n_obs
                )
                os_e = (
                    [float(v) if v is not None and not np.isnan(float(v)) else float("nan") for v in obs["os_event"].tolist()]
                    if "os_event" in obs.columns
                    else [float("nan")] * adata.n_obs
                )
                pfs_m = (
                    [float(v) if v is not None and not np.isnan(float(v)) else float("nan") for v in obs["pfs_months"].tolist()]
                    if "pfs_months" in obs.columns
                    else [float("nan")] * adata.n_obs
                )
                pfs_e = (
                    [float(v) if v is not None and not np.isnan(float(v)) else float("nan") for v in obs["pfs_event"].tolist()]
                    if "pfs_event" in obs.columns
                    else [float("nan")] * adata.n_obs
                )

                meta_df = pl.DataFrame(
                    {
                        "sample_id": prefixed_sids,
                        "raw_sample_id": raw_sids,
                        "cohort": [cid] * adata.n_obs,
                        "cancer_type": [COHORT_INDICATION_MAP[cid]] * adata.n_obs,
                        "sample_treatment": sample_trt,
                        "is_baseline_pre": [
                            str(t).strip().lower() in ("pre", "baseline", "nan", "none", "")
                            for t in sample_trt
                        ],
                        "response_binary": resp_bin,
                        "response_recist": resp_rec,
                        "os_months": os_m,
                        "os_event": os_e,
                        "pfs_months": pfs_m,
                        "pfs_event": pfs_e,
                    }
                )
                merged = meta_df.hstack(ssg_df.drop("sample_id"))
                raw_cohort_dfs.append(merged)
                print(f"Scored 29 Bagaev Fges on {cid:18s}: {adata.n_obs:3d} samples")

    all_raw_df = pl.concat(raw_cohort_dfs, how="vertical")
    all_raw_df.write_parquet(out_dir / "iatlas_29fges_raw_ssgsea.parquet")

    # 1. Within-cohort median/MAD scaling (primary) & indication-pooled scaling (sensitivity)
    scaled_cohort_res = median_mad_scale_signatures(
        all_raw_df, BAGAEV_MFP_SIGNATURE_ORDER, group_col=Some("cohort"), clip=Some(4.0)
    )
    scaled_ind_res = median_mad_scale_signatures(
        all_raw_df, BAGAEV_MFP_SIGNATURE_ORDER, group_col=Some("cancer_type"), clip=Some(4.0)
    )
    if isinstance(scaled_cohort_res, Failure) or isinstance(scaled_ind_res, Failure):
        raise RuntimeError("Median/MAD scaling failed on iAtlas.")
    scaled_df = scaled_cohort_res.unwrap()
    scaled_ind_df = scaled_ind_res.unwrap()

    # 2. Supervised K=35 TCGA k-NN: (a) Indication-matched TCGA & (b) Pan-Cancer TCGA
    knn_cfg = BagaevKnnConfig(k_neighbors=35, metric="correlation")
    ind_preds_parts: list[pl.DataFrame] = []
    for cid in IATLAS_COHORTS:
        sub_q = scaled_df.filter(pl.col("cohort") == cid)
        tcga_projs = list(IATLAS_TO_TCGA_PROJECTS[cid])
        sub_ref = tcga_ref.filter(pl.col("TCGA_project").is_in(tcga_projs))
        pred_res = classify_bagaev_knn(sub_q, sub_ref, config=knn_cfg)
        if isinstance(pred_res, Failure):
            raise RuntimeError(f"Indication k-NN failed for {cid}: {pred_res.failure()}")
        ind_preds_parts.append(pred_res.unwrap())
    ind_knn_df = pl.concat(ind_preds_parts, how="vertical")

    pan_knn_res = classify_bagaev_knn(scaled_df, tcga_ref, config=knn_cfg)
    ind_pool_knn_parts: list[pl.DataFrame] = []
    for cid in IATLAS_COHORTS:
        sub_q_ind = scaled_ind_df.filter(pl.col("cohort") == cid)
        sub_ref = tcga_ref.filter(pl.col("TCGA_project").is_in(list(IATLAS_TO_TCGA_PROJECTS[cid])))
        r_ind = classify_bagaev_knn(sub_q_ind, sub_ref, config=knn_cfg)
        if isinstance(r_ind, Success):
            ind_pool_knn_parts.append(
                r_ind.unwrap().select(
                    [pl.col("sample_id"), pl.col("bagaev_subtype").alias("subtype_knn_indication_pooled_scale")]
                )
            )
    if isinstance(pan_knn_res, Failure):
        raise RuntimeError(f"Pan-Cancer k-NN failed: {pan_knn_res.failure()}")
    pan_knn_df = pan_knn_res.unwrap().select(
        [pl.col("sample_id"), pl.col("bagaev_subtype").alias("subtype_knn_pancancer_tcga")]
    )

    # 3. Unsupervised de novo Louvain clustering (Pan-iAtlas & Per-Indication)
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

    pc1, pc2, umap1, umap2 = _compute_2d_embedding(scaled_df)

    # Assemble master iAtlas sample-level table
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

    master_df.write_parquet(out_dir / "iatlas_bagaev_sample_subtypes.parquet")
    master_df.select(
        [
            "sample_id", "cohort", "cancer_type", "is_baseline_pre",
            "response_binary", "response_recist", "bagaev_subtype",
            "subtype_knn_pancancer_tcga", "louvain_subtype",
            "prob_IE", "prob_IEF", "prob_F", "prob_D",
            "antitumor_immune_score", "stromal_fibrotic_score",
        ]
    ).write_csv(out_dir / "iatlas_bagaev_sample_subtypes.tsv", separator="\t")

    # 4. Validate against published Pan-Melanoma Bagaev MFP labels & inter-method concordance
    concordance_rows: list[dict[str, str | float]] = []
    c_ind_pan = compare_partition_concordance(
        master_df["bagaev_subtype"].to_list(),
        master_df["subtype_knn_pancancer_tcga"].to_list(),
    )
    concordance_rows.append({"comparison": "Indication_TCGA_kNN_vs_PanCancer_TCGA_kNN", **c_ind_pan})

    c_ind_louv = compare_partition_concordance(
        master_df["bagaev_subtype"].to_list(),
        master_df["louvain_subtype"].to_list(),
    )
    concordance_rows.append({"comparison": "Indication_TCGA_kNN_vs_DeNovo_Louvain", **c_ind_louv})

    if "subtype_knn_indication_pooled_scale" in master_df.columns:
        c_scale = compare_partition_concordance(
            master_df["bagaev_subtype"].to_list(),
            master_df["subtype_knn_indication_pooled_scale"].to_list(),
        )
        concordance_rows.append({"comparison": "Cohort_MAD_Scale_vs_Indication_Pooled_Scale", **c_scale})

    pan_mel = _load_panmelanoma_published_labels(lair_dir)
    if not pan_mel.is_empty():
        mel_map = {
            "Hugo-iAtlas": "Hugo",
            "Riaz-iAtlas": "Riaz",
            "Liu-iAtlas": "Liu",
            "Gide-iAtlas": "Gide",
        }
        matched_pub: list[str] = []
        matched_our: list[str] = []
        for iat_cid, bg_cid in mel_map.items():
            sub_pub = pan_mel.filter(
                (pl.col("Cohort") == bg_cid) & pl.col("MFP").is_in(list(CANONICAL_SUBTYPE_ORDER))
            )
            pub_lookup: dict[str, str] = {}
            for r in sub_pub.iter_rows(named=True):
                mfp_lbl = str(r["MFP"])
                samp_str = str(r.get("Sample") or "").strip()
                pat_str = str(r.get("Patient") or "").strip()
                tm_str = str(r.get("Time") or "").strip()
                pub_lookup[_normalize_sample_key(samp_str)] = mfp_lbl
                if pat_str:
                    pub_lookup[_normalize_sample_key(pat_str)] = mfp_lbl
                if bg_cid == "Riaz" and pat_str.startswith("Riaz_") and tm_str:
                    # Riaz_Pt1 + Pre -> PT1_PRE
                    pt_core = pat_str.removeprefix("Riaz_")
                    pub_lookup[_normalize_sample_key(f"{pt_core}_{tm_str}")] = mfp_lbl
                elif bg_cid == "Liu" and samp_str.startswith("Patient"):
                    # Patient1 -> LIU_SAMPLE1
                    num_core = samp_str.removeprefix("Patient")
                    pub_lookup[_normalize_sample_key(f"Liu_Sample{num_core}")] = mfp_lbl
                elif bg_cid == "Gide" and pat_str.startswith("Gide_") and tm_str:
                    # Gide_PD1_2 + Pre -> PD02_PRE ; Gide_ipiPD1_1 + Pre -> IP01_PRE
                    if pat_str.startswith("Gide_PD1_"):
                        num_s = pat_str.removeprefix("Gide_PD1_")
                        if num_s.isdigit():
                            pub_lookup[_normalize_sample_key(f"PD{int(num_s):02d}_{tm_str}")] = mfp_lbl
                    elif pat_str.startswith("Gide_ipiPD1_"):
                        num_s = pat_str.removeprefix("Gide_ipiPD1_")
                        if num_s.isdigit():
                            pub_lookup[_normalize_sample_key(f"iP{int(num_s):02d}_{tm_str}")] = mfp_lbl

            sub_iat = master_df.filter(pl.col("cohort") == iat_cid)
            c_pub: list[str] = []
            c_our: list[str] = []
            for r in sub_iat.iter_rows(named=True):
                raw_k = _normalize_sample_key(str(r["raw_sample_id"]))
                base_k = raw_k.split("_PRE")[0].split("_EDT")[0].split("_ON")[0]
                hit = pub_lookup.get(raw_k) or pub_lookup.get(base_k)
                if hit is not None:
                    c_pub.append(hit)
                    c_our.append(str(r["bagaev_subtype"]))
            if c_pub:
                m_c = compare_partition_concordance(c_pub, c_our)
                concordance_rows.append({"comparison": f"Published_Bagaev_vs_Our_kNN_{iat_cid}", **m_c})
                matched_pub.extend(c_pub)
                matched_our.extend(c_our)
        if matched_pub:
            m_all = compare_partition_concordance(matched_pub, matched_our)
            concordance_rows.append({"comparison": "Published_Bagaev_vs_Our_kNN_All_Melanoma", **m_all})

    # Also evaluate 29x29 correlation matrix conservation of iAtlas vs TCGA
    cons_iatlas = evaluate_fges_correlation_conservation(scaled_df, tcga_ref)
    if isinstance(cons_iatlas, Success):
        pl.DataFrame([{"dataset": "Pan_iAtlas_All9", **cons_iatlas.unwrap()}]).write_parquet(
            out_dir / "iatlas_vs_tcga_29fges_conservation.parquet"
        )

    concordance_df = pl.DataFrame(concordance_rows)
    concordance_df.write_parquet(out_dir / "iatlas_classification_concordance.parquet")
    concordance_df.write_csv(out_dir / "iatlas_classification_concordance.tsv", separator="\t")

    # 5. Response rate associations across Pan-Cancer & Per-Cohort (Pre-treatment & All samples)
    resp_tables: list[pl.DataFrame] = []
    for sample_filter_name, df_sub in (
        ("Baseline_Pre", master_df.filter(pl.col("is_baseline_pre"))),
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

    # RECIST 1.1 breakdown across subtypes
    recist_df = (
        master_df.filter(pl.col("is_baseline_pre") & (pl.col("response_recist") != "Unknown"))
        .group_by(["cohort", "bagaev_subtype", "response_recist"])
        .len()
        .rename({"len": "n_samples"})
        .sort(["cohort", "bagaev_subtype", "response_recist"])
    )
    recist_df.write_parquet(out_dir / "iatlas_subtype_recist_breakdown.parquet")

    # 6. Survival (OS & PFS) Cox HR for Favorable (IE + IE/F) vs Unfavorable (F + D) and IE vs rest
    surv_rows: list[dict[str, str | float | int]] = []
    pre_df = master_df.filter(pl.col("is_baseline_pre"))
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

            for contrast_name, indicator in (
                ("Favorable_IE_IEF_vs_FD", valid_s["bagaev_subtype"].is_in(["IE", "IE/F"]).to_numpy().astype(np.float64)),
                ("IE_vs_Rest", (valid_s["bagaev_subtype"] == "IE").to_numpy().astype(np.float64)),
            ):
                if np.std(indicator) <= 1e-8:
                    continue
                cox_res = fit_univariable_cox(times, events, indicator)
                if isinstance(cox_res, Success):
                    cx = cox_res.unwrap()
                    surv_rows.append(
                        {
                            "cohort": scope_id,
                            "endpoint": endpoint,
                            "contrast": contrast_name,
                            "n_samples": int(valid_s.height),
                            "n_events": int(np.sum(events)),
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

    print(f"Step 01 complete. Saved outputs to {out_dir}")


def main() -> None:
    paths = get_data_paths()
    default_lair = Path("/storage/halu/lair") if Path("/storage/halu/lair").exists() else paths.lair_dir
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
        default=default_lair,
        help="Directory containing DatasetBagaev-TCGA and DatasetBagaev-PanMelanoma.",
    )
    args = parser.parse_args()
    run_step01(out_dir=args.out_dir, lair_dir=args.lair_dir)


if __name__ == "__main__":
    main()
