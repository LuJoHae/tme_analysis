#!/usr/bin/env python3
"""Step 03: Tumor-Only Bagaev 29-Fges Pseudobulk Clustering & Response Analysis Across 23 SC Cohorts."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Final, Literal
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Nothing, Some
from returns.result import Failure, Success
from scipy import stats
from sklearn.metrics import roc_auc_score  # type: ignore[import-untyped]

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.bagaev_mfp_analysis.shared import assert_inside_data_root, eval_tagged_response, unwrap_or_raise, write_table
from tme_datasets import (
    BAGAEV_MFP_SIGNATURE_ORDER,
    TumorPseudobulkSpec,
    compute_ssgsea_matrix,
    deduplicate_patient_profiles,
    get_bagaev_mfp_collection,
    get_data_paths,
    load_ensembl_symbol_map,
    load_tcga_bagaev_reference,
    median_mad_scale_signatures,
    stream_cohort_tumor_pseudobulk,
)
from tme_response import (
    CANONICAL_SUBTYPE_ORDER,
    INDICATION_TO_TCGA_PROJECTS,
    BagaevKnnConfig,
    BagaevLouvainConfig,
    classify_bagaev_knn,
    cluster_bagaev_louvain,
    evaluate_fges_correlation_conservation,
)

CACHE_VERSION_PREFIX: Final[str] = "v2_tumor_only"


class SingleCellCohortSpec(BaseModel):
    """Immutable declarative specification for each of the 23 single-cell cohorts."""

    model_config = ConfigDict(frozen=True)

    cohort_id: str
    dataset_id: str
    display_name: str
    cohort_group: Literal["SC9_Milopy", "Atlas_14"]
    indication: str
    sorting_regime: Literal["Whole_Tumor", "FACS_Immune_Sorted"]
    sample_col: str
    has_icb_response: bool
    min_cells_per_sample: int = 25


SC23_SPECS: Final[tuple[SingleCellCohortSpec, ...]] = (
    SingleCellCohortSpec(cohort_id="GSE120575", dataset_id="GSE120575", display_name="Sade-Feldman 2018 (Melanoma)", cohort_group="SC9_Milopy", indication="Melanoma", sorting_regime="FACS_Immune_Sorted", sample_col="patient_geo", has_icb_response=True),
    SingleCellCohortSpec(cohort_id="CELLxGENE_7b20c613", dataset_id="CELLxGENE_7b20c613", display_name="Gondal 2025 (Pan-Cancer)", cohort_group="SC9_Milopy", indication="Pan-Cancer", sorting_regime="Whole_Tumor", sample_col="donor_id_pre_post", has_icb_response=True),
    SingleCellCohortSpec(cohort_id="CELLxGENE_05a8c945", dataset_id="CELLxGENE_05a8c945", display_name="Marteau 2026 (CRC)", cohort_group="SC9_Milopy", indication="CRC", sorting_regime="Whole_Tumor", sample_col="sample_id", has_icb_response=True),
    SingleCellCohortSpec(cohort_id="CELLxGENE_6f9de485", dataset_id="CELLxGENE_6f9de485", display_name="Yan 2026 (Breast)", cohort_group="SC9_Milopy", indication="Breast", sorting_regime="Whole_Tumor", sample_col="sample_id", has_icb_response=True),
    SingleCellCohortSpec(cohort_id="GSE207422", dataset_id="GSE207422", display_name="Hu 2023 (NSCLC)", cohort_group="SC9_Milopy", indication="NSCLC", sorting_regime="Whole_Tumor", sample_col="sample_id", has_icb_response=True),
    SingleCellCohortSpec(cohort_id="GSE243013", dataset_id="GSE243013", display_name="Liu 2025 (NSCLC)", cohort_group="SC9_Milopy", indication="NSCLC", sorting_regime="Whole_Tumor", sample_col="sample_id", has_icb_response=True),
    SingleCellCohortSpec(cohort_id="GSE233203", dataset_id="GSE233203", display_name="Hong 2025 (NSCLC)", cohort_group="SC9_Milopy", indication="NSCLC", sorting_regime="Whole_Tumor", sample_col="sample_id", has_icb_response=True, min_cells_per_sample=20),
    SingleCellCohortSpec(cohort_id="GSE200996", dataset_id="GSE200996", display_name="Luoma 2022 (HNSCC)", cohort_group="SC9_Milopy", indication="HNSCC", sorting_regime="FACS_Immune_Sorted", sample_col="sample_id", has_icb_response=True),
    SingleCellCohortSpec(cohort_id="GSE316195", dataset_id="GSE316195", display_name="Raufi 2026 (PDAC)", cohort_group="SC9_Milopy", indication="PDAC", sorting_regime="Whole_Tumor", sample_col="sample_id", has_icb_response=True),
    SingleCellCohortSpec(cohort_id="Azizi_2018_BRCA", dataset_id="GSE114727", display_name="Azizi 2018 (Breast)", cohort_group="Atlas_14", indication="Breast", sorting_regime="FACS_Immune_Sorted", sample_col="sample", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Becker_2022_COAD", dataset_id="GSE201349", display_name="Becker 2022 (COAD)", cohort_group="Atlas_14", indication="CRC", sorting_regime="Whole_Tumor", sample_col="sample", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Biermann_2022_BrainMet", dataset_id="GSE200218", display_name="Biermann 2022 (Mel BrainMet)", cohort_group="Atlas_14", indication="Melanoma", sorting_regime="Whole_Tumor", sample_col="orig.ident", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Borcherding_2021_ccRCC", dataset_id="GSE121638", display_name="Borcherding 2021 (ccRCC)", cohort_group="Atlas_14", indication="ccRCC", sorting_regime="FACS_Immune_Sorted", sample_col="sample", has_icb_response=False, min_cells_per_sample=20),
    SingleCellCohortSpec(cohort_id="Cheng_2021_PanMyeloid", dataset_id="GSE154763", display_name="Cheng 2021 (Pan-Myeloid)", cohort_group="Atlas_14", indication="Pan-Cancer", sorting_regime="FACS_Immune_Sorted", sample_col="patient", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Durante_2020_UVM", dataset_id="GSE139829", display_name="Durante 2020 (UVM)", cohort_group="Atlas_14", indication="UVM", sorting_regime="Whole_Tumor", sample_col="sample", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Khaliq_2022_CRC", dataset_id="GSE200997", display_name="Khaliq 2022 (CRC)", cohort_group="Atlas_14", indication="CRC", sorting_regime="Whole_Tumor", sample_col="samples", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Kim_2020_LUAD", dataset_id="GSE131907", display_name="Kim 2020 (LUAD)", cohort_group="Atlas_14", indication="LUAD", sorting_regime="Whole_Tumor", sample_col="Sample", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Leader_2021_NSCLC", dataset_id="GSE154826", display_name="Leader 2021 (NSCLC)", cohort_group="Atlas_14", indication="NSCLC", sorting_regime="FACS_Immune_Sorted", sample_col="sample", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Lu_2022_HCC", dataset_id="GSE149614", display_name="Lu 2022 (HCC)", cohort_group="Atlas_14", indication="HCC", sorting_regime="Whole_Tumor", sample_col="sample", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Pelka_2021_CRC", dataset_id="GSE178341", display_name="Pelka 2021 (CRC)", cohort_group="Atlas_14", indication="CRC", sorting_regime="Whole_Tumor", sample_col="PatientTypeID", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Pu_2021_THCA", dataset_id="GSE184362", display_name="Pu 2021 (THCA)", cohort_group="Atlas_14", indication="THCA", sorting_regime="Whole_Tumor", sample_col="sample", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Qian_2020_PanCancer", dataset_id="E-MTAB-8107", display_name="Qian 2020 (Pan-Cancer)", cohort_group="Atlas_14", indication="Pan-Cancer", sorting_regime="Whole_Tumor", sample_col="Comment[ENA_SAMPLE]", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Sharma_2020_HCC", dataset_id="GSE156625", display_name="Sharma 2020 (HCC)", cohort_group="Atlas_14", indication="HCC", sorting_regime="Whole_Tumor", sample_col="patient_tumorsection", has_icb_response=False),
)


def _enrich_sc9_audited_response_metadata(cohort_id: str, meta_df: pl.DataFrame, pb_cache_dir: Path) -> pl.DataFrame:
    """Backfill audited SC9 response_binary & treatment_status from companion sample_meta parquet."""
    comp_path = pb_cache_dir / f"{cohort_id}.sample_meta.parquet"
    if not comp_path.is_file():
        return meta_df
    comp = pl.read_parquet(comp_path)
    if not {"sample_id", "response_binary"}.issubset(comp.columns):
        return meta_df
    rows = comp.to_dicts()
    lookup_r = {str(r["sample_id"]).strip(): float(r["response_binary"]) for r in rows if r["response_binary"] in (0, 1, 0.0, 1.0)}
    lookup_tp = {str(r["sample_id"]).strip(): str(r.get("treatment_status") or "Unknown") for r in rows}

    def _res_r(r: dict[str, object]) -> float:
        sid, pid = str(r["sample_id"]).strip(), str(r.get("patient_id") or r["sample_id"]).strip()
        raw_v = r["response_binary"]
        cur = float(raw_v) if isinstance(raw_v, (int, float)) else float("nan")
        if cur in (0.0, 1.0):
            return cur
        return lookup_r.get(sid, lookup_r[pid] if pid not in ("", "NA", "NAN") and pid in lookup_r else cur)

    def _res_tp(r: dict[str, object]) -> str:
        sid, pid = str(r["sample_id"]).strip(), str(r.get("patient_id") or r["sample_id"]).strip()
        cur = str(r.get("biopsy_timepoint") or "Unknown")
        return lookup_tp.get(sid, lookup_tp.get(pid, cur)) if cur == "Unknown" else cur

    meta_rows = meta_df.to_dicts()
    return meta_df.with_columns(
        [
            pl.Series("response_binary", [_res_r(r) for r in meta_rows]),
            pl.Series("biopsy_timepoint", [_res_tp(r) for r in meta_rows]),
        ]
    )


def _matched_tcga_ref(tcga_ref: pl.DataFrame, ind_val: str) -> pl.DataFrame:
    projs = INDICATION_TO_TCGA_PROJECTS.get(ind_val)
    sub_ref = tcga_ref.filter(pl.col("TCGA_project").is_in(list(projs))) if projs else tcga_ref
    return sub_ref if sub_ref.height >= 35 else tcga_ref


def _classify_cohort_vs_matched_tcga(
    knn_in: pl.DataFrame, tcga_ref: pl.DataFrame, default_indication: str, knn_cfg: BagaevKnnConfig
) -> tuple[pl.DataFrame, pl.DataFrame]:
    ind_vals = sorted({str(v) for v in knn_in["sample_indication"].to_list()})
    if default_indication == "Pan-Cancer" and len(ind_vals) > 1:
        parts = [
            res.unwrap()
            for ind_val in ind_vals
            for res in (classify_bagaev_knn(knn_in.filter(pl.col("sample_indication") == ind_val), _matched_tcga_ref(tcga_ref, ind_val), config=knn_cfg),)
            if isinstance(res, Success)
        ]
        if not parts:
            raise RuntimeError("Pan-Cancer k-NN classification produced no outputs.")
        return pl.concat(parts, how="vertical"), tcga_ref
    ref_c = _matched_tcga_ref(tcga_ref, default_indication)
    return unwrap_or_raise(classify_bagaev_knn(knn_in, ref_c, config=knn_cfg), "TCGA k-NN failed"), ref_c


def _safe_auc(df: pl.DataFrame, score_expr: pl.Expr) -> float:
    if df.is_empty():
        return float("nan")
    y = df["response_binary"].to_numpy().astype(np.int64)
    if int(np.sum(y == 1)) < 2 or int(np.sum(y == 0)) < 2:
        return float("nan")
    s = df.select(score_expr.alias("s"))["s"].to_numpy().astype(np.float64)
    return float(roc_auc_score(y, s))


def run_step03(
    out_dir: Path,
    preprocessed_dir: Path,
    ensembl_dir: Path,
    pb_cache_dir: Path,
    sc9_milo_dir: Path,
) -> None:
    """Execute Step 03 across all 23 single-cell cohorts (9 Milopy SC9 + 14 Reference Atlases)."""
    paths = get_data_paths()
    for p in (preprocessed_dir, ensembl_dir, pb_cache_dir, paths.raw_dir):
        assert_inside_data_root(p, paths.data_root)

    out_dir.mkdir(parents=True, exist_ok=True)
    pb_fges_cache_dir = out_dir / "sc_pseudobulk_29fges_cache"
    pb_fges_cache_dir.mkdir(parents=True, exist_ok=True)

    tcga_ref = unwrap_or_raise(load_tcga_bagaev_reference())
    ens_map = load_ensembl_symbol_map(ensembl_dir, fallback_h5ad=Some(preprocessed_dir / "Hugo-iAtlas.h5ad"))
    collection = get_bagaev_mfp_collection()
    knn_cfg = BagaevKnnConfig(k_neighbors=35, metric="correlation")

    stromal_markers = {"COL1A1", "COL1A2", "DCN", "LUM", "FAP", "ACTA2", "PECAM1", "VWF", "CDH5", "KDR"}
    immune_markers = {"PTPRC", "CD3D", "CD3E", "CD8A", "CD4", "NKG7", "CD68", "CD14", "HLA-DRA", "MS4A1"}

    cohort_summary_rows: list[dict[str, str | int | float]] = []
    all_sample_tables: list[pl.DataFrame] = []
    emp_corr_rows: list[dict[str, str | float]] = []

    for spec in SC23_SPECS:
        cache_file = pb_fges_cache_dir / f"{CACHE_VERSION_PREFIX}_{spec.cohort_id}_raw_fges.parquet"
        if cache_file.is_file():
            raw_cohort_df = pl.read_parquet(cache_file)
            print(f"[Cached] Loaded {CACHE_VERSION_PREFIX} 29-Fges for {spec.cohort_id:24s}: {raw_cohort_df.height} tumor samples")
        else:
            print(f"[Streaming] Aggregating tumor-only pseudobulks for {spec.cohort_id} ({spec.display_name})...")
            pb_spec = TumorPseudobulkSpec(
                cohort_id=spec.cohort_id, sample_col=spec.sample_col, default_indication=spec.indication, min_cells_per_sample=spec.min_cells_per_sample
            )
            pb_out = unwrap_or_raise(
                stream_cohort_tumor_pseudobulk(spec.dataset_id, pb_spec, ens_map),
                f"Pseudobulk aggregation failed for {spec.cohort_id}",
            )
            meta_df = (
                _enrich_sc9_audited_response_metadata(spec.cohort_id, pb_out.meta_df, pb_cache_dir)
                if spec.cohort_group == "SC9_Milopy"
                else pb_out.meta_df
            )
            sym_upper = [g.strip().upper() for g in pb_out.gene_symbols]
            st_idx = [i for i, g in enumerate(sym_upper) if g in stromal_markers]
            im_idx = [i for i, g in enumerate(sym_upper) if g in immune_markers]
            mean_st = np.mean(pb_out.log2_cpm[:, st_idx], axis=1) if st_idx else np.zeros(len(pb_out.sample_ids))
            mean_im = np.mean(pb_out.log2_cpm[:, im_idx], axis=1) if im_idx else np.zeros(len(pb_out.sample_ids))
            ssg_df = unwrap_or_raise(
                compute_ssgsea_matrix(pb_out.log2_cpm, pb_out.sample_ids, pb_out.gene_symbols, collection),
                f"ssGSEA failed on {spec.cohort_id}",
            )
            raw_cohort_df = meta_df.with_columns(
                [
                    pl.lit(spec.cohort_id).alias("cohort"),
                    pl.lit(spec.display_name).alias("display_name"),
                    pl.lit(spec.cohort_group).alias("cohort_group"),
                    pl.lit(spec.indication).alias("indication"),
                    pl.lit(spec.sorting_regime).alias("sorting_regime"),
                    pl.lit(pb_out.scale_mode).alias("expression_scale_mode"),
                    pl.lit(pb_out.n_cells_pre_filter).alias("n_cells_raw"),
                    pl.lit(pb_out.n_cells_post_tumor_filter).alias("n_cells_tumor_gated"),
                    pl.lit(len(set(sym_upper))).alias("n_genes_detected"),
                    pl.Series("raw_stromal_log2cpm", mean_st.tolist()),
                    pl.Series("raw_immune_log2cpm", mean_im.tolist()),
                ]
            ).hstack(ssg_df.drop("sample_id"))
            raw_cohort_df.write_parquet(cache_file)
            print(f"  -> {raw_cohort_df.height} tumor samples | {len(set(sym_upper))} genes | scale={pb_out.scale_mode} | cells={pb_out.n_cells_post_tumor_filter:,}/{pb_out.n_cells_pre_filter:,}")

        if spec.cohort_group == "SC9_Milopy":
            raw_cohort_df = _enrich_sc9_audited_response_metadata(spec.cohort_id, raw_cohort_df, pb_cache_dir)

        sc_res = median_mad_scale_signatures(raw_cohort_df, BAGAEV_MFP_SIGNATURE_ORDER, group_col=Nothing, clip=Some(4.0))
        if isinstance(sc_res, Failure):
            continue
        scaled_c = sc_res.unwrap().with_columns((pl.col("cohort") + pl.lit("__") + pl.col("sample_id").cast(pl.String)).alias("global_sample_id"))

        knn_in = scaled_c.drop("sample_id").rename({"global_sample_id": "sample_id"})
        knn_df, ref_sub = _classify_cohort_vs_matched_tcga(knn_in, tcga_ref, spec.indication, knn_cfg)
        knn_df = knn_df.rename({"sample_id": "global_sample_id"})

        louv_sil, louv_tau = float("nan"), float("nan")
        if knn_in.height >= 8:
            louv_res = cluster_bagaev_louvain(knn_in, config=BagaevLouvainConfig(seed=42))
            if isinstance(louv_res, Success):
                l_out = louv_res.unwrap()
                louv_tau = l_out.selected_threshold
                sel_row = l_out.threshold_sweep_df.filter(pl.col("threshold") == louv_tau)
                louv_sil = float(sel_row["silhouette_corr"][0]) if not sel_row.is_empty() else louv_sil
                knn_df = knn_df.join(
                    l_out.assignments_df.select([pl.col("sample_id").alias("global_sample_id"), pl.col("louvain_subtype")]),
                    on="global_sample_id",
                    how="left",
                )
        if "louvain_subtype" not in knn_df.columns:
            knn_df = knn_df.with_columns(pl.col("bagaev_subtype").alias("louvain_subtype"))

        sample_full = scaled_c.join(knn_df, on="global_sample_id", how="left")
        all_sample_tables.append(sample_full)

        cons_res = evaluate_fges_correlation_conservation(scaled_c, ref_sub) if scaled_c.height >= 4 else Failure("n<4")
        cons_metrics = cons_res.unwrap() if isinstance(cons_res, Success) else {}
        st_counts = {st_n: int(sample_full.filter(pl.col("bagaev_subtype") == st_n).height) for st_n in CANONICAL_SUBTYPE_ORDER}
        stromal_ratio = float(np.mean(sample_full["raw_stromal_log2cpm"].to_numpy()) / max(0.01, float(np.mean(sample_full["raw_immune_log2cpm"].to_numpy()))))

        cohort_summary_rows.append(
            {
                "cohort_id": spec.cohort_id,
                "display_name": spec.display_name,
                "cohort_group": spec.cohort_group,
                "indication": spec.indication,
                "sorting_regime": spec.sorting_regime,
                "n_samples": int(sample_full.height),
                "n_unique_patients": int(sample_full["patient_id"].n_unique()),
                "n_genes_detected": int(sample_full["n_genes_detected"][0]),
                "median_cells_per_sample": float(np.median(sample_full["n_cells"].to_numpy())),
                "stromal_to_immune_expr_ratio": stromal_ratio,
                "n_IE": st_counts["IE"],
                "n_IEF": st_counts["IE/F"],
                "n_F": st_counts["F"],
                "n_D": st_counts["D"],
                "n_distinct_subtypes": sum(1 for v in st_counts.values() if v > 0),
                "louvain_selected_tau": louv_tau,
                "louvain_silhouette_corr": louv_sil,
                "offdiag_pearson_r_vs_tcga": float(cons_metrics.get("offdiag_pearson_r_vs_tcga", float("nan"))),
                "offdiag_spearman_rho_vs_tcga": float(cons_metrics.get("offdiag_spearman_rho_vs_tcga", float("nan"))),
                "block_coherence_delta": float(cons_metrics.get("block_coherence_delta", float("nan"))),
            }
        )

        for f_col, sig_col, lineage_label in (
            ("emp_frac_T_NK", "T_cells", "T_NK_vs_Bagaev_T_cells"),
            ("emp_frac_CAF", "CAF", "CAF_vs_Bagaev_CAF"),
            ("emp_frac_Endothelial", "Endothelium", "Endothelial_vs_Bagaev_Endothelium"),
            ("emp_frac_Myeloid", "Macrophages", "Myeloid_vs_Bagaev_Macrophages"),
        ):
            sub_v = sample_full.filter(pl.col(f_col).is_finite() & pl.col(sig_col).is_finite())
            if sub_v.height >= 6:
                x_f, y_s = sub_v[f_col].to_numpy().astype(np.float64), sub_v[sig_col].to_numpy().astype(np.float64)
                if np.std(x_f) > 1e-6 and np.std(y_s) > 1e-6:
                    sp = stats.spearmanr(x_f, y_s)
                    emp_corr_rows.append(
                        {"cohort_id": spec.cohort_id, "sorting_regime": spec.sorting_regime, "pair": lineage_label, "n_samples": float(sub_v.height), "spearman_rho": float(sp.statistic), "p_value": float(sp.pvalue)}
                    )

    master_sc_df = pl.concat(all_sample_tables, how="diagonal_relaxed")
    write_table(master_sc_df, out_dir, "sc23_pseudobulk_bagaev_samples", write_tsv=False)
    summary_df = pl.DataFrame(cohort_summary_rows)
    write_table(summary_df, out_dir, "sc23_pseudobulk_cohort_summary")
    if emp_corr_rows:
        write_table(pl.DataFrame(emp_corr_rows), out_dir, "sc23_empirical_fraction_vs_fges_correlations")

    wt_r = summary_df.filter(pl.col("sorting_regime") == "Whole_Tumor")["offdiag_pearson_r_vs_tcga"].drop_nans().to_numpy()
    facs_r = summary_df.filter(pl.col("sorting_regime") == "FACS_Immune_Sorted")["offdiag_pearson_r_vs_tcga"].drop_nans().to_numpy()
    if len(wt_r) >= 2 and len(facs_r) >= 2:
        _, mwu_p = stats.mannwhitneyu(wt_r, facs_r, alternative="greater")
        write_table(
            pl.DataFrame([{"n_whole_tumor_cohorts": len(wt_r), "n_facs_sorted_cohorts": len(facs_r), "median_r_whole_tumor_vs_tcga": float(np.median(wt_r)), "median_r_facs_sorted_vs_tcga": float(np.median(facs_r)), "mannwhitney_one_sided_p": float(mwu_p)}]),
            out_dir,
            "sc23_sorting_regime_conservation_test",
            write_tsv=False,
        )

    sc9_all_resp = master_sc_df.filter((pl.col("cohort_group") == "SC9_Milopy") & pl.col("response_binary").cast(pl.Float64).is_in([0.0, 1.0]))
    sc9_uniq_pref = deduplicate_patient_profiles(sc9_all_resp, mode="prefer_pre_one_per_patient", patient_col="patient_id", weight_col=Some("n_cells"))
    sc9_strict_pre = deduplicate_patient_profiles(sc9_all_resp, mode="strict_pre", patient_col="patient_id", weight_col=Some("n_cells"))
    sc9_specs = [s for s in SC23_SPECS if s.cohort_group == "SC9_Milopy"]

    sc9_resp_tables: list[pl.DataFrame] = []
    for scope_name, scope_df in (
        ("Unique_Patients_PreferPre", sc9_uniq_pref),
        ("Strict_Pre_Only", sc9_strict_pre),
        ("All_Tumor_Biopsies", sc9_all_resp),
    ):
        if scope_df.is_empty():
            continue
        for stratum_name, sub_df in (
            ("SC9_All_9_Pooled", scope_df),
            ("SC9_WholeTumor_7_Pooled", scope_df.filter(pl.col("sorting_regime") == "Whole_Tumor")),
            ("SC9_FACS_ImmuneSorted_2_Pooled", scope_df.filter(pl.col("sorting_regime") == "FACS_Immune_Sorted")),
        ):
            for m_col, m_name in (("bagaev_subtype", "TCGA_kNN"), ("louvain_subtype", "Cohort_Louvain")):
                sc9_resp_tables.extend(eval_tagged_response(sub_df, subtype_col=m_col, cohort_name=stratum_name, sample_scope=scope_name, method=m_name, stratify_by_cohort=True))
        for spec in sc9_specs:
            c_sub = scope_df.filter(pl.col("cohort") == spec.cohort_id)
            if c_sub.height >= 4:
                sc9_resp_tables.extend(eval_tagged_response(c_sub, subtype_col="bagaev_subtype", cohort_name=spec.cohort_id, sample_scope=scope_name, method="TCGA_kNN", stratify_by_cohort=False))
    if sc9_resp_tables:
        write_table(pl.concat(sc9_resp_tables, how="diagonal_relaxed"), out_dir, "sc9_milopy_subtype_response_rates")

    ie_expr = pl.col("prob_IE") + 0.5 * pl.col("prob_IEF")
    imm_fib_expr = pl.col("antitumor_immune_score") - 0.5 * pl.col("stromal_fibrotic_score")
    per_cohort_aucs = [
        {
            "cohort": spec.cohort_id,
            "display_name": spec.display_name,
            "sorting_regime": spec.sorting_regime,
            "n_evaluable_patients": int(c_sub.height),
            "n_responders": int((c_sub["response_binary"].to_numpy().astype(np.int64) == 1).sum()),
            "auc_bagaev_knn_ie_prob": ie_auc,
            "auc_bagaev_immune_minus_fibrotic": _safe_auc(c_sub, imm_fib_expr),
            "auc_strict_pre_only": _safe_auc(sc9_strict_pre.filter(pl.col("cohort") == spec.cohort_id), imm_fib_expr),
        }
        for spec in sc9_specs
        if (c_sub := sc9_uniq_pref.filter(pl.col("cohort") == spec.cohort_id)).height >= 4
        and np.isfinite(ie_auc := _safe_auc(c_sub, ie_expr))
    ]
    if per_cohort_aucs:
        auc_df = pl.DataFrame(per_cohort_aucs)
        milo_da_path = sc9_milo_dir / "milopy_cell_state_da_all_resolutions.parquet"
        if milo_da_path.is_file():
            milo_df = pl.read_parquet(milo_da_path)
            coh_col = next((c for c in ("cohort", "cohort_id", "dataset_id") if c in milo_df.columns), "")
            lfc_col = next((c for c in ("logFC", "mean_logFC", "milo_logfc") if c in milo_df.columns), "")
            if coh_col and lfc_col:
                auc_df = auc_df.join(
                    milo_df.group_by(coh_col).agg(pl.col(lfc_col).abs().max().alias("max_abs_milopy_state_logFC")).rename({coh_col: "cohort"}),
                    on="cohort",
                    how="left",
                )
        write_table(auc_df, out_dir, "sc9_milopy_cohort_bagaev_aucs")

    print(f"Step 03 complete ({master_sc_df.height} tumor samples across {summary_df.height} cohorts; {sc9_uniq_pref.height} unique SC9 evaluable patients).")


def main() -> None:
    paths = get_data_paths()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=paths.output_dir / "bagaev_mfp_analysis")
    parser.add_argument("--preprocessed-dir", type=Path, default=paths.preprocessed_dir)
    parser.add_argument("--ensembl-dir", type=Path, default=paths.ensembl_dir)
    parser.add_argument("--pb-cache-dir", type=Path, default=paths.pseudobulk_dir)
    parser.add_argument("--sc9-milo-dir", type=Path, default=paths.output_dir / "sc9_milo_deconv_validation" / "04_milopy")
    args = parser.parse_args()
    run_step03(args.out_dir, args.preprocessed_dir, args.ensembl_dir, args.pb_cache_dir, args.sc9_milo_dir)


if __name__ == "__main__":
    main()
