#!/usr/bin/env python3
"""Step 03: Tumor-Only Bagaev 29-Fges Pseudobulk Clustering & Response Analysis Across 23 SC Cohorts.

Evaluates:
- Group A (9 Milopy ICB Single-Cell Cohorts, SC9):
  GSE120575, CELLxGENE_7b20c613, CELLxGENE_05a8c945, CELLxGENE_6f9de485,
  GSE207422, GSE243013, GSE233203, GSE200996, GSE316195.
- Group B (14 Curated Reference Single-Cell Tumor Atlases):
  Azizi 2018 (BRCA), Becker 2022 (COAD), Biermann 2022 (BrainMet), Borcherding 2021 (ccRCC),
  Cheng 2021 (Pan-Myeloid), Durante 2020 (UVM), Khaliq 2022 (CRC), Kim 2020 (LUAD),
  Leader 2021 (NSCLC), Lu 2022 (HCC), Pelka 2021 (CRC), Pu 2021 (THCA),
  Qian 2020 (Pan-Cancer), Sharma 2020 (HCC).

Uses upstreamed tumor-only H5AD pseudobulk streaming (`sc_pseudobulk_stream`), accurate
`ln(CP10k+1)` vs `log2(TPM+1)` scale inversion, Ensembl release-111 HGNC symbol lookup
(resolving `Durante_2020_UVM`), sample-matched TCGA classification for Pan-Cancer cohorts,
and patient-deduplicated ICB response evaluation across `SC9`.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Final, Literal, assert_never
import numpy as np
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Nothing, Some
from returns.result import Failure, Success
from scipy import stats
from sklearn.metrics import roc_auc_score  # type: ignore[import-untyped]

from tme_datasets import (
    BAGAEV_MFP_SIGNATURE_ORDER,
    TumorPseudobulkResult,
    TumorPseudobulkSpec,
    aggregate_gzip_h5ad_tumor_pseudobulk,
    aggregate_multi_file_tumor_pseudobulk,
    compute_ssgsea_matrix,
    deduplicate_patient_profiles,
    get_bagaev_mfp_collection,
    get_data_paths,
    load_ensembl_symbol_map,
    median_mad_scale_signatures,
    stream_tumor_pseudobulk_h5ad,
)
from tme_response import (
    CANONICAL_SUBTYPE_ORDER,
    INDICATION_TO_TCGA_PROJECTS,
    BagaevKnnConfig,
    BagaevLouvainConfig,
    classify_bagaev_knn,
    cluster_bagaev_louvain,
    evaluate_fges_correlation_conservation,
    evaluate_subtype_response_associations,
)

SourceKind = Literal["preprocessed_h5ad", "lair_h5ad", "lair_multi_h5ad", "lair_gz_h5ad"]
CACHE_VERSION_PREFIX: Final[str] = "v2_tumor_only"


class SingleCellCohortSpec(BaseModel):
    """Immutable declarative specification for each of the 23 single-cell cohorts."""

    model_config = ConfigDict(frozen=True)

    cohort_id: str
    display_name: str
    cohort_group: Literal["SC9_Milopy", "Atlas_14"]
    indication: str
    sorting_regime: Literal["Whole_Tumor", "FACS_Immune_Sorted"]
    source_kind: SourceKind
    rel_path: str
    sample_col: str
    has_icb_response: bool
    min_cells_per_sample: int = 25


SC23_SPECS: Final[tuple[SingleCellCohortSpec, ...]] = (
    # --- Group A: 9 Milopy ICB Single-Cell Cohorts (SC9) ---
    SingleCellCohortSpec(cohort_id="GSE120575", display_name="Sade-Feldman 2018 (Melanoma)", cohort_group="SC9_Milopy", indication="Melanoma", sorting_regime="FACS_Immune_Sorted", source_kind="preprocessed_h5ad", rel_path="GSE120575.h5ad", sample_col="patient_geo", has_icb_response=True),
    SingleCellCohortSpec(cohort_id="CELLxGENE_7b20c613", display_name="Gondal 2025 (Pan-Cancer)", cohort_group="SC9_Milopy", indication="Pan-Cancer", sorting_regime="Whole_Tumor", source_kind="preprocessed_h5ad", rel_path="CELLxGENE_7b20c613.h5ad", sample_col="donor_id_pre_post", has_icb_response=True),
    SingleCellCohortSpec(cohort_id="CELLxGENE_05a8c945", display_name="Marteau 2026 (CRC)", cohort_group="SC9_Milopy", indication="CRC", sorting_regime="Whole_Tumor", source_kind="preprocessed_h5ad", rel_path="CELLxGENE_05a8c945.h5ad", sample_col="sample_id", has_icb_response=True),
    SingleCellCohortSpec(cohort_id="CELLxGENE_6f9de485", display_name="Yan 2026 (Breast)", cohort_group="SC9_Milopy", indication="Breast", sorting_regime="Whole_Tumor", source_kind="preprocessed_h5ad", rel_path="CELLxGENE_6f9de485.h5ad", sample_col="sample_id", has_icb_response=True),
    SingleCellCohortSpec(cohort_id="GSE207422", display_name="Hu 2023 (NSCLC)", cohort_group="SC9_Milopy", indication="NSCLC", sorting_regime="Whole_Tumor", source_kind="preprocessed_h5ad", rel_path="GSE207422.h5ad", sample_col="sample_id", has_icb_response=True),
    SingleCellCohortSpec(cohort_id="GSE243013", display_name="Liu 2025 (NSCLC)", cohort_group="SC9_Milopy", indication="NSCLC", sorting_regime="Whole_Tumor", source_kind="preprocessed_h5ad", rel_path="GSE243013.h5ad", sample_col="sample_id", has_icb_response=True),
    SingleCellCohortSpec(cohort_id="GSE233203", display_name="Hong 2025 (NSCLC)", cohort_group="SC9_Milopy", indication="NSCLC", sorting_regime="Whole_Tumor", source_kind="preprocessed_h5ad", rel_path="GSE233203.h5ad", sample_col="sample_id", has_icb_response=True, min_cells_per_sample=20),
    SingleCellCohortSpec(cohort_id="GSE200996", display_name="Luoma 2022 (HNSCC)", cohort_group="SC9_Milopy", indication="HNSCC", sorting_regime="FACS_Immune_Sorted", source_kind="preprocessed_h5ad", rel_path="GSE200996.h5ad", sample_col="sample_id", has_icb_response=True),
    SingleCellCohortSpec(cohort_id="GSE316195", display_name="Raufi 2026 (PDAC)", cohort_group="SC9_Milopy", indication="PDAC", sorting_regime="Whole_Tumor", source_kind="preprocessed_h5ad", rel_path="GSE316195.h5ad", sample_col="sample_id", has_icb_response=True),
    # --- Group B: 14 Curated Reference Single-Cell Tumor Atlases ---
    SingleCellCohortSpec(cohort_id="Azizi_2018_BRCA", display_name="Azizi 2018 (Breast)", cohort_group="Atlas_14", indication="Breast", sorting_regime="FACS_Immune_Sorted", source_kind="lair_multi_h5ad", rel_path="DatasetscRNASeqSignature-AziziSingleCellMapDiverse2018Adata", sample_col="file_sample", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Becker_2022_COAD", display_name="Becker 2022 (COAD)", cohort_group="Atlas_14", indication="CRC", sorting_regime="Whole_Tumor", source_kind="preprocessed_h5ad", rel_path="GSE201349.h5ad", sample_col="sample", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Biermann_2022_BrainMet", display_name="Biermann 2022 (Mel BrainMet)", cohort_group="Atlas_14", indication="Melanoma", sorting_regime="Whole_Tumor", source_kind="preprocessed_h5ad", rel_path="GSE200218.h5ad", sample_col="orig.ident", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Borcherding_2021_ccRCC", display_name="Borcherding 2021 (ccRCC)", cohort_group="Atlas_14", indication="ccRCC", sorting_regime="FACS_Immune_Sorted", source_kind="preprocessed_h5ad", rel_path="GSE121638.h5ad", sample_col="sample", has_icb_response=False, min_cells_per_sample=20),
    SingleCellCohortSpec(cohort_id="Cheng_2021_PanMyeloid", display_name="Cheng 2021 (Pan-Myeloid)", cohort_group="Atlas_14", indication="Pan-Cancer", sorting_regime="FACS_Immune_Sorted", source_kind="lair_h5ad", rel_path="DatasetscRNASeqSignature-ChengPancancerSinglecellTranscriptional2021Adata/adata.h5ad", sample_col="patient", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Durante_2020_UVM", display_name="Durante 2020 (UVM)", cohort_group="Atlas_14", indication="UVM", sorting_regime="Whole_Tumor", source_kind="lair_h5ad", rel_path="DatasetscRNASeqSignature-DuranteSinglecellAnalysisReveals2020Adata/adata.h5ad", sample_col="sample", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Khaliq_2022_CRC", display_name="Khaliq 2022 (CRC)", cohort_group="Atlas_14", indication="CRC", sorting_regime="Whole_Tumor", source_kind="preprocessed_h5ad", rel_path="GSE200997.h5ad", sample_col="samples", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Kim_2020_LUAD", display_name="Kim 2020 (LUAD)", cohort_group="Atlas_14", indication="LUAD", sorting_regime="Whole_Tumor", source_kind="preprocessed_h5ad", rel_path="GSE131907.h5ad", sample_col="Sample", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Leader_2021_NSCLC", display_name="Leader 2021 (NSCLC)", cohort_group="Atlas_14", indication="NSCLC", sorting_regime="FACS_Immune_Sorted", source_kind="preprocessed_h5ad", rel_path="GSE154826.h5ad", sample_col="sample", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Lu_2022_HCC", display_name="Lu 2022 (HCC)", cohort_group="Atlas_14", indication="HCC", sorting_regime="Whole_Tumor", source_kind="lair_h5ad", rel_path="DatasetscRNASeqSignature-LuSinglecellAtlasMulticellular2022Adata/adata.h5ad", sample_col="sample", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Pelka_2021_CRC", display_name="Pelka 2021 (CRC)", cohort_group="Atlas_14", indication="CRC", sorting_regime="Whole_Tumor", source_kind="preprocessed_h5ad", rel_path="GSE178341.h5ad", sample_col="PatientTypeID", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Pu_2021_THCA", display_name="Pu 2021 (THCA)", cohort_group="Atlas_14", indication="THCA", sorting_regime="Whole_Tumor", source_kind="preprocessed_h5ad", rel_path="GSE184362.h5ad", sample_col="sample", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Qian_2020_PanCancer", display_name="Qian 2020 (Pan-Cancer)", cohort_group="Atlas_14", indication="Pan-Cancer", sorting_regime="Whole_Tumor", source_kind="lair_h5ad", rel_path="DatasetscRNASeqSignature-QianPancancerBlueprintHeterogeneous2020aAdata/adata.h5ad", sample_col="Comment[ENA_SAMPLE]", has_icb_response=False),
    SingleCellCohortSpec(cohort_id="Sharma_2020_HCC", display_name="Sharma 2020 (HCC)", cohort_group="Atlas_14", indication="HCC", sorting_regime="Whole_Tumor", source_kind="lair_gz_h5ad", rel_path="DatasetscRNASeqSignature-SharmaOncofetalReprogrammingEndothelial2020/GSE156625_HCCFscanpyobj.h5ad.gz", sample_col="patient_tumorsection", has_icb_response=False),
)


def _enrich_sc9_audited_response_metadata(
    cohort_id: str,
    meta_df: pl.DataFrame,
    pb_cache_dir: Path,
) -> pl.DataFrame:
    """Backfill audited SC9 response_binary & treatment_status from companion sample_meta parquet."""
    comp_path = pb_cache_dir / f"{cohort_id}.sample_meta.parquet"
    if not comp_path.is_file():
        return meta_df

    comp = pl.read_parquet(comp_path)
    if "sample_id" not in comp.columns or "response_binary" not in comp.columns:
        return meta_df

    lookup_r: dict[str, float] = {
        str(r["sample_id"]).strip(): float(r["response_binary"])
        for r in comp.iter_rows(named=True)
        if r["response_binary"] in (0, 1, 0.0, 1.0)
    }
    lookup_tp: dict[str, str] = {
        str(r["sample_id"]).strip(): str(r.get("treatment_status") or "Unknown")
        for r in comp.iter_rows(named=True)
    }
    new_r: list[float] = []
    new_tp: list[str] = []
    for r in meta_df.iter_rows(named=True):
        sid = str(r["sample_id"]).strip()
        pid = str(r.get("patient_id") or sid).strip()
        cur_r = float(r["response_binary"]) if r["response_binary"] is not None else float("nan")
        if cur_r not in (0.0, 1.0):
            if sid in lookup_r:
                cur_r = lookup_r[sid]
            elif pid in lookup_r and pid not in ("", "NA", "NAN"):
                cur_r = lookup_r[pid]
        new_r.append(cur_r)

        cur_tp = str(r.get("biopsy_timepoint") or "Unknown")
        if cur_tp == "Unknown":
            cur_tp = lookup_tp.get(sid, lookup_tp.get(pid, cur_tp))
        new_tp.append(cur_tp)

    return meta_df.with_columns(
        [
            pl.Series("response_binary", new_r),
            pl.Series("biopsy_timepoint", new_tp),
        ]
    )


def _aggregate_spec(
    spec: SingleCellCohortSpec,
    preprocessed_dir: Path,
    lair_dir: Path,
    ens_map: dict[str, str],
) -> TumorPseudobulkResult:
    """Dispatch tumor-only pseudobulk aggregation for a single-cell cohort specification."""
    pb_spec = TumorPseudobulkSpec(
        cohort_id=spec.cohort_id,
        sample_col=spec.sample_col,
        default_indication=spec.indication,
        min_cells_per_sample=spec.min_cells_per_sample,
    )
    match spec.source_kind:
        case "preprocessed_h5ad":
            res = stream_tumor_pseudobulk_h5ad(preprocessed_dir / spec.rel_path, pb_spec, ens_map)
        case "lair_h5ad":
            res = stream_tumor_pseudobulk_h5ad(lair_dir / spec.rel_path, pb_spec, ens_map)
        case "lair_multi_h5ad":
            res = aggregate_multi_file_tumor_pseudobulk(lair_dir / spec.rel_path, pb_spec, ens_map)
        case "lair_gz_h5ad":
            res = aggregate_gzip_h5ad_tumor_pseudobulk(lair_dir / spec.rel_path, pb_spec, ens_map)
        case _ as unreachable:
            assert_never(unreachable)

    if isinstance(res, Failure):
        raise RuntimeError(f"Pseudobulk aggregation failed for {spec.cohort_id}: {res.failure()}")
    return res.unwrap()


def _classify_cohort_vs_matched_tcga(
    knn_in: pl.DataFrame,
    tcga_ref: pl.DataFrame,
    default_indication: str,
    knn_cfg: BagaevKnnConfig,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Classify samples against indication-matched TCGA projects (per-sample for Pan-Cancer cohorts)."""
    ind_vals = sorted({str(v) for v in knn_in["sample_indication"].to_list()})
    if default_indication == "Pan-Cancer" and len(ind_vals) > 1:
        parts: list[pl.DataFrame] = []
        for ind_val in ind_vals:
            sub_q = knn_in.filter(pl.col("sample_indication") == ind_val)
            projs = INDICATION_TO_TCGA_PROJECTS.get(ind_val)
            sub_ref = tcga_ref.filter(pl.col("TCGA_project").is_in(list(projs))) if projs else tcga_ref
            if sub_ref.height < 35:
                sub_ref = tcga_ref
            res = classify_bagaev_knn(sub_q, sub_ref, config=knn_cfg)
            if isinstance(res, Success):
                parts.append(res.unwrap())
        if not parts:
            raise RuntimeError("Pan-Cancer k-NN classification produced no outputs.")
        return pl.concat(parts, how="vertical"), tcga_ref

    projs_c = INDICATION_TO_TCGA_PROJECTS.get(default_indication)
    ref_c = tcga_ref.filter(pl.col("TCGA_project").is_in(list(projs_c))) if projs_c else tcga_ref
    if ref_c.height < 35:
        ref_c = tcga_ref
    knn_res = classify_bagaev_knn(knn_in, ref_c, config=knn_cfg)
    if isinstance(knn_res, Failure):
        raise RuntimeError(f"TCGA k-NN failed: {knn_res.failure()}")
    return knn_res.unwrap(), ref_c


def run_step03(
    out_dir: Path,
    preprocessed_dir: Path,
    ensembl_dir: Path,
    lair_dir: Path,
    pb_cache_dir: Path,
    sc9_milo_dir: Path,
) -> None:
    """Execute Step 03 across all 23 single-cell cohorts (9 Milopy SC9 + 14 Reference Atlases)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    pb_fges_cache_dir = out_dir / "sc_pseudobulk_29fges_cache"
    pb_fges_cache_dir.mkdir(parents=True, exist_ok=True)

    tcga_ref = pl.read_parquet(out_dir / "tcga_bagaev_reference_29fges.parquet")
    ens_map = load_ensembl_symbol_map(
        ensembl_dir, fallback_h5ad=Some(preprocessed_dir / "Hugo-iAtlas.h5ad")
    )
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
            pb_out = _aggregate_spec(spec, preprocessed_dir, lair_dir, ens_map)
            meta_df = pb_out.meta_df
            if spec.cohort_group == "SC9_Milopy":
                meta_df = _enrich_sc9_audited_response_metadata(spec.cohort_id, meta_df, pb_cache_dir)

            sym_upper = [g.strip().upper() for g in pb_out.gene_symbols]
            st_idx = [i for i, g in enumerate(sym_upper) if g in stromal_markers]
            im_idx = [i for i, g in enumerate(sym_upper) if g in immune_markers]
            mean_stromal_expr = np.mean(pb_out.log2_cpm[:, st_idx], axis=1) if st_idx else np.zeros(len(pb_out.sample_ids))
            mean_immune_expr = np.mean(pb_out.log2_cpm[:, im_idx], axis=1) if im_idx else np.zeros(len(pb_out.sample_ids))

            ssg_res = compute_ssgsea_matrix(pb_out.log2_cpm, pb_out.sample_ids, pb_out.gene_symbols, collection)
            if isinstance(ssg_res, Failure):
                raise RuntimeError(f"ssGSEA failed on {spec.cohort_id}: {ssg_res.failure()}")

            raw_cohort_df = (
                meta_df.with_columns(
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
                        pl.Series("raw_stromal_log2cpm", mean_stromal_expr.tolist()),
                        pl.Series("raw_immune_log2cpm", mean_immune_expr.tolist()),
                    ]
                )
                .hstack(ssg_res.unwrap().drop("sample_id"))
            )
            raw_cohort_df.write_parquet(cache_file)
            print(
                f"  -> {raw_cohort_df.height} tumor samples | {len(set(sym_upper))} genes | scale={pb_out.scale_mode} "
                f"| cells={pb_out.n_cells_post_tumor_filter:,}/{pb_out.n_cells_pre_filter:,}"
            )

        if spec.cohort_group == "SC9_Milopy":
            raw_cohort_df = _enrich_sc9_audited_response_metadata(spec.cohort_id, raw_cohort_df, pb_cache_dir)

        sc_res = median_mad_scale_signatures(raw_cohort_df, BAGAEV_MFP_SIGNATURE_ORDER, group_col=Nothing, clip=Some(4.0))
        if isinstance(sc_res, Failure):
            continue
        scaled_c = sc_res.unwrap().with_columns(
            (pl.col("cohort") + pl.lit("__") + pl.col("sample_id").cast(pl.String)).alias("global_sample_id")
        )

        knn_in = scaled_c.drop("sample_id").rename({"global_sample_id": "sample_id"})
        knn_df, ref_sub = _classify_cohort_vs_matched_tcga(knn_in, tcga_ref, spec.indication, knn_cfg)
        knn_df = knn_df.rename({"sample_id": "global_sample_id"})

        louv_sil = float("nan")
        louv_tau = float("nan")
        if knn_in.height >= 8:
            louv_res = cluster_bagaev_louvain(knn_in, config=BagaevLouvainConfig(seed=42))
            if isinstance(louv_res, Success):
                l_out = louv_res.unwrap()
                louv_tau = l_out.selected_threshold
                sel_row = l_out.threshold_sweep_df.filter(pl.col("threshold") == louv_tau)
                if not sel_row.is_empty():
                    louv_sil = float(sel_row["silhouette_corr"][0])
                knn_df = knn_df.join(
                    l_out.assignments_df.select(
                        [pl.col("sample_id").alias("global_sample_id"), pl.col("louvain_subtype")]
                    ),
                    on="global_sample_id",
                    how="left",
                )
        if "louvain_subtype" not in knn_df.columns:
            knn_df = knn_df.with_columns(pl.col("bagaev_subtype").alias("louvain_subtype"))

        sample_full = scaled_c.join(knn_df, on="global_sample_id", how="left")
        all_sample_tables.append(sample_full)

        # Correlation conservation (NaN cleanly when n_samples < 4, e.g. Borcherding N=3 tumor samples)
        cons_metrics: dict[str, float] = {}
        if scaled_c.height >= 4:
            cons_res = evaluate_fges_correlation_conservation(scaled_c, ref_sub)
            if isinstance(cons_res, Success):
                cons_metrics = cons_res.unwrap()

        st_counts = {
            st_n: int(sample_full.filter(pl.col("bagaev_subtype") == st_n).height)
            for st_n in CANONICAL_SUBTYPE_ORDER
        }
        stromal_ratio = float(
            np.mean(sample_full["raw_stromal_log2cpm"].to_numpy())
            / max(0.01, float(np.mean(sample_full["raw_immune_log2cpm"].to_numpy())))
        )

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
            sub_valid = sample_full.filter(pl.col(f_col).is_finite() & pl.col(sig_col).is_finite())
            if sub_valid.height >= 6:
                x_f = sub_valid[f_col].to_numpy().astype(np.float64)
                y_s = sub_valid[sig_col].to_numpy().astype(np.float64)
                if np.std(x_f) > 1e-6 and np.std(y_s) > 1e-6:
                    sp_res = stats.spearmanr(x_f, y_s)
                    emp_corr_rows.append(
                        {
                            "cohort_id": spec.cohort_id,
                            "sorting_regime": spec.sorting_regime,
                            "pair": lineage_label,
                            "n_samples": float(sub_valid.height),
                            "spearman_rho": float(sp_res.statistic),
                            "p_value": float(sp_res.pvalue),
                        }
                    )

    master_sc_df = pl.concat(all_sample_tables, how="diagonal_relaxed")
    master_sc_df.write_parquet(out_dir / "sc23_pseudobulk_bagaev_samples.parquet")

    summary_df = pl.DataFrame(cohort_summary_rows)
    summary_df.write_parquet(out_dir / "sc23_pseudobulk_cohort_summary.parquet")
    summary_df.write_csv(out_dir / "sc23_pseudobulk_cohort_summary.tsv", separator="\t")

    if emp_corr_rows:
        emp_df = pl.DataFrame(emp_corr_rows)
        emp_df.write_parquet(out_dir / "sc23_empirical_fraction_vs_fges_correlations.parquet")
        emp_df.write_csv(out_dir / "sc23_empirical_fraction_vs_fges_correlations.tsv", separator="\t")

    wt_r = summary_df.filter(pl.col("sorting_regime") == "Whole_Tumor")["offdiag_pearson_r_vs_tcga"].drop_nans().to_numpy()
    facs_r = summary_df.filter(pl.col("sorting_regime") == "FACS_Immune_Sorted")["offdiag_pearson_r_vs_tcga"].drop_nans().to_numpy()
    if len(wt_r) >= 2 and len(facs_r) >= 2:
        _, mwu_sort_p = stats.mannwhitneyu(wt_r, facs_r, alternative="greater")
        pl.DataFrame(
            [
                {
                    "n_whole_tumor_cohorts": len(wt_r),
                    "n_facs_sorted_cohorts": len(facs_r),
                    "median_r_whole_tumor_vs_tcga": float(np.median(wt_r)),
                    "median_r_facs_sorted_vs_tcga": float(np.median(facs_r)),
                    "mannwhitney_one_sided_p": float(mwu_sort_p),
                }
            ]
        ).write_parquet(out_dir / "sc23_sorting_regime_conservation_test.parquet")

    # Evaluable SC9 clinical response comparisons: patient-deduplicated & strict pre sensitivity
    sc9_all_resp = master_sc_df.filter(
        (pl.col("cohort_group") == "SC9_Milopy") & pl.col("response_binary").cast(pl.Float64).is_in([0.0, 1.0])
    )
    sc9_uniq_pref = deduplicate_patient_profiles(
        sc9_all_resp, mode="prefer_pre_one_per_patient", patient_col="patient_id", weight_col=Some("n_cells")
    )
    sc9_strict_pre = deduplicate_patient_profiles(
        sc9_all_resp, mode="strict_pre", patient_col="patient_id", weight_col=Some("n_cells")
    )

    sc9_resp_tables: list[pl.DataFrame] = []
    for scope_name, scope_df in (
        ("Unique_Patients_PreferPre", sc9_uniq_pref),
        ("Strict_Pre_Only", sc9_strict_pre),
        ("All_Tumor_Biopsies", sc9_all_resp),
    ):
        if scope_df.is_empty():
            continue
        strata_defs = [
            ("SC9_All_9_Pooled", scope_df),
            ("SC9_WholeTumor_7_Pooled", scope_df.filter(pl.col("sorting_regime") == "Whole_Tumor")),
            ("SC9_FACS_ImmuneSorted_2_Pooled", scope_df.filter(pl.col("sorting_regime") == "FACS_Immune_Sorted")),
        ]
        for stratum_name, sub_df in strata_defs:
            for m_col, m_name in (("bagaev_subtype", "TCGA_kNN"), ("louvain_subtype", "Cohort_Louvain")):
                r_assoc = evaluate_subtype_response_associations(sub_df, subtype_col=m_col, cohort_col=Some("cohort"))
                if isinstance(r_assoc, Success):
                    sc9_resp_tables.append(
                        r_assoc.unwrap().with_columns(
                            [
                                pl.lit(stratum_name).alias("cohort"),
                                pl.lit(scope_name).alias("sample_scope"),
                                pl.lit(m_name).alias("method"),
                            ]
                        )
                    )
        for spec in SC23_SPECS:
            if spec.cohort_group != "SC9_Milopy":
                continue
            c_sub = scope_df.filter(pl.col("cohort") == spec.cohort_id)
            if c_sub.height < 4:
                continue
            c_assoc = evaluate_subtype_response_associations(c_sub, subtype_col="bagaev_subtype", cohort_col=Nothing)
            if isinstance(c_assoc, Success):
                sc9_resp_tables.append(
                    c_assoc.unwrap().with_columns(
                        [
                            pl.lit(spec.cohort_id).alias("cohort"),
                            pl.lit(scope_name).alias("sample_scope"),
                            pl.lit("TCGA_kNN").alias("method"),
                        ]
                    )
                )

    per_cohort_aucs: list[dict[str, str | float | int]] = []
    for spec in SC23_SPECS:
        if spec.cohort_group != "SC9_Milopy":
            continue
        c_sub = sc9_uniq_pref.filter(pl.col("cohort") == spec.cohort_id)
        if c_sub.height < 4:
            continue
        yc = c_sub["response_binary"].to_numpy().astype(np.int64)
        if int(np.sum(yc == 1)) >= 2 and int(np.sum(yc == 0)) >= 2:
            ie_score = (c_sub["prob_IE"] + 0.5 * c_sub["prob_IEF"]).to_numpy().astype(np.float64)
            imm_minus_fib = (c_sub["antitumor_immune_score"] - 0.5 * c_sub["stromal_fibrotic_score"]).to_numpy().astype(np.float64)

            c_pre = sc9_strict_pre.filter(pl.col("cohort") == spec.cohort_id)
            pre_auc = float("nan")
            if not c_pre.is_empty():
                y_pre = c_pre["response_binary"].to_numpy().astype(np.int64)
                if int(np.sum(y_pre == 1)) >= 2 and int(np.sum(y_pre == 0)) >= 2:
                    pre_score = (c_pre["antitumor_immune_score"] - 0.5 * c_pre["stromal_fibrotic_score"]).to_numpy().astype(np.float64)
                    pre_auc = float(roc_auc_score(y_pre, pre_score))

            per_cohort_aucs.append(
                {
                    "cohort": spec.cohort_id,
                    "display_name": spec.display_name,
                    "sorting_regime": spec.sorting_regime,
                    "n_evaluable_patients": int(c_sub.height),
                    "n_responders": int(np.sum(yc == 1)),
                    "auc_bagaev_knn_ie_prob": float(roc_auc_score(yc, ie_score)),
                    "auc_bagaev_immune_minus_fibrotic": float(roc_auc_score(yc, imm_minus_fib)),
                    "auc_strict_pre_only": pre_auc,
                }
            )

    if sc9_resp_tables:
        sc9_resp_df = pl.concat(sc9_resp_tables, how="diagonal_relaxed")
        sc9_resp_df.write_parquet(out_dir / "sc9_milopy_subtype_response_rates.parquet")
        sc9_resp_df.write_csv(out_dir / "sc9_milopy_subtype_response_rates.tsv", separator="\t")

    if per_cohort_aucs:
        auc_df = pl.DataFrame(per_cohort_aucs)
        milo_da_path = sc9_milo_dir / "milopy_cell_state_da_all_resolutions.parquet"
        if milo_da_path.is_file():
            milo_df = pl.read_parquet(milo_da_path)
            coh_col = next((c for c in ("cohort", "cohort_id", "dataset_id") if c in milo_df.columns), "")
            lfc_col = next((c for c in ("logFC", "mean_logFC", "milo_logfc") if c in milo_df.columns), "")
            if coh_col and lfc_col:
                milo_agg = (
                    milo_df.group_by(coh_col)
                    .agg(pl.col(lfc_col).abs().max().alias("max_abs_milopy_state_logFC"))
                    .rename({coh_col: "cohort"})
                )
                auc_df = auc_df.join(milo_agg, on="cohort", how="left")
        auc_df.write_parquet(out_dir / "sc9_milopy_cohort_bagaev_aucs.parquet")
        auc_df.write_csv(out_dir / "sc9_milopy_cohort_bagaev_aucs.tsv", separator="\t")

    print(
        f"Step 03 complete ({master_sc_df.height} tumor samples across {summary_df.height} cohorts; "
        f"{sc9_uniq_pref.height} unique SC9 evaluable patients)."
    )


def main() -> None:
    paths = get_data_paths()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=paths.output_dir / "bagaev_mfp_analysis")
    parser.add_argument("--preprocessed-dir", type=Path, default=paths.preprocessed_dir)
    parser.add_argument("--ensembl-dir", type=Path, default=paths.ensembl_dir)
    parser.add_argument("--lair-dir", type=Path, default=paths.lair_dir)
    parser.add_argument("--pb-cache-dir", type=Path, default=paths.pseudobulk_dir)
    parser.add_argument(
        "--sc9-milo-dir",
        type=Path,
        default=paths.output_dir / "sc9_milo_deconv_validation" / "04_milopy",
    )
    args = parser.parse_args()
    run_step03(
        out_dir=args.out_dir,
        preprocessed_dir=args.preprocessed_dir,
        ensembl_dir=args.ensembl_dir,
        lair_dir=args.lair_dir,
        pb_cache_dir=args.pb_cache_dir,
        sc9_milo_dir=args.sc9_milo_dir,
    )


if __name__ == "__main__":
    main()
