#!/usr/bin/env python3
"""Step 03: Bagaev 29-Fges Pseudobulk Clustering & Response Analysis Across 23 Single-Cell Cohorts.

Evaluates:
- Group A (9 Milopy ICB Single-Cell Cohorts, SC9):
  GSE120575, CELLxGENE_7b20c613, CELLxGENE_05a8c945, CELLxGENE_6f9de485,
  GSE207422, GSE243013, GSE233203, GSE200996, GSE316195.
- Group B (14 Curated Reference Single-Cell Tumor Atlases):
  Azizi 2018 (BRCA), Becker 2022 (COAD), Biermann 2022 (BrainMet), Borcherding 2021 (ccRCC),
  Cheng 2021 (Pan-Myeloid), Durante 2020 (UVM), Khaliq 2022 (CRC), Kim 2020 (LUAD),
  Leader 2021 (NSCLC), Lu 2022 (HCC), Pelka 2021 (CRC), Pu 2021 (THCA),
  Qian 2020 (Pan-Cancer), Sharma 2020 (HCC).

Streams full-transcriptome patient/sample pseudobulks in memory-bounded CSR chunks, normalizes to
log2(CPM_1e6 + 1), scores all 29 Bagaev Fges via ssGSEA, compares 29x29 signature correlation
conservation vs TCGA (Whole-Tumor vs FACS Immune-Sorted), runs TCGA K=35 k-NN & Louvain clustering,
validates against empirical single-cell fractions, and tests ICB response correlation across SC9.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import gzip
from pathlib import Path
import shutil
import tempfile
from typing import Any, Final
import anndata as ad
from anndata.io import read_elem
import h5py
import numpy as np
import numpy.typing as npt
import pandas as pd
import polars as pl
from returns.maybe import Nothing, Some
from returns.result import Failure, Success
from scipy import stats
import scipy.sparse as sp
from sklearn.metrics import roc_auc_score

from tme_datasets import (
    BAGAEV_MFP_SIGNATURE_ORDER,
    compute_ssgsea_matrix,
    get_bagaev_mfp_collection,
    get_data_paths,
    median_mad_scale_signatures,
)
from tme_datasets.preprocessing.metadata import binarize_response, standardize_timepoint
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


@dataclass(frozen=True)
class SingleCellCohortSpec:
    """Declarative specification for each of the 23 single-cell cohorts."""

    cohort_id: str
    display_name: str
    cohort_group: str  # 'SC9_Milopy' or 'Atlas_14'
    indication: str
    sorting_regime: str  # 'Whole_Tumor' or 'FACS_Immune_Sorted'
    source_kind: str  # 'preprocessed_h5ad', 'lair_h5ad', 'lair_multi_h5ad', 'lair_gz_h5ad'
    rel_path: str
    sample_col: str
    has_icb_response: bool
    min_cells_per_sample: int = 25


SC23_SPECS: Final[tuple[SingleCellCohortSpec, ...]] = (
    # --- Group A: 9 Milopy ICB Single-Cell Cohorts (SC9) ---
    SingleCellCohortSpec("GSE120575", "Sade-Feldman 2018 (Melanoma)", "SC9_Milopy", "Melanoma", "FACS_Immune_Sorted", "preprocessed_h5ad", "GSE120575.h5ad", "patient_geo", True),
    SingleCellCohortSpec("CELLxGENE_7b20c613", "Gondal 2025 (Pan-Ca/Mel)", "SC9_Milopy", "Melanoma", "Whole_Tumor", "preprocessed_h5ad", "CELLxGENE_7b20c613.h5ad", "PMID_donor_id", True),
    SingleCellCohortSpec("CELLxGENE_05a8c945", "Marteau 2026 (CRC)", "SC9_Milopy", "CRC", "Whole_Tumor", "preprocessed_h5ad", "CELLxGENE_05a8c945.h5ad", "sample_id", True),
    SingleCellCohortSpec("CELLxGENE_6f9de485", "Yan 2026 (Breast)", "SC9_Milopy", "Breast", "Whole_Tumor", "preprocessed_h5ad", "CELLxGENE_6f9de485.h5ad", "sample_id", True),
    SingleCellCohortSpec("GSE207422", "Hu 2023 (NSCLC)", "SC9_Milopy", "NSCLC", "Whole_Tumor", "preprocessed_h5ad", "GSE207422.h5ad", "sample_id", True),
    SingleCellCohortSpec("GSE243013", "Liu 2025 (NSCLC)", "SC9_Milopy", "NSCLC", "Whole_Tumor", "preprocessed_h5ad", "GSE243013.h5ad", "sample_id", True),
    SingleCellCohortSpec("GSE233203", "Hong 2025 (NSCLC)", "SC9_Milopy", "NSCLC", "Whole_Tumor", "preprocessed_h5ad", "GSE233203.h5ad", "sample_id", True, min_cells_per_sample=20),
    SingleCellCohortSpec("GSE200996", "Luoma 2022 (HNSCC)", "SC9_Milopy", "HNSCC", "FACS_Immune_Sorted", "preprocessed_h5ad", "GSE200996.h5ad", "sample_id", True),
    SingleCellCohortSpec("GSE316195", "Raufi 2026 (PDAC)", "SC9_Milopy", "PDAC", "Whole_Tumor", "preprocessed_h5ad", "GSE316195.h5ad", "sample_id", True),
    # --- Group B: 14 Curated Reference Single-Cell Tumor Atlases ---
    SingleCellCohortSpec("Azizi_2018_BRCA", "Azizi 2018 (Breast)", "Atlas_14", "Breast", "FACS_Immune_Sorted", "lair_multi_h5ad", "DatasetscRNASeqSignature-AziziSingleCellMapDiverse2018Adata", "file_sample", False),
    SingleCellCohortSpec("Becker_2022_COAD", "Becker 2022 (COAD)", "Atlas_14", "CRC", "Whole_Tumor", "preprocessed_h5ad", "GSE201349.h5ad", "sample", False),
    SingleCellCohortSpec("Biermann_2022_BrainMet", "Biermann 2022 (Mel BrainMet)", "Atlas_14", "Melanoma", "Whole_Tumor", "preprocessed_h5ad", "GSE200218.h5ad", "orig.ident", False),
    SingleCellCohortSpec("Borcherding_2021_ccRCC", "Borcherding 2021 (ccRCC)", "Atlas_14", "ccRCC", "FACS_Immune_Sorted", "preprocessed_h5ad", "GSE121638.h5ad", "sample", False, min_cells_per_sample=20),
    SingleCellCohortSpec("Cheng_2021_PanMyeloid", "Cheng 2021 (Pan-Myeloid)", "Atlas_14", "Pan-Cancer", "FACS_Immune_Sorted", "lair_h5ad", "DatasetscRNASeqSignature-ChengPancancerSinglecellTranscriptional2021Adata/adata.h5ad", "patient", False),
    SingleCellCohortSpec("Durante_2020_UVM", "Durante 2020 (UVM)", "Atlas_14", "UVM", "Whole_Tumor", "lair_h5ad", "DatasetscRNASeqSignature-DuranteSinglecellAnalysisReveals2020Adata/adata.h5ad", "sample", False),
    SingleCellCohortSpec("Khaliq_2022_CRC", "Khaliq 2022 (CRC)", "Atlas_14", "CRC", "Whole_Tumor", "preprocessed_h5ad", "GSE200997.h5ad", "samples", False),
    SingleCellCohortSpec("Kim_2020_LUAD", "Kim 2020 (LUAD)", "Atlas_14", "LUAD", "Whole_Tumor", "preprocessed_h5ad", "GSE131907.h5ad", "Sample", False),
    SingleCellCohortSpec("Leader_2021_NSCLC", "Leader 2021 (NSCLC)", "Atlas_14", "NSCLC", "FACS_Immune_Sorted", "preprocessed_h5ad", "GSE154826.h5ad", "sample", False),
    SingleCellCohortSpec("Lu_2022_HCC", "Lu 2022 (HCC)", "Atlas_14", "HCC", "Whole_Tumor", "lair_h5ad", "DatasetscRNASeqSignature-LuSinglecellAtlasMulticellular2022Adata/adata.h5ad", "sample", False),
    SingleCellCohortSpec("Pelka_2021_CRC", "Pelka 2021 (CRC)", "Atlas_14", "CRC", "Whole_Tumor", "preprocessed_h5ad", "GSE178341.h5ad", "PatientTypeID", False),
    SingleCellCohortSpec("Pu_2021_THCA", "Pu 2021 (THCA)", "Atlas_14", "THCA", "Whole_Tumor", "preprocessed_h5ad", "GSE184362.h5ad", "sample", False),
    SingleCellCohortSpec("Qian_2020_PanCancer", "Qian 2020 (Pan-Cancer)", "Atlas_14", "Pan-Cancer", "Whole_Tumor", "lair_h5ad", "DatasetscRNASeqSignature-QianPancancerBlueprintHeterogeneous2020aAdata/adata.h5ad", "Comment[ENA_SAMPLE]", False),
    SingleCellCohortSpec("Sharma_2020_HCC", "Sharma 2020 (HCC)", "Atlas_14", "HCC", "Whole_Tumor", "lair_gz_h5ad", "DatasetscRNASeqSignature-SharmaOncofetalReprogrammingEndothelial2020/GSE156625_HCCFscanpyobj.h5ad.gz", "patient_tumorsection", False),
)

RESP_CANDIDATE_COLS: Final[tuple[str, ...]] = (
    "response_binary",
    "clinical_response",
    "response",
    "outcome",
    "Combined_outcome",
    "Path_response",
    "RECIST",
    "therapeutic_response",
)

TIMEPOINT_CANDIDATE_COLS: Final[tuple[str, ...]] = (
    "biopsy_timepoint",
    "treatment_status",
    "timepoint",
    "time_point",
    "pre_post",
)

CELLTYPE_CANDIDATE_COLS: Final[tuple[str, ...]] = (
    "cell_type",
    "cell_type_author",
    "cell_type_main",
    "celltype",
    "MajorCluster",
    "SINGLECELL_TYPE",
)


_ENSEMBL_TO_SYMBOL_CACHE: dict[str, str] = {}


def _get_ensembl_to_symbol_map(preprocessed_dir: Path) -> dict[str, str]:
    """Load cached Ensembl-ID-to-HGNC-symbol dictionary from preprocessed reference .var table."""
    if _ENSEMBL_TO_SYMBOL_CACHE:
        return _ENSEMBL_TO_SYMBOL_CACHE
    ref_h5 = preprocessed_dir / "Hugo-iAtlas.h5ad"
    if ref_h5.is_file():
        with h5py.File(ref_h5, "r") as f:
            vdf: pd.DataFrame = read_elem(f["var"])
            col_sym = next((c for c in ("gene_name", "original_id", "gene_symbol") if c in vdf.columns), None)
            if col_sym is not None:
                for eid, sym in zip(vdf.index.tolist(), vdf[col_sym].tolist(), strict=False):
                    clean_e = str(eid).split(".")[0].strip().upper()
                    clean_s = str(sym).strip().upper()
                    if clean_e.startswith("ENSG") and clean_s and not clean_s.startswith("ENSG"):
                        _ENSEMBL_TO_SYMBOL_CACHE[clean_e] = clean_s
    return _ENSEMBL_TO_SYMBOL_CACHE


def _extract_var_gene_symbols(var_df: pd.DataFrame) -> list[str]:
    """Resolve Hugo gene symbols from an AnnData .var DataFrame (including bare ENSG indices)."""
    for col in ("gene_name", "gene_symbol", "gene_symbols", "feature_name", "symbol"):
        if col in var_df.columns:
            vals = [str(v).strip().upper() for v in var_df[col].tolist()]
            if any(v and not v.startswith("ENSG") and v != "NAN" for v in vals[:50]):
                return vals
    idx_vals = [str(v).strip().upper() for v in var_df.index.tolist()]
    if "original_id" in var_df.columns and idx_vals and idx_vals[0].startswith("ENSG"):
        orig_vals = [str(v).strip().upper() for v in var_df["original_id"].tolist()]
        if any(v and not v.startswith("ENSG") for v in orig_vals[:50]):
            return orig_vals
    if idx_vals and sum(1 for v in idx_vals[:50] if v.startswith("ENSG")) >= 25:
        ens_map = _get_ensembl_to_symbol_map(get_data_paths().preprocessed_dir)
        if ens_map:
            return [ens_map.get(v.split(".")[0], v) for v in idx_vals]
    return idx_vals


def _classify_broad_lineage(ct_raw: str) -> str:
    """Classify heterogeneous single-cell author labels into broad TME compartments."""
    s = ct_raw.strip().lower()
    if any(k in s for k in ("fibro", "caf", "stroma", "mesench", "pericyte", "smooth muscle", "stellate")):
        return "Stromal_CAF"
    if any(k in s for k in ("endothel", "vascular", "lymphatic")):
        return "Endothelial"
    if any(k in s for k in ("cd8", "cd4", "t cell", "t_cell", "tcell", "treg", "nk", "lympho", "cytotox")):
        return "T_NK_Lymphoid"
    if any(k in s for k in ("b cell", "b_cell", "bcell", "plasma")):
        return "B_Plasma"
    if any(k in s for k in ("macro", "mono", "myeloid", "dc", "dendrit", "neutro", "mast", "tam", "granulo")):
        return "Myeloid"
    if any(k in s for k in ("malig", "tumor", "cancer", "epithel", "melanom", "hepatoc", "carcinoma", "neoplast")):
        return "Malignant_Epithelial"
    return "Other"


def _stream_pseudobulk_from_h5ad(
    h5_path: Path,
    sample_col: str,
    min_cells: int = 25,
    chunk_rows: int = 50000,
) -> tuple[npt.NDArray[np.float64], list[str], list[str], pl.DataFrame]:
    """Stream CSR/CSC/dense H5AD in bounded memory to build log2(CPM_1e6 + 1) sample pseudobulks."""
    with h5py.File(h5_path, "r") as f:
        obs: pd.DataFrame = read_elem(f["obs"])
        var: pd.DataFrame = read_elem(f["var"])
        gene_symbols = _extract_var_gene_symbols(var)
        n_cells = len(obs)
        n_genes = len(gene_symbols)

        actual_scol = sample_col
        if actual_scol not in obs.columns:
            actual_scol = next(
                (c for c in ("sample_id", "sample", "samples", "patient_id", "patient", "donor_id", "orig.ident") if c in obs.columns),
                obs.columns[0],
            )

        raw_sids = [str(v).strip() for v in obs[actual_scol].tolist()]
        valid_cell_mask = np.array(
            [bool(s) and s.lower() not in ("nan", "none", "unknown", "na") for s in raw_sids],
            dtype=bool,
        )

        # Pick counts layer if available, otherwise X
        x_node: Any = f["layers"]["counts"] if ("layers" in f and "counts" in f["layers"]) else f["X"]
        is_csr = isinstance(x_node, h5py.Group) and "indptr" in x_node

        # Determine if X is log-normalized (only when reading /X directly)
        needs_expm1 = False
        is_raw_int = True
        if not ("layers" in f and "counts" in f["layers"]):
            sample_data = (
                x_node["data"][: min(50000, len(x_node["data"]))]
                if is_csr and len(x_node["data"]) > 0
                else np.asarray(x_node[: min(200, n_cells)]).ravel()
            )
            nz = sample_data[sample_data != 0]
            if nz.size > 0:
                max_v = float(np.nanmax(nz))
                is_raw_int = bool(np.allclose(nz[:5000], np.round(nz[:5000]), atol=1e-3))
                if not is_raw_int and max_v <= 30.0:
                    needs_expm1 = True

        # Map valid samples with >= min_cells
        s_arr = np.asarray(raw_sids, dtype=object)
        u_sids, u_cnts = np.unique(s_arr[valid_cell_mask], return_counts=True)
        kept_samples = sorted(str(s) for s, c in zip(u_sids.tolist(), u_cnts.tolist(), strict=True) if int(c) >= min_cells)
        if not kept_samples:
            raise ValueError(f"No samples with >= {min_cells} cells in {h5_path.name} (col={actual_scol})")

        sid_to_idx = {s: i for i, s in enumerate(kept_samples)}
        cell_sample_idx = np.array([sid_to_idx.get(str(s), -1) for s in raw_sids], dtype=np.int64)

        pb_sums = np.zeros((len(kept_samples), n_genes), dtype=np.float64)
        kept_cell_counts = np.zeros(len(kept_samples), dtype=np.int64)

        if is_csr:
            indptr = np.asarray(x_node["indptr"][:], dtype=np.int64)
            is_csc = len(indptr) == (n_genes + 1) and len(indptr) != (n_cells + 1)
            if is_csc:
                mat_full: sp.spmatrix = read_elem(x_node).tocsr()
                for start in range(0, n_cells, chunk_rows):
                    end = min(start + chunk_rows, n_cells)
                    sub_csr = mat_full[start:end].tocsr()
                    if needs_expm1 and sub_csr.nnz > 0:
                        sub_csr = sp.csr_matrix(
                            (np.expm1(np.clip(sub_csr.data.astype(np.float64), 0.0, 25.0)), sub_csr.indices, sub_csr.indptr),
                            shape=sub_csr.shape,
                        )
                    row_sums = np.asarray(sub_csr.sum(axis=1)).ravel()
                    min_umi = 150.0 if is_raw_int else 1.0
                    c_sidx = cell_sample_idx[start:end]
                    ok = (c_sidx >= 0) & (row_sums >= min_umi)
                    if not np.any(ok):
                        continue
                    valid_local = np.flatnonzero(ok)
                    grp_ids = c_sidx[valid_local]
                    indicator = sp.csr_matrix(
                        (np.ones(len(valid_local), dtype=np.float64), (grp_ids, valid_local)),
                        shape=(len(kept_samples), end - start),
                    )
                    pb_sums += np.asarray((indicator @ sub_csr).toarray(), dtype=np.float64)
                    np.add.at(kept_cell_counts, grp_ids, 1)
            else:
                indices_ds = x_node["indices"]
                data_ds = x_node["data"]
                for start in range(0, n_cells, chunk_rows):
                    end = min(start + chunk_rows, n_cells)
                    c_sidx = cell_sample_idx[start:end]
                    if not np.any(c_sidx >= 0):
                        continue
                    p0 = int(indptr[start])
                    p1 = int(indptr[end])
                    if p1 <= p0:
                        continue
                    sub_indptr = indptr[start : end + 1] - p0
                    sub_indices = np.asarray(indices_ds[p0:p1], dtype=np.int32)
                    sub_data = np.asarray(data_ds[p0:p1], dtype=np.float64)
                    if needs_expm1:
                        sub_data = np.expm1(np.clip(sub_data, 0.0, 25.0))
                    sub_csr = sp.csr_matrix((sub_data, sub_indices, sub_indptr), shape=(end - start, n_genes))
                    row_sums = np.asarray(sub_csr.sum(axis=1)).ravel()
                    min_umi = 150.0 if is_raw_int else 1.0
                    ok = (c_sidx >= 0) & (row_sums >= min_umi)
                    if not np.any(ok):
                        continue
                    valid_local = np.flatnonzero(ok)
                    grp_ids = c_sidx[valid_local]
                    indicator = sp.csr_matrix(
                        (np.ones(len(valid_local), dtype=np.float64), (grp_ids, valid_local)),
                        shape=(len(kept_samples), end - start),
                    )
                    pb_sums += np.asarray((indicator @ sub_csr).toarray(), dtype=np.float64)
                    np.add.at(kept_cell_counts, grp_ids, 1)
        else:
            for start in range(0, n_cells, chunk_rows):
                end = min(start + chunk_rows, n_cells)
                c_sidx = cell_sample_idx[start:end]
                if not np.any(c_sidx >= 0):
                    continue
                block = np.asarray(x_node[start:end], dtype=np.float64)
                if needs_expm1:
                    block = np.expm1(np.clip(block, 0.0, 25.0))
                row_sums = np.nansum(block, axis=1)
                min_umi = 150.0 if is_raw_int else 1.0
                ok = (c_sidx >= 0) & (row_sums >= min_umi)
                for local_r in np.flatnonzero(ok):
                    g_id = int(c_sidx[local_r])
                    pb_sums[g_id] += np.nan_to_num(block[local_r], nan=0.0)
                    kept_cell_counts[g_id] += 1

    # Filter out samples that fell below min_cells after UMI quality filtering
    non_empty = kept_cell_counts >= max(10, min_cells // 2)
    pb_sums = pb_sums[non_empty]
    kept_cell_counts = kept_cell_counts[non_empty]
    final_samples = [s for s, keep in zip(kept_samples, non_empty.tolist(), strict=True) if keep]

    # Normalize each sample pseudobulk to CPM = 10^6 and transform via log2(CPM + 1)
    lib_sizes = np.sum(pb_sums, axis=1, keepdims=True)
    lib_sizes = np.where(lib_sizes > 0.0, lib_sizes, 1.0)
    cpm = (pb_sums / lib_sizes) * 1e6
    log2_cpm = np.log2(cpm + 1.0)

    # Extract sample-level clinical response, timepoint, and empirical lineage fractions from obs
    resp_col = next((c for c in RESP_CANDIDATE_COLS if c in obs.columns), None)
    tp_col = next((c for c in TIMEPOINT_CANDIDATE_COLS if c in obs.columns), None)
    ct_col = next((c for c in CELLTYPE_CANDIDATE_COLS if c in obs.columns), None)

    meta_rows: list[dict[str, str | int | float]] = []
    obs_sid_series = pd.Series(raw_sids, index=obs.index)
    for s_idx, sid in enumerate(final_samples):
        sub_obs = obs.loc[(obs_sid_series == sid).to_numpy()]
        r_bin = float("nan")
        if resp_col is not None and len(sub_obs) > 0:
            vals_r = [binarize_response(v) for v in sub_obs[resp_col].dropna().unique().tolist()]
            finite_r = [v for v in vals_r if not np.isnan(v)]
            if finite_r:
                r_bin = float(np.round(np.mean(finite_r)))

        tp_str = "Unknown"
        if tp_col is not None and len(sub_obs) > 0:
            tps = [standardize_timepoint(v) for v in sub_obs[tp_col].dropna().unique().tolist()]
            known_tps = [t for t in tps if t != "Unknown"]
            tp_str = known_tps[0] if known_tps else "Unknown"

        frac_t = float("nan")
        frac_mye = float("nan")
        frac_caf = float("nan")
        frac_endo = float("nan")
        frac_mal = float("nan")
        if ct_col is not None and len(sub_obs) > 0:
            cats = [_classify_broad_lineage(str(v)) for v in sub_obs[ct_col].tolist()]
            n_tot_c = float(len(cats))
            frac_t = sum(1 for c in cats if c == "T_NK_Lymphoid") / n_tot_c
            frac_mye = sum(1 for c in cats if c == "Myeloid") / n_tot_c
            frac_caf = sum(1 for c in cats if c == "Stromal_CAF") / n_tot_c
            frac_endo = sum(1 for c in cats if c == "Endothelial") / n_tot_c
            frac_mal = sum(1 for c in cats if c == "Malignant_Epithelial") / n_tot_c

        meta_rows.append(
            {
                "sample_id": sid,
                "n_cells": int(kept_cell_counts[s_idx]),
                "response_binary": r_bin,
                "biopsy_timepoint": tp_str,
                "emp_frac_T_NK": frac_t,
                "emp_frac_Myeloid": frac_mye,
                "emp_frac_CAF": frac_caf,
                "emp_frac_Endothelial": frac_endo,
                "emp_frac_Malignant": frac_mal,
            }
        )

    return log2_cpm, final_samples, gene_symbols, pl.DataFrame(meta_rows)


def _aggregate_azizi_multi_h5ad(lair_subdir: Path) -> tuple[npt.NDArray[np.float64], list[str], list[str], pl.DataFrame]:
    """Aggregate Azizi 2018 Breast tumor sample H5AD files into a unified log2(CPM+1) pseudobulk matrix."""
    tumor_files = sorted(p for p in lair_subdir.glob("*TUMOR*.h5ad") if not p.name.startswith("._"))
    if not tumor_files:
        tumor_files = sorted(p for p in lair_subdir.glob("*.h5ad") if not p.name.startswith("._"))

    sample_dicts: list[dict[str, float]] = []
    sample_names: list[str] = []
    cell_counts: list[int] = []
    all_genes_set: set[str] = set()

    for fpath in tumor_files:
        a = ad.read_h5ad(fpath)
        if a.n_obs < 25:
            continue
        syms = _extract_var_gene_symbols(a.var)
        X = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
        col_sums = np.asarray(np.sum(X, axis=0), dtype=np.float64).ravel()
        tot = float(np.sum(col_sums))
        if tot <= 0:
            continue
        cpm_vec = (col_sums / tot) * 1e6
        s_dict: dict[str, float] = {}
        for g, val in zip(syms, cpm_vec.tolist(), strict=True):
            if g and g != "NAN" and val > 0.0:
                s_dict[g] = max(s_dict.get(g, 0.0), float(val))
        all_genes_set.update(s_dict.keys())
        sample_dicts.append(s_dict)
        sample_names.append(fpath.stem.replace("_counts", "").replace("_matrix", ""))
        cell_counts.append(int(a.n_obs))

    genes_sorted = sorted(all_genes_set)
    g_idx = {g: i for i, g in enumerate(genes_sorted)}
    mat = np.zeros((len(sample_names), len(genes_sorted)), dtype=np.float64)
    for r_i, s_dict in enumerate(sample_dicts):
        for g, v in s_dict.items():
            mat[r_i, g_idx[g]] = np.log2(v + 1.0)

    meta_df = pl.DataFrame(
        {
            "sample_id": sample_names,
            "n_cells": cell_counts,
            "response_binary": [float("nan")] * len(sample_names),
            "biopsy_timepoint": ["Unknown"] * len(sample_names),
            "emp_frac_T_NK": [float("nan")] * len(sample_names),
            "emp_frac_Myeloid": [float("nan")] * len(sample_names),
            "emp_frac_CAF": [float("nan")] * len(sample_names),
            "emp_frac_Endothelial": [float("nan")] * len(sample_names),
            "emp_frac_Malignant": [float("nan")] * len(sample_names),
        }
    )
    return mat, sample_names, genes_sorted, meta_df


def _aggregate_sharma_gz_h5ad(gz_path: Path, sample_col: str) -> tuple[npt.NDArray[np.float64], list[str], list[str], pl.DataFrame]:
    """Decompress Sharma 2020 HCC scanpy H5AD (using .raw if broader) and aggregate tumor samples."""
    with tempfile.NamedTemporaryFile(suffix=".h5ad") as tmp:
        with gzip.open(gz_path, "rb") as gi:
            shutil.copyfileobj(gi, tmp)
        tmp.flush()
        adata = ad.read_h5ad(tmp.name)

    if "NormalvsTumor" in adata.obs.columns:
        tumor_mask = adata.obs["NormalvsTumor"].astype(str).isin(["Tumor", "F1", "F2"]).to_numpy()
        if int(np.sum(tumor_mask)) >= 500:
            adata = adata[tumor_mask].copy()

    src_X = adata.raw.X if (adata.raw is not None and adata.raw.n_vars > adata.n_vars) else adata.X
    src_var = adata.raw.var if (adata.raw is not None and adata.raw.n_vars > adata.n_vars) else adata.var
    syms = _extract_var_gene_symbols(src_var)

    X_csr = src_X.tocsr() if sp.issparse(src_X) else sp.csr_matrix(src_X)
    nz_sample = X_csr.data[: min(20000, len(X_csr.data))]
    if nz_sample.size > 0 and not np.allclose(nz_sample, np.round(nz_sample), atol=1e-3) and float(np.max(nz_sample)) <= 30.0:
        X_csr = sp.csr_matrix((np.expm1(np.clip(X_csr.data.astype(np.float64), 0.0, 25.0)), X_csr.indices, X_csr.indptr), shape=X_csr.shape)

    sids_raw = [str(v).strip() for v in adata.obs[sample_col].tolist()]
    u_s, cnt_s = np.unique(np.asarray(sids_raw, dtype=object), return_counts=True)
    kept = sorted(str(s) for s, c in zip(u_s.tolist(), cnt_s.tolist(), strict=True) if int(c) >= 25)
    s_map = {s: i for i, s in enumerate(kept)}
    grp = np.array([s_map.get(s, -1) for s in sids_raw], dtype=np.int64)
    valid_cells = np.flatnonzero(grp >= 0)
    ind = sp.csr_matrix(
        (np.ones(len(valid_cells), dtype=np.float64), (grp[valid_cells], valid_cells)),
        shape=(len(kept), adata.n_obs),
    )
    sums = np.asarray((ind @ X_csr).toarray(), dtype=np.float64)
    lib = np.maximum(np.sum(sums, axis=1, keepdims=True), 1.0)
    log2_cpm = np.log2((sums / lib) * 1e6 + 1.0)
    c_counts = [int(np.sum(grp == i)) for i in range(len(kept))]

    meta_df = pl.DataFrame(
        {
            "sample_id": kept,
            "n_cells": c_counts,
            "response_binary": [float("nan")] * len(kept),
            "biopsy_timepoint": ["Unknown"] * len(kept),
            "emp_frac_T_NK": [float("nan")] * len(kept),
            "emp_frac_Myeloid": [float("nan")] * len(kept),
            "emp_frac_CAF": [float("nan")] * len(kept),
            "emp_frac_Endothelial": [float("nan")] * len(kept),
            "emp_frac_Malignant": [float("nan")] * len(kept),
        }
    )
    return log2_cpm, kept, syms, meta_df


def _enrich_sc9_audited_response_metadata(
    cohort_id: str,
    meta_df: pl.DataFrame,
    pb_cache_dir: Path,
) -> pl.DataFrame:
    """Backfill audited SC9 response_binary & treatment_status from companion sample_meta parquet if needed."""
    valid_resp_count = meta_df.filter(pl.col("response_binary").is_in([0.0, 1.0])).height
    if valid_resp_count >= 5:
        return meta_df

    comp_path = pb_cache_dir / f"{cohort_id}.sample_meta.parquet"
    if not comp_path.is_file():
        return meta_df

    comp = pl.read_parquet(comp_path)
    if "sample_id" not in comp.columns or "response_binary" not in comp.columns:
        return meta_df

    lookup_r = {
        str(r["sample_id"]).strip(): float(r["response_binary"])
        for r in comp.iter_rows(named=True)
        if r["response_binary"] in (0, 1, 0.0, 1.0)
    }
    lookup_tp = {
        str(r["sample_id"]).strip(): str(r.get("treatment_status", "Unknown"))
        for r in comp.iter_rows(named=True)
    }
    new_r: list[float] = []
    new_tp: list[str] = []
    for r in meta_df.iter_rows(named=True):
        sid = str(r["sample_id"]).strip()
        cur_r = float(r["response_binary"]) if r["response_binary"] is not None else float("nan")
        new_r.append(lookup_r.get(sid, cur_r))
        cur_tp = str(r["biopsy_timepoint"])
        new_tp.append(lookup_tp.get(sid, cur_tp) if cur_tp == "Unknown" else cur_tp)

    return meta_df.with_columns(
        [
            pl.Series("response_binary", new_r),
            pl.Series("biopsy_timepoint", new_tp),
        ]
    )


def run_step03(
    out_dir: Path,
    preprocessed_dir: Path,
    lair_dir: Path,
    pb_cache_dir: Path,
    sc9_milo_dir: Path,
) -> None:
    """Execute Step 03 across all 23 single-cell cohorts (9 Milopy SC9 + 14 Reference Atlases)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    pb_fges_cache_dir = out_dir / "sc_pseudobulk_29fges_cache"
    pb_fges_cache_dir.mkdir(parents=True, exist_ok=True)

    tcga_ref = pl.read_parquet(out_dir / "tcga_bagaev_reference_29fges.parquet")
    collection = get_bagaev_mfp_collection()
    knn_cfg = BagaevKnnConfig(k_neighbors=35, metric="correlation")

    stromal_markers = {"COL1A1", "COL1A2", "DCN", "LUM", "FAP", "ACTA2", "PECAM1", "VWF", "CDH5", "KDR"}
    immune_markers = {"PTPRC", "CD3D", "CD3E", "CD8A", "CD4", "NKG7", "CD68", "CD14", "HLA-DRA", "MS4A1"}

    cohort_summary_rows: list[dict[str, str | int | float]] = []
    all_sample_tables: list[pl.DataFrame] = []
    emp_corr_rows: list[dict[str, str | float]] = []

    for spec in SC23_SPECS:
        cache_file = pb_fges_cache_dir / f"{spec.cohort_id}_raw_fges.parquet"
        if cache_file.is_file():
            raw_cohort_df = pl.read_parquet(cache_file)
            print(f"[Cached] Loaded 29-Fges pseudobulks for {spec.cohort_id:24s}: {raw_cohort_df.height} samples")
        else:
            print(f"[Streaming] Aggregating full-transcriptome pseudobulks for {spec.cohort_id} ({spec.display_name})...")
            if spec.source_kind == "preprocessed_h5ad":
                mat, sids, syms, meta_df = _stream_pseudobulk_from_h5ad(
                    preprocessed_dir / spec.rel_path, spec.sample_col, min_cells=spec.min_cells_per_sample
                )
            elif spec.source_kind == "lair_h5ad":
                mat, sids, syms, meta_df = _stream_pseudobulk_from_h5ad(
                    lair_dir / spec.rel_path, spec.sample_col, min_cells=spec.min_cells_per_sample
                )
            elif spec.source_kind == "lair_multi_h5ad":
                mat, sids, syms, meta_df = _aggregate_azizi_multi_h5ad(lair_dir / spec.rel_path)
            elif spec.source_kind == "lair_gz_h5ad":
                mat, sids, syms, meta_df = _aggregate_sharma_gz_h5ad(lair_dir / spec.rel_path, spec.sample_col)
            else:
                raise ValueError(f"Unknown source_kind: {spec.source_kind}")

            if spec.cohort_group == "SC9_Milopy":
                meta_df = _enrich_sc9_audited_response_metadata(spec.cohort_id, meta_df, pb_cache_dir)

            # Empirical pre-scaling stromal vs immune log2(CPM+1) diagnostic
            sym_upper = [g.strip().upper() for g in syms]
            st_idx = [i for i, g in enumerate(sym_upper) if g in stromal_markers]
            im_idx = [i for i, g in enumerate(sym_upper) if g in immune_markers]
            mean_stromal_expr = np.mean(mat[:, st_idx], axis=1) if st_idx else np.zeros(len(sids))
            mean_immune_expr = np.mean(mat[:, im_idx], axis=1) if im_idx else np.zeros(len(sids))

            ssg_res = compute_ssgsea_matrix(mat, sids, syms, collection)
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
                        pl.lit(len(set(sym_upper))).alias("n_genes_detected"),
                        pl.Series("raw_stromal_log2cpm", mean_stromal_expr.tolist()),
                        pl.Series("raw_immune_log2cpm", mean_immune_expr.tolist()),
                    ]
                )
                .hstack(ssg_res.unwrap().drop("sample_id"))
            )
            raw_cohort_df.write_parquet(cache_file)
            print(f"  -> Aggregated {raw_cohort_df.height} samples across {len(set(sym_upper))} unique genes.")

        if spec.cohort_group == "SC9_Milopy":
            raw_cohort_df = _enrich_sc9_audited_response_metadata(spec.cohort_id, raw_cohort_df, pb_cache_dir)

        # Within-cohort median/MAD scaling
        sc_res = median_mad_scale_signatures(raw_cohort_df, BAGAEV_MFP_SIGNATURE_ORDER, group_col=Nothing, clip=Some(4.0))
        if isinstance(sc_res, Failure):
            continue
        scaled_c = sc_res.unwrap().with_columns(
            (pl.col("cohort") + pl.lit("__") + pl.col("sample_id").cast(pl.String)).alias("global_sample_id")
        )

        # Supervised TCGA K=35 k-NN (indication-matched if available, else Pan-Cancer)
        tcga_projs = INDICATION_TO_TCGA_PROJECTS.get(spec.indication)
        ref_sub = tcga_ref.filter(pl.col("TCGA_project").is_in(list(tcga_projs))) if tcga_projs else tcga_ref
        if ref_sub.height < 35:
            ref_sub = tcga_ref

        knn_in = scaled_c.drop("sample_id").rename({"global_sample_id": "sample_id"})
        knn_res = classify_bagaev_knn(knn_in, ref_sub, config=knn_cfg)
        if isinstance(knn_res, Failure):
            continue
        knn_df = knn_res.unwrap().rename({"sample_id": "global_sample_id"})

        # De novo Louvain clustering if >= 8 samples
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

        # 29x29 correlation matrix conservation vs matched TCGA & Pan-Cancer TCGA
        cons_res = evaluate_fges_correlation_conservation(scaled_c, ref_sub)
        cons_metrics = cons_res.unwrap() if isinstance(cons_res, Success) else {}

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

        # Empirical single-cell fraction validation where cell-type annotations exist
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

    # Compare Whole-Tumor vs FACS-Immune-Sorted 29x29 TCGA conservation
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

    # Evaluate Milopy SC9 ICB response associations (Pooled All-9, Whole-Tumor 7, FACS 2, and Per-Cohort)
    sc9_samples = master_sc_df.filter(
        (pl.col("cohort_group") == "SC9_Milopy") & pl.col("response_binary").cast(pl.Float64).is_in([0.0, 1.0])
    )
    sc9_resp_tables: list[pl.DataFrame] = []

    strata_defs = [
        ("SC9_All_9_Pooled", sc9_samples),
        ("SC9_WholeTumor_7_Pooled", sc9_samples.filter(pl.col("sorting_regime") == "Whole_Tumor")),
        ("SC9_FACS_ImmuneSorted_2_Pooled", sc9_samples.filter(pl.col("sorting_regime") == "FACS_Immune_Sorted")),
    ]
    sc9_pre = sc9_samples.filter(pl.col("biopsy_timepoint").is_in(["Pre", "Unknown", "Combined"]))
    if sc9_pre.height >= 20:
        strata_defs.append(("SC9_Baseline_Pre_WholeTumor", sc9_pre.filter(pl.col("sorting_regime") == "Whole_Tumor")))

    for stratum_name, sub_df in strata_defs:
        for m_col, m_name in (("bagaev_subtype", "TCGA_kNN"), ("louvain_subtype", "Cohort_Louvain")):
            r_assoc = evaluate_subtype_response_associations(sub_df, subtype_col=m_col, cohort_col=Some("cohort"))
            if isinstance(r_assoc, Success):
                sc9_resp_tables.append(
                    r_assoc.unwrap().with_columns(
                        [pl.lit(stratum_name).alias("cohort"), pl.lit(m_name).alias("method")]
                    )
                )

    per_cohort_aucs: list[dict[str, str | float | int]] = []
    for spec in SC23_SPECS:
        if spec.cohort_group != "SC9_Milopy":
            continue
        c_sub = sc9_samples.filter(pl.col("cohort") == spec.cohort_id)
        if c_sub.height < 5:
            continue
        c_assoc = evaluate_subtype_response_associations(c_sub, subtype_col="bagaev_subtype", cohort_col=Nothing)
        if isinstance(c_assoc, Success):
            sc9_resp_tables.append(
                c_assoc.unwrap().with_columns(
                    [pl.lit(spec.cohort_id).alias("cohort"), pl.lit("TCGA_kNN").alias("method")]
                )
            )
        yc = c_sub["response_binary"].to_numpy().astype(np.int64)
        if int(np.sum(yc == 1)) >= 2 and int(np.sum(yc == 0)) >= 2:
            ie_score = (c_sub["prob_IE"] + 0.5 * c_sub["prob_IEF"]).to_numpy().astype(np.float64)
            imm_minus_fib = (c_sub["antitumor_immune_score"] - 0.5 * c_sub["stromal_fibrotic_score"]).to_numpy().astype(np.float64)
            per_cohort_aucs.append(
                {
                    "cohort": spec.cohort_id,
                    "display_name": spec.display_name,
                    "sorting_regime": spec.sorting_regime,
                    "n_evaluable_samples": int(c_sub.height),
                    "n_responders": int(np.sum(yc == 1)),
                    "auc_bagaev_knn_ie_prob": float(roc_auc_score(yc, ie_score)),
                    "auc_bagaev_immune_minus_fibrotic": float(roc_auc_score(yc, imm_minus_fib)),
                }
            )

    if sc9_resp_tables:
        sc9_resp_df = pl.concat(sc9_resp_tables, how="diagonal_relaxed")
        sc9_resp_df.write_parquet(out_dir / "sc9_milopy_subtype_response_rates.parquet")
        sc9_resp_df.write_csv(out_dir / "sc9_milopy_subtype_response_rates.tsv", separator="\t")

    if per_cohort_aucs:
        auc_df = pl.DataFrame(per_cohort_aucs)
        # Join Milopy single-cell DA benchmark summary if present for side-by-side comparison
        milo_da_path = sc9_milo_dir / "milopy_cell_state_da_all_resolutions.parquet"
        if milo_da_path.is_file():
            milo_df = pl.read_parquet(milo_da_path)
            coh_col = next((c for c in ("cohort", "cohort_id", "dataset_id") if c in milo_df.columns), None)
            lfc_col = next((c for c in ("logFC", "mean_logFC", "milo_logfc") if c in milo_df.columns), None)
            if coh_col and lfc_col:
                milo_agg = (
                    milo_df.group_by(coh_col)
                    .agg(pl.col(lfc_col).abs().max().alias("max_abs_milopy_state_logFC"))
                    .rename({coh_col: "cohort"})
                )
                auc_df = auc_df.join(milo_agg, on="cohort", how="left")
        auc_df.write_parquet(out_dir / "sc9_milopy_cohort_bagaev_aucs.parquet")
        auc_df.write_csv(out_dir / "sc9_milopy_cohort_bagaev_aucs.tsv", separator="\t")

    print(f"Step 03 complete. Saved 23-cohort single-cell pseudobulk Bagaev outputs to {out_dir}")


def main() -> None:
    paths = get_data_paths()
    default_lair = Path("/storage/halu/lair") if Path("/storage/halu/lair").exists() else paths.lair_dir
    local_pb = Path("data/pseudobulk")
    default_pb = local_pb if local_pb.is_dir() else paths.pseudobulk_dir
    local_milo = Path("output/sc9_milo_deconv_validation/04_milopy")
    default_milo = local_milo if local_milo.is_dir() else (paths.output_dir / "sc9_milo_deconv_validation/04_milopy")

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=paths.output_dir / "bagaev_mfp_analysis")
    parser.add_argument("--preprocessed-dir", type=Path, default=paths.preprocessed_dir)
    parser.add_argument("--lair-dir", type=Path, default=default_lair)
    parser.add_argument("--pb-cache-dir", type=Path, default=default_pb)
    parser.add_argument("--sc9-milo-dir", type=Path, default=default_milo)
    args = parser.parse_args()
    run_step03(
        out_dir=args.out_dir,
        preprocessed_dir=args.preprocessed_dir,
        lair_dir=args.lair_dir,
        pb_cache_dir=args.pb_cache_dir,
        sc9_milo_dir=args.sc9_milo_dir,
    )


if __name__ == "__main__":
    main()
