"""Memory-bounded tumor-only single-cell H5AD pseudobulk streaming and baseline deduplication."""

from __future__ import annotations

import gzip
from pathlib import Path
import re
import shutil
import tempfile
from typing import Any, Final, Literal, Mapping, assert_never, cast
import anndata as ad  # type: ignore[import-untyped]
from anndata.io import read_elem  # type: ignore[import-untyped]
import h5py  # type: ignore[import-untyped]
import numpy as np
import numpy.typing as npt
import pandas as pd
import polars as pl
from pydantic import BaseModel, ConfigDict
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success
import scipy.sparse as sp

from ..preprocessing.metadata import binarize_response, standardize_timepoint

ScaleMode = Literal["layer_counts", "raw_counts", "ln_cp10k", "log2_tpm"]
DeduplicateMode = Literal["strict_pre", "prefer_pre_one_per_patient"]

_NON_TUMOR_SAMPLE_RE: Final[re.Pattern[str]] = re.compile(
    r"(^|[_-])(pbmc|blood|normal|adjnorm)([_-]|$)", re.IGNORECASE
)
_GSE120575_PAT_RE: Final[re.Pattern[str]] = re.compile(r"^(Pre|Post)_([A-Za-z0-9]+?)(?:_\d+)?$", re.IGNORECASE)

_RESP_CANDIDATE_COLS: Final[tuple[str, ...]] = (
    "response_binary",
    "clinical_response",
    "response",
    "outcome",
    "Combined_outcome",
    "Path_response",
    "RECIST",
    "therapeutic_response",
)

_TIMEPOINT_CANDIDATE_COLS: Final[tuple[str, ...]] = (
    "biopsy_timepoint",
    "pre_post",
    "Resource",
    "Stage",
    "treatment_status",
    "timepoint",
    "time_point",
)

_CELLTYPE_CANDIDATE_COLS: Final[tuple[str, ...]] = (
    "clMidwayPr",
    "cell_type_coarse_crc_atlas",
    "author_cell_type",
    "cell_type_main",
    "Cell_type.refined",
    "Cell_type",
    "major_cell_type",
    "celltype",
    "MajorCluster",
    "cell_type",
    "cell_type_author",
    "SINGLECELL_TYPE",
)

_PATIENT_CANDIDATE_COLS: Final[tuple[str, ...]] = (
    "PMID_donor_id",
    "donor_id",
    "Patient_ID",
    "patient_id",
    "Patient",
    "PID",
    "Characteristics[individual]",
    "patientno",
    "patient",
)

_GONDAL_CANCER_MAP: Final[dict[str, str]] = {
    "TNBC": "Breast",
    "ER+": "Breast",
    "HER2+": "Breast",
    "ccRCC": "ccRCC",
    "HCC": "HCC",
    "iCCA": "CHOL",
    "Mel": "Melanoma",
    "MBM": "Melanoma",
    "BCC": "Melanoma",
}

_QIAN_DISEASE_MAP: Final[dict[str, str]] = {
    "breast cancer": "Breast",
    "colorectal cancer": "CRC",
    "ovarian cancer": "OV",
}


class TumorPseudobulkSpec(BaseModel):
    """Immutable specification for tumor-only single-cell pseudobulk aggregation."""

    model_config = ConfigDict(frozen=True)

    cohort_id: str
    sample_col: str
    default_indication: str = "Pan-Cancer"
    min_cells_per_sample: int = 25
    chunk_rows: int = 50000
    preferred_layer: str = "counts"


class TumorPseudobulkResult(BaseModel):
    """Immutable output container from streaming tumor-only H5AD pseudobulk aggregation."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    log2_cpm: npt.NDArray[np.float64]
    sample_ids: tuple[str, ...]
    gene_symbols: tuple[str, ...]
    meta_df: pl.DataFrame
    scale_mode: ScaleMode
    n_cells_pre_filter: int
    n_cells_post_tumor_filter: int


def load_ensembl_symbol_map(
    ensembl_dir: Path,
    fallback_h5ad: Maybe[Path] = Nothing,
) -> dict[str, str]:
    """Purely load Ensembl-ID-to-HGNC-symbol dictionary from Ensembl release cache or H5AD."""
    ens_map: dict[str, str] = {}
    cache_pq = ensembl_dir / "gene_mapping_cache_release_111.parquet"
    if cache_pq.is_file():
        df_map = pl.read_parquet(cache_pq).select(["gene_id", "gene_name"])
        for gid_raw, sym_raw in df_map.iter_rows():
            gid = str(gid_raw or "").split(".")[0].strip().upper()
            sym = str(sym_raw or "").strip().upper()
            if gid.startswith("ENSG") and sym and sym != "NAN" and not sym.startswith("ENSG"):
                ens_map.setdefault(gid, sym)

    match fallback_h5ad:
        case Some(fb_path) if fb_path.is_file():
            with h5py.File(fb_path, "r") as f:
                vdf = cast(pd.DataFrame, read_elem(f["var"]))
                col_sym = next(
                    (str(c) for c in ("gene_name", "original_id", "gene_symbol") if c in vdf.columns),
                    "",
                )
                if col_sym:
                    for eid, sym_val in zip(vdf.index.tolist(), vdf[col_sym].tolist(), strict=False):
                        clean_e = str(eid).split(".")[0].strip().upper()
                        clean_s = str(sym_val).strip().upper()
                        if clean_e.startswith("ENSG") and clean_s and clean_s != "NAN" and not clean_s.startswith("ENSG"):
                            ens_map.setdefault(clean_e, clean_s)
        case _:
            pass
    return ens_map


def extract_h5ad_var_symbols(
    var_df: pd.DataFrame,
    ens_map: Mapping[str, str],
) -> tuple[str, ...]:
    """Resolve upper-case Hugo gene symbols from an AnnData `.var` table (including bare ENSG indices)."""
    for col in ("gene_name", "gene_symbol", "gene_symbols", "feature_name", "symbol"):
        if col in var_df.columns:
            vals = [str(v).strip().upper() for v in var_df[col].tolist()]
            if any(v and not v.startswith("ENSG") and v != "NAN" for v in vals[:50]):
                return tuple(vals)

    idx_vals = [str(v).strip().upper() for v in var_df.index.tolist()]
    if "original_id" in var_df.columns and idx_vals and idx_vals[0].startswith("ENSG"):
        orig_vals = [str(v).strip().upper() for v in var_df["original_id"].tolist()]
        if any(v and not v.startswith("ENSG") and v != "NAN" for v in orig_vals[:50]):
            return tuple(orig_vals)

    if (
        idx_vals
        and ens_map
        and sum(1 for v in idx_vals[:50] if v.startswith("ENSG")) >= max(1, min(20, len(idx_vals[:50]) // 2))
    ):
        return tuple(ens_map.get(v.split(".")[0], v) for v in idx_vals)
    return tuple(idx_vals)


def classify_broad_lineage(ct_raw: str) -> str:
    """Classify heterogeneous single-cell author labels into broad TME compartments."""
    s = ct_raw.strip().lower()
    if any(k in s for k in ("fibro", "caf", "strom", "mesench", "pericyte", "peri", "smooth muscle", "smoothmuscle", "stellate")):
        return "Stromal_CAF"
    if any(k in s for k in ("endothel", "vascular", "lymphatic")) or s == "endo":
        return "Endothelial"
    if any(k in s for k in ("cd8", "cd4", "t cell", "t_cell", "tcell", "treg", "nk", "lympho", "cytotox", "tnkilc")) or s in {"t", "t/nk"}:
        return "T_NK_Lymphoid"
    if any(k in s for k in ("b cell", "b_cell", "bcell", "plasma")) or s == "b":
        return "B_Plasma"
    if any(k in s for k in ("macro", "mono", "myeloid", "dc", "dendrit", "neutro", "mast", "tam", "granulo")) or s == "mye":
        return "Myeloid"
    if any(k in s for k in ("malig", "tumor", "cancer", "epithel", "melanom", "hepatoc", "carcinoma", "neoplast")) or s in {"epi"}:
        return "Malignant_Epithelial"
    return "Other"


def build_tumor_cell_mask(
    obs: pd.DataFrame,
    cohort_id: str,
    sample_values: Maybe[tuple[str, ...]] = Nothing,
) -> npt.NDArray[np.bool_]:
    """Build strict tumor-only cell mask, excluding PBMC/blood, normal adjacent, benign, and uninvolved LN."""
    n = len(obs)
    mask = np.ones(n, dtype=bool)

    if cohort_id == "CELLxGENE_05a8c945" and "sample_type" in obs.columns:
        mask &= obs["sample_type"].astype(str).str.strip().str.lower().isin({"tumor", "metastasis"}).to_numpy()
    elif cohort_id == "Cheng_2021_PanMyeloid" and "tissue" in obs.columns:
        mask &= (obs["tissue"].astype(str).str.strip().str.upper() == "T").to_numpy()
    elif cohort_id == "Khaliq_2022_CRC" and "Condition" in obs.columns:
        mask &= (obs["Condition"].astype(str).str.strip().str.lower() == "tumor").to_numpy()
    elif cohort_id == "Kim_2020_LUAD" and "Sample_Origin" in obs.columns:
        mask &= obs["Sample_Origin"].astype(str).str.strip().isin({"tLung", "tL/B", "mBrain", "mLN"}).to_numpy()
    elif cohort_id == "Lu_2022_HCC" and "site" in obs.columns:
        mask &= obs["site"].astype(str).str.strip().isin({"Tumor", "PVTT"}).to_numpy()
    elif cohort_id == "Pelka_2021_CRC" and "SPECIMEN_TYPE" in obs.columns:
        mask &= (obs["SPECIMEN_TYPE"].astype(str).str.strip().str.upper() == "T").to_numpy()
    elif cohort_id == "Pu_2021_THCA" and "biopsy_site" in obs.columns:
        mask &= obs["biopsy_site"].astype(str).str.strip().isin({"T", "SC"}).to_numpy()
    elif cohort_id == "Qian_2020_PanCancer" and "Characteristics[sampling site]" in obs.columns:
        mask &= (
            obs["Characteristics[sampling site]"]
            .astype(str)
            .str.strip()
            .str.lower()
            .isin({"neoplasm", "neoplasm core", "neoplasm border"})
            .to_numpy()
        )
    elif cohort_id == "Sharma_2020_HCC" and "NormalvsTumor" in obs.columns:
        mask &= (obs["NormalvsTumor"].astype(str).str.strip() == "Tumor").to_numpy()
    elif cohort_id == "Borcherding_2021_ccRCC" and "sample" in obs.columns:
        mask &= obs["sample"].astype(str).str.strip().str.endswith("_T").to_numpy()

    if cohort_id == "Biermann_2022_BrainMet" and "cell_type_main" in obs.columns:
        mask &= (~obs["cell_type_main"].astype(str).isin({"Low-quality cells", "Contamination", "Doublets"})).to_numpy()

    match sample_values:
        case Some(svals) if len(svals) == n:
            sid_ok = np.array([not bool(_NON_TUMOR_SAMPLE_RE.search(s)) for s in svals], dtype=bool)
            mask &= sid_ok
        case _:
            pass

    return mask


def _detect_scale_mode(
    x_node: Any,
    has_preferred_layer: bool,
    is_csr_or_csc: bool,
    n_cells: int,
) -> ScaleMode:
    """Detect whether source matrix stores raw counts, ln(CP10k+1), or log2(TPM+1)."""
    if has_preferred_layer:
        return "layer_counts"
    sample_data = (
        np.asarray(x_node["data"][: min(50000, len(x_node["data"]))], dtype=np.float64)
        if is_csr_or_csc and len(x_node["data"]) > 0
        else np.asarray(x_node[: min(200, n_cells)], dtype=np.float64).ravel()
    )
    nz = sample_data[np.isfinite(sample_data) & (sample_data != 0.0)]
    if nz.size == 0:
        return "raw_counts"
    if bool(np.allclose(nz[:5000], np.round(nz[:5000]), atol=1e-3)):
        return "raw_counts"
    max_v = float(np.nanmax(nz))
    return "ln_cp10k" if max_v <= 12.5 else "log2_tpm"


def _invert_to_linear(arr: npt.NDArray[np.float64], mode: ScaleMode) -> npt.NDArray[np.float64]:
    """Convert expression values back to linear scale strictly matching the detected log base."""
    match mode:
        case "layer_counts" | "raw_counts":
            return arr
        case "ln_cp10k":
            return np.asarray(np.expm1(np.clip(arr, 0.0, 20.0)), dtype=np.float64)
        case "log2_tpm":
            return np.asarray(np.exp2(np.clip(arr, 0.0, 26.0)) - 1.0, dtype=np.float64)
        case _ as unreachable:
            assert_never(unreachable)


def _min_cell_total(mode: ScaleMode) -> float:
    """Minimum total cellular transcript mass per cell by scale mode."""
    match mode:
        case "layer_counts" | "raw_counts":
            return 150.0
        case "ln_cp10k" | "log2_tpm":
            return 1.0
        case _ as unreachable:
            assert_never(unreachable)


def _resolve_sample_column(obs: pd.DataFrame, preferred_col: str) -> str:
    """Choose sample identifier column from `.obs` deterministically."""
    if preferred_col in obs.columns:
        return preferred_col
    return next(
        (
            str(c)
            for c in ("sample_id", "sample", "samples", "patient_id", "patient", "donor_id", "orig.ident")
            if c in obs.columns
        ),
        str(obs.columns[0]),
    )


def _extract_patient_id(obs_sub: pd.DataFrame, cohort_id: str, sid: str) -> str:
    """Extract normalized patient/donor identifier for a single-cell sample."""
    if cohort_id == "GSE120575":
        m = _GSE120575_PAT_RE.match(sid)
        return m.group(2) if m is not None else sid
    if cohort_id == "Kim_2020_LUAD":
        # Note Kim_2020 obs['patient'] stores Sample_Origin ('tLung', 'mLN'); obs['Sample'] identifies specimen
        return str(obs_sub["Sample"].iloc[0]).strip() if "Sample" in obs_sub.columns else sid
    if cohort_id == "GSE200996" and "Patient_ID" in obs_sub.columns:
        vals = [str(v).strip() for v in obs_sub["Patient_ID"].dropna().unique().tolist() if str(v).strip() not in ("", "NA", "nan", "None")]
        return vals[0] if vals else "NA"

    pat_col = next((c for c in _PATIENT_CANDIDATE_COLS if c in obs_sub.columns), "")
    if pat_col:
        vals = [
            str(v).strip()
            for v in obs_sub[pat_col].dropna().unique().tolist()
            if str(v).strip() not in ("", "NA", "nan", "None", "Unknown", "unknown")
        ]
        if vals:
            return vals[0]
    return sid


def _extract_timepoint(obs_sub: pd.DataFrame, cohort_id: str, sid: str) -> str:
    """Extract harmonized biopsy timepoint {'Pre', 'On-Treatment', 'Post', 'Unknown'}."""
    if cohort_id == "GSE120575":
        s_low = sid.strip().lower()
        if s_low.startswith("pre_"):
            return "Pre"
        if s_low.startswith("post_"):
            return "Post"
    tp_col = next((c for c in _TIMEPOINT_CANDIDATE_COLS if c in obs_sub.columns), "")
    if tp_col and len(obs_sub) > 0:
        raw_vals = [str(v).strip() for v in obs_sub[tp_col].dropna().unique().tolist()]
        for raw in raw_vals:
            low = raw.lower()
            if low in ("naive", "pre", "pre-tx", "pre-treatment", "pre-treatment biopsy"):
                return "Pre"
            if low in ("treated", "post-treatment surgery"):
                return "Post"
            std = standardize_timepoint(raw)
            if std != "Unknown":
                return str(std)
    return "Unknown"


def _extract_indication(obs_sub: pd.DataFrame, cohort_id: str, default_indication: str) -> str:
    """Extract per-sample cancer indication for pan-cancer single-cell cohorts."""
    if cohort_id == "CELLxGENE_7b20c613" and "Cancer_type_update" in obs_sub.columns and len(obs_sub) > 0:
        raw_ct = str(obs_sub["Cancer_type_update"].iloc[0]).strip()
        return _GONDAL_CANCER_MAP.get(raw_ct, default_indication)
    if cohort_id == "Qian_2020_PanCancer" and "Characteristics[disease]" in obs_sub.columns and len(obs_sub) > 0:
        raw_dis = str(obs_sub["Characteristics[disease]"].iloc[0]).strip().lower()
        return _QIAN_DISEASE_MAP.get(raw_dis, default_indication)
    return default_indication


def _build_sample_metadata(
    obs: pd.DataFrame,
    raw_sids: tuple[str, ...],
    tumor_mask: npt.NDArray[np.bool_],
    final_samples: tuple[str, ...],
    kept_cell_counts: npt.NDArray[np.int64],
    spec: TumorPseudobulkSpec,
) -> pl.DataFrame:
    """Construct pure sample-level metadata (patient_id, response, timepoint, indication, cell fractions)."""
    resp_col = next((c for c in _RESP_CANDIDATE_COLS if c in obs.columns), "")
    ct_col = next((c for c in _CELLTYPE_CANDIDATE_COLS if c in obs.columns), "")

    obs_sid_arr = np.asarray(raw_sids, dtype=object)
    meta_rows: list[dict[str, str | int | float]] = []

    for s_idx, sid in enumerate(final_samples):
        row_mask = tumor_mask & (obs_sid_arr == sid)
        sub_obs = obs.loc[row_mask]

        r_bin = float("nan")
        if resp_col and len(sub_obs) > 0:
            vals_r = [binarize_response(v) for v in sub_obs[resp_col].dropna().unique().tolist()]
            finite_r = [float(v) for v in vals_r if not np.isnan(v)]
            if finite_r:
                r_bin = float(np.round(np.mean(finite_r)))

        pat_id = _extract_patient_id(sub_obs, spec.cohort_id, sid)
        tp_str = _extract_timepoint(sub_obs, spec.cohort_id, sid)
        ind_str = _extract_indication(sub_obs, spec.cohort_id, spec.default_indication)

        frac_t = float("nan")
        frac_mye = float("nan")
        frac_caf = float("nan")
        frac_endo = float("nan")
        frac_mal = float("nan")
        if ct_col and len(sub_obs) > 0:
            cats = [classify_broad_lineage(str(v)) for v in sub_obs[ct_col].tolist()]
            # Skip unannotated cohorts where all labels collapse to 'Other'
            if any(c != "Other" for c in cats):
                n_tot_c = float(len(cats))
                frac_t = sum(1 for c in cats if c == "T_NK_Lymphoid") / n_tot_c
                frac_mye = sum(1 for c in cats if c == "Myeloid") / n_tot_c
                frac_caf = sum(1 for c in cats if c == "Stromal_CAF") / n_tot_c
                frac_endo = sum(1 for c in cats if c == "Endothelial") / n_tot_c
                frac_mal = sum(1 for c in cats if c == "Malignant_Epithelial") / n_tot_c

        meta_rows.append(
            {
                "sample_id": sid,
                "patient_id": pat_id,
                "sample_indication": ind_str,
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

    return pl.DataFrame(meta_rows)


def stream_tumor_pseudobulk_h5ad(
    h5_path: Path,
    spec: TumorPseudobulkSpec,
    ens_map: Mapping[str, str],
) -> Result[TumorPseudobulkResult, str]:
    """Stream an H5AD in bounded chunks over tumor-only cells to build log2(CPM_1e6 + 1) pseudobulks."""
    if not h5_path.is_file():
        return Failure(f"H5AD file not found at {h5_path}.")

    with h5py.File(h5_path, "r") as f:
        obs = cast(pd.DataFrame, read_elem(f["obs"]))
        var = cast(pd.DataFrame, read_elem(f["var"]))
        gene_symbols = extract_h5ad_var_symbols(var, ens_map)
        n_cells = len(obs)
        n_genes = len(gene_symbols)

        actual_scol = _resolve_sample_column(obs, spec.sample_col)
        raw_sids = tuple(str(v).strip() for v in obs[actual_scol].tolist())
        tumor_mask = build_tumor_cell_mask(obs, spec.cohort_id, sample_values=Some(raw_sids))
        valid_sid_mask = np.array(
            [bool(s) and s.lower() not in ("nan", "none", "unknown", "na") for s in raw_sids],
            dtype=bool,
        )
        eligible_mask = tumor_mask & valid_sid_mask
        n_post_tumor = int(np.sum(eligible_mask))
        if n_post_tumor < spec.min_cells_per_sample:
            return Failure(
                f"Fewer than {spec.min_cells_per_sample} eligible tumor cells in {h5_path.name} ({n_post_tumor})."
            )

        has_pref_layer = "layers" in f and spec.preferred_layer in f["layers"]
        x_node: Any = f["layers"][spec.preferred_layer] if has_pref_layer else f["X"]
        is_sparse_grp = isinstance(x_node, h5py.Group) and "indptr" in x_node
        scale_mode = _detect_scale_mode(x_node, has_pref_layer, is_sparse_grp, n_cells)
        min_umi = _min_cell_total(scale_mode)

        s_arr = np.asarray(raw_sids, dtype=object)
        u_sids, u_cnts = np.unique(s_arr[eligible_mask], return_counts=True)
        kept_samples = tuple(
            sorted(
                str(s)
                for s, cnt in zip(u_sids.tolist(), u_cnts.tolist(), strict=True)
                if int(cnt) >= spec.min_cells_per_sample
            )
        )
        if not kept_samples:
            return Failure(
                f"No tumor samples with >= {spec.min_cells_per_sample} cells in {h5_path.name} (col={actual_scol})."
            )

        sid_to_idx = {s: i for i, s in enumerate(kept_samples)}
        cell_sample_idx = np.where(
            eligible_mask,
            np.array([sid_to_idx.get(s, -1) for s in raw_sids], dtype=np.int64),
            -1,
        )

        pb_sums = np.zeros((len(kept_samples), n_genes), dtype=np.float64)
        kept_cell_counts = np.zeros(len(kept_samples), dtype=np.int64)

        if is_sparse_grp:
            indptr = np.asarray(x_node["indptr"][:], dtype=np.int64)
            is_csc = len(indptr) == (n_genes + 1) and len(indptr) != (n_cells + 1)
            if is_csc:
                mat_full = read_elem(x_node).tocsr()
                for start in range(0, n_cells, spec.chunk_rows):
                    end = min(start + spec.chunk_rows, n_cells)
                    c_sidx = cell_sample_idx[start:end]
                    if not np.any(c_sidx >= 0):
                        continue
                    sub_csr = mat_full[start:end].tocsr()
                    if scale_mode in ("ln_cp10k", "log2_tpm") and sub_csr.nnz > 0:
                        lin_data = _invert_to_linear(np.asarray(sub_csr.data, dtype=np.float64), scale_mode)
                        sub_csr = sp.csr_matrix((lin_data, sub_csr.indices, sub_csr.indptr), shape=sub_csr.shape)
                    row_sums = np.asarray(sub_csr.sum(axis=1)).ravel()
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
                for start in range(0, n_cells, spec.chunk_rows):
                    end = min(start + spec.chunk_rows, n_cells)
                    c_sidx = cell_sample_idx[start:end]
                    if not np.any(c_sidx >= 0):
                        continue
                    p0 = int(indptr[start])
                    p1 = int(indptr[end])
                    if p1 <= p0:
                        continue
                    sub_indptr = indptr[start : end + 1] - p0
                    sub_indices = np.asarray(indices_ds[p0:p1], dtype=np.int32)
                    sub_data = _invert_to_linear(np.asarray(data_ds[p0:p1], dtype=np.float64), scale_mode)
                    sub_csr = sp.csr_matrix((sub_data, sub_indices, sub_indptr), shape=(end - start, n_genes))
                    row_sums = np.asarray(sub_csr.sum(axis=1)).ravel()
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
            for start in range(0, n_cells, spec.chunk_rows):
                end = min(start + spec.chunk_rows, n_cells)
                c_sidx = cell_sample_idx[start:end]
                if not np.any(c_sidx >= 0):
                    continue
                block = _invert_to_linear(np.asarray(x_node[start:end], dtype=np.float64), scale_mode)
                row_sums = np.nansum(block, axis=1)
                ok = (c_sidx >= 0) & (row_sums >= min_umi)
                for local_r in np.flatnonzero(ok):
                    g_id = int(c_sidx[local_r])
                    pb_sums[g_id] += np.nan_to_num(block[local_r], nan=0.0)
                    kept_cell_counts[g_id] += 1

    non_empty = kept_cell_counts >= max(10, spec.min_cells_per_sample // 2)
    if not np.any(non_empty):
        return Failure(f"All candidate tumor samples in {h5_path.name} failed UMI cell filtering.")

    pb_sums = pb_sums[non_empty]
    kept_cell_counts = kept_cell_counts[non_empty]
    final_samples = tuple(s for s, keep in zip(kept_samples, non_empty.tolist(), strict=True) if keep)

    lib_sizes = np.where((raw_ls := np.sum(pb_sums, axis=1, keepdims=True)) > 0.0, raw_ls, 1.0)
    log2_cpm = np.log2((pb_sums / lib_sizes) * 1e6 + 1.0)
    meta_df = _build_sample_metadata(
        obs=obs,
        raw_sids=raw_sids,
        tumor_mask=eligible_mask,
        final_samples=final_samples,
        kept_cell_counts=kept_cell_counts,
        spec=spec,
    )
    return Success(
        TumorPseudobulkResult(
            log2_cpm=log2_cpm,
            sample_ids=final_samples,
            gene_symbols=gene_symbols,
            meta_df=meta_df,
            scale_mode=scale_mode,
            n_cells_pre_filter=n_cells,
            n_cells_post_tumor_filter=n_post_tumor,
        )
    )


def aggregate_multi_file_tumor_pseudobulk(
    directory: Path,
    spec: TumorPseudobulkSpec,
    ens_map: Mapping[str, str],
) -> Result[TumorPseudobulkResult, str]:
    """Aggregate directory of single-sample tumor H5AD files (e.g. Azizi 2018 BRCA *TUMOR*.h5ad)."""
    if not directory.is_dir():
        return Failure(f"Multi-H5AD directory not found at {directory}.")

    tumor_files = sorted(p for p in directory.glob("*TUMOR*.h5ad") if not p.name.startswith("._"))
    if not tumor_files:
        tumor_files = sorted(
            p
            for p in directory.glob("*.h5ad")
            if not p.name.startswith("._") and not bool(_NON_TUMOR_SAMPLE_RE.search(p.stem))
        )
    if not tumor_files:
        return Failure(f"No tumor H5AD files found in {directory}.")

    sample_dicts: list[dict[str, float]] = []
    sample_names: list[str] = []
    cell_counts: list[int] = []
    all_genes_set: set[str] = set()
    total_pre = 0
    total_post = 0

    for fpath in tumor_files:
        a = ad.read_h5ad(fpath)
        total_pre += int(a.n_obs)
        if int(a.n_obs) < spec.min_cells_per_sample:
            continue
        syms = extract_h5ad_var_symbols(a.var, ens_map)
        x_raw = np.asarray(cast(Any, a.X).toarray() if sp.issparse(a.X) else a.X, dtype=np.float64)
        col_sums = np.sum(x_raw, axis=0).ravel()
        tot = float(np.sum(col_sums))
        if tot <= 0.0:
            continue
        total_post += int(a.n_obs)
        cpm_vec = (col_sums / tot) * 1e6
        s_dict: dict[str, float] = {}
        for g, val in zip(syms, cpm_vec.tolist(), strict=True):
            if g and g != "NAN" and val > 0.0:
                s_dict[g] = max(s_dict.get(g, 0.0), float(val))
        all_genes_set.update(s_dict.keys())
        sample_dicts.append(s_dict)
        sample_names.append(fpath.stem.replace("_counts", "").replace("_matrix", ""))
        cell_counts.append(int(a.n_obs))

    if not sample_names:
        return Failure(f"No non-empty tumor sample H5ADs in {directory}.")

    genes_sorted = tuple(sorted(all_genes_set))
    g_idx = {g: i for i, g in enumerate(genes_sorted)}
    mat = np.zeros((len(sample_names), len(genes_sorted)), dtype=np.float64)
    for r_i, s_dict in enumerate(sample_dicts):
        for g, v in s_dict.items():
            mat[r_i, g_idx[g]] = np.log2(v + 1.0)

    patient_ids = [re.sub(r"_(TUMOR\d*|T\d+)$", "", s, flags=re.IGNORECASE) for s in sample_names]
    meta_df = pl.DataFrame(
        {
            "sample_id": sample_names,
            "patient_id": patient_ids,
            "sample_indication": [spec.default_indication] * len(sample_names),
            "n_cells": cell_counts,
            "response_binary": [float("nan")] * len(sample_names),
            "biopsy_timepoint": ["Pre"] * len(sample_names),
            "emp_frac_T_NK": [float("nan")] * len(sample_names),
            "emp_frac_Myeloid": [float("nan")] * len(sample_names),
            "emp_frac_CAF": [float("nan")] * len(sample_names),
            "emp_frac_Endothelial": [float("nan")] * len(sample_names),
            "emp_frac_Malignant": [float("nan")] * len(sample_names),
        }
    )
    return Success(
        TumorPseudobulkResult(
            log2_cpm=mat,
            sample_ids=tuple(sample_names),
            gene_symbols=genes_sorted,
            meta_df=meta_df,
            scale_mode="raw_counts",
            n_cells_pre_filter=total_pre,
            n_cells_post_tumor_filter=total_post,
        )
    )


def _aggregate_in_memory_adata(
    adata: ad.AnnData,
    spec: TumorPseudobulkSpec,
    ens_map: Mapping[str, str],
) -> Result[TumorPseudobulkResult, str]:
    """Aggregate an in-memory AnnData object into tumor-only log2(CPM+1) pseudobulk profiles."""
    n_pre = int(adata.n_obs)
    if n_pre == 0:
        return Failure(f"Empty AnnData for {spec.cohort_id}.")
    scol = _resolve_sample_column(adata.obs, spec.sample_col)
    raw_sids_all = tuple(str(v).strip() for v in adata.obs[scol].tolist())
    tumor_mask = build_tumor_cell_mask(adata.obs, spec.cohort_id, sample_values=Some(raw_sids_all))
    if int(np.sum(tumor_mask)) < spec.min_cells_per_sample:
        return Failure(f"Fewer than {spec.min_cells_per_sample} tumor cells in {spec.cohort_id}.")

    sub_ad = adata[tumor_mask].copy()
    n_post = int(sub_ad.n_obs)
    use_raw = sub_ad.raw is not None and int(sub_ad.raw.n_vars) > int(sub_ad.n_vars)
    src_x = sub_ad.raw.X if use_raw and sub_ad.raw is not None else sub_ad.X
    src_var = sub_ad.raw.var if use_raw and sub_ad.raw is not None else sub_ad.var
    syms = extract_h5ad_var_symbols(src_var, ens_map)

    x_csr = cast(Any, src_x).tocsr() if sp.issparse(src_x) else sp.csr_matrix(src_x)
    nz = np.asarray(x_csr.data[: min(20000, len(x_csr.data))], dtype=np.float64)
    scale_mode: ScaleMode = "raw_counts"
    if nz.size > 0 and not np.allclose(nz, np.round(nz), atol=1e-3):
        scale_mode = "ln_cp10k" if float(np.nanmax(nz)) <= 12.5 else "log2_tpm"
        lin_data = _invert_to_linear(np.asarray(x_csr.data, dtype=np.float64), scale_mode)
        x_csr = sp.csr_matrix((lin_data, x_csr.indices, x_csr.indptr), shape=x_csr.shape)

    sids_raw = tuple(str(v).strip() for v in sub_ad.obs[scol].tolist())
    u_s, cnt_s = np.unique(np.asarray(sids_raw, dtype=object), return_counts=True)
    kept = tuple(
        sorted(
            str(s)
            for s, c in zip(u_s.tolist(), cnt_s.tolist(), strict=True)
            if int(c) >= spec.min_cells_per_sample and not bool(_NON_TUMOR_SAMPLE_RE.search(str(s)))
        )
    )
    if not kept:
        return Failure(f"No tumor samples with >= {spec.min_cells_per_sample} cells in {spec.cohort_id}.")

    s_map = {s: i for i, s in enumerate(kept)}
    grp = np.array([s_map.get(s, -1) for s in sids_raw], dtype=np.int64)
    valid_cells = np.flatnonzero(grp >= 0)
    ind = sp.csr_matrix(
        (np.ones(len(valid_cells), dtype=np.float64), (grp[valid_cells], valid_cells)),
        shape=(len(kept), n_post),
    )
    sums = np.asarray((ind @ x_csr).toarray(), dtype=np.float64)
    lib = np.maximum(np.sum(sums, axis=1, keepdims=True), 1.0)
    log2_cpm = np.log2((sums / lib) * 1e6 + 1.0)
    c_counts = np.array([int(np.sum(grp == i)) for i in range(len(kept))], dtype=np.int64)

    meta_df = _build_sample_metadata(
        obs=sub_ad.obs,
        raw_sids=sids_raw,
        tumor_mask=np.ones(n_post, dtype=bool),
        final_samples=kept,
        kept_cell_counts=c_counts,
        spec=spec,
    )
    return Success(
        TumorPseudobulkResult(
            log2_cpm=log2_cpm,
            sample_ids=kept,
            gene_symbols=syms,
            meta_df=meta_df,
            scale_mode=scale_mode,
            n_cells_pre_filter=n_pre,
            n_cells_post_tumor_filter=n_post,
        )
    )


def _is_complete_h5ad(h5_path: Path) -> bool:
    if not h5_path.is_file() or h5_path.stat().st_size == 0:
        return False
    try:
        with h5py.File(h5_path, "r") as f:
            return "X" in f and "obs" in f and "var" in f and len(f["obs"].keys()) > 0
    except Exception:
        return False


def stream_cohort_tumor_pseudobulk(
    dataset_id: str,
    spec: TumorPseudobulkSpec,
    ens_map: Mapping[str, str],
    *,
    auto_download: bool = False,
) -> Result[TumorPseudobulkResult, str]:
    """Stream tumor-only pseudobulk strictly from `TME_DATA_DIR/{preprocessed,raw}`."""
    from ..paths import get_data_paths
    from ..query import load_dataset

    paths = get_data_paths()
    h5_path = paths.preprocessed_dir / f"{dataset_id}.h5ad"
    if _is_complete_h5ad(h5_path):
        res = stream_tumor_pseudobulk_h5ad(h5_path, spec, ens_map)
        if isinstance(res, Success):
            return res

    raw_ds_dir = paths.raw_dir / dataset_id
    if dataset_id == "GSE156625":
        gz_h5 = raw_ds_dir / "GSE156625_HCCFscanpyobj.h5ad.gz"
        if gz_h5.is_file():
            return aggregate_gzip_h5ad_tumor_pseudobulk(gz_h5, spec, ens_map)

    ds_res = load_dataset(dataset_id, auto_download=auto_download)
    if isinstance(ds_res, Failure):
        return Failure(f"Failed loading {dataset_id} from {paths.data_root}: {ds_res.failure()}")
    return _aggregate_in_memory_adata(ds_res.unwrap(), spec, ens_map)


def aggregate_gzip_h5ad_tumor_pseudobulk(
    gz_path: Path,
    spec: TumorPseudobulkSpec,
    ens_map: Mapping[str, str],
) -> Result[TumorPseudobulkResult, str]:
    """Decompress gzipped H5AD (e.g. Sharma 2020 HCC) and aggregate tumor-only pseudobulk."""
    if not gz_path.is_file():
        return Failure(f"Gzipped H5AD file not found at {gz_path}.")

    with tempfile.NamedTemporaryFile(suffix=".h5ad") as tmp:
        with gzip.open(gz_path, "rb") as gi:
            shutil.copyfileobj(gi, tmp)
        tmp.flush()
        adata = ad.read_h5ad(tmp.name)

    n_pre = int(adata.n_obs)
    tumor_mask = build_tumor_cell_mask(adata.obs, spec.cohort_id)
    if int(np.sum(tumor_mask)) < spec.min_cells_per_sample:
        return Failure(f"Fewer than {spec.min_cells_per_sample} tumor cells in {gz_path.name}.")

    sub_ad = adata[tumor_mask].copy()
    n_post = int(sub_ad.n_obs)
    use_raw = sub_ad.raw is not None and int(sub_ad.raw.n_vars) > int(sub_ad.n_vars)
    src_x = sub_ad.raw.X if use_raw and sub_ad.raw is not None else sub_ad.X
    src_var = sub_ad.raw.var if use_raw and sub_ad.raw is not None else sub_ad.var
    syms = extract_h5ad_var_symbols(src_var, ens_map)

    x_csr = cast(Any, src_x).tocsr() if sp.issparse(src_x) else sp.csr_matrix(src_x)
    nz = np.asarray(x_csr.data[: min(20000, len(x_csr.data))], dtype=np.float64)
    scale_mode: ScaleMode = "raw_counts"
    if nz.size > 0 and not np.allclose(nz, np.round(nz), atol=1e-3):
        scale_mode = "ln_cp10k" if float(np.nanmax(nz)) <= 12.5 else "log2_tpm"
        lin_data = _invert_to_linear(np.asarray(x_csr.data, dtype=np.float64), scale_mode)
        x_csr = sp.csr_matrix((lin_data, x_csr.indices, x_csr.indptr), shape=x_csr.shape)

    scol = _resolve_sample_column(sub_ad.obs, spec.sample_col)
    sids_raw = tuple(str(v).strip() for v in sub_ad.obs[scol].tolist())
    u_s, cnt_s = np.unique(np.asarray(sids_raw, dtype=object), return_counts=True)
    kept = tuple(
        sorted(
            str(s)
            for s, c in zip(u_s.tolist(), cnt_s.tolist(), strict=True)
            if int(c) >= spec.min_cells_per_sample and not bool(_NON_TUMOR_SAMPLE_RE.search(str(s)))
        )
    )
    if not kept:
        return Failure(f"No tumor sections with >= {spec.min_cells_per_sample} cells in {gz_path.name}.")

    s_map = {s: i for i, s in enumerate(kept)}
    grp = np.array([s_map.get(s, -1) for s in sids_raw], dtype=np.int64)
    valid_cells = np.flatnonzero(grp >= 0)
    ind = sp.csr_matrix(
        (np.ones(len(valid_cells), dtype=np.float64), (grp[valid_cells], valid_cells)),
        shape=(len(kept), n_post),
    )
    sums = np.asarray((ind @ x_csr).toarray(), dtype=np.float64)
    lib = np.maximum(np.sum(sums, axis=1, keepdims=True), 1.0)
    log2_cpm = np.log2((sums / lib) * 1e6 + 1.0)
    c_counts = np.array([int(np.sum(grp == i)) for i in range(len(kept))], dtype=np.int64)

    meta_df = _build_sample_metadata(
        obs=sub_ad.obs,
        raw_sids=sids_raw,
        tumor_mask=np.ones(n_post, dtype=bool),
        final_samples=kept,
        kept_cell_counts=c_counts,
        spec=spec,
    )
    return Success(
        TumorPseudobulkResult(
            log2_cpm=log2_cpm,
            sample_ids=kept,
            gene_symbols=syms,
            meta_df=meta_df,
            scale_mode=scale_mode,
            n_cells_pre_filter=n_pre,
            n_cells_post_tumor_filter=n_post,
        )
    )


def deduplicate_patient_profiles(
    df: pl.DataFrame,
    mode: DeduplicateMode = "strict_pre",
    patient_col: str = "patient_id",
    cohort_col: str = "cohort",
    timepoint_col: str = "biopsy_timepoint",
    weight_col: Maybe[str] = Nothing,
) -> pl.DataFrame:
    """Select one independent tumor profile per `(cohort, patient_id)` without longitudinal leakage."""
    if df.is_empty():
        return df

    pat_expr = (
        pl.col(patient_col).cast(pl.String).str.strip_chars()
        if patient_col in df.columns
        else pl.col("sample_id").cast(pl.String)
    )
    coh_expr = (
        pl.col(cohort_col).cast(pl.String)
        if cohort_col in df.columns
        else pl.lit("Cohort")
    )
    tp_expr = (
        pl.col(timepoint_col).cast(pl.String).str.strip_chars().str.to_lowercase()
        if timepoint_col in df.columns
        else pl.lit("unknown")
    )
    id_col = "global_sample_id" if "global_sample_id" in df.columns else "sample_id"

    work = df.with_columns(
        [
            pat_expr.alias("_dedup_pat"),
            coh_expr.alias("_dedup_coh"),
            tp_expr.alias("_dedup_tp"),
        ]
    ).filter(
        pl.col("_dedup_pat").is_not_null()
        & (~pl.col("_dedup_pat").str.to_uppercase().is_in(["", "NA", "NAN", "NONE", "UNKNOWN"]))
    )

    match mode:
        case "strict_pre":
            # Exclude 'unknown' only inside cohorts where pre/post labels explicitly exist alongside 'unknown'
            cohort_has_explicit_tp = (
                pl.col("_dedup_tp").is_in(["pre", "post", "on-treatment", "on", "treated", "naive"]).any().over("_dedup_coh")
            )
            strict_mask = (
                pl.col("_dedup_tp").is_in(["pre", "baseline", "naive", "untreated"])
                | ((~cohort_has_explicit_tp) & (pl.col("_dedup_tp") == "unknown"))
            )
            work = work.filter(strict_mask).with_columns(pl.lit(0).alias("_tp_rank"))
        case "prefer_pre_one_per_patient":
            tp_rank = (
                pl.when(pl.col("_dedup_tp").is_in(["pre", "baseline", "naive", "untreated"]))
                .then(0)
                .when(pl.col("_dedup_tp") == "unknown")
                .then(1)
                .when(pl.col("_dedup_tp").is_in(["on-treatment", "on"]))
                .then(2)
                .otherwise(3)
            )
            work = work.with_columns(tp_rank.alias("_tp_rank"))
        case _ as unreachable:
            assert_never(unreachable)

    sort_cols = ["_dedup_coh", "_dedup_pat", "_tp_rank"]
    sort_desc = [False, False, False]
    match weight_col:
        case Some(w_c) if w_c in work.columns:
            work = work.with_columns(pl.col(w_c).cast(pl.Float64).fill_null(0.0).alias("_dedup_w"))
            sort_cols.extend(["_dedup_w", id_col])
            sort_desc.extend([True, False])
        case _:
            sort_cols.append(id_col)
            sort_desc.append(False)

    dedup = (
        work.sort(sort_cols, descending=sort_desc)
        .unique(subset=["_dedup_coh", "_dedup_pat"], keep="first", maintain_order=True)
        .drop(["_dedup_pat", "_dedup_coh", "_dedup_tp", "_tp_rank"] + (["_dedup_w"] if "_dedup_w" in work.columns else []))
    )
    return dedup
