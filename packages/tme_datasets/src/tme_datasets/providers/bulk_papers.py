"""Bulk ICI response dataset providers directly from primary publication deposits and GEO."""

from __future__ import annotations

import gzip
from pathlib import Path
import re
from typing import Any, Mapping
import anndata as ad  # type: ignore[import-untyped]
import numpy as np
import pandas as pd
from returns.result import Failure, Result, Success

from ..download.fetcher import download_single_file
from ..logging import get_logger
from ..preprocessing.metadata import binarize_response, standardize_recist

logger = get_logger("providers.bulk_papers")

PAPER_DATASETS = (
    "Auslander",
    "Chen-CTLA4",
    "Chen-PD1",
    "Cho",
    "Freeman",
    "Gide",
    "Hugo",
    "Jung",
    "Lauss",
    "Liu",
    "Prat",
    "Ravi",
    "Riaz",
    "Rose",
    "Snyder",
    "VanAllen",
)

RAVI_ESM_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-023-01355-5/"
    "MediaObjects/41588_2023_1355_MOESM3_ESM.xlsx"
)
RAVI_ESM_FILENAME = "41588_2023_1355_MOESM3_ESM.xlsx"

PRAT_SUPPL_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93157/suppl/GSE93157_raw_data_values.txt.gz"
)
PRAT_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93157/matrix/GSE93157_series_matrix.txt.gz"
)

JUNG_SUPPL_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"
)
JUNG_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz"
)

CHO_SUPPL_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz"
)
CHO_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/matrix/GSE126044_series_matrix.txt.gz"
)

# Jung et al. 2019 (Nat Commun, DOI: 10.1038/s41467-019-12159-9) Supplementary Table 3
# (Cross-verified with Cho et al. 2020 Supplemental Table S6A):
# (mutation_burden, global_methylation_level, aneuploidy_level, pfs_months, pfs_event, clinical_benefit)
_JUNG_SUPPL_TABLE_3: Mapping[str, tuple[int, float, float, float, int, str]] = {
    "378": (117, 0.414296576, 102187138.2, 13.73333333, 0, "DCB"),
    "573": (366, 0.410597593, 114018200.9, 1.5, 1, "NDB"),
    "825": (194, 0.251712972, 784194239.2, 1.433333333, 1, "NDB"),
    "947": (401, 0.387089023, 58799436.69, 20.6, 0, "DCB"),
    "990": (43, 0.325908608, 285011126.5, 3.166666667, 1, "NDB"),
    "1017": (82, 0.417538488, 286678444.4, 1.266666667, 1, "NDB"),
    "1066": (42, 0.395929718, 48979857.76, 0.766666667, 1, "NDB"),
    "1079": (101, 0.392672276, 128159753.6, 5.8, 1, "NDB"),
    "1104": (1100, 0.293149031, 895819215.0, 0.966666667, 1, "NDB"),
    "1145": (178, 0.421131362, 5860926.472, 2.433333333, 1, "NDB"),
    "1155": (273, 0.335438797, 443235770.9, 1.133333333, 1, "NDB"),
    "1164": (175, 0.293888515, 520968472.1, 0.366666667, 1, "NDB"),
    "1203": (137, 0.359713630, 437223965.6, 0.833333333, 1, "NDB"),
    "1327": (441, 0.426781560, 99766614.32, 6.833333333, 0, "DCB"),
    "1352": (1391, 0.220544740, 553574113.5, 8.333333333, 0, "DCB"),
    "1358": (328, 0.377759517, 649708156.1, 8.566666667, 1, "DCB"),
    "1401": (42, 0.436672231, 17494258.78, 0.933333333, 1, "NDB"),
    "1412": (202, 0.264443355, 545902993.4, 0.1, 1, "NDB"),
    "1425": (260, 0.303500801, 346192818.8, 1.966666667, 1, "NDB"),
    "1508": (18, 0.415487136, 7542228.411, 2.166666667, 1, "NDB"),
    "1510": (273, 0.377090295, 465391056.4, 10.8, 0, "DCB"),
    "1528": (146, 0.364790201, 374128096.1, 9.3, 0, "DCB"),
    "1554": (38, 0.400425652, 5947504.686, 1.233333333, 1, "NDB"),
    "1619": (13, 0.455254226, 7041421.973, 2.733333333, 1, "NDB"),
    "1708": (266, 0.450911862, 2076055.089, 5.6, 1, "DCB"),
    "1809": (156, 0.324772186, 400743320.7, 1.466666667, 1, "NDB"),
    "1873": (388, 0.307758327, 486476347.6, 1.233333333, 1, "NDB"),
}

# Cho et al. 2020 (Exp Mol Med, DOI: 10.1038/s12276-020-00493-8) Supplemental Table S2A:
# (age, sex, histology, pdl1_expression, responsiveness, drug, pfs_months, os_months, best_response)
_CHO_SUPPL_TABLE_S2A: Mapping[str, tuple[int, str, str, str, str, str, float, float, str]] = {
    "Dis_01": (72, "male", "Squamous cell carcinoma", "No", "Non-responder", "Nivolumab", 0.5, 1.1, "PD"),
    "Dis_02": (76, "male", "Squamous cell carcinoma", "Yes", "Responder", "Nivolumab", 13.8, 14.1, "PR"),
    "Dis_03": (65, "male", "Squamous cell carcinoma", "NA", "Non-responder", "Nivolumab", 2.3, 11.2, "PD"),
    "Dis_04": (55, "female", "Adenocarcinoma", "Yes", "Responder", "Nivolumab", 13.5, 18.7, "PR"),
    "Dis_05": (72, "male", "Adenocarcinoma", "No", "Non-responder", "Nivolumab", 2.8, 9.6, "PD"),
    "Dis_06": (54, "male", "Adenocarcinoma", "No", "Non-responder", "Nivolumab", 0.8, 0.8, "PD"),
    "Dis_07": (66, "male", "Squamous cell carcinoma", "No", "Non-responder", "Nivolumab", 1.1, 3.5, "PD"),
    "Dis_08": (34, "female", "Adenocarcinoma", "No", "Non-responder", "Nivolumab", 1.0, 1.0, "PD"),
    "Dis_09": (61, "male", "Squamous cell carcinoma", "NA", "Non-responder", "Nivolumab", 0.7, 1.0, "PD"),
    "Dis_10": (75, "male", "Squamous cell carcinoma", "Yes", "Responder", "Nivolumab", 14.2, 14.2, "PR"),
    "Dis_11": (47, "male", "Adenocarcinoma", "No", "Non-responder", "Nivolumab", 2.1, 3.1, "PD"),
    "Dis_12": (65, "male", "Squamous cell carcinoma", "No", "Non-responder", "Nivolumab", 3.4, 12.7, "PD"),
    "Dis_15": (58, "male", "Adenocarcinoma", "No", "Responder", "Nivolumab", 17.7, 18.1, "PR"),
    "Dis_16": (56, "male", "Adenocarcinoma", "Yes", "Non-responder", "Nivolumab", 3.6, 10.2, "PD"),
    "Dis_17": (64, "male", "Squamous cell carcinoma", "No", "Responder", "Nivolumab", 7.3, 12.9, "SD"),
    "Dis_18": (67, "male", "Squamous cell carcinoma", "Yes", "Non-responder", "Nivolumab", 0.5, 1.7, "PD"),
}


def _sanitize_adata_columns(adata: ad.AnnData) -> ad.AnnData:
    """Sanitize object and nullable columns on AnnData .obs and .var for h5py serialization."""
    for col in adata.obs.columns:
        if pd.api.types.is_bool_dtype(adata.obs[col]):
            continue
        if not pd.api.types.is_numeric_dtype(adata.obs[col]):
            adata.obs[col] = adata.obs[col].fillna("nan").astype(str)
    for col in adata.var.columns:
        if not pd.api.types.is_numeric_dtype(adata.var[col]):
            adata.var[col] = adata.var[col].fillna("").astype(str)
    return adata


def _parse_geo_series_matrix_metadata(matrix_path: Path) -> Result[pd.DataFrame, str]:
    """Parse !Sample_* metadata rows from a GEO series_matrix.txt(.gz) file into a DataFrame."""
    if not matrix_path.exists():
        return Failure(f"GEO series matrix file not found: {matrix_path}")

    try:
        open_fn = gzip.open if str(matrix_path).endswith(".gz") else open
        with open_fn(matrix_path, "rt", encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()

        meta: dict[str, list[str]] = {}
        for line in lines:
            if line.startswith(
                ("!Sample_title", "!Sample_geo_accession", "!Sample_source_name_ch1", "!Sample_description")
            ):
                parts = [p.strip().strip('"') for p in line.split("\t")]
                key = parts[0].lstrip("!")
                vals = parts[1:]
                if key not in meta:
                    meta[key] = vals
            elif line.startswith("!Sample_characteristics_ch1"):
                parts = [p.strip().strip('"') for p in line.split("\t")[1:]]
                first_with_colon = next((p for p in parts if ":" in p), "")
                if not first_with_colon:
                    continue
                char_key = first_with_colon.split(":", 1)[0].strip()
                char_vals = [p.split(":", 1)[1].strip() if ":" in p else p for p in parts]
                meta[char_key] = char_vals

        if not meta:
            return Failure(f"No !Sample_* metadata fields found in {matrix_path.name}")

        return Success(pd.DataFrame(meta))
    except Exception as exc:
        return Failure(f"Failed to parse GEO series matrix {matrix_path}: {exc}")


def _ici_agent_to_targets(ici_agent: object) -> set[str]:
    """Map ICI agent string from Ravi et al. 2023 to molecular targets."""
    s = str(ici_agent) if ici_agent is not None else ""
    targets: set[str] = set()
    if re.search(r"Nivolumab|Pembrolizumab", s, re.IGNORECASE):
        targets.add("PD-1")
    if re.search(r"Durvalumab|Avelumab|Atezolizumab", s, re.IGNORECASE):
        targets.add("PD-L1")
    if re.search(r"Ipilimumab|Tremelimumab", s, re.IGNORECASE):
        targets.add("CTLA-4")
    if re.search(r"Lirilumab", s, re.IGNORECASE):
        targets.add("KIR")
    if re.search(r"Urelumab", s, re.IGNORECASE):
        targets.add("CD137")
    if re.search(r"Epacadostat", s, re.IGNORECASE):
        targets.add("IDO")
    if re.search(r"LAG-3", s, re.IGNORECASE):
        targets.add("LAG-3")
    if re.search(r"Carboplatin|Pemetrexed|Gemcitabine", s, re.IGNORECASE):
        targets.add("Cytostasis")
    return targets


# -----------------------------------------------------------------------------
# 1. Ravi et al. 2023 (SU2C-MARK NSCLC Cohort, N = 152)
# -----------------------------------------------------------------------------
def download_ravi(
    raw_dir: Path,
    force_download: bool = False,
) -> Result[Path, str]:
    """Download Ravi et al. 2023 supplementary Excel workbook (41588_2023_1355_MOESM3_ESM.xlsx)."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    dest_xlsx = raw_dir / RAVI_ESM_FILENAME
    if force_download and dest_xlsx.exists():
        dest_xlsx.unlink(missing_ok=True)
    return download_single_file(RAVI_ESM_URL, dest_xlsx)


def preprocess_ravi(xlsx_path_or_dir: Path) -> Result[ad.AnnData, str]:
    """Parse Ravi et al. 2023 Excel workbook into a standardized AnnData object (152 samples)."""
    xlsx_path = (
        xlsx_path_or_dir / RAVI_ESM_FILENAME
        if xlsx_path_or_dir.is_dir()
        else xlsx_path_or_dir
    )
    if not xlsx_path.exists():
        return Failure(f"Ravi supplementary Excel file not found at {xlsx_path}")

    try:
        logger.info("Reading Ravi et al. 2023 RNA TPM table from %s...", xlsx_path.name)
        df_tpm = pd.read_excel(
            xlsx_path,
            sheet_name="Table_S13_RNA_TPM",
            index_col=0,
            skiprows=2,
        )
        descriptions = (
            df_tpm["Description"].astype(str).tolist()
            if "Description" in df_tpm.columns
            else [str(g) for g in df_tpm.index]
        )
        if "Description" in df_tpm.columns:
            df_tpm = df_tpm.drop(columns=["Description"])

        # Samples are columns, genes are index (Ensembl IDs with version suffix)
        raw_gene_ids = [str(g).split(".")[0] for g in df_tpm.index]
        sample_ids = [str(s) for s in df_tpm.columns]
        X_mat = df_tpm.to_numpy(dtype=np.float32).T

        logger.info("Reading Ravi et al. 2023 clinical annotations from %s...", xlsx_path.name)
        raw_clin = pd.read_excel(
            xlsx_path,
            sheet_name="Table_S1_Clinical_Annotations",
            skiprows=2,
        )
        raw_clin = raw_clin.dropna(subset=["Harmonized_SU2C_RNA_Tumor_Sample_ID_v2"]).copy()
        raw_clin["Harmonized_SU2C_RNA_Tumor_Sample_ID_v2"] = (
            raw_clin["Harmonized_SU2C_RNA_Tumor_Sample_ID_v2"].astype(str).str.strip()
        )
        raw_clin = raw_clin.set_index("Harmonized_SU2C_RNA_Tumor_Sample_ID_v2", drop=False)
        raw_clin = raw_clin.reindex(sample_ids)

        recist_to_resp = {"CR": "R", "PR": "R", "SD": "NR", "PD": "NR"}

        obs = pd.DataFrame(index=pd.Index(sample_ids, name="sample_id"))
        if "Patient_Age_at_Diagnosis" in raw_clin.columns:
            obs["age"] = pd.to_numeric(
                raw_clin["Patient_Age_at_Diagnosis"].map(lambda x: "89" if str(x).strip() == ">89" else x),
                errors="coerce",
            )
        else:
            obs["age"] = np.nan

        if "Patient_Sex" in raw_clin.columns:
            obs["sex"] = raw_clin["Patient_Sex"].map(
                lambda x: {"F": "female", "M": "male"}.get(str(x).strip().upper(), str(x))
            )
        else:
            obs["sex"] = "nan"

        if "Histology_Harmonized" in raw_clin.columns:
            obs["histology"] = raw_clin["Histology_Harmonized"].astype(str)
            obs["primary_tumor"] = raw_clin["Histology_Harmonized"].map(
                lambda x: "LUSC" if "squamous" in str(x).lower() and "non" not in str(x).lower() else "LUAD"
            )
        else:
            obs["primary_tumor"] = "LUAD"

        obs["cancer_type"] = "NSCLC"
        obs["overall_survival"] = (
            pd.to_numeric(raw_clin["Harmonized_OS_Days"], errors="coerce")
            if "Harmonized_OS_Days" in raw_clin.columns
            else np.nan
        )
        obs["os_months"] = obs["overall_survival"] / 30.4375
        if "Harmonized_OS_Event" in raw_clin.columns:
            obs["os_event"] = pd.to_numeric(raw_clin["Harmonized_OS_Event"], errors="coerce")

        if "Harmonized_PFS_Days" in raw_clin.columns:
            pfs_days = pd.to_numeric(raw_clin["Harmonized_PFS_Days"], errors="coerce")
            obs["pfs_days"] = pfs_days
            obs["pfs_months"] = pfs_days / 30.4375
        if "Harmonized_PFS_Event" in raw_clin.columns:
            obs["pfs_event"] = pd.to_numeric(raw_clin["Harmonized_PFS_Event"], errors="coerce")

        bor_series = (
            raw_clin["Harmonized_Confirmed_BOR"].astype(str).str.strip()
            if "Harmonized_Confirmed_BOR" in raw_clin.columns
            else pd.Series("Unknown", index=obs.index)
        )
        obs["RECIST"] = bor_series.map(lambda x: x if x in {"CR", "PR", "SD", "PD"} else "nan")
        obs["response"] = bor_series.map(lambda x: recist_to_resp.get(x, "nan"))
        obs["response_binary"] = bor_series.map(binarize_response)
        obs["response_recist"] = bor_series.map(standardize_recist)

        agents = (
            raw_clin["Agent_PD1"].tolist()
            if "Agent_PD1" in raw_clin.columns
            else ["PD-1"] * len(sample_ids)
        )
        obs["treatment.CTLA4"] = ["CTLA-4" in _ici_agent_to_targets(a) for a in agents]
        obs["treatment.PD1"] = ["PD-1" in _ici_agent_to_targets(a) for a in agents]
        obs["treatment.PDL1"] = ["PD-L1" in _ici_agent_to_targets(a) for a in agents]
        obs["treatment.KIR"] = ["KIR" in _ici_agent_to_targets(a) for a in agents]
        obs["treatment.CD137"] = ["CD137" in _ici_agent_to_targets(a) for a in agents]
        obs["treatment.IDO"] = ["IDO" in _ici_agent_to_targets(a) for a in agents]
        obs["treatment.LAG-3"] = ["LAG-3" in _ici_agent_to_targets(a) for a in agents]
        obs["treatment.ACT"] = False
        obs["treatment.CYTOSTASIS"] = ["Cytostasis" in _ici_agent_to_targets(a) for a in agents]
        obs["timing"] = "Pre"
        obs["original.patient_id"] = (
            raw_clin["Harmonized_SU2C_Participant_ID_v2"].astype(str)
            if "Harmonized_SU2C_Participant_ID_v2" in raw_clin.columns
            else sample_ids
        )
        obs["original.sample_id"] = sample_ids

        var = pd.DataFrame(index=pd.Index(raw_gene_ids, name="gene_id"))
        var["Description"] = descriptions
        var["gene_name"] = descriptions

        adata = ad.AnnData(X=X_mat, obs=obs, var=var)
        adata.uns["expression_type"] = "tpm"
        return Success(_sanitize_adata_columns(adata))
    except Exception as exc:
        msg = f"Failed to preprocess Ravi et al. 2023 cohort from {xlsx_path}: {exc}"
        logger.error(msg)
        return Failure(msg)


def load_ravi(
    raw_dir: Path,
    auto_download: bool = True,
    force_download: bool = False,
) -> Result[ad.AnnData, str]:
    """Download (if needed) and preprocess Ravi et al. 2023 NSCLC cohort."""
    xlsx_path = raw_dir / RAVI_ESM_FILENAME
    if (auto_download or force_download) and (force_download or not xlsx_path.exists()):
        dl_res = download_ravi(raw_dir, force_download=force_download)
        if isinstance(dl_res, Failure):
            return dl_res
    return preprocess_ravi(xlsx_path)


# -----------------------------------------------------------------------------
# 2. Prat et al. 2017 (GSE93157 Pan-Cancer / NSCLC PD-1 Cohort, N = 65)
# -----------------------------------------------------------------------------
def download_prat(
    raw_dir: Path,
    force_download: bool = False,
) -> Result[tuple[Path, Path], str]:
    """Download Prat et al. 2017 (GSE93157) raw expression and series matrix files from GEO."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    suppl_path = raw_dir / "GSE93157_raw_data_values.txt.gz"
    matrix_path = raw_dir / "GSE93157_series_matrix.txt.gz"
    if force_download:
        suppl_path.unlink(missing_ok=True)
        matrix_path.unlink(missing_ok=True)

    res_suppl = download_single_file(PRAT_SUPPL_URL, suppl_path)
    if isinstance(res_suppl, Failure):
        return res_suppl
    res_mat = download_single_file(PRAT_MATRIX_URL, matrix_path)
    if isinstance(res_mat, Failure):
        return res_mat
    return Success((suppl_path, matrix_path))


def preprocess_prat(raw_dir: Path) -> Result[ad.AnnData, str]:
    """Parse Prat et al. 2017 (GSE93157) files into a standardized AnnData object (65 samples)."""
    suppl_path = raw_dir / "GSE93157_raw_data_values.txt.gz"
    matrix_path = raw_dir / "GSE93157_series_matrix.txt.gz"
    legacy_xls = raw_dir / "00085472can163556-sup-176636_2_unknown_upload_3997093_ddtlnx.xls"

    try:
        if suppl_path.exists() and matrix_path.exists():
            logger.info("Parsing Prat et al. 2017 GEO files from %s...", raw_dir)
            open_fn = gzip.open if str(suppl_path).endswith(".gz") else open
            with open_fn(suppl_path, "rt", encoding="utf-8", errors="replace") as f:
                raw_lines = [
                    line.rstrip("\r\n")
                    for line in f
                    if line.strip() and not line.lstrip().startswith(("#", '"#'))
                ]

            header_cols = [c.strip() for c in raw_lines[0].split("\t") if c.strip()]
            sample_cols = header_cols[1:]  # SAMPLE 1 .. SAMPLE 65
            n_samples = len(sample_cols)

            genes: list[str] = []
            rows: list[list[float]] = []
            for line in raw_lines[1:]:
                parts = line.split("\t")
                gene_sym = parts[0].strip()
                if not gene_sym:
                    continue
                vals = [
                    float(parts[i]) if i < len(parts) and parts[i].strip() not in ("", "null", "NA") else 0.0
                    for i in range(1, n_samples + 1)
                ]
                genes.append(gene_sym)
                rows.append(vals)

            X_mat = np.array(rows, dtype=np.float32).T  # (65, 770)

            meta_res = _parse_geo_series_matrix_metadata(matrix_path)
            if isinstance(meta_res, Failure):
                return meta_res
            meta_df = meta_res.unwrap()

            # Align SAMPLE k with Patientk
            meta_df["sample_num"] = [
                int(re.sub(r"\D+", "", str(t))) if re.sub(r"\D+", "", str(t)) else idx + 1
                for idx, t in enumerate(meta_df.get("Sample_title", pd.Series(range(1, n_samples + 1))))
            ]
            meta_df = meta_df.set_index("sample_num", drop=False)
            sample_nums = [
                int(re.sub(r"\D+", "", s)) if re.sub(r"\D+", "", s) else idx + 1
                for idx, s in enumerate(sample_cols)
            ]
            meta_df = meta_df.reindex(sample_nums)

            patient_ids = [
                str(t) for t in meta_df.get("Sample_title", pd.Series([f"Patient{n}" for n in sample_nums]))
            ]
            obs = pd.DataFrame(index=pd.Index(patient_ids, name="sample_id"))
            obs["geo_accession"] = meta_df.get("Sample_geo_accession", pd.Series("", index=meta_df.index)).values
            obs["age"] = pd.to_numeric(meta_df.get("age", pd.Series(np.nan, index=meta_df.index)), errors="coerce").values

            sex_raw = meta_df.get("Sex", pd.Series("nan", index=meta_df.index))
            obs["sex"] = [
                {"F": "female", "M": "male"}.get(str(x).strip().upper(), str(x)) for x in sex_raw
            ]

            source_raw = meta_df.get("Sample_source_name_ch1", pd.Series("", index=meta_df.index))
            tumor_map = {
                "LUNG NON-SQUAMOUS CANCER": "LUAD",
                "SQUAMOUS LUNG CANCER": "LUSC",
                "MELANOMA": "SKCM",
                "HEADNECK": "HNSC",
            }
            cancer_map = {
                "LUNG NON-SQUAMOUS CANCER": "NSCLC",
                "SQUAMOUS LUNG CANCER": "NSCLC",
                "MELANOMA": "Melanoma",
                "HEADNECK": "HNSCC",
            }
            obs["primary_tumor"] = [tumor_map.get(str(x).strip().upper(), str(x)) for x in source_raw]
            obs["cancer_type"] = [cancer_map.get(str(x).strip().upper(), str(x)) for x in source_raw]
            obs["overall_survival"] = np.nan

            best_resp = [str(x).strip().upper() for x in meta_df.get("best.resp", pd.Series("Unknown", index=meta_df.index))]
            recist_to_resp = {"CR": "R", "PR": "R", "SD": "NR", "PD": "NR"}
            obs["RECIST"] = [r if r in recist_to_resp else "nan" for r in best_resp]
            obs["response"] = [recist_to_resp.get(r, "nan") for r in best_resp]
            obs["response_binary"] = [binarize_response(r) for r in best_resp]
            obs["response_recist"] = [standardize_recist(r) for r in best_resp]

            if "pfs" in meta_df.columns:
                obs["pfs_months"] = pd.to_numeric(meta_df["pfs"], errors="coerce").values
            if "pfse" in meta_df.columns:
                obs["pfs_event"] = pd.to_numeric(meta_df["pfse"], errors="coerce").values

            for extra_col in ("drug", "biopsy", "smoking", "ecog", "prev. lines", "braf", "nras", "egfr", "alk"):
                if extra_col in meta_df.columns:
                    obs[extra_col.replace(". ", "_")] = meta_df[extra_col].astype(str).values

            obs["treatment.CTLA4"] = False
            obs["treatment.PD1"] = True
            obs["treatment.PDL1"] = False
            obs["timing"] = "Pre"
            obs["sample_source"] = obs.get("biopsy", "ARCHIVAL")
            obs["original.patient_id"] = patient_ids
            obs["original.sample_id"] = patient_ids

            var = pd.DataFrame(index=pd.Index(genes, name="gene_name"))
            var["gene_name"] = genes
            adata = ad.AnnData(X=X_mat, obs=obs, var=var)
            return Success(_sanitize_adata_columns(adata))

        if legacy_xls.exists():
            logger.info("Parsing Prat et al. 2017 legacy Excel file from %s...", legacy_xls.name)
            data = pd.read_excel(legacy_xls, sheet_name="raw data", index_col=0).T
            obs_raw = pd.read_excel(legacy_xls, sheet_name="clinical data", index_col=0)
            var_df = data.loc[["Accession", "Class Name"]].T
            data = data.drop(["Accession", "Class Name"])
            obs_raw["original.sample_id"] = obs_raw.index.astype(str)
            obs_raw["original.patient_id"] = obs_raw.index.astype(str)
            obs_raw["response_binary"] = obs_raw["BEST.RESP"].map(binarize_response)
            obs_raw["response_recist"] = obs_raw["BEST.RESP"].map(standardize_recist)
            adata = ad.AnnData(X=data.to_numpy(dtype=np.float32), obs=obs_raw, var=var_df)
            return Success(_sanitize_adata_columns(adata))

        return Failure(f"Prat GSE93157 files not found in {raw_dir}")
    except Exception as exc:
        msg = f"Failed to preprocess Prat et al. 2017 (GSE93157) from {raw_dir}: {exc}"
        logger.error(msg)
        return Failure(msg)


def load_prat(
    raw_dir: Path,
    auto_download: bool = True,
    force_download: bool = False,
) -> Result[ad.AnnData, str]:
    """Download (if needed) and preprocess Prat et al. 2017 (GSE93157) cohort."""
    suppl_path = raw_dir / "GSE93157_raw_data_values.txt.gz"
    matrix_path = raw_dir / "GSE93157_series_matrix.txt.gz"
    if (auto_download or force_download) and (
        force_download or not suppl_path.exists() or not matrix_path.exists()
    ):
        dl_res = download_prat(raw_dir, force_download=force_download)
        if isinstance(dl_res, Failure):
            return dl_res
    return preprocess_prat(raw_dir)


# -----------------------------------------------------------------------------
# 3. Jung et al. 2019 (GSE135222 NSCLC Anti-PD-1/PD-L1 Cohort, N = 27)
# -----------------------------------------------------------------------------
def download_jung(
    raw_dir: Path,
    force_download: bool = False,
) -> Result[tuple[Path, Path], str]:
    """Download Jung et al. 2019 (GSE135222) RNA-seq TPM and series matrix files from GEO."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    suppl_path = raw_dir / "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"
    matrix_path = raw_dir / "GSE135222_series_matrix.txt.gz"
    if force_download:
        suppl_path.unlink(missing_ok=True)
        matrix_path.unlink(missing_ok=True)

    res_suppl = download_single_file(JUNG_SUPPL_URL, suppl_path)
    if isinstance(res_suppl, Failure):
        return res_suppl
    res_mat = download_single_file(JUNG_MATRIX_URL, matrix_path)
    if isinstance(res_mat, Failure):
        return res_mat
    return Success((suppl_path, matrix_path))


def preprocess_jung(raw_dir: Path) -> Result[ad.AnnData, str]:
    """Parse Jung et al. 2019 (GSE135222) into a standardized AnnData object (27 NSCLC samples)."""
    suppl_path = raw_dir / "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"
    matrix_path = raw_dir / "GSE135222_series_matrix.txt.gz"

    if not suppl_path.exists():
        return Failure(f"Jung GSE135222 expression matrix not found at {suppl_path}")

    try:
        logger.info("Parsing Jung et al. 2019 (GSE135222) RNA-seq expression from %s...", suppl_path.name)
        df_expr = pd.read_csv(suppl_path, sep="\t", index_col=0)
        raw_gene_ids = [str(g).split(".")[0] for g in df_expr.index]
        sample_ids = [str(s).strip() for s in df_expr.columns]
        X_mat = df_expr.to_numpy(dtype=np.float32).T

        obs = pd.DataFrame(index=pd.Index(sample_ids, name="sample_id"))
        patient_nums = [re.sub(r"^NSCLC\s*", "", s).strip() for s in sample_ids]
        obs["patient_number"] = patient_nums

        if matrix_path.exists():
            meta_res = _parse_geo_series_matrix_metadata(matrix_path)
            if isinstance(meta_res, Success):
                meta_df = meta_res.unwrap()
                meta_df["clean_id"] = [
                    str(t).replace(" ", "").strip()
                    for t in meta_df.get("Sample_title", pd.Series(sample_ids))
                ]
                meta_df = meta_df.set_index("clean_id", drop=False).reindex(sample_ids)
                obs["geo_accession"] = meta_df.get("Sample_geo_accession", pd.Series("", index=obs.index)).values
                obs["age"] = pd.to_numeric(meta_df.get("age", pd.Series(np.nan, index=obs.index)), errors="coerce").values
                obs["sex"] = [
                    str(g).strip().lower()
                    for g in meta_df.get("gender", pd.Series("nan", index=obs.index))
                ]
                if "pfs.time" in meta_df.columns:
                    obs["pfs_days"] = pd.to_numeric(meta_df["pfs.time"], errors="coerce").values
                if "progression-free survival (pfs)" in meta_df.columns:
                    obs["pfs_event"] = pd.to_numeric(
                        meta_df["progression-free survival (pfs)"], errors="coerce"
                    ).values

        # Enrich with Jung et al. 2019 Supplementary Table 3 clinical benefit (DCB vs NDB) and genomic metrics
        tmb_list: list[float] = []
        meth_list: list[float] = []
        aneu_list: list[float] = []
        pfs_mo_list: list[float] = []
        pfs_ev_list: list[float] = []
        benefit_list: list[str] = []

        for idx, pnum in enumerate(patient_nums):
            if pnum in _JUNG_SUPPL_TABLE_3:
                tmb, meth, aneu, pfs_mo, pfs_ev, benefit = _JUNG_SUPPL_TABLE_3[pnum]
                tmb_list.append(float(tmb))
                meth_list.append(float(meth))
                aneu_list.append(float(aneu))
                pfs_mo_list.append(float(pfs_mo))
                pfs_ev_list.append(float(pfs_ev))
                benefit_list.append(benefit)
            else:
                tmb_list.append(np.nan)
                meth_list.append(np.nan)
                aneu_list.append(np.nan)
                pfs_days_val = float(obs["pfs_days"].iloc[idx]) if "pfs_days" in obs.columns else np.nan
                pfs_mo_val = pfs_days_val / 30.0 if not np.isnan(pfs_days_val) else np.nan
                pfs_mo_list.append(pfs_mo_val)
                pfs_ev_val = float(obs["pfs_event"].iloc[idx]) if "pfs_event" in obs.columns else np.nan
                pfs_ev_list.append(pfs_ev_val)
                benefit_list.append("DCB" if not np.isnan(pfs_mo_val) and pfs_mo_val >= 6.0 else "NDB")

        obs["mutation_burden"] = tmb_list
        obs["global_methylation_level"] = meth_list
        obs["aneuploidy_level"] = aneu_list
        obs["pfs_months"] = pfs_mo_list
        if "pfs_event" not in obs.columns:
            obs["pfs_event"] = pfs_ev_list
        obs["clinical_benefit"] = benefit_list
        obs["response"] = ["R" if b == "DCB" else ("NR" if b == "NDB" else "nan") for b in benefit_list]
        obs["RECIST"] = "nan"
        obs["response_binary"] = [binarize_response(b) for b in benefit_list]
        obs["response_recist"] = "Unknown"
        obs["primary_tumor"] = "NSCLC"
        obs["cancer_type"] = "NSCLC"
        obs["treatment.PD1"] = True
        obs["treatment.PDL1"] = True
        obs["treatment.CTLA4"] = False
        obs["timing"] = "Pre"
        obs["original.patient_id"] = sample_ids
        obs["original.sample_id"] = sample_ids

        var = pd.DataFrame(index=pd.Index(raw_gene_ids, name="gene_id"))
        var["original_ensembl_id"] = [str(g) for g in df_expr.index]

        adata = ad.AnnData(X=X_mat, obs=obs, var=var)
        adata.uns["expression_type"] = "tpm"
        return Success(_sanitize_adata_columns(adata))
    except Exception as exc:
        msg = f"Failed to preprocess Jung et al. 2019 (GSE135222) from {raw_dir}: {exc}"
        logger.error(msg)
        return Failure(msg)


def load_jung(
    raw_dir: Path,
    auto_download: bool = True,
    force_download: bool = False,
) -> Result[ad.AnnData, str]:
    """Download (if needed) and preprocess Jung et al. 2019 (GSE135222) NSCLC cohort."""
    suppl_path = raw_dir / "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"
    matrix_path = raw_dir / "GSE135222_series_matrix.txt.gz"
    if (auto_download or force_download) and (
        force_download or not suppl_path.exists() or not matrix_path.exists()
    ):
        dl_res = download_jung(raw_dir, force_download=force_download)
        if isinstance(dl_res, Failure):
            return dl_res
    return preprocess_jung(raw_dir)


# -----------------------------------------------------------------------------
# 4. Cho et al. 2020 (GSE126044 NSCLC Anti-PD-1 Cohort, N = 16)
# -----------------------------------------------------------------------------
def download_cho(
    raw_dir: Path,
    force_download: bool = False,
) -> Result[tuple[Path, Path], str]:
    """Download Cho et al. 2020 (GSE126044) raw counts and series matrix files from GEO."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    suppl_path = raw_dir / "GSE126044_counts.txt.gz"
    matrix_path = raw_dir / "GSE126044_series_matrix.txt.gz"
    if force_download:
        suppl_path.unlink(missing_ok=True)
        matrix_path.unlink(missing_ok=True)

    res_suppl = download_single_file(CHO_SUPPL_URL, suppl_path)
    if isinstance(res_suppl, Failure):
        return res_suppl
    res_mat = download_single_file(CHO_MATRIX_URL, matrix_path)
    if isinstance(res_mat, Failure):
        return res_mat
    return Success((suppl_path, matrix_path))


def preprocess_cho(raw_dir: Path) -> Result[ad.AnnData, str]:
    """Parse Cho et al. 2020 (GSE126044) into a standardized AnnData object (16 NSCLC samples)."""
    suppl_path = raw_dir / "GSE126044_counts.txt.gz"
    matrix_path = raw_dir / "GSE126044_series_matrix.txt.gz"

    if not suppl_path.exists():
        return Failure(f"Cho GSE126044 count matrix not found at {suppl_path}")

    try:
        logger.info("Parsing Cho et al. 2020 (GSE126044) RNA-seq counts from %s...", suppl_path.name)
        df_counts = pd.read_csv(suppl_path, sep="\t", index_col=0)
        genes = [str(g).strip() for g in df_counts.index]
        sample_ids = [str(s).strip() for s in df_counts.columns]
        X_mat = df_counts.to_numpy(dtype=np.float32).T

        obs = pd.DataFrame(index=pd.Index(sample_ids, name="sample_id"))

        if matrix_path.exists():
            meta_res = _parse_geo_series_matrix_metadata(matrix_path)
            if isinstance(meta_res, Success):
                meta_df = meta_res.unwrap()
                meta_df["clean_id"] = [
                    re.sub(r"^RNA-seq_", "", str(t)).strip()
                    for t in meta_df.get("Sample_title", pd.Series(sample_ids))
                ]
                meta_df = meta_df.set_index("clean_id", drop=False).reindex(sample_ids)
                obs["geo_accession"] = meta_df.get("Sample_geo_accession", pd.Series("", index=obs.index)).values
                obs["patient_response_geo"] = meta_df.get(
                    "patient response", pd.Series("", index=obs.index)
                ).values
                obs["sample_source"] = meta_df.get("sample", pd.Series("", index=obs.index)).values

        ages: list[float] = []
        sexes: list[str] = []
        histologies: list[str] = []
        primary_tumors: list[str] = []
        pdl1_list: list[str] = []
        drugs: list[str] = []
        pfs_mo_list: list[float] = []
        os_mo_list: list[float] = []
        recist_list: list[str] = []
        resp_list: list[str] = []
        resp_bin_list: list[float] = []

        for idx, sid in enumerate(sample_ids):
            if sid in _CHO_SUPPL_TABLE_S2A:
                age, sex, hist, pdl1, resp_str, drug, pfs_mo, os_mo, bor = _CHO_SUPPL_TABLE_S2A[sid]
                ages.append(float(age))
                sexes.append(sex)
                histologies.append(hist)
                primary_tumors.append("LUSC" if "squamous" in hist.lower() else "LUAD")
                pdl1_list.append(pdl1)
                drugs.append(drug)
                pfs_mo_list.append(float(pfs_mo))
                os_mo_list.append(float(os_mo))
                recist_list.append(bor)
                resp_list.append("R" if resp_str.lower() == "responder" else "NR")
                resp_bin_list.append(binarize_response(resp_str))
            else:
                ages.append(np.nan)
                sexes.append("nan")
                histologies.append("NSCLC")
                primary_tumors.append("NSCLC")
                pdl1_list.append("NA")
                drugs.append("Nivolumab")
                pfs_mo_list.append(np.nan)
                os_mo_list.append(np.nan)
                recist_list.append("Unknown")
                geo_resp = (
                    str(obs["patient_response_geo"].iloc[idx])
                    if "patient_response_geo" in obs.columns
                    else ""
                )
                b_val = binarize_response(geo_resp)
                resp_bin_list.append(b_val)
                resp_list.append("R" if b_val == 1.0 else ("NR" if b_val == 0.0 else "nan"))

        obs["age"] = ages
        obs["sex"] = sexes
        obs["histology"] = histologies
        obs["primary_tumor"] = primary_tumors
        obs["cancer_type"] = "NSCLC"
        obs["pdl1_expression"] = pdl1_list
        obs["drug"] = drugs
        obs["pfs_months"] = pfs_mo_list
        obs["os_months"] = os_mo_list
        obs["overall_survival"] = [m * 30.4375 if not np.isnan(m) else np.nan for m in os_mo_list]
        obs["RECIST"] = recist_list
        obs["response"] = resp_list
        obs["response_binary"] = resp_bin_list
        obs["response_recist"] = [standardize_recist(r) for r in recist_list]
        obs["treatment.PD1"] = True
        obs["treatment.PDL1"] = False
        obs["treatment.CTLA4"] = False
        obs["timing"] = "Pre"
        obs["original.patient_id"] = sample_ids
        obs["original.sample_id"] = sample_ids

        var = pd.DataFrame(index=pd.Index(genes, name="gene_name"))
        var["gene_name"] = genes

        adata = ad.AnnData(X=X_mat, obs=obs, var=var)
        adata.uns["expression_type"] = "counts"
        return Success(_sanitize_adata_columns(adata))
    except Exception as exc:
        msg = f"Failed to preprocess Cho et al. 2020 (GSE126044) from {raw_dir}: {exc}"
        logger.error(msg)
        return Failure(msg)


def load_cho(
    raw_dir: Path,
    auto_download: bool = True,
    force_download: bool = False,
) -> Result[ad.AnnData, str]:
    """Download (if needed) and preprocess Cho et al. 2020 (GSE126044) NSCLC cohort."""
    suppl_path = raw_dir / "GSE126044_counts.txt.gz"
    matrix_path = raw_dir / "GSE126044_series_matrix.txt.gz"
    if (auto_download or force_download) and (
        force_download or not suppl_path.exists() or not matrix_path.exists()
    ):
        dl_res = download_cho(raw_dir, force_download=force_download)
        if isinstance(dl_res, Failure):
            return dl_res
    return preprocess_cho(raw_dir)


# -----------------------------------------------------------------------------
# Existing Pre-Built Paper H5AD & Genentech EGAD Loaders
# -----------------------------------------------------------------------------
def load_paper_h5ad(h5ad_path: Path) -> Result[ad.AnnData, str]:
    """Load an authentic paper H5AD dataset from publication archives or data/manual_download/."""
    if not h5ad_path.exists():
        msg = (
            f"Paper H5AD file not found at '{h5ad_path}'. "
            f"This dataset requires manual placement into 'data/manual_download/{h5ad_path.name}'. "
            f"If working on the cluster, run 'python scripts/setup_manual_downloads.py' (or 'make setup-manual-downloads') "
            f"to copy paper cohorts from cluster storage into data/manual_download/."
        )
        logger.error(msg)
        return Failure(msg)

    try:
        logger.info("Reading paper H5AD file from %s...", h5ad_path.name)
        adata = ad.read_h5ad(h5ad_path)
        # Harmonize response column if present
        resp_cols = [c for c in adata.obs.columns if any(k in c.lower() for k in ("response", "recist", "dcb"))]
        if resp_cols and "response_binary" not in adata.obs.columns:
            resp_col = resp_cols[0]
            logger.debug("Standardizing response column '%s' for %s", resp_col, h5ad_path.name)
            adata.obs["response_binary"] = adata.obs[resp_col].apply(binarize_response)
            adata.obs["response_recist"] = adata.obs[resp_col].apply(standardize_recist)

        logger.info("Successfully parsed %s: %d samples x %d genes", h5ad_path.name, adata.n_obs, adata.n_vars)
        return Success(adata)
    except Exception as exc:
        msg = f"Failed to load paper H5AD from {h5ad_path}: {exc}"
        logger.error(msg)
        return Failure(msg)


def load_genentech_egad(align_dir: Path) -> Result[ad.AnnData, str]:
    """Load Genentech IMvigor210 raw RNA-seq alignment read counts (EGAD00001006631)."""
    if not align_dir.exists():
        msg = (
            f"EGAD alignment directory not found at '{align_dir}'. "
            f"This dataset requires controlled-access / manual download. "
            f"Please ensure the alignment directory is placed at 'data/manual_download/EGAD00001006631-align' "
            f"(or run 'python scripts/setup_manual_downloads.py' on the cluster)."
        )
        logger.error(msg)
        return Failure(msg)

    try:
        logger.info("Scanning patient subdirectories in %s...", align_dir)
        all_counts = []
        for dirpath in align_dir.iterdir():
            if not dirpath.is_dir():
                continue
            tab_file = dirpath / "ReadsPerGene.out.tab"
            if tab_file.exists():
                patient_id = dirpath.name
                counts = pd.read_csv(tab_file, sep="\t", header=None, skiprows=4, index_col=0)
                # Select stranded_reverse (col 3)
                series = counts[3]
                series.name = patient_id
                all_counts.append(series)

        if not all_counts:
            msg = f"No ReadsPerGene.out.tab files found in {align_dir}"
            logger.error(msg)
            return Failure(msg)

        logger.info("Concatenating %d patient alignment count profiles...", len(all_counts))
        df_counts = pd.concat(all_counts, axis=1).T
        adata = ad.AnnData(
            X=df_counts.values.astype(np.float32),
            obs=pd.DataFrame(index=df_counts.index),
            var=pd.DataFrame(index=df_counts.columns),
        )
        return Success(adata)
    except Exception as exc:
        msg = f"Failed to load Genentech EGAD dataset: {exc}"
        logger.error(msg)
        return Failure(msg)
