"""Bulk ICI validation cohort providers via cBioPortal / iAtlas."""

from __future__ import annotations

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from returns.result import Failure, Result, Success

from ..logging import get_logger
from ..preprocessing.metadata import binarize_response, standardize_recist

logger = get_logger("providers.bulk_iatlas")

IATLAS_COHORTS = (
    "Hugo-iAtlas",
    "Riaz-iAtlas",
    "Liu-iAtlas",
    "Gide-iAtlas",
    "Rosenberg-iAtlas",
    "Padron-iAtlas",
    "Anders-iAtlas",
    "McDermott-iAtlas",
    "Choueiri-iAtlas",
)


def load_iatlas_cohort(
    cohort_dir: Path,
    cohort_name: str | None = None,
    auto_download: bool = True,
    force_download: bool = False,
) -> Result[ad.AnnData, str]:
    """Load a single cBioPortal / iAtlas bulk cohort directory into an AnnData object."""
    # cBioPortal tar.gz archive mapping
    cbioportal_archives = {
        "Hugo-iAtlas": "mel_iatlas_hugo_ucla_2016.tar.gz",
        "Riaz-iAtlas": "mel_iatlas_riaz_nivolumab_2017.tar.gz",
        "Liu-iAtlas": "mel_iatlas_liu_2019.tar.gz",
        "Gide-iAtlas": "mel_iatlas_gide_2019.tar.gz",
        "Rosenberg-iAtlas": "blca_iatlas_imvigor210_2017.tar.gz",
        "Padron-iAtlas": "paad_iatlas_prince_2022.tar.gz",
        "Anders-iAtlas": "brca_iatlas_anders_2022.tar.gz",
        "McDermott-iAtlas": "rcc_iatlas_immotion150_2018.tar.gz",
        "Choueiri-iAtlas": "ccrcc_iatlas_choueiri_2016.tar.gz",
        "Cloughesy-iAtlas": "gbm_iatlas_prins_2019.tar.gz",
    }

    expr_files = [
        "data_mrna_seq_tpm.txt",
        "data_mrna_seq_expression.txt",
        "data_mrna_seq_rpkm.txt",
        "data_mrna_seq_fpkm.txt",
    ]

    # Auto-download if missing or force_download=True
    if (auto_download or force_download) and cohort_name in cbioportal_archives:
        tar_name = cbioportal_archives[cohort_name]
        tar_path = cohort_dir / tar_name
        # Check if already extracted
        has_expr = cohort_dir.exists() and any(
            any(cohort_dir.rglob(f)) for f in expr_files
        )
        if force_download or not has_expr:
            from ..download.fetcher import download_single_file, unpack_tar

            url = f"https://datahub.assets.cbioportal.org/{tar_name}"
            if force_download and tar_path.exists():
                tar_path.unlink()
            download_single_file(url, tar_path)
            if tar_path.exists():
                unpack_tar(tar_path, cohort_dir)

    if not cohort_dir.is_dir():
        return Failure(f"Directory does not exist: {cohort_dir}")

    try:
        # Search directly in cohort_dir or any nested unpacked folder
        found_expr = None
        for f in expr_files:
            matches = list(cohort_dir.rglob(f))
            if matches:
                found_expr = matches[0]
                break

        if not found_expr:
            msg = f"No mRNA expression file found in {cohort_dir}"
            logger.error(msg)
            return Failure(msg)

        target_dir = found_expr.parent
        logger.info("Reading mRNA expression file from %s...", found_expr.name)

        df_expr = pd.read_csv(found_expr, sep="\t", index_col=0)
        # Drop redundant Hugo symbol or Entrez column if present
        if "Entrez_Gene_Id" in df_expr.columns:
            df_expr = df_expr.drop(columns=["Entrez_Gene_Id"])

        # Samples are columns, genes are index
        sample_ids = list(df_expr.columns)
        genes = list(df_expr.index)
        logger.info(
            "Parsed %d samples across %d genes from %s",
            len(sample_ids),
            len(genes),
            found_expr.name,
        )

        # 2. Load clinical annotations
        clinical_sample_file = next(
            (p for p in (target_dir / "data_clinical_sample.txt", cohort_dir / "data_clinical_sample.txt") if p.exists()),
            None,
        )
        if not clinical_sample_file:
            clin_matches = list(cohort_dir.rglob("data_clinical_sample.txt"))
            clinical_sample_file = clin_matches[0] if clin_matches else None

        df_clinical = pd.DataFrame(index=sample_ids)

        if clinical_sample_file and clinical_sample_file.exists():
            logger.info("Parsing clinical annotations from %s...", clinical_sample_file.name)
            # Skip comment lines beginning with '#'
            df_raw_clin = pd.read_csv(clinical_sample_file, sep="\t", comment="#")
            id_col = next((c for c in df_raw_clin.columns if any(k in c.upper() for k in ("SAMPLE_ID", "SAMPLEID"))), df_raw_clin.columns[0])
            df_raw_clin[id_col] = df_raw_clin[id_col].astype(str).str.strip()

            # Merge data_clinical_patient.txt if present (adds OS, PFS, and drug target metadata)
            clinical_patient_file = next(
                (p for p in (target_dir / "data_clinical_patient.txt", cohort_dir / "data_clinical_patient.txt") if p.exists()),
                None,
            )
            if not clinical_patient_file:
                pat_matches = list(cohort_dir.rglob("data_clinical_patient.txt"))
                clinical_patient_file = pat_matches[0] if pat_matches else None

            if clinical_patient_file and clinical_patient_file.exists():
                logger.info("Parsing patient annotations from %s...", clinical_patient_file.name)
                df_raw_pat = pd.read_csv(clinical_patient_file, sep="\t", comment="#")
                pat_id_col = next(
                    (c for c in df_raw_pat.columns if any(k in c.upper() for k in ("PATIENT_ID", "PATIENTID"))),
                    df_raw_pat.columns[0],
                )
                df_raw_pat[pat_id_col] = df_raw_pat[pat_id_col].astype(str).str.strip()

                clin_pat_col = next(
                    (c for c in df_raw_clin.columns if any(k in c.upper() for k in ("PATIENT_ID", "PATIENTID"))),
                    None,
                )
                if clin_pat_col:
                    df_raw_clin[clin_pat_col] = df_raw_clin[clin_pat_col].astype(str).str.strip()
                    df_raw_pat[pat_id_col] = df_raw_pat[pat_id_col].astype(str).str.strip()
                    df_raw_clin = df_raw_clin.merge(
                        df_raw_pat,
                        left_on=clin_pat_col,
                        right_on=pat_id_col,
                        how="left",
                        suffixes=("", "_patient"),
                    )
                else:
                    df_raw_clin[id_col] = df_raw_clin[id_col].astype(str).str.strip()
                    df_raw_pat[pat_id_col] = df_raw_pat[pat_id_col].astype(str).str.strip()
                    df_raw_clin = df_raw_clin.merge(
                        df_raw_pat,
                        left_on=id_col,
                        right_on=pat_id_col,
                        how="left",
                        suffixes=("", "_patient"),
                    )

            df_raw_clin = df_raw_clin.set_index(id_col)
            # Reindex to match expression samples
            df_clinical.index = df_clinical.index.astype(str).str.strip()
            common_samples = df_clinical.index.intersection(df_raw_clin.index)
            df_clinical.loc[common_samples, df_raw_clin.columns] = df_raw_clin.loc[common_samples]

        # 3. Standardize response
        resp_col = next((c for c in df_clinical.columns if any(k in c.upper() for k in ("RESPONSE", "RECIST", "CLINICAL_BENEFIT"))), None)
        if resp_col:
            df_clinical["response_binary"] = df_clinical[resp_col].apply(binarize_response)
            df_clinical["response_recist"] = df_clinical[resp_col].apply(standardize_recist)

        # 4. Standardize survival endpoints (OS & PFS)
        if "OS_STATUS" in df_clinical.columns:
            def _parse_os_status(val: object) -> float:
                if pd.isna(val):
                    return np.nan
                s = str(val).upper().strip()
                if "DECEASED" in s or "DEAD" in s or s.startswith("1:"):
                    return 1.0
                if "LIVING" in s or "ALIVE" in s or s.startswith("0:"):
                    return 0.0
                return np.nan
            df_clinical["os_event"] = df_clinical["OS_STATUS"].apply(_parse_os_status)

        if "OS_MONTHS" in df_clinical.columns:
            df_clinical["os_months"] = pd.to_numeric(df_clinical["OS_MONTHS"], errors="coerce")

        if "PFS_STATUS" in df_clinical.columns:
            def _parse_pfs_status(val: object) -> float:
                if pd.isna(val):
                    return np.nan
                s = str(val).upper().strip()
                if "PROGRESSED" in s or s.startswith("1:"):
                    return 1.0
                if "NOT_PROGRESSED" in s or "CENSORED" in s or "LIVING" in s or s.startswith("0:"):
                    return 0.0
                return np.nan
            df_clinical["pfs_event"] = df_clinical["PFS_STATUS"].apply(_parse_pfs_status)

        if "PFS_MONTHS" in df_clinical.columns:
            df_clinical["pfs_months"] = pd.to_numeric(df_clinical["PFS_MONTHS"], errors="coerce")

        # 5. Sanitize object columns for robust H5AD serialization
        for col in df_clinical.columns:
            if df_clinical[col].dtype == object:
                if any(isinstance(v, bool) for v in df_clinical[col].dropna()):
                    df_clinical[col] = df_clinical[col].astype("boolean")
                else:
                    df_clinical[col] = df_clinical[col].astype(str)

        X_mat = df_expr.values.T.astype(np.float32)
        adata = ad.AnnData(
            X=X_mat,
            obs=df_clinical,
            var=pd.DataFrame(index=genes),
        )
        return Success(adata)
    except Exception as exc:
        return Failure(f"Failed to load iAtlas cohort from {cohort_dir}: {exc}")
