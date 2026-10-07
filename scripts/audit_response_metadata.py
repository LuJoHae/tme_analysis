#!/usr/bin/env python3
"""Comprehensive Audit Engine for scRNA-seq Immunotherapy Clinical Response Labels.

Audits downloaded datasets across Tier 0 and Tier 1 on disk to empirically verify:
1. Does the dataset actually contain clinical response labels (RECIST CR/PR/SD/PD, R/NR, pCR/MPR)?
2. Where are the labels located (Cell-level AnnData .obs, supplementary Excel/TSV, or GEO Series Matrix)?
3. What are the exact column names, patient ID keys, and unique response categories?
4. Reclassifies cohorts lacking verified patient response into Tier 1 (ICB Treated) or Tier 2 (Baseline Atlas).
"""

from __future__ import annotations

import argparse
import gzip
import io
import os
import re
import sys
import tarfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import anndata as ad
import h5py
import pandas as pd
import polars as pl


@dataclass(frozen=True)
class AuditRecord:
    """Audit result for a single cohort."""

    accession: str
    indication: str
    original_tier: str
    passed_audit: bool
    veracity_class: str          # Class A, B, C, D, E
    response_source: str         # File name / source
    response_column: str         # Exact column name
    patient_id_column: str       # Exact patient/sample key
    response_categories: tuple[str, ...]
    endpoint_type: str           # RECIST, Pathologic, Survival, Setting_Only, None
    recommended_tier: str
    audit_notes: str


# Regex patterns for clinical response detection
RESPONSE_COL_PATTERNS = [
    r"recist",
    r"path.*resp",
    r"patholog.*resp",
    r"^response$",
    r"^responder$",
    r"^resp$",
    r"clinical.*resp",
    r"treatment.*resp",
    r"treatment.*response",
    r"overall.*resp",
    r"best.*resp",
    r"objective.*resp",
    r"response.*status",
    r"response.*group",
    r"radiological.*resp",
    r"^pcr$",
    r"^mpr$",
    r"pcr.*status",
    r"mpr.*status",
    r"benefit",
    r"clinical.*benefit",
    r"^outcome$",
    r"combined.*outcome",
    r"clinical.*outcome",
    r"donor.*outcome",
    r"^bor$",
    r"^irrc$",
    r"efficacy",
]

# Regex patterns for treatment setting / kinetics
TREATMENT_COL_PATTERNS = [
    r"treatment",
    r"therapy",
    r"timepoint",
    r"time.*point",
    r"treatment.*status",
    r"pre.*post",
    r"pre.*on",
    r"sample.*type",
    r"cohort",
    r"group",
    r"arm",
]

# Known response values / categories
VALID_RESPONSE_TOKENS = {
    "cr", "pr", "sd", "pd",
    "complete response", "partial response", "stable disease", "progressive disease",
    "responder", "non-responder", "non responder", "nr", "r",
    "response", "non-response",
    "pcr", "npcr", "non-pcr", "mpr", "nmpr", "non-mpr", "rd", "residual disease",
    "durable benefit", "non-durable benefit", "dcb", "ndb", "clinical benefit",
    "favourable", "unfavourable", "favorable", "unfavorable",
    "high", "medium", "low",
    "post-ici (resistant)", "ici_pr", "ici_sd", "or",
}


def search_response_in_dataframe(
    df: pd.DataFrame, source_name: str
) -> tuple[str, str, tuple[str, ...], str] | None:
    """Inspect a dataframe for response columns and return (col, patient_col, unique_vals, endpoint_type)."""
    # Find prospective patient/sample identifier column
    patient_col = ""
    for col in df.columns:
        c_lower = str(col).lower()
        if any(k in c_lower for k in ["patient", "donor", "subject", "sample", "orig.ident", "case", "donor_id"]):
            patient_col = str(col)
            break

    # 1. Search for explicit RECIST or Pathologic Response columns
    for col in df.columns:
        c_lower = str(col).strip().lower()
        if any(re.search(pat, c_lower) for pat in RESPONSE_COL_PATTERNS):
            vals = [str(v).strip() for v in df[col].dropna().unique() if str(v).strip().lower() not in {"nan", "na", "", "none"}]
            if not vals:
                continue

            vals_lower = {v.lower() for v in vals}
            # Check if values contain recognizable clinical response labels
            is_recist = any(v in VALID_RESPONSE_TOKENS for v in vals_lower) or any(
                k in " ".join(vals_lower) for k in ["responder", "response", "benefit", "progression", "stable", "favourable", "favorable"]
            )
            is_pathologic = any(v in {"pcr", "mpr", "npcr", "nmpr", "rd", "high", "medium", "low"} for v in vals_lower) and any(
                k in c_lower for k in ["path", "pcr", "mpr"]
            )

            if is_recist:
                return (str(col), patient_col, tuple(sorted(vals[:15])), "RECIST")
            elif is_pathologic:
                return (str(col), patient_col, tuple(sorted(vals[:15])), "Pathologic")
            elif len(vals) <= 6 and any(v in VALID_RESPONSE_TOKENS for v in vals_lower):
                # Potential discrete clinical categorization
                return (str(col), patient_col, tuple(sorted(vals[:15])), "Clinical_Categorical")

    # 2. Check for survival / PFS / OS time and event columns
    pfs_col = next((c for c in df.columns if any(k in str(c).lower() for k in ["pfs", "progression_free", "rfs", "recurrence_free"])), None)
    if pfs_col is not None:
        vals = [str(v).strip() for v in df[pfs_col].dropna().unique()][:10]
        return (str(pfs_col), patient_col, tuple(vals), "Survival_Only")

    return None

    return None


def inspect_h5ad_file(h5ad_path: Path) -> tuple[str, str, tuple[str, ...], str] | None:
    """Inspect AnnData .obs table for response columns using direct h5py to avoid loading gigabytes."""
    try:
        with h5py.File(h5ad_path, "r") as f:
            if "obs" not in f:
                return None
            obs = f["obs"]
            obs_keys = list(obs.keys())

            # Find prospective patient/donor column
            patient_col = ""
            for k in obs_keys:
                k_lower = str(k).lower()
                if any(x in k_lower for x in ["donor_id", "patient", "donor", "subject", "sample"]):
                    patient_col = str(k)
                    break

            # Search response columns
            for col in obs_keys:
                c_lower = str(col).strip().lower()
                if any(re.search(pat, c_lower) for pat in RESPONSE_COL_PATTERNS):
                    target = obs[col]
                    vals: list[str] = []
                    if isinstance(target, h5py.Group) and "categories" in target:
                        raw_cats = target["categories"][:]
                        vals = [
                            c.decode("utf-8", errors="replace").strip() if isinstance(c, bytes) else str(c).strip()
                            for c in raw_cats
                        ]
                    elif isinstance(target, h5py.Dataset):
                        sample = target[:5000]
                        vals = list({
                            s.decode("utf-8", errors="replace").strip() if isinstance(s, bytes) else str(s).strip()
                            for s in sample
                        })

                    vals = [v for v in vals if v.lower() not in {"nan", "na", "", "none"}]
                    if not vals:
                        continue

                    vals_lower = {v.lower() for v in vals}
                    is_recist = any(v in VALID_RESPONSE_TOKENS for v in vals_lower) or any(
                        k in " ".join(vals_lower) for k in ["responder", "response", "benefit", "progression", "stable", "favourable", "favorable"]
                    )
                    is_pathologic = any(v in {"pcr", "mpr", "npcr", "nmpr", "rd", "high", "medium", "low"} for v in vals_lower) and any(
                        k in c_lower for k in ["path", "pcr", "mpr"]
                    )
                    if is_recist:
                        return (str(col), patient_col, tuple(sorted(vals[:15])), "RECIST")
                    elif is_pathologic:
                        return (str(col), patient_col, tuple(sorted(vals[:15])), "Pathologic")
                    elif len(vals) <= 6 and any(v in VALID_RESPONSE_TOKENS for v in vals_lower):
                        return (str(col), patient_col, tuple(sorted(vals[:15])), "Clinical_Categorical")
    except Exception as exc:
        print(f"    [WARN] Failed to inspect H5AD {h5ad_path.name}: {exc}", file=sys.stderr)
        return None
    return None


def inspect_excel_file(xlsx_path: Path) -> tuple[str, str, tuple[str, ...], str, str] | None:
    """Inspect all sheets in an Excel clinical metadata file."""
    try:
        xl = pd.ExcelFile(xlsx_path)
        for sheet in xl.sheet_names:
            df = xl.parse(sheet, nrows=500)
            res = search_response_in_dataframe(df, source_name=f"{xlsx_path.name} [{sheet}]")
            if res is not None:
                col, pat_col, vals, endpoint = res
                return (col, pat_col, vals, endpoint, f"{xlsx_path.name} (sheet: {sheet})")
    except Exception as exc:
        print(f"    [WARN] Failed to read Excel {xlsx_path.name}: {exc}", file=sys.stderr)
    return None


def inspect_table_file(table_path: Path) -> tuple[str, str, tuple[str, ...], str] | None:
    """Inspect TSV, CSV, or TXT metadata tables."""
    try:
        sep = "\t" if (table_path.name.endswith(".tsv") or table_path.name.endswith(".tsv.gz") or "txt" in table_path.name) else ","
        df = pd.read_csv(
            table_path,
            sep=sep,
            nrows=500,
            on_bad_lines="skip",
            low_memory=False,
            encoding="utf8",
        )
        return search_response_in_dataframe(df, source_name=table_path.name)
    except Exception:
        # Try alternate separator
        try:
            alt_sep = "," if sep == "\t" else "\t"
            df = pd.read_csv(table_path, sep=alt_sep, nrows=500, on_bad_lines="skip", low_memory=False)
            return search_response_in_dataframe(df, source_name=table_path.name)
        except Exception:
            return None


def fetch_geo_series_matrix_characteristics(
    gse_id: str,
) -> tuple[str, str, tuple[str, ...], str] | None:
    """Query NCBI GEO Series Matrix to inspect Sample characteristics (characteristics_ch1)."""
    prefix = gse_id[:-3]
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}nnn/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Python tme_datasets)"})
        with urllib.request.urlopen(req, timeout=20) as resp:
            compressed = resp.read()

        text = gzip.decompress(compressed).decode("utf-8", errors="replace")
        char_lines = [line.strip() for line in text.splitlines() if line.startswith("!Sample_characteristics_ch1")]
        if not char_lines:
            return None

        # Look for response keys across sample characteristics lines
        for line in char_lines:
            tokens = [t.strip().strip('"') for t in line.split("\t")[1:]]
            line_str = " ".join(tokens).lower()
            if any(k in line_str for k in ["response:", "responder:", "recist:", "pfs:", "outcome:", "benefit:"]):
                # Found sample-level response annotation
                # Extract unique values
                parsed_vals: set[str] = set()
                col_name = "characteristics_ch1"
                for tok in tokens:
                    if ":" in tok:
                        k, v = tok.split(":", 1)
                        if any(rk in k.lower() for rk in ["resp", "recist", "pfs", "outcome", "benefit"]):
                            col_name = k.strip()
                            parsed_vals.add(v.strip())
                    else:
                        parsed_vals.add(tok)

                vals_tuple = tuple(sorted(list(parsed_vals)[:15]))
                return (col_name, "!Sample_geo_accession", vals_tuple, "GEO_Characteristics")
    except Exception:
        return None
    return None


def audit_cohort_directory(
    accession: str,
    indication: str,
    original_tier: str,
    cohort_dir: Path,
    use_online_fallback: bool = True,
) -> AuditRecord:
    """Audit all files within a downloaded cohort directory."""
    if not cohort_dir.exists():
        return AuditRecord(
            accession=accession,
            indication=indication,
            original_tier=original_tier,
            passed_audit=False,
            veracity_class="Class D",
            response_source="Missing Directory",
            response_column="",
            patient_id_column="",
            response_categories=(),
            endpoint_type="None",
            recommended_tier="Tier 2 (Baseline Atlas)",
            audit_notes="Cohort directory does not exist on disk.",
        )

    all_files = list(cohort_dir.glob("*"))

    # 1. Check for native H5AD files first (Direct Cell-Level: Class A)
    h5ad_files = [f for f in all_files if f.name.endswith(".h5ad")]
    for h5_f in h5ad_files:
        res = inspect_h5ad_file(h5_f)
        if res is not None:
            col, pat_col, vals, endpoint = res
            if endpoint in {"RECIST", "Pathologic", "Clinical_Categorical"}:
                return AuditRecord(
                    accession=accession,
                    indication=indication,
                    original_tier=original_tier,
                    passed_audit=True,
                    veracity_class="Class A",
                    response_source=h5_f.name,
                    response_column=col,
                    patient_id_column=pat_col,
                    response_categories=vals,
                    endpoint_type=endpoint,
                    recommended_tier=original_tier if "tier 0" in original_tier.lower() else "Tier 1 (ICB Response)",
                    audit_notes=f"Verified directly in cell AnnData .obs table ({col}: {vals})",
                )

    # 2. Check for Excel clinical metadata supplements (Class B)
    excel_files = [f for f in all_files if f.name.endswith((".xlsx", ".xls")) and not f.name.startswith("~$")]
    for xl_f in excel_files:
        res = inspect_excel_file(xl_f)
        if res is not None:
            col, pat_col, vals, endpoint, src_desc = res
            if endpoint in {"RECIST", "Pathologic", "Clinical_Categorical"}:
                return AuditRecord(
                    accession=accession,
                    indication=indication,
                    original_tier=original_tier,
                    passed_audit=True,
                    veracity_class="Class B",
                    response_source=src_desc,
                    response_column=col,
                    patient_id_column=pat_col,
                    response_categories=vals,
                    endpoint_type=endpoint,
                    recommended_tier=original_tier if "tier 0" in original_tier.lower() else "Tier 1 (ICB Response)",
                    audit_notes=f"Verified in clinical Excel supplement ({src_desc} -> col: {col}, values: {vals})",
                )

    # 3. Check for metadata TSV/CSV/TXT tables (Class A or Class B)
    meta_tables = [
        f for f in all_files
        if any(f.name.lower().endswith(ext) for ext in [".tsv.gz", ".csv.gz", ".txt.gz", ".tsv", ".csv", ".txt"])
        and not f.name.startswith("filelist")
        and not any(m in f.name.lower() for m in ["counts", "matrix", "features", "genes", "barcodes"])
    ]

    for tbl_f in meta_tables:
        res = inspect_table_file(tbl_f)
        if res is not None:
            col, pat_col, vals, endpoint = res
            if endpoint in {"RECIST", "Pathologic", "Clinical_Categorical"}:
                return AuditRecord(
                    accession=accession,
                    indication=indication,
                    original_tier=original_tier,
                    passed_audit=True,
                    veracity_class="Class B",
                    response_source=tbl_f.name,
                    response_column=col,
                    patient_id_column=pat_col,
                    response_categories=vals,
                    endpoint_type=endpoint,
                    recommended_tier=original_tier if "tier 0" in original_tier.lower() else "Tier 1 (ICB Response)",
                    audit_notes=f"Verified in metadata table ({tbl_f.name} -> col: {col}, values: {vals})",
                )

    # 4. Check for TAR archives: peek inside without full extraction
    tar_files = [f for f in all_files if f.name.endswith((".tar", ".tar.gz"))]
    for tar_f in tar_files:
        try:
            with tarfile.open(tar_f, "r:*") as tar:
                members = tar.getmembers()
                # Find metadata candidates inside tar
                meta_members = [
                    m for m in members
                    if any(k in m.name.lower() for k in ["meta", "clinical", "sample", "patient", "annot"])
                    and not any(k in m.name.lower() for k in ["matrix", "counts", "gene", "feature"])
                    and m.size < 50 * 1024 * 1024  # Less than 50MB
                ]
                for m in meta_members:
                    f_obj = tar.extractfile(m)
                    if f_obj is not None:
                        try:
                            # Read bytes into buffer
                            content = f_obj.read()
                            if m.name.endswith(".gz"):
                                content = gzip.decompress(content)
                            df = pd.read_csv(io.BytesIO(content), sep=None, engine="python", nrows=500)
                            res = search_response_in_dataframe(df, source_name=f"{tar_f.name} [{m.name}]")
                            if res is not None:
                                col, pat_col, vals, endpoint = res
                                if endpoint in {"RECIST", "Pathologic", "Clinical_Categorical"}:
                                    return AuditRecord(
                                        accession=accession,
                                        indication=indication,
                                        original_tier=original_tier,
                                        passed_audit=True,
                                        veracity_class="Class B",
                                        response_source=f"{tar_f.name} -> {m.name}",
                                        response_column=col,
                                        patient_id_column=pat_col,
                                        response_categories=vals,
                                        endpoint_type=endpoint,
                                        recommended_tier=original_tier if "tier 0" in original_tier.lower() else "Tier 1 (ICB Response)",
                                        audit_notes=f"Verified inside TAR archive ({tar_f.name}/{m.name} -> {col}: {vals})",
                                    )
                        except Exception:
                            continue
        except Exception:
            pass

    # 5. Online Fallback: Query NCBI GEO Series Matrix Characteristics
    if use_online_fallback and accession.startswith("GSE"):
        print(f"    [INFO] Local metadata lacked explicit response. Querying NCBI GEO Series Matrix for {accession}...", flush=True)
        geo_res = fetch_geo_series_matrix_characteristics(accession)
        if geo_res is not None:
            col, pat_col, vals, endpoint = geo_res
            vals_lower = {v.lower() for v in vals}
            if any(v in VALID_RESPONSE_TOKENS for v in vals_lower) or any(
                k in " ".join(vals_lower) for k in ["responder", "response", "benefit", "progression", "stable", "favourable", "favorable"]
            ):
                return AuditRecord(
                    accession=accession,
                    indication=indication,
                    original_tier=original_tier,
                    passed_audit=True,
                    veracity_class="Class B",
                    response_source="GEO Series Matrix (Sample Characteristics)",
                    response_column=col,
                    patient_id_column=pat_col,
                    response_categories=vals,
                    endpoint_type="RECIST_GEO",
                    recommended_tier=original_tier if "tier 0" in original_tier.lower() else "Tier 1 (ICB Response)",
                    audit_notes=f"Verified in GEO Sample Characteristics ({col}: {vals})",
                )

    # 6. Check if it is a treated cohort with longitudinal kinetics (Class C)
    # Check if files indicate treatment timing
    is_treated = "treated" in original_tier.lower() or any(
        k in f.name.lower() for f in all_files for k in ["pembro", "nivo", "treated", "pre_post", "longitudinal"]
    )
    if is_treated:
        return AuditRecord(
            accession=accession,
            indication=indication,
            original_tier=original_tier,
            passed_audit=False,
            veracity_class="Class C",
            response_source="Files present, no outcome breakdown",
            response_column="",
            patient_id_column="",
            response_categories=(),
            endpoint_type="Setting_Only",
            recommended_tier="Tier 1 (ICB Treated)",
            audit_notes="Cohort was treated with immunotherapy, but patient-level response outcome is not deposited.",
        )

    # 7. Default: Unverified / In vitro / Missing (Class D)
    return AuditRecord(
        accession=accession,
        indication=indication,
        original_tier=original_tier,
        passed_audit=False,
        veracity_class="Class D",
        response_source="None",
        response_column="",
        patient_id_column="",
        response_categories=(),
        endpoint_type="None",
        recommended_tier="Tier 2 (Baseline Atlas)",
        audit_notes="No patient-level clinical response labels found in public files; reclassified to Baseline Atlas.",
    )


def run_full_audit(
    registry_path: Path,
    raw_dir: Path,
    use_online_fallback: bool = True,
) -> list[AuditRecord]:
    """Run audit across all Tier 0, Tier 1, and downloaded cohorts on disk."""
    df = pl.read_parquet(registry_path)

    existing_dirs = {p.name for p in raw_dir.iterdir() if p.is_dir()} if raw_dir.exists() else set()

    # Target all cohorts in Tier 0, Tier 1, or downloaded on disk
    target_df = df.filter(
        pl.col("tier").str.starts_with("Tier 0")
        | pl.col("tier").str.starts_with("Tier 1")
        | pl.col("accession").is_in(list(existing_dirs))
    )

    records: list[AuditRecord] = []
    total = len(target_df)

    print(f"\n=======================================================", flush=True)
    print(f"  Starting Clinical Response Label Audit ({total} cohorts)", flush=True)
    print(f"  Raw Directory: {raw_dir}", flush=True)
    print(f"  Online GEO Matrix Fallback: {use_online_fallback}", flush=True)
    print(f"=======================================================\n", flush=True)

    for idx, row in enumerate(target_df.iter_rows(named=True), 1):
        acc = str(row["accession"])
        ind = str(row["indication"])
        tier = str(row["tier"])
        c_dir = raw_dir / acc

        rec = audit_cohort_directory(
            accession=acc,
            indication=ind,
            original_tier=tier,
            cohort_dir=c_dir,
            use_online_fallback=use_online_fallback,
        )
        records.append(rec)

        symbol = "✓" if rec.passed_audit else ("⚠" if rec.veracity_class == "Class C" else "✗")
        col_str = f"[{rec.response_column}]" if rec.response_column else "[-]"
        print(f"[{idx:3d}/{total:3d}] {symbol} {rec.accession:<18} | {rec.veracity_class:<7} | {rec.endpoint_type:<16} | {col_str:<25} | {rec.audit_notes[:50]}...", flush=True)

    return records


def generate_audit_report(records: list[AuditRecord], output_path: Path) -> None:
    """Generate comprehensive publication-grade markdown audit report."""
    md: list[str] = []
    md.append("# Immunotherapy Response Metadata Audit Report\n")
    md.append(
        f"This document provides the empirical verification results for all **{len(records)} single-cell cohorts** "
        "(Tier 0 Benchmark Core, Tier 1 ICB Response / Treated, and on-disk candidates) inspected on server `olm`.\n"
    )

    total = len(records)
    verified_cnt = sum(1 for r in records if r.passed_audit)
    class_a_cnt = sum(1 for r in records if r.veracity_class == "Class A")
    class_b_cnt = sum(1 for r in records if r.veracity_class == "Class B")
    class_c_cnt = sum(1 for r in records if r.veracity_class == "Class C")
    class_d_cnt = sum(1 for r in records if r.veracity_class == "Class D")

    md.append("## Executive Summary Statistics")
    md.append(f"- **Total Cohorts Audited**: {total}")
    md.append(f"- **Verified Response Ground Truth (Class A + B)**: **{verified_cnt} cohorts ({verified_cnt/total*100:.1f}%)**")
    md.append(f"  - **Class A (Direct Cell-Level Response)**: {class_a_cnt} cohorts")
    md.append(f"  - **Class B (Patient-Level Table / GEO Matrix)**: {class_b_cnt} cohorts")
    md.append(f"- **Class C (ICB Treated Only, No Response Breakdown)**: {class_c_cnt} cohorts")
    md.append(f"- **Class D (False Positive / Pre-Clinical / In Vitro)**: {class_d_cnt} cohorts\n")

    md.append("## Master Audit Table\n")
    md.append("| Accession | Indication | Original Tier | Veracity Class | Verified? | Response Column | Patient ID Key | Unique Response Categories | Source File |")
    md.append("| :--- | :--- | :--- | :--- | :---: | :--- | :--- | :--- | :--- |")

    for r in records:
        v_mark = "✓ YES" if r.passed_audit else "✗ NO"
        cat_str = ", ".join(r.response_categories[:6]) if r.response_categories else "None"
        md.append(
            f"| **{r.accession}** | {r.indication} | {r.original_tier} | **{r.veracity_class}** | {v_mark} | "
            f"`{r.response_column}` | `{r.patient_id_column}` | `{cat_str}` | {r.response_source} |"
        )

    md.append("\n---\n")
    md.append("## Detailed Cohort Profiles & Notes\n")

    for idx, r in enumerate(records, 1):
        status_str = "VERIFIED RESPONSE GROUND TRUTH" if r.passed_audit else "UNVERIFIED / SETTING ONLY"
        md.append(f"### {idx}. {r.accession} — {r.indication} ({r.veracity_class})")
        md.append(f"**Verification Status**: `{status_str}`\n")
        md.append(f"- **Original Tier**: {r.original_tier}")
        md.append(f"- **Recommended Tier**: **{r.recommended_tier}**")
        md.append(f"- **Veracity Class**: `{r.veracity_class}`")
        md.append(f"- **Endpoint Type**: `{r.endpoint_type}`")
        md.append(f"- **Response Column**: `{r.response_column}`")
        md.append(f"- **Patient / Sample Key**: `{r.patient_id_column}`")
        md.append(f"- **Response Categories Found**: `{list(r.response_categories)}`")
        md.append(f"- **Metadata Source**: `{r.response_source}`")
        md.append(f"- **Audit Notes**: {r.audit_notes}\n")

    output_path.write_text("\n".join(md))
    print(f"\nWritten audit report to {output_path}")


def apply_registry_updates(
    registry_path: Path,
    tsv_path: Path,
    records: list[AuditRecord],
) -> None:
    """Update Parquet and TSV registries with audited response flags and reclassifications."""
    df = pl.read_parquet(registry_path)
    rec_dict = {r.accession: r for r in records}

    new_tiers: list[str] = []
    annotated_flags: list[bool] = []
    verified_flags: list[bool] = []
    veracity_classes: list[str] = []
    response_cols: list[str] = []
    patient_cols: list[str] = []
    response_cats: list[str] = []
    audit_notes: list[str] = []

    for row in df.iter_rows(named=True):
        acc = str(row["accession"])
        curr_tier = str(row["tier"])
        curr_annot = bool(row["clinical_response_annotated"])

        rec = rec_dict.get(acc)
        if rec:
            final_tier = "Tier 0 (Benchmark Core)" if rec.passed_audit else rec.recommended_tier
            new_tiers.append(final_tier)
            annotated_flags.append(rec.passed_audit)
            verified_flags.append(rec.passed_audit)
            veracity_classes.append(rec.veracity_class)
            response_cols.append(rec.response_column)
            patient_cols.append(rec.patient_id_column)
            response_cats.append("; ".join(rec.response_categories))
            audit_notes.append(rec.audit_notes)
        else:
            new_tiers.append("Tier 2 (Baseline Atlas)")
            annotated_flags.append(False)
            verified_flags.append(False)
            veracity_classes.append("Class D (Untreated Atlas)")
            response_cols.append("")
            patient_cols.append("")
            response_cats.append("")
            audit_notes.append("Baseline tumor atlas; not in clinical immunotherapy pool.")

    updated_df = df.with_columns(
        pl.Series("tier", new_tiers),
        pl.Series("clinical_response_annotated", annotated_flags),
        pl.Series("audited_response_verified", verified_flags),
        pl.Series("response_veracity_class", veracity_classes),
        pl.Series("response_column_name", response_cols),
        pl.Series("patient_id_column", patient_cols),
        pl.Series("response_categories", response_cats),
        pl.Series("audit_notes", audit_notes),
    )

    updated_df.write_parquet(registry_path)
    updated_df.write_csv(tsv_path, separator="\t")
    print(f"Successfully updated registry: {registry_path} and {tsv_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit and verify clinical response labels across scRNA-seq cohorts.")
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=Path("/storage/halu/data-test/raw"),
        help="Path to downloaded raw cohort directories.",
    )
    parser.add_argument(
        "--registry",
        type=Path,
        default=Path("data/registry/discovered_solid_tumor_sc_datasets.parquet"),
        help="Path to Parquet registry.",
    )
    parser.add_argument(
        "--tsv",
        type=Path,
        default=Path("data/registry/discovered_solid_tumor_sc_datasets.tsv"),
        help="Path to TSV registry.",
    )
    parser.add_argument(
        "--output-report",
        type=Path,
        default=Path("docs/immunotherapy_response_metadata_audit.md"),
        help="Path to markdown output report.",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Apply audited reclassifications and columns to Parquet and TSV registries.",
    )
    parser.add_argument(
        "--no-online-fallback",
        action="store_true",
        help="Disable online GEO series matrix characteristics fetching.",
    )

    args = parser.parse_args()

    records = run_full_audit(
        registry_path=args.registry,
        raw_dir=args.raw_dir,
        use_online_fallback=not args.no_online_fallback,
    )

    generate_audit_report(records, args.output_report)

    if args.apply:
        apply_registry_updates(
            registry_path=args.registry,
            tsv_path=args.tsv,
            records=records,
        )


if __name__ == "__main__":
    main()
