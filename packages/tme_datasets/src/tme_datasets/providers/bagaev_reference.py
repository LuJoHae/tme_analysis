"""Pure loaders and ID-reconciliation helpers for BostonGene/Bagaev et al. (2021) reference tables."""

from __future__ import annotations

from pathlib import Path
from typing import Final
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

from ..download.fetcher import download_single_file
from ..genesets.models import GeneSet, GeneSetCollection
from ..genesets.parser import parse_gmt
from ..paths import get_data_paths

VALID_BAGAEV_SUBTYPES: Final[tuple[str, ...]] = ("IE", "IE/F", "F", "D")

BAGAEV_REFERENCE_URLS: Final[dict[str, dict[str, str]]] = {
    "Bagaev_Signature": {
        "gene_signatures.gmt": "https://raw.githubusercontent.com/BostonGene/MFP/refs/heads/master/signatures/gene_signatures.gmt",
        "gene_signatures_order.tsv": "https://raw.githubusercontent.com/BostonGene/MFP/refs/heads/master/signatures/gene_signatures_order.tsv",
    },
    "Bagaev_PanMelanoma": {
        "annotations.tsv": "https://science.bostongene.com/tumor-portrait/api/annotation_panmi/",
        "signature.tsv": "https://science.bostongene.com/tumor-portrait/api/signatures_panmi/",
    },
    "Bagaev_TCGA": {
        "annotations.tsv": "https://science.bostongene.com/tumor-portrait/api/annotation_tcga/",
        "signature.tsv": "https://science.bostongene.com/tumor-portrait/api/signatures_tcga/",
    },
}

_MELANOMA_COHORT_MAP: Final[dict[str, str]] = {
    "Hugo-iAtlas": "Hugo",
    "Riaz-iAtlas": "Riaz",
    "Liu-iAtlas": "Liu",
    "Gide-iAtlas": "Gide",
}


def normalize_sample_key(raw_id: str) -> str:
    """Normalize sample or patient identifiers for cross-table matching without delimiter sensitivity."""
    s = raw_id.strip()
    if "__" in s:
        s = s.split("__", 1)[1]
    return s.upper().replace("-", "_").replace(".", "_")


def _normalize_timepoint_token(raw_tp: str) -> str:
    """Normalize timepoint tokens into {'PRE', 'ON', 'UNK'} for biopsy-level matching."""
    t = raw_tp.strip().upper()
    if t in {"PRE", "BASELINE", "NAIVE", "PRE_TX", "PRE-TX"} or t.startswith("PRE"):
        return "PRE"
    if t in {"ON", "EDT", "POST", "ON_TREATMENT", "POST_TX", "POST-TX", "TREATED"} or t.startswith(("ON", "POST", "EDT")):
        return "ON"
    return "UNK"


def _resolve_raw_root(raw_root: Maybe[Path]) -> Path:
    return raw_root.unwrap() if isinstance(raw_root, Some) else get_data_paths().raw_dir


def ensure_bagaev_references(raw_root: Maybe[Path] = Nothing) -> Result[Path, str]:
    """Download Bagaev TCGA, PanMelanoma, and 29-Fges GMT files into `raw_dir/Bagaev_*`."""
    base_raw = _resolve_raw_root(raw_root)
    for subdir_name, files in BAGAEV_REFERENCE_URLS.items():
        target_dir = base_raw / subdir_name
        target_dir.mkdir(parents=True, exist_ok=True)
        for fname, url in files.items():
            dest = target_dir / fname
            if dest.is_file() and dest.stat().st_size > 0:
                continue
            dl_res = download_single_file(url, dest, resume=False)
            if isinstance(dl_res, Failure):
                return Failure(f"Failed downloading {subdir_name}/{fname}: {dl_res.failure()}")
    return Success(base_raw)


def _ensure_subdir(subdir_name: str, raw_root: Maybe[Path], auto_download: bool) -> Result[Path, str]:
    base_raw = _resolve_raw_root(raw_root)
    target_dir = base_raw / subdir_name
    req_files = BAGAEV_REFERENCE_URLS[subdir_name]
    missing = [fname for fname in req_files if not (target_dir / fname).is_file()]
    if missing and auto_download:
        target_dir.mkdir(parents=True, exist_ok=True)
        for fname in missing:
            dl_res = download_single_file(req_files[fname], target_dir / fname, resume=False)
            if isinstance(dl_res, Failure):
                return Failure(f"Failed downloading {subdir_name}/{fname}: {dl_res.failure()}")
    return Success(target_dir)


def load_bagaev_mfp_gmt_collection(
    raw_root: Maybe[Path] = Nothing,
    auto_download: bool = False,
) -> Result[GeneSetCollection, str]:
    """Load the 29-Fges GeneSetCollection from `raw/Bagaev_Signature/` inside `TME_DATA_DIR`."""
    dir_res = _ensure_subdir("Bagaev_Signature", raw_root, auto_download)
    if isinstance(dir_res, Failure):
        return dir_res
    sig_dir = dir_res.unwrap()
    gmt_path = sig_dir / "gene_signatures.gmt"
    if not gmt_path.is_file():
        return Failure(f"Missing Bagaev signature GMT at {gmt_path}.")

    gmt_res = parse_gmt(gmt_path)
    if isinstance(gmt_res, Failure):
        return Failure(gmt_res.failure())
    col = gmt_res.unwrap()

    order_path = sig_dir / "gene_signatures_order.tsv"
    if order_path.is_file():
        ordered_ids = [
            line.split("\t")[0].strip()
            for line in order_path.read_text(encoding="utf-8", errors="replace").splitlines()
            if line.strip() and line.split("\t")[0].strip() in col.gene_sets
        ]
        if ordered_ids:
            ordered_sets: dict[str, GeneSet] = {gid: col.gene_sets[gid] for gid in ordered_ids}
            for gid, gs in col.gene_sets.items():
                if gid not in ordered_sets:
                    ordered_sets[gid] = gs
            return Success(
                GeneSetCollection(
                    id="bagaev_mfp_29",
                    name="Bagaev MFP 29 Fges",
                    description="29 Functional Gene Expression Signatures from BostonGene/MFP",
                    gene_sets=ordered_sets,
                )
            )
    return Success(col)


def load_tcga_bagaev_reference(
    raw_root: Maybe[Path] = Nothing,
    auto_download: bool = False,
) -> Result[pl.DataFrame, str]:
    """Purely load BostonGene/Bagaev labeled TCGA 29-Fges reference profiles into Polars."""
    dir_res = _ensure_subdir("Bagaev_TCGA", raw_root, auto_download)
    if isinstance(dir_res, Failure):
        return dir_res
    tcga_dir = dir_res.unwrap()
    annot_path = tcga_dir / "annotations.tsv"
    sig_path = tcga_dir / "signature.tsv"
    if not annot_path.is_file() or not sig_path.is_file():
        return Failure(
            f"Missing Bagaev TCGA reference files in {tcga_dir}. Expected annotations.tsv and signature.tsv."
        )

    annot_lines = annot_path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(annot_lines) < 2:
        return Failure(f"Empty TCGA annotations file at {annot_path}.")

    header = ["sample_id", *annot_lines[0].split("\t")[1:]]
    annot_rows = [
        {header[i]: (parts[i].strip() if i < len(parts) else "") for i in range(len(header))}
        for line in annot_lines[1:]
        if line.strip() and (parts := line.split("\t"))
    ]
    req_annot_cols = ["sample_id", "TCGA_project", "HISTOLOGICAL_SUBTYPE", "MFP"]
    if not annot_rows or any(c not in header for c in req_annot_cols):
        return Failure(f"Unexpected schema in {annot_path}: found columns {header}.")
    annot_df = pl.DataFrame(annot_rows).select(req_annot_cols)

    sig_lines = sig_path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(sig_lines) < 2:
        return Failure(f"Empty TCGA signature matrix at {sig_path}.")

    sample_cols = [s.strip() for s in sig_lines[0].split("\t")[1:]]
    sig_matrix: dict[str, list[float]] = {
        parts[0].strip(): [
            float(x) if x.strip() not in ("", "NA", "nan") else 0.0
            for x in parts[1 : 1 + len(sample_cols)]
        ]
        for line in sig_lines[1:]
        if line.strip() and (parts := line.split("\t"))
    }
    if not sig_matrix:
        return Failure(f"No signature rows parsed from {sig_path}.")

    sig_df = pl.DataFrame({"sample_id": sample_cols, **sig_matrix})
    joined = annot_df.join(sig_df, on="sample_id", how="inner").filter(
        pl.col("MFP").is_in(list(VALID_BAGAEV_SUBTYPES))
    )
    if joined.is_empty():
        return Failure("Joined Bagaev TCGA reference is empty after subtype filtering.")
    return Success(joined)


def load_panmelanoma_bagaev_labels(
    raw_root: Maybe[Path] = Nothing,
    auto_download: bool = False,
) -> Result[pl.DataFrame, str]:
    """Purely load published Bagaev et al. Pan-Melanoma clinical trial MFP annotations."""
    dir_res = _ensure_subdir("Bagaev_PanMelanoma", raw_root, auto_download)
    if isinstance(dir_res, Failure):
        return dir_res
    pan_path = dir_res.unwrap() / "annotations.tsv"
    if not pan_path.is_file():
        return Failure(f"Pan-Melanoma annotations file not found at {pan_path}.")

    raw_df = pl.read_csv(pan_path, separator="\t", infer_schema_length=5000, ignore_errors=True)
    keep_cols = [c for c in ("Sample", "Cohort", "MFP", "Patient", "Response", "Time") if c in raw_df.columns]
    if "Sample" not in keep_cols or "Cohort" not in keep_cols or "MFP" not in keep_cols:
        return Failure(f"Missing required columns in {pan_path}: {raw_df.columns}.")
    return Success(raw_df.select(keep_cols))


def _derive_panmel_aliases(
    bg_cohort: str,
    samp_str: str,
    pat_str: str,
    tm_str: str,
) -> tuple[str, ...]:
    """Derive exact timepoint-preserving sample aliases for Bagaev Pan-Melanoma cohorts."""
    aliases: list[str] = [normalize_sample_key(samp_str)] if samp_str else []
    if bg_cohort == "Riaz" and pat_str.startswith("Riaz_") and tm_str:
        pt_core = pat_str.removeprefix("Riaz_")
        aliases.append(normalize_sample_key(f"{pt_core}_{tm_str}"))
    elif bg_cohort == "Liu" and samp_str.startswith("Patient"):
        num_core = samp_str.removeprefix("Patient")
        aliases.append(normalize_sample_key(f"Liu_Sample{num_core}"))
    elif bg_cohort == "Gide" and pat_str.startswith("Gide_") and tm_str:
        if pat_str.startswith("Gide_PD1_"):
            num_s = pat_str.removeprefix("Gide_PD1_")
            if num_s.isdigit():
                aliases.append(normalize_sample_key(f"PD{int(num_s):02d}_{tm_str}"))
        elif pat_str.startswith("Gide_ipiPD1_"):
            num_s = pat_str.removeprefix("Gide_ipiPD1_")
            if num_s.isdigit():
                aliases.append(normalize_sample_key(f"iP{int(num_s):02d}_{tm_str}"))
    return tuple(k for k in aliases if k)


def match_iatlas_to_panmelanoma_labels(
    iatlas_df: pl.DataFrame,
    panmel_df: pl.DataFrame,
) -> pl.DataFrame:
    """Attach published Pan-Melanoma `published_bagaev_mfp` labels preserving Pre/On biopsy specificity."""
    if iatlas_df.is_empty() or panmel_df.is_empty():
        return iatlas_df.with_columns(pl.lit("").alias("published_bagaev_mfp"))

    valid_pan = panmel_df.filter(pl.col("MFP").is_in(list(VALID_BAGAEV_SUBTYPES)))
    sample_lookup: dict[tuple[str, str], str] = {}
    patient_time_lookup: dict[tuple[str, str, str], str] = {}
    patient_mfp_sets: dict[tuple[str, str], set[str]] = {}

    for row in valid_pan.iter_rows(named=True):
        bg_coh = str(row.get("Cohort") or "").strip()
        mfp_lbl = str(row.get("MFP") or "").strip()
        samp_str = str(row.get("Sample") or "").strip()
        pat_str = str(row.get("Patient") or "").strip()
        tm_str = str(row.get("Time") or "").strip()
        tm_norm = _normalize_timepoint_token(tm_str)

        for alias in _derive_panmel_aliases(bg_coh, samp_str, pat_str, tm_str):
            sample_lookup.setdefault((bg_coh, alias), mfp_lbl)

        if pat_str:
            norm_pat = normalize_sample_key(pat_str)
            short_pat = normalize_sample_key(pat_str.removeprefix(f"{bg_coh}_"))
            for p_key in (norm_pat, short_pat):
                if p_key:
                    patient_mfp_sets.setdefault((bg_coh, p_key), set()).add(mfp_lbl)
                    if tm_norm != "UNK":
                        patient_time_lookup.setdefault((bg_coh, p_key, tm_norm), mfp_lbl)

    # Only allow timepoint-agnostic patient fallback when that patient has one unique label
    patient_unique_lookup: dict[tuple[str, str], str] = {
        k: next(iter(labels)) for k, labels in patient_mfp_sets.items() if len(labels) == 1
    }

    matched_labels: list[str] = []
    has_raw = "raw_sample_id" in iatlas_df.columns
    has_pat = "patient_id" in iatlas_df.columns
    has_trt = "sample_treatment" in iatlas_df.columns

    for row in iatlas_df.iter_rows(named=True):
        cid = str(row.get("cohort") or "")
        bg_coh = _MELANOMA_COHORT_MAP.get(cid, "")
        if not bg_coh:
            matched_labels.append("")
            continue

        raw_id = str(row.get("raw_sample_id") if has_raw else row.get("sample_id") or "")
        raw_k = normalize_sample_key(raw_id)
        pat_raw = str(row.get("patient_id") or "") if has_pat else ""
        pat_k = normalize_sample_key(pat_raw) if pat_raw else ""
        if not pat_k:
            pat_k = raw_k.split("_PRE")[0].split("_EDT")[0].split("_ON")[0]

        tm_raw = str(row.get("sample_treatment") or "") if has_trt else ""
        tm_norm = _normalize_timepoint_token(tm_raw)
        if tm_norm == "UNK":
            if "_PRE" in raw_k:
                tm_norm = "PRE"
            elif "_EDT" in raw_k or "_ON" in raw_k:
                tm_norm = "ON"

        hit = (
            sample_lookup.get((bg_coh, raw_k))
            or patient_time_lookup.get((bg_coh, pat_k, tm_norm))
            or patient_unique_lookup.get((bg_coh, pat_k))
            or ""
        )
        matched_labels.append(hit)

    return iatlas_df.with_columns(pl.Series("published_bagaev_mfp", matched_labels, dtype=pl.String))
