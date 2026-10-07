#!/usr/bin/env python3
"""Thin CLI orchestrator for preprocessing the 8 verified Tier 0 scRNA-seq ICB cohorts.

Standardizes raw cohort inputs into publication-grade AnnData (.h5ad) files containing:
- layers['counts']: Raw integer count matrix (sparse CSR float32).
- layers['tpm']: Linear TPM/CPM normalized matrix (10^6 per cell).
- layers['log1p']: Natural log-transformed normalized matrix log(1 + TPM).
- adata.X: Aligned with log1p for Scanpy / downstream tools.
- Standardized clinical observations in adata.obs (patient_id, clinical_response, etc.).
"""

from __future__ import annotations

import argparse
import gc
import os
import sys
import time
from pathlib import Path

import anndata as ad
from returns.result import Failure, Success

# Set scratch directory
os.environ.setdefault("TMPDIR", "/storage/halu/tmp")

from tme_datasets.logging import get_logger
from tme_datasets.models import QualityControlSpec, SingleCellProcessingSpec
from tme_datasets.preprocessing.pipeline import process_single_cell_dataset
from tme_datasets.query import load_dataset
from tme_datasets.registry import get_dataset_spec, resolve_dataset_id

logger = get_logger("preprocess_tier0")

# The 9 verified cohorts with 100% patient-traceable clinical response metadata
VERIFIED_TIER0_COHORTS: tuple[str, ...] = (
    "GSE120575",           # Melanoma (Sade-Feldman et al., Cell 2018)
    "CELLxGENE_7b20c613",  # Melanoma (Gondal et al., Sci Data 2025)
    "CELLxGENE_05a8c945",  # CRC (Marteau et al., Cancer Cell 2026)
    "CELLxGENE_6f9de485",  # TNBC (Navin et al., Nature 2026)
    "GSE207422",           # NSCLC (Hu et al., Genome Med 2023)
    "GSE243013",           # NSCLC (Liu et al., Cell 2025)
    "GSE233203",           # NSCLC (Hong et al., MedComm 2025)
    "GSE200996",           # HNSCC (Luoma et al., Cell 2022)
    "GSE316195",           # PDAC (Nature Comms 2026 snRNA-seq)
)



def preprocess_cohort(
    cohort_id: str,
    raw_dir: Path,
    output_dir: Path,
    overwrite: bool = False,
    dry_run: bool = False,
) -> dict[str, str | int]:
    """Preprocesses a single Tier 0 cohort, returning execution metrics."""
    canonical_id = resolve_dataset_id(cohort_id)
    cohort_raw_dir = raw_dir / canonical_id
    out_file = output_dir / f"{canonical_id}.h5ad"

    logger.info("=================================================================")
    logger.info("Processing Cohort: %s", canonical_id)
    logger.info("Raw Directory:     %s", cohort_raw_dir)
    logger.info("Target File:       %s", out_file)
    logger.info("=================================================================")

    if out_file.exists() and not overwrite and not dry_run:
        logger.info("Output file already exists: %s (skipping, use --overwrite to recompute)", out_file)
        adata_head = ad.read_h5ad(out_file, backed="r")
        n_cells = adata_head.n_obs
        n_vars = adata_head.n_vars
        n_pts = int(adata_head.obs["patient_id"].nunique()) if "patient_id" in adata_head.obs.columns else 0
        adata_head.file.close()
        del adata_head
        gc.collect()
        return {
            "cohort": canonical_id,
            "status": "ALREADY_EXISTS",
            "cells": n_cells,
            "genes": n_vars,
            "patients": n_pts,
            "file_size_mb": round(out_file.stat().st_size / (1024 * 1024), 2),
        }

    t0 = time.time()
    # Step 1: Load raw cohort via tme_datasets provider
    logger.info("[%s] Loading raw single-cell data...", canonical_id)
    load_res = load_dataset(
        canonical_id,
        raw_dir=cohort_raw_dir,
        auto_download=False,
        force_recompute=True,
        normalize_ensembl=False,
        cache_h5ad=False,
        apply_qc=False,
    )
    match load_res:
        case Failure(err):
            logger.error("[%s] Failed to load raw dataset: %s", canonical_id, err)
            return {"cohort": canonical_id, "status": f"FAILED: {err}", "cells": 0, "genes": 0, "patients": 0, "file_size_mb": 0}
        case Success(adata):
            pass

    raw_cells = adata.n_obs
    raw_genes = adata.n_vars
    n_patients = int(adata.obs["patient_id"].nunique()) if "patient_id" in adata.obs.columns else 0
    resp_counts = dict(adata.obs["clinical_response"].value_counts()) if "clinical_response" in adata.obs.columns else {}

    logger.info("[%s] Raw AnnData loaded: %d cells, %d genes, %d patients", canonical_id, raw_cells, raw_genes, n_patients)
    logger.info("[%s] Response breakdown: %s", canonical_id, resp_counts)

    if dry_run:
        logger.info("[%s] DRY-RUN completed successfully.", canonical_id)
        return {
            "cohort": canonical_id,
            "status": "DRY_RUN_OK",
            "cells": raw_cells,
            "genes": raw_genes,
            "patients": n_patients,
            "file_size_mb": 0,
        }

    # Step 2: Quality Control & Multi-layer Standardization
    # SnRNA-seq has tighter MT threshold (10%), Smart-seq2 has higher gene coverage and lower count floor
    is_sn = bool(adata.uns.get("is_single_nucleus", False)) or canonical_id == "GSE316195"
    is_smartseq = (
        bool(adata.uns.get("is_smartseq2", False))
        or canonical_id == "GSE120575"
        or adata.uns.get("expression_type") == "tpm"
    )
    qc_spec = QualityControlSpec(
        min_genes_per_cell=200,
        max_genes_per_cell=12000 if is_smartseq else 9000,
        min_counts_per_cell=100 if is_smartseq else 500,
        max_pct_mitochondrial=10.0 if is_sn else 20.0,
        min_cells_per_gene=3,
    )
    proc_spec = SingleCellProcessingSpec(
        qc=qc_spec,
        target_sum=1e6,
        use_adaptive_qc=False,
        generate_plots=False,
    )

    logger.info("[%s] Executing QC and standardization pipeline (MT max=%.1f%%)...", canonical_id, qc_spec.max_pct_mitochondrial)
    proc_res = process_single_cell_dataset(adata, spec=proc_spec, dataset_name=canonical_id)
    match proc_res:
        case Failure(err):
            logger.error("[%s] Preprocessing pipeline failed: %s", canonical_id, err)
            return {"cohort": canonical_id, "status": f"FAILED_QC: {err}", "cells": raw_cells, "genes": raw_genes, "patients": n_patients, "file_size_mb": 0}
        case Success(proc_result):
            processed_adata = proc_result.adata

    del adata
    gc.collect()

    # Step 3: Serialize to H5AD
    logger.info("[%s] Writing standardized .h5ad to: %s ...", canonical_id, out_file)
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_out = output_dir / f".tmp_{canonical_id}.h5ad"
    processed_adata.write_h5ad(temp_out, compression="gzip")
    temp_out.replace(out_file)

    elapsed = time.time() - t0
    final_size_mb = round(out_file.stat().st_size / (1024 * 1024), 2)
    post_cells = processed_adata.n_obs
    post_genes = processed_adata.n_vars
    logger.info("[%s] COMPLETED in %.1f seconds. Output size: %.1f MB (%d cells, %d genes)", canonical_id, elapsed, final_size_mb, post_cells, post_genes)

    del processed_adata
    gc.collect()

    return {
        "cohort": canonical_id,
        "status": "SUCCESS",
        "cells": raw_cells,
        "post_qc_cells": post_cells,
        "genes": post_genes,
        "patients": n_patients,
        "file_size_mb": final_size_mb,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Preprocess Tier 0 scRNA-seq ICB response cohorts into standardized .h5ad files.")
    parser.add_argument("--cohort", type=str, default="", help="Specific cohort to process (e.g. CELLxGENE_7b20c613, GSE207422).")
    parser.add_argument("--all", action="store_true", help="Process all 8 verified Tier 0 cohorts.")
    parser.add_argument("--raw-dir", type=Path, default=Path("/storage/halu/data-test/raw"), help="Path to raw dataset directory on disk.")
    parser.add_argument("--output-dir", type=Path, default=Path("/storage/halu/data-test/preprocessed"), help="Path to write preprocessed .h5ad files.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing preprocessed .h5ad files.")
    parser.add_argument("--dry-run", action="store_true", help="Inspect raw datasets and response columns without processing layers.")
    args = parser.parse_args()

    if not args.cohort and not args.all:
        parser.print_help()
        print("\nPlease specify either --cohort [ID] or --all.")
        sys.exit(1)

    cohorts_to_run = VERIFIED_TIER0_COHORTS if args.all else (args.cohort,)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    results: list[dict[str, str | int]] = []
    for c in cohorts_to_run:
        res = preprocess_cohort(
            cohort_id=c,
            raw_dir=args.raw_dir,
            output_dir=output_dir,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
        )
        results.append(res)

    print("\n" + "=" * 80)
    print("TIER 0 PREPROCESSING SUMMARY REPORT")
    print("=" * 80)
    print(f"{'Cohort':<22} {'Status':<16} {'Cells':<10} {'Genes':<10} {'Patients':<10} {'Size (MB)':<10}")
    print("-" * 80)
    for r in results:
        print(f"{r['cohort']:<22} {str(r['status'])[:15]:<16} {str(r.get('cells', 0)):<10} {str(r.get('genes', 0)):<10} {str(r.get('patients', 0)):<10} {str(r.get('file_size_mb', 0)):<10}")
    print("=" * 80)


if __name__ == "__main__":
    main()
