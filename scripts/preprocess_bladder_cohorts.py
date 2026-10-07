#!/usr/bin/env python3
"""CLI orchestrator for preprocessing Bladder Cancer single-cell & single-nucleus cohorts.

Target Cohorts:
- GSE222315: Primary Bladder Cancer vs. Normal Adjacent Tissue (13 samples, paired MTX+TSV)
- GSE302781: Metastatic Urothelial Carcinoma Rapid Autopsy snRNA-seq (15 samples, 10x triplets)
- GSE326225: Muscle-Invasive Bladder Cancer scRNA-seq (8 samples, 10x triplets)
- GSE301651: Single-cell sequencing of bladder cancer patients (13 samples, 10x triplets)

Standardizes raw cohort inputs into publication-grade AnnData (.h5ad) files containing:
- layers['counts']: Raw integer count matrix (sparse CSR float32).
- layers['tpm']: Linear TPM/CPM normalized matrix (10^6 per cell).
- layers['log1p']: Natural log-transformed normalized matrix log(1 + TPM).
- layers['log1p_norm']: Semantic alias for log1p.
- adata.X: Aligned with log1p for Scanpy / downstream tools.
- adata.raw: Post-QC raw counts snapshot.
- Canonical Ensembl Release 111 gene identifiers with enriched .var attributes.
- Standardized clinical observations in adata.obs (patient_id, sample_id, tissue_type, etc.).
"""

from __future__ import annotations

import argparse
import gc
import os
from pathlib import Path
import sys
import time

import anndata as ad
from returns.result import Failure, Success

# Set scratch directory
os.environ.setdefault("TMPDIR", "/storage/halu/tmp")

from tme_datasets.logging import get_logger
from tme_datasets.models import QualityControlSpec, SingleCellProcessingSpec
from tme_datasets.preprocessing.gene_normalization import normalize_dataset_to_ensembl
from tme_datasets.preprocessing.pipeline import process_single_cell_dataset
from tme_datasets.query import load_dataset
from tme_datasets.registry import resolve_dataset_id

logger = get_logger("preprocess_bladder")

BLADDER_COHORTS: tuple[str, ...] = (
    "GSE222315",  # Primary Bladder Cancer vs Normal Adjacent Tissue (13 samples)
    "GSE302781",  # Metastatic Urothelial Carcinoma snRNA-seq (15 samples)
    "GSE326225",  # Muscle-Invasive Bladder Cancer scRNA-seq (8 samples)
    "GSE301651",  # Bladder Cancer scRNA-seq (13 samples, primary + LN met + PBMC)
)


def preprocess_bladder_cohort(
    cohort_id: str,
    raw_dir: Path,
    output_dir: Path,
    overwrite: bool = False,
    dry_run: bool = False,
    skip_ensembl: bool = False,
) -> dict[str, object]:
    """Preprocesses a single Bladder Cancer cohort into a standardized .h5ad."""
    canonical_id = resolve_dataset_id(cohort_id)
    cohort_raw_dir = raw_dir / canonical_id
    out_file = output_dir / f"{canonical_id}.h5ad"

    logger.info("=================================================================")
    logger.info("Processing Bladder Cohort: %s", canonical_id)
    logger.info("Raw Directory:             %s", cohort_raw_dir)
    logger.info("Target File:               %s", out_file)
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
    logger.info("[%s] Loading raw single-cell / single-nucleus data...", canonical_id)
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
            return {"cohort": canonical_id, "status": f"FAILED_LOAD: {err}", "cells": 0, "genes": 0, "patients": 0, "file_size_mb": 0}
        case Success(adata):
            pass

    raw_cells = adata.n_obs
    raw_genes = adata.n_vars
    n_patients = int(adata.obs["patient_id"].nunique()) if "patient_id" in adata.obs.columns else 0
    n_samples = int(adata.obs["sample_id"].nunique()) if "sample_id" in adata.obs.columns else 0

    logger.info(
        "[%s] Raw AnnData loaded: %d cells, %d genes, %d patients, %d samples",
        canonical_id,
        raw_cells,
        raw_genes,
        n_patients,
        n_samples,
    )

    if dry_run:
        logger.info("[%s] DRY-RUN completed successfully.", canonical_id)
        return {
            "cohort": canonical_id,
            "status": "DRY_RUN_OK",
            "cells": raw_cells,
            "genes": raw_genes,
            "patients": n_patients,
            "samples": n_samples,
            "file_size_mb": 0,
        }

    # Step 2: Ensembl Release 111 Normalization
    if not skip_ensembl:
        logger.info("[%s] Normalizing gene identifiers to Ensembl Release 111...", canonical_id)
        norm_res = normalize_dataset_to_ensembl(
            adata,
            release=111,
            species="human",
            drop_unmapped=True,
            aggregation="sum",
        )
        match norm_res:
            case Failure(err):
                logger.error("[%s] Ensembl normalization failed: %s", canonical_id, err)
                return {"cohort": canonical_id, "status": f"FAILED_ENSEMBL: {err}", "cells": raw_cells, "genes": raw_genes, "patients": n_patients, "file_size_mb": 0}
            case Success(norm_adata):
                adata = norm_adata
        logger.info("[%s] Ensembl normalization complete: %d canonical genes retained.", canonical_id, adata.n_vars)

    # Step 3: Quality Control & Multi-layer Standardization
    # GSE302781 is snRNA-seq (10% MT max); GSE222315 and GSE326225 are scRNA-seq (20% MT max)
    is_sn = canonical_id == "GSE302781" or bool(adata.uns.get("is_single_nucleus", False))
    qc_spec = QualityControlSpec(
        min_genes_per_cell=200,
        max_genes_per_cell=9000,
        min_counts_per_cell=500,
        max_pct_mitochondrial=10.0 if is_sn else 20.0,
        min_cells_per_gene=3,
    )
    proc_spec = SingleCellProcessingSpec(
        qc=qc_spec,
        target_sum=1e6,
        use_adaptive_qc=False,
        generate_plots=False,
    )

    logger.info("[%s] Executing QC pipeline (snRNA=%s, max_pct_mt=%.1f%%)...", canonical_id, is_sn, qc_spec.max_pct_mitochondrial)
    proc_res = process_single_cell_dataset(adata, spec=proc_spec, dataset_name=canonical_id)
    del adata
    gc.collect()

    match proc_res:
        case Failure(err):
            logger.error("[%s] Preprocessing pipeline failed: %s", canonical_id, err)
            return {"cohort": canonical_id, "status": f"FAILED_QC: {err}", "cells": raw_cells, "genes": raw_genes, "patients": n_patients, "file_size_mb": 0}
        case Success(proc_result):
            processed_adata = proc_result.adata

    # Ensure index names do not clash with columns during H5AD serialization
    if processed_adata.obs_names.name is not None and processed_adata.obs_names.name in processed_adata.obs.columns:
        processed_adata.obs_names.name = None
    if processed_adata.var_names.name is not None and processed_adata.var_names.name in processed_adata.var.columns:
        processed_adata.var_names.name = None

    # Step 4: Serialize to H5AD
    logger.info("[%s] Serializing standardized .h5ad to: %s ...", canonical_id, out_file)
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_out = output_dir / f".tmp_{canonical_id}.h5ad"
    processed_adata.write_h5ad(temp_out, compression="gzip")
    temp_out.replace(out_file)

    elapsed = time.time() - t0
    final_size_mb = round(out_file.stat().st_size / (1024 * 1024), 2)
    post_cells = processed_adata.n_obs
    post_genes = processed_adata.n_vars
    logger.info(
        "[%s] COMPLETED in %.1f seconds. Output size: %.1f MB (%d cells, %d genes)",
        canonical_id,
        elapsed,
        final_size_mb,
        post_cells,
        post_genes,
    )

    del processed_adata
    gc.collect()

    return {
        "cohort": canonical_id,
        "status": "SUCCESS",
        "cells": raw_cells,
        "post_qc_cells": post_cells,
        "genes": post_genes,
        "patients": n_patients,
        "samples": n_samples,
        "file_size_mb": final_size_mb,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Preprocess Bladder Cancer single-cell & snRNA-seq cohorts into standardized .h5ad files.")
    parser.add_argument("--cohort", type=str, default="", help="Specific cohort to process (e.g. GSE222315, GSE302781, GSE326225).")
    parser.add_argument("--all", action="store_true", help="Process all 3 bladder cancer cohorts.")
    parser.add_argument("--raw-dir", type=Path, default=Path("/storage/halu/data-test/raw"), help="Path to raw dataset directory on disk.")
    parser.add_argument("--output-dir", type=Path, default=Path("/storage/halu/data-test/preprocessed"), help="Path to write preprocessed .h5ad files.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing preprocessed .h5ad files.")
    parser.add_argument("--dry-run", action="store_true", help="Inspect raw datasets and metadata without saving H5AD.")
    parser.add_argument("--skip-ensembl", action="store_true", help="Skip Ensembl Release 111 re-indexing.")
    args = parser.parse_args()

    if not args.cohort and not args.all:
        parser.print_help()
        print("\nPlease specify either --cohort [ID] or --all.")
        sys.exit(1)

    cohorts_to_run = BLADDER_COHORTS if args.all else (args.cohort,)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    results: list[dict[str, object]] = []
    for c in cohorts_to_run:
        res = preprocess_bladder_cohort(
            cohort_id=c,
            raw_dir=args.raw_dir,
            output_dir=output_dir,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
            skip_ensembl=args.skip_ensembl,
        )
        results.append(res)

    print("\n" + "=" * 90)
    print("BLADDER CANCER PREPROCESSING SUMMARY REPORT")
    print("=" * 90)
    print(f"{'Cohort':<14} {'Status':<16} {'Raw Cells':<11} {'Post-QC':<11} {'Genes':<9} {'Pts':<5} {'Smps':<6} {'Size (MB)':<10}")
    print("-" * 90)
    for r in results:
        print(
            f"{str(r['cohort']):<14} "
            f"{str(r['status'])[:15]:<16} "
            f"{str(r.get('cells', 0)):<11} "
            f"{str(r.get('post_qc_cells', 0)):<11} "
            f"{str(r.get('genes', 0)):<9} "
            f"{str(r.get('patients', 0)):<5} "
            f"{str(r.get('samples', 0)):<6} "
            f"{str(r.get('file_size_mb', 0)):<10}"
        )
    print("=" * 90)


if __name__ == "__main__":
    main()
