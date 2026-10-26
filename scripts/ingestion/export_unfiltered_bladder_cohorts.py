#!/usr/bin/env python3
"""Export unfiltered Bladder Cancer single-cell & single-nucleus AnnData (.h5ad) files.

Retains ALL cells (including low-quality droplets / failing cells) and annotates them
with explicit boolean quality control flags and covariates:
- obs['is_retained_qc']: True if cell passed all QC criteria, False otherwise.
- obs['qc_fail_low_genes']: True if detected genes < min_genes_per_cell (200).
- obs['qc_fail_high_genes']: True if detected genes > max_genes_per_cell (9000).
- obs['qc_fail_low_counts']: True if UMI counts < min_counts_per_cell (500).
- obs['qc_fail_high_mt']: True if mitochondrial % > threshold (10% for snRNA-seq, 20% for scRNA-seq).
- obs['total_counts']: Total UMI count per cell.
- obs['n_genes_by_counts']: Number of detected genes per cell.
- obs['pct_counts_mt']: Mitochondrial percentage per cell.
- obs['pct_counts_ribo']: Ribosomal percentage per cell.
- obs['pct_counts_hb']: Hemoglobin percentage per cell.
- var['is_retained_qc_gene']: True if gene expressed in >= 3 cells.
- uns['qc_spec']: The exact QC parameters applied.
- uns['sc_processing_summary']: Cell and gene counts before/after hypothetical filtering.

Target Cohorts:
- GSE326225 (Muscle-Invasive Bladder Cancer scRNA-seq)
- GSE302781 (Metastatic Urothelial Carcinoma Rapid Autopsy snRNA-seq)
- GSE222315 (Primary Bladder Cancer vs Normal Adjacent Tissue scRNA-seq)
- GSE301651 (Primary Bladder Cancer, Matched LN Met, and PBMC scRNA-seq)
"""

from __future__ import annotations

import argparse
import gc
import os
from pathlib import Path
import sys
import time
from typing import Mapping

import anndata as ad
import numpy as np
from pydantic import BaseModel, ConfigDict
from returns.result import Failure, Result, Success
import scipy.sparse as sp

# Set scratch directory to avoid /tmp out-of-space crashes during internal temp file operations
os.environ.setdefault("TMPDIR", "/ritom/scratch/cscs/lhaeuser/tmp")

from tme_datasets.logging import get_logger
from tme_datasets.models import QualityControlSpec
from tme_datasets.preprocessing.gene_normalization import normalize_dataset_to_ensembl
from tme_datasets.preprocessing.qc import compute_qc_covariates
from tme_datasets.query import load_dataset
from tme_datasets.registry import resolve_dataset_id

logger = get_logger("export_unfiltered_bladder")

DEFAULT_COHORTS: tuple[str, ...] = (
    "GSE326225",
    "GSE302781",
    "GSE222315",
    "GSE301651",
)


class ExportResult(BaseModel):
    model_config = ConfigDict(frozen=True)
    cohort: str
    total_cells: int
    passing_cells: int
    failing_cells: int
    pct_passing: float
    total_genes: int
    passing_genes: int
    output_path: str
    file_size_mb: float
    status: str


def export_unfiltered_cohort(
    cohort_id: str,
    raw_dir: Path,
    output_dir: Path,
    skip_ensembl: bool = False,
    overwrite: bool = False,
) -> Result[ExportResult, str]:
    """Loads raw cohort, computes QC flags for all cells, and saves unfiltered .h5ad."""
    canonical_id = resolve_dataset_id(cohort_id)
    cohort_raw_dir = raw_dir / canonical_id
    out_file = output_dir / f"{canonical_id}_unfiltered.h5ad"

    logger.info("=================================================================")
    logger.info("Exporting Unfiltered Cohort: %s", canonical_id)
    logger.info("Raw Directory:               %s", cohort_raw_dir)
    logger.info("Output Destination:          %s", out_file)
    logger.info("=================================================================")

    if out_file.exists() and not overwrite:
        logger.info("[%s] Output file already exists: %s (skipping)", canonical_id, out_file)
        adata_head = ad.read_h5ad(out_file, backed="r")
        n_total = adata_head.n_obs
        n_passing = int(adata_head.obs["is_retained_qc"].sum()) if "is_retained_qc" in adata_head.obs.columns else 0
        n_vars = adata_head.n_vars
        if hasattr(adata_head, "file") and adata_head.file is not None:
            adata_head.file.close()
        del adata_head
        gc.collect()

        return Success(
            ExportResult(
                cohort=canonical_id,
                total_cells=n_total,
                passing_cells=n_passing,
                failing_cells=n_total - n_passing,
                pct_passing=round((n_passing / max(n_total, 1)) * 100.0, 2),
                total_genes=n_vars,
                passing_genes=n_vars,
                output_path=str(out_file),
                file_size_mb=round(out_file.stat().st_size / (1024 * 1024), 2),
                status="ALREADY_EXISTS",
            )
        )

    t0 = time.time()

    # Step 1: Load raw data without QC or filtering
    logger.info("[%s] Loading raw single-cell/nucleus data...", canonical_id)
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
            return Failure(f"Failed to load raw dataset: {err}")
        case Success(raw_adata):
            adata = raw_adata

    raw_cells = adata.n_obs
    raw_genes = adata.n_vars
    logger.info("[%s] Raw AnnData loaded: %d cells x %d genes", canonical_id, raw_cells, raw_genes)

    # Step 2: Ensembl Release 111 normalization
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
                return Failure(f"Ensembl normalization failed: {err}")
            case Success(norm_adata):
                adata = norm_adata
        logger.info("[%s] Ensembl normalization complete: %d canonical genes retained.", canonical_id, adata.n_vars)

    # Step 3: Compute QC covariates on ALL cells
    logger.info("[%s] Computing QC covariates for all %d cells...", canonical_id, adata.n_obs)
    adata = compute_qc_covariates(adata)

    # Step 4: Apply QC criteria to establish boolean masks
    is_sn = canonical_id == "GSE302781" or bool(adata.uns.get("is_single_nucleus", False))
    max_pct_mt = 10.0 if is_sn else 20.0
    min_genes = 200
    max_genes = 9000
    min_counts = 500
    min_cells_per_gene = 3

    genes_detected = adata.obs["n_genes_by_counts"].to_numpy()
    total_counts = adata.obs["total_counts"].to_numpy()
    pct_mt = adata.obs["pct_counts_mt"].to_numpy()

    fail_low_genes = genes_detected < min_genes
    fail_high_genes = genes_detected > max_genes
    fail_low_counts = total_counts < min_counts
    fail_high_mt = pct_mt > max_pct_mt

    cell_passing_mask = (
        (~fail_low_genes)
        & (~fail_high_genes)
        & (~fail_low_counts)
        & (~fail_high_mt)
    )

    # Annotate cell-level flags
    adata.obs["is_retained_qc"] = cell_passing_mask
    adata.obs["qc_fail_low_genes"] = fail_low_genes
    adata.obs["qc_fail_high_genes"] = fail_high_genes
    adata.obs["qc_fail_low_counts"] = fail_low_counts
    adata.obs["qc_fail_high_mt"] = fail_high_mt

    # Detailed categoric QC verdict for easy inspection
    verdict = np.full(adata.n_obs, "Pass", dtype=object)
    verdict[fail_low_counts] = "Fail:LowCounts"
    verdict[fail_low_genes] = "Fail:LowGenes"
    verdict[fail_high_genes] = "Fail:HighGenes"
    verdict[fail_high_mt] = "Fail:HighMT"
    # Multiple failure reasons
    multi_fail = (
        fail_low_genes.astype(int)
        + fail_high_genes.astype(int)
        + fail_low_counts.astype(int)
        + fail_high_mt.astype(int)
    ) > 1
    verdict[multi_fail] = "Fail:MultiCriteria"
    adata.obs["qc_verdict"] = verdict

    # Gene-level QC mask: expressed in >= min_cells_per_gene
    X_mat = adata.X
    if sp.isspmatrix_csr(X_mat) or sp.isspmatrix_csc(X_mat):
        gene_expressed_cells = np.asarray((X_mat > 0).sum(axis=0)).ravel()
    else:
        gene_expressed_cells = np.sum(np.asarray(X_mat) > 0, axis=0)
    gene_passing_mask = gene_expressed_cells >= min_cells_per_gene
    adata.var["is_retained_qc_gene"] = gene_passing_mask
    adata.var["n_cells_by_counts"] = gene_expressed_cells

    # Step 5: Establish layers layout
    # Store raw integer counts in layers['counts'] as well as .X
    if not sp.isspmatrix_csr(adata.X):
        adata.X = sp.csr_matrix(adata.X, dtype=np.float32)
    else:
        adata.X = adata.X.astype(np.float32)
    adata.layers["counts"] = adata.X.copy()

    # Step 6: Annotate provenance metadata in uns
    n_total_cells = adata.n_obs
    n_passing_cells = int(np.sum(cell_passing_mask))
    n_failing_cells = n_total_cells - n_passing_cells
    pct_passing = (n_passing_cells / max(n_total_cells, 1)) * 100.0

    n_total_genes = adata.n_vars
    n_passing_genes = int(np.sum(gene_passing_mask))

    adata.uns["qc_spec"] = {
        "min_genes_per_cell": min_genes,
        "max_genes_per_cell": max_genes,
        "min_counts_per_cell": min_counts,
        "max_pct_mitochondrial": max_pct_mt,
        "min_cells_per_gene": min_cells_per_gene,
        "is_single_nucleus": is_sn,
    }
    adata.uns["sc_processing_summary"] = {
        "dataset_name": canonical_id,
        "n_cells_total": n_total_cells,
        "n_cells_passing_qc": n_passing_cells,
        "n_cells_failing_qc": n_failing_cells,
        "pct_cells_passing_qc": pct_passing,
        "n_genes_total": n_total_genes,
        "n_genes_passing_qc": n_passing_genes,
        "min_counts": min_counts,
        "min_genes": min_genes,
        "max_genes": max_genes,
        "max_pct_mt": max_pct_mt,
        "min_cells_per_gene": min_cells_per_gene,
    }

    # Clean index names for clean HDF5 serialization
    if adata.obs_names.name is not None and adata.obs_names.name in adata.obs.columns:
        adata.obs_names.name = None
    if adata.var_names.name is not None and adata.var_names.name in adata.var.columns:
        adata.var_names.name = None

    # Step 7: Serialize to output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_file = output_dir / f".tmp_{canonical_id}_unfiltered.h5ad"
    logger.info("[%s] Writing unfiltered .h5ad (%d cells, %d genes) to %s ...", canonical_id, n_total_cells, n_total_genes, out_file)
    adata.write_h5ad(temp_file, compression="gzip")
    temp_file.replace(out_file)

    elapsed = time.time() - t0
    final_size_mb = round(out_file.stat().st_size / (1024 * 1024), 2)
    logger.info(
        "[%s] COMPLETED in %.1f s. Output: %.1f MB (%d total cells, %d [%.1f%%] passing QC)",
        canonical_id,
        elapsed,
        final_size_mb,
        n_total_cells,
        n_passing_cells,
        pct_passing,
    )

    del adata
    gc.collect()

    return Success(
        ExportResult(
            cohort=canonical_id,
            total_cells=n_total_cells,
            passing_cells=n_passing_cells,
            failing_cells=n_failing_cells,
            pct_passing=round(pct_passing, 2),
            total_genes=n_total_genes,
            passing_genes=n_passing_genes,
            output_path=str(out_file),
            file_size_mb=final_size_mb,
            status="SUCCESS",
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export unfiltered AnnData (.h5ad) files with explicit QC boolean flags in .obs."
    )
    parser.add_argument(
        "--cohorts",
        nargs="+",
        default=list(DEFAULT_COHORTS),
        help="Cohort IDs to export (default: GSE326225 GSE302781 GSE222315 GSE301651).",
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=Path("/ritom/scratch/cscs/lhaeuser/data-test/raw"),
        help="Path to raw dataset directory on disk.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/tmp/bladder_unfiltered"),
        help="Output directory to save unfiltered .h5ad files (default: /tmp/bladder_unfiltered).",
    )
    parser.add_argument(
        "--skip-ensembl",
        action="store_true",
        help="Skip Ensembl Release 111 gene ID normalization.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )
    args = parser.parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    results: list[ExportResult] = []
    for c in args.cohorts:
        res = export_unfiltered_cohort(
            cohort_id=c,
            raw_dir=args.raw_dir,
            output_dir=output_dir,
            skip_ensembl=args.skip_ensembl,
            overwrite=args.overwrite,
        )
        match res:
            case Success(r):
                results.append(r)
            case Failure(err):
                logger.error("Failed cohort %s: %s", c, err)
                results.append(
                    ExportResult(
                        cohort=c,
                        total_cells=0,
                        passing_cells=0,
                        failing_cells=0,
                        pct_passing=0.0,
                        total_genes=0,
                        passing_genes=0,
                        output_path="",
                        file_size_mb=0.0,
                        status=f"FAILED: {err}",
                    )
                )

    print("\n" + "=" * 105)
    print("UNFILTERED BLADDER COHORTS EXPORT SUMMARY")
    print("=" * 105)
    print(
        f"{'Cohort':<14} {'Status':<16} {'Total Cells':<12} {'Passing':<10} {'Failing':<10} {'% Pass':<8} {'Genes':<8} {'Size (MB)':<10}"
    )
    print("-" * 105)
    for r in results:
        print(
            f"{r.cohort:<14} "
            f"{r.status[:15]:<16} "
            f"{r.total_cells:<12} "
            f"{r.passing_cells:<10} "
            f"{r.failing_cells:<10} "
            f"{r.pct_passing:<8.1f} "
            f"{r.total_genes:<8} "
            f"{r.file_size_mb:<10.1f}"
        )
    print("=" * 105)
    print(f"Files saved in: {output_dir}\n")


if __name__ == "__main__":
    main()
