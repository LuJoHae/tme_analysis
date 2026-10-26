#!/usr/bin/env python3
"""Batch execution orchestrator for Luecken & Theis (2019) single-cell pre-processing pipeline.

Processes all available single-cell transcriptomic datasets across the configured data hierarchy:
  1. Computes cell-level QC covariates (counts, genes, mito %, ribo %, hemoglobin %).
  2. Derives adaptive MAD thresholds (+/- 3 MADs) for outlier detection.
  3. Renders 5-panel Luecken & Theis QC diagnostic dashboards exported to vector SVG and 300 DPI PNG.
  4. Filters non-viable droplets, dying cells, and unexpressed genes (< 3 cells).
  5. Scales library sizes to TPM (10^6) in .layers['tpm'].
  6. Applies natural log1p transformation in .layers['log1p'] and .X.
  7. Archives raw integer counts in .layers['counts'] and .raw.
  8. Serializes standardized AnnData to data/preprocessed/{dataset_id}.h5ad with Ensembl release 111 IDs.
  9. Compiles an aggregate markdown and CSV summary audit report across all cohorts.
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import sys
import time
from typing import Any, Sequence

import anndata as ad
import polars as pl
from returns.maybe import Nothing, Some
from returns.result import Failure, Result, Success

# Enable writing nullable string columns in AnnData
ad.settings.allow_write_nullable_strings = True

from tme_datasets import (
    Modality,
    QualityControlSpec,
    SingleCellProcessingSpec,
    find_dataset_h5ad,
    get_data_paths,
    get_preprocessed_h5ad_path,
    list_registered_datasets,
    load_dataset,
    process_single_cell_dataset,
)
from tme_datasets.logging import get_logger

logger = get_logger("orchestrator.single_cell_qc")


def discover_available_single_cell_cohorts(
    explicit_cohorts: Sequence[str] | None = None,
) -> list[str]:
    """Identify which registered single-cell datasets have available data assets on disk."""
    registered = list_registered_datasets().filter(modality=Modality.SINGLE_CELL)
    candidate_ids = [s.id for s in registered]

    if explicit_cohorts:
        # Validate requested cohorts
        return [cid for cid in explicit_cohorts if cid in candidate_ids]

    available: list[str] = []
    for cid in candidate_ids:
        found = find_dataset_h5ad(cid)
        if isinstance(found, Some):
            available.append(cid)
        else:
            # Check raw directory
            cfg = get_data_paths()
            raw_dir = cfg.raw_dir / cid
            if raw_dir.exists() and any(raw_dir.iterdir()):
                available.append(cid)

    return available


def execute_cohort_pipeline(
    dataset_id: str,
    output_reports_dir: Path,
    use_adaptive_qc: bool = True,
    n_mads: float = 3.0,
    target_sum: float = 1e6,
    min_cells_per_gene: int = 3,
    force_recompute: bool = False,
) -> Result[dict[str, Any], str]:
    """Execute complete Luecken & Theis pipeline on a single cohort."""
    t0 = time.time()
    logger.info("=" * 70)
    logger.info("Processing single-cell cohort: %s", dataset_id)
    logger.info("=" * 70)

    plot_dir = output_reports_dir / dataset_id
    plot_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load dataset (normalize to Ensembl Release 111)
    logger.info("[%s] Loading dataset into memory...", dataset_id)
    load_res = load_dataset(
        dataset_id,
        force_recompute=force_recompute,
        apply_qc=False,  # QC is executed within pipeline
        normalize_ensembl=True,
    )
    match load_res:
        case Failure(err):
            msg = f"Failed to load dataset '{dataset_id}': {err}"
            logger.error(msg)
            return Failure(msg)
        case Success(adata):
            raw_adata = adata

    if raw_adata.n_obs == 0 or raw_adata.n_vars == 0:
        msg = f"Dataset '{dataset_id}' has empty dimensions: {raw_adata.shape}. Skipping."
        logger.warning(msg)
        return Failure(msg)

    logger.info("[%s] Ingested matrix: %d cells x %d genes", dataset_id, raw_adata.n_obs, raw_adata.n_vars)

    # 2. Build Pipeline Specification
    qc_spec = QualityControlSpec(
        min_cells_per_gene=min_cells_per_gene,
    )
    pipeline_spec = SingleCellProcessingSpec(
        qc=qc_spec,
        use_adaptive_qc=use_adaptive_qc,
        n_mads=n_mads,
        target_sum=target_sum,
        generate_plots=True,
        output_plot_dir=Some(plot_dir),
    )

    # 3. Execute Processing Pipeline
    logger.info("[%s] Running QC covariate calculation, outlier detection, and normalization...", dataset_id)
    pipe_res = process_single_cell_dataset(
        raw_adata,
        spec=pipeline_spec,
        dataset_name=dataset_id,
        output_dir=plot_dir,
        display_plots=False,
    )
    match pipe_res:
        case Failure(err):
            msg = f"Pipeline execution failed for '{dataset_id}': {err}"
            logger.error(msg)
            return Failure(msg)
        case Success(result):
            processed_adata = result.adata

    # 4. Serialize Processed H5AD
    target_h5ad = get_preprocessed_h5ad_path(dataset_id)
    target_h5ad.parent.mkdir(parents=True, exist_ok=True)
    logger.info("[%s] Writing processed AnnData to %s...", dataset_id, target_h5ad)
    processed_adata.write_h5ad(target_h5ad)
    size_mb = target_h5ad.stat().st_size / (1024 * 1024)

    duration = time.time() - t0
    logger.info(
        "[%s] Completed in %.1fs: %d -> %d cells (%.1f%% retained), %d -> %d genes, size=%.1f MB",
        dataset_id,
        duration,
        result.n_cells_pre_qc,
        result.n_cells_post_qc,
        result.pct_cells_retained,
        result.n_genes_pre_qc,
        result.n_genes_post_qc,
        size_mb,
    )

    record = {
        "dataset_id": dataset_id,
        "n_cells_pre_qc": result.n_cells_pre_qc,
        "n_cells_post_qc": result.n_cells_post_qc,
        "pct_cells_retained": round(result.pct_cells_retained, 2),
        "n_genes_pre_qc": result.n_genes_pre_qc,
        "n_genes_post_qc": result.n_genes_post_qc,
        "min_counts": result.resolved_qc_spec.min_counts_per_cell,
        "max_counts": result.resolved_qc_spec.max_counts_per_cell.value_or("None"),
        "min_genes": result.resolved_qc_spec.min_genes_per_cell,
        "max_genes": result.resolved_qc_spec.max_genes_per_cell,
        "max_pct_mt": round(result.resolved_qc_spec.max_pct_mitochondrial, 2),
        "min_cells_per_gene": result.resolved_qc_spec.min_cells_per_gene,
        "h5ad_path": str(target_h5ad),
        "file_size_mb": round(size_mb, 1),
        "duration_sec": round(duration, 1),
        "n_plots": len(result.plot_paths),
    }
    return Success(record)


def generate_summary_reports(
    records: list[dict[str, Any]],
    output_reports_dir: Path,
) -> None:
    """Generate aggregate CSV and markdown reports summarizing quality control filtering."""
    if not records:
        return

    df = pl.DataFrame(records)
    csv_path = output_reports_dir / "single_cell_qc_run_summary.csv"
    df.write_csv(csv_path)
    logger.info("Saved QC summary CSV: %s", csv_path)

    # Format Markdown Report
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    md_lines = [
        "# Single-Cell RNA-Seq Quality Control & Normalization Run Summary",
        f"\n**Execution Timestamp**: {timestamp}  ",
        f"**Pipeline Standard**: Luecken & Theis (2019) Best Practices  ",
        f"**Normalization**: Linear TPM (scale=$10^6$), Natural $\\log(1+p)$, Raw integer counts in `.layers['counts']`  ",
        f"**Gene Nomenclature**: Ensembl Release 111 Canonical IDs (`ENSG...`)  \n",
        "## Cohort Summary Table\n",
        "| Dataset ID | Pre-QC Cells | Post-QC Cells | Retained (%) | Pre-QC Genes | Post-QC Genes | Min Counts | Min Genes | Max Genes | Max MT (%) | Size (MB) | Time (s) |",
        "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
    ]

    for r in records:
        line = (
            f"| `{r['dataset_id']}` "
            f"| {r['n_cells_pre_qc']:,} "
            f"| {r['n_cells_post_qc']:,} "
            f"| {r['pct_cells_retained']:.1f}% "
            f"| {r['n_genes_pre_qc']:,} "
            f"| {r['n_genes_post_qc']:,} "
            f"| {r['min_counts']:,} "
            f"| {r['min_genes']:,} "
            f"| {r['max_genes']:,} "
            f"| {r['max_pct_mt']:.1f}% "
            f"| {r['file_size_mb']:.1f} "
            f"| {r['duration_sec']:.1f} |"
        )
        md_lines.append(line)

    md_lines.append("\n## Diagnostic Plot Locations\n")
    for r in records:
        md_lines.append(f"- **`{r['dataset_id']}`**: `output/reports/qc/{r['dataset_id']}/{r['dataset_id']}_qc_dashboard.svg` and `.png`")

    md_content = "\n".join(md_lines) + "\n"
    md_path = output_reports_dir / "single_cell_qc_run_summary.md"
    md_path.write_text(md_content)
    logger.info("Saved QC summary Markdown report: %s", md_path)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Execute Luecken & Theis (2019) QC & normalization pipeline across all available single-cell cohorts."
    )
    parser.add_argument(
        "--datasets",
        "-d",
        nargs="*",
        default=None,
        help="Optional list of dataset IDs to process. If omitted, auto-discovers all available single-cell cohorts.",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        type=Path,
        default=Path("output/reports/qc"),
        help="Output directory for reports and QC plots (default: output/reports/qc).",
    )
    parser.add_argument(
        "--mads",
        type=float,
        default=3.0,
        help="Number of Median Absolute Deviations for adaptive outlier detection (default: 3.0).",
    )
    parser.add_argument(
        "--target-sum",
        type=float,
        default=1e6,
        help="Target library size for TPM normalization (default: 1e6).",
    )
    parser.add_argument(
        "--min-cells-per-gene",
        type=int,
        default=3,
        help="Minimum number of expressing cells per gene (default: 3).",
    )
    parser.add_argument(
        "--force-recompute",
        action="store_true",
        help="Force recomputation from raw assets instead of cached H5AD.",
    )

    args = parser.parse_args()

    reports_dir = args.output_dir
    reports_dir.mkdir(parents=True, exist_ok=True)

    cohorts = discover_available_single_cell_cohorts(args.datasets)
    if not cohorts:
        logger.error("No matching single-cell cohorts found with available data on disk.")
        return 1

    logger.info("Found %d single-cell cohort(s) to process: %s", len(cohorts), cohorts)

    records: list[dict[str, Any]] = []
    failures: list[tuple[str, str]] = []

    for idx, cid in enumerate(cohorts, start=1):
        logger.info("[%d/%d] Beginning processing for '%s'...", idx, len(cohorts), cid)
        res = execute_cohort_pipeline(
            dataset_id=cid,
            output_reports_dir=reports_dir,
            use_adaptive_qc=True,
            n_mads=args.mads,
            target_sum=args.target_sum,
            min_cells_per_gene=args.min_cells_per_gene,
            force_recompute=args.force_recompute,
        )
        match res:
            case Success(record):
                records.append(record)
            case Failure(err):
                logger.error("Cohort '%s' failed: %s", cid, err)
                failures.append((cid, err))

    generate_summary_reports(records, reports_dir)

    print("\n" + "=" * 70)
    print(f"BATCH QC RUN COMPLETE: {len(records)} succeeded, {len(failures)} failed.")
    print(f"Summary Markdown: {reports_dir / 'single_cell_qc_run_summary.md'}")
    print(f"Summary CSV:      {reports_dir / 'single_cell_qc_run_summary.csv'}")
    print("=" * 70)

    if failures:
        print("\nFailures:")
        for cid, err in failures:
            print(f"  - {cid}: {err}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
