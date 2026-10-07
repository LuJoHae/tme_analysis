#!/usr/bin/env python3
"""CLI utility to execute best-practices single-cell pre-processing pipeline per Luecken & Theis (2019)."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from returns.maybe import Nothing, Some
from returns.result import Failure, Success

from tme_datasets import (
    QualityControlSpec,
    SingleCellProcessingSpec,
    get_preprocessed_h5ad_path,
    load_dataset,
    process_single_cell_dataset,
)
from tme_datasets.logging import get_logger

logger = get_logger("cli.preprocess_single_cell")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run best-practices scRNA-seq quality control, TPM normalization, and log1p pipeline (Luecken & Theis 2019)."
    )
    parser.add_argument(
        "--dataset",
        "-d",
        type=str,
        required=True,
        help="Registered dataset ID (e.g. GSE120575, GSE125449, Maynard_NSCLC).",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        type=Path,
        default=None,
        help="Output directory for SVG and PNG plots (default: output/reports/qc/<dataset>).",
    )
    parser.add_argument(
        "--output-h5ad",
        type=Path,
        default=None,
        help="Optional explicit path to write the processed H5AD file.",
    )
    parser.add_argument(
        "--no-adaptive",
        action="store_true",
        help="Disable adaptive MAD outlier filtering and rely strictly on fixed thresholds.",
    )
    parser.add_argument(
        "--mads",
        type=float,
        default=3.0,
        help="Number of Median Absolute Deviations for adaptive filtering (default: 3.0).",
    )
    parser.add_argument("--min-counts", type=int, default=500, help="Minimum counts per cell (default: 500).")
    parser.add_argument("--max-counts", type=int, default=None, help="Maximum counts per cell (doublet filter).")
    parser.add_argument("--min-genes", type=int, default=200, help="Minimum detected genes per cell (default: 200).")
    parser.add_argument("--max-genes", type=int, default=8000, help="Maximum detected genes per cell (default: 8000).")
    parser.add_argument("--max-mt", type=float, default=20.0, help="Maximum mitochondrial count percentage (default: 20.0).")
    parser.add_argument("--min-cells-per-gene", type=int, default=3, help="Minimum cells per gene (default: 3).")
    parser.add_argument("--target-sum", type=float, default=1e6, help="Target library size for TPM/CPM scaling (default: 1e6).")

    args = parser.parse_args()

    plot_dir = args.output_dir or Path(f"output/reports/qc/{args.dataset}")
    plot_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"Single-Cell Processing Pipeline: {args.dataset}")
    print("=" * 70)

    # 1. Load dataset (raw counts or preprocessed baseline)
    print(f"[1/3] Loading dataset '{args.dataset}'...")
    res = load_dataset(
        args.dataset,
        force_recompute=False,
        apply_qc=False,  # QC will be handled cleanly by the pipeline
        normalize_ensembl=True,
    )
    match res:
        case Failure(err):
            print(f"Error loading dataset: {err}", file=sys.stderr)
            return 1
        case Success(adata):
            raw_adata = adata

    print(f"   Loaded: {raw_adata.n_obs} cells x {raw_adata.n_vars} genes")

    # 2. Build specification
    qc_spec = QualityControlSpec(
        min_counts_per_cell=args.min_counts,
        max_counts_per_cell=Some(args.max_counts) if args.max_counts is not None else Nothing,
        min_genes_per_cell=args.min_genes,
        max_genes_per_cell=args.max_genes,
        max_pct_mitochondrial=args.max_mt,
        min_cells_per_gene=args.min_cells_per_gene,
    )
    pipeline_spec = SingleCellProcessingSpec(
        qc=qc_spec,
        use_adaptive_qc=not args.no_adaptive,
        n_mads=args.mads,
        target_sum=args.target_sum,
        generate_plots=True,
        output_plot_dir=Some(plot_dir),
    )

    # 3. Execute pipeline
    print(f"[2/3] Running QC diagnostics, filtering, TPM normalization, and log1p...")
    pipe_res = process_single_cell_dataset(
        raw_adata,
        spec=pipeline_spec,
        dataset_name=args.dataset,
        output_dir=plot_dir,
        display_plots=False,
    )

    match pipe_res:
        case Failure(err):
            print(f"Pipeline execution failed: {err}", file=sys.stderr)
            return 1
        case Success(result):
            processed_adata = result.adata

    print(f"[3/3] Execution completed successfully!")
    print("-" * 70)
    print(f"Summary for {args.dataset}:")
    print(f"  Cells pre-QC:         {result.n_cells_pre_qc:,}")
    print(f"  Cells post-QC:        {result.n_cells_post_qc:,} ({result.pct_cells_retained:.1f}% retained)")
    print(f"  Genes pre-QC:         {result.n_genes_pre_qc:,}")
    print(f"  Genes post-QC:        {result.n_genes_post_qc:,}")
    print(f"  Resolved min counts:  {result.resolved_qc_spec.min_counts_per_cell:,}")
    print(f"  Resolved min genes:   {result.resolved_qc_spec.min_genes_per_cell:,}")
    print(f"  Resolved max genes:   {result.resolved_qc_spec.max_genes_per_cell:,}")
    print(f"  Resolved max mito %:  {result.resolved_qc_spec.max_pct_mitochondrial:.1f}%")
    print(f"  Generated plots:")
    for p in result.plot_paths:
        print(f"    - {p}")

    target_h5ad = args.output_h5ad or get_preprocessed_h5ad_path(f"{args.dataset}_sc_processed")
    target_h5ad.parent.mkdir(parents=True, exist_ok=True)
    processed_adata.write_h5ad(target_h5ad)
    print(f"  Saved processed H5AD: {target_h5ad} ({target_h5ad.stat().st_size / (1024*1024):.1f} MB)")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
