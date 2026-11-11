"""Batch execution orchestrator for Luecken & Theis (2019) single-cell QC & normalization inside tme_datasets."""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path
import sys
import time
from typing import Any, Sequence

import anndata as ad  # type: ignore[import-untyped]
import polars as pl
from returns.maybe import Maybe, Nothing, Some
from returns.result import Failure, Result, Success

from ..logging import get_logger
from ..models import QualityControlSpec, SingleCellProcessingSpec
from ..paths import find_dataset_h5ad, get_data_paths, get_preprocessed_h5ad_path
from ..preprocessing.pipeline import process_single_cell_dataset
from ..query import load_dataset
from ..registry import list_registered_datasets
from ..types import Modality

ad.settings.allow_write_nullable_strings = True
logger = get_logger("pipeline.qc_orchestrator")


def discover_available_single_cell_cohorts(
    explicit_cohorts: Maybe[tuple[str, ...]] = Nothing,
) -> tuple[str, ...]:
    """Identify which registered single-cell datasets have available data assets on disk."""
    registered = list_registered_datasets().filter(modality=Modality.SINGLE_CELL)
    candidate_ids = tuple(s.id for s in registered)

    match explicit_cohorts:
        case Some(requested) if len(requested) > 0:
            return tuple(cid for cid in requested if cid in candidate_ids)
        case _:
            pass

    cfg = get_data_paths()
    available: list[str] = []
    for cid in candidate_ids:
        if isinstance(find_dataset_h5ad(cid), Some):
            available.append(cid)
        else:
            raw_dir = cfg.raw_dir / cid
            if raw_dir.exists() and any(raw_dir.iterdir()):
                available.append(cid)
    return tuple(available)


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
    plot_dir = output_reports_dir / dataset_id
    plot_dir.mkdir(parents=True, exist_ok=True)

    load_res = load_dataset(
        dataset_id,
        force_recompute=force_recompute,
        apply_qc=False,
        normalize_ensembl=True,
    )
    match load_res:
        case Failure(err):
            return Failure(f"Failed to load dataset '{dataset_id}': {err}")
        case Success(raw_adata):
            pass

    if raw_adata.n_obs == 0 or raw_adata.n_vars == 0:
        return Failure(f"Dataset '{dataset_id}' has empty dimensions: {raw_adata.shape}.")

    qc_spec = QualityControlSpec(min_cells_per_gene=min_cells_per_gene)
    pipeline_spec = SingleCellProcessingSpec(
        qc=qc_spec,
        use_adaptive_qc=use_adaptive_qc,
        n_mads=n_mads,
        target_sum=target_sum,
        generate_plots=True,
        output_plot_dir=Some(plot_dir),
    )

    pipe_res = process_single_cell_dataset(
        raw_adata,
        spec=pipeline_spec,
        dataset_name=dataset_id,
        output_dir=plot_dir,
        display_plots=False,
    )
    match pipe_res:
        case Failure(err):
            return Failure(f"Pipeline execution failed for '{dataset_id}': {err}")
        case Success(result):
            processed_adata = result.adata

    target_h5ad = get_preprocessed_h5ad_path(dataset_id)
    target_h5ad.parent.mkdir(parents=True, exist_ok=True)
    processed_adata.write_h5ad(target_h5ad)
    size_mb = float(target_h5ad.stat().st_size) / float(1024 * 1024)
    duration = time.time() - t0

    record: dict[str, Any] = {
        "dataset_id": dataset_id,
        "n_cells_pre_qc": result.n_cells_pre_qc,
        "n_cells_post_qc": result.n_cells_post_qc,
        "pct_cells_retained": round(result.pct_cells_retained, 2),
        "n_genes_pre_qc": result.n_genes_pre_qc,
        "n_genes_post_qc": result.n_genes_post_qc,
        "min_counts": result.resolved_qc_spec.min_counts_per_cell,
        "max_counts": str(result.resolved_qc_spec.max_counts_per_cell.value_or(0) or "Unbounded"),
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
    records: Sequence[dict[str, Any]],
    output_reports_dir: Path,
) -> None:
    """Generate aggregate CSV and markdown reports summarizing quality control filtering."""
    if not records:
        return

    df = pl.DataFrame(list(records))
    csv_path = output_reports_dir / "single_cell_qc_run_summary.csv"
    df.write_csv(csv_path)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    md_lines = [
        "# Single-Cell RNA-Seq Quality Control & Normalization Run Summary",
        f"\n**Execution Timestamp**: {timestamp}  ",
        "**Pipeline Standard**: Luecken & Theis (2019) Best Practices  ",
        "**Normalization**: Linear TPM (scale=$10^6$), Natural $\\log(1+p)$, Raw integer counts in `.layers['counts']`  ",
        "**Gene Nomenclature**: Ensembl Release 111 Canonical IDs (`ENSG...`)  \n",
        "## Cohort Summary Table\n",
        "| Dataset ID | Pre-QC Cells | Post-QC Cells | Retained (%) | Pre-QC Genes | Post-QC Genes | Min Counts | Min Genes | Max Genes | Max MT (%) | Size (MB) | Time (s) |",
        "| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |",
    ]
    for r in records:
        md_lines.append(
            f"| `{r['dataset_id']}` | {r['n_cells_pre_qc']:,} | {r['n_cells_post_qc']:,} "
            f"| {r['pct_cells_retained']:.1f}% | {r['n_genes_pre_qc']:,} | {r['n_genes_post_qc']:,} "
            f"| {r['min_counts']:,} | {r['min_genes']:,} | {r['max_genes']:,} "
            f"| {r['max_pct_mt']:.1f}% | {r['file_size_mb']:.1f} | {r['duration_sec']:.1f} |"
        )

    md_path = output_reports_dir / "single_cell_qc_run_summary.md"
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")


def run_single_cell_qc_pipeline(
    *,
    output_dir: Path = Path("output/tme_datasets/qc"),
    explicit_cohorts: Maybe[tuple[str, ...]] = Nothing,
    n_mads: float = 3.0,
    target_sum: float = 1e6,
    min_cells_per_gene: int = 3,
    force_recompute: bool = False,
) -> Result[pl.DataFrame, str]:
    """Run Luecken & Theis single-cell QC across discovered or explicit cohorts."""
    output_dir.mkdir(parents=True, exist_ok=True)
    cohorts = discover_available_single_cell_cohorts(explicit_cohorts)
    if not cohorts:
        return Failure("No matching single-cell cohorts found with available data on disk.")

    records: list[dict[str, Any]] = []
    failures: list[str] = []
    for cid in cohorts:
        match execute_cohort_pipeline(
            dataset_id=cid,
            output_reports_dir=output_dir,
            use_adaptive_qc=True,
            n_mads=n_mads,
            target_sum=target_sum,
            min_cells_per_gene=min_cells_per_gene,
            force_recompute=force_recompute,
        ):
            case Success(rec):
                records.append(rec)
            case Failure(err):
                failures.append(f"{cid}: {err}")

    generate_summary_reports(records, output_dir)
    return (
        Failure("; ".join(failures))
        if failures and not records
        else Success(pl.DataFrame(records))
    )


def main() -> int:
    """CLI entrypoint for tme-datasets-qc."""
    parser = argparse.ArgumentParser(
        description="Execute Luecken & Theis (2019) QC & normalization pipeline across single-cell cohorts."
    )
    parser.add_argument(
        "--datasets",
        "-d",
        nargs="*",
        default=[],
        help="Optional dataset IDs to process (default: auto-discover all available single-cell cohorts).",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        type=Path,
        default=Path("output/tme_datasets/qc"),
        help="Output directory for reports and QC plots.",
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
    explicit: Maybe[tuple[str, ...]] = (
        Some(tuple(str(d) for d in args.datasets)) if args.datasets else Nothing
    )
    outcome = run_single_cell_qc_pipeline(
        output_dir=Path(args.output_dir),
        explicit_cohorts=explicit,
        n_mads=float(args.mads),
        target_sum=float(args.target_sum),
        min_cells_per_gene=int(args.min_cells_per_gene),
        force_recompute=bool(args.force_recompute),
    )
    match outcome:
        case Success(_):
            return 0
        case Failure(err):
            sys.stderr.write(f"ERROR: {err}\n")
            return 1
        case _:
            return 1


if __name__ == "__main__":
    sys.exit(main())
