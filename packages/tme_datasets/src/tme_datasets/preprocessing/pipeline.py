"""End-to-end single-cell RNA-seq pre-processing pipeline conforming to Luecken & Theis (2019)."""

from __future__ import annotations

from pathlib import Path
import gc
import anndata as ad
from returns.maybe import Some
from returns.result import Failure, Result, Success

from ..logging import get_logger
from ..models import QualityControlSpec, SingleCellProcessingResult, SingleCellProcessingSpec
from .normalization import standardize_processed_layers
from .qc import (
    apply_quality_control,
    calculate_adaptive_thresholds,
    compute_qc_covariates,
    extract_qc_metrics_dataframe,
)
from .qc_plots import build_qc_dashboard, display_qc_plots_inline, export_qc_plots

logger = get_logger("preprocessing.pipeline")


def process_single_cell_dataset(
    adata: ad.AnnData,
    spec: SingleCellProcessingSpec | None = None,
    dataset_name: str = "single_cell",
    output_dir: Path | None = None,
    display_plots: bool = True,
) -> Result[SingleCellProcessingResult, str]:
    """Execute complete best-practices single-cell pre-processing pipeline per Luecken & Theis (2019).

    Pipeline stages:
        1. Calculate cellular QC covariates: count depth, gene detection, mitochondrial, ribosomal, and hemoglobin %.
        2. Determine QC thresholds via hybrid adaptive statistical detection (±3 MADs) or explicit thresholds.
        3. Construct and export the 5-panel Luecken & Theis Fig 2 QC dashboard to vector SVG / PNG and inline display.
        4. Filter non-viable droplets, dying cells, multiplets, and unexpressed genes (expressed in <3 cells).
        5. Scale cell library sizes to linear TPM/CPM (default 10^6) in .layers['tpm'].
        6. Apply natural log1p transformation in .layers['log1p'] and .X.
        7. Archive raw matrix in .layers['counts'] and .raw.

    Args:
        adata: Raw input AnnData object with counts in .X.
        spec: SingleCellProcessingSpec configuration.
        dataset_name: Dataset identifier for titles and plot filenames.
        output_dir: Directory where SVG and PNG plots will be exported.
        display_plots: Whether to render plots inline in Jupyter/IPython sessions.

    Returns:
        Success(SingleCellProcessingResult) or Failure(error_message).
    """
    cfg = spec or SingleCellProcessingSpec()

    # Step 1: Compute QC Covariates on raw data
    logger.info("[%s] Step 1/4: Computing QC covariates on %d cells...", dataset_name, adata.n_obs)
    with_covariates = compute_qc_covariates(adata)

    # Step 2: Determine QC Thresholds
    if cfg.use_adaptive_qc:
        logger.info("[%s] Step 2/4: Determining adaptive QC thresholds (%.1f MADs)...", dataset_name, cfg.n_mads)
        resolved_spec = calculate_adaptive_thresholds(with_covariates, n_mads=cfg.n_mads, base_spec=cfg.qc)
    else:
        resolved_spec = cfg.qc

    # Step 3: Diagnostic Visualizations
    plot_paths: tuple[Path, ...] = ()
    if cfg.generate_plots:
        logger.info("[%s] Step 3/4: Generating Luecken & Theis QC diagnostic plots...", dataset_name)
        metrics_df = extract_qc_metrics_dataframe(with_covariates)
        dashboard_chart = build_qc_dashboard(metrics_df, resolved_spec, dataset_name=dataset_name)

        if display_plots:
            display_qc_plots_inline(dashboard_chart)

        target_dir = output_dir or (cfg.output_plot_dir.unwrap() if isinstance(cfg.output_plot_dir, Some) else None)
        if target_dir is not None:
            export_res = export_qc_plots(dashboard_chart, target_dir, dataset_name=dataset_name)
            match export_res:
                case Success(paths):
                    plot_paths = paths
                case Failure(err):
                    logger.warning("Could not export QC plots: %s", err)

    # Step 4: Quality Control Filtering
    logger.info("[%s] Step 4/4: Filtering cells and genes...", dataset_name)
    qc_res = apply_quality_control(
        with_covariates,
        qc_spec=resolved_spec,
        use_adaptive_qc=False,  # Already resolved
        min_cells_per_gene=resolved_spec.min_cells_per_gene,
    )
    del with_covariates
    gc.collect()

    match qc_res:
        case Failure(err):
            return Failure(err)
        case Success(filtered_adata):
            pass

    # Step 5: TPM Normalization and Log1p Transformation
    logger.info("[%s] Standardizing multi-layer AnnData (TPM scale=%.0e, log1p)...", dataset_name, cfg.target_sum)
    norm_res = standardize_processed_layers(filtered_adata, target_sum=cfg.target_sum)
    del filtered_adata
    gc.collect()

    match norm_res:
        case Failure(err):
            return Failure(err)
        case Success(processed_adata):
            pass

    # Annotate provenance in uns
    n_pre_cells = adata.n_obs
    n_post_cells = processed_adata.n_obs
    n_pre_genes = adata.n_vars
    n_post_genes = processed_adata.n_vars
    pct_retained = (n_post_cells / max(n_pre_cells, 1)) * 100.0

    summary = {
        "dataset_name": dataset_name,
        "n_cells_pre_qc": n_pre_cells,
        "n_cells_post_qc": n_post_cells,
        "pct_cells_retained": pct_retained,
        "n_genes_pre_qc": n_pre_genes,
        "n_genes_post_qc": n_post_genes,
        "target_sum": cfg.target_sum,
        "min_counts": resolved_spec.min_counts_per_cell,
        "max_counts": resolved_spec.max_counts_per_cell.value_or(None),
        "min_genes": resolved_spec.min_genes_per_cell,
        "max_genes": resolved_spec.max_genes_per_cell,
        "max_pct_mt": resolved_spec.max_pct_mitochondrial,
        "min_cells_per_gene": resolved_spec.min_cells_per_gene,
        "plot_paths": [str(p) for p in plot_paths],
    }
    processed_adata.uns["sc_processing_summary"] = summary

    return Success(
        SingleCellProcessingResult(
            adata=processed_adata,
            resolved_qc_spec=resolved_spec,
            plot_paths=plot_paths,
            n_cells_pre_qc=n_pre_cells,
            n_cells_post_qc=n_post_cells,
            n_genes_pre_qc=n_pre_genes,
            n_genes_post_qc=n_post_genes,
            pct_cells_retained=pct_retained,
            summary=summary,
        )
    )


__all__ = ["process_single_cell_dataset"]
