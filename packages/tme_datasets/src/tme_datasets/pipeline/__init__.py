"""Dataset ingestion, batch preprocessing, and Luecken–Theis single-cell QC pipelines."""

from __future__ import annotations

from .batch_preprocess import (
    CohortPreprocessTask,
    MemoryGovernor,
    execute_preprocess_task,
    run_batch_preprocess_cohorts,
)
from .download_preprocess import (
    DatasetProcessResult,
    PipelineConfig,
    download_and_preprocess_all,
    run_parallel_download_and_preprocess,
)
from .qc_orchestrator import (
    discover_available_single_cell_cohorts,
    execute_cohort_pipeline,
    run_single_cell_qc_pipeline,
)

__all__ = [
    "CohortPreprocessTask",
    "DatasetProcessResult",
    "MemoryGovernor",
    "PipelineConfig",
    "discover_available_single_cell_cohorts",
    "download_and_preprocess_all",
    "execute_cohort_pipeline",
    "execute_preprocess_task",
    "run_batch_preprocess_cohorts",
    "run_parallel_download_and_preprocess",
    "run_single_cell_qc_pipeline",
]
