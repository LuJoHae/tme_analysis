"""Package single_cell_immuno_datasets: ICB single-cell data downloader, preprocessor, and report generator.

Deprecated: single_cell_immuno_datasets is superseded by tme_datasets.
Please use `from tme_datasets import load_dataset, query_datasets` instead.
"""

from .config import (
    DatasetSpec,
    QualityControlSpec,
    DataDirectories,
    TIER_1_DATASETS,
)
from .downloader import (
    download_single_file,
    download_dataset,
    download_all_tier1,
)
from .preprocessor import (
    check_transform_state,
    harmonize_metadata,
    preprocess_anndata,
)
from .gondal2025 import (
    download_gondal2025,
    preprocess_gondal2025,
)
from .gse120575 import (
    fetch_and_format_gse120575,
    download_gse120575,
    process_to_parquet,
)
from .report import (
    generate_summary_report,
    ALL_EVALUATED_DATASETS,
)
from .milopy_checker import (
    MilopyCheckResult,
    verify_milopy_compatibility,
)
from .milo_analysis import (
    MiloDatasetConfig,
    MiloRunSummary,
    READY_DATASET_CONFIGS,
    run_single_milo_pipeline,
    plot_cohort_summary,
    plot_cohort_percentage_summary,
    load_existing_summary,
    compute_and_save_umaps,
    ensure_umaps_for_dataset,
)

__all__ = [
    "DatasetSpec",
    "QualityControlSpec",
    "DataDirectories",
    "TIER_1_DATASETS",
    "download_single_file",
    "download_dataset",
    "download_all_tier1",
    "check_transform_state",
    "harmonize_metadata",
    "preprocess_anndata",
    "download_gondal2025",
    "preprocess_gondal2025",
    "fetch_and_format_gse120575",
    "download_gse120575",
    "process_to_parquet",
    "generate_summary_report",
    "ALL_EVALUATED_DATASETS",
    "MilopyCheckResult",
    "verify_milopy_compatibility",
    "MiloDatasetConfig",
    "MiloRunSummary",
    "READY_DATASET_CONFIGS",
    "run_single_milo_pipeline",
    "plot_cohort_summary",
    "plot_cohort_percentage_summary",
    "load_existing_summary",
    "compute_and_save_umaps",
    "ensure_umaps_for_dataset",
]
