#!/usr/bin/env python3
"""Compatibility entrypoint delegating to tme_datasets.pipeline.download_preprocess."""

from __future__ import annotations

import sys

from tme_datasets.pipeline.download_preprocess import (
    DatasetProcessResult,
    PipelineConfig,
    download_and_preprocess_all,
    main,
    run_parallel_download_and_preprocess,
)

__all__ = [
    "DatasetProcessResult",
    "PipelineConfig",
    "download_and_preprocess_all",
    "main",
    "run_parallel_download_and_preprocess",
]

if __name__ == "__main__":
    sys.exit(main())
