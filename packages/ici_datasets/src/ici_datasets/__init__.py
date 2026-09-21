"""Immune checkpoint inhibitor (ICI) response datasets.

Deprecated: ici_datasets is superseded by tme_datasets.
Please use `from tme_datasets import load_dataset, query_datasets` instead.
"""

from __future__ import annotations

__version__ = "0.1.0"

from tme_datasets import load_dataset, query_datasets
from tme_datasets.compat import CBioPortalDataset

__all__ = [
    "__version__",
    "bagaev_datasets",
    "cbioportal_datasets",
    "other_datasets"
]