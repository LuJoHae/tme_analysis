from __future__ import annotations

__version__ = "0.1.0"

import ici_datasets.bagaev_datasets
import ici_datasets.cbioportal_datasets
import ici_datasets.other_datasets

__all__ = [
    "__version__",
    "bagaev_datasets",
    "cbioportal_datasets",
    "other_datasets"
]