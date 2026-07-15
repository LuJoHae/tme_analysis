"""Processed datasets of scRNA-seq data from patient tumors"""

from __future__ import annotations

__version__ = "0.1.0"

from single_cell_datasets._single_cell_datasets import (
    SingleCellDataProcessStep01,
    SingleCellDataProcessStep02,
    SingleCellDataProcessStep03,
    SingleCellDataProcessStep04,
    SingleCellDataProcessStep05,
    SingleCellDataProcessStep06,
    SingleCellDataProcessStep07
)

from single_cell_datasets._single_cell_reference import (
    SingleCellReference,
    SingleCellClustered,
    SingleCellSubclustered,
    SingleCellClusterMeans,
    SingleCellDeconvolution,
    SingleCellBagaevDeconvolution
)

__all__ = [
    "__version__",
    "SingleCellDataProcessStep01",
    "SingleCellDataProcessStep02",
    "SingleCellDataProcessStep03",
    "SingleCellDataProcessStep04",
    "SingleCellDataProcessStep05",
    "SingleCellDataProcessStep06",
    "SingleCellDataProcessStep07",
    "SingleCellReference",
    "SingleCellClustered",
    "SingleCellSubclustered",
    "SingleCellClusterMeans",
    "SingleCellDeconvolution",
    "SingleCellBagaevDeconvolution"
]