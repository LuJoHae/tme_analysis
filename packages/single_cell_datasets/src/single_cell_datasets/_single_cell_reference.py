import datalair
from pathlib import Path
import singlecellrnasignature
from gene_utils import norm_genes
import anndata as ad
from scipy.sparse import csr_matrix
import scanpy as sc
import pandas as pd
import numpy as np
import single_cell_datasets
import gene_utils
from single_cell_datasets import SingleCellDataProcessStep07


class ReferenceSingleCell(datalair.Dataset):
    """Datalair Dataset class for all single cell datasets."""

    def __init__(self) -> None:
        """Initialize this dataset class as a datalair.Dataset class with namespace `DatasetSingleCell`."""
        super().__init__(namespace="DatasetSingleCell")


class SingleCellReference(ReferenceSingleCell):
    storage_path = Path("/storage/halu").resolve()

    def derive(self, lair: datalair.Lair) -> None:
        ad.settings.allow_write_nullable_strings = True
        output_dir = lair.get_path(self)
        ds = SingleCellDataProcessStep07()
        lair.get_dataset_filepaths(ds)
        adata = ad.read_h5ad(lair.get_dataset_filepaths(ds)["adata.h5ad"])



        adata.write_h5ad(output_dir.joinpath("adata.h5ad"))