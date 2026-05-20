import datalair
from pathlib import Path
import pandas as pd


class DatasetReferences(datalair.Dataset):
    """Datalair Dataset class for all Datasets containing gene_utils."""

    def __init__(self) -> None:
        """Initialize this dataset class as a datalair.Dataset class with namespace `DatasetReferences`."""
        super().__init__(namespace="DatasetGenentecAlign")


class Reference(DatasetReferences):
    storage_path = Path("/storage/halu").resolve()

    def derive(self, lair: datalair.Lair) -> None:
        output_dir = lair.get_path(self)

