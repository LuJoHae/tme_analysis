import datalair
from pathlib import Path
import pandas as pd


class DatasetGenentechAlign(datalair.Dataset):
    """Datalair Dataset class for all Datasets from Genentech."""

    def __init__(self) -> None:
        """Initialize this dataset class as a datalair.Dataset class with namespace `DatasetGenentecAlign`."""
        super().__init__(namespace="DatasetGenentecAlign")


class EGAD00001006631(DatasetGenentechAlign):
    storage_path = Path("/storage/halu").resolve()

    def derive(self, lair: datalair.Lair) -> None:
        output_dir = lair.get_path(self)
        all_counts = []
        assert self.storage_path.exists()
        for dirpath in (self.storage_path / "manual-download/EGAD00001006631-align").iterdir():
            patient_id = dirpath.name
            counts = pd.read_csv(dirpath / "ReadsPerGene.out.tab", sep="\t", header=None, skiprows=4, index_col=0)
            counts = counts.rename(columns={0: None, 1: "unstranded", 2: "stranded_forward", 3: "stranded_reverse"})
            counts.index.name = "gene_id"
            assert counts["stranded_reverse"].mean() >= counts["stranded_forward"].mean()
            # assert counts["stranded_reverse"].mean() >= counts["unstranded"].mean()
            counts = counts["stranded_reverse"]
            counts.name = patient_id
            all_counts.append(counts)
        all_counts = pd.concat(all_counts, axis=1).T
        all_counts = all_counts[sorted(list(all_counts.columns))]
        all_counts.to_hdf(output_dir / "counts.h5pd", key="counts")