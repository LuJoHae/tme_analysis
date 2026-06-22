import datalair
from gene_utils import download
import pandas as pd
import numpy as np


class DatasetBagaev(datalair.Dataset):
    """Datalair Dataset class for all single cell datasets."""

    def __init__(self) -> None:
        """Initialize this dataset class as a datalair.Dataset class with namespace `DatasetSingleCell`."""
        super().__init__(namespace="DatasetBagaev")


class Signature(DatasetBagaev):
    def derive(self, lair: datalair.Lair) -> None:
        output_dir = lair.get_path(self)
        download(
            "https://raw.githubusercontent.com/BostonGene/MFP/refs/heads/master/signatures/gene_signatures.gmt",
            output_dir.joinpath("gene_signatures.gmt"),
        )
        download(
            "https://raw.githubusercontent.com/BostonGene/MFP/refs/heads/master/signatures/gene_signatures_order.tsv",
            output_dir.joinpath("gene_signatures_order.tsv"),
        )


class PanMelanoma(DatasetBagaev):
    def derive(self, lair: datalair.Lair) -> None:
        output_dir = lair.get_path(self)
        download(
            "https://science.bostongene.com/tumor-portrait/api/annotation_panmi/",
            output_dir.joinpath("annotations.tsv"),
        )
        download(
            "https://science.bostongene.com/tumor-portrait/api/signatures_panmi/",
            output_dir.joinpath("signature.tsv"),
        )


class TCGA(DatasetBagaev):
    def derive(self, lair: datalair.Lair) -> None:
        output_dir = lair.get_path(self)
        download(
            "https://science.bostongene.com/tumor-portrait/api/annotation_tcga/",
            output_dir.joinpath("annotations.tsv"),
        )
        download(
            "https://science.bostongene.com/tumor-portrait/api/signatures_tcga/",
            output_dir.joinpath("signature.tsv"),
        )


def get_pan_melanoma_bagaev_dataset(lair):
    ds = PanMelanoma()
    lair.safe_derive(ds)
    filepaths = lair.get_dataset_filepaths(ds)
    annotations = pd.read_csv(filepaths["annotations.tsv"], sep="\t").set_index(
        "Sample")
    signature = pd.read_csv(filepaths["signature.tsv"], sep="\t", index_col=0)
    signature = signature.loc[np.isnan(signature).sum(axis=1) == 0]
    return annotations, signature