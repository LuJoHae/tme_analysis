import anndata as ad
import scanpy as sc
import tcga
import icir
import datalair
from gene_utils import read_gene_sets, ssgsea_formula, median_scale, \
    download_github_file
from ici_datasets.bagaev_datasets import Signature


def load_tcga(lair) -> ad.AnnData:
    adatas = []
    for name, path in lair.get_dataset_filepaths(tcga.AllProjectsAdata()).items():
        adata = ad.read_h5ad(path).T
        adata.obs["dataset"] = name.removesuffix(".h5ad")
        adatas.append(adata)
    adata = ad.concat(adatas, axis=0)
    return adata


def get_vanallen_bagaev_signature(lair):
    ds = Signature()
    lair.safe_derive(ds)
    filepaths = lair.get_dataset_filepaths(ds)
    bagaev_signature = read_gene_sets(filepaths["gene_signatures.gmt"])

    ds = icir.datasets.ImmuneCheckpointTherapyResponseProcessedGeneNormalizedClinicalDataNormalized()
    filepath = lair.get_dataset_filepaths(ds)["VanAllen.h5ad"]
    vanallen = ad.read_h5ad(filepath)

    sc.pp.log1p(vanallen)
    genes = set.union(*[x.genes for x in bagaev_signature.values()])
    vanallen = vanallen[:, vanallen.var["gene_name"].isin(genes)]
    vanallen.var.set_index("gene_name", inplace=True)
    vanallen.obs.set_index("original.patient_id", inplace=True)
    vanallen = vanallen.to_df()
    vanallen_signatures = ssgsea_formula(vanallen, bagaev_signature)
    vanallen_signatures_scaled = median_scale(vanallen_signatures, 4)
    return vanallen_signatures_scaled


class DatasetVanAllen(datalair.Dataset):
    """Datalair Dataset class for all single cell datasets."""

    def __init__(self) -> None:
        """Initialize this dataset class as a datalair.Dataset class with namespace `DatasetVanAllen`."""
        super().__init__(namespace="DatasetVanAllen")


class VanAllenCTLA4RNAseqTPM(DatasetVanAllen):
    def derive(self, lair: datalair.Lair) -> None:
        output_dir = lair.get_path(self)
        download_github_file(
            owner="vanallenlab",
            repo="VanAllen_CTLA4_Science_RNASeq_TPM",
            commit_hash="ca459ce4a3eaeffdcd47454064afc2cb319b74ce",  # First git commit
            file_path="TPM_RSEM_VAScience2015.txt",
            output_filename=output_dir / "TPM_RSEM_VAScience2015.txt"
        )
